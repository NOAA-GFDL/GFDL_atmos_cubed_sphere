module fv_ufs_restart_io_mod

  use mpp_mod,            only: mpp_pe, mpp_root_pe, mpp_chksum, mpp_pe
  use fv_arrays_mod,      only: fv_atmos_type
  use tracer_manager_mod,      only: tr_get_tracer_names=>get_tracer_names, &
                                     get_tracer_names, get_number_tracers, &
                                     set_tracer_profile, &
                                     get_tracer_index
  use field_manager_mod,       only: MODEL_ATMOS

  implicit none

  public fv_ufs_restart_register

  public fv_srf_wnd_restart_output
  public fv_srf_wnd_restart_bundle_setup

  public fv_tracer_restart_output
  public fv_tracer_restart_bundle_setup

  private

  ! fv_core.res

  ! fv_srf_wnd.res
  integer :: nvar2d_srf_wnd
  real, allocatable, target, dimension(:,:,:) :: srf_wnd_var2
  character(32), allocatable,dimension(:) :: srf_wnd_var2_names

  ! fv_tracers
  integer :: nvar3d_tracers
  integer :: ntprog, ntdiag, tracers_zsize
  real, allocatable, target, dimension(:,:,:,:) :: tracers_var3
  character(32), allocatable,dimension(:) :: tracers_var3_names

 contains

 subroutine fv_ufs_restart_register (Atm)

   ! this subroutine must allocate all data buffers and set the variable names
   ! for restart bundles

   implicit none

   type(fv_atmos_type), intent(inout) :: Atm

   integer :: isc, iec, jsc, jec, nx, ny
   integer :: num, nt
   character(len=64) :: tracer_name

   isc = Atm%bd%isc
   iec = Atm%bd%iec
   jsc = Atm%bd%jsc
   jec = Atm%bd%jec
   nx  = (iec - isc + 1)
   ny  = (jec - jsc + 1)

   if(mpp_pe() == mpp_root_pe()) write(0,*) ' in fv_ufs_restart_register '

   ! srf_wnd
   nvar2d_srf_wnd = 2
   allocate (srf_wnd_var2(nx,ny,nvar2d_srf_wnd), srf_wnd_var2_names(nvar2d_srf_wnd))
   srf_wnd_var2_names(1) = 'u_srf'
   srf_wnd_var2_names(2) = 'v_srf'

   ! tracers
   ntprog = size(Atm%q,4)
   ntdiag = size(Atm%qdiag,4)
   nvar3d_tracers = ntprog+ntdiag
   tracers_zsize = size(Atm%q,3)
   allocate (tracers_var3(nx,ny,tracers_zsize,nvar3d_tracers), tracers_var3_names(nvar3d_tracers))

   do nt = 1, ntprog
      call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
      tracers_var3_names(nt) = tracer_name
   enddo
   do nt = ntprog+1, nvar3d_tracers
      call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
      tracers_var3_names(nt) = tracer_name
   enddo

 end subroutine fv_ufs_restart_register

 subroutine fv_srf_wnd_restart_output(Atm)

   implicit none

   type(fv_atmos_type), intent(inout) :: Atm

   if(mpp_pe() == mpp_root_pe()) write(0,*) ' in fv_srf_wnd_restart_output '

   srf_wnd_var2(:,:,1) = Atm%u_srf
   srf_wnd_var2(:,:,2) = Atm%v_srf

 end subroutine fv_srf_wnd_restart_output

 subroutine fv_tracer_restart_output(Atm)

   implicit none

   type(fv_atmos_type), intent(inout) :: Atm

   integer :: nt
   integer :: isc, iec, jsc, jec

   if(mpp_pe() == mpp_root_pe()) write(0,*) ' in fv_tracer_restart_output '

   isc = Atm%bd%isc
   iec = Atm%bd%iec
   jsc = Atm%bd%jsc
   jec = Atm%bd%jec

   do nt = 1, ntprog
     tracers_var3(:,:,:,nt) = Atm%q(isc:iec,jsc:jec,:,nt)
   enddo
   do nt = ntprog+1, nvar3d_tracers
     tracers_var3(:,:,:,nt) = Atm%qdiag(isc:iec,jsc:jec,:,nt)
   enddo

 end subroutine fv_tracer_restart_output

 subroutine fv_srf_wnd_restart_bundle_setup(bundle, grid, rc)
!
!-------------------------------------------------------------
!*** set esmf bundle for phys restart fields
!------------------------------------------------------------
!
   use esmf

   implicit none

   type(ESMF_FieldBundle),intent(inout)        :: bundle
   type(ESMF_Grid),intent(inout)               :: grid
   integer,intent(out)                         :: rc

!*** local variables
   integer i, j, k, n
   character(128)    :: bdl_name
   type(ESMF_Field)  :: field
   character(128)    :: outputfile
   integer :: num
   real,dimension(:,:),pointer   :: temp_r2d
   real,dimension(:,:,:),pointer   :: temp_r3d

   if(mpp_pe() == mpp_root_pe()) write(0,*) ' in fv_srf_wnd_restart_bundle_setup '

   call ESMF_FieldBundleGet(bundle, name=bdl_name,rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

   outputfile = trim(bdl_name)

   do num = 1,nvar2d_srf_wnd
       temp_r2d => srf_wnd_var2(:,:,num)
       call create_2d_field_and_add_to_bundle(temp_r2d, trim(srf_wnd_var2_names(num)), trim(outputfile), grid, bundle)
   enddo

 end subroutine fv_srf_wnd_restart_bundle_setup

 subroutine fv_tracer_restart_bundle_setup(bundle, grid, rc)
!
!-------------------------------------------------------------
!*** set esmf bundle for phys restart fields
!------------------------------------------------------------
!
   use esmf

   implicit none

   type(ESMF_FieldBundle),intent(inout)        :: bundle
   type(ESMF_Grid),intent(inout)               :: grid
   integer,intent(out)                         :: rc

!*** local variables
   integer i, j, k, n
   character(128)    :: bdl_name
   type(ESMF_Field)  :: field
   character(128)    :: outputfile
   integer :: num
   real,dimension(:,:),pointer   :: temp_r2d
   real,dimension(:,:,:),pointer   :: temp_r3d

   if(mpp_pe() == mpp_root_pe()) write(0,*) ' in fv_tracer_restart_bundle_setup '

   call ESMF_FieldBundleGet(bundle, name=bdl_name,rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

   outputfile = trim(bdl_name)

   do num = 1,nvar3d_tracers
       temp_r3d => tracers_var3(:,:,:,num)
       call create_3d_field_and_add_to_bundle(temp_r3d, trim(tracers_var3_names(num)), "zaxis_1", tracers_zsize, trim(outputfile), grid, bundle)
   enddo

 end subroutine fv_tracer_restart_bundle_setup

 subroutine create_2d_field_and_add_to_bundle(temp_r2d, field_name, outputfile, grid, bundle)

   use esmf

   implicit none

   real, dimension(:,:),   pointer, intent(in)    :: temp_r2d
   character(len=*),                           intent(in)    :: field_name
   character(len=*),                           intent(in)    :: outputfile
   type(ESMF_Grid),                            intent(in)    :: grid
   type(ESMF_FieldBundle),                     intent(inout) :: bundle

   type(ESMF_Field) :: field

   integer :: rc, i

   field = ESMF_FieldCreate(grid, temp_r2d, datacopyflag=ESMF_DATACOPY_REFERENCE, &
                            name=trim(field_name), indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, file=__FILE__)) &
   call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", attrList=(/"output_file"/), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", name='output_file', value=trim(outputfile), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_FieldBundleAdd(bundle, (/field/), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  end subroutine create_2d_field_and_add_to_bundle

 subroutine create_3d_field_and_add_to_bundle(temp_r3d, field_name, axis_name, num_levels, outputfile, grid, bundle)

   use esmf

   implicit none

   real, dimension(:,:,:), pointer, intent(in)    :: temp_r3d
   character(len=*),                           intent(in)    :: field_name
   character(len=*),                           intent(in)    :: axis_name
   integer,                                    intent(in)    :: num_levels
   character(len=*),                           intent(in)    :: outputfile
   type(ESMF_Grid),                            intent(in)    :: grid
   type(ESMF_FieldBundle),                     intent(inout) :: bundle

   type(ESMF_Field) :: field

   integer :: rc, i

   ! if(mpp_pe() == mpp_root_pe()) write(0,*)'in add 3d field to bundle ', trim(field_name)
   field = ESMF_FieldCreate(grid, temp_r3d, datacopyflag=ESMF_DATACOPY_REFERENCE, &
                            name=trim(field_name), indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, file=__FILE__)) &
   call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", attrList=(/"output_file"/), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", name='output_file', value=trim(outputfile), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call add_zaxis_to_field(field, axis_name, num_levels)

   call ESMF_FieldBundleAdd(bundle, (/field/), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  end subroutine create_3d_field_and_add_to_bundle

 subroutine add_zaxis_to_field(field, axis_name, num_levels)

   use esmf

   implicit none

   type(ESMF_Field), intent(inout) :: field
   character(len=*), intent(in)    :: axis_name
   integer,          intent(in)    :: num_levels

   real(8), allocatable, dimension(:) :: buffer
   integer :: rc, i

   call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
                          attrList=(/"ESMF:ungridded_dim_labels"/), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
   call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
                          name="ESMF:ungridded_dim_labels", valueList=(/trim(axis_name)/), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

   call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3-dim",  &
                          attrList=(/trim(axis_name),trim(axis_name)//":cartesian_axis"/), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

   allocate( buffer(num_levels) )
   do i=1, num_levels
      buffer(i)=i
   end do
   call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3-dim", &
                          name=trim(axis_name), valueList=buffer, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
   deallocate(buffer)

   call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3-dim", &
                          name=trim(axis_name)//"cartesian_axis", value="Z", rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

 end subroutine add_zaxis_to_field


end module fv_ufs_restart_io_mod
