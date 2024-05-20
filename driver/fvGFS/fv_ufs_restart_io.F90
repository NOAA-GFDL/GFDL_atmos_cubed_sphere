module fv_ufs_restart_io_mod

  use esmf
  use mpp_mod,            only: mpp_pe, mpp_root_pe, mpp_chksum, mpp_npes, mpp_get_current_pelist
  use fv_arrays_mod,      only: fv_atmos_type
  use tracer_manager_mod, only: get_tracer_names
  use field_manager_mod,  only: MODEL_ATMOS
  use atmosphere_mod,     only: atmosphere_resolution
  use fv_io_mod,          only: fv_io_write_BCs

  implicit none

  public fv_dyn_restart_register
  public fv_dyn_restart_output

  public fv_core_restart_bundle_setup
  public fv_srf_wnd_restart_bundle_setup
  public fv_tracer_restart_bundle_setup

  private

  ! fv_core.res
  integer :: nvar3d_core_center, core_zsize
  real, allocatable, target, dimension(:,:,:) :: core_center_var2
  real, allocatable, target, dimension(:,:,:,:) :: core_center_var3, core_east_var3, core_north_var3
  character(len=32), allocatable,dimension(:)  :: core_center_var2_names, core_center_var3_names, core_east_var3_names, core_north_var3_names

  ! fv_srf_wnd.res
  integer :: nvar2d_srf_wnd
  real, allocatable, target, dimension(:,:,:) :: srf_wnd_var2
  character(len=32), allocatable, dimension(:) :: srf_wnd_var2_names

  ! fv_tracers
  integer :: nvar3d_tracers
  integer :: ntprog, ntdiag, tracers_zsize
  real, allocatable, target, dimension(:,:,:,:) :: tracers_var3
  character(len=32), allocatable, dimension(:) :: tracers_var3_names

  type(ESMF_FieldBundle) :: core_bundle, srf_wnd_bundle, tracer_bundle

 contains

 subroutine fv_dyn_restart_register (Atm)

   ! this subroutine must allocate all data buffers and set the variable names
   ! for restart bundles
   ! must be consistent with fv_io_register_restart in atmos_cubed_sphere/tools/fv_io.F90

   implicit none

   type(fv_atmos_type), intent(inout) :: Atm

   integer :: isc, iec, jsc, jec, nx, ny, nx_c, ny_c, nx_s, ny_s, nx_e, ny_e
   integer :: nlon, nlat, mlon, mlat
   integer :: num, nt, n
   character(len=64) :: tracer_name

   ! if(mpp_pe() == mpp_root_pe()) write(0,*) ' in fv_dyn_restart_register'

   isc = Atm%bd%isc
   iec = Atm%bd%iec
   jsc = Atm%bd%jsc
   jec = Atm%bd%jec
   nx  = (iec - isc + 1)
   ny  = (jec - jsc + 1)

   ! core
   ! 'u', 'v'
   nvar3d_core_center = 0
   if (.not.atm%flagstruct%hydrostatic) then
     nvar3d_core_center = nvar3d_core_center + 2 ! 'w', 'dz'
     if (atm%flagstruct%hybrid_z) then
       nvar3d_core_center = nvar3d_core_center + 1 ! 'ze0'
     endif
   endif
   nvar3d_core_center = nvar3d_core_center + 2 ! 't', 'delp'
   !--- include agrid winds in restarts for use in data assimilation
   if (atm%flagstruct%agrid_vel_rst) then
      nvar3d_core_center = nvar3d_core_center + 2 ! 'ua', 'va'
   endif
   core_zsize = size(Atm%u,3)

   call atmosphere_resolution (mlon, mlat, global=.true.)
   nx_c = nx
   ny_c = ny
   if (iec == mlon) then
      ! we are on task at the 'east' edge of the cubed sphere face or regional domain
      ! corner arrays should have one extra element in 'i' direction
      nx_c = nx + 1
   end if
   if (jec == mlat) then
      ! we are on task at the 'north' edge of the cubed sphere face or regional domain
      ! corner arrays should have one extra element in 'j' direction
      ny_c = ny + 1
   end if

   allocate(core_north_var3  (nx  , ny_c,core_zsize,1))
   allocate(core_east_var3   (nx_c, ny  ,core_zsize,1))
   allocate(core_center_var3 (nx  , ny  ,core_zsize,nvar3d_core_center))
   allocate(core_center_var2 (nx  , ny  ,1))

   core_north_var3  = 0.0
   core_east_var3   = 0.0
   core_center_var3 = 0.0
   core_center_var2 = 0.0

   allocate(core_north_var3_names  (1))
   allocate(core_east_var3_names   (1))
   allocate(core_center_var3_names (nvar3d_core_center))
   allocate(core_center_var2_names (1))

   core_north_var3_names(1) = 'u'
   core_east_var3_names(1) = 'v'
   n = 1
   if (.not.atm%flagstruct%hydrostatic) then
     core_center_var3_names(n) = 'W'; n=n+1
     core_center_var3_names(n) = 'DZ'; n=n+1
     if (atm%flagstruct%hybrid_z) then
       core_center_var3_names(n) = 'ze0'; n=n+1
     endif
   endif
   core_center_var3_names(n) = 'T'; n=n+1
   core_center_var3_names(n) = 'delp'; n=n+1
   !--- include agrid winds in restarts for use in data assimilation
   if (atm%flagstruct%agrid_vel_rst) then
     core_center_var3_names(n) = 'ua'; n=n+1
     core_center_var3_names(n) = 'va'; n=n+1
   endif

   core_center_var2_names(1) = 'phis'

   ! srf_wnd
   nvar2d_srf_wnd = 2
   allocate (srf_wnd_var2(nx,ny,nvar2d_srf_wnd), srf_wnd_var2_names(nvar2d_srf_wnd))
   srf_wnd_var2 = 0.0
   srf_wnd_var2_names(1) = 'u_srf'
   srf_wnd_var2_names(2) = 'v_srf'

   ! tracers
   ntprog = size(Atm%q,4)
   ntdiag = size(Atm%qdiag,4)
   nvar3d_tracers = ntprog+ntdiag
   tracers_zsize = size(Atm%q,3)
   allocate (tracers_var3(nx,ny,tracers_zsize,nvar3d_tracers), tracers_var3_names(nvar3d_tracers))
   tracers_var3 = 0.0

   do nt = 1, ntprog
      call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
      tracers_var3_names(nt) = tracer_name
   enddo
   do nt = ntprog+1, nvar3d_tracers
      call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
      tracers_var3_names(nt) = tracer_name
   enddo

 end subroutine fv_dyn_restart_register

 subroutine fv_dyn_restart_output(Atm, timestamp)

   implicit none

   type(fv_atmos_type), intent(inout) :: Atm
   character(len=*), intent(in) :: timestamp

   integer :: isc, iec, jsc, jec, nx, ny
   integer :: n, nt
   integer :: mlon, mlat, nx_c, ny_c

   isc = Atm%bd%isc
   iec = Atm%bd%iec
   jsc = Atm%bd%jsc
   jec = Atm%bd%jec
   nx  = (iec - isc + 1)
   ny  = (jec - jsc + 1)

   call atmosphere_resolution (mlon, mlat, global=.true.)
   nx_c = nx
   ny_c = ny
   if (iec == mlon) then
      ! we are on task at the 'east' edge of the cubed sphere face or regional domain
      ! corner arrays should have one extra element in 'i' direction
      nx_c = nx + 1
   end if
   if (jec == mlat) then
      ! we are on task at the 'north' edge of the cubed sphere face or regional domain
      ! corner arrays should have one extra element in 'j' direction
      ny_c = ny + 1
   end if

   ! ---- fv_core.res

   core_north_var3(:,:,:,1) = Atm%u(isc:iec,jsc:(jsc+ny_c-1),:)
   core_east_var3 (:,:,:,1) = Atm%v(isc:(isc+nx_c-1),jsc:jec,:)

   n = 1
   if (.not.atm%flagstruct%hydrostatic) then
      core_center_var3(:,:,:,n) = Atm%w(isc:iec,jsc:jec,:); n=n+1
      core_center_var3(:,:,:,n) = Atm%delz; n=n+1
      if (atm%flagstruct%hybrid_z) then
        core_center_var3(:,:,:,n) = Atm%ze0; n=n+1
      endif
   endif
   core_center_var3(:,:,:,n) = Atm%pt(isc:iec,jsc:jec,:);   n=n+1
   core_center_var3(:,:,:,n) = Atm%delp(isc:iec,jsc:jec,:); n=n+1
   !--- include agrid winds in restarts for use in data assimilation
   if (atm%flagstruct%agrid_vel_rst) then
     core_center_var3(:,:,:,n) = Atm%ua(isc:iec,jsc:jec,:); n=n+1
     core_center_var3(:,:,:,n) = Atm%va(isc:iec,jsc:jec,:); n=n+1
   endif

   core_center_var2(:,:,1) = Atm%phis(isc:iec,jsc:jec)

   ! ---- fv_srf_wnd.res
   srf_wnd_var2(:,:,1) = Atm%u_srf
   srf_wnd_var2(:,:,2) = Atm%v_srf

   ! ---- fv_tracer_res
   do nt = 1, ntprog
     tracers_var3(:,:,:,nt) = Atm%q(isc:iec,jsc:jec,:,nt)
   enddo
   do nt = ntprog+1, nvar3d_tracers
     tracers_var3(:,:,:,nt) = Atm%qdiag(isc:iec,jsc:jec,:,nt)
   enddo

   ! Instead of creating yet another esmf bundle just to write Atm%ak and Atm%bk, write them here synchronously
   call write_ak_bk(Atm, timestamp)

   ! Write bcs restarts
   if (Atm%neststruct%nested) then
     call fv_io_write_BCs(Atm)
   endif

 end subroutine fv_dyn_restart_output

 subroutine fv_core_restart_bundle_setup(bundle, grid, rc)
!
!-------------------------------------------------------------
!*** set esmf bundle for fv_core restart fields
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

   core_bundle = bundle

   call ESMF_FieldBundleGet(bundle, name=bdl_name,rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

   outputfile = trim(bdl_name)

   temp_r3d => core_north_var3(:,:,:,1)
   call create_3d_field_and_add_to_bundle(temp_r3d, trim(core_north_var3_names(1)), "zaxis_1", core_zsize, trim(outputfile), grid, bundle, staggerloc=ESMF_STAGGERLOC_EDGE2)

   temp_r3d => core_east_var3(:,:,:,1)
   call create_3d_field_and_add_to_bundle(temp_r3d, trim(core_east_var3_names(1)), "zaxis_1", core_zsize, trim(outputfile), grid, bundle, staggerloc=ESMF_STAGGERLOC_EDGE1)

   do n = 1, nvar3d_core_center
     temp_r3d => core_center_var3(:,:,:,n)
     call create_3d_field_and_add_to_bundle(temp_r3d, trim(core_center_var3_names(n)), "zaxis_1", core_zsize, trim(outputfile), grid, bundle)
   end do

   temp_r2d => core_center_var2(:,:,1)
   call create_2d_field_and_add_to_bundle(temp_r2d, trim(core_center_var2_names(1)), trim(outputfile), grid, bundle)

 end subroutine fv_core_restart_bundle_setup

 subroutine fv_srf_wnd_restart_bundle_setup(bundle, grid, rc)
!
!-------------------------------------------------------------
!*** set esmf bundle for fv_srf_wnd restart fields
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

   srf_wnd_bundle = bundle

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
!*** set esmf bundle for fv_tracer restart fields
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

   tracer_bundle = bundle

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

 subroutine create_3d_field_and_add_to_bundle(temp_r3d, field_name, axis_name, num_levels, outputfile, grid, bundle, staggerloc)

   use esmf

   implicit none

   real, dimension(:,:,:), pointer, intent(in)    :: temp_r3d
   character(len=*),                           intent(in)    :: field_name
   character(len=*),                           intent(in)    :: axis_name
   integer,                                    intent(in)    :: num_levels
   character(len=*),                           intent(in)    :: outputfile
   type(ESMF_Grid),                            intent(in)    :: grid
   type(ESMF_FieldBundle),                     intent(inout) :: bundle
   type(ESMF_StaggerLoc),optional,             intent(in)    :: staggerloc

   type(ESMF_Field) :: field
   type(ESMF_StaggerLoc) :: stagger

   integer :: rc, i

   stagger = ESMF_STAGGERLOC_CENTER
   if(present(staggerloc)) then
     stagger = staggerloc
   end if

   field = ESMF_FieldCreate(grid, temp_r3d, datacopyflag=ESMF_DATACOPY_REFERENCE, &
                            name=trim(field_name), indexFlag=ESMF_INDEX_DELOCAL, staggerloc=stagger, rc=rc)
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

   real, allocatable, dimension(:) :: buffer
   integer :: rc, i

   call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
                          attrList=(/"ESMF:ungridded_dim_labels"/), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
   call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
                          name="ESMF:ungridded_dim_labels", valueList=(/trim(axis_name)/), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

   call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3-dim",  &
                          attrList=(/trim(axis_name)//"               ",trim(axis_name)//":cartesian_axis"/), rc=rc)
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
                          name=trim(axis_name)//"axis", value="Z", rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

 end subroutine add_zaxis_to_field

 subroutine write_ak_bk(Atm, timestamp)


  use fms2_io_mod,             only: FmsNetcdfFile_t, &
                                     register_restart_field, register_axis, unlimited, &
                                     open_file, write_restart, &
                                     close_file, register_field, write_data, &
                                     register_variable_attribute
   implicit none

   type(fv_atmos_type), intent(in) :: Atm
   character(len=*), intent(in) :: timestamp

   character(len=120) :: fname
   integer, allocatable, dimension(:) :: pes
   type(FmsNetcdfFile_t) :: Fv_restart
   logical :: Fv_restart_is_open
   integer, dimension(:), allocatable :: buffer
   character(len=8), dimension(2)  :: dim_names_2d
   integer :: j

#ifdef OVERLOAD_R4
   character(len=5), parameter :: axis_type = 'float'
#else
   character(len=6), parameter :: axis_type = 'double'
#endif

   dim_names_2d(1) = "xaxis_1"
   dim_names_2d(2) = "Time"

   if (trim(timestamp) == '') then
      fname = 'RESTART/'//'fv_core.res.nc'
   else
      fname = 'RESTART/'//trim(timestamp)//'.fv_core.res.nc'
   endif

   allocate(pes(mpp_npes()))
   call mpp_get_current_pelist(pes)

   Fv_restart_is_open = open_file(Fv_restart, fname, "overwrite", is_restart=.true., pelist=pes)

   call register_axis(Fv_restart, "xaxis_1", size(Atm%ak(:), 1))
   call register_axis(Fv_restart, "Time", unlimited)
   call register_field(Fv_restart, "xaxis_1", axis_type, (/"xaxis_1"/))
   call register_variable_attribute(Fv_restart,"xaxis_1", "axis", "X", str_len=1)
   if (allocated(buffer)) deallocate(buffer)
   allocate(buffer(size(Atm%ak(:), 1)))
   do j = 1, size(Atm%ak(:), 1)
      buffer(j) = j
   end do
   call write_data(Fv_restart, "xaxis_1", buffer)
   deallocate(buffer)
   call register_field(Fv_restart, "Time", axis_type, (/"Time"/))
   call register_variable_attribute(Fv_restart, dim_names_2d(2), "cartesian_axis", "T", str_len=1)
   call register_variable_attribute(Fv_restart, dim_names_2d(2), "units", "time level", str_len=len("time level"))
   call register_variable_attribute(Fv_restart, dim_names_2d(2), "long_name", dim_names_2d(2), str_len=len(dim_names_2d(2)))
   call write_data(Fv_restart, "Time", 1)
   call register_restart_field (Fv_restart, 'ak', Atm%ak(:), dim_names_2d)
   call register_restart_field (Fv_restart, 'bk', Atm%bk(:), dim_names_2d)

   call write_restart(Fv_restart)
   call close_file(Fv_restart)
   deallocate(pes)

 end subroutine write_ak_bk

end module fv_ufs_restart_io_mod
