! $Id$
!
!----------------
! FV contro panel
!----------------

module fv_control_mod

   use constants_mod,       only: pi, kappa, radius, grav, rdgas
   use field_manager_mod,   only: MODEL_ATMOS
   use fms_mod,             only: write_version_number, open_namelist_file, &
                                  check_nml_error, close_file, file_exist
   use mpp_mod,             only: FATAL, mpp_error, mpp_pe, stdlog, &
                                  mpp_npes, mpp_get_current_pelist, &
                                  input_nml_file, get_unit, WARNING, &
                                  read_ascii_file, INPUT_STR_LENGTH, mpp_pe
   use mpp_domains_mod,     only: mpp_get_data_domain, mpp_get_compute_domain
   use tracer_manager_mod,  only: tm_get_number_tracers => get_number_tracers, &
                                  tm_get_tracer_index   => get_tracer_index,   &
                                  tm_get_tracer_indices => get_tracer_indices, &
                                  tm_set_tracer_profile => set_tracer_profile, &
                                  tm_get_tracer_names   => get_tracer_names,   &
                                  tm_check_if_prognostic=> check_if_prognostic,&
                                  tm_register_tracers   => register_tracers

   use fv_io_mod,           only: fv_io_exit
   use fv_restart_mod,      only: fv_restart_init, fv_restart_end
   use fv_arrays_mod,       only: fv_atmos_type, allocate_fv_atmos_type, deallocate_fv_atmos_type
   use fv_grid_utils_mod,   only: grid_utils_init, grid_utils_end, ptop_min
   use fv_eta_mod,          only: set_eta
   use fv_grid_tools_mod,   only: init_grid, read_grid, broadcast_aligned_nest
   use fv_mp_mod,           only: mp_start, mp_assign_gid, domain_decomp
   use fv_mp_mod,           only: ng, switch_current_Atm
   use fv_mp_mod,           only: broadcast_domains, mp_barrier, is_master, setup_master
!!! CLEANUP: should be replaced by a getter function?
   use test_cases_mod,      only: test_case, alpha
   use fv_timing_mod,       only: timing_on, timing_off, timing_init, timing_prt
   use mpp_domains_mod,     only: domain2D
   use mpp_domains_mod,     only: mpp_define_nest_domains, nest_domain_type, mpp_get_global_domain
   use mpp_domains_mod,     only: mpp_get_C2F_index, mpp_get_F2C_index, mpp_broadcast_domain
   use mpp_domains_mod,     only: CENTER, CORNER, NORTH, EAST, WEST, SOUTH
   use mpp_mod,             only: mpp_send, mpp_sync, mpp_transmit, mpp_set_current_pelist, mpp_declare_pelist, mpp_root_pe, mpp_recv, mpp_sync_self, mpp_broadcast, read_input_nml

   implicit none
   private

!-----------------------------------------------------------------------
! Grid descriptor file setup
!-----------------------------------------------------------------------
!------------------------------------------
! Model Domain parameters
! See fv_arrays.F90 for descriptions
!------------------------------------------
!CLEANUP module pointers
   character(len=80) , pointer :: grid_name
   character(len=120), pointer :: grid_file
   integer, pointer :: grid_type
   integer , pointer :: hord_mt 
   integer , pointer :: kord_mt 
   integer , pointer :: kord_wz 
   integer , pointer :: hord_vt 
   integer , pointer :: hord_tm 
   integer , pointer :: hord_dp 
   integer , pointer :: kord_tm 
   integer , pointer :: hord_tr 
   integer , pointer :: kord_tr 
   real    , pointer :: scale_z 
   real    , pointer :: w_max 
   real    , pointer :: z_min 

   integer , pointer :: nord
   real    , pointer :: dddmp 
   real    , pointer :: d2_bg 
   real    , pointer :: d4_bg 
   real    , pointer :: vtdm4 
   real    , pointer :: d2_bg_k1 
   real    , pointer :: d2_bg_k2 
   real    , pointer :: d2_divg_max_k1 
   real    , pointer :: d2_divg_max_k2 
   real    , pointer :: damp_k_k1 
   real    , pointer :: damp_k_k2 
   integer , pointer ::    n_zs_filter
   integer , pointer :: nord_zs_filter

   logical , pointer :: consv_am
   logical , pointer :: do_sat_adj
   logical , pointer :: do_f3d
   logical , pointer :: no_dycore 
   logical , pointer :: replace_w 
   logical , pointer :: convert_ke 
   logical , pointer :: do_vort_damp 
   logical , pointer :: use_old_omega 
! PG off centering:
   real    , pointer :: beta  
   integer , pointer :: n_sponge 
   real    , pointer :: d_ext 
   integer , pointer :: nwat  
   logical , pointer :: warm_start 
   logical , pointer :: inline_q 
   real , pointer :: shift_fac   
   logical , pointer :: do_schmidt 
   real , pointer :: stretch_fac 
   real , pointer :: target_lat  
   real , pointer :: target_lon  

   logical , pointer :: reset_eta 
   real    , pointer :: p_fac
   real    , pointer :: a_imp
   integer , pointer :: n_split 
                             ! Default 
   integer , pointer :: m_split 
   integer , pointer :: k_split 
   logical , pointer :: use_logp

   integer , pointer :: q_split 
   integer , pointer :: print_freq 

   integer , pointer :: npx           
   integer , pointer :: npy           
   integer , pointer :: npz           
   integer , pointer :: npz_rst 
                                      
   integer , pointer :: ncnst 
   integer , pointer :: pnats 
   integer , pointer :: dnats 
   integer , pointer :: ntiles        
   integer , pointer :: nf_omega  
   integer , pointer :: fv_sg_adj 
                                      
   integer , pointer :: na_init 
   real    , pointer :: p_ref 
   real    , pointer :: dry_mass 
   integer , pointer :: nt_prog 
   integer , pointer :: nt_phys 
   real    , pointer :: tau_h2o 

   real    , pointer :: d_con 
   real    , pointer :: consv_te 
   real    , pointer :: tau 
   real    , pointer :: rf_cutoff
   logical , pointer :: filter_phys 
   logical , pointer :: dwind_2d 
   logical , pointer :: breed_vortex_inline 
   logical , pointer :: range_warn 
   logical , pointer :: fill 
   logical , pointer :: fill_dp 
   logical , pointer :: fill_wz 
   logical , pointer :: check_negative
   logical , pointer :: non_ortho 
   logical , pointer :: adiabatic 
   logical , pointer :: moist_phys 
   logical , pointer :: do_Held_Suarez 
   logical , pointer :: do_reed_physics
   logical , pointer :: reed_cond_only
   logical , pointer :: reproduce_sum 
   logical , pointer :: adjust_dry_mass 
   logical , pointer :: fv_debug  
   logical , pointer :: srf_init  
   logical , pointer :: mountain  
   logical , pointer :: remap_t  
   logical , pointer :: z_tracer 

   logical , pointer :: old_divg_damp 
   logical , pointer :: fv_land 
   logical , pointer :: nudge 
   logical , pointer :: nudge_ic
   logical , pointer :: ncep_ic 
   logical , pointer :: fv_diag_ic 
   logical , pointer :: external_ic 
   character(len=128) , pointer :: res_latlon_dynamics
   character(len=128) , pointer :: res_latlon_tracers 
   logical , pointer :: hydrostatic 
   logical , pointer :: phys_hydrostatic 
   logical , pointer :: hybrid_z    
   logical , pointer :: Make_NH     
   logical , pointer :: make_hybrid_z  
   real,     pointer :: add_noise

   integer , pointer :: a2b_ord 
   integer , pointer :: c2l_ord 

   integer, pointer :: ndims

  real, pointer :: dx_const
  real, pointer :: dy_const
  real, pointer :: deglon_start, deglon_stop, &  ! boundaries of latlon patch
          deglat_start, deglat_stop
  real, pointer :: deglat

  logical, pointer :: nested, twowaynest
  integer, pointer :: parent_tile, refinement, nestbctype, nestupdate, nsponge, ioffset, joffset
  real, pointer :: s_weight

  integer, pointer :: layout(:), io_layout(:)

   integer :: ntilesMe                ! Number of tiles on this process =1 for now

   real    :: too_big  = 1.E35
   public :: fv_init, fv_end

   integer, public :: ngrids = 1
   integer, public, allocatable :: pelist_all(:)
   integer :: commID, max_refinement_of_global = 1.
   integer :: gid

!---- version number -----
   character(len=128) :: version = '$Id$'
   character(len=128) :: tagname = '$Name$'


   integer dumm
   real :: umax = 350.           ! max wave speed for grid_type>3
   integer :: parent_grid_num = -1

 contains

!-------------------------------------------------------------------------------
         
 subroutine fv_init(Atm, dt_atmos, grids_on_this_pe, p_split)

   type(fv_atmos_type), allocatable, intent(inout), target :: Atm(:)
   real,                intent(in)    :: dt_atmos
   logical, allocatable, intent(INOUT) :: grids_on_this_pe(:)
   integer, intent(INOUT) :: p_split

   integer :: i, j, k, n, p
   real :: sdt

! tracers
   integer :: num_family          ! output of register_tracers

   integer :: isd_p, ied_p, jsd_p, jed_p, isg, ieg, jsg, jeg
   integer :: ic, jc
   integer, allocatable :: p_ind(:,:,:)

   gid = mpp_pe()

   call init_nesting(Atm, grids_on_this_pe, p_split)

   !This call is needed to set up the pointers for fv_current_grid, even for a single-grid run
   !call switch_current_Atm(Atm(1), .false.)
   call setup_pointers(Atm(1))

! Start up MPI

   !call mp_assign_gid

    ! Initialize timing routines
      call timing_init
      call timing_on('TOTAL')

    ! Setup the run from namelist 
      ntilesMe = size(Atm(:)) !Full number of Atm arrays; one less than number of grids, if multiple grids

      call run_setup(Atm,dt_atmos, grids_on_this_pe, p_split)   ! initializes domain_decomp

      do n=1,ntilesMe

         !In a single-grid run this will still be needed to correctly set the domain
         call switch_current_Atm(Atm(n))
         call setup_pointers(Atm(n))
         
         target_lon = target_lon * pi/180.
         target_lat = target_lat * pi/180.

!--------------------------------------------------
! override number of tracers by reading field_table
!--------------------------------------------------

         !not sure if this works with multiple grids
         call tm_register_tracers (MODEL_ATMOS, ncnst, nt_prog, pnats, num_family)
         if(is_master()) then
            write(*,*) 'ncnst=', ncnst,' num_prog=',nt_prog,' pnats=',pnats,' dnats=',dnats,' num_family=',num_family         
            print*, ''
         endif

         if (grids_on_this_pe(n)) then
            call allocate_fv_atmos_type(Atm(n), Atm(n)%bd%isd, Atm(n)%bd%ied, Atm(n)%bd%jsd, Atm(n)%bd%jed, Atm(n)%bd%isc, Atm(n)%bd%iec, Atm(n)%bd%jsc, Atm(n)%bd%jec, &
                 npx, npy, npz, ndims, ncnst, ncnst-pnats, ng, .false., grids_on_this_pe(n))

            if (grids_on_this_pe(n)) then

               call switch_current_Atm(Atm(n))
               call setup_pointers(Atm(n))

               if ( (Atm(n)%bd%iec-Atm(n)%bd%isc+1).lt.4 .or. (Atm(n)%bd%jec-Atm(n)%bd%jsc+1).lt.4 ) then
                  if (is_master()) write(*,'(6I6)') Atm(n)%bd%isc, Atm(n)%bd%iec, Atm(n)%bd%jsc, Atm(n)%bd%jec, n
                  call mpp_error(FATAL,'Domain Decomposition:  Cubed Sphere compute domain has a &
                       &minium requirement of 4 points in X and Y, respectively')
               end if

            endif

            !!CLEANUP: Convenience pointers
            Atm(n)%gridstruct%nested    => Atm(n)%neststruct%nested
            Atm(n)%gridstruct%grid_type => Atm(n)%flagstruct%grid_type
            Atm(n)%flagstruct%grid_number => Atm(n)%grid_number

            ! Read Grid from GRID_FILE and setup grid descriptors
            ! needs modification for multiple tiles
            if(grid_type <0 .AND. trim(grid_file) == 'INPUT/grid_spec.nc') then
               call read_grid(Atm(n), grid_name, grid_file, npx, npy, npz, ndims, ntiles, ng)
            else
               call init_grid(Atm(n), grid_name, grid_file, npx, npy, npz, ndims, ntiles, ng)
            endif

            ! Initialize the SW (2D) part of the model
            !!!CLEANUP: this call could definitely use some cleaning up
            call grid_utils_init(Atm(n), npx, npy, npz, Atm(n)%gridstruct%grid, Atm(n)%gridstruct%agrid,   &
                 Atm(n)%gridstruct%area, Atm(n)%gridstruct%area_c, &
                 Atm(n)%gridstruct%cosa, Atm(n)%gridstruct%sina,   &
                 Atm(n)%gridstruct%dx, Atm(n)%gridstruct%dy,       &
                 Atm(n)%gridstruct%dxa, Atm(n)%gridstruct%dya,     &
                 Atm(n)%gridstruct%dxc, Atm(n)%gridstruct%dyc, non_ortho,   &
                 grid_type, c2l_ord)


            !!!CLEANUP: Are these correctly writing out on all pes?
            if ( is_master() ) then
               sdt =  dt_atmos/real(n_split*k_split*p_split)
               write(*,*) ' '
               write(*,*) 'Divergence damping Coefficients'
               write(*,*) 'For small dt=', sdt
               write(*,*) 'External mode del-2 (m**2/s)=',  d_ext*Atm(n)%gridstruct%da_min_c/sdt
               write(*,*) 'Internal mode del-2 SMAG dimensionless coeff=',  dddmp
               write(*,*) 'Internal mode del-2 background diff=', d2_bg*Atm(n)%gridstruct%da_min_c/sdt

               if (nord==1) then
                   write(*,*) 'Internal mode del-4 background diff=', d4_bg
                   write(*,*) 'Vorticity del-4 (m**4/s)=', (vtdm4*Atm(n)%gridstruct%da_min)**2/sdt*1.E-6
               endif
               if (nord==2) write(*,*) 'Internal mode del-6 background diff=', d4_bg
               if (nord==3) write(*,*) 'Internal mode del-8 background diff=', d4_bg

               write(*,*) 'Vorticity del-4 (m**4/s)=', (vtdm4*Atm(n)%gridstruct%da_min)**2/sdt*1.E-6
               write(*,*) 'beta=', beta
               write(*,*) ' '
            endif


            Atm(n)%ts   = 300.
            Atm(n)%phis = too_big
            ! The following statements are to prevent the phatom corner regions from
            ! growing instability
            Atm(n)%u  = 0.
            Atm(n)%v  = 0.
            Atm(n)%ua = too_big
            Atm(n)%va = too_big

         else !this grid is NOT defined on this pe

            !Allocate dummy arrays
            call allocate_fv_atmos_type(Atm(n),  Atm(n)%bd%isd, Atm(n)%bd%ied, Atm(n)%bd%jsd, Atm(n)%bd%jed, Atm(n)%bd%isc, Atm(n)%bd%iec, Atm(n)%bd%jsc, Atm(n)%bd%jec, &
                 npx, npy, npz, ndims, ncnst, ncnst-pnats, ng, .true., .false.)

            !Need to SEND grid_global to any child grids; this is received in setup_aligned_nest in fv_grid_tools
            if (Atm(n)%neststruct%nested) then

               call mpp_get_global_domain( Atm(n)%parent_grid%domain, &
                    isg, ieg, jsg, jeg)

               if (gid == Atm(n)%parent_grid%pelist(1)) then
                  do p=1,size(Atm(n)%pelist)
                     call mpp_send(Atm(n)%parent_grid%grid_global(isg-ng:ieg+1+ng,jsg-ng:jeg+1+ng,1:2,parent_tile), &
                          size(Atm(n)%parent_grid%grid_global(isg-ng:ieg+1+ng,jsg-ng:jeg+1+ng,1:2,parent_tile)), &
                          Atm(n)%pelist(p)) !send to p_ind in setup_aligned_nest
                     call mpp_sync_self
                  enddo
               endif

               if (Atm(n)%neststruct%twowaynest) then

                  !Also need to receive P_IND from child grids to set up ind_update_h
                  isd_p = Atm(n)%parent_grid%bd%isd
                  ied_p = Atm(n)%parent_grid%bd%ied
                  jsd_p = Atm(n)%parent_grid%bd%jsd
                  jed_p = Atm(n)%parent_grid%bd%jed
                  !               allocate(Atm(n)%ind_update_h(isd_p:ied_p+1,jsd_p:jed_p+1,2))
                  allocate(p_ind(1-ng:npx+ng,1-ng:npy+ng,1:2))
                  call mpp_recv(p_ind,size(p_ind),Atm(n)%pelist(1)) !receiving from p_ind setup_aligned_grids 
                  call mpp_sync_self

                  Atm(n)%neststruct%ind_update_h = 1000000

                  if (Atm(n)%parent_grid%tile == Atm(n)%neststruct%parent_tile) then
                     do j=1,npy
                        do i=1,npx

                           ic = p_ind(i,j,1)
                           jc = p_ind(i,j,2)

                           if (ic < isd_p .or. ic > ied_p .or. jc < jsd_p .or. jc > jed_p) cycle

                           if (i < Atm(n)%neststruct%ind_update_h(ic,jc,1) .and. &
                                j < Atm(n)%neststruct%ind_update_h(ic,jc,2) ) then
                              Atm(n)%neststruct%ind_update_h(ic,jc,:) = (/i, j/)
                           end if

                        end do
                     end do
                  end if

                  deallocate(p_ind)

               end if

            endif
         endif
      end do
      
    ! Initialize restart functions
      call fv_restart_init()

!     if ( reset_eta ) then
!         do n=1, ntilesMe
!            call set_eta(npz, Atm(n)%ks, ptop, Atm(n)%ak, Atm(n)%bk)
!         enddo
!         if(is_master()) write(*,*) "Hybrid sigma-p coordinate has been reset"
!     endif

      if (ntilesMe > 1) call switch_current_Atm(Atm(1))
      if (ntilesMe > 1) call setup_pointers(Atm(1))

 end subroutine fv_init
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
         
 subroutine fv_end(Atm, grids_on_this_pe)

    type(fv_atmos_type), intent(inout) :: Atm(:)
    logical, intent(INOUT) :: grids_on_this_pe(:)

    integer :: n

    call timing_off('TOTAL')
    call timing_prt( gid )

    call fv_restart_end(Atm, grids_on_this_pe)
    call fv_io_exit()

  ! Free temporary memory from sw_core routines

  ! Deallocate
    call grid_utils_end

    do n = 1, ntilesMe
       call deallocate_fv_atmos_type(Atm(n))
    end do


 end subroutine fv_end
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!
!     run_setup :: initialize run from namelist
!
      subroutine run_setup(Atm, dt_atmos, grids_on_this_pe, p_split)
      type(fv_atmos_type), intent(inout), target :: Atm(:)
      real, intent(in)                   :: dt_atmos
      logical, intent(INOUT) :: grids_on_this_pe(:)
      integer, intent(INOUT) :: p_split

      character(len=80) :: filename, tracerName, errString, nested_grid_filename
      integer :: ios, ierr, f_unit, unit
      logical :: exists

      real :: dim0 = 180.           ! base dimension
      real :: dt0  = 1800.          ! base time step
      real :: ns0  = 5.             ! base nsplit for base dimension 
                                    ! For cubed sphere 5 is better
      !real :: umax = 350.           ! max wave speed for grid_type>3 ! Now defined above
      real :: dimx, dl, dp, dxmin, dymin, d_fac

      integer :: n0split
      integer :: n, nn, i

      integer :: pe_counter

!     local version of these variables to allow PGI compiler to compile
      character(len=128) :: res_latlon_dynamics = ''
      character(len=128) :: res_latlon_tracers  = ''
      character(len=80)  :: grid_name = ''
      character(len=120) :: grid_file = ''

      namelist /fv_grid_nml/ grid_name, grid_file
      namelist /fv_core_nml/npx, npy, ntiles, npz, npz_rst, layout, io_layout, ncnst, nwat,  &
                         use_logp, p_fac, a_imp, k_split, n_split, m_split, q_split, print_freq, do_schmidt,      &
                         hord_mt, hord_vt, hord_tm, hord_dp, hord_tr, shift_fac, stretch_fac, target_lat, target_lon, &
                         kord_mt, kord_wz, kord_tm, kord_tr, fv_debug, fv_land, nudge, do_sat_adj, do_f3d, &
                         external_ic, ncep_ic, fv_diag_ic, res_latlon_dynamics, res_latlon_tracers, &
                         scale_z, w_max, z_min, dddmp, d2_bg, d4_bg, vtdm4, d_ext, beta, non_ortho, n_sponge, &
                         warm_start, adjust_dry_mass, mountain, d_con, nord, convert_ke, use_old_omega, &
                         dry_mass, grid_type, do_Held_Suarez, do_reed_physics, reed_cond_only, &
                         consv_te, fill, filter_phys, fill_dp, fill_wz, consv_am, &
                         range_warn, dwind_2d, inline_q, z_tracer, reproduce_sum, adiabatic, do_vort_damp, no_dycore,   &
                         replace_w, tau, tau_h2o, rf_cutoff, nf_omega, hydrostatic, fv_sg_adj, breed_vortex_inline,  &
                         na_init, hybrid_z, Make_NH, n_zs_filter, nord_zs_filter, reset_eta,         &
                         pnats, dnats, a2b_ord, remap_t, p_ref, d2_bg_k1, d2_bg_k2,  &
                         c2l_ord, dx_const, dy_const, umax, deglat,      &
                         deglon_start, deglon_stop, deglat_start, deglat_stop, &
                         phys_hydrostatic, make_hybrid_z, old_divg_damp, add_noise, &
                         nested, twowaynest, parent_grid_num, parent_tile, &
                         refinement, nestbctype, nestupdate, nsponge, s_weight, &
                         ioffset, joffset, check_negative, nudge_ic

      namelist /test_case_nml/test_case,alpha

      pe_counter = mpp_root_pe()

! Make alpha = 0 the default:
      alpha = 0.
      test_case = 11   ! (USGS terrain)

      filename = "input.nml"

      inquire(file=filename,exist=exists)
      if (.not. exists) then  ! This will be replaced with fv_error wrapper
        if(is_master()) write(*,*) "file ",trim(filename)," doesn't exist" 
        call mpp_error(FATAL,'FV core terminating 1')
     endif

#ifdef INTERNAL_FILE_NML
!      rewind (f_unit)
   ! Read Main namelist
      read (input_nml_file,fv_grid_nml,iostat=ios)
      ierr = check_nml_error(ios,'fv_grid_nml')
   ! Read Test_Case namelist
      read (input_nml_file,test_case_nml,iostat=ios)
      ierr = check_nml_error(ios,'test_case_nml')
#else
      f_unit=open_namelist_file()
      rewind (f_unit)
   ! Read Main namelist
      read (f_unit,fv_grid_nml,iostat=ios)
      ierr = check_nml_error(ios,'fv_grid_nml')
      rewind (f_unit)
#endif

      unit = stdlog()
      write(unit, nml=fv_grid_nml)

      do n=1,size(Atm)

         call switch_current_Atm(Atm(n), .false.)
         call setup_pointers(Atm(n))
         Atm(n)%grid_number = n

#ifdef INTERNAL_FILE_NML
         if (size(Atm) > 1) then
            call mpp_error(FATAL, "Nesting not implemented with INTERNAL_FILE_NML")
         endif
   ! Read FVCORE namelist 
      read (input_nml_file,fv_core_nml,iostat=ios)
      ierr = check_nml_error(ios,'fv_core_nml')
   ! Read Test_Case namelist
      read (input_nml_file,test_case_nml,iostat=ios)
      ierr = check_nml_error(ios,'test_case_nml')
#else
      if (size(Atm) == 1) then
         f_unit = open_namelist_file()
      else if (n == 1) then
         f_unit = open_namelist_file('input.nml')
      else 
         write(nested_grid_filename,'(A10, I2.2, A4)') 'input_nest', n, '.nml'
         f_unit = open_namelist_file(nested_grid_filename)
      endif

   ! Read FVCORE namelist 
      read (f_unit,fv_core_nml,iostat=ios)
      ierr = check_nml_error(ios,'fv_core_nml')

   ! Read Test_Case namelist
      rewind (f_unit)
      read (f_unit,test_case_nml,iostat=ios)
      ierr = check_nml_error(ios,'test_case_nml')

      call close_file(f_unit)
#endif         
          if (len_trim(grid_file) /= 0) Atm(n)%flagstruct%grid_file = grid_file
          if (len_trim(grid_name) /= 0) Atm(n)%flagstruct%grid_name = grid_name
          if (len_trim(res_latlon_dynamics) /= 0) Atm(n)%flagstruct%res_latlon_dynamics = res_latlon_dynamics
          if (len_trim(res_latlon_tracers)  /= 0) Atm(n)%flagstruct%res_latlon_tracers = res_latlon_tracers

         write(unit, nml=fv_core_nml)
         write(unit, nml=test_case_nml)

         !*** single tile for Cartesian grids
         if (grid_type>3) then
            ntiles=1
            non_ortho = .false.
            nf_omega = 0
         endif

         if (.not. nested) Atm(n)%neststruct%npx_global = npx

         ! Define n_split if not in namelist
         if (ntiles==6) then
            dimx = 4.0*(npx-1)
         if ( hydrostatic ) then
            if ( npx >= 120 ) ns0 = 6
         else
            if ( npx <= 45 ) then
                 ns0 = 6
            elseif ( npx <=90 ) then
                 ns0 = 7
            else
                 ns0 = 8
            endif
         endif
      else
         dimx = max ( npx, 2*(npy-1) )
      endif
          
      if (grid_type < 4) then
         n0split = nint ( ns0*abs(dt_atmos)*dimx/(dt0*dim0) + 0.49 )
      elseif (grid_type == 4 .or. grid_type == 7) then
         n0split = nint ( 2.*umax*dt_atmos/sqrt(dx_const**2 + dy_const**2) + 0.49 )
      elseif (grid_type == 5 .or. grid_type == 6) then
         if (grid_type == 6) then
            deglon_start = 0.; deglon_stop  = 360.
         endif
         dl = (deglon_stop-deglon_start)*pi/(180.*(npx-1))
         dp = (deglat_stop-deglat_start)*pi/(180.*(npy-1))

         dxmin=dl*radius*min(cos(deglat_start*pi/180.-ng*dp),   &
                             cos(deglat_stop *pi/180.+ng*dp))
         dymin=dp*radius
         n0split = nint ( 2.*umax*dt_atmos/sqrt(dxmin**2 + dymin**2) + 0.49 )
      endif
      n0split = max ( 1, n0split )

      if ( n_split == 0 ) then
           n_split = nint( real(n0split)/real(k_split*p_split) * stretch_fac + 0.5 )
           if(is_master()) write(*,*) 'For k_split (remapping)=', k_split
           if(is_master()) write(*,198) 'n_split is set to ', n_split, ' for resolution-dt=',npx,npy,ntiles,dt_atmos
      else
          if(is_master()) write(*,199) 'Using n_split from the namelist: ', n_split
      endif
      if (is_master() .and. n == 1 .and. p_split > 1) then
         write(*,199) 'Using p_split = ', p_split
      endif

      if (Atm(n)%neststruct%nested) then
         do i=1,n-1
            if (Atm(i)%grid_number == parent_grid_num) then
               Atm(n)%parent_grid => Atm(i)
               exit
            end if
         end do
         if (.not. associated(Atm(n)%parent_grid)) then
            write(errstring,'(2(A,I3))')  "Could not find parent grid #", parent_grid_num, ' for grid #', n
            call mpp_error(FATAL, errstring)
         end if

         !Note that if a gnomonic grid has a parent it is a NESTED gnomonic grid and therefore only has one tile
         if ( Atm(n)%parent_grid%flagstruct%grid_type < 3 .and. &
              .not. associated(Atm(n)%parent_grid%parent_grid)) then
            if (parent_tile > 6 .or. parent_tile < 1) then
               call mpp_error(FATAL, 'parent tile must be between 1 and 6 if the parent is a cubed-sphere grid')
            end if
         else
            if (parent_tile /= 1) then
               call mpp_error(FATAL, 'parent tile must be 1 if the parent is not a cubed-sphere grid')
            end if
         end if

         if ( refinement < 1 ) call mpp_error(FATAL, 'grid refinement must be positive')

         if (nestupdate == 1 .or. nestupdate == 2) then

            if (mod(npx-1,refinement) /= 0 .or. mod(npy-1,refinement) /= 0) then
               call mpp_error(WARNING, 'npx-1 or npy-1 is not evenly divisible by the refinement ratio; averaging update cannot be mass-conservative.')
            end if

         end if

         if ( consv_te > 0.) then
            call mpp_error(FATAL, 'The global energy fixer cannot be used on a nested grid. consv_te must be set to 0.')
         end if

         Atm(n)%neststruct%refinement_of_global = Atm(n)%neststruct%refinement * Atm(n)%parent_grid%neststruct%refinement_of_global
         max_refinement_of_global = max(Atm(n)%neststruct%refinement_of_global,max_refinement_of_global)
         Atm(n)%neststruct%npx_global = Atm(n)%neststruct%refinement * Atm(n)%parent_grid%neststruct%npx_global

      else
         Atm(n)%neststruct%ioffset                = -999
         Atm(n)%neststruct%joffset                = -999   
         Atm(n)%neststruct%parent_tile            = -1      
         Atm(n)%neststruct%refinement             = -1
      end if

      if (Atm(n)%neststruct%nested) then
         if (Atm(n)%flagstruct%grid_type >= 4 .and. Atm(n)%parent_grid%flagstruct%grid_type >= 4) then
            Atm(n)%flagstruct%dx_const = Atm(n)%parent_grid%flagstruct%dx_const / real(Atm(n)%neststruct%refinement)
            Atm(n)%flagstruct%dy_const = Atm(n)%parent_grid%flagstruct%dy_const / real(Atm(n)%neststruct%refinement)
         end if
      end if


!----------------------------------------
! Adjust divergence damping coefficients:
!----------------------------------------
!      d_fac = real(n0split)/real(n_split)
!      dddmp = dddmp * d_fac
!      d2_bg = d2_bg * d_fac
!      d4_bg = d4_bg * d_fac
!      d_ext = d_ext * d_fac
!      vtdm4 = vtdm4 * d_fac
      if (old_divg_damp) then
        if (is_master()) write(*,*) " fv_control: using original values for divergence damping "
        d2_bg_k1 = 6.         ! factor for d2_bg (k=1)  - default(4.)
        d2_bg_k2 = 4.         ! factor for d2_bg (k=2)  - default(2.)
        d2_divg_max_k1 = 0.02 ! d2_divg max value (k=1) - default(0.05)
        d2_divg_max_k2 = 0.01 ! d2_divg max value (k=2) - default(0.02)
        damp_k_k1 = 0.        ! damp_k value (k=1)      - default(0.05)
        damp_k_k2 = 0.        ! damp_k value (k=2)      - default(0.025)
      elseif (n_sponge == 0 ) then
        if ( d2_bg_k1 > 1. ) d2_bg_k1 = 0.20
        if ( d2_bg_k2 > 1. ) d2_bg_k2 = 0.015
      endif

!     if ( beta < 1.e-5 ) beta = 0.   ! beta < 0 is used for non-hydrostatic "one_grad_p"

      if ( .not.hydrostatic ) then
         if ( m_split==0 ) then
              m_split = 1. + abs(dt_atmos)/real(k_split*n_split*p_split)
         endif
         if(is_master()) write(*,198) 'm_split is set to ', m_split
         if(is_master()) then
            write(*,*) 'Off center implicit scheme param=', a_imp
            write(*,*) ' p_fac=', p_fac
         endif
      endif

      if(is_master()) then
         write(*,199) 'Using n_sponge : ', n_sponge
         write(*,197) 'Using non_ortho : ', non_ortho
      endif

 197  format(A,l7)
 198  format(A,i2.2,A,i4.4,'x',i4.4,'x',i1.1,'-',f9.3)
 199  format(A,i2.2)

      if (.not. nested) alpha = alpha*pi


      allocate(Atm(n)%neststruct%child_grids(size(Atm)))
      Atm(N)%neststruct%child_grids = .false.

      !Broadcast data

      !Check layout

   enddo

   !Set pelists
   do n=1,size(Atm)
      if (ANY(Atm(n)%pelist == gid)) then
            call mpp_set_current_pelist(Atm(n)%pelist)
            call mpp_get_current_pelist(Atm(n)%pelist, commID=commID)
            call mp_start(commID)
      endif

      if (Atm(n)%neststruct%nested) then
         Atm(n)%neststruct%parent_proc = ANY(Atm(n)%parent_grid%pelist == gid)
         Atm(n)%neststruct%child_proc = ANY(Atm(n)%pelist == gid)
      endif
   enddo

   do n=1,size(Atm)
      
      call switch_current_Atm(Atm(n),.false.)
      call setup_pointers(Atm(n))
      !! CLEANUP: WARNING not sure what changes to domain_decomp may cause
      call domain_decomp(npx,npy,ntiles,grid_type,nested,Atm(n),layout,io_layout)
   enddo

   !!! CLEANUP: This sets the pelist to ALL, which is also
   !!! required for the define_nest_domains step in the next loop.
   !!! Later the pelist must be reset to the 'local' pelist.
   call broadcast_domains(Atm)

   do n=1,size(Atm)
      call switch_current_Atm(Atm(n))
      call setup_pointers(Atm(n))

      if (nested) then
         if (mod(npx-1 , refinement) /= 0 .or. mod(npy-1, refinement) /= 0) &
              call mpp_error(FATAL, 'npx or npy not an even refinement of its coarse grid.')
         
         !Pelist needs to be set to ALL (which should have been done
         !in broadcast_domains) to get this to work
         call mpp_define_nest_domains(Atm(n)%neststruct%nest_domain, Atm(n)%domain, Atm(parent_grid_num)%domain, &
              7, parent_tile, &
              1, npx-1, 1, npy-1,                  & !Grid cells, not points
              ioffset, ioffset + (npx-1)/refinement - 1, &
              joffset, joffset + (npy-1)/refinement - 1,         &
              (/ (i,i=0,mpp_npes()-1)  /), 0, name="nest_domain") !What pelist to use?

         Atm(parent_grid_num)%neststruct%child_grids(n) = .true.

         if (Atm(n)%neststruct%nestbctype > 1) then

            call mpp_error(FATAL, 'nestbctype > 1 not yet implemented')

            !This check is due to a bug which has not yet been identified. Beware.
!            if (Atm(n)%parent_grid%flagstruct%hord_tr == 7) &
!                 call mpp_error(FATAL, "Flux-form nested  BCs (nestbctype > 1) should not use hord_tr == 7 (on parent grid), since there is no guarantee of tracer mass conservation with this option.")

!!$            if (Atm(n)%flagstruct%q_split > 0 .and. Atm(n)%parent_grid%flagstruct%q_split > 0) then
!!$               if (mod(Atm(n)%flagstruct%q_split,Atm(n)%parent_grid%flagstruct%q_split) /= 0) call mpp_error(FATAL, &
!!$                    "Flux-form nested BCs (nestbctype > 1) require q_split on the nested grid to be evenly divisible by that on the coarse grid.")
!!$            endif
!!$            if (mod((Atm(n)%npx-1),Atm(n)%neststruct%refinement) /= 0 .or. mod((Atm(n)%npy-1),Atm(n)%neststruct%refinement) /= 0) call mpp_error(FATAL, &
!!$                 "Flux-form nested BCs (nestbctype > 1) requires npx and npy to be one more than a multiple of the refinement ratio.")
!!$            Atm(n)%parent_grid%neststruct%do_flux_BCs = .true.
!!$            if (Atm(n)%neststruct%nestbctype == 3 .or. Atm(n)%neststruct%nestbctype == 4) Atm(n)%parent_grid%neststruct%do_2way_flux_BCs = .true.
         endif

      end if

      do nn=1,size(Atm)
         if (n == 1) allocate(Atm(nn)%neststruct%nest_domain_all(size(Atm)))
         Atm(nn)%neststruct%nest_domain_all(n) = Atm(n)%neststruct%nest_domain
      enddo

   end do

   do n=1,size(Atm)
      if (ANY(Atm(n)%pelist == gid)) then
         call mpp_set_current_pelist(Atm(n)%pelist)
      endif
   enddo

  end subroutine run_setup

  subroutine init_nesting(Atm, grids_on_this_pe, p_split)
    
    type(fv_atmos_type), intent(inout), allocatable :: Atm(:)
   logical, allocatable, intent(INOUT) :: grids_on_this_pe(:)
    integer, intent(INOUT) :: p_split
    character(100) :: pe_list_name
    integer :: nest_pes(100)
    integer :: n, npes, ntiles, pecounter, i
    integer, allocatable :: pelist(:)
    integer :: f_unit, ios, ierr

    !This is an OPTIONAL namelist, that needs to be read before everything else
    namelist /nest_nml/ ngrids, ntiles, nest_pes, p_split

    call mp_assign_gid

    nest_pes = 0
    ntiles = -999

#ifdef INTERNAL_FILE_NML
      read (input_nml_file,nest_nml,iostat=ios)
      ierr = check_nml_error(ios,'nest_nml')
#else
      f_unit=open_namelist_file()
      rewind (f_unit)
      read (f_unit,nest_nml,iostat=ios)
      ierr = check_nml_error(ios,'nest_nml')
      call close_file(f_unit)
#endif

      if (ntiles /= -999) ngrids = ntiles
      if (ngrids > 10) call mpp_error(FATAL, "More than 10 nested grids not supported")

      allocate(Atm(ngrids))
    
      allocate(grids_on_this_pe(ngrids))
      grids_on_this_pe = .false. !initialization

      npes = mpp_npes()

      ! Need to get a global pelist to send data around later?
      allocate( pelist_all(npes) )
      pelist_all = (/ (i,i=0,npes-1) /)
      pelist_all = pelist_all + mpp_root_pe()

      if (ngrids == 1) then

         !Set up the single pelist
         allocate(Atm(1)%pelist(npes))
         Atm(1)%pelist = (/(i, i=0, npes-1)/)
         Atm(1)%pelist = Atm(1)%pelist + mpp_root_pe()
         call mpp_declare_pelist(Atm(1)%pelist)
         call mpp_set_current_pelist(Atm(1)%pelist)
         !Now set in domain_decomp
         !masterproc = Atm(1)%pelist(1)
         call setup_master(Atm(1)%pelist)
         grids_on_this_pe(1) = .true.
         Atm(1)%npes_this_grid = npes

      else

         pecounter = mpp_root_pe()
         do n=1,ngrids
            if (n == 1) then
               pe_list_name = ''
            else
               write(pe_list_name,'(A4, I2.2)') 'nest', n
            endif

            if (nest_pes(n) == 0) then
               if (n < ngrids) call mpp_error(FATAL, 'Only nest_pes(ngrids) in nest_nml can be zero; preceeding values must be nonzero.')
               allocate(Atm(n)%pelist(npes-pecounter))
               Atm(n)%pelist = (/(i, i=pecounter, npes-1)/)
               if (n > 1) then
                  call mpp_declare_pelist(Atm(n)%pelist, trim(pe_list_name))
                  !Make sure nested-grid input file exists
                  if (.not. file_exist('input_'//trim(pe_list_name)//'.nml')) then
                     call mpp_error(FATAL, "Could not find nested grid namelist input_"//trim(pe_list_name)//".nml")
                  endif
               endif
               exit
            else
               allocate(Atm(n)%pelist(nest_pes(n)))
               Atm(n)%pelist = (/ (i, i=pecounter, pecounter+nest_pes(n)-1) /)
               if (Atm(n)%pelist(nest_pes(n)) >= npes) then
                  call mpp_error(FATAL, 'PEs assigned by nest_pes in nest_nml exceeds number of available PEs.')
               endif

               call mpp_declare_pelist(Atm(n)%pelist, trim(pe_list_name))
               !Make sure nested-grid input file exists
               if (n > 1) then
                  if (.not. file_exist('input_'//trim(pe_list_name)//'.nml')) then
                     call mpp_error(FATAL, "Could not find nested grid namelist input_"//trim(pe_list_name)//".nml")
                  endif
               endif
               pecounter = pecounter+nest_pes(n)
            endif
         enddo

         !Set pelists
         do n=1,ngrids
            Atm(n)%npes_this_grid = size(Atm(n)%pelist)
            if (ANY(gid == Atm(n)%pelist)) then
                  call mpp_set_current_pelist(Atm(n)%pelist)
                  !now set in domain_decomp
                  !masterproc = Atm(n)%pelist(1)
                  call setup_master(Atm(n)%pelist)
                  grids_on_this_pe(n) = .true.
#if defined (INTERNAL_FILE_NML)
                  if (n > 1) call read_input_nml
#else
                  !Namelist file read in fv_control.F90
#endif
                  exit
               endif
         enddo

         if (pecounter /= npes) then
            call mpp_error(FATAL, 'nest_pes in nest_nml does not assign all of the available PEs.')
         endif
      endif

      !Layout is checked later, in fv_control

  end subroutine init_nesting

  subroutine setup_pointers(Atm)

    type(fv_atmos_type), intent(INOUT), target :: Atm

    !This routine associates the MODULE flag pointers with the ARRAY flag variables for the grid active on THIS pe so the flags can be read in from the namelist.

     res_latlon_dynamics           => Atm%flagstruct%res_latlon_dynamics
     res_latlon_tracers            => Atm%flagstruct%res_latlon_tracers

     grid_type                     => Atm%flagstruct%grid_type
     grid_name                     => Atm%flagstruct%grid_name
     grid_file                     => Atm%flagstruct%grid_file
     hord_mt                       => Atm%flagstruct%hord_mt
     kord_mt                       => Atm%flagstruct%kord_mt
     kord_wz                       => Atm%flagstruct%kord_wz
     hord_vt                       => Atm%flagstruct%hord_vt
     hord_tm                       => Atm%flagstruct%hord_tm
     hord_dp                       => Atm%flagstruct%hord_dp
     kord_tm                       => Atm%flagstruct%kord_tm
     hord_tr                       => Atm%flagstruct%hord_tr
     kord_tr                       => Atm%flagstruct%kord_tr
     scale_z                       => Atm%flagstruct%scale_z
     w_max                         => Atm%flagstruct%w_max
     z_min                         => Atm%flagstruct%z_min
     nord                          => Atm%flagstruct%nord
     dddmp                         => Atm%flagstruct%dddmp
     d2_bg                         => Atm%flagstruct%d2_bg
     d4_bg                         => Atm%flagstruct%d4_bg
     vtdm4                         => Atm%flagstruct%vtdm4
     d2_bg_k1                      => Atm%flagstruct%d2_bg_k1
     d2_bg_k2                      => Atm%flagstruct%d2_bg_k2
     d2_divg_max_k1                => Atm%flagstruct%d2_divg_max_k1
     d2_divg_max_k2                => Atm%flagstruct%d2_divg_max_k2
     damp_k_k1                     => Atm%flagstruct%damp_k_k1
     damp_k_k2                     => Atm%flagstruct%damp_k_k2
     n_zs_filter                   => Atm%flagstruct%n_zs_filter
     nord_zs_filter                => Atm%flagstruct%nord_zs_filter
     consv_am                      => Atm%flagstruct%consv_am
     do_sat_adj                    => Atm%flagstruct%do_sat_adj
     do_f3d                        => Atm%flagstruct%do_f3d
     no_dycore                     => Atm%flagstruct%no_dycore
     replace_w                     => Atm%flagstruct%replace_w
     convert_ke                    => Atm%flagstruct%convert_ke
     do_vort_damp                  => Atm%flagstruct%do_vort_damp
     use_old_omega                 => Atm%flagstruct%use_old_omega
     beta                          => Atm%flagstruct%beta
     n_sponge                      => Atm%flagstruct%n_sponge
     d_ext                         => Atm%flagstruct%d_ext
     nwat                          => Atm%flagstruct%nwat
     use_logp                      => Atm%flagstruct%use_logp
     warm_start                    => Atm%flagstruct%warm_start
     inline_q                      => Atm%flagstruct%inline_q
     shift_fac                     => Atm%flagstruct%shift_fac
     do_schmidt                    => Atm%flagstruct%do_schmidt
     stretch_fac                   => Atm%flagstruct%stretch_fac
     target_lat                    => Atm%flagstruct%target_lat
     target_lon                    => Atm%flagstruct%target_lon
     reset_eta                     => Atm%flagstruct%reset_eta
     p_fac                         => Atm%flagstruct%p_fac
     a_imp                         => Atm%flagstruct%a_imp
     n_split                       => Atm%flagstruct%n_split
     m_split                       => Atm%flagstruct%m_split
     k_split                       => Atm%flagstruct%k_split
     use_logp                      => Atm%flagstruct%use_logp
     q_split                       => Atm%flagstruct%q_split
     print_freq                    => Atm%flagstruct%print_freq
     npx                           => Atm%flagstruct%npx
     npy                           => Atm%flagstruct%npy
     npz                           => Atm%flagstruct%npz
     npz_rst                       => Atm%flagstruct%npz_rst
     ncnst                         => Atm%flagstruct%ncnst
     pnats                         => Atm%flagstruct%pnats
     dnats                         => Atm%flagstruct%dnats
     ntiles                        => Atm%flagstruct%ntiles
     nf_omega                      => Atm%flagstruct%nf_omega
     fv_sg_adj                     => Atm%flagstruct%fv_sg_adj
     na_init                       => Atm%flagstruct%na_init
     p_ref                         => Atm%flagstruct%p_ref
     dry_mass                      => Atm%flagstruct%dry_mass
     nt_prog                       => Atm%flagstruct%nt_prog
     nt_phys                       => Atm%flagstruct%nt_phys
     tau_h2o                       => Atm%flagstruct%tau_h2o
     d_con                         => Atm%flagstruct%d_con
     consv_te                      => Atm%flagstruct%consv_te
     tau                           => Atm%flagstruct%tau
     rf_cutoff                     => Atm%flagstruct%rf_cutoff
     filter_phys                   => Atm%flagstruct%filter_phys
     dwind_2d                      => Atm%flagstruct%dwind_2d
     breed_vortex_inline           => Atm%flagstruct%breed_vortex_inline
     range_warn                    => Atm%flagstruct%range_warn
     fill                          => Atm%flagstruct%fill
     fill_dp                       => Atm%flagstruct%fill_dp
     fill_wz                       => Atm%flagstruct%fill_wz
     check_negative                => Atm%flagstruct%check_negative
     non_ortho                     => Atm%flagstruct%non_ortho
     adiabatic                     => Atm%flagstruct%adiabatic
     moist_phys                    => Atm%flagstruct%moist_phys
     do_Held_Suarez                => Atm%flagstruct%do_Held_Suarez
     do_reed_physics               => Atm%flagstruct%do_reed_physics
     reed_cond_only                => Atm%flagstruct%reed_cond_only
     reproduce_sum                 => Atm%flagstruct%reproduce_sum
     adjust_dry_mass               => Atm%flagstruct%adjust_dry_mass
     fv_debug                      => Atm%flagstruct%fv_debug
     srf_init                      => Atm%flagstruct%srf_init
     mountain                      => Atm%flagstruct%mountain
     remap_t                       => Atm%flagstruct%remap_t
     z_tracer                      => Atm%flagstruct%z_tracer
     old_divg_damp                 => Atm%flagstruct%old_divg_damp
     fv_land                       => Atm%flagstruct%fv_land
     nudge                         => Atm%flagstruct%nudge
     nudge_ic                      => Atm%flagstruct%nudge_ic
     ncep_ic                       => Atm%flagstruct%ncep_ic
     fv_diag_ic                    => Atm%flagstruct%fv_diag_ic
     external_ic                   => Atm%flagstruct%external_ic

     hydrostatic                   => Atm%flagstruct%hydrostatic
     phys_hydrostatic              => Atm%flagstruct%phys_hydrostatic
     hybrid_z                      => Atm%flagstruct%hybrid_z
     Make_NH                       => Atm%flagstruct%Make_NH
     make_hybrid_z                 => Atm%flagstruct%make_hybrid_z
     add_noise                     => Atm%flagstruct%add_noise
     a2b_ord                       => Atm%flagstruct%a2b_ord
     c2l_ord                       => Atm%flagstruct%c2l_ord
     ndims                         => Atm%flagstruct%ndims
    
     dx_const                      => Atm%flagstruct%dx_const
     dy_const                      => Atm%flagstruct%dy_const
     deglon_start                  => Atm%flagstruct%deglon_start
     deglon_stop                   => Atm%flagstruct%deglon_stop
     deglat_start                  => Atm%flagstruct%deglat_start
     deglat_stop                   => Atm%flagstruct%deglat_stop

     deglat                        => Atm%flagstruct%deglat

     nested                        => Atm%neststruct%nested
     twowaynest                    => Atm%neststruct%twowaynest
     parent_tile                   => Atm%neststruct%parent_tile
     refinement                    => Atm%neststruct%refinement
     nestbctype                    => Atm%neststruct%nestbctype
     nestupdate                    => Atm%neststruct%nestupdate
     nsponge                       => Atm%neststruct%nsponge
     s_weight                      => Atm%neststruct%s_weight
     ioffset                       => Atm%neststruct%ioffset
     joffset                       => Atm%neststruct%joffset

     layout                        => Atm%layout
     io_layout                     => Atm%io_layout
  end subroutine setup_pointers

       
end module fv_control_mod
