! $Id: fv_control.F90,v 17.0.2.9.2.20.2.3 2012/06/12 18:03:35 Lucas.Harris Exp $
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
                                  input_nml_file, get_unit, WARNING, read_ascii_file, INPUT_STR_LENGTH
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
   use fv_grid_utils_mod,   only: grid_utils_init, grid_utils_end, ptop_min, deglat, &
                                 da_min_c, da_min, ptop, ks
   use fv_eta_mod,          only: set_eta
   use fv_grid_tools_mod,   only: init_grid, cosa, sina, area, area_c, dx, dy, dxa, dya, &
                                 dxc, dyc, grid_type, dx_const, dy_const,                         &
                                 deglon_start, deglon_stop, deglat_start, deglat_stop, &
                                 read_grid
   use fv_mp_mod,           only: mp_start, mp_assign_gid, domain_decomp, domain, &
                                 ng, tile, gid, switch_current_Atm, &
                                 npes_x, npes_y, npes_all => npes, io_domain_layout
   use test_cases_mod,      only: test_case, alpha
   use fv_timing_mod,       only: timing_on, timing_off, timing_init, timing_prt
   use fv_current_grid_mod, only: grid_name, grid_file, hord_mt,&
                                  kord_mt, kord_wz, hord_vt, hord_tm, hord_dp, kord_tm,&
                                  hord_tr, kord_tr, scale_z, w_max, z_min, nord,&
                                  dddmp, d2_bg, d4_bg, vtdm4, d2_bg_k1, d2_bg_k2,&
                                  d2_divg_max_k1, d2_divg_max_k2, damp_k_k1, damp_k_k2,&
                                  n_zs_filter, nord_zs_filter, no_dycore, replace_w,&
                                  convert_ke, do_vort_damp, use_old_omega, beta, n_sponge,&
                                  d_ext, nwat, warm_start, inline_q, n_sponge, d_ext, nwat,&
                                  warm_start, inline_q, shift_fac, do_schmidt, stretch_fac,&
                                  target_lat, target_lon, reset_eta, n_split,&
                                  m_split, k_split, q_split, print_freq, npx, npy, npz,&
                                  npz_rst, ncnst, pnats, ntiles,&
                                  nf_omega, fv_sg_adj, na_init, p_ref, dry_mass, p_ref, dry_mass, nt_prog, nt_phys,&
                                  tau_h2o, d_con, consv_te, tau, rf_center,&
                                  tq_filter, filter_phys, dwind_2d, breed_vortex_inline,&
                                  range_warn, fill, fill_dp, fill_wz, non_ortho, adiabatic,&
                                  moist_phys, do_Held_Suarez, reproduce_sum, adjust_dry_mass,&
                                  fv_debug, srf_init, mountain, remap_t,&
                                  z_tracer, old_divg_damp, master, fv_land,&
                                  nudge, ncep_ic, fv_diag_ic, external_ic,&
                                  res_latlon_dynamics, res_latlon_tracers, hydrostatic,&
                                  phys_hydrostatic, hybrid_z, Make_NH, make_hybrid_z,&
                                  a2b_ord, c2l_ord, ndims, &
                                  nested, twowaynest, parent_grid, parent_tile, &
                                  refinement, nestbctype, nestupdate, nsponge, s_weight, ioffset, joffset, &
                                  isd, ied, jsd, jed, isc, iec, jsc, jec, is, ie, js, je, layout, io_layout, a_imp, p_fac
   use fv_grid_tools_mod,   only: grid, agrid
#ifdef MARS_GCM
   use fv_current_grid_mod, only: reference_sfc_pres, sponge_damp
#endif
   use fv_current_grid_mod, only: switch_current_grid_pointers
   use mpp_domains_mod,     only: domain2D
   use fv_mp_mod,           only: concurrent, grids_on_this_pe, masterproc, broadcast_domains, npes, mp_barrier
   use mpp_domains_mod,     only: mpp_define_nest_domains, nest_domain_type, mpp_get_global_domain
   use mpp_domains_mod,     only : mpp_get_C2F_index, mpp_get_F2C_index, mpp_broadcast_domain
   use mpp_domains_mod,     only: CENTER, CORNER, NORTH, EAST, WEST, SOUTH
   use fv_current_grid_mod, only: nest_domain, child_grids, pelist, grid_global, npes_this_grid
   use mpp_mod,             only: mpp_send, mpp_sync, mpp_transmit, mpp_set_current_pelist, mpp_declare_pelist, mpp_root_pe, mpp_recv, mpp_sync_self

   implicit none
   private

!-----------------------------------------------------------------------
! Grid descriptor file setup
!-----------------------------------------------------------------------
!------------------------------------------
! Model Domain parameters
! See fv_arrays.F90 for descriptions
!------------------------------------------
   integer :: ntilesMe                ! Number of tiles on this process =1 for now

   real    :: too_big  = 1.E35
   public :: shift_fac, do_schmidt, stretch_fac, target_lat, target_lon 
   public :: npx,npy,npz, npz_rst, ntiles, ncnst, pnats, nwat
   public :: hord_mt, hord_vt, kord_mt, kord_wz, hord_tm, hord_dp, kord_tm, hord_tr, kord_tr
   public :: nord, fill_dp, fill_wz, inline_q, breed_vortex_inline, dwind_2d, filter_phys, tq_filter 
   public :: p_fac, a_imp, k_split, n_split, m_split, q_split, master
   public :: scale_z, w_max, z_min, dddmp, d2_bg, d4_bg, d_ext, vtdm4, beta
   public :: n_sponge, p_ref, mountain
   public :: remap_t,  z_tracer, fv_debug
   public :: external_ic, ncep_ic, fv_diag_ic, res_latlon_dynamics, res_latlon_tracers, fv_land
   public :: na_init, fv_sg_adj, tau, tau_h2o, rf_center, d_con
   public :: fv_init, fv_end
   public :: domain
   public :: adiabatic, nf_omega, moist_phys, range_warn, do_vort_damp
   public :: no_dycore, replace_w, convert_ke, use_old_omega
   public :: hydrostatic, phys_hydrostatic,  hybrid_z, a2b_ord
   public :: nt_prog, nt_phys
   public :: d2_bg_k1, d2_bg_k2, d2_divg_max_k1, d2_divg_max_k2, damp_k_k1, damp_k_k2

#ifdef MARS_GCM
   public :: reference_sfc_pres, sponge_damp
#endif MARS

   integer, allocatable :: pelist_all(:)
   integer :: commID, max_refinement_of_global = 1.

!---- version number -----
   character(len=128) :: version = '$Id: fv_control.F90,v 17.0.2.9.2.20.2.3 2012/06/12 18:03:35 Lucas.Harris Exp $'
   character(len=128) :: tagname = '$Name: siena_201211 $'

 contains

!-------------------------------------------------------------------------------
         
 subroutine fv_init(Atm, dt_atmos)

   type(fv_atmos_type), intent(inout) :: Atm(:)
   real,                intent(in)    :: dt_atmos

   integer :: i, j, k, n, p
   real :: sdt

! tracers
   integer :: num_family          ! output of register_tracers

   integer :: isd_p, ied_p, jsd_p, jed_p, isg, ieg, jsg, jeg
   integer :: ic, jc
   integer, allocatable :: p_ind(:,:,:)

   !This call is needed to set up the pointers for fv_current_grid, even for a single-grid run
   call switch_current_Atm(Atm(1), .false.)

! Start up MPI

   call mp_assign_gid

    ! Initialize timing routines
      call timing_init
      call timing_on('TOTAL')

    ! Setup the run from namelist 
      ntilesMe = size(Atm(:)) !Full number of Atm arrays; one less than number of grids, if multiple grids

      allocate(grids_on_this_pe(size(Atm)))

      call run_setup(Atm,dt_atmos)   ! initializes domain_decomp

      do n=1,ntilesMe

         !In a single-grid run this will still be needed to correctly set the domain
         call switch_current_Atm(Atm(n))
         
         target_lon = target_lon * pi/180.
         target_lat = target_lat * pi/180.

!--------------------------------------------------
! override number of tracers by reading field_table
!--------------------------------------------------

         !not sure if this works with multiple grids
         call tm_register_tracers (MODEL_ATMOS, ncnst, nt_prog, pnats, num_family)
         if(master) write(*,*) 'ncnst=', ncnst,' num_prog=',nt_prog,' pnats=',pnats,' num_family=',num_family

         
         if (gid == 0) print*, ''
         if (grids_on_this_pe(n)) then
            call allocate_fv_atmos_type(Atm(n), isd, ied, jsd, jed, isc, iec, jsc, jec, &
                 npx, npy, npz, ndims, ncnst, ng, .false., concurrent, grids_on_this_pe(n))

            if (grids_on_this_pe(n)) then

               call switch_current_Atm(Atm(n))

               !if (master) print*, isc, iec, jsc, jec
               if ( (iec-isc+1).lt.4 .or. (jec-jsc+1).lt.4 ) &
                    call mpp_error(FATAL,'Domain Decomposition:  Cubed Sphere compute domain has a &
                    &minium requirement of 4 points in X and Y, respectively')

            endif

            ! Read Grid from GRID_FILE and setup grid descriptors
            ! needs modification for multiple tiles
            if(grid_type <0 .AND. trim(grid_file) == 'INPUT/grid_spec.nc') then
               call read_grid(Atm(n), grid_name, grid_file, npx, npy, npz, ndims, ntiles, ng)
            else
               call init_grid(Atm(n), grid_name, grid_file, npx, npy, npz, ndims, ntiles, ng)
            endif

            ! Initialize the SW (2D) part of the model
            call grid_utils_init(Atm(n), npx, npy, npz, grid, agrid,   &
                 area, area_c, cosa, sina, dx, dy, dxa, dya, dxc, dyc, non_ortho,   &
                 grid_type, c2l_ord)

            if ( master ) then
               sdt =  dt_atmos/real(n_split*k_split)
               write(*,*) ' '
               write(*,*) 'Divergence damping Coefficients * 1.E6:'
               write(*,*) 'For small dt=', sdt
               write(*,*) 'External mode del-2 (m**2/s)=',  d_ext*da_min_c     /sdt*1.E-6
               write(*,*) 'Internal mode del-2 SMAG dimensionless coeff=',  dddmp
               write(*,*) 'Internal mode del-2 background diff=', d2_bg

               if (nord==1) write(*,*) 'Internal mode del-4 background diff=', d4_bg
               if (nord==2) write(*,*) 'Internal mode del-6 background diff=', d4_bg
               if (nord==3) write(*,*) 'Internal mode del-8 background diff=', d4_bg

               write(*,*) 'Vorticity del-4 (m**4/s)=', (vtdm4*da_min)**2/sdt*1.E-6
               write(*,*) ' '
            endif


            Atm(n)%ts   = 0.
            Atm(n)%phis = too_big
            ! The following statements are to prevent the phatom corner regions from
            ! growing instability
            Atm(n)%u  = 0.
            Atm(n)%v  = 0.
            Atm(n)%ua = too_big
            Atm(n)%va = too_big

         else !this grid is NOT defined on this pe

            !Allocate dummy arrays
            call allocate_fv_atmos_type(Atm(n),  isd, ied, jsd, jed, isc, iec, jsc, jec, &
                 npx, npy, npz, ndims, ncnst, ng, .true., concurrent, .false.)

            !Need to SEND grid_global to any child grids; this is received in setup_aligned_nest in fv_grid_tools
            if (nested) then

               call mpp_get_global_domain( parent_grid%domain, &
                    isg, ieg, jsg, jeg)

               if (gid == parent_grid%pelist(1)) then
                  do p=1,size(pelist)
                     call mpp_send(parent_grid%grid_global(isg-ng:ieg+1+ng,jsg-ng:jeg+1+ng,1:2,parent_tile), &
                          size(parent_grid%grid_global(isg-ng:ieg+1+ng,jsg-ng:jeg+1+ng,1:2,parent_tile)), &
                          pelist(p)) !send to p_ind in setup_aligned_nest
                     call mpp_sync_self
                  enddo
               endif

               if (twowaynest) then

                  !Also need to receive P_IND from child grids to set up ind_update_h
                  isd_p = Atm(n)%parent_grid%isd
                  ied_p = Atm(n)%parent_grid%ied
                  jsd_p = Atm(n)%parent_grid%jsd
                  jed_p = Atm(n)%parent_grid%jed
                  !               allocate(Atm(n)%ind_update_h(isd_p:ied_p+1,jsd_p:jed_p+1,2))
                  allocate(p_ind(1-ng:npx+ng,1-ng:npy+ng,1:2))
                  call mpp_recv(p_ind,size(p_ind),pelist(1)) !receiving from p_ind setup_aligned_grids 
                  call mpp_sync

                  Atm(n)%ind_update_h = 1000000

                  if (Atm(n)%parent_grid%tile == Atm(n)%parent_tile) then
                     do j=1,npy
                        do i=1,npx

                           ic = p_ind(i,j,1)
                           jc = p_ind(i,j,2)

                           if (ic < isd_p .or. ic > ied_p .or. jc < jsd_p .or. jc > jed_p) cycle

                           if (i < Atm(n)%ind_update_h(ic,jc,1) .and. &
                                j < Atm(n)%ind_update_h(ic,jc,2) ) then
                              Atm(n)%ind_update_h(ic,jc,:) = (/i, j/)
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
!            call set_eta(npz, ks, ptop, Atm(n)%ak, Atm(n)%bk)
!         enddo
!         if(master) write(6,*) "Hybrid sigma-p coordinate has been reset"
!     endif

      if (ntilesMe > 1) call switch_current_Atm(Atm(1))

 end subroutine fv_init
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
         
 subroutine fv_end(Atm)

    type(fv_atmos_type), intent(inout) :: Atm(:)

    integer :: n

    call timing_off('TOTAL')
    call timing_prt( mpp_pe() )

    call fv_restart_end(Atm)
    call fv_io_exit()

  ! Free temporary memory from sw_core routines

  ! Deallocate
    call grid_utils_end

    do n = 1, ntilesMe
       call deallocate_fv_atmos_type(Atm(n))
    end do


 end subroutine fv_end
!-------------------------------------------------------------------------------

!!!! WARNING: New logic in this section may cause problems with nesting code

!-------------------------------------------------------------------------------
!
!     run_setup :: initialize run from namelist
!
      subroutine run_setup(Atm, dt_atmos)
      type(fv_atmos_type), intent(inout), target :: Atm(:)
      real, intent(in)                   :: dt_atmos

      character(len=80) :: filename, tracerName, errString, nested_nml_filename
      integer :: ios, ierr, f_unit, unit
      logical :: exists

      real :: dim0 = 180.           ! base dimension
      real :: dt0  = 1800.          ! base time step
      real :: ns0  = 5.             ! base nsplit for base dimension 
                                    ! For cubed sphere 5 is better
      real :: umax = 350.           ! max wave speed for grid_type>3
      real :: dimx, dl, dp, dxmin, dymin, d_fac

      integer :: n0split
      integer :: n, i, parent_grid_num = -1

      integer :: pe_counter = 0

!     local version of these variables to allow PGI compiler to compile
      character(len=128) :: res_latlon_dynamics = ''
      character(len=128) :: res_latlon_tracers  = ''
      character(len=80)  :: grid_name = ''
      character(len=120) :: grid_file = ''

      character(len=INPUT_STR_LENGTH), dimension(:), allocatable :: nested_nml_file

      namelist /mpi_nml/npes_x,npes_y  ! Use of this namelist is deprecated
      namelist /fv_grid_nml/ grid_name, grid_file
      namelist /fv_core_nml/npx, npy, ntiles, npz, npz_rst, layout, io_layout, ncnst, nwat,  &
                            p_fac, a_imp, k_split, n_split, m_split, q_split, print_freq, do_schmidt,      &
                            hord_mt, hord_vt, hord_tm, hord_dp, hord_tr, shift_fac, stretch_fac, target_lat, target_lon, &
                            kord_mt, kord_wz, kord_tm, kord_tr, fv_debug, fv_land, nudge,  &
                            external_ic, ncep_ic, fv_diag_ic, res_latlon_dynamics, res_latlon_tracers, &
                            scale_z, w_max, z_min, dddmp, d2_bg, d4_bg, vtdm4, d_ext, beta, non_ortho, n_sponge, &
                            warm_start, adjust_dry_mass, mountain, d_con, nord, convert_ke, use_old_omega, &
                            dry_mass, grid_type, do_Held_Suarez, consv_te, fill, tq_filter, filter_phys, fill_dp, fill_wz, &
                            range_warn, dwind_2d, inline_q, z_tracer, reproduce_sum, adiabatic, do_vort_damp, no_dycore,   &
                            replace_w, tau, tau_h2o, rf_center, nf_omega, hydrostatic, fv_sg_adj, breed_vortex_inline,  &
                            na_init, hybrid_z, Make_NH, n_zs_filter, nord_zs_filter, reset_eta,         &
                            a2b_ord, remap_t, p_ref, d2_bg_k1, d2_bg_k2,  &
#ifdef MARS_GCM
                            sponge_damp, reference_sfc_pres,                     &
#endif
                            c2l_ord, dx_const, dy_const, umax, deglat,      &
                            deglon_start, deglon_stop, deglat_start, deglat_stop, &
                            phys_hydrostatic, make_hybrid_z, old_divg_damp, &
                            nested, twowaynest, parent_grid_num, parent_tile, &
                            refinement, nestbctype, nestupdate, nsponge, s_weight, &
                            ioffset, joffset

      namelist /test_case_nml/test_case,alpha

! Make alpha = 0 the default:
      alpha = 0.
      test_case = 11   ! (USGS terrain)

      filename = "input.nml"

      inquire(file=filename,exist=exists)
      if (.not. exists) then  ! This will be replaced with fv_error wrapper
        if(master) write(6,*) "file ",trim(filename)," doesn't exist" 
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
   ! Look for deprecated mpi_nml
      read (input_nml_file,mpi_nml,iostat=ios)
      if (ios == 0) call mpp_error(FATAL,'mpi_nml is deprecated. Use layout in fv_core_nml')
!      rewind (f_unit)
#else
      f_unit=open_namelist_file()
      rewind (f_unit)
   ! Read Main namelist
      read (f_unit,fv_grid_nml,iostat=ios)
      ierr = check_nml_error(ios,'fv_grid_nml')
   ! Read Test_Case namelist
      rewind (f_unit)
      read (f_unit,test_case_nml,iostat=ios)
      ierr = check_nml_error(ios,'test_case_nml')
   ! Look for deprecated mpi_nml
      rewind (f_unit)
      read (f_unit,mpi_nml,iostat=ios)
      if (ios == 0) call mpp_error(FATAL,'mpi_nml is deprecated. Use layout in fv_core_nml')
      rewind (f_unit)
#endif

      unit = stdlog()
      write(unit, nml=fv_grid_nml)
      write(unit, nml=test_case_nml)

      do n=1,size(Atm)

         call switch_current_Atm(Atm(n), .false.)
         Atm(n)%grid_number = n

         if (n == 1) then
#ifdef INTERNAL_FILE_NML
   ! Read FVCORE namelist 
      read (input_nml_file,fv_core_nml,iostat=ios)
      ierr = check_nml_error(ios,'fv_core_nml')
#else
   ! Read FVCORE namelist 
      read (f_unit,fv_core_nml,iostat=ios)
      ierr = check_nml_error(ios,'fv_core_nml')
      call close_file(f_unit)
#endif         
          if (len_trim(grid_file) /= 0) Atm(n)%grid_file = grid_file
          if (len_trim(grid_name) /= 0) Atm(n)%grid_name = grid_name
          if (len_trim(res_latlon_dynamics) /= 0) Atm(n)%res_latlon_dynamics = res_latlon_dynamics
          if (len_trim(res_latlon_tracers)  /= 0) Atm(n)%res_latlon_tracers = res_latlon_tracers
         else
            write(nested_nml_filename, '(A, I1, A)') 'input', n, '.nml'
            call read_ascii_file(nested_nml_filename, INPUT_STR_LENGTH, nested_nml_file)
            read(nested_nml_file,fv_core_nml,iostat=ios)
            ierr = check_nml_error(ios,'fv_core_nml')
            deallocate (nested_nml_file)
            Atm(n)%grid_file = Atm(1)%grid_file
            Atm(n)%grid_name = Atm(1)%grid_name
            Atm(n)%res_latlon_dynamics = Atm(1)%res_latlon_dynamics
            Atm(n)%res_latlon_tracers = Atm(1)%res_latlon_tracers
         endif

         write(unit, nml=fv_core_nml)

         if (n > 1) then
            test_case = Atm(1)%test_case
            alpha     = Atm(1)%alpha

            !*** single tile for Cartesian grids
            if (grid_type>3) then
               ntiles=1
               non_ortho = .false.
               nf_omega = 0
            endif
         end if

         npes_x = layout(1)
         npes_y = layout(2)
         io_domain_layout = io_layout

         if (.not. nested) Atm(n)%npx_global = npx

         ! Define n_split if not in namelist
         if (ntiles==6) then
            dimx = 4.0*(npx-1)
#ifdef MARS_GCM
         ns0 = 8
#else
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
#endif
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
           n_split = nint( real(n0split)/real(k_split) * stretch_fac + 0.5 )
           if(master) write(6,*) 'For k_split (remapping)=', k_split
           if(master) write(6,198) 'n_split is set to ', n_split, ' for resolution-dt=',npx,npy,ntiles,dt_atmos
      else
          if(master) write(6,199) 'Using n_split from the namelist: ', n_split
      endif

      if (Atm(n)%nested) then
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
         if ( Atm(n)%parent_grid%grid_type < 3 .and. &
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

!!$         !A nested grid should have the same test case as its parent
!!$         Atm(n)%alpha = Atm(n)%parent_grid%alpha
!!$         Atm(n)%test_case = Atm(n)%parent_grid%test_case

         if (nestupdate == 1 .or. nestupdate == 2) then

            if (mod(npx-1,refinement) /= 0 .or. mod(npy-1,refinement) /= 0) then
               call mpp_error(WARNING, 'npx-1 or npy-1 is not evenly divisible by the refinement ratio; averaging update cannot be mass-conservative.')
            end if

         end if

         if ( consv_te > 0.) then
            call mpp_error(FATAL, 'The global energy fixer cannot be used on a nested grid. consv_te must be set to 0.')
         end if
         if ( d_con > 0. ) then
            call mpp_error(FATAL, 'The divergence damping heating cannot be used on a nested grid. d_con must be set to 0.')
         end if

         Atm(n)%refinement_of_global = Atm(n)%refinement * Atm(n)%parent_grid%refinement_of_global
         max_refinement_of_global = max(Atm(n)%refinement_of_global,max_refinement_of_global)
         Atm(n)%npx_global = Atm(n)%refinement * Atm(n)%parent_grid%npx_global

      else
         Atm(n)%ioffset                = -999
         Atm(n)%joffset                = -999   
         Atm(n)%parent_tile            = -1      
         Atm(n)%refinement             = -1
      end if

      if (Atm(n)%nested) then
         if (Atm(n)%grid_type >= 4 .and. Atm(n)%parent_grid%grid_type >= 4) then
            Atm(n)%dx_const = Atm(n)%parent_grid%dx_const / real(Atm(n)%refinement)
            Atm(n)%dy_const = Atm(n)%parent_grid%dy_const / real(Atm(n)%refinement)
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
        if (master) write(6,*) " fv_control: using original values for divergence damping "
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

      if ( .not.hydrostatic ) then
         if ( m_split==0 ) then
              m_split = 1. + abs(dt_atmos)/real(k_split*n_split)
         endif
         if(master) write(*,198) 'm_split is set to ', m_split
         if ( a_imp > 0.5 ) then
              if(master) write(*,*) 'Off center implicit scheme param=', a_imp, ' p_fac=', p_fac
         endif
      endif

      if(master) then
         write(6,199) 'Using n_sponge : ', n_sponge
         write(6,197) 'Using non_ortho : ', non_ortho
      endif

 197  format(A,l7)
 198  format(A,i2.2,A,i4.4,'x',i4.4,'x',i1.1,'-',f9.3)
 199  format(A,i2.2)

      if (.not. nested) alpha = alpha*pi


      allocate(Atm(n)%child_grids(size(Atm)))
      Atm(N)%child_grids = .false.

      !Determine concurrency; if first grid uses up all tiles then we are NOT concurrent.
      npes_this_grid =  npes_x*npes_y*ntiles
      if (n == 1) concurrent =  (npes_this_grid < mpp_npes())

      !Set up PElists now, before calling domain_decomp
      if (concurrent) then
         if (pe_counter+npes_this_grid-1 >= mpp_npes()) then
            if (gid == mpp_root_pe()) print*, 'pes for this grid     = ', npes_this_grid
            if (gid == mpp_root_pe()) print*, 'pes available         = ', mpp_npes()
            if (gid == mpp_root_pe()) print*, 'pes already assigned  = ', pe_counter
            call mpp_error(FATAL, 'DOMAIN_DECOMP: grid for concurrent nesting requires' // &
                 'more PEs than are available.')
         endif
      else
         npes_this_grid = mpp_npes()
      endif
      
      allocate(Atm(n)%pelist(npes_this_grid))
      Atm(n)%pelist=(/ (i,i=pe_counter,pe_counter+npes_this_grid-1)    /)
      call mpp_declare_pelist(Atm(n)%pelist)
      if (concurrent) then
         pe_counter = pe_counter + npes_this_grid
      endif
            
   enddo

   !Set pelists
   do n=1,size(Atm)
      if (ANY(Atm(n)%pelist == gid)) then
         if (concurrent) then
            call mpp_set_current_pelist(Atm(n)%pelist)
         endif
         masterproc = Atm(n)%pelist(1)
         master = (masterproc == gid)
         grids_on_this_pe(n) = .true.

         if (.not. allocated(pelist_all)) then
            !The next call is needed to get commID
            allocate( pelist_all(size(Atm(n)%pelist)) )
            call mpp_get_current_pelist(pelist_all, commID=commID)
            call mp_start(commID)
         endif
      else
         grids_on_this_pe(n) = .false.
      endif

      if (Atm(n)%nested) then
         Atm(n)%parent_proc = ANY(Atm(n)%parent_grid%pelist == gid)
         Atm(n)%child_proc = ANY(Atm(n)%pelist == gid)
      endif
   enddo

   do n=1,size(Atm)
      
      call switch_current_Atm(Atm(n),.false.)
      call domain_decomp(npx,npy,ntiles,ng,grid_type,Atm(n)%nested)
   enddo

   call broadcast_domains(Atm)

   do n=1,size(Atm)
      call switch_current_Atm(Atm(n))

      if (nested) then
         if (mod(npx-1 , refinement) /= 0 .or. mod(npy-1, refinement) /= 0) &
              call mpp_error(FATAL, 'npx or npy not an even refinement of its coarse grid.')
         
         !Pelist needs to be set to ALL (which should have been done
         !in broadcast_domains) to get this to work
         call mpp_define_nest_domains(Atm(n)%nest_domain, Atm(n)%domain, Atm(parent_grid_num)%domain, &
              7, parent_tile, &
              1, npx-1, 1, npy-1,                  & !Grid cells, not points
              ioffset, ioffset + (npx-1)/refinement - 1, &
              joffset, joffset + (npy-1)/refinement - 1,         &
              (/ (i,i=0,npes-1)  /), 0, name="nest_domain") !What pelist to use?

         Atm(parent_grid_num)%child_grids(n) = .true.

         if (Atm(n)%nestbctype > 1) then

            !This check is due to a bug which has not yet been identified. Beware.
!!$            if (Atm(n)%parent_grid%hord_tr == 7) &
!!$                 call mpp_error(FATAL, "Flux-form nested  BCs (nestbctype > 1) should not use hord_tr == 7 (on parent grid), since there is no guarantee of tracer mass conservation with this option.")

            if (.not. concurrent) call mpp_error(FATAL,"Flux-form nested BCs (nestbctype > 1) only compatible with concurrent nesting.")
            if (Atm(n)%q_split > 0 .and. Atm(n)%parent_grid%q_split > 0) then
               if (mod(Atm(n)%q_split,Atm(n)%parent_grid%q_split) /= 0) call mpp_error(FATAL, &
                    "Flux-form nested BCs (nestbctype > 1) require q_split on the nested grid to be evenly divisible by that on the coarse grid.")
            endif
            if (mod((Atm(n)%npx-1),Atm(n)%refinement) /= 0 .or. mod((Atm(n)%npy-1),Atm(n)%refinement) /= 0) call mpp_error(FATAL, &
                 "Flux-form nested BCs (nestbctype > 1) requires npx and npy to be one more than a multiple of the refinement ratio.")
            Atm(n)%parent_grid%do_flux_BCs = .true.
            if (Atm(n)%nestbctype == 3 .or. Atm(n)%nestbctype == 4) Atm(n)%parent_grid%do_2way_flux_BCs = .true.
         endif
      end if

   end do

      if (concurrent) then
         do n=1,size(Atm)
            if (ANY(Atm(n)%pelist == gid)) then
               call mpp_set_current_pelist(Atm(n)%pelist)
            endif
         enddo
      endif

  end subroutine run_setup


       
end module fv_control_mod
