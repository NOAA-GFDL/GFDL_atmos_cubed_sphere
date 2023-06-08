!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be
!* useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

module fv_restart_mod

  !<OVERVIEW>
  ! Restart facilities for FV core
  !</OVERVIEW>
  !<DESCRIPTION>
  ! This module writes and reads restart files for the FV core. Additionally
  ! it provides setup and calls routines necessary to provide a complete restart
  ! for the model.
  !</DESCRIPTION>

  use constants_mod,       only: kappa, pi=>pi_8, rdgas, grav, rvgas, cp_air
  use fv_arrays_mod,       only: radius, omega ! scaled for small earth
  use fv_arrays_mod,       only: fv_atmos_type, fv_nest_type, fv_grid_bounds_type, R_GRID
  use fv_io_mod,           only: fv_io_init, fv_io_read_restart, fv_io_write_restart, &
                                 remap_restart, fv_io_write_BCs, fv_io_read_BCs
  use fv_grid_utils_mod,   only: ptop_min, fill_ghost, g_sum, &
                                 make_eta_level, cubed_to_latlon, great_circle_dist
  use fv_diagnostics_mod,  only: prt_maxmin, gn
  use init_hydro_mod,      only: p_var
  use mpp_domains_mod,     only: mpp_update_domains, domain2d, DGRID_NE
  use mpp_domains_mod,     only: mpp_get_compute_domain, mpp_get_data_domain, mpp_get_global_domain
  use mpp_domains_mod,     only: CENTER, CORNER, NORTH, EAST,  mpp_get_C2F_index, WEST, SOUTH
  use mpp_domains_mod,     only: mpp_global_field
  use mpp_mod,             only: mpp_chksum, stdout, mpp_error, FATAL, NOTE
  use mpp_mod,             only: get_unit, mpp_sum, mpp_broadcast, mpp_max
  use mpp_mod,             only: mpp_get_current_pelist, mpp_npes, mpp_set_current_pelist
  use mpp_mod,             only: mpp_send, mpp_recv, mpp_sync_self, mpp_pe, mpp_sync
  use fms2_io_mod,         only: file_exists, set_filename_appendix, FmsNetcdfFile_t, open_file, close_file
  use test_cases_mod,      only: alpha, init_case, init_double_periodic!, init_latlon
  use fv_mp_mod,           only: is_master, mp_reduce_min, mp_reduce_max, corners_YDir => YDir, fill_corners, tile_fine, global_nest_domain
  use fv_surf_map_mod,     only: sgh_g, oro_g
  use tracer_manager_mod,  only: get_tracer_index, get_tracer_names, set_tracer_profile
  use field_manager_mod,   only: MODEL_ATMOS
  use external_ic_mod,     only: get_external_ic
  use fv_eta_mod,          only: compute_dz_var, compute_dz_L32, set_hybrid_z
  use fv_surf_map_mod,     only: del2_cubed_sphere, del4_cubed_sphere
  use boundary_mod,        only: fill_nested_grid, nested_grid_BC, update_coarse_grid
  use fv_timing_mod,       only: timing_on, timing_off
  use fv_treat_da_inc_mod, only: read_da_inc
  use fv_regional_mod,     only: write_full_fields
  use coarse_grained_restart_files_mod, only: fv_io_write_restart_coarse

  implicit none
  private

  public :: fv_restart_init, fv_restart_end, fv_restart, fv_write_restart

  real(kind=R_GRID), parameter :: cnst_0p20=0.20d0
  !--- private data type
  logical                       :: module_is_initialized = .FALSE.

contains

  !#####################################################################
  ! <SUBROUTINE NAME="fv_restart_init">
  !
  ! <DESCRIPTION>
  ! Initialize the fv core restart facilities
  ! </DESCRIPTION>
  !
  subroutine fv_restart_init()

    call fv_io_init()
    module_is_initialized = .TRUE.

  end subroutine fv_restart_init
  ! </SUBROUTINE> NAME="fv_restart_init"


    !#####################################################################
  ! <SUBROUTINE NAME="fv_restart">
  !
  ! <DESCRIPTION>
  ! The fv core restart facility
  ! </DESCRIPTION>
  !
  subroutine fv_restart(fv_domain, Atm, dt_atmos, seconds, days, cold_start, grid_type, &
                        this_grid)
    type(domain2d),      intent(inout) :: fv_domain
    type(fv_atmos_type), intent(inout) :: Atm(:)
    real,                intent(in)    :: dt_atmos
    integer,             intent(out)   :: seconds
    integer,             intent(out)   :: days
    logical,             intent(inout)    :: cold_start
    integer,             intent(in)    :: grid_type, this_grid

    integer :: i, j, k, l, m, n, ntileMe, nt, iq
    integer :: isc, iec, jsc, jec, ncnst, ntprog, ntdiag
    integer :: isd, ied, jsd, jed, npz
    integer isd_p, ied_p, jsd_p, jed_p, isc_p, iec_p, jsc_p, jec_p, isg, ieg, jsg,jeg, npx_p, npy_p
    real, allocatable :: g_dat(:,:,:)

    integer :: unit
    real, allocatable :: dz1(:)
    real rgrav, f00, ztop, pertn, ph
    logical :: hybrid
    character(len=128):: tname, errstring, fname, tracer_name
    character(len=120):: fname_ne, fname_sw
    character(len=3) :: gn
    character(len=10) :: inputdir
    character(len=6) :: gnn

    integer :: npts, sphum, aero_id
    integer, allocatable :: pelist(:), global_pelist(:), smoothed_topo(:)
    real    :: sumpertn
    real    :: zvir, nbg_inv

    type(FmsNetcdfFile_t) :: fileobj
    logical :: do_read_restart = .false.
    logical :: do_read_restart_bc = .false.
    integer, allocatable :: ideal_test_case(:), new_nest_topo(:)
    integer :: nest_level
    integer, allocatable :: BC_remap_level(:)

    rgrav = 1. / grav

    if(.not.module_is_initialized) call mpp_error(FATAL, 'You must call fv_restart_init.')

    ntileMe = size(Atm(:))
    allocate(smoothed_topo(ntileme))
    smoothed_topo(:) = 0
    allocate(ideal_test_case(ntileme))
    ideal_test_case(:) = 0
    allocate(new_nest_topo(ntileme))
    new_nest_topo(:) = 0

    do n = 1, ntileMe

       isd = Atm(n)%bd%isd
       ied = Atm(n)%bd%ied
       jsd = Atm(n)%bd%jsd
       jed = Atm(n)%bd%jed
       isc = Atm(n)%bd%isc
       iec = Atm(n)%bd%iec
       jsc = Atm(n)%bd%jsc
       jec = Atm(n)%bd%jec
       ncnst = Atm(n)%ncnst
       if( is_master() ) write(*,*) 'in fv_restart ncnst=', ncnst
       npz = Atm(n)%npz
       ntprog = size(Atm(n)%q,4)
       ntdiag = size(Atm(n)%qdiag,4)

       !1. sort out restart, external_ic, and cold-start (idealized) plus initialize tracers
       if (Atm(n)%neststruct%nested) then
          write(fname,   '(A, I2.2, A)') 'INPUT/fv_core.res.nest', Atm(n)%grid_number, '.nc'
          write(fname_ne,'(A, I2.2, A)') 'INPUT/fv_BC_ne.res.nest', Atm(n)%grid_number, '.nc'
          write(fname_sw,'(A, I2.2, A)') 'INPUT/fv_BC_sw.res.nest', Atm(n)%grid_number, '.nc'
          if (is_master()) print*, 'Searching for nested grid BC files ', trim(fname_ne), ' ', trim (fname_sw)
          do_read_restart = file_exists(fname)
          do_read_restart_bc = file_exists(fname_ne) .and. file_exists(fname_sw)
          if (is_master()) then
             print*, 'FV_RESTART: ', n, do_read_restart, do_read_restart_bc
             if (.not. do_read_restart_bc) write(*,*) 'BC files not found, re-generating nested grid boundary conditions'
          endif
          Atm(N)%neststruct%first_step = .not. do_read_restart_bc
       else
          fname='INPUT/fv_core.res.nc'
          do_read_restart = open_file(fileobj, fname, "read", is_restart=.true.)
          if (do_read_restart) call close_file(fileobj)
          if (is_master()) print*, 'FV_RESTART: ', n, do_read_restart, do_read_restart_bc
       endif

       !initialize tracers
       do nt = 1, ntprog
          call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
          ! set all tracers to an initial profile value
          call set_tracer_profile (MODEL_ATMOS, nt, Atm(n)%q(:,:,:,nt))
       enddo
       do nt = ntprog+1, ntprog+ntdiag
          call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
          ! set all tracers to an initial profile value
          call set_tracer_profile (MODEL_ATMOS, nt, Atm(n)%qdiag(:,:,:,nt))
       enddo

       !2. Register restarts
       !No longer need to register restarts in fv_restart_mod with fms2_io implementation

       ! The two calls are needed until everything uses fms2io
       if (Atm(n)%neststruct%nested .and. n==this_grid) then
          write(gnn,'(A4, I2.2)') "nest", Atm(n)%grid_number
          call set_filename_appendix(gnn)
       endif

       !3preN. Topography BCs for nest, including setup for blending

       if (Atm(n)%neststruct%nested) then
          if (.not. allocated(pelist)) then
             allocate(pelist(0:mpp_npes()-1))
             call mpp_get_current_pelist(pelist)
          endif
          call mpp_set_current_pelist() !global
          call mpp_broadcast(Atm(n)%flagstruct%external_ic,Atm(n)%pelist(1))
          call mpp_sync()
          call mpp_set_current_pelist(pelist)
          if ( ( smoothed_topo(Atm(n)%parent_grid%grid_number) > 0 .or. &
                 .not. do_read_restart_bc .or. &
                 Atm(n)%flagstruct%external_ic  ) ) then
             new_nest_topo(n) = 1
             if (n==this_grid .or. this_grid==Atm(n)%parent_grid%grid_number) then
                call fill_nested_grid_topo(Atm(n), n==this_grid)
             endif

! This sets the nest BCs. Ideally, it should be done outside of this loop but
! it is moved here -before ex_ic- to avoid repro issues on the nest
             do nest_level=1,Atm(this_grid)%neststruct%num_nest_level

                if (Atm(this_grid)%neststruct%nested .AND. Atm(this_grid)%neststruct%nlevel==nest_level)then
                   call nested_grid_BC(Atm(this_grid)%ps, Atm(this_grid)%parent_grid%ps, global_nest_domain, &
                        Atm(this_grid)%neststruct%ind_h, Atm(this_grid)%neststruct%wt_h, 0, 0, &
                        Atm(this_grid)%npx, Atm(this_grid)%npy, Atm(this_grid)%bd, 1, Atm(this_grid)%npx-1, 1,&
                        Atm(this_grid)%npy-1,nest_level=Atm(this_grid)%neststruct%nlevel)
                   call nested_grid_BC(Atm(this_grid)%phis, Atm(this_grid)%parent_grid%phis, global_nest_domain, &
                        Atm(this_grid)%neststruct%ind_h, Atm(this_grid)%neststruct%wt_h, 0, 0, &
                        Atm(this_grid)%npx, Atm(this_grid)%npy, Atm(this_grid)%bd, 1, Atm(this_grid)%npx-1, 1, &
                        Atm(this_grid)%npy-1,nest_level=Atm(this_grid)%neststruct%nlevel)
               endif

                if (ANY (Atm(this_grid)%neststruct%child_grids) .AND. Atm(this_grid)%neststruct%nlevel==nest_level-1) then
                   call nested_grid_BC(Atm(this_grid)%ps, global_nest_domain, 0, 0, nest_level=Atm(this_grid)%neststruct%nlevel+1)
                   call nested_grid_BC(Atm(this_grid)%phis, global_nest_domain, 0, 0, nest_level=Atm(this_grid)%neststruct%nlevel+1)
                endif

            enddo
          endif
       endif

       !This call still appears to be necessary to get isd, etc. correct
       !call switch_current_Atm(Atm(n)) !TODO should NOT be necessary now that we manually set isd, etc.

       !--- call fv_io_register_restart to register restart field to be written out in fv_io_write_restart
       !if (n==this_grid) call fv_io_register_restart(Atm(n)%domain,Atm(n:n))
       !if (Atm(n)%neststruct%nested) call fv_io_register_restart_BCs(Atm(n)) !TODO put into fv_io_register_restart

       if (n==this_grid) then

          !3. External_ic
          if (Atm(n)%flagstruct%external_ic) then

             if( is_master() ) write(*,*) 'Calling get_external_ic'
             call get_external_ic(Atm(n), .not. do_read_restart)
             if( is_master() ) write(*,*) 'IC generated from the specified external source'

          !4. Restart
          elseif (do_read_restart) then

             if ( Atm(n)%flagstruct%npz_rst /= 0 .and. Atm(n)%flagstruct%npz_rst /= Atm(n)%npz ) then
                !Remap vertically the prognostic variables for the chosen vertical resolution
                if( is_master() ) then
                   write(*,*) ' '
                   write(*,*) '***** Important Note from FV core ********************'
                   write(*,*) 'Remapping dynamic IC from', Atm(n)%flagstruct%npz_rst, 'levels to ', Atm(n)%npz,'levels'
                   write(*,*) '***** End Note from FV core **************************'
                   write(*,*) ' '
                endif
                call remap_restart( Atm(n:n) )
                if( is_master() ) write(*,*) 'Done remapping dynamical IC'
             else
                if( is_master() ) write(*,*) 'Warm starting, calling fv_io_restart'
                call fv_io_read_restart(Atm(n)%domain_for_read,Atm(n:n))
                !====== PJP added DA functionality ======
                if (Atm(n)%flagstruct%read_increment) then
                   ! print point in middle of domain for a sanity check
                   i = (Atm(n)%bd%isc + Atm(n)%bd%iec)/2
                   j = (Atm(n)%bd%jsc + Atm(n)%bd%jec)/2
                   k = Atm(n)%npz/2
                   if( is_master() ) write(*,*) 'Calling read_da_inc',Atm(n)%pt(i,j,k)
                   call read_da_inc(Atm(n), Atm(n)%domain)
                   if( is_master() ) write(*,*) 'Back from read_da_inc',Atm(n)%pt(i,j,k)
                endif
                !====== end PJP added DA functionailty======
             endif

             seconds = 0; days = 0   ! Restart needs to be modified to record seconds and days.

             if (Atm(n)%neststruct%nested) then
                if ( Atm(n)%flagstruct%npz_rst /= 0 .and. Atm(n)%flagstruct%npz_rst /= npz ) then
                   call mpp_error(FATAL, "Remap-restart not implemented for nests.")
                endif
                if (do_read_restart_BC) call fv_io_read_BCs(Atm(n))
                call mpp_update_domains(Atm(n)%u, Atm(n)%v, Atm(n)%domain, gridtype=DGRID_NE, complete=.true.)
             endif

             if ( Atm(n)%flagstruct%mountain ) then
                ! !!! Additional terrain filter -- should not be called repeatedly !!!
                if ( Atm(n)%flagstruct%n_zs_filter > 0 ) then
                   if ( Atm(n)%flagstruct%nord_zs_filter == 2 ) then
                      !!! TODO: move this block into its own routine or CLEAN UP these subroutine calls
                      call del2_cubed_sphere(Atm(n)%npx, Atm(n)%npy, Atm(n)%phis, &
                           Atm(n)%gridstruct%area_64, Atm(n)%gridstruct%dx, Atm(n)%gridstruct%dy,   &
                           Atm(n)%gridstruct%dxc, Atm(n)%gridstruct%dyc, Atm(n)%gridstruct%sin_sg, &
                           Atm(n)%flagstruct%n_zs_filter, cnst_0p20*Atm(n)%gridstruct%da_min, &
                           .false., oro_g, Atm(n)%gridstruct%bounded_domain, Atm(n)%domain, Atm(n)%bd)
                      if ( is_master() ) write(*,*) 'Warning !!! del-2 terrain filter has been applied ', &
                           Atm(n)%flagstruct%n_zs_filter, ' times'
                   else if( Atm(n)%flagstruct%nord_zs_filter == 4 ) then
                      call del4_cubed_sphere(Atm(n)%npx, Atm(n)%npy, Atm(n)%phis, Atm(n)%gridstruct%area_64, &
                           Atm(n)%gridstruct%dx, Atm(n)%gridstruct%dy,   &
                           Atm(n)%gridstruct%dxc, Atm(n)%gridstruct%dyc, Atm(n)%gridstruct%sin_sg, &
                           Atm(n)%flagstruct%n_zs_filter, .false., oro_g, Atm(n)%gridstruct%bounded_domain, &
                           Atm(n)%domain, Atm(n)%bd)
                      if ( is_master() ) write(*,*) 'Warning !!! del-4 terrain filter has been applied ', &
                           Atm(n)%flagstruct%n_zs_filter, ' times'
                   endif
                   if ( Atm(n)%flagstruct%fv_land .and. allocated(sgh_g) .and. allocated(oro_g) ) then
                      do j=jsc,jec
                         do i=isc,iec
                            Atm(n)%sgh(i,j) = sgh_g(i,j)
                            Atm(n)%oro(i,j) = oro_g(i,j)
                         enddo
                      enddo
                   endif
                endif
                call mpp_update_domains( Atm(n)%phis, Atm(n)%domain, complete=.true. )
             else
                Atm(n)%phis = 0.
                if( is_master() ) write(*,*) 'phis set to zero'
             endif !mountain

             !5. Idealized test case
          elseif (Atm(n)%flagstruct%is_ideal_case) then

             ideal_test_case(n) = 1

             if ( Atm(n)%flagstruct%make_hybrid_z ) then
                hybrid = .false.
             else
                hybrid = Atm(n)%flagstruct%hybrid_z
             endif
             if (grid_type < 4) then
                if ( .not. Atm(n)%flagstruct%external_ic ) then
                   call init_case(Atm(n)%u,Atm(n)%v,Atm(n)%w,Atm(n)%pt,Atm(n)%delp,Atm(n)%q, &
                        Atm(n)%phis, Atm(n)%ps,Atm(n)%pe, Atm(n)%peln,Atm(n)%pk,Atm(n)%pkz, &
                        Atm(n)%uc,Atm(n)%vc, Atm(n)%ua,Atm(n)%va,        &
                        Atm(n)%ak, Atm(n)%bk, Atm(n)%gridstruct, Atm(n)%flagstruct,&
                        Atm(n)%npx, Atm(n)%npy, npz, Atm(n)%ng, &
                        ncnst, Atm(n)%flagstruct%nwat,  &
                        Atm(n)%flagstruct%ndims, Atm(n)%flagstruct%ntiles, &
                        Atm(n)%flagstruct%dry_mass, &
                        Atm(n)%flagstruct%mountain,       &
                        Atm(n)%flagstruct%moist_phys, Atm(n)%flagstruct%hydrostatic, &
                        hybrid, Atm(n)%delz, Atm(n)%ze0, &
                        Atm(n)%flagstruct%adiabatic, Atm(n)%ks, Atm(n)%neststruct%npx_global, &
                        Atm(n)%ptop, Atm(n)%domain, Atm(n)%tile_of_mosaic, Atm(n)%bd)
                endif
             elseif (grid_type == 4) then
                call init_double_periodic(Atm(n)%u,Atm(n)%v,Atm(n)%w,Atm(n)%pt, &
                     Atm(n)%delp,Atm(n)%q,Atm(n)%phis, Atm(n)%ps,Atm(n)%pe, &
                     Atm(n)%peln,Atm(n)%pk,Atm(n)%pkz, &
                     Atm(n)%uc,Atm(n)%vc, Atm(n)%ua,Atm(n)%va,        &
                     Atm(n)%ak, Atm(n)%bk, &
                     Atm(n)%gridstruct, Atm(n)%flagstruct, &
                     Atm(n)%npx, Atm(n)%npy, npz, Atm(n)%ng, &
                     ncnst, Atm(n)%flagstruct%nwat,  &
                     Atm(n)%flagstruct%ndims, Atm(n)%flagstruct%ntiles, &
                     Atm(n)%flagstruct%dry_mass, Atm(n)%flagstruct%mountain, &
                     Atm(n)%flagstruct%moist_phys, Atm(n)%flagstruct%hydrostatic, &
                     hybrid, Atm(n)%delz, Atm(n)%ze0, Atm(n)%ks, Atm(n)%ptop, &
                     Atm(n)%domain, Atm(n)%tile_of_mosaic, Atm(n)%bd)
                if( is_master() ) write(*,*) 'Doubly Periodic IC generated'
             elseif (grid_type == 5 .or. grid_type == 6) then
                call mpp_error(FATAL, "Idealized test cases for grid_type == 5,6 (global lat-lon) grid not supported")
             endif

             !Turn this off on the nested grid if you are just interpolating topography from the coarse grid!
             !These parameters are needed in LM3/LM4, and are communicated through restart files
             if ( Atm(n)%flagstruct%fv_land  .and. allocated(sgh_g) .and. allocated(oro_g)) then
                do j=jsc,jec
                   do i=isc,iec
                      Atm(n)%sgh(i,j) = sgh_g(i,j)
                      Atm(n)%oro(i,j) = oro_g(i,j)
                   enddo
                enddo
             endif

             Atm(n)%u0 = Atm(n)%u
             Atm(n)%v0 = Atm(n)%v

          else

                call mpp_error(FATAL, "If there is no restart file, either external_ic or is_ideal_case must be set true.")

          endif !external_ic vs. restart vs. idealized


       endif !n==this_grid


          !!!! NOT NEEDED??
          !Currently even though we do fill in the nested-grid IC from
          ! init_case or external_ic we appear to overwrite it using
          !  coarse-grid data

!       endif  !end cold_start check

       !5n. Nesting setup (part I)

       !Broadcast data for nesting
       if (ntileMe > 1) then
          if (.not. allocated(pelist)) then
             allocate(pelist(0:mpp_npes()-1))
             call mpp_get_current_pelist(pelist)
          endif

          call mpp_set_current_pelist()!global
          !for remap BCs
          call mpp_broadcast(Atm(n)%ptop,Atm(n)%pelist(1))
          call mpp_broadcast(Atm(n)%ak,Atm(n)%npz+1,Atm(n)%pelist(1))
          call mpp_broadcast(Atm(n)%bk,Atm(n)%npz+1,Atm(n)%pelist(1))
          !smoothed_topo
          call mpp_broadcast(smoothed_topo(n),Atm(n)%pelist(1))

          call mpp_sync()
          call mpp_set_current_pelist(pelist)


          if (Atm(n)%neststruct%nested) then
             Atm(n)%neststruct%do_remap_BC(ntileMe) = .false.

             if (Atm(n)%npz /= Atm(n)%parent_grid%npz) then
                Atm(n)%neststruct%do_remap_BC(n) = .true.
             else
                do k=1,Atm(n)%npz+1
                   if (Atm(n)%ak(k) /= Atm(n)%parent_grid%ak(k)) then
                      Atm(n)%neststruct%do_remap_BC(n) = .true.
                      exit
                   endif
                   if (Atm(n)%bk(k) /= Atm(n)%parent_grid%bk(k)) then
                      Atm(n)%neststruct%do_remap_BC(n) = .true.
                      exit
                   endif
                enddo
             endif

             Atm(n)%parent_grid%neststruct%do_remap_BC(n) = Atm(n)%neststruct%do_remap_BC(n)
             if (is_master() .and. n==this_grid) then
                if (Atm(n)%neststruct%do_remap_BC(n)) then
                   print*, ' Remapping BCs ENABLED on grid', n
                else
                   print*, ' Remapping BCs DISABLED (not necessary) on grid', n
                endif
                write(*,'(A, I3, A, F8.2, A)') ' Nested grid ', n, ',  ptop = ', Atm(n)%ak(1), ' Pa'
                write(*,'(A, I3, A, F8.2, A)') ' Parent grid ', n, ',  ptop = ', Atm(n)%parent_grid%ak(1), ' Pa'
                if (Atm(n)%ak(1) < Atm(n)%parent_Grid%ak(1)) then
                   print*, ' WARNING nested grid top above parent grid top. May have problems with remapping BCs.'
                endif
             endif
          endif

       endif

    end do !break cycling loop to finish nesting setup

    ! The following section is simply to set up the logical do_remap_BC_level
    ! do_remap_BC_level is true if the BCs of any grid need remapping at a certain nest level
    ! This is to accomodate the BC communications in fv_nesting which is done by level
    ! and does not mean that all nests at a certain level are undergoing a BC remapping
    ! remapping is actually happening when do_remap_BC=.true.
    if (ntileMe>1) then
       if (.not. allocated (BC_remap_level))then
          allocate (BC_remap_level(Atm(this_grid)%neststruct%num_nest_level))
          BC_remap_level(:)=0
       endif

       do nest_level=1,Atm(this_grid)%neststruct%num_nest_level
          if (Atm(this_grid)%neststruct%nlevel==nest_level .AND. (Atm(this_grid)%neststruct%do_remap_BC(this_grid))) then
             BC_remap_level(nest_level) = 1
          endif
       enddo

       call mpp_set_current_pelist() !global

       if (.not. allocated (global_pelist))     allocate(global_pelist(mpp_npes()))

       call mpp_get_current_pelist(global_pelist)
       call mpp_max(BC_remap_level,Atm(this_grid)%neststruct%num_nest_level,global_pelist)
       call mpp_set_current_pelist(pelist)

       do nest_level=1,Atm(this_grid)%neststruct%num_nest_level
          Atm(this_grid)%neststruct%do_remap_BC_level(nest_level) = (BC_remap_level(nest_level) == 1 )
       enddo
    endif

    ! Topo twoway update
    do n = ntileMe,1,-1
       if (atm(n)%neststruct%twowaynest) then
          if (new_nest_topo(n) > 0) then
             if (Atm(n)%parent_grid%grid_number==this_grid) then    !only parent?!
                call twoway_topo_update(Atm(n), n==this_grid)
             elseif (n==this_grid .or. Atm(this_grid)%neststruct%nlevel==Atm(n)%neststruct%nlevel) then
                call twoway_topo_update(Atm(this_grid), n==this_grid)
             endif
          endif
       endif
    end do

    !6. Data Setup
    do n = 1, ntileMe

       if (n/=this_grid) cycle

       isd = Atm(n)%bd%isd
       ied = Atm(n)%bd%ied
       jsd = Atm(n)%bd%jsd
       jed = Atm(n)%bd%jed
       isc = Atm(n)%bd%isc
       iec = Atm(n)%bd%iec
       jsc = Atm(n)%bd%jsc
       jec = Atm(n)%bd%jec
       ncnst = Atm(n)%ncnst
       if( is_master() ) write(*,*) 'in fv_restart ncnst=', ncnst
       npz = Atm(n)%npz
       ntprog = size(Atm(n)%q,4)
       ntdiag = size(Atm(n)%qdiag,4)

       if (ideal_test_case(n) == 0) then
#ifdef SW_DYNAMICS
          Atm(n)%pt(:,:,:)=1.
#else
          if ( .not.Atm(n)%flagstruct%hybrid_z ) then
             if(Atm(n)%ptop/=Atm(n)%ak(1)) call mpp_error(FATAL,'FV restart: ptop not equal Atm(n)%ak(1)')
          else
             Atm(n)%ptop = Atm(n)%ak(1);  Atm(n)%ks = 0
          endif
          call p_var(npz,         isc,         iec,       jsc,     jec,   Atm(n)%ptop,     ptop_min,  &
               Atm(n)%delp, Atm(n)%delz, Atm(n)%pt, Atm(n)%ps, Atm(n)%pe, Atm(n)%peln,   &
               Atm(n)%pk,   Atm(n)%pkz, kappa, Atm(n)%q, Atm(n)%ng, &
               ncnst,  Atm(n)%gridstruct%area_64, Atm(n)%flagstruct%dry_mass,  &
               Atm(n)%flagstruct%adjust_dry_mass,  Atm(n)%flagstruct%mountain, &
               Atm(n)%flagstruct%moist_phys,  Atm(n)%flagstruct%hydrostatic, &
               Atm(n)%flagstruct%nwat, Atm(n)%domain, Atm(1)%flagstruct%adiabatic, Atm(n)%flagstruct%make_nh)
#endif
          if ( grid_type < 7 .and. grid_type /= 4 ) then
             ! Fill big values in the non-existing corner regions:
             !          call fill_ghost(Atm(n)%phis, Atm(n)%npx, Atm(n)%npy, big_number)
             do j=jsd,jed+1
                do i=isd,ied+1
                   Atm(n)%gridstruct%fc(i,j) = 2.*omega*( -cos(Atm(n)%gridstruct%grid(i,j,1))*cos(Atm(n)%gridstruct%grid(i,j,2))*sin(alpha) + &
                        sin(Atm(n)%gridstruct%grid(i,j,2))*cos(alpha) )
                enddo
             enddo
             do j=jsd,jed
                do i=isd,ied
                   Atm(n)%gridstruct%f0(i,j) = 2.*omega*( -cos(Atm(n)%gridstruct%agrid(i,j,1))*cos(Atm(n)%gridstruct%agrid(i,j,2))*sin(alpha) + &
                        sin(Atm(n)%gridstruct%agrid(i,j,2))*cos(alpha) )
                enddo
             enddo
          else
             f00 = 2.*omega*sin(Atm(n)%flagstruct%deglat/180.*pi)
             do j=jsd,jed+1
                do i=isd,ied+1
                   Atm(n)%gridstruct%fc(i,j) = f00
                enddo
             enddo
             do j=jsd,jed
                do i=isd,ied
                   Atm(n)%gridstruct%f0(i,j) = f00
                enddo
             enddo
          endif
          call mpp_update_domains( Atm(n)%gridstruct%f0, Atm(n)%domain )
          if ( Atm(n)%gridstruct%cubed_sphere .and. (.not. Atm(n)%gridstruct%bounded_domain))then
             call fill_corners(Atm(n)%gridstruct%f0, Atm(n)%npx, Atm(n)%npy, Corners_YDir)
          endif
       endif


!---------------------------------------------------------------------------------------------
! Transform the (starting) Eulerian vertical coordinate from sigma-p to hybrid_z
     if ( Atm(n)%flagstruct%hybrid_z ) then
       if ( Atm(n)%flagstruct%make_hybrid_z ) then
          allocate ( dz1(npz) )
          if( npz==32 ) then
              call compute_dz_L32(npz, ztop, dz1)
          else
              ztop = 45.E3
              call compute_dz_var(npz, ztop, dz1)
          endif
          call set_hybrid_z(isc, iec, jsc, jec, Atm(n)%ng, npz, ztop, dz1, rgrav,  &
                            Atm(n)%phis, Atm(n)%ze0)
          deallocate ( dz1 )
!         call prt_maxmin('ZE0', Atm(n)%ze0,  isc, iec, jsc, jec, 0, npz, 1.E-3)
!         call prt_maxmin('DZ0', Atm(n)%delz, isc, iec, jsc, jec, 0, npz, 1.   )
       endif
!      call make_eta_level(npz, Atm(n)%pe, area, Atm(n)%ks, Atm(n)%ak, Atm(n)%bk, Atm(n)%ptop)
     endif
!---------------------------------------------------------------------------------------------

     if (Atm(n)%flagstruct%do_aerosol) then
       aero_id = get_tracer_index(MODEL_ATMOS, 'aerosol')
       if (aero_id .gt. 0) then
         Atm(n)%q(isc:iec,jsc:jec,:,aero_id) = 0.0
       endif
     endif

     if (Atm(n)%flagstruct%add_noise > 0.) then
        write(errstring,'(A, E16.9)') "Adding thermal noise of amplitude ", Atm(n)%flagstruct%add_noise
        call mpp_error(NOTE, errstring)
        call random_seed
        npts = 0
        sumpertn = 0.
        do k=1,npz
        do j=jsc,jec
        do i=isc,iec
           call random_number(pertn)
           Atm(n)%pt(i,j,k) = Atm(n)%pt(i,j,k) + pertn*Atm(n)%flagstruct%add_noise
           npts = npts + 1
           sumpertn = sumpertn + pertn*Atm(n)%flagstruct%add_noise ** 2
        enddo
        enddo
        enddo
        call mpp_update_domains(Atm(n)%pt, Atm(n)%domain)
        call mpp_sum(sumpertn)
        call mpp_sum(npts)
        write(errstring,'(A, E16.9)') "RMS added noise: ", sqrt(sumpertn/npts)
        call mpp_error(NOTE, errstring)
     endif

     if (Atm(n)%flagstruct%fv_sg_adj > 0 .and. Atm(n)%flagstruct%sg_cutoff > 0) then
        !Choose n_sponge from first reference level above sg_cutoff
        do k=1,npz
           ph = Atm(n)%ak(k+1) +  Atm(n)%bk(k+1)*Atm(n)%flagstruct%p_ref
           if (ph > Atm(n)%flagstruct%sg_cutoff) exit
        enddo
        Atm(n)%flagstruct%n_sponge = min(k,npz)
        write(errstring,'(A, I3, A)') ' Override n_sponge: applying 2dz filter to ', k , ' levels'
        call mpp_error(NOTE, errstring)
     endif

      if (Atm(n)%grid_number > 1) then
         write(gn,'(A2, I1)') " g", Atm(n)%grid_number
      else
         gn = ''
      end if

      unit = stdout()
      !!!NOTE: Checksums not yet working in stand-alone regional model!!
      write(unit,*)
      write(unit,*) 'fv_restart u   ', trim(gn),' = ', mpp_chksum(Atm(n)%u(isc:iec,jsc:jec,:))
      write(unit,*) 'fv_restart v   ', trim(gn),' = ', mpp_chksum(Atm(n)%v(isc:iec,jsc:jec,:))
      write(unit,*) 'fv_restart ua   ', trim(gn),' = ', mpp_chksum(Atm(n)%ua(isc:iec,jsc:jec,:))
      write(unit,*) 'fv_restart va   ', trim(gn),' = ', mpp_chksum(Atm(n)%va(isc:iec,jsc:jec,:))
      if ( .not.Atm(n)%flagstruct%hydrostatic )   &
        write(unit,*) 'fv_restart w   ', trim(gn),' = ', mpp_chksum(Atm(n)%w(isc:iec,jsc:jec,:))
      write(unit,*) 'fv_restart delp', trim(gn),' = ', mpp_chksum(Atm(n)%delp(isc:iec,jsc:jec,:))
      write(unit,*) 'fv_restart phis', trim(gn),' = ', mpp_chksum(Atm(n)%phis(isc:iec,jsc:jec))

#ifdef SW_DYNAMICS
      call prt_maxmin('H ', Atm(n)%delp, isc, iec, jsc, jec, Atm(n)%ng, 1, rgrav)
#else
      write(unit,*) 'fv_restart pt  ', trim(gn),' = ', mpp_chksum(Atm(n)%pt(isc:iec,jsc:jec,:))
      if (ntprog>0) &
           write(unit,*) 'fv_restart q(prog) nq  ', trim(gn),' =',ntprog, mpp_chksum(Atm(n)%q(isc:iec,jsc:jec,:,:))
      if (ntdiag>0) &
           write(unit,*) 'fv_restart q(diag) nq  ', trim(gn),' =',ntdiag, mpp_chksum(Atm(n)%qdiag(isc:iec,jsc:jec,:,:))
      do iq=1,min(17, ntprog)     ! Check up to 17 tracers
        call get_tracer_names(MODEL_ATMOS, iq, tracer_name)
        write(unit,*) 'fv_restart '//trim(tracer_name)//' = ', mpp_chksum(Atm(n)%q(isc:iec,jsc:jec,:,iq))
      enddo

!---------------
! Check Min/Max:
!---------------
      call pmaxmn_g('ZS', Atm(n)%phis, isc, iec, jsc, jec, 1, rgrav, Atm(n)%gridstruct%area_64, Atm(n)%domain)
      call pmaxmn_g('PS', Atm(n)%ps,   isc, iec, jsc, jec, 1, 0.01,  Atm(n)%gridstruct%area_64, Atm(n)%domain)
      call pmaxmn_g('T ', Atm(n)%pt,   isc, iec, jsc, jec, npz, 1.,  Atm(n)%gridstruct%area_64, Atm(n)%domain)

! Check tracers:
      do i=1, ntprog
          call get_tracer_names ( MODEL_ATMOS, i, tname )
          call pmaxmn_g(trim(tname), Atm(n)%q(isd:ied,jsd:jed,1:npz,i:i), isc, iec, jsc, jec, npz, &
                        1., Atm(n)%gridstruct%area_64, Atm(n)%domain)
      enddo
#endif
      call prt_maxmin('U (local) ', Atm(n)%u(isc:iec,jsc:jec,1:npz), isc, iec, jsc, jec, 0, npz, 1.)
      call prt_maxmin('V (local) ', Atm(n)%v(isc:iec,jsc:jec,1:npz), isc, iec, jsc, jec, 0, npz, 1.)
      ! compute ua, va
      call cubed_to_latlon(Atm(n)%u, Atm(n)%v, Atm(n)%ua, Atm(n)%va, &
           Atm(n)%gridstruct, &
           Atm(n)%npx, Atm(n)%npy, npz, 1, &
           Atm(n)%gridstruct%grid_type, Atm(n)%domain, &
           Atm(n)%gridstruct%bounded_domain, Atm(n)%flagstruct%c2l_ord, Atm(n)%bd)
      call prt_maxmin('UA ', Atm(n)%ua, isc, iec, jsc, jec, Atm(n)%ng, npz, 1.)
      call prt_maxmin('VA ', Atm(n)%va, isc, iec, jsc, jec, Atm(n)%ng, npz, 1.)

      if ( (.not.Atm(n)%flagstruct%hydrostatic) .and. Atm(n)%flagstruct%make_nh ) then
         call mpp_error(NOTE, "  Initializing w to 0")
         Atm(n)%w = 0.
         sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
         if ( .not.Atm(n)%flagstruct%hybrid_z ) then
            if (Atm(n)%flagstruct%adiabatic .or. sphum < 0) then
             zvir = 0.
            else
             zvir = rvgas/rdgas - 1.
            endif
             do k=1,npz
                do j=jsc,jec
                   do i=isc,iec
                      Atm(n)%delz(i,j,k) = (rdgas*rgrav)*Atm(n)%pt(i,j,k)*(1.+zvir*Atm(n)%q(i,j,k,sphum))*(Atm(n)%peln(i,k,j)-Atm(n)%peln(i,k+1,j))
                   enddo
                enddo
             enddo
         endif
      endif

      if ( .not.Atm(n)%flagstruct%hydrostatic )   &
      call pmaxmn_g('W ', Atm(n)%w, isc, iec, jsc, jec, npz, 1., Atm(n)%gridstruct%area_64, Atm(n)%domain)

      if (is_master()) write(unit,*)

!--------------------------------------------
! Initialize surface winds for flux coupler:
!--------------------------------------------
      if ( .not. Atm(n)%flagstruct%srf_init ) then
         call cubed_to_latlon(Atm(n)%u, Atm(n)%v, Atm(n)%ua, Atm(n)%va, &
              Atm(n)%gridstruct, &
              Atm(n)%npx, Atm(n)%npy, npz, 1, &
              Atm(n)%gridstruct%grid_type, Atm(n)%domain, &
              Atm(n)%gridstruct%bounded_domain, Atm(n)%flagstruct%c2l_ord, Atm(n)%bd)
         do j=jsc,jec
            do i=isc,iec
               Atm(n)%u_srf(i,j) = Atm(n)%ua(i,j,npz)
               Atm(n)%v_srf(i,j) = Atm(n)%va(i,j,npz)
            enddo
         enddo
         Atm(n)%flagstruct%srf_init = .true.
      endif

    end do   ! n_tile

  end subroutine fv_restart
  ! </SUBROUTINE> NAME="fv_restart"


  subroutine fill_nested_grid_topo_halo(Atm, proc_in)

    type(fv_atmos_type), intent(INOUT) :: Atm
    logical, intent(IN), OPTIONAL :: proc_in
    integer :: isd, ied, jsd, jed

    if (.not. Atm%neststruct%nested) return

    call mpp_get_data_domain( Atm%parent_grid%domain, &
         isd, ied, jsd, jed)

    !This is 2D and doesn't need remapping
    if (is_master()) print*, '  FILLING NESTED GRID HALO WITH INTERPOLATED TERRAIN'
    call nested_grid_BC(Atm%phis, Atm%parent_grid%phis, global_nest_domain, &
         Atm%neststruct%ind_h, Atm%neststruct%wt_h, 0, 0, &
         Atm%npx, Atm%npy, Atm%bd, isd, ied, jsd, jed, proc_in=proc_in, nest_level=Atm%grid_number-1)

  end subroutine fill_nested_grid_topo_halo

!!! We call this routine to fill the nested grid with topo so that we can do the boundary smoothing.
!!! Interior topography is then over-written in get_external_ic.
!!! Input grid is the nest; use Atm%parent_grid% to reference parent
  subroutine fill_nested_grid_topo(Atm, proc_in)

    type(fv_atmos_type), intent(INOUT) :: Atm
    logical, intent(IN), OPTIONAL :: proc_in
    real, allocatable :: g_dat(:,:,:)
    integer :: p, sending_proc
    integer :: isd_p, ied_p, jsd_p, jed_p
    integer :: isg, ieg, jsg,jeg

    logical :: process

    process = .true.
    if (present(proc_in)) then
       process = proc_in
    else
       process = .true.
    endif

    call mpp_get_global_domain( Atm%parent_grid%domain, &
         isg, ieg, jsg, jeg)
    call mpp_get_data_domain( Atm%parent_grid%domain, &
         isd_p,  ied_p,  jsd_p,  jed_p  )

    allocate(g_dat( isg:ieg, jsg:jeg, 1) )

    call timing_on('COMM_TOTAL')

    !!! FIXME: For whatever reason this code CRASHES if the lower-left corner
    !!!        of the nested grid lies within the first PE of a grid tile.

    if (is_master() .and. .not. Atm%flagstruct%external_ic ) print*, ' FILLING NESTED GRID INTERIOR WITH INTERPOLATED TERRAIN'

    sending_proc = (Atm%parent_grid%pelist(1)) + &
         (Atm%neststruct%parent_tile-tile_fine(Atm%parent_grid%grid_number)+Atm%parent_grid%flagstruct%ntiles-1)*Atm%parent_grid%npes_per_tile
    if (Atm%neststruct%parent_tile == Atm%parent_grid%global_tile) then
    !if (Atm%neststruct%parent_proc .and. Atm%neststruct%parent_tile == Atm%parent_grid%global_tile) then
       call mpp_global_field( &
            Atm%parent_grid%domain, &
            Atm%parent_grid%phis(isd_p:ied_p,jsd_p:jed_p), g_dat(isg:,jsg:,1), position=CENTER)
       if (mpp_pe() == sending_proc) then
          do p=1,size(Atm%pelist)
             call mpp_send(g_dat,size(g_dat),Atm%pelist(p))
          enddo
       endif
    endif

    if (ANY(Atm%pelist == mpp_pe())) then
       call mpp_recv(g_dat, size(g_dat), sending_proc)
    endif

    call timing_off('COMM_TOTAL')

    if (process) call fill_nested_grid(Atm%phis, g_dat(isg:,jsg:,1), &
         Atm%neststruct%ind_h, Atm%neststruct%wt_h, &
         0, 0,  isg, ieg, jsg, jeg, Atm%bd)

    call mpp_sync_self

    deallocate(g_dat)


  end subroutine fill_nested_grid_topo

  !This will still probably be needed for moving nests
  !NOTE: this has NOT been maintained and so %global_tile is now meaningless if not referring to data on the current PE
  !      needs to be re-coded to follow method in fill_nested_grid_Topo
  subroutine fill_nested_grid_data(Atm, proc_in)

    type(fv_atmos_type), intent(INOUT) :: Atm(:) !Only intended to be one element; needed for cubed_sphere_terrain
    logical, intent(IN), OPTIONAL :: proc_in
    real, allocatable :: g_dat(:,:,:), pt_coarse(:,:,:)
    integer :: i,j,k,nq, sphum, ncnst, istart, iend, npz
    integer :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer :: isd_p, ied_p, jsd_p, jed_p, isc_p, iec_p, jsc_p, jec_p
    integer :: isg, ieg, jsg,jeg, npx_p, npy_p
    integer :: isg_n, ieg_n, jsg_n, jeg_n, npx_n, npy_n
    real zvir, gh0, p1(2), p2(2), r, r0

    integer :: p, sending_proc, gid
    logical process

    call mpp_error(FATAL, " FILL_NESTED_GRID_DATA not yet updated for remap BCs")

    if (present(proc_in)) then
       process = proc_in
    else
       process = .true.
    endif

    isd = Atm(1)%bd%isd
    ied = Atm(1)%bd%ied
    jsd = Atm(1)%bd%jsd
    jed = Atm(1)%bd%jed
    ncnst = Atm(1)%ncnst
    isc = Atm(1)%bd%isc; iec = Atm(1)%bd%iec; jsc = Atm(1)%bd%jsc; jec = Atm(1)%bd%jec
    npz     = Atm(1)%npz

    gid = mpp_pe()

    sending_proc = Atm(1)%parent_grid%pelist(1) + (Atm(1)%neststruct%parent_tile-1)*Atm(1)%parent_grid%npes_per_tile

       call mpp_get_data_domain( Atm(1)%parent_grid%domain, &
            isd_p,  ied_p,  jsd_p,  jed_p  )
       call mpp_get_compute_domain( Atm(1)%parent_grid%domain, &
            isc_p,  iec_p,  jsc_p,  jec_p  )
    call mpp_get_global_domain( Atm(1)%parent_grid%domain, &
         isg, ieg, jsg, jeg, xsize=npx_p, ysize=npy_p)

    if (process) then

       call mpp_error(NOTE, "FILLING NESTED GRID DATA")

    else

       call mpp_error(NOTE, "SENDING TO FILL NESTED GRID DATA")

    endif

    !delp

    allocate(g_dat( isg:ieg, jsg:jeg, npz) )

    call timing_on('COMM_TOTAL')

    !Call mpp_global_field on the procs that have the required data.
    !Then broadcast from the head PE to the receiving PEs
    if (Atm(1)%neststruct%parent_proc .and. Atm(1)%neststruct%parent_tile == Atm(1)%parent_grid%global_tile) then
       call mpp_global_field( &
            Atm(1)%parent_grid%domain, &
            Atm(1)%parent_grid%delp(isd_p:ied_p,jsd_p:jed_p,:), g_dat, position=CENTER)
       if (gid == sending_proc) then !crazy logic but what we have for now
          do p=1,size(Atm(1)%pelist)
             call mpp_send(g_dat,size(g_dat),Atm(1)%pelist(p))
          enddo
       endif
    endif
    if (ANY(Atm(1)%pelist == gid)) then
       call mpp_recv(g_dat, size(g_dat), sending_proc)
    endif

    call timing_off('COMM_TOTAL')

    if (process) call fill_nested_grid(Atm(1)%delp, g_dat, &
         Atm(1)%neststruct%ind_h, Atm(1)%neststruct%wt_h, &
         0, 0,  isg, ieg, jsg, jeg, npz, Atm(1)%bd)

    call mpp_sync_self

    !tracers
    do nq=1,ncnst

       call timing_on('COMM_TOTAL')

       if (ANY(Atm(1)%parent_grid%pelist == gid) .and. Atm(1)%neststruct%parent_tile == Atm(1)%parent_grid%global_tile) then
          call mpp_global_field( &
            Atm(1)%parent_grid%domain, &
            Atm(1)%parent_grid%q(isd_p:ied_p,jsd_p:jed_p,:,nq), g_dat, position=CENTER)
          if (gid == sending_proc) then
             do p=1,size(Atm(1)%pelist)
                call mpp_send(g_dat,size(g_dat),Atm(1)%pelist(p))
             enddo
          endif
       endif
       if (ANY(Atm(1)%pelist == gid)) then
          call mpp_recv(g_dat, size(g_dat), sending_proc)
       endif

       call timing_off('COMM_TOTAL')

       if (process) call fill_nested_grid(Atm(1)%q(isd:ied,jsd:jed,:,nq), g_dat, &
            Atm(1)%neststruct%ind_h, Atm(1)%neststruct%wt_h, &
            0, 0,  isg, ieg, jsg, jeg, npz, Atm(1)%bd)

       call mpp_sync_self

    end do

    !Note that we do NOT fill in phis (surface geopotential), which should
    !be computed exactly instead of being interpolated.


#ifndef SW_DYNAMICS
    !pt --- actually temperature

    call timing_on('COMM_TOTAL')

    if (ANY(Atm(1)%parent_grid%pelist == gid) .and. Atm(1)%neststruct%parent_tile == Atm(1)%parent_grid%global_tile) then
          call mpp_global_field( &
            Atm(1)%parent_grid%domain, &
            Atm(1)%parent_grid%pt(isd_p:ied_p,jsd_p:jed_p,:), g_dat, position=CENTER)
       if (gid == sending_proc) then
          do p=1,size(Atm(1)%pelist)
             call mpp_send(g_dat,size(g_dat),Atm(1)%pelist(p))
          enddo
       endif
    endif
    if (ANY(Atm(1)%pelist == gid)) then
       call mpp_recv(g_dat, size(g_dat), sending_proc)
    endif

    call mpp_sync_self

    call timing_off('COMM_TOTAL')

    if (process) call fill_nested_grid(Atm(1)%pt, g_dat, &
         Atm(1)%neststruct%ind_h, Atm(1)%neststruct%wt_h, &
         0, 0,  isg, ieg, jsg, jeg, npz, Atm(1)%bd)


    if ( Atm(1)%flagstruct%nwat > 0 ) then
       sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
    else
       sphum = 1
    endif
    if ( Atm(1)%parent_grid%flagstruct%adiabatic .or. Atm(1)%parent_grid%flagstruct%do_Held_Suarez ) then
       zvir = 0.         ! no virtual effect
    else
       zvir = rvgas/rdgas - 1.
    endif

    call timing_on('COMM_TOTAL')

    if (ANY(Atm(1)%parent_grid%pelist == gid) .and. Atm(1)%neststruct%parent_tile == Atm(1)%parent_grid%global_tile) then
          call mpp_global_field( &
            Atm(1)%parent_grid%domain, &
            Atm(1)%parent_grid%pkz(isc_p:iec_p,jsc_p:jec_p,:), g_dat, position=CENTER)
       if (gid == sending_proc) then
          do p=1,size(Atm(1)%pelist)
             call mpp_send(g_dat,size(g_dat),Atm(1)%pelist(p))
          enddo
       endif
    endif
    if (ANY(Atm(1)%pelist == gid)) then
       call mpp_recv(g_dat, size(g_dat), sending_proc)
    endif

    call mpp_sync_self

    call timing_off('COMM_TOTAL')

    if (process) then
       allocate(pt_coarse(isd:ied,jsd:jed,npz))
       call fill_nested_grid(pt_coarse, g_dat, &
            Atm(1)%neststruct%ind_h, Atm(1)%neststruct%wt_h, &
            0, 0,  isg, ieg, jsg, jeg, npz, Atm(1)%bd)

       if (Atm(1)%bd%is == 1) then
          do k=1,npz
             do j=Atm(1)%bd%jsd,Atm(1)%bd%jed
                do i=Atm(1)%bd%isd,0
                   Atm(1)%pt(i,j,k) = cp_air*Atm(1)%pt(i,j,k)/pt_coarse(i,j,k)*(1.+zvir*Atm(1)%q(i,j,k,sphum))
                end do
             end do
          end do
       end if

       if (Atm(1)%bd%js == 1) then
          if (Atm(1)%bd%is == 1) then
             istart = Atm(1)%bd%is
          else
             istart = Atm(1)%bd%isd
          end if
          if (Atm(1)%bd%ie == Atm(1)%npx-1) then
             iend = Atm(1)%bd%ie
          else
             iend = Atm(1)%bd%ied
          end if

          do k=1,npz
             do j=Atm(1)%bd%jsd,0
                do i=istart,iend
                   Atm(1)%pt(i,j,k) = cp_air*Atm(1)%pt(i,j,k)/pt_coarse(i,j,k)*(1.+zvir*Atm(1)%q(i,j,k,sphum))
                end do
             end do
          end do
       end if

       if (Atm(1)%bd%ie == Atm(1)%npx-1) then
          do k=1,npz
             do j=Atm(1)%bd%jsd,Atm(1)%bd%jed
                do i=Atm(1)%npx,Atm(1)%bd%ied
                   Atm(1)%pt(i,j,k) = cp_air*Atm(1)%pt(i,j,k)/pt_coarse(i,j,k)*(1.+zvir*Atm(1)%q(i,j,k,sphum))
                end do
             end do
          end do
       end if

       if (Atm(1)%bd%je == Atm(1)%npy-1) then
          if (Atm(1)%bd%is == 1) then
             istart = Atm(1)%bd%is
          else
             istart = Atm(1)%bd%isd
          end if
          if (Atm(1)%bd%ie == Atm(1)%npx-1) then
             iend = Atm(1)%bd%ie
          else
             iend = Atm(1)%bd%ied
          end if

          do k=1,npz
             do j=Atm(1)%npy,Atm(1)%bd%jed
                do i=istart,iend
                   Atm(1)%pt(i,j,k) = cp_air*Atm(1)%pt(i,j,k)/pt_coarse(i,j,k)*(1.+zvir*Atm(1)%q(i,j,k,sphum))
                end do
             end do
          end do
       end if

       deallocate(pt_coarse)

    end if

    if (.not. Atm(1)%flagstruct%hydrostatic) then

       !delz
       call timing_on('COMM_TOTAL')

       if (ANY(Atm(1)%parent_grid%pelist == gid) .and. Atm(1)%neststruct%parent_tile == Atm(1)%parent_grid%global_tile) then
          call mpp_global_field( &
            Atm(1)%parent_grid%domain, &
            Atm(1)%parent_grid%delz(isd_p:ied_p,jsd_p:jed_p,:), g_dat, position=CENTER)
          if (gid == sending_proc) then
             do p=1,size(Atm(1)%pelist)
                call mpp_send(g_dat,size(g_dat),Atm(1)%pelist(p))
             enddo
          endif
       endif
       if (ANY(Atm(1)%pelist == gid)) then
          call mpp_recv(g_dat, size(g_dat), sending_proc)
       endif

       call mpp_sync_self

       call timing_off('COMM_TOTAL')

       if (process) call fill_nested_grid(Atm(1)%delz, g_dat, &
            Atm(1)%neststruct%ind_h, Atm(1)%neststruct%wt_h, &
            0, 0,  isg, ieg, jsg, jeg, npz, Atm(1)%bd)

       !w

       call timing_on('COMM_TOTAL')

       if (ANY(Atm(1)%parent_grid%pelist == gid) .and. Atm(1)%neststruct%parent_tile == Atm(1)%parent_grid%global_tile) then
          call mpp_global_field( &
            Atm(1)%parent_grid%domain, &
            Atm(1)%parent_grid%w(isd_p:ied_p,jsd_p:jed_p,:), g_dat, position=CENTER)
          if (gid == sending_proc) then
             do p=1,size(Atm(1)%pelist)
                call mpp_send(g_dat,size(g_dat),Atm(1)%pelist(p))
             enddo
          endif
       endif
       if (ANY(Atm(1)%pelist == gid)) then
          call mpp_recv(g_dat, size(g_dat), sending_proc)
       endif

       call mpp_sync_self

       call timing_off('COMM_TOTAL')

       if (process) call fill_nested_grid(Atm(1)%w, g_dat, &
            Atm(1)%neststruct%ind_h, Atm(1)%neststruct%wt_h, &
            0, 0,  isg, ieg, jsg, jeg, npz, Atm(1)%bd)
       !

    end if

#endif
    deallocate(g_dat)

    !u

    allocate(g_dat( isg:ieg, jsg:jeg+1, npz) )
    g_dat = 1.e25

    call timing_on('COMM_TOTAL')

    if (ANY(Atm(1)%parent_grid%pelist == gid) .and. Atm(1)%neststruct%parent_tile == Atm(1)%parent_grid%global_tile) then
          call mpp_global_field( &
            Atm(1)%parent_grid%domain, &
            Atm(1)%parent_grid%u(isd_p:ied_p,jsd_p:jed_p+1,:), g_dat, position=NORTH)
       if (gid == sending_proc) then
          do p=1,size(Atm(1)%pelist)
             call mpp_send(g_dat,size(g_dat),Atm(1)%pelist(p))
          enddo
       endif
    endif
    if (ANY(Atm(1)%pelist == gid)) then
       call mpp_recv(g_dat, size(g_dat), sending_proc)
    endif

    call mpp_sync_self

    call timing_off('COMM_TOTAL')

    call mpp_sync_self
    if (process) call fill_nested_grid(Atm(1)%u, g_dat, &
         Atm(1)%neststruct%ind_u, Atm(1)%neststruct%wt_u, &
         0, 1,  isg, ieg, jsg, jeg, npz, Atm(1)%bd)
    deallocate(g_dat)

    !v

    allocate(g_dat( isg:ieg+1, jsg:jeg, npz) )
    g_dat = 1.e25

    call timing_on('COMM_TOTAL')

    if (ANY(Atm(1)%parent_grid%pelist == gid) .and. Atm(1)%neststruct%parent_tile == Atm(1)%parent_grid%global_tile) then
          call mpp_global_field( &
            Atm(1)%parent_grid%domain, &
            Atm(1)%parent_grid%v(isd_p:ied_p+1,jsd_p:jed_p,:), g_dat, position=EAST)
       if (gid == sending_proc) then
          do p=1,size(Atm(1)%pelist)
             call mpp_send(g_dat,size(g_dat),Atm(1)%pelist(p))
          enddo
       endif
    endif
    if (ANY(Atm(1)%pelist == gid)) then
       call mpp_recv(g_dat, size(g_dat), sending_proc)
    endif

    call mpp_sync_self

    call timing_off('COMM_TOTAL')

    if (process) call fill_nested_grid(Atm(1)%v, g_dat, &
         Atm(1)%neststruct%ind_v, Atm(1)%neststruct%wt_v, &
         1, 0,  isg, ieg, jsg, jeg, npz, Atm(1)%bd)

    deallocate(g_dat)

  end subroutine fill_nested_grid_data

  !This routine actually sets up the coarse-grid TOPOGRAPHY.
  subroutine twoway_topo_update(Atm, proc_in)

    type(fv_atmos_type), intent(INOUT) :: Atm
    logical, intent(IN), OPTIONAL :: proc_in
    integer :: i,j,k,nq, sphum, ncnst, istart, iend, npz
    integer :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer :: isd_p, ied_p, jsd_p, jed_p, isc_p, iec_p, jsc_p, jec_p
    integer :: isg, ieg, jsg,jeg, npx_p, npy_p
    real zvir

    integer :: p , sending_proc
    logical :: process

    if (present(proc_in)) then
       process = proc_in
    else
       process = .true.
    endif

    isd = Atm%bd%isd
    ied = Atm%bd%ied
    jsd = Atm%bd%jsd
    jed = Atm%bd%jed
    ncnst = Atm%ncnst
    isc = Atm%bd%isc; iec = Atm%bd%iec; jsc = Atm%bd%jsc; jec = Atm%bd%jec
    npz     = Atm%npz

    isd_p = Atm%parent_grid%bd%isd
    ied_p = Atm%parent_grid%bd%ied
    jsd_p = Atm%parent_grid%bd%jsd
    jed_p = Atm%parent_grid%bd%jed
    isc_p = Atm%parent_grid%bd%isc
    iec_p = Atm%parent_grid%bd%iec
    jsc_p = Atm%parent_grid%bd%jsc
    jec_p = Atm%parent_grid%bd%jec
    !sending_proc = Atm%parent_grid%pelist(1) + (Atm%neststruct%parent_tile-1)*Atm%parent_grid%npes_per_tile

    call mpp_get_global_domain( Atm%parent_grid%domain, &
         isg, ieg, jsg, jeg, xsize=npx_p, ysize=npy_p)


    !NOW: what we do is to update the nested-grid terrain to the coarse grid,
    !to ensure consistency between the two grids.
    if ( process ) call mpp_update_domains(Atm%phis, Atm%domain, complete=.true.)
 !   if (Atm%neststruct%twowaynest) then
       if (ANY(Atm%parent_grid%pelist == mpp_pe()) .or. Atm%neststruct%child_proc) then
          call update_coarse_grid(Atm%parent_grid%phis, &
               Atm%phis, global_nest_domain, &
               Atm%gridstruct%dx, Atm%gridstruct%dy, Atm%gridstruct%area, &
               Atm%bd, isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, &
               Atm%neststruct%isu, Atm%neststruct%ieu, Atm%neststruct%jsu, Atm%neststruct%jeu, &
               Atm%npx, Atm%npy, 0, 0, &
               Atm%neststruct%refinement, Atm%neststruct%nestupdate, 0, 0, &
               ANY(Atm%parent_grid%pelist == mpp_pe()), Atm%neststruct%child_proc, Atm%parent_grid, Atm%neststruct%nlevel)
          Atm%parent_grid%neststruct%parent_of_twoway = .true.
          !NOTE: mpp_update_nest_coarse (and by extension, update_coarse_grid) does **NOT** pass data
          !allowing a two-way update into the halo of the coarse grid. It only passes data so that the INTERIOR
          ! can have the two-way update. Thus, on the nest's cold start, if this update_domains call is not done,
          ! the coarse grid will have the wrong topography in the halo, which will CHANGE when a restart is done!!
         !if (Atm%neststruct%parent_proc) call mpp_update_domains(Atm%parent_grid%phis, Atm%parent_grid%domain)
         if (ANY(Atm%parent_grid%pelist == mpp_pe())) call mpp_update_domains(Atm%parent_grid%phis, Atm%parent_grid%domain)
       end if

  !  end if


#ifdef SW_DYNAMICS
#else
    !Reset p_var after updating topography
    if (process) call p_var(npz, isc, iec, jsc, jec, Atm%ptop, ptop_min, Atm%delp, &
         Atm%delz, Atm%pt, Atm%ps,   &
         Atm%pe, Atm%peln, Atm%pk, Atm%pkz, kappa, Atm%q, &
         Atm%ng, ncnst, Atm%gridstruct%area_64, Atm%flagstruct%dry_mass, .false., Atm%flagstruct%mountain, &
         Atm%flagstruct%moist_phys, .true., Atm%flagstruct%nwat, Atm%domain, Atm%flagstruct%adiabatic)
#endif



  end subroutine twoway_topo_update


  !#######################################################################
  ! <SUBROUTINE NAME="fv_write_restart">
  ! <DESCRIPTION>
  !  Write out restart files registered through register_restart_file
  ! </DESCRIPTION>
  subroutine fv_write_restart(Atm, timestamp)
    type(fv_atmos_type), intent(inout) :: Atm
    character(len=*),    intent(in)    :: timestamp

    if (Atm%coarse_graining%write_coarse_restart_files) then
       call fv_io_write_restart_coarse(Atm, timestamp)
       if (.not. Atm%coarse_graining%write_only_coarse_intermediate_restarts) then
          call fv_io_write_restart(Atm, prefix=timestamp)
       endif
    else
       call fv_io_write_restart(Atm, prefix=timestamp)
    endif

    if (Atm%neststruct%nested) then
       call fv_io_write_BCs(Atm)
    endif

  end subroutine fv_write_restart
  ! </SUBROUTINE>



  !#####################################################################
  ! <SUBROUTINE NAME="fv_restart_end">
  !
  ! <DESCRIPTION>
  ! Initialize the fv core restart facilities
  ! </DESCRIPTION>
  !
  subroutine fv_restart_end(Atm)
    type(fv_atmos_type), intent(inout) :: Atm

    integer :: isc, iec, jsc, jec
    integer :: iq, ncnst, ntprog, ntdiag
    integer :: isd, ied, jsd, jed, npz
    integer :: unit
    integer :: file_unit
    integer, allocatable :: pelist(:)
    character(len=128):: tracer_name
    character(len=3):: gn


    call mpp_set_current_pelist(Atm%pelist)

    isc = Atm%bd%isc; iec = Atm%bd%iec; jsc = Atm%bd%jsc; jec = Atm%bd%jec

    isd = Atm%bd%isd
    ied = Atm%bd%ied
    jsd = Atm%bd%jsd
    jed = Atm%bd%jed
    npz = Atm%npz
    ncnst = Atm%ncnst
    ntprog = size(Atm%q,4)
    ntdiag = size(Atm%qdiag,4)

    if (Atm%grid_number > 1) then
       write(gn,'(A2, I1)') " g", Atm%grid_number
    else
       gn = ''
    end if

    unit = stdout()
    write(unit,*)
    write(unit,*) 'fv_restart_end u   ', trim(gn),' = ', mpp_chksum(Atm%u(isc:iec,jsc:jec,:))
    write(unit,*) 'fv_restart_end v   ', trim(gn),' = ', mpp_chksum(Atm%v(isc:iec,jsc:jec,:))
    write(unit,*) 'fv_restart_end ua   ', trim(gn),' = ', mpp_chksum(Atm%ua(isc:iec,jsc:jec,:))
    write(unit,*) 'fv_restart_end va   ', trim(gn),' = ', mpp_chksum(Atm%va(isc:iec,jsc:jec,:))
    if ( .not. Atm%flagstruct%hydrostatic )    &
         write(unit,*) 'fv_restart_end w   ', trim(gn),' = ', mpp_chksum(Atm%w(isc:iec,jsc:jec,:))
    write(unit,*) 'fv_restart_end delp', trim(gn),' = ', mpp_chksum(Atm%delp(isc:iec,jsc:jec,:))
    write(unit,*) 'fv_restart_end phis', trim(gn),' = ', mpp_chksum(Atm%phis(isc:iec,jsc:jec))
#ifndef SW_DYNAMICS
    write(unit,*) 'fv_restart_end pt  ', trim(gn),' = ', mpp_chksum(Atm%pt(isc:iec,jsc:jec,:))
    if (ntprog>0) &
         write(unit,*) 'fv_restart_end q(prog) nq  ', trim(gn),' =',ntprog, mpp_chksum(Atm%q(isc:iec,jsc:jec,:,:))
    if (ntdiag>0) &
         write(unit,*) 'fv_restart_end q(diag) nq  ', trim(gn),' =',ntdiag, mpp_chksum(Atm%qdiag(isc:iec,jsc:jec,:,:))
    do iq=1,min(17, ntprog)     ! Check up to 17 tracers
       call get_tracer_names(MODEL_ATMOS, iq, tracer_name)
       write(unit,*) 'fv_restart_end '//trim(tracer_name)// trim(gn),' = ', mpp_chksum(Atm%q(isc:iec,jsc:jec,:,iq))
    enddo

    !---------------
    ! Check Min/Max:
    !---------------
    !     call prt_maxmin('ZS', Atm%phis, isc, iec, jsc, jec, Atm%ng, 1, 1./grav)
    call pmaxmn_g('ZS', Atm%phis, isc, iec, jsc, jec, 1, 1./grav, Atm%gridstruct%area_64, Atm%domain)
    call pmaxmn_g('PS ', Atm%ps,   isc, iec, jsc, jec, 1, 0.01   , Atm%gridstruct%area_64, Atm%domain)
    call prt_maxmin('PS*', Atm%ps, isc, iec, jsc, jec, Atm%ng, 1, 0.01)
    call prt_maxmin('U (local) ', Atm%u(isc:iec,jsc:jec,1:npz), isc, iec, jsc, jec, 0, npz, 1.)
    call prt_maxmin('V (local) ', Atm%v(isc:iec,jsc:jec,1:npz), isc, iec, jsc, jec, 0, npz, 1.)
    call prt_maxmin('UA ', Atm%ua, isc, iec, jsc, jec, Atm%ng, npz, 1.)
    call prt_maxmin('VA ', Atm%va, isc, iec, jsc, jec, Atm%ng, npz, 1.)
    if ( .not. Atm%flagstruct%hydrostatic )    &
         call prt_maxmin('W ', Atm%w , isc, iec, jsc, jec, Atm%ng, npz, 1.)
    call prt_maxmin('T ', Atm%pt, isc, iec, jsc, jec, Atm%ng, npz, 1.)
    do iq=1, ntprog
       call get_tracer_names ( MODEL_ATMOS, iq, tracer_name )
       call pmaxmn_g(trim(tracer_name), Atm%q(isd:ied,jsd:jed,1:npz,iq:iq), isc, iec, jsc, jec, npz, &
            1., Atm%gridstruct%area_64, Atm%domain)
    enddo
    ! Write4 energy correction term
#endif

       call fv_io_write_restart(Atm)
       if (Atm%coarse_graining%write_coarse_restart_files) then
          call fv_io_write_restart_coarse(Atm)
       endif

 if (Atm%neststruct%nested) call fv_io_write_BCs(Atm)

 if (Atm%flagstruct%write_restart_with_bcs) call write_full_fields(Atm)

 module_is_initialized = .FALSE.

#ifdef EFLUX_OUT
 if( is_master() ) then
    write(*,*) steps, 'Mean equivalent Heat flux for this integration period=',Atm(1)%idiag%efx_sum/real(max(1,Atm(1)%idiag%steps)), &
         'Mean nesting-related flux for this integration period=',Atm(1)%idiag%efx_sum_nest/real(max(1,Atm(1)%idiag%steps)), &
         'Mean mountain torque=',Atm(1)%idiag%mtq_sum/real(max(1,Atm(1)%idiag%steps))
    file_unit = get_unit()
    open (unit=file_unit, file='e_flux.data', form='unformatted',status='unknown', access='sequential')
    do n=1,steps
       write(file_unit) Atm(1)%idiag%efx(n)
       write(file_unit) Atm(1)%idiag%mtq(n)    ! time series global mountain torque
       !write(file_unit) Atm(1)%idiag%efx_nest(n)
    enddo
    close(unit=file_unit)
 endif
#endif

  end subroutine fv_restart_end
  ! </SUBROUTINE> NAME="fv_restart_end"


subroutine pmaxmn_g(qname, q, is, ie, js, je, km, fac, area, domain)
      character(len=*), intent(in)::  qname
      integer, intent(in):: is, ie, js, je
      integer, intent(in):: km
      real, intent(in)::    q(is-3:ie+3, js-3:je+3, km)
      real, intent(in)::    fac
      real(kind=R_GRID), intent(IN)::    area(is-3:ie+3, js-3:je+3)
      type(domain2d), intent(INOUT) :: domain
!
      real qmin, qmax, gmean
      integer i,j,k
      character(len=8) :: display_name

      logical, SAVE :: first_time = .true.

      qmin = q(is,js,1)
      qmax = qmin

      do k=1,km
      do j=js,je
         do i=is,ie
            !if ( (q(i,j,k) >= 1e30) .eqv. (q(i,j,k) < 1e30) ) then !NAN checking
            !   print*, ' NAN found for ', qname, mpp_pe(), i,j,k
            !else
            if( q(i,j,k) < qmin) then
                qmin = q(i,j,k)
            elseif( q(i,j,k) > qmax ) then
                qmax = q(i,j,k)
            endif
          enddo
      enddo
      enddo

      call mp_reduce_min(qmin)
      call mp_reduce_max(qmax)

      gmean = g_sum(domain, q(is:ie,js:je,km), is, ie, js, je, 3, area, 1, .true.)

      if(is_master()) then
         j = min(len(trim(qname)),8)
         display_name = qname(1:j)
         write(6,*) display_name, trim(gn), qmax*fac, qmin*fac, gmean*fac
      endif

end subroutine pmaxmn_g

end module fv_restart_mod
