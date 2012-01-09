module atmosphere_mod

!-----------------------------------------------------------------------
!
! Interface for Cubed_Sphere fv dynamical core and physics
!
!-----------------------------------------------------------------------

!-----------------
! FMS modules:
!-----------------
use constants_mod,      only: cp_air, rdgas, grav, rvgas, kappa, pstd_mks
use time_manager_mod,   only: time_type, get_time, set_time, operator(+)
use fms_mod,            only: file_exist, open_namelist_file,    &
                              close_file, error_mesg, FATAL,     &
                              check_nml_error, stdlog,           &
                              write_version_number,              &
                              mpp_pe, mpp_root_pe, set_domain,   &
                              mpp_clock_id, mpp_clock_begin,     &
                              mpp_clock_end, CLOCK_SUBCOMPONENT, &
                              clock_flag_default, nullify_domain
use mpp_mod,            only: mpp_error, FATAL, NOTE, input_nml_file
use mpp_domains_mod,    only: domain2d
use xgrid_mod,          only: grid_box_type
!miz
use diag_manager_mod,   only: diag_axis_init, register_diag_field, &
                              register_static_field, send_data
!miz
use field_manager_mod,  only: MODEL_ATMOS
use tracer_manager_mod, only: get_tracer_index,&
                              get_number_tracers, &
                              get_tracer_names

!-----------------
! FV core modules:
!-----------------
use fv_grid_tools_mod,  only: area, grid_type, dx, dy, area
use fv_grid_utils_mod,  only: edge_w, edge_e, edge_s, edge_n, en1, en2, vlon, vlat
use fv_arrays_mod,      only: fv_atmos_type
use fv_control_mod,     only: fv_init, domain, fv_end, p_ref
use fv_io_mod,          only: fv_io_register_nudge_restart
use fv_mp_mod,          only: domain_for_coupler
use fv_dynamics_mod,    only: fv_dynamics
use fv_diagnostics_mod, only: fv_diag_init, fv_diag, fv_time
use fv_restart_mod,     only: fv_restart, fv_write_restart
use fv_timing_mod,      only: timing_on, timing_off
use fv_physics_mod,     only: fv_physics_down, fv_physics_up,  &
                              fv_physics_init, fv_physics_end, &
                              surf_diff_type, fv_physics_restart
#if defined (ATMOS_NUDGE)
use atmos_nudge_mod,      only: atmos_nudge_init, atmos_nudge_end
#elif defined (CLIMATE_NUDGE)
use fv_climate_nudge_mod, only: fv_climate_nudge_init,fv_climate_nudge_end
#else
use fv_nwp_nudge_mod,     only: fv_nwp_nudge_init, fv_nwp_nudge_end
#endif

implicit none
private

public  atmosphere_down,       atmosphere_up,       &
        atmosphere_init,       atmosphere_end,      &
        atmosphere_resolution, atmosphere_boundary, &
        get_atmosphere_axes,   atmosphere_domain,   &
        get_bottom_mass,       get_bottom_wind,     &
        atmosphere_cell_area,  atmosphere_restart,  &
        get_stock_pe,          surf_diff_type

!-----------------------------------------------------------------------

character(len=128) :: version = '$Id: atmosphere.F90,v 19.0 2012/01/06 19:55:50 fms Exp $'
character(len=128) :: tagname = '$Name: siena $'
character(len=7)   :: mod_name = 'atmos'

!---- namelist (saved in file input.nml) ----
!
! physics_window  The number of "i" by "j" rows processed each time
!                 the modular physics is called. To process the entire
!                 domain use physics_window = (/0,0/).
!                   [integer, default: physics_window = 0,0]

  integer, dimension(2) :: physics_window = (/0,0/)
  namelist /atmosphere_nml/ physics_window

!---- private data ----
  type (time_type) :: Time_step_atmos
  type (fv_atmos_type), allocatable :: Atm(:)
  public Atm

  real    :: dt_atmos
  real    :: zvir
  integer :: npx, npy, npz, ncnst, pnats
  integer :: isc, iec, jsc, jec
  integer :: nq                       ! transported tracers
  integer :: sec, seconds, days
  integer :: atmos_axes(4)
  integer :: ntiles=1
  integer :: id_dynam, id_phys_down, id_phys_up, id_fv_diag
  logical :: cold_start = .false.       ! read in initial condition

  integer, dimension(:), allocatable :: id_tracerdt_dyn
  integer :: num_tracers = 0
!miz
  integer :: id_tdt_dyn, id_qdt_dyn, id_qldt_dyn, id_qidt_dyn, id_qadt_dyn
  logical :: used
  character(len=64) :: field
  real, allocatable :: ttend(:,:,:)
  real, allocatable :: qtendyyf(:,:,:,:)
  real, allocatable :: qtend(:,:,:,:)
  real              :: mv = -1.e10
!miz

contains



 subroutine atmosphere_init (Time_init, Time, Time_step, Surf_diff, Grid_box)
   type (time_type),     intent(in)    :: Time_init, Time, Time_step
   type(surf_diff_type), intent(inout) :: Surf_diff
   type(grid_box_type),  intent(inout) :: Grid_box

   integer :: unit, ierr, io, i
   integer :: itrac
   logical :: do_atmos_nudge
   character(len=32) :: tracer_name, tracer_units

   call get_number_tracers(MODEL_ATMOS, num_prog= num_tracers)

   zvir = rvgas/rdgas - 1.

!----- read namelist -----
#ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=atmosphere_nml, iostat=io)
   ierr = check_nml_error (io, 'atmosphere_nml')
#else
   if ( file_exist('input.nml') ) then
       unit = open_namelist_file ( )
       ierr=1
       do while (ierr /= 0)
          read (unit, nml=atmosphere_nml, iostat=io, end=10)
          ierr = check_nml_error (io, 'atmosphere_nml')
       enddo
 10    call close_file (unit)
   endif
#endif

!----- write version and namelist to log file -----

   unit = stdlog()
   call write_version_number ( version, tagname )
   if ( mpp_pe() == mpp_root_pe() ) write (unit, nml=atmosphere_nml)

!---- compute physics/atmos time step in seconds ----

   Time_step_atmos = Time_step
   call get_time (Time_step_atmos, sec)
   dt_atmos = real(sec)

!----- initialize FV dynamical core -----
   cold_start = (.not.file_exist('INPUT/fv_core.res.nc'))

   allocate( Atm(ntiles) )

   call fv_init( Atm(:), dt_atmos )  ! allocates Atm components

   npx   = Atm(1)%npx
   npy   = Atm(1)%npy
   npz   = Atm(1)%npz
   ncnst = Atm(1)%ncnst
   pnats = Atm(1)%pnats

   isc = Atm(1)%isc
   iec = Atm(1)%iec
   jsc = Atm(1)%jsc
   jec = Atm(1)%jec

   ! Allocate grid variables to be used to calculate gradient in 2nd order flux exchange
   allocate(Grid_box%dx    (   isc:iec  , jsc:jec+1))
   allocate(Grid_box%dy    (   isc:iec+1, jsc:jec  ))
   allocate(Grid_box%area  (   isc:iec  , jsc:jec  ))
   allocate(Grid_box%edge_w(              jsc:jec+1))
   allocate(Grid_box%edge_e(              jsc:jec+1))
   allocate(Grid_box%edge_s(   isc:iec+1           ))
   allocate(Grid_box%edge_n(   isc:iec+1           ))
   allocate(Grid_box%en1   (3, isc:iec  , jsc:jec+1))
   allocate(Grid_box%en2   (3, isc:iec+1, jsc:jec  ))
   allocate(Grid_box%vlon  (3, isc:iec  , jsc:jec  ))
   allocate(Grid_box%vlat  (3, isc:iec  , jsc:jec  ))
   Grid_box%dx    (   isc:iec  , jsc:jec+1) = dx    (   isc:iec,   jsc:jec+1)
   Grid_box%dy    (   isc:iec+1, jsc:jec  ) = dy    (   isc:iec+1, jsc:jec  )
   Grid_box%area  (   isc:iec  , jsc:jec  ) = area  (   isc:iec  , jsc:jec  )
   Grid_box%edge_w(              jsc:jec+1) = edge_w(              jsc:jec+1)
   Grid_box%edge_e(              jsc:jec+1) = edge_e(              jsc:jec+1)
   Grid_box%edge_s(   isc:iec+1           ) = edge_s(   isc:iec+1)
   Grid_box%edge_n(   isc:iec+1           ) = edge_n(   isc:iec+1)
   Grid_box%en1   (:, isc:iec  , jsc:jec+1) = en1   (:, isc:iec  , jsc:jec+1)
   Grid_box%en2   (:, isc:iec+1, jsc:jec  ) = en2   (:, isc:iec+1, jsc:jec  )
   if (allocated(vlon) .and. allocated(vlat)) then
      do i = 1, 3
         Grid_box%vlon  (i, isc:iec  , jsc:jec  ) = vlon  (isc:iec ,  jsc:jec, i  )
         Grid_box%vlat  (i, isc:iec  , jsc:jec  ) = vlat  (isc:iec ,  jsc:jec, i  )
      end do
   else
      do i = 1, 3
         Grid_box%vlon  (i, isc:iec  , jsc:jec  ) = 0.
         Grid_box%vlat  (i, isc:iec  , jsc:jec  ) = 0.
      end do
   endif
   nq = ncnst-pnats

   call fv_restart(domain, Atm, dt_atmos, seconds, days, cold_start, grid_type)

   fv_time = Time

!----- initialize atmos_axes and fv_dynamics diagnostics

   call fv_diag_init(Atm, atmos_axes, Time, npx, npy, npz, p_ref)

!----- initialize physics interface -----
!----- initialize domains for reading global physics data -----

   call set_domain ( domain )

   call fv_physics_init (Atm, atmos_axes, Time, physics_window, Surf_diff)

   call nullify_domain ( )

!--- initialize nudging module ---
#if defined (ATMOS_NUDGE)
    call atmos_nudge_init ( Time, atmos_axes(1:3), flag=do_atmos_nudge )
    if ( do_atmos_nudge .and. Atm(1)%nudge ) then
         call mpp_error(NOTE, 'Code compiled with atmospheric nudging, but fv_core_nml nudge is also set to .true.')
    elseif ( do_atmos_nudge) then
         call mpp_error(NOTE, 'Code compiled with and using atmospheric nudging')
    endif
    Atm(1)%nudge = do_atmos_nudge
#elif defined (CLIMATE_NUDGE)
    call fv_climate_nudge_init ( Time, atmos_axes(1:3), flag=do_atmos_nudge )
    if ( do_atmos_nudge .and. Atm(1)%nudge ) then
         call mpp_error(NOTE, 'Code compiled with climate nudging, but fv_core_nml nudge is also set to .true.')
    elseif ( do_atmos_nudge ) then
         call mpp_error(NOTE, 'Code compiled with and using climate nudging')
    endif
    Atm(1)%nudge = do_atmos_nudge
#else
   if ( Atm(1)%nudge ) then
        call fv_nwp_nudge_init( Time, atmos_axes, npz, zvir, Atm(1)%ak, Atm(1)%bk, Atm(1)%ts, Atm(1)%phis)
        call mpp_error(NOTE, 'NWP nudging is active')
   endif
#endif

! This call needs to be separate from the register nudging restarts after initialization
   call fv_io_register_nudge_restart ( Atm )

!miz
   if( Atm(1)%ncep_ic ) Surf_diff%sst_miz(:,:) = Atm(1)%ts(isc:iec, jsc:jec)

   id_tdt_dyn =register_diag_field(mod_name,'tdt_dyn',  atmos_axes(1:3),Time,'tdt_dyn', 'K/s', missing_value=mv)
   id_qdt_dyn =register_diag_field(mod_name,'qdt_dyn',  atmos_axes(1:3),Time,'qdt_dyn', 'kg/kg/s', missing_value=mv)
   id_qldt_dyn=register_diag_field(mod_name,'qldt_dyn', atmos_axes(1:3),Time,'qldt_dyn','kg/kg/s', missing_value=mv)
   id_qidt_dyn=register_diag_field(mod_name,'qidt_dyn', atmos_axes(1:3),Time,'qidt_dyn','kg/kg/s', missing_value=mv)
   id_qadt_dyn=register_diag_field(mod_name,'qadt_dyn', atmos_axes(1:3),Time,'qadt_dyn','1/s', missing_value=mv)

!yyf---allocate id_tracer_dyn 
   allocate (id_tracerdt_dyn    (num_tracers))
!yyf---loop for tracers
   do itrac = 1, num_tracers
     call get_tracer_names (MODEL_ATMOS, itrac, name = tracer_name, &
                                                  units = tracer_units)
     if (get_tracer_index(MODEL_ATMOS,tracer_name)>0) &
         id_tracerdt_dyn(itrac) = register_diag_field  &
             (mod_name, TRIM(tracer_name)//'dt_dyn', atmos_axes(1:3), &
             Time, TRIM(tracer_name)//' total tendency from advection',&
             TRIM(tracer_units)//'/s', missing_value = mv)
   enddo
   if (any(id_tracerdt_dyn(:)>0))   &
                     allocate(qtendyyf(isc:iec, jsc:jec,1:npz,num_tracers))
!yyf---end loop

   if ( id_tdt_dyn>0 ) allocate(ttend (isc:iec, jsc:jec, 1:npz))
   if ( id_qdt_dyn>0 .or. id_qldt_dyn>0 .or. id_qidt_dyn>0 .or. id_qadt_dyn>0 )   &
   allocate(qtend (isc:iec, jsc:jec, 1:npz, 4))
!miz

!  --- initialize clocks for dynamics, physics_down and physics_up
   id_dynam     = mpp_clock_id ('FV dynamical core',   &
          flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
   id_phys_down = mpp_clock_id ('Physics_down',   &
          flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
   id_phys_up   = mpp_clock_id ('Physics_up',   &
          flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
   id_fv_diag   = mpp_clock_id ('FV Diag',   &
          flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )

 end subroutine atmosphere_init



 subroutine atmosphere_down ( Time,    frac_land,          &
                           t_surf,  albedo,                &
                           albedo_vis_dir, albedo_nir_dir, &
                           albedo_vis_dif, albedo_nir_dif, &
                           rough_mom,                      &
                           u_star,  b_star, q_star,        &
                           dtau_du, dtau_dv, tau_x, tau_y, &
                           frac_open_sea,                  &
                           gust, coszen, flux_sw,          &
                           flux_sw_dir, flux_sw_dif,       &
                           flux_sw_down_vis_dir,           &
                           flux_sw_down_vis_dif,           &
                           flux_sw_down_total_dir,         &
                           flux_sw_down_total_dif,         &
                           flux_sw_vis,                    &
                           flux_sw_vis_dir,                &
                           flux_sw_vis_dif,                &
                           flux_lw,                        &
                           Surf_diff                       )
!
!        Time = time at the current time level
!
   type(time_type),intent(in)      :: Time
   real, intent(in), dimension(:,:):: frac_land, t_surf, albedo,      &
                                      albedo_vis_dir, albedo_nir_dir, &
                                      albedo_vis_dif, albedo_nir_dif, &
                                      rough_mom, u_star, b_star,      &
                                      q_star, dtau_du, dtau_dv,       &
                                      frac_open_sea
   real, intent(inout), dimension(:,:):: tau_x,  tau_y
   real, intent(out),   dimension(:,:):: gust, coszen, flux_sw,    &
                                         flux_sw_dir, flux_sw_dif, &
                                         flux_sw_down_vis_dir,     &
                                         flux_sw_down_total_dir,   &
                                         flux_sw_down_vis_dif,     &
                                         flux_sw_down_total_dif,   &
                                         flux_sw_vis,              &
                                         flux_sw_vis_dir,          &
                                         flux_sw_vis_dif, flux_lw
   type(surf_diff_type), intent(inout):: Surf_diff
   type(time_type) :: Time_prev, Time_next
   integer         :: itrac

   Time_prev = Time                       ! two time-level scheme
   Time_next = Time + Time_step_atmos

!---- Call FV dynamics -----

   call mpp_clock_begin (id_dynam)
                    call timing_on('fv_dynamics')

!miz
   if ( id_tdt_dyn>0 ) ttend(:, :, :) = Atm(1)%pt(isc:iec, jsc:jec, :)
   if ( id_qdt_dyn>0 .or. id_qldt_dyn>0 .or. id_qidt_dyn>0 .or. id_qadt_dyn>0 )   &
   qtend(:, :, :, :) = Atm(1)%q (isc:iec, jsc:jec, :, :)
!miz
   do itrac = 1, num_tracers
     if (id_tracerdt_dyn (itrac) >0 )   &
            qtendyyf(:,:,:,itrac) = Atm(1)%q(isc:iec, jsc:jec, :,itrac)
   enddo

   call fv_dynamics(npx, npy, npz, nq, Atm(1)%ng, dt_atmos, Atm(1)%consv_te,         &
                    Atm(1)%fill,  Atm(1)%reproduce_sum, kappa, cp_air, zvir,         &
                    Atm(1)%ks,    ncnst,        Atm(1)%n_split,   Atm(1)%q_split,    &
                    Atm(1)%u, Atm(1)%v, Atm(1)%w, Atm(1)%delz, Atm(1)%hydrostatic,   & 
                    Atm(1)%pt, Atm(1)%delp, Atm(1)%q, Atm(1)%ps,                     &
                    Atm(1)%pe, Atm(1)%pk, Atm(1)%peln, Atm(1)%pkz, Atm(1)%phis,      &
                    Atm(1)%omga, Atm(1)%ua, Atm(1)%va, Atm(1)%uc, Atm(1)%vc,         &
                    Atm(1)%ak, Atm(1)%bk, Atm(1)%mfx, Atm(1)%mfy,                    &
                    Atm(1)%cx, Atm(1)%cy, Atm(1)%ze0, Atm(1)%hybrid_z)

                    call timing_off('fv_dynamics')
!miz
   if ( id_tdt_dyn>0 ) then
        ttend = (Atm(1)%pt(isc:iec, jsc:jec, :)   - ttend(:, :, :   ))/dt_atmos
         used = send_data(id_tdt_dyn,  ttend(:,:,:),   Time)
   endif

   if ( id_qdt_dyn>0 .or. id_qldt_dyn>0 .or. id_qidt_dyn>0 .or. id_qadt_dyn>0 ) then
        qtend = (Atm(1)%q (isc:iec, jsc:jec, :, :)- qtend(:, :, :, :))/dt_atmos
        used = send_data(id_qdt_dyn,  qtend(:,:,:,1), Time)
        used = send_data(id_qldt_dyn, qtend(:,:,:,2), Time)
        used = send_data(id_qidt_dyn, qtend(:,:,:,3), Time)
        used = send_data(id_qadt_dyn, qtend(:,:,:,4), Time)
   endif
!miz

   do itrac = 1, num_tracers
     if(id_tracerdt_dyn(itrac)>0) then
       qtendyyf(:,:,:,itrac) = (Atm(1)%q (isc:iec, jsc:jec, :,itrac)-  &
                                        qtendyyf(:,:,:,itrac))/dt_atmos
       used = send_data(id_tracerdt_dyn(itrac), qtendyyf(:,:,:,itrac), &
                                                           Time)
     endif
   enddo

   call mpp_clock_end (id_dynam)


   call set_domain ( domain )
   call mpp_clock_begin (id_phys_down)
                         call timing_on('fv_physics_down')
   call fv_physics_down (Atm, dt_atmos, Time_prev, Time, Time_next,     &
                         frac_land, albedo,              &
                         albedo_vis_dir, albedo_nir_dir, &
                         albedo_vis_dif, albedo_nir_dif, &
                         rough_mom,  t_surf,             &
                         u_star,  b_star, q_star,        &
                         dtau_du, dtau_dv, tau_x, tau_y, &
                         flux_sw, flux_sw_dir,           &
                         flux_sw_dif,                    &
                         flux_sw_down_vis_dir,           &
                         flux_sw_down_vis_dif,           &
                         flux_sw_down_total_dir,         &
                         flux_sw_down_total_dif,         &
                         flux_sw_vis, flux_sw_vis_dir,   &
                         flux_sw_vis_dif, flux_lw,       &
                         coszen, gust, Surf_diff,        &
                         frac_open_sea )
                         call timing_off('fv_physics_down')
   call mpp_clock_end (id_phys_down)
   call nullify_domain ( )

 end subroutine atmosphere_down



 subroutine atmosphere_up ( Time,  frac_land, Surf_diff, lprec, fprec, gust, &
                            u_star, b_star, q_star )

   type(time_type),intent(in)         :: Time
   type(surf_diff_type), intent(inout):: Surf_diff
   real, intent(in),  dimension(:,:)  :: frac_land
   real, intent(inout), dimension(:,:):: gust
   real, intent(out), dimension(:,:)  :: lprec,   fprec
   real, intent(in), dimension(:,:)   :: u_star, b_star, q_star

   type(time_type) :: Time_prev, Time_next

   Time_prev = Time                       ! two time-level scheme
   Time_next = Time + Time_step_atmos

   call set_domain ( domain )
   call mpp_clock_begin (id_phys_up)
!-----------------------------------------------------------------------
                       call timing_on('fv_physics_up')
   call fv_physics_up( Atm, dt_atmos, Time_prev, Time, Time_next,      &
                       frac_land, Surf_diff, lprec, fprec, gust ,      &
                       u_star, b_star, q_star   )
                       call timing_off('fv_physics_up')
!-----------------------------------------------------------------------
   call mpp_clock_end (id_phys_up)

   call mpp_clock_begin(id_fv_diag)

   fv_time = Time_next
   call get_time (fv_time, seconds,  days)
!-----------------------------------------------------------------------
                call timing_on('FV_DIAG')
   call fv_diag( Atm, zvir, fv_time, Atm(1)%print_freq )
                call timing_off('FV_DIAG')     
!-----------------------------------------------------------------------
   call mpp_clock_end(id_fv_diag)
   call nullify_domain ( )


 end subroutine atmosphere_up



 subroutine atmosphere_end (Time, Grid_box)
   type (time_type),       intent(in) :: Time
   type(grid_box_type), intent(inout) :: Grid_box

  ! initialize domains for writing global physics data
   call set_domain ( domain )

   call get_time (Time, seconds,  days)

   call fv_physics_end(Atm, Time)

!--- end nudging module ---
#if defined (ATMOS_NUDGE)
   if ( Atm(1)%nudge ) call atmos_nudge_end
#elif defined (CLIMATE_NUDGE)
   if ( Atm(1)%nudge ) call fv_climate_nudge_end
#else
   if ( Atm(1)%nudge ) call fv_nwp_nudge_end
#endif

   call nullify_domain ( )
   call fv_end(Atm)
   deallocate (Atm)

   deallocate(Grid_box%dx)
   deallocate(Grid_box%dy)
   deallocate(Grid_box%area)
   deallocate(Grid_box%edge_w)
   deallocate(Grid_box%edge_e)
   deallocate(Grid_box%edge_s)
   deallocate(Grid_box%edge_n)
   deallocate(Grid_box%en1)
   deallocate(Grid_box%en2)
   deallocate(Grid_box%vlon)
   deallocate(Grid_box%vlat)

   if ( id_tdt_dyn>0 ) deallocate(ttend)
   if ( id_qdt_dyn>0 .or. id_qldt_dyn>0 .or. id_qidt_dyn>0 .or. id_qadt_dyn>0 ) deallocate(qtend)

   if (allocated(qtendyyf)) deallocate (qtendyyf)

 end subroutine atmosphere_end



  !#######################################################################
  ! <SUBROUTINE NAME="atmosphere_restart">
  ! <DESCRIPTION>
  !  Write out restart files registered through register_restart_file
  ! </DESCRIPTION>
  subroutine atmosphere_restart(timestamp)
    character(len=*),  intent(in) :: timestamp

    call fv_physics_restart(timestamp)
    call fv_write_restart(Atm, timestamp)

  end subroutine atmosphere_restart
  ! </SUBROUTINE>



 subroutine atmosphere_resolution (i_size, j_size, global)
   integer, intent(out)          :: i_size, j_size
   logical, intent(in), optional :: global
   logical :: local

   local = .true.
   if( PRESENT(global) ) local = .NOT.global

   if( local ) then
       i_size = iec - isc + 1
       j_size = jec - jsc + 1
   else
       i_size = npx - 1
       j_size = npy - 1
   end if

 end subroutine atmosphere_resolution



 subroutine atmosphere_cell_area  (area_out)
   real, dimension(:,:),  intent(out)          :: area_out       

   area_out(1:iec-isc+1, 1:jec-jsc+1) =  area (isc:iec,jsc:jec)                        

 end subroutine atmosphere_cell_area 





 subroutine atmosphere_boundary (blon, blat, global)
!---------------------------------------------------------------
!    returns the longitude and latitude grid box edges
!    for either the local PEs grid (default) or the global grid
!---------------------------------------------------------------
    real,    intent(out) :: blon(:,:), blat(:,:)   ! Unit: radian
    logical, intent(in), optional :: global
! Local data:
    integer i,j

    if( PRESENT(global) ) then
      if (global) call mpp_error(FATAL, '==> global grid is no longer available &
                               & in the Cubed Sphere')
    endif

    do j=jsc,jec+1
       do i=isc,iec+1
          blon(i-isc+1,j-jsc+1) = Atm(1)%grid(i,j,1)
          blat(i-isc+1,j-jsc+1) = Atm(1)%grid(i,j,2)
       enddo
    end do

 end subroutine atmosphere_boundary



 subroutine atmosphere_domain ( fv_domain )
   type(domain2d), intent(out) :: fv_domain
!  returns the domain2d variable associated with the coupling grid
!  note: coupling is done using the mass/temperature grid with no halos

   fv_domain = domain_for_coupler

 end subroutine atmosphere_domain



 subroutine get_atmosphere_axes ( axes )
   integer, intent(out) :: axes (:)

!----- returns the axis indices for the atmospheric (mass) grid -----
   if ( size(axes(:)) < 0 .or. size(axes(:)) > 4 ) call error_mesg (    &
                               'get_atmosphere_axes in atmosphere_mod', &
                               'size of argument is incorrect', FATAL   )

   axes (1:size(axes(:))) = atmos_axes (1:size(axes(:)))
 
 end subroutine get_atmosphere_axes




 subroutine get_bottom_mass ( t_bot, tr_bot, p_bot, z_bot, p_surf, slp )
!--------------------------------------------------------------
! returns temp, sphum, pres, height at the lowest model level
! and surface pressure
!--------------------------------------------------------------
   real, intent(out), dimension(isc:iec,jsc:jec):: t_bot, p_bot, z_bot, p_surf
   real, intent(out), optional, dimension(isc:iec,jsc:jec):: slp
   real, intent(out), dimension(isc:iec,jsc:jec,nq):: tr_bot
   integer :: i, j, m, k, kr
   real    :: rrg, sigtop, sigbot
   real, dimension(isc:iec,jsc:jec) :: tref
   real, parameter :: tlaps = 6.5e-3

   rrg  = rdgas / grav

   do j=jsc,jec
      do i=isc,iec
         p_surf(i,j) = Atm(1)%ps(i,j)
         t_bot(i,j) = Atm(1)%pt(i,j,npz)
         p_bot(i,j) = Atm(1)%delp(i,j,npz)/(Atm(1)%peln(i,npz+1,j)-Atm(1)%peln(i,npz,j))
         z_bot(i,j) = rrg*t_bot(i,j)*(1.+zvir*Atm(1)%q(i,j,npz,1)) *  &
                      (1. - Atm(1)%pe(i,npz,j)/p_bot(i,j))
      enddo
   enddo

   if ( present(slp) ) then
     ! determine 0.8 sigma reference level
     sigtop = Atm(1)%ak(1)/pstd_mks+Atm(1)%bk(1)
     do k = 1, npz 
        sigbot = Atm(1)%ak(k+1)/pstd_mks+Atm(1)%bk(k+1)
        if (sigbot+sigtop > 1.6) then
           kr = k  
           exit    
        endif   
        sigtop = sigbot
     enddo
     do j=jsc,jec
        do i=isc,iec
           ! sea level pressure
           tref(i,j) = Atm(1)%pt(i,j,kr) * (Atm(1)%delp(i,j,kr)/ &
                            ((Atm(1)%peln(i,kr+1,j)-Atm(1)%peln(i,kr,j))*Atm(1)%ps(i,j)))**(-rrg*tlaps)
           slp(i,j) = Atm(1)%ps(i,j)*(1.+tlaps*Atm(1)%phis(i,j)/(tref(i,j)*grav))**(1./(rrg*tlaps))
        enddo
     enddo
   endif

! Copy tracers
   do m=1,nq
      do j=jsc,jec
         do i=isc,iec
            tr_bot(i,j,m) = Atm(1)%q(i,j,npz,m)
         enddo
      enddo
   enddo

 end subroutine get_bottom_mass



 subroutine get_bottom_wind ( u_bot, v_bot )
!-----------------------------------------------------------
! returns u and v on the mass grid at the lowest model level
!-----------------------------------------------------------
   real, intent(out), dimension(isc:iec,jsc:jec):: u_bot, v_bot
   integer i, j

   do j=jsc,jec
      do i=isc,iec
         u_bot(i,j) = Atm(1)%u_srf(i,j)
         v_bot(i,j) = Atm(1)%v_srf(i,j)
      enddo
   enddo

 end subroutine get_bottom_wind



 subroutine get_stock_pe(index, value)
   integer, intent(in) :: index
   real,   intent(out) :: value

#ifdef USE_STOCK
   include 'stock.inc' 
#endif

   real wm(isc:iec,jsc:jec)
   integer i,j,k
   
   select case (index)

#ifdef USE_STOCK
   case (ISTOCK_WATER)
#else
   case (1)
#endif
     
!----------------------
! Perform vertical sum:
!----------------------
     wm = 0.
     do j=jsc,jec
        do k=1,npz
           do i=isc,iec
! Warning: the following works only with AM2 physics: water vapor; cloud water, cloud ice.
              wm(i,j) = wm(i,j) + Atm(1)%delp(i,j,k) * ( Atm(1)%q(i,j,k,1) +    &
                                                         Atm(1)%q(i,j,k,2) +    &
                                                         Atm(1)%q(i,j,k,3) )
           enddo
        enddo
     enddo

!----------------------
! Horizontal sum:
!----------------------
     value = 0.
     do j=jsc,jec
        do i=isc,iec
           value = value + wm(i,j)*area(i,j)
        enddo
     enddo
     value = value/grav

   case default
     value = 0.0
   end select

 end subroutine get_stock_pe 

end module atmosphere_mod
