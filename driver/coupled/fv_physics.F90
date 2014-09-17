module fv_physics_mod

!-----------------------------------------------------------------------
!
!  Interface for Cubed_sphere FV dynamics with GFDL atmospheric physics
!  History: modified by SJL based on Memphis release
!           tq_filter option removed, 2012/08/06
!-----------------------------------------------------------------------

!-----------------
! FMS modules:
!-----------------
use atmos_co2_mod,         only: atmos_co2_rad, co2_radiation_override
use constants_mod,         only: rdgas, grav, rvgas, WTMAIR, WTMCO2, cp_air
use time_manager_mod,      only: time_type, get_time, operator(-)
use fms_mod,               only: error_mesg, FATAL, NOTE, WARNING, write_version_number,clock_flag_default
use physics_driver_mod,    only: physics_driver_init, physics_driver_end,   &
                                 physics_driver_moist_init, &
                                 physics_driver_moist_end, &
                                 physics_driver_down, physics_driver_up, surf_diff_type, &
                                 physics_driver_down_time_vary, &
                                 physics_driver_up_time_vary,  &
                                 physics_driver_down_endts,  &
                                 physics_driver_up_endts, &
                                 physics_driver_restart
use field_manager_mod,     only: MODEL_ATMOS
use tracer_manager_mod,    only: get_tracer_index, NO_TRACER
use mpp_domains_mod,       only: mpp_global_sum, BITWISE_EXACT_SUM, NON_BITWISE_EXACT_SUM
use mpp_mod,               only: mpp_error, mpp_clock_id,  mpp_clock_begin,  &
                                 mpp_clock_end, CLOCK_MODULE_DRIVER, mpp_pe

!-----------------
! FV core modules:
!-----------------
use fv_grid_utils_mod,     only: g_sum
use fv_arrays_mod,         only: fv_atmos_type
use fv_eta_mod,            only: get_eta_level
use fv_update_phys_mod,    only: fv_update_phys, del2_phys
use fv_sg_mod,             only: fv_dry_conv
use fv_timing_mod,         only: timing_on, timing_off
use mpp_domains_mod,       only: domain2d !domain data type

#ifdef _OPENMP
use omp_lib
#endif

implicit none
private

public  fv_physics_down, fv_physics_up, fv_physics_init, fv_physics_end
public  surf_diff_type, fv_physics_restart

!-----------------------------------------------------------------------

   real, allocatable, dimension(:,:,:)   :: u_dt, v_dt, t_dt
   real, allocatable, dimension(:,:,:,:) :: q_dt  ! potentially a huge array
   real, allocatable, dimension(:,:,:)   :: p_full, z_full, p_half, z_half
   real    :: zvir, rrg, ginv
   integer :: id_fv_physics_down, id_fv_physics_up, id_fv_update_phys
   integer :: isc, iec, jsc, jec, ngc, nt_prog, nt_dyn
   integer :: isd, ied, jsd, jed
   integer :: isw, iew, jsw, jew  ! window start/end in global index space
   integer :: nx_win, ny_win      ! iew-isw+1, jew-jsw+1 (window sizes)
   integer :: nx_dom, ny_dom      ! ie-is+1, je-js+1 (compute domain sizes)
   integer :: sphum, w_diff
   integer, allocatable, dimension(:)  :: physics_window_x, physics_window_y
   integer :: ny_per_thread, num_phys_windows
   real, parameter:: w0_big = 60.  ! to prevent negative w-tracer diffusion

   integer :: numthreads
!   integer, EXTERNAL :: omp_get_thread_num, omp_get_num_threads      

!---- version number -----
   character(len=128) :: version = '$Id: fv_physics.F90,v 20.0.4.1.2.1 2014/09/04 19:06:26 Rusty.Benson Exp $'
   character(len=128) :: tagname = '$Name: tikal_201409 $'

contains


  subroutine fv_physics_init (Atm, axes, Time, window, Surf_diff)
    type(fv_atmos_type), intent(inout) :: Atm(:)
!-----------------------------------------------------------------------
!   axes      = array of axis indices for diagnostics (x,y,pf,ph)
!   Time      = current time (time_type)
!-----------------------------------------------------------------------
    integer,               intent(in)    :: axes(4)
    integer,               intent(in)    :: window(2)
    type (time_type),      intent(in)    :: Time
    type (surf_diff_type), intent(inout) :: Surf_diff
!-----------------------------------------------------------------------
! Local:
    character(len=132)  :: text
    real, allocatable   :: p_edge(:,:,:)  ! used by atmos_tracer_driver_init

    real    ::  pref(Atm(1)%npz+1,2)
    real    :: phalf(Atm(1)%npz+1)
    real    :: ps1, ps2
    integer :: i, j, k
    integer :: ios
    character(len=80) evalue

! All tracers are prognostic
    nt_prog = Atm(1)%ncnst - Atm(1)%flagstruct%pnats
   
    isc = Atm(1)%bd%isc
    iec = Atm(1)%bd%iec
    jsc = Atm(1)%bd%jsc
    jec = Atm(1)%bd%jec

    ngc = Atm(1)%ng

    isd = isc - ngc
    ied = iec + ngc
    jsd = jsc - ngc
    jed = jec + ngc

    zvir = rvgas/rdgas - 1.
    ginv = 1./ grav
    rrg  = rdgas / grav        

    numthreads = 1
!$OMP PARALLEL
!$OMP MASTER
!$  numthreads = omp_get_num_threads()
!$OMP END MASTER
!$OMP END PARALLEL

!----- write version to logfile --------
    call write_version_number (version,tagname)

! Specific humidity is assumed to be q(:,:,:,1)
    sphum = get_tracer_index (MODEL_ATMOS, 'sphum' )
    if(sphum /= 1) call error_mesg('fv_physics_init:','sphum /= 1', FATAL)

    w_diff = get_tracer_index (MODEL_ATMOS, 'w_diff' )
    if ( w_diff /= NO_TRACER ) then
         nt_dyn = nt_prog - 1
    else
         nt_dyn = nt_prog
    endif

!---------- reference profile -----------
    ps1 = 101325.
    ps2 =  81060.
    pref(Atm(1)%npz+1,1) = ps1
    pref(Atm(1)%npz+1,2) = ps2
    call get_eta_level ( Atm(1)%npz, ps1, pref(1,1), phalf, Atm(1)%ak, Atm(1)%bk )
    call get_eta_level ( Atm(1)%npz, ps2, pref(1,2), phalf, Atm(1)%ak, Atm(1)%bk )

    allocate (p_edge(isc:iec,jsc:jec, Atm(1)%npz+1))
!------- pressure at model layer interfaces -----
    do k=1,Atm(1)%npz+1
       do j=jsc,jec
          do i=isc,iec
             p_edge(i,j,k) = Atm(1)%pe(i,k,j)
          enddo
       enddo
    enddo
!---------- initialize physics -------

    call physics_driver_init(Time, Atm(1)%gridstruct%grid(isc:iec+1,jsc:jec+1,1),             &
                                   Atm(1)%gridstruct%grid(isc:iec+1,jsc:jec+1,2),             &
                                   Atm(1)%gridstruct%agrid(isc:iec,jsc:jec,1),                &
                                   Atm(1)%gridstruct%agrid(isc:iec,jsc:jec,2),                &
                             axes, pref, Atm(1)%q(isc:iec,jsc:jec,1:Atm(1)%npz,1:Atm(1)%ncnst),  &
                             Surf_diff,  p_edge )
    deallocate ( p_edge )

    allocate (  u_dt(isd:ied,jsd:jed, Atm(1)%npz) )
    allocate (  v_dt(isd:ied,jsd:jed, Atm(1)%npz) )
    allocate (  t_dt(isc:iec,jsc:jec, Atm(1)%npz) )
    allocate (  q_dt(isc:iec,jsc:jec, Atm(1)%npz, nt_prog) )

! physics window
    nx_dom = iec - isc + 1
    ny_dom = jec - jsc + 1

    nx_win = window(1)
    ny_win = window(2)

    if( nx_win.LE.0) nx_win = nx_dom
    if( ny_win.LE.0) ny_win = ny_dom

!!$    if( nx_win.LE.0 .or. Atm(1)%neststruct%nested) nx_win = nx_dom
!!$    if( ny_win.LE.0 .or. Atm(1)%neststruct%nested) ny_win = ny_dom

! Consistency check:
    if( mod(nx_dom,nx_win).NE.0 )then
        write( text,'(a,i5,a,i5)' )'FV_PHYSICS_INIT: atmosphere_nml problem,'// &
             ' physics_window must divide domain size for X-dimension - window:'&
             ,nx_win, ' domain:', nx_dom
        call mpp_error( FATAL, text )
    end if
    if( mod(ny_dom,ny_win).NE.0 )then
        write( text,'(a,i5,a,i5)' )'FV_PHYSICS_INIT: atmosphere_nml problem,'// &
             ' physics_window must divide domain size for Y-dimension - window:'&
             , ny_win, ' domain: ',ny_dom
        call mpp_error( FATAL, text )
    end if

    allocate( p_full(isc:iec,jsc:jec,Atm(1)%npz) )
    allocate( z_full(isc:iec,jsc:jec,Atm(1)%npz) )
    allocate( p_half(isc:iec,jsc:jec,Atm(1)%npz+1) )
    allocate( z_half(isc:iec,jsc:jec,Atm(1)%npz+1) )

!MPP clocks
    id_fv_physics_down = mpp_clock_id( 'FV_PHYSICS_DOWN', &
         flags=clock_flag_default, grain=CLOCK_MODULE_DRIVER )
    id_fv_physics_up = mpp_clock_id( 'FV_PHYSICS_UP', &
         flags=clock_flag_default, grain=CLOCK_MODULE_DRIVER )
    id_fv_update_phys = mpp_clock_id( 'FV_UPDATE_PHYS', &
         flags=clock_flag_default, grain=CLOCK_MODULE_DRIVER )

    ny_per_thread = max(1,ny_win/numthreads)
!    if (mod(ny_win, numthreads ) /= 0) then
!      call error_mesg ('physics_driver_down', &
!         'The number of OpenMP threads must be an integral multiple &
!                  &of the number of rows in the physics window', FATAL)
!    endif

    num_phys_windows = (nx_dom/nx_win)*(ny_dom/ny_win) 
    write(text,'(a,2i4)') 'num_phys_windows, numthreads',num_phys_windows,numthreads
    call error_mesg ('fv_physics_init', trim(text), WARNING) ! NOTE)
    allocate(physics_window_x(num_phys_windows))
    allocate(physics_window_y(num_phys_windows))
    i = 1 
    do jsw = jsc,jec,ny_win
       do isw = isc,iec,nx_win
          physics_window_x(i) =isw
          physics_window_y(i) =jsw
          i = i + 1
       enddo
    enddo
    if (numthreads > num_phys_windows ) then
      call error_mesg ('fv_physics_init', &
         'The number of OpenMP threads is greater than the number of physics windows. &
                  &Please use more physics windows via atmosphere_nml.', FATAL)
    endif
    
  end subroutine fv_physics_init



  subroutine fv_physics_down(Atm, dt_phys, Time_prev, Time, Time_next, &
                             frac_land,   albedo,            &
                             albedo_vis_dir, albedo_nir_dir, &
                             albedo_vis_dif, albedo_nir_dif, &
                             rough_vel,   t_surf,            &
                             t_ref, q_ref,                   &  ! cjg: PBL depth mods
                             u_star, b_star, q_star,         &
                             dtau_du, dtau_dv, tau_x, tau_y, &
                             flux_sw,                        &
                             flux_sw_dir, flux_sw_dif,       &
                             flux_sw_down_vis_dir,           &
                             flux_sw_down_vis_dif,           &
                             flux_sw_down_total_dir,         &
                             flux_sw_down_total_dif,         &
                             flux_sw_vis, flux_sw_vis_dir,   &
                             flux_sw_vis_dif,                &
                             flux_lw, coszen,                &
                             gust, Surf_diff, frac_open_sea )
!-----------------------------------------------------------------------
!
!   Time_prev =  time at the previous time level, tau-1 (time_type)
!   Time      =  time at the current time level,  tau   (time_type)
!   Time_next =  time at the next time level,     tau+1 (time_type)
!
!   NOTE: for a two time level scheme (e.g., forward-backward scheme)
!         Time_prev = Time.
!
!-----------------------------------------------------------------------
    type(time_type),     intent(in) :: Time_prev, Time, Time_next
    type(fv_atmos_type), intent(inout) :: Atm
    real,                intent(in) :: dt_phys
    real, intent(in), dimension(isc:iec,jsc:jec):: frac_land,  albedo, &
                                       albedo_vis_dir, albedo_nir_dir, &
                                       albedo_vis_dif, albedo_nir_dif, &
                                       rough_vel, t_surf,              &
                                       t_ref, q_ref,                   &  ! cjg: PBL depth mods
                                       u_star, b_star, q_star,         &
                                       dtau_du, dtau_dv,               &
                                       frac_open_sea

    type(surf_diff_type), intent(inout) :: Surf_diff
    real, intent(inout), dimension(isc:iec,jsc:jec):: tau_x, tau_y
    real, intent(out),   dimension(isc:iec,jsc:jec):: flux_sw, &
                                   flux_sw_dir, flux_sw_dif,   &
                                   flux_sw_down_vis_dir,       &
                                   flux_sw_down_vis_dif,       &
                                   flux_sw_down_total_dir,     &
                                   flux_sw_down_total_dif,     &
                                   flux_sw_vis,                &
                                   flux_sw_vis_dir,            &
                                   flux_sw_vis_dif, flux_lw, coszen, gust
!-----------------------------------------------------------------------
    real :: gavg_rrv(nt_prog)
    real :: gavg_rrvp(nt_prog+1)
    real :: atm_mass
    real, dimension(isc:iec,jsc:jec,1):: sum1
    integer:: iq, idx, phys_loop
    integer:: i, j, k
    real    :: dt 
    integer :: sec, day

!----------------------------------------------------------------------
! obtain pressure-weighted global mean co2 dry volume mixing ratio for
! use by radiation package.
!----------------------------------------------------------------------
    gavg_rrv = 0.
! check to override predicted global pressure-weighted rad co2
    idx = get_tracer_index(MODEL_ATMOS, 'co2')
    if(idx /= NO_TRACER .and. co2_radiation_override) then
      call atmos_co2_rad(Time, gavg_rrv(idx))
    elseif (idx /= NO_TRACER) then
      call compute_g_avg(gavg_rrv, 'co2', Atm%pe, Atm%q, Atm%gridstruct%area, Atm%domain, Atm%npz, Atm%ncnst)
    endif

!------------------------------------------------------------------------
!  compute global mean atmospheric mass to be used later in diagnostic.
!------------------------------------------------------------------------
    sum1 = 0.
    do j=jsc,jec
      do i=isc,iec
        sum1(i,j,1) = Atm%pe(i,Atm%npz+1,j)*Atm%gridstruct%area(i,j)
      end do
    end do
    atm_mass = REAL(mpp_global_sum(Atm%domain, sum1,   &
                                    flags=NON_BITWISE_EXACT_SUM),KIND=4)

!---------------------------------------------------------------------
! compute the physics time step (from tau-1 to tau+1).
!---------------------------------------------------------------------
    call get_time (Time_next-Time_prev, sec, day)
    dt = real(sec+day*86400)
 
!---------------------------------------------------------------------
! call physics_driver_down_time_vary to do the time-dependent, spatially
! independent calculations before entering windows / threads loop. 
!--------------------------------------------------------------------- 
    call physics_driver_down_time_vary (Time, Time_next, gavg_rrv, dt)


    if ( Atm%flagstruct%fv_sg_adj > 0 ) then
         call fv_dry_conv( isd, ied, jsd, jed, isc, iec, jsc, jec, Atm%npz, nt_dyn, dt_phys,   &
                           Atm%flagstruct%fv_sg_adj, Atm%delp, Atm%pe, Atm%peln,  &
                           Atm%pkz, Atm%pt, Atm%q, Atm%ua, Atm%va,  &
                           Atm%flagstruct%hydrostatic, Atm%w, Atm%delz, u_dt, v_dt, t_dt, q_dt )
    else
         u_dt(:,:,:) = 0.
         v_dt(:,:,:) = 0.
         v_dt(:,:,:) = 0.
         t_dt(:,:,:) = 0.
         q_dt(:,:,:,:) = 0.
    endif

    if ( .not. Atm%flagstruct%hydrostatic .and. w_diff /= NO_TRACER ) then
!$OMP parallel do default(shared) private(i, j, k)
       do k=1, Atm%npz
          do j=jsc, jec
             do i=isc, iec
                Atm%q(i,j,k,w_diff) = Atm%w(i,j,k) + w0_big
                    q_dt(i,j,k,w_diff) = 0.
             enddo
          enddo
        enddo
    endif

    call mpp_clock_begin(id_fv_physics_down)


    call compute_p_z(Atm%npz, isc, jsc, nx_dom, ny_dom, Atm%phis, Atm%pt, &
                     Atm%q, Atm%delp, Atm%pe, Atm%peln,         &
                     Atm%delz, Atm%flagstruct%phys_hydrostatic)

    gavg_rrvp(1:nt_prog) = gavg_rrv(:)
    gavg_rrvp(nt_prog+1) = atm_mass   

!$OMP parallel do schedule(dynamic) default(shared) private(phys_loop, isw, iew, jsw, jew)
    do phys_loop = 1, size(physics_window_y)!num_phys_windows
       jsw = physics_window_y(phys_loop)
       jew = jsw + ny_win - 1
       isw = physics_window_x(phys_loop)
       iew = isw + nx_win - 1

          call physics_driver_down( isw-isc+1, iew-isc+1, jsw-jsc+1, jew-jsc+1, &
                   Time_prev, Time, Time_next                              , &
                   Atm%gridstruct%agrid(isw:iew,jsw:jew,2)                         , &
                   Atm%gridstruct%agrid(isw:iew,jsw:jew,1)                         , &
                   Atm%gridstruct%area(isw:iew,jsw:jew)                 , &
                   p_half(isw:iew,jsw:jew,:), p_full(isw:iew,jsw:jew,:)    , &
                   z_half(isw:iew,jsw:jew,:), z_full(isw:iew,jsw:jew,:)    , &
                   p_half(isw:isw,jsw:jsw,:)                               , &
                   Atm%ua(isw:iew,jsw:jew,:)                            , &
                   Atm%va(isw:iew,jsw:jew,:)                            , &
                   Atm%pt(isw:iew,jsw:jew,:)                            , &
                   Atm%q (isw:iew,jsw:jew,:,1)                          , &
                   Atm%q (isw:iew,jsw:jew,:,:)                          , &
                   Atm%ua(isw:iew,jsw:jew,:)                            , &
                   Atm%va(isw:iew,jsw:jew,:)                            , &
                   Atm%pt(isw:iew,jsw:jew,:)                            , &
                   Atm%q (isw:iew,jsw:jew,:,1)                          , &
                   Atm%q (isw:iew,jsw:jew,:,:)                          , &
                   frac_land(isw:iew,jsw:jew), rough_vel(isw:iew,jsw:jew)  , &
                   frac_open_sea(isw:iew,jsw:jew),                           &
                   albedo   (isw:iew,jsw:jew)                              , &
                   albedo_vis_dir(isw:iew,jsw:jew)                         , &
                   albedo_nir_dir(isw:iew,jsw:jew)                         , &
                   albedo_vis_dif(isw:iew,jsw:jew)                         , &
                   albedo_nir_dif(isw:iew,jsw:jew)                         , &
                   t_surf  (isw:iew,jsw:jew)                               , &
                   t_ref(isw:iew,jsw:jew), q_ref(isw:iew,jsw:jew)          , & ! cjg: PBL depth mods  
                   u_star(isw:iew,jsw:jew)                                 , &
                   b_star  (isw:iew,jsw:jew),  q_star(isw:iew,jsw:jew)     , &
                   dtau_du (isw:iew,jsw:jew), dtau_dv(isw:iew,jsw:jew)     , &
                   tau_x   (isw:iew,jsw:jew),   tau_y(isw:iew,jsw:jew)     , &
                   u_dt    (isw:iew,jsw:jew,:),  v_dt(isw:iew,jsw:jew,:)   , &
                   t_dt    (isw:iew,jsw:jew,:),  q_dt(isw:iew,jsw:jew,:,1) , &
                   q_dt    (isw:iew,jsw:jew,:,1:nt_prog)                   , &
                   flux_sw               (isw:iew,jsw:jew)                 , &
                   flux_sw_dir           (isw:iew,jsw:jew)                 , &
                   flux_sw_dif           (isw:iew,jsw:jew)                 , &
                   flux_sw_down_vis_dir  (isw:iew,jsw:jew)                 , &
                   flux_sw_down_vis_dif  (isw:iew,jsw:jew)                 , &
                   flux_sw_down_total_dir(isw:iew,jsw:jew)                 , &
                   flux_sw_down_total_dif(isw:iew,jsw:jew)                 , &
                   flux_sw_vis           (isw:iew,jsw:jew)                 , &
                   flux_sw_vis_dir       (isw:iew,jsw:jew)                 , &
                   flux_sw_vis_dif       (isw:iew,jsw:jew)                 , &
                   flux_lw               (isw:iew,jsw:jew)                 , &
                   coszen                (isw:iew,jsw:jew)                 , &
                   gust                  (isw:iew,jsw:jew)                 , &
                   Surf_diff,   gavg_rrvp )
    enddo


    call physics_driver_down_endts (1, 1)


    call mpp_clock_end(id_fv_physics_down)

  end subroutine fv_physics_down



  subroutine fv_physics_up( Atm, dt_phys, Time_prev, Time, Time_next, &
                            frac_land, Surf_diff, lprec, fprec, gust, &
                            u_star, b_star, q_star )
!-----------------------------------------------------------------------
!
!   Time_prev =  time at the previous time level, tau-1 (time_type)
!   Time      =  time at the current time level,  tau   (time_type)
!   Time_next =  time at the next time level,     tau+1 (time_type)
!
!   NOTE: for a two time level scheme (e.g., forward-backward scheme)
!         Time_prev = Time.
!
!-----------------------------------------------------------------------
    type(time_type),     intent(in)    :: Time_prev, Time, Time_next
    type(fv_atmos_type), intent(inout) :: Atm
    type(surf_diff_type),intent(inout) :: Surf_diff
    real,                intent(in)    :: dt_phys
    real, intent(in),    dimension(isc:iec,jsc:jec) :: frac_land
    real, intent(inout), dimension(isc:iec,jsc:jec) :: gust
    real, intent(out),   dimension(isc:iec,jsc:jec) :: lprec, fprec
    real, intent(in),    dimension(isc:iec,jsc:jec) :: u_star, b_star, q_star
    integer seconds, days
    integer :: i, j, k
    real::  rcp
    integer :: phys_loop
    integer :: sec, day
    real    :: dt
    type(time_type) :: Time_step

    rcp = 1. / cp_air

    call mpp_clock_begin(id_fv_physics_up)

!---------------------------------------------------------------------
!    compute the physics time step (from tau-1 to tau+1).
!---------------------------------------------------------------------
    call get_time (Time_next-Time_prev, sec, day)
    dt = real(sec+day*86400)
 
    call physics_driver_up_time_vary (Time, Time_next, dt)

    call compute_p_z(Atm%npz, isc , jsc , nx_dom, ny_dom, Atm%phis, Atm%pt,  &
                     Atm%q, Atm%delp, Atm%pe, Atm%peln,      &
                     Atm%delz, Atm%flagstruct%phys_hydrostatic)

    call physics_driver_moist_init (nx_dom, ny_dom,  Atm%npz, nt_prog, Atm%ncnst) 

!$OMP parallel  do default(shared) private(phys_loop, isw, iew, jsw, jew)
    do phys_loop = 1, size(physics_window_y)!num_phys_windows
       jsw = physics_window_y(phys_loop)
       jew = jsw + ny_win - 1  
       isw = physics_window_x(phys_loop)
       iew = isw + nx_win - 1

          call physics_driver_up( isw-isc+1, iew-isc+1, jsw-jsc+1, jew-jsc+1, &
                                  Time_prev, Time, Time_next             , &
                                  Atm%gridstruct%agrid(isw:iew,jsw:jew,2)        , &
                                  Atm%gridstruct%agrid(isw:iew,jsw:jew,1)        , &
                                  Atm%gridstruct%area(isw:iew,jsw:jew), &
                                  p_half(isw:iew,jsw:jew,:)              , &
                                  p_full(isw:iew,jsw:jew,:)              , &
                                  z_half(isw:iew,jsw:jew,:)              , &
                                  z_full(isw:iew,jsw:jew,:)              , &
                                  Atm%omga(isw:iew,jsw:jew,:)         , &
                                  Atm%ua(isw:iew,jsw:jew,:)           , &
                                  Atm%va(isw:iew,jsw:jew,:)           , &
                                  Atm%pt(isw:iew,jsw:jew,:)           , &
                                  Atm%q(isw:iew,jsw:jew,:,1)          , &
                                  Atm%q(isw:iew,jsw:jew,:,:)          , &   ! cjg: pass all tracers
                                  Atm%ua(isw:iew,jsw:jew,:)           , &
                                  Atm%va(isw:iew,jsw:jew,:)           , &
                                  Atm%pt(isw:iew,jsw:jew,:)           , &
                                  Atm%q(isw:iew,jsw:jew,:,1)          , &
                                  Atm%q(isw:iew,jsw:jew,:,:)          , &   ! cjg: pass all tracers
                                  frac_land(isw:iew,jsw:jew)             , &
                                  u_star   (isw:iew,jsw:jew)             , &
                                  b_star   (isw:iew,jsw:jew)             , &
                                  q_star   (isw:iew,jsw:jew)             , &
                                  u_dt     (isw:iew,jsw:jew,:)           , &
                                  v_dt     (isw:iew,jsw:jew,:)           , &
                                  t_dt     (isw:iew,jsw:jew,:)           , &
                                  q_dt     (isw:iew,jsw:jew,:,1)         , &
                                  q_dt     (isw:iew,jsw:jew,:,1:nt_prog) , &
                                  Surf_diff                              , &
                                  lprec    (isw:iew,jsw:jew)             , &
                                  fprec    (isw:iew,jsw:jew)             , &
                                  gust     (isw:iew,jsw:jew)             , &
                                  hydrostatic=Atm%flagstruct%hydrostatic , &
                                  phys_hydrostatic=Atm%flagstruct%phys_hydrostatic )
    enddo

    if ( .not.Atm%flagstruct%hydrostatic .and. w_diff /= NO_TRACER ) then
!$OMP parallel do default(shared) private(i, j, k)
       do k=1, Atm%npz
          do j=jsc, jec
             do i=isc, iec
                Atm%q(i,j,k,w_diff) = q_dt(i,j,k,w_diff) ! w tendency due to phys
                t_dt(i,j,k) = t_dt(i,j,k) - q_dt(i,j,k,w_diff)*rcp*(Atm%w(i,j,k)+0.5*dt_phys*q_dt(i,j,k,w_diff))
                Atm%w(i,j,k) = Atm%w(i,j,k) + dt_phys*Atm%q(i,j,k,w_diff)
! Heating due to loss of KE (vertical diffusion of w)
             enddo
          enddo
        enddo
    endif

    call physics_driver_moist_end

    call physics_driver_up_endts (1, 1)

    call mpp_clock_end(id_fv_physics_up)

    call mpp_clock_begin(id_fv_update_phys)
                                                            call timing_on('update_fv')
    call fv_update_phys( dt_phys,   isc,        iec,         jsc,    jec,   isd,                   &
                         ied,       jsd,        jed,         ngc,       nt_dyn,                    &
                         Atm%u,  Atm%v,   Atm%delp, Atm%pt, Atm%q,                                 &
                         Atm%ua, Atm%va,  Atm%ps,   Atm%pe, Atm%peln,                              &
                         Atm%pk, Atm%pkz, Atm%ak,   Atm%bk, Atm%phis,                              &
                         Atm%u_srf, Atm%v_srf, Atm%ts, Atm%delz, Atm%flagstruct%hydrostatic,       &
                         u_dt, v_dt, t_dt, q_dt, .true., Time_next, Atm%flagstruct%nudge,          &
                         Atm%gridstruct, Atm%gridstruct%agrid(:,:,1), Atm%gridstruct%agrid(:,:,2), &
                         Atm%npx, Atm%npy, Atm%npz, Atm%flagstruct, Atm%neststruct, Atm%bd,        &
                         Atm%domain, Atm%ptop)
                                                            call timing_off('update_fv')
    call mpp_clock_end(id_fv_update_phys)

  end subroutine fv_physics_up



  subroutine fv_physics_end (Atm, Time)
    type(fv_atmos_type), intent(inout) :: Atm(:)
    type(time_type), intent(in) :: Time
!                                 NOTE: this is not the dynamics time
    call physics_driver_end (Time)

    deallocate ( u_dt )
    deallocate ( v_dt )
    deallocate ( t_dt )
    deallocate ( q_dt )
    deallocate ( p_full )
    deallocate ( z_full )
    deallocate ( p_half )
    deallocate ( z_half )

  end subroutine fv_physics_end



  !#######################################################################
  ! <SUBROUTINE NAME="fv_physics_restart">
  ! <DESCRIPTION>
  !  Write out restart files registered through register_restart_file
  ! </DESCRIPTION>
  subroutine fv_physics_restart(timestamp)
    character(len=*),  intent(in) :: timestamp

    call physics_driver_restart(timestamp)

  end subroutine fv_physics_restart
  ! </SUBROUTINE>  



  subroutine compute_p_z (nlev, istart, jstart, isiz, jsiz, phis, pt, q,   &
                          delp, pe, peln, delz, hydrostatic)
    integer, intent(in):: nlev
    integer, intent(in):: istart, jstart, isiz, jsiz
    real,    intent(in):: phis(isd:ied,jsd:jed)
    real,    intent(in)::   pt(isd:ied,jsd:jed,nlev)
    real,    intent(in)::    q(isd:ied,jsd:jed,nlev,sphum)
    real,    intent(in):: delp(isd:ied,jsd:jed,nlev)
    real,    intent(in)::   pe(isc-1:iec+1,nlev+1,jsc-1:jec+1)
    real,    intent(in):: peln(isc  :iec,  nlev+1,jsc  :jec)
    real,    intent(in):: delz(isd:ied,jsd:jed,nlev)
    logical, intent(in):: hydrostatic
! local
    integer i,j,k,id,jd
    real    tvm

!----------------------------------------------------
! Compute pressure and height at full and half levels
!----------------------------------------------------

!$OMP parallel do default(shared) private(jd,id)
    do j=1,jsiz
       jd = j + jstart - 1
       do i=1,isiz
          id = i + istart - 1
          z_half(id,jd,nlev+1) = phis(id,jd) * ginv
       enddo
    end do

!$OMP parallel do default(shared) private(jd,id)
    do k=1,nlev+1
       do j=1,jsiz
          jd = j + jstart - 1
          do i=1,isiz
             id = i + istart - 1
             p_half(id,jd,k) = pe(id,k,jd)
          enddo
       enddo
    enddo

    if ( hydrostatic ) then
      do k=nlev,1,-1
!Rusty recommends explicitly making loop indices private
!Also here and in the next loop OMP needs to be inside the k loop
! since there is a k+1 dependency
!$OMP parallel do default(shared) private(jd,id, tvm, i, j)
         do j=1,jsiz
            jd = j + jstart - 1
            do i=1,isiz
               id = i + istart - 1
               tvm = rrg*pt(id,jd,k)*(1.+zvir*q(id,jd,k,sphum))
               p_full(id,jd,k) = delp(id,jd,k)/(peln(id,k+1,jd)-peln(id,k,jd))
               z_full(id,jd,k) = z_half(id,jd,k+1) + tvm*(1.-p_half(id,jd,k)/p_full(id,jd,k))
               z_half(id,jd,k) = z_half(id,jd,k+1) + tvm*(peln(id,k+1,jd)-peln(id,k,jd))
            enddo
         enddo
      enddo
    else
!--------- Non-Hydrostatic option ------------------------------------------
      do k=nlev,1,-1
!$OMP parallel do default(shared) private(jd,id,j,i)
         do j=1,jsiz
            jd = j + jstart - 1
            do i=1,isiz
               id = i + istart - 1
               p_full(id,jd,k) = delp(id,jd,k)/(peln(id,k+1,jd)-peln(id,k,jd))
               z_half(id,jd,k) = z_half(id,jd,k+1) - delz(id,jd,k)
               z_full(id,jd,k) = 0.5*(z_half(id,jd,k) + z_half(id,jd,k+1))
            enddo
         enddo
      enddo
!--------- Non-Hydrostatic option ------------------------------------------
    endif

  end subroutine compute_p_z



  subroutine compute_g_avg(rrv, tracer_name, pe, q, area, domain, npz, ncnst )
    real,          intent(inout) :: rrv(nt_prog)
    character(len=*), intent(in) :: tracer_name
    real, intent(in):: pe(isc-1:iec+1,npz+1,jsc-1:jec+1)
    real, intent(in)::  q(isd:ied,jsd:jed,npz, ncnst)
    real, intent(in)::area(isd:ied,jsd:jed)
    type(domain2D), intent(INOUT)  :: domain
    integer, intent(IN) :: npz, ncnst
!------------------------------------------------------------
    real psfc_sum(isc:iec,jsc:jec,1), qp_sum(isc:iec,jsc:jec,1)
    real qp, s1, s2
    integer j, i, k, idx

    psfc_sum = 0.
    qp_sum = 0.
    idx = get_tracer_index(MODEL_ATMOS, trim(tracer_name))

    if(idx /= NO_TRACER) then
   
!$OMP parallel do default(shared) private(qp)
       do j=jsc,jec
          do i=isc,iec
             psfc_sum(i,j,1) = pe(i,npz+1,j)*area(i,j)
!---------------------------------------------------------------------
!  define pressure-weighted column mean value of dry mass mixing
!  ratio  for tracer idx. assumption is that the tracer field q_phys
!  is a moist mass mixing ratio. convert to dry mass mixing ratio by
!  dividing by (1 - qh2o).
!---------------------------------------------------------------------
             qp = 0.0
             do k=1,npz
! SJL 2012/08/06
                qp = qp + q(i,j,k,idx)*(pe(i,k+1,j)-pe(i,k,j))/(1.-q(i,j,k,sphum))
             enddo
             qp_sum(i,j,1) = qp * area(i,j)
          enddo
       enddo
       s1 = REAL(mpp_global_sum(domain, psfc_sum, flags=NON_BITWISE_EXACT_SUM),KIND=4)
       s2 = REAL(mpp_global_sum(domain, qp_sum,   flags=NON_BITWISE_EXACT_SUM),KIND=4)
       rrv(idx) = s2 / s1
!---------------------------------------------------------------------
!    convert the tracer dry mass mixing ratio to the dry volume
!    mixing ratio.
!---------------------------------------------------------------------
       if (trim(tracer_name).eq.'co2') then
          rrv(idx) = rrv(idx)*WTMAIR/WTMCO2
       end if
    endif

  end subroutine compute_g_avg

end module fv_physics_mod
