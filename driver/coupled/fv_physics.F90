module fv_physics_mod

!-----------------------------------------------------------------------
!
!  Interface for Cubed_sphere FV dynamics with GFDL atmospheric physics
!  History: modified by SJL based on Memphis release
!-----------------------------------------------------------------------

!-----------------
! FMS modules:
!-----------------
use constants_mod,         only: rdgas, grav, rvgas
use time_manager_mod,      only: time_type
use fms_mod,               only: error_mesg, FATAL, write_version_number,clock_flag_default
use physics_driver_mod,    only: physics_driver_init, physics_driver_end,   &
                                 physics_driver_down, physics_driver_up, surf_diff_type
use field_manager_mod,     only: MODEL_ATMOS
use tracer_manager_mod,    only: get_tracer_index, NO_TRACER
use mpp_domains_mod,       only: mpp_global_sum, BITWISE_EXACT_SUM
use mpp_mod,               only: mpp_error, mpp_clock_id,  mpp_clock_begin,  &
                                 mpp_clock_end, CLOCK_MODULE_DRIVER, mpp_pe
#ifdef ATMOS_NUDGE
use atmos_nudge_mod,       only: atmos_nudge_init, atmos_nudge_end
#endif

!-----------------
! FV core modules:
!-----------------
use fv_grid_tools_mod,            only: area
use fv_arrays_mod,         only: fv_atmos_type
use fv_control_mod,            only: npx, npy, npz, ncnst, pnats, domain
use fv_eta_mod,               only: get_eta_level
use fv_update_phys_mod,    only: fv_update_phys
use fv_sg_mod,             only: fv_sg_conv
use fv_timing_mod,          only: timing_on, timing_off

implicit none
private

public  fv_physics_down, fv_physics_up, fv_physics_init, fv_physics_end
public  surf_diff_type

!-----------------------------------------------------------------------
character(len=128) :: version = '$Id: fv_physics.F90,v 1.1.2.8.4.2.2.1.2.1 2007/11/09 19:19:46 sjl Exp $'
character(len=128) :: tag = '$Name: omsk_2007_12 $'
!-----------------------------------------------------------------------

   real, allocatable, dimension(:,:,:)   :: u_dt, v_dt, t_dt
   real, allocatable, dimension(:,:,:,:) :: q_dt  ! potentially a huge array
   real, allocatable, dimension(:,:,:)   :: p_full, z_full, p_half, z_half
   logical :: do_atmos_nudge
   real    :: zvir, rrg, ginv
   integer :: id_fv_physics_down, id_fv_physics_up, id_fv_update_phys
   integer :: isc, iec, jsc, jec, ngc, nt_prog
   integer :: isd, ied, jsd, jed
   integer :: isw, iew, jsw, jew  ! window start/end in global index space
   integer :: nx_win, ny_win      ! iew-isw+1, jew-jsw+1 (window sizes)
   integer :: nx_dom, ny_dom      ! ie-is+1, je-js+1 (compute domain sizes)
   integer :: sphum

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

    real    ::  pref(npz+1,2)
    real    :: phalf(npz+1)
    real    :: ps1, ps2
    integer :: i, j, k

! All tracers are prognostic
    nt_prog = ncnst - pnats

    isc = Atm(1)%isc
    iec = Atm(1)%iec
    jsc = Atm(1)%jsc
    jec = Atm(1)%jec

    ngc = Atm(1)%ng

    isd = isc - ngc
    ied = iec + ngc
    jsd = jsc - ngc
    jed = jec + ngc

    zvir = rvgas/rdgas - 1.
    ginv = 1./ grav
    rrg  = rdgas / grav        

!----- write version to logfile --------
    call write_version_number (version,tag)

! Specific humidity is assumed to be q(:,:,:,1)
    sphum = get_tracer_index (MODEL_ATMOS, 'sphum' )
    if(sphum /= 1) call error_mesg('fv_physics_init:','sphum /= 1', FATAL)

!---------- reference profile -----------
    ps1 = 101325.
    ps2 =  81060.
    pref(npz+1,1) = ps1
    pref(npz+1,2) = ps2
    call get_eta_level ( npz, ps1, pref(1,1), phalf, Atm(1)%ak, Atm(1)%bk )
    call get_eta_level ( npz, ps2, pref(1,2), phalf, Atm(1)%ak, Atm(1)%bk )

    allocate (  u_dt(isd:ied,jsd:jed, npz) )
    allocate (  v_dt(isd:ied,jsd:jed, npz) )
    allocate (  t_dt(isc:iec,jsc:jec, npz) )
    allocate (  q_dt(isc:iec,jsc:jec, npz, nt_prog) )
    allocate (p_edge(isc:iec,jsc:jec, npz+1))

!------- pressure at model layer interfaces -----
    do k=1,npz+1
       do j=jsc,jec
          do i=isc,iec
             p_edge(i,j,k) = Atm(1)%pe(i,k,j)
          enddo
       enddo
    enddo
!---------- initialize physics -------

    call physics_driver_init(Time, Atm(1)%grid(isc:iec+1,jsc:jec+1,1),             &
                                   Atm(1)%grid(isc:iec+1,jsc:jec+1,2),             &
                             axes, pref, Atm(1)%q(isc:iec,jsc:jec,1:npz,1:ncnst),  &
                             Surf_diff,  p_edge )
    deallocate ( p_edge )

!--- initialize nudging module ---
#ifdef ATMOS_NUDGE
    call atmos_nudge_init ( Time, axes(1:3), flag=do_atmos_nudge )
#endif

! physics window
    nx_dom = iec - isc + 1
    ny_dom = jec - jsc + 1

    nx_win = window(1)
    ny_win = window(2)

    if( nx_win.LE.0 ) nx_win = nx_dom
    if( ny_win.LE.0 ) ny_win = ny_dom

! Consistency check:
    if( mod(nx_dom,nx_win).NE.0 )then
        write( text,'(a,2i4)' )'FV_PHYSICS_INIT: atmosphere_nml problem,'// &
             ' physics_window must divide domain size: ',  nx_win, nx_dom
        call mpp_error( FATAL, text )
    end if
    if( mod(ny_dom,ny_win).NE.0 )then
        write( text,'(a,2i4)' )'FV_PHYSICS_INIT: atmosphere_nml problem,'// &
             ' physics_window must divide domain size: ',  ny_win, ny_dom
        call mpp_error( FATAL, text )
    end if

    allocate( p_full(nx_win,ny_win,npz) )
    allocate( z_full(nx_win,ny_win,npz) )
    allocate( p_half(nx_win,ny_win,npz+1) )
    allocate( z_half(nx_win,ny_win,npz+1) )

!MPP clocks
    id_fv_physics_down = mpp_clock_id( 'FV_PHYSICS_DOWN', &
         flags=clock_flag_default, grain=CLOCK_MODULE_DRIVER )
    id_fv_physics_up = mpp_clock_id( 'FV_PHYSICS_UP', &
         flags=clock_flag_default, grain=CLOCK_MODULE_DRIVER )
    id_fv_update_phys = mpp_clock_id( 'FV_UPDATE_PHYS', &
         flags=clock_flag_default, grain=CLOCK_MODULE_DRIVER )

  end subroutine fv_physics_init


  subroutine fv_physics_down(Atm, dt_phys, Time_prev, Time, Time_next, &
                             frac_land,   albedo,            &
                             albedo_vis_dir, albedo_nir_dir, &
                             albedo_vis_dif, albedo_nir_dif, &
                             rough_vel,   t_surf,            &
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
                             gust, Surf_diff )
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
    type(fv_atmos_type), intent(inout) :: Atm(:)
    real,                intent(in) :: dt_phys
    real, intent(in), dimension(isc:iec,jsc:jec):: frac_land,  albedo, &
                                       albedo_vis_dir, albedo_nir_dir, &
                                       albedo_vis_dif, albedo_nir_dif, &
                                       rough_vel, t_surf, u_star,      &
                                       b_star, q_star, dtau_du, dtau_dv

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
    real :: rdt
    integer :: i, j, k, m
!---------------------------- do physics -------------------------------

    rdt = 1./ dt_phys
    gavg_rrv = 0.
    call compute_g_avg(gavg_rrv, 'co2', Atm(1)%pe, Atm(1)%q)
    call mpp_clock_begin(id_fv_physics_down)

    do jsw = jsc, jec, ny_win
       jew = jsw + ny_win - 1
       do isw = isc, iec, nx_win
          iew = isw + nx_win - 1

          if ( Atm(1)%fv_sg_adj>0 ) then
             call fv_sg_conv(isw, iew, jsw, jew, isd, ied, jsd, jed,    &
                             isc, iec, jsc, jec, npz, nt_prog, dt_phys, &
                             Atm(1)%fv_sg_adj, Atm(1)%delp, Atm(1)%pe, Atm(1)%peln,  &
                             Atm(1)%pt, Atm(1)%q, Atm(1)%ua, Atm(1)%va,       &
                             u_dt, v_dt, t_dt, q_dt, Atm(1)%ak, Atm(1)%bk) 
             do k=1,npz
                do j=jsw,jew
                   do i=isw,iew
                      u_dt(i,j,k) = (u_dt(i,j,k) - Atm(1)%ua(i,j,k)) * rdt
                      v_dt(i,j,k) = (v_dt(i,j,k) - Atm(1)%va(i,j,k)) * rdt
                      t_dt(i,j,k) = (t_dt(i,j,k) - Atm(1)%pt(i,j,k)) * rdt
                   enddo
                enddo
             enddo
             do m=1,nt_prog
             do k=1,npz
                do j=jsw,jew
                   do i=isw,iew
                      q_dt(i,j,k,m) = (q_dt(i,j,k,m) - Atm(1)%q(i,j,k,m)) * rdt
                   enddo
                enddo
             enddo
             enddo
          else
! Initialize tendencies due to parameterizations:
             do k=1,npz
                do j=jsw,jew
                   do i=isw,iew
                      u_dt(i,j,k) = 0.
                      v_dt(i,j,k) = 0.
                      t_dt(i,j,k) = 0.
                   enddo
                enddo
             enddo
             do m=1,nt_prog
             do k=1,npz
                do j=jsw,jew
                   do i=isw,iew
                      q_dt(i,j,k,m) = 0.
                   enddo
                enddo
             enddo
             enddo
          endif

          call compute_p_z(npz, isw, jsw, nx_win, ny_win, Atm(1)%phis, Atm(1)%pt, &
                           Atm(1)%q, Atm(1)%delp, Atm(1)%pe, Atm(1)%peln,         &
                           Atm(1)%delz,  Atm(1)%hydrostatic)

          call physics_driver_down( isw-isc+1, iew-isc+1, jsw-jsc+1, jew-jsc+1, &
                   Time_prev, Time, Time_next                              , &
                   Atm(1)%agrid(isw:iew,jsw:jew,2)                         , & 
                   Atm(1)%agrid(isw:iew,jsw:jew,1)                         , & 
                   area(isw:iew,jsw:jew), p_half,  p_full, z_half,  z_full , &
                   Atm(1)%ua(isw:iew,jsw:jew,:)                            , &
                   Atm(1)%va(isw:iew,jsw:jew,:)                            , &
                   Atm(1)%pt(isw:iew,jsw:jew,:)                            , &
                   Atm(1)%q (isw:iew,jsw:jew,:,1)                          , &
                   Atm(1)%q (isw:iew,jsw:jew,:,:)                          , &
                   Atm(1)%ua(isw:iew,jsw:jew,:)                            , &
                   Atm(1)%va(isw:iew,jsw:jew,:)                            , &
                   Atm(1)%pt(isw:iew,jsw:jew,:)                            , &
                   Atm(1)%q (isw:iew,jsw:jew,:,1)                          , &
                   Atm(1)%q (isw:iew,jsw:jew,:,:)                          , &
                   frac_land(isw:iew,jsw:jew), rough_vel(isw:iew,jsw:jew)  , &
                   albedo   (isw:iew,jsw:jew)                              , &
                   albedo_vis_dir(isw:iew,jsw:jew)                         , &
                   albedo_nir_dir(isw:iew,jsw:jew)                         , &
                   albedo_vis_dif(isw:iew,jsw:jew)                         , &
                   albedo_nir_dif(isw:iew,jsw:jew)                         , &
                   t_surf  (isw:iew,jsw:jew),  u_star(isw:iew,jsw:jew)     , &
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
                   Surf_diff,   gavg_rrv )
       enddo
    enddo
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
    type(fv_atmos_type), intent(inout) :: Atm(:)
    type(surf_diff_type),intent(inout) :: Surf_diff
    real,                intent(in)    :: dt_phys
    real, intent(in),    dimension(isc:iec,jsc:jec) :: frac_land
    real, intent(inout), dimension(isc:iec,jsc:jec) :: gust
    real, intent(out),   dimension(isc:iec,jsc:jec) :: lprec, fprec
    real, intent(in),    dimension(isc:iec,jsc:jec) :: u_star, b_star, q_star

    call mpp_clock_begin(id_fv_physics_up)

    do jsw = jsc, jec, ny_win
       jew = jsw + ny_win - 1
       do isw = isc,iec,nx_win
          iew = isw + nx_win - 1

          call compute_p_z(npz, isw, jsw, nx_win, ny_win, Atm(1)%phis, Atm(1)%pt,   &
                           Atm(1)%q, Atm(1)%delp, Atm(1)%pe, Atm(1)%peln,         &
                           Atm(1)%delz,  Atm(1)%hydrostatic)

          call physics_driver_up( isw-isc+1, iew-isc+1, jsw-jsc+1, jew-jsc+1, &
                                  Time_prev, Time, Time_next             , &
                                  Atm(1)%agrid(isw:iew,jsw:jew,2)        , &
                                  Atm(1)%agrid(isw:iew,jsw:jew,1)        , &
                                  area(isw:iew,jsw:jew), p_half, p_full  , &
                                  z_half,  z_full                        , &
                                  Atm(1)%omga(isw:iew,jsw:jew,:)         , &
                                  Atm(1)%ua(isw:iew,jsw:jew,:)           , &
                                  Atm(1)%va(isw:iew,jsw:jew,:)           , &
                                  Atm(1)%pt(isw:iew,jsw:jew,:)           , &
                                  Atm(1)%q(isw:iew,jsw:jew,:,1)          , &
                                  Atm(1)%q(isw:iew,jsw:jew,:,1:nt_prog)  , &
                                  Atm(1)%ua(isw:iew,jsw:jew,:)           , &
                                  Atm(1)%va(isw:iew,jsw:jew,:)           , &
                                  Atm(1)%pt(isw:iew,jsw:jew,:)           , &
                                  Atm(1)%q(isw:iew,jsw:jew,:,1)          , &
                                  Atm(1)%q(isw:iew,jsw:jew,:,1:nt_prog)  , &
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
                                  gust     (isw:iew,jsw:jew) )
       enddo
    enddo
    call mpp_clock_end(id_fv_physics_up)

    call mpp_clock_begin(id_fv_update_phys)
                                                            call timing_on('update_fv')
    call fv_update_phys( dt_phys,   isc,        iec,         jsc,    jec,   isd,       &
                         ied,       jsd,        jed,         ngc,       nt_prog,       &
                         Atm(1)%u,  Atm(1)%v,   Atm(1)%delp, Atm(1)%pt, Atm(1)%q,      &
                         Atm(1)%ua, Atm(1)%va,  Atm(1)%ps,   Atm(1)%pe, Atm(1)%peln,   &
                         Atm(1)%pk, Atm(1)%pkz, Atm(1)%ak,   Atm(1)%bk, u_dt,          &
                         v_dt, t_dt, q_dt, Atm(1)%u_srf, Atm(1)%v_srf,  Atm(1)%delz,   &
                         Atm(1)%hydrostatic, .true., Time_next )
                                                            call timing_off('update_fv')
    call mpp_clock_end(id_fv_update_phys)

  end subroutine fv_physics_up



  subroutine fv_physics_end (Time)
   type(time_type), intent(in) :: Time
!                                 NOTE: this is not the dynamics time
    call physics_driver_end (Time)
#ifdef ATMOS_NUDEG
    call atmos_nudge_end
#endif

    deallocate ( u_dt )
    deallocate ( v_dt )
    deallocate ( t_dt )
    deallocate ( q_dt )
    deallocate ( p_full )
    deallocate ( z_full )
    deallocate ( p_half )
    deallocate ( z_half )

  end subroutine fv_physics_end



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
  real,    intent(in):: delz(isc:iec,jsc:jec,nlev)
  logical, intent(in):: hydrostatic
! local
  integer i,j,k,id,jd
  real    tvm

!----------------------------------------------------
! Compute pressure and height at full and half levels
!----------------------------------------------------
  do j=1,jsiz
     jd = j + jstart - 1
     do i=1,isiz
        id = i + istart - 1
        z_half(i,j,nlev+1) = phis(id,jd) * ginv
     enddo
  end do

  do k=1,nlev+1
     do j=1,jsiz
        jd = j + jstart - 1
        do i=1,isiz
           id = i + istart - 1
           p_half(i,j,k) = pe(id,k,jd)
        enddo
     enddo
  enddo

#ifdef TEST_NH
  if ( hydrostatic ) then
#endif
  do k=nlev,1,-1
     do j=1,jsiz
        jd = j + jstart - 1
        do i=1,isiz
           id = i + istart - 1
           tvm = rrg*pt(id,jd,k)*(1.+zvir*q(id,jd,k,sphum))
           p_full(i,j,k) = delp(id,jd,k)/(peln(id,k+1,jd)-peln(id,k,jd))
           z_full(i,j,k) = z_half(i,j,k+1) + tvm*(1.-p_half(i,j,k)/p_full(i,j,k))
           z_half(i,j,k) = z_half(i,j,k+1) + tvm*(peln(id,k+1,jd)-peln(id,k,jd))
        enddo
     enddo
  enddo
#ifdef TEST_NH
  else
  do k=nlev,1,-1
     do j=1,jsiz
        jd = j + jstart - 1
        do i=1,isiz
           id = i + istart - 1
! Consistency between p_half and P_full? obtain p_half from Riemman solver?
! OR, use hydrostatic values?
           p_full(i,j,k) = -rrg*delp(id,jd,k)/delz(id,jd,k)*pt(id,jd,k)*(1.+zvir*q(id,jd,k,sphum))
           z_half(i,j,k) = z_half(i,j,k+1) - delz(i,j,k)
           z_full(i,j,k) = 0.5*(z_half(i,j,k) + z_half(i,j,k+1))
        enddo
     enddo
  enddo
  endif
#endif

  end subroutine compute_p_z



  subroutine compute_g_avg(rrv, tracer_name, pe, q)

  real,          intent(inout) :: rrv(nt_prog)
  character(len=*), intent(in) :: tracer_name
  real, intent(in):: pe(isc-1:iec+1,npz+1,jsc-1:jec+1)
  real, intent(in)::  q(isd:ied,jsd:jed,npz, ncnst)
!------------------------------------------------------------
  real psfc_sum(isc:iec,jsc:jec,1), qp_sum(isc:iec,jsc:jec,1)
  real qp, s1, s2
  integer j, i, k, idx

  psfc_sum = 0.
  qp_sum = 0.
  idx = get_tracer_index(MODEL_ATMOS, trim(tracer_name))

  if(idx /= NO_TRACER) then
     do j=jsc,jec
        do i=isc,iec
           psfc_sum(i,j,1) = pe(i,npz+1,j)*area(i,j)
           qp = 0.0
           do k=1,npz
              qp = qp + q(i,j,k,idx)*(pe(i,k+1,j) - pe(i,k,j))
           enddo
           qp_sum(i,j,1) = qp * area(i,j)
        enddo
     enddo
     s1 = mpp_global_sum(domain, psfc_sum, flags=BITWISE_EXACT_SUM)
     s2 = mpp_global_sum(domain, qp_sum,   flags=BITWISE_EXACT_SUM)
     rrv(idx) = s2 / s1
  endif
  
  end subroutine compute_g_avg

end module fv_physics_mod
