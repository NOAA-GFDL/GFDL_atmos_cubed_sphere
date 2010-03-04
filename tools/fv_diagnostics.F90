module fv_diagnostics_mod

 use constants_mod,    only: grav, rdgas, pi, radius, kappa
 use fms_io_mod,       only: set_domain, nullify_domain
 use time_manager_mod, only: time_type, get_date, get_time
 use mpp_domains_mod,  only: domain2d, mpp_update_domains, DGRID_NE
 use diag_manager_mod, only: diag_axis_init, register_diag_field, &
                             register_static_field, send_data, diag_grid_init
 use fv_arrays_mod,    only: fv_atmos_type
 use fv_mapz_mod,      only: E_Flux
 use fv_mp_mod,        only: domain, gid, masterproc, &
                             mp_reduce_sum, mp_reduce_min, mp_reduce_max
 use fv_eta_mod,        only: get_eta_level, gw_1d
 use fv_grid_tools_mod, only: dx, dy, rdxa, rdya, area, rarea
 use fv_grid_utils_mod, only: f0, cosa_s, g_sum, sina_u, sina_v, en1, en2, vlon
 use a2b_edge_mod,     only: a2b_ord4
 use fv_surf_map_mod,  only: zs_g
 use fv_sg_mod,        only: qsmith

 use tracer_manager_mod, only: get_tracer_names, get_number_tracers, get_tracer_index
 use field_manager_mod,  only: MODEL_ATMOS
 use mpp_mod,            only: mpp_error, FATAL, stdlog

#if defined(MARS_GCM) && defined(MARS_SURFACE)
 use mars_surface_mod,     only:  sfc_snow, sfc_frost
#  ifdef DUST_SOURCE
 use dust_source_mod,      only:  sfc_dust, odcol
!!!  use aerosol_mod,          only:  ndust_bins, nice_bins, nice_moms, aerosol_bins, &
!!!                                   dust_indx, ice_bin_indx, ice_mom_indx
#  endif DUST_SOURCE

#  ifdef WATER_CYCLE
 use cloud_physics_mod,   only:  cldcol, wcol
#  endif WATER_CYCLE
#endif


 implicit none
 private

 integer ::id_ps, id_slp, id_ua, id_va, id_pt, id_omga, id_vort,  &
           id_tm, id_pv, id_zsurf, id_oro, id_sgh, id_divg, id_w, &
           id_te, id_zs, id_ze, id_mq, id_vorts, id_us, id_vs,    &
           id_tq, id_rh, id_c15, id_c25, id_c35, id_c45,          &
                         id_f15, id_f25, id_f35, id_f45,          &
           id_ppt, id_ts

! Selected p-level fields from 3D variables:
 integer ::id_h300, id_h500, id_vort850, id_u850, id_v850, id_w850,  &
           id_u200, id_v200, id_w200, id_s200, id_sl12, id_sl13
 integer ::id_t200, id_t500, id_t700, id_t850,              &
           id_q200, id_q500, id_q850,                       &
           id_omg200, id_omg500, id_omg850, id_h200,        &
           id_h50, id_t50, id_q50, id_u50, id_v50
! Ned diag for HFIP experiments:
 integer :: id_h850, id_u700, id_v700, id_h700, id_u500, id_v500

#ifdef MARS_GCM
 integer ::  id_t05
 integer ::  id_tdust, id_sfc_dust
#endif MARS_GCM

 integer, parameter:: max_step = 1000
 integer steps
 real(kind=4):: efx(max_step), mtq(max_step)
 real(kind=4):: efx_sum,       mtq_sum
! For initial conditions:
 integer ic_ps, ic_ua, ic_va, ic_ppt
 integer, allocatable :: id_tracer(:)

 integer  :: ncnst
 real :: missing_value = -1.e10
 real :: ginv
 real, allocatable :: phalf(:)
 real, allocatable :: zsurf(:,:)
 real, allocatable :: zxg(:,:)
 real, allocatable :: pt1(:)
 real :: pk0
 logical master

 type(time_type) :: fv_time

 logical :: module_is_initialized=.false.
 logical :: moist_phys
 integer  sphum, liq_wat, ice_wat       ! GFDL physics
 integer  rainwat, snowwat, graupel
 real    :: ptop
 real    :: rad2deg
! tracers
 character(len=128)   :: tname
 character(len=256)   :: tlongname, tunits

 public :: fv_diag_init, fv_time, fv_diag, prt_maxmin, range_check, id_divg, id_te
 public :: efx, efx_sum, mtq, mtq_sum, steps

contains

 subroutine fv_diag_init(Atm, axes, Time, npx, npy, npz, p_ref)
    type(fv_atmos_type), intent(inout) :: Atm(:)
    integer, intent(out) :: axes(4)
    type(time_type), intent(in) :: Time
    integer,         intent(in) :: npx, npy, npz
    real, intent(in):: p_ref

    real, allocatable :: grid_xt(:), grid_yt(:), grid_xe(:), grid_ye(:), grid_xn(:), grid_yn(:)
    real, allocatable :: grid_x(:),  grid_y(:)
    real              :: vrange(2), vsrange(2), wrange(2), trange(2), slprange(2), rhrange(2)
    real, allocatable :: a3(:,:,:)
    real              :: pfull(npz)
    real              :: hyam(npz), hybm(npz)

    integer :: id_bk, id_pk, id_area, id_lon, id_lat, id_lont, id_latt, id_phalf, id_pfull
    integer :: id_hyam, id_hybm
    integer :: i, j, k, n, ntileMe, id_xt, id_yt, id_x, id_y, id_xe, id_ye, id_xn, id_yn
    integer :: isc, iec, jsc, jec

    logical :: used

    character(len=64) :: field
    integer              :: ntprog
    integer              :: unit

    rad2deg = 180./pi

! For total energy diagnostics:
    steps = 0
    efx = 0.;       efx_sum = 0.
    mtq = 0.;       mtq_sum = 0.

    ncnst = Atm(1)%ncnst
    moist_phys = Atm(1)%moist_phys

    call set_domain(Atm(1)%domain)  ! Set domain so that diag_manager can access tile information

    if ( Atm(1)%nwat>=3 ) then
         sphum   = get_tracer_index (MODEL_ATMOS, 'sphum')
         liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
         ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
    endif

    if ( Atm(1)%nwat==6 ) then
        rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat')
        snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat')
        graupel = get_tracer_index (MODEL_ATMOS, 'graupel')
    else
         sphum   = 1
     endif

! valid range for some fields

!!!  This will need mods for more than 1 tile per pe  !!!

    vsrange = (/ -200.,  200. /)  ! surface (lowest layer) winds

    vrange = (/ -330.,  330. /)  ! winds
    wrange = (/  -80.,   50. /)  ! vertical wind
   rhrange = (/  -10.,  150. /)  ! RH

#if defined(MARS_GCM)
    slprange = (/0.,  100./)  ! sea-level-pressure
    trange = (/  50., 360. /)  ! temperature
#elif defined(VENUS_GCM)
    trange = (/  100.,  900. /)  ! temperature
    slprange = (/80.E3,  98.E3/)  ! sea-level-pressure
#else
    trange = (/  100.,  350. /)  ! temperature
    slprange = (/800.,  1200./)  ! sea-level-pressure
#endif

    ginv = 1./GRAV
    fv_time = Time

    allocate ( phalf(npz+1) )
    call get_eta_level(Atm(1)%npz, p_ref, pfull, phalf, Atm(1)%ak, Atm(1)%bk, 0.01)

!   allocate(grid_xt(npx-1), grid_yt(npy-1), grid_xe(npx), grid_ye(npy-1), grid_xn(npx-1), grid_yn(npy))
    allocate(grid_xt(npx-1), grid_yt(npy-1))
    grid_xt = (/ (i, i=1,npx-1) /)
    grid_yt = (/ (j, j=1,npy-1) /)
!   grid_xe = (/ (i, i=1,npx) /)
!   grid_ye = (/ (j, j=1,npy-1) /)
!   grid_xn = (/ (i, i=1,npx-1) /)
!   grid_yn = (/ (j, j=1,npy) /)

    allocate(grid_x(npx), grid_y(npy))
    grid_x = (/ (i, i=1,npx) /)
    grid_y = (/ (j, j=1,npy) /)

    n=1
    isc = Atm(n)%isc; iec = Atm(n)%iec
    jsc = Atm(n)%jsc; jec = Atm(n)%jec

    ! Send diag_manager the grid informtaion
    call diag_grid_init(DOMAIN=Atm(n)%domain, &
         &              GLO_LON=rad2deg*Atm(n)%grid(isc:iec+1,jsc:jec+1,1), &
         &              GLO_LAT=rad2deg*Atm(n)%grid(isc:iec+1,jsc:jec+1,2))

    ntileMe = size(Atm(:))
    do n = 1, ntileMe
       field = 'grid'

       id_xt = diag_axis_init('grid_xt',grid_xt,'degrees_E','x','T-cell longitude', &
                           set_name=trim(field),Domain2=Domain, tile_count=n)
       id_yt = diag_axis_init('grid_yt',grid_yt,'degrees_N','y','T-cell latitude',  &
                           set_name=trim(field), Domain2=Domain, tile_count=n)
!  Don't need these right now
!      id_xe = diag_axis_init ('grid_xe',grid_xe,'degrees_E','x','E-cell longitude', &
!                              set_name=trim(field),Domain2=Domain, tile_count=n)
!      id_ye = diag_axis_init ('grid_ye',grid_ye,'degrees_N','y','E-cell latitude',  &
!                              set_name=trim(field), Domain2=Domain, tile_count=n)
!      id_xn = diag_axis_init ('grid_xn',grid_xn,'degrees_E','x','N-cell longitude', &
!                              set_name=trim(field),Domain2=Domain, aux='geolon_n, geolat_n', tile_count=n)
!      id_yn = diag_axis_init ('grid_yn',grid_yn,'degrees_N','y','N-cell latitude',  &
!                              set_name=trim(field), Domain2=Domain, tile_count=n)

       id_x = diag_axis_init('grid_x',grid_x,'degrees_E','x','cell corner longitude', &
                           set_name=trim(field),Domain2=Domain, tile_count=n)
       id_y = diag_axis_init('grid_y',grid_y,'degrees_N','y','cell corner latitude',  &
                           set_name=trim(field), Domain2=Domain, tile_count=n)

    end do
!   deallocate(grid_xt, grid_yt, grid_xe, grid_ye, grid_xn, grid_yn)
    deallocate(grid_xt, grid_yt)
    deallocate(grid_x,  grid_y )

    id_phalf = diag_axis_init('phalf', phalf, 'mb', 'z', &
            'ref half pressure level', direction=-1, set_name="dynamics")
    id_pfull = diag_axis_init('pfull', pfull, 'mb', 'z', &
            'ref full pressure level', direction=-1, set_name="dynamics", edges=id_phalf)

!---- register static fields -------

    id_bk    = register_static_field ( "dynamics", 'bk', (/id_phalf/), &
         'vertical coordinate sigma value', 'none' )

    id_pk    = register_static_field ( "dynamics", 'pk', (/id_phalf/), &
         'pressure part of the hybrid coordinate', 'pascal' )

    id_hyam    = register_static_field ( "dynamics", 'hyam', (/id_pfull/), &
         'vertical coordinate A value', '1E-5 Pa' )

    id_hybm    = register_static_field ( "dynamics", 'hybm', (/id_pfull/), &
         'vertical coordinate B value', 'none' )

!--- Send static data

    if ( id_bk > 0 )    used = send_data ( id_bk,Atm(1)%bk, Time )
    if ( id_pk > 0 )    used = send_data ( id_pk,Atm(1)%ak, Time )
    if ( id_hyam > 0 ) then
         do k=1,npz
            hyam(k) = 0.5 * ( Atm(1)%ak(k) + Atm(1)%ak(k+1) ) * 1.E-5
         enddo
         used = send_data ( id_hyam, hyam, Time )
    endif
    if ( id_hybm > 0 ) then
         do k=1,npz
            hybm(k) = 0.5 * ( Atm(1)%bk(k) + Atm(1)%bk(k+1) )
         enddo
         used = send_data ( id_hybm, hybm, Time )
    endif

!   Approach will need modification if we wish to write values on other than A grid.
    axes(1) = id_xt
    axes(2) = id_yt
    axes(3) = id_pfull
    axes(4) = id_phalf

!---- register time independent fields -------

    do n = 1, ntileMe
       field= 'dynamics'
       id_lon  = register_static_field ( trim(field), 'grid_lon', (/id_x,id_y/),  &
                                         'longitude', 'degrees_E' )
       id_lat  = register_static_field ( trim(field), 'grid_lat', (/id_x,id_y/),  &
                                         'latitude', 'degrees_N' )
       id_lont = register_static_field ( trim(field), 'grid_lont', (/id_xt,id_yt/),  &
                                         'longitude', 'degrees_E' )
       id_latt = register_static_field ( trim(field), 'grid_latt', (/id_xt,id_yt/),  &
                                         'latitude', 'degrees_N' )
       id_area = register_static_field ( trim(field), 'area', axes(1:2),  &
                                         'cell area', 'm**2' )
#ifndef DYNAMICS_ZS
       id_zsurf = register_static_field ( trim(field), 'zsurf', axes(1:2),  &
                                         'surface height', 'm' )
#endif
       id_zs = register_static_field ( trim(field), 'zs', axes(1:2),  &
                                        'Original Mean Terrain', 'm' )
! 3D hybrid_z fields:
       id_ze = register_static_field ( trim(field), 'ze', axes(1:3),  &
                                        'Hybrid_Z_surface', 'm' )
! For mountain torque in zonal dir:
       id_oro = register_static_field ( trim(field), 'oro', axes(1:2),  &
                                        'Land/Water Mask', 'none' )
       id_sgh = register_static_field ( trim(field), 'sgh', axes(1:2),  &
                                        'Terrain Standard deviation', 'm' )
       id_ts = register_static_field ( trim(field), 'ts', axes(1:2),  &
                                        'Skin temperature', 'K' )

!--------------------
! Initial conditions:
!--------------------
       ic_ps  = register_static_field ( trim(field), 'ps_ic', axes(1:2),  &
                                         'initial surface pressure', 'Pa' )
       ic_ua = register_static_field ( trim(field), 'ua_ic', axes(1:3),        &
            'zonal wind', 'm/sec' )
       ic_va = register_static_field ( trim(field), 'va_ic', axes(1:3),        &
            'meridional wind', 'm/sec' )
       ic_ppt= register_static_field ( trim(field), 'ppt_ic', axes(1:3),        &
            'potential temperature perturbation', 'K' )

    end do

    master = (gid == masterproc)

    n = 1 

    allocate ( zsurf(isc:iec,jsc:jec) )

    do j=jsc,jec
       do i=isc,iec
          zsurf(i,j) = ginv * Atm(n)%phis(i,j)
       enddo
    enddo

!--- Send time independent data

    do n = 1, ntileMe
       isc = Atm(n)%isc; iec = Atm(n)%iec
       jsc = Atm(n)%jsc; jec = Atm(n)%jec

       if (id_lon  > 0) used = send_data(id_lon,  180./pi*Atm(n)%grid(isc:iec+1,jsc:jec+1,1), Time)
       if (id_lat  > 0) used = send_data(id_lat,  180./pi*Atm(n)%grid(isc:iec+1,jsc:jec+1,2), Time)
       if (id_lont > 0) used = send_data(id_lont, 180./pi*Atm(n)%agrid(isc:iec,jsc:jec,1), Time)
       if (id_latt > 0) used = send_data(id_latt, 180./pi*Atm(n)%agrid(isc:iec,jsc:jec,2), Time)
       if (id_area > 0) used = send_data(id_area, area(isc:iec,jsc:jec), Time)
#ifndef DYNAMICS_ZS
       if (id_zsurf > 0) used = send_data(id_zsurf, zsurf, Time)
#endif
       if ( Atm(n)%fv_land ) then
         if (id_zs  > 0) used = send_data(id_zs , zs_g, Time)
         if (id_oro > 0) used = send_data(id_oro, Atm(n)%oro(isc:iec,jsc:jec), Time)
         if (id_sgh > 0) used = send_data(id_sgh, Atm(n)%sgh(isc:iec,jsc:jec), Time)
       endif

       if ( Atm(n)%ncep_ic ) then
         if (id_ts > 0) used = send_data(id_ts, Atm(n)%ts(isc:iec,jsc:jec), Time)
       endif

       if ( Atm(n)%hybrid_z .and. id_ze > 0 ) &
                      used = send_data(id_ze, Atm(n)%ze0(isc:iec,jsc:jec,1:npz), Time)

       if (ic_ps > 0) used = send_data(ic_ps, Atm(n)%ps(isc:iec,jsc:jec)*ginv, Time)

       if(ic_ua > 0) used=send_data(ic_ua, Atm(n)%ua(isc:iec,jsc:jec,:), Time)
       if(ic_va > 0) used=send_data(ic_va, Atm(n)%va(isc:iec,jsc:jec,:), Time)

       pk0 = 1000.E2 ** kappa
       if(ic_ppt> 0) then
! Potential temperature
          allocate ( pt1(npz) )
          allocate ( a3(isc:iec,jsc:jec,npz) )
#ifdef TEST_GWAVES
          call gw_1d(npz, 1000.E2, Atm(n)%ak, Atm(n)%ak, Atm(n)%ak(1), 10.E3, pt1)
#else
          pt1 = 0.
#endif
          do k=1,npz
          do j=jsc,jec
             do i=isc,iec
                a3(i,j,k) =  (Atm(n)%pt(i,j,k)/Atm(n)%pkz(i,j,k) - pt1(k)) * pk0
             enddo
          enddo
          enddo
          used=send_data(ic_ppt, a3, Time)
          deallocate ( a3 )
          deallocate ( pt1 )
       endif
    end do

!--------------------------------------------------------------
! Register main prognostic fields: ps, (u,v), t, omega (dp/dt)
!--------------------------------------------------------------

    allocate(id_tracer(ncnst))

    do n = 1, ntileMe
       field= 'dynamics'

#ifdef DYNAMICS_ZS
       id_zsurf = register_diag_field ( trim(field), 'zsurf', axes(1:2), Time,           &
                                       'surface height', 'm')
#endif
!-------------------
! Surface pressure
!-------------------
       id_ps = register_diag_field ( trim(field), 'ps', axes(1:2), Time,           &
            'surface pressure', 'Pa', missing_value=missing_value )

!-------------------
! Mountain torque
!-------------------
       id_mq = register_diag_field ( trim(field), 'mq', axes(1:2), Time,           &
            'mountain torque', 'Hadleys per unit area', missing_value=missing_value )

!--------------
! 50 mb Height
!--------------
      id_h50 = register_diag_field (trim(field), 'h50', axes(1:2),  Time,   &
                                     '50-mb hght', 'm', missing_value=missing_value )
!--------------
! 200 mb Height
!--------------
      id_h200 = register_diag_field (trim(field), 'h200', axes(1:2),  Time,   &
                                     '200-mb hght', 'm', missing_value=missing_value )
!--------------
! 300 mb Height
!--------------
      id_h300 = register_diag_field (trim(field), 'h300', axes(1:2),  Time,   &
                                     '300-mb hght', 'm', missing_value=missing_value )
!--------------
! 500 mb Height
!--------------
      id_h500 = register_diag_field (trim(field), 'h500', axes(1:2),  Time,   &
                                     '500-mb hght', 'm', missing_value=missing_value )
!--------------
! 700 mb Height
!--------------
      id_h700 = register_diag_field (trim(field), 'h700', axes(1:2),  Time,   &
                                     '700-mb hght', 'm', missing_value=missing_value )
!--------------
! 850 mb Height
!--------------
      id_h850 = register_diag_field (trim(field), 'h850', axes(1:2),  Time,   &
                                     '850-mb hght', 'm', missing_value=missing_value )
!-----------------------------
! mean temp between 300-500 mb
!-----------------------------
      id_tm = register_diag_field (trim(field), 'tm', axes(1:2),  Time,   &
                                   'mean 300-500 mb temp', 'K', missing_value=missing_value )

!-------------------
! Sea-level-pressure
!-------------------
       id_slp = register_diag_field (trim(field), 'slp', axes(1:2),  Time,   &
                                     'sea-level pressure', 'mb', missing_value=missing_value,  &
                                      range=slprange )
!-------------------
! Hurricane scales:
!-------------------
! Net effects: ~ intensity * freq
       id_c15 = register_diag_field (trim(field), 'cat15', axes(1:2),  Time,   &
                                     'de-pression < 1000', 'mb', missing_value=missing_value)
       id_c25 = register_diag_field (trim(field), 'cat25', axes(1:2),  Time,   &
                                     'de-pression < 980', 'mb', missing_value=missing_value)
       id_c35 = register_diag_field (trim(field), 'cat35', axes(1:2),  Time,   &
                                     'de-pression < 964', 'mb', missing_value=missing_value)
       id_c45 = register_diag_field (trim(field), 'cat45', axes(1:2),  Time,   &
                                     'de-pression < 944', 'mb', missing_value=missing_value)
! Frequency:
       id_f15 = register_diag_field (trim(field), 'f15', axes(1:2),  Time,   &
                                     'Cat15 frequency', 'none', missing_value=missing_value)
       id_f25 = register_diag_field (trim(field), 'f25', axes(1:2),  Time,   &
                                     'Cat25 frequency', 'none', missing_value=missing_value)
       id_f35 = register_diag_field (trim(field), 'f35', axes(1:2),  Time,   &
                                     'Cat35 frequency', 'none', missing_value=missing_value)
       id_f45 = register_diag_field (trim(field), 'f45', axes(1:2),  Time,   &
                                     'Cat45 frequency', 'none', missing_value=missing_value)
!-------------------
! A grid winds (lat-lon)
!-------------------
       id_ua = register_diag_field ( trim(field), 'ucomp', axes(1:3), Time,        &
            'zonal wind', 'm/sec', missing_value=missing_value, range=vrange )
       id_va = register_diag_field ( trim(field), 'vcomp', axes(1:3), Time,        &
            'meridional wind', 'm/sec', missing_value=missing_value, range=vrange)

       id_w = register_diag_field ( trim(field), 'w', axes(1:3), Time,        &
            'vertical wind', 'm/sec', missing_value=missing_value, range=wrange )

       id_pt   = register_diag_field ( trim(field), 'temp', axes(1:3), Time,       &
            'temperature', 'K', missing_value=missing_value, range=trange )
       id_ppt  = register_diag_field ( trim(field), 'ppt', axes(1:3), Time,       &
            'potential temperature perturbation', 'K', missing_value=missing_value )
       id_omga = register_diag_field ( trim(field), 'omega', axes(1:3), Time,      &
            'omega', 'Pa/s', missing_value=missing_value )
       id_divg  = register_diag_field ( trim(field), 'divg', axes(1:3), Time,      &
            'mean divergence', '1/s', missing_value=missing_value )

       id_rh = register_diag_field ( trim(field), 'rh', axes(1:3), Time,        &
            'Relative Humidity', '%', missing_value=missing_value, range=rhrange )
! Total energy (only when moist_phys = .T.)
       id_te    = register_diag_field ( trim(field), 'te', axes(1:2), Time,      &
            'Total Energy', 'J/kg', missing_value=missing_value )
!--------------------
! Relative vorticity
!--------------------
       id_vort = register_diag_field ( trim(field), 'vort', axes(1:3), Time,       &
            'vorticity', '1/s', missing_value=missing_value )
!--------------------
! Potential vorticity
!--------------------
       id_pv = register_diag_field ( trim(field), 'pv', axes(1:3), Time,       &
            'potential vorticity', '1/s', missing_value=missing_value )

#ifdef MARS_GCM
!--------------------------
! Extra Martian diagnostics:
!--------------------------

       id_t05 = register_diag_field ( trim(field), 't05', axes(1:2), Time,       &
               '0.5-mb temperature', 'K', missing_value=missing_value )
!!       id_sfc_dust = register_diag_field ( trim(field), 'sfc_dust', axes(1:2), Time,        &
!!             'Total sfc dust', 'kg/m**2', missing_value=missing_value )
!!        id_tdust = register_diag_field ( trim(field), 'odcol', axes(1:2), Time,        &
!!             'Total dust column', 'kg/m**2', missing_value=missing_value )
#endif MARS_GCM

!--------------------------
! Extra surface diagnistics:
!--------------------------
! Surface (lowest layer) vorticity: for tropical cyclones diag.
       id_vorts = register_diag_field ( trim(field), 'vorts', axes(1:2), Time,       &
            'surface vorticity', '1/s', missing_value=missing_value )
       id_us = register_diag_field ( trim(field), 'us', axes(1:2), Time,        &
            'surface u-wind', 'm/sec', missing_value=missing_value, range=vsrange )
       id_vs = register_diag_field ( trim(field), 'vs', axes(1:2), Time,        &
            'surface v-wind', 'm/sec', missing_value=missing_value, range=vsrange )
#ifdef OLD_TQ
       id_tq = register_diag_field ( trim(field), 'tq', axes(1:2), Time,        &
            'Total water vapor', 'kg/m**2', missing_value=missing_value )
#else
       id_tq = register_diag_field ( trim(field), 'tq', axes(1:2), Time,        &
            'Total water path', 'kg/m**2', missing_value=missing_value )
#endif

!--------------------------
! 850-mb vorticity
!--------------------------
       id_vort850 = register_diag_field ( trim(field), 'vort850', axes(1:2), Time,       &
                           '850-mb vorticity', '1/s', missing_value=missing_value )

!--------------------------
! 50-mb winds:
!--------------------------
       id_u50 = register_diag_field ( trim(field), 'u50', axes(1:2), Time,       &
                           '50-mb u-wind', '1/s', missing_value=missing_value )
       id_v50 = register_diag_field ( trim(field), 'v50', axes(1:2), Time,       &
                           '50-mb v-wind', '1/s', missing_value=missing_value )
!--------------------------
! 200-mb winds:
!--------------------------
       id_u200 = register_diag_field ( trim(field), 'u200', axes(1:2), Time,       &
                           '200-mb u-wind', '1/s', missing_value=missing_value )
       id_v200 = register_diag_field ( trim(field), 'v200', axes(1:2), Time,       &
                           '200-mb v-wind', '1/s', missing_value=missing_value )
       id_w200 = register_diag_field ( trim(field), 'w200', axes(1:2), Time,       &
                           '200-mb w-wind', '1/s', missing_value=missing_value )
! s200: wind speed for computing KE spectrum
! Cubed_2_latlon interpolation is more accurate, particularly near the poles, using
! winds speed (a scalar), rather than wind vectors or kinetic energy directly.
       id_s200 = register_diag_field ( trim(field), 's200', axes(1:2), Time,       &
                           '200-mb wind_speed', 'm/s', missing_value=missing_value )
       id_sl12 = register_diag_field ( trim(field), 'sl12', axes(1:2), Time,       &
                           '12th L wind_speed', 'm/s', missing_value=missing_value )
       id_sl13 = register_diag_field ( trim(field), 'sl13', axes(1:2), Time,       &
                           '13th L wind_speed', 'm/s', missing_value=missing_value )
!--------------------------
! 500-mb winds:
!--------------------------
       id_u500 = register_diag_field ( trim(field), 'u500', axes(1:2), Time,       &
                           '500-mb u-wind', '1/s', missing_value=missing_value )
       id_v500 = register_diag_field ( trim(field), 'v500', axes(1:2), Time,       &
                           '500-mb v-wind', '1/s', missing_value=missing_value )
!--------------------------
! 700-mb winds:
!--------------------------
       id_u700 = register_diag_field ( trim(field), 'u700', axes(1:2), Time,       &
                           '700-mb u-wind', '1/s', missing_value=missing_value )
       id_v700 = register_diag_field ( trim(field), 'v700', axes(1:2), Time,       &
                           '700-mb v-wind', '1/s', missing_value=missing_value )
!--------------------------
! 850-mb winds:
!--------------------------
       id_u850 = register_diag_field ( trim(field), 'u850', axes(1:2), Time,       &
                           '850-mb u-wind', '1/s', missing_value=missing_value )
       id_v850 = register_diag_field ( trim(field), 'v850', axes(1:2), Time,       &
                           '850-mb v-wind', '1/s', missing_value=missing_value )
       id_w850 = register_diag_field ( trim(field), 'w850', axes(1:2), Time,       &
                           '850-mb w-wind', '1/s', missing_value=missing_value )
!--------------------------
! temperature:
!--------------------------
       id_t50 = register_diag_field ( trim(field), 't50', axes(1:2), Time,       &
                           '50-mb temperature', 'K', missing_value=missing_value )
       id_t200 = register_diag_field ( trim(field), 't200', axes(1:2), Time,       &
                           '200-mb temperature', 'K', missing_value=missing_value )
       id_t500 = register_diag_field ( trim(field), 't500', axes(1:2), Time,       &
                           '500-mb temperature', 'K', missing_value=missing_value )
       id_t700 = register_diag_field ( trim(field), 't700', axes(1:2), Time,       &
                           '700-mb temperature', 'K', missing_value=missing_value )
       id_t850 = register_diag_field ( trim(field), 't850', axes(1:2), Time,       &
                           '850-mb temperature', 'K', missing_value=missing_value )
!--------------------------
! specific humidity:
!--------------------------
       id_q50 = register_diag_field ( trim(field), 'q50', axes(1:2), Time,       &
                           '50-mb specific humidity', 'kg/kg', missing_value=missing_value )
       id_q200 = register_diag_field ( trim(field), 'q200', axes(1:2), Time,       &
                           '200-mb specific humidity', 'kg/kg', missing_value=missing_value )
       id_q500 = register_diag_field ( trim(field), 'q500', axes(1:2), Time,       &
                           '500-mb specific humidity', 'kg/kg', missing_value=missing_value )
       id_q850 = register_diag_field ( trim(field), 'q850', axes(1:2), Time,       &
                           '850-mb specific humidity', 'kg/kg', missing_value=missing_value )
!--------------------------
! specific humidity:
!--------------------------
       id_omg200 = register_diag_field ( trim(field), 'omg200', axes(1:2), Time,       &
                           '200-mb omega', 'Pa/s', missing_value=missing_value )
       id_omg500 = register_diag_field ( trim(field), 'omg500', axes(1:2), Time,       &
                           '500-mb omega', 'Pa/s', missing_value=missing_value )
       id_omg850 = register_diag_field ( trim(field), 'omg850', axes(1:2), Time,       &
                           '850-mb omega', 'Pa/s', missing_value=missing_value )

!--------------------
! Tracer diagnostics:
!--------------------
       do i=1, ncnst
           call get_tracer_names ( MODEL_ATMOS, i, tname, tlongname, tunits )
           id_tracer(i) = register_diag_field ( field, trim(tname),  &
                axes(1:3), Time, trim(tlongname), &
                trim(tunits), missing_value=missing_value)
           if (master) then
               if (id_tracer(i) > 0) then
                   unit = stdlog()
                   write(unit,'(a,a,a,a)') &
                        & 'Diagnostics available for tracer ',tname, &
                        ' in module ', field
               end if
           endif
       enddo

       if ( id_mq > 0 )  then
            allocate ( zxg(isc:iec,jsc:jec) )
! Initialize gradient of terrain for mountain torque computation:
            call init_mq(Atm(n)%phis, Atm(n)%agrid(isc:iec,jsc:jec,2), npx, npy, isc, iec, jsc, jec, Atm(n)%ng)
       endif

    end do

    call nullify_domain()  ! Nullify  set_domain info

    module_is_initialized=.true.
 end subroutine fv_diag_init


 subroutine init_mq(phis, rlat, npx, npy, is, ie, js, je, ng)
    integer, intent(in):: npx, npy, is, ie, js, je, ng
    real, intent(in):: phis(is-ng:ie+ng, js-ng:je+ng)
    real, intent(in):: rlat(is:ie, js:je)  ! latitude (radian)
! local:
    real zs(is-ng:ie+ng, js-ng:je+ng)
    real zb(is-ng:ie+ng, js-ng:je+ng)
    real pdx(3,is:ie,js:je+1)
    real pdy(3,is:ie+1,js:je)
    integer i, j, n

!   do j=js,je
!      do i=is,ie
    do j=js-ng,je+ng
       do i=is-ng,ie+ng
          zs(i,j) = phis(i,j) / grav
       enddo
    enddo
!   call mpp_update_domains( zs, domain )

    call a2b_ord4(zs, zb, npx, npy, is, ie, js, je, ng)

    do j=js,je+1
       do i=is,ie
          do n=1,3
             pdx(n,i,j) = 0.5*(zb(i,j)+zb(i+1,j))*dx(i,j)*en1(n,i,j)
          enddo
       enddo
    enddo
    do j=js,je
       do i=is,ie+1
          do n=1,3
             pdy(n,i,j) = 0.5*(zb(i,j)+zb(i,j+1))*dy(i,j)*en2(n,i,j)
          enddo
       enddo
    enddo

! Compute gradient by Green's theorem
    do j=js,je
       do i=is,ie
          zxg(i,j) = vlon(i,j,1)*(pdx(1,i,j+1)-pdx(1,i,j)-pdy(1,i,j)+pdy(1,i+1,j))  &
                   + vlon(i,j,2)*(pdx(2,i,j+1)-pdx(2,i,j)-pdy(2,i,j)+pdy(2,i+1,j))  &
                   + vlon(i,j,3)*(pdx(3,i,j+1)-pdx(3,i,j)-pdy(3,i,j)+pdy(3,i+1,j))
! Times surface pressure to get Hadleys per unit area
! Unit Hadley = 1.E18 kg m**2 / s**2
          zxg(i,j) = -zxg(i,j) * radius * cos(rlat(i,j)) * rarea(i,j) * 1.E-18
       enddo
    enddo

 end subroutine init_mq

 subroutine fv_diag(Atm, zvir, Time, print_freq)

    type(fv_atmos_type), intent(inout) :: Atm(:)
    type(time_type),     intent(in) :: Time
    real,                intent(in):: zvir
    integer,             intent(in):: print_freq

    integer :: isc, iec, jsc, jec, n, ntileMe
    integer :: isd, ied, jsd, jed, npz, itrac
    integer :: ngc, nwater

    real, allocatable :: a2(:,:),a3(:,:,:), wk(:,:,:), wz(:,:,:), ucoor(:,:,:), vcoor(:,:,:)
    real, allocatable :: slp(:,:), depress(:,:), ws_max(:,:), tc_count(:,:)
    real, allocatable :: u2(:,:), v2(:,:)
    real height(2)
    real plevs(6)
    real tot_mq, tmp
    logical :: used
    logical :: bad_range
    logical :: prt_minmax
    integer i,j,k, yr, mon, dd, hr, mn, days, seconds
    character(len=128)   :: tname
    real, parameter:: ws_0 = 16.   ! minimum max_wind_speed within the 7x7 search box
    real, parameter:: ws_1 = 20.
    real, parameter:: vort_c0= 2.2e-5 
    logical, allocatable :: storm(:,:), cat_crt(:,:)

#ifdef MARS_GCM
    real  ::   atm_mass,  sfc_mass, atm_cloud
    real  ::   tsfc_dust, tcol_dust
#endif

! cat15: SLP<1000; srf_wnd>ws_0; vort>vort_c0
! cat25: SLP< 980; srf_wnd>ws_1; vort>vort_c0
! cat35: SLP< 964; srf_wnd>ws_1; vort>vort_c0
! cat45: SLP< 944; srf_wnd>ws_1; vort>vort_c0

    height(1) = 5.E3      ! for computing 5-km "pressure"
    height(2) = 0.        ! for sea-level pressure

    ntileMe = size(Atm(:))
    n = 1
    isc = Atm(n)%isc; iec = Atm(n)%iec
    jsc = Atm(n)%jsc; jec = Atm(n)%jec
    ngc = Atm(n)%ng
    npz = Atm(n)%npz
    ptop = Atm(n)%ak(1)

    isd = Atm(n)%isd; ied = Atm(n)%ied
    jsd = Atm(n)%jsd; jed = Atm(n)%jed


    if( id_c15>0 ) then
        allocate (   storm(isc:iec,jsc:jec) )
        allocate ( depress(isc:iec,jsc:jec) )
        allocate (  ws_max(isc:iec,jsc:jec) )
        allocate ( cat_crt(isc:iec,jsc:jec) )
        allocate (tc_count(isc:iec,jsc:jec) )
    endif

    fv_time = Time
    call set_domain(Atm(1)%domain)

    if ( moist_phys ) then
#if defined(MARS_GCM) || defined(VENUS_GCM)
         call get_time (fv_time, seconds,  days)
         mn= 0
         hr= 0
         mon= 0
#else
         call get_date(fv_time, yr, mon, dd, hr, mn, seconds)
#endif 
         if( print_freq == 0 ) then
                 prt_minmax = .false.
         elseif( print_freq < 0 ) then
                 prt_minmax = .true.
         else
                 prt_minmax = mod(hr, print_freq) == 0 .and. mn==0 .and. seconds==0
         endif
     else
         call get_time (fv_time, seconds,  days)
         if( print_freq == 0 ) then
                 prt_minmax = .false.
         elseif( print_freq < 0 ) then
                 prt_minmax = .true.
         else
                 prt_minmax = mod(seconds, 3600*print_freq) == 0
         endif
     endif

     if(prt_minmax) then
#if defined(MARS_GCM) || defined(VENUS_GCM)
        if(master) write(6,*) Days, seconds
#else
         if ( moist_phys ) then
              if(master) write(6,*) yr, mon, dd, hr, mn, seconds
         else
              if(master) write(6,*) Days, seconds
         endif
#endif
     endif

    if( prt_minmax ) then

        call prt_maxmin('ZS', zsurf,     isc, iec, jsc, jec, 0,   1, 1.0,  master)
        call prt_maxmin('PS', Atm(n)%ps, isc, iec, jsc, jec, ngc, 1, 0.01, master)

        call prt_mass(npz, ncnst, isc, iec, jsc, jec, ngc, Atm(n)%nwat,    &
                      Atm(n)%ps, Atm(n)%delp, Atm(n)%q, master)
#ifndef SW_DYNAMICS
             steps = steps + 1
           efx_sum = efx_sum + E_Flux
        if ( steps <= max_step ) efx(steps) = E_Flux
        if (master)  then
            write(6,*) 'ENG Deficit (W/m**2)=', E_Flux
        endif
#endif
        call prt_maxmin('UA', Atm(n)%ua, isc, iec, jsc, jec, ngc, npz, 1., master)
        call prt_maxmin('VA', Atm(n)%va, isc, iec, jsc, jec, ngc, npz, 1., master)

        if ( .not. Atm(n)%hydrostatic ) then
          call prt_maxmin('W ', Atm(n)%w , isc, iec, jsc, jec, ngc, npz, 1., master)
          if ( Atm(n)%hybrid_z ) call prt_maxmin('Hybrid_ZTOP (km)', Atm(n)%ze0(isc:iec,jsc:jec,1), &
                                                 isc, iec, jsc, jec, 0, 1, 1.E-3, master)
          call prt_maxmin('Bottom DZ (m)', Atm(n)%delz(isc:iec,jsc:jec,npz),    &
                          isc, iec, jsc, jec, 0, 1, 1., master)
        endif

#ifndef SW_DYNAMICS
        call prt_maxmin('TA', Atm(n)%pt,   isc, iec, jsc, jec, ngc, npz, 1., master)
        call prt_maxmin('OM', Atm(n)%omga, isc, iec, jsc, jec, ngc, npz, 1., master)
#endif

#if defined(MARS_GCM) && defined(MARS_SURFACE)
        atm_mass  = g_sum( Atm(n)%ps(isc:iec,jsc:jec), isc, iec, jsc, jec, ngc, area,mode=1)
        sfc_mass  = g_sum( sfc_snow,isc, iec, jsc, jec, ngc, area,mode=1)
        sfc_mass= sfc_mass*grav   !   Conversion to pressure units

        if(master) write(*,*) 'Atmospheric CO2 (mb) =', atm_mass*0.01
        if(master) write(*,*) 'CO2 sfc frost   (mb) =', sfc_mass*0.01
        if(master) write(*,*) 'Total CO2 Inventory  =', (atm_mass+sfc_mass)*0.01

#ifdef WATER_CYCLE
        sfc_mass  = g_sum( sfc_frost, isc, iec, jsc, jec, ngc, area,mode=1)
        atm_mass  = g_sum(      wcol, isc, iec, jsc, jec, ngc, area,mode=1)
        atm_cloud = g_sum(    cldcol, isc, iec, jsc, jec, ngc, area,mode=1)
        sfc_mass= sfc_mass - 3.7   !  Arbitrary offset

        if(master) write(*,*) 'Atmospheric H2o vapor (kg/m**2) =', atm_mass
        if(master) write(*,*) 'Atmospheric H2o cloud (kg/m**2) =', atm_cloud
        if(master) write(*,*) 'Total Atmospheric H2o ', atm_cloud + atm_mass

        if(master) write(*,*) 'H2O surface frost (kg/m**2) ==', sfc_mass
        if(master) write(*,*) 'Total H2O inventory =', atm_mass+sfc_mass+atm_cloud
#endif WATER_CYCLE

#ifdef DUST_SOURCE
        tsfc_dust  = g_sum( sfc_dust(:,:,1),isc, iec, jsc, jec, ngc, area,mode=1)
        tcol_dust  = g_sum( odcol   (:,:,1),isc, iec, jsc, jec, ngc, area,mode=1)

        if(master) write(*,*) 'Surface dust inventory (kg/m**2) =', tsfc_dust - 30.0
        if(master) write(*,*) 'Atmospheric dust (kg/m**2) =', tcol_dust
        if(master) write(*,*) 'Total dust inventory ', tsfc_dust - 30.0 + tcol_dust
#endif DUST_SOURCE
#endif

    elseif ( Atm(n)%range_warn ) then
         call range_check('DELP', Atm(n)%delp, isc, iec, jsc, jec, ngc, npz, Atm(n)%agrid,    &
                           master, 0.1*ptop, 200.E2, bad_range)
         call range_check('UA', Atm(n)%ua, isc, iec, jsc, jec, ngc, npz, Atm(n)%agrid,   &
                           master, -220., 250., bad_range)
         call range_check('VA', Atm(n)%ua, isc, iec, jsc, jec, ngc, npz, Atm(n)%agrid,   &
                           master, -220., 220., bad_range)
#ifndef SW_DYNAMICS
         call range_check('TA', Atm(n)%pt, isc, iec, jsc, jec, ngc, npz, Atm(n)%agrid,   &
                           master, 150., 350., bad_range)
#endif

    endif

    allocate ( a2(isc:iec,jsc:jec) )
    allocate ( wk(isc:iec,jsc:jec,npz) )

    do n = 1, ntileMe

#ifdef DYNAMICS_ZS
       if(id_zsurf > 0)  used=send_data(id_zsurf, zsurf, Time)
#endif
       if(id_ps > 0) used=send_data(id_ps, Atm(n)%ps(isc:iec,jsc:jec), Time)

       if(id_c15>0 .or. id_c25>0 .or. id_c35>0 .or. id_c45>0) then
          call wind_max(isc, iec, jsc, jec ,isd, ied, jsd, jed, Atm(n)%ua(isc:iec,jsc:jec,npz),   &
                        Atm(n)%va(isc:iec,jsc:jec,npz), ws_max)
          do j=jsc,jec
             do i=isc,iec
                if( abs(Atm(n)%agrid(i,j,2)*rad2deg)<60.0 .and.     &
                    Atm(n)%phis(i,j)*ginv<500.0 .and. ws_max(i,j)>ws_0 ) then
                    storm(i,j) = .true.
                else
                    storm(i,j) = .false.
                endif
             enddo
          enddo
       endif

       if ( id_vort850>0 .or. id_vorts>0 .or. id_vort>0 .or. id_pv>0 .or. id_rh>0 ) then
          call get_vorticity(isc, iec, jsc, jec, isd, ied, jsd, jed, npz, Atm(n)%u, Atm(n)%v, wk)

          if(id_vort >0) used=send_data(id_vort,  wk, Time)
          if(id_vorts>0) used=send_data(id_vorts, wk(isc:iec,jsc:jec,npz), Time)

          if(id_c15>0) then
             do j=jsc,jec
                do i=isc,iec
                   if ( storm(i,j) )    &
                   storm(i,j) = (Atm(n)%agrid(i,j,2)>0. .and. wk(i,j,npz)> vort_c0) .or. &
                                (Atm(n)%agrid(i,j,2)<0. .and. wk(i,j,npz)<-vort_c0) 
                enddo
             enddo
          endif


          if(id_vort850>0) then
             call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                       850.e2, Atm(n)%peln, wk, a2)
             used=send_data(id_vort850, a2, Time)

             if(id_c15>0) then
             do j=jsc,jec
                do i=isc,iec
                   if ( storm(i,j) )    &
                     storm(i,j) = (Atm(n)%agrid(i,j,2)>0. .and. a2(i,j)> vort_c0) .or.     &
                                  (Atm(n)%agrid(i,j,2)<0. .and. a2(i,j)<-vort_c0) 
                enddo
             enddo
             endif

          endif

          if ( id_pv > 0 ) then
! Note: this is expensive computation.
              call pv_entropy(isc, iec, jsc, jec, ngc, npz, wk,    &
                              f0, Atm(n)%pt, Atm(n)%pkz, Atm(n)%delp, grav)
              used = send_data ( id_pv, wk, Time )
          endif

! Relative Humidity
          if ( id_rh > 0 ) then
! Compute FV mean pressure
               do k=1,npz
                  do j=jsc,jec
                     do i=isc,iec
                        a2(i,j) = Atm(n)%delp(i,j,k)/(Atm(n)%peln(i,k+1,j)-Atm(n)%peln(i,k,j))
                     enddo
                  enddo
                  call qsmith(iec-isc+1, jec-jsc+1, 1, Atm(n)%pt(isc:iec,jsc:jec,k),   &
                              a2, Atm(n)%q(isc:iec,jsc:jec,k,sphum), wk(isc,jsc,k))
                  do j=jsc,jec
                     do i=isc,iec
                        wk(i,j,k) = 100.*Atm(n)%q(i,j,k,sphum)/wk(i,j,k)
                     enddo
                  enddo
               enddo
               used = send_data ( id_rh, wk, Time )
               if(prt_minmax) then
                  call prt_maxmin('RH_sf (%)', wk(isc:iec,jsc:jec,npz), isc, iec, jsc, jec, 0,   1, 1., master)
                  call prt_maxmin('RH_3D (%)', wk, isc, iec, jsc, jec, 0, npz, 1., master)
               endif
          endif

       endif

       if(id_c25>0 .or. id_c35>0 .or. id_c45>0) then
          do j=jsc,jec
             do i=isc,iec
                if ( storm(i,j) .and. ws_max(i,j)>ws_1 ) then
                     cat_crt(i,j) = .true.
                else
                     cat_crt(i,j) = .false.
                endif
             enddo
          enddo
       endif



       if( id_slp>0 .or. id_tm>0 .or. id_h50>0 .or. id_h200>0 .or. id_h300>0 .or. id_h500>0 .or. &
           id_h700>0 .or. id_h850>0 .or. id_c15>0 ) then

          allocate ( wz(isc:iec,jsc:jec,npz+1) )
          call get_height_field(isc, iec, jsc, jec, ngc, npz, wz, Atm(n)%pt, Atm(n)%q, Atm(n)%peln, zvir)
          if( prt_minmax )   &
          call prt_maxmin('ZTOP', wz(isc:iec,jsc:jec,1), isc, iec, jsc, jec, 0, 1, 1.E-3, master)

          if(id_slp > 0) then
! Cumpute SLP (pressure at height=0)
          allocate ( slp(isc:iec,jsc:jec) )
          call get_pressure_given_height(isc, iec, jsc, jec, ngc, npz, wz, 1, height(2),   &
                                        Atm(n)%pt(:,:,npz), Atm(n)%peln, slp, 0.01)
          used = send_data (id_slp, slp, Time)
          if( prt_minmax )   &
             call prt_maxmin('SLP', slp, isc, iec, jsc, jec, 0, 1, 1., master)
          endif

! Compute H3000 and/or H500
          if( id_tm>0 .or. id_h50>0 .or. id_h200>0 .or. id_h300>0 .or. id_h500>0 .or. id_h700>0 .or. id_h850>0 .or. id_ppt>0) then

              allocate( a3(isc:iec,jsc:jec,6) )
              plevs(1) = log( 5000. )
              plevs(2) = log( 20000. )
              plevs(3) = log( 30000. )
              plevs(4) = log( 50000. )
              plevs(5) = log( 70000. )
              plevs(6) = log( 85000. )

             call get_height_given_pressure(isc, iec, jsc, jec, ngc, npz, wz, 6, plevs, Atm(n)%peln, a3)
             if(id_h50>0)  used = send_data ( id_h50,  a3(isc:iec,jsc:jec,1), Time )
             if(id_h200>0) used = send_data ( id_h200, a3(isc:iec,jsc:jec,2), Time )
             if(id_h300>0) used = send_data ( id_h300, a3(isc:iec,jsc:jec,3), Time )
             if(id_h500>0) used = send_data ( id_h500, a3(isc:iec,jsc:jec,4), Time )
             if(id_h700>0) used = send_data ( id_h700, a3(isc:iec,jsc:jec,5), Time )
             if(id_h850>0) used = send_data ( id_h850, a3(isc:iec,jsc:jec,6), Time )

             if( id_tm>0 ) then
                 do j=jsc,jec
                    do i=isc,iec
                       a2(i,j) = grav*(a3(i,j,3)-a3(i,j,4))/(rdgas*(plevs(4)-plevs(3)))
                    enddo
                 enddo
                 used = send_data ( id_tm, a2, Time )
             endif

            if(id_c15>0 .or. id_c25>0 .or. id_c35>0 .or. id_c45>0) then
             do j=jsc,jec
                do i=isc,iec
! Minimum warm core:
                   if ( storm(i,j) ) then
                        if( a2(i,j)<254.0 .or. Atm(n)%pt(i,j,npz)<281.0 ) Then
                              storm(i,j) = .false.
                            cat_crt(i,j) = .false.
                        endif
                   endif
                enddo
             enddo
! Cat 1-5:
             do j=jsc,jec
                do i=isc,iec
                   if ( storm(i,j) .and. slp(i,j)<1000.0 ) then
                         depress(i,j) = 1000. - slp(i,j)
                        tc_count(i,j) = 1.
                   else
                         depress(i,j) = 0.
                        tc_count(i,j) = 0.
                   endif
                enddo
             enddo
             used = send_data(id_c15, depress, Time)
             if(id_f15>0) used = send_data(id_f15, tc_count, Time)
             if(prt_minmax) then
!               tmp = g_sum(depress, isc, iec, jsc, jec, ngc, area, 1) 
!               if(master) write(*,*) 'Mean Tropical Cyclone depression (mb)=', tmp
                call prt_maxmin('Depression', depress, isc, iec, jsc, jec, 0,   1, 1., master)
             endif
            endif

! Cat 2-5:
            if(id_c25>0) then
             do j=jsc,jec
                do i=isc,iec
                   if ( cat_crt(i,j) .and. slp(i,j)<980.0 ) then
                        depress(i,j) = 980. - slp(i,j)
                       tc_count(i,j) = 1.
                   else
                        depress(i,j) = 0.
                       tc_count(i,j) = 0.
                   endif
                enddo
             enddo
             used = send_data(id_c25, depress, Time)
             if(id_f25>0) used = send_data(id_f25, tc_count, Time)
            endif

! Cat 3-5:
            if(id_c35>0) then
             do j=jsc,jec
                do i=isc,iec
                   if ( cat_crt(i,j) .and. slp(i,j)<964.0 ) then
                        depress(i,j) = 964. - slp(i,j)
                       tc_count(i,j) = 1.
                   else
                        depress(i,j) = 0.
                       tc_count(i,j) = 0.
                   endif
                enddo
             enddo
             used = send_data(id_c35, depress, Time)
             if(id_f35>0) used = send_data(id_f35, tc_count, Time)
            endif

! Cat 4-5:
            if(id_c45>0) then
             do j=jsc,jec
                do i=isc,iec
                   if ( cat_crt(i,j) .and. slp(i,j)<944.0 ) then
                        depress(i,j) = 944. - slp(i,j)
                       tc_count(i,j) = 1.
                   else
                        depress(i,j) = 0.
                       tc_count(i,j) = 0.
                   endif
                enddo
             enddo
             used = send_data(id_c45, depress, Time)
             if(id_f45>0) used = send_data(id_f45, tc_count, Time)
            endif

            if (id_c15>0) then
                deallocate(depress)
                deallocate(cat_crt)
                deallocate(storm)
                deallocate(ws_max)
                deallocate(tc_count)
            endif

            if(id_slp>0 )  deallocate( slp )

            deallocate( a3 )
          endif

         deallocate ( wz )
       endif


       if(id_mq > 0)  then
          do j=jsc,jec
             do i=isc,iec
                a2(i,j) = Atm(n)%ps(i,j)*zxg(i,j)
             enddo
          enddo
          used = send_data(id_mq, a2, Time)
          if( prt_minmax ) then
              tot_mq  = g_sum( a2, isc, iec, jsc, jec, ngc, area, 0) 
              mtq_sum = mtq_sum + tot_mq
              if ( steps <= max_step ) mtq(steps) = tot_mq
              if(master) write(*,*) 'Total (global) mountain torque (Hadleys)=', tot_mq
          endif
       endif

#ifdef MARS_GCM
       if ( id_t05>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      0.5e2, Atm(n)%peln, Atm(n)%pt(isc:iec,jsc:jec,:), a2)
            used=send_data(id_t05, a2, Time)
       endif
#  ifdef WATER_CYCLE
       if ( id_tq>0 ) then
          itrac= get_tracer_index (MODEL_ATMOS, 'h2o_vapor')
          a2 = 0.
          do k=1,npz
          do j=jsc,jec
             do i=isc,iec
                a2(i,j) = a2(i,j) + Atm(n)%q(i,j,k,itrac)*Atm(n)%delp(i,j,k)
             enddo
          enddo
          enddo
          used = send_data(id_tq, a2*ginv, Time)
       endif
#  endif WATER_CYCLE
#else
       if ( id_tq>0 ) then
          nwater = Atm(1)%nwat
          a2 = 0.
          do k=1,npz
          do j=jsc,jec
             do i=isc,iec
#ifdef OLD_TQ
                a2(i,j) = a2(i,j) + Atm(n)%q(i,j,k,1)*Atm(n)%delp(i,j,k)
#else
                a2(i,j) = a2(i,j) + sum(Atm(n)%q(i,j,k,1:nwater))*Atm(n)%delp(i,j,k)
#endif
             enddo
          enddo
          enddo
          used = send_data(id_tq, a2*ginv, Time)
       endif
#endif MARS_GCM

       if(id_us > 0) used=send_data(id_us, Atm(n)%ua(isc:iec,jsc:jec,npz), Time)
       if(id_vs > 0) used=send_data(id_vs, Atm(n)%va(isc:iec,jsc:jec,npz), Time)

       if(id_ua > 0) used=send_data(id_ua, Atm(n)%ua(isc:iec,jsc:jec,:), Time)
       if(id_va > 0) used=send_data(id_va, Atm(n)%va(isc:iec,jsc:jec,:), Time)

! 50-mb
       if ( id_u50>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      50.e2, Atm(n)%peln, Atm(n)%ua(isc:iec,jsc:jec,:), a2)
            used=send_data(id_u50, a2, Time)
       endif
       if ( id_v50>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      50.e2, Atm(n)%peln, Atm(n)%va(isc:iec,jsc:jec,:), a2)
            used=send_data(id_v50, a2, Time)
       endif
       if ( id_t50>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      50.e2, Atm(n)%peln, Atm(n)%pt(isc:iec,jsc:jec,:), a2)
            used=send_data(id_t50, a2, Time)
       endif
       if ( id_q50>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      50.e2, Atm(n)%peln, Atm(n)%q(isc:iec,jsc:jec,:,sphum), a2)
            used=send_data(id_q50, a2, Time)
       endif
! 200-mb
       if ( id_u200>0 .or. id_s200>0 ) then
            allocate( u2(isc:iec,jsc:jec) )
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      200.e2, Atm(n)%peln, Atm(n)%ua(isc:iec,jsc:jec,:), u2)
            if( id_u200>0 ) used=send_data(id_u200, u2, Time)
       endif
       if ( id_v200>0 .or. id_s200>0 ) then
            allocate( v2(isc:iec,jsc:jec) )
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      200.e2, Atm(n)%peln, Atm(n)%va(isc:iec,jsc:jec,:), v2)
            if( id_v200>0 ) used=send_data(id_v200, v2, Time)
       endif
       if ( id_s200>0 ) then
            do j=jsc,jec
               do i=isc,iec
                  a2(i,j) = sqrt(u2(i,j)**2 + v2(i,j)**2)
               enddo
            enddo
            used=send_data(id_s200, a2, Time)
       endif
       if ( id_sl12>0 ) then   ! 13th level wind speed (~ 222 mb for the 32L setup)
            do j=jsc,jec
               do i=isc,iec
                  a2(i,j) = sqrt(Atm(n)%ua(i,j,12)**2 + Atm(n)%va(i,j,12)**2)
               enddo
            enddo
            used=send_data(id_sl12, a2, Time)
       endif
       if ( id_sl13>0 ) then   ! 13th level wind speed (~ 222 mb for the 32L setup)
            do j=jsc,jec
               do i=isc,iec
                  a2(i,j) = sqrt(Atm(n)%ua(i,j,13)**2 + Atm(n)%va(i,j,13)**2)
               enddo
            enddo
            used=send_data(id_sl13, a2, Time)
       endif
       if ( allocated (u2) )  deallocate ( u2 )
       if ( allocated (v2) )  deallocate ( v2 )

       if ( id_w200>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      200.e2, Atm(n)%peln, Atm(n)%w(isc:iec,jsc:jec,:), a2)
            used=send_data(id_w200, a2, Time)
       endif
       if ( id_t200>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      200.e2, Atm(n)%peln, Atm(n)%pt(isc:iec,jsc:jec,:), a2)
            used=send_data(id_t200, a2, Time)
       endif
       if ( id_q200>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      200.e2, Atm(n)%peln, Atm(n)%q(isc:iec,jsc:jec,:,sphum), a2)
            used=send_data(id_q200, a2, Time)
       endif
       if ( id_omg200>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      200.e2, Atm(n)%peln, Atm(n)%omga(isc:iec,jsc:jec,:), a2)
            used=send_data(id_omg200, a2, Time)
       endif
! 500-mb
       if ( id_u500>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      500.e2, Atm(n)%peln, Atm(n)%ua(isc:iec,jsc:jec,:), a2)
            used=send_data(id_u500, a2, Time)
       endif
       if ( id_v500>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      500.e2, Atm(n)%peln, Atm(n)%va(isc:iec,jsc:jec,:), a2)
            used=send_data(id_v500, a2, Time)
       endif
       if ( id_t500>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      500.e2, Atm(n)%peln, Atm(n)%pt(isc:iec,jsc:jec,:), a2)
            used=send_data(id_t500, a2, Time)
       endif
       if ( id_q500>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      500.e2, Atm(n)%peln, Atm(n)%q(isc:iec,jsc:jec,:,sphum), a2)
            used=send_data(id_q500, a2, Time)
       endif
       if ( id_omg500>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      500.e2, Atm(n)%peln, Atm(n)%omga(isc:iec,jsc:jec,:), a2)
            used=send_data(id_omg500, a2, Time)
       endif
! 700-mb
       if ( id_u700>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      700.e2, Atm(n)%peln, Atm(n)%ua(isc:iec,jsc:jec,:), a2)
            used=send_data(id_u700, a2, Time)
       endif
       if ( id_v700>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      700.e2, Atm(n)%peln, Atm(n)%va(isc:iec,jsc:jec,:), a2)
            used=send_data(id_v700, a2, Time)
       endif
       if ( id_t700>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      700.e2, Atm(n)%peln, Atm(n)%pt(isc:iec,jsc:jec,:), a2)
            used=send_data(id_t700, a2, Time)
       endif
! 850-mb
       if ( id_u850>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      850.e2, Atm(n)%peln, Atm(n)%ua(isc:iec,jsc:jec,:), a2)
            used=send_data(id_u850, a2, Time)
       endif
       if ( id_v850>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      850.e2, Atm(n)%peln, Atm(n)%va(isc:iec,jsc:jec,:), a2)
            used=send_data(id_v850, a2, Time)
       endif
       if ( id_w850>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      850.e2, Atm(n)%peln, Atm(n)%w(isc:iec,jsc:jec,:), a2)
            used=send_data(id_w850, a2, Time)
       endif
       if ( id_t850>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      850.e2, Atm(n)%peln, Atm(n)%pt(isc:iec,jsc:jec,:), a2)
            used=send_data(id_t850, a2, Time)
       endif
       if ( id_q850>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      850.e2, Atm(n)%peln, Atm(n)%q(isc:iec,jsc:jec,:,sphum), a2)
            used=send_data(id_q850, a2, Time)
       endif
       if ( id_omg850>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      850.e2, Atm(n)%peln, Atm(n)%omga(isc:iec,jsc:jec,:), a2)
            used=send_data(id_omg850, a2, Time)
       endif

       if ( .not.Atm(n)%hydrostatic .and. id_w>0  )     &
                 used=send_data(id_w, Atm(n)%w(isc:iec,jsc:jec,:), Time)

       if(id_pt   > 0) used=send_data(id_pt  , Atm(n)%pt  (isc:iec,jsc:jec,:), Time)
       if(id_omga > 0) used=send_data(id_omga, Atm(n)%omga(isc:iec,jsc:jec,:), Time)

       if(id_ppt> 0) then
! Potential temperature perturbation for gravity wave test_case
          allocate ( pt1(npz) )
#ifdef TEST_GWAVES
          call gw_1d(npz, 1000.E2, Atm(n)%ak, Atm(n)%ak, Atm(n)%ak(1), 10.E3, pt1)
#else
          pt1 = 0. 
#endif
          do k=1,npz
          do j=jsc,jec
             do i=isc,iec
                wk(i,j,k) =  (Atm(n)%pt(i,j,k)/Atm(n)%pkz(i,j,k) - pt1(k)) * pk0
             enddo
          enddo
          enddo
          used=send_data(id_ppt, wk, Time)

          if( prt_minmax ) then
              call prt_maxmin('PoTemp', wk, isc, iec, jsc, jec, 0, npz, 1., master)
          endif

          deallocate ( pt1 )
       endif


        do itrac=1, ncnst
          if (id_tracer(itrac) > 0) &
               & used = send_data (id_tracer(itrac), Atm(n)%q(isc:iec,jsc:jec,:,itrac), Time )
          if( prt_minmax ) then
              call get_tracer_names ( MODEL_ATMOS, itrac, tname )
#ifndef SW_DYNAMICS
              call prt_maxmin(trim(tname), Atm(n)%q(:,:,1,itrac), &
                              isc, iec, jsc, jec, ngc, npz, 1., master)
#endif
          endif
        enddo

    enddo

    deallocate ( a2 )
    deallocate ( wk )

    call nullify_domain()


 end subroutine fv_diag

 subroutine wind_max(isc, iec, jsc, jec ,isd, ied, jsd, jed, us, vs, ws_max)
 integer isc, iec, jsc, jec
 integer isd, ied, jsd, jed
 real, intent(in), dimension(isc:iec,jsc:jec):: us, vs
 real, intent(out) :: ws_max(isc:iec,jsc:jec)
! Local
 real :: wx(isc:iec,jsd:jed), ws(isd:ied,jsd:jed)
 integer:: i,j

 ws = 0.   ! fill corners with zeros
 do j=jsc,jec
    do i=isc,iec
       ws(i,j) = sqrt(us(i,j)**2 + vs(i,j)**2)
    enddo
 enddo

 call mpp_update_domains( ws, domain )

 do j=jsd,jed
    do i=isc,iec
       wx(i,j) = max(ws(i-3,j), ws(i-2,j), ws(i-1,j), ws(i,j), ws(i+1,j), ws(i+2,j), ws(i+3,j))
    enddo
 enddo

 do j=jsc,jec
    do i=isc,iec
       ws_max(i,j) = max(wx(i,j-3), wx(i,j-2), wx(i,j-1), wx(i,j), wx(i,j+1), wx(i,j+2), wx(i,j+3))
    enddo
 enddo

 end subroutine wind_max


 subroutine get_vorticity(isc, iec, jsc, jec ,isd, ied, jsd, jed, npz, u, v, vort)
 integer isd, ied, jsd, jed, npz
 integer isc, iec, jsc, jec
 real, intent(in)  :: u(isd:ied, jsd:jed+1, npz), v(isd:ied+1, jsd:jed, npz)
 real, intent(out) :: vort(isc:iec, jsc:jec, npz)
! Local
 real :: utmp(isc:iec, jsc:jec+1), vtmp(isc:iec+1, jsc:jec)
 integer :: i,j,k

      do k=1,npz
         do j=jsc,jec+1
            do i=isc,iec
               utmp(i,j) = u(i,j,k)*dx(i,j)
            enddo
         enddo
         do j=jsc,jec
            do i=isc,iec+1
               vtmp(i,j) = v(i,j,k)*dy(i,j)
            enddo
         enddo

         do j=jsc,jec
            do i=isc,iec
               vort(i,j,k) = rarea(i,j)*(utmp(i,j)-utmp(i,j+1)-vtmp(i,j)+vtmp(i+1,j))
            enddo
         enddo
      enddo

 end subroutine get_vorticity


 subroutine get_height_field(is, ie, js, je, ng, km, wz, pt, q, peln, zvir)
  integer, intent(in):: is, ie, js, je, km, ng
  real, intent(in):: peln(is:ie,km+1,js:je)
  real, intent(in):: pt(is-ng:ie+ng,js-ng:je+ng,km)
  real, intent(in)::  q(is-ng:ie+ng,js-ng:je+ng,km,*) ! water vapor
  real, intent(in):: zvir
  real, intent(out):: wz(is:ie,js:je,km+1)
!
  integer i,j,k
  real gg

      gg  = rdgas * ginv

      do j=js,je
         do i=is,ie
            wz(i,j,km+1) = zsurf(i,j)
         enddo
         do k=km,1,-1
            do i=is,ie
               wz(i,j,k) = wz(i,j,k+1) + gg*pt(i,j,k)*(1.+zvir*q(i,j,k,sphum))  &
                          *(peln(i,k+1,j)-peln(i,k,j))
            enddo
         enddo
      enddo

 end subroutine get_height_field

 subroutine range_check(qname, q, is, ie, js, je, n_g, km, pos, master, q_low, q_hi, bad_range)
      character(len=*), intent(in)::  qname
      integer, intent(in):: is, ie, js, je
      integer, intent(in):: n_g, km
      real, intent(in)::    q(is-n_g:ie+n_g, js-n_g:je+n_g, km)
      real, intent(in):: pos(is-n_g:ie+n_g, js-n_g:je+n_g,2)
      real, intent(in):: q_low, q_hi
      logical, intent(in):: master
      logical, optional, intent(out):: bad_range
!
      real qmin, qmax
      integer i,j,k

      if ( present(bad_range) ) bad_range = .false. 
      qmin = q(is,js,1)
      qmax = qmin

      do k=1,km
      do j=js,je
         do i=is,ie
            if( q(i,j,k) < qmin ) then
                qmin = q(i,j,k)
            elseif( q(i,j,k) > qmax ) then
                qmax = q(i,j,k)
            endif
          enddo
      enddo
      enddo

      call mp_reduce_min(qmin)
      call mp_reduce_max(qmax)

      if( qmin<q_low .or. qmax>q_hi ) then
          if(master) write(6,*) 'Range_check Warning:', qname, ' max = ', qmax, ' min = ', qmin
          if ( present(bad_range) ) then
               bad_range = .true. 
          endif
      endif

      if ( present(bad_range) ) then
! Print out where the bad value(s) is (are)
         if ( bad_range .EQV. .false. ) return      ! kerr
         do k=1,km
            do j=js,je
               do i=is,ie
                  if( q(i,j,k)<q_low .or. q(i,j,k)>q_hi ) then
!                     write(*,*) 'gid=', gid, k,i,j, pos(i,j,1)*rad2deg, pos(i,j,2)*rad2deg, q(i,j,k)
                      write(*,*) 'k=',k,' (i,j)=',i,j, pos(i,j,1)*rad2deg, pos(i,j,2)*rad2deg, q(i,j,k)
                      if ( k/= 1 ) write(*,*) k-1, q(i,j,k-1)
                      if ( k/=km ) write(*,*) k+1, q(i,j,k+1)
                  endif
               enddo
            enddo
         enddo
         call mpp_error(FATAL,'==> Error from range_check: data out of bound')
      endif

 end subroutine range_check

 subroutine prt_maxmin(qname, q, is, ie, js, je, n_g, km, fac, master)
      character(len=*), intent(in)::  qname
      integer, intent(in):: is, ie, js, je
      integer, intent(in):: n_g, km
      real, intent(in)::    q(is-n_g:ie+n_g, js-n_g:je+n_g, km)
      real, intent(in)::    fac
      logical, intent(in):: master
!
      real qmin, qmax
      integer i,j,k

      qmin = q(is,js,1)
      qmax = qmin

      do k=1,km
      do j=js,je
         do i=is,ie
!           qmin = min(qmin, q(i,j,k))
!           qmax = max(qmax, q(i,j,k))
            if( q(i,j,k) < qmin ) then
                qmin = q(i,j,k)
            elseif( q(i,j,k) > qmax ) then
                qmax = q(i,j,k)
            endif
          enddo
      enddo
      enddo

      call mp_reduce_min(qmin)
      call mp_reduce_max(qmax)

      if(master) write(6,*) qname, ' max = ', qmax*fac, ' min = ', qmin*fac

 end subroutine prt_maxmin


 subroutine prt_mass(km, nq, is, ie, js, je, n_g, nwat, ps, delp, q, master)

 integer, intent(in):: is, ie, js, je
 integer, intent(in):: nq, n_g, km, nwat
 real, intent(in)::   ps(is-n_g:ie+n_g, js-n_g:je+n_g)
 real, intent(in):: delp(is-n_g:ie+n_g, js-n_g:je+n_g, km)
 real, intent(in)::  q(is-n_g:ie+n_g, js-n_g:je+n_g, km, nq)
 logical, intent(in):: master
! Local:
 real psq(is:ie,js:je,nwat)
 real q_strat(is:ie,js:je)
 real qtot(nwat)
 real psmo, totw, psdry
 integer k, n, kstrat

#if defined(MARS_GCM) || defined(VENUS_GCM)
 psmo = g_sum(ps(is:ie,js:je), is, ie, js, je, n_g, area, 1)
 totw  = 0.0
 psdry = psmo - totw
 if ( nwat > 0 ) then
    qtot(:)= 0.0
 endif

 if( master ) then
     write(6,*) 'Total surface pressure (mb) = ',  0.01*psmo
!!!     write(6,*) 'mean dry surface pressure = ',    0.01*psdry
  endif
#else

 if ( nwat==0 ) then
      psmo = g_sum(ps(is:ie,js:je), is, ie, js, je, n_g, area, 1) 
      if( master ) write(6,*) 'Total surface pressure (mb) = ',  0.01*psmo
      return
 endif

 call z_sum(is, ie, js, je, km, n_g, delp, q(is-n_g,js-n_g,1,sphum  ), psq(is,js,sphum  )) 
 call z_sum(is, ie, js, je, km, n_g, delp, q(is-n_g,js-n_g,1,liq_wat), psq(is,js,liq_wat))
 call z_sum(is, ie, js, je, km, n_g, delp, q(is-n_g,js-n_g,1,ice_wat), psq(is,js,ice_wat))

! Mean water vapor in the "stratosphere" (75 mb and above):
 kstrat = 1
 do k=1,km
    if ( phalf(k+1) > 75. ) exit
    kstrat = k
 enddo

 call z_sum(is, ie, js, je, kstrat, n_g, delp, q(is-n_g,js-n_g,1,sphum), q_strat(is,js)) 
 psmo = g_sum(q_strat(is,js), is, ie, js, je, n_g, area, 1) * 1.e6           &
      / p_sum(is, ie, js, je, kstrat, n_g, delp)

 if(master) write(*,*) 'Mean specific humidity (mg/kg) above 75 mb=', psmo, kstrat



 if ( nwat==6 ) then
     call z_sum(is, ie, js, je, km, n_g, delp, q(is-n_g,js-n_g,1,rainwat), psq(is,js,rainwat))
     call z_sum(is, ie, js, je, km, n_g, delp, q(is-n_g,js-n_g,1,snowwat), psq(is,js,snowwat))
     call z_sum(is, ie, js, je, km, n_g, delp, q(is-n_g,js-n_g,1,graupel), psq(is,js,graupel))
 endif

!-------------------
! Check global means
!-------------------
 psmo = g_sum(ps(is:ie,js:je), is, ie, js, je, n_g, area, 1) 

 do n=1,nwat
    qtot(n) = g_sum(psq(is,js,n), is, ie, js, je, n_g, area, 1) 
 enddo

 totw  = sum(qtot(1:nwat))
 psdry = psmo - totw

 if( master ) then
     write(6,*) 'Total surface pressure (mb) = ',  0.01*psmo
     write(6,*) 'mean dry surface pressure = ',    0.01*psdry
     write(6,*) 'Total Water Vapor (kg/m**2) =',  qtot(sphum)*ginv
     if ( nwat==6 ) then
          write(6,*) '--- Micro Phys water substances (kg/m**2) ---'
          write(6,*) 'Total cloud water=', qtot(liq_wat)*ginv
          write(6,*) 'Total rain  water=', qtot(rainwat)*ginv
          write(6,*) 'Total cloud ice  =', qtot(ice_wat)*ginv
          write(6,*) 'Total snow       =', qtot(snowwat)*ginv
          write(6,*) 'Total graupel    =', qtot(graupel)*ginv
          write(6,*) '---------------------------------------------'
     endif
  endif

#endif MARS_GCM

 end subroutine prt_mass


 subroutine z_sum(is, ie, js, je, km, n_g, delp, q, sum2)
 integer, intent(in):: is, ie, js, je,  n_g, km
 real, intent(in):: delp(is-n_g:ie+n_g, js-n_g:je+n_g, km)
 real, intent(in)::    q(is-n_g:ie+n_g, js-n_g:je+n_g, km)
 real, intent(out):: sum2(is:ie,js:je)
 integer i,j,k

 do j=js,je
    do i=is,ie
       sum2(i,j) = delp(i,j,1)*q(i,j,1)
    enddo
    do k=2,km
       do i=is,ie
          sum2(i,j) = sum2(i,j) + delp(i,j,k)*q(i,j,k)
       enddo
    enddo
 enddo

 end subroutine z_sum


 real function p_sum(is, ie, js, je, km, n_g, delp)
 integer, intent(in):: is, ie, js, je,  n_g, km
 real, intent(in):: delp(is-n_g:ie+n_g, js-n_g:je+n_g, km)
 real :: sum2(is:ie,js:je)
 integer i,j,k

 do j=js,je
    do i=is,ie
       sum2(i,j) = delp(i,j,1)
    enddo
    do k=2,km
       do i=is,ie
          sum2(i,j) = sum2(i,j) + delp(i,j,k)
       enddo
    enddo
 enddo
 p_sum = g_sum(sum2, is, ie, js, je, n_g, area, 1)

 end function p_sum



 subroutine get_pressure_given_height(is, ie, js, je, ng, km, wz, kd, height,   &
                                      ts, peln, a2, fac)

 integer,  intent(in):: is, ie, js, je, km, ng
 integer,  intent(in):: kd           ! vertical dimension of the ouput height
 real, intent(in):: wz(is:ie,js:je,km+1)
 real, intent(in):: ts(is-ng:ie+ng,js-ng:je+ng)
 real, intent(in):: peln(is:ie,km+1,js:je)
 real, intent(in):: height(kd)   ! must be monotonically decreasing with increasing k
 real, intent(out):: a2(is:ie,js:je,kd)      ! pressure (pa)
 real, optional, intent(in):: fac

! local:
 integer n, i,j,k
 real ptmp, tm


 do n=1,kd

!$omp parallel do default(shared) private(i, j, k, ptmp, tm)
    do j=js,je

       do 1000 i=is,ie

         if ( height(n) >= wz(i,j,km+1) ) then
!---------------------
! Search from top down
!---------------------
          do k=1,km
             if( height(n) < wz(i,j,k) .and. height(n) >= wz(i,j,k+1) ) then
! Found it!
                 ptmp = peln(i,k,j) + (peln(i,k+1,j)-peln(i,k,j)) *   &
                       (wz(i,j,k)-height(n)) / (wz(i,j,k)-wz(i,j,k+1))
                 a2(i,j,n) = exp(ptmp)
                 go to 500
             endif
          enddo

         else
!-----------------------------------------
! xtrapolation: mean laspe rate 6.5 deg/km
!-----------------------------------------
                tm = rdgas*ginv*(ts(i,j) + 3.25E-3*(wz(i,j,km)-height(n)))
          a2(i,j,n) = exp( peln(i,km+1,j) + (wz(i,j,km+1) - height(n))/tm )
         endif
500      if ( present(fac) ) a2(i,j,n) = fac * a2(i,j,n)
1000   continue
    enddo
 enddo

 end subroutine get_pressure_given_height


 subroutine get_height_given_pressure(is, ie, js, je, ng, km, wz, kd, log_p,   &
                                      peln, a2)
 integer,  intent(in):: is, ie, js, je, ng, km
 integer,  intent(in):: kd           ! vertical dimension of the ouput height
 real, intent(in):: log_p(kd)    ! must be monotonically decreasing with increasing k
                                 ! log (p)
 real, intent(in):: wz(is:ie,js:je,km+1)
 real, intent(in):: peln(is:ie,km+1,js:je)
 real, intent(out):: a2(is:ie,js:je,kd)      ! height (m)

! local:
 integer n, i,j,k

 do n=1,kd

!$omp parallel do default(shared) private(i, j, k)
    do j=js,je
       do 1000 i=is,ie
          do k=1,km
             if( log_p(n) <= peln(i,k+1,j) .and. log_p(n) >= peln(i,k,j) ) then
! Found it!
                 a2(i,j,n) = wz(i,j,k)  +  (wz(i,j,k+1) - wz(i,j,k)) *   &
                       (log_p(n)-peln(i,k,j)) / (peln(i,k+1,j)-peln(i,k,j) )
                 go to 1000
             endif
          enddo
!                a2(i,j,n) = missing_value
! Extrapolation into ground:  use wz(km-1:km+1)
                 a2(i,j,n) = wz(i,j,km+1) + (wz(i,j,km+1) - wz(i,j,km-1)) *   &
                       (log_p(n)-peln(i,km+1,j)) / (peln(i,km+1,j)-peln(i,km-1,j) )
1000   continue
    enddo
 enddo

 end subroutine get_height_given_pressure


 subroutine interpolate_vertical(is, ie, js, je, km, plev, peln, a3, a2)

 integer,  intent(in):: is, ie, js, je, km
 real, intent(in):: peln(is:ie,km+1,js:je)
 real, intent(in):: a3(is:ie,js:je,km)
 real, intent(in):: plev
 real, intent(out):: a2(is:ie,js:je)
! local:
 real pm(km)
 real logp
 integer i,j,k

 logp = log(plev)

 do j=js,je
    do 1000 i=is,ie

       do k=1,km
          pm(k) = 0.5*(peln(i,k,j)+peln(i,k+1,j))
       enddo

       if( logp <= pm(1) ) then
           a2(i,j) = a3(i,j,1)
       elseif ( logp >= pm(km) ) then
           a2(i,j) = a3(i,j,km)
       else 
           do k=1,km-1
              if( logp <= pm(k+1) .and. logp >= pm(k) ) then
                  a2(i,j) = a3(i,j,k) + (a3(i,j,k+1)-a3(i,j,k))*(logp-pm(k))/(pm(k+1)-pm(k))
                  go to 1000
              endif
           enddo
       endif
1000   continue
 enddo

 end subroutine interpolate_vertical



 subroutine pv_entropy(is, ie, js, je, ng, km, vort, f_d, pt, pkz, delp, grav)

! !INPUT PARAMETERS:
   integer, intent(in)::  is, ie, js, je, ng, km
   real, intent(in):: grav
   real, intent(in):: pt(is-ng:ie+ng,js-ng:je+ng,km) 
   real, intent(in):: pkz(is:ie,js:je,km) 
   real, intent(in):: delp(is-ng:ie+ng,js-ng:je+ng,km)
   real, intent(in):: f_d(is-ng:ie+ng,js-ng:je+ng) 

! vort is relative vorticity as input. Becomes PV on output
      real, intent(inout):: vort(is:ie,js:je,km)

! !DESCRIPTION:
!        EPV = 1/r * (vort+f_d) * d(S)/dz; where S is a conservative scalar
!        r the fluid density, and S is chosen to be the entropy here: S = log(pt)
!        pt == potential temperature.
! Computation are performed assuming the input is on "x-y-z" Cartesian coordinate.
! The approximation made here is that the relative vorticity computed on constant
! z-surface is not that different from the hybrid sigma-p coordinate.
! See page 39, Pedlosky 1979: Geophysical Fluid Dynamics
!
! The follwoing simplified form is strictly correct only if vort is computed on 
! constant z surfaces. In addition hydrostatic approximation is made.
!        EPV = - GRAV * (vort+f_d) / del(p) * del(pt) / pt 
! where del() is the vertical difference operator.
!
! programmer: S.-J. Lin, shian-jiann.lin@noaa.gov
!
!EOP
!---------------------------------------------------------------------
!BOC
      real w3d(is:ie,js:je,km)
      real te(is:ie,js:je,km+1), t2(is:ie,km), delp2(is:ie,km)
      real te2(is:ie,km+1)
      integer i, j, k

#ifdef SW_DYNAMICS
        do j=js,je
          do i=is,ie
            vort(i,j,1) = grav * (vort(i,j,1)+f_d(i,j)) / delp(i,j,1)
          enddo
        enddo
#else
! Compute PT at layer edges.
     do j=js,je
        do k=1,km
          do i=is,ie
               t2(i,k) = pt(i,j,k) / pkz(i,j,k)
              w3d(i,j,k) = t2(i,k)
            delp2(i,k) = delp(i,j,k)
          enddo
        enddo

        call ppme(t2, te2, delp2, ie-is+1, km)

        do k=1,km+1
           do i=is,ie
              te(i,j,k) = te2(i,k)
           enddo
        enddo
     enddo

     do k=1,km
        do j=js,je
          do i=is,ie
! Entropy is the thermodynamic variable in the following form
            vort(i,j,k) = (vort(i,j,k)+f_d(i,j)) * ( te(i,j,k)-te(i,j,k+1) )  &
                          / ( w3d(i,j,k)*delp(i,j,k) )  * grav
          enddo
        enddo
     enddo
#endif

 end subroutine pv_entropy


 subroutine ppme(p,qe,delp,im,km)

  integer, intent(in):: im, km
  real, intent(in)::    p(im,km)
  real, intent(in):: delp(im,km)
  real, intent(out)::qe(im,km+1)

! local arrays.
      real dc(im,km),delq(im,km), a6(im,km)
      real c1, c2, c3, tmp, qmax, qmin
      real a1, a2, s1, s2, s3, s4, ss3, s32, s34, s42
      real a3, b2, sc, dm, d1, d2, f1, f2, f3, f4
      real qm, dq
      integer i, k, km1

      km1 = km - 1

      do 500 k=2,km
      do 500 i=1,im
500   a6(i,k) = delp(i,k-1) + delp(i,k)

      do 1000 k=1,km1
      do 1000 i=1,im
      delq(i,k) = p(i,k+1) - p(i,k)
1000  continue

      do 1220 k=2,km1
      do 1220 i=1,im
      c1 = (delp(i,k-1)+0.5*delp(i,k))/a6(i,k+1)
      c2 = (delp(i,k+1)+0.5*delp(i,k))/a6(i,k)
      tmp = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /    &
                                    (a6(i,k)+delp(i,k+1))
      qmax = max(p(i,k-1),p(i,k),p(i,k+1)) - p(i,k)
      qmin = p(i,k) - min(p(i,k-1),p(i,k),p(i,k+1))
      dc(i,k) = sign(min(abs(tmp),qmax,qmin), tmp)
1220  continue

!****6***0*********0*********0*********0*********0*********0**********72
! 4th order interpolation of the provisional cell edge value
!****6***0*********0*********0*********0*********0*********0**********72

   do k=3,km1
      do i=1,im
         c1 = delq(i,k-1)*delp(i,k-1) / a6(i,k)
         a1 = a6(i,k-1) / (a6(i,k) + delp(i,k-1))
         a2 = a6(i,k+1) / (a6(i,k) + delp(i,k))
         qe(i,k) = p(i,k-1) + c1 + 2./(a6(i,k-1)+a6(i,k+1)) *        &
                   ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -         &
                                delp(i,k-1)*a1*dc(i,k  ) )
      enddo
   enddo

! three-cell parabolic subgrid distribution at model top

   do i=1,im
! three-cell PP-distribution
! Compute a,b, and c of q = aP**2 + bP + c using cell averages and delp
! a3 = a / 3
! b2 = b / 2
      s1 = delp(i,1)
      s2 = delp(i,2) + s1
!
      s3 = delp(i,2) + delp(i,3)
      s4 = s3 + delp(i,4)
      ss3 =  s3 + s1
      s32 = s3*s3
      s42 = s4*s4
      s34 = s3*s4
! model top
      a3 = (delq(i,2) - delq(i,1)*s3/s2) / (s3*ss3)
!
      if(abs(a3) .gt. 1.e-14) then
         b2 =  delq(i,1)/s2 - a3*(s1+s2)
         sc = -b2/(3.*a3)
         if(sc .lt. 0. .or. sc .gt. s1) then
             qe(i,1) = p(i,1) - s1*(a3*s1 + b2)
         else
             qe(i,1) = p(i,1) - delq(i,1)*s1/s2
         endif
      else
! Linear
         qe(i,1) = p(i,1) - delq(i,1)*s1/s2
      endif
      dc(i,1) = p(i,1) - qe(i,1)
! compute coef. for the off-centered area preserving cubic poly.
      dm = delp(i,1) / (s34*ss3*(delp(i,2)+s3)*(s4+delp(i,1)))
      f1 = delp(i,2)*s34 / ( s2*ss3*(s4+delp(i,1)) )
      f2 = (delp(i,2)+s3) * (ss3*(delp(i,2)*s3+s34+delp(i,2)*s4)   &
            + s42*(delp(i,2)+s3+s32/s2))
      f3 = -delp(i,2)*( ss3*(s32*(s3+s4)/(s4-delp(i,2))            &
            + (delp(i,2)*s3+s34+delp(i,2)*s4))                     &
            + s42*(delp(i,2)+s3) )
      f4 = ss3*delp(i,2)*s32*(delp(i,2)+s3) / (s4-delp(i,2))
      qe(i,2) = f1*p(i,1)+(f2*p(i,2)+f3*p(i,3)+f4*p(i,4))*dm
   enddo

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
   do i=1,im
      d1 = delp(i,km)
      d2 = delp(i,km1)
      qm = (d2*p(i,km)+d1*p(i,km1)) / (d1+d2)
      dq = 2.*(p(i,km1)-p(i,km)) / (d1+d2)
      c1 = (qe(i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
      c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1**2)
      qe(i,km  ) = qm - c1*d1*d2*(d2+3.*d1)
      qe(i,km+1) = d1*(8.*c1*d1**2-c3) + qe(i,km)
   enddo

 end subroutine ppme


end module fv_diagnostics_mod
