module fv_diagnostics_mod

 use constants_mod,    only: grav, rdgas, rvgas, pi, radius, kappa, WTMAIR, WTMCO2, omega
 use fms_io_mod,       only: set_domain, nullify_domain
 use time_manager_mod, only: time_type, get_date, get_time
 use mpp_domains_mod,  only: domain2d, mpp_update_domains, DGRID_NE
 use diag_manager_mod, only: diag_axis_init, register_diag_field, &
                             register_static_field, send_data, diag_grid_init
 use fv_arrays_mod,    only: fv_atmos_type, fv_grid_type, fv_diag_type
 !!! CLEANUP needs removal?
 use fv_mapz_mod,      only: E_Flux, E_Flux_nest
 use fv_mp_mod,        only: mp_reduce_sum, mp_reduce_min, mp_reduce_max
 use fv_eta_mod,        only: get_eta_level, gw_1d
 use fv_grid_utils_mod, only: g_sum
 use a2b_edge_mod,     only: a2b_ord2, a2b_ord4
 use fv_surf_map_mod,  only: zs_g
 use fv_sg_mod,        only: qsmith

 use tracer_manager_mod, only: get_tracer_names, get_number_tracers, get_tracer_index
 use field_manager_mod,  only: MODEL_ATMOS
 use mpp_mod,            only: mpp_error, FATAL, stdlog, mpp_pe, mpp_root_pe
 use sat_vapor_pres_mod, only: compute_qs, lookup_es

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

 use fv_arrays_mod, only: max_step !Parameter; should be OK

 implicit none
 private


 real :: missing_value = -1.e10
 real :: ginv
 real :: pk0
 logical master

! private (to this module) diag:

 type(time_type) :: fv_time
 type(fv_diag_type), pointer :: idiag
 character(len=3) :: gn = ''

 logical :: module_is_initialized=.false.
 logical :: prt_minmax =.false.
 logical :: m_calendar
 integer  sphum, liq_wat, ice_wat       ! GFDL physics
 integer  rainwat, snowwat, graupel
 integer :: istep
 real    :: ptop
 real    :: rad2deg
! tracers
 character(len=128)   :: tname
 character(len=256)   :: tlongname, tunits

 public :: fv_diag_init, fv_time, fv_diag, prt_mxm, prt_maxmin, range_check!, id_divg, id_te
 public :: prt_mass, prt_minmax, ppme

!---- version number -----
 character(len=128) :: version = '$Id: fv_diagnostics.F90,v 20.0.2.1.2.1.2.1.2.1.4.1 2014/11/24 19:16:18 Lucas.Harris Exp $'
 character(len=128) :: tagname = '$Name: testing $'

contains

 subroutine fv_diag_init(Atm, axes, Time, npx, npy, npz, p_ref)
    type(fv_atmos_type), intent(inout), target :: Atm(:)
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

    !These id_* are not needed later since they are for static data which is not used elsewhere
    integer :: id_bk, id_pk, id_area, id_lon, id_lat, id_lont, id_latt, id_phalf, id_pfull
    integer :: id_hyam, id_hybm
    integer :: i, j, k, n, ntileMe, id_xt, id_yt, id_x, id_y, id_xe, id_ye, id_xn, id_yn
    integer :: isc, iec, jsc, jec

    logical :: used

    character(len=64) :: field
    integer              :: ntprog
    integer              :: unit

    integer :: ncnst

    if (Atm(1)%grid_Number > 1) then
       write(gn,"(A2,I1)") " g", Atm(1)%grid_number
    else
       gn = ""
    end if

    idiag => Atm(1)%idiag

    rad2deg = 180./pi

! For total energy diagnostics:
    idiag%steps = 0
    idiag%efx = 0.;       idiag%efx_sum = 0.
    idiag%mtq = 0.;       idiag%mtq_sum = 0.

    ncnst = Atm(1)%ncnst
    m_calendar = Atm(1)%flagstruct%moist_phys

    call set_domain(Atm(1)%domain)  ! Set domain so that diag_manager can access tile information

         sphum   = get_tracer_index (MODEL_ATMOS, 'sphum')
    if ( Atm(1)%flagstruct%nwat>=3 ) then
         liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
         ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
    endif

    if ( Atm(1)%flagstruct%nwat==6 ) then
        rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat')
        snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat')
        graupel = get_tracer_index (MODEL_ATMOS, 'graupel')
     endif

! valid range for some fields

!!!  This will need mods for more than 1 tile per pe  !!!

    vsrange = (/ -200.,  200. /)  ! surface (lowest layer) winds

    vrange = (/ -330.,  330. /)  ! winds
    wrange = (/ -100.,  100. /)  ! vertical wind
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
     if (Atm(1)%grid_number == 1) fv_time = Time

    allocate ( idiag%phalf(npz+1) )
    call get_eta_level(Atm(1)%npz, p_ref, pfull, idiag%phalf, Atm(1)%ak, Atm(1)%bk, 0.01)

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
    isc = Atm(n)%bd%isc; iec = Atm(n)%bd%iec
    jsc = Atm(n)%bd%jsc; jec = Atm(n)%bd%jec

    ! Send diag_manager the grid informtaion
    call diag_grid_init(DOMAIN=Atm(n)%domain, &
         &              GLO_LON=rad2deg*Atm(n)%gridstruct%grid(isc:iec+1,jsc:jec+1,1), &
         &              GLO_LAT=rad2deg*Atm(n)%gridstruct%grid(isc:iec+1,jsc:jec+1,2), &
         &              AGLO_LON=rad2deg*Atm(n)%gridstruct%agrid(isc-1:iec+1,jsc-1:jec+1,1), &
         &              AGLO_LAT=rad2deg*Atm(n)%gridstruct%agrid(isc-1:iec+1,jsc-1:jec+1,2))

    ntileMe = size(Atm(:))
    if (ntileMe > 1) call mpp_error(FATAL, "fv_diag_init can only be called with one grid at a time.")

!    do n = 1, ntileMe
    n = 1
       field = 'grid'

       id_xt = diag_axis_init('grid_xt',grid_xt,'degrees_E','x','T-cell longitude', &
                           set_name=trim(field),Domain2=Atm(n)%Domain, tile_count=n)
       id_yt = diag_axis_init('grid_yt',grid_yt,'degrees_N','y','T-cell latitude',  &
                           set_name=trim(field), Domain2=Atm(n)%Domain, tile_count=n)
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
                           set_name=trim(field),Domain2=Atm(n)%Domain, tile_count=n)
       id_y = diag_axis_init('grid_y',grid_y,'degrees_N','y','cell corner latitude',  &
                           set_name=trim(field), Domain2=Atm(n)%Domain, tile_count=n)

!    end do
!   deallocate(grid_xt, grid_yt, grid_xe, grid_ye, grid_xn, grid_yn)
    deallocate(grid_xt, grid_yt)
    deallocate(grid_x,  grid_y )

    id_phalf = diag_axis_init('phalf', idiag%phalf, 'mb', 'z', &
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

!    do n = 1, ntileMe
    n = 1
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
       idiag%id_zsurf = register_static_field ( trim(field), 'zsurf', axes(1:2),  &
                                         'surface height', 'm' )
#endif
       idiag%id_zs = register_static_field ( trim(field), 'zs', axes(1:2),  &
                                        'Original Mean Terrain', 'm' )
! 3D hybrid_z fields:
       idiag%id_ze = register_static_field ( trim(field), 'ze', axes(1:3),  &
                                        'Hybrid_Z_surface', 'm' )
       idiag%id_oro = register_static_field ( trim(field), 'oro', axes(1:2),  &
                                        'Land/Water Mask', 'none' )
       idiag%id_sgh = register_static_field ( trim(field), 'sgh', axes(1:2),  &
                                        'Terrain Standard deviation', 'm' )
!       idiag%id_ts = register_static_field ( trim(field), 'ts', axes(1:2),  &
!                                        'Skin temperature', 'K' )

!--------------------
! Initial conditions:
!--------------------
       idiag%ic_ps  = register_static_field ( trim(field), 'ps_ic', axes(1:2),  &
                                         'initial surface pressure', 'Pa' )
       idiag%ic_ua = register_static_field ( trim(field), 'ua_ic', axes(1:3),        &
            'zonal wind', 'm/sec' )
       idiag%ic_va = register_static_field ( trim(field), 'va_ic', axes(1:3),        &
            'meridional wind', 'm/sec' )
       idiag%ic_ppt= register_static_field ( trim(field), 'ppt_ic', axes(1:3),        &
            'potential temperature perturbation', 'K' )
#ifdef LASPRAT
       idiag%ic_sphum  = register_static_field ( trim(field), 'sphum_ic', axes(1:2),  &
                                         'initial surface pressure', 'Pa' )
#endif

!    end do

    master = (mpp_pe()==mpp_root_pe())

    n=1
    isc = Atm(n)%bd%isc; iec = Atm(n)%bd%iec
    jsc = Atm(n)%bd%jsc; jec = Atm(n)%bd%jec

    allocate ( idiag%zsurf(isc:iec,jsc:jec) )

    do j=jsc,jec
       do i=isc,iec
          idiag%zsurf(i,j) = ginv * Atm(n)%phis(i,j)
       enddo
    enddo

!--- Send time independent data

!    do n = 1, ntileMe
    n = 1
       isc = Atm(n)%bd%isc; iec = Atm(n)%bd%iec
       jsc = Atm(n)%bd%jsc; jec = Atm(n)%bd%jec
       if (id_lon  > 0) used = send_data(id_lon,  180./pi*Atm(n)%gridstruct%grid(isc:iec+1,jsc:jec+1,1), Time)
       if (id_lat  > 0) used = send_data(id_lat,  180./pi*Atm(n)%gridstruct%grid(isc:iec+1,jsc:jec+1,2), Time)
       if (id_lont > 0) used = send_data(id_lont, 180./pi*Atm(n)%gridstruct%agrid(isc:iec,jsc:jec,1), Time)
       if (id_latt > 0) used = send_data(id_latt, 180./pi*Atm(n)%gridstruct%agrid(isc:iec,jsc:jec,2), Time)
       if (id_area > 0) used = send_data(id_area, Atm(n)%gridstruct%area(isc:iec,jsc:jec), Time)
#ifndef DYNAMICS_ZS
       if (idiag%id_zsurf > 0) used = send_data(idiag%id_zsurf, idiag%zsurf, Time)
#endif
       if ( Atm(n)%flagstruct%fv_land ) then
         if (idiag%id_zs  > 0) used = send_data(idiag%id_zs , zs_g, Time)
         if (idiag%id_oro > 0) used = send_data(idiag%id_oro, Atm(n)%oro(isc:iec,jsc:jec), Time)
         if (idiag%id_sgh > 0) used = send_data(idiag%id_sgh, Atm(n)%sgh(isc:iec,jsc:jec), Time)
       endif

       if ( Atm(n)%flagstruct%ncep_ic ) then
         if (idiag%id_ts > 0) used = send_data(idiag%id_ts, Atm(n)%ts(isc:iec,jsc:jec), Time)
       endif

       if ( Atm(n)%flagstruct%hybrid_z .and. idiag%id_ze > 0 ) &
                      used = send_data(idiag%id_ze, Atm(n)%ze0(isc:iec,jsc:jec,1:npz), Time)

       if (idiag%ic_ps > 0) used = send_data(idiag%ic_ps, Atm(n)%ps(isc:iec,jsc:jec)*ginv, Time)

       if(idiag%ic_ua > 0) used=send_data(idiag%ic_ua, Atm(n)%ua(isc:iec,jsc:jec,:), Time)
       if(idiag%ic_va > 0) used=send_data(idiag%ic_va, Atm(n)%va(isc:iec,jsc:jec,:), Time)

       pk0 = 1000.E2 ** kappa
       if(idiag%ic_ppt> 0) then
! Potential temperature
          allocate ( idiag%pt1(npz) )
          allocate ( a3(isc:iec,jsc:jec,npz) )
#ifdef TEST_GWAVES
          call gw_1d(npz, 1000.E2, Atm(n)%ak, Atm(n)%ak, Atm(n)%ak(1), 10.E3, idiag%pt1)
#else
          idiag%pt1 = 0.
#endif
          do k=1,npz
          do j=jsc,jec
             do i=isc,iec
                a3(i,j,k) =  (Atm(n)%pt(i,j,k)/Atm(n)%pkz(i,j,k) - idiag%pt1(k)) * pk0
             enddo
          enddo
          enddo
          used=send_data(idiag%ic_ppt, a3, Time)
          deallocate ( a3 )
          deallocate ( idiag%pt1 )
       endif
!    end do

!--------------------------------------------------------------
! Register main prognostic fields: ps, (u,v), t, omega (dp/dt)
!--------------------------------------------------------------

    allocate(idiag%id_tracer(ncnst))
    allocate(idiag%id_tracer_dmmr(ncnst))
    allocate(idiag%id_tracer_dvmr(ncnst))
    allocate(idiag%w_mr(ncnst))
    idiag%id_tracer(:)      = 0
    idiag%id_tracer_dmmr(:) = 0 
    idiag%id_tracer_dvmr(:) = 0 
    idiag%w_mr(:) = 0.E0

!    do n = 1, ntileMe
    n = 1
       field= 'dynamics'

#ifdef DYNAMICS_ZS
       idiag%id_zsurf = register_diag_field ( trim(field), 'zsurf', axes(1:2), Time,           &
                                       'surface height', 'm')
#endif
!-------------------
! Surface pressure
!-------------------
       idiag%id_ps = register_diag_field ( trim(field), 'ps', axes(1:2), Time,           &
            'surface pressure', 'Pa', missing_value=missing_value )

!-------------------
! Mountain torque
!-------------------
       idiag%id_mq = register_diag_field ( trim(field), 'mq', axes(1:2), Time,           &
            'mountain torque', 'Hadleys per unit area', missing_value=missing_value )
!-------------------
! Angular momentum
!-------------------
       idiag%id_aam = register_diag_field ( trim(field), 'aam', axes(1:2), Time,           &
            'angular momentum', 'kg*m^2/s', missing_value=missing_value )
       idiag%id_amdt = register_diag_field ( trim(field), 'amdt', axes(1:2), Time,           &
            'angular momentum error', 'kg*m^2/s^2', missing_value=missing_value )


!--------------
! 10 mb Height
!--------------
      idiag%id_h10 = register_diag_field (trim(field), 'h10', axes(1:2),  Time,   &
                                     '10-mb hght', 'm', missing_value=missing_value )
!--------------
! 50 mb Height
!--------------
      idiag%id_h50 = register_diag_field (trim(field), 'h50', axes(1:2),  Time,   &
                                     '50-mb hght', 'm', missing_value=missing_value )
!--------------
! 100 mb Height
!--------------
      idiag%id_h100 = register_diag_field (trim(field), 'h100', axes(1:2),  Time,   &
                                     '100-mb hght', 'm', missing_value=missing_value )
!--------------
! 200 mb Height
!--------------
      idiag%id_h200 = register_diag_field (trim(field), 'h200', axes(1:2),  Time,   &
                                     '200-mb hght', 'm', missing_value=missing_value )
!--------------
! 250 mb Height
!--------------
      idiag%id_h250 = register_diag_field (trim(field), 'h250', axes(1:2),  Time,   &
                                     '250-mb hght', 'm', missing_value=missing_value )
!--------------
! 300 mb Height
!--------------
      idiag%id_h300 = register_diag_field (trim(field), 'h300', axes(1:2),  Time,   &
                                     '300-mb hght', 'm', missing_value=missing_value )
!--------------
! 500 mb Height
!--------------
      idiag%id_h500 = register_diag_field (trim(field), 'h500', axes(1:2),  Time,   &
                                     '500-mb hght', 'm', missing_value=missing_value )
!--------------
! 700 mb Height
!--------------
      idiag%id_h700 = register_diag_field (trim(field), 'h700', axes(1:2),  Time,   &
                                     '700-mb hght', 'm', missing_value=missing_value )
!--------------
! 850 mb Height
!--------------
      idiag%id_h850 = register_diag_field (trim(field), 'h850', axes(1:2),  Time,   &
                                     '850-mb hght', 'm', missing_value=missing_value )
!--------------
! 1000 mb Height
!--------------
      idiag%id_h1000= register_diag_field (trim(field), 'h1000', axes(1:2),  Time,   &
                                     '1000-mb hght', 'm', missing_value=missing_value )

      ! flag for calculation of geopotential
      if ( idiag%id_h10>0  .or. idiag%id_h50>0  .or. idiag%id_h100>0 .or. idiag%id_h200>0 .or. idiag%id_h250>0 .or. &
           idiag%id_h300>0 .or. idiag%id_h500>0 .or. idiag%id_h700>0 .or. idiag%id_h850>0 .or. idiag%id_h1000>0 ) then
           idiag%id_hght = 1
      else
           idiag%id_hght = 0
      endif
!-----------------------------
! mean temp between 300-500 mb
!-----------------------------
      idiag%id_tm = register_diag_field (trim(field), 'tm', axes(1:2),  Time,   &
                                   'mean 300-500 mb temp', 'K', missing_value=missing_value )

!-------------------
! Sea-level-pressure
!-------------------
       idiag%id_slp = register_diag_field (trim(field), 'slp', axes(1:2),  Time,   &
                                     'sea-level pressure', 'mb', missing_value=missing_value,  &
                                      range=slprange )
!----------------------------------
! Bottom level pressure for masking
!----------------------------------
       idiag%id_pmask = register_diag_field (trim(field), 'pmask', axes(1:2),  Time,   &
                                     'masking pressure at lowest level', 'mb',   &
                                      missing_value=missing_value )
!------------------------------------------
! Fix for Bottom level pressure for masking
!------------------------------------------
       idiag%id_pmaskv2 = register_diag_field(TRIM(field), 'pmaskv2', axes(1:2), Time,&
            & 'masking pressure at lowest level', 'mb', missing_value=missing_value)
                                     
!-------------------
! Hurricane scales:
!-------------------
! Net effects: ~ intensity * freq
       idiag%id_c15 = register_diag_field (trim(field), 'cat15', axes(1:2),  Time,   &
                                     'de-pression < 1000', 'mb', missing_value=missing_value)
       idiag%id_c25 = register_diag_field (trim(field), 'cat25', axes(1:2),  Time,   &
                                     'de-pression < 980', 'mb', missing_value=missing_value)
       idiag%id_c35 = register_diag_field (trim(field), 'cat35', axes(1:2),  Time,   &
                                     'de-pression < 964', 'mb', missing_value=missing_value)
       idiag%id_c45 = register_diag_field (trim(field), 'cat45', axes(1:2),  Time,   &
                                     'de-pression < 944', 'mb', missing_value=missing_value)
! Frequency:
       idiag%id_f15 = register_diag_field (trim(field), 'f15', axes(1:2),  Time,   &
                                     'Cat15 frequency', 'none', missing_value=missing_value)
       idiag%id_f25 = register_diag_field (trim(field), 'f25', axes(1:2),  Time,   &
                                     'Cat25 frequency', 'none', missing_value=missing_value)
       idiag%id_f35 = register_diag_field (trim(field), 'f35', axes(1:2),  Time,   &
                                     'Cat35 frequency', 'none', missing_value=missing_value)
       idiag%id_f45 = register_diag_field (trim(field), 'f45', axes(1:2),  Time,   &
                                     'Cat45 frequency', 'none', missing_value=missing_value)
!-------------------
! A grid winds (lat-lon)
!-------------------
       idiag%id_ua = register_diag_field ( trim(field), 'ucomp', axes(1:3), Time,        &
            'zonal wind', 'm/sec', missing_value=missing_value, range=vrange )
       idiag%id_va = register_diag_field ( trim(field), 'vcomp', axes(1:3), Time,        &
            'meridional wind', 'm/sec', missing_value=missing_value, range=vrange)

       idiag%id_w = register_diag_field ( trim(field), 'w', axes(1:3), Time,        &
            'vertical wind', 'm/sec', missing_value=missing_value, range=wrange )

       idiag%id_pt   = register_diag_field ( trim(field), 'temp', axes(1:3), Time,       &
            'temperature', 'K', missing_value=missing_value, range=trange )
       idiag%id_ppt  = register_diag_field ( trim(field), 'ppt', axes(1:3), Time,       &
            'potential temperature perturbation', 'K', missing_value=missing_value )
       idiag%id_omga = register_diag_field ( trim(field), 'omega', axes(1:3), Time,      &
            'omega', 'Pa/s', missing_value=missing_value )
       idiag%id_divg  = register_diag_field ( trim(field), 'divg', axes(1:3), Time,      &
            'mean divergence', '1/s', missing_value=missing_value )

       idiag%id_rh = register_diag_field ( trim(field), 'rh', axes(1:3), Time,        &
            'Relative Humidity', '%', missing_value=missing_value )
!            'Relative Humidity', '%', missing_value=missing_value, range=rhrange )
! Total energy (only when moist_phys = .T.)
       idiag%id_te    = register_diag_field ( trim(field), 'te', axes(1:2), Time,      &
            'Total Energy', 'J/kg', missing_value=missing_value )
! Total Kinetic energy
       idiag%id_ke    = register_diag_field ( trim(field), 'ke', axes(1:2), Time,      &
            'Total KE', 'm^2/s^2', missing_value=missing_value )
       idiag%id_delp = register_diag_field ( trim(field), 'delp', axes(1:3), Time,        &
            'pressure thickness', 'pa', missing_value=missing_value )
       idiag%id_delz = register_diag_field ( trim(field), 'delz', axes(1:3), Time,        &
            'height thickness', 'm', missing_value=missing_value )
       idiag%id_zratio = register_diag_field ( trim(field), 'zratio', axes(1:3), Time,        &
            'nonhydro_ratio', 'n/a', missing_value=missing_value )
       idiag%id_ws     = register_diag_field ( trim(field), 'ws', axes(1:2), Time,        &
            'Terrain W', 'm/s', missing_value=missing_value )
!--------------------
! 3D Condensate
!--------------------
       idiag%id_qn = register_diag_field ( trim(field), 'qn', axes(1:3), Time,       &
               'cloud condensate', 'g/g', missing_value=missing_value )
       idiag%id_qp = register_diag_field ( trim(field), 'qp', axes(1:3), Time,       &
               'precip condensate', 'g/g', missing_value=missing_value )
! fast moist phys tendencies:
       idiag%id_mdt = register_diag_field ( trim(field), 'mdt', axes(1:3), Time,       &
               'DT/Dt: fast moist phys', 'deg/sec', missing_value=missing_value )
!--------------------
! Relative vorticity
!--------------------
       idiag%id_vort = register_diag_field ( trim(field), 'vort', axes(1:3), Time,       &
            'vorticity', '1/s', missing_value=missing_value )
!--------------------
! Potential vorticity
!--------------------
       idiag%id_pv = register_diag_field ( trim(field), 'pv', axes(1:3), Time,       &
            'potential vorticity', '1/s', missing_value=missing_value )

#ifdef MARS_GCM
!--------------------------
! Extra Martian diagnostics:
!--------------------------

       idiag%id_t05 = register_diag_field ( trim(field), 't05', axes(1:2), Time,       &
               '0.5-mb temperature', 'K', missing_value=missing_value )
!!       idiag%id_sfc_dust = register_diag_field ( trim(field), 'sfc_dust', axes(1:2), Time,        &
!!             'Total sfc dust', 'kg/m**2', missing_value=missing_value )
!!        idiag%id_tdust = register_diag_field ( trim(field), 'odcol', axes(1:2), Time,        &
!!             'Total dust column', 'kg/m**2', missing_value=missing_value )
#endif MARS_GCM

!--------------------------
! Extra surface diagnistics:
!--------------------------
! Surface (lowest layer) vorticity: for tropical cyclones diag.
       idiag%id_vorts = register_diag_field ( trim(field), 'vorts', axes(1:2), Time,       &
            'surface vorticity', '1/s', missing_value=missing_value )
       idiag%id_us = register_diag_field ( trim(field), 'us', axes(1:2), Time,        &
            'surface u-wind', 'm/sec', missing_value=missing_value, range=vsrange )
       idiag%id_vs = register_diag_field ( trim(field), 'vs', axes(1:2), Time,        &
            'surface v-wind', 'm/sec', missing_value=missing_value, range=vsrange )
       idiag%id_tq = register_diag_field ( trim(field), 'tq', axes(1:2), Time,        &
            'Total water path', 'kg/m**2', missing_value=missing_value )
       idiag%id_iw = register_diag_field ( trim(field), 'iw', axes(1:2), Time,        &
            'Ice water path', 'kg/m**2', missing_value=missing_value )
       idiag%id_lw = register_diag_field ( trim(field), 'lw', axes(1:2), Time,        &
            'Liquid water path', 'kg/m**2', missing_value=missing_value )
       idiag%id_ts = register_diag_field ( trim(field), 'ts', axes(1:2), Time,  &
                                        'Skin temperature', 'K' )
       idiag%id_tb = register_diag_field ( trim(field), 'tb', axes(1:2), Time,  &
                                        'lowest layer temperature', 'K' )

!--------------------------
! 850-mb vorticity
!--------------------------
       idiag%id_vort850 = register_diag_field ( trim(field), 'vort850', axes(1:2), Time,       &
                           '850-mb vorticity', '1/s', missing_value=missing_value )

!--------------------------
! 10-mb winds:
!--------------------------
       idiag%id_u10 = register_diag_field ( trim(field), 'u10', axes(1:2), Time,       &
                           '10-mb u-wind', '1/s', missing_value=missing_value )
       idiag%id_v10 = register_diag_field ( trim(field), 'v10', axes(1:2), Time,       &
                           '10-mb v-wind', '1/s', missing_value=missing_value )
!--------------------------
! 50-mb winds:
!--------------------------
       idiag%id_u50 = register_diag_field ( trim(field), 'u50', axes(1:2), Time,       &
                           '50-mb u-wind', '1/s', missing_value=missing_value )
       idiag%id_v50 = register_diag_field ( trim(field), 'v50', axes(1:2), Time,       &
                           '50-mb v-wind', '1/s', missing_value=missing_value )
!--------------------------
! 100-mb winds:
!--------------------------
       idiag%id_u100 = register_diag_field ( trim(field), 'u100', axes(1:2), Time,       &
                           '100-mb u-wind', '1/s', missing_value=missing_value )
       idiag%id_v100 = register_diag_field ( trim(field), 'v100', axes(1:2), Time,       &
                           '100-mb v-wind', '1/s', missing_value=missing_value )
!--------------------------
! 200-mb winds:
!--------------------------
       idiag%id_u200 = register_diag_field ( trim(field), 'u200', axes(1:2), Time,       &
                           '200-mb u-wind', '1/s', missing_value=missing_value )
       idiag%id_v200 = register_diag_field ( trim(field), 'v200', axes(1:2), Time,       &
                           '200-mb v-wind', '1/s', missing_value=missing_value )
       idiag%id_w200 = register_diag_field ( trim(field), 'w200', axes(1:2), Time,       &
                           '200-mb w-wind', '1/s', missing_value=missing_value )
! s200: wind speed for computing KE spectrum
! Cubed_2_latlon interpolation is more accurate, particularly near the poles, using
! winds speed (a scalar), rather than wind vectors or kinetic energy directly.
       idiag%id_s200 = register_diag_field ( trim(field), 's200', axes(1:2), Time,       &
                           '200-mb wind_speed', 'm/s', missing_value=missing_value )
       idiag%id_sl12 = register_diag_field ( trim(field), 'sl12', axes(1:2), Time,       &
                           '12th L wind_speed', 'm/s', missing_value=missing_value )
       idiag%id_sl13 = register_diag_field ( trim(field), 'sl13', axes(1:2), Time,       &
                           '13th L wind_speed', 'm/s', missing_value=missing_value )
!--------------------------
! 250-mb winds:
!--------------------------
       idiag%id_u250 = register_diag_field ( trim(field), 'u250', axes(1:2), Time,       &
                           '250-mb u-wind', '1/s', missing_value=missing_value )
       idiag%id_v250 = register_diag_field ( trim(field), 'v250', axes(1:2), Time,       &
                           '250-mb v-wind', '1/s', missing_value=missing_value )
!--------------------------
! 300-mb winds:
!--------------------------
       idiag%id_u300 = register_diag_field ( trim(field), 'u300', axes(1:2), Time,       &
                           '300-mb u-wind', '1/s', missing_value=missing_value )
       idiag%id_v300 = register_diag_field ( trim(field), 'v300', axes(1:2), Time,       &
                           '300-mb v-wind', '1/s', missing_value=missing_value )
!--------------------------
! 500-mb winds:
!--------------------------
       idiag%id_u500 = register_diag_field ( trim(field), 'u500', axes(1:2), Time,       &
                           '500-mb u-wind', '1/s', missing_value=missing_value )
       idiag%id_v500 = register_diag_field ( trim(field), 'v500', axes(1:2), Time,       &
                           '500-mb v-wind', '1/s', missing_value=missing_value )
!--------------------------
! 700-mb winds:
!--------------------------
       idiag%id_u700 = register_diag_field ( trim(field), 'u700', axes(1:2), Time,       &
                           '700-mb u-wind', '1/s', missing_value=missing_value )
       idiag%id_v700 = register_diag_field ( trim(field), 'v700', axes(1:2), Time,       &
                           '700-mb v-wind', '1/s', missing_value=missing_value )
!--------------------------
! 850-mb winds:
!--------------------------
       idiag%id_u850 = register_diag_field ( trim(field), 'u850', axes(1:2), Time,       &
                           '850-mb u-wind', '1/s', missing_value=missing_value )
       idiag%id_v850 = register_diag_field ( trim(field), 'v850', axes(1:2), Time,       &
                           '850-mb v-wind', '1/s', missing_value=missing_value )
       idiag%id_w850 = register_diag_field ( trim(field), 'w850', axes(1:2), Time,       &
                           '850-mb w-wind', '1/s', missing_value=missing_value )
       idiag%id_w5km = register_diag_field ( trim(field), 'w5km', axes(1:2), Time,       &
                           '5-km w-wind', '1/s', missing_value=missing_value )
! helicity
       idiag%id_x850 = register_diag_field ( trim(field), 'x850', axes(1:2), Time,       &
                           '850-mb vertical comp. of helicity', 'm/s**2', missing_value=missing_value )

! Storm Relative Helicity
       idiag%id_srh = register_diag_field ( trim(field), 'srh', axes(1:2), Time,       &
                           'Storm Relative Helicity', 'm/s**2', missing_value=missing_value )
!--------------------------
! 1000-mb winds:
!--------------------------
       idiag%id_u1000 = register_diag_field ( trim(field), 'u1000', axes(1:2), Time,       &
                           '1000-mb u-wind', '1/s', missing_value=missing_value )
       idiag%id_v1000 = register_diag_field ( trim(field), 'v1000', axes(1:2), Time,       &
                           '1000-mb v-wind', '1/s', missing_value=missing_value )
!--------------------------
! temperature:
!--------------------------
       idiag%id_t10 = register_diag_field ( trim(field), 't10', axes(1:2), Time,       &
                           '10-mb temperature', 'K', missing_value=missing_value )
       idiag%id_t50 = register_diag_field ( trim(field), 't50', axes(1:2), Time,       &
                           '50-mb temperature', 'K', missing_value=missing_value )
       idiag%id_t100 = register_diag_field ( trim(field), 't100', axes(1:2), Time,       &
                           '100-mb temperature', 'K', missing_value=missing_value )
       idiag%id_t200 = register_diag_field ( trim(field), 't200', axes(1:2), Time,       &
                           '200-mb temperature', 'K', missing_value=missing_value )
       idiag%id_t250 = register_diag_field ( trim(field), 't250', axes(1:2), Time,       &
                           '250-mb temperature', 'K', missing_value=missing_value )
       idiag%id_t300 = register_diag_field ( trim(field), 't300', axes(1:2), Time,       &
                           '300-mb temperature', 'K', missing_value=missing_value )
       idiag%id_t500 = register_diag_field ( trim(field), 't500', axes(1:2), Time,       &
                           '500-mb temperature', 'K', missing_value=missing_value )
       idiag%id_t700 = register_diag_field ( trim(field), 't700', axes(1:2), Time,       &
                           '700-mb temperature', 'K', missing_value=missing_value )
       idiag%id_t850 = register_diag_field ( trim(field), 't850', axes(1:2), Time,       &
                           '850-mb temperature', 'K', missing_value=missing_value )
       idiag%id_t1000 = register_diag_field ( trim(field), 't1000', axes(1:2), Time,       &
                           '1000-mb temperature', 'K', missing_value=missing_value )
!--------------------------
! specific humidity:
!--------------------------
       idiag%id_q10 = register_diag_field ( trim(field), 'q10', axes(1:2), Time,       &
                           '10-mb specific humidity', 'kg/kg', missing_value=missing_value )
       idiag%id_q50 = register_diag_field ( trim(field), 'q50', axes(1:2), Time,       &
                           '50-mb specific humidity', 'kg/kg', missing_value=missing_value )
       idiag%id_q100 = register_diag_field ( trim(field), 'q100', axes(1:2), Time,       &
                           '100-mb specific humidity', 'kg/kg', missing_value=missing_value )
       idiag%id_q200 = register_diag_field ( trim(field), 'q200', axes(1:2), Time,       &
                           '200-mb specific humidity', 'kg/kg', missing_value=missing_value )
       idiag%id_q250 = register_diag_field ( trim(field), 'q250', axes(1:2), Time,       &
                           '250-mb specific humidity', 'kg/kg', missing_value=missing_value )
       idiag%id_q300 = register_diag_field ( trim(field), 'q300', axes(1:2), Time,       &
                           '300-mb specific humidity', 'kg/kg', missing_value=missing_value )
       idiag%id_q500 = register_diag_field ( trim(field), 'q500', axes(1:2), Time,       &
                           '500-mb specific humidity', 'kg/kg', missing_value=missing_value )
       idiag%id_q700 = register_diag_field ( trim(field), 'q700', axes(1:2), Time,       &
                           '700-mb specific humidity', 'kg/kg', missing_value=missing_value )
       idiag%id_q850 = register_diag_field ( trim(field), 'q850', axes(1:2), Time,       &
                           '850-mb specific humidity', 'kg/kg', missing_value=missing_value )
       idiag%id_q1000 = register_diag_field ( trim(field), 'q1000', axes(1:2), Time,       &
                           '1000-mb specific humidity', 'kg/kg', missing_value=missing_value )
!--------------------------
! relative humidity (physics definition):
!--------------------------
       idiag%id_rh10 = register_diag_field ( trim(field), 'rh10', axes(1:2), Time,       &
                           '10-mb relative humidity', '%', missing_value=missing_value )
       idiag%id_rh50 = register_diag_field ( trim(field), 'rh50', axes(1:2), Time,       &
                           '50-mb relative humidity', '%', missing_value=missing_value )
       idiag%id_rh100 = register_diag_field ( trim(field), 'rh100', axes(1:2), Time,       &
                           '100-mb relative humidity', '%', missing_value=missing_value )
       idiag%id_rh200 = register_diag_field ( trim(field), 'rh200', axes(1:2), Time,       &
                           '200-mb relative humidity', '%', missing_value=missing_value )
       idiag%id_rh250 = register_diag_field ( trim(field), 'rh250', axes(1:2), Time,       &
                           '250-mb relative humidity', '%', missing_value=missing_value )
       idiag%id_rh300 = register_diag_field ( trim(field), 'rh300', axes(1:2), Time,       &
                           '300-mb relative humidity', '%', missing_value=missing_value )
       idiag%id_rh500 = register_diag_field ( trim(field), 'rh500', axes(1:2), Time,       &
                           '500-mb relative humidity', '%', missing_value=missing_value )
       idiag%id_rh700 = register_diag_field ( trim(field), 'rh700', axes(1:2), Time,       &
                           '700-mb relative humidity', '%', missing_value=missing_value )
       idiag%id_rh850 = register_diag_field ( trim(field), 'rh850', axes(1:2), Time,       &
                           '850-mb relative humidity', '%', missing_value=missing_value )
       idiag%id_rh1000 = register_diag_field ( trim(field), 'rh1000', axes(1:2), Time,       &
                           '1000-mb relative humidity', '%', missing_value=missing_value )
!--------------------------
! relative humidity (CMIP definition):
!--------------------------
       idiag%id_rh10_cmip = register_diag_field ( trim(field), 'rh10_cmip', axes(1:2), Time,       &
                           '10-mb relative humidity (CMIP)', '%', missing_value=missing_value )
       idiag%id_rh50_cmip = register_diag_field ( trim(field), 'rh50_cmip', axes(1:2), Time,       &
                           '50-mb relative humidity (CMIP)', '%', missing_value=missing_value )
       idiag%id_rh100_cmip = register_diag_field ( trim(field), 'rh100_cmip', axes(1:2), Time,       &
                           '100-mb relative humidity (CMIP)', '%', missing_value=missing_value )
       idiag%id_rh250_cmip = register_diag_field ( trim(field), 'rh250_cmip', axes(1:2), Time,       &
                           '250-mb relative humidity (CMIP)', '%', missing_value=missing_value )
       idiag%id_rh300_cmip = register_diag_field ( trim(field), 'rh300_cmip', axes(1:2), Time,       &
                           '300-mb relative humidity (CMIP)', '%', missing_value=missing_value )
       idiag%id_rh500_cmip = register_diag_field ( trim(field), 'rh500_cmip', axes(1:2), Time,       &
                           '500-mb relative humidity (CMIP)', '%', missing_value=missing_value )
       idiag%id_rh700_cmip = register_diag_field ( trim(field), 'rh700_cmip', axes(1:2), Time,       &
                           '700-mb relative humidity (CMIP)', '%', missing_value=missing_value )
       idiag%id_rh850_cmip = register_diag_field ( trim(field), 'rh850_cmip', axes(1:2), Time,       &
                           '850-mb relative humidity (CMIP)', '%', missing_value=missing_value )
       idiag%id_rh1000_cmip = register_diag_field ( trim(field), 'rh1000_cmip', axes(1:2), Time,       &
                           '1000-mb relative humidity (CMIP)', '%', missing_value=missing_value )
!--------------------------
! specific humidity:
!--------------------------
       idiag%id_omg10 = register_diag_field ( trim(field), 'omg10', axes(1:2), Time,       &
                           '10-mb omega', 'Pa/s', missing_value=missing_value )
       idiag%id_omg50 = register_diag_field ( trim(field), 'omg50', axes(1:2), Time,       &
                           '50-mb omega', 'Pa/s', missing_value=missing_value )
       idiag%id_omg100 = register_diag_field ( trim(field), 'omg100', axes(1:2), Time,       &
                           '100-mb omega', 'Pa/s', missing_value=missing_value )
       idiag%id_omg200 = register_diag_field ( trim(field), 'omg200', axes(1:2), Time,       &
                           '200-mb omega', 'Pa/s', missing_value=missing_value )
       idiag%id_omg250 = register_diag_field ( trim(field), 'omg250', axes(1:2), Time,       &
                           '250-mb omega', 'Pa/s', missing_value=missing_value )
       idiag%id_omg300 = register_diag_field ( trim(field), 'omg300', axes(1:2), Time,       &
                           '300-mb omega', 'Pa/s', missing_value=missing_value )
       idiag%id_omg500 = register_diag_field ( trim(field), 'omg500', axes(1:2), Time,       &
                           '500-mb omega', 'Pa/s', missing_value=missing_value )
       idiag%id_omg700 = register_diag_field ( trim(field), 'omg700', axes(1:2), Time,       &
                           '700-mb omega', 'Pa/s', missing_value=missing_value )
       idiag%id_omg850 = register_diag_field ( trim(field), 'omg850', axes(1:2), Time,       &
                           '850-mb omega', 'Pa/s', missing_value=missing_value )
       idiag%id_omg1000 = register_diag_field ( trim(field), 'omg1000', axes(1:2), Time,       &
                           '1000-mb omega', 'Pa/s', missing_value=missing_value )

       do i=1, ncnst
!--------------------
! Tracer diagnostics:
!--------------------
           call get_tracer_names ( MODEL_ATMOS, i, tname, tlongname, tunits )
           idiag%id_tracer(i) = register_diag_field ( field, trim(tname),  &
                axes(1:3), Time, trim(tlongname), &
                trim(tunits), missing_value=missing_value)
           if (master) then
               if (idiag%id_tracer(i) > 0) then
                   unit = stdlog()
                   write(unit,'(a,a,a,a)') &
                        & 'Diagnostics available for tracer ',trim(tname), &
                        ' in module ', trim(field)
               end if
           endif
!----------------------------------
! ESM Tracer dmmr/dvmr diagnostics:
!   for specific elements only
!----------------------------------
!---co2
           if (trim(tname).eq.'co2') then
               idiag%w_mr(:) = WTMCO2
               idiag%id_tracer_dmmr(i) = register_diag_field ( field, trim(tname)//'_dmmr',  &
                    axes(1:3), Time, trim(tlongname)//" (dry mmr)",           &
                    trim(tunits), missing_value=missing_value)
               idiag%id_tracer_dvmr(i) = register_diag_field ( field, trim(tname)//'_dvmr',  &
                    axes(1:3), Time, trim(tlongname)//" (dry vmr)",           &
                    'mol/mol', missing_value=missing_value)
               if (master) then
                   unit = stdlog()
                   if (idiag%id_tracer_dmmr(i) > 0) then
                       write(unit,'(a,a,a,a)') 'Diagnostics available for '//trim(tname)//' dry mmr ', &
                              trim(tname)//'_dmmr', ' in module ', trim(field)
                   end if
                   if (idiag%id_tracer_dvmr(i) > 0) then
                       write(unit,'(a,a,a,a)') 'Diagnostics available for '//trim(tname)//' dry vmr ', &
                            trim(tname)//'_dvmr', ' in module ', trim(field)
                   end if
               endif
           endif
!---end co2

       enddo

       if ( Atm(1)%flagstruct%consv_am .or. idiag%id_mq > 0 .or. idiag%id_amdt > 0 )  then
            allocate ( idiag%zxg(isc:iec,jsc:jec) )
! Initialize gradient of terrain for mountain torque computation:
            call init_mq(Atm(n)%phis, Atm(n)%gridstruct, &
                 npx, npy, isc, iec, jsc, jec, Atm(n)%ng)
       endif

!    end do

    call nullify_domain()  ! Nullify  set_domain info

    module_is_initialized=.true.
    istep = 0
 end subroutine fv_diag_init


 subroutine init_mq(phis, gridstruct, npx, npy, is, ie, js, je, ng)
    integer, intent(in):: npx, npy, is, ie, js, je, ng
    real, intent(in):: phis(is-ng:ie+ng, js-ng:je+ng)
    type(fv_grid_type), intent(IN), target :: gridstruct

! local:
    real zs(is-ng:ie+ng, js-ng:je+ng)
    real zb(is-ng:ie+ng, js-ng:je+ng)
    real pdx(3,is:ie,js:je+1)
    real pdy(3,is:ie+1,js:je)
    integer i, j, n

    real, pointer :: rarea(:,:)
    real, pointer, dimension(:,:) :: dx, dy
    real, pointer, dimension(:,:,:) :: en1, en2, agrid, vlon, vlat

    rarea => gridstruct%rarea
    dx    => gridstruct%dx
    dy    => gridstruct%dy
    en1   => gridstruct%en1
    en2   => gridstruct%en2
    agrid => gridstruct%agrid
    vlon  => gridstruct%vlon
    vlat  => gridstruct%vlat

!   do j=js,je
!      do i=is,ie
    do j=js-ng,je+ng
       do i=is-ng,ie+ng
          zs(i,j) = phis(i,j) / grav
       enddo
    enddo
!   call mpp_update_domains( zs, domain )

!   call a2b_ord2(zs, zb, gridstruct, npx, npy, is, ie, js, je, ng)
    call a2b_ord4(zs, zb, gridstruct, npx, npy, is, ie, js, je, ng)

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

! Compute "volume-mean" gradient by Green's theorem
    do j=js,je
       do i=is,ie
          idiag%zxg(i,j) = vlon(1,i,j)*(pdx(1,i,j+1)-pdx(1,i,j)-pdy(1,i,j)+pdy(1,i+1,j))  &
                         + vlon(2,i,j)*(pdx(2,i,j+1)-pdx(2,i,j)-pdy(2,i,j)+pdy(2,i+1,j))  &
                         + vlon(3,i,j)*(pdx(3,i,j+1)-pdx(3,i,j)-pdy(3,i,j)+pdy(3,i+1,j))
! dF/d(lamda) = radius*cos(agrid(i,j,2)) * dF/dx, F is a scalar
!                                                       ________________________
          idiag%zxg(i,j) =  idiag%zxg(i,j)*rarea(i,j) * radius*cos(agrid(i,j,2))
!                                                       ^^^^^^^^^^^^^^^^^^^^^^^^
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
    real, allocatable :: u2(:,:), v2(:,:), x850(:,:)
    real, allocatable :: dmmr(:,:,:), dvmr(:,:,:)
    real height(2)
    real plevs(10)
    real    :: tot_mq, tmp, sar, slon, slat
    logical :: used
    logical :: bad_range
    integer i,j,k, yr, mon, dd, hr, mn, days, seconds, nq
    character(len=128)   :: tname
    real, parameter:: ws_0 = 16.   ! minimum max_wind_speed within the 7x7 search box
    real, parameter:: ws_1 = 20.
    real, parameter:: vort_c0= 2.2e-5 
    logical, allocatable :: storm(:,:), cat_crt(:,:)

#ifdef MARS_GCM
    real  ::   atm_mass,  sfc_mass, atm_cloud
    real  ::   tsfc_dust, tcol_dust
#endif

    real :: tmp2, pvsum

    !!! CLEANUP: does it really make sense to have this routine loop over Atm% anymore? We assume n=1 below anyway

! cat15: SLP<1000; srf_wnd>ws_0; vort>vort_c0
! cat25: SLP< 980; srf_wnd>ws_1; vort>vort_c0
! cat35: SLP< 964; srf_wnd>ws_1; vort>vort_c0
! cat45: SLP< 944; srf_wnd>ws_1; vort>vort_c0

    height(1) = 5.E3      ! for computing 5-km "pressure"
    height(2) = 0.        ! for sea-level pressure

    ntileMe = size(Atm(:))
    n = 1
    isc = Atm(n)%bd%isc; iec = Atm(n)%bd%iec
    jsc = Atm(n)%bd%jsc; jec = Atm(n)%bd%jec
    ngc = Atm(n)%ng
    npz = Atm(n)%npz
    ptop = Atm(n)%ak(1)
    nq = size (Atm(n)%q,4)

    isd = Atm(n)%bd%isd; ied = Atm(n)%bd%ied
    jsd = Atm(n)%bd%jsd; jed = Atm(n)%bd%jed


    if( idiag%id_c15>0 ) then
        allocate (   storm(isc:iec,jsc:jec) )
        allocate ( depress(isc:iec,jsc:jec) )
        allocate (  ws_max(isc:iec,jsc:jec) )
        allocate ( cat_crt(isc:iec,jsc:jec) )
        allocate (tc_count(isc:iec,jsc:jec) )
    endif

    if( idiag%id_x850>0 ) then
        allocate ( x850(isc:iec,jsc:jec) )
    endif

    fv_time = Time
    call set_domain(Atm(1)%domain)

    if ( m_calendar ) then
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
                 istep = istep + 1
                 prt_minmax = mod(istep, -print_freq) == 0
         else
                 prt_minmax = mod(hr, print_freq) == 0 .and. mn==0 .and. seconds==0
         endif
     else
         call get_time (fv_time, seconds,  days)
         if( print_freq == 0 ) then
                 prt_minmax = .false.
         elseif( print_freq < 0 ) then
                 istep = istep + 1
                 prt_minmax = mod(istep, -print_freq) == 0
         else
                 prt_minmax = mod(seconds, 3600*print_freq) == 0
         endif
     endif

     if(prt_minmax) then
#if defined(MARS_GCM) || defined(VENUS_GCM)
        if(master) write(*,*) Days, seconds
#else
         if ( m_calendar ) then
              if(master) write(*,*) yr, mon, dd, hr, mn, seconds
         else
              if(master) write(*,*) Days, seconds
         endif
#endif
     endif

    if( prt_minmax ) then

        call prt_mxm('ZS', idiag%zsurf,     isc, iec, jsc, jec, 0,   1, 1.0, Atm(n)%gridstruct%area, Atm(n)%domain)
        call prt_maxmin('PS', Atm(n)%ps, isc, iec, jsc, jec, ngc, 1, 0.01)

#ifdef TEST_TRACER
        call prt_mass(npz, nq, isc, iec, jsc, jec, ngc, max(1,Atm(n)%flagstruct%nwat),    &
                      Atm(n)%ps, Atm(n)%delp, Atm(n)%q, Atm(n)%gridstruct%area, Atm(n)%domain)
#else
        call prt_mass(npz, nq, isc, iec, jsc, jec, ngc, Atm(n)%flagstruct%nwat,    &
                      Atm(n)%ps, Atm(n)%delp, Atm(n)%q, Atm(n)%gridstruct%area, Atm(n)%domain)
#endif

#ifndef SW_DYNAMICS
             idiag%steps = idiag%steps + 1
           idiag%efx_sum = idiag%efx_sum + E_Flux
        if ( idiag%steps <= max_step ) idiag%efx(idiag%steps) = E_Flux
        if (master)  then
            write(*,*) 'ENG Deficit (W/m**2)', trim(gn), '=', E_Flux
        endif

           idiag%efx_sum_nest = idiag%efx_sum_nest + E_Flux_nest
        if ( idiag%steps <= max_step ) idiag%efx_nest(idiag%steps) = E_Flux_nest
        if (master .and. abs(E_Flux_nest) > 0.)  then
            write(*,*) 'ENG Deficit from two-way update (W/m**2)', trim(gn), '=', E_Flux_nest
        endif
#endif
        call prt_maxmin('UA', Atm(n)%ua, isc, iec, jsc, jec, ngc, npz, 1.)
        call prt_maxmin('VA', Atm(n)%va, isc, iec, jsc, jec, ngc, npz, 1.)

        if ( .not. Atm(n)%flagstruct%hydrostatic ) then
          call prt_maxmin('W ', Atm(n)%w , isc, iec, jsc, jec, ngc, npz, 1.)
          if ( Atm(n)%flagstruct%hybrid_z ) call prt_maxmin('Hybrid_ZTOP (km)', Atm(n)%ze0(isc:iec,jsc:jec,1), &
                                                 isc, iec, jsc, jec, 0, 1, 1.E-3)
          call prt_maxmin('DZ (m)', Atm(n)%delz(isc:iec,jsc:jec,1:npz),    &
                          isc, iec, jsc, jec, 0, npz, 1.)
          call prt_maxmin('Bottom DZ (m)', Atm(n)%delz(isc:iec,jsc:jec,npz),    &
                          isc, iec, jsc, jec, 0, 1, 1.)
!         call prt_maxmin('Top DZ (m)', Atm(n)%delz(isc:iec,jsc:jec,1),    &
!                         isc, iec, jsc, jec, 0, 1, 1.)
        endif

#ifndef SW_DYNAMICS
        call prt_maxmin('TA', Atm(n)%pt,   isc, iec, jsc, jec, ngc, npz, 1.)
        call prt_maxmin('OM', Atm(n)%omga, isc, iec, jsc, jec, ngc, npz, 1.)
#endif

#if defined(MARS_GCM) && defined(MARS_SURFACE)
        atm_mass  = g_sum( Atm(n)%domain, Atm(n)%ps(isc:iec,jsc:jec), isc, iec, jsc, jec, ngc, area,mode=1)
        sfc_mass  = g_sum( Atm(n)%domain, sfc_snow,isc, iec, jsc, jec, ngc, area,mode=1)
        sfc_mass= sfc_mass*grav   !   Conversion to pressure units

        if(master) write(*,*) 'Atmospheric CO2 (mb)', trim(gn), ' =', atm_mass*0.01
        if(master) write(*,*) 'CO2 sfc frost   (mb)', trim(gn), ' =', sfc_mass*0.01
        if(master) write(*,*) 'Total CO2 Inventory ', trim(gn), ' =', (atm_mass+sfc_mass)*0.01

#ifdef WATER_CYCLE
        sfc_mass  = g_sum( Atm(n)%domain, sfc_frost, isc, iec, jsc, jec, ngc, area,mode=1)
        atm_mass  = g_sum( Atm(n)%domain,      wcol, isc, iec, jsc, jec, ngc, area,mode=1)
        atm_cloud = g_sum( Atm(n)%domain,   cldcol, isc, iec, jsc, jec, ngc, area,mode=1)
        sfc_mass= sfc_mass - 3.7   !  Arbitrary offset

        if(master) write(*,*) 'Atmospheric H2o vapor (kg/m**2)', trim(gn), ' =', atm_mass
        if(master) write(*,*) 'Atmospheric H2o cloud (kg/m**2)', trim(gn), ' =', atm_cloud
        if(master) write(*,*) 'Total Atmospheric H2o', trim(gn), ' ', atm_cloud + atm_mass

        if(master) write(*,*) 'H2O surface frost (kg/m**2)', trim(gn), ' ==', sfc_mass
        if(master) write(*,*) 'Total H2O inventory', trim(gn), ' =', atm_mass+sfc_mass+atm_cloud
#endif WATER_CYCLE

#ifdef DUST_SOURCE
        tsfc_dust  = g_sum( Atm(n)%domain, sfc_dust(:,:,1),isc, iec, jsc, jec, ngc, area,mode=1)
        tcol_dust  = g_sum( Atm(n)%domain, odcol   (:,:,1),isc, iec, jsc, jec, ngc, area,mode=1)

        if(master) write(*,*) 'Surface dust inventory (kg/m**2)', trim(gn), ' =', tsfc_dust - 30.0
        if(master) write(*,*) 'Atmospheric dust (kg/m**2)', trim(gn), ' =', tcol_dust
        if(master) write(*,*) 'Total dust inventory', trim(gn), ' ', tsfc_dust - 30.0 + tcol_dust
#endif DUST_SOURCE
#endif

    elseif ( Atm(n)%flagstruct%range_warn ) then
         call range_check('DELP', Atm(n)%delp, isc, iec, jsc, jec, ngc, npz, Atm(n)%gridstruct%agrid,    &
                           0.01*ptop, 200.E2, bad_range)
         call range_check('UA', Atm(n)%ua, isc, iec, jsc, jec, ngc, npz, Atm(n)%gridstruct%agrid,   &
                           -220., 250., bad_range)
         call range_check('VA', Atm(n)%ua, isc, iec, jsc, jec, ngc, npz, Atm(n)%gridstruct%agrid,   &
                           -220., 220., bad_range)
#ifndef SW_DYNAMICS
         call range_check('TA', Atm(n)%pt, isc, iec, jsc, jec, ngc, npz, Atm(n)%gridstruct%agrid,   &
                           150., 350., bad_range)
#endif

    endif

    allocate ( a2(isc:iec,jsc:jec) )
    allocate ( u2(isc:iec,jsc:jec) )
    allocate ( v2(isc:iec,jsc:jec) )
    allocate ( wk(isc:iec,jsc:jec,npz) )
    if ( any(idiag%id_tracer_dmmr > 0) .or. any(idiag%id_tracer_dvmr > 0) ) then
        allocate ( dmmr(isc:iec,jsc:jec,1:npz) )
        allocate ( dvmr(isc:iec,jsc:jec,1:npz) )
    endif

!    do n = 1, ntileMe
    n = 1

#ifdef DYNAMICS_ZS
       if(idiag%id_zsurf > 0)  used=send_data(idiag%id_zsurf, idiag%zsurf, Time)
#endif
       if(idiag%id_ps > 0) used=send_data(idiag%id_ps, Atm(n)%ps(isc:iec,jsc:jec), Time)

       if(idiag%id_c15>0 .or. idiag%id_c25>0 .or. idiag%id_c35>0 .or. idiag%id_c45>0) then
          call wind_max(isc, iec, jsc, jec ,isd, ied, jsd, jed, Atm(n)%ua(isc:iec,jsc:jec,npz),   &
                        Atm(n)%va(isc:iec,jsc:jec,npz), ws_max, Atm(n)%domain)
          do j=jsc,jec
             do i=isc,iec
                if( abs(Atm(n)%gridstruct%agrid(i,j,2)*rad2deg)<45.0 .and.     &
                    Atm(n)%phis(i,j)*ginv<500.0 .and. ws_max(i,j)>ws_0 ) then
                    storm(i,j) = .true.
                else
                    storm(i,j) = .false.
                endif
             enddo
          enddo
       endif

       if ( idiag%id_vort850>0 .or. idiag%id_vorts>0 .or. idiag%id_vort>0 .or. idiag%id_pv>0 .or. idiag%id_rh>0 .or. idiag%id_x850>0 ) then
          call get_vorticity(isc, iec, jsc, jec, isd, ied, jsd, jed, npz, Atm(n)%u, Atm(n)%v, wk, &
               Atm(n)%gridstruct%dx, Atm(n)%gridstruct%dy, Atm(n)%gridstruct%rarea)

          if(idiag%id_vort >0) used=send_data(idiag%id_vort,  wk, Time)
          if(idiag%id_vorts>0) used=send_data(idiag%id_vorts, wk(isc:iec,jsc:jec,npz), Time)

          if(idiag%id_c15>0) then
             do j=jsc,jec
                do i=isc,iec
                   if ( storm(i,j) )    &
                   storm(i,j) = (Atm(n)%gridstruct%agrid(i,j,2)>0. .and. wk(i,j,npz)> vort_c0) .or. &
                                (Atm(n)%gridstruct%agrid(i,j,2)<0. .and. wk(i,j,npz)<-vort_c0) 
                enddo
             enddo
          endif


          if(idiag%id_vort850>0 .or. idiag%id_c15>0 .or. idiag%id_x850>0) then
             call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                       850.e2, Atm(n)%peln, wk, a2)
             used=send_data(idiag%id_vort850, a2, Time)
             if ( idiag%id_x850>0 ) x850(:,:) = a2(:,:) 

             if(idiag%id_c15>0) then
             do j=jsc,jec
                do i=isc,iec
                   if ( storm(i,j) )    &
                     storm(i,j) = (Atm(n)%gridstruct%agrid(i,j,2)>0. .and. a2(i,j)> vort_c0) .or.     &
                                  (Atm(n)%gridstruct%agrid(i,j,2)<0. .and. a2(i,j)<-vort_c0) 
                enddo
             enddo
             endif

          endif

          if ( idiag%id_srh > 0 ) then
             call helicity_relative(isc, iec, jsc, jec, ngc, npz, zvir, sphum, a2, &
                  Atm(n)%ua, Atm(n)%va, Atm(n)%delz, Atm(n)%q,   &
                  Atm(n)%flagstruct%hydrostatic, Atm(n)%pt, Atm(n)%peln, Atm(n)%phis, grav)
             used = send_data ( idiag%id_srh, a2, Time )
             if(prt_minmax) then
                do j=jsc,jec
                   do i=isc,iec
                      tmp = rad2deg * Atm(n)%gridstruct%agrid(i,j,1)
                      tmp2 = rad2deg * Atm(n)%gridstruct%agrid(i,j,2)
                      if (  tmp2<25. .or. tmp2>50.    &
                           .or. tmp<235. .or. tmp>300. ) then
                         a2(i,j) = 0.
                      endif
                   enddo
                enddo
                call prt_maxmin('SRH over CONUS', a2, isc, iec, jsc, jec, 0,   1, 1.)
             endif
          endif

          if ( idiag%id_pv > 0 ) then
! Note: this is expensive computation.
              call pv_entropy(isc, iec, jsc, jec, ngc, npz, wk,    &
                              Atm(n)%gridstruct%f0, Atm(n)%pt, Atm(n)%pkz, Atm(n)%delp, grav)
              used = send_data ( idiag%id_pv, wk, Time )
          endif

! Relative Humidity
          if ( idiag%id_rh > 0 ) then
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
               used = send_data ( idiag%id_rh, wk, Time )
               if(prt_minmax) then
                  call prt_maxmin('RH_sf (%)', wk(isc:iec,jsc:jec,npz), isc, iec, jsc, jec, 0,   1, 1.)
                  call prt_maxmin('RH_3D (%)', wk, isc, iec, jsc, jec, 0, npz, 1.)
               endif
          endif

       endif

       ! rel hum from physics at selected press levels (for IPCC)
       if (idiag%id_rh50>0  .or. idiag%id_rh100>0 .or. idiag%id_rh200>0 .or. idiag%id_rh250>0 .or. &
           idiag%id_rh500>0 .or. idiag%id_rh700>0 .or. idiag%id_rh850>0 .or. idiag%id_rh1000>0) then
           ! compute mean pressure
           do k=1,npz
               do j=jsc,jec
               do i=isc,iec
                   a2(i,j) = Atm(n)%delp(i,j,k)/(Atm(n)%peln(i,k+1,j)-Atm(n)%peln(i,k,j))
               enddo
               enddo
               call rh_calc (a2, Atm(n)%pt(isc:iec,jsc:jec,k), &
                             Atm(n)%q(isc:iec,jsc:jec,k,sphum), wk(isc:iec,jsc:jec,k))
           enddo
           if (idiag%id_rh50>0) then
               call interpolate_vertical(isc, iec, jsc, jec, npz, 50.e2, Atm(n)%peln, wk(isc:iec,jsc:jec,:), a2)
               used=send_data(idiag%id_rh50, a2, Time)
           endif
           if (idiag%id_rh100>0) then
               call interpolate_vertical(isc, iec, jsc, jec, npz, 100.e2, Atm(n)%peln, wk(isc:iec,jsc:jec,:), a2)
               used=send_data(idiag%id_rh100, a2, Time)
           endif
           if (idiag%id_rh200>0) then
               call interpolate_vertical(isc, iec, jsc, jec, npz, 200.e2, Atm(n)%peln, wk(isc:iec,jsc:jec,:), a2)
               used=send_data(idiag%id_rh200, a2, Time)
           endif
           if (idiag%id_rh250>0) then
               call interpolate_vertical(isc, iec, jsc, jec, npz, 250.e2, Atm(n)%peln, wk(isc:iec,jsc:jec,:), a2)
               used=send_data(idiag%id_rh250, a2, Time)
           endif
           if (idiag%id_rh500>0) then
               call interpolate_vertical(isc, iec, jsc, jec, npz, 500.e2, Atm(n)%peln, wk(isc:iec,jsc:jec,:), a2)
               used=send_data(idiag%id_rh500, a2, Time)
           endif
           if (idiag%id_rh700>0) then
               call interpolate_vertical(isc, iec, jsc, jec, npz, 700.e2, Atm(n)%peln, wk(isc:iec,jsc:jec,:), a2)
               used=send_data(idiag%id_rh700, a2, Time)
           endif
           if (idiag%id_rh850>0) then
               call interpolate_vertical(isc, iec, jsc, jec, npz, 850.e2, Atm(n)%peln, wk(isc:iec,jsc:jec,:), a2)
               used=send_data(idiag%id_rh850, a2, Time)
           endif
           if (idiag%id_rh1000>0) then
               call interpolate_vertical(isc, iec, jsc, jec, npz, 1000.e2, Atm(n)%peln, wk(isc:iec,jsc:jec,:), a2)
               used=send_data(idiag%id_rh1000, a2, Time)
           endif
       endif
       ! rel hum (CMIP definition) at selected press levels  (for IPCC)
       if (idiag%id_rh10_cmip>0 .or. idiag%id_rh50_cmip>0  .or. idiag%id_rh100_cmip>0 .or. &
           idiag%id_rh250_cmip>0 .or. idiag%id_rh500_cmip>0 .or. idiag%id_rh700_cmip>0 .or. &
           idiag%id_rh850_cmip>0 .or. idiag%id_rh1000_cmip>0) then
           ! compute mean pressure
           do k=1,npz
               do j=jsc,jec
               do i=isc,iec
                   a2(i,j) = Atm(n)%delp(i,j,k)/(Atm(n)%peln(i,k+1,j)-Atm(n)%peln(i,k,j))
               enddo
               enddo
               call rh_calc (a2, Atm(n)%pt(isc:iec,jsc:jec,k), &
                             Atm(n)%q(isc:iec,jsc:jec,k,sphum), wk(isc:iec,jsc:jec,k), do_cmip=.true.)
           enddo
           if (idiag%id_rh10_cmip>0) then
               call interpolate_vertical(isc, iec, jsc, jec, npz, 10.e2, Atm(n)%peln, wk(isc:iec,jsc:jec,:), a2)
               used=send_data(idiag%id_rh10_cmip, a2, Time)
           endif
           if (idiag%id_rh50_cmip>0) then
               call interpolate_vertical(isc, iec, jsc, jec, npz, 50.e2, Atm(n)%peln, wk(isc:iec,jsc:jec,:), a2)
               used=send_data(idiag%id_rh50_cmip, a2, Time)
           endif
           if (idiag%id_rh100_cmip>0) then
               call interpolate_vertical(isc, iec, jsc, jec, npz, 100.e2, Atm(n)%peln, wk(isc:iec,jsc:jec,:), a2)
               used=send_data(idiag%id_rh100_cmip, a2, Time)
           endif
           if (idiag%id_rh250_cmip>0) then
               call interpolate_vertical(isc, iec, jsc, jec, npz, 250.e2, Atm(n)%peln, wk(isc:iec,jsc:jec,:), a2)
               used=send_data(idiag%id_rh250_cmip, a2, Time)
           endif
           if (idiag%id_rh500_cmip>0) then
               call interpolate_vertical(isc, iec, jsc, jec, npz, 500.e2, Atm(n)%peln, wk(isc:iec,jsc:jec,:), a2)
               used=send_data(idiag%id_rh500_cmip, a2, Time)
           endif
           if (idiag%id_rh700_cmip>0) then
               call interpolate_vertical(isc, iec, jsc, jec, npz, 700.e2, Atm(n)%peln, wk(isc:iec,jsc:jec,:), a2)
               used=send_data(idiag%id_rh700_cmip, a2, Time)
           endif
           if (idiag%id_rh850_cmip>0) then
               call interpolate_vertical(isc, iec, jsc, jec, npz, 850.e2, Atm(n)%peln, wk(isc:iec,jsc:jec,:), a2)
               used=send_data(idiag%id_rh850_cmip, a2, Time)
           endif
           if (idiag%id_rh1000_cmip>0) then
               call interpolate_vertical(isc, iec, jsc, jec, npz, 1000.e2, Atm(n)%peln, wk(isc:iec,jsc:jec,:), a2)
               used=send_data(idiag%id_rh1000_cmip, a2, Time)
           endif
       endif

       if(idiag%id_c25>0 .or. idiag%id_c35>0 .or. idiag%id_c45>0) then
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



       if( idiag%id_slp>0 .or. idiag%id_tm>0 .or. idiag%id_hght>0 .or. idiag%id_c15>0 ) then

          allocate ( wz(isc:iec,jsc:jec,npz+1) )
          call get_height_field(isc, iec, jsc, jec, ngc, npz, wz, Atm(n)%pt, Atm(n)%q, Atm(n)%peln, zvir)
          if( prt_minmax )   &
          call prt_maxmin('ZTOP', wz(isc:iec,jsc:jec,1), isc, iec, jsc, jec, 0, 1, 1.E-3)

          if(idiag%id_slp > 0) then
! Cumpute SLP (pressure at height=0)
          allocate ( slp(isc:iec,jsc:jec) )
          call get_pressure_given_height(isc, iec, jsc, jec, ngc, npz, wz, 1, height(2),   &
                                        Atm(n)%pt(:,:,npz), Atm(n)%peln, slp, 0.01)
          used = send_data (idiag%id_slp, slp, Time)
          if( prt_minmax )   &
             call prt_maxmin('SLP', slp, isc, iec, jsc, jec, 0, 1, 1.)
          endif

! Compute H3000 and/or H500
          if( idiag%id_tm>0 .or. idiag%id_hght>0 .or. idiag%id_ppt>0) then

              allocate( a3(isc:iec,jsc:jec,size(plevs,1)) )
              plevs(1)  = log( 1000. )
              plevs(2)  = log( 5000. )
              plevs(3)  = log( 10000. )
              plevs(4)  = log( 20000. )
              plevs(5)  = log( 25000. )
              plevs(6)  = log( 30000. )
              plevs(7)  = log( 50000. )
              plevs(8)  = log( 70000. )
              plevs(9)  = log( 85000. )
              plevs(10) = log(100000. )

             call get_height_given_pressure(isc, iec, jsc, jec, ngc, npz, wz, 10, plevs, Atm(n)%peln, a3)
             if(idiag%id_h10>0)   used = send_data ( idiag%id_h10,   a3(isc:iec,jsc:jec,1),  Time )
             if(idiag%id_h50>0)   used = send_data ( idiag%id_h50,   a3(isc:iec,jsc:jec,2),  Time )
             if(idiag%id_h100>0)  used = send_data ( idiag%id_h100,  a3(isc:iec,jsc:jec,3),  Time )
             if(idiag%id_h200>0)  used = send_data ( idiag%id_h200,  a3(isc:iec,jsc:jec,4),  Time )
             if(idiag%id_h250>0)  used = send_data ( idiag%id_h250,  a3(isc:iec,jsc:jec,5),  Time )
             if(idiag%id_h300>0)  used = send_data ( idiag%id_h300,  a3(isc:iec,jsc:jec,6),  Time )
             if(idiag%id_h500>0)  used = send_data ( idiag%id_h500,  a3(isc:iec,jsc:jec,7),  Time )
             if(idiag%id_h700>0)  used = send_data ( idiag%id_h700,  a3(isc:iec,jsc:jec,8),  Time )
             if(idiag%id_h850>0)  used = send_data ( idiag%id_h850,  a3(isc:iec,jsc:jec,9),  Time )
             if(idiag%id_h1000>0) used = send_data ( idiag%id_h1000, a3(isc:iec,jsc:jec,10), Time )

             ! mean temp 300mb to 500mb
             if( idiag%id_tm>0 ) then
                 do j=jsc,jec
                    do i=isc,iec
                       a2(i,j) = grav*(a3(i,j,6)-a3(i,j,7))/(rdgas*(plevs(7)-plevs(6)))
                    enddo
                 enddo
                 used = send_data ( idiag%id_tm, a2, Time )
             endif

            if(idiag%id_c15>0 .or. idiag%id_c25>0 .or. idiag%id_c35>0 .or. idiag%id_c45>0) then
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
             used = send_data(idiag%id_c15, depress, Time)
             if(idiag%id_f15>0) used = send_data(idiag%id_f15, tc_count, Time)
             if(prt_minmax) then
                do j=jsc,jec
                   do i=isc,iec
                      slon = rad2deg*Atm(n)%gridstruct%agrid(i,j,1)
                      slat = rad2deg*Atm(n)%gridstruct%agrid(i,j,2)
! Western Pac: negative; positive elsewhere
                      if ( slat>0. .and. slat<40. .and. slon>110. .and. slon<180. ) then
                           depress(i,j) = -depress(i,j)
                      endif
                   enddo
                enddo
                call prt_maxmin('Depress', depress, isc, iec, jsc, jec, 0,   1, 1.)
                do j=jsc,jec
                   do i=isc,iec
                      if ( Atm(n)%gridstruct%agrid(i,j,2)<0.) then
! Excluding the SH cyclones
                           depress(i,j) = 0.
                      endif
                   enddo
                enddo
                call prt_maxmin('NH Deps', depress, isc, iec, jsc, jec, 0,   1, 1.)

! ATL basin cyclones
                do j=jsc,jec
                   do i=isc,iec
                      tmp = rad2deg * Atm(n)%gridstruct%agrid(i,j,1)
                      if ( tmp<280. ) then
                           depress(i,j) = 0.
                      endif
                   enddo
                enddo
                call prt_maxmin('ATL Deps', depress, isc, iec, jsc, jec, 0,   1, 1.)
             endif
            endif

! Cat 2-5:
            if(idiag%id_c25>0) then
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
             used = send_data(idiag%id_c25, depress, Time)
             if(idiag%id_f25>0) used = send_data(idiag%id_f25, tc_count, Time)
            endif

! Cat 3-5:
            if(idiag%id_c35>0) then
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
             used = send_data(idiag%id_c35, depress, Time)
             if(idiag%id_f35>0) used = send_data(idiag%id_f35, tc_count, Time)
            endif

! Cat 4-5:
            if(idiag%id_c45>0) then
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
             used = send_data(idiag%id_c45, depress, Time)
             if(idiag%id_f45>0) used = send_data(idiag%id_f45, tc_count, Time)
            endif

            if (idiag%id_c15>0) then
                deallocate(depress)
                deallocate(cat_crt)
                deallocate(storm)
                deallocate(ws_max)
                deallocate(tc_count)
            endif

            if(idiag%id_slp>0 )  deallocate( slp )

            deallocate( a3 )
          endif

         deallocate ( wz )
       endif


       if(idiag%id_mq > 0)  then
          do j=jsc,jec
             do i=isc,iec
! zxg * surface pressure * 1.e-18--> Hadleys per unit area
! Unit Hadley = 1.E18 kg m**2 / s**2
                a2(i,j) = -1.e-18 * Atm(n)%ps(i,j)*idiag%zxg(i,j)
             enddo
          enddo
          used = send_data(idiag%id_mq, a2, Time)
          if( prt_minmax ) then
              tot_mq  = g_sum( Atm(n)%domain, a2, isc, iec, jsc, jec, ngc, Atm(n)%gridstruct%area, 0) 
              idiag%mtq_sum = idiag%mtq_sum + tot_mq
              if ( idiag%steps <= max_step ) idiag%mtq(idiag%steps) = tot_mq
              if(master) write(*,*) 'Total (global) mountain torque (Hadleys)=', tot_mq
          endif
       endif

       if (idiag%id_ts > 0) used = send_data(idiag%id_ts, Atm(n)%ts(isc:iec,jsc:jec), Time)

#ifdef MARS_GCM
       if ( idiag%id_t05>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      0.5e2, Atm(n)%peln, Atm(n)%pt(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_t05, a2, Time)
       endif
#  ifdef WATER_CYCLE
       if ( idiag%id_tq>0 ) then
          itrac= get_tracer_index (MODEL_ATMOS, 'h2o_vapor')
          a2 = 0.
          do k=1,npz
          do j=jsc,jec
             do i=isc,iec
                a2(i,j) = a2(i,j) + Atm(n)%q(i,j,k,itrac)*Atm(n)%delp(i,j,k)
             enddo
          enddo
          enddo
          used = send_data(idiag%id_tq, a2*ginv, Time)
       endif
#  endif WATER_CYCLE
#else
       if ( idiag%id_tq>0 ) then
          nwater = Atm(1)%flagstruct%nwat
          a2 = 0.
          do k=1,npz
          do j=jsc,jec
             do i=isc,iec
!                a2(i,j) = a2(i,j) + Atm(n)%q(i,j,k,1)*Atm(n)%delp(i,j,k)
                a2(i,j) = a2(i,j) + sum(Atm(n)%q(i,j,k,1:nwater))*Atm(n)%delp(i,j,k)
             enddo
          enddo
          enddo
          used = send_data(idiag%id_tq, a2*ginv, Time)
       endif
       if ( idiag%id_iw>0 ) then
          a2 = 0.
          do k=1,npz
          do j=jsc,jec
             do i=isc,iec
                a2(i,j) = a2(i,j) + Atm(n)%delp(i,j,k) * (     &
                                    Atm(n)%q(i,j,k,ice_wat) +  &
                                    Atm(n)%q(i,j,k,snowwat) +  &
                                    Atm(n)%q(i,j,k,graupel) )
             enddo
          enddo
          enddo
          used = send_data(idiag%id_iw, a2*ginv, Time)
       endif
       if ( idiag%id_lw>0 ) then
          a2 = 0.
          do k=1,npz
          do j=jsc,jec
             do i=isc,iec
                a2(i,j) = a2(i,j) + (Atm(n)%q(i,j,k,liq_wat)+Atm(n)%q(i,j,k,rainwat)) *    &
                                     Atm(n)%delp(i,j,k)
             enddo
          enddo
          enddo
          used = send_data(idiag%id_lw, a2*ginv, Time)
       endif
! Total 3D clouds
       if ( idiag%id_qn>0 ) then
          do k=1,npz
          do j=jsc,jec
             do i=isc,iec
                wk(i,j,k) = Atm(n)%q(i,j,k,liq_wat) + Atm(n)%q(i,j,k,ice_wat)
             enddo
          enddo
          enddo
          used = send_data(idiag%id_qn, wk, Time)
       endif
! Total 3D condensates
       if ( idiag%id_qp>0 ) then
          do k=1,npz
          do j=jsc,jec
             do i=isc,iec
                wk(i,j,k) = Atm(n)%q(i,j,k,snowwat) + Atm(n)%q(i,j,k,rainwat) +   &
                            Atm(n)%q(i,j,k,graupel)
             enddo
          enddo
          enddo
          used = send_data(idiag%id_qp, wk, Time)
       endif
#endif MARS_GCM

       if(idiag%id_us > 0 .and. idiag%id_vs > 0) then
          u2(:,:) = Atm(n)%ua(isc:iec,jsc:jec,npz)
          v2(:,:) = Atm(n)%va(isc:iec,jsc:jec,npz)
          do j=jsc,jec
             do i=isc,iec
                a2(i,j) = sqrt(u2(i,j)**2 + v2(i,j)**2)
             enddo
          enddo
          used=send_data(idiag%id_us, u2, Time)
          used=send_data(idiag%id_vs, v2, Time)
          if(prt_minmax) call prt_maxmin('Surf_wind_speed', a2, isc, iec, jsc, jec, 0, 1, 1.)
       endif
       if(idiag%id_tb > 0) then
          a2(:,:) = Atm(n)%pt(isc:iec,jsc:jec,npz)
          used=send_data(idiag%id_tb, a2, Time)
          if( prt_minmax )   &
          call prt_mxm('T_bot:', a2, isc, iec, jsc, jec, 0, 1, 1., Atm(n)%gridstruct%area, Atm(n)%domain)
       endif

       if(idiag%id_ua > 0) used=send_data(idiag%id_ua, Atm(n)%ua(isc:iec,jsc:jec,:), Time)
       if(idiag%id_va > 0) used=send_data(idiag%id_va, Atm(n)%va(isc:iec,jsc:jec,:), Time)

       if(idiag%id_ke > 0) then
          a2(:,:) = 0.
          do k=1,npz
          do j=jsc,jec
             do i=isc,iec
                a2(i,j) = a2(i,j) + Atm(n)%delp(i,j,k)*(Atm(n)%ua(i,j,k)**2+Atm(n)%va(i,j,k)**2)
             enddo
          enddo
          enddo
! Mass weighted KE 
          do j=jsc,jec
             do i=isc,iec
                a2(i,j) = 0.5*a2(i,j)/(Atm(n)%ps(i,j)-ptop)
             enddo
          enddo
          used=send_data(idiag%id_ke, a2, Time)
          if(prt_minmax) then
             tot_mq  = g_sum( Atm(n)%domain, a2, isc, iec, jsc, jec, ngc, Atm(n)%gridstruct%area, 1) 
             if (master) write(*,*) 'SQRT(2.*KE; m/s)=', sqrt(2.*tot_mq)  
          endif
       endif

       if(idiag%id_delp > 0) used=send_data(idiag%id_delp, Atm(n)%delp(isc:iec,jsc:jec,:), Time)
       if(idiag%id_delz > 0) then
          do k=1,npz
            do j=jsc,jec
            do i=isc,iec
               wk(i,j,k) = -Atm(n)%delz(i,j,k)
            enddo
            enddo
          enddo
          used=send_data(idiag%id_delz, wk, Time)
       endif

! pressure for masking p-level fields
! incorrectly defines a2 to be ps (in mb).
       if (idiag%id_pmask>0) then
            do j=jsc,jec
            do i=isc,iec
                a2(i,j) = exp((Atm(n)%peln(i,npz+1,j)+Atm(n)%peln(i,npz+1,j))*0.5)*0.01
               !a2(i,j) = Atm(n)%delp(i,j,k)/(Atm(n)%peln(i,k+1,j)-Atm(n)%peln(i,k,j))*0.01
            enddo
            enddo
            used=send_data(idiag%id_pmask, a2, Time)
       endif
! fix for pressure for masking p-level fields
! based on lowest-level pfull
! define pressure at lowest level the same as interpolate_vertical (in mb)
       if (idiag%id_pmaskv2>0) then
            do j=jsc,jec
            do i=isc,iec
                a2(i,j) = exp((Atm(n)%peln(i,npz,j)+Atm(n)%peln(i,npz+1,j))*0.5)*0.01
            enddo
            enddo
            used=send_data(idiag%id_pmaskv2, a2, Time)
       endif

! 10-mb
       if ( idiag%id_u10>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      10.e2, Atm(n)%peln, Atm(n)%ua(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_u10, a2, Time)
       endif
       if ( idiag%id_v10>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      10.e2, Atm(n)%peln, Atm(n)%va(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_v10, a2, Time)
       endif
       if ( idiag%id_t10>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      10.e2, Atm(n)%peln, Atm(n)%pt(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_t10, a2, Time)
       endif
       if ( idiag%id_q10>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      10.e2, Atm(n)%peln, Atm(n)%q(isc:iec,jsc:jec,:,sphum), a2)
            used=send_data(idiag%id_q10, a2, Time)
       endif
       if ( idiag%id_omg10>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      10.e2, Atm(n)%peln, Atm(n)%omga(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_omg10, a2, Time)
       endif
! 50-mb
       if ( idiag%id_u50>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      50.e2, Atm(n)%peln, Atm(n)%ua(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_u50, a2, Time)
       endif
       if ( idiag%id_v50>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      50.e2, Atm(n)%peln, Atm(n)%va(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_v50, a2, Time)
       endif
       if ( idiag%id_t50>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      50.e2, Atm(n)%peln, Atm(n)%pt(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_t50, a2, Time)
       endif
       if ( idiag%id_q50>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      50.e2, Atm(n)%peln, Atm(n)%q(isc:iec,jsc:jec,:,sphum), a2)
            used=send_data(idiag%id_q50, a2, Time)
       endif
       if ( idiag%id_omg50>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      50.e2, Atm(n)%peln, Atm(n)%omga(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_omg50, a2, Time)
       endif
! 100-mb
       if ( idiag%id_u100>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      100.e2, Atm(n)%peln, Atm(n)%ua(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_u100, a2, Time)
       endif
       if ( idiag%id_v100>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      100.e2, Atm(n)%peln, Atm(n)%va(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_v100, a2, Time)
       endif
       if ( idiag%id_t100>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      100.e2, Atm(n)%peln, Atm(n)%pt(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_t100, a2, Time)
            if( prt_minmax ) then
            call prt_mxm('T100:', a2, isc, iec, jsc, jec, 0, 1, 1., Atm(n)%gridstruct%area, Atm(n)%domain)
#ifdef WPAC
! Compute mean temp at 100 mb over the WPAC main convective region:
            tmp = 0.
            sar = 0.
            do j=jsc,jec
               do i=isc,iec
                  slon = Atm(n)%gridstruct%agrid(i,j,1)*rad2deg
                  slat = Atm(n)%gridstruct%agrid(i,j,2)*rad2deg
                  if( (slat>0.0 .and. slat<25.) .and. (slon>105.0 .and. slon<140.) ) then
                      sar = sar + Atm(n)%gridstruct%area(i,j)
                      tmp = tmp + a2(i,j)*Atm(n)%gridstruct%area(i,j)
                  endif
               enddo
            enddo
            call mp_reduce_sum(sar)
            if ( sar > 0. ) then
                 call mp_reduce_sum(tmp)
                 if (master) write(*,*) 'WPac: mean T100 =', tmp/sar
            endif
            endif
#else
! Compute mean temp at 100 mb over the [-10,10]
            tmp = 0.
            sar = 0.
            do j=jsc,jec
               do i=isc,iec
                  slat = Atm(n)%gridstruct%agrid(i,j,2)*rad2deg
                  if( (slat>-10.0 .and. slat<10.) ) then
                      sar = sar + Atm(n)%gridstruct%area(i,j)
                      tmp = tmp + a2(i,j)*Atm(n)%gridstruct%area(i,j)
                  endif
               enddo
            enddo
            call mp_reduce_sum(sar)
            if ( sar > 0. ) then
                 call mp_reduce_sum(tmp)
                 if (master) write(*,*) 'Tropical [10s,10n] mean T100 =', tmp/sar
            endif
            endif
#endif
       endif
       if ( idiag%id_q100>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      100.e2, Atm(n)%peln, Atm(n)%q(isc:iec,jsc:jec,:,sphum), a2)
            used=send_data(idiag%id_q100, a2, Time)
       endif
       if ( idiag%id_omg100>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      100.e2, Atm(n)%peln, Atm(n)%omga(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_omg100, a2, Time)
       endif
! 200-mb
       if ( idiag%id_u200>0 .or. idiag%id_s200>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      200.e2, Atm(n)%peln, Atm(n)%ua(isc:iec,jsc:jec,:), u2)
            if( idiag%id_u200>0 ) used=send_data(idiag%id_u200, u2, Time)
       endif
       if ( idiag%id_v200>0 .or. idiag%id_s200>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      200.e2, Atm(n)%peln, Atm(n)%va(isc:iec,jsc:jec,:), v2)
            if( idiag%id_v200>0 ) used=send_data(idiag%id_v200, v2, Time)
       endif
       if ( idiag%id_s200>0 ) then
            do j=jsc,jec
               do i=isc,iec
                  a2(i,j) = sqrt(u2(i,j)**2 + v2(i,j)**2)
               enddo
            enddo
            used=send_data(idiag%id_s200, a2, Time)
       endif
       if ( idiag%id_sl12>0 ) then   ! 13th level wind speed (~ 222 mb for the 32L setup)
            do j=jsc,jec
               do i=isc,iec
                  a2(i,j) = sqrt(Atm(n)%ua(i,j,12)**2 + Atm(n)%va(i,j,12)**2)
               enddo
            enddo
            used=send_data(idiag%id_sl12, a2, Time)
       endif
       if ( idiag%id_sl13>0 ) then   ! 13th level wind speed (~ 222 mb for the 32L setup)
            do j=jsc,jec
               do i=isc,iec
                  a2(i,j) = sqrt(Atm(n)%ua(i,j,13)**2 + Atm(n)%va(i,j,13)**2)
               enddo
            enddo
            used=send_data(idiag%id_sl13, a2, Time)
       endif

       if ( idiag%id_w200>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      200.e2, Atm(n)%peln, Atm(n)%w(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_w200, a2, Time)
       endif
       if ( idiag%id_t200>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      200.e2, Atm(n)%peln, Atm(n)%pt(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_t200, a2, Time)
            if( prt_minmax ) then
            call prt_mxm('T200:', a2, isc, iec, jsc, jec, 0, 1, 1., Atm(n)%gridstruct%area, Atm(n)%domain)
! Compute mean temp at 200 mb over the WPAC main convective region:
            tmp = 0.
            sar = 0.
            do j=jsc,jec
               do i=isc,iec
                  slon = Atm(n)%gridstruct%agrid(i,j,1)*rad2deg
                  slat = Atm(n)%gridstruct%agrid(i,j,2)*rad2deg
                  if( (slat>0.0 .and. slat<25.) .and. (slon>105.0 .and. slon<140.) ) then
                      sar = sar + Atm(n)%gridstruct%area(i,j)
                      tmp = tmp + a2(i,j)*Atm(n)%gridstruct%area(i,j)
                  endif
               enddo
            enddo
            call mp_reduce_sum(sar)
            if ( sar > 0. ) then
                 call mp_reduce_sum(tmp)
                 if (master) write(*,*) 'WPac: mean T200 =', tmp/sar
            endif
            endif
       endif
       if ( idiag%id_q200>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      200.e2, Atm(n)%peln, Atm(n)%q(isc:iec,jsc:jec,:,sphum), a2)
            used=send_data(idiag%id_q200, a2, Time)
       endif
       if ( idiag%id_omg200>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      200.e2, Atm(n)%peln, Atm(n)%omga(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_omg200, a2, Time)
       endif
! 250-mb
       if ( idiag%id_u250>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      250.e2, Atm(n)%peln, Atm(n)%ua(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_u250, a2, Time)
       endif
       if ( idiag%id_v250>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      250.e2, Atm(n)%peln, Atm(n)%va(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_v250, a2, Time)
       endif
       if ( idiag%id_t250>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      250.e2, Atm(n)%peln, Atm(n)%pt(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_t250, a2, Time)
       endif
       if ( idiag%id_q250>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      250.e2, Atm(n)%peln, Atm(n)%q(isc:iec,jsc:jec,:,sphum), a2)
            used=send_data(idiag%id_q250, a2, Time)
       endif
       if ( idiag%id_omg250>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      250.e2, Atm(n)%peln, Atm(n)%omga(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_omg250, a2, Time)
       endif
! 500-mb
       if ( idiag%id_u500>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      500.e2, Atm(n)%peln, Atm(n)%ua(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_u500, a2, Time)
       endif
       if ( idiag%id_v500>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      500.e2, Atm(n)%peln, Atm(n)%va(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_v500, a2, Time)
       endif
       if ( idiag%id_t500>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      500.e2, Atm(n)%peln, Atm(n)%pt(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_t500, a2, Time)
       endif
       if ( idiag%id_q500>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      500.e2, Atm(n)%peln, Atm(n)%q(isc:iec,jsc:jec,:,sphum), a2)
            used=send_data(idiag%id_q500, a2, Time)
       endif
       if ( idiag%id_omg500>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      500.e2, Atm(n)%peln, Atm(n)%omga(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_omg500, a2, Time)
       endif
! 700-mb
       if ( idiag%id_u700>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      700.e2, Atm(n)%peln, Atm(n)%ua(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_u700, a2, Time)
       endif
       if ( idiag%id_v700>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      700.e2, Atm(n)%peln, Atm(n)%va(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_v700, a2, Time)
       endif
       if ( idiag%id_t700>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      700.e2, Atm(n)%peln, Atm(n)%pt(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_t700, a2, Time)
       endif
       if ( idiag%id_q700>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      700.e2, Atm(n)%peln, Atm(n)%q(isc:iec,jsc:jec,:,sphum), a2)
            used=send_data(idiag%id_q700, a2, Time)
       endif
       if ( idiag%id_omg700>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      700.e2, Atm(n)%peln, Atm(n)%omga(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_omg700, a2, Time)
       endif
! 850-mb
       if ( idiag%id_u850>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      850.e2, Atm(n)%peln, Atm(n)%ua(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_u850, a2, Time)
       endif
       if ( idiag%id_v850>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      850.e2, Atm(n)%peln, Atm(n)%va(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_v850, a2, Time)
       endif
       if ( idiag%id_w850>0 .or. idiag%id_x850>0) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      850.e2, Atm(n)%peln, Atm(n)%w(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_w850, a2, Time)

            if ( idiag%id_x850>0 .and. idiag%id_vort850>0 ) then
                 x850(:,:) = x850(:,:)*a2(:,:) 
                 used=send_data(idiag%id_x850, x850, Time)
                 deallocate ( x850 )
            endif
       endif
       if ( idiag%id_w5km>0 ) then
          if (.not.allocated(wz)) allocate ( wz(isc:iec,jsc:jec,npz+1) )
!$omp parallel do default(shared) private(i, j, k)
            do j=jsc,jec
               do i=isc,iec
                  wz(i,j,npz+1) = Atm(n)%phis(i,j)/grav
               enddo
               do k=npz,1,-1
                  do i=isc,iec
                     wz(i,j,k) = wz(i,j,k+1) - Atm(n)%delz(i,j,k)
                  enddo
               enddo
            enddo
            call interpolate_z(isc, iec, jsc, jec, npz, 5.e3, wz, Atm(n)%w(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_w5km, a2, Time)
            if(prt_minmax) call prt_maxmin('W5km', a2, isc, iec, jsc, jec, 0, 1, 1.)
           deallocate (wz)
       endif
       if ( idiag%id_t850>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      850.e2, Atm(n)%peln, Atm(n)%pt(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_t850, a2, Time)
       endif
       if ( idiag%id_q850>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      850.e2, Atm(n)%peln, Atm(n)%q(isc:iec,jsc:jec,:,sphum), a2)
            used=send_data(idiag%id_q850, a2, Time)
       endif
       if ( idiag%id_omg850>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      850.e2, Atm(n)%peln, Atm(n)%omga(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_omg850, a2, Time)
       endif
! 1000-mb
       if ( idiag%id_u1000>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      1000.e2, Atm(n)%peln, Atm(n)%ua(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_u1000, a2, Time)
       endif
       if ( idiag%id_v1000>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      1000.e2, Atm(n)%peln, Atm(n)%va(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_v1000, a2, Time)
       endif
       if ( idiag%id_t1000>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      1000.e2, Atm(n)%peln, Atm(n)%pt(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_t1000, a2, Time)
       endif
       if ( idiag%id_q1000>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      1000.e2, Atm(n)%peln, Atm(n)%q(isc:iec,jsc:jec,:,sphum), a2)
            used=send_data(idiag%id_q1000, a2, Time)
       endif
       if ( idiag%id_omg1000>0 ) then
            call interpolate_vertical(isc, iec, jsc, jec, npz,   &
                                      1000.e2, Atm(n)%peln, Atm(n)%omga(isc:iec,jsc:jec,:), a2)
            used=send_data(idiag%id_omg1000, a2, Time)
       endif

       if ( .not.Atm(n)%flagstruct%hydrostatic .and. idiag%id_w>0  )     &
                 used=send_data(idiag%id_w, Atm(n)%w(isc:iec,jsc:jec,:), Time)

       if(idiag%id_pt   > 0) used=send_data(idiag%id_pt  , Atm(n)%pt  (isc:iec,jsc:jec,:), Time)
       if(idiag%id_omga > 0) used=send_data(idiag%id_omga, Atm(n)%omga(isc:iec,jsc:jec,:), Time)

       if(idiag%id_ppt> 0) then
! Potential temperature perturbation for gravity wave test_case
          allocate ( idiag%pt1(npz) )
          if( .not. allocated(a3) ) allocate ( a3(isc:iec,jsc:jec,npz) )
#ifdef TEST_GWAVES
          call gw_1d(npz, 1000.E2, Atm(n)%ak, Atm(n)%ak, Atm(n)%ak(1), 10.E3, idiag%pt1)
#else
          idiag%pt1 = 0. 
#endif
          do k=1,npz
          do j=jsc,jec
             do i=isc,iec
!               wk(i,j,k) =  (Atm(n)%pt(i,j,k)-300.)/Atm(n)%pkz(i,j,k) * pk0
                wk(i,j,k) =  (Atm(n)%pt(i,j,k)/Atm(n)%pkz(i,j,k) - idiag%pt1(k)) * pk0
             enddo
          enddo
          enddo
          used=send_data(idiag%id_ppt, wk, Time)

          if( prt_minmax ) then
              call prt_maxmin('PoTemp', wk, isc, iec, jsc, jec, 0, npz, 1.)
          endif

          if( allocated(a3) ) deallocate ( a3 )
          deallocate ( idiag%pt1 )
       endif


        do itrac=1, Atm(n)%ncnst
          call get_tracer_names (MODEL_ATMOS, itrac, tname)
          if (idiag%id_tracer(itrac) > 0 .and. itrac .gt. nq) then
            used = send_data (idiag%id_tracer(itrac), Atm(n)%qdiag(isc:iec,jsc:jec,:,itrac), Time )
            if( prt_minmax ) call prt_maxmin(trim(tname), Atm(n)%qdiag(:,:,1,itrac), &
                              isc, iec, jsc, jec, ngc, npz, 1.)
          else
            used = send_data (idiag%id_tracer(itrac), Atm(n)%q(isc:iec,jsc:jec,:,itrac), Time )
            if( prt_minmax ) call prt_maxmin(trim(tname), Atm(n)%q(:,:,1,itrac), &
                              isc, iec, jsc, jec, ngc, npz, 1.)
          endif
!-------------------------------
! ESM TRACER diagnostics output:
! jgj: per SJ email (jul 17 2008): q_dry = q_moist/(1-sphum)
! mass mixing ratio: q_dry = mass_tracer/mass_dryair = mass_tracer/(mass_air - mass_water) ~ q_moist/(1-sphum)
! co2_mmr = (wco2/wair) * co2_vmr
! Note: There is a check to ensure tracer number one is sphum

          if (idiag%id_tracer_dmmr(itrac) > 0 .or. idiag%id_tracer_dvmr(itrac) > 0) then
              if (itrac .gt. nq) then
                dmmr(:,:,:) = Atm(n)%qdiag(isc:iec,jsc:jec,1:npz,itrac)  &
                              /(1.0-Atm(n)%q(isc:iec,jsc:jec,1:npz,1))
              else
                dmmr(:,:,:) = Atm(n)%q(isc:iec,jsc:jec,1:npz,itrac)  &
                              /(1.0-Atm(n)%q(isc:iec,jsc:jec,1:npz,1))
              endif
              dvmr(:,:,:) = dmmr(isc:iec,jsc:jec,1:npz) * WTMAIR/idiag%w_mr(itrac)
              used = send_data (idiag%id_tracer_dmmr(itrac), dmmr, Time )
              used = send_data (idiag%id_tracer_dvmr(itrac), dvmr, Time )
              if( prt_minmax ) then
                 call prt_maxmin(trim(tname)//'_dmmr', dmmr, &
                    isc, iec, jsc, jec, 0, npz, 1.)
                 call prt_maxmin(trim(tname)//'_dvmr', dvmr, & 
                    isc, iec, jsc, jec, 0, npz, 1.)
            endif
          endif
        enddo



   ! enddo  ! end ntileMe do-loop

    deallocate ( a2 )
    deallocate ( u2 )
    deallocate ( v2 )
    deallocate ( wk )

    if ( any(idiag%id_tracer_dmmr > 0) .or. any(idiag%id_tracer_dvmr > 0) ) then
        deallocate (dmmr)
        deallocate ( dvmr)
    endif

    call nullify_domain()


 end subroutine fv_diag

 subroutine wind_max(isc, iec, jsc, jec ,isd, ied, jsd, jed, us, vs, ws_max, domain)
 integer isc, iec, jsc, jec
 integer isd, ied, jsd, jed
 real, intent(in), dimension(isc:iec,jsc:jec):: us, vs
 real, intent(out) :: ws_max(isc:iec,jsc:jec)
 type(domain2d), intent(INOUT) :: domain
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


 subroutine get_vorticity(isc, iec, jsc, jec ,isd, ied, jsd, jed, npz, u, v, vort, dx, dy, rarea)
 integer isd, ied, jsd, jed, npz
 integer isc, iec, jsc, jec
 real, intent(in)  :: u(isd:ied, jsd:jed+1, npz), v(isd:ied+1, jsd:jed, npz)
 real, intent(out) :: vort(isc:iec, jsc:jec, npz)
 real, intent(IN) :: dx(isd:ied,jsd:jed+1)
 real, intent(IN) :: dy(isd:ied+1,jsd:jed)
 real, intent(IN) :: rarea(isd:ied,jsd:jed)
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
            wz(i,j,km+1) = idiag%zsurf(i,j)
         enddo
         do k=km,1,-1
            do i=is,ie
               wz(i,j,k) = wz(i,j,k+1) + gg*pt(i,j,k)*(1.+zvir*q(i,j,k,sphum))  &
                          *(peln(i,k+1,j)-peln(i,k,j))
            enddo
         enddo
      enddo

 end subroutine get_height_field

 subroutine range_check(qname, q, is, ie, js, je, n_g, km, pos, q_low, q_hi, bad_range)
      character(len=*), intent(in)::  qname
      integer, intent(in):: is, ie, js, je
      integer, intent(in):: n_g, km
      real, intent(in)::    q(is-n_g:ie+n_g, js-n_g:je+n_g, km)
      real, intent(in):: pos(is-n_g:ie+n_g, js-n_g:je+n_g,2)
      real, intent(in):: q_low, q_hi
      logical, optional, intent(out):: bad_range
!
      real qmin, qmax
      integer i,j,k

      master = (mpp_pe()==mpp_root_pe())
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
          if(master) write(*,*) 'Range_check Warning:', qname, ' max = ', qmax, ' min = ', qmin
          if ( present(bad_range) ) then
               bad_range = .true. 
          endif
      endif

      if ( present(bad_range) ) then
! Print out where the bad value(s) is (are)
         if ( bad_range .EQV. .false. ) return
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

 subroutine prt_maxmin(qname, q, is, ie, js, je, n_g, km, fac)
      character(len=*), intent(in)::  qname
      integer, intent(in):: is, ie, js, je
      integer, intent(in):: n_g, km
      real, intent(in)::    q(is-n_g:ie+n_g, js-n_g:je+n_g, km)
      real, intent(in)::    fac

      real qmin, qmax
      integer i,j,k
      master = (mpp_pe()==mpp_root_pe())

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

      if(master) then
            write(*,*) qname//trim(gn), ' max = ', qmax*fac, ' min = ', qmin*fac
      endif

 end subroutine prt_maxmin

 subroutine prt_mxm(qname, q, is, ie, js, je, n_g, km, fac, area, domain)
      character(len=*), intent(in)::  qname
      integer, intent(in):: is, ie, js, je
      integer, intent(in):: n_g, km
      real, intent(in)::    q(is-n_g:ie+n_g, js-n_g:je+n_g, km)
      real, intent(in)::    fac
! BUG !!!
!     real, intent(IN)::    area(is-n_g:ie+n_g, js-n_g:je+n_g, km)
      real, intent(IN)::    area(is-3:ie+3, js-3:je+3)
      type(domain2d), intent(INOUT) :: domain
!
      real qmin, qmax, gmean
      integer i,j,k

      master = (mpp_pe()==mpp_root_pe())
      qmin = q(is,js,1)
      qmax = qmin
      gmean = 0.

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

! SJL: BUG!!!
!     gmean = g_sum(domain, q(is,js,km), is, ie, js, je, 3, area, 1) 
      gmean = g_sum(domain, q(is:ie,js:je,km), is, ie, js, je, 3, area, 1) 

      if(master) write(6,*) qname, qmax*fac, qmin*fac, gmean*fac

 end subroutine prt_mxm

 !Added nwat == 1 case for water vapor diagnostics
 subroutine prt_mass(km, nq, is, ie, js, je, n_g, nwat, ps, delp, q, area, domain)

 integer, intent(in):: is, ie, js, je
 integer, intent(in):: nq, n_g, km, nwat
 real, intent(in)::   ps(is-n_g:ie+n_g, js-n_g:je+n_g)
 real, intent(in):: delp(is-n_g:ie+n_g, js-n_g:je+n_g, km)
 real, intent(in)::  q(is-n_g:ie+n_g, js-n_g:je+n_g, km, nq)
 real, intent(IN):: area(is-n_g:ie+n_g,js-n_g:je+n_g)
 type(domain2d), intent(INOUT) :: domain
! Local:
 real psq(is:ie,js:je,nwat), psqv(is:ie,js:je)
 real q_strat(is:ie,js:je)
 real qtot(nwat), qwat
 real psmo, totw, psdry
 integer k, n, kstrat
 master = (mpp_pe()==mpp_root_pe())

#if defined(MARS_GCM) || defined(VENUS_GCM)
 psmo = g_sum(domain, ps(is:ie,js:je), is, ie, js, je, n_g, area, 1)
 totw  = 0.0
 psdry = psmo - totw
 if ( nwat > 0 ) then
    qtot(:)= 0.0
 endif

 if( master ) then
     write(*,*) 'Total surface pressure (mb)', trim(gn), ' = ',  0.01*psmo
!!!     write(*,*) 'mean dry surface pressure = ',    0.01*psdry
  endif
#else

 if ( nwat==0 ) then
      psmo = g_sum(domain, ps(is:ie,js:je), is, ie, js, je, n_g, area, 1) 
      if( master ) write(*,*) 'Total surface pressure (mb)', trim(gn), ' = ',  0.01*psmo
      call z_sum(is, ie, js, je, km, n_g, delp, q(is-n_g,js-n_g,1,1  ), psqv(is,js), domain) 
!     qwat = g_sum(domain, psqv(is,js), is, ie, js, je, n_g, area, 1) 
!     if( master ) write(*,*) 'Total Water Vapor (kg/m**2)', trim(gn), ' =',  qwat*ginv
      return
 endif

 call z_sum(is, ie, js, je, km, n_g, delp, q(is-n_g,js-n_g,1,sphum  ), psq(is,js,sphum  ), domain) 
 if (nwat > 1) then
 call z_sum(is, ie, js, je, km, n_g, delp, q(is-n_g,js-n_g,1,liq_wat), psq(is,js,liq_wat), domain)
 call z_sum(is, ie, js, je, km, n_g, delp, q(is-n_g,js-n_g,1,ice_wat), psq(is,js,ice_wat), domain)
 endif

! Mean water vapor in the "stratosphere" (75 mb and above):
 if ( idiag%phalf(2)< 75. ) then
 kstrat = 1
 do k=1,km
    if ( idiag%phalf(k+1) > 75. ) exit
    kstrat = k
 enddo
 call z_sum(is, ie, js, je, kstrat, n_g, delp, q(is-n_g,js-n_g,1,sphum), q_strat(is,js), domain) 
 psmo = g_sum(domain, q_strat(is,js), is, ie, js, je, n_g, area, 1) * 1.e6           &
      / p_sum(is, ie, js, je, kstrat, n_g, delp, area, domain)
 if(master) write(*,*) 'Mean specific humidity (mg/kg) above 75 mb', trim(gn), '=', psmo
 endif



 if ( nwat==6 ) then
     call z_sum(is, ie, js, je, km, n_g, delp, q(is-n_g,js-n_g,1,rainwat), psq(is,js,rainwat), domain)
     call z_sum(is, ie, js, je, km, n_g, delp, q(is-n_g,js-n_g,1,snowwat), psq(is,js,snowwat), domain)
     call z_sum(is, ie, js, je, km, n_g, delp, q(is-n_g,js-n_g,1,graupel), psq(is,js,graupel), domain)
 endif

!-------------------
! Check global means
!-------------------
 psmo = g_sum(domain, ps(is:ie,js:je), is, ie, js, je, n_g, area, 1) 

 do n=1,nwat
    qtot(n) = g_sum(domain, psq(is,js,n), is, ie, js, je, n_g, area, 1) 
 enddo

 totw  = sum(qtot(1:nwat))
 psdry = psmo - totw

 if( master ) then
     write(*,*) 'Total surface pressure (mb)', trim(gn), ' = ',  0.01*psmo
     write(*,*) 'mean dry surface pressure', trim(gn), ' = ',    0.01*psdry
     write(*,*) 'Total Water Vapor (kg/m**2)', trim(gn), ' =',  qtot(sphum)*ginv
     if ( nwat==6 ) then
          write(*,*) '--- Micro Phys water substances (kg/m**2) ---'
          write(*,*) 'Total cloud water', trim(gn), '=', qtot(liq_wat)*ginv
          write(*,*) 'Total rain  water', trim(gn), '=', qtot(rainwat)*ginv
          write(*,*) 'Total cloud ice  ', trim(gn), '=', qtot(ice_wat)*ginv
          write(*,*) 'Total snow       ', trim(gn), '=', qtot(snowwat)*ginv
          write(*,*) 'Total graupel    ', trim(gn), '=', qtot(graupel)*ginv
          write(*,*) '---------------------------------------------'
     endif
  endif

#endif MARS_GCM

 end subroutine prt_mass


 subroutine z_sum(is, ie, js, je, km, n_g, delp, q, sum2, domain)
 integer, intent(in):: is, ie, js, je,  n_g, km
 real, intent(in):: delp(is-n_g:ie+n_g, js-n_g:je+n_g, km)
 real, intent(in)::    q(is-n_g:ie+n_g, js-n_g:je+n_g, km)
 real, intent(out):: sum2(is:ie,js:je)
 type(domain2d), intent(INOUT) :: domain

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


 real function p_sum(is, ie, js, je, km, n_g, delp, area, domain)
 integer, intent(in):: is, ie, js, je,  n_g, km
 real, intent(in):: delp(is-n_g:ie+n_g, js-n_g:je+n_g, km)
 real, intent(IN) :: area(is-n_g:ie+n_g, js-n_g:je+n_g)
 real :: sum2(is:ie,js:je)
 integer i,j,k
 type(domain2d), intent(INOUT) :: domain

!$omp parallel do default(shared)
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
 p_sum = g_sum(domain, sum2, is, ie, js, je, n_g, area, 1)

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

!$omp parallel do default(shared) private(i, j, k, pm)
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
 subroutine interpolate_z(is, ie, js, je, km, zl, hght, a3, a2)

 integer,  intent(in):: is, ie, js, je, km
 real, intent(in):: hght(is:ie,js:je,km+1)  ! hght(k) > hght(k+1)
 real, intent(in):: a3(is:ie,js:je,km)
 real, intent(in):: zl
 real, intent(out):: a2(is:ie,js:je)
! local:
 real zm(km)
 integer i,j,k


!$omp parallel do default(shared) private(i, j, k, zm)
 do j=js,je
    do 1000 i=is,ie
       do k=1,km
          zm(k) = 0.5*(hght(i,j,k)+hght(i,j,k+1))
       enddo
       if( zl >= zm(1) ) then
           a2(i,j) = a3(i,j,1)
       elseif ( zl <= zm(km) ) then
           a2(i,j) = a3(i,j,km)
       else 
           do k=1,km-1
              if( zl <= zm(k) .and. zl >= zm(k+1) ) then
                  a2(i,j) = a3(i,j,k) + (a3(i,j,k+1)-a3(i,j,k))*(zm(k)-zl)/(zm(k)-zm(k+1))
                  go to 1000
              endif
           enddo
       endif
1000   continue
 enddo

 end subroutine interpolate_z

 subroutine helicity_relative(is, ie, js, je, ng, km, zvir, sphum, srh,   &
                              ua, va, delz, q, hydrostatic, pt, peln, phis, grav)
! !INPUT PARAMETERS:
   integer, intent(in):: is, ie, js, je, ng, km, sphum
   real, intent(in):: grav, zvir
   real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng,km):: pt, ua, va
   real, intent(in):: delz(is-ng:ie+ng,js-ng:je+ng,km)
   real, intent(in):: q(is-ng:ie+ng,js-ng:je+ng,km,*)
   real, intent(in):: phis(is-ng:ie+ng,js-ng:je+ng)
   real, intent(in):: peln(is:ie,km+1,js:je) 
   logical, intent(in):: hydrostatic
   real, intent(out):: srh(is:ie,js:je)   ! unit: (m/s)**2
   real, parameter:: z_crit = 3.e3   ! lowest 3-km
!---------------------------------------------------------------------------------
! SRH = 150-299 ... supercells possible with weak tornadoes
! SRH = 300-449 ... very favourable to supercells development and strong tornadoes
! SRH > 450 ... violent tornadoes
!---------------------------------------------------------------------------------
! if z_crit = 1E3, the threshold for supercells is 100 (m/s)**2
! Coded by S.-J. Lin for CONUS regional climate simulations
!
   real:: rdg
   real, dimension(is:ie):: zh, uc, vc, dz
   integer i, j, k, k0

   rdg = rdgas / grav

!$omp parallel do default(shared), private(i,j,k,zh,uc,vc,dz,k0)
   do j=js,je

      do i=is,ie
         uc(i) = 0.
         vc(i) = 0.
         zh(i) = 0.
         srh(i,j) = 0.

!        if ( phis(i,j)/grav < 1.E3 ) then
         do k=km,1,-1
            if ( hydrostatic ) then
                 dz(i) = rdg*pt(i,j,k)*(1.+zvir*q(i,j,k,sphum))*(peln(i,k+1,j)-peln(i,k,j))
            else
                 dz(i) = - delz(i,j,k)
            endif
            zh(i) = zh(i) + dz(i)
! Compute mean winds below z_crit
            if ( zh(i) < z_crit ) then
                uc(i) = uc(i) + ua(i,j,k)*dz(i)
                vc(i) = vc(i) + va(i,j,k)*dz(i)
                k0 = k
            else
                uc(i) = uc(i) / (zh(i)-dz(i))
                vc(i) = vc(i) / (zh(i)-dz(i))
                goto 123
            endif
         enddo
123      continue

! Lowest layer wind shear computed betw top edge and mid-layer
         k = km
         srh(i,j) = 0.5*(va(i,j,km)-vc(i))*(ua(i,j,km-1)-ua(i,j,km))  -  &
                    0.5*(ua(i,j,km)-uc(i))*(va(i,j,km-1)-va(i,j,km))
         do k=k0, km-1
            srh(i,j) = srh(i,j) + 0.5*(va(i,j,k)-vc(i))*(ua(i,j,k-1)-ua(i,j,k+1)) -  &
                                  0.5*(ua(i,j,k)-uc(i))*(va(i,j,k-1)-va(i,j,k+1))
         enddo
!        endif
      enddo  ! i-loop
   enddo   ! j-loop

 end subroutine helicity_relative



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
!$omp parallel do default(shared), private(i, j) 
        do j=js,je
          do i=is,ie
            vort(i,j,1) = grav * (vort(i,j,1)+f_d(i,j)) / delp(i,j,1)
          enddo
        enddo
#else
! Compute PT at layer edges.
!$omp parallel do default(shared), private(i, j, k, t2, delp2) 
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

!$omp parallel do default(shared), private(i, j, k)
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

!#######################################################################

subroutine rh_calc (pfull, t, qv, rh, do_cmip)
  
   real, intent (in),  dimension(:,:) :: pfull, t, qv
   real, intent (out), dimension(:,:) :: rh
   real, dimension(size(t,1),size(t,2)) :: esat
   logical, intent(in), optional :: do_cmip

   real, parameter :: d622 = rdgas/rvgas
   real, parameter :: d378 = 1.-d622

   logical :: do_simple = .false.

! because Betts-Miller uses a simplified scheme for calculating the relative humidity

     if (do_simple) then
        call lookup_es (t, esat)
        rh(:,:) = pfull(:,:)
        rh(:,:) = MAX(rh(:,:),esat(:,:))  !limit where pfull ~ esat
        rh(:,:)=100.*qv(:,:)/(d622*esat(:,:)/rh(:,:))
     else
        if (present(do_cmip)) then
         call compute_qs (t, pfull, rh, q=qv, es_over_liq_and_ice = .true.)
         rh(:,:)=100.*qv(:,:)/rh(:,:)
        else
        call compute_qs (t, pfull, rh, q=qv)
        rh(:,:)=100.*qv(:,:)/rh(:,:)
     endif
     endif

end subroutine rh_calc

!#######################################################################

end module fv_diagnostics_mod
