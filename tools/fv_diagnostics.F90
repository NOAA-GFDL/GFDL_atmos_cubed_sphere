module fv_diagnostics_mod

 use constants_mod,    only: grav, rdgas, pi, radius
 use fms_io_mod,       only: set_domain, nullify_domain
 use time_manager_mod, only: time_type, get_date, get_time
 use mpp_domains_mod,  only: domain2d, mpp_update_domains, DGRID_NE
 use diag_manager_mod, only: diag_axis_init, register_diag_field, &
                             register_static_field, send_data
#ifdef USE_STDOUT
 use mpp_mod,          only: stdout
#endif
 use fv_arrays_mod,    only: fv_atmos_type
 use mapz_module,      only: E_Flux
 use mp_mod,           only: domain, gid, masterproc, &
                             mp_reduce_sum, mp_reduce_min, mp_reduce_max
 use eta_mod,          only: get_eta_level
 use grid_tools,       only: dx, dy, dxa, dya, area, rarea
 use grid_utils,       only: f0, cosa_s, cubed_to_latlon, g_sum, sina_u, sina_v,   &
                             en1, en2, vlon
 use a2b_edge_mod,     only: a2b_ord4
 use sw_core,          only: d2a2c_vect
 use surf_map,         only: zs_g

 use tracer_manager_mod, only: get_tracer_names, get_number_tracers
 use field_manager_mod,  only: MODEL_ATMOS
 use fms_mod,            only: error_mesg, FATAL, stdlog

 implicit none
 private

 logical master
 integer ::id_ps, id_slp, id_h500, id_ua, id_va, id_pt, id_omga, id_vort,  &
           id_pv, id_zsurf, id_oro, id_sgh, id_prec, id_divg, id_u, id_v, id_w,     &
           id_te, id_zs, id_mq
 integer, parameter:: max_step = 1000
 integer steps
 real*4:: efx(max_step), mtq(max_step)
 real*4:: efx_sum,       mtq_sum
! For initial conditions:
 integer ic_ps, ic_ua, ic_va
 integer, allocatable :: id_tracer(:), id_tracer_tend(:)

 integer  :: ncnst
 real :: missing_value = -1.e10
 real :: ginv
 real, allocatable :: phalf(:)
 real, allocatable :: zsurf(:,:)
 real, allocatable :: zxg(:,:)

 type(time_type) :: fv_time

 logical :: module_is_initialized=.false.
 logical :: full_phys
 real    :: ptop

 public :: fv_diag_init, fv_time, fv_diag, prt_maxmin, id_prec, id_divg, id_te
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
    real              :: vrange(2), wrange(2), trange(2), slprange(2)
    real              :: pfull(npz)

    integer :: id_bk, id_pk, id_lon, id_lat, id_lont, id_latt, id_phalf, id_pfull
    integer :: i, j, n, ntileMe, id_xt, id_yt, id_x, id_y, id_xe, id_ye, id_xn, id_yn
    integer :: isc, iec, jsc, jec

    logical :: used

    character(len=64) :: field

! tracers
    character(len=128)   :: tname
    character(len=256)   :: tlongname, tunits
    integer              :: ntprog

! For total energy diagnostics:
    steps = 0
    efx = 0.;       efx_sum = 0.
    mtq = 0.;       mtq_sum = 0.

    ncnst = Atm(1)%ncnst

    call set_domain(Atm(1)%domain)  ! Set domain so that diag_manager can access tile information

! valid range for some fields

!!!  This will need mods for more than 1 tile per pe  !!!

    vrange = (/ -330.,  330. /)  ! winds
    wrange = (/  -50.,   50. /)  ! vertical wind

#ifdef MARS_GCM
    slprange = (/0.,  100./)  ! sea-level-pressure
    trange = (/  50., 360. /)  ! temperature
#else
    trange = (/  100.,  350. /)  ! temperature
    slprange = (/800.,  1200./)  ! sea-level-pressure
#endif


    ginv = 1./GRAV
    fv_time = Time  ! Fix the fv_time is used

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

!--- Send static data

    if ( id_bk > 0 )    used = send_data ( id_bk,Atm(1)%bk, Time )
    if ( id_pk > 0 )    used = send_data ( id_pk,Atm(1)%ak, Time )

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

#ifndef DYNAMICS_ZS
       id_zsurf = register_static_field ( trim(field), 'zsurf', axes(1:2),  &
                                         'surface height', 'm' )
#endif
       id_zs = register_static_field ( trim(field), 'zs', axes(1:2),  &
                                        'Original Mean Terrain', 'm' )
! For mountain torque in zonal dir:
       id_oro = register_static_field ( trim(field), 'oro', axes(1:2),  &
                                        'Land/Water Mask', 'none' )
       id_sgh = register_static_field ( trim(field), 'sgh', axes(1:2),  &
                                        'Terrain Standard deviation', 'm' )
       ic_ps  = register_static_field ( trim(field), 'ps_ic', axes(1:2),  &
                                         'initial surface pressure', 'Pa' )

       ic_ua = register_static_field ( trim(field), 'ua_ic', axes(1:3),        &
            'zonal wind', 'm/sec' )
       ic_va = register_static_field ( trim(field), 'va_ic', axes(1:3),        &
            'meridional wind', 'm/sec' )

    end do

    master = (gid == masterproc)

    n=1
    isc = Atm(n)%isc; iec = Atm(n)%iec
    jsc = Atm(n)%jsc; jec = Atm(n)%jec

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
#ifndef DYNAMICS_ZS
       if (id_zsurf > 0) used = send_data(id_zsurf, zsurf, Time)
#endif
#ifdef FV_LAND
       if (id_zs  > 0) used = send_data(id_zs , zs_g, Time)
       if (id_oro > 0) used = send_data(id_oro, Atm(n)%oro(isc:iec,jsc:jec), Time)
       if (id_sgh > 0) used = send_data(id_sgh, Atm(n)%sgh(isc:iec,jsc:jec), Time)
#endif
#ifdef SW_DYNAMICS
       Atm(n)%ps(isc:iec,jsc:jec) = ginv * Atm(n)%delp(isc:iec,jsc:jec,1)
#endif
       if (ic_ps  > 0) used = send_data(ic_ps, Atm(n)%ps(isc:iec,jsc:jec), Time)

#ifdef SW_DYNAMICS
       if (ic_ua>0 .or. ic_va>0) then
          call cubed_to_latlon(Atm(n)%u, Atm(n)%v, Atm(n)%ua, Atm(n)%va,   &
                               dx, dy, dxa, dya, 1)
       endif
#endif
       if(ic_ua > 0) used=send_data(ic_ua, Atm(n)%ua(isc:iec,jsc:jec,:), Time)
       if(ic_va > 0) used=send_data(ic_va, Atm(n)%va(isc:iec,jsc:jec,:), Time)
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
! 500 mb Height
!--------------
      id_h500 = register_diag_field (trim(field), 'h500', axes(1:2),  Time,   &
                                     '500-mb hght', 'm', missing_value=missing_value )

!-------------------
! Sea-level-pressure
!-------------------
       id_slp = register_diag_field (trim(field), 'slp', axes(1:2),  Time,   &
                                     'sea-level pressure', 'mb', missing_value=missing_value,  &
                                      range=slprange )
!----------------------------
! Precip from simple physics
!----------------------------
      id_prec = register_diag_field (trim(field), 'prec', axes(1:2),  Time,   &
                                     'Precip_Simple', 'mm/day', missing_value=missing_value )
!-------------------
! A grid winds (lat-lon)
!-------------------
       id_u  = register_diag_field ( trim(field), 'ucoor', axes(1:3), Time,        &
            'u wind', 'm/sec', missing_value=missing_value, range=vrange )
       id_v  = register_diag_field ( trim(field), 'vcoor', axes(1:3), Time,        &
            'v wind', 'm/sec', missing_value=missing_value, range=vrange)

       id_ua = register_diag_field ( trim(field), 'ucomp', axes(1:3), Time,        &
            'zonal wind', 'm/sec', missing_value=missing_value, range=vrange )
       id_va = register_diag_field ( trim(field), 'vcomp', axes(1:3), Time,        &
            'meridional wind', 'm/sec', missing_value=missing_value, range=vrange)

       id_w = register_diag_field ( trim(field), 'w', axes(1:3), Time,        &
            'vertical wind', 'm/sec', missing_value=missing_value, range=wrange )

       id_pt   = register_diag_field ( trim(field), 'temp', axes(1:3), Time,       &
            'temperature', 'K', missing_value=missing_value, range=trange )
       id_omga = register_diag_field ( trim(field), 'omega', axes(1:3), Time,      &
            'omega', 'Pa/s', missing_value=missing_value )
       id_divg  = register_diag_field ( trim(field), 'divg', axes(1:3), Time,      &
            'mean divergence', '1/s', missing_value=missing_value )
! Total energy (moist only when full_phys = .T.)
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

        do i = 1, ncnst
           call get_tracer_names ( MODEL_ATMOS, i, tname, tlongname, tunits )
           id_tracer(i) = register_diag_field ( field, trim(tname),  &
                axes(1:3), Time, trim(tlongname), &
                trim(tunits), missing_value=missing_value)
           if (master) then
               if (id_tracer(i) > 0) then
                   write(stdlog(),'(a,a,a,a)') &
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
    integer :: ngc

    real, allocatable :: a2(:,:), wk(:,:,:), wz(:,:,:), ucoor(:,:,:), vcoor(:,:,:)
    real height(2)
    real tot_mq 
    logical :: used
    logical :: prt_minmax
    integer i,j,  yr, mon, dd, hr, mn, days, seconds
    character(len=128)   :: tname


    height(1) = 5.E3      ! for computing 5-km "pressure"
    height(2) = 0.        ! for sea-level pressure

    ntileMe = size(Atm(:))
    n = 1
    isc = Atm(n)%isc; iec = Atm(n)%iec
    jsc = Atm(n)%jsc; jec = Atm(n)%jec
    ngc = Atm(n)%ng
    npz = Atm(n)%npz
    full_phys = Atm(n)%full_phys
    ptop = Atm(n)%ak(1)

    allocate( a2(isc:iec,jsc:jec) )


    fv_time = Time
    call set_domain(Atm(1)%domain)

    if ( full_phys ) then
#ifdef MARS_GCM
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
         if ( full_phys ) then
#ifdef USE_STDOUT
              write(stdout(),*) yr, mon, dd, hr, mn, seconds
#else
              if(master) write(6,*) yr, mon, dd, hr, mn, seconds
#endif
         else
#ifdef USE_STDOUT
              write(stdout(),*) Days, seconds
#else
              if(master) write(6,*) Days, seconds
#endif
         endif
     endif

     if(id_u>0 .or. id_v>0) then
        isd = Atm(n)%isd; ied = Atm(n)%ied
        jsd = Atm(n)%jsd; jed = Atm(n)%jed
        allocate(ucoor(isc:iec,jsc:jec,npz), vcoor(isc:iec,jsc:jec,npz) )
        call d2a_vect(Atm(n)%u, Atm(n)%v, ucoor, vcoor)
        if(id_u > 0) used=send_data(id_u, ucoor, Time)
        if(id_v > 0) used=send_data(id_v, vcoor, Time)
        deallocate(ucoor, vcoor)
     end if
 
#ifdef SW_DYNAMICS
       if (id_ua>0 .or. id_va>0) then
          isd = Atm(n)%isd; ied = Atm(n)%ied
          jsd = Atm(n)%jsd; jed = Atm(n)%jed
          call cubed_to_latlon(Atm(n)%u, Atm(n)%v, Atm(n)%ua, Atm(n)%va,   &
                               dx, dy, dxa, dya, 1)
       endif
#endif

    if( prt_minmax ) then
        call prt_mass(npz, ncnst, isc, iec, jsc, jec, ngc, Atm(n)%ps, Atm(n)%delp, Atm(n)%q, master)
#ifdef SW_DYNAMICS
        call prt_maxmin('PS', Atm(n)%ps, isc, iec, jsc, jec, ngc, 1, 1./grav, master)
#else
             steps = steps + 1
           efx_sum = efx_sum + E_Flux
        if ( steps <= max_step ) efx(steps) = E_Flux
        if (master)  then
            write(6,*) 'ENG Deficit (W/m**2)=', E_Flux
        endif

        call prt_maxmin('PS', Atm(n)%ps, isc, iec, jsc, jec, ngc, 1, 0.01, master)
#endif
!       call prt_maxmin('ZS', zsurf, isc, iec, jsc, jec, 0, 1, 1., master)
        call prt_maxmin('UA', Atm(n)%ua, isc, iec, jsc, jec, ngc, npz, 1., master)
        call prt_maxmin('VA', Atm(n)%va, isc, iec, jsc, jec, ngc, npz, 1., master)

        if ( .not. Atm(n)%hydrostatic )   &
        call prt_maxmin('W ', Atm(n)%w , isc, iec, jsc, jec, ngc, npz, 1., master)

#ifndef SW_DYNAMICS
        call prt_maxmin('TA', Atm(n)%pt,   isc, iec, jsc, jec, ngc, npz, 1., master)
        call prt_maxmin('OM', Atm(n)%omga, isc, iec, jsc, jec, ngc, npz, 1., master)

! Tracers:
! Have moved this to further down 
!       if ( ncnst >= 1 )    &
!       call prt_maxmin('Q1', Atm(n)%q(isc-ngc,jsc-ngc,1,1), isc, iec, jsc, jec, ngc, npz, 1., master)
!       if ( ncnst >= 2 )    &
!       call prt_maxmin('Q2', Atm(n)%q(isc-ngc,jsc-ngc,1,2), isc, iec, jsc, jec, ngc, npz, 1., master)
!       if ( ncnst >= 3 )    &
!       call prt_maxmin('Q3', Atm(n)%q(isc-ngc,jsc-ngc,1,3), isc, iec, jsc, jec, ngc, npz, 1., master)
!       if ( ncnst >= 4 )    &
!       call prt_maxmin('Q4', Atm(n)%q(isc-ngc,jsc-ngc,1,4), isc, iec, jsc, jec, ngc, npz, 1., master)
#endif
    endif


    do n = 1, ntileMe

       isc = Atm(n)%isc; iec = Atm(n)%iec
       jsc = Atm(n)%jsc; jec = Atm(n)%jec

#ifdef DYNAMICS_ZS
       if(id_zsurf > 0)  used=send_data(id_zsurf, zsurf, Time)
#endif
       if(id_ps > 0) used=send_data(id_ps, Atm(n)%ps(isc:iec,jsc:jec), Time)


       if(id_slp > 0 .or. id_h500 > 0 ) then

          allocate ( wz(isc:iec,jsc:jec,npz+1) )
          call get_height_field(isc, iec, jsc, jec, ngc, npz, wz, Atm(n)%pt, Atm(n)%q, Atm(n)%peln, zvir)
          if( prt_minmax )   &
          call prt_maxmin('ZTOP', wz(isc,jsc,1), isc, iec, jsc, jec, 0, 1, 1.E-3, master)

          if(id_slp > 0) then
! Cumpute SLP (pressure at height=0)
          call get_pressure_given_height(isc, iec, jsc, jec, ngc, npz, wz, 1, height(2),   &
                                        Atm(n)%pt(:,:,npz), Atm(n)%peln, a2, 0.01)
          used = send_data (id_slp, a2, Time)
          if( prt_minmax )   &
             call prt_maxmin('SLP', a2, isc, iec, jsc, jec, 0, 1, 1., master)
          endif

! Compute H500
          if(id_h500 > 0) then
             height(1) = log( 50000. )
             call get_height_given_pressure(isc, iec, jsc, jec, ngc, npz, wz, 1, height(1),   &
                                            Atm(n)%peln, a2)
             used = send_data ( id_h500, a2, Time )
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
              tot_mq  = g_sum( a2, isc, iec, jsc, jec, ngc, area) 
              mtq_sum = mtq_sum + tot_mq
              if ( steps <= max_step ) mtq(steps) = tot_mq
              if(master) write(*,*) 'Total (global) mountain torque (Hadleys)=', tot_mq
          endif
       endif


       if(id_ua > 0) used=send_data(id_ua, Atm(n)%ua(isc:iec,jsc:jec,:), Time)
       if(id_va > 0) used=send_data(id_va, Atm(n)%va(isc:iec,jsc:jec,:), Time)

       if ( .not. Atm(n)%hydrostatic .and. id_w>0  )     &
                     used=send_data(id_w, Atm(n)%w(isc:iec,jsc:jec,:), Time)

       if(id_pt   > 0) used=send_data(id_pt  , Atm(n)%pt  (isc:iec,jsc:jec,:), Time)
       if(id_omga > 0) used=send_data(id_omga, Atm(n)%omga(isc:iec,jsc:jec,:), Time)

       if ( id_vort>0 .or. id_pv>0 ) then
          allocate ( wk(isc:iec,jsc:jec,npz) )
          isd = Atm(n)%isd; ied = Atm(n)%ied
          jsd = Atm(n)%jsd; jed = Atm(n)%jed
          npz = Atm(n)%npz

          if(id_vort >0 .or. id_pv>0) then
             call get_vorticity(Atm(n)%u, Atm(n)%v, wk)
             used=send_data(id_vort, wk, Time)
          endif

         if ( id_pv > 0 ) then
! Note: this is expensive computation.
              call pv_entropy(isc, iec, jsc, jec, ngc, npz, wk,    &
                              f0, Atm(n)%pt, Atm(n)%pkz, Atm(n)%delp, grav)
              used = send_data ( id_pv, wk, Time )
         endif

         deallocate ( wk )
       endif

        do itrac=1, ncnst
          if (id_tracer(itrac) > 0) &
               & used = send_data (id_tracer(itrac), Atm(n)%q(isc:iec,jsc:jec,:,itrac), Time )
          if( prt_minmax ) then
              call get_tracer_names ( MODEL_ATMOS, itrac, tname )
              call prt_maxmin(trim(tname), Atm(n)%q(isc-ngc,jsc-ngc,1,itrac), &
                   isc, iec, jsc, jec, ngc, npz, 1., master)
          endif
        enddo
       deallocate ( a2 )
    enddo

    call nullify_domain()

  contains

    subroutine get_vorticity(u, v, vort)

      real, intent(in)  :: u(isd:ied, jsd:jed+1, npz), v(isd:ied+1, jsd:jed, npz)
      real, intent(out) :: vort(isc:iec, jsc:jec, npz)

      real    :: utmp(isc:iec, jsc:jec+1), vtmp(isc:iec+1, jsc:jec)
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

    subroutine d2a_vect(u, v, ua, va)

      real, intent(in)  :: u(isd:ied, jsd:jed+1, npz), v(isd:ied+1, jsd:jed, npz)
      real, intent(out) :: ua(isc:iec, jsc:jec, npz), va(isc:iec, jsc:jec, npz)

      real :: ut(isc:iec+1, jsc:jec), vt(isc:iec, jsc:jec+1)
      real :: utmp, vtmp, rsin2
      integer :: i,j,k

      do k=1,npz
         do j=jsc,jec+1
            do i=isc,iec
               vt(i,j) = u(i,j,k)*dx(i,j)
            enddo
         enddo
         do j=jsc,jec
            do i=isc,iec+1
               ut(i,j) = v(i,j,k)*dy(i,j)
            enddo
         enddo
         do j=jsc,jec
            do i=isc,iec
               utmp = 0.5*(vt(i,j) + vt(i,j+1)) / dxa(i,j)
               vtmp = 0.5*(ut(i,j) + ut(i+1,j)) / dya(i,j)
               rsin2 = 1. / (1.-cosa_s(i,j)**2)
               ua(i,j,k) = (utmp-vtmp*cosa_s(i,j)) * rsin2
               va(i,j,k) = (vtmp-utmp*cosa_s(i,j)) * rsin2
            enddo
         enddo
      enddo

    end subroutine d2a_vect

 end subroutine fv_diag


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
               wz(i,j,k) = wz(i,j,k+1) + gg*pt(i,j,k)*(1.+zvir*q(i,j,k,1))   &
                          *(peln(i,k+1,j)-peln(i,k,j))
            enddo
         enddo
      enddo

 end subroutine get_height_field

 subroutine prt_maxmin(qname, q, is, ie, js, je, n_g, km, fac, master)
      character*(*), intent(in)::  qname
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

#ifdef USE_STDOUT
      write(stdout(),*) qname, ' max = ', qmax*fac, ' min = ', qmin*fac
#else
      if(master) write(6,*) qname, ' max = ', qmax*fac, ' min = ', qmin*fac
#endif

 end subroutine prt_maxmin

 subroutine prt_mass(km, nq, is, ie, js, je, n_g, ps, delp, q, master)

      integer, intent(in):: is, ie, js, je
      integer, intent(in):: nq, n_g, km
      real, intent(in)::   ps(is-n_g:ie+n_g, js-n_g:je+n_g)
      real, intent(in):: delp(is-n_g:ie+n_g, js-n_g:je+n_g, km)
      real, intent(in)::  q(is-n_g:ie+n_g, js-n_g:je+n_g, km, nq)
      logical, intent(in):: master
! Local:
      real psq(is:ie,js:je)
      real psmo, psdry, qtot
      integer i,j,k, ip

      ip = min(3,size(q,4))

      do 1000 j=js,je
         do i=is,ie
            psq(i,j) = 0.
         enddo

        if( full_phys ) then
          do k=1,km
             do i=is,ie
                psq(i,j) = psq(i,j) + sum( q(i,j,k,1:ip) ) * delp(i,j,k)
            enddo
          enddo
        else
          do k=1,km
             do i=is,ie
                psq(i,j) = psq(i,j) + delp(i,j,k)*q(i,j,k,1)
             enddo
          enddo
        endif
1000  continue

! Check global means
       psmo  = g_sum( ps(is:ie,js:je), is, ie, js, je, n_g, area, mode=1) 
       qtot  = g_sum( psq, is, ie, js, je, n_g, area, mode=1) 
       psdry = psmo - qtot

#ifdef USE_STDOUT
       write(stdout(),*) 'Total surface pressure (mb) = ', 0.01*psmo
       if ( full_phys ) then
            write(stdout(),*) 'mean dry surface pressure = ', 0.01*psdry
            write(stdout(),*) 'Total Column Water (kg/m**2) =', qtot/grav
       else
            write(stdout(),*) 'Total Tracer-1 Mass=', qtot/grav
       endif
#else
       if(master) then
          write(6,*) 'Total surface pressure (mb) = ', 0.01*psmo
       if ( full_phys ) then
            write(6,*) 'mean dry surface pressure = ', 0.01*psdry
            write(6,*) 'Total Column Water (kg/m**2) =', qtot/grav
       else
            write(6,*) 'Total Tracer-1 Mass=', qtot/grav
       endif
       endif
#endif

 end subroutine prt_mass


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

!$omp parallel do private(i,j,k, ptmp, tm)
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

!$omp parallel do private(i,j,k)
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
                 a2(i,j,n) = missing_value
1000   continue
    enddo
 enddo

 end subroutine get_height_given_pressure


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
