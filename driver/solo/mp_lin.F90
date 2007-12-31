module mp_lin_mod
 use fv_diagnostics_mod,only: prt_maxmin
 use fv_grid_tools_mod, only: area
 use fv_grid_utils_mod, only: g_sum, ptop
 use fv_timing_mod,     only: timing_on, timing_off
 use fv_mp_mod,         only: gid, ng
 use mpp_mod,           only: stdlog
 use diag_manager_mod,  only: register_diag_field, send_data
 use time_manager_mod,  only: time_type, get_date, get_time
 use constants_mod,     only: grav, rdgas, rvgas, cp_air, hlv, hlf
 use fms_mod,           only: write_version_number, open_namelist_file, &
                              check_nml_error, file_exist, close_file,  &
                               error_mesg, FATAL 

 implicit none
 private

 character(len=128) :: version = ''
 character(len=128) :: tagname = ''

 public  micro_phys_driver, micro_phys_init, micro_phys_end
 real             :: missing_value = -1.e10
 logical          :: module_is_initialized = .false.
 character(len=7) :: mod_name = 'mp_lin'

!==== fms constants ====================
 real, parameter :: rair  = rdgas        ! gas constant for dry air           (j/kg/k)
 real, parameter :: cp    = cp_air       ! heat capacity at constant pressure (j/kg/k)
 real, parameter :: eps   = rdgas/rvgas
 real, parameter :: latv  = hlv          ! = 2.500e6
 real, parameter :: lati  = hlf          ! = 3.34e5
 real, parameter :: rvapr = rvgas
!==== fms constants ====================
 real, parameter :: qvmin  = 1.e-20
 real, parameter :: qrmin  = 1.e-8
 real, parameter :: sfcrho = 1.20         ! surface air density
 real, parameter :: vmin = 1.e-3          ! minimum fall speed for rain/graupel
 real, parameter ::  tice = 273.16
 real, parameter ::   so3 = 7./3.

 real :: cracs,csacr,cgacr,cgacs,acco(3,4),csacw,               &
         craci,csaci,cgacw,cgaci,cracw,cssub(5),cgsub(5), &
         crevp(5),cgfr(2),csmlt(5),cgmlt(5)
 real :: rmi50, rmi40
 real :: es0, ces0, c1brg, c2brg

! Derived variables:
 real :: dts, rdts, pie
 real :: lcp, icp, tcp, rgrav
 real :: ql0_crt, raut_rate, c_praut
 real :: mp_count = 0.

 logical :: do_setup=.true.
 logical :: master 

 real, allocatable:: vt_r(:,:,:), vt_s(:,:,:), vt_g(:,:,:), vt_i(:,:,:)
 real, allocatable:: prec0(:,:), rain0(:,:), snow0(:,:), ice0(:,:), graupel0(:,:)
 real, allocatable:: prec1(:,:), prec_mp(:,:), cond(:,:) 
 real, allocatable:: table(:), table2(:), tablew(:),  des(:), des2(:), desw(:)

 integer:: isc, iec, jsc, jec
 integer:: id_vtr, id_vts,  id_vtg, id_vti, id_rain, id_snow, id_graupel, &
           id_ice, id_prec, id_cond

!----------------------
! namelist  parameters:
!----------------------
 real :: mp_time = 180.   ! maximum micro-physics time step (sec)
 real ::  tau_i  =  10.   ! (sec) cloud ice melt real quick
 real ::  tau_g  = 180.   ! graupel melt
 real ::  tau_s  = 120.   ! snow melt
 real ::  tau_r  = 180.   ! freezing rain time scale
 real ::  tau_f  =  60.   ! freezing scale for cloud_wat-->cld_ic
 real ::  tau_evap  = 300.   ! 
 real ::  rh_l   = 0.95   ! RH threshold over land
 real ::  rh_o   = 1.00   ! RH threshold over ocean
 real :: evap_eff = 1.0   ! rain re-evaporation efficiency
 real :: evap_rhc = 0.90  ! rain re-evaporation threshold

!-------------------------------------------------------------
! Cld_wat --> rain autoconversion "inverse" rate
! ql0_crt controls the total amount of cloud water 
! *** lower resolution needs smaller value of ql0_crt ***

! * keeping ql0_crt = const; higher c_praut reduces total cloud water
! * Keeping c_praut const; higher ql0_crt produces more cloud water
! * Keeping ql0_crt*c_praut=const; higher ql0_crt produces more cloud water
! *** Higher cloud water usually associated with more TPW
! For land, uses a higher qol_crt while keeping the same c_praut
!real :: ql0_crt = 2.5e-4   ! cld_wat--> rain; was 7.e-4 in Purdue_Lin/WRF
!-------------------------------------------------------------
#ifdef LFO_AUTO
 real :: ql0_ocean = 1.4e-4  !               over ocean
 real :: ql0_land  = 7.0e-4  !  critical ql0 over land
 real :: c_praut_o = 5.0e-3 ! autoconversion: cld_wat --> rain was 1.0e-3
 real :: c_praut_l = 1.0e-3 ! autoconversion: cld_wat --> rain was 1.0e-3
#else
 real :: ql0_land  = 4.0e-4  !  critical ql0 over land
 real :: ql0_ocean = 4.0e-4
 real :: c_praut_l = 6.78
 real :: c_praut_o = 6.78
#endif
!-------------------------------------------------------------
 real :: qs0_crt = 1.0e-3   ! snow -> graupel threshold
 real :: qi0_crt = 1.0e-3   ! ice  -> snow threshold; was 0.6e-3
 real :: qr0_crt = 1.0e-4   ! rain --> snow or graupel threshold

 real :: c_psaut = 1.0e-3   ! autoconversion: cloud_ice --> snow
 real :: c_psaci = 1.0e-2   ! accretion: cloud ice --> snow; was 0.1 (Zetac)
 real :: c_piacr = 2.e2     ! accretion: rain --> ice
 real :: c_pgacs = 0.05     ! snow--> graupel eff. (original setting wss 0.1)

! fall velocity tunning constants:
 real :: vr_fac = 1.
 real :: vs_fac = 1.
 real :: vg_fac = 1.
 real :: vi_fac = 0.

 logical :: use_cld_prof = .false.

! Khairoutdinov and Kogan 2000: autoconversion of cloud water to rain
 logical :: auto_kk2000  = .true.     ! c_praut~6.774 in WSM6/KK2000
                                       ! Larger value produces lower LWP
 logical :: falling_rain_mp = .true.
 logical :: s2g_accr = .true.
 logical :: s2g_auto = .true.
 logical :: mp_debug = .false.
 logical :: mp_print = .true.

 namelist /mp_lin_nml/mp_time, tau_r, tau_evap, tau_s, tau_i, tau_g,  &
                      vr_fac, vs_fac, vg_fac, vi_fac, use_cld_prof,  &
                      auto_kk2000, qs0_crt, qi0_crt, qr0_crt,   &
                      ql0_land, ql0_ocean, rh_l, rh_o,  &
                      c_praut_l, c_praut_o, evap_rhc, evap_eff,   &
                      tau_f, s2g_accr, s2g_auto, c_piacr,           &
                      c_psaut, c_psaci, falling_rain_mp, mp_debug, mp_print

 contains
 

  subroutine micro_phys_driver(qv,    ql,    qr,    qi,    qs,    qg,    qa,    &
                               qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, qa_dt, & 
                               pt_dt, pt, p3, dz,  delp, dt_in,                 &
                               land,  rain, snow, ice, graupel,                 &
                               is,ie, js,je, ks,ke, ktop, kbot, time)

  type(time_type), intent(in):: time
  integer,         intent(in):: is,ie, js,je, ks,ke
  integer,         intent(in):: ktop, kbot  ! vertical compute domain
  real,            intent(in):: dt_in

  real, intent(in), dimension(is:ie,js:je,ks:ke) :: p3, delp
  real, intent(in), dimension(is:ie,js:je):: land  !land fraction

  real, intent(in   ), dimension(is:ie,js:je,ks:ke):: pt, qv, ql, qr, qi, qs, qg, qa
  real, intent(inout), dimension(is:ie,js:je,ks:ke):: pt_dt,  qa_dt, dz
  real, intent(inout), dimension(is:ie,js:je,ks:ke):: qv_dt, ql_dt, qr_dt, qi_dt,  &
                                                      qs_dt, qg_dt
  real, intent(out), dimension(is:ie,js:je):: rain, snow, ice, graupel

! local:
  real :: mpdt, rdt, convt, tot_prec
  integer :: i,j,k
  integer :: seconds, days, ntimes
  logical used

                                        call timing_on (" split_mp")

! tendency zero out for am moist processes should be done outside the driver

     mpdt = min(dt_in, mp_time)
      rdt = 1. / dt_in
   ntimes = nint( dt_in/mpdt )
! small time step:
      dts = dt_in / real(ntimes)
     rdts = 1./dts

  call get_time (time, seconds, days)

  if ( do_setup ) then
      call setup_con (is, ie, js, je, ks, ke)
      call setupm
      do_setup = .false.
  endif

  if ( mp_debug ) then
       call prt_maxmin('T_b_mp',    pt, is, ie, js, je, 0, kbot, 1., master)
       call prt_maxmin('qg_dt_b_mp',  qg_dt, is, ie, js, je, 0, kbot, 1., master)
  endif


  do j=js, je
     do i=is, ie
        graupel(i,j) = 0.
           rain(i,j) = 0.
           snow(i,j) = 0.
            ice(i,j) = 0.
           cond(i,j) = 0.
     enddo
  enddo
  do j=js,je
     call mpdrv( delp, pt, p3, qv, ql, qr, qi, qs, qg, qa, dz,  &
                  is, ie, js, je, ks, ke, ktop, kbot, j, dt_in,  & 
                  ntimes, rain(is,j), snow(is,j), graupel(is,j), &
                  ice(is,j), cond(is,j), land(is,j),  &
                  pt_dt, qv_dt, ql_dt, qr_dt, qi_dt,    &
                  qs_dt, qg_dt, qa_dt )
  enddo

! no clouds allowed above ktop
   if ( ks < ktop ) then
      do k=ks, ktop
         do j=js,je
            do i=is,ie
!              qa(i,j,k) = 0.
               qa_dt(i,j,k) = -qa(i,j,k) * rdt
            enddo
         enddo
      enddo
   endif

   if ( id_vtr> 0 ) used=send_data(id_vtr, vt_r, time)
   if ( id_vts> 0 ) used=send_data(id_vts, vt_s, time)
   if ( id_vtg> 0 ) used=send_data(id_vtg, vt_g, time)
   if ( id_vts> 0 ) used=send_data(id_vti, vt_i, time)

   convt = 86400.*rdt*rgrav
   do j=js,je
      do i=is,ie
            rain(i,j) =    rain(i,j) * convt
            snow(i,j) =    snow(i,j) * convt
             ice(i,j) =     ice(i,j) * convt
         graupel(i,j) = graupel(i,j) * convt
         prec_mp(i,j) =    rain(i,j) + snow(i,j) + ice(i,j) + graupel(i,j)
      enddo
   enddo

   if ( id_cond>0 ) then
        do j=js,je
           do i=is,ie
              cond(i,j) = cond(i,j)*rgrav
           enddo
        enddo
        used=send_data(id_cond, cond, time)
   endif

   if ( id_snow>0 ) then
        used=send_data(id_snow,    snow,    time)
        if ( seconds==0 ) then
!            tot_prec = g_sum(snow, is, ie, js, je, ng, area, 1) 
!            if(master) write(*,*) 'mean snow=', tot_prec
        endif
        snow0(:,:) = snow0(:,:) + snow(:,:)
   endif

   if ( id_graupel>0 ) then
        used=send_data(id_graupel, graupel, time)
        if ( seconds==0 ) then
!            tot_prec = g_sum(graupel, is, ie, js, je, ng, area, 1) 
!            if(master) write(*,*) 'mean graupel=', tot_prec
        endif
        graupel0(:,:) = graupel0(:,:) + graupel(:,:)
   endif

   if ( id_ice>0 ) then
        used=send_data(id_ice, ice, time)
        if ( seconds==0 ) then
!            tot_prec = g_sum(ice, is, ie, js, je, ng, area, 1) 
!            if(master) write(*,*) 'mean ice_mp=', tot_prec
        endif
        ice0(:,:) = ice0(:,:) + ice(:,:)
   endif

   if ( id_rain>0 ) then
        used=send_data(id_rain,    rain,    time)
        if ( seconds==0 ) then
!            tot_prec = g_sum(rain, is, ie, js, je, ng, area, 1) 
!            if(master) write(*,*) 'mean rain=', tot_prec
        endif
        rain0(:,:) = rain0(:,:) + rain(:,:)
   endif
   

   if ( id_prec>0 ) then
        used=send_data(id_prec, prec_mp, time)
   endif

!----------------------------------------------------------------------------

        prec0(:,:) = prec0(:,:) + prec_mp(:,:)
        prec1(:,:) = prec1(:,:) + prec_mp(:,:)
        mp_count = mp_count + 1.

        if ( seconds==0 .and. mp_print ) then
             tot_prec = g_sum(prec1*dt_in/86400., is, ie, js, je, ng, area, 1) 
             if(master) write(*,*) 'Daily prec_mp=', tot_prec
             prec1(:,:) = 0.
        endif
!----------------------------------------------------------------------------


  if ( mp_debug ) then
       call prt_maxmin('T_a_mp',    pt, is, ie, js, je, 0, kbot, 1., master)
       call prt_maxmin('qg_dt_a_mp',  qg_dt, is, ie, js, je, 0, kbot, 1., master)
       call prt_maxmin('prec', prec_mp, is, ie, js, je, 0,    1, 1., master)
  endif

                                        call timing_off(" split_mp")

 end subroutine micro_phys_driver



 subroutine mpdrv( delp, pt, p3, qv, ql, qr, qi, qs, qg, qa, dz,     &
                    is, ie, js, je, ks, ke, ktop, kbot, j, dt_in, ntimes,  & 
                    rain, snow, graupel, ice, &
                    cond, land, pt_dt, qv_dt, ql_dt, qr_dt, qi_dt,    &
                    qs_dt, qg_dt, qa_dt )

!-------------------------------------------------------------------
!  lin et al., 1983, jam, 1065-1092, and
!  rutledge and hobbs, 1984, jas, 2949-2972
!-------------------------------------------------------------------
! terminal fall is handled lagrangianly by conservative fv algorithm
!
! pt: temperature (k)
! 6 water species:
! 1) qv: water vapor (kg/kg)
! 2) ql: cloud water (kg/kg)
! 3) qr: rain        (kg/kg)
! 4) qi: cloud ice   (kg/kg)
! 5) qs: snow        (kg/kg)
! 6) qg: graupel     (kg/kg)

  integer,         intent(in):: j, is,ie, js,je, ks,ke
  integer,         intent(in):: ntimes, ktop, kbot
  real,            intent(in):: dt_in

  real, intent(in), dimension(is:ie,js:je,ks:ke) :: p3, delp
  real, intent(in), dimension(is:ie):: land  !land fraction
  real, intent(in   ), dimension(is:ie,js:je,ks:ke):: pt, qv, ql, qr, qi, qs, qg, qa
  real, intent(inout), dimension(is:ie,js:je,ks:ke):: pt_dt,  qa_dt, dz,   &
                                   qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt
  real, intent(out), dimension(is:ie):: rain, snow, ice, graupel, cond
!----------
! local var
!----------
  real, dimension(ktop:kbot):: qvz, qlz, qrz, &
                               qiz, qsz, qgz, qaz, &
                               vtiz, vtsz, vtgz, vtrz, &
                               dp1, qv0, ql0, qr0, qi0, qs0, qg0, qa0, t0, den, &
                               tz, p1, dz0, dz1, rhofac, rh_p, rh
  real :: r1, s1, i1, g1, rdt, omq, qaz1
  real :: qpz, clouds, ts, esk, adj_fra
  real :: dt5, t_fac, p_fac, rh_t1, rh_t2
  integer :: i,k,n

   rdt = 1. / dt_in

   dt5 = 0.5*dts

   do 2000 i=is, ie

   do k=ktop, kbot
      dp1(k) = delp(i,j,k)
       t0(k) = pt(i,j,k)
       tz(k) = t0(k) 
       p1(k) = p3(i,j,k)
!-----------------------------------
      qvz(k) = max(qvmin, qv(i,j,k))
!-----------------------------------
      qlz(k) = ql(i,j,k)
      qrz(k) = qr(i,j,k)
      qiz(k) = qi(i,j,k)
      qsz(k) = qs(i,j,k)
      qgz(k) = qg(i,j,k)
      qa0(k) = qa(i,j,k)
!------------------------------
! convert to dry mixing ratios:
!------------------------------
         omq = 1. - (qvz(k)+qlz(k)+qrz(k)+qiz(k)+qsz(k)+qgz(k))
      dp1(k) = dp1(k) * omq
         omq = 1. / omq
      qvz(k) = qvz(k)*omq
      qv0(k) = qvz(k)
      qlz(k) = qlz(k)*omq
      ql0(k) = qlz(k)
      qrz(k) = qrz(k)*omq
      qr0(k) = qrz(k)
      qiz(k) = qiz(k)*omq
      qi0(k) = qiz(k)
      qsz(k) = qsz(k)*omq
      qs0(k) = qsz(k)
      qgz(k) = qgz(k)*omq
      qg0(k) = qgz(k)
      qaz(k) = 0.  
      dz0(k) = dz(i,j,k)
      dz1(k) = dz0(k)
   enddo

     ql0_crt =  ql0_land*land(i) + ql0_ocean*(1.-land(i))
     c_praut = c_praut_l*land(i) + c_praut_o*(1.-land(i))

   if ( use_cld_prof ) then
        t_fac = max(0.,min(1.,(tz(kbot)-tice)/(30.)))
        !For land:
        rh_t1 = 1. - (1.-rh_l)*t_fac
        !For ocean
        rh_t2 = 1. - (1.-rh_o)*t_fac
   endif

   do k=ktop, kbot
      if ( use_cld_prof ) then
           if ( p1(k) > 900.e2 ) then 
                rh_p(k) = 1.
           elseif ( p1(k) > 800.e2 ) then 
                p_fac= (p1(k) - 800.e2) / 100.e2
                p_fac= p_fac ** 2
                rh_p(k) = (rh_t1 + (1.-rh_t1)*p_fac)*land(i) +  &
                          (rh_t2 + (1.-rh_t2)*p_fac)*(1.-land(i))
           else
                rh_p(k) = rh_t1*land(i) + rh_t2*(1.-land(i)) 
           endif
      else
           rh_p(k) = 1.
      endif
   enddo

   do k=ktop, kbot
         den(k) = p1(k)/(rair*tz(k))
      rhofac(k) = sqrt(sfcrho/den(k))
   enddo

 do 1000 n=1,ntimes

   call neg_adj(ktop, kbot, p1, tz, dp1, qvz, qlz, qrz, qiz, qsz, qgz)

   call fall_speed(ktop, kbot, den, qsz, qiz, qgz, qrz, vtsz, vtiz, vtgz, vtrz)

   call terminal_fall (dts, ktop, kbot, tz, qvz, qlz, qrz, qgz, qsz, qiz, p1, &
                       dz1, dp1, den, rhofac, vtrz, vtgz, vtsz, vtiz, rh_p,   &
                       r1, g1, s1, i1)

      rain(i) = rain(i)    + r1
      snow(i) = snow(i)    + s1
   graupel(i) = graupel(i) + g1
       ice(i) = ice(i)     + i1

   adj_fra = 0.5

   call warm_rain(dts, ktop, kbot, tz, qvz, qlz, qrz, p1, den, rhofac, rh_p, adj_fra)

!-------------------------------------------------
! * ice-phase microphysics & saturation adjustment
!-------------------------------------------------

   call icloud( ktop, kbot, tz, p1, qvz, qlz, qrz, qiz, qsz, qgz,  &
                den, rhofac, vtsz, vtgz, vtrz, rh_p, rh )

! Fractional clouds are only possible in sub-cycle mode
   do k=ktop, kbot
      if( rh(k) >= rh_p(k) ) then
          qaz(k) = qaz(k) + 1.
      elseif ( p1(k)<400.e2 .and. qiz(k)>1.e-4 ) then
! High ice clouds:
          qaz(k) = qaz(k) + dim(rh(k), 0.5)
      endif
   enddo

! update delta-height:
   do k=ktop, kbot
      dz1(k) = dz0(k)*tz(k)/t0(k) 
   enddo
 
1000  continue

   do k = ktop, kbot
      pt_dt(i,j,k) = pt_dt(i,j,k) + rdt*(tz(k)- t0(k))
               omq = dp1(k) / delp(i,j,k)
      qv_dt(i,j,k) = qv_dt(i,j,k) + rdt*(qvz(k)-qv0(k))*omq
      ql_dt(i,j,k) = ql_dt(i,j,k) + rdt*(qlz(k)-ql0(k))*omq
      qr_dt(i,j,k) = qr_dt(i,j,k) + rdt*(qrz(k)-qr0(k))*omq
      qi_dt(i,j,k) = qi_dt(i,j,k) + rdt*(qiz(k)-qi0(k))*omq
      qs_dt(i,j,k) = qs_dt(i,j,k) + rdt*(qsz(k)-qs0(k))*omq
      qg_dt(i,j,k) = qg_dt(i,j,k) + rdt*(qgz(k)-qg0(k))*omq
      qa_dt(i,j,k) = qa_dt(i,j,k) + rdt*(qaz(k)/real(ntimes)-qa0(k))
         dz(i,j,k) = dz1(k)
   enddo

!-----------------
! fms diagnostics:
!-----------------
   if ( id_cond>0 ) then
     do k=ktop,kbot                   ! total condensate
        cond(i) = cond(i) + dp1(k)*(qlz(k)+qrz(k)+qsz(k)+qiz(k)+qgz(k))
     enddo
   endif

   if ( id_vtr> 0 ) then
        do k=ktop, kbot
           vt_r(i,j,k) = vtrz(k)
        enddo
   endif
   if ( id_vts> 0 ) then
        do k=ktop, kbot
           vt_s(i,j,k) = vtsz(k)
        enddo
   endif
   if ( id_vtg> 0 ) then
        do k=ktop, kbot
           vt_g(i,j,k) = vtgz(k)
        enddo
   endif
   if ( id_vts> 0 ) then
        do k=ktop, kbot
           vt_i(i,j,k) = vtiz(k)
        enddo
   endif

2000  continue

 end subroutine mpdrv


 subroutine warm_rain( dt, ktop, kbot, tz, qv, ql, qr, pm,  &
                       den, denfac, rh_p, adj_fra )

 real,    intent(in):: dt                    ! time step (s)
 real,    intent(in):: adj_fra
 integer, intent(in):: ktop, kbot
 real,    intent(in), dimension(ktop:kbot):: pm, rh_p, den, denfac
 real,    intent(inout), dimension(ktop:kbot):: tz, qv, ql, qr
! local:
 real:: lcpk(ktop:kbot)
 real:: qsat, dqsdt, wsat, dwsdt, evap, cond, factor, qden, tsq, sink
 real:: clouds, qpz, ts, rh
 integer k

  do k=ktop,kbot
     lcpk(k) = latv / (cp - rdgas*ptop/pm(k))  ! Total energy conserving 
  enddo

 if ( auto_kk2000 ) then
! * KK2000 form
      raut_rate = dt*c_praut

  do k=ktop,kbot
         ts = tz(k) - lcpk(k)*ql(k)
       qsat = qs1d(ts, pm(k), dqsdt)   ! pure water table
        qpz = qv(k) + ql(k)
         rh = qpz/qsat
     if ( rh<rh_p(k) ) then
          tz(k) = ts
          qv(k) = qpz
          ql(k) = 0.    ! evaporate all cloud water before rain evap
          if ( qr(k) > qrmin ) then
               tsq = tz(k)**2
              qden = qr(k)*den(k)
              evap = crevp(1)*tsq*dim(rh_p(k)*qsat, qv(k)) *     &
                    (crevp(2)*sqrt(qden)+crevp(3)*qden**0.725)  &
                   /(crevp(4)*tsq + crevp(5)*qsat*den(k)) * evap_eff
              evap = min( qr(k), dt*evap,                &
                    evap_rhc*dim(rh_p(k)*qsat, qv(k))/(1.+lcpk(k)*dqsdt) )
             qv(k) = qv(k) + evap
             qr(k) = qr(k) - evap
             tz(k) = tz(k) - evap*lcpk(k)
          endif
     else
#ifdef TEST_ADJ
! partial temperature range adjustment (full adjustment done in icloud)
         if( tz(k) > tice ) then
             wsat = ws1d(tz(k), pm(k), dwsdt)
             cond = adj_fra*(qv(k)-rh_p(k)*wsat)/(1.+lcpk(k)*dwsdt)
             cond = max( -ql(k), cond )  ! prevent negative ql
            ql(k) = ql(k) + cond
            qv(k) = qv(k) - cond
            tz(k) = tz(k) + cond*lcpk(k)
         endif
#endif

         if ( ql(k) > qrmin ) then
! * Accretion
             factor = dt*denfac(k)*cracw*(max(qrmin,qr(k))*den(k))**0.95
               sink = factor/(1.+factor)*ql(k)
!            sink = min(ql(k), dt*67.*(ql(k)*max(qrmin,qr(k)))**1.15)
! * autoconversion
             if ( ql(k) > ql0_crt )   &
                  sink = min(ql(k), sink + raut_rate*ql(k)**so3)
            ql(k) = ql(k) - sink
            qr(k) = qr(k) + sink
         endif
     endif

  enddo

  else

  raut_rate = dt*c_praut/(1.+dt*c_praut)

  do k=ktop,kbot
         ts = tz(k) - lcpk(k)*ql(k)
       qsat = qs1d(ts, pm(k), dqsdt)
        qpz = qv(k) + ql(k)
         rh = qpz/qsat
     if ( rh<rh_p(k) ) then
          tz(k) = ts
          qv(k) = qpz
          ql(k) = 0.
          if ( qr(k) > qrmin ) then
               tsq = tz(k)**2
              qden = qr(k)*den(k)
              evap = crevp(1)*tsq*dim(rh_p(k)*qsat, qv(k)) *     &
                    (crevp(2)*sqrt(qden)+crevp(3)*qden**0.725)  &
                   /(crevp(4)*tsq + crevp(5)*qsat*den(k)) * evap_eff
              evap = min( qr(k), dt*evap,                &
                    evap_rhc*dim(rh_p(k)*qsat, qv(k))/(1.+lcpk(k)*dqsdt) )
             qv(k) = qv(k) + evap
             qr(k) = qr(k) - evap
             tz(k) = tz(k) - evap*lcpk(k)
          endif
#ifdef TEST_ADJ
     elseif( tz(k) > tice ) then
! partial adjustment (full adjustment done in icloud)
          wsat = ws1d(tz(k), pm(k), dwsdt)
          cond = adj_fra*(qv(k)-rh_p(k)*wsat)/(1.+lcpk(k)*dwsdt)
          cond = max( -ql(k), cond )  ! prevent negative
          ql(k) = ql(k) + cond
          qv(k) = qv(k) - cond
          tz(k) = tz(k) + cond*lcpk(k)
#else
     else
#endif
          if ( ql(k) > qrmin ) then
! * Accretion & autoconversion
             factor = dt*denfac(k)*cracw*(max(qrmin, qr(k))*den(k))**0.95
               sink = min(ql(k), factor/(1.+factor)*ql(k)+raut_rate*dim(ql(k), ql0_crt))
              ql(k) = ql(k) - sink
              qr(k) = qr(k) + sink
          endif
     endif
  enddo

 endif

 end subroutine warm_rain



 subroutine icloud(ktop, kbot, tzk, p1, qvk, qlk, qrk, qik, qsk, qgk, den,    &
                   denfac, vts, vtg, vtr, rh_p, rh)

!----------------------------------------------------
! Simplified cloud micro-physics; processes splitting
! with some un-split sub-grouping
! Time implicit (when possible) accretion and autoconversion
! Author: Shian-Jiann Lin, GFDL
!-------------------------------------------------------

 integer, intent(in) :: ktop, kbot
 real, intent(in), dimension(ktop:kbot):: p1, den, denfac,  &
                                          vts, vtg, vtr, rh_p
 real, intent(inout), dimension(ktop:kbot):: tzk, qvk, qlk, qrk, qik, qsk, qgk
 real, intent(out  ), dimension(ktop:kbot):: rh
! local:
 real, dimension(ktop:kbot) :: lcpk, icpk, tcpk
 real :: tz, qv, ql, qr, qi, qs, qg
 real :: c3racs(3), c3sacr(3), c3gacr(3), c3gacs(3)
 real :: praut, pracw, pracs, psacw, pgacw, pgmlt,   &
         psmlt, prevp, psacr, pgacr, pgfr,  pgacs,   &
         pgaut, pgaci, praci, psaut, psaci, pssub,   &
         pgsub, psfw,  psfi,  pidw,  piacr
 real :: tc, tsq, dqs0, si, qsi
 real :: factor, frez, sink, ratio
 real :: esk, qsat, qden
 real :: qsiz, qswz, qpz
 real :: clouds, rql, tmp1, tmp2, rqi, ts
 real :: ql0, qi0, fsat
 real :: dqsdt, dwsdt
 real, parameter:: t40 = tice - 40.
 integer, parameter:: n_qsat=12
 integer :: i, j, k, n
 real:: acco1(3,4)
 equivalence (c3racs(1),acco1(1,1)),(c3sacr(1),acco1(1,2)),  &
             (c3gacr(1),acco1(1,3)),(c3gacs(1),acco1(1,4))

 do j=1,4
    do i=1,3
       acco1(i,j) = acco(i,j)
    enddo
 enddo

 do k=ktop,kbot
    tmp1 = cp - rdgas*ptop/p1(k)
    lcpk(k) =  latv / tmp1
    icpk(k) =  lati / tmp1
    tcpk(k) = lcpk(k) + icpk(k)
 enddo

 do k=ktop, kbot
    if( tzk(k) < t40 ) then
! * pihom: freezing of cloud water
        frez = min(qlk(k), (tice-tzk(k))/icpk(k))
    elseif( tzk(k) > tice ) then
! * pimlt: melting of cloud ice 
        frez = -min(qik(k), (tzk(k)-tice)/icpk(k))
    else
       tmp1 = tice - tzk(k)
       frez = min(qlk(k), tmp1/icpk(k)) * min(1., tmp1*dts/(40.*tau_f))
    endif
    qlk(k) = qlk(k) - frez
    qik(k) = qik(k) + frez
    tzk(k) = tzk(k) + frez*icpk(k)
 enddo


 do 1000 k=ktop, kbot

   tz = tzk(k)
   qv = qvk(k)
   ql = qlk(k)
   qi = qik(k)
   qr = qrk(k)
   qs = qsk(k)
   qg = qgk(k)

  clouds = ql + qi
     qpz = qv + clouds
      ts = tz - (lcpk(k)*clouds + icpk(k)*qi)
     esk = es2_table(ts)     ! over pure ice if T<TICE
    qsat = eps*esk / max(esk, p1(k)-esk)
    rh(k) = qpz/qsat
       
   if ( rh(k)<rh_p(k) ) then
        tz = ts
        qv = qpz
        qi = 0.
        ql = 0.
        if ( (qr+qs+qg)<qrmin ) then
              tzk(k) = tz
              qvk(k) = qv
              qlk(k) = 0.
              qik(k) = 0.
              go to 1000
        endif
   endif

!--------------------------------------
! *** Split-micro_physics_processes ***
!--------------------------------------
! Zetac: excluded (from LFO83) term: psdep
! pgwet removed by SJL (reducing complexity)
! piacr added by SJL to increase ice production

   pgacr = 0.
   pgacw = 0.
   tc = tz-tice

if ( tc >= 0.0 ) then

     dqs0 = ces0/(p1(k)-es0) - qv

     if( qs>qrmin ) then ! melting of snow into rain
! * accretion: cloud water --> snow
! only rate is used (for snow melt) since tc > 0.
        if( ql>qrmin ) then
            factor = dts*denfac(k)*csacw*(qs*den(k))**0.8125
             psacw = rdts * factor/(1.+factor)*ql     ! rate
        else
             psacw = 0.
        endif

! * accretion: melted snow --> rain:
        if ( qr>qrmin ) then
             psacr = min(acr3d(vts(k), vtr(k), qr, qs, csacr, c3sacr, denfac(k), den(k)), qr*rdts)
        else
             psacr = 0.
        endif
! * Snow melt (due to rain accretion): snow --> rain
        psmlt = max(0., smlt(tc, dqs0, qs*den(k), psacw, psacr, csmlt, den(k), denfac(k)))

! * accretion: snow --> rain
        if( qr>qrmin ) then
            pracs = acr3d(vtr(k), vts(k), qs, qr, cracs, c3racs, denfac(k), den(k))
        else
            pracs = 0.
        endif

! Total snow sink:
        sink = min(qs, dts*(psmlt+pracs), tc/icpk(k))
        qs = qs - sink
        qr = qr + sink
        tz = tz - sink*icpk(k)
     endif   ! snow existed

#ifndef NO_GRAUPEL
     if ( qg>qrmin ) then
         if ( qr>qrmin ) then
! * accretion: rain --> graupel
              pgacr = min(acr3d(vtg(k), vtr(k), qr, qg, cgacr, c3gacr, denfac(k), den(k)), rdts*qr)
         endif

         if( ql>qrmin ) then
! * accretion: cloud water --> graupel
             factor = cgacw*(qg*den(k))**0.875 / sqrt(den(k))
              pgacw = factor/(1.+dts*factor) * ql  ! rate
         endif

! * melting: graupel --> rain
            tc = tz-tice
         pgmlt = dts*gmlt(tc, dqs0, qg*den(k), pgacw, pgacr, cgmlt, den(k))
         pgmlt = min( max(0., pgmlt), qg, max(0.,tc/icpk(k)) )
            qg = qg - pgmlt 
            qr = qr + pgmlt 
            tz = tz - pgmlt*icpk(k)
     endif   ! graupel existed
#endif

!----------------------------
else        ! tc < 0
!----------------------------

! Cloud water sink:
  if( ql>qrmin ) then

!#11 * pidw;     cloud water --> ice
      pidw = dts*realidw(tc, qi, den(k))   ! cld_wat --> cld_ice

!#5B * accretion: cloud water --> snow
      if( qs>qrmin ) then
          factor = dts*denfac(k)*csacw*(qs*den(k))**0.8125
          psacw = factor/(1.+factor)*ql
      else
          psacw = 0.
      endif

      sink = pidw + psacw    ! cloud water sink
      if ( sink>qrmin ) then
           ratio = pidw/sink
           sink = min(ql, sink, -tc/icpk(k))
           qi = qi + sink*ratio
           qs = qs + sink*(1.-ratio)
           ql = ql - sink
           tz = tz + sink*icpk(k)
           tc = tz - tice
      endif

   endif  ! cloud water existed
 
  if ( qi>qrmin ) then

!#9 * accretion: cloud ice --> snow
     if ( qs>qrmin )  then
! Note the c_psaci=0.1 (from Zetac):
! check Lin Eq. 22 has temperature dependency
          factor = dts*c_psaci*denfac(k)*csaci*(qs*den(k))**0.8125
           psaci = factor/(1.+factor) * qi
     else
           psaci = 0.
     endif

!#10 * autoconversion: cloud ice --> snow
! (according to LFO 1983: Eq. 21 solved implicitly)
     if ( qi>qi0_crt ) then
          factor = dts*c_psaut*exp(0.025*tc)
           psaut = factor/(1.+factor)*(qi-qi0_crt)
     else
           psaut = 0.
     endif
!------------------------------------------
! Bergeron processes -- psfw and psfi
!------------------------------------------
!#12: psfi: cloud ice   --> snow
!#13: psfw: cloud water --> snow

      call bergrn(tc, ql, qi, ql*den(k), psfw, psfi)
      psfw = min(dts*psfw, ql, max(0.,-tc/icpk(k)))
      psfi = min(dts*psfi, qi)

      sink = min(qi, psaci+psaut+psfi)
        qi = qi - sink
        qs = qs + sink + psfw
        ql = ql - psfw
        tz = tz + psfw*icpk(k)

#ifndef NO_GRAUPEL
!#15 * accretion: cloud ice --> graupel
! do graupel last!
      if ( qg>qrmin ) then
           factor = dts*cgaci*(qg*den(k))**0.875 / sqrt(den(k))
            pgaci = factor/(1.+factor)*qi
               qi = qi - pgaci
               qg = qg + pgaci
      endif
#endif

  endif  ! cloud ice existed
 
  if ( qr>qrmin ) then

! rain to ice, snow, graupel processes: it is ordered this way for better
! radiative balance in GFDL model.

!#26 * piacr: accretion of rain by cloud ice [lfo 26]
! rain --> ice  (factor somewhat arbitrary; totally tunable)
! The value of "factor" needs to be near order(1) to be really effective

       if ( qi>0. ) then   ! if ice existed: rain --> ice
            factor = dts*denfac(k)*qi * c_piacr  ! no effect if no cie
            piacr = factor/(1.+factor)*qr
            piacr = min( piacr, dim(tice,tz)/icpk(k) )
               qr = qr - piacr
               qi = qi + piacr
               tz = tz + piacr*icpk(k)
       endif

!#17 * accretion of rain by snow; rain --> snow
       if ( qs>0.0 .and. qr>0. ) then   ! if snow
            psacr = dts*acr3d(vts(k), vtr(k), qr, qs, csacr, c3sacr, denfac(k), den(k))
            psacr = min( psacr, qr,  dim(tice,tz)/icpk(k) )
               qr = qr - psacr
               qs = qs + psacr
               tz = tz + psacr*icpk(k)
       endif
 
!#14 * accretion: accretion of cloud ice by rain to produce snow or graupel
! (LFO: produces snow or graupel; cloud ice sink.. via psacr & pgfr)
       tc = tz-tice
       if ( qr>qrmin .and. tc<0. )  then
            factor = dts*denfac(k)*craci*(qr*den(k))**0.95
             praci = factor/(1.+factor)*qi   ! check praci form
#ifndef NO_GRAUPEL
             if ( qr > qr0_crt ) then
                  qg = qg + praci
             else
                  qs = qs + praci
             endif
#else
             qs = qs + praci
#endif
             qi = qi - praci
       endif

!#18 * freezing of remaining rain --> graupel/hail
       if ( qr>qrmin .and. tc<0. )  then
            pgfr = dts*cgfr(1)*(exp(-cgfr(2)*tc)-1.)*(qr*den(k))**1.75/den(k)
            pgfr = min(pgfr, -tc/icpk(k), qr)
              qr = qr - pgfr
#ifndef NO_GRAUPEL
              qg = qg + pgfr
#else
              qs = qs + pgfr
#endif
              tz = tz + pgfr*icpk(k)
       endif
  endif
 
  if ( qs > 0.0 ) then   ! snow <--> vapor
       esk = es2_table(tz)
       qsi = eps*esk / max(esk, p1(k)-esk)
        si = (qv+qi)/qsi
!#19 * Sublimation of snow:
       qden = qs*den(k)
        tsq = tz**2
      pssub = cssub(1)*tsq*qsi*(si-rh_p(k))*(cssub(2)*sqrt(qden) + &
              cssub(3)*(qden)**0.65625*sqrt(denfac(k)))/(cssub(4)*tsq+cssub(5)*qsi*den(k))
      pssub = min(dts*pssub, qv)
      pssub = max(pssub, -qs)
         qs = qs + pssub 
         qv = qv - pssub 
         tz = tz + pssub*tcpk(k)

#ifndef NO_GRAUPEL
    if ( s2g_auto .and. qs>qs0_crt ) then
!#21 * Snow --> graupel auto-conversion
         factor = dts * 1.e-3*exp(0.09*(tz-tice))
          pgaut = factor/(1.+factor)*(qs-qs0_crt)
             qs = qs - pgaut
             qg = qg + pgaut
    endif
#endif

  endif

 
#ifndef NO_GRAUPEL
! Graupel production terms:
  if ( qg>qrmin ) then
     if ( qr>qrmin ) then 
!#23 * accretion: rain --> graupel
          pgacr = min(dts*acr3d(vtg(k), vtr(k), qr, qg, cgacr, c3gacr, denfac(k), den(k)), qr)
             qr = qr - pgacr
     endif

     if( ql>qrmin ) then
!#24 * accretion: cloud water --> graupel
         factor = dts*cgacw*(qg*den(k))**0.875 / sqrt(den(k))
          pgacw = factor/(1.+factor)*ql
             ql = ql - pgacw
     endif

!#20 * accretion: snow --> graupel
     if( s2g_accr .and. qs>qrmin  ) then
         pgacs = min(dts*acr3d(vtg(k), vts(k), qs, qg, cgacs*c_pgacs, c3gacs, denfac(k),den(k)), qs)
         qs = qs - pgacs
         qg = qg + pgacs
     endif

     qg = qg +  pgacr+pgacw
     tz = tz + (pgacr+pgacw)*icpk(k)

     esk = es2_table(tz)
     qsi = eps*esk / max(esk, p1(k)-esk)
     si = (qv+qi)/qsi
     if ( (ql+qi)<qrmin .and. si<rh_p(k) ) then
           qden = qg*den(k)
            tsq = tz**2
!#22 * graupel sublimation (unsaturated)
          pgsub = cgsub(1)*tsq*qsi*(rh_p(k)-si)*(cgsub(2)*sqrt(qden) + &
                  cgsub(3)*(qden**0.6875)/den(k)**0.25)/(cgsub(4)*tsq+cgsub(5)*qsi*den(k))
          pgsub = min(dts*pgsub, qg)
             qg = qg - pgsub
             qv = qv + pgsub 
             tz = tz - pgsub*tcpk(k)
     endif
 endif    ! graupel existed
#endif
endif   ! end if tc < 0


!----------------------
! saturation adjustment
!----------------------
  clouds = ql + qi
     qpz = qv + clouds

  if( clouds > 1.e-5 ) then
      rqi = qi/clouds
  else
      rqi = max(0., min(1., (tice-tz)/(tice-248.16)))  ! to -25 deg
  end if
  rql = 1. - rqi

   ql0 = ql;  qi0 = qi
  tmp1 = lcpk(k) + rqi*icpk(k)
  tmp2 = tz  + lcpk(k)*qv + icpk(k)*(rqi*qpz - qi)

  fsat = 1.;  ts = tz

  do n=1,n_qsat
     qsiz = qs1d(ts, p1(k), dqsdt)
     if ( ts<tice .and. ts>t40 ) then
          qswz = ws1d(ts, p1(k), dwsdt)
     else
          qswz = qsiz
         dwsdt = dqsdt   ! pure water or pure ice
     endif
     qsat = rql*qswz + rqi*qsiz
     if( abs(fsat) < 0.01 ) go to 555
     fsat = ts + tmp1*qsat - tmp2
      ts = ts - fsat / (1.+tmp1*(rql*dwsdt + rqi*dqsdt))
  enddo
  call error_mesg ('Micro_physics', 'qsat adjustment not converging', FATAL) 

555 continue

  rh(k) = qpz/qsat
  if( rh(k) >= rh_p(k) ) then
      qv = qsat * rh_p(k)
      ql = rql * (qpz-qv)
      qi = qpz - (qv+ql)
  else
      qi = qi * dim(rh(k), 0.5)   ! ice sink
      qv = qpz - qi
      ql = 0.
  end if

  tzk(k) = tz + lcpk(k)*(ql-ql0) + tcpk(k)*(qi-qi0)
  qvk(k) = qv
  qlk(k) = ql
  qik(k) = qi
  qrk(k) = qr
  qsk(k) = qs
  qgk(k) = qg

1000 continue

 end subroutine icloud



 subroutine terminal_fall(dtm, ktop, kbot, tz, qv, ql, qr, qg, qs, qi, pm, dz, dp,  &
                          den, denfac, vtr, vtg, vts, vti, rh_p, r1, g1, s1, i1)

! lagrangian control-volume method:

 real,    intent(in):: dtm                    ! time step (s)
 integer, intent(in):: ktop, kbot
 real,    intent(in), dimension(ktop:kbot):: dp, vtr, vtg, vts, vti, pm, rh_p,  &
                                             den, denfac
 real,    intent(inout), dimension(ktop:kbot):: dz, qv, ql, qr, qg, qs, qi, tz
 real,    intent(out):: r1, g1, s1, i1
! local:
 real, parameter:: dz_min = 1.e-2
 real, dimension(ktop:kbot+1):: ze, zt
 real:: qsat, dqsdt, dt5, melt, evap, dtime
 real:: factor, frac, frozen_p
 real:: qden, accr, tsq, tmp1
 real, dimension(ktop:kbot):: freezing_rain, l2r, lcpk, icpk
 real:: zs = 0.
 integer k, k0, m

  do k=ktop,kbot
        tmp1 = cp - rdgas*ptop/pm(k)
     lcpk(k) = latv / tmp1
     icpk(k) = lati / tmp1
  enddo

  dt5 = 0.5*dtm

  ze(kbot+1) = zs
  do k=kbot, ktop, -1
     ze(k) = ze(k+1) - dz(k)  ! dz<0
  enddo
  zt(ktop) = ze(ktop)


 do k=ktop+1,kbot
    zt(k) = ze(k) - dt5*(vtr(k-1)+vtr(k))
 enddo
 zt(kbot+1) = zs - dtm*vtr(kbot)

 do k=ktop,kbot
    if( zt(k+1)>=zt(k) ) zt(k+1) = zt(k) - dz_min
    freezing_rain(k) = 0.
    l2r(k) = 0.
 enddo


 frozen_p = 0.

 if ( falling_rain_mp ) then
  do k=kbot-1,ktop,-1
     if ( qr(k) > qrmin ) then
          do m=k+1, kbot
             if ( zt(k+1)>=ze(m) ) exit
!            if ( zt(k)<ze(m) ) then    ! the top of the c-v is in the layer below
             if ( zt(k)<ze(m+1) ) then
! In case of small time step the following process should have no effect
! There is possibility of double count; but it is better than no count at all.
                dtime = min( dtm, (ze(m)-ze(m+1))/vtr(k) )
                 qden = qr(k)*den(k)
                if ( ql(m)>qrmin ) then
#ifndef TEST_ACC
                     factor = dtime*denfac(m)*cracw*qden**0.95
! Accretion: linear form of KK2000: 3.7*qr*qr; solved implicitly
!                    factor = 3.7*dtime*qr(k)*den(k)/den(m)
                       accr = factor/(1.+factor)*ql(m)
                      ql(m) =  ql(m) - accr
                     l2r(m) = l2r(m) + accr
#endif
                else
! * no clouds * evaporation of falling rain:
                    qsat = qs1d(tz(m), pm(m), dqsdt)
                    evap = dim( rh_p(m)*qsat, qv(m) ) *             &
                           min( dtime/(tau_evap*2.0e-2)*sqrt(qden)*(tz(m)/300.),   &
                                0.9/(1.+lcpk(m)*dqsdt) )
                    evap = min ( qr(k)*dp(k)/dp(m), evap )
                   qv(m) = qv(m) + evap
                   tz(m) = tz(m) - evap*lcpk(m)
                   qr(k) = qr(k) - evap*dp(m)/dp(k)
                endif

                if ( tz(m)<tice ) then 
                     dtime = min(1., (ze(m)-ze(m+1))/(vtr(k)*tau_r) )
                      melt = min(qr(k)*dp(k)/dp(m), dtime*dim(tice, tz(m))/icpk(m))
                     tz(m) = tz(m) + melt*icpk(m)
                     qr(k) = qr(k) - melt*dp(m)/dp(k)
                     if ( zt(k)<zs ) then
                          frozen_p = frozen_p + melt*dp(m)
                     else
                          freezing_rain(m) = freezing_rain(m) + melt
                     endif
                endif
             endif
             if ( qr(k) < qrmin ) exit 
          enddo
     endif
  enddo
 endif

 call lagrangian_fall(ktop, kbot, zs, ze, zt, dp, qr, r1)

#ifndef TEST_ACC
 do k=ktop,kbot
    qr(k) = qr(k) + l2r(k)
 enddo
#endif

! find freezing level
  k0 = kbot
  do k=ktop, kbot-1
     if ( tz(k) > tice ) then
          k0 = k
          go to 11
     endif
  enddo
11  continue

!-----
! ice:
!-----
  if ( vi_fac < 1.e-5 ) then
     i1 = 0.
  else
  do k=ktop+1,kbot
     zt(k) = ze(k) - dt5*(vti(k-1)+vti(k))
  enddo
  zt(kbot+1) = zs - dtm*vti(kbot)

  do k=ktop,kbot
     if( zt(k+1)>=zt(k) ) zt(k+1) = zt(k) - dz_min
  enddo

  if ( k0 < kbot ) then
  do k=kbot-1,k0,-1
     if ( qi(k) > qrmin ) then
          do m=k+1, kbot
             if ( zt(k+1)>=ze(m) ) exit
             if ( zt(k)<ze(m) ) then    ! the top of the c-v is in the layer below
                  dtime = min(1., (ze(m)-ze(m+1)) / (max(vmin,vti(k))*tau_i))
                   melt = min(qi(k)*dp(k)/dp(m), dtime*dim(tz(m),tice)/icpk(m))
                  ql(m) = ql(m) + melt      ! melt into local cloud water
                  tz(m) = tz(m) - melt*icpk(m)
                  qi(k) = qi(k) - melt*dp(m)/dp(k)
             endif
          enddo
     endif
  enddo
  endif
  call lagrangian_fall(ktop, kbot, zs, ze, zt, dp, qi, i1)
  endif

!--------------------------------------------
! melting of falling snow (qs) into rain(qr)
!--------------------------------------------
  do k=ktop+1,kbot
     zt(k) = ze(k) - dt5*(vts(k-1)+vts(k))
  enddo
  zt(kbot+1) = zs - dtm*vts(kbot)

  do k=ktop,kbot
     if( zt(k+1)>=zt(k) ) zt(k+1) = zt(k) - dz_min
  enddo

  if ( k0 < kbot ) then
  do k=kbot-1,k0,-1
     if ( qs(k) > qrmin ) then
          do m=k+1, kbot
             if ( zt(k+1)>=ze(m) ) exit
             dtime = min( dtm, (ze(m)-ze(m+1))/(vmin+vts(k)) )
             if ( zt(k)<ze(m+1) ) then
!            if ( zt(k)<ze(m) ) then    ! the top of the c-v is in the layer below
                  dtime = min(1., dtime/tau_s)
                   melt = min(qs(k)*dp(k)/dp(m), dtime*dim(tz(m),tice)/icpk(m))
                  tz(m) = tz(m) - melt*icpk(m)
                  qs(k) = qs(k) - melt*dp(m)/dp(k)
                  if ( zt(k)<zs ) then
                       r1 = r1 + melt*dp(m)
                  else
                       qr(m) = qr(m) + melt   ! qr source here will fall next time step
                  endif
             endif
             if ( qs(k) < qrmin ) exit
          enddo
     endif
  enddo
  endif
  call lagrangian_fall(ktop, kbot, zs, ze, zt, dp, qs, s1)

#ifndef NO_GRAUPEL
!----------------------------------------------
! melting of falling graupel (qg) into rain(qr)
!----------------------------------------------
  do k=ktop+1,kbot
     zt(k) = ze(k) - dt5*(vtg(k-1)+vtg(k))
  enddo
  zt(kbot+1) = zs - dtm*vtg(kbot)

  do k=ktop,kbot
     if( zt(k+1)>=zt(k) ) zt(k+1) = zt(k) - dz_min
  enddo

  if ( k0 < kbot ) then
  do k=kbot-1,k0,-1
     if ( qg(k) > qrmin ) then
          do m=k+1, kbot
             if ( zt(k+1)>=ze(m) ) exit
             dtime = min( dtm, (ze(m)-ze(m+1))/vtg(k) )
             if ( zt(k)<ze(m+1) ) then
!            if ( zt(k)<ze(m) ) then    ! the top of the c-v is in the layer below
                  dtime = min(1., dtime/tau_g)
                   melt = min(qg(k)*dp(k)/dp(m), dtime*dim(tz(m),tice)/icpk(m))
                  tz(m) = tz(m) - melt*icpk(m)
                  qg(k) = qg(k) -  melt*dp(m)/dp(k)
                  if ( zt(k)<zs ) then
                       r1 = r1 + melt*dp(m)
                  else
                       qr(m) = qr(m) + melt
                  endif
             endif
             if ( qg(k) < qrmin ) exit
           enddo
     endif
  enddo
  endif
  call lagrangian_fall(ktop, kbot, zs, ze, zt, dp, qg, g1)

  g1 = g1 + frozen_p
  do k=ktop,kbot
     qg(k) = qg(k) + freezing_rain(k)
  enddo

#else
  g1 = 0.
  s1 = s1 + frozen_p
  do k=ktop,kbot
     qs(k) = qs(k) + freezing_rain(k)
  enddo
#endif


 end subroutine terminal_fall



 subroutine lagrangian_rain(ktop, kbot, zs, ze, zt, dp, qr, ql, precip)
!
! Perform precip and accretion/evaporation
!
 real,    intent(in):: zs
 integer, intent(in):: ktop, kbot
 real,    intent(in), dimension(ktop:kbot):: dp
 real,    intent(in), dimension(ktop:kbot+1):: ze, zt
 real,    intent(inout), dimension(ktop:kbot):: qr, ql
 real,    intent(out):: precip
! local:
 real, dimension(ktop:kbot):: rden, qm1, qm2
 integer k, k0, n, m

! density of rain
  do k=ktop,kbot
     rden(k) = qr(k)*dp(k)
      qm1(k) = rden(k) / (zt(k)-zt(k+1))
  enddo

  k0 = ktop

  do k=kbot, ktop, -1
 
     if ( zt(k)<=zs ) exit

     do m=k0,kbot

! Search where the top edge of the cv is located
        if ( zt(k)<ze(m) .and. zt(k)>ze(m+1) ) then

             if( zt(k+1)>=ze(m+1) ) then
             endif

! Search where the bottom of the cv is located
             do n=m,kbot
                if( zt(k+1)>ze(n) ) then
                endif               
             enddo

        endif

     enddo
  enddo

  precip = 0.

! Determine precip
  do k=ktop,kbot
     if ( zt(k+1)<zs ) then
          precip = qr(k)*dp(k)*(zs-zt(k+1))/(zt(k)-zt(k+1)) 
          if ( k/=kbot ) then
             do m=k+1,kbot
                precip = precip + qr(m)*dp(m)
             enddo
          endif
     endif
  enddo

1000 continue

  do k=ktop,kbot
     qr(k) = qm2(k) / dp(k)
  enddo

 end subroutine lagrangian_rain



 subroutine lagrangian_fall(ktop, kbot, zs, ze, zt, dp, q, precip)
 real,    intent(in):: zs
 integer, intent(in):: ktop, kbot
 real,    intent(in), dimension(ktop:kbot):: dp
 real,    intent(in), dimension(ktop:kbot+1):: ze, zt
 real,    intent(inout), dimension(ktop:kbot):: q
 real,    intent(out):: precip
! local:
 real, dimension(ktop:kbot):: qm1, qm2
 integer k, k0, n, m

! density:
  do k=ktop,kbot
     qm1(k) = q(k)*dp(k) / (zt(k)-zt(k+1))
  enddo

   k0 = ktop
   do k=ktop,kbot
      do n=k0,kbot
      if(ze(k) <= zt(n) .and. ze(k) >= zt(n+1)) then
         if(ze(k+1) >= zt(n+1)) then
!                          entire new grid is within the original grid
            qm2(k) = qm1(n)*(ze(k)-ze(k+1))
            k0 = n
            goto 555
         else
            qm2(k) = qm1(n)*(ze(k)-zt(n+1))    ! fractional area
            do m=n+1,kbot
!                                        locate the bottom edge: ze(k+1)
               if(ze(k+1) < zt(m+1) ) then
                  qm2(k) = qm2(k) + q(m)*dp(m)
               else
                  qm2(k) = qm2(k) + qm1(m)*(zt(m)-ze(k+1))
                  k0 = m
                  goto 555
               endif
            enddo
            goto 555
         endif
      endif
      enddo
555 continue
   enddo

     precip = 0.
! direct algorithm (prevent small negatives)
     do k=ktop,kbot
        if ( zt(k+1) < zs ) then
             precip = qm1(k)*(zs-zt(k+1)) 
             if ( (k+1) > kbot ) goto 777
                  do m=k+1,kbot
                     precip = precip + q(m)*dp(m)
                  enddo
             goto 777
        endif
     enddo
777  continue

   do k=ktop,kbot
      q(k) = qm2(k) / dp(k)
   enddo

 end subroutine lagrangian_fall


 subroutine fall_speed(ktop, kbot, den, qs, qi, qg, qr, vts, vti, vtg, vtr)
 integer, intent(in)                     :: ktop, kbot
 real, intent(in ), dimension(ktop:kbot) :: den, qs, qi, qg, qr
 real, intent(out), dimension(ktop:kbot) :: vts, vti, vtg, vtr
! fall velocity constants:
 real :: vconr = 2503.23638966667,   vcons = 6.6280504,        &
         vcong = 87.2382675,         vconi = 3.29
 real, parameter :: thr=1.e-9, thg = 1.e-9, ths=1.e-9, thi=1.e-9
 real:: rhof, rhos
 integer:: k
!-----------------------------------------------------------------------
! marshall-palmer constants
!-----------------------------------------------------------------------
 real :: normr = 25132741228.7183,   norms = 942477796.076938, &
         normg =  5026548245.74367

!-----------------------------------------------------------------------
! marshall-palmer formula
!-----------------------------------------------------------------------

! sjl: try the local air density -- for global model; the true value could be
! much smaller than sfcrho over high mountains
!  rhos = den(kbot) 
   rhos = 1.2

   do k=ktop, kbot
      rhof = sqrt( min(50., max(1., rhos/den(k))) )
! rain:
      vtr(k) = vr_fac*max(vmin, vconr*rhof*(dim(qr(k)*den(k),thr)/normr)**0.2)
! snow:
      vts(k) = vs_fac*vcons*rhof*(dim(qs(k)*den(k),ths)/norms)**0.0625

#ifndef NO_GRAUPEL
! graupel:
      vtg(k) = vg_fac*max(vmin, vcong*rhof*(dim(qg(k)*den(k),thg)/normg)**0.125)
#endif

! ice:
      vti(k) = vi_fac*vconi*rhof*dim(qi(k)*den(k),thi)**0.16
   enddo

 end subroutine fall_speed

 real function realidw(tc,qik,rho )
      real, intent(in):: tc, qik,  rho
      real :: fnc, fmi, a1, a2, tc1, a1t, a2t, rminuc
      dimension a1t(0:31),a2t(0:31)
! fu *******************************************************************
!
!     note -- a2t is identical to a2t in subroutine bergrn, but a1t
!     is multiplied by a factor of 1.e-5.
!
! fu *******************************************************************
      data a1t/.0,.7939e-12,.7841e-11,.3369e-10,0.4336e-10,0.5285e-10, &
              0.3728e-10,0.1852e-10,0.2991e-11,0.4248e-11,0.7434e-11,&
              0.1812e-10,0.4394e-10,0.9145e-10,0.1725e-9 ,0.3348e-9 ,&
              0.1725e-9 ,0.9175e-10,0.4412e-10,0.2252e-10,0.9115e-11,&
              0.4876e-11,0.3473e-11,0.4758e-11,0.6306e-11,0.8573e-11,&
              0.7868e-11,0.7192e-11,0.6513e-11,0.5956e-11,0.5333e-11,&
              0.4834e-11/
      data a2t/0.0,0.4006,0.4831,0.5320,0.5307,0.5319, &
              0.5249,0.4888,0.3894,0.4047,0.4318, &
              0.4771,0.5183,0.5463,0.5651,0.5813, &
              0.5655,0.5478,0.5203,0.4906,0.4447, &
              0.4126,0.3960,0.4149,0.4320,0.4506, &
              0.4483,0.4460,0.4433,0.4413,0.4382, &
              0.4361/
      data rminuc/1.05e-15/

      tc1 = max(tc,-30.0)
      a1  = (a1t(-int(tc1))-a1t(-int(tc1)+1))*(tc1-int(tc1)+1.0)+ &
           a1t(-int(tc1)+1)
      a2  = (a2t(-int(tc1))-a2t(-int(tc1)+1))*(tc1-int(tc1)+1.0)+ &
           a2t(-int(tc1)+1)
      fnc = 0.01 * exp ( - 0.6 * tc )
      fmi = rho * qik / fnc * 1000.0
      realidw = exp(-0.6*tc)*a1*fmi**a2/rho

 end function realidw



 subroutine setupm

 real:: gcon, cd, scm3, pisq, act(8), acc(3)
 real :: vdifu, tcond
 real :: rdrya,visk
 real :: cwacs, ch2o, hltf
 real ::  hlts, hltc, ri50

 real :: gam263, gam275, gam290,                                &
         gam325, gam350, gam380,                                &
         gam425, gam450, gam480,                                &
         gam625, gam680

 data  gam263/1.456943/,   gam275/1.608355/,  gam290/1.827363/  &
       gam325/2.54925/,    gam350/3.323363/,  gam380/4.694155/  &
       gam425/8.285063/,   gam450/11.631769/, gam480/17.837789/ &
       gam625/184.860962/, gam680/496.604067/
!
!     physical constants (mks)
!     lin's constants(mks) except rmi50,rmi40 (cgs)
!
 real :: alin, clin
 real :: rnzr, rnzs, rnzg, rhor, rhos, rhog
 data alin, clin  /842.0, 4.80/
 data rnzr /8.0e6/  ! lin83
 data rnzs /3.0e6/  ! lin83
 data rnzg /4.0e6/  ! rh84
 data rhor /1.0e3/  ! lin83
 data rhos /0.1e3/  ! lin83
 data rhog /0.4e3/  ! rh84
 data acc/5.0,2.0,0.5/

 integer :: k, i

      pie = 4.*atan(1.0)

      vdifu=2.11e-5
      tcond=2.36e-2

!     rdrya=2.87e2

      visk=1.259e-5
      hlts=2.8336e6
      hltc=2.5e6
      hltf=3.336e5

      ch2o=4.1855e3
      rmi50=3.84e-6
      rmi40=2.46e-7
      ri50=1.e-4

!     parameters passed from model which uses adams-bashforth scheme
!     cplc = hltc/cp
!     cpls = hlts/cp
!     cplf = hltf/cp
      pisq = pie*pie
      scm3 = (visk/vdifu)**(1./3.)
!
!     acr3:  four lead constants required, three for each summation
!            four separate processes:  racs,sacr,gacr,gacs
!
      cracs = pisq*rnzr*rnzs*rhos
      csacr = pisq*rnzr*rnzs*rhor
      cgacr = pisq*rnzr*rnzg*rhor
      cgacs = pisq*rnzg*rnzs*rhos
!
!     act:  1-2:racs(s-r); 3-4:sacr(r-s);
!           5-6:gacr(r-g); 7-8:gacs(s-g)
!
      act(1) = pie * rnzs * rhos
      act(2) = pie * rnzr * rhor
      act(6) = pie * rnzg * rhog
      act(3) = act(2)
      act(4) = act(1)
      act(5) = act(2)
      act(7) = act(1)
      act(8) = act(6)

      do i=1,3
         do k=1,4
            acco(i,k) = acc(i)/(act(2*k-1)**((7-i)*0.25)*act(2*k)**(i*0.25))
         enddo
      enddo
!
      gcon  = 40.74 * sqrt( sfcrho )   ! 44.628
!
!     acr1:  single constant required
!     five separate processes:  sacw,wacs,iacr,raci,saci
!
      csacw = pie*rnzs*clin*gam325/(4.*act(1)**0.8125)
!     cwacs = pisq*rhos*rnzs*clin*gam625/(1.0056e-10*act(1)**1.5625)
!     ciacr = pisq*rhor*rnzr*alin*gam680/(1.0056e-11*act(2)**1.7)
      craci = pie*rnzr*alin*gam380/(4.*act(2)**0.95)
      csaci = csacw
!
!     cagci is for dry growth
!
      cgacw = pie*rnzg*gam350*gcon/(4.*act(6)**0.875)
      cgaci = cgacw*0.1
!
!     racw
!
      cracw = craci            ! cracw= 3.27206196043822
!
!     subl and revp:  five constants for three separate processes
!
      cssub(1) = 2.*pie*vdifu*tcond*rvapr*rnzs
      cgsub(1) = 2.*pie*vdifu*tcond*rvapr*rnzg
      crevp(1) = 2.*pie*vdifu*tcond*rvapr*rnzr
      cssub(2) = 0.78/sqrt(act(1))
      cgsub(2) = 0.78/sqrt(act(6))
      crevp(2) = 0.78/sqrt(act(2))
      cssub(3) = 0.31*scm3*gam263*sqrt(clin/visk)/act(1)**0.65625
      cgsub(3) = 0.31*scm3*gam275*sqrt(gcon/visk)/act(6)**0.6875
      crevp(3) = 0.31*scm3*gam290*sqrt(alin/visk)/act(2)**0.725
      cssub(4) = tcond*rvapr
      cssub(5) = hlts**2*vdifu
      cgsub(4) = cssub(4)
      crevp(4) = cssub(4)
      cgsub(5) = cssub(5)
      crevp(5) = hltc**2*vdifu
!
!     gfr:  two constants
!
      cgfr(1) = 20.e2*pisq*rnzr*rhor/act(2)**1.75
      cgfr(2) = 0.66
!
!sk   smlt:  five constants ( lin et al. 1983 )
!
!sk ********************************************************************
      csmlt(1) = 2.*pie*tcond*rnzs/hltf
      csmlt(2) = 2.*pie*vdifu*rnzs*hltc/hltf
      csmlt(3) = cssub(2)
      csmlt(4) = cssub(3)
      csmlt(5) = ch2o/hltf
!sk ********************************************************************
!
!     gmlt:  five constants
!
      cgmlt(1) = 2.*pie*tcond*rnzg/hltf
      cgmlt(2) = 2.*pie*vdifu*rnzg*hltc/hltf
!sk ********************************************************************
      cgmlt(3) = cgsub(2)
      cgmlt(4) = cgsub(3)
!sk ********************************************************************
      cgmlt(5) = ch2o/hltf
!
!     gwet:  two constants plus calc. of sat. vapor pressure at tc=0.0
!     ventilation coefficient of 10 included
!
!     cgwet(1) = cgmlt(1)
!     cgwet(2) = cgmlt(2)
!     cgwet(3) = cgmlt(3)
!     cgwet(4) = cgmlt(4)
      es0 = 6.107799961e2
      ces0 = eps*es0
!
!     bergrn: two constants
!     c2brg has conversion factor of 10**3
!
      c1brg = dts/rmi50
!lin  c2brg = ri50**2*1.e3 ! error
      c2brg = pie*ri50**2*1.e3
!
      do_setup = .false.

 end subroutine setupm


 subroutine micro_phys_init(axes, time)
    integer,         intent(in) :: axes(4)
    type(time_type), intent(in) :: time
    
    integer   :: unit, io, ierr
    logical   :: flag

    if( file_exist( 'input.nml' ) ) then
       unit = open_namelist_file ()
       io = 1
       do while ( io .ne. 0 )
          read( unit, nml = mp_lin_nml, iostat = io, end = 10 )
          ierr = check_nml_error(io,'mp_lin_nml')
       end do
10     call close_file ( unit )
    end if
    call write_version_number (version, tagname)
    write( stdlog(), nml = mp_lin_nml )
 
    id_vtr = register_diag_field ( mod_name, 'vt_r', axes(1:3), time,        &
         'rain fall speed', 'm/sec', missing_value=missing_value )
    id_vts = register_diag_field ( mod_name, 'vt_s', axes(1:3), time,        &
         'snow fall speed', 'm/sec', missing_value=missing_value )
    id_vtg = register_diag_field ( mod_name, 'vt_g', axes(1:3), time,        &
         'graupel fall speed', 'm/sec', missing_value=missing_value )
    id_vti = register_diag_field ( mod_name, 'vt_i', axes(1:3), time,        &
         'ice fall speed', 'm/sec', missing_value=missing_value )

    id_rain = register_diag_field ( mod_name, 'rain_lin', axes(1:2), time,        &
         'rain_lin', 'mm/day', missing_value=missing_value )
    id_snow = register_diag_field ( mod_name, 'snow_lin', axes(1:2), time,        &
         'snow_lin', 'mm/day', missing_value=missing_value )
    id_graupel = register_diag_field ( mod_name, 'graupel_lin', axes(1:2), time,  &
         'graupel_lin', 'mm/day', missing_value=missing_value )
    id_ice = register_diag_field ( mod_name, 'ice_lin', axes(1:2), time,        &
         'ice_lin', 'mm/day', missing_value=missing_value )
    id_prec = register_diag_field ( mod_name, 'prec_lin', axes(1:2), time,     &
         'prec_lin', 'mm/day', missing_value=missing_value )

    id_cond = register_diag_field ( mod_name, 'cond_lin', axes(1:2), time,     &
         'total condensate', 'kg/m**2', missing_value=missing_value )


    module_is_initialized = .true.
    
 end subroutine micro_phys_init



 subroutine micro_phys_end
 real gmp

  if ( mp_print ) then
! the g_sum call does not work if physics window is used *****
   if ( id_ice> 0 ) then
        gmp = g_sum(ice0, isc, iec, jsc, jec, ng, area, 1) 
        if(master) write(*,*) 'total ice=', gmp/mp_count
   endif
   if ( id_graupel> 0 ) then
        gmp = g_sum(graupel0, isc, iec, jsc, jec, ng, area, 1) 
        if(master) write(*,*) 'total graupel=', gmp/mp_count
   endif
   if ( id_snow> 0 ) then
        gmp = g_sum(snow0, isc, iec, jsc, jec, ng, area, 1) 
        if(master) write(*,*) 'total snow=', gmp/mp_count
   endif
   if ( id_rain> 0 ) then
        gmp = g_sum(rain0, isc, iec, jsc, jec, ng, area, 1) 
        if(master) write(*,*) 'total rain=', gmp/mp_count
   endif
!  if ( id_prec> 0 ) then
        gmp = g_sum(prec0, isc, iec, jsc, jec, ng, area, 1) 
        if(master) write(*,*) 'total prec=', gmp/mp_count
!  endif
  endif

   if ( id_vtr> 0 ) then
        deallocate ( vt_r )
   endif
   if ( id_vts> 0 ) then
        deallocate ( vt_s )
   endif
   if ( id_vti> 0 ) then
        deallocate ( vt_i )
   endif
   if ( id_vtg> 0 ) then
        deallocate ( vt_g )
   endif

   deallocate (  prec_mp  )
   deallocate (  prec0    )
   deallocate (  prec1    )
   deallocate (  rain0    )
   deallocate (  snow0    )
   deallocate (  ice0     )
   deallocate (  graupel0 )
   deallocate (  cond )

   deallocate ( table  )
   deallocate ( table2 )
   deallocate ( tablew )
   deallocate ( des )
   deallocate ( des2 )
   deallocate ( desw )
    
 end subroutine micro_phys_end



 subroutine setup_con( is, ie, js, je, ks, ke )
 integer, intent(in) :: is,ie, js,je, ks, ke

  master = gid==0

  isc = is;   iec = ie
  jsc = js;   jec = je

  lcp = latv / cp
  icp = lati / cp
  tcp = (latv+lati) / cp

  rgrav = 1./ grav

  call qsmith_init

! fall speed diagnostics:
      if ( id_vtr> 0 ) then
           allocate ( vt_r(is:ie, js:je, ks:ke) )
           vt_r = 0.
      endif
      if ( id_vts> 0 ) then
           allocate ( vt_s(is:ie, js:je, ks:ke) )
           vt_s = 0.
      endif
      if ( id_vtg> 0 ) then
           allocate ( vt_g(is:ie, js:je, ks:ke) )
           vt_g = 0.
      endif
      if ( id_vti> 0 ) then
           allocate ( vt_i(is:ie, js:je, ks:ke) )
           vt_i = 0.
      endif

      allocate (     cond(is:ie, js:je) )
      allocate (  prec_mp(is:ie, js:je) )
      allocate (    prec0(is:ie, js:je) )
      allocate (    prec1(is:ie, js:je) )
      allocate (    rain0(is:ie, js:je) )
      allocate (    snow0(is:ie, js:je) )
      allocate (     ice0(is:ie, js:je) )
      allocate ( graupel0(is:ie, js:je) )

      prec0 = 0.
      prec1 = 0.
      rain0 = 0.
      snow0 = 0.
       ice0 = 0.
   graupel0 = 0.
 

 end subroutine setup_con


 subroutine bergrn(tc,ql,qi,qlrho,psfw,psfi)
      real, intent(inout) :: psfi, psfw
      real, intent(in):: qlrho, qi, ql, tc
      real :: dt1, a21, a2,  a1, tc1, a1t, a2t
      dimension a1t(0:31),a2t(0:31)
      data a1t/0.,0.7939e-7,0.7841e-6,0.3369e-5,0.4336e-5,0.5285e-5, &
              0.3728e-5,0.1852e-5,0.2991e-6,0.4248e-6,0.7434e-6, &
              0.1812e-5,0.4394e-5,0.9145e-5,0.1725e-4,0.3348e-4, &
              0.1725e-4,0.9175e-5,0.4412e-5,0.2252e-5,0.9115e-6, &
              0.4876e-6,0.3473e-6,0.4758e-6,0.6306e-6,0.8573e-6, &
              0.7868e-6,0.7192e-6,0.6513e-6,0.5956e-6,0.5333e-6, &
              0.4834e-6/
      data a2t/0.,0.4006,0.4831,0.5320,0.5307,0.5319, &
              0.5249,0.4888,0.3894,0.4047,0.4318, &
              0.4771,0.5183,0.5463,0.5651,0.5813, &
              0.5655,0.5478,0.5203,0.4906,0.4447, &
              0.4126,0.3960,0.4149,0.4320,0.4506, &
              0.4483,0.4460,0.4433,0.4413,0.4382, &
              0.4361/
 
      tc1 = max(tc, -30.0)
      if ( tc1 .gt. -1.0 ) then
      a1 = a1t(1)
      a2 = a2t(1)
      else
      a1  = (a1t(-int(tc1))-a1t(-int(tc1)+1))*(tc1-int(tc1)+1.0)+ &
           a1t(-int(tc1)+1)
      a2  = (a2t(-int(tc1))-a2t(-int(tc1)+1))*(tc1-int(tc1)+1.0)+ &
           a2t(-int(tc1)+1)
      endif
      a21 = 1.0 - a2
      dt1 = (rmi50**a21 - rmi40**a21) / (a1*a21)
!     note:  mks units, ui50=1.0 m/sec, eiw=1.0
      psfw = c1brg*qi/dt1*(a1*rmi50**a2+qlrho*c2brg)
      psfi = qi/dt1        ! sjl

 end subroutine bergrn


 real function acr3d(v1, v2, q1, q2, c, cac, denfac, rho)
 real, intent(in) :: v1,v2,c, denfac,rho
 real, intent(in) :: q1, q2    ! mixing ratio!!!
 real, intent(in) :: cac(3)
 real :: a
 integer :: k
      a=0.0
      do k=1,3
         a = a + cac(k)*( (q1*rho)**((7-k)*0.25) * (q2*rho)**(k*0.25) )
      enddo
      acr3d = c * abs(v1-v2) * a / rho
 end function acr3d


 real function smlt(tc, dqs, qsrho,psacw,psacr,c,rho, rhofac)
 real, intent(in)::  tc,dqs,qsrho,psacw,psacr,c(5),rho, rhofac
     
 smlt = (c(1)*tc/rho-c(2)*dqs) * (c(3)*sqrt(qsrho)+ &
         c(4)*qsrho**0.65625*sqrt(rhofac)) + c(5)*tc*(psacw+psacr)

 end function smlt
 

 real function gmlt(tc, dqs,qgrho,pgacw,pgacr,c, rho)
 real, intent(in)::  tc,dqs,qgrho,pgacw,pgacr,c(5),rho
     
!     note:  pgacw and pgacr must be calc before gmlt is called
!
 gmlt = (c(1)*tc/rho-c(2)*dqs) * (c(3)*sqrt(qgrho)+ &
         c(4)*qgrho**0.6875/rho**0.25) + c(5)*tc*(pgacw+pgacr)
 end function gmlt


 subroutine qsmith_init
  integer, parameter:: length=2621 
  integer i

  if( .not. allocated(table) ) then
!                            generate es table (dt = 0.1 deg. c)
       allocate ( table( length) )
       allocate ( table2(length) )
       allocate ( tablew(length) )
       allocate (   des (length) )
       allocate (   des2(length) )
       allocate (   desw(length) )

       call qs_table (length )
       call qs_table2(length )
       call qs_tablew(length )

       do i=1,length-1
           des(i) = max(0.,  table(i+1) -  table(i))
          des2(i) = max(0., table2(i+1) - table2(i))
          desw(i) = max(0., tablew(i+1) - tablew(i))
       enddo
        des(length) =  des(length-1)
       des2(length) = des2(length-1)
       desw(length) = desw(length-1)
  endif
 
 end subroutine qsmith_init
 

 real function qs1d(ta, pa, dqdt)
! 2-phase tabel
  real, intent(in):: ta, pa
  real, intent(out):: dqdt
! local:
  real es, ap1
  real, parameter:: tmin=tice - 160.
  real, parameter:: eps10 = 10.*eps
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = table2(it) + (ap1-it)*des2(it)
      qs1d = eps*es/max(es, pa-es)
        it = ap1 - 0.5
      dqdt = eps10*(des2(it) + (ap1-it)*(des2(it+1)-des2(it)))/max(es, pa-es)

 end function qs1d


 real function ws1d(ta, pa, dqdt)
! Pure water phase
  real, intent(in):: ta, pa
  real, intent(out):: dqdt
! local:
  real es, ap1
  real, parameter:: tmin=tice - 160.
  real, parameter:: eps10 = 10.*eps
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = tablew(it) + (ap1-it)*desw(it)
      ws1d = eps*es/max(es, pa-es)
        it = ap1 - 0.5
      dqdt = eps10*(desw(it) + (ap1-it)*(desw(it+1)-desw(it)))/max(es, pa-es)

 end function ws1d


 real function es2_table(ta)
! two-phase table
  real, intent(in):: ta
  real, parameter:: tmin=tice - 160.
  real  ap1
  integer it
       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
      es2_table = table2(it) + (ap1-it)*des2(it)
 end function es2_table


 real function esw_table(ta)
! pure water phase table
  real, intent(in):: ta
  real, parameter:: tmin=tice - 160.
  real  ap1
  integer it
       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
      esw_table = tablew(it) + (ap1-it)*desw(it)
 end function esw_table


 subroutine qs_tablew(n)
! 2-phase table
      integer, intent(in):: n
      real:: delt=0.1
      real esbasw, tbasw, esbasi, tbasi, tmin, tem, aa, b, c, d, e
      integer i

! constants
      esbasw = 1013246.0
       tbasw =     373.16
      esbasi =    6107.1
       tbasi =     273.16
        tmin = tbasi - 160.

     do i=1,n
        tem = tmin+delt*real(i-1)
!  compute es over water
!  see smithsonian meteorological tables page 350.
        aa  = -7.90298*(tbasw/tem-1.)
        b   =  5.02808*alog10(tbasw/tem)
        c   = -1.3816e-07*(10**((1.-tem/tbasw)*11.344)-1.)
        d   =  8.1328e-03*(10**((tbasw/tem-1.)*(-3.49149))-1.)
        e   = alog10(esbasw)
        tablew(i) = 0.1 * 10**(aa+b+c+d+e)
     enddo

 end subroutine qs_tablew


 subroutine qs_table2(n)
! 2-phase table
  integer, intent(in):: n
  real:: delt=0.1
  real esbasw, tbasw, esbasi, tbasi, tmin, tem, aa, b, c, d, e
  integer :: i0, i1
  real :: tem0, tem1
  integer i

! constants
      esbasw = 1013246.0
       tbasw =     373.16
      esbasi =    6107.1
       tbasi =     273.16
      tmin = tbasi - 160.

     do i=1,n
        tem = tmin+delt*real(i-1)
        if ( i<= 1600 ) then
!  compute es over ice between -160c and 0 c.
!  see smithsonian meteorological tables page 350.
              aa  = -9.09718 *(tbasi/tem-1.)
              b   = -3.56654 *alog10(tbasi/tem)
              c   =  0.876793*(1.-tem/tbasi)
              e   = alog10(esbasi)
             table2(i) = 0.1 * 10**(aa+b+c+e)
        else
!  compute es over water between 0c and 102c.
!  see smithsonian meteorological tables page 350.
             aa  = -7.90298*(tbasw/tem-1.)
             b   =  5.02808*alog10(tbasw/tem)
             c   = -1.3816e-07*(10**((1.-tem/tbasw)*11.344)-1.)
             d   =  8.1328e-03*(10**((tbasw/tem-1.)*(-3.49149))-1.)
             e   = alog10(esbasw)
             table2(i) = 0.1 * 10**(aa+b+c+d+e)
        endif
     enddo

!----------
! smoother
!----------
      i0 = 1600;  i1 = 1601
      tem0 = 0.25*(table2(i0-1) + 2.*table(i0) + table2(i0+1))
      tem1 = 0.25*(table2(i1-1) + 2.*table(i1) + table2(i1+1))
      table2(i0) = tem0
      table2(i1) = tem1

 end subroutine qs_table2



 real function qs1d_blend(t, p, dqdt)
! this routine is for use in mp_lin; which uses dry mixing ratio
! Blended
  real, intent(in):: t, p
  real, intent(out):: dqdt
! local:
  real es, ap1
  real, parameter:: tmin=tice - 160.
  integer it

       ap1 = 10.*dim(t, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = table(it) + (ap1-it)*des(it)
      qs1d_blend = eps*es/max(es, p-es)
        it = ap1 - 0.5
      dqdt = 10.*eps*(des(it) + (ap1-it)*(des(it+1)-des(it)))/max(es, p-es)

 end function qs1d_blend

 subroutine qs_table(n)
      integer, intent(in):: n
      real esupc(200)
      real:: delt=0.1
      real esbasw, tbasw, esbasi, tbasi, tmin, tem, aa, b, c, d, e, esh20 
      real wice, wh2o
      integer i

! constants
      esbasw = 1013246.0
       tbasw =     373.16
      esbasi =    6107.1
       tbasi =     273.16

!  compute es over ice between -160c and 0 c.
      tmin = tbasi - 160.
!  see smithsonian meteorological tables page 350.
      do i=1,1600
         tem = tmin+delt*real(i-1)
         aa  = -9.09718 *(tbasi/tem-1.)
         b   = -3.56654 *alog10(tbasi/tem)
         c   =  0.876793*(1.-tem/tbasi)
         e   = alog10(esbasi)
         table(i)=10**(aa+b+c+e)
      enddo

!  compute es over water between -20c and 102c.
!  see smithsonian meteorological tables page 350.
      do  i=1,1221
          tem = 253.16+delt*real(i-1)
          aa  = -7.90298*(tbasw/tem-1.)
          b   =  5.02808*alog10(tbasw/tem)
          c   = -1.3816e-07*(10**((1.-tem/tbasw)*11.344)-1.)
          d   =  8.1328e-03*(10**((tbasw/tem-1.)*(-3.49149))-1.)
          e   = alog10(esbasw)
          esh20  = 10**(aa+b+c+d+e)
          if (i <= 200) then
              esupc(i) = esh20
          else
              table(i+1400) = esh20
          endif
      enddo

!  derive blended es over ice and supercooled water between -20c and 0c
      do i=1,200
         tem  = 253.16+delt*real(i-1)
         wice = 0.05*(273.16-tem)
         wh2o = 0.05*(tem-253.16)
         table(i+1400) = wice*table(i+1400)+wh2o*esupc(i)
      enddo

      do i=1,n
         table(i) = table(i)*0.1
      enddo

 end subroutine qs_table


  subroutine qsmith(im, km, k1, t, p, q, dqdt)
! input t in deg k; p (pa)
  integer, intent(in):: im, km, k1
  real, intent(in),dimension(im,km):: t, p
  real, intent(out),dimension(im,km):: q
  real, intent(out), optional:: dqdt(im,km)
! local:
  real es(im,km)
  real ap1
  real tmin, oms
  integer i, k, it

  tmin = tice-160.

  oms = 1. - eps
 
  if( .not. allocated(table) ) then
       call  qsmith_init
  endif
 
      do k=k1,km
         do i=1,im
            ap1 = 10.*dim(t(i,k), tmin) + 1.
            ap1 = min(2621., ap1)
            it = ap1
            es(i,k) = table(it) + (ap1-it)*des(it)
            q(i,k) = eps*es(i,k) / max(es(i,k), p(i,k)-oms*es(i,k))
            q(i,k) = min(1., q(i,k)) 
         enddo
      enddo

      if ( present(dqdt) ) then
      do k=k1,km
           do i=1,im
              ap1 = 10.*dim(t(i,k), tmin) + 1.
              ap1 = min(2621., ap1) - 0.5
              it  = ap1
              dqdt(i,k) = 10.*eps*(des(it) + (ap1-it)*(des(it+1)-des(it)))/p(i,k)
           enddo
      enddo
      endif
 
 end subroutine qsmith



 subroutine neg_adj(ktop, kbot, p1, pt, dp, qv, ql, qr, qi, qs, qg)
! 1d version:
! this is designed for 6-class micro-physics schemes
 integer, intent(in):: ktop, kbot
 real, intent(in):: dp(ktop:kbot), p1(ktop:kbot)
 real, intent(inout), dimension(ktop:kbot)::    &
                                pt, qv, ql, qr, qi, qs, qg
! local:
 real lcpk(ktop:kbot), icpk(ktop:kbot)
 real dq, tmp1
 integer k

 do k=ktop,kbot
       tmp1 = cp - rdgas*ptop/p1(k)
    lcpk(k) = latv / tmp1
    icpk(k) = lati / tmp1
 enddo

 do k=ktop, kbot
!-----------
! ice-phase:
!-----------
! if ice<0 borrow from snow
          if( qi(k) < 0. ) then
              qs(k) = qs(k) + qi(k)
              qi(k) = 0.
          endif
#ifdef NO_GRAUPEL
    qs(k) = qs(k) + qg(k)
    qg(k) = 0.
! if snow < 0 then borrow from rain
          if ( qs(k) < 0. ) then
               qr(k) = qr(k) + qs(k)
               pt(k) = pt(k) - qs(k)*icpk(k)   ! heating
               qs(k) = 0.
          endif
#else
! if snow<0 borrow from graupel
          if( qs(k) < 0. ) then
              qg(k) = qg(k) + qs(k)
              qs(k) = 0.
          endif
! if graupel < 0 then borrow from rain
          if ( qg(k) < 0. ) then
               qr(k) = qr(k) + qg(k)
               pt(k) = pt(k) - qg(k)*icpk(k)   ! heating
               qg(k) = 0.
          endif
#endif

! liquid phase:
! fix negative rain by borrowing from cloud water
          if ( qr(k) < 0. ) then
               ql(k) = ql(k) + qr(k)
               qr(k) = 0.
          endif
! fix negative cloud water with vapor
          if ( ql(k) < 0. ) then
               qv(k) = qv(k) + ql(k)
               pt(k) = pt(k) - ql(k)*lcpk(k)
               ql(k) = 0.
          endif
 enddo

!-----------------------------------
! fix water vapor; borrow from below
!-----------------------------------
 do k=ktop,kbot-1
    if( qv(k) < 0. ) then
        qv(k+1) = qv(k+1) + qv(k)*dp(k)/dp(k+1)
        qv(k  ) = 0.
    endif
 enddo
 
! bottom layer; borrow from above
 if( qv(kbot) < 0. .and. qv(kbot-1)>0.) then
             dq = min(-qv(kbot)*dp(kbot), qv(kbot-1)*dp(kbot-1))
     qv(kbot-1) = qv(kbot-1) - dq/dp(kbot-1) 
     qv(kbot  ) = qv(kbot  ) + dq/dp(kbot  ) 
 endif
! if qv is still < 0

 end subroutine neg_adj


 subroutine neg_adj3(is, ie, js, je, ng, kbot,      &
                     pt, dp, qv, ql, qr, qi, qs, qg)

! this is designed for 6-class micro-physics schemes
 integer, intent(in):: is, ie, js, je, ng, kbot
 real, intent(in):: dp(is-ng:ie+ng,js-ng:je+ng,kbot)
 real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,kbot)::    &
                                 pt, qv, ql, qr, qi, qs, qg
! local:
 real dq
 integer i, j, k

 do k=1, kbot
    do j=js, je
       do i=is, ie
!-----------
! ice-phase:
!-----------
! if ice<0 borrow from snow
          if( qi(i,j,k) < 0. ) then
              qs(i,j,k) = qs(i,j,k) + qi(i,j,k)
              qi(i,j,k) = 0.
          endif
! if snow<0 borrow from graupel
          if( qs(i,j,k) < 0. ) then
              qg(i,j,k) = qg(i,j,k) + qs(i,j,k)
              qs(i,j,k) = 0.
          endif
! if graupel < 0 then borrow from rain
          if ( qg(i,j,k) < 0. ) then
               qr(i,j,k) = qr(i,j,k) + qg(i,j,k)
               pt(i,j,k) = pt(i,j,k) - qg(i,j,k)*icp
               qg(i,j,k) = 0.
          endif

! liquid phase:
! fix negative rain by borrowing from cloud water
          if ( qr(i,j,k) < 0. ) then
               ql(i,j,k) = ql(i,j,k) + qr(i,j,k)
               qr(i,j,k) = 0.
          endif
! fix negative cloud water with vapor
          if ( ql(i,j,k) < 0. ) then
               qv(i,j,k) = qv(i,j,k) + ql(i,j,k)
               pt(i,j,k) = pt(i,j,k) - ql(i,j,k)*lcp
               ql(i,j,k) = 0.
          endif
     enddo
   enddo
 enddo

!-----------------------------------
! fix water vapor; borrow from below
!-----------------------------------
 do k=1,kbot-1
    do j=js, je
       do i=is, ie
          if( qv(i,j,k) < 0. ) then
              qv(i,j,k+1) = qv(i,j,k+1) + qv(i,j,k)*dp(i,j,k)/dp(i,j,k+1)
              qv(i,j,k  ) = 0.
          endif
     enddo
   enddo
 enddo
 
! bottom layer; borrow from above
  do j=js, je
     do i=is, ie
        if( qv(i,j,kbot) < 0. .and. qv(i,j,kbot-1)>0.) then
            dq = min(-qv(i,j,kbot)*dp(i,j,kbot), qv(i,j,kbot-1)*dp(i,j,kbot-1))
            qv(i,j,kbot-1) = qv(i,j,kbot-1) - dq/dp(i,j,kbot-1) 
            qv(i,j,kbot  ) = qv(i,j,kbot  ) + dq/dp(i,j,kbot  ) 
        endif
! if qv is still < 0
   enddo
 enddo

 end subroutine neg_adj3


end module mp_lin_mod
