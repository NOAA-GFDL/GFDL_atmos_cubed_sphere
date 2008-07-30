module mp_lin_mod
 use fv_diagnostics_mod,only: prt_maxmin
 use fv_grid_tools_mod, only: area
 use fv_grid_utils_mod, only: g_sum, ptop
 use fv_timing_mod,     only: timing_on, timing_off
 use fv_mp_mod,         only: gid, ng, masterproc
 use mpp_mod,           only: stdlog
 use diag_manager_mod,  only: register_diag_field, send_data
 use time_manager_mod,  only: time_type, get_date, get_time
 use constants_mod,     only: grav, rdgas, rvgas, cp_air, hlv, hlf, kappa
 use fms_mod,           only: write_version_number, open_namelist_file, &
                              check_nml_error, file_exist, close_file,  &
                              error_mesg, FATAL 
 use fv_control_mod,    only: hydrostatic, phys_hydrostatic

 implicit none
 private

 character(len=128) :: version = ''
 character(len=128) :: tagname = ''

 public  micro_phys_driver, micro_phys_init, micro_phys_end, sg_conv
 public  qsmith_init, qsmith, es2_table1d, es3_table1d, esw_table1d
 real             :: missing_value = -1.e10
 logical          :: module_is_initialized = .false.
 character(len=7) :: mod_name = 'mp_lin'

!==== fms constants ====================
 real, parameter :: cp    = cp_air       ! heat capacity at constant pressure (j/kg/k)
 real, parameter :: eps   = rdgas/rvgas  ! = 0.621971831
 real, parameter :: zvir  = rvgas/rdgas-1.    ! = 0.607789855
 real, parameter :: latv  = hlv          ! = 2.500e6
 real, parameter :: lati  = hlf          ! = 3.34e5
 real, parameter :: lats  = hlv+hlf
!==== fms constants ====================
 real, parameter :: qvmin  = 1.e-22
 real, parameter :: qrmin  = 1.e-9
 real, parameter :: sfcrho = 1.20         ! surface air density
 real, parameter :: vmin   = 1.e-2        ! minimum fall speed for rain/graupel
 real, parameter :: tice   = 273.16
 real, parameter :: tice2  = 275.16
 real, parameter :: rhor   = 1.0e3  ! LFO83

 real :: cracs,csacr,cgacr,cgacs,acco(3,4),csacw,               &
         craci, csaci,cgacw,cgaci,cracw,cssub(5),cgsub(5), &
         crevp(5), cgfr(2), csmlt(5), cgmlt(5)
 real :: rmi50, es0, ces0, c1brg, c2brg

 real, parameter:: dz_min = 1.e-2
! Derived variables:
 real :: dts, rdts, pie, f_exp
 real :: lcp, icp, tcp, rgrav
 real :: c_praut
 real :: fac_rc
 real :: mp_count = 0.

 logical :: do_setup=.true.
 logical :: master 

 real, allocatable:: vt_r(:,:,:), vt_s(:,:,:), vt_g(:,:,:), vt_i(:,:,:)
 real, allocatable:: prec0(:,:), rain0(:,:), snow0(:,:), ice0(:,:), graupel0(:,:)
 real, allocatable:: prec1(:,:), prec_mp(:,:), cond(:,:) 
 real, allocatable:: table(:), table2(:), table3(:), tablew(:), des(:), des2(:), des3(:), desw(:)

 integer:: isc, iec, jsc, jec
 integer:: id_vtr, id_vts,  id_vtg, id_vti, id_rain, id_snow, id_graupel, &
           id_ice, id_prec, id_cond

 real :: qc_crt = 1.0e-6

!----------------------
! namelist  parameters:
!----------------------
 real :: mp_time = 150.   ! maximum micro-physics time step (sec)
 real ::  tau_i  =  10.   ! (sec) cloud ice melt
 real ::  tau_s  = 120.   ! snow melt
 real ::  tau_r  = 150.   ! freezing rain time scale; not used !!!
 real ::  tau_g  = 180.   ! graupel melt
 real ::  tau_evap = 900.  ! not used!!!
 real ::  tau_wsub = 150. ! cloud water erosion/evaporation time scale
 real ::  rh_l   = 0.80   ! RH threshold over land
 real ::  rh_o   = 0.80   ! RH threshold over ocean
 real ::  cfac   = 20.0   ! 
 real :: evap_rhc = 1.0   ! rain re-evaporation threshold

! NCAR CAM settings (near surafce): 400 (land), 150 (ocean), 75 (sea ice) 
 real :: ccn   = 150.     !
 real :: ccn_o = 100.    
 real :: ccn_l = 300.    
 real :: rthresh = 8.0e-6     ! critical cloud drop radius (micro m)

!-------------------------------------------------------------
 real :: qi0_crt = 1.2e-4   ! ice  -> snow autocon threshold (density)
!real :: qi0_crt = 6.0e-4   ! ice  -> snow autocon threshold (mixing ratio)

 real :: qr0_crt = 1.0e-4   ! rain --> snow or graupel threshold
 real :: c_psaut = 1.0e-3   ! autoconversion rate: cloud_ice -> snow
 real :: c_psaci = 0.1      ! accretion: cloud ice --> snow (was 0.1 in Zetac)
 real :: c_piacr = 250.     ! accretion: rain --> ice (was 2.e2)
 real :: c_cracw = 1.0      ! rain accretion factor

!-----------------
! Graupel control:
!-----------------
 real :: qs0_crt = 2.5e-3   ! snow --> graupel threshold (was 1.e-3)
 real :: c_pgacs = 3.0e-3   ! snow --> graupel "accretion" eff. (was 0.1 in Zetac)

! fall velocity tuning constants:
 real :: den_ref = 1.2      ! Reference (surface) density for fall speed
                            ! Larger value produce larger fall speed
 real :: vr_fac = 1.
 real :: vs_fac = 1.
 real :: vg_fac = 1.
 real :: vi_fac = 1.
 real :: thi    = 1.0e-9    ! cloud ice *density* threshold for terminal fall
 real :: i2r_crt = 2.5e-3   ! ice --> rain/cloud water melting threshold
 real :: p_crt   = 200.E2   ! 

 logical :: use_ccn  = .false.
 logical :: use_ppm  = .true.
 logical :: falling_rain_mp = .true.
 logical :: use_cld_prof = .false.
 logical :: mp_debug = .false.
 logical :: mp_print = .true.

 integer :: k_moist = 99

 namelist /mp_lin_nml/mp_time, tau_r, tau_evap, tau_s, tau_i, tau_g,  &
                      vr_fac, vs_fac, vg_fac, vi_fac, use_cld_prof,  &
                      qs0_crt, qi0_crt, qr0_crt,   &
                      rh_l, rh_o, cfac, den_ref, use_ccn, &
                      rthresh, ccn_l, ccn_o,  evap_rhc,     &
                      thi, c_piacr, tau_wsub,       &
                      c_psaut, c_psaci, c_pgacs, falling_rain_mp, &
                      c_cracw, i2r_crt, k_moist, p_crt,     &
                      use_ppm, mp_debug, mp_print

 contains
 

  subroutine micro_phys_driver(qv,    ql,    qr,    qi,    qs,    qg,    qa,    &
                               qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, qa_dt, & 
                               pt_dt, pt, p3, dz,  delp, dt_in,                 &
                               land,  rain, snow, ice, graupel,                 &
                               is,ie, js,je, ks,ke, ktop, kbot, time)

  type(time_type), intent(in):: time
  integer,         intent(in):: is,ie, js,je  ! physics window
  integer,         intent(in):: ks,ke         ! vertical dimension
  integer,         intent(in):: ktop, kbot    ! vertical compute domain
  real,            intent(in):: dt_in

  real, intent(in   ), dimension(is:ie,js:je)      :: land  !land fraction
  real, intent(out  ), dimension(is:ie,js:je)      :: rain, snow, ice, graupel
  real, intent(in   ), dimension(is:ie,js:je,ks:ke):: p3, delp    ! p3 not used
  real, intent(in   ), dimension(is:ie,js:je,ks:ke):: pt, qv, ql, qr, qi, qs, qg, qa
  real, intent(inout), dimension(is:ie,js:je,ks:ke):: pt_dt,  qa_dt, dz
  real, intent(inout), dimension(is:ie,js:je,ks:ke):: qv_dt, ql_dt, qr_dt, qi_dt,  &
                                                      qs_dt, qg_dt

! local:
  logical used
  real    :: mpdt, rdt, convt, tot_prec
  integer :: i,j,k
  integer :: seconds, days, ntimes

                                        call timing_on (" split_mp")

! tendency zero out for am moist processes should be done outside the driver

     mpdt = min(dt_in, mp_time)
      rdt = 1. / dt_in
   ntimes = nint( dt_in/mpdt )
! small time step:
      dts = dt_in / real(ntimes)
     rdts = 1./dts
    f_exp = 1. - exp(-0.5*dts/tau_wsub)

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
     call mpdrv( delp, pt, qv, ql, qr, qi, qs, qg, qa, dz,  &
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

#ifdef SIM_PHYS
   if ( id_vtr> 0 ) used=send_data(id_vtr, vt_r, time)
   if ( id_vts> 0 ) used=send_data(id_vts, vt_s, time)
   if ( id_vtg> 0 ) used=send_data(id_vtg, vt_g, time)
   if ( id_vts> 0 ) used=send_data(id_vti, vt_i, time)
#else
   if ( id_vtr> 0 ) used=send_data(id_vtr, vt_r, time, is, js)
   if ( id_vts> 0 ) used=send_data(id_vts, vt_s, time, is, js)
   if ( id_vtg> 0 ) used=send_data(id_vtg, vt_g, time, is, js)
   if ( id_vts> 0 ) used=send_data(id_vti, vt_i, time, is, js)
#endif

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
#ifdef SIM_PHYS
        used=send_data(id_cond, cond, time)
#else
        used=send_data(id_cond, cond, time, is, js)
#endif
   endif

   if ( id_snow>0 ) then
#ifdef SIM_PHYS
        used=send_data(id_snow,    snow,    time)
#else
        used=send_data(id_snow,    snow,    time, is, js)
#endif
        if ( seconds==0 ) then
             tot_prec = g_sum(snow, is, ie, js, je, ng, area, 1) 
             if(master) write(*,*) 'mean snow=', tot_prec
        endif
        snow0(:,:) = snow0(:,:) + snow(:,:)
   endif

   if ( id_graupel>0 ) then
#ifdef SIM_PHYS
        used=send_data(id_graupel, graupel, time)
#else
        used=send_data(id_graupel, graupel, time, is, js)
#endif
        if ( seconds==0 ) then
             tot_prec = g_sum(graupel, is, ie, js, je, ng, area, 1) 
             if(master) write(*,*) 'mean graupel=', tot_prec
        endif
        graupel0(:,:) = graupel0(:,:) + graupel(:,:)
   endif

   if ( id_ice>0 ) then
#ifdef SIM_PHYS
        used=send_data(id_ice, ice, time)
#else
        used=send_data(id_ice, ice, time, is, js)
#endif
        if ( seconds==0 ) then
             tot_prec = g_sum(ice, is, ie, js, je, ng, area, 1) 
             if(master) write(*,*) 'mean ice_mp=', tot_prec
        endif
        ice0(:,:) = ice0(:,:) + ice(:,:)
   endif

   if ( id_rain>0 ) then
#ifdef SIM_PHYS
        used=send_data(id_rain,    rain,    time)
#else
        used=send_data(id_rain,    rain,    time, is, js)
#endif
        if ( seconds==0 ) then
!            tot_prec = g_sum(rain, is, ie, js, je, ng, area, 1) 
!            if(master) write(*,*) 'mean rain=', tot_prec
        endif
        rain0(:,:) = rain0(:,:) + rain(:,:)
   endif
   

   if ( id_prec>0 ) then
#ifdef SIM_PHYS
        used=send_data(id_prec, prec_mp, time)
#else
        used=send_data(id_prec, prec_mp, time,is, js)
#endif
   endif

!----------------------------------------------------------------------------

        prec0(:,:) = prec0(:,:) + prec_mp(:,:)
        prec1(:,:) = prec1(:,:) + prec_mp(:,:)
        mp_count = mp_count + 1.

        if ( seconds==0 .and. mp_print ) then
             tot_prec = g_sum(prec1*dt_in/86400., is, ie, js, je, ng, area, 1) 
             if(master) write(*,*) 'Daily prec_mp=', tot_prec
!            call prt_maxmin('prec_mp', prec1*dt_in/86400., is, ie, js, je, 0, 1, 1., master)
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



 subroutine mpdrv( delp, pt, qv, ql, qr, qi, qs, qg, qa, dz,     &
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

  real, intent(in), dimension(is:ie,js:je,ks:ke) :: delp
  real, intent(in), dimension(is:ie):: land  !land fraction
  real, intent(in   ), dimension(is:ie,js:je,ks:ke):: pt, qv, ql, qr, qi, qs, qg, qa
  real, intent(inout), dimension(is:ie,js:je,ks:ke):: pt_dt,  qa_dt, dz,   &
                                            qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt
  real, intent(out), dimension(is:ie):: rain, snow, ice, graupel, cond
!----------
! local var
!----------
  real, dimension(ktop:kbot):: qvz, qlz, qrz, qiz, qsz, qgz, qaz, &
                               vtiz, vtsz, vtgz, vtrz, &
                               dp1, qv0, ql0, qr0, qi0, qs0, qg0, qa0, t0, den, &
                               den0, tz, p1, dz0, dz1, denfac, rh_p, rh_c, rh

  real :: r1, s1, i1, g1, rdt, omq, qaz1
  real :: qpz, clouds, cpaut
  real :: t_fac, p_fac, rh_t1, rh_t2
  real :: dt_rain
  integer :: i,k,n
! real:: x, pexp
! pexp(x) = 1.+x*(1.+x*(0.5+x/6.*(1.+x*(0.25+0.05*x))))

   dt_rain = dts * 0.5

   rdt = 1. / dt_in

   cpaut = 0.55*0.104*grav/1.717e-5

   do 2000 i=is, ie

   do k=ktop, kbot
       t0(k) = pt(i,j,k)
       tz(k) = t0(k) 
!-----------------------------------
      qvz(k) = max(qvmin, qv(i,j,k))
!-----------------------------------
      qlz(k) = ql(i,j,k)
      qrz(k) = qr(i,j,k)
      qiz(k) = qi(i,j,k)
      qsz(k) = qs(i,j,k)
      qgz(k) = qg(i,j,k)
      qa0(k) = qa(i,j,k)
      qaz(k) = 0.  
      dz0(k) = dz(i,j,k)
!--------------------------
         omq = 1. - (qvz(k)+qlz(k)+qrz(k)+qiz(k)+qsz(k)+qgz(k))
      dp1(k) = delp(i,j,k) * omq         ! dry air mass * grav
     den0(k) = -dp1(k)/(grav*dz0(k))     ! density of dry air
       p1(k) = den0(k)*rdgas*t0(k)       ! dry pressure
!------------------------------
! convert to dry mixing ratios:
!------------------------------
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
   enddo

! Compute dry pressure for non-hydrostatic case
!-----------------------------------------------
!  if ( .not. phys_hydrostatic ) then
!      do k=ktop, kbot
!         p1(k) = den0(k)*rdgas*t0(k)
!      enddo
!  endif
!-----------------------------------------------

! Based on Klein Eq. 15
   ccn = (ccn_l*land(i) + ccn_o*(1.-land(i))) * 1.e6
   if ( use_ccn ) then
!      CCN = ccn_surf * (den/den_surf)
       ccn = ccn * rdgas*tz(kbot)/p1(kbot)
   endif
   c_praut =  cpaut * (ccn*rhor)**(-1./3.)

   t_fac = max(0.,min(1.,(tz(kbot)-tice)/(30.)))

   !For land:
   rh_t1 = 1. - (1.-rh_l)*t_fac
   !For ocean
   rh_t2 = 1. - (1.-rh_o)*t_fac

   do k=ktop, kbot
      if ( p1(k) > 950.e2 ) then 
           rh_c(k) = 1.
      elseif ( p1(k) > 800.e2 ) then 
           p_fac= (p1(k) - 800.e2) / 150.e2
           p_fac= p_fac ** 2
           rh_c(k) = (rh_t1 + (1.-rh_t1)*p_fac)*land(i) +  &
                     (rh_t2 + (1.-rh_t2)*p_fac)*(1.-land(i))
      else
           rh_c(k) = rh_t1*land(i) + rh_t2*(1.-land(i)) 
      endif
   enddo

   if ( use_cld_prof ) then
      do k=ktop, kbot
         rh_p(k) = rh_c(k)
      enddo
   else
      do k=ktop, kbot
         rh_p(k) = 1.
      enddo
   endif

!-------------------------
! * fix all negatives
!-------------------------

 call neg_adj(ktop, kbot, p1, tz, dp1, qvz, qlz, qrz, qiz, qsz, qgz)

 do 1000 n=1,ntimes

   do k=ktop, kbot
         dz1(k) = dz0(k)*tz(k)/t0(k) 
         den(k) = den0(k)*dz0(k)/dz1(k)
      denfac(k) = sqrt(sfcrho/den(k))
   enddo

!-------------------------------------------
! Time-split warm rain processes: first pass
!-------------------------------------------
!                                       call timing_on (" warm_rain")
   call warm_rain(dt_rain, ktop, kbot, dp1, dz1, tz, qvz, qlz, qrz, p1, den, denfac, vtrz, r1)
!                                       call timing_off(" warm_rain")
   rain(i) = rain(i) + r1

!------------------------------------------------
! * sedimentation of cloud ice, snow, and graupel
!------------------------------------------------
!                                       call timing_on (" terminal_fall")
   call fall_speed(ktop, kbot, den, qsz, qiz, qgz, vtsz, vtiz, vtgz)

   call terminal_fall ( dts, ktop, kbot, tz, qvz, qlz, qrz, qgz, qsz, qiz, p1, &
                        dz1, dp1, den, denfac, vtgz, vtsz, vtiz, rh_p,   &
                        r1, g1, s1, i1 )
!                                       call timing_off(" terminal_fall")

      rain(i) = rain(i)    + r1  ! from melted snow & ice that reached the ground
      snow(i) = snow(i)    + s1
   graupel(i) = graupel(i) + g1
       ice(i) = ice(i)     + i1

!-------------------------------------------
! Time-split warm rain processes: 2nd pass
!-------------------------------------------
!                                       call timing_on (" warm_rain")
   call warm_rain(dt_rain, ktop, kbot, dp1, dz1, tz, qvz, qlz, qrz, p1, den, denfac, vtrz, r1)
!                                       call timing_off(" warm_rain")
   rain(i) = rain(i) + r1

!-------------------------
! * ice-phase microphysics
!-------------------------

!                                       call timing_on (" icloud")
   call icloud( ktop, kbot, tz, p1, qvz, qlz, qrz, qiz, qsz, qgz,  &
                den, denfac, vtsz, vtgz, vtrz, rh_p, rh )
!                                       call timing_off(" icloud")

!-------------------------
! * cloud fraction
!-------------------------

   do k=ktop, kbot
      if ( rh(k) >= 1. ) then
! Full (100%) cloud cover
          qaz(k) = qaz(k) + 1.
      elseif ( rh(k)>rh_c(k) .and. (qlz(k)+qiz(k)+qsz(k)) > qc_crt ) then
! Partial cloudiness 
          qaz(k) = qaz(k) +  exp( -(cfac*rh_c(k)*(1.-rh(k)))**2 )
      endif
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
!        dz(i,j,k) = dz1(k)
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



 subroutine warm_rain( dt, ktop, kbot, dp, dz, tz, qv, ql, qr, pm,  &
                       den, denfac, vtr, r1)

 integer, intent(in):: ktop, kbot
 real,    intent(in):: dt                    ! time step (s)
 real,    intent(in),    dimension(ktop:kbot):: dp, dz, pm, den, denfac
 real, intent(inout), dimension(ktop:kbot):: tz, qv, ql, qr, vtr
 real, intent(out):: r1
! local:
 real, parameter:: fo3 = 4./3.
 real, parameter:: so3 = 7./3.
 real:: lcpk(ktop:kbot)
 real, dimension(ktop:kbot+1):: ze, zt
 real:: qsat, dqsdt, evap, factor, tsq
 real:: sink, dq
 real:: rho0, qden
 real:: zs = 0.
!-----------------------------------------------------------------------
! fall velocity constants:
!-----------------------------------------------------------------------
 real, parameter :: vconr = 2503.23638966667
 real, parameter :: normr = 25132741228.7183
 real, parameter :: thr=1.e-9
 real:: dt5
 integer k
 logical no_fall

! Terminal speed of falling rain:

  call check_column(ktop, kbot, qr, no_fall)
  if ( no_fall ) then
       vtr = vmin
       r1 = 0.
       go to 999
  endif

  if ( den_ref < 0. ) then
       rho0 = -den_ref*den(kbot) 
  else
       rho0 = den_ref   ! default=1.2
  endif

  do k=ktop, kbot
     qden = qr(k)*den(k)
     vtr(k) = max(vmin, vconr*sqrt(min(100., rho0/den(k)))      &
                             *exp(0.2*log(dim(qden,thr)/normr)))
  enddo

  dt5 = 0.5*dt
  ze(kbot+1) = zs
  do k=kbot, ktop, -1
     ze(k) = ze(k+1) - dz(k)  ! dz<0
  enddo
  zt(ktop) = ze(ktop)


 do k=ktop+1,kbot
    zt(k) = ze(k) - dt5*(vtr(k-1)+vtr(k))
 enddo
 zt(kbot+1) = zs - dt5*vtr(kbot)

 do k=ktop,kbot
    if( zt(k+1)>=zt(k) ) zt(k+1) = zt(k) - dz_min
 enddo

! call lagrangian_fall_ppm(ktop, kbot, zs, ze, zt, dp, qr, r1)
  call lagrangian_fall_pcm(ktop, kbot, zs, ze, zt, dp, qr, r1)

999  continue

!---------------------
! warm-rain processes:
!---------------------

  do k=ktop,kbot
!    lcpk(k) = latv / (cp - rdgas*ptop/pm(k))  ! Total energy conserving 
     lcpk(k) = lcp
  enddo

  do k=ktop,kbot

!-------------------
! * Rain evaporation
!-------------------
     if ( qr(k) > qrmin ) then
         qsat = ws1d(tz(k), pm(k), dqsdt)   ! water-phase table
         if ( qv(k) < qsat ) then
               tsq = tz(k)**2
              qden = qr(k)*den(k)
              evap = crevp(1)*tsq*(qsat-qv(k)) *     &
                    (crevp(2)*sqrt(qden)+crevp(3)*exp(0.725*log(qden)))   &
                   /(crevp(4)*tsq + crevp(5)*qsat*den(k))
              evap = min(qr(k), dt*evap, (qsat-qv(k))/(1.+lcpk(k)*dqsdt))
             qv(k) = qv(k) + evap
             qr(k) = qr(k) - evap
             tz(k) = tz(k) - evap*lcpk(k)
         endif
     endif             ! rain existed

     if ( ql(k) > qrmin ) then
!-------------------
! * Accretion: pracc
!-------------------
        if ( qr(k) > qrmin ) then
               qden = qr(k)*den(k)
             factor = dt*denfac(k)*cracw*exp(0.95*log(qden))
             sink = factor/(1.+factor)*ql(k)
        else
             sink = 0.
        endif

!-------------------
! * autoconversion
!-------------------
#ifdef USE_CCN
        if ( use_ccn ) then
! As in Klein's GFDL AM2 stratiform scheme.
! But CCN is formulted as CCN = CCN_surface * (den/den_surface)
             dq = ql(k) - fac_rc*ccn
             if ( dq > 0. ) then
                  sink = min( ql(k), sink + min(dq,     &
                              dt*c_praut*den(k)*ql(k)**so3 ))
             endif
        else
#endif
! Klein's GFDL AM2 stratiform scheme.
             dq = ql(k) - fac_rc*ccn/den(k)
             if ( dq > 0. ) then
! Note: ql is factored out of the (7/3) power computation
                  sink = min( ql(k), sink + min(dq,     &
                              dt*c_praut*ql(k)*exp(fo3*log(den(k)*ql(k)))) )
             endif
#ifdef USE_CCN
        endif
#endif
        ql(k) = ql(k) - sink
        qr(k) = qr(k) + sink
     endif    ! cloud water existed

  enddo


 end subroutine warm_rain





 subroutine icloud(ktop, kbot, tzk, p1, qvk, qlk, qrk, qik, qsk, qgk, den,    &
                   denfac, vts, vtg, vtr, rh_p, rh)

!----------------------------------------------------
! Bulk cloud micro-physics; processes splitting
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
!real, parameter:: tice1 = 274.16
 real, parameter:: t40 = tice - 40.
 real, parameter:: rhos = 0.1e3    ! snow density (1/10 of water)

 real, dimension(ktop:kbot) :: lcpk, icpk, tcpk
 real :: tz, qv, ql, qr, qi, qs, qg
 real :: praut, pracw, pracs, psacw, pgacw, pgmlt,   &
         psmlt, prevp, psacr, pgacr, pgfr,  pgacs,   &
         pgaut, pgaci, praci, psaut, psaci, pssub,   &
         pgsub, psfw,  psfi,  pidw,  piacr
 real :: tc, tsq, dqs0, qden, qi_tho
 real :: factor, frez, melt, sink
 real :: clouds, tmp1, rqi, tin, qsw, qsi, qpz
 real :: dqsdt, dwsdt, dq, qimin, pidep, ai
 real :: n0s, lamda
 integer :: i, j, k
 
! Taylor series expansion of A**B
!------------------------------------------------------------------------
!  A**B = 1 +  x + (1/2!)*x**2 + (1/3!)*x**3 + (1/4!)*x**4 + (1/5!)*x**5
!  where x = B*log(A)   this is useful if x << 1
!------------------------------------------------------------------------
! real:: x, pexp
! pexp(x) = 1.+x*(1.+x*(0.5+x/6.*(1.+x*(0.25+0.05*x))))

 real:: c3racs(3), c3sacr(3), c3gacr(3), c3gacs(3)
 real:: acco1(3,4)
 equivalence (c3racs(1),acco1(1,1)),(c3sacr(1),acco1(1,2)),  &
             (c3gacr(1),acco1(1,3)),(c3gacs(1),acco1(1,4))

 do j=1,4
    do i=1,3
       acco1(i,j) = acco(i,j)
    enddo
 enddo

 do k=ktop,kbot
!--------------------------------------
!      tmp1 = cp - rdgas*ptop/p1(k)
!   lcpk(k) =  latv / tmp1
!   icpk(k) =  lati / tmp1
!   tcpk(k) = lcpk(k) + icpk(k)
!--------------------------------------
    lcpk(k) = lcp
    icpk(k) = icp
    tcpk(k) = tcp
 enddo

 do k=ktop, kbot

! * pihom: homogeneous freezing of cloud water
    if( qlk(k) > qrmin ) then
        if( tzk(k) < t40 ) then
              frez = min(qlk(k), (t40-tzk(k))/icpk(k))
            qlk(k) = qlk(k) - frez
            qik(k) = qik(k) + frez
            tzk(k) = tzk(k) + frez*icpk(k)
        endif

! Biggs 1953: The supercooling water. Proc. Phys. Soc. London, B66, 688-694.
        if( qlk(k)>qrmin .and. tzk(k) < tice ) then
! This is a very small term. It may not be needed as it is doing similar 
! job as realidw
! pihtf = dt*1.44e-12;  for tmp1 = 20, den=0.8 and ql = 1.e-4
! max(pihtf) ~ dt*1.e-6;  for tmp1 = 40, den=1.0 and ql = 1.e-3
             tmp1 = tice - tzk(k)
            frez = 3.3333e-10*(exp(0.66*tmp1) - 1.)*den(k)*qlk(k)*qlk(k) 
            frez = min(qlk(k), tmp1/icpk(k), dts*frez)
           qlk(k) = qlk(k) - frez
           qik(k) = qik(k) + frez
           tzk(k) = tzk(k) + frez*icpk(k)
        endif
    endif

! * pimlt: instant melting of cloud ice
    if( tzk(k) > tice .and. qik(k) > 0. ) then
          melt = min(qik(k), (tzk(k)-tice)/icpk(k))
        qlk(k) = qlk(k) + melt      ! cloud ice --> cloud water
        qik(k) = qik(k) - melt
        tzk(k) = tzk(k) - melt*icpk(k)
!--------------------------------------------------------------
! Melting starts if T >1 C
!       melt = min(qik(k), dim(tzk(k),tice1)/icpk(k))
!
!       if( qik(k) > i2r_crt ) then
!           qrk(k) = qrk(k) + melt      ! ice --> rain
!       else
!           qlk(k) = qlk(k) + melt      ! ice --> cloud water
!       endif 
! At this stage it is possible to have cloud ice & T > Tice
! Convert remaining ice to snow
!       qsk(k) = qsk(k) + qik(k) - melt ! ice --> snow
!       qik(k) = 0.
!--------------------------------------------------------------
    endif

 enddo

 do 1000 k=ktop, kbot

   tz = tzk(k)
   qv = qvk(k)
   ql = qlk(k)
   qi = qik(k)
   qr = qrk(k)
   qs = qsk(k)
   qg = qgk(k)

!--------------------------------------
! *** Split-micro_physics_processes ***
!--------------------------------------
! Zetac: excluded (from LFO83) term: psdep
! pgwet removed by SJL (reducing complexity)
! piacr added to increase ice production

   pgacr = 0.
   pgacw = 0.
   tc = tz-tice

if ( tc > 1.0 ) then

!---------------------------------------------
! * Melting of snow (> 1 C) and graupel (> 2C)
!---------------------------------------------

     dqs0 = ces0/p1(k) - qv

     if( qs>qrmin ) then ! melting of snow into rain
! * accretion: cloud water --> snow
! only rate is used (for snow melt) since tc > 0.
        if( ql>qrmin ) then
!           factor = denfac(k)*csacw*(qs*den(k))**0.8125
            factor = denfac(k)*csacw*exp(0.8125*log(qs*den(k)))
             psacw = factor/(1.+dts*factor)*ql     ! rate
        else
             psacw = 0.
        endif

        if ( qr>qrmin ) then
! * accretion: melted snow --> rain:
             psacr = min(acr3d(vts(k), vtr(k), qr, qs, csacr, c3sacr, den(k)), qr*rdts)
! * accretion: snow --> rain
             pracs = acr3d(vtr(k), vts(k), qs, qr, cracs, c3racs, den(k))
        else
             psacr = 0.
             pracs = 0.
        endif

! * Snow melt (due to rain accretion): snow --> rain
        psmlt = max(0., smlt(tc, dqs0, qs*den(k), psacw, psacr, csmlt, den(k), denfac(k)))

! Total snow sink:
        sink = min(qs, dts*(psmlt+pracs), tc/icpk(k))
        qs = qs - sink
        qr = qr + sink
        tz = tz - sink*icpk(k)
     endif   ! snow existed

     tc = tz-tice
     if ( qg>qrmin .and. tc>2. ) then
         if ( qr>qrmin ) then
! * accretion: rain --> graupel
              pgacr = min(acr3d(vtg(k), vtr(k), qr, qg, cgacr, c3gacr, den(k)), rdts*qr)
         endif

         qden = qg*den(k)
         if( ql>qrmin ) then
! * accretion: cloud water --> graupel
!            factor = cgacw/sqrt(den(k))*(qg*den(k))**0.875
             factor = cgacw*qden/sqrt(den(k)*sqrt(sqrt(qden)))
              pgacw = factor/(1.+dts*factor) * ql  ! rate
         endif

! * melting: graupel --> rain
         pgmlt = dts*gmlt(tc, dqs0, qden, pgacw, pgacr, cgmlt, den(k))
         pgmlt = min( max(0., pgmlt), qg, max(0.,tc/icpk(k)) )
            qg = qg - pgmlt 
            qr = qr + pgmlt 
            tz = tz - pgmlt*icpk(k)
     endif   ! graupel existed



elseif( tc < 0. ) then 

!------------------
! Cloud water sink:
!------------------

  if( ql>qrmin ) then
! cloud water --> ice
      if ( qi>qrmin ) then
           pidw = min(ql, -tc/icpk(k), dts*realidw(tc, qi, den(k)))
           qi = qi + pidw
           ql = ql - pidw
           tz = tz + pidw*icpk(k)
           tc = tz - tice
      endif
! cloud water --> Snow
      if( qs>qrmin ) then
!         factor = dts*denfac(k)*csacw*(qs*den(k))**0.8125
          factor = dts*denfac(k)*csacw*exp(0.8125*log(qs*den(k)))
          psacw = min( factor/(1.+factor)*ql, -tc/icpk(k) )
          qs = qs + psacw
          ql = ql - psacw
          tz = tz + psacw*icpk(k)
          tc = tz - tice
      endif
  endif  ! (significant) cloud water existed
 


!------------------
! Cloud ice proc:
!------------------

  if ( qi>qrmin ) then

! * accretion (pacr): cloud ice --> snow
     if ( qs>qrmin )  then
#ifdef ZETAC_PSACI
! check Lin Eq. 22 has temperature dependency
! The following is from the "Lin Micro-physics" in Zetac
          factor = dts*denfac(k)*csaci*(qs*den(k))**0.8125
#else
! Eq(10) in HDC 2004, MWR
!  Gamma(3.41) = 3.0133  (combined into the constant: 27.737)
! Optimized form:
            n0s = 2.e6*exp(-0.12*tc)
          lamda = (qs*den(k))/(pie*rhos*n0s)
!        factor = dts*denfac(k)*27.737*n0s*exp(0.05*tc)*lamda**0.8525
         factor = dts*denfac(k)*27.737*n0s*exp(0.05*tc + 0.8525*log(lamda))
#endif
          psaci = factor/(1.+factor) * qi
     else
          psaci = 0.
     endif

! * autoconversion: cloud ice --> snow
! (according to LFO 1983: Eq. 21 solved implicitly)
! Other than the Bergeron process this is the only sink to cloud ice if snow does not pre-exist
! and the air is saturated.

   qi_tho = qi0_crt / den(k)

   if ( qi>qi_tho ) then
        factor = dts*c_psaut*exp(0.025*tc)
         psaut = factor/(1.+factor)*(qi-qi_tho)
   else
         psaut = 0.
   endif

!---------------------------------------------
! "Bergeron" processes -- psfw and psfi
!---------------------------------------------
! psfi: cloud ice   --> snow
! psfw: cloud water --> snow

     if ( ql>qrmin ) then
          call bergrn(dts, tc, ql, qi, den(k), psfw, psfi)
     else
          psfw = 0.
          psfi = 0.
     endif

     sink = min(qi, psaci+psaut+psfi)
       qi = qi - sink
       ql = ql - psfw
       qs = qs + sink + psfw
       tz = tz + psfw*icpk(k)

! * accretion: cloud ice --> graupel
      if ( qg>qrmin .and. qi>qrmin ) then
!          factor = dts*cgaci/sqrt(den(k))*(qg*den(k))**0.875
             qden = qg*den(k)
           factor = dts*cgaci*qden/sqrt(den(k)*sqrt(sqrt(qden)))
            pgaci = factor/(1.+factor)*qi
               qi = qi - pgaci
               qg = qg + pgaci
      endif

  endif  ! cloud ice existed
 

!----------------
! Cold-Rain proc:
!----------------
! rain to ice, snow, graupel processes:

  if ( qr>qrmin ) then

! * piacr: accretion of rain by cloud ice [lfo 26]
! rain --> ice  (factor somewhat arbitrary; totally tunable)
! The value of "factor" needs to be near order(1) to be really effective

       if ( qi>qrmin ) then   ! if ice existed: rain --> ice
            factor = dts*denfac(k)*qi * c_piacr
            piacr = factor/(1.+factor)*qr
            piacr = min( piacr, dim(tice,tz)/icpk(k) )
               qr = qr - piacr
               qi = qi + piacr
               tz = tz + piacr*icpk(k)
       endif

! * accretion of rain by snow; rain --> snow
       if ( qs>qrmin .and. qr>qrmin ) then   ! if snow
            psacr = dts*acr3d(vts(k), vtr(k), qr, qs, csacr, c3sacr, den(k))
            psacr = min( psacr, qr,  dim(tice,tz)/icpk(k) )
               qr = qr - psacr
               qs = qs + psacr
               tz = tz + psacr*icpk(k)
       endif
 
       tc = tz-tice
       if ( qr>qrmin .and. tc<0. )  then

! * accretion: accretion of cloud ice by rain to produce snow or graupel
! (LFO: produces snow or graupel; cloud ice sink.. via psacr & pgfr)
!           factor = dts*denfac(k)*craci*(qr*den(k))**0.95
            factor = dts*denfac(k)*craci*exp(0.95*log(qr*den(k)))
             praci = factor/(1.+factor)*qi   ! check praci form
             if ( qr > qr0_crt ) then
                  qg = qg + praci
             else
                  qs = qs + praci
             endif
             qi = qi - praci

! * rain freezing --> graupel/hail
!           pgfr = dts*cgfr(1)*(exp(-cgfr(2)*tc)-1.)*(qr*den(k))**1.75/den(k)
            qden = qr*den(k)
            pgfr = dts*cgfr(1)*(exp(-cgfr(2)*tc)-1.)*qden*qden/(sqrt(sqrt(qden))*den(k))
            pgfr = min(pgfr, -tc/icpk(k), qr)
              qr = qr - pgfr
              qg = qg + pgfr
              tz = tz + pgfr*icpk(k)
       endif
  endif
 


  if ( qs>qrmin ) then

!-----------------------
! * sublimation of snow:
!-----------------------
! snow <--> vapor  two-way conversion
        qsi = qs1d(tz, p1(k), dqsdt)
       qden = qs*den(k)
       tmp1 = exp(0.65625*log(qden))
        tsq = tz*tz
      pssub = cssub(1)*tsq*(cssub(2)*sqrt(qden) + &
              cssub(3)*tmp1*sqrt(denfac(k)))/(cssub(4)*tsq+cssub(5)*qsi*den(k))

      pssub = (qsi-qv)*min(dts*pssub, 1./(1.+tcpk(k)*dqsdt))

      if ( pssub > 0. ) then
           pssub = min(pssub, qs)
      else
           pssub = max(pssub, (tz-tice)/tcpk(k))
      endif

      qs = qs - pssub 
      qv = qv + pssub 
      tz = tz - pssub*tcpk(k)

!-----------------------------------
! * Autoconversion Snow --> graupel
!-----------------------------------
     if ( qs>qs0_crt ) then
          factor = dts*1.e-3*exp(0.09*(tz-tice))
           pgaut = factor/(1.+factor)*(qs-qs0_crt)
              qs = qs - pgaut
              qg = qg + pgaut
     endif

  endif    ! snow existed

 
!--------------------------
! Graupel production terms:
!--------------------------

  if ( qg>qrmin ) then

! * accretion: rain --> graupel
     if ( qr>qrmin ) then 
          pgacr = min(dts*acr3d(vtg(k), vtr(k), qr, qg, cgacr, c3gacr, den(k)), qr)
             qr = qr - pgacr
     endif

! * accretion: cloud water --> graupel
     if( ql>qrmin ) then
!        factor = dts*cgacw/sqrt(den(k))*(qg*den(k))**0.875
           qden = qg*den(k)
         factor = dts*cgacw*qden/sqrt(den(k)*sqrt(sqrt(qden)))
          pgacw = factor/(1.+factor)*ql
             ql = ql - pgacw
     endif

! * accretion: snow --> graupel
     if( qs>qrmin ) then
         pgacs = min(dts*acr3d(vtg(k), vts(k), qs, qg, cgacs, c3gacs, den(k)), qs)
         qs = qs - pgacs
         qg = qg + pgacs
     endif

     qg = qg +  pgacr+pgacw
     tz = tz + (pgacr+pgacw)*icpk(k)

! * graupel sublimation (unsaturated only)
! make it 2-way?
     qsi = qs1d(tz, p1(k), dqsdt)
      dq = (qsi-qv)/(1.+tcpk(k)*dqsdt)

     if ( dq > 0. ) then
           qden = qg*den(k)
           tmp1 = exp(0.6875*log(qden))
            tsq = tz*tz
          pgsub = cgsub(1)*tsq*dq*( cgsub(2)*sqrt(qden) +  &
                  cgsub(3)*tmp1/sqrt(sqrt(den(k))) ) / (cgsub(4)*tsq + cgsub(5)*qsi*den(k))
          pgsub = min(dts*pgsub, qg, dq)
             qg = qg - pgsub
             qv = qv + pgsub 
             tz = tz - pgsub*tcpk(k)
     endif
 endif    ! graupel existed

endif   ! end ice-physics 

!-----------------------------
! * Soft saturation adjustment
!-----------------------------
     qvk(k) = qv
     qlk(k) = ql
     qik(k) = qi
     tzk(k) = tz

!-------------------
! * Pre-conditioner:
!-------------------
   clouds = ql + qi
      qpz = qv + clouds
      tin = tz - ( lcpk(k)*clouds + icpk(k)*qi )  ! minimum possible temperature
    rh(k) = qpz*p1(k)/(eps*es2_table(tin))        ! ice-phase

! Instant evaporation of all clouds is allowed if RH<0.85
! Otherwise, clouds < qc_crt is required
    if ( rh(k)<0.85 .or. (clouds<qc_crt .and. rh(k)<1.0) ) then  
         tzk(k) = tin
         qvk(k) = qpz
         qlk(k) = 0.
         qik(k) = 0.
         goto 123
    endif


! * Cloud water <---> water vapor
! Time split the faster liquid-phase process: perform half time step before the ice-phase
   if ( tzk(k) > t40 ) then
         qsw = ws1d(tzk(k), p1(k), dwsdt)
        sink = f_exp*(qsw-qvk(k))/(1.+lcpk(k)*dwsdt)
        sink = min(qlk(k), sink)
      qvk(k) = qvk(k) + sink
      qlk(k) = qlk(k) - sink
      tzk(k) = tzk(k) - sink*lcpk(k)
   endif

! * Cloud ice  <----> water vapor

   tc = tzk(k) - tice

   if ( tc < 0.0 ) then

       qsi = qs1d(tzk(k), p1(k), dqsdt)
      sink = (qsi-qvk(k))/(1.+tcpk(k)*dqsdt)

!------------------------------------------------
! Initiation of cloud ice (pigen): following WSM6
!------------------------------------------------
      if ( sink<0.0 .and. qi<1.e-5 .and. tc<-1.0 ) then
!     tmp1 = 4.92e-11 * (1.e3*exp(-0.1*tc))**1.33
!     tmp1 = 4.92e-11 * exp(1.33*log(1.e3*exp(-0.1*tc)))
         qimin = 4.92e-11*exp(9.1873145 - 0.133*tc)/den(k)
            dq = qik(k) - qimin
      else
            dq = 1.
      endif

      if ( dq < 0.0 ) then
           sink = max(dq, sink, tc/tcpk(k))
      else

! pidep: Deposition/Sublimation rate of ice
        if ( qik(k) > qrmin ) then
!          pidep = 47.6*sqrt(5.38e7)*(qsi-qvk(k))*(den(k)*qik(k))**(7./8.)
! Optimization:
             qden = qik(k)*den(k)
            pidep = 349138.78*(qsi-qvk(k))*qden/sqrt(sqrt(sqrt(qden)))
! Dudhia 1989: page 3103 Eq (B7) and (B8)
               ai = den(k)*lats*lats/(0.0243*rvgas*tzk(k)**2)
            pidep = dts*pidep/(ai*qsi + 4.42478e4)
            if ( sink > 0. ) then
                 sink = min(qik(k), sink, pidep)       ! sublimation of ice
            else
                 sink = max(sink, pidep, tc/tcpk(k))   ! growth of ice c. by deposition
            endif
        else
            sink = 0.
        endif
      endif
      qvk(k) = qvk(k) + sink
      qik(k) = qik(k) - sink
      tzk(k) = tzk(k) - sink*tcpk(k)
   endif

! Finish the remaining half time step:
! * Cloud water <---> water vapor

   if ( tzk(k) > t40 ) then
         qsw = ws1d(tzk(k), p1(k), dwsdt)
        sink = f_exp*(qsw-qvk(k))/(1.+lcpk(k)*dwsdt)
        sink = min(qlk(k), sink)
      qvk(k) = qvk(k) + sink
      qlk(k) = qlk(k) - sink
      tzk(k) = tzk(k) - sink*lcpk(k)
   endif

!*****************
! * RH computation
!*****************
   clouds = qlk(k) + qik(k) + qs  ! snow field is included for radiation
   qpz    = clouds + qvk(k)

   if ( tzk(k) <= tice ) qsi = iqsat(tzk(k), p1(k))
   if ( tzk(k) >= t40  ) qsw = wqsat(tzk(k), p1(k))

   if( tzk(k) <= t40 ) then
       rh(k) = qpz / qsi
   elseif ( tzk(k) >= tice ) then
       rh(k) = qpz / qsw
   else
! mixed phase:
       if( clouds > 1.e-5 ) then
           rqi = (qik(k)+qs) / clouds
       else
           rqi = max(0., min(1., (tice-tz)/30.))
       endif
       rh(k) = qpz / (rqi*qsi + (1.-rqi)*qsw)
   end if

123   continue

 qrk(k) = qr
 qsk(k) = qs
 qgk(k) = qg

1000 continue

 end subroutine icloud



 subroutine terminal_fall(dtm, ktop, kbot, tz, qv, ql, qr, qg, qs, qi, pm, dz, dp,  &
                          den, denfac, vtg, vts, vti, rh_p, r1, g1, s1, i1)

! lagrangian control-volume method:

 real,    intent(in):: dtm                    ! time step (s)
 integer, intent(in):: ktop, kbot
 real,    intent(in), dimension(ktop:kbot):: dp, vtg, vts, vti, pm, rh_p,  &
                                             den, denfac
 real,    intent(inout), dimension(ktop:kbot):: dz, qv, ql, qr, qg, qs, qi, tz
 real,    intent(out):: r1, g1, s1, i1
! local:
 real, dimension(ktop:kbot+1):: ze, zt
 real:: qsat, dqsdt, dt5, melt, evap, dtime
 real:: factor, frac
 real:: qden, accr, tsq, tmp1, precip
 real, dimension(ktop:kbot):: lcpk, icpk
 real:: zs = 0.
 integer k, k0, m
 logical no_fall

  do k=ktop,kbot
!       tmp1 = cp - rdgas*ptop/pm(k)
!    lcpk(k) = latv / tmp1
!    icpk(k) = lati / tmp1
     lcpk(k) = lcp
     icpk(k) = icp
  enddo

  dt5 = 0.5*dtm

! find melting level; tice+2.
  k0 = kbot
  do k=ktop, kbot-1
     if ( tz(k) > tice2 ) then
          k0 = k
          go to 11
     endif
  enddo
11  continue

!-----
! ice:
!-----

  ze(kbot+1) = zs
  do k=kbot, ktop, -1
     ze(k) = ze(k+1) - dz(k)  ! dz<0
  enddo

  zt(ktop) = ze(ktop)

  call check_column(ktop, kbot, qi, no_fall)

  if ( vi_fac < 1.e-5 .or. no_fall ) then
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
             if ( zt(k)<ze(m) .and. tz(m)>tice2 ) then
                  dtime = min(1., (ze(m)-ze(m+1))/(max(vmin,vti(k))*tau_i))
                   melt = min(qi(k)*dp(k)/dp(m), dtime*(tz(m)-tice2)/icpk(m))
                  ql(m) = ql(m) + melt         ! melt into cloud water
                  tz(m) = tz(m) - melt*icpk(m)
                  qi(k) = qi(k) - melt*dp(m)/dp(k)
!            elseif ( zt(k+1)<ze(m) .and. zt(k+1)>ze(m+1) ) then
!                  frac = (ze(m)-zt(k+1))/(zt(k)-zt(k+1))
!                  melt = min( frac*qi(k)*dp(k)/dp(m), dim(tz(m),tice2)/icpk(m) )
!                 ql(m) = ql(m) + melt          ! melt into cloud water
!                 tz(m) = tz(m) - melt*icpk(m)
!                 qi(k) = qi(k) - melt*dp(m)/dp(k)
             endif
          enddo
     endif
  enddo
  endif

  if ( use_ppm ) then
       call lagrangian_fall_ppm(ktop, kbot, zs, ze, zt, dp, qi, i1)
  else
       call lagrangian_fall_pcm(ktop, kbot, zs, ze, zt, dp, qi, i1)
  endif

  endif

!--------------------------------------------
! melting of falling snow (qs) into rain(qr)
!--------------------------------------------
  r1 = 0.

  call check_column(ktop, kbot, qs, no_fall)

  if ( no_fall ) then
       s1 = 0.
  else

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
!            if ( zt(k)<ze(m) ) then    ! the top of the c-v is in the layer below
             if ( zt(k)<ze(m+1) .and. tz(m)>tice2 ) then
                  dtime = min(1., dtime/tau_s)
                   melt = min(qs(k)*dp(k)/dp(m), dtime*(tz(m)-tice2)/icpk(m))
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

  if ( use_ppm ) then
       call lagrangian_fall_ppm(ktop, kbot, zs, ze, zt, dp, qs, s1)
  else
       call lagrangian_fall_pcm(ktop, kbot, zs, ze, zt, dp, qs, s1)
  endif
  endif

!----------------------------------------------
! melting of falling graupel (qg) into rain(qr)
!----------------------------------------------
  call check_column(ktop, kbot, qg, no_fall)

  if ( no_fall ) then
       g1 = 0.
  else
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
             if ( zt(k)<ze(m+1) .and. tz(m)>tice2 ) then
                  dtime = min(1., dtime/tau_g)
                   melt = min(qg(k)*dp(k)/dp(m), dtime*(tz(m)-tice2)/icpk(m))
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

  if ( use_ppm ) then
       call lagrangian_fall_ppm(ktop, kbot, zs, ze, zt, dp, qg, g1)
  else
       call lagrangian_fall_pcm(ktop, kbot, zs, ze, zt, dp, qg, g1)
  endif
  endif


 end subroutine terminal_fall


 subroutine check_column(ktop, kbot, q, no_fall)
 integer, intent(in):: ktop, kbot
 real,    intent(in):: q(ktop:kbot)
 logical, intent(out):: no_fall
! local:
 integer k

 no_fall = .true.
 do k=ktop, kbot
    if ( q(k) > qrmin ) then
         no_fall = .false.
         exit
    endif
 enddo

 end subroutine check_column


 subroutine lagrangian_fall_pcm(ktop, kbot, zs, ze, zt, dp, q, precip)
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
     qm2(k) = 0.
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

 end subroutine lagrangian_fall_pcm



 subroutine lagrangian_fall_ppm(ktop, kbot, zs, ze, zt, dp, q, precip)
 integer, intent(in):: ktop, kbot
 real,    intent(in):: zs
 real,    intent(in), dimension(ktop:kbot):: dp
 real,    intent(in), dimension(ktop:kbot+1):: ze, zt
 real,    intent(inout), dimension(ktop:kbot):: q
 real,    intent(out):: precip
! local:
 real, dimension(ktop:kbot):: qm0, qm1, qm2, dz
 real a4(4,ktop:kbot)
 real pl, pr, delz, esl
 integer k, k0, n, m
 real, parameter:: r3 = 1./3., r23 = 2./3.

! density:
  do k=ktop,kbot
      dz(k) = zt(k) - zt(k+1)      ! note: dz is positive
     qm0(k) = q(k)*dp(k)
     qm1(k) = qm0(k) / dz(k)
     qm2(k) = 0.
     a4(1,k) = qm1(k)
  enddo

! Construct qm1 profile with zt as coordinate

   call c1_profile(a4(1,ktop), dz(ktop), kbot-ktop+1, .true.)

   k0 = ktop
   do k=ktop,kbot
      do n=k0,kbot
      if(ze(k) <= zt(n) .and. ze(k) >= zt(n+1)) then
         pl = (zt(n)-ze(k)) / dz(n)
         if( zt(n+1) <= ze(k+1) ) then
!                          entire new grid is within the original grid
                pr = (zt(n)-ze(k+1)) / dz(n)
            qm2(k) = a4(2,n) + 0.5*(a4(4,n)+a4(3,n)-a4(2,n))*(pr+pl) -  &
                     a4(4,n)*r3*(pr*(pr+pl)+pl**2)
            qm2(k) = qm2(k)*(ze(k)-ze(k+1))
            k0 = n
            goto 555
         else
            qm2(k) = (ze(k)-zt(n+1)) * (a4(2,n)+0.5*(a4(4,n)+   &
                      a4(3,n)-a4(2,n))*(1.+pl) - a4(4,n)*( r3*(1.+pl*(1.+pl))) )
            if ( n<kbot ) then
               do m=n+1,kbot
!                                        locate the bottom edge: ze(k+1)
                  if( ze(k+1) < zt(m+1) ) then
                     qm2(k) = qm2(k) + q(m)*dp(m)
                  else
                     delz = zt(m) - ze(k+1)
                      esl = delz / dz(m)
                     qm2(k) = qm2(k) + delz*( a4(2,m) + 0.5*esl*        &
                             (a4(3,m)-a4(2,m)+a4(4,m)*(1.-r23*esl)) )
                     k0 = m
                     goto 555
                  endif
               enddo
            endif
            goto 555
         endif
      endif
      enddo
555 continue
   enddo

   precip = 0.

   do k=ktop,kbot
      precip = precip + qm0(k) - qm2(k)
   enddo
!  precip = max(0., precip)

   do k=ktop,kbot
      q(k) = qm2(k) / dp(k)
   enddo

 end subroutine lagrangian_fall_ppm


 subroutine c1_profile(a4, del, km, do_huynh)
 integer, intent(in):: km      ! vertical dimension
 real   , intent(in):: del(km)
 logical, intent(in):: do_huynh
 real , intent(inout):: a4(4,km)
!-----------------------------------------------------------------------
 real  gam(km)
 real  q(km+1)
 real   d4, bet, a_bot, grat, pmp, lac
 integer k

     grat = del(2) / del(1)   ! grid ratio
      bet = grat*(grat+0.5)
     q(1) = (2.*grat*(grat+1.)*a4(1,1)+a4(1,2)) / bet
   gam(1) = ( 1. + grat*(grat+1.5) ) / bet

  do k=2,km
      d4 = del(k-1) / del(k)
     bet =  2. + 2.*d4 - gam(k-1)
     q(k) = (3.*(a4(1,k-1)+d4*a4(1,k))-q(k-1))/bet
     gam(k) = d4 / bet
  enddo
 
       a_bot = 1. + d4*(d4+1.5)
     q(km+1) = (2.*d4*(d4+1.)*a4(1,km)+a4(1,km-1)-a_bot*q(km))  &
             / ( d4*(d4+0.5) - a_bot*gam(km) )

  do k=km,1,-1
     q(k) = q(k) - gam(k)*q(k+1)
  enddo

!------------------
! Apply constraints
!------------------
  do k=2,km
     gam(k) = a4(1,k) - a4(1,k-1)
  enddo

! Apply large-scale constraints to ALL fields if not local max/min
! Top:

  if ( (q(2)-q(1))*(q(3)-q(2))>0. ) then
       q(2) = min( q(2), max(a4(1,1), a4(1,2)) )
       q(2) = max( q(2), min(a4(1,1), a4(1,2)) )
  else
       q(2) = max(0., q(2))
  endif

! Interior:
  do k=3,km-1
     if ( gam(k-1)*gam(k+1)>0. ) then
          q(k) = min( q(k), max(a4(1,k-1),a4(1,k)) )
          q(k) = max( q(k), min(a4(1,k-1),a4(1,k)) )
     else
          q(k) = max(0., q(k))
     endif
  enddo

! Bottom:
  if ( (q(km)-q(km-1))*(q(km+1)-q(km))>0. ) then
       q(km) = min( q(km), max(a4(1,km-1), a4(1,km)) )
       q(km) = max( q(km), min(a4(1,km-1), a4(1,km)) )
  else
       q(km) = max(0., q(km))
  endif

  do k=1,km
     a4(2,k) = q(k  )
     a4(3,k) = q(k+1)
  enddo

!-----------------------------------------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
!-----------------------------------------------------------
! Top & bot surfaces
  a4(2, 1) = max(0., a4(2, 1))
  a4(3,km) = max(0., a4(3,km))

!-------------------------------------
! Huynh's 2nd constraint for interior:
!-------------------------------------
  if ( do_huynh ) then
     do k=3,km-2
! Left  edges
        pmp = a4(1,k) - 2.*gam(k+1)
        lac = pmp + 1.5*gam(k+2)
        a4(2,k) = min( max(a4(2,k), min(a4(1,k), pmp, lac)),   &
                                    max(a4(1,k), pmp, lac) )
! Right edges
        pmp = a4(1,k) + 2.*gam(k)
        lac = pmp - 1.5*gam(k-1)
        a4(3,k) = min( max(a4(3,k), min(a4(1,k), pmp, lac)),    &
                                    max(a4(1,k), pmp, lac) )
     enddo
  endif

  do k=1,km
     a4(4,k) = 3.*(2.*a4(1,k) - (a4(2,k)+a4(3,k)))
  enddo

  call c1_limiters(km, a4)

 end subroutine c1_profile



 subroutine c1_limiters(km, a4)
 integer, intent(in) :: km
 real, intent(inout) :: a4(4,km)   ! PPM array
! !LOCAL VARIABLES:
 real  da1, da2, a6da, fmin
 real, parameter:: r12 = 1./12.
 integer k

! Positive definite constraint

 do k=1,km
 if( abs(a4(3,k)-a4(2,k)) < -a4(4,k) ) then
     fmin = a4(1,k)+0.25*(a4(3,k)-a4(2,k))**2/a4(4,k)+a4(4,k)*r12
     if( fmin < 0. ) then
         if( a4(1,k)<a4(3,k) .and. a4(1,k)<a4(2,k) ) then
             a4(3,k) = a4(1,k)
             a4(2,k) = a4(1,k)
             a4(4,k) = 0.
         elseif( a4(3,k) > a4(2,k) ) then
             a4(4,k) = 3.*(a4(2,k)-a4(1,k))
             a4(3,k) = a4(2,k) - a4(4,k)
         else
             a4(4,k) = 3.*(a4(3,k)-a4(1,k))
             a4(2,k) = a4(3,k) - a4(4,k)
         endif
     endif
 endif
 enddo

 end subroutine c1_limiters


 subroutine fall_speed(ktop, kbot, den, qs, qi, qg, vts, vti, vtg)
 integer, intent(in)                     :: ktop, kbot
 real, intent(in ), dimension(ktop:kbot) :: den, qs, qi, qg
 real, intent(out), dimension(ktop:kbot) :: vts, vti, vtg
! fall velocity constants:
 real :: vcons = 6.6280504, vcong = 87.2382675,         vconi = 3.29
 real, parameter :: thg = 1.e-9, ths=1.e-9
!-----------------------------------------------------------------------
! marshall-palmer constants
!-----------------------------------------------------------------------
 real :: norms = 942477796.076938, &
         normg =  5026548245.74367
 real :: rhof, rho0, qden
 integer:: k
!-----------------------------------------------------------------------
! marshall-palmer formula
!-----------------------------------------------------------------------

! try the local air density -- for global model; the true value could be
! much smaller than sfcrho over high mountains

  if ( den_ref < 0. ) then
       rho0 = -den_ref*den(kbot) 
  else
       rho0 = den_ref   ! default=1.2
  endif

   do k=ktop, kbot
        rhof = sqrt( min(100., rho0/den(k)) )
! snow:
!     vts(k) = vs_fac*vcons*rhof*(dim(qs(k)*den(k),ths)/norms)**0.0625
      vts(k) = vs_fac*vcons*rhof*exp(0.0625*log(dim(qs(k)*den(k),ths)/norms))
! graupel:
!     vtg(k) = vg_fac*max(vmin, vcong*rhof*(dim(qg(k)*den(k),thg)/normg)**0.125)
      vtg(k) = vg_fac*max(vmin, vcong*rhof*sqrt(sqrt(sqrt(dim(qg(k)*den(k),thg)/normg))) )
! ice:
!     vti(k) = vi_fac*vconi*rhof*dim(qi(k)*den(k),thi)**0.16
      vti(k) = vi_fac*vconi*rhof*exp(0.16*log(dim(qi(k)*den(k),thi)))
   enddo

 end subroutine fall_speed



 real function realidw(tc, qik, rho )
! Production of ice from super cooled cloud water
      real, intent(in):: tc, qik,  rho
      real :: fmi, a1, a2, tc1, a1t, a2t, expt
      dimension a1t(0:31),a2t(0:31)
! fu *******************************************************************
!     note -- a2t is identical to a2t in subroutine bergrn, but a1t
!     is multiplied by a factor of 1.e-5.
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

      tc1 = max(tc,-30.0)

      a1  = (a1t(-int(tc1))-a1t(-int(tc1)+1))*(tc1-int(tc1)+1.0)+ &
             a1t(-int(tc1)+1)
      a2  = (a2t(-int(tc1))-a2t(-int(tc1)+1))*(tc1-int(tc1)+1.0)+ &
             a2t(-int(tc1)+1)

      expt = exp(-0.6*tc)
       fmi = rho*qik/expt * 1.E5
      realidw = expt*a1*fmi**a2/rho

 end function realidw


 subroutine bergrn1(dtbg, tc, ql, qi, rho, psfw, psfi)
      real, intent(in):: dtbg, rho, qi, ql, tc
      real, intent(out):: psfi, psfw
      real, parameter:: rmi40=2.46e-7
      real:: a21, a2,  a1, a1t, a2t
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
 
  if ( tc > -1.0 ) then
       a1 = a1t(1)
       a2 = a2t(1)
       a21 = 1. - a2
       psfi = min(qi, dtbg*qi*a1*a21/(rmi50**a21 - rmi40**a21))
       psfw = min(ql, psfi*c1brg*(a1*rmi50**a2 + ql*rho*c2brg))
  elseif( tc > -30. ) then
       a1 = (a1t(-int(tc))-a1t(-int(tc)+1))*(tc-int(tc)+1.0)+ &
             a1t(-int(tc)+1)
       a2 = (a2t(-int(tc))-a2t(-int(tc)+1))*(tc-int(tc)+1.0)+ &
             a2t(-int(tc)+1)
       a21 = 1. - a2
! psfi: cloud ice --> snow
      psfi = min(qi, dtbg*qi*a1*a21/(rmi50**a21 - rmi40**a21))
! psfw: cloud wat --> snow
      psfw = min(ql, psfi*c1brg*(a1*rmi50**a2 + ql*rho*c2brg))
  else
! No snow production if T < -30C
      psfi = 0.
      psfw = 0.
  endif


 end subroutine bergrn1


 subroutine bergrn(dtbg, tc, ql, qi, rho, psfw, psfi)
      real, intent(in):: dtbg, rho, qi, ql, tc
      real, intent(out):: psfi, psfw
      real, parameter:: rmi40=2.46e-7
      real:: a21, a2,  a1, a1t, a2t, tc1
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
 
      tc1 = max(tc, -30.)

      if ( tc1 > -1.0 ) then
           a1 = a1t(1)
           a2 = a2t(1)
      else
           a1 = (a1t(-int(tc1))-a1t(-int(tc1)+1))*(tc1-int(tc1)+1.0)+ &
                 a1t(-int(tc1)+1)
           a2 = (a2t(-int(tc1))-a2t(-int(tc1)+1))*(tc1-int(tc1)+1.0)+ &
                 a2t(-int(tc1)+1)
      endif
!     note:  mks units, ui50=1.0 m/sec, eiw=1.0
      a21 = 1. - a2
! The snow production here is independent of amount of ql available?
! psfi: cloud ice --> snow
      psfi = dtbg*qi*a1*a21/(rmi50**a21 - rmi40**a21)
!----------------------------------------
! ramp to zero linearly from -30 to -40 C
!----------------------------------------
      if ( tc<-30. )  psfi = psfi * max(0., 1.+(30.+tc)/10.)
      psfi = min(qi, psfi)
! psfw: cloud water --> snow
      psfw = min(ql, psfi*c1brg*(a1*rmi50**a2 + ql*rho*c2brg))

 end subroutine bergrn


 subroutine setupm

 real :: gcon, cd, scm3, pisq, act(8), acc(3)
 real :: vdifu, tcond
 real :: visk
 real :: ch2o, hltf
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
 real :: rnzr, rnzs, rnzg, rhos, rhog
 data alin, clin  /842.0, 4.80/
 data rnzr /8.0e6/  ! lin83
 data rnzs /3.0e6/  ! lin83
 data rnzg /4.0e6/  ! rh84
 data rhos /0.1e3/  ! lin83    (snow density; 1/10 of water)
 data rhog /0.4e3/  ! rh84     (graupel density)
 data acc/5.0,2.0,0.5/

 real den_rc
 integer :: k, i

      pie = 4.*atan(1.0)

! S. Klein's formular (EQ 16) from AM2
      fac_rc = (4./3.)*pie*rhor*rthresh**3
      den_rc = fac_rc * ccn_o*1.e6
      if(master) write(*,*) 'MP: rthresh=', rthresh, 'vi_fac=', vi_fac
      if(master) write(*,*) 'MP: for ccn_o=', ccn_o, 'ql_rc=', den_rc
      den_rc = fac_rc * ccn_l*1.e6
      if(master) write(*,*) 'MP: for ccn_l=', ccn_l, 'ql_rc=', den_rc

      vdifu=2.11e-5
      tcond=2.36e-2

      visk=1.259e-5
      hlts=2.8336e6
      hltc=2.5e6
      hltf=3.336e5

      ch2o=4.1855e3
      rmi50=3.84e-6
!     rmi40=2.46e-7
      ri50=1.e-4

      pisq = pie*pie
      scm3 = (visk/vdifu)**(1./3.)
!
      cracs = pisq*rnzr*rnzs*rhos
      csacr = pisq*rnzr*rnzs*rhor
      cgacr = pisq*rnzr*rnzg*rhor
      cgacs = pisq*rnzg*rnzs*rhos
      cgacs = cgacs*c_pgacs
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
      csacw = pie*rnzs*clin*gam325/(4.*act(1)**0.8125)
!     ciacr = pisq*rhor*rnzr*alin*gam680/(1.0056e-11*act(2)**1.7)
      craci = pie*rnzr*alin*gam380/(4.*act(2)**0.95)
      csaci = csacw * c_psaci
!
      cgacw = pie*rnzg*gam350*gcon/(4.*act(6)**0.875)
      cgaci = cgacw*0.1
!
      cracw = craci            ! cracw= 3.27206196043822
!     cracw = c_cracw * cracw
!
!     subl and revp:  five constants for three separate processes
!
      cssub(1) = 2.*pie*vdifu*tcond*rvgas*rnzs
      cgsub(1) = 2.*pie*vdifu*tcond*rvgas*rnzg
      crevp(1) = 2.*pie*vdifu*tcond*rvgas*rnzr
      cssub(2) = 0.78/sqrt(act(1))
      cgsub(2) = 0.78/sqrt(act(6))
      crevp(2) = 0.78/sqrt(act(2))
      cssub(3) = 0.31*scm3*gam263*sqrt(clin/visk)/act(1)**0.65625
      cgsub(3) = 0.31*scm3*gam275*sqrt(gcon/visk)/act(6)**0.6875
      crevp(3) = 0.31*scm3*gam290*sqrt(alin/visk)/act(2)**0.725
      cssub(4) = tcond*rvgas
      cssub(5) = hlts**2*vdifu
      cgsub(4) = cssub(4)
      crevp(4) = cssub(4)
      cgsub(5) = cssub(5)
      crevp(5) = hltc**2*vdifu
!
      cgfr(1) = 20.e2*pisq*rnzr*rhor/act(2)**1.75
      cgfr(2) = 0.66
!
!sk ********************************************************************
!sk   smlt:  five constants ( lin et al. 1983 )
      csmlt(1) = 2.*pie*tcond*rnzs/hltf
      csmlt(2) = 2.*pie*vdifu*rnzs*hltc/hltf
      csmlt(3) = cssub(2)
      csmlt(4) = cssub(3)
      csmlt(5) = ch2o/hltf
!sk ********************************************************************
!     gmlt:  five constants
      cgmlt(1) = 2.*pie*tcond*rnzg/hltf
      cgmlt(2) = 2.*pie*vdifu*rnzg*hltc/hltf
      cgmlt(3) = cgsub(2)
      cgmlt(4) = cgsub(3)
      cgmlt(5) = ch2o/hltf
!sk ********************************************************************
      es0 = 6.107799961e2   ! ~6.1 mb
      ces0 = eps*es0
!
!     c2brg has conversion factor of 10**3
      c1brg = dts/rmi50
!lin  c2brg = ri50**2*1.e3 ! error
      c2brg = pie*ri50**2*1.e3
 
      do_setup = .false.

 end subroutine setupm


 subroutine micro_phys_init(axes, time)
    integer,         intent(in) :: axes(4)
    type(time_type), intent(in) :: time
    
    integer   :: unit, io, ierr
    logical   :: flag

    master = gid == 0

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
    if (gid==masterproc) write( stdlog(), nml = mp_lin_nml )
 
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
!   if ( master ) write(*,*) 'prec_lin diagnostics initialized.', id_prec

    id_cond = register_diag_field ( mod_name, 'cond_lin', axes(1:2), time,     &
         'total condensate', 'kg/m**2', missing_value=missing_value )

    call qsmith_init

    if ( master ) write(*,*) 'mp_lin diagnostics initialized.'

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
   deallocate ( table3 )
   deallocate ( tablew )
   deallocate ( des )
   deallocate ( des2 )
   deallocate ( des3 )
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



 real function acr3d(v1, v2, q1, q2, c, cac, rho)
 real, intent(in) :: v1, v2, c, rho
 real, intent(in) :: q1, q2    ! mixing ratio!!!
 real, intent(in) :: cac(3)
 real :: t1, s1, s2
!integer :: k
! real:: a
!     a=0.0
!     do k=1,3
!        a = a + cac(k)*( (q1*rho)**((7-k)*0.25) * (q2*rho)**(k*0.25) )
!     enddo
!     acr3d = c * abs(v1-v2) * a/rho
!----------
! Optimized
!----------
      t1 = sqrt(q1*rho)
      s1 = sqrt(q2*rho)
      s2 = sqrt(s1)       ! s1 = s2**2
      acr3d = c*abs(v1-v2)*q1*s2*(cac(1)*t1 + cac(2)*sqrt(t1)*s2 + cac(3)*s1)

 end function acr3d




 real function smlt(tc, dqs, qsrho,psacw,psacr,c,rho, rhofac)
 real, intent(in):: tc,dqs,qsrho,psacw,psacr,c(5),rho, rhofac
     
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
       allocate ( table3(length) )
       allocate ( tablew(length) )
       allocate (   des (length) )
       allocate (   des2(length) )
       allocate (   des3(length) )
       allocate (   desw(length) )

       call qs_table (length )
       call qs_table2(length )
       call qs_table3(length )
       call qs_tablew(length )

       do i=1,length-1
           des(i) = max(0.,  table(i+1) -  table(i))
          des2(i) = max(0., table2(i+1) - table2(i))
          des3(i) = max(0., table3(i+1) - table3(i))
          desw(i) = max(0., tablew(i+1) - tablew(i))
       enddo
        des(length) =  des(length-1)
       des2(length) = des2(length-1)
       des3(length) = des3(length-1)
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
      qs1d = eps*es/pa
        it = ap1 - 0.5
      dqdt = eps10*(des2(it) + (ap1-it)*(des2(it+1)-des2(it)))/pa

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
      ws1d = eps*es/pa
        it = ap1 - 0.5
      dqdt = eps10*(desw(it) + (ap1-it)*(desw(it+1)-desw(it)))/pa

 end function ws1d


 real function wqsat(ta, pa)
! Pure water phase
  real, intent(in):: ta, pa
! local:
  real es, ap1
  real, parameter:: tmin=tice - 160.
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = tablew(it) + (ap1-it)*desw(it)
     wqsat = eps*es/pa

 end function wqsat

 real function iqsat(ta, pa)
! Pure water phase
  real, intent(in):: ta, pa
! local:
  real es, ap1
  real, parameter:: tmin=tice - 160.
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = table2(it) + (ap1-it)*des2(it)
     iqsat = eps*es/pa

 end function iqsat



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


 subroutine esw_table1d(ta, es, n)
  integer, intent(in):: n
! For waterphase only
  real, intent(in)::  ta(n)
  real, intent(out):: es(n)
  real, parameter:: tmin=tice - 160.
  real  ap1
  integer i, it

  do i=1, n
       ap1 = 10.*dim(ta(i), tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
     es(i) = tablew(it) + (ap1-it)*desw(it)
  enddo
 end subroutine esw_table1d



 subroutine es2_table1d(ta, es, n)
  integer, intent(in):: n
! two-phase table with -2C as the transition point for ice-water phase
! For sea ice model
  real, intent(in)::  ta(n)
  real, intent(out):: es(n)
  real, parameter:: tmin=tice - 160.
  real  ap1
  integer i, it

  do i=1, n
       ap1 = 10.*dim(ta(i), tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
     es(i) = table2(it) + (ap1-it)*des2(it)
  enddo
 end subroutine es2_table1d


 subroutine es3_table1d(ta, es, n)
  integer, intent(in):: n
! two-phase table with -2C as the transition point for ice-water phase
  real, intent(in)::  ta(n)
  real, intent(out):: es(n)
  real, parameter:: tmin=tice - 160.
  real  ap1
  integer i, it

  do i=1, n
       ap1 = 10.*dim(ta(i), tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
     es(i) = table3(it) + (ap1-it)*des3(it)
  enddo
 end subroutine es3_table1d



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



 subroutine qs_table3(n)
! 2-phase table with "-2 C" as the transition point
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
!       if ( i<= 1600 ) then
        if ( i<= 1580 ) then  ! to -2 C
!  compute es over ice between -160c and 0 c.
!  see smithsonian meteorological tables page 350.
              aa  = -9.09718 *(tbasi/tem-1.)
              b   = -3.56654 *alog10(tbasi/tem)
              c   =  0.876793*(1.-tem/tbasi)
              e   = alog10(esbasi)
             table3(i) = 0.1 * 10**(aa+b+c+e)
        else
!  compute es over water between -2c and 102c.
!  see smithsonian meteorological tables page 350.
             aa  = -7.90298*(tbasw/tem-1.)
             b   =  5.02808*alog10(tbasw/tem)
             c   = -1.3816e-07*(10**((1.-tem/tbasw)*11.344)-1.)
             d   =  8.1328e-03*(10**((tbasw/tem-1.)*(-3.49149))-1.)
             e   = alog10(esbasw)
             table3(i) = 0.1 * 10**(aa+b+c+d+e)
        endif
     enddo

!----------
! smoother
!----------
      i0 = 1580
      tem0 = 0.25*(table3(i0-1) + 2.*table(i0) + table3(i0+1))
      i1 = 1581
      tem1 = 0.25*(table3(i1-1) + 2.*table(i1) + table3(i1+1))
      table3(i0) = tem0
      table3(i1) = tem1

 end subroutine qs_table3


 real function qs1d_blend(t, p, q)
! Note: this routine is based on "moist" mixing ratio
! Blended mixed phase table
  real, intent(in):: t, p, q
  real es, ap1
  real, parameter:: tmin=tice - 160.
  integer it

       ap1 = 10.*dim(t, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = table(it) + (ap1-it)*des(it)
      qs1d_blend = eps*es*(1.+zvir*q)/p

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


 subroutine qsmith(im, km, ks, t, p, q, qs, dqdt)
! input t in deg k; p (pa) : moist pressure
  integer, intent(in):: im, km, ks
  real, intent(in),dimension(im,km):: t, p, q
  real, intent(out),dimension(im,km):: qs
  real, intent(out), optional:: dqdt(im,km)
! local:
  real, parameter:: eps10 = 10.*eps
  real es(im,km)
  real ap1
  real tmin
  integer i, k, it

  tmin = tice-160.

  if( .not. allocated(table) ) then
       call  qsmith_init
  endif
 
      do k=ks,km
         do i=1,im
            ap1 = 10.*dim(t(i,k), tmin) + 1.
            ap1 = min(2621., ap1)
            it = ap1
            es(i,k) = table(it) + (ap1-it)*des(it)
            qs(i,k) = eps*es(i,k)*(1.+zvir*q(i,k))/p(i,k)
         enddo
      enddo

      if ( present(dqdt) ) then
      do k=ks,km
           do i=1,im
              ap1 = 10.*dim(t(i,k), tmin) + 1.
              ap1 = min(2621., ap1) - 0.5
              it  = ap1
              dqdt(i,k) = eps10*(des(it)+(ap1-it)*(des(it+1)-des(it)))*(1.+zvir*q(i,k))/p(i,k)
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
!      tmp1 = cp - rdgas*ptop/p1(k)
!   lcpk(k) = latv / tmp1
!   icpk(k) = lati / tmp1
    lcpk(k) = latv / cp
    icpk(k) = lati / cp
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




 subroutine sg_conv(is, ie, js, je, isd, ied, jsd, jed,    &
                    isc, iec, jsc, jec,  km, nq, dt, tau,  &
                    delp, phalf, pm, zfull, zhalf, ta, qa, ua, va, w, &
                    u_dt, v_dt, t_dt, q_dt, mcond, land, pblht, nqv, nql, nqi)
! Non-precipitating sub-grid scale convective adjustment-mixing
!-------------------------------------------
      integer, intent(in):: is, ie, js, je, km, nq
      integer, intent(in):: mcond
      integer, intent(in):: isc, iec, jsc, jec
      integer, intent(in):: isd, ied, jsd, jed
      integer, intent(in):: tau            ! Relaxation time scale
      integer, intent(in):: nqv, nql, nqi  ! vapor, liquid, ice
      real, intent(in):: dt             ! model time step
      real, intent(in):: phalf(is:ie,js:je,km+1) 
      real, intent(in):: pm(is:ie,js:je,km)
      real, intent(in):: zfull(is:ie,js:je,km)
      real, intent(in):: zhalf(is:ie,js:je,km+1)
      real, intent(in):: delp(isd:ied,jsd:jed,km)      ! Delta p at each model level
      real, intent(in)::   ta(isd:ied,jsd:jed,km)      ! Temperature
      real, intent(in)::   qa(isd:ied,jsd:jed,km,nq)   ! Specific humidity & tracers
      real, intent(in)::   ua(isd:ied,jsd:jed,km)
      real, intent(in)::   va(isd:ied,jsd:jed,km)
      real, intent(inout):: w(isd:ied,jsd:jed,km)
      real, intent(in)::  land(is:ie,js:je)
      real, intent(in):: pblht(is:ie,js:je)     ! depth (m) of the PBL
! Output:
! Updated fields:
      real, intent(out):: u_dt(isd:ied,jsd:jed,km)   ! updated u-wind field
      real, intent(out):: v_dt(isd:ied,jsd:jed,km)   !         v-wind
      real, intent(out):: t_dt(isc:iec,jsc:jec,km)   !         temperature
      real, intent(out):: q_dt(isc:iec,jsc:jec,km,nq) !
!---------------------------Local variables-----------------------------
      real, dimension(is:ie,km):: tvm, u0, v0, w0, t0, gz, hd, pkz
      real, dimension(is:ie,km+1):: pk, peln
      real q0(is:ie,km,nq) 
      real gzh(is:ie)
      real pbot, ri, pt1, pt2, lf, ratio
      real rdt, dh, dh0, dhs, dq, tv, h0, mc, mx,  fra, rk, rz, rcp
      real qs1, detn
      real clouds, rqi
      integer kcond
      integer i, j, k, n, m, iq, kk, ik
      real, parameter:: ustar2 = 1.E-8
      real, parameter:: dh_min = 1.E-4

      if ( nqv /= 1 .or. nql/=2 ) then
           call error_mesg ('sg_conv', 'Tracer indexing error', FATAL) 
      endif

!     call prt_maxmin('PBL_HT', pblht, is, ie, js, je, 0, 1, 1., master)


    rz = rvgas - rdgas          ! rz = zvir * rdgas
    rk = cp_air/rdgas + 1.
   rcp = 1./cp_air

    m = 4
    rdt = 1. / dt
    fra = dt/real(tau)

!------------
! Compute gz: center 
!------------
  do 1000 j=js,je       ! this main loop can be OpneMPed in j

    do k=mcond,km+1
       do i=is,ie
          peln(i,k) = log(phalf(i,j,k))
!           pk(i,k) = phalf(i,j,k)**kappa
            pk(i,k) = exp(kappa*peln(i,k))
       enddo
    enddo

    do k=mcond,km
       do i=is,ie
          u0(i,k) = ua(i,j,k)
          v0(i,k) = va(i,j,k)
          t0(i,k) = ta(i,j,k)
         pkz(i,k) = (pk(i,k+1)-pk(i,k))/(kappa*(peln(i,k+1)-peln(i,k)))
       enddo
    enddo

    if ( .not.hydrostatic ) then
       do k=mcond,km
          do i=is,ie
             w0(i,k) = w(i,j,k)
          enddo
       enddo
    endif

    do iq=1,nq
       do k=mcond,km
          do i=is,ie
             q0(i,k,iq) = qa(i,j,k,iq)
          enddo
       enddo
    enddo


!-----------------
! K-H instability:
!-----------------
   kcond = mcond

   do n=1,m
      ratio = real(n)/real(m)

    if( phys_hydrostatic ) then
       do i=is,ie
          gzh(i) = 0.
       enddo
       do k=km, mcond,-1
          do i=is,ie
           tvm(i,k) = t0(i,k)*(1.+zvir*q0(i,k,nqv))
                tv  = rdgas*tvm(i,k)
            gz(i,k) = gzh(i) + tv*(1.-phalf(i,j,k)/pm(i,j,k))
            hd(i,k) = cp_air*tvm(i,k)+gz(i,k)+0.5*(u0(i,k)**2+v0(i,k)**2)
             gzh(i) = gzh(i) + tv*(peln(i,k+1)-peln(i,k))
          enddo
       enddo
       do i=is,ie
          gzh(i) = 0.
       enddo
    else
       do k=mcond,km
          do i=is,ie
             gz(i,k) = grav*zfull(i,j,k)
             hd(i,k) = cp_air*t0(i,k)+gz(i,k)+0.5*(u0(i,k)**2+v0(i,k)**2+w0(i,k)**2)
          enddo
       enddo
    endif

      do k=km,kcond+1,-1
         do i=is,ie
! Richardson number at interface: g*delz * (del_theta/theta) / (del_u**2 + del_v**2)
            pt1 = t0(i,k-1)/pkz(i,k-1)
            pt2 = t0(i,k  )/pkz(i,k  )
             ri = (gz(i,k-1)-gz(i,k))*(pt1-pt2)/( 0.5*(pt1+pt2)*        &
                 ((u0(i,k-1)-u0(i,k))**2+(v0(i,k-1)-v0(i,k))**2+ustar2) )
! Dry convective mixing for K-H instability & CAT (Clear Air Turbulence):
! Compute equivalent mass flux: mc
            if ( ri < 0.25 ) then
                 mc = ratio * (1.-max(0.0,4.*ri)) ** 2
                 mc = mc*delp(i,j,k-1)*delp(i,j,k)/(delp(i,j,k-1)+delp(i,j,k))
                 do iq=1,nq
                              h0 = mc*(q0(i,k,iq)-q0(i,k-1,iq))
                    q0(i,k-1,iq) = q0(i,k-1,iq) + h0/delp(i,j,k-1)
                    q0(i,k  ,iq) = q0(i,k  ,iq) - h0/delp(i,j,k  )
                 enddo
! u:
                        h0 = mc*(u0(i,k)-u0(i,k-1))
                 u0(i,k-1) = u0(i,k-1) + h0/delp(i,j,k-1)
                 u0(i,k  ) = u0(i,k  ) - h0/delp(i,j,k  )
! v:
                        h0 = mc*(v0(i,k)-v0(i,k-1))
                 v0(i,k-1) = v0(i,k-1) + h0/delp(i,j,k-1)
                 v0(i,k  ) = v0(i,k  ) - h0/delp(i,j,k  )
! h:
                          h0 = mc*(hd(i,k)-hd(i,k-1))
                   hd(i,k-1) = hd(i,k-1) + h0/delp(i,j,k-1)
                   hd(i,k  ) = hd(i,k  ) - h0/delp(i,j,k  )
                if ( .not.hydrostatic ) then
                           h0 = mc*(w0(i,k)-w0(i,k-1))
                    w0(i,k-1) = w0(i,k-1) + h0/delp(i,j,k-1)
                    w0(i,k  ) = w0(i,k  ) - h0/delp(i,j,k  )
                endif
            endif
         enddo
!--------------
! Retrive Temp:
!--------------
      if ( phys_hydrostatic ) then
         kk = k
         do i=is,ie
            t0(i,kk) = (hd(i,kk)-gzh(i)-0.5*(u0(i,kk)**2+v0(i,kk)**2))  &
                     / ( rk - phalf(i,j,kk)/pm(i,j,kk) )
              gzh(i) = gzh(i) + t0(i,kk)*(peln(i,kk+1)-peln(i,kk))
            t0(i,kk) = t0(i,kk) / ( rdgas + rz*q0(i,kk,nqv) )
         enddo
         kk = k-1
         do i=is,ie
            t0(i,kk) = (hd(i,kk)-gzh(i)-0.5*(u0(i,kk)**2+v0(i,kk)**2))  &
                     / ((rk-phalf(i,j,kk)/pm(i,j,kk))*(rdgas+rz*q0(i,kk,nqv)))
         enddo
      else
! Non-hydrostatic under constant volume heating/cooling
         do kk=k-1,k
            do i=is,ie
               t0(i,kk) = rcp*(hd(i,kk)-gz(i,kk)-0.5*(u0(i,kk)**2+v0(i,kk)**2+w0(i,kk)**2))
            enddo
         enddo
      endif
      enddo
   enddo       ! n-loop


!-------------------------
! Moist adjustment/mixing:
!-------------------------
 m = 4

 if( km>k_moist+1 ) then
   do n=1,m

    ratio = real(n)/real(m)

    if ( phys_hydrostatic ) then
       do i=is,ie
          gzh(i) = 0.
       enddo
    endif

    do k=km,max(kcond,k_moist)+1,-1
       do i=is,ie
          if ( phalf(i,j,k) > p_crt ) then
!--------------------------------------------------------------------
!           qs1 = qs1d_blend(t0(i,k-1), pm(i,j,k-1), q0(i,k-1,nqv))
!            lf = hlv + hlf*min(1.0, max(0.0, (tice-t0(i,k-1))/30.))
!--------------------------------------------------------------------
            clouds = q0(i,k-1,nql) + q0(i,k-1,nqi)
            if( clouds > 1.e-5 ) then
                rqi = q0(i,k-1,nqi) / clouds
            else
                rqi = max(0., min(1., (tice-t0(i,k-1))/30.))
            end if
            qs1 = rqi*es2_table(t0(i,k-1)) + (1.-rqi)*esw_table(t0(i,k-1))
            qs1 = eps*qs1*(1.+zvir*q0(i,k-1,nqv))/pm(i,j,k-1)
             lf = hlv + rqi*hlf

              dh0 = hd(i,k) - hd(i,k-1)
              dhs = dh0 + lf*(q0(i,k,nqv)-qs1        )
              dh  = dh0 + lf*(q0(i,k,nqv)-q0(i,k-1,nqv))
!             if ( dhs>0.0 .and. dh>dh_min ) then
              if ( dhs>0.0 .and. dh>dh_min .and. q0(i,k,nqv)>q0(i,k-1,nqv) ) then
                   mc = ratio*min(1.0, 0.5*dhs/dh)*    &
                        delp(i,j,k-1)*delp(i,j,k)/(delp(i,j,k-1)+delp(i,j,k))
                          h0 = mc*dh0
                   hd(i,k-1) = hd(i,k-1) + h0/delp(i,j,k-1)
                   hd(i,k  ) = hd(i,k  ) - h0/delp(i,j,k  )
! Perform local mixing of all advected tracers:
#ifdef DET_CON
                 if ( zhalf(i,j,k) > (pblht(i,j)+zhalf(i,j,km+1)) ) then
                      detn = min(1., zhalf(i,j,k)/7.e3)
! specific humidity:
                              h0 = mc*(q0(i,k,nqv)-q0(i,k-1,nqv))
                              dq = h0/delp(i,j,k-1)
                   q0(i,k-1,nqv) = q0(i,k-1,nqv) + dq*(1.-detn)
                   q0(i,k  ,nqv) = q0(i,k  ,nqv) - h0/delp(i,j,k  )
                   do iq=2,nq
                                h0 = mc*(q0(i,k,iq)-q0(i,k-1,iq))
                      q0(i,k-1,iq) = q0(i,k-1,iq) + h0/delp(i,j,k-1)
                      q0(i,k  ,iq) = q0(i,k  ,iq) - h0/delp(i,j,k  )
                   enddo
!--------------
! Condensation:
!--------------
                   dq = dq * detn
                   q0(i,k-1,nql) = q0(i,k-1,nql) + dq*(1.-rqi)
                   q0(i,k-1,nqi) = q0(i,k-1,nqi) + dq*rqi
                   hd(i,k-1) = hd(i,k-1) + dq*lf

                 else
                   do iq=1,nq
                                h0 = mc*(q0(i,k,iq)-q0(i,k-1,iq))
                      q0(i,k-1,iq) = q0(i,k-1,iq) + h0/delp(i,j,k-1)
                      q0(i,k  ,iq) = q0(i,k  ,iq) - h0/delp(i,j,k  )
                   enddo
                 endif
#else
                   do iq=1,nq
                                h0 = mc*(q0(i,k,iq)-q0(i,k-1,iq))
                      q0(i,k-1,iq) = q0(i,k-1,iq) + h0/delp(i,j,k-1)
                      q0(i,k  ,iq) = q0(i,k  ,iq) - h0/delp(i,j,k  )
                   enddo
#endif
! u:
                          h0 = mc*(u0(i,k)-u0(i,k-1))
                   u0(i,k-1) = u0(i,k-1) + h0/delp(i,j,k-1)
                   u0(i,k  ) = u0(i,k  ) - h0/delp(i,j,k  )
! v:
                          h0 = mc*(v0(i,k)-v0(i,k-1))
                   v0(i,k-1) = v0(i,k-1) + h0/delp(i,j,k-1)
                   v0(i,k  ) = v0(i,k  ) - h0/delp(i,j,k  )
! *** Non-hydrostatic:
                  if ( .not.hydrostatic ) then
                          h0 = mc*(w0(i,k)-w0(i,k-1))
                   w0(i,k-1) = w0(i,k-1) + h0/delp(i,j,k-1)
                   w0(i,k  ) = w0(i,k  ) - h0/delp(i,j,k  )
                  endif
! ***
              endif  ! dh check
            endif    ! p_crt check
         enddo
!--------------
! Retrive Temp:
!--------------
       if ( phys_hydrostatic ) then
         kk = k
         do i=is,ie
            t0(i,kk) = (hd(i,kk)-gzh(i)-0.5*(u0(i,kk)**2+v0(i,kk)**2))  &
                     / ( rk - phalf(i,j,kk)/pm(i,j,kk) )
              gzh(i) = gzh(i) + t0(i,kk)*(peln(i,kk+1)-peln(i,kk))
            t0(i,kk) = t0(i,kk) / ( rdgas + rz*q0(i,kk,nqv) )
         enddo
         kk = k-1
         do i=is,ie
            t0(i,kk) = (hd(i,kk)-gzh(i)-0.5*(u0(i,kk)**2+v0(i,kk)**2))  &
                     / ((rk-phalf(i,j,kk)/pm(i,j,kk))*(rdgas+rz*q0(i,kk,nqv)))
         enddo
       else
! Non-hydrostatic under constant volume heating/cooling
         do kk=k-1,k
            do i=is,ie
               t0(i,kk) = rcp*(hd(i,kk)-gz(i,kk)-0.5*(u0(i,kk)**2+v0(i,kk)**2+w0(i,kk)**2))
            enddo
         enddo
       endif
      enddo
   enddo       ! n-loop
 endif      ! k_moist check

   if ( fra < 1. ) then
      do k=mcond,km
         do i=is,ie
            t0(i,k) = ta(i,j,k) + (t0(i,k) - ta(i,j,k))*fra
            u0(i,k) = ua(i,j,k) + (u0(i,k) - ua(i,j,k))*fra
            v0(i,k) = va(i,j,k) + (v0(i,k) - va(i,j,k))*fra
         enddo
      enddo

      if ( .not.hydrostatic ) then
      do k=mcond,km
         do i=is,ie
            w0(i,k) = w(i,j,k) + (w0(i,k) - w(i,j,k))*fra
         enddo
      enddo
      endif

      do iq=1,nq
         do k=mcond,km
            do i=is,ie
               q0(i,k,iq) = qa(i,j,k,iq) + (q0(i,k,iq) - qa(i,j,k,iq))*fra
            enddo
         enddo
      enddo
   endif

!--------------------
! Update fields:
!--------------------
   do k=1,mcond-1
      do i=is,ie
         u_dt(i,j,k) = ua(i,j,k)
         v_dt(i,j,k) = va(i,j,k)
         t_dt(i,j,k) = ta(i,j,k)
      enddo
   enddo
   do k=mcond,km
      do i=is,ie
         u_dt(i,j,k) = u0(i,k)
         v_dt(i,j,k) = v0(i,k)
         t_dt(i,j,k) = t0(i,k)
      enddo
   enddo

   if ( .not.hydrostatic ) then
      do k=mcond,km
         do i=is,ie
            w(i,j,k) = w0(i,k)
         enddo
      enddo
   endif

   do iq=1,nq
      do k=1,mcond-1
         do i=is,ie
            q_dt(i,j,k,iq) = qa(i,j,k,iq)
         enddo
      enddo
      do k=mcond,km
         do i=is,ie
            q_dt(i,j,k,iq) = q0(i,k,iq)
         enddo
      enddo
   enddo

1000 continue


 end subroutine sg_conv

end module mp_lin_mod
