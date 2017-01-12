!
! Cloud micro-physics package for GFDL global cloud resolving model
! The algorithms are originally derived from Lin et al 1983. Most of the key
! elements have been simplified/improved. This code at this stage bears little
! to no similarity to the original Lin MP in Zeta. Therefore, it is best to be called
! GFDL Micro-Physics (GFDL MP).
! Developer: Shian-Jiann Lin
!
module lin_cld_microphys_mod
! use mpp_mod,           only: stdlog, mpp_pe, mpp_root_pe, mpp_clock_id, &
!                              mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE, &
!                              input_nml_file, mpp_max
! use diag_manager_mod,  only: register_diag_field, send_data
! use time_manager_mod,  only: time_type, get_time
! use constants_mod,     only: grav, rdgas, rvgas, cp_air, cp_vapor, hlv, hlf, kappa, pi=>pi_8
! use fms_mod,           only: write_version_number, open_namelist_file, &
!                              check_nml_error, file_exist, close_file,  &
!                              error_mesg, FATAL

 implicit none
 private

 public  lin_cld_microphys_driver, lin_cld_microphys_init, lin_cld_microphys_end, wqs1, wqs2, qs_blend
 public  qsmith_init, qsmith, es2_table1d, es3_table1d, esw_table1d, wqsat_moist, wqsat2_moist
 public  setup_con, wet_bulb
 public  cloud_diagnosis
 real             :: missing_value = -1.e10
 logical          :: module_is_initialized = .false.
 logical          :: qsmith_tables_initialized = .false.
 character(len=17) :: mod_name = 'lin_cld_microphys'

!==== constants_mod ====
integer, public, parameter :: R_GRID=8
real, parameter :: grav = 9.80665_R_GRID
real, parameter :: rdgas = 287.05_R_GRID
real, parameter :: rvgas = 461.50_R_GRID
real, parameter :: cp_air = 1004.6_R_GRID
real, parameter :: cp_vapor = 4.0_R_GRID*RVGAS
real, parameter :: hlv = 2.5e6_R_GRID
real, parameter :: hlf = 3.3358e5_R_GRID
real, parameter :: kappa = rdgas/cp_air
real, parameter :: pi = 3.1415926535897931_R_GRID
!==== constants_mod ====

!==== fms constants ====================
!!! real, parameter :: latv  = hlv             ! = 2.500e6
!!!  real, parameter:: cv_air = 717.56   ! Satoh value
!!! real, parameter :: lati  = hlf             ! = 3.34e5
!!! real, parameter :: lats  = latv+lati       ! = 2.834E6
! rdgas = 287.04;   rvgas = 461.50
! cp_air =rdgas * 7./2. = 1006.64           ! heat capacity at constant pressure (j/kg/k)
! The following two are from Emanuel's book "Atmospheric Convection"
!!! real, parameter :: c_liq = 4190.           ! heat capacity of water at 0C
                                            !
 real, parameter :: eps   = rdgas/rvgas     ! = 0.621971831
 real, parameter :: zvir  = rvgas/rdgas-1.  ! = 0.607789855
 real, parameter :: table_ice  = 273.16  ! freezing point for qs table
 real, parameter :: cv_air =  cp_air - rdgas ! = rdgas * (7/2-1) = 2.5*rdgas=717.68

   real, parameter:: e00 = 610.71  ! saturation vapor pressure at T0
!   real, parameter:: c_liq = 4190.       ! heat capacity of water at 0C
! GFS value
   real, parameter:: c_liq = 4.1855e+3       ! heat capacity of water at 0C
   real, parameter:: c_ice = 2106.           ! heat capacity of ice at 0C: c=c_ice+7.3*(T-Tice) 
   real, parameter:: cp_vap = cp_vapor   ! 1846.
!  real, parameter:: cv_vap = 1410.0     ! Emanuel value
! For consistency, cv_vap derived FMS constants:
   real, parameter:: cv_vap = cp_vap - rvgas  ! 1384.5
   real, parameter:: dc_vap = cp_vap - c_liq     ! = -2344.    isobaric heating/cooling
   real, parameter:: dc_ice =  c_liq - c_ice     ! =  2084

! Values at 0 Deg C
! GFS value
   real, parameter:: hlv0 = 2.5e6
!  real, parameter:: hlv0 = 2.501e6   ! Emanuel Appendix-2
! GFS value
   real, parameter:: hlf0 = 3.3358e5
!  real, parameter:: hlf0 = 3.337e5   ! Emanuel
   real, parameter:: t_ice = 273.16
! Latent heat at absolute zero:
   real, parameter:: Lv0 =  hlv0 - dc_vap*t_ice   ! = 3.141264e6
   real, parameter:: li00 = hlf0 - dc_ice*t_ice   ! = -2.355446e5

   real, parameter:: d2ice = cp_vap - c_ice
   real, parameter:: Li2 = hlv0+hlf0 - d2ice*t_ice

 real, parameter :: qrmin  = 1.e-8
 real, parameter :: qvmin  = 1.e-20      ! min value for water vapor (treated as zero)
 real, parameter :: qcmin  = 1.e-12      ! min value for cloud condensates
 real, parameter :: sfcrho = 1.2         ! surface air density
 real, parameter :: vr_min = 1.e-3       ! minimum fall speed for rain/graupel
 real, parameter :: vf_min = 1.0E-5
 real, parameter :: rhor   = 1.0e3  ! LFO83
 real, parameter :: dz_min = 1.e-2
 real :: cracs, csacr, cgacr, cgacs, acco(3,4), csacw,          &
         craci, csaci, cgacw, cgaci, cracw, cssub(5), cgsub(5), &
         crevp(5), cgfr(2), csmlt(5), cgmlt(5)
 real :: es0, ces0
 real :: pie, rgrav, fac_rc
 real :: lcp, icp, tcp
 real :: lv00, c_air, c_vap

 logical :: de_ice = .false.     !
 logical :: sedi_transport = .true.     !
 logical :: do_sedi_w = .false.
 logical :: do_sedi_heat = .true.      !
 logical :: prog_ccn = .false.     ! do prognostic CCN (Yi Ming's method)
 logical :: do_qa   = .true.       ! do inline cloud fraction
 logical :: rad_snow =.false.
 logical :: rad_graupel =.false.
 logical :: rad_rain =.false.
 logical :: fix_negative =.false.
 logical :: do_setup=.true.
 logical :: master
 logical :: p_nonhydro = .false.

 real, allocatable:: table(:), table2(:), table3(:), tablew(:), des(:), des2(:), des3(:), desw(:)
 logical :: tables_are_initialized = .false.

 integer:: id_rh, id_vtr, id_vts,  id_vtg, id_vti, id_rain, id_snow, id_graupel, &
           id_ice, id_prec, id_cond, id_var, id_droplets
 real:: lati, latv, lats

 real, parameter :: dt_fr = 8.       ! homogeneous freezing of all cloud water at t_wfr - dt_fr
                                     ! minimum temperature water can exist (Moore & Molinero Nov. 2011, Nature)
                                     ! dt_fr can be considered as the error bar
 integer :: lin_cld_mp_clock   ! clock for timing of driver routine

 real :: t_snow_melt = 16.      ! snow melt tempearture scale factor
 real :: t_grau_melt = 32.      ! graupel melt tempearture scale factor
 real :: p_min = 100.    ! minimum pressure (Pascal) for MP to operate

! For cloud-resolving: 1-5 km
!                       qi0_crt = 0.8E-4
!                       qs0_crt = 0.6E-3
!                       c_psaci = 0.1
!                       c_pgacs = 0.1
!----------------------
! namelist  parameters:
!----------------------
 real :: cld_min = 0.05
 real :: tice  = 273.16  ! set tice = 165. to trun off ice-phase phys (Kessler emulator)

 real :: qc_crt  = 5.0e-8  ! minimum condensate mixing ratio to allow partial cloudiness
 real :: t_min   = 178.  ! Min temp to freeze-dry all water vapor
 real :: t_sub   = 184.  ! Min temp for sublimation of cloud ice
 real :: mp_time = 150.  ! maximum micro-physics time step (sec)

 real :: rh_inc = 0.25   ! rh increment for complete evap of ql and qi
 real :: rh_inr = 0.25
 real :: rh_ins = 0.25   ! rh increment for sublimation of snow

! The following 3 time scales are for melting during terminal falls
 real :: tau_s   = 600.    ! snow melt
 real :: tau_g   = 600.    ! graupel melt
 real :: tau_mlt = 600.    ! ice melting time-scale

! Fast MP:
 real :: tau_i2s = 1000.  ! ice2snow auto-conversion time scale (sec) 
! cloud water
 real :: tau_v2l = 150.   ! vapor --> cloud water (condensation)  time scale
 real :: tau_l2v = 300.   ! cloud water --> vapor (evaporation)  time scale
! Graupel
 real :: tau_g2v = 900.   ! Grapuel sublimation time scale
 real :: tau_v2g = 21600. ! Grapuel deposition -- make it a slow process

 real :: dw_land  = 0.20  ! base value for subgrid deviation/variability over land
 real :: dw_ocean = 0.15  ! base value for ocean
 real :: ccn_o =  90.
 real :: ccn_l = 270.
 real :: rthresh = 10.0e-6     ! critical cloud drop radius (micro m)

!-------------------------------------------------------------
! WRF/WSM6 scheme: qi_gen = 4.92e-11 * (1.e3*exp(0.1*tmp))**1.33
!  optimized:      qi_gen = 4.92e-11 * exp( 1.33*log(1.e3*exp(0.1*tmp)) )
! qi_gen ~ 4.808e-7 at 0 C; 1.818e-6 at -10 C, 9.82679e-5 at -40C
! the following value is constructed such that qc_crt = 0 at zero C and @ -10C matches
! WRF/WSM6 ice initiation scheme; qi_crt = qi_gen*min(qi_lim, 0.1*tmp) / den
!
 real :: qi_gen  = 1.82E-6
 real :: qi_lim  = 10.
 real :: ql_mlt  = 2.0e-3    ! max value of cloud water allowed from melted cloud ice
 real :: ql_gen  = 1.0e-3    ! max ql generation during remapping step if fast_sat_adj = .T.
 real :: sat_adj0 = 0.99     ! adjustment factor (0: no, 1: full) during fast_sat_adj

! Cloud condensate upper bounds: "safety valves" for ql & qi
 real :: ql0_max = 2.0e-3    ! max ql value (auto converted to rain)
 real :: qi0_max = 5.0e-6    ! max qi value (by other sources)

 real :: qi0_crt = 1.0e-4    ! ice  --> snow autocon threshold (was 1.E-4)
                             ! qi0_crt is highly dependent on horizontal resolution
 real :: qr0_crt = 1.0e-4    ! rain --> snow or graupel/hail threshold
                             ! LFO used *mixing ratio* = 1.E-4 (hail in LFO)
 real :: c_paut  = 0.55     ! autoconversion ql --> qr  (use 0.5 to reduce autoconversion)
 real :: c_psaci = 0.02     ! accretion: cloud ice --> snow (was 0.1 in Zetac)
 real :: c_piacr = 5.0      ! accretion: rain --> ice:
 real :: c_cracw = 0.9      ! rain accretion efficiency

! Decreasing  clin to reduce csacw (so as to reduce cloud water ---> snow)
 real:: alin = 842.0
 real:: clin = 4.8      ! 4.8 --> 6. (to ehance ql--> qs)

!-----------------
! Graupel control:
!-----------------
 real :: qs0_crt = 1.0e-3   ! snow --> graupel density threshold (0.6e-3 in Purdue Lin scheme)
 real :: c_pgacs = 2.0e-3   ! snow --> graupel "accretion" eff. (was 0.1 in Zetac)

! fall velocity tuning constants:
 logical :: const_vi  = .false.    ! If .T. the constants are specified by v*_fac
 logical :: const_vs  = .false.
 logical :: const_vg  = .false.
 logical :: const_vr  = .false.
                     ! Good values:
 real :: vi_fac = 1.      ! If const_vi: 1/3
 real :: vs_fac = 1.      ! If const_vs: 1.
 real :: vg_fac = 1.      ! If const_vg: 2.
 real :: vr_fac = 1.      ! If const_vr: 4.
! Upper bounds of fall speed (with variable speed option)
 real :: vi_max = 0.5   !  max fall speed for ice
 real :: vs_max = 5.0   !  max fall speed for snow
 real :: vg_max = 8.0   !  max fall speed for graupel
 real :: vr_max = 12.   !  max fall speed for rain

 logical :: fast_sat_adj = .false.
 logical :: z_slope_liq  = .true.          !  use linear mono slope for autocconversions
 logical :: z_slope_ice  = .false.          !  use linear mono slope for autocconversions
 logical :: use_ccn      = .false.
 logical :: use_ppm      = .false.
 logical :: mono_prof = .true.          ! perform terminal fall with mono ppm scheme
 logical :: mp_print = .false.

 real:: global_area = -1.

 real:: tice0, t_wfr
 real:: log_10

 public mp_time, t_min, t_sub, tau_s, tau_g, dw_land, dw_ocean,  &
        vi_fac, vr_fac, vs_fac, vg_fac, ql_mlt, do_qa, fix_negative, &
        vi_max, vs_max, vg_max, vr_max,        &
        qs0_crt, qi_gen, ql0_max, qi0_max, qi0_crt, qr0_crt, fast_sat_adj, &
        rh_inc, rh_ins, rh_inr, const_vi, const_vs, const_vg, const_vr,    &
        use_ccn, rthresh, ccn_l, ccn_o, qc_crt, tau_g2v, tau_v2g, sat_adj0,    &
        c_piacr, tau_mlt, tau_v2l, tau_l2v, tau_i2s, qi_lim, ql_gen,  &
        c_paut, c_psaci, c_pgacs, z_slope_liq, z_slope_ice, prog_ccn,  &
        c_cracw, alin, clin, tice, rad_snow, rad_graupel, rad_rain,   &
        cld_min, use_ppm, mono_prof, do_sedi_heat, sedi_transport,   &
        do_sedi_w, de_ice, mp_print

!---- version number -----
 character(len=128) :: version = '$Id: lin_cloud_microphys.F90,v 21.0.2.1 2014/12/18 21:14:54 Lucas.Harris Exp $'
 character(len=128) :: tagname = '$Name:  $'

 contains


!  subroutine lin_cld_microphys_driver(qv, ql, qr, qi, qs, qg, qa, qn,                &
!                               qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, qa_dt,      &
!                               pt_dt, pt, w, uin, vin, udt, vdt, dz, delp, area, dt_in, &
!                               land,  rain, snow, ice, graupel,                      &
!                               hydrostatic, phys_hydrostatic,                        &
!                               iis,iie, jjs,jje, kks,kke, ktop, kbot, time)
  subroutine lin_cld_microphys_driver(qv, ql, qr, qi, qs, qg, qa, qn,                &
                               qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, qa_dt,      &
                               pt_dt, pt, w, uin, vin, udt, vdt, dz, delp, area, dt_in, &
                               land,  rain, snow, ice, graupel,                      &
                               hydrostatic, phys_hydrostatic,                        &
                               iis,iie, jjs,jje, kks,kke, ktop, kbot, seconds)
! kks == 1; kke == kbot == npz
  logical,         intent(in):: hydrostatic, phys_hydrostatic
  integer,         intent(in):: iis,iie, jjs,jje  ! physics window
  integer,         intent(in):: kks,kke           ! vertical dimension
  integer,         intent(in):: ktop, kbot        ! vertical compute domain
  integer,         intent(in):: seconds
  real,            intent(in):: dt_in

  real, intent(in   ), dimension(:,:)  :: area
  real, intent(in   ), dimension(:,:)  :: land  !land fraction
  real, intent(out  ), dimension(:,:)  :: rain, snow, ice, graupel
  real, intent(in   ), dimension(:,:,:):: delp, dz, uin, vin
  real, intent(in   ), dimension(:,:,:):: pt, qv, ql, qr, qg, qa, qn 
  real, intent(inout), dimension(:,:,:):: qi, qs
  real, intent(inout), dimension(:,:,:):: pt_dt,  qa_dt, udt, vdt, w
  real, intent(inout), dimension(:,:,:):: qv_dt, ql_dt, qr_dt, qi_dt,  &
                                          qs_dt, qg_dt


! local:
  logical used
  real    :: mpdt, rdt, dts, convt, tot_prec
  integer :: i,j,k
  integer :: is,ie, js,je  ! physics window
  integer :: ks,ke         ! vertical dimension
  integer :: days, ntimes

  real, dimension(iie-iis+1,jje-jjs+1)           :: prec_mp, prec1, cond, w_var, rh0
  real, dimension(iie-iis+1,jje-jjs+1,kke-kks+1) :: vt_r, vt_s, vt_g, vt_i, qn2
  real, dimension(size(pt,1), size(pt,3)) :: m2_rain, m2_sol

  real :: allmax

  is = 1
  js = 1
  ks = 1
  ie = iie-iis+1
  je = jje-jjs+1
  ke = kke-kks+1

!  call mpp_clock_begin (lin_cld_mp_clock)

  if ( phys_hydrostatic .or. hydrostatic ) then
        lv00 = Lv0
      c_air = cp_air
      c_vap = cp_vapor
      p_nonhydro = .false.
  else
        lv00 = Lv0
      c_air = cv_air
      c_vap = cv_vap
      p_nonhydro = .true.
  endif
  if ( hydrostatic ) do_sedi_w = .false.
  latv = hlv
  lati = hlf
  lats = latv + lati
  lcp = latv / cp_air
  icp = lati / cp_air
  tcp = (latv+lati) / cp_air

! tendency zero out for am moist processes should be done outside the driver
     mpdt = min(dt_in, mp_time)
      rdt = 1. / dt_in
   ntimes = nint( dt_in/mpdt )
! small time step:
      dts = dt_in / real(ntimes)

!  call get_time (time, seconds, days)


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
     call mpdrv( hydrostatic, uin, vin, w, delp, pt, qv, ql, qr, qi, qs, qg, qa, qn, dz,  &
                 is, ie, js, je, ks, ke, ktop, kbot, j, dt_in,  &
                 ntimes, rain(:,j), snow(:,j), graupel(:,j), &
                 ice(:,j), m2_rain, m2_sol, cond(:,j), area(:,j), land(:,j),  &
                 udt, vdt, pt_dt, qv_dt, ql_dt, qr_dt, qi_dt,    &
                 qs_dt, qg_dt, qa_dt, &
                 w_var, vt_r, vt_s, vt_g, vt_i, qn2 )

  enddo

! no clouds allowed above ktop
   if ( ks < ktop ) then
      do k=ks, ktop
        if ( do_qa ) then
           do j=js,je
              do i=is,ie
                 qa_dt(i,j,k) = 0.
              enddo
           enddo
        else
           do j=js,je
              do i=is,ie
!!!              qa_dt(i,j,k) = -qa(i,j,k) * rdt
                 qa_dt(i,j,k) = 0.   ! GFS
              enddo
           enddo
        endif
      enddo
   endif

   if ( id_vtr> 0 ) then
!     used=send_data(id_vtr, vt_r, time, is_in=iis, js_in=jjs)
   endif

   if ( id_vts> 0 ) then
!     used=send_data(id_vts, vt_s, time, is_in=iis, js_in=jjs)
   endif

   if ( id_vtg> 0 ) then
!     used=send_data(id_vtg, vt_g, time, is_in=iis, js_in=jjs)
   endif

   if ( id_vti> 0 ) then
!     used=send_data(id_vti, vt_i, time, is_in=iis, js_in=jjs)
   endif

   if ( id_droplets> 0 ) then
!     used=send_data(id_droplets, qn2, time, is_in=iis, js_in=jjs)
   endif

   if ( id_var> 0 ) then
!     used=send_data(id_var, w_var, time, is_in=iis, js_in=jjs)
   endif

! Convert to mm/day
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
!        used=send_data(id_cond, cond, time, is_in=iis, js_in=jjs)
   endif

   if ( id_snow>0 ) then
!       used=send_data(id_snow,    snow,    time, iis, jjs)
!        used=send_data(id_snow, snow, time, is_in=iis, js_in=jjs)
        if ( mp_print .and. seconds==0 ) then
!             tot_prec = g_sum(snow, is, ie, js, je, area, 1)
!             if(master) write(*,*) 'mean snow=', tot_prec
        endif
   endif

   if ( id_graupel>0 ) then
!       used=send_data(id_graupel, graupel, time, iis, jjs)
!        used=send_data(id_graupel, graupel, time, is_in=iis, js_in=jjs)
        if ( mp_print .and. seconds==0 ) then
!             tot_prec = g_sum(graupel, is, ie, js, je, area, 1)
!             if(master) write(*,*) 'mean graupel=', tot_prec
        endif
   endif

   if ( id_ice>0 ) then
!       used=send_data(id_ice, ice, time, iis, jjs)
!        used=send_data(id_ice, ice, time, is_in=iis, js_in=jjs)
        if ( mp_print .and. seconds==0 ) then
!             tot_prec = g_sum(ice, is, ie, js, je, area, 1)
!             if(master) write(*,*) 'mean ice_mp=', tot_prec
        endif
   endif

   if ( id_rain>0 ) then
!       used=send_data(id_rain,    rain,    time, iis, jjs)
!        used=send_data(id_rain, rain, time, is_in=iis, js_in=jjs)
        if ( mp_print .and. seconds==0 ) then
!             tot_prec = g_sum(rain, is, ie, js, je, area, 1)
!             if(master) write(*,*) 'mean rain=', tot_prec
        endif
   endif

   if ( id_rh>0 ) then !Not used?
!       used=send_data(id_rh,  rh0,   time, iis, jjs)
!        used=send_data(id_rh, rh0, time, is_in=iis, js_in=jjs)
   endif


   if ( id_prec>0 ) then
!       used=send_data(id_prec, prec_mp, time, iis, jjs)
!        used=send_data(id_prec, prec_mp, time, is_in=iis, js_in=jjs)
   endif

!----------------------------------------------------------------------------

   if ( mp_print ) then
        prec1(:,:) = prec1(:,:) + prec_mp(:,:)
        if ( seconds==0 ) then
             prec1(:,:) = prec1(:,:)*dt_in/86400.
!             tot_prec = g_sum(prec1, is, ie, js, je, area, 1)
!             if(master) write(*,*) 'Daily prec_mp=', tot_prec
             prec1(:,:) = 0.
        endif
   endif
!----------------------------------------------------------------------------


!   call mpp_clock_end (lin_cld_mp_clock)

 end subroutine lin_cld_microphys_driver



 subroutine mpdrv( hydrostatic, uin, vin, w, delp, pt, qv, ql, qr, qi, qs, qg, qa, qn, dz,     &
                   is, ie, js, je, ks, ke, ktop, kbot, j, dt_in, ntimes,  &
                   rain, snow, graupel, ice, m2_rain, m2_sol, &
                   cond, area1, land, u_dt, v_dt, pt_dt, qv_dt, ql_dt, qr_dt, qi_dt, &
                   qs_dt, qg_dt, qa_dt, &
                   w_var, vt_r, vt_s, vt_g, vt_i, qn2 )

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

  logical,         intent(in):: hydrostatic
  integer,         intent(in):: j, is,ie, js,je, ks,ke
  integer,         intent(in):: ntimes, ktop, kbot
  real,            intent(in):: dt_in

  real, intent(in), dimension(is:):: area1, land
  real, intent(in), dimension(is:,js:,ks:):: uin, vin, delp, pt, qv, ql, qr, qg, qa, qn, dz
  real, intent(inout), dimension(is:,js:,ks:):: qi, qs
  real, intent(inout), dimension(is:,js:,ks:):: u_dt, v_dt, w, pt_dt,  qa_dt,  &
                                            qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt
  real, intent(inout), dimension(is:):: rain, snow, ice, graupel, cond
  real, intent(out), dimension(is:,js:)       :: w_var
  real, intent(out), dimension(is:,js:,ks:) :: vt_r, vt_s, vt_g, vt_i, qn2
  real, intent(out), dimension(is:,ks:):: m2_rain,  m2_sol

!----------
! local var
!----------
  real, dimension(ktop:kbot):: qvz, qlz, qrz, qiz, qsz, qgz, qaz, &
                               vtiz, vtsz, vtgz, vtrz, dp0, &
                               dp1, qv0, ql0, qr0, qi0, qs0, qg0, qa0, t0, den, &
                               den0, tz, p1, dz0, dz1, denfac
  real, dimension(ktop:kbot):: ccn, c_praut, m1_rain, m1_sol, m1, u0, v0, u1, v1, w1
  real :: cpaut,  rh_adj, rh_rain
  real :: r1, s1, i1, g1, rdt, ccn0
  real :: dt_rain, dts
  real :: s_leng, t_land, t_ocean, h_var
  real :: cvm, tmp, omq
  real:: dqi, qio, qin
  integer :: i,k,n

       dts = dt_in / real(ntimes)
   dt_rain = dts * 0.5
   rdt = 1. / dt_in

   cpaut = c_paut*0.104*grav/1.717e-5

   do 2000 i=is, ie
      do k=ktop, kbot
         qiz(k) = qi(i,j,k)
         qsz(k) = qs(i,j,k)
      enddo
      if ( de_ice ) then
! This is to prevent excessive build-up of cloud ice from external sources
        do k=ktop, kbot
           qio = qiz(k) - dt_in*qi_dt(i,j,k)   ! original qi before phys
           qin = max(qio, qi0_max)             ! Adjusted value
           if ( qiz(k) > qin ) then
                qsz(k) = qsz(k) + qiz(k) - qin
                qiz(k) = qin
                dqi = (qin - qio)*rdt  ! modified qi tendency
                qs_dt(i,j,k) = qs_dt(i,j,k) + qi_dt(i,j,k) - dqi
                qi_dt(i,j,k) = dqi
                qi(i,j,k) = qiz(k)
                qs(i,j,k) = qsz(k)
           endif
        enddo
      endif

   do k=ktop, kbot
       t0(k) = pt(i,j,k)
       tz(k) = t0(k)
      dp1(k) = delp(i,j,k)
      dp0(k) = dp1(k)       ! moist air mass * grav
!-- Moist mixing ratios ------------
      qvz(k) = qv(i,j,k)
      qlz(k) = ql(i,j,k)
      qrz(k) = qr(i,j,k)
      qgz(k) = qg(i,j,k)
!-- dp1: Dry air_mass
!      dp1(k) = dp1(k)*(1. -(qvz(k)+qlz(k)+qrz(k)+qiz(k)+qsz(k)+qgz(k)))
      dp1(k) = dp1(k)*(1. -qvz(k))   ! GFS
!-- Convert to dry mixing ratios
         omq = dp0(k) / dp1(k)
      qvz(k) = qvz(k)*omq
      qlz(k) = qlz(k)*omq
      qrz(k) = qrz(k)*omq
      qiz(k) = qiz(k)*omq
      qsz(k) = qsz(k)*omq
      qgz(k) = qgz(k)*omq
!-----------------------------------
      qa0(k) = qa(i,j,k)
      qaz(k) = 0.
      dz0(k) = dz(i,j,k)
!--------------------------
     den0(k) = -dp1(k)/(grav*dz0(k))  ! density of dry air
       p1(k) = den0(k)*rdgas*t0(k)    ! dry air pressure
!------------------------------
! Save a copy of old value for computing tendencies
      qv0(k) = qvz(k)
      ql0(k) = qlz(k)
      qr0(k) = qrz(k)
      qi0(k) = qiz(k)
      qs0(k) = qsz(k)
      qg0(k) = qgz(k)
! *** For sedi_momentum ***
       m1(k) = 0.
       u0(k) = uin(i,j,k)
       v0(k) = vin(i,j,k)
       u1(k) = u0(k)
       v1(k) = v0(k)
! *** For sedi_momentum ***
   enddo

   if ( do_sedi_w ) then
     do k=ktop, kbot
        w1(k) = w(i,j,k)
     enddo
   endif
! Based on Klein Eq. 15

 if ( prog_ccn ) then
   do k=ktop, kbot
!convert #/cc to #/m^3
      ccn(k) = qn(i,j,k) * 1.e6
      c_praut(k) = cpaut * (ccn(k)*rhor)**(-1./3.)
   enddo
   use_ccn = .false.
 else
   ccn0 = (ccn_l*land(i) + ccn_o*(1.-land(i))) * 1.e6
   if ( use_ccn ) then
! *** CCN is formulted as CCN = CCN_surface * (den/den_surface) ***
        ccn0 = ccn0 * rdgas*tz(kbot)/p1(kbot)
   endif
   tmp = cpaut * (ccn0*rhor)**(-1./3.)
   do k=ktop, kbot
      c_praut(k) = tmp
          ccn(k) = ccn0
   enddo
  endif

!--------------------------------------------------------
! Total water subgrid deviation in horizontal direction
!--------------------------------------------------------
!       default area dependent form: use dx ~ 100 km as the base
   s_leng  = sqrt(sqrt(area1(i)/1.E10))
   t_land  = dw_land  * s_leng
   t_ocean = dw_ocean * s_leng
   h_var = t_land*land(i) + t_ocean*(1.-land(i))
   h_var = min(0.20, max(0.01, h_var))            ! cap:
   if ( id_var>0 ) w_var(i,j) = h_var

   rh_adj  = 1. - h_var - rh_inc
   rh_rain = max(0.35, rh_adj - rh_inr)  ! rh_inr = 0.25

!-------------------------
! * fix all negatives
!-------------------------

 if ( fix_negative )  &
 call neg_adj(ktop, kbot, tz, dp1, qvz, qlz, qrz, qiz, qsz, qgz)

 m2_rain(:,:) = 0.
 m2_sol(:,:) = 0.
 do 1000 n=1,ntimes

   if ( p_nonhydro ) then
       do k=ktop, kbot
          dz1(k) = dz0(k)
          den(k) = den0(k)    ! dry air density remains the same
          denfac(k) = sqrt(sfcrho/den(k))
       enddo
   else
       do k=ktop, kbot
          dz1(k) = dz0(k)*tz(k)/t0(k)    ! hydrostatic balance
          den(k) = den0(k)*dz0(k)/dz1(k)
          denfac(k) = sqrt(sfcrho/den(k))
       enddo
   endif

!-------------------------------------------
! Time-split warm rain processes: first pass
!-------------------------------------------
   call warm_rain(dt_rain, ktop, kbot, dp1, dz1, tz, qvz, qlz, qrz, qiz, qsz, qgz, den, denfac, ccn,  &
                  c_praut, rh_rain, vtrz, r1, m1_rain, w1, h_var)
   rain(i) = rain(i) + r1

   do k=ktop, kbot
      m2_rain(i,k) = m2_rain(i,k) + m1_rain(k)
      m1(k) = m1(k) + m1_rain(k)
   enddo

!------------------------------------------------
! * sedimentation of cloud ice, snow, and graupel
!------------------------------------------------
   call fall_speed(ktop, kbot, den, qsz, qiz, qgz, qlz, tz, vtsz, vtiz, vtgz)

   call terminal_fall ( dts, ktop, kbot, tz, qvz, qlz, qrz, qgz, qsz, qiz,  &
                        dz1, dp1, den, vtgz, vtsz, vtiz,    &
                        r1, g1, s1, i1, m1_sol, w1 )
      rain(i) = rain(i)    + r1  ! from melted snow & ice that reached the ground
      snow(i) = snow(i)    + s1
   graupel(i) = graupel(i) + g1
       ice(i) = ice(i)     + i1

   if ( do_sedi_heat )  &
   call sedi_heat(ktop, kbot, dp1, m1_sol, dz1, tz, qvz, qlz, qrz, qiz, qsz, qgz, c_ice)

!-------------------------------------------
! Time-split warm rain processes: 2nd pass
!-------------------------------------------
   call warm_rain(dt_rain, ktop, kbot, dp1, dz1, tz, qvz, qlz, qrz, qiz, qsz, qgz, den, denfac, ccn,  &
                  c_praut, rh_rain, vtrz, r1, m1_rain, w1, h_var)
   rain(i) = rain(i) + r1
   do k=ktop, kbot
      m2_rain(i,k) = m2_rain(i,k) + m1_rain(k)
       m2_sol(i,k) = m2_sol(i,k) + m1_sol(k)
      m1(k) = m1(k) + m1_rain(k) + m1_sol(k)
   enddo

!-------------------------
! * ice-phase microphysics
!-------------------------

   call icloud( ktop, kbot, tz, p1, qvz, qlz, qrz, qiz, qsz, qgz,  &
                dp1, den, denfac, vtsz, vtgz, vtrz, qaz, rh_adj, rh_rain, dts, h_var )

1000  continue  ! sub-cycle

   if ( sedi_transport ) then
! Note: dp1 is dry mass; dp0 is the old moist (total) mass
       do k = ktop+1, kbot
          u1(k) = (dp0(k)*u1(k) + m1(k-1)*u1(k-1))/(dp0(k)+m1(k-1))
          v1(k) = (dp0(k)*v1(k) + m1(k-1)*v1(k-1))/(dp0(k)+m1(k-1))
          u_dt(i,j,k) = u_dt(i,j,k) + (u1(k)-u0(k))*rdt
          v_dt(i,j,k) = v_dt(i,j,k) + (v1(k)-v0(k))*rdt
       enddo
   endif
   if ( do_sedi_w ) then
        do k = ktop, kbot
           w(i,j,k) = w1(k)
        enddo
   endif

   do k = ktop, kbot
! update moist air mass (actually hydrostatic pressure)
!-- Convert to dry mixing ratios
         omq = dp1(k) / dp0(k)
      qv_dt(i,j,k) = qv_dt(i,j,k) + rdt*(qvz(k)-qv0(k))*omq
      ql_dt(i,j,k) = ql_dt(i,j,k) + rdt*(qlz(k)-ql0(k))*omq
      qr_dt(i,j,k) = qr_dt(i,j,k) + rdt*(qrz(k)-qr0(k))*omq
      qi_dt(i,j,k) = qi_dt(i,j,k) + rdt*(qiz(k)-qi0(k))*omq
      qs_dt(i,j,k) = qs_dt(i,j,k) + rdt*(qsz(k)-qs0(k))*omq
      qg_dt(i,j,k) = qg_dt(i,j,k) + rdt*(qgz(k)-qg0(k))*omq
      cvm = c_air + qvz(k)*c_vap + (qrz(k)+qlz(k))*c_liq + (qiz(k)+qsz(k)+qgz(k))*c_ice
           pt_dt(i,j,k) = pt_dt(i,j,k) + rdt*(tz(k)- t0(k))*cvm/cp_air
   enddo

   do k = ktop, kbot
      if ( do_qa ) then
           qa_dt(i,j,k) = 0.
      else
           qa_dt(i,j,k) = qa_dt(i,j,k) + rdt*( qaz(k)/real(ntimes)-qa0(k))
      endif
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

   if ( id_droplets> 0 ) then
        do k=ktop, kbot
           qn2(i,j,k) = ccn(k)
        enddo
   endif

2000  continue

 end subroutine mpdrv



 subroutine sedi_heat(ktop, kbot, dm, m1, dz, tz, qv, ql, qr, qi, qs, qg, cw)
 integer, intent(in):: ktop, kbot
! Input q fields are dry mixing ratios, and dm is dry air mass
 real, intent(in   ), dimension(ktop:kbot):: dm, m1, dz, qv, ql, qr, qi, qs, qg
 real, intent(inout), dimension(ktop:kbot):: tz
 real, intent(in):: cw    ! heat capacity
!
 real,dimension(ktop:kbot):: dgz, cvm, cvn
 real:: tmp
 integer:: k

   do k=ktop, kbot
      dgz(k) = -0.5*grav*dz(k)    ! > 0
      cvn(k) = dm(k)*(cv_air + qv(k)*cv_vap + (qr(k)+ql(k))*c_liq + (qi(k)+qs(k)+qg(k))*c_ice)
   enddo

! SJL July 2014
! Assumption: the KE in the falling condensates is negligible compared to the potential energy
! that was unaccounted for. Local thermal equilibrium is assumed, and the loss in PE is transformed
! into internal energy (to heat the whole grid box)
! Backward time-implicit upwind transport scheme:
! dm here is dry air mass

   k = ktop
     tmp = cvn(k) + m1(k)*cw
   tz(k) = (tmp*tz(k) + m1(k)*dgz(k) ) / tmp

! Implicit algorithm: can't be vectorized
! Needs an inner i-loop for vectorization
   do k=ktop+1, kbot
      tz(k) = ( (cvn(k)+cw*(m1(k)-m1(k-1)))*tz(k) + m1(k-1)*cw*tz(k-1) + dgz(k)*(m1(k-1)+m1(k)) )    &
            / (  cvn(k)+cw*m1(k) )
   enddo

 end subroutine sedi_heat


 subroutine warm_rain( dt, ktop, kbot, dp, dz, tz, qv, ql, qr, qi, qs, qg, &
                       den, denfac, ccn, c_praut, rh_rain, vtr, r1, m1_rain, w1, h_var)

 integer, intent(in):: ktop, kbot
 real,    intent(in):: dt                    ! time step (s)
 real,    intent(in), dimension(ktop:kbot):: dp, dz, den, denfac, ccn, c_praut
 real,    intent(in):: rh_rain, h_var
 real, intent(inout), dimension(ktop:kbot):: tz, qv, ql, qr, qi, qs, qg, vtr
 real, intent(out):: r1
 real, intent(inout), dimension(ktop:kbot):: m1_rain, w1

! local:
 real, parameter:: so3 = 7./3.
 real, dimension(ktop:kbot):: dl, dm
 real, dimension(ktop:kbot+1):: ze, zt
 real:: sink, dq, qc0, qc
 real:: qden
 real:: zs = 0.
 real:: dt5
 integer k
!-----------------------------------------------------------------------
! fall velocity constants:
!-----------------------------------------------------------------------
 real, parameter :: vconr = 2503.23638966667
 real, parameter :: normr = 25132741228.7183
 real, parameter :: thr=1.e-8
 logical no_fall

!---------------------
! warm-rain processes:
!---------------------

  dt5 = 0.5*dt

!------------------------
! Terminal speed of rain:
!------------------------

  m1_rain(:) = 0.
  call check_column(ktop, kbot, qr, no_fall)
  if ( no_fall ) then
       vtr(:) = vf_min
       r1 = 0.
       go to 999   ! jump to auto-conversion
  endif

  if ( const_vr ) then
     vtr(:) = vr_fac   ! IFS_2016: 4.0
   else
     do k=ktop, kbot
        qden = qr(k)*den(k)
        if ( qr(k) < thr ) then
             vtr(k) = vr_min
        else
             vtr(k) = vr_fac*vconr*sqrt(min(10., sfcrho/den(k)))*exp(0.2*log(qden/normr))
             vtr(k) = min(vr_max, max(vr_min, vtr(k)))
        endif
     enddo
   endif

  ze(kbot+1) = zs
  do k=kbot, ktop, -1
     ze(k) = ze(k+1) - dz(k)  ! dz<0
  enddo

! Evap_acc of rain for 1/2 time step
!!!  if ( .not. fast_sat_adj )  &
  call revap_racc( ktop, kbot, dt5, tz, qv, ql, qr, qi, qs, qg, den, denfac, rh_rain, h_var )

  if ( do_sedi_w ) then
      do k=ktop, kbot
         dm(k) = dp(k)*(1.+qv(k)+ql(k)+qr(k)+qi(k)+qs(k)+qg(k))
      enddo
  endif

  if ( use_ppm ) then
       zt(ktop) = ze(ktop)
       do k=ktop+1,kbot
          zt(k) = ze(k) - dt5*(vtr(k-1)+vtr(k))
       enddo
       zt(kbot+1) = zs - dt*vtr(kbot)

       do k=ktop,kbot
          if( zt(k+1)>=zt(k) ) zt(k+1) = zt(k) - dz_min
       enddo
       call lagrangian_fall_ppm(ktop, kbot, zs, ze, zt, dp, qr, r1, m1_rain, mono_prof)
  else
       call implicit_fall(dt, ktop, kbot, ze, vtr, dp, qr, r1, m1_rain)
  endif

  if ( do_sedi_w ) then
      w1(ktop) = (dm(ktop)*w1(ktop)+m1_rain(ktop)*vtr(ktop)) /(dm(ktop)-m1_rain(ktop))
      do k=ktop+1, kbot
         w1(k) = (dm(k)*w1(k)-m1_rain(k-1)*vtr(k-1)+m1_rain(k)*vtr(k))   &
               / (dm(k)      +m1_rain(k-1)         -m1_rain(k))
      enddo
  endif
  if ( do_sedi_heat )  &
  call sedi_heat(ktop, kbot, dp, m1_rain, dz, tz, qv, ql, qr, qi, qs, qg, c_liq)
! Finish the remaing 1/2 time step
  call revap_racc( ktop, kbot, dt5, tz, qv, ql, qr, qi, qs, qg, den, denfac, rh_rain, h_var )

999  continue

!-------------------
! * auto-conversion
!-------------------
! Assuming linear subgrid vertical distribution of cloud water
! following Lin et al. 1994, MWR

  call linear_prof( kbot-ktop+1, ql(ktop), dl(ktop), z_slope_liq, h_var )

  do k=ktop,kbot
    qc0 = fac_rc*ccn(k)
    if ( tz(k) > t_wfr + dt_fr ) then
      dl(k) = min(max(1.e-6, dl(k)), 0.5*ql(k))

!--------------------------------------------------------------------
!  As in Klein's GFDL AM2 stratiform scheme (with subgrid variations)
!--------------------------------------------------------------------
      if ( use_ccn ) then
!  CCN is formulted as CCN = CCN_surface * (den/den_surface)
           qc = qc0
      else
           qc = qc0/den(k)
      endif

      dq = 0.5*(ql(k) + dl(k) - qc)    ! dq = dl  if qc==q_minus = ql - dl
                                       ! dq = 0   if qc==q_plus  = ql + dl
      if ( dq > 0. ) then    ! q_plus > qc
! Revised continuous form: linearly decays (with subgrid dl) to zero at qc == ql + dl
            sink = min(1.,dq/dl(k)) * dt*c_praut(k)*den(k)*exp(so3*log(ql(k)))
           ql(k) = ql(k) - sink
           qr(k) = qr(k) + sink
      endif
    endif
  enddo


 end subroutine warm_rain

 subroutine revap_racc( ktop, kbot, dt, tz, qv, ql, qr, qi, qs, qg, den, denfac, rh_rain, h_var)
 integer, intent(in):: ktop, kbot
 real,    intent(in):: dt                 ! time step (s)
 real,    intent(in), dimension(ktop:kbot):: den, denfac
 real,    intent(in)                      :: rh_rain, h_var
 real, intent(inout), dimension(ktop:kbot):: tz, qv, qr, ql, qi, qs, qg
! local:
 real:: dqv, qsat, dqsdt, evap, t2, qden, q_plus, q_minus, sink
 real:: qpz, dq, dqh, tin, lcpk
 integer k

  do k=ktop,kbot
     if ( tz(k) > t_wfr .and. qr(k) > qrmin ) then

      lcpk = (lv00+dc_vap*tz(k)) / (c_air + qv(k)*c_vap + (qr(k)+ql(k))*c_liq + &
                (qi(k)+qs(k)+qg(k))*c_ice)

             tin = tz(k) - lcpk*ql(k) ! presence of clouds suppresses the rain evap
             qpz = qv(k) + ql(k)
            qsat = wqs2(tin, den(k), dqsdt)
             dqh = max(ql(k), h_var*max(qpz, qcmin))
             dqh = min( dqh, 0.2*qpz ) ! New limiter
             dqv = qsat - qv(k)   ! use this to prevent super-sat the gird box
         q_minus = qpz - dqh
         q_plus  = qpz + dqh

! qsat must be > q_minus to activate evaporation
! qsat must be < q_plus  to activate accretion

!-------------------
! * Rain evaporation
!-------------------
         if ( dqv>qvmin .and. qsat > q_minus ) then
              if ( qsat > q_plus ) then
                   dq = qsat - qpz
              else
! q_minus < qsat < q_plus
! dq == dqh if qsat == q_minus
                  dq = 0.25*(q_minus-qsat)**2 / dqh
              endif
              qden = qr(k)*den(k)
               t2 = tin*tin
              evap =  crevp(1)*t2*dq*(crevp(2)*sqrt(qden)+crevp(3)*exp(0.725*log(qden)))   &
                   / (crevp(4)*t2 + crevp(5)*qsat*den(k))
              evap = min( qr(k), dt*evap, dqv/(1.+lcpk*dqsdt) )
! Alternative Minimum Evap in dry environmental air
              sink = min( qr(k), dim(rh_rain*qsat, qv(k))/(1.+lcpk*dqsdt) )
              evap = max( evap, sink )
             qr(k) = qr(k) - evap
             qv(k) = qv(k) + evap
             tz(k) = tz(k) - evap*lcpk
         endif

!        if ( qr(k)>qrmin .and. ql(k)>1.E-7  .and.  qsat<q_plus ) then
         if ( qr(k)>qrmin .and. ql(k)>1.E-6 .and. qsat<q_minus ) then
!-------------------
! * Accretion: pracc
!-------------------
               sink = dt*denfac(k)*cracw * exp(0.95*log(qr(k)*den(k)))
               sink = sink/(1.+sink)*ql(k)
              ql(k) = ql(k) - sink
              qr(k) = qr(k) + sink
         endif

     endif   ! warm-rain
  enddo

 end subroutine revap_racc


 subroutine linear_prof(km, q, dm, z_var, h_var)
! Used for cloud ice and cloud water autoconversion
! qi --> ql  & ql --> qr
! Edges: qE == qbar +/- dm
 integer, intent(in):: km
 real, intent(in ):: q(km), h_var
 real, intent(out):: dm(km)
 logical, intent(in):: z_var
!
 real:: dq(km)
 integer:: k

 if ( z_var ) then
    do k=2,km
       dq(k) = 0.5*(q(k) - q(k-1))
    enddo
    dm(1) = 0.
    do k=2, km-1
! Use twice the strength of the positive definiteness limiter (Lin et al 1994)
       dm(k) = 0.5*min(abs(dq(k) + dq(k+1)), 0.5*q(k))
       if ( dq(k)*dq(k+1) <= 0. ) then
            if ( dq(k) > 0. ) then   ! Local max
                 dm(k) = min( dm(k), dq(k), -dq(k+1) )
            else
                 dm(k) = 0.
            endif
       endif
    enddo
    dm(km) = 0.
! impose a presumed background horizontal variability that is proportional to the value itself
    do k=1, km
       dm(k) = max( dm(k), qvmin, h_var*q(k) )
    enddo
 else
    do k=1, km
       dm(k) = max( qvmin, h_var*q(k) )
    enddo
 endif

 end subroutine linear_prof


 subroutine icloud(ktop, kbot, tzk, p1, qvk, qlk, qrk, qik, qsk, qgk, dp1, &
                   den, denfac, vts, vtg, vtr, qak, rh_adj, rh_rain, dts, h_var)

!----------------------------------------------------
! Bulk cloud micro-physics; processes splitting
! with some un-split sub-grouping
! Time implicit (when possible) accretion and autoconversion
! Author: Shian-Jiann Lin, GFDL
!-------------------------------------------------------

 integer, intent(in) :: ktop, kbot
 real, intent(in),    dimension(ktop:kbot):: p1, dp1, den, denfac, vts, vtg, vtr
 real, intent(inout), dimension(ktop:kbot):: tzk, qvk, qlk, qrk, qik, qsk, qgk, qak
 real, intent(in) :: rh_adj, rh_rain, dts, h_var
! local:
 real, dimension(ktop:kbot) :: lcpk, icpk, tcpk, di
 real:: rdts, fac_g2v, fac_v2g, fac_i2s, fac_mlt
 real:: tz, qv, ql, qr, qi, qs, qg, melt
 real:: pracs, psacw, pgacw, pgmlt,    &
        psmlt, psacr, pgacr, pgfr,     &
        pgaut, pgaci, praci, psaut, psaci, pgsub
 real:: tc, tsq, dqs0, qden, qim, qsm, pssub
 real:: dt5, factor, sink, qi_crt
 real:: tmp1, qsw, qsi, dqsdt, dq
 real:: dtmp, qc, q_plus, q_minus, cvm
 integer:: k

 fac_i2s = 1. - exp( -dts/tau_i2s )
 fac_g2v = 1. - exp( -dts/tau_g2v )
 fac_v2g = 1. - exp( -dts/tau_v2g )
 dt5 = 0.5*dts
 fac_mlt = 1. - exp( -dt5/tau_mlt )

 rdts = 1./dts

 do k=ktop,kbot
    icpk(k) = (li00+dc_ice*tzk(k)) /        &
              (c_air+qvk(k)*c_vap+(qlk(k)+qrk(k))*c_liq+(qik(k)+qsk(k)+qgk(k))*c_ice)
 enddo

! Sources of cloud ice: pihom, cold rain, and the sat_adj
! (initiation plus deposition)

! Sources of snow: cold rain, auto conversion + accretion (from cloud ice)
! sat_adj (deposition; requires pre-existing snow); initial snow comes from auto conversion

 do k=ktop, kbot
    if( tzk(k) > tice .and. qik(k) > qcmin ) then
!--------------------------------------
! * pimlt: instant melting of cloud ice
!--------------------------------------
        melt = min( qik(k), fac_mlt*(tzk(k)-tice)/icpk(k) )
        tmp1 = min( melt, dim(ql_mlt, qlk(k)) )   ! max ql amount
        qlk(k) = qlk(k) + tmp1
        qrk(k) = qrk(k) + melt - tmp1
        qik(k) = qik(k) - melt
        tzk(k) = tzk(k) - melt*icpk(k)
    elseif (tzk(k) < t_wfr .and. qlk(k) > qcmin ) then
!--------------------------------------------------------------
! * pihom * homogeneous Freezing of cloud water into cloud ice:
!--------------------------------------------------------------
! This is the 1st occurance of liquid water freezing in the split MP process
          dtmp = t_wfr - tzk(k)
        factor = min( 1., dtmp/dt_fr )
          sink = min( qlk(k)*factor, dtmp/icpk(k) )
! DE_ICE_GFS: prevent large ice build-up from SAS
          qi_crt = qi_gen*min(qi_lim, 0.1*(tice-tzk(k))) / den(k)
          tmp1 = min(sink, dim(qi_crt, qik(k)))
            qlk(k) = qlk(k) - sink
            qsk(k) = qsk(k) + sink - tmp1
            qik(k) = qik(k) + tmp1
            tzk(k) = tzk(k) + sink*icpk(k)
    endif
 enddo

 call linear_prof( kbot-ktop+1, qik(ktop), di(ktop), z_slope_ice, h_var )

 do k=ktop,kbot
    cvm = c_air + qvk(k)*c_vap + (qlk(k)+qrk(k))*c_liq + (qik(k)+qsk(k)+qgk(k))*c_ice
    lcpk(k) = (lv00+dc_vap*tzk(k)) / cvm
    icpk(k) = (li00+dc_ice*tzk(k)) / cvm
    tcpk(k) = lcpk(k) + icpk(k)
 enddo

 do 3000 k=ktop, kbot

   if( tzk(k) < t_min .or. p1(k) < p_min ) goto 3000

   tz = tzk(k)
   qv = qvk(k)
   ql = qlk(k)
   qi = qik(k)
   qr = qrk(k)
   qs = qsk(k)
   qg = qgk(k)

   pgacr = 0.
   pgacw = 0.
   tc = tz-tice

if ( tc .ge. 0. ) then

!-----------------------------
!* Melting of snow and graupel
!-----------------------------
     dqs0 = ces0/p1(k) - qv

     if( qs>qcmin ) then

! * accretion: cloud water --> snow
! only rate is used (for snow melt) since tc > 0.
        if( ql>qrmin ) then
            factor = denfac(k)*csacw*exp(0.8125*log(qs*den(k)))
             psacw = factor/(1.+dts*factor)*ql     ! rate
        else
             psacw = 0.
        endif

        if ( qr>qrmin ) then
! * accretion: melted snow --> rain:
             psacr = min(acr3d(vts(k), vtr(k), qr, qs, csacr, acco(1,2), den(k)), qr*rdts)
! * accretion: snow --> rain
             pracs = acr3d(vtr(k), vts(k), qs, qr, cracs, acco(1,1), den(k))
        else
             psacr = 0.
             pracs = 0.
        endif

! Total snow sink:
! * Snow melt (due to rain accretion): snow --> rain
        psmlt = max(0., smlt(tc, dqs0, qs*den(k), psacw, psacr, csmlt, den(k), denfac(k)))
         sink = min(qs, dts*(psmlt+pracs), tc/icpk(k))

        qs = qs - sink
        qr = qr + sink
        tz = tz - sink*icpk(k)    ! cooling due to snow melting
        tc = tz-tice
     endif

     if ( qg>qcmin .and. tc>0. ) then

! * accretion: rain --> graupel
         if ( qr>qrmin )   &
         pgacr = min(acr3d(vtg(k), vtr(k), qr, qg, cgacr, acco(1,3), den(k)), rdts*qr)

         qden = qg*den(k)
         if( ql>qrmin ) then
! * accretion: cloud water --> graupel
             factor = cgacw*qden/sqrt(den(k)*sqrt(sqrt(qden)))
              pgacw = factor/(1.+dts*factor) * ql  ! rate
         endif

! * melting: graupel --> rain
         pgmlt = dts*gmlt(tc, dqs0, qden, pgacw, pgacr, cgmlt, den(k))
         pgmlt = min( max(0., pgmlt), qg, tc/icpk(k) )
            qg = qg - pgmlt
            qr = qr + pgmlt
            tz = tz - pgmlt*icpk(k)
     endif   ! graupel existed

else

!------------------
! Cloud ice proc:
!------------------
  if ( qi>3.E-7 ) then  ! cloud ice sink terms

!----------------------------------------
! * accretion (pacr): cloud ice --> snow
!----------------------------------------
     if ( qs>1.E-7 )  then
! SJL added (following Lin Eq. 23) the temperature dependency
! To reduce accretion, use Esi = exp(0.05*tc) as in Hong et al 2004
          factor = dts*denfac(k)*csaci*exp(0.05*tc + 0.8125*log(qs*den(k)))
          psaci = factor/(1.+factor) * qi
     else
          psaci = 0.
     endif

!-------------------------------------
! * autoconversion: cloud ice --> snow
!-------------------------------------
! Similar to LFO 1983: Eq. 21 solved implicitly
! Threshold from WSM6 scheme, Hong et al 2004, Eq(13) : qi0_crt ~0.8E-4
    qim = qi0_crt / den(k)
! Assuming linear subgrid vertical distribution of cloud ice
! The mismatch computation following Lin et al. 1994, MWR

    if ( const_vi ) then
        tmp1 = fac_i2s
    else
        tmp1 = fac_i2s * exp(0.025*tc)
    endif

    di(k) = max( di(k), qrmin )
    q_plus = qi + di(k)
    if ( q_plus > (qim+qrmin) ) then
         if ( qim > (qi - di(k)) ) then
              dq = (0.25*(q_plus-qim)**2) / di(k)
         else
              dq = qi - qim
         endif
         psaut  = tmp1 * dq
    else
         psaut = 0.
    endif
! sink is no greater than 50% of qi
    sink = min( 0.5*qi, psaci+psaut )
      qi = qi - sink
      qs = qs + sink

!-----------------------------------
! * accretion: cloud ice --> graupel
!-----------------------------------
    if ( qg>1.E-6 ) then
!        factor = dts*cgaci/sqrt(den(k))*exp(0.05*tc + 0.875*log(qg*den(k)))
! Simplified form: remove temp dependency & set the exponent "0.875" --> 1
         factor = dts*cgaci*sqrt(den(k))*qg
          pgaci = factor/(1.+factor)*qi
             qi = qi - pgaci
             qg = qg + pgaci
    endif

  endif  ! significant cloud ice existed

!----------------
! Cold-Rain proc:
!----------------
! rain to ice, snow, graupel processes:

  tc = tz-tice

  if ( qr>1.e-7 .and. tc < 0. ) then

! *sink* terms to qr: psacr + pgfr
! source terms to qs: psacr
! source terms to qg: pgfr

!----------------------------
! * accretion of rain by snow
!----------------------------
      if ( qs > 1.E-7 ) then   ! if snow exists
           psacr = dts*acr3d(vts(k), vtr(k), qr, qs, csacr, acco(1,2), den(k))
      else
           psacr = 0.
      endif

!----------------------------
! * rain freezing --> graupel
!----------------------------
     pgfr = dts*cgfr(1)/den(k)*(exp(-cgfr(2)*tc)-1.)*exp( 1.75*log(qr*den(k)) )

!--- Total sink to qr
       sink = psacr + pgfr
     factor = min( sink, qr, -tc/icpk(k) ) / max( sink, qrmin )

      psacr = factor * psacr
      pgfr  = factor * pgfr

      sink = psacr + pgfr
        tz = tz + sink*icpk(k)
        qr = qr - sink
        qs = qs + psacr
        qg = qg + pgfr
  endif  ! qr existed

!--------------------------
! Graupel production terms:
!--------------------------

  if( qs > 1.e-7 ) then
! * accretion: snow --> graupel
      if ( qg > qrmin ) then
           sink = dts*acr3d(vtg(k), vts(k), qs, qg, cgacs, acco(1,4), den(k))
      else
           sink = 0.
      endif

      qsm = qs0_crt / den(k)
      if ( qs > qsm ) then
! * Autoconversion Snow --> graupel
           factor = dts*1.e-3*exp(0.09*(tz-tice))
             sink = sink + factor/(1.+factor)*(qs-qsm)
      endif
      sink = min( qs, sink )
        qs = qs - sink
        qg = qg + sink

  endif   ! snow existed

  if ( qg>1.E-7 .and. tz < tice0 ) then

! * accretion: cloud water --> graupel
     if( ql>1.E-6 ) then
           qden = qg*den(k)
         factor = dts*cgacw*qden/sqrt(den(k)*sqrt(sqrt(qden)))
          pgacw = factor/(1.+factor)*ql
     else
          pgacw = 0.
     endif

! * accretion: rain --> graupel
     if ( qr>1.E-6 ) then
          pgacr = min(dts*acr3d(vtg(k), vtr(k), qr, qg, cgacr, acco(1,3), den(k)), qr)
     else
          pgacr = 0.
     endif

       sink = pgacr + pgacw
     factor = min( sink, dim(tice,tz)/icpk(k) ) / max( sink, qrmin )
      pgacr = factor * pgacr
      pgacw = factor * pgacw

     sink = pgacr + pgacw
       tz = tz + sink*icpk(k)
       qg = qg + sink
       qr = qr - pgacr
       ql = ql - pgacw
  endif    ! graupel existed

endif   ! end ice-physics

     tzk(k) = tz
     qvk(k) = qv
     qlk(k) = ql
     qik(k) = qi
     qrk(k) = qr
     qsk(k) = qs
     qgk(k) = qg

3000 continue   ! k-loop

   call subgrid_z_proc(ktop, kbot, p1, den, denfac, dts, rh_adj, tzk, qvk, qlk, qrk, qik, qsk, qgk, qak, h_var, rh_rain)

 end subroutine icloud


 subroutine subgrid_z_proc(ktop, kbot, p1, den, denfac, dts, rh_adj, tz, qv, ql, qr, qi, qs, qg, qa, h_var, rh_rain)

! Temperature sentive high vertical resolution processes:

 integer, intent(in):: ktop, kbot
 real, intent(in),    dimension(ktop:kbot):: p1, den, denfac
 real, intent(in)                         :: dts, rh_adj, h_var, rh_rain
 real, intent(inout), dimension(ktop:kbot):: tz, qv, ql, qr, qi, qs, qg, qa
! local:
 real, dimension(ktop:kbot):: lcpk, icpk, tcpk, tcp3
 real:: fac_v2l, fac_l2v
 real:: pidep, qi_crt
! qstar over water may be accurate only down to -80 C with ~10% uncertainty
! must not be too large to allow PSC
 real:: rh, clouds,  rqi, tin, qsw, qsi, qpz, qstar
 real:: dqsdt, dwsdt, dq, dq0, factor, tmp, cvm
 real:: q_plus, q_minus, dt_evap, dt_pisub
 real:: evap, sink, tc, pisub, iwt, q_adj, dtmp, lat2
 real:: pssub, pgsub, tsq, qden, fac_g2v, fac_v2g
 integer :: k

 if ( fast_sat_adj ) then
    dt_evap = 0.5*dts
 else
    dt_evap = dts
 endif

 fac_v2l = 1. - exp( -dt_evap/tau_v2l )        !
 fac_l2v = 1. - exp( -dt_evap/tau_l2v )        !

 fac_g2v = 1. - exp( -dts/tau_g2v )
 fac_v2g = 1. - exp( -dts/tau_v2g )

 lat2 = lats * lats

 do k=ktop, kbot
! Moist heat capacity
    cvm = c_air + qv(k)*c_vap + (ql(k)+qr(k))*c_liq + (qi(k)+qs(k)+qg(k))*c_ice
    lcpk(k) = (lv00+dc_vap*tz(k)) / cvm
    icpk(k) = (li00+dc_ice*tz(k)) / cvm
    tcpk(k) = lcpk(k) + icpk(k)
    tcp3(k) = lcpk(k) + icpk(k)*min(1., dim(tice,tz(k))/(tice-t_wfr))
 enddo

 do 4000 k=ktop,kbot

   if ( p1(k) < p_min ) goto 4000

   if ( tz(k) < t_min ) then
        sink = dim(qv(k), 1.e-7)
        qv(k) = qv(k) - sink
        qi(k) = qi(k) + sink
        tz(k) = tz(k) + sink*tcpk(k)
        if ( .not.do_qa ) qa(k) = qa(k) + 1.  ! Air fully saturated; 100 % cloud cover
        goto 4000
    endif

! Instant evaporation/sublimation of all clouds if RH<rh_adj --> cloud free
      iwt = qi(k)
   clouds = ql(k) + iwt

   tin = tz(k) - ( lcpk(k)*clouds + icpk(k)*iwt )  ! minimum  temperature
   if ( tin>t_sub+6. ) then
   qpz = qv(k) + clouds
       rh = qpz / iqs1(tin, den(k))

    if ( rh < rh_adj ) then  ! qpz / rh_adj < qs
         tz(k) = tin
         qv(k) = qpz
         ql(k) = 0.
         qi(k) = 0.
         goto 4000            ! cloud free
    endif
    endif

! cloud water <--> vapor adjustment:
   qsw = wqs2(tz(k), den(k), dwsdt)
   dq0 = qsw - qv(k)
   if ( dq0 > 0. ) then
        evap =  min( ql(k), fac_l2v*dq0/(1.+tcp3(k)*dwsdt) )
   else   ! condensate all excess vapor into cloud water
        evap = fac_v2l*dq0/(1.+tcp3(k)*dwsdt)
   endif
   qv(k) = qv(k) + evap
   ql(k) = ql(k) - evap
   tz(k) = tz(k) - evap*lcpk(k)

! Enforce complete freezing below -48 C
   dtmp = t_wfr - tz(k)   ! [-40,-48]
   if( dtmp>0. .and. ql(k) > qcmin ) then
       sink = min( ql(k),  ql(k)*dtmp*0.125, dtmp/icpk(k) )
      ql(k) = ql(k) - sink
      qi(k) = qi(k) + sink
      tz(k) = tz(k) + sink*icpk(k)
   endif

  if ( fast_sat_adj ) then
      dt_pisub = 0.5*dts
  else
      dt_pisub = dts
      tc = tice - tz(k)
      if( ql(k) > qrmin .and. tc> 0. ) then
! Bigg mechanism
         sink = dts*3.3333e-10*(exp(0.66*tc)-1.)*den(k)*ql(k)*ql(k)
         sink = min(ql(k), tc/icpk(k), sink)
         ql(k) = ql(k) - sink
         qi(k) = qi(k) + sink
         tz(k) = tz(k) + sink*icpk(k)
      endif ! significant ql existed
  endif

   cvm = c_air + qv(k)*c_vap + (ql(k)+qr(k))*c_liq + (qi(k)+qs(k)+qg(k))*c_ice
   lcpk(k) = (lv00+dc_vap*tz(k)) / cvm
   icpk(k) = (li00+dc_ice*tz(k)) / cvm
   tcpk(k) = lcpk(k) + icpk(k)

!------------------------------------------
! * pidep: sublimation/deposition of ice:
!------------------------------------------
    if ( tz(k) < tice ) then
         qsi = iqs2(tz(k), den(k), dqsdt)
          dq = qv(k) - qsi
        sink = dq/(1.+tcpk(k)*dqsdt)

        if ( qi(k) > qrmin ) then
! Eq 9, Hong et al. 2004, MWR
! For A and B, see Dudhia 1989: page 3103 Eq (B7) and (B8)
             pidep = dt_pisub*dq*349138.78*exp(0.875*log(qi(k)*den(k)))     &
                   / (qsi*den(k)*lat2/(0.0243*rvgas*tz(k)**2) + 4.42478e4)
        else
             pidep = 0.
        endif

        if ( dq > 0. ) then   ! vapor -> ice
             tmp = tice - tz(k)
! 20160912: the following should produce more ice at higher altitude
!            qi_crt = 4.92e-11*exp( 1.33*log(1.e3*exp(0.1*tmp)) ) / den(k)
             qi_crt = qi_gen*min(qi_lim, 0.1*tmp) / den(k)
             sink = min(sink, max(qi_crt-qi(k),pidep), tmp/tcpk(k))
        else  ! ice --> vapor
             pidep = pidep * min(1., dim(tz(k),t_sub)*0.2)
             sink = max(pidep, sink, -qi(k))
        endif
        qv(k) = qv(k) - sink
        qi(k) = qi(k) + sink
        tz(k) = tz(k) + sink*tcpk(k)
    endif

!----------------------------------
! * sublimation/deposition of snow:
!----------------------------------
! This process happens for all temp rage
  if ( qs(k)>qrmin ) then
       qsi = iqs2(tz(k), den(k), dqsdt)
          qden = qs(k)*den(k)
           tmp = exp(0.65625*log(qden))
           tsq = tz(k)*tz(k)
            dq = (qsi-qv(k))/(1.+tcpk(k)*dqsdt)
         pssub =  cssub(1)*tsq*(cssub(2)*sqrt(qden) + cssub(3)*tmp*sqrt(denfac(k)))  &
               / (cssub(4)*tsq+cssub(5)*qsi*den(k))
         pssub = (qsi-qv(k))*dts*pssub
         if ( pssub > 0. ) then     ! qs --> qv,  sublimation
              pssub = min(pssub*min(1.,dim(tz(k),t_sub)*0.2), qs(k))
! Minimum sublimation:
!             pssub = max(pssub, min(qs(k), dim(0.35*qsi,qv(k))/(1.+tcpk(k)*dqsdt)) )
         else
              if ( tz(k)>tice ) then
                   pssub = 0.  ! no deposition
              else
                   pssub = max( pssub, 0.5*dq, (tz(k)-tice)/tcpk(k) )
              endif
         endif
         qs(k) = qs(k) - pssub
         qv(k) = qv(k) + pssub
         tz(k) = tz(k) - pssub*tcpk(k)
  endif

!------------------------------------------------------------
! * Simplified 2-way grapuel sublimation-deposition mechanism
!------------------------------------------------------------
  if ( qg(k)>qrmin ) then
         qsi = iqs2(tz(k), den(k), dqsdt)
          dq = (qv(k)-qsi)/(1.+tcpk(k)*dqsdt)
       pgsub = (qv(k)/qsi-1.) * qg(k)
       if ( pgsub > 0. ) then        ! deposition
            if ( tz(k)>tice ) then
                 pgsub = 0.   ! no deposition
            else
                 pgsub = min(fac_v2g*pgsub, 0.2*dq, ql(k)+qr(k), (tice-tz(k))/tcpk(k))
            endif
       else                   ! submilation
            pgsub = max( fac_g2v*pgsub, dq )*min(1.,dim(tz(k),t_sub)*0.1)
! Minimum sublimation:
!           pgsub = min(pgsub, -min(qg(k), dim(0.3*qsi,qv(k))/(1.+tcpk(k)*dqsdt)) )
       endif
       qg(k) = qg(k) + pgsub
       qv(k) = qv(k) - pgsub
       tz(k) = tz(k) + pgsub*tcpk(k)
  endif

!------------------------------------------------
! * Minimum Evap of rain in dry environmental air
!------------------------------------------------
  if( qr(k)>qcmin) then
      qsw = wqs2(tz(k), den(k), dqsdt)
     sink = min(qr(k), dim(rh_rain*qsw, qv(k))/(1.+lcpk(k)*dqsdt))
      qv(k) = qv(k) + sink
      qr(k) = qr(k) - sink
      tz(k) = tz(k) - sink*lcpk(k)
  endif

   if ( do_qa ) goto 4000

   if ( rad_snow ) then
        iwt = qi(k) + qs(k)
   else
        iwt = qi(k)
   endif
   if ( rad_rain ) then
        clouds = ql(k) + qr(k) + iwt
   else
        clouds = ql(k) + iwt
   endif

   qpz = qv(k) + clouds     ! qpz is conserved
   tin = tz(k) - ( lcpk(k)*clouds + icpk(k)*iwt )  ! minimum  temperature

!--------------------
! * determine qstar
!--------------------
! Using the "liquid-frozen water temperature": tin
   if( tin <= t_wfr ) then
       qstar = iqs1(tin, den(k))
   elseif ( tin >= tice ) then
       qstar = wqs1(tin, den(k))
   else
! mixed phase:
       qsi = iqs1(tin, den(k))
       qsw = wqs1(tin, den(k))
       if( clouds > 3.E-6 ) then
           rqi = iwt / clouds
       else
! Mostly liquid water clouds at initial cloud development stage
           rqi = (tice-tin)/(tice-t_wfr)
       endif
       qstar = rqi*qsi + (1.-rqi)*qsw
   endif

!-------------------------
! * cloud fraction
!-------------------------
! Assuming subgrid linear distribution in horizontal; this is effectively a smoother for the
! binary cloud scheme

   if ( qpz > qrmin ) then
! Partial cloudiness by PDF:
            dq = max(qcmin, h_var*qpz)
       q_plus  = qpz + dq        ! cloud free if qstar > q_plus
       q_minus = qpz - dq
       if ( qstar < q_minus ) then
            qa(k) = qa(k) + 1.       ! Air fully saturated; 100 % cloud cover
       elseif ( qstar<q_plus .and. clouds>qc_crt ) then
            qa(k) = qa(k) + (q_plus-qstar)/(dq+dq)       ! partial cloud cover
       endif
   endif

4000 continue

 end subroutine subgrid_z_proc


 subroutine revap_rac1(hydrostatic, is, ie, dt, tz, qv, ql, qr, qi, qs, qg, den, hvar)
 logical, intent(in):: hydrostatic
 integer, intent(in):: is, ie
 real, intent(in):: dt       ! time step (s)
 real, intent(in),    dimension(is:ie):: den, hvar, qi, qs, qg
 real, intent(inout), dimension(is:ie):: tz, qv, qr, ql
! local:
 real, dimension(is:ie):: lcp2, denfac
 real:: dqv, qsat, dqsdt, evap, qden, q_plus, q_minus, sink
 real:: tin, t2, qpz, dq, dqh, q_liq, q_sol
 integer i

  if ( hydrostatic ) then
  do i=is, ie
     q_liq = ql(i) + qr(i)
     q_sol = qi(i) + qs(i) + qg(i)
     lcp2(i) = (Lv0+dc_vap*tz(i))/((1.-(qv(i)+q_liq+q_sol))*cp_air+qv(i)*cp_vap+q_liq*c_liq+q_sol*c_ice)
!    denfac(i) = sqrt(sfcrho/den(i))
  enddo
  else
  do i=is, ie
     q_liq = ql(i) + qr(i)
     q_sol = qi(i) + qs(i) + qg(i)
     lcp2(i) = (Lv0+dc_vap*tz(i))/((1.-(qv(i)+q_liq+q_sol))*cv_air+qv(i)*cv_vap+q_liq*c_liq+q_sol*c_ice)
!    denfac(i) = sqrt(sfcrho/den(i))
  enddo
  endif

  do i=is, ie
     if ( qr(i) > qrmin .and. tz(i) > t_wfr ) then
          qpz = qv(i) + ql(i)
          tin = tz(i) - lcp2(i)*ql(i) ! presence of clouds suppresses the rain evap
         qsat = wqs2(tin, den(i), dqsdt)
          dqh = max(ql(i), hvar(i)*max(qpz, qcmin))
          dqv = qsat - qv(i)
         q_minus = qpz - dqh
         q_plus  = qpz + dqh

! qsat must be > q_minus to activate evaporation
! qsat must be < q_plus  to activate accretion

!-------------------
! * Rain evaporation
!-------------------
         if ( dqv > qvmin  .and. qsat > q_minus ) then
              if ( qsat > q_plus ) then
                   dq = qsat - qpz
              else
! q_minus < qsat < q_plus
! dq == dqh if qsat == q_minus
                  dq = 0.25*(q_minus-qsat)**2 / dqh
              endif
              qden = qr(i)*den(i)
              t2 = tin * tin
              evap =  crevp(1)*t2*dq*(crevp(2)*sqrt(qden)+crevp(3)*exp(0.725*log(qden)))   &
                   / (crevp(4)*t2  +  crevp(5)*qsat*den(i))
              evap = min( qr(i), dt*evap, dqv/(1.+lcp2(i)*dqsdt) )
             qr(i) = qr(i) - evap
             qv(i) = qv(i) + evap
             tz(i) = tz(i) - evap*lcp2(i)
         endif

         if ( qr(i)>qrmin .and. ql(i)>1.E-8 .and. qsat<q_plus ) then
!-------------------
! * Accretion: pracc
!-------------------
               denfac(i) = sqrt(sfcrho/den(i))
               sink = dt*denfac(i)*cracw*exp(0.95*log(qr(i)*den(i)))
               sink = sink/(1.+sink)*ql(i)
              ql(i) = ql(i) - sink
              qr(i) = qr(i) + sink
         endif
     endif
  enddo

 end subroutine revap_rac1

 subroutine terminal_fall(dtm, ktop, kbot, tz, qv, ql, qr, qg, qs, qi, dz, dp,  &
                          den, vtg, vts, vti, r1, g1, s1, i1, m1_sol, w1)

! lagrangian control-volume method:

 real,    intent(in):: dtm                    ! time step (s)
 integer, intent(in):: ktop, kbot
 real,    intent(in), dimension(ktop:kbot):: vtg, vts, vti, den, dp, dz
 real,    intent(inout), dimension(ktop:kbot):: qv, ql, qr, qg, qs, qi, tz, m1_sol, w1
 real,    intent(out):: r1, g1, s1, i1
! local:
 real, dimension(ktop:kbot+1):: ze, zt
 real:: qsat, dqsdt, dt5, melt, evap, dtime
 real:: factor, frac
 real:: tmp1, precip, tc, sink
 real, dimension(ktop:kbot):: lcpk, icpk
 real, dimension(ktop:kbot):: m1, dm
 real:: zs = 0.
 real:: cvm, fac_mlt
 integer k, k0, m
 logical no_fall

 
  do k=ktop,kbot
     m1_sol(k) = 0.
     cvm = c_air + qv(k)*c_vap + (qr(k)+ql(k))*c_liq + (qi(k)+qs(k)+qg(k))*c_ice
     lcpk(k) = (lv00+dc_vap*tz(k)) / cvm
     icpk(k) = (li00+dc_ice*tz(k)) / cvm
  enddo

  dt5 = 0.5*dtm
  fac_mlt = 1. - exp( -dt5/tau_mlt )


! Melting of cloud_ice and snow (before fall):

! find significant melting level
  k0 = kbot
  do k=ktop, kbot-1
     if ( tz(k) > tice ) then
          k0 = k
          go to 11
     endif
  enddo
11  continue

 do k=k0, kbot
!------
! * ice
!------
    tc = tz(k) - tice
    if( qi(k) > qcmin .and. tc>0. ) then
        melt = min( qi(k), fac_mlt*tc/icpk(k) )
        tmp1 = min( melt, dim(ql_mlt, ql(k)) )
        ql(k) = ql(k) + tmp1
        qr(k) = qr(k) + melt - tmp1
        qi(k) = qi(k) - melt
        tz(k) = tz(k) - melt*icpk(k)
           tc = tz(k) - tice
    endif
 enddo

  if ( dtm < 60. ) k0 = kbot
  k0 = kbot

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
             if ( zt(k)<ze(m+1) .and. tz(m)>tice ) then
                  dtime = min( 1.0, (ze(m)-ze(m+1))/(max(vr_min,vti(k))*tau_mlt) )
                   melt = min( qi(k)*dp(k)/dp(m), dtime*(tz(m)-tice)/icpk(m) )
                   tmp1 = min( melt, dim(ql_mlt, ql(m)) )
                  ql(m) = ql(m) + tmp1
                  qr(m) = qr(m) - tmp1 + melt
                  tz(m) = tz(m) - melt*icpk(m)
                  qi(k) = qi(k) - melt*dp(m)/dp(k)
             endif
          enddo
     endif
  enddo
  endif

  if ( do_sedi_w ) then
      do k=ktop, kbot
         dm(k) = dp(k)*(1.+qv(k)+ql(k)+qr(k)+qi(k)+qs(k)+qg(k))
      enddo
  endif
  if ( use_ppm ) then
       call lagrangian_fall_ppm(ktop, kbot, zs, ze, zt, dp, qi, i1, m1_sol, mono_prof)
  else
       call implicit_fall(dtm, ktop, kbot, ze, vti, dp, qi, i1, m1_sol)
  endif

  if ( do_sedi_w ) then
     w1(ktop) = (dm(ktop)*w1(ktop)+m1_sol(ktop)*vti(ktop)) /(dm(ktop)-m1_sol(ktop))
     do k=ktop+1, kbot
        w1(k) = (dm(k)*w1(k)-m1_sol(k-1)*vti(k-1)+m1_sol(k)*vti(k))   &
              / (dm(k)      +m1_sol(k-1)         -m1_sol(k))
     enddo
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
                  dtime = min( dtm, (ze(m)-ze(m+1))/(vr_min+vts(k)) )
             if ( zt(k)<ze(m+1) .and. tz(m)>tice ) then
                  dtime = min(1., dtime/tau_s)
                   melt = min(qs(k)*dp(k)/dp(m), dtime*(tz(m)-tice)/icpk(m))
                  tz(m) = tz(m) - melt*icpk(m)
                  qs(k) = qs(k) - melt*dp(m)/dp(k)
                  if ( zt(k)<zs ) then
                       r1 = r1 + melt*dp(m)   ! precip as rain
                  else
!                      qr source here will fall next time step (therefore, can evap)
                       qr(m) = qr(m) + melt
                  endif
             endif
             if ( qs(k) < qrmin ) exit
          enddo
     endif
  enddo
  endif

  if ( do_sedi_w ) then
      do k=ktop, kbot
         dm(k) = dp(k)*(1.+qv(k)+ql(k)+qr(k)+qi(k)+qs(k)+qg(k))
      enddo
  endif
  if ( use_ppm ) then
       call lagrangian_fall_ppm(ktop, kbot, zs, ze, zt, dp, qs, s1, m1, mono_prof)
  else
       call implicit_fall(dtm, ktop, kbot, ze, vts, dp, qs, s1, m1)
  endif
  do k=ktop,kbot
     m1_sol(k) = m1_sol(k) + m1(k)
  enddo
  if ( do_sedi_w ) then
     w1(ktop) = (dm(ktop)*w1(ktop)+m1(ktop)*vts(ktop)) /(dm(ktop)-m1(ktop))
     do k=ktop+1, kbot
        w1(k) = (dm(k)*w1(k) - m1(k-1)*vts(k-1) + m1(k)*vts(k))   &
              / (dm(k)       + m1(k-1)          - m1(k))
     enddo
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
             if ( zt(k)<ze(m+1) .and. tz(m)>tice ) then
                  dtime = min(1., dtime/tau_g)
                   melt = min(qg(k)*dp(k)/dp(m), dtime*(tz(m)-tice)/icpk(m))
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

  if ( do_sedi_w ) then
      do k=ktop, kbot
         dm(k) = dp(k)*(1.+qv(k)+ql(k)+qr(k)+qi(k)+qs(k)+qg(k))
      enddo
  endif
  if ( use_ppm ) then
       call lagrangian_fall_ppm(ktop, kbot, zs, ze, zt, dp, qg, g1, m1, mono_prof)
  else
       call implicit_fall(dtm, ktop, kbot, ze, vtg, dp, qg, g1, m1)
  endif
  do k=ktop,kbot
     m1_sol(k) = m1_sol(k) + m1(k)
  enddo
  if ( do_sedi_w ) then
     w1(ktop) = (dm(ktop)*w1(ktop)+m1(ktop)*vtg(ktop)) /(dm(ktop)-m1(ktop))
     do k=ktop+1, kbot
        w1(k) = (dm(k)*w1(k) - m1(k-1)*vtg(k-1) + m1(k)*vtg(k))   &
              / (dm(k)       + m1(k-1)          - m1(k))
     enddo
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


 subroutine implicit_fall(dt, ktop, kbot, ze, vt, dp, q, precip, m1)
! Developed by SJ Lin, 2016
 real,    intent(in):: dt
 integer, intent(in):: ktop, kbot
 real,    intent(in), dimension(ktop:kbot+1):: ze
 real,    intent(in), dimension(ktop:kbot):: vt, dp
 real,    intent(inout), dimension(ktop:kbot):: q
 real,    intent(out), dimension(ktop:kbot):: m1
 real,    intent(out):: precip
! local:
 real, dimension(ktop:kbot):: dz, qm, dd
 integer:: k

  do k=ktop,kbot
     dz(k) = ze(k) - ze(k+1)
     dd(k) = dt * vt(k)
      q(k) = q(k)*dp(k)
  enddo

! Sedimentation: non-vectorizable loop
  qm(ktop) = q(ktop) / (dz(ktop)+dd(ktop))
  do k=ktop+1,kbot
     qm(k) = (q(k) + dd(k-1)*qm(k-1)) / (dz(k) + dd(k))
  enddo
! qm is density at this stage

  do k=ktop,kbot
     qm(k) = qm(k)*dz(k)
  enddo

! Output mass fluxes: non-vectorizable loop
  m1(ktop) = q(ktop) - qm(ktop)
  do k=ktop+1,kbot
     m1(k) = m1(k-1) + q(k) - qm(k)
  enddo
  precip = m1(kbot)

! Update:
  do k=ktop,kbot
     q(k) = qm(k) / dp(k)
  enddo

 end subroutine implicit_fall



 subroutine lagrangian_fall_ppm(ktop, kbot, zs, ze, zt, dp, q, precip, m1, mono)
 integer, intent(in):: ktop, kbot
 real,    intent(in):: zs
 logical, intent(in):: mono
 real,    intent(in), dimension(ktop:kbot+1):: ze, zt
 real,    intent(in), dimension(ktop:kbot):: dp
! m1: flux
 real,    intent(inout), dimension(ktop:kbot):: q, m1
 real,    intent(out):: precip
! local:
 real, dimension(ktop:kbot):: qm, dz
 real a4(4,ktop:kbot)
 real pl, pr, delz, esl
 integer k, k0, n, m
 real, parameter:: r3 = 1./3., r23 = 2./3.

! density:
  do k=ktop,kbot
      dz(k) = zt(k) - zt(k+1)      ! note: dz is positive
       q(k) = q(k) * dp(k)
     a4(1,k) = q(k) / dz(k)
     qm(k) = 0.
  enddo

! Construct vertical profile with zt as coordinate

   call cs_profile(a4(1,ktop), dz(ktop), kbot-ktop+1, mono)

   k0 = ktop
   do k=ktop,kbot
      do n=k0,kbot
      if(ze(k) <= zt(n) .and. ze(k) >= zt(n+1)) then
         pl = (zt(n)-ze(k)) / dz(n)
         if( zt(n+1) <= ze(k+1) ) then
!                          entire new grid is within the original grid
                pr = (zt(n)-ze(k+1)) / dz(n)
            qm(k) = a4(2,n) + 0.5*(a4(4,n)+a4(3,n)-a4(2,n))*(pr+pl) -  &
                     a4(4,n)*r3*(pr*(pr+pl)+pl**2)
            qm(k) = qm(k)*(ze(k)-ze(k+1))
            k0 = n
            goto 555
         else
            qm(k) = (ze(k)-zt(n+1)) * (a4(2,n)+0.5*(a4(4,n)+   &
                      a4(3,n)-a4(2,n))*(1.+pl) - a4(4,n)*( r3*(1.+pl*(1.+pl))) )
            if ( n<kbot ) then
               do m=n+1,kbot
!                                        locate the bottom edge: ze(k+1)
                  if( ze(k+1) < zt(m+1) ) then
                     qm(k) = qm(k) + q(m)
                  else
                     delz = zt(m) - ze(k+1)
                      esl = delz / dz(m)
                     qm(k) = qm(k) + delz*( a4(2,m) + 0.5*esl*        &
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

  m1(ktop) = q(ktop) - qm(ktop)
  do k=ktop+1,kbot
     m1(k) = m1(k-1) + q(k) - qm(k)
  enddo
  precip = m1(kbot)

! Convert back to *dry* mixing ratio:
! dp must be dry air_mass (because moist air mass will be changed due to terminal fall).
   do k=ktop,kbot
     q(k) = qm(k)/dp(k)
   enddo

 end subroutine lagrangian_fall_ppm


 subroutine cs_profile(a4, del, km, do_mono)
 integer, intent(in):: km      ! vertical dimension
 real   , intent(in):: del(km)
 logical, intent(in):: do_mono
 real , intent(inout):: a4(4,km)
!-----------------------------------------------------------------------
 real, parameter:: qp_min = 1.e-6
 real  gam(km)
 real  q(km+1)
 real   d4, bet, a_bot, grat, pmp, lac
 real   pmp_1, lac_1, pmp_2, lac_2
 real  da1, da2, a6da
 integer k
 logical extm(km)

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
  q(1) = max( q(1), 0. )
  q(2) = min( q(2), max(a4(1,1), a4(1,2)) )
  q(2) = max( q(2), min(a4(1,1), a4(1,2)), 0. )

! Interior:
  do k=3,km-1
     if ( gam(k-1)*gam(k+1)>0. ) then
          q(k) = min( q(k), max(a4(1,k-1),a4(1,k)) )
          q(k) = max( q(k), min(a4(1,k-1),a4(1,k)) )
     else
          if ( gam(k-1) > 0. ) then
! There exists a local max
               q(k) = max( q(k), min(a4(1,k-1),a4(1,k)) )
          else
! There exists a local min
               q(k) = min( q(k), max(a4(1,k-1),a4(1,k)) )
               q(k) = max( q(k), 0.0 )
          endif
     endif
  enddo

  q(km  ) = min( q(km), max(a4(1,km-1), a4(1,km)) )
  q(km  ) = max( q(km), min(a4(1,km-1), a4(1,km)), 0. )
! q(km+1) = max( q(km+1), 0.)

!-----------------------------------------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
!-----------------------------------------------------------
  do k=1,km-1
     a4(2,k) = q(k  )
     a4(3,k) = q(k+1)
  enddo

  do k=2,km-1
     if ( gam(k)*gam(k+1) > 0.0 ) then
          extm(k) = .false.
     else
          extm(k) = .true.
     endif
  enddo

  if ( do_mono ) then
     do k=3,km-2
        if ( extm(k) ) then
! positive definite constraint ONLY if true local extrema
           if ( a4(1,k)<qp_min .or. extm(k-1) .or. extm(k+1) ) then
               a4(2,k) = a4(1,k)
               a4(3,k) = a4(1,k)
           endif
        else
           a4(4,k) = 6.*a4(1,k) - 3.*(a4(2,k)+a4(3,k))
           if( abs(a4(4,k)) > abs(a4(2,k)-a4(3,k)) ) then
! Check within the smooth region if subgrid profile is non-monotonic
                pmp_1 = a4(1,k) - 2.0*gam(k+1)
                lac_1 = pmp_1   + 1.5*gam(k+2)
              a4(2,k) = min( max(a4(2,k), min(a4(1,k), pmp_1, lac_1)),  &
                                          max(a4(1,k), pmp_1, lac_1) )
                pmp_2 = a4(1,k) + 2.0*gam(k)
                lac_2 = pmp_2   - 1.5*gam(k-1)
              a4(3,k) = min( max(a4(3,k), min(a4(1,k), pmp_2, lac_2)),  &
                                          max(a4(1,k), pmp_2, lac_2) )
           endif
        endif
     enddo
  else
     do k=3,km-2
        if ( extm(k) .and. (extm(k-1).or.extm(k+1).or.a4(1,k)<qp_min) ) then
             a4(2,k) = a4(1,k)
             a4(3,k) = a4(1,k)
        endif
     enddo
  endif

  do k=1,km-1
     a4(4,k) = 6.*a4(1,k) - 3.*(a4(2,k)+a4(3,k))
  enddo

  k = km-1
  if( extm(k) ) then
      a4(2,k) = a4(1,k)
      a4(3,k) = a4(1,k)
      a4(4,k) = 0.
  else
      da1  = a4(3,k) - a4(2,k)
      da2  = da1**2
      a6da = a4(4,k)*da1
      if(a6da < -da2) then
         a4(4,k) = 3.*(a4(2,k)-a4(1,k))
         a4(3,k) = a4(2,k) - a4(4,k)
      elseif(a6da > da2) then
         a4(4,k) = 3.*(a4(3,k)-a4(1,k))
         a4(2,k) = a4(3,k) - a4(4,k)
      endif
  endif

  call cs_limiters(km-1, a4)

! Bottom layer:
  a4(2,km) = a4(1,km)
  a4(3,km) = a4(1,km)
  a4(4,km) = 0.

 end subroutine cs_profile



 subroutine cs_limiters(km, a4)
 integer, intent(in) :: km
 real, intent(inout) :: a4(4,km)   ! PPM array
! !LOCAL VARIABLES:
 real, parameter:: r12 = 1./12.
 integer k

! Positive definite constraint

 do k=1,km
 if( abs(a4(3,k)-a4(2,k)) < -a4(4,k) ) then
     if( (a4(1,k)+0.25*(a4(3,k)-a4(2,k))**2/a4(4,k)+a4(4,k)*r12) < 0. ) then
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

 end subroutine cs_limiters



 subroutine fall_speed(ktop, kbot, den, qs, qi, qg, ql, tk, vts, vti, vtg)
 integer, intent(in)                     :: ktop, kbot
 real, intent(in ), dimension(ktop:kbot) :: den, qs, qi, qg, ql, tk
 real, intent(out), dimension(ktop:kbot) :: vts, vti, vtg
! fall velocity constants:
 real, parameter:: thi = 1.0e-8 ! cloud ice threshold for terminal fall
 real, parameter:: thg = 1.0e-8
 real, parameter:: ths = 1.0e-8
 real, parameter:: aa = -4.14122e-5, bb = -0.00538922, cc = -0.0516344, dd = 0.00216078, ee = 1.9714
! marshall-palmer constants
 real, parameter:: vcons = 6.6280504, vcong = 87.2382675
 real, parameter:: norms = 942477796.076938, normg =  5026548245.74367
 real, dimension(ktop:kbot) :: qden, tc, rhof
 real :: vi0
 integer:: k
!-----------------------------------------------------------------------
! marshall-palmer formula
!-----------------------------------------------------------------------

! try the local air density -- for global model; the true value could be
! much smaller than sfcrho over high mountains

  do k=ktop, kbot
     rhof(k) = sqrt( min(10.,  sfcrho/den(k)) )
  enddo

! ice:
   if ( const_vi ) then
       vti(:) = vi_fac
   else
! use Deng and Mace (2008, GRL), which gives smaller fall speed than HD90 formula
       vi0 = 0.01*vi_fac
       do k=ktop, kbot
          if ( qi(k) < thi ) then  ! this is needed as the fall-speed maybe problematic for small qi
               vti(k) = vf_min
          else
             tc(k) = tk(k) - tice
            vti(k) = (3.+log10(qi(k)*den(k)))*(tc(k)*(aa*tc(k)+bb)+cc)+dd*tc(k)+ee
            vti(k) = vi0*exp(log_10*vti(k))
            vti(k) = min(vi_max, max(vf_min, vti(k)))
          endif
       enddo
   endif

   if ( const_vs ) then
      vts(:) = vs_fac     ! 1.   IFS_2016
   else
      do k=ktop, kbot
! snow:
         if ( qs(k) < ths ) then
           vts(k) = vf_min
         else
           vts(k) = vs_fac*vcons*rhof(k)*exp(0.0625*log(qs(k)*den(k)/norms))
           vts(k) = min(vs_max, max(vf_min, vts(k)))
         endif
      enddo
   endif

! graupel:
   if ( const_vg ) then
      vtg(:) = vg_fac     ! 2.
   else
     do k=ktop, kbot
        if ( qg(k) < thg ) then
             vtg(k) = vf_min
        else
             vtg(k) = vg_fac*vcong*rhof(k)*sqrt(sqrt(sqrt(qg(k)*den(k)/normg)))
             vtg(k) = min(vg_max, max(vf_min, vtg(k)))
        endif
     enddo
   endif

 end subroutine fall_speed


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
!
 real :: rnzr, rnzs, rnzg, rhos, rhog
 !Intercept parameters
 data rnzr /8.0e6/  ! lin83
 data rnzs /3.0e6/  ! lin83
 data rnzg /4.0e6/  ! rh84
 !Density parameters
 data rhos /0.1e3/  ! lin83    (snow density; 1/10 of water)
 data rhog /0.4e3/  ! rh84     (graupel density)
 data acc/5.0,2.0,0.5/

 real den_rc
 integer :: k, i

      pie = 4.*atan(1.0)

! S. Klein's formular (EQ 16) from AM2
      fac_rc = (4./3.)*pie*rhor*rthresh**3

   if ( prog_ccn ) then
!      if(master) write(*,*) 'prog_ccn option is .T.'
   else
      den_rc = fac_rc * ccn_o*1.e6
!      if(master) write(*,*) 'MP: rthresh=', rthresh, 'vi_fac=', vi_fac
!      if(master) write(*,*) 'MP: for ccn_o=', ccn_o, 'ql_rc=', den_rc
      den_rc = fac_rc * ccn_l*1.e6
!      if(master) write(*,*) 'MP: for ccn_l=', ccn_l, 'ql_rc=', den_rc
   endif

      vdifu=2.11e-5
      tcond=2.36e-2

      visk=1.259e-5
      hlts=2.8336e6
      hltc=2.5e6
      hltf=3.336e5

      ch2o=4.1855e3
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
! Decreasing  csacw to reduce cloud water ---> snow

      craci = pie*rnzr*alin*gam380/(4.*act(2)**0.95)
      csaci = csacw * c_psaci
!
      cgacw = pie*rnzg*gam350*gcon/(4.*act(6)**0.875)
!     cgaci = cgacw*0.1
! SJL, May 28, 2012
      cgaci = cgacw*0.05
!
      cracw = craci            ! cracw= 3.27206196043822
      cracw = c_cracw * cracw
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

 end subroutine setupm


! subroutine lin_cld_microphys_init(id, jd, kd, axes, time)
 subroutine lin_cld_microphys_init

!    integer,         intent(in) :: id, jd, kd
!    integer,         intent(in) :: axes(4)
!    type(time_type), intent(in) :: time

    integer   :: unit, io, ierr, k, logunit
    logical   :: flag
    real :: tmp, q1, q2

    log_10 = log(10.)

!    master = (mpp_pe().eq.mpp_root_pe())

!#ifdef INTERNAL_FILE_NML
!    read( input_nml_file, nml = lin_cld_microphys_nml, iostat = io )
!    ierr = check_nml_error(io,'lin_cloud_microphys_nml')
!#else
!    if( file_exist( 'input.nml' ) ) then
!       unit = open_namelist_file ()
!       io = 1
!       do while ( io .ne. 0 )
!          read( unit, nml = lin_cld_microphys_nml, iostat = io, end = 10 )
!          ierr = check_nml_error(io,'lin_cloud_microphys_nml')
!       end do
!10     call close_file ( unit )
!    end if
!#endif
!    call write_version_number (version, tagname)
!    logunit = stdlog()

    if ( do_setup ) then
      call setup_con
      call setupm
      do_setup = .false.
    endif

    tice0 = tice - 0.01
    t_wfr = tice - 40. ! supercooled water can exist down to -48 C, which is the "absolute"

!    if (master) write( logunit, nml = lin_cld_microphys_nml )

!    id_vtr = register_diag_field ( mod_name, 'vt_r', axes(1:3), time,        &
!         'rain fall speed', 'm/sec', missing_value=missing_value )
    id_vtr = 1
!    id_vts = register_diag_field ( mod_name, 'vt_s', axes(1:3), time,        &
!         'snow fall speed', 'm/sec', missing_value=missing_value )
    id_vts = 1
!    id_vtg = register_diag_field ( mod_name, 'vt_g', axes(1:3), time,        &
!         'graupel fall speed', 'm/sec', missing_value=missing_value )
    id_vtg = 1
!    id_vti = register_diag_field ( mod_name, 'vt_i', axes(1:3), time,        &
!         'ice fall speed', 'm/sec', missing_value=missing_value )
    id_vti = 1
!    id_droplets = register_diag_field ( mod_name, 'droplets', axes(1:3), time,        &
!         'droplet number concentration', '#/m3', missing_value=missing_value )
    id_droplets = 1
!    id_rh = register_diag_field ( mod_name, 'rh_lin', axes(1:2), time,        &
!         'relative humidity', 'n/a', missing_value=missing_value )
    id_rh = 1

!    id_rain = register_diag_field ( mod_name, 'rain_lin', axes(1:2), time,        &
!         'rain_lin', 'mm/day', missing_value=missing_value )
    id_rain = 1
!    id_snow = register_diag_field ( mod_name, 'snow_lin', axes(1:2), time,        &
!         'snow_lin', 'mm/day', missing_value=missing_value )
    id_snow = 1
!    id_graupel = register_diag_field ( mod_name, 'graupel_lin', axes(1:2), time,  &
!         'graupel_lin', 'mm/day', missing_value=missing_value )
    id_graupel = 1
!    id_ice = register_diag_field ( mod_name, 'ice_lin', axes(1:2), time,        &
!         'ice_lin', 'mm/day', missing_value=missing_value )
    id_ice = 1
!    id_prec = register_diag_field ( mod_name, 'prec_lin', axes(1:2), time,     &
!         'prec_lin', 'mm/day', missing_value=missing_value )
    id_prec = 1
!   if ( master ) write(*,*) 'prec_lin diagnostics initialized.', id_prec

!    id_cond = register_diag_field ( mod_name, 'cond_lin', axes(1:2), time,     &
!         'total condensate', 'kg/m**2', missing_value=missing_value )
    id_cond = 1

!    id_var = register_diag_field ( mod_name, 'var_lin', axes(1:2), time,     &
!         'subgrid variance', 'n/a',  missing_value=missing_value )
    id_var = 1

!   call qsmith_init

! TESTING the water vapor tables
!   if ( mp_debug .and. master ) then
!        write(*,*) 'TESTING water vapor tables in lin_cld_microphys'
!        tmp = tice - 90.
!   do k=1,25
!      q1 = wqsat_moist(tmp, 0., 1.E5)
!      q2 = qs1d_m(tmp, 0., 1.E5)
!      write(*,*) NINT(tmp-tice), q1, q2, 'dq=', q1-q2
!      tmp = tmp + 5.
!   enddo
!   endif

!   if ( master ) write(*,*) 'lin_cld_micrphys diagnostics initialized.'

!   lin_cld_mp_clock = mpp_clock_id('Lin_cld_microphys', grain=CLOCK_ROUTINE)

   module_is_initialized = .true.

 end subroutine lin_cld_microphys_init



 subroutine lin_cld_microphys_end

   deallocate ( table  )
   deallocate ( table2 )
   deallocate ( table3 )
   deallocate ( tablew )
   deallocate ( des )
   deallocate ( des2 )
   deallocate ( des3 )
   deallocate ( desw )

   tables_are_initialized = .false.

 end subroutine lin_cld_microphys_end



 subroutine setup_con

!  master = (mpp_pe().eq.mpp_root_pe())
  rgrav = 1./ grav

  if ( .not. qsmith_tables_initialized ) call qsmith_init
  qsmith_tables_initialized = .true.

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

  if( .not. tables_are_initialized ) then

!    master = (mpp_pe().eq.mpp_root_pe())
!    if (master) print*, ' lin MP: initializing qs tables'
!!! DEBUG CODE
!    print*, mpp_pe(), allocated(table), allocated(table2), allocated(table3), allocated(tablew), allocated(des), allocated(des2), allocated(des3), allocated(desw)
!!! END DEBUG CODE

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

       tables_are_initialized = .true.
  endif

 end subroutine qsmith_init

 real function wqs1(ta, den)
! Pure water phase; universal dry/moist formular using air density
! Input "den" can be either dry or moist air density
  real, intent(in):: ta, den
! local:
  real es, ap1
  real, parameter:: tmin=table_ice - 160.
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = tablew(it) + (ap1-it)*desw(it)
      wqs1 = es / (rvgas*ta*den)

 end function wqs1

 real function wqs2(ta, den, dqdt)
! Pure water phase; universal dry/moist formular using air density
! Input "den" can be either dry or moist air density
  real, intent(in):: ta, den
  real, intent(out):: dqdt
! local:
  real es, ap1
  real, parameter:: tmin=table_ice - 160.
  integer it

  if (.not. tables_are_initialized) call qsmith_init

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = tablew(it) + (ap1-it)*desw(it)
      wqs2 = es / (rvgas*ta*den)
        it = ap1 - 0.5
! Finite diff, del_T = 0.1:
      dqdt = 10.*(desw(it) + (ap1-it)*(desw(it+1)-desw(it))) / (rvgas*ta*den)

 end function wqs2

 real function wet_bulb(q, t, den)
! Liquid phase only
 real, intent(in):: t, q, den
 real:: qs, tp, dqdt

  wet_bulb = t
  qs = wqs2(wet_bulb, den, dqdt)
  tp = 0.5*(qs-q)/(1.+lcp*dqdt)*lcp
  wet_bulb = wet_bulb - tp 
! tp is negative if super-saturated
  if ( tp > 0.01 ) then
       qs = wqs2(wet_bulb, den, dqdt)
       tp = (qs-q)/(1.+lcp*dqdt)*lcp
       wet_bulb = wet_bulb - tp 
  endif

 end function wet_bulb
 real function iqs1(ta, den)
! water-ice phase; universal dry/moist formular using air density
! Input "den" can be either dry or moist air density
  real, intent(in):: ta, den
! local:
  real es, ap1
  real, parameter:: tmin=table_ice - 160.
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = table2(it) + (ap1-it)*des2(it)
      iqs1 = es / (rvgas*ta*den)

 end function iqs1

 real function iqs2(ta, den, dqdt)
! water-ice phase; universal dry/moist formular using air density
! Input "den" can be either dry or moist air density
  real, intent(in):: ta, den
  real, intent(out):: dqdt
! local:
  real es, ap1
  real, parameter:: tmin=table_ice - 160.
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = table2(it) + (ap1-it)*des2(it)
      iqs2 = es / (rvgas*ta*den)
        it = ap1 - 0.5
      dqdt = 10.*(des2(it) + (ap1-it)*(des2(it+1)-des2(it))) / (rvgas*ta*den)

 end function iqs2

 real function qs1d_moist(ta, qv, pa, dqdt)
! 2-phase tabel
  real, intent(in):: ta, pa, qv
  real, intent(out):: dqdt
! local:
  real es, ap1
  real, parameter:: tmin=table_ice - 160.
  real, parameter:: eps10 = 10.*eps
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = table2(it) + (ap1-it)*des2(it)
      qs1d_moist = eps*es*(1.+zvir*qv)/pa
        it = ap1 - 0.5
      dqdt = eps10*(des2(it) + (ap1-it)*(des2(it+1)-des2(it)))*(1.+zvir*qv)/pa

 end function qs1d_moist

 real function wqsat2_moist(ta, qv, pa, dqdt)
! Pure water phase
  real, intent(in):: ta, pa, qv
  real, intent(out):: dqdt
! local:
  real es, ap1
  real, parameter:: tmin=table_ice - 160.
  real, parameter:: eps10 = 10.*eps
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = tablew(it) + (ap1-it)*desw(it)
     wqsat2_moist = eps*es*(1.+zvir*qv)/pa
     dqdt = eps10*(desw(it) + (ap1-it)*(desw(it+1)-desw(it)))*(1.+zvir*qv)/pa

 end function wqsat2_moist

 real function wqsat_moist(ta, qv, pa)
! Pure water phase
  real, intent(in):: ta, pa, qv
! local:
  real es, ap1
  real, parameter:: tmin=table_ice - 160.
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = tablew(it) + (ap1-it)*desw(it)
     wqsat_moist = eps*es*(1.+zvir*qv)/pa

 end function wqsat_moist

 real function qs1d_m(ta, qv, pa)
! 2-phase tabel
  real, intent(in):: ta, pa, qv
! local:
  real es, ap1
  real, parameter:: tmin=table_ice - 160.
  real, parameter:: eps10 = 10.*eps
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = table2(it) + (ap1-it)*des2(it)
      qs1d_m = eps*es*(1.+zvir*qv)/pa

 end function qs1d_m

 real function d_sat(ta)
! Computes the difference in saturation vapor *density* between water and ice
  real, intent(in):: ta
  real, parameter:: tmin=table_ice - 160.
  real es_w, es_i, ap1
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
! over Water:
       es_w = tablew(it) + (ap1-it)*desw(it)
! over Ice:
       es_i = table2(it) + (ap1-it)*des2(it)
      d_sat = dim(es_w, es_i)/(rvgas*ta)  ! Take positive difference

 end function d_sat


 real function esw_table(ta)
! pure water phase table
  real, intent(in):: ta
  real, parameter:: tmin=table_ice - 160.
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
  real, parameter:: tmin=table_ice - 160.
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
  real, parameter:: tmin=table_ice - 160.
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
  real, parameter:: tmin=table_ice - 160.
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
  real, parameter:: tmin=table_ice - 160.
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
! Over water
      integer, intent(in):: n
      real:: delt=0.1
      real esbasw, tbasw, esbasi, tbasi, tmin, tem
      integer i

! constants
      esbasw = 1013246.0
       tbasw = table_ice + 100.  !   373.16
      esbasi =    6107.1
       tbasi = table_ice
        tmin = tbasi - 160.

     do i=1,n
        tem = tmin+delt*real(i-1)
!  compute es over water
        tablew(i) = e00*exp((dc_vap*log(tem/t_ice)+Lv0*(tem-t_ice)/(tem*t_ice))/rvgas)
     enddo

 end subroutine qs_tablew


 subroutine qs_table2(n)
! 2-phase table
  integer, intent(in):: n
  real:: delt=0.1
  real esbasw, tbasw, esbasi, tbasi, tmin, tem
  integer :: i0, i1
  real :: tem0, tem1
  integer i

! constants
      esbasw = 1013246.0
       tbasw = table_ice + 100. !    373.16
      esbasi =    6107.1
       tbasi = table_ice
      tmin = tbasi - 160.

     do i=1,n
        tem = tmin+delt*real(i-1)
        if ( i<= 1600 ) then
!  compute es over ice between -160c and 0 c.
             table2(i) = e00*exp((d2ice*log(tem/t_ice)+Li2*(tem-t_ice)/(tem*t_ice))/rvgas)
        else
!  compute es over water between 0c and 102c.
             table2(i) = e00*exp((dc_vap*log(tem/t_ice)+Lv0*(tem-t_ice)/(tem*t_ice))/rvgas)
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
       tbasw = table_ice + 100. !     373.16
      esbasi =    6107.1
       tbasi = table_ice !     273.16
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


 real function qs_blend(t, p, q)
! Note: this routine is based on "moist" mixing ratio
! Blended mixed phase table
  real, intent(in):: t, p, q
  real es, ap1
  real, parameter:: tmin=table_ice - 160.
  integer it

       ap1 = 10.*dim(t, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = table(it) + (ap1-it)*des(it)
      qs_blend = eps*es*(1.+zvir*q)/p

 end function qs_blend

 subroutine qs_table(n)
      integer, intent(in):: n
      real esupc(200)
      real:: delt=0.1
      real esbasw, tbasw, esbasi, tbasi, tmin, tem, esh20
      real wice, wh2o
      integer i

! constants
      esbasw = 1013246.0
       tbasw = table_ice + 100.  !  373.16
      esbasi =    6107.1
       tbasi = table_ice !     273.16

!  compute es over ice between -160c and 0 c.
      tmin = tbasi - 160.
!  see smithsonian meteorological tables page 350.
      do i=1,1600
         tem = tmin+delt*real(i-1)
         table(i) = e00*exp((d2ice*log(tem/t_ice)+Li2*(tem-t_ice)/(tem*t_ice))/rvgas)
      enddo

!  compute es over water between -20c and 102c.
!  see smithsonian meteorological tables page 350.
      do  i=1,1221
          tem = 253.16+delt*real(i-1)
          esh20 = e00*exp((dc_vap*log(tem/t_ice)+Lv0*(tem-t_ice)/(tem*t_ice))/rvgas)
          if (i <= 200) then
              esupc(i) = esh20
          else
              table(i+1400) = esh20
          endif
      enddo

!  derive blended es over ice and supercooled water between -20c and 0c
      do i=1,200
         tem  = 253.16+delt*real(i-1)
         wice = 0.05*(table_ice-tem)
         wh2o = 0.05*(tem-253.16)
         table(i+1400) = wice*table(i+1400)+wh2o*esupc(i)
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
  real, parameter:: tmin=table_ice - 160.
  integer i, k, it

  if( .not. tables_are_initialized ) then
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


 subroutine neg_adj(ktop, kbot, pt, dp, qv, ql, qr, qi, qs, qg)
! 1d version:
! this is designed for 6-class micro-physics schemes
 integer, intent(in):: ktop, kbot
 real, intent(in):: dp(ktop:kbot)
 real, intent(inout), dimension(ktop:kbot)::    &
                                pt, qv, ql, qr, qi, qs, qg
! local:
 real lcpk(ktop:kbot), icpk(ktop:kbot)
 real dq, tmp1, cvm
 integer k

 do k=ktop,kbot
    cvm = c_air + qv(k)*c_vap + (qr(k)+ql(k))*c_liq + (qi(k)+qs(k)+qg(k))*c_ice
    lcpk(k) = (lv00+dc_vap*pt(k)) / cvm
    icpk(k) = (li00+dc_ice*pt(k)) / cvm
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


! real function g_sum(p, ifirst, ilast, jfirst, jlast, area, mode)
!!-------------------------
!! Quick local sum algorithm
!!-------------------------
! use mpp_mod,           only: mpp_sum
! integer, intent(IN) :: ifirst, ilast
! integer, intent(IN) :: jfirst, jlast
! integer, intent(IN) :: mode  ! if ==1 divided by area
! real, intent(IN) :: p(ifirst:ilast,jfirst:jlast)      ! field to be summed
! real, intent(IN) :: area(ifirst:ilast,jfirst:jlast)
! integer :: i,j
! real gsum
!
!   if( global_area < 0. ) then
!       global_area = 0.
!       do j=jfirst,jlast
!          do i=ifirst,ilast
!             global_area = global_area + area(i,j)
!          enddo
!       enddo
!       call mpp_sum(global_area)
!   end if
!
!   gsum = 0.
!   do j=jfirst,jlast
!      do i=ifirst,ilast
!         gsum = gsum + p(i,j)*area(i,j)
!      enddo
!   enddo
!   call mpp_sum(gsum)
!
!   if ( mode==1 ) then
!        g_sum = gsum / global_area
!   else
!        g_sum = gsum
!   endif
!
! end function g_sum

 subroutine interpolate_z(is, ie, js, je, km, zl, hght, a3, a2)

 integer,  intent(in):: is, ie, js, je, km
 real, intent(in):: hght(is:ie,js:je,km+1)  ! hght(k) > hght(k+1)
 real, intent(in):: a3(is:ie,js:je,km)
 real, intent(in):: zl
 real, intent(out):: a2(is:ie,js:je)
! local:
 real zm(km)
 integer i,j,k


!$OMP parallel do default(none) shared(is,ie,js,je,km,hght,zl,a2,a3) private(zm)
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

 subroutine cloud_diagnosis(is, ie, js, je, den, qw, qi, qr, qs, qg, T, qcw, qci, qcr, qcs, qcg, rew, rei, rer, res, reg)

   implicit none

   integer, intent(in) :: is, ie, js, je
   real, dimension(is:ie,js:je), intent(in) :: den, T
   real, dimension(is:ie,js:je), intent(in) :: qw, qi, qr, qs, qg        ! units: kg/kg
   real, dimension(is:ie,js:je), intent(out) :: qcw, qci, qcr, qcs, qcg  ! units: kg/m^3
   real, dimension(is:ie,js:je), intent(out) :: rew, rei, rer, res, reg  ! units: micron

   integer :: i, j
   real :: lambdar, lambdas, lambdag

   real :: rhow = 1.0E3, rhor = 1.0E3, rhos = 1.0E2, rhog = 4.0E2
   real :: n0r = 8.0E6, n0s = 3.0E6, n0g = 4.0E6
   real :: alphar = 0.8, alphas = 0.25, alphag = 0.5
   real :: gammar = 17.837789, gammas = 8.2850630, gammag = 11.631769
   real :: qmin = 1.0E-5, ccn = 1.0E8, beta = 1.22

!   real :: rewmin = 1.0, rewmax = 25.0
!   real :: reimin = 10.0, reimax = 300.0
!   real :: rermin = 25.0, rermax = 225.0
!   real :: resmin = 300, resmax = 1000.0
!   real :: regmin = 1000.0, regmax = 1.0E5
   real :: rewmin = 5.0, rewmax = 10.0
   real :: reimin = 10.0, reimax = 150.0
   real :: rermin = 0.0, rermax = 10000.0
   real :: resmin = 0.0, resmax = 10000.0
   real :: regmin = 0.0, regmax = 10000.0

   do j = js, je
     do i = is, ie

! cloud water (Martin et al., 1994)
       if (qw(i,j) .gt. qmin) then
         qcw(i,j) = den(i,j) * qw(i,j)
         rew(i,j) = exp(1.0 / 3.0 * log((3 * qcw(i,j)) / (4 * pi * rhow * ccn))) * 1.0E6
         rew(i,j) = max(rewmin, min(rewmax, rew(i,j)))
       else
         qcw(i,j) = 0.0
         rew(i,j) = rewmin
       end if

! cloud ice (Heymsfield and McFarquhar, 1996)
       if (qi(i,j) .gt. qmin) then
         qci(i,j) = den(i,j) * qi(i,j)
         if (T(i,j) - tice .lt. -50) then
           rei(i,j) = beta / 9.917 * exp((1 - 0.891) * log(1.0E3 * qci(i,j))) * 1.0E3
         elseif (T(i,j) - tice .lt. -40) then
           rei(i,j) = beta / 9.337 * exp((1 - 0.920) * log(1.0E3 * qci(i,j))) * 1.0E3
         elseif (T(i,j) - tice .lt. -30) then
           rei(i,j) = beta / 9.208 * exp((1 - 0.945) * log(1.0E3 * qci(i,j))) * 1.0E3
         else
           rei(i,j) = beta / 9.387 * exp((1 - 0.969) * log(1.0E3 * qci(i,j))) * 1.0E3
         end if
         rei(i,j) = max(reimin, min(reimax, rei(i,j)))
       else
         qci(i,j) = 0.0
         rei(i,j) = reimin
       end if

! rain (Lin et al., 1983)
       if (qr(i,j) .gt. qmin) then
         qcr(i,j) = den(i,j) * qr(i,j)
         lambdar = exp(0.25 * log(pi * rhor * n0r / qcr(i,j)))
         rer(i,j) = 0.5 * exp(log(gammar / 6) / alphar) / lambdar * 1.0E6
         rer(i,j) = max(rermin, min(rermax, rer(i,j)))
       else
         qcr(i,j) = 0.0
         rer(i,j) = rermin
       end if

! snow (Lin et al., 1983)
       if (qs(i,j) .gt. qmin) then
         qcs(i,j) = den(i,j) * qs(i,j)
         lambdas = exp(0.25 * log(pi * rhos * n0s / qcs(i,j)))
         res(i,j) = 0.5 * exp(log(gammas / 6) / alphas) / lambdas * 1.0E6
         res(i,j) = max(resmin, min(resmax, res(i,j)))
       else
         qcs(i,j) = 0.0
         res(i,j) = resmin
       end if

! graupel (Lin et al., 1983)
       if (qg(i,j) .gt. qmin) then
         qcg(i,j) = den(i,j) * qg(i,j)
         lambdag = exp(0.25 * log(pi * rhog * n0g / qcg(i,j)))
         reg(i,j) = 0.5 * exp(log(gammag / 6) / alphag) / lambdag * 1.0E6
         reg(i,j) = max(regmin, min(regmax, reg(i,j)))
       else
         qcg(i,j) = 0.0
         reg(i,j) = regmin
       end if

     end do
   end do

 end subroutine cloud_diagnosis

end module lin_cld_microphys_mod
