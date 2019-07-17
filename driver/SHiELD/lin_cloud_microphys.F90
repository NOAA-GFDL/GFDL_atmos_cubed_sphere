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
 public  cracw
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
! real, parameter:: e00 = 610.71  ! saturation vapor pressure at T0
   real, parameter:: e00 = 611.21  ! IFS: saturation vapor pressure at T0
   real, parameter:: c_liq = 4.1855e+3       ! heat capacity of water at 0C
!  real, parameter:: c_liq = 4218.        ! ECMWF-IFS
!  real, parameter:: c_ice = 2106.           ! heat capacity of ice at 0C: c=c_ice+7.3*(T-Tice) 
   real, parameter:: c_ice = 1972.           ! heat capacity of ice at -15 C
   real, parameter:: cp_vap = cp_vapor   ! 1846.
!  real, parameter:: cv_vap = 1410.0     ! Emanuel value
! For consistency, cv_vap derived FMS constants:
   real, parameter:: cv_vap = cp_vap - rvgas  ! 1384.5
   real, parameter:: dc_ice =  c_liq - c_ice     ! =  2084
   real, parameter:: dc_vap = cp_vap - c_liq     ! = -2344.    isobaric heating/cooling
! Values at 0 Deg C
! GFS value
   real, parameter:: hlv0 = 2.5e6
!  real, parameter:: hlv0 = 2.501e6   ! Emanuel Appendix-2
! GFS value
   real, parameter:: hlf0 = 3.3358e5
!  real, parameter:: hlf0 = 3.337e5   ! Emanuel
   real, parameter:: t_ice = 273.16
! Latent heat at absolute zero:
   real, parameter:: li00 = hlf0 - dc_ice*t_ice   ! = -2.355446e5
   real, parameter:: Lv0 =  hlv0 - dc_vap*t_ice   ! = 3.141264e6

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
 real :: lv00, d0_vap, c_air, c_vap

 logical :: de_ice = .false.     !
 logical :: sedi_transport = .true.     !
 logical :: do_sedi_w = .false.
 logical :: do_sedi_heat = .true.      !
 logical :: prog_ccn = .false.     ! do prognostic CCN (Yi Ming's method)
 logical :: do_qa   = .true.       ! do inline cloud fraction
 logical :: rad_snow =.true.
 logical :: rad_graupel =.true.
 logical :: rad_rain =.true.
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
 real :: tau_r   = 900.    ! rain freezing time scale during fast_sat
 real :: tau_s   = 900.    ! snow melt
 real :: tau_g   = 600.    ! graupel melt
 real :: tau_mlt = 600.    ! ice melting time-scale

! Fast MP:
 real :: tau_i2s = 1000.  ! ice2snow auto-conversion time scale (sec) 
 real :: tau_l2r = 900.
! cloud water
 real :: tau_v2l = 150.   ! vapor --> cloud water (condensation)  time scale
 real :: tau_l2v = 300.   ! cloud water --> vapor (evaporation)  time scale
! Graupel
 real :: tau_g2v = 900.   ! Grapuel sublimation time scale
 real :: tau_v2g = 21600. ! Grapuel deposition -- make it a slow process

 real :: dw_land  = 0.20  ! base value for subgrid deviation/variability over land
 real :: dw_ocean = 0.10  ! base value for ocean
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
 real :: qi_lim  = 1.
 real :: ql_mlt  = 2.0e-3    ! max value of cloud water allowed from melted cloud ice
 real :: ql_gen  = 1.0e-3    ! max ql generation during remapping step if fast_sat_adj = .T.
 real :: sat_adj0 = 0.90     ! adjustment factor (0: no, 1: full) during fast_sat_adj

! Cloud condensate upper bounds: "safety valves" for ql & qi
 real :: ql0_max = 2.0e-3    ! max ql value (auto converted to rain)
 real :: qi0_max = 1.0e-4    ! max qi value (by other sources)

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

 public mp_time, t_min, t_sub, tau_r, tau_s, tau_g, dw_land, dw_ocean,  &
        vi_fac, vr_fac, vs_fac, vg_fac, ql_mlt, do_qa, fix_negative, &
        vi_max, vs_max, vg_max, vr_max,        &
        qs0_crt, qi_gen, ql0_max, qi0_max, qi0_crt, qr0_crt, fast_sat_adj, &
        rh_inc, rh_ins, rh_inr, const_vi, const_vs, const_vg, const_vr,    &
        use_ccn, rthresh, ccn_l, ccn_o, qc_crt, tau_g2v, tau_v2g, sat_adj0,    &
        c_piacr, tau_mlt, tau_v2l, tau_l2v, tau_i2s, tau_l2r, qi_lim, ql_gen,  &
        c_paut, c_psaci, c_pgacs, z_slope_liq, z_slope_ice, prog_ccn,  &
        c_cracw, alin, clin, tice, rad_snow, rad_graupel, rad_rain,   &
        cld_min, use_ppm, mono_prof, do_sedi_heat, sedi_transport,   &
        do_sedi_w, de_ice, mp_print

!---- version number -----
 character(len=128) :: version = '$Id: lin_cloud_microphys.F90,v 21.0.2.1 2014/12/18 21:14:54 Lucas.Harris Exp $'
 character(len=128) :: tagname = '$Name:  $'

 contains


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


  end subroutine lin_cld_microphys_driver


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


 subroutine lin_cld_microphys_init

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
    lcpk(k) = (lv00+d0_vap*pt(k)) / cvm
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
