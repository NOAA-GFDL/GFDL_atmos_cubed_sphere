!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
! =======================================================================
! cloud micro - physics package for gfdl global cloud resolving model
! the algorithms are originally derived from lin et al 1983. most of the
! key elements have been simplified / improved. this code at this stage
! bears little to no similarity to the original lin mp in zetac.
! therefore, it is best to be called gfdl micro - physics (gfdl mp) .
! developer: shian - jiann lin, linjiong zhou
! revision: inline gfdl cloud microphysics, 9 / 8 / 2017
! =======================================================================

module gfdl_mp_mod

    ! use mpp_mod, only: stdlog, mpp_pe, mpp_root_pe, mpp_clock_id, &
    ! mpp_clock_begin, mpp_clock_end, clock_routine, &
    ! input_nml_file
    ! use time_manager_mod, only: time_type
    ! use constants_mod, only: grav, rdgas, rvgas, cp_air, hlv, hlf, pi => pi_8
    ! use fms_mod, only: write_version_number, open_namelist_file, &
    ! check_nml_error, file_exist, close_file
    use fv_arrays_mod, only: r_grid

    implicit none

    private

    public gfdl_mp_driver, gfdl_mp_init, gfdl_mp_end, wqs1, do_hail, wqsat_moist, wqs2, iqs1, iqs2, qsmith_init, c_liq

    real :: missing_value = - 1.e10

    logical :: module_is_initialized = .false.
    logical :: qsmith_tables_initialized = .false.

    character (len = 17) :: mod_name = 'gfdl_mp'

    real, parameter :: grav = 9.80665 ! gfs: acceleration due to gravity
    real, parameter :: rdgas = 287.05 ! gfs: gas constant for dry air
    real, parameter :: rvgas = 461.50 ! gfs: gas constant for water vapor
    real, parameter :: cp_air = 1004.6 ! gfs: heat capacity of dry air at constant pressure
    real, parameter :: hlv = 2.5e6 ! gfs: latent heat of evaporation
    real, parameter :: hlf = 3.3358e5 ! gfs: latent heat of fusion
    real, parameter :: pi = 3.1415926535897931 ! gfs: ratio of circle circumference to diameter
    ! real, parameter :: cp_air = rdgas * 7. / 2. ! 1004.675, heat capacity of dry air at constant pressure
    real, parameter :: cp_vap = 4.0 * rvgas ! 1846.0, heat capacity of water vapore at constnat pressure
    ! real, parameter :: cv_air = 717.56 ! satoh value
    real, parameter :: cv_air = cp_air - rdgas ! 717.55, heat capacity of dry air at constant volume
    ! real, parameter :: cv_vap = 1410.0 ! emanuel value
    real, parameter :: cv_vap = 3.0 * rvgas ! 1384.5, heat capacity of water vapor at constant volume

#ifdef TEST_ICE0
    real, parameter :: c_ice = 1972. ! gfdl: heat capacity of ice at - 15 deg c
    real, parameter :: c_liq = 4.1855e+3 ! gfs: heat capacity of water at 15 c
    ! c_liq - c_ice = 2213
#else
    real, parameter :: c_ice = 2106. ! heat capacity of ice at 0. deg c
    ! ifs documentation:
    real, parameter :: c_liq = 4218. ! c_liq - c_ice = 2112
    ! emanual's book:
    ! real, parameter :: c_liq = 4190.0 ! heat capacity of water at 0 deg c
#endif

    real, parameter :: eps = rdgas / rvgas ! 0.6219934995
    real, parameter :: zvir = rvgas / rdgas - 1. ! 0.6077338443

    real, parameter :: t_ice = 273.16 ! freezing temperature
    real, parameter :: table_ice = 273.16 ! freezing point for qs table

    ! real, parameter :: e00 = 610.71 ! gfdl: saturation vapor pressure at 0 deg c
    real (kind = r_grid), parameter :: e00 = 611.21 ! ifs: saturation vapor pressure at 0 deg c

    real, parameter :: dc_vap = cp_vap - c_liq ! - 2339.5, isobaric heating / cooling
    real, parameter :: dc_ice = c_liq - c_ice ! 2213.5, isobaric heating / colling

    real, parameter :: hlv0 = hlv ! gfs: evaporation latent heat coefficient at 0 deg c
    ! real, parameter :: hlv0 = 2.501e6 ! emanuel appendix - 2
    real, parameter :: hlf0 = hlf ! gfs: fussion latent heat coefficient at 0 deg c
    ! real, parameter :: hlf0 = 3.337e5 ! emanuel

    real, parameter :: lv0 = hlv0 - dc_vap * t_ice! 3.13905782e6, evaporation latent heat coefficient at 0 deg k
    real, parameter :: li0 = hlf0 - dc_ice * t_ice! - 2.7105966e5, fussion latend heat coefficient at 0 deg k

    ! real (kind = r_grid), parameter :: d2ice = dc_vap + dc_ice ! - 126, isobaric heating / cooling
    real (kind = r_grid), parameter :: d2ice = cp_vap - c_ice
    ! d2ice = cp_vap - c_ice
    real (kind = r_grid), parameter :: li2 = lv0 + li0 ! 2.86799816e6, sublimation latent heat coefficient at 0 deg k

    real, parameter :: qrmin = 1.e-8 ! min value for ???
    real, parameter :: qvmin = 1.e-20 ! min value for water vapor (treated as zero)
    real, parameter :: qcmin = 1.e-12 ! min value for cloud condensates

    real, parameter :: vr_min = 1.e-3 ! min fall speed for rain
    real, parameter :: vf_min = 1.e-5 ! min fall speed for cloud ice, snow, graupel

    real, parameter :: dz_min = 1.e-2 ! use for correcting flipped height

    real, parameter :: sfcrho = 1.2 ! surface air density

    ! intercept parameters

    real, parameter :: rnzr = 8.0e6 ! lin83
    real, parameter :: rnzs = 3.0e6 ! lin83
    real, parameter :: rnzg = 4.0e6 ! rh84
    real, parameter :: rnzh = 4.0e4 ! lin83 --- lmh 29 sep 17

    ! density parameters

    real, parameter :: rhor = 1.e3 ! density of rain water, lin83
    real, parameter :: rhos = 0.1e3 ! lin83 (snow density; 1 / 10 of water)
    real, parameter :: rhog = 0.4e3 ! rh84 (graupel density)
    real, parameter :: rhoh = 0.917e3 ! lin83 --- lmh 29 sep 17

    real :: cracs, csacr, cgacr, cgacs, csacw, craci, csaci, cgacw, cgaci, cracw ! constants for accretions
    real :: acco (3, 4) ! constants for accretions
    real :: cssub (5), cgsub (5), crevp (5), cgfr (2), csmlt (5), cgmlt (5)

    real :: es0, ces0
    real :: pie, rgrav, fac_rc
    real :: c_air, c_vap

    real :: lat2, lcp, icp, tcp ! used in bigg mechanism and wet bulk

    real :: d0_vap ! the same as dc_vap, except that cp_vap can be cp_vap or cv_vap
    real (kind = r_grid) :: lv00, li00, li20
    ! scaled constants:
    real (kind = r_grid) :: d1_vap, d1_ice, c1_vap, c1_liq, c1_ice
    real (kind = r_grid), parameter :: one_r8 = 1.

    integer :: ntimes = 1 ! cloud microphysics sub cycles

    ! cloud microphysics switchers

    integer :: icloud_f = 0 ! cloud scheme
    integer :: irain_f = 0 ! cloud water to rain auto conversion scheme

    logical :: de_ice = .false. ! to prevent excessive build - up of cloud ice from external sources
    logical :: sedi_transport = .true. ! transport of momentum in sedimentation
    logical :: do_sedi_w = .true. ! transport of vertical momentum during sedimentation
    logical :: do_sedi_heat = .true. ! transport of heat in sedimentation
    logical :: prog_ccn = .false. ! do prognostic ccn (yi ming's method)
    logical :: do_qa = .true. ! do inline cloud fraction
    logical :: rad_snow = .true. ! consider snow in cloud fraciton calculation
    logical :: rad_graupel = .true. ! consider graupel in cloud fraction calculation
    logical :: rad_rain = .true. ! consider rain in cloud fraction calculation
    logical :: fix_negative = .false. ! fix negative water species
    logical :: do_setup = .true. ! setup constants and parameters
    logical :: disp_heat = .false. ! dissipative heating due to sedimentation
    logical :: do_cond_timescale = .false. ! whether to apply a timescale to condensation

    real, allocatable :: table (:), table2 (:), table3 (:), tablew (:)
    real, allocatable :: des (:), des2 (:), des3 (:), desw (:)

    logical :: tables_are_initialized = .false.

    ! logical :: master
    ! integer :: id_rh, id_vtr, id_vts, id_vtg, id_vti, id_rain, id_snow, id_graupel, &
    ! id_ice, id_prec, id_cond, id_var, id_droplets
    ! integer :: gfdl_mp_clock ! clock for timing of driver routine

    real, parameter :: dt_fr = 8. ! homogeneous freezing of all cloud water at t_wfr - dt_fr
    ! minimum temperature water can exist (moore & molinero nov. 2011, nature)
    ! dt_fr can be considered as the error bar

    real, parameter :: p0_min = 100. ! minimum pressure (pascal) for mp to operate
    real :: p_min

    ! -----------------------------------------------------------------------
    ! namelist parameters
    ! -----------------------------------------------------------------------

    real :: cld_fac = 1.0 ! multiplication factor for cloud fraction
    real :: cld_min = 0.05 ! minimum cloud fraction
    real :: tice = 273.16 ! set tice = 165. to trun off ice - phase phys (kessler emulator)

    ! real :: t_min = 178. ! min temp to freeze - dry all water vapor
    ! sjl 20181123
    real :: t_min = 170. ! min temp to freeze - dry all water vapor
    real :: t_sub = 184. ! min temp for sublimation of cloud ice

    ! relative humidity increment

    real :: rh_inc = 0.25 ! rh increment for complete evaporation of cloud water and cloud ice
    real :: rh_inr = 0.1 ! rh increment for minimum evaporation of rain (not used---originally for "alternative minimum evaporation")
    real :: rh_ins = 0.1 ! rh increment for sublimation of snow (not used)

    ! conversion time scale

    real :: tau_r2g = 900. ! rain freezing during fast_sat
    real :: tau_smlt = 900. ! snow melting
    real :: tau_g2r = 600. ! graupel melting to rain
    real :: tau_imlt = 600. ! cloud ice melting
    real :: tau_i2s = 1000. ! cloud ice to snow auto - conversion
    real :: tau_l2r = 900. ! cloud water to rain auto - conversion
    real :: tau_v2l = 150. ! water vapor to cloud water (condensation)
    real :: tau_l2v = 300. ! cloud water to water vapor (evaporation)
    real :: tau_g2v = 900. ! grapuel sublimation
    real :: tau_v2g = 21600. ! grapuel deposition -- make it a slow process

    ! horizontal subgrid variability

    real :: dw_land = 0.20 ! base value for subgrid deviation / variability over land
    real :: dw_ocean = 0.10 ! base value for ocean

    ! prescribed ccn

    real :: ccn_o = 90. ! ccn over ocean (cm^ - 3)
    real :: ccn_l = 270. ! ccn over land (cm^ - 3)

    real :: rthresh = 10.0e-6 ! critical cloud drop radius (micro m)

    ! -----------------------------------------------------------------------
    ! wrf / wsm6 scheme: qi_gen = 4.92e-11 * (1.e3 * exp (0.1 * tmp)) ** 1.33
    ! optimized: qi_gen = 4.92e-11 * exp (1.33 * log (1.e3 * exp (0.1 * tmp)))
    ! qi_gen ~ 4.808e-7 at 0 c; 1.818e-6 at - 10 c, 9.82679e-5 at - 40c
    ! the following value is constructed such that qc_crt = 0 at zero c and @ - 10c matches
    ! wrf / wsm6 ice initiation scheme; qi_crt = qi_gen * min (qi_lim, 0.1 * tmp) / den
    ! -----------------------------------------------------------------------

    real :: sat_adj0 = 0.90 ! adjustment factor (0: no, 1: full) during fast_sat_adj

    real :: qc_crt = 5.0e-8 ! mini condensate mixing ratio to allow partial cloudiness

    real :: qi_lim = 1. ! cloud ice limiter to prevent large ice build up

    real :: ql_mlt = 2.0e-3 ! max value of cloud water allowed from melted cloud ice
    real :: qs_mlt = 1.0e-6 ! max cloud water due to snow melt

    real :: ql_gen = 1.0e-3 ! max cloud water generation during remapping step if fast_sat_adj = .t.
    real :: qi_gen = 1.82e-6 ! max cloud ice generation during remapping step

    ! cloud condensate upper bounds: "safety valves" for ql & qi

    real :: ql0_max = 2.0e-3 ! max cloud water value (auto converted to rain)
    real :: qi0_max = 1.0e-4 ! max cloud ice value (by other sources)
    real :: qi0_crt = 1.0e-4 ! cloud ice to snow autoconversion threshold (was 1.e-4)
    ! qi0_crt if negative, its magnitude is used as the mixing ration threshold; otherwise, used as density
    real :: qr0_crt = 1.0e-4 ! rain to snow or graupel / hail threshold (not used)
    ! lfo used * mixing ratio * = 1.e-4 (hail in lfo)
    real :: qs0_crt = 1.0e-3 ! snow to graupel density threshold (0.6e-3 in purdue lin scheme)

    real :: c_paut = 0.55 ! autoconversion cloud water to rain (use 0.5 to reduce autoconversion)
    real :: c_psaci = 0.02 ! accretion: cloud ice to snow (was 0.1 in zetac)
    real :: c_piacr = 5.0 ! accretion: rain to ice: (not used)
    real :: c_cracw = 0.9 ! rain accretion efficiency
    real :: c_pgacs = 2.0e-3 ! snow to graupel "accretion" eff. (was 0.1 in zetac)

    ! decreasing clin to reduce csacw (so as to reduce cloud water --- > snow)

    real :: alin = 842.0 ! "a" in lin1983
    real :: clin = 4.8 ! "c" in lin 1983, 4.8 -- > 6. (to ehance ql -- > qs)

    ! fall velocity tuning constants:

    logical :: const_vi = .false. ! if .t. the constants are specified by v * _fac
    logical :: const_vs = .false. ! if .t. the constants are specified by v * _fac
    logical :: const_vg = .false. ! if .t. the constants are specified by v * _fac
    logical :: const_vr = .false. ! if .t. the constants are specified by v * _fac

    ! good values:

    real :: vi_fac = 1. ! if const_vi: 1 / 3
    real :: vs_fac = 1. ! if const_vs: 1.
    real :: vg_fac = 1. ! if const_vg: 2.
    real :: vr_fac = 1. ! if const_vr: 4.

    ! upper bounds of fall speed (with variable speed option)

    real :: vi_max = 0.5 ! max fall speed for ice
    real :: vs_max = 5.0 ! max fall speed for snow
    real :: vg_max = 8.0 ! max fall speed for graupel
    real :: vr_max = 12. ! max fall speed for rain

    ! cloud microphysics switchers

    ! this should be removed with the inline code
    logical :: fast_sat_adj = .false. ! has fast saturation adjustments
    logical :: z_slope_liq = .true. ! use linear mono slope for autocconversions
    logical :: z_slope_ice = .false. ! use linear mono slope for autocconversions
    logical :: use_ccn = .false. ! must be true when prog_ccn is false
    logical :: use_ppm = .false. ! use ppm fall scheme
    logical :: use_ppm_ice = .false. ! use ppm fall scheme for cloud ice
    logical :: mono_prof = .true. ! perform terminal fall with mono ppm scheme
    logical :: mp_print = .false. ! cloud microphysics debugging printout
    logical :: do_hail = .false. ! use hail parameters instead of graupel

    ! real :: global_area = - 1.

    real :: g2, log_10, tice0, t_wfr

    ! -----------------------------------------------------------------------
    ! namelist
    ! -----------------------------------------------------------------------

    namelist / gfdl_mp_nml / &
        t_min, t_sub, tau_r2g, tau_smlt, tau_g2r, dw_land, dw_ocean, &
        vi_fac, vr_fac, vs_fac, vg_fac, ql_mlt, do_qa, fix_negative, vi_max, &
        vs_max, vg_max, vr_max, qs_mlt, qs0_crt, qi_gen, ql0_max, qi0_max, &
        qi0_crt, fast_sat_adj, rh_inc, rh_ins, rh_inr, const_vi, &
        const_vs, const_vg, const_vr, use_ccn, rthresh, ccn_l, ccn_o, qc_crt, &
        tau_g2v, tau_v2g, sat_adj0, tau_imlt, tau_v2l, tau_l2v, &
        tau_i2s, tau_l2r, qi_lim, ql_gen, c_paut, c_psaci, c_pgacs, &
        z_slope_liq, z_slope_ice, prog_ccn, c_cracw, alin, clin, tice, &
        rad_snow, rad_graupel, rad_rain, cld_fac, cld_min, use_ppm, use_ppm_ice, mono_prof, &
        do_sedi_heat, sedi_transport, do_sedi_w, de_ice, icloud_f, irain_f, mp_print, &
        ntimes, disp_heat, do_hail, do_cond_timescale

contains

! -----------------------------------------------------------------------
! the driver of the gfdl cloud microphysics
! -----------------------------------------------------------------------

subroutine gfdl_mp_driver (qv, ql, qr, qi, qs, qg, qa, qn, &
        pt, w, ua, va, dz, delp, gsize, dts, hs, rain, snow, ice, &
        graupel, hydrostatic, phys_hydrostatic, is, ie, ks, ke, q_con, cappa, consv_te, &
        te, last_step)

    implicit none

    logical, intent (in) :: hydrostatic, phys_hydrostatic
    logical, intent (in) :: last_step
    logical, intent (in) :: consv_te

    integer, intent (in) :: is, ie ! physics window
    integer, intent (in) :: ks, ke ! vertical dimension

    real, intent (in) :: dts ! physics time step

    real, intent (in), dimension (is:ie) :: hs, gsize

    real, intent (in), dimension (is:ie, ks:ke) :: dz
    real, intent (in), dimension (is:ie, ks:ke) :: qn

    real, intent (inout), dimension (is:ie, ks:ke) :: delp
    real, intent (inout), dimension (is:ie, ks:ke) :: qv, ql, qr, qi, qs, qg, qa
    real, intent (inout), dimension (is:ie, ks:ke) :: pt, ua, va, w
    real, intent (inout), dimension (is:ie, ks:ke) :: q_con, cappa
    real, intent (inout), dimension (is:ie) :: rain, snow, ice, graupel

    real, intent (out), dimension (is:ie, ks:ke) :: te
    ! logical :: used
    real, dimension (is:ie) :: w_var
    real, dimension (is:ie, ks:ke) :: vt_r, vt_s, vt_g, vt_i, qn2
    real, dimension (is:ie, ks:ke) :: m2_rain, m2_sol

    ! call mpp_clock_begin (gfdl_mp_clock)

    if (last_step) then
        p_min = p0_min ! final clean - up
    else
        p_min = 30.e2 ! time saving trick
    endif

    ! -----------------------------------------------------------------------
    ! define heat capacity of dry air and water vapor based on hydrostatical property
    ! -----------------------------------------------------------------------

    if (hydrostatic .or. phys_hydrostatic) then
        c_air = cp_air
        c_vap = cp_vap
        if (hydrostatic) do_sedi_w = .false.
    else
        c_air = cv_air
        c_vap = cv_vap
    endif
    d0_vap = c_vap - c_liq

    ! scaled constants (to reduce fp errors for 32 - bit) :
    d1_vap = d0_vap / c_air
    d1_ice = dc_ice / c_air

    ! lv0 = hlv0 - (c_vap - c_liq) * t_ice! 3.13905782e6, evaporation latent heat coefficient at 0 deg k
    lv00 = (hlv0 - d0_vap * t_ice) / c_air
    li00 = (hlf0 - dc_ice * t_ice) / c_air
    li20 = lv00 + li00

    c1_vap = c_vap / c_air
    c1_liq = c_liq / c_air
    c1_ice = c_ice / c_air

    ! -----------------------------------------------------------------------
    ! define latent heat coefficient used in wet bulb and bigg mechanism
    ! -----------------------------------------------------------------------

    lat2 = (hlv + hlf) ** 2

    lcp = hlv / cp_air
    icp = hlf / cp_air
    tcp = (hlv + hlf) / cp_air

    ! tendency zero out for am moist processes should be done outside the driver

    ! -----------------------------------------------------------------------
    ! major cloud microphysics
    ! -----------------------------------------------------------------------

    call mpdrv (hydrostatic, ua, va, w, delp, pt, qv, ql, qr, qi, qs, qg, &
        qa, qn, dz, is, ie, ks, ke, dts, &
        rain, snow, graupel, ice, m2_rain, m2_sol, gsize, hs, &
        w_var, vt_r, vt_s, vt_g, vt_i, qn2, q_con, cappa, consv_te, te, &
        last_step)

    ! call mpp_clock_end (gfdl_mp_clock)

end subroutine gfdl_mp_driver

! -----------------------------------------------------------------------
! gfdl cloud microphysics, major program
! lin et al., 1983, jam, 1065 - 1092, and
! rutledge and hobbs, 1984, jas, 2949 - 2972
! terminal fall is handled lagrangianly by conservative fv algorithm
! pt: temperature (k)
! 6 water species:
! 1) qv: water vapor (kg / kg)
! 2) ql: cloud water (kg / kg)
! 3) qr: rain (kg / kg)
! 4) qi: cloud ice (kg / kg)
! 5) qs: snow (kg / kg)
! 6) qg: graupel (kg / kg)
! -----------------------------------------------------------------------

subroutine mpdrv (hydrostatic, ua, va, w, delp, pt, qv, ql, qr, qi, qs, &
        qg, qa, qn, dz, is, ie, ks, ke, dt_in, &
        rain, snow, graupel, ice, m2_rain, m2_sol, gsize, hs, &
        w_var, vt_r, vt_s, vt_g, vt_i, qn2, q_con, cappa, consv_te, te, &
        last_step)

    implicit none

    logical, intent (in) :: hydrostatic
    logical, intent (in) :: last_step
    logical, intent (in) :: consv_te
    integer, intent (in) :: is, ie, ks, ke
    real, intent (in) :: dt_in
    real, intent (in), dimension (is:ie) :: gsize
    real, intent (in), dimension (is:ie) :: hs
    real, intent (in), dimension (is:ie, ks:ke) :: dz
    real, intent (in), dimension (is:ie, ks:ke) :: qn

    real, intent (inout), dimension (is:ie, ks:ke) :: delp
    real, intent (inout), dimension (is:ie, ks:ke) :: qv, ql, qr, qi, qs, qg, qa
    real, intent (inout), dimension (is:ie, ks:ke) :: pt, ua, va, w
    real, intent (inout), dimension (is:ie, ks:ke) :: q_con, cappa
    real, intent (inout), dimension (is:ie) :: rain, snow, ice, graupel

    real, intent (out), dimension (is:ie) :: w_var
    real, intent (out), dimension (is:ie, ks:ke) :: vt_r, vt_s, vt_g, vt_i, qn2
    real, intent (out), dimension (is:ie, ks:ke) :: m2_rain, m2_sol
    real, intent (out), dimension (is:ie, ks:ke) :: te
    ! local:
    real, dimension (ks:ke) :: q_liq, q_sol
    real, dimension (ks:ke) :: qvz, qlz, qrz, qiz, qsz, qgz, qaz
    real, dimension (ks:ke) :: vtiz, vtsz, vtgz, vtrz
    real, dimension (ks:ke) :: dp1, dz1
    real, dimension (ks:ke) :: den, p1, denfac
    real, dimension (ks:ke) :: ccn, c_praut, m1_rain, m1_sol, m1
    real, dimension (ks:ke) :: u0, v0, u1, v1, w1

    real :: cpaut, rh_adj, rh_rain
    real :: r1, s1, i1, g1, rdt, ccn0
    real :: dt_rain
    real :: s_leng, t_land, t_ocean, h_var, tmp
    real (kind = r_grid), dimension (ks:ke) :: dp0, tz, cvm
    real (kind = r_grid) :: con_r8, c8
    real :: convt
    real :: dts

    integer :: i, k, n

    dts = dt_in / real (ntimes)

    dt_rain = dts * 0.5
    rdt = one_r8 / dts

    ! convert to mm / day
    convt = 86400. * rdt * rgrav

    ! -----------------------------------------------------------------------
    ! use local variables
    ! -----------------------------------------------------------------------

    do i = is, ie

        do k = ks, ke
#ifdef MOIST_CAPPA
            tz (k) = pt (i, k) / ((1. + zvir * qv (i, k)) * (1. - (ql (i, k) + qr (i, k) + qi (i, k) + qs (i, k) + qg (i, k))))
#else
            tz (k) = pt (i, k) / (1. + zvir * qv (i, k))
#endif
            dp0 (k) = delp (i, k)
            ! -----------------------------------------------------------------------
            ! convert moist mixing ratios to dry mixing ratios
            ! -----------------------------------------------------------------------
            qvz (k) = qv (i, k)
            qlz (k) = ql (i, k)
            qrz (k) = qr (i, k)
            qiz (k) = qi (i, k)
            qsz (k) = qs (i, k)
            qgz (k) = qg (i, k)
            ! save moist ratios for te:
            q_liq (k) = qlz (k) + qrz (k)
            q_sol (k) = qiz (k) + qsz (k) + qgz (k)
            q_con (i, k) = q_liq (k) + q_sol (k)
            qaz (k) = 0.
            dz1 (k) = dz (i, k)
            con_r8 = one_r8 - (qvz (k) + q_con (i, k))
            ! dp1 is dry mass (no change during mp)
            dp1 (k) = dp0 (k) * con_r8
            con_r8 = one_r8 / con_r8
            qvz (k) = qvz (k) * con_r8
            qlz (k) = qlz (k) * con_r8
            qrz (k) = qrz (k) * con_r8
            qiz (k) = qiz (k) * con_r8
            qsz (k) = qsz (k) * con_r8
            qgz (k) = qgz (k) * con_r8

            den (k) = - dp1 (k) / (grav * dz1 (k)) ! density of dry air
            p1 (k) = den (k) * rdgas * tz (k) ! dry air pressure

            ! -----------------------------------------------------------------------
            ! for sedi_momentum transport:
            ! -----------------------------------------------------------------------

            m1 (k) = 0.
            u0 (k) = ua (i, k)
            v0 (k) = va (i, k)
            w1 (k) = w (i, k)
            u1 (k) = u0 (k)
            v1 (k) = v0 (k)
            denfac (k) = sqrt (sfcrho / den (k))
        enddo

        ! -----------------------------------------------------------------------
        ! fix energy conservation
        ! -----------------------------------------------------------------------

        if (consv_te) then
            if (hydrostatic) then
                do k = ks, ke
                    te (i, k) = - c_air * tz (k) * delp (i, k)
                enddo
            else
                do k = ks, ke
#ifdef MOIST_CAPPA
                    q_liq (k) = ql (i, k) + qr (i, k)
                    q_sol (k) = qi (i, k) + qs (i, k) + qg (i, k)
                    q_con (i, k) = q_liq (k) + q_sol (k)
                    cvm (k) = (one_r8 - (qv (i, k) + q_con (i, k))) * c_air + qv (i, k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
                    te (i, k) = - cvm (k) * tz (k) * delp (i, k)
#else
                    te (i, k) = - c_air * tz (k) * delp (i, k)
#endif
                enddo
            endif
        endif

        ! -----------------------------------------------------------------------
        ! calculate cloud condensation nuclei (ccn)
        ! the following is based on klein eq. 15
        ! -----------------------------------------------------------------------

        cpaut = c_paut * 0.104 * grav / 1.717e-5

        if (prog_ccn) then
            do k = ks, ke
                ! convert # / cc to # / m^3
                ccn (k) = qn (i, k) * 1.e6
                c_praut (k) = cpaut * (ccn (k) * rhor) ** (- 1. / 3.)
            enddo
            use_ccn = .false.
        else
            ccn0 = (ccn_l * min (1., abs (hs (i)) / (10. * grav)) + ccn_o * (1. - min (1., abs (hs (i)) / (10. * grav)))) * 1.e6
            if (use_ccn) then
                ! -----------------------------------------------------------------------
                ! ccn is formulted as ccn = ccn_surface * (den / den_surface)
                ! -----------------------------------------------------------------------
                ccn0 = ccn0 * rdgas * tz (ke) / p1 (ke)
            endif
            tmp = cpaut * (ccn0 * rhor) ** (- 1. / 3.)
            do k = ks, ke
                c_praut (k) = tmp
                ccn (k) = ccn0
            enddo
        endif

        ! -----------------------------------------------------------------------
        ! calculate horizontal subgrid variability
        ! total water subgrid deviation in horizontal direction
        ! default area dependent form: use dx ~ 100 km as the base
        ! -----------------------------------------------------------------------

        s_leng = sqrt (gsize (i) / 1.e5)
        t_land = dw_land * s_leng
        t_ocean = dw_ocean * s_leng
        tmp = min (1., abs (hs (i)) / (10. * grav))
        h_var = t_land * tmp + t_ocean * (1. - tmp)
        h_var = min (0.20, max (0.01, h_var))

        ! -----------------------------------------------------------------------
        ! relative humidity thresholds
        ! -----------------------------------------------------------------------

        rh_adj = 1. - h_var - rh_inc
        rh_rain = max (0.6, rh_adj - rh_inr) ! rh_inr = 0.2

        ! -----------------------------------------------------------------------
        ! fix all negative water species
        ! -----------------------------------------------------------------------

        if (fix_negative) &
            call neg_adj (ks, ke, tz, dp1, qvz, qlz, qrz, qiz, qsz, qgz)

        m2_rain (i, :) = 0.
        m2_sol (i, :) = 0.

        do n = 1, ntimes

            ! -----------------------------------------------------------------------
            ! time - split warm rain processes: 1st pass
            ! -----------------------------------------------------------------------

            call warm_rain (dt_rain, ks, ke, dp1, dz1, tz, qvz, qlz, qrz, qiz, qsz, &
                qgz, den, denfac, ccn, c_praut, rh_rain, vtrz, r1, m1_rain, w1, h_var)

            rain (i) = rain (i) + r1 * convt

            do k = ks, ke
                m2_rain (i, k) = m2_rain (i, k) + m1_rain (k)
                m1 (k) = m1 (k) + m1_rain (k)
            enddo

            ! -----------------------------------------------------------------------
            ! sedimentation of cloud ice, snow, and graupel
            ! -----------------------------------------------------------------------

            call fall_speed (ks, ke, den, qsz, qiz, qgz, qlz, tz, vtsz, vtiz, vtgz)

            call terminal_fall (dts, ks, ke, tz, qvz, qlz, qrz, qgz, qsz, qiz, &
                dz1, dp1, den, vtgz, vtsz, vtiz, r1, g1, s1, i1, m1_sol, w1)

            rain (i) = rain (i) + r1 * convt ! from melted snow & ice that reached the ground
            snow (i) = snow (i) + s1 * convt
            graupel (i) = graupel (i) + g1 * convt
            ice (i) = ice (i) + i1 * convt

            ! -----------------------------------------------------------------------
            ! heat transportation during sedimentation
            ! -----------------------------------------------------------------------

            if (do_sedi_heat) &
                call sedi_heat (ks, ke, dp1, m1_sol, dz1, tz, qvz, qlz, qrz, qiz, &
                qsz, qgz, c_ice)

            ! -----------------------------------------------------------------------
            ! time - split warm rain processes: 2nd pass
            ! -----------------------------------------------------------------------

            call warm_rain (dt_rain, ks, ke, dp1, dz1, tz, qvz, qlz, qrz, qiz, qsz, &
                qgz, den, denfac, ccn, c_praut, rh_rain, vtrz, r1, m1_rain, w1, h_var)

            rain (i) = rain (i) + r1 * convt

            do k = ks, ke
                m2_rain (i, k) = m2_rain (i, k) + m1_rain (k)
                m2_sol (i, k) = m2_sol (i, k) + m1_sol (k)
                m1 (k) = m1 (k) + m1_rain (k) + m1_sol (k)
            enddo

            ! -----------------------------------------------------------------------
            ! ice - phase microphysics
            ! -----------------------------------------------------------------------

            call icloud (ks, ke, tz, p1, qvz, qlz, qrz, qiz, qsz, qgz, dp1, den, &
                denfac, vtsz, vtgz, vtrz, qaz, rh_adj, rh_rain, dts, h_var, last_step)

        enddo

        ! -----------------------------------------------------------------------
        ! momentum transportation during sedimentation
        ! note: dp1 is dry mass; dp0 is the old moist (total) mass
        ! -----------------------------------------------------------------------

        if (sedi_transport) then
            do k = ks + 1, ke
                u1 (k) = (dp0 (k) * u1 (k) + m1 (k - 1) * u1 (k - 1)) / (dp0 (k) + m1 (k - 1))
                v1 (k) = (dp0 (k) * v1 (k) + m1 (k - 1) * v1 (k - 1)) / (dp0 (k) + m1 (k - 1))
                ua (i, k) = u1 (k)
                va (i, k) = v1 (k)
            enddo
            ! sjl modify tz due to ke loss:
            ! seperate loop (vectorize better with no k - dependency)
            if (disp_heat) then
                do k = ks + 1, ke
#ifdef MOIST_CAPPA
                    c8 = c_air + qvz (k) * c_vap + (qrz (k) + qlz (k)) * c_liq + (qiz (k) + qsz (k) + qgz (k)) * c_ice
                    tz (k) = tz (k) + 0.5 * (u0 (k) ** 2 + v0 (k) ** 2 - (u1 (k) ** 2 + v1 (k) ** 2)) / c8
#else
                    tz (k) = tz (k) + 0.5 * (u0 (k) ** 2 + v0 (k) ** 2 - (u1 (k) ** 2 + v1 (k) ** 2)) / c_air
#endif
                enddo
            endif
        endif

        if (do_sedi_w) then
            ! conserve local te
            !#ifdef disp_w
            if (disp_heat) then
                do k = ks, ke
#ifdef MOIST_CAPPA
                    c8 = c_air + qvz (k) * c_vap + (qrz (k) + qlz (k)) * c_liq + (qiz (k) + qsz (k) + qgz (k)) * c_ice
                    tz (k) = tz (k) + 0.5 * (w (i, k) ** 2 - w1 (k) * w1 (k)) / c8
#else
                    tz (k) = tz (k) + 0.5 * (w (i, k) ** 2 - w1 (k) * w1 (k)) / c_air
#endif
                enddo
            endif
            !#endif
            do k = ks, ke
                w (i, k) = w1 (k)
            enddo
        endif

        ! -----------------------------------------------------------------------
        ! update moist air mass (actually hydrostatic pressure)
        ! convert to dry mixing ratios
        ! -----------------------------------------------------------------------

        do k = ks, ke
            ! total mass changed due to sedimentation !!!
            con_r8 = one_r8 + qvz (k) + qlz (k) + qrz (k) + qiz (k) + qsz (k) + qgz (k)
            delp (i, k) = dp1 (k) * con_r8
            ! convert back to moist mixing ratios
            con_r8 = one_r8 / con_r8
            qvz (k) = qvz (k) * con_r8
            qlz (k) = qlz (k) * con_r8
            qrz (k) = qrz (k) * con_r8
            qiz (k) = qiz (k) * con_r8
            qsz (k) = qsz (k) * con_r8
            qgz (k) = qgz (k) * con_r8
            ! all are moist mixing ratios at this point on:
            qv (i, k) = qvz (k)
            ql (i, k) = qlz (k)
            qr (i, k) = qrz (k)
            qi (i, k) = qiz (k)
            qs (i, k) = qsz (k)
            qg (i, k) = qgz (k)
            q_liq (k) = qlz (k) + qrz (k)
            q_sol (k) = qiz (k) + qsz (k) + qgz (k)
            q_con (i, k) = q_liq (k) + q_sol (k)
            cvm (k) = (one_r8 - (qvz (k) + q_con (i, k))) * c_air + qvz (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
            tmp = rdgas * (1. + zvir * qvz (k))
            cappa (i, k) = tmp / (tmp + cvm (k))
#ifdef MOIST_CAPPA
            pt (i, k) = tz (k) * (1. + zvir * qvz (k)) * (1. - q_con (i, k))
#else
            pt (i, k) = tz (k) * (1. + zvir * qvz (k))
#endif
        enddo

        ! -----------------------------------------------------------------------
        ! fix energy conservation
        ! -----------------------------------------------------------------------

        if (consv_te) then
            if (hydrostatic) then
                do k = ks, ke
                    te (i, k) = te (i, k) + c_air * tz (k) * delp (i, k)
                enddo
            else
                do k = ks, ke
#ifdef MOIST_CAPPA
                    te (i, k) = te (i, k) + cvm (k) * tz (k) * delp (i, k)
#else
                    te (i, k) = te (i, k) + c_air * tz (k) * delp (i, k)
#endif
                enddo
            endif
        endif

        ! -----------------------------------------------------------------------
        ! update cloud fraction tendency
        ! -----------------------------------------------------------------------

        do k = ks, ke
            qa (i, k) = qaz (k)
        enddo

    enddo

end subroutine mpdrv

! -----------------------------------------------------------------------
! sedimentation of heat
! -----------------------------------------------------------------------

subroutine sedi_heat (ks, ke, dm, m1, dz, tz, qv, ql, qr, qi, qs, qg, cw)
    ! revised with a precise energy conserving form: s. - j. lin, jan 22, 2018
    ! input q fields are dry mixing ratios, and dm is dry air mass
    implicit none
    integer, intent (in) :: ks, ke
    real, intent (in), dimension (ks:ke) :: dm, m1, dz, qv, ql, qr, qi, qs, qg
    real (kind = r_grid), intent (inout), dimension (ks:ke) :: tz
    real, intent (in) :: cw ! heat capacity
    ! local:
    real, dimension (ks:ke) :: dgz, cv0
    integer :: k

    ! this is the vectorized loop
    do k = ks + 1, ke
        dgz (k) = - g2 * (dz (k - 1) + dz (k))
        cv0 (k) = dm (k) * (cv_air + qv (k) * cv_vap + (qr (k) + ql (k)) * c_liq + &
             (qi (k) + qs (k) + qg (k)) * c_ice) + cw * (m1 (k) - m1 (k - 1))
        ! cvm_new + cw * m1 (k) = cvm_old + cw * m1 (k - 1)
    enddo
    ! -----------------------------------------------------------------------
    ! implicit algorithm: can't be vectorized
    ! needs an inner i - loop for vectorization
    ! -----------------------------------------------------------------------
    ! top layer: cv0 = cvn + cw * m1 (k)
    ! tz (k) = cv0 (k) * tz (k) / (cvn (k) + cw * m1 (k)) = tz (k) -- > no change
    do k = ks + 1, ke
        tz (k) = (cv0 (k) * tz (k) + m1 (k - 1) * (cw * tz (k - 1) + dgz (k))) / (cv0 (k) + cw * m1 (k - 1))
    enddo

end subroutine sedi_heat

! -----------------------------------------------------------------------
! warm rain cloud microphysics
! -----------------------------------------------------------------------

subroutine warm_rain (dt, ks, ke, dp, dz, tz, qv, ql, qr, qi, qs, qg, &
        den, denfac, ccn, c_praut, rh_rain, vtr, r1, m1_rain, w1, h_var)

    implicit none

    integer, intent (in) :: ks, ke
    real, intent (in) :: dt ! time step (s)
    real, intent (in) :: rh_rain, h_var
    real, intent (in), dimension (ks:ke) :: dp, dz, den
    real, intent (in), dimension (ks:ke) :: denfac, ccn, c_praut

    real (kind = r_grid), intent (inout), dimension (ks:ke) :: tz
    real, intent (inout), dimension (ks:ke) :: vtr, qv, ql, qr, qi, qs, qg, m1_rain, w1
    real, intent (out) :: r1
    real, parameter :: so3 = 7. / 3.
    ! fall velocity constants:
    real, parameter :: vconr = 2503.23638966667
    real, parameter :: normr = 25132741228.7183
    real, parameter :: thr = 1.e-8

    real, dimension (ks:ke) :: dl, dm
    real, dimension (ks:ke + 1) :: ze, zt
    real :: sink, dq, qc0, qc
    real :: qden
    real :: zs = 0.
    real :: dt5
    integer :: k

    logical :: no_fall

    dt5 = 0.5 * dt

    ! -----------------------------------------------------------------------
    ! terminal speed of rain
    ! -----------------------------------------------------------------------

    m1_rain (:) = 0.

    call check_column (ks, ke, qr, no_fall)

    if (no_fall) then
        vtr (:) = vf_min
        r1 = 0.
    else

        ! -----------------------------------------------------------------------
        ! fall speed of rain
        ! -----------------------------------------------------------------------

        if (const_vr) then
            vtr (:) = vr_fac ! ifs_2016: 4.0
        else
            do k = ks, ke
                qden = qr (k) * den (k)
                if (qr (k) < thr) then
                    vtr (k) = vr_min
                else
                    vtr (k) = vr_fac * vconr * sqrt (min (10., sfcrho / den (k))) * &
                        exp (0.2 * log (qden / normr))
                    vtr (k) = min (vr_max, max (vr_min, vtr (k)))
                endif
            enddo
        endif

        ze (ke + 1) = zs
        do k = ke, ks, - 1
            ze (k) = ze (k + 1) - dz (k) ! dz < 0
        enddo

        ! -----------------------------------------------------------------------
        ! evaporation and accretion of rain for the first 1 / 2 time step
        ! -----------------------------------------------------------------------

        call revap_racc (ks, ke, dt5, tz, qv, ql, qr, qi, qs, qg, den, denfac, rh_rain, h_var)

        if (do_sedi_w) then
            do k = ks, ke
                dm (k) = dp (k) * (1. + qv (k) + ql (k) + qr (k) + qi (k) + qs (k) + qg (k))
            enddo
        endif

        ! -----------------------------------------------------------------------
        ! mass flux induced by falling rain
        ! -----------------------------------------------------------------------

        if (use_ppm) then
            zt (ks) = ze (ks)
            do k = ks + 1, ke
                zt (k) = ze (k) - dt5 * (vtr (k - 1) + vtr (k))
            enddo
            zt (ke + 1) = zs - dt * vtr (ke)

            do k = ks, ke
                if (zt (k + 1) >= zt (k)) zt (k + 1) = zt (k) - dz_min
            enddo
            call lagrangian_fall_ppm (ks, ke, zs, ze, zt, dp, qr, r1, m1_rain, mono_prof)
        else
            call implicit_fall (dt, ks, ke, ze, vtr, dp, qr, r1, m1_rain)
        endif

        ! -----------------------------------------------------------------------
        ! vertical velocity transportation during sedimentation
        ! -----------------------------------------------------------------------

        if (do_sedi_w) then
            ! conservation of vertical momentum:
            w1 (ks) = w1 (ks) + m1_rain (ks) * vtr (ks) / dm (ks)
            do k = ks + 1, ke
                w1 (k) = (dm (k) * w1 (k) + m1_rain (k - 1) * (w1 (k - 1) - vtr (k - 1)) + m1_rain (k) * vtr (k)) &
                     / (dm (k) + m1_rain (k - 1))
            enddo
        endif

        ! -----------------------------------------------------------------------
        ! heat exchanges during sedimentation
        ! -----------------------------------------------------------------------

        if (do_sedi_heat) &
            call sedi_heat (ks, ke, dp, m1_rain, dz, tz, qv, ql, qr, qi, qs, qg, c_liq)

        ! -----------------------------------------------------------------------
        ! evaporation and accretion of rain for the remaing 1 / 2 time step
        ! -----------------------------------------------------------------------

        call revap_racc (ks, ke, dt5, tz, qv, ql, qr, qi, qs, qg, den, denfac, rh_rain, h_var)

    endif

    ! -----------------------------------------------------------------------
    ! auto - conversion
    ! assuming linear subgrid vertical distribution of cloud water
    ! following lin et al. 1994, mwr
    ! -----------------------------------------------------------------------

    if (irain_f /= 0) then

        ! -----------------------------------------------------------------------
        ! no subgrid varaibility
        ! -----------------------------------------------------------------------

        do k = ks, ke
            qc0 = fac_rc * ccn (k)
            if (tz (k) > t_wfr) then
                if (use_ccn) then
                    ! -----------------------------------------------------------------------
                    ! ccn is formulted as ccn = ccn_surface * (den / den_surface)
                    ! -----------------------------------------------------------------------
                    qc = qc0
                else
                    qc = qc0 / den (k)
                endif
                dq = ql (k) - qc
                if (dq > 0.) then
                    sink = min (dq, dt * c_praut (k) * den (k) * exp (so3 * log (ql (k))))
                    ql (k) = ql (k) - sink
                    qr (k) = qr (k) + sink
                endif
            endif
        enddo

    else

        ! -----------------------------------------------------------------------
        ! with subgrid varaibility
        ! -----------------------------------------------------------------------

        call linear_prof (ke - ks + 1, ql (ks), dl (ks), z_slope_liq, h_var)

        do k = ks, ke
            qc0 = fac_rc * ccn (k)
            if (tz (k) > t_wfr + dt_fr) then
                dl (k) = min (max (1.e-6, dl (k)), 0.5 * ql (k))
                ! --------------------------------------------------------------------
                ! as in klein's gfdl am2 stratiform scheme (with subgrid variations)
                ! --------------------------------------------------------------------
                if (use_ccn) then
                    ! --------------------------------------------------------------------
                    ! ccn is formulted as ccn = ccn_surface * (den / den_surface)
                    ! --------------------------------------------------------------------
                    qc = qc0
                else
                    qc = qc0 / den (k)
                endif
                dq = 0.5 * (ql (k) + dl (k) - qc)
                ! --------------------------------------------------------------------
                ! dq = dl if qc == q_minus = ql - dl
                ! dq = 0 if qc == q_plus = ql + dl
                ! --------------------------------------------------------------------
                if (dq > 0.) then ! q_plus > qc
                    ! --------------------------------------------------------------------
                    ! revised continuous form: linearly decays (with subgrid dl) to zero at qc == ql + dl
                    ! --------------------------------------------------------------------
                    sink = min (1., dq / dl (k)) * dt * c_praut (k) * den (k) * exp (so3 * log (ql (k)))
                    ql (k) = ql (k) - sink
                    qr (k) = qr (k) + sink
                endif
            endif
        enddo
    endif

end subroutine warm_rain

! -----------------------------------------------------------------------
! evaporation of rain
! -----------------------------------------------------------------------

subroutine revap_racc (ks, ke, dt, tz, qv, ql, qr, qi, qs, qg, den, denfac, rh_rain, h_var)

    implicit none

    integer, intent (in) :: ks, ke
    real, intent (in) :: dt ! time step (s)
    real, intent (in) :: rh_rain, h_var
    real, intent (in), dimension (ks:ke) :: den, denfac
    real (kind = r_grid), intent (inout), dimension (ks:ke) :: tz
    real, intent (inout), dimension (ks:ke) :: qv, qr, ql, qi, qs, qg
    ! local:
    real (kind = r_grid), dimension (ks:ke) :: cvm
    real, dimension (ks:ke) :: q_liq, q_sol, lcpk
    real :: dqv, qsat, dqsdt, evap, t2, qden, q_plus, q_minus, sink
    real :: qpz, dq, dqh, tin

    integer :: k

    do k = ks, ke

        if (tz (k) > t_wfr .and. qr (k) > qrmin) then

            ! -----------------------------------------------------------------------
            ! define heat capacity and latent heat coefficient
            ! -----------------------------------------------------------------------

            q_liq (k) = ql (k) + qr (k)
            q_sol (k) = qi (k) + qs (k) + qg (k)

            cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
            lcpk (k) = (lv00 + d1_vap * tz (k)) / cvm (k)
            tin = (tz (k) * cvm (k) - lv00 * ql (k)) / (1. + (qv (k) + ql (k)) * c1_vap + qr (k) * c1_liq + q_sol (k) * c1_ice)
            !
            qpz = qv (k) + ql (k)
            qsat = wqs2 (tin, den (k), dqsdt)
            dqh = max (ql (k), h_var * max (qpz, qcmin))
            dqh = min (dqh, 0.2 * qpz) ! new limiter
            dqv = qsat - qv (k) ! use this to prevent super - sat the gird box
            q_minus = qpz - dqh
            q_plus = qpz + dqh

            ! -----------------------------------------------------------------------
            ! qsat must be > q_minus to activate evaporation
            ! qsat must be < q_plus to activate accretion
            ! -----------------------------------------------------------------------

            ! -----------------------------------------------------------------------
            ! rain evaporation
            ! -----------------------------------------------------------------------

            if (dqv > 0. .and. qsat > q_minus) then
                if (qsat > q_plus) then
                    dq = qsat - qpz
                else
                    ! -----------------------------------------------------------------------
                    ! q_minus < qsat < q_plus
                    ! dq == dqh if qsat == q_minus
                    ! -----------------------------------------------------------------------
                    dq = 0.25 * (qsat - q_minus) ** 2 / dqh
                endif
                qden = qr (k) * den (k)
                t2 = tin * tin
                evap = crevp (1) * t2 * dq * (crevp (2) * sqrt (qden) + crevp (3) * &
                    exp (0.725 * log (qden))) / (crevp (4) * t2 + crevp (5) * qsat * den (k))
                evap = min (qr (k), dt * evap, dqv / (1. + lcpk (k) * dqsdt))
                ! -----------------------------------------------------------------------
                ! alternative minimum evap in dry environmental air
                ! sjl 20180831:
                sink = min (qr (k), dim (rh_rain * qsat, qv (k)) / (1. + lcpk (k) * dqsdt))
                evap = max (evap, sink)
                ! -----------------------------------------------------------------------
                qr (k) = qr (k) - evap
                qv (k) = qv (k) + evap
                q_liq (k) = q_liq (k) - evap
                tz (k) = (cvm (k) * tz (k) - lv00 * evap) / (one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice)
            endif

            ! -----------------------------------------------------------------------
            ! accretion: pracc
            ! -----------------------------------------------------------------------

            ! if (qr (k) > qrmin .and. ql (k) > 1.e-7 .and. qsat < q_plus) then
            if (qr (k) > 1.e-6 .and. ql (k) > 2.e-6 .and. qsat < q_minus) then
                sink = dt * denfac (k) * cracw * exp (0.95 * log (qr (k) * den (k)))
                sink = sink / (1. + sink) * ql (k)
                ql (k) = ql (k) - sink
                qr (k) = qr (k) + sink
            endif

        endif ! warm - rain
    enddo

end subroutine revap_racc

! -----------------------------------------------------------------------
! definition of vertical subgrid variability
! used for cloud ice and cloud water autoconversion
! qi -- > ql & ql -- > qr
! edges: qe == qbar + / - dm
! -----------------------------------------------------------------------

subroutine linear_prof (km, q, dm, z_var, h_var)

    implicit none

    integer, intent (in) :: km
    real, intent (in) :: q (km), h_var
    real, intent (out) :: dm (km)
    logical, intent (in) :: z_var
    real :: dq (km)
    integer :: k

    if (z_var) then
        do k = 2, km
            dq (k) = 0.5 * (q (k) - q (k - 1))
        enddo
        dm (1) = 0.

        ! -----------------------------------------------------------------------
        ! use twice the strength of the positive definiteness limiter (lin et al 1994)
        ! -----------------------------------------------------------------------

        do k = 2, km - 1
            dm (k) = 0.5 * min (abs (dq (k) + dq (k + 1)), 0.5 * q (k))
            if (dq (k) * dq (k + 1) <= 0.) then
                if (dq (k) > 0.) then ! local max
                    dm (k) = min (dm (k), dq (k), - dq (k + 1))
                else
                    dm (k) = 0.
                endif
            endif
        enddo
        dm (km) = 0.

        ! -----------------------------------------------------------------------
        ! impose a presumed background horizontal variability that is proportional to the value itself
        ! -----------------------------------------------------------------------

        do k = 1, km
            dm (k) = max (dm (k), qvmin, h_var * q (k))
        enddo
    else
        do k = 1, km
            dm (k) = max (qvmin, h_var * q (k))
        enddo
    endif

end subroutine linear_prof

! =======================================================================
! ice cloud microphysics processes
! bulk cloud micro - physics; processes splitting
! with some un - split sub - grouping
! time implicit (when possible) accretion and autoconversion
! author: shian - jiann lin, gfdl
! =======================================================================

subroutine icloud (ks, ke, tzk, p1, qvk, qlk, qrk, qik, qsk, qgk, dp1, &
        den, denfac, vts, vtg, vtr, qak, rh_adj, rh_rain, dts, h_var, last_step)

    implicit none

    logical, intent (in) :: last_step
    integer, intent (in) :: ks, ke
    real, intent (in), dimension (ks:ke) :: p1, dp1, den, denfac, vts, vtg, vtr
    real (kind = r_grid), intent (inout), dimension (ks:ke) :: tzk
    real, intent (inout), dimension (ks:ke) :: qvk, qlk, qrk, qik, qsk, qgk, qak
    real, intent (in) :: rh_adj, rh_rain, dts, h_var
    ! local:
    real, dimension (ks:ke) :: icpk, di
    real, dimension (ks:ke) :: q_liq, q_sol
    real (kind = r_grid), dimension (ks:ke) :: cvm, te8
    real (kind = r_grid) :: tz
    real :: rdts, fac_g2v, fac_v2g, fac_i2s, fac_imlt
    real :: qv, ql, qr, qi, qs, qg, melt
    real :: pracs, psacw, pgacw, psacr, pgacr, pgaci, praci, psaci
    real :: pgmlt, psmlt, pgfr, psaut
    real :: tc, dqs0, qden, qim, qsm
    real :: dt5, factor, sink, qi_crt
    real :: tmp, qsw, qsi, dqsdt, dq
    real :: dtmp, qc, q_plus, q_minus
    integer :: k

    dt5 = 0.5 * dts
    rdts = 1. / dts

    ! -----------------------------------------------------------------------
    ! define conversion scalar / factor
    ! -----------------------------------------------------------------------

    fac_i2s = 1. - exp (- dts / tau_i2s)
    fac_g2v = 1. - exp (- dts / tau_g2v)
    fac_v2g = 1. - exp (- dts / tau_v2g)
    fac_imlt = 1. - exp (- dt5 / tau_imlt)

    ! -----------------------------------------------------------------------
    ! define heat capacity and latend heat coefficient
    ! -----------------------------------------------------------------------

    do k = ks, ke
        q_liq (k) = qlk (k) + qrk (k)
        q_sol (k) = qik (k) + qsk (k) + qgk (k)
        cvm (k) = one_r8 + qvk (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
        te8 (k) = cvm (k) * tzk (k) + lv00 * qvk (k) - li00 * q_sol (k)
        icpk (k) = (li00 + d1_ice * tzk (k)) / cvm (k)
    enddo

    ! -----------------------------------------------------------------------
    ! sources of cloud ice: pihom, cold rain, and the sat_adj
    ! (initiation plus deposition)
    ! sources of snow: cold rain, auto conversion + accretion (from cloud ice)
    ! sat_adj (deposition; requires pre - existing snow) ; initial snow comes from auto conversion
    ! -----------------------------------------------------------------------

    do k = ks, ke
        if (tzk (k) > tice .and. qik (k) > qcmin) then

            ! -----------------------------------------------------------------------
            ! pimlt: instant melting of cloud ice
            ! -----------------------------------------------------------------------

            melt = min (qik (k), fac_imlt * (tzk (k) - tice) / icpk (k))
            tmp = min (melt, dim (ql_mlt, qlk (k))) ! max ql amount
            qlk (k) = qlk (k) + tmp
            qrk (k) = qrk (k) + melt - tmp
            qik (k) = qik (k) - melt
            q_liq (k) = q_liq (k) + melt
            q_sol (k) = q_sol (k) - melt
        elseif (tzk (k) < t_wfr .and. qlk (k) > qcmin) then

            ! -----------------------------------------------------------------------
            ! pihom: homogeneous freezing of cloud water into cloud ice
            ! -----------------------------------------------------------------------

            dtmp = t_wfr - tzk (k)
            factor = min (1., dtmp / dt_fr)
            sink = min (qlk (k) * factor, dtmp / icpk (k))
            qi_crt = qi_gen * min (qi_lim, 0.1 * (tice - tzk (k))) / den (k)
            tmp = min (sink, dim (qi_crt, qik (k)))
            qlk (k) = qlk (k) - sink
            qsk (k) = qsk (k) + sink - tmp
            qik (k) = qik (k) + tmp
            q_liq (k) = q_liq (k) - sink
            q_sol (k) = q_sol (k) + sink
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! vertical subgrid variability
    ! -----------------------------------------------------------------------

    call linear_prof (ke - ks + 1, qik (ks), di (ks), z_slope_ice, h_var)

    ! -----------------------------------------------------------------------
    ! update capacity heat and latend heat coefficient
    ! -----------------------------------------------------------------------

    do k = ks, ke
        cvm (k) = one_r8 + qvk (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
        tzk (k) = (te8 (k) - lv00 * qvk (k) + li00 * q_sol (k)) / cvm (k)
        icpk (k) = (li00 + d1_ice * tzk (k)) / cvm (k)
    enddo

    do k = ks, ke

        ! -----------------------------------------------------------------------
        ! do nothing above p_min
        ! -----------------------------------------------------------------------

        if (p1 (k) < p_min) cycle

        tz = tzk (k)
        qv = qvk (k)
        ql = qlk (k)
        qi = qik (k)
        qr = qrk (k)
        qs = qsk (k)
        qg = qgk (k)

        pgacr = 0.
        pgacw = 0.
        tc = tz - tice

        if (tc .ge. 0.) then

            ! -----------------------------------------------------------------------
            ! melting of snow
            ! -----------------------------------------------------------------------

            dqs0 = ces0 / p1 (k) - qv ! not sure if this is correct; check again

            if (qs > qcmin) then

                ! -----------------------------------------------------------------------
                ! psacw: accretion of cloud water by snow
                ! only rate is used (for snow melt) since tc > 0.
                ! -----------------------------------------------------------------------

                if (ql > qrmin) then
                    factor = denfac (k) * csacw * exp (0.8125 * log (qs * den (k)))
                    psacw = factor / (1. + dts * factor) * ql ! rate
                else
                    psacw = 0.
                endif

                ! -----------------------------------------------------------------------
                ! psacr: accretion of rain by melted snow
                ! pracs: accretion of snow by rain
                ! -----------------------------------------------------------------------

                if (qr > qrmin) then
                    psacr = min (acr3d (vts (k), vtr (k), qr, qs, csacr, acco (1, 2), &
                        den (k)), qr * rdts)
                    pracs = acr3d (vtr (k), vts (k), qs, qr, cracs, acco (1, 1), den (k))
                else
                    psacr = 0.
                    pracs = 0.
                endif

                ! -----------------------------------------------------------------------
                ! total snow sink:
                ! psmlt: snow melt (due to rain accretion)
                ! -----------------------------------------------------------------------

                psmlt = max (0., smlt (tc, dqs0, qs * den (k), psacw, psacr, csmlt, &
                    den (k), denfac (k)))
                sink = min (qs, dts * (psmlt + pracs), tc / icpk (k))
                qs = qs - sink

                ! melt all snow if t > 12 c
                if (qs > qcmin .and. tz > tice + 12.) then
                    sink = sink + qs
                    qs = 0.
                endif

                tmp = min (sink, dim (qs_mlt, ql)) ! max ql due to snow melt
                ql = ql + tmp
                qr = qr + sink - tmp
                q_liq (k) = q_liq (k) + sink
                q_sol (k) = q_sol (k) - sink

                cvm (k) = one_r8 + qv * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
                tz = (te8 (k) - lv00 * qv + li00 * q_sol (k)) / cvm (k)
                tc = tz - tice
                icpk (k) = (li00 + d1_ice * tz) / cvm (k)

            endif

            ! -----------------------------------------------------------------------
            ! update capacity heat and latend heat coefficient
            ! -----------------------------------------------------------------------


            ! -----------------------------------------------------------------------
            ! melting of graupel
            ! -----------------------------------------------------------------------

            if (qg > qcmin .and. tc > 0.) then

                ! -----------------------------------------------------------------------
                ! pgacr: accretion of rain by graupel
                ! -----------------------------------------------------------------------

                if (qr > qrmin) &
                    pgacr = min (acr3d (vtg (k), vtr (k), qr, qg, cgacr, acco (1, 3), &
                    den (k)), rdts * qr)

                ! -----------------------------------------------------------------------
                ! pgacw: accretion of cloud water by graupel
                ! -----------------------------------------------------------------------

                qden = qg * den (k)
                if (ql > qrmin) then
                    factor = cgacw * qden / sqrt (den (k) * sqrt (sqrt (qden)))
                    pgacw = factor / (1. + dts * factor) * ql ! rate
                endif

                ! -----------------------------------------------------------------------
                ! pgmlt: graupel melt
                ! -----------------------------------------------------------------------

                pgmlt = dts * gmlt (tc, dqs0, qden, pgacw, pgacr, cgmlt, den (k))
                pgmlt = min (max (0., pgmlt), qg, tc / icpk (k))
                qg = qg - pgmlt
                qr = qr + pgmlt
                q_liq (k) = q_liq (k) + pgmlt
                q_sol (k) = q_sol (k) - pgmlt
                cvm (k) = one_r8 + qv * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
                tz = (te8 (k) - lv00 * qv + li00 * q_sol (k)) / cvm (k)
            endif

        else

            ! -----------------------------------------------------------------------
            ! cloud ice proc:
            ! -----------------------------------------------------------------------

            ! -----------------------------------------------------------------------
            ! psaci: accretion of cloud ice by snow
            ! -----------------------------------------------------------------------

            if (qi > 1.e-6) then ! cloud ice sink terms
                if (qs > 1.e-6) then
                    ! -----------------------------------------------------------------------
                    ! sjl added (following lin eq. 23) the temperature dependency
                    ! to reduce accretion, use esi = exp (0.05 * tc) as in hong et al 2004
                    ! -----------------------------------------------------------------------
                    factor = dts * denfac (k) * csaci * exp (0.05 * tc + 0.8125 * log (qs * den (k)))
                    psaci = factor / (1. + factor) * qi
                else
                    psaci = 0.
                endif

                ! -----------------------------------------------------------------------
                ! pasut: autoconversion: cloud ice -- > snow
                ! -----------------------------------------------------------------------

                ! -----------------------------------------------------------------------
                ! similar to lfo 1983: eq. 21 solved implicitly
                ! threshold from wsm6 scheme, hong et al 2004, eq (13) : qi0_crt ~0.8e-4
                ! -----------------------------------------------------------------------

                if (qi0_crt < 0.) then
                    qim = - qi0_crt
                else
                    qim = qi0_crt / den (k)
                endif

                ! -----------------------------------------------------------------------
                ! assuming linear subgrid vertical distribution of cloud ice
                ! the mismatch computation following lin et al. 1994, mwr
                ! -----------------------------------------------------------------------

                if (const_vi) then
                    tmp = fac_i2s
                else
                    tmp = fac_i2s * exp (0.025 * tc)
                endif

                di (k) = max (di (k), qrmin)
                q_plus = qi + di (k)
                if (q_plus > (qim + qrmin)) then
                    if (qim > (qi - di (k))) then
                        dq = (0.25 * (q_plus - qim) ** 2) / di (k)
                    else
                        dq = qi - qim
                    endif
                    psaut = tmp * dq
                else
                    psaut = 0.
                endif
                ! -----------------------------------------------------------------------
                ! sink is no greater than 75% of qi
                ! -----------------------------------------------------------------------
                sink = min (0.75 * qi, psaci + psaut)
                qi = qi - sink
                qs = qs + sink

                ! -----------------------------------------------------------------------
                ! pgaci: accretion of cloud ice by graupel
                ! -----------------------------------------------------------------------

                if (qg > 3.e-6) then
                    ! -----------------------------------------------------------------------
                    ! factor = dts * cgaci / sqrt (den (k)) * exp (0.05 * tc + 0.875 * log (qg * den (k)))
                    ! simplified form: remove temp dependency & set the exponent "0.875" -- > 1
                    ! -----------------------------------------------------------------------
                    factor = dts * cgaci * sqrt (den (k)) * qg
                    pgaci = factor / (1. + factor) * qi
                    qi = qi - pgaci
                    qg = qg + pgaci
                endif

            endif

            ! -----------------------------------------------------------------------
            ! cold - rain proc:
            ! -----------------------------------------------------------------------

            ! -----------------------------------------------------------------------
            ! rain to ice, snow, graupel processes:
            ! -----------------------------------------------------------------------

            tc = tz - tice

            if (qr > 1.e-6 .and. tc < 0.) then

                ! -----------------------------------------------------------------------
                ! * sink * terms to qr: psacr + pgfr
                ! source terms to qs: psacr
                ! source terms to qg: pgfr
                ! -----------------------------------------------------------------------

                ! -----------------------------------------------------------------------
                ! psacr accretion of rain by snow
                ! -----------------------------------------------------------------------

                if (qs > 1.e-6) then ! if snow exists
                    psacr = dts * acr3d (vts (k), vtr (k), qr, qs, csacr, acco (1, 2), den (k))
                else
                    psacr = 0.
                endif

                ! -----------------------------------------------------------------------
                ! pgfr: rain freezing -- > graupel
                ! -----------------------------------------------------------------------

                pgfr = dts * cgfr (1) / den (k) * (exp (- cgfr (2) * tc) - 1.) * &
                    exp (1.75 * log (qr * den (k)))

                ! -----------------------------------------------------------------------
                ! total sink to qr
                ! -----------------------------------------------------------------------

                sink = psacr + pgfr
                factor = min (sink, qr, - tc / icpk (k)) / max (sink, qrmin)

                psacr = factor * psacr
                pgfr = factor * pgfr

                sink = psacr + pgfr
                qr = qr - sink
                qs = qs + psacr
                qg = qg + pgfr
                q_liq (k) = q_liq (k) - sink
                q_sol (k) = q_sol (k) + sink

                cvm (k) = one_r8 + qv * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
                tz = (te8 (k) - lv00 * qv + li00 * q_sol (k)) / cvm (k)
                icpk (k) = (li00 + d1_ice * tz) / cvm (k)
            endif

            ! -----------------------------------------------------------------------
            ! graupel production terms:
            ! -----------------------------------------------------------------------

            if (qs > 3.e-6) then

                ! -----------------------------------------------------------------------
                ! accretion: snow -- > graupel
                ! -----------------------------------------------------------------------

                if (qg > qrmin) then
                    sink = dts * acr3d (vtg (k), vts (k), qs, qg, cgacs, acco (1, 4), den (k))
                else
                    sink = 0.
                endif

                ! -----------------------------------------------------------------------
                ! autoconversion snow -- > graupel
                ! -----------------------------------------------------------------------

                qsm = qs0_crt / den (k)
                if (qs > qsm) then
                    factor = dts * 1.e-3 * exp (0.09 * (tz - tice))
                    sink = sink + factor / (1. + factor) * (qs - qsm)
                endif
                sink = min (qs, sink)
                qs = qs - sink
                qg = qg + sink

            endif ! snow existed

            if (qg > 1.e-6 .and. tz < tice0) then

                ! -----------------------------------------------------------------------
                ! pgacw: accretion of cloud water by graupel
                ! -----------------------------------------------------------------------

                if (ql > 1.e-6) then
                    qden = qg * den (k)
                    factor = dts * cgacw * qden / sqrt (den (k) * sqrt (sqrt (qden)))
                    pgacw = factor / (1. + factor) * ql
                else
                    pgacw = 0.
                endif

                ! -----------------------------------------------------------------------
                ! pgacr: accretion of rain by graupel
                ! -----------------------------------------------------------------------

                if (qr > 1.e-6) then
                    pgacr = min (dts * acr3d (vtg (k), vtr (k), qr, qg, cgacr, acco (1, 3), &
                        den (k)), qr)
                else
                    pgacr = 0.
                endif

                sink = pgacr + pgacw
                factor = min (sink, dim (tice, tz) / icpk (k)) / max (sink, qrmin)
                pgacr = factor * pgacr
                pgacw = factor * pgacw

                sink = pgacr + pgacw
                qg = qg + sink
                qr = qr - pgacr
                ql = ql - pgacw

                q_liq (k) = q_liq (k) - sink
                q_sol (k) = q_sol (k) + sink
                tz = (te8 (k) - lv00 * qv + li00 * q_sol (k)) / (one_r8 + qv * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice)
            endif

        endif

        tzk (k) = tz
        qvk (k) = qv
        qlk (k) = ql
        qik (k) = qi
        qrk (k) = qr
        qsk (k) = qs
        qgk (k) = qg

    enddo

    call subgrid_z_proc (ks, ke, p1, den, denfac, dts, rh_adj, tzk, qvk, &
        qlk, qrk, qik, qsk, qgk, qak, h_var, rh_rain, te8, last_step)

end subroutine icloud

! =======================================================================
! temperature sentive high vertical resolution processes
! =======================================================================

subroutine subgrid_z_proc (ks, ke, p1, den, denfac, dts, rh_adj, tz, qv, &
        ql, qr, qi, qs, qg, qa, h_var, rh_rain, te8, last_step)

    implicit none

    integer, intent (in) :: ks, ke
    real, intent (in) :: dts, rh_adj, h_var, rh_rain
    real, intent (in), dimension (ks:ke) :: p1, den, denfac
    real (kind = r_grid), intent (in), dimension (ks:ke) :: te8
    real (kind = r_grid), intent (inout), dimension (ks:ke) :: tz
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, qa
    logical, intent (in) :: last_step
    ! local:
    real, dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    real, dimension (ks:ke) :: q_liq, q_sol, q_cond
    real (kind = r_grid), dimension (ks:ke) :: cvm
    real :: pidep, qi_crt
    ! -----------------------------------------------------------------------
    ! qstar over water may be accurate only down to - 80 deg c with ~10% uncertainty
    ! must not be too large to allow psc
    ! -----------------------------------------------------------------------
    real :: rh, rqi, tin, qsw, qsi, qpz, qstar
    real :: dqsdt, dwsdt, dq, dq0, factor, tmp, liq, ice
    real :: q_plus, q_minus
    real :: evap, sink, tc, dtmp
    real :: pssub, pgsub, tsq, qden
    real :: fac_l2v, fac_v2l, fac_g2v, fac_v2g
    integer :: k

    ! -----------------------------------------------------------------------
    ! define conversion scalar / factor
    ! -----------------------------------------------------------------------

    fac_l2v = 1. - exp (- dts / tau_l2v)
    fac_v2l = 1. - exp (- dts / tau_v2l)
    fac_g2v = 1. - exp (- dts / tau_g2v)
    fac_v2g = 1. - exp (- dts / tau_v2g)

    ! -----------------------------------------------------------------------
    ! define heat capacity and latend heat coefficient
    ! -----------------------------------------------------------------------

    do k = ks, ke
        q_liq (k) = ql (k) + qr (k)
        q_sol (k) = qi (k) + qs (k) + qg (k)
        cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
        lcpk (k) = (lv00 + d1_vap * tz (k)) / cvm (k)
        icpk (k) = (li00 + d1_ice * tz (k)) / cvm (k)
        tcp3 (k) = lcpk (k) + icpk (k) * min (1., dim (tice, tz (k)) / (tice - t_wfr))
    enddo

    do k = ks, ke

        if (p1 (k) < p_min) cycle

        ! -----------------------------------------------------------------------
        ! instant deposit all water vapor to cloud ice when temperature is super low
        ! -----------------------------------------------------------------------

        if (tz (k) < t_min) then
            sink = dim (qv (k), 1.e-7)
            qv (k) = qv (k) - sink
            qi (k) = qi (k) + sink
            q_sol (k) = q_sol (k) + sink
            tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / (one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice)
            if (do_qa) qa (k) = 1. ! air fully saturated; 100 % cloud cover
            cycle
        endif

        ! -----------------------------------------------------------------------
        ! instant evaporation / sublimation of all clouds if rh < rh_adj -- > cloud free
        ! -----------------------------------------------------------------------
        ! rain water is handled in warm - rain process.
        qpz = qv (k) + ql (k) + qi (k) + qs (k)
        tin = (te8 (k) - lv00 * qpz + li00 * qg (k)) / (one_r8 + qpz * c1_vap + qr (k) * c1_liq + qg (k) * c1_ice)
        if (tin > t_sub + 6.) then
            rh = qpz / iqs1 (tin, den (k))
            if (rh < rh_adj) then ! qpz / rh_adj < qs
                tz (k) = tin
                qv (k) = qpz
                ql (k) = 0.
                qi (k) = 0.
                qs (k) = 0.
                cycle ! cloud free
            endif
        endif

        ! -----------------------------------------------------------------------
        ! cloud water < -- > vapor adjustment:
        ! -----------------------------------------------------------------------

        tin = tz (k)
        qsw = wqs2 (tin, den (k), dwsdt)
        dq0 = qsw - qv (k)
        if (dq0 > 0.) then ! evaporation
            factor = min (1., fac_l2v * (10. * dq0 / qsw)) ! the rh dependent factor = 1 at 90%
            evap = min (ql (k), factor * dq0 / (1. + tcp3 (k) * dwsdt))
        elseif (do_cond_timescale) then
           factor = min ( 1., fac_v2l * ( 10. * (-dq0) / qsw ))
           evap = - min ( qv (k), factor * -dq0 / (1. + tcp3 (k) * dwsdt))
        else ! condensate all excess vapor into cloud water
            evap = dq0 / (1. + tcp3 (k) * dwsdt)
        endif
        ! sjl on jan 23 2018: reversible evap / condensation:
        qv (k) = qv (k) + evap
        ql (k) = ql (k) - evap
        q_liq (k) = q_liq (k) - evap

        cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
        tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / cvm (k)

        ! -----------------------------------------------------------------------
        ! update heat capacity and latend heat coefficient
        ! -----------------------------------------------------------------------

        icpk (k) = (li00 + d1_ice * tz (k)) / cvm (k)

        ! -----------------------------------------------------------------------
        ! enforce complete freezing below - 48 c
        ! -----------------------------------------------------------------------

        dtmp = t_wfr - tz (k) ! [ - 40, - 48]
        if (dtmp > 0. .and. ql (k) > qcmin) then
            sink = min (ql (k), ql (k) * dtmp * 0.125, dtmp / icpk (k))
            ql (k) = ql (k) - sink
            qi (k) = qi (k) + sink
            q_liq (k) = q_liq (k) - sink
            q_sol (k) = q_sol (k) + sink
            cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
            tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / cvm (k)
            icpk (k) = (li00 + d1_ice * tz (k)) / cvm (k)
        endif

        ! -----------------------------------------------------------------------
        ! bigg mechanism
        ! -----------------------------------------------------------------------

        tc = tice - tz (k)
        if (ql (k) > qrmin .and. tc > 0.1) then
            sink = 3.3333e-10 * dts * (exp (0.66 * tc) - 1.) * den (k) * ql (k) * ql (k)
            sink = min (ql (k), tc / icpk (k), sink)
            ql (k) = ql (k) - sink
            qi (k) = qi (k) + sink
            q_liq (k) = q_liq (k) - sink
            q_sol (k) = q_sol (k) + sink
            cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
            tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / cvm (k)
        endif ! significant ql existed

        ! -----------------------------------------------------------------------
        ! update capacity heat and latend heat coefficient
        ! -----------------------------------------------------------------------

        tcpk (k) = (li20 + (d1_vap + d1_ice) * tz (k)) / cvm (k)

        ! -----------------------------------------------------------------------
        ! sublimation / deposition of ice
        ! -----------------------------------------------------------------------

        if (tz (k) < tice) then
            qsi = iqs2 (tz (k), den (k), dqsdt)
            dq = qv (k) - qsi
            sink = dq / (1. + tcpk (k) * dqsdt)
            if (qi (k) > qrmin) then
                ! eq 9, hong et al. 2004, mwr
                ! for a and b, see dudhia 1989: page 3103 eq (b7) and (b8)
                pidep = dts * dq * 349138.78 * exp (0.875 * log (qi (k) * den (k))) &
                     / (qsi * den (k) * lat2 / (0.0243 * rvgas * tz (k) ** 2) + 4.42478e4)
            else
                pidep = 0.
            endif
            if (dq > 0.) then ! vapor - > ice
                tmp = tice - tz (k)
                ! 20160912: the following should produce more ice at higher altitude
                ! qi_crt = 4.92e-11 * exp (1.33 * log (1.e3 * exp (0.1 * tmp))) / den (k)
                qi_crt = qi_gen * min (qi_lim, 0.1 * tmp) / den (k)
                sink = min (sink, max (qi_crt - qi (k), pidep), tmp / tcpk (k))
            else ! ice -- > vapor
                pidep = pidep * min (1., dim (tz (k), t_sub) * 0.2)
                sink = max (pidep, sink, - qi (k))
            endif
            qv (k) = qv (k) - sink
            qi (k) = qi (k) + sink
            q_sol (k) = q_sol (k) + sink
            cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
            tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / cvm (k)
        endif

        ! -----------------------------------------------------------------------
        ! update capacity heat and latend heat coefficient
        ! -----------------------------------------------------------------------

        tcpk (k) = (li20 + (d1_vap + d1_ice) * tz (k)) / cvm (k)

        ! -----------------------------------------------------------------------
        ! sublimation / deposition of snow
        ! this process happens for all temp rage
        ! -----------------------------------------------------------------------

        if (qs (k) > qrmin) then
            qsi = iqs2 (tz (k), den (k), dqsdt)
            qden = qs (k) * den (k)
            tmp = exp (0.65625 * log (qden))
            tsq = tz (k) * tz (k)
            dq = (qsi - qv (k)) / (1. + tcpk (k) * dqsdt)
            pssub = cssub (1) * tsq * (cssub (2) * sqrt (qden) + cssub (3) * tmp * &
                sqrt (denfac (k))) / (cssub (4) * tsq + cssub (5) * qsi * den (k))
            pssub = (qsi - qv (k)) * dts * pssub
            if (pssub > 0.) then ! qs -- > qv, sublimation
                pssub = min (pssub * min (1., dim (tz (k), t_sub) * 0.2), qs (k))
            else
                if (tz (k) > tice) then
                    pssub = 0. ! no deposition
                else
                    pssub = max (pssub, dq, (tz (k) - tice) / tcpk (k))
                endif
            endif

            ! *******************************
            ! evap all snow if tz (k) > 12. c
            !s ******************************
            if (tz (k) > tice + 12.) then
                tmp = qs (k) - pssub
                if (tmp > 0.) pssub = pssub + tmp
            endif

            qs (k) = qs (k) - pssub
            qv (k) = qv (k) + pssub
            q_sol (k) = q_sol (k) - pssub
            cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
            tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / cvm (k)
            tcpk (k) = (li20 + (d1_vap + d1_ice) * tz (k)) / cvm (k)
        endif

        ! -----------------------------------------------------------------------
        ! simplified 2 - way grapuel sublimation - deposition mechanism
        ! -----------------------------------------------------------------------
        if (qg (k) > qrmin) then
            qsi = iqs2 (tz (k), den (k), dqsdt)
            dq = (qv (k) - qsi) / (1. + tcpk (k) * dqsdt)
            pgsub = (qv (k) / qsi - 1.) * qg (k)
            if (pgsub > 0.) then ! deposition
                if (tz (k) > tice .or. qg (k) < 1.e-6) then
                    pgsub = 0. ! no deposition
                else
                    pgsub = min (fac_v2g * pgsub, 0.2 * dq, ql (k) + qr (k), &
                         (tice - tz (k)) / tcpk (k))
                endif
            else ! submilation
                pgsub = max (fac_g2v * pgsub, dq) * min (1., dim (tz (k), t_sub) * 0.1)
            endif
            qg (k) = qg (k) + pgsub
            qv (k) = qv (k) - pgsub
            q_sol (k) = q_sol (k) + pgsub
            cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
            tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / cvm (k)
        endif

        ! -----------------------------------------------------------------------
        ! update capacity heat and latend heat coefficient
        ! -----------------------------------------------------------------------

        ! lcpk (k) = (lv00 + d1_vap * tz (k)) / cvm (k)
        ! icpk (k) = (li00 + d1_ice * tz (k)) / cvm (k)

        ! -----------------------------------------------------------------------
        ! compute cloud fraction
        ! -----------------------------------------------------------------------

        ! -----------------------------------------------------------------------
        ! combine water species
        ! -----------------------------------------------------------------------

        if (.not. (do_qa .and. last_step)) cycle

        ice = q_sol (k)
        if (rad_snow) then
            if (rad_graupel) then
                q_sol (k) = qi (k) + qs (k) + qg (k)
            else
                q_sol (k) = qi (k) + qs (k)
            endif
        else
            q_sol (k) = qi (k)
        endif
        liq = q_liq (k)
        if (rad_rain) then
            q_liq (k) = ql (k) + qr (k)
        else
            q_liq (k) = ql (k)
        endif

        q_cond (k) = q_liq (k) + q_sol (k)
        qpz = qv (k) + q_cond (k)

        ! -----------------------------------------------------------------------
        ! use the "liquid - frozen water temperature" (tin) to compute saturated specific humidity
        ! -----------------------------------------------------------------------
        ! tin = tz (k) - (lcpk (k) * q_cond (k) + icpk (k) * q_sol (k)) ! minimum temperature
        !! tin = (tz (k) * cvm (i) + li00 * q_sol (k) - lv00 * q_cond (k)) / &
        !! (one_r8 + (qv (k) + q_cond (k)) * c1_vap)
        ice = ice - q_sol (k)
        liq = liq - q_liq (k)
        tin = (te8 (k) - lv00 * qpz + li00 * ice) / (one_r8 + qpz * c1_vap + liq * c1_liq + ice * c1_ice)
        ! -----------------------------------------------------------------------
        ! determine saturated specific humidity
        ! -----------------------------------------------------------------------

        if (tin <= t_wfr) then
            ! ice phase:
            qstar = iqs1 (tin, den (k))
        elseif (tin >= tice) then
            ! liquid phase:
            qstar = wqs1 (tin, den (k))
        else
            ! mixed phase:
            qsi = iqs1 (tin, den (k))
            qsw = wqs1 (tin, den (k))
            if (q_cond (k) > 3.e-6) then
                rqi = q_sol (k) / q_cond (k)
            else
                ! -----------------------------------------------------------------------
                ! mostly liquid water q_cond (k) at initial cloud development stage
                ! -----------------------------------------------------------------------
                rqi = (tice - tin) / (tice - t_wfr)
            endif
            qstar = rqi * qsi + (1. - rqi) * qsw
        endif

        ! -----------------------------------------------------------------------
        ! assuming subgrid linear distribution in horizontal; this is effectively a smoother for the
        ! binary cloud scheme
        ! -----------------------------------------------------------------------


        ! -----------------------------------------------------------------------
        ! partial cloudiness by pdf:
        ! assuming subgrid linear distribution in horizontal; this is effectively a smoother for the
        ! binary cloud scheme; qa = 0.5 if qstar == qpz
        ! -----------------------------------------------------------------------

        qpz = cld_fac * qpz
        rh = qpz / qstar

        ! -----------------------------------------------------------------------
        ! icloud_f = 0: bug - fixed
        ! icloud_f = 1: old fvgfs gfdl) mp implementation
        ! icloud_f = 2: binary cloud scheme (0 / 1)
        ! -----------------------------------------------------------------------

        if (rh > 0.80 .and. qpz > 1.e-6) then

            dq = h_var * qpz
            q_plus = qpz + dq
            q_minus = qpz - dq

            if (icloud_f == 2) then
                if (qstar < qpz) then
                    qa (k) = 1.
                else
                    qa (k) = 0.
                endif
            else
                if (qstar < q_minus) then
                    qa (k) = 1.
                else
                    if (qstar < q_plus) then
                        if (icloud_f == 0) then
                            qa (k) = (q_plus - qstar) / (dq + dq)
                        else
                            qa (k) = (q_plus - qstar) / (2. * dq * (1. - q_cond (k)))
                        endif
                    else
                        qa (k) = 0.
                    endif
                    ! impose minimum cloudiness if substantial q_cond (k) exist
                    if (q_cond (k) > 1.e-6) then
                        qa (k) = max (cld_min, qa (k))
                    endif
                    qa (k) = min (1., qa (k))
                endif
            endif
        else
            qa (k) = 0.
        endif

    enddo

end subroutine subgrid_z_proc

! =======================================================================
! compute terminal fall speed
! consider cloud ice, snow, and graupel's melting during fall
! =======================================================================

subroutine terminal_fall (dtm, ks, ke, tz, qv, ql, qr, qg, qs, qi, dz, dp, &
        den, vtg, vts, vti, r1, g1, s1, i1, m1_sol, w1)

    implicit none

    integer, intent (in) :: ks, ke
    real, intent (in) :: dtm ! time step (s)
    real, intent (in), dimension (ks:ke) :: vtg, vts, vti, den, dp, dz
    real (kind = r_grid), intent (inout), dimension (ks:ke) :: tz
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qg, qs, qi, m1_sol, w1
    real, intent (out) :: r1, g1, s1, i1
    ! local:
    real, dimension (ks:ke + 1) :: ze, zt
    real :: qsat, dqsdt, dt5, evap, dtime
    real :: factor, frac
    real :: tmp, precip, tc, sink
    real, dimension (ks:ke) :: lcpk, icpk, cvm, q_liq, q_sol
    real, dimension (ks:ke) :: m1, dm
    real :: zs = 0.
    real :: fac_imlt

    integer :: k, k0, m
    logical :: no_fall

    dt5 = 0.5 * dtm
    fac_imlt = 1. - exp (- dt5 / tau_imlt)

    ! -----------------------------------------------------------------------
    ! define heat capacity and latend heat coefficient
    ! -----------------------------------------------------------------------

    do k = ks, ke
        m1_sol (k) = 0.
        q_liq (k) = ql (k) + qr (k)
        q_sol (k) = qi (k) + qs (k) + qg (k)
        cvm (k) = 1. + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
        lcpk (k) = (lv00 + d1_vap * tz (k)) / cvm (k)
        icpk (k) = (li00 + d1_ice * tz (k)) / cvm (k)
    enddo

    ! -----------------------------------------------------------------------
    ! find significant melting level
    ! -----------------------------------------------------------------------

    k0 = ke
    do k = ks, ke - 1
        if (tz (k) > tice) then
            k0 = k
            exit
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! melting of cloud_ice (before fall) :
    ! -----------------------------------------------------------------------

    do k = k0, ke
        tc = tz (k) - tice
        if (qi (k) > qcmin .and. tc > 0.) then
            sink = min (qi (k), fac_imlt * tc / icpk (k))
            tmp = min (sink, dim (ql_mlt, ql (k)))
            ql (k) = ql (k) + tmp
            qr (k) = qr (k) + sink - tmp
            qi (k) = qi (k) - sink
            q_liq (k) = q_liq (k) + sink
            q_sol (k) = q_sol (k) - sink
            tz (k) = tz (k) * cvm (k) - li00 * sink
            cvm (k) = 1. + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
            tz (k) = tz (k) / cvm (k)
            tc = tz (k) - tice
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! turn off melting when cloud microphysics time step is small
    ! -----------------------------------------------------------------------

    ! sjl, turn off melting of falling cloud ice, snow and graupel
    ! if (dtm < 60.) k0 = ke
    k0 = ke
    ! sjl, turn off melting of falling cloud ice, snow and graupel

    ze (ke + 1) = zs
    do k = ke, ks, - 1
        ze (k) = ze (k + 1) - dz (k) ! dz < 0
    enddo

    zt (ks) = ze (ks)

    ! -----------------------------------------------------------------------
    ! update capacity heat and latend heat coefficient
    ! -----------------------------------------------------------------------

    do k = k0, ke
        icpk (k) = (li00 + d1_ice * tz (k)) / cvm (k)
    enddo

    ! -----------------------------------------------------------------------
    ! melting of falling cloud ice into rain
    ! -----------------------------------------------------------------------

    call check_column (ks, ke, qi, no_fall)

    if (vi_fac < 1.e-5 .or. no_fall) then
        i1 = 0.
    else

        do k = ks + 1, ke
            zt (k) = ze (k) - dt5 * (vti (k - 1) + vti (k))
        enddo
        zt (ke + 1) = zs - dtm * vti (ke)

        do k = ks, ke
            if (zt (k + 1) >= zt (k)) zt (k + 1) = zt (k) - dz_min
        enddo

        if (k0 < ke) then
            do k = ke - 1, k0, - 1
                if (qi (k) > qrmin) then
                    do m = k + 1, ke
                        if (zt (k + 1) >= ze (m)) exit
                        if (zt (k) < ze (m + 1) .and. tz (m) > tice) then
                            dtime = min (1.0, (ze (m) - ze (m + 1)) / (max (vr_min, vti (k)) * tau_imlt))
                            sink = min (qi (k) * dp (k) / dp (m), dtime * (tz (m) - tice) / icpk (m))
                            tmp = min (sink, dim (ql_mlt, ql (m)))
                            ql (m) = ql (m) + tmp
                            qr (m) = qr (m) - tmp + sink
                            qi (k) = qi (k) - sink * dp (m) / dp (k)
                            tz (m) = (tz (m) * cvm (m) - li00 * sink) / &
                                 (1. + qv (m) * c1_vap + (ql (m) + qr (m)) * c1_liq + (qi (m) + qs (m) + qg (m)) * c1_ice)
                        endif
                    enddo
                endif
            enddo
        endif

        if (do_sedi_w) then
            do k = ks, ke
                dm (k) = dp (k) * (1. + qv (k) + ql (k) + qr (k) + qi (k) + qs (k) + qg (k))
            enddo
        endif

        if (use_ppm_ice) then
            call lagrangian_fall_ppm (ks, ke, zs, ze, zt, dp, qi, i1, m1_sol, mono_prof)
        else
            call implicit_fall (dtm, ks, ke, ze, vti, dp, qi, i1, m1_sol)
        endif

        if (do_sedi_w) then
            w1 (ks) = w1 (ks) + m1_sol (ks) * vti (ks) / dm (ks)
            do k = ks + 1, ke
                w1 (k) = (dm (k) * w1 (k) + m1_sol (k - 1) * (w1 (k - 1) - vti (k - 1)) + m1_sol (k) * vti (k)) &
                     / (dm (k) + m1_sol (k - 1))
            enddo
        endif

    endif

    ! -----------------------------------------------------------------------
    ! melting of falling snow into rain
    ! -----------------------------------------------------------------------

    r1 = 0.

    call check_column (ks, ke, qs, no_fall)

    if (no_fall) then
        s1 = 0.
    else

        do k = ks + 1, ke
            zt (k) = ze (k) - dt5 * (vts (k - 1) + vts (k))
        enddo
        zt (ke + 1) = zs - dtm * vts (ke)

        do k = ks, ke
            if (zt (k + 1) >= zt (k)) zt (k + 1) = zt (k) - dz_min
        enddo

        if (k0 < ke) then
            do k = ke - 1, k0, - 1
                if (qs (k) > qrmin) then
                    do m = k + 1, ke
                        if (zt (k + 1) >= ze (m)) exit
                        dtime = min (dtm, (ze (m) - ze (m + 1)) / (vr_min + vts (k)))
                        if (zt (k) < ze (m + 1) .and. tz (m) > tice) then
                            dtime = min (1.0, dtime / tau_smlt)
                            sink = min (qs (k) * dp (k) / dp (m), dtime * (tz (m) - tice) / icpk (m))
                            tz (m) = tz (m) - sink * icpk (m)
                            qs (k) = qs (k) - sink * dp (m) / dp (k)
                            if (zt (k) < zs) then
                                r1 = r1 + sink * dp (m) ! precip as rain
                            else
                                ! qr source here will fall next time step (therefore, can evap)
                                qr (m) = qr (m) + sink
                            endif
                        endif
                        if (qs (k) < qrmin) exit
                    enddo
                endif
            enddo
        endif

        if (do_sedi_w) then
            do k = ks, ke
                dm (k) = dp (k) * (1. + qv (k) + ql (k) + qr (k) + qi (k) + qs (k) + qg (k))
            enddo
        endif

        if (use_ppm) then
            call lagrangian_fall_ppm (ks, ke, zs, ze, zt, dp, qs, s1, m1, mono_prof)
        else
            call implicit_fall (dtm, ks, ke, ze, vts, dp, qs, s1, m1)
        endif

        do k = ks, ke
            m1_sol (k) = m1_sol (k) + m1 (k)
        enddo

        if (do_sedi_w) then
            w1 (ks) = w1 (ks) + m1 (ks) * vts (ks) / dm (ks)
            do k = ks + 1, ke
                w1 (k) = (dm (k) * w1 (k) + m1 (k - 1) * (w1 (k - 1) - vts (k - 1)) + m1 (k) * vts (k)) &
                     / (dm (k) + m1 (k - 1))
            enddo
        endif

    endif

    ! ----------------------------------------------
    ! melting of falling graupel into rain
    ! ----------------------------------------------

    call check_column (ks, ke, qg, no_fall)

    if (no_fall) then
        g1 = 0.
    else

        do k = ks + 1, ke
            zt (k) = ze (k) - dt5 * (vtg (k - 1) + vtg (k))
        enddo
        zt (ke + 1) = zs - dtm * vtg (ke)

        do k = ks, ke
            if (zt (k + 1) >= zt (k)) zt (k + 1) = zt (k) - dz_min
        enddo

        if (k0 < ke) then
            do k = ke - 1, k0, - 1
                if (qg (k) > qrmin) then
                    do m = k + 1, ke
                        if (zt (k + 1) >= ze (m)) exit
                        dtime = min (dtm, (ze (m) - ze (m + 1)) / vtg (k))
                        if (zt (k) < ze (m + 1) .and. tz (m) > tice) then
                            dtime = min (1., dtime / tau_g2r)
                            sink = min (qg (k) * dp (k) / dp (m), dtime * (tz (m) - tice) / icpk (m))
                            tz (m) = tz (m) - sink * icpk (m)
                            qg (k) = qg (k) - sink * dp (m) / dp (k)
                            if (zt (k) < zs) then
                                r1 = r1 + sink * dp (m)
                            else
                                qr (m) = qr (m) + sink
                            endif
                        endif
                        if (qg (k) < qrmin) exit
                    enddo
                endif
            enddo
        endif

        if (do_sedi_w) then
            do k = ks, ke
                dm (k) = dp (k) * (1. + qv (k) + ql (k) + qr (k) + qi (k) + qs (k) + qg (k))
            enddo
        endif

        if (use_ppm) then
            call lagrangian_fall_ppm (ks, ke, zs, ze, zt, dp, qg, g1, m1, mono_prof)
        else
            call implicit_fall (dtm, ks, ke, ze, vtg, dp, qg, g1, m1)
        endif

        do k = ks, ke
            m1_sol (k) = m1_sol (k) + m1 (k)
        enddo

        if (do_sedi_w) then
            w1 (ks) = w1 (ks) + m1 (ks) * vtg (ks) / dm (ks)
            do k = ks + 1, ke
                w1 (k) = (dm (k) * w1 (k) + m1 (k - 1) * (w1 (k - 1) - vtg (k - 1)) + m1 (k) * vtg (k)) &
                     / (dm (k) + m1 (k - 1))
            enddo
        endif

    endif

end subroutine terminal_fall

! =======================================================================
! check if water species large enough to fall
! =======================================================================

subroutine check_column (ks, ke, q, no_fall)

    implicit none

    integer, intent (in) :: ks, ke
    real, intent (in) :: q (ks:ke)
    logical, intent (out) :: no_fall
    integer :: k

    no_fall = .true.

    do k = ks, ke
        if (q (k) > qrmin) then
            no_fall = .false.
            exit
        endif
    enddo

end subroutine check_column

! =======================================================================
! time - implicit monotonic scheme
! developed by sj lin, 2016
! =======================================================================

subroutine implicit_fall (dt, ks, ke, ze, vt, dp, q, precip, m1)

    implicit none

    integer, intent (in) :: ks, ke
    real, intent (in) :: dt
    real, intent (in), dimension (ks:ke + 1) :: ze
    real, intent (in), dimension (ks:ke) :: vt, dp
    real, intent (inout), dimension (ks:ke) :: q
    real, intent (out), dimension (ks:ke) :: m1
    real, intent (out) :: precip
    real, dimension (ks:ke) :: dz, qm, dd
    integer :: k

    do k = ks, ke
        dz (k) = ze (k) - ze (k + 1)
        dd (k) = dt * vt (k)
        q (k) = q (k) * dp (k)
    enddo

    ! -----------------------------------------------------------------------
    ! sedimentation: non - vectorizable loop
    ! -----------------------------------------------------------------------

    qm (ks) = q (ks) / (dz (ks) + dd (ks))
    do k = ks + 1, ke
        qm (k) = (q (k) + dd (k - 1) * qm (k - 1)) / (dz (k) + dd (k))
    enddo

    ! -----------------------------------------------------------------------
    ! qm is density at this stage
    ! -----------------------------------------------------------------------

    do k = ks, ke
        qm (k) = qm (k) * dz (k)
    enddo

    ! -----------------------------------------------------------------------
    ! output mass fluxes: non - vectorizable loop
    ! -----------------------------------------------------------------------

    m1 (ks) = q (ks) - qm (ks)
    do k = ks + 1, ke
        m1 (k) = m1 (k - 1) + q (k) - qm (k)
    enddo
    precip = m1 (ke)

    ! -----------------------------------------------------------------------
    ! update:
    ! -----------------------------------------------------------------------

    do k = ks, ke
        q (k) = qm (k) / dp (k)
    enddo

end subroutine implicit_fall

! =======================================================================
! lagrangian scheme
! developed by sj lin, around 2006
! =======================================================================

subroutine lagrangian_fall_ppm (ks, ke, zs, ze, zt, dp, q, precip, m1, mono)

    implicit none

    integer, intent (in) :: ks, ke
    real, intent (in) :: zs
    logical, intent (in) :: mono
    real, intent (in), dimension (ks:ke + 1) :: ze, zt
    real, intent (in), dimension (ks:ke) :: dp

    ! m1: flux
    real, intent (inout), dimension (ks:ke) :: q, m1
    real, intent (out) :: precip
    real, dimension (ks:ke) :: qm, dz

    real :: a4 (4, ks:ke)
    real :: pl, pr, delz, esl
    integer :: k, k0, n, m
    real, parameter :: r3 = 1. / 3., r23 = 2. / 3.

    ! -----------------------------------------------------------------------
    ! density:
    ! -----------------------------------------------------------------------

    do k = ks, ke
        dz (k) = zt (k) - zt (k + 1) ! note: dz is positive
        q (k) = q (k) * dp (k)
        a4 (1, k) = q (k) / dz (k)
        qm (k) = 0.
    enddo

    ! -----------------------------------------------------------------------
    ! construct vertical profile with zt as coordinate
    ! -----------------------------------------------------------------------

    call cs_profile (a4 (1, ks), dz (ks), ke - ks + 1, mono)

    k0 = ks
    do k = ks, ke
        do n = k0, ke
            if (ze (k) <= zt (n) .and. ze (k) >= zt (n + 1)) then
                pl = (zt (n) - ze (k)) / dz (n)
                if (zt (n + 1) <= ze (k + 1)) then
                    ! entire new grid is within the original grid
                    pr = (zt (n) - ze (k + 1)) / dz (n)
                    qm (k) = a4 (2, n) + 0.5 * (a4 (4, n) + a4 (3, n) - a4 (2, n)) * (pr + pl) - &
                        a4 (4, n) * r3 * (pr * (pr + pl) + pl ** 2)
                    qm (k) = qm (k) * (ze (k) - ze (k + 1))
                    k0 = n
                    goto 555
                else
                    qm (k) = (ze (k) - zt (n + 1)) * (a4 (2, n) + 0.5 * (a4 (4, n) + &
                        a4 (3, n) - a4 (2, n)) * (1. + pl) - a4 (4, n) * (r3 * (1. + pl * (1. + pl))))
                    if (n < ke) then
                        do m = n + 1, ke
                            ! locate the bottom edge: ze (k + 1)
                            if (ze (k + 1) < zt (m + 1)) then
                                qm (k) = qm (k) + q (m)
                            else
                                delz = zt (m) - ze (k + 1)
                                esl = delz / dz (m)
                                qm (k) = qm (k) + delz * (a4 (2, m) + 0.5 * esl * &
                                     (a4 (3, m) - a4 (2, m) + a4 (4, m) * (1. - r23 * esl)))
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

    m1 (ks) = q (ks) - qm (ks)
    do k = ks + 1, ke
        m1 (k) = m1 (k - 1) + q (k) - qm (k)
    enddo
    precip = m1 (ke)

    ! convert back to * dry * mixing ratio:
    ! dp must be dry air_mass (because moist air mass will be changed due to terminal fall) .

    do k = ks, ke
        q (k) = qm (k) / dp (k)
    enddo

end subroutine lagrangian_fall_ppm

subroutine cs_profile (a4, del, km, do_mono)

    implicit none

    integer, intent (in) :: km ! vertical dimension
    real, intent (in) :: del (km)
    logical, intent (in) :: do_mono
    real, intent (inout) :: a4 (4, km)
    real, parameter :: qp_min = 1.e-6
    real :: gam (km)
    real :: q (km + 1)
    real :: d4, bet, a_bot, grat, pmp, lac
    real :: pmp_1, lac_1, pmp_2, lac_2
    real :: da1, da2, a6da

    integer :: k

    logical extm (km)

    grat = del (2) / del (1) ! grid ratio
    bet = grat * (grat + 0.5)
    q (1) = (2. * grat * (grat + 1.) * a4 (1, 1) + a4 (1, 2)) / bet
    gam (1) = (1. + grat * (grat + 1.5)) / bet

    do k = 2, km
        d4 = del (k - 1) / del (k)
        bet = 2. + 2. * d4 - gam (k - 1)
        q (k) = (3. * (a4 (1, k - 1) + d4 * a4 (1, k)) - q (k - 1)) / bet
        gam (k) = d4 / bet
    enddo

    a_bot = 1. + d4 * (d4 + 1.5)
    q (km + 1) = (2. * d4 * (d4 + 1.) * a4 (1, km) + a4 (1, km - 1) - a_bot * q (km)) &
         / (d4 * (d4 + 0.5) - a_bot * gam (km))

    do k = km, 1, - 1
        q (k) = q (k) - gam (k) * q (k + 1)
    enddo

    ! -----------------------------------------------------------------------
    ! apply constraints
    ! -----------------------------------------------------------------------

    do k = 2, km
        gam (k) = a4 (1, k) - a4 (1, k - 1)
    enddo

    ! -----------------------------------------------------------------------
    ! apply large - scale constraints to all fields if not local max / min
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! top:
    ! -----------------------------------------------------------------------

    q (1) = max (q (1), 0.)
    q (2) = min (q (2), max (a4 (1, 1), a4 (1, 2)))
    q (2) = max (q (2), min (a4 (1, 1), a4 (1, 2)), 0.)

    ! -----------------------------------------------------------------------
    ! interior:
    ! -----------------------------------------------------------------------

    do k = 3, km - 1
        if (gam (k - 1) * gam (k + 1) > 0.) then
            q (k) = min (q (k), max (a4 (1, k - 1), a4 (1, k)))
            q (k) = max (q (k), min (a4 (1, k - 1), a4 (1, k)))
        else
            if (gam (k - 1) > 0.) then
                ! there exists a local max
                q (k) = max (q (k), min (a4 (1, k - 1), a4 (1, k)))
            else
                ! there exists a local min
                q (k) = min (q (k), max (a4 (1, k - 1), a4 (1, k)))
                q (k) = max (q (k), 0.0)
            endif
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! bottom :
    ! -----------------------------------------------------------------------

    q (km) = min (q (km), max (a4 (1, km - 1), a4 (1, km)))
    q (km) = max (q (km), min (a4 (1, km - 1), a4 (1, km)), 0.)
    ! q (km + 1) = max (q (km + 1), 0.)

    ! -----------------------------------------------------------------------
    ! f (s) = al + s * [ (ar - al) + a6 * (1 - s) ] (0 <= s <= 1)
    ! -----------------------------------------------------------------------

    do k = 1, km - 1
        a4 (2, k) = q (k)
        a4 (3, k) = q (k + 1)
    enddo

    do k = 2, km - 1
        if (gam (k) * gam (k + 1) > 0.0) then
            extm (k) = .false.
        else
            extm (k) = .true.
        endif
    enddo

    if (do_mono) then
        do k = 3, km - 2
            if (extm (k)) then
                ! positive definite constraint only if true local extrema
                if (a4 (1, k) < qp_min .or. extm (k - 1) .or. extm (k + 1)) then
                    a4 (2, k) = a4 (1, k)
                    a4 (3, k) = a4 (1, k)
                endif
            else
                a4 (4, k) = 6. * a4 (1, k) - 3. * (a4 (2, k) + a4 (3, k))
                if (abs (a4 (4, k)) > abs (a4 (2, k) - a4 (3, k))) then
                    ! check within the smooth region if subgrid profile is non - monotonic
                    pmp_1 = a4 (1, k) - 2.0 * gam (k + 1)
                    lac_1 = pmp_1 + 1.5 * gam (k + 2)
                    a4 (2, k) = min (max (a4 (2, k), min (a4 (1, k), pmp_1, lac_1)), &
                        max (a4 (1, k), pmp_1, lac_1))
                    pmp_2 = a4 (1, k) + 2.0 * gam (k)
                    lac_2 = pmp_2 - 1.5 * gam (k - 1)
                    a4 (3, k) = min (max (a4 (3, k), min (a4 (1, k), pmp_2, lac_2)), &
                        max (a4 (1, k), pmp_2, lac_2))
                endif
            endif
        enddo
    else
        do k = 3, km - 2
            if (extm (k)) then
                if (a4 (1, k) < qp_min .or. extm (k - 1) .or. extm (k + 1)) then
                    a4 (2, k) = a4 (1, k)
                    a4 (3, k) = a4 (1, k)
                endif
            endif
        enddo
    endif

    do k = 1, km - 1
        a4 (4, k) = 6. * a4 (1, k) - 3. * (a4 (2, k) + a4 (3, k))
    enddo

    k = km - 1
    if (extm (k)) then
        a4 (2, k) = a4 (1, k)
        a4 (3, k) = a4 (1, k)
        a4 (4, k) = 0.
    else
        da1 = a4 (3, k) - a4 (2, k)
        da2 = da1 ** 2
        a6da = a4 (4, k) * da1
        if (a6da < - da2) then
            a4 (4, k) = 3. * (a4 (2, k) - a4 (1, k))
            a4 (3, k) = a4 (2, k) - a4 (4, k)
        elseif (a6da > da2) then
            a4 (4, k) = 3. * (a4 (3, k) - a4 (1, k))
            a4 (2, k) = a4 (3, k) - a4 (4, k)
        endif
    endif

    call cs_limiters (km - 1, a4)

    ! -----------------------------------------------------------------------
    ! bottom layer:
    ! -----------------------------------------------------------------------

    a4 (2, km) = a4 (1, km)
    a4 (3, km) = a4 (1, km)
    a4 (4, km) = 0.

end subroutine cs_profile

subroutine cs_limiters (km, a4)

    implicit none

    integer, intent (in) :: km

    real, intent (inout) :: a4 (4, km) ! ppm array

    real, parameter :: r12 = 1. / 12.

    integer :: k

    ! -----------------------------------------------------------------------
    ! positive definite constraint
    ! -----------------------------------------------------------------------

    do k = 1, km
        if (abs (a4 (3, k) - a4 (2, k)) < - a4 (4, k)) then
            if ((a4 (1, k) + 0.25 * (a4 (3, k) - a4 (2, k)) ** 2 / a4 (4, k) + a4 (4, k) * r12) < 0.) then
                if (a4 (1, k) < a4 (3, k) .and. a4 (1, k) < a4 (2, k)) then
                    a4 (3, k) = a4 (1, k)
                    a4 (2, k) = a4 (1, k)
                    a4 (4, k) = 0.
                elseif (a4 (3, k) > a4 (2, k)) then
                    a4 (4, k) = 3. * (a4 (2, k) - a4 (1, k))
                    a4 (3, k) = a4 (2, k) - a4 (4, k)
                else
                    a4 (4, k) = 3. * (a4 (3, k) - a4 (1, k))
                    a4 (2, k) = a4 (3, k) - a4 (4, k)
                endif
            endif
        endif
    enddo

end subroutine cs_limiters

! =======================================================================
! calculation of vertical fall speed
! =======================================================================

subroutine fall_speed (ks, ke, den, qs, qi, qg, ql, tk, vts, vti, vtg)

    implicit none

    integer, intent (in) :: ks, ke

    real (kind = r_grid), intent (in), dimension (ks:ke) :: tk
    real, intent (in), dimension (ks:ke) :: den, qs, qi, qg, ql
    real, intent (out), dimension (ks:ke) :: vts, vti, vtg

    ! fall velocity constants:

    real, parameter :: thi = 1.0e-8 ! cloud ice threshold for terminal fall
    real, parameter :: thg = 1.0e-8
    real, parameter :: ths = 1.0e-8

    real, parameter :: aa = - 4.14122e-5
    real, parameter :: bb = - 0.00538922
    real, parameter :: cc = - 0.0516344
    real, parameter :: dd = 0.00216078
    real, parameter :: ee = 1.9714

    ! marshall - palmer constants

    real, parameter :: vcons = 6.6280504
    real, parameter :: vcong = 87.2382675
    real, parameter :: vconh = vcong * sqrt (rhoh / rhog) ! 132.087495104005
    real, parameter :: norms = 942477796.076938
    real, parameter :: normg = 5026548245.74367
    real, parameter :: normh = pi * rhoh * rnzh ! 115233618.533674

    real, dimension (ks:ke) :: qden, tc, rhof

    real :: vi0

    integer :: k

    ! -----------------------------------------------------------------------
    ! marshall - palmer formula
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! try the local air density -- for global model; the true value could be
    ! much smaller than sfcrho over high mountains
    ! -----------------------------------------------------------------------

    do k = ks, ke
        rhof (k) = sqrt (min (10., sfcrho / den (k)))
    enddo

    ! -----------------------------------------------------------------------
    ! ice:
    ! -----------------------------------------------------------------------

    if (const_vi) then
        vti (:) = vi_fac
    else
        ! -----------------------------------------------------------------------
        ! use deng and mace (2008, grl), which gives smaller fall speed than hd90 formula
        ! -----------------------------------------------------------------------
        vi0 = 0.01 * vi_fac
        do k = ks, ke
            if (qi (k) < thi) then ! this is needed as the fall - speed maybe problematic for small qi
                vti (k) = vf_min
            else
                tc (k) = tk (k) - tice
                vti (k) = (3. + log10 (qi (k) * den (k))) * (tc (k) * (aa * tc (k) + bb) + cc) + dd * tc (k) + ee
                vti (k) = vi0 * exp (log_10 * vti (k))
                vti (k) = min (vi_max, max (vf_min, vti (k)))
            endif
        enddo
    endif

    ! -----------------------------------------------------------------------
    ! snow:
    ! -----------------------------------------------------------------------

    if (const_vs) then
        vts (:) = vs_fac ! 1. ifs_2016
    else
        do k = ks, ke
            if (qs (k) < ths) then
                vts (k) = vf_min
            else
                vts (k) = vs_fac * vcons * rhof (k) * exp (0.0625 * log (qs (k) * den (k) / norms))
                vts (k) = min (vs_max, max (vf_min, vts (k)))
            endif
        enddo
    endif

    ! -----------------------------------------------------------------------
    ! graupel:
    ! -----------------------------------------------------------------------

    if (const_vg) then
        vtg (:) = vg_fac ! 2.
    else
        if (do_hail) then
            do k = ks, ke
                if (qg (k) < thg) then
                    vtg (k) = vf_min
                else
                    vtg (k) = vg_fac * vconh * rhof (k) * sqrt (sqrt (sqrt (qg (k) * den (k) / normh)))
                    vtg (k) = min (vg_max, max (vf_min, vtg (k)))
                endif
            enddo
        else
            do k = ks, ke
                if (qg (k) < thg) then
                    vtg (k) = vf_min
                else
                    vtg (k) = vg_fac * vcong * rhof (k) * sqrt (sqrt (sqrt (qg (k) * den (k) / normg)))
                    vtg (k) = min (vg_max, max (vf_min, vtg (k)))
                endif
            enddo
        endif
    endif

end subroutine fall_speed

! =======================================================================
! setup gfdl cloud microphysics parameters
! =======================================================================

subroutine setupm

    implicit none

    real :: gcon, cd, scm3, pisq, act (8)
    real :: vdifu, tcond
    real :: visk
    real :: ch2o, hltf
    real :: hlts, hltc, ri50

    real, parameter :: gam263 = 1.456943, gam275 = 1.608355, gam290 = 1.827363, &
        gam325 = 2.54925, gam350 = 3.323363, gam380 = 4.694155, &
        gam425 = 8.285063, gam450 = 11.631769, gam480 = 17.837789, &
        gam625 = 184.860962, gam680 = 496.604067

    real, parameter :: acc (3) = (/ 5.0, 2.0, 0.5 /)

    real den_rc

    integer :: i, k

    pie = 4. * atan (1.0)

    ! s. klein's formular (eq 16) from am2

    fac_rc = (4. / 3.) * pie * rhor * rthresh ** 3

    if (prog_ccn) then
        ! if (master) write (*, *) 'prog_ccn option is .t.'
    else
        den_rc = fac_rc * ccn_o * 1.e6
        ! if (master) write (*, *) 'mp: for ccn_o = ', ccn_o, 'ql_rc = ', den_rc
        den_rc = fac_rc * ccn_l * 1.e6
        ! if (master) write (*, *) 'mp: for ccn_l = ', ccn_l, 'ql_rc = ', den_rc
    endif

    vdifu = 2.11e-5
    tcond = 2.36e-2

    visk = 1.259e-5
    hlts = 2.8336e6
    hltc = 2.5e6
    hltf = 3.336e5

    ch2o = 4.1855e3
    ri50 = 1.e-4

    pisq = pie * pie
    scm3 = (visk / vdifu) ** (1. / 3.)

    cracs = pisq * rnzr * rnzs * rhos
    csacr = pisq * rnzr * rnzs * rhor
    if (do_hail) then
        cgacr = pisq * rnzr * rnzh * rhor
        cgacs = pisq * rnzh * rnzs * rhos
    else
        cgacr = pisq * rnzr * rnzg * rhor
        cgacs = pisq * rnzg * rnzs * rhos
    endif
    cgacs = cgacs * c_pgacs

    ! act: 1 - 2:racs (s - r) ; 3 - 4:sacr (r - s) ;
    ! 5 - 6:gacr (r - g) ; 7 - 8:gacs (s - g)

    act (1) = pie * rnzs * rhos
    act (2) = pie * rnzr * rhor
    if (do_hail) then
        act (6) = pie * rnzh * rhoh
    else
        act (6) = pie * rnzg * rhog
    endif
    act (3) = act (2)
    act (4) = act (1)
    act (5) = act (2)
    act (7) = act (1)
    act (8) = act (6)

    do i = 1, 3
        do k = 1, 4
            acco (i, k) = acc (i) / (act (2 * k - 1) ** ((7 - i) * 0.25) * act (2 * k) ** (i * 0.25))
        enddo
    enddo

    gcon = 40.74 * sqrt (sfcrho) ! 44.628

    csacw = pie * rnzs * clin * gam325 / (4. * act (1) ** 0.8125)
    ! decreasing csacw to reduce cloud water --- > snow

    craci = pie * rnzr * alin * gam380 / (4. * act (2) ** 0.95)
    csaci = csacw * c_psaci

    if (do_hail) then
        cgacw = pie * rnzh * gam350 * gcon / (4. * act (6) ** 0.875)
    else
        cgacw = pie * rnzg * gam350 * gcon / (4. * act (6) ** 0.875)
    endif
    ! cgaci = cgacw * 0.1

    ! sjl, may 28, 2012
    cgaci = cgacw * 0.05
    ! sjl, may 28, 2012

    cracw = craci ! cracw = 3.27206196043822
    cracw = c_cracw * cracw

    ! subl and revp: five constants for three separate processes

    cssub (1) = 2. * pie * vdifu * tcond * rvgas * rnzs
    if (do_hail) then
        cgsub (1) = 2. * pie * vdifu * tcond * rvgas * rnzh
    else
        cgsub (1) = 2. * pie * vdifu * tcond * rvgas * rnzg
    endif
    crevp (1) = 2. * pie * vdifu * tcond * rvgas * rnzr
    cssub (2) = 0.78 / sqrt (act (1))
    cgsub (2) = 0.78 / sqrt (act (6))
    crevp (2) = 0.78 / sqrt (act (2))
    cssub (3) = 0.31 * scm3 * gam263 * sqrt (clin / visk) / act (1) ** 0.65625
    cgsub (3) = 0.31 * scm3 * gam275 * sqrt (gcon / visk) / act (6) ** 0.6875
    crevp (3) = 0.31 * scm3 * gam290 * sqrt (alin / visk) / act (2) ** 0.725
    cssub (4) = tcond * rvgas
    cssub (5) = hlts ** 2 * vdifu
    cgsub (4) = cssub (4)
    crevp (4) = cssub (4)
    cgsub (5) = cssub (5)
    crevp (5) = hltc ** 2 * vdifu

    cgfr (1) = 20.e2 * pisq * rnzr * rhor / act (2) ** 1.75
    cgfr (2) = 0.66

    ! smlt: five constants (lin et al. 1983)

    csmlt (1) = 2. * pie * tcond * rnzs / hltf
    csmlt (2) = 2. * pie * vdifu * rnzs * hltc / hltf
    csmlt (3) = cssub (2)
    csmlt (4) = cssub (3)
    csmlt (5) = ch2o / hltf

    ! gmlt: five constants

    if (do_hail) then
        cgmlt (1) = 2. * pie * tcond * rnzh / hltf
        cgmlt (2) = 2. * pie * vdifu * rnzh * hltc / hltf
    else
        cgmlt (1) = 2. * pie * tcond * rnzg / hltf
        cgmlt (2) = 2. * pie * vdifu * rnzg * hltc / hltf
    endif
    cgmlt (3) = cgsub (2)
    cgmlt (4) = cgsub (3)
    cgmlt (5) = ch2o / hltf

    es0 = 6.107799961e2 ! ~6.1 mb
    ces0 = eps * es0

end subroutine setupm

! =======================================================================
! initialization of gfdl cloud microphysics
! =======================================================================

!subroutine gfdl_mp_init (id, jd, kd, axes, time)
subroutine gfdl_mp_init (me, master, nlunit, input_nml_file, logunit, fn_nml)

    implicit none

    integer, intent (in) :: me
    integer, intent (in) :: master
    integer, intent (in) :: nlunit
    integer, intent (in) :: logunit

    character (len = 64), intent (in) :: fn_nml
    character (len = *), intent (in) :: input_nml_file (:)

    integer :: ios
    logical :: exists

    ! integer, intent (in) :: id, jd, kd
    ! integer, intent (in) :: axes (4)
    ! type (time_type), intent (in) :: time

    ! integer :: unit, io, ierr, k, logunit
    ! logical :: flag
    ! real :: tmp, q1, q2

    ! master = (mpp_pe () .eq.mpp_root_pe ())

    !#ifdef internal_file_nml
    ! read (input_nml_file, nml = gfdl_mp_nml, iostat = io)
    ! ierr = check_nml_error (io, 'gfdl_mp_nml')
    !#else
    ! if (file_exist ('input.nml')) then
    ! unit = open_namelist_file ()
    ! io = 1
    ! do while (io .ne. 0)
    ! read (unit, nml = gfdl_mp_nml, iostat = io, end = 10)
    ! ierr = check_nml_error (io, 'gfdl_mp_nml')
    ! enddo
    !10 call close_file (unit)
    ! endif
    !#endif
    ! call write_version_number ('gfdl_mp_mod', version)
    ! logunit = stdlog ()

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml = gfdl_mp_nml)
#else
    inquire (file = trim (fn_nml), exist = exists)
    if (.not. exists) then
        write (6, *) 'gfdl - mp :: namelist file: ', trim (fn_nml), ' does not exist'
        stop
    else
        open (unit = nlunit, file = fn_nml, readonly, status = 'old', iostat = ios)
    endif
    rewind (nlunit)
    read (nlunit, nml = gfdl_mp_nml)
    close (nlunit)
#endif

    ! write version number and namelist to log file

    if (me == master) then
        write (logunit, *) " ================================================================== "
        write (logunit, *) "gfdl_mp_mod"
        write (logunit, nml = gfdl_mp_nml)
    endif

    if (do_setup) then
        call setup_con
        call setupm
        do_setup = .false.
    endif

    g2 = 0.5 * grav
    log_10 = log (10.)

    tice0 = tice - 0.01
    t_wfr = tice - 40.0 ! supercooled water can exist down to - 48 c, which is the "absolute"

    ! if (master) write (logunit, nml = gfdl_mp_nml)

    ! if (master) write (*, *) 'prec_lin diagnostics initialized.', id_prec

    ! call qsmith_init

    ! testing the water vapor tables

    ! if (mp_debug .and. master) then
    ! write (*, *) 'testing water vapor tables in gfdl_mp'
    ! tmp = tice - 90.
    ! do k = 1, 25
    ! q1 = wqsat_moist (tmp, 0., 1.e5)
    ! q2 = qs1d_m (tmp, 0., 1.e5)
    ! write (*, *) nint (tmp - tice), q1, q2, 'dq = ', q1 - q2
    ! tmp = tmp + 5.
    ! enddo
    ! endif

    ! if (master) write (*, *) 'gfdl_cloud_micrphys diagnostics initialized.'

    ! gfdl_mp_clock = mpp_clock_id ('gfdl_mp', grain = clock_routine)

    module_is_initialized = .true.

end subroutine gfdl_mp_init

! =======================================================================
! end of gfdl cloud microphysics
! =======================================================================

subroutine gfdl_mp_end

    implicit none

    deallocate (table)
    deallocate (table2)
    deallocate (table3)
    deallocate (tablew)
    deallocate (des)
    deallocate (des2)
    deallocate (des3)
    deallocate (desw)

    tables_are_initialized = .false.

end subroutine gfdl_mp_end

! =======================================================================
! qsmith table initialization
! =======================================================================

subroutine setup_con

    implicit none

    ! master = (mpp_pe () .eq.mpp_root_pe ())

    rgrav = 1. / grav

    if (.not. qsmith_tables_initialized) call qsmith_init

    qsmith_tables_initialized = .true.

end subroutine setup_con

! =======================================================================
! accretion function (lin et al. 1983)
! =======================================================================

real function acr3d (v1, v2, q1, q2, c, cac, rho)

    implicit none

    real, intent (in) :: v1, v2, c, rho
    real, intent (in) :: q1, q2 ! mixing ratio!!!
    real, intent (in) :: cac (3)

    real :: t1, s1, s2

    ! integer :: k
    !
    ! real :: a
    !
    ! a = 0.0
    ! do k = 1, 3
    ! a = a + cac (k) * ((q1 * rho) ** ((7 - k) * 0.25) * (q2 * rho) ** (k * 0.25))
    ! enddo
    ! acr3d = c * abs (v1 - v2) * a / rho

    ! optimized

    t1 = sqrt (q1 * rho)
    s1 = sqrt (q2 * rho)
    s2 = sqrt (s1) ! s1 = s2 ** 2
    acr3d = c * abs (v1 - v2) * q1 * s2 * (cac (1) * t1 + cac (2) * sqrt (t1) * s2 + cac (3) * s1)

end function acr3d

! =======================================================================
! melting of snow function (lin et al. 1983)
! note: psacw and psacr must be calc before smlt is called
! =======================================================================

real function smlt (tc, dqs, qsrho, psacw, psacr, c, rho, rhofac)

    implicit none

    real, intent (in) :: tc, dqs, qsrho, psacw, psacr, c (5), rho, rhofac

    smlt = (c (1) * tc / rho - c (2) * dqs) * (c (3) * sqrt (qsrho) + &
        c (4) * qsrho ** 0.65625 * sqrt (rhofac)) + c (5) * tc * (psacw + psacr)

end function smlt

! =======================================================================
! melting of graupel function (lin et al. 1983)
! note: pgacw and pgacr must be calc before gmlt is called
! =======================================================================

real function gmlt (tc, dqs, qgrho, pgacw, pgacr, c, rho)

    implicit none

    real, intent (in) :: tc, dqs, qgrho, pgacw, pgacr, c (5), rho

    gmlt = (c (1) * tc / rho - c (2) * dqs) * (c (3) * sqrt (qgrho) + &
        c (4) * qgrho ** 0.6875 / rho ** 0.25) + c (5) * tc * (pgacw + pgacr)

end function gmlt

! =======================================================================
! initialization
! prepare saturation water vapor pressure tables
! =======================================================================

subroutine qsmith_init

    implicit none

    integer, parameter :: length = 2621

    integer :: i

    if (.not. tables_are_initialized) then

        ! master = (mpp_pe () .eq. mpp_root_pe ())
        ! if (master) print *, ' gfdl mp: initializing qs tables'

        ! debug code
        ! print *, mpp_pe (), allocated (table), allocated (table2), &
        ! allocated (table3), allocated (tablew), allocated (des), &
        ! allocated (des2), allocated (des3), allocated (desw)
        ! end debug code

        ! generate es table (dt = 0.1 deg. c)

        allocate (table (length))
        allocate (table2 (length))
        allocate (table3 (length))
        allocate (tablew (length))
        allocate (des (length))
        allocate (des2 (length))
        allocate (des3 (length))
        allocate (desw (length))

        call qs_table (length)
        call qs_table2 (length)
        call qs_table3 (length)
        call qs_tablew (length)

        do i = 1, length - 1
            des (i) = max (0., table (i + 1) - table (i))
            des2 (i) = max (0., table2 (i + 1) - table2 (i))
            des3 (i) = max (0., table3 (i + 1) - table3 (i))
            desw (i) = max (0., tablew (i + 1) - tablew (i))
        enddo
        des (length) = des (length - 1)
        des2 (length) = des2 (length - 1)
        des3 (length) = des3 (length - 1)
        desw (length) = desw (length - 1)

        tables_are_initialized = .true.

    endif

end subroutine qsmith_init

! =======================================================================
! compute the saturated specific humidity for table ii
! =======================================================================

real function wqsat_moist (ta, qv, pa)

    implicit none

    real, intent (in) :: ta, pa, qv

    real :: es, ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = tablew (it) + (ap1 - it) * desw (it)
    wqsat_moist = eps * es * (1. + zvir * qv) / pa

end function wqsat_moist


real function wqs1 (ta, den)

    implicit none

    ! pure water phase; universal dry / moist formular using air density
    ! input "den" can be either dry or moist air density

    real, intent (in) :: ta, den

    real :: es, ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    !NOTE: a crash here usually means NaN
    !if (it < 1 .or. it > 2621) then
    !   write(*,*), 'WQS1: table range violation', it, ta, tmin, den
    !endif
    es = tablew (it) + (ap1 - it) * desw (it)
    wqs1 = es / (rvgas * ta * den)

end function wqs1

! =======================================================================
! compute the gradient of saturated specific humidity for table ii
! =======================================================================

real function wqs2 (ta, den, dqdt)

    implicit none

    ! pure water phase; universal dry / moist formular using air density
    ! input "den" can be either dry or moist air density

    real, intent (in) :: ta, den
    real, intent (out) :: dqdt
    real :: es, ap1, tmin
    integer :: it

    tmin = table_ice - 160.

    if (.not. tables_are_initialized) call qsmith_init

    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    !NOTE: a crash here usually means NaN
    !if (it < 1 .or. it > 2621) then
    !   write(*,*), 'WQS2: table range violation', it, ta, tmin, den
    !endif
    es = tablew (it) + (ap1 - it) * desw (it)
    wqs2 = es / (rvgas * ta * den)
    it = ap1 - 0.5
    ! finite diff, del_t = 0.1:
    dqdt = 10. * (desw (it) + (ap1 - it) * (desw (it + 1) - desw (it))) / (rvgas * ta * den)

end function wqs2

! =======================================================================
! compute wet buld temperature
! =======================================================================

real function wet_bulb (q, t, den)

    implicit none

    real, intent (in) :: t, q, den

    real :: qs, tp, dqdt

    wet_bulb = t
    qs = wqs2 (wet_bulb, den, dqdt)
    tp = 0.5 * (qs - q) / (1. + lcp * dqdt) * lcp
    wet_bulb = wet_bulb - tp

    ! tp is negative if super - saturated
    if (tp > 0.01) then
        qs = wqs2 (wet_bulb, den, dqdt)
        tp = (qs - q) / (1. + lcp * dqdt) * lcp
        wet_bulb = wet_bulb - tp
    endif

end function wet_bulb

! =======================================================================
! compute the saturated specific humidity for table iii
! =======================================================================

real function iqs1 (ta, den)

    implicit none

    ! water - ice phase; universal dry / moist formular using air density
    ! input "den" can be either dry or moist air density

    real, intent (in) :: ta, den

    real :: es, ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = table2 (it) + (ap1 - it) * des2 (it)
    iqs1 = es / (rvgas * ta * den)

end function iqs1

! =======================================================================
! compute the gradient of saturated specific humidity for table iii
! =======================================================================

real function iqs2 (ta, den, dqdt)

    implicit none

    ! water - ice phase; universal dry / moist formular using air density
    ! input "den" can be either dry or moist air density

    real (kind = r_grid), intent (in) :: ta
    real, intent (in) :: den
    real, intent (out) :: dqdt
    real (kind = r_grid) :: tmin, es, ap1
    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = table2 (it) + (ap1 - it) * des2 (it)
    iqs2 = es / (rvgas * ta * den)
    it = ap1 - 0.5
    dqdt = 10. * (des2 (it) + (ap1 - it) * (des2 (it + 1) - des2 (it))) / (rvgas * ta * den)

end function iqs2


! =======================================================================
! saturation water vapor pressure table ii
! 1 - phase table
! =======================================================================

subroutine qs_tablew (n)

    implicit none

    integer, intent (in) :: n

    real (kind = r_grid) :: delt = 0.1
    real (kind = r_grid) :: tmin, tem, fac0, fac1, fac2

    integer :: i

    tmin = table_ice - 160.

    ! -----------------------------------------------------------------------
    ! compute es over water
    ! -----------------------------------------------------------------------

    do i = 1, n
        tem = tmin + delt * real (i - 1)
        fac0 = (tem - t_ice) / (tem * t_ice)
        fac1 = fac0 * lv0
        fac2 = (dc_vap * log (tem / t_ice) + fac1) / rvgas
        tablew (i) = e00 * exp (fac2)
    enddo

end subroutine qs_tablew

! =======================================================================
! saturation water vapor pressure table iii
! 2 - phase table
! =======================================================================

subroutine qs_table2 (n)

    implicit none

    integer, intent (in) :: n

    real (kind = r_grid) :: delt = 0.1
    real (kind = r_grid) :: tmin, tem0, tem1, fac0, fac1, fac2

    integer :: i, i0, i1

    tmin = table_ice - 160.

    do i = 1, n
        tem0 = tmin + delt * real (i - 1)
        fac0 = (tem0 - t_ice) / (tem0 * t_ice)
        if (i <= 1600) then
            ! -----------------------------------------------------------------------
            ! compute es over ice between - 160 deg c and 0 deg c.
            ! -----------------------------------------------------------------------
            fac1 = fac0 * li2
            fac2 = (d2ice * log (tem0 / t_ice) + fac1) / rvgas
        else
            ! -----------------------------------------------------------------------
            ! compute es over water between 0 deg c and 102 deg c.
            ! -----------------------------------------------------------------------
            fac1 = fac0 * lv0
            fac2 = (dc_vap * log (tem0 / t_ice) + fac1) / rvgas
        endif
        table2 (i) = e00 * exp (fac2)
    enddo

    ! -----------------------------------------------------------------------
    ! smoother around 0 deg c
    ! -----------------------------------------------------------------------

    i0 = 1600
    i1 = 1601
    tem0 = 0.25 * (table2 (i0 - 1) + 2. * table (i0) + table2 (i0 + 1))
    tem1 = 0.25 * (table2 (i1 - 1) + 2. * table (i1) + table2 (i1 + 1))
    table2 (i0) = tem0
    table2 (i1) = tem1

end subroutine qs_table2

! =======================================================================
! saturation water vapor pressure table iv
! 2 - phase table with " - 2 c" as the transition point
! =======================================================================

subroutine qs_table3 (n)

    implicit none

    integer, intent (in) :: n

    real (kind = r_grid) :: delt = 0.1
    real (kind = r_grid) :: esbasw, tbasw, esbasi, tmin, tem, aa, b, c, d, e
    real (kind = r_grid) :: tem0, tem1

    integer :: i, i0, i1

    esbasw = 1013246.0
    tbasw = table_ice + 100.
    esbasi = 6107.1
    tmin = table_ice - 160.

    do i = 1, n
        tem = tmin + delt * real (i - 1)
        ! if (i <= 1600) then
        if (i <= 1580) then ! change to - 2 c
            ! -----------------------------------------------------------------------
            ! compute es over ice between - 160 deg c and 0 deg c.
            ! see smithsonian meteorological tables page 350.
            ! -----------------------------------------------------------------------
            aa = - 9.09718 * (table_ice / tem - 1.)
            b = - 3.56654 * log10 (table_ice / tem)
            c = 0.876793 * (1. - tem / table_ice)
            e = log10 (esbasi)
            table3 (i) = 0.1 * 10 ** (aa + b + c + e)
        else
            ! -----------------------------------------------------------------------
            ! compute es over water between - 2 deg c and 102 deg c.
            ! see smithsonian meteorological tables page 350.
            ! -----------------------------------------------------------------------
            aa = - 7.90298 * (tbasw / tem - 1.)
            b = 5.02808 * log10 (tbasw / tem)
            c = - 1.3816e-7 * (10 ** ((1. - tem / tbasw) * 11.344) - 1.)
            d = 8.1328e-3 * (10 ** ((tbasw / tem - 1.) * (- 3.49149)) - 1.)
            e = log10 (esbasw)
            table3 (i) = 0.1 * 10 ** (aa + b + c + d + e)
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! smoother around - 2 deg c
    ! -----------------------------------------------------------------------

    i0 = 1580
    i1 = 1581
    tem0 = 0.25 * (table3 (i0 - 1) + 2. * table (i0) + table3 (i0 + 1))
    tem1 = 0.25 * (table3 (i1 - 1) + 2. * table (i1) + table3 (i1 + 1))
    table3 (i0) = tem0
    table3 (i1) = tem1

end subroutine qs_table3


! =======================================================================
! saturation water vapor pressure table i
! 3 - phase table
! =======================================================================

subroutine qs_table (n)

    implicit none

    integer, intent (in) :: n

    real (kind = r_grid) :: delt = 0.1
    real (kind = r_grid) :: tmin, tem, esh20
    real (kind = r_grid) :: wice, wh2o, fac0, fac1, fac2
    real (kind = r_grid) :: esupc (200)

    integer :: i

    tmin = table_ice - 160.

    ! -----------------------------------------------------------------------
    ! compute es over ice between - 160 deg c and 0 deg c.
    ! -----------------------------------------------------------------------

    do i = 1, 1600
        tem = tmin + delt * real (i - 1)
        fac0 = (tem - t_ice) / (tem * t_ice)
        fac1 = fac0 * li2
        fac2 = (d2ice * log (tem / t_ice) + fac1) / rvgas
        table (i) = e00 * exp (fac2)
    enddo

    ! -----------------------------------------------------------------------
    ! compute es over water between - 20 deg c and 102 deg c.
    ! -----------------------------------------------------------------------

    do i = 1, 1221
        tem = 253.16 + delt * real (i - 1)
        fac0 = (tem - t_ice) / (tem * t_ice)
        fac1 = fac0 * lv0
        fac2 = (dc_vap * log (tem / t_ice) + fac1) / rvgas
        esh20 = e00 * exp (fac2)
        if (i <= 200) then
            esupc (i) = esh20
        else
            table (i + 1400) = esh20
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! derive blended es over ice and supercooled water between - 20 deg c and 0 deg c
    ! -----------------------------------------------------------------------

    do i = 1, 200
        tem = 253.16 + delt * real (i - 1)
        wice = 0.05 * (table_ice - tem)
        wh2o = 0.05 * (tem - 253.16)
        table (i + 1400) = wice * table (i + 1400) + wh2o * esupc (i)
    enddo

end subroutine qs_table


! =======================================================================
! fix negative water species
! this is designed for 6 - class micro - physics schemes
! =======================================================================

subroutine neg_adj (ks, ke, pt, dp, qv, ql, qr, qi, qs, qg)

    implicit none

    integer, intent (in) :: ks, ke
    real, intent (in), dimension (ks:ke) :: dp
    real (kind = r_grid), intent (inout), dimension (ks:ke) :: pt
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg

    real, dimension (ks:ke) :: lcpk, icpk

    real :: dq, cvm

    integer :: k

    ! -----------------------------------------------------------------------
    ! define heat capacity and latent heat coefficient
    ! -----------------------------------------------------------------------

    do k = ks, ke
        cvm = 1. + qv (k) * c1_vap + (qr (k) + ql (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
        lcpk (k) = (lv00 + d1_vap * pt (k)) / cvm
        icpk (k) = (li00 + d1_ice * pt (k)) / cvm
    enddo

    do k = ks, ke

        ! -----------------------------------------------------------------------
        ! ice phase:
        ! -----------------------------------------------------------------------

        ! if cloud ice < 0, borrow from snow
        if (qi (k) < 0.) then
            qs (k) = qs (k) + qi (k)
            qi (k) = 0.
        endif
        ! if snow < 0, borrow from graupel
        if (qs (k) < 0.) then
            qg (k) = qg (k) + qs (k)
            qs (k) = 0.
        endif
        ! if graupel < 0, borrow from rain
#ifdef HIGH_NEG_HT
        if (qg (k) < 0.) then
            qr (k) = qr (k) + qg (k)
            pt (k) = pt (k) - qg (k) * icpk (k) ! heating
            qg (k) = 0.
        endif
#endif

        ! -----------------------------------------------------------------------
        ! liquid phase:
        ! -----------------------------------------------------------------------

        ! if rain < 0, borrow from cloud water
        if (qr (k) < 0.) then
            ql (k) = ql (k) + qr (k)
            qr (k) = 0.
        endif

    enddo

end subroutine neg_adj

! =======================================================================
! compute global sum
! quick local sum algorithm
! =======================================================================

!real function g_sum (p, ifirst, ilast, jfirst, jlast, area, mode)
!
! use mpp_mod, only: mpp_sum
!
! implicit none
!
! integer, intent (in) :: ifirst, ilast, jfirst, jlast
! integer, intent (in) :: mode ! if == 1 divided by area
!
! real, intent (in), dimension (ifirst:ilast, jfirst:jlast) :: p, area
!
! integer :: i, j
!
! real :: gsum
!
! if (global_area < 0.) then
! global_area = 0.
! do j = jfirst, jlast
! do i = ifirst, ilast
! global_area = global_area + area (i, j)
! enddo
! enddo
! call mpp_sum (global_area)
! endif
!
! gsum = 0.
! do j = jfirst, jlast
! do i = ifirst, ilast
! gsum = gsum + p (i, j) * area (i, j)
! enddo
! enddo
! call mpp_sum (gsum)
!
! if (mode == 1) then
! g_sum = gsum / global_area
! else
! g_sum = gsum
! endif
!
!end function g_sum

end module gfdl_mp_mod
