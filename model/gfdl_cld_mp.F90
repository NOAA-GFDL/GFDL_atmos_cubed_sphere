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
!* useful, but WITHOUT ANY WARRANTY; without even the implied warranty
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
! =======================================================================

module gfdl_cld_mp_mod

#ifdef GFS_PHYS
    use machine, only: r_grid => kind_phys
#endif

    implicit none

    private

    public gfdl_cld_mp_driver, gfdl_cld_mp_init, gfdl_cld_mp_end
    public wqs1, wqs2, iqs1, iqs2, mpdrv, sedi_heat, warm_rain, revap_racc, &
        linear_prof, icloud, subgrid_z_proc, terminal_fall, check_column, implicit_fall, &
        lagrangian_fall_ppm, cs_profile, cs_limiters, fall_speed, setupm, setup_con, &
        qsmith_init, qs_tablew, qs_table2, qs_table3, qs_table, neg_adj, acr3d, smlt, gmlt, &
        wet_bulb, qsmith, qs_blend, es3_table1d, es2_table1d, esw_table1d, es2_table, &
        esw_table, d_sat, qs1d_m, wqsat_moist, wqsat2_moist, qs1d_moist, revap_rac1, &
        wqs2_vect, rhow, rhor, rhos, rhog, rhoh, rnzr, rnzs, rnzg, rnzh, rvgas, rdgas, &
        grav, hlv, hlf, cp_air, cp_vap, cv_air, cv_vap, c_ice, c_liq, dc_vap, dc_ice, &
        t_ice, t_wfr, e00, pi, zvir, rgrav

#ifndef GFS_PHYS
    integer, parameter :: r_grid = 8
#endif

    logical :: module_is_initialized = .false.
    logical :: qsmith_tables_initialized = .false.

    real, parameter :: grav = 9.80665 ! gfs: acceleration due to gravity
    real, parameter :: rdgas = 287.05 ! gfs: gas constant for dry air
    real, parameter :: rvgas = 461.50 ! gfs: gas constant for water vapor
    real, parameter :: cp_air = 1.0046e3 ! gfs: heat capacity of dry air at constant pressure
    real, parameter :: hlv = 2.5e6 ! gfs: latent heat of evaporation
    real, parameter :: hlf = 3.3358e5 ! gfs: latent heat of fusion
    real, parameter :: pi = 3.1415926535897931 ! gfs: ratio of circle circumference to diameter

    ! real, parameter :: cp_air = rdgas * 7. / 2. ! 1004.675, heat capacity of dry air at constant pressure
    real, parameter :: cp_vap = 4.0 * rvgas ! 1846.0, heat capacity of water vapore at constnat pressure
    ! real, parameter :: cv_air = 717.56 ! satoh value, heat capacity of dry air at constant volume
    real, parameter :: cv_air = cp_air - rdgas ! 717.55, heat capacity of dry air at constant volume
    ! real, parameter :: cv_vap = 1410.0 ! emanuel value, heat capacity of water vapor at constant volume
    real, parameter :: cv_vap = 3.0 * rvgas ! 1384.5, heat capacity of water vapor at constant volume

    ! http: // www.engineeringtoolbox.com / ice - thermal - properties - d_576.html
    ! c_ice = 2050.0 at 0 deg c
    ! c_ice = 2000.0 at - 10 deg c
    ! c_ice = 1943.0 at - 20 deg c
    ! c_ice = 1882.0 at - 30 deg c
    ! c_ice = 1818.0 at - 40 deg c

    ! https: // www.engineeringtoolbox.com / specific - heat - capacity - water - d_660.html
    ! c_liq = 4219.9 at 0.01 deg c
    ! c_liq = 4195.5 at 10 deg c
    ! c_liq = 4184.4 at 20 deg c
    ! c_liq = 4180.1 at 30 deg c
    ! c_liq = 4179.6 at 40 deg c

    ! the following two are from emanuel's book "atmospheric convection"
    ! real, parameter :: c_ice = 2.106e3 ! heat capacity of ice at 0 deg c: c = c_ice + 7.3 * (t - tice)
    ! real, parameter :: c_liq = 4.190e3 ! heat capacity of water at 0 deg c
    ! real, parameter :: c_ice = 1.972e3 ! gfdl: heat capacity of ice at - 15 deg c
    ! real, parameter :: c_liq = 4.1855e3 ! gfdl: heat capacity of water at 15 deg c
    ! real, parameter :: c_ice = 2.106e3 ! gfs: heat capacity of ice at 0 deg c
    ! real, parameter :: c_liq = 4.1855e3 ! gfs: heat capacity of liquid at 15 deg c
    real, parameter :: c_ice = 2.106e3 ! ifs: heat capacity of ice at 0 deg c
    real, parameter :: c_liq = 4.218e3 ! ifs: heat capacity of water at 0 deg c

    real, parameter :: eps = rdgas / rvgas ! 0.6219934995
    real, parameter :: zvir = rvgas / rdgas - 1. ! 0.6077338443

    real, parameter :: dc_vap = cp_vap - c_liq ! - 2.372e3, isobaric heating / cooling
    real, parameter :: dc_ice = c_liq - c_ice ! 2.112e3, isobaric heating / colling

    real, parameter :: t_ice = 273.16 ! freezing temperature
    real, parameter :: table_ice = 273.16 ! freezing point for qs table
    real :: t_wfr ! complete freezing temperature

    real (kind = r_grid), parameter :: e00 = 611.21 ! ifs: saturation vapor pressure at 0 deg c
    ! real (kind = r_grid), parameter :: e00 = 610.71 ! gfdl: saturation vapor pressure at 0 deg c

    real, parameter :: hlv0 = hlv ! gfs: evaporation latent heat coefficient at 0 deg c
    ! real, parameter :: hlv0 = 2.501e6 ! emanuel value
    real, parameter :: hlf0 = hlf ! gfs: fussion latent heat coefficient at 0 deg c
    ! real, parameter :: hlf0 = 3.337e5 ! emanuel value

    real, parameter :: lv0 = hlv0 - dc_vap * t_ice ! 3.14893552e6, evaporation latent heat coefficient at 0 deg k
    real, parameter :: li0 = hlf0 - dc_ice * t_ice ! - 2.2691392e5, fussion latend heat coefficient at 0 deg k

    real (kind = r_grid), parameter :: d2ice = cp_vap - c_ice ! - 260.0, isobaric heating / cooling
    real (kind = r_grid), parameter :: li2 = lv0 + li0 ! 2.9220216e6, sublimation latent heat coefficient at 0 deg k

    real, parameter :: qrmin = 1.e-8 ! min value for cloud condensates
    real, parameter :: qvmin = 1.e-20 ! min value for water vapor (treated as zero)
    real, parameter :: qcmin = 1.e-12 ! min value for cloud condensates

    real, parameter :: vr_min = 1.e-3 ! min fall speed for rain
    real, parameter :: vf_min = 1.e-5 ! min fall speed for cloud ice, snow, graupel

    real, parameter :: dz_min = 1.e-2 ! used for correcting flipped height

    real, parameter :: sfcrho = 1.2 ! surface air density

    real, parameter :: rnzr = 8.0e6 ! lin et al. 1983
    real, parameter :: rnzs = 3.0e6 ! lin et al. 1983
    real, parameter :: rnzg = 4.0e6 ! rutledge and hobbs 1984
    ! lmh, 20170929
    real, parameter :: rnzh = 4.0e4 ! lin et al. 1983

    real, parameter :: rhow = 1.0e3 ! density of cloud water
    real, parameter :: rhor = 1.0e3 ! lin et al. 1983
    real, parameter :: rhos = 0.1e3 ! lin et al. 1983
    real, parameter :: rhog = 0.4e3 ! rutledge and hobbs 1984
    ! lmh, 20170929
    real, parameter :: rhoh = 0.917e3 ! lin et al. 1983

    real, parameter :: rgrav = 1. / grav

    real :: cracs, csacr, cgacr, cgacs, csacw, craci, csaci, cgacw, cgaci, cracw ! constants for accretions
    real :: acco (3, 4) ! constants for accretions
    ! constants for sublimation / deposition, freezing / melting, condensation / evaporation
    real :: cssub (5), cgsub (5), crevp (5), cgfr (2), csmlt (5), cgmlt (5)

    real :: es0, ces0
    real :: pie, fac_rc
    real :: c_air, c_vap

    real :: lat2, lcp, icp, tcp ! used in bigg mechanism and wet bulk

    real :: d0_vap ! the same as dc_vap, except that cp_vap can be cp_vap or cv_vap
    real (kind = r_grid) :: lv00, li00, li20
    real (kind = r_grid) :: d1_vap, d1_ice, c1_vap, c1_liq, c1_ice
    real (kind = r_grid), parameter :: one_r8 = 1.

    real, allocatable :: table (:), table2 (:), table3 (:), tablew (:)
    real, allocatable :: des (:), des2 (:), des3 (:), desw (:)

    logical :: tables_are_initialized = .false.

    real, parameter :: dt_fr = 8. ! homogeneous freezing of all cloud water at t_wfr - dt_fr
    ! minimum temperature water can exist (moore & molinero nov. 2011, nature)
    ! dt_fr can be considered as the error bar

    real, parameter :: p0_min = 100. ! minimum pressure (pascal) for mp to operate
    real :: p_min

    ! -----------------------------------------------------------------------
    ! namelist parameters
    ! -----------------------------------------------------------------------

    integer :: ntimes = 1 ! cloud microphysics sub cycles

    integer :: icloud_f = 0 ! cloud scheme
    integer :: irain_f = 0 ! cloud water to rain auto conversion scheme

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

    real :: cld_fac = 1.0 ! multiplication factor for cloud fraction
    real :: cld_min = 0.05 ! minimum cloud fraction
    real :: tice = 273.16 ! set tice = 165. to trun off ice - phase phys (kessler emulator)
    real :: tice_mlt = 273.16 ! set ice melting temperature to 268.0 based on observation (kay et al., 2016, jc)

    real :: t_min = 178. ! min temp to freeze - dry all water vapor
    real :: t_sub = 184. ! min temp for sublimation of cloud ice
    real :: mp_time = 150. ! maximum micro - physics time step (sec)

    real :: rh_inc = 0.25 ! rh increment for complete evaporation of cloud water and cloud ice
    real :: rh_inr = 0.25 ! rh increment for minimum evaporation of rain
    real :: rh_ins = 0.25 ! rh increment for sublimation of snow

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
    real :: tau_revp = 0. ! rain evaporation

    real :: dw_land = 0.20 ! base value for subgrid deviation / variability over land
    real :: dw_ocean = 0.10 ! base value for ocean

    real :: ccn_o = 90. ! ccn over ocean (cm^ - 3)
    real :: ccn_l = 270. ! ccn over land (cm^ - 3)

    real :: rthresh = 10.0e-6 ! critical cloud drop radius (micron)

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

    real :: ql_gen = 1.0e-3 ! max cloud water generation during remapping step if do_sat_adj = .t.
    real :: qi_gen = 1.82e-6 ! max cloud ice generation during remapping step

    ! cloud condensate upper bounds: "safety valves" for ql & qi

    real :: ql0_max = 2.0e-3 ! max cloud water value (auto converted to rain)
    real :: qi0_max = 1.0e-4 ! max cloud ice value (by other sources)

    real :: qi0_crt = 1.0e-4 ! cloud ice to snow autoconversion threshold (was 1.e-4)
    ! qi0_crt if negative, its magnitude is used as the mixing ration threshold; otherwise, used as density
    real :: qr0_crt = 1.0e-4 ! rain to snow or graupel / hail threshold
    ! lin et al. (1983) used * mixing ratio * = 1.e-4 (hail)
    real :: qs0_crt = 1.0e-3 ! snow to graupel density threshold (0.6e-3 in purdue lin scheme)

    real :: c_paut = 0.55 ! autoconversion cloud water to rain (use 0.5 to reduce autoconversion)
    real :: c_psaci = 0.02 ! accretion: cloud ice to snow (was 0.1 in zetac)
    real :: c_piacr = 5.0 ! accretion: rain to ice: (not used)
    real :: c_cracw = 0.9 ! rain accretion efficiency
    real :: c_pgacs = 2.0e-3 ! snow to graupel "accretion" eff. (was 0.1 in zetac)

    ! decreasing clin to reduce csacw (so as to reduce cloud water --- > snow)

    real :: alin = 842.0 ! "a" in lin et al. (1983)
    real :: clin = 4.8 ! "c" in lin et al. (1983), 4.8 -- > 6. (to ehance ql -- > qs)

    logical :: const_vi = .false. ! if .t. the constants are specified by v * _fac
    logical :: const_vs = .false. ! if .t. the constants are specified by v * _fac
    logical :: const_vg = .false. ! if .t. the constants are specified by v * _fac
    logical :: const_vr = .false. ! if .t. the constants are specified by v * _fac

    real :: vi_fac = 1. ! ifs: if const_vi: 1 / 3
    real :: vs_fac = 1. ! ifs: if const_vs: 1.
    real :: vg_fac = 1. ! ifs: if const_vg: 2.
    real :: vr_fac = 1. ! ifs: if const_vr: 4.

    real :: vi_max = 0.5 ! max fall speed for ice
    real :: vs_max = 5.0 ! max fall speed for snow
    real :: vg_max = 8.0 ! max fall speed for graupel
    real :: vr_max = 12. ! max fall speed for rain

    real :: xr_a = 0.25 ! p value in xu and randall, 1996
    real :: xr_b = 100. ! alpha_0 value in xu and randall, 1996
    real :: xr_c = 0.49 ! gamma value in xu and randall, 1996

    real :: te_err = 1.e-14 ! 64bit: 1.e-14, 32bit: 1.e-7

    logical :: do_sat_adj = .false. ! has fast saturation adjustments
    logical :: z_slope_liq = .true. ! use linear mono slope for autocconversions
    logical :: z_slope_ice = .false. ! use linear mono slope for autocconversions
    logical :: use_ccn = .false. ! must be true when prog_ccn is false
    logical :: use_ppm = .false. ! use ppm fall scheme
    logical :: use_ppm_ice = .false. ! use ppm fall scheme for cloud ice
    logical :: mono_prof = .true. ! perform terminal fall with mono ppm scheme
    logical :: do_hail = .false. ! use hail parameters instead of graupel
    logical :: hd_icefall = .false. ! use heymsfield and donner, 1990's fall speed of cloud ice
    logical :: use_xr_cloud = .false. ! use xu and randall, 1996's cloud diagnosis
    logical :: use_park_cloud = .false. ! park et al. 2016
    logical :: use_gi_cloud = .false. ! gultepe and isaac (2007, grl)
    logical :: use_rhc_cevap = .false. ! cap of rh for cloud water evaporation
    logical :: use_rhc_revap = .false. ! cap of rh for rain evaporation
    logical :: consv_checker = .false. ! turn on energy and water conservation checker
    logical :: do_warm_rain_mp = .false. ! do warm rain cloud microphysics only
    ! turn off to save time, turn on only in c48 64bit

    real :: g2, log_10

    real :: rh_thres = 0.75
    real :: rhc_cevap = 0.85 ! cloud water
    real :: rhc_revap = 0.85 ! cloud water

    real :: f_dq_p = 1.0
    real :: f_dq_m = 1.0
    logical :: do_cld_adj = .false.

    integer :: inflag = 1 ! ice nucleation scheme
    ! 1: hong et al., 2004
    ! 2: meyers et al., 1992
    ! 3: meyers et al., 1992
    ! 4: cooper, 1986
    ! 5: flecther, 1962

    ! -----------------------------------------------------------------------
    ! namelist
    ! -----------------------------------------------------------------------

    namelist / gfdl_mp_nml / &
        t_min, t_sub, tau_r2g, tau_smlt, tau_g2r, dw_land, dw_ocean, &
        vi_fac, vr_fac, vs_fac, vg_fac, ql_mlt, do_qa, fix_negative, vi_max, &
        vs_max, vg_max, vr_max, qs_mlt, qs0_crt, qi_gen, ql0_max, qi0_max, &
        qi0_crt, do_sat_adj, rh_inc, rh_ins, rh_inr, const_vi, &
        const_vs, const_vg, const_vr, use_ccn, rthresh, ccn_l, ccn_o, qc_crt, &
        tau_g2v, tau_v2g, sat_adj0, tau_imlt, tau_v2l, tau_l2v, &
        tau_i2s, tau_l2r, qi_lim, ql_gen, c_paut, c_psaci, c_pgacs, &
        z_slope_liq, z_slope_ice, prog_ccn, c_cracw, alin, clin, tice, &
        rad_snow, rad_graupel, rad_rain, cld_fac, cld_min, use_ppm, use_ppm_ice, mono_prof, &
        do_sedi_heat, sedi_transport, do_sedi_w, icloud_f, irain_f, &
        ntimes, disp_heat, do_hail, use_xr_cloud, xr_a, xr_b, xr_c, tau_revp, tice_mlt, hd_icefall, &
        do_cond_timescale, mp_time, consv_checker, te_err, use_park_cloud, &
        use_gi_cloud, use_rhc_cevap, use_rhc_revap, inflag, do_warm_rain_mp, &
        rh_thres, f_dq_p, f_dq_m, do_cld_adj

    public &
        t_min, t_sub, tau_r2g, tau_smlt, tau_g2r, dw_land, dw_ocean, &
        vi_fac, vr_fac, vs_fac, vg_fac, ql_mlt, do_qa, fix_negative, vi_max, &
        vs_max, vg_max, vr_max, qs_mlt, qs0_crt, qi_gen, ql0_max, qi0_max, &
        qi0_crt, do_sat_adj, rh_inc, rh_ins, rh_inr, const_vi, &
        const_vs, const_vg, const_vr, use_ccn, rthresh, ccn_l, ccn_o, qc_crt, &
        tau_g2v, tau_v2g, sat_adj0, tau_imlt, tau_v2l, tau_l2v, &
        tau_i2s, tau_l2r, qi_lim, ql_gen, c_paut, c_psaci, c_pgacs, &
        z_slope_liq, z_slope_ice, prog_ccn, c_cracw, alin, clin, tice, &
        rad_snow, rad_graupel, rad_rain, cld_fac, cld_min, use_ppm, use_ppm_ice, mono_prof, &
        do_sedi_heat, sedi_transport, do_sedi_w, icloud_f, irain_f, &
        ntimes, disp_heat, do_hail, use_xr_cloud, xr_a, xr_b, xr_c, tau_revp, tice_mlt, hd_icefall, &
        do_cond_timescale, mp_time, consv_checker, te_err, use_park_cloud, &
        use_gi_cloud, use_rhc_cevap, use_rhc_revap, inflag, do_warm_rain_mp, &
        rh_thres, f_dq_p, f_dq_m, do_cld_adj

contains

! -----------------------------------------------------------------------
! the driver of the gfdl cloud microphysics
! -----------------------------------------------------------------------

subroutine gfdl_cld_mp_driver (qv, ql, qr, qi, qs, qg, qa, qnl, qni, &
        pt, w, ua, va, dz, delp, gsize, dts, hs, rain, snow, ice, &
        graupel, hydrostatic, is, ie, ks, ke, q_con, cappa, consv_te, &
        te, condensation, deposition, evaporation, sublimation, last_step, do_inline_mp)

    implicit none

    logical, intent (in) :: hydrostatic
    logical, intent (in) :: last_step
    logical, intent (in) :: consv_te
    logical, intent (in) :: do_inline_mp

    integer, intent (in) :: is, ie ! physics window
    integer, intent (in) :: ks, ke ! vertical dimension

    real, intent (in) :: dts ! physics time step

    real, intent (in), dimension (is:ie) :: hs, gsize

    real, intent (in), dimension (is:ie, ks:ke) :: dz
    real, intent (in), dimension (is:ie, ks:ke) :: qnl, qni

    real, intent (inout), dimension (is:ie, ks:ke) :: delp
    real, intent (inout), dimension (is:ie, ks:ke) :: qv, ql, qr, qi, qs, qg, qa
    real, intent (inout), dimension (is:ie, ks:ke) :: pt, ua, va, w
    real, intent (inout), dimension (is:, ks:) :: q_con, cappa
    real, intent (inout), dimension (is:ie) :: rain, snow, ice, graupel
    real, intent (inout), dimension (is:ie) :: condensation, deposition
    real, intent (inout), dimension (is:ie) :: evaporation, sublimation

    real, intent (inout), dimension (is:ie, ks:ke) :: te
    ! logical :: used
    real, dimension (is:ie) :: w_var
    real, dimension (is:ie, ks:ke) :: vt_r, vt_s, vt_g, vt_i
    real, dimension (is:ie, ks:ke) :: m2_rain, m2_sol

    if (last_step) then
        p_min = p0_min ! final clean - up
    else
        p_min = 30.e2 ! time saving trick
    endif

    ! -----------------------------------------------------------------------
    ! define heat capacity of dry air and water vapor based on hydrostatical property
    ! -----------------------------------------------------------------------

    if (hydrostatic) then
        c_air = cp_air
        c_vap = cp_vap
        do_sedi_w = .false.
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
        qa, qnl, qni, dz, is, ie, ks, ke, dts, &
        rain, snow, graupel, ice, m2_rain, m2_sol, gsize, hs, &
        w_var, vt_r, vt_s, vt_g, vt_i, q_con, cappa, consv_te, te, &
        condensation, deposition, evaporation, sublimation, last_step, do_inline_mp)

end subroutine gfdl_cld_mp_driver

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
        qg, qa, qnl, qni, dz, is, ie, ks, ke, dt_in, &
        rain, snow, graupel, ice, m2_rain, m2_sol, gsize, hs, &
        w_var, vt_r, vt_s, vt_g, vt_i, q_con, cappa, consv_te, te, &
        condensation, deposition, evaporation, sublimation, last_step, do_inline_mp)

    implicit none

    logical, intent (in) :: hydrostatic
    logical, intent (in) :: last_step
    logical, intent (in) :: consv_te
    logical, intent (in) :: do_inline_mp
    integer, intent (in) :: is, ie, ks, ke
    real, intent (in) :: dt_in
    real, intent (in), dimension (is:ie) :: gsize
    real, intent (in), dimension (is:ie) :: hs
    real, intent (in), dimension (is:ie, ks:ke) :: dz
    real, intent (in), dimension (is:ie, ks:ke) :: qnl, qni

    real, intent (inout), dimension (is:ie, ks:ke) :: delp
    real, intent (inout), dimension (is:ie, ks:ke) :: qv, ql, qr, qi, qs, qg, qa
    real, intent (inout), dimension (is:ie, ks:ke) :: pt, ua, va, w
    real, intent (inout), dimension (is:, ks:) :: q_con, cappa
    real, intent (inout), dimension (is:ie) :: rain, snow, ice, graupel
    real, intent (inout), dimension (is:ie) :: condensation, deposition
    real, intent (inout), dimension (is:ie) :: evaporation, sublimation

    real, intent (out), dimension (is:ie) :: w_var
    real, intent (out), dimension (is:ie, ks:ke) :: vt_r, vt_s, vt_g, vt_i
    real, intent (out), dimension (is:ie, ks:ke) :: m2_rain, m2_sol
    real, intent (out), dimension (is:ie, ks:ke) :: te
    ! local:
    real, dimension (ks:ke) :: q_liq, q_sol
    real, dimension (ks:ke) :: qvz, qlz, qrz, qiz, qsz, qgz, qaz
    real, dimension (ks:ke) :: vtiz, vtsz, vtgz, vtrz
    real, dimension (ks:ke) :: dp1, dz1
    real, dimension (ks:ke) :: den, p1, denfac
    real, dimension (ks:ke) :: ccn, cin, c_praut, m1_rain, m1_sol, m1
    real, dimension (ks:ke) :: u0, v0, u1, v1, w1

    real (kind = r_grid), dimension (is:ie, ks:ke) :: te_beg, te_end, tw_beg, tw_end
    real (kind = r_grid), dimension (is:ie, ks:ke) :: te_beg_0, te_end_0, tw_beg_0, tw_end_0
    real (kind = r_grid), dimension (is:ie) :: te_b_beg, te_b_end, tw_b_beg, tw_b_end, dte, te_loss
    real (kind = r_grid), dimension (is:ie) :: te_b_beg_0, te_b_end_0, tw_b_beg_0, tw_b_end_0
    real (kind = r_grid), dimension (ks:ke) :: te1, te2

    real :: cpaut, rh_adj, rh_rain
    real :: r1, s1, i1, g1, rdt, ccn0
    real :: dt_rain
    real :: s_leng, t_land, t_ocean, h_var, tmp
    real (kind = r_grid), dimension (ks:ke) :: dp0, tz, cvm
    real (kind = r_grid) :: con_r8, c8
    real :: convt
    real :: dts, q_cond
    real :: cond, dep, reevap, sub

    integer :: i, k, n

    ntimes = max (ntimes, int (dt_in / min (dt_in, mp_time)))
    dts = dt_in / real (ntimes)

    dt_rain = dts * 0.5
    rdt = one_r8 / dts

    dte = 0.0

    ! convert to mm / day
    convt = 86400. * rdt * rgrav
    cond = 0.0

    ! -----------------------------------------------------------------------
    ! use local variables
    ! -----------------------------------------------------------------------

    do i = is, ie

        do k = ks, ke
            if (do_inline_mp) then
#ifdef MOIST_CAPPA
                tz (k) = pt (i, k) / ((1. + zvir * qv (i, k)) * (1. - (ql (i, k) + qr (i, k) + qi (i, k) + qs (i, k) + qg (i, k))))
#else
                tz (k) = pt (i, k) / (1. + zvir * qv (i, k))
#endif
            else
                tz (k) = pt (i, k)
            endif
        enddo

        ! -----------------------------------------------------------------------
        ! total energy checker
        ! -----------------------------------------------------------------------

        if (consv_checker) then
            do k = ks, ke
                q_liq (k) = ql (i, k) + qr (i, k)
                q_sol (k) = qi (i, k) + qs (i, k) + qg (i, k)
                cvm (k) = c_air * (1.0 - qv (i, k) - q_liq (k) - q_sol (k)) + &
                    qv (i, k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
                te_beg_0 (i, k) = cvm (k) * tz (k) + lv00 * c_air * qv (i, k) - li00 * c_air * q_sol (k)
                if (hydrostatic) then
                    te_beg_0 (i, k) = te_beg_0 (i, k) + 0.5 * (ua (i, k) ** 2 + va (i, k) ** 2)
                else
                    te_beg_0 (i, k) = te_beg_0 (i, k) + 0.5 * (ua (i, k) ** 2 + va (i, k) ** 2 + w (i, k) ** 2)
                endif
                te_beg_0 (i, k) = rgrav * te_beg_0 (i, k) * delp (i, k) * gsize (i) ** 2.0
                tw_beg_0 (i, k) = rgrav * (qv (i, k) + q_liq (k) + q_sol (k)) * delp (i, k) * gsize (i) ** 2.0
            enddo
            te_b_beg_0 (i) = (dte (i) - li00 * c_air * (ice (i) + snow (i) + graupel (i)) * dt_in / 86400) * gsize (i) ** 2.0
            tw_b_beg_0 (i) = (rain (i) + ice (i) + snow (i) + graupel (i)) * dt_in / 86400 * gsize (i) ** 2.0
        endif

        do k = ks, ke
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
            q_cond = q_liq (k) + q_sol (k)
            qaz (k) = 0.
            dz1 (k) = dz (i, k)
            con_r8 = one_r8 - (qvz (k) + q_cond)
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
            if (.not. hydrostatic) then
                w1 (k) = w (i, k)
            endif
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
                    q_cond = q_liq (k) + q_sol (k)
                    cvm (k) = (one_r8 - (qv (i, k) + q_cond)) * c_air + &
                        qv (i, k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
                    te (i, k) = - cvm (k) * tz (k) * delp (i, k)
#else
                    te (i, k) = - c_air * tz (k) * delp (i, k)
#endif
                enddo
            endif
        endif

        ! -----------------------------------------------------------------------
        ! total energy checker
        ! -----------------------------------------------------------------------

        if (consv_checker) then
            do k = ks, ke
                q_liq (k) = qlz (k) + qrz (k)
                q_sol (k) = qiz (k) + qsz (k) + qgz (k)
                cvm (k) = c_air + qvz (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
                te_beg (i, k) = cvm (k) * tz (k) + lv00 * c_air * qvz (k) - li00 * c_air * q_sol (k)
                if (hydrostatic) then
                    te_beg (i, k) = te_beg (i, k) + 0.5 * (u1 (k) ** 2 + v1 (k) ** 2)
                else
                    te_beg (i, k) = te_beg (i, k) + 0.5 * (u1 (k) ** 2 + v1 (k) ** 2 + w1 (k) ** 2)
                endif
                te_beg (i, k) = rgrav * te_beg (i, k) * dp1 (k) * gsize (i) ** 2.0
                tw_beg (i, k) = rgrav * (qvz (k) + q_liq (k) + q_sol (k)) * dp1 (k) * gsize (i) ** 2.0
            enddo
            te_b_beg (i) = (dte (i) - li00 * c_air * (ice (i) + snow (i) + graupel (i)) * dt_in / 86400) * gsize (i) ** 2.0
            tw_b_beg (i) = (rain (i) + ice (i) + snow (i) + graupel (i)) * dt_in / 86400 * gsize (i) ** 2.0
        endif

        ! -----------------------------------------------------------------------
        ! calculate cloud condensation nuclei (ccn)
        ! the following is based on klein eq. 15
        ! -----------------------------------------------------------------------

        cpaut = c_paut * 0.104 * grav / 1.717e-5

        if (prog_ccn) then
            do k = ks, ke
                ! convert # / cm^3 to # / m^3
                ccn (k) = max (10.0, qnl (i, k)) * 1.e6
                cin (k) = max (10.0, qni (i, k)) * 1.e6
                ccn (k) = ccn (k) / den (k)
                c_praut (k) = cpaut * (ccn (k) * rhor) ** (- 1. / 3.)
            enddo
        else
            ! convert # / cm^3 to # / m^3
            ccn0 = (ccn_l * min (1., abs (hs (i)) / (10. * grav)) + &
                ccn_o * (1. - min (1., abs (hs (i)) / (10. * grav)))) * 1.e6
            do k = ks, ke
                ccn (k) = ccn0 / den (k)
                c_praut (k) = cpaut * (ccn (k) * rhor) ** (- 1. / 3.)
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
        rh_rain = max (0.35, rh_adj - rh_inr) ! rh_inr = 0.25

        ! -----------------------------------------------------------------------
        ! fix all negative water species
        ! -----------------------------------------------------------------------

        if (fix_negative) &
            call neg_adj (ks, ke, tz, dp1, qvz, qlz, qrz, qiz, qsz, qgz, cond)

        condensation (i) = condensation (i) + cond * convt * ntimes

        m2_rain (i, :) = 0.
        m2_sol (i, :) = 0.

        do n = 1, ntimes

            ! -----------------------------------------------------------------------
            ! time - split warm rain processes: 1st pass
            ! -----------------------------------------------------------------------

            call warm_rain (dt_rain, ks, ke, dp1, dz1, tz, qvz, qlz, qrz, qiz, qsz, &
                qgz, den, denfac, ccn, c_praut, rh_rain, vtrz, r1, m1_rain, w1, h_var, reevap, dte (i))

            evaporation (i) = evaporation (i) + reevap * convt
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
                dz1, dp1, den, vtgz, vtsz, vtiz, r1, g1, s1, i1, m1_sol, w1, dte (i))

            rain (i) = rain (i) + r1 * convt ! from melted snow & ice that reached the ground
            snow (i) = snow (i) + s1 * convt
            graupel (i) = graupel (i) + g1 * convt
            ice (i) = ice (i) + i1 * convt

            ! -----------------------------------------------------------------------
            ! energy loss during sedimentation heating
            ! -----------------------------------------------------------------------

            if (consv_checker) then
                do k = ks, ke
                    te1 (k) = one_r8 + qvz (k) * c1_vap + (qlz (k) + qrz (k)) * c1_liq + (qiz (k) + qsz (k) + qgz (k)) * c1_ice
                    te1 (k) = rgrav * te1 (k) * c_air * tz (k) * dp1 (k)
                enddo
            endif

            ! -----------------------------------------------------------------------
            ! heat transportation during sedimentation
            ! -----------------------------------------------------------------------

            if (do_sedi_heat) then
                call sedi_heat (ks, ke, dp1, m1_sol, dz1, tz, qvz, qlz, qrz, qiz, &
                    qsz, qgz, c_ice)
            endif

            ! -----------------------------------------------------------------------
            ! energy loss during sedimentation heating
            ! -----------------------------------------------------------------------

            if (consv_checker) then
                do k = ks, ke
                    te2 (k) = one_r8 + qvz (k) * c1_vap + (qlz (k) + qrz (k)) * c1_liq + (qiz (k) + qsz (k) + qgz (k)) * c1_ice
                    te2 (k) = rgrav * te2 (k) * c_air * tz (k) * dp1 (k)
                enddo
                dte (i) = dte (i) + sum (te1) - sum (te2)
            endif

            ! -----------------------------------------------------------------------
            ! time - split warm rain processes: 2nd pass
            ! -----------------------------------------------------------------------

            call warm_rain (dt_rain, ks, ke, dp1, dz1, tz, qvz, qlz, qrz, qiz, qsz, &
                qgz, den, denfac, ccn, c_praut, rh_rain, vtrz, r1, m1_rain, w1, h_var, reevap, dte (i))

            evaporation (i) = evaporation (i) + reevap * convt
            rain (i) = rain (i) + r1 * convt

            do k = ks, ke
                m2_rain (i, k) = m2_rain (i, k) + m1_rain (k)
                m2_sol (i, k) = m2_sol (i, k) + m1_sol (k)
                m1 (k) = m1 (k) + m1_rain (k) + m1_sol (k)
            enddo

            ! -----------------------------------------------------------------------
            ! ice - phase microphysics
            ! -----------------------------------------------------------------------

            call icloud (ks, ke, tz, p1, qvz, qlz, qrz, qiz, qsz, qgz, dp1, den, ccn, &
                cin, denfac, vtsz, vtgz, vtrz, qaz, rh_adj, rh_rain, dts, h_var, gsize (i), &
                cond, dep, reevap, sub, last_step)

            condensation (i) = condensation (i) + cond * convt
            deposition (i) = deposition (i) + dep * convt
            evaporation (i) = evaporation (i) + reevap * convt
            sublimation (i) = sublimation (i) + sub * convt

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
                    tz (k) = tz (k) + 0.5 * (w (i, k) ** 2 - w1 (k) ** 2) / c8
#else
                    tz (k) = tz (k) + 0.5 * (w (i, k) ** 2 - w1 (k) ** 2) / c_air
#endif
                enddo
            endif
            !#endif
            do k = ks, ke
                w (i, k) = w1 (k)
            enddo
        endif

        ! -----------------------------------------------------------------------
        ! total energy checker
        ! -----------------------------------------------------------------------

        if (consv_checker) then
            do k = ks, ke
                q_liq (k) = qlz (k) + qrz (k)
                q_sol (k) = qiz (k) + qsz (k) + qgz (k)
                cvm (k) = c_air + qvz (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
                te_end (i, k) = cvm (k) * tz (k) + lv00 * c_air * qvz (k) - li00 * c_air * q_sol (k)
                if (hydrostatic) then
                    te_end (i, k) = te_end (i, k) + 0.5 * (u1 (k) ** 2 + v1 (k) ** 2)
                else
                    te_end (i, k) = te_end (i, k) + 0.5 * (u1 (k) ** 2 + v1 (k) ** 2 + w1 (k) ** 2)
                endif
                te_end (i, k) = rgrav * te_end (i, k) * dp1 (k) * gsize (i) ** 2.0
                tw_end (i, k) = rgrav * (qvz (k) + q_liq (k) + q_sol (k)) * dp1 (k) * gsize (i) ** 2.0
            enddo
            te_b_end (i) = (dte (i) - li00 * c_air * (ice (i) + snow (i) + graupel (i)) * dt_in / 86400) * gsize (i) ** 2.0
            tw_b_end (i) = (rain (i) + ice (i) + snow (i) + graupel (i)) * dt_in / 86400 * gsize (i) ** 2.0
            ! total energy loss due to sedimentation and its heating
            te_loss (i) = dte (i) * gsize (i) ** 2.0
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
            q_cond = q_liq (k) + q_sol (k)
            cvm (k) = (one_r8 - (qvz (k) + q_cond)) * c_air + qvz (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
#ifdef MOIST_CAPPA
            q_con (i, k) = q_cond
            tmp = rdgas * (1. + zvir * qvz (k))
            cappa (i, k) = tmp / (tmp + cvm (k))
#endif
            if (do_inline_mp) then
#ifdef MOIST_CAPPA
                pt (i, k) = tz (k) * (1. + zvir * qvz (k)) * (1. - q_cond)
#else
                pt (i, k) = tz (k) * (1. + zvir * qvz (k))
#endif
            else
                pt (i, k) = pt (i, k) + (tz (k) - pt (i, k)) * cvm (k) / cp_air
            endif
        enddo

        ! -----------------------------------------------------------------------
        ! total energy checker
        ! -----------------------------------------------------------------------

        if (consv_checker) then
            do k = ks, ke
                q_liq (k) = ql (i, k) + qr (i, k)
                q_sol (k) = qi (i, k) + qs (i, k) + qg (i, k)
                cvm (k) = c_air * (1.0 - qv (i, k) - q_liq (k) - q_sol (k)) + &
                    qv (i, k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
                te_end_0 (i, k) = cvm (k) * tz (k) + lv00 * c_air * qv (i, k) - li00 * c_air * q_sol (k)
                te_end_0 (i, k) = te_end_0 (i, k) + 0.5 * (ua (i, k) ** 2 + va (i, k) ** 2 + w (i, k) ** 2)
                te_end_0 (i, k) = rgrav * te_end_0 (i, k) * delp (i, k) * gsize (i) ** 2.0
                tw_end_0 (i, k) = rgrav * (qv (i, k) + q_liq (k) + q_sol (k)) * delp (i, k) * gsize (i) ** 2.0
            enddo
            te_b_end_0 (i) = (dte (i) - li00 * c_air * (ice (i) + snow (i) + graupel (i)) * dt_in / 86400) * gsize (i) ** 2.0
            tw_b_end_0 (i) = (rain (i) + ice (i) + snow (i) + graupel (i)) * dt_in / 86400 * gsize (i) ** 2.0
        endif

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

    ! -----------------------------------------------------------------------
    ! total energy checker
    ! -----------------------------------------------------------------------

    if (consv_checker) then
        if (abs (sum (te_end) + sum (te_b_end) - sum (te_beg) - sum (te_b_beg)) / (sum (te_beg) + sum (te_b_beg)) .gt. te_err) then
            print *, "gfdl_cld_mp te: ", sum (te_beg) / sum (gsize ** 2) + sum (te_b_beg) / sum (gsize ** 2), &
                sum (te_end) / sum (gsize ** 2) + sum (te_b_end) / sum (gsize ** 2), &
                 (sum (te_end) + sum (te_b_end) - sum (te_beg) - sum (te_b_beg)) / (sum (te_beg) + sum (te_b_beg))
        endif
        if (abs (sum (tw_end) + sum (tw_b_end) - sum (tw_beg) - sum (tw_b_beg)) / (sum (tw_beg) + sum (tw_b_beg)) .gt. te_err) then
            print *, "gfdl_cld_mp tw: ", sum (tw_beg) / sum (gsize ** 2) + sum (tw_b_beg) / sum (gsize ** 2), &
                sum (tw_end) / sum (gsize ** 2) + sum (tw_b_end) / sum (gsize ** 2), &
                 (sum (tw_end) + sum (tw_b_end) - sum (tw_beg) - sum (tw_b_beg)) / (sum (tw_beg) + sum (tw_b_beg))
        endif
        ! print *, "gfdl_cld_mp te loss (%) : ", sum (te_loss) / (sum (te_beg) + sum (te_b_beg)) * 100.0
    endif

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
        den, denfac, ccn, c_praut, rh_rain, vtr, r1, m1_rain, w1, h_var, reevap, dte)

    implicit none

    integer, intent (in) :: ks, ke
    real, intent (in) :: dt ! time step (s)
    real, intent (in) :: rh_rain, h_var
    real, intent (in), dimension (ks:ke) :: dp, dz, den
    real, intent (in), dimension (ks:ke) :: denfac, ccn, c_praut

    real (kind = r_grid), intent (inout), dimension (ks:ke) :: tz
    real, intent (inout), dimension (ks:ke) :: vtr, qv, ql, qr, qi, qs, qg, m1_rain, w1
    real (kind = r_grid), intent (inout) :: dte
    real, intent (out) :: r1
    real, intent (out) :: reevap
    real, parameter :: so3 = 7. / 3.
    ! fall velocity constants:
    real, parameter :: vconr = 2503.23638966667
    real, parameter :: normr = 25132741228.7183
    real, parameter :: thr = 1.e-8

    real, dimension (ks:ke) :: dl, dm
    real (kind = r_grid), dimension (ks:ke) :: te1, te2
    real, dimension (ks:ke + 1) :: ze, zt
    real :: sink, dq, qc
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

    reevap = 0

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

        call revap_racc (ks, ke, dt5, tz, qv, ql, qr, qi, qs, qg, den, denfac, rh_rain, h_var, dp, reevap)

        if (do_sedi_w) then
            do k = ks, ke
                dm (k) = dp (k) * (1. + qv (k) + ql (k) + qr (k) + qi (k) + qs (k) + qg (k))
            enddo
        endif

        ! -----------------------------------------------------------------------
        ! energy loss during sedimentation
        ! -----------------------------------------------------------------------

        if (consv_checker) then
            do k = ks, ke
                te1 (k) = one_r8 + qv (k) * c1_vap + (ql (k) + qr (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
                te1 (k) = rgrav * te1 (k) * c_air * tz (k) * dp (k)
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
        ! energy loss during sedimentation
        ! -----------------------------------------------------------------------

        if (consv_checker) then
            do k = ks, ke
                te2 (k) = one_r8 + qv (k) * c1_vap + (ql (k) + qr (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
                te2 (k) = rgrav * te2 (k) * c_air * tz (k) * dp (k)
            enddo
            dte = dte + sum (te1) - sum (te2)
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
        ! energy loss during sedimentation heating
        ! -----------------------------------------------------------------------

        if (consv_checker) then
            do k = ks, ke
                te1 (k) = one_r8 + qv (k) * c1_vap + (ql (k) + qr (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
                te1 (k) = rgrav * te1 (k) * c_air * tz (k) * dp (k)
            enddo
        endif

        ! -----------------------------------------------------------------------
        ! heat exchanges during sedimentation
        ! -----------------------------------------------------------------------

        if (do_sedi_heat) then
            call sedi_heat (ks, ke, dp, m1_rain, dz, tz, qv, ql, qr, qi, qs, qg, c_liq)
        endif

        ! -----------------------------------------------------------------------
        ! energy loss during sedimentation heating
        ! -----------------------------------------------------------------------

        if (consv_checker) then
            do k = ks, ke
                te2 (k) = one_r8 + qv (k) * c1_vap + (ql (k) + qr (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
                te2 (k) = rgrav * te2 (k) * c_air * tz (k) * dp (k)
            enddo
            dte = dte + sum (te1) - sum (te2)
        endif


        ! -----------------------------------------------------------------------
        ! evaporation and accretion of rain for the remaing 1 / 2 time step
        ! -----------------------------------------------------------------------

        call revap_racc (ks, ke, dt5, tz, qv, ql, qr, qi, qs, qg, den, denfac, rh_rain, h_var, dp, reevap)

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
            qc = fac_rc * ccn (k)
            if (tz (k) > t_wfr) then
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
            qc = fac_rc * ccn (k)
            if (tz (k) > t_wfr + dt_fr) then
                dl (k) = min (max (1.e-6, dl (k)), 0.5 * ql (k))
                ! --------------------------------------------------------------------
                ! as in klein's gfdl am2 stratiform scheme (with subgrid variations)
                ! --------------------------------------------------------------------
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

subroutine revap_racc (ks, ke, dt, tz, qv, ql, qr, qi, qs, qg, den, denfac, rh_rain, h_var, dp, reevap)

    implicit none

    integer, intent (in) :: ks, ke
    real, intent (in) :: dt ! time step (s)
    real, intent (in) :: rh_rain, h_var
    real, intent (in), dimension (ks:ke) :: den, denfac, dp
    real (kind = r_grid), intent (inout), dimension (ks:ke) :: tz
    real, intent (inout), dimension (ks:ke) :: qv, qr, ql, qi, qs, qg
    real, intent (out) :: reevap
    ! local:
    real (kind = r_grid), dimension (ks:ke) :: cvm
    real, dimension (ks:ke) :: q_liq, q_sol, lcpk
    real :: dqv, qsat, dqsdt, evap, t2, qden, q_plus, q_minus, sink
    real :: qpz, dq, dqh, tin
    real :: fac_revp, rh_tem

    integer :: k

    if (tau_revp .gt. 1.e-6) then
        fac_revp = 1. - exp (- dt / tau_revp)
    else
        fac_revp = 1.
    endif

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

            rh_tem = qpz / iqs1 (tin, den (k))

            if (dqv > qvmin .and. qsat > q_minus) then
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
                if (use_rhc_revap) then
                    evap = 0.0
                    if (rh_tem < rhc_revap) then
                        evap = crevp (1) * t2 * dq * (crevp (2) * sqrt (qden) + crevp (3) * &
                            exp (0.725 * log (qden)) * sqrt (denfac (k))) / (crevp (4) * t2 + crevp (5) * qsat * den (k))
                        evap = min (qr (k), dt * fac_revp * evap, dqv / (1. + lcpk (k) * dqsdt))
                    endif
                else
                    evap = crevp (1) * t2 * dq * (crevp (2) * sqrt (qden) + crevp (3) * &
                        exp (0.725 * log (qden))) / (crevp (4) * t2 + crevp (5) * qsat * den (k))
                    evap = min (qr (k), dt * fac_revp * evap, dqv / (1. + lcpk (k) * dqsdt))
                endif
                reevap = reevap + evap * dp (k)

                ! -----------------------------------------------------------------------
                ! alternative minimum evap in dry environmental air
                ! sink = min (qr (k), dim (rh_rain * qsat, qv (k)) / (1. + lcpk (k) * dqsdt))
                ! evap = max (evap, sink)
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
            if (qr (k) > qrmin .and. ql (k) > 1.e-6 .and. qsat < q_minus) then
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

subroutine icloud (ks, ke, tzk, p1, qvk, qlk, qrk, qik, qsk, qgk, dp1, den, &
        ccn, cin, denfac, vts, vtg, vtr, qak, rh_adj, rh_rain, dts, h_var, &
        gsize, cond, dep, reevap, sub, last_step)

    implicit none

    logical, intent (in) :: last_step
    integer, intent (in) :: ks, ke
    real, intent (in), dimension (ks:ke) :: p1, dp1, den, denfac, vts, vtg, vtr, ccn
    real (kind = r_grid), intent (inout), dimension (ks:ke) :: tzk
    real, intent (inout), dimension (ks:ke) :: qvk, qlk, qrk, qik, qsk, qgk, qak
    real, intent (inout), dimension (ks:ke) :: cin
    real, intent (in) :: rh_adj, rh_rain, dts, h_var, gsize
    real, intent (out) :: cond, dep, reevap, sub
    ! local:
    real, dimension (ks:ke) :: icpk, di, qim
    real, dimension (ks:ke) :: q_liq, q_sol
    real (kind = r_grid), dimension (ks:ke) :: cvm, te8
    real (kind = r_grid) :: tz
    real :: rdts, fac_g2v, fac_v2g, fac_i2s, fac_imlt
    real :: qv, ql, qr, qi, qs, qg, melt
    real :: pracs, psacw, pgacw, psacr, pgacr, pgaci, praci, psaci
    real :: pgmlt, psmlt, pgfr, psaut
    real :: tc, dqs0, qden, qsm
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
    ! similar to lfo 1983: eq. 21 solved implicitly
    ! threshold from wsm6 scheme, hong et al 2004, eq (13) : qi0_crt ~0.8e-4
    ! -----------------------------------------------------------------------

    do k = ks, ke
        if (qi0_crt < 0.) then
            qim (k) = - qi0_crt
        else
            qim (k) = qi0_crt / den (k)
        endif
    enddo

    if (.not. do_warm_rain_mp) then

        ! -----------------------------------------------------------------------
        ! sources of cloud ice: pihom, cold rain, and the sat_adj
        ! (initiation plus deposition)
        ! sources of snow: cold rain, auto conversion + accretion (from cloud ice)
        ! sat_adj (deposition; requires pre - existing snow) ; initial snow comes from auto conversion
        ! -----------------------------------------------------------------------

        do k = ks, ke
            if (tzk (k) > tice_mlt .and. qik (k) > qcmin) then

                ! -----------------------------------------------------------------------
                ! pimlt: instant melting of cloud ice
                ! -----------------------------------------------------------------------

                melt = min (qik (k), fac_imlt * (tzk (k) - tice_mlt) / icpk (k))
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
                tmp = min (sink, dim (qim (k), qik (k)))
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

                if (qi > 3.e-7) then ! cloud ice sink terms

                    if (qs > 1.e-7) then
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
                    if (q_plus > (qim (k) + qrmin)) then
                        if (qim (k) > (qi - di (k))) then
                            dq = (0.25 * (q_plus - qim (k)) ** 2) / di (k)
                        else
                            dq = qi - qim (k)
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

                    if (qg > 1.e-6) then
                        ! -----------------------------------------------------------------------
                        ! factor = dts * cgaci / sqrt (den (k)) * exp (0.05 * tc + 0.875 * log (qg * den (k)))
                        ! simplified form: remove temp dependency & set the exponent "0.875" -- > 1
                        ! -----------------------------------------------------------------------
                        factor = dts * cgaci / sqrt (den (k)) * exp (0.875 * log (qg * den (k)))
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

                if (qr > 1.e-7 .and. tc < 0.) then

                    ! -----------------------------------------------------------------------
                    ! * sink * terms to qr: psacr + pgfr
                    ! source terms to qs: psacr
                    ! source terms to qg: pgfr
                    ! -----------------------------------------------------------------------

                    ! -----------------------------------------------------------------------
                    ! psacr accretion of rain by snow
                    ! -----------------------------------------------------------------------

                    if (qs > 1.e-7) then ! if snow exists
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

                if (qs > 1.e-7) then

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

                if (qg > 1.e-7 .and. tz < tice) then

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

    endif

    call subgrid_z_proc (ks, ke, p1, den, denfac, dts, rh_adj, tzk, qvk, qlk, &
        qrk, qik, qsk, qgk, qak, dp1, h_var, rh_rain, te8, ccn, cin, gsize, &
        cond, dep, reevap, sub, last_step)

end subroutine icloud

! =======================================================================
! temperature sentive high vertical resolution processes
! =======================================================================

subroutine subgrid_z_proc (ks, ke, p1, den, denfac, dts, rh_adj, tz, qv, ql, qr, &
        qi, qs, qg, qa, dp1, h_var, rh_rain, te8, ccn, cin, gsize, cond, dep, reevap, sub, last_step)

    implicit none

    integer, intent (in) :: ks, ke
    real, intent (in) :: dts, rh_adj, h_var, rh_rain, gsize
    real, intent (in), dimension (ks:ke) :: p1, den, denfac, ccn, dp1
    real (kind = r_grid), intent (in), dimension (ks:ke) :: te8
    real (kind = r_grid), intent (inout), dimension (ks:ke) :: tz
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, qa
    real, intent (inout), dimension (ks:ke) :: cin
    logical, intent (in) :: last_step
    real, intent (out) :: cond, dep, reevap, sub
    ! local:
    real, dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    real, dimension (ks:ke) :: q_liq, q_sol, q_cond
    real (kind = r_grid), dimension (ks:ke) :: cvm
    real :: pidep, qi_crt
    real :: sigma, gam
    ! -----------------------------------------------------------------------
    ! qstar over water may be accurate only down to - 80 deg c with ~10% uncertainty
    ! must not be too large to allow psc
    ! -----------------------------------------------------------------------
    real :: rh, rqi, tin, qsw, qsi, qpz, qstar, rh_tem
    real :: dqsdt, dwsdt, dq, dq0, factor, tmp, liq, ice
    real :: q_plus, q_minus, dt_evap, dt_pisub
    real :: evap, sink, tc, dtmp, qa10, qa100
    real :: pssub, pgsub, tsq, qden
    real :: fac_l2v, fac_v2l, fac_g2v, fac_v2g
    integer :: k

    if (do_sat_adj) then
        dt_evap = 0.5 * dts
    else
        dt_evap = dts
    endif

    ! -----------------------------------------------------------------------
    ! define conversion scalar / factor
    ! -----------------------------------------------------------------------

    fac_l2v = 1. - exp (- dt_evap / tau_l2v)
    fac_v2l = 1. - exp (- dt_evap / tau_v2l)
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

    cond = 0
    dep = 0
    reevap = 0
    sub = 0

    do k = ks, ke

        if (p1 (k) < p_min) cycle

        if (.not. do_warm_rain_mp) then

            ! -----------------------------------------------------------------------
            ! instant deposit all water vapor to cloud ice when temperature is super low
            ! -----------------------------------------------------------------------

            if (tz (k) < t_min) then
                sink = dim (qv (k), 1.e-7)
                dep = dep + sink * dp1 (k)
                qv (k) = qv (k) - sink
                qi (k) = qi (k) + sink
                q_sol (k) = q_sol (k) + sink
                tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / &
                     (one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice)
                if (do_qa) qa (k) = 1. ! air fully saturated; 100 % cloud cover
                cycle
            endif

            ! -----------------------------------------------------------------------
            ! instant evaporation / sublimation of all clouds if rh < rh_adj -- > cloud free
            ! -----------------------------------------------------------------------
            ! rain water is handled in warm - rain process.
            qpz = qv (k) + ql (k) + qi (k)
            tin = (te8 (k) - lv00 * qpz + li00 * (qs (k) + qg (k))) / &
                 (one_r8 + qpz * c1_vap + qr (k) * c1_liq + (qs (k) + qg (k)) * c1_ice)
            if (tin > t_sub + 6.) then
                rh = qpz / iqs1 (tin, den (k))
                if (rh < rh_adj) then ! qpz / rh_adj < qs
                    reevap = reevap + ql (k) * dp1 (k)
                    sub = sub + qi (k) * dp1 (k)
                    tz (k) = tin
                    qv (k) = qpz
                    ql (k) = 0.
                    qi (k) = 0.
                    cycle ! cloud free
                endif
            endif

        endif

        ! -----------------------------------------------------------------------
        ! cloud water < -- > vapor adjustment:
        ! -----------------------------------------------------------------------

        tin = tz (k)
        rh_tem = qpz / iqs1 (tin, den (k))
        qsw = wqs2 (tin, den (k), dwsdt)
        dq0 = qsw - qv (k)
        if (use_rhc_cevap) then
            evap = 0.
            if (rh_tem .lt. rhc_cevap) then
                if (dq0 > 0.) then ! evaporation
                    factor = min (1., fac_l2v * (10. * dq0 / qsw)) ! the rh dependent factor = 1 at 90%
                    evap = min (ql (k), factor * dq0 / (1. + tcp3 (k) * dwsdt))
                    reevap = reevap + evap * dp1 (k)
                elseif (do_cond_timescale) then
                    factor = min (1., fac_v2l * (10. * (- dq0) / qsw))
                    evap = - min (qv (k), factor * (- dq0) / (1. + tcp3 (k) * dwsdt))
                    cond = cond - evap * dp1 (k)
                else ! condensate all excess vapor into cloud water
                    evap = dq0 / (1. + tcp3 (k) * dwsdt)
                    cond = cond - evap * dp1 (k)
                endif
            endif
        else
            if (dq0 > 0.) then ! evaporation
                factor = min (1., fac_l2v * (10. * dq0 / qsw)) ! the rh dependent factor = 1 at 90%
                evap = min (ql (k), factor * dq0 / (1. + tcp3 (k) * dwsdt))
                reevap = reevap + evap * dp1 (k)
            elseif (do_cond_timescale) then
                factor = min (1., fac_v2l * (10. * (- dq0) / qsw))
                evap = - min (qv (k), factor * (- dq0) / (1. + tcp3 (k) * dwsdt))
                cond = cond - evap * dp1 (k)
            else ! condensate all excess vapor into cloud water
                evap = dq0 / (1. + tcp3 (k) * dwsdt)
                cond = cond - evap * dp1 (k)
            endif
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

        if (.not. do_warm_rain_mp) then

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

            if (do_sat_adj) then
                dt_pisub = 0.5 * dts
            else
                dt_pisub = dts
                tc = tice - tz (k)
                if (ql (k) > qrmin .and. tc > 0.1) then
                    sink = 100. / (rhow * ccn (k)) * dts * (exp (0.66 * tc) - 1.) * ql (k) ** 2
                    sink = min (ql (k), tc / icpk (k), sink)
                    ql (k) = ql (k) - sink
                    qi (k) = qi (k) + sink
                    q_liq (k) = q_liq (k) - sink
                    q_sol (k) = q_sol (k) + sink
                    cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
                    tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / cvm (k)
                endif ! significant ql existed
            endif

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
                    if (.not. prog_ccn) then
                        if (inflag .eq. 1) &
                            ! hong et al., 2004
                            cin (k) = 5.38e7 * exp (0.75 * log (qi (k) * den (k)))
                        if (inflag .eq. 2) &
                            ! meyers et al., 1992
                            cin (k) = exp (-2.80 + 0.262 * (tice - tz (k))) * 1000.0 ! convert from L^-1 to m^-3
                        if (inflag .eq. 3) &
                            ! meyers et al., 1992
                            cin (k) = exp (-0.639 + 12.96 * (qv (k) / qsi - 1.0)) * 1000.0 ! convert from L^-1 to m^-3
                        if (inflag .eq. 4) &
                            ! cooper, 1986
                            cin (k) = 5.e-3 * exp (0.304 * (tice - tz (k))) * 1000.0 ! convert from L^-1 to m^-3
                        if (inflag .eq. 5) &
                            ! flecther, 1962
                            cin (k) = 1.e-5 * exp (0.5 * (tice - tz (k))) * 1000.0 ! convert from L^-1 to m^-3
                    endif
                    pidep = dt_pisub * dq * 4.0 * 11.9 * exp (0.5 * log (qi (k) * den (k) * cin (k))) &
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
                    dep = dep + sink * dp1 (k)
                else ! ice -- > vapor
                    pidep = pidep * min (1., dim (tz (k), t_sub) * 0.2)
                    sink = max (pidep, sink, - qi (k))
                    sub = sub - sink * dp1 (k)
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
                    sub = sub + pssub * dp1 (k)
                else
                    if (tz (k) > tice) then
                        pssub = 0. ! no deposition
                    else
                        pssub = max (pssub, dq, (tz (k) - tice) / tcpk (k))
                    endif
                    dep = dep - pssub * dp1 (k)
                endif
                qs (k) = qs (k) - pssub
                qv (k) = qv (k) + pssub
                q_sol (k) = q_sol (k) - pssub
                cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
                tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / cvm (k)
                tcpk (k) = (li20 + (d1_vap + d1_ice) * tz (k)) / cvm (k)
            endif

            ! -----------------------------------------------------------------------
            ! sublimation / deposition of graupel
            ! this process happens for all temp rage
            ! -----------------------------------------------------------------------

            if (qg (k) > qrmin) then
                qsi = iqs2 (tz (k), den (k), dqsdt)
                qden = qg (k) * den (k)
                tmp = exp (0.6875 * log (qden))
                tsq = tz (k) * tz (k)
                dq = (qsi - qv (k)) / (1. + tcpk (k) * dqsdt)
                pgsub = cgsub (1) * tsq * (cgsub (2) * sqrt (qden) + cgsub (3) * tmp / &
                    sqrt (sqrt (den (k)))) / (cgsub (4) * tsq + cgsub (5) * qsi * den (k))
                pgsub = (qsi - qv (k)) * dts * pgsub
                if (pgsub > 0.) then ! qs -- > qv, sublimation
                    pgsub = min (pgsub * min (1., dim (tz (k), t_sub) * 0.2), qg (k))
                    sub = sub + pgsub * dp1 (k)
                else
                    if (tz (k) > tice) then
                        pgsub = 0. ! no deposition
                    else
                        pgsub = max (pgsub, dq, (tz (k) - tice) / tcpk (k))
                    endif
                    dep = dep - pgsub * dp1 (k)
                endif
                qg (k) = qg (k) - pgsub
                qv (k) = qv (k) + pgsub
                q_sol (k) = q_sol (k) - pgsub
                cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
                tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / cvm (k)
                tcpk (k) = (li20 + (d1_vap + d1_ice) * tz (k)) / cvm (k)
            endif

        endif

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
        ! icloud_f = 3: revision of icloud = 0
        ! -----------------------------------------------------------------------

        if (use_xr_cloud) then ! xu and randall cloud scheme (1996)
            if (rh >= 1.0) then
                qa (k) = 1.0
            elseif (rh > rh_thres .and. q_cond (k) > 1.e-6) then
                qa (k) = rh ** xr_a * (1.0 - exp (- xr_b * max (0.0, q_cond (k)) / &
                    max (1.e-5, (max (1.e-10, 1.0 - rh) * qstar) ** xr_c)))
                qa (k) = max (0.0, min (1., qa (k)))
            else
                qa (k) = 0.0
            endif
        elseif (use_park_cloud) then ! park et al. 2016 (mon. wea. review)
            if (q_cond (k) > 1.e-6) then
                qa (k) = 1. / 50. * (5.77 * (100. - gsize / 1000.) * max (0.0, q_cond (k) * 1000.) ** 1.07 + &
                    4.82 * (gsize / 1000. - 50.) * max (0.0, q_cond (k) * 1000.) ** 0.94)
                qa (k) = qa (k) * (0.92 / 0.96 * q_liq (k) / q_cond (k) + 1.0 / 0.96 * q_sol (k) / q_cond (k))
                qa (k) = max (0.0, min (1., qa (k)))
            else
                qa (k) = 0.0
            endif
        elseif (use_gi_cloud) then ! gultepe and isaac (2007)
            sigma = 0.28 + max (0.0, q_cond (k) * 1000.) ** 0.49
            gam = max (0.0, q_cond (k) * 1000.) / sigma
            if (gam < 0.18) then
                qa10 = 0.
            elseif (gam > 2.0) then
                qa10 = 1.0
            else
                qa10 = - 0.1754 + 0.9811 * gam - 0.2223 * gam ** 2 + 0.0104 * gam ** 3
                qa10 = max (0.0, min (1., qa10))
            endif
            if (gam < 0.12) then
                qa100 = 0.
            elseif (gam > 1.85) then
                qa100 = 1.0
            else
                qa100 = - 0.0913 + 0.7213 * gam + 0.1060 * gam ** 2 - 0.0946 * gam ** 3
                qa100 = max (0.0, min (1., qa100))
            endif
            qa (k) = qa10 + (log10 (gsize / 1000.) - 1) * (qa100 - qa10)
            qa (k) = max (0.0, min (1., qa (k)))
        else
            if (rh > rh_thres .and. qpz > 1.e-6) then

                dq = h_var * qpz
                if (do_cld_adj) then
                    q_plus = qpz + dq * f_dq_p * min(1.0, max(0.0, (p1 (k) - 200.e2) / (1000.e2 - 200.e2)))
                else
                    q_plus = qpz + dq * f_dq_p
                endif
                q_minus = qpz - dq * f_dq_m

                if (icloud_f .eq. 2) then
                    if (qstar < qpz) then
                        qa (k) = 1.
                    else
                        qa (k) = 0.
                    endif
                elseif (icloud_f .eq. 3) then
                    if (qstar < qpz) then
                        qa (k) = 1.
                    else
                        if (qstar < q_plus) then
                            qa (k) = (q_plus - qstar) / (dq * f_dq_p)
                        else
                            qa (k) = 0.
                        endif
                        ! impose minimum cloudiness if substantial q_cond (k) exist
                        if (q_cond (k) > 1.e-6) then
                            qa (k) = max (cld_min, qa (k))
                        endif
                        qa (k) = min (1., qa (k))
                    endif
                else
                    if (qstar < q_minus) then
                        qa (k) = 1.
                    else
                        if (qstar < q_plus) then
                            if (icloud_f .eq. 0) then
                                qa (k) = (q_plus - qstar) / (dq * f_dq_p + dq * f_dq_m)
                            else
                                qa (k) = (q_plus - qstar) / ((dq * f_dq_p + dq * f_dq_m) * (1. - q_cond (k)))
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
        endif

    enddo

end subroutine subgrid_z_proc

! =======================================================================
! rain evaporation
! =======================================================================

subroutine revap_rac1 (hydrostatic, is, ie, dt, tz, qv, ql, qr, qi, qs, qg, den, hvar)

    implicit none

    logical, intent (in) :: hydrostatic

    integer, intent (in) :: is, ie

    real, intent (in) :: dt ! time step (s)

    real, intent (in), dimension (is:ie) :: den, hvar, qi, qs, qg

    real, intent (inout), dimension (is:ie) :: tz, qv, qr, ql

    real, dimension (is:ie) :: lcp2, denfac, q_liq, q_sol, cvm, lhl

    real :: dqv, qsat, dqsdt, evap, qden, q_plus, q_minus, sink
    real :: tin, t2, qpz, dq, dqh

    integer :: i

    ! -----------------------------------------------------------------------
    ! define latend heat coefficient
    ! -----------------------------------------------------------------------

    do i = is, ie
        lhl (i) = lv00 + d0_vap * tz (i)
        q_liq (i) = ql (i) + qr (i)
        q_sol (i) = qi (i) + qs (i) + qg (i)
        cvm (i) = c_air + qv (i) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
        lcp2 (i) = lhl (i) / cvm (i)
        ! denfac (i) = sqrt (sfcrho / den (i))
    enddo

    do i = is, ie
        if (qr (i) > qrmin .and. tz (i) > t_wfr) then
            qpz = qv (i) + ql (i)
            tin = tz (i) - lcp2 (i) * ql (i) ! presence of clouds suppresses the rain evap
            qsat = wqs2 (tin, den (i), dqsdt)
            dqh = max (ql (i), hvar (i) * max (qpz, qcmin))
            dqv = qsat - qv (i)
            q_minus = qpz - dqh
            q_plus = qpz + dqh

            ! -----------------------------------------------------------------------
            ! qsat must be > q_minus to activate evaporation
            ! qsat must be < q_plus to activate accretion
            ! -----------------------------------------------------------------------

            ! -----------------------------------------------------------------------
            ! rain evaporation
            ! -----------------------------------------------------------------------

            if (dqv > qvmin .and. qsat > q_minus) then
                if (qsat > q_plus) then
                    dq = qsat - qpz
                else
                    ! q_minus < qsat < q_plus
                    ! dq == dqh if qsat == q_minus
                    dq = 0.25 * (q_minus - qsat) ** 2 / dqh
                endif
                qden = qr (i) * den (i)
                t2 = tin * tin
                evap = crevp (1) * t2 * dq * (crevp (2) * sqrt (qden) + crevp (3) * exp (0.725 * log (qden))) &
                     / (crevp (4) * t2 + crevp (5) * qsat * den (i))
                evap = min (qr (i), dt * evap, dqv / (1. + lcp2 (i) * dqsdt))
                qr (i) = qr (i) - evap
                qv (i) = qv (i) + evap
                q_liq (i) = q_liq (i) - evap
                cvm (i) = c_air + qv (i) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
                tz (i) = tz (i) - evap * lhl (i) / cvm (i)
            endif

            ! -----------------------------------------------------------------------
            ! accretion: pracc
            ! -----------------------------------------------------------------------

            if (qr (i) > qrmin .and. ql (i) > 1.e-8 .and. qsat < q_plus) then
                denfac (i) = sqrt (sfcrho / den (i))
                sink = dt * denfac (i) * cracw * exp (0.95 * log (qr (i) * den (i)))
                sink = sink / (1. + sink) * ql (i)
                ql (i) = ql (i) - sink
                qr (i) = qr (i) + sink
            endif
        endif
    enddo

end subroutine revap_rac1

! =======================================================================
! compute terminal fall speed
! consider cloud ice, snow, and graupel's melting during fall
! =======================================================================

subroutine terminal_fall (dtm, ks, ke, tz, qv, ql, qr, qg, qs, qi, dz, dp, &
        den, vtg, vts, vti, r1, g1, s1, i1, m1_sol, w1, dte)

    implicit none

    integer, intent (in) :: ks, ke
    real, intent (in) :: dtm ! time step (s)
    real, intent (in), dimension (ks:ke) :: vtg, vts, vti, den, dp, dz
    real (kind = r_grid), intent (inout), dimension (ks:ke) :: tz
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qg, qs, qi, m1_sol, w1
    real (kind = r_grid), intent (inout) :: dte
    real, intent (out) :: r1, g1, s1, i1
    ! local:
    real, dimension (ks:ke + 1) :: ze, zt
    real :: qsat, dqsdt, dt5, evap, dtime
    real :: factor, frac
    real :: tmp, precip, tc, sink
    real, dimension (ks:ke) :: lcpk, icpk, cvm, q_liq, q_sol
    real, dimension (ks:ke) :: m1, dm
    real (kind = r_grid), dimension (ks:ke) :: te1, te2
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
    ! melting of falling cloud ice into cloud water and rain
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

        ! -----------------------------------------------------------------------
        ! energy loss during sedimentation
        ! -----------------------------------------------------------------------

        if (consv_checker) then
            do k = ks, ke
                te1 (k) = one_r8 + qv (k) * c1_vap + (ql (k) + qr (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
                te1 (k) = rgrav * te1 (k) * c_air * tz (k) * dp (k)
            enddo
        endif

        if (use_ppm_ice) then
            call lagrangian_fall_ppm (ks, ke, zs, ze, zt, dp, qi, i1, m1_sol, mono_prof)
        else
            call implicit_fall (dtm, ks, ke, ze, vti, dp, qi, i1, m1_sol)
        endif

        ! -----------------------------------------------------------------------
        ! energy loss during sedimentation
        ! -----------------------------------------------------------------------

        if (consv_checker) then
            do k = ks, ke
                te2 (k) = one_r8 + qv (k) * c1_vap + (ql (k) + qr (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
                te2 (k) = rgrav * te2 (k) * c_air * tz (k) * dp (k)
            enddo
            dte = dte + sum (te1) - sum (te2)
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

        ! -----------------------------------------------------------------------
        ! energy loss during sedimentation
        ! -----------------------------------------------------------------------

        if (consv_checker) then
            do k = ks, ke
                te1 (k) = one_r8 + qv (k) * c1_vap + (ql (k) + qr (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
                te1 (k) = rgrav * te1 (k) * c_air * tz (k) * dp (k)
            enddo
        endif

        if (use_ppm) then
            call lagrangian_fall_ppm (ks, ke, zs, ze, zt, dp, qs, s1, m1, mono_prof)
        else
            call implicit_fall (dtm, ks, ke, ze, vts, dp, qs, s1, m1)
        endif

        ! -----------------------------------------------------------------------
        ! energy loss during sedimentation
        ! -----------------------------------------------------------------------

        if (consv_checker) then
            do k = ks, ke
                te2 (k) = one_r8 + qv (k) * c1_vap + (ql (k) + qr (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
                te2 (k) = rgrav * te2 (k) * c_air * tz (k) * dp (k)
            enddo
            dte = dte + sum (te1) - sum (te2)
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

        ! -----------------------------------------------------------------------
        ! energy loss during sedimentation
        ! -----------------------------------------------------------------------

        if (consv_checker) then
            do k = ks, ke
                te1 (k) = one_r8 + qv (k) * c1_vap + (ql (k) + qr (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
                te1 (k) = rgrav * te1 (k) * c_air * tz (k) * dp (k)
            enddo
        endif

        if (use_ppm) then
            call lagrangian_fall_ppm (ks, ke, zs, ze, zt, dp, qg, g1, m1, mono_prof)
        else
            call implicit_fall (dtm, ks, ke, ze, vtg, dp, qg, g1, m1)
        endif

        ! -----------------------------------------------------------------------
        ! energy loss during sedimentation
        ! -----------------------------------------------------------------------

        if (consv_checker) then
            do k = ks, ke
                te2 (k) = one_r8 + qv (k) * c1_vap + (ql (k) + qr (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
                te2 (k) = rgrav * te2 (k) * c_air * tz (k) * dp (k)
            enddo
            dte = dte + sum (te1) - sum (te2)
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
                if (hd_icefall) then
                    ! heymsfield and donner, 1990, jas
                    vti (k) = vi_fac * 3.29 * (qi (k) * den (k)) ** 0.16
                else
                    ! deng and mace, 2008, grl
                    vti (k) = (3. + log10 (qi (k) * den (k))) * (tc (k) * (aa * tc (k) + bb) + cc) + dd * tc (k) + ee
                    vti (k) = vi0 * exp (log_10 * vti (k))
                endif
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
    hlts = hlv + hlf
    hltc = hlv
    hltf = hlf

    ch2o = c_liq
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

    es0 = e00
    ces0 = eps * es0

end subroutine setupm

! =======================================================================
! initialization of gfdl cloud microphysics
! =======================================================================

subroutine gfdl_cld_mp_init (input_nml_file, logunit)

    implicit none

    character (len = *), intent (in) :: input_nml_file (:)
    integer, intent (in) :: logunit

    logical :: exists

    read (input_nml_file, nml = gfdl_mp_nml)

    ! write version number and namelist to log file
    write (logunit, *) " ================================================================== "
    write (logunit, *) "gfdl_mp_mod"
    write (logunit, nml = gfdl_mp_nml)

    if (do_setup) then
        call setup_con
        call setupm
        do_setup = .false.
    endif

    g2 = 0.5 * grav
    log_10 = log (10.)

    if (do_warm_rain_mp) then
        t_wfr = t_min
    else
        t_wfr = t_ice - 40.0
    endif

    module_is_initialized = .true.

end subroutine gfdl_cld_mp_init

! =======================================================================
! end of gfdl cloud microphysics
! =======================================================================

subroutine gfdl_cld_mp_end

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

end subroutine gfdl_cld_mp_end

! =======================================================================
! qsmith table initialization
! =======================================================================

subroutine setup_con

    implicit none

    ! master = (mpp_pe () .eq.mpp_root_pe ())

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
    es = tablew (it) + (ap1 - it) * desw (it)
    wqs2 = es / (rvgas * ta * den)
    it = ap1 - 0.5
    ! finite diff, del_t = 0.1:
    dqdt = 10. * (desw (it) + (ap1 - it) * (desw (it + 1) - desw (it))) / (rvgas * ta * den)

end function wqs2

! =======================================================================
! compute the gradient of saturated specific humidity for table ii
! it is the same as "wqs2", but written as vector function
! =======================================================================

subroutine wqs2_vect (is, ie, ta, den, wqsat, dqdt)

    implicit none

    ! pure water phase; universal dry / moist formular using air density
    ! input "den" can be either dry or moist air density

    integer, intent (in) :: is, ie

    real, intent (in), dimension (is:ie) :: ta, den

    real, intent (out), dimension (is:ie) :: wqsat, dqdt

    real :: es, ap1, tmin

    integer :: i, it

    tmin = t_ice - 160.

    do i = is, ie
        ap1 = 10. * dim (ta (i), tmin) + 1.
        ap1 = min (2621., ap1)
        it = ap1
        es = tablew (it) + (ap1 - it) * desw (it)
        wqsat (i) = es / (rvgas * ta (i) * den (i))
        it = ap1 - 0.5
        ! finite diff, del_t = 0.1:
        dqdt (i) = 10. * (desw (it) + (ap1 - it) * (desw (it + 1) - desw (it))) / (rvgas * ta (i) * den (i))
    enddo

end subroutine wqs2_vect

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
! compute the gradient of saturated specific humidity for table iii
! =======================================================================

real function qs1d_moist (ta, qv, pa, dqdt)

    implicit none

    real, intent (in) :: ta, pa, qv

    real, intent (out) :: dqdt

    real :: es, ap1, tmin, eps10

    integer :: it

    tmin = table_ice - 160.
    eps10 = 10. * eps
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = table2 (it) + (ap1 - it) * des2 (it)
    qs1d_moist = eps * es * (1. + zvir * qv) / pa
    it = ap1 - 0.5
    dqdt = eps10 * (des2 (it) + (ap1 - it) * (des2 (it + 1) - des2 (it))) * (1. + zvir * qv) / pa

end function qs1d_moist

! =======================================================================
! compute the gradient of saturated specific humidity for table ii
! =======================================================================

real function wqsat2_moist (ta, qv, pa, dqdt)

    implicit none

    real, intent (in) :: ta, pa, qv

    real, intent (out) :: dqdt

    real :: es, ap1, tmin, eps10

    integer :: it

    tmin = table_ice - 160.
    eps10 = 10. * eps
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = tablew (it) + (ap1 - it) * desw (it)
    wqsat2_moist = eps * es * (1. + zvir * qv) / pa
    it = ap1 - 0.5
    dqdt = eps10 * (desw (it) + (ap1 - it) * (desw (it + 1) - desw (it))) * (1. + zvir * qv) / pa

end function wqsat2_moist

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

! =======================================================================
! compute the saturated specific humidity for table iii
! =======================================================================

real function qs1d_m (ta, qv, pa)

    implicit none

    real, intent (in) :: ta, pa, qv

    real :: es, ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = table2 (it) + (ap1 - it) * des2 (it)
    qs1d_m = eps * es * (1. + zvir * qv) / pa

end function qs1d_m

! =======================================================================
! computes the difference in saturation vapor * density * between water and ice
! =======================================================================

real function d_sat (ta, den)

    implicit none

    real, intent (in) :: ta, den

    real :: es_w, es_i, ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es_w = tablew (it) + (ap1 - it) * desw (it)
    es_i = table2 (it) + (ap1 - it) * des2 (it)
    d_sat = dim (es_w, es_i) / (rvgas * ta * den) ! take positive difference

end function d_sat

! =======================================================================
! compute the saturated water vapor pressure for table ii
! =======================================================================

real function esw_table (ta)

    implicit none

    real, intent (in) :: ta

    real :: ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    esw_table = tablew (it) + (ap1 - it) * desw (it)

end function esw_table

! =======================================================================
! compute the saturated water vapor pressure for table iii
! =======================================================================

real function es2_table (ta)

    implicit none

    real, intent (in) :: ta

    real :: ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es2_table = table2 (it) + (ap1 - it) * des2 (it)

end function es2_table

! =======================================================================
! compute the saturated water vapor pressure for table ii
! =======================================================================

subroutine esw_table1d (ta, es, n)

    implicit none

    integer, intent (in) :: n

    real, intent (in) :: ta (n)

    real, intent (out) :: es (n)

    real :: ap1, tmin

    integer :: i, it

    tmin = table_ice - 160.

    do i = 1, n
        ap1 = 10. * dim (ta (i), tmin) + 1.
        ap1 = min (2621., ap1)
        it = ap1
        es (i) = tablew (it) + (ap1 - it) * desw (it)
    enddo

end subroutine esw_table1d

! =======================================================================
! compute the saturated water vapor pressure for table iii
! =======================================================================

subroutine es2_table1d (ta, es, n)

    implicit none

    integer, intent (in) :: n

    real, intent (in) :: ta (n)

    real, intent (out) :: es (n)

    real :: ap1, tmin

    integer :: i, it

    tmin = table_ice - 160.

    do i = 1, n
        ap1 = 10. * dim (ta (i), tmin) + 1.
        ap1 = min (2621., ap1)
        it = ap1
        es (i) = table2 (it) + (ap1 - it) * des2 (it)
    enddo

end subroutine es2_table1d

! =======================================================================
! compute the saturated water vapor pressure for table iv
! =======================================================================

subroutine es3_table1d (ta, es, n)

    implicit none

    integer, intent (in) :: n

    real, intent (in) :: ta (n)

    real, intent (out) :: es (n)

    real :: ap1, tmin

    integer :: i, it

    tmin = table_ice - 160.

    do i = 1, n
        ap1 = 10. * dim (ta (i), tmin) + 1.
        ap1 = min (2621., ap1)
        it = ap1
        es (i) = table3 (it) + (ap1 - it) * des3 (it)
    enddo

end subroutine es3_table1d

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
! compute the saturated specific humidity for table
! note: this routine is based on "moist" mixing ratio
! =======================================================================

real function qs_blend (t, p, q)

    implicit none

    real, intent (in) :: t, p, q

    real :: es, ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (t, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = table (it) + (ap1 - it) * des (it)
    qs_blend = eps * es * (1. + zvir * q) / p

end function qs_blend

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
! compute the saturated specific humidity and the gradient of saturated specific humidity
! input t in deg k, p in pa; p = rho rdry tv, moist pressure
! =======================================================================

subroutine qsmith (im, km, ks, t, p, q, qs, dqdt)

    implicit none

    integer, intent (in) :: im, km, ks

    real, intent (in), dimension (im, km) :: t, p, q

    real, intent (out), dimension (im, km) :: qs

    real, intent (out), dimension (im, km), optional :: dqdt

    real :: eps10, ap1, tmin

    real, dimension (im, km) :: es

    integer :: i, k, it

    tmin = table_ice - 160.
    eps10 = 10. * eps

    if (.not. tables_are_initialized) then
        call qsmith_init
    endif

    do k = ks, km
        do i = 1, im
            ap1 = 10. * dim (t (i, k), tmin) + 1.
            ap1 = min (2621., ap1)
            it = ap1
            es (i, k) = table (it) + (ap1 - it) * des (it)
            qs (i, k) = eps * es (i, k) * (1. + zvir * q (i, k)) / p (i, k)
        enddo
    enddo

    if (present (dqdt)) then
        do k = ks, km
            do i = 1, im
                ap1 = 10. * dim (t (i, k), tmin) + 1.
                ap1 = min (2621., ap1) - 0.5
                it = ap1
                dqdt (i, k) = eps10 * (des (it) + (ap1 - it) * (des (it + 1) - des (it))) * (1. + zvir * q (i, k)) / p (i, k)
            enddo
        enddo
    endif

end subroutine qsmith

! =======================================================================
! fix negative water species
! this is designed for 6 - class micro - physics schemes
! =======================================================================

subroutine neg_adj (ks, ke, pt, dp, qv, ql, qr, qi, qs, qg, cond)

    implicit none

    integer, intent (in) :: ks, ke
    real, intent (in), dimension (ks:ke) :: dp
    real (kind = r_grid), intent (inout), dimension (ks:ke) :: pt
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (out) :: cond

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

    cond = 0

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
        if (qg (k) < 0.) then
            qr (k) = qr (k) + qg (k)
            pt (k) = pt (k) - qg (k) * icpk (k) ! heating
            qg (k) = 0.
        endif

        ! -----------------------------------------------------------------------
        ! liquid phase:
        ! -----------------------------------------------------------------------

        ! if rain < 0, borrow from cloud water
        if (qr (k) < 0.) then
            ql (k) = ql (k) + qr (k)
            qr (k) = 0.
        endif
        ! if cloud water < 0, borrow from water vapor
        if (ql (k) < 0.) then
            cond = cond - ql (k) * dp (k)
            qv (k) = qv (k) + ql (k)
            pt (k) = pt (k) - ql (k) * lcpk (k) ! heating
            ql (k) = 0.
        endif

    enddo

    ! -----------------------------------------------------------------------
    ! fix water vapor; borrow from below
    ! -----------------------------------------------------------------------

    do k = ks, ke - 1
        if (qv (k) < 0.) then
            qv (k + 1) = qv (k + 1) + qv (k) * dp (k) / dp (k + 1)
            qv (k) = 0.
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! bottom layer; borrow from above
    ! -----------------------------------------------------------------------

    if (qv (ke) < 0. .and. qv (ke - 1) > 0.) then
        dq = min (- qv (ke) * dp (ke), qv (ke - 1) * dp (ke - 1))
        qv (ke - 1) = qv (ke - 1) - dq / dp (ke - 1)
        qv (ke) = qv (ke) + dq / dp (ke)
    endif

end subroutine neg_adj

end module gfdl_cld_mp_mod
