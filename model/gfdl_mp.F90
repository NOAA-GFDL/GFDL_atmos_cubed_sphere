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
! GFDL Cloud Microphysics Package (GFDL MP) Version 3
! The algorithms are originally derived from Lin et al. (1983).
! Most of the key elements have been simplified / improved.
! This code at this stage bears little to no similarity to the original Lin MP in ZETAC.
! Developers: Linjiong Zhou and the GFDL FV3 Team
! References:
! Version 0: Chen and Lin (2011 doi: 10.1029/2011GL047629, 2013 doi: 10.1175/JCLI-D-12-00061.1)
! Version 1: Zhou et al. (2019 doi: 10.1175/BAMS-D-17-0246.1)
! Version 2: Harris et al. (2020 doi: 10.1029/2020MS002223), Zhou et al. (2022 doi: 10.25923/pz3c-8b96)
! Version 3: Zhou et al. (2022 doi: 10.1029/2021MS002971)
! =======================================================================

module gfdl_mp_mod
    
    implicit none
    
    private
    
    ! -----------------------------------------------------------------------
    ! interface functions
    ! -----------------------------------------------------------------------
    
    interface wqs
        procedure wes_t
        procedure wqs_trho
        procedure wqs_ptqv
    end interface wqs
    
    interface mqs
        procedure mes_t
        procedure mqs_trho
        procedure mqs_ptqv
    end interface mqs
    
    interface iqs
        procedure ies_t
        procedure iqs_trho
        procedure iqs_ptqv
    end interface iqs
    
    interface mhc
        procedure mhc3
        procedure mhc4
        procedure mhc6
    end interface mhc
    
    interface wet_bulb
        procedure wet_bulb_dry
        procedure wet_bulb_moist
    end interface wet_bulb
    
    ! -----------------------------------------------------------------------
    ! public subroutines, functions, and variables
    ! -----------------------------------------------------------------------
    
    public :: gfdl_mp_init
    public :: gfdl_mp_driver
    public :: gfdl_mp_end
    public :: fast_sat_adj, cld_eff_rad, rad_ref
    public :: qs_init, wqs, mqs, mqs3d
    public :: c_liq, c_ice, rhow, wet_bulb
    public :: cv_air, cv_vap, mtetw
    public :: hlv, hlf, tice
    
    ! -----------------------------------------------------------------------
    ! precision definition
    ! -----------------------------------------------------------------------
    
    integer, parameter :: r8 = 8 ! double precision
    
    ! -----------------------------------------------------------------------
    ! initialization conditions
    ! -----------------------------------------------------------------------
    
    logical :: tables_are_initialized = .false. ! initialize satuation tables
    
    ! -----------------------------------------------------------------------
    ! physics constants
    ! -----------------------------------------------------------------------
    
    real, parameter :: grav = 9.80665 ! acceleration due to gravity (m/s^2), ref: IFS
    
    real, parameter :: rgrav = 1.0 / grav ! inversion of gravity acceleration (s^2/m)
    
    real, parameter :: pi = 4.0 * atan (1.0) ! ratio of circle circumference to diameter

    real, parameter :: boltzmann = 1.38064852e-23 ! boltzmann constant (J/K)
    real, parameter :: avogadro = 6.02214076e23 ! avogadro number (1/mol)
    real, parameter :: runiver = avogadro * boltzmann ! 8.314459727525675, universal gas constant (J/K/mol)
    real, parameter :: mmd = 2.89644e-2 ! dry air molar mass (kg/mol), ref: IFS
    real, parameter :: mmv = 1.80153e-2 ! water vapor molar mass (kg/mol), ref: IFS
    
    real, parameter :: rdgas = 287.05 ! gas constant for dry air (J/kg/K): ref: GFDL, GFS
    real, parameter :: rvgas = 461.50 ! gas constant for water vapor (J/kg/K): ref: GFDL, GFS
    !real, parameter :: rdgas = runiver / mmd ! 287.0578961596192, gas constant for dry air (J/kg/K)
    !real, parameter :: rvgas = runiver / mmv ! 461.52213549181386, gas constant for water vapor (J/kg/K)

    real, parameter :: zvir = rvgas / rdgas - 1. ! 0.6077667316114637
    real, parameter :: eps = rdgas / rvgas ! 0.6219934994582882
    real, parameter :: epsm1 = rdgas / rvgas - 1. ! -0.3780065005417118
    
    real, parameter :: tice = 273.15 ! freezing temperature (K): ref: GFDL, GFS
    !real, parameter :: tice = 273.16 ! freezing temperature (K), ref: IFS
    
    real, parameter :: cp_air = 1004.6 ! heat capacity of dry air at constant pressure (J/kg/K): ref: GFDL, GFS
    real, parameter :: cv_air = cp_air - rdgas ! 717.55, heat capacity of dry air at constant volume (J/kg/K): ref: GFDL, GFS
    !real, parameter :: cp_air = 7. / 2. * rdgas ! 1004.7026365586671, heat capacity of dry air at constant pressure (J/kg/K)
    !real, parameter :: cv_air = 5. / 2. * rdgas ! 717.644740399048, heat capacity of dry air at constant volume (J/kg/K)
    real, parameter :: cp_vap = 4.0 * rvgas ! 1846.0885419672554, heat capacity of water vapor at constnat pressure (J/kg/K)
    real, parameter :: cv_vap = 3.0 * rvgas ! 1384.5664064754415, heat capacity of water vapor at constant volume (J/kg/K)
    
    real, parameter :: c_ice = 2.106e3 ! heat capacity of ice at 0 deg C (J/kg/K), ref: IFS
    real, parameter :: c_liq = 4.218e3 ! heat capacity of water at 0 deg C (J/kg/K), ref: IFS
    
    real, parameter :: dc_vap = cp_vap - c_liq ! - 2371.9114580327446, isobaric heating / cooling (J/kg/K)
    real, parameter :: dc_ice = c_liq - c_ice ! 2112.0, isobaric heating / colling (J/kg/K)
    real, parameter :: d2_ice = cp_vap - c_ice ! - 259.9114580327446, isobaric heating / cooling (J/kg/K)
    
    real, parameter :: hlv = 2.5e6 ! latent heat of evaporation at 0 deg C (J/kg): ref: GFDL, GFS
    real, parameter :: hlf = 3.3358e5 ! latent heat of fusion at 0 deg C (J/kg): ref: GFDL, GFS
    !real, parameter :: hlv = 2.5008e6 ! latent heat of evaporation at 0 deg C (J/kg), ref: IFS
    !real, parameter :: hlf = 3.345e5 ! latent heat of fusion at 0 deg C (J/kg), ref: IFS
    
    real, parameter :: visd = 1.717e-5 ! dynamics viscosity of air at 0 deg C and 1000 hPa (Mason, 1971) (kg/m/s)
    real, parameter :: visk = 1.35e-5 ! kinematic viscosity of air at 0 deg C  and 1000 hPa (Mason, 1971) (m^2/s)
    real, parameter :: vdifu = 2.25e-5 ! diffusivity of water vapor in air at 0 deg C  and 1000 hPa (Mason, 1971) (m^2/s)
    real, parameter :: tcond = 2.40e-2 ! thermal conductivity of air at 0 deg C  and 1000 hPa (Mason, 1971) (J/m/s/K)

    real, parameter :: rho0 = 1.0 ! reference air density (kg/m^3), ref: IFS
    real, parameter :: cdg = 3.15121 ! drag coefficient of graupel (Locatelli and Hobbs, 1974)
    real, parameter :: cdh = 0.5 ! drag coefficient of hail (Heymsfield and Wright, 2014)
    
    real (kind = r8), parameter :: lv0 = hlv - dc_vap * tice ! 3148711.3338762247, evaporation latent heat coeff. at 0 deg K (J/kg)
    real (kind = r8), parameter :: li0 = hlf - dc_ice * tice ! - 242413.92000000004, fussion latent heat coeff. at 0 deg K (J/kg)
    real (kind = r8), parameter :: li2 = lv0 + li0 ! 2906297.413876225, sublimation latent heat coeff. at 0 deg K (J/kg)
    
    real (kind = r8), parameter :: e00 = 611.21 ! saturation vapor pressure at 0 deg C (Pa), ref: IFS
    
    ! -----------------------------------------------------------------------
    ! predefined parameters
    ! -----------------------------------------------------------------------
    
    integer, parameter :: length = 2621 ! length of the saturation table
    
    real, parameter :: qcmin = 1.0e-15 ! min value for cloud condensates (kg/kg)
    real, parameter :: qfmin = 1.0e-8 ! min value for sedimentation (kg/kg)
    
    real, parameter :: dz_min = 1.0e-2 ! used for correcting flipped height (m)
    
    real, parameter :: rhow = 1.0e3 ! density of cloud water (kg/m^3)
    real, parameter :: rhoi = 9.17e2 ! density of cloud ice (kg/m^3)
    real, parameter :: rhor = 1.0e3 ! density of rain (Lin et al. 1983) (kg/m^3)
    real, parameter :: rhos = 1.0e2 ! density of snow (Lin et al. 1983) (kg/m^3)
    real, parameter :: rhog = 4.0e2 ! density of graupel (Rutledge and Hobbs 1984) (kg/m^3)
    real, parameter :: rhoh = 9.17e2 ! density of hail (Lin et al. 1983) (kg/m^3)
    
    real, parameter :: dt_fr = 8.0 ! t_wfr - dt_fr: minimum temperature water can exist (Moore and Molinero 2011)
    
    real (kind = r8), parameter :: one_r8 = 1.0 ! constant 1
    
    ! -----------------------------------------------------------------------
    ! namelist parameters
    ! -----------------------------------------------------------------------
    
    integer :: ntimes = 1 ! cloud microphysics sub cycles
    
    integer :: cfflag = 1 ! cloud fraction scheme
    ! 1: GFDL cloud scheme
    ! 2: Xu and Randall (1996)
    ! 3: Park et al. (2016)
    ! 4: Gultepe and Isaac (2007)
    
    integer :: icloud_f = 0 ! GFDL cloud scheme
    ! 0: subgrid variability based scheme
    ! 1: same as 0, but for old fvgfs implementation
    ! 2: binary cloud scheme
    ! 3: extension of 0
    
    integer :: irain_f = 0 ! cloud water to rain auto conversion scheme
    ! 0: subgrid variability based scheme
    ! 1: no subgrid varaibility
    
    integer :: inflag = 1 ! ice nucleation scheme
    ! 1: Hong et al. (2004)
    ! 2: Meyers et al. (1992)
    ! 3: Meyers et al. (1992)
    ! 4: Cooper (1986)
    ! 5: Fletcher (1962)
    
    integer :: igflag = 3 ! ice generation scheme
    ! 1: WSM6
    ! 2: WSM6 with 0 at 0 C
    ! 3: WSM6 with 0 at 0 C and fixed value at - 10 C
    ! 4: combination of 1 and 3
    
    integer :: ifflag = 1 ! ice fall scheme
    ! 1: Deng and Mace (2008)
    ! 2: Heymsfield and Donner (1990)
    
    integer :: rewflag = 1 ! cloud water effective radius scheme
    ! 1: Martin et al. (1994)
    ! 2: Martin et al. (1994), GFDL revision
    ! 3: Kiehl et al. (1994)
    ! 4: effective radius
    
    integer :: reiflag = 5 ! cloud ice effective radius scheme
    ! 1: Heymsfield and Mcfarquhar (1996)
    ! 2: Donner et al. (1997)
    ! 3: Fu (2007)
    ! 4: Kristjansson et al. (2000)
    ! 5: Wyser (1998)
    ! 6: Sun and Rikus (1999), Sun (2001)
    ! 7: effective radius
    
    integer :: rerflag = 1 ! rain effective radius scheme
    ! 1: effective radius
    
    integer :: resflag = 1 ! snow effective radius scheme
    ! 1: effective radius
    
    integer :: regflag = 1 ! graupel effective radius scheme
    ! 1: effective radius
    
    integer :: radr_flag = 1 ! radar reflectivity for rain
    ! 1: Mark Stoelinga (2005)
    ! 2: Smith et al. (1975), Tong and Xue (2005)
    ! 3: Marshall-Palmer formula (https://en.wikipedia.org/wiki/DBZ_(meteorology))
    
    integer :: rads_flag = 1 ! radar reflectivity for snow
    ! 1: Mark Stoelinga (2005)
    ! 2: Smith et al. (1975), Tong and Xue (2005)
    ! 3: Marshall-Palmer formula (https://en.wikipedia.org/wiki/DBZ_(meteorology))
    
    integer :: radg_flag = 1 ! radar reflectivity for graupel
    ! 1: Mark Stoelinga (2005)
    ! 2: Smith et al. (1975), Tong and Xue (2005)
    ! 3: Marshall-Palmer formula (https://en.wikipedia.org/wiki/DBZ_(meteorology))

    integer :: sedflag = 1 ! sedimentation scheme
    ! 1: implicit scheme
    ! 2: explicit scheme
    ! 3: lagrangian scheme
    ! 4: combined implicit and lagrangian scheme
    
    integer :: vdiffflag = 1 ! wind difference scheme in accretion
    ! 1: Wisner et al. (1972)
    ! 2: Mizuno (1990)
    ! 3: Murakami (1990)
    
    logical :: do_sedi_uv = .true. ! transport of horizontal momentum in sedimentation
    logical :: do_sedi_w = .true. ! transport of vertical momentum in sedimentation
    logical :: do_sedi_heat = .true. ! transport of heat in sedimentation
    logical :: do_sedi_melt = .true. ! melt cloud ice, snow, and graupel during sedimentation
    
    logical :: do_qa = .true. ! do inline cloud fraction
    logical :: rad_snow = .true. ! include snow in cloud fraciton calculation
    logical :: rad_graupel = .true. ! include graupel in cloud fraction calculation
    logical :: rad_rain = .true. ! include rain in cloud fraction calculation
    logical :: do_cld_adj = .false. ! do cloud fraction adjustment
    
    logical :: z_slope_liq = .true. ! use linear mono slope for autocconversions
    logical :: z_slope_ice = .true. ! use linear mono slope for autocconversions
    
    logical :: use_rhc_cevap = .false. ! cap of rh for cloud water evaporation
    logical :: use_rhc_revap = .false. ! cap of rh for rain evaporation
    
    logical :: const_vw = .false. ! if .ture., the constants are specified by v * _fac
    logical :: const_vi = .false. ! if .ture., the constants are specified by v * _fac
    logical :: const_vs = .false. ! if .ture., the constants are specified by v * _fac
    logical :: const_vg = .false. ! if .ture., the constants are specified by v * _fac
    logical :: const_vr = .false. ! if .ture., the constants are specified by v * _fac
    
    logical :: liq_ice_combine = .false. ! combine all liquid water, combine all solid water
    logical :: snow_grauple_combine = .true. ! combine snow and graupel
    
    logical :: prog_ccn = .false. ! do prognostic ccn (Yi Ming's method)
    
    logical :: fix_negative = .true. ! fix negative water species
    
    logical :: do_cond_timescale = .false. ! whether to apply a timescale to condensation
    
    logical :: do_hail = .false. ! use hail parameters instead of graupel
    
    logical :: consv_checker = .false. ! turn on energy and water conservation checker
    
    logical :: do_warm_rain_mp = .false. ! do warm rain cloud microphysics only
    
    logical :: do_wbf = .false. ! do Wegener Bergeron Findeisen process

    logical :: do_psd_water_fall = .false. ! calculate cloud water terminal velocity based on PSD
    logical :: do_psd_ice_fall = .false. ! calculate cloud ice terminal velocity based on PSD

    logical :: do_psd_water_num = .false. ! calculate cloud water number concentration based on PSD
    logical :: do_psd_ice_num = .false. ! calculate cloud ice number concentration based on PSD

    logical :: do_new_acc_water = .false. ! perform the new accretion for cloud water
    logical :: do_new_acc_ice = .false. ! perform the new accretion for cloud ice
    
    logical :: cp_heating = .false. ! update temperature based on constant pressure
    
    real :: mp_time = 150.0 ! maximum microphysics time step (s)
    
    real :: n0w_sig = 1.1 ! intercept parameter (significand) of cloud water (Lin et al. 1983) (1/m^4) (Martin et al. 1994)
    !real :: n0w_sig = 1.4 ! intercept parameter (significand) of cloud water (Lin et al. 1983) (1/m^4) (Martin et al. 1994)
    real :: n0i_sig = 1.3 ! intercept parameter (significand) of cloud ice (Lin et al. 1983) (1/m^4) (McFarquhar et al. 2015)
    !real :: n0i_sig = 9.4 ! intercept parameter (significand) of cloud ice (Lin et al. 1983) (1/m^4) (McFarquhar et al. 2015)
    real :: n0r_sig = 8.0 ! intercept parameter (significand) of rain (Lin et al. 1983) (1/m^4) (Marshall and Palmer 1948)
    real :: n0s_sig = 3.0 ! intercept parameter (significand) of snow (Lin et al. 1983) (1/m^4) (Gunn and Marshall 1958)
    real :: n0g_sig = 4.0 ! intercept parameter (significand) of graupel (Rutledge and Hobbs 1984) (1/m^4) (Houze et al. 1979)
    real :: n0h_sig = 4.0 ! intercept parameter (significand) of hail (Lin et al. 1983) (1/m^4) (Federer and Waldvogel 1975)
    
    real :: n0w_exp = 41 ! intercept parameter (exponent) of cloud water (Lin et al. 1983) (1/m^4) (Martin et al. 1994)
    !real :: n0w_exp = 91 ! intercept parameter (exponent) of cloud water (Lin et al. 1983) (1/m^4) (Martin et al. 1994)
    real :: n0i_exp = 18 ! intercept parameter (exponent) of cloud ice (Lin et al. 1983) (1/m^4) (McFarquhar et al. 2015)
    !real :: n0i_exp = 17 ! intercept parameter (exponent) of cloud ice (Lin et al. 1983) (1/m^4) (McFarquhar et al. 2015)
    real :: n0r_exp = 6 ! intercept parameter (exponent) of rain (Lin et al. 1983) (1/m^4) (Marshall and Palmer 1948)
    real :: n0s_exp = 6 ! intercept parameter (exponent) of snow (Lin et al. 1983) (1/m^4) (Gunn and Marshall 1958)
    real :: n0g_exp = 6 ! intercept parameter (exponent) of graupel (Rutledge and Hobbs 1984) (1/m^4) (Houze et al. 1979)
    real :: n0h_exp = 4 ! intercept parameter (exponent) of hail (Lin et al. 1983) (1/m^4) (Federer and Waldvogel 1975)
    
    real :: muw = 6.0 ! shape parameter of cloud water in Gamma distribution (Martin et al. 1994)
    !real :: muw = 16.0 ! shape parameter of cloud water in Gamma distribution (Martin et al. 1994)
    real :: mui = 3.35 ! shape parameter of cloud ice in Gamma distribution (McFarquhar et al. 2015)
    !real :: mui = 3.54 ! shape parameter of cloud ice in Gamma distribution (McFarquhar et al. 2015)
    real :: mur = 1.0 ! shape parameter of rain in Gamma distribution (Marshall and Palmer 1948)
    real :: mus = 1.0 ! shape parameter of snow in Gamma distribution (Gunn and Marshall 1958)
    real :: mug = 1.0 ! shape parameter of graupel in Gamma distribution (Houze et al. 1979)
    real :: muh = 1.0 ! shape parameter of hail in Gamma distribution (Federer and Waldvogel 1975)
    
    real :: alinw = 3.e7 ! "a" in Lin et al. (1983) for cloud water (Ikawa and Saito 1990)
    real :: alini = 7.e2 ! "a" in Lin et al. (1983) for cloud ice (Ikawa and Saita 1990)
    real :: alinr = 842.0 ! "a" in Lin et al. (1983) for rain (Liu and Orville 1969)
    real :: alins = 4.8 ! "a" in Lin et al. (1983) for snow (straka 2009)
    real :: aling = 1.0 ! "a" in Lin et al. (1983), similar to a, but for graupel (Pruppacher and Klett 2010)
    real :: alinh = 1.0 ! "a" in Lin et al. (1983), similar to a, but for hail (Pruppacher and Klett 2010)

    real :: blinw = 2.0 ! "b" in Lin et al. (1983) for cloud water (Ikawa and Saito 1990)
    real :: blini = 1.0 ! "b" in Lin et al. (1983) for cloud ice (Ikawa and Saita 1990)
    real :: blinr = 0.8 ! "b" in Lin et al. (1983) for rain (Liu and Orville 1969)
    real :: blins = 0.25 ! "b" in Lin et al. (1983) for snow (straka 2009)
    real :: bling = 0.5 ! "b" in Lin et al. (1983), similar to b, but for graupel (Pruppacher and Klett 2010)
    real :: blinh = 0.5 ! "b" in Lin et al. (1983), similar to b, but for hail (Pruppacher and Klett 2010)
    
    real :: tice_mlt = 273.16 ! can set ice melting temperature to 268 based on observation (Kay et al. 2016) (K)
    
    real :: t_min = 178.0 ! minimum temperature to freeze - dry all water vapor (K)
    real :: t_sub = 184.0 ! minimum temperature for sublimation of cloud ice (K)
    
    real :: rh_inc = 0.25 ! rh increment for complete evaporation of cloud water and cloud ice
    real :: rh_inr = 0.25 ! rh increment for minimum evaporation of rain
    real :: rh_ins = 0.25 ! rh increment for sublimation of snow
    
    real :: tau_r2g = 900.0 ! rain freezing to graupel time scale (s)
    real :: tau_i2s = 1000.0 ! cloud ice to snow autoconversion time scale (s)
    real :: tau_l2r = 900.0 ! cloud water to rain autoconversion time scale (s)
    real :: tau_v2l = 150.0 ! water vapor to cloud water condensation time scale (s)
    real :: tau_l2v = 300.0 ! cloud water to water vapor evaporation time scale (s)
    real :: tau_revp = 0.0 ! rain evaporation time scale (s)
    real :: tau_imlt = 1200.0 ! cloud ice melting time scale (s)
    real :: tau_smlt = 900.0 ! snow melting time scale (s)
    real :: tau_gmlt = 600.0 ! graupel melting time scale (s)
    real :: tau_wbf = 300.0 ! graupel melting time scale (s)
    
    real :: dw_land = 0.20 ! base value for subgrid deviation / variability over land
    real :: dw_ocean = 0.10 ! base value for subgrid deviation / variability over ocean
    
    real :: ccn_o = 90.0 ! ccn over ocean (1/cm^3)
    real :: ccn_l = 270.0 ! ccn over land (1/cm^3)
    
    real :: rthresh = 10.0e-6 ! critical cloud drop radius (micron) for autoconversion
    
    real :: cld_min = 0.05 ! minimum cloud fraction
    
    real :: qi_lim = 1.0 ! cloud ice limiter (0: no, 1: full, >1: extra) to prevent large ice build up
    
    real :: ql_mlt = 2.0e-3 ! maximum cloud water allowed from melted cloud ice (kg/kg)
    real :: qs_mlt = 1.0e-6 ! maximum cloud water allowed from melted snow (kg/kg)
    
    real :: ql_gen = 1.0e-3 ! maximum cloud water generation during remapping step (kg/kg)
    
    real :: ql0_max = 2.0e-3 ! maximum cloud water value (autoconverted to rain) (kg/kg)
    real :: qi0_max = 1.0e-4 ! maximum cloud ice value (autoconverted to snow) (kg/m^3)
    
    real :: qi0_crt = 1.0e-4 ! cloud ice to snow autoconversion threshold (kg/m^3)
    real :: qs0_crt = 1.0e-3 ! snow to graupel autoconversion threshold (0.6e-3 in Purdue Lin scheme) (kg/m^3)
    
    real :: c_paut = 0.55 ! cloud water to rain autoconversion efficiency
    real :: c_psacw = 1.0 ! cloud water to snow accretion efficiency
    real :: c_psaci = 0.05 ! cloud ice to snow accretion efficiency (was 0.1 in ZETAC)
    real :: c_pracw = 0.8 ! cloud water to rain accretion efficiency
    real :: c_praci = 1.0 ! cloud ice to rain accretion efficiency
    real :: c_pgacw = 1.0 ! cloud water to graupel accretion efficiency
    real :: c_pgaci = 0.05 ! cloud ice to graupel accretion efficiency (was 0.1 in ZETAC)
    real :: c_pracs = 1.0 ! snow to rain accretion efficiency
    real :: c_psacr = 1.0 ! rain to snow accretion efficiency
    real :: c_pgacr = 1.0 ! rain to graupel accretion efficiency
    real :: c_pgacs = 0.01 ! snow to graupel accretion efficiency (was 0.1 in ZETAC)

    real :: is_fac = 0.2 ! cloud ice sublimation temperature factor
    real :: ss_fac = 0.2 ! snow sublimation temperature factor
    real :: gs_fac = 0.2 ! graupel sublimation temperature factor

    real :: rh_fac = 10.0 ! cloud water condensation / evaporation relative humidity factor

    real :: sed_fac = 1.0 ! coefficient for sedimentation fall, scale from 1.0 (implicit) to 0.0 (lagrangian)
    
    real :: vw_fac = 1.0
    real :: vi_fac = 1.0 ! IFS: if const_vi: 1 / 3
    real :: vs_fac = 1.0 ! IFS: if const_vs: 1.
    real :: vg_fac = 1.0 ! IFS: if const_vg: 2.
    real :: vr_fac = 1.0 ! IFS: if const_vr: 4.
    
    real :: vw_max = 0.01 ! maximum fall speed for cloud water (m/s)
    real :: vi_max = 0.5 ! maximum fall speed for cloud ice (m/s)
    real :: vs_max = 5.0 ! maximum fall speed for snow (m/s)
    real :: vg_max = 8.0 ! maximum fall speed for graupel (m/s)
    real :: vr_max = 12.0 ! maximum fall speed for rain (m/s)
    
    real :: xr_a = 0.25 ! p value in Xu and Randall (1996)
    real :: xr_b = 100.0 ! alpha_0 value in Xu and Randall (1996)
    real :: xr_c = 0.49 ! gamma value in Xu and Randall (1996)
    
    real :: te_err = 1.e-5 ! 64bit: 1.e-14, 32bit: 1.e-7; turn off to save computer time
    real :: tw_err = 1.e-8 ! 64bit: 1.e-14, 32bit: 1.e-7; turn off to save computer time
    
    real :: rh_thres = 0.75 ! minimum relative humidity for cloud fraction
    real :: rhc_cevap = 0.85 ! maximum relative humidity for cloud water evaporation
    real :: rhc_revap = 0.85 ! maximum relative humidity for rain evaporation
    
    real :: f_dq_p = 1.0 ! cloud fraction adjustment for supersaturation
    real :: f_dq_m = 1.0 ! cloud fraction adjustment for undersaturation
    
    real :: fi2s_fac = 1.0 ! maximum sink of cloud ice to form snow: 0-1
    real :: fi2g_fac = 1.0 ! maximum sink of cloud ice to form graupel: 0-1
    real :: fs2g_fac = 1.0 ! maximum sink of snow to form graupel: 0-1
    
    real :: beta = 1.22 ! defined in Heymsfield and Mcfarquhar (1996)
    
    real :: rewmin = 5.0, rewmax = 15.0 ! minimum and maximum effective radius for cloud water (micron)
    real :: reimin = 10.0, reimax = 150.0 ! minimum and maximum effective radius for cloud ice (micron)
    real :: rermin = 15.0, rermax = 10000.0 ! minimum and maximum effective radius for rain (micron)
    real :: resmin = 150.0, resmax = 10000.0 ! minimum and maximum effective radius for snow (micron)
    real :: regmin = 150.0, regmax = 10000.0 ! minimum and maximum effective radius for graupel
    !real :: rewmax = 15.0, rermin = 15.0 ! Kokhanovsky (2004)

    real :: rewfac = 1.0 ! this is a tuning parameter to compromise the inconsistency between
                         ! GFDL MP's PSD and cloud water radiative property's PSD assumption.
                         ! after the cloud water radiative property's PSD is rebuilt,
                         ! this parameter should be 1.0.
    real :: reifac = 1.0 ! this is a tuning parameter to compromise the inconsistency between
                         ! GFDL MP's PSD and cloud ice radiative property's PSD assumption.
                         ! after the cloud ice radiative property's PSD is rebuilt,
                         ! this parameter should be 1.0.
    
    ! -----------------------------------------------------------------------
    ! local shared variables
    ! -----------------------------------------------------------------------
    
    real :: acco (3, 10), acc (20)
    real :: cracs, csacr, cgacr, cgacs, csacw, craci, csaci, cgacw, cgaci, cracw
    real :: cssub (5), cgsub (5), crevp (5), cgfr (2), csmlt (4), cgmlt (4)
    
    real :: t_wfr, fac_rc, c_air, c_vap, d0_vap
    
    real (kind = r8) :: lv00, li00, li20, cpaut
    real (kind = r8) :: d1_vap, d1_ice, c1_vap, c1_liq, c1_ice
    real (kind = r8) :: normw, normr, normi, norms, normg, normh
    real (kind = r8) :: expow, expor, expoi, expos, expog, expoh
    real (kind = r8) :: pcaw, pcar, pcai, pcas, pcag, pcah
    real (kind = r8) :: pcbw, pcbr, pcbi, pcbs, pcbg, pcbh
    real (kind = r8) :: edaw, edar, edai, edas, edag, edah
    real (kind = r8) :: edbw, edbr, edbi, edbs, edbg, edbh
    real (kind = r8) :: oeaw, oear, oeai, oeas, oeag, oeah
    real (kind = r8) :: oebw, oebr, oebi, oebs, oebg, oebh
    real (kind = r8) :: rraw, rrar, rrai, rras, rrag, rrah
    real (kind = r8) :: rrbw, rrbr, rrbi, rrbs, rrbg, rrbh
    real (kind = r8) :: tvaw, tvar, tvai, tvas, tvag, tvah
    real (kind = r8) :: tvbw, tvbr, tvbi, tvbs, tvbg, tvbh
    
    real, allocatable :: table0 (:), table1 (:), table2 (:), table3 (:), table4 (:)
    real, allocatable :: des0 (:), des1 (:), des2 (:), des3 (:), des4 (:)
    
    ! -----------------------------------------------------------------------
    ! namelist
    ! -----------------------------------------------------------------------
    
    namelist / gfdl_mp_nml / &
        t_min, t_sub, tau_r2g, tau_smlt, tau_gmlt, dw_land, dw_ocean, vw_fac, vi_fac, &
        vr_fac, vs_fac, vg_fac, ql_mlt, do_qa, fix_negative, vw_max, vi_max, vs_max, &
        vg_max, vr_max, qs_mlt, qs0_crt, ql0_max, qi0_max, qi0_crt, ifflag, &
        rh_inc, rh_ins, rh_inr, const_vw, const_vi, const_vs, const_vg, const_vr, rthresh, &
        ccn_l, ccn_o, igflag, c_paut, tau_imlt, tau_v2l, tau_l2v, tau_i2s, &
        tau_l2r, qi_lim, ql_gen, do_hail, inflag, c_psacw, c_psaci, c_pracs, &
        c_psacr, c_pgacr, c_pgacs, c_pgacw, c_pgaci, z_slope_liq, z_slope_ice, &
        prog_ccn, c_pracw, c_praci, rad_snow, rad_graupel, rad_rain, cld_min, &
        sedflag, sed_fac, do_sedi_uv, do_sedi_w, do_sedi_heat, icloud_f, &
        irain_f, xr_a, xr_b, xr_c, ntimes, tau_revp, tice_mlt, do_cond_timescale, &
        mp_time, consv_checker, te_err, tw_err, use_rhc_cevap, use_rhc_revap, tau_wbf, &
        do_warm_rain_mp, rh_thres, f_dq_p, f_dq_m, do_cld_adj, rhc_cevap, &
        rhc_revap, beta, liq_ice_combine, rewflag, reiflag, rerflag, resflag, &
        regflag, rewmin, rewmax, reimin, reimax, rermin, rermax, resmin, &
        resmax, regmin, regmax, fs2g_fac, fi2s_fac, fi2g_fac, do_sedi_melt, &
        radr_flag, rads_flag, radg_flag, do_wbf, do_psd_water_fall, do_psd_ice_fall, &
        n0w_sig, n0i_sig, n0r_sig, n0s_sig, n0g_sig, n0h_sig, n0w_exp, n0i_exp, &
        n0r_exp, n0s_exp, n0g_exp, n0h_exp, muw, mui, mur, mus, mug, muh, &
        alinw, alini, alinr, alins, aling, alinh, blinw, blini, blinr, blins, bling, blinh, &
        do_new_acc_water, do_new_acc_ice, is_fac, ss_fac, gs_fac, rh_fac, &
        snow_grauple_combine, do_psd_water_num, do_psd_ice_num, vdiffflag, rewfac, reifac, &
        cp_heating
    
contains

! =======================================================================
! GFDL cloud microphysics initialization
! =======================================================================

subroutine gfdl_mp_init (input_nml_file, logunit, hydrostatic)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    logical, intent (in) :: hydrostatic
    
    integer, intent (in) :: logunit
    
    character (len = *), intent (in) :: input_nml_file (:)
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    logical :: exists
    
    ! -----------------------------------------------------------------------
    ! read namelist
    ! -----------------------------------------------------------------------
    
    read (input_nml_file, nml = gfdl_mp_nml)
    
    ! -----------------------------------------------------------------------
    ! write namelist to log file
    ! -----------------------------------------------------------------------
    
    write (logunit, *) " ================================================================== "
    write (logunit, *) "gfdl_mp_mod"
    write (logunit, nml = gfdl_mp_nml)
    
    ! -----------------------------------------------------------------------
    ! initialize microphysics variables
    ! -----------------------------------------------------------------------
    
    if (.not. tables_are_initialized) call qs_init
    
    call setup_mp
    
    ! -----------------------------------------------------------------------
    ! define various heat capacities and latent heat coefficients at 0 deg K
    ! -----------------------------------------------------------------------
    
    call setup_mhc_lhc (hydrostatic)
    
end subroutine gfdl_mp_init

! =======================================================================
! GFDL cloud microphysics driver
! =======================================================================

subroutine gfdl_mp_driver (qv, ql, qr, qi, qs, qg, qa, qnl, qni, pt, wa, &
        ua, va, delz, delp, gsize, dtm, hs, water, rain, ice, snow, graupel, &
        hydrostatic, is, ie, ks, ke, q_con, cappa, consv_te, adj_vmr, te, dte, &
        pcw, edw, oew, rrw, tvw, pci, edi, oei, rri, tvi, pcr, edr, oer, rrr, tvr, &
        pcs, eds, oes, rrs, tvs, pcg, edg, oeg, rrg, tvg, &
        prefluxw, prefluxr, prefluxi, prefluxs, prefluxg, condensation, &
        deposition, evaporation, sublimation, last_step, do_inline_mp)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: is, ie, ks, ke
    
    logical, intent (in) :: hydrostatic, last_step, consv_te, do_inline_mp
    
    real, intent (in) :: dtm
    
    real, intent (in), dimension (is:ie) :: hs, gsize
    
    real, intent (in), dimension (is:ie, ks:ke) :: qnl, qni
    
    real, intent (inout), dimension (is:ie, ks:ke) :: delp, delz, pt, ua, va, wa, te
    real, intent (inout), dimension (is:ie, ks:ke) :: qv, ql, qr, qi, qs, qg, qa
    real, intent (inout), dimension (is:ie, ks:ke) :: prefluxw, prefluxr, prefluxi, prefluxs, prefluxg
    
    real, intent (inout), dimension (is:, ks:) :: q_con, cappa
    
    real, intent (inout), dimension (is:ie) :: water, rain, ice, snow, graupel
    real, intent (inout), dimension (is:ie) :: condensation, deposition
    real, intent (inout), dimension (is:ie) :: evaporation, sublimation
    
    real, intent (out), dimension (is:ie, ks:ke) :: adj_vmr
    real, intent (out), dimension (is:ie, ks:ke) :: pcw, edw, oew, rrw, tvw
    real, intent (out), dimension (is:ie, ks:ke) :: pci, edi, oei, rri, tvi
    real, intent (out), dimension (is:ie, ks:ke) :: pcr, edr, oer, rrr, tvr
    real, intent (out), dimension (is:ie, ks:ke) :: pcs, eds, oes, rrs, tvs
    real, intent (out), dimension (is:ie, ks:ke) :: pcg, edg, oeg, rrg, tvg

    real (kind = r8), intent (out), dimension (is:ie) :: dte

    ! -----------------------------------------------------------------------
    ! major cloud microphysics driver
    ! -----------------------------------------------------------------------
    
    call mpdrv (hydrostatic, ua, va, wa, delp, pt, qv, ql, qr, qi, qs, qg, qa, &
        qnl, qni, delz, is, ie, ks, ke, dtm, water, rain, ice, snow, graupel, &
        gsize, hs, q_con, cappa, consv_te, adj_vmr, te, dte, pcw, edw, oew, rrw, tvw, &
        pci, edi, oei, rri, tvi, pcr, edr, oer, rrr, tvr, pcs, eds, oes, rrs, tvs, &
        pcg, edg, oeg, rrg, tvg, prefluxw, prefluxr, prefluxi, &
        prefluxs, prefluxg, condensation, deposition, evaporation, sublimation, &
        last_step, do_inline_mp, .false., .true.)
    
end subroutine gfdl_mp_driver

! =======================================================================
! GFDL cloud microphysics end
! =======================================================================

subroutine gfdl_mp_end
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! free up memory
    ! -----------------------------------------------------------------------
    
    deallocate (table0)
    deallocate (table1)
    deallocate (table2)
    deallocate (table3)
    deallocate (table4)
    deallocate (des0)
    deallocate (des1)
    deallocate (des2)
    deallocate (des3)
    deallocate (des4)
    
    tables_are_initialized = .false.
    
end subroutine gfdl_mp_end

! =======================================================================
! setup cloud microphysics parameters
! =======================================================================

subroutine setup_mp
    
    implicit none
    
    integer :: i, k
    
    real :: gcon, hcon, scm3, pisq, act (20), ace (20), occ (3), aone
    
    ! -----------------------------------------------------------------------
    ! complete freezing temperature
    ! -----------------------------------------------------------------------
    
    if (do_warm_rain_mp) then
        t_wfr = t_min
    else
        t_wfr = tice - 40.0
    endif
    
    ! -----------------------------------------------------------------------
    ! cloud water autoconversion, Hong et al. (2004)
    ! -----------------------------------------------------------------------
    
    fac_rc = (4. / 3.) * pi * rhow * rthresh ** 3
    
    aone = 2. / 9. * (3. / 4.) ** (4. / 3.) / pi ** (1. / 3.)
    cpaut = c_paut * aone * grav / visd
    
    ! -----------------------------------------------------------------------
    ! terminal velocities parameters, Lin et al. (1983)
    ! -----------------------------------------------------------------------
    
    gcon = (4. * grav * rhog / (3. * cdg * rho0)) ** 0.5
    hcon = (4. * grav * rhoh / (3. * cdh * rho0)) ** 0.5
    
    ! -----------------------------------------------------------------------
    ! part of the slope parameters
    ! -----------------------------------------------------------------------
    
    normw = pi * rhow * n0w_sig * gamma (muw + 3)
    normi = pi * rhoi * n0i_sig * gamma (mui + 3)
    normr = pi * rhor * n0r_sig * gamma (mur + 3)
    norms = pi * rhos * n0s_sig * gamma (mus + 3)
    normg = pi * rhog * n0g_sig * gamma (mug + 3)
    normh = pi * rhoh * n0h_sig * gamma (muh + 3)

    expow = exp (n0w_exp / (muw + 3) * log (10.))
    expoi = exp (n0i_exp / (mui + 3) * log (10.))
    expor = exp (n0r_exp / (mur + 3) * log (10.))
    expos = exp (n0s_exp / (mus + 3) * log (10.))
    expog = exp (n0g_exp / (mug + 3) * log (10.))
    expoh = exp (n0h_exp / (muh + 3) * log (10.))

    ! -----------------------------------------------------------------------
    ! parameters for particle concentration (pc), effective diameter (ed),
    ! optical extinction (oe), radar reflectivity factor (rr), and
    ! mass-weighted terminal velocity (tv)
    ! -----------------------------------------------------------------------
    
    pcaw = exp (3 / (muw + 3) * log (n0w_sig)) * gamma (muw) * exp (3 * n0w_exp / (muw + 3) * log (10.))
    pcai = exp (3 / (mui + 3) * log (n0i_sig)) * gamma (mui) * exp (3 * n0i_exp / (mui + 3) * log (10.))
    pcar = exp (3 / (mur + 3) * log (n0r_sig)) * gamma (mur) * exp (3 * n0r_exp / (mur + 3) * log (10.))
    pcas = exp (3 / (mus + 3) * log (n0s_sig)) * gamma (mus) * exp (3 * n0s_exp / (mus + 3) * log (10.))
    pcag = exp (3 / (mug + 3) * log (n0g_sig)) * gamma (mug) * exp (3 * n0g_exp / (mug + 3) * log (10.))
    pcah = exp (3 / (muh + 3) * log (n0h_sig)) * gamma (muh) * exp (3 * n0h_exp / (muh + 3) * log (10.))

    pcbw = exp (muw / (muw + 3) * log (pi * rhow * gamma (muw + 3)))
    pcbi = exp (mui / (mui + 3) * log (pi * rhoi * gamma (mui + 3)))
    pcbr = exp (mur / (mur + 3) * log (pi * rhor * gamma (mur + 3)))
    pcbs = exp (mus / (mus + 3) * log (pi * rhos * gamma (mus + 3)))
    pcbg = exp (mug / (mug + 3) * log (pi * rhog * gamma (mug + 3)))
    pcbh = exp (muh / (muh + 3) * log (pi * rhoh * gamma (muh + 3)))
    
    edaw = exp (- 1. / (muw + 3) * log (n0w_sig)) * (muw + 2) * exp (- n0w_exp / (muw + 3) * log (10.))
    edai = exp (- 1. / (mui + 3) * log (n0i_sig)) * (mui + 2) * exp (- n0i_exp / (mui + 3) * log (10.))
    edar = exp (- 1. / (mur + 3) * log (n0r_sig)) * (mur + 2) * exp (- n0r_exp / (mur + 3) * log (10.))
    edas = exp (- 1. / (mus + 3) * log (n0s_sig)) * (mus + 2) * exp (- n0s_exp / (mus + 3) * log (10.))
    edag = exp (- 1. / (mug + 3) * log (n0g_sig)) * (mug + 2) * exp (- n0g_exp / (mug + 3) * log (10.))
    edah = exp (- 1. / (muh + 3) * log (n0h_sig)) * (muh + 2) * exp (- n0h_exp / (muh + 3) * log (10.))

    edbw = exp (1. / (muw + 3) * log (pi * rhow * gamma (muw + 3)))
    edbi = exp (1. / (mui + 3) * log (pi * rhoi * gamma (mui + 3)))
    edbr = exp (1. / (mur + 3) * log (pi * rhor * gamma (mur + 3)))
    edbs = exp (1. / (mus + 3) * log (pi * rhos * gamma (mus + 3)))
    edbg = exp (1. / (mug + 3) * log (pi * rhog * gamma (mug + 3)))
    edbh = exp (1. / (muh + 3) * log (pi * rhoh * gamma (muh + 3)))
    
    oeaw = exp (1. / (muw + 3) * log (n0w_sig)) * pi * gamma (muw + 2) * &
        exp (n0w_exp / (muw + 3) * log (10.))
    oeai = exp (1. / (mui + 3) * log (n0i_sig)) * pi * gamma (mui + 2) * &
        exp (n0i_exp / (mui + 3) * log (10.))
    oear = exp (1. / (mur + 3) * log (n0r_sig)) * pi * gamma (mur + 2) * &
        exp (n0r_exp / (mur + 3) * log (10.))
    oeas = exp (1. / (mus + 3) * log (n0s_sig)) * pi * gamma (mus + 2) * &
        exp (n0s_exp / (mus + 3) * log (10.))
    oeag = exp (1. / (mug + 3) * log (n0g_sig)) * pi * gamma (mug + 2) * &
        exp (n0g_exp / (mug + 3) * log (10.))
    oeah = exp (1. / (muh + 3) * log (n0h_sig)) * pi * gamma (muh + 2) * &
        exp (n0h_exp / (muh + 3) * log (10.))

    oebw = 2 * exp ((muw + 2) / (muw + 3) * log (pi * rhow * gamma (muw + 3)))
    oebi = 2 * exp ((mui + 2) / (mui + 3) * log (pi * rhoi * gamma (mui + 3)))
    oebr = 2 * exp ((mur + 2) / (mur + 3) * log (pi * rhor * gamma (mur + 3)))
    oebs = 2 * exp ((mus + 2) / (mus + 3) * log (pi * rhos * gamma (mus + 3)))
    oebg = 2 * exp ((mug + 2) / (mug + 3) * log (pi * rhog * gamma (mug + 3)))
    oebh = 2 * exp ((muh + 2) / (muh + 3) * log (pi * rhoh * gamma (muh + 3)))
    
    rraw = exp (- 3 / (muw + 3) * log (n0w_sig)) * gamma (muw + 6) * &
        exp (- 3 * n0w_exp / (muw + 3) * log (10.))
    rrai = exp (- 3 / (mui + 3) * log (n0i_sig)) * gamma (mui + 6) * &
        exp (- 3 * n0i_exp / (mui + 3) * log (10.))
    rrar = exp (- 3 / (mur + 3) * log (n0r_sig)) * gamma (mur + 6) * &
        exp (- 3 * n0r_exp / (mur + 3) * log (10.))
    rras = exp (- 3 / (mus + 3) * log (n0s_sig)) * gamma (mus + 6) * &
        exp (- 3 * n0s_exp / (mus + 3) * log (10.))
    rrag = exp (- 3 / (mug + 3) * log (n0g_sig)) * gamma (mug + 6) * &
        exp (- 3 * n0g_exp / (mug + 3) * log (10.))
    rrah = exp (- 3 / (muh + 3) * log (n0h_sig)) * gamma (muh + 6) * &
        exp (- 3 * n0h_exp / (muh + 3) * log (10.))

    rrbw = exp ((muw + 6) / (muw + 3) * log (pi * rhow * gamma (muw + 3)))
    rrbi = exp ((mui + 6) / (mui + 3) * log (pi * rhoi * gamma (mui + 3)))
    rrbr = exp ((mur + 6) / (mur + 3) * log (pi * rhor * gamma (mur + 3)))
    rrbs = exp ((mus + 6) / (mus + 3) * log (pi * rhos * gamma (mus + 3)))
    rrbg = exp ((mug + 6) / (mug + 3) * log (pi * rhog * gamma (mug + 3)))
    rrbh = exp ((muh + 6) / (muh + 3) * log (pi * rhoh * gamma (muh + 3)))
    
    tvaw = exp (- blinw / (muw + 3) * log (n0w_sig)) * alinw * gamma (muw + blinw + 3) * &
        exp (- blinw * n0w_exp / (muw + 3) * log (10.))
    tvai = exp (- blini / (mui + 3) * log (n0i_sig)) * alini * gamma (mui + blini + 3) * &
        exp (- blini * n0i_exp / (mui + 3) * log (10.))
    tvar = exp (- blinr / (mur + 3) * log (n0r_sig)) * alinr * gamma (mur + blinr + 3) * &
        exp (- blinr * n0r_exp / (mur + 3) * log (10.))
    tvas = exp (- blins / (mus + 3) * log (n0s_sig)) * alins * gamma (mus + blins + 3) * &
        exp (- blins * n0s_exp / (mus + 3) * log (10.))
    tvag = exp (- bling / (mug + 3) * log (n0g_sig)) * aling * gamma (mug + bling + 3) * &
        exp (- bling * n0g_exp / (mug + 3) * log (10.)) * gcon
    tvah = exp (- blinh / (muh + 3) * log (n0h_sig)) * alinh * gamma (muh + blinh + 3) * &
        exp (- blinh * n0h_exp / (muh + 3) * log (10.)) * hcon

    tvbw = exp (blinw / (muw + 3) * log (pi * rhow * gamma (muw + 3))) * gamma (muw + 3)
    tvbi = exp (blini / (mui + 3) * log (pi * rhoi * gamma (mui + 3))) * gamma (mui + 3)
    tvbr = exp (blinr / (mur + 3) * log (pi * rhor * gamma (mur + 3))) * gamma (mur + 3)
    tvbs = exp (blins / (mus + 3) * log (pi * rhos * gamma (mus + 3))) * gamma (mus + 3)
    tvbg = exp (bling / (mug + 3) * log (pi * rhog * gamma (mug + 3))) * gamma (mug + 3)
    tvbh = exp (blinh / (muh + 3) * log (pi * rhoh * gamma (muh + 3))) * gamma (muh + 3)
    
    ! -----------------------------------------------------------------------
    ! Schmidt number, Sc ** (1 / 3) in Lin et al. (1983)
    ! -----------------------------------------------------------------------
    
    scm3 = exp (1. / 3. * log (visk / vdifu))
    
    pisq = pi * pi
    
    ! -----------------------------------------------------------------------
    ! accretion between cloud water, cloud ice, rain, snow, and graupel or hail, Lin et al. (1983)
    ! -----------------------------------------------------------------------

    cracw = pi * n0r_sig * alinr * gamma (2 + mur + blinr) / &
         (4. * exp ((2 + mur + blinr) / (mur + 3) * log (normr))) * &
         exp ((1 - blinr) * log (expor))
    craci = pi * n0r_sig * alinr * gamma (2 + mur + blinr) / &
         (4. * exp ((2 + mur + blinr) / (mur + 3) * log (normr))) * &
         exp ((1 - blinr) * log (expor))
    csacw = pi * n0s_sig * alins * gamma (2 + mus + blins) / &
         (4. * exp ((2 + mus + blins) / (mus + 3) * log (norms))) * &
         exp ((1 - blins) * log (expos))
    csaci = pi * n0s_sig * alins * gamma (2 + mus + blins) / &
         (4. * exp ((2 + mus + blins) / (mus + 3) * log (norms))) * &
         exp ((1 - blins) * log (expos))
    if (do_hail) then
        cgacw = pi * n0h_sig * alinh * gamma (2 + muh + blinh) * hcon / &
             (4. * exp ((2 + muh + blinh) / (muh + 3) * log (normh))) * &
             exp ((1 - blinh) * log (expoh))
        cgaci = pi * n0h_sig * alinh * gamma (2 + muh + blinh) * hcon / &
             (4. * exp ((2 + muh + blinh) / (muh + 3) * log (normh))) * &
             exp ((1 - blinh) * log (expoh))
    else
        cgacw = pi * n0g_sig * aling * gamma (2 + mug + bling) * gcon / &
             (4. * exp ((2 + mug + bling) / (mug + 3) * log (normg))) * &
             exp ((1 - bling) * log (expog))
        cgaci = pi * n0g_sig * aling * gamma (2 + mug + bling) * gcon / &
             (4. * exp ((2 + mug + bling) / (mug + 3) * log (normg))) * &
             exp ((1 - bling) * log (expog))
    endif

    if (do_new_acc_water) then
    
        cracw = pisq * n0r_sig * n0w_sig * rhow / 24.
        csacw = pisq * n0s_sig * n0w_sig * rhow / 24.
        if (do_hail) then
            cgacw = pisq * n0h_sig * n0w_sig * rhow / 24.
        else
            cgacw = pisq * n0g_sig * n0w_sig * rhow / 24.
        endif

    endif

    if (do_new_acc_ice) then
    
        craci = pisq * n0r_sig * n0i_sig * rhoi / 24.
        csaci = pisq * n0s_sig * n0i_sig * rhoi / 24.
        if (do_hail) then
            cgaci = pisq * n0h_sig * n0i_sig * rhoi / 24.
        else
            cgaci = pisq * n0g_sig * n0i_sig * rhoi / 24.
        endif

    endif

    cracw = cracw * c_pracw
    craci = craci * c_praci
    csacw = csacw * c_psacw
    csaci = csaci * c_psaci
    cgacw = cgacw * c_pgacw
    cgaci = cgaci * c_pgaci

    ! -----------------------------------------------------------------------
    ! accretion between cloud water, cloud ice, rain, snow, and graupel or hail, Lin et al. (1983)
    ! -----------------------------------------------------------------------
    
    cracs = pisq * n0r_sig * n0s_sig * rhos / 24.
    csacr = pisq * n0s_sig * n0r_sig * rhor / 24.
    if (do_hail) then
        cgacr = pisq * n0h_sig * n0r_sig * rhor / 24.
        cgacs = pisq * n0h_sig * n0s_sig * rhos / 24.
    else
        cgacr = pisq * n0g_sig * n0r_sig * rhor / 24.
        cgacs = pisq * n0g_sig * n0s_sig * rhos / 24.
    endif
    
    cracs = cracs * c_pracs
    csacr = csacr * c_psacr
    cgacr = cgacr * c_pgacr
    cgacs = cgacs * c_pgacs
    
    ! act / ace / acc:
    !  1 -  2: racs (s - r)
    !  3 -  4: sacr (r - s)
    !  5 -  6: gacr (r - g)
    !  7 -  8: gacs (s - g)
    !  9 - 10: racw (w - r)
    ! 11 - 12: raci (i - r)
    ! 13 - 14: sacw (w - s)
    ! 15 - 16: saci (i - s)
    ! 17 - 18: sacw (w - g)
    ! 19 - 20: saci (i - g)
    
    act (1) = norms
    act (2) = normr
    act (3) = act (2)
    act (4) = act (1)
    act (5) = act (2)
    if (do_hail) then
        act (6) = normh
    else
        act (6) = normg
    endif
    act (7) = act (1)
    act (8) = act (6)
    act (9) = normw
    act (10) = act (2)
    act (11) = normi
    act (12) = act (2)
    act (13) = act (9)
    act (14) = act (1)
    act (15) = act (11)
    act (16) = act (1)
    act (17) = act (9)
    act (18) = act (6)
    act (19) = act (11)
    act (20) = act (6)
    
    ace (1) = expos
    ace (2) = expor
    ace (3) = ace (2)
    ace (4) = ace (1)
    ace (5) = ace (2)
    if (do_hail) then
        ace (6) = expoh
    else
        ace (6) = expog
    endif
    ace (7) = ace (1)
    ace (8) = ace (6)
    ace (9) = expow
    ace (10) = ace (2)
    ace (11) = expoi
    ace (12) = ace (2)
    ace (13) = ace (9)
    ace (14) = ace (1)
    ace (15) = ace (11)
    ace (16) = ace (1)
    ace (17) = ace (9)
    ace (18) = ace (6)
    ace (19) = ace (11)
    ace (20) = ace (6)
    
    acc (1) = mus
    acc (2) = mur
    acc (3) = acc (2)
    acc (4) = acc (1)
    acc (5) = acc (2)
    if (do_hail) then
        acc (6) = muh
    else
        acc (6) = mug
    endif
    acc (7) = acc (1)
    acc (8) = acc (6)
    acc (9) = muw
    acc (10) = acc (2)
    acc (11) = mui
    acc (12) = acc (2)
    acc (13) = acc (9)
    acc (14) = acc (1)
    acc (15) = acc (11)
    acc (16) = acc (1)
    acc (17) = acc (9)
    acc (18) = acc (6)
    acc (19) = acc (11)
    acc (20) = acc (6)
    
    occ (1) = 1.
    occ (2) = 2.
    occ (3) = 1.
    
    do i = 1, 3
        do k = 1, 10
            acco (i, k) = occ (i) * gamma (6 + acc (2 * k - 1) - i) * gamma (acc (2 * k) + i - 1) / &
                 (exp ((6 + acc (2 * k - 1) - i) / (acc (2 * k - 1) + 3) * log (act (2 * k - 1))) * &
                exp ((acc (2 * k) + i - 1) / (acc (2 * k) + 3) * log (act (2 * k)))) * &
                exp ((i - 3) * log (ace (2 * k - 1))) * exp ((4 - i) * log (ace (2 * k)))
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! rain evaporation, snow sublimation, and graupel or hail sublimation, Lin et al. (1983)
    ! -----------------------------------------------------------------------
    
    crevp (1) = 2. * pi * vdifu * tcond * rvgas * n0r_sig * gamma (1 + mur) / &
        exp ((1 + mur) / (mur + 3) * log (normr)) * exp (2.0 * log (expor))
    crevp (2) = 0.78
    crevp (3) = 0.31 * scm3 * sqrt (alinr / visk) * gamma ((3 + 2 * mur + blinr) / 2) / &
        exp ((3 + 2 * mur + blinr) / (mur + 3) / 2 * log (normr)) * &
        exp ((1 + mur) / (mur + 3) * log (normr)) / gamma (1 + mur) * &
        exp ((- 1 - blinr) / 2. * log (expor))
    crevp (4) = tcond * rvgas
    crevp (5) = vdifu
    
    cssub (1) = 2. * pi * vdifu * tcond * rvgas * n0s_sig * gamma (1 + mus) / &
        exp ((1 + mus) / (mus + 3) * log (norms)) * exp (2.0 * log (expos))
    cssub (2) = 0.78
    cssub (3) = 0.31 * scm3 * sqrt (alins / visk) * gamma ((3 + 2 * mus + blins) / 2) / &
        exp ((3 + 2 * mus + blins) / (mus + 3) / 2 * log (norms)) * &
        exp ((1 + mus) / (mus + 3) * log (norms)) / gamma (1 + mus) * &
        exp ((- 1 - blins) / 2. * log (expos))
    cssub (4) = tcond * rvgas
    cssub (5) = vdifu
    
    if (do_hail) then
        cgsub (1) = 2. * pi * vdifu * tcond * rvgas * n0h_sig * gamma (1 + muh) / &
            exp ((1 + muh) / (muh + 3) * log (normh)) * exp (2.0 * log (expoh))
        cgsub (2) = 0.78
        cgsub (3) = 0.31 * scm3 * sqrt (alinh * hcon / visk) * gamma ((3 + 2 * muh + blinh) / 2) / &
            exp (1. / (muh + 3) * (3 + 2 * muh + blinh) / 2 * log (normh)) * &
            exp (1. / (muh + 3) * (1 + muh) * log (normh)) / gamma (1 + muh) * &
            exp ((- 1 - blinh) / 2. * log (expoh))
    else
        cgsub (1) = 2. * pi * vdifu * tcond * rvgas * n0g_sig * gamma (1 + mug) / &
            exp ((1 + mug) / (mug + 3) * log (normg)) * exp (2.0 * log (expog))
        cgsub (2) = 0.78
        cgsub (3) = 0.31 * scm3 * sqrt (aling * gcon / visk) * gamma ((3 + 2 * mug + bling) / 2) / &
            exp ((3 + 2 * mug + bling) / (mug + 3) / 2 * log (normg)) * &
            exp ((1 + mug) / (mug + 3) * log (normg)) / gamma (1 + mug) * &
            exp ((- 1 - bling) / 2. * log (expog))
    endif
    cgsub (4) = tcond * rvgas
    cgsub (5) = vdifu
    
    ! -----------------------------------------------------------------------
    ! snow melting, Lin et al. (1983)
    ! -----------------------------------------------------------------------
    
    csmlt (1) = 2. * pi * tcond * n0s_sig * gamma (1 + mus) / &
        exp ((1 + mus) / (mus + 3) * log (norms)) * exp (2.0 * log (expos))
    csmlt (2) = 2. * pi * vdifu * n0s_sig * gamma (1 + mus) / &
        exp ((1 + mus) / (mus + 3) * log (norms)) * exp (2.0 * log (expos))
    csmlt (3) = cssub (2)
    csmlt (4) = cssub (3)
    
    ! -----------------------------------------------------------------------
    ! graupel or hail melting, Lin et al. (1983)
    ! -----------------------------------------------------------------------
    
    if (do_hail) then
        cgmlt (1) = 2. * pi * tcond * n0h_sig * gamma (1 + muh) / &
            exp ((1 + muh) / (muh + 3) * log (normh)) * exp (2.0 * log (expoh))
        cgmlt (2) = 2. * pi * vdifu * n0h_sig * gamma (1 + muh) / &
            exp ((1 + muh) / (muh + 3) * log (normh)) * exp (2.0 * log (expoh))
    else
        cgmlt (1) = 2. * pi * tcond * n0g_sig * gamma (1 + mug) / &
            exp ((1 + mug) / (mug + 3) * log (normg)) * exp (2.0 * log (expog))
        cgmlt (2) = 2. * pi * vdifu * n0g_sig * gamma (1 + mug) / &
            exp ((1 + mug) / (mug + 3) * log (normg)) * exp (2.0 * log (expog))
    endif
    cgmlt (3) = cgsub (2)
    cgmlt (4) = cgsub (3)
    
    ! -----------------------------------------------------------------------
    ! rain freezing, Lin et al. (1983)
    ! -----------------------------------------------------------------------
    
    cgfr (1) = 1.e2 / 36 * pisq * n0r_sig * rhor * gamma (6 + mur) / &
        exp ((6 + mur) / (mur + 3) * log (normr)) * exp (- 3.0 * log (expor))
    cgfr (2) = 0.66
    
end subroutine setup_mp

! =======================================================================
! define various heat capacities and latent heat coefficients at 0 deg K
! =======================================================================

subroutine setup_mhc_lhc (hydrostatic)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    logical, intent (in) :: hydrostatic
    
    if (hydrostatic) then
        c_air = cp_air
        c_vap = cp_vap
        do_sedi_w = .false.
    else
        c_air = cv_air
        c_vap = cv_vap
    endif
    d0_vap = c_vap - c_liq
    
    ! scaled constants (to reduce float point errors for 32-bit)
    
    d1_vap = d0_vap / c_air
    d1_ice = dc_ice / c_air
    
    lv00 = (hlv - d0_vap * tice) / c_air
    li00 = (hlf - dc_ice * tice) / c_air
    li20 = lv00 + li00
    
    c1_vap = c_vap / c_air
    c1_liq = c_liq / c_air
    c1_ice = c_ice / c_air
    
end subroutine setup_mhc_lhc

! =======================================================================
! major cloud microphysics driver
! =======================================================================

subroutine mpdrv (hydrostatic, ua, va, wa, delp, pt, qv, ql, qr, qi, qs, qg, &
        qa, qnl, qni, delz, is, ie, ks, ke, dtm, water, rain, ice, snow, graupel, &
        gsize, hs, q_con, cappa, consv_te, adj_vmr, te, dte, pcw, edw, oew, rrw, tvw, &
        pci, edi, oei, rri, tvi, pcr, edr, oer, rrr, tvr, pcs, eds, oes, rrs, tvs, &
        pcg, edg, oeg, rrg, tvg, prefluxw, prefluxr, prefluxi, prefluxs, prefluxg, &
        condensation, deposition, evaporation, sublimation, last_step, do_inline_mp, &
        do_mp_fast, do_mp_full)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: is, ie, ks, ke
    
    logical, intent (in) :: hydrostatic, last_step, consv_te, do_inline_mp
    logical, intent (in) :: do_mp_fast, do_mp_full
    
    real, intent (in) :: dtm
    
    real, intent (in), dimension (is:ie) :: gsize, hs
    
    real, intent (in), dimension (is:ie, ks:ke) :: qnl, qni
    
    real, intent (inout), dimension (is:ie, ks:ke) :: delp, delz, pt, ua, va, wa
    real, intent (inout), dimension (is:ie, ks:ke) :: qv, ql, qr, qi, qs, qg, qa
    real, intent (inout), dimension (is:ie, ks:ke) :: prefluxw, prefluxr, prefluxi, prefluxs, prefluxg
    
    real, intent (inout), dimension (is:, ks:) :: q_con, cappa
    
    real, intent (inout), dimension (is:ie) :: water, rain, ice, snow, graupel
    real, intent (inout), dimension (is:ie) :: condensation, deposition
    real, intent (inout), dimension (is:ie) :: evaporation, sublimation
    
    real, intent (out), dimension (is:ie, ks:ke) :: te, adj_vmr
    real, intent (out), dimension (is:ie, ks:ke) :: pcw, edw, oew, rrw, tvw
    real, intent (out), dimension (is:ie, ks:ke) :: pci, edi, oei, rri, tvi
    real, intent (out), dimension (is:ie, ks:ke) :: pcr, edr, oer, rrr, tvr
    real, intent (out), dimension (is:ie, ks:ke) :: pcs, eds, oes, rrs, tvs
    real, intent (out), dimension (is:ie, ks:ke) :: pcg, edg, oeg, rrg, tvg
    
    real (kind = r8), intent (out), dimension (is:ie) :: dte

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: i, k, n
    
    real :: rh_adj, rh_rain, ccn0, cin0, cond, q1, q2
    real :: convt, dts, q_cond, t_lnd, t_ocn, h_var, tmp, nl, ni
    
    real, dimension (ks:ke) :: q_liq, q_sol, dp, dz, dp0
    real, dimension (ks:ke) :: qvz, qlz, qrz, qiz, qsz, qgz, qaz
    real, dimension (ks:ke) :: den, pz, denfac, ccn, cin
    real, dimension (ks:ke) :: u, v, w
    
    real (kind = r8) :: con_r8, c8, cp8
    
    real (kind = r8), dimension (is:ie, ks:ke) :: te_beg_d, te_end_d, tw_beg_d, tw_end_d
    real (kind = r8), dimension (is:ie, ks:ke) :: te_beg_m, te_end_m, tw_beg_m, tw_end_m
    
    real (kind = r8), dimension (is:ie) :: te_b_beg_d, te_b_end_d, tw_b_beg_d, tw_b_end_d, te_loss
    real (kind = r8), dimension (is:ie) :: te_b_beg_m, te_b_end_m, tw_b_beg_m, tw_b_end_m
    
    real (kind = r8), dimension (ks:ke) :: tz, tzuv, tzw
    
    ! -----------------------------------------------------------------------
    ! time steps
    ! -----------------------------------------------------------------------
    
    ntimes = max (ntimes, int (dtm / min (dtm, mp_time)))
    dts = dtm / real (ntimes)
    
    ! -----------------------------------------------------------------------
    ! initialization of total energy difference and condensation diag
    ! -----------------------------------------------------------------------
    
    dte = 0.0
    cond = 0.0
    adj_vmr = 1.0
    
    ! -----------------------------------------------------------------------
    ! unit convert to mm/day
    ! -----------------------------------------------------------------------
    
    convt = 86400. * rgrav / dts
    
    do i = is, ie
        
        ! -----------------------------------------------------------------------
        ! conversion of temperature
        ! -----------------------------------------------------------------------
        
        if (do_inline_mp) then
            do k = ks, ke
                q_cond = ql (i, k) + qr (i, k) + qi (i, k) + qs (i, k) + qg (i, k)
                tz (k) = pt (i, k) / ((1. + zvir * qv (i, k)) * (1. - q_cond))
            enddo
        else
            do k = ks, ke
                tz (k) = pt (i, k)
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! calculate base total energy
        ! -----------------------------------------------------------------------
        
        if (consv_te) then
            if (hydrostatic) then
                do k = ks, ke
                    te (i, k) = - c_air * tz (k) * delp (i, k)
                enddo
            else
                do k = ks, ke
                    te (i, k) = - mte (qv (i, k), ql (i, k), qr (i, k), qi (i, k), &
                        qs (i, k), qg (i, k), tz (k), delp (i, k), .true.) * grav
                enddo
            endif
        endif
        
        ! -----------------------------------------------------------------------
        ! total energy checker
        ! -----------------------------------------------------------------------
        
        if (consv_checker) then
            call mtetw (ks, ke, qv (i, :), ql (i, :), qr (i, :), qi (i, :), &
                qs (i, :), qg (i, :), tz, ua (i, :), va (i, :), wa (i, :), &
                delp (i, :), gsize (i), dte (i), 0.0, water (i), rain (i), &
                ice (i), snow (i), graupel (i), 0.0, 0.0, dtm, te_beg_m (i, :), &
                tw_beg_m (i, :), te_b_beg_m (i), tw_b_beg_m (i), .true., hydrostatic)
        endif
        
        do k = ks, ke
            
            ! -----------------------------------------------------------------------
            ! convert specific ratios to mass mixing ratios
            ! -----------------------------------------------------------------------
            
            qvz (k) = qv (i, k)
            qlz (k) = ql (i, k)
            qrz (k) = qr (i, k)
            qiz (k) = qi (i, k)
            qsz (k) = qs (i, k)
            qgz (k) = qg (i, k)
            qaz (k) = qa (i, k)
            
            if (do_inline_mp) then
                q_cond = qlz (k) + qrz (k) + qiz (k) + qsz (k) + qgz (k)
                con_r8 = one_r8 - (qvz (k) + q_cond)
            else
                con_r8 = one_r8 - qvz (k)
            endif            

            dp0 (k) = delp (i, k)
            dp (k) = delp (i, k) * con_r8
            con_r8 = one_r8 / con_r8
            qvz (k) = qvz (k) * con_r8
            qlz (k) = qlz (k) * con_r8
            qrz (k) = qrz (k) * con_r8
            qiz (k) = qiz (k) * con_r8
            qsz (k) = qsz (k) * con_r8
            qgz (k) = qgz (k) * con_r8
            
            ! -----------------------------------------------------------------------
            ! dry air density and layer-mean pressure thickness
            ! -----------------------------------------------------------------------
            
            dz (k) = delz (i, k)
            den (k) = - dp (k) / (grav * dz (k))
            pz (k) = den (k) * rdgas * tz (k)
            
            ! -----------------------------------------------------------------------
            ! for sedi_momentum transport
            ! -----------------------------------------------------------------------
            
            u (k) = ua (i, k)
            v (k) = va (i, k)
            if (.not. hydrostatic) then
                w (k) = wa (i, k)
            endif
            
        enddo
        
        do k = ks, ke
            denfac (k) = sqrt (den (ke) / den (k))
        enddo
        
        ! -----------------------------------------------------------------------
        ! total energy checker
        ! -----------------------------------------------------------------------
        
        if (consv_checker) then
            call mtetw (ks, ke, qvz, qlz, qrz, qiz, qsz, qgz, tz, u, v, w, &
                dp, gsize (i), dte (i), 0.0, water (i), rain (i), ice (i), &
                snow (i), graupel (i), 0.0, 0.0, dtm, te_beg_d (i, :), tw_beg_d (i, :), &
                te_b_beg_d (i), tw_b_beg_d (i), .false., hydrostatic)
        endif
        
        ! -----------------------------------------------------------------------
        ! cloud condensation nuclei (CCN), cloud ice nuclei (CIN)
        ! -----------------------------------------------------------------------
        
        if (prog_ccn) then
            do k = ks, ke
                ! boucher and lohmann (1995)
                nl = min (1., abs (hs (i)) / (10. * grav)) * &
                     (10. ** 2.24 * (qnl (i, k) * den (k) * 1.e9) ** 0.257) + &
                     (1. - min (1., abs (hs (i)) / (10. * grav))) * &
                     (10. ** 2.06 * (qnl (i, k) * den (k) * 1.e9) ** 0.48)
                ni = qni (i, k)
                ccn (k) = max (10.0, nl) * 1.e6
                cin (k) = max (10.0, ni) * 1.e6
                ccn (k) = ccn (k) / den (k)
                cin (k) = cin (k) / den (k)
            enddo
        else
            ccn0 = (ccn_l * min (1., abs (hs (i)) / (10. * grav)) + &
                ccn_o * (1. - min (1., abs (hs (i)) / (10. * grav)))) * 1.e6
            cin0 = 0.0
            do k = ks, ke
                ccn (k) = ccn0 / den (k)
                cin (k) = cin0 / den (k)
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! subgrid deviation in horizontal direction
        ! default area dependent form: use dx ~ 100 km as the base
        ! -----------------------------------------------------------------------
        
        t_lnd = dw_land * sqrt (gsize (i) / 1.e5)
        t_ocn = dw_ocean * sqrt (gsize (i) / 1.e5)
        tmp = min (1., abs (hs (i)) / (10. * grav))
        h_var = t_lnd * tmp + t_ocn * (1. - tmp)
        h_var = min (0.20, max (0.01, h_var))
        
        ! -----------------------------------------------------------------------
        ! relative humidity thresholds
        ! -----------------------------------------------------------------------
        
        rh_adj = 1. - h_var - rh_inc
        rh_rain = max (0.35, rh_adj - rh_inr)
        
        ! -----------------------------------------------------------------------
        ! fix negative water species from outside
        ! -----------------------------------------------------------------------
        
        if (fix_negative) &
            call neg_adj (ks, ke, tz, dp, qvz, qlz, qrz, qiz, qsz, qgz, cond)
        
        condensation (i) = condensation (i) + cond * convt * ntimes
        
        ! -----------------------------------------------------------------------
        ! fast microphysics loop
        ! -----------------------------------------------------------------------
        
        if (do_mp_fast) then
            
            call mp_fast (ks, ke, tz, qvz, qlz, qrz, qiz, qsz, qgz, dtm, dp, den, &
                ccn, cin, condensation (i), deposition (i), evaporation (i), &
                sublimation (i), convt)
            
        endif
        
        ! -----------------------------------------------------------------------
        ! full microphysics loop
        ! -----------------------------------------------------------------------
        
        if (do_mp_full) then
            
            call mp_full (ks, ke, ntimes, tz, qvz, qlz, qrz, qiz, qsz, qgz, dp, dz, &
                u, v, w, den, denfac, ccn, cin, dts, rh_adj, rh_rain, h_var, dte (i), &
                water (i), rain (i), ice (i), snow (i), graupel (i), prefluxw (i, :), &
                prefluxr (i, :), prefluxi (i, :), prefluxs (i, :), prefluxg (i, :), &
                condensation (i), deposition (i), evaporation (i), sublimation (i), convt)
            
        endif
        
        ! -----------------------------------------------------------------------
        ! cloud fraction diagnostic
        ! -----------------------------------------------------------------------
        
        if (do_qa .and. last_step) then
            call cloud_fraction (ks, ke, pz, den, qvz, qlz, qrz, qiz, qsz, qgz, qaz, &
                tz, h_var, gsize (i))
        endif
        
        ! =======================================================================
        ! calculation of particle concentration (pc), effective diameter (ed),
        ! optical extinction (oe), radar reflectivity factor (rr), and
        ! mass-weighted terminal velocity (tv)
        ! =======================================================================

        pcw (i, :) = 0.0
        edw (i, :) = 0.0
        oew (i, :) = 0.0
        rrw (i, :) = 0.0
        tvw (i, :) = 0.0
        pci (i, :) = 0.0
        edi (i, :) = 0.0
        oei (i, :) = 0.0
        rri (i, :) = 0.0
        tvi (i, :) = 0.0
        pcr (i, :) = 0.0
        edr (i, :) = 0.0
        oer (i, :) = 0.0
        rrr (i, :) = 0.0
        tvr (i, :) = 0.0
        pcs (i, :) = 0.0
        eds (i, :) = 0.0
        oes (i, :) = 0.0
        rrs (i, :) = 0.0
        tvs (i, :) = 0.0
        pcg (i, :) = 0.0
        edg (i, :) = 0.0
        oeg (i, :) = 0.0
        rrg (i, :) = 0.0
        tvg (i, :) = 0.0
        
        do k = ks, ke
            if (qlz (k) .gt. qcmin) then
                call cal_pc_ed_oe_rr_tv (qlz (k), den (k), blinw, muw, pcaw, pcbw, pcw (i, k), &
                    edaw, edbw, edw (i, k), oeaw, oebw, oew (i, k), rraw, rrbw, rrw (i, k), &
                    tvaw, tvbw, tvw (i, k))
            endif
            if (qiz (k) .gt. qcmin) then
                call cal_pc_ed_oe_rr_tv (qiz (k), den (k), blini, mui, pcai, pcbi, pci (i, k), &
                    edai, edbi, edi (i, k), oeai, oebi, oei (i, k), rrai, rrbi, rri (i, k), &
                    tvai, tvbi, tvi (i, k))
            endif
            if (qrz (k) .gt. qcmin) then
                call cal_pc_ed_oe_rr_tv (qrz (k), den (k), blinr, mur, pcar, pcbr, pcr (i, k), &
                    edar, edbr, edr (i, k), oear, oebr, oer (i, k), rrar, rrbr, rrr (i, k), &
                    tvar, tvbr, tvr (i, k))
            endif
            if (qsz (k) .gt. qcmin) then
                call cal_pc_ed_oe_rr_tv (qsz (k), den (k), blins, mus, pcas, pcbs, pcs (i, k), &
                    edas, edbs, eds (i, k), oeas, oebs, oes (i, k), rras, rrbs, rrs (i, k), &
                    tvas, tvbs, tvs (i, k))
            endif
            if (do_hail) then
                if (qgz (k) .gt. qcmin) then
                    call cal_pc_ed_oe_rr_tv (qgz (k), den (k), blinh, muh, pcah, pcbh, pcg (i, k), &
                        edah, edbh, edg (i, k), oeah, oebh, oeg (i, k), rrah, rrbh, rrg (i, k), &
                        tvah, tvbh, tvg (i, k))
                endif
            else
                if (qgz (k) .gt. qcmin) then
                    call cal_pc_ed_oe_rr_tv (qgz (k), den (k), bling, mug, pcag, pcbg, pcg (i, k), &
                        edag, edbg, edg (i, k), oeag, oebg, oeg (i, k), rrag, rrbg, rrg (i, k), &
                        tvag, tvbg, tvg (i, k))
                endif
            endif
        enddo

        ! -----------------------------------------------------------------------
        ! momentum transportation during sedimentation
        ! update temperature before delp and q update
        ! -----------------------------------------------------------------------
        
        if (do_sedi_uv) then
            do k = ks, ke
                c8 = mhc (qvz (k), qlz (k), qrz (k), qiz (k), qsz (k), qgz (k)) * c_air
                tzuv (k) = 0.5 * (ua (i, k) ** 2 + va (i, k) ** 2 - (u (k) ** 2 + v (k) ** 2)) / c8
                tz (k) = tz (k) + tzuv (k)
            enddo
        endif
        
        if (do_sedi_w) then
            do k = ks, ke
                c8 = mhc (qvz (k), qlz (k), qrz (k), qiz (k), qsz (k), qgz (k)) * c_air
                tzw (k) = 0.5 * (wa (i, k) ** 2 - w (k) ** 2) / c8
                tz (k) = tz (k) + tzw (k)
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! total energy checker
        ! -----------------------------------------------------------------------
        
        if (consv_checker) then
            call mtetw (ks, ke, qvz, qlz, qrz, qiz, qsz, qgz, tz, u, v, w, &
                dp, gsize (i), dte (i), 0.0, water (i), rain (i), ice (i), &
                snow (i), graupel (i), 0.0, 0.0, dtm, te_end_d (i, :), tw_end_d (i, :), &
                te_b_end_d (i), tw_b_end_d (i), .false., hydrostatic, te_loss (i))
        endif
        
        do k = ks, ke
            
            ! -----------------------------------------------------------------------
            ! convert mass mixing ratios back to specific ratios
            ! -----------------------------------------------------------------------
            
            if (do_inline_mp) then
                q_cond = qlz (k) + qrz (k) + qiz (k) + qsz (k) + qgz (k)
                con_r8 = one_r8 + qvz (k) + q_cond
            else
                con_r8 = one_r8 + qvz (k)
            endif
            
            delp (i, k) = dp (k) * con_r8
            con_r8 = one_r8 / con_r8
            qvz (k) = qvz (k) * con_r8
            qlz (k) = qlz (k) * con_r8
            qrz (k) = qrz (k) * con_r8
            qiz (k) = qiz (k) * con_r8
            qsz (k) = qsz (k) * con_r8
            qgz (k) = qgz (k) * con_r8
            
            q1 = qv (i, k) + ql (i, k) + qr (i, k) + qi (i, k) + qs (i, k) + qg (i, k)
            q2 = qvz (k) + qlz (k) + qrz (k) + qiz (k) + qsz (k) + qgz (k)
            adj_vmr (i, k) = ((one_r8 - q1) / (one_r8 - q2)) / (one_r8 + q2 - q1)

            qv (i, k) = qvz (k)
            ql (i, k) = qlz (k)
            qr (i, k) = qrz (k)
            qi (i, k) = qiz (k)
            qs (i, k) = qsz (k)
            qg (i, k) = qgz (k)
            qa (i, k) = qaz (k)
            
            ! -----------------------------------------------------------------------
            ! calculate some more variables needed outside
            ! -----------------------------------------------------------------------
            
            q_liq (k) = qlz (k) + qrz (k)
            q_sol (k) = qiz (k) + qsz (k) + qgz (k)
            q_cond = q_liq (k) + q_sol (k)
            con_r8 = one_r8 - (qvz (k) + q_cond)
            c8 = mhc (con_r8, qvz (k), q_liq (k), q_sol (k)) * c_air
            
#ifdef USE_COND
            q_con (i, k) = q_cond
#endif
#ifdef MOIST_CAPPA
            tmp = rdgas * (1. + zvir * qvz (k))
            cappa (i, k) = tmp / (tmp + c8)
#endif
            
        enddo
        
        ! -----------------------------------------------------------------------
        ! momentum transportation during sedimentation
        ! update temperature after delp and q update
        ! -----------------------------------------------------------------------
        
        if (do_sedi_uv) then
            do k = ks, ke
                tz (k) = tz (k) - tzuv (k)
                q_liq (k) = qlz (k) + qrz (k)
                q_sol (k) = qiz (k) + qsz (k) + qgz (k)
                q_cond = q_liq (k) + q_sol (k)
                con_r8 = one_r8 - (qvz (k) + q_cond)
                c8 = mhc (con_r8, qvz (k), q_liq (k), q_sol (k)) * c_air
                tzuv (k) = (0.5 * (ua (i, k) ** 2 + va (i, k) ** 2) * dp0 (k) - &
                    0.5 * (u (k) ** 2 + v (k) ** 2) * delp (i, k)) / c8 / delp (i, k)
                tz (k) = tz (k) + tzuv (k)
            enddo
            do k = ks, ke
                ua (i, k) = u (k)
                va (i, k) = v (k)
            enddo
        endif
        
        if (do_sedi_w) then
            do k = ks, ke
                tz (k) = tz (k) - tzw (k)
                q_liq (k) = qlz (k) + qrz (k)
                q_sol (k) = qiz (k) + qsz (k) + qgz (k)
                q_cond = q_liq (k) + q_sol (k)
                con_r8 = one_r8 - (qvz (k) + q_cond)
                c8 = mhc (con_r8, qvz (k), q_liq (k), q_sol (k)) * c_air
                tzw (k) = (0.5 * (wa (i, k) ** 2) * dp0 (k) - &
                    0.5 * (w (k) ** 2) * delp (i, k)) / c8 / delp (i, k)
                tz (k) = tz (k) + tzw (k)
            enddo
            do k = ks, ke
                wa (i, k) = w (k)
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! total energy checker
        ! -----------------------------------------------------------------------
        
        if (consv_checker) then
            call mtetw (ks, ke, qv (i, :), ql (i, :), qr (i, :), qi (i, :), &
                qs (i, :), qg (i, :), tz, ua (i, :), va (i, :), wa (i, :), &
                delp (i, :), gsize (i), dte (i), 0.0, water (i), rain (i), &
                ice (i), snow (i), graupel (i), 0.0, 0.0, dtm, te_end_m (i, :), &
                tw_end_m (i, :), te_b_end_m (i), tw_b_end_m (i), .true., hydrostatic)
        endif
        
        ! -----------------------------------------------------------------------
        ! calculate total energy loss or gain
        ! -----------------------------------------------------------------------
        
        if (consv_te) then
            if (hydrostatic) then
                do k = ks, ke
                    te (i, k) = te (i, k) + c_air * tz (k) * delp (i, k)
                enddo
            else
                do k = ks, ke
                    te (i, k) = te (i, k) + mte (qv (i, k), ql (i, k), qr (i, k), qi (i, k), &
                        qs (i, k), qg (i, k), tz (k), delp (i, k), .true.) * grav
                enddo
            endif
        endif
        
        ! -----------------------------------------------------------------------
        ! conversion of temperature
        ! -----------------------------------------------------------------------
        
        if (do_inline_mp) then
            do k = ks, ke
                q_cond = qlz (k) + qrz (k) + qiz (k) + qsz (k) + qgz (k)
                if (cp_heating) then
                    con_r8 = one_r8 - (qvz (k) + q_cond)
                    c8 = mhc (con_r8, qvz (k), q_liq (k), q_sol (k)) * c_air
                    cp8 = con_r8 * cp_air + qvz (k) * cp_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
                    delz (i, k) = delz (i, k) / pt (i, k)
                    pt (i, k) = pt (i, k) + (tz (k) * ((1. + zvir * qvz (k)) * (1. - q_cond)) - pt (i, k)) * c8 / cp8
                    delz (i, k) = delz (i, k) * pt (i, k)
                else
                    pt (i, k) = tz (k) * ((1. + zvir * qvz (k)) * (1. - q_cond))
                endif
            enddo
        else
            do k = ks, ke
                q_liq (k) = qlz (k) + qrz (k)
                q_sol (k) = qiz (k) + qsz (k) + qgz (k)
                q_cond = q_liq (k) + q_sol (k)
                con_r8 = one_r8 - (qvz (k) + q_cond)
                c8 = mhc (con_r8, qvz (k), q_liq (k), q_sol (k)) * c_air
                pt (i, k) = pt (i, k) + (tz (k) - pt (i, k)) * c8 / cp_air
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! total energy checker
        ! -----------------------------------------------------------------------
        
        if (consv_checker) then
            if (abs (sum (te_end_d (i, :)) + te_b_end_d (i) - sum (te_beg_d (i, :)) - te_b_beg_d (i)) / &
                 (sum (te_beg_d (i, :)) + te_b_beg_d (i)) .gt. te_err) then
                print*, "GFDL-MP-DRY TE: ", &
                    !(sum (te_beg_d (i, :)) + te_b_beg_d (i)) / (gsize (i) ** 2), &
                    !(sum (te_end_d (i, :)) + te_b_end_d (i)) / (gsize (i) ** 2), &
                    (sum (te_end_d (i, :)) + te_b_end_d (i) - sum (te_beg_d (i, :)) - te_b_beg_d (i)) / &
                    (sum (te_beg_d (i, :)) + te_b_beg_d (i))
            endif
            if (abs (sum (tw_end_d (i, :)) + tw_b_end_d (i) - sum (tw_beg_d (i, :)) - tw_b_beg_d (i)) / &
                 (sum (tw_beg_d (i, :)) + tw_b_beg_d (i)) .gt. tw_err) then
                print*, "GFDL-MP-DRY TW: ", &
                    !(sum (tw_beg_d (i, :)) + tw_b_beg_d (i)) / (gsize (i) ** 2), &
                    !(sum (tw_end_d (i, :)) + tw_b_end_d (i)) / (gsize (i) ** 2), &
                    (sum (tw_end_d (i, :)) + tw_b_end_d (i) - sum (tw_beg_d (i, :)) - tw_b_beg_d (i)) / &
                    (sum (tw_beg_d (i, :)) + tw_b_beg_d (i))
            endif
            !print*, "GFDL MP TE DRY LOSS (%) : ", te_loss (i) / (sum (te_beg_d (i, :)) + te_b_beg_d (i)) * 100.0
            if (abs (sum (te_end_m (i, :)) + te_b_end_m (i) - sum (te_beg_m (i, :)) - te_b_beg_m (i)) / &
                 (sum (te_beg_m (i, :)) + te_b_beg_m (i)) .gt. te_err) then
                print*, "GFDL-MP-WET TE: ", &
                    !(sum (te_beg_m (i, :)) + te_b_beg_m (i)) / (gsize (i) ** 2), &
                    !(sum (te_end_m (i, :)) + te_b_end_m (i)) / (gsize (i) ** 2), &
                    (sum (te_end_m (i, :)) + te_b_end_m (i) - sum (te_beg_m (i, :)) - te_b_beg_m (i)) / &
                    (sum (te_beg_m (i, :)) + te_b_beg_m (i))
            endif
            if (abs (sum (tw_end_m (i, :)) + tw_b_end_m (i) - sum (tw_beg_m (i, :)) - tw_b_beg_m (i)) / &
                 (sum (tw_beg_m (i, :)) + tw_b_beg_m (i)) .gt. tw_err) then
                print*, "GFDL-MP-WET TW: ", &
                    !(sum (tw_beg_m (i, :)) + tw_b_beg_m (i)) / (gsize (i) ** 2), &
                    !(sum (tw_end_m (i, :)) + tw_b_end_m (i)) / (gsize (i) ** 2), &
                    (sum (tw_end_m (i, :)) + tw_b_end_m (i) - sum (tw_beg_m (i, :)) - tw_b_beg_m (i)) / &
                    (sum (tw_beg_m (i, :)) + tw_b_beg_m (i))
            endif
            !print*, "GFDL MP TE WET LOSS (%) : ", te_loss_0 (i) / (sum (te_beg_m (i, :)) + te_b_beg_m (i)) * 100.0
        endif
    
    enddo ! i loop
    
end subroutine mpdrv

! =======================================================================
! fix negative water species
! =======================================================================

subroutine neg_adj (ks, ke, tz, dp, qv, ql, qr, qi, qs, qg, cond)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in), dimension (ks:ke) :: dp
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    
    real, intent (out) :: cond
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: dq, sink
    
    real, dimension (ks:ke) :: q_liq, q_sol, lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), dimension (ks:ke) :: cvm, te8
    
    ! -----------------------------------------------------------------------
    ! initialization
    ! -----------------------------------------------------------------------
    
    cond = 0
    
    ! -----------------------------------------------------------------------
    ! calculate moist heat capacity and latent heat coefficients
    ! -----------------------------------------------------------------------
    
    call cal_mhc_lhc (ks, ke, qv, ql, qr, qi, qs, qg, q_liq, q_sol, cvm, te8, tz, &
        lcpk, icpk, tcpk, tcp3)
    
    do k = ks, ke
        
        ! -----------------------------------------------------------------------
        ! fix negative solid-phase hydrometeors
        ! -----------------------------------------------------------------------
        
        ! if cloud ice < 0, borrow from snow
        if (qi (k) .lt. 0.) then
            sink = min (- qi (k), max (0., qs (k)))
            call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., 0., 0., sink, - sink, 0.)
        endif
        
        ! if snow < 0, borrow from graupel
        if (qs (k) .lt. 0.) then
            sink = min (- qs (k), max (0., qg (k)))
            call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., 0., 0., 0., sink, - sink)
        endif
        
        ! if graupel < 0, borrow from rain
        if (qg (k) .lt. 0.) then
            sink = min (- qg (k), max (0., qr (k)))
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., 0., - sink, 0., 0., sink, te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
        endif
        
        ! -----------------------------------------------------------------------
        ! fix negative liquid-phase hydrometeors
        ! -----------------------------------------------------------------------
        
        ! if rain < 0, borrow from cloud water
        if (qr (k) .lt. 0.) then
            sink = min (- qr (k), max (0., ql (k)))
            call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., - sink, sink, 0., 0., 0.)
        endif
        
        ! if cloud water < 0, borrow from water vapor
        if (ql (k) .lt. 0.) then
            sink = min (- ql (k), max (0., qv (k)))
            cond = cond + sink * dp (k)
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                 - sink, sink, 0., 0., 0., 0., te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
        endif
        
    enddo
    
    ! -----------------------------------------------------------------------
    ! fix negative water vapor
    ! -----------------------------------------------------------------------
    
    ! if water vapor < 0, borrow water vapor from below
    do k = ks, ke - 1
        if (qv (k) .lt. 0.) then
            qv (k + 1) = qv (k + 1) + qv (k) * dp (k) / dp (k + 1)
            qv (k) = 0.
        endif
    enddo
    
    ! if water vapor < 0, borrow water vapor from above
    if (qv (ke) .lt. 0. .and. qv (ke - 1) .gt. 0.) then
        dq = min (- qv (ke) * dp (ke), qv (ke - 1) * dp (ke - 1))
        qv (ke - 1) = qv (ke - 1) - dq / dp (ke - 1)
        qv (ke) = qv (ke) + dq / dp (ke)
    endif
    
end subroutine neg_adj

! =======================================================================
! full microphysics loop
! =======================================================================

subroutine mp_full (ks, ke, ntimes, tz, qv, ql, qr, qi, qs, qg, dp, dz, u, v, w, &
        den, denfac, ccn, cin, dts, rh_adj, rh_rain, h_var, dte, water, rain, ice, &
        snow, graupel, prefluxw, prefluxr, prefluxi, prefluxs, prefluxg, &
        condensation, deposition, evaporation, sublimation, convt)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke, ntimes
    
    real, intent (in) :: dts, rh_adj, rh_rain, h_var, convt
    
    real, intent (in), dimension (ks:ke) :: dp, dz, den, denfac
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, u, v, w, ccn, cin
    real, intent (inout), dimension (ks:ke) :: prefluxw, prefluxr, prefluxi, prefluxs, prefluxg
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    real, intent (inout) :: water, rain, ice, snow, graupel
    real, intent (inout) :: condensation, deposition
    real, intent (inout) :: evaporation, sublimation
    
    real (kind = r8), intent (inout) :: dte
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: n
    
    real :: w1, r1, i1, s1, g1, cond, dep, reevap, sub
    
    real, dimension (ks:ke) :: vtw, vtr, vti, vts, vtg, pfw, pfr, pfi, pfs, pfg
    
    do n = 1, ntimes
        
        ! -----------------------------------------------------------------------
        ! sedimentation of cloud ice, snow, graupel or hail, and rain
        ! -----------------------------------------------------------------------
        
        call sedimentation (dts, ks, ke, tz, qv, ql, qr, qi, qs, qg, &
            dz, dp, vtw, vtr, vti, vts, vtg, w1, r1, i1, s1, g1, pfw, pfr, pfi, pfs, pfg, &
            u, v, w, den, denfac, dte)
        
        water = water + w1 * convt
        rain = rain + r1 * convt
        ice = ice + i1 * convt
        snow = snow + s1 * convt
        graupel = graupel + g1 * convt
        
        prefluxw = prefluxw + pfw * convt
        prefluxr = prefluxr + pfr * convt
        prefluxi = prefluxi + pfi * convt
        prefluxs = prefluxs + pfs * convt
        prefluxg = prefluxg + pfg * convt
        
        ! -----------------------------------------------------------------------
        ! warm rain cloud microphysics
        ! -----------------------------------------------------------------------
        
        call warm_rain (dts, ks, ke, dp, dz, tz, qv, ql, qr, qi, qs, qg, &
            den, denfac, vtw, vtr, ccn, rh_rain, h_var, reevap)
        
        evaporation = evaporation + reevap * convt
        
        ! -----------------------------------------------------------------------
        ! ice cloud microphysics
        ! -----------------------------------------------------------------------
        
        call ice_cloud (ks, ke, tz, qv, ql, qr, qi, qs, qg, den, &
            denfac, vtw, vtr, vti, vts, vtg, dts, h_var)
        
        ! -----------------------------------------------------------------------
        ! temperature sentive high vertical resolution processes
        ! -----------------------------------------------------------------------
        
        call subgrid_z_proc (ks, ke, den, denfac, dts, rh_adj, tz, qv, ql, &
            qr, qi, qs, qg, dp, ccn, cin, cond, dep, reevap, sub)
        
        condensation = condensation + cond * convt
        deposition = deposition + dep * convt
        evaporation = evaporation + reevap * convt
        sublimation = sublimation + sub * convt
        
    enddo
    
end subroutine mp_full

! =======================================================================
! fast microphysics loop
! =======================================================================

subroutine mp_fast (ks, ke, tz, qv, ql, qr, qi, qs, qg, dtm, dp, den, &
        ccn, cin, condensation, deposition, evaporation, sublimation, convt)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dtm, convt
    
    real, intent (in), dimension (ks:ke) :: dp, den
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, ccn, cin
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    real, intent (inout) :: condensation, deposition
    real, intent (inout) :: evaporation, sublimation
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    real :: cond, dep, reevap, sub
    
    real, dimension (ks:ke) :: q_liq, q_sol, lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), dimension (ks:ke) :: cvm, te8
    
    ! -----------------------------------------------------------------------
    ! initialization
    ! -----------------------------------------------------------------------
    
    cond = 0
    dep = 0
    reevap = 0
    sub = 0
    
    ! -----------------------------------------------------------------------
    ! calculate heat capacities and latent heat coefficients
    ! -----------------------------------------------------------------------
    
    call cal_mhc_lhc (ks, ke, qv, ql, qr, qi, qs, qg, q_liq, q_sol, cvm, te8, tz, &
        lcpk, icpk, tcpk, tcp3)
    
    if (.not. do_warm_rain_mp) then
        
        ! -----------------------------------------------------------------------
        ! cloud ice melting to form cloud water and rain
        ! -----------------------------------------------------------------------
        
        call pimlt (ks, ke, dtm, qv, ql, qr, qi, qs, qg, tz, cvm, te8, &
            lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! enforce complete freezing below t_wfr
        ! -----------------------------------------------------------------------
        
        call pcomp (ks, ke, qv, ql, qr, qi, qs, qg, tz, cvm, te8, &
            lcpk, icpk, tcpk, tcp3)
        
    endif
    
    ! -----------------------------------------------------------------------
    ! cloud water condensation and evaporation
    ! -----------------------------------------------------------------------
    
    call pcond_pevap (ks, ke, dtm, qv, ql, qr, qi, qs, qg, tz, dp, cvm, te8, den, &
        lcpk, icpk, tcpk, tcp3, cond, reevap)
    
    condensation = condensation + cond * convt
    evaporation = evaporation + reevap * convt
        
    if (.not. do_warm_rain_mp) then
 
        ! -----------------------------------------------------------------------
        ! cloud water freezing to form cloud ice and snow
        ! -----------------------------------------------------------------------
        
        call pifr (ks, ke, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, &
            lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! Wegener Bergeron Findeisen process
        ! -----------------------------------------------------------------------
        
        call pwbf (ks, ke, dtm, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, &
            lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! Bigg freezing mechanism
        ! -----------------------------------------------------------------------
        
        call pbigg (ks, ke, dtm, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, ccn, &
            lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! rain freezing to form graupel
        ! -----------------------------------------------------------------------
        
        call pgfr_simp (ks, ke, dtm, qv, ql, qr, qi, qs, qg, tz, cvm, te8, &
            lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! snow melting to form cloud water and rain
        ! -----------------------------------------------------------------------
        
        call psmlt_simp (ks, ke, dtm, qv, ql, qr, qi, qs, qg, tz, cvm, te8, &
            lcpk, icpk, tcpk, tcp3)
        
    endif
    
    ! -----------------------------------------------------------------------
    ! cloud water to rain autoconversion
    ! -----------------------------------------------------------------------
    
    call praut_simp (ks, ke, dtm, tz, qv, ql, qr, qi, qs, qg)
    
    if (.not. do_warm_rain_mp) then
        
        ! -----------------------------------------------------------------------
        ! cloud ice deposition and sublimation
        ! -----------------------------------------------------------------------
        
        call pidep_pisub (ks, ke, dtm, qv, ql, qr, qi, qs, qg, tz, dp, cvm, te8, den, &
            lcpk, icpk, tcpk, tcp3, cin, dep, sub)
        
        deposition = deposition + dep * convt
        sublimation = sublimation + sub * convt
        
        ! -----------------------------------------------------------------------
        ! cloud ice to snow autoconversion
        ! -----------------------------------------------------------------------
        
        call psaut_simp (ks, ke, dtm, qv, ql, qr, qi, qs, qg, tz, den)
        
    endif
    
end subroutine mp_fast

! =======================================================================
! sedimentation of cloud ice, snow, graupel or hail, and rain
! =======================================================================

subroutine sedimentation (dts, ks, ke, tz, qv, ql, qr, qi, qs, qg, dz, dp, &
        vtw, vtr, vti, vts, vtg, w1, r1, i1, s1, g1, pfw, pfr, pfi, pfs, pfg, &
        u, v, w, den, denfac, dte)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: dp, dz, den, denfac
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, u, v, w
    
    real, intent (out) :: w1, r1, i1, s1, g1
    
    real, intent (out), dimension (ks:ke) :: vtw, vtr, vti, vts, vtg, pfw, pfr, pfi, pfs, pfg
    
    real (kind = r8), intent (inout) :: dte
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real, dimension (ks:ke) :: q_liq, q_sol, lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), dimension (ks:ke) :: te8, cvm
    
    w1 = 0.
    r1 = 0.
    i1 = 0.
    s1 = 0.
    g1 = 0.
    
    vtw = 0.
    vtr = 0.
    vti = 0.
    vts = 0.
    vtg = 0.
    
    pfw = 0.
    pfr = 0.
    pfi = 0.
    pfs = 0.
    pfg = 0.

    ! -----------------------------------------------------------------------
    ! calculate heat capacities and latent heat coefficients
    ! -----------------------------------------------------------------------
    
    call cal_mhc_lhc (ks, ke, qv, ql, qr, qi, qs, qg, q_liq, q_sol, cvm, te8, tz, &
        lcpk, icpk, tcpk, tcp3)
    
    ! -----------------------------------------------------------------------
    ! terminal fall and melting of falling cloud ice into rain
    ! -----------------------------------------------------------------------
    
    if (do_psd_ice_fall) then
        call term_rsg (ks, ke, qi, den, denfac, vi_fac, blini, mui, tvai, tvbi, vi_max, const_vi, vti)
    else
        call term_ice (ks, ke, tz, qi, den, vi_fac, vi_max, const_vi, vti)
    endif
    
    if (do_sedi_melt) then
        call sedi_melt (dts, ks, ke, tz, qv, ql, qr, qi, qs, qg, dz, dp, &
            vti, r1, tau_imlt, icpk, "qi")
    endif
    
    call terminal_fall (dts, ks, ke, tz, qv, ql, qr, qi, qs, qg, dz, dp, &
        vti, i1, pfi, u, v, w, dte, "qi")
    
    pfi (ks) = max (0.0, pfi (ks))
    do k = ke, ks + 1, -1
        pfi (k) = max (0.0, pfi (k) - pfi (k - 1))
    enddo

    ! -----------------------------------------------------------------------
    ! terminal fall and melting of falling snow into rain
    ! -----------------------------------------------------------------------
    
    call term_rsg (ks, ke, qs, den, denfac, vs_fac, blins, mus, tvas, tvbs, vs_max, const_vs, vts)
    
    if (do_sedi_melt) then
        call sedi_melt (dts, ks, ke, tz, qv, ql, qr, qi, qs, qg, dz, dp, &
            vts, r1, tau_smlt, icpk, "qs")
    endif
    
    call terminal_fall (dts, ks, ke, tz, qv, ql, qr, qi, qs, qg, dz, dp, &
        vts, s1, pfs, u, v, w, dte, "qs")
    
    pfs (ks) = max (0.0, pfs (ks))
    do k = ke, ks + 1, -1
        pfs (k) = max (0.0, pfs (k) - pfs (k - 1))
    enddo

    ! -----------------------------------------------------------------------
    ! terminal fall and melting of falling graupel into rain
    ! -----------------------------------------------------------------------
    
    if (do_hail) then
        call term_rsg (ks, ke, qg, den, denfac, vg_fac, blinh, muh, tvah, tvbh, vg_max, const_vg, vtg)
    else
        call term_rsg (ks, ke, qg, den, denfac, vg_fac, bling, mug, tvag, tvbg, vg_max, const_vg, vtg)
    endif
    
    if (do_sedi_melt) then
        call sedi_melt (dts, ks, ke, tz, qv, ql, qr, qi, qs, qg, dz, dp, &
            vtg, r1, tau_gmlt, icpk, "qg")
    endif
    
    call terminal_fall (dts, ks, ke, tz, qv, ql, qr, qi, qs, qg, dz, dp, &
        vtg, g1, pfg, u, v, w, dte, "qg")
    
    pfg (ks) = max (0.0, pfg (ks))
    do k = ke, ks + 1, -1
        pfg (k) = max (0.0, pfg (k) - pfg (k - 1))
    enddo
    
    ! -----------------------------------------------------------------------
    ! terminal fall of cloud water
    ! -----------------------------------------------------------------------
    
    if (do_psd_water_fall) then

        call term_rsg (ks, ke, ql, den, denfac, vw_fac, blinw, muw, tvaw, tvbw, vw_max, const_vw, vtw)
    
        call terminal_fall (dts, ks, ke, tz, qv, ql, qr, qi, qs, qg, dz, dp, &
            vtw, w1, pfw, u, v, w, dte, "ql")

        pfw (ks) = max (0.0, pfw (ks))
        do k = ke, ks + 1, -1
            pfw (k) = max (0.0, pfw (k) - pfw (k - 1))
        enddo

    endif
    
    ! -----------------------------------------------------------------------
    ! terminal fall of rain
    ! -----------------------------------------------------------------------
    
    call term_rsg (ks, ke, qr, den, denfac, vr_fac, blinr, mur, tvar, tvbr, vr_max, const_vr, vtr)
    
    call terminal_fall (dts, ks, ke, tz, qv, ql, qr, qi, qs, qg, dz, dp, &
        vtr, r1, pfr, u, v, w, dte, "qr")
    
    pfr (ks) = max (0.0, pfr (ks))
    do k = ke, ks + 1, -1
        pfr (k) = max (0.0, pfr (k) - pfr (k - 1))
    enddo

end subroutine sedimentation

! =======================================================================
! terminal velocity for cloud ice
! =======================================================================

subroutine term_ice (ks, ke, tz, q, den, v_fac, v_max, const_v, vt)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    logical, intent (in) :: const_v
    
    real, intent (in) :: v_fac, v_max
    
    real, intent (in), dimension (ks:ke) :: q, den
    
    real (kind = r8), intent (in), dimension (ks:ke) :: tz
    
    real, intent (out), dimension (ks:ke) :: vt
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: qden
    
    real, parameter :: aa = - 4.14122e-5
    real, parameter :: bb = - 0.00538922
    real, parameter :: cc = - 0.0516344
    real, parameter :: dd = 0.00216078
    real, parameter :: ee = 1.9714
    
    real, dimension (ks:ke) :: tc
    
    if (const_v) then
        vt (:) = v_fac
    else
        do k = ks, ke
            qden = q (k) * den (k)
            if (q (k) .lt. qfmin) then
                vt (k) = 0.0
            else
                tc (k) = tz (k) - tice
                if (ifflag .eq. 1) then
                    vt (k) = (3. + log10 (qden)) * (tc (k) * (aa * tc (k) + bb) + cc) + &
                        dd * tc (k) + ee
                    vt (k) = 0.01 * v_fac * exp (vt (k) * log (10.))
                endif
                if (ifflag .eq. 2) &
                    vt (k) = v_fac * 3.29 * exp (0.16 * log (qden))
                vt (k) = min (v_max, max (0.0, vt (k)))
            endif
        enddo
    endif
    
end subroutine term_ice

! =======================================================================
! terminal velocity for rain, snow, and graupel, Lin et al. (1983)
! =======================================================================

subroutine term_rsg (ks, ke, q, den, denfac, v_fac, blin, mu, tva, tvb, v_max, const_v, vt)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    logical, intent (in) :: const_v
    
    real, intent (in) :: v_fac, blin, v_max, mu
    
    real (kind = r8), intent (in) :: tva, tvb
    
    real, intent (in), dimension (ks:ke) :: q, den, denfac
    
    real, intent (out), dimension (ks:ke) :: vt
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    if (const_v) then
        vt (:) = v_fac
    else
        do k = ks, ke
            if (q (k) .lt. qfmin) then
                vt (k) = 0.0
            else
                call cal_pc_ed_oe_rr_tv (q (k), den (k), blin, mu, &
                    tva = tva, tvb = tvb, tv = vt (k))
                vt (k) = v_fac * vt (k) * denfac (k)
                vt (k) = min (v_max, max (0.0, vt (k)))
            endif
        enddo
    endif
    
end subroutine term_rsg

! =======================================================================
! melting during sedimentation
! =======================================================================

subroutine sedi_melt (dts, ks, ke, tz, qv, ql, qr, qi, qs, qg, dz, dp, &
        vt, r1, tau_mlt, icpk, qflag)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts, tau_mlt
    
    real, intent (in), dimension (ks:ke) :: vt, dp, dz, icpk
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    
    real, intent (inout) :: r1
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    character (len = 2), intent (in) :: qflag
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k, m
    
    real :: dtime, sink, zs
    
    real, dimension (ks:ke) :: q
    
    real, dimension (ks:ke + 1) :: ze, zt
    
    real (kind = r8), dimension (ks:ke) :: cvm
    
    call zezt (ks, ke, dts, zs, dz, vt, ze, zt)
    
    select case (qflag)
        case ("qi")
            q = qi
        case ("qs")
            q = qs
        case ("qg")
            q = qg
        case default
            print *, "gfdl_mp: qflag error!"
    end select
    
    ! -----------------------------------------------------------------------
    ! melting to rain
    ! -----------------------------------------------------------------------
    
    do k = ke - 1, ks, - 1
        if (vt (k) .lt. 1.e-10) cycle
        if (q (k) .gt. qcmin) then
            do m = k + 1, ke
                if (zt (k + 1) .ge. ze (m)) exit
                if (zt (k) .lt. ze (m + 1) .and. tz (m) .gt. tice) then
                    cvm (k) = mhc (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k))
                    cvm (m) = mhc (qv (m), ql (m), qr (m), qi (m), qs (m), qg (m))
                    dtime = min (dts, (ze (m) - ze (m + 1)) / vt (k))
                    dtime = min (1.0, dtime / tau_mlt)
                    sink = min (q (k) * dp (k) / dp (m), dtime * (tz (m) - tice) / icpk (m))
                    q (k) = q (k) - sink * dp (m) / dp (k)
                    if (zt (k) .lt. zs) then
                        r1 = r1 + sink * dp (m)
                    else
                        qr (m) = qr (m) + sink
                    endif
                    select case (qflag)
                        case ("qi")
                            qi (k) = q (k)
                        case ("qs")
                            qs (k) = q (k)
                        case ("qg")
                            qg (k) = q (k)
                        case default
                            print *, "gfdl_mp: qflag error!"
                    end select
                    tz (k) = (tz (k) * cvm (k) - li00 * sink * dp (m) / dp (k)) / &
                        mhc (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k))
                    tz (m) = (tz (m) * cvm (m)) / &
                        mhc (qv (m), ql (m), qr (m), qi (m), qs (m), qg (m))
                endif
                if (q (k) .lt. qcmin) exit
            enddo
        endif
    enddo
    
end subroutine sedi_melt

! =======================================================================
! melting during sedimentation
! =======================================================================

subroutine terminal_fall (dts, ks, ke, tz, qv, ql, qr, qi, qs, qg, dz, dp, &
        vt, x1, m1, u, v, w, dte, qflag)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: vt, dp, dz
    
    character (len = 2), intent (in) :: qflag
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, u, v, w
    
    real, intent (inout) :: x1
    
    real (kind = r8), intent (inout) :: dte
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    real, intent (out), dimension (ks:ke) :: m1
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    logical :: no_fall
    
    real :: zs
    
    real, dimension (ks:ke) :: dm, q
    
    real, dimension (ks:ke + 1) :: ze, zt
    
    real (kind = r8), dimension (ks:ke) :: te1, te2

    m1 = 0.0
    
    call zezt (ks, ke, dts, zs, dz, vt, ze, zt)
    
    select case (qflag)
        case ("ql")
            q = ql
        case ("qr")
            q = qr
        case ("qi")
            q = qi
        case ("qs")
            q = qs
        case ("qg")
            q = qg
        case default
            print *, "gfdl_mp: qflag error!"
    end select
    
    call check_column (ks, ke, q, no_fall)
    
    if (no_fall) return
    
    ! -----------------------------------------------------------------------
    ! momentum transportation during sedimentation
    ! -----------------------------------------------------------------------
    
    if (do_sedi_w) then
        do k = ks, ke
            dm (k) = dp (k) * (1. + qv (k) + ql (k) + qr (k) + qi (k) + qs (k) + qg (k))
        enddo
    endif
    
    ! -----------------------------------------------------------------------
    ! energy change during sedimentation
    ! -----------------------------------------------------------------------
    
    do k = ks, ke
        te1 (k) = mte (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), tz (k), dp (k), .false.)
    enddo
    
    ! -----------------------------------------------------------------------
    ! sedimentation
    ! -----------------------------------------------------------------------
    
    select case (qflag)
        case ("ql")
            q = ql
        case ("qr")
            q = qr
        case ("qi")
            q = qi
        case ("qs")
            q = qs
        case ("qg")
            q = qg
        case default
            print *, "gfdl_mp: qflag error!"
    end select

    if (sedflag .eq. 1) &
        call implicit_fall (dts, ks, ke, ze, vt, dp, q, x1, m1)
    if (sedflag .eq. 2) &
        call explicit_fall (dts, ks, ke, ze, vt, dp, q, x1, m1)
    if (sedflag .eq. 3) &
        call lagrangian_fall (ks, ke, zs, ze, zt, dp, q, x1, m1)
    if (sedflag .eq. 4) &
        call implicit_lagrangian_fall (dts, ks, ke, zs, ze, zt, vt, dp, q, &
            x1, m1, sed_fac)
    
    select case (qflag)
        case ("ql")
            ql = q
        case ("qr")
            qr = q
        case ("qi")
            qi = q
        case ("qs")
            qs = q
        case ("qg")
            qg = q
        case default
            print *, "gfdl_mp: qflag error!"
    end select
    
    ! -----------------------------------------------------------------------
    ! energy change during sedimentation
    ! -----------------------------------------------------------------------
    
    do k = ks, ke
        te2 (k) = mte (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), tz (k), dp (k), .false.)
    enddo
    dte = dte + sum (te1) - sum (te2)
    
    ! -----------------------------------------------------------------------
    ! momentum transportation during sedimentation
    ! -----------------------------------------------------------------------
    
    if (do_sedi_uv) then
        call sedi_uv (ks, ke, m1, dp, u, v)
    endif
    
    if (do_sedi_w) then
        call sedi_w (ks, ke, m1, w, vt, dm)
    endif
    
    ! -----------------------------------------------------------------------
    ! energy change during sedimentation heating
    ! -----------------------------------------------------------------------
    
    do k = ks, ke
        te1 (k) = mte (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), tz (k), dp (k), .false.)
    enddo
    
    ! -----------------------------------------------------------------------
    ! heat exchanges during sedimentation
    ! -----------------------------------------------------------------------
    
    if (do_sedi_heat) then
        call sedi_heat (ks, ke, dp, m1, dz, tz, qv, ql, qr, qi, qs, qg, c_ice)
    endif
    
    ! -----------------------------------------------------------------------
    ! energy change during sedimentation heating
    ! -----------------------------------------------------------------------
    
    do k = ks, ke
        te2 (k) = mte (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), tz (k), dp (k), .false.)
    enddo
    dte = dte + sum (te1) - sum (te2)
    
end subroutine terminal_fall

! =======================================================================
! calculate ze zt for sedimentation
! =======================================================================

subroutine zezt (ks, ke, dts, zs, dz, vt, ze, zt)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: dz, vt
    
    real, intent (out) :: zs
    
    real, intent (out), dimension (ks:ke + 1) :: ze, zt
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: dt5
    
    dt5 = 0.5 * dts
    zs = 0.0
    ze (ke + 1) = zs
    do k = ke, ks, - 1
        ze (k) = ze (k + 1) - dz (k)
    enddo
    zt (ks) = ze (ks)
    do k = ks + 1, ke
        zt (k) = ze (k) - dt5 * (vt (k - 1) + vt (k))
    enddo
    zt (ke + 1) = zs - dts * vt (ke)
    do k = ks, ke
        if (zt (k + 1) .ge. zt (k)) zt (k + 1) = zt (k) - dz_min
    enddo
    
end subroutine zezt

! =======================================================================
! check if water species is large enough to fall
! =======================================================================

subroutine check_column (ks, ke, q, no_fall)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: q (ks:ke)
    
    logical, intent (out) :: no_fall
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    no_fall = .true.
    
    do k = ks, ke
        if (q (k) .gt. qfmin) then
            no_fall = .false.
            exit
        endif
    enddo
    
end subroutine check_column

! =======================================================================
! warm rain cloud microphysics
! =======================================================================

subroutine warm_rain (dts, ks, ke, dp, dz, tz, qv, ql, qr, qi, qs, qg, &
        den, denfac, vtw, vtr, ccn, rh_rain, h_var, reevap)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts, rh_rain, h_var
    
    real, intent (in), dimension (ks:ke) :: dp, dz, den, denfac, vtw, vtr
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, ccn
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    real, intent (out) :: reevap
    
    ! -----------------------------------------------------------------------
    ! initialization
    ! -----------------------------------------------------------------------
    
    reevap = 0
    
    ! -----------------------------------------------------------------------
    ! rain evaporation to form water vapor
    ! -----------------------------------------------------------------------
    
    call prevp (ks, ke, dts, tz, qv, ql, qr, qi, qs, qg, den, denfac, rh_rain, h_var, dp, reevap)
    
    ! -----------------------------------------------------------------------
    ! rain accretion with cloud water
    ! -----------------------------------------------------------------------
    
    call pracw (ks, ke, dts, tz, qv, ql, qr, qi, qs, qg, den, denfac, vtw, vtr)
    
    ! -----------------------------------------------------------------------
    ! cloud water to rain autoconversion
    ! -----------------------------------------------------------------------
    
    call praut (ks, ke, dts, tz, qv, ql, qr, qi, qs, qg, den, ccn, h_var)
    
end subroutine warm_rain

! =======================================================================
! rain evaporation to form water vapor, Lin et al. (1983)
! =======================================================================

subroutine prevp (ks, ke, dts, tz, qv, ql, qr, qi, qs, qg, den, denfac, rh_rain, h_var, dp, reevap)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts, rh_rain, h_var
    
    real, intent (in), dimension (ks:ke) :: den, denfac, dp
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    real, intent (inout), dimension (ks:ke) :: qv, qr, ql, qi, qs, qg
    
    real, intent (out) :: reevap
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: dqv, qsat, dqdt, tmp, t2, qden, q_plus, q_minus, sink
    real :: qpz, dq, dqh, tin, fac_revp, rh_tem
    
    real, dimension (ks:ke) :: q_liq, q_sol, lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), dimension (ks:ke) :: cvm, te8
    
    ! -----------------------------------------------------------------------
    ! initialization
    ! -----------------------------------------------------------------------
    
    reevap = 0
    
    ! -----------------------------------------------------------------------
    ! time-scale factor
    ! -----------------------------------------------------------------------
    
    fac_revp = 1.
    if (tau_revp .gt. 1.e-6) then
        fac_revp = 1. - exp (- dts / tau_revp)
    endif
    
    ! -----------------------------------------------------------------------
    ! calculate heat capacities and latent heat coefficients
    ! -----------------------------------------------------------------------
    
    call cal_mhc_lhc (ks, ke, qv, ql, qr, qi, qs, qg, q_liq, q_sol, cvm, te8, tz, &
        lcpk, icpk, tcpk, tcp3)
    
    do k = ks, ke
        
        tin = (tz (k) * cvm (k) - lv00 * ql (k)) / mhc (qv (k) + ql (k), qr (k), q_sol (k))
        
        ! -----------------------------------------------------------------------
        ! calculate supersaturation and subgrid variability of water
        ! -----------------------------------------------------------------------
        
        qpz = qv (k) + ql (k)
        qsat = wqs (tin, den (k), dqdt)
        dqv = qsat - qv (k)
        
        dqh = max (ql (k), h_var * max (qpz, qcmin))
        dqh = min (dqh, 0.2 * qpz)
        q_minus = qpz - dqh
        q_plus = qpz + dqh
        
        ! -----------------------------------------------------------------------
        ! rain evaporation
        ! -----------------------------------------------------------------------
        
        rh_tem = qpz / qsat
            
        if (tz (k) .gt. t_wfr .and. qr (k) .gt. qcmin .and. dqv .gt. 0.0 .and. qsat .gt. q_minus) then

            if (qsat .gt. q_plus) then
                dq = qsat - qpz
            else
                dq = 0.25 * (qsat - q_minus) ** 2 / dqh
            endif
            qden = qr (k) * den (k)
            t2 = tin * tin
            sink = psub (t2, dq, qden, qsat, crevp, den (k), denfac (k), blinr, mur, lcpk (k), cvm (k))
            sink = min (qr (k), dts * fac_revp * sink, dqv / (1. + lcpk (k) * dqdt))
            if (use_rhc_revap .and. rh_tem .ge. rhc_revap) then
                sink = 0.0
            endif
            
            ! -----------------------------------------------------------------------
            ! alternative minimum evaporation in dry environmental air
            ! -----------------------------------------------------------------------
            ! tmp = min (qr (k), dim (rh_rain * qsat, qv (k)) / (1. + lcpk (k) * dqdt))
            ! sink = max (sink, tmp)
            
            reevap = reevap + sink * dp (k)
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                sink, 0., - sink, 0., 0., 0., te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo ! k loop
    
end subroutine prevp

! =======================================================================
! rain accretion with cloud water, Lin et al. (1983)
! =======================================================================

subroutine pracw (ks, ke, dts, tz, qv, ql, qr, qi, qs, qg, den, denfac, vtw, vtr)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den, denfac, vtw, vtr
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    real, intent (inout), dimension (ks:ke) :: qv, qr, ql, qi, qs, qg
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: qden, sink
    
    do k = ks, ke
        
        if (tz (k) .gt. t_wfr .and. qr (k) .gt. qcmin .and. ql (k) .gt. qcmin) then
            
            qden = qr (k) * den (k)
            if (do_new_acc_water) then
                sink = dts * acr3d (vtr (k), vtw (k), ql (k), qr (k), cracw, acco (:, 5), &
                    acc (9), acc (10), den (k))
            else
                sink = dts * acr2d (qden, cracw, denfac (k), blinr, mur)
                sink = sink / (1. + sink) * ql (k)
            endif
            
            call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., - sink, sink, 0., 0., 0.)
            
        endif
        
    enddo
    
end subroutine pracw

! =======================================================================
! cloud water to rain autoconversion, Hong et al. (2004)
! =======================================================================

subroutine praut (ks, ke, dts, tz, qv, ql, qr, qi, qs, qg, den, ccn, h_var)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts, h_var
    
    real, intent (in), dimension (ks:ke) :: den
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, ccn
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    real, parameter :: so3 = 7.0 / 3.0
    real, parameter :: so1 = - 1.0 / 3.0
    
    integer :: k
    
    real :: sink, dq, qc
    
    real, dimension (ks:ke) :: dl, c_praut

    if (irain_f .eq. 0) then
        
        call linear_prof (ke - ks + 1, ql (ks), dl (ks), z_slope_liq, h_var)
        
        do k = ks, ke
            
            if (tz (k) .gt. t_wfr .and. ql (k) .gt. qcmin) then
                
                if (do_psd_water_num) then
                    call cal_pc_ed_oe_rr_tv (ql (k), den (k), blinw, muw, &
                        pca = pcaw, pcb = pcbw, pc = ccn (k))
                    ccn (k) = ccn (k) / den (k)
                endif

                qc = fac_rc * ccn (k)
                dl (k) = min (max (qcmin, dl (k)), 0.5 * ql (k))
                dq = 0.5 * (ql (k) + dl (k) - qc)
                
                if (dq .gt. 0.) then
                    
                    c_praut (k) = cpaut * exp (so1 * log (ccn (k) * rhow))
                    sink = min (1., dq / dl (k)) * dts * c_praut (k) * den (k) * &
                        exp (so3 * log (ql (k)))
                    sink = min (ql (k), sink)
                    
                    call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                        0., - sink, sink, 0., 0., 0.)
                    
                endif
                
            endif
            
        enddo
        
    endif
    
    if (irain_f .eq. 1) then
        
        do k = ks, ke
            
            if (tz (k) .gt. t_wfr .and. ql (k) .gt. qcmin) then
                
                if (do_psd_water_num) then
                    call cal_pc_ed_oe_rr_tv (ql (k), den (k), blinw, muw, &
                        pca = pcaw, pcb = pcbw, pc = ccn (k))
                    ccn (k) = ccn (k) / den (k)
                endif

                qc = fac_rc * ccn (k)
                dq = ql (k) - qc
                
                if (dq .gt. 0.) then
                    
                    c_praut (k) = cpaut * exp (so1 * log (ccn (k) * rhow))
                    sink = min (dq, dts * c_praut (k) * den (k) * exp (so3 * log (ql (k))))
                    sink = min (ql (k), sink)
                    
                    call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                        0., - sink, sink, 0., 0., 0.)
                    
                endif
                
            endif
            
        enddo
        
    endif

end subroutine praut
    
! =======================================================================
! ice cloud microphysics
! =======================================================================

subroutine ice_cloud (ks, ke, tz, qv, ql, qr, qi, qs, qg, den, &
        denfac, vtw, vtr, vti, vts, vtg, dts, h_var)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts, h_var
    
    real, intent (in), dimension (ks:ke) :: den, denfac, vtw, vtr, vti, vts, vtg
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    real, dimension (ks:ke) :: di, q_liq, q_sol, lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), dimension (ks:ke) :: cvm, te8
    
    ! -----------------------------------------------------------------------
    ! calculate heat capacities and latent heat coefficients
    ! -----------------------------------------------------------------------
    
    call cal_mhc_lhc (ks, ke, qv, ql, qr, qi, qs, qg, q_liq, q_sol, cvm, te8, tz, &
        lcpk, icpk, tcpk, tcp3)
    
    if (.not. do_warm_rain_mp) then
        
        ! -----------------------------------------------------------------------
        ! cloud ice melting to form cloud water and rain
        ! -----------------------------------------------------------------------
        
        call pimlt (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! cloud water freezing to form cloud ice and snow
        ! -----------------------------------------------------------------------
        
        call pifr (ks, ke, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! vertical subgrid variability
        ! -----------------------------------------------------------------------
        
        call linear_prof (ke - ks + 1, qi, di, z_slope_ice, h_var)
        
        ! -----------------------------------------------------------------------
        ! snow melting (includes snow accretion with cloud water and rain) to form cloud water and rain
        ! -----------------------------------------------------------------------
        
        call psmlt (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, denfac, &
            vtw, vtr, vts, lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! graupel melting (includes graupel accretion with cloud water and rain) to form rain
        ! -----------------------------------------------------------------------
        
        call pgmlt (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, denfac, &
            vtw, vtr, vtg, lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! snow accretion with cloud ice
        ! -----------------------------------------------------------------------
        
        call psaci (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, den, denfac, vti, vts)
        
        ! -----------------------------------------------------------------------
        ! cloud ice to snow autoconversion
        ! -----------------------------------------------------------------------
        
        call psaut (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, den, di)
        
        ! -----------------------------------------------------------------------
        ! graupel accretion with cloud ice
        ! -----------------------------------------------------------------------
        
        call pgaci (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, den, denfac, vti, vtg)
        
        ! -----------------------------------------------------------------------
        ! snow accretion with rain and rain freezing to form graupel
        ! -----------------------------------------------------------------------
        
        call psacr_pgfr (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, denfac, &
            vtr, vts, lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! graupel accretion with snow
        ! -----------------------------------------------------------------------
        
        call pgacs (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, den, vts, vtg)
        
        ! -----------------------------------------------------------------------
        ! snow to graupel autoconversion
        ! -----------------------------------------------------------------------
        
        call pgaut (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, den)
        
        ! -----------------------------------------------------------------------
        ! graupel accretion with cloud water and rain
        ! -----------------------------------------------------------------------
        
        call pgacw_pgacr (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, denfac, &
            vtr, vtg, lcpk, icpk, tcpk, tcp3)
        
    endif ! do_warm_rain_mp
    
end subroutine ice_cloud

! =======================================================================
! cloud ice melting to form cloud water and rain, Lin et al. (1983)
! =======================================================================

subroutine pimlt (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, tmp, sink, fac_imlt
    
    fac_imlt = 1. - exp (- dts / tau_imlt)
    
    do k = ks, ke
        
        tc = tz (k) - tice_mlt
        
        if (tc .gt. 0 .and. qi (k) .gt. qcmin) then
            
            sink = fac_imlt * tc / icpk (k)
            sink = min (qi (k), sink)
            tmp = min (sink, dim (ql_mlt, ql (k)))
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., tmp, sink - tmp, - sink, 0., 0., te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine pimlt

! =======================================================================
! cloud water freezing to form cloud ice and snow, Lin et al. (1983)
! =======================================================================

subroutine pifr (ks, ke, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in), dimension (ks:ke) :: den
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, tmp, sink, qim
    
    do k = ks, ke
        
        tc = t_wfr - tz (k)
        
        if (tc .gt. 0. .and. ql (k) .gt. qcmin) then
            
            sink = ql (k) * tc / dt_fr
            sink = min (ql (k), sink, tc / icpk (k))
            qim = qi0_crt / den (k)
            tmp = min (sink, dim (qim, qi (k)))
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., - sink, 0., tmp, sink - tmp, 0., te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine pifr

! =======================================================================
! snow melting (includes snow accretion with cloud water and rain) to form cloud water and rain
! Lin et al. (1983)
! =======================================================================

subroutine psmlt (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, denfac, &
        vtw, vtr, vts, lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den, denfac, vtw, vtr, vts
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, factor, tmp, sink, qden, dqdt, tin, dq, qsi
    real :: psacw, psacr, pracs
    
    do k = ks, ke
        
        tc = tz (k) - tice
        
        if (tc .ge. 0. .and. qs (k) .gt. qcmin) then
            
            psacw = 0.
            qden = qs (k) * den (k)
            if (ql (k) .gt. qcmin) then
                if (do_new_acc_water) then
                    psacw = acr3d (vts (k), vtw (k), ql (k), qs (k), csacw, acco (:, 7), &
                        acc (13), acc (14), den (k))
                else
                    factor = acr2d (qden, csacw, denfac (k), blins, mus)
                    psacw = factor / (1. + dts * factor) * ql (k)
                endif
            endif
            
            psacr = 0.
            pracs = 0.
            if (qr (k) .gt. qcmin) then
                psacr = min (acr3d (vts (k), vtr (k), qr (k), qs (k), csacr, acco (:, 2), &
                    acc (3), acc (4), den (k)), qr (k) / dts)
                pracs = acr3d (vtr (k), vts (k), qs (k), qr (k), cracs, acco (:, 1), &
                    acc (1), acc (2), den (k))
            endif
            
            tin = tz (k)
            qsi = iqs (tin, den (k), dqdt)
            dq = qsi - qv (k)
            sink = max (0., pmlt (tc, dq, qden, psacw, psacr, csmlt, den (k), denfac (k), blins, mus, &
                lcpk (k), icpk (k), cvm (k)))
            
            sink = min (qs (k), (sink + pracs) * dts, tc / icpk (k))
            tmp = min (sink, dim (qs_mlt, ql (k)))
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., tmp, sink - tmp, 0., - sink, 0., te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine psmlt

! =======================================================================
! graupel melting (includes graupel accretion with cloud water and rain) to form rain
! Lin et al. (1983)
! =======================================================================

subroutine pgmlt (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, denfac, &
        vtw, vtr, vtg, lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den, denfac, vtw, vtr, vtg
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, factor, sink, qden, dqdt, tin, dq, qsi
    real :: pgacw, pgacr
    
    do k = ks, ke
        
        tc = tz (k) - tice
        
        if (tc .ge. 0. .and. qg (k) .gt. qcmin) then
            
            pgacw = 0.
            qden = qg (k) * den (k)
            if (ql (k) .gt. qcmin) then
                if (do_new_acc_water) then
                    pgacw = acr3d (vtg (k), vtw (k), ql (k), qg (k), cgacw, acco (:, 9), &
                        acc (17), acc (18), den (k))
                else
                    if (do_hail) then
                        factor = acr2d (qden, cgacw, denfac (k), blinh, muh)
                    else
                        factor = acr2d (qden, cgacw, denfac (k), bling, mug)
                    endif
                    pgacw = factor / (1. + dts * factor) * ql (k)
                endif
            endif
            
            pgacr = 0.
            if (qr (k) .gt. qcmin) then
                pgacr = min (acr3d (vtg (k), vtr (k), qr (k), qg (k), cgacr, acco (:, 3), &
                    acc (5), acc (6), den (k)), qr (k) / dts)
            endif
            
            tin = tz (k)
            qsi = iqs (tin, den (k), dqdt)
            dq = qsi - qv (k)
            if (do_hail) then
                sink = max (0., pmlt (tc, dq, qden, pgacw, pgacr, cgmlt, den (k), denfac (k), &
                    blinh, muh, lcpk (k), icpk (k), cvm (k)))
            else
                sink = max (0., pmlt (tc, dq, qden, pgacw, pgacr, cgmlt, den (k), denfac (k), &
                    bling, mug, lcpk (k), icpk (k), cvm (k)))
            endif
            
            sink = min (qg (k), sink * dts, tc / icpk (k))
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., 0., sink, 0., 0., - sink, te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine pgmlt

! =======================================================================
! snow accretion with cloud ice, Lin et al. (1983)
! =======================================================================

subroutine psaci (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, den, denfac, vti, vts)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den, denfac, vti, vts
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, factor, sink, qden
    
    do k = ks, ke
        
        tc = tz (k) - tice
        
        if (tc .lt. 0. .and. qi (k) .gt. qcmin) then
            
            sink = 0.
            qden = qs (k) * den (k)
            if (qs (k) .gt. qcmin) then
                if (do_new_acc_ice) then
                    sink = dts * acr3d (vts (k), vti (k), qi (k), qs (k), csaci, acco (:, 8), &
                        acc (15), acc (16), den (k))
                else
                    factor = dts * acr2d (qden, csaci, denfac (k), blins, mus)
                    sink = factor / (1. + factor) * qi (k)
                endif
            endif
            
            sink = min (fi2s_fac * qi (k), sink)
            
            call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., 0., 0., - sink, sink, 0.)
            
        endif
        
    enddo
    
end subroutine psaci

! =======================================================================
! cloud ice to snow autoconversion, Lin et al. (1983)
! =======================================================================

subroutine psaut (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, den, di)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, di
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, sink, fac_i2s, q_plus, qim, dq, tmp
    
    fac_i2s = 1. - exp (- dts / tau_i2s)
    
    do k = ks, ke
        
        tc = tz (k) - tice
        
        if (tc .lt. 0. .and. qi (k) .gt. qcmin) then
            
            sink = 0.
            tmp = fac_i2s * exp (0.025 * tc)
            di (k) = max (di (k), qcmin)
            q_plus = qi (k) + di (k)
            qim = qi0_crt / den (k)
            if (q_plus .gt. (qim + qcmin)) then
                if (qim .gt. (qi (k) - di (k))) then
                    dq = (0.25 * (q_plus - qim) ** 2) / di (k)
                else
                    dq = qi (k) - qim
                endif
                sink = tmp * dq
            endif
            
            sink = min (fi2s_fac * qi (k), sink)
            
            call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., 0., 0., - sink, sink, 0.)
            
        endif
        
    enddo
    
end subroutine psaut

! =======================================================================
! graupel accretion with cloud ice, Lin et al. (1983)
! =======================================================================

subroutine pgaci (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, den, denfac, vti, vtg)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den, denfac, vti, vtg
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, factor, sink, qden
    
    do k = ks, ke
        
        tc = tz (k) - tice
        
        if (tc .lt. 0. .and. qi (k) .gt. qcmin) then
            
            sink = 0.
            qden = qg (k) * den (k)
            if (qg (k) .gt. qcmin) then
                if (do_new_acc_ice) then
                    sink = dts * acr3d (vtg (k), vti (k), qi (k), qg (k), cgaci, acco (:, 10), &
                        acc (19), acc (20), den (k))
                else
                    if (do_hail) then
                        factor = dts * acr2d (qden, cgaci, denfac (k), blinh, muh)
                    else
                        factor = dts * acr2d (qden, cgaci, denfac (k), bling, mug)
                    endif
                    sink = factor / (1. + factor) * qi (k)
                endif
            endif
            
            sink = min (fi2g_fac * qi (k), sink)
            
            call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., 0., 0., - sink, 0., sink)
            
        endif
        
    enddo
    
end subroutine pgaci

! =======================================================================
! snow accretion with rain and rain freezing to form graupel, Lin et al. (1983)
! =======================================================================

subroutine psacr_pgfr (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, denfac, &
        vtr, vts, lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den, denfac, vtr, vts
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, factor, sink
    real :: psacr, pgfr
    
    do k = ks, ke
        
        tc = tz (k) - tice
        
        if (tc .lt. 0. .and. qr (k) .gt. qcmin) then
            
            psacr = 0.
            if (qs (k) .gt. qcmin) then
                psacr = dts * acr3d (vts (k), vtr (k), qr (k), qs (k), csacr, acco (:, 2), &
                    acc (3), acc (4), den (k))
            endif
            
            pgfr = dts * cgfr (1) / den (k) * (exp (- cgfr (2) * tc) - 1.) * &
                exp ((6 + mur) / (mur + 3) * log (6 * qr (k) * den (k)))
            
            sink = psacr + pgfr
            factor = min (sink, qr (k), - tc / icpk (k)) / max (sink, qcmin)
            psacr = factor * psacr
            pgfr = factor * pgfr
            
            sink = min (qr (k), psacr + pgfr)
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., 0., - sink, 0., psacr, pgfr, te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine psacr_pgfr

! =======================================================================
! graupel accretion with snow, Lin et al. (1983)
! =======================================================================

subroutine pgacs (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, den, vts, vtg)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den, vts, vtg
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: sink
    
    do k = ks, ke
        
        if (tz (k) .lt. tice .and. qs (k) .gt. qcmin .and. qg (k) .gt. qcmin) then
            
            sink = dts * acr3d (vtg (k), vts (k), qs (k), qg (k), cgacs, acco (:, 4), &
                acc (7), acc (8), den (k))
            sink = min (fs2g_fac * qs (k), sink)
            
            call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., 0., 0., 0., - sink, sink)
            
        endif
        
    enddo
    
end subroutine pgacs

! =======================================================================
! snow to graupel autoconversion, Lin et al. (1983)
! =======================================================================

subroutine pgaut (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, den)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, factor, sink, qsm
    
    do k = ks, ke
        
        tc = tz (k) - tice
        
        if (tc .lt. 0. .and. qs (k) .gt. qcmin) then
            
            sink = 0
            qsm = qs0_crt / den (k)
            if (qs (k) .gt. qsm) then
                factor = dts * 1.e-3 * exp (0.09 * (tz (k) - tice))
                sink = factor / (1. + factor) * (qs (k) - qsm)
            endif
            
            sink = min (fs2g_fac * qs (k), sink)
            
            call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., 0., 0., 0., - sink, sink)
            
        endif
        
    enddo
    
end subroutine pgaut

! =======================================================================
! graupel accretion with cloud water and rain, Lin et al. (1983)
! =======================================================================

subroutine pgacw_pgacr (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, denfac, &
        vtr, vtg, lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den, denfac, vtr, vtg
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, factor, sink, qden
    real :: pgacw, pgacr
    
    do k = ks, ke
        
        tc = tz (k) - tice
        
        if (tc .lt. 0. .and. qg (k) .gt. qcmin) then
            
            pgacw = 0.
            if (ql (k) .gt. qcmin) then
                qden = qg (k) * den (k)
                if (do_hail) then
                    factor = dts * acr2d (qden, cgacw, denfac (k), blinh, muh)
                else
                    factor = dts * acr2d (qden, cgacw, denfac (k), bling, mug)
                endif
                pgacw = factor / (1. + factor) * ql (k)
            endif
            
            pgacr = 0.
            if (qr (k) .gt. qcmin) then
                pgacr = min (dts * acr3d (vtg (k), vtr (k), qr (k), qg (k), cgacr, acco (:, 3), &
                    acc (5), acc (6), den (k)), qr (k))
            endif
            
            sink = pgacr + pgacw
            factor = min (sink, dim (tice, tz (k)) / icpk (k)) / max (sink, qcmin)
            pgacr = factor * pgacr
            pgacw = factor * pgacw
            
            sink = pgacr + pgacw
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., - pgacw, - pgacr, 0., 0., sink, te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine pgacw_pgacr

! =======================================================================
! temperature sentive high vertical resolution processes
! =======================================================================

subroutine subgrid_z_proc (ks, ke, den, denfac, dts, rh_adj, tz, qv, ql, qr, &
        qi, qs, qg, dp, ccn, cin, cond, dep, reevap, sub)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts, rh_adj
    
    real, intent (in), dimension (ks:ke) :: den, denfac, dp
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, ccn, cin
    
    real, intent (out) :: cond, dep, reevap, sub
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    real, dimension (ks:ke) :: q_liq, q_sol, q_cond, lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), dimension (ks:ke) :: cvm, te8
    
    ! -----------------------------------------------------------------------
    ! initialization
    ! -----------------------------------------------------------------------
    
    cond = 0
    dep = 0
    reevap = 0
    sub = 0
    
    ! -----------------------------------------------------------------------
    ! calculate heat capacities and latent heat coefficients
    ! -----------------------------------------------------------------------
    
    call cal_mhc_lhc (ks, ke, qv, ql, qr, qi, qs, qg, q_liq, q_sol, cvm, te8, tz, &
        lcpk, icpk, tcpk, tcp3)
    
    ! -----------------------------------------------------------------------
    ! instant processes (include deposition, evaporation, and sublimation)
    ! -----------------------------------------------------------------------
    
    if (.not. do_warm_rain_mp) then
        
        call pinst (ks, ke, qv, ql, qr, qi, qs, qg, tz, dp, cvm, te8, den, &
            lcpk, icpk, tcpk, tcp3, rh_adj, dep, sub, reevap)
        
    endif
    
    ! -----------------------------------------------------------------------
    ! cloud water condensation and evaporation
    ! -----------------------------------------------------------------------
    
    call pcond_pevap (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, dp, cvm, te8, den, &
        lcpk, icpk, tcpk, tcp3, cond, reevap)
    
    if (.not. do_warm_rain_mp) then
        
        ! -----------------------------------------------------------------------
        ! enforce complete freezing below t_wfr
        ! -----------------------------------------------------------------------
        
        call pcomp (ks, ke, qv, ql, qr, qi, qs, qg, tz, cvm, te8, lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! Wegener Bergeron Findeisen process
        ! -----------------------------------------------------------------------
        
        call pwbf (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! Bigg freezing mechanism
        ! -----------------------------------------------------------------------
        
        call pbigg (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, ccn, lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! cloud ice deposition and sublimation
        ! -----------------------------------------------------------------------
        
        call pidep_pisub (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, dp, cvm, te8, den, &
            lcpk, icpk, tcpk, tcp3, cin, dep, sub)
        
        ! -----------------------------------------------------------------------
        ! snow deposition and sublimation
        ! -----------------------------------------------------------------------
        
        call psdep_pssub (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, dp, cvm, te8, den, &
            denfac, lcpk, icpk, tcpk, tcp3, dep, sub)
        
        ! -----------------------------------------------------------------------
        ! graupel deposition and sublimation
        ! -----------------------------------------------------------------------
        
        call pgdep_pgsub (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, dp, cvm, te8, den, &
            denfac, lcpk, icpk, tcpk, tcp3, dep, sub)
        
    endif
    
end subroutine subgrid_z_proc

! =======================================================================
! instant processes (include deposition, evaporation, and sublimation)
! =======================================================================

subroutine pinst (ks, ke, qv, ql, qr, qi, qs, qg, tz, dp, cvm, te8, den, &
        lcpk, icpk, tcpk, tcp3, rh_adj, dep, sub, reevap)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: rh_adj
    
    real, intent (in), dimension (ks:ke) :: den, dp
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    real, intent (out) :: dep, reevap, sub
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: sink, tin, qpz, rh, dqdt, tmp, qsi
    
    do k = ks, ke
        
        ! -----------------------------------------------------------------------
        ! instant deposit all water vapor to cloud ice when temperature is super low
        ! -----------------------------------------------------------------------
        
        if (tz (k) .lt. t_min) then
            
            sink = dim (qv (k), qcmin)
            dep = dep + sink * dp (k)
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                 - sink, 0., 0., sink, 0., 0., te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
        ! -----------------------------------------------------------------------
        ! instant evaporation / sublimation of all clouds when rh < rh_adj
        ! -----------------------------------------------------------------------
        
        qpz = qv (k) + ql (k) + qi (k)
        tin = (te8 (k) - lv00 * qpz + li00 * (qs (k) + qg (k))) / &
            mhc (qpz, qr (k), qs (k) + qg (k))
        
        if (tin .gt. t_sub + 6.) then
            
            qsi = iqs (tin, den (k), dqdt)
            rh = qpz / qsi
            if (rh .lt. rh_adj) then
                
                sink = ql (k)
                tmp = qi (k)
                
                reevap = reevap + sink * dp (k)
                sub = sub + tmp * dp (k)
                
                call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                    sink + tmp, - sink, 0., - tmp, 0., 0., te8 (k), cvm (k), tz (k), &
                    lcpk (k), icpk (k), tcpk (k), tcp3 (k))
                
            endif
            
        endif
        
    enddo
    
end subroutine pinst

! =======================================================================
! cloud water condensation and evaporation, Hong and Lim (2006)
! =======================================================================

subroutine pcond_pevap (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, dp, cvm, te8, den, &
        lcpk, icpk, tcpk, tcp3, cond, reevap)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den, dp
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    real, intent (out) :: cond, reevap
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: sink, tin, qpz, dqdt, qsw, rh_tem, dq, factor, fac_l2v, fac_v2l
    
    fac_l2v = 1. - exp (- dts / tau_l2v)
    fac_v2l = 1. - exp (- dts / tau_v2l)
    
    do k = ks, ke
        
        tin = tz (k)
        qsw = wqs (tin, den (k), dqdt)
        qpz = qv (k) + ql (k) + qi (k)
        rh_tem = qpz / qsw
        dq = qsw - qv (k)
        if (dq .gt. 0.) then
            factor = min (1., fac_l2v * (rh_fac * dq / qsw))
            sink = min (ql (k), factor * dq / (1. + tcp3 (k) * dqdt))
            if (use_rhc_cevap .and. rh_tem .ge. rhc_cevap) then
                sink = 0.
            endif
            reevap = reevap + sink * dp (k)
        elseif (do_cond_timescale) then
            factor = min (1., fac_v2l * (rh_fac * (- dq) / qsw))
            sink = - min (qv (k), factor * (- dq) / (1. + tcp3 (k) * dqdt))
            cond = cond - sink * dp (k)
        else
            sink = - min (qv (k), - dq / (1. + tcp3 (k) * dqdt))
            cond = cond - sink * dp (k)
        endif
        
        call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
            sink, - sink, 0., 0., 0., 0., te8 (k), cvm (k), tz (k), &
            lcpk (k), icpk (k), tcpk (k), tcp3 (k))
        
    enddo
    
end subroutine pcond_pevap

! =======================================================================
! enforce complete freezing below t_wfr, Lin et al. (1983)
! =======================================================================

subroutine pcomp (ks, ke, qv, ql, qr, qi, qs, qg, tz, cvm, te8, lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, sink
    
    do k = ks, ke
        
        tc = t_wfr - tz (k)
        
        if (tc .gt. 0. .and. ql (k) .gt. qcmin) then
            
            sink = ql (k) * tc / dt_fr
            sink = min (ql (k), sink, tc / icpk (k))
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., - sink, 0., sink, 0., 0., te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine pcomp

! =======================================================================
! Wegener Bergeron Findeisen process, Storelvmo and Tan (2015)
! =======================================================================

subroutine pwbf (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, tin, sink, dqdt, qsw, qsi, qim, tmp, fac_wbf

    if (.not. do_wbf) return
    
    fac_wbf = 1. - exp (- dts / tau_wbf)
    
    do k = ks, ke
        
        tc = tice - tz (k)
        
        tin = tz (k)
        qsw = wqs (tin, den (k), dqdt)
        qsi = iqs (tin, den (k), dqdt)

        if (tc .gt. 0. .and. ql (k) .gt. qcmin .and. qi (k) .gt. qcmin .and. &
            qv (k) .gt. qsi .and. qv (k) .lt. qsw) then

            sink = min (fac_wbf * ql (k), tc / icpk (k))
            qim = qi0_crt / den (k)
            tmp = min (sink, dim (qim, qi (k)))
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., - sink, 0., tmp, sink - tmp, 0., te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine pwbf

! =======================================================================
! Bigg freezing mechanism, Bigg (1953)
! =======================================================================

subroutine pbigg (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, ccn, lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, ccn
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: sink, tc
    
    do k = ks, ke
        
        tc = tice - tz (k)
        
        if (tc .gt. 0 .and. ql (k) .gt. qcmin) then
            
            if (do_psd_water_num) then
                call cal_pc_ed_oe_rr_tv (ql (k), den (k), blinw, muw, &
                    pca = pcaw, pcb = pcbw, pc = ccn (k))
                ccn (k) = ccn (k) / den (k)
            endif

            sink = 100. / (rhow * ccn (k)) * dts * (exp (0.66 * tc) - 1.) * ql (k) ** 2
            sink = min (ql (k), sink, tc / icpk (k))
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., - sink, 0., sink, 0., 0., te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo

end subroutine pbigg
        
! =======================================================================
! cloud ice deposition and sublimation, Hong et al. (2004)
! =======================================================================

subroutine pidep_pisub (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, dp, cvm, te8, den, &
        lcpk, icpk, tcpk, tcp3, cin, dep, sub)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den, dp
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, cin
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    real, intent (out) :: dep, sub
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: sink, tin, dqdt, qsi, dq, pidep, tmp, tc, qi_gen, qi_crt
    
    do k = ks, ke
        
        if (tz (k) .lt. tice) then
            
            pidep = 0.
            tin = tz (k)
            qsi = iqs (tin, den (k), dqdt)
            dq = qv (k) - qsi
            tmp = dq / (1. + tcpk (k) * dqdt)
            
            if (qi (k) .gt. qcmin) then
                if (.not. prog_ccn) then
                    if (inflag .eq. 1) &
                        cin (k) = 5.38e7 * exp (0.75 * log (qi (k) * den (k)))
                    if (inflag .eq. 2) &
                        cin (k) = exp (- 2.80 + 0.262 * (tice - tz (k))) * 1000.0
                    if (inflag .eq. 3) &
                        cin (k) = exp (- 0.639 + 12.96 * (qv (k) / qsi - 1.0)) * 1000.0
                    if (inflag .eq. 4) &
                        cin (k) = 5.e-3 * exp (0.304 * (tice - tz (k))) * 1000.0
                    if (inflag .eq. 5) &
                        cin (k) = 1.e-5 * exp (0.5 * (tice - tz (k))) * 1000.0
                endif
                if (do_psd_ice_num) then
                    call cal_pc_ed_oe_rr_tv (qi (k), den (k), blini, mui, &
                        pca = pcai, pcb = pcbi, pc = cin (k))
                    cin (k) = cin (k) / den (k)
                endif
                pidep = dts * dq * 4.0 * 11.9 * exp (0.5 * log (qi (k) * den (k) * cin (k))) / &
                     (qsi * den (k) * (tcpk (k) * cvm (k)) ** 2 / (tcond * rvgas * tz (k) ** 2) + &
                    1. / vdifu)
            endif
            
            if (dq .gt. 0.) then
                tc = tice - tz (k)
                qi_gen = 4.92e-11 * exp (1.33 * log (1.e3 * exp (0.1 * tc)))
                if (igflag .eq. 1) &
                    qi_crt = qi_gen / den (k)
                if (igflag .eq. 2) &
                    qi_crt = qi_gen * min (qi_lim, 0.1 * tc) / den (k)
                if (igflag .eq. 3) &
                    qi_crt = 1.82e-6 * min (qi_lim, 0.1 * tc) / den (k)
                if (igflag .eq. 4) &
                    qi_crt = max (qi_gen, 1.82e-6) * min (qi_lim, 0.1 * tc) / den (k)
                sink = min (tmp, max (qi_crt - qi (k), pidep), tc / tcpk (k))
                dep = dep + sink * dp (k)
            else
                pidep = pidep * min (1., dim (tz (k), t_sub) * is_fac)
                sink = max (pidep, tmp, - qi (k))
                sub = sub - sink * dp (k)
            endif
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                 - sink, 0., 0., sink, 0., 0., te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine pidep_pisub

! =======================================================================
! snow deposition and sublimation, Lin et al. (1983)
! =======================================================================

subroutine psdep_pssub (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, dp, cvm, te8, den, &
        denfac, lcpk, icpk, tcpk, tcp3, dep, sub)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den, dp, denfac
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    real, intent (out) :: dep, sub
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: sink, tin, dqdt, qsi, qden, t2, dq, pssub
    
    do k = ks, ke
        
        if (qs (k) .gt. qcmin) then
            
            tin = tz (k)
            qsi = iqs (tin, den (k), dqdt)
            qden = qs (k) * den (k)
            t2 = tz (k) * tz (k)
            dq = qsi - qv (k)
            pssub = psub (t2, dq, qden, qsi, cssub, den (k), denfac (k), blins, mus, tcpk (k), cvm (k))
            pssub = dts * pssub
            dq = dq / (1. + tcpk (k) * dqdt)
            if (pssub .gt. 0.) then
                sink = min (pssub * min (1., dim (tz (k), t_sub) * ss_fac), qs (k))
                sub = sub + sink * dp (k)
            else
                sink = 0.
                if (tz (k) .le. tice) then
                    sink = max (pssub, dq, (tz (k) - tice) / tcpk (k))
                endif
                dep = dep - sink * dp (k)
            endif
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                sink, 0., 0., 0., - sink, 0., te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine psdep_pssub

! =======================================================================
! graupel deposition and sublimation, Lin et al. (1983)
! =======================================================================

subroutine pgdep_pgsub (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, dp, cvm, te8, den, &
        denfac, lcpk, icpk, tcpk, tcp3, dep, sub)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den, dp, denfac
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    real, intent (out) :: dep, sub
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: sink, tin, dqdt, qsi, qden, t2, dq, pgsub
    
    do k = ks, ke
        
        if (qg (k) .gt. qcmin) then
            
            tin = tz (k)
            qsi = iqs (tin, den (k), dqdt)
            qden = qg (k) * den (k)
            t2 = tz (k) * tz (k)
            dq = qsi - qv (k)
            if (do_hail) then
                pgsub = psub (t2, dq, qden, qsi, cgsub, den (k), denfac (k), &
                    blinh, muh, tcpk (k), cvm (k))
            else
                pgsub = psub (t2, dq, qden, qsi, cgsub, den (k), denfac (k), &
                    bling, mug, tcpk (k), cvm (k))
            endif
            pgsub = dts * pgsub
            dq = dq / (1. + tcpk (k) * dqdt)
            if (pgsub .gt. 0.) then
                sink = min (pgsub * min (1., dim (tz (k), t_sub) * gs_fac), qg (k))
                sub = sub + sink * dp (k)
            else
                sink = 0.
                if (tz (k) .le. tice) then
                    sink = max (pgsub, dq, (tz (k) - tice) / tcpk (k))
                endif
                dep = dep - sink * dp (k)
            endif
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                sink, 0., 0., 0., 0., - sink, te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine pgdep_pgsub

! =======================================================================
! cloud fraction diagnostic
! =======================================================================

subroutine cloud_fraction (ks, ke, pz, den, qv, ql, qr, qi, qs, qg, qa, tz, h_var, gsize)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: h_var, gsize
    
    real, intent (in), dimension (ks:ke) :: pz, den
    
    real (kind = r8), intent (in), dimension (ks:ke) :: tz
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, qa
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: q_plus, q_minus
    real :: rh, rqi, tin, qsw, qsi, qpz, qstar, sigma, gam
    real :: dqdt, dq, liq, ice
    real :: qa10, qa100
    
    real, dimension (ks:ke) :: q_liq, q_sol, q_cond, lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), dimension (ks:ke) :: cvm, te8
    
    ! -----------------------------------------------------------------------
    ! calculate heat capacities and latent heat coefficients
    ! -----------------------------------------------------------------------
    
    call cal_mhc_lhc (ks, ke, qv, ql, qr, qi, qs, qg, q_liq, q_sol, cvm, te8, tz, &
        lcpk, icpk, tcpk, tcp3)
    
    do k = ks, ke
        
        ! combine water species
        
        ice = q_sol (k)
        q_sol (k) = qi (k)
        if (rad_snow) then
            q_sol (k) = qi (k) + qs (k)
            if (rad_graupel) then
                q_sol (k) = qi (k) + qs (k) + qg (k)
            endif
        endif
        
        liq = q_liq (k)
        q_liq (k) = ql (k)
        if (rad_rain) then
            q_liq (k) = ql (k) + qr (k)
        endif
        
        q_cond (k) = q_liq (k) + q_sol (k)
        qpz = qv (k) + q_cond (k)
        
        ! use the "liquid - frozen water temperature" (tin) to compute saturated specific humidity
        
        ice = ice - q_sol (k)
        liq = liq - q_liq (k)
        tin = (te8 (k) - lv00 * qpz + li00 * ice) / mhc (qpz, liq, ice)
        
        ! calculate saturated specific humidity
        
        if (tin .le. t_wfr) then
            qstar = iqs (tin, den (k), dqdt)
        elseif (tin .ge. tice) then
            qstar = wqs (tin, den (k), dqdt)
        else
            qsi = iqs (tin, den (k), dqdt)
            qsw = wqs (tin, den (k), dqdt)
            if (q_cond (k) .gt. qcmin) then
                rqi = q_sol (k) / q_cond (k)
            else
                rqi = (tice - tin) / (tice - t_wfr)
            endif
            qstar = rqi * qsi + (1. - rqi) * qsw
        endif
        
        ! cloud schemes
        
        rh = qpz / qstar
        
        if (cfflag .eq. 1) then
            if (rh .gt. rh_thres .and. qpz .gt. qcmin) then
                
                dq = h_var * qpz
                if (do_cld_adj) then
                    q_plus = qpz + dq * f_dq_p * min (1.0, max (0.0, (pz (k) - 200.e2) / &
                         (1000.e2 - 200.e2)))
                else
                    q_plus = qpz + dq * f_dq_p
                endif
                q_minus = qpz - dq * f_dq_m
                
                if (icloud_f .eq. 2) then
                    if (qstar .lt. qpz) then
                        qa (k) = 1.
                    else
                        qa (k) = 0.
                    endif
                elseif (icloud_f .eq. 3) then
                    if (qstar .lt. qpz) then
                        qa (k) = 1.
                    else
                        if (qstar .lt. q_plus) then
                            qa (k) = (q_plus - qstar) / (dq * f_dq_p)
                        else
                            qa (k) = 0.
                        endif
                        if (q_cond (k) .gt. qcmin) then
                            qa (k) = max (cld_min, qa (k))
                        endif
                        qa (k) = min (1., qa (k))
                    endif
                else
                    if (qstar .lt. q_minus) then
                        qa (k) = 1.
                    else
                        if (qstar .lt. q_plus) then
                            if (icloud_f .eq. 0) then
                                qa (k) = (q_plus - qstar) / (dq * f_dq_p + dq * f_dq_m)
                            else
                                qa (k) = (q_plus - qstar) / ((dq * f_dq_p + dq * f_dq_m) * &
                                     (1. - q_cond (k)))
                            endif
                        else
                            qa (k) = 0.
                        endif
                        if (q_cond (k) .gt. qcmin) then
                            qa (k) = max (cld_min, qa (k))
                        endif
                        qa (k) = min (1., qa (k))
                    endif
                endif
            else
                qa (k) = 0.
            endif
        endif
        
        if (cfflag .eq. 2) then
            if (rh .ge. 1.0) then
                qa (k) = 1.0
            elseif (rh .gt. rh_thres .and. q_cond (k) .gt. qcmin) then
                qa (k) = exp (xr_a * log (rh)) * (1.0 - exp (- xr_b * max (0.0, q_cond (k)) / &
                    max (1.e-5, exp (xr_c * log (max (1.e-10, 1.0 - rh) * qstar)))))
                qa (k) = max (0.0, min (1., qa (k)))
            else
                qa (k) = 0.0
            endif
        endif
        
        if (cfflag .eq. 3) then
            if (q_cond (k) .gt. qcmin) then
                qa (k) = 1. / 50. * (5.77 * (100. - gsize / 1000.) * &
                    exp (1.07 * log (max (qcmin * 1000., q_cond (k) * 1000.))) + &
                    4.82 * (gsize / 1000. - 50.) * &
                    exp (0.94 * log (max (qcmin * 1000., q_cond (k) * 1000.))))
                qa (k) = qa (k) * (0.92 / 0.96 * q_liq (k) / q_cond (k) + &
                    1.0 / 0.96 * q_sol (k) / q_cond (k))
                qa (k) = max (0.0, min (1., qa (k)))
            else
                qa (k) = 0.0
            endif
        endif
        
        if (cfflag .eq. 4) then
            sigma = 0.28 + exp (0.49 * log (max (qcmin * 1000., q_cond (k) * 1000.)))
            gam = max (0.0, q_cond (k) * 1000.) / sigma
            if (gam .lt. 0.18) then
                qa10 = 0.
            elseif (gam .gt. 2.0) then
                qa10 = 1.0
            else
                qa10 = - 0.1754 + 0.9811 * gam - 0.2223 * gam ** 2 + 0.0104 * gam ** 3
                qa10 = max (0.0, min (1., qa10))
            endif
            if (gam .lt. 0.12) then
                qa100 = 0.
            elseif (gam .gt. 1.85) then
                qa100 = 1.0
            else
                qa100 = - 0.0913 + 0.7213 * gam + 0.1060 * gam ** 2 - 0.0946 * gam ** 3
                qa100 = max (0.0, min (1., qa100))
            endif
            qa (k) = qa10 + (log10 (gsize / 1000.) - 1) * (qa100 - qa10)
            qa (k) = max (0.0, min (1., qa (k)))
        endif
        
    enddo
    
end subroutine cloud_fraction

! =======================================================================
! piecewise parabolic lagrangian scheme
! this subroutine is the same as map1_q2 in fv_mapz_mod.
! =======================================================================

subroutine lagrangian_fall (ks, ke, zs, ze, zt, dp, q, precip, m1)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: zs
    
    real, intent (in), dimension (ks:ke + 1) :: ze, zt
    
    real, intent (in), dimension (ks:ke) :: dp
    
    real, intent (inout), dimension (ks:ke) :: q
    
    real, intent (inout) :: precip
    
    real, intent (out), dimension (ks:ke) :: m1
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k, k0, n, m
    
    real :: a4 (4, ks:ke), pl, pr, delz, esl
    
    real, parameter :: r3 = 1. / 3., r23 = 2. / 3.
    
    real, dimension (ks:ke) :: qm, dz
    
    ! -----------------------------------------------------------------------
    ! density:
    ! -----------------------------------------------------------------------
    
    do k = ks, ke
        dz (k) = zt (k) - zt (k + 1)
        q (k) = q (k) * dp (k)
        a4 (1, k) = q (k) / dz (k)
        qm (k) = 0.
    enddo
    
    ! -----------------------------------------------------------------------
    ! construct vertical profile with zt as coordinate
    ! -----------------------------------------------------------------------
    
    call cs_profile (a4 (1, ks), dz (ks), ke - ks + 1)
    
    k0 = ks
    do k = ks, ke
        do n = k0, ke
            if (ze (k) .le. zt (n) .and. ze (k) .ge. zt (n + 1)) then
                pl = (zt (n) - ze (k)) / dz (n)
                if (zt (n + 1) .le. ze (k + 1)) then
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
                    if (n .lt. ke) then
                        do m = n + 1, ke
                            ! locate the bottom edge: ze (k + 1)
                            if (ze (k + 1) .lt. zt (m + 1)) then
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
    precip = precip + m1 (ke)
    
    ! -----------------------------------------------------------------------
    ! convert back to * dry * mixing ratio:
    ! dp must be dry air_mass (because moist air mass will be changed due to terminal fall) .
    ! -----------------------------------------------------------------------
    
    do k = ks, ke
        q (k) = qm (k) / dp (k)
    enddo
    
end subroutine lagrangian_fall

! =======================================================================
! vertical profile reconstruction
! this subroutine is the same as cs_profile in fv_mapz_mod where iv = 0 and kord = 9
! =======================================================================

subroutine cs_profile (a4, del, km)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: km
    
    real, intent (in) :: del (km)
    
    real, intent (inout) :: a4 (4, km)
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    logical :: extm (km)
    
    real :: gam (km), q (km + 1), d4, bet, a_bot, grat, pmp, lac
    real :: pmp_1, lac_1, pmp_2, lac_2, da1, da2, a6da
    
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
    ! top:
    ! -----------------------------------------------------------------------
    
    q (1) = max (q (1), 0.)
    q (2) = min (q (2), max (a4 (1, 1), a4 (1, 2)))
    q (2) = max (q (2), min (a4 (1, 1), a4 (1, 2)), 0.)
    
    ! -----------------------------------------------------------------------
    ! interior:
    ! -----------------------------------------------------------------------
    
    do k = 3, km - 1
        if (gam (k - 1) * gam (k + 1) .gt. 0.) then
            ! apply large - scale constraints to all fields if not local max / min
            q (k) = min (q (k), max (a4 (1, k - 1), a4 (1, k)))
            q (k) = max (q (k), min (a4 (1, k - 1), a4 (1, k)))
        else
            if (gam (k - 1) .gt. 0.) then
                ! there exists a local max
                q (k) = max (q (k), min (a4 (1, k - 1), a4 (1, k)))
            else
                ! there exists a local min
                q (k) = min (q (k), max (a4 (1, k - 1), a4 (1, k)))
                ! positive-definite
                q (k) = max (q (k), 0.0)
            endif
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! bottom:
    ! -----------------------------------------------------------------------
    
    q (km) = min (q (km), max (a4 (1, km - 1), a4 (1, km)))
    q (km) = max (q (km), min (a4 (1, km - 1), a4 (1, km)), 0.)
    q (km + 1) = max (q (km + 1), 0.)
    
    do k = 1, km
        a4 (2, k) = q (k)
        a4 (3, k) = q (k + 1)
    enddo
    
    do k = 1, km
        if (k .eq. 1 .or. k .eq. km) then
            extm (k) = (a4 (2, k) - a4 (1, k)) * (a4 (3, k) - a4 (1, k)) .gt. 0.
        else
            extm (k) = gam (k) * gam (k + 1) .lt. 0.
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! apply constraints
    ! f (s) = al + s * [ (ar - al) + a6 * (1 - s) ] (0 <= s <= 1)
    ! always use monotonic mapping
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! top:
    ! -----------------------------------------------------------------------

    a4 (2, 1) = max (0., a4 (2, 1))
    
    ! -----------------------------------------------------------------------
    ! Huynh's 2nd constraint for interior:
    ! -----------------------------------------------------------------------

    do k = 3, km - 2
        if (extm (k)) then
            ! positive definite constraint only if true local extrema
            if (a4 (1, k) .lt. qcmin .or. extm (k - 1) .or. extm (k + 1)) then
                a4 (2, k) = a4 (1, k)
                a4 (3, k) = a4 (1, k)
            endif
        else
            a4 (4, k) = 6. * a4 (1, k) - 3. * (a4 (2, k) + a4 (3, k))
            if (abs (a4 (4, k)) .gt. abs (a4 (2, k) - a4 (3, k))) then
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
        if (a6da .lt. - da2) then
            a4 (4, k) = 3. * (a4 (2, k) - a4 (1, k))
            a4 (3, k) = a4 (2, k) - a4 (4, k)
        elseif (a6da .gt. da2) then
            a4 (4, k) = 3. * (a4 (3, k) - a4 (1, k))
            a4 (2, k) = a4 (3, k) - a4 (4, k)
        endif
    endif
    
    call cs_limiters (km - 1, a4)
    
    ! -----------------------------------------------------------------------
    ! bottom:
    ! -----------------------------------------------------------------------
    
    a4 (2, km) = a4 (1, km)
    a4 (3, km) = a4 (1, km)
    a4 (4, km) = 0.
    
end subroutine cs_profile

! =======================================================================
! cubic spline (cs) limiters or boundary conditions
! a positive-definite constraint (iv = 0) is applied to tracers in every layer, 
! adjusting the top-most and bottom-most interface values to enforce positive.
! this subroutine is the same as cs_limiters in fv_mapz_mod where iv = 0.
! =======================================================================

subroutine cs_limiters (km, a4)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: km
    
    real, intent (inout) :: a4 (4, km) ! ppm array
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real, parameter :: r12 = 1. / 12.
    
    do k = 1, km
        if (a4 (1, k) .le. 0.) then
            a4 (2, k) = a4 (1, k)
            a4 (3, k) = a4 (1, k)
            a4 (4, k) = 0.
        else
            if (abs (a4 (3, k) - a4 (2, k)) .lt. - a4 (4, k)) then
                if ((a4 (1, k) + 0.25 * (a4 (3, k) - a4 (2, k)) ** 2 / a4 (4, k) + &
                    a4 (4, k) * r12) .lt. 0.) then
                    ! local minimum is negative
                    if (a4 (1, k) .lt. a4 (3, k) .and. a4 (1, k) .lt. a4 (2, k)) then
                        a4 (3, k) = a4 (1, k)
                        a4 (2, k) = a4 (1, k)
                        a4 (4, k) = 0.
                    elseif (a4 (3, k) .gt. a4 (2, k)) then
                        a4 (4, k) = 3. * (a4 (2, k) - a4 (1, k))
                        a4 (3, k) = a4 (2, k) - a4 (4, k)
                    else
                        a4 (4, k) = 3. * (a4 (3, k) - a4 (1, k))
                        a4 (2, k) = a4 (3, k) - a4 (4, k)
                    endif
                endif
            endif
        endif
    enddo
    
end subroutine cs_limiters

! =======================================================================
! time-implicit monotonic scheme
! =======================================================================

subroutine implicit_fall (dts, ks, ke, ze, vt, dp, q, precip, m1)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke + 1) :: ze
    
    real, intent (in), dimension (ks:ke) :: vt, dp
    
    real, intent (inout), dimension (ks:ke) :: q
    
    real, intent (inout) :: precip
    
    real, intent (out), dimension (ks:ke) :: m1
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real, dimension (ks:ke) :: dz, qm, dd
    
    do k = ks, ke
        dz (k) = ze (k) - ze (k + 1)
        dd (k) = dts * vt (k)
        q (k) = q (k) * dp (k)
    enddo
    
    qm (ks) = q (ks) / (dz (ks) + dd (ks))
    do k = ks + 1, ke
        qm (k) = (q (k) + qm (k - 1) * dd (k - 1)) / (dz (k) + dd (k))
    enddo
    
    do k = ks, ke
        qm (k) = qm (k) * dz (k)
    enddo
    
    m1 (ks) = q (ks) - qm (ks)
    do k = ks + 1, ke
        m1 (k) = m1 (k - 1) + q (k) - qm (k)
    enddo
    precip = precip + m1 (ke)
    
    do k = ks, ke
        q (k) = qm (k) / dp (k)
    enddo
    
end subroutine implicit_fall

! =======================================================================
! time-explicit monotonic scheme
! =======================================================================

subroutine explicit_fall (dts, ks, ke, ze, vt, dp, q, precip, m1)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke + 1) :: ze
    
    real, intent (in), dimension (ks:ke) :: vt, dp
    
    real, intent (inout), dimension (ks:ke) :: q
    
    real, intent (inout) :: precip
    
    real, intent (out), dimension (ks:ke) :: m1
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: n, k, nstep
    
    real, dimension (ks:ke) :: dz, qm, q0, dd
    
    do k = ks, ke
        dz (k) = ze (k) - ze (k + 1)
        dd (k) = dts * vt (k)
        q0 (k) = q (k) * dp (k)
    enddo
    
    nstep = 1 + int (maxval (dd / dz))
    do k = ks, ke
        dd (k) = dd (k) / nstep
        q (k) = q0 (k)
    enddo
    
    do n = 1, nstep
        qm (ks) = q (ks) - q (ks) * dd (ks) / dz (ks)
        do k = ks + 1, ke
            qm (k) = q (k) - q (k) * dd (k) / dz (k) + q (k - 1) * dd (k - 1) / dz (k - 1)
        enddo
        q = qm
    enddo
    
    m1 (ks) = q0 (ks) - qm (ks)
    do k = ks + 1, ke
        m1 (k) = m1 (k - 1) + q0 (k) - qm (k)
    enddo
    precip = precip + m1 (ke)
    
    do k = ks, ke
        q (k) = qm (k) / dp (k)
    enddo
    
end subroutine explicit_fall

! =======================================================================
! combine time-implicit monotonic scheme with the piecewise parabolic lagrangian scheme
! =======================================================================

subroutine implicit_lagrangian_fall (dts, ks, ke, zs, ze, zt, vt, dp, q, &
        precip, flux, sed_fac)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: zs, dts, sed_fac
    
    real, intent (in), dimension (ks:ke + 1) :: ze, zt
    
    real, intent (in), dimension (ks:ke) :: vt, dp
    
    real, intent (inout), dimension (ks:ke) :: q
    
    real, intent (inout) :: precip
    
    real, intent (out), dimension (ks:ke) :: flux

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    real :: pre0, pre1
    
    real, dimension (ks:ke) :: q0, q1, m0, m1

    q0 = q
    pre0 = precip
    
    call implicit_fall (dts, ks, ke, ze, vt, dp, q0, pre0, m0)

    q1 = q
    pre1 = precip

    call lagrangian_fall (ks, ke, zs, ze, zt, dp, q1, pre1, m1)

    q = q0 * sed_fac + q1 * (1.0 - sed_fac)
    flux = m0 * sed_fac + m1 * (1.0 - sed_fac)
    precip = pre0 * sed_fac + pre1 * (1.0 - sed_fac)
    
end subroutine implicit_lagrangian_fall
    
! =======================================================================
! vertical subgrid variability used for cloud ice and cloud water autoconversion
! edges: qe == qbar + / - dm
! =======================================================================

subroutine linear_prof (km, q, dm, z_var, h_var)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: km
    
    logical, intent (in) :: z_var
    
    real, intent (in) :: q (km), h_var
    
    real, intent (out) :: dm (km)
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: dq (km)
    
    if (z_var) then
        do k = 2, km
            dq (k) = 0.5 * (q (k) - q (k - 1))
        enddo
        dm (1) = 0.
        ! -----------------------------------------------------------------------
        ! use twice the strength of the positive definiteness limiter (Lin et al. 1994)
        ! -----------------------------------------------------------------------
        do k = 2, km - 1
            dm (k) = 0.5 * min (abs (dq (k) + dq (k + 1)), 0.5 * q (k))
            if (dq (k) * dq (k + 1) .le. 0.) then
                if (dq (k) .gt. 0.) then
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
            dm (k) = max (dm (k), 0.0, h_var * q (k))
        enddo
    else
        do k = 1, km
            dm (k) = max (0.0, h_var * q (k))
        enddo
    endif
    
end subroutine linear_prof

! =======================================================================
! accretion function, Lin et al. (1983)
! =======================================================================

function acr2d (qden, c, denfac, blin, mu)
    
    implicit none
    
    real :: acr2d
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: qden, c, denfac, blin, mu
    
    acr2d = denfac * c * exp ((2 + mu + blin) / (mu + 3) * log (6 * qden))
    
end function acr2d

! =======================================================================
! accretion function, Lin et al. (1983)
! =======================================================================

function acr3d (v1, v2, q1, q2, c, acco, acc1, acc2, den)
    
    implicit none
    
    real :: acr3d
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: v1, v2, c, den, q1, q2, acco (3), acc1, acc2
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: i
    
    real :: t1, t2, tmp, vdiff
    
    t1 = exp (1. / (acc1 + 3) * log (6 * q1 * den))
    t2 = exp (1. / (acc2 + 3) * log (6 * q2 * den))

    if (vdiffflag .eq. 1) vdiff = abs (v1 - v2)
    if (vdiffflag .eq. 2) vdiff = sqrt ((1.20 * v1 - 0.95 * v2) ** 2. + 0.08 * v1 * v2)
    if (vdiffflag .eq. 3) vdiff = sqrt ((1.00 * v1 - 1.00 * v2) ** 2. + 0.04 * v1 * v2)
    
    acr3d = c * vdiff / den
    
    tmp = 0
    do i = 1, 3
        tmp = tmp + acco (i) * exp ((6 + acc1 - i) * log (t1)) * exp ((acc2 + i - 1) * log (t2))
    enddo
    
    acr3d = acr3d * tmp
    
end function acr3d

! =======================================================================
! ventilation coefficient, Lin et al. (1983)
! =======================================================================

function vent_coeff (qden, c1, c2, denfac, blin, mu)
    
    implicit none
    
    real :: vent_coeff
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: qden, c1, c2, denfac, blin, mu
    
    vent_coeff = c1 + c2 * exp ((3 + 2 * mu + blin) / (mu + 3) / 2 * log (6 * qden)) * &
        sqrt (denfac) / exp ((1 + mu) / (mu + 3) * log (6 * qden))
    
end function vent_coeff

! =======================================================================
! sublimation or evaporation function, Lin et al. (1983)
! =======================================================================

function psub (t2, dq, qden, qsat, c, den, denfac, blin, mu, cpk, cvm)
    
    implicit none
    
    real :: psub
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: t2, dq, qden, qsat, c (5), den, denfac, blin, cpk, mu
    
    real (kind = r8), intent (in) :: cvm
    
    psub = c (1) * t2 * dq * exp ((1 + mu) / (mu + 3) * log (6 * qden)) * &
        vent_coeff (qden, c (2), c (3), denfac, blin, mu) / &
         (c (4) * t2 + c (5) * (cpk * cvm) ** 2 * qsat * den)
    
end function psub

! =======================================================================
! melting function, Lin et al. (1983)
! =======================================================================

function pmlt (tc, dq, qden, pxacw, pxacr, c, den, denfac, blin, mu, lcpk, icpk, cvm)
    
    implicit none
    
    real :: pmlt
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: tc, dq, qden, pxacw, pxacr, c (4), den, denfac, blin, lcpk, icpk, mu
    
    real (kind = r8), intent (in) :: cvm
    
    pmlt = (c (1) / (icpk * cvm) * tc / den - c (2) * lcpk / icpk * dq) * &
        exp ((1 + mu) / (mu + 3) * log (6 * qden)) * &
        vent_coeff (qden, c (3), c (4), denfac, blin, mu) + &
        c_liq / (icpk * cvm) * tc * (pxacw + pxacr)
    
end function pmlt

! =======================================================================
! sedimentation of horizontal momentum
! =======================================================================

subroutine sedi_uv (ks, ke, m1, dp, u, v)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in), dimension (ks:ke) :: m1, dp
    
    real, intent (inout), dimension (ks:ke) :: u, v
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    do k = ks + 1, ke
        u (k) = (dp (k) * u (k) + m1 (k - 1) * u (k - 1)) / (dp (k) + m1 (k - 1))
        v (k) = (dp (k) * v (k) + m1 (k - 1) * v (k - 1)) / (dp (k) + m1 (k - 1))
    enddo
    
end subroutine sedi_uv

! =======================================================================
! sedimentation of vertical momentum
! =======================================================================

subroutine sedi_w (ks, ke, m1, w, vt, dm)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in), dimension (ks:ke) :: m1, vt, dm
    
    real, intent (inout), dimension (ks:ke) :: w
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    w (ks) = w (ks) + m1 (ks) * vt (ks) / dm (ks)
    do k = ks + 1, ke
        w (k) = (dm (k) * w (k) + m1 (k - 1) * (w (k - 1) - vt (k - 1)) + m1 (k) * vt (k)) / &
             (dm (k) + m1 (k - 1))
    enddo
    
end subroutine sedi_w

! =======================================================================
! sedimentation of heat
! =======================================================================

subroutine sedi_heat (ks, ke, dm, m1, dz, tz, qv, ql, qr, qi, qs, qg, cw)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: cw
    
    real, intent (in), dimension (ks:ke) :: dm, m1, dz, qv, ql, qr, qi, qs, qg
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real, dimension (ks:ke) :: dgz, cv0
    
    do k = ks + 1, ke
        dgz (k) = - 0.5 * grav * (dz (k - 1) + dz (k))
        cv0 (k) = dm (k) * (cv_air + qv (k) * cv_vap + (qr (k) + ql (k)) * c_liq + &
             (qi (k) + qs (k) + qg (k)) * c_ice) + cw * (m1 (k) - m1 (k - 1))
    enddo
    
    do k = ks + 1, ke
        tz (k) = (cv0 (k) * tz (k) + m1 (k - 1) * (cw * tz (k - 1) + dgz (k))) / &
             (cv0 (k) + cw * m1 (k - 1))
    enddo
    
end subroutine sedi_heat

! =======================================================================
! fast saturation adjustments
! =======================================================================

subroutine fast_sat_adj (dtm, is, ie, ks, ke, hydrostatic, consv_te, &
        adj_vmr, te, dte, qv, ql, qr, qi, qs, qg, qa, qnl, qni, hs, delz, &
        pt, delp, q_con, cappa, gsize, last_step, condensation, &
        evaporation, deposition, sublimation, do_sat_adj)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: is, ie, ks, ke
    
    logical, intent (in) :: hydrostatic, last_step, consv_te, do_sat_adj
    
    real, intent (in) :: dtm
    
    real, intent (in), dimension (is:ie) :: hs, gsize
    
    real, intent (in), dimension (is:ie, ks:ke) :: qnl, qni
    
    real, intent (inout), dimension (is:ie, ks:ke) :: delp, delz, pt, te
    real, intent (inout), dimension (is:ie, ks:ke) :: qv, ql, qr, qi, qs, qg, qa
    
    real, intent (inout), dimension (is:, ks:) :: q_con, cappa
    
    real, intent (inout), dimension (is:ie) :: condensation, deposition
    real, intent (inout), dimension (is:ie) :: evaporation, sublimation
    
    real, intent (out), dimension (is:ie, ks:ke) :: adj_vmr

    real (kind = r8), intent (out), dimension (is:ie) :: dte

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    real, dimension (is:ie, ks:ke) :: ua, va, wa, prefluxw, prefluxr, prefluxi, prefluxs, prefluxg
    
    real, dimension (is:ie) :: water, rain, ice, snow, graupel
    
    real, dimension (is:ie, ks:ke) :: pcw, edw, oew, rrw, tvw
    real, dimension (is:ie, ks:ke) :: pci, edi, oei, rri, tvi
    real, dimension (is:ie, ks:ke) :: pcr, edr, oer, rrr, tvr
    real, dimension (is:ie, ks:ke) :: pcs, eds, oes, rrs, tvs
    real, dimension (is:ie, ks:ke) :: pcg, edg, oeg, rrg, tvg

    ! -----------------------------------------------------------------------
    ! initialization
    ! -----------------------------------------------------------------------
    
    ua = 0.0
    va = 0.0
    wa = 0.0
    
    water = 0.0
    rain = 0.0
    ice = 0.0
    snow = 0.0
    graupel = 0.0
    
    prefluxw = 0.0
    prefluxr = 0.0
    prefluxi = 0.0
    prefluxs = 0.0
    prefluxg = 0.0
    
    ! -----------------------------------------------------------------------
    ! major cloud microphysics driver
    ! -----------------------------------------------------------------------
    
    call mpdrv (hydrostatic, ua, va, wa, delp, pt, qv, ql, qr, qi, qs, qg, qa, &
        qnl, qni, delz, is, ie, ks, ke, dtm, water, rain, ice, snow, graupel, &
        gsize, hs, q_con, cappa, consv_te, adj_vmr, te, dte, pcw, edw, oew, rrw, tvw, &
        pci, edi, oei, rri, tvi, pcr, edr, oer, rrr, tvr, pcs, eds, oes, rrs, tvs, &
        pcg, edg, oeg, rrg, tvg, prefluxw, prefluxr, prefluxi, &
        prefluxs, prefluxg, condensation, deposition, evaporation, sublimation, &
        last_step, .true., do_sat_adj, .false.)
    
end subroutine fast_sat_adj

! =======================================================================
! rain freezing to form graupel, simple version
! =======================================================================

subroutine pgfr_simp (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, &
        lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, sink, fac_r2g
    
    fac_r2g = 1. - exp (- dts / tau_r2g)
    
    do k = ks, ke
        
        tc = tz (k) - tice
        
        if (tc .lt. 0. .and. qr (k) .gt. qcmin) then
            
            sink = (- tc * 0.025) ** 2 * qr (k)
            sink = min (qr (k), sink, - fac_r2g * tc / icpk (k))
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., 0., - sink, 0., 0., sink, te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine pgfr_simp

! =======================================================================
! snow melting to form cloud water and rain, simple version
! =======================================================================

subroutine psmlt_simp (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, &
        lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, tmp, sink, fac_smlt
    
    fac_smlt = 1. - exp (- dts / tau_smlt)
    
    do k = ks, ke
        
        tc = tz (k) - tice
        
        if (tc .ge. 0. .and. qs (k) .gt. qcmin) then
            
            sink = (tc * 0.1) ** 2 * qs (k)
            sink = min (qs (k), sink, fac_smlt * tc / icpk (k))
            tmp = min (sink, dim (qs_mlt, ql (k)))
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., tmp, sink - tmp, 0., - sink, 0., te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine psmlt_simp

! =======================================================================
! cloud water to rain autoconversion, simple version
! =======================================================================

subroutine praut_simp (ks, ke, dts, tz, qv, ql, qr, qi, qs, qg)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, sink, fac_l2r
    
    fac_l2r = 1. - exp (- dts / tau_l2r)
    
    do k = ks, ke
        
        tc = tz (k) - t_wfr
        
        if (tc .gt. 0 .and. ql (k) .gt. ql0_max) then
            
            sink = fac_l2r * (ql (k) - ql0_max)
            
            call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., - sink, sink, 0., 0., 0.)
            
        endif
        
    enddo

end subroutine praut_simp
    
! =======================================================================
! cloud ice to snow autoconversion, simple version
! =======================================================================

subroutine psaut_simp (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, den)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, sink, fac_i2s, qim
    
    fac_i2s = 1. - exp (- dts / tau_i2s)
    
    do k = ks, ke
        
        tc = tz (k) - tice
        
        qim = qi0_max / den (k)
        
        if (tc .lt. 0. .and. qi (k) .gt. qim) then
            
            sink = fac_i2s * (qi (k) - qim)
            
            call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., 0., 0., - sink, sink, 0.)
            
        endif
        
    enddo
    
end subroutine psaut_simp

! =======================================================================
! cloud radii diagnosis built for gfdl cloud microphysics
! =======================================================================

subroutine cld_eff_rad (is, ie, ks, ke, lsm, p, delp, t, qv, qw, qi, qr, qs, qg, qa, &
        qcw, qci, qcr, qcs, qcg, rew, rei, rer, res, reg, cld, cloud, snowd, &
        cnvw, cnvi, cnvc)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: is, ie, ks, ke
    
    real, intent (in), dimension (is:ie) :: lsm, snowd
    
    real, intent (in), dimension (is:ie, ks:ke) :: delp, t, p, cloud
    real, intent (in), dimension (is:ie, ks:ke) :: qv, qw, qi, qr, qs, qg, qa
    
    real, intent (in), dimension (is:ie, ks:ke), optional :: cnvw, cnvi, cnvc
    
    real, intent (inout), dimension (is:ie, ks:ke) :: qcw, qci, qcr, qcs, qcg
    real, intent (inout), dimension (is:ie, ks:ke) :: rew, rei, rer, res, reg
    real, intent (inout), dimension (is:ie, ks:ke) :: cld
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: i, k, ind
    
    real, dimension (is:ie, ks:ke) :: qmw, qmr, qmi, qms, qmg
    
    real :: dpg, rho, ccnw, mask, cor, tc, bw
    real :: lambdaw, lambdar, lambdai, lambdas, lambdag, rei_fac
    
    real :: retab (138) = (/ &
        0.05000, 0.05000, 0.05000, 0.05000, 0.05000, 0.05000, &
        0.05500, 0.06000, 0.07000, 0.08000, 0.09000, 0.10000, &
        0.20000, 0.30000, 0.40000, 0.50000, 0.60000, 0.70000, &
        0.80000, 0.90000, 1.00000, 1.10000, 1.20000, 1.30000, &
        1.40000, 1.50000, 1.60000, 1.80000, 2.00000, 2.20000, &
        2.40000, 2.60000, 2.80000, 3.00000, 3.20000, 3.50000, &
        3.80000, 4.10000, 4.40000, 4.70000, 5.00000, 5.30000, &
        5.60000, 5.92779, 6.26422, 6.61973, 6.99539, 7.39234, &
        7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930, &
        10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319, &
        15.0351, 15.8799, 16.7674, 17.6986, 18.6744, 19.6955, &
        20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125, &
        27.7895, 28.6450, 29.4167, 30.1088, 30.7306, 31.2943, &
        31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601, &
        34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078, &
        38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635, &
        42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221, &
        50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898, &
        65.5604, 71.2885, 75.4113, 79.7368, 84.2351, 88.8833, &
        93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424, &
        124.954, 130.630, 136.457, 142.446, 148.608, 154.956, &
        161.503, 168.262, 175.248, 182.473, 189.952, 197.699, &
        205.728, 214.055, 222.694, 231.661, 240.971, 250.639 /)
    
    qmw = qw
    qmi = qi
    qmr = qr
    qms = qs
    qmg = qg
    cld = cloud
    
    ! -----------------------------------------------------------------------
    ! merge convective cloud to total cloud
    ! -----------------------------------------------------------------------
    
    if (present (cnvw)) then
        qmw = qmw + cnvw
    endif
    if (present (cnvi)) then
        qmi = qmi + cnvi
    endif
    if (present (cnvc)) then
        cld = cnvc + (1 - cnvc) * cld
    endif
    
    ! -----------------------------------------------------------------------
    ! combine liquid and solid phases
    ! -----------------------------------------------------------------------
    
    if (liq_ice_combine) then
        do i = is, ie
            do k = ks, ke
                qmw (i, k) = qmw (i, k) + qmr (i, k)
                qmr (i, k) = 0.0
                qmi (i, k) = qmi (i, k) + qms (i, k) + qmg (i, k)
                qms (i, k) = 0.0
                qmg (i, k) = 0.0
            enddo
        enddo
    endif
    
    ! -----------------------------------------------------------------------
    ! combine snow and graupel
    ! -----------------------------------------------------------------------
    
    if (snow_grauple_combine) then
        do i = is, ie
            do k = ks, ke
                qms (i, k) = qms (i, k) + qmg (i, k)
                qmg (i, k) = 0.0
            enddo
        enddo
    endif
    
        
    do i = is, ie

        do k = ks, ke
            
            qmw (i, k) = max (qmw (i, k), qcmin)
            qmi (i, k) = max (qmi (i, k), qcmin)
            qmr (i, k) = max (qmr (i, k), qcmin)
            qms (i, k) = max (qms (i, k), qcmin)
            qmg (i, k) = max (qmg (i, k), qcmin)
            
            cld (i, k) = min (max (cld (i, k), 0.0), 1.0)
            
            mask = min (max (lsm (i), 0.0), 2.0)
            
            dpg = abs (delp (i, k)) / grav
            rho = p (i, k) / (rdgas * t (i, k) * (1. + zvir * qv (i, k)))
            
            tc = t (i, k) - tice
            
            if (rewflag .eq. 1) then
                
                ! -----------------------------------------------------------------------
                ! cloud water (Martin et al. 1994)
                ! -----------------------------------------------------------------------
                
                if (prog_ccn) then
                    ! boucher and lohmann (1995)
                    ccnw = (1.0 - abs (mask - 1.0)) * &
                         (10. ** 2.24 * (qa (i, k) * rho * 1.e9) ** 0.257) + &
                        abs (mask - 1.0) * &
                         (10. ** 2.06 * (qa (i, k) * rho * 1.e9) ** 0.48)
                else
                    ccnw = ccn_o * abs (mask - 1.0) + ccn_l * (1.0 - abs (mask - 1.0))
                endif
                
                if (qmw (i, k) .gt. qcmin) then
                    qcw (i, k) = dpg * qmw (i, k) * 1.0e3
                    rew (i, k) = exp (1.0 / 3.0 * log ((3.0 * qmw (i, k) * rho) / &
                         (4.0 * pi * rhow * ccnw))) * 1.0e4
                    rew (i, k) = max (rewmin, min (rewmax, rew (i, k)))
                else
                    qcw (i, k) = 0.0
                    rew (i, k) = rewmin
                endif
                
            endif
            
            if (rewflag .eq. 2) then
                
                ! -----------------------------------------------------------------------
                ! cloud water (Martin et al. 1994, gfdl revision)
                ! -----------------------------------------------------------------------
                
                if (prog_ccn) then
                    ! boucher and lohmann (1995)
                    ccnw = (1.0 - abs (mask - 1.0)) * &
                         (10. ** 2.24 * (qa (i, k) * rho * 1.e9) ** 0.257) + &
                        abs (mask - 1.0) * &
                         (10. ** 2.06 * (qa (i, k) * rho * 1.e9) ** 0.48)
                else
                    ccnw = 1.077 * ccn_o * abs (mask - 1.0) + 1.143 * ccn_l * (1.0 - abs (mask - 1.0))
                endif
                
                if (qmw (i, k) .gt. qcmin) then
                    qcw (i, k) = dpg * qmw (i, k) * 1.0e3
                    rew (i, k) = exp (1.0 / 3.0 * log ((3.0 * qmw (i, k) * rho) / &
                         (4.0 * pi * rhow * ccnw))) * 1.0e4
                    rew (i, k) = max (rewmin, min (rewmax, rew (i, k)))
                else
                    qcw (i, k) = 0.0
                    rew (i, k) = rewmin
                endif
                
            endif
            
            if (rewflag .eq. 3) then
                
                ! -----------------------------------------------------------------------
                ! cloud water (Kiehl et al. 1994)
                ! -----------------------------------------------------------------------
                
                if (qmw (i, k) .gt. qcmin) then
                    qcw (i, k) = dpg * qmw (i, k) * 1.0e3
                    rew (i, k) = 14.0 * abs (mask - 1.0) + &
                         (8.0 + (14.0 - 8.0) * min (1.0, max (0.0, - tc / 30.0))) * &
                         (1.0 - abs (mask - 1.0))
                    rew (i, k) = rew (i, k) + (14.0 - rew (i, k)) * &
                        min (1.0, max (0.0, snowd (i) / 1000.0))
                    rew (i, k) = max (rewmin, min (rewmax, rew (i, k)))
                else
                    qcw (i, k) = 0.0
                    rew (i, k) = rewmin
                endif
                
            endif
            
            if (rewflag .eq. 4) then
                
                ! -----------------------------------------------------------------------
                ! cloud water derived from PSD
                ! -----------------------------------------------------------------------
                
                if (qmw (i, k) .gt. qcmin) then
                    qcw (i, k) = dpg * qmw (i, k) * 1.0e3
                    call cal_pc_ed_oe_rr_tv (qmw (i, k), rho, blinw, muw, &
                        eda = edaw, edb = edbw, ed = rew (i, k))
                    rew (i, k) = rewfac * 0.5 * rew (i, k) * 1.0e6
                    rew (i, k) = max (rewmin, min (rewmax, rew (i, k)))
                else
                    qcw (i, k) = 0.0
                    rew (i, k) = rewmin
                endif
                
            endif
            
            if (reiflag .eq. 1) then
                
                ! -----------------------------------------------------------------------
                ! cloud ice (Heymsfield and Mcfarquhar 1996)
                ! -----------------------------------------------------------------------
                
                if (qmi (i, k) .gt. qcmin) then
                    qci (i, k) = dpg * qmi (i, k) * 1.0e3
                    rei_fac = log (1.0e3 * qmi (i, k) * rho)
                    if (tc .lt. - 50) then
                        rei (i, k) = beta / 9.917 * exp (0.109 * rei_fac) * 1.0e3
                    elseif (tc .lt. - 40) then
                        rei (i, k) = beta / 9.337 * exp (0.080 * rei_fac) * 1.0e3
                    elseif (tc .lt. - 30) then
                        rei (i, k) = beta / 9.208 * exp (0.055 * rei_fac) * 1.0e3
                    else
                        rei (i, k) = beta / 9.387 * exp (0.031 * rei_fac) * 1.0e3
                    endif
                    rei (i, k) = max (reimin, min (reimax, rei (i, k)))
                else
                    qci (i, k) = 0.0
                    rei (i, k) = reimin
                endif
                
            endif
            
            if (reiflag .eq. 2) then
                
                ! -----------------------------------------------------------------------
                ! cloud ice (Donner et al. 1997)
                ! -----------------------------------------------------------------------
                
                if (qmi (i, k) .gt. qcmin) then
                    qci (i, k) = dpg * qmi (i, k) * 1.0e3
                    if (tc .le. - 55) then
                        rei (i, k) = 15.41627
                    elseif (tc .le. - 50) then
                        rei (i, k) = 16.60895
                    elseif (tc .le. - 45) then
                        rei (i, k) = 32.89967
                    elseif (tc .le. - 40) then
                        rei (i, k) = 35.29989
                    elseif (tc .le. - 35) then
                        rei (i, k) = 55.65818
                    elseif (tc .le. - 30) then
                        rei (i, k) = 85.19071
                    elseif (tc .le. - 25) then
                        rei (i, k) = 72.35392
                    else
                        rei (i, k) = 92.46298
                    endif
                    rei (i, k) = max (reimin, min (reimax, rei (i, k)))
                else
                    qci (i, k) = 0.0
                    rei (i, k) = reimin
                endif
                
            endif
            
            if (reiflag .eq. 3) then
                
                ! -----------------------------------------------------------------------
                ! cloud ice (Fu 2007)
                ! -----------------------------------------------------------------------
                
                if (qmi (i, k) .gt. qcmin) then
                    qci (i, k) = dpg * qmi (i, k) * 1.0e3
                    rei (i, k) = 47.05 + tc * (0.6624 + 0.001741 * tc)
                    rei (i, k) = max (reimin, min (reimax, rei (i, k)))
                else
                    qci (i, k) = 0.0
                    rei (i, k) = reimin
                endif
                
            endif
            
            if (reiflag .eq. 4) then
                
                ! -----------------------------------------------------------------------
                ! cloud ice (Kristjansson et al. 2000)
                ! -----------------------------------------------------------------------
                
                if (qmi (i, k) .gt. qcmin) then
                    qci (i, k) = dpg * qmi (i, k) * 1.0e3
                    ind = min (max (int (t (i, k) - 136.0), 44), 138 - 1)
                    cor = t (i, k) - int (t (i, k))
                    rei (i, k) = retab (ind) * (1. - cor) + retab (ind + 1) * cor
                    rei (i, k) = max (reimin, min (reimax, rei (i, k)))
                else
                    qci (i, k) = 0.0
                    rei (i, k) = reimin
                endif
                
            endif
            
            if (reiflag .eq. 5) then
                
                ! -----------------------------------------------------------------------
                ! cloud ice (Wyser 1998)
                ! -----------------------------------------------------------------------
                
                if (qmi (i, k) .gt. qcmin) then
                    qci (i, k) = dpg * qmi (i, k) * 1.0e3
                    bw = - 2. + 1.e-3 * log10 (rho * qmi (i, k) / 50.e-3) * &
                        exp (1.5 * log (max (1.e-10, - tc)))
                    rei (i, k) = 377.4 + bw * (203.3 + bw * (37.91 + 2.3696 * bw))
                    rei (i, k) = max (reimin, min (reimax, rei (i, k)))
                else
                    qci (i, k) = 0.0
                    rei (i, k) = reimin
                endif
                
            endif
            
            if (reiflag .eq. 6) then
                
                ! -----------------------------------------------------------------------
                ! cloud ice (Sun and Rikus 1999, Sun 2001)
                ! -----------------------------------------------------------------------
                
                if (qmi (i, k) .gt. qcmin) then
                    qci (i, k) = dpg * qmi (i, k) * 1.0e3
                    rei_fac = log (1.0e3 * qmi (i, k) * rho)
                    rei (i, k) = 45.8966 * exp (0.2214 * rei_fac) + &
                        0.7957 * exp (0.2535 * rei_fac) * (tc + 190.0)
                    rei (i, k) = (1.2351 + 0.0105 * tc) * rei (i, k)
                    rei (i, k) = max (reimin, min (reimax, rei (i, k)))
                else
                    qci (i, k) = 0.0
                    rei (i, k) = reimin
                endif
                
            endif
            
            if (reiflag .eq. 7) then
                
                ! -----------------------------------------------------------------------
                ! cloud ice derived from PSD
                ! -----------------------------------------------------------------------
                
                if (qmi (i, k) .gt. qcmin) then
                    qci (i, k) = dpg * qmi (i, k) * 1.0e3
                    call cal_pc_ed_oe_rr_tv (qmi (i, k), rho, blini, mui, &
                        eda = edai, edb = edbi, ed = rei (i, k))
                    rei (i, k) = reifac * 0.5 * rei (i, k) * 1.0e6
                    rei (i, k) = max (reimin, min (reimax, rei (i, k)))
                else
                    qci (i, k) = 0.0
                    rei (i, k) = reimin
                endif
                
            endif
            
            if (rerflag .eq. 1) then
                
                ! -----------------------------------------------------------------------
                ! rain derived from PSD
                ! -----------------------------------------------------------------------
                
                if (qmr (i, k) .gt. qcmin) then
                    qcr (i, k) = dpg * qmr (i, k) * 1.0e3
                    call cal_pc_ed_oe_rr_tv (qmr (i, k), rho, blinr, mur, &
                        eda = edar, edb = edbr, ed = rer (i, k))
                    rer (i, k) = 0.5 * rer (i, k) * 1.0e6
                    rer (i, k) = max (rermin, min (rermax, rer (i, k)))
                else
                    qcr (i, k) = 0.0
                    rer (i, k) = rermin
                endif
                
            endif
            
            if (resflag .eq. 1) then
                
                ! -----------------------------------------------------------------------
                ! snow derived from PSD
                ! -----------------------------------------------------------------------
                
                if (qms (i, k) .gt. qcmin) then
                    qcs (i, k) = dpg * qms (i, k) * 1.0e3
                    call cal_pc_ed_oe_rr_tv (qms (i, k), rho, blins, mus, &
                        eda = edas, edb = edbs, ed = res (i, k))
                    res (i, k) = 0.5 * res (i, k) * 1.0e6
                    res (i, k) = max (resmin, min (resmax, res (i, k)))
                else
                    qcs (i, k) = 0.0
                    res (i, k) = resmin
                endif
                
            endif
            
            if (regflag .eq. 1) then
                
                ! -----------------------------------------------------------------------
                ! graupel derived from PSD
                ! -----------------------------------------------------------------------
                
                if (qmg (i, k) .gt. qcmin) then
                    qcg (i, k) = dpg * qmg (i, k) * 1.0e3
                    if (do_hail) then
                        call cal_pc_ed_oe_rr_tv (qmg (i, k), rho, blinh, muh, &
                            eda = edah, edb = edbh, ed = reg (i, k))
                    else
                        call cal_pc_ed_oe_rr_tv (qmg (i, k), rho, bling, mug, &
                            eda = edag, edb = edbg, ed = reg (i, k))
                    endif
                    reg (i, k) = 0.5 * reg (i, k) * 1.0e6
                    reg (i, k) = max (regmin, min (regmax, reg (i, k)))
                else
                    qcg (i, k) = 0.0
                    reg (i, k) = regmin
                endif
                
            endif
            
        enddo
        
    enddo
    
end subroutine cld_eff_rad

! =======================================================================
! radar reflectivity
! =======================================================================

subroutine rad_ref (is, ie, js, je, isd, ied, jsd, jed, q, pt, delp, peln, &
        delz, dbz, maxdbz, allmax, npz, ncnst, hydrostatic, zvir, &
        do_inline_mp, sphum, liq_wat, ice_wat, rainwat, snowwat, graupel, mp_top)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    logical, intent (in) :: hydrostatic, do_inline_mp
    
    integer, intent (in) :: is, ie, js, je, isd, ied, jsd, jed
    integer, intent (in) :: npz, ncnst, mp_top
    integer, intent (in) :: sphum, liq_wat, ice_wat, rainwat, snowwat, graupel
    
    real, intent (in) :: zvir
    
    real, intent (in), dimension (is:, js:, 1:) :: delz
    
    real, intent (in), dimension (isd:ied, jsd:jed, npz) :: pt, delp
    
    real, intent (in), dimension (isd:ied, jsd:jed, npz, ncnst) :: q
    
    real, intent (in), dimension (is:ie, npz + 1, js:je) :: peln
    
    real, intent (out) :: allmax
    
    real, intent (out), dimension (is:ie, js:je) :: maxdbz
    
    real, intent (out), dimension (is:ie, js:je, npz) :: dbz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: i, j, k
    
    real, parameter :: alpha = 0.224, mp_const = 200 * exp (1.6 * log (3.6e6))
    
    real (kind = r8) :: qden, z_e
    real :: fac_r, fac_s, fac_g
    
    real, dimension (npz) :: den, denfac, qmr, qms, qmg, vtr, vts, vtg
    
    ! -----------------------------------------------------------------------
    ! return if the microphysics scheme doesn't include rain
    ! -----------------------------------------------------------------------
    
    if (rainwat .lt. 1) return
    
    ! -----------------------------------------------------------------------
    ! initialization
    ! -----------------------------------------------------------------------
    
    dbz = - 20.
    maxdbz = - 20.
    allmax = - 20.
    
    ! -----------------------------------------------------------------------
    ! calculate radar reflectivity
    ! -----------------------------------------------------------------------
    
    do j = js, je
        do i = is, ie
            
            ! -----------------------------------------------------------------------
            ! air density
            ! -----------------------------------------------------------------------
            
            do k = 1, npz
                if (hydrostatic) then
                    den (k) = delp (i, j, k) / ((peln (i, k + 1, j) - peln (i, k, j)) * &
                        rdgas * pt (i, j, k) * (1. + zvir * q (i, j, k, sphum)))
                else
                    den (k) = - delp (i, j, k) / (grav * delz (i, j, k))
                endif
                qmr (k) = max (qcmin, q (i, j, k, rainwat))
                qms (k) = max (qcmin, q (i, j, k, snowwat))
                qmg (k) = max (qcmin, q (i, j, k, graupel))
            enddo
            
            do k = 1, npz
                denfac (k) = sqrt (den (npz) / den (k))
            enddo
            
            ! -----------------------------------------------------------------------
            ! fall speed
            ! -----------------------------------------------------------------------
            
            if (radr_flag .eq. 3) then
                call term_rsg (1, npz, qmr, den, denfac, vr_fac, blinr, &
                    mur, tvar, tvbr, vr_max, const_vr, vtr)
                vtr = vtr / rhor
            endif
            
            if (rads_flag .eq. 3) then
                call term_rsg (1, npz, qms, den, denfac, vs_fac, blins, &
                    mus, tvas, tvbs, vs_max, const_vs, vts)
                vts = vts / rhos
            endif
            
            if (radg_flag .eq. 3) then
                if (do_hail .and. .not. do_inline_mp) then
                    call term_rsg (1, npz, qmg, den, denfac, vg_fac, blinh, &
                        muh, tvah, tvbh, vg_max, const_vg, vtg)
                    vtg = vtg / rhoh
                else
                    call term_rsg (1, npz, qmg, den, denfac, vg_fac, bling, &
                        mug, tvag, tvbg, vg_max, const_vg, vtg)
                    vtg = vtg / rhog
                endif
            endif
            
            ! -----------------------------------------------------------------------
            ! radar reflectivity
            ! -----------------------------------------------------------------------
            
            do k = mp_top + 1, npz
                z_e = 0.
                
                if (rainwat .gt. 0) then
                    qden = den (k) * qmr (k)
                    if (qmr (k) .gt. qcmin) then
                        call cal_pc_ed_oe_rr_tv (qmr (k), den (k), blinr, mur, &
                            rra = rrar, rrb = rrbr, rr = fac_r)
                    else
                        fac_r = 0.0
                    endif
                    if (radr_flag .eq. 1 .or. radr_flag .eq. 2) then
                        z_e = z_e + fac_r * 1.e18
                    endif
                    if (radr_flag .eq. 3) then
                        z_e = z_e + mp_const * exp (1.6 * log (qden * vtr (k)))
                    endif
                endif
                
                if (snowwat .gt. 0) then
                    qden = den (k) * qms (k)
                    if (qms (k) .gt. qcmin) then
                        call cal_pc_ed_oe_rr_tv (qms (k), den (k), blins, mus, &
                            rra = rras, rrb = rrbs, rr = fac_s)
                    else
                        fac_s = 0.0
                    endif
                    if (rads_flag .eq. 1) then
                        if (pt (i, j, k) .lt. tice) then
                            z_e = z_e + fac_s * 1.e18 * alpha * (rhos / rhor) ** 2
                        else
                            z_e = z_e + fac_s * 1.e18 * alpha * (rhos / rhor) ** 2 / alpha
                        endif
                    endif
                    if (rads_flag .eq. 2) then
                        if (pt (i, j, k) .lt. tice) then
                            z_e = z_e + fac_s * 1.e18 * alpha * (rhos / rhoi) ** 2
                        else
                            z_e = z_e + fac_s * 1.e18
                        endif
                    endif
                    if (rads_flag .eq. 3) then
                        z_e = z_e + mp_const * exp (1.6 * log (qden * vts (k)))
                    endif
                endif
                
                if (graupel .gt. 0) then
                    qden = den (k) * qmg (k)
                    if (do_hail .and. .not. do_inline_mp) then
                        if (qmg (k) .gt. qcmin) then
                            call cal_pc_ed_oe_rr_tv (qmg (k), den (k), blinh, muh, &
                                rra = rrah, rrb = rrbh, rr = fac_g)
                        else
                            fac_g = 0.0
                        endif
                        if (radg_flag .eq. 1) then
                            if (pt (i, j, k) .lt. tice) then
                                z_e = z_e + fac_g * 1.e18 * alpha * (rhoh / rhor) ** 2
                            else
                                z_e = z_e + fac_g * 1.e18 * alpha * (rhoh / rhor) ** 2 / alpha
                            endif
                        endif
                        if (radg_flag .eq. 2) then
                            z_e = z_e + fac_g * 1.e18
                        endif
                    else
                        if (qmg (k) .gt. qcmin) then
                            call cal_pc_ed_oe_rr_tv (qmg (k), den (k), bling, mug, &
                                rra = rrag, rrb = rrbg, rr = fac_g)
                        else
                            fac_g = 0.0
                        endif
                        if (radg_flag .eq. 1) then
                            if (pt (i, j, k) .lt. tice) then
                                z_e = z_e + fac_g * 1.e18 * alpha * (rhog / rhor) ** 2
                            else
                                z_e = z_e + fac_g * 1.e18 * alpha * (rhog / rhor) ** 2 / alpha
                            endif
                        endif
                        if (radg_flag .eq. 2) then
                            z_e = z_e + fac_g * 1.e18
                        endif
                    endif
                    if (radg_flag .eq. 3) then
                        z_e = z_e + mp_const * exp (1.6 * log (qden * vtg (k)))
                    endif
                endif
                
                dbz (i, j, k) = 10. * log10 (max (0.01, z_e))
            enddo
            
            do k = mp_top + 1, npz
                maxdbz (i, j) = max (dbz (i, j, k), maxdbz (i, j))
            enddo
            
            allmax = max (maxdbz (i, j), allmax)
            
        enddo
    enddo
    
end subroutine rad_ref

! =======================================================================
! moist heat capacity, 3 input variables
! =======================================================================

function mhc3 (qv, q_liq, q_sol)
    
    implicit none
    
    real (kind = r8) :: mhc3
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: qv, q_liq, q_sol
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    mhc3 = one_r8 + qv * c1_vap + q_liq * c1_liq + q_sol * c1_ice
    
end function mhc3

! =======================================================================
! moist heat capacity, 4 input variables
! =======================================================================

function mhc4 (qd, qv, q_liq, q_sol)
    
    implicit none
    
    real (kind = r8) :: mhc4
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: qv, q_liq, q_sol
    
    real (kind = r8), intent (in) :: qd
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    mhc4 = qd + qv * c1_vap + q_liq * c1_liq + q_sol * c1_ice
    
end function mhc4

! =======================================================================
! moist heat capacity, 6 input variables
! =======================================================================

function mhc6 (qv, ql, qr, qi, qs, qg)
    
    implicit none
    
    real (kind = r8) :: mhc6
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: qv, ql, qr, qi, qs, qg
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    real :: q_liq, q_sol
    
    q_liq = ql + qr
    q_sol = qi + qs + qg
    mhc6 = mhc (qv, q_liq, q_sol)
    
end function mhc6

! =======================================================================
! moist total energy
! =======================================================================

function mte (qv, ql, qr, qi, qs, qg, tk, dp, moist_q)
    
    implicit none
    
    real (kind = r8) :: mte
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    logical, intent (in) :: moist_q
    
    real, intent (in) :: qv, ql, qr, qi, qs, qg, dp
    
    real (kind = r8), intent (in) :: tk
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    real :: q_liq, q_sol, q_cond
    
    real (kind = r8) :: cvm, con_r8
    
    q_liq = ql + qr
    q_sol = qi + qs + qg
    q_cond = q_liq + q_sol
    con_r8 = one_r8 - (qv + q_cond)
    if (moist_q) then
        cvm = mhc (con_r8, qv, q_liq, q_sol)
    else
        cvm = mhc (qv, q_liq, q_sol)
    endif
    mte = rgrav * cvm * c_air * tk * dp
    
end function mte

! =======================================================================
! moist total energy and total water
! =======================================================================

subroutine mtetw (ks, ke, qv, ql, qr, qi, qs, qg, tz, ua, va, wa, delp, &
        gsize, dte, vapor, water, rain, ice, snow, graupel, sen, stress, dts, &
        te, tw, te_b, tw_b, moist_q, hydrostatic, te_loss)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    logical, intent (in) :: moist_q, hydrostatic
    
    real, intent (in) :: gsize, vapor, water, rain, ice, snow, graupel, dts, sen, stress
    
    real, intent (in), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, ua, va, wa, delp
    
    real (kind = r8), intent (in) :: dte
    
    real (kind = r8), intent (in), dimension (ks:ke) :: tz
    
    real (kind = r8), intent (out) :: te_b, tw_b
    
    real (kind = r8), intent (out), optional :: te_loss
    
    real (kind = r8), intent (out), dimension (ks:ke) :: te, tw
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: q_cond
    
    real (kind = r8) :: con_r8
    
    real, dimension (ks:ke) :: q_liq, q_sol
    
    real (kind = r8), dimension (ks:ke) :: cvm

     do k = ks, ke
         q_liq (k) = ql (k) + qr (k)
         q_sol (k) = qi (k) + qs (k) + qg (k)
         q_cond = q_liq (k) + q_sol (k)
         con_r8 = one_r8 - (qv (k) + q_cond)
         if (moist_q) then
             cvm (k) = mhc (con_r8, qv (k), q_liq (k), q_sol (k))
         else
             cvm (k) = mhc (qv (k), q_liq (k), q_sol (k))
         endif
         te (k) = (cvm (k) * tz (k) + lv00 * qv (k) - li00 * q_sol (k)) * c_air
         if (hydrostatic) then
             te (k) = te (k) + 0.5 * (ua (k) ** 2 + va (k) ** 2)
         else
             te (k) = te (k) + 0.5 * (ua (k) ** 2 + va (k) ** 2 + wa (k) ** 2)
         endif
         te (k) = rgrav * te (k) * delp (k) * gsize ** 2.0
         tw (k) = rgrav * (qv (k) + q_cond) * delp (k) * gsize ** 2.0
     enddo
     te_b = (dte + (lv00 * c_air * vapor - li00 * c_air * (ice + snow + graupel)) * dts / 86400 + sen * dts + stress * dts) * gsize ** 2.0
     tw_b = (vapor + water + rain + ice + snow + graupel) * dts / 86400 * gsize ** 2.0

     if (present (te_loss)) then
          ! total energy change due to sedimentation and its heating
          te_loss = dte * gsize ** 2.0
     endif

end subroutine mtetw
    
! =======================================================================
! calculate heat capacities and latent heat coefficients
! =======================================================================

subroutine cal_mhc_lhc (ks, ke, qv, ql, qr, qi, qs, qg, q_liq, q_sol, &
        cvm, te8, tz, lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    
    real (kind = r8), intent (in), dimension (ks:ke) :: tz
    
    real, intent (out), dimension (ks:ke) :: q_liq, q_sol, lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (out), dimension (ks:ke) :: cvm, te8
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    do k = ks, ke
        q_liq (k) = ql (k) + qr (k)
        q_sol (k) = qi (k) + qs (k) + qg (k)
        cvm (k) = mhc (qv (k), q_liq (k), q_sol (k))
        te8 (k) = cvm (k) * tz (k) + lv00 * qv (k) - li00 * q_sol (k)
        lcpk (k) = (lv00 + d1_vap * tz (k)) / cvm (k)
        icpk (k) = (li00 + d1_ice * tz (k)) / cvm (k)
        tcpk (k) = (li20 + (d1_vap + d1_ice) * tz (k)) / cvm (k)
        tcp3 (k) = lcpk (k) + icpk (k) * min (1., dim (tice, tz (k)) / (tice - t_wfr))
    enddo

end subroutine cal_mhc_lhc
    
! =======================================================================
! update hydrometeors
! =======================================================================

subroutine update_qq (qv, ql, qr, qi, qs, qg, dqv, dql, dqr, dqi, dqs, dqg)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: dqv, dql, dqr, dqi, dqs, dqg
    
    real, intent (inout) :: qv, ql, qr, qi, qs, qg
    
    qv = qv + dqv
    ql = ql + dql
    qr = qr + dqr
    qi = qi + dqi
    qs = qs + dqs
    qg = qg + dqg
    
end subroutine update_qq

! =======================================================================
! update hydrometeors and temperature
! =======================================================================

subroutine update_qt (qv, ql, qr, qi, qs, qg, dqv, dql, dqr, dqi, dqs, dqg, te8, &
        cvm, tk, lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: dqv, dql, dqr, dqi, dqs, dqg
    
    real (kind = r8), intent (in) :: te8
    
    real, intent (inout) :: qv, ql, qr, qi, qs, qg
    
    real, intent (out) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (out) :: cvm, tk
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    qv = qv + dqv
    ql = ql + dql
    qr = qr + dqr
    qi = qi + dqi
    qs = qs + dqs
    qg = qg + dqg
    
    cvm = mhc (qv, ql, qr, qi, qs, qg)
    tk = (te8 - lv00 * qv + li00 * (qi + qs + qg)) / cvm
    
    lcpk = (lv00 + d1_vap * tk) / cvm
    icpk = (li00 + d1_ice * tk) / cvm
    tcpk = (li20 + (d1_vap + d1_ice) * tk) / cvm
    tcp3 = lcpk + icpk * min (1., dim (tice, tk) / (tice - t_wfr))
    
end subroutine update_qt

! =======================================================================
! calculation of particle concentration (pc), effective diameter (ed),
! optical extinction (oe), radar reflectivity factor (rr), and
! mass-weighted terminal velocity (tv)
! =======================================================================

subroutine cal_pc_ed_oe_rr_tv (q, den, blin, mu, pca, pcb, pc, eda, edb, ed, &
        oea, oeb, oe, rra, rrb, rr, tva, tvb, tv)

    implicit none

    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------

    real, intent (in) :: blin, mu

    real, intent (in) :: q, den

    real (kind = r8), intent (in), optional :: pca, pcb, eda, edb, oea, oeb, rra, rrb, tva, tvb

    real, intent (out), optional :: pc, ed, oe, rr, tv

    if (present (pca) .and. present (pcb) .and. present (pc)) then
        pc = pca / pcb * exp (mu / (mu + 3) * log (6 * den * q))
    endif
    if (present (eda) .and. present (edb) .and. present (ed)) then
        ed = eda / edb * exp (1. / (mu + 3) * log (6 * den * q))
    endif
    if (present (oea) .and. present (oeb) .and. present (oe)) then
        oe = oea / oeb * exp ((mu + 2) / (mu + 3) * log (6 * den * q))
    endif
    if (present (rra) .and. present (rrb) .and. present (rr)) then
        rr = rra / rrb * exp ((mu + 6) / (mu + 3) * log (6 * den * q))
    endif
    if (present (tva) .and. present (tvb) .and. present (tv)) then
        tv = tva / tvb * exp (blin / (mu + 3) * log (6 * den * q))
    endif

end subroutine cal_pc_ed_oe_rr_tv

! =======================================================================
! prepare saturation water vapor pressure tables
! =======================================================================

subroutine qs_init
    
    implicit none
    
    integer :: i
    
    if (.not. tables_are_initialized) then
        
        allocate (table0 (length))
        allocate (table1 (length))
        allocate (table2 (length))
        allocate (table3 (length))
        allocate (table4 (length))
        
        allocate (des0 (length))
        allocate (des1 (length))
        allocate (des2 (length))
        allocate (des3 (length))
        allocate (des4 (length))
        
        call qs_table0 (length)
        call qs_table1 (length)
        call qs_table2 (length)
        call qs_table3 (length)
        call qs_table4 (length)
        
        do i = 1, length - 1
            des0 (i) = max (0., table0 (i + 1) - table0 (i))
            des1 (i) = max (0., table1 (i + 1) - table1 (i))
            des2 (i) = max (0., table2 (i + 1) - table2 (i))
            des3 (i) = max (0., table3 (i + 1) - table3 (i))
            des4 (i) = max (0., table4 (i + 1) - table4 (i))
        enddo
        des0 (length) = des0 (length - 1)
        des1 (length) = des1 (length - 1)
        des2 (length) = des2 (length - 1)
        des3 (length) = des3 (length - 1)
        des4 (length) = des4 (length - 1)
        
        tables_are_initialized = .true.

    endif
    
end subroutine qs_init

! =======================================================================
! saturation water vapor pressure table, core function
! =======================================================================

subroutine qs_table_core (n, n_blend, do_smith_table, table)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: n, n_blend
    
    logical, intent (in) :: do_smith_table
    
    real, intent (out), dimension (n) :: table
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: i
    integer, parameter :: n_min = 1600
    
    real (kind = r8) :: delt = 0.1
    real (kind = r8) :: tmin, tem, esh
    real (kind = r8) :: wice, wh2o, fac0, fac1, fac2
    real (kind = r8) :: esbasw, tbasw, esbasi, a, b, c, d, e
    real (kind = r8) :: esupc (n_blend)
    
    esbasw = 1013246.0
    tbasw = tice + 100.
    esbasi = 6107.1
    tmin = tice - n_min * delt
    
    ! -----------------------------------------------------------------------
    ! compute es over ice between - (n_min * delt) deg C and 0 deg C
    ! -----------------------------------------------------------------------
    
    if (do_smith_table) then
        do i = 1, n_min
            tem = tmin + delt * real (i - 1)
            a = - 9.09718 * (tice / tem - 1.)
            b = - 3.56654 * log10 (tice / tem)
            c = 0.876793 * (1. - tem / tice)
            e = log10 (esbasi)
            table (i) = 0.1 * exp ((a + b + c + e) * log (10.))
        enddo
    else
        do i = 1, n_min
            tem = tmin + delt * real (i - 1)
            fac0 = (tem - tice) / (tem * tice)
            fac1 = fac0 * li2
            fac2 = (d2_ice * log (tem / tice) + fac1) / rvgas
            table (i) = e00 * exp (fac2)
        enddo
    endif
    
    ! -----------------------------------------------------------------------
    ! compute es over water between - (n_blend * delt) deg C and [ (n - n_min - 1) * delt] deg C
    ! -----------------------------------------------------------------------
    
    if (do_smith_table) then
        do i = 1, n - n_min + n_blend
            tem = tice + delt * (real (i - 1) - n_blend)
            a = - 7.90298 * (tbasw / tem - 1.)
            b = 5.02808 * log10 (tbasw / tem)
            c = - 1.3816e-7 * (exp ((1. - tem / tbasw) * 11.344 * log (10.)) - 1.)
            d = 8.1328e-3 * (exp ((tbasw / tem - 1.) * (- 3.49149) * log (10.)) - 1.)
            e = log10 (esbasw)
            esh = 0.1 * exp ((a + b + c + d + e) * log (10.))
            if (i .le. n_blend) then
                esupc (i) = esh
            else
                table (i + n_min - n_blend) = esh
            endif
        enddo
    else
        do i = 1, n - n_min + n_blend
            tem = tice + delt * (real (i - 1) - n_blend)
            fac0 = (tem - tice) / (tem * tice)
            fac1 = fac0 * lv0
            fac2 = (dc_vap * log (tem / tice) + fac1) / rvgas
            esh = e00 * exp (fac2)
            if (i .le. n_blend) then
                esupc (i) = esh
            else
                table (i + n_min - n_blend) = esh
            endif
        enddo
    endif
    
    ! -----------------------------------------------------------------------
    ! derive blended es over ice and supercooled water between - (n_blend * delt) deg C and 0 deg C
    ! -----------------------------------------------------------------------
    
    do i = 1, n_blend
        tem = tice + delt * (real (i - 1) - n_blend)
        wice = 1.0 / (delt * n_blend) * (tice - tem)
        wh2o = 1.0 / (delt * n_blend) * (tem - tice + delt * n_blend)
        table (i + n_min - n_blend) = wice * table (i + n_min - n_blend) + wh2o * esupc (i)
    enddo
    
end subroutine qs_table_core

! =======================================================================
! saturation water vapor pressure table 0, water only
! useful for idealized experiments
! it can also be used in warm rain microphyscis only
! =======================================================================

subroutine qs_table0 (n)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: n
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: i
    
    real (kind = r8) :: delt = 0.1
    real (kind = r8) :: tmin, tem, fac0, fac1, fac2
    
    tmin = tice - 160.
    
    ! -----------------------------------------------------------------------
    ! compute es over water only
    ! -----------------------------------------------------------------------
    
    do i = 1, n
        tem = tmin + delt * real (i - 1)
        fac0 = (tem - tice) / (tem * tice)
        fac1 = fac0 * lv0
        fac2 = (dc_vap * log (tem / tice) + fac1) / rvgas
        table0 (i) = e00 * exp (fac2)
    enddo
    
end subroutine qs_table0

! =======================================================================
! saturation water vapor pressure table 1, water and ice
! blended between -20 deg C and 0 deg C
! the most realistic saturation water vapor pressure for the full temperature range
! =======================================================================

subroutine qs_table1 (n)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: n
    
    call qs_table_core (n, 200, .false., table1)
    
end subroutine qs_table1

! =======================================================================
! saturation water vapor pressure table 2, water and ice
! same as table 1, but the blending is replaced with smoothing around 0 deg C
! it is not designed for mixed-phase cloud microphysics
! used for ice microphysics (< 0 deg C) or warm rain microphysics (> 0 deg C)
! =======================================================================

subroutine qs_table2 (n)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: n
    
    call qs_table_core (n, 0, .false., table2)
    
end subroutine qs_table2

! =======================================================================
! saturation water vapor pressure table 3, water and ice
! blended between -20 deg C and 0 deg C
! the same as table 1, but from smithsonian meteorological tables page 350
! =======================================================================

subroutine qs_table3 (n)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: n
    
    call qs_table_core (n, 200, .true., table3)
    
end subroutine qs_table3

! =======================================================================
! saturation water vapor pressure table 4, water and ice
! same as table 3, but the blending is replaced with smoothing around 0 deg C
! the same as table 2, but from smithsonian meteorological tables page 350
! =======================================================================

subroutine qs_table4 (n)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: n
    
    call qs_table_core (n, 0, .true., table4)
    
end subroutine qs_table4

! =======================================================================
! compute the saturated water pressure, core function
! =======================================================================

function es_core (length, tk, table, des)
    
    implicit none
    
    real :: es_core
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------

    integer, intent (in) :: length
    
    real, intent (in) :: tk
    
    real, intent (in), dimension (length) :: table, des
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: it
    
    real :: ap1, tmin
    
    if (.not. tables_are_initialized) call qs_init
    
    tmin = tice - 160.
    ap1 = 10. * dim (tk, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es_core = table (it) + (ap1 - it) * des (it)
    
end function es_core

! =======================================================================
! compute the saturated specific humidity, core function
! =======================================================================

function qs_core (length, tk, den, dqdt, table, des)
    
    implicit none
    
    real :: qs_core
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------

    integer, intent (in) :: length
    
    real, intent (in) :: tk, den
    
    real, intent (in), dimension (length) :: table, des
    
    real, intent (out) :: dqdt
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: it
    
    real :: ap1, tmin
    
    tmin = tice - 160.
    ap1 = 10. * dim (tk, tmin) + 1.
    ap1 = min (2621., ap1)
    qs_core = es_core (length, tk, table, des) / (rvgas * tk * den)
    it = ap1 - 0.5
    dqdt = 10. * (des (it) + (ap1 - it) * (des (it + 1) - des (it))) / (rvgas * tk * den)
    
end function qs_core

! =======================================================================
! compute the saturated water pressure based on table 0, water only
! useful for idealized experiments
! it can also be used in warm rain microphyscis only
! =======================================================================

function wes_t (tk)
    
    implicit none
    
    real :: wes_t
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: tk
    
    wes_t = es_core (length, tk, table0, des0)
    
end function wes_t

! =======================================================================
! compute the saturated water pressure based on table 1, water and ice
! the most realistic saturation water vapor pressure for the full temperature range
! =======================================================================

function mes_t (tk)
    
    implicit none
    
    real :: mes_t
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: tk
    
    mes_t = es_core (length, tk, table1, des1)
    
end function mes_t

! =======================================================================
! compute the saturated water pressure based on table 2, water and ice
! it is not designed for mixed-phase cloud microphysics
! used for ice microphysics (< 0 deg C) or warm rain microphysics (> 0 deg C)
! =======================================================================

function ies_t (tk)
    
    implicit none
    
    real :: ies_t
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: tk
    
    ies_t = es_core (length, tk, table2, des2)
    
end function ies_t

! =======================================================================
! compute the saturated specific humidity based on table 0, water only
! useful for idealized experiments
! it can also be used in warm rain microphyscis only
! =======================================================================

function wqs_trho (tk, den, dqdt)
    
    implicit none
    
    real :: wqs_trho
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: tk, den
    
    real, intent (out) :: dqdt
    
    wqs_trho = qs_core (length, tk, den, dqdt, table0, des0)
    
end function wqs_trho

! =======================================================================
! compute the saturated specific humidity based on table 1, water and ice
! the most realistic saturation water vapor pressure for the full temperature range
! =======================================================================

function mqs_trho (tk, den, dqdt)
    
    implicit none
    
    real :: mqs_trho
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: tk, den
    
    real, intent (out) :: dqdt
    
    mqs_trho = qs_core (length, tk, den, dqdt, table1, des1)
    
end function mqs_trho

! =======================================================================
! compute the saturated specific humidity based on table 2, water and ice
! it is not designed for mixed-phase cloud microphysics
! used for ice microphysics (< 0 deg C) or warm rain microphysics (> 0 deg C)
! =======================================================================

function iqs_trho (tk, den, dqdt)
    
    implicit none
    
    real :: iqs_trho
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: tk, den
    
    real, intent (out) :: dqdt
    
    iqs_trho = qs_core (length, tk, den, dqdt, table2, des2)
    
end function iqs_trho

! =======================================================================
! compute the saturated specific humidity based on table 0, water only
! useful for idealized experiments
! it can also be used in warm rain microphyscis only
! =======================================================================

function wqs_ptqv (tk, pa, qv, dqdt)
    
    implicit none
    
    real :: wqs_ptqv
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: tk, pa, qv
    
    real, intent (out) :: dqdt
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    real :: den
    
    den = pa / (rdgas * tk * (1. + zvir * qv))
    
    wqs_ptqv = wqs (tk, den, dqdt)
    
end function wqs_ptqv

! =======================================================================
! compute the saturated specific humidity based on table 1, water and ice
! the most realistic saturation water vapor pressure for the full temperature range
! =======================================================================

function mqs_ptqv (tk, pa, qv, dqdt)
    
    implicit none
    
    real :: mqs_ptqv
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: tk, pa, qv
    
    real, intent (out) :: dqdt
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    real :: den
    
    den = pa / (rdgas * tk * (1. + zvir * qv))
    
    mqs_ptqv = mqs (tk, den, dqdt)
    
end function mqs_ptqv

! =======================================================================
! compute the saturated specific humidity based on table 2, water and ice
! it is not designed for mixed-phase cloud microphysics
! used for ice microphysics (< 0 deg C) or warm rain microphysics (> 0 deg C)
! =======================================================================

function iqs_ptqv (tk, pa, qv, dqdt)
    
    implicit none
    
    real :: iqs_ptqv
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: tk, pa, qv
    
    real, intent (out) :: dqdt
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    real :: den
    
    den = pa / (rdgas * tk * (1. + zvir * qv))
    
    iqs_ptqv = iqs (tk, den, dqdt)
    
end function iqs_ptqv

! =======================================================================
! compute the saturated specific humidity based on table 1, water and ice
! the most realistic saturation water vapor pressure for the full temperature range
! it is the 3d version of "mqs"
! =======================================================================

subroutine mqs3d (im, km, ks, tk, pa, qv, qs, dqdt)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: im, km, ks
    
    real, intent (in), dimension (im, ks:km) :: tk, pa, qv
    
    real, intent (out), dimension (im, ks:km) :: qs
    
    real, intent (out), dimension (im, ks:km), optional :: dqdt
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: i, k
    
    real :: dqdt0
    
    if (present (dqdt)) then
        do k = ks, km
            do i = 1, im
                qs (i, k) = mqs (tk (i, k), pa (i, k), qv (i, k), dqdt (i, k))
            enddo
        enddo
    else
        do k = ks, km
            do i = 1, im
                qs (i, k) = mqs (tk (i, k), pa (i, k), qv (i, k), dqdt0)
            enddo
        enddo
    endif
    
end subroutine mqs3d

! =======================================================================
! compute wet buld temperature, core function
! Knox et al. (2017)
! =======================================================================

function wet_bulb_core (qv, tk, den, lcp)
    
    implicit none
    
    real :: wet_bulb_core
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: qv, tk, den, lcp
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    logical :: do_adjust = .false.
    
    real :: factor = 1. / 3.
    real :: qsat, tp, dqdt
    
    wet_bulb_core = tk
    qsat = wqs (wet_bulb_core, den, dqdt)
    tp = factor * (qsat - qv) / (1. + lcp * dqdt) * lcp
    wet_bulb_core = wet_bulb_core - tp
    
    if (do_adjust .and. tp .gt. 0.0) then
        qsat = wqs (wet_bulb_core, den, dqdt)
        tp = (qsat - qv) / (1. + lcp * dqdt) * lcp
        wet_bulb_core = wet_bulb_core - tp
    endif
    
end function wet_bulb_core

! =======================================================================
! compute wet buld temperature, dry air case
! =======================================================================

function wet_bulb_dry (qv, tk, den)
    
    implicit none
    
    real :: wet_bulb_dry
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: qv, tk, den
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    real :: lcp
    
    lcp = hlv / cp_air
    
    wet_bulb_dry = wet_bulb_core (qv, tk, den, lcp)
    
end function wet_bulb_dry

! =======================================================================
! compute wet buld temperature, moist air case
! =======================================================================

function wet_bulb_moist (qv, ql, qi, qr, qs, qg, tk, den)
    
    implicit none
    
    real :: wet_bulb_moist
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: qv, ql, qi, qr, qs, qg, tk, den
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    real :: lcp, q_liq, q_sol
    
    real (kind = r8) :: cvm
    
    q_liq = ql + qr
    q_sol = qi + qs + qg
    cvm = mhc (qv, q_liq, q_sol)
    lcp = (lv00 + d1_vap * tk) / cvm
    
    wet_bulb_moist = wet_bulb_core (qv, tk, den, lcp)
    
end function wet_bulb_moist

end module gfdl_mp_mod
