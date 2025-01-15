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
! Scale-Aware Aerosol-Aware Mass-Flux (SA-AAMP) Convection Scheme
! This code was originally from GFSv16. It was later rewritten as an inline scheme.
! Developers: Jongil Han, Linjiong Zhou, and the GFDL FV3 Team
! References: Han and Pan (2011), Han et al. (2017), Han and Bretherton (2019)
! =======================================================================

module sa_aamf_mod

    use fms_mod, only: check_nml_error
    use gfdl_mp_mod, only: mqs

    implicit none

    private

    ! -----------------------------------------------------------------------
    ! public subroutines, functions, and variables
    ! -----------------------------------------------------------------------

    public :: sa_aamf_init
    public :: sa_aamf_deep
    public :: sa_aamf_shal

    ! -----------------------------------------------------------------------
    ! physics constants
    ! -----------------------------------------------------------------------

    real, parameter :: grav = 9.80665 ! acceleration due to gravity (m/s^2), ref: IFS

    real, parameter :: rdgas = 287.05 ! gas constant for dry air (J/kg/K): ref: GFDL, GFS
    real, parameter :: rvgas = 461.50 ! gas constant for water vapor (J/kg/K): ref: GFDL, GFS

    real, parameter :: zvir = rvgas / rdgas - 1. ! 0.6077667316114637
    real, parameter :: eps = rdgas / rvgas ! 0.6219934994582882
    real, parameter :: epsm1 = rdgas / rvgas - 1. ! -0.3780065005417118

    real, parameter :: tice = 273.15 ! freezing temperature (K): ref: GFDL, GFS

    real, parameter :: cp_air = 1004.6 ! heat capacity of dry air at constant pressure (J/kg/K): ref: GFDL, GFS
    real, parameter :: cp_vap = 4.0 * rvgas ! 1846.0885419672554, heat capacity of water vapor at constnat pressure (J/kg/K)

    real, parameter :: c_liq = 4.218e3 ! heat capacity of water at 0 deg C (J/kg/K), ref: IFS

    real, parameter :: hlv = 2.5e6 ! latent heat of evaporation at 0 deg C (J/kg): ref: GFDL, GFS

    real, parameter :: qcmin = 1.0e-15 ! min value for cloud condensates (kg/kg)

    ! -----------------------------------------------------------------------
    ! namelist parameters
    ! -----------------------------------------------------------------------

    ! mass flux deep convection

    real :: clam_deep    = 0.1   ! c_e for deep convection (Han and Pan, 2011, eq(6))
    real :: c0s_deep     = 0.002 ! conversion parameter of detrainment from liquid water into convetive precipitaiton
    real :: c1_deep      = 0.002 ! conversion parameter of detrainment from liquid water into grid-scale cloud water
    real :: pgcon_deep   = 0.55  ! control the reduction in momentum transport
                                 ! 0.7 : Gregory et al. (1997, QJRMS)
                                 ! 0.55: Zhang & Wu (2003, JAS)
    real :: asolfac_deep = 0.89  ! aerosol-aware parameter based on Lim & Hong (2012)
                                 ! asolfac_deep= cx / c0s_deep(=.002)
                                 ! cx = min([-0.7 ln(Nccn) + 24]*1.e-4, c0s_deep)
                                 ! Nccn: CCN number concentration in cm^(-3)
                                 ! Until a realistic Nccn is provided, typical Nccns are assumed
                                 ! as Nccn=100 for sea and Nccn=7000 for land
    real :: evfact_deep  = 0.3   ! evaporation factor
    real :: evfactl_deep = 0.3   ! evaporation factor over land
    real :: betal_deep   = 0.05  ! downdraft heat flux contribution over land
    real :: betas_deep   = 0.05  ! downdraft heat flux contribution over ocean
    real :: dxcrtas_deep = 8.e3  ! the threshold value (unit: m) for the quasi-equilibrium assumption of Arakawa-Schubert
    real :: cthk_deep    = 200.  ! min cloud top for deep convection
    real :: betaw_deep   = 0.03  ! ratio between cloud base mass flux and mean updraft (eq 6 in Han et al 2017)

    ! mass flux shallow convection

    real :: clam_shal    = 0.3   ! c_e for shallow convection (Han and Pan, 2011, eq(6))
    real :: c0s_shal     = 0.002 ! conversion parameter of detrainment from liquid water into convetive precipitaiton
    real :: c1_shal      = 5.e-4 ! conversion parameter of detrainment from liquid water into grid-scale cloud water
    real :: pgcon_shal   = 0.55  ! control the reduction in momentum transport
                                 ! 0.7 : Gregory et al. (1997, QJRMS)
                                 ! 0.55: Zhang & Wu (2003, JAS)
    real :: asolfac_shal = 0.89  ! aerosol-aware parameter based on Lim & Hong (2012)
                                 ! asolfac_shal= cx / c0s_shal(=.002)
                                 ! cx = min([-0.7 ln(Nccn) + 24]*1.e-4, c0s_shal)
                                 ! Nccn: CCN number concentration in cm^(-3)
                                 ! Until a realistic Nccn is provided, typical Nccns are assumed
                                 ! as Nccn=100 for sea and Nccn=7000 for land
    real :: evfact_shal  = 0.3   ! rain evaporation efficiency over the ocean
    real :: evfactl_shal = 0.3   ! rain evaporation efficiency over the land
    real :: cthk_shal    = 200.  ! max cloud top for shallow convection
    real :: top_shal     = 0.7   ! max cloud height for shallow convection (P/Ps < top_shal)
    real :: betaw_shal   = 0.03  ! ratio between cloud base mass flux and mean updraft (eq 6 in Han et al 2017)
    real :: dxcrt_shal   = 15.e3 ! critical resolution for calculating scale-aware cloud base mass flux

    ! for both convections

    logical :: use_tke_conv    = .false. ! flag for adjusting entrainment/detrainment rates in conv
    logical :: use_shear_conv  = .false. ! flag for considering shear effect for wu/wd in conv
    logical :: limit_shal_conv = .false. ! flag for constraining shal conv based on diagnosed cloud depth/top


    ! -----------------------------------------------------------------------
    ! namelist
    ! -----------------------------------------------------------------------

    namelist / sa_aamf_nml / &
        clam_deep, c0s_deep, c1_deep, pgcon_deep, asolfac_deep, evfact_deep, evfactl_deep, &
        betal_deep, betas_deep, dxcrtas_deep, cthk_deep, betaw_deep, clam_shal, c0s_shal, &
        c1_shal, pgcon_shal, asolfac_shal, evfact_shal, evfactl_shal, cthk_shal, top_shal, &
        betaw_shal, dxcrt_shal, use_tke_conv, use_shear_conv, limit_shal_conv

contains

! =======================================================================
! SAMP initialization
! =======================================================================

subroutine sa_aamf_init (input_nml_file, logunit)

    implicit none

    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------

    integer, intent (in) :: logunit

    character (len = *), intent (in) :: input_nml_file (:)

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: ios, ierr

    ! -----------------------------------------------------------------------
    ! read namelist
    ! -----------------------------------------------------------------------

    read (input_nml_file, nml = sa_aamf_nml, iostat = ios)
    ierr = check_nml_error (ios, 'sa_aamf_nml')

    ! -----------------------------------------------------------------------
    ! write namelist to log file
    ! -----------------------------------------------------------------------

    write (logunit, *) " ================================================================== "
    write (logunit, *) "sa_aamf_mod"
    write (logunit, nml = sa_aamf_nml)

end subroutine sa_aamf_init

! =======================================================================
! Scale-Aware Aerosol-Aware Mass-Flux Deep Convection

! The Scale-Aware Mass-Flux (SAMF) deep convection scheme is an updated version of the previous
! Simplified Arakawa-Schubert (SAS) scheme with scale and aerosol awareness and parameterizes the
! effect of deep convection on the environment (represented by the model state variables) in the
! following way.
!
! First, a simple cloud model is used to determine the change in model state variables due to one
! entraining/detraining cloud type, per unit cloud-base mass flux. Next, the total change in state
! variables is retrieved by determining the actual cloud base mass flux using the quasi-equilibrium
! assumption (for grid sizes larger than a threshold value [currently set to 8 km]) or a mean
! updraft velocity (for grid sizes smaller than the threshold value). With a scale-aware
! parameterization, the cloud mass flux decreases with increasing grid resolution. A simple
! aerosol-aware parameterization is employed, where rain conversion in the convective updraft is
! modified by aerosol number concentration. The name SAS is replaced with SAMF as for the smaller
! grid sizes, the parameterization does not use Arakawa-Schubert's quasi-equilibrium assumption any
! longer where the cloud work function (interpreted as entrainment-moderated Convective Available
! Potential Energy [CAPE]) by the large scale dynamics is in balance with the consumption of the
! cloud work function by the convection.
!
! The SAS scheme uses the working concepts put forth in Arakawa and Schubert (1974) but includes
! modifications and simplifications from Grell (1993) such as saturated downdrafts and only one
! cloud type (the deepest possible), rather than a spectrum based on cloud top heights or assumed
! entrainment rates. The scheme was implemented for the GFS in 1995 by Pan and Wu, with further
! modifications discussed in Han and Pan (2011), including the calculation of cloud top, a greater
! CFL-criterion-based maximum cloud base mass flux, updated cloud model entrainment and detrainment,
! improved convective transport of horizontal momentum, a more general triggering function, and
! the inclusion of convective overshooting.
!
! The SAMF scheme updates the SAS scheme with scale- and aerosol-aware parameterizations from
! Han et al. (2017) based on the studies by Arakawa and Wu (2013) and Grell and Freitas (2014) for
! scale awareness and by Lim (2011) for aerosol awareness. The ratio of advective time to
! convective turnover time is also taken into account for the scale-aware parameterization.
! Along with the scale- and aerosol-aware parameterizations, more changes are made to the SAMF
! scheme. The cloud base mass-flux computation is modified to use convective turnover time as the
! convective adjustment time scale. The rain conversion rate is modified to decrease with
! decreasing air temperature above the freezing level. Convective inhibition in the sub-cloud layer
! is used as an additional trigger condition. Convective cloudiness is enhanced by considering
! suspended cloud condensate in the updraft. The lateral entrainment is also enhanced to more
! strongly suppress convection in a drier environment.
!
! In further update for FY19 GFS implementation, interaction with Turbulent Kinetic Energy (TKE),
! which is a prognostic variable used in a scale-aware TKE-based moist EDMF vertical turbulent
! mixing scheme, is included. Entrainment rates in updrafts and downdrafts are proportional to sub-
! cloud mean TKE. TKE is transported by cumulus convection. TKE contribution from cumulus
! convection is deduced from cumulus mass flux. On the other hand, tracers such as ozone and
! aerosol are also transported by cumulus convection.
!
! Occasional model crashes have been occurred when stochastic physics is on, due to too much
! convective cooling and heating tendencies near the cumulus top which are amplified by stochastic
! physics. To reduce too much convective cooling at the cloud top, the convection schemes have been
! modified for the rain conversion rate, entrainment and detrainment rates, overshooting layers,
! and maximum allowable cloudbase mass flux (as of JUNE 2018) .
!
! For grid sizes larger than threshold value, as in Grell (1993) , the SAMF deep convection scheme
! can be described in terms of three types of "controls": static, dynamic, and feedback. The static
! control component consists of the simple entraining/detraining updraft/downdraft cloud model and
! is used to determine the cloud properties, convective precipitation, as well as the convective
! cloud top height. The dynamic control is the determination of the potential energy available for
! convection to "consume", or how primed the large-scale environment is for convection to occur due
! to changes by the dyanmics of the host model. The feedback control is the determination of how the
! parameterized convection changes the large-scale environment (the host model state variables)
! given the changes to the state variables per unit cloud base mass flux calculated in the static
! control portion and the deduced cloud base mass flux determined from the dynamic control.
!
! For grid sizes smaller than threshold value, the cloud base mass flux in the SAMF scheme is
! determined by the cumulus updraft velocity averaged ove the whole cloud depth (Han et al., 2017),
! which in turn, determines changes of the large-scale environment due to the cumulus convection.
!
! \param[in] IM      number of used points
! \param[in] KM      vertical layer dimension
! \param[in] DELT    physics time step in seconds
! \param[in] NTK     index for tke
! \param[in] NTR     total number of tracers including tke
! \param[in] DELP    pressure difference between level k and k + 1 (pa)
! \param[in] PRSLP   mean layer presure (pa)
! \param[in] PSP     surface pressure (pa)
! \param[in] PHIL    layer geopotential (\f$m^2 / s^2\f$)
! \param[in] QTR     tracer array including cloud condensate (\f$kg / kg\f$)
! \param[inout] QL   cloud water or ice (kg / kg)
! \param[inout] Q1   updated tracers (kg / kg)
! \param[inout] T1   updated temperature (k)
! \param[inout] U1   updated zonal wind (\f$m s^{ - 1}\f$)
! \param[inout] V1   updated meridional wind (\f$m s^{ - 1}\f$)
! \param[out] RN     convective rain (m)
! \param[out] KBOT   index for cloud base
! \param[out] KTOP   index for cloud top
! \param[out] KCNV   flag to denote deep convection (0 = no, 1 = yes)
! \param[in] ISLIMSK sea / land / ice mask (= 0 / 1 / 2)
! \param[in] GSIZE   size of grid box (\f$m\f$)
! \param[in] DOT     layer mean vertical velocity (pa / s)
! \param[in] NCLOUD  number of cloud species
! \param[out] UD_MF  updraft mass flux multiplied by time step (\f$kg / m^2\f$)
! \param[out] DD_MF  downdraft mass flux multiplied by time step (\f$kg / m^2\f$)
! \param[out] DT_MF  ud_mf at cloud top (\f$kg / m^2\f$)
! \param[out] CNVW   convective cloud water (kg / kg)
! \param[out] CNVC   convective cloud cover (unitless)
!
! General Algorithm
! # Compute preliminary quantities needed for static, dynamic, and feedback control portions of the algorithm.
! # Perform calculations related to the updraft of the entraining/detraining cloud model ("static control") .
! # Perform calculations related to the downdraft of the entraining/detraining cloud model ("static control") .
!
! # For grid sizes larger than the threshold value (currently 8 km) :
! 1) Using the updated temperature and moisture profiles that were modified by the convection on a short time-scale,
!    recalculate the total cloud work function to determine the change in the cloud work function due to convection,
!    or the stabilizing effect of the cumulus.
! 2) For the "dynamic control", using a reference cloud work function, estimate the change in cloud work function
!    due to the large-scale dynamics. following the quasi-equilibrium assumption, calculate the cloud base mass flux
!    required to keep the large-scale convective destabilization in balance with the stabilization effect of the convection.
!
! # For grid sizes smaller than the threshold value (currently 8 km) :
! 1) Compute the cloud base mass flux using the cumulus updraft velocity averaged ove the whole cloud depth.
! # For scale awareness, the updraft fraction (sigma) is obtained as a function of cloud base entrainment.
!   Then, the final cloud base mass flux is obtained by the original mass flux multiplied by the (1sigma) 2
! # For the "feedback control", calculate updated values of the state variables by multiplying the cloud base
!   mass flux and the tendencies calculated per unit cloud base mass flux from the static control.
! =======================================================================

subroutine sa_aamf_deep (im, km, delt, itc, ntc, ntw, nti, ntk, ntr, delp, &
        prslp, psp, phil, qtr, q1, t1, u1, v1, qr, fscav, rn, kbot, ktop, &
        kcnv, islimsk, gsize, dot, ncloud, ud_mf, dd_mf, dt_mf, cnvw, cnvc)

    implicit none

    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------

    integer, intent (in) :: im, km, itc, ntc, ntw, nti, ntk, ntr, ncloud, islimsk (im)

    real, intent (in) :: delt
    real, intent (in) :: psp (im), delp (im, km), &
        prslp (im, km), gsize (im), dot (im, km), phil (im, km)
    real, intent (in) :: fscav (ntc)

    integer, intent (inout) :: kcnv (im)

    real, intent (inout) :: qtr (im, km, ntr + 2), &
        q1 (im, km), t1 (im, km), u1 (im, km), v1 (im, km)

    integer, intent (out) :: kbot (im), ktop (im)

    real, intent (out) :: rn (im), qr (im, km)
    real, intent (out), optional :: cnvw (im, km), cnvc (im, km), &
        ud_mf (im, km), dd_mf (im, km), dt_mf (im, km)

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: i, indx, jmn, k, kk, km1, n

    real :: clamd, tkemx, tkemn, dtke, &
        dbeta, betamx, betamn, &
        cxlame, cxlamd, &
        xlamde, xlamdd, &
        crtlame, crtlamd

    ! real :: detad

    real :: adw, aup, aafac, beta, d0, &
        dellat, delta, desdt, dg, &
        dh, dhh, dp, &
        dq, dqsdp, dqsdt, dt, &
        dt2, dtmax, dtmin, &
        dxcrtuf, &
        dv1h, dv2h, dv3h, &
        dv1q, dv2q, dv3q, &
        dz, dz1, e1, edtmax, &
        edtmaxl, edtmaxs, el2orc, elocp, &
        es, etah, &
        dthk, &
        evef, fact1, &
        fact2, factor, &
        g, gamma, pprime, cm, &
        qlk, qrch, qs, &
        rain, rfact, shear, tfac, &
        val, val1, val2, &
        w1, w1l, w1s, w2, &
        w2l, w2s, w3, w3l, &
        w3s, w4, w4l, w4s, &
        rho, &
        xdby, xpw, xpwd, &
        xqrch, tem, tem1, tem2, &
        ptem, ptem1, ptem2

    integer :: kb (im), kbcon (im), kbcon1 (im), &
        ktcon (im), ktcon1 (im), ktconn (im), &
        jmin (im), lmin (im), kbmax (im), &
        kbm (im), kmax (im)

    ! real :: acrt (im), acrtfct (im),

    real :: aa1 (im), tkemean (im), clamt (im), &
        umean (im), tauadv (im), &
        delhbar (im), delq (im), delq2 (im), &
        delqbar (im), delqev (im), deltbar (im), &
        deltv (im), dtconv (im), edt (im), &
        edto (im), edtx (im), fld (im), &
        hcdo (im, km), hmax (im), hmin (im), &
        ucdo (im, km), vcdo (im, km), aa2 (im), &
        ecdo (im, km, ntr), &
        pdot (im), po (im, km), &
        pwavo (im), pwevo (im), mbdt (im), &
        qcdo (im, km), qcond (im), qevap (im), &
        rntot (im), vshear (im), xaa0 (im), &
        xk (im), xlamd (im), cina (im), &
        xmb (im), xmbmax (im), xpwav (im), &
        xpwev (im), delebar (im, ntr), &
        ! xlamx (im), &
        delubar (im), delvbar (im), &
        xlamdet (im), xlamddt (im), &
        cxlamet (im), cxlamdt (im)

    real :: c0 (im)

    real :: cinpcr, cinpcrmx, cinpcrmn, &
        cinacr, cinacrmx, cinacrmn

    ! parameters for updraft velocity calculation
    real :: bet1, cd1, f1, gam1, &
        bb1, bb2, wucb, csmf, tkcrt, cmxfac

    ! physical parameters
    parameter (g = grav)
    parameter (elocp = hlv / cp_air, el2orc = hlv * hlv / (rvgas * cp_air))
    parameter (d0 = .001)

    ! asolfac_deep: aerosol - aware parameter based on lim & hong (2012)
    ! asolfac_deep = cx / c0s_deep (= .002)
    ! cx = min ([ - 0.7 ln (nccn) + 24] * 1.e-4, c0s_deep)
    ! nccn: ccn number concentration in cm^ (- 3)
    ! until a realistic nccn is provided, typical nccns are assumed
    ! as nccn = 100 for sea and nccn = 1000 for land

    parameter (cm = 1.0, delta = zvir)
    parameter (fact1 = (cp_vap - c_liq) / rvgas, fact2 = hlv / rvgas - fact1 * tice)
    parameter (clamd = 0.03, tkemx = 0.65, tkemn = 0.05)
    parameter (dtke = tkemx - tkemn)
    parameter (dbeta = 0.1)
    parameter (dthk = 25.)
    parameter (cinpcrmx = 180., cinpcrmn = 120.)
    parameter (cinacrmx = - 120., cinacrmn = - 80.)
    parameter (bet1 = 1.875, cd1 = .506, f1 = 2.0, gam1 = .5)
    parameter (dxcrtuf = 15.e3)
    parameter (bb1 = 4.0, bb2 = 0.8, csmf = 0.2)
    parameter (tkcrt = 2., cmxfac = 15.)

    ! local variables and arrays
    real :: pfld (im, km), to (im, km), qo (im, km), &
        uo (im, km), vo (im, km), qeso (im, km), &
        ctr (im, km, ntr), ctro (im, km, ntr)

    ! for aerosol transport
    real :: qaero (im, km, ntc)

    ! for updraft velocity calculation
    real :: wu2 (im, km), buo (im, km), drag (im, km), wush (im, km)
    real :: wc (im), scaldfunc (im), sigmagfm (im)

    real :: qlko_ktcon (im), dellal (im, km), tvo (im, km), &
        dbyo (im, km), zo (im, km), &
        xlamue (im, km), xlamud (im, km), &
        fent1 (im, km), fent2 (im, km), frh (im, km), &
        heo (im, km), heso (im, km), &
        qrcd (im, km), dellah (im, km), dellaq (im, km), &
        dellae (im, km, ntr), &
        dellau (im, km), dellav (im, km), hcko (im, km), &
        ucko (im, km), vcko (im, km), qcko (im, km), &
        ecko (im, km, ntr), &
        eta (im, km), etad (im, km), zi (im, km), &
        qrcko (im, km), qrcdo (im, km), &
        pwo (im, km), pwdo (im, km), c0t (im, km), &
        tx1 (im), sumx (im), cnvwt (im, km)
    ! rhbar (im)

    logical :: do_aerosols, totflg, cnvflg (im), asqecflg (im), flg (im)

    ! asqecflg: flag for the quasi - equilibrium assumption of arakawa - schubert

    ! real :: pcrit (15), acritt (15), acrit (15)
    ! save pcrit, acritt
    ! data pcrit / 850., 800., 750., 700., 650., 600., 550., 500., 450., 400., &
    ! 350., 300., 250., 200., 150. /
    ! data acritt / .0633, .0445, .0553, .0664, .075, .1082, .1521, .2216, &
    ! .3151, .3677, .41, .5255, .7663, 1.1686, 1.6851 /
    ! gdas derived acrit
    ! data acritt / .203, .515, .521, .566, .625, .665, .659, .688, &
    ! .743, .813, .886, .947, 1.138, 1.377, 1.896 /

    real :: tf, tcr, tcrf
    parameter (tf = 233.16, tcr = 263.16, tcrf = 1.0 / (tcr - tf))

    ! -----------------------------------------------------------------------
    ! determine whether to perform aerosol transport
    ! -----------------------------------------------------------------------

    do_aerosols = (itc > 0) .and. (ntc > 0) .and. (ntr > 0)
    if (do_aerosols) do_aerosols = (ntr >= itc + ntc - 3)

    ! -----------------------------------------------------------------------
    ! compute preliminary quantities needed for static, dynamic, and feedback control portions of the algorithm.
    ! convert input pressure terms to centibar units.
    ! convert input pa terms to cb terms -- moorthi
    ! -----------------------------------------------------------------------

    km1 = km - 1

    ! -----------------------------------------------------------------------
    ! initialize column - integrated and other single - value - per - column variable arrays.
    ! initialize arrays
    ! -----------------------------------------------------------------------

    do i = 1, im
        cnvflg (i) = .true.
        rn (i) = 0.
        mbdt (i) = 10.
        kbot (i) = km + 1
        ktop (i) = 0
        kbcon (i) = km
        ktcon (i) = 1
        ktconn (i) = 1
        dtconv (i) = 3600.
        pdot (i) = 0.
        lmin (i) = 1
        jmin (i) = 1
        qlko_ktcon (i) = 0.
        edt (i) = 0.
        edto (i) = 0.
        edtx (i) = 0.
        ! acrt (i) = 0.
        ! acrtfct (i) = 1.
        aa1 (i) = 0.
        aa2 (i) = 0.
        xaa0 (i) = 0.
        cina (i) = 0.
        pwavo (i) = 0.
        pwevo (i) = 0.
        xpwav (i) = 0.
        xpwev (i) = 0.
        vshear (i) = 0.
    enddo

    ! -----------------------------------------------------------------------
    ! determine aerosol - aware rain conversion parameter over land
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (islimsk (i) == 1) then
            c0 (i) = c0s_deep * asolfac_deep
        else
            c0 (i) = c0s_deep
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! determine rain conversion parameter above the freezing level which exponentially
    ! decreases with decreasing temperature from han et al.'s (2017) equation 8.
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (t1 (i, k) > 273.16) then
                c0t (i, k) = c0 (i)
            else
                tem = d0 * (t1 (i, k) - 273.16)
                tem1 = exp (tem)
                c0t (i, k) = c0 (i) * tem1
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! initialize convective cloud water and cloud cover to zero.
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (present (cnvw)) cnvw (i, k) = 0.
            if (present (cnvc)) cnvc (i, k) = 0.
        enddo
    enddo

    do k = 1, km
        do i = 1, im
            qr (i, k) = 0.
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! hchuang code change
    ! initialize updraft and downdraft mass fluxes to zero.
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (present (ud_mf)) ud_mf (i, k) = 0.
            if (present (dd_mf)) dd_mf (i, k) = 0.
            if (present (dt_mf)) dt_mf (i, k) = 0.
        enddo
    enddo

    ! do k = 1, 15
    ! acrit (k) = acritt (k) * (975. - pcrit (k))
    ! enddo

    dt2 = delt
    ! val = 1200.
    val = 600.
    dtmin = max (dt2, val)
    ! val = 5400.
    val = 10800.
    dtmax = max (dt2, val)

    ! -----------------------------------------------------------------------
    ! model tunable parameters are all here
    ! -----------------------------------------------------------------------

    edtmaxl = .3
    edtmaxs = .3
    ! clam_deep = .1
    aafac = .05
    ! betal_deep = .15
    ! betas_deep = .15
    ! betal_deep = .05
    ! betas_deep = .05
    ! evef = 0.07
    ! evfact_deep = 0.3
    ! evfactl_deep = 0.3

    crtlame = 1.0e-4
    crtlamd = 1.0e-4

    ! cxlame = 1.0e-3
    cxlame = 1.0e-4
    cxlamd = 1.0e-4
    xlamde = 1.0e-4
    xlamdd = 1.0e-4

    ! pgcon_deep = 0.7 ! gregory et al. (1997, qjrms)
    ! pgcon_deep = 0.55 ! zhang & wu (2003, jas)

    w1l = - 8.e-3
    w2l = - 4.e-2
    w3l = - 5.e-3
    w4l = - 5.e-4
    w1s = - 2.e-4
    w2s = - 2.e-3
    w3s = - 1.e-3
    w4s = - 2.e-5

    ! -----------------------------------------------------------------------
    ! define top layer for search of the downdraft originating layer
    ! and the maximum thetae for updraft
    ! determine maximum indices for the parcel starting point (kbm), lfc (kbmax), and cloud top (kmax) .
    ! -----------------------------------------------------------------------

    do i = 1, im
        kbmax (i) = km
        kbm (i) = km
        kmax (i) = km
        tx1 (i) = 1.0 / psp (i)
    enddo

    do k = 1, km
        do i = 1, im
            if (prslp (i, k) * tx1 (i) > 0.04) kmax (i) = k + 1
            if (prslp (i, k) * tx1 (i) > 0.45) kbmax (i) = k + 1
            if (prslp (i, k) * tx1 (i) > 0.70) kbm (i) = k + 1
        enddo
    enddo

    do i = 1, im
        kmax (i) = min (km, kmax (i))
        kbmax (i) = min (kbmax (i), kmax (i))
        kbm (i) = min (kbm (i), kmax (i))
    enddo

    ! -----------------------------------------------------------------------
    ! hydrostatic height assume zero terr and initially assume
    ! updraft entrainment rate as an inverse function of height
    ! calculate hydrostatic height at layer centers assuming a flat surface (no terrain) from the geopotential.
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            zo (i, k) = phil (i, k) / g
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! calculate interface height
    ! -----------------------------------------------------------------------

    do k = 1, km1
        do i = 1, im
            zi (i, k) = 0.5 * (zo (i, k) + zo (i, k + 1))
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! convert surface pressure to mb from cb
    ! convert prslp from centibar to millibar, set normalized mass fluxes to 1, cloud properties to 0, and save model state variables (after advection / turbulence) .
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (k <= kmax (i)) then
                pfld (i, k) = prslp (i, k) * 0.01
                eta (i, k) = 1.
                fent1 (i, k) = 1.
                fent2 (i, k) = 1.
                frh (i, k) = 0.
                hcko (i, k) = 0.
                qcko (i, k) = 0.
                qrcko (i, k) = 0.
                ucko (i, k) = 0.
                vcko (i, k) = 0.
                etad (i, k) = 1.
                hcdo (i, k) = 0.
                qcdo (i, k) = 0.
                ucdo (i, k) = 0.
                vcdo (i, k) = 0.
                qrcd (i, k) = 0.
                qrcdo (i, k) = 0.
                dbyo (i, k) = 0.
                pwo (i, k) = 0.
                pwdo (i, k) = 0.
                dellal (i, k) = 0.
                to (i, k) = t1 (i, k)
                qo (i, k) = q1 (i, k)
                uo (i, k) = u1 (i, k)
                vo (i, k) = v1 (i, k)
                ! uo (i, k) = u1 (i, k) * rcs (i)
                ! vo (i, k) = v1 (i, k) * rcs (i)
                wu2 (i, k) = 0.
                buo (i, k) = 0.
                drag (i, k) = 0.
                cnvwt (i, k) = 0.
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! initialize tracer variables
    ! -----------------------------------------------------------------------

    kk = 0
    do n = 1, ntr + 2
        if (n .eq. ntw .or. n .eq. nti) cycle
        kk = kk + 1
        do k = 1, km
            do i = 1, im
                if (k <= kmax (i)) then
                    ctr (i, k, kk) = qtr (i, k, n)
                    ctro (i, k, kk) = qtr (i, k, n)
                    ecko (i, k, kk) = 0.
                    ecdo (i, k, kk) = 0.
                endif
            enddo
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! column variables
    ! p is pressure of the layer (mb)
    ! t is temperature at t - dt (k) ..tn
    ! q is mixing ratio at t - dt (kg / kg) ..qn
    ! to is temperature at t + dt (k) ... this is after advection and turbulan
    ! qo is mixing ratio at t + dt (kg / kg) ..q1
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! calculate saturation specific humidity and enforce minimum moisture values.
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (k <= kmax (i)) then
                qeso (i, k) = 0.01 * mqs (to (i, k)) ! mqs is in pa
                qeso (i, k) = eps * qeso (i, k) / (pfld (i, k) + epsm1 * qeso (i, k))
                val1 = 1.e-8
                qeso (i, k) = max (qeso (i, k), val1)
                val2 = 1.e-10
                qo (i, k) = max (qo (i, k), val2)
                ! qo (i, k) = min (qo (i, k), qeso (i, k))
                ! tvo (i, k) = to (i, k) + delta * to (i, k) * qo (i, k)
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute moist static energy
    ! calculate moist static energy (heo) and saturation moist static energy (heso) .
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (k <= kmax (i)) then
                ! tem = g * zo (i, k) + cp_air * to (i, k)
                tem = phil (i, k) + cp_air * to (i, k)
                heo (i, k) = tem + hlv * qo (i, k)
                heso (i, k) = tem + hlv * qeso (i, k)
                ! heo (i, k) = min (heo (i, k), heso (i, k))
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! determine level with largest moist static energy
    ! this is the level where updraft starts
    ! perform calculations related to the updraft of the entraining / detraining cloud model ("static control") .
    ! search below index "kbm" for the level of maximum moist static energy.
    ! -----------------------------------------------------------------------

    do i = 1, im
        hmax (i) = heo (i, 1)
        kb (i) = 1
    enddo

    do k = 2, km
        do i = 1, im
            if (k <= kbm (i)) then
                if (heo (i, k) > hmax (i)) then
                    kb (i) = k
                    hmax (i) = heo (i, k)
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! calculate the temperature, specific humidity, and pressure at interface levels.
    ! -----------------------------------------------------------------------

    do k = 1, km1
        do i = 1, im
            if (k <= kmax (i) - 1) then
                dz = .5 * (zo (i, k + 1) - zo (i, k))
                dp = .5 * (pfld (i, k + 1) - pfld (i, k))
                es = 0.01 * mqs (to (i, k + 1)) ! mqs is in pa
                pprime = pfld (i, k + 1) + epsm1 * es
                qs = eps * es / pprime
                dqsdp = - qs / pprime
                desdt = es * (fact1 / to (i, k + 1) + fact2 / (to (i, k + 1) ** 2))
                dqsdt = qs * pfld (i, k + 1) * desdt / (es * pprime)
                gamma = el2orc * qeso (i, k + 1) / (to (i, k + 1) ** 2)
                dt = (g * dz + hlv * dqsdp * dp) / (cp_air * (1. + gamma))
                dq = dqsdt * dt + dqsdp * dp
                to (i, k) = to (i, k + 1) + dt
                qo (i, k) = qo (i, k + 1) + dq
                po (i, k) = .5 * (pfld (i, k) + pfld (i, k + 1))
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! recalculate saturation specific humidity, moist static energy, saturation moist static energy, and horizontal momentum on interface levels. enforce minimum specific humidity and calculate \f$ (1 - rh) \f$.
    ! -----------------------------------------------------------------------

    do k = 1, km1
        do i = 1, im
            if (k <= kmax (i) - 1) then
                qeso (i, k) = 0.01 * mqs (to (i, k)) ! mqs is in pa
                qeso (i, k) = eps * qeso (i, k) / (po (i, k) + epsm1 * qeso (i, k))
                val1 = 1.e-8
                qeso (i, k) = max (qeso (i, k), val1)
                val2 = 1.e-10
                qo (i, k) = max (qo (i, k), val2)
                ! qo (i, k) = min (qo (i, k), qeso (i, k))
                frh (i, k) = 1. - min (qo (i, k) / qeso (i, k), 1.)
                heo (i, k) = .5 * g * (zo (i, k) + zo (i, k + 1)) + &
                    cp_air * to (i, k) + hlv * qo (i, k)
                heso (i, k) = .5 * g * (zo (i, k) + zo (i, k + 1)) + &
                    cp_air * to (i, k) + hlv * qeso (i, k)
                uo (i, k) = .5 * (uo (i, k) + uo (i, k + 1))
                vo (i, k) = .5 * (vo (i, k) + vo (i, k + 1))
            endif
        enddo
    enddo

    do n = 1, ntr
        do k = 1, km1
            do i = 1, im
                if (k <= kmax (i) - 1) then
                    ctro (i, k, n) = .5 * (ctro (i, k, n) + ctro (i, k + 1, n))
                endif
            enddo
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! look for the level of free convection as cloud base
    ! search below the index "kbmax" for the level of free convection (lfc) where the condition \f$h_b > h^ * \f$ is first met, where \f$h_b, h^ * \f$ are the state moist static energy at the parcel's starting level and saturation moist static energy, respectively. set "kbcon" to the index of the lfc.
    ! -----------------------------------------------------------------------

    do i = 1, im
        flg (i) = .true.
        kbcon (i) = kmax (i)
    enddo

    do k = 1, km1
        do i = 1, im
            if (flg (i) .and. k <= kbmax (i)) then
                if (k > kb (i) .and. heo (i, kb (i)) > heso (i, k)) then
                    kbcon (i) = k
                    flg (i) = .false.
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! if no lfc, return to the calling routine without modifying state variables.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (kbcon (i) == kmax (i)) cnvflg (i) = .false.
    enddo

    totflg = .true.
    do i = 1, im
        totflg = totflg .and. (.not. cnvflg (i))
    enddo

    if (totflg) return

    ! -----------------------------------------------------------------------
    ! determine the vertical pressure velocity at the lfc. after han and pan (2011), determine the maximum pressure thickness between a parcel's starting level and the lfc. if a parcel doesn't reach the lfc within the critical thickness, then the convective inhibition is deemed too great for convection to be triggered, and the subroutine returns to the calling routine without modifying the state variables.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            ! pdot (i) = 10. * dot (i, kbcon (i))
            pdot (i) = 0.01 * dot (i, kbcon (i)) ! now dot is in pa / s
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! turn off convection if pressure depth between parcel source level
    ! and cloud base is larger than a critical value, cinpcr
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            if (islimsk (i) == 1) then
                w1 = w1l
                w2 = w2l
                w3 = w3l
                w4 = w4l
            else
                w1 = w1s
                w2 = w2s
                w3 = w3s
                w4 = w4s
            endif
            if (pdot (i) <= w4) then
                tem = (pdot (i) - w4) / (w3 - w4)
            elseif (pdot (i) >= - w4) then
                tem = - (pdot (i) + w4) / (w4 - w3)
            else
                tem = 0.
            endif
            val1 = - 1.
            tem = max (tem, val1)
            val2 = 1.
            tem = min (tem, val2)
            ptem = 1. - tem
            ptem1 = .5 * (cinpcrmx - cinpcrmn)
            cinpcr = cinpcrmx - ptem * ptem1
            tem1 = pfld (i, kb (i)) - pfld (i, kbcon (i))
            if (tem1 > cinpcr) then
                cnvflg (i) = .false.
            endif
        endif
    enddo

    totflg = .true.
    do i = 1, im
        totflg = totflg .and. (.not. cnvflg (i))
    enddo

    if (totflg) return

    ! -----------------------------------------------------------------------
    ! turbulent entrainment rate assumed to be proportional
    ! to subcloud mean tke
    ! -----------------------------------------------------------------------

    if (ntk > 0) then

        do i = 1, im
            if (cnvflg (i)) then
                sumx (i) = 0.
                tkemean (i) = 0.
            endif
        enddo
        do k = 1, km1
            do i = 1, im
                if (cnvflg (i)) then
                    if (k >= kb (i) .and. k < kbcon (i)) then
                        dz = zo (i, k + 1) - zo (i, k)
                        tem = 0.5 * (qtr (i, k, ntk) + qtr (i, k + 1, ntk))
                        tkemean (i) = tkemean (i) + tem * dz
                        sumx (i) = sumx (i) + dz
                    endif
                endif
            enddo
        enddo
        do i = 1, im
            if (cnvflg (i)) then
                tkemean (i) = tkemean (i) / sumx (i)
                if (tkemean (i) > tkemx) then
                    clamt (i) = clam_deep + clamd
                else if (tkemean (i) < tkemn) then
                    clamt (i) = clam_deep - clamd
                else
                    tem = tkemx - tkemean (i)
                    tem1 = 1. - 2. * tem / dtke
                    clamt (i) = clam_deep + clamd * tem1
                endif
            endif
        enddo
        ! kgao 12 / 08 / 2023: adjust ent / det rates based on tke
        if (use_tke_conv) then
            do i = 1, im
                if (cnvflg (i)) then
                    xlamdet (i) = xlamde
                    xlamddt (i) = xlamdd
                    cxlamet (i) = cxlame
                    cxlamdt (i) = cxlamd
                    if (tkemean (i) > tkcrt) then
                        tem = 1. + tkemean (i) / tkcrt
                        tem1 = min (tem, cmxfac)
                        clamt (i) = tem1 * clam_deep
                        xlamdet (i) = tem1 * xlamdet (i)
                        xlamddt (i) = tem1 * xlamddt (i)
                        cxlamet (i) = tem1 * cxlamet (i)
                        cxlamdt (i) = tem1 * cxlamdt (i)
                    endif
                endif
            enddo
        endif
    else
        do i = 1, im
            if (cnvflg (i)) then
                clamt (i) = clam_deep
            endif
        enddo
    endif

    ! -----------------------------------------------------------------------
    ! also initially assume updraft entrainment rate
    ! is an inverse function of height
    ! -----------------------------------------------------------------------

    do k = 1, km1
        do i = 1, im
            if (cnvflg (i)) then
                xlamue (i, k) = clamt (i) / zi (i, k)
                xlamue (i, k) = max (xlamue (i, k), crtlame)
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! assume that updraft entrainment rate above cloud base is
    ! same as that at cloud base
    ! -----------------------------------------------------------------------

    ! calculate the entrainment rate according to han and pan (2011), equation 8, after bechtold et al. (2008), equation 2 given by:
    ! \f[
    ! \epsilon = \epsilon_0f_0 + d_1\left (1 - rh\right) f_1
    ! \f]
    ! where \f$\epsilon_0\f$ is the cloud base entrainment rate, \f$d_1\f$ is a tunable constant, and \f$f_0 = \left (\frac{q_s}{q_{s, b}}\right) ^2\f$ and \f$f_1 = \left (\frac{q_s}{q_{s, b}}\right) ^3\f$ where \f$q_s\f$ and \f$q_{s, b}\f$ are the saturation specific humidities at a given level and cloud base, respectively. the detrainment rate in the cloud is assumed to be equal to the entrainment rate at cloud base.
    ! do i = 1, im
    ! if (cnvflg (i)) then
    ! xlamx (i) = xlamue (i, kbcon (i))
    ! endif
    ! enddo
    ! do k = 2, km1
    ! do i = 1, im
    ! if (cnvflg (i) .and. (k > kbcon (i) .and. k < kmax (i))) then
    ! xlamue (i, k) = xlamx (i)
    ! endif
    ! enddo
    ! enddo

    ! -----------------------------------------------------------------------
    ! specify detrainment rate for the updrafts
    ! (the updraft detrainment rate is set constant and equal to the entrainment rate at cloud base.)
    ! the updraft detrainment rate is vertically constant and proportional to clamt
    ! -----------------------------------------------------------------------

    do k = 1, km1
        do i = 1, im
            if (cnvflg (i) .and. k < kmax (i)) then
                ! xlamud (i, k) = xlamx (i)
                ! xlamud (i, k) = crtlamd
                xlamud (i, k) = 0.001 * clamt (i)
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! functions rapidly decreasing with height, mimicking a cloud ensemble
    ! entrainment functions decreasing with height (fent),
    ! mimicking a cloud ensemble
    ! (bechtold et al., 2008)
    ! -----------------------------------------------------------------------

    do k = 2, km1
        do i = 1, im
            if (cnvflg (i) .and. (k > kbcon (i) .and. k < kmax (i))) then
                tem = qeso (i, k) / qeso (i, kbcon (i))
                fent1 (i, k) = tem ** 2
                fent2 (i, k) = tem ** 3
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! final entrainment and detrainment rates as the sum of turbulent part and
    ! organized one depending on the environmental relative humidity
    ! (bechtold et al., 2008; derbyshire et al., 2011)
    ! -----------------------------------------------------------------------

    ! kgao 12 / 21 / 2023
    if (use_tke_conv) then
        ! new code
        do k = 2, km1
            do i = 1, im
                if (cnvflg (i) .and. (k > kbcon (i) .and. k < kmax (i))) then
                    tem = cxlamet (i) * frh (i, k) * fent2 (i, k)
                    xlamue (i, k) = xlamue (i, k) * fent1 (i, k) + tem
                    tem1 = cxlamdt (i) * frh (i, k)
                    xlamud (i, k) = xlamud (i, k) + tem1
                endif
            enddo
        enddo
    else
        ! ori code
        do k = 2, km1
            do i = 1, im
                if (cnvflg (i) .and. (k > kbcon (i) .and. k < kmax (i))) then
                    tem = cxlame * frh (i, k) * fent2 (i, k)
                    xlamue (i, k) = xlamue (i, k) * fent1 (i, k) + tem
                    tem1 = cxlamd * frh (i, k)
                    xlamud (i, k) = xlamud (i, k) + tem1
                endif
            enddo
        enddo
    endif ! end of use_tke_conv

    ! -----------------------------------------------------------------------
    ! determine updraft mass flux for the subcloud layers
    ! calculate the normalized mass flux for subcloud and in - cloud layers according to pan and wu (1995) equation 1:
    ! \f[
    ! \frac{1}{\eta}\frac{\partial \eta}{\partial z} = \lambda_e - \lambda_d
    ! \f]
    ! where \f$\eta\f$ is the normalized mass flux, \f$\lambda_e\f$ is the entrainment rate and \f$\lambda_d\f$ is the detrainment rate.
    ! -----------------------------------------------------------------------

    do k = km1, 1, - 1
        do i = 1, im
            if (cnvflg (i)) then
                if (k < kbcon (i) .and. k >= kb (i)) then
                    dz = zi (i, k + 1) - zi (i, k)
                    tem = 0.5 * (xlamud (i, k) + xlamud (i, k + 1))
                    ptem = 0.5 * (xlamue (i, k) + xlamue (i, k + 1)) - tem
                    eta (i, k) = eta (i, k + 1) / (1. + ptem * dz)
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute mass flux above cloud base
    ! -----------------------------------------------------------------------

    do i = 1, im
        flg (i) = cnvflg (i)
    enddo

    do k = 2, km1
        do i = 1, im
            if (flg (i)) then
                if (k > kbcon (i) .and. k < kmax (i)) then
                    dz = zi (i, k) - zi (i, k - 1)
                    tem = 0.5 * (xlamud (i, k) + xlamud (i, k - 1))
                    ptem = 0.5 * (xlamue (i, k) + xlamue (i, k - 1)) - tem
                    eta (i, k) = eta (i, k - 1) * (1 + ptem * dz)
                    if (eta (i, k) <= 0.) then
                        kmax (i) = k
                        ktconn (i) = k
                        flg (i) = .false.
                    endif
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute updraft cloud properties
    ! set cloud properties equal to the state variables at updraft starting level (kb) .
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            indx = kb (i)
            hcko (i, indx) = heo (i, indx)
            ucko (i, indx) = uo (i, indx)
            vcko (i, indx) = vo (i, indx)
            pwavo (i) = 0.
        endif
    enddo

    ! for tracers
    do n = 1, ntr
        do i = 1, im
            if (cnvflg (i)) then
                indx = kb (i)
                ecko (i, indx, n) = ctro (i, indx, n)
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! cloud property is modified by the entrainment process
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! cloud property is modified by the entrainment process
    ! cm is an enhancement factor in entrainment rates for momentum
    ! calculate the cloud properties as a parcel ascends, modified by entrainment and detrainment. discretization follows appendix b of grell (1993). following han and pan (2006), the convective momentum transport is reduced by the convection - induced pressure gradient force by the constant "pgcon_deep", currently set to 0.55 after zhang and wu (2003).
    ! -----------------------------------------------------------------------

    do k = 2, km1
        do i = 1, im
            if (cnvflg (i)) then
                if (k > kb (i) .and. k < kmax (i)) then
                    dz = zi (i, k) - zi (i, k - 1)
                    tem = 0.5 * (xlamue (i, k) + xlamue (i, k - 1)) * dz
                    tem1 = 0.25 * (xlamud (i, k) + xlamud (i, k - 1)) * dz
                    factor = 1. + tem - tem1
                    hcko (i, k) = ((1. - tem1) * hcko (i, k - 1) + tem * 0.5 * &
                         (heo (i, k) + heo (i, k - 1))) / factor
                    dbyo (i, k) = hcko (i, k) - heso (i, k)

                    tem = 0.5 * cm * tem
                    factor = 1. + tem
                    ptem = tem + pgcon_deep
                    ptem1 = tem - pgcon_deep
                    ucko (i, k) = ((1. - tem) * ucko (i, k - 1) + ptem * uo (i, k) &
                         + ptem1 * uo (i, k - 1)) / factor
                    vcko (i, k) = ((1. - tem) * vcko (i, k - 1) + ptem * vo (i, k) &
                         + ptem1 * vo (i, k - 1)) / factor
                endif
            endif
        enddo
    enddo

    do n = 1, ntr
        do k = 2, km1
            do i = 1, im
                if (cnvflg (i)) then
                    if (k > kb (i) .and. k < kmax (i)) then
                        dz = zi (i, k) - zi (i, k - 1)
                        tem = 0.25 * (xlamue (i, k) + xlamue (i, k - 1)) * dz
                        factor = 1. + tem
                        ecko (i, k, n) = ((1. - tem) * ecko (i, k - 1, n) + tem * &
                            (ctro (i, k, n) + ctro (i, k - 1, n))) / factor
                    endif
                endif
            enddo
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! taking account into convection inhibition due to existence of
    ! dry layers below cloud base
    ! with entrainment, recalculate the lfc as the first level where buoyancy is positive. the difference in pressure levels between lfcs calculated with / without entrainment must be less than a threshold (currently 25 hpa) . otherwise, convection is inhibited and the scheme returns to the calling routine without modifying the state variables. this is the subcloud dryness trigger modification discussed in han and pan (2011).
    ! -----------------------------------------------------------------------

    do i = 1, im
        flg (i) = cnvflg (i)
        kbcon1 (i) = kmax (i)
    enddo

    do k = 2, km1
        do i = 1, im
            if (flg (i) .and. k < kmax (i)) then
                if (k >= kbcon (i) .and. dbyo (i, k) > 0.) then
                    kbcon1 (i) = k
                    flg (i) = .false.
                endif
            endif
        enddo
    enddo

    do i = 1, im
        if (cnvflg (i)) then
            if (kbcon1 (i) == kmax (i)) cnvflg (i) = .false.
        endif
    enddo

    do i = 1, im
        if (cnvflg (i)) then
            tem = pfld (i, kbcon (i)) - pfld (i, kbcon1 (i))
            if (tem > dthk) then
                cnvflg (i) = .false.
            endif
        endif
    enddo

    totflg = .true.
    do i = 1, im
        totflg = totflg .and. (.not. cnvflg (i))
    enddo

    if (totflg) return

    ! -----------------------------------------------------------------------
    ! calculate convective inhibition
    ! calculate additional trigger condition of the convective inhibition (cin) according to han et al.'s (2017) equation 13.
    ! -----------------------------------------------------------------------

    do k = 2, km1
        do i = 1, im
            if (cnvflg (i)) then
                if (k > kb (i) .and. k < kbcon1 (i)) then
                    dz1 = zo (i, k + 1) - zo (i, k)
                    gamma = el2orc * qeso (i, k) / (to (i, k) ** 2)
                    rfact = 1. + delta * cp_air * gamma &
                         * to (i, k) / hlv
                    cina (i) = cina (i) + &
                        ! dz1 * eta (i, k) * (g / (cp_air * to (i, k))) &
                        dz1 * (g / (cp_air * to (i, k))) &
                         * dbyo (i, k) / (1. + gamma) &
                         * rfact
                    val = 0.
                    cina (i) = cina (i) + &
                        ! dz1 * eta (i, k) * g * delta * &
                        dz1 * g * delta * &
                        max (val, (qeso (i, k) - qo (i, k)))
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! turn off convection if the cin is less than a critical value (cinacr) which is inversely proportional to the large - scale vertical velocity.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then

            if (islimsk (i) == 1) then
                w1 = w1l
                w2 = w2l
                w3 = w3l
                w4 = w4l
            else
                w1 = w1s
                w2 = w2s
                w3 = w3s
                w4 = w4s
            endif
            if (pdot (i) <= w4) then
                tem = (pdot (i) - w4) / (w3 - w4)
            elseif (pdot (i) >= - w4) then
                tem = - (pdot (i) + w4) / (w4 - w3)
            else
                tem = 0.
            endif

            val1 = - 1.
            tem = max (tem, val1)
            val2 = 1.
            tem = min (tem, val2)
            tem = 1. - tem
            tem1 = .5 * (cinacrmx - cinacrmn)
            cinacr = cinacrmx - tem * tem1

            ! cinacr = cinacrmx
            if (cina (i) < cinacr) cnvflg (i) = .false.
        endif
    enddo

    totflg = .true.
    do i = 1, im
        totflg = totflg .and. (.not. cnvflg (i))
    enddo

    if (totflg) return

    ! -----------------------------------------------------------------------
    ! determine first guess cloud top as the level of zero buoyancy
    ! calculate the cloud top as the first level where parcel buoyancy becomes negative. if the thickness of the calculated convection is less than a threshold (currently 200 hpa), then convection is inhibited, and the scheme returns to the calling routine.
    ! -----------------------------------------------------------------------

    do i = 1, im
        flg (i) = cnvflg (i)
        ktcon (i) = 1
    enddo

    do k = 2, km1
        do i = 1, im
            if (flg (i) .and. k < kmax (i)) then
                if (k > kbcon1 (i) .and. dbyo (i, k) < 0.) then
                    ktcon (i) = k
                    flg (i) = .false.
                endif
            endif
        enddo
    enddo

    do i = 1, im
        if (cnvflg (i)) then
            if (ktcon (i) == 1 .and. ktconn (i) > 1) then
                ktcon (i) = ktconn (i)
            endif
            tem = pfld (i, kbcon (i)) - pfld (i, ktcon (i))
            if (tem < cthk_deep) cnvflg (i) = .false.
        endif
    enddo

    totflg = .true.
    do i = 1, im
        totflg = totflg .and. (.not. cnvflg (i))
    enddo

    if (totflg) return

    ! -----------------------------------------------------------------------
    ! search for downdraft originating level above theta - e minimum
    ! to originate the downdraft, search for the level above the minimum in moist static energy. return to the calling routine without modification if this level is determined to be outside of the convective cloud layers.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            hmin (i) = heo (i, kbcon1 (i))
            lmin (i) = kbmax (i)
            jmin (i) = kbmax (i)
        endif
    enddo

    do k = 2, km1
        do i = 1, im
            if (cnvflg (i) .and. k <= kbmax (i)) then
                if (k > kbcon1 (i) .and. heo (i, k) < hmin (i)) then
                    lmin (i) = k + 1
                    hmin (i) = heo (i, k)
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! make sure that jmin (i) is within the cloud
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            jmin (i) = min (lmin (i), ktcon (i) - 1)
            jmin (i) = max (jmin (i), kbcon1 (i) + 1)
            if (jmin (i) >= ktcon (i)) cnvflg (i) = .false.
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! specify upper limit of mass flux at cloud base
    ! calculate the maximum value of the cloud base mass flux using the cfl - criterion - based formula of han and pan (2011), equation 7.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            ! xmbmax (i) = .1

            k = kbcon (i)
            dp = delp (i, k)
            xmbmax (i) = dp / (2. * g * dt2)
            ! xmbmax (i) = dp / (g * dt2)

            ! mbdt (i) = 0.1 * dp / g

            ! tem = dp / (g * dt2)
            ! xmbmax (i) = min (tem, xmbmax (i))
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! compute cloud moisture property and precipitation
    ! set cloud moisture property equal to the enviromental moisture at updraft starting level (kb) .
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            ! aa1 (i) = 0.
            qcko (i, kb (i)) = qo (i, kb (i))
            qrcko (i, kb (i)) = qo (i, kb (i))
            ! rhbar (i) = 0.
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! calculate the moisture content of the entraining / detraining parcel (qcko) and the value it would have if just saturated (qrch), according to equation a.14 in grell (1993). their difference is the amount of convective cloud water (qlk = rain + condensate) . determine the portion of convective cloud water that remains suspended and the portion that is converted into convective precipitation (pwo) . calculate and save the negative cloud work function (aa1) due to water loading. the liquid water in the updraft layer is assumed to be detrained from the layers above the level of the minimum moist static energy into the grid - scale cloud water (dellal) .
    ! -----------------------------------------------------------------------

    do k = 2, km1
        do i = 1, im
            if (cnvflg (i)) then
                if (k > kb (i) .and. k < ktcon (i)) then

                    dz = zi (i, k) - zi (i, k - 1)
                    gamma = el2orc * qeso (i, k) / (to (i, k) ** 2)
                    qrch = qeso (i, k) &
                         + gamma * dbyo (i, k) / (hlv * (1. + gamma))
                    tem = 0.5 * (xlamue (i, k) + xlamue (i, k - 1)) * dz
                    tem1 = 0.25 * (xlamud (i, k) + xlamud (i, k - 1)) * dz
                    factor = 1. + tem - tem1
                    qcko (i, k) = ((1. - tem1) * qcko (i, k - 1) + tem * 0.5 * &
                         (qo (i, k) + qo (i, k - 1))) / factor
                    qrcko (i, k) = qcko (i, k)
                    dq = eta (i, k) * (qcko (i, k) - qrch)

                    ! rhbar (i) = rhbar (i) + qo (i, k) / qeso (i, k)

                    ! -----------------------------------------------------------------------
                    ! check if there is excess moisture to release latent heat
                    ! -----------------------------------------------------------------------

                    if (k >= kbcon (i) .and. dq > 0.) then
                        etah = .5 * (eta (i, k) + eta (i, k - 1))
                        dp = delp (i, k)
                        if (ncloud > 0 .and. k > jmin (i)) then
                            ptem = c0t (i, k) + c1_deep
                            qlk = dq / (eta (i, k) + etah * ptem * dz)
                            dellal (i, k) = etah * c1_deep * dz * qlk * g / dp
                        else
                            qlk = dq / (eta (i, k) + etah * c0t (i, k) * dz)
                        endif
                        ! aa1 (i) = aa1 (i) - dz * g * qlk * etah
                        ! aa1 (i) = aa1 (i) - dz * g * qlk
                        buo (i, k) = buo (i, k) - g * qlk
                        qcko (i, k) = qlk + qrch
                        pwo (i, k) = etah * c0t (i, k) * dz * qlk
                        pwavo (i) = pwavo (i) + pwo (i, k)
                        ! cnvwt (i, k) = (etah * qlk + pwo (i, k)) * g / dp
                        cnvwt (i, k) = etah * qlk * g / dp
                    endif

                    ! -----------------------------------------------------------------------
                    ! compute buoyancy and drag for updraft velocity
                    ! -----------------------------------------------------------------------

                    if (k >= kbcon (i)) then
                        rfact = 1. + delta * cp_air * gamma &
                             * to (i, k) / hlv
                        buo (i, k) = buo (i, k) + (g / (cp_air * to (i, k))) &
                             * dbyo (i, k) / (1. + gamma) &
                             * rfact
                        val = 0.
                        buo (i, k) = buo (i, k) + g * delta * &
                            max (val, (qeso (i, k) - qo (i, k)))
                        drag (i, k) = max (xlamue (i, k), xlamud (i, k))
                        ! kgao 12 / 18 / 2023: considers shear effect
                        tem = ((uo (i, k) - uo (i, k - 1)) / dz) ** 2
                        tem = tem + ((vo (i, k) - vo (i, k - 1)) / dz) ** 2
                        wush (i, k) = csmf * sqrt (tem)
                    endif

                endif
            endif
        enddo
    enddo

    ! do i = 1, im
    ! if (cnvflg (i)) then
    ! indx = ktcon (i) - kb (i) - 1
    ! rhbar (i) = rhbar (i) / float (indx)
    ! endif
    ! enddo

    ! -----------------------------------------------------------------------
    ! calculate cloud work function
    ! -----------------------------------------------------------------------

    ! do k = 2, km1
    ! do i = 1, im
    ! if (cnvflg (i)) then
    ! if (k >= kbcon (i) .and. k < ktcon (i)) then
    ! dz1 = zo (i, k + 1) - zo (i, k)
    ! gamma = el2orc * qeso (i, k) / (to (i, k) ** 2)
    ! rfact = 1. + delta * cp_air * gamma &
    ! * to (i, k) / hlv
    ! aa1 (i) = aa1 (i) + &
    ! dz1 * eta (i, k) * (g / (cp_air * to (i, k))) &
    ! dz1 * (g / (cp_air * to (i, k))) &
    ! * dbyo (i, k) / (1. + gamma) &
    ! * rfact
    ! val = 0.
    ! aa1 (i) = aa1 (i) + &
    ! dz1 * eta (i, k) * g * delta * &
    ! dz1 * g * delta * &
    ! max (val, (qeso (i, k) - qo (i, k)))
    ! endif
    ! endif
    ! enddo
    ! enddo

    ! -----------------------------------------------------------------------
    ! calculate cloud work function
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! calculate the cloud work function according to pan and wu (1995) equation 4:
    ! \f[
    ! a_u = \int_{z_0}^{z_t}\frac{g}{c_pt (z) }\frac{\eta}{1 + \gamma}[h (z) - h^ * (z) ]dz
    ! \f]
    ! (discretized according to grell (1993) equation b.10 using b.2 and b.3 of arakawa and schubert (1974) and assuming \f$\eta = 1\f$) where \f$a_u\f$ is the updraft cloud work function, \f$z_0\f$ and \f$z_t\f$ are cloud base and cloud top, respectively, \f$\gamma = \frac{l}{c_p}\left (\frac{\partial \overline{q_s}}{\partial t}\right) _p\f$ and other quantities are previously defined.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            aa1 (i) = 0.
        endif
    enddo

    do k = 2, km1
        do i = 1, im
            if (cnvflg (i)) then
                if (k >= kbcon (i) .and. k < ktcon (i)) then
                    dz1 = zo (i, k + 1) - zo (i, k)
                    ! aa1 (i) = aa1 (i) + buo (i, k) * dz1 * eta (i, k)
                    aa1 (i) = aa1 (i) + buo (i, k) * dz1
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! if the updraft cloud work function is negative, convection does not occur, and the scheme returns to the calling routine.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i) .and. aa1 (i) <= 0.) cnvflg (i) = .false.
    enddo

    totflg = .true.
    do i = 1, im
        totflg = totflg .and. (.not. cnvflg (i))
    enddo

    if (totflg) return

    ! -----------------------------------------------------------------------
    ! estimate the onvective overshooting as the level
    ! where the [aafac * cloud work function] becomes zero,
    ! which is the final cloud top
    ! continue calculating the cloud work function past the point of neutral buoyancy to represent overshooting according to han and pan (2011). convective overshooting stops when \f$ ca_u < 0\f$ where \f$c\f$ is currently 10%, or when 10% of the updraft cloud work function has been consumed by the stable buoyancy force.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            aa2 (i) = aafac * aa1 (i)
        endif
    enddo

    do i = 1, im
        flg (i) = cnvflg (i)
        ktcon1 (i) = kmax (i)
    enddo

    do k = 2, km1
        do i = 1, im
            if (flg (i)) then
                if (k >= ktcon (i) .and. k < kmax (i)) then
                    dz1 = zo (i, k + 1) - zo (i, k)
                    gamma = el2orc * qeso (i, k) / (to (i, k) ** 2)
                    rfact = 1. + delta * cp_air * gamma &
                         * to (i, k) / hlv
                    aa2 (i) = aa2 (i) + &
                        ! dz1 * eta (i, k) * (g / (cp_air * to (i, k))) &
                        dz1 * (g / (cp_air * to (i, k))) &
                         * dbyo (i, k) / (1. + gamma) &
                         * rfact
                    ! val = 0.
                    ! aa2 (i) = aa2 (i) + &
                    ! dz1 * eta (i, k) * g * delta * &
                    ! dz1 * g * delta * &
                    ! max (val, (qeso (i, k) - qo (i, k)))
                    if (aa2 (i) < 0.) then
                        ktcon1 (i) = k
                        flg (i) = .false.
                    endif
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute cloud moisture property, detraining cloud water
    ! and precipitation in overshooting layers
    ! for the overshooting convection, calculate the moisture content of the entraining / detraining parcel as before. partition convective cloud water and precipitation and detrain convective cloud water above the mimimum in moist static energy.
    ! -----------------------------------------------------------------------

    do k = 2, km1
        do i = 1, im
            if (cnvflg (i)) then
                if (k >= ktcon (i) .and. k < ktcon1 (i)) then
                    dz = zi (i, k) - zi (i, k - 1)
                    gamma = el2orc * qeso (i, k) / (to (i, k) ** 2)
                    qrch = qeso (i, k) &
                         + gamma * dbyo (i, k) / (hlv * (1. + gamma))
                    tem = 0.5 * (xlamue (i, k) + xlamue (i, k - 1)) * dz
                    tem1 = 0.25 * (xlamud (i, k) + xlamud (i, k - 1)) * dz
                    factor = 1. + tem - tem1
                    qcko (i, k) = ((1. - tem1) * qcko (i, k - 1) + tem * 0.5 * &
                         (qo (i, k) + qo (i, k - 1))) / factor
                    qrcko (i, k) = qcko (i, k)
                    dq = eta (i, k) * (qcko (i, k) - qrch)

                    ! -----------------------------------------------------------------------
                    ! check if there is excess moisture to release latent heat
                    ! -----------------------------------------------------------------------

                    if (dq > 0.) then
                        etah = .5 * (eta (i, k) + eta (i, k - 1))
                        dp = delp (i, k)
                        if (ncloud > 0) then
                            ptem = c0t (i, k) + c1_deep
                            qlk = dq / (eta (i, k) + etah * ptem * dz)
                            dellal (i, k) = etah * c1_deep * dz * qlk * g / dp
                        else
                            qlk = dq / (eta (i, k) + etah * c0t (i, k) * dz)
                        endif
                        qcko (i, k) = qlk + qrch
                        pwo (i, k) = etah * c0t (i, k) * dz * qlk
                        pwavo (i) = pwavo (i) + pwo (i, k)
                        ! cnvwt (i, k) = (etah * qlk + pwo (i, k)) * g / dp
                        cnvwt (i, k) = etah * qlk * g / dp
                    endif
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute updraft velocity square (wu2)
    ! calculate updraft velocity square (wu2) according to han et al.'s (2017) equation 7.
    ! -----------------------------------------------------------------------

    ! bb1 = 2. * (1. + bet1 * cd1)
    ! bb2 = 2. / (f1 * (1. + gam1))
    !
    ! bb1 = 3.9
    ! bb2 = 0.67
    !
    ! bb1 = 2.0
    ! bb2 = 4.0

    ! bb1 = 4.0
    ! bb2 = 0.8

    ! do i = 1, im
    !     if (cnvflg (i)) then
    !         k = kbcon1 (i)
    !         tem = po (i, k) / (rdgas * to (i, k))
    !         wucb = - 0.01 * dot (i, k) / (tem * g)
    !         if (wucb > 0.) then
    !             wu2 (i, k) = wucb * wucb
    !         else
    !             wu2 (i, k) = 0.
    !         endif
    !     endif
    ! enddo

    do k = 2, km1
        do i = 1, im
            if (cnvflg (i)) then
                if (k > kbcon1 (i) .and. k < ktcon (i)) then
                    dz = zi (i, k) - zi (i, k - 1)
                    tem = 0.25 * bb1 * (drag (i, k) + drag (i, k - 1)) * dz
                    tem1 = 0.5 * bb2 * (buo (i, k) + buo (i, k - 1)) * dz
                    ! kgao 12 / 18 / 2023
                    if (use_shear_conv) then
                        tem2 = wush (i, k) * sqrt (wu2 (i, k - 1))
                        tem2 = (tem1 - tem2) * dz
                        ptem = (1. - tem) * wu2 (i, k - 1)
                        ptem1 = 1. + tem
                        wu2 (i, k) = (ptem + tem2) / ptem1
                    else
                        ptem = (1. - tem) * wu2 (i, k - 1)
                        ptem1 = 1. + tem
                        wu2 (i, k) = (ptem + tem1) / ptem1
                    endif
                    wu2 (i, k) = max (wu2 (i, k), 0.)
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute updraft velocity average over the whole cumulus
    ! calculate the mean updraft velocity within the cloud (wc) .
    ! -----------------------------------------------------------------------

    do i = 1, im
        wc (i) = 0.
        sumx (i) = 0.
    enddo

    do k = 2, km1
        do i = 1, im
            if (cnvflg (i)) then
                if (k > kbcon1 (i) .and. k < ktcon (i)) then
                    dz = zi (i, k) - zi (i, k - 1)
                    tem = 0.5 * (sqrt (wu2 (i, k)) + sqrt (wu2 (i, k - 1)))
                    wc (i) = wc (i) + tem * dz
                    sumx (i) = sumx (i) + dz
                endif
            endif
        enddo
    enddo

    do i = 1, im
        if (cnvflg (i)) then
            if (sumx (i) == 0.) then
                cnvflg (i) = .false.
            else
                wc (i) = wc (i) / sumx (i)
            endif
            val = 1.e-4
            if (wc (i) < val) cnvflg (i) = .false.
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! exchange ktcon with ktcon1
    ! swap the indices of the convective cloud top (ktcon) and the overshooting convection top (ktcon1) to use the same cloud top level in the calculations of \f$a^ + \f$ and \f$a^ * \f$.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            kk = ktcon (i)
            ktcon (i) = ktcon1 (i)
            ktcon1 (i) = kk
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! this section is ready for cloud water
    ! separate the total updraft cloud water at cloud top into vapor and condensate.
    ! -----------------------------------------------------------------------

    if (ncloud > 0) then

        ! -----------------------------------------------------------------------
        ! compute liquid and vapor separation at cloud top
        ! -----------------------------------------------------------------------

        do i = 1, im
            if (cnvflg (i)) then
                k = ktcon (i) - 1
                gamma = el2orc * qeso (i, k) / (to (i, k) ** 2)
                qrch = qeso (i, k) &
                     + gamma * dbyo (i, k) / (hlv * (1. + gamma))
                dq = qcko (i, k) - qrch

                ! -----------------------------------------------------------------------
                ! check if there is excess moisture to release latent heat
                ! -----------------------------------------------------------------------

                if (dq > 0.) then
                    qlko_ktcon (i) = dq
                    qcko (i, k) = qrch
                endif
            endif
        enddo
    endif

    ! -----------------------------------------------------------------------
    ! downdraft calculations
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! compute precipitation efficiency in terms of windshear
    ! perform calculations related to the downdraft of the entraining / detraining cloud model ("static control") .
    ! first, in order to calculate the downdraft mass flux (as a fraction of the updraft mass flux), calculate the wind shear and precipitation efficiency according to equation 58 in fritsch and chappell (1980):
    ! \f[
    ! e = 1.591 - 0.639\frac{\delta v}{\delta z} + 0.0953\left (\frac{\delta v}{\delta z}\right) ^2 - 0.00496\left (\frac{\delta v}{\delta z}\right) ^3
    ! \f]
    ! where \f$\delta v\f$ is the integrated horizontal shear over the cloud depth, \f$\delta z\f$, (the ratio is converted to units of \f$10^{ - 3} s^{ - 1}\f$) . the variable "edto" is \f$1 - e\f$ and is constrained to the range \f$[0, 0.9]\f$.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            vshear (i) = 0.
        endif
    enddo

    do k = 2, km
        do i = 1, im
            if (cnvflg (i)) then
                if (k > kb (i) .and. k <= ktcon (i)) then
                    shear = sqrt ((uo (i, k) - uo (i, k - 1)) ** 2 &
                         + (vo (i, k) - vo (i, k - 1)) ** 2)
                    vshear (i) = vshear (i) + shear
                endif
            endif
        enddo
    enddo

    do i = 1, im
        if (cnvflg (i)) then
            vshear (i) = 1.e3 * vshear (i) / (zi (i, ktcon (i)) - zi (i, kb (i)))
            e1 = 1.591 - .639 * vshear (i) &
                 + .0953 * (vshear (i) ** 2) - .00496 * (vshear (i) ** 3)
            edt (i) = 1. - e1
            val = .9
            edt (i) = min (edt (i), val)
            val = .0
            edt (i) = max (edt (i), val)
            edto (i) = edt (i)
            edtx (i) = edt (i)
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! determine detrainment rate between 1 and kbcon
    ! next, calculate the variable detrainment rate between the surface and the lfc according to:
    ! \f[
    ! \lambda_d = \frac{1 - \beta^{\frac{1}{k_{lfc}}}}{\overline{\delta z}}
    ! \f]
    ! \f$\lambda_d\f$ is the detrainment rate, \f$\beta\f$ is a constant currently set to 0.05, implying that only 5% of downdraft mass flux at lfc reaches the ground surface due to detrainment, \f$k_{lfc}\f$ is the vertical index of the lfc level, and \f$\overline{\delta z}\f$ is the average vertical grid spacing below the lfc.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            sumx (i) = 0.
        endif
    enddo

    do k = 1, km1
        do i = 1, im
            if (cnvflg (i)) then
                if (k >= 1 .and. k < kbcon (i)) then
                    dz = zi (i, k + 1) - zi (i, k)
                    sumx (i) = sumx (i) + dz
                endif
            endif
        enddo
    enddo

    do i = 1, im
        if (cnvflg (i)) then
            betamn = betas_deep
            if (islimsk (i) == 1) betamn = betal_deep
            if (ntk > 0) then
                betamx = betamn + dbeta
                if (tkemean (i) > tkemx) then
                    beta = betamn
                else if (tkemean (i) < tkemn) then
                    beta = betamx
                else
                    tem = (betamx - betamn) * (tkemean (i) - tkemn)
                    beta = betamx - tem / dtke
                endif
            else
                beta = betamn
            endif
            dz = (sumx (i) + zi (i, 1)) / float (kbcon (i))
            tem = 1. / float (kbcon (i))
            xlamd (i) = (1. - beta ** tem) / dz
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! determine downdraft mass flux
    ! calculate the normalized downdraft mass flux from equation 1 of pan and wu (1995). downdraft entrainment and detrainment rates are constants from the downdraft origination to the lfc.
    ! -----------------------------------------------------------------------

    do k = km1, 1, - 1
        do i = 1, im
            if (cnvflg (i) .and. k <= kmax (i) - 1) then
                if (k < jmin (i) .and. k >= kbcon (i)) then
                    dz = zi (i, k + 1) - zi (i, k)
                    ! kgao 12 / 18 / 2023
                    if (use_tke_conv) then
                        ptem = xlamddt (i) - xlamdet (i)
                    else
                        ptem = xlamdd - xlamde
                    endif
                    etad (i, k) = etad (i, k + 1) * (1. - ptem * dz)
                else if (k < kbcon (i)) then
                    dz = zi (i, k + 1) - zi (i, k)
                    ! kgao 12 / 18 / 2023
                    if (use_tke_conv) then
                        ptem = xlamd (i) + xlamddt (i) - xlamdet (i)
                    else
                        ptem = xlamd (i) + xlamdd - xlamde
                    endif
                    etad (i, k) = etad (i, k + 1) * (1. - ptem * dz)
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! downdraft moisture properties
    ! set initial cloud downdraft properties equal to the state variables at the downdraft origination level.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            jmn = jmin (i)
            hcdo (i, jmn) = heo (i, jmn)
            qcdo (i, jmn) = qo (i, jmn)
            qrcdo (i, jmn) = qo (i, jmn)
            ucdo (i, jmn) = uo (i, jmn)
            vcdo (i, jmn) = vo (i, jmn)
            pwevo (i) = 0.
        endif
    enddo

    ! for tracers
    do n = 1, ntr
        do i = 1, im
            if (cnvflg (i)) then
                jmn = jmin (i)
                ecdo (i, jmn, n) = ctro (i, jmn, n)
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! calculate the cloud properties as a parcel descends, modified by entrainment and detrainment. discretization follows appendix b of grell (1993).
    ! -----------------------------------------------------------------------

    do k = km1, 1, - 1
        do i = 1, im
            if (cnvflg (i) .and. k < jmin (i)) then
                dz = zi (i, k + 1) - zi (i, k)
                if (k >= kbcon (i)) then
                    ! kgao 12 / 18 / 2023
                    if (use_tke_conv) then
                        tem = xlamdet (i) * dz
                        tem1 = 0.5 * xlamddt (i) * dz
                    else
                        tem = xlamde * dz
                        tem1 = 0.5 * xlamdd * dz
                    endif
                else
                    ! kgao 12 / 18 / 2023
                    if (use_tke_conv) then
                        tem = xlamdet (i) * dz
                        tem1 = 0.5 * (xlamd (i) + xlamddt (i)) * dz
                    else
                        tem = xlamde * dz
                        tem1 = 0.5 * (xlamd (i) + xlamdd) * dz
                    endif
                endif
                factor = 1. + tem - tem1
                hcdo (i, k) = ((1. - tem1) * hcdo (i, k + 1) + tem * 0.5 * &
                     (heo (i, k) + heo (i, k + 1))) / factor
                dbyo (i, k) = hcdo (i, k) - heso (i, k)

                tem = 0.5 * cm * tem
                factor = 1. + tem
                ptem = tem - pgcon_deep
                ptem1 = tem + pgcon_deep
                ucdo (i, k) = ((1. - tem) * ucdo (i, k + 1) + ptem * uo (i, k + 1) &
                     + ptem1 * uo (i, k)) / factor
                vcdo (i, k) = ((1. - tem) * vcdo (i, k + 1) + ptem * vo (i, k + 1) &
                     + ptem1 * vo (i, k)) / factor
            endif
        enddo
    enddo

    do n = 1, ntr
        do k = km1, 1, - 1
            do i = 1, im
                if (cnvflg (i) .and. k < jmin (i)) then
                    dz = zi (i, k + 1) - zi (i, k)
                    ! kgao 12 / 18 / 2023
                    if (use_tke_conv) then
                        tem = 0.5 * xlamdet (i) * dz
                    else
                        tem = 0.5 * xlamde * dz
                    endif
                    factor = 1. + tem
                    ecdo (i, k, n) = ((1. - tem) * ecdo (i, k + 1, n) + tem * &
                        (ctro (i, k, n) + ctro (i, k + 1, n))) / factor
                endif
            enddo
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute the amount of moisture that is necessary to keep the downdraft saturated.
    ! -----------------------------------------------------------------------

    do k = km1, 1, - 1
        do i = 1, im
            if (cnvflg (i) .and. k < jmin (i)) then
                gamma = el2orc * qeso (i, k) / (to (i, k) ** 2)
                qrcdo (i, k) = qeso (i, k) + &
                     (1. / hlv) * (gamma / (1. + gamma)) * dbyo (i, k)
                ! detad = etad (i, k + 1) - etad (i, k)

                dz = zi (i, k + 1) - zi (i, k)
                ! kgao 12 / 18 / 2023
                if (use_tke_conv) then
                    if (k >= kbcon (i)) then
                        tem = xlamdet (i) * dz
                        tem1 = 0.5 * xlamddt (i) * dz
                    else
                        tem = xlamdet (i) * dz
                        tem1 = 0.5 * (xlamd (i) + xlamddt (i)) * dz
                    endif
                else
                    if (k >= kbcon (i)) then
                        tem = xlamde * dz
                        tem1 = 0.5 * xlamdd * dz
                    else
                        tem = xlamde * dz
                        tem1 = 0.5 * (xlamd (i) + xlamdd) * dz
                    endif
                endif
                factor = 1. + tem - tem1
                qcdo (i, k) = ((1. - tem1) * qrcdo (i, k + 1) + tem * 0.5 * &
                     (qo (i, k) + qo (i, k + 1))) / factor

                ! pwdo (i, k) = etad (i, k + 1) * qcdo (i, k + 1) - &
                ! etad (i, k) * qrcdo (i, k)
                ! pwdo (i, k) = pwdo (i, k) - detad * &
                ! .5 * (qrcdo (i, k) + qrcdo (i, k + 1))

                pwdo (i, k) = etad (i, k) * (qcdo (i, k) - qrcdo (i, k))
                pwevo (i) = pwevo (i) + pwdo (i, k)
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! final downdraft strength dependent on precip
    ! efficiency (edt), normalized condensate (pwav), and
    ! evaporate (pwev)
    ! update the precipitation efficiency (edto) based on the ratio of normalized cloud condensate (pwavo) to normalized cloud evaporate (pwevo) .
    ! -----------------------------------------------------------------------

    do i = 1, im
        edtmax = edtmaxl
        if (islimsk (i) == 0) edtmax = edtmaxs
        if (cnvflg (i)) then
            if (pwevo (i) < 0.) then
                edto (i) = - edto (i) * pwavo (i) / pwevo (i)
                edto (i) = min (edto (i), edtmax)
            else
                edto (i) = 0.
            endif
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! downdraft cloudwork functions
    ! calculate downdraft cloud work function (\f$a_d\f$) according to equation a.42 (discretized by b.11) in grell (1993). add it to the updraft cloud work function, \f$a_u\f$.
    ! -----------------------------------------------------------------------

    do k = km1, 1, - 1
        do i = 1, im
            if (cnvflg (i) .and. k < jmin (i)) then
                gamma = el2orc * qeso (i, k) / to (i, k) ** 2
                dhh = hcdo (i, k)
                dt = to (i, k)
                dg = gamma
                dh = heso (i, k)
                dz = - 1. * (zo (i, k + 1) - zo (i, k))
                ! aa1 (i) = aa1 (i) + edto (i) * dz * etad (i, k)
                aa1 (i) = aa1 (i) + edto (i) * dz &
                     * (g / (cp_air * dt)) * ((dhh - dh) / (1. + dg)) &
                     * (1. + delta * cp_air * dg * dt / hlv)
                val = 0.
                ! aa1 (i) = aa1 (i) + edto (i) * dz * etad (i, k)
                aa1 (i) = aa1 (i) + edto (i) * dz &
                     * g * delta * max (val, (qeso (i, k) - qo (i, k)))
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! check for negative total cloud work function; if found, return to calling routine without modifying state variables.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i) .and. aa1 (i) <= 0.) then
            cnvflg (i) = .false.
        endif
    enddo

    totflg = .true.
    do i = 1, im
        totflg = totflg .and. (.not. cnvflg (i))
    enddo

    if (totflg) return

    ! -----------------------------------------------------------------------
    ! what would the change be, that a cloud with unit mass
    ! will do to the environment?
    ! calculate the change in moist static energy, moisture mixing ratio, and horizontal winds per unit cloud base mass flux near the surface using equations b.18 and b.19 from grell (1993), for all layers below cloud top from equations b.14 and b.15, and for the cloud top from b.16 and b.17.
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (cnvflg (i) .and. k <= kmax (i)) then
                dellah (i, k) = 0.
                dellaq (i, k) = 0.
                dellau (i, k) = 0.
                dellav (i, k) = 0.
            endif
        enddo
    enddo

    do n = 1, ntr
        do k = 1, km
            do i = 1, im
                if (cnvflg (i) .and. k <= kmax (i)) then
                    dellae (i, k, n) = 0.
                endif
            enddo
        enddo
    enddo

    do i = 1, im
        if (cnvflg (i)) then
            dp = delp (i, 1)
            dellah (i, 1) = edto (i) * etad (i, 1) * (hcdo (i, 1) &
                 - heo (i, 1)) * g / dp
            dellaq (i, 1) = edto (i) * etad (i, 1) * (qrcdo (i, 1) &
                 - qo (i, 1)) * g / dp
            dellau (i, 1) = edto (i) * etad (i, 1) * (ucdo (i, 1) &
                 - uo (i, 1)) * g / dp
            dellav (i, 1) = edto (i) * etad (i, 1) * (vcdo (i, 1) &
                 - vo (i, 1)) * g / dp
        endif
    enddo

    do n = 1, ntr
        do i = 1, im
            if (cnvflg (i)) then
                dp = delp (i, 1)
                dellae (i, 1, n) = edto (i) * etad (i, 1) * (ecdo (i, 1, n) &
                    - ctro (i, 1, n)) * g / dp
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! changed due to subsidence and entrainment
    ! -----------------------------------------------------------------------

    do k = 2, km1
        do i = 1, im
            if (cnvflg (i) .and. k < ktcon (i)) then
                aup = 1.
                if (k <= kb (i)) aup = 0.
                adw = 1.
                if (k > jmin (i)) adw = 0.
                dp = delp (i, k)
                dz = zi (i, k) - zi (i, k - 1)

                dv1h = heo (i, k)
                dv2h = .5 * (heo (i, k) + heo (i, k - 1))
                dv3h = heo (i, k - 1)
                dv1q = qo (i, k)
                dv2q = .5 * (qo (i, k) + qo (i, k - 1))
                dv3q = qo (i, k - 1)

                tem = 0.5 * (xlamue (i, k) + xlamue (i, k - 1))
                tem1 = 0.5 * (xlamud (i, k) + xlamud (i, k - 1))

                ! kgao 12 / 18 / 2023
                if (use_tke_conv) then
                    if (k <= kbcon (i)) then
                        ptem = xlamdet (i)
                        ptem1 = xlamd (i) + xlamddt (i)
                    else
                        ptem = xlamdet (i)
                        ptem1 = xlamddt (i)
                    endif
                else
                    if (k <= kbcon (i)) then
                        ptem = xlamde
                        ptem1 = xlamd (i) + xlamdd
                    else
                        ptem = xlamde
                        ptem1 = xlamdd
                    endif
                endif

                dellah (i, k) = dellah (i, k) + &
                     ((aup * eta (i, k) - adw * edto (i) * etad (i, k)) * dv1h &
                     - (aup * eta (i, k - 1) - adw * edto (i) * etad (i, k - 1)) * dv3h &
                     - (aup * tem * eta (i, k - 1) + adw * edto (i) * ptem * etad (i, k)) * dv2h * dz &
                     + aup * tem1 * eta (i, k - 1) * .5 * (hcko (i, k) + hcko (i, k - 1)) * dz &
                     + adw * edto (i) * ptem1 * etad (i, k) * .5 * (hcdo (i, k) + hcdo (i, k - 1)) * dz &
                    ) * g / dp

                dellaq (i, k) = dellaq (i, k) + &
                     ((aup * eta (i, k) - adw * edto (i) * etad (i, k)) * dv1q &
                     - (aup * eta (i, k - 1) - adw * edto (i) * etad (i, k - 1)) * dv3q &
                     - (aup * tem * eta (i, k - 1) + adw * edto (i) * ptem * etad (i, k)) * dv2q * dz &
                     + aup * tem1 * eta (i, k - 1) * .5 * (qrcko (i, k) + qcko (i, k - 1)) * dz &
                     + adw * edto (i) * ptem1 * etad (i, k) * .5 * (qrcdo (i, k) + qcdo (i, k - 1)) * dz &
                    ) * g / dp

                tem1 = eta (i, k) * (uo (i, k) - ucko (i, k))
                tem2 = eta (i, k - 1) * (uo (i, k - 1) - ucko (i, k - 1))
                ptem1 = etad (i, k) * (uo (i, k) - ucdo (i, k))
                ptem2 = etad (i, k - 1) * (uo (i, k - 1) - ucdo (i, k - 1))
                dellau (i, k) = dellau (i, k) + &
                     (aup * (tem1 - tem2) - adw * edto (i) * (ptem1 - ptem2)) * g / dp

                tem1 = eta (i, k) * (vo (i, k) - vcko (i, k))
                tem2 = eta (i, k - 1) * (vo (i, k - 1) - vcko (i, k - 1))
                ptem1 = etad (i, k) * (vo (i, k) - vcdo (i, k))
                ptem2 = etad (i, k - 1) * (vo (i, k - 1) - vcdo (i, k - 1))
                dellav (i, k) = dellav (i, k) + &
                     (aup * (tem1 - tem2) - adw * edto (i) * (ptem1 - ptem2)) * g / dp

            endif
        enddo
    enddo

    do n = 1, ntr
        do k = 2, km1
            do i = 1, im
                if (cnvflg (i) .and. k < ktcon (i)) then
                    aup = 1.
                    if (k <= kb (i)) aup = 0.
                    adw = 1.
                    if (k > jmin (i)) adw = 0.
                    dp = delp (i, k)
                    tem1 = eta (i, k) * (ctro (i, k, n) - ecko (i, k, n))
                    tem2 = eta (i, k - 1) * (ctro (i, k - 1, n) - ecko (i, k - 1, n))
                    ptem1 = etad (i, k) * (ctro (i, k, n) - ecdo (i, k, n))
                    ptem2 = etad (i, k - 1) * (ctro (i, k - 1, n) - ecdo (i, k - 1, n))
                    dellae (i, k, n) = dellae (i, k, n) + &
                        (aup * (tem1 - tem2) - adw * edto (i) * (ptem1 - ptem2)) * g / dp
                endif
            enddo
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! cloud top
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            indx = ktcon (i)
            dp = delp (i, indx)
            dv1h = heo (i, indx - 1)
            dellah (i, indx) = eta (i, indx - 1) * &
                 (hcko (i, indx - 1) - dv1h) * g / dp
            dv1q = qo (i, indx - 1)
            dellaq (i, indx) = eta (i, indx - 1) * &
                 (qcko (i, indx - 1) - dv1q) * g / dp
            dellau (i, indx) = eta (i, indx - 1) * &
                 (ucko (i, indx - 1) - uo (i, indx - 1)) * g / dp
            dellav (i, indx) = eta (i, indx - 1) * &
                 (vcko (i, indx - 1) - vo (i, indx - 1)) * g / dp

            ! -----------------------------------------------------------------------
            ! cloud water
            ! -----------------------------------------------------------------------

            dellal (i, indx) = eta (i, indx - 1) * &
                qlko_ktcon (i) * g / dp
        endif
    enddo

    do n = 1, ntr
        do i = 1, im
            if (cnvflg (i)) then
                indx = ktcon (i)
                dp = delp (i, indx)
                dellae (i, indx, n) = eta (i, indx - 1) * &
                    (ecko (i, indx - 1, n) - ctro (i, indx - 1, n)) * g / dp
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! final changed variable per unit mass flux
    ! if grid size is less than a threshold value (dxcrtas_deep: currently 8km), the quasi - equilibrium assumption of arakawa - schubert is not used any longer.
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! if grid size is less than a threshold value (dxcrtas_deep),
    ! the quasi - equilibrium assumption of arakawa - schubert is not
    ! used any longer.
    ! -----------------------------------------------------------------------

    do i = 1, im
        asqecflg (i) = cnvflg (i)
        if (asqecflg (i) .and. gsize (i) < dxcrtas_deep) then
            asqecflg (i) = .false.
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! if grid size is larger than the threshold value (i.e., asqecflg = .true.), the quasi - equilibrium assumption is used to obtain the cloud base mass flux. to begin with, calculate the change in the temperature and moisture profiles per unit cloud base mass flux.
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (asqecflg (i) .and. k <= kmax (i)) then
                if (k > ktcon (i)) then
                    qo (i, k) = q1 (i, k)
                    to (i, k) = t1 (i, k)
                endif
                if (k <= ktcon (i)) then
                    qo (i, k) = dellaq (i, k) * mbdt (i) + q1 (i, k)
                    dellat = (dellah (i, k) - hlv * dellaq (i, k)) / cp_air
                    to (i, k) = dellat * mbdt (i) + t1 (i, k)
                    val = 1.e-10
                    qo (i, k) = max (qo (i, k), val)
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! the above changed environment is now used to calulate the
    ! effect the arbitrary cloud (with unit mass flux)
    ! would have on the stability,
    ! which then is used to calculate the real mass flux,
    ! necessary to keep this change in balance with the large - scale
    ! destabilization.
    ! environmental conditions again, first heights
    ! using the updated temperature and moisture profiles that were modified by the convection on a short time - scale, recalculate the total cloud work function to determine the change in the cloud work function due to convection, or the stabilizing effect of the cumulus.
    ! using notation from pan and wu (1995), the previously calculated cloud work function is denoted by \f$a^ + \f$. now, it is necessary to use the entraining / detraining cloud model ("static control") to determine the cloud work function of the environment after the stabilization of the arbitrary convective element (per unit cloud base mass flux) has been applied, denoted by \f$a^ * \f$.
    ! recalculate saturation specific humidity.
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (asqecflg (i) .and. k <= kmax (i)) then
                qeso (i, k) = 0.01 * mqs (to (i, k)) ! mqs is in pa
                qeso (i, k) = eps * qeso (i, k) / (pfld (i, k) + epsm1 * qeso (i, k))
                val = 1.e-8
                qeso (i, k) = max (qeso (i, k), val)
                ! tvo (i, k) = to (i, k) + delta * to (i, k) * qo (i, k)
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! moist static energy
    ! recalculate moist static energy and saturation moist static energy.
    ! -----------------------------------------------------------------------

    do k = 1, km1
        do i = 1, im
            if (asqecflg (i) .and. k <= kmax (i) - 1) then
                dz = .5 * (zo (i, k + 1) - zo (i, k))
                dp = .5 * (pfld (i, k + 1) - pfld (i, k))
                es = 0.01 * mqs (to (i, k + 1)) ! mqs is in pa
                pprime = pfld (i, k + 1) + epsm1 * es
                qs = eps * es / pprime
                dqsdp = - qs / pprime
                desdt = es * (fact1 / to (i, k + 1) + fact2 / (to (i, k + 1) ** 2))
                dqsdt = qs * pfld (i, k + 1) * desdt / (es * pprime)
                gamma = el2orc * qeso (i, k + 1) / (to (i, k + 1) ** 2)
                dt = (g * dz + hlv * dqsdp * dp) / (cp_air * (1. + gamma))
                dq = dqsdt * dt + dqsdp * dp
                to (i, k) = to (i, k + 1) + dt
                qo (i, k) = qo (i, k + 1) + dq
                po (i, k) = .5 * (pfld (i, k) + pfld (i, k + 1))
            endif
        enddo
    enddo

    do k = 1, km1
        do i = 1, im
            if (asqecflg (i) .and. k <= kmax (i) - 1) then
                qeso (i, k) = 0.01 * mqs (to (i, k)) ! mqs is in pa
                qeso (i, k) = eps * qeso (i, k) / (po (i, k) + epsm1 * qeso (i, k))
                val1 = 1.e-8
                qeso (i, k) = max (qeso (i, k), val1)
                val2 = 1.e-10
                qo (i, k) = max (qo (i, k), val2)
                ! qo (i, k) = min (qo (i, k), qeso (i, k))
                heo (i, k) = .5 * g * (zo (i, k) + zo (i, k + 1)) + &
                    cp_air * to (i, k) + hlv * qo (i, k)
                heso (i, k) = .5 * g * (zo (i, k) + zo (i, k + 1)) + &
                    cp_air * to (i, k) + hlv * qeso (i, k)
            endif
        enddo
    enddo

    do i = 1, im
        if (asqecflg (i)) then
            k = kmax (i)
            heo (i, k) = g * zo (i, k) + cp_air * to (i, k) + hlv * qo (i, k)
            heso (i, k) = g * zo (i, k) + cp_air * to (i, k) + hlv * qeso (i, k)
            ! heo (i, k) = min (heo (i, k), heso (i, k))
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! static control
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! moisture and cloud work functions
    ! as before, recalculate the updraft cloud work function.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (asqecflg (i)) then
            xaa0 (i) = 0.
            xpwav (i) = 0.
        endif
    enddo

    do i = 1, im
        if (asqecflg (i)) then
            indx = kb (i)
            hcko (i, indx) = heo (i, indx)
            qcko (i, indx) = qo (i, indx)
        endif
    enddo

    do k = 2, km1
        do i = 1, im
            if (asqecflg (i)) then
                if (k > kb (i) .and. k <= ktcon (i)) then
                    dz = zi (i, k) - zi (i, k - 1)
                    tem = 0.5 * (xlamue (i, k) + xlamue (i, k - 1)) * dz
                    tem1 = 0.25 * (xlamud (i, k) + xlamud (i, k - 1)) * dz
                    factor = 1. + tem - tem1
                    hcko (i, k) = ((1. - tem1) * hcko (i, k - 1) + tem * 0.5 * &
                         (heo (i, k) + heo (i, k - 1))) / factor
                endif
            endif
        enddo
    enddo

    do k = 2, km1
        do i = 1, im
            if (asqecflg (i)) then
                if (k > kb (i) .and. k < ktcon (i)) then
                    dz = zi (i, k) - zi (i, k - 1)
                    gamma = el2orc * qeso (i, k) / (to (i, k) ** 2)
                    xdby = hcko (i, k) - heso (i, k)
                    xqrch = qeso (i, k) &
                         + gamma * xdby / (hlv * (1. + gamma))

                    tem = 0.5 * (xlamue (i, k) + xlamue (i, k - 1)) * dz
                    tem1 = 0.25 * (xlamud (i, k) + xlamud (i, k - 1)) * dz
                    factor = 1. + tem - tem1
                    qcko (i, k) = ((1. - tem1) * qcko (i, k - 1) + tem * 0.5 * &
                         (qo (i, k) + qo (i, k - 1))) / factor

                    dq = eta (i, k) * (qcko (i, k) - xqrch)

                    if (k >= kbcon (i) .and. dq > 0.) then
                        etah = .5 * (eta (i, k) + eta (i, k - 1))
                        if (ncloud > 0 .and. k > jmin (i)) then
                            ptem = c0t (i, k) + c1_deep
                            qlk = dq / (eta (i, k) + etah * ptem * dz)
                        else
                            qlk = dq / (eta (i, k) + etah * c0t (i, k) * dz)
                        endif
                        if (k < ktcon1 (i)) then
                            ! xaa0 (i) = xaa0 (i) - dz * g * qlk * etah
                            xaa0 (i) = xaa0 (i) - dz * g * qlk
                        endif
                        qcko (i, k) = qlk + xqrch
                        xpw = etah * c0t (i, k) * dz * qlk
                        xpwav (i) = xpwav (i) + xpw
                    endif
                endif
                if (k >= kbcon (i) .and. k < ktcon1 (i)) then
                    dz1 = zo (i, k + 1) - zo (i, k)
                    gamma = el2orc * qeso (i, k) / (to (i, k) ** 2)
                    rfact = 1. + delta * cp_air * gamma &
                         * to (i, k) / hlv
                    xaa0 (i) = xaa0 (i) &
                        ! + dz1 * eta (i, k) * (g / (cp_air * to (i, k))) &
                         + dz1 * (g / (cp_air * to (i, k))) &
                         * xdby / (1. + gamma) &
                         * rfact
                    val = 0.
                    xaa0 (i) = xaa0 (i) + &
                        ! dz1 * eta (i, k) * g * delta * &
                        dz1 * g * delta * &
                        max (val, (qeso (i, k) - qo (i, k)))
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! downdraft calculations
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! downdraft moisture properties
    ! as before, recalculate the downdraft cloud work function.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (asqecflg (i)) then
            jmn = jmin (i)
            hcdo (i, jmn) = heo (i, jmn)
            qcdo (i, jmn) = qo (i, jmn)
            qrcd (i, jmn) = qo (i, jmn)
            xpwev (i) = 0.
        endif
    enddo

    do k = km1, 1, - 1
        do i = 1, im
            if (asqecflg (i) .and. k < jmin (i)) then
                dz = zi (i, k + 1) - zi (i, k)
                ! kgao 12 / 18 / 2023
                if (use_tke_conv) then
                    if (k >= kbcon (i)) then
                        tem = xlamdet (i) * dz
                        tem1 = 0.5 * xlamddt (i) * dz
                    else
                        tem = xlamdet (i) * dz
                        tem1 = 0.5 * (xlamd (i) + xlamddt (i)) * dz
                    endif
                else
                    if (k >= kbcon (i)) then
                        tem = xlamde * dz
                        tem1 = 0.5 * xlamdd * dz
                    else
                        tem = xlamde * dz
                        tem1 = 0.5 * (xlamd (i) + xlamdd) * dz
                    endif
                endif
                factor = 1. + tem - tem1
                hcdo (i, k) = ((1. - tem1) * hcdo (i, k + 1) + tem * 0.5 * &
                     (heo (i, k) + heo (i, k + 1))) / factor
            endif
        enddo
    enddo

    do k = km1, 1, - 1
        do i = 1, im
            if (asqecflg (i) .and. k < jmin (i)) then
                dq = qeso (i, k)
                dt = to (i, k)
                gamma = el2orc * dq / dt ** 2
                dh = hcdo (i, k) - heso (i, k)
                qrcd (i, k) = dq + (1. / hlv) * (gamma / (1. + gamma)) * dh
                ! detad = etad (i, k + 1) - etad (i, k)

                dz = zi (i, k + 1) - zi (i, k)
                ! kgao 12 / 18 / 2023
                if (use_tke_conv) then
                    if (k >= kbcon (i)) then
                        tem = xlamdet (i) * dz
                        tem1 = 0.5 * xlamddt (i) * dz
                    else
                        tem = xlamdet (i) * dz
                        tem1 = 0.5 * (xlamd (i) + xlamddt (i)) * dz
                    endif
                else
                    if (k >= kbcon (i)) then
                        tem = xlamde * dz
                        tem1 = 0.5 * xlamdd * dz
                    else
                        tem = xlamde * dz
                        tem1 = 0.5 * (xlamd (i) + xlamdd) * dz
                    endif
                endif
                factor = 1. + tem - tem1
                qcdo (i, k) = ((1. - tem1) * qrcd (i, k + 1) + tem * 0.5 * &
                     (qo (i, k) + qo (i, k + 1))) / factor

                ! xpwd = etad (i, k + 1) * qcdo (i, k + 1) - &
                ! etad (i, k) * qrcd (i, k)
                ! xpwd = xpwd - detad * &
                ! .5 * (qrcd (i, k) + qrcd (i, k + 1))

                xpwd = etad (i, k) * (qcdo (i, k) - qrcd (i, k))
                xpwev (i) = xpwev (i) + xpwd
            endif
        enddo
    enddo

    do i = 1, im
        edtmax = edtmaxl
        if (islimsk (i) == 0) edtmax = edtmaxs
        if (asqecflg (i)) then
            if (xpwev (i) >= 0.) then
                edtx (i) = 0.
            else
                edtx (i) = - edtx (i) * xpwav (i) / xpwev (i)
                edtx (i) = min (edtx (i), edtmax)
            endif
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! downdraft cloudwork functions
    ! -----------------------------------------------------------------------

    do k = km1, 1, - 1
        do i = 1, im
            if (asqecflg (i) .and. k < jmin (i)) then
                gamma = el2orc * qeso (i, k) / to (i, k) ** 2
                dhh = hcdo (i, k)
                dt = to (i, k)
                dg = gamma
                dh = heso (i, k)
                dz = - 1. * (zo (i, k + 1) - zo (i, k))
                ! xaa0 (i) = xaa0 (i) + edtx (i) * dz * etad (i, k)
                xaa0 (i) = xaa0 (i) + edtx (i) * dz &
                     * (g / (cp_air * dt)) * ((dhh - dh) / (1. + dg)) &
                     * (1. + delta * cp_air * dg * dt / hlv)
                val = 0.
                ! xaa0 (i) = xaa0 (i) + edtx (i) * dz * etad (i, k)
                xaa0 (i) = xaa0 (i) + edtx (i) * dz &
                     * g * delta * max (val, (qeso (i, k) - qo (i, k)))
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! calculate critical cloud work function
    ! -----------------------------------------------------------------------

    ! do i = 1, im
    ! if (cnvflg (i)) then
    ! if (pfld (i, ktcon (i)) < pcrit (15)) then
    ! acrt (i) = acrit (15) * (975. - pfld (i, ktcon (i))) &
    ! / (975. - pcrit (15))
    ! else if (pfld (i, ktcon (i)) > pcrit (1)) then
    ! acrt (i) = acrit (1)
    ! else
    ! k = int ((850. - pfld (i, ktcon (i))) / 50.) + 2
    ! k = min (k, 15)
    ! k = max (k, 2)
    ! acrt (i) = acrit (k) + (acrit (k - 1) - acrit (k)) * &
    ! (pfld (i, ktcon (i)) - pcrit (k)) / (pcrit (k - 1) - pcrit (k))
    ! endif
    ! endif
    ! enddo
    ! do i = 1, im
    ! if (cnvflg (i)) then
    ! if (islimsk (i) == 1) then
    ! w1 = w1l
    ! w2 = w2l
    ! w3 = w3l
    ! w4 = w4l
    ! else
    ! w1 = w1s
    ! w2 = w2s
    ! w3 = w3s
    ! w4 = w4s
    ! endif

    ! -----------------------------------------------------------------------
    ! modify critical cloud workfunction by cloud base vertical velocity
    ! -----------------------------------------------------------------------

    ! if (pdot (i) <= w4) then
    ! acrtfct (i) = (pdot (i) - w4) / (w3 - w4)
    ! elseif (pdot (i) >= - w4) then
    ! acrtfct (i) = - (pdot (i) + w4) / (w4 - w3)
    ! else
    ! acrtfct (i) = 0.
    ! endif
    ! val1 = - 1.
    ! acrtfct (i) = max (acrtfct (i), val1)
    ! val2 = 1.
    ! acrtfct (i) = min (acrtfct (i), val2)
    ! acrtfct (i) = 1. - acrtfct (i)

    ! -----------------------------------------------------------------------
    ! modify acrtfct (i) by colume mean rh if rhbar (i) is greater than 80 percent
    ! -----------------------------------------------------------------------

    ! if (rhbar (i) >= .8) then
    ! acrtfct (i) = acrtfct (i) * (.9 - min (rhbar (i), .9)) * 10.
    ! endif

    ! -----------------------------------------------------------------------
    ! modify adjustment time scale by cloud base vertical velocity
    ! -----------------------------------------------------------------------

    ! dtconv (i) = dt2 + max ((1800. - dt2), 0.) * &
    ! (pdot (i) - w2) / (w1 - w2)
    ! dtconv (i) = max (dtconv (i), dt2)
    ! dtconv (i) = 1800. * (pdot (i) - w2) / (w1 - w2)

    ! dtconv (i) = max (dtconv (i), dtmin)
    ! dtconv (i) = min (dtconv (i), dtmax)

    ! endif
    ! enddo

    ! -----------------------------------------------------------------------
    ! compute convective turn - over time
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! following bechtold et al. (2008), the convective adjustment time (dtconv) is set to be proportional to the convective turnover time, which is computed using the mean updraft velocity (wc) and the cloud depth. it is also proportional to the grid size (gsize) .
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            tem = zi (i, ktcon1 (i)) - zi (i, kbcon1 (i))
            dtconv (i) = tem / wc (i)
            tfac = 1. + gsize (i) / 75000.
            dtconv (i) = tfac * dtconv (i)
            dtconv (i) = max (dtconv (i), dtmin)
            dtconv (i) = min (dtconv (i), dtmax)
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! calculate advective time scale (tauadv) using a mean cloud layer wind speed.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            sumx (i) = 0.
            umean (i) = 0.
        endif
    enddo

    do k = 2, km1
        do i = 1, im
            if (cnvflg (i)) then
                if (k >= kbcon1 (i) .and. k < ktcon1 (i)) then
                    dz = zi (i, k) - zi (i, k - 1)
                    tem = sqrt (u1 (i, k) * u1 (i, k) + v1 (i, k) * v1 (i, k))
                    umean (i) = umean (i) + tem * dz
                    sumx (i) = sumx (i) + dz
                endif
            endif
        enddo
    enddo

    do i = 1, im
        if (cnvflg (i)) then
            umean (i) = umean (i) / sumx (i)
            umean (i) = max (umean (i), 1.)
            tauadv (i) = gsize (i) / umean (i)
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! from han et al.'s (2017) equation 6, calculate cloud base mass flux as a function of the mean updraft velcoity for the grid sizes where the quasi - equilibrium assumption of arakawa - schubert is not valid any longer.
    ! as discussed in han et al. (2017), when dtconv is larger than tauadv, the convective mixing is not fully conducted before the cumulus cloud is advected out of the grid cell. in this case, therefore, the cloud base mass flux is further reduced in proportion to the ratio of tauadv to dtconv.
    ! compute cloud base mass flux as a function of the mean
    ! updraft velcoity for the grid sizes where
    ! the quasi - equilibrium assumption of arakawa - schubert is not
    ! valid any longer.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i) .and. .not.asqecflg (i)) then
            k = kbcon (i)
            rho = po (i, k) * 100. / (rdgas * to (i, k))
            tfac = tauadv (i) / dtconv (i)
            tfac = min (tfac, 1.)
            xmb (i) = tfac * betaw_deep * rho * wc (i)
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! compute cloud base mass flux using
    ! the quasi - equilibrium assumption of arakawa - schubert
    ! for the cases where the quasi - equilibrium assumption of arakawa - schubert is valid, first calculate the large scale destabilization as in equation 5 of pan and wu (1995):
    ! \f[
    ! \frac{\partial a}{\partial t}_{ls} = \frac{a^ + - ca^0}{\delta t_{ls}}
    ! \f]
    ! here \f$a^0\f$ is set to zero following han et al.'s (2017), implying that the instability is completely eliminated after the convective adjustment time, \f$\delta t_{ls}\f$.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (asqecflg (i)) then
            ! fld (i) = (aa1 (i) - acrt (i) * acrtfct (i)) / dtconv (i)
            fld (i) = aa1 (i) / dtconv (i)
            if (fld (i) <= 0.) then
                asqecflg (i) = .false.
                cnvflg (i) = .false.
            endif
        endif

        ! -----------------------------------------------------------------------
        ! calculate the stabilization effect of the convection (per unit cloud base mass flux) as in equation 6 of pan and wu (1995):
        ! \f[
        ! \frac{\partial a}{\partial t}_{cu} = \frac{a^ * - a^ + }{\delta t_{cu}}
        ! \f]
        ! \f$\delta t_{cu}\f$ is the short timescale of the convection.
        ! -----------------------------------------------------------------------

        if (asqecflg (i)) then
            ! xaa0 (i) = max (xaa0 (i), 0.)
            xk (i) = (xaa0 (i) - aa1 (i)) / mbdt (i)
            if (xk (i) >= 0.) then
                asqecflg (i) = .false.
                cnvflg (i) = .false.
            endif
        endif

        ! -----------------------------------------------------------------------
        ! kernel, cloud base mass flux
        ! -----------------------------------------------------------------------

        ! -----------------------------------------------------------------------
        ! the cloud base mass flux (xmb) is then calculated from equation 7 of pan and wu (1995)
        ! \f[
        ! m_c = \frac{ - \frac{\partial a}{\partial t}_{ls}}{\frac{\partial a}{\partial t}_{cu}}
        ! \f]
        ! -----------------------------------------------------------------------
        ! again when dtconv is larger than tauadv, the cloud base mass flux is further reduced in proportion to the ratio of tauadv to dtconv.
        ! -----------------------------------------------------------------------

        if (asqecflg (i)) then
            tfac = tauadv (i) / dtconv (i)
            tfac = min (tfac, 1.)
            xmb (i) = - tfac * fld (i) / xk (i)
            ! xmb (i) = min (xmb (i), xmbmax (i))
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! if the large scale destabilization is less than zero, or the stabilization by the convection is greater than zero, then the scheme returns to the calling routine without modifying the state variables.
    ! -----------------------------------------------------------------------

    totflg = .true.
    do i = 1, im
        totflg = totflg .and. (.not. cnvflg (i))
    enddo

    if (totflg) return

    ! -----------------------------------------------------------------------
    ! modified grell & freitas' (2014) updraft fraction which uses
    ! actual entrainment rate at cloud base
    ! for scale - aware parameterization, the updraft fraction (sigmagfm) is first computed as a function of the lateral entrainment rate at cloud base (see han et al.'s (2017) equation 4 and 5), following the study by grell and freitas (2014).
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            tem = min (max (xlamue (i, kbcon (i)), 7.e-5), 3.e-4)
            tem = 0.2 / tem
            tem1 = 3.14 * tem * tem
            sigmagfm (i) = tem1 / (gsize (i) ** 2.0)
            sigmagfm (i) = max (sigmagfm (i), 0.001)
            sigmagfm (i) = min (sigmagfm (i), 0.999)
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! compute scale - aware function based on arakawa & wu (2013)
    ! then, calculate the reduction factor (scaldfunc) of the vertical convective eddy transport of mass flux as a function of updraft fraction from the studies by arakawa and wu (2013) (also see han et al.'s (2017) equation 1 and 2) . the final cloud base mass flux with scale - aware parameterization is obtained from the mass flux when sigmagfm < < 1, multiplied by the reduction factor (han et al.'s (2017) equation 2) .
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            if (gsize (i) < dxcrtuf) then
                scaldfunc (i) = (1. - sigmagfm (i)) * (1. - sigmagfm (i))
                scaldfunc (i) = max (min (scaldfunc (i), 1.0), 0.)
            else
                scaldfunc (i) = 1.0
            endif
            xmb (i) = xmb (i) * scaldfunc (i)
            xmb (i) = min (xmb (i), xmbmax (i))
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! transport aerosols if present
    ! -----------------------------------------------------------------------

    if (do_aerosols) &
        call sa_aamf_deep_aero (im, km, itc, ntc, ntr, delt, &
        xlamde, xlamdd, cnvflg, jmin, kb, kmax, kbcon, ktcon, fscav, &
        edto, xlamd, xmb, c0t, eta, etad, zi, xlamue, xlamud, delp, &
        qtr, qaero)

    ! -----------------------------------------------------------------------
    ! restore to, qo, uo, vo to t1, q1, u1, v1 in case convection stops
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (cnvflg (i) .and. k <= kmax (i)) then
                to (i, k) = t1 (i, k)
                qo (i, k) = q1 (i, k)
                uo (i, k) = u1 (i, k)
                vo (i, k) = v1 (i, k)
                qeso (i, k) = 0.01 * mqs (t1 (i, k)) ! mqs is in pa
                qeso (i, k) = eps * qeso (i, k) / (pfld (i, k) + epsm1 * qeso (i, k))
                val = 1.e-8
                qeso (i, k) = max (qeso (i, k), val)
            endif
        enddo
    enddo

    do n = 1, ntr
        do k = 1, km
            do i = 1, im
                if (cnvflg (i) .and. k <= kmax (i)) then
                    ctro (i, k, n) = ctr (i, k, n)
                endif
            enddo
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! feedback: simply the changes from the cloud with unit mass flux
    ! multiplied by the mass flux necessary to keep the
    ! equilibrium with the larger - scale.
    ! for the "feedback" control, calculate updated values of the state variables by multiplying the cloud base mass flux and the tendencies calculated per unit cloud base mass flux from the static control.
    ! calculate the temperature tendency from the moist static energy and specific humidity tendencies.
    ! update the temperature, specific humidity, and horiztonal wind state variables by multiplying the cloud base mass flux - normalized tendencies by the cloud base mass flux.
    ! accumulate column - integrated tendencies.
    ! -----------------------------------------------------------------------

    do i = 1, im
        delhbar (i) = 0.
        delqbar (i) = 0.
        deltbar (i) = 0.
        delubar (i) = 0.
        delvbar (i) = 0.
        qcond (i) = 0.
    enddo

    do n = 1, ntr
        do i = 1, im
            delebar (i, n) = 0.
        enddo
    enddo

    do k = 1, km
        do i = 1, im
            if (cnvflg (i) .and. k <= kmax (i)) then
                if (k <= ktcon (i)) then
                    dellat = (dellah (i, k) - hlv * dellaq (i, k)) / cp_air
                    t1 (i, k) = t1 (i, k) + dellat * xmb (i) * dt2
                    q1 (i, k) = q1 (i, k) + dellaq (i, k) * xmb (i) * dt2
                    ! tem = 1. / rcs (i)
                    ! u1 (i, k) = u1 (i, k) + dellau (i, k) * xmb (i) * dt2 * tem
                    ! v1 (i, k) = v1 (i, k) + dellav (i, k) * xmb (i) * dt2 * tem
                    u1 (i, k) = u1 (i, k) + dellau (i, k) * xmb (i) * dt2
                    v1 (i, k) = v1 (i, k) + dellav (i, k) * xmb (i) * dt2
                    dp = delp (i, k)
                    delhbar (i) = delhbar (i) + dellah (i, k) * xmb (i) * dp / g
                    delqbar (i) = delqbar (i) + dellaq (i, k) * xmb (i) * dp / g
                    deltbar (i) = deltbar (i) + dellat * xmb (i) * dp / g
                    delubar (i) = delubar (i) + dellau (i, k) * xmb (i) * dp / g
                    delvbar (i) = delvbar (i) + dellav (i, k) * xmb (i) * dp / g
                endif
            endif
        enddo
    enddo

    kk = 0
    do n = 1, ntr + 2
        if (n .eq. ntw .or. n .eq. nti) cycle
        kk = kk + 1
        do k = 1, km
            do i = 1, im
                if (cnvflg (i) .and. k <= kmax (i)) then
                    if (k <= ktcon (i)) then
                        ctr (i, k, kk) = ctr (i, k, kk) + dellae (i, k, kk) * xmb (i) * dt2
                        delebar (i, kk) = delebar (i, kk) + dellae (i, k, kk) * xmb (i) * dp / g
                        qtr (i, k, n) = ctr (i, k, kk)
                    endif
                endif
            enddo
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! recalculate saturation specific humidity using the updated temperature.
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (cnvflg (i) .and. k <= kmax (i)) then
                if (k <= ktcon (i)) then
                    qeso (i, k) = 0.01 * mqs (t1 (i, k)) ! mqs is in pa
                    qeso (i, k) = eps * qeso (i, k) / (pfld (i, k) + epsm1 * qeso (i, k))
                    val = 1.e-8
                    qeso (i, k) = max (qeso (i, k), val)
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! add up column - integrated convective precipitation by multiplying the normalized value by the cloud base mass flux.
    ! -----------------------------------------------------------------------

    do i = 1, im
        rntot (i) = 0.
        delqev (i) = 0.
        delq2 (i) = 0.
        flg (i) = cnvflg (i)
    enddo

    do k = km, 1, - 1
        do i = 1, im
            if (cnvflg (i) .and. k <= kmax (i)) then
                if (k < ktcon (i)) then
                    aup = 1.
                    if (k <= kb (i)) aup = 0.
                    adw = 1.
                    if (k >= jmin (i)) adw = 0.
                    rain = aup * pwo (i, k) + adw * edto (i) * pwdo (i, k)
                    rntot (i) = rntot (i) + rain * xmb (i) * .001 * dt2
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! determine the evaporation of the convective precipitation and update the integrated convective precipitation.
    ! update state temperature and moisture to account for evaporation of convective precipitation.
    ! update column - integrated tendencies to account for evaporation of convective precipitation.
    ! -----------------------------------------------------------------------

    do k = km, 1, - 1
        do i = 1, im
            if (k <= kmax (i)) then
                deltv (i) = 0.
                delq (i) = 0.
                qevap (i) = 0.
                if (cnvflg (i) .and. k < ktcon (i)) then
                    aup = 1.
                    if (k <= kb (i)) aup = 0.
                    adw = 1.
                    if (k >= jmin (i)) adw = 0.
                    rain = aup * pwo (i, k) + adw * edto (i) * pwdo (i, k)
                    rn (i) = rn (i) + rain * xmb (i) * .001 * dt2
                    qr (i, k) = qr (i, k) + rain * xmb (i) * .001 * dt2
                endif
                if (flg (i) .and. k < ktcon (i)) then
                    evef = edt (i) * evfact_deep
                    if (islimsk (i) == 1) evef = edt (i) * evfactl_deep
                    ! if (islimsk (i) == 1) evef = .07
                    ! if (islimsk (i) == 1) evef = 0.
                    qcond (i) = evef * (q1 (i, k) - qeso (i, k)) &
                         / (1. + el2orc * qeso (i, k) / t1 (i, k) ** 2)
                    dp = delp (i, k)
                    if (rn (i) > 0. .and. qcond (i) < 0.) then
                        qevap (i) = - qcond (i) * (1. - exp (- .32 * sqrt (dt2 * rn (i))))
                        qevap (i) = min (qevap (i), rn (i) * 1000. * g / dp)
                        delq2 (i) = delqev (i) + .001 * qevap (i) * dp / g
                    endif
                    if (rn (i) > 0. .and. qcond (i) < 0. .and. delq2 (i) > rntot (i)) then
                        qevap (i) = 1000. * g * (rntot (i) - delqev (i)) / dp
                        flg (i) = .false.
                    endif
                    if (rn (i) > 0. .and. qevap (i) > 0.) then
                        q1 (i, k) = q1 (i, k) + qevap (i)
                        t1 (i, k) = t1 (i, k) - elocp * qevap (i)
                        rn (i) = rn (i) - .001 * qevap (i) * dp / g
                        qr (i, k) = qr (i, k) - .001 * qevap (i) * dp / g
                        deltv (i) = - elocp * qevap (i) / dt2
                        delq (i) = + qevap (i) / dt2
                        delqev (i) = delqev (i) + .001 * dp * qevap (i) / g
                    endif
                    delqbar (i) = delqbar (i) + delq (i) * dp / g
                    deltbar (i) = deltbar (i) + deltv (i) * dp / g
                endif
            endif
        enddo
    enddo

    ! do i = 1, im
    ! if (me == 31 .and. cnvflg (i)) then
    ! if (cnvflg (i)) then
    ! print *, ' deep delhbar, delqbar, deltbar = ', &
    ! delhbar (i), hlv * delqbar (i), cp_air * deltbar (i)
    ! print *, ' deep delubar, delvbar = ', delubar (i), delvbar (i)
    ! print *, ' precip = ', hlv * rn (i) * 1000. / dt2
    ! print *, 'pdif = ', pfld (i, kbcon (i)) - pfld (i, ktcon (i))
    ! endif
    ! enddo

    ! -----------------------------------------------------------------------
    ! precipitation rate converted to actual precip
    ! in unit of m instead of kg
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then

            ! -----------------------------------------------------------------------
            ! in the event of upper level rain evaporation and lower level downdraft
            ! moistening, rn can become negative, in this case, we back out of the
            ! heating and the moistening
            ! -----------------------------------------------------------------------

            if (rn (i) < 0. .and. .not.flg (i)) rn (i) = 0.
            if (rn (i) <= 0.) then
                rn (i) = 0.
            else
                ktop (i) = ktcon (i)
                kbot (i) = kbcon (i)
                kcnv (i) = 1
            endif
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! convective cloud water
    ! calculate convective cloud water.
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (present (cnvw) .and. cnvflg (i) .and. rn (i) > 0.) then
                if (k >= kbcon (i) .and. k < ktcon (i)) then
                    cnvw (i, k) = cnvwt (i, k) * xmb (i) * dt2
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! convective cloud cover
    ! calculate convective cloud cover, which is used when pdf - based cloud fraction is used (i.e., pdfcld = .true.) .
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (present (cnvc) .and. cnvflg (i) .and. rn (i) > 0.) then
                if (k >= kbcon (i) .and. k < ktcon (i)) then
                    cnvc (i, k) = 0.04 * log (1. + 675. * eta (i, k) * xmb (i))
                    cnvc (i, k) = min (cnvc (i, k), 0.6)
                    cnvc (i, k) = max (cnvc (i, k), 0.0)
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! cloud water
    ! separate detrained cloud water into liquid and ice species as a function of temperature only.
    ! -----------------------------------------------------------------------

    if (ncloud > 0) then

        do k = 1, km
            do i = 1, im
                if (cnvflg (i) .and. rn (i) > 0.) then
                    ! if (k > kb (i) .and. k <= ktcon (i)) then
                    if (k >= kbcon (i) .and. k <= ktcon (i)) then
                        tem = dellal (i, k) * xmb (i) * dt2
                        qtr (i, k, ntw) = qtr (i, k, ntw) + tem
                    endif
                endif
            enddo
        enddo

    endif

    ! -----------------------------------------------------------------------
    ! if convective precipitation is zero or negative, reset the updated state variables back to their original values (negating convective changes) .
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (cnvflg (i) .and. rn (i) <= 0.) then
                if (k <= kmax (i)) then
                    t1 (i, k) = to (i, k)
                    q1 (i, k) = qo (i, k)
                    u1 (i, k) = uo (i, k)
                    v1 (i, k) = vo (i, k)
                endif
            endif
        enddo
    enddo

    kk = 0
    do n = 1, ntr + 2
        if (n .eq. ntw .or. n .eq. nti) cycle
        kk = kk + 1
        do k = 1, km
            do i = 1, im
                if (cnvflg (i) .and. rn (i) <= 0.) then
                    if (k <= kmax (i)) then
                        ctr (i, k, kk) = ctro (i, k, kk)
                        qtr (i, k, n) = ctr (i, k, kk)
                    endif
                endif
            enddo
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! store aerosol concentrations if present
    ! -----------------------------------------------------------------------

    if (do_aerosols) then
        do n = 1, ntc
            kk = n + itc - 1
            do k = 1, km
                do i = 1, im
                    if (cnvflg (i) .and. rn (i) > 0.) then
                        if (k <= kmax (i)) qtr (i, k, kk) = qaero (i, k, n)
                    endif
                enddo
            enddo
        enddo
    endif

    ! -----------------------------------------------------------------------
    ! hchuang code change
    ! calculate and retain the updraft and downdraft mass fluxes for dust transport by cumulus convection.
    ! calculate the updraft convective mass flux.
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (present (ud_mf) .and. cnvflg (i) .and. rn (i) > 0.) then
                if (k >= kb (i) .and. k < ktop (i)) then
                    ud_mf (i, k) = eta (i, k) * xmb (i) * dt2
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! save the updraft convective mass flux at cloud top.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (present (dt_mf) .and. present (ud_mf) .and. cnvflg (i) .and. rn (i) > 0.) then
            k = ktop (i) - 1
            dt_mf (i, k) = ud_mf (i, k)
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! calculate the downdraft convective mass flux.
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (present (dd_mf) .and. cnvflg (i) .and. rn (i) > 0.) then
                if (k >= 1 .and. k <= jmin (i)) then
                    dd_mf (i, k) = edto (i) * etad (i, k) * xmb (i) * dt2
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! include tke contribution from deep convection
    ! -----------------------------------------------------------------------

    if (ntk > 0) then

        do k = 2, km1
            do i = 1, im
                if (cnvflg (i) .and. rn (i) > 0.) then
                    if (k > kb (i) .and. k < ktop (i)) then
                        tem = 0.5 * (eta (i, k - 1) + eta (i, k)) * xmb (i)
                        tem1 = pfld (i, k) * 100. / (rdgas * t1 (i, k))
                        sigmagfm (i) = max (sigmagfm (i), betaw_deep)
                        ptem = tem / (sigmagfm (i) * tem1)
                        qtr (i, k, ntk) = qtr (i, k, ntk) + 0.5 * sigmagfm (i) * ptem * ptem
                    endif
                endif
            enddo
        enddo

        do k = 2, km1
            do i = 1, im
                if (cnvflg (i) .and. rn (i) > 0.) then
                    if (k > 1 .and. k <= jmin (i)) then
                        tem = 0.5 * edto (i) * (etad (i, k - 1) + etad (i, k)) * xmb (i)
                        tem1 = pfld (i, k) * 100. / (rdgas * t1 (i, k))
                        sigmagfm (i) = max (sigmagfm (i), betaw_deep)
                        ptem = tem / (sigmagfm (i) * tem1)
                        qtr (i, k, ntk) = qtr (i, k, ntk) + 0.5 * sigmagfm (i) * ptem * ptem
                    endif
                endif
            enddo
        enddo

    endif

end subroutine sa_aamf_deep

! =======================================================================
! Scale-Aware Aerosol-Aware Mass-Flux Shallow Convection
!
! The Scale-Aware Mass-Flux shallow (SAMF_shal) convection scheme is an updated version of the
! previous mass-flux shallow convection scheme with scale and aerosol awareness and parameterizes
! the effect of shallow convection on the environment. The SAMF_shal scheme is similar to the SAMF
! deep convection scheme but with a few key differences.
!
! First, no quasi-equilibrium assumption is used for any grid size and the shallow cloud base mass
! flux is parameterized using a mean updraft velocity. Further, there are no convective downdrafts,
! the entrainment rate is greater than for deep convection, and the shallow convection is limited
! to not extend over the level where \f$p = 0.7p_{sfc}\f$. The paramerization of scale and aerosol
! awareness follows that of the samf deep convection scheme.
!
! The previous version of the shallow convection scheme (shalcnv.f) is described in Han and Pan
! (2011) and differences between the shallow and deep convection schemes are presented in Han and
! Pan (2011) and Han et al. (2017). Details of scale- and aerosol-aware parameterizations are
! described in Han et al. (2017).
!
! In further update for FY19 GFS implementation, interaction with Turbulent Kinetic Energy (TKE),
! which is a prognostic variable used in a scale-aware tke-based moist EDMF vertical turbulent
! mixing scheme, is included. Entrainment rates in updrafts are proportional to sub-cloud mean TKE.
! TKE is transported by cumulus convection. TKE contribution from cumulus convection is deduced
! from cumulus mass flux. On the other hand, tracers such as ozone and aerosol are also transported
! by cumulus convection.
!
! To reduce too much convective cooling at the cloud top, the convection schemes have been modified
! for the rain conversion rate, entrainment and detrainment rates, overshooting layers, and maximum
! allowable cloudbase mass flux (as of JUNE 2018) .
!
! contains the entire samf shallow convection scheme.
!
! This routine follows the SAMF deep scheme quite closely, although it can be interpreted as only
! having the "static" and "feedback" control portions, since the "dynamic" control is not necessary
! to find the cloud base mass flux. The algorithm is simplified from SAMF deep convection by
! excluding convective downdrafts and being confined to operate below \f$p = 0.7p_{sfc}\f$. Also,
! entrainment is both simpler and stronger in magnitude compared to the deep scheme.
!
! \param[in] IM      number of used points
! \param[in] KM      vertical layer dimension
! \param[in] DELT    physics time step in seconds
! \param[in] NTK     index for tke
! \param[in] NTR     total number of tracers including tke
! \param[in] DELP    pressure difference between level k and k + 1 (pa)
! \param[in] PRSLP   mean layer presure (pa)
! \param[in] PSP     surface pressure (pa)
! \param[in] PHIL    layer geopotential (\f$m^s / s^2\f$)
! \param[in] QTR     tracer array including cloud condensate (\f$kg / kg\f$)
! \param[inout] QL   cloud water or ice (kg / kg)
! \param[inout] Q1   updated tracers (kg / kg)
! \param[inout] T1   updated temperature (k)
! \param[inout] U1   updated zonal wind (\f$m s^{ - 1}\f$)
! \param[inout] V1   updated meridional wind (\f$m s^{ - 1}\f$)
! \param[out] RN     convective rain (m)
! \param[out] KBOT   index for cloud base
! \param[out] KTOP   index for cloud top
! \param[out] KCNV   flag to denote deep convection (0 = no, 1 = yes)
! \param[in] ISLIMSK sea / land / ice mask (= 0 / 1 / 2)
! \param[in] DOT     layer mean vertical velocity (pa / s)
! \param[in] NCLOUD  number of cloud species
! \param[in] HPBL    pbl height (m)
! \param[in] HEAT    surface sensible heat flux (k m / s)
! \param[in] EVAP    surface latent heat flux (kg / kg m / s)
! \param[out] UD_MF  updraft mass flux multiplied by time step (\f$kg / m^2\f$)
! \param[out] DT_MF  ud_mf at cloud top (\f$kg / m^2\f$)
! \param[out] CNVW   convective cloud water (kg / kg)
! \param[out] CNVC   convective cloud cover (unitless)
!
! General Algorithm
! # Compute preliminary quantities needed for the static and feedback control portions of the algorithm.
! # Perform calculations related to the updraft of the entraining / detraining cloud model ("static control") .
! # The cloud base mass flux is obtained using the cumulus updraft velocity averaged ove the whole cloud depth.
! # Calculate the tendencies of the state variables (per unit cloud base mass flux) and the cloud base mass flux.
! # For the "feedback control", calculate updated values of the state variables by multiplying the cloud base
!   mass flux and the tendencies calculated per unit cloud base mass flux from the static control.
! =======================================================================

subroutine sa_aamf_shal (im, km, delt, itc, ntc, ntw, nti, ntk, ntr, delp, &
        prslp, psp, phil, qtr, q1, t1, u1, v1, qr, fscav, rn, kbot, ktop, &
        kcnv, islimsk, gsize, dot, ncloud, hpbl, ud_mf, dt_mf, cnvw, cnvc)

    implicit none

    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------

    integer, intent (in) :: im, km, itc, ntc, ntw, nti, ntk, ntr, ncloud, islimsk (im)

    real, intent (in) :: delt
    real, intent (in) :: psp (im), delp (im, km), &
        prslp (im, km), gsize (im), hpbl (im), dot (im, km), phil (im, km)
    real, intent (in) :: fscav (ntc)

    integer, intent (inout) :: kbot (im), ktop (im), kcnv (im)

    real, intent (inout) :: qtr (im, km, ntr + 2), q1 (im, km), t1 (im, km), &
        u1 (im, km), v1 (im, km)

    real, intent (out) :: rn (im), qr (im, km)
    real, intent (out), optional :: cnvw (im, km), cnvc (im, km), &
        ! hchuang code change mass flux output
        ud_mf (im, km), dt_mf (im, km)

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: i, j, indx, k, kk, km1, n
    integer :: kpbl (im)

    real :: clamd, tkemx, tkemn, dtke

    real :: dellat, delta, &
        c0l, d0, &
        desdt, dp, &
        dq, dqsdp, dqsdt, dt, &
        dt2, dtmax, dtmin, &
        dv1h, dv2h, dv3h, &
        dv1q, dv2q, dv3q, &
        dz, dz1, e1, &
        el2orc, elocp, aafac, cm, &
        es, etah, h1, &
        evef, fact1, &
        fact2, factor, dthk, &
        g, gamma, pprime, &
        qlk, qrch, qs, &
        rfact, shear, tfac, &
        val, val1, val2, &
        w1, w1l, w1s, w2, &
        w2l, w2s, w3, w3l, &
        w3s, w4, w4l, w4s, &
        rho, tem, tem1, tem2, &
        ptem, ptem1

    integer :: kb (im), kbcon (im), kbcon1 (im), &
        ktcon (im), ktcon1 (im), ktconn (im), &
        kbm (im), kmax (im)

    real :: aa1 (im), cina (im), &
        tkemean (im), clamt (im), &
        umean (im), tauadv (im), &
        delhbar (im), delq (im), delq2 (im), &
        delqbar (im), delqev (im), deltbar (im), &
        deltv (im), dtconv (im), edt (im), &
        pdot (im), po (im, km), &
        qcond (im), qevap (im), hmax (im), &
        rntot (im), vshear (im), &
        xlamud (im), xmb (im), xmbmax (im), &
        delebar (im, ntr), &
        delubar (im), delvbar (im)

    real :: c0 (im)

    real :: crtlamd

    real :: cinpcr, cinpcrmx, cinpcrmn, &
        cinacr, cinacrmx, cinacrmn

    ! parameters for updraft velocity calculation
    real :: bet1, cd1, f1, gam1, &
        bb1, bb2, tkcrt, cmxfac, csmf

    ! physical parameters
    parameter (g = grav)
    parameter (elocp = hlv / cp_air, &
        el2orc = hlv * hlv / (rvgas * cp_air))
    parameter (d0 = .001)

    ! asolfac_shal: aerosol - aware parameter based on lim & hong (2012)
    ! asolfac_shal = cx / c0s_shal (= .002)
    ! cx = min ([ - 0.7 ln (nccn) + 24] * 1.e-4, c0s_shal)
    ! nccn: ccn number concentration in cm^ (- 3)
    ! until a realistic nccn is provided, typical nccns are assumed
    ! as nccn = 100 for sea and nccn = 1000 for land

    parameter (cm = 1.0, delta = zvir)
    parameter (fact1 = (cp_vap - c_liq) / rvgas, fact2 = hlv / rvgas - fact1 * tice)
    parameter (clamd = 0.1, tkemx = 0.65, tkemn = 0.05)
    parameter (dtke = tkemx - tkemn)
    parameter (dthk = 25.)
    parameter (cinpcrmx = 180., cinpcrmn = 120.)
    parameter (cinacrmx = - 120., cinacrmn = - 80.)
    parameter (crtlamd = 3.e-4)
    parameter (dtmax = 10800., dtmin = 600.)
    parameter (bet1 = 1.875, cd1 = .506, f1 = 2.0, gam1 = .5)
    parameter (h1 = 0.33333333)
    parameter (bb1 = 4.0, bb2 = 0.8, csmf = 0.2)
    parameter (tkcrt = 2., cmxfac = 15.)

    ! local variables and arrays
    real :: pfld (im, km), to (im, km), qo (im, km), &
        uo (im, km), vo (im, km), qeso (im, km), &
        ctr (im, km, ntr), ctro (im, km, ntr)

    ! for aerosol transport
    real :: qaero (im, km, ntc)

    ! for updraft velocity calculation
    real :: wu2 (im, km), buo (im, km), drag (im, km), wush (im, km)
    real :: wc (im), scaldfunc (im), sigmagfm (im)

    ! cloud water
    ! real :: tvo (im, km),
    real :: qlko_ktcon (im), dellal (im, km), &
        dbyo (im, km), zo (im, km), xlamue (im, km), &
        heo (im, km), heso (im, km), &
        dellah (im, km), dellaq (im, km), &
        dellae (im, km, ntr), &
        dellau (im, km), dellav (im, km), hcko (im, km), &
        ucko (im, km), vcko (im, km), qcko (im, km), &
        qrcko (im, km), eta (im, km), &
        ecko (im, km, ntr), &
        zi (im, km), pwo (im, km), c0t (im, km), &
        sumx (im), tx1 (im), cnvwt (im, km)

    logical :: do_aerosols, totflg, cnvflg (im), flg (im)

    real :: tf, tcr, tcrf
    parameter (tf = 233.16, tcr = 263.16, tcrf = 1.0 / (tcr - tf))

    ! -----------------------------------------------------------------------
    ! determine whether to perform aerosol transport
    ! -----------------------------------------------------------------------

    do_aerosols = (itc > 0) .and. (ntc > 0) .and. (ntr > 0)
    if (do_aerosols) do_aerosols = (ntr >= itc + ntc - 3)

    ! -----------------------------------------------------------------------
    ! convert input pa terms to cb terms -- moorthi
    ! compute preliminary quantities needed for the static and feedback control portions of the algorithm.
    ! convert input pressure terms to centibar units.
    ! -----------------------------------------------------------------------

    km1 = km - 1

    ! -----------------------------------------------------------------------
    ! initialize arrays
    ! initialize column - integrated and other single - value - per - column variable arrays.
    ! -----------------------------------------------------------------------

    do i = 1, im
        cnvflg (i) = .true.
        if (kcnv (i) == 1) cnvflg (i) = .false.
        if (cnvflg (i)) then
            kbot (i) = km + 1
            ktop (i) = 0
        endif
        rn (i) = 0.
        kbcon (i) = km
        ktcon (i) = 1
        ktconn (i) = 1
        kb (i) = km
        pdot (i) = 0.
        qlko_ktcon (i) = 0.
        edt (i) = 0.
        aa1 (i) = 0.
        cina (i) = 0.
        vshear (i) = 0.
    enddo

    ! -----------------------------------------------------------------------
    ! return to the calling routine if deep convection is present or the surface buoyancy flux is negative.
    ! -----------------------------------------------------------------------

    totflg = .true.
    do i = 1, im
        totflg = totflg .and. (.not. cnvflg (i))
    enddo

    if (totflg) return

    ! -----------------------------------------------------------------------
    ! determine aerosol - aware rain conversion parameter over land
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (islimsk (i) == 1) then
            c0 (i) = c0s_shal * asolfac_shal
        else
            c0 (i) = c0s_shal
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! determine rain conversion parameter above the freezing level which exponentially decreases with decreasing temperature from han et al.'s (2017) equation 8.
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (t1 (i, k) > 273.16) then
                c0t (i, k) = c0 (i)
            else
                tem = d0 * (t1 (i, k) - 273.16)
                tem1 = exp (tem)
                c0t (i, k) = c0 (i) * tem1
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! initialize convective cloud water and cloud cover to zero.
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (present (cnvw)) cnvw (i, k) = 0.
            if (present (cnvc)) cnvc (i, k) = 0.
        enddo
    enddo

    do k = 1, km
        do i = 1, im
            qr (i, k) = 0.
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! hchuang code change
    ! initialize updraft mass fluxes to zero.
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (present (ud_mf)) ud_mf (i, k) = 0.
            if (present (dt_mf)) dt_mf (i, k) = 0.
        enddo
    enddo

    dt2 = delt

    ! -----------------------------------------------------------------------
    ! model tunable parameters are all here
    ! -----------------------------------------------------------------------

    ! clam_shal = .3
    aafac = .05
    ! evef = 0.07
    ! evfact_shal = 0.3
    ! evfactl_shal = 0.3

    ! pgcon_shal = 0.7 ! gregory et al. (1997, qjrms)
    ! pgcon_shal = 0.55 ! zhang & wu (2003, jas)

    w1l = - 8.e-3
    w2l = - 4.e-2
    w3l = - 5.e-3
    w4l = - 5.e-4
    w1s = - 2.e-4
    w2s = - 2.e-3
    w3s = - 1.e-3
    w4s = - 2.e-5

    ! -----------------------------------------------------------------------
    ! define top layer for search of the downdraft originating layer
    ! and the maximum thetae for updraft
    ! determine maximum indices for the parcel starting point (kbm) and cloud top (kmax) .
    ! -----------------------------------------------------------------------

    do i = 1, im
        kbm (i) = km
        kmax (i) = km
        tx1 (i) = 1.0 / psp (i)
    enddo

    do k = 1, km
        do i = 1, im
            if (prslp (i, k) * tx1 (i) > 0.70) kbm (i) = k + 1
            if (prslp (i, k) * tx1 (i) > 0.60) kmax (i) = k + 1
        enddo
    enddo

    do i = 1, im
        kbm (i) = min (kbm (i), kmax (i))
    enddo

    ! -----------------------------------------------------------------------
    ! hydrostatic height assume zero terr and compute
    ! updraft entrainment rate as an inverse function of height
    ! calculate hydrostatic height at layer centers assuming a flat surface (no terrain) from the geopotential.
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            zo (i, k) = phil (i, k) / g
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! calculate interface height
    ! -----------------------------------------------------------------------

    do k = 1, km1
        do i = 1, im
            zi (i, k) = 0.5 * (zo (i, k) + zo (i, k + 1))
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! pbl height
    ! find the index for the pbl top using the pbl height; enforce that it is lower than the maximum parcel starting level.
    ! -----------------------------------------------------------------------

    do i = 1, im
        flg (i) = cnvflg (i)
        kpbl (i) = 1
    enddo

    do k = 2, km1
        do i = 1, im
            if (flg (i) .and. zo (i, k) <= hpbl (i)) then
                kpbl (i) = k
            else
                flg (i) = .false.
            endif
        enddo
    enddo

    do i = 1, im
        kpbl (i) = min (kpbl (i), kbm (i))
    enddo

    ! -----------------------------------------------------------------------
    ! convert surface pressure to mb from cb
    ! convert prslp from centibar to millibar, set normalized mass flux to 1, cloud properties to 0, and save model state variables (after advection / turbulence) .
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (cnvflg (i) .and. k <= kmax (i)) then
                pfld (i, k) = prslp (i, k) * 0.01
                eta (i, k) = 1.
                hcko (i, k) = 0.
                qcko (i, k) = 0.
                qrcko (i, k) = 0.
                ucko (i, k) = 0.
                vcko (i, k) = 0.
                dbyo (i, k) = 0.
                pwo (i, k) = 0.
                dellal (i, k) = 0.
                to (i, k) = t1 (i, k)
                qo (i, k) = q1 (i, k)
                uo (i, k) = u1 (i, k)
                vo (i, k) = v1 (i, k)
                ! uo (i, k) = u1 (i, k) * rcs (i)
                ! vo (i, k) = v1 (i, k) * rcs (i)
                wu2 (i, k) = 0.
                buo (i, k) = 0.
                drag (i, k) = 0.
                cnvwt (i, k) = 0.
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! initialize tracer variables
    ! -----------------------------------------------------------------------

    kk = 0
    do n = 1, ntr + 2
        if (n .eq. ntw .or. n .eq. nti) cycle
        kk = kk + 1
        do k = 1, km
            do i = 1, im
                if (cnvflg (i) .and. k <= kmax (i)) then
                    ctr (i, k, kk) = qtr (i, k, n)
                    ctro (i, k, kk) = qtr (i, k, n)
                    ecko (i, k, kk) = 0.
                endif
            enddo
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! column variables
    ! p is pressure of the layer (mb)
    ! t is temperature at t - dt (k) ..tn
    ! q is mixing ratio at t - dt (kg / kg) ..qn
    ! to is temperature at t + dt (k) ... this is after advection and turbulan
    ! qo is mixing ratio at t + dt (kg / kg) ..q1
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! calculate saturation specific humidity and enforce minimum moisture values.
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (cnvflg (i) .and. k <= kmax (i)) then
                qeso (i, k) = 0.01 * mqs (to (i, k)) ! mqs is in pa
                qeso (i, k) = eps * qeso (i, k) / (pfld (i, k) + epsm1 * qeso (i, k))
                val1 = 1.e-8
                qeso (i, k) = max (qeso (i, k), val1)
                val2 = 1.e-10
                qo (i, k) = max (qo (i, k), val2)
                ! qo (i, k) = min (qo (i, k), qeso (i, k))
                ! tvo (i, k) = to (i, k) + delta * to (i, k) * qo (i, k)
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute moist static energy
    ! calculate moist static energy (heo) and saturation moist static energy (heso) .
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (cnvflg (i) .and. k <= kmax (i)) then
                ! tem = g * zo (i, k) + cp_air * to (i, k)
                tem = phil (i, k) + cp_air * to (i, k)
                heo (i, k) = tem + hlv * qo (i, k)
                heso (i, k) = tem + hlv * qeso (i, k)
                ! heo (i, k) = min (heo (i, k), heso (i, k))
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! determine level with largest moist static energy within pbl
    ! this is the level where updraft starts
    ! perform calculations related to the updraft of the entraining / detraining cloud model ("static control") .
    ! search in the pbl for the level of maximum moist static energy to start the ascending parcel.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            hmax (i) = heo (i, 1)
            kb (i) = 1
        endif
    enddo

    do k = 2, km
        do i = 1, im
            if (cnvflg (i) .and. k <= kpbl (i)) then
                if (heo (i, k) > hmax (i)) then
                    kb (i) = k
                    hmax (i) = heo (i, k)
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! calculate the temperature, water vapor mixing ratio, and pressure at interface levels.
    ! -----------------------------------------------------------------------

    do k = 1, km1
        do i = 1, im
            if (cnvflg (i) .and. k <= kmax (i) - 1) then
                dz = .5 * (zo (i, k + 1) - zo (i, k))
                dp = .5 * (pfld (i, k + 1) - pfld (i, k))
                es = 0.01 * mqs (to (i, k + 1)) ! mqs is in pa
                pprime = pfld (i, k + 1) + epsm1 * es
                qs = eps * es / pprime
                dqsdp = - qs / pprime
                desdt = es * (fact1 / to (i, k + 1) + fact2 / (to (i, k + 1) ** 2))
                dqsdt = qs * pfld (i, k + 1) * desdt / (es * pprime)
                gamma = el2orc * qeso (i, k + 1) / (to (i, k + 1) ** 2)
                dt = (g * dz + hlv * dqsdp * dp) / (cp_air * (1. + gamma))
                dq = dqsdt * dt + dqsdp * dp
                to (i, k) = to (i, k + 1) + dt
                qo (i, k) = qo (i, k + 1) + dq
                po (i, k) = .5 * (pfld (i, k) + pfld (i, k + 1))
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! recalculate saturation specific humidity, moist static energy, saturation moist static energy, and horizontal momentum on interface levels. enforce minimum specific humidity.
    ! -----------------------------------------------------------------------

    do k = 1, km1
        do i = 1, im
            if (cnvflg (i) .and. k <= kmax (i) - 1) then
                qeso (i, k) = 0.01 * mqs (to (i, k)) ! mqs is in pa
                qeso (i, k) = eps * qeso (i, k) / (po (i, k) + epsm1 * qeso (i, k))
                val1 = 1.e-8
                qeso (i, k) = max (qeso (i, k), val1)
                val2 = 1.e-10
                qo (i, k) = max (qo (i, k), val2)
                ! qo (i, k) = min (qo (i, k), qeso (i, k))
                heo (i, k) = .5 * g * (zo (i, k) + zo (i, k + 1)) + &
                    cp_air * to (i, k) + hlv * qo (i, k)
                heso (i, k) = .5 * g * (zo (i, k) + zo (i, k + 1)) + &
                    cp_air * to (i, k) + hlv * qeso (i, k)
                uo (i, k) = .5 * (uo (i, k) + uo (i, k + 1))
                vo (i, k) = .5 * (vo (i, k) + vo (i, k + 1))
            endif
        enddo
    enddo

    do n = 1, ntr
        do k = 1, km1
            do i = 1, im
                if (cnvflg (i) .and. k <= kmax (i) - 1) then
                    ctro (i, k, n) = .5 * (ctro (i, k, n) + ctro (i, k + 1, n))
                endif
            enddo
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! look for the level of free convection as cloud base
    ! search below the index "kbm" for the level of free convection (lfc) where the condition \f$h_b > h^ * \f$ is first met, where \f$h_b, h^ * \f$ are the state moist static energy at the parcel's starting level and saturation moist static energy, respectively. set "kbcon" to the index of the lfc.
    ! -----------------------------------------------------------------------

    do i = 1, im
        flg (i) = cnvflg (i)
        if (flg (i)) kbcon (i) = kmax (i)
    enddo

    do k = 2, km1
        do i = 1, im
            if (flg (i) .and. k < kbm (i)) then
                if (k > kb (i) .and. heo (i, kb (i)) > heso (i, k)) then
                    kbcon (i) = k
                    flg (i) = .false.
                endif
            endif
        enddo
    enddo

    do i = 1, im
        if (cnvflg (i)) then
            if (kbcon (i) == kmax (i)) cnvflg (i) = .false.
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! if no lfc, return to the calling routine without modifying state variables.
    ! -----------------------------------------------------------------------

    totflg = .true.
    do i = 1, im
        totflg = totflg .and. (.not. cnvflg (i))
    enddo

    if (totflg) return

    ! -----------------------------------------------------------------------
    ! determine the vertical pressure velocity at the lfc. after han and pan (2011), determine the maximum pressure thickness between a parcel's starting level and the lfc. if a parcel doesn't reach the lfc within the critical thickness, then the convective inhibition is deemed too great for convection to be triggered, and the subroutine returns to the calling routine without modifying the state variables.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            ! pdot (i) = 10. * dot (i, kbcon (i))
            pdot (i) = 0.01 * dot (i, kbcon (i)) ! now dot is in pa / s
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! turn off convection if pressure depth between parcel source level
    ! and cloud base is larger than a critical value, cinpcr
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            if (islimsk (i) == 1) then
                w1 = w1l
                w2 = w2l
                w3 = w3l
                w4 = w4l
            else
                w1 = w1s
                w2 = w2s
                w3 = w3s
                w4 = w4s
            endif
            if (pdot (i) <= w4) then
                tem = (pdot (i) - w4) / (w3 - w4)
            elseif (pdot (i) >= - w4) then
                tem = - (pdot (i) + w4) / (w4 - w3)
            else
                tem = 0.
            endif
            val1 = - 1.
            tem = max (tem, val1)
            val2 = 1.
            tem = min (tem, val2)
            ptem = 1. - tem
            ptem1 = .5 * (cinpcrmx - cinpcrmn)
            cinpcr = cinpcrmx - ptem * ptem1
            tem1 = pfld (i, kb (i)) - pfld (i, kbcon (i))
            if (tem1 > cinpcr) then
                cnvflg (i) = .false.
            endif
        endif
    enddo

    totflg = .true.
    do i = 1, im
        totflg = totflg .and. (.not. cnvflg (i))
    enddo

    if (totflg) return

    ! -----------------------------------------------------------------------
    ! turbulent entrainment rate assumed to be proportional
    ! to subcloud mean tke
    ! -----------------------------------------------------------------------

    if (ntk > 0) then
        do i = 1, im
            if (cnvflg (i)) then
                sumx (i) = 0.
                tkemean (i) = 0.
            endif
        enddo
        do k = 1, km1
            do i = 1, im
                if (cnvflg (i)) then
                    if (k >= kb (i) .and. k < kbcon (i)) then
                        dz = zo (i, k + 1) - zo (i, k)
                        tem = 0.5 * (qtr (i, k, ntk) + qtr (i, k + 1, ntk))
                        tkemean (i) = tkemean (i) + tem * dz
                        sumx (i) = sumx (i) + dz
                    endif
                endif
            enddo
        enddo
        do i = 1, im
            if (cnvflg (i)) then
                tkemean (i) = tkemean (i) / sumx (i)
                if (tkemean (i) > tkemx) then
                    clamt (i) = clam_shal + clamd
                else if (tkemean (i) < tkemn) then
                    clamt (i) = clam_shal - clamd
                else
                    tem = tkemx - tkemean (i)
                    tem1 = 1. - 2. * tem / dtke
                    clamt (i) = clam_shal + clamd * tem1
                endif
            endif
        enddo
        ! kgao 12 / 18 / 2023 - adjust entrainment rate based on tke
        if (use_tke_conv) then
            do i = 1, im
                if (cnvflg (i)) then
                    if (tkemean (i) > tkcrt) then
                        tem = 1. + tkemean (i) / tkcrt
                        tem1 = min (tem, cmxfac)
                        clamt (i) = tem1 * clam_shal
                    endif
                endif
            enddo
        endif
    else
        do i = 1, im
            if (cnvflg (i)) then
                clamt (i) = clam_shal
            endif
        enddo
    endif

    ! -----------------------------------------------------------------------
    ! assume updraft entrainment rate
    ! is an inverse function of height
    ! -----------------------------------------------------------------------

    do k = 1, km1
        do i = 1, im
            if (cnvflg (i)) then
                xlamue (i, k) = clamt (i) / zi (i, k)
            endif
        enddo
    enddo

    do i = 1, im
        if (cnvflg (i)) then
            xlamue (i, km) = xlamue (i, km1)
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! specify the detrainment rate for the updrafts
    ! (the updraft detrainment rate is set constant and equal to the entrainment rate at cloud base.)
    ! the updraft detrainment rate is vertically constant and proportional to clamt
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            ! xlamud (i) = xlamue (i, kbcon (i))
            ! xlamud (i) = crtlamd
            xlamud (i) = 0.001 * clamt (i)
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! determine updraft mass flux for the subcloud layers
    ! calculate the normalized mass flux for subcloud and in - cloud layers according to pan and wu (1995) equation 1:
    ! \f[
    ! \frac{1}{\eta}\frac{\partial \eta}{\partial z} = \lambda_e - \lambda_d
    ! \f]
    ! where \f$\eta\f$ is the normalized mass flux, \f$\lambda_e\f$ is the entrainment rate and \f$\lambda_d\f$ is the detrainment rate. the normalized mass flux increases upward below the cloud base and decreases upward above.
    ! -----------------------------------------------------------------------

    do k = km1, 1, - 1
        do i = 1, im
            if (cnvflg (i)) then
                if (k < kbcon (i) .and. k >= kb (i)) then
                    dz = zi (i, k + 1) - zi (i, k)
                    ptem = 0.5 * (xlamue (i, k) + xlamue (i, k + 1)) - xlamud (i)
                    eta (i, k) = eta (i, k + 1) / (1. + ptem * dz)
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute mass flux above cloud base
    ! -----------------------------------------------------------------------

    do i = 1, im
        flg (i) = cnvflg (i)
    enddo

    do k = 2, km1
        do i = 1, im
            if (flg (i)) then
                if (k > kbcon (i) .and. k < kmax (i)) then
                    dz = zi (i, k) - zi (i, k - 1)
                    ptem = 0.5 * (xlamue (i, k) + xlamue (i, k - 1)) - xlamud (i)
                    eta (i, k) = eta (i, k - 1) * (1 + ptem * dz)
                    if (eta (i, k) <= 0.) then
                        kmax (i) = k
                        ktconn (i) = k
                        kbm (i) = min (kbm (i), kmax (i))
                        flg (i) = .false.
                    endif
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute updraft cloud property
    ! set cloud properties equal to the state variables at updraft starting level (kb) .
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            indx = kb (i)
            hcko (i, indx) = heo (i, indx)
            ucko (i, indx) = uo (i, indx)
            vcko (i, indx) = vo (i, indx)
        endif
    enddo

    ! for tracers
    do n = 1, ntr
        do i = 1, im
            if (cnvflg (i)) then
                indx = kb (i)
                ecko (i, indx, n) = ctro (i, indx, n)
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! cm is an enhancement factor in entrainment rates for momentum
    ! calculate the cloud properties as a parcel ascends, modified by entrainment and detrainment. discretization follows appendix b of grell (1993). following han and pan (2006), the convective momentum transport is reduced by the convection - induced pressure gradient force by the constant "pgcon_shal", currently set to 0.55 after zhang and wu (2003).
    ! -----------------------------------------------------------------------

    do k = 2, km1
        do i = 1, im
            if (cnvflg (i)) then
                if (k > kb (i) .and. k < kmax (i)) then
                    dz = zi (i, k) - zi (i, k - 1)
                    tem = 0.5 * (xlamue (i, k) + xlamue (i, k - 1)) * dz
                    tem1 = 0.5 * xlamud (i) * dz
                    factor = 1. + tem - tem1
                    hcko (i, k) = ((1. - tem1) * hcko (i, k - 1) + tem * 0.5 * &
                         (heo (i, k) + heo (i, k - 1))) / factor
                    dbyo (i, k) = hcko (i, k) - heso (i, k)

                    tem = 0.5 * cm * tem
                    factor = 1. + tem
                    ptem = tem + pgcon_shal
                    ptem1 = tem - pgcon_shal
                    ucko (i, k) = ((1. - tem) * ucko (i, k - 1) + ptem * uo (i, k) &
                         + ptem1 * uo (i, k - 1)) / factor
                    vcko (i, k) = ((1. - tem) * vcko (i, k - 1) + ptem * vo (i, k) &
                         + ptem1 * vo (i, k - 1)) / factor
                endif
            endif
        enddo
    enddo

    do n = 1, ntr
        do k = 2, km1
            do i = 1, im
                if (cnvflg (i)) then
                    if (k > kb (i) .and. k < kmax (i)) then
                        dz = zi (i, k) - zi (i, k - 1)
                        tem = 0.25 * (xlamue (i, k) + xlamue (i, k - 1)) * dz
                        factor = 1. + tem
                        ecko (i, k, n) = ((1. - tem) * ecko (i, k - 1, n) + tem * &
                            (ctro (i, k, n) + ctro (i, k - 1, n))) / factor
                    endif
                endif
            enddo
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! taking account into convection inhibition due to existence of
    ! dry layers below cloud base
    ! with entrainment, recalculate the lfc as the first level where buoyancy is positive. the difference in pressure levels between lfcs calculated with / without entrainment must be less than a threshold (currently 25 hpa) . otherwise, convection is inhibited and the scheme returns to the calling routine without modifying the state variables. this is the subcloud dryness trigger modification discussed in han and pan (2011).
    ! -----------------------------------------------------------------------

    do i = 1, im
        flg (i) = cnvflg (i)
        kbcon1 (i) = kmax (i)
    enddo

    do k = 2, km1
        do i = 1, im
            if (flg (i) .and. k < kbm (i)) then
                if (k >= kbcon (i) .and. dbyo (i, k) > 0.) then
                    kbcon1 (i) = k
                    flg (i) = .false.
                endif
            endif
        enddo
    enddo

    do i = 1, im
        if (cnvflg (i)) then
            if (kbcon1 (i) == kmax (i)) cnvflg (i) = .false.
        endif
    enddo

    do i = 1, im
        if (cnvflg (i)) then
            tem = pfld (i, kbcon (i)) - pfld (i, kbcon1 (i))
            if (tem > dthk) then
                cnvflg (i) = .false.
            endif
        endif
    enddo

    totflg = .true.
    do i = 1, im
        totflg = totflg .and. (.not. cnvflg (i))
    enddo

    if (totflg) return

    ! -----------------------------------------------------------------------
    ! calculate convective inhibition
    ! calculate additional trigger condition of the convective inhibition (cin) according to han et al.'s (2017) equation 13.
    ! -----------------------------------------------------------------------

    do k = 2, km1
        do i = 1, im
            if (cnvflg (i)) then
                if (k > kb (i) .and. k < kbcon1 (i)) then
                    dz1 = zo (i, k + 1) - zo (i, k)
                    gamma = el2orc * qeso (i, k) / (to (i, k) ** 2)
                    rfact = 1. + delta * cp_air * gamma &
                         * to (i, k) / hlv
                    cina (i) = cina (i) + &
                        ! dz1 * eta (i, k) * (g / (cp_air * to (i, k))) &
                        dz1 * (g / (cp_air * to (i, k))) &
                         * dbyo (i, k) / (1. + gamma) &
                         * rfact
                    val = 0.
                    cina (i) = cina (i) + &
                        ! dz1 * eta (i, k) * g * delta * &
                        dz1 * g * delta * &
                        max (val, (qeso (i, k) - qo (i, k)))
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! turn off convection if the cin is less than a critical value (cinacr) which is inversely proportional to the large - scale vertical velocity.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then

            if (islimsk (i) == 1) then
                w1 = w1l
                w2 = w2l
                w3 = w3l
                w4 = w4l
            else
                w1 = w1s
                w2 = w2s
                w3 = w3s
                w4 = w4s
            endif
            if (pdot (i) <= w4) then
                tem = (pdot (i) - w4) / (w3 - w4)
            elseif (pdot (i) >= - w4) then
                tem = - (pdot (i) + w4) / (w4 - w3)
            else
                tem = 0.
            endif

            val1 = - 1.
            tem = max (tem, val1)
            val2 = 1.
            tem = min (tem, val2)
            tem = 1. - tem
            tem1 = .5 * (cinacrmx - cinacrmn)
            cinacr = cinacrmx - tem * tem1

            ! cinacr = cinacrmx
            if (cina (i) < cinacr) cnvflg (i) = .false.
        endif
    enddo

    totflg = .true.
    do i = 1, im
        totflg = totflg .and. (.not. cnvflg (i))
    enddo

    if (totflg) return

    ! -----------------------------------------------------------------------
    ! determine first guess cloud top as the level of zero buoyancy
    ! limited to the level of p / ps = 0.7
    ! calculate the cloud top as the first level where parcel buoyancy becomes negative; the maximum possible value is at \f$p = 0.7p_{sfc}\f$.
    ! -----------------------------------------------------------------------

    do i = 1, im
        flg (i) = cnvflg (i)
        if (flg (i)) ktcon (i) = kbm (i)
    enddo

    do k = 2, km1
        do i = 1, im
            if (flg (i) .and. k < kbm (i)) then
                if (k > kbcon1 (i) .and. dbyo (i, k) < 0.) then
                    ktcon (i) = k
                    flg (i) = .false.
                endif
            endif
        enddo
    enddo

    ! kg change: turn off shal conv based on diagnosed cloud depth or top
    ! the idea here is that if the cloud is too deep or too high, it should not be
    ! handled by shal conv
    do i = 1, im
        if (cnvflg (i) .and. limit_shal_conv) then
            ! a) cloud depth criterion as in deep conv
            tem = pfld (i, kbcon (i)) - pfld (i, ktcon (i))
            if (tem >= cthk_shal) cnvflg (i) = .false.
            ! b) cloud top criterion
            if (prslp (i, ktcon (i)) * tx1 (i) < top_shal) cnvflg (i) = .false.
            !if (ktcon (i) > kmax (i)) cnvflg (i) = .false.
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! specify upper limit of mass flux at cloud base
    ! calculate the maximum value of the cloud base mass flux using the cfl - criterion - based formula of han and pan (2011), equation 7.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            ! xmbmax (i) = .1
            !
            k = kbcon (i)
            dp = delp (i, k)
            xmbmax (i) = dp / (2. * g * dt2)
            ! xmbmax (i) = dp / (g * dt2)
            !
            ! tem = dp / (g * dt2)
            ! xmbmax (i) = min (tem, xmbmax (i))
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! compute cloud moisture property and precipitation
    ! set cloud moisture property equal to the enviromental moisture at updraft starting level (kb) .
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            aa1 (i) = 0.
            qcko (i, kb (i)) = qo (i, kb (i))
            qrcko (i, kb (i)) = qo (i, kb (i))
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! calculate the moisture content of the entraining / detraining parcel (qcko) and the value it would have if just saturated (qrch), according to equation a.14 in grell (1993). their difference is the amount of convective cloud water (qlk = rain + condensate) . determine the portion of convective cloud water that remains suspended and the portion that is converted into convective precipitation (pwo) . calculate and save the negative cloud work function (aa1) due to water loading. above the level of minimum moist static energy, some of the cloud water is detrained into the grid - scale cloud water from every cloud layer with a rate of 0.0005 \f$m^{ - 1}\f$ (dellal) .
    ! -----------------------------------------------------------------------

    do k = 2, km1
        do i = 1, im
            if (cnvflg (i)) then
                if (k > kb (i) .and. k < ktcon (i)) then
                    dz = zi (i, k) - zi (i, k - 1)
                    gamma = el2orc * qeso (i, k) / (to (i, k) ** 2)
                    qrch = qeso (i, k) &
                         + gamma * dbyo (i, k) / (hlv * (1. + gamma))

                    tem = 0.5 * (xlamue (i, k) + xlamue (i, k - 1)) * dz
                    tem1 = 0.5 * xlamud (i) * dz
                    factor = 1. + tem - tem1
                    qcko (i, k) = ((1. - tem1) * qcko (i, k - 1) + tem * 0.5 * &
                         (qo (i, k) + qo (i, k - 1))) / factor
                    qrcko (i, k) = qcko (i, k)

                    dq = eta (i, k) * (qcko (i, k) - qrch)

                    ! rhbar (i) = rhbar (i) + qo (i, k) / qeso (i, k)

                    ! -----------------------------------------------------------------------
                    ! below lfc check if there is excess moisture to release latent heat
                    ! -----------------------------------------------------------------------

                    if (k >= kbcon (i) .and. dq > 0.) then
                        etah = .5 * (eta (i, k) + eta (i, k - 1))
                        dp = delp (i, k)
                        if (ncloud > 0) then
                            ptem = c0t (i, k) + c1_shal
                            qlk = dq / (eta (i, k) + etah * ptem * dz)
                            dellal (i, k) = etah * c1_shal * dz * qlk * g / dp
                        else
                            qlk = dq / (eta (i, k) + etah * c0t (i, k) * dz)
                        endif
                        buo (i, k) = buo (i, k) - g * qlk
                        qcko (i, k) = qlk + qrch
                        pwo (i, k) = etah * c0t (i, k) * dz * qlk
                        cnvwt (i, k) = etah * qlk * g / dp
                    endif

                    ! -----------------------------------------------------------------------
                    ! compute buoyancy and drag for updraft velocity
                    ! -----------------------------------------------------------------------

                    if (k >= kbcon (i)) then
                        rfact = 1. + delta * cp_air * gamma &
                             * to (i, k) / hlv
                        buo (i, k) = buo (i, k) + (g / (cp_air * to (i, k))) &
                             * dbyo (i, k) / (1. + gamma) &
                             * rfact
                        val = 0.
                        buo (i, k) = buo (i, k) + g * delta * &
                            max (val, (qeso (i, k) - qo (i, k)))
                        drag (i, k) = max (xlamue (i, k), xlamud (i))
                        ! kgao 12 / 18 / 2023
                        tem = ((uo (i, k) - uo (i, k - 1)) / dz) ** 2
                        tem = tem + ((vo (i, k) - vo (i, k - 1)) / dz) ** 2
                        wush (i, k) = csmf * sqrt (tem)
                    endif

                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! calculate cloud work function
    ! -----------------------------------------------------------------------

    ! do k = 2, km1
    ! do i = 1, im
    ! if (cnvflg (i)) then
    ! if (k >= kbcon (i) .and. k < ktcon (i)) then
    ! dz1 = zo (i, k + 1) - zo (i, k)
    ! gamma = el2orc * qeso (i, k) / (to (i, k) ** 2)
    ! rfact = 1. + delta * cp_air * gamma &
    ! * to (i, k) / hlv
    ! aa1 (i) = aa1 (i) + &
    ! dz1 * eta (i, k) * (g / (cp_air * to (i, k))) &
    ! dz1 * (g / (cp_air * to (i, k))) &
    ! * dbyo (i, k) / (1. + gamma) &
    ! * rfact
    ! val = 0.
    ! aa1 (i) = aa1 (i) + &
    ! dz1 * eta (i, k) * g * delta * &
    ! dz1 * g * delta * &
    ! max (val, (qeso (i, k) - qo (i, k)))
    ! endif
    ! endif
    ! enddo
    ! enddo
    ! do i = 1, im
    ! if (cnvflg (i) .and. aa1 (i) <= 0.) cnvflg (i) = .false.
    ! enddo

    ! -----------------------------------------------------------------------
    ! calculate cloud work function
    ! -----------------------------------------------------------------------


    ! -----------------------------------------------------------------------
    ! calculate the cloud work function according to pan and wu (1995) equation 4:
    ! \f[
    ! a_u = \int_{z_0}^{z_t}\frac{g}{c_pt (z) }\frac{\eta}{1 + \gamma}[h (z) - h^ * (z) ]dz
    ! \f]
    ! (discretized according to grell (1993) equation b.10 using b.2 and b.3 of arakawa and schubert (1974) and assuming \f$\eta = 1\f$) where \f$a_u\f$ is the updraft cloud work function, \f$z_0\f$ and \f$z_t\f$ are cloud base and cloud top, respectively, \f$\gamma = \frac{l}{c_p}\left (\frac{\partial \overline{q_s}}{\partial t}\right) _p\f$ and other quantities are previously defined.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            aa1 (i) = 0.
        endif
    enddo

    do k = 2, km1
        do i = 1, im
            if (cnvflg (i)) then
                if (k >= kbcon (i) .and. k < ktcon (i)) then
                    dz1 = zo (i, k + 1) - zo (i, k)
                    aa1 (i) = aa1 (i) + buo (i, k) * dz1
                endif
            endif
        enddo
    enddo

    do i = 1, im
        if (cnvflg (i) .and. aa1 (i) <= 0.) cnvflg (i) = .false.
    enddo

    ! -----------------------------------------------------------------------
    ! if the updraft cloud work function is negative, convection does not occur, and the scheme returns to the calling routine.
    ! -----------------------------------------------------------------------

    totflg = .true.
    do i = 1, im
        totflg = totflg .and. (.not. cnvflg (i))
    enddo

    if (totflg) return

    ! -----------------------------------------------------------------------
    ! estimate the onvective overshooting as the level
    ! where the [aafac * cloud work function] becomes zero,
    ! which is the final cloud top
    ! limited to the level of p / ps = 0.7
    ! continue calculating the cloud work function past the point of neutral buoyancy to represent overshooting according to han and pan (2011). convective overshooting stops when \f$ ca_u < 0\f$ where \f$c\f$ is currently 10%, or when 10% of the updraft cloud work function has been consumed by the stable buoyancy force. overshooting is also limited to the level where \f$p = 0.7p_{sfc}\f$.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            aa1 (i) = aafac * aa1 (i)
        endif
    enddo

    do i = 1, im
        flg (i) = cnvflg (i)
        ktcon1 (i) = kbm (i)
    enddo

    do k = 2, km1
        do i = 1, im
            if (flg (i)) then
                if (k >= ktcon (i) .and. k < kbm (i)) then
                    dz1 = zo (i, k + 1) - zo (i, k)
                    gamma = el2orc * qeso (i, k) / (to (i, k) ** 2)
                    rfact = 1. + delta * cp_air * gamma &
                         * to (i, k) / hlv
                    aa1 (i) = aa1 (i) + &
                        ! dz1 * eta (i, k) * (g / (cp_air * to (i, k))) &
                        dz1 * (g / (cp_air * to (i, k))) &
                         * dbyo (i, k) / (1. + gamma) &
                         * rfact
                    ! val = 0.
                    ! aa1 (i) = aa1 (i) + &
                    ! dz1 * eta (i, k) * g * delta * &
                    ! dz1 * g * delta * &
                    ! max (val, (qeso (i, k) - qo (i, k)))
                    if (aa1 (i) < 0.) then
                        ktcon1 (i) = k
                        flg (i) = .false.
                    endif
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute cloud moisture property, detraining cloud water
    ! and precipitation in overshooting layers
    ! for the overshooting convection, calculate the moisture content of the entraining / detraining parcel as before. partition convective cloud water and precipitation and detrain convective cloud water in the overshooting layers.
    ! -----------------------------------------------------------------------

    do k = 2, km1
        do i = 1, im
            if (cnvflg (i)) then
                if (k >= ktcon (i) .and. k < ktcon1 (i)) then
                    dz = zi (i, k) - zi (i, k - 1)
                    gamma = el2orc * qeso (i, k) / (to (i, k) ** 2)
                    qrch = qeso (i, k) &
                         + gamma * dbyo (i, k) / (hlv * (1. + gamma))

                    tem = 0.5 * (xlamue (i, k) + xlamue (i, k - 1)) * dz
                    tem1 = 0.5 * xlamud (i) * dz
                    factor = 1. + tem - tem1
                    qcko (i, k) = ((1. - tem1) * qcko (i, k - 1) + tem * 0.5 * &
                         (qo (i, k) + qo (i, k - 1))) / factor
                    qrcko (i, k) = qcko (i, k)

                    dq = eta (i, k) * (qcko (i, k) - qrch)

                    ! -----------------------------------------------------------------------
                    ! check if there is excess moisture to release latent heat
                    ! -----------------------------------------------------------------------

                    if (dq > 0.) then
                        etah = .5 * (eta (i, k) + eta (i, k - 1))
                        dp = delp (i, k)
                        if (ncloud > 0) then
                            ptem = c0t (i, k) + c1_shal
                            qlk = dq / (eta (i, k) + etah * ptem * dz)
                            dellal (i, k) = etah * c1_shal * dz * qlk * g / dp
                        else
                            qlk = dq / (eta (i, k) + etah * c0t (i, k) * dz)
                        endif
                        qcko (i, k) = qlk + qrch
                        pwo (i, k) = etah * c0t (i, k) * dz * qlk
                        cnvwt (i, k) = etah * qlk * g / dp
                    endif
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute updraft velocity square (wu2)
    ! calculate updraft velocity square (wu2) according to han et al.'s (2017) equation 7.
    ! -----------------------------------------------------------------------

    ! bb1 = 2. * (1. + bet1 * cd1)
    ! bb2 = 2. / (f1 * (1. + gam1))

    ! bb1 = 3.9
    ! bb2 = 0.67

    ! bb1 = 2.0
    ! bb2 = 4.0

    ! bb1 = 4.0
    ! bb2 = 0.8

    ! do i = 1, im
    !     if (cnvflg (i)) then
    !         k = kbcon1 (i)
    !         tem = po (i, k) / (rdgas * to (i, k))
    !         wucb = - 0.01 * dot (i, k) / (tem * g)
    !         if (wucb > 0.) then
    !             wu2 (i, k) = wucb * wucb
    !         else
    !             wu2 (i, k) = 0.
    !         endif
    !     endif
    ! enddo

    do k = 2, km1
        do i = 1, im
            if (cnvflg (i)) then
                if (k > kbcon1 (i) .and. k < ktcon (i)) then
                    dz = zi (i, k) - zi (i, k - 1)
                    tem = 0.25 * bb1 * (drag (i, k) + drag (i, k - 1)) * dz
                    tem1 = 0.5 * bb2 * (buo (i, k) + buo (i, k - 1)) * dz
                    ! kgao 12 / 18 / 2023 - considers shear effect on updraft
                    if (use_shear_conv) then
                        tem2 = wush (i, k) * sqrt (wu2 (i, k - 1))
                        tem2 = (tem1 - tem2) * dz
                        ptem = (1. - tem) * wu2 (i, k - 1)
                        ptem1 = 1. + tem
                        wu2 (i, k) = (ptem + tem2) / ptem1
                    else
                        ptem = (1. - tem) * wu2 (i, k - 1)
                        ptem1 = 1. + tem
                        wu2 (i, k) = (ptem + tem1) / ptem1
                    endif
                    wu2 (i, k) = max (wu2 (i, k), 0.)
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute updraft velocity averaged over the whole cumulus
    ! calculate the mean updraft velocity within the cloud (wc) .
    ! -----------------------------------------------------------------------

    do i = 1, im
        wc (i) = 0.
        sumx (i) = 0.
    enddo

    do k = 2, km1
        do i = 1, im
            if (cnvflg (i)) then
                if (k > kbcon1 (i) .and. k < ktcon (i)) then
                    dz = zi (i, k) - zi (i, k - 1)
                    tem = 0.5 * (sqrt (wu2 (i, k)) + sqrt (wu2 (i, k - 1)))
                    wc (i) = wc (i) + tem * dz
                    sumx (i) = sumx (i) + dz
                endif
            endif
        enddo
    enddo

    do i = 1, im
        if (cnvflg (i)) then
            if (sumx (i) == 0.) then
                cnvflg (i) = .false.
            else
                wc (i) = wc (i) / sumx (i)
            endif
            val = 1.e-4
            if (wc (i) < val) cnvflg (i) = .false.
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! exchange ktcon with ktcon1
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            kk = ktcon (i)
            ktcon (i) = ktcon1 (i)
            ktcon1 (i) = kk
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! this section is ready for cloud water
    ! -----------------------------------------------------------------------

    if (ncloud > 0) then

        ! -----------------------------------------------------------------------
        ! compute liquid and vapor separation at cloud top
        ! separate the total updraft cloud water at cloud top into vapor and condensate.
        ! -----------------------------------------------------------------------

        do i = 1, im
            if (cnvflg (i)) then
                k = ktcon (i) - 1
                gamma = el2orc * qeso (i, k) / (to (i, k) ** 2)
                qrch = qeso (i, k) &
                     + gamma * dbyo (i, k) / (hlv * (1. + gamma))
                dq = qcko (i, k) - qrch

                ! -----------------------------------------------------------------------
                ! check if there is excess moisture to release latent heat
                ! -----------------------------------------------------------------------

                if (dq > 0.) then
                    qlko_ktcon (i) = dq
                    qcko (i, k) = qrch
                endif
            endif
        enddo
    endif

    ! -----------------------------------------------------------------------
    ! compute precipitation efficiency in terms of windshear
    ! calculate the wind shear and precipitation efficiency according to equation 58 in fritsch and chappell (1980):
    ! \f[
    ! e = 1.591 - 0.639\frac{\delta v}{\delta z} + 0.0953\left (\frac{\delta v}{\delta z}\right) ^2 - 0.00496\left (\frac{\delta v}{\delta z}\right) ^3
    ! \f]
    ! where \f$\delta v\f$ is the integrated horizontal shear over the cloud depth, \f$\delta z\f$, (the ratio is converted to units of \f$10^{ - 3} s^{ - 1}\f$) . the variable "edt" is \f$1 - e\f$ and is constrained to the range \f$[0, 0.9]\f$.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            vshear (i) = 0.
        endif
    enddo

    do k = 2, km
        do i = 1, im
            if (cnvflg (i)) then
                if (k > kb (i) .and. k <= ktcon (i)) then
                    shear = sqrt ((uo (i, k) - uo (i, k - 1)) ** 2 &
                         + (vo (i, k) - vo (i, k - 1)) ** 2)
                    vshear (i) = vshear (i) + shear
                endif
            endif
        enddo
    enddo

    do i = 1, im
        if (cnvflg (i)) then
            vshear (i) = 1.e3 * vshear (i) / (zi (i, ktcon (i)) - zi (i, kb (i)))
            e1 = 1.591 - .639 * vshear (i) &
                 + .0953 * (vshear (i) ** 2) - .00496 * (vshear (i) ** 3)
            edt (i) = 1. - e1
            val = .9
            edt (i) = min (edt (i), val)
            val = .0
            edt (i) = max (edt (i), val)
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! what would the change be, that a cloud with unit mass
    ! will do to the environment?
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! calculate the tendencies of the state variables (per unit cloud base mass flux) and the cloud base mass flux.
    ! calculate the change in moist static energy, moisture mixing ratio, and horizontal winds per unit cloud base mass flux for all layers below cloud top from equations b.14 and b.15 from grell (1993), and for the cloud top from b.16 and b.17.
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (cnvflg (i) .and. k <= kmax (i)) then
                dellah (i, k) = 0.
                dellaq (i, k) = 0.
                dellau (i, k) = 0.
                dellav (i, k) = 0.
            endif
        enddo
    enddo

    do n = 1, ntr
        do k = 1, km
            do i = 1, im
                if (cnvflg (i) .and. k <= kmax (i)) then
                    dellae (i, k, n) = 0.
                endif
            enddo
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! changed due to subsidence and entrainment
    ! -----------------------------------------------------------------------

    do k = 2, km1
        do i = 1, im
            if (cnvflg (i)) then
                if (k > kb (i) .and. k < ktcon (i)) then
                    dp = delp (i, k)
                    dz = zi (i, k) - zi (i, k - 1)

                    dv1h = heo (i, k)
                    dv2h = .5 * (heo (i, k) + heo (i, k - 1))
                    dv3h = heo (i, k - 1)
                    dv1q = qo (i, k)
                    dv2q = .5 * (qo (i, k) + qo (i, k - 1))
                    dv3q = qo (i, k - 1)

                    tem = 0.5 * (xlamue (i, k) + xlamue (i, k - 1))
                    tem1 = xlamud (i)

                    dellah (i, k) = dellah (i, k) + &
                         (eta (i, k) * dv1h - eta (i, k - 1) * dv3h &
                         - tem * eta (i, k - 1) * dv2h * dz &
                         + tem1 * eta (i, k - 1) * .5 * (hcko (i, k) + hcko (i, k - 1)) * dz &
                        ) * g / dp

                    dellaq (i, k) = dellaq (i, k) + &
                         (eta (i, k) * dv1q - eta (i, k - 1) * dv3q &
                         - tem * eta (i, k - 1) * dv2q * dz &
                         + tem1 * eta (i, k - 1) * .5 * (qrcko (i, k) + qcko (i, k - 1)) * dz &
                        ) * g / dp

                    tem1 = eta (i, k) * (uo (i, k) - ucko (i, k))
                    tem2 = eta (i, k - 1) * (uo (i, k - 1) - ucko (i, k - 1))
                    dellau (i, k) = dellau (i, k) + (tem1 - tem2) * g / dp

                    tem1 = eta (i, k) * (vo (i, k) - vcko (i, k))
                    tem2 = eta (i, k - 1) * (vo (i, k - 1) - vcko (i, k - 1))
                    dellav (i, k) = dellav (i, k) + (tem1 - tem2) * g / dp

                endif
            endif
        enddo
    enddo

    do n = 1, ntr
        do k = 2, km1
            do i = 1, im
                if (cnvflg (i)) then
                    if (k > kb (i) .and. k < ktcon (i)) then
                        dp = delp (i, k)
                        tem1 = eta (i, k) * (ctro (i, k, n) - ecko (i, k, n))
                        tem2 = eta (i, k - 1) * (ctro (i, k - 1, n) - ecko (i, k - 1, n))
                        dellae (i, k, n) = dellae (i, k, n) + (tem1 - tem2) * g / dp
                    endif
                endif
            enddo
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! cloud top
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            indx = ktcon (i)
            dp = delp (i, indx)
            dv1h = heo (i, indx - 1)
            dellah (i, indx) = eta (i, indx - 1) * &
                 (hcko (i, indx - 1) - dv1h) * g / dp
            dv1q = qo (i, indx - 1)
            dellaq (i, indx) = eta (i, indx - 1) * &
                 (qcko (i, indx - 1) - dv1q) * g / dp
            dellau (i, indx) = eta (i, indx - 1) * &
                 (ucko (i, indx - 1) - uo (i, indx - 1)) * g / dp
            dellav (i, indx) = eta (i, indx - 1) * &
                 (vcko (i, indx - 1) - vo (i, indx - 1)) * g / dp

            ! -----------------------------------------------------------------------
            ! cloud water
            ! -----------------------------------------------------------------------

            dellal (i, indx) = eta (i, indx - 1) * &
                qlko_ktcon (i) * g / dp
        endif
    enddo

    do n = 1, ntr
        do i = 1, im
            if (cnvflg (i)) then
                indx = ktcon (i)
                dp = delp (i, indx)
                dellae (i, indx, n) = eta (i, indx - 1) * &
                    (ecko (i, indx - 1, n) - ctro (i, indx - 1, n)) * g / dp
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute convective turn - over time
    ! following bechtold et al. (2008), calculate the convective turnover time using the mean updraft velocity (wc) and the cloud depth. it is also proportional to the grid size (gsize) .
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            tem = zi (i, ktcon1 (i)) - zi (i, kbcon1 (i))
            dtconv (i) = tem / wc (i)
            tfac = 1. + gsize (i) / 75000.
            dtconv (i) = tfac * dtconv (i)
            dtconv (i) = max (dtconv (i), dtmin)
            dtconv (i) = max (dtconv (i), dt2)
            dtconv (i) = min (dtconv (i), dtmax)
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! calculate advective time scale (tauadv) using a mean cloud layer wind speed.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            sumx (i) = 0.
            umean (i) = 0.
        endif
    enddo

    do k = 2, km1
        do i = 1, im
            if (cnvflg (i)) then
                if (k >= kbcon1 (i) .and. k < ktcon1 (i)) then
                    dz = zi (i, k) - zi (i, k - 1)
                    tem = sqrt (u1 (i, k) * u1 (i, k) + v1 (i, k) * v1 (i, k))
                    umean (i) = umean (i) + tem * dz
                    sumx (i) = sumx (i) + dz
                endif
            endif
        enddo
    enddo

    do i = 1, im
        if (cnvflg (i)) then
            umean (i) = umean (i) / sumx (i)
            umean (i) = max (umean (i), 1.)
            tauadv (i) = gsize (i) / umean (i)
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! compute cloud base mass flux as a function of the mean
    ! updraft velcoity
    ! from han et al.'s (2017) equation 6, calculate cloud base mass flux as a function of the mean updraft velcoity.
    ! as discussed in han et al. (2017), when dtconv is larger than tauadv, the convective mixing is not fully conducted before the cumulus cloud is advected out of the grid cell. in this case, therefore, the cloud base mass flux is further reduced in proportion to the ratio of tauadv to dtconv.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            k = kbcon (i)
            rho = po (i, k) * 100. / (rdgas * to (i, k))
            tfac = tauadv (i) / dtconv (i)
            tfac = min (tfac, 1.)
            xmb (i) = tfac * betaw_shal * rho * wc (i)
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! for scale - aware parameterization, the updraft fraction (sigmagfm) is first computed as a function of the lateral entrainment rate at cloud base (see han et al.'s (2017) equation 4 and 5), following the study by grell and freitas (2014).
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            tem = min (max (xlamue (i, kbcon (i)), 2.e-4), 6.e-4)
            tem = 0.2 / tem
            tem1 = 3.14 * tem * tem
            sigmagfm (i) = tem1 / (gsize (i) ** 2.0)
            sigmagfm (i) = max (sigmagfm (i), 0.001)
            sigmagfm (i) = min (sigmagfm (i), 0.999)
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! then, calculate the reduction factor (scaldfunc) of the vertical convective eddy transport of mass flux as a function of updraft fraction from the studies by arakawa and wu (2013) (also see han et al.'s (2017) equation 1 and 2) . the final cloud base mass flux with scale - aware parameterization is obtained from the mass flux when sigmagfm < < 1, multiplied by the reduction factor (han et al.'s (2017) equation 2) .
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            if (gsize (i) < dxcrt_shal) then
                scaldfunc (i) = (1. - sigmagfm (i)) * (1. - sigmagfm (i))
                scaldfunc (i) = max (min (scaldfunc (i), 1.0), 0.)
            else
                scaldfunc (i) = 1.0
            endif
            xmb (i) = xmb (i) * scaldfunc (i)
            xmb (i) = min (xmb (i), xmbmax (i))
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! transport aerosols if present
    ! -----------------------------------------------------------------------

    if (do_aerosols) &
        call sa_aamf_shal_aero (im, km, itc, ntc, ntr, delt, &
        cnvflg, kb, kmax, kbcon, ktcon, fscav, &
        xmb, c0t, eta, zi, xlamue, xlamud, delp, &
        qtr, qaero)

    ! -----------------------------------------------------------------------
    ! for the "feedback control", calculate updated values of the state variables by multiplying the cloud base mass flux and the tendencies calculated per unit cloud base mass flux from the static control.
    ! - recalculate saturation specific humidity.
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (cnvflg (i) .and. k <= kmax (i)) then
                qeso (i, k) = 0.01 * mqs (t1 (i, k)) ! mqs is in pa
                qeso (i, k) = eps * qeso (i, k) / (pfld (i, k) + epsm1 * qeso (i, k))
                val = 1.e-8
                qeso (i, k) = max (qeso (i, k), val)
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! calculate the temperature tendency from the moist static energy and specific humidity tendencies.
    ! update the temperature, specific humidity, and horiztonal wind state variables by multiplying the cloud base mass flux - normalized tendencies by the cloud base mass flux.
    ! accumulate column - integrated tendencies.
    ! -----------------------------------------------------------------------

    do i = 1, im
        delhbar (i) = 0.
        delqbar (i) = 0.
        deltbar (i) = 0.
        delubar (i) = 0.
        delvbar (i) = 0.
        qcond (i) = 0.
    enddo

    do n = 1, ntr
        do i = 1, im
            delebar (i, n) = 0.
        enddo
    enddo

    do k = 1, km
        do i = 1, im
            if (cnvflg (i)) then
                if (k > kb (i) .and. k <= ktcon (i)) then
                    dellat = (dellah (i, k) - hlv * dellaq (i, k)) / cp_air
                    t1 (i, k) = t1 (i, k) + dellat * xmb (i) * dt2
                    q1 (i, k) = q1 (i, k) + dellaq (i, k) * xmb (i) * dt2
                    ! tem = 1. / rcs (i)
                    ! u1 (i, k) = u1 (i, k) + dellau (i, k) * xmb (i) * dt2 * tem
                    ! v1 (i, k) = v1 (i, k) + dellav (i, k) * xmb (i) * dt2 * tem
                    u1 (i, k) = u1 (i, k) + dellau (i, k) * xmb (i) * dt2
                    v1 (i, k) = v1 (i, k) + dellav (i, k) * xmb (i) * dt2
                    dp = delp (i, k)
                    delhbar (i) = delhbar (i) + dellah (i, k) * xmb (i) * dp / g
                    delqbar (i) = delqbar (i) + dellaq (i, k) * xmb (i) * dp / g
                    deltbar (i) = deltbar (i) + dellat * xmb (i) * dp / g
                    delubar (i) = delubar (i) + dellau (i, k) * xmb (i) * dp / g
                    delvbar (i) = delvbar (i) + dellav (i, k) * xmb (i) * dp / g
                endif
            endif
        enddo
    enddo

    kk = 0
    do n = 1, ntr + 2
        if (n .eq. ntw .or. n .eq. nti) cycle
        kk = kk + 1
        do k = 1, km
            do i = 1, im
                if (cnvflg (i) .and. k <= kmax (i)) then
                    if (k <= ktcon (i)) then
                        ctr (i, k, kk) = ctr (i, k, kk) + dellae (i, k, kk) * xmb (i) * dt2
                        delebar (i, kk) = delebar (i, kk) + dellae (i, k, kk) * xmb (i) * dp / g
                        qtr (i, k, n) = ctr (i, k, kk)
                    endif
                endif
            enddo
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! recalculate saturation specific humidity using the updated temperature.
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (cnvflg (i)) then
                if (k > kb (i) .and. k <= ktcon (i)) then
                    qeso (i, k) = 0.01 * mqs (t1 (i, k)) ! mqs is in pa
                    qeso (i, k) = eps * qeso (i, k) / (pfld (i, k) + epsm1 * qeso (i, k))
                    val = 1.e-8
                    qeso (i, k) = max (qeso (i, k), val)
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! add up column - integrated convective precipitation by multiplying the normalized value by the cloud base mass flux.
    ! -----------------------------------------------------------------------

    do i = 1, im
        rntot (i) = 0.
        delqev (i) = 0.
        delq2 (i) = 0.
        flg (i) = cnvflg (i)
    enddo

    do k = km, 1, - 1
        do i = 1, im
            if (cnvflg (i)) then
                if (k < ktcon (i) .and. k > kb (i)) then
                    rntot (i) = rntot (i) + pwo (i, k) * xmb (i) * .001 * dt2
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! evaporating rain
    ! determine the evaporation of the convective precipitation and update the integrated convective precipitation.
    ! update state temperature and moisture to account for evaporation of convective precipitation.
    ! update column - integrated tendencies to account for evaporation of convective precipitation.
    ! -----------------------------------------------------------------------

    do k = km, 1, - 1
        do i = 1, im
            if (k <= kmax (i)) then
                deltv (i) = 0.
                delq (i) = 0.
                qevap (i) = 0.
                if (cnvflg (i)) then
                    if (k < ktcon (i) .and. k > kb (i)) then
                        rn (i) = rn (i) + pwo (i, k) * xmb (i) * .001 * dt2
                        qr (i, k) = qr (i, k) + pwo (i, k) * xmb (i) * .001 * dt2
                    endif
                endif
                if (flg (i) .and. k < ktcon (i)) then
                    evef = edt (i) * evfact_shal
                    if (islimsk (i) == 1) evef = edt (i) * evfactl_shal
                    ! if (islimsk (i) == 1) evef = .07
                    ! if (islimsk (i) == 1) evef = 0.
                    qcond (i) = evef * (q1 (i, k) - qeso (i, k)) &
                         / (1. + el2orc * qeso (i, k) / t1 (i, k) ** 2)
                    dp = delp (i, k)
                    if (rn (i) > 0. .and. qcond (i) < 0.) then
                        qevap (i) = - qcond (i) * (1. - exp (- .32 * sqrt (dt2 * rn (i))))
                        qevap (i) = min (qevap (i), rn (i) * 1000. * g / dp)
                        delq2 (i) = delqev (i) + .001 * qevap (i) * dp / g
                    endif
                    if (rn (i) > 0. .and. qcond (i) < 0. .and. delq2 (i) > rntot (i)) then
                        qevap (i) = 1000. * g * (rntot (i) - delqev (i)) / dp
                        flg (i) = .false.
                    endif
                    if (rn (i) > 0. .and. qevap (i) > 0.) then
                        tem = .001 * dp / g
                        tem1 = qevap (i) * tem
                        if (tem1 > rn (i)) then
                            qevap (i) = rn (i) / tem
                            rn (i) = 0.
                        else
                            rn (i) = rn (i) - tem1
                        endif
                        qr (i, k) = qr (i, k) - qevap (i) * tem
                        q1 (i, k) = q1 (i, k) + qevap (i)
                        t1 (i, k) = t1 (i, k) - elocp * qevap (i)
                        deltv (i) = - elocp * qevap (i) / dt2
                        delq (i) = + qevap (i) / dt2
                        delqev (i) = delqev (i) + .001 * dp * qevap (i) / g
                    endif
                    delqbar (i) = delqbar (i) + delq (i) * dp / g
                    deltbar (i) = deltbar (i) + deltv (i) * dp / g
                endif
            endif
        enddo
    enddo

    ! do i = 1, im
    ! if (me == 31 .and. cnvflg (i)) then
    ! if (cnvflg (i)) then
    ! print *, ' shallow delhbar, delqbar, deltbar = ', &
    ! delhbar (i), hlv * delqbar (i), cp_air * deltbar (i)
    ! print *, ' shallow delubar, delvbar = ', delubar (i), delvbar (i)
    ! print *, ' precip = ', hlv * rn (i) * 1000. / dt2
    ! print *, 'pdif = ', pfld (i, kbcon (i)) - pfld (i, ktcon (i))
    ! endif
    ! enddo
    ! do n = 1, ntr
    ! do i = 1, im
    ! if (me == 31 .and. cnvflg (i)) then
    ! if (cnvflg (i)) then
    ! print *, ' tracer delebar = ', delebar (i, n)
    ! endif
    ! enddo
    ! enddo

    do i = 1, im
        if (cnvflg (i)) then
            if (rn (i) < 0. .or. .not.flg (i)) rn (i) = 0.
            ktop (i) = ktcon (i)
            kbot (i) = kbcon (i)
            kcnv (i) = 2
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! convective cloud water
    ! calculate shallow convective cloud water.
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (present (cnvw) .and. cnvflg (i)) then
                if (k >= kbcon (i) .and. k < ktcon (i)) then
                    cnvw (i, k) = cnvwt (i, k) * xmb (i) * dt2
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! convective cloud cover
    ! calculate convective cloud cover, which is used when pdf - based cloud fraction is used (i.e., pdfcld = .true.) .
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (present (cnvc) .and. cnvflg (i)) then
                if (k >= kbcon (i) .and. k < ktcon (i)) then
                    cnvc (i, k) = 0.04 * log (1. + 675. * eta (i, k) * xmb (i))
                    cnvc (i, k) = min (cnvc (i, k), 0.2)
                    cnvc (i, k) = max (cnvc (i, k), 0.0)
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! cloud water
    ! separate detrained cloud water into liquid and ice species as a function of temperature only.
    ! -----------------------------------------------------------------------

    if (ncloud > 0) then

        do k = 1, km1
            do i = 1, im
                if (cnvflg (i)) then
                    ! if (k > kb (i) .and. k <= ktcon (i)) then
                    if (k >= kbcon (i) .and. k <= ktcon (i)) then
                        tem = dellal (i, k) * xmb (i) * dt2
                        qtr (i, k, ntw) = qtr (i, k, ntw) + tem
                    endif
                endif
            enddo
        enddo

    endif

    ! -----------------------------------------------------------------------
    ! store aerosol concentrations if present
    ! -----------------------------------------------------------------------

    if (do_aerosols) then
        do n = 1, ntc
            kk = n + itc - 1
            do k = 1, km
                do i = 1, im
                    if (cnvflg (i) .and. rn (i) > 0.) then
                        if (k <= kmax (i)) qtr (i, k, kk) = qaero (i, k, n)
                    endif
                enddo
            enddo
        enddo
    endif

    ! -----------------------------------------------------------------------
    ! hchuang code change
    ! calculate and retain the updraft mass flux for dust transport by cumulus convection.
    ! calculate the updraft convective mass flux.
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (present (ud_mf) .and. cnvflg (i)) then
                if (k >= kb (i) .and. k < ktop (i)) then
                    ud_mf (i, k) = eta (i, k) * xmb (i) * dt2
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! save the updraft convective mass flux at cloud top.
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (present (dt_mf) .and. present (ud_mf) .and. cnvflg (i)) then
            k = ktop (i) - 1
            dt_mf (i, k) = ud_mf (i, k)
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! include tke contribution from shallow convection
    ! -----------------------------------------------------------------------

    if (ntk > 0) then

        do k = 2, km1
            do i = 1, im
                if (cnvflg (i)) then
                    if (k > kb (i) .and. k < ktop (i)) then
                        tem = 0.5 * (eta (i, k - 1) + eta (i, k)) * xmb (i)
                        tem1 = pfld (i, k) * 100. / (rdgas * t1 (i, k))
                        sigmagfm (i) = max (sigmagfm (i), betaw_shal)
                        ptem = tem / (sigmagfm (i) * tem1)
                        qtr (i, k, ntk) = qtr (i, k, ntk) + 0.5 * sigmagfm (i) * ptem * ptem
                    endif
                endif
            enddo
        enddo

    endif

end subroutine sa_aamf_shal

! =======================================================================
! Aerosol Transportation in Deep Convection
! =======================================================================

subroutine sa_aamf_deep_aero (im, km, itc, ntc, ntr, delt, &
        xlamde, xlamdd, cnvflg, jmin, kb, kmax, kbcon, ktcon, fscav, &
        edto, xlamd, xmb, c0t, eta, etad, zi, xlamue, xlamud, delp, &
        qtr, qaero)

    implicit none

    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------

    integer, intent (in) :: im, km, itc, ntc, ntr

    real, intent (in) :: delt, xlamde, xlamdd

    logical, dimension (im), intent (in) :: cnvflg

    integer, dimension (im), intent (in) :: jmin, kb, kmax, kbcon, ktcon

    real, dimension (im), intent (in) :: edto,xlamd, xmb

    real, dimension (ntc), intent (in) :: fscav

    real, dimension (im, km), intent (in) :: c0t, eta, etad, zi, xlamue, xlamud

    real, dimension (im, km), intent (in) :: delp

    real, dimension (im, km, ntr + 2), intent (in) :: qtr

    real, dimension (im, km, ntc), intent (out) :: qaero

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    ! general variables

    integer :: i, indx, it, k, kk, km1, kp1, n

    real :: adw, aup, dtime_max, dv1q, dv2q, dv3q, dtovdz, dz, factor, ptem, ptem1, qamax, tem, tem1

    real, dimension (im, km) :: xmbp

    ! chemical transport variables

    real, dimension (im, km, ntc) :: ctro2, ecko2, ecdo2, dellae2

    ! additional variables for tracers for wet deposition,

    real, dimension (im, km, ntc) :: chem_c, chem_pw, wet_dep

    ! if reevaporation is enabled, uncomment lines below

    ! real, dimension (im, ntc) :: pwav
    ! real, dimension (im, km) :: pwdper
    ! real, dimension (im, km, ntr) :: chem_pwd

    !  additional variables for fct

    real, dimension (im, km) :: flx_lo, totlout, clipout

    real, parameter :: one = 1.0
    real, parameter :: half = 0.5
    real, parameter :: quarter = 0.25
    real, parameter :: zero = 0.0
    real, parameter :: epsil = 1.e-22 ! prevent division by zero

    ! -----------------------------------------------------------------------
    ! begin
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! check if aerosols are present
    ! -----------------------------------------------------------------------

    if (ntc <= 0 .or. itc <= 0 .or. ntr <= 0) return
    if (ntr < itc + ntc - 3) return

    ! -----------------------------------------------------------------------
    ! initialize work variables
    ! -----------------------------------------------------------------------

    km1 = km - 1

    chem_c = zero
    chem_pw = zero
    ctro2 = zero
    dellae2 = zero
    ecdo2 = zero
    ecko2 = zero
    qaero = zero

    ! -----------------------------------------------------------------------
    ! set work arrays
    ! -----------------------------------------------------------------------

    do n = 1, ntc
        it = n + itc - 1
        do k = 1, km
            do i = 1, im
                if (k <= kmax (i)) qaero (i, k, n) = max (qcmin, qtr (i, k, it))
            enddo
        enddo
    enddo

    do k = 1, km
        do i = 1, im
            xmbp (i, k) = grav * xmb (i) / delp (i, k)
        enddo
    enddo

    do n = 1, ntc

        ! -----------------------------------------------------------------------
        ! interface level
        ! -----------------------------------------------------------------------

        do k = 1, km1
            kp1 = k + 1
            do i = 1, im
                if (kp1 <= kmax (i)) ctro2 (i, k, n) = &
                    half * (qaero (i, k, n) + qaero (i, kp1, n))
            enddo
        enddo

        ! -----------------------------------------------------------------------
        ! top level
        ! -----------------------------------------------------------------------

        do i = 1, im
            ctro2 (i, kmax (i), n) = qaero (i, kmax (i), n)
        enddo
    enddo

    do n = 1, ntc
        do k = 1, km
            do i = 1, im
                if (cnvflg (i) .and. (k <= kb (i))) &
                    ecko2 (i, k, n) = ctro2 (i, k, n)
            enddo
        enddo
    enddo

    do n = 1, ntc
        do i = 1, im
            if (cnvflg (i)) ecdo2 (i, jmin (i), n) = ctro2 (i, jmin (i), n)
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! do chemical tracers, first need to know how much reevaporates
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! aerosol re - evaporation is set to zero for now
    ! uncomment and edit the following code to enable re - evaporation
    ! chem_pwd = zero
    ! pwdper = zero
    ! pwav = zero
    ! do i = 1, im
    ! do k = 1, jmin (i)
    ! pwdper (i, k) = - edto (i) * pwdo (i, k) / pwavo (i)
    ! enddo
    ! enddo
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! calculate include mixing ratio (ecko2), how much goes into
    ! rainwater to be rained out (chem_pw), and total scavenged,
    ! if not reevaporated (pwav)
    ! -----------------------------------------------------------------------

    do n = 1, ntc
        do k = 2, km1
            kk = k - 1
            do i = 1, im
                if (cnvflg (i)) then
                    if ((k > kb (i)) .and. (k < ktcon (i))) then
                        dz = zi (i, k) - zi (i, kk)
                        tem = half * (xlamue (i, k) + xlamue (i, kk)) * dz
                        tem1 = quarter * (xlamud (i, k) + xlamud (i, kk)) * dz
                        factor = one + tem - tem1

                        ! -----------------------------------------------------------------------
                        ! if conserved (not scavenging) then
                        ! -----------------------------------------------------------------------

                        ecko2 (i, k, n) = ((one - tem1) * ecko2 (i, kk, n) &
                            + half * tem * (ctro2 (i, k, n) + ctro2 (i, kk, n))) / factor

                        ! -----------------------------------------------------------------------
                        ! how much will be scavenged
                        ! -----------------------------------------------------------------------

                        ! -----------------------------------------------------------------------
                        ! this choice was used in gf, and is also described in a
                        ! successful implementation into cesm in grl (yu et al. 2019),
                        ! it uses dimesnsionless scavenging coefficients (fscav),
                        ! but includes henry coeffs with gas phase chemistry
                        ! -----------------------------------------------------------------------

                        ! -----------------------------------------------------------------------
                        ! fraction fscav is going into liquid
                        ! -----------------------------------------------------------------------

                        chem_c (i, k, n) = fscav (n) * ecko2 (i, k, n)

                        ! -----------------------------------------------------------------------
                        ! of that part is going into rain out (chem_pw)
                        ! -----------------------------------------------------------------------

                        tem = chem_c (i, k, n) / (one + c0t (i, k) * dz)
                        chem_pw (i, k, n) = c0t (i, k) * dz * tem * eta (i, kk) !etah
                        ecko2 (i, k, n) = tem + ecko2 (i, k, n) - chem_c (i, k, n)

                        ! -----------------------------------------------------------------------
                        ! pwav needed fo reevaporation in downdraft
                        ! if including reevaporation, please uncomment code below
                        ! pwav (i, n) = pwav (i, n) + chem_pw (i, k, n)
                        ! -----------------------------------------------------------------------

                    endif
                endif
            enddo
        enddo
        do k = 1, km1
            do i = 1, im
                if (k >= ktcon (i)) ecko2 (i, k, n) = ctro2 (i, k, n)
            enddo
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! reevaporation of some, pw and pwd terms needed later for dellae2
    ! -----------------------------------------------------------------------

    do n = 1, ntc
        do k = km1, 1, - 1
            kp1 = k + 1
            do i = 1, im
                if (cnvflg (i) .and. (k < jmin (i))) then
                    dz = zi (i, kp1) - zi (i, k)
                    if (k >= kbcon (i)) then
                        tem = xlamde * dz
                        tem1 = half * xlamdd * dz
                    else
                        tem = xlamde * dz
                        tem1 = half * (xlamd (i) + xlamdd) * dz
                    endif
                    factor = one + tem - tem1
                    ecdo2 (i, k, n) = ((one - tem1) * ecdo2 (i, kp1, n) &
                        + half * tem * (ctro2 (i, k, n) + ctro2 (i, kp1, n))) / factor

                    ! -----------------------------------------------------------------------
                    ! if including reevaporation, please uncomment code below
                    ! -----------------------------------------------------------------------

                    ! ecdo2 (i, k, n) = ecdo2 (i, k, n) + pwdper (i, kp1) * pwav (i, n)
                    ! chem_pwd (i, k, n) = max (zero, pwdper (i, kp1) * pwav (i, n))
                endif
            enddo
        enddo
    enddo

    do n = 1, ntc
        do i = 1, im
            if (cnvflg (i)) then

                ! -----------------------------------------------------------------------
                ! subsidence term treated in fct routine
                ! -----------------------------------------------------------------------

                dellae2 (i, 1, n) = edto (i) * etad (i, 1) * ecdo2 (i, 1, n) * xmbp (i, 1)
            endif
        enddo
    enddo

    do n = 1, ntc
        do i = 1, im
            if (cnvflg (i)) then
                k = ktcon (i)
                kk = k - 1

                ! -----------------------------------------------------------------------
                ! for the subsidence term already is considered
                ! -----------------------------------------------------------------------

                dellae2 (i, k, n) = eta (i, kk) * ecko2 (i, kk, n) * xmbp (i, k)
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! for updraft & downdraft vertical transport
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! initialize maximum allowed timestep for upstream difference approach
    ! -----------------------------------------------------------------------

    dtime_max = delt
    do k = 2, km1
        kk = k - 1
        do i = 1, im
            if (kk < ktcon (i)) dtime_max = min (dtime_max, half * delp (i, kk))
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! now for every chemistry tracer
    ! -----------------------------------------------------------------------

    do n = 1, ntc
        do k = 2, km1
            kk = k - 1
            do i = 1, im
                if (cnvflg (i) .and. (k < ktcon (i))) then
                    dz = zi (i, k) - zi (i, kk)
                    aup = one
                    if (k <= kb (i)) aup = zero
                    adw = one
                    if (k > jmin (i)) adw = zero

                    dv1q = half * (ecko2 (i, k, n) + ecko2 (i, kk, n))
                    dv2q = half * (ctro2 (i, k, n) + ctro2 (i, kk, n))
                    dv3q = half * (ecdo2 (i, k, n) + ecdo2 (i, kk, n))

                    tem = half * (xlamue (i, k) + xlamue (i, kk))
                    tem1 = half * (xlamud (i, k) + xlamud (i, kk))

                    if (k <= kbcon (i)) then
                        ptem = xlamde
                        ptem1 = xlamd (i) + xlamdd
                    else
                        ptem = xlamde
                        ptem1 = xlamdd
                    endif
                    dellae2 (i, k, n) = dellae2 (i, k, n) + &

                    ! -----------------------------------------------------------------------
                    ! detrainment from updraft
                    ! -----------------------------------------------------------------------

                        (aup * tem1 * eta (i, kk) * dv1q &

                    ! -----------------------------------------------------------------------
                    ! entrainement into up and downdraft
                    ! -----------------------------------------------------------------------

                        - (aup * tem * eta (i, kk) + adw * edto (i) * ptem * etad (i, k)) * dv2q &

                    ! -----------------------------------------------------------------------
                    ! detrainment from downdraft
                    ! -----------------------------------------------------------------------

                        + (adw * edto (i) * ptem1 * etad (i, k) * dv3q)) * dz * xmbp (i, k)

                    wet_dep (i, k, n) = chem_pw (i, k, n) * grav / delp (i, k)

                    ! -----------------------------------------------------------------------
                    ! sinks from where updraft and downdraft start
                    ! -----------------------------------------------------------------------

                    if (k == jmin (i) + 1) then
                        dellae2 (i, k, n) = dellae2 (i, k, n) &
                            - edto (i) * etad (i, kk) * ctro2 (i, kk, n) * xmbp (i, k)
                    endif
                    if (k == kb (i)) then
                        dellae2 (i, k, n) = dellae2 (i, k, n) &
                            - eta (i, k) * ctro2 (i, k, n) * xmbp (i, k)
                    endif
                endif
            enddo
        enddo

        do i = 1, im
            if (cnvflg (i)) then
                if (kb (i) == 1) then
                    k = kb (i)
                    dellae2 (i, k, n) = dellae2 (i, k, n) &
                        - eta (i, k) * ctro2 (i, k, n) * xmbp (i, k)
                endif
            endif
        enddo

    enddo

    ! -----------------------------------------------------------------------
    ! for every tracer...
    ! -----------------------------------------------------------------------

    do n = 1, ntc
        flx_lo = zero
        totlout = zero
        clipout = zero

        ! -----------------------------------------------------------------------
        ! compute low - order mass flux, upstream
        ! -----------------------------------------------------------------------

        do k = 2, km1
            kk = k - 1
            do i = 1, im
                if (cnvflg (i) .and. (kk < ktcon (i))) then
                    tem = zero
                    if (kk >= kb (i)) tem = eta (i, kk)
                    if (kk <= jmin (i)) tem = tem - edto (i) * etad (i, kk)

                    ! -----------------------------------------------------------------------
                    ! low - order flux, upstream
                    ! -----------------------------------------------------------------------

                    if (tem > zero) then
                        flx_lo (i, k) = - xmb (i) * tem * qaero (i, k, n)
                    elseif (tem < zero) then
                        flx_lo (i, k) = - xmb (i) * tem * qaero (i, kk, n)
                    endif
                endif
            enddo
        enddo

        ! -----------------------------------------------------------------------
        ! make sure low - ord fluxes don't violate positive - definiteness
        ! -----------------------------------------------------------------------

        do k = 1, km1
            kp1 = k + 1
            do i = 1, im
                if (cnvflg (i) .and. (k <= ktcon (i))) then

                    ! -----------------------------------------------------------------------
                    ! time step / grid spacing
                    ! -----------------------------------------------------------------------

                    dtovdz = grav * dtime_max / abs (delp (i, k))

                    ! -----------------------------------------------------------------------
                    ! total flux out
                    ! -----------------------------------------------------------------------

                    totlout (i, k) = max (zero, flx_lo (i, kp1)) - min (zero, flx_lo (i, k))
                    clipout (i, k) = min (one, qaero (i, k, n) / max (epsil, totlout (i, k)) &
                        / (1.0001 * dtovdz))
                endif
            enddo
        enddo

        ! -----------------------------------------------------------------------
        ! recompute upstream mass fluxes
        ! -----------------------------------------------------------------------

        do k = 2, km1
            kk = k - 1
            do i = 1, im
                if (cnvflg (i) .and. (kk < ktcon (i))) then
                    tem = zero
                    if (kk >= kb (i)) tem = eta (i, kk)
                    if (kk <= jmin (i)) tem = tem - edto (i) * etad (i, kk)
                    if (tem > zero) then
                        flx_lo (i, k) = flx_lo (i, k) * clipout (i, k)
                    elseif (tem < zero) then
                        flx_lo (i, k) = flx_lo (i, k) * clipout (i, kk)
                    endif
                endif
            enddo
        enddo

        ! -----------------------------------------------------------------------
        ! a positive - definite low - order (diffusive) solution for the subsidnce fluxes
        ! -----------------------------------------------------------------------

        do k = 1, km1
            kp1 = k + 1
            do i = 1, im
                if (cnvflg (i) .and. (k <= ktcon (i))) then
                    dtovdz = grav * dtime_max / abs (delp (i, k)) ! time step / grid spacing
                    dellae2 (i, k, n) = dellae2 (i, k, n) &
                        - (flx_lo (i, kp1) - flx_lo (i, k)) * dtovdz / dtime_max
                endif
            enddo
        enddo

    enddo ! ctr

    ! -----------------------------------------------------------------------
    ! convert wet deposition to total mass deposited over dt and dp
    ! -----------------------------------------------------------------------

    do n = 1, ntc
        do k = 1, km
            do i = 1, im
                if (cnvflg (i) .and. (k < ktcon (i))) &
                    wet_dep (i, k, n) = wet_dep (i, k, n) * xmb (i) * delt * delp (i, k)
            enddo
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute final aerosol concentrations
    ! -----------------------------------------------------------------------

    do n = 1, ntc
        do k = 1, km
            do i = 1, im
                if (cnvflg (i) .and. (k <= min (kmax (i), ktcon (i)))) then
                    qaero (i, k, n) = qaero (i, k, n) + dellae2 (i, k, n) * delt
                    if (qaero (i, k, n) < zero) then

                        ! -----------------------------------------------------------------------
                        ! add negative mass to wet deposition
                        ! -----------------------------------------------------------------------

                        wet_dep (i, k, n) = wet_dep (i, k, n) - qaero (i, k, n) * delp (i, k)
                        qaero (i, k, n) = qcmin
                    endif
                endif
            enddo
        enddo
    enddo

end subroutine sa_aamf_deep_aero

! =======================================================================
! Aerosol Transportation in Shallow Convection
! =======================================================================

subroutine sa_aamf_shal_aero (im, km, itc, ntc, ntr, delt, &
        cnvflg, kb, kmax, kbcon, ktcon, fscav, &
        xmb, c0t, eta, zi, xlamue, xlamud, delp, &
        qtr, qaero)

    implicit none

    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------

    integer, intent (in) :: im, km, itc, ntc, ntr

    real, intent (in) :: delt

    ! real, intent (in) :: xlamde, xlamdd

    logical, dimension (im), intent (in) :: cnvflg

    ! integer, dimension (im), intent (in) :: jmin

    integer, dimension (im), intent (in) :: kb, kmax, kbcon, ktcon

    real, dimension (im), intent (in) :: xmb, xlamud

    real, dimension (ntc), intent (in) :: fscav

    real, dimension (im, km), intent (in) :: c0t, eta, zi, xlamue

    real, dimension (im, km), intent (in) :: delp

    real, dimension (im, km, ntr + 2), intent (in) :: qtr

    real, dimension (im, km, ntc), intent (out) :: qaero

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    ! general variables

    integer :: i, indx, it, k, kk, km1, kp1, n

    ! real :: adw, aup, dtime_max, dv1q, dv2q, dv3q,

    real :: aup, dtime_max, dv1q, dv2q, dv3q, dtovdz, dz, factor, ptem, ptem1, qamax, tem, tem1

    real, dimension (im, km) :: xmbp

    ! chemical transport variables

    real, dimension (im, km, ntc) :: ctro2, ecko2, dellae2

    ! additional variables for tracers for wet deposition,

    real, dimension (im, km, ntc) :: chem_c, chem_pw, wet_dep

    ! if reevaporation is enabled, uncomment lines below

    ! real, dimension (im, ntc) :: pwav
    ! real, dimension (im, km) :: pwdper
    ! real, dimension (im, km, ntr) :: chem_pwd

    ! additional variables for fct

    real, dimension (im, km) :: flx_lo, totlout, clipout

    real, parameter :: one = 1.0
    real, parameter :: half = 0.5
    real, parameter :: quarter = 0.25
    real, parameter :: zero = 0.0
    real, parameter :: epsil = 1.e-22 ! prevent division by zero
    real, parameter :: escav = 0.8 ! wet scavenging efficiency

    ! -----------------------------------------------------------------------
    ! begin
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! check if aerosols are present
    ! -----------------------------------------------------------------------

    if (ntc <= 0 .or. itc <= 0 .or. ntr <= 0) return
    if (ntr < itc + ntc - 3) return

    ! -----------------------------------------------------------------------
    ! initialize work variables
    ! -----------------------------------------------------------------------

    km1 = km - 1

    chem_c = zero
    chem_pw = zero
    ctro2 = zero
    dellae2 = zero
    !ecdo2 = zero
    ecko2 = zero
    qaero = zero

    ! -----------------------------------------------------------------------
    ! set work arrays
    ! -----------------------------------------------------------------------

    do n = 1, ntc
        it = n + itc - 1
        do k = 1, km
            do i = 1, im
                if (k <= kmax (i)) qaero (i, k, n) = max (qcmin, qtr (i, k, it))
            enddo
        enddo
    enddo

    do k = 1, km
        do i = 1, im
            xmbp (i, k) = grav * xmb (i) / delp (i, k)
        enddo
    enddo

    do n = 1, ntc

        ! -----------------------------------------------------------------------
        ! interface level
        ! -----------------------------------------------------------------------

        do k = 1, km1
            kp1 = k + 1
            do i = 1, im
                if (kp1 <= kmax (i)) ctro2 (i, k, n) = &
                    half * (qaero (i, k, n) + qaero (i, kp1, n))
            enddo
        enddo

        ! -----------------------------------------------------------------------
        ! top level
        ! -----------------------------------------------------------------------

        do i = 1, im
            ctro2 (i, kmax (i), n) = qaero (i, kmax (i), n)
        enddo
    enddo

    do n = 1, ntc
        do k = 1, km
            do i = 1, im
                if (cnvflg (i) .and. (k <= kb (i))) &
                    ecko2 (i, k, n) = ctro2 (i, k, n)
            enddo
        enddo
    enddo

    !do n = 1, ntc
    ! do i = 1, im
    ! if (cnvflg (i)) ecdo2 (i, jmin (i), n) = ctro2 (i, jmin (i), n)
    ! enddo
    !enddo

    ! -----------------------------------------------------------------------
    ! do chemical tracers, first need to know how much reevaporates
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! aerosol re - evaporation is set to zero for now
    ! uncomment and edit the following code to enable re - evaporation
    ! chem_pwd = zero
    ! pwdper = zero
    ! pwav = zero
    ! do i = 1, im
    ! do k = 1, jmin (i)
    ! pwdper (i, k) = - edto (i) * pwdo (i, k) / pwavo (i)
    ! enddo
    ! enddo
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! calculate include mixing ratio (ecko2), how much goes into
    ! rainwater to be rained out (chem_pw), and total scavenged,
    ! if not reevaporated (pwav)
    ! -----------------------------------------------------------------------

    do n = 1, ntc
        do k = 2, km1
            kk = k - 1
            do i = 1, im
                if (cnvflg (i)) then
                    if ((k > kb (i)) .and. (k < ktcon (i))) then
                        dz = zi (i, k) - zi (i, kk)
                        tem = half * (xlamue (i, k) + xlamue (i, kk)) * dz
                        ! tem1 = quarter * (xlamud (i, k) + xlamud (i, kk)) * dz
                        tem1 = quarter * (xlamud (i) + xlamud (i)) * dz
                        factor = one + tem - tem1

                        ! -----------------------------------------------------------------------
                        ! if conserved (not scavenging) then
                        ! -----------------------------------------------------------------------

                        ecko2 (i, k, n) = ((one - tem1) * ecko2 (i, kk, n) &
                            + half * tem * (ctro2 (i, k, n) + ctro2 (i, kk, n))) / factor

                        ! -----------------------------------------------------------------------
                        ! how much will be scavenged
                        ! -----------------------------------------------------------------------

                        ! -----------------------------------------------------------------------
                        ! this choice was used in gf, and is also described in a
                        ! successful implementation into cesm in grl (yu et al. 2019),
                        ! it uses dimesnsionless scavenging coefficients (fscav),
                        ! but includes henry coeffs with gas phase chemistry
                        ! -----------------------------------------------------------------------

                        ! -----------------------------------------------------------------------
                        ! fraction fscav is going into liquid
                        ! -----------------------------------------------------------------------

                        chem_c (i, k, n) = escav * fscav (n) * ecko2 (i, k, n)

                        ! -----------------------------------------------------------------------
                        ! of that part is going into rain out (chem_pw)
                        ! -----------------------------------------------------------------------

                        tem = chem_c (i, k, n) / (one + c0t (i, k) * dz)
                        chem_pw (i, k, n) = c0t (i, k) * dz * tem * eta (i, kk) !etah
                        ecko2 (i, k, n) = tem + ecko2 (i, k, n) - chem_c (i, k, n)

                        ! -----------------------------------------------------------------------
                        ! pwav needed fo reevaporation in downdraft
                        ! if including reevaporation, please uncomment code below
                        ! pwav (i, n) = pwav (i, n) + chem_pw (i, k, n)
                        ! -----------------------------------------------------------------------

                    endif
                endif
            enddo
        enddo
        do k = 1, km1
            do i = 1, im
                if (k >= ktcon (i)) ecko2 (i, k, n) = ctro2 (i, k, n)
            enddo
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! reevaporation of some, pw and pwd terms needed later for dellae2
    ! -----------------------------------------------------------------------

    ! do n = 1, ntc
    ! do k = km1, 1, - 1
    ! kp1 = k + 1
    ! do i = 1, im
    ! if (cnvflg (i) .and. (k < jmin (i))) then
    ! dz = zi (i, kp1) - zi (i, k)
    ! if (k >= kbcon (i)) then
    ! tem = xlamde * dz
    ! tem1 = half * xlamdd * dz
    ! else
    ! tem = xlamde * dz
    ! tem1 = half * (xlamd (i) + xlamdd) * dz
    ! endif
    ! factor = one + tem - tem1
    ! ecdo2 (i, k, n) = ((one - tem1) * ecdo2 (i, kp1, n) &
    ! + half * tem * (ctro2 (i, k, n) + ctro2 (i, kp1, n))) / factor
    ! if including reevaporation, please uncomment code below
    ! ecdo2 (i, k, n) = ecdo2 (i, k, n) + pwdper (i, kp1) * pwav (i, n)
    ! chem_pwd (i, k, n) = max (zero, pwdper (i, kp1) * pwav (i, n))
    ! endif
    ! enddo
    ! enddo
    ! enddo

    ! do n = 1, ntc
    ! do i = 1, im
    ! if (cnvflg (i)) then
    ! subsidence term treated in fct routine
    ! dellae2 (i, 1, n) = edto (i) * etad (i, 1) * ecdo2 (i, 1, n) * xmbp (i, 1)
    ! endif
    ! enddo
    ! enddo

    do n = 1, ntc
        do i = 1, im
            if (cnvflg (i)) then
                k = ktcon (i)
                kk = k - 1

                ! -----------------------------------------------------------------------
                ! for the subsidence term already is considered
                ! -----------------------------------------------------------------------

                dellae2 (i, k, n) = eta (i, kk) * ecko2 (i, kk, n) * xmbp (i, k)
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! for updraft & downdraft vertical transport
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! initialize maximum allowed timestep for upstream difference approach
    ! -----------------------------------------------------------------------

    dtime_max = delt
    do k = 2, km1
        kk = k - 1
        do i = 1, im
            if (kk < ktcon (i)) dtime_max = min (dtime_max, half * delp (i, kk))
        enddo
    enddo

    ! now for every chemistry tracer
    do n = 1, ntc
        do k = 2, km1
            kk = k - 1
            do i = 1, im
                if (cnvflg (i) .and. (k < ktcon (i))) then
                    dz = zi (i, k) - zi (i, kk)
                    aup = one
                    if (k <= kb (i)) aup = zero
                    ! adw = one
                    ! if (k > jmin (i)) adw = zero

                    dv1q = half * (ecko2 (i, k, n) + ecko2 (i, kk, n))
                    dv2q = half * (ctro2 (i, k, n) + ctro2 (i, kk, n))
                    ! dv3q = half * (ecdo2 (i, k, n) + ecdo2 (i, kk, n))

                    tem = half * (xlamue (i, k) + xlamue (i, kk))
                    !tem1 = half * (xlamud (i, k) + xlamud (i, kk))
                    tem1 = half * (xlamud (i) + xlamud (i))

                    ! if (k <= kbcon (i)) then
                    ! ptem = xlamde
                    ! ptem1 = xlamd (i) + xlamdd
                    ! else
                    ! ptem = xlamde
                    ! ptem1 = xlamdd
                    ! endif

                    dellae2 (i, k, n) = dellae2 (i, k, n) + &

                    ! -----------------------------------------------------------------------
                    ! detrainment from updraft
                    ! -----------------------------------------------------------------------

                        (aup * tem1 * eta (i, kk) * dv1q &

                    ! -----------------------------------------------------------------------
                    ! entrainement into up and downdraft
                    ! -----------------------------------------------------------------------

                        ! - (aup * tem * eta (i, kk) + adw * edto (i) * ptem * etad (i, k)) * dv2q &
                        - (aup * tem * eta (i, kk)) * dv2q &

                    ! -----------------------------------------------------------------------
                    ! detrainment from downdraft
                    ! -----------------------------------------------------------------------

                        ! + (adw * edto (i) * ptem1 * etad (i, k) * dv3q) &
                        ) * dz * xmbp (i, k)

                    wet_dep (i, k, n) = chem_pw (i, k, n) * grav / delp (i, k)

                    ! -----------------------------------------------------------------------
                    ! sinks from where updraft and downdraft start
                    ! -----------------------------------------------------------------------

                    ! if (k == jmin (i) + 1) then
                        ! dellae2 (i, k, n) = dellae2 (i, k, n) &
                            ! - edto (i) * etad (i, kk) * ctro2 (i, kk, n) * xmbp (i, k)
                    ! endif

                    if (k == kb (i)) then
                        dellae2 (i, k, n) = dellae2 (i, k, n) &
                            - eta (i, k) * ctro2 (i, k, n) * xmbp (i, k)
                    endif
                endif
            enddo
        enddo

        do i = 1, im
            if (cnvflg (i)) then
                if (kb (i) == 1) then
                    k = kb (i)
                    dellae2 (i, k, n) = dellae2 (i, k, n) &
                        - eta (i, k) * ctro2 (i, k, n) * xmbp (i, k)
                endif
            endif
        enddo

    enddo

    ! -----------------------------------------------------------------------
    ! for every tracer...
    ! -----------------------------------------------------------------------

    do n = 1, ntc
        flx_lo = zero
        totlout = zero
        clipout = zero

        ! -----------------------------------------------------------------------
        ! compute low - order mass flux, upstream
        ! -----------------------------------------------------------------------

        do k = 2, km1
            kk = k - 1
            do i = 1, im
                if (cnvflg (i) .and. (kk < ktcon (i))) then
                    tem = zero
                    if (kk >= kb (i)) tem = eta (i, kk)
                    ! if (kk <= jmin (i)) tem = tem - edto (i) * etad (i, kk)

                    ! -----------------------------------------------------------------------
                    ! low - order flux, upstream
                    ! -----------------------------------------------------------------------

                    if (tem > zero) then
                        flx_lo (i, k) = - xmb (i) * tem * qaero (i, k, n)
                    elseif (tem < zero) then
                        flx_lo (i, k) = - xmb (i) * tem * qaero (i, kk, n)
                    endif
                endif
            enddo
        enddo

        ! -----------------------------------------------------------------------
        ! make sure low - ord fluxes don't violate positive - definiteness
        ! -----------------------------------------------------------------------

        do k = 1, km1
            kp1 = k + 1
            do i = 1, im
                if (cnvflg (i) .and. (k <= ktcon (i))) then

                    ! -----------------------------------------------------------------------
                    ! time step / grid spacing
                    ! -----------------------------------------------------------------------

                    dtovdz = grav * dtime_max / abs (delp (i, k))

                    ! -----------------------------------------------------------------------
                    ! total flux out
                    ! -----------------------------------------------------------------------

                    totlout (i, k) = max (zero, flx_lo (i, kp1)) - min (zero, flx_lo (i, k))
                    clipout (i, k) = min (one, qaero (i, k, n) / max (epsil, totlout (i, k)) &
                        / (1.0001 * dtovdz))
                endif
            enddo
        enddo

        ! -----------------------------------------------------------------------
        ! recompute upstream mass fluxes
        ! -----------------------------------------------------------------------

        do k = 2, km1
            kk = k - 1
            do i = 1, im
                if (cnvflg (i) .and. (kk < ktcon (i))) then
                    tem = zero
                    if (kk >= kb (i)) tem = eta (i, kk)
                    ! if (kk <= jmin (i)) tem = tem - edto (i) * etad (i, kk)
                    if (tem > zero) then
                        flx_lo (i, k) = flx_lo (i, k) * clipout (i, k)
                    elseif (tem < zero) then
                        flx_lo (i, k) = flx_lo (i, k) * clipout (i, kk)
                    endif
                endif
            enddo
        enddo

        ! -----------------------------------------------------------------------
        ! a positive - definite low - order (diffusive) solution for the subsidnce fluxes
        ! -----------------------------------------------------------------------

        do k = 1, km1
            kp1 = k + 1
            do i = 1, im
                if (cnvflg (i) .and. (k <= ktcon (i))) then
                    dtovdz = grav * dtime_max / abs (delp (i, k)) ! time step / grid spacing
                    dellae2 (i, k, n) = dellae2 (i, k, n) &
                        - (flx_lo (i, kp1) - flx_lo (i, k)) * dtovdz / dtime_max
                endif
            enddo
        enddo

    enddo ! ctr

    ! -----------------------------------------------------------------------
    ! convert wet deposition to total mass deposited over dt and dp
    ! -----------------------------------------------------------------------

    do n = 1, ntc
        do k = 1, km
            do i = 1, im
                if (cnvflg (i) .and. (k < ktcon (i))) &
                    wet_dep (i, k, n) = wet_dep (i, k, n) * xmb (i) * delt * delp (i, k)
            enddo
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute final aerosol concentrations
    ! -----------------------------------------------------------------------

    do n = 1, ntc
        do k = 1, km
            do i = 1, im
                if (cnvflg (i) .and. (k <= min (kmax (i), ktcon (i)))) then
                    qaero (i, k, n) = qaero (i, k, n) + dellae2 (i, k, n) * delt
                    if (qaero (i, k, n) < zero) then

                        ! -----------------------------------------------------------------------
                        ! add negative mass to wet deposition
                        ! -----------------------------------------------------------------------

                        wet_dep (i, k, n) = wet_dep (i, k, n) - qaero (i, k, n) * delp (i, k)
                        qaero (i, k, n) = qcmin
                    endif
                endif
            enddo
        enddo
    enddo

end subroutine sa_aamf_shal_aero

end module sa_aamf_mod
