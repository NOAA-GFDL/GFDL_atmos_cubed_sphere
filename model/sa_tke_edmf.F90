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
! Scale-Aware Turbulent-Kinetic-Energy based Moist-Eddy-Diffusivity-Mass-Flux
! (SA-TKE-EDMF) Subgrid Vertical Turbulence Mixing Scheme
! For the convective boundary layer, the scheme adopts EDMF parameterization
! (Siebesma et al. 2007) to take into account non-local transport by
! large eddies (mfpblt.f).
! A new mass-flux parameterizaiton for stratocumulus-top-induced turbulence
! mixing has been introduced (previously, it was eddy diffusion form) (mfscu.f).
! For local turbulence mixing, a TKE closure model is used.
! Developers: Jongil Han, Kun Gao, Linjiong Zhou, and the GFDL FV3 Team
! References: Han et al. (2016), Han and Bretherton (2019)
! =======================================================================

! =======================================================================
! Updates at GFDL:
! 1) Jul 2019 by Kun Gao
!    goal: to allow for tke advection
!    change: rearange tracers (q1g)
!    TKE no longer needs to be the last tracer
! 2) Nov 2019 by Kun Gao
!    turn off non-local mixing for hydrometers to avoid unphysical negative values
! 3) Jun 2020 by Kun Gao
!    a) add option for turning off upper-limter on background diff. in inversion layer
!       over land/ice points (cap_k0_land)
!    b) use different xkzm_m, xkzm_h for land, ocean and sea ice points
!    c) add option for turning off hb19 formula for surface backgroud diff. (do_dk_hb19)
! 4) May 2022 by Linjiong Zhou
!    put it into the FV3 dynamical core and revise accordingly
! =======================================================================

module sa_tke_edmf_mod

    use fms_mod, only: check_nml_error
    use gfdl_mp_mod, only: mqs

    implicit none

    private

    ! -----------------------------------------------------------------------
    ! public subroutines, functions, and variables
    ! -----------------------------------------------------------------------

    public :: sa_tke_edmf_init
    public :: sa_tke_edmf_pbl
    public :: sa_tke_edmf_sfc

    ! -----------------------------------------------------------------------
    ! physics constants
    ! -----------------------------------------------------------------------

    real, parameter :: grav = 9.80665 ! acceleration due to gravity (m/s^2), ref: IFS

    real, parameter :: sbc = 5.670400e-8 ! Stefan-Boltzmann constant (kg/s^3/K^4)

    real, parameter :: rdgas = 287.05 ! gas constant for dry air (J/kg/K): ref: GFDL, GFS
    real, parameter :: rvgas = 461.50 ! gas constant for water vapor (J/kg/K): ref: GFDL, GFS

    real, parameter :: zvir = rvgas / rdgas - 1. ! 0.6077667316114637
    real, parameter :: eps = rdgas / rvgas ! 0.6219934994582882
    real, parameter :: epsm1 = rdgas / rvgas - 1. ! -0.3780065005417118

    real, parameter :: t0ice = 273.15 ! freezing temperature (K): ref: GFDL, GFS
    real, parameter :: tgice = 271.2 ! freezing temperature at sea (K)

    real, parameter :: cp_air = 1004.6 ! heat capacity of dry air at constant pressure (J/kg/K): ref: GFDL, GFS
    real, parameter :: cp_vap = 4.0 * rvgas ! 1846.0885419672554, heat capacity of water vapor at constnat pressure (J/kg/K)

    real, parameter :: c_liq = 4.218e3 ! heat capacity of water at 0 deg C (J/kg/K), ref: IFS

    real, parameter :: hlv = 2.5e6 ! latent heat of evaporation at 0 deg C (J/kg): ref: GFDL, GFS
    real, parameter :: hlf = 3.3358e5 ! latent heat of fusion at 0 deg C (J/kg): ref: GFDL, GFS

    ! -----------------------------------------------------------------------
    ! namelist parameters
    ! -----------------------------------------------------------------------

    logical :: cap_k0_land = .true.  ! flag for applying limter on background diff in inversion
    logical :: do_dk_hb19  = .false. ! flag for using hb19 formula for background diff
    logical :: dspheat     = .false. ! flag for tke dissipative heating
    logical :: sfc_gfdl    = .false. ! flag for using updated sfc layer scheme

    logical :: redrag              = .false. ! flag for reduced drag coeff. over sea
    logical :: do_z0_moon          = .false. ! flag for using z0 scheme in Moon et al. 2007
    logical :: do_z0_hwrf15        = .false. ! flag for using z0 scheme in 2015 HWRF
    logical :: do_z0_hwrf17        = .false. ! flag for using z0 scheme in 2017 HWRF
    logical :: do_z0_hwrf17_hwonly = .false. ! flag for using z0 scheme in 2017 HWRF only under high wind

    integer :: ivegsrc = 2 ! ivegsrc = 0   => USGS,
                           ! ivegsrc = 1   => IGBP (20 category)
                           ! ivegsrc = 2   => UMD  (13 category)

    real :: xkzm_mo  = 1.0   ! bkgd_vdif_m background vertical diffusion for momentum over ocean
    real :: xkzm_ho  = 1.0   ! bkgd_vdif_h background vertical diffusion for heat q over ocean
    real :: xkzm_ml  = 1.0   ! bkgd_vdif_m background vertical diffusion for momentum over land
    real :: xkzm_hl  = 1.0   ! bkgd_vdif_h background vertical diffusion for heat q over land
    real :: xkzm_mi  = 1.0   ! bkgd_vdif_m background vertical diffusion for momentum over ice
    real :: xkzm_hi  = 1.0   ! bkgd_vdif_h background vertical diffusion for heat q over ice
    real :: xkzm_s   = 1.0   ! bkgd_vdif_s sigma threshold for background mom. diffusion
    real :: xkzm_lim = 0.01  ! background vertical diffusion limit
    real :: xkzm_fac = 1.0   ! background vertical diffusion factor
    real :: xkzinv   = 0.15  ! diffusivity in inversion layers
    real :: xkgdx    = 25.e3 ! background vertical diffusion threshold
    real :: rlmn     = 30.   ! lower-limter on asymtotic mixing length in satmedmfdiff.f
    real :: rlmx     = 300.  ! upper-limter on asymtotic mixing length in satmedmfdiff.f

    real :: czilc        = 0.8     ! Zilintkivitch constant
    real :: z0s_max      = .317e-2 ! a limiting value for z0 under high windskk
    real :: wind_th_hwrf = 33.     ! wind speed threshold when z0 level off as in HWRF

    real :: ck0 = 0.4  ! proportionality coefficient for momentum in PBL
    real :: ck1 = 0.15 ! proportionality coefficient for momentum above PBL
    real :: ch0 = 0.4  ! proportionality coefficient for heat & q in PBL
    real :: ch1 = 0.15 ! proportionality coefficient for heat & q above PBL

    ! -----------------------------------------------------------------------
    ! namelist
    ! -----------------------------------------------------------------------

    namelist / sa_tke_edmf_nml / &
        xkzm_mo, xkzm_ho, xkzm_ml, xkzm_hl, xkzm_mi, xkzm_hi, xkzm_s, &
        xkzm_lim, xkzm_fac, xkzinv, xkgdx, rlmn, rlmx, sfc_gfdl, &
        cap_k0_land, do_dk_hb19, dspheat, redrag, do_z0_moon, &
        do_z0_hwrf15, do_z0_hwrf17, do_z0_hwrf17_hwonly, czilc, &
        z0s_max, wind_th_hwrf, ivegsrc, ck0, ck1, ch0, ch1

contains

! =======================================================================
! SA-TKE-EDMF initialization
! =======================================================================

subroutine sa_tke_edmf_init (input_nml_file, logunit)

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

    read (input_nml_file, nml = sa_tke_edmf_nml, iostat = ios)
    ierr = check_nml_error (ios, 'sa_tke_edmf_nml')

    ! -----------------------------------------------------------------------
    ! write namelist to log file
    ! -----------------------------------------------------------------------

    write (logunit, *) " ================================================================== "
    write (logunit, *) "sa_tke_edmf_mod"
    write (logunit, nml = sa_tke_edmf_nml)

end subroutine sa_tke_edmf_init

! =======================================================================
! SA-TKE-EDMF scheme
! =======================================================================

subroutine sa_tke_edmf_pbl (im, km, ntrac, ntcw, ntiw, ntke, &
        delt, u1, v1, t1, q1, gsize, islimsk, &
        radh, rbsoil, zorl, u10m, v10m, fm, fh, &
        tsea, heat, evap, stress, spd1, kinver, &
        psk, del, prsi, prsl, prslk, phii, phil, &
        hpbl, kpbl, dusfc, dvsfc, dtsfc, dqsfc, dkt_out)

    implicit none

    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------

    integer, intent (in) :: im, km, ntrac, ntcw, ntiw, ntke
    integer, intent (in) :: kinver (im), islimsk (im)

    real, intent (in) :: delt
    real, intent (in) :: radh (im, km), gsize (im), &
        psk (im), rbsoil (im), &
        zorl (im), tsea (im), &
        u10m (im), v10m (im), &
        fm (im), fh (im), &
        evap (im), heat (im), &
        stress (im), spd1 (im), &
        prsi (im, km + 1), del (im, km), &
        prsl (im, km), prslk (im, km), &
        phii (im, km + 1), phil (im, km)

    real, intent (inout) :: u1 (im, km), v1 (im, km), &
        t1 (im, km), q1 (im, km, ntrac)

    integer, intent (out) :: kpbl (im)

    real, intent (out) :: hpbl (im)

    real, intent (out), optional :: dusfc (im), dvsfc (im), dtsfc (im), dqsfc (im), &
        dkt_out (im, km)

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: i, is, k, kk, n, km1, kmpbl, kmscu, ntrac1, ntcw_new
    integer :: lcld (im), kcld (im), krad (im), mrad (im)
    integer :: kx1 (im), kpblx (im)

    real :: tke (im, km), tkeh (im, km - 1)

    real :: theta (im, km), thvx (im, km), thlvx (im, km), &
        qlx (im, km), thetae (im, km), thlx (im, km), &
        slx (im, km), svx (im, km), qtx (im, km), &
        tvx (im, km), pix (im, km), radx (im, km - 1), &
        dku (im, km - 1), dkt (im, km - 1), dkq (im, km - 1), &
        cku (im, km - 1), ckt (im, km - 1), q1g (im, km, ntrac), &
        vdt (im, km), udt (im, km), tdt (im, km), qdt (im, km)

    real :: plyr (im, km), rhly (im, km), cfly (im, km), &
        qstl (im, km)

    real :: dtdz1 (im), gdx (im), &
        phih (im), phim (im), prn (im, km - 1), &
        rbdn (im), rbup (im), thermal (im), &
        ustar (im), wstar (im), hpblx (im), &
        ust3 (im), wst3 (im), &
        z0 (im), crb (im), &
        hgamt (im), hgamq (im), &
        wscale (im), vpert (im), &
        zol (im), sflux (im), radj (im), &
        tx1 (im), tx2 (im)

    real :: radmin (im)

    real :: zi (im, km + 1), zl (im, km), zm (im, km), &
        xkzo (im, km - 1), xkzmo (im, km - 1), &
        xkzm_hx (im), xkzm_mx (im), &
        rdzt (im, km - 1), &
        al (im, km - 1), ad (im, km), au (im, km - 1), &
        f1 (im, km), f2 (im, km * (ntrac - 1))

    real :: elm (im, km), ele (im, km), rle (im, km - 1), &
        ckz (im, km), chz (im, km), &
        diss (im, km - 1), prod (im, km - 1), &
        bf (im, km - 1), shr2 (im, km - 1), &
        xlamue (im, km - 1), xlamde (im, km - 1), &
        gotvx (im, km), rlam (im, km - 1)

    ! variables for updrafts (thermals)
    real :: tcko (im, km), qcko (im, km, ntrac), &
        ucko (im, km), vcko (im, km), &
        buou (im, km), xmf (im, km)

    ! variables for stratocumulus - top induced downdrafts
    real :: tcdo (im, km), qcdo (im, km, ntrac), &
        ucdo (im, km), vcdo (im, km), &
        buod (im, km), xmfd (im, km)

    logical :: pblflg (im), sfcflg (im), flg (im)
    logical :: scuflg (im), pcnvflg (im)
    logical :: mlenflg

    ! pcnvflg: true for unstable pbl
    real :: aphi16, aphi5, &
        wfac, cfac, &
        gamcrt, gamcrq, sfcfrac, &
        conq, cont, conw, &
        dsdz2, dsdzt, dkmax, &
        dsig, dt2, dtodsd, &
        dtodsu, g, factor, dz, &
        gocp, gravi, zol1, zolcru, &
        buop, shrp, dtn, cdtn, &
        prnum, prmax, prmin, prtke, &
        prscu, dw2, dw2min, zk, &
        elmfac, elefac, dspmax, &
        alp, clwt, cql, &
        f0, robn, crbmin, crbmax, &
        es, qs, value, onemrh, &
        cfh, gamma, elocp, el2orc, &
        epsi, beta, chx, cqx, &
        rdt, rdz, qmin, qlmin, &
        ri, rimin, &
        rbcr, rbint, tdzmin, &
        elmx, &
        ttend, utend, vtend, qtend, &
        zfac, zfmin, vk, spdk2, &
        tkmin, dspfac, &
        zlup, zldn, bsum, &
        tem, tem1, tem2, &
        ptem, ptem0, ptem1, ptem2

    real :: ce0, rchck

    real :: qlcr, zstblmax

    real :: h1

    parameter (gravi = 1.0 / grav)
    parameter (g = grav)
    parameter (gocp = g / cp_air)
    parameter (cont = cp_air / g, conq = hlv / g, conw = 1.0 / g)
    parameter (elocp = hlv / cp_air, el2orc = hlv * hlv / (rvgas * cp_air))
    parameter (wfac = 7.0, cfac = 4.5)
    parameter (gamcrt = 3., gamcrq = 0., sfcfrac = 0.1)
    parameter (vk = 0.4, rimin = - 100.)
    parameter (rbcr = 0.25, zolcru = - 0.02, tdzmin = 1.e-3)
    parameter (prmin = 0.25, prmax = 4.0, prtke = 1.0, prscu = 0.67)
    parameter (f0 = 1.e-4, crbmin = 0.15, crbmax = 0.35)
    parameter (tkmin = 1.e-9, dspfac = 0.5, dspmax = 10.0)
    parameter (qmin = 1.e-8, qlmin = 1.e-12, zfmin = 1.e-8)
    parameter (aphi5 = 5., aphi16 = 16.)
    parameter (elmfac = 1.0, elefac = 1.0, cql = 100.)
    parameter (dw2min = 1.e-4, dkmax = 1000.)
    parameter (qlcr = 3.5e-5, zstblmax = 2500.)
    parameter (h1 = 0.33333333)
    parameter (ce0 = 0.4)
    parameter (rchck = 1.5, cdtn = 25.)

    elmx = rlmx

    ! -----------------------------------------------------------------------
    ! kgao note (jul 2019)
    ! the code was originally written assuming ntke = ntrac
    ! in this version ntke does not need to be equal to ntrac
    ! in the following we rearrange q1g so that tke is the last tracer
    ! -----------------------------------------------------------------------

    !if (ntrac >= 3) then
    if (ntke == ntrac) then ! tke is the last tracer
        q1g (:, :, :) = q1 (:, :, :)
    else ! tke is not
        do kk = 1, ntke - 1
            q1g (:, :, kk) = q1 (:, :, kk)
        enddo
        do kk = ntke + 1, ntrac
            q1g (:, :, kk - 1) = q1 (:, :, kk)
        enddo
        q1g (:, :, ntrac) = q1 (:, :, ntke)
    endif
    !endif

    dt2 = delt
    rdt = 1. / dt2

    ntrac1 = ntrac - 1
    km1 = km - 1
    kmpbl = km / 2
    kmscu = km / 2

    do k = 1, km
        do i = 1, im
            zi (i, k) = phii (i, k) * gravi
            zl (i, k) = phil (i, k) * gravi
            xmf (i, k) = 0.
            xmfd (i, k) = 0.
            buou (i, k) = 0.
            buod (i, k) = 0.
            ckz (i, k) = ck1
            chz (i, k) = ch1
        enddo
    enddo

    do i = 1, im
        zi (i, km + 1) = phii (i, km + 1) * gravi
    enddo
    do k = 1, km
        do i = 1, im
            zm (i, k) = zi (i, k + 1)
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! horizontal grid size
    ! -----------------------------------------------------------------------

    do i = 1, im
        gdx (i) = gsize (i)
    enddo

    do k = 1, km
        do i = 1, im
            tke (i, k) = max (q1 (i, k, ntke), tkmin) ! tke at layer centers
        enddo
    enddo
    do k = 1, km1
        do i = 1, im
            tkeh (i, k) = 0.5 * (tke (i, k) + tke (i, k + 1)) ! tke at interfaces
        enddo
    enddo

    do k = 1, km1
        do i = 1, im
            rdzt (i, k) = 1.0 / (zl (i, k + 1) - zl (i, k))
            prn (i, k) = 1.0
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! han and bretherton, 2019
    ! set background diffusivities as a function of
    ! horizontal grid size with xkzm_h & xkzm_m for gdx >= xkgdx
    ! and 0.01 for gdx = 5m, i.e.,
    ! xkzm_hx = 0.01 + (xkzm_h - 0.01) / (xkgdx - 5.) * (gdx - 5.)
    ! xkzm_mx = 0.01 + (xkzm_h - 0.01) / (xkgdx - 5.) * (gdx - 5.)
    ! -----------------------------------------------------------------------

    do i = 1, im
        kx1 (i) = 1
        tx1 (i) = 1.0 / prsi (i, 1)
        tx2 (i) = tx1 (i)

        ! -----------------------------------------------------------------------
        ! kgao change - set surface value of background diff (dk) below
        ! -----------------------------------------------------------------------

        if (do_dk_hb19) then ! use eq43 in hb2019

            if (gdx (i) >= xkgdx) then ! resolution coarser than xkgdx
                if (islimsk (i) == 1) then ! land points
                    xkzm_hx (i) = xkzm_hl
                    xkzm_mx (i) = xkzm_ml
                elseif (islimsk (i) == 2) then! sea ice points
                    xkzm_hx (i) = xkzm_hi
                    xkzm_mx (i) = xkzm_mi
                else ! ocean points
                    xkzm_hx (i) = xkzm_ho
                    xkzm_mx (i) = xkzm_mo
                endif
            else ! resolution finer than xkgdx
                tem = 1. / (xkgdx - 5.)
                if (islimsk (i) == 1) then ! land points
                    tem1 = (xkzm_hl - xkzm_lim) * tem
                    tem2 = (xkzm_ml - xkzm_lim) * tem
                elseif (islimsk (i) == 2) then! sea ice points
                    tem1 = (xkzm_hi - xkzm_lim) * tem
                    tem2 = (xkzm_mi - xkzm_lim) * tem
                else ! ocean points
                    tem1 = (xkzm_ho - xkzm_lim) * tem
                    tem2 = (xkzm_mo - xkzm_lim) * tem
                endif
                ptem = gdx (i) - 5.
                xkzm_hx (i) = xkzm_lim + tem1 * ptem
                xkzm_mx (i) = xkzm_lim + tem2 * ptem
            endif

        else ! use values in the namelist; no res dependency

            if (islimsk (i) == 1) then ! land points
                xkzm_hx (i) = xkzm_hl
                xkzm_mx (i) = xkzm_ml
            elseif (islimsk (i) == 2) then ! sea ice points
                xkzm_hx (i) = xkzm_hi
                xkzm_mx (i) = xkzm_mi
            else ! ocean points
                xkzm_hx (i) = xkzm_ho
                xkzm_mx (i) = xkzm_mo
            endif

        endif
    enddo

    do k = 1, km1
        do i = 1, im
            xkzo (i, k) = 0.0
            xkzmo (i, k) = 0.0
            if (k < kinver (i)) then
                ! -----------------------------------------------------------------------
                ! vertical background diffusivity
                ! -----------------------------------------------------------------------
                ptem = prsi (i, k + 1) * tx1 (i)
                tem1 = 1.0 - ptem
                tem1 = tem1 * tem1 * 10.0
                xkzo (i, k) = xkzm_hx (i) * min (1.0, exp (- tem1))
                ! -----------------------------------------------------------------------
                ! vertical background diffusivity for momentum
                ! -----------------------------------------------------------------------
                if (ptem >= xkzm_s) then
                    xkzmo (i, k) = xkzm_mx (i)
                    kx1 (i) = k + 1
                else
                    if (k == kx1 (i) .and. k > 1) tx2 (i) = 1.0 / prsi (i, k)
                    tem1 = 1.0 - prsi (i, k + 1) * tx2 (i)
                    tem1 = tem1 * tem1 * 5.0
                    xkzmo (i, k) = xkzm_mx (i) * min (1.0, exp (- tem1))
                endif
            endif
        enddo
    enddo

    do i = 1, im
        z0 (i) = 0.01 * zorl (i)
        if (present (dusfc)) dusfc (i) = 0.
        if (present (dvsfc)) dvsfc (i) = 0.
        if (present (dtsfc)) dtsfc (i) = 0.
        if (present (dqsfc)) dqsfc (i) = 0.
        kpbl (i) = 1
        hpbl (i) = 0.
        kpblx (i) = 1
        hpblx (i) = 0.
        pblflg (i) = .true.
        sfcflg (i) = .true.
        if (rbsoil (i) > 0.) sfcflg (i) = .false.
        pcnvflg (i) = .false.
        scuflg (i) = .true.
        if (scuflg (i)) then
            radmin (i) = 0.
            mrad (i) = km1
            krad (i) = 1
            lcld (i) = km1
            kcld (i) = km1
        endif
    enddo

    do k = 1, km
        do i = 1, im
            pix (i, k) = psk (i) / prslk (i, k)
            theta (i, k) = t1 (i, k) * pix (i, k)
            if (ntiw > 0) then
                tem = max (q1 (i, k, ntcw), qlmin)
                tem1 = max (q1 (i, k, ntiw), qlmin)
                qlx (i, k) = tem + tem1
                ptem = hlv * tem + (hlv + hlf) * tem1
                slx (i, k) = cp_air * t1 (i, k) + phil (i, k) - ptem
            else
                qlx (i, k) = max (q1 (i, k, ntcw), qlmin)
                slx (i, k) = cp_air * t1 (i, k) + phil (i, k) - hlv * qlx (i, k)
            endif
            tem2 = 1. + zvir * max (q1g (i, k, 1), qmin) - qlx (i, k)
            thvx (i, k) = theta (i, k) * tem2
            tvx (i, k) = t1 (i, k) * tem2
            qtx (i, k) = max (q1g (i, k, 1), qmin) + qlx (i, k)
            thlx (i, k) = theta (i, k) - pix (i, k) * elocp * qlx (i, k)
            thlvx (i, k) = thlx (i, k) * (1. + zvir * qtx (i, k))
            svx (i, k) = cp_air * tvx (i, k)
            ptem1 = elocp * pix (i, k) * max (q1g (i, k, 1), qmin)
            thetae (i, k) = theta (i, k) + ptem1
            gotvx (i, k) = g / tvx (i, k)
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! the background vertical diffusivities in the inversion layers are limited
    ! to be less than or equal to xkzminv
    ! -----------------------------------------------------------------------

    do k = 1, km1
        do i = 1, im
            tem1 = (tvx (i, k + 1) - tvx (i, k)) * rdzt (i, k)

            if (cap_k0_land) then
                if (tem1 > 1.e-5) then
                    xkzo (i, k) = min (xkzo (i, k), xkzinv)
                    xkzmo (i, k) = min (xkzmo (i, k), xkzinv)
                endif
            else
                ! -----------------------------------------------------------------------
                ! kgao note: do not apply upper - limiter over land and sea ice points
                ! (consistent with change in satmedmfdifq.f in jun 2020)
                ! -----------------------------------------------------------------------
                if (tem1 > 0. .and. islimsk (i) == 0) then
                    xkzo (i, k) = min (xkzo (i, k), xkzinv)
                    xkzmo (i, k) = min (xkzmo (i, k), xkzinv)
                endif
            endif

        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute an empirical cloud fraction based on
    ! xu & randall's (1996, jas) study
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            plyr (i, k) = 0.01 * prsl (i, k) ! pa to mb (hpa)
            ! compute relative humidity
            es = 0.01 * mqs (t1 (i, k)) ! mqs in pa
            ! revise it to a stable format -- Linjiong Zhou, 7/19/2022
            ! qs = max (qmin, eps * es / (plyr (i, k) + epsm1 * es))
            qs = max (qmin, es / plyr (i, k) * eps * (1 + zvir * q1g (i, k, 1)))
            rhly (i, k) = max (0.0, min (1.0, max (qmin, q1g (i, k, 1)) / qs))
            qstl (i, k) = qs
        enddo
    enddo

    do k = 1, km
        do i = 1, im
            cfly (i, k) = 0.
            clwt = 1.0e-6 * (plyr (i, k) * 0.001)
            if (qlx (i, k) > clwt) then
                onemrh = max (1.e-10, 1.0 - rhly (i, k))
                tem1 = min (max ((onemrh * qstl (i, k)) ** 0.49, 0.0001), 1.0)
                tem1 = cql / tem1
                value = max (min (tem1 * qlx (i, k), 50.0), 0.0)
                tem2 = sqrt (sqrt (rhly (i, k)))
                cfly (i, k) = min (max (tem2 * (1.0 - exp (- value)), 0.0), 1.0)
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute buoyancy modified by clouds
    ! -----------------------------------------------------------------------

    do k = 1, km1
        do i = 1, im
            tem = 0.5 * (svx (i, k) + svx (i, k + 1))
            tem1 = 0.5 * (t1 (i, k) + t1 (i, k + 1))
            tem2 = 0.5 * (qstl (i, k) + qstl (i, k + 1))
            cfh = min (cfly (i, k + 1), 0.5 * (cfly (i, k) + cfly (i, k + 1)))
            alp = g / tem
            gamma = el2orc * tem2 / (tem1 ** 2)
            epsi = tem1 / elocp
            beta = (1. + gamma * epsi * (1. + zvir)) / (1. + gamma)
            chx = cfh * alp * beta + (1. - cfh) * alp
            cqx = cfh * alp * hlv * (beta - epsi)
            cqx = cqx + (1. - cfh) * zvir * g
            ptem1 = (slx (i, k + 1) - slx (i, k)) * rdzt (i, k)
            ptem2 = (qtx (i, k + 1) - qtx (i, k)) * rdzt (i, k)
            bf (i, k) = chx * ptem1 + cqx * ptem2
        enddo
    enddo

    do k = 1, km1
        do i = 1, im
            dku (i, k) = 0.
            dkt (i, k) = 0.
            dkq (i, k) = 0.
            cku (i, k) = 0.
            ckt (i, k) = 0.
            tem = zi (i, k + 1) - zi (i, k)
            radx (i, k) = tem * radh (i, k)
        enddo
    enddo

    do i = 1, im
        sflux (i) = heat (i) + evap (i) * zvir * theta (i, 1)
        if (.not.sfcflg (i) .or. sflux (i) <= 0.) pblflg (i) = .false.
    enddo

    ! -----------------------------------------------------------------------
    ! compute critical bulk richardson number
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (pblflg (i)) then
            ! thermal (i) = thvx (i, 1)
            thermal (i) = thlvx (i, 1)
            crb (i) = rbcr
        else
            thermal (i) = tsea (i) * (1. + zvir * max (q1g (i, 1, 1), qmin))
            tem = sqrt (u10m (i) ** 2 + v10m (i) ** 2)
            tem = max (tem, 1.)
            robn = tem / (f0 * z0 (i))
            tem1 = 1.e-7 * robn
            crb (i) = 0.16 * (tem1 ** (- 0.18))
            crb (i) = max (min (crb (i), crbmax), crbmin)
        endif
    enddo

    do i = 1, im
        dtdz1 (i) = dt2 / (zi (i, 2) - zi (i, 1))
    enddo

    do i = 1, im
        ustar (i) = sqrt (stress (i))
    enddo

    ! -----------------------------------------------------------------------
    ! compute buoyancy (bf) and winshear square
    ! -----------------------------------------------------------------------

    do k = 1, km1
        do i = 1, im
            rdz = rdzt (i, k)
            ! bf (i, k) = gotvx (i, k) * (thvx (i, k + 1) - thvx (i, k)) * rdz
            dw2 = (u1 (i, k) - u1 (i, k + 1)) ** 2 + &
                (v1 (i, k) - v1 (i, k + 1)) ** 2
            shr2 (i, k) = max (dw2, dw2min) * rdz * rdz
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! find pbl height based on bulk richardson number (mrf pbl scheme)
    ! and also for diagnostic purpose
    ! -----------------------------------------------------------------------

    do i = 1, im
        flg (i) = .false.
        rbup (i) = rbsoil (i)
    enddo

    do k = 1, kmpbl
        do i = 1, im
            if (.not.flg (i)) then
                rbdn (i) = rbup (i)
                spdk2 = max ((u1 (i, k) ** 2 + v1 (i, k) ** 2), 1.)
                ! rbup (i) = (thvx (i, k) - thermal (i)) * &
                ! (g * zl (i, k) / thvx (i, 1)) / spdk2
                rbup (i) = (thlvx (i, k) - thermal (i)) * &
                    (g * zl (i, k) / thlvx (i, 1)) / spdk2
                kpblx (i) = k
                flg (i) = rbup (i) > crb (i)
            endif
        enddo
    enddo

    do i = 1, im
        if (kpblx (i) > 1) then
            k = kpblx (i)
            if (rbdn (i) >= crb (i)) then
                rbint = 0.
            elseif (rbup (i) <= crb (i)) then
                rbint = 1.
            else
                rbint = (crb (i) - rbdn (i)) / (rbup (i) - rbdn (i))
            endif
            hpblx (i) = zl (i, k - 1) + rbint * (zl (i, k) - zl (i, k - 1))
            if (hpblx (i) < zi (i, kpblx (i))) kpblx (i) = kpblx (i) - 1
        else
            hpblx (i) = zl (i, 1)
            kpblx (i) = 1
        endif
        hpbl (i) = hpblx (i)
        kpbl (i) = kpblx (i)
        if (kpbl (i) <= 1) pblflg (i) = .false.
    enddo

    ! -----------------------------------------------------------------------
    ! compute similarity parameters
    ! -----------------------------------------------------------------------

    do i = 1, im
        zol (i) = max (rbsoil (i) * fm (i) * fm (i) / fh (i), rimin)
        if (sfcflg (i)) then
            zol (i) = min (zol (i), - zfmin)
        else
            zol (i) = max (zol (i), zfmin)
        endif

        zol1 = zol (i) * sfcfrac * hpbl (i) / zl (i, 1)
        if (sfcflg (i)) then
            tem = 1.0 / (1. - aphi16 * zol1)
            phih (i) = sqrt (tem)
            phim (i) = sqrt (phih (i))
        else
            phim (i) = 1. + aphi5 * zol1
            phih (i) = phim (i)
        endif
    enddo

    do i = 1, im
        if (pblflg (i)) then
            if (zol (i) < zolcru) then
                pcnvflg (i) = .true.
            endif
            wst3 (i) = gotvx (i, 1) * sflux (i) * hpbl (i)
            wstar (i) = wst3 (i) ** h1
            ust3 (i) = ustar (i) ** 3.
            wscale (i) = (ust3 (i) + wfac * vk * wst3 (i) * sfcfrac) ** h1
            ptem = ustar (i) / aphi5
            wscale (i) = max (wscale (i), ptem)
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! compute a thermal excess
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (pcnvflg (i)) then
            hgamt (i) = heat (i) / wscale (i)
            hgamq (i) = evap (i) / wscale (i)
            vpert (i) = hgamt (i) + hgamq (i) * zvir * theta (i, 1)
            vpert (i) = max (vpert (i), 0.)
            tem = min (cfac * vpert (i), gamcrt)
            thermal (i) = thermal (i) + tem !jih jul2020
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! enhance the pbl height by considering the thermal excess
    ! (overshoot pbl top) -- jih jul2020
    ! -----------------------------------------------------------------------

    do i = 1, im
        flg (i) = .true.
        if (pcnvflg (i)) then
            flg (i) = .false.
            rbup (i) = rbsoil (i)
        endif
    enddo
    do k = 2, kmpbl
        do i = 1, im
            if (.not.flg (i)) then
                rbdn (i) = rbup (i)
                spdk2 = max ((u1 (i, k) ** 2 + v1 (i, k) ** 2), 1.)
                rbup (i) = (thlvx (i, k) - thermal (i)) * &
                    (g * zl (i, k) / thlvx (i, 1)) / spdk2
                kpbl (i) = k
                flg (i) = rbup (i) > crb (i)
            endif
        enddo
    enddo

    do i = 1, im
        if (pcnvflg (i)) then
            k = kpbl (i)
            if (rbdn (i) >= crb (i)) then
                rbint = 0.
            elseif (rbup (i) <= crb (i)) then
                rbint = 1.
            else
                rbint = (crb (i) - rbdn (i)) / (rbup (i) - rbdn (i))
            endif
            hpbl (i) = zl (i, k - 1) + rbint * (zl (i, k) - zl (i, k - 1))
            if (hpbl (i) < zi (i, kpbl (i))) then
                kpbl (i) = kpbl (i) - 1
            endif
            if (kpbl (i) <= 1) then
                pcnvflg (i) = .false.
                pblflg (i) = .false.
            endif
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! look for stratocumulus
    ! -----------------------------------------------------------------------

    do i = 1, im
        flg (i) = scuflg (i)
    enddo
    do k = 1, km1
        do i = 1, im
            if (flg (i) .and.zl (i, k) >= zstblmax) then
                lcld (i) = k
                flg (i) = .false.
            endif
        enddo
    enddo
    do i = 1, im
        flg (i) = scuflg (i)
    enddo
    do k = kmscu, 1, - 1
        do i = 1, im
            if (flg (i) .and. k <= lcld (i)) then
                if (qlx (i, k) >= qlcr) then
                    kcld (i) = k
                    flg (i) = .false.
                endif
            endif
        enddo
    enddo
    do i = 1, im
        if (scuflg (i) .and. kcld (i) == km1) scuflg (i) = .false.
    enddo

    do i = 1, im
        flg (i) = scuflg (i)
    enddo
    do k = kmscu, 1, - 1
        do i = 1, im
            if (flg (i) .and. k <= kcld (i)) then
                if (qlx (i, k) >= qlcr) then
                    if (radx (i, k) < radmin (i)) then
                        radmin (i) = radx (i, k)
                        krad (i) = k
                    endif
                else
                    flg (i) = .false.
                endif
            endif
        enddo
    enddo
    do i = 1, im
        if (scuflg (i) .and. krad (i) <= 1) scuflg (i) = .false.
        if (scuflg (i) .and. radmin (i) >= 0.) scuflg (i) = .false.
    enddo

    ! -----------------------------------------------------------------------
    ! compute components for mass flux mixing by large thermals
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (pcnvflg (i)) then
                tcko (i, k) = t1 (i, k)
                ucko (i, k) = u1 (i, k)
                vcko (i, k) = v1 (i, k)
            endif
            if (scuflg (i)) then
                tcdo (i, k) = t1 (i, k)
                ucdo (i, k) = u1 (i, k)
                vcdo (i, k) = v1 (i, k)
            endif
        enddo
    enddo
    do kk = 1, ntrac1
        do k = 1, km
            do i = 1, im
                if (pcnvflg (i)) then
                    qcko (i, k, kk) = q1g (i, k, kk)
                endif
                if (scuflg (i)) then
                    qcdo (i, k, kk) = q1g (i, k, kk)
                endif
            enddo
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! kgao note - change ntcw if q1g is rearranged
    ! -----------------------------------------------------------------------

    if (ntke > ntcw) then
        ntcw_new = ntcw
    else
        ntcw_new = ntcw - 1
    endif

    ! -----------------------------------------------------------------------
    ! edmf parameterization siebesma et al. (2007)
    ! -----------------------------------------------------------------------

    call mfpblt (im, km, kmpbl, ntcw_new, ntrac1, dt2, &
        pcnvflg, zl, zm, q1g, t1, u1, v1, plyr, pix, thlx, thvx, &
        gdx, hpbl, kpbl, vpert, buou, xmf, &
        tcko, qcko, ucko, vcko, xlamue)

    ! -----------------------------------------------------------------------
    ! mass - flux parameterization for stratocumulus - top - induced turbulence mixing
    ! -----------------------------------------------------------------------

    call mfscu (im, km, kmscu, ntcw_new, ntrac1, dt2, &
        scuflg, zl, zm, q1g, t1, u1, v1, plyr, pix, &
        thlx, thvx, thlvx, gdx, thetae, radj, &
        krad, mrad, radmin, buod, xmfd, &
        tcdo, qcdo, ucdo, vcdo, xlamde)

    ! -----------------------------------------------------------------------
    ! compute prandtl number and exchange coefficient varying with height
    ! -----------------------------------------------------------------------

    do k = 1, kmpbl
        do i = 1, im
            if (k < kpbl (i)) then
                tem = phih (i) / phim (i)
                ptem = - 3. * (max (zi (i, k + 1) - sfcfrac * hpbl (i), 0.)) ** 2. &
                     / hpbl (i) ** 2.
                if (pcnvflg (i)) then
                    prn (i, k) = 1. + (tem - 1.) * exp (ptem)
                else
                    prn (i, k) = tem
                endif
                prn (i, k) = min (prn (i, k), prmax)
                prn (i, k) = max (prn (i, k), prmin)

                ckz (i, k) = ck1 + (ck0 - ck1) * exp (ptem)
                ckz (i, k) = min (ckz (i, k), ck0)
                ckz (i, k) = max (ckz (i, k), ck1)
                chz (i, k) = ch1 + (ch0 - ch1) * exp (ptem)
                chz (i, k) = min (chz (i, k), ch0)
                chz (i, k) = max (chz (i, k), ch1)
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute an asymtotic mixing length
    ! -----------------------------------------------------------------------

    do k = 1, km1
        do i = 1, im
            zlup = 0.0
            bsum = 0.0
            mlenflg = .true.
            do n = k, km1
                if (mlenflg) then
                    dz = zl (i, n + 1) - zl (i, n)
                    ptem = gotvx (i, n) * (thvx (i, n + 1) - thvx (i, k)) * dz
                    ! ptem = gotvx (i, n) * (thlvx (i, n + 1) - thlvx (i, k)) * dz
                    bsum = bsum + ptem
                    zlup = zlup + dz
                    if (bsum >= tke (i, k)) then
                        if (ptem >= 0.) then
                            tem2 = max (ptem, zfmin)
                        else
                            tem2 = min (ptem, - zfmin)
                        endif
                        ptem1 = (bsum - tke (i, k)) / tem2
                        zlup = zlup - ptem1 * dz
                        zlup = max (zlup, 0.)
                        mlenflg = .false.
                    endif
                endif
            enddo
            zldn = 0.0
            bsum = 0.0
            mlenflg = .true.
            do n = k, 1, - 1
                if (mlenflg) then
                    if (n == 1) then
                        dz = zl (i, 1)
                        tem1 = tsea (i) * (1. + zvir * max (q1g (i, 1, 1), qmin))
                    else
                        dz = zl (i, n) - zl (i, n - 1)
                        tem1 = thvx (i, n - 1)
                        ! tem1 = thlvx (i, n - 1)
                    endif
                    ptem = gotvx (i, n) * (thvx (i, k) - tem1) * dz
                    ! ptem = gotvx (i, n) * (thlvx (i, k) - tem1) * dz
                    bsum = bsum + ptem
                    zldn = zldn + dz
                    if (bsum >= tke (i, k)) then
                        if (ptem >= 0.) then
                            tem2 = max (ptem, zfmin)
                        else
                            tem2 = min (ptem, - zfmin)
                        endif
                        ptem1 = (bsum - tke (i, k)) / tem2
                        zldn = zldn - ptem1 * dz
                        zldn = max (zldn, 0.)
                        mlenflg = .false.
                    endif
                endif
            enddo

            tem = 0.5 * (zi (i, k + 1) - zi (i, k))
            tem1 = min (tem, rlmn)

            ptem2 = min (zlup, zldn)
            rlam (i, k) = elmfac * ptem2
            rlam (i, k) = max (rlam (i, k), tem1)
            rlam (i, k) = min (rlam (i, k), rlmx)

            ptem2 = sqrt (zlup * zldn)
            ele (i, k) = elefac * ptem2
            ele (i, k) = max (ele (i, k), tem1)
            ele (i, k) = min (ele (i, k), elmx)

        enddo
    enddo

    do k = 1, km1
        do i = 1, im
            tem = vk * zl (i, k)
            if (zol (i) < 0.) then
                ptem = 1. - 100. * zol (i)
                ptem1 = ptem ** 0.2
                zk = tem * ptem1
            elseif (zol (i) >= 1.) then
                zk = tem / 3.7
            else
                ptem = 1. + 2.7 * zol (i)
                zk = tem / ptem
            endif
            elm (i, k) = zk * rlam (i, k) / (rlam (i, k) + zk)

            dz = zi (i, k + 1) - zi (i, k)
            tem = max (gdx (i), dz)
            elm (i, k) = min (elm (i, k), tem)
            ele (i, k) = min (ele (i, k), tem)

        enddo
    enddo
    do i = 1, im
        elm (i, km) = elm (i, km1)
        ele (i, km) = ele (i, km1)
    enddo

    ! -----------------------------------------------------------------------
    ! compute eddy diffusivities
    ! -----------------------------------------------------------------------

    do k = 1, km1
        do i = 1, im
            tem = 0.5 * (elm (i, k) + elm (i, k + 1))
            tem = tem * sqrt (tkeh (i, k))
            if (k < kpbl (i)) then
                if (pblflg (i)) then
                    dku (i, k) = ckz (i, k) * tem
                    dkt (i, k) = dku (i, k) / prn (i, k)
                else
                    dkt (i, k) = chz (i, k) * tem
                    dku (i, k) = dkt (i, k) * prn (i, k)
                endif
            else
                ri = max (bf (i, k) / shr2 (i, k), rimin)
                if (ri < 0.) then ! unstable regime
                    dku (i, k) = ck1 * tem
                    dkt (i, k) = rchck * dku (i, k)
                else ! stable regime
                    dkt (i, k) = ch1 * tem
                    prnum = 1.0 + 2.1 * ri
                    prnum = min (prnum, prmax)
                    dku (i, k) = dkt (i, k) * prnum
                endif
            endif

            if (scuflg (i)) then
                if (k >= mrad (i) .and. k < krad (i)) then
                    tem1 = ckz (i, k) * tem
                    ptem1 = tem1 / prscu
                    dku (i, k) = max (dku (i, k), tem1)
                    dkt (i, k) = max (dkt (i, k), ptem1)
                endif
            endif

            dkq (i, k) = prtke * dkt (i, k)

            dkt (i, k) = min (dkt (i, k), dkmax)
            dkt (i, k) = max (dkt (i, k), xkzo (i, k))
            dkq (i, k) = min (dkq (i, k), dkmax)
            dkq (i, k) = max (dkq (i, k), xkzo (i, k))
            dku (i, k) = min (dku (i, k), dkmax)
            dku (i, k) = max (dku (i, k), xkzmo (i, k))

        enddo
    enddo

    do i = 1, im
        if (scuflg (i)) then
            k = krad (i)
            tem = bf (i, k) / gotvx (i, k)
            tem1 = max (tem, tdzmin)
            ptem = radj (i) / tem1
            dkt (i, k) = dkt (i, k) + ptem
            dku (i, k) = dku (i, k) + ptem
            dkq (i, k) = dkq (i, k) + ptem
        endif
    enddo

    if (present (dkt_out)) then
        do k = 1, km1
            do i = 1, im
                dkt_out (i, k) = dkt (i, k)
            enddo
        enddo
    endif

    ! -----------------------------------------------------------------------
    ! compute buoyancy and shear productions of tke
    ! -----------------------------------------------------------------------

    do k = 1, km1
        do i = 1, im
            if (k == 1) then
                tem = - dkt (i, 1) * bf (i, 1)
                ! if (pcnvflg (i)) then
                ! ptem1 = xmf (i, 1) * buou (i, 1)
                ! else
                ptem1 = 0.
                ! endif
                if (scuflg (i) .and. mrad (i) == 1) then
                    ptem2 = xmfd (i, 1) * buod (i, 1)
                else
                    ptem2 = 0.
                endif
                tem = tem + ptem1 + ptem2
                buop = 0.5 * (gotvx (i, 1) * sflux (i) + tem)

                tem1 = dku (i, 1) * shr2 (i, 1)

                tem = (u1 (i, 2) - u1 (i, 1)) * rdzt (i, 1)
                ! if (pcnvflg (i)) then
                ! ptem = xmf (i, 1) * tem
                ! ptem1 = 0.5 * ptem * (u1 (i, 2) - ucko (i, 2))
                ! else
                ptem1 = 0.
                ! endif
                if (scuflg (i) .and. mrad (i) == 1) then
                    ptem = ucdo (i, 1) + ucdo (i, 2) - u1 (i, 1) - u1 (i, 2)
                    ptem = 0.5 * tem * xmfd (i, 1) * ptem
                else
                    ptem = 0.
                endif
                ptem1 = ptem1 + ptem

                tem = (v1 (i, 2) - v1 (i, 1)) * rdzt (i, 1)
                ! if (pcnvflg (i)) then
                ! ptem = xmf (i, 1) * tem
                ! ptem2 = 0.5 * ptem * (v1 (i, 2) - vcko (i, 2))
                ! else
                ptem2 = 0.
                ! endif
                if (scuflg (i) .and. mrad (i) == 1) then
                    ptem = vcdo (i, 1) + vcdo (i, 2) - v1 (i, 1) - v1 (i, 2)
                    ptem = 0.5 * tem * xmfd (i, 1) * ptem
                else
                    ptem = 0.
                endif
                ptem2 = ptem2 + ptem

                ! tem2 = stress (i) * spd1 (i) / zl (i, 1)
                tem2 = stress (i) * ustar (i) * phim (i) / (vk * zl (i, 1))
                shrp = 0.5 * (tem1 + ptem1 + ptem2 + tem2)
            else
                tem1 = - dkt (i, k - 1) * bf (i, k - 1)
                tem2 = - dkt (i, k) * bf (i, k)
                tem = 0.5 * (tem1 + tem2)
                if (pcnvflg (i) .and. k <= kpbl (i)) then
                    ptem = 0.5 * (xmf (i, k - 1) + xmf (i, k))
                    ptem1 = ptem * buou (i, k)
                else
                    ptem1 = 0.
                endif
                if (scuflg (i)) then
                    if (k >= mrad (i) .and. k < krad (i)) then
                        ptem0 = 0.5 * (xmfd (i, k - 1) + xmfd (i, k))
                        ptem2 = ptem0 * buod (i, k)
                    else
                        ptem2 = 0.
                    endif
                else
                    ptem2 = 0.
                endif
                buop = tem + ptem1 + ptem2

                tem1 = dku (i, k - 1) * shr2 (i, k - 1)
                tem2 = dku (i, k) * shr2 (i, k)
                tem = 0.5 * (tem1 + tem2)
                tem1 = (u1 (i, k + 1) - u1 (i, k)) * rdzt (i, k)
                tem2 = (u1 (i, k) - u1 (i, k - 1)) * rdzt (i, k - 1)
                if (pcnvflg (i) .and. k <= kpbl (i)) then
                    ptem = xmf (i, k) * tem1 + xmf (i, k - 1) * tem2
                    ptem1 = 0.5 * ptem * (u1 (i, k) - ucko (i, k))
                else
                    ptem1 = 0.
                endif
                if (scuflg (i)) then
                    if (k >= mrad (i) .and. k < krad (i)) then
                        ptem0 = xmfd (i, k) * tem1 + xmfd (i, k - 1) * tem2
                        ptem2 = 0.5 * ptem0 * (ucdo (i, k) - u1 (i, k))
                    else
                        ptem2 = 0.
                    endif
                else
                    ptem2 = 0.
                endif
                shrp = tem + ptem1 + ptem2
                tem1 = (v1 (i, k + 1) - v1 (i, k)) * rdzt (i, k)
                tem2 = (v1 (i, k) - v1 (i, k - 1)) * rdzt (i, k - 1)
                if (pcnvflg (i) .and. k <= kpbl (i)) then
                    ptem = xmf (i, k) * tem1 + xmf (i, k - 1) * tem2
                    ptem1 = 0.5 * ptem * (v1 (i, k) - vcko (i, k))
                else
                    ptem1 = 0.
                endif
                if (scuflg (i)) then
                    if (k >= mrad (i) .and. k < krad (i)) then
                        ptem0 = xmfd (i, k) * tem1 + xmfd (i, k - 1) * tem2
                        ptem2 = 0.5 * ptem0 * (vcdo (i, k) - v1 (i, k))
                    else
                        ptem2 = 0.
                    endif
                else
                    ptem2 = 0.
                endif
                shrp = shrp + ptem1 + ptem2
            endif
            prod (i, k) = buop + shrp
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! first predict tke due to tke production & dissipation (diss)
    ! -----------------------------------------------------------------------

    do k = 1, km1
        do i = 1, im
            rle (i, k) = ce0 / ele (i, k)
        enddo
    enddo
    kk = max (nint (dt2 / cdtn), 1)
    dtn = dt2 / float (kk)
    do n = 1, kk
        do k = 1, km1
            do i = 1, im
                tem = sqrt (tke (i, k))
                diss (i, k) = rle (i, k) * tke (i, k) * tem
                tem1 = prod (i, k) + tke (i, k) / dtn
                diss (i, k) = max (min (diss (i, k), tem1), 0.)
                tke (i, k) = tke (i, k) + dtn * (prod (i, k) - diss (i, k)) ! no diffusion yet
                tke (i, k) = max (tke (i, k), tkmin)
            enddo
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute updraft & downdraft properties for tke
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            if (pcnvflg (i)) then
                ! kgao change
                ! qcko (i, k, ntke) = tke (i, k)
                qcko (i, k, ntrac) = tke (i, k)
            endif
            if (scuflg (i)) then
                ! kgao change
                ! qcdo (i, k, ntke) = tke (i, k)
                qcdo (i, k, ntrac) = tke (i, k)
            endif
        enddo
    enddo
    do k = 2, kmpbl
        do i = 1, im
            if (pcnvflg (i) .and. k <= kpbl (i)) then
                dz = zl (i, k) - zl (i, k - 1)
                tem = 0.5 * xlamue (i, k - 1) * dz
                factor = 1. + tem
                ! kgao change
                ! qcko (i, k, ntke) = ((1. - tem) * qcko (i, k - 1, ntke) + tem * &
                ! (tke (i, k) + tke (i, k - 1))) / factor
                qcko (i, k, ntrac) = ((1. - tem) * qcko (i, k - 1, ntrac) + tem * &
                    (tke (i, k) + tke (i, k - 1))) / factor
            endif
        enddo
    enddo
    do k = kmscu, 1, - 1
        do i = 1, im
            if (scuflg (i) .and. k < krad (i)) then
                if (k >= mrad (i)) then
                    dz = zl (i, k + 1) - zl (i, k)
                    tem = 0.5 * xlamde (i, k) * dz
                    factor = 1. + tem
                    ! kgao change
                    ! qcdo (i, k, ntke) = ((1. - tem) * qcdo (i, k + 1, ntke) + tem * &
                    ! (tke (i, k) + tke (i, k + 1))) / factor
                    qcdo (i, k, ntrac) = ((1. - tem) * qcdo (i, k + 1, ntrac) + tem * &
                        (tke (i, k) + tke (i, k + 1))) / factor
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute tridiagonal matrix elements for turbulent kinetic energy
    ! -----------------------------------------------------------------------

    do i = 1, im
        ad (i, 1) = 1.0
        f1 (i, 1) = tke (i, 1)
    enddo

    do k = 1, km1
        do i = 1, im
            dtodsd = dt2 / del (i, k)
            dtodsu = dt2 / del (i, k + 1)
            dsig = prsl (i, k) - prsl (i, k + 1)
            rdz = rdzt (i, k)
            tem1 = dsig * dkq (i, k) * rdz
            dsdz2 = tem1 * rdz
            au (i, k) = - dtodsd * dsdz2
            al (i, k) = - dtodsu * dsdz2
            ad (i, k) = ad (i, k) - au (i, k)
            ad (i, k + 1) = 1. - al (i, k)
            tem2 = dsig * rdz

            if (pcnvflg (i) .and. k < kpbl (i)) then
                ptem = 0.5 * tem2 * xmf (i, k)
                ptem1 = dtodsd * ptem
                ptem2 = dtodsu * ptem
                tem = tke (i, k) + tke (i, k + 1)
                ! kgao change
                ! ptem = qcko (i, k, ntke) + qcko (i, k + 1, ntke)
                ptem = qcko (i, k, ntrac) + qcko (i, k + 1, ntrac)
                f1 (i, k) = f1 (i, k) - (ptem - tem) * ptem1
                f1 (i, k + 1) = tke (i, k + 1) + (ptem - tem) * ptem2
            else
                f1 (i, k + 1) = tke (i, k + 1)
            endif

            if (scuflg (i)) then
                if (k >= mrad (i) .and. k < krad (i)) then
                    ptem = 0.5 * tem2 * xmfd (i, k)
                    ptem1 = dtodsd * ptem
                    ptem2 = dtodsu * ptem
                    tem = tke (i, k) + tke (i, k + 1)
                    ! kgao change
                    ! ptem = qcdo (i, k, ntke) + qcdo (i, k + 1, ntke)
                    ptem = qcdo (i, k, ntrac) + qcdo (i, k + 1, ntrac)
                    f1 (i, k) = f1 (i, k) + (ptem - tem) * ptem1
                    f1 (i, k + 1) = f1 (i, k + 1) - (ptem - tem) * ptem2
                endif
            endif

        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! solve tridiagonal problem for tke
    ! -----------------------------------------------------------------------

    call tridit (im, km, 1, al, ad, au, f1, au, f1)

    ! -----------------------------------------------------------------------
    ! recover tendency of tke
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            ! fix negative tke
            f1 (i, k) = max (f1 (i, k), tkmin)
            q1g (i, k, ntrac) = f1 (i, k)
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute tridiagonal matrix elements for heat and moisture (and other tracers, except tke)
    ! -----------------------------------------------------------------------

    do i = 1, im
        ad (i, 1) = 1.
        f1 (i, 1) = t1 (i, 1) + dtdz1 (i) * heat (i)
        f2 (i, 1) = q1g (i, 1, 1) + dtdz1 (i) * evap (i)
    enddo
    if (ntrac1 >= 2) then
        do kk = 2, ntrac1
            is = (kk - 1) * km
            do i = 1, im
                f2 (i, 1 + is) = q1g (i, 1, kk)
            enddo
        enddo
    endif

    do k = 1, km1
        do i = 1, im
            dtodsd = dt2 / del (i, k)
            dtodsu = dt2 / del (i, k + 1)
            dsig = prsl (i, k) - prsl (i, k + 1)
            rdz = rdzt (i, k)
            tem1 = dsig * dkt (i, k) * rdz
            dsdzt = tem1 * gocp
            dsdz2 = tem1 * rdz
            au (i, k) = - dtodsd * dsdz2
            al (i, k) = - dtodsu * dsdz2
            ad (i, k) = ad (i, k) - au (i, k)
            ad (i, k + 1) = 1. - al (i, k)
            tem2 = dsig * rdz

            if (pcnvflg (i) .and. k < kpbl (i)) then
                ptem = 0.5 * tem2 * xmf (i, k)
                ptem1 = dtodsd * ptem
                ptem2 = dtodsu * ptem
                tem = t1 (i, k) + t1 (i, k + 1)
                ptem = tcko (i, k) + tcko (i, k + 1)
                f1 (i, k) = f1 (i, k) + dtodsd * dsdzt - (ptem - tem) * ptem1
                f1 (i, k + 1) = t1 (i, k + 1) - dtodsu * dsdzt + (ptem - tem) * ptem2
                tem = q1g (i, k, 1) + q1g (i, k + 1, 1)
                ptem = qcko (i, k, 1) + qcko (i, k + 1, 1)
                f2 (i, k) = f2 (i, k) - (ptem - tem) * ptem1
                f2 (i, k + 1) = q1g (i, k + 1, 1) + (ptem - tem) * ptem2
            else
                f1 (i, k) = f1 (i, k) + dtodsd * dsdzt
                f1 (i, k + 1) = t1 (i, k + 1) - dtodsu * dsdzt
                f2 (i, k + 1) = q1g (i, k + 1, 1)
            endif

            if (scuflg (i)) then
                if (k >= mrad (i) .and. k < krad (i)) then
                    ptem = 0.5 * tem2 * xmfd (i, k)
                    ptem1 = dtodsd * ptem
                    ptem2 = dtodsu * ptem
                    ptem = tcdo (i, k) + tcdo (i, k + 1)
                    tem = t1 (i, k) + t1 (i, k + 1)
                    f1 (i, k) = f1 (i, k) + (ptem - tem) * ptem1
                    f1 (i, k + 1) = f1 (i, k + 1) - (ptem - tem) * ptem2
                    tem = q1g (i, k, 1) + q1g (i, k + 1, 1)
                    ptem = qcdo (i, k, 1) + qcdo (i, k + 1, 1)
                    f2 (i, k) = f2 (i, k) + (ptem - tem) * ptem1
                    f2 (i, k + 1) = f2 (i, k + 1) - (ptem - tem) * ptem2
                endif
            endif
        enddo
    enddo

    if (ntrac1 >= 2) then
        do kk = 2, ntrac1
            is = (kk - 1) * km
            do k = 1, km1
                do i = 1, im
                    if (pcnvflg (i) .and. k < kpbl (i)) then
                        dtodsd = dt2 / del (i, k)
                        dtodsu = dt2 / del (i, k + 1)
                        dsig = prsl (i, k) - prsl (i, k + 1)
                        tem = dsig * rdzt (i, k)
                        ptem = 0.5 * tem * xmf (i, k)
                        ptem1 = dtodsd * ptem
                        ptem2 = dtodsu * ptem
                        tem1 = qcko (i, k, kk) + qcko (i, k + 1, kk)
                        tem2 = q1g (i, k, kk) + q1g (i, k + 1, kk)
                        ! kgao note - turn off non - local mixing
                        f2 (i, k + is) = f2 (i, k + is) ! - (tem1 - tem2) * ptem1
                        f2 (i, k + 1 + is) = q1g (i, k + 1, kk) ! + (tem1 - tem2) * ptem2
                    else
                        f2 (i, k + 1 + is) = q1g (i, k + 1, kk)
                    endif

                    if (scuflg (i)) then
                        if (k >= mrad (i) .and. k < krad (i)) then
                            dtodsd = dt2 / del (i, k)
                            dtodsu = dt2 / del (i, k + 1)
                            dsig = prsl (i, k) - prsl (i, k + 1)
                            tem = dsig * rdzt (i, k)
                            ptem = 0.5 * tem * xmfd (i, k)
                            ptem1 = dtodsd * ptem
                            ptem2 = dtodsu * ptem
                            tem1 = qcdo (i, k, kk) + qcdo (i, k + 1, kk)
                            tem2 = q1g (i, k, kk) + q1g (i, k + 1, kk)
                            ! kgao note - turn off non - local mixing
                            f2 (i, k + is) = f2 (i, k + is) ! + (tem1 - tem2) * ptem1
                            f2 (i, k + 1 + is) = f2 (i, k + 1 + is) ! - (tem1 - tem2) * ptem2
                        endif
                    endif

                enddo
            enddo
        enddo
    endif

    ! -----------------------------------------------------------------------
    ! solve tridiagonal problem for heat and moisture
    ! -----------------------------------------------------------------------

    call tridin (im, km, ntrac1, al, ad, au, f1, f2, au, f1, f2)

    ! -----------------------------------------------------------------------
    ! recover tendencies of heat and moisture
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            tdt (i, k) = (f1 (i, k) - t1 (i, k)) * rdt
            qdt (i, k) = (f2 (i, k) - q1g (i, k, 1)) * rdt
            if (present (dtsfc)) dtsfc (i) = dtsfc (i) + cont * del (i, k) * tdt (i, k)
            if (present (dqsfc)) dqsfc (i) = dqsfc (i) + conq * del (i, k) * qdt (i, k)
            t1 (i, k) = f1 (i, k)
            q1g (i, k, 1) = f2 (i, k)
        enddo
    enddo

    if (ntrac1 >= 2) then
        do kk = 2, ntrac1
            is = (kk - 1) * km
            do k = 1, km
                do i = 1, im
                    q1g (i, k, kk) = f2 (i, k + is)
                enddo
            enddo
        enddo
    endif

    ! -----------------------------------------------------------------------
    ! kgao note - rearrange tracer tendencies
    ! -----------------------------------------------------------------------

    !if (ntrac >= 3) then
    if (ntke == ntrac) then ! tke is the last tracer
        q1 (:, :, :) = q1g (:, :, :)
    else ! tke is not
        do kk = 1, ntke - 1
            q1 (:, :, kk) = q1g (:, :, kk)
        enddo
        q1 (:, :, ntke) = q1g (:, :, ntrac)
        do kk = ntke + 1, ntrac
            q1 (:, :, kk) = q1g (:, :, kk - 1)
        enddo
    endif
    !endif

    ! -----------------------------------------------------------------------
    ! add tke dissipative heating to temperature tendency
    ! -----------------------------------------------------------------------

    if (dspheat) then
        do k = 1, km1
            do i = 1, im
                ! tem = min (diss (i, k), dspmax)
                ! ttend = tem / cp_air
                ttend = diss (i, k) / cp_air
                t1 (i, k) = t1 (i, k) + dspfac * ttend * dt2
            enddo
        enddo
    endif

    ! -----------------------------------------------------------------------
    ! compute tridiagonal matrix elements for momentum
    ! -----------------------------------------------------------------------

    do i = 1, im
        ad (i, 1) = 1.0 + dtdz1 (i) * stress (i) / spd1 (i)
        f1 (i, 1) = u1 (i, 1)
        f2 (i, 1) = v1 (i, 1)
    enddo

    do k = 1, km1
        do i = 1, im
            dtodsd = dt2 / del (i, k)
            dtodsu = dt2 / del (i, k + 1)
            dsig = prsl (i, k) - prsl (i, k + 1)
            rdz = rdzt (i, k)
            tem1 = dsig * dku (i, k) * rdz
            dsdz2 = tem1 * rdz
            au (i, k) = - dtodsd * dsdz2
            al (i, k) = - dtodsu * dsdz2
            ad (i, k) = ad (i, k) - au (i, k)
            ad (i, k + 1) = 1. - al (i, k)
            tem2 = dsig * rdz

            if (pcnvflg (i) .and. k < kpbl (i)) then
                ptem = 0.5 * tem2 * xmf (i, k)
                ptem1 = dtodsd * ptem
                ptem2 = dtodsu * ptem
                tem = u1 (i, k) + u1 (i, k + 1)
                ptem = ucko (i, k) + ucko (i, k + 1)
                f1 (i, k) = f1 (i, k) - (ptem - tem) * ptem1
                f1 (i, k + 1) = u1 (i, k + 1) + (ptem - tem) * ptem2
                tem = v1 (i, k) + v1 (i, k + 1)
                ptem = vcko (i, k) + vcko (i, k + 1)
                f2 (i, k) = f2 (i, k) - (ptem - tem) * ptem1
                f2 (i, k + 1) = v1 (i, k + 1) + (ptem - tem) * ptem2
            else
                f1 (i, k + 1) = u1 (i, k + 1)
                f2 (i, k + 1) = v1 (i, k + 1)
            endif

            if (scuflg (i)) then
                if (k >= mrad (i) .and. k < krad (i)) then
                    ptem = 0.5 * tem2 * xmfd (i, k)
                    ptem1 = dtodsd * ptem
                    ptem2 = dtodsu * ptem
                    tem = u1 (i, k) + u1 (i, k + 1)
                    ptem = ucdo (i, k) + ucdo (i, k + 1)
                    f1 (i, k) = f1 (i, k) + (ptem - tem) * ptem1
                    f1 (i, k + 1) = f1 (i, k + 1) - (ptem - tem) * ptem2
                    tem = v1 (i, k) + v1 (i, k + 1)
                    ptem = vcdo (i, k) + vcdo (i, k + 1)
                    f2 (i, k) = f2 (i, k) + (ptem - tem) * ptem1
                    f2 (i, k + 1) = f2 (i, k + 1) - (ptem - tem) * ptem2
                endif
            endif

        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! solve tridiagonal problem for momentum
    ! -----------------------------------------------------------------------

    call tridi2 (im, km, al, ad, au, f1, f2, au, f1, f2)

    ! -----------------------------------------------------------------------
    ! recover tendencies of momentum
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, im
            udt (i, k) = (f1 (i, k) - u1 (i, k)) * rdt
            vdt (i, k) = (f2 (i, k) - v1 (i, k)) * rdt
            if (present (dusfc)) dusfc (i) = dusfc (i) + conw * del (i, k) * udt (i, k)
            if (present (dvsfc)) dvsfc (i) = dvsfc (i) + conw * del (i, k) * vdt (i, k)
            u1 (i, k) = f1 (i, k)
            v1 (i, k) = f2 (i, k)
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! pbl height for diagnostic purpose
    ! -----------------------------------------------------------------------

    do i = 1, im
        hpbl (i) = hpblx (i)
        kpbl (i) = kpblx (i)
    enddo

    return

end subroutine sa_tke_edmf_pbl

! =======================================================================
! subroutine to calcualte surface variables for PBL
! =======================================================================

subroutine sa_tke_edmf_sfc (im, lsoil, ps, u1, v1, t1, q1, &
        delt, tsurf, prsl1, prslki, evap, hflx, fm, fh, &
        z1, snwdph, zorl, ztrl, islimsk, ustar, sigmaf, &
        vegtype, shdmax, sfcemis, dlwflx, sfcnsw, &
        sfcdsw, srflag, hice, fice, tice, weasd, &
        tprcp, stc, qsurf, cmm, chh, gflux, ep, &
        u10m_out, v10m_out, t2m_out, q2m_out, &
        cm_out, ch_out, rb_out, stress_out, wind_out)

    implicit none

    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------

    integer, intent (in) :: im, lsoil

    integer, intent (in) :: islimsk (im), vegtype (im)

    real, intent (in) :: delt

    real, intent (in) :: ps (im), u1 (im), v1 (im), t1 (im), q1 (im), &
        prslki (im), z1 (im), prsl1 (im), sigmaf (im), shdmax (im), &
        sfcemis (im), dlwflx (im), sfcnsw (im), sfcdsw (im), srflag (im)

    real, intent (inout) :: fm (im), fh (im), zorl (im), ztrl (im), ustar (im), snwdph (im), &
        hice (im), fice (im), tice (im), weasd (im), tprcp (im), stc (im, lsoil), &
        evap (im), hflx (im), tsurf (im), qsurf (im), cmm (im), chh (im), &
        gflux (im), ep (im)

    real, intent (out), optional :: u10m_out (im), v10m_out (im), &
        t2m_out (im), q2m_out (im), cm_out (im), ch_out (im), rb_out (im), &
        stress_out (im), wind_out (im)

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    logical :: mom4ice = .false.

    integer :: lsm = 1

    real :: fm10 (im), fh2 (im), u10m (im), v10m (im), t2m (im), q2m (im), &
        cm (im), ch (im), rb (im), stress (im), wind (im), snowmt (im)

    ! -----------------------------------------------------------------------
    ! calculate surface exchange coefficients and near-surface wind
    ! -----------------------------------------------------------------------

    if ( sfc_gfdl ) then

        call sfc_exch_gfdl (im, ps, u1, v1, t1, q1, z1, &
            snwdph, tsurf, zorl, ztrl, cm, ch, rb, &
            prsl1, prslki, islimsk, stress, fm, fh, &
            ustar, wind, fm10, fh2, sigmaf, vegtype, shdmax)

    else

        call sfc_exch (im, ps, u1, v1, t1, q1, z1, &
            snwdph, tsurf, zorl, cm, ch, rb, &
            prsl1, prslki, islimsk, stress, fm, fh, &
            ustar, wind, fm10, fh2, sigmaf, vegtype, shdmax)

    endif

    ! -----------------------------------------------------------------------
    ! surface energy balance over ocean
    ! -----------------------------------------------------------------------

    call sfc_ocea (im, ps, u1, v1, t1, q1, tsurf, cm, ch, &
        prsl1, prslki, islimsk, qsurf, cmm, chh, gflux, evap, hflx, ep)

    ! -----------------------------------------------------------------------
    ! surface energy balance over land
    ! -----------------------------------------------------------------------

    ! TBD

    ! -----------------------------------------------------------------------
    ! surface energy balance over seaice
    ! -----------------------------------------------------------------------

    call sfc_seai (im, lsoil, ps, u1, v1, t1, q1, delt, &
        sfcemis, dlwflx, sfcnsw, sfcdsw, srflag, &
        cm, ch, prsl1, prslki, islimsk, mom4ice, lsm, &
        hice, fice, tice, weasd, tsurf, tprcp, stc, ep, &
        snwdph, qsurf, snowmt, gflux, cmm, chh, evap, hflx)

    ! -----------------------------------------------------------------------
    ! update near surface fields
    ! -----------------------------------------------------------------------

    call sfc_updt (im, ps, u1, v1, t1, q1, &
        tsurf, qsurf, u10m, v10m, t2m, q2m, &
        prslki, evap, fm, fh, fm10, fh2)

    ! -----------------------------------------------------------------------
    ! optional output
    ! -----------------------------------------------------------------------

    if (present (u10m_out)) u10m_out = u10m
    if (present (v10m_out)) v10m_out = v10m
    if (present (t2m_out)) t2m_out = t2m
    if (present (q2m_out)) q2m_out = q2m
    if (present (cm_out)) cm_out = cm
    if (present (ch_out)) ch_out = ch
    if (present (rb_out)) rb_out = rb
    if (present (stress_out)) stress_out = stress
    if (present (wind_out)) wind_out = wind

end subroutine sa_tke_edmf_sfc

! =======================================================================
! subroutine to calculate surface exchange coefficients and near-surface wind
! =======================================================================

subroutine sfc_exch (im, ps, u1, v1, t1, q1, z1, &
        snwdph, tsurf, zorl, cm, ch, rb, &
        prsl1, prslki, islimsk, stress, fm, fh, &
        ustar, wind, fm10, fh2, sigmaf, vegtype, shdmax)

    implicit none

    integer im
    real, dimension (im) :: ps, u1, v1, t1, q1, z1, &
        tsurf, zorl, cm, ch, rb, &
        prsl1, prslki, stress, &
		fm, fh, ustar, wind, ddvel, &
        fm10, fh2, sigmaf, shdmax, &
        snwdph
    integer, dimension (im) :: vegtype, islimsk

    logical flag_iter (im) ! added by s.lu

    ! locals

    integer i

    real :: aa, aa0, bb, bb0, dtv, adtv, qs1, &
        hl1, hl12, pm, ph, pm10, ph2, rat, &
        thv1, tvs, z1i, z0, z0max, ztmax, &
        fms, fhs, hl0, hl0inf, hlinf, &
        hl110, hlt, hltinf, olinf, &
        restar, tem1, tem2, ztmax1, &
        z0_adj, wind_th_moon, ustar_th, a, b, c, & !kgao
    u10m, v10m, ws10m !kgao

    real, parameter :: &
        charnock = .014, ca = .4, & ! ca - von karman constant
    alpha = 5., a0 = - 3.975, a1 = 12.32, alpha4 = 4.0 * alpha, &
        b1 = - 7.755, b2 = 6.041, alpha2 = alpha + alpha, beta = 1.0, &
        a0p = - 7.941, a1p = 24.75, b1p = - 8.705, b2p = 7.899, &
        vis = 1.4e-5, rnu = 1.51e-5, visi = 1.0 / vis, &
        log01 = log (0.01), log05 = log (0.05), log07 = log (0.07), &
        ztmin1 = - 999.0, &
        ! following is added by kgao
    bs0 = - 8.367276172397277e-12, &
        bs1 = 1.7398510865876079e-09, &
        bs2 = - 1.331896578363359e-07, &
        bs3 = 4.507055294438727e-06, &
        bs4 = - 6.508676881906914e-05, &
        bs5 = 0.00044745137674732834, &
        bs6 = - 0.0010745704660847233, &
        cf0 = 2.1151080765239772e-13, &
        cf1 = - 3.2260663894433345e-11, &
        cf2 = - 3.329705958751961e-10, &
        cf3 = 1.7648562021709124e-07, &
        cf4 = 7.107636825694182e-06, &
        cf5 = - 0.0013914681964973246, &
        cf6 = 0.0406766967657759, &
        p13 = - 1.296521881682694e-02, &
        p12 = 2.855780863283819e-01, &
        p11 = - 1.597898515251717e+00, &
        p10 = - 8.396975715683501e+00, &
        p25 = 3.790846746036765e-10, &
        p24 = 3.281964357650687e-09, &
        p23 = 1.962282433562894e-07, &
        p22 = - 1.240239171056262e-06, &
        p21 = 1.739759082358234e-07, &
        p20 = 2.147264020369413e-05, &
        p35 = 1.840430200185075e-07, &
        p34 = - 2.793849676757154e-05, &
        p33 = 1.735308193700643e-03, &
        p32 = - 6.139315534216305e-02, &
        p31 = 1.255457892775006e+00, &
        p30 = - 1.663993561652530e+01, &
        p40 = 4.579369142033410e-04

    ! parameter (charnock = .014, ca = .4) !c ca is the von karman constant
    ! parameter (alpha = 5., a0 = - 3.975, a1 = 12.32, b1 = - 7.755, b2 = 6.041)
    ! parameter (a0p = - 7.941, a1p = 24.75, b1p = - 8.705, b2p = 7.899, vis = 1.4e-5)

    ! real aa1, bb1, bb2, cc, cc1, cc2, arnu
    ! parameter (aa1 = - 1.076, bb1 = .7045, cc1 = - .05808)
    ! parameter (bb2 = - .1954, cc2 = .009999)
    ! parameter (arnu = .135 * rnu)

    ! z0s_max = .196e-2 for u10_crit = 25 m / s
    ! z0s_max = .317e-2 for u10_crit = 30 m / s
    ! z0s_max = .479e-2 for u10_crit = 35 m / s

    ! mbek -- toga - coare flux algorithm
    ! parameter (rnu = 1.51e-5, arnu = 0.11 * rnu)

    ! initialize variables. all units are supposedly m.k.s. unless specified
    ! ps is in pascals, wind is wind speed,
    ! surface roughness length is converted to m from cm

    ddvel = 0.0
    flag_iter = .true.

    do i = 1, im
        if (flag_iter (i)) then
            wind (i) = max (sqrt (u1 (i) * u1 (i) + v1 (i) * v1 (i)) &
                 + max (0.0, min (ddvel (i), 30.0)), 1.0)
            tem1 = 1.0 + zvir * max (q1 (i), 1.e-8)
            thv1 = t1 (i) * prslki (i) * tem1
            tvs = 0.5 * (tsurf (i) + tsurf (i)) * tem1
            qs1 = mqs (t1 (i))
            qs1 = max (1.0e-8, eps * qs1 / (prsl1 (i) + epsm1 * qs1))

            z0 = 0.01 * zorl (i)
            z0max = max (1.0e-6, min (z0, z1 (i)))
            z1i = 1.0 / z1 (i)

            ! compute stability dependent exchange coefficients
            ! this portion of the code is presently suppressed


            if (islimsk (i) == 0) then ! over ocean
                ustar (i) = sqrt (grav * z0 / charnock)

                ! ** test xubin's new z0

                ! ztmax = z0max

                restar = max (ustar (i) * z0max * visi, 0.000001)

                ! restar = log (restar)
                ! restar = min (restar, 5.)
                ! restar = max (restar, - 5.)
                ! rat = aa1 + (bb1 + cc1 * restar) * restar
                ! rat = rat / (1. + (bb2 + cc2 * restar) * restar))
                ! rat taken from zeng, zhao and dickinson 1997

                rat = min (7.0, 2.67 * sqrt (sqrt (restar)) - 2.57)
                ztmax = z0max * exp (- rat)

            else ! over land and sea ice
                ! ** xubin's new z0 over land and sea ice
                tem1 = 1.0 - shdmax (i)
                tem2 = tem1 * tem1
                tem1 = 1.0 - tem2

                if (ivegsrc == 1) then

                    if (vegtype (i) == 10) then
                        z0max = exp (tem2 * log01 + tem1 * log07)
                    elseif (vegtype (i) == 6) then
                        z0max = exp (tem2 * log01 + tem1 * log05)
                    elseif (vegtype (i) == 7) then
                        ! z0max = exp (tem2 * log01 + tem1 * log01)
                        z0max = 0.01
                    elseif (vegtype (i) == 16) then
                        ! z0max = exp (tem2 * log01 + tem1 * log01)
                        z0max = 0.01
                    else
                        z0max = exp (tem2 * log01 + tem1 * log (z0max))
                    endif

                elseif (ivegsrc == 2) then

                    if (vegtype (i) == 7) then
                        z0max = exp (tem2 * log01 + tem1 * log07)
                    elseif (vegtype (i) == 8) then
                        z0max = exp (tem2 * log01 + tem1 * log05)
                    elseif (vegtype (i) == 9) then
                        ! z0max = exp (tem2 * log01 + tem1 * log01)
                        z0max = 0.01
                    elseif (vegtype (i) == 11) then
                        ! z0max = exp (tem2 * log01 + tem1 * log01)
                        z0max = 0.01
                    else
                        z0max = exp (tem2 * log01 + tem1 * log (z0max))
                    endif

                endif
                z0max = max (z0max, 1.0e-6)

                ! czilc = 10.0 ** (- (0.40 / 0.07) * z0) ! fei's canopy height dependance of czil
                ! czilc = 0.8

                tem1 = 1.0 - sigmaf (i)
                ztmax = z0max * exp (- tem1 * tem1 &
                     * czilc * ca * sqrt (ustar (i) * (0.01 / 1.5e-05)))

            endif ! end of if (islimsk (i) == 0) then

            ztmax = max (ztmax, 1.0e-6)
            tem1 = z0max / z1 (i)
            if (abs (1.0 - tem1) > 1.0e-6) then
                ztmax1 = - beta * log (tem1) / (alpha2 * (1. - tem1))
            else
                ztmax1 = 99.0
            endif
            if (z0max < 0.05 .and. snwdph (i) < 10.0) ztmax1 = 99.0


            ! compute stability indices (rb and hlinf)

            dtv = thv1 - tvs
            adtv = max (abs (dtv), 0.001)
            dtv = sign (1., dtv) * adtv
            rb (i) = max (- 5000.0, (grav + grav) * dtv * z1 (i) &
                 / ((thv1 + tvs) * wind (i) * wind (i)))
            tem1 = 1.0 / z0max
            tem2 = 1.0 / ztmax
            fm (i) = log ((z0max + z1 (i)) * tem1)
            fh (i) = log ((ztmax + z1 (i)) * tem2)
            fm10 (i) = log ((z0max + 10.) * tem1)
            fh2 (i) = log ((ztmax + 2.) * tem2)
            hlinf = rb (i) * fm (i) * fm (i) / fh (i)
            hlinf = min (max (hlinf, ztmin1), ztmax1)

            ! stable case

            if (dtv >= 0.0) then
                hl1 = hlinf
                if (hlinf > .25) then
                    tem1 = hlinf * z1i
                    hl0inf = z0max * tem1
                    hltinf = ztmax * tem1
                    aa = sqrt (1. + alpha4 * hlinf)
                    aa0 = sqrt (1. + alpha4 * hl0inf)
                    bb = aa
                    bb0 = sqrt (1. + alpha4 * hltinf)
                    pm = aa0 - aa + log ((aa + 1.) / (aa0 + 1.))
                    ph = bb0 - bb + log ((bb + 1.) / (bb0 + 1.))
                    fms = fm (i) - pm
                    fhs = fh (i) - ph
                    hl1 = fms * fms * rb (i) / fhs
                    hl1 = min (max (hl1, ztmin1), ztmax1)
                endif

                ! second iteration

                tem1 = hl1 * z1i
                hl0 = z0max * tem1
                hlt = ztmax * tem1
                aa = sqrt (1. + alpha4 * hl1)
                aa0 = sqrt (1. + alpha4 * hl0)
                bb = aa
                bb0 = sqrt (1. + alpha4 * hlt)
                pm = aa0 - aa + log ((1.0 + aa) / (1.0 + aa0))
                ph = bb0 - bb + log ((1.0 + bb) / (1.0 + bb0))
                hl110 = hl1 * 10. * z1i
                hl110 = min (max (hl110, ztmin1), ztmax1)
                aa = sqrt (1. + alpha4 * hl110)
                pm10 = aa0 - aa + log ((1.0 + aa) / (1.0 + aa0))
                hl12 = (hl1 + hl1) * z1i
                hl12 = min (max (hl12, ztmin1), ztmax1)
                ! aa = sqrt (1. + alpha4 * hl12)
                bb = sqrt (1. + alpha4 * hl12)
                ph2 = bb0 - bb + log ((1.0 + bb) / (1.0 + bb0))

                ! unstable case - check for unphysical obukhov length

            else ! dtv < 0 case
                olinf = z1 (i) / hlinf
                tem1 = 50.0 * z0max
                if (abs (olinf) <= tem1) then
                    hlinf = - z1 (i) / tem1
                    hlinf = min (max (hlinf, ztmin1), ztmax1)
                endif

                ! get pm and ph

                if (hlinf >= - 0.5) then
                    hl1 = hlinf
                    pm = (a0 + a1 * hl1) * hl1 / (1. + (b1 + b2 * hl1) * hl1)
                    ph = (a0p + a1p * hl1) * hl1 / (1. + (b1p + b2p * hl1) * hl1)
                    hl110 = hl1 * 10. * z1i
                    hl110 = min (max (hl110, ztmin1), ztmax1)
                    pm10 = (a0 + a1 * hl110) * hl110 / (1. + (b1 + b2 * hl110) * hl110)
                    hl12 = (hl1 + hl1) * z1i
                    hl12 = min (max (hl12, ztmin1), ztmax1)
                    ph2 = (a0p + a1p * hl12) * hl12 / (1. + (b1p + b2p * hl12) * hl12)
                else ! hlinf < 0.05
                    hl1 = - hlinf
                    tem1 = 1.0 / sqrt (hl1)
                    pm = log (hl1) + 2. * sqrt (tem1) - .8776
                    ph = log (hl1) + .5 * tem1 + 1.386
                    ! pm = log (hl1) + 2.0 * hl1 ** (- .25) - .8776
                    ! ph = log (hl1) + 0.5 * hl1 ** (- .5) + 1.386
                    hl110 = hl1 * 10. * z1i
                    hl110 = min (max (hl110, ztmin1), ztmax1)
                    pm10 = log (hl110) + 2.0 / sqrt (sqrt (hl110)) - .8776
                    ! pm10 = log (hl110) + 2. * hl110 ** (- .25) - .8776
                    hl12 = (hl1 + hl1) * z1i
                    hl12 = min (max (hl12, ztmin1), ztmax1)
                    ph2 = log (hl12) + 0.5 / sqrt (hl12) + 1.386
                    ! ph2 = log (hl12) + .5 * hl12 ** (- .5) + 1.386
                endif

            endif ! end of if (dtv >= 0) then loop

            ! finish the exchange coefficient computation to provide fm and fh

            fm (i) = fm (i) - pm
            fh (i) = fh (i) - ph
            fm10 (i) = fm10 (i) - pm10
            fh2 (i) = fh2 (i) - ph2
            cm (i) = ca * ca / (fm (i) * fm (i))
            ch (i) = ca * ca / (fm (i) * fh (i))
            tem1 = 0.00001 / z1 (i)
            cm (i) = max (cm (i), tem1)
            ch (i) = max (ch (i), tem1)
            stress (i) = cm (i) * wind (i) * wind (i)
            ustar (i) = sqrt (stress (i))

            ! update z0 over ocean

            if (islimsk (i) == 0) then

                z0 = (charnock / grav) * ustar (i) * ustar (i)

                ! mbek -- toga - coare flux algorithm
                ! z0 = (charnock / grav) * ustar (i) * ustar (i) + arnu / ustar (i)
                ! new implementation of z0
                ! cc = ustar (i) * z0 / rnu
                ! pp = cc / (1. + cc)
                ! ff = grav * arnu / (charnock * ustar (i) ** 3)
                ! z0 = arnu / (ustar (i) * ff ** pp)

                ! -------------------------- modify z0 by kgao

                ! diagnose 10m wind (same as sfc_diag.f)

                u10m = u1 (i) * fm10 (i) / fm (i)
                v10m = v1 (i) * fm10 (i) / fm (i)
                ws10m = sqrt (u10m * u10m + v10m * v10m)

                ! option - uri / gfdl (hwrf 2015)
                ! note there is discontinuity at 10m / s in original formulation
                ! needs to be fixed

                if (do_z0_hwrf15) then
                    if (ws10m <= 5.0) then
                        z0 = 0.0185 / 9.8 * (7.59e-4 * ws10m ** 2 + 2.46e-2 * ws10m) ** 2
                    elseif (ws10m > 5.0 .and. ws10m <= 10.) then
                        z0 = 0.00000235 * (ws10m ** 2 - 25.) + 3.805129199617346e-05
                    elseif (ws10m > 10.0 .and. ws10m <= 60.) then
                        z0 = bs6 + bs5 * ws10m + bs4 * ws10m ** 2 + bs3 * ws10m ** 3 &
                             + bs2 * ws10m ** 4 + bs1 * ws10m ** 5 + bs0 * ws10m ** 6
                    else
                        z0 = cf6 + cf5 * ws10m + cf4 * ws10m ** 2 + cf3 * ws10m ** 3 &
                             + cf2 * ws10m ** 4 + cf1 * ws10m ** 5 + cf0 * ws10m ** 6
                    endif
                endif

                ! option - hwrf 2017

                if (do_z0_hwrf17) then
                    if (ws10m <= 6.5) then
                        z0 = exp (p10 + p11 * ws10m + p12 * ws10m ** 2 + p13 * ws10m ** 3)
                    elseif (ws10m > 6.5 .and. ws10m <= 15.7) then
                        z0 = p25 * ws10m ** 5 + p24 * ws10m ** 4 + p23 * ws10m ** 3 &
                             + p22 * ws10m ** 2 + p21 * ws10m + p20
                    elseif (ws10m > 15.7 .and. ws10m <= 53.) then
                        z0 = exp (p35 * ws10m ** 5 + p34 * ws10m ** 4 + p33 * ws10m ** 3 &
                             + p32 * ws10m ** 2 + p31 * ws10m + p30)
                    else
                        z0 = p40
                    endif
                endif

                ! option - gfs (low wind) + hwrf 2017 (high wind)

                if (do_z0_hwrf17_hwonly) then

                    if (ws10m > wind_th_hwrf .and. ws10m <= 53.) then
                        z0 = exp (p35 * ws10m ** 5 + p34 * ws10m ** 4 + p33 * ws10m ** 3 &
                             + p32 * ws10m ** 2 + p31 * ws10m + p30)
                    elseif (ws10m > 53.) then
                        z0 = p40
                    endif

                endif

                ! option - gfs (low wind) + moon et al (high wind)

                if (do_z0_moon) then
                    wind_th_moon = 20.
                    a = 0.56
                    b = - 20.255
                    c = wind_th_moon - 2.458
                    ustar_th = (- b - sqrt (b * b - 4 * a * c)) / (2 * a)

                    z0_adj = 0.001 * (0.085 * wind_th_moon - 0.58) - &
                         (charnock / grav) * ustar_th * ustar_th

                    ws10m = 2.458 + ustar (i) * (20.255 - 0.56 * ustar (i)) ! eq (7) moon et al. 2007
                    if (ws10m > wind_th_moon) then ! no modification in low wind conditions
                        z0 = 0.001 * (0.085 * ws10m - 0.58) - z0_adj ! eq (8b) moon et al. 2007 modified by kgao
                    endif
                endif

                ! ---------------------------- modify z0 end

                if (redrag) then
                    zorl (i) = 100.0 * max (min (z0, z0s_max), 1.e-7)
                else
                    zorl (i) = 100.0 * max (min (z0, .1), 1.e-7)
                endif
            endif
        endif ! end of if (flagiter) loop
    enddo

end subroutine sfc_exch

! =======================================================================
! subroutine to calculate surface exchange coefficients and near-surface wind
! Oct 2019 - a clean and updated version by Kun Gao at GFDL (kun.gao@noaa.gov)
! =======================================================================

subroutine sfc_exch_gfdl (im, ps, u1, v1, t1, q1, z1, &
        snwdph, tsurf, zorl, ztrl, cm, ch, rb, &
        prsl1, prslki, islimsk, &
        stress, fm, fh, &
        ustar, wind, fm10, fh2, &
        sigmaf, vegtype, shdmax)

    implicit none

    ! --- input / output

    integer im

    real, dimension (im) :: ps, u1, v1, t1, q1, z1, &
        zorl, ztrl, cm, ch, rb, &
        prsl1, prslki, stress, &
        fm, fh, ustar, wind, ddvel, &
        fm10, fh2, sigmaf, shdmax, &
        tsurf, snwdph

    integer, dimension (im) :: vegtype, islimsk

    logical flag_iter (im)

    ! --- local

    integer i

    real :: aa, aa0, bb, bb0, dtv, adtv, qs1, &
        hl1, hl12, pm, ph, pm10, ph2, rat, &
        thv1, tvs, z1i, z0, zt, z0max, ztmax, &
        fms, fhs, hl0, hl0inf, hlinf, &
        hl110, hlt, hltinf, olinf, &
        restar, czilc, tem1, tem2, &
        u10m, v10m, ws10m, ws10m_moon, &
        z0_1, zt_1, fm1, fh1, ustar_1, ztmax_1

    real, parameter :: &
        charnock = .014, ca = .4, &
        vis = 1.4e-5, rnu = 1.51e-5, visi = 1.0 / vis, &
        log01 = log (0.01), log05 = log (0.05), log07 = log (0.07), &
        ztmin1 = - 999.0

    ! ================================================
    ! main program starts here
    ! ================================================

    ddvel = 0.0
    flag_iter = .true.

    do i = 1, im

        if (flag_iter (i)) then

            ! --- get variables at model lowest layer and surface (water / ice / land)

            wind (i) = max (sqrt (u1 (i) * u1 (i) + v1 (i) * v1 (i)) &
                + max (0.0, min (ddvel (i), 30.0)), 1.0)
            tem1 = 1.0 + zvir * max (q1 (i), 1.e-8)
            thv1 = t1 (i) * prslki (i) * tem1
            tvs = 0.5 * (tsurf (i) + tsurf (i)) * tem1
            qs1 = mqs (t1 (i))
            qs1 = max (1.0e-8, eps * qs1 / (prsl1 (i) + epsm1 * qs1))

            ! (sea / land / ice mask = 0 / 1 / 2)
            if (islimsk (i) == 1 .or. islimsk (i) == 2) then ! over land or sea ice

                ! ================================================
                ! if over land or sea ice:
                ! step 1 - get z0 / zt
                ! step 2 - call similarity
                ! ================================================

                ! --- get surface roughness for momentum (z0)

                z0 = 0.01 * zorl (i)
                z0max = max (1.0e-6, min (z0, z1 (i)))

                !xubin's new z0 over land and sea ice
                tem1 = 1.0 - shdmax (i) ! shdmax is max vegetation area fraction
                tem2 = tem1 * tem1
                tem1 = 1.0 - tem2

                if (ivegsrc == 1) then

                    if (vegtype (i) == 10) then
                        z0max = exp (tem2 * log01 + tem1 * log07)
                    elseif (vegtype (i) == 6) then
                        z0max = exp (tem2 * log01 + tem1 * log05)
                    elseif (vegtype (i) == 7) then
                        ! z0max = exp (tem2 * log01 + tem1 * log01)
                        z0max = 0.01
                    elseif (vegtype (i) == 16) then
                        ! z0max = exp (tem2 * log01 + tem1 * log01)
                        z0max = 0.01
                    else
                        z0max = exp (tem2 * log01 + tem1 * log (z0max))
                    endif

                elseif (ivegsrc == 2) then

                    if (vegtype (i) == 7) then
                        z0max = exp (tem2 * log01 + tem1 * log07)
                    elseif (vegtype (i) == 8) then
                        z0max = exp (tem2 * log01 + tem1 * log05)
                    elseif (vegtype (i) == 9) then
                        ! z0max = exp (tem2 * log01 + tem1 * log01)
                        z0max = 0.01
                    elseif (vegtype (i) == 11) then
                        ! z0max = exp (tem2 * log01 + tem1 * log01)
                        z0max = 0.01
                    else
                        z0max = exp (tem2 * log01 + tem1 * log (z0max))
                    endif

                    z0max = max (z0max, 1.0e-6)

                endif

                ! --- get surface roughness for heat (zt)

                ! czilc = 10.0 ** (- (0.40 / 0.07) * z0) ! let czilc depend on canopy height
                czilc = 0.8

                tem1 = 1.0 - sigmaf (i)
                ztmax = z0max * exp (- tem1 * tem1 * &
                    czilc * ca * sqrt (ustar (i) * (0.01 / 1.5e-05)))

                ztmax = max (ztmax, 1.0e-6)

                ! --- call similarity

                call monin_obukhov_similarity &
                    (z1 (i), snwdph (i), thv1, wind (i), z0max, ztmax, tvs, &
                    rb (i), fm (i), fh (i), fm10 (i), fh2 (i), &
                    cm (i), ch (i), stress (i), ustar (i))

            elseif (islimsk (i) == 0) then ! over water

                ! ================================================
                ! if over water (redesigned by kun gao)
                ! iteration 1
                ! step 1 get z0 / zt from previous step
                ! step 2 call similarity
                ! iteration 2
                ! step 1 update z0 / zt
                ! step 2 call similarity
                ! ================================================

                ! === iteration 1

                ! --- get z0 / zt
                z0 = 0.01 * zorl (i)
                zt = 0.01 * ztrl (i)

                z0max = max (1.0e-6, min (z0, z1 (i)))
                ztmax = max (zt, 1.0e-6)

                ! --- call similarity
                call monin_obukhov_similarity &
                    (z1 (i), snwdph (i), thv1, wind (i), z0max, ztmax, tvs, &
                    rb (i), fm (i), fh (i), fm10 (i), fh2 (i), &
                    cm (i), ch (i), stress (i), ustar (i))

                ! === iteration 2

                ! --- get z0 / zt following the old sfc_diff.f
                z0 = (charnock / grav) * ustar (i) * ustar (i)
                if (redrag) then
                    z0 = max (min (z0, z0s_max), 1.e-7)
                else
                    z0 = max (min (z0, .1), 1.e-7)
                endif

                ! zt calculations copied from old sfc_diff.f
                !ustar (i) = sqrt (grav * z0 / charnock)
                !restar = max (ustar (i) * z0max * visi, 0.000001)
                !rat = min (7.0, 2.67 * sqrt (sqrt (restar)) - 2.57)
                !ztmax = z0max * exp (- rat)

                ustar_1 = sqrt (grav * z0 / charnock)
                restar = max (ustar_1 * z0max * visi, 0.000001)
                rat = min (7.0, 2.67 * sqrt (sqrt (restar)) - 2.57)
                zt = z0max * exp (- rat) ! zeng, zhao and dickinson 1997 (eq 25)

                ! --- update z0 / zt with new options
                ! only z0 options in the following
                ! will add zt options in the future

                u10m = u1 (i) * fm10 (i) / fm (i)
                v10m = v1 (i) * fm10 (i) / fm (i)
                ws10m = sqrt (u10m * u10m + v10m * v10m)

                if (do_z0_hwrf15) then
                    ! option 1: hwrf15, originally developed by uri / gfdl
                    call cal_z0_hwrf15 (ws10m, z0)
                    call cal_zt_hwrf15 (ws10m, zt)

                elseif (do_z0_hwrf17) then
                    ! option 2: hwrf17
                    call cal_z0_hwrf17 (ws10m, z0)
                    call cal_zt_hwrf17 (ws10m, zt)

                elseif (do_z0_hwrf17_hwonly) then
                    ! option 3: hwrf17 under high wind only
                    if (ws10m > wind_th_hwrf) then
                        call cal_z0_hwrf17 (ws10m, z0)
                        z0 = max (min (z0, z0s_max), 1.e-7) ! must apply limiter here
                    endif

                elseif (do_z0_moon) then
                    ! option 4: moon et al 2007 under high winds (same as in hiram)
                    ws10m_moon = 2.458 + ustar (i) * (20.255 - 0.56 * ustar (i)) ! eq (7) moon et al. 2007
                    if (ws10m_moon > 20.) then
                        call cal_z0_moon (ws10m_moon, z0)
                        z0 = max (min (z0, z0s_max), 1.e-7) ! must apply limiter here
                    endif
                endif

                z0max = max (z0, 1.0e-6)
                ztmax = max (zt, 1.0e-6)

                ! --- call similarity
                call monin_obukhov_similarity &
                    (z1 (i), snwdph (i), thv1, wind (i), z0max, ztmax, tvs, &
                    rb (i), fm (i), fh (i), fm10 (i), fh2 (i), &
                    cm (i), ch (i), stress (i), ustar (i))

                zorl (i) = 100.0 * z0max
                ztrl (i) = 100.0 * ztmax

            endif ! end of if (islimsk) loop
        endif ! end of if (flagiter) loop
    enddo ! end of do i = 1, im loop

    return

end subroutine sfc_exch_gfdl

! =======================================================================
! Originally developed by URI/GFDL
! Coded by Kun Gao (kun.gao@noaa.gov)
! =======================================================================

subroutine cal_z0_hwrf15 (ws10m, z0)

    real :: ws10m, z0

    real, parameter :: &
        a0 = - 8.367276172397277e-12, &
        a1 = 1.7398510865876079e-09, &
        a2 = - 1.331896578363359e-07, &
        a3 = 4.507055294438727e-06, &
        a4 = - 6.508676881906914e-05, &
        a5 = 0.00044745137674732834, &
        a6 = - 0.0010745704660847233, &
        b0 = 2.1151080765239772e-13, &
        b1 = - 3.2260663894433345e-11, &
        b2 = - 3.329705958751961e-10, &
        b3 = 1.7648562021709124e-07, &
        b4 = 7.107636825694182e-06, &
        b5 = - 0.0013914681964973246, &
        b6 = 0.0406766967657759

    if (ws10m <= 5.0) then
        z0 = 0.0185 / 9.8 * (7.59e-4 * ws10m ** 2 + 2.46e-2 * ws10m) ** 2
    elseif (ws10m > 5.0 .and. ws10m <= 10.) then
        z0 = 0.00000235 * (ws10m ** 2 - 25.) + 3.805129199617346e-05
    elseif (ws10m > 10.0 .and. ws10m <= 60.) then
        z0 = a6 + a5 * ws10m + a4 * ws10m ** 2 + a3 * ws10m ** 3 + &
            a2 * ws10m ** 4 + a1 * ws10m ** 5 + a0 * ws10m ** 6
    else
        z0 = b6 + b5 * ws10m + b4 * ws10m ** 2 + b3 * ws10m ** 3 + &
            b2 * ws10m ** 4 + b1 * ws10m ** 5 + b0 * ws10m ** 6
    endif

end subroutine cal_z0_hwrf15

! =======================================================================
! Originally developed by URI/GFDL
! Coded by Kun Gao (kun.gao@noaa.gov)
! =======================================================================

subroutine cal_zt_hwrf15 (ws10m, zt)

    real :: ws10m, zt

    real, parameter :: &
        a0 = 2.51715926619e-09, &
        a1 = - 1.66917514012e-07, &
        a2 = 4.57345863551e-06, &
        a3 = - 6.64883696932e-05, &
        a4 = 0.00054390175125, &
        a5 = - 0.00239645231325, &
        a6 = 0.00453024927761, &
        b0 = - 1.72935914649e-14, &
        b1 = 2.50587455802e-12, &
        b2 = - 7.90109676541e-11, &
        b3 = - 4.40976353607e-09, &
        b4 = 3.68968179733e-07, &
        b5 = - 9.43728336756e-06, &
        b6 = 8.90731312383e-05, &
        c0 = 4.68042680888e-14, &
        c1 = - 1.98125754931e-11, &
        c2 = 3.41357133496e-09, &
        c3 = - 3.05130605309e-07, &
        c4 = 1.48243563819e-05, &
        c5 = - 0.000367207751936, &
        c6 = 0.00357204479347

    if (ws10m <= 7.0) then
        zt = 0.0185 / 9.8 * (7.59e-4 * ws10m ** 2 + 2.46e-2 * ws10m) ** 2
    elseif (ws10m > 7.0 .and. ws10m <= 15.) then
        zt = a6 + a5 * ws10m + a4 * ws10m ** 2 + a3 * ws10m ** 3 + &
            a2 * ws10m ** 4 + a1 * ws10m ** 5 + a0 * ws10m ** 6
    elseif (ws10m > 15.0 .and. ws10m <= 60.) then
        zt = b6 + b5 * ws10m + b4 * ws10m ** 2 + b3 * ws10m ** 3 + &
            b2 * ws10m ** 4 + b1 * ws10m ** 5 + b0 * ws10m ** 6
    else
        zt = c6 + c5 * ws10m + c4 * ws10m ** 2 + c3 * ws10m ** 3 + &
            c2 * ws10m ** 4 + c1 * ws10m ** 5 + c0 * ws10m ** 6
    endif

end subroutine cal_zt_hwrf15

! =======================================================================
! Coded by Kun Gao (kun.gao@noaa.gov)
! =======================================================================

subroutine cal_z0_hwrf17 (ws10m, z0)

    real :: ws10m, z0

    real, parameter :: &
        p13 = - 1.296521881682694e-02, &
        p12 = 2.855780863283819e-01, &
        p11 = - 1.597898515251717e+00, &
        p10 = - 8.396975715683501e+00, &
        p25 = 3.790846746036765e-10, &
        p24 = 3.281964357650687e-09, &
        p23 = 1.962282433562894e-07, &
        p22 = - 1.240239171056262e-06, &
        p21 = 1.739759082358234e-07, &
        p20 = 2.147264020369413e-05, &
        p35 = 1.840430200185075e-07, &
        p34 = - 2.793849676757154e-05, &
        p33 = 1.735308193700643e-03, &
        p32 = - 6.139315534216305e-02, &
        p31 = 1.255457892775006e+00, &
        p30 = - 1.663993561652530e+01, &
        p40 = 4.579369142033410e-04

    if (ws10m <= 6.5) then
        z0 = exp (p10 + p11 * ws10m + p12 * ws10m ** 2 + p13 * ws10m ** 3)
    elseif (ws10m > 6.5 .and. ws10m <= 15.7) then
        z0 = p25 * ws10m ** 5 + p24 * ws10m ** 4 + p23 * ws10m ** 3 + &
            p22 * ws10m ** 2 + p21 * ws10m + p20
    elseif (ws10m > 15.7 .and. ws10m <= 53.) then
        z0 = exp (p35 * ws10m ** 5 + p34 * ws10m ** 4 + p33 * ws10m ** 3 + &
            p32 * ws10m ** 2 + p31 * ws10m + p30)
    else
        z0 = p40
    endif

end subroutine cal_z0_hwrf17

! =======================================================================
! Coded by Kun Gao (kun.gao@noaa.gov)
! =======================================================================

subroutine cal_zt_hwrf17 (ws10m, zt)

    real :: ws10m, zt

    real, parameter :: p00 = 1.100000000000000e-04, &
        p15 = - 9.144581627678278e-10, p14 = 7.020346616456421e-08, &
        p13 = - 2.155602086883837e-06, p12 = 3.333848806567684e-05, &
        p11 = - 2.628501274963990e-04, p10 = 8.634221567969181e-04, &
        p25 = - 8.654513012535990e-12, p24 = 1.232380050058077e-09, &
        p23 = - 6.837922749505057e-08, p22 = 1.871407733439947e-06, &
        p21 = - 2.552246987137160e-05, p20 = 1.428968311457630e-04, &
        p35 = 3.207515102100162e-12, p34 = - 2.945761895342535e-10, &
        p33 = 8.788972147364181e-09, p32 = - 3.814457439412957e-08, &
        p31 = - 2.448983648874671e-06, p30 = 3.436721779020359e-05, &
        p45 = - 3.530687797132211e-11, p44 = 3.939867958963747e-09, &
        p43 = - 1.227668406985956e-08, p42 = - 1.367469811838390e-05, &
        p41 = 5.988240863928883e-04, p40 = - 7.746288511324971e-03, &
        p56 = - 1.187982453329086e-13, p55 = 4.801984186231693e-11, &
        p54 = - 8.049200462388188e-09, p53 = 7.169872601310186e-07, &
        p52 = - 3.581694433758150e-05, p51 = 9.503919224192534e-04, &
        p50 = - 1.036679430885215e-02, &
        p60 = 4.751256171799112e-05

    if (ws10m >= 0.0 .and. ws10m < 5.9) then
        zt = p00
    elseif (ws10m >= 5.9 .and. ws10m <= 15.4) then
        zt = p10 + ws10m * (p11 + ws10m * (p12 + ws10m * (p13 + &
            ws10m * (p14 + ws10m * p15))))
    elseif (ws10m > 15.4 .and. ws10m <= 21.6) then
        zt = p20 + ws10m * (p21 + ws10m * (p22 + ws10m * (p23 + &
            ws10m * (p24 + ws10m * p25))))
    elseif (ws10m > 21.6 .and. ws10m <= 42.2) then
        zt = p30 + ws10m * (p31 + ws10m * (p32 + ws10m * (p33 + &
            ws10m * (p34 + ws10m * p35))))
    elseif (ws10m > 42.2 .and. ws10m <= 53.3) then
        zt = p40 + ws10m * (p41 + ws10m * (p42 + ws10m * (p43 + &
            ws10m * (p44 + ws10m * p45))))
    elseif (ws10m > 53.3 .and. ws10m <= 80.0) then
        zt = p50 + ws10m * (p51 + ws10m * (p52 + ws10m * (p53 + &
            ws10m * (p54 + ws10m * (p55 + ws10m * p56)))))
    elseif (ws10m > 80.0) then
        zt = p60
    endif

end subroutine cal_zt_hwrf17

! =======================================================================
! Coded by Kun Gao (kun.gao@noaa.gov)
! =======================================================================

subroutine cal_z0_moon (ws10m, z0)

    real :: ws10m, z0

    real :: ustar_th, z0_adj

    real, parameter :: &
        charnock = .014, &
        wind_th_moon = 20., &
        a = 0.56, &
        b = - 20.255, &
        c = wind_th_moon - 2.458

    ustar_th = (- b - sqrt (b * b - 4 * a * c)) / (2 * a)

    z0_adj = 0.001 * (0.085 * wind_th_moon - 0.58) - &
        (charnock / grav) * ustar_th * ustar_th

    z0 = 0.001 * (0.085 * ws10m - 0.58) - z0_adj ! eq (8b) moon et al. 2007 modified by kgao

end subroutine cal_z0_moon

! =======================================================================
! Monin Obukhov Similarity
! =======================================================================

subroutine monin_obukhov_similarity ( &
        z1, snwdph, thv1, wind, z0max, ztmax, tvs, &
        rb, fm, fh, fm10, fh2, cm, ch, stress, ustar)

    ! --- input
    ! z1 - lowest model level height
    ! snwdph - surface snow thickness
    ! wind - wind speed at lowest model layer
    ! thv1 - virtual potential temp at lowest model layer
    ! tvs - surface temp
    ! z0max - surface roughness length for momentum
    ! ztmax - surface roughness length for heat
    !
    ! --- output
    ! rb - a bulk richardson number
    ! fm, fh - similarity function defined at lowest model layer
    ! fm10, fh2 - similarity function defined at 10m (for momentum) and 2m (for heat)
    ! cm, ch - surface exchange coefficients for momentum and heat
    ! stress - surface wind stress
    ! ustar - surface frictional velocity

    ! --- inputs:
    real, intent (in) :: z1, snwdph, thv1, wind, z0max, ztmax, tvs

    ! --- outputs:
    real, intent (out) :: rb, fm, fh, fm10, fh2, cm, ch, stress, ustar

    ! --- locals:

    real, parameter :: alpha = 5., a0 = - 3.975, &
        a1 = 12.32, alpha4 = 4.0 * alpha, &
        b1 = - 7.755, b2 = 6.041, alpha2 = alpha + alpha, beta = 1.0, &
        a0p = - 7.941, a1p = 24.75, b1p = - 8.705, b2p = 7.899, &
        ztmin1 = - 999.0, ca = .4

    real :: aa, aa0, bb, bb0, dtv, adtv, &
        hl1, hl12, pm, ph, pm10, ph2, &
        z1i, &
        fms, fhs, hl0, hl0inf, hlinf, &
        hl110, hlt, hltinf, olinf, &
        tem1, tem2, ztmax1

    z1i = 1.0 / z1

    tem1 = z0max / z1
    if (abs (1.0 - tem1) > 1.0e-6) then
        ztmax1 = - beta * log (tem1) / (alpha2 * (1. - tem1))
    else
        ztmax1 = 99.0
    endif
    if (z0max < 0.05 .and. snwdph < 10.0) ztmax1 = 99.0

    !
    ! compute stability indices (rb and hlinf)
    !
    dtv = thv1 - tvs
    adtv = max (abs (dtv), 0.001)
    dtv = sign (1., dtv) * adtv
    rb = max (- 5000.0, (grav + grav) * dtv * z1 / &
        ((thv1 + tvs) * wind * wind))
    tem1 = 1.0 / z0max
    tem2 = 1.0 / ztmax
    fm = log ((z0max + z1) * tem1)
    fh = log ((ztmax + z1) * tem2)
    fm10 = log ((z0max + 10.) * tem1)
    fh2 = log ((ztmax + 2.) * tem2)
    hlinf = rb * fm * fm / fh
    hlinf = min (max (hlinf, ztmin1), ztmax1)
    !
    ! stable case
    !
    if (dtv >= 0.0) then
        hl1 = hlinf
        if (hlinf > .25) then
            tem1 = hlinf * z1i
            hl0inf = z0max * tem1
            hltinf = ztmax * tem1
            aa = sqrt (1. + alpha4 * hlinf)
            aa0 = sqrt (1. + alpha4 * hl0inf)
            bb = aa
            bb0 = sqrt (1. + alpha4 * hltinf)
            pm = aa0 - aa + log ((aa + 1.) / (aa0 + 1.))
            ph = bb0 - bb + log ((bb + 1.) / (bb0 + 1.))
            fms = fm - pm
            fhs = fh - ph
            hl1 = fms * fms * rb / fhs
            hl1 = min (max (hl1, ztmin1), ztmax1)
        endif
        !
        ! second iteration
        !
        tem1 = hl1 * z1i
        hl0 = z0max * tem1
        hlt = ztmax * tem1
        aa = sqrt (1. + alpha4 * hl1)
        aa0 = sqrt (1. + alpha4 * hl0)
        bb = aa
        bb0 = sqrt (1. + alpha4 * hlt)
        pm = aa0 - aa + log ((1.0 + aa) / (1.0 + aa0))
        ph = bb0 - bb + log ((1.0 + bb) / (1.0 + bb0))
        hl110 = hl1 * 10. * z1i
        hl110 = min (max (hl110, ztmin1), ztmax1)
        aa = sqrt (1. + alpha4 * hl110)
        pm10 = aa0 - aa + log ((1.0 + aa) / (1.0 + aa0))
        hl12 = (hl1 + hl1) * z1i
        hl12 = min (max (hl12, ztmin1), ztmax1)
        ! aa = sqrt (1. + alpha4 * hl12)
        bb = sqrt (1. + alpha4 * hl12)
        ph2 = bb0 - bb + log ((1.0 + bb) / (1.0 + bb0))
        !
        ! unstable case - check for unphysical obukhov length
        !
    else ! dtv < 0 case
        olinf = z1 / hlinf
        tem1 = 50.0 * z0max
        if (abs (olinf) <= tem1) then
            hlinf = - z1 / tem1
            hlinf = min (max (hlinf, ztmin1), ztmax1)
        endif
        !
        ! get pm and ph
        !
        if (hlinf >= - 0.5) then
            hl1 = hlinf
            pm = (a0 + a1 * hl1) * hl1 / (1. + (b1 + b2 * hl1) * hl1)
            ph = (a0p + a1p * hl1) * hl1 / (1. + (b1p + b2p * hl1) * hl1)
            hl110 = hl1 * 10. * z1i
            hl110 = min (max (hl110, ztmin1), ztmax1)
            pm10 = (a0 + a1 * hl110) * hl110 / (1. + (b1 + b2 * hl110) * hl110)
            hl12 = (hl1 + hl1) * z1i
            hl12 = min (max (hl12, ztmin1), ztmax1)
            ph2 = (a0p + a1p * hl12) * hl12 / (1. + (b1p + b2p * hl12) * hl12)
        else ! hlinf < 0.05
            hl1 = - hlinf
            tem1 = 1.0 / sqrt (hl1)
            pm = log (hl1) + 2. * sqrt (tem1) - .8776
            ph = log (hl1) + .5 * tem1 + 1.386
            ! pm = log (hl1) + 2.0 * hl1 ** (- .25) - .8776
            ! ph = log (hl1) + 0.5 * hl1 ** (- .5) + 1.386
            hl110 = hl1 * 10. * z1i
            hl110 = min (max (hl110, ztmin1), ztmax1)
            pm10 = log (hl110) + 2.0 / sqrt (sqrt (hl110)) - .8776
            ! pm10 = log (hl110) + 2. * hl110 ** (- .25) - .8776
            hl12 = (hl1 + hl1) * z1i
            hl12 = min (max (hl12, ztmin1), ztmax1)
            ph2 = log (hl12) + 0.5 / sqrt (hl12) + 1.386
            ! ph2 = log (hl12) + .5 * hl12 ** (- .5) + 1.386
        endif

    endif ! end of if (dtv >= 0) then loop
    !
    ! finish the exchange coefficient computation to provide fm and fh
    !
    fm = fm - pm
    fh = fh - ph
    fm10 = fm10 - pm10
    fh2 = fh2 - ph2
    cm = ca * ca / (fm * fm)
    ch = ca * ca / (fm * fh)
    tem1 = 0.00001 / z1
    cm = max (cm, tem1)
    ch = max (ch, tem1)
    stress = cm * wind * wind
    ustar = sqrt (stress)

    return

end subroutine monin_obukhov_similarity

! =======================================================================
! subroutine to surface energy balance over ocean
!
! program history log:
! 2005 -- created from the original progtm to account for ocean only
! oct 2006 -- h. wei added cmm and chh to the output
! apr 2009 -- y. - t. hou modified to match the modified gbphys.f
! reformatted the code and added program documentation
! sep 2009 -- s. moorthi removed rcl and made pa as pressure unit
! and furthur reformatted the code
!
! inputs: size
! im - integer, horizontal dimension 1
! ps - real, surface pressure im
! u1, v1 - real, u / v component of surface layer wind im
! t1 - real, surface layer mean temperature (k) im
! q1 - real, surface layer mean specific humidity im
! tsurf - real, ground surface skin temperature (k) im
! cm - real, surface exchange coeff for momentum (m / s) im
! ch - real, surface exchange coeff heat & moisture (m / s) im
! prsl1 - real, surface layer mean pressure im
! prslki - real, im
! islimsk - integer, sea / land / ice mask (= 0 / 1 / 2) im
! ddvel - real, wind enhancement due to convection (m / s) im
! flag_iter - logical, im
!
! outputs: size
! qsurf - real, specific humidity at sfc im
! cmm - real, im
! chh - real, im
! gflux - real, ground heat flux (zero for ocean) im
! evap - real, evaporation from latent heat flux im
! hflx - real, sensible heat flux im
! ep - real, potential evaporation im
! =======================================================================

subroutine sfc_ocea (im, ps, u1, v1, t1, q1, tsurf, cm, ch, &
    prsl1, prslki, islimsk, qsurf, cmm, chh, gflux, evap, hflx, ep)

    implicit none

    ! --- constant parameters:
    real, parameter :: cpinv = 1.0 / cp_air, &
        hvapi = 1.0 / hlv, &
        elocp = hlv / cp_air

    ! --- inputs:
    integer, intent (in) :: im

    real, dimension (im), intent (in) :: ps, u1, v1, &
        t1, q1, tsurf, cm, ch, prsl1, prslki
    integer, dimension (im), intent (in) :: islimsk

    ! --- outputs:
    real, dimension (im), intent (inout) :: qsurf, &
        cmm, chh, gflux, evap, hflx, ep

    ! --- locals:

    real :: q0, qss, rch, rho, wind, tem, ddvel (im)

    integer :: i

    logical :: flag (im), flag_iter (im)
    !
    ! ===> ... begin here
    !
    ! --- ... flag for open water

    ddvel = 0.0
    flag_iter = .true.

    do i = 1, im
        flag (i) = (islimsk (i) == 0 .and. flag_iter (i))

        ! --- ... initialize variables. all units are supposedly m.k.s. unless specified
        ! ps is in pascals, wind is wind speed,
        ! rho is density, qss is sat. hum. at surface

        if (flag (i)) then

            wind = max (sqrt (u1 (i) * u1 (i) + v1 (i) * v1 (i)) &
                 + max (0.0, min (ddvel (i), 30.0)), 1.0)

            q0 = max (q1 (i), 1.0e-8)
            rho = prsl1 (i) / (rdgas * t1 (i) * (1.0 + zvir * q0))

            qss = mqs (tsurf (i))
            qss = eps * qss / (ps (i) + epsm1 * qss)

            evap (i) = 0.0
            hflx (i) = 0.0
            ep (i) = 0.0
            gflux (i) = 0.0

            ! --- ... rcp = rho cp_air ch v

            rch = rho * cp_air * ch (i) * wind
            cmm (i) = cm (i) * wind
            chh (i) = rho * ch (i) * wind

            ! --- ... sensible and latent heat flux over open water

            hflx (i) = rch * (tsurf (i) - t1 (i) * prslki (i))

            evap (i) = elocp * rch * (qss - q0)
            qsurf (i) = qss

            tem = 1.0 / rho
            hflx (i) = hflx (i) * tem * cpinv
            evap (i) = evap (i) * tem * hvapi
        endif
    enddo

end subroutine sfc_ocea

! =======================================================================
! subroutine to surface energy balance over land
! =======================================================================

! =======================================================================
! subroutine to surface energy balance over seaice
!
! program history log:
! 2005 -- xingren wu created from original progtm and added
! two - layer ice model
! 200x -- sarah lu added flag_iter
! oct 2006 -- h. wei added cmm and chh to output
! 2007 -- x. wu modified for mom4 coupling (i.e. mom4ice)
! 2007 -- s. moorthi micellaneous changes
! may 2009 -- y. - t. hou modified to include surface emissivity
! effect on lw radiation. replaced the confusing
! slrad with sfc net sw sfcnsw (dn - up) . reformatted
! the code and add program documentation block.
! sep 2009 -- s. moorthi removed rcl, changed pressure units and
! further optimized
! jan 2015 -- x. wu change "cimin = 0.15" for both
! uncoupled and coupled case
!
! inputs: size
! im, km - integer, horiz dimension and num of soil layers 1
! ps - real, surface pressure im
! u1, v1 - real, u / v component of surface layer wind im
! t1 - real, surface layer mean temperature (k) im
! q1 - real, surface layer mean specific humidity im
! delt - real, time interval (second) 1
! sfcemis - real, sfc lw emissivity (fraction) im
! dlwflx - real, total sky sfc downward lw flux (w / m ** 2) im
! sfcnsw - real, total sky sfc netsw flx into ground (w / m ** 2) im
! sfcdsw - real, total sky sfc downward sw flux (w / m ** 2) im
! srflag - real, snow / rain flag for precipitation im
! cm - real, surface exchange coeff for momentum (m / s) im
! ch - real, surface exchange coeff heat & moisture (m / s) im
! prsl1 - real, surface layer mean pressure im
! prslki - real, im
! islimsk - integer, sea / land / ice mask (= 0 / 1 / 2) im
! ddvel - real, im
! flag_iter - logical, im
! mom4ice - logical, im
! islimsk - integer, flag for land surface model scheme 1
! = 0: use osu scheme; = 1: use noah scheme
!
! input / outputs:
! hice - real, sea - ice thickness im
! fice - real, sea - ice concentration im
! tice - real, sea - ice surface temperature im
! weasd - real, water equivalent accumulated snow depth (mm) im
! tsurf - real, ground surface skin temperature (k) im
! tprcp - real, total precipitation im
! stc - real, soil temp (k) im, km
! ep - real, potential evaporation im
!
! outputs:
! snwdph - real, water equivalent snow depth (mm) im
! qsurf - real, specific humidity at sfc im
! snowmt - real, snow melt (m) im
! gflux - real, soil heat flux (w / m ** 2) im
! cmm - real, im
! chh - real, im
! evap - real, evaperation from latent heat flux im
! hflx - real, sensible heat flux im
! =======================================================================

subroutine sfc_seai (im, km, ps, u1, v1, t1, q1, delt, &
    sfcemis, dlwflx, sfcnsw, sfcdsw, srflag, &
    cm, ch, prsl1, prslki, islimsk, mom4ice, lsm, &
    hice, fice, tice, weasd, tsurf, tprcp, stc, ep, &
    snwdph, qsurf, snowmt, gflux, cmm, chh, evap, hflx)

    implicit none

    ! --- constant parameters:
    integer, parameter :: kmi = 2 ! 2 - layer of ice
    real, parameter :: cpinv = 1.0 / cp_air
    real, parameter :: hvapi = 1.0 / hlv
    real, parameter :: elocp = hlv / cp_air
    real, parameter :: himax = 8.0 ! maximum ice thickness allowed
    real, parameter :: himin = 0.1 ! minimum ice thickness required
    real, parameter :: hsmax = 2.0 ! maximum snow depth allowed
    real, parameter :: timin = 173.0 ! minimum temperature allowed for snow / ice
    real, parameter :: albfw = 0.06 ! albedo for lead
    real, parameter :: dsi = 1.0 / 0.33

    ! --- inputs:
    integer, intent (in) :: im, km, lsm

    real, dimension (im), intent (in) :: ps, u1, v1, &
        t1, q1, sfcemis, dlwflx, sfcnsw, sfcdsw, srflag, cm, ch, &
        prsl1, prslki

    integer, dimension (im), intent (in) :: islimsk
    real, intent (in) :: delt

    logical, intent (in) :: mom4ice

    ! --- input / outputs:
    real, dimension (im), intent (inout) :: hice, &
        fice, tice, weasd, tsurf, tprcp, ep

    real, dimension (im, km), intent (inout) :: stc

    ! --- outputs:
    real, dimension (im), intent (inout) :: snwdph, &
        qsurf, snowmt, gflux, cmm, chh, evap, hflx

    ! --- locals:
    real, dimension (im) :: ffw, evapi, evapw, &
        sneti, snetw, hfd, hfi, &
        ! hflxi, hflxw, sneti, snetw, qssi, qssw, hfd, hfi, hfw, &
    focn, snof, hi_save, hs_save, rch, rho, &
        snowd, theta1, ddvel

    real :: t12, t14, tem, stsice (im, kmi), &
        hflxi, hflxw, q0, qs1, wind, qssi, qssw
    real, parameter :: cimin = 0.15 ! --- minimum ice concentration

    integer :: i, k

    logical :: flag (im), flag_iter (im)
    !
    ! ===> ... begin here
    !
    ! --- ... set flag for sea - ice

    ddvel = 0.0
    flag_iter = .true.

    do i = 1, im
        flag (i) = (islimsk (i) >= 2) .and. flag_iter (i)
        if (flag_iter (i) .and. islimsk (i) < 2) then
            hice (i) = 0.0
            fice (i) = 0.0
        endif
    enddo

    ! --- ... update sea ice temperature

    do k = 1, kmi
        do i = 1, im
            if (flag (i)) then
                stsice (i, k) = stc (i, k)
            endif
        enddo
    enddo
    !
    if (mom4ice) then
        do i = 1, im
            if (flag (i)) then
                hi_save (i) = hice (i)
                hs_save (i) = weasd (i) * 0.001
            endif
        enddo
    elseif (lsm > 0) then ! --- ... snow - rain detection
        do i = 1, im
            if (flag (i)) then
                if (srflag (i) == 1.0) then
                    ep (i) = 0.0
                    weasd (i) = weasd (i) + 1.e3 * tprcp (i)
                    tprcp (i) = 0.0
                endif
            endif
        enddo
    endif

    ! --- ... initialize variables. all units are supposedly m.k.s. unless specifie
    ! psurf is in pascals, wind is wind speed, theta1 is adiabatic surface
    ! temp from level 1, rho is density, qs1 is sat. hum. at level1 and qss
    ! is sat. hum. at surface
    ! convert slrad to the civilized unit from langley minute - 1 k - 4

    do i = 1, im
        if (flag (i)) then
            ! psurf (i) = 1000.0 * ps (i)
            ! ps1 (i) = 1000.0 * prsl1 (i)

            ! dlwflx has been given a negative sign for downward longwave
            ! sfcnsw is the net shortwave flux (direction: dn - up)

            wind = max (sqrt (u1 (i) * u1 (i) + v1 (i) * v1 (i)) &
                 + max (0.0, min (ddvel (i), 30.0)), 1.0)

            q0 = max (q1 (i), 1.0e-8)
            theta1 (i) = t1 (i) * prslki (i)
            rho (i) = prsl1 (i) / (rdgas * t1 (i) * (1.0 + zvir * q0))
            qs1 = mqs (t1 (i))
            qs1 = max (eps * qs1 / (prsl1 (i) + epsm1 * qs1), 1.e-8)
            q0 = min (qs1, q0)

            ffw (i) = 1.0 - fice (i)
            if (fice (i) < cimin) then
                print *, 'warning: ice fraction is low:', fice (i)
                fice (i) = cimin
                ffw (i) = 1.0 - fice (i)
                tice (i) = tgice
                tsurf (i) = tgice
                print *, 'fix ice fraction: reset it to:', fice (i)
            endif

            qssi = mqs (tice (i))
            qssi = eps * qssi / (ps (i) + epsm1 * qssi)
            qssw = mqs (tgice)
            qssw = eps * qssw / (ps (i) + epsm1 * qssw)

            ! --- ... snow depth in water equivalent is converted from mm to m unit

            if (mom4ice) then
                snowd (i) = weasd (i) * 0.001 / fice (i)
            else
                snowd (i) = weasd (i) * 0.001
            endif
            ! flagsnw (i) = .false.

            ! --- ... when snow depth is less than 1 mm, a patchy snow is assumed and
            ! soil is allowed to interact with the atmosphere.
            ! we should eventually move to a linear combination of soil and
            ! snow under the condition of patchy snow.

            ! --- ... rcp = rho cp_air ch v

            cmm (i) = cm (i) * wind
            chh (i) = rho (i) * ch (i) * wind
            rch (i) = chh (i) * cp_air

            ! --- ... sensible and latent heat flux over open water & sea ice

            evapi (i) = elocp * rch (i) * (qssi - q0)
            evapw (i) = elocp * rch (i) * (qssw - q0)
            ! evap (i) = fice (i) * evapi (i) + ffw (i) * evapw (i)

            ! if (lprnt) write (0, *) ' tice = ', tice (ipr)

            snetw (i) = sfcdsw (i) * (1.0 - albfw)
            snetw (i) = min (3.0 * sfcnsw (i) / (1.0 + 2.0 * ffw (i)), snetw (i))
            sneti (i) = (sfcnsw (i) - ffw (i) * snetw (i)) / fice (i)

            t12 = tice (i) * tice (i)
            t14 = t12 * t12

            ! --- ... hfi = net non - solar and upir heat flux @ ice surface

            hfi (i) = - dlwflx (i) + sfcemis (i) * sbc * t14 + evapi (i) &
                 + rch (i) * (tice (i) - theta1 (i))
            hfd (i) = 4.0 * sfcemis (i) * sbc * tice (i) * t12 &
                 + (1.0 + elocp * eps * hlv * qs1 / (rdgas * t12)) * rch (i)

            t12 = tgice * tgice
            t14 = t12 * t12

            ! --- ... hfw = net heat flux @ water surface (within ice)

            ! hfw (i) = - dlwflx (i) + sfcemis (i) * sbc * t14 + evapw (i) &
            ! + rch (i) * (tgice - theta1 (i)) - snetw (i)

            focn (i) = 2.0 ! heat flux from ocean - should be from ocn model
            snof (i) = 0.0 ! snowfall rate - snow accumulates in gbphys

            hice (i) = max (min (hice (i), himax), himin)
            snowd (i) = min (snowd (i), hsmax)

            if (snowd (i) > (2.0 * hice (i))) then
                print *, 'warning: too much snow :', snowd (i)
                snowd (i) = hice (i) + hice (i)
                print *, 'fix: decrease snow depth to:', snowd (i)
            endif
        endif
    enddo

    ! if (lprnt) write (0, *) ' tice2 = ', tice (ipr)
    call ice3lay &
        ! --- inputs: !
        (im, kmi, fice, flag, hfi, hfd, sneti, focn, delt, &
        ! --- outputs: !
        snowd, hice, stsice, tice, snof, snowmt, gflux) !

    ! if (lprnt) write (0, *) ' tice3 = ', tice (ipr)
    if (mom4ice) then
        do i = 1, im
            if (flag (i)) then
                hice (i) = hi_save (i)
                snowd (i) = hs_save (i)
            endif
        enddo
    endif

    do i = 1, im
        if (flag (i)) then
            if (tice (i) < timin) then
                print *, 'warning: snow / ice temperature is too low:', tice (i), ' i = ', i
                tice (i) = timin
                print *, 'fix snow / ice temperature: reset it to:', tice (i)
            endif

            if (stsice (i, 1) < timin) then
                print *, 'warning: layer 1 ice temp is too low:', stsice (i, 1), ' i = ', i
                stsice (i, 1) = timin
                print *, 'fix layer 1 ice temp: reset it to:', stsice (i, 1)
            endif

            if (stsice (i, 2) < timin) then
                print *, 'warning: layer 2 ice temp is too low:', stsice (i, 2)
                stsice (i, 2) = timin
                print *, 'fix layer 2 ice temp: reset it to:', stsice (i, 2)
            endif

            tsurf (i) = tice (i) * fice (i) + tgice * ffw (i)
        endif
    enddo

    do k = 1, kmi
        do i = 1, im
            if (flag (i)) then
                stc (i, k) = min (stsice (i, k), t0ice)
            endif
        enddo
    enddo

    do i = 1, im
        if (flag (i)) then
            ! --- ... calculate sensible heat flux (& evap over sea ice)

            hflxi = rch (i) * (tice (i) - theta1 (i))
            hflxw = rch (i) * (tgice - theta1 (i))
            hflx (i) = fice (i) * hflxi + ffw (i) * hflxw
            evap (i) = fice (i) * evapi (i) + ffw (i) * evapw (i)
            !
            ! --- ... the rest of the output

            qsurf (i) = q1 (i) + evap (i) / (elocp * rch (i))

            ! --- ... convert snow depth back to mm of water equivalent

            weasd (i) = snowd (i) * 1000.0
            snwdph (i) = weasd (i) * dsi ! snow depth in mm

            tem = 1.0 / rho (i)
            hflx (i) = hflx (i) * tem * cpinv
            evap (i) = evap (i) * tem * hvapi
        endif
    enddo

end subroutine sfc_seai

! =======================================================================
! Three-Layer Sea Ice Vertical Thermodynamics
!
! based on:  m. winton, "a reformulated three-layer sea ice model",
! journal of atmospheric and oceanic technology, 2000
!
!
!        -> +---------+ <- tice - diagnostic surface temperature ( <= 0c )
!       /   |         |
!   snowd   |  snow   | <- 0-heat capacity snow layer
!       \   |         |
!        => +---------+
!       /   |         |
!      /    |         | <- t1 - upper 1/2 ice temperature; this layer has
!     /     |         |         a variable (t/s dependent) heat capacity
!   hice    |...ice...|
!     \     |         |
!      \    |         | <- t2 - lower 1/2 ice temp. (fixed heat capacity)
!       \   |         |
!        -> +---------+ <- base of ice fixed at seawater freezing temp.
!
! inputs: size
! im, kmi - integer, horiz dimension and num of ice layers 1
! fice - real, sea - ice concentration im
! flag - logical, ice mask flag 1
! hfi - real, net non - solar and heat flux @ surface (w / m^2) im
! hfd - real, heat flux derivatice @ sfc (w / m^2 / deg - c) im
! sneti - real, net solar incoming at top (w / m^2) im
! focn - real, heat flux from ocean (w / m^2) im
! delt - real, timestep (sec) 1
!
! input / outputs:
! snowd - real, surface pressure im
! hice - real, sea - ice thickness im
! stsice - real, temp @ midpt of ice levels (deg c) im, km
! tice - real, surface temperature (deg c) im
! snof - real, snowfall rate (m / sec) im
!
! outputs:
! snowmt - real, snow melt during delt (m) im
! gflux - real, conductive heat flux (w / m^2) im
!
! locals:
! hdi - real, ice - water interface (m)
! hsni - real, snow - ice (m)
!
! =======================================================================

subroutine ice3lay &
    !...................................
    ! --- inputs:
    (im, kmi, fice, flag, hfi, hfd, sneti, focn, delt, &
    ! --- input / outputs:
    snowd, hice, stsice, tice, snof, &
    ! --- outputs:
    snowmt, gflux)

    implicit none

    ! --- constant parameters: (properties of ice, snow, and seawater)
    real, parameter :: ds = 330.0 ! snow (ov sea ice) density (kg / m^3)
    real, parameter :: dw = 1000.0 ! fresh water density (kg / m^3)
    real, parameter :: dsdw = ds / dw
    real, parameter :: dwds = dw / ds
    real, parameter :: ks = 0.31 ! conductivity of snow (w / mk)
    real, parameter :: i0 = 0.3 ! ice surface penetrating solar fraction
    real, parameter :: ki = 2.03 ! conductivity of ice (w / mk)
    real, parameter :: di = 917.0 ! density of ice (kg / m^3)
    real, parameter :: didw = di / dw
    real, parameter :: dsdi = ds / di
    real, parameter :: ci = 2054.0 ! heat capacity of fresh ice (j / kg / k)
    real, parameter :: li = 3.34e5 ! latent heat of fusion (j / kg - ice)
    real, parameter :: si = 1.0 ! salinity of sea ice
    real, parameter :: mu = 0.054 ! relates freezing temp to salinity
    real, parameter :: tfi = - mu * si ! sea ice freezing temp = - mu * salinity
    real, parameter :: tfw = - 1.8 ! tfw - seawater freezing temp (c)
    real, parameter :: tfi0 = tfi - 0.0001
    real, parameter :: dici = di * ci
    real, parameter :: dili = di * li
    real, parameter :: dsli = ds * li
    real, parameter :: ki4 = ki * 4.0

    ! --- inputs:
    integer, intent (in) :: im, kmi

    real, dimension (im), intent (in) :: fice, hfi, hfd, sneti, focn

    real, intent (in) :: delt

    logical, dimension (im), intent (in) :: flag

    ! --- input / outputs:
    real, dimension (im), intent (inout) :: snowd, hice, tice, snof

    real, dimension (im, kmi), intent (inout) :: stsice

    ! --- outputs:
    real, dimension (im), intent (out) :: snowmt, gflux

    ! --- locals:

    real :: dt2, dt4, dt6, h1, h2, dh, wrk, wrk1, &
        dt2i, hdi, hsni, ai, bi, a1, b1, a10, b10, &
        c1, ip, k12, k32, tsf, f1, tmelt, bmelt

    integer :: i
    !
    ! ===> ... begin here
    !
    dt2 = 2.0 * delt
    dt4 = 4.0 * delt
    dt6 = 6.0 * delt
    dt2i = 1.0 / dt2

    do i = 1, im
        if (flag (i)) then
            snowd (i) = snowd (i) * dwds
            hdi = (dsdw * snowd (i) + didw * hice (i))

            if (hice (i) < hdi) then
                snowd (i) = snowd (i) + hice (i) - hdi
                hsni = (hdi - hice (i)) * dsdi
                hice (i) = hice (i) + hsni
            endif

            snof (i) = snof (i) * dwds
            tice (i) = tice (i) - t0ice
            stsice (i, 1) = min (stsice (i, 1) - t0ice, tfi0) ! degc
            stsice (i, 2) = min (stsice (i, 2) - t0ice, tfi0) ! degc

            ip = i0 * sneti (i) ! ip + v (in winton ip = - i0 * sneti as sol - v)
            if (snowd (i) > 0.0) then
                tsf = 0.0
                ip = 0.0
            else
                tsf = tfi
                ip = i0 * sneti (i) ! ip + v here (in winton ip = - i0 * sneti)
            endif
            tice (i) = min (tice (i), tsf)

            ! --- ... compute ice temperature

            bi = hfd (i)
            ai = hfi (i) - sneti (i) + ip - tice (i) * bi ! + v sol input here
            k12 = ki4 * ks / (ks * hice (i) + ki4 * snowd (i))
            k32 = (ki + ki) / hice (i)

            wrk = 1.0 / (dt6 * k32 + dici * hice (i))
            a10 = dici * hice (i) * dt2i + k32 * (dt4 * k32 + dici * hice (i)) * wrk
            b10 = - di * hice (i) * (ci * stsice (i, 1) + li * tfi / stsice (i, 1)) &
                 * dt2i - ip &
                 - k32 * (dt4 * k32 * tfw + dici * hice (i) * stsice (i, 2)) * wrk

            wrk1 = k12 / (k12 + bi)
            a1 = a10 + bi * wrk1
            b1 = b10 + ai * wrk1
            c1 = dili * tfi * dt2i * hice (i)

            stsice (i, 1) = - (sqrt (b1 * b1 - 4.0 * a1 * c1) + b1) / (a1 + a1)
            tice (i) = (k12 * stsice (i, 1) - ai) / (k12 + bi)

            if (tice (i) > tsf) then
                a1 = a10 + k12
                b1 = b10 - k12 * tsf
                stsice (i, 1) = - (sqrt (b1 * b1 - 4.0 * a1 * c1) + b1) / (a1 + a1)
                tice (i) = tsf
                tmelt = (k12 * (stsice (i, 1) - tsf) - (ai + bi * tsf)) * delt
            else
                tmelt = 0.0
                snowd (i) = snowd (i) + snof (i) * delt
            endif

            stsice (i, 2) = (dt2 * k32 * (stsice (i, 1) + tfw + tfw) &
                 + dici * hice (i) * stsice (i, 2)) * wrk

            bmelt = (focn (i) + ki4 * (stsice (i, 2) - tfw) / hice (i)) * delt

            ! --- ... resize the ice ...

            h1 = 0.5 * hice (i)
            h2 = 0.5 * hice (i)

            ! --- ... top ...

            if (tmelt <= snowd (i) * dsli) then
                snowmt (i) = tmelt / dsli
                snowd (i) = snowd (i) - snowmt (i)
            else
                snowmt (i) = snowd (i)
                h1 = h1 - (tmelt - snowd (i) * dsli) &
                     / (di * (ci - li / stsice (i, 1)) * (tfi - stsice (i, 1)))
                snowd (i) = 0.0
            endif

            ! --- ... and bottom

            if (bmelt < 0.0) then
                dh = - bmelt / (dili + dici * (tfi - tfw))
                stsice (i, 2) = (h2 * stsice (i, 2) + dh * tfw) / (h2 + dh)
                h2 = h2 + dh
            else
                h2 = h2 - bmelt / (dili + dici * (tfi - stsice (i, 2)))
            endif

            ! --- ... if ice remains, even up 2 layers, else, pass negative energy back in snow

            hice (i) = h1 + h2

            if (hice (i) > 0.0) then
                if (h1 > 0.5 * hice (i)) then
                    f1 = 1.0 - (h2 + h2) / hice (i)
                    stsice (i, 2) = f1 * (stsice (i, 1) + li * tfi / (ci * stsice (i, 1))) &
                         + (1.0 - f1) * stsice (i, 2)

                    if (stsice (i, 2) > tfi) then
                        hice (i) = hice (i) - h2 * ci * (stsice (i, 2) - tfi) / (li * delt)
                        stsice (i, 2) = tfi
                    endif
                else
                    f1 = (h1 + h1) / hice (i)
                    stsice (i, 1) = f1 * (stsice (i, 1) + li * tfi / (ci * stsice (i, 1))) &
                         + (1.0 - f1) * stsice (i, 2)
                    stsice (i, 1) = (stsice (i, 1) - sqrt (stsice (i, 1) * stsice (i, 1) &
                         - 4.0 * tfi * li / ci)) * 0.5
                endif

                k12 = ki4 * ks / (ks * hice (i) + ki4 * snowd (i))
                gflux (i) = k12 * (stsice (i, 1) - tice (i))
            else
                snowd (i) = snowd (i) + (h1 * (ci * (stsice (i, 1) - tfi) &
                     - li * (1.0 - tfi / stsice (i, 1))) &
                     + h2 * (ci * (stsice (i, 2) - tfi) - li)) / li

                hice (i) = max (0.0, snowd (i) * dsdi)
                snowd (i) = 0.0
                stsice (i, 1) = tfw
                stsice (i, 2) = tfw
                gflux (i) = 0.0
            endif ! endif_hice_block

            gflux (i) = fice (i) * gflux (i)
            snowmt (i) = snowmt (i) * dsdw
            snowd (i) = snowd (i) * dsdw
            tice (i) = tice (i) + t0ice
            stsice (i, 1) = stsice (i, 1) + t0ice
            stsice (i, 2) = stsice (i, 2) + t0ice
        endif ! endif_flag_block
    enddo ! enddo_i_loop

end subroutine ice3lay

! =======================================================================
! subroutine to update near surface fields
! =======================================================================

subroutine sfc_updt (im, ps, u1, v1, t1, q1, &
        tsurf, qsurf, u10m, v10m, t2m, q2m, &
        prslki, evap, fm, fh, fm10, fh2)

    implicit none

    integer im
    real, dimension (im) :: ps, u1, v1, t1, q1, tsurf, qsurf, &
        u10m, v10m, t2m, q2m, prslki, evap, &
        fm, fh, fm10, fh2

    ! locals

    real, parameter :: qmin = 1.0e-8
    integer k, i

    real :: fhi, qss, wrk, f10m (im)
    ! real :: sig2k, fhi, qss

    ! real, parameter :: g = grav

    ! estimate sigma ** k at 2 m

    ! sig2k = 1. - 4. * g * 2. / (cp_air * 280.)

    ! initialize variables. all units are supposedly m.k.s. unless specified
    ! ps is in pascals

    do i = 1, im
        f10m (i) = fm10 (i) / fm (i)
        ! f10m (i) = min (f10m (i), 1.)
        u10m (i) = f10m (i) * u1 (i)
        v10m (i) = f10m (i) * v1 (i)
        fhi = fh2 (i) / fh (i)
        ! t2m (i) = tsurf (i) * (1. - fhi) + t1 (i) * prslki (i) * fhi
        ! sig2k = 1. - (grav + grav) / (cp_air * t2m (i))
        ! t2m (i) = t2m (i) * sig2k
        wrk = 1.0 - fhi

        t2m (i) = tsurf (i) * wrk + t1 (i) * prslki (i) * fhi - (grav + grav) / cp_air

        if (evap (i) >= 0.) then ! for evaporation > 0, use inferred qsurf to deduce q2m
            q2m (i) = qsurf (i) * wrk + max (qmin, q1 (i)) * fhi
        else ! for dew formation, use saturated q at tsurf
            qss = mqs (tsurf (i))
            qss = eps * qss / (ps (i) + epsm1 * qss)
            q2m (i) = qss * wrk + max (qmin, q1 (i)) * fhi
        endif
        qss = mqs (t2m (i))
        qss = eps * qss / (ps (i) + epsm1 * qss)
        q2m (i) = min (q2m (i), qss)
    enddo

end subroutine sfc_updt

! =======================================================================
! solve tridiagonal problem for tke
! =======================================================================

subroutine tridit (l, n, nt, cl, cm, cu, rt, au, at)

    implicit none

    integer :: is, k, kk, n, nt, l, i

    real :: fk (l)

    real :: cl (l, 2:n), cm (l, n), cu (l, n - 1), &
        rt (l, n * nt), &
        au (l, n - 1), at (l, n * nt), &
        fkk (l, 2:n - 1)

    do i = 1, l
        fk (i) = 1. / cm (i, 1)
        au (i, 1) = fk (i) * cu (i, 1)
    enddo
    do k = 1, nt
        is = (k - 1) * n
        do i = 1, l
            at (i, 1 + is) = fk (i) * rt (i, 1 + is)
        enddo
    enddo
    do k = 2, n - 1
        do i = 1, l
            fkk (i, k) = 1. / (cm (i, k) - cl (i, k) * au (i, k - 1))
            au (i, k) = fkk (i, k) * cu (i, k)
        enddo
    enddo
    do kk = 1, nt
        is = (kk - 1) * n
        do k = 2, n - 1
            do i = 1, l
                at (i, k + is) = fkk (i, k) * (rt (i, k + is) - cl (i, k) * at (i, k + is - 1))
            enddo
        enddo
    enddo
    do i = 1, l
        fk (i) = 1. / (cm (i, n) - cl (i, n) * au (i, n - 1))
    enddo
    do k = 1, nt
        is = (k - 1) * n
        do i = 1, l
            at (i, n + is) = fk (i) * (rt (i, n + is) - cl (i, n) * at (i, n + is - 1))
        enddo
    enddo
    do kk = 1, nt
        is = (kk - 1) * n
        do k = n - 1, 1, - 1
            do i = 1, l
                at (i, k + is) = at (i, k + is) - au (i, k) * at (i, k + is + 1)
            enddo
        enddo
    enddo

end subroutine tridit

! =======================================================================
! edmf parameterization siebesma et al. (2007)
! =======================================================================

subroutine mfpblt (im, km, kmpbl, ntcw, ntrac1, delt, &
        cnvflg, zl, zm, q1, t1, u1, v1, plyr, pix, thlx, thvx, &
        gdx, hpbl, kpbl, vpert, buo, xmf, &
        tcko, qcko, ucko, vcko, xlamue)

    implicit none

    integer, intent (in) :: im, km, kmpbl, ntcw, ntrac1
    integer :: kpbl (im)

    logical :: cnvflg (im)

    real :: delt
    real :: q1 (im, km, ntrac1), &
        t1 (im, km), u1 (im, km), v1 (im, km), &
        plyr (im, km), pix (im, km), thlx (im, km), &
        thvx (im, km), zl (im, km), zm (im, km), &
        gdx (im), &
        hpbl (im), vpert (im), &
        buo (im, km), xmf (im, km), &
        tcko (im, km), qcko (im, km, ntrac1), &
        ucko (im, km), vcko (im, km), &
        xlamue (im, km - 1)

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: i, j, k, n, ndc
    integer :: kpblx (im), kpbly (im)

    real :: dt2, dz, ce0, cm, &
        factor, gocp, &
        g, b1, f1, &
        bb1, bb2, &
        alp, a1, pgcon, &
        qmin, qlmin, xmmx, rbint, &
        tem, tem1, tem2, &
        ptem, ptem1, ptem2

    real :: elocp, el2orc, qs, es, &
        tlu, gamma, qlu, &
        thup, thvu, dq

    real :: rbdn (im), rbup (im), hpblx (im), &
        xlamuem (im, km - 1)

    real :: wu2 (im, km), thlu (im, km), &
        qtx (im, km), qtu (im, km)

    real :: xlamavg (im), sigma (im), &
        scaldfunc (im), sumx (im)

    logical :: totflg, flg (im)

    ! physical parameters
    parameter (g = grav)
    parameter (gocp = g / cp_air)
    parameter (elocp = hlv / cp_air, el2orc = hlv * hlv / (rvgas * cp_air))
    parameter (ce0 = 0.4, cm = 1.0)
    parameter (qmin = 1.e-8, qlmin = 1.e-12)
    parameter (alp = 1.0, pgcon = 0.55)
    parameter (a1 = 0.13, b1 = 0.5, f1 = 0.15)

    totflg = .true.
    do i = 1, im
        totflg = totflg .and. (.not. cnvflg (i))
    enddo
    if (totflg) return

    dt2 = delt

    do k = 1, km
        do i = 1, im
            if (cnvflg (i)) then
                buo (i, k) = 0.
                wu2 (i, k) = 0.
                qtx (i, k) = q1 (i, k, 1) + q1 (i, k, ntcw)
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute thermal excess
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            ptem = alp * vpert (i)
            ptem = min (ptem, 3.0)
            thlu (i, 1) = thlx (i, 1) + ptem
            qtu (i, 1) = qtx (i, 1)
            buo (i, 1) = g * ptem / thvx (i, 1)
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! compute entrainment rate
    ! -----------------------------------------------------------------------

    do k = 1, kmpbl
        do i = 1, im
            if (cnvflg (i)) then
                dz = zl (i, k + 1) - zl (i, k)
                if (k < kpbl (i)) then
                    ptem = 1. / (zm (i, k) + dz)
                    tem = max ((hpbl (i) - zm (i, k) + dz), dz)
                    ptem1 = 1. / tem
                    xlamue (i, k) = ce0 * (ptem + ptem1)
                else
                    xlamue (i, k) = ce0 / dz
                endif
                xlamuem (i, k) = cm * xlamue (i, k)
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute buoyancy for updraft air parcel
    ! -----------------------------------------------------------------------

    do k = 2, kmpbl
        do i = 1, im
            if (cnvflg (i)) then
                dz = zl (i, k) - zl (i, k - 1)
                tem = 0.5 * xlamue (i, k - 1) * dz
                factor = 1. + tem

                thlu (i, k) = ((1. - tem) * thlu (i, k - 1) + tem * &
                    (thlx (i, k - 1) + thlx (i, k))) / factor
                qtu (i, k) = ((1. - tem) * qtu (i, k - 1) + tem * &
                     (qtx (i, k - 1) + qtx (i, k))) / factor

                tlu = thlu (i, k) / pix (i, k)
                es = 0.01 * mqs (tlu) ! mqs in pa
                qs = max (qmin, eps * es / (plyr (i, k) + epsm1 * es))
                dq = qtu (i, k) - qs

                if (dq > 0.) then
                    gamma = el2orc * qs / (tlu ** 2)
                    qlu = dq / (1. + gamma)
                    qtu (i, k) = qs + qlu
                    tem1 = 1. + zvir * qs - qlu
                    thup = thlu (i, k) + pix (i, k) * elocp * qlu
                    thvu = thup * tem1
                else
                    tem1 = 1. + zvir * qtu (i, k)
                    thvu = thlu (i, k) * tem1
                endif
                buo (i, k) = g * (thvu / thvx (i, k) - 1.)

            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute updraft velocity square (wu2)
    ! -----------------------------------------------------------------------

    ! tem = 1. - 2. * f1
    ! bb1 = 2. * b1 / tem
    ! bb2 = 2. / tem
    ! from soares et al. (2004, qjrms)
    ! bb1 = 2.
    ! bb2 = 4.

    ! from bretherton et al. (2004, mwr)
    ! bb1 = 4.
    ! bb2 = 2.

    ! from our tuning
    bb1 = 2.0
    bb2 = 4.0

    do i = 1, im
        if (cnvflg (i)) then
            dz = zm (i, 1)
            tem = 0.5 * bb1 * xlamue (i, 1) * dz
            tem1 = bb2 * buo (i, 1) * dz
            ptem1 = 1. + tem
            wu2 (i, 1) = tem1 / ptem1
        endif
    enddo
    do k = 2, kmpbl
        do i = 1, im
            if (cnvflg (i)) then
                dz = zm (i, k) - zm (i, k - 1)
                tem = 0.25 * bb1 * (xlamue (i, k) + xlamue (i, k - 1)) * dz
                tem1 = bb2 * buo (i, k) * dz
                ptem = (1. - tem) * wu2 (i, k - 1)
                ptem1 = 1. + tem
                wu2 (i, k) = (ptem + tem1) / ptem1
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! update pbl height as the height where updraft velocity vanishes
    ! -----------------------------------------------------------------------

    do i = 1, im
        flg (i) = .true.
        kpbly (i) = kpbl (i)
        if (cnvflg (i)) then
            flg (i) = .false.
            rbup (i) = wu2 (i, 1)
        endif
    enddo
    do k = 2, kmpbl
        do i = 1, im
            if (.not.flg (i)) then
                rbdn (i) = rbup (i)
                rbup (i) = wu2 (i, k)
                kpblx (i) = k
                flg (i) = rbup (i) .le.0.
            endif
        enddo
    enddo
    do i = 1, im
        if (cnvflg (i)) then
            k = kpblx (i)
            if (rbdn (i) <= 0.) then
                rbint = 0.
            elseif (rbup (i) >= 0.) then
                rbint = 1.
            else
                rbint = rbdn (i) / (rbdn (i) - rbup (i))
            endif
            hpblx (i) = zm (i, k - 1) + rbint * (zm (i, k) - zm (i, k - 1))
        endif
    enddo

    do i = 1, im
        if (cnvflg (i)) then
            if (kpbl (i) > kpblx (i)) then
                kpbl (i) = kpblx (i)
                hpbl (i) = hpblx (i)
            endif
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! update entrainment rate
    ! -----------------------------------------------------------------------

    do k = 1, kmpbl
        do i = 1, im
            if (cnvflg (i) .and. kpbly (i) > kpblx (i)) then
                dz = zl (i, k + 1) - zl (i, k)
                if (k < kpbl (i)) then
                    ptem = 1. / (zm (i, k) + dz)
                    tem = max ((hpbl (i) - zm (i, k) + dz), dz)
                    ptem1 = 1. / tem
                    xlamue (i, k) = ce0 * (ptem + ptem1)
                else
                    xlamue (i, k) = ce0 / dz
                endif
                xlamuem (i, k) = cm * xlamue (i, k)
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute entrainment rate averaged over the whole pbl
    ! -----------------------------------------------------------------------

    do i = 1, im
        xlamavg (i) = 0.
        sumx (i) = 0.
    enddo
    do k = 1, kmpbl
        do i = 1, im
            if (cnvflg (i) .and. k < kpbl (i)) then
                dz = zl (i, k + 1) - zl (i, k)
                xlamavg (i) = xlamavg (i) + xlamue (i, k) * dz
                sumx (i) = sumx (i) + dz
            endif
        enddo
    enddo
    do i = 1, im
        if (cnvflg (i)) then
            xlamavg (i) = xlamavg (i) / sumx (i)
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! updraft mass flux as a function of updraft velocity profile
    ! -----------------------------------------------------------------------

    do k = 1, kmpbl
        do i = 1, im
            if (cnvflg (i) .and. k < kpbl (i)) then
                if (wu2 (i, k) > 0.) then
                    tem = sqrt (wu2 (i, k))
                else
                    tem = 0.
                endif
                xmf (i, k) = a1 * tem
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute updraft fraction as a function of mean entrainment rate
    ! (grell & freitas, 2014)
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            tem = 0.2 / xlamavg (i)
            tem1 = 3.14 * tem * tem
            sigma (i) = tem1 / (gdx (i) * gdx (i))
            sigma (i) = max (sigma (i), 0.001)
            sigma (i) = min (sigma (i), 0.999)
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! compute scale - aware function based on arakawa & wu (2013)
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            if (sigma (i) > a1) then
                scaldfunc (i) = (1. - sigma (i)) * (1. - sigma (i))
                scaldfunc (i) = max (min (scaldfunc (i), 1.0), 0.)
            else
                scaldfunc (i) = 1.0
            endif
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! final scale - aware updraft mass flux
    ! -----------------------------------------------------------------------

    do k = 1, kmpbl
        do i = 1, im
            if (cnvflg (i) .and. k < kpbl (i)) then
                xmf (i, k) = scaldfunc (i) * xmf (i, k)
                dz = zl (i, k + 1) - zl (i, k)
                xmmx = dz / dt2
                xmf (i, k) = min (xmf (i, k), xmmx)
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute updraft property using updated entranment rate
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            thlu (i, 1) = thlx (i, 1)
        endif
    enddo

    ! do i = 1, im
    ! if (cnvflg (i)) then
    ! ptem1 = max (qcko (i, 1, ntcw), 0.)
    ! tlu = thlu (i, 1) / pix (i, 1)
    ! tcko (i, 1) = tlu + elocp * ptem1
    ! endif
    ! enddo

    do k = 2, kmpbl
        do i = 1, im
            if (cnvflg (i) .and. k <= kpbl (i)) then
                dz = zl (i, k) - zl (i, k - 1)
                tem = 0.5 * xlamue (i, k - 1) * dz
                factor = 1. + tem

                thlu (i, k) = ((1. - tem) * thlu (i, k - 1) + tem * &
                    (thlx (i, k - 1) + thlx (i, k))) / factor
                qtu (i, k) = ((1. - tem) * qtu (i, k - 1) + tem * &
                    (qtx (i, k - 1) + qtx (i, k))) / factor

                tlu = thlu (i, k) / pix (i, k)
                es = 0.01 * mqs (tlu) ! mqs in pa
                qs = max (qmin, eps * es / (plyr (i, k) + epsm1 * es))
                dq = qtu (i, k) - qs

                if (dq > 0.) then
                    gamma = el2orc * qs / (tlu ** 2)
                    qlu = dq / (1. + gamma)
                    qtu (i, k) = qs + qlu
                    qcko (i, k, 1) = qs
                    qcko (i, k, ntcw) = qlu
                    tcko (i, k) = tlu + elocp * qlu
                else
                    qcko (i, k, 1) = qtu (i, k)
                    qcko (i, k, ntcw) = 0.
                    tcko (i, k) = tlu
                endif

            endif
        enddo
    enddo

    do k = 2, kmpbl
        do i = 1, im
            if (cnvflg (i) .and. k <= kpbl (i)) then
                dz = zl (i, k) - zl (i, k - 1)
                tem = 0.5 * xlamuem (i, k - 1) * dz
                factor = 1. + tem
                ptem = tem + pgcon
                ptem1 = tem - pgcon
                ucko (i, k) = ((1. - tem) * ucko (i, k - 1) + ptem * u1 (i, k) + &
                    ptem1 * u1 (i, k - 1)) / factor
                vcko (i, k) = ((1. - tem) * vcko (i, k - 1) + ptem * v1 (i, k) + &
                    ptem1 * v1 (i, k - 1)) / factor
            endif
        enddo
    enddo

    if (ntcw > 2) then

        do n = 2, ntcw - 1
            do k = 2, kmpbl
                do i = 1, im
                    if (cnvflg (i) .and. k <= kpbl (i)) then
                        dz = zl (i, k) - zl (i, k - 1)
                        tem = 0.5 * xlamue (i, k - 1) * dz
                        factor = 1. + tem

                        qcko (i, k, n) = ((1. - tem) * qcko (i, k - 1, n) + tem * &
                            (q1 (i, k, n) + q1 (i, k - 1, n))) / factor
                    endif
                enddo
            enddo
        enddo

    endif

    ndc = ntrac1 - ntcw

    if (ndc > 0) then

        do n = ntcw + 1, ntrac1
            do k = 2, kmpbl
                do i = 1, im
                    if (cnvflg (i) .and. k <= kpbl (i)) then
                        dz = zl (i, k) - zl (i, k - 1)
                        tem = 0.5 * xlamue (i, k - 1) * dz
                        factor = 1. + tem

                        qcko (i, k, n) = ((1. - tem) * qcko (i, k - 1, n) + tem * &
                            (q1 (i, k, n) + q1 (i, k - 1, n))) / factor
                    endif
                enddo
            enddo
        enddo

    endif

    return

end subroutine mfpblt

! =======================================================================
! mass - flux parameterization for stratocumulus - top - induced turbulence mixing
! =======================================================================

subroutine mfscu (im, km, kmscu, ntcw, ntrac1, delt, &
        cnvflg, zl, zm, q1, t1, u1, v1, plyr, pix, &
        thlx, thvx, thlvx, gdx, thetae, radj, &
        krad, mrad, radmin, buo, xmfd, &
        tcdo, qcdo, ucdo, vcdo, xlamde)

    implicit none

    integer, intent (in) :: im, km, kmscu, ntcw, ntrac1
    integer :: krad (im), mrad (im)

    logical :: cnvflg (im)

    real :: delt
    real :: q1 (im, km, ntrac1), t1 (im, km), &
        u1 (im, km), v1 (im, km), &
        plyr (im, km), pix (im, km), &
        thlx (im, km), &
        thvx (im, km), thlvx (im, km), &
        gdx (im), radj (im), &
        zl (im, km), zm (im, km), &
        thetae (im, km), radmin (im), &
        buo (im, km), xmfd (im, km), &
        tcdo (im, km), qcdo (im, km, ntrac1), &
        ucdo (im, km), vcdo (im, km), &
        xlamde (im, km - 1)

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: i, j, indx, k, n, kk, ndc

    integer :: krad1 (im), mradx (im), mrady (im)

    real :: dt2, dz, ce0, cm, &
        gocp, factor, g, tau, &
        b1, f1, bb1, bb2, &
        a1, a2, a11, a22, &
        cteit, pgcon, &
        qmin, qlmin, &
        xmmx, tem, tem1, tem2, &
        ptem, ptem1, ptem2

    real :: elocp, el2orc, qs, es, &
        tld, gamma, qld, thdn, &
        thvd, dq

    real :: wd2 (im, km), thld (im, km), &
        qtx (im, km), qtd (im, km), &
        thlvd (im), hrad (im), &
        xlamdem (im, km - 1), ra1 (im), ra2 (im)

    real :: xlamavg (im), sigma (im), &
        scaldfunc (im), sumx (im)

    logical :: totflg, flg (im)

    real :: actei, cldtime

    ! physical parameters
    parameter (g = grav)
    parameter (gocp = g / cp_air)
    parameter (elocp = hlv / cp_air, el2orc = hlv * hlv / (rvgas * cp_air))
    parameter (ce0 = 0.4, cm = 1.0, pgcon = 0.55)
    parameter (qmin = 1.e-8, qlmin = 1.e-12)
    parameter (b1 = 0.45, f1 = 0.15)
    parameter (a1 = 0.12, a2 = 0.5)
    parameter (a11 = 0.2, a22 = 1.0)
    parameter (cldtime = 500.)
    parameter (actei = 0.7)
    ! parameter (actei = 0.23)

    totflg = .true.
    do i = 1, im
        totflg = totflg .and. (.not. cnvflg (i))
    enddo
    if (totflg) return

    dt2 = delt

    do k = 1, km
        do i = 1, im
            if (cnvflg (i)) then
                buo (i, k) = 0.
                wd2 (i, k) = 0.
                qtx (i, k) = q1 (i, k, 1) + q1 (i, k, ntcw)
            endif
        enddo
    enddo

    do i = 1, im
        if (cnvflg (i)) then
            hrad (i) = zm (i, krad (i))
            krad1 (i) = krad (i) - 1
        endif
    enddo

    do i = 1, im
        if (cnvflg (i)) then
            k = krad (i)
            tem = zm (i, k + 1) - zm (i, k)
            tem1 = cldtime * radmin (i) / tem
            tem1 = max (tem1, - 3.0)
            thld (i, k) = thlx (i, k) + tem1
            qtd (i, k) = qtx (i, k)
            thlvd (i) = thlvx (i, k) + tem1
            buo (i, k) = - g * tem1 / thvx (i, k)
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! specify downdraft fraction
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            ra1 (i) = a1
            ra2 (i) = a11
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! if the condition for cloud - top instability is met,
    ! increase downdraft fraction
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            k = krad (i)
            tem = thetae (i, k) - thetae (i, k + 1)
            tem1 = qtx (i, k) - qtx (i, k + 1)
            if (tem > 0. .and. tem1 > 0.) then
                cteit = cp_air * tem / (hlv * tem1)
                if (cteit > actei) then
                    ra1 (i) = a2
                    ra2 (i) = a22
                endif
            endif
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! compute radiative flux jump at stratocumulus top
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            radj (i) = - ra2 (i) * radmin (i)
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! first - quess level of downdraft extension (mrad)
    ! -----------------------------------------------------------------------

    do i = 1, im
        flg (i) = cnvflg (i)
        mrad (i) = krad (i)
    enddo
    do k = kmscu, 1, - 1
        do i = 1, im
            if (flg (i) .and. k < krad (i)) then
                if (thlvd (i) <= thlvx (i, k)) then
                    mrad (i) = k
                else
                    flg (i) = .false.
                endif
            endif
        enddo
    enddo
    do i = 1, im
        if (cnvflg (i)) then
            kk = krad (i) - mrad (i)
            if (kk < 1) cnvflg (i) = .false.
        endif
    enddo

    totflg = .true.
    do i = 1, im
        totflg = totflg .and. (.not. cnvflg (i))
    enddo
    if (totflg) return

    ! -----------------------------------------------------------------------
    ! compute entrainment rate
    ! -----------------------------------------------------------------------

    do k = 1, kmscu
        do i = 1, im
            if (cnvflg (i)) then
                dz = zl (i, k + 1) - zl (i, k)
                if (k >= mrad (i) .and. k < krad (i)) then
                    if (mrad (i) == 1) then
                        ptem = 1. / (zm (i, k) + dz)
                    else
                        ptem = 1. / (zm (i, k) - zm (i, mrad (i) - 1) + dz)
                    endif
                    tem = max ((hrad (i) - zm (i, k) + dz), dz)
                    ptem1 = 1. / tem
                    xlamde (i, k) = ce0 * (ptem + ptem1)
                else
                    xlamde (i, k) = ce0 / dz
                endif
                xlamdem (i, k) = cm * xlamde (i, k)
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute buoyancy for downdraft air parcel
    ! -----------------------------------------------------------------------

    do k = kmscu, 1, - 1
        do i = 1, im
            if (cnvflg (i) .and. k < krad (i)) then
                dz = zl (i, k + 1) - zl (i, k)
                tem = 0.5 * xlamde (i, k) * dz
                factor = 1. + tem

                thld (i, k) = ((1. - tem) * thld (i, k + 1) + tem * &
                    (thlx (i, k) + thlx (i, k + 1))) / factor
                qtd (i, k) = ((1. - tem) * qtd (i, k + 1) + tem * &
                    (qtx (i, k) + qtx (i, k + 1))) / factor

                tld = thld (i, k) / pix (i, k)
                es = 0.01 * mqs (tld) ! mqs in pa
                qs = max (qmin, eps * es / (plyr (i, k) + epsm1 * es))
                dq = qtd (i, k) - qs

                if (dq > 0.) then
                    gamma = el2orc * qs / (tld ** 2)
                    qld = dq / (1. + gamma)
                    qtd (i, k) = qs + qld
                    tem1 = 1. + zvir * qs - qld
                    thdn = thld (i, k) + pix (i, k) * elocp * qld
                    thvd = thdn * tem1
                else
                    tem1 = 1. + zvir * qtd (i, k)
                    thvd = thld (i, k) * tem1
                endif
                buo (i, k) = g * (1. - thvd / thvx (i, k))

            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute downdraft velocity square (wd2)
    ! -----------------------------------------------------------------------

    ! tem = 1. - 2. * f1
    ! bb1 = 2. * b1 / tem
    ! bb2 = 2. / tem
    ! from soares et al. (2004, qjrms)
    ! bb1 = 2.
    ! bb2 = 4.

    ! from bretherton et al. (2004, mwr)
    ! bb1 = 4.
    ! bb2 = 2.

    ! from our tuning
    bb1 = 2.0
    bb2 = 4.0

    do i = 1, im
        if (cnvflg (i)) then
            k = krad1 (i)
            dz = zm (i, k + 1) - zm (i, k)
            ! tem = 0.25 * bb1 * (xlamde (i, k) + xlamde (i, k + 1)) * dz
            tem = 0.5 * bb1 * xlamde (i, k) * dz
            tem1 = bb2 * buo (i, k + 1) * dz
            ptem1 = 1. + tem
            wd2 (i, k) = tem1 / ptem1
        endif
    enddo
    do k = kmscu, 1, - 1
        do i = 1, im
            if (cnvflg (i) .and. k < krad1 (i)) then
                dz = zm (i, k + 1) - zm (i, k)
                tem = 0.25 * bb1 * (xlamde (i, k) + xlamde (i, k + 1)) * dz
                tem1 = bb2 * buo (i, k + 1) * dz
                ptem = (1. - tem) * wd2 (i, k + 1)
                ptem1 = 1. + tem
                wd2 (i, k) = (ptem + tem1) / ptem1
            endif
        enddo
    enddo

    do i = 1, im
        flg (i) = cnvflg (i)
        mrady (i) = mrad (i)
        if (flg (i)) mradx (i) = krad (i)
    enddo
    do k = kmscu, 1, - 1
        do i = 1, im
            if (flg (i) .and. k < krad (i)) then
                if (wd2 (i, k) > 0.) then
                    mradx (i) = k
                else
                    flg (i) = .false.
                endif
            endif
        enddo
    enddo

    do i = 1, im
        if (cnvflg (i)) then
            if (mrad (i) < mradx (i)) then
                mrad (i) = mradx (i)
            endif
        endif
    enddo

    do i = 1, im
        if (cnvflg (i)) then
            kk = krad (i) - mrad (i)
            if (kk < 1) cnvflg (i) = .false.
        endif
    enddo

    totflg = .true.
    do i = 1, im
        totflg = totflg .and. (.not. cnvflg (i))
    enddo
    if (totflg) return

    ! -----------------------------------------------------------------------
    ! update entrainment rate
    ! -----------------------------------------------------------------------

    do k = 1, kmscu
        do i = 1, im
            if (cnvflg (i) .and. mrady (i) < mradx (i)) then
                dz = zl (i, k + 1) - zl (i, k)
                if (k >= mrad (i) .and. k < krad (i)) then
                    if (mrad (i) == 1) then
                        ptem = 1. / (zm (i, k) + dz)
                    else
                        ptem = 1. / (zm (i, k) - zm (i, mrad (i) - 1) + dz)
                    endif
                    tem = max ((hrad (i) - zm (i, k) + dz), dz)
                    ptem1 = 1. / tem
                    xlamde (i, k) = ce0 * (ptem + ptem1)
                else
                    xlamde (i, k) = ce0 / dz
                endif
                xlamdem (i, k) = cm * xlamde (i, k)
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute entrainment rate averaged over the whole downdraft layers
    ! -----------------------------------------------------------------------

    do i = 1, im
        xlamavg (i) = 0.
        sumx (i) = 0.
    enddo
    do k = kmscu, 1, - 1
        do i = 1, im
            if (cnvflg (i) .and. (k >= mrad (i) .and. k < krad (i))) then
                dz = zl (i, k + 1) - zl (i, k)
                xlamavg (i) = xlamavg (i) + xlamde (i, k) * dz
                sumx (i) = sumx (i) + dz
            endif
        enddo
    enddo
    do i = 1, im
        if (cnvflg (i)) then
            xlamavg (i) = xlamavg (i) / sumx (i)
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! compute downdraft mass flux
    ! -----------------------------------------------------------------------

    do k = kmscu, 1, - 1
        do i = 1, im
            if (cnvflg (i) .and. (k >= mrad (i) .and. k < krad (i))) then
                if (wd2 (i, k) > 0.) then
                    tem = sqrt (wd2 (i, k))
                else
                    tem = 0.
                endif
                xmfd (i, k) = ra1 (i) * tem
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute downdraft fraction as a function of mean entrainment rate
    ! (grell & freitas, 2014)
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            tem = 0.2 / xlamavg (i)
            tem1 = 3.14 * tem * tem
            sigma (i) = tem1 / (gdx (i) * gdx (i))
            sigma (i) = max (sigma (i), 0.001)
            sigma (i) = min (sigma (i), 0.999)
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! compute scale - aware function based on arakawa & wu (2013)
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            if (sigma (i) > ra1 (i)) then
                scaldfunc (i) = (1. - sigma (i)) * (1. - sigma (i))
                scaldfunc (i) = max (min (scaldfunc (i), 1.0), 0.)
            else
                scaldfunc (i) = 1.0
            endif
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! final scale - aware downdraft mass flux
    ! -----------------------------------------------------------------------

    do k = kmscu, 1, - 1
        do i = 1, im
            if (cnvflg (i) .and. (k >= mrad (i) .and. k < krad (i))) then
                xmfd (i, k) = scaldfunc (i) * xmfd (i, k)
                dz = zl (i, k + 1) - zl (i, k)
                xmmx = dz / dt2
                xmfd (i, k) = min (xmfd (i, k), xmmx)
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! compute downdraft property using updated entranment rate
    ! -----------------------------------------------------------------------

    do i = 1, im
        if (cnvflg (i)) then
            k = krad (i)
            thld (i, k) = thlx (i, k)
        endif
    enddo

    ! do i = 1, im
    ! if (cnvflg (i)) then
    ! k = krad (i)
    ! ptem1 = max (qcdo (i, k, ntcw), 0.)
    ! tld = thld (i, k) / pix (i, k)
    ! tcdo (i, k) = tld + elocp * ptem1
    ! qcdo (i, k, 1) = qcdo (i, k, 1) + 0.2 * qcdo (i, k, 1)
    ! qcdo (i, k, ntcw) = qcdo (i, k, ntcw) + 0.2 * qcdo (i, k, ntcw)
    ! endif
    ! enddo

    do k = kmscu, 1, - 1
        do i = 1, im
            if (cnvflg (i) .and. (k >= mrad (i) .and. k < krad (i))) then
                dz = zl (i, k + 1) - zl (i, k)
                tem = 0.5 * xlamde (i, k) * dz
                factor = 1. + tem

                thld (i, k) = ((1. - tem) * thld (i, k + 1) + tem * &
                     (thlx (i, k) + thlx (i, k + 1))) / factor
                qtd (i, k) = ((1. - tem) * qtd (i, k + 1) + tem * &
                     (qtx (i, k) + qtx (i, k + 1))) / factor

                tld = thld (i, k) / pix (i, k)
                es = 0.01 * mqs (tld) ! mqs in pa
                qs = max (qmin, eps * es / (plyr (i, k) + epsm1 * es))
                dq = qtd (i, k) - qs

                if (dq > 0.) then
                    gamma = el2orc * qs / (tld ** 2)
                    qld = dq / (1. + gamma)
                    qtd (i, k) = qs + qld
                    qcdo (i, k, 1) = qs
                    qcdo (i, k, ntcw) = qld
                    tcdo (i, k) = tld + elocp * qld
                else
                    qcdo (i, k, 1) = qtd (i, k)
                    qcdo (i, k, ntcw) = 0.
                    tcdo (i, k) = tld
                endif

            endif
        enddo
    enddo

    do k = kmscu, 1, - 1
        do i = 1, im
            if (cnvflg (i) .and. k < krad (i)) then
                if (k >= mrad (i)) then
                    dz = zl (i, k + 1) - zl (i, k)
                    tem = 0.5 * xlamdem (i, k) * dz
                    factor = 1. + tem
                    ptem = tem - pgcon
                    ptem1 = tem + pgcon

                    ucdo (i, k) = ((1. - tem) * ucdo (i, k + 1) + ptem * u1 (i, k + 1) &
                         + ptem1 * u1 (i, k)) / factor
                    vcdo (i, k) = ((1. - tem) * vcdo (i, k + 1) + ptem * v1 (i, k + 1) &
                         + ptem1 * v1 (i, k)) / factor
                endif
            endif
        enddo
    enddo

    if (ntcw > 2) then

        do n = 2, ntcw - 1
            do k = kmscu, 1, - 1
                do i = 1, im
                    if (cnvflg (i) .and. k < krad (i)) then
                        if (k >= mrad (i)) then
                            dz = zl (i, k + 1) - zl (i, k)
                            tem = 0.5 * xlamde (i, k) * dz
                            factor = 1. + tem

                            qcdo (i, k, n) = ((1. - tem) * qcdo (i, k + 1, n) + tem * &
                                (q1 (i, k, n) + q1 (i, k + 1, n))) / factor
                        endif
                    endif
                enddo
            enddo
        enddo

    endif

    ndc = ntrac1 - ntcw

    if (ndc > 0) then

        do n = ntcw + 1, ntrac1
            do k = kmscu, 1, - 1
                do i = 1, im
                    if (cnvflg (i) .and. k < krad (i)) then
                        if (k >= mrad (i)) then
                            dz = zl (i, k + 1) - zl (i, k)
                            tem = 0.5 * xlamde (i, k) * dz
                            factor = 1. + tem

                            qcdo (i, k, n) = ((1. - tem) * qcdo (i, k + 1, n) + tem * &
                                (q1 (i, k, n) + q1 (i, k + 1, n))) / factor
                        endif
                    endif
                enddo
            enddo
        enddo

    endif

    return

end subroutine mfscu

! =======================================================================
! routine to solve the tridiagonal system to calculate temperature and
! moisture at \f$ t + \delta t \f$; part of two - part process to
! calculate time tendencies due to vertical diffusion.
! =======================================================================

subroutine tridi2 (l, n, cl, cm, cu, r1, r2, au, a1, a2)

    implicit none

    integer :: k, n, l, i

    real :: fk

    real :: cl (l, 2:n), cm (l, n), cu (l, n - 1), r1 (l, n), r2 (l, n), &
        au (l, n - 1), a1 (l, n), a2 (l, n)

    do i = 1, l
        fk = 1. / cm (i, 1)
        au (i, 1) = fk * cu (i, 1)
        a1 (i, 1) = fk * r1 (i, 1)
        a2 (i, 1) = fk * r2 (i, 1)
    enddo
    do k = 2, n - 1
        do i = 1, l
            fk = 1. / (cm (i, k) - cl (i, k) * au (i, k - 1))
            au (i, k) = fk * cu (i, k)
            a1 (i, k) = fk * (r1 (i, k) - cl (i, k) * a1 (i, k - 1))
            a2 (i, k) = fk * (r2 (i, k) - cl (i, k) * a2 (i, k - 1))
        enddo
    enddo
    do i = 1, l
        fk = 1. / (cm (i, n) - cl (i, n) * au (i, n - 1))
        a1 (i, n) = fk * (r1 (i, n) - cl (i, n) * a1 (i, n - 1))
        a2 (i, n) = fk * (r2 (i, n) - cl (i, n) * a2 (i, n - 1))
    enddo
    do k = n - 1, 1, - 1
        do i = 1, l
            a1 (i, k) = a1 (i, k) - au (i, k) * a1 (i, k + 1)
            a2 (i, k) = a2 (i, k) - au (i, k) * a2 (i, k + 1)
        enddo
    enddo

end subroutine tridi2

! =======================================================================
! routine to solve the tridiagonal system to calculate u - and v -
! momentum at \f$ t + \delta t \f$; part of two - part process to
! calculate time tendencies due to vertical diffusion.
! =======================================================================

subroutine tridin (l, n, nt, cl, cm, cu, r1, r2, au, a1, a2)

    implicit none

    integer :: is, k, kk, n, nt, l, i

    real :: fk (l)

    real :: cl (l, 2:n), cm (l, n), cu (l, n - 1), &
        r1 (l, n), r2 (l, n * nt), &
        au (l, n - 1), a1 (l, n), a2 (l, n * nt), &
        fkk (l, 2:n - 1)

    do i = 1, l
        fk (i) = 1. / cm (i, 1)
        au (i, 1) = fk (i) * cu (i, 1)
        a1 (i, 1) = fk (i) * r1 (i, 1)
    enddo
    do k = 1, nt
        is = (k - 1) * n
        do i = 1, l
            a2 (i, 1 + is) = fk (i) * r2 (i, 1 + is)
        enddo
    enddo
    do k = 2, n - 1
        do i = 1, l
            fkk (i, k) = 1. / (cm (i, k) - cl (i, k) * au (i, k - 1))
            au (i, k) = fkk (i, k) * cu (i, k)
            a1 (i, k) = fkk (i, k) * (r1 (i, k) - cl (i, k) * a1 (i, k - 1))
        enddo
    enddo
    do kk = 1, nt
        is = (kk - 1) * n
        do k = 2, n - 1
            do i = 1, l
                a2 (i, k + is) = fkk (i, k) * (r2 (i, k + is) - cl (i, k) * a2 (i, k + is - 1))
            enddo
        enddo
    enddo
    do i = 1, l
        fk (i) = 1. / (cm (i, n) - cl (i, n) * au (i, n - 1))
        a1 (i, n) = fk (i) * (r1 (i, n) - cl (i, n) * a1 (i, n - 1))
    enddo
    do k = 1, nt
        is = (k - 1) * n
        do i = 1, l
            a2 (i, n + is) = fk (i) * (r2 (i, n + is) - cl (i, n) * a2 (i, n + is - 1))
        enddo
    enddo
    do k = n - 1, 1, - 1
        do i = 1, l
            a1 (i, k) = a1 (i, k) - au (i, k) * a1 (i, k + 1)
        enddo
    enddo
    do kk = 1, nt
        is = (kk - 1) * n
        do k = n - 1, 1, - 1
            do i = 1, l
                a2 (i, k + is) = a2 (i, k + is) - au (i, k) * a2 (i, k + is + 1)
            enddo
        enddo
    enddo

end subroutine tridin

end module sa_tke_edmf_mod
