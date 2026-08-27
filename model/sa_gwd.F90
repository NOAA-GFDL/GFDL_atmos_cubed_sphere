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
! Scale-Aware Gravity Wave Drag (SA-GWD) Package
! This package includes orographic gravity wave drag, mountain blokcing,
! and convective gravity wave drag
! Developer:
! References:
! =======================================================================

module sa_gwd_mod

    use fms_mod, only: check_nml_error

    implicit none

    private

    ! -----------------------------------------------------------------------
    ! public subroutines, functions, and variables
    ! -----------------------------------------------------------------------

    public :: sa_gwd_init
    public :: sa_gwd_oro
    public :: sa_gwd_cnv

    ! -----------------------------------------------------------------------
    ! physics constants
    ! -----------------------------------------------------------------------

    real, parameter :: grav = 9.80665 ! acceleration due to gravity (m/s^2), ref: IFS

    real, parameter :: rerth = 6.3712e6 ! radius of earth (m)

    real, parameter :: pi = 4.0 * atan (1.0) ! ratio of circle circumference to diameter

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

    ! -----------------------------------------------------------------------
    ! namelist parameters
    ! -----------------------------------------------------------------------

    integer :: nmtvr = 14 ! number of topographic variables such as variance etc

    real :: cdmbgwd (2) = (/2.0, 0.25/) ! multiplication factors for cdmb and gwd
    real :: p_crit = 0.                 ! Optional level above which GWD stress decays with height
    real :: cgwf (2) = (/0.5, 0.05/)    !< multiplication factor for convective GWD

    ! -----------------------------------------------------------------------
    ! namelist
    ! -----------------------------------------------------------------------

    namelist / sa_gwd_nml / nmtvr, cdmbgwd, p_crit, cgwf

contains

! =======================================================================
! GWD initialization
! =======================================================================

subroutine sa_gwd_init (input_nml_file, logunit)

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

    read (input_nml_file, nml = sa_gwd_nml, iostat = ios)
    ierr = check_nml_error (ios, 'sa_gwd_nml')

    ! -----------------------------------------------------------------------
    ! write namelist to log file
    ! -----------------------------------------------------------------------

    write (logunit, *) " ================================================================== "
    write (logunit, *) "gwd_mod"
    write (logunit, nml = sa_gwd_nml)

end subroutine sa_gwd_init

! =======================================================================
! This subroutine is the parameterization of orographic gravity wave
! drag and mountain blocking.
!
! At present, global models must be run with horizontal resolutions
! that cannot typically resolve atmospheric phenomena shorter than
! ~10-100 km or greater for weather prediction and ~100-1000 km or
! greater for climate predicition. Many atmospheric processes have
! shorter horizontal scales than these "subgrid-scale" processes
! interact with and affect the larger-scale atmosphere in important
! ways.
!
! Atmospheric gravity waves are one such unresolved processes. These
! waves are generated by lower atmospheric sources. E.g., flow over
! irregularities at the Earth's surface such as mountains and valleys,
! uneven distribution of diabatic heat sources asscociated with
! convective systems, and highly dynamic atmospheric processes such
! as jet streams and fronts. The dissipation of these waves produces
! synoptic-scale body forces on the atmospheric flow, known as
! "gravity wave drag" (GWD), which affects both short-term evolution
! of weather systems and long-term climate. However, the spatial
! scales of these waves (in the range of ~5-500 km horizontally) are
! too short to be fully captured in models, and so GWD must be
! parameterized. In addition, the role of GWD in driving the global
! middle atmosphere circulation and thus global mean wind/temperature
! structures is well established. Thus, GWD parametrizations are now
! critical components of virtually all large-scale atmospheric models.
! GFS physics includes parameterizations of gravity waves from two
! important sources: mountains and convection.
!
! Atmospheric flow is significantly influenced by orography creating
! lift and frictional forces. The representation of orography and its
! influence in numerical weather prediction models are necessarily
! divided into the resolvable scales of motion and treated by
! primitive equations, the remaining sub-grid scales to be treated by
! parameterization. In terms of large scale NWP models, mountain
! blocking of wind flow around sub-grid scale orograph is a process
! that retards motion at various model vertical levels near or in the
! boundary layer. Flow around the mountain encounters larger
! frictional forces by being in contact with the mountain surfaces
! for longer time as well as the interaction of the atmospheric
! environment with vortex shedding which occurs in numerous
! observations. Lott and Miller (1997) \cite lott_and_miller_1997,
! incorporated the dividing streamline and mountain blocking in
! conjunction with sub-grid scale vertically propagating gravity wave
! parameterization in the context of NWP. The dividing streamline is
! seen as a source of gravity waves to the atmosphere above and
! nonlinear subgrid low-level mountain drag effect below.
!
! In a review paper on gravity waves in the middle atmosphere, Fritts
! (1984) \cite fritts_1984 showed that a large portion of observed
! gravity wave momentum flux has higher frequencies than those of
! stationary mountain waves. This phenomenon was explained by cumulus
! convection, which is an additional source of tropospheric gravity
! waves, and is particularly important in summertime. When the surface
! wind and stability are weak, the magnitude of the surface drag and
! the resultant influence of orographically-induced gravity wave drag
! on the large-scale flow are relatively small compared with those in
! wintertime (Palmer et al. 1986 \cite palmer_et_al_1986). In this
! situation, the relative importance of cumulus convection as a source
! of gravity waves is larger. in addition, in the tropical regions
! where persistent convection exists, deep cumulus clouds impinging on
! the stable stratosphere can generate gravity waves that influence
! the large-scale flow.
!
!-----------------------------------------------------------------------
! Outlines GWD parameterization in GFS:
!
! - Gravity-wave drag is simulated as described by Alpert et al.
! (1988) \cite alpert_et_al_1988. The parameterization includes
! determination of the momentum flux due to gravity waves at the
! surface, as well as upper levels. The surface stress is a nonlinear
! function of the surface wind speed and the local froude number,
! following Pierrehumbert (1987) \cite pierrehumbert_1987. Vertical
! variations in the momentum flux occur when the local richardson
! number is less than 0.25 (the stress vanishes), or when wave
! breaking occurs (local froude number becomes critical); in the
! latter case, the momentum flux is reduced according to the
! Lindzen (1981) \cite lindzen_1981 wave saturation hypothesis.
! Modifications are made to avoid instability when the critical layer
! is near the surface, since the time scale for gravity-wave drag is
! shorter than the model time step.
!
! - The treatment of the GWD in the lower troposphere is enhanced
! according to Kim and Arakawa (1995) \cite kim_and_arakawa_1995.
! Orographic Std Dev (HPRIME), Convexity (OC), Asymmetry (OA4) and Lx
! (CLX4) are input topographic statistics needed (see Appendix in Kim
! and Arakawa (1995) \cite kim_and_arakawa_1995).
!
! - Mountain blocking influences are incorporated following the Lott
! and Miller (1997) \cite lott_and_miller_1997 parameterization with
! minor changes, including their dividing streamline concept. The
! model subgrid scale orography is represented by four parameters,
! after Baines and Palmer (1990) \cite baines_and_palmer_1990, the
! standard deviation (HPRIME), the anisotropy (GAMMA), the slope
! (SIGMA) and the geographical orientation of the orography (THETA).
! These are calculated off-line as a function of model resolution in
! the fortran code ml01rg2.f, with script mlb2.sh (see Appendix:
! specification of subgrid-scale orography in Lott and Miller (1997)
! \cite lott_and_miller_1997).
!
! - The orographic GWD parameterizations automatically scales
! with model resolution. For example, the T574L64 version of GFS uses
! four times stronger mountain blocking and one half the strength of
! gravity wave drag than the T383L64 version.
!
! - The parameterization of stationary convectively-forced GWD follows
! the development of Chun and Baik (1998) \cite chun_and_baik_1998,
! which was tested in GCMs by Chun et al. (2001, 2004)
! \cite chun_et_al_2001 \cite chun_et_al_2004 was implemented in GFS
! by Ake Johansson (2008) and the work of the GCWMB staff. Modest
! positive effects from using the parameterization are seen in the
! tropical upper troposphere and lower stratosphere.
!
!-----------------------------------------------------------------------
! intra_gwdps Intraphysics Communication
!
! - Routine gwdps (\ref orographic) is called from gbphys after call
! to moninedmf
! - Routine gwdc (\ref convective) is called from gbphys after call
! to sascnvn
!
! The time tendencies of zonal and meridional wind are altered to
! include the effect of mountain induced gravity wave drag from
! subgrid scale orography including convective breaking, shear
! breaking and the presence of critical levels.
!
!-----------------------------------------------------------------------
! \param[in] IM      horizontal number of used pts
! \param[in] KM      vertical layer dimension
! \param[in, out] VTGWD  non-linear tendency for v wind component
! \param[in, out] UTGWD  non-linear tendency for u wind component
! \param[in, out] TTGWD  non-linear tendency for temperature (not used)
! \param[in] U1      zonal wind component of model layer wind (m/s)
! \param[in] V1      meridional wind component of model layer wind (m/s)
! \param[in] T1      model layer mean temperature (K)
! \param[in] Q1      model layer mean specific humidity
! \param[in] KPBL    index for the PBL top layer
! \param[in] PRSI    pressure at layer interfaces
! \param[in] DEL     positive increment of p/psfc across layer
! \param[in] PRSL    mean layer pressure
! \param[in] PRSLK   exner function at layer
! \param[in] PHII    interface geopotential (\f$m^2/s^2\f$)
! \param[in] PHIL    layer geopotential (\f$m^2/s^2\f$)
! \param[in] DELT    physics time step in seconds
! \param[in] HPRIME  orographic standard deviation (m) (mtnvar (:, 1))
! \param[in] OC      orographic convexity (mtnvar (:, 2))
! \param[in] OA4     orographic asymmetry (mtnvar (:, 3:6))
! \param[in] CLX4    Lx, the fractional area covered by the subgrid-scale
!                    orography higher than a critical height for a grid box
!                    with the interval \f$ \triangle x \f$ (mtnvar (:, 7:10))
! \param[in] THETA   the angle of the mtn with that to the east (x) axis (mtnvar (:, 11))
! \param[in] SIGMA   orographic slope (mtnvar (:, 13))
! \param[in] GAMMA   orographic anisotropy (mtnvar (:, 12))
! \param[in] ELVMAX  orographic maximum (mtnvar (:, 14))
! \param[out] DUSFC  u component of surface stress
! \param[out] DVSFC  v component of surface stress
! \param[in] IMX     number of longitude points
! \param[in] NMTVR   number of topographic variables such as variance etc
!                    used in the GWD parameterization, current operational, nmtvr = 14
! \param[in] CDMBGWD multiplication factors for cdmb and gwd
!
!-----------------------------------------------------------------------
! Implementation Versions
!-----------------------------------------------------------------------
! -- not in this code -- history of GWDP at NCEP --
!
! Version 3 modified for gravity waves, location: .fr30 (v3gwd) * j *
! 3.1 includes variable saturation flux profile cf isigst
! 3.g includes ps combined w / ph (GLAS and GFDL)
!     also included is ri smooth over a thick lower layer
!     also included is decrease in de-acc at top by 1 / 2
!     the NMC GWD incorporating both GLAS (P&S) and GFDL (MIGWD)
!     mountain induced gravity wave drag
!     code from .fr30 (v3monnx) for monin3
!     this version (06 Mar 1987)
!     this version (26 Apr 1987) 3.g
!     this version (01 May 1987) 3.9
!     change to fortran 77 (Feb 1989) --- Hann-Ming Henry Juang
!     20070601 elvmax bug fix (* j *)
!
!-----------------------------------------------------------------------
! Version 4 -- this code
!
! Modified to implement the enhanced low tropospheric gravity
! wave drag developed by Kim and Arakawa (JAS, 1995) .
!
! Orographic Std Dev (HPRIME), Convexity (OC), Asymmetry (OA4)
! and Lx (CLX4) are input topographic statistics needed.
!
! Programmed and debugged by Hong, Alpert and Kim --- Jan 1996.
! Debugged again by Moorthi and Iredell --- May 1998.
!
! Further cleanup, optimization and modification by Moorthi --- May 98, Mar 99.
!
! Modified for usgs orography data (NCEP office note 424)
! and with several bugs fixed by Moorthi and Hong --- Jul 1999.
!
! Modified & implemented into NRL NOGAPS by Young and Joon Kim --- Jul 2000.
!
!-----------------------------------------------------------------------
! Version LM MB (6) : oz fix 8 / 2003 -- this code
!
! Changed to include the Lott and Miller MTN blocking
! with some modifications by (* j *) 4 / 02
! from a principal coordinate calculation using the
! hi res 8 minute orography, the angle of the
! mtn with that to the east (x) axis is theta, the slope
! parameter sigma. the anisotropy is in gamma - all are input
! topographic statistics needed. these are calculated off - line
! as a function of model resolution in the fortran code ml01rg2.f,
! with script mlb2.sh. (* j *)
!
! gwdps_mb.f version (following lmi) elvmax < hncrit (* j *)
! mb3a expt to enhance elvmax mtn hgt see sigfac & hncrit
! gwdps_gwdfix_v6.f fixgwd gf6.0 20070608 sigfac = 4.
!-----------------------------------------------------------------------
!
! USE
! routine is called from gbphys (after call to monnin)
!
! PURPOSE
! using the gwd parameterizations of ps - GLAS and ph -
! GFDL technique. the time tendencies of u v
! are altered to include the effect of mountain induced
! gravity wave drag from sub - grid scale orography including
! convective breaking, shear breaking and the presence of
! critical levels
!
! INPUT
! VTGWD (im, km) non-lin tendency for v wind component
! UTGWD (im, km) non-lin tendency for u wind component
! TTGWD (im, km) non-lin tendency for temperature
! U1 (im, km) zonal wind m / sec at t0 - dt
! V1 (im, km) meridional wind m / sec at t0 - dt
! T1 (im, km) temperature deg k at t0 - dt
! Q1 (im, km) specific humidity at t0 - dt
!
! delt time step secs
! SI (n) p / psfc at base of layer n
! SL (n) p / psfc at middle of layer n
! DEL (n) positive increment of p / psfc across layer n
! KPBL (im) is the index of the top layer of the pbl
!
! OUTPUT
! VTGWD, UTGWD as augmented by tendency due to gwdps
! other input variables unmodified.
!
!-----------------------------------------------------------------------
! REVISION LOG:
! May 2013 J. Wang change cleff back to opn setting
! Jan 2014 J. Wang merge henry and fangin's dissipation heat in GFS to nems
! =======================================================================

subroutine sa_gwd_oro (im, km, u1, v1, t1, q1, delt, gsize, &
        kpbl, prsi, del, prsl, prslk, phii, phil, &
        hprime, oc, oa4, clx4, theta, sigma, gamma, elvmax, &
        utgwd, vtgwd, ttgwd, dusfc, dvsfc, rdxzb)

    implicit none

    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------

    integer, intent (in) :: im, km
    integer, intent (in) :: kpbl (im) ! index for the pbl top layer

    real, intent (in) :: delt
    real, intent (in) :: gsize (im)
    real, intent (in) :: prsl (im, km), prsi (im, km + 1), del (im, km), &
        prslk (im, km), phil (im, km), phii (im, km + 1)
    real, intent (in) :: oc (im), oa4 (im, 4), clx4 (im, 4), hprime (im)
    real, intent (in) :: elvmax (im), theta (im), sigma (im), gamma (im)

    real, intent (inout) :: u1 (im, km), v1 (im, km), t1 (im, km), q1 (im, km)

    real, intent (out), optional :: utgwd (im, km), vtgwd (im, km), ttgwd (im, km)

    real, intent (out), optional :: dusfc (im), dvsfc (im), rdxzb (im)

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    ! for lm mtn blocking
    real :: wk (im), emax (im)
    real :: bnv2lm (im, km), pe (im), ek (im), zbk (im), up (im)
    real :: db (im, km), ang (im, km), uds (im, km)
    real :: zlen, dbtmp, r, phiang, cdmb (im), dbim
    real :: eng0, eng1

    ! some constants

    real :: dw2min, rimin, ric, bnv2min, efmin, &
        efmax, hpmax, hpmin, rad_to_deg, deg_to_rad
    parameter (rad_to_deg = 180.0 / pi, deg_to_rad = pi / 180.0)
    parameter (dw2min = 1., rimin = - 100., ric = 0.25, bnv2min = 1.0e-5)
    ! parameter (efmin = 0.0, efmax = 10.0, hpmax = 200.0)
    parameter (efmin = 0.0, efmax = 10.0, hpmax = 2400.0, hpmin = 1.0)
    ! parameter (p_crit = 30.e2)

    real :: frc, ce, ceofrc, frmax, cg, gmax, &
        veleps, factop, rlolev, rdi
    ! critac, veleps, factop, rlolev, rdi
    parameter (frc = 1.0, ce = 0.8, ceofrc = ce / frc, frmax = 100., cg = 0.5)
    parameter (gmax = 1.0, veleps = 1.0, factop = 0.5)
    ! parameter (gmax = 1.0, critac = 5.0e-4, veleps = 1.0, factop = 0.5)
    parameter (rlolev = 50000.0)
    ! parameter (rlolev = 500.0)
    ! parameter (rlolev = 0.5)

    real :: dpmin, hminmt, hncrit, minwnd, sigfac
    ! --- for lm mtn blocking
    ! parameter (cdmb = 1.0) ! non - dim sub grid mtn drag amp (* j *)
    parameter (hncrit = 8000.) ! max value in meters for emax (* j *)
    ! hncrit set to 8000m and sigfac added to enhance emax mtn hgt
    parameter (sigfac = 4.0) ! mb3a expt test for emax factor (* j *)
    parameter (hminmt = 50.) ! min mtn height (* j *)
    parameter (minwnd = 0.1) ! min wind component (* j *)

    ! parameter (dpmin = 00.0) ! minimum thickness of the reference layer
    !! parameter (dpmin = 05.0) ! minimum thickness of the reference layer
    ! parameter (dpmin = 20.0) ! minimum thickness of the reference layer
    ! in centibars
    parameter (dpmin = 5000.0) ! minimum thickness of the reference layer
    ! in pa

    real :: fdir
    integer :: mdir
    parameter (mdir = 8, fdir = mdir / (pi + pi))
    integer :: nwdir (mdir)
    data nwdir / 6, 7, 5, 8, 2, 3, 1, 4 /
    save nwdir

    logical :: icrilv (im)

    ! ---- mountain induced gravity wave drag

    real :: taub (im), xn (im), yn (im), ubar (im), &
        vbar (im), ulow (im), oa (im), clx (im), &
        roll (im), uloi (im), &
        dtfac (im), xlinv (im), delks (im), delks1 (im)

    real :: bnv2 (im, km), taup (im, km + 1), ri_n (im, km), &
        taud (im, km), ro (im, km), vtk (im, km), &
        vtj (im, km), scor (im), velco (im, km - 1), &
        bnv2bar (im)

    ! real :: velko (km - 1)
    integer :: kref (im), kint (im), iwk (im), ipt (im)
    ! for lm mtn blocking
    integer :: kreflm (im), iwklm (im)
    integer :: idxzb (im), ktrial, klevm1

    real :: gor, gocp, gr2, bnv, fr, &
        brvf, cleff (im), tem, tem1, tem2, temc, temv, &
        wdir, ti, rdz, dw2, shr2, bvf2, &
        rdelks, efact, coefm, gfobnv, &
        scork, rscor, hd, fro, rim, sira, &
        dtaux, dtauy, pkp1log, pklog
    integer :: kmm1, kmm2, lcap, lcapp1, kbps, kbpsp1, kbpsm1, &
        kmps, idir, nwd, i, j, k, klcap, kp1, kmpbl, npt, &
        kmll
    ! kmll, kmds, ihit, jhit

    ! parameter (cdmb = 1.0) ! non - dim sub grid mtn drag amp (* j *)
    ! non - dim sub grid mtn drag amp (* j *)
    ! cdmb = 1.0 / float (imx / 192)
    ! cdmb = 192.0 / float (imx)
    ! cdmb = 4.0 * 192.0 / float (imx)
    do i = 1, im
        cdmb (i) = 2.e-5 * gsize (i)
        if (cdmbgwd (1) >= 0.0) cdmb (i) = cdmb (i) * cdmbgwd (1)
    enddo

    do i = 1, im
        emax (i) = elvmax (i)
    enddo

    do k = 1, km
        do i = 1, im
            if (present (utgwd)) utgwd (i, k) = 0.
            if (present (vtgwd)) vtgwd (i, k) = 0.
            if (present (ttgwd)) ttgwd (i, k) = 0.
        enddo
    enddo

    do i = 1, im
        if (present (dusfc)) dusfc (i) = 0.
        if (present (dvsfc)) dvsfc (i) = 0.
    enddo

    do k = 1, km
        do i = 1, im
            db (i, k) = 0.
            ang (i, k) = 0.
            uds (i, k) = 0.
        enddo
    enddo

    rdi = 1.0 / rdgas
    gor = grav / rdgas
    gr2 = grav * gor
    gocp = grav / cp_air

    ! ncnt = 0
    kmm1 = km - 1
    kmm2 = km - 2
    lcap = km
    lcapp1 = lcap + 1


    if (nmtvr .eq. 14) then
        ! ---- for lm and gwd calculation points
        if (present (rdxzb)) rdxzb (:) = 0.
        ipt = 0
        npt = 0
        do i = 1, im
            if ((emax (i) .gt. hminmt) .and. (hprime (i) .gt. hpmin)) then
                npt = npt + 1
                ipt (npt) = i
            endif
        enddo
        if (npt .eq. 0) return ! no gwd / mb calculation done

        ! if (lprnt) print *, ' npt = ', npt, ' npr = ', npr, ' ipr = ', ipr, ' im = ', im, &
        ! ' ipt (npt) = ', ipt (npt)

        ! --- iwklm is the level above the height of the of the mountain.
        ! --- idxzb is the level of the dividing streamline.
        ! initialize dividing streamline (ds) control vector

        do i = 1, npt
            iwklm (i) = 2
            idxzb (i) = 0
            kreflm (i) = 0
        enddo
        ! if (lprnt) print *, ' in gwdps_lm.f npt, im, km, me = ', npt, im, km, me


        ! > --- subgrid mountain blocking section

        !..............................
        !..............................

        ! (* j *) 11 / 03: test upper limit on kmll = km - 1
        ! then do not need hncrit -- test with large hncrit first.
        ! kmll = km / 2 ! maximum mtnlm height : # of vertical levels / 2
        kmll = kmm1
        ! --- no mtn should be as high as kmll (so we do not have to start at
        ! --- the top of the model but could do calc for all levels) .

        do i = 1, npt
            j = ipt (i)
            emax (j) = min (emax (j) + sigfac * hprime (j), hncrit)
        enddo

        do k = 1, kmll
            do i = 1, npt
                j = ipt (i)
                ! --- interpolate to max mtn height for index, iwklm (i) wk[gz]
                ! --- emax is limited to hncrit because to hi res topo30 orog.
                pkp1log = phil (j, k + 1) / grav
                pklog = phil (j, k) / grav
                !!! ------- emax (j) = min (emax (j) + sigfac * hprime (j), hncrit)
                if ((emax (j) .le. pkp1log) .and. (emax (j) .ge. pklog)) then
                    ! print *, ' in gwdps_lm.f 1 = ', k, emax (j), pklog, pkp1log, me
                    ! --- wk for diags but can be saved and reused.
                    wk (i) = grav * emax (j) / (phil (j, k + 1) - phil (j, k))
                    iwklm (i) = max (iwklm (i), k + 1)
                    ! print *, ' in gwdps_lm.f 2 npt = ', npt, i, j, wk (i), iwklm (i), me
                endif

                ! --- find at prsl levels large scale environment variables
                ! --- these cover all possible mtn max heights
                vtj (i, k) = t1 (j, k) * (1. + zvir * q1 (j, k))
                vtk (i, k) = vtj (i, k) / prslk (j, k)
                ro (i, k) = rdi * prsl (j, k) / vtj (i, k) ! density kg / m ** 3
            enddo
        enddo

        ! testing for highest model level of mountain top

        ! ihit = 2
        ! jhit = 0
        ! do i = 1, npt
        ! j = ipt (i)
        ! if (iwklm (i) .gt. ihit) then
        ! ihit = iwklm (i)
        ! jhit = j
        ! endif
        ! enddo
        ! print *, ' mb: kdt, max (iwklm), jhit, phil, me = ', &
        ! kdt, ihit, jhit, phil (jhit, ihit), me

        klevm1 = kmll - 1
        do k = 1, klevm1
            do i = 1, npt
                j = ipt (i)
                rdz = grav / (phil (j, k + 1) - phil (j, k))
                ! --- brunt - vaisala frequency
                ! > - compute brunt - vaisala frequency \f$n\f$.
                bnv2lm (i, k) = (grav + grav) * rdz * (vtk (i, k + 1) - vtk (i, k)) &
                     / (vtk (i, k + 1) + vtk (i, k))
                bnv2lm (i, k) = max (bnv2lm (i, k), bnv2min)
            enddo
        enddo
        ! print *, ' in gwdps_lm.f 3 npt = ', npt, j, rdz, me

        do i = 1, npt
            j = ipt (i)
            delks (i) = 1.0 / (prsi (j, 1) - prsi (j, iwklm (i)))
            delks1 (i) = 1.0 / (prsl (j, 1) - prsl (j, iwklm (i)))
            ubar (i) = 0.0
            vbar (i) = 0.0
            roll (i) = 0.0
            pe (i) = 0.0
            ek (i) = 0.0
            bnv2bar (i) = (prsl (j, 1) - prsl (j, 2)) * delks1 (i) * bnv2lm (i, 1)
        enddo

        ! --- find the dividing stream line height
        ! --- starting from the level above the max mtn downward
        ! --- iwklm (i) is the k - index of mtn emax elevation
        ! > - find the dividing streamline height starting from the level above
        !! the maximum mountain height and processing downward.
        do ktrial = kmll, 1, - 1
            do i = 1, npt
                if (ktrial .lt. iwklm (i) .and. kreflm (i) .eq. 0) then
                    kreflm (i) = ktrial
                endif
            enddo
        enddo
        ! print *, ' in gwdps_lm.f 4 npt = ', npt, kreflm (npt), me

        ! --- in the layer kreflm (i) to 1 find pe (which needs n, emax)
        ! --- make averages, guess dividing stream (ds) line layer.
        ! --- this is not used in the first cut except for testing and
        ! --- is the vert ave of quantities from the surface to mtn top.

        do i = 1, npt
            do k = 1, kreflm (i)
                j = ipt (i)
                rdelks = del (j, k) * delks (i)
                ubar (i) = ubar (i) + rdelks * u1 (j, k) ! trial mean u below
                vbar (i) = vbar (i) + rdelks * v1 (j, k) ! trial mean v below
                roll (i) = roll (i) + rdelks * ro (i, k) ! trial mean ro below
                rdelks = (prsl (j, k) - prsl (j, k + 1)) * delks1 (i)
                bnv2bar (i) = bnv2bar (i) + bnv2lm (i, k) * rdelks
                ! --- these vert ave are for diags, testing and gwd to follow (* j *) .
            enddo
        enddo
        ! print *, ' in gwdps_lm.f 5 = ', i, kreflm (npt), bnv2bar (npt), me

        ! --- integrate to get pe in the trial layer.
        ! --- need the first layer where pe > ek - as soon as
        ! --- idxzb is not 0 we have a hit and zb is found.

        do i = 1, npt
            j = ipt (i)
            do k = iwklm (i), 1, - 1
                phiang = atan2 (v1 (j, k), u1 (j, k)) * rad_to_deg
                ang (i, k) = (theta (j) - phiang)
                if (ang (i, k) .gt. 90.) ang (i, k) = ang (i, k) - 180.
                if (ang (i, k) .lt. - 90.) ang (i, k) = ang (i, k) + 180.
                ang (i, k) = ang (i, k) * deg_to_rad

                ! > - compute wind speed uds
                !!\f[
                !! uds = \max (\sqrt{u1^2 + v1^2}, minwnd)
                !!\f]
                !! where \f$ minwnd = 0.1 \f$, \f$u1\f$ and \f$v1\f$ are zonal and
                !! meridional wind components of model layer wind.
                uds (i, k) = &
                    max (sqrt (u1 (j, k) * u1 (j, k) + v1 (j, k) * v1 (j, k)), minwnd)
                ! --- test to see if we found zb previously
                if (idxzb (i) .eq. 0) then
                    pe (i) = pe (i) + bnv2lm (i, k) * &
                         (grav * emax (j) - phil (j, k)) * &
                         (phii (j, k + 1) - phii (j, k)) / (grav * grav)
                    ! --- ke
                    ! --- wind projected on the line perpendicular to mtn range, u (zb (k)) .
                    ! --- kenetic energy is at the layer zb
                    ! --- theta ranges from - + 90deg |_ to the mtn "largest topo variations"
                    up (i) = uds (i, k) * cos (ang (i, k))
                    ek (i) = 0.5 * up (i) * up (i)

                    ! --- dividing stream lime is found when pe = exceeds ek.
                    if (pe (i) .ge. ek (i)) then
                        idxzb (i) = k
                        if (present (rdxzb)) rdxzb (j) = real (k)
                    endif
                    ! --- then mtn blocked flow is between zb = k (idxzb (i)) and surface

                    ! > - the dividing streamline height (idxzb), of a subgrid scale
                    !! obstable, is found by comparing the potential (pe) and kinetic
                    !! energies (ek) of the upstream large scale wind and subgrid scale air
                    !! parcel movements. the dividing streamline is found when
                    !! \f$pe\geq ek\f$. mountain - blocked flow is defined to exist between
                    !! the surface and the dividing streamline height (\f$h_d\f$), which
                    !! can be found by solving an integral equation for \f$h_d\f$:
                    !!\f[
                    !! \frac{u^{2} (h_{d}) }{2} = \int_{h_{d}}^{h} n^{2} (z) (h - z) dz
                    !!\f]
                    !! where \f$h\f$ is the maximum subgrid scale elevation within the grid
                    !! box of actual orography, \f$h\f$, obtained from the gtopo30 dataset
                    !! from the u.s. geological survey.
                endif
            enddo
        enddo

        ! print *, ' in gwdps_lm.f 6 = ', phiang, theta (ipt (npt)), me
        ! print *, ' in gwdps_lm.f 7 = ', idxzb (npt), pe (npt)

        ! if (lprnt .and. npr .gt. 0) then
        ! print *, ' bnv2bar, bnv2lm = ', bnv2bar (npr), bnv2lm (npr, 1:klevm1)
        ! print *, ' npr, idxzb, uds = ', npr, idxzb (npr), uds (npr, :)
        ! print *, ' pe, up, ek = ', pe (npr), up (npr), ek (npr)
        ! endif

        do i = 1, npt
            j = ipt (i)
            ! --- calc if n constant in layers (zb guess) - a diagnostic only.
            zbk (i) = emax (j) &
                 - sqrt (ubar (i) * ubar (i) + vbar (i) * vbar (i)) / bnv2bar (i)
        enddo

        ! if (lprnt .and. npr .gt. 0) then
        ! print *, ' iwklm, zbk = ', iwklm (npr), zbk (npr), idxzb (npr)
        ! print *, ' zb = ', phil (ipr), idxzb (npr)) / grav
        ! print *, ' in gwdps_lm.f 8 npt = ', npt, zbk (npt), up (npt), me
        ! endif

        ! --- the drag for mtn blocked flow

        do i = 1, npt
            j = ipt (i)
            zlen = 0.
            ! print *, ' in gwdps_lm.f 9 = ', i, j, idxzb (i), me
            if (idxzb (i) .gt. 0) then
                do k = idxzb (i), 1, - 1
                    if (phil (j, idxzb (i)) .gt. phil (j, k)) then

                        ! > - calculate \f$zlen\f$, which sums up a number of contributions of
                        !! elliptic obstables.
                        !!\f[
                        !! zlen = \sqrt{[\frac{h_{d} - z}{z + h'}]}
                        !!\f]
                        !! where \f$z\f$ is the height, \f$h'\f$ is the orographic standard
                        !! deviation (hprime) .
                        zlen = sqrt ((phil (j, idxzb (i)) - phil (j, k)) / &
                             (phil (j, k) + grav * hprime (j)))
                        ! --- lm eq 14:
                        ! > - calculate the drag coefficient to vary with the aspect ratio of
                        !! the obstable as seen by the incident flow (see eq.14 in lott and
                        !! miller (1997) \cite lott_and_miller_1997)
                        !!\f[
                        !! r = \frac{\cos^{2}\psi + \gamma\sin^{2}\psi}{\gamma\cos^{2}\psi + \sin^{2}\psi}
                        !!\f]
                        !! where \f$\psi\f$, which is derived from theta, is the angle between
                        !! the incident flow direction and the normal ridge direcion.
                        !! \f$\gamma\f$ is the orographic anisotropy (gamma) .
                        r = cos (ang (i, k)) ** 2 + gamma (j) * sin (ang (i, k)) ** 2
                        if (abs (r) .lt. 1.e-20) then
                            db (i, k) = 0.0
                        else
                            r = (gamma (j) * cos (ang (i, k)) ** 2 + sin (ang (i, k)) ** 2) / r
                            ! --- (negitive of db -- see sign at tendency)
                            ! > - in each model layer below the dividing streamlines, a drag from
                            !! the blocked flow is exerted by the obstacle on the large scale flow.
                            !! the drag per unit area and per unit height is written (eq.15 in
                            !! lott and miller (1997) \cite lott_and_miller_1997) :
                            !!\f[
                            !! d_{b} (z) = - c_{d}\max (2 - \frac{1}{r}, 0) \rho\frac{\sigma}{2h'}zlen\max (\cos\psi, \gamma\sin\psi) \frac{uds}{2}
                            !!\f]
                            !! where \f$c_{d}\f$ is a specified constant, \f$\sigma\f$ is the
                            !! orographic slope.

                            dbtmp = 0.25 * cdmb (j) * &
                                max (2. - r, 0.) * sigma (j) * &
                                max (cos (ang (i, k)), gamma (j) * sin (ang (i, k))) * &
                                zlen / hprime (j)
                            db (i, k) = dbtmp * uds (i, k)
                        endif

                        ! if (lprnt .and. i .eq. npr) then
                        ! print *, ' in gwdps_lmi.f 10 npt = ', npt, i, j, idxzb (i), &
                        ! dbtmp, r' ang = ', ang (i, k), ' gamma = ', gamma (j), ' k = ', k
                        ! print *, ' in gwdps_lmi.f 11 k = ', k, zlen, cos (ang (i, k))
                        ! print *, ' in gwdps_lmi.f 12 db = ', db (i, k), sin (ang (i, k))
                        ! endif
                    endif
                enddo
                ! if (lprnt) print *, ' @k = 1, zlen, dbtmp = ', k, zlen, dbtmp
            endif
        enddo

        !.............................
        !.............................
        ! end mtn blocking section

    elseif (nmtvr .ne. 14) then
        ! ---- for mb not present and gwd (nmtvr .ne .14)
        ipt = 0
        npt = 0
        do i = 1, im
            if (hprime (i) .gt. hpmin) then
                npt = npt + 1
                ipt (npt) = i
            endif
        enddo
        if (npt .eq. 0) return ! no gwd / mb calculation done

        ! if (lprnt) print *, ' npr = ', npr, ' npt = ', npt, ' ipr = ', ipr, &
        ! ' ipt (npt) = ', ipt (npt)

        do i = 1, npt
            idxzb (i) = 0
            if (present (rdxzb)) rdxzb (i) = 0.
        enddo
    endif

    !.............................
    !.............................

    ! > --- orographic gravity wave drag section
    kmpbl = km / 2 ! maximum pbl height : # of vertical levels / 2

    ! scale cleff between im = 384 * 2 and 192 * 2 for t126 / t170 and t62

    !if (imx .gt. 0) then
    !cleff = 1.0e-5 * sqrt (float (imx) / 384.0) ! this is inverse of cleff
    !cleff = 1.0e-5 * sqrt (float (imx) / 192.0) ! this is inverse of cleff
    !cleff = 0.5e-5 * sqrt (float (imx) / 192.0) ! this is inverse of cleff
    !cleff = 1.0e-5 * sqrt (float (imx) / 192) / float (imx / 192)
    !cleff = 1.0e-5 / sqrt (float (imx) / 192.0) ! this is inverse of cleff
    !cleff = 0.5e-5 / sqrt (float (imx) / 192.0) ! this is inverse of cleff
    !hmhj for ndsl
    !jw cleff = 0.1e-5 / sqrt (float (imx) / 192.0) ! this is inverse of cleff
    !cleff = 2.0e-5 * sqrt (float (imx) / 192.0) ! this is inverse of cleff
    !cleff = 2.5e-5 * sqrt (float (imx) / 192.0) ! this is inverse of cleff
    !ndif
    do i = 1, im
        cleff (i) = 0.25e-5 * sqrt (2.e-5 * gsize (i))
        if (cdmbgwd (2) >= 0.0) cleff (i) = cleff (i) * cdmbgwd (2)
    enddo

    do k = 1, km
        do i = 1, npt
            j = ipt (i)
            vtj (i, k) = t1 (j, k) * (1. + zvir * q1 (j, k))
            vtk (i, k) = vtj (i, k) / prslk (j, k)
            ro (i, k) = rdi * prsl (j, k) / vtj (i, k) ! density tons / m ** 3
            taup (i, k) = 0.0
        enddo
    enddo
    do k = 1, kmm1
        do i = 1, npt
            j = ipt (i)
            ti = 2.0 / (t1 (j, k) + t1 (j, k + 1))
            tem = ti / (prsl (j, k) - prsl (j, k + 1))
            rdz = grav / (phil (j, k + 1) - phil (j, k))
            tem1 = u1 (j, k) - u1 (j, k + 1)
            tem2 = v1 (j, k) - v1 (j, k + 1)
            dw2 = tem1 * tem1 + tem2 * tem2
            shr2 = max (dw2, dw2min) * rdz * rdz
            bvf2 = grav * (gocp + rdz * (vtj (i, k + 1) - vtj (i, k))) * ti
            ri_n (i, k) = max (bvf2 / shr2, rimin) ! richardson number
            ! brunt - vaisala frequency
            ! tem = gr2 * (prsl (j, k) + prsl (j, k + 1)) * tem
            ! bnv2 (i, k) = tem * (vtk (i, k + 1) - vtk (i, k)) / (vtk (i, k + 1) + vtk (i, k))
            bnv2 (i, k) = (grav + grav) * rdz * (vtk (i, k + 1) - vtk (i, k)) &
                 / (vtk (i, k + 1) + vtk (i, k))
            bnv2 (i, k) = max (bnv2 (i, k), bnv2min)
        enddo
    enddo
    ! print *, ' in gwdps_lm.f gwd:14 = ', npt, kmm1, bnv2 (npt, kmm1)

    ! apply 3 point smoothing on bnv2

    ! do k = 1, km
    ! do i = 1, im
    ! vtk (i, k) = bnv2 (i, k)
    ! enddo
    ! enddo
    ! do k = 2, kmm1
    ! do i = 1, im
    ! bnv2 (i, k) = 0.25 * (vtk (i, k - 1) + vtk (i, k + 1)) + 0.5 * vtk (i, k)
    ! enddo
    ! enddo

    ! finding the first interface index above 50 hpa level

    do i = 1, npt
        iwk (i) = 2
    enddo
    do k = 3, kmpbl
        do i = 1, npt
            j = ipt (i)
            tem = (prsi (j, 1) - prsi (j, k))
            if (tem .lt. dpmin) iwk (i) = k
        enddo
    enddo

    ! > - calculate the reference level index: kref = max (2, kpbl + 1) . where
    !! kpbl is the index for the pbl top layer.
    kbps = 1
    kmps = km
    do i = 1, npt
        j = ipt (i)
        kref (i) = max (iwk (i), kpbl (j) + 1) ! reference level
        delks (i) = 1.0 / (prsi (j, 1) - prsi (j, kref (i)))
        delks1 (i) = 1.0 / (prsl (j, 1) - prsl (j, kref (i)))
        ubar (i) = 0.0
        vbar (i) = 0.0
        roll (i) = 0.0
        kbps = max (kbps, kref (i))
        kmps = min (kmps, kref (i))

        bnv2bar (i) = (prsl (j, 1) - prsl (j, 2)) * delks1 (i) * bnv2 (i, 1)
    enddo
    ! print *, ' in gwdps_lm.f gwd:15 = ', kbps, kmps
    kbpsp1 = kbps + 1
    kbpsm1 = kbps - 1
    do k = 1, kbps
        do i = 1, npt
            if (k .lt. kref (i)) then
                j = ipt (i)
                rdelks = del (j, k) * delks (i)
                ubar (i) = ubar (i) + rdelks * u1 (j, k) ! mean u below kref
                vbar (i) = vbar (i) + rdelks * v1 (j, k) ! mean v below kref

                roll (i) = roll (i) + rdelks * ro (i, k) ! mean ro below kref
                rdelks = (prsl (j, k) - prsl (j, k + 1)) * delks1 (i)
                bnv2bar (i) = bnv2bar (i) + bnv2 (i, k) * rdelks
            endif
        enddo
    enddo
    ! print *, ' in gwdps_lm.f gwd:15b = ', bnv2bar (npt)

    ! figure out low - level horizontal wind direction and find 'oa'

    ! nwd 1 2 3 4 5 6 7 8
    ! wd w s sw nw e n ne se

    ! > - calculate low - level horizontal wind direction, the derived
    !! orographic asymmetry parameter (oa), and the derived lx (clx) .
    do i = 1, npt
        j = ipt (i)
        wdir = atan2 (ubar (i), vbar (i)) + pi
        idir = mod (nint (fdir * wdir), mdir) + 1
        nwd = nwdir (idir)
        oa (i) = (1 - 2 * int ((nwd - 1) / 4)) * oa4 (j, mod (nwd - 1, 4) + 1)
        clx (i) = clx4 (j, mod (nwd - 1, 4) + 1)
    enddo

    ! ----- xn, yn "low - level" wind projections in zonal &
    ! meridional directions
    ! ----- ulow "low - level" wind magnitude - (= u)
    ! ----- bnv2 bnv2 = n ** 2
    ! ----- taub base momentum flux
    ! ----- = - (ro * u ** 3 / (n * xl) * gf (fr) for n ** 2 > 0
    ! ----- = 0. for n ** 2 < 0
    ! ----- fr froude = n * hprime / u
    ! ----- g gmax * fr ** 2 / (fr ** 2 + cg / oc)

    ! ----- initialize some arrays

    do i = 1, npt
        xn (i) = 0.0
        yn (i) = 0.0
        taub (i) = 0.0
        ulow (i) = 0.0
        dtfac (i) = 1.0
        icrilv (i) = .false. ! initialize critical level control vector


        ! ---- compute the "low level" wind magnitude (m / s)

        ulow (i) = max (sqrt (ubar (i) * ubar (i) + vbar (i) * vbar (i)), 1.0)
        uloi (i) = 1.0 / ulow (i)
    enddo

    do k = 1, kmm1
        do i = 1, npt
            j = ipt (i)
            velco (i, k) = 0.5 * ((u1 (j, k) + u1 (j, k + 1)) * ubar (i) &
                 + (v1 (j, k) + v1 (j, k + 1)) * vbar (i))
            velco (i, k) = velco (i, k) * uloi (i)
            ! if ((velco (i, k) .lt.veleps) .and. (velco (i, k) .gt.0.)) then
            ! velco (i, k) = veleps
            ! endif
        enddo
    enddo


    ! find the interface level of the projected wind where
    ! low levels & upper levels meet above pbl

    ! do i = 1, npt
    ! kint (i) = km
    ! enddo
    ! do k = 1, kmm1
    ! do i = 1, npt
    ! if (k .gt. kref (i)) then
    ! if (velco (i, k) .lt. veleps .and. kint (i) .eq. km) then
    ! kint (i) = k + 1
    ! endif
    ! endif
    ! enddo
    ! enddo
    ! warning kint = kref !!!!!!!!
    do i = 1, npt
        kint (i) = kref (i)
    enddo

    ! if (lprnt) print *, ' ubar = ', ubar, &
    ! ' vbar = ', vbar, ' ulow = ', ulow, ' veleps = ', veleps

    do i = 1, npt
        j = ipt (i)
        bnv = sqrt (bnv2bar (i))
        fr = bnv * uloi (i) * min (hprime (j), hpmax)
        fr = min (fr, frmax)
        xn (i) = ubar (i) * uloi (i)
        yn (i) = vbar (i) * uloi (i)

        ! compute the base level stress and store it in taub
        ! calculate enhancement factor, number of mountains & aspect
        ! ratio const. use simplified relationship between standard
        ! deviation & critical hgt

        ! > - calculate enhancement factor (e), number of mountans (m') and
        !! aspect ratio constant.
        !!\n as in eq. (4.9), (4.10), (4.11) in kim and arakawa (1995)
        !! \cite kim_and_arakawa_1995, we define m' and e in such a way that they
        !! depend on the geometry and location of the subgrid - scale orography
        !! through oa and the nonlinearity of flow above the orography through
        !! fr. oc, which is the orographic convexity, and statistically
        !! determine how protruded (sharp) the subgrid - scale orography is, is
        !! included in the saturation flux g' in such a way that g' is
        !! proportional to oc. the forms of e, m' and g' are:
        !!\f[
        !! e (oa, f_{r_{0}}) = (oa + 2) ^{\delta}
        !!\f]
        !!\f[
        !! \delta = c_{e}f_{r_{0}} / f_{r_{c}}
        !!\f]
        !!\f[
        !! m' (oa, clx) = c_{m}\triangle x (1 + clx) ^{oa + 1}
        !!\f]
        !!\f[
        !! g' (oc, f_{r_{0}}) = \frac{f_{r_{0}}^2}{f_{r_{0}}^2 + a^{2}}
        !!\f]
        !!\f[
        !! a^{2} = c_{g}oc^{ - 1}
        !!\f]
        !! where \f$f_{r_{c}} (= 1) \f$ is the critical froude number,
        !! \f$f_{r_{0}}\f$ is the froude number. \f$c_{e}\f$, \f$c_{m}\f$,
        !! \f$c_{g}\f$ are constants.

        ! > - calculate the reference - level drag \f$\tau_{0}\f$ (eq. (4.8) in
        !! kim and arakawa (1995) \cite kim_and_arakawa_1995) :
        !!\f[
        !! \tau_0 = e\frac{m'}{\triangle x}\frac{\rho_{0}u_0^3}{n_{0}}g'
        !!\f]
        !! where \f$e\f$, \f$m'\f$, and \f$g'\f$ are the enhancement factor,
        !! "the number of mountains", and the flux function defined above,
        !! respectively.

        efact = (oa (i) + 2.) ** (ceofrc * fr)
        efact = min (max (efact, efmin), efmax)

        coefm = (1. + clx (i)) ** (oa (i) + 1.)

        xlinv (i) = coefm * cleff (j)

        tem = fr * fr * oc (j)
        gfobnv = gmax * tem / ((tem + cg) * bnv) ! g / n0

        taub (i) = xlinv (i) * roll (i) * ulow (i) * ulow (i) &
             * ulow (i) * gfobnv * efact ! base flux tau0

        ! tem = min (hprime (i), hpmax)
        ! taub (i) = xlinv (i) * roll (i) * ulow (i) * bnv * tem * tem

        k = max (1, kref (i) - 1)
        tem = max (velco (i, k) * velco (i, k), 0.1)
        scor (i) = bnv2 (i, k) / tem ! scorer parameter below ref level
    enddo
    ! if (lprnt) print *, ' taub = ', taub

    ! ---- set up bottom values of stress

    do k = 1, kbps
        do i = 1, npt
            if (k .le. kref (i)) taup (i, k) = taub (i)
        enddo
    enddo

    ! now compute vertical structure of the stress.

    do k = kmps, kmm1 ! vertical level k loop
        kp1 = k + 1
        do i = 1, npt

            ! ----- unstable layer if ri < ric
            ! ----- unstable layer if upper air vel comp along surf vel <= 0 (crit lay)
            ! ---- at (u - c) = 0. crit layer exists and bit vector should be set (.le.)

            if (k .ge. kref (i)) then
                icrilv (i) = icrilv (i) .or. (ri_n (i, k) .lt. ric) &
                    .or. (velco (i, k) .le. 0.0)
            endif
        enddo

        ! > - compute the drag above the reference level (\f$k\geq kref\f$) :
        !! - calculate the ratio of the scorer parameter (\f$r_{scor}\f$) .
        !! \n from a series of experiments, kim and arakawa (1995)
        !! \cite kim_and_arakawa_1995 found that the magnitude of drag divergence
        !! tends to be underestimated by the revised scheme in low - level
        !! downstream regions with wave breaking. therefore, at low levels when
        !! oa > 0 (i.e., in the "downstream" region) the saturation hypothesis
        !! is replaced by the following formula based on the ratio of the
        !! the scorer parameter:
        !!\f[
        !! r_{scor} = \min \left[\frac{\tau_i}{\tau_{i + 1}}, 1\right]
        !!\f]
        do i = 1, npt
            if (k .ge. kref (i)) then
                if (.not.icrilv (i) .and. taup (i, k) .gt. 0.0) then
                    temv = 1.0 / max (velco (i, k), 0.01)
                    ! if (oa (i) .gt. 0. .and. prsi (ipt (i), kp1) .gt.rlolev) then
                    if (oa (i) .gt.0. .and. kp1 .lt. kint (i)) then
                        scork = bnv2 (i, k) * temv * temv
                        rscor = min (1.0, scork / scor (i))
                        scor (i) = scork
                    else
                        rscor = 1.
                    endif

                    ! > - the drag above the reference level is expressed as:
                    !!\f[
                    !! \tau = \frac{m'}{\triangle x}\rho nuh_d^2
                    !!\f]
                    !! where \f$h_{d}\f$ is the displacement wave amplitude. in the absence
                    !! of wave breaking, the displacement amplitude for the \f$i^{th}\f$
                    !! layer can be expressed using the drag for the layer immediately
                    !! below. thus, assuming \f$\tau_i = \tau_{i + 1}\f$, we can get:
                    !!\f[
                    !! h_{d_i}^2 = \frac{\triangle x}{m'}\frac{\tau_{i + 1}}{\rho_{i}n_{i}u_{i}}
                    !!\f]

                    brvf = sqrt (bnv2 (i, k)) ! brunt - vaisala frequency
                    ! tem1 = xlinv (i) * (ro (i, kp1) + ro (i, k)) * brvf * velco (i, k) * 0.5
                    tem1 = xlinv (i) * (ro (i, kp1) + ro (i, k)) * brvf * 0.5&
                         * max (velco (i, k), 0.01)
                    hd = sqrt (taup (i, k) / tem1)
                    fro = brvf * hd * temv

                    ! rim is the minimum - richardson number by shutts (1985)

                    ! > - the minimum richardson number (\f$ri_{m}\f$) or local
                    !! wave - modified richardson number, which determines the onset of wave
                    !! breaking, is expressed in terms of \f$r_{i}\f$ and
                    !! \f$f_{r_{d}} = nh_{d} / u\f$:
                    !!\f[
                    !! ri_{m} = \frac{ri (1 - fr_{d}) }{ (1 + \sqrt{ri}\cdot fr_{d}) ^{2}}
                    !!\f]
                    !! see eq. (4.6) in kim and arakawa (1995) \cite kim_and_arakawa_1995.

                    tem2 = sqrt (ri_n (i, k))
                    tem = 1. + tem2 * fro
                    rim = ri_n (i, k) * (1. - fro) / (tem * tem)

                    ! check stability to employ the 'saturation hypothesis'
                    ! of lindzen (1981) except at tropospheric downstream regions

                    ! > - check stability to employ the 'saturation hypothesis' of lindzen
                    !! (1981) \cite lindzen_1981 except at tropospheric downstream regions.
                    !! \n wave breaking occurs when \f$ri_{m} < ri_{c} = 0.25\f$. then
                    !! lindzen's wave saturation hypothesis resets the displacement
                    !! amplitude \f$h_{d}\f$ to that corresponding to \f$ri_{m} = 0.25\f$,
                    !! we obtain the critical \f$h_{d}\f$ (or \f$h_{c}\f$) expressed in
                    !! terms of the mean values of \f$u\f$, \f$n\f$, and \f$ri\f$ (
                    !! eq. (4.7) in kim and arakawa (1995) \cite kim_and_arakawa_1995) :
                    !!\f[
                    !! h_{c} = \frac{u}{n}\left\{2 (2 + \frac{1}{\sqrt{ri}}) ^{1 / 2} - (2 + \frac{1}{\sqrt{ri}}) \right\}
                    !!\f]
                    !! if \f$ri_{m}\leq ri_{c}\f$, obtain \f$\tau\f$ from the drag above
                    !! the reference level by using \f$h_{c}\f$ computed above; otherwise
                    !! \f$\tau\f$ is unchanged (note: scaled by the ratio of the scorer
                    !! paramter) .
                    ! ----------------------
                    if (rim .le. ric .and. (oa (i) .le. 0. .or. kp1 .ge. kint (i))) then
                        ! & if (rim .le. ric .and. (oa (i) .le. 0. .or. prsi (ipt (i), kp1) .le.rlolev)) then
                        temc = 2.0 + 1.0 / tem2
                        hd = velco (i, k) * (2. * sqrt (temc) - temc) / brvf
                        taup (i, kp1) = tem1 * hd * hd
                    else
                        taup (i, kp1) = taup (i, k) * rscor
                    endif
                    taup (i, kp1) = min (taup (i, kp1), taup (i, k))
                endif
            endif
        enddo
    enddo

    ! do i = 1, im
    ! taup (i, km + 1) = taup (i, km)
    ! enddo

    if (lcap .le. km) then
        do klcap = lcapp1, km + 1
            do i = 1, npt
                sira = prsi (ipt (i), klcap) / prsi (ipt (i), lcap)
                taup (i, klcap) = sira * taup (i, lcap)
            enddo
        enddo
    endif
    ! sjl: linear decay above p_crit, becoming constant at 1 mb
    ! angular momentum conservation is ensured, except the top leakage
    ! ----------------------- sjl mod ------------------------------
    if (p_crit > 1.e-10) then
        do i = 1, npt
            j = ipt (i)
            do k = km / 2, km + 1
                if (prsi (j, k) < p_crit) then ! scale it to zero @ top
                    taup (i, k) = taup (i, k) * (prsi (j, k) - prsi (j, km + 1)) / &
                         (p_crit - prsi (j, km + 1))
                elseif (prsi (j, k) < 1.e2) then
                    taup (i, k) = taup (i, k - 1) ! constant stress - > zero drag
                endif
            enddo
        enddo
    endif
    ! ----------------------- sjl mod ------------------------------


    ! calculate - (grav / p *) * d (tau) / d (sigma) and decel terms dtaux, dtauy

    do k = 1, km
        do i = 1, npt
            taud (i, k) = grav * (taup (i, k + 1) - taup (i, k)) / del (ipt (i), k)
        enddo
    enddo

    ! ------ limit de - acceleration (momentum deposition) at top to 1 / 2 value
    ! ------ the idea is some stuff must go out the 'top'

    if (p_crit <= 1.e-10) then
        do klcap = lcap, km
            do i = 1, npt
                taud (i, klcap) = taud (i, klcap) * factop
            enddo
        enddo
    endif

    ! ------ if the gravity wave drag would force a critical line in the
    ! ------ layers below sigma = rlolev during the next delt timestep,
    ! ------ then only apply drag until that critical line is reached.

    do k = 1, kmm1
        do i = 1, npt
            if (k .gt. kref (i) .and. prsi (ipt (i), k) .ge. rlolev) then
                if (taud (i, k) .ne.0.) then
                    tem = delt * taud (i, k)
                    dtfac (i) = min (dtfac (i), abs (velco (i, k) / tem))
                endif
            endif
        enddo
    enddo

    ! if (lprnt .and. npr .gt. 0) then
    ! print *, ' before a = ', a (npr, :)
    ! print *, ' before b = ', b (npr, :)
    ! endif

    ! > - calculate outputs: a, b, dusfc, dvsfc (see parameter description) .
    !! - below the dividing streamline height (k < idxzb), mountain
    !! blocking (\f$d_{b}\f$) is applied.
    !! - otherwise (k >= idxzb), orographic gwd (\f$\tau\f$) is applied.
    do k = 1, km
        do i = 1, npt
            j = ipt (i)
            taud (i, k) = taud (i, k) * dtfac (i)
            dtaux = taud (i, k) * xn (i)
            dtauy = taud (i, k) * yn (i)
            eng0 = 0.5 * (u1 (j, k) ** 2 + v1 (j, k) ** 2)
            ! --- lm mb (* j *) changes overwrite gwd
            if (k .lt. idxzb (i) .and. idxzb (i) .ne. 0) then
                dbim = db (i, k) / (1. + db (i, k) * delt)
                if (present (vtgwd)) vtgwd (j, k) = - dbim * v1 (j, k)
                if (present (utgwd)) utgwd (j, k) = - dbim * u1 (j, k)
                if (present (dusfc)) dusfc (j) = dusfc (j) - dbim * u1 (j, k) * del (j, k)
                if (present (dvsfc)) dvsfc (j) = dvsfc (j) - dbim * v1 (j, k) * del (j, k)
                v1 (j, k) = v1 (j, k) - dbim * v1 (j, k) * delt
                u1 (j, k) = u1 (j, k) - dbim * u1 (j, k) * delt
                !eng1 = eng0 * (1.0 - dbim * delt) * (1.0 - dbim * delt)
                eng1 = 0.5 * (u1 (j, k) ** 2 + v1 (j, k) ** 2)
            else
                if (present (vtgwd)) vtgwd (j, k) = dtauy
                if (present (utgwd)) utgwd (j, k) = dtaux
                if (present (dusfc)) dusfc (j) = dusfc (j) + dtaux * del (j, k)
                if (present (dvsfc)) dvsfc (j) = dvsfc (j) + dtauy * del (j, k)
                v1 (j, k) = v1 (j, k) + dtauy * delt
                u1 (j, k) = u1 (j, k) + dtaux * delt
                eng1 = 0.5 * (u1 (j, k) ** 2 + v1 (j, k) ** 2)
            endif
            if (present (ttgwd)) ttgwd (j, k) = (eng0 - eng1) / cp_air / delt
            t1 (j, k) = t1 (j, k) + (eng0 - eng1) / cp_air
        enddo
    enddo
    ! if (lprnt) then
    ! print *, ' in gwdps_lm.f after a = ', a (ipr, :)
    ! print *, ' in gwdps_lm.f after b = ', b (ipr, :)
    ! print *, ' db = ', db (ipr, :)
    ! endif
    tem = - 1.0 / grav
    do i = 1, npt
        j = ipt (i)
        ! tem = (- 1.e3 / grav)
        if (present (dusfc)) dusfc (j) = tem * dusfc (j)
        if (present (dvsfc)) dvsfc (j) = tem * dvsfc (j)
    enddo

    ! monitor for excessive gravity wave drag tendencies if ncnt > 0

    ! if (ncnt.gt.0) then
    ! if (lat.ge.38.and.lat.le.42) then
    !cmic$ guard 37
    ! do 92 i = 1, im
    ! if (ikount.gt.ncnt) go to 92
    ! if (i.lt.319.or.i.gt.320) go to 92
    ! do 91 k = 1, km
    ! if (abs (taud (i, k)) .gt. critac) then
    ! if (i.le.im) then
    ! ikount = ikount + 1
    ! print 123, i, lat, kdt
    ! print 124, taub (i), bnv (i), ulow (i),
    ! 1 gf (i), fr (i), roll (i), hprime (i), xn (i), yn (i)
    ! print 124, (taud (i, kk), kk = 1, km)
    ! print 124, (taup (i, kk), kk = 1, km + 1)
    ! print 124, (ri_n (i, kk), kk = 1, km)
    ! do 93 kk = 1, kmm1
    ! velko (kk) =
    ! 1 0.5 * ((u1 (i, kk) + u1 (i, kk + 1)) * ubar (i) +
    ! 2 (v1 (i, kk) + v1 (i, kk + 1)) * vbar (i)) * uloi (i)
    !93 continue
    ! print 124, (velko (kk), kk = 1, kmm1)
    ! print 124, (a (i, kk), kk = 1, km)
    ! print 124, (dtauy (i, kk), kk = 1, km)
    ! print 124, (b (i, kk), kk = 1, km)
    ! print 124, (dtaux (i, kk), kk = 1, km)
    ! go to 92
    ! endif
    ! endif
    !91 continue
    !92 continue
    !cmic$ end guard 37
    !123 format (' *** migwd print *** i = ', i3, ' lat = ', i3, ' kdt = ', i3)
    !124 format (2x, 10e13.6)
    ! endif
    ! endif

    ! print *, ' in gwdps_lm.f 18 = ', a (ipt (1), idxzb (1)), &
    ! b (ipt (1), idxzb (1)), me

end subroutine sa_gwd_oro

! =======================================================================
! Stationary convection forced gravity wave drag based on chun and
! baik (1998) \cite chun_and_baik_1998
!
! This subroutine is the parameterization of convective gravity wave
! drag based on the theory given by Chun and Baik (1998)
! \cite chun_and_baik_1998 modified for implementation into the
! GFS / CFS by Ake Johansson (Aug 2005).
!
! Parameterizing subgrid-scale convection-induced gravity wave
! momentum flux for use in large-scale models inherently requires
! some information from subgrid-scale cumulus parameterization.
! The methodology for parameterizing the zonal momentum flux induced
! by thermal forcing can be summarized as follows. From the cloud-base
! to cloud-top height, the effect of the momentum flux
! induced by subgrid-scale diabatic forcing is not considered because
! subgrid-scale cumulus convection in large-scale models is only
! activated in a conditionally unstable atmosphere. Below the cloud
! base, the momentum flux is also not considered because of the wave
! momentum cancellation. At the cloud top, the momentum flux is
! obtained by eq. (18) and (19) in Chun and Baik (1998)
! \cite chun_and_baik_1998. Above the cloud top, there are two ways to
! construct the momentum flux profile. One way is to specify a
! vertical structure of the momentum flux normalized by the cloud-top
! value, similar to what has been done for mountain drag
! parameterization. The other way is to apply the wave saturation
! hypothesis in order to find wave breaking levels in terms of the
! richardon number criterion using the nonlinearity factor of
! thermally induced waves.
!
!-----------------------------------------------------------------------
! \param[in] IM      horizontal number of used pts
! \param[in] KM      vertical layer dimension
! \param[in] U1      u component of layer wind
! \param[in] V1      v component of layer wind
! \param[in] T1      layer mean temperature (k)
! \param[in] Q1      layer mean tracer concentration
! \param[in] PRSL    mean layer pressure
! \param[in] PRSI    pressure at layer interfaces
! \param[in] DEL     mean layer delta p
! \param[in] QMAX    maximum convective heating rate (k / s) in a
!                    horizontal grid point calculated
!                    from cumulus parameterization
! \param[in] KTOP    vertical level index for cloud top
! \param[in] KBOT    vertical level index for cloud bottom
! \param[in] KCNV    (0, 1) dependent on whether convection occur or not
! \param[in] CLDF    deep convective cloud fraction at the cloud top
! \param[in] DLENGTH grid spacing in the direction of basic wind at the cloud top
! \param[out] UTGWC  zonal wind tendency
! \param[out] VTGWC  meridional wind tendency
! \param[out] TAUCTX wave stress at the cloud top projected in the east
! \param[out] TAUCTY wave stress at the cloud top projected in the north
!
!-----------------------------------------------------------------------
! Aug 2005 Ake Johansson - original code for parameterization of convectively forced
! gravity wave drag from Yonsei university, Korea
! based on the theory given by Chun and Baik (JAS, 1998)
! modified for implementation into the GFS / CFSD by
! 2013 S. Moorthi - updated and optimized code for T1534 GFS implementation
! ??? ?? 2015 J. Alpert - reducing the magnitude of tauctmax to fix blow up in L64 GFS
! S. Kar & M. Young
! Aug 15 2016 - S. Moorthi - fix for exessive dissipation which led to blow up in
! 128 level runs with NEMS / GSM
!
!-----------------------------------------------------------------------
! ARGUMENTS
!
! input variables
!
! U       : midpoint zonal wind
! V       : midpoint meridional wind
! T       : midpoint temperatures
! PMID    : midpoint pressures
! PINT    : interface pressures
! DPMID   : midpoint delta p (pi (k) - pi (k - 1))
! QMAX    : deep convective heating
! KCLDTOP : vertical level index for cloud top (mid level)
! KCLDBOT : vertical level index for cloud bottom (mid level)
! KCNV    : (0, 1) dependent on whether convection occur or not
!
! output variables
!
! UTGWC : zonal wind tendency
! VTGWC : meridional wind tendency
!
!-----------------------------------------------------------------------
! LOCAL WORKSPACE
!
! i, k     : loop index
! kk       : loop index
! cldf     : deep convective cloud fraction at the cloud top.
! ugwdc    : zonal wind after gwdc paramterization
! vgwdc    : meridional wind after gwdc parameterization
! plnmid   : log (pmid) (mid level)
! plnint   : log (pint) (interface level)
! dpint    : delta pmid (interface level)
! tauct    : wave stress at the cloud top calculated using basic - wind
!            parallel to the wind vector at the cloud top (mid level)
! tauctx   : wave stress at the cloud top projected in the east
! taucty   : wave stress at the cloud top projected in the north
! qmax     : maximum deep convective heating rate (k s - 1) in a
!            horizontal grid point calculated from cumulus para -
!            meterization. (mid level)
! wtgwc    : wind tendency in direction to the wind vector at the cloud top level
!            due to convectively generated gravity waves (mid level)
! utgwcl   : zonal wind tendency due to convectively generated
!            gravity waves (mid level)
! vtgwcl   : meridional wind tendency due to convectively generated
!            gravity waves (mid level)
! taugwci  : profile of wave stress calculated using basic - wind
!            parallel to the wind vector at the cloud top
! taugwcxi : profile of zonal component of gravity wave stress
! taugwcyi : profile of meridional component of gravity wave stress
!
! taugwci, taugwcxi, and taugwcyi are defined at the interface level
!
! bruni    : brunt - vaisala frequency (interface level)
! brunm    : brunt - vaisala frequency (mid level)
! rhoi     : air density (interface level)
! rhom     : air density (mid level)
! ti       : temperature (interface level)
! basicum  : basic - wind profile. basic - wind is parallel to the wind
!            vector at the cloud top level. (mid level)
! basicui  : basic - wind profile. basic - wind is parallel to the wind
!            vector at the cloud top level. (interface level)
! riloc    : local richardson number (interface level)
! rimin    : minimum richardson number including both the basic - state
!            and gravity wave effects (interface level)
! gwdcloc  : horizontal location where the gwdc scheme is activated.
! break    : horizontal location where wave breaking is occurred.
! critic   : horizontal location where critical level filtering is
!            occurred.
! dogwdc   : logical flag whether the gwdc parameterization is
!            calculated at a grid point or not.
!
! dogwdc is used in order to lessen cpu time for gwdc calculation.
!
!-----------------------------------------------------------------------
! Local Variables
!
! ucltop    : zonal wind at the cloud top (mid level)
! vcltop    : meridional wind at the cloud top (mid level)
! windcltop : wind speed at the cloud top (mid level)
! shear     : vertical shear of basic wind
! cosphi    : cosine of angle of wind vector at the cloud top
! sinphi    : sine of angle of wind vector at the cloud top
! c1        : tunable parameter
! c2        : tunable parameter
! dlength   : grid spacing in the direction of basic wind at the cloud top
! nonlinct  : nonlinear parameter at the cloud top
! nonlin    : nonlinear parameter above the cloud top
! nonlins   : saturation nonlinear parameter
! taus      : saturation gravity wave drag == taugwci (i, k)
! n2        : square of brunt - vaisala frequency
! dtdp      : dt / dp
! xstress   : vertically integrated zonal momentum change due to gwdc
! ystress   : vertically integrated meridional momentum change due to gwdc
! crit1     : variable 1 for checking critical level
! crit2     : variable 2 for checking critical level
! =======================================================================

subroutine sa_gwd_cnv (im, km, u1, v1, t1, q1, delt, gsize, qmax, &
        prsl, prsi, del, ktop, kbot, kcnv, &
        utgwc, vtgwc, ttgwc, tauctx, taucty)

    implicit none

    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------

    integer, intent (in) :: im, km
    integer, intent (in) :: ktop (im), kbot (im), kcnv (im)

    real, intent (in) :: delt
    real, intent (in) :: gsize (im), qmax (im)
    real, intent (in) :: prsl (im, km), prsi (im, km + 1), del (im, km)

    real, intent (inout) :: u1 (im, km), v1 (im, km), t1 (im, km), q1 (im, km)

    real, intent (out), optional :: utgwc (im, km), vtgwc (im, km), ttgwc (im, km)

    real, intent (out), optional :: tauctx (im), taucty (im)

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    ! cumchr1 (im, km)

    integer :: i, ii, k, k1, kk, kb, ilev, npt, kcb, kcldm
    integer :: ipt (im)

    real :: cldf (im), dlength (im)

    real :: tem, tem1, tem2, qtem, wtgwc, tauct, &
        windcltop, shear, nonlinct, nonlin, nonlins, &
        n2, dtdp, crit1, crit2, p1, p2, &
        ! n2, dtdp, crit1, crit2, p1, p2, &
        gsqr, onebg, eng0, eng1
        ! taus, n2, dtdp, crit1, crit2, p1, p2

    integer, allocatable :: kcldtop (:), kcldbot (:)
    logical, allocatable :: do_gwc (:)
    real, allocatable :: tauctxl (:), tauctyl (:), &
        gwdcloc (:), break (:), &
        ! critic (:), &
    ! critic (:), angle (:), &
    cosphi (:), sinphi (:), &
        xstress (:), ystress (:), &
        ucltop (:), vcltop (:), &
        wrk (:), dtfac (:), &
        dlen (:), gqmcldlen (:)
    ! real, allocatable :: plnint (:, :), dpint (:, :), &
    ! taugwci (:, :), taugwcxi (:, :), &
    ! taugwcyi (:, :), bruni (:, :), &
    ! taugwcyi (:, :), bruni (:, :),
    real, allocatable :: plnint (:, :), velco (:, :), &
        taugwci (:, :), bruni (:, :), &
        rhoi (:, :), basicui (:, :), &
        ti (:, :), riloc (:, :), &
        rimin (:, :), pint (:, :)
    ! real, allocatable :: ugwdc (:, :), vgwdc (:, :),
    real, allocatable :: &
        ! plnmid (:, :), wtgwc (:, :), &
    plnmid (:, :), taugw (:, :), &
        utgwcl (:, :), vtgwcl (:, :), &
        basicum (:, :), u (:, :), v (:, :), &
        t (:, :), spfh (:, :), &
        pmid (:, :), dpmid (:, :), &
        ! pmid (:, :), cumchr (:, :), &
    brunm (:, :), rhom (:, :)

    ! copy from gsmphys/physcons.f90
    integer, parameter :: max_lon = 5000, max_lat = 2000, min_lon = 192, min_lat = 94
    real :: dxmax, dxmin, dxinv, work1 (im), work2 (im)

    real, parameter :: &
        c1 = 1.41, c2 = - 0.38, ricrit = 0.25, &
        n2min = 1.e-32, zero = 0.0, one = 1.0, &
        taumin = 1.0e-20, tauctmax = - 20., &
        ! taumin = 1.0e-20, tauctmax = - 5., &
    qmin = 1.0e-10, shmin = 1.0e-20, &
        rimax = 1.0e+20, rimaxm = 0.99e+20, &
        rimaxp = 1.01e+20, rilarge = 0.9e+20, &
        riminx = - 1.0e+20, riminm = - 1.01e+20, &
        riminp = - 0.99e+20, rismall = - 0.9e+20


    npt = 0
    do i = 1, im
        ipt (i) = 0
        if (kcnv (i) /= 0 .and. qmax (i) > zero) then
            npt = npt + 1
            ipt (npt) = i
        endif
    enddo
    do k = 1, km
        do i = 1, im
            if (present (utgwc)) utgwc (i, k) = 0.0
            if (present (vtgwc)) vtgwc (i, k) = 0.0
            if (present (ttgwc)) ttgwc (i, k) = 0.0
            ! brunm (i, k) = 0.0
            ! rhom (i, k) = 0.0
        enddo
    enddo
    do i = 1, im
        if (present (tauctx)) tauctx (i) = 0.0
        if (present (taucty)) taucty (i) = 0.0
    enddo
    if (npt == 0) return ! no gwdc calculation done

    !-------------------------------------------------------------------
    ! calculate deep convective cloud fraction at the cloud top.
    tem = rerth * rerth * (pi + pi) * pi
    dxmax = log (tem / (max_lon * max_lat))
    dxmin = log (tem / (min_lon * min_lat))
    dxinv = 1.0 / (dxmax - dxmin)
    do i = 1, im
        work1 (i) = (2.0 * log (gsize (i)) - dxmin) * dxinv
        work1 (i) = max (0.0, min (1.0, work1 (i)))
        work2 (i) = 1.0 - work1 (i)
        cldf (i) = cgwf (1) * work1 (i) + cgwf (2) * work2 (i)
    enddo
    !-------------------------------------------------------------------
    ! calculate grid spacing in the direction of basic wind at the cloud top
    do i = 1, im
        dlength(i) = sqrt (2. * gsize (i) * gsize (i))
    enddo

    ! ***********************************************************************

    ! begin gwdc

    ! ***********************************************************************

    ! -----------------------------------------------------------------------
    ! write out incoming variables
    ! -----------------------------------------------------------------------

    ! fhourpr = zero
    ! if (lprnt) then
    ! if (fhour >= fhourpr) then
    ! print *, ' '
    ! write (*, *) 'inside gwdc raw input start print at fhour = ', &
    ! fhour
    ! write (*, *) 'im km ', im, km
    ! write (*, *) 'kbot ktop qmax dlength kcnv ',
    ! + kbot (ipr), ktop (ipr), qmax (ipr), dlength (ipr), kcnv (ipr)
    ! write (*, *) 'grav cp_air rdgas ', grav, cp_air, rdgas

    ! -------- pressure levels ----------
    ! write (*, 9100)
    ! ilev = km + 1
    ! write (*, 9110) ilev, (10. * prsi (ipr, ilev))
    ! do ilev = km, 1, - 1
    ! write (*, 9120) ilev, (10. * prsl (ipr, ilev)), &
    ! (10. * del (ipr, ilev))
    ! write (*, 9110) ilev, (10. * prsi (ipr, ilev))
    ! enddo

    ! -------- u1 v1 t1 ----------
    ! write (*, 9130)
    ! do ilev = km, 1, - 1
    ! write (*, 9140) ilev, u1 (ipr, ilev), v1 (ipr, ilev), t1 (ipr, ilev)
    ! enddo

    ! print *, ' '
    ! print *, ' inside gwdc raw input end print'
    ! endif
    ! endif

    !9100 format (//, 14x, 'pressure levels', //,
    ! + ' ilev', 6x, 'prsi', 7x, 'prsl', 6x, 'del', /)
    !9110 format (i4, 2x, f10.3)
    !9120 format (i4, 12x, 2 (2x, f10.3))
    !9130 format (//, ' ilev', 7x, 'u1', 10x, 'v1', 10x, 't1', /)
    !9140 format (i4, 3 (2x, f10.3))

    ! allocate local arrays

    allocate (kcldtop (npt), kcldbot (npt), do_gwc (npt))
    allocate (tauctxl (npt), tauctyl (npt), dtfac (npt), &
        gwdcloc (npt), break (npt), cosphi (npt), &
        ! gwdcloc (npt), break (npt), critic (npt), cosphi (npt), &
    sinphi (npt), xstress (npt), ystress (npt), wrk (npt), &
        ucltop (npt), vcltop (npt), dlen (npt), gqmcldlen (npt))

    ! allocate (plnint (npt, 2:km + 1), dpint (npt, km + 1), &
    ! taugwci (npt, km + 1), taugwcxi (npt, km + 1), &
    ! taugwcyi (npt, km + 1), bruni (npt, km + 1),
    allocate (plnint (npt, 2:km + 1), &
        taugwci (npt, km + 1), bruni (npt, km + 1), &
        rhoi (npt, km + 1), basicui (npt, km + 1), &
        ti (npt, km + 1), riloc (npt, km + 1), &
        rimin (npt, km + 1), pint (npt, km + 1))

    ! allocate (ugwdc (npt, km), vgwdc (npt, km),
    allocate &
        ! (plnmid (npt, km), wtgwc (npt, km), &
     (plnmid (npt, km), velco (npt, km), &
        utgwcl (npt, km), vtgwcl (npt, km), &
        basicum (npt, km), u (npt, km), v (npt, km), &
        t (npt, km), spfh (npt, km), pmid (npt, km), &
        dpmid (npt, km), taugw (npt, km), &
        ! dpmid (npt, km), cumchr (npt, km), &
    brunm (npt, km), rhom (npt, km))

    ! -----------------------------------------------------------------------
    ! > - # create local arrays with reversed vertical indices
    !! and initialize local variables
    ! -----------------------------------------------------------------------
    gsqr = grav * grav
    onebg = one / grav

    ! if (lprnt) then
    ! npr = 1
    ! do i = 1, npt
    ! if (ipr == ipt (i)) then
    ! npr = i
    ! exit
    ! endif
    ! enddo
    ! endif

    do k = 1, km
        k1 = km - k + 1
        do i = 1, npt
            ii = ipt (i)
            u (i, k) = u1 (ii, k1)
            v (i, k) = v1 (ii, k1)
            t (i, k) = t1 (ii, k1)
            spfh (i, k) = max (q1 (ii, k1), qmin)
            pmid (i, k) = prsl (ii, k1)
            dpmid (i, k) = del (ii, k1) * onebg
            ! cumchr (i, k) = cumchr1 (ii, k1)

            rhom (i, k) = pmid (i, k) / (rdgas * t (i, k) * (1.0 + zvir * spfh (i, k)))
            plnmid (i, k) = log (pmid (i, k))
            utgwcl (i, k) = zero
            vtgwcl (i, k) = zero
            ! ugwdc (i, k) = zero
            ! vgwdc (i, k) = zero
            brunm (i, k) = zero
            basicum (i, k) = zero
        enddo
    enddo

    do k = 1, km + 1
        k1 = km - k + 2
        do i = 1, npt
            ii = ipt (i)
            pint (i, k) = prsi (ii, k1)
            taugwci (i, k) = zero
            bruni (i, k) = zero
            rhoi (i, k) = zero
            ti (i, k) = zero
            basicui (i, k) = zero
            riloc (i, k) = zero
            rimin (i, k) = zero
        enddo
    enddo
    do k = 2, km + 1
        do i = 1, npt
            plnint (i, k) = log (pint (i, k))
        enddo
    enddo

    do i = 1, npt
        ii = ipt (i)
        kcldtop (i) = km - ktop (ii) + 1
        kcldbot (i) = km - kbot (ii) + 1
        dlen (i) = dlength (ii)
        ! (grav * qmax (ii) * cldf (ii) * dlength (ii))
        gqmcldlen (i) = grav * qmax (ii) * cldf (ii) * dlen (i)
    enddo
    ! if (lprnt) write (7000, *) ' ktop = ', ktop (ipr), ' kbot = ', kbot (ipr), &
    ! ' kcldtop = ', kcldtop (npr), ' kcldbot = ', kcldbot (npr), &
    ! ' dlength = ', dlength (ipr), ' qmax = ', qmax (ipr), ' cldf = ', cldf (ipr)

    ! if (lprnt) then
    ! if (fhour.ge.fhourpr) then
    ! write (*, 9200)
    ! do i = 1, im
    ! write (*, 9201) kcnv (i), kcldbot (i), kcldtop (i)
    ! enddo
    ! endif
    ! endif

    !9200 format (//, ' inside gwdc local variables start print', //,
    ! + 2x, 'kcnv', 2x, 'kcldbot', 2x, 'kcldtop', //)
    !9201 format (i4, 2x, i5, 4x, i5)

    ! ***********************************************************************

    ! -----------------------------------------------------------------------

    ! pressure variables

    ! interface 1 ======== pint (1) *********
    ! mid - level 1 -------- pmid (1) dpmid (1)
    ! 2 ======== pint (2) dpint (2)
    ! 2 -------- pmid (2) dpmid (2)
    ! 3 ======== pint (3) dpint (3)
    ! 3 -------- pmid (3) dpmid (3)
    ! 4 ======== pint (4) dpint (4)
    ! 4 -------- pmid (4) dpmid (4)
    ! ........
    ! 17 ======== pint (17) dpint (17)
    ! 17 -------- pmid (17) dpmid (17)
    ! 18 ======== pint (18) dpint (18)
    ! 18 -------- pmid (18) dpmid (18)
    ! 19 ======== pint (19) *********

    ! -----------------------------------------------------------------------

    do i = 1, npt
        tauctxl (i) = zero
        tauctyl (i) = zero

        ! -----------------------------------------------------------------------
        ! thermal variables

        ! interface 1 ======== ti (1) rhoi (1) bruni (1)
        ! 1 -------- t (1) rhom (1) brunm (1)
        ! 2 ======== ti (2) rhoi (2) bruni (2)
        ! 2 -------- t (2) rhom (2) brunm (2)
        ! 3 ======== ti (3) rhoi (3) bruni (3)
        ! 3 -------- t (3) rhom (3) brunm (3)
        ! 4 ======== ti (4) rhoi (4) bruni (4)
        ! 4 -------- t (4) rhom (4) brunm (4)
        ! ........
        ! 17 ========
        ! 17 -------- t (17) rhom (17) brunm (17)
        ! 18 ======== ti (18) rhoi (18) bruni (18)
        ! 18 -------- t (18) rhom (18) brunm (18)
        ! 19 ======== ti (19) rhoi (19) bruni (19)



        ! > - the top interface temperature, density, and brunt - vaisala
        !! frequencies (\f$n\f$) are calculated assuming an isothermal
        !! atmosphere above the top mid level.

        ti (i, 1) = t (i, 1)
        rhoi (i, 1) = pint (i, 1) / (rdgas * ti (i, 1))
        bruni (i, 1) = sqrt (gsqr / (cp_air * ti (i, 1)))

        ! > - the bottom interface temperature, density, and brunt - vaisala
        !! frequencies (\f$n\f$) are calculated assuming an isothermal
        !! atmosphere below the bottom mid level.

        ti (i, km + 1) = t (i, km)
        rhoi (i, km + 1) = pint (i, km + 1) / (rdgas * ti (i, km + 1) * (1.0 + zvir * spfh (i, km)))
        bruni (i, km + 1) = sqrt (gsqr / (cp_air * ti (i, km + 1)))
    enddo

    ! -----------------------------------------------------------------------

    ! > - the interface level temperature, density, and brunt - vaisala
    !! frequencies (\f$n\f$) are calculated based on linear interpolation
    !! of temperature in ln (p) .

    ! -----------------------------------------------------------------------

    do k = 2, km
        do i = 1, npt
            tem1 = (plnmid (i, k) - plnint (i, k)) / (plnmid (i, k) - plnmid (i, k - 1))
            tem2 = one - tem1
            ti (i, k) = t (i, k - 1) * tem1 + t (i, k) * tem2
            qtem = spfh (i, k - 1) * tem1 + spfh (i, k) * tem2
            rhoi (i, k) = pint (i, k) / (rdgas * ti (i, k) * (1.0 + zvir * qtem))
            dtdp = (t (i, k) - t (i, k - 1)) / (pmid (i, k) - pmid (i, k - 1))
            n2 = gsqr / ti (i, k) * (1. / cp_air - rhoi (i, k) * dtdp)
            bruni (i, k) = sqrt (max (n2min, n2))
        enddo
    enddo

    deallocate (spfh)
    ! -----------------------------------------------------------------------

    ! > - the mid - level brunt - vaisala frequencies (\f$n\f$) are calculated
    !! based on interpolated interface temperatures.
    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, npt
            dtdp = (ti (i, k + 1) - ti (i, k)) / (pint (i, k + 1) - pint (i, k))
            n2 = gsqr / t (i, k) * (1. / cp_air - rhom (i, k) * dtdp)
            brunm (i, k) = sqrt (max (n2min, n2))
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! printout
    ! -----------------------------------------------------------------------

    ! if (lprnt) then
    ! if (fhour.ge.fhourpr) then

    ! -------- pressure levels ----------
    ! write (*, 9101)
    ! do ilev = 1, km
    ! write (*, 9111) ilev, (0.01 * pint (ipr, ilev)), &
    ! (0.01 * dpint (ipr, ilev)), plnint (ipr, ilev)
    ! write (*, 9121) ilev, (0.01 * pmid (ipr, ilev)), &
    ! (0.01 * dpmid (ipr, ilev)), plnmid (ipr, ilev)
    ! enddo
    ! ilev = km + 1
    ! write (*, 9111) ilev, (0.01 * pint (ipr, ilev)), &
    ! (0.01 * dpint (ipr, ilev)), plnint (ipr, ilev)

    ! 2
    ! -------- u v t n ----------
    ! write (*, 9102)
    ! do ilev = 1, km
    ! write (*, 9112) ilev, ti (ipr, ilev), (100. * bruni (ipr, ilev))
    ! write (*, 9122) ilev, u (ipr, ilev), v (ipr, ilev),
    ! + t (ipr, ilev), (100. * brunm (ipr, ilev))
    ! enddo
    ! ilev = km + 1
    ! write (*, 9112) ilev, ti (ipr, ilev), (100. * bruni (ipr, ilev))

    ! endif
    ! endif

    !9101 format (//, 14x, 'pressure levels', //,
    ! + ' ilev', 4x, 'pint', 4x, 'pmid', 4x, 'dpint', 3x, 'dpmid', 5x, 'lnp', /)
    !9111 format (i4, 1x, f8.2, 9x, f8.2, 9x, f8.2)
    !9121 format (i4, 9x, f8.2, 9x, f8.2, 1x, f8.2)
    !9102 format (// ' ilev', 5x, 'u', 7x, 'v', 5x, 'ti', 7x, 't',
    ! + 5x, 'bruni', 3x, 'brunm', //)
    !9112 format (i4, 16x, f8.2, 8x, f8.3)
    !9122 format (i4, 2f8.2, 8x, f8.2, 8x, f8.3)


    ! ***********************************************************************

    ! big loop over grid points only done if kcnv = 1

    ! ***********************************************************************

    kcldm = 1
    do i = 1, npt
        kk = kcldtop (i)
        kb = kcldbot (i)
        kcldm = max (kcldm, kk)

        ! -----------------------------------------------------------------------

        ! > - # calculate the cloud top wind components and speed.
        !! here, ucltop, vcltop, and windcltop are wind components and
        !! wind speed at mid - level cloud top index

        ! -----------------------------------------------------------------------

        ucltop (i) = u (i, kk)
        vcltop (i) = v (i, kk)
        ! windcltop = sqrt (ucltop (i) * ucltop (i) + vcltop (i) * vcltop (i))
        windcltop = 1.0 / sqrt (ucltop (i) * ucltop (i) &
             + vcltop (i) * vcltop (i))
        cosphi (i) = ucltop (i) * windcltop
        sinphi (i) = vcltop (i) * windcltop
        ! angle (i) = acos (cosphi) * 180. / pi
    enddo

    ! -----------------------------------------------------------------------

    ! > - # calculate the basic state wind projected in the direction of the
    !! cloud top wind at mid level and interface level (u, ui), where:
    !! \n u : basic - wind speed profile. basic - wind is parallel to the wind
    !! vector at the cloud top level. (mid level)
    !! \n ui: basic - wind speed profile. basic - wind is parallel to the wind
    !! vector at the cloud top level. (interface level)
    ! input u (i, k) and v (i, k) is defined at mid level

    ! -----------------------------------------------------------------------

    do k = 1, km
        do i = 1, npt
            basicum (i, k) = u (i, k) * cosphi (i) + v (i, k) * sinphi (i)
        enddo
    enddo

    ! -----------------------------------------------------------------------

    ! basic state wind at interface level is also calculated
    ! based on linear interpolation in ln (pressure)

    ! in the top and bottom boundaries, basic - state wind at interface level
    ! is assumed to be vertically uniform.

    ! -----------------------------------------------------------------------

    do i = 1, npt
        basicui (i, 1) = basicum (i, 1)
        basicui (i, km + 1) = basicum (i, km)
    enddo
    do k = 2, km
        do i = 1, npt
            tem1 = (plnmid (i, k) - plnint (i, k)) / (plnmid (i, k) - plnmid (i, k - 1))
            tem2 = one - tem1
            basicui (i, k) = basicum (i, k) * tem2 + basicum (i, k - 1) * tem1
        enddo
    enddo

    ! -----------------------------------------------------------------------

    ! > - # calculate the local richardson number
    !! \f[
    !! ri = n^2 / \eta^2
    !! \f]
    !! where \f$\eta\f$ is the vertical shear (\f$du / dz\f$) .

    ! basicum : u at mid level
    ! basicui : ui at interface level

    ! interface 1 ======== ui (1) rhoi (1) bruni (1) riloc (1)
    ! mid - level 1 -------- u (1)
    ! 2 ======== ui (2) dpint (2) rhoi (2) bruni (2) riloc (2)
    ! 2 -------- u (2)
    ! 3 ======== ui (3) dpint (3) rhoi (3) bruni (3) riloc (3)
    ! 3 -------- u (3)
    ! 4 ======== ui (4) dpint (4) rhoi (4) bruni (4) riloc (4)
    ! 4 -------- u (4)
    ! ........
    ! 17 ======== ui (17) dpint (17) rhoi (17) bruni (17) riloc (17)
    ! 17 -------- u (17)
    ! 18 ======== ui (18) dpint (18) rhoi (18) bruni (18) riloc (18)
    ! 18 -------- u (18)
    ! 19 ======== ui (19) rhoi (19) bruni (19) riloc (19)

    ! -----------------------------------------------------------------------

    do k = 2, km
        do i = 1, npt
            shear = grav * rhoi (i, k) * (basicum (i, k) - basicum (i, k - 1)) &
                 / (pmid (i, k) - pmid (i, k - 1))
            if (abs (shear) < shmin) then
                riloc (i, k) = rimax
            else
                tem = bruni (i, k) / shear
                riloc (i, k) = tem * tem
                if (riloc (i, k) >= rimax) riloc (i, k) = rilarge
            endif
        enddo
    enddo

    do i = 1, npt
        riloc (i, 1) = riloc (i, 2)
        riloc (i, km + 1) = riloc (i, km)
    enddo

    ! if (lprnt.and. (i.eq.ipr)) then
    ! if (fhour.ge.fhourpr) then
    ! write (*, 9104) ucltop, vcltop, windcltop, angle, kk
    ! do ilev = 1, km
    ! write (*, 9114) ilev, basicui (ipr, ilev), dpint (ipr, ilev),
    ! + rhoi (ipr, ilev), (100. * bruni (ipr, ilev)), riloc (ilev)
    ! write (*, 9124) ilev, (basicum (ipr, ilev))
    ! enddo
    ! ilev = km + 1
    ! write (*, 9114) ilev, basicui (ipr, ilev), dpint (ipr, ilev),
    ! + rhoi (ipr, ilev), (100. * bruni (ipr, ilev)), riloc (ilev)
    ! endif
    ! endif

    !9104 format (//, 'wind vector at cloudtop = (', f6.2, ', ', f6.2, ') = ',
    ! + f6.2, ' in direction ', f6.2, 4x, 'kk = ', i2, //,
    ! + ' ilev', 2x, 'basicum', 2x, 'basicui', 4x, 'dpint', 6x, 'rhoi', 5x,
    ! + 'bruni', 6x, 'ri', /)
    !9114 format (i4, 10x, f8.2, 4 (2x, f8.2))
    !9124 format (i4, 1x, f8.2)

    ! -----------------------------------------------------------------------

    ! > - # calculate the gravity wave stress at the interface level cloud top.

    ! kcldtopi : the interface level cloud top index
    ! kcldtop : the midlevel cloud top index
    ! kcldbot : the midlevel cloud bottom index

    ! a : find deep convective heating rate maximum

    ! if kcldtop (i) is less than kcldbot (i) in a horizontal grid point,
    ! it can be thought that there is deep convective cloud. however,
    ! deep convective heating between kcldbot and kcldtop is sometimes
    ! zero in spite of kcldtop less than kcldbot. in this case,
    ! maximum deep convective heating is assumed to be 1.e-30.

    ! b : kk is the vertical index for interface level cloud top

    ! c : total convective fractional cover (cldf) is used as the
    ! convective cloud cover for gwdc calculation instead of
    ! convective cloud cover in each layer (concld) .
    ! a1 = cldf * dlength
    ! you can see the difference between cldf (i) and concld (i)
    ! in (4.a.2) in description of the ncar community climate
    ! model (ccm3) .
    ! in ncar ccm3, cloud fractional cover in each layer in a deep
    ! cumulus convection is determined assuming total convective
    ! cloud cover is randomly overlapped in each layer in the
    ! cumulus convection.

    ! d : wave stress at cloud top is calculated when the atmosphere
    ! is dynamically stable at the cloud top

    ! e : cloud top wave stress and nonlinear parameter are calculated
    ! using density, temperature, and wind that are defined at mid
    ! level just below the interface level in which cloud top wave
    ! stress is defined.
    ! nonlinct is defined at the interface level.

    ! f : if the atmosphere is dynamically unstable at the cloud top,
    ! gwdc calculation in current horizontal grid is skipped.

    ! g : if mean wind at the cloud top is less than zero, gwdc

    ! > - wave stress at cloud top is calculated when the atmosphere
    !! is dynamically stable at the cloud top
    !
    ! > - the cloud top wave stress and nonlinear parameter are calculated
    !! using density, temperature, and wind that are defined at mid
    !! level just below the interface level in which cloud top wave
    !! stress is defined.
    !! the parameter \f$\mu\f$ is the nonlinearity factor of thermally
    !! induced internal gravity waves defined by eq. (17) in chun and
    !! baik, 1998 \cite chun_and_baik_1998
    !! \f[
    !! \mu = \frac{gq_{0}a_{1}}{c_{p}t_{0}nu^{2}}
    !! \f]
    !! where \f$q_{0}\f$ is the maximum deep convective heating rate in a
    !! horizontal grid point calculated from cumulus parameterization.
    !! \f$a_{1}\f$ is the half - width of
    !! the forcing function.\f$g\f$ is gravity. \f$c_{p}\f$ is specific
    !! heat at constant pressure. \f$t_{0}\f$ is the layer mean
    !! temperature (t1) . as eqs. (18) and (19) \cite chun_and_baik_1998,
    !! the zonal momentum flux is given by
    !! \f[
    !! \tau_{x} = - [\rho u^{3} / (n\triangle x) ]g (\mu)
    !! \f]
    !! where
    !! \f[
    !! g (\mu) = c_{1}c_2^2 \mu^{2}
    !! \f]
    !! wher \f$\rho\f$ is the local density.
    !! the tunable parameter \f$c_1\f$ is related to the horizontal
    !! structure of thermal forcing. the tunable parameter \f$c_2\f$ is
    !! related to the basic - state wind and stability and the bottom and
    !! top heights of thermal forcing. if the atmosphere is dynamically
    !! unstable at the cloud top, the convective gwd calculation is
    !! skipped at that grid point.
    !
    ! - if mean wind at the cloud top is less than zero, gwdc
    ! calculation in current horizontal grid is skipped.


    ! > - the stress is capped at tauctmax = - 5\f$n / m^2\f$
    !! in order to prevent numerical instability.

    ! -----------------------------------------------------------------------
    !d
    do i = 1, npt
        kk = kcldtop (i)
        if (abs (basicui (i, kk)) > zero .and. riloc (i, kk) > ricrit) then
            !e
            tem = basicum (i, kk)
            tem1 = tem * tem
            nonlinct = gqmcldlen (i) / (bruni (i, kk) * t (i, kk) * tem1) ! mu
            tem2 = c2 * nonlinct
            ! rhou^3c1 (c2mu) ^2 / ndx
            tauct = - rhom (i, kk) * tem * tem1 * c1 * tem2 * tem2&
                 / (bruni (i, kk) * dlen (i))

            tauct = max (tauctmax, tauct)
            tauctxl (i) = tauct * cosphi (i) ! x stress at cloud top
            tauctyl (i) = tauct * sinphi (i) ! y stress at cloud top
            taugwci (i, kk) = tauct ! * 1
            do_gwc (i) = .true.
        else
            !f
            tauctxl (i) = zero
            tauctyl (i) = zero
            do_gwc (i) = .false.
        endif
        !h
    enddo

    ! if (lprnt.and. (i.eq.ipr)) then
    ! if (fhour.ge.fhourpr) then
    ! write (*, 9210) tauctx (ipr), taucty (ipr), tauct (ipr), angle, kk
    ! endif
    ! endif

    !9210 format (/, 5x, 'stress vector = (', f8.3, ', ', f8.3, ') = ', f8.3,
    ! + ' in direction ', f6.2, 4x, 'kk = ', i2, /)

    ! -----------------------------------------------------------------------

    ! at this point, mean wind at the cloud top is larger than zero and
    ! local ri at the cloud top is larger than ricrit (= 0.25)

    ! calculate minimum of richardson number including both basic - state
    ! condition and wave effects.

    ! g * q_0 * alpha * dx ri_loc * (1 - mu * |c2|)
    ! mu = ---------------- ri_min = -----------------------------
    ! c_p * n * t * u^2 (1 + mu * ri_loc^ (0.5) * |c2|) ^2

    ! minimum ri is calculated for the following two cases

    ! (1) riloc < 1.e+20
    ! (2) riloc = 1.e+20 ---- > vertically uniform basic - state wind

    ! riloc cannot be smaller than zero because n^2 becomes 1.e-32 in the
    ! case of n^2 < 0.. thus the sign of rinum is determined by
    ! 1 - nonlin * |c2|.

    ! -----------------------------------------------------------------------
    ! > - # calculate the minimum richardson number including both the
    !! basic - state condition and wave effects.
    !!\f[
    !! ri_{min}\approx\frac{ri (1 - \mu|c_{2}|) }{ (1 + \mu ri^{1 / 2}|c_{2}|) ^{2}}
    !!\f]

    do k = kcldm, 1, - 1

        do i = 1, npt
            if (do_gwc (i)) then
                kk = kcldtop (i)
                if (k > kk) cycle
                if (k /= 1) then
                    tem1 = (u (i, k) + u (i, k - 1)) * 0.5
                    tem2 = (v (i, k) + v (i, k - 1)) * 0.5
                    crit1 = ucltop (i) * tem1
                    crit2 = vcltop (i) * tem2
                    velco (i, k) = tem1 * cosphi (i) + tem2 * sinphi (i)
                else
                    crit1 = ucltop (i) * u (i, 1)
                    crit2 = vcltop (i) * v (i, 1)
                    velco (i, 1) = u (i, 1) * cosphi (i) + v (i, 1) * sinphi (i)
                endif
                ! if (lprnt .and. i == npr) write (7000, *) ' k = ', k, ' crit1 = ', &
                ! crit1, ' crit2 = ', crit2, ' basicui = ', basicui (i, k)

                if (abs (basicui (i, k)) > zero .and. crit1 > zero .and. crit2 > zero) then
                    tem = basicui (i, k) * basicui (i, k)
                    nonlin = gqmcldlen (i) / (bruni (i, k) * ti (i, k) * tem)
                    tem = nonlin * abs (c2)
                    if (riloc (i, k) < rimaxm) then
                        tem1 = 1 + tem * sqrt (riloc (i, k))
                        rimin (i, k) = riloc (i, k) * (1 - tem) / (tem1 * tem1)
                    else if ((riloc (i, k) > rimaxm) .and. (riloc (i, k) < rimaxp)) then
                        rimin (i, k) = (1 - tem) / (tem * tem)
                    endif
                    if (rimin (i, k) <= riminx) then
                        rimin (i, k) = rismall
                    endif
                else
                    rimin (i, k) = riminx
                endif
                ! if (lprnt .and. i == npr) write (7000, *) ' rimin = ', rimin (i, k)

                ! -----------------------------------------------------------------------

                ! if the minimum \f$r_{i}\f$ at interface cloud top is less than or equal to 1 / 4,
                ! the convective gwd calculation is skipped at that grid point.

                ! -----------------------------------------------------------------------

                ! -----------------------------------------------------------------------

                ! > - # calculate the gravity wave stress profile using the wave
                !! saturation hypothesis of lindzen (1981) \cite lindzen_1981.

                ! assuming kcldtop (i) = 10 and kcldbot = 16,

                ! taugwci riloc rimin utgwc

                ! interface 1 ======== - 0.001 - 1.e20
                ! 1 -------- 0.000
                ! 2 ======== - 0.001 - 1.e20
                ! 2 -------- 0.000
                ! 3 ======== - 0.001 - 1.e20
                ! 3 --------- .xxx
                ! 4 ======== - 0.001 2.600 2.000
                ! 4 -------- 0.000
                ! 5 ======== - 0.001 2.500 2.000
                ! 5 -------- 0.000
                ! 6 ======== - 0.001 1.500 0.110
                ! 6 -------- + .xxx
                ! 7 ======== - 0.005 2.000 3.000
                ! 7 -------- 0.000
                ! 8 ======== - 0.005 1.000 0.222
                ! 8 -------- + .xxx
                ! 9 ======== - 0.010 1.000 2.000
                ! 9 -------- 0.000
                ! kcldtopi 10 ======== $$$ - 0.010
                ! kcldtop 10 -------- $$$ yyyyy
                ! 11 ======== $$$ 0
                ! 11 -------- $$$
                ! 12 ======== $$$ 0
                ! 12 -------- $$$
                ! 13 ======== $$$ 0
                ! 13 -------- $$$
                ! 14 ======== $$$ 0
                ! 14 -------- $$$
                ! 15 ======== $$$ 0
                ! 15 -------- $$$
                ! 16 ======== $$$ 0
                ! kcldbot 16 -------- $$$
                ! 17 ======== 0
                ! 17 --------
                ! 18 ======== 0
                ! 18 --------
                ! 19 ======== 0

                ! -----------------------------------------------------------------------

                ! even though the cloud top level obtained in deep convective para -
                ! meterization is defined in mid - level, the cloud top level for
                ! the gwdc calculation is assumed to be the interface level just
                ! above the mid - level cloud top vertical level index.

                ! -----------------------------------------------------------------------

                ! > - when \f$ri_{min}\f$ is set to 1 / 4 based on lindzen's (1981)
                !! \cite lindzen_1981 saturation hypothesis, the nonlinearity factor
                !! for wave saturation can be derived by
                !! \f[
                !! \mu_{s} = \frac{1}{|c_{2}|}[2\sqrt{2 + \frac{1}{\sqrt{ri}}} - (2 + \frac{1}{\sqrt{ri}}) ]
                !! \f]
                !! then the saturation zonal momentum flux is given by
                !! \f[
                !! \tau_{s} = - [\rho u^{3} / (n\triangle x) ]c_{1}c_2^2\mu_s^2
                !! \f]

                if (k < kk .and. k > 1) then
                    if (abs (taugwci (i, k + 1)) > taumin) then ! taugwci
                        if (riloc (i, k) > ricrit) then ! riloc
                            if (rimin (i, k) > ricrit) then ! rimin
                                taugwci (i, k) = taugwci (i, k + 1)
                            elseif (rimin (i, k) > riminp) then
                                tem = 2.0 + 1.0 / sqrt (riloc (i, k))
                                nonlins = (1.0 / abs (c2)) * (2. * sqrt (tem) - tem)
                                tem1 = basicui (i, k)
                                tem2 = c2 * nonlins * tem1
                                taugwci (i, k) = - rhoi (i, k) * c1 * tem1 * tem2 * tem2&
                                     / (bruni (i, k) * dlen (i))
                            elseif (rimin (i, k) > riminm) then
                                taugwci (i, k) = zero
                                ! taugwci (i, k) = taugwci (i, k + 1)
                            endif ! rimin
                        else

                            ! > - if the minimum \f$r_{i}\f$ at interface cloud top is less than
                            !! or equal to 1 / 4, the convective gwd calculation is skipped at that
                            !! grid point.

                            taugwci (i, k) = zero
                        endif ! riloc
                    else
                        taugwci (i, k) = zero
                    endif ! taugwci

                    if ((basicum (i, k + 1) * basicum (i, k)) < 0.) then
                        taugwci (i, k + 1) = zero
                        taugwci (i, k) = zero
                    endif

                    if (abs (taugwci (i, k)) > abs (taugwci (i, k + 1))) then
                        taugwci (i, k) = taugwci (i, k + 1)
                    endif

                elseif (k == 1) then

                    ! > - as an upper boundary condition, upward propagation of gravity
                    !! wave energy is permitted.

                    taugwci (i, 1) = taugwci (i, 2)
                endif

                ! if (lprnt .and. i == npr) then
                ! write (7000, *) 'k = ', k, ' taugwci = ', taugwci (i, k), &
                ! 'riloc', riloc (i, k), 'riminp = ', riminp, ' ricrit = ', ricrit, &
                ! 'bruni (i, k) = ', bruni (i, k), ' deln = ', bruni (i, k), &
                ! 'basicui (i, k) = ', basicui (i, k), ' rimin = ', rimin (i, k), &
                ! ' dlen = ', dlen (i), ' rhoi = ', rhoi (i, k)
                ! endif

            endif
        enddo ! end of i = 1, npt loop
    enddo ! end of k = kcldm, 1, - 1 loop

    do i = 1, npt
        dtfac (i) = 1.0
    enddo
    do k = 1, km
        do i = 1, npt
            if (do_gwc (i)) then
                kk = kcldtop (i)
                if (k < kk) then
                    taugw (i, k) = (taugwci (i, k + 1) - taugwci (i, k)) / dpmid (i, k)
                    if (taugw (i, k) /= 0.0) then
                        tem = delt * taugw (i, k)
                        dtfac (i) = min (dtfac (i), abs (velco (i, k) / tem))
                    endif
                else
                    taugw (i, k) = 0.0
                endif
            else
                taugw (i, k) = 0.0
            endif
        enddo
    enddo

    !!!!!! vertical differentiation
    !!!!!
    ! > - # calculate wind tendency in direction to the wind vector, zonal
    !! wind tendency and meridional wind tendency above the cloud top
    !! level due to convectively generated gravity waves.

    do k = 1, km
        do i = 1, npt
            if (do_gwc (i)) then
                kk = kcldtop (i)
                if (k < kk) then
                    ! wtgwc = (taugwci (i, k + 1) - taugwci (i, k)) / dpmid (i, k)
                    wtgwc = taugw (i, k) * dtfac (i)
                    utgwcl (i, k) = wtgwc * cosphi (i)
                    vtgwcl (i, k) = wtgwc * sinphi (i)
                else
                    utgwcl (i, k) = zero
                    vtgwcl (i, k) = zero
                endif
                ! if (lprnt .and. i == npr) then
                ! write (7000, *) 'k = ', k, ' wtgwc = ', wtgwc, ' taugwci = ', taugwci (i, k), &
                ! taugwci (i, k + 1), ' dpmid = ', dpmid (i, k), ' cosphi = ', cosphi (i), &
                ! ' sinphi = ', sinphi (i), ' utgwcl = ', utgwcl (i, k), &
                ! 'vtgwcl = ', vtgwcl (i, k), ' dtfac = ', dtfac (i)
                ! endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------

    ! calculate momentum flux = stress deposited above cloup top
    ! apply equal amount with opposite sign within cloud

    ! -----------------------------------------------------------------------

    do i = 1, npt
        xstress (i) = zero
        ystress (i) = zero
    enddo
    do k = 1, kcldm
        do i = 1, npt
            if (do_gwc (i)) then
                xstress (i) = xstress (i) + utgwcl (i, k) * dpmid (i, k)
                ystress (i) = ystress (i) + vtgwcl (i, k) * dpmid (i, k)
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! alt 1 only uppermost layer
    ! -----------------------------------------------------------------------

    ! kk = kcldtop (i)
    ! tem1 = grav / dpmid (i, kk)
    ! utgwc (i, kk) = - tem1 * xstress
    ! vtgwc (i, kk) = - tem1 * ystress

    ! -----------------------------------------------------------------------
    ! alt 2 sin (kt - kb)
    ! -----------------------------------------------------------------------

    do i = 1, npt
        if (do_gwc (i)) then
            wrk (i) = 0.5 * pi / (pint (i, kcldbot (i) + 1) - pint (i, kcldtop (i)))
        endif
    enddo
    do k = 1, km
        do i = 1, npt
            if (do_gwc (i)) then
                kk = kcldtop (i)
                if (k >= kk .and. k <= kcldbot (i)) then
                    p1 = sin (wrk (i) * (pint (i, k) - pint (i, kk)))
                    p2 = sin (wrk (i) * (pint (i, k + 1) - pint (i, kk)))
                    tem = - (p2 - p1) / dpmid (i, k)
                    utgwcl (i, k) = tem * xstress (i)
                    vtgwcl (i, k) = tem * ystress (i)
                endif
            endif
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! alt 3 from kt to kb proportional to conv heating
    ! -----------------------------------------------------------------------

    ! do k = kcldtop (i), kcldbot (i)
    ! p1 = cumchr (i, k)
    ! p2 = cumchr (i, k + 1)
    ! utgwcl (i, k) = - grav * xstress * (p1 - p2) / dpmid (i, k)
    ! enddo

    ! -----------------------------------------------------------------------

    ! the gwdc should accelerate the zonal and meridional wind in the
    ! opposite direction of the previous zonal and meridional wind,
    ! respectively

    ! -----------------------------------------------------------------------

    ! do k = 1, kcldtop (i) - 1

    ! if (utgwcl (i, k) * u (i, k) .gt. 0.0) then

    ! -------------------- x - component -------------------

    ! write (6, ' (a) ')
    ! + ' (gwdc) warning: the gwdc should accelerate the zonal wind '
    ! write (6, ' (a, a, i3, a, i3) ')
    ! + 'in the opposite direction of the previous zonal wind',
    ! + ' at i = ', i, ' and j = ', lat
    ! write (6, ' (4 (1x, e17.10)) ') u (i, kk), v (i, kk), u (i, k), v (i, k)
    ! write (6, ' (a, 1x, e17.10)) ') 'vcld . v = ',
    ! + u (i, kk) * u (i, k) + v (i, kk) * v (i, k)

    ! if (u (i, kcldtop (i)) * u (i, k) + v (i, kcldtop (i)) * v (i, k) .gt.0.0) then
    ! do k1 = 1, km
    ! write (6, ' (i2, 36x, 2 (1x, e17.10)) ')
    ! + k1, taugwcxi (i, k1), taugwci (i, k1)
    ! write (6, ' (i2, 2 (1x, e17.10)) ') k1, utgwcl (i, k1), u (i, k1)
    ! enddo
    ! write (6, ' (i2, 36x, 1x, e17.10) ') (km + 1), taugwcxi (i, km + 1)
    ! endif

    ! -------------------- along wind at cloud top -----

    ! do k1 = 1, km
    ! write (6, ' (i2, 36x, 2 (1x, e17.10)) ')
    ! + k1, taugwci (i, k1)
    ! write (6, ' (i2, 2 (1x, e17.10)) ') k1, wtgwc (i, k1), basicum (i, k1)
    ! enddo
    ! write (6, ' (i2, 36x, 1x, e17.10) ') (km + 1), taugwci (i, km + 1)

    ! endif

    ! if (vtgwc (i, k) * v (i, k) .gt. 0.0) then
    ! write (6, ' (a) ')
    ! + ' (gwdc) warning: the gwdc should accelerate the meridional wind'
    ! write (6, ' (a, a, i3, a, i3) ')
    ! + 'in the opposite direction of the previous meridional wind',
    ! + ' at i = ', i, ' and j = ', lat
    ! write (6, ' (4 (1x, e17.10)) ') u (i, kcldtop (i)), v (i, kcldtop (i)),
    ! + u (i, k), v (i, k)
    ! write (6, ' (a, 1x, e17.10)) ') 'vcld . v = ',
    ! + u (i, kcldtop (i)) * u (i, k) + v (i, kcldtop (i)) * v (i, k)
    ! if (u (i, kcldtop (i)) * u (i, k) + v (i, kcldtop (i)) * v (i, k) .gt.0.0) then
    ! do k1 = 1, km
    ! write (6, ' (i2, 36x, 2 (1x, e17.10)) ')
    ! + k1, taugwcyi (i, k1), taugwci (i, k1)
    ! write (6, ' (i2, 2 (1x, e17.10)) ') k1, vtgwc (i, k1), v (i, k1)
    ! enddo
    ! write (6, ' (i2, 36x, 1x, e17.10) ') (km + 1), taugwcyi (i, km + 1)
    ! endif
    ! endif

    ! enddo

    !1000 continue


    ! ***********************************************************************

    ! if (lprnt) then
    ! if (fhour.ge.fhourpr) then
    ! -------- utgwc vtgwc ----------
    ! write (*, 9220)
    ! do ilev = 1, km
    ! write (*, 9221) ilev, (86400. * utgwcl (ipr, ilev)),
    ! + (86400. * vtgwcl (ipr, ilev))
    ! enddo
    ! endif
    ! endif

    !9220 format (//, 14x, 'tendency due to gwdc', //,
    ! + ' ilev', 6x, 'utgwc', 7x, 'vtgwc', /)
    !9221 format (i4, 2 (2x, f10.3))

    ! -----------------------------------------------------------------------

    ! for gwdc performance analysis

    ! -----------------------------------------------------------------------

    ! do k = 1, kk - 1
    ! do i = 1, nct

    ! kk = kcldtop (i)

    ! if ((abs (taugwci (i, kk)) > taumin)) then

    ! gwdcloc (i) = one

    ! if (abs (taugwci (i, k) - taugwci (i, kk)) > taumin) then
    ! break (i) = 1.0
    ! go to 2000
    ! endif
    ! enddo
    !2000 continue

    ! do k = 1, kk - 1

    ! if ((abs (taugwci (i, k)) .lt.taumin) .and. (abs (taugwci (i, k + 1)) .gt.taumin) .and. (basicum (i, k + 1) * basicum (i, k) .lt. 0.)) then
    ! critic (i) = 1.0
    ! print *, i, k, ' inside gwdc taugwci (k) = ', taugwci (i, k)
    ! print *, i, k + 1, ' inside gwdc taugwci (k + 1) = ', taugwci (i, k + 1)
    ! print *, i, k, ' inside gwdc basicum (k) = ', basicum (i, k)
    ! print *, i, k + 1, ' inside gwdc basicum (k + 1) = ', basicum (i, k + 1)
    ! print *, i, ' inside gwdc critic = ', critic (i)
    ! goto 2010
    ! endif
    ! enddo
    !2010 continue

    ! endif

    ! enddo

    ! -----------------------------------------------------------------------
    ! > - # convert back local convective gwd tendency arrays to gfs model
    !! vertical indices.
    ! outgoing (fu1, fv1) = (utgwc, vtgwc)
    ! -----------------------------------------------------------------------

    do k = 1, km
        k1 = km - k + 1
        do i = 1, npt
            ii = ipt (i)
            eng0 = 0.5 * (u1 (ii, k1) ** 2 + v1 (ii, k1) ** 2)
            if (present (utgwc)) utgwc (ii, k1) = utgwcl (i, k)
            if (present (vtgwc)) vtgwc (ii, k1) = vtgwcl (i, k)

            ! brunm (ii, kk) = brunm (i, k)
            ! brunm (i, k) = tem

            ! rhom (ii, kk) = rhom (i, k)
            u1 (ii, k1) = u1 (ii, k1) + utgwcl (i, k) * delt
            v1 (ii, k1) = v1 (ii, k1) + vtgwcl (i, k) * delt
            eng1 = 0.5 * (u1 (ii, k1) ** 2 + v1 (ii, k1) ** 2)
            if (present (ttgwc)) ttgwc (ii, k1) = (eng0 - eng1) / cp_air / delt
            t1 (ii, k1) = t1 (ii, k1) + (eng0 - eng1) / cp_air
        enddo
        ! if (lprnt) write (7000, *) ' k = ', k, ' k1 = ', k1, ' utgwc = ', &
        ! utgwc (ipr, k1), ' vtgwc = ', vtgwc (ipr, k1)
    enddo
    do i = 1, npt
        ii = ipt (i)
        if (present (tauctx)) tauctx (ii) = tauctxl (i)
        if (present (taucty)) taucty (ii) = tauctyl (i)
    enddo

    ! if (lprnt) then
    ! if (fhour.ge.fhourpr) then
    ! -------- utgwc vtgwc ----------
    ! write (*, 9225)
    ! do ilev = km, 1, - 1
    ! write (*, 9226) ilev, (86400. * fu1 (ipr, ilev)),
    ! + (86400. * fv1 (ipr, ilev))
    ! enddo
    ! endif
    ! endif

    !9225 format (//, 14x, 'tendency due to gwdc - to gbphys', //,
    ! + ' ilev', 6x, 'utgwc', 7x, 'vtgwc', /)
    !9226 format (i4, 2 (2x, f10.3))

    deallocate (kcldtop, kcldbot, do_gwc)
    deallocate (tauctxl, tauctyl, dtfac, &
        ! gwdcloc, break, critic, cosphi, &
        gwdcloc, break, cosphi, &
        sinphi, xstress, ystress, &
        dlen, ucltop, vcltop, gqmcldlen, wrk)

    deallocate (plnint, taugwci, velco, &
        bruni, rhoi, basicui, &
        ti, riloc, rimin, pint)

    deallocate (plnmid, utgwcl, vtgwcl, basicum, u, v, t, &
        pmid, dpmid, brunm, rhom, taugw)

end subroutine sa_gwd_cnv

end module sa_gwd_mod
