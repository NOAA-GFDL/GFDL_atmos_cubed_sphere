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
! cloud radii diagnosis built for gfdl cloud microphysics
! authors: linjiong zhou and shian - jiann lin
! =======================================================================
module cloud_diagnosis_mod

    implicit none

    private

    public cloud_diagnosis, cloud_diagnosis_init

    real, parameter :: grav = 9.80665 ! gfs: acceleration due to gravity
    real, parameter :: rdgas = 287.05 ! gfs: gas constant for dry air
    real, parameter :: rvgas = 461.50 ! gfs: gas constant for water vapor
    real, parameter :: pi = 3.1415926535897931 ! gfs: ratio of circle circumference to diameter

    real, parameter :: zvir = rvgas / rdgas - 1. ! 0.6077338443

    real :: tice = 273.16 ! set tice = 165. to trun off ice - phase phys (kessler emulator)

    real :: ql0_max = 2.0e-3 ! max cloud water value (auto converted to rain)
    real :: qi0_max = 2.0e-4 ! max cloud ice value (by other sources)
    real :: qi0_rei = 0.8e-4 ! max cloud ice value (by other sources)

    real :: ccn_o = 100. ! ccn over ocean (cm^ - 3)
    real :: ccn_l = 300. ! ccn over land (cm^ - 3)

    ! cloud diagnosis

    real :: qmin = 1.0e-12 ! minimum mass mixing ratio (kg / kg)
    ! real :: beta = 1.22 ! defined in heymsfield and mcfarquhar, 1996
    real :: beta = 1.
    ! real :: beta = 0.5 ! testing

    ! real :: rewmin = 1.0, rewmax = 25.0
    ! real :: reimin = 10.0, reimax = 300.0
    ! real :: rermin = 25.0, rermax = 225.0
    ! real :: resmin = 300, resmax = 1000.0
    ! real :: regmin = 1000.0, regmax = 1.0e5
    ! lz
    ! real :: rewmin = 5.0, rewmax = 10.0
    ! real :: reimin = 10.0, reimax = 150.0
    ! real :: rermin = 0.0, rermax = 10000.0
    ! real :: resmin = 0.0, resmax = 10000.0
    ! real :: regmin = 0.0, regmax = 10000.0
    ! sjl
    !!! real :: reimin = 10.0, reimax = 150.0
    real :: rewmin = 4.0, rewmax = 10.0
    real :: reimin = 4.0, reimax = 250.0
    real :: rermin = 5.0, rermax = 2000.0
    real :: resmin = 5.0, resmax = 2000.0
    real :: regmin = 5.0, regmax = 2000.0

    real :: betaw = 1.0
    real :: betai = 1.0
    real :: betar = 1.0
    real :: betas = 1.0
    real :: betag = 1.0

    logical :: liq_ice_combine = .true.

    integer :: rewflag = 1
    ! 1: martin et al., 1994
    ! 2: martin et al., 1994, gfdl revision
    ! 3: kiehl et al., 1994
    integer :: reiflag = 1
    ! 1: heymsfield and mcfarquhar, 1996
    ! 2: donner et al., 1997
    ! 3: fu, 2007
    ! 4: kristjansson et al., 2000
    ! 5: wyser, 1998

    namelist / cloud_diagnosis_nml / &
        ql0_max, qi0_max, qi0_rei, ccn_o, ccn_l, qmin, beta, liq_ice_combine, rewflag, reiflag, &
        rewmin, rewmax, reimin, reimax, rermin, rermax, resmin, resmax, regmin, regmax, &
        betaw, betai, betar, betas, betag

contains

! =======================================================================
! radius of cloud species diagnosis
! =======================================================================

subroutine cloud_diagnosis (is, ie, ks, ke, lsm, p, delp, t, qw, qi, qr, qs, qg, &
        qcw, qci, qcr, qcs, qcg, rew, rei, rer, res, reg, &
        cld, cloud, snowd, cnvw, cnvi, cnvc)

    implicit none

    integer, intent (in) :: is, ie
    integer, intent (in) :: ks, ke

    real, intent (in), dimension (is:ie) :: lsm ! land sea mask, 0: ocean, 1: land, 2: sea ice
    real, intent (in), dimension (is:ie) :: snowd ! snow depth (mm)

    real, intent (in), dimension (is:ie, ks:ke) :: delp, t, p
    real, intent (in), dimension (is:ie, ks:ke) :: cloud ! cloud fraction
    real, intent (in), dimension (is:ie, ks:ke) :: qw, qi, qr, qs, qg ! mass mixing ratio (kg / kg)

    real, intent (in), dimension (is:ie, ks:ke), optional :: cnvw, cnvi ! convective cloud water, cloud ice mass mixing ratio (kg / kg)
    real, intent (in), dimension (is:ie, ks:ke), optional :: cnvc ! convective cloud fraction

    real, intent (out), dimension (is:ie, ks:ke) :: qcw, qci, qcr, qcs, qcg ! units: g / m^2
    real, intent (out), dimension (is:ie, ks:ke) :: rew, rei, rer, res, reg ! radii (micron)
    real, intent (out), dimension (is:ie, ks:ke) :: cld ! total cloud fraction

    ! local variables

    integer :: i, k, ind

    real, dimension (is:ie, ks:ke) :: qmw, qmr, qmi, qms, qmg ! mass mixing ratio (kg / kg)

    real :: dpg ! dp / g
    real :: rho ! density (kg / m^3)
    real :: ccnw ! cloud condensate nuclei for cloud water (cm^ - 3)
    real :: mask
    real :: cor
    real :: tc0
    real :: bw

    real :: lambdar, lambdas, lambdag
    real :: rei_fac

    real :: rhow = 1.0e3, rhor = 1.0e3, rhos = 1.0e2, rhog = 4.0e2 ! density (kg / m^3)
    real :: n0r = 8.0e6, n0s = 3.0e6, n0g = 4.0e6 ! intercept parameters (m^ - 4)
    real :: alphar = 0.8, alphas = 0.25, alphag = 0.5 ! parameters in terminal equation in lin et al., 1983
    real :: gammar = 17.837789, gammas = 8.2850630, gammag = 11.631769 ! gamma values as a result of different alpha
    real, parameter :: rho_0 = 50.e-3

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

    if (present (cnvw)) then
        qmw = qmw + cnvw
    endif
    if (present (cnvi)) then
        qmi = qmi + cnvi
    endif
    if (present (cnvc)) then
        cld = cnvc + (1 - cnvc) * cld
    endif

    if (liq_ice_combine) then
        do k = ks, ke
            do i = is, ie

                ! frozen condensates:
                ! cloud ice treated as snow above freezing and graupel exists
                if (t (i, k) > tice) then
                    qms (i, k) = qmi (i, k) + qms (i, k)
                    qmi (i, k) = 0.
                else
                    qmi (i, k) = qmi (i, k) + qms (i, k)
                    if (qmi (i, k) .gt. qi0_max) then
                        qms (i, k) = qmi (i, k) - qi0_max + qmg (i, k)
                        qmi (i, k) = qi0_max
                    else
                        qms (i, k) = qmg (i, k)
                    endif
                    qmg (i, k) = 0. ! treating all graupel as "snow"
                endif
            enddo
        enddo
    else
        ! treating snow as ice, graupel as snow
        ! qmi (:, :) = qmi (:, :) + qms (:, :)
        ! qms (:, :) = qmg (:, :)
        ! qmg (:, :) = 0. ! treating all graupel as "snow"
        do k = ks, ke
            do i = is, ie
                ! step - 1: combine cloud ice & snow
                qmi (i, k) = qmi (i, k) + qms (i, k)
                ! step - 2: auto - convert cloud ice if > qi0_max
                qms (i, k) = qmi (i, k) - qi0_max
                if (qms (i, k) .gt. 0.) then
                    qmi (i, k) = qi0_max
                else
                    qms (i, k) = 0.0
                endif
            enddo
        enddo
    endif

    ! liquid condensates:

    ! sjl: 20180825
#ifdef COMBINE_QR
    do k = ks, ke
        do i = is, ie
            ! step - 1: combine cloud water & rain
            qmw (i, k) = qmw (i, k) + qmr (i, k)
            ! step - 2: auto - convert cloud wat if > ql0_max
            qmr (i, k) = qmw (i, k) - ql0_max
            if (qmr (i, k) .gt. 0.) then
                qmw (i, k) = ql0_max
            else
                qmr (i, k) = 0.0
            endif
        enddo
    enddo
#endif

    do k = ks, ke

        do i = is, ie

            qmw (i, k) = max (qmw (i, k), 0.0)
            qmi (i, k) = max (qmi (i, k), 0.0)
            qmr (i, k) = max (qmr (i, k), 0.0)
            qms (i, k) = max (qms (i, k), 0.0)
            qmg (i, k) = max (qmg (i, k), 0.0)

            cld (i, k) = min (max (cld (i, k), 0.0), 1.0)

            mask = min (max (lsm (i), 0.0), 2.0)

            dpg = abs (delp (i, k)) / grav
            ! sjl:
            ! rho = p (i, k) / (rdgas * t (i, k) * (1. + zvir * qv)) ! needs qv
            rho = p (i, k) / (rdgas * t (i, k))
            ! use rho = dpg / delz ! needs delz

            tc0 = t (i, k) - tice

            if (rewflag .eq. 1) then

                ! -----------------------------------------------------------------------
                ! cloud water (martin et al., 1994)
                ! -----------------------------------------------------------------------

                ccnw = 0.80 * (- 1.15e-3 * (ccn_o ** 2) + 0.963 * ccn_o + 5.30) * abs (mask - 1.0) + &
                0.67 * (- 2.10e-4 * (ccn_l ** 2) + 0.568 * ccn_l - 27.9) * (1.0 - abs (mask - 1.0))

                if (qmw (i, k) .gt. qmin) then
                    qcw (i, k) = betaw * dpg * qmw (i, k) * 1.0e3
                    rew (i, k) = exp (1.0 / 3.0 * log ((3.0 * qmw (i, k) * rho) / (4.0 * pi * rhow * ccnw))) * 1.0e4
                    rew (i, k) = max (rewmin, min (rewmax, rew (i, k)))
                else
                    qcw (i, k) = 0.0
                    rew (i, k) = rewmin
                endif

            endif

            if (rewflag .eq. 2) then

                ! -----------------------------------------------------------------------
                ! cloud water (martin et al., 1994, gfdl revision)
                ! -----------------------------------------------------------------------

                ccnw = 1.077 * ccn_o * abs (mask - 1.0) + 1.143 * ccn_l * (1.0 - abs (mask - 1.0))

                if (qmw (i, k) .gt. qmin) then
                    qcw (i, k) = betaw * dpg * qmw (i, k) * 1.0e3
                    rew (i, k) = exp (1.0 / 3.0 * log ((3.0 * qmw (i, k) / cld (i, k) * rho) / (4.0 * pi * rhow * ccnw))) * 1.0e4
                    rew (i, k) = max (rewmin, min (rewmax, rew (i, k)))
                else
                    qcw (i, k) = 0.0
                    rew (i, k) = rewmin
                endif

            endif

            if (rewflag .eq. 3) then

                ! -----------------------------------------------------------------------
                ! cloud water (kiehl et al., 1994)
                ! -----------------------------------------------------------------------

                if (qmw (i, k) .gt. qmin) then
                    qcw (i, k) = betaw * dpg * qmw (i, k) * 1.0e3
                    rew (i, k) = 14.0 * abs (mask - 1.0) + &
                     (8.0 + (14.0 - 8.0) * min (1.0, max (0.0, - tc0 / 30.0))) * (1.0 - abs (mask - 1.0))
                    rew (i, k) = rew (i, k) + (14.0 - rew (i, k)) * min (1.0, max (0.0, snowd (i) / 1000.0))
                    rew (i, k) = max (rewmin, min (rewmax, rew (i, k)))
                else
                    qcw (i, k) = 0.0
                    rew (i, k) = rewmin
                endif

            endif

            if (reiflag .eq. 1) then

                ! -----------------------------------------------------------------------
                ! cloud ice (heymsfield and mcfarquhar, 1996)
                ! -----------------------------------------------------------------------

                if (qmi (i, k) .gt. qmin) then
                    qci (i, k) = betai * dpg * qmi (i, k) * 1.0e3
                    ! sjl
                    ! rei_fac = log (1.0e3 * qmi (i, k) * rho)
                    rei_fac = log (1.0e3 * min (qi0_rei, qmi (i, k)) * rho)
                    if (tc0 .lt. - 50) then
                        ! rei (i, k) = beta / 9.917 * exp ((1. - 0.891) * rei_fac) * 1.0e3
                        rei (i, k) = beta / 9.917 * exp (0.109 * rei_fac) * 1.0e3
                    elseif (tc0 .lt. - 40) then
                        ! rei (i, k) = beta / 9.337 * exp ((1. - 0.920) * rei_fac) * 1.0e3
                        rei (i, k) = beta / 9.337 * exp (0.08 * rei_fac) * 1.0e3
                    elseif (tc0 .lt. - 30) then
                        ! rei (i, k) = beta / 9.208 * exp ((1. - 0.945) * rei_fac) * 1.0e3
                        rei (i, k) = beta / 9.208 * exp (0.055 * rei_fac) * 1.0e3
                    else
                        ! rei (i, k) = beta / 9.387 * exp ((1. - 0.969) * rei_fac) * 1.0e3
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
                ! cloud ice (donner et al., 1997)
                ! -----------------------------------------------------------------------

                if (qmi (i, k) .gt. qmin) then
                    qci (i, k) = betai * dpg * qmi (i, k) * 1.0e3
                    if (tc0 .le. - 55) then
                        rei (i, k) = 15.41627
                    elseif (tc0 .le. - 50) then
                        rei (i, k) = 16.60895
                    elseif (tc0 .le. - 45) then
                        rei (i, k) = 32.89967
                    elseif (tc0 .le. - 40) then
                        rei (i, k) = 35.29989
                    elseif (tc0 .le. - 35) then
                        rei (i, k) = 55.65818
                    elseif (tc0 .le. - 30) then
                        rei (i, k) = 85.19071
                    elseif (tc0 .le. - 25) then
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
                ! cloud ice (fu, 2007)
                ! -----------------------------------------------------------------------

                if (qmi (i, k) .gt. qmin) then
                    qci (i, k) = betai * dpg * qmi (i, k) * 1.0e3
                    ! use fu2007 form below - 10 c
                    if (tc0 > - 10) then
                        ! tc = - 10, rei = 40.6
                        rei (i, k) = 100.0 + tc0 * 5.94
                    else
                        rei (i, k) = 47.05 + tc0 * (0.6624 + 0.001741 * tc0)
                    endif
                    ! rei (i, k) = max (reimin, min (reimax, rei (i, k)))
                    rei (i, k) = max (reimin, rei (i, k))
                else
                    qci (i, k) = 0.0
                    rei (i, k) = reimin
                endif

            endif

            if (reiflag .eq. 4) then

                ! -----------------------------------------------------------------------
                ! cloud ice (kristjansson et al., 2000)
                ! -----------------------------------------------------------------------

                if (qmi (i, k) .gt. qmin) then
                    qci (i, k) = betai * dpg * qmi (i, k) * 1.0e3
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
                ! cloud ice (wyser, 1998)
                ! -----------------------------------------------------------------------

                if (qmi (i, k) .gt. qmin) then
                    qci (i, k) = betai * dpg * qmi (i, k) * 1.0e3
                    bw = - 2. + 1.e-3 * log10 (rho * qmi (i, k) / rho_0) * exp (1.5 * log (- min (- 1.e-6, tc0)))
                    rei (i, k) = 377.4 + bw * (203.3 + bw * (37.91 + 2.3696 * bw))
                    ! rei (i, k) = max (reimin, min (reimax, rei (i, k)))
                    rei (i, k) = max (reimin, rei (i, k))
                else
                    qci (i, k) = 0.0
                    rei (i, k) = reimin
                endif

            endif

            ! -----------------------------------------------------------------------
            ! rain (lin et al., 1983)
            ! -----------------------------------------------------------------------

            if (qmr (i, k) .gt. qmin) then
                qcr (i, k) = betar * dpg * qmr (i, k) * 1.0e3
                lambdar = exp (0.25 * log (pi * rhor * n0r / qmr (i, k) / rho))
                rer (i, k) = 0.5 * exp (log (gammar / 6) / alphar) / lambdar * 1.0e6
                rer (i, k) = max (rermin, min (rermax, rer (i, k)))
            else
                qcr (i, k) = 0.0
                rer (i, k) = rermin
            endif

            ! -----------------------------------------------------------------------
            ! snow (lin et al., 1983)
            ! -----------------------------------------------------------------------

            if (qms (i, k) .gt. qmin) then
                qcs (i, k) = betas * dpg * qms (i, k) * 1.0e3
                lambdas = exp (0.25 * log (pi * rhos * n0s / qms (i, k) / rho))
                res (i, k) = 0.5 * exp (log (gammas / 6) / alphas) / lambdas * 1.0e6
                res (i, k) = max (resmin, min (resmax, res (i, k)))
            else
                qcs (i, k) = 0.0
                res (i, k) = resmin
            endif

            ! -----------------------------------------------------------------------
            ! graupel (lin et al., 1983)
            ! -----------------------------------------------------------------------

            if (qmg (i, k) .gt. qmin) then
                qcg (i, k) = betag * dpg * qmg (i, k) * 1.0e3
                lambdag = exp (0.25 * log (pi * rhog * n0g / qmg (i, k) / rho))
                reg (i, k) = 0.5 * exp (log (gammag / 6) / alphag) / lambdag * 1.0e6
                reg (i, k) = max (regmin, min (regmax, reg (i, k)))
            else
                qcg (i, k) = 0.0
                reg (i, k) = regmin
            endif

        enddo

    enddo

end subroutine cloud_diagnosis

subroutine cloud_diagnosis_init (nlunit, input_nml_file, logunit, fn_nml)

    implicit none

    integer, intent (in) :: nlunit
    integer, intent (in) :: logunit

    character (len = 64), intent (in) :: fn_nml
    character (len = *), intent (in) :: input_nml_file (:)

    integer :: ios
    logical :: exists

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml = cloud_diagnosis_nml, iostat = ios)
#else
    inquire (file = trim (fn_nml), exist = exists)
    if (.not. exists) then
        write (6, *) 'cloud_diagnosis :: namelist file: ', trim (fn_nml), ' does not exist'
        stop
    else
        open (unit = nlunit, file = fn_nml, readonly, status = 'old', iostat = ios)
    endif
    rewind (nlunit)
    read (nlunit, nml = cloud_diagnosis_nml)
    close (nlunit)
#endif

end subroutine cloud_diagnosis_init

end module cloud_diagnosis_mod
