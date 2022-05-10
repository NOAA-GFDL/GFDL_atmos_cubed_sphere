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

module rad_ref_mod

    use constants_mod, only: grav, rdgas, pi => pi_8
    use fv_arrays_mod, only: fv_grid_bounds_type, r_grid
    use gfdl_mp_mod, only: do_hail, rhor, rhos, rhog, rhoh, rnzr, rnzs, rnzg, rnzh
    use gfdl_mp_mod, only: do_hail_inline => do_hail ! assuming same densities and numbers in both inline and traditional gfdl mp

contains

subroutine rad_ref (q, pt, delp, peln, delz, dbz, maxdbz, allmax, bd, &
        npz, ncnst, hydrostatic, zvir, in0r, in0s, in0g, iliqskin, do_inline_mp, &
        sphum, liq_wat, ice_wat, rainwat, snowwat, graupel, mp_top)

    ! code from mark stoelinga's dbzcalc.f from the rip package.
    ! currently just using values taken directly from that code, which is
    ! consistent for the mm5 reisner - 2 microphysics. from that file:

    ! this routine computes equivalent reflectivity factor (in dbz) at
    ! each model grid point. in calculating ze, the rip algorithm makes
    ! assumptions consistent with those made in an early version
    ! (ca. 1996) of the bulk mixed - phase microphysical scheme in the mm5
    ! model (i.e., the scheme known as "resiner - 2") . for each species:
    !
    ! 1. particles are assumed to be spheres of constant density. the
    ! densities of rain drops, snow particles, and graupel particles are
    ! taken to be rho_r = rho_l = 1000 kg m^ - 3, rho_s = 100 kg m^ - 3, and
    ! rho_g = 400 kg m^ - 3, respectively. (l refers to the density of
    ! liquid water.)
    !
    ! 2. the size distribution (in terms of the actual diameter of the
    ! particles, rather than the melted diameter or the equivalent solid
    ! ice sphere diameter) is assumed to follow an exponential
    ! distribution of the form n (d) = n_0 * exp (lambda * d) .
    !
    ! 3. if in0x = 0, the intercept parameter is assumed constant (as in
    ! early reisner - 2), with values of 8x10^6, 2x10^7, and 4x10^6 m^ - 4,
    ! for rain, snow, and graupel, respectively. various choices of
    ! in0x are available (or can be added) . currently, in0x = 1 gives the
    ! variable intercept for each species that is consistent with
    ! thompson, rasmussen, and manning (2004, monthly weather review,
    ! vol. 132, no. 2, pp. 519 - 542.)
    !
    ! 4. if iliqskin = 1, frozen particles that are at a temperature above
    ! freezing are assumed to scatter as a liquid particle.
    !
    ! more information on the derivation of simulated reflectivity in rip
    ! can be found in stoelinga (2005, unpublished write - up) . contact
    ! mark stoelinga (stoeling@atmos.washington.edu) for a copy.

    ! 22sep16: modifying to use the gfdl mp parameters. if doing so remember
    ! that the gfdl mp assumes a constant intercept (in0x = .false.)
    ! ferrier - aligo has an option for fixed slope (rather than fixed intercept) .
    ! thompson presumably is an extension of reisner mp.

    implicit none

    type (fv_grid_bounds_type), intent (in) :: bd

    logical, intent (in) :: hydrostatic, in0r, in0s, in0g, iliqskin, do_inline_mp

    integer, intent (in) :: npz, ncnst, mp_top
    integer, intent (in) :: sphum, liq_wat, ice_wat, rainwat, snowwat, graupel

    real, intent (in), dimension (bd%isd:bd%ied, bd%jsd:bd%jed, npz) :: pt, delp
    real, intent (in), dimension (bd%is:, bd%js:, 1:) :: delz
    real, intent (in), dimension (bd%isd:bd%ied, bd%jsd:bd%jed, npz, ncnst) :: q
    real, intent (in), dimension (bd%is :bd%ie, npz + 1, bd%js:bd%je) :: peln
    real, intent (out), dimension (bd%is :bd%ie, bd%js :bd%je, npz) :: dbz
    real, intent (out), dimension (bd%is :bd%ie, bd%js :bd%je) :: maxdbz

    real, intent (in) :: zvir
    real, intent (out) :: allmax

    ! parameters for constant intercepts (in0[rsg] = .false.)
    ! using gfdl mp values

    real (kind = r_grid), parameter :: vconr = 2503.23638966667
    real (kind = r_grid), parameter :: vcong = 87.2382675
    real (kind = r_grid), parameter :: vcons = 6.6280504
    real (kind = r_grid), parameter :: vconh = vcong
    real (kind = r_grid), parameter :: normr = 25132741228.7183
    real (kind = r_grid), parameter :: normg = 5026548245.74367
    real (kind = r_grid), parameter :: normh = pi * rhoh * rnzh
    real (kind = r_grid), parameter :: norms = 942477796.076938

    ! constants for variable intercepts
    ! will need to be changed based on mp scheme

    real, parameter :: r1 = 1.e-15
    real, parameter :: ron = 8.e6
    real, parameter :: ron2 = 1.e10
    real, parameter :: son = 2.e7
    real, parameter :: gon = 5.e7
    real, parameter :: ron_min = 8.e6
    real, parameter :: ron_qr0 = 0.00010
    real, parameter :: ron_delqr0 = 0.25 * ron_qr0
    real, parameter :: ron_const1r = (ron2 - ron_min) * 0.5
    real, parameter :: ron_const2r = (ron2 + ron_min) * 0.5

    ! other constants

    real, parameter :: gamma_seven = 720.
    real, parameter :: alpha = 0.224
    real (kind = r_grid), parameter :: factor_s = gamma_seven * 1.e18 * (1. / (pi * rhos)) ** 1.75 &
         * (rhos / rhor) ** 2 * alpha
    real, parameter :: qmin = 1.e-12
    real, parameter :: tice = 273.16

    ! double precision

    real (kind = r_grid), dimension (bd%is:bd%ie) :: rhoair, denfac, z_e
    real (kind = r_grid) :: qr1, qs1, qg1, t1, t2, t3, rwat, vtr, vtg, vts
    real (kind = r_grid) :: factorb_s, factorb_g
    real (kind = r_grid) :: temp_c, pres, sonv, gonv, ronv

    real :: rhogh, vcongh, normgh

    integer :: i, j, k
    integer :: is, ie, js, je

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je

    if (rainwat < 1) return

    dbz (:, :, 1:mp_top) = - 20.
    maxdbz (:, :) = - 20. ! minimum value
    allmax = - 20.

    if ((do_hail .and. .not. do_inline_mp) .or. (do_hail_inline .and. do_inline_mp)) then
        rhogh = rhoh
        vcongh = vconh
        normgh = normh
    else
        rhogh = rhog
        vcongh = vcong
        normgh = normg
    endif

    !$omp parallel do default (shared) private (rhoair, t1, t2, t3, denfac, vtr, vtg, vts, z_e)
    do k = mp_top + 1, npz
        do j = js, je
            if (hydrostatic) then
                do i = is, ie
                    rhoair (i) = delp (i, j, k) / ((peln (i, k + 1, j) - peln (i, k, j)) * &
                        rdgas * pt (i, j, k) * (1. + zvir * q (i, j, k, sphum)))
                    denfac (i) = sqrt (min (10., 1.2 / rhoair (i)))
                    z_e (i) = 0.
                enddo
            else
                do i = is, ie
                    rhoair (i) = - delp (i, j, k) / (grav * delz (i, j, k)) ! moist air density
                    denfac (i) = sqrt (min (10., 1.2 / rhoair (i)))
                    z_e (i) = 0.
                enddo
            endif
            if (rainwat > 0) then
                do i = is, ie
                    ! the following form vectorizes better & more consistent with gfdl_mp
                    ! sjl notes: marshall - palmer, dbz = 200 * precip ** 1.6, precip = 3.6e6 * t1 / rhor * vtr ! [mm / hr]
                    ! gfdl_mp terminal fall speeds are used
                    ! date modified 20170701
                    ! account for excessively high cloud water - > autoconvert (diag only) excess cloud water
                    t1 = rhoair (i) * max (qmin, q (i, j, k, rainwat) + dim (q (i, j, k, liq_wat), 1.0e-3))
                    vtr = max (1.e-3, vconr * denfac (i) * exp (0.2 * log (t1 / normr)))
                    z_e (i) = 200. * exp (1.6 * log (3.6e6 * t1 / rhor * vtr))
                    ! z_e = 200. * (exp (1.6 * log (3.6e6 * t1 / rhor * vtr)) + &
                    ! exp (1.6 * log (3.6e6 * t3 / rhogh * vtg)) + &
                    ! exp (1.6 * log (3.6e6 * t2 / rhos * vts)))
                enddo
            endif
            if (graupel > 0) then
                do i = is, ie
                    t3 = rhoair (i) * max (qmin, q (i, j, k, graupel))
                    vtg = max (1.e-3, vcongh * denfac (i) * exp (0.125 * log (t3 / normgh)))
                    z_e (i) = z_e (i) + 200. * exp (1.6 * log (3.6e6 * t3 / rhogh * vtg))
                enddo
            endif
            if (snowwat > 0) then
                do i = is, ie
                    t2 = rhoair (i) * max (qmin, q (i, j, k, snowwat))
                    ! vts = max (1.e-3, vcons * denfac * exp (0.0625 * log (t2 / norms)))
                    z_e (i) = z_e (i) + (factor_s / alpha) * t2 * exp (0.75 * log (t2 / rnzs))
                    ! z_e = 200. * (exp (1.6 * log (3.6e6 * t1 / rhor * vtr)) + &
                    ! exp (1.6 * log (3.6e6 * t3 / rhogh * vtg)) + &
                    ! exp (1.6 * log (3.6e6 * t2 / rhos * vts)))
                enddo
            endif
            do i = is, ie
                dbz (i, j, k) = 10. * log10 (max (0.01, z_e (i)))
            enddo
        enddo
    enddo

    !$omp parallel do default (shared)
    do j = js, je
        do k = mp_top + 1, npz
            do i = is, ie
                maxdbz (i, j) = max (dbz (i, j, k), maxdbz (i, j))
            enddo
        enddo
    enddo

    do j = js, je
        do i = is, ie
            allmax = max (maxdbz (i, j), allmax)
        enddo
    enddo

end subroutine rad_ref

end module rad_ref_mod
