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
! Intermediate Physics Interface
! Developer: Linjiong Zhou
! Last Update: 5/19/2022
! =======================================================================

module intermediate_phys_mod

    use constants_mod, only: rdgas, grav
    use fv_grid_utils_mod, only: cubed_to_latlon, update_dwinds_phys
    use fv_arrays_mod, only: fv_grid_type, fv_grid_bounds_type, inline_mp_type
    use mpp_domains_mod, only: domain2d, mpp_update_domains
    use tracer_manager_mod, only: get_tracer_index, get_tracer_names
    use field_manager_mod, only: model_atmos
    use gfdl_mp_mod, only: gfdl_mp_driver, fast_sat_adj, mtetw
    
    implicit none
    
    private

    real, parameter :: consv_min = 0.001

    public :: intermediate_phys

    ! -----------------------------------------------------------------------
    ! precision definition
    ! -----------------------------------------------------------------------
    
    integer, parameter :: r8 = 8 ! double precision
    
contains

subroutine intermediate_phys (is, ie, js, je, isd, ied, jsd, jed, km, npx, npy, nq, nwat, &
               c2l_ord, mdt, consv, akap, ptop, pfull, hs, te0_2d, u, v, w, pt, &
               delp, delz, q_con, cappa, q, pkz, r_vir, te_err, tw_err, inline_mp, &
               gridstruct, domain, bd, hydrostatic, do_adiabatic_init, &
               do_inline_mp, do_sat_adj, last_step, do_fast_phys, consv_checker, adj_mass_vmr)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------

    integer, intent (in) :: is, ie, js, je, isd, ied, jsd, jed, km, npx, npy, nq, c2l_ord, nwat

    logical, intent (in) :: hydrostatic, do_adiabatic_init, do_inline_mp, consv_checker
    logical, intent (in) :: do_sat_adj, last_step, do_fast_phys
    integer, intent (in) :: adj_mass_vmr

    real, intent (in) :: consv, mdt, akap, r_vir, ptop, te_err, tw_err

    real, intent (in), dimension (km) :: pfull

    real, intent (in), dimension (isd:ied, jsd:jed) :: hs

    real, intent (inout), dimension (is:, js:, 1:) :: delz
    
    real, intent (inout), dimension (isd:, jsd:, 1:) :: q_con, cappa, w
    
    real, intent (inout), dimension (is:ie, js:je) :: te0_2d

    real, intent (inout), dimension (isd:ied, jsd:jed, km) :: pt, delp

    real, intent (inout), dimension (isd:ied, jsd:jed, km, *) :: q

    real, intent (inout), dimension (isd:ied, jsd:jed+1, km) :: u

    real, intent (inout), dimension (isd:ied+1, jsd:jed, km) :: v

    real, intent (out), dimension (is:ie, js:je, km) :: pkz

    type (fv_grid_type), intent (in), target :: gridstruct

    type (fv_grid_bounds_type), intent (in) :: bd

    type (domain2d), intent (inout) :: domain

    type (inline_mp_type), intent (inout) :: inline_mp

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    logical, allocatable, dimension (:) :: conv_vmr_mmr

    integer :: i, j, k, m, kmp, sphum, liq_wat, ice_wat
    integer :: rainwat, snowwat, graupel, cld_amt, ccn_cm3, cin_cm3, aerosol

    real :: rrg

    real, dimension (is:ie) :: gsize

    real, dimension (is:ie, km) :: q2, q3, qliq, qsol, adj_vmr

    real, dimension (is:ie, km+1) :: phis, pe, peln

    real, dimension (isd:ied, jsd:jed, km) :: te, ua, va

    real, allocatable, dimension (:) :: wz

    real, allocatable, dimension (:,:) :: dz, wa

    real, allocatable, dimension (:,:,:) :: u_dt, v_dt, dp0, u0, v0
    
    real (kind = r8), allocatable, dimension (:) :: tz

    real (kind = r8), dimension (is:ie) :: te_b_beg, te_b_end, tw_b_beg, tw_b_end, dte, te_loss

    real (kind = r8), dimension (is:ie, 1:km) :: te_beg, te_end, tw_beg, tw_end
    
    character (len = 32) :: tracer_units, tracer_name

    sphum = get_tracer_index (model_atmos, 'sphum')
    liq_wat = get_tracer_index (model_atmos, 'liq_wat')
    ice_wat = get_tracer_index (model_atmos, 'ice_wat')
    rainwat = get_tracer_index (model_atmos, 'rainwat')
    snowwat = get_tracer_index (model_atmos, 'snowwat')
    graupel = get_tracer_index (model_atmos, 'graupel')
    cld_amt = get_tracer_index (model_atmos, 'cld_amt')
    ccn_cm3 = get_tracer_index (model_atmos, 'ccn_cm3')
    cin_cm3 = get_tracer_index (model_atmos, 'cin_cm3')
    aerosol = get_tracer_index (model_atmos, 'aerosol')

    rrg = - rdgas / grav

    ! time saving trick
    if (last_step) then
        kmp = 1
    else
        do k = 1, km
            kmp = k
            if (pfull (k) .gt. 50.E2) exit
        enddo
    endif

    ! decide which tracer needs adjustment
    if (.not. allocated (conv_vmr_mmr)) allocate (conv_vmr_mmr (nq))
    conv_vmr_mmr (:) = .false.
    if (adj_mass_vmr .gt. 0) then
        do m = 1, nq
            call get_tracer_names (model_atmos, m, name = tracer_name, units = tracer_units)
            if (trim (tracer_units) .eq. 'vmr') then
                conv_vmr_mmr (m) = .true.
            else
                conv_vmr_mmr (m) = .false.
            endif
        enddo
    endif

    !-----------------------------------------------------------------------
    ! Fast Saturation Adjustment >>>
    !-----------------------------------------------------------------------

    ! Note: pt at this stage is T_v
    if ((do_adiabatic_init .or. (.not. do_inline_mp) .or. do_sat_adj) .and. nwat .eq. 6) then

        allocate (dz (is:ie, kmp:km))

        allocate (tz (kmp:km))
        allocate (wz (kmp:km))

!$OMP parallel do default (none) shared (is, ie, js, je, isd, jsd, kmp, km, te, ptop, w, &
!$OMP                                    delp, hydrostatic, hs, pt, delz, rainwat, ua, va, &
!$OMP                                    liq_wat, ice_wat, snowwat, graupel, q_con, r_vir, &
!$OMP                                    sphum, pkz, last_step, consv, te0_2d, gridstruct, &
!$OMP                                    q, mdt, cld_amt, cappa, rrg, akap, ccn_cm3, &
!$OMP                                    cin_cm3, aerosol, do_sat_adj, adj_mass_vmr, &
!$OMP                                    conv_vmr_mmr, nq, consv_checker, te_err, tw_err) &
!$OMP                           private (q2, q3, gsize, dz, pe, peln, adj_vmr, qliq, qsol, &
!$OMP                                    tz, wz, dte, te_beg, tw_beg, te_b_beg, tw_b_beg, &
!$OMP                                    te_end, tw_end, te_b_end, tw_b_end, te_loss)

        do j = js, je

            ! grid size
            gsize (is:ie) = sqrt (gridstruct%area_64 (is:ie, j))

            ! aerosol
            if (aerosol .gt. 0) then
                q2 (is:ie, kmp:km) = q (is:ie, j, kmp:km, aerosol)
            elseif (ccn_cm3 .gt. 0) then
                q2 (is:ie, kmp:km) = q (is:ie, j, kmp:km, ccn_cm3)
            else
                q2 (is:ie, kmp:km) = 0.0
            endif
            if (cin_cm3 .gt. 0) then
                q3 (is:ie, kmp:km) = q (is:ie, j, kmp:km, cin_cm3)
            else
                q3 (is:ie, kmp:km) = 0.0
            endif
 
            ! total energy checker
            if (consv_checker) then
                qliq (is:ie, kmp:km) = q (is:ie, j, kmp:km, liq_wat) + q (is:ie, j, kmp:km, rainwat)
                qsol (is:ie, kmp:km) = q (is:ie, j, kmp:km, ice_wat) + q (is:ie, j, kmp:km, snowwat) + q (is:ie, j, kmp:km, graupel)
                te_beg (is:ie, kmp:km) = 0.0
                tw_beg (is:ie, kmp:km) = 0.0
                te_b_beg (is:ie) = 0.0
                tw_b_beg (is:ie) = 0.0
                do i = is, ie
                    tz (kmp:km) = pt (i, j, kmp:km) / ((1. + r_vir * q (i, j, kmp:km, sphum)) * (1. - (qliq (i, kmp:km) + qsol (i, kmp:km))))
                    dte (i) = 0.0
                    wz (kmp:km) = 0.0
                    ua (i, j, kmp:km) = 0.0
                    va (i, j, kmp:km) = 0.0
                    call mtetw (kmp, km, q (i, j, kmp:km, sphum), q (i, j, kmp:km, liq_wat), &
                        q (i, j, kmp:km, rainwat), q (i, j, kmp:km, ice_wat), q (i, j, kmp:km, snowwat), &
                        q (i, j, kmp:km, graupel), tz (kmp:km), ua (i, j, kmp:km), va (i, j, kmp:km), wz (kmp:km), &
                        delp (i, j, kmp:km), dte (i), 0.0, 0.0, 0.0, 0.0, 0.0, &
                        0.0, 0.0, 0.0, abs (mdt), te_beg (i, kmp:km), tw_beg (i, kmp:km), &
                        te_b_beg (i), tw_b_beg (i), .true., hydrostatic)
                enddo
            endif

            ! calculate pe, peln
            pe (is:ie, 1) = ptop
            peln (is:ie, 1) = log (ptop)
            do k = 2, km + 1
                pe (is:ie, k) = pe (is:ie, k-1) + delp (is:ie, j, k-1)
                peln (is:ie, k) = log (pe (is:ie, k))
            enddo

            ! layer thickness
            if (.not. hydrostatic) then
                dz (is:ie, kmp:km) = delz (is:ie, j, kmp:km)
            else
                dz (is:ie, kmp:km) = (peln (is:ie, kmp+1:km+1) - peln (is:ie, kmp:km)) * &
                    rrg * pt (is:ie, j, kmp:km)
            endif

            ! fast saturation adjustment
            call fast_sat_adj (abs (mdt), is, ie, kmp, km, hydrostatic, consv .gt. consv_min, &
                     adj_vmr (is:ie, kmp:km), te (is:ie, j, kmp:km), dte (is:ie), q (is:ie, j, kmp:km, sphum), &
                     q (is:ie, j, kmp:km, liq_wat), q (is:ie, j, kmp:km, rainwat), &
                     q (is:ie, j, kmp:km, ice_wat), q (is:ie, j, kmp:km, snowwat), &
                     q (is:ie, j, kmp:km, graupel), q (is:ie, j, kmp:km, cld_amt), &
                     q2 (is:ie, kmp:km), q3 (is:ie, kmp:km), hs (is:ie, j), &
                     dz (is:ie, kmp:km), pt (is:ie, j, kmp:km), delp (is:ie, j, kmp:km), &
#ifdef USE_COND
                     q_con (is:ie, j, kmp:km), &
#else
                     q_con (isd:, jsd, 1:), &
#endif
#ifdef MOIST_CAPPA
                     cappa (is:ie, j, kmp:km), &
#else
                     cappa (isd:, jsd, 1:), &
#endif
                     gsize, last_step, do_sat_adj)

            ! update non-microphyiscs tracers due to mass change
            if (adj_mass_vmr .gt. 0) then
                do m = 1, nq
                    if (conv_vmr_mmr (m)) then
                        q (is:ie, j, kmp:km, m) = q (is:ie, j, kmp:km, m) * adj_vmr (is:ie, kmp:km)
                    endif
                enddo
            endif

            ! update pkz
            if (.not. hydrostatic) then
#ifdef MOIST_CAPPA
                pkz (is:ie, j, kmp:km) = exp (cappa (is:ie, j, kmp:km) * &
                    log (rrg * delp (is:ie, j, kmp:km) / &
                    delz (is:ie, j, kmp:km) * pt (is:ie, j, kmp:km)))
#else
                pkz (is:ie, j, kmp:km) = exp (akap * log (rrg * delp (is:ie, j, kmp:km) / &
                    delz (is:ie, j, kmp:km) * pt (is:ie, j, kmp:km)))
#endif
            endif
 
            ! total energy checker
            if (consv_checker) then
                qliq (is:ie, kmp:km) = q (is:ie, j, kmp:km, liq_wat) + q (is:ie, j, kmp:km, rainwat)
                qsol (is:ie, kmp:km) = q (is:ie, j, kmp:km, ice_wat) + q (is:ie, j, kmp:km, snowwat) + q (is:ie, j, kmp:km, graupel)
                te_end (is:ie, kmp:km) = 0.0
                tw_end (is:ie, kmp:km) = 0.0
                te_b_end (is:ie) = 0.0
                tw_b_end (is:ie) = 0.0
                do i = is, ie
                    tz (kmp:km) = pt (i, j, kmp:km) / ((1. + r_vir * q (i, j, kmp:km, sphum)) * (1. - (qliq (i, kmp:km) + qsol (i, kmp:km))))
                    wz (kmp:km) = 0.0
                    ua (i, j, kmp:km) = 0.0
                    va (i, j, kmp:km) = 0.0
                    call mtetw (kmp, km, q (i, j, kmp:km, sphum), q (i, j, kmp:km, liq_wat), &
                        q (i, j, kmp:km, rainwat), q (i, j, kmp:km, ice_wat), q (i, j, kmp:km, snowwat), &
                        q (i, j, kmp:km, graupel), tz (kmp:km), ua (i, j, kmp:km), va (i, j, kmp:km), wz (kmp:km), &
                        delp (i, j, kmp:km), dte (i), 0.0, 0.0, 0.0, 0.0, 0.0, &
                        0.0, 0.0, 0.0, abs (mdt), te_end (i, kmp:km), tw_end (i, kmp:km), &
                        te_b_end (i), tw_b_end (i), .true., hydrostatic, te_loss (i))
                enddo
            endif

            ! add total energy change to te0_2d
            if (consv .gt. consv_min) then
                do i = is, ie
                    do k = kmp, km
                        te0_2d (i, j) = te0_2d (i, j) + te (i, j, k)
                    enddo
                enddo
            endif

            ! total energy checker
            if (consv_checker) then
                do i = is, ie
                    if (abs (sum (te_end (i, kmp:km)) + te_b_end (i) - sum (te_beg (i, kmp:km)) - te_b_beg (i)) / &
                         (sum (te_beg (i, kmp:km)) + te_b_beg (i)) .gt. te_err) then
                        print*, "FAST_SAT_ADJ TE: ", &
                            !(sum (te_beg (i, kmp:km)) + te_b_beg (i)), &
                            !(sum (te_end (i, kmp:km)) + te_b_end (i)), &
                            (sum (te_end (i, kmp:km)) + te_b_end (i) - sum (te_beg (i, kmp:km)) - te_b_beg (i)) / &
                            (sum (te_beg (i, kmp:km)) + te_b_beg (i))
                    endif
                    if (abs (sum (tw_end (i, kmp:km)) + tw_b_end (i) - sum (tw_beg (i, kmp:km)) - tw_b_beg (i)) / &
                         (sum (tw_beg (i, kmp:km)) + tw_b_beg (i)) .gt. tw_err) then
                        print*, "FAST_SAT_ADJ TW: ", &
                            !(sum (tw_beg (i, kmp:km)) + tw_b_beg (i)), &
                            !(sum (tw_end (i, kmp:km)) + tw_b_end (i)), &
                            (sum (tw_end (i, kmp:km)) + tw_b_end (i) - sum (tw_beg (i, kmp:km)) - tw_b_beg (i)) / &
                            (sum (tw_beg (i, kmp:km)) + tw_b_beg (i))
                    endif
                    !print*, "FAST_SAT_ADJ LOSS (%) : ", te_loss (i) / (sum (te_beg (i, kmp:km)) + te_b_beg (i)) * 100.0
                enddo
            endif

        enddo

        deallocate (dz)

        deallocate (tz)
        deallocate (wz)

    endif

    !-----------------------------------------------------------------------
    ! <<< Fast Saturation Adjustment
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! Inline GFDL MP >>>
    !-----------------------------------------------------------------------

    if ((.not. do_adiabatic_init) .and. do_inline_mp .and. nwat .eq. 6) then

        allocate (u_dt (isd:ied, jsd:jed, km))
        allocate (v_dt (isd:ied, jsd:jed, km))

        allocate (tz (kmp:km))
        allocate (wz (kmp:km))

        ! initialize wind tendencies
        do k = 1, km
            do j = jsd, jed
                do i = isd, ied
                    u_dt (i, j, k) = 0.
                    v_dt (i, j, k) = 0.
                enddo
            enddo
        enddo

        ! save D grid u and v
        if (consv .gt. consv_min) then
            allocate (u0 (isd:ied, jsd:jed+1, km))
            allocate (v0 (isd:ied+1, jsd:jed, km))
            u0 = u
            v0 = v
        endif

        ! D grid wind to A grid wind remap
        call cubed_to_latlon (u, v, ua, va, gridstruct, npx, npy, km, 1, gridstruct%grid_type, &
                 domain, gridstruct%bounded_domain, c2l_ord, bd)

        ! save delp
        if (consv .gt. consv_min) then
            allocate (dp0 (isd:ied, jsd:jed, km))
            dp0 = delp
        endif

        allocate (dz (is:ie, kmp:km))
        allocate (wa (is:ie, kmp:km))

!$OMP parallel do default (none) shared (is, ie, js, je, isd, jsd, kmp, km, ua, va, &
!$OMP                                    te, delp, hydrostatic, hs, pt, delz, ptop, &
!$OMP                                    rainwat, liq_wat, ice_wat, snowwat, graupel, q_con, &
!$OMP                                    sphum, w, pkz, last_step, consv, te0_2d, r_vir, &
!$OMP                                    gridstruct, q, mdt, cld_amt, cappa, rrg, akap, &
!$OMP                                    ccn_cm3, cin_cm3, inline_mp, do_inline_mp, consv_checker, &
!$OMP                                    u_dt, v_dt, aerosol, adj_mass_vmr, conv_vmr_mmr, nq, &
!$OMP                                    te_err, tw_err) &
!$OMP                           private (q2, q3, gsize, dz, wa, pe, peln, adj_vmr, qliq, qsol, &
!$OMP                                    tz, wz, dte, te_beg, tw_beg, te_b_beg, tw_b_beg, &
!$OMP                                    te_end, tw_end, te_b_end, tw_b_end, te_loss)

        do j = js, je

            ! grid size
            gsize (is:ie) = sqrt (gridstruct%area_64 (is:ie, j))

            ! aerosol
            if (aerosol .gt. 0) then
                q2 (is:ie, kmp:km) = q (is:ie, j, kmp:km, aerosol)
            elseif (ccn_cm3 .gt. 0) then
                q2 (is:ie, kmp:km) = q (is:ie, j, kmp:km, ccn_cm3)
            else
                q2 (is:ie, kmp:km) = 0.0
            endif
            if (cin_cm3 .gt. 0) then
                q3 (is:ie, kmp:km) = q (is:ie, j, kmp:km, cin_cm3)
            else
                q3 (is:ie, kmp:km) = 0.0
            endif
 
            ! note: ua and va are A-grid variables
            ! note: pt is virtual temperature at this point
            ! note: w is vertical velocity (m/s)
            ! note: delz is negative, delp is positive, delz doesn't change in constant volume situation
            ! note: hs is geopotential height (m^2/s^2)
            ! note: the unit of q2 or q3 is #/cm^3
            ! note: the unit of area is m^2
            ! note: the unit of prew, prer, prei, pres, preg is mm/day
            ! note: the unit of prefluxw, prefluxr, prefluxi, prefluxs, prefluxg is mm/day

            ! save ua, va for wind tendency calculation
            u_dt (is:ie, j, kmp:km) = ua (is:ie, j, kmp:km)
            v_dt (is:ie, j, kmp:km) = va (is:ie, j, kmp:km)

            ! initialize tendencies diagnostic
            if (allocated (inline_mp%liq_wat_dt)) inline_mp%liq_wat_dt (is:ie, j, kmp:km) = &
                inline_mp%liq_wat_dt (is:ie, j, kmp:km) - q (is:ie, j, kmp:km, liq_wat)
            if (allocated (inline_mp%ice_wat_dt)) inline_mp%ice_wat_dt (is:ie, j, kmp:km) = &
                inline_mp%ice_wat_dt (is:ie, j, kmp:km) - q (is:ie, j, kmp:km, ice_wat)
            if (allocated (inline_mp%qv_dt)) inline_mp%qv_dt (is:ie, j, kmp:km) = &
                inline_mp%qv_dt (is:ie, j, kmp:km) - q (is:ie, j, kmp:km, sphum)
            if (allocated (inline_mp%ql_dt)) inline_mp%ql_dt (is:ie, j, kmp:km) = &
                inline_mp%ql_dt (is:ie, j, kmp:km) - (q (is:ie, j, kmp:km, liq_wat) + &
                q (is:ie, j, kmp:km, rainwat))
            if (allocated (inline_mp%qi_dt)) inline_mp%qi_dt (is:ie, j, kmp:km) = &
                inline_mp%qi_dt (is:ie, j, kmp:km) - (q (is:ie, j, kmp:km, ice_wat) + &
                q (is:ie, j, kmp:km, snowwat) + q (is:ie, j, kmp:km, graupel))
            if (allocated (inline_mp%qr_dt)) inline_mp%qr_dt (is:ie, j, kmp:km) = &
                inline_mp%qr_dt (is:ie, j, kmp:km) - q (is:ie, j, kmp:km, rainwat)
            if (allocated (inline_mp%qs_dt)) inline_mp%qs_dt (is:ie, j, kmp:km) = &
                inline_mp%qs_dt (is:ie, j, kmp:km) - q (is:ie, j, kmp:km, snowwat)
            if (allocated (inline_mp%qg_dt)) inline_mp%qg_dt (is:ie, j, kmp:km) = &
                inline_mp%qg_dt (is:ie, j, kmp:km) - q (is:ie, j, kmp:km, graupel)
            if (allocated (inline_mp%t_dt)) inline_mp%t_dt (is:ie, j, kmp:km) = &
                inline_mp%t_dt (is:ie, j, kmp:km) - pt (is:ie, j, kmp:km)
            if (allocated (inline_mp%u_dt)) inline_mp%u_dt (is:ie, j, kmp:km) = &
                inline_mp%u_dt (is:ie, j, kmp:km) - ua (is:ie, j, kmp:km)
            if (allocated (inline_mp%v_dt)) inline_mp%v_dt (is:ie, j, kmp:km) = &
                inline_mp%v_dt (is:ie, j, kmp:km) - va (is:ie, j, kmp:km)

            ! total energy checker
            if (consv_checker) then
                qliq (is:ie, kmp:km) = q (is:ie, j, kmp:km, liq_wat) + q (is:ie, j, kmp:km, rainwat)
                qsol (is:ie, kmp:km) = q (is:ie, j, kmp:km, ice_wat) + q (is:ie, j, kmp:km, snowwat) + q (is:ie, j, kmp:km, graupel)
                te_beg (is:ie, kmp:km) = 0.0
                tw_beg (is:ie, kmp:km) = 0.0
                te_b_beg (is:ie) = 0.0
                tw_b_beg (is:ie) = 0.0
                do i = is, ie
                    tz (kmp:km) = pt (i, j, kmp:km) / ((1. + r_vir * q (i, j, kmp:km, sphum)) * (1. - (qliq (i, kmp:km) + qsol (i, kmp:km))))
                    if (hydrostatic) then
                        wz (kmp:km) = 0.0
                    else
                        wz (kmp:km) = w (i, j, kmp:km)
                    endif
                    dte (i) = 0.0
                    call mtetw (kmp, km, q (i, j, kmp:km, sphum), q (i, j, kmp:km, liq_wat), &
                        q (i, j, kmp:km, rainwat), q (i, j, kmp:km, ice_wat), q (i, j, kmp:km, snowwat), &
                        q (i, j, kmp:km, graupel), tz (kmp:km), ua (i, j, kmp:km), va (i, j, kmp:km), wz (kmp:km), &
                        delp (i, j, kmp:km), dte (i), 0.0, inline_mp%prew (i, j), &
                        inline_mp%prer (i, j), inline_mp%prei (i, j), inline_mp%pres (i, j), &
                        inline_mp%preg (i, j), 0.0, 0.0, abs (mdt), te_beg (i, kmp:km), tw_beg (i, kmp:km), &
                        te_b_beg (i), tw_b_beg (i), .true., hydrostatic)
                enddo
            endif

            ! calculate pe, peln
            pe (is:ie, 1) = ptop
            peln (is:ie, 1) = log (ptop)
            do k = 2, km + 1
                pe (is:ie, k) = pe (is:ie, k-1) + delp (is:ie, j, k-1)
                peln (is:ie, k) = log (pe (is:ie, k))
            enddo

            ! vertical velocity and layer thickness
            if (.not. hydrostatic) then
                wa (is:ie, kmp:km) = w (is:ie, j, kmp:km)
                dz (is:ie, kmp:km) = delz (is:ie, j, kmp:km)
            else
                dz (is:ie, kmp:km) = (peln (is:ie, kmp+1:km+1) - peln (is:ie, kmp:km)) * &
                    rrg * pt (is:ie, j, kmp:km)
            endif

            ! GFDL cloud microphysics main program
            call gfdl_mp_driver (q (is:ie, j, kmp:km, sphum), q (is:ie, j, kmp:km, liq_wat), &
                     q (is:ie, j, kmp:km, rainwat), q (is:ie, j, kmp:km, ice_wat), &
                     q (is:ie, j, kmp:km, snowwat), q (is:ie, j, kmp:km, graupel), &
                     q (is:ie, j, kmp:km, cld_amt), q2 (is:ie, kmp:km), &
                     q3 (is:ie, kmp:km), pt (is:ie, j, kmp:km), wa (is:ie, kmp:km), &
                     ua (is:ie, j, kmp:km), va (is:ie, j, kmp:km), dz (is:ie, kmp:km), &
                     delp (is:ie, j, kmp:km), gsize, abs (mdt), hs (is:ie, j), &
                     inline_mp%prew (is:ie, j), inline_mp%prer (is:ie, j), &
                     inline_mp%prei (is:ie, j), inline_mp%pres (is:ie, j), &
                     inline_mp%preg (is:ie, j), hydrostatic, is, ie, kmp, km, &
#ifdef USE_COND
                     q_con (is:ie, j, kmp:km), &
#else
                     q_con (isd:, jsd, 1:), &
#endif
#ifdef MOIST_CAPPA
                     cappa (is:ie, j, kmp:km), &
#else
                     cappa (isd:, jsd, 1:), &
#endif
                     consv .gt. consv_min, adj_vmr (is:ie, kmp:km), te (is:ie, j, kmp:km), dte (is:ie), &
                     inline_mp%prefluxw(is:ie, j, kmp:km), &
                     inline_mp%prefluxr(is:ie, j, kmp:km), inline_mp%prefluxi(is:ie, j, kmp:km), &
                     inline_mp%prefluxs(is:ie, j, kmp:km), inline_mp%prefluxg(is:ie, j, kmp:km), &
                     last_step, do_inline_mp)

            ! update non-microphyiscs tracers due to mass change
            if (adj_mass_vmr .gt. 0) then
                do m = 1, nq
                    if (conv_vmr_mmr (m)) then
                        q (is:ie, j, kmp:km, m) = q (is:ie, j, kmp:km, m) * adj_vmr (is:ie, kmp:km)
                    endif
                enddo
            endif

            ! update vertical velocity
            if (.not. hydrostatic) then
                w (is:ie, j, kmp:km) = wa (is:ie, kmp:km)
            endif

            ! compute wind tendency at A grid fori D grid wind update
            u_dt (is:ie, j, kmp:km) = (ua (is:ie, j, kmp:km) - u_dt (is:ie, j, kmp:km)) / abs (mdt)
            v_dt (is:ie, j, kmp:km) = (va (is:ie, j, kmp:km) - v_dt (is:ie, j, kmp:km)) / abs (mdt)

            ! update layer thickness
            if (.not. hydrostatic) then
                delz (is:ie, j, kmp:km) = dz (is:ie, kmp:km)
            endif

            ! tendencies diagnostic
            if (allocated (inline_mp%liq_wat_dt)) inline_mp%liq_wat_dt (is:ie, j, kmp:km) = &
                inline_mp%liq_wat_dt (is:ie, j, kmp:km) + q (is:ie, j, kmp:km, liq_wat)
            if (allocated (inline_mp%ice_wat_dt)) inline_mp%ice_wat_dt (is:ie, j, kmp:km) = &
                inline_mp%ice_wat_dt (is:ie, j, kmp:km) + q (is:ie, j, kmp:km, ice_wat)
            if (allocated (inline_mp%qv_dt)) inline_mp%qv_dt (is:ie, j, kmp:km) = &
                inline_mp%qv_dt (is:ie, j, kmp:km) + q (is:ie, j, kmp:km, sphum)
            if (allocated (inline_mp%ql_dt)) inline_mp%ql_dt (is:ie, j, kmp:km) = &
                inline_mp%ql_dt (is:ie, j, kmp:km) + (q (is:ie, j, kmp:km, liq_wat) + &
                q (is:ie, j, kmp:km, rainwat))
            if (allocated (inline_mp%qi_dt)) inline_mp%qi_dt (is:ie, j, kmp:km) = &
                inline_mp%qi_dt (is:ie, j, kmp:km) + (q (is:ie, j, kmp:km, ice_wat) + &
                q (is:ie, j, kmp:km, snowwat) + q (is:ie, j, kmp:km, graupel))
            if (allocated (inline_mp%qr_dt)) inline_mp%qr_dt (is:ie, j, kmp:km) = &
                inline_mp%qr_dt (is:ie, j, kmp:km) + q (is:ie, j, kmp:km, rainwat)
            if (allocated (inline_mp%qs_dt)) inline_mp%qs_dt (is:ie, j, kmp:km) = &
                inline_mp%qs_dt (is:ie, j, kmp:km) + q (is:ie, j, kmp:km, snowwat)
            if (allocated (inline_mp%qg_dt)) inline_mp%qg_dt (is:ie, j, kmp:km) = &
                inline_mp%qg_dt (is:ie, j, kmp:km) + q (is:ie, j, kmp:km, graupel)
            if (allocated (inline_mp%t_dt)) inline_mp%t_dt (is:ie, j, kmp:km) = &
                inline_mp%t_dt (is:ie, j, kmp:km) + pt (is:ie, j, kmp:km)
            if (allocated (inline_mp%u_dt)) inline_mp%u_dt (is:ie, j, kmp:km) = &
                inline_mp%u_dt (is:ie, j, kmp:km) + ua (is:ie, j, kmp:km)
            if (allocated (inline_mp%v_dt)) inline_mp%v_dt (is:ie, j, kmp:km) = &
                inline_mp%v_dt (is:ie, j, kmp:km) + va (is:ie, j, kmp:km)

            ! update pkz
            if (.not. hydrostatic) then
#ifdef MOIST_CAPPA
                pkz (is:ie, j, kmp:km) = exp (cappa (is:ie, j, kmp:km) * &
                    log (rrg * delp (is:ie, j, kmp:km) / &
                    delz (is:ie, j, kmp:km) * pt (is:ie, j, kmp:km)))
#else
                pkz (is:ie, j, kmp:km) = exp (akap * log (rrg * delp (is:ie, j, kmp:km) / &
                    delz (is:ie, j, kmp:km) * pt (is:ie, j, kmp:km)))
#endif
            endif
 
            ! total energy checker
            if (consv_checker) then
                qliq (is:ie, kmp:km) = q (is:ie, j, kmp:km, liq_wat) + q (is:ie, j, kmp:km, rainwat)
                qsol (is:ie, kmp:km) = q (is:ie, j, kmp:km, ice_wat) + q (is:ie, j, kmp:km, snowwat) + q (is:ie, j, kmp:km, graupel)
                te_end (is:ie, kmp:km) = 0.0
                tw_end (is:ie, kmp:km) = 0.0
                te_b_end (is:ie) = 0.0
                tw_b_end (is:ie) = 0.0
                do i = is, ie
                    tz (kmp:km) = pt (i, j, kmp:km) / ((1. + r_vir * q (i, j, kmp:km, sphum)) * (1. - (qliq (i, kmp:km) + qsol (i, kmp:km))))
                    if (hydrostatic) then
                        wz (kmp:km) = 0.0
                    else
                        wz (kmp:km) = w (i, j, kmp:km)
                    endif
                    call mtetw (kmp, km, q (i, j, kmp:km, sphum), q (i, j, kmp:km, liq_wat), &
                        q (i, j, kmp:km, rainwat), q (i, j, kmp:km, ice_wat), q (i, j, kmp:km, snowwat), &
                        q (i, j, kmp:km, graupel), tz (kmp:km), ua (i, j, kmp:km), va (i, j, kmp:km), wz (kmp:km), &
                        delp (i, j, kmp:km), dte (i), 0.0, inline_mp%prew (i, j), &
                        inline_mp%prer (i, j), inline_mp%prei (i, j), inline_mp%pres (i, j), &
                        inline_mp%preg (i, j), 0.0, 0.0, abs (mdt), te_end (i, kmp:km), tw_end (i, kmp:km), &
                        te_b_end (i), tw_b_end (i), .true., hydrostatic, te_loss (i))
                enddo
            endif

            ! add total energy change to te0_2d
            if (consv .gt. consv_min) then
                do i = is, ie
                    do k = kmp, km
                        te0_2d (i, j) = te0_2d (i, j) + te (i, j, k)
                    enddo
                enddo
            endif

            ! total energy checker
            if (consv_checker) then
                do i = is, ie
                    if (abs (sum (te_end (i, kmp:km)) + te_b_end (i) - sum (te_beg (i, kmp:km)) - te_b_beg (i)) / &
                         (sum (te_beg (i, kmp:km)) + te_b_beg (i)) .gt. te_err) then
                        print*, "GFDL-MP-INTM TE: ", &
                            !(sum (te_beg (i, kmp:km)) + te_b_beg (i)), &
                            !(sum (te_end (i, kmp:km)) + te_b_end (i)), &
                            (sum (te_end (i, kmp:km)) + te_b_end (i) - sum (te_beg (i, kmp:km)) - te_b_beg (i)) / &
                            (sum (te_beg (i, kmp:km)) + te_b_beg (i))
                    endif
                    if (abs (sum (tw_end (i, kmp:km)) + tw_b_end (i) - sum (tw_beg (i, kmp:km)) - tw_b_beg (i)) / &
                         (sum (tw_beg (i, kmp:km)) + tw_b_beg (i)) .gt. tw_err) then
                        print*, "GFDL-MP-INTM TW: ", &
                            !(sum (tw_beg (i, kmp:km)) + tw_b_beg (i)), &
                            !(sum (tw_end (i, kmp:km)) + tw_b_end (i)), &
                            (sum (tw_end (i, kmp:km)) + tw_b_end (i) - sum (tw_beg (i, kmp:km)) - tw_b_beg (i)) / &
                            (sum (tw_beg (i, kmp:km)) + tw_b_beg (i))
                    endif
                    !print*, "GFDL-MP-INTM LOSS (%) : ", te_loss (i) / (sum (te_beg (i, kmp:km)) + te_b_beg (i)) * 100.0
                enddo
            endif

        enddo

        deallocate (dz)
        deallocate (wa)

        ! Note: (ua, va) are *lat-lon* wind tendenies on cell centers
        if ( gridstruct%square_domain ) then
            call mpp_update_domains (u_dt, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.false.)
            call mpp_update_domains (v_dt, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
        else
            call mpp_update_domains (u_dt, domain, complete=.false.)
            call mpp_update_domains (v_dt, domain, complete=.true.)
        endif
        ! update u_dt and v_dt in halo
        call mpp_update_domains (u_dt, v_dt, domain)

        ! update D grid wind
        call update_dwinds_phys (is, ie, js, je, isd, ied, jsd, jed, abs (mdt), u_dt, v_dt, u, v, &
                 gridstruct, npx, npy, km, domain)

        deallocate (u_dt)
        deallocate (v_dt)

        deallocate (tz)
        deallocate (wz)

        ! update dry total energy
        if (consv .gt. consv_min) then
!$OMP parallel do default (none) shared (is, ie, js, je, km, te0_2d, hydrostatic, delp, &
!$OMP                                    gridstruct, u, v, dp0, u0, v0, hs, delz, w) &
!$OMP                           private (phis)
            do j = js, je
                if (hydrostatic) then
                    do k = 1, km
                        do i = is, ie
                            te0_2d (i, j) = te0_2d (i, j) + delp (i, j, k) * &
                                (0.25 * gridstruct%rsin2 (i, j) * (u (i, j, k) ** 2 + &
                                u (i, j+1, k) ** 2 + v (i, j, k) ** 2 + v (i+1, j, k) ** 2 - &
                                (u (i, j, k) + u (i, j+1, k)) * (v (i, j, k) + v (i+1, j, k)) * &
                                gridstruct%cosa_s (i, j))) - dp0 (i, j, k) * &
                                (0.25 * gridstruct%rsin2 (i, j) * (u0 (i, j, k) ** 2 + &
                                u0 (i, j+1, k) ** 2 + v0 (i, j, k) ** 2 + v0 (i+1, j, k) ** 2 - &
                                (u0 (i, j, k) + u0 (i, j+1, k)) * (v0 (i, j, k) + v0 (i+1, j, k)) * &
                                gridstruct%cosa_s (i, j)))
                        enddo
                    enddo
                else
                    do i = is, ie
                        phis (i, km+1) = hs (i, j)
                    enddo
                    do k = km, 1, -1
                        do i = is, ie
                            phis (i, k) = phis (i, k+1) - grav * delz (i, j, k)
                        enddo
                    enddo
                    do k = 1, km
                        do i = is, ie
                            te0_2d (i, j) = te0_2d (i, j) + delp (i, j, k) * &
                                (0.5 * (phis (i, k) + phis (i, k+1) + w (i, j, k) ** 2 + 0.5 * &
                                gridstruct%rsin2 (i, j) * (u (i, j, k) ** 2 + u (i, j+1, k) ** 2 + &
                                v (i, j, k) ** 2 + v (i+1, j, k) ** 2 - (u (i, j, k) + &
                                u (i, j+1, k)) * (v (i, j, k) + v (i+1, j, k)) * &
                                gridstruct%cosa_s (i, j)))) - dp0 (i, j, k) * &
                                (0.5 * (phis (i, k) + phis (i, k+1) + w (i, j, k) ** 2 + &
                                0.5 * gridstruct%rsin2 (i, j) * (u0 (i, j, k) ** 2 + &
                                u0 (i, j+1, k) ** 2 + v0 (i, j, k) ** 2 + v0 (i+1, j, k) ** 2 - &
                                (u0 (i, j, k) + u0 (i, j+1, k)) * (v0 (i, j, k) + v0 (i+1, j, k)) * &
                                gridstruct%cosa_s (i, j))))
                        enddo
                    enddo
                endif
            enddo
        end if

        if (consv .gt. consv_min) then
            deallocate (u0)
            deallocate (v0)
            deallocate (dp0)
        endif

    endif

    !-----------------------------------------------------------------------
    ! <<< Inline GFDL MP
    !-----------------------------------------------------------------------

end subroutine intermediate_phys

end module intermediate_phys_mod
