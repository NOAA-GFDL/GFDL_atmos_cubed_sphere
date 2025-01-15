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

#ifdef OVERLOAD_R4
    use constantsR4_mod, only: rdgas, rvgas, grav, kappa, cp_air
#else
    use constants_mod, only: rdgas, rvgas, grav, kappa, cp_air
#endif
    use fv_grid_utils_mod, only: cubed_to_latlon, update_dwinds_phys
    use fv_arrays_mod, only: fv_grid_type, fv_grid_bounds_type, inline_mp_type
    use fv_arrays_mod, only: inline_pbl_type, inline_cnv_type, inline_gwd_type
    use fv_arrays_mod, only: fv_thermo_type
    use mpp_domains_mod, only: domain2d, mpp_update_domains
    use tracer_manager_mod, only: get_tracer_index, get_tracer_names
    use field_manager_mod, only: model_atmos
    use gfdl_mp_mod, only: gfdl_mp_driver, fast_sat_adj, c_liq, c_ice, cv_air, &
                           cv_vap, hlv, hlf, mtetw, tice
    use sa_tke_edmf_mod, only: sa_tke_edmf_sfc, sa_tke_edmf_pbl
    use sa_tke_edmf_new_mod, only: sa_tke_edmf_new_sfc, sa_tke_edmf_new_pbl
    use sa_sas_mod, only: sa_sas_deep, sa_sas_shal
    use sa_aamf_mod, only: sa_aamf_deep, sa_aamf_shal
    use sa_gwd_mod, only: sa_gwd_oro, sa_gwd_cnv
    use fv_timing_mod, only: timing_on, timing_off

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
               mdt, consv, akap, ptop, pfull, hs, te0_2d, u, v, w, omga, pt, &
               delp, delz, q_con, cappa, q, pkz, r_vir, te_err, tw_err, &
               inline_mp, inline_pbl, inline_cnv, inline_gwd, &
               gridstruct, thermostruct, domain, bd, hydrostatic, do_adiabatic_init, &
               do_inline_mp, do_inline_pbl, do_inline_cnv, do_inline_gwd, &
               do_sat_adj, last_step, do_fast_phys, consv_checker, adj_mass_vmr, &
               inline_pbl_flag, inline_cnv_flag)

    implicit none

    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------

    integer, intent (in) :: is, ie, js, je, isd, ied, jsd, jed, km, npx, npy, nq, nwat, inline_pbl_flag, inline_cnv_flag, adj_mass_vmr

    logical, intent (in) :: hydrostatic, do_adiabatic_init, do_inline_mp, do_inline_pbl, consv_checker
    logical, intent (in) :: do_inline_cnv, do_inline_gwd, do_sat_adj, last_step, do_fast_phys

    real, intent (in) :: consv, mdt, akap, r_vir, ptop, te_err, tw_err

    real, intent (in), dimension (km) :: pfull

    real, intent (in), dimension (isd:ied, jsd:jed) :: hs

    real, intent (in), dimension (isd:ied, jsd:jed, km) :: omga

    real, intent (inout), dimension (is:, js:, 1:) :: delz

    real, intent (inout), dimension (isd:, jsd:, 1:) :: q_con, cappa, w

    real, intent (inout), dimension (is:ie, js:je) :: te0_2d

    real, intent (inout), dimension (isd:ied, jsd:jed, km) :: pt, delp

    real, intent (inout), dimension (isd:ied, jsd:jed, km, *) :: q

    real, intent (inout), dimension (isd:ied, jsd:jed+1, km) :: u

    real, intent (inout), dimension (isd:ied+1, jsd:jed, km) :: v

    real, intent (out), dimension (is:ie, js:je, km) :: pkz

    type (fv_grid_type), intent (in), target :: gridstruct

    type (fv_thermo_type), intent (in), target :: thermostruct

    type (fv_grid_bounds_type), intent (in) :: bd

    type (domain2d), intent (inout) :: domain

    type (inline_mp_type), intent (inout) :: inline_mp

    type (inline_pbl_type), intent (inout) :: inline_pbl

    type (inline_cnv_type), intent (inout) :: inline_cnv

    type (inline_gwd_type), intent (inout) :: inline_gwd

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    logical :: safety_check = .true.

    logical, allocatable, dimension (:) :: conv_vmr_mmr

    integer :: i, j, k, m, n, kr, kmp, ncld, ntke, sphum, liq_wat, ice_wat, lsoil
    integer :: rainwat, snowwat, graupel, cld_amt, ccn_cm3, cin_cm3, aerosol
    integer :: ios, ntchm, ntchs
    integer :: k_con, k_cappa

    real :: rrg, tem

    real, dimension (is:ie) :: gsize, dqv, dql, dqi, dqr, dqs, dqg, ps_dt, q_liq, q_sol, c_moist, k1, k2

    real, dimension (is:ie, km) :: q2, q3, qliq, qsol, cvm, adj_vmr

    real, dimension (is:ie, km+1) :: phis, pe, peln

    real, dimension (isd:ied, jsd:jed, km) :: te, ua, va

    integer, allocatable, dimension (:) :: kinver, vegtype

    real, allocatable, dimension (:) :: rn, rb, u10m, v10m, sigmaf, stress, wind, tmp, wz, fscav

    real, allocatable, dimension (:) :: dtsfc, dqvsfc, dqlsfc, dqisfc, dqrsfc, dqssfc, dqgsfc

    real, allocatable, dimension (:,:) :: dz, zm, zi, wa, dp, pm, pi, pmk, pik, qv, ql, qr, ta, uu, vv, ww, radh

    real, allocatable, dimension (:,:,:) :: u_dt, v_dt, dp0, u0, v0, qa

    real (kind = r8), allocatable, dimension (:) :: tz

    real (kind = r8), dimension (is:ie) :: te_b_beg, te_b_end, tw_b_beg, tw_b_end, dte, te_loss

    real (kind = r8), dimension (is:ie, 1:km) :: te_beg, te_end, tw_beg, tw_end, te8, dte8

    character (len = 32) :: tracer_units, tracer_name

    character (len = 20) :: fscav_aero (20) = 'default'

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
    ntke = get_tracer_index (model_atmos, 'sgs_tke')

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

    if (thermostruct%use_cond) then
       k_con = kmp
    else
       k_con = 1
    endif
    if (thermostruct%moist_kappa) then
       k_cappa = kmp
    else
       k_cappa = 1
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
!$OMP                                    conv_vmr_mmr, nq, consv_checker, te_err, tw_err, &
!$OMP                                    inline_mp,k_con,k_cappa,thermostruct) &
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
                     q_con (is:ie, j, k_con:), cappa (is:ie, j, k_cappa:), &
                     gsize, inline_mp%mppcw (is:ie, j), inline_mp%mppew (is:ie, j), inline_mp%mppe1 (is:ie, j), &
                     inline_mp%mpper (is:ie, j), inline_mp%mppdi (is:ie, j), inline_mp%mppd1 (is:ie, j), &
                     inline_mp%mppds (is:ie, j), inline_mp%mppdg (is:ie, j), inline_mp%mppsi (is:ie, j), &
                     inline_mp%mpps1 (is:ie, j), inline_mp%mppss (is:ie, j), inline_mp%mppsg (is:ie, j), &
                     inline_mp%mppfw (is:ie, j), inline_mp%mppfr (is:ie, j), inline_mp%mppmi (is:ie, j), &
                     inline_mp%mppms (is:ie, j), inline_mp%mppmg (is:ie, j), inline_mp%mppm1 (is:ie, j), &
                     inline_mp%mppm2 (is:ie, j), inline_mp%mppm3 (is:ie, j), inline_mp%mppar (is:ie, j), &
                     inline_mp%mppas (is:ie, j), inline_mp%mppag (is:ie, j), inline_mp%mpprs (is:ie, j), &
                     inline_mp%mpprg (is:ie, j), inline_mp%mppxr (is:ie, j), inline_mp%mppxs (is:ie, j), &
                     inline_mp%mppxg (is:ie, j), last_step, do_sat_adj, &
                     thermostruct%use_cond, thermostruct%moist_kappa)

            ! update non-microphyiscs tracers due to mass change
            if (adj_mass_vmr .gt. 0) then
                do m = 1, nq
                    if (conv_vmr_mmr (m)) then
                        q (is:ie, j, kmp:km, m) = q (is:ie, j, kmp:km, m) * adj_vmr (is:ie, kmp:km)
                    endif
                enddo
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

            ! update pkz
            if (.not. hydrostatic) then
               if (thermostruct%moist_kappa) then
                  pkz (is:ie, j, kmp:km) = exp (cappa (is:ie, j, kmp:km) * &
                       log (rrg * delp (is:ie, j, kmp:km) / &
                       delz (is:ie, j, kmp:km) * pt (is:ie, j, kmp:km)))
               else
                  pkz (is:ie, j, kmp:km) = exp (akap * log (rrg * delp (is:ie, j, kmp:km) / &
                       delz (is:ie, j, kmp:km) * pt (is:ie, j, kmp:km)))
               endif
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
                    !if (abs (sum (te_end (i, kmp:km)) + te_b_end (i) - sum (te_beg (i, kmp:km)) - te_b_beg (i)) / &
                    !     (sum (te_beg (i, kmp:km)) + te_b_beg (i)) .gt. te_err) then
                    !    print*, "FAST_SAT_ADJ TE: ", &
                    !        !(sum (te_beg (i, kmp:km)) + te_b_beg (i)), &
                    !        !(sum (te_end (i, kmp:km)) + te_b_end (i)), &
                    !        (sum (te_end (i, kmp:km)) + te_b_end (i) - sum (te_beg (i, kmp:km)) - te_b_beg (i)) / &
                    !        (sum (te_beg (i, kmp:km)) + te_b_beg (i))
                    !endif
                    inline_mp%fast_te_a_chg (i, j) = sum (te_end (i, :)) - sum (te_beg (i, :))
                    inline_mp%fast_te_b_chg (i, j) = te_b_end (i) - te_b_beg (i)
                    !if (abs (sum (tw_end (i, kmp:km)) + tw_b_end (i) - sum (tw_beg (i, kmp:km)) - tw_b_beg (i)) / &
                    !     (sum (tw_beg (i, kmp:km)) + tw_b_beg (i)) .gt. tw_err) then
                    !    print*, "FAST_SAT_ADJ TW: ", &
                    !        !(sum (tw_beg (i, kmp:km)) + tw_b_beg (i)), &
                    !        !(sum (tw_end (i, kmp:km)) + tw_b_end (i)), &
                    !        (sum (tw_end (i, kmp:km)) + tw_b_end (i) - sum (tw_beg (i, kmp:km)) - tw_b_beg (i)) / &
                    !        (sum (tw_beg (i, kmp:km)) + tw_b_beg (i))
                    !endif
                    inline_mp%fast_tw_a_chg (i, j) = sum (tw_end (i, :)) - sum (tw_beg (i, :))
                    inline_mp%fast_tw_b_chg (i, j) = tw_b_end (i) - tw_b_beg (i)
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
    ! Inline Planetary Boundary Layer >>>
    !-----------------------------------------------------------------------

    if ((.not. do_adiabatic_init) .and. do_inline_pbl .and. (.not. do_fast_phys)) then

        allocate (kinver (is:ie))

        allocate (dz (is:ie, 1:km))
        allocate (zm (is:ie, 1:km))
        allocate (zi (is:ie, 1:km+1))
        allocate (dp (is:ie, 1:km))
        allocate (pm (is:ie, 1:km))
        allocate (pi (is:ie, 1:km+1))
        allocate (pmk (is:ie, 1:km))
        allocate (pik (is:ie, 1:km+1))

        allocate (ta (is:ie, 1:km))
        allocate (uu (is:ie, 1:km))
        allocate (vv (is:ie, 1:km))
        allocate (qa (is:ie, 1:km, 1:nq))

        allocate (radh (is:ie, 1:km))
        allocate (rb (is:ie))
        allocate (u10m (is:ie))
        allocate (v10m (is:ie))
        allocate (stress (is:ie))
        allocate (wind (is:ie))
        allocate (sigmaf (is:ie))
        allocate (vegtype (is:ie))

        allocate (dtsfc (is:ie))
        allocate (dqvsfc (is:ie))
        allocate (dqlsfc (is:ie))
        allocate (dqisfc (is:ie))
        allocate (dqrsfc (is:ie))
        allocate (dqssfc (is:ie))
        allocate (dqgsfc (is:ie))

        allocate (tz (1:km))
        allocate (wz (1:km))

        allocate (u_dt (isd:ied, jsd:jed, km))
        allocate (v_dt (isd:ied, jsd:jed, km))

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
                 domain, gridstruct%bounded_domain, 4, bd)

        ! save delp
        if (consv .gt. consv_min) then
            allocate (dp0 (isd:ied, jsd:jed, km))
            dp0 = delp
        endif

!$OMP parallel do default (none) shared (is, ie, js, je, isd, jsd, km, nq, ua, va, w, &
!$OMP                                    te, delp, hydrostatic, hs, pt, delz, q_con, &
!$OMP                                    rainwat, liq_wat, ice_wat, snowwat, graupel, &
!$OMP                                    sphum, pkz, consv, te0_2d, gridstruct, q, &
!$OMP                                    mdt, cappa, rrg, akap, r_vir, u_dt, v_dt, &
!$OMP                                    ptop, ntke, inline_pbl, safety_check, nwat, &
!$OMP                                    adj_mass_vmr, conv_vmr_mmr, consv_checker, &
!$OMP                                    te_err, tw_err, inline_pbl_flag, thermostruct) &
!$OMP                           private (gsize, dz, zi, pi, pik, pmk, lsoil, pe, &
!$OMP                                    zm, dp, pm, ta, uu, vv, qliq, qsol, qa, adj_vmr, &
!$OMP                                    radh, rb, u10m, v10m, sigmaf, vegtype, q_liq, &
!$OMP                                    stress, wind, kinver, q_sol, c_moist, peln, &
!$OMP                                    cvm, kr, dqv, dql, dqi, dqr, dqs, dqg, ps_dt, &
!$OMP                                    tz, wz, dte, te_beg, tw_beg, te_b_beg, tw_b_beg, &
!$OMP                                    te_end, tw_end, te_b_end, tw_b_end, te_loss, &
!$OMP                                    dtsfc, dqvsfc, dqlsfc, dqisfc, dqrsfc, dqssfc, dqgsfc)

        do j = js, je

            ! grid size
            gsize (is:ie) = sqrt (gridstruct%area_64 (is:ie, j))

            ! save ua, va for wind tendency calculation
            u_dt (is:ie, j, 1:km) = ua (is:ie, j, 1:km)
            v_dt (is:ie, j, 1:km) = va (is:ie, j, 1:km)

            kinver = km
            lsoil = 4

            ! total energy before parameterization
            if (consv .gt. consv_min) then
                qliq = q (is:ie, j, 1:km, liq_wat) + q (is:ie, j, 1:km, rainwat)
                qsol = q (is:ie, j, 1:km, ice_wat) + q (is:ie, j, 1:km, snowwat) + q (is:ie, j, 1:km, graupel)
                cvm = (1 - (q (is:ie, j, 1:km, sphum) + qliq + qsol)) * cv_air + &
                    q (is:ie, j, 1:km, sphum) * cv_vap + qliq * c_liq + qsol * c_ice
                te (is:ie, j, 1:km) = - cvm * pt (is:ie, j, 1:km) / ((1. + r_vir * q (is:ie, j, 1:km, sphum)) * &
                    (1. - (qliq + qsol))) * delp (is:ie, j, 1:km)
            endif

            ! total energy checker
            if (consv_checker) then
                qliq = q (is:ie, j, 1:km, liq_wat) + q (is:ie, j, 1:km, rainwat)
                qsol = q (is:ie, j, 1:km, ice_wat) + q (is:ie, j, 1:km, snowwat) + q (is:ie, j, 1:km, graupel)
                te_beg (is:ie, 1:km) = 0.0
                tw_beg (is:ie, 1:km) = 0.0
                te_b_beg (is:ie) = 0.0
                tw_b_beg (is:ie) = 0.0
                do i = is, ie
                    tz = pt (i, j, 1:km) / ((1. + r_vir * q (i, j, 1:km, sphum)) * (1. - (qliq (i, 1:km) + qsol (i, 1:km))))
                    if (hydrostatic) then
                        wz = 0.0
                    else
                        wz = w (i, j, 1:km)
                    endif
                    dte (i) = 0.0
                    call mtetw (1, km, q (i, j, 1:km, sphum), q (i, j, 1:km, liq_wat), &
                        q (i, j, 1:km, rainwat), q (i, j, 1:km, ice_wat), q (i, j, 1:km, snowwat), &
                        q (i, j, 1:km, graupel), tz, ua (i, j, 1:km), va (i, j, 1:km), wz, &
                        delp (i, j, 1:km), dte (i), 0.0, 0.0, 0.0, 0.0, 0.0, &
                        0.0, 0.0, 0.0, abs (mdt), te_beg (i, 1:km), tw_beg (i, 1:km), &
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

            ! vertical index flip over
            zi (is:ie, 1) = 0.0
            pi (is:ie, 1) = pe (is:ie, km+1)
            pik (is:ie, 1) = exp (kappa * log (pi (is:ie, 1) * 1.e-5))
            inline_pbl%dtsfc (is:ie, j) = 0.0
            inline_pbl%dqsfc (is:ie, j) = 0.0
            dtsfc (is:ie) = 0.0
            dqvsfc (is:ie) = 0.0
            dqlsfc (is:ie) = 0.0
            dqisfc (is:ie) = 0.0
            dqrsfc (is:ie) = 0.0
            dqssfc (is:ie) = 0.0
            dqgsfc (is:ie) = 0.0
            inline_pbl%dusfc (is:ie, j) = 0.0
            inline_pbl%dvsfc (is:ie, j) = 0.0
            inline_pbl%dksfc (is:ie, j) = 0.0
            do k = 1, km
                kr = km - k + 1
                dp (is:ie, k) = delp (is:ie, j, kr)
                pi (is:ie, k+1) = pe (is:ie, kr)
                pik (is:ie, k+1) = exp (kappa * log (pi (is:ie, k+1) * 1.e-5))
                if (.not. hydrostatic) then
                    pm (is:ie, k) = dp (is:ie, k) / delz (is:ie, j, kr) * &
                        rrg * pt (is:ie, j, kr)
                    dz (is:ie, k) = delz (is:ie, j, kr)
                    ! ensure subgrid monotonicity of pressure
                    do i = is, ie
                        pm (i, k) = min (pm (i, k), pi (i, k) - 0.01 * pm (i, k))
                        pm (i, k) = max (pm (i, k), pi (i, k+1) + 0.01 * pm (i, k))
                    enddo
                else
                    pm (is:ie, k) = dp (is:ie, k) / (peln (is:ie, kr+1) - peln (is:ie, kr))
                    dz (is:ie, k) = (peln (is:ie, kr+1) - peln (is:ie, kr)) * &
                        rrg * pt (is:ie, j, kr)
                endif
                pmk (is:ie, k) = exp (kappa * log (pm (is:ie, k) * 1.e-5))
                zi (is:ie, k+1) = zi (is:ie, k) - dz (is:ie, k) * grav
                if (k .eq. 1) then
                    zm (is:ie, k) = - 0.5 * dz (is:ie, k) * grav
                else
                    zm (is:ie, k) = zm (is:ie, k-1) - 0.5 * (dz (is:ie, k-1) + dz (is:ie, k)) * grav
                endif
                q_liq = q (is:ie, j, kr, liq_wat) + q (is:ie, j, kr, rainwat)
                q_sol = q (is:ie, j, kr, ice_wat) + q (is:ie, j, kr, snowwat) + q (is:ie, j, kr, graupel)
                ta (is:ie, k) = pt (is:ie, j, kr) / ((1. + r_vir * q (is:ie, j, kr, sphum)) * &
                    (1. - (q_liq + q_sol)))
                uu (is:ie, k) = ua (is:ie, j, kr)
                vv (is:ie, k) = va (is:ie, j, kr)
                qa (is:ie, k, 1:nq) = q (is:ie, j, kr, 1:nq)
                radh (is:ie, k) = inline_pbl%radh (is:ie, j, kr)
                c_moist = (1 - (q (is:ie, j, kr, sphum) + q_liq + q_sol)) * cv_air + &
                    q (is:ie, j, kr, sphum) * cv_vap + q_liq * c_liq + q_sol * c_ice
                inline_pbl%dtsfc (is:ie, j) = inline_pbl%dtsfc (is:ie, j) - cp_air * ta (is:ie, k) * delp (is:ie, j, kr) / grav / abs (mdt)
                inline_pbl%dqsfc (is:ie, j) = inline_pbl%dqsfc (is:ie, j) - (hlv - rvgas * tice + (cv_vap - c_liq) * (ta (is:ie, k) - tice)) * q (is:ie, j, kr, sphum) * delp (is:ie, j, kr) / grav / abs (mdt)
                dtsfc (is:ie) = dtsfc (is:ie) - c_moist * ta (is:ie, k) * delp (is:ie, j, kr) / grav / abs (mdt)
                dqvsfc (is:ie) = dqvsfc (is:ie) - q (is:ie, j, kr, sphum) * delp (is:ie, j, kr) / grav / abs (mdt)
                dqlsfc (is:ie) = dqlsfc (is:ie) - q (is:ie, j, kr, liq_wat) * delp (is:ie, j, kr) / grav / abs (mdt)
                dqisfc (is:ie) = dqisfc (is:ie) - q (is:ie, j, kr, ice_wat) * delp (is:ie, j, kr) / grav / abs (mdt)
                dqrsfc (is:ie) = dqrsfc (is:ie) - q (is:ie, j, kr, rainwat) * delp (is:ie, j, kr) / grav / abs (mdt)
                dqssfc (is:ie) = dqssfc (is:ie) - q (is:ie, j, kr, snowwat) * delp (is:ie, j, kr) / grav / abs (mdt)
                dqgsfc (is:ie) = dqgsfc (is:ie) - q (is:ie, j, kr, graupel) * delp (is:ie, j, kr) / grav / abs (mdt)
                inline_pbl%dusfc (is:ie, j) = inline_pbl%dusfc (is:ie, j) - ua (is:ie, j, kr) * delp (is:ie, j, kr) / grav / abs (mdt)
                inline_pbl%dvsfc (is:ie, j) = inline_pbl%dvsfc (is:ie, j) - va (is:ie, j, kr) * delp (is:ie, j, kr) / grav / abs (mdt)
                inline_pbl%dksfc (is:ie, j) = inline_pbl%dksfc (is:ie, j) - 0.5 * (ua (is:ie, j, kr) ** 2 + va (is:ie, j, kr) ** 2 + w (is:ie, j, kr) ** 2) * delp (is:ie, j, kr) / grav / abs (mdt)
            enddo

            do i = is, ie
                sigmaf (i) = max (inline_pbl%vfrac (i, j), 0.01)
                vegtype (i) = int (inline_pbl%vtype (i, j) + 0.5)
            enddo

            ! check if pressure or height cross over
            if (safety_check) then
                do k = 1, km
                    do i = is, ie
                        if (k .lt. km) then
                            if (pm (i, k) .le. pm (i, k+1)) then
                                print*, "Warning: inline edmf pressure layer cross over", k, pm (i, k), pm (i, k+1)
                            endif
                            if (zm (i, k) .ge. zm (i, k+1)) then
                                print*, "Warning: inline edmf height layer cross over", k, zm (i, k), zm (i, k+1)
                            endif
                        endif
                        if (pi (i, k) .le. pi (i, k+1)) then
                            print*, "Warning: inline edmf pressure interface cross over", k, pi (i, k), pi (i, k+1)
                        endif
                        if (zi (i, k) .ge. zi (i, k+1)) then
                            print*, "Warning: inline edmf height interface cross over", k, zi (i, k), zi (i, k+1)
                        endif
                    enddo
                enddo
            endif

            if (inline_pbl_flag .eq. 1) &
                ! diagnose surface variables for PBL parameterization
                call sa_tke_edmf_sfc (ie-is+1, lsoil, pi (is:ie, 1), uu (is:ie, 1), &
                    vv (is:ie, 1), ta (is:ie, 1), qa (is:ie, 1, sphum), &
                    abs (mdt), inline_pbl%tsfc (is:ie, j), pm (is:ie, 1), &
                    pik (is:ie, 1) / pmk (is:ie, 1), inline_pbl%evap (is:ie, j), &
                    inline_pbl%hflx (is:ie, j), inline_pbl%ffmm (is:ie, j), &
                    inline_pbl%ffhh (is:ie, j), zm (is:ie, 1) / grav, &
                    inline_pbl%snowd (is:ie, j), inline_pbl%zorl (is:ie, j), inline_pbl%ztrl (is:ie, j), &
                    inline_pbl%lsm (is:ie, j), inline_pbl%uustar (is:ie, j), sigmaf, vegtype, &
                    inline_pbl%shdmax (is:ie, j), inline_pbl%sfcemis (is:ie, j), &
                    inline_pbl%dlwflx (is:ie, j), inline_pbl%sfcnsw (is:ie, j), &
                    inline_pbl%sfcdsw (is:ie, j), inline_pbl%srflag (is:ie, j), &
                    inline_pbl%hice (is:ie, j), inline_pbl%fice (is:ie, j), &
                    inline_pbl%tice (is:ie, j), inline_pbl%weasd (is:ie, j), &
                    inline_pbl%tprcp (is:ie, j), inline_pbl%stc (is:ie, j, :), &
                    inline_pbl%qsurf (is:ie, j), inline_pbl%cmm (is:ie, j), &
                    inline_pbl%chh (is:ie, j), inline_pbl%gflux (is:ie, j), &
                    inline_pbl%ep (is:ie, j), u10m_out = u10m, v10m_out = v10m, &
                    rb_out = rb, stress_out = stress, wind_out = wind)

                ! SA-TKE-EDMF main program
                call sa_tke_edmf_pbl (ie-is+1, km, nq, liq_wat, ice_wat, ntke, &
                    abs (mdt), uu, vv, ta, qa, gsize, inline_pbl%lsm (is:ie, j), &
                    radh, rb, inline_pbl%zorl (is:ie, j), u10m, v10m, &
                    inline_pbl%ffmm (is:ie, j), inline_pbl%ffhh (is:ie, j), &
                    inline_pbl%tsfc (is:ie, j), inline_pbl%hflx (is:ie, j), &
                    inline_pbl%evap (is:ie, j), stress, wind, kinver, &
                    pik (is:ie, 1), dp, pi, pm, pmk, zi, zm, &
                    inline_pbl%hpbl (is:ie, j), inline_pbl%kpbl (is:ie, j))
                    !inline_pbl%dusfc (is:ie, j), inline_pbl%dvsfc (is:ie, j), &
                    !inline_pbl%dtsfc (is:ie, j), inline_pbl%dqsfc (is:ie, j))

            if (inline_pbl_flag .eq. 2) &
                ! diagnose surface variables for PBL parameterization
                call sa_tke_edmf_new_sfc (ie-is+1, lsoil, pi (is:ie, 1), uu (is:ie, 1), &
                    vv (is:ie, 1), ta (is:ie, 1), qa (is:ie, 1, sphum), &
                    abs (mdt), inline_pbl%tsfc (is:ie, j), pm (is:ie, 1), &
                    pik (is:ie, 1) / pmk (is:ie, 1), inline_pbl%evap (is:ie, j), &
                    inline_pbl%hflx (is:ie, j), inline_pbl%ffmm (is:ie, j), &
                    inline_pbl%ffhh (is:ie, j), zm (is:ie, 1) / grav, &
                    inline_pbl%snowd (is:ie, j), inline_pbl%zorl (is:ie, j), inline_pbl%ztrl (is:ie, j), &
                    inline_pbl%lsm (is:ie, j), inline_pbl%uustar (is:ie, j), sigmaf, vegtype, &
                    inline_pbl%shdmax (is:ie, j), inline_pbl%sfcemis (is:ie, j), &
                    inline_pbl%dlwflx (is:ie, j), inline_pbl%sfcnsw (is:ie, j), &
                    inline_pbl%sfcdsw (is:ie, j), inline_pbl%srflag (is:ie, j), &
                    inline_pbl%hice (is:ie, j), inline_pbl%fice (is:ie, j), &
                    inline_pbl%tice (is:ie, j), inline_pbl%weasd (is:ie, j), &
                    inline_pbl%tprcp (is:ie, j), inline_pbl%stc (is:ie, j, :), &
                    inline_pbl%qsurf (is:ie, j), inline_pbl%cmm (is:ie, j), &
                    inline_pbl%chh (is:ie, j), inline_pbl%gflux (is:ie, j), &
                    inline_pbl%ep (is:ie, j), u10m_out = u10m, v10m_out = v10m, &
                    rb_out = rb, stress_out = stress, wind_out = wind)

                ! SA-TKE-EDMF main program
                call sa_tke_edmf_new_pbl (ie-is+1, km, nq, liq_wat, ice_wat, ntke, &
                    abs (mdt), uu, vv, ta, qa, gsize, inline_pbl%lsm (is:ie, j), &
                    radh, rb, sigmaf, inline_pbl%zorl (is:ie, j), u10m, v10m, &
                    inline_pbl%ffmm (is:ie, j), inline_pbl%ffhh (is:ie, j), &
                    inline_pbl%tsfc (is:ie, j), inline_pbl%hflx (is:ie, j), &
                    inline_pbl%evap (is:ie, j), stress, wind, kinver, &
                    pik (is:ie, 1), dp, pi, pm, pmk, zi, zm, &
                    inline_pbl%hpbl (is:ie, j), inline_pbl%kpbl (is:ie, j))
                    !inline_pbl%dusfc (is:ie, j), inline_pbl%dvsfc (is:ie, j), &
                    !inline_pbl%dtsfc (is:ie, j), inline_pbl%dqsfc (is:ie, j))

            ! update u, v, T, q, and delp, vertical index flip over
            do k = 1, km
                kr = km - k + 1
                q (is:ie, j, kr, nwat+1:nq) = qa (is:ie, k, nwat+1:nq)
                dqv = qa (is:ie, k, sphum) - q (is:ie, j, kr, sphum)
                dql = qa (is:ie, k, liq_wat) - q (is:ie, j, kr, liq_wat)
                dqi = qa (is:ie, k, ice_wat) - q (is:ie, j, kr, ice_wat)
                dqr = qa (is:ie, k, rainwat) - q (is:ie, j, kr, rainwat)
                dqs = qa (is:ie, k, snowwat) - q (is:ie, j, kr, snowwat)
                dqg = qa (is:ie, k, graupel) - q (is:ie, j, kr, graupel)
                ps_dt = 1 + dqv + dql + dqi + dqr + dqs + dqg
                adj_vmr (is:ie, kr) = (ps_dt - (qa (is:ie, k, sphum) + &
                    qa (is:ie, k, liq_wat) + qa (is:ie, k, ice_wat) + &
                    qa (is:ie, k, rainwat) + qa (is:ie, k, snowwat) + &
                    qa (is:ie, k, graupel))) / (1. - (qa (is:ie, k, sphum) + &
                    qa (is:ie, k, liq_wat) + qa (is:ie, k, ice_wat) + &
                    qa (is:ie, k, rainwat) + qa (is:ie, k, snowwat) + &
                    qa (is:ie, k, graupel))) / ps_dt
                q (is:ie, j, kr, sphum) = qa (is:ie, k, sphum) / ps_dt
                q (is:ie, j, kr, liq_wat) = qa (is:ie, k, liq_wat) / ps_dt
                q (is:ie, j, kr, ice_wat) = qa (is:ie, k, ice_wat) / ps_dt
                q (is:ie, j, kr, rainwat) = qa (is:ie, k, rainwat) / ps_dt
                q (is:ie, j, kr, snowwat) = qa (is:ie, k, snowwat) / ps_dt
                q (is:ie, j, kr, graupel) = qa (is:ie, k, graupel) / ps_dt
                delp (is:ie, j, kr) = delp (is:ie, j, kr) * ps_dt
                q_liq = q (is:ie, j, kr, liq_wat) + q (is:ie, j, kr, rainwat)
                q_sol = q (is:ie, j, kr, ice_wat) + q (is:ie, j, kr, snowwat) + q (is:ie, j, kr, graupel)
                if (thermostruct%use_cond) then
                    q_con (is:ie, j, kr) = q_liq + q_sol
                endif
                c_moist = (1 - (q (is:ie, j, kr, sphum) + q_liq + q_sol)) * cv_air + &
                    q (is:ie, j, kr, sphum) * cv_vap + q_liq * c_liq + q_sol * c_ice
                if (thermostruct%moist_kappa) then
                    cappa (is:ie, j, kr) = rdgas / (rdgas + c_moist / (1. + r_vir * q (is:ie, j, kr, sphum)))
                endif
                pt (is:ie, j, kr) = pt (is:ie, j, kr) + (ta (is:ie, k) * &
                    ((1. + r_vir * q (is:ie, j, kr, sphum)) * (1. - (q_liq + q_sol))) - &
                    pt (is:ie, j, kr)) * cp_air / c_moist
                ua (is:ie, j, kr) = uu (is:ie, k)
                va (is:ie, j, kr) = vv (is:ie, k)
                inline_pbl%dtsfc (is:ie, j) = inline_pbl%dtsfc (is:ie, j) + cp_air * ta (is:ie, k) * delp (is:ie, j, kr) / ps_dt / grav / abs (mdt)
                inline_pbl%dqsfc (is:ie, j) = inline_pbl%dqsfc (is:ie, j) + (hlv - rvgas * tice + (cv_vap - c_liq) * (ta (is:ie, k) - tice)) * q (is:ie, j, kr, sphum) * delp (is:ie, j, kr) / grav / abs (mdt)
                dtsfc (is:ie) = dtsfc (is:ie) + c_moist * (pt (is:ie, j, kr) / ((1. + r_vir * q (is:ie, j, kr, sphum)) * (1. - (q_liq + q_sol)))) * delp (is:ie, j, kr) / grav / abs (mdt)
                dqvsfc (is:ie) = dqvsfc (is:ie) + q (is:ie, j, kr, sphum) * delp (is:ie, j, kr) / grav / abs (mdt)
                dqlsfc (is:ie) = dqlsfc (is:ie) + q (is:ie, j, kr, liq_wat) * delp (is:ie, j, kr) / grav / abs (mdt)
                dqisfc (is:ie) = dqisfc (is:ie) + q (is:ie, j, kr, ice_wat) * delp (is:ie, j, kr) / grav / abs (mdt)
                dqrsfc (is:ie) = dqrsfc (is:ie) + q (is:ie, j, kr, rainwat) * delp (is:ie, j, kr) / grav / abs (mdt)
                dqssfc (is:ie) = dqssfc (is:ie) + q (is:ie, j, kr, snowwat) * delp (is:ie, j, kr) / grav / abs (mdt)
                dqgsfc (is:ie) = dqgsfc (is:ie) + q (is:ie, j, kr, graupel) * delp (is:ie, j, kr) / grav / abs (mdt)
                inline_pbl%dusfc (is:ie, j) = inline_pbl%dusfc (is:ie, j) + ua (is:ie, j, kr) * delp (is:ie, j, kr) / grav / abs (mdt)
                inline_pbl%dvsfc (is:ie, j) = inline_pbl%dvsfc (is:ie, j) + va (is:ie, j, kr) * delp (is:ie, j, kr) / grav / abs (mdt)
                inline_pbl%dksfc (is:ie, j) = inline_pbl%dksfc (is:ie, j) + 0.5 * (ua (is:ie, j, kr) ** 2 + va (is:ie, j, kr) ** 2 + w (is:ie, j, kr) ** 2) * delp (is:ie, j, kr) / grav / abs (mdt)
            enddo

            ! update non-microphyiscs tracers due to mass change
            if (adj_mass_vmr .gt. 0) then
                do m = 1, nq
                    if (conv_vmr_mmr (m)) then
                        q (is:ie, j, 1:km, m) = q (is:ie, j, 1:km, m) * adj_vmr (is:ie, 1:km)
                    endif
                enddo
            endif

            ! compute wind tendency at A grid fori D grid wind update
            u_dt (is:ie, j, 1:km) = (ua (is:ie, j, 1:km) - u_dt (is:ie, j, 1:km)) / abs (mdt)
            v_dt (is:ie, j, 1:km) = (va (is:ie, j, 1:km) - v_dt (is:ie, j, 1:km)) / abs (mdt)

            ! update pkz
            if (.not. hydrostatic) then
                if (thermostruct%moist_kappa) then
                    pkz (is:ie, j, 1:km) = exp (cappa (is:ie, j, 1:km) * &
                        log (rrg * delp (is:ie, j, 1:km) / &
                        delz (is:ie, j, 1:km) * pt (is:ie, j, 1:km)))
                else
                    pkz (is:ie, j, 1:km) = exp (akap * log (rrg * delp (is:ie, j, 1:km) / &
                        delz (is:ie, j, 1:km) * pt (is:ie, j, 1:km)))
                endif
            endif

            ! total energy checker
            if (consv_checker) then
                qliq = q (is:ie, j, 1:km, liq_wat) + q (is:ie, j, 1:km, rainwat)
                qsol = q (is:ie, j, 1:km, ice_wat) + q (is:ie, j, 1:km, snowwat) + q (is:ie, j, 1:km, graupel)
                te_end (is:ie, 1:km) = 0.0
                tw_end (is:ie, 1:km) = 0.0
                te_b_end (is:ie) = 0.0
                tw_b_end (is:ie) = 0.0
                do i = is, ie
                    tz = pt (i, j, 1:km) / ((1. + r_vir * q (i, j, 1:km, sphum)) * (1. - (qliq (i, 1:km) + qsol (i, 1:km))))
                    if (hydrostatic) then
                        wz = 0.0
                    else
                        wz = w (i, j, 1:km)
                    endif
                    dte (i) = 0.0
                    call mtetw (1, km, q (i, j, 1:km, sphum), q (i, j, 1:km, liq_wat), &
                        q (i, j, 1:km, rainwat), q (i, j, 1:km, ice_wat), q (i, j, 1:km, snowwat), &
                        q (i, j, 1:km, graupel), tz, ua (i, j, 1:km), va (i, j, 1:km), wz, &
                        delp (i, j, 1:km), dte (i), - dqvsfc (i) * 86400, - dqlsfc (i) * 86400, &
                        - dqrsfc (i) * 86400, - dqisfc (i) * 86400, - dqssfc (i) * 86400, - dqgsfc (i) * 86400, &
                        - dtsfc (i), - inline_pbl%dksfc (i, j), abs (mdt), te_end (i, 1:km), tw_end (i, 1:km), &
                        te_b_end (i), tw_b_end (i), .true., hydrostatic, te_loss (i))
                enddo
            endif

            ! total energy after parameterization, add total energy change to te0_2d
            if (consv .gt. consv_min) then
                qliq = q (is:ie, j, 1:km, liq_wat) + q (is:ie, j, 1:km, rainwat)
                qsol = q (is:ie, j, 1:km, ice_wat) + q (is:ie, j, 1:km, snowwat) + q (is:ie, j, 1:km, graupel)
                cvm = (1 - (q (is:ie, j, 1:km, sphum) + qliq + qsol)) * cv_air + &
                    q (is:ie, j, 1:km, sphum) * cv_vap + qliq * c_liq + qsol * c_ice
                te (is:ie, j, 1:km) = te (is:ie, j, 1:km) + &
                    cvm * pt (is:ie, j, 1:km) / ((1. + r_vir * q (is:ie, j, 1:km, sphum)) * &
                    (1. - (qliq + qsol))) * delp (is:ie, j, 1:km)
                do k = 1, km
                    te0_2d (is:ie, j) = te0_2d (is:ie, j) + te (is:ie, j, k)
                enddo
            endif

            ! total energy checker
            if (consv_checker) then
                do i = is, ie
                    !if (abs (sum (te_end (i, :)) + te_b_end (i) - sum (te_beg (i, :)) - te_b_beg (i)) / &
                    !     (sum (te_beg (i, :)) + te_b_beg (i)) .gt. te_err) then
                    !    print*, "PBL-INTM TE: ", &
                    !        !(sum (te_beg (i, :)) + te_b_beg (i)), &
                    !        !(sum (te_end (i, :)) + te_b_end (i)), &
                    !        (sum (te_end (i, :)) + te_b_end (i) - sum (te_beg (i, :)) - te_b_beg (i)) / &
                    !        (sum (te_beg (i, :)) + te_b_beg (i))
                    !endif
                    inline_pbl%intm_te_a_chg (i, j) = sum (te_end (i, :)) - sum (te_beg (i, :))
                    inline_pbl%intm_te_b_chg (i, j) = te_b_end (i) - te_b_beg (i)
                    !if (abs (sum (tw_end (i, :)) + tw_b_end (i) - sum (tw_beg (i, :)) - tw_b_beg (i)) / &
                    !     (sum (tw_beg (i, :)) + tw_b_beg (i)) .gt. tw_err) then
                    !    print*, "PBL-INTM TW: ", &
                    !        !(sum (tw_beg (i, :)) + tw_b_beg (i)), &
                    !        !(sum (tw_end (i, :)) + tw_b_end (i)), &
                    !        (sum (tw_end (i, :)) + tw_b_end (i) - sum (tw_beg (i, :)) - tw_b_beg (i)) / &
                    !        (sum (tw_beg (i, :)) + tw_b_beg (i))
                    !endif
                    inline_pbl%intm_tw_a_chg (i, j) = sum (tw_end (i, :)) - sum (tw_beg (i, :))
                    inline_pbl%intm_tw_b_chg (i, j) = tw_b_end (i) - tw_b_beg (i)
                    !print*, "PBL-INTM LOSS (%) : ", te_loss (i) / (sum (te_beg (i, :)) + te_b_beg (i)) * 100.0
                enddo
            endif

        enddo

        deallocate (kinver)

        deallocate (dz)
        deallocate (zm)
        deallocate (zi)
        deallocate (dp)
        deallocate (pm)
        deallocate (pi)
        deallocate (pmk)
        deallocate (pik)

        deallocate (ta)
        deallocate (uu)
        deallocate (vv)
        deallocate (qa)

        deallocate (radh)
        deallocate (rb)
        deallocate (u10m)
        deallocate (v10m)
        deallocate (stress)
        deallocate (wind)
        deallocate (sigmaf)
        deallocate (vegtype)

        deallocate (dtsfc)
        deallocate (dqvsfc)
        deallocate (dqlsfc)
        deallocate (dqisfc)
        deallocate (dqrsfc)
        deallocate (dqssfc)
        deallocate (dqgsfc)

        deallocate (tz)
        deallocate (wz)

        ! Note: (ua, va) are *lat-lon* wind tendenies on cell centers
        call timing_on('COMM_TOTAL')
        if ( gridstruct%square_domain ) then
            call mpp_update_domains (u_dt, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.false.)
            call mpp_update_domains (v_dt, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
        else
            call mpp_update_domains (u_dt, domain, complete=.false.)
            call mpp_update_domains (v_dt, domain, complete=.true.)
        endif
        call timing_off('COMM_TOTAL')

        ! update D grid wind
        call update_dwinds_phys (is, ie, js, je, isd, ied, jsd, jed, abs (mdt), u_dt, v_dt, u, v, &
                 gridstruct, npx, npy, km, domain)

        deallocate (u_dt)
        deallocate (v_dt)

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
    ! <<< Inline Planetary Boundary Layer
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! Inline Convection >>>
    !-----------------------------------------------------------------------

    if ((.not. do_adiabatic_init) .and. do_inline_cnv) then

        allocate (rn (is:ie))
        allocate (tmp (is:ie))

        allocate (dz (is:ie, 1:km))
        allocate (zm (is:ie, 1:km))
        allocate (dp (is:ie, 1:km))
        allocate (pm (is:ie, 1:km))
        allocate (pi (is:ie, 1:km+1))

        allocate (ta (is:ie, 1:km))
        allocate (qv (is:ie, 1:km))
        allocate (qr (is:ie, 1:km))
        allocate (uu (is:ie, 1:km))
        allocate (vv (is:ie, 1:km))
        allocate (ww (is:ie, 1:km))

        if (inline_cnv_flag .eq. 1) allocate (ql (is:ie, 1:km))
        if (inline_cnv_flag .eq. 2) allocate (qa (is:ie, 1:km, 1:nq))

        allocate (u_dt (isd:ied, jsd:jed, km))
        allocate (v_dt (isd:ied, jsd:jed, km))

        allocate (tz (1:km))
        allocate (wz (1:km))

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
                 domain, gridstruct%bounded_domain, 4, bd)

        ! save delp
        if (consv .gt. consv_min) then
            allocate (dp0 (isd:ied, jsd:jed, km))
            dp0 = delp
        endif

        if (inline_cnv_flag .eq. 2) then
            ntchm = 0 ! number of chemical tracers
            ntchs = get_tracer_index (model_atmos, 'so2') ! tracer index for first chemical tracer
            if (ntchs .gt. 0) then
                ntchm = get_tracer_index (model_atmos, 'pp10')
                if (ntchm .gt. 0) then
                    ntchm = ntchm - ntchs + 1
                endif
            endif
            ! setup aerosol scavenging factors
            allocate (fscav (ntchm))
            if (ntchm .gt. 0) then
                ! initialize to default
                fscav = 0.6
                n = get_tracer_index (model_atmos, 'seas1') - ntchs + 1
                if (n .gt. 0) fscav (n) = 1.0
                n = get_tracer_index (model_atmos, 'seas2') - ntchs + 1
                if (n .gt. 0) fscav (n) = 1.0
                n = get_tracer_index (model_atmos, 'seas3') - ntchs + 1
                if (n .gt. 0) fscav (n) = 1.0
                n = get_tracer_index (model_atmos, 'seas4') - ntchs + 1
                if (n .gt. 0) fscav (n) = 1.0
                n = get_tracer_index (model_atmos, 'seas5') - ntchs + 1
                if (n .gt. 0) fscav (n) = 1.0
                ! read factors from namelist
                do i = 1, size (fscav_aero)
                    j = index (fscav_aero (i), ":")
                    if (j .gt. 1) then
                        read (fscav_aero (i) (j+1:), *, iostat = ios) tem
                        if (ios .ne. 0) cycle
                        if (adjustl (fscav_aero (i) (:j-1)) .eq. "*") then
                            fscav = tem
                            exit
                        else
                            n = get_tracer_index (model_atmos, adjustl(fscav_aero (i) (:j-1))) - ntchs + 1
                            if (n .gt. 0) fscav (n) = tem
                        endif
                    endif
                enddo
            endif
        endif

!$OMP parallel do default (none) shared (is, ie, js, je, isd, jsd, km, ua, va, q_con, w, &
!$OMP                                    te, delp, hydrostatic, hs, pt, delz, omga, &
!$OMP                                    rainwat, liq_wat, ice_wat, snowwat, graupel, &
!$OMP                                    sphum, pkz, consv, te0_2d, gridstruct, q, &
!$OMP                                    mdt, cappa, rrg, akap, r_vir, inline_cnv, &
!$OMP                                    u_dt, v_dt, inline_pbl, safety_check, ptop, &
!$OMP                                    adj_mass_vmr, conv_vmr_mmr, nq, consv_checker, &
!$OMP                                    te_err, tw_err, inline_cnv_flag, fscav, ntchs, &
!$OMP                                    ntchm, ntke, nwat, thermostruct) &
!$OMP                           private (gsize, dz, pi, rn, tmp, q_liq, q_sol, pe, peln, qa, &
!$OMP                                    zm, dp, pm, qv, ql, qr, ta, uu, vv, ww, ncld, qliq, qsol, &
!$OMP                                    cvm, kr, dqv, dql, dqi, dqr, dqs, dqg, ps_dt, c_moist, &
!$OMP                                    adj_vmr, k1, k2, tz, wz, dte, te_beg, tw_beg, te_b_beg, tw_b_beg, &
!$OMP                                    te_end, tw_end, te_b_end, tw_b_end, te_loss, te8, dte8)

        do j = js, je

            ! grid size
            gsize (is:ie) = sqrt (gridstruct%area_64 (is:ie, j))

            ! save ua, va for wind tendency calculation
            u_dt (is:ie, j, 1:km) = ua (is:ie, j, 1:km)
            v_dt (is:ie, j, 1:km) = va (is:ie, j, 1:km)

            rn = 0.0
            ncld = 1
            inline_cnv%ktop (is:ie, j) = 1
            inline_cnv%kbot (is:ie, j) = km
            inline_cnv%kcnv (is:ie, j) = 0
            inline_cnv%cumabs (is:ie, j) = 0

            ! total energy before parameterization
            if (consv .gt. consv_min) then
                qliq = q (is:ie, j, 1:km, liq_wat) + q (is:ie, j, 1:km, rainwat)
                qsol = q (is:ie, j, 1:km, ice_wat) + q (is:ie, j, 1:km, snowwat) + q (is:ie, j, 1:km, graupel)
                cvm = (1 - (q (is:ie, j, 1:km, sphum) + qliq + qsol)) * cv_air + &
                    q (is:ie, j, 1:km, sphum) * cv_vap + qliq * c_liq + qsol * c_ice
                te (is:ie, j, 1:km) = - cvm * pt (is:ie, j, 1:km) / ((1. + r_vir * q (is:ie, j, 1:km, sphum)) * &
                    (1. - (qliq + qsol))) * delp (is:ie, j, 1:km)
            endif

            ! total energy checker
            if (consv_checker) then
                qliq = q (is:ie, j, 1:km, liq_wat) + q (is:ie, j, 1:km, rainwat)
                qsol = q (is:ie, j, 1:km, ice_wat) + q (is:ie, j, 1:km, snowwat) + q (is:ie, j, 1:km, graupel)
                te_beg (is:ie, 1:km) = 0.0
                tw_beg (is:ie, 1:km) = 0.0
                te_b_beg (is:ie) = 0.0
                tw_b_beg (is:ie) = 0.0
                do i = is, ie
                    tz = pt (i, j, 1:km) / ((1. + r_vir * q (i, j, 1:km, sphum)) * (1. - (qliq (i, 1:km) + qsol (i, 1:km))))
                    if (hydrostatic) then
                        wz = 0.0
                    else
                        wz = w (i, j, 1:km)
                    endif
                    dte (i) = 0.0
                    call mtetw (1, km, q (i, j, 1:km, sphum), q (i, j, 1:km, liq_wat), &
                        q (i, j, 1:km, rainwat), q (i, j, 1:km, ice_wat), q (i, j, 1:km, snowwat), &
                        q (i, j, 1:km, graupel), tz, ua (i, j, 1:km), va (i, j, 1:km), wz, &
                        delp (i, j, 1:km), dte (i), 0.0, 0.0, inline_cnv%prec (i, j) / abs (mdt) * 1.e3 * 86400, 0.0, 0.0, &
                        0.0, 0.0, 0.0, abs (mdt), te_beg (i, 1:km), tw_beg (i, 1:km), &
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

            ! vertical index flip over
            pi (is:ie, 1) = pe (is:ie, km+1)
            do k = 1, km
                kr = km - k + 1
                dp (is:ie, k) = delp (is:ie, j, kr)
                pi (is:ie, k+1) = pe (is:ie, kr)
                if (.not. hydrostatic) then
                    pm (is:ie, k) = dp (is:ie, k) / delz (is:ie, j, kr) * &
                        rrg * pt (is:ie, j, kr)
                    dz (is:ie, k) = delz (is:ie, j, kr)
                    ! ensure subgrid monotonicity of pressure
                    do i = is, ie
                        pm (i, k) = min (pm (i, k), pi (i, k) - 0.01 * pm (i, k))
                        pm (i, k) = max (pm (i, k), pi (i, k+1) + 0.01 * pm (i, k))
                    enddo
                else
                    pm (is:ie, k) = dp (is:ie, k) / (peln (is:ie, kr+1) - peln (is:ie, kr))
                    dz (is:ie, k) = (peln (is:ie, kr+1) - peln (is:ie, kr)) * &
                        rrg * pt (is:ie, j, kr)
                endif
                if (k .eq. 1) then
                    zm (is:ie, k) = - 0.5 * dz (is:ie, k) * grav
                else
                    zm (is:ie, k) = zm (is:ie, k-1) - 0.5 * (dz (is:ie, k-1) + dz (is:ie, k)) * grav
                endif
                qv (is:ie, k) = q (is:ie, j, kr, sphum)
                if (inline_cnv_flag .eq. 1) ql (is:ie, k) = q (is:ie, j, kr, liq_wat)
                q_liq = q (is:ie, j, kr, liq_wat) + q (is:ie, j, kr, rainwat)
                q_sol = q (is:ie, j, kr, ice_wat) + q (is:ie, j, kr, snowwat) + q (is:ie, j, kr, graupel)
                ta (is:ie, k) = pt (is:ie, j, kr) / ((1. + r_vir * q (is:ie, j, kr, sphum)) * &
                    (1. - (q_liq + q_sol)))
                uu (is:ie, k) = ua (is:ie, j, kr)
                vv (is:ie, k) = va (is:ie, j, kr)
                ww (is:ie, k) = omga (is:ie, j, kr)
                if (inline_cnv_flag .eq. 2) qa (is:ie, k, 1:nq) = q (is:ie, j, kr, 1:nq)
                c_moist = (1 - (q (is:ie, j, kr, sphum) + q_liq + q_sol)) * cv_air + &
                    q (is:ie, j, kr, sphum) * cv_vap + q_liq * c_liq + q_sol * c_ice
                te8 (is:ie, k) = (c_moist * ta (is:ie, k) + &
                    (hlv - rvgas * tice - (cv_vap - c_liq) * tice) * q (is:ie, j, kr, sphum) - &
                    (hlf - (c_liq - c_ice) * tice) * q_sol) * delp (is:ie, j, kr) / grav
                dte8 (is:ie, k) = 0.0
            enddo

            ! check if pressure or height cross over
            if (safety_check) then
                do k = 1, km
                    do i = is, ie
                        if (k .lt. km) then
                            if (pm (i, k) .le. pm (i, k+1)) then
                                print*, "Warning: inline sas pressure layer cross over", k, pm (i, k), pm (i, k+1)
                            endif
                        if (pi (i, k) .le. pi (i, k+1)) then
                            print*, "Warning: inline sas pressure interface cross over", k, pi (i, k), pi (i, k+1)
                        endif
                            if (zm (i, k) .ge. zm (i, k+1)) then
                                print*, "Warning: inline sas height layer cross over", k, zm (i, k), zm (i, k+1)
                            endif
                        endif
                    enddo
                enddo
            endif

            if (inline_cnv_flag .eq. 1) &
                ! SA-SAS deep convection main program
                call sa_sas_deep (ie-is+1, km, abs (mdt), dp, pm, pi (is:ie, 1), zm, ql, &
                    qv, ta, uu, vv, qr, rn, inline_cnv%kbot (is:ie, j), inline_cnv%ktop (is:ie, j), &
                    inline_cnv%kcnv (is:ie, j), inline_pbl%lsm (is:ie, j), gsize, ww, ncld)

            if (inline_cnv_flag .eq. 2) &
                ! SA-AAMF deep convection main program
                call sa_aamf_deep (ie-is+1, km, abs (mdt), ntchs, ntchm, liq_wat, ice_wat, ntke, nq - 2, dp, &
                    pm, pi (is:ie, 1), zm, qa, qv, ta, uu, vv, qr, fscav, rn, inline_cnv%kbot (is:ie, j), &
                    inline_cnv%ktop (is:ie, j), inline_cnv%kcnv (is:ie, j), inline_pbl%lsm (is:ie, j), gsize, ww, ncld)

            ! convective precipitation accumulation
            inline_cnv%prec (is:ie, j) = inline_cnv%prec (is:ie, j) + rn

            if (inline_cnv_flag .eq. 1) &
                ! SA-SAS shallow convection main program
                call sa_sas_shal (ie-is+1, km, abs (mdt), dp, pm, pi (is:ie, 1), zm, ql, &
                    qv, ta, uu, vv, qr, rn, inline_cnv%kbot (is:ie, j), inline_cnv%ktop (is:ie, j), &
                    inline_cnv%kcnv (is:ie, j), inline_pbl%lsm (is:ie, j), gsize, ww, ncld, &
                    inline_pbl%hpbl (is:ie, j))

            if (inline_cnv_flag .eq. 2) &
                ! SA-AAMF shallow convection main program
                call sa_aamf_shal (ie-is+1, km, abs (mdt), ntchs, ntchm, liq_wat, ice_wat, ntke, nq - 2, dp, &
                    pm, pi (is:ie, 1), zm, qa, qv, ta, uu, vv, qr, fscav, rn, inline_cnv%kbot (is:ie, j), &
                    inline_cnv%ktop (is:ie, j), inline_cnv%kcnv (is:ie, j), inline_pbl%lsm (is:ie, j), gsize, ww, ncld, &
                    inline_pbl%hpbl (is:ie, j))

            ! convective precipitation accumulation
            inline_cnv%prec (is:ie, j) = inline_cnv%prec (is:ie, j) + rn

            ! convective heating for convective gravity wave drag parameterization
            tmp (is:ie) = 0.0
            do k = 1, km
                kr = km - k + 1
                do i = is, ie
                    if (k .ge. inline_cnv%kbot (i, j) .and. k .le. inline_cnv%ktop (i, j)) then
                        inline_cnv%cumabs (i, j) = inline_cnv%cumabs (i, j) + &
                            (ta (i, k) - pt (i, j, kr)) * dp (i, k)
                        tmp (i) = tmp (i) + dp (i, k)
                    endif
                enddo
            enddo
            do i = is, ie
                if (tmp (i) .gt. 0.0) inline_cnv%cumabs (i, j) = inline_cnv%cumabs (i, j) / (abs (mdt) * tmp (i))
            enddo

            ! update u, v, T, q, and delp, vertical index flip over
            do k = 1, km
                kr = km - k + 1
                k1 = 0.5 * (ua (is:ie, j, kr) ** 2 + va (is:ie, j, kr) ** 2 + w (is:ie, j, kr) ** 2) * delp (is:ie, j, kr)
                if (inline_cnv_flag .eq. 1) then
                    dqv = qv (is:ie, k) - q (is:ie, j, kr, sphum)
                    dql = ql (is:ie, k) - q (is:ie, j, kr, liq_wat)
                    ps_dt = 1 + dqv + dql
                    adj_vmr (is:ie, kr) = (ps_dt - (qv (is:ie, k) + ql (is:ie, k) + &
                        q (is:ie, j, kr, ice_wat) + q (is:ie, j, kr, rainwat) + &
                        q (is:ie, j, kr, snowwat) + q (is:ie, j, kr, graupel))) / &
                        (1. - (qv (is:ie, k) + ql (is:ie, k) + q (is:ie, j, kr, ice_wat) + &
                        q (is:ie, j, kr, rainwat) + q (is:ie, j, kr, snowwat) + &
                        q (is:ie, j, kr, graupel))) / ps_dt
                    q (is:ie, j, kr, sphum) = qv (is:ie, k) / ps_dt
                    q (is:ie, j, kr, liq_wat) = ql (is:ie, k) / ps_dt
                    q (is:ie, j, kr, ice_wat) = q (is:ie, j, kr, ice_wat) / ps_dt
                    q (is:ie, j, kr, rainwat) = q (is:ie, j, kr, rainwat) / ps_dt
                    q (is:ie, j, kr, snowwat) = q (is:ie, j, kr, snowwat) / ps_dt
                    q (is:ie, j, kr, graupel) = q (is:ie, j, kr, graupel) / ps_dt
                endif
                if (inline_cnv_flag .eq. 2) then
                    q (is:ie, j, kr, nwat+1:nq) = qa (is:ie, k, nwat+1:nq)
                    dqv = qv (is:ie, k) - q (is:ie, j, kr, sphum)
                    dql = qa (is:ie, k, liq_wat) - q (is:ie, j, kr, liq_wat)
                    dqi = qa (is:ie, k, ice_wat) - q (is:ie, j, kr, ice_wat)
                    dqr = qa (is:ie, k, rainwat) - q (is:ie, j, kr, rainwat)
                    dqs = qa (is:ie, k, snowwat) - q (is:ie, j, kr, snowwat)
                    dqg = qa (is:ie, k, graupel) - q (is:ie, j, kr, graupel)
                    ps_dt = 1 + dqv + dql + dqi + dqr + dqs + dqg
                    adj_vmr (is:ie, kr) = (ps_dt - (qv (is:ie, k) + &
                        qa (is:ie, k, liq_wat) + qa (is:ie, k, ice_wat) + &
                        qa (is:ie, k, rainwat) + qa (is:ie, k, snowwat) + &
                        qa (is:ie, k, graupel))) / (1. - (qv (is:ie, k) + &
                        qa (is:ie, k, liq_wat) + qa (is:ie, k, ice_wat) + &
                        qa (is:ie, k, rainwat) + qa (is:ie, k, snowwat) + &
                        qa (is:ie, k, graupel))) / ps_dt
                    q (is:ie, j, kr, sphum) = qv (is:ie, k) / ps_dt
                    q (is:ie, j, kr, liq_wat) = qa (is:ie, k, liq_wat) / ps_dt
                    q (is:ie, j, kr, ice_wat) = qa (is:ie, k, ice_wat) / ps_dt
                    q (is:ie, j, kr, rainwat) = qa (is:ie, k, rainwat) / ps_dt
                    q (is:ie, j, kr, snowwat) = qa (is:ie, k, snowwat) / ps_dt
                    q (is:ie, j, kr, graupel) = qa (is:ie, k, graupel) / ps_dt
                endif
                delp (is:ie, j, kr) = delp (is:ie, j, kr) * ps_dt
                q_liq = q (is:ie, j, kr, liq_wat) + q (is:ie, j, kr, rainwat)
                q_sol = q (is:ie, j, kr, ice_wat) + q (is:ie, j, kr, snowwat) + q (is:ie, j, kr, graupel)
                if (thermostruct%use_cond) then
                    q_con (is:ie, j, kr) = q_liq + q_sol
                endif
                c_moist = (1 - (q (is:ie, j, kr, sphum) + q_liq + q_sol)) * cv_air + &
                    q (is:ie, j, kr, sphum) * cv_vap + q_liq * c_liq + q_sol * c_ice
                if (thermostruct%moist_kappa) then
                    cappa (is:ie, j, kr) = rdgas / (rdgas + c_moist / (1. + r_vir * q (is:ie, j, kr, sphum)))
                endif
                pt (is:ie, j, kr) = pt (is:ie, j, kr) + (ta (is:ie, k) * &
                    ((1. + r_vir * q (is:ie, j, kr, sphum)) * (1. - (q_liq + q_sol))) - &
                    pt (is:ie, j, kr)) * cp_air / c_moist
                dte8 (is:ie, k) = te8 (is:ie, k) - (c_moist * pt (is:ie, j, kr) / &
                    ((1. + r_vir * q (is:ie, j, kr, sphum)) * (1. - (q_liq + q_sol))) + &
                    (hlv - rvgas * tice - (cv_vap - c_liq) * tice) * q (is:ie, j, kr, sphum) - &
                    (hlf - (c_liq - c_ice) * tice) * q_sol) * delp (is:ie, j, kr) / grav
                ua (is:ie, j, kr) = uu (is:ie, k)
                va (is:ie, j, kr) = vv (is:ie, k)
                k2 = 0.5 * (ua (is:ie, j, kr) ** 2 + va (is:ie, j, kr) ** 2 + w (is:ie, j, kr) ** 2) * delp (is:ie, j, kr)
                pt (is:ie, j, kr) = pt (is:ie, j, kr) + (k1 - k2) / c_moist / delp (is:ie, j, kr) * &
                    ((1. + r_vir * q (is:ie, j, kr, sphum)) * (1. - (q_liq + q_sol)))
            enddo

            ! update non-microphyiscs tracers due to mass change
            if (adj_mass_vmr .gt. 0) then
                do m = 1, nq
                    if (conv_vmr_mmr (m)) then
                        q (is:ie, j, 1:km, m) = q (is:ie, j, 1:km, m) * adj_vmr (is:ie, 1:km)
                    endif
                enddo
            endif

            ! compute wind tendency at A grid fori D grid wind update
            u_dt (is:ie, j, 1:km) = (ua (is:ie, j, 1:km) - u_dt (is:ie, j, 1:km)) / abs (mdt)
            v_dt (is:ie, j, 1:km) = (va (is:ie, j, 1:km) - v_dt (is:ie, j, 1:km)) / abs (mdt)

            ! update pkz
            if (.not. hydrostatic) then
                if (thermostruct%moist_kappa) then
                    pkz (is:ie, j, 1:km) = exp (cappa (is:ie, j, 1:km) * &
                        log (rrg * delp (is:ie, j, 1:km) / &
                        delz (is:ie, j, 1:km) * pt (is:ie, j, 1:km)))
                else
                    pkz (is:ie, j, 1:km) = exp (akap * log (rrg * delp (is:ie, j, 1:km) / &
                        delz (is:ie, j, 1:km) * pt (is:ie, j, 1:km)))
                endif
            endif

            ! total energy checker
            if (consv_checker) then
                qliq = q (is:ie, j, 1:km, liq_wat) + q (is:ie, j, 1:km, rainwat)
                qsol = q (is:ie, j, 1:km, ice_wat) + q (is:ie, j, 1:km, snowwat) + q (is:ie, j, 1:km, graupel)
                te_end (is:ie, 1:km) = 0.0
                tw_end (is:ie, 1:km) = 0.0
                te_b_end (is:ie) = 0.0
                tw_b_end (is:ie) = 0.0
                do i = is, ie
                    tz = pt (i, j, 1:km) / ((1. + r_vir * q (i, j, 1:km, sphum)) * (1. - (qliq (i, 1:km) + qsol (i, 1:km))))
                    if (hydrostatic) then
                        wz = 0.0
                    else
                        wz = w (i, j, 1:km)
                    endif
                    dte (i) = sum (dte8 (i, 1:km))
                    call mtetw (1, km, q (i, j, 1:km, sphum), q (i, j, 1:km, liq_wat), &
                        q (i, j, 1:km, rainwat), q (i, j, 1:km, ice_wat), q (i, j, 1:km, snowwat), &
                        q (i, j, 1:km, graupel), tz, ua (i, j, 1:km), va (i, j, 1:km), wz, &
                        delp (i, j, 1:km), dte (i), 0.0, 0.0, inline_cnv%prec (i, j) / abs (mdt) * 1.e3 * 86400, 0.0, 0.0, &
                        0.0, 0.0, 0.0, abs (mdt), te_end (i, 1:km), tw_end (i, 1:km), &
                        te_b_end (i), tw_b_end (i), .true., hydrostatic, te_loss (i))
                enddo
            endif

            ! total energy after parameterization, add total energy change to te0_2d
            if (consv .gt. consv_min) then
                qliq = q (is:ie, j, 1:km, liq_wat) + q (is:ie, j, 1:km, rainwat)
                qsol = q (is:ie, j, 1:km, ice_wat) + q (is:ie, j, 1:km, snowwat) + q (is:ie, j, 1:km, graupel)
                cvm = (1 - (q (is:ie, j, 1:km, sphum) + qliq + qsol)) * cv_air + &
                    q (is:ie, j, 1:km, sphum) * cv_vap + qliq * c_liq + qsol * c_ice
                te (is:ie, j, 1:km) = te (is:ie, j, 1:km) + &
                    cvm * pt (is:ie, j, 1:km) / ((1. + r_vir * q (is:ie, j, 1:km, sphum)) * &
                    (1. - (qliq + qsol))) * delp (is:ie, j, 1:km)
                do k = 1, km
                    te0_2d (is:ie, j) = te0_2d (is:ie, j) + te (is:ie, j, k)
                enddo
            endif

            ! total energy checker
            if (consv_checker) then
                do i = is, ie
                    !if (abs (sum (te_end (i, :)) + te_b_end (i) - sum (te_beg (i, :)) - te_b_beg (i)) / &
                    !     (sum (te_beg (i, :)) + te_b_beg (i)) .gt. te_err) then
                    !    print*, "CNV-INTM TE: ", &
                    !        !(sum (te_beg (i, :)) + te_b_beg (i)), &
                    !        !(sum (te_end (i, :)) + te_b_end (i)), &
                    !        (sum (te_end (i, :)) + te_b_end (i) - sum (te_beg (i, :)) - te_b_beg (i)) / &
                    !        (sum (te_beg (i, :)) + te_b_beg (i))
                    !endif
                    inline_cnv%intm_te_a_chg (i, j) = sum (te_end (i, :)) - sum (te_beg (i, :))
                    inline_cnv%intm_te_b_chg (i, j) = te_b_end (i) - te_b_beg (i)
                    !if (abs (sum (tw_end (i, :)) + tw_b_end (i) - sum (tw_beg (i, :)) - tw_b_beg (i)) / &
                    !     (sum (tw_beg (i, :)) + tw_b_beg (i)) .gt. tw_err) then
                    !    print*, "CNV-INTM TW: ", &
                    !        !(sum (tw_beg (i, :)) + tw_b_beg (i)), &
                    !        !(sum (tw_end (i, :)) + tw_b_end (i)), &
                    !        (sum (tw_end (i, :)) + tw_b_end (i) - sum (tw_beg (i, :)) - tw_b_beg (i)) / &
                    !        (sum (tw_beg (i, :)) + tw_b_beg (i))
                    !endif
                    inline_cnv%intm_tw_a_chg (i, j) = sum (tw_end (i, :)) - sum (tw_beg (i, :))
                    inline_cnv%intm_tw_b_chg (i, j) = tw_b_end (i) - tw_b_beg (i)
                    !print*, "CNV-INTM LOSS (%) : ", te_loss (i) / (sum (te_beg (i, :)) + te_b_beg (i)) * 100.0
                enddo
            endif

        enddo

        if (inline_cnv_flag .eq. 2) then
            deallocate (fscav)
        endif

        deallocate (rn)
        deallocate (tmp)

        deallocate (dz)
        deallocate (zm)
        deallocate (dp)
        deallocate (pm)
        deallocate (pi)

        deallocate (ta)
        deallocate (qv)
        deallocate (qr)
        deallocate (uu)
        deallocate (vv)
        deallocate (ww)

        if (inline_cnv_flag .eq. 1) deallocate (ql)
        if (inline_cnv_flag .eq. 2) deallocate (qa)

        deallocate (tz)
        deallocate (wz)

        ! Note: (ua, va) are *lat-lon* wind tendenies on cell centers
        call timing_on('COMM_TOTAL')
        if ( gridstruct%square_domain ) then
            call mpp_update_domains (u_dt, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.false.)
            call mpp_update_domains (v_dt, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
        else
            call mpp_update_domains (u_dt, domain, complete=.false.)
            call mpp_update_domains (v_dt, domain, complete=.true.)
        endif
        call timing_off('COMM_TOTAL')

        ! update D grid wind
        call update_dwinds_phys (is, ie, js, je, isd, ied, jsd, jed, abs (mdt), u_dt, v_dt, u, v, &
                 gridstruct, npx, npy, km, domain)

        deallocate (u_dt)
        deallocate (v_dt)

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
    ! <<< Inline Convection
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! Inline Gravity Wave Drag >>>
    !-----------------------------------------------------------------------

    if ((.not. do_adiabatic_init) .and. do_inline_gwd) then

        allocate (dz (is:ie, 1:km))
        allocate (zm (is:ie, 1:km))
        allocate (zi (is:ie, 1:km+1))
        allocate (dp (is:ie, 1:km))
        allocate (pm (is:ie, 1:km))
        allocate (pi (is:ie, 1:km+1))
        allocate (pmk (is:ie, 1:km))

        allocate (ta (is:ie, 1:km))
        allocate (qv (is:ie, 1:km))
        allocate (uu (is:ie, 1:km))
        allocate (vv (is:ie, 1:km))

        allocate (u_dt (isd:ied, jsd:jed, km))
        allocate (v_dt (isd:ied, jsd:jed, km))

        allocate (tz (1:km))
        allocate (wz (1:km))

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
                 domain, gridstruct%bounded_domain, 4, bd)

        ! save delp
        if (consv .gt. consv_min) then
            allocate (dp0 (isd:ied, jsd:jed, km))
            dp0 = delp
        endif

!$OMP parallel do default (none) shared (is, ie, js, je, isd, jsd, km, ua, va, w, &
!$OMP                                    te, delp, hydrostatic, pt, delz, q_con, &
!$OMP                                    rainwat, liq_wat, ice_wat, snowwat, graupel, &
!$OMP                                    sphum, pkz, consv, te0_2d, gridstruct, q, &
!$OMP                                    mdt, cappa, rrg, akap, r_vir, inline_gwd, &
!$OMP                                    ptop, inline_pbl, inline_cnv, u_dt, v_dt, &
!$OMP                                    safety_check, do_fast_phys, &
!$OMP                                    conv_vmr_mmr, nq, consv_checker, &
!$OMP                                    te_err, tw_err, thermostruct) &
!$OMP                           private (gsize, dz, pi, pmk, zi, q_liq, q_sol, pe, &
!$OMP                                    zm, dp, pm, qv, ta, uu, vv, qliq, qsol, &
!$OMP                                    cvm, kr, c_moist, peln, &
!$OMP                                    tz, wz, dte, te_beg, tw_beg, te_b_beg, tw_b_beg, &
!$OMP                                    te_end, tw_end, te_b_end, tw_b_end, te_loss)

        do j = js, je

            ! grid size
            gsize (is:ie) = sqrt (gridstruct%area_64 (is:ie, j))

            ! save ua, va for wind tendency calculation
            u_dt (is:ie, j, 1:km) = ua (is:ie, j, 1:km)
            v_dt (is:ie, j, 1:km) = va (is:ie, j, 1:km)

            ! total energy before parameterization
            if (consv .gt. consv_min) then
                qliq = q (is:ie, j, 1:km, liq_wat) + q (is:ie, j, 1:km, rainwat)
                qsol = q (is:ie, j, 1:km, ice_wat) + q (is:ie, j, 1:km, snowwat) + q (is:ie, j, 1:km, graupel)
                cvm = (1 - (q (is:ie, j, 1:km, sphum) + qliq + qsol)) * cv_air + &
                    q (is:ie, j, 1:km, sphum) * cv_vap + qliq * c_liq + qsol * c_ice
                te (is:ie, j, 1:km) = - cvm * pt (is:ie, j, 1:km) / ((1. + r_vir * q (is:ie, j, 1:km, sphum)) * &
                    (1. - (qliq + qsol))) * delp (is:ie, j, 1:km)
            endif

            ! total energy checker
            if (consv_checker) then
                qliq = q (is:ie, j, 1:km, liq_wat) + q (is:ie, j, 1:km, rainwat)
                qsol = q (is:ie, j, 1:km, ice_wat) + q (is:ie, j, 1:km, snowwat) + q (is:ie, j, 1:km, graupel)
                te_beg (is:ie, 1:km) = 0.0
                tw_beg (is:ie, 1:km) = 0.0
                te_b_beg (is:ie) = 0.0
                tw_b_beg (is:ie) = 0.0
                do i = is, ie
                    tz = pt (i, j, 1:km) / ((1. + r_vir * q (i, j, 1:km, sphum)) * (1. - (qliq (i, 1:km) + qsol (i, 1:km))))
                    if (hydrostatic) then
                        wz = 0.0
                    else
                        wz = w (i, j, 1:km)
                    endif
                    dte (i) = 0.0
                    call mtetw (1, km, q (i, j, 1:km, sphum), q (i, j, 1:km, liq_wat), &
                        q (i, j, 1:km, rainwat), q (i, j, 1:km, ice_wat), q (i, j, 1:km, snowwat), &
                        q (i, j, 1:km, graupel), tz, ua (i, j, 1:km), va (i, j, 1:km), wz, &
                        delp (i, j, 1:km), dte (i), 0.0, 0.0, 0.0, 0.0, 0.0, &
                        0.0, 0.0, 0.0, abs (mdt), te_beg (i, 1:km), tw_beg (i, 1:km), &
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

            ! vertical index flip over
            zi (is:ie, 1) = 0.0
            pi (is:ie, 1) = pe (is:ie, km+1)
            do k = 1, km
                kr = km - k + 1
                dp (is:ie, k) = delp (is:ie, j, kr)
                pi (is:ie, k+1) = pe (is:ie, kr)
                if (.not. hydrostatic) then
                    pm (is:ie, k) = dp (is:ie, k) / delz (is:ie, j, kr) * &
                        rrg * pt (is:ie, j, kr)
                    dz (is:ie, k) = delz (is:ie, j, kr)
                    ! ensure subgrid monotonicity of pressure
                    do i = is, ie
                        pm (i, k) = min (pm (i, k), pi (i, k) - 0.01 * pm (i, k))
                        pm (i, k) = max (pm (i, k), pi (i, k+1) + 0.01 * pm (i, k))
                    enddo
                else
                    pm (is:ie, k) = dp (is:ie, k) / (peln (is:ie, kr+1) - peln (is:ie, kr))
                    dz (is:ie, k) = (peln (is:ie, kr+1) - peln (is:ie, kr)) * &
                        rrg * pt (is:ie, j, kr)
                endif
                pmk (is:ie, k) = exp (kappa * log (pm (is:ie, k) * 1.e-5))
                zi (is:ie, k+1) = zi (is:ie, k) - dz (is:ie, k) * grav
                if (k .eq. 1) then
                    zm (is:ie, k) = - 0.5 * dz (is:ie, k) * grav
                else
                    zm (is:ie, k) = zm (is:ie, k-1) - 0.5 * (dz (is:ie, k-1) + dz (is:ie, k)) * grav
                endif
                qv (is:ie, k) = q (is:ie, j, kr, sphum)
                q_liq = q (is:ie, j, kr, liq_wat) + q (is:ie, j, kr, rainwat)
                q_sol = q (is:ie, j, kr, ice_wat) + q (is:ie, j, kr, snowwat) + q (is:ie, j, kr, graupel)
                ta (is:ie, k) = pt (is:ie, j, kr) / ((1. + r_vir * q (is:ie, j, kr, sphum)) * &
                    (1. - (q_liq + q_sol)))
                uu (is:ie, k) = ua (is:ie, j, kr)
                vv (is:ie, k) = va (is:ie, j, kr)
            enddo

            ! check if pressure or height cross over
            if (safety_check) then
                do k = 1, km
                    do i = is, ie
                        if (k .lt. km) then
                            if (pm (i, k) .le. pm (i, k+1)) then
                                print*, "Warning: inline gwd pressure layer cross over", k, pm (i, k), pm (i, k+1)
                            endif
                            if (zm (i, k) .ge. zm (i, k+1)) then
                                print*, "Warning: inline gwd height layer cross over", k, zm (i, k), zm (i, k+1)
                            endif
                        endif
                        if (pi (i, k) .le. pi (i, k+1)) then
                            print*, "Warning: inline gwd pressure interface cross over", k, pi (i, k), pi (i, k+1)
                        endif
                        if (zi (i, k) .ge. zi (i, k+1)) then
                            print*, "Warning: inline gwd height interface cross over", k, zi (i, k), zi (i, k+1)
                        endif
                    enddo
                enddo
            endif

            if (.not. do_fast_phys) then

                ! orographic gravity wave drag and mountain blocking main program
                call sa_gwd_oro (ie-is+1, km, uu, vv, ta, qv, abs (mdt), gsize, &
                    inline_pbl%kpbl (is:ie, j), pi, dp, pm, pmk, zi, zm, &
                    inline_gwd%hprime (is:ie, j), inline_gwd%oc (is:ie, j), inline_gwd%oa (is:ie, j, :), &
                    inline_gwd%ol (is:ie, j, :), inline_gwd%theta (is:ie, j), inline_gwd%sigma (is:ie, j), &
                    inline_gwd%gamma (is:ie, j), inline_gwd%elvmax (is:ie, j))

            endif

            ! convective gravity wave drag main program
            call sa_gwd_cnv (ie-is+1, km, uu, vv, ta, qv, abs (mdt), gsize, inline_cnv%cumabs (is:ie, j), &
                pm, pi, dp, inline_cnv%ktop (is:ie, j), inline_cnv%kbot (is:ie, j), inline_cnv%kcnv (is:ie, j))

            ! update u, v, T, q, and delp, vertical index flip over
            do k = 1, km
                kr = km - k + 1
                q_liq = q (is:ie, j, kr, liq_wat) + q (is:ie, j, kr, rainwat)
                q_sol = q (is:ie, j, kr, ice_wat) + q (is:ie, j, kr, snowwat) + q (is:ie, j, kr, graupel)
                if (thermostruct%use_cond) then
                    q_con (is:ie, j, kr) = q_liq + q_sol
                endif
                c_moist = (1 - (q (is:ie, j, kr, sphum) + q_liq + q_sol)) * cv_air + &
                    q (is:ie, j, kr, sphum) * cv_vap + q_liq * c_liq + q_sol * c_ice
                if (thermostruct%moist_kappa) then
                    cappa (is:ie, j, kr) = rdgas / (rdgas + c_moist / (1. + r_vir * q (is:ie, j, kr, sphum)))
                endif
                pt (is:ie, j, kr) = pt (is:ie, j, kr) + (ta (is:ie, k) * &
                    ((1. + r_vir * q (is:ie, j, kr, sphum)) * (1. - (q_liq + q_sol))) - &
                    pt (is:ie, j, kr)) * cp_air / c_moist
                ua (is:ie, j, kr) = uu (is:ie, k)
                va (is:ie, j, kr) = vv (is:ie, k)
            enddo

            ! compute wind tendency at A grid fori D grid wind update
            u_dt (is:ie, j, 1:km) = (ua (is:ie, j, 1:km) - u_dt (is:ie, j, 1:km)) / abs (mdt)
            v_dt (is:ie, j, 1:km) = (va (is:ie, j, 1:km) - v_dt (is:ie, j, 1:km)) / abs (mdt)

            ! update pkz
            if (.not. hydrostatic) then
                if (thermostruct%moist_kappa) then
                    pkz (is:ie, j, 1:km) = exp (cappa (is:ie, j, 1:km) * &
                        log (rrg * delp (is:ie, j, 1:km) / &
                        delz (is:ie, j, 1:km) * pt (is:ie, j, 1:km)))
                else
                    pkz (is:ie, j, 1:km) = exp (akap * log (rrg * delp (is:ie, j, 1:km) / &
                        delz (is:ie, j, 1:km) * pt (is:ie, j, 1:km)))
                endif
            endif

            ! total energy checker
            if (consv_checker) then
                qliq = q (is:ie, j, 1:km, liq_wat) + q (is:ie, j, 1:km, rainwat)
                qsol = q (is:ie, j, 1:km, ice_wat) + q (is:ie, j, 1:km, snowwat) + q (is:ie, j, 1:km, graupel)
                te_end (is:ie, 1:km) = 0.0
                tw_end (is:ie, 1:km) = 0.0
                te_b_end (is:ie) = 0.0
                tw_b_end (is:ie) = 0.0
                do i = is, ie
                    tz = pt (i, j, 1:km) / ((1. + r_vir * q (i, j, 1:km, sphum)) * (1. - (qliq (i, 1:km) + qsol (i, 1:km))))
                    if (hydrostatic) then
                        wz = 0.0
                    else
                        wz = w (i, j, 1:km)
                    endif
                    dte (i) = 0.0
                    call mtetw (1, km, q (i, j, 1:km, sphum), q (i, j, 1:km, liq_wat), &
                        q (i, j, 1:km, rainwat), q (i, j, 1:km, ice_wat), q (i, j, 1:km, snowwat), &
                        q (i, j, 1:km, graupel), tz, ua (i, j, 1:km), va (i, j, 1:km), wz, &
                        delp (i, j, 1:km), dte (i), 0.0, 0.0, 0.0, 0.0, 0.0, &
                        0.0, 0.0, 0.0, abs (mdt), te_end (i, 1:km), tw_end (i, 1:km), &
                        te_b_end (i), tw_b_end (i), .true., hydrostatic, te_loss (i))
                enddo
            endif

            ! total energy after parameterization, add total energy change to te0_2d
            if (consv .gt. consv_min) then
                qliq = q (is:ie, j, 1:km, liq_wat) + q (is:ie, j, 1:km, rainwat)
                qsol = q (is:ie, j, 1:km, ice_wat) + q (is:ie, j, 1:km, snowwat) + q (is:ie, j, 1:km, graupel)
                cvm = (1 - (q (is:ie, j, 1:km, sphum) + qliq + qsol)) * cv_air + &
                    q (is:ie, j, 1:km, sphum) * cv_vap + qliq * c_liq + qsol * c_ice
                te (is:ie, j, 1:km) = te (is:ie, j, 1:km) + &
                    cvm * pt (is:ie, j, 1:km) / ((1. + r_vir * q (is:ie, j, 1:km, sphum)) * &
                    (1. - (qliq + qsol))) * delp (is:ie, j, 1:km)
                do k = 1, km
                    te0_2d (is:ie, j) = te0_2d (is:ie, j) + te (is:ie, j, k)
                enddo
            endif

            ! total energy checker
            if (consv_checker) then
                do i = is, ie
                    !if (abs (sum (te_end (i, :)) + te_b_end (i) - sum (te_beg (i, :)) - te_b_beg (i)) / &
                    !     (sum (te_beg (i, :)) + te_b_beg (i)) .gt. te_err) then
                    !    print*, "GWD-INTM TE: ", &
                    !        !(sum (te_beg (i, :)) + te_b_beg (i)), &
                    !        !(sum (te_end (i, :)) + te_b_end (i)), &
                    !        (sum (te_end (i, :)) + te_b_end (i) - sum (te_beg (i, :)) - te_b_beg (i)) / &
                    !        (sum (te_beg (i, :)) + te_b_beg (i))
                    !endif
                    inline_gwd%intm_te_a_chg (i, j) = sum (te_end (i, :)) - sum (te_beg (i, :))
                    inline_gwd%intm_te_b_chg (i, j) = te_b_end (i) - te_b_beg (i)
                    !if (abs (sum (tw_end (i, :)) + tw_b_end (i) - sum (tw_beg (i, :)) - tw_b_beg (i)) / &
                    !     (sum (tw_beg (i, :)) + tw_b_beg (i)) .gt. tw_err) then
                    !    print*, "GWD-INTM TW: ", &
                    !        !(sum (tw_beg (i, :)) + tw_b_beg (i)), &
                    !        !(sum (tw_end (i, :)) + tw_b_end (i)), &
                    !        (sum (tw_end (i, :)) + tw_b_end (i) - sum (tw_beg (i, :)) - tw_b_beg (i)) / &
                    !        (sum (tw_beg (i, :)) + tw_b_beg (i))
                    !endif
                    inline_gwd%intm_tw_a_chg (i, j) = sum (tw_end (i, :)) - sum (tw_beg (i, :))
                    inline_gwd%intm_tw_b_chg (i, j) = tw_b_end (i) - tw_b_beg (i)
                    !print*, "GWD-INTM LOSS (%) : ", te_loss (i) / (sum (te_beg (i, :)) + te_b_beg (i)) * 100.0
                enddo
            endif

        enddo

        deallocate (dz)
        deallocate (zm)
        deallocate (dp)
        deallocate (pm)
        deallocate (pi)
        deallocate (pmk)

        deallocate (ta)
        deallocate (qv)
        deallocate (uu)
        deallocate (vv)

        deallocate (tz)
        deallocate (wz)

        ! Note: (ua, va) are *lat-lon* wind tendenies on cell centers
        call timing_on('COMM_TOTAL')
        if ( gridstruct%square_domain ) then
            call mpp_update_domains (u_dt, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.false.)
            call mpp_update_domains (v_dt, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
        else
            call mpp_update_domains (u_dt, domain, complete=.false.)
            call mpp_update_domains (v_dt, domain, complete=.true.)
        endif
        call timing_off('COMM_TOTAL')

        ! update D grid wind
        call update_dwinds_phys (is, ie, js, je, isd, ied, jsd, jed, abs (mdt), u_dt, v_dt, u, v, &
                 gridstruct, npx, npy, km, domain)

        deallocate (u_dt)
        deallocate (v_dt)

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
    ! <<< Inline Gravity Wave Drag
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! Inline Microphysics >>>
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
                 domain, gridstruct%bounded_domain, 4, bd)

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
!$OMP                                    te_err, tw_err, k_con, k_cappa, thermostruct) &
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
                     q_con (is:ie, j, k_con:), cappa (is:ie, j, k_cappa:), &
                     consv .gt. consv_min, adj_vmr (is:ie, kmp:km), te (is:ie, j, kmp:km), dte (is:ie), &
                     inline_mp%prefluxw(is:ie, j, kmp:km), &
                     inline_mp%prefluxr(is:ie, j, kmp:km), inline_mp%prefluxi(is:ie, j, kmp:km), &
                     inline_mp%prefluxs(is:ie, j, kmp:km), inline_mp%prefluxg(is:ie, j, kmp:km), &
                     inline_mp%mppcw (is:ie, j), inline_mp%mppew (is:ie, j), inline_mp%mppe1 (is:ie, j), &
                     inline_mp%mpper (is:ie, j), inline_mp%mppdi (is:ie, j), inline_mp%mppd1 (is:ie, j), &
                     inline_mp%mppds (is:ie, j), inline_mp%mppdg (is:ie, j), inline_mp%mppsi (is:ie, j), &
                     inline_mp%mpps1 (is:ie, j), inline_mp%mppss (is:ie, j), inline_mp%mppsg (is:ie, j), &
                     inline_mp%mppfw (is:ie, j), inline_mp%mppfr (is:ie, j), inline_mp%mppmi (is:ie, j), &
                     inline_mp%mppms (is:ie, j), inline_mp%mppmg (is:ie, j), inline_mp%mppm1 (is:ie, j), &
                     inline_mp%mppm2 (is:ie, j), inline_mp%mppm3 (is:ie, j), inline_mp%mppar (is:ie, j), &
                     inline_mp%mppas (is:ie, j), inline_mp%mppag (is:ie, j), inline_mp%mpprs (is:ie, j), &
                     inline_mp%mpprg (is:ie, j), inline_mp%mppxr (is:ie, j), inline_mp%mppxs (is:ie, j), &
                     inline_mp%mppxg (is:ie, j), last_step, do_inline_mp, &
                     thermostruct%use_cond, thermostruct%moist_kappa)

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
               if (thermostruct%moist_kappa) then
                  pkz (is:ie, j, kmp:km) = exp (cappa (is:ie, j, kmp:km) * &
                       log (rrg * delp (is:ie, j, kmp:km) / &
                       delz (is:ie, j, kmp:km) * pt (is:ie, j, kmp:km)))
               else
                  pkz (is:ie, j, kmp:km) = exp (akap * log (rrg * delp (is:ie, j, kmp:km) / &
                       delz (is:ie, j, kmp:km) * pt (is:ie, j, kmp:km)))
               endif
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
                    !if (abs (sum (te_end (i, kmp:km)) + te_b_end (i) - sum (te_beg (i, kmp:km)) - te_b_beg (i)) / &
                    !     (sum (te_beg (i, kmp:km)) + te_b_beg (i)) .gt. te_err) then
                    !    print*, "MP-INTM TE: ", &
                    !        !(sum (te_beg (i, kmp:km)) + te_b_beg (i)), &
                    !        !(sum (te_end (i, kmp:km)) + te_b_end (i)), &
                    !        (sum (te_end (i, kmp:km)) + te_b_end (i) - sum (te_beg (i, kmp:km)) - te_b_beg (i)) / &
                    !        (sum (te_beg (i, kmp:km)) + te_b_beg (i))
                    !endif
                    inline_mp%intm_te_a_chg (i, j) = sum (te_end (i, :)) - sum (te_beg (i, :))
                    inline_mp%intm_te_b_chg (i, j) = te_b_end (i) - te_b_beg (i)
                    !if (abs (sum (tw_end (i, kmp:km)) + tw_b_end (i) - sum (tw_beg (i, kmp:km)) - tw_b_beg (i)) / &
                    !     (sum (tw_beg (i, kmp:km)) + tw_b_beg (i)) .gt. tw_err) then
                    !    print*, "MP-INTM TW: ", &
                    !        !(sum (tw_beg (i, kmp:km)) + tw_b_beg (i)), &
                    !        !(sum (tw_end (i, kmp:km)) + tw_b_end (i)), &
                    !        (sum (tw_end (i, kmp:km)) + tw_b_end (i) - sum (tw_beg (i, kmp:km)) - tw_b_beg (i)) / &
                    !        (sum (tw_beg (i, kmp:km)) + tw_b_beg (i))
                    !endif
                    inline_mp%intm_tw_a_chg (i, j) = sum (tw_end (i, :)) - sum (tw_beg (i, :))
                    inline_mp%intm_tw_b_chg (i, j) = tw_b_end (i) - tw_b_beg (i)
                    !print*, "MP-INTM LOSS (%) : ", te_loss (i) / (sum (te_beg (i, kmp:km)) + te_b_beg (i)) * 100.0
                enddo
            endif

        enddo

        deallocate (dz)
        deallocate (wa)

        ! Note: (ua, va) are *lat-lon* wind tendenies on cell centers
        call timing_on('COMM_TOTAL')
        if ( gridstruct%square_domain ) then
            call mpp_update_domains (u_dt, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.false.)
            call mpp_update_domains (v_dt, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
        else
            call mpp_update_domains (u_dt, domain, complete=.false.)
            call mpp_update_domains (v_dt, domain, complete=.true.)
        endif
        call timing_off('COMM_TOTAL')

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
    ! <<< Inline Microphysics
    !-----------------------------------------------------------------------

end subroutine intermediate_phys

end module intermediate_phys_mod
