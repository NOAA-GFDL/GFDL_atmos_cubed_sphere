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
! Fast Physics Interface
! Developer: Linjiong Zhou
! Last Update: 5/19/2022
! =======================================================================

module fast_phys_mod

#ifdef OVERLOAD_R4
    use constantsR4_mod, only: rdgas, rvgas, grav, kappa, cp_air
#else
    use constants_mod, only: rdgas, rvgas, grav, kappa, cp_air
#endif
    use fv_grid_utils_mod, only: cubed_to_latlon, update_dwinds_phys
    use fv_arrays_mod, only: fv_grid_type, fv_grid_bounds_type, fv_thermo_type
    use fv_arrays_mod, only: inline_pbl_type, inline_gwd_type
    use mpp_domains_mod, only: domain2d, mpp_update_domains
    use tracer_manager_mod, only: get_tracer_index, get_tracer_names
    use field_manager_mod, only: model_atmos
    use gfdl_mp_mod, only: c_liq, c_ice, cv_air, cv_vap, hlv, mtetw, tice
    use sa_tke_edmf_mod, only: sa_tke_edmf_sfc, sa_tke_edmf_pbl
    use sa_tke_edmf_new_mod, only: sa_tke_edmf_new_sfc, sa_tke_edmf_new_pbl
    use sa_gwd_mod, only: sa_gwd_oro
    use fv_timing_mod, only: timing_on, timing_off

    implicit none

    private

    real, parameter :: consv_min = 0.001

    public :: fast_phys

    ! -----------------------------------------------------------------------
    ! precision definition
    ! -----------------------------------------------------------------------
    
    integer, parameter :: r8 = 8 ! double precision

contains

subroutine fast_phys (is, ie, js, je, isd, ied, jsd, jed, km, npx, npy, nq, nwat, &
               mdt, consv, akap, ptop, hs, te0_2d, u, v, w, pt, &
               delp, delz, q_con, cappa, q, pkz, r_vir, te_err, tw_err, inline_pbl, inline_gwd, &
               gridstruct, thermostruct, domain, bd, hydrostatic, do_adiabatic_init, &
               do_inline_pbl, do_inline_gwd, consv_checker, adj_mass_vmr, inline_pbl_flag)

    implicit none

    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------

    integer, intent (in) :: is, ie, js, je, isd, ied, jsd, jed, km, npx, npy, nq, nwat, adj_mass_vmr, inline_pbl_flag

    logical, intent (in) :: hydrostatic, do_adiabatic_init, do_inline_pbl, do_inline_gwd
    logical, intent (in) :: consv_checker

    real, intent (in) :: consv, mdt, akap, r_vir, ptop, te_err, tw_err

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

    type (fv_thermo_type), intent (in), target :: thermostruct

    type (fv_grid_bounds_type), intent (in) :: bd

    type (domain2d), intent (inout) :: domain

    type (inline_pbl_type), intent (inout) :: inline_pbl

    type (inline_gwd_type), intent (inout) :: inline_gwd

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    logical :: safety_check = .true.

    logical, allocatable, dimension (:) :: conv_vmr_mmr

    integer :: i, j, k, m, kr, kmp, ncld, ntke, sphum, liq_wat, ice_wat, lsoil
    integer :: rainwat, snowwat, graupel, cld_amt, ccn_cm3, cin_cm3, aerosol

    real :: rrg

    real, dimension (is:ie) :: gsize, dqv, dql, dqi, dqr, dqs, dqg, ps_dt, q_liq, q_sol, c_moist

    real, dimension (is:ie, km) :: q2, q3, qliq, qsol, cvm, adj_vmr

    real, dimension (is:ie, km+1) :: phis, pe, peln

    real, dimension (isd:ied, jsd:jed, km) :: te, ua, va

    integer, allocatable, dimension (:) :: kinver, vegtype

    real, allocatable, dimension (:) :: rn, rb, u10m, v10m, sigmaf, stress, wind, tmp, wz

    real, allocatable, dimension (:) :: dtsfc, dqvsfc, dqlsfc, dqisfc, dqrsfc, dqssfc, dqgsfc

    real, allocatable, dimension (:,:) :: dz, zm, zi, wa, dp, pm, pi, pmk, pik, qv, ql, ta, uu, vv, ww, radh

    real, allocatable, dimension (:,:,:) :: u_dt, v_dt, dp0, u0, v0, qa

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
    ntke = get_tracer_index (model_atmos, 'sgs_tke')

    rrg = - rdgas / grav

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
    ! pt conversion
    !-----------------------------------------------------------------------

    if (thermostruct%moist_kappa) then
       do k = 1, km
           do j = js, je
               do i = is, ie
                   pt (i, j, k) = pt (i, j, k) * exp (cappa (i, j, k) / (1. - cappa (i, j, k)) * &
                       log (rrg * delp (i, j, k) / delz (i, j, k) * pt (i, j, k)))
               enddo
           enddo
       enddo
    else
       do k = 1, km
           do j = js, je
               do i = is, ie
                   pt (i, j, k) = pt (i, j, k) * exp (akap / (1 - akap) * &
                       log (rrg * delp (i, j, k) / delz (i, j, k) * pt (i, j, k)))
               enddo
           enddo
       enddo
    endif

    !-----------------------------------------------------------------------
    ! Inline Planetary Boundary Layer >>>
    !-----------------------------------------------------------------------

    if ((.not. do_adiabatic_init) .and. do_inline_pbl) then

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

            if (inline_pbl_flag .eq. 1) then

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

            endif

            if (inline_pbl_flag .eq. 2) then

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

            endif

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
                    !    print*, "PBL-FAST TE: ", &
                    !        !(sum (te_beg (i, :)) + te_b_beg (i)), &
                    !        !(sum (te_end (i, :)) + te_b_end (i)), &
                    !        (sum (te_end (i, :)) + te_b_end (i) - sum (te_beg (i, :)) - te_b_beg (i)) / &
                    !        (sum (te_beg (i, :)) + te_b_beg (i))
                    !endif
                    inline_pbl%fast_te_a_chg (i, j) = sum (te_end (i, :)) - sum (te_beg (i, :))
                    inline_pbl%fast_te_b_chg (i, j) = te_b_end (i) - te_b_beg (i)
                    !if (abs (sum (tw_end (i, :)) + tw_b_end (i) - sum (tw_beg (i, :)) - tw_b_beg (i)) / &
                    !     (sum (tw_beg (i, :)) + tw_b_beg (i)) .gt. tw_err) then
                    !    print*, "PBL-FAST TW: ", &
                    !        !(sum (tw_beg (i, :)) + tw_b_beg (i)), &
                    !        !(sum (tw_end (i, :)) + tw_b_end (i)), &
                    !        (sum (tw_end (i, :)) + tw_b_end (i) - sum (tw_beg (i, :)) - tw_b_beg (i)) / &
                    !        (sum (tw_beg (i, :)) + tw_b_beg (i))
                    !endif
                    inline_pbl%fast_tw_a_chg (i, j) = sum (tw_end (i, :)) - sum (tw_beg (i, :))
                    inline_pbl%fast_tw_b_chg (i, j) = tw_b_end (i) - tw_b_beg (i)
                    !print*, "PBL-FAST LOSS (%) : ", te_loss (i) / (sum (te_beg (i, :)) + te_b_beg (i)) * 100.0
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
!$OMP                                    ptop, inline_pbl, u_dt, v_dt, safety_check, &
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

            ! orographic gravity wave drag and mountain blocking main program
            call sa_gwd_oro (ie-is+1, km, uu, vv, ta, qv, abs (mdt), gsize, &
                inline_pbl%kpbl (is:ie, j), pi, dp, pm, pmk, zi, zm, &
                inline_gwd%hprime (is:ie, j), inline_gwd%oc (is:ie, j), inline_gwd%oa (is:ie, j, :), &
                inline_gwd%ol (is:ie, j, :), inline_gwd%theta (is:ie, j), inline_gwd%sigma (is:ie, j), &
                inline_gwd%gamma (is:ie, j), inline_gwd%elvmax (is:ie, j))

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
                    !    print*, "GWD-FAST TE: ", &
                    !        !(sum (te_beg (i, :)) + te_b_beg (i)), &
                    !        !(sum (te_end (i, :)) + te_b_end (i)), &
                    !        (sum (te_end (i, :)) + te_b_end (i) - sum (te_beg (i, :)) - te_b_beg (i)) / &
                    !        (sum (te_beg (i, :)) + te_b_beg (i))
                    !endif
                    inline_gwd%fast_te_a_chg (i, j) = sum (te_end (i, :)) - sum (te_beg (i, :))
                    inline_gwd%fast_te_b_chg (i, j) = te_b_end (i) - te_b_beg (i)
                    !if (abs (sum (tw_end (i, :)) + tw_b_end (i) - sum (tw_beg (i, :)) - tw_b_beg (i)) / &
                    !     (sum (tw_beg (i, :)) + tw_b_beg (i)) .gt. tw_err) then
                    !    print*, "GWD-FAST TW: ", &
                    !        !(sum (tw_beg (i, :)) + tw_b_beg (i)), &
                    !        !(sum (tw_end (i, :)) + tw_b_end (i)), &
                    !        (sum (tw_end (i, :)) + tw_b_end (i) - sum (tw_beg (i, :)) - tw_b_beg (i)) / &
                    !        (sum (tw_beg (i, :)) + tw_b_beg (i))
                    !endif
                    inline_gwd%fast_tw_a_chg (i, j) = sum (tw_end (i, :)) - sum (tw_beg (i, :))
                    inline_gwd%fast_tw_b_chg (i, j) = tw_b_end (i) - tw_b_beg (i)
                    !print*, "GWD-FAST LOSS (%) : ", te_loss (i) / (sum (te_beg (i, :)) + te_b_beg (i)) * 100.0
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
    ! pt conversion
    !-----------------------------------------------------------------------

    if (thermostruct%moist_kappa) then
       do k = 1, km
           do j = js, je
               do i = is, ie
                   pkz (i, j, k) = exp (cappa (i, j, k) * &
                       log (rrg * delp (i, j, k) / &
                       delz (i, j, k) * pt (i, j, k)))
                   pt (i, j, k) = pt (i, j, k) / pkz (i, j, k)
               enddo
           enddo
       enddo
    else
       do k = 1, km
           do j = js, je
               do i = is, ie
                   pkz (i, j, k) = exp (akap * &
                       log (rrg * delp (i, j, k) / &
                       delz (i, j, k) * pt (i, j, k)))
                   pt (i, j, k) = pt (i, j, k) / pkz (i, j, k)
               enddo
           enddo
       enddo
    endif

end subroutine fast_phys

end module fast_phys_mod
