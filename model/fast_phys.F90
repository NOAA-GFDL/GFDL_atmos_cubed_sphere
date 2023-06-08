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

    use constants_mod, only: rdgas, grav
    use fv_grid_utils_mod, only: cubed_to_latlon, update_dwinds_phys
    use fv_arrays_mod, only: fv_grid_type, fv_grid_bounds_type
    use mpp_domains_mod, only: domain2d, mpp_update_domains
    use tracer_manager_mod, only: get_tracer_index, get_tracer_names
    use field_manager_mod, only: model_atmos
    use gfdl_mp_mod, only: mtetw
    
    implicit none
    
    private

    real, parameter :: consv_min = 0.001

    public :: fast_phys

    ! -----------------------------------------------------------------------
    ! precision definition
    ! -----------------------------------------------------------------------
    
    integer, parameter :: r8 = 8 ! double precision
    
contains

subroutine fast_phys (is, ie, js, je, isd, ied, jsd, jed, km, npx, npy, nq, &
               c2l_ord, mdt, consv, akap, ptop, hs, te0_2d, u, v, w, pt, &
               delp, delz, q_con, cappa, q, pkz, r_vir, te_err, tw_err, &
               gridstruct, domain, bd, hydrostatic, do_adiabatic_init, &
               consv_checker, adj_mass_vmr)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------

    integer, intent (in) :: is, ie, js, je, isd, ied, jsd, jed, km, npx, npy, nq, c2l_ord, adj_mass_vmr

    logical, intent (in) :: hydrostatic, do_adiabatic_init, consv_checker

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

    type (fv_grid_bounds_type), intent (in) :: bd

    type (domain2d), intent (inout) :: domain


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

    do k = 1, km
        do j = js, je
            do i = is, ie
#ifdef MOIST_CAPPA
                pt (i, j, k) = pt (i, j, k) * exp (cappa (i, j, k) / (1. - cappa (i, j, k)) * &
                    log (rrg * delp (i, j, k) / delz (i, j, k) * pt (i, j, k)))
#else
                pt (i, j, k) = pt (i, j, k) * exp (akap / (1 - akap) * &
                    log (rrg * delp (i, j, k) / delz (i, j, k) * pt (i, j, k)))
#endif
            enddo
        enddo
    enddo

    !-----------------------------------------------------------------------
    ! pt conversion
    !-----------------------------------------------------------------------

    do k = 1, km
        do j = js, je
            do i = is, ie
#ifdef MOIST_CAPPA
                pkz (i, j, k) = exp (cappa (i, j, k) * &
                    log (rrg * delp (i, j, k) / &
                    delz (i, j, k) * pt (i, j, k)))
#else
                pkz (i, j, k) = exp (akap * &
                    log (rrg * delp (i, j, k) / &
                    delz (i, j, k) * pt (i, j, k)))
#endif
                pt (i, j, k) = pt (i, j, k) / pkz (i, j, k)
            enddo
        enddo
    enddo

end subroutine fast_phys

end module fast_phys_mod
