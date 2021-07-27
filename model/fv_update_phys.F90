
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

!>@brief The module 'fv_update_phys' applies physics tendencies consistent with
!! the FV3 discretization and definition of the prognostic variables.

module fv_update_phys_mod

! Modules Included:
! <table>
!   <tr>
!     <th>Module Name</th>
!     <th>Functions Included</th>
!   </tr>
!   <tr>
!     <td>boundary_mod</td>
!     <td>nested_grid_BC, extrapolation_BC</td>
!   </tr>
!   <tr>
!     <td>constants_mod</td>
!     <td>kappa, rdgas, rvgas, grav, cp_air, cp_vapor, pi=>pi_8, radius</td>
!   </tr>
!   <tr>
!     <td>field_manager_mod</td>
!     <td>MODEL_ATMOS</td>
!   </tr>
!   <tr>
!     <td>fv_arrays_mod</td>
!     <td>fv_flags_type, fv_nest_type, R_GRID</td>
!   </tr>
!   <tr>
!     <td>fv_diagnostics_mod</td>
!     <td>prt_maxmin</td>
!   </tr>
!   <tr>
!     <td>fv_eta_mod</td>
!     <td>get_eta_level</td>
!   </tr>
!   <tr>
!     <td>fv_mapz_mod</td>
!     <td>moist_cv, moist_cp</td>
!   </tr>
!   <tr>
!     <td>fv_mp_mod</td>
!     <td>start_group_halo_update, complete_group_halo_update,group_halo_update_type</td>
!   </tr>
!   <tr>
!     <td>fv_timing_mod</td>
!     <td>timing_on, timing_off</td>
!   </tr>
!   <tr>
!     <td>mpp_mod</td>
!     <td>FATAL, mpp_error,NOTE, WARNING</td>
!   </tr>
!   <tr>
!     <td>mpp_domains_mod</td>
!     <td>mpp_update_domains, domain2d</td>
!   </tr>
!   <tr>
!     <td>mpp_parameter_mod</td>
!     <td>AGRID_PARAM=>AGRID</td>
!   </tr>
!   <tr>
!     <td>time_manager_mod</td>
!     <td>time_type</td>
!   </tr>
!   <tr>
!     <td>tracer_manager_mod</td>
!     <td>get_tracer_index, adjust_mass, get_tracer_names</td>
!   </tr>
! </table>

  use constants_mod,      only: kappa, rdgas, rvgas, grav, cp_air, cp_vapor, pi=>pi_8, radius, TFREEZE
  use field_manager_mod,  only: MODEL_ATMOS
  use mpp_domains_mod,    only: mpp_update_domains, domain2d
  use mpp_parameter_mod,  only: AGRID_PARAM=>AGRID
  use mpp_mod,            only: FATAL, mpp_error
  use mpp_mod,            only: mpp_error, NOTE, WARNING, mpp_pe
  use time_manager_mod,   only: time_type
  use tracer_manager_mod, only: get_tracer_index, adjust_mass, get_tracer_names
  use fv_mp_mod,          only: start_group_halo_update, complete_group_halo_update
  use fv_mp_mod,          only: group_halo_update_type
  use fv_arrays_mod,      only: fv_flags_type, fv_nest_type, R_GRID, phys_diag_type, nudge_diag_type
  use boundary_mod,       only: nested_grid_BC
  use boundary_mod,       only: extrapolation_BC
  use fv_eta_mod,         only: get_eta_level
  use fv_timing_mod,      only: timing_on, timing_off
  use fv_diagnostics_mod, only: prt_maxmin, range_check
  use fv_mapz_mod,        only: moist_cv, moist_cp
#if defined (ATMOS_NUDGE)
  use atmos_nudge_mod,    only: get_atmos_nudge, do_ps
#elif defined (CLIMATE_NUDGE)
  use fv_climate_nudge_mod, only: fv_climate_nudge, do_ps
#elif defined (ADA_NUDGE)
  use fv_ada_nudge_mod,   only: fv_ada_nudge
#else
  use fv_nwp_nudge_mod,   only: fv_nwp_nudge
#endif
  use fv_arrays_mod,      only: fv_grid_type, fv_nest_type, fv_grid_bounds_type
  use fv_grid_utils_mod,  only: cubed_to_latlon, update_dwinds_phys, update2d_dwinds_phys
  use fv_nesting_mod,     only: set_physics_BCs
  use sat_vapor_pres_mod, only: tcmin, tcmax
#ifdef MULTI_GASES
  use multi_gases_mod,    only: virq, virqd, vicpqd, vicvqd, num_gas
#endif

  implicit none

  public :: fv_update_phys, del2_phys
  real,parameter:: con_cp  = cp_air
  real, parameter :: tmax = 330

  contains

  subroutine fv_update_phys ( dt, is, ie, js, je, isd, ied, jsd, jed, ng, nq,     &
                              u, v, w, delp, pt, q, qdiag, ua, va, ps, pe,  peln, pk, pkz,  &
                              ak, bk, phis, u_srf, v_srf, ts, delz, hydrostatic,  &
                              u_dt, v_dt, t_dt, moist_phys, Time, nudge,    &
                              gridstruct, lona, lata, npx, npy, npz, flagstruct,  &
                              neststruct, bd, domain, ptop, phys_diag, &
                              nudge_diag, q_dt)
    real, intent(in)   :: dt, ptop
    integer, intent(in):: is,  ie,  js,  je, ng
    integer, intent(in):: isd, ied, jsd, jed
    integer, intent(in):: nq            ! tracers modified by physics
                                        ! ncnst is the total number of tracers
    logical, intent(in):: moist_phys
    logical, intent(in):: hydrostatic
    logical, intent(in):: nudge

    type (time_type), intent(in) :: Time

    real, intent(in), dimension(npz+1):: ak, bk
    real, intent(in) :: phis(isd:ied,jsd:jed)
    real, intent(inout):: delz(is:,js:,1:)

! optional arguments for atmospheric nudging
    real, intent(in), dimension(isd:ied,jsd:jed), optional ::   &
                                lona, lata   !< A-grid (physics) lon and lat

! Winds on lat-lon grid:
    real, intent(inout), dimension(isd:ied,jsd:jed,npz):: ua, va
    real, intent(inout), dimension(isd:   ,jsd:   ,1: ):: w

! Tendencies from Physics:
    real, intent(inout), dimension(isd:ied,jsd:jed,npz):: u_dt, v_dt
    real, intent(inout):: t_dt(is:ie,js:je,npz)
    real, intent(inout), optional :: q_dt(is:ie,js:je,npz,nq)
    type(phys_diag_type), intent(inout) :: phys_diag
    type(nudge_diag_type), intent(inout) :: nudge_diag

! Saved Bottom winds for GFDL Physics Interface
    real, intent(out), dimension(is:ie,js:je):: u_srf, v_srf, ts

    type(fv_flags_type) :: flagstruct
    type(fv_grid_bounds_type), intent(IN)  :: bd
    type(domain2d), intent(INOUT) :: domain

    real, intent(inout):: u(isd:ied  ,jsd:jed+1,npz)  !< D grid zonal wind (m/s)
    real, intent(inout):: v(isd:ied+1,jsd:jed  ,npz)  !< D grid meridional wind (m/s)
    real, intent(inout), dimension(isd:ied,jsd:jed,npz):: pt, delp
    real, intent(inout):: q(isd:ied,jsd:jed,npz,nq)   !< specific humidity and constituents
    real, intent(inout):: qdiag(isd:ied,jsd:jed,npz,nq+1:flagstruct%ncnst) !< diagnostic tracers

!-----------------------------------------------------------------------
! Auxilliary pressure arrays:
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
    real, intent(inout):: ps  (isd:ied  ,jsd:jed)           !< Surface pressure (pascal)
    real, intent(inout):: pe  (is-1:ie+1, npz+1,js-1:je+1)  !< edge pressure (pascal)
    real, intent(inout):: pk  (is:ie,js:je  , npz+1)        !< pe**cappa
    real, intent(inout):: peln(is:ie,npz+1,js:je)           !< ln(pe)
    real, intent(inout):: pkz (is:ie,js:je,npz)             !< finite-volume mean pk
    real, parameter:: tice = 273.16

    type(fv_grid_type) :: gridstruct
    type(fv_nest_type) :: neststruct

    real :: q_dt_nudge(is:ie,js:je,npz,nq)

    integer, intent(IN) :: npx, npy, npz

!***********
! Haloe Data
!***********
    real, parameter::    q1_h2o = 2.2E-6
    real, parameter::    q7_h2o = 3.8E-6
    real, parameter::  q100_h2o = 3.8E-6
    real, parameter:: q1000_h2o = 3.1E-6
    real, parameter:: q2000_h2o = 2.8E-6
    real, parameter:: q3000_h2o = 3.0E-6

! Local arrays:
    real  ps_dt(is:ie,js:je)
    real  cvm(is:ie), qc(is:ie)
    real  phalf(npz+1), pfull(npz)

#ifdef MULTI_GASES
    integer :: nn, nm
#endif

    type(group_halo_update_type), save :: i_pack(2)
    integer  i, j, k, m, n, nwat
    integer  sphum, liq_wat, ice_wat, cld_amt   ! GFDL AM physics
    integer  rainwat, snowwat, graupel          ! GFDL Cloud Microphysics
    integer  w_diff                             ! w-tracer for PBL diffusion
    real:: qstar, dbk, rdg, zvir, p_fac, cv_air, gama_dt, tbad
    logical :: bad_range

!f1p
!account for change in air molecular weight because of H2O change
    logical, dimension(nq) :: conv_vmr_mmr
    real                   :: adj_vmr(is:ie,js:je,npz)
    character(len=32)      :: tracer_units, tracer_name

    cv_air = cp_air - rdgas ! = rdgas * (7/2-1) = 2.5*rdgas=717.68

    rdg = -rdgas / grav

    nwat = flagstruct%nwat

    if ( moist_phys .or. nwat/=0 ) then
           zvir = rvgas/rdgas - 1.
    else
           zvir = 0.
    endif

!f1p
    conv_vmr_mmr(1:nq) = .false.
    if (flagstruct%adj_mass_vmr) then
    do m=1,nq
       call get_tracer_names (MODEL_ATMOS, m, name = tracer_name,  &
            units = tracer_units)
       if ( trim(tracer_units) .eq. 'vmr' ) then
          conv_vmr_mmr(m) = .true.
       else
          conv_vmr_mmr(m) = .false.
       end if
    end do
    end if

    sphum   = get_tracer_index (MODEL_ATMOS, 'sphum')
    liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
    ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
    rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat')
    snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat')
    graupel = get_tracer_index (MODEL_ATMOS, 'graupel')
    cld_amt = get_tracer_index (MODEL_ATMOS, 'cld_amt')

    if ( .not. hydrostatic ) then
        w_diff = get_tracer_index (MODEL_ATMOS, 'w_diff')
!       if ( w_diff<8 ) call mpp_error(FATAL,'W_tracer index error')
    else
        w_diff = 0
    endif

    if ( .not. hydrostatic .and. .not. flagstruct%phys_hydrostatic .and. nwat == 0 ) then
       gama_dt = dt*cp_air/cv_air
    else
       gama_dt = -1.e24
    endif

    if ( flagstruct%fv_debug ) then
      call prt_maxmin('delp_b_update', delp, is, ie, js,  je, ng, npz, 0.01)
      if (present(q_dt)) then
        do m=1,nq
          call prt_maxmin('q_dt', q_dt(is,js,1,m), is, ie, js, je, 0, npz, 1.)
        enddo
      endif
      call prt_maxmin('u_dt', u_dt, is, ie, js,  je, ng, npz, 1.)
      call prt_maxmin('v_dt', v_dt, is, ie, js,  je, ng, npz, 1.)
      call prt_maxmin('T_dt', t_dt, is, ie, js,  je, 0, npz, 1.)
    endif

    call get_eta_level(npz, 1.0E5, pfull, phalf, ak, bk)
#ifdef MULTI_GASES
        nm = max(           nwat,num_gas)
        nn = max(flagstruct%nwat,num_gas)
#endif

    if (size(neststruct%child_grids) > 1) then
       call set_physics_BCs(ps, u_dt, v_dt, flagstruct, gridstruct, neststruct, npx, npy, npz, ng, ak, bk, bd)
    endif

    if (allocated(phys_diag%phys_t_dt)) phys_diag%phys_t_dt = pt(is:ie,js:je,:)
    if (present(q_dt)) then
       if (allocated(phys_diag%phys_qv_dt)) phys_diag%phys_qv_dt = q(is:ie,js:je,:,sphum)
       if (allocated(phys_diag%phys_ql_dt)) then
          if (liq_wat < 0) call mpp_error(FATAL, " phys_ql_dt needs at least one liquid water tracer defined")
          phys_diag%phys_ql_dt = q(is:ie,js:je,:,liq_wat)
          if (rainwat > 0) phys_diag%phys_ql_dt = q(is:ie,js:je,:,rainwat) + phys_diag%phys_ql_dt
       endif
       if (allocated(phys_diag%phys_qi_dt)) then
          if (ice_wat < 0) then
             call mpp_error(WARNING, " phys_qi_dt needs at least one ice water tracer defined")
             phys_diag%phys_qi_dt = 0.
          endif
          phys_diag%phys_qi_dt = q(is:ie,js:je,:,ice_wat)
          if (snowwat > 0) phys_diag%phys_qi_dt = q(is:ie,js:je,:,snowwat) + phys_diag%phys_qi_dt
          if (graupel > 0) phys_diag%phys_qi_dt = q(is:ie,js:je,:,graupel) + phys_diag%phys_qi_dt
       endif
    endif

!$OMP parallel do default(none) &
!$OMP             shared(is,ie,js,je,npz,flagstruct,pfull,q_dt,sphum,q,qdiag,  &
!$OMP                    nq,w_diff,dt,nwat,liq_wat,rainwat,ice_wat,snowwat,    &
!$OMP                    graupel,delp,cld_amt,hydrostatic,pt,t_dt,delz,adj_vmr,&
!$OMP                    gama_dt,cv_air,ua,u_dt,va,v_dt,isd,ied,jsd,jed,       &
#ifdef MULTI_GASES
!$OMP                    nn, nm,                                               &
#endif
!$OMP                    conv_vmr_mmr,pe,ptop,gridstruct,phys_diag)                                         &
!$OMP             private(cvm, qc, qstar, ps_dt, p_fac,tbad)
    do k=1, npz

      if (present(q_dt)) then

        if (flagstruct%tau_h2o<0.0 .and. pfull(k) < 100.E2 ) then
! Wipe the stratosphere clean:
! This should only be used for initialization from a bad model state
           p_fac = -flagstruct%tau_h2o*86400.
           do j=js,je
              do i=is,ie
                 q_dt(i,j,k,sphum) = q_dt(i,j,k,sphum) + (3.E-6-q(i,j,k,sphum))/p_fac
              enddo
           enddo
        elseif ( flagstruct%tau_h2o>0.0 .and. pfull(k) < 3000. ) then
! Do idealized Ch4 chemistry

           if ( pfull(k) < 1. ) then
               qstar = q1_h2o
               p_fac = 0.2 * flagstruct%tau_h2o*86400.
           elseif ( pfull(k) <   7. .and. pfull(k) >=    1. ) then
               qstar = q1_h2o + (q7_h2o-q1_h2o)*log(pfull(k)/1.)/log(7.)
               p_fac = 0.3 * flagstruct%tau_h2o*86400.
           elseif ( pfull(k) <  100. .and. pfull(k) >=    7. ) then
               qstar = q7_h2o + (q100_h2o-q7_h2o)*log(pfull(k)/7.)/log(100./7.)
               p_fac = 0.4 * flagstruct%tau_h2o*86400.
           elseif ( pfull(k) < 1000. .and. pfull(k) >=  100. ) then
               qstar = q100_h2o + (q1000_h2o-q100_h2o)*log(pfull(k)/1.E2)/log(10.)
               p_fac = 0.5 * flagstruct%tau_h2o*86400.
           elseif ( pfull(k) < 2000. .and. pfull(k) >= 1000. ) then
               qstar = q1000_h2o + (q2000_h2o-q1000_h2o)*log(pfull(k)/1.E3)/log(2.)
               p_fac = 0.75 * flagstruct%tau_h2o*86400.
           else
               qstar = q3000_h2o
               p_fac = flagstruct%tau_h2o*86400.
           endif

           do j=js,je
              do i=is,ie
                 q_dt(i,j,k,sphum) = q_dt(i,j,k,sphum) + (qstar-q(i,j,k,sphum))/p_fac
              enddo
           enddo
        endif

!----------------
! Update tracers:
!----------------
        do m=1,nq
          if( m /= w_diff ) then
            do j=js,je
              do i=is,ie
                q(i,j,k,m) = q(i,j,k,m) + dt*q_dt(i,j,k,m)
              enddo
            enddo
          endif
        enddo

!--------------------------------------------------------
! Adjust total air mass due to changes in water substance
!--------------------------------------------------------
        do j=js,je
          do i=is,ie
#ifdef MULTI_GASES
            ps_dt(i,j)  = 1. + dt*sum(q_dt(i,j,k,1:nm))
#else
            ps_dt(i,j)  = 1. + dt*sum(q_dt(i,j,k,1:nwat))
#endif
            delp(i,j,k) = delp(i,j,k) * ps_dt(i,j)
            if (flagstruct%adj_mass_vmr) then
#ifdef MULTI_GASES
               adj_vmr(i,j,k) =                          &
                    (ps_dt(i,j) - sum(q(i,j,k,1:nn))) /  &
                    (1.d0       - sum(q(i,j,k,1:nn)))
#else
               adj_vmr(i,j,k) =                          &
                    (ps_dt(i,j) - sum(q(i,j,k,1:flagstruct%nwat))) /  &
                    (1.d0       - sum(q(i,j,k,1:flagstruct%nwat)))
#endif
            end if
          enddo
        enddo

!-----------------------------------------
! Adjust mass mixing ratio of all tracers
!-----------------------------------------
        if ( nwat /=0 ) then
          do m=1,flagstruct%ncnst
!-- check to query field_table to determine if tracer needs mass adjustment
            if( m /= cld_amt .and. m /= w_diff .and. adjust_mass(MODEL_ATMOS,m)) then
              if (m <= nq)  then
                q(is:ie,js:je,k,m) = q(is:ie,js:je,k,m) / ps_dt(is:ie,js:je)
                if (conv_vmr_mmr(m)) &
                   q(is:ie,js:je,k,m) = q(is:ie,js:je,k,m) * adj_vmr(is:ie,js:je,k)
              else
                qdiag(is:ie,js:je,k,m) = qdiag(is:ie,js:je,k,m) / ps_dt(is:ie,js:je)
              endif
            endif
          enddo
        endif

      endif ! present(q_dt)

      if ( hydrostatic ) then
         do j=js,je
            call moist_cp(is,ie,isd,ied,jsd,jed, npz, j, k, nwat, sphum, liq_wat, rainwat,    &
                          ice_wat, snowwat, graupel, q, qc, cvm, pt(is:ie,j,k) )
            do i=is,ie
               pt(i,j,k) = pt(i,j,k) + t_dt(i,j,k)*dt*con_cp/cvm(i)
            enddo
         enddo
      else
         if ( flagstruct%phys_hydrostatic ) then
! Constant pressure
             do j=js,je
                call moist_cp(is,ie,isd,ied,jsd,jed, npz, j, k, nwat, sphum, liq_wat, rainwat,    &
                              ice_wat, snowwat, graupel, q, qc, cvm, pt(is:ie,j,k) )
                do i=is,ie
                   delz(i,j,k) = delz(i,j,k) / pt(i,j,k)
                   pt(i,j,k) = pt(i,j,k) + t_dt(i,j,k)*dt*con_cp/cvm(i)
                   delz(i,j,k) = delz(i,j,k) * pt(i,j,k)
                enddo
             enddo
         else
            !NOTE: only works for either no physics or GFDL MP
            if (nwat == 0) then
               do j=js,je
                  do i=is,ie
                     pt(i,j,k) = pt(i,j,k) + t_dt(i,j,k)*gama_dt
                  enddo
               enddo
            else
               do j=js,je
                  call moist_cv(is,ie,isd,ied,jsd,jed, npz, j, k, nwat, sphum, liq_wat, rainwat,    &
                                ice_wat, snowwat, graupel, q, qc, cvm, pt(is:ie,j,k))
                  do i=is,ie
                     pt(i,j,k) = pt(i,j,k) + t_dt(i,j,k)*dt*con_cp/cvm(i)
                  enddo
               enddo
            endif
         endif
      endif

#ifndef GFS_PHYS
      do j=js,je
         do i=is,ie
            ua(i,j,k) = ua(i,j,k) + dt*u_dt(i,j,k)
            va(i,j,k) = va(i,j,k) + dt*v_dt(i,j,k)
         enddo
      enddo
#endif

   enddo ! k-loop

   if (allocated(phys_diag%phys_t_dt)) phys_diag%phys_t_dt = (pt(is:ie,js:je,:) - phys_diag%phys_t_dt) / dt
   if (present(q_dt)) then
      if (allocated(phys_diag%phys_qv_dt)) phys_diag%phys_qv_dt = (q(is:ie,js:je,:,sphum) - phys_diag%phys_qv_dt) / dt
      if (allocated(phys_diag%phys_ql_dt)) then
         if (liq_wat < 0) call mpp_error(FATAL, " phys_ql_dt needs at least one liquid water tracer defined")
         phys_diag%phys_ql_dt = q(is:ie,js:je,:,liq_wat) - phys_diag%phys_qv_dt
         if (rainwat > 0) phys_diag%phys_ql_dt = q(is:ie,js:je,:,rainwat) + phys_diag%phys_ql_dt
         phys_diag%phys_ql_dt = phys_diag%phys_ql_dt / dt
      endif
      if (allocated(phys_diag%phys_qi_dt)) then
         if (ice_wat < 0) then
            call mpp_error(WARNING, " phys_qi_dt needs at least one ice water tracer defined")
            phys_diag%phys_qi_dt = 0.
         endif
         phys_diag%phys_qi_dt = q(is:ie,js:je,:,ice_wat) - phys_diag%phys_qi_dt
         if (snowwat > 0) phys_diag%phys_qi_dt = q(is:ie,js:je,:,snowwat) + phys_diag%phys_qi_dt
         if (graupel > 0) phys_diag%phys_qi_dt = q(is:ie,js:je,:,graupel) + phys_diag%phys_qi_dt
         phys_diag%phys_qi_dt = phys_diag%phys_qi_dt / dt
      endif
   endif

   if ( flagstruct%range_warn ) then
      call range_check('PT UPDATE', pt, is, ie, js, je, ng, npz, gridstruct%agrid,    &
                          tcmin+TFREEZE, tcmax+TFREEZE, bad_range, Time)
      if (bad_range) then
         do k=1,npz
         do j=js,je
         do i=is,ie
            if (pt(i,j,k) < tcmin+TFREEZE .or. pt(i,j,k) > tcmax+TFREEZE) then
               write(*,*) 'PT UPDATE: ', t_dt(i,j,k)*dt, i,j,k, gridstruct%agrid(i,j,:)
            endif
         enddo
         enddo
         enddo
      endif
   endif

! [delp, (ua, va), pt, q] updated. Perform nudging if requested
!------- nudging of atmospheric variables toward specified data --------

    ps_dt(:,:) = 0.

    if ( nudge ) then
#if defined (ATMOS_NUDGE)
!--------------------------------------------
! All fields will be updated; tendencies added
!--------------------------------------------
        call get_atmos_nudge ( Time, dt, is, ie, js, je,    &
             npz, ng, ps(is:ie,js:je), ua(is:ie, js:je,:),  &
             va(is:ie,js:je,:), pt(is:ie,js:je,:),          &
             q(is:ie,js:je,:,:), ps_dt(is:ie,js:je), u_dt(is:ie,js:je,:),  &
             v_dt(is:ie,js:je,:), t_dt(is:ie,js:je,:), &
             q_dt_nudge(is:ie,js:je,:,:) )

!--------------
! Update delp
!--------------
        if (do_ps) then
!$omp parallel do default(none) shared(is,ie,js,je,npz,bk,delp,ps_dt) private(dbk)
            do k=1,npz
               dbk = dt * (bk(k+1) - bk(k))
               do j=js,je
                  do i=is,ie
                     delp(i,j,k) = delp(i,j,k) + dbk*ps_dt(i,j)
                  enddo
               enddo
            enddo
#elif defined (CLIMATE_NUDGE)
!--------------------------------------------
! All fields will be updated; tendencies added
!--------------------------------------------
             lona(is:ie,js:je), lata(is:ie,js:je), phis(is:ie,js:je), &
             ptop, ak, bk, &
             ps(is:ie,js:je), ua(is:ie,js:je,:), va(is:ie,js:je,:), &
             pt(is:ie,js:je,:), q(is:ie,js:je,:,sphum:sphum),   &
             ps_dt(is:ie,js:je), u_dt(is:ie,js:je,:),  &
             v_dt(is:ie,js:je,:), t_dt(is:ie,js:je,:), &
             q_dt_nudge(is:ie,js:je,:,sphum:sphum) )

!--------------
! Update delp
!--------------
        if (do_ps) then
!$omp parallel do default(none) shared(is,ie,js,je,npz,dt,bk,delp,ps_dt) private(dbk)
            do k=1,npz
              dbk = dt * (bk(k+1) - bk(k))
               do j=js,je
                  do i=is,ie
                     delp(i,j,k) = delp(i,j,k) + dbk*ps_dt(i,j)
                   enddo
               enddo
            enddo
        endif
#elif defined (ADA_NUDGE)
! All fields will be updated except winds; wind tendencies added
!$omp parallel do default(shared)
        do j=js,je
          do k=2,npz+1
            do i=is,ie
              pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
            enddo
          enddo
          do i=is,ie
            ps(i,j) = pe(i,npz+1,j)
          enddo
        enddo
        call fv_ada_nudge ( Time, dt, npx, npy, npz,  ps_dt, u_dt, v_dt, t_dt, q_dt_nudge,   &
                            zvir, ptop, ak, bk, ts, ps, delp, ua, va, pt,    &
                            nwat, q,  phis, gridstruct, bd, domain )
#else
! All fields will be updated except winds; wind tendencies added
!$OMP parallel do default(none) shared(is,ie,js,je,npz,pe,delp,ps)
        do j=js,je
          do k=2,npz+1
            do i=is,ie
              pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
            enddo
          enddo
          do i=is,ie
            ps(i,j) = pe(i,npz+1,j)
          enddo
        enddo
        call fv_nwp_nudge ( Time, dt, npx, npy, npz,  ps_dt, u_dt, v_dt, t_dt, q_dt_nudge,   &
                            zvir, ptop, ak, bk, ts, ps, delp, ua, va, pt,    &
                            nwat, q,  phis, gridstruct, bd, domain )
#endif

  endif         ! end nudging

  if ( .not.flagstruct%dwind_2d ) then

                                                            call timing_on('COMM_TOTAL')
      if ( gridstruct%square_domain ) then
           call start_group_halo_update(i_pack(1), u_dt, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.false.)
           call start_group_halo_update(i_pack(1), v_dt, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
      else
           call start_group_halo_update(i_pack(1), u_dt, domain, complete=.false.)
           call start_group_halo_update(i_pack(1), v_dt, domain, complete=.true.)
      endif
                                                           call timing_off('COMM_TOTAL')
  endif

!----------------------------------------
! Update pe, peln, pkz, and surface winds
!----------------------------------------
  if ( flagstruct%fv_debug ) then
       call prt_maxmin('PS_b_update',     ps, is, ie, js,  je, ng,   1, 0.01)
       call prt_maxmin('delp_a_update', delp, is, ie, js,  je, ng, npz, 0.01)
  endif

!$OMP parallel do default(none) shared(is,ie,js,je,npz,pe,delp,peln,pk,ps,u_srf,v_srf, &
#ifdef MULTI_GASES
!$OMP                                  q,                                              &
#endif
!$OMP                                  ua,va,pkz,hydrostatic)
   do j=js,je
      do k=2,npz+1
         do i=is,ie
              pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
            peln(i,k,j) = log( pe(i,k,j) )
              pk(i,j,k) = exp( kappa*peln(i,k,j) )
         enddo
      enddo

      do i=is,ie
            ps(i,j) = pe(i,npz+1,j)
         u_srf(i,j) = ua(i,j,npz)
         v_srf(i,j) = va(i,j,npz)
      enddo

      if ( hydrostatic ) then
         do k=1,npz
            do i=is,ie
               pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k))/(kappa*(peln(i,k+1,j)-peln(i,k,j)))
#ifdef MULTI_GASES
               pkz(i,j,k) = exp( virqd(q(i,j,k,:))/vicpqd(q(i,j,k,:)) * log( pkz(i,j,k) ) )
#endif
            enddo
         enddo
      endif
   enddo      ! j-loop

                                                    call timing_on(' Update_dwinds')
  if ( flagstruct%dwind_2d ) then
    call update2d_dwinds_phys(is, ie, js, je, isd, ied, jsd, jed, dt, u_dt, v_dt, u, v, gridstruct, &
         npx,npy,npz,domain)
  else

     !I have not seen dwind_2d be used for anything; so we will only handle nesting assuming dwind_2d == .false.

    call timing_on('COMM_TOTAL')

    call complete_group_halo_update(i_pack(1), domain)

    call timing_off('COMM_TOTAL')
!
! for regional grid need to set values for u_dt and v_dt at the edges.
! Note from Lucas:The physics only operates on the compute domain.
! One snag is that in fv_update_phys.F90 u_dt and v_dt from the physics need to be interpolated to the D-grids,
! which requires BCs for u_dt and v_dt. For the nested grid I can simply get the BCs from the coarse grid, but
! in your case I would recommend just setting the boundary conditions to 0 or to constant values (ie the value
! of the cell closest to the boundary).
    if (gridstruct%regional) then
     if (is == 1) then
      do k=1,npz
       do j = js,je
        u_dt(is-1,j,k) = u_dt(is,j,k)
        v_dt(is-1,j,k) = v_dt(is,j,k)
       enddo
      enddo
     endif
     if (ie == npx) then
      do k=1,npz
       do j = js,je
        u_dt(ie+1,j,k) = u_dt(ie,j,k)
        v_dt(ie+1,j,k) = v_dt(ie,j,k)
       enddo
      enddo
     endif
     if (js == 1) then
      do k=1,npz
       do i = is,ie
        u_dt(i,js-1,k) = u_dt(i,js,k)
        v_dt(i,js-1,k) = v_dt(i,js,k)
       enddo
      enddo
     endif
     if (je == npy) then
      do k=1,npz
       do i = is,ie
        u_dt(i,je+1,k) = u_dt(i,je,k)
        v_dt(i,je+1,k) = v_dt(i,je,k)
       enddo
      enddo
     endif
!
! corners
!
     do k=1,npz
      if (is == 1 .and. js == 1) then
       u_dt(is-1,js-1,k) = u_dt(is,js,k)
       v_dt(is-1,js-1,k) = v_dt(is,js,k)
      elseif (is == 1 .and. je == npy) then
       u_dt(is-1,je+1,k) = u_dt(is,je,k)
       v_dt(is-1,je+1,k) = v_dt(is,je,k)
      elseif (ie == npx .and. js == 1) then
       u_dt(ie+1,js-1,k) = u_dt(ie,je,k)
       v_dt(ie+1,js-1,k) = v_dt(ie,je,k)
      elseif (ie == npx .and. je == npy) then
       u_dt(ie+1,je+1,k) = u_dt(ie,je,k)
       v_dt(ie+1,je+1,k) = v_dt(ie,je,k)
      endif
     enddo
    endif !regional
!
    call update_dwinds_phys(is, ie, js, je, isd, ied, jsd, jed, dt, u_dt, v_dt, u, v, gridstruct, npx, npy, npz, domain)
 endif
                                                    call timing_off(' Update_dwinds')
#ifdef GFS_PHYS
    call cubed_to_latlon(u, v, ua, va, gridstruct, &
         npx, npy, npz, 1, gridstruct%grid_type, domain, gridstruct%bounded_domain, flagstruct%c2l_ord, bd)
#endif

  if ( flagstruct%fv_debug ) then
       call prt_maxmin('PS_a_update', ps, is, ie, js, je, ng,   1, 0.01)
  endif

  if (allocated(phys_diag%phys_u_dt)) then
     do k=1,npz
     do j=js,je
     do i=is,ie
        phys_diag%phys_u_dt(i,j,k) = u_dt(i,j,k)
     enddo
     enddo
     enddo
  endif
  if (allocated(phys_diag%phys_v_dt)) then
     do k=1,npz
     do j=js,je
     do i=is,ie
        phys_diag%phys_v_dt(i,j,k) = v_dt(i,j,k)
     enddo
     enddo
     enddo
  endif

  end subroutine fv_update_phys

!>@brief The subroutine 'del2_phys' is for filtering the physics tendency.
  subroutine del2_phys(qdt, delp, gridstruct, cd, npx, npy, km, is, ie, js, je, &
                       isd, ied, jsd, jed, ngc, domain)
! This routine is for filtering the physics tendency
   integer, intent(in):: npx, npy, km
   integer, intent(in):: is, ie, js, je, isd, ied, jsd, jed, ngc
   real,    intent(in):: cd            !< cd = K * da_min;   0 < K < 0.25
   real, intent(in   ):: delp(isd:ied,jsd:jed,km)
   real, intent(inout):: qdt(is-ngc:ie+ngc,js-ngc:je+ngc,km)
   type(fv_grid_type), intent(IN), target :: gridstruct
   type(domain2d), intent(INOUT) :: domain

   real, pointer, dimension(:,:) :: rarea, dx, dy, sina_u, sina_v, rdxc, rdyc
   real, pointer, dimension(:,:,:) :: sin_sg
!
   real :: q(isd:ied,jsd:jed,km)
   real :: fx(is:ie+1,js:je), fy(is:ie,js:je+1)
   real :: mask(is:ie+1,js:je+1)
   real :: f1(is:ie+1), f2(js:je+1)
   real :: damp
   integer i,j,k

   rarea  => gridstruct%rarea
   dx     => gridstruct%dx
   dy     => gridstruct%dy
   sina_u => gridstruct%sina_u
   sina_v => gridstruct%sina_v
   rdxc   => gridstruct%rdxc
   rdyc   => gridstruct%rdyc
   sin_sg => gridstruct%sin_sg

! Applying mask to cd, the damping coefficient?
   damp = 0.25 * cd * gridstruct%da_min

! Mask defined at corners

!$OMP parallel do default(none) shared(is,ie,f1,npx)
   do i=is,ie+1
      f1(i) = (1. - sin(real(i-1)/real(npx-1)*pi))**2
   enddo

!$OMP parallel do default(none) shared(is,ie,js,je,f1,f2,npy,mask,damp)
   do j=js,je+1
      f2(j) = (1. - sin(real(j-1)/real(npy-1)*pi))**2
      do i=is,ie+1
         mask(i,j) = damp * (f1(i) + f2(j))
      enddo
   enddo

! mass weighted tendency from physics is filtered

!$OMP parallel do default(none) shared(is,ie,js,je,km,q,qdt,delp)
   do k=1,km
      do j=js,je
         do i=is,ie
            q(i,j,k) = qdt(i,j,k)*delp(i,j,k)
         enddo
      enddo
   enddo
                     call timing_on('COMM_TOTAL')
   call mpp_update_domains(q, domain, complete=.true.)
                     call timing_off('COMM_TOTAL')

!$OMP parallel do default(none) shared(is,ie,js,je,km,mask,dy,sina_u,q,rdxc,gridstruct, &
!$OMP                                  sin_sg,npx,dx,npy,rdyc,sina_v,qdt,rarea,delp)    &
!$OMP                          private(fx, fy)
   do k=1,km
      do j=js,je
         do i=is,ie+1
            fx(i,j) = &
                 (mask(i,j)+mask(i,j+1))*dy(i,j)*sina_u(i,j)* &
                 (q(i-1,j,k)-q(i,j,k))*rdxc(i,j)
         enddo
         if (is == 1 .and. .not. gridstruct%bounded_domain)   fx(i,j) = &
              (mask(is,j)+mask(is,j+1))*dy(is,j)*(q(is-1,j,k)-q(is,j,k))*rdxc(is,j)* &
            0.5*(sin_sg(1,j,1) + sin_sg(0,j,3))
         if (ie+1==npx .and. .not. gridstruct%bounded_domain) fx(i,j) = &
              (mask(ie+1,j)+mask(ie+1,j+1))*dy(ie+1,j)*(q(ie,j,k)-q(ie+1,j,k))*rdxc(ie+1,j)* &
            0.5*(sin_sg(npx,j,1) + sin_sg(npx-1,j,3))
      enddo
      do j=js,je+1
         if ((j == 1 .OR. j == npy) .and. .not. gridstruct%bounded_domain) then
            do i=is,ie
               fy(i,j) = (mask(i,j)+mask(i+1,j))*dx(i,j)*&
                    (q(i,j-1,k)-q(i,j,k))*rdyc(i,j) &
                    *0.5*(sin_sg(i,j,2) + sin_sg(i,j-1,4) )
            enddo
         else
            do i=is,ie
               fy(i,j) = (mask(i,j)+mask(i+1,j))*dx(i,j)*sina_v(i,j)*&
                    (q(i,j-1,k)-q(i,j,k))*rdyc(i,j)
            enddo
         end if
      enddo
      do j=js,je
         do i=is,ie
            qdt(i,j,k) = qdt(i,j,k) + rarea(i,j)*(fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))/delp(i,j,k)
         enddo
      enddo
   enddo

  end subroutine del2_phys

end module fv_update_phys_mod
