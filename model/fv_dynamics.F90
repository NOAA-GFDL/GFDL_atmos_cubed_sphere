module fv_dynamics_mod
   use constants_mod,       only: grav, pi, radius, hlv, rdgas    ! latent heat of water vapor
   use dyn_core_mod,        only: dyn_core, del2_cubed
   use fv_mapz_mod,         only: compute_total_energy, Lagrangian_to_Eulerian
   use fv_tracer2d_mod,     only: tracer_2d, tracer_2d_1L, tracer_2d_nested
   use fv_grid_utils_mod,   only: cubed_to_latlon, c2l_ord2
   use fv_mp_mod,           only: is_master
   use fv_timing_mod,       only: timing_on, timing_off
   use diag_manager_mod,    only: send_data
   use fv_diagnostics_mod,  only: fv_time, prt_maxmin, range_check
   use fv_diagnostics_mod,  only: prt_mass
   use mpp_domains_mod,     only: DGRID_NE, CGRID_NE, mpp_update_domains, domain2D
   use mpp_domains_mod,     only: mpp_start_update_domains, mpp_complete_update_domains
   use field_manager_mod,   only: MODEL_ATMOS
   use tracer_manager_mod,  only: get_tracer_index
   use fv_sg_mod,           only: neg_adj3
   use fv_nesting_mod,      only: setup_nested_grid_BCs
   use fv_arrays_mod,       only: fv_grid_type, fv_flags_type, fv_atmos_type, fv_nest_type, fv_diag_type, fv_grid_bounds_type

implicit none
   logical :: RF_initialized = .false.
   logical :: bad_range
   real, allocatable ::  rf(:), rw(:)
   integer :: kmax=1
   !Arrays for global grid total energy, used for grid nesting
   real, allocatable :: te_2d_coarse(:,:)
   real, allocatable :: dp1_coarse(:,:,:)
private
public :: fv_dynamics

!---- version number -----
   character(len=128) :: version = '$Id: fv_dynamics.F90,v 20.0 2013/12/13 23:04:25 fms Exp $'
   character(len=128) :: tagname = '$Name: tikal_201409 $'

contains

!-----------------------------------------------------------------------
!     fv_dynamics :: FV dynamical core driver
!-----------------------------------------------------------------------
 
  subroutine fv_dynamics(npx, npy, npz, nq,  ng, bdt, consv_te, fill,               &
                        reproduce_sum, kappa, cp_air, zvir, ptop, ks, ncnst, n_split,     &
                        q_split, u, v, w, delz, hydrostatic, pt, delp, q,   &
#ifdef PKC
                        ps, pe, pk, peln, pkz, pkc, phis, omga, ua, va, uc, vc,          &
#else
                        ps, pe, pk, peln, pkz, phis, omga, ua, va, uc, vc,          &
#endif
                        ak, bk, mfx, mfy, cx, cy, ze0, hybrid_z, &
                        gridstruct, flagstruct, neststruct, idiag, bd, &
                        parent_grid, domain, time_total)

    real, intent(IN) :: bdt  ! Large time-step
    real, intent(IN) :: consv_te
    real, intent(IN) :: kappa, cp_air
    real, intent(IN) :: zvir, ptop
    real, intent(IN), optional :: time_total

    integer, intent(IN) :: npx
    integer, intent(IN) :: npy
    integer, intent(IN) :: npz
    integer, intent(IN) :: nq             ! transported tracers
    integer, intent(IN) :: ng
    integer, intent(IN) :: ks
    integer, intent(IN) :: ncnst
    integer, intent(IN) :: n_split        ! small-step horizontal dynamics
    integer, intent(IN) :: q_split        ! tracer
    logical, intent(IN) :: fill
    logical, intent(IN) :: reproduce_sum
    logical, intent(IN) :: hydrostatic
    logical, intent(IN) :: hybrid_z       ! Using hybrid_z for remapping

    type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(inout), dimension(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) :: u ! D grid zonal wind (m/s)
    real, intent(inout), dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) :: v ! D grid meridional wind (m/s)
    real, intent(inout) :: w(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  !  W (m/s)
    real, intent(inout) :: pt(  bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! temperature (K)
    real, intent(inout) :: delp(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! pressure thickness (pascal)
    real, intent(inout) :: q(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz, ncnst) ! specific humidity and constituents
    real, intent(inout) :: delz(bd%isd:bd%ied,bd%jsd:bd%jed,npz)   ! delta-height (m); non-hydrostatic only
    real, intent(inout) ::  ze0(bd%is:bd%ie,bd%js:bd%je,npz+1) ! height at edges (m); non-hydrostatic

!-----------------------------------------------------------------------
! Auxilliary pressure arrays:    
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
    real, intent(inout) :: ps  (bd%isd:bd%ied  ,bd%jsd:bd%jed)           ! Surface pressure (pascal)
    real, intent(inout) :: pe  (bd%is-1:bd%ie+1, npz+1,bd%js-1:bd%je+1)  ! edge pressure (pascal)
    real, intent(inout) :: pk  (bd%is:bd%ie,bd%js:bd%je, npz+1)          ! pe**cappa
    real, intent(inout) :: peln(bd%is:bd%ie,npz+1,bd%js:bd%je)           ! ln(pe)
    real, intent(inout) :: pkz (bd%is:bd%ie,bd%js:bd%je,npz)             ! finite-volume mean pk
#ifdef PKC
    real, intent(inout) :: pkc (bd%isd:bd%ied,bd%jsd:bd%jed,npz+1)             ! finite-volume edge pk
#endif
    
!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
    real, intent(inout) :: phis(bd%isd:bd%ied,bd%jsd:bd%jed)       ! Surface geopotential (g*Z_surf)
    real, intent(inout) :: omga(bd%isd:bd%ied,bd%jsd:bd%jed,npz)   ! Vertical pressure velocity (pa/s)
    real, intent(inout) :: uc(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) ! (uc,vc) mostly used as the C grid winds
    real, intent(inout) :: vc(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz)

    real, intent(inout), dimension(bd%isd:bd%ied ,bd%jsd:bd%jed ,npz):: ua, va
    real, intent(in),    dimension(npz+1):: ak, bk

! Accumulated Mass flux arrays: the "Flux Capacitor"
    real, intent(inout) ::  mfx(bd%is:bd%ie+1, bd%js:bd%je,   npz)
    real, intent(inout) ::  mfy(bd%is:bd%ie  , bd%js:bd%je+1, npz)
! Accumulated Courant number arrays
    real, intent(inout) ::  cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    real, intent(inout) ::  cy(bd%isd:bd%ied ,bd%js:bd%je+1, npz)

    type(fv_grid_type),  intent(inout), target :: gridstruct
    type(fv_flags_type), intent(INOUT) :: flagstruct
    type(fv_nest_type),  intent(INOUT) :: neststruct
    type(domain2d), intent(INOUT) :: domain
    type(fv_atmos_type), intent(INOUT) :: parent_grid
    type(fv_diag_type), intent(IN) :: idiag

! Local Arrays
      real:: ws(bd%is:bd%ie,bd%js:bd%je)
      real:: te_2d(bd%is:bd%ie,bd%js:bd%je)
      real::   teq(bd%is:bd%ie,bd%js:bd%je)
      real:: pfull(npz)
      real:: gz(bd%is:bd%ie)
      real, allocatable :: dp1(:,:,:)
      real:: akap, rg, rdg, ph1, ph2, mdt
      integer:: kord_tracer(ncnst)
      integer :: i,j,k, n, iq, n_map
      integer :: sphum, liq_wat, ice_wat      ! GFDL physics
      integer :: rainwat, snowwat, graupel, cld_amt
      logical used, last_step
      integer, parameter :: max_packs=10
      integer:: i_pack(max_packs)

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

                                                  call timing_on('FV_DYN_START')
      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

      if ( flagstruct%no_dycore ) goto 911

      allocate ( dp1(is:ie, js:je, 1:npz) )
      dp1 = 0.
#ifdef SW_DYNAMICS
      akap  = 1.
      pfull(1) = 0.5*flagstruct%p_ref
#ifdef TEST_TRACER
      sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
#endif
#else


#if defined(MARS_GCM) || defined(VENUS_GCM) 
           sphum = 1
           cld_amt = -1   ! to cause trouble if (mis)used
#else
      if ( flagstruct%nwat==6 ) then
             sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
           liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
           ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
           rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat')
           snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat')
           graupel = get_tracer_index (MODEL_ATMOS, 'graupel')
           cld_amt = get_tracer_index (MODEL_ATMOS, 'cld_amt')
      elseif (flagstruct%nwat==3 ) then
             sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
           liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
           ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
           cld_amt = get_tracer_index (MODEL_ATMOS, 'cld_amt')
      else
           sphum = 1
           cld_amt = -1   ! to cause trouble if (mis)used
      endif
#endif

      akap  = kappa
      rg = kappa*cp_air

!$omp parallel do default(shared) private(ph1, ph2)
      do k=1,npz
         ph1 = ak(k  ) + bk(k  )*flagstruct%p_ref
         ph2 = ak(k+1) + bk(k+1)*flagstruct%p_ref
         pfull(k) = (ph2 - ph1) / log(ph2/ph1)
      enddo

      if ( flagstruct%fv_debug ) then
         call prt_maxmin('T_dyn_b',   pt, is, ie, js, je, ng, npz, 1.)
         call prt_maxmin('delp_b ', delp, is, ie, js, je, ng, npz, 0.01)
         call prt_maxmin('pk_b',    pk, is, ie, js, je, 0, npz+1, 1.)
         call prt_maxmin('pkz_b',  pkz, is, ie, js, je, 0, npz, 1.)
      endif

!$omp parallel do default(shared) 
      do k=1,npz
         do j=js,je
            do i=is,ie
               dp1(i,j,k) = zvir*q(i,j,k,sphum)
            enddo
         enddo
      enddo

! Re-compute pkz to enforce restart reproducibility when in solo mode
      if ( .not. hydrostatic ) then
         rdg = -rdgas / grav
!$omp parallel do default(shared)
         do k=1,npz
            do j=js,je
               do i=is,ie
                  pkz(i,j,k) = exp( kappa*log(rdg*delp(i,j,k)*pt(i,j,k)*    &
                       (1.+dp1(i,j,k))/delz(i,j,k)) )
               enddo
            enddo
         enddo
      endif

!---------------------
! Compute Total Energy
!---------------------
      if ( consv_te > 0. ) then
           call compute_total_energy(is, ie, js, je, isd, ied, jsd, jed, npz,        &
                                     u, v, w, delz, pt, delp, q, dp1, pe, peln, phis, &
                                     gridstruct%rsin2, gridstruct%cosa_s, &
                                     zvir, cp_air, rg, hlv, te_2d, ua, va, teq,            &
                                     flagstruct%moist_phys, sphum, liq_wat, ice_wat, hydrostatic, idiag%id_te)
           if( idiag%id_te>0 ) then
               used = send_data(idiag%id_te, teq, fv_time)
!              te_den=1.E-9*g_sum(teq, is, ie, js, je, ng, area, 0)/(grav*4.*pi*radius**2)
!              if(is_master())  write(*,*) 'Total Energy Density (Giga J/m**2)=',te_den
           endif
      endif

      if( flagstruct%tau > 0. ) then
        if ( flagstruct%n_sponge == 0 ) then
             call Rayleigh_Super(abs(bdt), npx, npy, npz, ks, pfull, flagstruct%tau, flagstruct%rf_center, u, v, w, pt,  &
                  ua, va, delz, gridstruct%agrid, cp_air, rg, ptop, hydrostatic, .true., flagstruct%rf_cutoff, gridstruct, domain, bd)
        else
             call Rayleigh_Friction(abs(bdt), npx, npy, npz, ks, pfull, flagstruct%tau, flagstruct%rf_center, u, v, w, pt,  &
                  ua, va, delz, cp_air, rg, ptop, hydrostatic, .true., flagstruct%rf_cutoff, gridstruct, domain, bd)
        endif
      endif

#endif

      !We call this BEFORE converting pt to virtual potential temperature, 
      !since we interpolate on (regular) temperature rather than theta.
      if (gridstruct%nested .or. ANY(neststruct%child_grids)) then
                                           call timing_on('NEST_BCs')
         call setup_nested_grid_BCs(npx, npy, npz, cp_air, zvir, ncnst, sphum,     &
              u, v, w, pt, delp, delz, q, uc, vc, pkz, &
              neststruct%nested, flagstruct%inline_q, flagstruct%make_nh, ng, &
              gridstruct, flagstruct, neststruct, &
              neststruct%nest_timestep, neststruct%tracer_nest_timestep, domain, bd)
                                           call timing_off('NEST_BCs')
      endif

#ifndef SW_DYNAMICS
! Convert pt to virtual potential temperature * CP
!$omp parallel do default(shared) 
  do k=1,npz
     do j=js,je
        do i=is,ie
           pt(i,j,k) = cp_air*pt(i,j,k)/pkz(i,j,k) * (1.+dp1(i,j,k))
        enddo
     enddo
  enddo
#endif

  last_step = .false.
  mdt = bdt / real(flagstruct%k_split)

                                                  call timing_off('FV_DYN_START')
                                                  call timing_on('FV_DYN_LOOP')
  do n_map=1, flagstruct%k_split   ! first level of time-split
                                           call timing_on('COMM_TOTAL')
      i_pack(1) = mpp_start_update_domains(delp, domain)
      i_pack(2) = mpp_start_update_domains(pt,   domain)
      i_pack(8) = mpp_start_update_domains(u, v, domain, gridtype=DGRID_NE)
                                           call timing_off('COMM_TOTAL')
!$omp parallel do default(shared) 
      do k=1,npz
         do j=js,je
            do i=is,ie
               dp1(i,j,k) = delp(i,j,k)
            enddo
         enddo
      enddo

      if ( n_map==flagstruct%k_split ) last_step = .true.

                                           call timing_on('DYN_CORE')
      call dyn_core(npx, npy, npz, ng, sphum, nq, mdt, n_split, zvir, cp_air, akap, grav, hydrostatic, &
                    u, v, w, delz, pt, q, delp, pe, pk, phis, ws, omga, ptop, pfull, ua, va,           & 
#ifdef PKC
                    uc, vc, mfx, mfy, cx, cy, pkz, pkc, peln, ak, bk, ks, &
#else
                    uc, vc, mfx, mfy, cx, cy, pkz, peln, ak, bk, ks, &
#endif
                    gridstruct, flagstruct, neststruct, idiag, bd, &
                    domain, n_map==1, i_pack, last_step, time_total)
                                           call timing_off('DYN_CORE')

  if ( flagstruct%fv_debug ) then
       call prt_maxmin('delp_a1',  delp, is, ie, js, je, ng, npz, 0.01)
       call prt_maxmin('T_dyn_a1',  pt, is, ie, js, je, ng, npz, 1.)
       call prt_maxmin('pk_a1',   pk, is, ie, js, je, 0, npz+1, 1.)
       call prt_maxmin('pkz_a1',  pkz, is, ie, js, je, 0, npz, 1.)
  endif

#ifdef SW_DYNAMICS
!$omp parallel do default(shared)
      do j=js,je
         do i=is,ie
            ps(i,j) = delp(i,j,1) / grav
         enddo
      enddo
#else
      if( .not. flagstruct%inline_q .and. nq /= 0 ) then    
!--------------------------------------------------------
! Perform large-time-step scalar transport using the accumulated CFL and
! mass fluxes
                                              call timing_on('tracer_2d')
       !!! CLEANUP: merge these two calls?
       if (gridstruct%nested .or. ANY(neststruct%child_grids)) then
         call tracer_2d_nested(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, npx, npy, npz, nq,    &
                        flagstruct%hord_tr, q_split, mdt, idiag%id_divg, i_pack(10), &
                        flagstruct%z_tracer, flagstruct%k_split, neststruct, parent_grid)          
       else
         call tracer_2d(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, npx, npy, npz, nq,    &
                        flagstruct%hord_tr, q_split, mdt, idiag%id_divg, i_pack(10), &
                        flagstruct%z_tracer, flagstruct%k_split)
       endif
                                             call timing_off('tracer_2d')

         if( last_step .and. idiag%id_divg>0 ) then
             used = send_data(idiag%id_divg, dp1, fv_time) 
             if(flagstruct%fv_debug) call prt_maxmin('divg',  dp1, is, ie, js, je, 0, npz, 1.)
         endif
      endif

      if ( npz > 4 ) then
!------------------------------------------------------------------------
! Peroform vertical remapping from Lagrangian control-volume to
! the Eulerian coordinate as specified by the routine set_eta.
! Note that this finite-volume dycore is otherwise independent of the vertical
! Eulerian coordinate.
!------------------------------------------------------------------------

         do iq=1,nq
                                kord_tracer(iq) = flagstruct%kord_tr
            if ( iq==cld_amt )  kord_tracer(iq) = 9      ! monotonic
         enddo
                                                  call timing_on('Remapping')

         call Lagrangian_to_Eulerian(last_step, consv_te, ps, pe, delp,          &
                     pkz, pk, mdt, bdt, npz, is,ie,js,je, isd,ied,jsd,jed,       &
                     nq, flagstruct%nwat, sphum, u,  v, w, delz, pt, q, phis,    &
                     zvir, cp_air, akap, flagstruct%kord_mt, flagstruct%kord_wz, &
                     kord_tracer, flagstruct%kord_tm, peln, te_2d,               &
                     ng, ua, va, omga, dp1, ws, fill, reproduce_sum,             &
                     ptop, ak, bk, gridstruct, domain, ze0, flagstruct%remap_t,  &
                     flagstruct%do_sat_adj, hydrostatic, hybrid_z, last_step)

                                                  call timing_off('Remapping')
         if( last_step )  then
!--------------------------
! Filter omega for physics:
!--------------------------
            if(flagstruct%nf_omega>0) call del2_cubed(omga, 0.20*gridstruct%da_min, gridstruct, domain, npx, npy, npz, flagstruct%nf_omega, bd)
            if(.not.hydrostatic .and. flagstruct%replace_w) then
!$omp parallel do default(shared) 
               do k=1,npz
                  do j=js,je
                     do i=is,ie
!                                          dz/dt = - omega / (g*density)
                        w(i,j,k) = delz(i,j,k)/delp(i,j,k) * omga(i,j,k)
                     enddo
                  enddo
               enddo
               if( is_master() ) write(6,*) 'Warning: W-wind replaced by scaled Omega'
               !This is the one place where we allow a member of flagstruct% to be 
               !changed outside of fv_control.
               flagstruct%replace_w = .false.
            endif
         endif
      end if

  enddo    ! n_map loop
                                                  call timing_off('FV_DYN_LOOP')

                                                  call timing_on('FV_DYN_END')
! Convert back to temperature
! if remap_t = .True. (which is so by default) and kord_tm < 0
! (which is not a default) then pt is converted in fv_mapz
#ifndef OLD_PT_TO_T
  if ( .not. (flagstruct%remap_t .and. flagstruct%kord_tm < 0) ) then
#endif
!$omp parallel do default(shared) 
     do k=1,npz
        do j=js,je
           do i=is,ie
              pt(i,j,k) = pt(i,j,k)*pkz(i,j,k)/(cp_air*(1.+zvir*q(i,j,k,sphum)))
           enddo
        enddo
     enddo
#ifndef OLD_PT_TO_T
  endif
#endif
#endif

      if ( flagstruct%fv_debug ) then
         call prt_maxmin('delp_a',  delp, is, ie, js, je, ng, npz, 0.01)
         call prt_maxmin('T_dyn_a',  pt, is, ie, js, je, ng, npz, 1.)
         call prt_maxmin('pk_a',   pk, is, ie, js, je, 0, npz+1, 1.)
         call prt_maxmin('pkz_a',  pkz, is, ie, js, je, 0, npz, 1.)
      endif

      deallocate ( dp1 )

  if( flagstruct%nwat==6 ) then
      call neg_adj3(is, ie, js, je, ng, npz,        &
                    flagstruct%hydrostatic,         &
                    peln, delz,                     &
                    pt, delp, q(isd,jsd,1,sphum),   &
                              q(isd,jsd,1,liq_wat), &
                              q(isd,jsd,1,rainwat), &
                              q(isd,jsd,1,ice_wat), &
                              q(isd,jsd,1,snowwat), &
                              q(isd,jsd,1,graupel), &
                              q(isd,jsd,1,cld_amt)  )
     if ( flagstruct%fv_debug ) then
       call prt_maxmin('SPHUM_dyn',   q(isd,jsd,1,sphum  ), is, ie, js, je, ng, npz, 1.)
       call prt_maxmin('liq_wat_dyn', q(isd,jsd,1,liq_wat), is, ie, js, je, ng, npz, 1.)
       call prt_maxmin('ice_wat_dyn', q(isd,jsd,1,ice_wat), is, ie, js, je, ng, npz, 1.)
       call prt_maxmin('snowwat_dyn', q(isd,jsd,1,snowwat), is, ie, js, je, ng, npz, 1.)
       call prt_maxmin('graupel_dyn', q(isd,jsd,1,graupel), is, ie, js, je, ng, npz, 1.)
     endif
  endif

911  call cubed_to_latlon(u, v, ua, va, gridstruct, &
          npx, npy, npz, 1, gridstruct%grid_type, domain, gridstruct%nested, flagstruct%c2l_ord, bd)
     if ( flagstruct%fv_debug ) then
       call prt_maxmin('UA', ua, is, ie, js, je, ng, npz, 1.)
       call prt_maxmin('VA', va, is, ie, js, je, ng, npz, 1.)
     endif

  if ( flagstruct%range_warn ) then
       call range_check('UA_dyn', ua, is, ie, js, je, ng, npz, gridstruct%agrid,   &
                         -220., 260., bad_range)
       call range_check('VA_dyn', ua, is, ie, js, je, ng, npz, gridstruct%agrid,   &
                         -220., 220., bad_range)
#ifndef SW_DYNAMICS
       call range_check('TA_dyn', pt, is, ie, js, je, ng, npz, gridstruct%agrid,   &
                         160., 330., bad_range)
       if ( .not. hydrostatic ) &
       call range_check('W_dyn', w, is, ie, js, je, ng, npz, gridstruct%agrid,   &
                         -20., 20., bad_range)
#endif
  endif

                                                  call timing_off('FV_DYN_END')

  end subroutine fv_dynamics


 subroutine Rayleigh_Super(dt, npx, npy, npz, ks, pm, tau, p_c, u, v, w, pt,  &
                           ua, va, delz, agrid, cp, rg, ptop, hydrostatic, conserve, rf_cutoff, gridstruct, domain, bd)
    real, intent(in):: dt
    real, intent(in):: tau              ! time scale (days)
    real, intent(in):: p_c
    real, intent(in):: cp, rg, ptop, rf_cutoff
    real, intent(in),  dimension(npz):: pm
    integer, intent(in):: npx, npy, npz, ks
    logical, intent(in):: hydrostatic
    logical, intent(in):: conserve
    type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(inout):: u(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) ! D grid zonal wind (m/s)
    real, intent(inout):: v(bd%isd:bd%ied+1,bd%jsd:bd%jed,npz) ! D grid meridional wind (m/s)
    real, intent(inout)::  w(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! cell center vertical wind (m/s)
    real, intent(inout):: pt(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! temp
    real, intent(inout):: ua(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! 
    real, intent(inout):: va(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! 
    real, intent(inout):: delz(bd%isd:bd%ied,bd%jsd:bd%jed,npz)   ! delta-height (m); non-hydrostatic only
    real,   intent(in) :: agrid(bd%isd:bd%ied,  bd%jsd:bd%jed,2)
    type(fv_grid_type), intent(IN) :: gridstruct
    type(domain2d), intent(INOUT) :: domain
!
    real, allocatable ::  u2f(:,:,:)
    real, parameter:: u0   = 60.   ! scaling velocity
    real, parameter:: sday = 86400.
    real c1, pc, rcv
    integer i, j, k

    integer :: is,  ie,  js,  je
    integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

    rcv = 1. / (cp - rg)

     if ( .not. RF_initialized ) then
          allocate( rf(npz) )

          if ( p_c <= 0. ) then
               pc = pm(1)
          else
               pc = p_c
          endif

          if( is_master() ) write(6,*) 'Rayleigh friction E-folding time [days]:'
          c1 = 1. / (tau*sday)
          do k=1, npz
             if ( pm(k) < rf_cutoff ) then
                  rf(k) = c1*(1.+tanh(log10(pc/pm(k))))
                  if( is_master() ) write(6,*) k, 0.01*pm(k), 1./(rf(k)*sday)
                  kmax = k
             else
                  exit
             endif
          enddo
          RF_initialized = .true.
     endif

    call c2l_ord2(u, v, ua, va, gridstruct, npz, gridstruct%grid_type, bd)

    allocate( u2f(isd:ied,jsd:jed,kmax) )
    u2f = 1.
!$omp parallel do default(shared)
    do k=1,kmax
       if ( pm(k) < rf_cutoff ) then
       do j=js,je
        if ( hydrostatic ) then
          do i=is,ie
             if ( abs(ua(i,j,k)) > 35.*cos(agrid(i,j,2)) )  then
                  u2f(i,j,k) = 1./(1.+dt*rf(k)*sqrt(ua(i,j,k)**2+va(i,j,k)**2)/u0)
             endif
          enddo
        else
          do i=is,ie
             if ( abs(ua(i,j,k)) > 35.*cos(agrid(i,j,2)) .or. abs(w(i,j,k))>7.5 )  then
                  u2f(i,j,k) = 1./(1.+dt*rf(k)*sqrt(ua(i,j,k)**2+va(i,j,k)**2+w(i,j,k)**2)/u0)
             endif
          enddo
        endif
       enddo
       endif ! p check
    enddo
                                        call timing_on('COMM_TOTAL')
    call mpp_update_domains(u2f, domain)
                                        call timing_off('COMM_TOTAL')


!$omp parallel do default(shared)
     do k=1,kmax
        if ( pm(k) < rf_cutoff ) then
! Add heat so as to conserve TE
          if ( conserve ) then
             if ( hydrostatic ) then
               do j=js,je
                  do i=is,ie
                     pt(i,j,k) = pt(i,j,k) + 0.5*(ua(i,j,k)**2+va(i,j,k)**2)*(1.-u2f(i,j,k)**2)/(cp-rg*ptop/pm(k))
                  enddo
               enddo
             else
               do j=js,je
                  do i=is,ie
                     pt(i,j,k) = pt(i,j,k) + 0.5*(ua(i,j,k)**2+va(i,j,k)**2+w(i,j,k)**2)*(1.-u2f(i,j,k)**2)*rcv
                  enddo
               enddo
             endif
          endif
             do j=js,je+1
                do i=is,ie
                   u(i,j,k) = 0.5*(u2f(i,j-1,k)+u2f(i,j,k))*u(i,j,k)
                enddo
             enddo
             do j=js,je
                do i=is,ie+1
                   v(i,j,k) = 0.5*(u2f(i-1,j,k)+u2f(i,j,k))*v(i,j,k)
                enddo
             enddo
          if ( .not. hydrostatic ) then
             do j=js,je
                do i=is,ie
                   w(i,j,k) = u2f(i,j,k)*w(i,j,k)
                enddo
             enddo
          endif
        endif
     enddo

     deallocate ( u2f )

 end subroutine Rayleigh_Super


 subroutine Rayleigh_Friction(dt, npx, npy, npz, ks, pm, tau, p_c, u, v, w, pt,  &
                              ua, va, delz, cp, rg, ptop, hydrostatic, conserve, &
                              rf_cutoff, gridstruct, domain, bd)
    real, intent(in):: dt
    real, intent(in):: tau              ! time scale (days)
    real, intent(in):: p_c
    real, intent(in):: cp, rg, ptop, rf_cutoff
    real, intent(in),  dimension(npz):: pm
    integer, intent(in):: npx, npy, npz, ks
    logical, intent(in):: hydrostatic
    logical, intent(in):: conserve
    type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(inout):: u(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) ! D grid zonal wind (m/s)
    real, intent(inout):: v(bd%isd:bd%ied+1,bd%jsd:bd%jed,npz) ! D grid meridional wind (m/s)
    real, intent(inout)::  w(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! cell center vertical wind (m/s)
    real, intent(inout):: pt(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! temp
    real, intent(inout):: ua(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! 
    real, intent(inout):: va(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! 
    real, intent(inout):: delz(bd%isd:bd%ied,bd%jsd:bd%jed,npz)   ! delta-height (m); non-hydrostatic only
    type(fv_grid_type), intent(IN) :: gridstruct
    type(domain2d), intent(INOUT) :: domain
! local:
    real, allocatable ::  u2f(:,:,:)
    real, parameter:: sday = 86400.
    real, parameter:: u000 = 4900.   ! scaling velocity  **2
    real c1, pc, rcv
    integer i, j, k

    integer :: is,  ie,  js,  je
    integer :: isd, ied, jsd, jed

    
    is  = bd%is
    ie  = bd%ie
    js  = bd%js
    je  = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed


    rcv = 1. / (cp - rg)

    if ( .not. RF_initialized ) then
          allocate( rf(npz) )
          allocate( rw(npz) )

          if ( p_c <= 0. ) then
               pc = pm(1)
          else
               pc = p_c
          endif

          if( is_master() ) write(6,*) 'Rayleigh friction E-folding time [days]:'
          c1 = 1. / (tau*sday)

          kmax = 1
          do k=1,npz
             if ( pm(k) < 40.E2 ) then
                  rf(k) = c1*(1.+tanh(log10(pc/pm(k))))
                  kmax = k
                  if( is_master() ) write(6,*) k, 0.01*pm(k), 1./(rf(k)*sday)
             else
                exit
             endif
          enddo
          if( is_master() ) write(6,*) 'Rayleigh Friction kmax=', kmax

          RF_initialized = .true.
    endif

    allocate( u2f(isd:ied,jsd:jed,kmax) )

    call c2l_ord2(u, v, ua, va, gridstruct, npz, gridstruct%grid_type, bd)
    u2f = 0.
!$omp parallel do default(shared)
    do k=1,kmax
        if ( hydrostatic ) then
           do j=js,je
              do i=is,ie
                 u2f(i,j,k) = ua(i,j,k)**2 + va(i,j,k)**2
              enddo
           enddo
        else
           do j=js,je
              do i=is,ie
                 u2f(i,j,k) = ua(i,j,k)**2 + va(i,j,k)**2 + w(i,j,k)**2
              enddo
           enddo
        endif
    enddo
                                        call timing_on('COMM_TOTAL')
    call mpp_update_domains(u2f, domain)
                                        call timing_off('COMM_TOTAL')

!$omp parallel do default(shared)
     do k=1,kmax

        if ( conserve ) then
           if ( hydrostatic ) then
             do j=js,je
                do i=is,ie
                   pt(i,j,k) = pt(i,j,k) + 0.5*u2f(i,j,k)/(cp-rg*ptop/pm(k))      &
                             * ( 1. - 1./(1.+dt*rf(k)*sqrt(u2f(i,j,k)/u000))**2 )
                enddo
             enddo
           else
             do j=js,je
                do i=is,ie
                   delz(i,j,k) = delz(i,j,k) / pt(i,j,k)
                   pt(i,j,k) = pt(i,j,k) + 0.5*u2f(i,j,k) * rcv      &
                             * ( 1. - 1./(1.+dt*rf(k)*sqrt(u2f(i,j,k)/u000))**2 )
                   delz(i,j,k) = delz(i,j,k) * pt(i,j,k)
                enddo
             enddo
           endif
        endif

        do j=js-1,je+1
           do i=is-1,ie+1
              u2f(i,j,k) = dt*rf(k)*sqrt(u2f(i,j,k)/u000)
           enddo
        enddo

        do j=js,je+1
           do i=is,ie
              u(i,j,k) = u(i,j,k) / (1.+0.5*(u2f(i,j-1,k)+u2f(i,j,k)))
           enddo
        enddo
        do j=js,je
           do i=is,ie+1
              v(i,j,k) = v(i,j,k) / (1.+0.5*(u2f(i-1,j,k)+u2f(i,j,k)))
           enddo
        enddo

        if ( .not. hydrostatic ) then
              do j=js,je
                 do i=is,ie
                    w(i,j,k) = w(i,j,k) / (1.+u2f(i,j,k))
                 enddo
              enddo
        endif

     enddo

     deallocate ( u2f )

 end subroutine Rayleigh_Friction


end module fv_dynamics_mod
