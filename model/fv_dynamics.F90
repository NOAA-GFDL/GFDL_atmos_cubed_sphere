module fv_dynamics_mod
   use constants_mod,       only: grav, pi, radius, hlv, rdgas, omega
   use dyn_core_mod,        only: dyn_core, del2_cubed
   use fv_mapz_mod,         only: compute_total_energy, Lagrangian_to_Eulerian
   use fv_tracer2d_mod,     only: tracer_2d, tracer_2d_1L, tracer_2d_nested
   use fv_grid_utils_mod,   only: cubed_to_latlon, c2l_ord2, g_sum
   use fv_mp_mod,           only: is_master
   use fv_timing_mod,       only: timing_on, timing_off
   use diag_manager_mod,    only: send_data
   use fv_diagnostics_mod,  only: fv_time, prt_mxm, range_check, prt_minmax
   use mpp_domains_mod,     only: DGRID_NE, CGRID_NE, mpp_update_domains, domain2D
   use mpp_domains_mod,     only: mpp_start_update_domains, mpp_complete_update_domains
   use field_manager_mod,   only: MODEL_ATMOS
   use tracer_manager_mod,  only: get_tracer_index
   use fv_sg_mod,           only: neg_adj3
   use fv_nesting_mod,      only: setup_nested_grid_BCs
   use fv_arrays_mod,       only: fv_grid_type, fv_flags_type, fv_atmos_type, fv_nest_type, fv_diag_type, fv_grid_bounds_type
   use fv_nwp_nudge_mod,    only: do_adiabatic_init

implicit none
   logical :: RF_initialized = .false.
   logical :: bad_range
   real, allocatable ::  rf(:)
   integer :: kmax=1
   real :: agrav
#ifdef HIWPP
   real, allocatable:: u00(:,:,:), v00(:,:,:)
#endif
private
public :: fv_dynamics

!---- version number -----
   character(len=128) :: version = '$Id$'
   character(len=128) :: tagname = '$Name$'

contains

!-----------------------------------------------------------------------
!     fv_dynamics :: FV dynamical core driver
!-----------------------------------------------------------------------
 
  subroutine fv_dynamics(npx, npy, npz, nq_tot,  ng, bdt, consv_te, fill,               &
                        reproduce_sum, kappa, cp_air, zvir, ptop, ks, ncnst, n_split,     &
                        q_split, u, v, w, delz, hydrostatic, pt, delp, q,   &
                        ps, pe, pk, peln, pkz, phis, q_con, omga, ua, va, uc, vc,          &
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
    integer, intent(IN) :: nq_tot             ! transported tracers
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
    real, intent(inout) :: w(   bd%isd:  ,bd%jsd:  ,1:)  !  W (m/s)
    real, intent(inout) :: pt(  bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! temperature (K)
    real, intent(inout) :: delp(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! pressure thickness (pascal)
    real, intent(inout) :: q(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz, ncnst) ! specific humidity and constituents
    real, intent(inout) :: delz(bd%isd:,bd%jsd:,1:)   ! delta-height (m); non-hydrostatic only
    real, intent(inout) ::  ze0(bd%is:, bd%js: ,1:) ! height at edges (m); non-hydrostatic

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
    real, intent(inout):: q_con(bd%isd:, bd%jsd:, 1:)
    
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
      real:: ps2(bd%isd:bd%ied,bd%jsd:bd%jed)
      real:: m_fac(bd%is:bd%ie,bd%js:bd%je)
      real:: pfull(npz)
      real:: gz(bd%is:bd%ie)
      real, allocatable :: dp1(:,:,:), dtdt_m(:,:,:)
      real:: akap, rg, rdg, ph1, ph2, mdt, gam, amdt, u0
      integer:: kord_tracer(ncnst)
      integer :: i,j,k, n, iq, n_map, nq, nwat, k_split
      integer :: sphum, liq_wat, ice_wat      ! GFDL physics
      integer :: rainwat, snowwat, graupel, cld_amt
      logical used, last_step, consv_am, do_omega
      integer, parameter :: max_packs=11
      integer:: i_pack(max_packs)
      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed
      real :: rcv, dt2, consv_fac
      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed
 
 
      agrav = 1. / grav
        dt2 = 0.5*bdt
      consv_am = flagstruct%consv_am
      k_split = flagstruct%k_split
      nwat = flagstruct%nwat
      nq = nq_tot - flagstruct%dnats
      allocate ( dp1(is:ie, js:je, 1:npz) )

      if ( flagstruct%no_dycore ) goto 911

#ifdef SW_DYNAMICS
      akap  = 1.
      pfull(1) = 0.5*flagstruct%p_ref
#ifdef TEST_TRACER
      sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
#endif
#else

      if (nwat>=3 ) then
             sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
           liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
           ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
           cld_amt = get_tracer_index (MODEL_ATMOS, 'cld_amt')
      endif
      if ( nwat==6 ) then
           rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat')
           snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat')
           graupel = get_tracer_index (MODEL_ATMOS, 'graupel')
      else
           sphum = 1
           cld_amt = -1   ! to cause trouble if (mis)used
      endif

      akap  = kappa
      rg = kappa*cp_air
      rcv = 1./(cp_air-rg)

!$omp parallel do default(shared) private(ph1, ph2)
      do k=1,npz
         ph1 = ak(k  ) + bk(k  )*flagstruct%p_ref
         ph2 = ak(k+1) + bk(k+1)*flagstruct%p_ref
         pfull(k) = (ph2 - ph1) / log(ph2/ph1)
      enddo

!$omp parallel do default(shared) 
      do k=1,npz
         do j=js,je
            do i=is,ie
               dp1(i,j,k) = zvir*q(i,j,k,sphum)
#ifdef USE_COND
#ifdef USE_NWAT3
               q_con(i,j,k) = q(i,j,k,liq_wat) + q(i,j,k,ice_wat)
#else
               q_con(i,j,k) = q(i,j,k,liq_wat) + q(i,j,k,rainwat) + q(i,j,k,ice_wat)  &
                            + q(i,j,k,snowwat) + q(i,j,k,graupel)
#endif
#endif
            enddo
         enddo
      enddo

      if ( .not. hydrostatic ) then
         rdg = -rdgas * agrav
!$omp parallel do default(shared)
         do k=1,npz
            do j=js,je
               do i=is,ie
#ifdef USE_COND
                  pkz(i,j,k) = exp( kappa*log(rdg*delp(i,j,k)*pt(i,j,k)*    &
                               (1.+dp1(i,j,k)-q_con(i,j,k))/delz(i,j,k)) )
#else
                  pkz(i,j,k) = exp( kappa*log(rdg*delp(i,j,k)*pt(i,j,k)*    &
                               (1.+dp1(i,j,k))/delz(i,j,k)) )
#endif
               enddo
            enddo
         enddo
      endif

      if ( flagstruct%fv_debug ) then
         call prt_mxm('PS',        ps, is, ie, js, je, ng,   1, 0.01, gridstruct%area, domain)
         call prt_mxm('T_dyn_b',   pt, is, ie, js, je, ng, npz, 1.,   gridstruct%area, domain)
         if (.not. hydrostatic) then
            call prt_mxm('delz',  delz, is, ie, js, je, ng, npz, 1., gridstruct%area, domain)
         endif
         call prt_mxm('delp_b ', delp, is, ie, js, je, ng, npz, 0.01, gridstruct%area, domain)
         call prt_mxm('pk_b',    pk, is, ie, js, je, 0, npz+1, 1.,gridstruct%area, domain)
         call prt_mxm('pkz_b',   pkz,is, ie, js, je, 0, npz,   1.,gridstruct%area, domain)
      endif

!---------------------
! Compute Total Energy
!---------------------
      if ( consv_te > 0.  .and. (.not.do_adiabatic_init) ) then
           call compute_total_energy(is, ie, js, je, isd, ied, jsd, jed, npz,        &
                                     u, v, w, delz, pt, delp, q, dp1, pe, peln, phis, &
                                     gridstruct%rsin2, gridstruct%cosa_s, &
                                     zvir, cp_air, rg, hlv, te_2d, ua, va, teq,            &
                                     flagstruct%moist_phys, sphum, liq_wat, rainwat,   &
                                     ice_wat, snowwat, graupel, hydrostatic, idiag%id_te)
           if( idiag%id_te>0 ) then
               used = send_data(idiag%id_te, teq, fv_time)
!              te_den=1.E-9*g_sum(teq, is, ie, js, je, ng, area, 0)/(grav*4.*pi*radius**2)
!              if(is_master())  write(*,*) 'Total Energy Density (Giga J/m**2)=',te_den
           endif
      endif

      if( (consv_am.or.idiag%id_amdt>0) .and. (.not.do_adiabatic_init) ) then
          call compute_aam(npz, is, ie, js, je, isd, ied, jsd, jed, gridstruct, bd,   &
                           ptop, ua, va, u, v, delp, teq, ps2, m_fac)
      endif

      if( flagstruct%tau > 0. ) then
        if ( gridstruct%grid_type<4 ) then
             call Rayleigh_Super(abs(bdt), npx, npy, npz, ks, pfull, phis, flagstruct%tau, u, v, w, pt,  &
                  ua, va, delz, gridstruct%agrid, cp_air, rg, ptop, hydrostatic, .true., flagstruct%rf_cutoff, gridstruct, domain, bd)
        else
             call Rayleigh_Friction(abs(bdt), npx, npy, npz, ks, pfull, flagstruct%tau, u, v, w, pt,  &
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
#ifdef USE_COND
           pt(i,j,k) = cp_air*pt(i,j,k)*(1.+dp1(i,j,k)-q_con(i,j,k))/pkz(i,j,k)
#else
           pt(i,j,k) = cp_air*pt(i,j,k)*(1.+dp1(i,j,k))/pkz(i,j,k)
#endif
        enddo
     enddo
  enddo
#endif

  last_step = .false.
  mdt = bdt / real(k_split)

  if ( idiag%id_mdt > 0 ) then
       allocate ( dtdt_m(is:ie,js:je,npz) )
!$omp parallel do default(shared) 
       do k=1,npz
          do j=js,je
             do i=is,ie
                dtdt_m(i,j,k) = 0.
             enddo
          enddo
       enddo
  endif
                                                  call timing_on('FV_DYN_LOOP')
  do n_map=1, k_split   ! first level of time-split
                                           call timing_on('COMM_TOTAL')
      i_pack(1) = mpp_start_update_domains(delp, domain)
      i_pack(2) = mpp_start_update_domains(pt,   domain)
#ifndef ROT3
      i_pack(8) = mpp_start_update_domains(u, v, domain, gridtype=DGRID_NE)
#endif
#ifdef USE_COND
      i_pack(11) = mpp_start_update_domains(q_con, domain)
#endif
                                           call timing_off('COMM_TOTAL')
!$omp parallel do default(shared) 
      do k=1,npz
         do j=js,je
            do i=is,ie
               dp1(i,j,k) = delp(i,j,k)
            enddo
         enddo
      enddo

      if ( n_map==k_split ) last_step = .true.
#ifdef USE_COND
                                           call timing_on('COMM_TOTAL')
     call mpp_complete_update_domains(i_pack(11), q_con, domain)
                                           call timing_off('COMM_TOTAL')
#endif

                                           call timing_on('DYN_CORE')
      call dyn_core(npx, npy, npz, ng, sphum, nq, mdt, n_split, zvir, cp_air, akap, grav, hydrostatic, &
                    u, v, w, delz, pt, q, delp, pe, pk, phis, ws, omga, ptop, pfull, ua, va,           & 
                    uc, vc, mfx, mfy, cx, cy, pkz, peln, q_con, ak, bk, ks, &
                    gridstruct, flagstruct, neststruct, idiag, bd, &
                    domain, n_map==1, i_pack, last_step, time_total)
                                           call timing_off('DYN_CORE')

  if ( flagstruct%fv_debug ) then
       call prt_mxm('delp_a1',  delp, is, ie, js, je, ng, npz, 0.01, gridstruct%area, domain)
       call prt_mxm('PT_dyn_a1',  pt, is, ie, js, je, ng, npz, 1., gridstruct%area, domain)
       call prt_mxm('pk_a1',   pk, is, ie, js, je, 0, npz+1, 1., gridstruct%area, domain)
  endif

#ifdef SW_DYNAMICS
!$omp parallel do default(shared)
      do j=js,je
         do i=is,ie
            ps(i,j) = delp(i,j,1) * agrav
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
                        flagstruct%z_tracer, k_split, neststruct, parent_grid)          
       else
         call tracer_2d(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, npx, npy, npz, nq,    &
                        flagstruct%hord_tr, q_split, mdt, idiag%id_divg, i_pack(10), &
                        flagstruct%z_tracer, k_split)
       endif
                                             call timing_off('tracer_2d')

         if( last_step .and. idiag%id_divg>0 ) then
             used = send_data(idiag%id_divg, dp1, fv_time) 
             if(flagstruct%fv_debug) call prt_mxm('divg',  dp1, is, ie, js, je, 0, npz, 1.,gridstruct%area, domain)
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
                     nq, nwat, sphum, q_con, u,  v, w, delz, pt, q, phis,    &
                     zvir, cp_air, akap, flagstruct%kord_mt, flagstruct%kord_wz, &
                     kord_tracer, flagstruct%kord_tm, peln, te_2d,               &
                     ng, ua, va, omga, dp1, ws, fill, reproduce_sum,             &
                     idiag%id_mdt>0, dtdt_m, &
                     ptop, ak, bk, gridstruct, domain, ze0, flagstruct%remap_t,  &
                     flagstruct%do_sat_adj, hydrostatic, hybrid_z, last_step, do_adiabatic_init)

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

! Convert back to temperature
            if ( .not. flagstruct%remap_t ) then
              if ( hydrostatic ) then
!$omp parallel do default(shared) 
               do k=1,npz
                  do j=js,je
                  do i=is,ie
                     pt(i,j,k) = pt(i,j,k)*pkz(i,j,k)/(cp_air*(1.+zvir*q(i,j,k,sphum)))
                  enddo
                  enddo
               enddo
              else
!$omp parallel do default(shared) 
               do k=1,npz
                  do j=js,je
                  do i=is,ie
#ifdef USE_COND
                     pt(i,j,k) = pt(i,j,k)*pkz(i,j,k)/(cp_air*(1.+zvir*q(i,j,k,sphum)-q_con(i,j,k)))
#else
                     pt(i,j,k) = pt(i,j,k)*pkz(i,j,k)/(cp_air*(1.+zvir*q(i,j,k,sphum)))
#endif
                  enddo
                  enddo
               enddo
            endif
            endif

         endif
      end if
#endif
  enddo    ! n_map loop
                                                  call timing_off('FV_DYN_LOOP')

  if ( idiag%id_mdt > 0 .and. (.not.do_adiabatic_init) ) then
! Output temperature tendency due to inline moist physics:
!$omp parallel do default(shared) 
       do k=1,npz
          do j=js,je
             do i=is,ie
                dtdt_m(i,j,k) = dtdt_m(i,j,k) / bdt
             enddo
          enddo
       enddo
       used = send_data(idiag%id_mdt, dtdt_m, fv_time)
       deallocate ( dtdt_m )
  endif


  if( nwat==6 ) then
     if (cld_amt > 0) then
      call neg_adj3(is, ie, js, je, ng, npz,        &
                    flagstruct%hydrostatic,         &
                    peln, delz,                     &
                    pt, delp, q(isd,jsd,1,sphum),   &
                              q(isd,jsd,1,liq_wat), &
                              q(isd,jsd,1,rainwat), &
                              q(isd,jsd,1,ice_wat), &
                              q(isd,jsd,1,snowwat), &
                              q(isd,jsd,1,graupel), &
                              q(isd,jsd,1,cld_amt), flagstruct%check_negative)
     else
        call neg_adj3(is, ie, js, je, ng, npz,        &
                      flagstruct%hydrostatic,         &
                      peln, delz,                     &
                      pt, delp, q(isd,jsd,1,sphum),   &
                                q(isd,jsd,1,liq_wat), &
                                q(isd,jsd,1,rainwat), &
                                q(isd,jsd,1,ice_wat), &
                                q(isd,jsd,1,snowwat), &
                                q(isd,jsd,1,graupel), check_negative=flagstruct%check_negative)
     endif
     if ( flagstruct%fv_debug ) then
       call prt_mxm('T_dyn_a3',    pt, is, ie, js, je, ng, npz, 1., gridstruct%area, domain)
       call prt_mxm('SPHUM_dyn',   q(isd,jsd,1,sphum  ), is, ie, js, je, ng, npz, 1.,gridstruct%area, domain)
       call prt_mxm('liq_wat_dyn', q(isd,jsd,1,liq_wat), is, ie, js, je, ng, npz, 1.,gridstruct%area, domain)
       call prt_mxm('rainwat_dyn', q(isd,jsd,1,rainwat), is, ie, js, je, ng, npz, 1.,gridstruct%area, domain)
       call prt_mxm('ice_wat_dyn', q(isd,jsd,1,ice_wat), is, ie, js, je, ng, npz, 1.,gridstruct%area, domain)
       call prt_mxm('snowwat_dyn', q(isd,jsd,1,snowwat), is, ie, js, je, ng, npz, 1.,gridstruct%area, domain)
       call prt_mxm('graupel_dyn', q(isd,jsd,1,graupel), is, ie, js, je, ng, npz, 1.,gridstruct%area, domain)
     endif
  endif

  if( consv_am .or. idiag%id_amdt>0 .or. idiag%id_aam>0 .and. (.not.do_adiabatic_init)  ) then
      call compute_aam(npz, is, ie, js, je, isd, ied, jsd, jed, gridstruct, bd,   &
                       ptop, ua, va, u, v, delp, te_2d, ps, m_fac)
      if( idiag%id_aam>0 ) then
          used = send_data(idiag%id_aam, te_2d, fv_time)
          if ( prt_minmax ) then
             gam = g_sum( domain, te_2d, is, ie, js, je, ng, gridstruct%area, 0) 
             if( is_master() ) write(6,*) 'Total AAM =', gam
          endif
      endif
  endif

  if( consv_am .or. idiag%id_amdt>0 .and. (.not.do_adiabatic_init)  ) then
!$omp parallel do default(shared) 
      do j=js,je
         do i=is,ie
! Note: the mountain torque computation contains also numerical error
! The numerical error is mostly from the zonal gradient of the terrain (zxg)
            te_2d(i,j) = te_2d(i,j)-teq(i,j) + dt2*(ps2(i,j)+ps(i,j))*idiag%zxg(i,j)
         enddo
      enddo
      if( idiag%id_amdt>0 ) used = send_data(idiag%id_amdt, te_2d/bdt, fv_time)

      if ( consv_am .or. prt_minmax ) then
         amdt = g_sum( domain, te_2d, is, ie, js, je, ng, gridstruct%area, 0) 
         u0 = -radius*amdt/g_sum( domain, m_fac, is, ie, js, je, ng, gridstruct%area, 0)
         u0 = real(u0, 4)    ! truncate to enforce reproducibility
         if(is_master() .and. prt_minmax)         &
!        write(6,*) 'Dynamic AM tendency =', amdt/(bdt*1.e18), 'del-u (per yr)=', u0*365.*86400./bdt
         write(6,*) 'Dynamic Angular Momentum tendency (Hadleys)=', amdt/(bdt*1.e18)
      endif

      if( consv_am ) then
!$omp parallel do default(shared)
      do j=js,je
         do i=is,ie
            m_fac(i,j) = u0*cos(gridstruct%agrid(i,j,2))
         enddo
      enddo
!$omp parallel do default(shared)
      do k=1,npz
      if ( hydrostatic ) then
         do j=js,je
         do i=is,ie
            pt(i,j,k) = pt(i,j,k) - m_fac(i,j)*(0.5*m_fac(i,j)+ua(i,j,k))/cp_air
         enddo
         enddo
      else
         do j=js,je
         do i=is,ie
            pt(i,j,k) = pt(i,j,k) - m_fac(i,j)*(0.5*m_fac(i,j)+ua(i,j,k))*rcv
         enddo
         enddo
      endif
      do j=js,je+1
         do i=is,ie
            u(i,j,k) = u(i,j,k) + u0*gridstruct%l2c_u(i,j)
         enddo
      enddo
      do j=js,je
         do i=is,ie+1
            v(i,j,k) = v(i,j,k) + u0*gridstruct%l2c_v(i,j)
         enddo
      enddo
      enddo
      endif   !  consv_am
  endif

911  call cubed_to_latlon(u, v, ua, va, gridstruct, &
          npx, npy, npz, 1, gridstruct%grid_type, domain, gridstruct%nested, flagstruct%c2l_ord, bd)

     if ( flagstruct%fv_debug ) then
       call prt_mxm('UA', ua, is, ie, js, je, ng, npz, 1., gridstruct%area, domain)
       call prt_mxm('VA', va, is, ie, js, je, ng, npz, 1., gridstruct%area, domain)
       call prt_mxm('TA', pt, is, ie, js, je, ng, npz, 1., gridstruct%area, domain)
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
  deallocate ( dp1 )

  end subroutine fv_dynamics


 subroutine Rayleigh_Super(dt, npx, npy, npz, ks, pm, phis, tau, u, v, w, pt,  &
                           ua, va, delz, agrid, cp, rg, ptop, hydrostatic, conserve, rf_cutoff, gridstruct, domain, bd)
    real, intent(in):: dt
    real, intent(in):: tau              ! time scale (days)
    real, intent(in):: cp, rg, ptop, rf_cutoff
    real, intent(in),  dimension(npz):: pm
    integer, intent(in):: npx, npy, npz, ks
    logical, intent(in):: hydrostatic
    logical, intent(in):: conserve
    type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(inout):: u(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) ! D grid zonal wind (m/s)
    real, intent(inout):: v(bd%isd:bd%ied+1,bd%jsd:bd%jed,npz) ! D grid meridional wind (m/s)
    real, intent(inout)::  w(bd%isd:      ,bd%jsd:      ,1: ) ! cell center vertical wind (m/s)
    real, intent(inout):: pt(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! temp
    real, intent(inout):: ua(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! 
    real, intent(inout):: va(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! 
    real, intent(inout):: delz(bd%isd:      ,bd%jsd:      ,1: )   ! delta-height (m); non-hydrostatic only
    real,   intent(in) :: agrid(bd%isd:bd%ied,  bd%jsd:bd%jed,2)
    real, intent(in) :: phis(bd%isd:bd%ied,bd%jsd:bd%jed)       ! Surface geopotential (g*Z_surf)
    type(fv_grid_type), intent(IN) :: gridstruct
    type(domain2d), intent(INOUT) :: domain
!
    real, allocatable ::  u2f(:,:,:)
    real, parameter:: u0   = 60.   ! scaling velocity
    real, parameter:: sday = 86400.
    real rcv, tau0
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
#ifdef HIWPP
          allocate ( u00(is:ie,  js:je+1,npz) )
          allocate ( v00(is:ie+1,js:je  ,npz) )
!$omp parallel do default(shared)
          do k=1,npz
             do j=js,je+1
                do i=is,ie
                   u00(i,j,k) = u(i,j,k)
                enddo
             enddo
             do j=js,je
                do i=is,ie+1
                   v00(i,j,k) = v(i,j,k)
                enddo
             enddo
          enddo
#endif
#ifdef SMALL_EARTH
          tau0 = tau
#else
          tau0 = tau * sday
#endif
          allocate( rf(npz) )
          rf(:) = 0.

          if( is_master() ) write(6,*) 'Rayleigh friction E-folding time (days):'
          do k=1, npz
             if ( pm(k) < rf_cutoff ) then
                  rf(k) = dt/tau0*sin(0.5*pi*log(rf_cutoff/pm(k))/log(rf_cutoff/ptop))**2
                  if( is_master() ) write(6,*) k, 0.01*pm(k), dt/(rf(k)*sday)
                  kmax = k
             else
                  exit
             endif
          enddo
          RF_initialized = .true.
     endif

    call c2l_ord2(u, v, ua, va, gridstruct, npz, gridstruct%grid_type, bd)

    allocate( u2f(isd:ied,jsd:jed,kmax) )

!$omp parallel do default(shared)
    do k=1,kmax
       if ( pm(k) < rf_cutoff ) then
       do j=js,je
        if ( hydrostatic ) then
          do i=is,ie
             if ( abs(ua(i,j,k)) > 35.*cos(agrid(i,j,2)) )  then
                  u2f(i,j,k) = 1./(1.+rf(k)*sqrt(ua(i,j,k)**2+va(i,j,k)**2)/u0)
             else
                  u2f(i,j,k) = 1.
             endif
          enddo
        else
          do i=is,ie
             if ( abs(ua(i,j,k)) > 35.*cos(agrid(i,j,2)) .or. abs(w(i,j,k))>7.5 )  then
                  u2f(i,j,k) = 1./(1.+rf(k)*sqrt(ua(i,j,k)**2+va(i,j,k)**2+w(i,j,k)**2)/u0)
             else
                  u2f(i,j,k) = 1.
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
#ifdef HIWPP
             do j=js,je
                do i=is,ie
                   w(i,j,k) = w(i,j,k)/(1.+rf(k))
                enddo
             enddo
             do j=js,je+1
                do i=is,ie
                   u(i,j,k) = (u(i,j,k)+rf(k)*u00(i,j,k))/(1.+rf(k))
                enddo
             enddo
             do j=js,je
                do i=is,ie+1
                   v(i,j,k) = (v(i,j,k)+rf(k)*v00(i,j,k))/(1.+rf(k))
                enddo
             enddo
#else
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
#endif
        endif
     enddo

     deallocate ( u2f )

 end subroutine Rayleigh_Super


 subroutine Rayleigh_Friction(dt, npx, npy, npz, ks, pm, tau, u, v, w, pt,  &
                              ua, va, delz, cp, rg, ptop, hydrostatic, conserve, &
                              rf_cutoff, gridstruct, domain, bd)
    real, intent(in):: dt
    real, intent(in):: tau              ! time scale (days)
    real, intent(in):: cp, rg, ptop, rf_cutoff
    real, intent(in),  dimension(npz):: pm
    integer, intent(in):: npx, npy, npz, ks
    logical, intent(in):: hydrostatic
    logical, intent(in):: conserve
    type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(inout):: u(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) ! D grid zonal wind (m/s)
    real, intent(inout):: v(bd%isd:bd%ied+1,bd%jsd:bd%jed,npz) ! D grid meridional wind (m/s)
    real, intent(inout)::  w(bd%isd:      ,bd%jsd:      ,1: ) ! cell center vertical wind (m/s)
    real, intent(inout):: pt(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! temp
    real, intent(inout):: ua(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! 
    real, intent(inout):: va(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! 
    real, intent(inout):: delz(bd%isd:      ,bd%jsd:      ,1: )   ! delta-height (m); non-hydrostatic only
    type(fv_grid_type), intent(IN) :: gridstruct
    type(domain2d), intent(INOUT) :: domain
! local:
    real, allocatable ::  u2f(:,:,:)
    real, parameter:: sday = 86400.
    real, parameter:: u000 = 4900.   ! scaling velocity  **2
    real  rcv
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
          if( is_master() ) write(6,*) 'Rayleigh friction E-folding time (days):'
          do k=1, npz
             if ( pm(k) < rf_cutoff ) then
                  rf(k) = dt/(tau*sday)*sin(0.5*pi*log(rf_cutoff/pm(k))/log(rf_cutoff/ptop))**2
                  if( is_master() ) write(6,*) k, 0.01*pm(k), dt/(rf(k)*sday)
                  kmax = k
             else
                  exit
             endif
          enddo
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
                             * ( 1. - 1./(1.+rf(k)*sqrt(u2f(i,j,k)/u000))**2 )
                enddo
             enddo
           else
             do j=js,je
                do i=is,ie
                   delz(i,j,k) = delz(i,j,k) / pt(i,j,k)
                   pt(i,j,k) = pt(i,j,k) + 0.5*u2f(i,j,k) * rcv      &
                             * ( 1. - 1./(1.+rf(k)*sqrt(u2f(i,j,k)/u000))**2 )
                   delz(i,j,k) = delz(i,j,k) * pt(i,j,k)
                enddo
             enddo
           endif
        endif

        do j=js-1,je+1
           do i=is-1,ie+1
              u2f(i,j,k) = rf(k)*sqrt(u2f(i,j,k)/u000)
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

 subroutine compute_aam(npz, is, ie, js, je, isd, ied, jsd, jed, gridstruct, bd,   &
                        ptop, ua, va, u, v, delp, aam, ps, m_fac)
! Compute vertically (mass) integrated Atmospheric Angular Momentum
    integer, intent(in):: npz
    integer, intent(in):: is,  ie,  js,  je
    integer, intent(in):: isd, ied, jsd, jed
    real, intent(in):: ptop
    real, intent(inout):: u(isd:ied  ,jsd:jed+1,npz) ! D grid zonal wind (m/s)
    real, intent(inout):: v(isd:ied+1,jsd:jed,npz) ! D grid meridional wind (m/s)
    real, intent(inout):: delp(isd:ied,jsd:jed,npz)
    real, intent(inout), dimension(isd:ied,jsd:jed, npz):: ua, va
    real, intent(out):: aam(is:ie,js:je)
    real, intent(out):: m_fac(is:ie,js:je)
    real, intent(out):: ps(isd:ied,jsd:jed)
    type(fv_grid_bounds_type), intent(IN) :: bd
    type(fv_grid_type), intent(IN) :: gridstruct
! local:
    real, dimension(is:ie):: r1, r2, dm
    integer i, j, k

  call c2l_ord2(u, v, ua, va, gridstruct, npz, gridstruct%grid_type, bd)
    
!$omp parallel do default(shared) private(r1, r2, dm)
  do j=js,je
     do i=is,ie
        r1(i) = radius*cos(gridstruct%agrid(i,j,2))
        r2(i) = r1(i)*r1(i)
        aam(i,j) = 0.
        m_fac(i,j) = 0.
        ps(i,j) = ptop
     enddo
     do k=1,npz
        do i=is,ie
           dm(i) = delp(i,j,k)
           ps(i,j) = ps(i,j) + dm(i)
           dm(i) = dm(i)*agrav
           aam(i,j) = aam(i,j) + (r2(i)*omega + r1(i)*ua(i,j,k)) * dm(i)
           m_fac(i,j) = m_fac(i,j) + dm(i)*r2(i)
        enddo
     enddo
  enddo

 end subroutine compute_aam

end module fv_dynamics_mod
