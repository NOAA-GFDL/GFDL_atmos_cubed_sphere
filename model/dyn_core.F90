module dyn_core_mod

  use constants_mod,      only: rdgas
  use mpp_mod,            only: mpp_pe 
  use mpp_domains_mod,    only: CGRID_NE, DGRID_NE, mpp_get_boundary, mpp_update_domains,  &
                                mpp_start_update_domains, mpp_complete_update_domains, domain2d
  use mpp_parameter_mod,  only: CORNER
  use fv_mp_mod,          only: is_master
  use sw_core_mod,        only: c_sw, d_sw
  use a2b_edge_mod,       only: a2b_ord2, a2b_ord4
  use nh_core_mod,        only: Riem_Solver3, Riem_Solver_C, update_dz_c, update_dz_d, geopk_halo_nh
  use tp_core_mod,        only: copy_corners
  use fv_timing_mod,      only: timing_on, timing_off
  use fv_diagnostics_mod, only: prt_maxmin, fv_time
#if defined (ADA_NUDGE)
  use fv_ada_nudge_mod,   only: breed_slp_inline_ada
#else
  use fv_nwp_nudge_mod,   only: breed_slp_inline
#endif
  use diag_manager_mod,   only: send_data
  use fv_arrays_mod,      only: fv_grid_type, fv_flags_type, fv_nest_type, fv_diag_type, fv_grid_bounds_type

  use boundary_mod,         only: extrapolation_BC,  nested_grid_BC_apply, nested_grid_BC_apply_intT

#ifdef SW_DYNAMICS
  use test_cases_mod,      only: test_case, case9_forcing1, case9_forcing2
#endif

implicit none
private

public :: dyn_core, del2_cubed

  real :: ptk, rgrav
  real :: d3_damp
  real, allocatable, dimension(:,:,:) ::  ut, vt, crx, cry, xfx, yfx, divgd, &
#ifdef PKC
                                          zh, du, dv, delpc, pk3, ptc, gz
#else
                                          zh, du, dv, pkc, delpc, pk3, ptc, gz
#endif

!---- version number -----
  character(len=128) :: version = '$Id: dyn_core.F90,v 20.0 2013/12/13 23:04:18 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal_201409 $'

contains

!-----------------------------------------------------------------------
!     dyn_core :: FV Lagrangian dynamics driver
!-----------------------------------------------------------------------
 
 subroutine dyn_core(npx, npy, npz, ng, sphum, nq, bdt, n_split, zvir, cp, akap, grav, hydrostatic,  &
                     u,  v,  w, delz, pt, q, delp, pe, pk, phis, ws, omga, ptop, pfull, ua, va, & 
#ifdef PKC
                     uc, vc, mfx, mfy, cx, cy, pkz, pkc, peln, ak, bk, &
#else
                     uc, vc, mfx, mfy, cx, cy, pkz, peln, ak, bk, &
#endif
                     ks, gridstruct, flagstruct, neststruct, idiag, bd, domain, &
                     init_step, i_pack, end_step, time_total)
    integer, intent(IN) :: npx
    integer, intent(IN) :: npy
    integer, intent(IN) :: npz
    integer, intent(IN) :: ng, nq, sphum
    integer, intent(IN) :: n_split
    real   , intent(IN) :: bdt
    real   , intent(IN) :: zvir, cp, akap, grav
    real   , intent(IN) :: ptop
    logical, intent(IN) :: hydrostatic
    logical, intent(IN) :: init_step, end_step
    real, intent(in) :: pfull(npz)
    real, intent(in),     dimension(npz+1) :: ak, bk
    integer, intent(IN) :: ks
    integer, intent(inout) :: i_pack(*)
    type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(inout), dimension(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz):: u  ! D grid zonal wind (m/s)
    real, intent(inout), dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz):: v  ! D grid meridional wind (m/s)
    real, intent(inout) :: w(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! vertical vel. (m/s)
    real, intent(inout) :: delz(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! delta-height (m)
    real, intent(inout) :: pt(  bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! temperature (K)
    real, intent(inout) :: delp(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! pressure thickness (pascal)
    real, intent(inout) :: q(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz, nq)  ! 
    real, intent(in), optional:: time_total  ! total time (seconds) since start

!-----------------------------------------------------------------------
! Auxilliary pressure arrays:    
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
    real, intent(inout):: phis(bd%isd:bd%ied,bd%jsd:bd%jed)      ! Surface geopotential (g*Z_surf)
    real, intent(inout):: pe(bd%is-1:bd%ie+1, npz+1,bd%js-1:bd%je+1)  ! edge pressure (pascal)
    real, intent(inout):: peln(bd%is:bd%ie,npz+1,bd%js:bd%je)          ! ln(pe)
    real, intent(inout):: pk(bd%is:bd%ie,bd%js:bd%je, npz+1)        ! pe**kappa

!-----------------------------------------------------------------------
! Others:
    real,    parameter:: near0 = 1.E-8
    real,    parameter:: huge_r = 1.E40
    real,    parameter:: air_viscosity = 1.E-5   ! [m**2/sec] for T ~ 260 K
!-----------------------------------------------------------------------
    real, intent(out  ):: ws(bd%is:bd%ie,bd%js:bd%je)        ! w at surface
    real, intent(inout):: omga(bd%isd:bd%ied,bd%jsd:bd%jed,npz)    ! Vertical pressure velocity (pa/s)
    real, intent(inout):: uc(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz)  ! (uc, vc) are mostly used as the C grid winds
    real, intent(inout):: vc(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz)
    real, intent(inout), dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz):: ua, va

! The Flux capacitors: accumulated Mass flux arrays
    real, intent(inout)::  mfx(bd%is:bd%ie+1, bd%js:bd%je,   npz)
    real, intent(inout)::  mfy(bd%is:bd%ie  , bd%js:bd%je+1, npz)
! Accumulated Courant number arrays
    real, intent(inout)::  cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    real, intent(inout)::  cy(bd%isd:bd%ied ,bd%js:bd%je+1, npz)
! Work:
    real, intent(inout):: pkz(bd%is:bd%ie,bd%js:bd%je,npz)  ! 
#ifdef PKC
    real, intent(inout):: pkc(bd%isd:bd%ied,bd%jsd:bd%jed,npz+1)  ! 
#endif

    type(fv_grid_type),  intent(INOUT), target :: gridstruct
    type(fv_flags_type), intent(IN),    target :: flagstruct
    type(fv_nest_type),  intent(INOUT)         :: neststruct
    type(fv_diag_type),  intent(IN)            :: idiag
    type(domain2d),      intent(INOUT)         :: domain

    real, allocatable, dimension(:,:,:):: pem, heat_source
! Auto 1D & 2D arrays:
    real, dimension(bd%isd:bd%ied,bd%jsd:bd%jed):: ws3
    real:: dp_ref(npz)
    real:: zs(bd%isd:bd%ied,bd%jsd:bd%jed)        ! surface height (m)
    real:: p1d(bd%is:bd%ie)
    real:: om2d(bd%is:bd%ie,npz)
    real wbuffer(npy+2,npz)
    real ebuffer(npy+2,npz)
    real nbuffer(npx+2,npz)
    real sbuffer(npx+2,npz)
! ----   For external mode:
    real divg2(bd%is:bd%ie+1,bd%js:bd%je+1)
    real wk(bd%isd:bd%ied,bd%jsd:bd%jed)
    real fz(bd%is: bd%ie+1,bd%js: bd%je+1)
    real heat_s(bd%is:bd%ie,bd%js:bd%je)
    real damp_vt(npz+1)
    integer nord_v(npz+1)
!-------------------------------------
    integer :: hord_m, hord_v, hord_t, hord_p, nord_k
    integer :: w_pack
    integer :: ms
!---------------------------------------
    integer :: i,j,k, it, iq, n_con
    integer :: iep1, jep1
    real    :: beta_d, damp_k, d_con_k
    real    :: dt, dt2, rdt
    real    :: d2_divg, dd_divg, d3_divg
    real    :: diff_z0, k1k, kapag, tmcp
    logical :: last_step, remap_step
    logical used
    real :: split_timestep_bc

    real, pointer, dimension(:,:) :: rarea

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

    ptk = ptop ** akap
    dt  = bdt / real(n_split)
    dt2 = 0.5*dt
    rdt = 1.0/dt
    ms = max(1, flagstruct%m_split/2)

! Indexes:
    iep1 = ie + 1
    jep1 = je + 1

!!! Initialize pointers
    rarea => gridstruct%rarea

    if ( .not.hydrostatic ) then

         rgrav = 1.0/grav
           k1k =  akap / (1.-akap)    ! rg/Cv=0.4
         kapag = -akap / grav
       diff_z0 = 0.25 * flagstruct%scale_z**2

!$omp parallel do default(shared)
       do k=1,npz
          dp_ref(k) = ak(k+1)-ak(k) + (bk(k+1)-bk(k))*1.E5  
       enddo

!$omp parallel do default(shared)
       do j=jsd,jed
          do i=isd,ied
             zs(i,j) = phis(i,j) * rgrav
          enddo
       enddo
    endif

      if ( init_step ) then  ! Start of the big dynamic time stepping

           allocate(    gz(isd:ied, jsd:jed ,npz+1) )
             call init_ijk_mem(isd,ied, jsd,jed, npz+1, gz, huge_r)
#ifndef PKC
           allocate(   pkc(isd:ied, jsd:jed ,npz+1) )
#endif
           allocate(   ptc(isd:ied, jsd:jed ,npz ) )
           allocate( crx(is :ie+1, jsd:jed,  npz) )
           allocate( xfx(is :ie+1, jsd:jed,  npz) )
           allocate( cry(isd:ied,  js :je+1, npz) )
           allocate( yfx(isd:ied,  js :je+1, npz) )
           allocate( divgd(isd:ied+1,jsd:jed+1,npz) )
           allocate( delpc(isd:ied, jsd:jed  ,npz  ) )
                     call init_ijk_mem(isd,ied, jsd,jed, npz, delpc, 0.)

           allocate( ut(isd:ied, jsd:jed, npz) )
                     call init_ijk_mem(isd,ied, jsd,jed, npz, ut, 0.)
           allocate( vt(isd:ied, jsd:jed, npz) )
                     call init_ijk_mem(isd,ied, jsd,jed, npz, vt, 0.)

          if ( .not. hydrostatic ) then
               allocate( zh(isd:ied, jsd:jed, npz+1) )
               call init_ijk_mem(isd,ied, jsd,jed, npz+1, zh, huge_r )
               allocate ( pk3(isd:ied,jsd:jed,npz+1) )
          endif
          if ( flagstruct%beta > near0 ) then
               allocate( du(isd:ied,  jsd:jed+1,npz) )
               call init_ijk_mem(isd,ied,   jsd,jed+1, npz, du, 0.)
               allocate( dv(isd:ied+1,jsd:jed,  npz) )
               call init_ijk_mem(isd,ied+1, jsd,jed  , npz, dv, 0.)
          endif
      endif    ! end init_step

      !moved to fv_control
      !if ( flagstruct%beta < 1.e-5 ) flagstruct%beta = 0.

! Empty the "flux capacitors"
    call init_ijk_mem(is, ie+1, js,  je,   npz, mfx, 0.)
    call init_ijk_mem(is, ie  , js,  je+1, npz, mfy, 0.)
    call init_ijk_mem(is, ie+1, jsd, jed,  npz, cx, 0.)
    call init_ijk_mem(isd, ied, js,  je+1, npz, cy, 0.)

    if ( flagstruct%d_con > 1.0E-5 ) then
         allocate( heat_source(isd:ied, jsd:jed, npz) )
         call init_ijk_mem(isd, ied, jsd, jed, npz, heat_source, 0.)
    endif

    if ( flagstruct%convert_ke .or. flagstruct%vtdm4> 1.E-3 ) then
         n_con = npz
    else
         if ( flagstruct%d2_bg_k1 < 1.E-3 ) then
              n_con = 0
         else
              if ( flagstruct%d2_bg_k2 < 1.E-3 ) then
                   n_con = 1
              else
                   n_con = 2
              endif
         endif
    endif


     if (gridstruct%nested) then
        split_timestep_bc = real(n_split*flagstruct%k_split+neststruct%nest_timestep)
     endif

!-----------------------------------------------------
  do it=1,n_split
!-----------------------------------------------------
     if ( flagstruct%breed_vortex_inline .or. it==n_split ) then
          remap_step = .true.
     else
          remap_step = .false.
     endif

     if ( flagstruct%fv_debug ) then
          if(is_master()) write(*,*) 'n_split loop, it=', it
     endif

     if ( nq > 0 ) then
                                    call timing_on('COMM_TOTAL')
                                        call timing_on('COMM_TRACER')
         if ( flagstruct%inline_q ) then
                      i_pack(10) = mpp_start_update_domains(q, domain)
         elseif ( it==n_split ) then
                      i_pack(10) = mpp_start_update_domains(q, domain)
         endif
                                       call timing_off('COMM_TRACER')
                                   call timing_off('COMM_TOTAL')
     endif

     if ( .not. hydrostatic ) then

      if ( it==1 ) then
                             call timing_on('COMM_TOTAL')
         w_pack = mpp_start_update_domains(w, domain)
                             call timing_off('COMM_TOTAL')
!$omp parallel do default(shared)
         do j=jsd,jed
            do i=isd,ied
               gz(i,j,npz+1) = zs(i,j)
            enddo
            do k=npz,1,-1
               do i=isd,ied
                  gz(i,j,k) = gz(i,j,k+1) - delz(i,j,k)
               enddo
            enddo
         enddo
                             call timing_on('COMM_TOTAL')
         i_pack(5) = mpp_start_update_domains(gz,  domain)
                             call timing_off('COMM_TOTAL')
      endif

     endif


#ifdef SW_DYNAMICS
     if (test_case>1) then
     if (test_case==9) call case9_forcing1(phis, time_total)
#endif

     if ( it==1 ) then
                                       call timing_on('COMM_TOTAL')
          call mpp_complete_update_domains(i_pack(1), delp, domain)
          call mpp_complete_update_domains(i_pack(2), pt,   domain)
                                      call timing_off('COMM_TOTAL')
          beta_d = 0.
     else
          beta_d = flagstruct%beta
     endif

     if ( it==n_split .and. end_step ) then
       if ( flagstruct%use_old_omega ) then
          allocate ( pem(is-1:ie+1,npz+1,js-1:je+1) )
!$omp parallel do default(shared)
         do j=js-1,je+1
            do i=is-1,ie+1
               pem(i,1,j) = ptop
            enddo
            do k=1,npz
               do i=is-1,ie+1
                  pem(i,k+1,j) = pem(i,k,j) + delp(i,j,k)
               enddo
            enddo
         enddo
       endif
          last_step = .true.
     else
          last_step = .false.
     endif
       
                                                     call timing_on('COMM_TOTAL')
     call mpp_complete_update_domains(i_pack(8), u, v, domain, gridtype=DGRID_NE)
     if( (.not. hydrostatic) .and.  it == 1) then
        call mpp_complete_update_domains(w_pack, w, domain)
     endif
                                                     call timing_off('COMM_TOTAL')

                                                     call timing_on('c_sw')
!$omp parallel do default(shared)
      do k=1,npz
         call c_sw(delpc(isd,jsd,k), delp(isd,jsd,k),  ptc(isd,jsd,k),    &
                      pt(isd,jsd,k),    u(isd,jsd,k),    v(isd,jsd,k),    &
                       w(isd,jsd,k),   uc(isd,jsd,k),   vc(isd,jsd,k),    &
                      ua(isd,jsd,k),   va(isd,jsd,k), omga(isd,jsd,k),    &
                      ut(isd,jsd,k),   vt(isd,jsd,k), divgd(isd,jsd,k),   &
                      flagstruct%nord,   dt2,  hydrostatic,  .true., bd,  &
                      gridstruct, flagstruct)
      enddo
                                                     call timing_off('c_sw')
      if ( flagstruct%nord > 0 ) then
                                                   call timing_on('COMM_TOTAL')
          i_pack(3) = mpp_start_update_domains(divgd, domain, position=CORNER)
                                                  call timing_off('COMM_TOTAL')
      endif

!     if( flagstruct%fill_dp ) call mix_dp(hydrostatic, omga, delpc, ptc, npz, ak, bk, .true., bd)

      if ( hydrostatic ) then
           call geopk(ptop, pe, peln, delpc, pkc, gz, phis, ptc, pkz, npz, akap, .true., gridstruct%nested, .false., npx, npy, flagstruct%a2b_ord, bd)
      else
           if ( it == 1 ) then

                                      call timing_on('COMM_TOTAL')
              call mpp_complete_update_domains(i_pack(5), gz, domain)
                                     call timing_off('COMM_TOTAL')

!$omp parallel do default(shared)
           do k=1,npz+1
              do j=jsd,jed
                 do i=isd,ied
! Save edge heights for update_dz_d
                    zh(i,j,k) = gz(i,j,k)
                 enddo
              enddo
           enddo

        else 
!$omp parallel do default(shared)
           do k=1, npz+1
              do j=jsd,jed
                 do i=isd,ied
                    gz(i,j,k) = zh(i,j,k)
                 enddo
              enddo
           enddo
        endif
                                            call timing_on('UPDATE_DZ_C')
         call update_dz_c(is, ie, js, je, npz, ng, dt2, dp_ref, zs, gridstruct%area, ut, vt, gz, ws3, &
             npx, npy, gridstruct%sw_corner, gridstruct%se_corner, &
             gridstruct%ne_corner, gridstruct%nw_corner, bd, gridstruct%grid_type)
                                            call timing_off('UPDATE_DZ_C')

                                               call timing_on('Riem_Solver')
                                                         call timing_on('Riem_Solver')
           call Riem_Solver_C( ms, dt2,   is,  ie,   js,   je,   npz,   ng,   &
                               akap,  cp,  ptop, phis, omga, ptc,  &
                                delpc, gz,  pkc, ws3, flagstruct%p_fac, &
                                flagstruct%a_imp, flagstruct%scale_z )
                                               call timing_off('Riem_Solver')

!!$                                                                   call timing_on('COMM_TOTAL')
!!$           if ( gridstruct%square_domain ) then
!!$             i_pack(4) = mpp_start_update_domains(pkc, domain,whalo=1,ehalo=1,shalo=1, nhalo=1)
!!$             i_pack(5) = mpp_start_update_domains(gz , domain,whalo=1,ehalo=1,shalo=1, nhalo=1)
!!$             call mpp_complete_update_domains(i_pack(4), pkc,domain,whalo=1,ehalo=1,shalo=1, nhalo=1)
!!$             call mpp_complete_update_domains(i_pack(5), gz ,domain,whalo=1,ehalo=1,shalo=1, nhalo=1)
!!$!            call mpp_update_domains(pkc,domain,whalo=1,ehalo=1,shalo=1, nhalo=1, complete=.false.)
!!$!            call mpp_update_domains(gz ,domain,whalo=1,ehalo=1,shalo=1, nhalo=1, complete=.true.)
!!$           else   
!!$             i_pack(4) = mpp_start_update_domains(pkc, domain)
!!$             i_pack(5) = mpp_start_update_domains(gz,  domain)
!!$             call mpp_complete_update_domains(i_pack(4), pkc, domain)
!!$             call mpp_complete_update_domains(i_pack(5), gz , domain)
!!$!            call mpp_update_domains(pkc, domain, complete=.false.)
!!$!            call mpp_update_domains(gz , domain, complete=.true.)
!!$           endif
!!$                                                                   call timing_off('COMM_TOTAL')
           if (gridstruct%nested) then
                 call nested_grid_BC_apply_intT(delpc, &
                      0, 0, npx, npy, npz, bd, split_timestep_BC+0.5, real(n_split*flagstruct%k_split), &
                      var_east_t0=neststruct%delp_BC%east_t0, &
                      var_west_t0=neststruct%delp_BC%west_t0, &
                      var_north_t0=neststruct%delp_BC%north_t0, &
                      var_south_t0=neststruct%delp_BC%south_t0, &
                      var_east_t1=neststruct%delp_BC%east_t1, &
                      var_west_t1=neststruct%delp_BC%west_t1, &
                      var_north_t1=neststruct%delp_BC%north_t1, &
                      var_south_t1=neststruct%delp_BC%south_t1, &
                      bctype=neststruct%nestbctype, &
                      nsponge=neststruct%nsponge, s_weight=neststruct%s_weight   )
                 call nested_grid_BC_apply_intT(ptc, &
                      0, 0, npx, npy, npz, bd, split_timestep_BC+0.5, real(n_split*flagstruct%k_split), &
                      var_east_t0=neststruct%pt_BC%east_t0, &
                      var_west_t0=neststruct%pt_BC%west_t0, &
                      var_north_t0=neststruct%pt_BC%north_t0, &
                      var_south_t0=neststruct%pt_BC%south_t0, &
                      var_east_t1=neststruct%pt_BC%east_t1, &
                      var_west_t1=neststruct%pt_BC%west_t1, &
                      var_north_t1=neststruct%pt_BC%north_t1, &
                      var_south_t1=neststruct%pt_BC%south_t1, &
                      bctype=neststruct%nestbctype, &
                      nsponge=neststruct%nsponge, s_weight=neststruct%s_weight   )
                 call nested_grid_BC_apply_intT(delz, &
                      0, 0, npx, npy, npz, bd, split_timestep_BC+0.5, real(n_split*flagstruct%k_split), &
                      var_east_t0=neststruct%delz_BC%east_t0, &
                      var_west_t0=neststruct%delz_BC%west_t0, &
                      var_north_t0=neststruct%delz_BC%north_t0, &
                      var_south_t0=neststruct%delz_BC%south_t0, &
                      var_east_t1=neststruct%delz_BC%east_t1, &
                      var_west_t1=neststruct%delz_BC%west_t1, &
                      var_north_t1=neststruct%delz_BC%north_t1, &
                      var_south_t1=neststruct%delz_BC%south_t1, &
                      bctype=neststruct%nestbctype, &
                      nsponge=neststruct%nsponge, s_weight=neststruct%s_weight   )

!!$         do k=1,npz
!!$            call extrapolation_BC( delz(:,:,k), 0, 0, npx, npy, bd)
!!$            call extrapolation_BC(delpc(:,:,k), 0, 0, npx, npy, bd)
!!$            call extrapolation_BC(  ptc(:,:,k), 0, 0, npx, npy, bd)
!!$         end do


              !Compute gz/pkc
              !NOTE: nominally only need to compute quantities one out in the halo for p_grad_c
              !(instead of entire halo)
              call geopk_halo_nh(ptop, grav, akap, cp, delpc, delz, ptc, phis, pkc, gz, pk3, &
                   npx, npy, npz, gridstruct%nested, .false., .false., .false., bd)

           endif

      endif   ! end hydro check

      call p_grad_c(dt2, npz, delpc, pkc, gz, uc, vc, bd, gridstruct%rdxc, gridstruct%rdyc, hydrostatic)

                                                                   call timing_on('COMM_TOTAL')
      i_pack(9) = mpp_start_update_domains(uc, vc, domain, gridtype=CGRID_NE)
                                                     call timing_off('COMM_TOTAL')
#ifdef SW_DYNAMICS
      if (test_case==9) call case9_forcing2(phis)
      endif !test_case>1
#endif

                                                                   call timing_on('COMM_TOTAL')
    if (flagstruct%inline_q .and. nq>0) call mpp_complete_update_domains(i_pack(10), q, domain)
    if (flagstruct%nord > 0) call mpp_complete_update_domains(i_pack(3), divgd,  domain, position=CORNER)
                             call mpp_complete_update_domains(i_pack(9), uc, vc, domain, gridtype=CGRID_NE)

                                                                   call timing_off('COMM_TOTAL')
      if (gridstruct%nested) then
         !On a nested grid we have to do SOMETHING with uc and vc in
         ! the boundary halo, particularly at the corners of the
         ! domain and of each processor element. We must either
         ! apply an interpolated BC, or extrapolate into the
         ! boundary halo
         ! NOTE: 
         !The update_domains calls for uc and vc need to go BEFORE the BCs to ensure cross-restart
         !bitwise-consistent solutions when doing the spatial extrapolation; should not make a
         !difference for interpolated BCs from the coarse grid.


!!$         do k=1,npz
!!$            call extrapolation_BC(uc(:,:,k), 1, 0, npx, npy, bd)
!!$            call extrapolation_BC(vc(:,:,k), 0, 1, npx, npy, bd)
!!$         end do
!!$


         !vc
            call nested_grid_BC_apply_intT(vc, &
                 0, 1, npx, npy, npz, bd, split_timestep_bc+0.5, real(n_split*flagstruct%k_split), & 
                 var_east_t0=neststruct%vc_BC%east_t0, &
                 var_west_t0=neststruct%vc_BC%west_t0, &
                 var_north_t0=neststruct%vc_BC%north_t0, &
                 var_south_t0=neststruct%vc_BC%south_t0, &
                 var_east_t1=neststruct%vc_BC%east_t1, &
                 var_west_t1=neststruct%vc_BC%west_t1, &
                 var_north_t1=neststruct%vc_BC%north_t1, &
                 var_south_t1=neststruct%vc_BC%south_t1, &
                 bctype=neststruct%nestbctype, &
                 nsponge=neststruct%nsponge, s_weight=neststruct%s_weight  )

            !uc
            call nested_grid_BC_apply_intT(uc, &
                 1, 0, npx, npy, npz, bd, split_timestep_bc+0.5, real(n_split*flagstruct%k_split), &
                 var_east_t0=neststruct%uc_BC%east_t0, &
                 var_west_t0=neststruct%uc_BC%west_t0, &
                 var_north_t0=neststruct%uc_BC%north_t0, &
                 var_south_t0=neststruct%uc_BC%south_t0, &
                 var_east_t1=neststruct%uc_BC%east_t1, &
                 var_west_t1=neststruct%uc_BC%west_t1, &
                 var_north_t1=neststruct%uc_BC%north_t1, &
                 var_south_t1=neststruct%uc_BC%south_t1, &
                 bctype=neststruct%nestbctype, &
                 nsponge=neststruct%nsponge, s_weight=neststruct%s_weight   )

      end if

    if ( flagstruct%inline_q ) then
         if (gridstruct%nested) then
            do iq=1,nq
                  call nested_grid_BC_apply_intT(q(isd:ied,jsd:jed,:,iq), &
                       0, 0, npx, npy, npz, bd, split_timestep_BC+1, real(n_split*flagstruct%k_split), &
                       var_east_t0=neststruct%q_BC(iq)%east_t0, &
                       var_west_t0=neststruct%q_BC(iq)%west_t0, &
                       var_north_t0=neststruct%q_BC(iq)%north_t0, &
                       var_south_t0=neststruct%q_BC(iq)%south_t0, &
                       var_east_t1=neststruct%q_BC(iq)%east_t1, &
                       var_west_t1=neststruct%q_BC(iq)%west_t1, &
                       var_north_t1=neststruct%q_BC(iq)%north_t1, &
                       var_south_t1=neststruct%q_BC(iq)%south_t1, &
                       bctype=neststruct%nestbctype, &
                       nsponge=neststruct%nsponge, s_weight=neststruct%s_weight   )
            end do
         end if
      endif

                                                     call timing_on('d_sw')
!$omp parallel do default(shared) private(nord_k, damp_k, d2_divg, dd_divg, hord_m, hord_v, hord_t, hord_p, wk, heat_s, d_con_k)
    do k=1,npz
       hord_m = flagstruct%hord_mt
       hord_t = flagstruct%hord_tm
       hord_v = flagstruct%hord_vt
       hord_p = flagstruct%hord_dp
       nord_k = flagstruct%nord
       nord_v(k) = min(2, flagstruct%nord)
       damp_k = flagstruct%dddmp
       d2_divg = min(0.20, flagstruct%d2_bg*(1.-3.*tanh(0.1*log(pfull(k)/pfull(npz)))))
       dd_divg = flagstruct%d4_bg
       d_con_k = flagstruct%d_con
       if ( flagstruct%do_vort_damp ) then
            damp_vt(k) = flagstruct%vtdm4
       else
            damp_vt(k) = 0.
       endif
       if ( npz==1 .or. flagstruct%n_sponge<0 ) then
           d2_divg = flagstruct%d2_bg
       elseif ( flagstruct%n_sponge==0 ) then
! New Del-2 Sponge layer: formulation
              if ( k==1 ) then
                   damp_k = 0.
                   nord_k = 0
                   d2_divg = max(0.21, flagstruct%d2_bg_k1, flagstruct%d2_bg, 0.05)
                   nord_v(k) = 0
                   damp_vt(k) = flagstruct%d2_bg_k1
              elseif ( k==2 ) then
                   damp_k = 0.
                   nord_k = 0
                   d2_divg = min(0.21, max(1.25*flagstruct%d2_bg_k2, flagstruct%d2_bg, 0.01))
                   nord_v(k) = 0
                   damp_vt(k) = flagstruct%d2_bg_k2
              elseif ( k==3 .and. flagstruct%d2_bg_k2>0.5 ) then
                   damp_k = 0.
                   nord_k = 0
                   d2_divg = max(0.2*flagstruct%d2_bg_k2, flagstruct%d2_bg, 0.01)
                   nord_v(k) = 0
                   damp_vt(k) = 0.2*flagstruct%d2_bg_k2
              endif

              if ( damp_vt(k) < 0.01 .and. nord_k>0 ) d_con_k = 0.
              if ( nord_v(k)==0 .and. damp_vt(k)>0.01 ) then
                   hord_t = 6
                   hord_v = 6
              endif
       else
           if( k <= flagstruct%n_sponge .and. npz>16 ) then
! Apply first order scheme for damping the sponge layer
               hord_m = 1
               hord_v = 1
               hord_t = 1
               hord_p = 1
               nord_k = 0
               damp_k = flagstruct%damp_k_k1
               d2_divg = min(0.20, flagstruct%d2_bg_k1*flagstruct%d2_bg)   ! 0.25 is the stability limit
               d2_divg = max(flagstruct%d2_divg_max_k1, d2_divg)
           elseif( k == flagstruct%n_sponge+1 .and. npz>24 ) then
               hord_v = 2
               hord_t = 2
               hord_p = 2
               nord_k = max(0, flagstruct%nord-1)
               d2_divg = min(0.20, flagstruct%d2_bg_k2*flagstruct%d2_bg)
               d2_divg = max(flagstruct%d2_divg_max_k2, d2_divg)
               if ( flagstruct%nord > 1 ) then
                    damp_k = 0.
               else
                    damp_k = flagstruct%damp_k_k2
               endif
           endif
       endif
       damp_k = max(damp_k, flagstruct%dddmp)

       if( .not. flagstruct%use_old_omega .and. last_step ) then
! Average horizontal "convergence" to cell center
            do j=js,je
               do i=is,ie
                  omga(i,j,k) = delp(i,j,k)
               enddo
            enddo
       endif

!--- external mode divergence damping ---
       if ( flagstruct%d_ext > 0. )  &
            call a2b_ord2(delp(isd,jsd,k), wk, gridstruct, npx, npy, is,    &
                          ie, js, je, ng, .false.)

       call d_sw(vt(isd,jsd,k), delp(isd,jsd,k), ptc(isd,jsd,k),  pt(isd,jsd,k),      &
                  u(isd,jsd,k),    v(isd,jsd,k),   w(isd,jsd,k),  uc(isd,jsd,k),      &
                  vc(isd,jsd,k),   ua(isd,jsd,k),  va(isd,jsd,k), divgd(isd,jsd,k),   &
                  mfx(is, js, k),  mfy(is, js, k),  cx(is, jsd,k),  cy(isd,js, k),    &
                  crx(is, jsd,k),  cry(isd,js, k), xfx(is, jsd,k), yfx(isd,js, k),    &
                  heat_s, zvir, sphum, nq,  q,  k,  npz, flagstruct%inline_q,  dt,               &
                  flagstruct%hord_tr, hord_m, hord_v, hord_t, hord_p, nord_k, nord_v(k), damp_k, &
                  d2_divg, dd_divg, damp_vt(k), d_con_k, hydrostatic, gridstruct, flagstruct, bd)

       if( .not. flagstruct%use_old_omega .and. last_step ) then
! Average horizontal "convergence" to cell center
            do j=js,je
               do i=is,ie
                  omga(i,j,k) = omga(i,j,k)*(xfx(i,j,k)-xfx(i+1,j,k)+yfx(i,j,k)-yfx(i,j+1,k))*rarea(i,j)*rdt
               enddo
            enddo
       endif

       if ( flagstruct%d_ext > 0. ) then
            do j=js,jep1
               do i=is,iep1
                  ptc(i,j,k) = wk(i,j)    ! delp at cell corners
               enddo
            enddo
       endif
       if ( d_con_k > 1.0E-5 ) then
! Average horizontal "convergence" to cell center
            do j=js,je
               do i=is,ie
                  heat_source(i,j,k) = heat_source(i,j,k) + heat_s(i,j)
               enddo
            enddo
       endif
    enddo           ! end openMP k-loop
                                                     call timing_off('d_sw')

    if( flagstruct%fill_dp ) call mix_dp(hydrostatic, w, delp, pt, npz, ak, bk, .false., flagstruct%fv_debug, bd)

                                                             call timing_on('COMM_TOTAL')
    i_pack(1) = mpp_start_update_domains(delp, domain)
    i_pack(2) = mpp_start_update_domains(pt,   domain)
    if ( .not. hydrostatic )  then
       w_pack = mpp_start_update_domains(w, domain)
    endif
                                                             call timing_off('COMM_TOTAL')

    if ( flagstruct%d_ext > 0. ) then
          d2_divg = flagstruct%d_ext * gridstruct%da_min_c
!$omp parallel do default(shared)
          do j=js,jep1
              do i=is,iep1
                    wk(i,j) = ptc(i,j,1)
                 divg2(i,j) = wk(i,j)*vt(i,j,1)
              enddo
              do k=2,npz
                 do i=is,iep1
                       wk(i,j) =    wk(i,j) + ptc(i,j,k)
                    divg2(i,j) = divg2(i,j) + ptc(i,j,k)*vt(i,j,k)
                 enddo
              enddo
              do i=is,iep1
                 divg2(i,j) = d2_divg*divg2(i,j)/wk(i,j)
              enddo
          enddo
    else
        divg2(:,:) = 0.
    endif

                                       call timing_on('COMM_TOTAL')
     call mpp_complete_update_domains(i_pack(1), delp, domain)
     call mpp_complete_update_domains(i_pack(2), pt,   domain)
                                       call timing_off('COMM_TOTAL')
     if ( hydrostatic ) then
          call geopk(ptop, pe, peln, delp, pkc, gz, phis, pt, pkz, npz, akap, .false., gridstruct%nested, .true., npx, npy, flagstruct%a2b_ord, bd)
          !Want to move this block into the hydro/nonhydro branch above and merge the two if structures
          if (gridstruct%nested) then

             call nested_grid_BC_apply_intT(delp, &
                  0, 0, npx, npy, npz, bd, split_timestep_BC+1, real(n_split*flagstruct%k_split), &
                  var_east_t0=neststruct%delp_BC%east_t0, &
                  var_west_t0=neststruct%delp_BC%west_t0, &
                  var_north_t0=neststruct%delp_BC%north_t0, &
                  var_south_t0=neststruct%delp_BC%south_t0, &
                  var_east_t1=neststruct%delp_BC%east_t1, &
                  var_west_t1=neststruct%delp_BC%west_t1, &
                  var_north_t1=neststruct%delp_BC%north_t1, &
                  var_south_t1=neststruct%delp_BC%south_t1, &
                  bctype=neststruct%nestbctype, &
                  nsponge=neststruct%nsponge, s_weight=neststruct%s_weight   )

#ifndef SW_DYNAMICS

             call nested_grid_BC_apply_intT(pt, &
                  0, 0, npx, npy, npz, bd, split_timestep_BC+1, real(n_split*flagstruct%k_split), &
                  var_east_t0=neststruct%pt_BC%east_t0, &
                  var_west_t0=neststruct%pt_BC%west_t0, &
                  var_north_t0=neststruct%pt_BC%north_t0, &
                  var_south_t0=neststruct%pt_BC%south_t0, &
                  var_east_t1=neststruct%pt_BC%east_t1, &
                  var_west_t1=neststruct%pt_BC%west_t1, &
                  var_north_t1=neststruct%pt_BC%north_t1, &
                  var_south_t1=neststruct%pt_BC%south_t1, &
                  bctype=neststruct%nestbctype, &
                  nsponge=neststruct%nsponge, s_weight=neststruct%s_weight   )

#endif

          end if
          !Is this extra geopk call necessary?
          call geopk(ptop, pe, peln, delp, pkc, gz, phis, pt, pkz, npz, akap, &
               .false., gridstruct%nested, .true., npx, npy, flagstruct%a2b_ord, bd)
     else
                                            call timing_on('UPDATE_DZ')
        call update_dz_d(nord_v, damp_vt, flagstruct%hord_tm, is, ie, js, je, npz, ng, npx, npy, gridstruct%area,  &
                         gridstruct%rarea, dp_ref, zs, zh, crx, cry, xfx, yfx, delz, ws, rdt, gridstruct, bd)
                                            call timing_off('UPDATE_DZ')

        if (idiag%id_ws>0 .and. last_step) then
!           call prt_maxmin('WS', ws, is, ie, js, je, 0, 1, 1., master)
            used=send_data(idiag%id_ws, ws, fv_time)
        endif

                                call timing_on('COMM_TOTAL')
        i_pack(5) = mpp_start_update_domains(zh,  domain)
        call mpp_complete_update_domains(i_pack(5), zh,  domain)
                                call timing_off('COMM_TOTAL')

        !Want to move this block into the hydro/nonhydro branch above and merge the two if structures
        if (gridstruct%nested) then

           call nested_grid_BC_apply_intT(delp, &
                0, 0, npx, npy, npz, bd, split_timestep_BC+1, real(n_split*flagstruct%k_split), &
                var_east_t0=neststruct%delp_BC%east_t0, &
                var_west_t0=neststruct%delp_BC%west_t0, &
                var_north_t0=neststruct%delp_BC%north_t0, &
                var_south_t0=neststruct%delp_BC%south_t0, &
                var_east_t1=neststruct%delp_BC%east_t1, &
                var_west_t1=neststruct%delp_BC%west_t1, &
                var_north_t1=neststruct%delp_BC%north_t1, &
                var_south_t1=neststruct%delp_BC%south_t1, &
                bctype=neststruct%nestbctype, &
                nsponge=neststruct%nsponge, s_weight=neststruct%s_weight   )

#ifndef SW_DYNAMICS

           call nested_grid_BC_apply_intT(pt, &
                0, 0, npx, npy, npz, bd, split_timestep_BC+1, real(n_split*flagstruct%k_split), &
                var_east_t0=neststruct%pt_BC%east_t0, &
                var_west_t0=neststruct%pt_BC%west_t0, &
                var_north_t0=neststruct%pt_BC%north_t0, &
                var_south_t0=neststruct%pt_BC%south_t0, &
                var_east_t1=neststruct%pt_BC%east_t1, &
                var_west_t1=neststruct%pt_BC%west_t1, &
                var_north_t1=neststruct%pt_BC%north_t1, &
                var_south_t1=neststruct%pt_BC%south_t1, &
                bctype=neststruct%nestbctype, &
                nsponge=neststruct%nsponge, s_weight=neststruct%s_weight   )

#endif

        end if

!-----------------------------------------------------------
! pkc is non-hydrostatic perturbation pressure
!-----------------------------------------------------------
!$omp parallel do default(shared)
        do j=jsd, jed
           do i=isd, ied
              ws3(i,j) = ( zs(i,j) - zh(i,j,npz+1) ) * rdt
           enddo
        enddo


                                call timing_on('COMM_TOTAL')
        call mpp_complete_update_domains(w_pack, w, domain)
                                call timing_off('COMM_TOTAL')

                                                         call timing_on('Riem_Solver')
        call Riem_Solver3(flagstruct%m_split, dt,  is,  ie,   js,   je, npz, ng,     &
                         isd, ied, jsd, jed, &
                         akap, cp,  ptop, phis, w,  delz, pt, delp, zh,    &
                         gz,   pkc, pk3, pk, pe, peln, ws3, &
                         flagstruct%p_fac, flagstruct%a_imp, &
                         flagstruct%scale_z, flagstruct%use_logp, remap_step)
                                                         call timing_off('Riem_Solver')

     endif    ! end hydro check

     if (gridstruct%nested) then
!!! How to best handle this?
                                       call timing_on('COMM_TOTAL')
        call mpp_update_domains(pkc, domain, complete=.false.)
        call mpp_update_domains(gz , domain, complete=.true.)
                                       call timing_off('COMM_TOTAL')
     end if

#ifdef SW_DYNAMICS
      if (test_case > 1) then
#else
      if ( hydrostatic .and. remap_step ) then
!$omp parallel do default(shared)
           do k=1,npz+1
              do j=js,je
                 do i=is,ie
                    pk(i,j,k) = pkc(i,j,k)
                 enddo
              enddo
           enddo
      endif
#endif

!----------------------------
! Compute pressure gradient:
!----------------------------
    if ( hydrostatic ) then
       if ( flagstruct%beta > 0. ) then
          call grad1_p_update(divg2, u, v, du, dv, pkc, gz, dt, ng, gridstruct, bd, npx, npy, npz, ptop, beta_d, flagstruct%a2b_ord)
       else
          call one_grad_p(u, v, pkc, gz, divg2, delp, dt, ng, gridstruct, bd, npx, npy, npz, ptop, hydrostatic, flagstruct%a2b_ord, flagstruct%d_ext)
       endif

    else

       if (gridstruct%nested) then
           call nested_grid_BC_apply_intT(delz, &
                    0, 0, npx, npy, npz, bd, split_timestep_BC+1., real(n_split*flagstruct%k_split), &
                    var_east_t0=neststruct%delz_BC%east_t0, &
                    var_west_t0=neststruct%delz_BC%west_t0, &
                    var_north_t0=neststruct%delz_BC%north_t0, &
                    var_south_t0=neststruct%delz_BC%south_t0, &
                    var_east_t1=neststruct%delz_BC%east_t1, &
                    var_west_t1=neststruct%delz_BC%west_t1, &
                    var_north_t1=neststruct%delz_BC%north_t1, &
                    var_south_t1=neststruct%delz_BC%south_t1, &
                    bctype=neststruct%nestbctype, &
                    nsponge=neststruct%nsponge, s_weight=neststruct%s_weight   )

            !Compute gz/pkc/pk3; note that now pkc should be nonhydro pert'n pressure
            call geopk_halo_nh(ptop, grav, akap, cp, delp, delz, pt, phis, pkc, gz, pk3, npx, npy, npz, gridstruct%nested, .true., .true., .true., bd)

       endif

       if ( flagstruct%beta > 0. ) then
          call split_p_grad(du, dv, u, v, pkc, gz, delp, pk3, divg2, beta_d, dt, ng, gridstruct, bd, npx, npy, npz)
       else
          call nh_p_grad(u, v, pkc, gz, delp, pk3, divg2, dt, ng, gridstruct, bd, npx, npy, npz)
       endif

   endif

!-------------------------------------------------------------------------------------------------------
    if ( flagstruct%breed_vortex_inline ) then
#if defined (ADA_NUDGE)
         call breed_slp_inline_ada( it, dt, npz, ak, bk, phis, pe, pk, peln, pkz,     &
                                delp, u, v, pt, q, flagstruct%nwat, zvir, gridstruct, ks, domain, bd )
#else
         call breed_slp_inline( it, dt, npz, ak, bk, phis, pe, pk, peln, pkz,     &
                                delp, u, v, pt, q, flagstruct%nwat, zvir, gridstruct, ks, domain, bd )

#endif
    endif
!-------------------------------------------------------------------------------------------------------

                                                     call timing_on('COMM_TOTAL')
    if( it==n_split .and. gridstruct%grid_type<4 .and. .not. gridstruct%nested) then
! Prevent accumulation of rounding errors at overlapped domain edges:
          call mpp_get_boundary(u, v, domain, wbuffery=wbuffer, ebuffery=ebuffer,  &
                            sbufferx=sbuffer, nbufferx=nbuffer, gridtype=DGRID_NE )
!$omp parallel do default(shared)
          do k=1,npz
             do i=is,ie
                u(i,je+1,k) = nbuffer(i-is+1,k)
             enddo
             do j=js,je
                v(ie+1,j,k) = ebuffer(j-js+1,k)
             enddo
          enddo

    endif

    if ( it/=n_split)   &
         i_pack(8) = mpp_start_update_domains(u, v, domain, gridtype=DGRID_NE)
                                                     call timing_off('COMM_TOTAL')

#ifdef SW_DYNAMICS
    endif
#endif
      if ( gridstruct%nested ) then
         neststruct%nest_timestep = neststruct%nest_timestep + 1
      endif

#ifdef SW_DYNAMICS
#else
    if ( last_step ) then
      if ( flagstruct%use_old_omega ) then
!$omp parallel do default(shared) private(om2d)
         do k=1,npz
            do j=js,je
               do i=is,ie
                  omga(i,j,k) = (pe(i,k+1,j) - pem(i,k+1,j)) * rdt
               enddo
            enddo
         enddo
!------------------------------
! Compute the "advective term"
!------------------------------
         call adv_pe(ua, va, pem, omga, gridstruct, bd, npx, npy,  npz, ng)
      else
!$omp parallel do default(shared) private(om2d)
         do j=js,je
            do k=1,npz
               do i=is,ie
                  om2d(i,k) = omga(i,j,k)
               enddo
            enddo
            do k=2,npz
               do i=is,ie
                  om2d(i,k) = om2d(i,k-1) + omga(i,j,k)
               enddo
            enddo
            do k=2,npz
               do i=is,ie
                  omga(i,j,k) = om2d(i,k)
               enddo
            enddo
         enddo
      endif
      if (idiag%id_ws>0 .and. hydrostatic) then
!$omp parallel do default(shared) 
          do j=js,je
             do i=is,ie
                ws(i,j) = delz(i,j,npz)/delp(i,j,npz) * omga(i,j,npz)
             enddo
          enddo
          used=send_data(idiag%id_ws, ws, fv_time)
      endif
    endif
#endif

    if (gridstruct%nested) then
#ifdef SW_DYNAMICS
#else



         if (.not. hydrostatic) then
               call nested_grid_BC_apply_intT(w, &
                    0, 0, npx, npy, npz, bd, split_timestep_BC, real(n_split*flagstruct%k_split), &
                    var_east_t0=neststruct%w_BC%east_t0, &
                    var_west_t0=neststruct%w_BC%west_t0, &
                    var_north_t0=neststruct%w_BC%north_t0, &
                    var_south_t0=neststruct%w_BC%south_t0, &
                    var_east_t1=neststruct%w_BC%east_t1, &
                    var_west_t1=neststruct%w_BC%west_t1, &
                    var_north_t1=neststruct%w_BC%north_t1, &
                    var_south_t1=neststruct%w_BC%south_t1, &
                    bctype=neststruct%nestbctype, &
                    nsponge=neststruct%nsponge, s_weight=neststruct%s_weight   )
         end if
#endif
         
         
            call nested_grid_BC_apply_intT(u, &
                 0, 1, npx, npy, npz, bd, split_timestep_BC, real(n_split*flagstruct%k_split), &
                 var_east_t0=neststruct%u_BC%east_t0, &
                 var_west_t0=neststruct%u_BC%west_t0, &
                 var_north_t0=neststruct%u_BC%north_t0, &
                 var_south_t0=neststruct%u_BC%south_t0, &
                 var_east_t1=neststruct%u_BC%east_t1, &
                 var_west_t1=neststruct%u_BC%west_t1, &
                 var_north_t1=neststruct%u_BC%north_t1, &
                 var_south_t1=neststruct%u_BC%south_t1, &
                 bctype=neststruct%nestbctype, &
                 nsponge=neststruct%nsponge, s_weight=neststruct%s_weight   )

         !v
            call nested_grid_BC_apply_intT(v, &
                 1, 0, npx, npy, npz, bd, split_timestep_BC, real(n_split*flagstruct%k_split), &
                 var_east_t0=neststruct%v_BC%east_t0, &
                 var_west_t0=neststruct%v_BC%west_t0, &
                 var_north_t0=neststruct%v_BC%north_t0, &
                 var_south_t0=neststruct%v_BC%south_t0, &
                 var_east_t1=neststruct%v_BC%east_t1, &
                 var_west_t1=neststruct%v_BC%west_t1, &
                 var_north_t1=neststruct%v_BC%north_t1, &
                 var_south_t1=neststruct%v_BC%south_t1, &
                 bctype=neststruct%nestbctype, &
                 nsponge=neststruct%nsponge, s_weight=neststruct%s_weight   )

      end if

!-----------------------------------------------------
  enddo   ! time split loop
!-----------------------------------------------------

  if ( flagstruct%fv_debug ) then
       if(is_master()) write(*,*) 'End of n_split loop'
  endif


  if ( n_con/=0 .and. flagstruct%d_con > 1.e-5 ) then
       call del2_cubed(heat_source, 0.20*gridstruct%da_min, gridstruct, domain, npx, npy, npz, 3, bd)

! Note: pt here is cp*(Virtual_Temperature/pkz)
    if ( hydrostatic ) then
!
! del(Cp*T) = - del(KE)
!
!$omp parallel do default(shared)
       do k=1,n_con
          do j=js,je
             do i=is,ie
                pt(i,j,k) = pt(i,j,k) + heat_source(i,j,k)/(delp(i,j,k)*pkz(i,j,k))
             enddo
          enddo
       enddo
    else
!$omp parallel do default(shared) private(tmcp)
       do k=1,n_con
          do j=js,je
             do i=is,ie
                pkz(i,j,k) = exp( k1k*log(kapag*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)) )
! tmcp = termperature_v * cp
                tmcp = pt(i,j,k) * pkz(i,j,k)
                delz(i,j,k) = delz(i,j,k) / tmcp
                tmcp = tmcp + heat_source(i,j,k)/delp(i,j,k)
                pt(i,j,k) = tmcp / pkz(i,j,k)
                delz(i,j,k) = delz(i,j,k) * tmcp
             enddo
          enddo
       enddo
    endif
       !deallocate( heat_source )    
    endif
    if (allocated(heat_source)) deallocate( heat_source ) !If ncon == 0 but d_con > 1.e-5, this would not be deallocated in earlier versions of the code


  if ( end_step ) then
    deallocate(    gz )
    deallocate(   ptc )
    deallocate(   crx )
    deallocate(   xfx )
    deallocate(   cry )
    deallocate(   yfx )
    deallocate( divgd )
#ifndef PKC
    deallocate(   pkc )
#endif
    deallocate( delpc )

    if( allocated(ut))   deallocate( ut )
    if( allocated(vt))   deallocate( vt )
    if ( allocated (du) ) deallocate( du )
    if ( allocated (dv) ) deallocate( dv )
    if ( .not. hydrostatic ) then
         deallocate( zh )
         if( allocated(pk3) )   deallocate ( pk3 )
    endif

  endif
  if( allocated(pem) )   deallocate ( pem )

if ( flagstruct%fv_debug ) then
   if(is_master()) write(*,*) 'End of dyn_core'
endif

end subroutine dyn_core


subroutine adv_pe(ua, va, pem, om, gridstruct, bd, npx, npy, npz, ng)

integer, intent(in) :: npx, npy, npz, ng
type(fv_grid_bounds_type), intent(IN) :: bd
! Contra-variant wind components:
real, intent(in), dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz):: ua, va
! Pressure at edges:
real, intent(in) :: pem(bd%is-1:bd%ie+1,1:npz+1,bd%js-1:bd%je+1)
real, intent(inout) :: om(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
type(fv_grid_type), intent(INOUT), target :: gridstruct

! Local:
real, dimension(bd%is:bd%ie,bd%js:bd%je):: up, vp
real v3(3,bd%is:bd%ie,bd%js:bd%je)

real pin(bd%isd:bd%ied,bd%jsd:bd%jed)
real  pb(bd%isd:bd%ied,bd%jsd:bd%jed)

real grad(3,bd%is:bd%ie,bd%js:bd%je)
real pdx(3,bd%is:bd%ie,bd%js:bd%je+1)
real pdy(3,bd%is:bd%ie+1,bd%js:bd%je)
real, pointer, dimension(:,:) :: rarea, dx, dy
real, pointer, dimension(:,:,:) :: ec1, ec2, en1, en2

integer :: i,j,k, n

integer :: is,  ie,  js,  je

rarea => gridstruct%rarea
dx    => gridstruct%dx
dy    => gridstruct%dy
ec1   => gridstruct%ec1
ec2   => gridstruct%ec2
en1   => gridstruct%en1
en2   => gridstruct%en2

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je

!$omp parallel do default(shared) private(i, j, k, n, pdx, pdy, pin, pb, up, vp, grad, v3)
do k=1,npz
   if ( k==npz ) then
      do j=js,je
         do i=is,ie
            up(i,j) = ua(i,j,npz)
            vp(i,j) = va(i,j,npz)
         enddo
      enddo
   else
      do j=js,je
         do i=is,ie
            up(i,j) = 0.5*(ua(i,j,k)+ua(i,j,k+1))
            vp(i,j) = 0.5*(va(i,j,k)+va(i,j,k+1))
         enddo
      enddo
   endif

   ! Compute Vect wind:
   do j=js,je
      do i=is,ie
         do n=1,3
            v3(n,i,j) = up(i,j)*ec1(n,i,j) + vp(i,j)*ec2(n,i,j) 
         enddo
      enddo
   enddo

   do j=js-1,je+1
      do i=is-1,ie+1
         pin(i,j) = pem(i,k+1,j)
      enddo
   enddo

   ! Compute pe at 4 cell corners:
   call a2b_ord2(pin, pb, gridstruct, npx, npy, is, ie, js, je, ng)


   do j=js,je+1
      do i=is,ie
         do n=1,3
            pdx(n,i,j) = (pb(i,j)+pb(i+1,j))*dx(i,j)*en1(n,i,j)
         enddo
      enddo
   enddo
   do j=js,je
      do i=is,ie+1
         do n=1,3
            pdy(n,i,j) = (pb(i,j)+pb(i,j+1))*dy(i,j)*en2(n,i,j)
         enddo
      enddo
   enddo

   ! Compute grad (pe) by Green's theorem
   do j=js,je
      do i=is,ie
         do n=1,3
            grad(n,i,j) = pdx(n,i,j+1) - pdx(n,i,j) - pdy(n,i,j) + pdy(n,i+1,j)
         enddo
      enddo
   enddo

   ! Compute inner product: V3 * grad (pe)
   do j=js,je
      do i=is,ie
         om(i,j,k) = om(i,j,k) + 0.5*rarea(i,j)*(v3(1,i,j)*grad(1,i,j) +   &
              v3(2,i,j)*grad(2,i,j) + v3(3,i,j)*grad(3,i,j))
      enddo
   enddo
enddo

end subroutine adv_pe




subroutine p_grad_c(dt2, npz, delpc, pkc, gz, uc, vc, bd, rdxc, rdyc, hydrostatic)

integer, intent(in):: npz
real,    intent(in):: dt2
type(fv_grid_bounds_type), intent(IN) :: bd
real, intent(in), dimension(bd%isd:, bd%jsd: ,:  ):: delpc
! pkc is pe**cappa     if hydrostatic
! pkc is full pressure if non-hydrostatic
real, intent(in), dimension(bd%isd:bd%ied, bd%jsd:bd%jed ,npz+1):: pkc, gz
real, intent(inout):: uc(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz)
real, intent(inout):: vc(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz)
real, intent(IN) :: rdxc(bd%isd:bd%ied+1,bd%jsd:bd%jed+1)
real, intent(IN) :: rdyc(bd%isd:bd%ied  ,bd%jsd:bd%jed)
logical, intent(in):: hydrostatic
! Local:
real:: wk(bd%is-1:bd%ie+1,bd%js-1:bd%je+1)
integer:: i,j,k

integer :: is,  ie,  js,  je

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je


!$omp parallel do default(shared) private(wk)
do k=1,npz

   if ( hydrostatic ) then
      do j=js-1,je+1
         do i=is-1,ie+1
            wk(i,j) = pkc(i,j,k+1) - pkc(i,j,k)
         enddo
      enddo
   else
      do j=js-1,je+1
         do i=is-1,ie+1
            wk(i,j) = delpc(i,j,k)
         enddo
      enddo
   endif

   do j=js,je
      do i=is,ie+1
         uc(i,j,k) = uc(i,j,k) + dt2*rdxc(i,j) / (wk(i-1,j)+wk(i,j)) *   &
              ( (gz(i-1,j,k+1)-gz(i,j,k  ))*(pkc(i,j,k+1)-pkc(i-1,j,k))  &
              + (gz(i-1,j,k) - gz(i,j,k+1))*(pkc(i-1,j,k+1)-pkc(i,j,k)) )
      enddo
   enddo
   do j=js,je+1
      do i=is,ie
         vc(i,j,k) = vc(i,j,k) + dt2*rdyc(i,j) / (wk(i,j-1)+wk(i,j)) *   &
              ( (gz(i,j-1,k+1)-gz(i,j,k  ))*(pkc(i,j,k+1)-pkc(i,j-1,k))  &
              + (gz(i,j-1,k) - gz(i,j,k+1))*(pkc(i,j-1,k+1)-pkc(i,j,k)) )
      enddo
   enddo
enddo

end subroutine p_grad_c


subroutine nh_p_grad(u, v, pp, gz, delp, pk, divg2, dt, ng, gridstruct, bd, npx, npy, npz)
integer, intent(IN) :: ng, npx, npy, npz
real,    intent(IN) :: dt
type(fv_grid_bounds_type), intent(IN) :: bd
real,    intent(in) :: divg2(bd%is:bd%ie+1, bd%js:bd%je+1)
real, intent(inout) ::  delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
real, intent(inout) ::    pp(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)  ! perturbation pressure
real, intent(inout) ::    pk(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)  ! p**kappa
real, intent(inout) ::    gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)  ! g * h
real, intent(inout) ::     u(bd%isd:bd%ied,  bd%jsd:bd%jed+1,npz) 
real, intent(inout) ::     v(bd%isd:bd%ied+1,bd%jsd:bd%jed,  npz)
    type(fv_grid_type), intent(INOUT), target :: gridstruct
! Local:
real wk1(bd%isd:bd%ied, bd%jsd:bd%jed)
real  wk(bd%is: bd%ie+1,bd%js: bd%je+1)
real du, dv
integer i,j,k
real, pointer, dimension(:,:) :: rdx,rdy

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

      rdx => gridstruct%rdx
      rdy => gridstruct%rdy

!Remember that not all compilers set pp to zero by default
!$omp parallel do default(shared)
do j=js,je+1
   do i=is,ie+1
      pp(i,j,1) = 0.
      pk(i,j,1) = ptk
   enddo
enddo

!$omp parallel do default(shared) private(wk1)
do k=1,npz+1
   if ( k/=1 ) then
      call a2b_ord4(pp(isd,jsd,k), wk1, gridstruct, npx, npy, is, ie, js, je, ng, .true.)
      call a2b_ord4(pk(isd,jsd,k), wk1, gridstruct, npx, npy, is, ie, js, je, ng, .true.)
   endif
   call a2b_ord4( gz(isd,jsd,k), wk1, gridstruct, npx, npy, is, ie, js, je, ng, .true.)
enddo

!$omp parallel do default(shared) private(wk1, wk, du, dv)
do k=1,npz
   call a2b_ord4(delp(isd,jsd,k), wk1, gridstruct, npx, npy, is, ie, js, je, ng)
   do j=js,je+1
      do i=is,ie+1
         wk(i,j) = pk(i,j,k+1) - pk(i,j,k)
      enddo
   enddo
   do j=js,je+1
      do i=is,ie
         ! hydrostatic contributions from past time-step already added in the "beta" part
         ! Current gradient from "hydrostatic" components:
         du =  dt / (wk(i,j)+wk(i+1,j)) *   &
              ((gz(i,j,k+1)-gz(i+1,j,k))*(pk(i+1,j,k+1)-pk(i,j,k)) +  &
              (gz(i,j,k)-gz(i+1,j,k+1))*(pk(i,j,k+1)-pk(i+1,j,k)))
         ! Non-hydrostatic contribution for half time-step
         u(i,j,k) = (u(i,j,k) + divg2(i,j)-divg2(i+1,j) + du +  &
              dt/(wk1(i,j)+wk1(i+1,j)) *  &
              ((gz(i,j,k+1)-gz(i+1,j,k))*(pp(i+1,j,k+1)-pp(i,j,k))    &
              + (gz(i,j,k)-gz(i+1,j,k+1))*(pp(i,j,k+1)-pp(i+1,j,k))))*rdx(i,j)
      enddo
   enddo
   do j=js,je
      do i=is,ie+1
         ! Current gradient from "hydrostatic" components:
         dv = dt / (wk(i,j)+wk(i,j+1)) *   &
              ((gz(i,j,k+1)-gz(i,j+1,k))*(pk(i,j+1,k+1)-pk(i,j,k)) +  &
              (gz(i,j,k)-gz(i,j+1,k+1))*(pk(i,j,k+1)-pk(i,j+1,k)))
         ! Non-hydrostatic contribution for half time-step
         v(i,j,k) = (v(i,j,k)  + divg2(i,j)-divg2(i,j+1) + dv + &
              dt/(wk1(i,j)+wk1(i,j+1)) *  &
              ((gz(i,j,k+1)-gz(i,j+1,k))*(pp(i,j+1,k+1)-pp(i,j,k))   &
              + (gz(i,j,k)-gz(i,j+1,k+1))*(pp(i,j,k+1)-pp(i,j+1,k))))*rdy(i,j)
      enddo
   enddo

enddo    ! end k-loop
end subroutine nh_p_grad


subroutine split_p_grad(du, dv, u, v, pp, gz, delp, pk, divg2, beta, dt, ng, gridstruct, bd, npx, npy, npz)
integer, intent(IN) :: ng, npx, npy, npz
real,    intent(IN) :: beta, dt
type(fv_grid_bounds_type), intent(IN) :: bd
real,    intent(in) :: divg2(bd%is:bd%ie+1, bd%js:bd%je+1)
real, intent(inout) ::  delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
real, intent(inout) ::    pp(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)  ! perturbation pressure
real, intent(inout) ::    pk(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)  ! p**kappa
real, intent(inout) ::    gz(bd%isd:bd%ied, bd%jsd:bd%jed, npz+1)  ! g * h
real, intent(inout) ::    du(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) 
real, intent(inout) ::    dv(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz)
real, intent(inout) ::     u(bd%isd:bd%ied,  bd%jsd:bd%jed+1,npz) 
real, intent(inout) ::     v(bd%isd:bd%ied+1,bd%jsd:bd%jed,  npz)
type(fv_grid_type), intent(INOUT), target :: gridstruct
! Local:
real wk1(bd%isd:bd%ied, bd%jsd:bd%jed)
real  wk(bd%is: bd%ie+1,bd%js: bd%je+1)
real  alpha
integer i,j,k

real, pointer, dimension(:,:) :: rdx,rdy

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

      rdx => gridstruct%rdx
      rdy => gridstruct%rdy

alpha = 1. - beta

!$omp parallel do default(shared)
do j=js,je+1
   do i=is,ie+1
      pp(i,j,1) = 0.
      pk(i,j,1) = ptk
   enddo
enddo

!$omp parallel do default(shared) private(wk1)
do k=1,npz+1
   if ( k/=1 ) then
      call a2b_ord4(pp(isd,jsd,k), wk1, gridstruct, npx, npy, is, ie, js, je, ng, .true.)
      call a2b_ord4(pk(isd,jsd,k), wk1, gridstruct, npx, npy, is, ie, js, je, ng, .true.)
   endif
   call a2b_ord4( gz(isd,jsd,k), wk1, gridstruct, npx, npy, is, ie, js, je, ng, .true.)
enddo

!$omp parallel do default(shared) private(wk1, wk)
do k=1,npz
   call a2b_ord4(delp(isd,jsd,k), wk1, gridstruct, npx, npy, is, ie, js, je, ng)

   do j=js,je+1
      do i=is,ie+1
         wk(i,j) = pk(i,j,k+1) - pk(i,j,k)
      enddo
   enddo

   do j=js,je+1
      do i=is,ie
         u(i,j,k) = u(i,j,k) + beta*du(i,j,k)
         ! hydrostatic contributions from past time-step already added in the "beta" part
         ! Current gradient from "hydrostatic" components:
         !---------------------------------------------------------------------------------
         du(i,j,k) =  dt / (wk(i,j)+wk(i+1,j)) *   &
              ((gz(i,j,k+1)-gz(i+1,j,k))*(pk(i+1,j,k+1)-pk(i,j,k)) +  &
              (gz(i,j,k)-gz(i+1,j,k+1))*(pk(i,j,k+1)-pk(i+1,j,k)))
         !---------------------------------------------------------------------------------
         ! Non-hydrostatic contribution for half time-step
!!$#ifdef GEOPK_CHECK
!!$         if (abs(pp(i,j,k) ) > 1.e5 ) then
!!$            print*, mpp_pe(),i,j,k, 'P Pert:', pp(i,j,k)
!!$         endif
!!$#endif
         u(i,j,k) = (u(i,j,k) + divg2(i,j)-divg2(i+1,j) + alpha*du(i,j,k) +  &
              dt/(wk1(i,j)+wk1(i+1,j)) *  &
              ((gz(i,j,k+1)-gz(i+1,j,k))*(pp(i+1,j,k+1)-pp(i,j,k))    &
              + (gz(i,j,k)-gz(i+1,j,k+1))*(pp(i,j,k+1)-pp(i+1,j,k))))*rdx(i,j)
      enddo
   enddo
   do j=js,je
      do i=is,ie+1
         v(i,j,k) = v(i,j,k) + beta*dv(i,j,k)
         ! Current gradient from "hydrostatic" components:
         !---------------------------------------------------------------------------------
         dv(i,j,k) = dt / (wk(i,j)+wk(i,j+1)) *   &
              ((gz(i,j,k+1)-gz(i,j+1,k))*(pk(i,j+1,k+1)-pk(i,j,k)) +  &
              (gz(i,j,k)-gz(i,j+1,k+1))*(pk(i,j,k+1)-pk(i,j+1,k)))
         !---------------------------------------------------------------------------------
         ! Non-hydrostatic contribution for half time-step
         v(i,j,k) = (v(i,j,k)  + divg2(i,j)-divg2(i,j+1) + alpha*dv(i,j,k) + &
              dt/(wk1(i,j)+wk1(i,j+1)) *  &
              ((gz(i,j,k+1)-gz(i,j+1,k))*(pp(i,j+1,k+1)-pp(i,j,k))   &
              + (gz(i,j,k)-gz(i,j+1,k+1))*(pp(i,j,k+1)-pp(i,j+1,k))))*rdy(i,j)
      enddo
   enddo

enddo    ! end k-loop


end subroutine split_p_grad



subroutine one_grad_p(u, v, pk, gz, divg2, delp, dt, ng, gridstruct, bd, npx, npy, npz,  &
   ptop, hydrostatic, a2b_ord, d_ext)  

integer, intent(IN) :: ng, npx, npy, npz, a2b_ord
real,    intent(IN) :: dt, ptop, d_ext
logical, intent(in) :: hydrostatic
type(fv_grid_bounds_type), intent(IN) :: bd
real,    intent(in) :: divg2(bd%is:bd%ie+1,bd%js:bd%je+1)
real, intent(inout) ::    pk(bd%isd:bd%ied,  bd%jsd:bd%jed  ,npz+1)
real, intent(inout) ::    gz(bd%isd:bd%ied,  bd%jsd:bd%jed  ,npz+1)
real, intent(inout) ::  delp(bd%isd:bd%ied,  bd%jsd:bd%jed  ,npz)
real, intent(inout) ::     u(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) 
real, intent(inout) ::     v(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz)
type(fv_grid_type), intent(INOUT), target :: gridstruct
! Local:
real, dimension(bd%isd:bd%ied,bd%jsd:bd%jed):: wk
real:: wk1(bd%is:bd%ie+1,bd%js:bd%je)
real:: wk2(bd%is:bd%ie,bd%js:bd%je+1)
real top_value
integer i,j,k

real, pointer, dimension(:,:) :: rdx,rdy

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

      rdx => gridstruct%rdx
      rdy => gridstruct%rdy

if ( hydrostatic ) then
   ! pk is pe**kappa if hydrostatic
   top_value = ptk
else
   ! pk is full pressure if non-hydrostatic
   top_value = ptop
endif

!$omp parallel do default(shared) 
do j=js,je+1
   do i=is,ie+1
      pk(i,j,1) = top_value
   enddo
enddo

!$omp parallel do default(shared) private(wk)
do k=2,npz+1
   if ( a2b_ord==4 ) then
      call a2b_ord4(pk(isd,jsd,k), wk, gridstruct, npx, npy, is, ie, js, je, ng, .true.)
   else
      call a2b_ord2(pk(isd,jsd,k), wk, gridstruct, npx, npy, is, ie, js, je, ng, .true.)
   endif
enddo

!$omp parallel do default(shared) private(wk)
do k=1,npz+1
   if ( a2b_ord==4 ) then
      call a2b_ord4( gz(isd,jsd,k), wk, gridstruct, npx, npy, is, ie, js, je, ng, .true.)
   else
      call a2b_ord2( gz(isd,jsd,k), wk, gridstruct, npx, npy, is, ie, js, je, ng, .true.)
   endif
enddo

if ( d_ext > 0. ) then

   !$omp parallel do default(shared) 
   do j=js,je+1
      do i=is,ie
         wk2(i,j) = divg2(i,j)-divg2(i+1,j)
      enddo
   enddo

   !$omp parallel do default(shared) 
   do j=js,je
      do i=is,ie+1
         wk1(i,j) = divg2(i,j)-divg2(i,j+1)
      enddo
   enddo

else

   !$omp parallel do default(shared) 
   do j=js,je+1
      do i=is,ie
         wk2(i,j) = 0.
      enddo
      do i=is,ie+1
         wk1(i,j) = 0.
      enddo
   enddo

endif

!$omp parallel do default(shared) private(wk)
do k=1,npz

   if ( hydrostatic ) then
      do j=js,je+1
         do i=is,ie+1
            wk(i,j) = pk(i,j,k+1) - pk(i,j,k)
         enddo
      enddo
   else
      if ( a2b_ord==4 ) then
         call a2b_ord4(delp(isd,jsd,k), wk, gridstruct, npx, npy, is, ie, js, je, ng)
      else
         call a2b_ord2(delp(isd,jsd,k), wk, gridstruct, npx, npy, is, ie, js, je, ng)
      endif
   endif

   do j=js,je+1
      do i=is,ie
         u(i,j,k) = rdx(i,j)*(wk2(i,j)+u(i,j,k) + dt/(wk(i,j)+wk(i+1,j)) * &
              ((gz(i,j,k+1)-gz(i+1,j,k))*(pk(i+1,j,k+1)-pk(i,j,k)) &
              + (gz(i,j,k)-gz(i+1,j,k+1))*(pk(i,j,k+1)-pk(i+1,j,k))))
      enddo
   enddo
   do j=js,je
      do i=is,ie+1
         v(i,j,k) = rdy(i,j)*(wk1(i,j)+v(i,j,k) + dt/(wk(i,j)+wk(i,j+1)) * &
              ((gz(i,j,k+1)-gz(i,j+1,k))*(pk(i,j+1,k+1)-pk(i,j,k)) &
              + (gz(i,j,k)-gz(i,j+1,k+1))*(pk(i,j,k+1)-pk(i,j+1,k))))
      enddo
   enddo
enddo    ! end k-loop

end subroutine one_grad_p



subroutine grad1_p_update(divg2, u, v, delu, delv, pk, gz, dt, ng, gridstruct, bd, npx, npy, npz, ptop, beta, a2b_ord)

integer, intent(in) :: ng, npx, npy, npz, a2b_ord
real,    intent(in) :: dt, ptop, beta
type(fv_grid_bounds_type), intent(IN) :: bd
real, intent(in):: divg2(bd%is:bd%ie+1,bd%js:bd%je+1)
real, intent(inout) ::    pk(bd%isd:bd%ied,  bd%jsd:bd%jed  ,npz+1)
real, intent(inout) ::    gz(bd%isd:bd%ied,  bd%jsd:bd%jed  ,npz+1)
real, intent(inout) ::     u(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) 
real, intent(inout) ::     v(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz)

real, intent(inout) ::    delu(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) 
real, intent(inout) ::    delv(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz)
type(fv_grid_type), intent(INOUT), target :: gridstruct

! Local:
real:: wk(bd%isd:bd%ied,bd%jsd:bd%jed)
real top_value, alpha
integer i,j,k
real, pointer, dimension(:,:) :: rdx,rdy

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd  = bd%isd
      ied  = bd%ied
      jsd  = bd%jsd
      jed  = bd%jed

      rdx => gridstruct%rdx
      rdy => gridstruct%rdy

alpha = 1. - beta

! pk is pe**kappa if hydrostatic
top_value = ptk

!$omp parallel do default(shared)
do j=js,je+1
   do i=is,ie+1
      pk(i,j,1) = top_value
   enddo
enddo
!$omp parallel do default(shared) private(wk)
do k=2,npz+1
   if ( a2b_ord==4 ) then
      call a2b_ord4(pk(isd,jsd,k), wk, gridstruct, npx, npy, is, ie, js, je, ng, .true.)
   else
      call a2b_ord2(pk(isd,jsd,k), wk, gridstruct, npx, npy, is, ie, js, je, ng, .true.)
   endif
enddo

!$omp parallel do default(shared) private(wk)
do k=1,npz+1
   if ( a2b_ord==4 ) then
      call a2b_ord4( gz(isd,jsd,k), wk, gridstruct, npx, npy, is, ie, js, je, ng, .true.)
   else
      call a2b_ord2( gz(isd,jsd,k), wk, gridstruct, npx, npy, is, ie, js, je, ng, .true.)
   endif
enddo

!$omp parallel do default(shared) private(wk)
do k=1,npz

   do j=js,je+1
      do i=is,ie+1
         wk(i,j) = pk(i,j,k+1) - pk(i,j,k)
      enddo
   enddo

   do j=js,je+1
      do i=is,ie
         u(i,j,k) = u(i,j,k) + beta*delu(i,j,k)
         delu(i,j,k) = dt/(wk(i,j)+wk(i+1,j)) *  &
              ((gz(i,j,k+1)-gz(i+1,j,k))*(pk(i+1,j,k+1)-pk(i,j,k)) &
              + (gz(i,j,k)-gz(i+1,j,k+1))*(pk(i,j,k+1)-pk(i+1,j,k)))
         u(i,j,k) = (u(i,j,k) + divg2(i,j)-divg2(i+1,j) + alpha*delu(i,j,k))*rdx(i,j)
      enddo
   enddo
   do j=js,je
      do i=is,ie+1
         v(i,j,k) = v(i,j,k) + beta*delv(i,j,k)
         delv(i,j,k) = dt/(wk(i,j)+wk(i,j+1)) *  &
              ((gz(i,j,k+1)-gz(i,j+1,k))*(pk(i,j+1,k+1)-pk(i,j,k)) &
              + (gz(i,j,k)-gz(i,j+1,k+1))*(pk(i,j,k+1)-pk(i,j+1,k)))
         v(i,j,k) = (v(i,j,k) + divg2(i,j)-divg2(i,j+1) + alpha*delv(i,j,k))*rdy(i,j)
      enddo
   enddo
enddo    ! end k-loop

end subroutine grad1_p_update


subroutine mix_wz(dt, w, delp, delz, pt, pkc, ng, km, zratio, last_step, grav, w_max, z_min, id_zratio, bd)
integer, intent(in):: km, ng, id_zratio
real, intent(in):: dt, grav, w_max, z_min
type(fv_grid_bounds_type), intent(IN) :: bd
real, intent(in):: pkc(bd%is-ng:bd%ie+ng, bd%js-ng:bd%je+ng, km+1)
real, intent(in), dimension(bd%isd:bd%ied,bd%jsd:bd%jed,km):: delp
real, intent(inout), dimension(bd%isd:bd%ied,bd%jsd:bd%jed,km):: pt
real, intent(inout) :: delz(bd%isd:bd%ied, bd%jsd:bd%jed, km)  ! delta-height (m)
real, intent(inout), dimension(bd%isd:bd%ied,bd%jsd:bd%jed,km):: w
real, intent(out):: zratio(bd%is:bd%ie, bd%js:bd%je, km)
logical, intent(in):: last_step
! Local:
real, parameter:: tau_w = 300.
real:: wk(bd%is:bd%ie,km)
real:: wm, tm, dz1, dz2, dz3, rat, relx
integer::ip(bd%js:bd%je)
integer i, j, k, ip_sum
logical used

      integer :: is,  ie,  js,  je

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je


relx = dt/tau_w

!$omp parallel do default(shared) private(wk, wm, tm, dz1, dz2, dz3, rat)
do j=js,je

        ip(j) = 0

#ifdef TEST_WZ

   do k=1,km
      do i=is,ie
         wk(i,k) = pt(i,j,k)*(pkc(i,j,k+1)-pkc(i,j,k))*rgrav
      enddo
   enddo

   ! Bottom: 
   k = km
   do i=is,ie
      dz2 = wk(i,k)
      zratio(i,j,k) = -delz(i,j,k)/dz2 
      if( abs(w(i,j,k)) > w_max .or. zratio(i,j,k)<z_min ) then

         !               tm = (wk(i,k-1)+wk(i,k))/(pkc(i,j,k+1)-pkc(i,j,k-1))*grav
         !               pt(i,j,k-1) = tm
         !               pt(i,j,k  ) = tm

         wm =  (w(i,j,k-1)*delp(i,j,k-1)+w(i,j,k)*delp(i,j,k))/(delp(i,j,k-1)+delp(i,j,k))
         w(i,j,k-1) = wm
         w(i,j,k  ) = wm
         dz1 = wk(i,k-1)
         rat = (delz(i,j,k-1) + delz(i,j,k)) / (dz1 + dz2)
         delz(i,j,k-1) = rat * dz1
         delz(i,j,k  ) = rat * dz2
         ip(j) = ip(j) + 1
      endif
   enddo

   do k=km-1,2,-1
      do i=is,ie
         dz2 = wk(i,k)
         zratio(i,j,k) = -delz(i,j,k)/dz2 
         if( abs(w(i,j,k)) > w_max .or. zratio(i,j,k)<z_min ) then

            !               tm = (wk(i,k-1)+wk(i,k)+wk(i,k+1))/(pkc(i,j,k+2)-pkc(i,j,k-1))*grav
            !               pt(i,j,k-1) = tm
            !               pt(i,j,k  ) = tm
            !               pt(i,j,k+1) = tm

            ! totally mix w
            wm = (w(i,j,k-1)*delp(i,j,k-1)+w(i,j,k)*delp(i,j,k)+w(i,j,k+1)*delp(i,j,k+1))   &
                 / (delp(i,j,k-1) + delp(i,j,k) + delp(i,j,k+1))
            w(i,j,k-1) = wm
            w(i,j,k  ) = wm
            w(i,j,k+1) = wm
            ! local hydrostatic adjustment of dz
            dz1 = wk(i,k-1)
            dz3 = wk(i,k+1)
            rat = (delz(i,j,k-1)+delz(i,j,k)+delz(i,j,k+1)) / (dz1 + dz2 + dz3)
            delz(i,j,k-1) = rat * dz1
            delz(i,j,k  ) = rat * dz2
            delz(i,j,k+1) = rat * dz3
            ip(j) = ip(j) + 1
         endif
      enddo
   enddo
#endif
   ! Last resort check
   do k=1,km
      do i=is,ie
         if( abs(w(i,j,k)) > w_max ) then
            ! Relax to w_max
            !                w(i,j,k) = (w(i,j,k)+ relx*sign(w_max, w(i,j,k)))/(1.+relx)
            w(i,j,k) = sign(w_max, w(i,j,k))
         endif
      enddo
   enddo
enddo   ! j-openMP loop

! ip_sum = sum ( ip(js:je) ) 
! if ( ip_sum>km*(je-js+1)*(ie-is+1)/100 ) write(*,*) 'Warning: Mix_wz for GID=', mpp_pe(), ' total points=', ip_sum 

if ( id_zratio>0 .and. last_step ) then
   zratio(is:ie,js:je,1) = 1.
   call prt_maxmin('DZ_ratio', zratio, is, ie, js, je, 0, km, 1.)
   used=send_data(id_zratio, zratio, fv_time)
endif

end subroutine  mix_wz


subroutine mix_dp(hydrostatic, w, delp, pt, km, ak, bk, CG, fv_debug, bd)
integer, intent(IN) :: km
real   , intent(IN) :: ak(km+1), bk(km+1)
type(fv_grid_bounds_type), intent(IN) :: bd
real, intent(INOUT), dimension(bd%isd:bd%ied,bd%jsd:bd%jed,km):: pt, delp, w
logical, intent(IN) :: hydrostatic, CG, fv_debug
! Local:
real dp, dpmin
integer i, j, k, ip
integer ifirst, ilast
integer jfirst, jlast

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


if ( CG ) then
   ifirst = is-1; ilast = ie+1
   jfirst = js-1; jlast = je+1
else
   ifirst = is; ilast = ie
   jfirst = js; jlast = je
endif


!$omp parallel do default(shared) private(ip, dpmin, dp)
do 1000 j=jfirst,jlast

   ip = 0

   do k=1, km-1
      dpmin = 0.01 * ( ak(k+1)-ak(k) + (bk(k+1)-bk(k))*1.E5 )
      do i=ifirst, ilast
         if(delp(i,j,k) < dpmin) then
            ! Remap from below and mix pt
            dp = dpmin - delp(i,j,k)
            pt(i,j,k) = (pt(i,j,k)*delp(i,j,k) + pt(i,j,k+1)*dp) / dpmin
            if ( .not.hydrostatic ) w(i,j,k) = (w(i,j,k)*delp(i,j,k) + w(i,j,k+1)*dp) / dpmin
            delp(i,j,k) = dpmin
            delp(i,j,k+1) = delp(i,j,k+1) - dp
            ip = ip + 1
         endif
      enddo
   enddo

   ! Bottom (k=km):
   dpmin = 0.01 * ( ak(km+1)-ak(km) + (bk(km+1)-bk(km))*1.E5 )
   do i=ifirst, ilast
      if(delp(i,j,km) < dpmin) then
         ! Remap from above and mix pt
         dp = dpmin - delp(i,j,km)
         pt(i,j,km) = (pt(i,j,km)*delp(i,j,km) + pt(i,j,km-1)*dp)/dpmin
         if ( .not.hydrostatic ) w(i,j,km) = (w(i,j,km)*delp(i,j,km) + w(i,j,km-1)*dp) / dpmin
         delp(i,j,km) = dpmin
         delp(i,j,km-1) = delp(i,j,km-1) - dp
         ip = ip + 1
      endif
   enddo
   if ( fv_debug .and. ip/=0 ) write(*,*) 'Warning: Mix_dp', mpp_pe(), j, ip 
   !      if ( ip/=0 ) write(*,*) 'Warning: Mix_dp', mpp_pe(), j, ip 
1000 continue

 end subroutine  mix_dp



 subroutine geopk(ptop, pe, peln, delp, pk, gz, hs, pt, pkz, km, akap, CG, nested, computehalo, npx, npy, a2b_ord, bd)

   integer, intent(IN) :: km, npx, npy, a2b_ord
   real   , intent(IN) :: akap, ptop
   type(fv_grid_bounds_type), intent(IN) :: bd
   real   , intent(IN) :: hs(bd%isd:bd%ied,bd%jsd:bd%jed)
   real, intent(IN), dimension(bd%isd:bd%ied,bd%jsd:bd%jed,km):: pt, delp
   logical, intent(IN) :: CG, nested, computehalo
   ! !OUTPUT PARAMETERS
   real, intent(OUT), dimension(bd%isd:bd%ied,bd%jsd:bd%jed,km+1):: gz, pk
   real, intent(OUT) :: pe(bd%is-1:bd%ie+1,km+1,bd%js-1:bd%je+1)
   real, intent(out) :: peln(bd%is:bd%ie,km+1,bd%js:bd%je)          ! ln(pe)
   real, intent(out) :: pkz(bd%is:bd%ie,bd%js:bd%je,km)
   ! !DESCRIPTION:
   !    Calculates geopotential and pressure to the kappa.
   ! Local:
   real p1d(bd%isd:bd%ied)
   real logp(bd%isd:bd%ied)
   real peln1
   integer i, j, k
   integer ifirst, ilast
   integer jfirst, jlast

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

#ifndef SW_DYNAMICS
   peln1 = log(ptop)
#endif

   if ( (.not. CG .and. a2b_ord==4) .or. (nested .and. .not. CG) ) then   ! D-Grid
      ifirst = is-2; ilast = ie+2
      jfirst = js-2; jlast = je+2
   else
      ifirst = is-1; ilast = ie+1
      jfirst = js-1; jlast = je+1
   endif

   if (nested .and. computehalo) then
      if (is == 1)     ifirst = isd
      if (ie == npx-1) ilast  = ied
      if (js == 1)     jfirst = jsd
      if (je == npy-1) jlast  = jed
   end if

   !$omp parallel do default(shared) private(p1d, logp)
   do 2000 j=jfirst,jlast

      do i=ifirst, ilast
         p1d(i) = ptop
         pk(i,j,1) = ptk
         gz(i,j,km+1) = hs(i,j)
      enddo

#ifndef SW_DYNAMICS
      if( j>=js .and. j<=je) then
         do i=is,ie
            peln(i,1,j) = peln1
         enddo
      endif
#endif

      if( j>(js-2) .and. j<(je+2) ) then
         do i=max(ifirst,is-1), min(ilast,ie+1) 
            pe(i,1,j) = ptop
         enddo
      endif

      ! Top down
      do k=2,km+1
         do i=ifirst, ilast
            p1d(i)  = p1d(i) + delp(i,j,k-1)
#ifdef GEOPK_CHECK
            if (p1d(i) < 1.e-10) then
               print*, 'NEGATIVE P1D: ', i, j, k
               print*, i, j, k, p1d(i), delp(i,j,k-1)
            end if
#endif
            !             pk(i,j,k) = p1d(i) ** akap
            ! Optimized form:
            logp(i) = log(p1d(i))
            pk(i,j,k) = exp( akap*logp(i) ) 
         enddo

         if( j>(js-2) .and. j<(je+2) ) then
            do i=max(ifirst,is-1), min(ilast,ie+1) 
               pe(i,k,j) = p1d(i)
            enddo
            if( j>=js .and. j<=je) then
               do i=is,ie
                  peln(i,k,j) = logp(i)
               enddo
            endif
         endif

      enddo

      ! Bottom up
      do k=km,1,-1
         do i=ifirst, ilast
            gz(i,j,k) = gz(i,j,k+1) + pt(i,j,k)*(pk(i,j,k+1)-pk(i,j,k))
         enddo
      enddo

      if ( .not. CG .and. j .ge. js .and. j .le. je ) then
         do k=1,km
            do i=is,ie
               pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k))/(akap*(peln(i,k+1,j)-peln(i,k,j)))
            enddo
         enddo
      endif

2000  continue
    end subroutine geopk

    subroutine grad2_scalar(grad_x, grad_y, qc, gridstruct, domain, npx, npy, ng, ghosted, bd)
      !
      ! Utilities Routine to compute gradients of a scalar on the native cubed-sphere grid using Green's theorem
      ! 2D version
      !
      integer, intent(in):: npx, npy, ng
      logical, intent(in):: ghosted
      type(fv_grid_bounds_type), intent(IN) :: bd
      ! Scalar at cell centers:
      real, intent(inout) ::  qc(bd%isd:bd%ied,bd%jsd:bd%jed)
      ! Spatial gradient in lat_lon directions
      real, intent(out), dimension(bd%is:bd%ie,bd%js:bd%je):: grad_x, grad_y
      type(fv_grid_type), intent(INOUT), target :: gridstruct
      type(domain2d), intent(INOUT) :: domain

      ! Local:
      real  qb(bd%isd:bd%ied,bd%jsd:bd%jed) 
      real grad3(3,bd%is:bd%ie,bd%js:bd%je)
      real   pdx(3,bd%is:bd%ie,bd%js:bd%je+1)
      real   pdy(3,bd%is:bd%ie+1,bd%js:bd%je)
      real, pointer, dimension(:,:) :: rarea, dx, dy
      real, pointer, dimension(:,:,:) :: vlon, vlat, en1, en2
      integer :: i,j, n

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      rarea => gridstruct%rarea
      dx    => gridstruct%dx
      dy    => gridstruct%dy

      vlon  => gridstruct%vlon
      vlat  => gridstruct%vlat
      en1   => gridstruct%en1
      en2   => gridstruct%en2

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

      if( .not. ghosted ) call mpp_update_domains( qc, domain, complete=.true. )

      ! Compute qb at 4 cell corners:

      call a2b_ord2(qc(isd,jsd), qb, gridstruct, npx, npy, is, ie, js, je, ng)

      !$omp parallel do default(shared)
      do j=js,je+1
         do i=is,ie
            do n=1,3
               pdx(n,i,j) = (qb(i,j)+qb(i+1,j))*dx(i,j)*en1(n,i,j)
            enddo
         enddo
      enddo
      !$omp parallel do default(shared)
      do j=js,je
         do i=is,ie+1
            do n=1,3
               pdy(n,i,j) = (qb(i,j)+qb(i,j+1))*dy(i,j)*en2(n,i,j)
            enddo
         enddo
      enddo

      ! Compute gradient by Green's theorem
      !$omp parallel do default(shared)
      do j=js,je
         do i=is,ie
            do n=1,3
               grad3(n,i,j) = pdx(n,i,j+1) - pdx(n,i,j) - pdy(n,i,j) + pdy(n,i+1,j)
            enddo
         enddo
      enddo

      !$omp parallel do default(shared)
      do j=js,je
         do i=is,ie
            grad_x(i,j) = 0.5*rarea(i,j)*(grad3(1,i,j)*vlon(i,j,1) + grad3(2,i,j)*vlon(i,j,2) + grad3(3,i,j)*vlon(i,j,3))
            grad_y(i,j) = 0.5*rarea(i,j)*(grad3(1,i,j)*vlat(i,j,1) + grad3(2,i,j)*vlat(i,j,2) + grad3(3,i,j)*vlat(i,j,3))
         enddo
      enddo

    end subroutine grad2_scalar


    subroutine grad_scalar(grad_x, grad_y, qc, gridstruct, domain, npx, npy, npz, ng, ghosted, bd)
      !
      ! Utilities Routine to compute gradients of a scalar on the native cubed-sphere grid using Green's theorem
      !
      integer, intent(in):: npx, npy, npz, ng
      logical, intent(in):: ghosted
      type(fv_grid_bounds_type), intent(IN) :: bd
      type(fv_grid_type), intent(IN), target :: gridstruct
      type(domain2d), intent(INOUT) :: domain
      ! Scalar at cell centers:
      real, intent(inout) ::  qc(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
      ! Spatial gradient in lat_lon directions
      real, intent(out), dimension(bd%is:bd%ie,bd%js:bd%je,npz):: grad_x, grad_y

      ! Local:
      real  qb(bd%isd:bd%ied,bd%jsd:bd%jed) 
      real grad3(3,bd%is:bd%ie,bd%js:bd%je)
      real   pdx(3,bd%is:bd%ie,bd%js:bd%je+1)
      real   pdy(3,bd%is:bd%ie+1,bd%js:bd%je)
      real, pointer, dimension(:,:) :: rarea, dx, dy
      real, pointer, dimension(:,:,:) :: vlon, vlat, en1, en2
      integer :: i,j,k, n

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      rarea => gridstruct%rarea
      dx    => gridstruct%dx
      dy    => gridstruct%dy

      vlon  => gridstruct%vlon
      vlat  => gridstruct%vlat
      en1   => gridstruct%en1
      en2   => gridstruct%en2

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed


      if( .not. ghosted ) call mpp_update_domains( qc, domain, complete=.true. )

      !$omp parallel do default(shared) private(pdx, pdy, qb, grad3)
      do k=1,npz

         ! Compute qb at 4 cell corners:

         call a2b_ord2(qc(isd,jsd,k), qb, gridstruct, npx, npy, is, ie, js, je, ng)

         do j=js,je+1
            do i=is,ie
               do n=1,3
                  pdx(n,i,j) = (qb(i,j)+qb(i+1,j))*dx(i,j)*en1(n,i,j)
               enddo
            enddo
         enddo
         do j=js,je
            do i=is,ie+1
               do n=1,3
                  pdy(n,i,j) = (qb(i,j)+qb(i,j+1))*dy(i,j)*en2(n,i,j)
               enddo
            enddo
         enddo

         ! Compute gradient by Green's theorem
         do j=js,je
            do i=is,ie
               do n=1,3
                  grad3(n,i,j) = pdx(n,i,j+1) - pdx(n,i,j) - pdy(n,i,j) + pdy(n,i+1,j)
               enddo
            enddo
         enddo

         do j=js,je
            do i=is,ie
               grad_x(i,j,k) = 0.5*rarea(i,j)*(grad3(1,i,j)*vlon(i,j,1) + grad3(2,i,j)*vlon(i,j,2) + grad3(3,i,j)*vlon(i,j,3))
               grad_y(i,j,k) = 0.5*rarea(i,j)*(grad3(1,i,j)*vlat(i,j,1) + grad3(2,i,j)*vlat(i,j,2) + grad3(3,i,j)*vlat(i,j,3))
            enddo
         enddo
      enddo

    end subroutine grad_scalar

    subroutine del2_cubed(q, cd, gridstruct, domain, npx, npy, km, nmax, bd)
      !---------------------------------------------------------------
      ! This routine is for filtering the omega field for the physics
      !---------------------------------------------------------------
      integer, intent(in):: npx, npy, km, nmax
      real,    intent(in):: cd            ! cd = K * da_min;   0 < K < 0.25
      type(fv_grid_bounds_type), intent(IN) :: bd
      real, intent(inout):: q(bd%isd:bd%ied,bd%jsd:bd%jed,km)
      type(fv_grid_type), intent(IN), target :: gridstruct
      type(domain2d), intent(INOUT) :: domain
      real, parameter:: r3  = 1./3.
      real :: fx(bd%isd:bd%ied+1,bd%jsd:bd%jed), fy(bd%isd:bd%ied,bd%jsd:bd%jed+1)
      real :: q2(bd%isd:bd%ied,bd%jsd:bd%jed)
      integer i,j,k, n, nt, ntimes

      logical :: do_nullify = .false.

      !Local routine pointers
      real, pointer, dimension(:,:) :: rarea
      real, pointer, dimension(:,:) :: sina_u, sina_v
      real, pointer, dimension(:,:) :: rdxc, rdyc, dx, dy
      logical, pointer :: sw_corner, se_corner, ne_corner, nw_corner

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

      rarea => gridstruct%rarea
      rdxc => gridstruct%rdxc
      rdyc => gridstruct%rdyc
      dx => gridstruct%dx
      dy => gridstruct%dy
      sina_u => gridstruct%sina_u
      sina_v => gridstruct%sina_v
      
      sw_corner => gridstruct%sw_corner
      nw_corner => gridstruct%nw_corner
      se_corner => gridstruct%se_corner
      ne_corner => gridstruct%ne_corner

      ntimes = min(3, nmax)

      call timing_on('COMM_TOTAL')
      call mpp_update_domains(q, domain, complete=.true.)
      call timing_off('COMM_TOTAL')


      do n=1,ntimes
         nt = ntimes - n

         !$omp parallel do default(shared) private(fx, fy)
         do k=1,km

            if ( sw_corner ) then
               q(1,1,k) = (q(1,1,k)+q(0,1,k)+q(1,0,k)) * r3
               q(0,1,k) =  q(1,1,k)
               q(1,0,k) =  q(1,1,k)
            endif
            if ( se_corner ) then
               q(ie, 1,k) = (q(ie,1,k)+q(npx,1,k)+q(ie,0,k)) * r3
               q(npx,1,k) =  q(ie,1,k)
               q(ie, 0,k) =  q(ie,1,k)
            endif
            if ( ne_corner ) then
               q(ie, je,k) = (q(ie,je,k)+q(npx,je,k)+q(ie,npy,k)) * r3
               q(npx,je,k) =  q(ie,je,k)
               q(ie,npy,k) =  q(ie,je,k)
            endif
            if ( nw_corner ) then
               q(1, je,k) = (q(1,je,k)+q(0,je,k)+q(1,npy,k)) * r3
               q(0, je,k) =  q(1,je,k)
               q(1,npy,k) =  q(1,je,k)
            endif

            if(nt>0) call copy_corners(q(isd,jsd,k), npx, npy, 1, gridstruct%nested, bd, &
                 sw_corner, se_corner, nw_corner, ne_corner )
            do j=js-nt,je+nt
               do i=is-nt,ie+1+nt
                  fx(i,j) = dy(i,j)*sina_u(i,j)*(q(i-1,j,k)-q(i,j,k))*rdxc(i,j)
               enddo
            enddo

            if(nt>0) call copy_corners(q(isd,jsd,k), npx, npy, 2, gridstruct%nested, bd, &
                 sw_corner, se_corner, nw_corner, ne_corner)
            do j=js-nt,je+1+nt
               do i=is-nt,ie+nt
                  fy(i,j) = dx(i,j)*sina_v(i,j)*(q(i,j-1,k)-q(i,j,k))*rdyc(i,j)
               enddo
            enddo

            do j=js-nt,je+nt
               do i=is-nt,ie+nt
                  q(i,j,k) = q(i,j,k) + cd*rarea(i,j)*(fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))
               enddo
            enddo
         enddo
      enddo

    end subroutine del2_cubed

    subroutine init_ijk_mem(i1, i2, j1, j2, km, array, var)
      integer, intent(in):: i1, i2, j1, j2, km
      real, intent(out):: array(i1:i2,j1:j2,km)
      real, intent(in):: var
      integer:: i, j, k

      !$omp parallel do default(shared)
      do k=1,km
         do j=j1,j2
            do i=i1,i2
               array(i,j,k) = var
            enddo
         enddo
      enddo

    end subroutine init_ijk_mem


  end module dyn_core_mod
