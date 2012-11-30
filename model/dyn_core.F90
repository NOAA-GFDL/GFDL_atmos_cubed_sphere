module dyn_core_mod

  use mpp_mod,            only: mpp_pe 
  use mpp_domains_mod,    only: CGRID_NE, DGRID_NE, mpp_get_boundary, mpp_update_domains,  &
                                mpp_start_update_domains, mpp_complete_update_domains
  use mpp_parameter_mod,  only: CORNER
  use fv_mp_mod,          only: domain, isd, ied, jsd, jed, is, ie, js, je, gid, square_domain
  use fv_control_mod,     only: hord_mt, hord_vt, hord_tm, hord_dp, hord_tr, n_sponge,  &
                                dddmp, d2_bg, d4_bg, d_ext, vtdm4, beta,        &
                                a2b_ord, master, fv_debug, d_con, scale_z, do_vort_damp, nord, &
                                fill_dp, nwat, inline_q, breed_vortex_inline,    &
                                d2_bg_k1, d2_bg_k2, d2_divg_max_k1, d2_divg_max_k2, damp_k_k1, damp_k_k2, &
                                m_split, w_max, z_min, fill_wz, convert_ke, use_old_omega
  use sw_core_mod,        only: c_sw, d_sw
  use a2b_edge_mod,       only: a2b_ord2, a2b_ord4
  use nh_core_mod,        only: Riem_Solver, Riem_Solver_C, update_dz_c, update_dz_d
  use fv_grid_tools_mod,  only: rdx, rdy, rdxc, rdyc, dx, dy, area, rarea, grid_type
  use fv_grid_utils_mod,  only: edge_vect_w, edge_vect_e, edge_vect_s, edge_vect_n,  &
                                ec1, ec2, en1, en2, vlon, vlat, da_min, da_min_c
  use fv_grid_utils_mod,  only: sina_u, sina_v, cosa_u, cosa_v,          &
                                sw_corner, se_corner, ne_corner, nw_corner
  use tp_core_mod,        only: copy_corners
  use fv_timing_mod,      only: timing_on, timing_off
  use fv_diagnostics_mod, only: prt_maxmin, idiag, fv_time
  use fv_nwp_nudge_mod,   only: breed_slp_inline
  use diag_manager_mod,   only: send_data
  use fv_mp_mod,          only: concurrent

  use boundary_mod,         only: extrapolation_BC,  nested_grid_BC_apply, nested_grid_BC_apply_intT
   use fv_current_grid_mod, only: nest_timestep, nsponge, nested, nestbctype, s_weight, k_split
   use fv_current_grid_mod, only: h_east, h_west, h_south, h_north
   use fv_current_grid_mod, only: q_east, q_west, q_south, q_north
#ifndef SW_DYNAMICS
   use fv_current_grid_mod, only: pt_east, pt_west, pt_south, pt_north
   use fv_current_grid_mod, only: w_east, w_west, w_south, w_north
#endif
   use fv_current_grid_mod, only: u_east, u_west, u_south, u_north
   use fv_current_grid_mod, only: v_east, v_west, v_south, v_north
   use fv_current_grid_mod, only: h_east_t0, h_west_t0, h_south_t0, h_north_t0
   use fv_current_grid_mod, only: q_east_t0, q_west_t0, q_south_t0, q_north_t0
#ifndef SW_DYNAMICS
   use fv_current_grid_mod, only: pt_east_t0, pt_west_t0, pt_south_t0, pt_north_t0
   use fv_current_grid_mod, only: w_east_t0, w_west_t0, w_south_t0, w_north_t0
#endif
   use fv_current_grid_mod, only: u_east_t0, u_west_t0, u_south_t0, u_north_t0
   use fv_current_grid_mod, only: v_east_t0, v_west_t0, v_south_t0, v_north_t0
   use fv_current_grid_mod, only: uc_east_t1, uc_west_t1, uc_south_t1, uc_north_t1
   use fv_current_grid_mod, only: vc_east_t1, vc_west_t1, vc_south_t1, vc_north_t1
   use fv_current_grid_mod, only: uc_east_t0, uc_west_t0, uc_south_t0, uc_north_t0
   use fv_current_grid_mod, only: vc_east_t0, vc_west_t0, vc_south_t0, vc_north_t0

#ifdef SW_DYNAMICS
  use test_cases_mod,      only: test_case, case9_forcing1, case9_forcing2
#endif

implicit none
private

public :: dyn_core, del2_cubed

  real :: ptk, rgrav
  real :: d3_damp
  real, allocatable, dimension(:,:,:) ::  ut, vt, crx, cry, xfx, yfx, divgd, &
                                          zh, du, dv, pkc, delpc, pk3, ptc, gz

!---- version number -----
  character(len=128) :: version = '$Id: dyn_core.F90,v 17.0.2.10.2.18 2012/05/10 20:41:48 Lucas.Harris Exp $'
  character(len=128) :: tagname = '$Name: siena_201211 $'

contains

!-----------------------------------------------------------------------
!     dyn_core :: FV Lagrangian dynamics driver
!-----------------------------------------------------------------------
 
 subroutine dyn_core(npx, npy, npz, ng, sphum, nq, bdt, n_split, zvir, cp, akap, grav, hydrostatic,  &
                     u,  v,  w, delz, pt, q, delp, pe, pk, phis, ws, omga, ptop, pfull, ua, va, & 
                     uc, vc, mfx, mfy, cx, cy, pkz, peln, ak, bk, init_step, i_pack, end_step, time_total)
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
    integer, intent(inout) :: i_pack(*)
    real, intent(inout), dimension(isd:ied  ,jsd:jed+1,npz):: u  ! D grid zonal wind (m/s)
    real, intent(inout), dimension(isd:ied+1,jsd:jed  ,npz):: v  ! D grid meridional wind (m/s)
    real, intent(inout) :: w(   isd:ied  ,jsd:jed  ,npz)  ! vertical vel. (m/s)
    real, intent(inout) :: delz(is :ie   ,js :je   ,npz)  ! delta-height (m)
    real, intent(inout) :: pt(  isd:ied  ,jsd:jed  ,npz)  ! temperature (K)
    real, intent(inout) :: delp(isd:ied  ,jsd:jed  ,npz)  ! pressure thickness (pascal)
    real, intent(inout) :: q(   isd:ied  ,jsd:jed  ,npz, nq)  ! 
    real, intent(in), optional:: time_total  ! total time (seconds) since start

!-----------------------------------------------------------------------
! Auxilliary pressure arrays:    
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
    real, intent(inout):: phis(isd:ied,jsd:jed)      ! Surface geopotential (g*Z_surf)
    real, intent(inout):: pe(is-1:ie+1, npz+1,js-1:je+1)  ! edge pressure (pascal)
    real, intent(inout):: peln(is:ie,npz+1,js:je)          ! ln(pe)
    real, intent(inout):: pk(is:ie,js:je, npz+1)        ! pe**kappa

!-----------------------------------------------------------------------
! Others:
    real,    parameter:: near0 = 1.E-8
    real,    parameter:: huge_r = 1.E40
    real,    parameter:: air_viscosity = 1.E-5   ! [m**2/sec] for T ~ 260 K
!-----------------------------------------------------------------------
    real, intent(out  ):: ws(is:ie,js:je)        ! w at surface
    real, intent(inout):: omga(isd:ied,jsd:jed,npz)    ! Vertical pressure velocity (pa/s)
    real, intent(inout):: uc(isd:ied+1,jsd:jed  ,npz)  ! (uc, vc) are mostly used as the C grid winds
    real, intent(inout):: vc(isd:ied  ,jsd:jed+1,npz)
    real, intent(inout), dimension(isd:ied,jsd:jed,npz):: ua, va

! The Flux capacitors: accumulated Mass flux arrays
    real, intent(inout)::  mfx(is:ie+1, js:je,   npz)
    real, intent(inout)::  mfy(is:ie  , js:je+1, npz)
! Accumulated Courant number arrays
    real, intent(inout)::  cx(is:ie+1, jsd:jed, npz)
    real, intent(inout)::  cy(isd:ied ,js:je+1, npz)
! Work:
    real, intent(inout):: pkz(is:ie,js:je,npz)  ! 

    real, allocatable, dimension(:,:,:):: pem, heat_source
! Auto 1D & 2D arrays:
    real:: dp_ref(npz)
    real:: zs(is:ie,js:je)        ! surface height (m)
    real:: p1d(is:ie)
    real:: om2d(is:ie,npz)
    real wbuffer(npy+2,npz)
    real ebuffer(npy+2,npz)
    real nbuffer(npx+2,npz)
    real sbuffer(npx+2,npz)
! ----   For external mode:
    real divg2(is:ie+1,js:je+1)
    real wk(isd:ied,jsd:jed)
    real fz(is: ie+1,js: je+1)
    real heat_s(is:ie,js:je)
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
    logical :: last_step
    logical used
    real :: split_timestep_bc

    ptk = ptop ** akap
    dt  = bdt / real(n_split)
    dt2 = 0.5*dt
    rdt = 1.0/dt
    rgrav = 1.0/grav

#ifdef SW_DYNAMICS
    k1k = 0.
#else
      k1k =  akap / (1.-akap)    ! rg/Cv=0.4
#endif
    kapag = -akap / grav

     diff_z0 = 0.25 * scale_z**2

     do k=1,npz
        dp_ref(k) = ak(k+1)-ak(k) + (bk(k+1)-bk(k))*1.E5  
     enddo

! Indexes:
      iep1 = ie + 1
      jep1 = je + 1
      if ( init_step ) then  ! Start of the big dynamic time stepping

           allocate(    gz(isd:ied, jsd:jed ,npz+1) )
           allocate(   pkc(isd:ied, jsd:jed ,npz+1) )
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
          if ( beta > near0 ) then
               allocate( du(isd:ied,  jsd:jed+1,npz) )
               call init_ijk_mem(isd,ied,   jsd,jed+1, npz, du, 0.)
               allocate( dv(isd:ied+1,jsd:jed,  npz) )
               call init_ijk_mem(isd,ied+1, jsd,jed  , npz, dv, 0.)
          endif
      endif    ! end init_step

      if ( beta < 1.e-5 ) beta = 0.

! Empty the "flux capacitors"
    call init_ijk_mem(is, ie+1, js,  je,   npz, mfx, 0.)
    call init_ijk_mem(is, ie  , js,  je+1, npz, mfy, 0.)
    call init_ijk_mem(is, ie+1, jsd, jed,  npz, cx, 0.)
    call init_ijk_mem(isd, ied, js,  je+1, npz, cy, 0.)

    if ( d_con > 1.0E-5 ) then
         allocate( heat_source(isd:ied, jsd:jed, npz) )
         call init_ijk_mem(isd, ied, jsd, jed, npz, heat_source, 0.)
    endif

    if ( .not. hydrostatic ) then
!$omp parallel do default(shared)
        do j=js,je
           do i=is,ie
              zs(i,j) = phis(i,j) * rgrav
           enddo
        enddo
    endif

    if ( convert_ke .or. vtdm4> 1.E-3 ) then
         n_con = npz
    else
         if ( d2_bg_k1 < 1.E-3 ) then
              n_con = 0
         else
              if ( d2_bg_k2 < 1.E-3 ) then
                   n_con = 1
              else
                   n_con = 2
              endif
         endif
    endif


     if (nested) then
        if (concurrent) then
           split_timestep_bc = real(n_split*k_split+nest_timestep)
        else
           split_timestep_bc = real(nest_timestep)
        endif
     endif

!-----------------------------------------------------
  do it=1,n_split
!-----------------------------------------------------

     if ( nq > 0 ) then
                                    call timing_on('COMM_TOTAL')
                                        call timing_on('COMM_TRACER')
         if ( inline_q ) then
                      i_pack(10) = mpp_start_update_domains(q, domain)
         elseif ( it==n_split ) then
                      i_pack(10) = mpp_start_update_domains(q, domain)
         endif
                                       call timing_off('COMM_TRACER')
                                   call timing_off('COMM_TOTAL')
     endif

     if ( .not. hydrostatic ) then
                             call timing_on('COMM_TOTAL')
         w_pack = mpp_start_update_domains(w, domain)
                             call timing_off('COMM_TOTAL')
!$omp parallel do default(shared)
         do j=js,je
            do i=is,ie
               gz(i,j,npz+1) = zs(i,j)
            enddo
            do k=npz,1,-1
               do i=is,ie
                  gz(i,j,k) = gz(i,j,k+1) - delz(i,j,k)
               enddo
            enddo
         enddo
                             call timing_on('COMM_TOTAL')
         i_pack(5) = mpp_start_update_domains(gz,  domain)
                             call timing_off('COMM_TOTAL')
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
          beta_d = beta
     endif

     if ( it==n_split .and. end_step ) then
       if ( use_old_omega ) then
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
     if( .not. hydrostatic ) call mpp_complete_update_domains(w_pack, w, domain)
                                                     call timing_off('COMM_TOTAL')

                                                     call timing_on('c_sw')
!$omp parallel do default(shared)
      do k=1,npz
         call c_sw(delpc(isd,jsd,k), delp(isd,jsd,k),  ptc(isd,jsd,k),    &
                      pt(isd,jsd,k),    u(isd,jsd,k),    v(isd,jsd,k),    &
                       w(isd,jsd,k),   uc(isd,jsd,k),   vc(isd,jsd,k),    &
                      ua(isd,jsd,k),   va(isd,jsd,k), omga(isd,jsd,k),    &
                      ut(isd,jsd,k),   vt(isd,jsd,k), divgd(isd,jsd,k),   &
                      nord,   dt2,  hydrostatic,  .true.)
      enddo
                                                     call timing_off('c_sw')
      if ( nord > 0 ) then
                                                   call timing_on('COMM_TOTAL')
          i_pack(3) = mpp_start_update_domains(divgd, domain, position=CORNER)
                                                  call timing_off('COMM_TOTAL')
      endif

!     if( fill_dp ) call mix_dp(hydrostatic, omga, delpc, ptc, npz, ak, bk, .true.)

      if ( hydrostatic ) then
           call geopk(ptop, pe, peln, delpc, pkc, gz, phis, ptc, pkz, npz, akap, .true., .false., npx, npy)
      else
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

                                            call timing_on('UPDATE_DZ_C')
           call update_dz_c(is, ie, js, je, npz, ng, dt2, dp_ref, zs, area, ut, vt, gz, ws)
                                            call timing_off('UPDATE_DZ_C')

                                               call timing_on('Riem_Solver_C')
           ms = max(1, m_split/2)
           call Riem_Solver_C( ms, dt2,   is,  ie,   js,   je,   npz,   ng,   &
                               akap,  cp,  ptop, phis, omga, ptc,  &
                               delpc, gz,  pkc, ws )
                                               call timing_off('Riem_Solver_C')

                                                                   call timing_on('COMM_TOTAL')
           if ( square_domain ) then
             i_pack(4) = mpp_start_update_domains(pkc, domain,whalo=1,ehalo=1,shalo=1, nhalo=1)
             i_pack(5) = mpp_start_update_domains(gz , domain,whalo=1,ehalo=1,shalo=1, nhalo=1)
             call mpp_complete_update_domains(i_pack(4), pkc,domain,whalo=1,ehalo=1,shalo=1, nhalo=1)
             call mpp_complete_update_domains(i_pack(5), gz ,domain,whalo=1,ehalo=1,shalo=1, nhalo=1)
!            call mpp_update_domains(pkc,domain,whalo=1,ehalo=1,shalo=1, nhalo=1, complete=.false.)
!            call mpp_update_domains(gz ,domain,whalo=1,ehalo=1,shalo=1, nhalo=1, complete=.true.)
           else   
             i_pack(4) = mpp_start_update_domains(pkc, domain)
             i_pack(5) = mpp_start_update_domains(gz,  domain)
             call mpp_complete_update_domains(i_pack(4), pkc, domain)
             call mpp_complete_update_domains(i_pack(5), gz , domain)
!            call mpp_update_domains(pkc, domain, complete=.false.)
!            call mpp_update_domains(gz , domain, complete=.true.)
           endif
                                                                   call timing_off('COMM_TOTAL')
      endif   ! end hydro check

      call p_grad_c(dt2, npz, delpc, pkc, gz, uc, vc, hydrostatic)

                                                                   call timing_on('COMM_TOTAL')
      i_pack(9) = mpp_start_update_domains(uc, vc, domain, gridtype=CGRID_NE)
                                                     call timing_off('COMM_TOTAL')
#ifdef SW_DYNAMICS
      if (test_case==9) call case9_forcing2(phis)
      endif !test_case>1
#endif

                                                                   call timing_on('COMM_TOTAL')
    if (inline_q .and. nq>0) call mpp_complete_update_domains(i_pack(10), q, domain)
    if ( nord > 0 )          call mpp_complete_update_domains(i_pack(3), divgd,  domain, position=CORNER)
                             call mpp_complete_update_domains(i_pack(9), uc, vc, domain, gridtype=CGRID_NE)

                                                                   call timing_off('COMM_TOTAL')
      if (nested) then
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
!!$            call extrapolation_BC(uc(:,:,k), 1, 0, npx, npy)
!!$            call extrapolation_BC(vc(:,:,k), 0, 1, npx, npy)
!!$         end do
!!$


         !vc
            call nested_grid_BC_apply_intT(vc, &
                 0, 1, npx, npy, npz, split_timestep_bc+0.5, real(n_split*k_split), & 
                 var_east_t0=vc_east_t0, &
                 var_west_t0=vc_west_t0, &
                 var_north_t0=vc_north_t0, &
                 var_south_t0=vc_south_t0, &
                 var_east_t1=vc_east_t1, &
                 var_west_t1=vc_west_t1, &
                 var_north_t1=vc_north_t1, &
                 var_south_t1=vc_south_t1, &
                 bctype=nestbctype, &
                 nsponge=nsponge, s_weight=s_weight  )

            !uc
            call nested_grid_BC_apply_intT(uc, &
                 1, 0, npx, npy, npz, split_timestep_bc+0.5, real(n_split*k_split), &
                 var_east_t0=uc_east_t0, &
                 var_west_t0=uc_west_t0, &
                 var_north_t0=uc_north_t0, &
                 var_south_t0=uc_south_t0, &
                 var_east_t1=uc_east_t1, &
                 var_west_t1=uc_west_t1, &
                 var_north_t1=uc_north_t1, &
                 var_south_t1=uc_south_t1, &
                 bctype=nestbctype, &
                 nsponge=nsponge, s_weight=s_weight   )

      end if

    if ( inline_q ) then
         if (nested) then
            do iq=1,nq
               if (concurrent) then
                  call nested_grid_BC_apply_intT(q(isd:ied,jsd:jed,:,iq), &
                       0, 0, npx, npy, npz, split_timestep_BC+1, real(n_split*k_split), &
                       var_east_t0=q_east_t0(:,:,:,iq), &
                       var_west_t0=q_west_t0(:,:,:,iq), &
                       var_north_t0=q_north_t0(:,:,:,iq), &
                       var_south_t0=q_south_t0(:,:,:,iq), &
                       var_east_t1=q_east(:,:,:,iq), &
                       var_west_t1=q_west(:,:,:,iq), &
                       var_north_t1=q_north(:,:,:,iq), &
                       var_south_t1=q_south(:,:,:,iq), &
                       bctype=nestbctype, &
                       nsponge=nsponge, s_weight=s_weight   )
               else
                  call nested_grid_BC_apply(q(isd:ied,jsd:jed,:,iq), &
                       0, 0, npx, npy, npz, nest_timestep+1, n_split*k_split, &
                       var_east=q_east(:,:,:,iq), &
                       var_west=q_west(:,:,:,iq), &
                       var_north=q_north(:,:,:,iq), &
                       var_south=q_south(:,:,:,iq), bctype=nestbctype, &
                       nsponge=nsponge, s_weight=s_weight   )
               endif
            end do
         end if
      endif

                                                     call timing_on('d_sw')
!$omp parallel do default(shared) private(nord_k, damp_k, d2_divg, dd_divg, hord_m, hord_v, hord_t, hord_p, wk, heat_s, d_con_k)
    do k=1,npz
       hord_m = hord_mt
       hord_t = hord_tm
       hord_v = hord_vt
       hord_p = hord_dp
       nord_k = nord
       nord_v(k) = min(2, nord)
       damp_k = dddmp
       d2_divg = min(0.20, d2_bg*(1.-3.*tanh(0.1*log(pfull(k)/pfull(npz)))))
       dd_divg = d4_bg
       d_con_k = d_con
       if ( do_vort_damp ) then
            damp_vt(k) = vtdm4
       else
            damp_vt(k) = 0.
       endif
       if ( npz==1 .or. n_sponge<0 ) then
           d2_divg = d2_bg
       elseif ( n_sponge==0 ) then
! New Del-2 Sponge layer: formulation
              if ( k==1 ) then
                   damp_k = 0.
                   nord_k = 0
                   d2_divg = max(0.21, d2_bg_k1, d2_bg, 0.05)
                   nord_v(k) = 0
                   damp_vt(k) = d2_bg_k1
              elseif ( k==2 ) then
                   damp_k = 0.
                   nord_k = 0
                   d2_divg = min(0.21, max(1.25*d2_bg_k2, d2_bg, 0.01))
                   nord_v(k) = 0
                   damp_vt(k) = d2_bg_k2
              elseif ( k==3 .and. d2_bg_k2>0.5 ) then
                   damp_k = 0.
                   nord_k = 0
                   d2_divg = max(0.2*d2_bg_k2, d2_bg, 0.01)
                   nord_v(k) = 0
                   damp_vt(k) = 0.2*d2_bg_k2
              endif

              if ( damp_vt(k) < 0.01 .and. nord_k>0 ) d_con_k = 0.
              if ( nord_v(k)==0 .and. damp_vt(k)>0.01 ) then
                   hord_t = 6
                   hord_v = 6
              endif
       else
           if( k <= n_sponge .and. npz>16 ) then
! Apply first order scheme for damping the sponge layer
               hord_m = 1
               hord_v = 1
               hord_t = 1
               hord_p = 1
               nord_k = 0
               damp_k = damp_k_k1
               d2_divg = min(0.20, d2_bg_k1*d2_bg)   ! 0.25 is the stability limit
               d2_divg = max(d2_divg_max_k1, d2_divg)
           elseif( k == n_sponge+1 .and. npz>24 ) then
               hord_v = 2
               hord_t = 2
               hord_p = 2
               nord_k = max(0, nord-1)
               d2_divg = min(0.20, d2_bg_k2*d2_bg)
               d2_divg = max(d2_divg_max_k2, d2_divg)
               if ( nord > 1 ) then
                    damp_k = 0.
               else
                    damp_k = damp_k_k2
               endif
           endif
       endif
       damp_k = max(damp_k, dddmp)

       if( .not. use_old_omega .and. last_step ) then
! Average horizontal "convergence" to cell center
            do j=js,je
               do i=is,ie
                  omga(i,j,k) = delp(i,j,k)
               enddo
            enddo
       endif

!--- external mode divergence damping ---
       if ( d_ext > 0. )  &
            call a2b_ord2(delp(isd,jsd,k), wk, npx, npy, is,    &
                          ie, js, je, ng, .false.)

       call d_sw(vt(isd,jsd,k), delp(isd,jsd,k), ptc(isd,jsd,k),  pt(isd,jsd,k),      &
                  u(isd,jsd,k),    v(isd,jsd,k),   w(isd,jsd,k),  uc(isd,jsd,k),      &
                  vc(isd,jsd,k),   ua(isd,jsd,k),  va(isd,jsd,k), divgd(isd,jsd,k),   &
                  mfx(is, js, k),  mfy(is, js, k),  cx(is, jsd,k),  cy(isd,js, k),    &
                  crx(is, jsd,k),  cry(isd,js, k), xfx(is, jsd,k), yfx(isd,js, k),    &
                  heat_s, zvir, sphum, nq,  q,  k,  npz, inline_q,  dt,               &
                  hord_tr, hord_m, hord_v, hord_t, hord_p, nord_k, nord_v(k), damp_k, &
                  d2_divg, dd_divg, damp_vt(k), d_con_k, hydrostatic)

       if( .not. use_old_omega .and. last_step ) then
! Average horizontal "convergence" to cell center
            do j=js,je
               do i=is,ie
                  omga(i,j,k) = omga(i,j,k)*(xfx(i,j,k)-xfx(i+1,j,k)+yfx(i,j,k)-yfx(i,j+1,k))*rarea(i,j)*rdt
               enddo
            enddo
       endif

       if ( d_ext > 0. ) then
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
    enddo         
                                                     call timing_off('d_sw')

    if( fill_dp ) call mix_dp(hydrostatic, w, delp, pt, npz, ak, bk, .false.)

                                                             call timing_on('COMM_TOTAL')
                                        i_pack(1) = mpp_start_update_domains(delp, domain)
    if ( hydrostatic .or. it/=n_split ) i_pack(2) = mpp_start_update_domains(pt,   domain)
!                                       i_pack(2) = mpp_start_update_domains(pt,   domain)
                                                             call timing_off('COMM_TOTAL')

    if ( d_ext > 0. ) then
          d2_divg = d_ext * da_min_c
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

     if ( hydrostatic ) then
                                       call timing_on('COMM_TOTAL')
          call mpp_complete_update_domains(i_pack(1), delp, domain)
          call mpp_complete_update_domains(i_pack(2), pt,   domain)
                                       call timing_off('COMM_TOTAL')
     else
                                            call timing_on('UPDATE_DZ')
        call update_dz_d(nord_v, damp_vt, hord_tm, is, ie, js, je, npz, ng, npx, npy, area,  &
                         dp_ref, zs, zh, crx, cry, xfx, yfx, delz, ws, rdt)
                                            call timing_off('UPDATE_DZ')

        if (idiag%id_ws>0 .and. last_step) then
            used=send_data(idiag%id_ws, ws, fv_time)
!           call prt_maxmin('WS', ws, is, ie, js, je, 0, 1, 1., master)
        endif
        if(fill_wz) call mix_wz(w, zh, delp, delz, pt, pk, 0, npz, pkz, last_step, grav)

                                       call timing_on('COMM_TOTAL')
        call mpp_complete_update_domains(i_pack(1), delp, domain)
                                       call timing_off('COMM_TOTAL')

     endif    ! end hydro check

    if (nested) then

       if (concurrent) then
          !!! Not sure whether to use +1 or not for the timestep for this variable.
          call nested_grid_BC_apply_intT(delp, &
               0, 0, npx, npy, npz, split_timestep_BC+1, real(n_split*k_split), &
               var_east_t0=h_east_t0, &
               var_west_t0=h_west_t0, &
               var_north_t0=h_north_t0, &
               var_south_t0=h_south_t0, &
               var_east_t1=h_east, &
               var_west_t1=h_west, &
               var_north_t1=h_north, &
               var_south_t1=h_south, &
               bctype=nestbctype, &
               nsponge=nsponge, s_weight=s_weight   )
       else
          call nested_grid_BC_apply(delp, &
               0, 0, npx, npy, npz, nest_timestep+1, n_split*k_split, &
               var_east=h_east, &
               var_west=h_west, &
               var_north=h_north, &
               var_south=h_south, bctype=nestbctype, &
               nsponge=nsponge, s_weight=s_weight  )
       endif


#ifndef SW_DYNAMICS

       if (concurrent) then
          call nested_grid_BC_apply_intT(pt, &
               0, 0, npx, npy, npz, split_timestep_BC+1, real(n_split*k_split), &
               var_east_t0=pt_east_t0, &
               var_west_t0=pt_west_t0, &
               var_north_t0=pt_north_t0, &
               var_south_t0=pt_south_t0, &
               var_east_t1=pt_east, &
               var_west_t1=pt_west, &
               var_north_t1=pt_north, &
               var_south_t1=pt_south, &
               bctype=nestbctype, &
               nsponge=nsponge, s_weight=s_weight   )
       else
          call nested_grid_BC_apply(pt, &
               0, 0, npx, npy, npz, nest_timestep+1, n_split*k_split, &
               var_east=pt_east, &
               var_west=pt_west, &
               var_north=pt_north, &
               var_south=pt_south, bctype=nestbctype, &
               nsponge=nsponge, s_weight=s_weight  )
       endif
#endif

    end if

     if ( hydrostatic ) then
          call geopk(ptop, pe, peln, delp, pkc, gz, phis, pt, pkz, npz, akap, .false., .true., npx, npy)
     else
!-----------------------------------------------------------
! pkc is non-hydrostatic perturbation pressure
!-----------------------------------------------------------
                                                         call timing_on('Riem_Solver')
        call Riem_Solver(m_split, dt,  is,  ie,   js,   je, npz,  ng,      &
                         akap, cp,  ptop, phis, w,  delz, pt, delp, zh,    &
                         gz,   pkc, pk, pe, peln, ws, da_min,    &
                         diff_z0,  it==n_split)
                                                         call timing_off('Riem_Solver')

                                                         call timing_on('COMM_TOTAL')
        i_pack(4) = mpp_start_update_domains(pkc, domain)
        i_pack(5) = mpp_start_update_domains(gz,  domain)
!$omp parallel do default(shared)
        do k=1,npz+1
           do j=js,je
              do i=is,ie
                 pk3(i,j,k) = pk(i,j,k)
              enddo
           enddo
        enddo
        i_pack(6) = mpp_start_update_domains(pk3, domain)
                                                         call timing_off('COMM_TOTAL')
     endif    ! end hydro check
     if (nested) then
                                       call timing_on('COMM_TOTAL')
        call mpp_update_domains(pkc, domain, complete=.false.)
        call mpp_update_domains(gz , domain, complete=.true.)
                                       call timing_off('COMM_TOTAL')
     end if

#ifdef SW_DYNAMICS
      if (test_case > 1) then !! Not sure whether to keep this check or not---does this code still work for case 1?
#else
      if ( hydrostatic .and. (breed_vortex_inline .or. it==n_split) ) then
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
    if ( .not. hydrostatic ) then
                                                                call timing_on('COMM_TOTAL')
         call mpp_complete_update_domains(i_pack(4), pkc, domain)
         call mpp_complete_update_domains(i_pack(5), gz,  domain)
         call mpp_complete_update_domains(i_pack(6), pk3, domain)
                                                                call timing_off('COMM_TOTAL')

       if ( beta > 0. ) then
         call split_p_grad(du, dv, u, v, pkc, gz, delp, pk3, divg2, beta_d, dt, ng, npx, npy, npz)
       else
         call nh_p_grad(u, v, pkc, gz, delp, pk3, divg2, dt, ng, npx, npy, npz)
       endif

    else
      if ( beta > 0. ) then
         call grad1_p_update(divg2, u, v, du, dv, pkc, gz, dt, ng, npx, npy, npz, ptop, beta_d)
      else
         call one_grad_p(u, v, pkc, gz, divg2, delp, dt, ng, npx, npy, npz, ptop, hydrostatic)
      endif
    endif
                                                                       call timing_on('COMM_TOTAL')
    if (.not. hydrostatic .and. it/=n_split) call mpp_complete_update_domains(i_pack(2), pt, domain)
!   if (.not. hydrostatic ) call mpp_complete_update_domains(i_pack(2), pt, domain)
                                                                       call timing_off('COMM_TOTAL')

!-------------------------------------------------------------------------------------------------------
    if ( breed_vortex_inline ) then
         call breed_slp_inline( it, dt, npz, ak, bk, phis, pe, pk, peln, pkz,     &
                                delp, u, v, pt, q, nwat, zvir )
    endif
!-------------------------------------------------------------------------------------------------------

                                                     call timing_on('COMM_TOTAL')
    if( it==n_split .and. grid_type<4 .and. .not. nested) then
! Prevent accumulation of rounding errors at overlapped domain edges:
          call mpp_get_boundary(u, v, domain, wbuffery=wbuffer, ebuffery=ebuffer,  &
                            sbufferx=sbuffer, nbufferx=nbuffer, gridtype=DGRID_NE )
#ifdef TEST_BUFF
          u(is:ie,je+1,1:npz) = nbuffer
          v(ie+1,js:je,1:npz) = ebuffer
#else
!$omp parallel do default(shared)
          do k=1,npz
             do i=is,ie
                u(i,je+1,k) = nbuffer(i-is+1,k)
             enddo
          enddo

!$omp parallel do default(shared)
          do k=1,npz
             do j=js,je
                v(ie+1,j,k) = ebuffer(j-js+1,k)
             enddo
          enddo
#endif
    endif

    if ( it/=n_split)   &
         i_pack(8) = mpp_start_update_domains(u, v, domain, gridtype=DGRID_NE)
                                                     call timing_off('COMM_TOTAL')

#ifdef SW_DYNAMICS
    endif
#endif
      if ( nested ) then
         nest_timestep = nest_timestep + 1
      endif

#ifdef SW_DYNAMICS
#else
    if ( last_step ) then
      if ( use_old_omega ) then
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
         call adv_pe(ua, va, pem, omga, npx, npy,  npz, ng)
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

    if (nested) then
#ifdef SW_DYNAMICS
#else



         if (.not. hydrostatic) then
            if (concurrent) then
               call nested_grid_BC_apply_intT(w, &
                    0, 0, npx, npy, npz, split_timestep_BC, real(n_split*k_split), &
                    var_east_t0=w_east_t0, &
                    var_west_t0=w_west_t0, &
                    var_north_t0=w_north_t0, &
                    var_south_t0=w_south_t0, &
                    var_east_t1=w_east, &
                    var_west_t1=w_west, &
                    var_north_t1=w_north, &
                    var_south_t1=w_south, &
                    bctype=nestbctype, &
                    nsponge=nsponge, s_weight=s_weight   )
            else
               call nested_grid_BC_apply(w, &
                    0, 0, npx, npy, npz, nest_timestep, n_split*k_split, &
                    var_east=w_east, &
                    var_west=w_west, &
                    var_north=w_north, &
                    var_south=w_south, bctype=nestbctype, &
                    nsponge=nsponge, s_weight=s_weight  )
            endif

         end if
#endif
         
         
         if (concurrent) then
            call nested_grid_BC_apply_intT(u, &
                 0, 1, npx, npy, npz, split_timestep_BC, real(n_split*k_split), &
                 var_east_t0=u_east_t0, &
                 var_west_t0=u_west_t0, &
                 var_north_t0=u_north_t0, &
                 var_south_t0=u_south_t0, &
                 var_east_t1=u_east, &
                 var_west_t1=u_west, &
                 var_north_t1=u_north, &
                 var_south_t1=u_south, &
                 bctype=nestbctype, &
                 nsponge=nsponge, s_weight=s_weight   )
         else
            call nested_grid_BC_apply(u, &
                 0, 1, npx, npy, npz, nest_timestep, n_split*k_split, &
                 var_east=u_east, &
                 var_west=u_west, &
                 var_north=u_north, &
                 var_south=u_south, bctype=nestbctype, &
                 nsponge=nsponge, s_weight=s_weight  )
         endif

         !v
         if (concurrent) then
            call nested_grid_BC_apply_intT(v, &
                 1, 0, npx, npy, npz, split_timestep_BC, real(n_split*k_split), &
                 var_east_t0=v_east_t0, &
                 var_west_t0=v_west_t0, &
                 var_north_t0=v_north_t0, &
                 var_south_t0=v_south_t0, &
                 var_east_t1=v_east, &
                 var_west_t1=v_west, &
                 var_north_t1=v_north, &
                 var_south_t1=v_south, &
                 bctype=nestbctype, &
                 nsponge=nsponge, s_weight=s_weight   )
         else
            call nested_grid_BC_apply(v, &
                 1, 0, npx, npy, npz, nest_timestep, n_split*k_split, &
                 var_east=v_east, &
                 var_west=v_west, &
                 var_north=v_north, &
                 var_south=v_south, bctype=nestbctype, &
                 nsponge=nsponge, s_weight=s_weight   )
         endif

      end if
!-----------------------------------------------------
  enddo   ! time split loop
!-----------------------------------------------------

  if ( fv_debug ) then
       if(master) write(*,*) 'End of n_split loop'
  endif


  if ( n_con/=0 .and. d_con > 1.e-5 ) then
       call del2_cubed(heat_source, 0.20*da_min, npx, npy, npz, 3)

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
    deallocate(   pkc )
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

  if ( fv_debug ) then
       if(master) write(*,*) 'End of dyn_core'
  endif

 end subroutine dyn_core


 subroutine adv_pe(ua, va, pem, om, npx, npy, npz, ng)

 integer, intent(in) :: npx, npy, npz, ng
! Contra-variant wind components:
 real, intent(in), dimension(isd:ied,jsd:jed,npz):: ua, va
! Pressure at edges:
 real, intent(in) :: pem(is-1:ie+1,1:npz+1,js-1:je+1)
 real, intent(inout) :: om(isd:ied,jsd:jed,npz)

! Local:
 real, dimension(is:ie,js:je):: up, vp
 real v3(3,is:ie,js:je)

 real pin(isd:ied,jsd:jed)
 real  pb(isd:ied,jsd:jed)

 real grad(3,is:ie,js:je)
 real pdx(3,is:ie,js:je+1)
 real pdy(3,is:ie+1,js:je)
 integer :: i,j,k, n

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
    call a2b_ord2(pin, pb, npx, npy, is, ie, js, je, ng)


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




 subroutine p_grad_c(dt2, npz, delpc, pkc, gz, uc, vc, hydrostatic)

 integer, intent(in):: npz
 real,    intent(in):: dt2
 real, intent(in), dimension(isd:ied, jsd:jed ,npz  ):: delpc
! pkc is pe**cappa     if hydrostatic
! pkc is full pressure if non-hydrostatic
 real, intent(in), dimension(isd:ied, jsd:jed ,npz+1):: pkc, gz
 real, intent(inout):: uc(isd:ied+1,jsd:jed  ,npz)
 real, intent(inout):: vc(isd:ied  ,jsd:jed+1,npz)
 logical, intent(in):: hydrostatic
! Local:
 real:: wk(is-1:ie+1,js-1:je+1)
 integer:: i,j,k

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


 subroutine nh_p_grad(u, v, pp, gh, delp, pk, divg2, dt, ng, npx, npy, npz)
    integer, intent(IN) :: ng, npx, npy, npz
    real,    intent(IN) :: dt
    real,    intent(in) :: divg2(is:ie+1, js:je+1)
    real, intent(inout) ::  delp(isd:ied, jsd:jed, npz)
    real, intent(inout) ::    pp(isd:ied, jsd:jed, npz+1)  ! perturbation pressure
    real, intent(inout) ::    pk(isd:ied, jsd:jed, npz+1)  ! p**kappa
    real, intent(inout) ::    gh(isd:ied, jsd:jed, npz+1)  ! g * h
    real, intent(inout) ::     u(isd:ied,  jsd:jed+1,npz) 
    real, intent(inout) ::     v(isd:ied+1,jsd:jed,  npz)
! Local:
    real wk1(isd:ied, jsd:jed)
    real  wk(is: ie+1,js: je+1)
    real du, dv
    integer i,j,k

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
           call a2b_ord4(pp(isd,jsd,k), wk1, npx, npy, is, ie, js, je, ng, .true.)
           call a2b_ord4(pk(isd,jsd,k), wk1, npx, npy, is, ie, js, je, ng, .true.)
       endif
       call a2b_ord4( gh(isd,jsd,k), wk1, npx, npy, is, ie, js, je, ng, .true.)
    enddo

!$omp parallel do default(shared) private(wk1, wk, du, dv)
    do k=1,npz
       call a2b_ord4(delp(isd,jsd,k), wk1, npx, npy, is, ie, js, je, ng)
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
                  ((gh(i,j,k+1)-gh(i+1,j,k))*(pk(i+1,j,k+1)-pk(i,j,k)) +  &
                  (gh(i,j,k)-gh(i+1,j,k+1))*(pk(i,j,k+1)-pk(i+1,j,k)))
! Non-hydrostatic contribution for half time-step
             u(i,j,k) = (u(i,j,k) + divg2(i,j)-divg2(i+1,j) + du +  &
                        dt/(wk1(i,j)+wk1(i+1,j)) *  &
                       ((gh(i,j,k+1)-gh(i+1,j,k))*(pp(i+1,j,k+1)-pp(i,j,k))    &
                      + (gh(i,j,k)-gh(i+1,j,k+1))*(pp(i,j,k+1)-pp(i+1,j,k))))*rdx(i,j)
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
! Current gradient from "hydrostatic" components:
             dv = dt / (wk(i,j)+wk(i,j+1)) *   &
                 ((gh(i,j,k+1)-gh(i,j+1,k))*(pk(i,j+1,k+1)-pk(i,j,k)) +  &
                 (gh(i,j,k)-gh(i,j+1,k+1))*(pk(i,j,k+1)-pk(i,j+1,k)))
! Non-hydrostatic contribution for half time-step
             v(i,j,k) = (v(i,j,k)  + divg2(i,j)-divg2(i,j+1) + dv + &
                        dt/(wk1(i,j)+wk1(i,j+1)) *  &
                       ((gh(i,j,k+1)-gh(i,j+1,k))*(pp(i,j+1,k+1)-pp(i,j,k))   &
                      + (gh(i,j,k)-gh(i,j+1,k+1))*(pp(i,j,k+1)-pp(i,j+1,k))))*rdy(i,j)
          enddo
       enddo

    enddo    ! end k-loop
 end subroutine nh_p_grad


 subroutine split_p_grad(du, dv, u, v, pp, gh, delp, pk, divg2, beta, dt, ng, npx, npy, npz)
    integer, intent(IN) :: ng, npx, npy, npz
    real,    intent(IN) :: beta, dt
    real,    intent(in) :: divg2(is:ie+1, js:je+1)
    real, intent(inout) ::  delp(isd:ied, jsd:jed, npz)
    real, intent(inout) ::    pp(isd:ied, jsd:jed, npz+1)  ! perturbation pressure
    real, intent(inout) ::    pk(isd:ied, jsd:jed, npz+1)  ! p**kappa
    real, intent(inout) ::    gh(isd:ied, jsd:jed, npz+1)  ! g * h
    real, intent(inout) ::    du(isd:ied  ,jsd:jed+1,npz) 
    real, intent(inout) ::    dv(isd:ied+1,jsd:jed  ,npz)
    real, intent(inout) ::     u(isd:ied,  jsd:jed+1,npz) 
    real, intent(inout) ::     v(isd:ied+1,jsd:jed,  npz)
! Local:
    real wk1(isd:ied, jsd:jed)
    real  wk(is: ie+1,js: je+1)
    real  alpha
    integer i,j,k

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
           call a2b_ord4(pp(isd,jsd,k), wk1, npx, npy, is, ie, js, je, ng, .true.)
           call a2b_ord4(pk(isd,jsd,k), wk1, npx, npy, is, ie, js, je, ng, .true.)
       endif
       call a2b_ord4( gh(isd,jsd,k), wk1, npx, npy, is, ie, js, je, ng, .true.)
    enddo

!$omp parallel do default(shared) private(wk1, wk)
    do k=1,npz
       call a2b_ord4(delp(isd,jsd,k), wk1, npx, npy, is, ie, js, je, ng)

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
                        ((gh(i,j,k+1)-gh(i+1,j,k))*(pk(i+1,j,k+1)-pk(i,j,k)) +  &
                         (gh(i,j,k)-gh(i+1,j,k+1))*(pk(i,j,k+1)-pk(i+1,j,k)))
!---------------------------------------------------------------------------------
! Non-hydrostatic contribution for half time-step
             u(i,j,k) = (u(i,j,k) + divg2(i,j)-divg2(i+1,j) + alpha*du(i,j,k) +  &
                        dt/(wk1(i,j)+wk1(i+1,j)) *  &
                       ((gh(i,j,k+1)-gh(i+1,j,k))*(pp(i+1,j,k+1)-pp(i,j,k))    &
                      + (gh(i,j,k)-gh(i+1,j,k+1))*(pp(i,j,k+1)-pp(i+1,j,k))))*rdx(i,j)
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             v(i,j,k) = v(i,j,k) + beta*dv(i,j,k)
! Current gradient from "hydrostatic" components:
!---------------------------------------------------------------------------------
             dv(i,j,k) = dt / (wk(i,j)+wk(i,j+1)) *   &
                        ((gh(i,j,k+1)-gh(i,j+1,k))*(pk(i,j+1,k+1)-pk(i,j,k)) +  &
                         (gh(i,j,k)-gh(i,j+1,k+1))*(pk(i,j,k+1)-pk(i,j+1,k)))
!---------------------------------------------------------------------------------
! Non-hydrostatic contribution for half time-step
             v(i,j,k) = (v(i,j,k)  + divg2(i,j)-divg2(i,j+1) + alpha*dv(i,j,k) + &
                        dt/(wk1(i,j)+wk1(i,j+1)) *  &
                       ((gh(i,j,k+1)-gh(i,j+1,k))*(pp(i,j+1,k+1)-pp(i,j,k))   &
                      + (gh(i,j,k)-gh(i,j+1,k+1))*(pp(i,j,k+1)-pp(i,j+1,k))))*rdy(i,j)
          enddo
       enddo

    enddo    ! end k-loop


 end subroutine split_p_grad



 subroutine one_grad_p(u, v, pk, gh, divg2, delp, dt, ng, npx, npy, npz,  &
                       ptop, hydrostatic)  

    integer, intent(IN) :: ng, npx, npy, npz
    real,    intent(IN) :: dt, ptop
    logical, intent(in) :: hydrostatic
    real,    intent(in) :: divg2(is:ie+1,js:je+1)
    real, intent(inout) ::    pk(isd:ied,  jsd:jed  ,npz+1)
    real, intent(inout) ::    gh(isd:ied,  jsd:jed  ,npz+1)
    real, intent(inout) ::  delp(isd:ied,  jsd:jed  ,npz)
    real, intent(inout) ::     u(isd:ied  ,jsd:jed+1,npz) 
    real, intent(inout) ::     v(isd:ied+1,jsd:jed  ,npz)
! Local:
    real, dimension(isd:ied,jsd:jed):: wk
    real:: wk1(is:ie+1,js:je)
    real:: wk2(is:ie,js:je+1)
    real top_value
    integer i,j,k

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
         call a2b_ord4(pk(isd,jsd,k), wk, npx, npy, is, ie, js, je, ng, .true.)
       else
         call a2b_ord2(pk(isd,jsd,k), wk, npx, npy, is, ie, js, je, ng, .true.)
       endif
    enddo

!$omp parallel do default(shared) private(wk)
    do k=1,npz+1
       if ( a2b_ord==4 ) then
         call a2b_ord4( gh(isd,jsd,k), wk, npx, npy, is, ie, js, je, ng, .true.)
       else
         call a2b_ord2( gh(isd,jsd,k), wk, npx, npy, is, ie, js, je, ng, .true.)
       endif
    enddo

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
            call a2b_ord4(delp(isd,jsd,k), wk, npx, npy, is, ie, js, je, ng)
         else
            call a2b_ord2(delp(isd,jsd,k), wk, npx, npy, is, ie, js, je, ng)
         endif
       endif

       do j=js,je+1
          do i=is,ie
             u(i,j,k) = rdx(i,j)*(wk2(i,j)+u(i,j,k) + dt/(wk(i,j)+wk(i+1,j)) * &
                        ((gh(i,j,k+1)-gh(i+1,j,k))*(pk(i+1,j,k+1)-pk(i,j,k)) &
                       + (gh(i,j,k)-gh(i+1,j,k+1))*(pk(i,j,k+1)-pk(i+1,j,k))))
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             v(i,j,k) = rdy(i,j)*(wk1(i,j)+v(i,j,k) + dt/(wk(i,j)+wk(i,j+1)) * &
                        ((gh(i,j,k+1)-gh(i,j+1,k))*(pk(i,j+1,k+1)-pk(i,j,k)) &
                       + (gh(i,j,k)-gh(i,j+1,k+1))*(pk(i,j,k+1)-pk(i,j+1,k))))
          enddo
       enddo
    enddo    ! end k-loop

 end subroutine one_grad_p



 subroutine grad1_p_update(divg2, u, v, delu, delv, pk, gh, dt, ng, npx, npy, npz, ptop, beta)

    integer, intent(in) :: ng, npx, npy, npz
    real,    intent(in) :: dt, ptop, beta
    real, intent(in):: divg2(is:ie+1,js:je+1)
    real, intent(inout) ::    pk(isd:ied,  jsd:jed  ,npz+1)
    real, intent(inout) ::    gh(isd:ied,  jsd:jed  ,npz+1)
    real, intent(inout) ::     u(isd:ied  ,jsd:jed+1,npz) 
    real, intent(inout) ::     v(isd:ied+1,jsd:jed  ,npz)

    real, intent(inout) ::    delu(isd:ied  ,jsd:jed+1,npz) 
    real, intent(inout) ::    delv(isd:ied+1,jsd:jed  ,npz)
! Local:
    real:: wk(isd:ied,jsd:jed)
    real top_value, alpha
    integer i,j,k

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
         call a2b_ord4(pk(isd,jsd,k), wk, npx, npy, is, ie, js, je, ng, .true.)
       else
         call a2b_ord2(pk(isd,jsd,k), wk, npx, npy, is, ie, js, je, ng, .true.)
       endif
    enddo

!$omp parallel do default(shared) private(wk)
    do k=1,npz+1
       if ( a2b_ord==4 ) then
         call a2b_ord4( gh(isd,jsd,k), wk, npx, npy, is, ie, js, je, ng, .true.)
       else
         call a2b_ord2( gh(isd,jsd,k), wk, npx, npy, is, ie, js, je, ng, .true.)
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
                         ((gh(i,j,k+1)-gh(i+1,j,k))*(pk(i+1,j,k+1)-pk(i,j,k)) &
                      + (gh(i,j,k)-gh(i+1,j,k+1))*(pk(i,j,k+1)-pk(i+1,j,k)))
             u(i,j,k) = (u(i,j,k) + divg2(i,j)-divg2(i+1,j) + alpha*delu(i,j,k))*rdx(i,j)
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             v(i,j,k) = v(i,j,k) + beta*delv(i,j,k)
             delv(i,j,k) = dt/(wk(i,j)+wk(i,j+1)) *  &
                         ((gh(i,j,k+1)-gh(i,j+1,k))*(pk(i,j+1,k+1)-pk(i,j,k)) &
                      + (gh(i,j,k)-gh(i,j+1,k+1))*(pk(i,j,k+1)-pk(i,j+1,k)))
             v(i,j,k) = (v(i,j,k) + divg2(i,j)-divg2(i,j+1) + alpha*delv(i,j,k))*rdy(i,j)
          enddo
       enddo
    enddo    ! end k-loop

 end subroutine grad1_p_update


 subroutine mix_wz(w, zh, delp, delz, pt, pkc, ng, km, zratio, last_step, grav)
  integer, intent(in):: km, ng
  real, intent(in):: grav
  real, intent(in):: pkc(is-ng:ie+ng, js-ng:je+ng, km+1)
  real, intent(in), dimension(isd:ied,jsd:jed,km):: delp
  real, intent(inout), dimension(isd:ied,jsd:jed,km+1):: zh
  real, intent(inout), dimension(isd:ied,jsd:jed,km):: pt
  real, intent(inout) :: delz(is:ie, js:je, km)  ! delta-height (m)
  real, intent(inout), dimension(isd:ied,jsd:jed,km):: w
  real, intent(out):: zratio(is:ie, js:je, km)
  logical, intent(in):: last_step
! Local:
  real:: wk(is:ie,km)
  real:: wm, tm, dz1, dz2
  integer::ip(js:je)
  integer i, j, k, ip_sum
  logical used


!$omp parallel do default(shared) private(wk, wm, tm, dz1, dz2)
     do j=js,je

        do k=1,km
           do i=is,ie
              wk(i,k) = pt(i,j,k)*(pkc(i,j,k+1)-pkc(i,j,k))*rgrav
           enddo
        enddo

        ip(j) = 0
        do k=km,2,-1
          do i=is,ie
             dz2 = wk(i,k)
             zratio(i,j,k) = (zh(i,j,k)-zh(i,j,k+1))/dz2 
             if( abs(w(i,j,k)) > w_max .or. zratio(i,j,k)<z_min ) then
                tm = (wk(i,k-1)+wk(i,k))/(pkc(i,j,k+1)-pkc(i,j,k-1))*grav
                pt(i,j,k-1) = tm
                pt(i,j,k  ) = tm
                wm =  (w(i,j,k-1)*delp(i,j,k-1)+w(i,j,k)*delp(i,j,k))/(delp(i,j,k-1)+delp(i,j,k))
                w(i,j,k-1) = wm
                w(i,j,k  ) = wm
                dz1 = wk(i,k-1)
                zh(i,j,k) = zh(i,j,k+1) + (zh(i,j,k-1)-zh(i,j,k+1))*dz2/(dz1+dz2)
                ip(j) = ip(j) + 1
             endif
          enddo
        enddo
     enddo

! ip_sum = sum ( ip(js:je) ) 
! if ( ip_sum>km*(je-js+1)*(ie-is+1)/100 ) write(*,*) 'Warning: Mix_wz for GID=', gid, ' total points=', ip_sum 

  if ( idiag%id_zratio>0 .and. last_step ) then
       zratio(is:ie,js:je,1) = 1.
       call prt_maxmin('DZ_ratio', zratio, is, ie, js, je, 0, km, 1., master)
       used=send_data(idiag%id_zratio, zratio, fv_time)
  endif

 end subroutine  mix_wz


 subroutine mix_dp(hydrostatic, w, delp, pt, km, ak, bk, CG)
  integer, intent(IN) :: km
  real   , intent(IN) :: ak(km+1), bk(km+1)
  real, intent(INOUT), dimension(isd:ied,jsd:jed,km):: pt, delp, w
  logical, intent(IN) :: hydrostatic, CG
! Local:
     real dp, dpmin
     integer i, j, k, ip
     integer ifirst, ilast
     integer jfirst, jlast


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
       if ( fv_debug .and. ip/=0 ) write(*,*) 'Warning: Mix_dp', gid, j, ip 
!      if ( ip/=0 ) write(*,*) 'Warning: Mix_dp', gid, j, ip 
1000   continue

 end subroutine  mix_dp



 subroutine geopk(ptop, pe, peln, delp, pk, gh, hs, pt, pkz, km, akap, CG, computehalo, npx, npy)

     integer, intent(IN) :: km, npx, npy
     real   , intent(IN) :: akap, ptop
     real   , intent(IN) :: hs(isd:ied,jsd:jed)
     real, intent(IN), dimension(isd:ied,jsd:jed,km):: pt, delp
     logical, intent(IN) :: CG, computehalo
! !OUTPUT PARAMETERS
     real, intent(OUT), dimension(isd:ied,jsd:jed,km+1):: gh, pk
     real, intent(OUT) :: pe(is-1:ie+1,km+1,js-1:je+1)
     real, intent(out) :: peln(is:ie,km+1,js:je)          ! ln(pe)
     real, intent(out) :: pkz(is:ie,js:je,km)
! !DESCRIPTION:
!    Calculates geopotential and pressure to the kappa.
! Local:
     real p1d(isd:ied)
     real logp(isd:ied)
     real peln1
     integer i, j, k
     integer ifirst, ilast
     integer jfirst, jlast

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
           gh(i,j,km+1) = hs(i,j)
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
              gh(i,j,k) = gh(i,j,k+1) + pt(i,j,k)*(pk(i,j,k+1)-pk(i,j,k))
           enddo
        enddo

      if ( .not. CG .and. j .ge. js .and. j .le. je ) then
! This is for hydrostatic only
         do k=1,km
            do i=is,ie
               pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k))/(akap*(peln(i,k+1,j)-peln(i,k,j)))
            enddo
         enddo
      endif

2000  continue
 end subroutine geopk



 subroutine grad2_scalar(grad_x, grad_y, qc, npx, npy, ng, ghosted)
!
! Utilities Routine to compute gradients of a scalar on the native cubed-sphere grid using Green's theorem
! 2D version
!
 integer, intent(in):: npx, npy, ng
 logical, intent(in):: ghosted
! Scalar at cell centers:
 real, intent(inout) ::  qc(isd:ied,jsd:jed)
! Spatial gradient in lat_lon directions
 real, intent(out), dimension(is:ie,js:je):: grad_x, grad_y

! Local:
 real  qb(isd:ied,jsd:jed) 
 real grad3(3,is:ie,js:je)
 real   pdx(3,is:ie,js:je+1)
 real   pdy(3,is:ie+1,js:je)
 integer :: i,j, n

 if( .not. ghosted ) call mpp_update_domains( qc, domain, complete=.true. )

! Compute qb at 4 cell corners:

    call a2b_ord2(qc(isd,jsd), qb, npx, npy, is, ie, js, je, ng)

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

 
 subroutine grad_scalar(grad_x, grad_y, qc, npx, npy, npz, ng, ghosted)
!
! Utilities Routine to compute gradients of a scalar on the native cubed-sphere grid using Green's theorem
!
 integer, intent(in):: npx, npy, npz, ng
 logical, intent(in):: ghosted
! Scalar at cell centers:
 real, intent(inout) ::  qc(isd:ied,jsd:jed,npz)
! Spatial gradient in lat_lon directions
 real, intent(out), dimension(is:ie,js:je,npz):: grad_x, grad_y

! Local:
 real  qb(isd:ied,jsd:jed) 
 real grad3(3,is:ie,js:je)
 real   pdx(3,is:ie,js:je+1)
 real   pdy(3,is:ie+1,js:je)
 integer :: i,j,k, n


 if( .not. ghosted ) call mpp_update_domains( qc, domain, complete=.true. )

!$omp parallel do default(shared) private(pdx, pdy, qb, grad3)
 do k=1,npz

! Compute qb at 4 cell corners:
  
    call a2b_ord2(qc(isd,jsd,k), qb, npx, npy, is, ie, js, je, ng)

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

 subroutine del2_cubed(q, cd, npx, npy, km, nmax)
!---------------------------------------------------------------
! This routine is for filtering the omega field for the physics
!---------------------------------------------------------------
   integer, intent(in):: npx, npy, km, nmax
   real,    intent(in):: cd            ! cd = K * da_min;   0 < K < 0.25
   real, intent(inout):: q(isd:ied,jsd:jed,km)
   real, parameter:: r3  = 1./3.
   real :: fx(isd:ied+1,jsd:jed), fy(isd:ied,jsd:jed+1)
   real :: q2(isd:ied,jsd:jed)
   integer i,j,k, n, nt, ntimes

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

      if(nt>0) call copy_corners(q(isd,jsd,k), npx, npy, 1)
      do j=js-nt,je+nt
         do i=is-nt,ie+1+nt
            fx(i,j) = dy(i,j)*sina_u(i,j)*(q(i-1,j,k)-q(i,j,k))*rdxc(i,j)
         enddo
      enddo

      if(nt>0) call copy_corners(q(isd,jsd,k), npx, npy, 2)
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
