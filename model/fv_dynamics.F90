module fv_dynamics_mod
   use constants_mod,   only: grav, pi, radius, hlv    ! latent heat of water vapor
   use dyn_core_mod,    only: dyn_core
   use fv_mapz_mod,     only: compute_total_energy, Lagrangian_to_Eulerian
   use fv_tracer2d_mod,   only: tracer_2d, tracer_2d_1L
   use fv_control_mod,  only: hord_mt, hord_vt, hord_tm, hord_tr, &
                              kord_mt, kord_tm, kord_tr, full_phys, &
                              z_tracer, tau, rf_center, nf_omega,   &
                              uniform_ppm, remap_t,  k_top, p_ref,  &
                              nwat, fv_debug
   use fv_grid_utils_mod,      only: sina_u, sina_v, sw_corner, se_corner, &
                              ne_corner, nw_corner, da_min, ptop,   &
                              cubed_to_latlon, c2l_ord2
   use fv_grid_tools_mod,      only: dx, dy, rdxa, rdya, rdxc, rdyc, area, rarea
   use fv_mp_mod,          only: is,js,ie,je, isd,jsd,ied,jed, gid, domain
   use fv_timing_mod,    only: timing_on, timing_off

   use diag_manager_mod,    only: send_data
   use fv_diagnostics_mod,  only: id_divg, id_te, fv_time, prt_maxmin
   use mpp_domains_mod, only: mpp_update_domains, DGRID_NE
   use field_manager_mod,  only: MODEL_ATMOS
   use tracer_manager_mod, only: get_tracer_index
   use fv_sg_mod,          only: neg_adj3

implicit none
   logical :: RF_initialized = .false.
   real, allocatable ::  rf(:), rw(:)
private
public :: fv_dynamics


contains

!-----------------------------------------------------------------------
!     fv_dynamics :: FV dynamical core driver
!-----------------------------------------------------------------------
 
  subroutine fv_dynamics(npx, npy, npz, nq,  ng, bdt, consv_te, fill,               &
                        reproduce_sum, kappa, cp_air, zvir, ks, ncnst, n_split,     &
                        q_split, u, v, w, delz, hydrostatic, pt, delp, q,           &
                        ps, pe, pk, peln, pkz, phis, omga, ua, va, uc, vc,          &
                        ak, bk, mfx, mfy, cx, cy, ze0, hybrid_z, time_total)

    real, intent(IN) :: bdt  ! Large time-step
    real, intent(IN) :: consv_te
    real, intent(IN) :: kappa, cp_air
    real, intent(IN) :: zvir
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

    real, intent(inout) :: u(   isd:ied  ,jsd:jed+1,npz)  ! D grid zonal wind (m/s)
    real, intent(inout) :: v(   isd:ied+1,jsd:jed  ,npz)  ! D grid meridional wind (m/s)
    real, intent(inout) :: w(   isd:ied  ,jsd:jed  ,npz)  !  W (m/s)
    real, intent(inout) :: pt(  isd:ied  ,jsd:jed  ,npz)  ! temperature (K)
    real, intent(inout) :: delp(isd:ied  ,jsd:jed  ,npz)  ! pressure thickness (pascal)
    real, intent(inout) :: q(   isd:ied  ,jsd:jed  ,npz, ncnst) ! specific humidity and constituents
    real, intent(inout) :: delz(is:ie,js:je,npz)   ! delta-height (m); non-hydrostatic only
    real, intent(inout) ::  ze0(is:ie,js:je,npz+1) ! height at edges (m); non-hydrostatic

!-----------------------------------------------------------------------
! Auxilliary pressure arrays:    
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
    real, intent(inout) :: ps  (isd:ied  ,jsd:jed)           ! Surface pressure (pascal)
    real, intent(inout) :: pe  (is-1:ie+1, npz+1,js-1:je+1)  ! edge pressure (pascal)
    real, intent(inout) :: pk  (is:ie,js:je, npz+1)          ! pe**cappa
    real, intent(inout) :: peln(is:ie,npz+1,js:je)           ! ln(pe)
    real, intent(inout) :: pkz (is:ie,js:je,npz)             ! finite-volume mean pk
    
!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
    real, intent(inout) :: phis(isd:ied,jsd:jed)       ! Surface geopotential (g*Z_surf)
    real, intent(inout) :: omga(isd:ied,jsd:jed,npz)   ! Vertical pressure velocity (pa/s)
    real, intent(inout) :: uc(isd:ied+1,jsd:jed  ,npz) ! (uc,vc) mostly used as the C grid winds
    real, intent(inout) :: vc(isd:ied  ,jsd:jed+1,npz)

    real, intent(inout), dimension(isd:ied ,jsd:jed ,npz):: ua, va
    real, intent(in),    dimension(npz+1):: ak, bk

! Accumulated Mass flux arrays: the "Flux Capacitor"
    real, intent(inout) ::  mfx(is:ie+1, js:je,   npz)
    real, intent(inout) ::  mfy(is:ie  , js:je+1, npz)
! Accumulated Courant number arrays
    real, intent(inout) ::  cx(is:ie+1, jsd:jed, npz)
    real, intent(inout) ::  cy(isd:ied ,js:je+1, npz)


! Local Arrays
      real:: q2(isd:ied,jsd:jed,nq)
      real:: te_2d(is:ie,js:je)
      real::   teq(is:ie,js:je)
      real:: pfull(npz)
      real, allocatable :: dp1(:,:,:)
      real, allocatable :: pem(:,:,:)
      real:: akap, rg, ph1, ph2
      integer :: i,j,k, iq
      integer :: sphum, liq_wat, ice_wat      ! GFDL physics
      integer :: rainwat, snowwat, graupel
      logical used
      real te_den

      allocate ( dp1(is:ie, js:je, 1:npz) )
      allocate ( pem(is-1:ie+1, 1:npz+1, js-1:je+1) )

#ifdef SW_DYNAMICS
      akap  = 1.
                                                  call timing_on('COMM_TOTAL')
      call mpp_update_domains(u, v, domain, gridtype=DGRID_NE, complete=.true.)
                                                 call timing_off('COMM_TOTAL')
#else
      if ( nwat==6 ) then
             sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
           liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
           ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
           rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat')
           snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat')
           graupel = get_tracer_index (MODEL_ATMOS, 'graupel')
      else
           sphum = 1
      endif

      akap  = kappa
      rg = kappa*cp_air

      do k=1,npz
         ph1 = ak(k  ) + bk(k  )*p_ref
         ph2 = ak(k+1) + bk(k+1)*p_ref
         pfull(k) = (ph2 - ph1) / log(ph2/ph1)
      enddo

      if( tau > 0. )      &
      call Rayleigh_Friction(bdt, npz, ks, pfull, tau, rf_center, u, v, w, pt,  &
                             ua, va, cp_air, rg,  hydrostatic, .true.)

                                                  call timing_on('COMM_TOTAL')
      call mpp_update_domains(u, v, domain, gridtype=DGRID_NE, complete=.true.)
                                                 call timing_off('COMM_TOTAL')
!---------------------
! Compute Total Energy
!---------------------
      if ( consv_te > 0. ) then
           call compute_total_energy(is, ie, js, je, isd, ied, jsd, jed, npz,  &
                                     u, v, pt, delp, q, pe, peln, phis, zvir,  &
                                     cp_air, rg, hlv, te_2d, ua, va, teq,      &
                                     full_phys, sphum, id_te)
           if( id_te>0 ) then
               used = send_data(id_te, teq, fv_time)
!              te_den=1.E-9*g_sum(teq, is, ie, js, je, ng, area, 0)/(grav*4.*pi*radius**2)
!              if(gid==0)  write(*,*) 'Total Energy Density (Giga J/m**2)=',te_den
           endif
      endif

 
! Convert pt to virtual potential temperature * CP
!$omp parallel do private (i, j, k)
      do k=1,npz
         do j=js,je
            do i=is,ie
                pt(i,j,k) = cp_air*pt(i,j,k)/pkz(i,j,k)*(1.+zvir*q(i,j,k,sphum))
            enddo
         enddo
      enddo
#endif

!$omp parallel do private (i, j, k)
      do k=1,npz
         do j=js,je
            do i=is,ie
               dp1(i,j,k) = delp(i,j,k)
            enddo
         enddo
      enddo

      call dyn_core(npx, npy, npz, ng, bdt, n_split, cp_air, akap, grav, hydrostatic, &
                    u, v, w, delz, pt, delp, pe, pk, phis, omga, ptop, pfull, ua, va, & 
                    uc, vc, mfx, mfy, cx, cy, pem, pkz, uniform_ppm, time_total)

#ifdef SW_DYNAMICS
      do j=js,je
         do i=is,ie
            ps(i,j) = delp(i,j,1) / grav
         enddo
      enddo
#else
      if(nq /= 0) then    
!--------------------------------------------------------
! Perform large-time-step scalar transport using the accumulated CFL and
! mass fluxes
         call timing_on('tracer_2d')
       if ( z_tracer ) then
         do k=1,npz
            do iq=1,nq
            do j=js,je
               do i=is,ie                   ! To_do list:
                  q2(i,j,iq) = q(i,j,k,iq)  ! The data copying can be avoided if q is
                                            ! re-dimensioned as q(i,j,nq,k)
               enddo
            enddo
            enddo
         call tracer_2d_1L(q2, dp1(is,js,k), mfx(is,js,k), mfy(is,js,k), &
                           cx(is,jsd,k),  cy(isd,js,k), npx, npy, npz,   &
                           nq, hord_tr, q_split, k, q, bdt, uniform_ppm, id_divg)
         enddo
       else
         call tracer_2d(q, dp1, mfx, mfy, cx, cy, npx, npy, npz, nq, &
                        hord_tr, q_split, bdt, uniform_ppm, id_divg)
       endif
         call timing_off('tracer_2d')
         if( id_divg>0 ) used = send_data(id_divg, dp1, fv_time) 
      endif


      if ( npz > 4 ) then
!------------------------------------------------------------------------
! Peroform vertical remapping from Lagrangian control-volume to
! the Eulerian coordinate as specified by the routine set_eta.
! Note that this finite-volume dycore is otherwise independent of the vertical
! Eulerian coordinate.
!------------------------------------------------------------------------
                                                  call timing_on('Remapping')

         call Lagrangian_to_Eulerian(consv_te, ps, pe, delp,  &
                     pkz, pk, bdt, npz, is,ie,js,je, isd,ied,jsd,jed, &
                     nq, sphum, u,  v, w, delz, pt, q, phis, grav, zvir, cp_air,   &
                     akap, kord_mt, kord_tr, kord_tm, peln, te_2d,  &
                     ng, ua, va, omga, dp1, pem, fill, reproduce_sum, &
                     ak, bk, ks, ze0, remap_t, hydrostatic, hybrid_z, k_top)

                                                  call timing_off('Remapping')
!--------------------------
! Filter omega for physics:
!--------------------------
         if( nf_omega>0 )  then
                                           call timing_on('OMEGA_DEL2')
            call del2_cubed(omga, 0.25*da_min, npx, npy, npz, nf_omega)
                                           call timing_off('OMEGA_DEL2')
         endif
      endif
#endif

      deallocate ( dp1 )
      deallocate ( pem )

  if ( fv_debug ) then
       call prt_maxmin('PS_dyn', ps, is, ie, js, je, ng,   1, 0.01, gid==0)
       call prt_maxmin('T_dyn',  pt, is, ie, js, je, ng, npz, 1., gid==0)
  endif

  if( nwat==6 ) then
      call neg_adj3(is, ie, js, je, ng, npz,        &
                    pt, delp, q(isd,jsd,1,sphum),   &
                              q(isd,jsd,1,liq_wat), &
                              q(isd,jsd,1,rainwat), &
                              q(isd,jsd,1,ice_wat), &
                              q(isd,jsd,1,snowwat), &
                              q(isd,jsd,1,graupel)  )
     if ( fv_debug ) then
       call prt_maxmin('SPHUM_dyn',   q(isd,jsd,1,sphum  ), is, ie, js, je, ng, npz, 1., gid==0)
       call prt_maxmin('liq_wat_dyn', q(isd,jsd,1,liq_wat), is, ie, js, je, ng, npz, 1., gid==0)
       call prt_maxmin('ice_wat_dyn', q(isd,jsd,1,ice_wat), is, ie, js, je, ng, npz, 1., gid==0)
       call prt_maxmin('snowwat_dyn', q(isd,jsd,1,snowwat), is, ie, js, je, ng, npz, 1., gid==0)
       call prt_maxmin('graupel_dyn', q(isd,jsd,1,graupel), is, ie, js, je, ng, npz, 1., gid==0)
!      call prt_maxmin('cld_amt_dyn', q(isd,jsd,1,cld_amt), is, ie, js, je, ng, npz, 1., gid==0)
     endif
  endif

  end subroutine fv_dynamics


 subroutine del2_cubed(q, cd, npx, npy, km, ntimes)
!---------------------------------------------------------------
! This routine is for filtering the omega field for the physics
!---------------------------------------------------------------
   integer, intent(in):: npx, npy, km, ntimes
   real,    intent(in):: cd            ! cd = K * da_min;   0 < K < 0.25
   real, intent(inout):: q(isd:ied,jsd:jed,km)
   real, parameter:: r3  = 1./3.
   real :: fx(is:ie+1,js:je), fy(is:ie,js:je+1)
   integer i,j,k, n

   do n=1,ntimes
                     call timing_on('COMM_TOTAL')
   call mpp_update_domains(q, domain, complete=.true.)
                     call timing_off('COMM_TOTAL')
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

      do j=js,je
         do i=is,ie+1
            fx(i,j) = cd*dy(i,j)*sina_u(i,j)*(q(i-1,j,k)-q(i,j,k))*rdxc(i,j)
         enddo
      enddo

      do j=js,je+1
         do i=is,ie
            fy(i,j) = cd*dx(i,j)*sina_v(i,j)*(q(i,j-1,k)-q(i,j,k))*rdyc(i,j)
         enddo
      enddo

      do j=js,je
         do i=is,ie
            q(i,j,k) = q(i,j,k) + rarea(i,j)*(fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))
         enddo
      enddo
   enddo
   enddo

 end subroutine del2_cubed



 subroutine Rayleigh_Friction(dt, npz, ks, pm, tau, p_c, u, v, w, pt,  &
                              ua, va, cp, rg, hydrostatic, conserve)
    real, intent(in):: dt
    real, intent(in):: tau              ! time scale (days)
    real, intent(in):: p_c
    real, intent(in):: cp, rg
    real, intent(in),  dimension(npz):: pm
    integer, intent(in):: npz, ks
    logical, intent(in):: hydrostatic
    logical, intent(in):: conserve
    real, intent(inout):: u(isd:ied  ,jsd:jed+1,npz) ! D grid zonal wind (m/s)
    real, intent(inout):: v(isd:ied+1,jsd:jed,npz) ! D grid meridional wind (m/s)
    real, intent(inout)::  w(isd:ied,jsd:jed,npz) ! cell center vertical wind (m/s)
    real, intent(inout):: pt(isd:ied,jsd:jed,npz) ! temp
    real, intent(inout):: ua(isd:ied,jsd:jed,npz) ! 
    real, intent(inout):: va(isd:ied,jsd:jed,npz) ! 
    real, parameter:: sday = 86400.
    real, parameter:: wfac = 10.     ! factor to amplify the drag on w
    real c1, pc, fac
    integer i, j, k, kmax

     kmax = max(npz/3+1, ks)

     if ( .not. RF_initialized ) then
          allocate( rf(npz) )
          allocate( rw(npz) )

          if ( p_c <= 0. ) then
               pc = pm(1)
          else
               pc = p_c
          endif

          if( gid==0 ) write(6,*) 'Rayleigh friction E-folding time [days]:'
          c1 = 1. / (tau*sday)
          do k=1,kmax
             if ( pm(k) < 30.E2 ) then
                  rf(k) = c1*(1.+tanh(log10(pc/pm(k))))
                  if( gid==0 ) write(6,*) k, 0.01*pm(k), 1./(rf(k)*sday)
                  rf(k) = 1./(1.+dt*rf(k))
                  rw(k) = 1./(1.+dt*rf(k)*wfac)
             endif
          enddo
          RF_initialized = .true.
     endif

     if(conserve) call c2l_ord2(u, v, ua, va, dx, dy, rdxa, rdya, npz)

     do k=1,kmax
        if ( pm(k) < 30.E2 ) then
! Add heat so as to conserve TE
          if ( conserve ) then
               fac = 0.5*(1.-rf(k)**2) / (cp-rg*ptop/pm(k))
               do j=js,je
                  do i=is,ie
                     pt(i,j,k) = pt(i,j,k) + fac*(ua(i,j,k)**2 + va(i,j,k)**2)
                  enddo
               enddo
          endif
             do j=js,je+1
                do i=is,ie
                   u(i,j,k) = u(i,j,k)*rf(k)
                enddo
             enddo
             do j=js,je
                do i=is,ie+1
                   v(i,j,k) = v(i,j,k)*rf(k)
                enddo
             enddo
          if ( .not. hydrostatic ) then
             do j=js,je
                do i=is,ie
                   w(i,j,k) = w(i,j,k)*rw(k)
                enddo
             enddo
          endif
        endif
     enddo

 end subroutine Rayleigh_Friction

end module fv_dynamics_mod
