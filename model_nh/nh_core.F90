module nh_core_mod

! Notes:
! Using k_top=2 to treat the top layer hydrostatically so that delz will
! be computed using hydrostatic balance (instead of the update by
! advection of height using extrapolated winds at the model top)
!
! To do list:
! include moisture effect in pt
!------------------------------

   use constants_mod,  only: rdgas, grav
   use fv_control_mod, only: m_split, quick_p_c, quick_p_d, k_top, m_riem, master
   use tp_core_mod,    only: fv_tp_2d, copy_corners
!  use fv_timing_mod,  only: timing_on, timing_off

   implicit none
   private

   public Riem_Solver, Riem_Solver_C, update_dz_c, update_dz_d
   real, parameter:: dz_max = -0.5               ! (meters)

CONTAINS 

  subroutine update_dz_c(is, ie, js, je, km, ng, area,   &
                         zh, ut, vt, dz_in, dz_out, wk)
! !INPUT PARAMETERS:
  integer, intent(in):: is, ie, js, je, ng, km
  real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng,km):: ut, vt, zh
  real, intent(in ):: area(is-ng:ie+ng,js-ng:je+ng)
  real, intent(in ):: dz_in (is:ie,js:je,km) 
  real, intent(out):: dz_out(is:ie,js:je,km) 
  real, intent(out):: wk(is-ng:ie+ng,js-ng:je+ng,km+1)  ! work array
! Local Work array:
  real, dimension(is:ie+1,js:je  ):: xfx, fx
  real, dimension(is:ie  ,js:je+1):: yfx, fy
  integer  i, j, k

  do 6000 k=k_top,km

!-------------------------------------------------------------
     if ( k==1 ) then
#ifdef TOP_WIND
        do j=js,je
           do i=is,ie+1
              xfx(i,j) = ut(i,j,k) 
           enddo
        enddo
        do j=js,je+1
           do i=is,ie
              yfx(i,j) = vt(i,j,k)
           enddo
        enddo
#else
        call top_edge(ut(is:ie+1,js:je,  1:3), xfx, is, ie+1, js, je, 3)
        call top_edge(vt(is:ie  ,js:je+1,1:3), yfx, is, ie, js, je+1, 3)
#endif
     else
        do j=js,je
           do i=is,ie+1
              xfx(i,j) = 0.5*(ut(i,j,k-1) + ut(i,j,k))
           enddo
        enddo
        do j=js,je+1
           do i=is,ie
              yfx(i,j) = 0.5*(vt(i,j,k-1) + vt(i,j,k))
           enddo
        enddo
     endif
!-------------------------------------------------------------

     do j=js,je
        do i=is,ie+1
           if( xfx(i,j) > 0. ) then
               fx(i,j) = zh(i-1,j,k)
           else
               fx(i,j) = zh(i  ,j,k)
           endif
           fx(i,j) = xfx(i,j)*fx(i,j)
        enddo
     enddo
     do j=js,je+1
        do i=is,ie
           if( yfx(i,j) > 0. ) then
               fy(i,j) = zh(i,j-1,k)
           else
               fy(i,j) = zh(i,j  ,k)
           endif
           fy(i,j) = yfx(i,j)*fy(i,j)
        enddo
     enddo

! height increments due to horizontal advection:
     do j=js,je
        do i=is,ie
           wk(i,j,k) = (zh(i,j,k)*area(i,j) +  fx(i,j)- fx(i+1,j)+ fy(i,j)- fy(i,j+1)) &
                     / (          area(i,j) + xfx(i,j)-xfx(i+1,j)+yfx(i,j)-yfx(i,j+1)) &
                     - zh(i,j,k)
        enddo
     enddo
6000 continue

  do k=k_top,km-1
     do j=js,je
        do i=is,ie
           dz_out(i,j,k) = dz_in(i,j,k) - wk(i,j,k) + wk(i,j,k+1)
        enddo
     enddo
  enddo

! Bottom layer: k=km
  do j=js,je
     do i=is,ie
        dz_out(i,j,km) = dz_in(i,j,km) - wk(i,j,km) 
     enddo
  enddo
  call fix_dz(is, ie, js, je, km, 0, dz_out)

  end subroutine update_dz_c



  subroutine update_dz_d(hord, is, ie, js, je, km, ng, npx, npy, area,    &
                         zh, crx, cry, xfx, yfx, delz, wk, delp, n_sponge)

  integer, intent(in):: is, ie, js, je, ng, km, npx, npy
  integer, intent(in):: hord, n_sponge
  real, intent(in)   :: area(is-ng:ie+ng,js-ng:je+ng)
  real, intent(inout) ::  zh(is-ng:ie+ng,js-ng:je+ng,km)
  real, intent(inout) ::delz(is:ie,js:je,km)
  real, intent(inout) ::delp(is-ng:ie+ng,js-ng:je+ng,km)
  real, intent(inout), dimension(is:ie+1,js-ng:je+ng,km):: crx, xfx
  real, intent(inout), dimension(is-ng:ie+ng,js:je+1,km):: cry, yfx
  real, intent(  out) ::   wk(is:ie,js:je,km)  ! work array
!-----------------------------------------------------
! Local array:
  real, dimension(is:   ie+1, js-ng:je+ng):: crx_adv, xfx_adv
  real, dimension(is-ng:ie+ng,js:   je+1 ):: cry_adv, yfx_adv
  real, dimension(is:ie+1,js:je  ):: fx
  real, dimension(is:ie  ,js:je+1):: fy
  real  delx(is:ie+1,km), dely(is-ng:ie+ng,km)
  real :: ra_x(is:ie,js-ng:je+ng)
  real :: ra_y(is-ng:ie+ng,js:je)
!--------------------------------------------------------------------
  integer  i, j, k, iord, isd, ied, jsd, jed, lm

  isd = is - ng;  ied = ie + ng
  jsd = js - ng;  jed = je + ng

#ifdef NO_CS_PROFILE
  do k=k_top,km
     if( k <= n_sponge .and. km>16 ) then
         iord = 1
     else
         iord = hord
     endif

     if ( k==1 ) then
          call top_edge(crx,                      crx_adv, is,  ie+1, jsd, jed, km)
          call top_edge(xfx,                      xfx_adv, is,  ie+1, jsd, jed, km)
          call top_edge(cry(isd:ied,js:je+1,1:3), cry_adv, isd, ied,  js,  je+1, 3)
          call top_edge(yfx(isd:ied,js:je+1,1:3), yfx_adv, isd, ied,  js,  je+1, 3)
     else
        do j=jsd,jed
           do i=is,ie+1
              crx_adv(i,j) = 0.5*(crx(i,j,k-1) + crx(i,j,k))
              xfx_adv(i,j) = 0.5*(xfx(i,j,k-1) + xfx(i,j,k))
           enddo
        enddo
        do j=js,je+1
           do i=isd,ied
              cry_adv(i,j) = 0.5*(cry(i,j,k-1) + cry(i,j,k))
              yfx_adv(i,j) = 0.5*(yfx(i,j,k-1) + yfx(i,j,k))
           enddo
        enddo
     endif

     do j=jsd,jed
        do i=is,ie
           ra_x(i,j) = area(i,j) + xfx_adv(i,j) - xfx_adv(i+1,j)
        enddo
     enddo
     do j=js,je
        do i=isd,ied
           ra_y(i,j) = area(i,j) + yfx_adv(i,j) - yfx_adv(i,j+1)
        enddo
     enddo

     call fv_tp_2d(zh(isd,jsd,k), crx_adv, cry_adv, npx, npy, iord, &
                   fx, fy, xfx_adv, yfx_adv, area, ra_x, ra_y)
     do j=js,je
        do i=is,ie
           wk(i,j,k) = (zh(i,j,k)*area(i,j)+fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))       &
                     / (area(i,j)+xfx_adv(i,j)-xfx_adv(i+1,j)+yfx_adv(i,j)-yfx_adv(i,j+1)) &
                     - zh(i,j,k)
        enddo
     enddo
  enddo    ! k-loop

#else
! lm = 1   ! interior limiter
  lm = 0   ! no limiter

#ifndef UNIFORM_CS
  do k=1,km
     call copy_corners(delp(isd,jsd,k), npx, npy, 1)
  enddo

  do j=jsd,jed
     do k=1,km
        do i=is,ie+1
           delx(i,k) = 0.5*(delp(i-1,j,k)+delp(i,j,k))
        enddo
     enddo
     call edge_profile(crx, xfx, is, ie+1, jsd, jed, j, km, lm, delx)
  enddo

  do k=1,km
     call copy_corners(delp(isd,jsd,k), npx, npy, 2)
  enddo

  do j=js,je+1
     do k=1,km
        do i=isd,ied
           dely(i,k) = 0.5*(delp(i,j-1,k)+delp(i,j,k))
        enddo
     enddo
     call edge_profile(cry, yfx, isd, ied, js, je+1, j, km, lm, dely)
  enddo

#else
  do j=jsd,jed
     call edge_profile(crx, xfx, is,  ie+1, jsd, jed, j, km, lm)
  enddo
  do j=js,je+1
     call edge_profile(cry, yfx, isd, ied,  js, je+1, j, km, lm)
  enddo
#endif

  do k=k_top,km

     do j=jsd,jed
        do i=is,ie
           ra_x(i,j) = area(i,j) + xfx(i,j,k) - xfx(i+1,j,k)
        enddo
     enddo
     do j=js,je
        do i=isd,ied
           ra_y(i,j) = area(i,j) + yfx(i,j,k) - yfx(i,j+1,k)
        enddo
     enddo

     call fv_tp_2d(zh(isd,jsd,k), crx(is,jsd,k), cry(isd,js,k), npx,  npy, iord, &
                   fx, fy, xfx(is,jsd,k), yfx(isd,js,k), area, ra_x, ra_y)
     do j=js,je
        do i=is,ie
           wk(i,j,k) = (zh(i,j,k)*area(i,j)+fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))   &
                     / (area(i,j)+xfx(i,j,k)-xfx(i+1,j,k)+yfx(i,j,k)-yfx(i,j+1,k)) &
                     - zh(i,j,k)
        enddo
     enddo
  enddo
#endif

  do k=k_top,km-1
     do j=js,je
        do i=is,ie
           delz(i,j,k) = delz(i,j,k) - wk(i,j,k) + wk(i,j,k+1)
        enddo
     enddo
  enddo
  do j=js,je
     do i=is,ie
        delz(i,j,km) = delz(i,j,km) - wk(i,j,km) 
     enddo
  enddo

  call fix_dz(is, ie, js, je, km, 0, delz)

  end subroutine update_dz_d


  subroutine Riem_Solver_C(dt,   is,  ie,   js, je, km,   ng,  &
                           akap, cp,  ptop, hs, w,  delz, pt,  &
                           delp, gz,  pk,   ip)

   integer, intent(in):: is, ie, js, je, ng, km
   integer, intent(in):: ip       ! ip==1 pk is full pressure
   real, intent(in):: dt,  akap, cp, ptop
   real, intent(in):: delz(is:ie,js:je,km)
   real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng,km):: pt, delp
   real, intent(in)::       hs(is-ng:ie+ng,js-ng:je+ng)
   real, intent(inout):: w(is-ng:ie+ng,js-ng:je+ng,km)
! OUTPUT PARAMETERS 
   real, intent(out), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: gz, pk
! Local:
  real, dimension(is:ie,km  ):: pm, dm, dz2
  real, dimension(is:ie,km+1):: pem, pk2
  real gama, rgrav, ptk
  integer i, j, k
  integer m_split_c


   m_split_c = max(1, m_split/2)
    gama = 1./(1.-akap)
   rgrav = 1./grav

   ptk = ptop ** akap

!$omp parallel do default(shared) private(i, j, k, dm, dz2, pem, pm, pk2)
   do 2000 j=js,je

      do k=1,km
         do i=is,ie
            dm(i,k) = delp(i,j,k)
         enddo
      enddo

      do i=is,ie
         pem(i,1) = ptop
         pk2(i,1) = ptk
      enddo

      do k=2,km+1
         do i=is,ie
            pem(i,k) = pem(i,k-1) + dm(i,k-1)
            pk2(i,k) = exp( akap*log(pem(i,k)) )
         enddo
      enddo

      do k=k_top,km
         do i=is,ie
            dz2(i,k) = delz(i,j,k)
             pm(i,k) = exp( gama*log(akap*dm(i,k)/(pk2(i,k+1)-pk2(i,k))) )
             dm(i,k) = dm(i,k) * rgrav
         enddo
      enddo

      call Riem_3D( m_split_c, dt, is, ie, js, je, ng, j, km, cp, gama, akap, &
                    pk, dm, pm, w, dz2, pt, quick_p_c, .true., k_top, m_riem )

      do i=is,ie
         pk(i,j,1) = ptop                     ! full pressure at top
      enddo
      do k=2,km+1
         do i=is,ie
            pk(i,j,k) = pk(i,j,k) + pem(i,k)  ! add hydrostatic component
         enddo
      enddo

!---------------------------------------
! Compute dz2 hydrostatically if k_top>1
!---------------------------------------
   if ( k_top>1 ) then
      do k=1,k_top-1
         do i=is,ie
            dz2(i,k) = pt(i,j,k)*(pk2(i,k)-pk2(i,k+1))*rgrav
         enddo
      enddo
   endif

! Compute Height * grav (for p-gradient computation)
      do i=is,ie
         gz(i,j,km+1) = hs(i,j)
      enddo

      do k=km,1,-1
         do i=is,ie
            gz(i,j,k) = gz(i,j,k+1) - dz2(i,k)*grav
         enddo
      enddo

2000  continue

  end subroutine Riem_Solver_C


  subroutine Riem_Solver(dt,   is,   ie,   js, je, km, ng,    &
                         akap, cp,   ptop, hs, peln, w,  delz, pt,  &
                         delp, gz,   pkc, pk, pe, last_call, ip)
!--------------------------------------------
! !OUTPUT PARAMETERS
! Ouput: gz: grav*height at edges
!        pe: full     hydrostatic pressure
!       pkc: full non-hydrostatic pressure
!--------------------------------------------
   integer, intent(in):: is, ie, js, je, km, ng
   integer, intent(in):: ip      ! ip==0 pkc is perturbation pressure
   real, intent(in):: dt         ! the BIG horizontal Lagrangian time step
   real, intent(in):: akap, cp, ptop
   real, intent(in):: hs(is-ng:ie+ng,js-ng:je+ng)
   logical, intent(in):: last_call
   real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,km):: w, delp, pt
   real, intent(inout):: delz(is:ie,js:je,km)
   real, intent(out), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: gz, pkc
   real, intent(out):: pk(is:ie,js:je,km+1)
   real, intent(out):: pe(is-1:ie+1,km+1,js-1:je+1)
   real, intent(out):: peln(is:ie,km+1,js:je)           ! ln(pe)
! Local:
  real, dimension(is:ie,km):: pm, dm, dz2
  real :: pem(is:ie,km+1)
  real gama, rgrav, ptk
  integer i, j, k

    gama = 1./(1.-akap)
   rgrav = 1./grav
     ptk = ptop ** akap

!$omp parallel do default(shared) private(i, j, k, dm, dz2, pem, pm)
   do 2000 j=js,je

      do k=1,km
         do i=is,ie
            dm(i,k) = delp(i,j,k)
         enddo
      enddo

      do i=is,ie
         pem(i,1) = ptop
         pk(i,j,1) = ptk
      enddo

      do k=2,km+1
         do i=is,ie
               pem(i,k) = pem(i,k-1) + dm(i,k-1)
            peln(i,k,j) = log(pem(i,k))
              pk(i,j,k) = exp(akap*peln(i,k,j))
         enddo
      enddo

      do k=k_top,km
         do i=is,ie
            dz2(i,k) = delz(i,j,k)
! hydrostatic pressure:
!           pm(i,k) = (akap*dm(i,k)/(pk(i,j,k+1)-pk(i,j,k))) ** gama
            pm(i,k) = exp( gama*log(akap*dm(i,k)/(pk(i,j,k+1)-pk(i,j,k))) )
            dm(i,k) = dm(i,k) * rgrav
         enddo
      enddo

      call Riem_3D(m_split, dt, is, ie, js, je, ng, j, km, cp, gama, akap,  &
                   pkc, dm, pm, w, dz2, pt, quick_p_d, .false., k_top, m_riem)
      if ( ip==1 ) then
           do i=is,ie
              pkc(i,j,1) = ptop
           enddo
           do k=2,km+1
              do i=is,ie
                 pkc(i,j,k) = pkc(i,j,k) + pem(i,k)
              enddo
           enddo
      endif

!---------------------------------------
! Compute dz2 hydrostatically if k_top>1
!---------------------------------------
   if ( k_top>1 ) then
      do k=1,k_top-1
         do i=is,ie
            dz2(i,k) = pt(i,j,k)*(pk(i,j,k)-pk(i,j,k+1))*rgrav
         enddo
      enddo
   endif

! Compute Height * grav (for p-gradient computation)
      do i=is,ie
         gz(i,j,km+1) = hs(i,j)
      enddo

      do k=km,1,-1
         do i=is,ie
            gz(i,j,k) = gz(i,j,k+1) - dz2(i,k)*grav
         enddo
      enddo

      do k=1,km
         do i=is,ie
            delz(i,j,k) = dz2(i,k)
         enddo
      enddo

2000  continue

  if ( last_call ) then

    do j=js-1,je+1
       do i=is-1,ie+1
          pe(i,1,j) = ptop
       enddo
       do k=2,km+1
          do i=is-1,ie+1
             pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
          enddo
       enddo
    enddo

  endif

  end subroutine Riem_Solver


  subroutine Riem_3D(ns, bdt, is, ie, js, je, ng, j, km, cp, gama, cappa, p3, dm2,    &
                     pm2, w, dz2, pt, quick_p, c_core, ktop, iad)

  integer, intent(in):: ns, is, ie, js, je, ng,  km, j
  integer, intent(in):: iad      ! time step scheme 
  integer, intent(in):: ktop     ! starting layer for non-hydrostatic dynamics
                                 ! 1: All non-hydrostatic
                                 ! 2: top sponge layer is hydrostatic
  real,    intent(in):: bdt, cp, gama, cappa
  real,    intent(in), dimension(is:ie,km):: dm2, pm2
  logical, intent(in):: quick_p       ! fast algorithm for pressure
  logical, intent(in):: c_core
  real, intent(in  ) :: pt (is-ng:ie+ng,js-ng:je+ng,km)
! IN/OUT:
  real, intent(inout):: dz2(is:ie,km)
  real, intent(inout)::   w(is-ng:ie+ng,js-ng:je+ng,km)
  real, intent(out  )::  p3(is-ng:ie+ng,js-ng:je+ng,km+1)
! --- Local 1D copyies -----
#ifdef USE_2D
  real, dimension(km,is:ie):: t2, p2, pt2
#else
  real, dimension(km):: c2, p2, pt2
#endif
  real, dimension(km):: r_p, r_n, rden, dz, dm, wm, dts, pdt
  real, dimension(km+1):: m_bot, m_top, r_bot, r_top, time_left, pe1, pbar, wbar

  real, parameter:: dzmx = 0.5*dz_max
  real    :: dt, rdt, grg, z_frac, t_left
  real    :: a1, b1, g2, rcp
  real    :: seq(ns)       ! time stepping sequence
  integer :: k2(km+1)
  integer :: i, k, n, ke, kt, k0, k1, k3

  call time_sequence( iad, ns, bdt, seq )

  grg = gama * rdgas  
  rcp = 1. / cp
  rdt = 1. / bdt

  if ( quick_p ) then
       a1 = 0.5               ! a1=1 fully implicit
       b1 = 1. - a1
       g2 = -2.*gama
  endif

#ifdef USE_2D
  do i=is,ie
     do k=ktop,km
        rden(k) = -rdgas*dm2(i,k)/dz2(i,k)
        pt2(k,i) = pt(i,j,k) * rcp
         p2(k,i) = ( rden(k)*pt2(k,i) )**gama
         t2(k,i) = p2(k,i) / rden(k)
     enddo
  enddo
#endif

  do k=1,km+1
     wbar(k) = 0.
     pbar(k) = 0.
  enddo

 do 6000 i=is,ie

    do k=ktop,km
#ifdef USE_2D
       dz(k) = dz2(i,k)
       dm(k) = dm2(i,k)
       wm(k) = w(i,j,k)*dm(k)
#else
       dz(k) = dz2(i,k)
       dm(k) = dm2(i,k)
       wm(k) = w(i,j,k)*dm(k)
       rden(k) = -rdgas*dm(k)/dz(k)
       pt2(k) = pt(i,j,k) * rcp   ! virtual effect is included in pt
!      p2(k) = ( rden(k)*pt2(k) )**gama
       p2(k) = exp( gama*log(rden(k)*pt2(k)) )
       c2(k) = sqrt( grg*p2(k)/rden(k) )
#endif
    enddo

    do k=1,km+1
       pe1(k) = 0.
    enddo

 do 5000 n=1,ns

   dt = seq(n)

   do k=ktop,km
#ifdef USE_2D
       dts(k) = -dz(k)/(sqrt(grg*t2(k,i)))
       pdt(k) = dts(k)*(p2(k,i)-pm2(i,k))
#else
       dts(k) = -dz(k) / c2(k)
       pdt(k) = dts(k)*(p2(k)-pm2(i,k))
#endif
       r_p(k) = wm(k) + pdt(k)
       r_n(k) = wm(k) - pdt(k)
   enddo

!--------------------------------------------------
! Compute r_top from bottom up: dm/dt > 0
!----------------------------------------------------
   do k=ktop+1,km+1
      k2(k) = k-1
      m_top(k) = 0.
      r_top(k) = 0.
      time_left(k) = dt
   enddo
  
   do 444 ke=km+1,ktop+1,-1
        kt=k2(ke)
     do k=kt,ktop,-1
        z_frac = time_left(ke)/dts(k)
        if ( z_frac <= 1. ) then
            if ( (ke-k) > 2 ) then
               k1 = ke-1
               k2(k1) = k
               m_top(k1) = m_top(ke) - dm(k1)
               r_top(k1) = r_top(ke) - r_n(k1)
               time_left(k1) = time_left(ke) + dts(k1)
            endif
            m_top(ke) = m_top(ke) + z_frac*dm(k)
            r_top(ke) = r_top(ke) + z_frac*r_n(k)
            go to 444
        else
            time_left(ke) = time_left(ke) - dts(k)
            m_top(ke) = m_top(ke) + dm(k)
            r_top(ke) = r_top(ke) + r_n(k)
        endif 
     enddo
! wave from ke already left the top
     if ( ke == ktop+1 ) exit
     do k=ke-1,ktop+1,-1
        m_top(k) = m_top(k+1) - dm(k)
        r_top(k) = r_top(k+1) - r_n(k)
     enddo
     exit
444 continue

!--------------------------------------------------
! Compute r_bot from top down: dm/dt < 0
!----------------------------------------------------
   do k=ktop,km
        k2(k) = k
     m_bot(k) = 0.
     r_bot(k) = 0.
    time_left(k) = dt
   enddo

  do 4000 ke=ktop,km
        kt = k2(ke)
     do k=kt,km
        z_frac = time_left(ke)/dts(k)
        if ( z_frac <= 1. ) then
             if ( (k-ke)>1 ) then
                time_left(ke+1) = time_left(ke) + dts(ke)
                m_bot(ke+1) =  m_bot(ke) - dm(ke)
                r_bot(ke+1) =  r_bot(ke) - r_p(ke)
                k2(ke+1) = k
             endif
                m_bot(ke) = m_bot(ke) + z_frac*dm(k)
                r_bot(ke) = r_bot(ke) + z_frac*r_p(k)
             if( ke==km ) go to 7777      ! All done
             go to 4000      ! to next interface
        else 
             time_left(ke) = time_left(ke) - dts(k)
             m_bot(ke) =  m_bot(ke) + dm(k)
             r_bot(ke) =  r_bot(ke) + r_p(k)
        endif
     enddo
!----------------------------------------
! Ray from edge-ke already hit the ground.
!----------------------------------------
     k3 = ke
     t_left = time_left(ke)
     exit
4000 continue


!---------------------------------
! Perfect reflection at the bottom
!---------------------------------
   k1 = km
   do kt=k3,km
     k0 = k1
     do k=k0,ktop,-1
        z_frac = t_left/dts(k)
        if ( z_frac <= 1. ) then
!-- next interface -------------------------------------
!          if ( kt /= km ) then
                    k1 = k
                 t_left = t_left + dts(kt)
                 m_bot(kt+1) = m_bot(kt) - dm(kt)
                 r_bot(kt+1) = r_bot(kt) - r_p(kt)
!          endif
!-------------------------------------------------------
           m_bot(kt) = m_bot(kt) + z_frac*dm(k)
           r_bot(kt) = r_bot(kt) - z_frac*r_n(k)
           exit            ! goto next interface
        else 
           m_bot(kt) = m_bot(kt) + dm(k)
           r_bot(kt) = r_bot(kt) - r_n(k)
              t_left = t_left - dts(k)
        endif
     enddo
   enddo

7777  continue

  pbar(ktop) = 0.
  wbar(ktop) = r_bot(ktop) / m_bot(ktop)
  do k=ktop+1,km
     wbar(k) = (r_bot(k)+r_top(k)) / (m_top(k)+m_bot(k))
     pbar(k) =  m_top(k)*wbar(k) - r_top(k)
  enddo
  pbar(km+1) = -r_top(km+1)
! wbar(km+1) = 0.

   do k=ktop+1,km+1
      pe1(k) = pe1(k) + pbar(k)
   enddo

   if ( n==ns ) then
      if ( c_core ) then
          do k=ktop,km
             dz2(i,k) = min(dzmx, dz(k) + dt*(wbar(k+1)-wbar(k)) )
          enddo
      else
          do k=ktop,km
             dz2(i,k) = min(dzmx, dz(k) + dt*(wbar(k+1)-wbar(k)) )
             w(i,j,k) = ( wm(k) + pbar(k+1) - pbar(k) ) / dm(k)
          enddo
      endif
   else
     if ( quick_p ) then
        do k=ktop,km
           wm(k) = wm(k) + pbar(k+1) - pbar(k)
           rden(k) = dt*(wbar(k+1)-wbar(k))   ! dz tendency
           pdt(k) = dz(k)                     ! old dz
           dz(k) = dz(k) + rden(k)            ! updated dz
           pdt(k) = g2*rden(k) / (pdt(k)+dz(k))
#ifdef USE_2D
           p2(k,i) = max(0., p2(k,i)*(1.+b1*pdt(k))/(1.-a1*pdt(k)))
           t2(k,i) = -p2(k,i)*dz(k)/(rdgas*dm(k))
#else
           p2(k) = max(0., p2(k)*(1.+b1*pdt(k))/(1.-a1*pdt(k)))
           c2(k) = sqrt( -gama*p2(k)*dz(k)/dm(k) )
#endif
        enddo
      else
        do k=ktop,km
           wm(k) = wm(k) + pbar(k+1) - pbar(k)
           dz(k) = min(dzmx, dz(k) + dt*(wbar(k+1)-wbar(k)) )
           rden(k) = -rdgas*dm(k)/dz(k)
#ifdef USE_2D
           p2(k,i) = (rden(k)*pt2(k,i))**gama
           t2(k,i) = p2(k,i)/rden(k)
#else
!          p2(k) = (rden(k)*pt2(k))**gama
           p2(k) = exp( gama*log(rden(k)*pt2(k)) )
           c2(k) = sqrt( grg*p2(k)/rden(k) )
#endif
        enddo
      endif
   endif

5000  continue

!---------------------------------------------------------------------
! Note: time-mean full pressure at edges and time-level-(n+1) DZ are
! used for computation of p-gradient
! Could use cell center full pressure at time-level-(n+1)
!---------------------------------------------------------------------
!        p3(i,j,1:ktop) = 0.
      do k=1,ktop
         p3(i,j,k) = 0.
      enddo
      do k=ktop+1,km+1
         p3(i,j,k) = pe1(k)*rdt
      enddo

6000  continue

 end subroutine Riem_3D

 subroutine time_sequence ( iad, ns, bdt, tseq )
 integer, intent(in) :: iad, ns
 real,    intent(in) :: bdt
 real, intent(out):: tseq(ns)
! local
  integer :: seq(ns)
  integer :: n, nstep
  real :: sdt

! Default uniform time stepping (iad=0)
  do n=1,ns
     seq(n) = 1
  enddo

! Note: iad=2 or 4 appear to be more stable than other options
  if ( ns>3 ) then
    if ( iad==1 ) then
                     ! 1, 1, 2, 2, ...
         do n=3,ns
            seq(n) = 2
         enddo
    elseif ( iad==2 ) then
                     ! 1, 2, 2, 2, ...
         do n=2,ns
            seq(n) = 2
         enddo
    elseif ( iad==3 ) then
                     ! 1, 2, 3, 3, 3. ...
         seq(2) = 2
         do n=3,ns
            seq(n) = 3
         enddo
    elseif ( iad==4 ) then
                     ! 1, 2, 4, 4, 4, ...
         seq(2) = 2
         do n=3,ns
            seq(n) = 4
         enddo
    elseif( iad==5 ) then
!---------------------
! Fibonacci sequence:
!---------------------
                     ! 1, 1, 2, 3, 5, 8, 13, 21, 34
         do n=3,ns
            seq(n) = seq(n-2) + seq(n-1)
         enddo
    endif
  endif

  nstep = 1
  if ( ns>1 ) then
       do n=2,ns
          nstep = nstep + seq(n) 
       enddo
  endif
  sdt = bdt / real(nstep)

  do n=1,ns
     tseq(n) = sdt * real ( seq(n) )
  enddo


 end subroutine time_sequence


 subroutine edge_profile(q1, q2, i1, i2, j1, j2, j, km, limiter, delp)
! Optimized for wind profile reconstruction:
! Developer: S.-J. Lin, NOAA/GFDL
 integer, intent(in):: i1, i2, j1, j2
 integer, intent(in):: j, km
 integer, intent(in):: limiter
 real, intent(in), optional::  delp(i1:i2,km)  ! layer thickness
 real, intent(inout), dimension(i1:i2,j1:j2,km):: q1, q2
!-----------------------------------------------------------------------
 real, dimension(i1:i2,km+1):: qe1, qe2  ! edge values
 real   d4(i1:i2)
 real  gam(i1:i2,km)
 real  gak(km)
 real  bet, a_bot, gratio, r2o3, r4o3, xt1, xt2
 integer i, k

 if ( present(delp) ) then
  do i=i1,i2
        gratio = delp(i,2) / delp(i,1)   ! grid ratio
           xt1 = 2.*gratio*(gratio+1. )
           bet =    gratio*(gratio+0.5)
      qe1(i,1) = ( xt1*q1(i,j,1) + q1(i,j,2) ) / bet
      qe2(i,1) = ( xt1*q2(i,j,1) + q2(i,j,2) ) / bet
      gam(i,1) = ( 1. + gratio*(gratio+1.5) ) / bet
  enddo

  do k=2,km
     do i=i1,i2
           d4(i) = delp(i,k-1) / delp(i,k)
             bet =  2. + 2.*d4(i) - gam(i,k-1)
        qe1(i,k) = ( 3.*(q1(i,j,k-1)+d4(i)*q1(i,j,k)) - qe1(i,k-1) ) / bet
        qe2(i,k) = ( 3.*(q2(i,j,k-1)+d4(i)*q2(i,j,k)) - qe2(i,k-1) ) / bet
        gam(i,k) = d4(i) / bet
     enddo
  enddo
 
  do i=i1,i2
           a_bot = 1. + d4(i)*(d4(i)+1.5)
             xt1 = 2.*d4(i)*(d4(i) + 1.)
             xt2 = d4(i)*(d4(i)+0.5) - a_bot*gam(i,km)
     qe1(i,km+1) = ( xt1*q1(i,j,km) + q1(i,j,km-1) - a_bot*qe1(i,km) ) / xt2
     qe2(i,km+1) = ( xt1*q2(i,j,km) + q2(i,j,km-1) - a_bot*qe2(i,km) ) / xt2
  enddo

  do k=km,1,-1
     do i=i1,i2
        qe1(i,k) = qe1(i,k) - gam(i,k)*qe1(i,k+1)
        qe2(i,k) = qe2(i,k) - gam(i,k)*qe2(i,k+1)
     enddo
  enddo

 else
!------------------------------------------------
! Optimized coding for uniform grid: SJL Apr 2007
!------------------------------------------------
     r2o3 = 2./3.
     r4o3 = 4./3.
     do i=i1,i2
        qe1(i,1) = r4o3*q1(i,j,1) + r2o3*q1(i,j,2)
        qe2(i,1) = r4o3*q2(i,j,1) + r2o3*q2(i,j,2)
     enddo

        gak(1) = 7./3.
     do k=2,km
        gak(k) =  1. / (4. - gak(k-1))
        do i=i1,i2
           qe1(i,k) = (3.*(q1(i,j,k-1) + q1(i,j,k)) - qe1(i,k-1)) * gak(k)
           qe2(i,k) = (3.*(q2(i,j,k-1) + q2(i,j,k)) - qe2(i,k-1)) * gak(k)
        enddo
     enddo

     bet = 1. / (1.5 - 3.5*gak(km))
     do i=i1,i2
        qe1(i,km+1) = (4.*q1(i,j,km) + q1(i,j,km-1) - 3.5*qe1(i,km)) * bet
        qe2(i,km+1) = (4.*q2(i,j,km) + q2(i,j,km-1) - 3.5*qe2(i,km)) * bet
     enddo

     do k=km,1,-1
        do i=i1,i2
           qe1(i,k) = qe1(i,k) - gak(k)*qe1(i,k+1)
           qe2(i,k) = qe2(i,k) - gak(k)*qe2(i,k+1)
        enddo
     enddo
 endif

!------------------
! Apply constraints
!------------------
    if ( limiter==2 ) then   ! limit the top winds
         do i=i1,i2
            if ( q1(i,j,1)*qe1(i,1) <= 0. ) qe1(i,1) = 0.
            if ( q2(i,j,1)*qe2(i,1) <= 0. ) qe2(i,1) = 0.
         enddo
    endif

    if ( limiter/=0 ) then
         do k=2,km
            do i=i1,i2
               qe1(i,k) = min( qe1(i,k), max(q1(i,j,k-1), q1(i,j,k)) )
               qe1(i,k) = max( qe1(i,k), min(q1(i,j,k-1), q1(i,j,k)) )
               qe2(i,k) = min( qe2(i,k), max(q2(i,j,k-1), q2(i,j,k)) )
               qe2(i,k) = max( qe2(i,k), min(q2(i,j,k-1), q2(i,j,k)) )
            enddo
         enddo
    endif

! mean values replaced by edge winds
    do k=1,km
       do i=i1,i2
          q1(i,j,k) = qe1(i,k)
          q2(i,j,k) = qe2(i,k)
       enddo
    enddo

 end subroutine edge_profile


 subroutine top_edge(p, qe, is, ie, js, je, km)
! Constant grid spacing:
  integer, intent(in) :: is, ie, js, je, km
  real,    intent(in) ::  p(is:ie,js:je,km)
  real, intent(out)::    qe(is:ie,js:je)
!----------------------------------------
  real dq1, dq2
  real a3, b2, sc
  integer i,j

! three-cell parabolic subgrid distribution at model top
! three-cell PP-distribution
! Compute a,b, and c of q = aP**2 + bP + c using cell averages and delp
! a3 = a / 3
! b2 = b / 2

  do j=js,je
     do i=is,ie
         dq1 = p(i,j,2) - p(i,j,1)
         dq2 = p(i,j,3) - p(i,j,2)
         a3 = (dq2 - 1.5*dq1) / 6.
         if( abs(a3) > 1.e-15 ) then
             b2 =  0.5*dq1 - 3.*a3
             sc = -b2/(3.*a3)
            if(sc < 0. .or. sc > 1.) then
               qe(i,j) = p(i,j,1) - (a3 + b2)
            else
               qe(i,j) = p(i,j,1) - 0.5*dq1
            endif
         else
! Linear profile:
            qe(i,j) = p(i,j,1) - 0.5*dq1
         endif
         if ( qe(i,j)*p(i,j,1) < 0. ) qe(i,j) = 0.
     enddo
  enddo

 end subroutine top_edge

 subroutine fix_dz(is, ie, js, je, km, ng, dz)
   integer,  intent(in):: is, ie, km
   integer,  intent(in):: js, je, ng
   real,  intent(inout):: dz(is:ie, js-ng:je+ng,km)
   integer i, j, k
   logical modified

   modified = .false.

!$omp parallel do default(shared) private(i, j, k)
   do j=js,je
      do k=km,k_top+1,-1
         do i=is,ie
            if( dz(i,j,k) > dz_max ) then
                dz(i,j,k-1) = dz(i,j,k-1) + dz(i,j,k) - dz_max
                dz(i,j,k  ) = dz_max
!               modified = .true.
            endif
         enddo
      enddo
! Top layer
      do i=is,ie
         if( dz(i,j,k_top) > dz_max ) then
             dz(i,j,k_top+1) = dz(i,j,k_top+1) + dz(i,j,k_top) - dz_max
             dz(i,j,k_top) = dz_max
!            modified = .true.
         endif
      enddo
   enddo

#ifdef NH_DEBUG
   if ( modified) write(*,*) 'DZ modified'
#endif

end subroutine fix_dz
!-----------------------------------------------------------------------

#ifdef DEEP_ATM
  subroutine rotate_uvw(dt, im, jm, km, jfirst, jlast, ng_d, ng_s, g_d,  &
                        ua, va, u, v, w, du)
  integer, intent(in):: im, jm, km
  integer, intent(in):: ng_d, ng_s
  integer, intent(in):: jfirst       ! first latitude of the subdomain
  integer, intent(in):: jlast        ! last latitude of the subdomain
  real, intent(in)::  dt
  real, intent(in)::  g_d(im,jfirst:jlast)

  real, intent(inout):: ua(im,jfirst:jlast,km)  ! a-grid u-Wind (m/s)
  real, intent(inout):: va(im,jfirst:jlast,km)  ! a-grid v-Wind (m/s)
  real, intent(inout)::  u(im,jfirst-ng_d:jlast+ng_s,km)  ! u-Wind (m/s)
  real, intent(inout)::  v(im,jfirst-ng_d:jlast+ng_d,km)  ! v-Wind (m/s)
  real, intent(inout)::  w(im,jfirst-ng_d:jlast+ng_d,km)  ! w-Wind (m/s)
  real, intent(out):: du(im,jfirst-1:jlast,km)

! Local:
  real wo(0:im)
  integer i,j,k
  
!-------------------------
! Fully explicit algorithm
!-------------------------
! Need to add mean height to radius!!
  do k=1,km
     do j=js,je+1
        do i=is,ie
!  A grid tendency:
           du(i,j,k) = -dt*w(i,j,k)*( g_d(i,j) + ua(i,j,k)/(radius+mz(i,j,k)) )
        enddo
     enddo
     do j=js,je
        do i=is,ie+1
           v(i,j,k) = v(i,j,k) * ( 1. - 0.5*dt*(wo(i-1)+wo(i))/radius ) 
        enddo
     enddo
     do j=js,je
        do i=is,ie
           w(i,j,k) = wo(i) + dt * ( g_d(i,j)*ua(i,j,k) +    &
                             (ua(i,j,k)**2+va(i,j,k)**2)/radius )
        enddo
     enddo
  enddo

  end subroutine rotate_uvw
#endif

end module nh_core_mod
