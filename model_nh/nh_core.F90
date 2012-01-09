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
   use fv_control_mod, only: m_split, k_top, m_riem, master, kd3, fill_wz
   use tp_core_mod,    only: fv_tp_2d, copy_corners

   implicit none
   private

   public Riem_Solver, update_gz_c, update_dz_d
   real, parameter:: dz_max = -2.5               ! (meters)
   real, parameter:: df_max = dz_max * grav
   real, parameter:: near0 = 1.E-8

   real, allocatable:: density(:,:,:)
   real:: damp3

!---- version number -----
   character(len=128) :: version = '$Id: nh_core.F90,v 19.0 2012/01/06 19:57:52 fms Exp $'
   character(len=128) :: tagname = '$Name: siena $'

contains

  subroutine update_gz_c(is, ie, js, je, km, ng, area, ut, vt, gz)
! !INPUT PARAMETERS:
  integer, intent(in):: is, ie, js, je, ng, km
  real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng,km):: ut, vt
  real, intent(in):: area(is-ng:ie+ng,js-ng:je+ng)
  real, intent(inout):: gz(is-ng:ie+ng,js-ng:je+ng,km+1)  ! work array
! Local Work array:
  real, dimension(is:ie+1,js:je  ):: xfx, fx
  real, dimension(is:ie  ,js:je+1):: yfx, fy
  integer  i, j, k

!$omp parallel do default(shared) private(xfx, yfx, fx, fy)
  do 6000 k=k_top,km
! Use winds below the height
     if ( k==1 ) then
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
     do j=js,je
        do i=is,ie+1
           if( xfx(i,j) > 0. ) then
               fx(i,j) = gz(i-1,j,k)
           else
               fx(i,j) = gz(i  ,j,k)
           endif
           fx(i,j) = xfx(i,j)*fx(i,j)
        enddo
     enddo
     do j=js,je+1
        do i=is,ie
           if( yfx(i,j) > 0. ) then
               fy(i,j) = gz(i,j-1,k)
           else
               fy(i,j) = gz(i,j  ,k)
           endif
           fy(i,j) = yfx(i,j)*fy(i,j)
        enddo
     enddo

     do j=js,je
        do i=is,ie
           gz(i,j,k) = (gz(i,j,k)*area(i,j) +  fx(i,j)- fx(i+1,j)+ fy(i,j)- fy(i,j+1)) &
                     / (          area(i,j) + xfx(i,j)-xfx(i+1,j)+yfx(i,j)-yfx(i,j+1))
        enddo
     enddo

#ifdef FIX_BOTTOM_LAYER
     if ( k==km ) then
        do j=js,je
           do i=is,ie
! Limit gz at top edge of the bottom layer; this essentially limited the mountain forcing on ws
! Prevents negative layer thickness
                  gz(i,j,km) = max(gz(i,j,km+1)-df_max, gz(i,j,km)) 
           enddo
        enddo
     endif
#endif

6000 continue

  end subroutine update_gz_c



  subroutine update_dz_d(ndif, damp, hord, is, ie, js, je, km, ng, npx, npy, area,    &
                         zs, zh, crx, cry, xfx, yfx, delz, ws, rdt, id_ws)

  integer, intent(in):: is, ie, js, je, ng, km, npx, npy
  integer, intent(in):: ndif(km)
  integer, intent(in):: hord
  integer, intent(in):: id_ws
  real, intent(in)   :: rdt
  real, intent(in)   :: damp(km)
  real, intent(in)   :: area(is-ng:ie+ng,js-ng:je+ng)
  real, intent(inout) ::  zs(is-ng:ie+ng,js-ng:je+ng)
  real, intent(inout) ::  zh(is-ng:ie+ng,js-ng:je+ng,km)
  real, intent(inout) ::delz(is:ie,js:je,km)
  real, intent(inout), dimension(is:ie+1,js-ng:je+ng,km):: crx, xfx
  real, intent(inout), dimension(is-ng:ie+ng,js:je+1,km):: cry, yfx
  real, intent(inout)   :: ws(is:ie,js:je)
!-----------------------------------------------------
  integer, parameter:: kint = 0 ! top layer(s) using delz transport
! Local array:
  real, dimension(is:   ie+1, js-ng:je+ng):: crx_adv, xfx_adv
  real, dimension(is-ng:ie+ng,js:   je+1 ):: cry_adv, yfx_adv
  real, dimension(is:ie+1,js:je  ):: fx
  real, dimension(is:ie  ,js:je+1):: fy
  real:: dz1(is-ng:ie+ng,js-ng:je+ng)
  real:: ra_x(is:ie,js-ng:je+ng)
  real:: ra_y(is-ng:ie+ng,js:je)
  real, parameter:: surf_fac = 0.9  ! lowest layer to "surface" reduction factor
!--------------------------------------------------------------------
  integer  i, j, k, isd, ied, jsd, jed, kmax

  isd = is - ng;  ied = ie + ng
  jsd = js - ng;  jed = je + ng

  if ( id_ws > 0 ) then
       kmax = km + 1
  else
       kmax = km
  endif

!$omp parallel do default(shared) private(ra_x, ra_y, fx, fy, crx_adv, xfx_adv, cry_adv, yfx_adv, dz1)
  do k=1,kmax
     if ( k <= kint ) then
        do j=jsd,jed
           do i=isd,ied
              dz1(i,j) = zh(i,j,k+1) - zh(i,j,k)
           enddo
        enddo
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

        call fv_tp_2d(dz1,    crx(is,jsd,k), cry(isd,js,k), npx,  npy, hord, &
                      fx, fy, xfx(is,jsd,k), yfx(isd,js,k), area, ra_x, ra_y, nord=ndif(k), damp_c=damp(k))
        do j=js,je
           do i=is,ie
              delz(i,j,k) = (dz1(i,j)*area(i,j)+fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))   &
                          / (ra_x(i,j) + yfx(i,j,k)-yfx(i,j+1,k)) 
           enddo
        enddo

     elseif ( k == km+1 ) then
! Compute ws:
       if ( id_ws > 0 ) then
        do j=jsd,jed 
           do i=is,ie+1
              crx_adv(i,j) = crx(i,j,km) * surf_fac
              xfx_adv(i,j) = xfx(i,j,km) * surf_fac
           enddo
        enddo
        do j=js,je+1
           do i=isd,ied
              cry_adv(i,j) = cry(i,j,km) * surf_fac
              yfx_adv(i,j) = yfx(i,j,km) * surf_fac
           enddo
        enddo
        call fv_tp_2d(zs, crx_adv, cry_adv, npx, npy, hord, &
                      fx, fy, xfx_adv, yfx_adv, area, ra_x, ra_y, nord=ndif(km), damp_c=damp(km))
        do j=js,je
           do i=is,ie
              ws(i,j) = rdt*(zs(i,j) - (zs(i,j)*area(i,j)+fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))   &
                                     / (ra_x(i,j) + yfx_adv(i,j)-yfx_adv(i,j+1)))
           enddo
        enddo
       endif

     else            

      if ( k==1 ) then
        do j=jsd,jed 
           do i=is,ie+1
              crx_adv(i,j) = crx(i,j,k)
              xfx_adv(i,j) = xfx(i,j,k)
           enddo
        enddo
        do j=js,je+1
           do i=isd,ied
              cry_adv(i,j) = cry(i,j,k)
              yfx_adv(i,j) = yfx(i,j,k)
           enddo
        enddo
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

        call fv_tp_2d(zh(isd,jsd,k), crx_adv, cry_adv, npx,  npy, hord, &
                      fx, fy, xfx_adv, yfx_adv, area, ra_x, ra_y, nord=ndif(k), damp_c=damp(k))
        do j=js,je
           do i=is,ie
              zh(i,j,k) = (zh(i,j,k)*area(i,j)+fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))   &
                        / (ra_x(i,j) + yfx_adv(i,j)-yfx_adv(i,j+1))
           enddo
        enddo

        if ( k==km ) then
           if ( fill_wz ) then
               do j=js,je
                  do i=is,ie
                     delz(i,j,km) = zs(i,j) - zh(i,j,km)
                  enddo
               enddo
           else
               do j=js,je
                  do i=is,ie
                     zh(i,j,km) = max(zs(i,j)-dz_max, zh(i,j,km)) 
                     delz(i,j,km) = zs(i,j) - zh(i,j,km)
                  enddo
               enddo
           endif
        endif

     endif
  enddo

  if ( kint < km-1 ) then
!$omp parallel do default(shared)
     do k=kint+1, km-1
        do j=js,je
           do i=is,ie
              delz(i,j,k) = zh(i,j,k+1) - zh(i,j,k)
           enddo
        enddo
     enddo
  endif


  end subroutine update_dz_d


  subroutine Riem_Solver(dt,   is,   ie,   js, je, km, ng,    &
                         akap, cp,   ptop, hs, w,  delz, pt,  &
                         delp, gz,   pkc, pk, pe, peln,       &
                         da_min, diff_z0, first_call, last_call, compute_pvar)
!--------------------------------------------
! !OUTPUT PARAMETERS
! Ouput: gz: grav*height at edges
!        pe: full     hydrostatic pressure
!       pkc: non-hydrostatic "perturbation" pressure
!--------------------------------------------
   integer, intent(in):: is, ie, js, je, km, ng
   real, intent(in):: dt         ! the BIG horizontal Lagrangian time step
   real, intent(in):: akap, cp, ptop
   real, intent(in):: da_min, diff_z0
   real, intent(in):: hs(is-ng:ie+ng,js-ng:je+ng)
   logical, intent(in):: first_call, last_call, compute_pvar
   real, intent(inout):: peln(is:ie,km+1,js:je)          ! ln(pe)
   real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,km):: w, delp, pt
   real, intent(inout):: delz(is:ie,js:je,km)
   real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: pkc
   real, intent(out):: pe(is-1:ie+1,km+1,js-1:je+1)
   real, intent(out):: pk(is:ie,js:je,km+1)
   real, intent(out), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: gz
! Local:
  real, dimension(is:ie,km):: pm, dm, dz2, w2
  real :: pem(is:ie,km+1)
  real gama, rgrav, ptk, peln1
  integer i, j, k

    gama = 1./(1.-akap)
   rgrav = 1./grav
     ptk = ptop ** akap
   peln1 = log(ptop)

   damp3 = (kd3*da_min/300.) / (dt/real(m_split))

   if ( first_call .and. kd3 > near0 )  allocate ( density(is:ie,1:km,js:je) )

!$omp parallel do default(shared) private(dm, dz2, pem, pm)
   do 2000 j=js,je

      do k=1,km
         do i=is,ie
            dm(i,k) = delp(i,j,k)
         enddo
      enddo
      do i=is,ie
         pem(i,1) = ptop
      enddo
      do k=2,km+1
         do i=is,ie
            pem(i,k) = pem(i,k-1) + dm(i,k-1)
         enddo
      enddo

      if ( compute_pvar ) then
           do i=is,ie
              pk(i,j,1) = ptk
              peln(i,1,j) = peln1
           enddo
           do k=2,km+1
              do i=is,ie
                 peln(i,k,j) = log(pem(i,k))
                 pk(i,j,k) = exp( akap*peln(i,k,j) )
              enddo
           enddo
      endif

      do k=k_top,km
         do i=is,ie
! Non-hydrostatic pressure:
!      p2(k) = ( -R*dm/dz*pt2(k) )**gama
! hydrostatic pressure:
            pm(i,k) = exp( gama*log(akap*dm(i,k)/(pk(i,j,k+1)-pk(i,j,k))) )
!           pm(i,k) = dm(i,k) / (peln(i,k+1,j)-peln(i,k,j))
            dm(i,k) = dm(i,k) * rgrav
            dz2(i,k) = delz(i,j,k)
             w2(i,k) = w(i,j,k)
         enddo
      enddo

      call Riem_3D(m_split, dt, is, ie, js, je, ng, j, km, cp, gama, akap,  &
                   pkc, dm, pm, pem, w2, dz2, pt, density(is,1,j), first_call, k_top, m_riem)

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

      if ( diff_z0 > near0 ) then
            call imp_diff_w(j, is, ie, js, je, ng, km, diff_z0, dz2, w2, w)
      else
          do k=1,km
             do i=is,ie
                w(i,j,k) = w2(i,k)
             enddo
          enddo
      endif

2000  continue

  if ( last_call ) then
    if( kd3 > near0 ) deallocate ( density )

!$omp parallel do default(shared) 
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


  subroutine imp_diff_w(j, is, ie, js, je, ng, km, cd, delz, w, w3)
  integer, intent(in) :: j, is, ie, js, je, km, ng
  real, intent(in) :: cd
  real, intent(in) :: delz(is:ie, km)  ! delta-height (m)
  real, intent(in) :: w(is:ie, km)  ! vertical vel. (m/s)
  real, intent(out) :: w3(is-ng:ie+ng,js-ng:je+ng,km)
! Local:
  real, dimension(is:ie,km):: c, gam, dz, wt
  real:: bet(is:ie)
  real:: a
  integer:: i, k

     do k=2,km
        do i=is,ie
           dz(i,k) = 0.5*(delz(i,k-1)+delz(i,k))
        enddo
     enddo
     do k=1,km-1
        do i=is,ie
           c(i,k) = -cd/(dz(i,k+1)*delz(i,k))
        enddo
     enddo

! model top:
     do i=is,ie
         bet(i) = 1. - c(i,1)      ! bet(i) = b
        wt(i,1) = w(i,1) / bet(i)
     enddo

! Interior:
     do k=2,km-1
        do i=is,ie
           gam(i,k) = c(i,k-1)/bet(i)
                  a = cd/(dz(i,k)*delz(i,k))
             bet(i) = (1.+a-c(i,k)) + a*gam(i,k)
            wt(i,k) = (w(i,k) + a*wt(i,k-1)) / bet(i)
        enddo
     enddo

! Bottom: k = km
     do i=is,ie
        gam(i,km) = c(i,km-1) / bet(i)
                a = cd/(dz(i,km)*delz(i,km))
         wt(i,km) = (w(i,km) + a*wt(i,km-1))/(1. + a + (cd+cd)/delz(i,km)**2 + a*gam(i,km))
     enddo
 
     do k=km-1,1,-1
        do i=is,ie
           wt(i,k) = wt(i,k) - gam(i,k+1)*wt(i,k+1)
        enddo
     enddo

     do k=1,km
        do i=is,ie
           w3(i,j,k) = wt(i,k)
        enddo
     enddo

  end subroutine imp_diff_w



  subroutine Riem_3D(ns, bdt, is, ie, js, je, ng, j, km, cp, gama, cappa, p3, dm2,    &
                     pm2, pe2, w2, dz2, pt, den2, first_call, ktop, iad)

  integer, intent(in):: ns, is, ie, js, je, ng,  km, j
  integer, intent(in):: iad      ! time step scheme 
  integer, intent(in):: ktop     ! starting layer for non-hydrostatic dynamics
                                 ! 1: All non-hydrostatic
                                 ! 2: top sponge layer is hydrostatic
  real,    intent(in):: bdt, cp, gama, cappa
  real,    intent(in), dimension(is:ie,km):: dm2, pm2
  logical, intent(in):: first_call
  real, intent(in  ) :: pt (is-ng:ie+ng,js-ng:je+ng,km)
  real, intent(in  ) :: pe2(is:ie,km+1)
! IN/OUT:
  real, intent(inout):: den2(is:ie,km)
  real, intent(inout):: dz2(is:ie,km)
  real, intent(inout)::  w2(is:ie,km)
  real, intent(out  )::  p3(is-ng:ie+ng,js-ng:je+ng,km+1)
! --- Local 1D copyies -----
  real, dimension(km):: c2, p2, pt2, r_n, rden, dz, wm
!-------------------------------------
! Halo into the ground:
  real, dimension(km+km+1):: r_p, dm, dts
!-------------------------------------
  real, dimension(km+1):: m_bot, m_top, r_bot, r_top, time_left, pe1, pbar, wbar
  real, parameter:: dzmx = -5.0
  real    :: dt, rdt, grg, z_frac, ptmp1, den
  real    :: rcp
  real    :: seq(ns)       ! time stepping sequence
  real, parameter:: p_fac = 0.9
  real:: pm_fac
  integer :: k2(km+1)
  integer :: i, k, n, ke, kt, k0, k1, k3, kmax

  call time_sequence( iad, ns, bdt, seq )

  grg = gama * rdgas  
  rcp = 1. / cp
  rdt = 1. / bdt

  r_n(:) = 1.E40
   dm(:) = 1.E40
  r_p(:) = 1.E40
  dts(:) = 1.E40


  do k=1,km+1
     wbar(k) = 0.
     pbar(k) = 0.
  enddo

 do 6000 i=is,ie

    do k=ktop,km
       dz(k) = dz2(i,k)
       dm(k) = dm2(i,k)
       wm(k) = w2(i,k)*dm(k)
       rden(k) = -rdgas*dm(k)/dz(k)
       pt2(k) = pt(i,j,k) * rcp   ! virtual effect is included in pt
!      p2(k) = ( rden(k)*pt2(k) )**gama
       p2(k) = exp( gama*log(rden(k)*pt2(k)) )      ! full pressure
       c2(k) = sqrt( grg*p2(k)/rden(k) )
    enddo
! Initialize "density"
    if ( kd3 > near0 .and. first_call ) then
         do k=ktop,km
            den2(i,k) = -dm(k)/dz(k)
         enddo
    endif

    do k=1,km+1
       pe1(k) = 0.
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! extend dm into the ground
    kmax = 2*km + 1
    do k=km+1, km+km
       dm(k) = dm(kmax-k)
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 do 5000 n=1,ns

   dt = seq(n)
   pm_fac = -p_fac*dt

   if ( kd3 > near0 ) then
! Apply time filter to p_prime !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do k=ktop,km
         dts(k) = -dz(k) / c2(k)
         den = -dm(k)/dz(k)
         ptmp1 = dts(k) * ( p2(k) - pm2(i,k) + damp3*(den2(i,k)-den) )
         den2(i,k) = den     !  save current density
         r_p(k) = wm(k) + ptmp1
         r_n(k) = wm(k) - ptmp1
      enddo
   else
      do k=ktop,km
         dts(k) = -dz(k) / c2(k)
         ptmp1 = dts(k)*(p2(k)-pm2(i,k))
         r_p(k) = wm(k) + ptmp1
         r_n(k) = wm(k) - ptmp1
      enddo
   endif

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
            if ( k < ke-2 ) then ! original
!           if ( k < ke-1 ) then
               k1 = ke-1
               k2(k1) = k
               m_top(k1) = m_top(ke) - dm(k1)
               r_top(k1) = r_top(ke) - r_n(k1)
               time_left(k1) = time_left(ke) + dts(k1)
            endif
            m_top(ke) = m_top(ke) + z_frac*dm(k)
            r_top(ke) = r_top(ke) + z_frac*r_n(k)
            go to 444     ! All done, next level
        else
            time_left(ke) = time_left(ke) - dts(k)
            m_top(ke) = m_top(ke) + dm(k)
            r_top(ke) = r_top(ke) + r_n(k)
        endif 
     enddo
! wave from ke already left the top
     if ( ke == ktop+1 ) exit
! Finish all interfaces above ke
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Extend r_p into the ground
    do k=km+1, km+km
       kt = kmax - k       !   [km,1,-1]
       r_p(k) = -r_n(kt)   ! w changes sign
       dts(k) =  dts(kt)
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do 4000 ke=ktop,km
        kt = k2(ke)
     do k=kt,km+km
        z_frac = time_left(ke)/dts(k)
        if ( z_frac <= 1. ) then
!            partial volume
             if ( k > ke+1 ) then
                time_left(ke+1) = time_left(ke) + dts(ke)
                m_bot(ke+1) =  m_bot(ke) - dm(ke)
                r_bot(ke+1) =  r_bot(ke) - r_p(ke)
                k2(ke+1) = k
             endif
             m_bot(ke) = m_bot(ke) + z_frac*dm(k)
             r_bot(ke) = r_bot(ke) + z_frac*r_p(k)
             go to 4000      ! to next interface
!            exit
        else 
!            whole volume
             time_left(ke) = time_left(ke) - dts(k)
                 m_bot(ke) = m_bot(ke) + dm(k)
                 r_bot(ke) = r_bot(ke) + r_p(k)
        endif
     enddo
4000 continue

! Top BC:
  pbar(ktop) = 0.
  wbar(ktop) = r_bot(ktop) / m_bot(ktop)
  do k=ktop+1,km
     wbar(k) = (r_bot(k)+r_top(k)) / (m_top(k)+m_bot(k))
  enddo
!---------------------------------------
! Apply limiter to wbar from bot to top
!---------------------------------------
  do k=km,ktop+1,-1
     wbar(k) = max( wbar(k), (dz(k)-dzmx)/dt + wbar(k+1) )
  enddo

! Bot BC:
! wbar(km+1) = 0.
! pbar(km+1) = -r_top(km+1)
  do k=ktop+1,km+1
     pbar(k) =  m_top(k)*wbar(k) - r_top(k)
!-- pbar limiter: total p > 0 --------------
     pbar(k) =  max( pbar(k), pm_fac*pe2(i,k) )
! pbar limiter: ----------------------------
      pe1(k) = pe1(k) + pbar(k)
  enddo

   if ( n==ns ) then
        do k=ktop,km
           dz2(i,k) = dz(k) + dt*(wbar(k+1)-wbar(k))
           w2(i,k) = ( wm(k) + pbar(k+1) - pbar(k) ) / dm(k)
        enddo
   else
        do k=ktop,km
           wm(k) = wm(k) + pbar(k+1) - pbar(k)
           dz(k) = dz(k) + dt*(wbar(k+1)-wbar(k))
           rden(k) = -rdgas*dm(k)/dz(k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          p2(k) = (rden(k)*pt2(k))**gama
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           p2(k) = exp( gama*log(rden(k)*pt2(k)) )
           c2(k) = sqrt( grg*p2(k)/rden(k) )
        enddo
   endif

5000  continue

!---------------------------------------------------------------------
! Note: time-mean full pressure at edges and time-level-(n+1) DZ are
! used for computation of p-gradient
! Could use cell center full pressure at time-level-(n+1)
!---------------------------------------------------------------------
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

end module nh_core_mod
