module nh_core_mod
! Developer: S.-J. Lin, NOAA/GFDL
! To do list:
! include moisture effect in pt
!------------------------------

#ifdef GEOPK_CHECK
   use mpp_mod,        only: mpp_pe  
#endif
   use constants_mod,  only: rdgas, grav
   use tp_core_mod,    only: fv_tp_2d, copy_corners
   use fv_arrays_mod,  only: fv_grid_type, fv_grid_bounds_type

   implicit none
   private

   public Riem_Solver, Riem_Solver_c, update_dz_c, update_dz_d, geopk_halo_nh
   real, parameter:: dz_min = 2.

CONTAINS 

  subroutine update_dz_c(is, ie, js, je, km, ng, dt, dp0, zs, area, ut, vt, gz, ws)
! !INPUT PARAMETERS:
  integer, intent(in):: is, ie, js, je, ng, km
  real, intent(in):: dt
  real, intent(in):: dp0(km)
  real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng,km):: ut, vt
  real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng):: area
  real, intent(inout):: gz(is-ng:ie+ng,js-ng:je+ng,km+1)
  real, intent(in   ):: zs(is-ng:ie+ng, js-ng:je+ng)
  real, intent(  out):: ws(is:ie, js:je)
! Local Work array:
  real, dimension(is:ie+1,js:je  ):: xfx, fx
  real, dimension(is:ie  ,js:je+1):: yfx, fy
  real, parameter:: r14 = 1./14.
  integer  i, j, k
  real:: rdt, top_ratio, bot_ratio, int_ratio
!--------------------------------------------------------------------

  rdt = 1. / dt

  top_ratio = dp0(1 ) / (dp0(   1)+dp0(2 ))
  bot_ratio = dp0(km) / (dp0(km-1)+dp0(km))

!$omp parallel do default(shared) private(xfx, yfx, fx, fy, int_ratio)
  do 6000 k=1,km+1

     if ( k==1 ) then
        do j=js,je
           do i=is,ie+1
!             xfx(i,j) = 1.5*ut(i,j,1)-0.5*ut(i,j,2)
              xfx(i,j) = ut(i,j,1) + (ut(i,j,1)-ut(i,j,2))*top_ratio
           enddo
        enddo
        do j=js,je+1
           do i=is,ie
!             yfx(i,j) = 1.5*vt(i,j,1)-0.5*vt(i,j,2)
              yfx(i,j) = vt(i,j,1) + (vt(i,j,1)-vt(i,j,2))*top_ratio
           enddo
        enddo
     elseif ( k==km+1 ) then
! Bottom extrapolation
        do j=js,je
           do i=is,ie+1
!             xfx(i,j) = 1.5*ut(i,j,km)-0.5*ut(i,j,km-1)
              xfx(i,j) = ut(i,j,km) + (ut(i,j,km)-ut(i,j,km-1))*bot_ratio
!             xfx(i,j) = r14*(3.*ut(i,j,km-2)-13.*ut(i,j,km-1)+24.*ut(i,j,km))
!             if ( xfx(i,j)*ut(i,j,km)<0. ) xfx(i,j) = 0.
           enddo
        enddo
        do j=js,je+1
           do i=is,ie
!             yfx(i,j) = 1.5*vt(i,j,km)-0.5*vt(i,j,km-1)
              yfx(i,j) = vt(i,j,km) + (vt(i,j,km)-vt(i,j,km-1))*bot_ratio
!             yfx(i,j) = r14*(3.*vt(i,j,km-2)-13.*vt(i,j,km-1)+24.*vt(i,j,km))
!             if ( yfx(i,j)*vt(i,j,km)<0. ) yfx(i,j) = 0.
           enddo
        enddo
     else
        int_ratio = 1./(dp0(k-1)+dp0(k))
        do j=js,je
           do i=is,ie+1
!             xfx(i,j) = 0.5*(ut(i,j,k-1) + ut(i,j,k))
              xfx(i,j) = (dp0(k)*ut(i,j,k-1)+dp0(k-1)*ut(i,j,k))*int_ratio
           enddo
        enddo
        do j=js,je+1
           do i=is,ie
!             yfx(i,j) = 0.5*(vt(i,j,k-1) + vt(i,j,k))
              yfx(i,j) = (dp0(k)*vt(i,j,k-1)+dp0(k-1)*vt(i,j,k))*int_ratio
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

     if ( k==km+1 ) then
        do j=js,je
           do i=is,ie
              ws(i,j) = ( zs(i,j) - gz(i,j,k) ) * rdt
           enddo
        enddo
     endif

6000 continue

! Enforce monotonicity of height to prevent blowup
!$omp parallel do default(shared)
  do j=js, je
     do k=km, 1, -1
        do i=is, ie
           gz(i,j,k) = max( gz(i,j,k), gz(i,j,k+1) + dz_min )
        enddo
     enddo
  enddo

  end subroutine update_dz_c


  subroutine update_dz_d(ndif, damp, hord, is, ie, js, je, km, ng, npx, npy, area, rarea,   &
                         dp0, zs, zh, crx, cry, xfx, yfx, delz, ws, rdt, gridstruct, bd)

  type(fv_grid_bounds_type), intent(IN) :: bd
  integer, intent(in):: is, ie, js, je, ng, km, npx, npy
  integer, intent(in):: hord
  real, intent(in)   :: rdt
  real, intent(in)   :: dp0(km)
  real, intent(in)   :: area(is-ng:ie+ng,js-ng:je+ng)
  real, intent(in)   :: rarea(is-ng:ie+ng,js-ng:je+ng)
  real,    intent(inout):: damp(km+1)
  integer, intent(inout):: ndif(km+1)
  real, intent(in   ) ::  zs(is-ng:ie+ng,js-ng:je+ng)
  real, intent(inout) ::  zh(is-ng:ie+ng,js-ng:je+ng,km+1)
  real, intent(  out) ::delz(is-ng:ie+ng,js-ng:je+ng,km)
  real, intent(inout), dimension(is:ie+1,js-ng:je+ng,km):: crx, xfx
  real, intent(inout), dimension(is-ng:ie+ng,js:je+1,km):: cry, yfx
  real, intent(inout)   :: ws(is:ie,js:je)

  type(fv_grid_type), intent(IN) :: gridstruct
!-----------------------------------------------------
! Local array:
  real, dimension(is:   ie+1, js-ng:je+ng,km+1):: crx_adv, xfx_adv
  real, dimension(is-ng:ie+ng,js:   je+1,km+1 ):: cry_adv, yfx_adv
  real, dimension(is:ie+1,js:je  ):: fx
  real, dimension(is:ie  ,js:je+1):: fy
  real:: ra_x(is:ie,js-ng:je+ng)
  real:: ra_y(is-ng:ie+ng,js:je)
!--------------------------------------------------------------------
  integer  i, j, k, isd, ied, jsd, jed

  damp(km+1) = damp(km)
  ndif(km+1) = ndif(km)
  
  isd = is - ng;  ied = ie + ng
  jsd = js - ng;  jed = je + ng

!$omp parallel do default(shared)
  do j=jsd,jed
     call edge_profile(crx, xfx, crx_adv, xfx_adv, is,  ie+1, jsd, jed, j, km, &
                            dp0, .true., 0)
     if ( j<=je+1 .and. j>=js )      &
     call edge_profile(cry, yfx, cry_adv, yfx_adv, isd, ied,  js, je+1, j, km, &
                            dp0, .true., 0)
  enddo

!$omp parallel do default(shared) private(ra_x, ra_y, fx, fy)
  do k=1,km+1

     do j=jsd,jed
        do i=is,ie
           ra_x(i,j) = area(i,j) + xfx_adv(i,j,k) - xfx_adv(i+1,j,k)
        enddo
     enddo
     do j=js,je
        do i=isd,ied
           ra_y(i,j) = area(i,j) + yfx_adv(i,j,k) - yfx_adv(i,j+1,k)
        enddo
     enddo

     call fv_tp_2d(zh(isd,jsd,k), crx_adv(is,jsd,k), cry_adv(isd,js,k), npx,  npy, hord, &
                   fx, fy, xfx_adv(is,jsd,k), yfx_adv(isd,js,k), gridstruct, bd, ra_x, ra_y,       &
                   nord=ndif(k), damp_c=damp(k))

     do j=js,je
        do i=is,ie
           zh(i,j,k) = (zh(i,j,k)*area(i,j)+fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))   &
                     / (ra_x(i,j) + yfx_adv(i,j,k)-yfx_adv(i,j+1,k))
        enddo
     enddo

     if ( k==km+1 ) then
         do j=js,je
            do i=is,ie
               ws(i,j) = (zs(i,j) - zh(i,j,k)) * rdt
            enddo
         enddo
     endif
  enddo

! Enforce monotonicity of height to prevent blowup
!$omp parallel do default(shared)
  do j=js, je
     do k=km, 1, -1
        do i=is, ie
! Enforce monotonicity of height to prevent blowup
           zh(i,j,k) = max( zh(i,j,k), zh(i,j,k+1) + dz_min )
        enddo
     enddo
  enddo

  end subroutine update_dz_d

  subroutine Riem_Solver_c(ms,   dt,  is,   ie,   js, je, km,   ng,  &
                           akap, cp,  ptop, hs, w3,  pt,  &
                           delp, gz,  pef,  ws, p_fac, a_imp)

   integer, intent(in):: is, ie, js, je, ng, km
   integer, intent(in):: ms
   real, intent(in):: dt,  akap, cp, ptop, p_fac, a_imp
   real, intent(in):: ws(is:ie,js:je)
   real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng,km):: pt, delp
   real, intent(in)::   hs(is-ng:ie+ng,js-ng:je+ng)
   real, intent(inout):: w3(is-ng:ie+ng,js-ng:je+ng,km)
! OUTPUT PARAMETERS 
   real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: gz
   real, intent(  out), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: pef !full NH pressure
! Local:
  real, dimension(is:ie,km  ):: pm, dm, dz2, w2, pt2
  real, dimension(is:ie,km+1):: pem, peln, pe2
  real gama, rgrav, peln1, rcp
  integer i, j, k

    gama = 1./(1.-akap)
   rgrav = 1./grav
     rcp = 1. / cp

   peln1 = log(ptop)

!$omp parallel do default(shared) private(dm, dz2, w2, pt2, pe2, pem, pm, peln)
   do 2000 j=js,je

      do k=1,km
         do i=is,ie
            dm(i,k) = delp(i,j,k)
         enddo
      enddo

      do i=is,ie
         pem(i,1) = ptop
         peln(i,1) = peln1
      enddo

      do k=2,km+1
         do i=is,ie
            pem(i,k) = pem(i,k-1) + dm(i,k-1)
#ifdef GEOPK_CHECK
            if (pem(i,k) < 0.) then
               print*, mpp_pe(),i,j,k, 'PEM = ', pem(i,k), dm(i,k-1)
            endif
#endif
            peln(i,k) = log( pem(i,k) )
         enddo
      enddo

      do k=1,km
         do i=is,ie
            dz2(i,k) = gz(i,j,k+1) - gz(i,j,k)
             pm(i,k) = dm(i,k) / (peln(i,k+1)-peln(i,k))
             dm(i,k) = dm(i,k) * rgrav
             w2(i,k) = w3(i,j,k)
            pt2(i,k) = pt(i,j,k) * rcp
         enddo
      enddo

      if ( a_imp <= 0.5 ) then
           call RIM_2D(ms, dt, is, ie, km, rdgas, gama, akap, pe2, &
                        dm, pm, pem, w2, dz2, pt2, ws(is,j), .true.)
      else
           call SIMPLEST_solver(dt, is, ie, km, rdgas, gama, akap, pe2,  &
                                dm, pm, pem, w2, dz2, pt2, ws(is,j), p_fac)
      endif

      do i=is,ie
         pef(i,j,1) = ptop                     ! full pressure at top
      enddo
      do k=2,km+1
         do i=is,ie
            pef(i,j,k) = pe2(i,k) + pem(i,k)  ! add hydrostatic component
         enddo
      enddo

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



  subroutine Riem_Solver(ms, dt,   is,   ie,   js, je, km, ng,    &
                         akap, cp,   ptop, hs, w,  delz, pt,  &
                         delp, zh, gz,  ppe, pk, pe, peln, &
                         ws, da_min, diff_z0, p_fac, a_imp, last_call)
!--------------------------------------------
! !OUTPUT PARAMETERS
! Ouput: gz: grav*height at edges
!        pe: full     hydrostatic pressure
!       ppe: non-hydrostatic pressure perturbation
!--------------------------------------------
   integer, intent(in):: ms, is, ie, js, je, km, ng
   real, intent(in):: dt         ! the BIG horizontal Lagrangian time step
   real, intent(in):: akap, cp, ptop
   real, intent(in):: da_min, diff_z0, p_fac, a_imp
   real, intent(in):: hs(is-ng:ie+ng,js-ng:je+ng)
   logical, intent(in):: last_call
   real, intent(in):: ws(is:ie,js:je)
   real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: zh
   real, intent(inout):: peln(is:ie,km+1,js:je)          ! ln(pe)
   real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,km):: w, delp, pt
   real, intent(out), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: ppe ! NH pert'n
   real, intent(out):: delz(is-ng:ie+ng,js-ng:je+ng,km)
   real, intent(out):: pe(is-1:ie+1,km+1,js-1:je+1)
   real, intent(out):: pk(is:ie,js:je,km+1)
   real, intent(out), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: gz
! Local:
  real, dimension(is:ie,km):: pm, dm, dz2, w2, pt2
  real, dimension(is:ie,km+1)::pem, pe2
  real gama, rgrav, ptk, peln1, rcp
  integer i, j, k

    gama = 1./(1.-akap)
   rgrav = 1./grav
     ptk = ptop ** akap
   peln1 = log(ptop)
     rcp = 1./cp

!$omp parallel do default(shared) private(dm, dz2, w2, pt2, pem, pm, pe2)
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

      do k=1,km
         do i=is,ie
            pm(i,k) = dm(i,k) / (peln(i,k+1,j)-peln(i,k,j))
            dm(i,k) = dm(i,k) * rgrav
            dz2(i,k) = zh(i,j,k+1) - zh(i,j,k)
             w2(i,k) = w(i,j,k)
            pt2(i,k) = pt(i,j,k) * rcp
         enddo
      enddo

      if ( a_imp <= 0.5 ) then
           call RIM_2D(ms, dt, is, ie, km, rdgas, gama, akap, pe2,   &
                        dm, pm, pem, w2, dz2, pt2, ws(is,j), .false.)
      elseif ( a_imp > 0.999 ) then
           call SIMPLEST_solver(dt, is, ie, km, rdgas, gama, akap, pe2, dm,   &
                                pm, pem, w2, dz2, pt2, ws(is,j), p_fac)
      else
           call SIM_solver(dt, is, ie, km, rdgas, gama, akap, pe2, dm,  &
                                pm, pem, w2, dz2, pt2, ws(is,j), p_fac, a_imp)
      endif

      do k=1,km+1
         do i=is,ie
            ppe(i,j,k) = pe2(i,k)
         enddo
      enddo

      do k=1,km
         do i=is,ie
            delz(i,j,k) = dz2(i,k)
         enddo
      enddo

      do i=is,ie
         gz(i,j,km+1) = hs(i,j)
      enddo

      do k=km,1,-1
         do i=is,ie
            gz(i,j,k) = gz(i,j,k+1) - dz2(i,k)*grav
         enddo
      enddo

      if ( diff_z0 > 1.E-3 ) then
            call imp_diff_w(j, is, ie, js, je, ng, km, diff_z0, dz2, ws(is,j), w2, w)
      else
          do k=1,km
             do i=is,ie
                w(i,j,k) = w2(i,k)
             enddo
          enddo
      endif

2000  continue

  if ( last_call ) then

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


  subroutine imp_diff_w(j, is, ie, js, je, ng, km, cd, delz, ws, w, w3)
  integer, intent(in) :: j, is, ie, js, je, km, ng
  real, intent(in) :: cd
  real, intent(in) :: delz(is:ie, km)  ! delta-height (m)
  real, intent(in) :: w(is:ie, km)  ! vertical vel. (m/s)
  real, intent(in) :: ws(is:ie)
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

! Bottom:
     do i=is,ie
        gam(i,km) = c(i,km-1) / bet(i)
                a = cd/(dz(i,km)*delz(i,km))
         wt(i,km) = (w(i,km) + 2.*ws(i)*cd/delz(i,km)**2                        &
                  +  a*wt(i,km-1))/(1. + a + (cd+cd)/delz(i,km)**2 + a*gam(i,km))
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


  subroutine RIM_2D(ms, bdt, is, ie, km, rgas, gama, cappa, p3, &
                     dm2, pm2, pe2, w2, dz2, pt2, ws, c_core )

  integer, intent(in):: ms, is, ie, km
  real,    intent(in):: bdt, gama, cappa, rgas
  real,    intent(in), dimension(is:ie,km):: dm2, pm2
  logical, intent(in):: c_core
  real, intent(in  ) :: pt2(is:ie,km)
  real, intent(in  ) :: pe2(is:ie,km+1)
  real, intent(in  ) :: ws(is:ie)
! IN/OUT:
  real, intent(inout):: dz2(is:ie,km)
  real, intent(inout)::  w2(is:ie,km)
  real, intent(out  )::  p3(is:ie,km+1)
! Local:
  real:: ws2(is:ie)
  real, dimension(km+1):: m_bot, m_top, r_bot, r_top, pe1, pbar, wbar
  real, dimension(km):: r_hi, r_lo, dz, wm, dm, dts
  real, dimension(km):: pf1, wc, cm , pp, pt1
  real:: dt, rdt, grg, z_frac, ptmp1, rden, pf, time_left
  real:: m_surf
  integer:: i, k, n, ke, kt1, ktop
  integer:: ks0, ks1

  grg = gama * rgas  
  rdt = 1. / bdt
  dt = bdt / real(ms)

    pbar(:) = 0.
    wbar(:) = 0.

    do i=is,ie
       ws2(i) = 2.*ws(i)
    enddo

 do 6000 i=is,ie

    do k=1,km
       dz(k) = dz2(i,k)
       dm(k) = dm2(i,k)
       wm(k) = w2(i,k)*dm(k)
      pt1(k) = pt2(i,k)
    enddo

    pe1(:) = 0.
    wbar(km+1) = ws(i)

    ks0 = 1
    if ( ms > 1 .and. ms < 8 ) then
! Continuity of (pbar, wbar) is maintained
         do k=1, km
              rden = -rgas*dm(k)/dz(k)
            pf1(k) = exp( gama*log(rden*pt1(k)) )
            dts(k) = -dz(k)/sqrt(grg*pf1(k)/rden)
            if ( bdt > dts(k) ) then
                 ks0 = k-1 
                 goto 222
            endif
         enddo
         ks0 = km
222      if ( ks0==1 ) goto 244

         do k=1, ks0
            cm(k) = dm(k) / dts(k)
            wc(k) = wm(k) / dts(k)
            pp(k) = pf1(k) - pm2(i,k)
         enddo

         wbar(1) = (wc(1)+pp(1)) / cm(1)
         do k=2, ks0
            wbar(k) = (wc(k-1)+wc(k) + pp(k)-pp(k-1)) / (cm(k-1)+cm(k))
            pbar(k) = bdt*(cm(k-1)*wbar(k) - wc(k-1) + pp(k-1))
             pe1(k) = pbar(k)
         enddo

         if ( ks0 == km ) then
              pbar(km+1) = bdt*(cm(km)*wbar(km+1) - wc(km) + pp(km))
              if ( c_core ) then
                   do k=1,km
                      dz2(i,k) = dz(k) + bdt*(wbar(k+1) - wbar(k))
                   enddo
              else
                   do k=1,km
                      dz2(i,k) = dz(k) + bdt*(wbar(k+1) - wbar(k))
                       w2(i,k) = (wm(k)+pbar(k+1)-pbar(k))/dm(k)
                   enddo
              endif
              p3(i,1) = 0.
              do k=2,km+1
                 p3(i,k) = pbar(k)*rdt
              enddo
              goto 6000  ! next i
         else
              if ( c_core ) then
                do k=1, ks0-1
                   dz2(i,k) = dz(k) + bdt*(wbar(k+1) - wbar(k))
                enddo
              else
                do k=1, ks0-1
                   dz2(i,k) = dz(k) + bdt*(wbar(k+1) - wbar(k))
                    w2(i,k) = (wm(k)+pbar(k+1)-pbar(k))/dm(k)
                enddo
              endif
              pbar(ks0) = pbar(ks0) / real(ms)
         endif
    endif
244 ks1 = ks0

 do 5000 n=1, ms

   do k=ks1, km
        rden = -rgas*dm(k)/dz(k)
          pf = exp( gama*log(rden*pt1(k)) )
      dts(k) = -dz(k) /  sqrt( grg*pf/rden )
       ptmp1 = dts(k)*(pf - pm2(i,k))
      r_lo(k) = wm(k) + ptmp1
      r_hi(k) = wm(k) - ptmp1
   enddo

   ktop = ks1
   do k=ks1, km
      if( dt > dts(k) ) then
          ktop = k-1
          goto 333
      endif
   enddo
   ktop = km
333   continue

 if ( ktop >= ks1 ) then
   do k=ks1, ktop
          z_frac = dt/dts(k)
      r_bot(k  ) = z_frac*r_lo(k)
      r_top(k+1) = z_frac*r_hi(k)
      m_bot(k  ) = z_frac* dm(k)
      m_top(k+1) = m_bot(k)
   enddo
   if ( ktop == km ) goto 666
 endif

   do k=ktop+2, km+1
      m_top(k) = 0.
      r_top(k) = 0.
   enddo

   kt1 = max(1, ktop)
   do 444 ke=km+1, ktop+2, -1
      time_left = dt
     do k=ke-1, kt1, -1
        if ( time_left     > dts(k) ) then
             time_left = time_left - dts(k)
             m_top(ke) = m_top(ke) +  dm(k)
             r_top(ke) = r_top(ke) + r_hi(k)
        else
               z_frac = time_left/dts(k)
            m_top(ke) = m_top(ke) + z_frac*dm(k)
            r_top(ke) = r_top(ke) + z_frac*r_hi(k)
            go to 444     ! next level
        endif 
     enddo
444 continue

  do k=ktop+1, km
     m_bot(k) = 0.
     r_bot(k) = 0.
  enddo

  do 4000 ke=ktop+1, km
     time_left = dt
     do k=ke, km
        if ( time_left > dts(k) ) then
             time_left = time_left -  dts(k)
             m_bot(ke) = m_bot(ke) +   dm(k)
             r_bot(ke) = r_bot(ke) + r_lo(k)
        else 
                z_frac = time_left/dts(k)
             m_bot(ke) = m_bot(ke) + z_frac*  dm(k)
             r_bot(ke) = r_bot(ke) + z_frac*r_lo(k)
             go to 4000      ! next interface
        endif
     enddo
     m_surf = m_bot(ke)
     do k=km, kt1, -1
        if ( time_left > dts(k) ) then
             time_left = time_left - dts(k)
             m_bot(ke) = m_bot(ke) +  dm(k)
             r_bot(ke) = r_bot(ke) - r_hi(k)
        else
                z_frac = time_left/dts(k)
             m_bot(ke) = m_bot(ke) + z_frac*  dm(k)
             r_bot(ke) = r_bot(ke) - z_frac*r_hi(k) + (m_bot(ke)-m_surf)*ws2(i)
             go to 4000      ! next interface
        endif
     enddo
4000 continue

666  if ( ks1==1 ) wbar(1) = r_bot(1) / m_bot(1)
     do k=ks1+1, km
        wbar(k) = (r_bot(k)+r_top(k)) / (m_top(k)+m_bot(k))
     enddo
! pbar here is actually dt*pbar
     do k=ks1+1, km+1
        pbar(k) = m_top(k)*wbar(k) - r_top(k)
         pe1(k) = pe1(k) + pbar(k)
     enddo

  if ( n==ms ) then
     if ( c_core ) then
       do k=ks1, km
          dz2(i,k) = dz(k) + dt*(wbar(k+1)-wbar(k))
       enddo
     else
       do k=ks1, km
          dz2(i,k) = dz(k) + dt*(wbar(k+1)-wbar(k))
           w2(i,k) = ( wm(k) + pbar(k+1) - pbar(k) ) / dm(k)
       enddo
     endif
  else
     do k=ks1, km
        dz(k) = dz(k) + dt*(wbar(k+1)-wbar(k))
        wm(k) = wm(k) + pbar(k+1) - pbar(k)
     enddo
  endif

5000  continue
      p3(i,1) = 0.
   do k=2,km+1
      p3(i,k) = pe1(k)*rdt
   enddo

6000  continue        ! end i-loop

 end subroutine RIM_2D

 subroutine SIMPLEST_solver(dt,  is,  ie, km, rgas,  gama, cappa, pe, dm2,   &
                            pm2, pem, w2, dz2, pt2, ws, p_fac)
   integer, intent(in):: is, ie, km
   real,    intent(in):: dt, rgas, gama, cappa, p_fac
   real, intent(in   ), dimension(is:ie,km):: dm2, pt2, pm2
   real, intent(in )::  ws(is:ie)
   real, intent(in ):: pem(is:ie,km+1)
   real, intent(out)::  pe(is:ie,km+1)
   real, intent(inout), dimension(is:ie,km):: dz2, w2
! Local
   real, parameter:: r3 = 1./3.
   real, dimension(is:ie,km  ):: aa, bb, dd, w1, g_rat, gam
   real, dimension(is:ie,km+1):: pp
   real, dimension(is:ie):: p1, bet
   real t1g, rdt, capa1
   integer i, k

      t1g = gama * 2.*dt*dt
      rdt = 1. / dt
    capa1 = cappa - 1.

    do k=1,km
       do i=is, ie
          w1(i,k) = w2(i,k)
          pe(i,k) = exp(gama*log(-dm2(i,k)/dz2(i,k)*rgas*pt2(i,k))) - pm2(i,k)
       enddo
    enddo

    do k=1,km-1
       do i=is, ie
          g_rat(i,k) = dm2(i,k)/dm2(i,k+1)
             bb(i,k) = 2.*(1.+g_rat(i,k))
             dd(i,k) = 3.*(pe(i,k) + g_rat(i,k)*pe(i,k+1))
       enddo
    enddo

    do i=is, ie
         bet(i) = bb(i,1)
        pp(i,1) = 0.
        pp(i,2) = dd(i,1) / bet(i)
       bb(i,km) = 2.
       dd(i,km) = 3.*pe(i,km)
    enddo

    do k=2,km
      do i=is, ie
          gam(i,k) =  g_rat(i,k-1) / bet(i)
            bet(i) =  bb(i,k) - gam(i,k)
         pp(i,k+1) = (dd(i,k) - pp(i,k) ) / bet(i)
      enddo
    enddo

    do k=km, 2, -1
      do i=is, ie
         pp(i,k) = pp(i,k) - gam(i,k)*pp(i,k+1)
      enddo
    enddo

! Start the w-solver
    do k=2, km
       do i=is, ie
          aa(i,k) = t1g/(dz2(i,k-1)+dz2(i,k)) * (pem(i,k)+pp(i,k))
       enddo
    enddo
    do i=is, ie
       bet(i)  = dm2(i,1) - aa(i,2)
       w2(i,1) = (dm2(i,1)*w1(i,1) + dt*pp(i,2)) / bet(i)
    enddo
    do k=2,km-1
       do i=is, ie
          gam(i,k) = aa(i,k) / bet(i)
            bet(i) =  dm2(i,k) - (aa(i,k) + aa(i,k+1) + aa(i,k)*gam(i,k))
           w2(i,k) = (dm2(i,k)*w1(i,k)+dt*(pp(i,k+1)-pp(i,k))-aa(i,k)*w2(i,k-1)) / bet(i)
       enddo
    enddo
    do i=is, ie
           p1(i) = t1g/dz2(i,km)*(pem(i,km+1)+pp(i,km+1))
       gam(i,km) = aa(i,km) / bet(i)
          bet(i) =  dm2(i,km) - (aa(i,km)+p1(i) + aa(i,km)*gam(i,km))
        w2(i,km) = (dm2(i,km)*w1(i,km)+dt*(pp(i,km+1)-pp(i,km))-p1(i)*ws(i)-aa(i,km)*w2(i,km-1))/bet(i)
    enddo
    do k=km-1, 1, -1
       do i=is, ie
          w2(i,k) = w2(i,k) - gam(i,k+1)*w2(i,k+1)
       enddo
    enddo

    do i=is, ie
       pe(i,1) = 0.
    enddo
    do k=1,km
       do i=is, ie
          pe(i,k+1) = pe(i,k) + dm2(i,k)*(w2(i,k)-w1(i,k))*rdt
       enddo
    enddo

    do i=is, ie
           p1(i) = ( pe(i,km) + 2.*pe(i,km+1) )*r3
       dz2(i,km) = -dm2(i,km)*rgas*pt2(i,km)*exp(capa1*log(max(p_fac*pm2(i,km),p1(i)+pm2(i,km))))
    enddo

    do k=km-1, 1, -1
       do i=is, ie
          p1(i) = (pe(i,k) + bb(i,k)*pe(i,k+1) + g_rat(i,k)*pe(i,k+2))*r3 - g_rat(i,k)*p1(i)
          dz2(i,k) = -dm2(i,k)*rgas*pt2(i,k)*exp(capa1*log(max(p_fac*pm2(i,k),p1(i)+pm2(i,k))))
       enddo
    enddo

 end subroutine SIMPLEST_solver


 subroutine SIM_solver(dt,  is,  ie, km, rgas, gama, cappa, pe2, dm2,   &
                       pm2, pem, w2, dz2, pt2, ws, p_fac, alpha)
   integer, intent(in):: is, ie, km
   real, intent(in):: dt, rgas, gama, cappa, p_fac, alpha
   real, intent(in   ), dimension(is:ie,km):: dm2, pt2, pm2
   real, intent(in )::  ws(is:ie)
   real, intent(in ):: pem(is:ie,km+1)
   real, intent(out):: pe2(is:ie,km+1)
   real, intent(inout), dimension(is:ie,km):: dz2, w2
! Local
   real, parameter:: r3 = 1./3.
   real, dimension(is:ie,km  ):: aa, bb, dd, w1, wk, g_rat, gam
   real, dimension(is:ie,km+1):: pp
   real, dimension(is:ie):: p1, wk1, bet
   real  beta, t2, t1g, rdt, ra, capa1
   integer i, k

    beta = 1. - alpha
      ra = 1. / alpha
      t2 = beta / alpha
     t1g = gama * 2.*(alpha*dt)**2
     rdt = 1. / dt
   capa1 = cappa - 1.

   do k=1,km
      do i=is, ie
          w1(i,k) = w2(i,k)
         pe2(i,k) = exp(gama*log(-dm2(i,k)/dz2(i,k)*rgas*pt2(i,k))) - pm2(i,k)
      enddo
   enddo

   do k=1,km-1
      do i=is, ie
         g_rat(i,k) = dm2(i,k)/dm2(i,k+1)
            bb(i,k) = 2.*(1.+g_rat(i,k))
            dd(i,k) = 3.*(pe2(i,k) + g_rat(i,k)*pe2(i,k+1))
      enddo
   enddo

   do i=is, ie
       bet(i) = bb(i,1)
      pp(i,1) = 0.
      pp(i,2) = dd(i,1) / bet(i)
      bb(i,km) = 2.
      dd(i,km) = 3.*pe2(i,km)
    enddo

    do k=2,km
       do i=is, ie
          gam(i,k) =  g_rat(i,k-1) / bet(i)
            bet(i) =  bb(i,k) - gam(i,k)
         pp(i,k+1) = (dd(i,k) - pp(i,k) ) / bet(i)
      enddo
    enddo

    do k=km, 2, -1
       do i=is, ie
          pp(i,k) = pp(i,k) - gam(i,k)*pp(i,k+1)
       enddo
    enddo

    do k=1, km+1
       do i=is, ie
          pe2(i,k) = pem(i,k) + pp(i,k)
       enddo
    enddo

    do k=2, km
       do i=is, ie
          aa(i,k) = t1g/(dz2(i,k-1)+dz2(i,k))*pe2(i,k)
          wk(i,k) = t2*aa(i,k)*(w1(i,k-1)-w1(i,k))
       enddo
    enddo
    do i=is, ie
       bet(i)  =  dm2(i,1) - aa(i,2)
       w2(i,1) = (dm2(i,1)*w1(i,1) + dt*pp(i,2) + wk(i,2)) / bet(i)
    enddo
    do k=2,km-1
       do i=is, ie
          gam(i,k) = aa(i,k) / bet(i)
            bet(i) =  dm2(i,k) - (aa(i,k)+aa(i,k+1) + aa(i,k)*gam(i,k))
           w2(i,k) = (dm2(i,k)*w1(i,k) + dt*(pp(i,k+1)-pp(i,k)) + wk(i,k+1)-wk(i,k)  &
                    - aa(i,k)*w2(i,k-1)) / bet(i)
       enddo
    enddo
    do i=is, ie
          wk1(i) = t1g/dz2(i,km)*pe2(i,km+1)
       gam(i,km) = aa(i,km) / bet(i)
          bet(i) =  dm2(i,km) - (aa(i,km)+wk1(i) + aa(i,km)*gam(i,km))
        w2(i,km) = (dm2(i,km)*w1(i,km) + dt*(pp(i,km+1)-pp(i,km))-wk1(i)*ws(i)*ra - &
                     wk(i,km) + t2*wk1(i)*w1(i,km) -aa(i,km)*w2(i,km-1)) / bet(i)
    enddo
    do k=km-1, 1, -1
      do i=is, ie
         w2(i,k) = w2(i,k) - gam(i,k+1)*w2(i,k+1)
      enddo
    enddo

    do i=is, ie
       pe2(i,1) = 0.
    enddo
    do k=1,km
       do i=is, ie
          pe2(i,k+1) = pe2(i,k) + ( dm2(i,k)*(w2(i,k)-w1(i,k))*rdt   &
                                  - beta*(pp(i,k+1)-pp(i,k)) ) * ra
       enddo
    enddo

    do i=is, ie
           p1(i) = (pe2(i,km)+ 2.*pe2(i,km+1))*r3
       dz2(i,km) = -dm2(i,km)*rgas*pt2(i,km)*exp(capa1*log(max(p_fac*pm2(i,km),p1(i)+pm2(i,km))))
    enddo

    do k=km-1, 1, -1
       do i=is, ie
             p1(i) = (pe2(i,k)+bb(i,k)*pe2(i,k+1)+g_rat(i,k)*pe2(i,k+2))*r3 - g_rat(i,k)*p1(i)
          dz2(i,k) = -dm2(i,k)*rgas*pt2(i,k)*exp(capa1*log(max(p_fac*pm2(i,k),p1(i)+pm2(i,k))))
       enddo
    enddo

    do k=1, km+1
       do i=is, ie
          pe2(i,k) = pe2(i,k) + beta*(pp(i,k)-pe2(i,k))
       enddo
    enddo

 end subroutine SIM_solver


 subroutine edge_scalar(q1, qe, i1, i2, km, id)
! Optimized for wind profile reconstruction:
! Developer: S.-J. Lin, NOAA/GFDL
 integer, intent(in):: i1, i2, km
 integer, intent(in):: id                       ! 0: pp 1: wind
 real, intent(in ), dimension(i1:i2,km):: q1
 real, intent(out), dimension(i1:i2,km+1):: qe
!-----------------------------------------------------------------------
 real, parameter:: r2o3 = 2./3.
 real, parameter:: r4o3 = 4./3. 
 real  gak(km)
 real  bet
 integer i, k

!------------------------------------------------
! Optimized coding for uniform grid: SJL Apr 2007
!------------------------------------------------

   if ( id==1 ) then
      do i=i1,i2
         qe(i,1) = r4o3*q1(i,1) + r2o3*q1(i,2)
      enddo
   else
      do i=i1,i2
         qe(i,1) = 1.E30
      enddo
   endif

   gak(1) = 7./3.
   do k=2,km
      gak(k) =  1. / (4. - gak(k-1))
      do i=i1,i2
         qe(i,k) = (3.*(q1(i,k-1) + q1(i,k)) - qe(i,k-1)) * gak(k)
      enddo
   enddo

   bet = 1. / (1.5 - 3.5*gak(km))
   do i=i1,i2
      qe(i,km+1) = (4.*q1(i,km) + q1(i,km-1) - 3.5*qe(i,km)) * bet
   enddo

   do k=km,1,-1
      do i=i1,i2
         qe(i,k) = qe(i,k) - gak(k)*qe(i,k+1)
      enddo
   enddo

 end subroutine edge_scalar



 subroutine edge_profile(q1, q2, q1e, q2e, i1, i2, j1, j2, j, km, dp0, uniform_grid, limiter)
! Optimized for wind profile reconstruction:
! Developer: S.-J. Lin, NOAA/GFDL
 integer, intent(in):: i1, i2, j1, j2
 integer, intent(in):: j, km
 integer, intent(in):: limiter
 logical, intent(in):: uniform_grid
 real, intent(in):: dp0(km)
 real, intent(in),  dimension(i1:i2,j1:j2,km):: q1, q2
 real, intent(out), dimension(i1:i2,j1:j2,km+1):: q1e, q2e
!-----------------------------------------------------------------------
 real, dimension(i1:i2,km+1):: qe1, qe2, gam  ! edge values
 real  gak(km)
 real  bet, r2o3, r4o3
 real  g0, gk, xt1, xt2, a_bot
 integer i, k

 if ( uniform_grid ) then
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
 else
! Assuming grid varying in vertical only
   g0 = dp0(2) / dp0(1)
  xt1 = 2.*g0*(g0+1. )
  bet =    g0*(g0+0.5)
  do i=i1,i2
      qe1(i,1) = ( xt1*q1(i,j,1) + q1(i,j,2) ) / bet
      qe2(i,1) = ( xt1*q2(i,j,1) + q2(i,j,2) ) / bet
      gam(i,1) = ( 1. + g0*(g0+1.5) ) / bet
  enddo

  do k=2,km
     gk = dp0(k-1) / dp0(k)
     do i=i1,i2
             bet =  2. + 2.*gk - gam(i,k-1)
        qe1(i,k) = ( 3.*(q1(i,j,k-1)+gk*q1(i,j,k)) - qe1(i,k-1) ) / bet
        qe2(i,k) = ( 3.*(q2(i,j,k-1)+gk*q2(i,j,k)) - qe2(i,k-1) ) / bet
        gam(i,k) = gk / bet
     enddo
  enddo
 
  a_bot = 1. + gk*(gk+1.5)
    xt1 =   2.*gk*(gk+1.)
  do i=i1,i2
             xt2 = gk*(gk+0.5) - a_bot*gam(i,km)
     qe1(i,km+1) = ( xt1*q1(i,j,km) + q1(i,j,km-1) - a_bot*qe1(i,km) ) / xt2
     qe2(i,km+1) = ( xt1*q2(i,j,km) + q2(i,j,km-1) - a_bot*qe2(i,km) ) / xt2
  enddo

  do k=km,1,-1
     do i=i1,i2
        qe1(i,k) = qe1(i,k) - gam(i,k)*qe1(i,k+1)
        qe2(i,k) = qe2(i,k) - gam(i,k)*qe2(i,k+1)
     enddo
  enddo
 endif

!------------------
! Apply constraints
!------------------
    if ( limiter/=0 ) then   ! limit the top & bottom winds
         do i=i1,i2
! Top
            if ( q1(i,j,1)*qe1(i,1) < 0. ) qe1(i,1) = 0.
            if ( q2(i,j,1)*qe2(i,1) < 0. ) qe2(i,1) = 0.
! Surface:
            if ( q1(i,j,km)*qe1(i,km+1) < 0. ) qe1(i,km+1) = 0.
            if ( q2(i,j,km)*qe2(i,km+1) < 0. ) qe2(i,km+1) = 0.
         enddo
    endif

    do k=1,km+1
       do i=i1,i2
          q1e(i,j,k) = qe1(i,k)
          q2e(i,j,k) = qe2(i,k)
       enddo
    enddo

 end subroutine edge_profile

    subroutine geopk_halo_nh(ptop, grav, kappa, cp, delp, delz, pt, phis, pkc, gz, pk3, &
         npx, npy, npz, nested, pkc_pertn, computepk3, fullhalo, bd)

      !INPUT: delp, delz, pt
      !OUTPUT: gz, pkc, pk3 (optional)
      integer, intent(IN) :: npx, npy, npz
      logical, intent(IN) :: pkc_pertn, computepk3, fullhalo, nested
      real, intent(IN) :: ptop, kappa, cp, grav
      type(fv_grid_bounds_type), intent(IN) :: bd
      real, intent(IN) :: phis(bd%isd:bd%ied,bd%jsd:bd%jed)
      real, intent(IN),  dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz):: pt, delp, delz
      real, intent(INOUT), dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz+1):: gz, pkc, pk3

      integer :: i,j,k
      real :: gama !'gamma'
      real :: rcp
      real :: ptk, rgrav, rkap
      real :: peln1

      real, dimension(bd%isd:bd%ied, npz+1, bd%jsd:bd%jed ) :: pe, peln
      real, dimension(bd%isd:bd%ied, npz) :: gam, bb, dd, pkz
      real, dimension(bd%isd:bd%ied, npz-1) :: g_rat
      real, dimension(bd%isd:bd%ied) :: bet
      real :: pm

      integer :: ifirst, ilast, jfirst, jlast

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

      if (.not. nested) return

!!$   if (fullhalo) then
!!$      ifirst = is
!!$      if (is == 1) ifirst = isd
!!$      ilast =  ie
!!$      if (ie == npx-1) ilast  = ied
!!$      jfirst = js
!!$      if (js == 1) jfirst = jsd
!!$      jlast =  je
!!$      if (je == npy-1) jlast  = jed
!!$   else
!!$      ifirst = is
!!$      if (is == 1) ifirst = is-1
!!$      ilast =  ie
!!$      if (ie == npx-1) ilast  = ie+1
!!$      jfirst = js
!!$      if (js == 1) jfirst = js-1
!!$      jlast =  je
!!$      if (je == npy-1) jlast  = je+1
!!$   endif
      ifirst = isd
      jfirst = jsd
      ilast = ied
      jlast = jed

      !Remember we want to compute these in the HALO. Note also this routine
      !requires an appropriate

      rgrav = 1./grav
      gama = 1./(1.-kappa)
      rcp = 1./cp !not R/cp
      ptk = ptop ** kappa
      rkap = 1./kappa
      peln1 = log(ptop)

!!$      do j=jfirst,jlast
!!$
!!$         !GZ
!!$         do i=ifirst,ilast
!!$            gz(i,j,npz+1) = phis(i,j)
!!$         enddo
!!$         do k=npz,1,-1
!!$         do i=ifirst,ilast
!!$               gz(i,j,k) = gz(i,j,k+1) - delz(i,j,k)*grav
!!$         enddo
!!$         enddo
!!$         
!!$         !Hydrostatic interface pressure
!!$         do i=ifirst,ilast
!!$            pe(i,1,j) = ptop
!!$            peln(i,1,j) = peln1
!!$         enddo
!!$         do k=2,npz+1
!!$         do i=ifirst,ilast
!!$            pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
!!$            peln(i,k,j) = log(pe(i,k,j))
!!$         enddo
!!$         enddo
!!$
!!$         !Perturbation nonhydro layer-mean pressure (NOT to the kappa)
!!$         do k=1,npz
!!$         do i=ifirst,ilast
!!$            !Full p
!!$            pkz(i,k) = exp(gama*log(-delp(i,j,k)*rgrav/delz(i,j,k)*rdgas*pt(i,j,k)*rcp))
!!$            !hydro
!!$            pm = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
!!$            !Remove hydro cell-mean pressure
!!$            pkz(i,k) = pkz(i,k) - pm
!!$         enddo
!!$         enddo
!!$
!!$         !Reversible interpolation on layer NH pressure perturbation
!!$         !                 to recover  lastge NH pressure perturbation
!!$         do k=1,npz-1
!!$         do i=ifirst,ilast
!!$            g_rat(i,k) = delp(i,j,k)/delp(i,j,k+1)
!!$            bb(i,k) = 2.*(1. + g_rat(i,k))
!!$            dd(i,k) = 3.*(pkz(i,k) + g_rat(i,k)*pkz(i,k+1))
!!$         enddo
!!$         enddo
!!$
!!$         do i=ifirst,ilast
!!$            bet(i) = bb(i,1)
!!$            pkc(i,j,1) = 0.
!!$            pkc(i,j,2) = dd(i,1)/bet(i)
!!$            bb(i,npz) = 2.
!!$            dd(i,npz) = 3.*pkz(i,npz)
!!$         enddo
!!$         do k=2,npz
!!$         do i=ifirst,ilast
!!$            gam(i,k) = g_rat(i,k-1)/bet(i)
!!$            bet(i) = bb(i,k) - gam(i,k)
!!$            pkc(i,j,k+1) = (dd(i,k) - pkc(i,j,k))/bet(i)
!!$         enddo
!!$         enddo
!!$         do k=npz,2,-1
!!$         do i=ifirst,ilast
!!$            pkc(i,j,k) = pkc(i,j,k) - gam(i,k)*pkc(i,j,k+1)
!!$         enddo
!!$         enddo
!!$
!!$      enddo
!!$
!!$      do j=jfirst,jlast
!!$
!!$         if (.not. pkc_pertn) then
!!$            do k=npz+1,1,-1
!!$            do i=ifirst,ilast
!!$               pkc(i,j,k) = pkc(i,j,k) + pe(i,k,j)
!!$            enddo
!!$            enddo
!!$         endif
!!$
!!$         !pk3 if necessary
!!$         if (computepk3) then
!!$            do i=ifirst,ilast
!!$               pk3(i,j,1) = ptk
!!$            enddo
!!$            do k=2,npz+1
!!$               do i=ifirst,ilast
!!$                  pk3(i,j,k) = exp(kappa*log(pe(i,k,j)))
!!$               enddo
!!$            enddo
!!$         endif
!!$
!!$      enddo


      !NOTE: Compiler does NOT like this sort of nested-grid BC code. Is it trying to do some ugly optimization?

      if (is == 1) then

         do j=jfirst,jlast

            !GZ
            do i=ifirst,0
               gz(i,j,npz+1) = phis(i,j)
            enddo
            do k=npz,1,-1
               do i=ifirst,0
                  gz(i,j,k) = gz(i,j,k+1) - delz(i,j,k)*grav
               enddo
            enddo

            !Hydrostatic interface pressure
            do i=ifirst,0
               pe(i,1,j) = ptop
               peln(i,1,j) = peln1
            enddo
            do k=2,npz+1
               do i=ifirst,0
                  pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
                  peln(i,k,j) = log(pe(i,k,j))
               enddo
            enddo

            !Perturbation nonhydro layer-mean pressure (NOT to the kappa)
            do k=1,npz
               do i=ifirst,0
                  !Full p
                  pkz(i,k) = exp(gama*log(-delp(i,j,k)*rgrav/delz(i,j,k)*rdgas*pt(i,j,k)*rcp))
                  !hydro
                  pm = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
                  !Remove hydro cell-mean pressure
                  pkz(i,k) = pkz(i,k) - pm
               enddo
            enddo

            !Reversible interpolation on layer NH pressure perturbation
            !                 to recover  edge NH pressure perturbation
            do k=1,npz-1
               do i=ifirst,0
                  g_rat(i,k) = delp(i,j,k)/delp(i,j,k+1)
                  bb(i,k) = 2.*(1. + g_rat(i,k))
                  dd(i,k) = 3.*(pkz(i,k) + g_rat(i,k)*pkz(i,k+1))
               enddo
            enddo

            do i=ifirst,0
               bet(i) = bb(i,1)
               pkc(i,j,1) = 0.
               pkc(i,j,2) = dd(i,1)/bet(i)
               bb(i,npz) = 2.
               dd(i,npz) = 3.*pkz(i,npz)
            enddo
            do k=2,npz
               do i=ifirst,0
                  gam(i,k) = g_rat(i,k-1)/bet(i)
                  bet(i) = bb(i,k) - gam(i,k)
                  pkc(i,j,k+1) = (dd(i,k) - pkc(i,j,k))/bet(i)
               enddo
            enddo
            do k=npz,2,-1
               do i=ifirst,0
                  pkc(i,j,k) = pkc(i,j,k) - gam(i,k)*pkc(i,j,k+1)
#ifdef NHNEST_DEBUG
                  if (abs(pkc(i,j,k)) > 1.e5) then
                     print*, mpp_pe(), i,j,k, 'PKC: ', pkc(i,j,k)
                  endif
#endif
               enddo
            enddo

         enddo

         do j=jfirst,jlast

            if (.not. pkc_pertn) then
               do k=npz+1,1,-1
                  do i=ifirst,0
                     pkc(i,j,k) = pkc(i,j,k) + pe(i,k,j)
                  enddo
               enddo
            endif

            !pk3 if necessary
            if (computepk3) then
               do i=ifirst,0
                  pk3(i,j,1) = ptk
               enddo
               do k=2,npz+1
                  do i=ifirst,0
                     pk3(i,j,k) = exp(kappa*log(pe(i,k,j)))
                  enddo
               enddo
            endif

         enddo

      endif

      if (ie == npx-1) then

         do j=jfirst,jlast

            !GZ
            do i=npx,ilast
               gz(i,j,npz+1) = phis(i,j)
            enddo
            do k=npz,1,-1
               do i=npx,ilast
                  gz(i,j,k) = gz(i,j,k+1) - delz(i,j,k)*grav
               enddo
            enddo

            !Hydrostatic interface pressure
            do i=npx,ilast
               pe(i,1,j) = ptop
               peln(i,1,j) = peln1
            enddo
            do k=2,npz+1
               do i=npx,ilast
                  pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
                  peln(i,k,j) = log(pe(i,k,j))
               enddo
            enddo

            !Perturbation nonhydro layer-mean pressure (NOT to the kappa)
            do k=1,npz
               do i=npx,ilast
                  !Full p
                  pkz(i,k) = exp(gama*log(-delp(i,j,k)*rgrav/delz(i,j,k)*rdgas*pt(i,j,k)*rcp))
                  !hydro
                  pm = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
                  !Remove hydro cell-mean pressure
                  pkz(i,k) = pkz(i,k) - pm
               enddo
            enddo

            !Reversible interpolation on layer NH pressure perturbation
            !                 to recover  edge NH pressure perturbation
            do k=1,npz-1
               do i=npx,ilast
                  g_rat(i,k) = delp(i,j,k)/delp(i,j,k+1)
                  bb(i,k) = 2.*(1. + g_rat(i,k))
                  dd(i,k) = 3.*(pkz(i,k) + g_rat(i,k)*pkz(i,k+1))
               enddo
            enddo

            do i=npx,ilast
               bet(i) = bb(i,1)
               pkc(i,j,1) = 0.
               pkc(i,j,2) = dd(i,1)/bet(i)
               bb(i,npz) = 2.
               dd(i,npz) = 3.*pkz(i,npz)
            enddo
            do k=2,npz
               do i=npx,ilast
                  gam(i,k) = g_rat(i,k-1)/bet(i)
                  bet(i) = bb(i,k) - gam(i,k)
                  pkc(i,j,k+1) = (dd(i,k) - pkc(i,j,k))/bet(i)
               enddo
            enddo
            do k=npz,2,-1
               do i=npx,ilast
                  pkc(i,j,k) = pkc(i,j,k) - gam(i,k)*pkc(i,j,k+1)
               enddo
            enddo


         enddo

         do j=jfirst,jlast

            if (.not. pkc_pertn) then
               do k=npz+1,1,-1
                  do i=npx,ilast
                     pkc(i,j,k) = pkc(i,j,k) + pe(i,k,j)
                  enddo
               enddo
            endif

            !pk3 if necessary
            if (computepk3) then
               do i=npx,ilast
                  pk3(i,j,1) = ptk
               enddo
               do k=2,npz+1
                  do i=npx,ilast
                     pk3(i,j,k) = exp(kappa*log(pe(i,k,j)))
                  enddo
               enddo
            endif

         enddo

      endif

      if (js == 1) then

         do j=jfirst,0

            !GZ
            do i=ifirst,ilast
               gz(i,j,npz+1) = phis(i,j)
            enddo
            do k=npz,1,-1
               do i=ifirst,ilast
                  gz(i,j,k) = gz(i,j,k+1) - delz(i,j,k)*grav
               enddo
            enddo

            !Hydrostatic interface pressure
            do i=ifirst,ilast
               pe(i,1,j) = ptop
               peln(i,1,j) = peln1
            enddo
            do k=2,npz+1
               do i=ifirst,ilast
                  pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
                  peln(i,k,j) = log(pe(i,k,j))
               enddo
            enddo

            !Perturbation nonhydro layer-mean pressure (NOT to the kappa)
            do k=1,npz
               do i=ifirst,ilast
                  !Full p
                  pkz(i,k) = exp(gama*log(-delp(i,j,k)*rgrav/delz(i,j,k)*rdgas*pt(i,j,k)*rcp))
                  !hydro
                  pm = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
                  !Remove hydro cell-mean pressure
                  pkz(i,k) = pkz(i,k) - pm
               enddo
            enddo

            !Reversible interpolation on layer NH pressure perturbation
            !                 to recover  edge NH pressure perturbation
            do k=1,npz-1
               do i=ifirst,ilast
                  g_rat(i,k) = delp(i,j,k)/delp(i,j,k+1)
                  bb(i,k) = 2.*(1. + g_rat(i,k))
                  dd(i,k) = 3.*(pkz(i,k) + g_rat(i,k)*pkz(i,k+1))
               enddo
            enddo

            do i=ifirst,ilast
               bet(i) = bb(i,1)
               pkc(i,j,1) = 0.
               pkc(i,j,2) = dd(i,1)/bet(i)
               bb(i,npz) = 2.
               dd(i,npz) = 3.*pkz(i,npz)
            enddo
            do k=2,npz
               do i=ifirst,ilast
                  gam(i,k) = g_rat(i,k-1)/bet(i)
                  bet(i) = bb(i,k) - gam(i,k)
                  pkc(i,j,k+1) = (dd(i,k) - pkc(i,j,k))/bet(i)
               enddo
            enddo
            do k=npz,2,-1
               do i=ifirst,ilast
                  pkc(i,j,k) = pkc(i,j,k) - gam(i,k)*pkc(i,j,k+1)
#ifdef NHNEST_DEBUG
                  if (abs(pkc(i,j,k)) > 1.e5) then
                     print*, mpp_pe(), i,j,k, 'PKC: ', pkc(i,j,k)
                  endif
#endif
               enddo
            enddo

         enddo

         do j=jfirst,0

            if (.not. pkc_pertn) then
               do k=npz+1,1,-1
                  do i=ifirst,ilast
                     pkc(i,j,k) = pkc(i,j,k) + pe(i,k,j)
                  enddo
               enddo
            endif

            !pk3 if necessary
            if (computepk3) then
               do i=ifirst,ilast
                  pk3(i,j,1) = ptk
               enddo
               do k=2,npz+1
                  do i=ifirst,ilast
                     pk3(i,j,k) = exp(kappa*log(pe(i,k,j)))
                  enddo
               enddo
            endif

         enddo

      endif

      if (je == npy-1) then

         do j=npy,jlast

            !GZ
            do i=ifirst,ilast
               gz(i,j,npz+1) = phis(i,j)
            enddo
            do k=npz,1,-1
               do i=ifirst,ilast
                  gz(i,j,k) = gz(i,j,k+1) - delz(i,j,k)*grav
               enddo
            enddo

            !Hydrostatic interface pressure
            do i=ifirst,ilast
               pe(i,1,j) = ptop
               peln(i,1,j) = peln1
            enddo
            do k=2,npz+1
               do i=ifirst,ilast
                  pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
                  peln(i,k,j) = log(pe(i,k,j))
               enddo
            enddo

            !Perturbation nonhydro layer-mean pressure (NOT to the kappa)
            do k=1,npz
               do i=ifirst,ilast
                  !Full p
                  pkz(i,k) = exp(gama*log(-delp(i,j,k)*rgrav/delz(i,j,k)*rdgas*pt(i,j,k)*rcp))
                  !hydro
                  pm = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
                  !Remove hydro cell-mean pressure
                  pkz(i,k) = pkz(i,k) - pm
               enddo
            enddo

            !Reversible interpolation on layer NH pressure perturbation
            !                 to recover  lastge NH pressure perturbation
            do k=1,npz-1
               do i=ifirst,ilast
                  g_rat(i,k) = delp(i,j,k)/delp(i,j,k+1)
                  bb(i,k) = 2.*(1. + g_rat(i,k))
                  dd(i,k) = 3.*(pkz(i,k) + g_rat(i,k)*pkz(i,k+1))
               enddo
            enddo

            do i=ifirst,ilast
               bet(i) = bb(i,1)
               pkc(i,j,1) = 0.
               pkc(i,j,2) = dd(i,1)/bet(i)
               bb(i,npz) = 2.
               dd(i,npz) = 3.*pkz(i,npz)
            enddo
            do k=2,npz
               do i=ifirst,ilast
                  gam(i,k) = g_rat(i,k-1)/bet(i)
                  bet(i) = bb(i,k) - gam(i,k)
                  pkc(i,j,k+1) = (dd(i,k) - pkc(i,j,k))/bet(i)
               enddo
            enddo
            do k=npz,2,-1
               do i=ifirst,ilast
                  pkc(i,j,k) = pkc(i,j,k) - gam(i,k)*pkc(i,j,k+1)
               enddo
            enddo


         enddo

         do j=npy,jlast

            if (.not. pkc_pertn) then
               do k=npz+1,1,-1
                  do i=ifirst,ilast
                     pkc(i,j,k) = pkc(i,j,k) + pe(i,k,j)
                  enddo
               enddo
            endif

            !pk3 if necessary
            if (computepk3) then
               do i=ifirst,ilast
                  pk3(i,j,1) = ptk
               enddo
               do k=2,npz+1
                  do i=ifirst,ilast
                     pk3(i,j,k) = exp(kappa*log(pe(i,k,j)))
                  enddo
               enddo
            endif

         enddo

      endif



    end subroutine geopk_halo_nh


end module nh_core_mod
