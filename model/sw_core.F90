 module sw_core_mod

 use fv_mp_mod,         only: ng, is,js,ie,je, isd,jsd,ied,jed,  &
                              mp_corner_comm, domain
 use fv_grid_tools_mod, only: npx=>npx_g,npy=>npy_g, cosa, sina,  &
                              rdxc, rdyc, dx,dy, dxc,dyc, dxa,dya,  &
                              rdxa, rdya, area, area_c, rarea, rarea_c, rdx, rdy
 use fv_grid_tools_mod, only: grid_type
 use tp_core_mod,       only: fv_tp_2d, pert_ppm, copy_corners
 use fv_grid_utils_mod, only: edge_vect_s,edge_vect_n,edge_vect_w,edge_vect_e,  &
                              sw_corner, se_corner, ne_corner, nw_corner,       &
                              cosa_u, cosa_v, cosa_s, sina_s, sina_u, sina_v,   &
                              rsin_u, rsin_v, rsina, ec1, ec2, ew, es,          &
                              big_number, da_min_c, da_min, fC, f0, stretched_grid, &
                              rsin2, divg_u, divg_v, gnomonic_grid, cos_sg, sin_sg
 use fv_mp_mod, only: fill_corners, XDir, YDir

#ifdef SW_DYNAMICS
 use test_cases_mod,   only: test_case
#endif

 implicit none

  real, parameter:: r3 =   1./3.
  real, parameter:: t11=27./28., t12=-13./28., t13=3./7., t14=6./7., t15=3./28.
  real, parameter:: s11=11./14., s13=-13./14., s14=4./7., s15=3./14.
  real, parameter:: near_zero = 1.E-9     ! for KE limiter
!----------------------
! PPM volume mean form:
!----------------------
  real, parameter:: p1 =  7./12.     ! 0.58333333
  real, parameter:: p2 = -1./12.
!----------------------------
! 4-pt Lagrange interpolation
!----------------------------
  real, parameter:: a1 =  0.5625
  real, parameter:: a2 = -0.0625
!----------------------------------------------
! volume-conserving cubic with 2nd drv=0 at end point:
  real, parameter:: c1 = -2./14.
  real, parameter:: c2 = 11./14.
  real, parameter:: c3 =  5./14.
! 3-pt off-center intp formular:
! real, parameter:: c1 = -0.125
! real, parameter:: c2 =  0.75
! real, parameter:: c3 =  0.375
!----------------------------------------------
! scheme 2.1: perturbation form
  real, parameter:: b1 =   1./30.
  real, parameter:: b2 = -13./60.
  real, parameter:: b3 = -13./60.
  real, parameter:: b4 =  0.45
  real, parameter:: b5 = -0.05

  private
  public :: c_sw, d_sw

!---- version number -----
  character(len=128) :: version = '$Id: sw_core.F90,v 19.0 2012/01/06 19:57:54 fms Exp $'
  character(len=128) :: tagname = '$Name: siena $'

  contains

 
   subroutine c_sw(delpc, delp, ptc, pt, u,v, uc,vc, ua,va,  &
                   ut, vt, divg_d, nord, dt2, hydrostatic, dord4)
      real, intent(INOUT), dimension(isd:ied,  jsd:jed+1):: u, vc
      real, intent(INOUT), dimension(isd:ied+1,jsd:jed  ):: v, uc
      real, intent(INOUT), dimension(isd:ied, jsd:jed):: delp,  pt,  ua, va, ut, vt
      real, intent(OUT  ), dimension(isd:ied, jsd:jed):: delpc, ptc
      real, intent(OUT), dimension(isd:ied+1,jsd:jed+1):: divg_d
      integer, intent(IN) :: nord
      real,    intent(IN) :: dt2
      logical, intent(IN) :: hydrostatic
      logical, intent(IN) :: dord4
! Local:
      real, dimension(is-1:ie+1,js-1:je+1):: vort, ke
      real, dimension(is-1:ie+2,js-1:je+1):: fx, fx1
      real, dimension(is-1:ie+1,js-1:je+2):: fy, fy1
      real :: dt4
      integer :: i,j, is2, ie1
      integer iep1, jep1

      iep1 = ie+1; jep1 = je+1

      call d2a2c_vect(u, v, ua, va, uc, vc, ut, vt, dord4)  
!     call d2a2c_vect_v2(u, v, ua, va, uc, vc, ut, vt)

      if( nord > 0 ) call divergence_corner(u, v, ua, va, divg_d)

      do j=js-1,jep1
         do i=is-1,iep1+1
#ifndef TEST_UP
            if (i == 1 .OR. i == npx) then
               if (ut(i,j) > 0.) then
                  ut(i,j) = dt2*ut(i,j)*dy(i,j)*sin_sg(i,j,1) 
               else
                  ut(i,j) = dt2*ut(i,j)*dy(i,j)*sin_sg(i-1,j,3)
               end if
            else
               ut(i,j) = dt2*ut(i,j)*dy(i,j)*sina_u(i,j)
            end if
#else
            if (ut(i,j) > 0.) then
                ut(i,j) = dt2*ut(i,j)*dy(i,j)*sin_sg(i-1,j,3) 
            else
                ut(i,j) = dt2*ut(i,j)*dy(i,j)*sin_sg(i,j,1)
            end if
#endif
         enddo
      enddo
      do j=js-1,jep1+1
         do i=is-1,iep1
#ifndef TEST_UP
            if (j == 1 .OR. j == npy) then
               if (vt(i,j) > 0.) then
                  vt(i,j) = dt2*vt(i,j)*dx(i,j)*sin_sg(i,j,2) 
               else
                  vt(i,j) = dt2*vt(i,j)*dx(i,j)*sin_sg(i,j-1,4)
               end if
            else
               vt(i,j) = dt2*vt(i,j)*dx(i,j)*sina_v(i,j)
            end if
#else
            if (vt(i,j) > 0.) then
                vt(i,j) = dt2*vt(i,j)*dx(i,j)*sin_sg(i,j-1,4) 
            else
                vt(i,j) = dt2*vt(i,j)*dx(i,j)*sin_sg(i,j,  2)
            end if
#endif
         enddo
      enddo

!----------------
! Transport delp:
!----------------
! Xdir:
      if (grid_type < 3) call fill2_4corners(delp, pt, 1)
      if ( hydrostatic ) then
#ifdef SW_DYNAMICS
           do j=js-1,jep1
              do i=is-1,iep1+1      
                 if ( ut(i,j) > 0. ) then
                      fx1(i,j) = delp(i-1,j)
                 else
                      fx1(i,j) = delp(i,j)
                 endif
                 fx1(i,j) =  ut(i,j)*fx1(i,j)
              enddo
           enddo
#else
           do j=js-1,jep1
              do i=is-1,iep1+1
                 if ( ut(i,j) > 0. ) then
                      fx1(i,j) = delp(i-1,j)
                       fx(i,j) =   pt(i-1,j)
                 else
                      fx1(i,j) = delp(i,j)
                       fx(i,j) =   pt(i,j)
                 endif
                 fx1(i,j) =  ut(i,j)*fx1(i,j)
                  fx(i,j) = fx1(i,j)* fx(i,j)
              enddo
           enddo
#endif
      else
           do j=js-1,je+1
              do i=is-1,ie+2      
                 if ( ut(i,j) > 0. ) then
                      fx1(i,j) = delp(i-1,j)
                       fx(i,j) =   pt(i-1,j)
                 else
                      fx1(i,j) = delp(i,j)
                       fx(i,j) =   pt(i,j)
                 endif
                 fx1(i,j) =  ut(i,j)*fx1(i,j)
                  fx(i,j) = fx1(i,j)* fx(i,j)
              enddo
           enddo
      endif

! Ydir:
      if (grid_type < 3) call fill2_4corners(delp, pt, 2)
      if ( hydrostatic ) then
           do j=js-1,jep1+1
              do i=is-1,iep1      
                 if ( vt(i,j) > 0. ) then
                      fy1(i,j) = delp(i,j-1)
                       fy(i,j) =   pt(i,j-1)
                 else
                      fy1(i,j) = delp(i,j)
                       fy(i,j) =   pt(i,j)
                 endif
                 fy1(i,j) =  vt(i,j)*fy1(i,j)
                  fy(i,j) = fy1(i,j)* fy(i,j)
              enddo
           enddo
           do j=js-1,jep1
              do i=is-1,iep1    
                 delpc(i,j) = delp(i,j) + (fx1(i,j)-fx1(i+1,j)+fy1(i,j)-fy1(i,j+1))*rarea(i,j)
#ifdef SW_DYNAMICS
                   ptc(i,j) = pt(i,j)
#else
                   ptc(i,j) = (pt(i,j)*delp(i,j) +   &
                              (fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))*rarea(i,j))/delpc(i,j)
#endif
              enddo
           enddo
      else
           do j=js-1,je+2
              do i=is-1,ie+1      
                 if ( vt(i,j) > 0. ) then
                      fy1(i,j) = delp(i,j-1)
                       fy(i,j) =   pt(i,j-1)
                 else
                      fy1(i,j) = delp(i,j)
                       fy(i,j) =   pt(i,j)
                 endif
                 fy1(i,j) =  vt(i,j)*fy1(i,j)
                  fy(i,j) = fy1(i,j)* fy(i,j)
              enddo
           enddo
           do j=js-1,je+1
              do i=is-1,ie+1    
                 delpc(i,j) = delp(i,j) + (fx1(i,j)-fx1(i+1,j)+fy1(i,j)-fy1(i,j+1))*rarea(i,j)
                   ptc(i,j) = (pt(i,j)*delp(i,j) +   &
                              (fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))*rarea(i,j))/delpc(i,j)
              enddo
           enddo
      endif

!------------
! Compute KE:
!------------

!Since uc = u*, i.e. the covariant wind perpendicular to the face edge, if we want to compute kinetic energy we will need the true coordinate-parallel covariant wind, computed through u = uc*sina + v*cosa. 
!Use the alpha for the cell KE is being computed in.
      do j=js-1,jep1
         do i=is-1,iep1
            if ( ua(i,j) > 0. ) then
                 if ( i==1 ) then
                    ke(1,j) = uc(1,  j)*sin_sg(1,j,1)+v(1,j)*cos_sg(1,j,1)
                 elseif ( i==npx ) then
                    ke(i,j) = uc(npx,j)*sin_sg(npx,j,1)+v(npx,j)*cos_sg(npx,j,1)
                 else
                    ke(i,j) = uc(i,j)
                 endif
            else
                 if ( i==0  ) then
                    ke(0,j) = uc(1,  j)*sin_sg(0,j,3)+v(1,j)*cos_sg(0,j,3)
                 elseif ( i==(npx-1)  ) then
                    ke(i,j) = uc(npx,j)*sin_sg(npx-1,j,3)+v(npx,j)*cos_sg(npx-1,j,3)
                 else
                    ke(i,j) = uc(i+1,j)
                 endif
            endif
         enddo
      enddo
      do j=js-1,jep1
         do i=is-1,iep1
            if ( va(i,j) > 0. ) then
               if ( j==1  ) then
                  vort(i,1) = vc(i,  1)*sin_sg(i,1,2)+u(i,  1)*cos_sg(i,1,2)
               elseif ( j==npy  ) then
                  vort(i,j) = vc(i,npy)*sin_sg(i,npy,2)+u(i,npy)*cos_sg(i,npy,2)
               else
                  vort(i,j) = vc(i,j)
               endif
            else
               if ( j==0  ) then
                  vort(i,0) = vc(i,  1)*sin_sg(i,0,4)+u(i,  1)*cos_sg(i,0,4)
               elseif ( j==(npy-1) ) then
                  vort(i,j) = vc(i,npy)*sin_sg(i,npy-1,4)+u(i,npy)*cos_sg(i,npy-1,4)
               else
                  vort(i,j) = vc(i,j+1)
               endif
            endif
         enddo
      enddo

      dt4 = 0.5*dt2
      do j=js-1,jep1
         do i=is-1,iep1
            ke(i,j) = dt4*(ua(i,j)*ke(i,j) + va(i,j)*vort(i,j)) 
         enddo
      enddo

!------------------------------
! Compute circulation on C grid
!------------------------------
! To consider using true co-variant winds at face edges?
#ifdef TEST_EDGE
      do j=js-1,je+1
         do i=is,ie+1
            fx(i,j) = uc(i,j) * dxc(i,j)
         enddo
      enddo
      do j=js,je+1
         do i=is-1,ie+1
            fy(i,j) = vc(i,j) * dyc(i,j)
         enddo
      enddo
#else
      is2 = max(2,is); ie1 = min(npx-1,ie+1)
      do j=js-1,je+1
         do i=is2,ie1
            fx(i,j) = uc(i,j)*dxc(i,j)
         enddo
         !Multiplying in sin increases the error slightly.
         if(  is   ==  1 ) fx(1,  j) = uc(1,  j)*dxc(1,  j)!* &
              !0.5*( sin_sg(1,j,1) + sin_sg(0,j,3)  )
         if( (ie+1)==npx ) fx(npx,j) = uc(npx,j)*dxc(npx,j)!* &
              !0.5*( sin_sg(npx,j,1) + sin_sg(npx-1,j,3)  )
      enddo

      do j=js,je+1
         if( j==1 .or. j==npy ) then
           do i=is-1,ie+1
              fy(i,j) = vc(i,j)*dyc(i,j)!* &
                   !0.5*( sin_sg(i,j,2) + sin_sg(i,j-1,4)  )
           enddo
         else
           do i=is-1,ie+1
              fy(i,j) = vc(i,j)*dyc(i,j)
           enddo
         endif
      enddo
#endif
      do j=js,je+1
         do i=is,ie+1
            vort(i,j) =  fx(i,j-1) - fx(i,j) - fy(i-1,j) + fy(i,j)
         enddo
      enddo

! Remove the extra term at the corners:
      if ( sw_corner ) vort(1,    1) = vort(1,    1) + fy(0,   1)
      if ( se_corner ) vort(npx  ,1) = vort(npx,  1) - fy(npx, 1)
      if ( ne_corner ) vort(npx,npy) = vort(npx,npy) - fy(npx,npy)
      if ( nw_corner ) vort(1,  npy) = vort(1,  npy) + fy(0,  npy)

!----------------------------
! Compute absolute vorticity
!----------------------------
      do j=js,je+1
         do i=is,ie+1
            vort(i,j) = fC(i,j) + rarea_c(i,j) * vort(i,j)
         enddo
      enddo

!----------------------------------
! Transport absolute vorticity:
!----------------------------------
!To go from v to contravariant v at the edges, we divide by sin_sg;
! but we then must multiply by sin_sg to get the proper flux.
! These cancel, leaving us with fy1 = dt2*v at the edges.
! (For the same reason we only divide by sin instead of sin**2 in the interior)

      do j=js,je
         do i=is,iep1
            if ( i==1 .or. i==npx ) then
                 fy1(i,j) = dt2*v(i,j)
            else
                 fy1(i,j) = dt2*(v(i,j)-uc(i,j)*cosa_u(i,j))/sina_u(i,j)
            endif
            if ( fy1(i,j) > 0. ) then
                 fy(i,j) = vort(i,j)
            else
                 fy(i,j) = vort(i,j+1)
            endif
          enddo
      enddo

      do j=js,jep1
         if ( j==1 .or. j==npy ) then
            do i=is,ie
               fx1(i,j) = dt2*u(i,j)
               if ( fx1(i,j) > 0. ) then
                    fx(i,j) = vort(i,j)
               else
                    fx(i,j) = vort(i+1,j)
               endif
            enddo
         else
            do i=is,ie
               fx1(i,j) = dt2*(u(i,j)-vc(i,j)*cosa_v(i,j))/sina_v(i,j)
               if ( fx1(i,j) > 0. ) then
                    fx(i,j) = vort(i,j)
               else
                    fx(i,j) = vort(i+1,j)
               endif
            enddo
         endif
      enddo

! Update time-centered winds on the C-Grid
      do j=js,je
         do i=is,iep1
            uc(i,j) = uc(i,j) + fy1(i,j)*fy(i,j) + rdxc(i,j)*(ke(i-1,j)-ke(i,j))
         enddo
      enddo
      do j=js,jep1
         do i=is,ie
            vc(i,j) = vc(i,j) - fx1(i,j)*fx(i,j) + rdyc(i,j)*(ke(i,j-1)-ke(i,j))
         enddo
      enddo

   end subroutine c_sw



!-------------------------------------------------------------------------------
!
!     d_sw :: D-Grid Shallow Water Routine
!
   subroutine d_sw(delpc, delp,  ptc,   pt, u,  v, w, uc,vc, &
                   ua, va, divg_d, xflux, yflux, cx, cy,              &
                   crx_adv, cry_adv,  xfx_adv, yfx_adv, heat_source,    &
                   zvir, sphum, nq, q, k, km, inline_q,  &
                   dt, hord_tr, hord_mt, hord_vt, hord_tm, hord_dp, nord,   &
                   nord_v, dddmp, d2_bg, d4_bg, vtdm4, d_con, hydrostatic)

      integer, intent(IN):: hord_tr, hord_mt, hord_vt, hord_tm, hord_dp
      integer, intent(IN):: nord   ! nord=1 (del-4) or 3 (del-8)
      integer, intent(IN):: nord_v
      integer, intent(IN):: sphum, nq, k, km
      real   , intent(IN):: dt, dddmp, d2_bg, d4_bg, vtdm4, d_con
      real   , intent(IN):: zvir
      real, intent(inout):: divg_d(isd:ied+1,jsd:jed+1) ! divergence
      real, intent(INOUT), dimension(isd:ied,  jsd:jed):: delp, pt, ua, va, w
      real, intent(INOUT), dimension(isd:ied  ,jsd:jed+1):: u, vc
      real, intent(INOUT), dimension(isd:ied+1,jsd:jed  ):: v, uc
      real, intent(INOUT):: q(isd:ied,jsd:jed,km,nq)
      real, intent(OUT),   dimension(isd:ied,  jsd:jed)  :: delpc, ptc
      real, intent(OUT),   dimension(is:ie,js:je):: heat_source
! The flux capacitors:
      real, intent(INOUT):: xflux(is:ie+1,js:je  )
      real, intent(INOUT):: yflux(is:ie  ,js:je+1)
!------------------------
      real, intent(INOUT)::    cx(is:ie+1,jsd:jed  )
      real, intent(INOUT)::    cy(isd:ied,js:je+1)
      logical, intent(IN):: hydrostatic
      logical, intent(IN):: inline_q
      real, intent(OUT), dimension(is:ie+1,jsd:jed):: crx_adv, xfx_adv
      real, intent(OUT), dimension(isd:ied,js:je+1):: cry_adv, yfx_adv
! Local:
      real :: ut(isd:ied+1,jsd:jed)
      real :: vt(isd:ied,  jsd:jed+1)
!---
      real :: fx2(isd:ied+1,jsd:jed)
      real :: fy2(isd:ied,  jsd:jed+1)
      real :: dw(is:ie,js:je) !  work array
!---
      real, dimension(is:ie+1,js:je+1):: ub, vb
      real :: wk(isd:ied,jsd:jed) !  work array
      real :: ke(isd:ied+1,jsd:jed+1) !  needs this for corner_comm
      real :: vort(isd:ied,jsd:jed)     ! Vorticity
      real ::   fx(is:ie+1,js:je  )  ! 1-D X-direction Fluxes
      real ::   fy(is:ie  ,js:je+1)  ! 1-D Y-direction Fluxes
      real :: ra_x(is:ie,jsd:jed)
      real :: ra_y(isd:ied,js:je)
      real :: gx(is:ie+1,js:je  ) 
      real :: gy(is:ie  ,js:je+1)  ! work Y-dir flux array
      logical :: fill_c

      real :: dt4, dt5, dt6
      real :: damp, damp2, damp4, dd8, u2, v2, du2, dv2
      real :: damp_v, damp_t, damp_m
      integer :: i,j, is2, ie1, js2, je1, n, nt, n2, iq

      damp_v = vtdm4     ! vorticity & w 
!      if ( hydrostatic ) then
!           damp_t = 0.
!           damp_m = 0.
!      else
           damp_t = vtdm4     ! potential temperature
           damp_m = vtdm4     ! delp (mass)
!      endif

#ifdef SW_DYNAMICS 
      if ( test_case == 1 ) then
        do j=jsd,jed
           do i=is,ie+1
              xfx_adv(i,j) = dt * uc(i,j) / sina_u(i,j)
              if (xfx_adv(i,j) > 0.) then
                  crx_adv(i,j) = xfx_adv(i,j) * rdxa(i-1,j)
              else
                  crx_adv(i,j) = xfx_adv(i,j) * rdxa(i,j)
              endif
              xfx_adv(i,j) = dy(i,j)*xfx_adv(i,j)*sina_u(i,j)
           enddo
        enddo

        do j=js,je+1
           do i=isd,ied
              yfx_adv(i,j) = dt * vc(i,j) / sina_v(i,j)
              if (yfx_adv(i,j) > 0.) then
                 cry_adv(i,j) = yfx_adv(i,j) * rdya(i,j-1)
              else
                 cry_adv(i,j) = yfx_adv(i,j) * rdya(i,j)
              endif
              yfx_adv(i,j) = dx(i,j)*yfx_adv(i,j)*sina_v(i,j)
           enddo
        enddo
      else
#endif
     if ( grid_type < 3 ) then

        do j=jsd,jed
           if(j/=0 .and. j/=1 .and. j/=(npy-1) .and. j/=npy) then
             do i=is-1,ie+2
                ut(i,j) = ( uc(i,j) - 0.25 * cosa_u(i,j) *     &
                    (vc(i-1,j)+vc(i,j)+vc(i-1,j+1)+vc(i,j+1)))*rsin_u(i,j)
             enddo
           endif
        enddo
        do j=js-1,je+2
           if( j/=1 .and. j/=npy ) then
              do i=isd,ied
                 vt(i,j) = ( vc(i,j) - 0.25 * cosa_v(i,j) *     &
                    (uc(i,j-1)+uc(i+1,j-1)+uc(i,j)+uc(i+1,j)))*rsin_v(i,j)
              enddo
           endif
        enddo

! West edge:
       if ( is==1 ) then
          do j=jsd,jed
             if ( uc(1,j)*dt > 0. ) then
                ut(1,j) = uc(1,j) / sin_sg(0,j,3)
             else
                ut(1,j) = uc(1,j) / sin_sg(1,j,1)
             endif
          enddo
          do j=max(3,js), min(npy-2,je+1)
             vt(0,j) = vc(0,j) - 0.25*cosa_v(0,j)*   &
                  (ut(0,j-1)+ut(1,j-1)+ut(0,j)+ut(1,j))
             vt(1,j) = vc(1,j) - 0.25*cosa_v(1,j)*   &
                  (ut(1,j-1)+ut(2,j-1)+ut(1,j)+ut(2,j))
          enddo
       endif   ! West face

! East edge:
       if ( (ie+1)==npx ) then
          do j=jsd,jed
             if ( uc(npx,j)*dt > 0. ) then
                ut(npx,j) = uc(npx,j) / sin_sg(npx-1,j,3)
             else
                ut(npx,j) = uc(npx,j) / sin_sg(npx,j,1)
             endif
          enddo

           do j=max(3,js), min(npy-2,je+1)
              vt(npx-1,j) = vc(npx-1,j) - 0.25*cosa_v(npx-1,j)*   &
                           (ut(npx-1,j-1)+ut(npx,j-1)+ut(npx-1,j)+ut(npx,j))
              vt(npx,j) = vc(npx,j) - 0.25*cosa_v(npx,j)*   &
                         (ut(npx,j-1)+ut(npx+1,j-1)+ut(npx,j)+ut(npx+1,j))
           enddo
       endif

! South (Bottom) edge:
       if ( js==1 ) then

           do i=isd,ied
              if ( vc(i,1)*dt > 0. ) then
                   vt(i,1) = vc(i,1) / sin_sg(i,0,4)
              else
                   vt(i,1) = vc(i,1) / sin_sg(i,1,2)
              endif
           enddo

           do i=max(3,is),min(npx-2,ie+1)
              ut(i,0) = uc(i,0) - 0.25*cosa_u(i,0)*   &
                       (vt(i-1,0)+vt(i,0)+vt(i-1,1)+vt(i,1))
              ut(i,1) = uc(i,1) - 0.25*cosa_u(i,1)*   &
                       (vt(i-1,1)+vt(i,1)+vt(i-1,2)+vt(i,2))
           enddo
       endif

! North edge:
       if ( (je+1)==npy ) then
           do i=isd,ied
              if ( vc(i,npy)*dt > 0. ) then
                   vt(i,npy) = vc(i,npy) / sin_sg(i,npy-1,4)
              else
                   vt(i,npy) = vc(i,npy) / sin_sg(i,npy,2)
              endif
           enddo
           do i=max(3,is),min(npx-2,ie+1)
              ut(i,npy-1) = uc(i,npy-1) - 0.25*cosa_u(i,npy-1)*   &
                           (vt(i-1,npy-1)+vt(i,npy-1)+vt(i-1,npy)+vt(i,npy))
              ut(i,npy) = uc(i,npy) - 0.25*cosa_u(i,npy)*   &
                         (vt(i-1,npy)+vt(i,npy)+vt(i-1,npy+1)+vt(i,npy+1))
           enddo
       endif

! The following code solves a 2x2 system to get the interior parallel-to-edge uc,vc values 
! near the corners (ex: for the sw corner ut(2,1) and vt(1,2) are solved for simultaneously). 
! It then computes the halo uc, vc values so as to be consistent with the computations on 
! the facing panel. 

       !The system solved is:
       !  ut(2,1) = uc(2,1) - avg(vt)*cosa_u(2,1)
       !  vt(1,2) = vc(1,2) - avg(ut)*cosa_v(1,2)
       ! in which avg(vt) includes vt(1,2) and avg(ut) includes ut(2,1)

        if( sw_corner ) then
            damp = 1. / (1.-0.0625*cosa_u(2,0)*cosa_v(1,0))
            ut(2,0) = (uc(2,0)-0.25*cosa_u(2,0)*(vt(1,1)+vt(2,1)+vt(2,0) +vc(1,0) -   &
                      0.25*cosa_v(1,0)*(ut(1,0)+ut(1,-1)+ut(2,-1))) ) * damp
            damp = 1. / (1.-0.0625*cosa_u(0,1)*cosa_v(0,2))
            vt(0,2) = (vc(0,2)-0.25*cosa_v(0,2)*(ut(1,1)+ut(1,2)+ut(0,2)+uc(0,1) -   &
                      0.25*cosa_u(0,1)*(vt(0,1)+vt(-1,1)+vt(-1,2))) ) * damp

            damp = 1. / (1.-0.0625*cosa_u(2,1)*cosa_v(1,2))
            ut(2,1) = (uc(2,1)-0.25*cosa_u(2,1)*(vt(1,1)+vt(2,1)+vt(2,2)+vc(1,2) -   &
                      0.25*cosa_v(1,2)*(ut(1,1)+ut(1,2)+ut(2,2))) ) * damp

            vt(1,2) = (vc(1,2)-0.25*cosa_v(1,2)*(ut(1,1)+ut(1,2)+ut(2,2)+uc(2,1) -   &
                      0.25*cosa_u(2,1)*(vt(1,1)+vt(2,1)+vt(2,2))) ) * damp
!!$
!!$           !*********DEBUG CODE************
!!$           ut(2,0) = uc(2,0)
!!$           vt(0,2) = vc(0,2)
!!$           ut(2,1) = uc(2,1)
!!$           vt(1,2) = vc(1,2)



        endif

        if( se_corner ) then
            damp = 1. / (1. - 0.0625*cosa_u(npx-1,0)*cosa_v(npx-1,0))
            ut(npx-1,0) = ( uc(npx-1,0)-0.25*cosa_u(npx-1,0)*(   &
                            vt(npx-1,1)+vt(npx-2,1)+vt(npx-2,0)+vc(npx-1,0) -   &
                      0.25*cosa_v(npx-1,0)*(ut(npx,0)+ut(npx,-1)+ut(npx-1,-1))) ) * damp
            damp = 1. / (1. - 0.0625*cosa_u(npx+1,1)*cosa_v(npx,2))
            vt(npx,  2) = ( vc(npx,2)-0.25*cosa_v(npx,2)*(  &
                            ut(npx,1)+ut(npx,2)+ut(npx+1,2)+uc(npx+1,1) -   &
                      0.25*cosa_u(npx+1,1)*(vt(npx,1)+vt(npx+1,1)+vt(npx+1,2))) ) * damp

            damp = 1. / (1. - 0.0625*cosa_u(npx-1,1)*cosa_v(npx-1,2))
            ut(npx-1,1) = ( uc(npx-1,1)-0.25*cosa_u(npx-1,1)*(  &
                            vt(npx-1,1)+vt(npx-2,1)+vt(npx-2,2)+vc(npx-1,2) -   &
                      0.25*cosa_v(npx-1,2)*(ut(npx,1)+ut(npx,2)+ut(npx-1,2))) ) * damp
            vt(npx-1,2) = ( vc(npx-1,2)-0.25*cosa_v(npx-1,2)*(  &
                            ut(npx,1)+ut(npx,2)+ut(npx-1,2)+uc(npx-1,1) -   &
                      0.25*cosa_u(npx-1,1)*(vt(npx-1,1)+vt(npx-2,1)+vt(npx-2,2))) ) * damp

!!$           !*********DEBUG CODE************
!!$           ut(npx-1,0) = uc(npx-1,0)
!!$           vt(npx,2) = vc(npx,2)
!!$           ut(npx-1,1) = uc(npx-1,1)
!!$           vt(npx-1,2) = vc(npx-1,2)

           

        endif

        if( ne_corner ) then
            damp = 1. / (1. - 0.0625*cosa_u(npx-1,npy)*cosa_v(npx-1,npy+1))
            ut(npx-1,npy) = ( uc(npx-1,npy)-0.25*cosa_u(npx-1,npy)*(   &
                              vt(npx-1,npy)+vt(npx-2,npy)+vt(npx-2,npy+1)+vc(npx-1,npy+1) -   &
                0.25*cosa_v(npx-1,npy+1)*(ut(npx,npy)+ut(npx,npy+1)+ut(npx-1,npy+1))) ) * damp
            damp = 1. / (1. - 0.0625*cosa_u(npx+1,npy-1)*cosa_v(npx,npy-1))
            vt(npx,  npy-1) = ( vc(npx,npy-1)-0.25*cosa_v(npx,npy-1)*(   &
                                ut(npx,npy-1)+ut(npx,npy-2)+ut(npx+1,npy-2)+uc(npx+1,npy-1) -   &
                0.25*cosa_u(npx+1,npy-1)*(vt(npx,npy)+vt(npx+1,npy)+vt(npx+1,npy-1))) ) * damp

            damp = 1. / (1. - 0.0625*cosa_u(npx-1,npy-1)*cosa_v(npx-1,npy-1))
            ut(npx-1,npy-1) = ( uc(npx-1,npy-1)-0.25*cosa_u(npx-1,npy-1)*(  &
                                vt(npx-1,npy)+vt(npx-2,npy)+vt(npx-2,npy-1)+vc(npx-1,npy-1) -  &
                0.25*cosa_v(npx-1,npy-1)*(ut(npx,npy-1)+ut(npx,npy-2)+ut(npx-1,npy-2))) ) * damp
            vt(npx-1,npy-1) = ( vc(npx-1,npy-1)-0.25*cosa_v(npx-1,npy-1)*(  &
                                ut(npx,npy-1)+ut(npx,npy-2)+ut(npx-1,npy-2)+uc(npx-1,npy-1) -  &
                0.25*cosa_u(npx-1,npy-1)*(vt(npx-1,npy)+vt(npx-2,npy)+vt(npx-2,npy-1))) ) * damp

!!$           !*********DEBUG CODE************
!!$           ut(npx-1, npy) = uc(npx-1,npy)
!!$           vt(npx, npy-1) = vc(npx,npy-1)
!!$           ut(npx-1,npy-1) = uc(npx-1,npy-1)
!!$           vt(npx-1,npy-1) = vc(npx-1,npy-1)


        endif

        if( nw_corner ) then
            damp = 1. / (1. - 0.0625*cosa_u(2,npy)*cosa_v(1,npy+1))
            ut(2,npy) = ( uc(2,npy)-0.25*cosa_u(2,npy)*(   &
                          vt(1,npy)+vt(2,npy)+vt(2,npy+1)+vc(1,npy+1) -   &
                      0.25*cosa_v(1,npy+1)*(ut(1,npy)+ut(1,npy+1)+ut(2,npy+1))) ) * damp
            damp = 1. / (1. - 0.0625*cosa_u(0,npy-1)*cosa_v(0,npy-1))
            vt(0,npy-1) = ( vc(0,npy-1)-0.25*cosa_v(0,npy-1)*(  &
                            ut(1,npy-1)+ut(1,npy-2)+ut(0,npy-2)+uc(0,npy-1) -   &
                      0.25*cosa_u(0,npy-1)*(vt(0,npy)+vt(-1,npy)+vt(-1,npy-1))) ) * damp

            damp = 1. / (1. - 0.0625*cosa_u(2,npy-1)*cosa_v(1,npy-1))
            ut(2,npy-1) = ( uc(2,npy-1)-0.25*cosa_u(2,npy-1)*(  &
                            vt(1,npy)+vt(2,npy)+vt(2,npy-1)+vc(1,npy-1) -   &
                      0.25*cosa_v(1,npy-1)*(ut(1,npy-1)+ut(1,npy-2)+ut(2,npy-2))) ) * damp

            vt(1,npy-1) = ( vc(1,npy-1)-0.25*cosa_v(1,npy-1)*(  &
                            ut(1,npy-1)+ut(1,npy-2)+ut(2,npy-2)+uc(2,npy-1) -   &
                      0.25*cosa_u(2,npy-1)*(vt(1,npy)+vt(2,npy)+vt(2,npy-1))) ) * damp
!!$
!!$           !*********DEBUG CODE************
!!$           ut(2,npy) = uc(2,npy)
!!$           vt(0,npy-1) = vc(0,npy-1)
!!$           ut(2,npy-1) = uc(2,npy-1)
!!$           vt(1,npy-1) = vc(1,npy-1)

        endif
           

     else
! grid_type >= 3
        do j=jsd,jed
           do i=is,ie+1
              ut(i,j) =  uc(i,j)
           enddo
        enddo
        
        do j=js,je+1
           do i=isd,ied
              vt(i,j) = vc(i,j) 
           enddo
        enddo
     endif      ! end grid_type choices

        do j=jsd,jed
           do i=is,ie+1
              xfx_adv(i,j) = dt*ut(i,j)
           enddo
        enddo

        do j=js,je+1
           do i=isd,ied
              yfx_adv(i,j) = dt*vt(i,j)
           enddo
        enddo

! Explanation of the following code:
!    xfx_adv = dt*ut*dy   
!    crx_adv = dt*ut/dx



        do j=jsd,jed
           do i=is,ie+1
              if ( xfx_adv(i,j) > 0. ) then
                   crx_adv(i,j) = xfx_adv(i,j) * rdxa(i-1,j)
                   !if (i == 1 .OR. i == npx) then
                   !   xfx_adv(i,j) = dy(i,j)*xfx_adv(i,j)
                   !else
                      xfx_adv(i,j) = dy(i,j)*xfx_adv(i,j)*sin_sg(i-1,j,3)
                   !end if
              else
                   crx_adv(i,j) = xfx_adv(i,j) * rdxa(i,j)
                   !if (i == 1 .OR. i == npx) then
                   !   xfx_adv(i,j) = dy(i,j)*xfx_adv(i,j)
                   !else
                      xfx_adv(i,j) = dy(i,j)*xfx_adv(i,j)*sin_sg(i,j,1)
                   !end if
             endif
           enddo
        enddo
        do j=js,je+1
           do i=isd,ied
              if ( yfx_adv(i,j) > 0. ) then
                   cry_adv(i,j) = yfx_adv(i,j) * rdya(i,j-1)
                   !if (j == 1 .OR. j == npy) then
                   !   yfx_adv(i,j) = dx(i,j)*yfx_adv(i,j)
                   !else
                      yfx_adv(i,j) = dx(i,j)*yfx_adv(i,j)*sin_sg(i,j-1,4)
                   !end if
              else
                   cry_adv(i,j) = yfx_adv(i,j) * rdya(i,j)
                   !if (j == 1 .OR. j == npy) then
                   !   yfx_adv(i,j) = dx(i,j)*yfx_adv(i,j)
                   !else
                      yfx_adv(i,j) = dx(i,j)*yfx_adv(i,j)*sin_sg(i,j,2) 
                   !end if
              endif
           enddo
        enddo

#ifdef SW_DYNAMICS
      endif
#endif

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


      call fv_tp_2d(delp, crx_adv, cry_adv, npx, npy, hord_dp, fx, fy,  &
                    xfx_adv,yfx_adv, area, ra_x, ra_y, nord=nord_v, damp_c=damp_m)

#ifdef SW_DYNAMICS
        do j=js,je
           do i=is,ie
              delp(i,j) = delp(i,j) +    &
                         (fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))*rarea(i,j)
              ptc(i,j) = pt(i,j)
           enddo
        enddo
#else

! <<< Save the mass fluxes to the "Flux Capacitor" for tracer transport >>>
        do j=jsd,jed
            do i=is,ie+1
              cx(i,j) = cx(i,j) + crx_adv(i,j)
           enddo
        enddo       
        do j=js,je
           do i=is,ie+1
              xflux(i,j) = xflux(i,j) + fx(i,j)
           enddo
        enddo       

        do j=js,je+1
           do i=isd,ied
              cy(i,j) = cy(i,j) + cry_adv(i,j)
           enddo
           do i=is,ie
              yflux(i,j) = yflux(i,j) + fy(i,j)
           enddo
        enddo 

        do j=js,je
           do i=is,ie
              heat_source(i,j) = 0.
           enddo
        enddo

        if ( .not. hydrostatic ) then
            if ( damp_v>1.E-5 ) then
                 damp4 = (damp_v*da_min_c)**(nord_v+1)
                 call del6_vt_flux(nord_v, npx, npy, damp4, w, wk, fx2, fy2)
                do j=js,je
                   do i=is,ie
                      dw(i,j) = (fx2(i,j)-fx2(i+1,j)+fy2(i,j)-fy2(i,j+1))*rarea(i,j)
! 0.5 * [ (w+dw)**2 - w**2 ] = w*dw + 0.5*dw*dw
                    heat_source(i,j) = -d_con*dw(i,j)*(w(i,j)+0.5*dw(i,j))
                   enddo
                enddo
            endif
            call fv_tp_2d(w, crx_adv,cry_adv, npx, npy, hord_vt, gx, gy, xfx_adv, yfx_adv, &
                          area, ra_x, ra_y, mfx=fx, mfy=fy)
            do j=js,je
               do i=is,ie
                  w(i,j) = delp(i,j)*w(i,j) + (gx(i,j)-gx(i+1,j)+gy(i,j)-gy(i,j+1))*rarea(i,j)
               enddo
            enddo
        endif

     if ( inline_q ) then
        do j=jsd,jed
           do i=isd,ied
              pt(i,j) = pt(i,j)/(1.+zvir*q(i,j,k,sphum))
           enddo
        enddo
     endif

        call fv_tp_2d(pt, crx_adv,cry_adv, npx, npy, hord_tm, gx, gy,  &
                      xfx_adv,yfx_adv, area, ra_x, ra_y, mfx=fx, mfy=fy, mass=delp, nord=nord_v, damp_c=damp_t)

     if ( inline_q ) then
        do j=js,je
           do i=is,ie
                wk(i,j) = delp(i,j)
              delp(i,j) = wk(i,j) + (fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))*rarea(i,j)
              pt(i,j) = (pt(i,j)*wk(i,j) +               &
                        (gx(i,j)-gx(i+1,j)+gy(i,j)-gy(i,j+1))*rarea(i,j))/delp(i,j)
           enddo
        enddo
        do iq=1,nq
           call fv_tp_2d(q(isd,jsd,k,iq), crx_adv,cry_adv, npx, npy, hord_tr, gx, gy,  &
                         xfx_adv,yfx_adv, area, ra_x, ra_y, mfx=fx, mfy=fy)
           do j=js,je
              do i=is,ie
                 q(i,j,k,iq) = (q(i,j,k,iq)*wk(i,j) +               &
                         (gx(i,j)-gx(i+1,j)+gy(i,j)-gy(i,j+1))*rarea(i,j))/delp(i,j)
              enddo
           enddo
        enddo
        do j=js,je
           do i=is,ie
              pt(i,j) = pt(i,j)*(1.+zvir*q(i,j,k,sphum))
           enddo
        enddo
     else
        do j=js,je
           do i=is,ie
              pt(i,j) = pt(i,j)*delp(i,j) +               &
                         (gx(i,j)-gx(i+1,j)+gy(i,j)-gy(i,j+1))*rarea(i,j)
              delp(i,j) = delp(i,j) +                     &
                         (fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))*rarea(i,j)
              pt(i,j) = pt(i,j) / delp(i,j)
           enddo
        enddo
     endif

        if ( .not. hydrostatic ) then
          if ( damp_v>1.E-5 ) then
            do j=js,je
               do i=is,ie
                  w(i,j) = w(i,j) / delp(i,j) + dw(i,j)
               enddo
            enddo
          else
            do j=js,je
               do i=is,ie
                  w(i,j) = w(i,j) / delp(i,j)
               enddo
            enddo
          endif
        endif
#endif

#ifdef SW_DYNAMICS
      if (test_case > 1) then
#endif

!----------------------
! Kinetic Energy Fluxes
!----------------------
! Compute B grid contra-variant components for KE:

      dt5 = 0.5 *dt
      dt4 = 0.25*dt

      is2 = max(2,is); ie1 = min(npx-1,ie+1)
      js2 = max(2,js); je1 = min(npy-1,je+1)

      if (grid_type < 3) then

         if ( js==1 ) then
            do i=is,ie+1
               vb(i,1) = dt5*(vt(i-1,1)+vt(i,1))       ! corner values are incorrect
            enddo
         endif
         
         do j=js2,je1
            do i=is2,ie1
               vb(i,j) = dt5*(vc(i-1,j)+vc(i,j)-(uc(i,j-1)+uc(i,j))*cosa(i,j))*rsina(i,j)
            enddo
            if ( is==1 ) then
!               vb(1,j) = dt5*(vt(0,j)+vt(1,j)) 
! 2-pt extrapolation from both sides:
               vb(1,j) = dt4*(-vt(-1,j) + 3.*(vt(0,j)+vt(1,j)) - vt(2,j))
            endif
            if ( (ie+1)==npx ) then
!               vb(npx,j) = dt5*(vt(npx-1,j)+vt(npx,j))
               ! 2-pt extrapolation from both sides:
               vb(npx,j) = dt4*(-vt(npx-2,j) + 3.*(vt(npx-1,j)+vt(npx,j)) - vt(npx+1,j))
            endif
         enddo

         if ( (je+1)==npy ) then
            do i=is,ie+1
               vb(i,npy) = dt5*(vt(i-1,npy)+vt(i,npy)) ! corner values are incorrect
            enddo
         endif
         
      else
         do j=js,je+1
            do i=is,ie+1
               vb(i,j) = dt5*(vc(i-1,j)+vc(i,j))
            enddo
         enddo
      endif

      call ytp_v(vb, u, v, ub, hord_mt)

      do j=js,je+1
         do i=is,ie+1
            ke(i,j) = vb(i,j)*ub(i,j)
         enddo
      enddo

      if (grid_type < 3) then
         if ( is==1 ) then
            do j=js,je+1
               ub(1,j) = dt5*(ut(1,j-1)+ut(1,j))       ! corner values are incorrect
            enddo
         endif
         
         do j=js,je+1
            if ( j==1 .or. j==npy ) then
               do i=is2,ie1
!                  ub(i,j) = dt5*(ut(i,j-1)+ut(i,j))
! 2-pt extrapolation from both sides:
                  ub(i,j) = dt4*(-ut(i,j-2) + 3.*(ut(i,j-1)+ut(i,j)) - ut(i,j+1))
               enddo
            else
               do i=is2,ie1
                  ub(i,j) = dt5*(uc(i,j-1)+uc(i,j)-(vc(i-1,j)+vc(i,j))*cosa(i,j))*rsina(i,j)
               enddo
            endif
         enddo
         
         if ( (ie+1)==npx ) then
            do j=js,je+1
               ub(npx,j) = dt5*(ut(npx,j-1)+ut(npx,j))       ! corner values are incorrect
            enddo
         endif
         
      else
         do j=js,je+1
            do i=is,ie+1
               ub(i,j) = dt5*(uc(i,j-1)+uc(i,j))
            enddo
         enddo
      endif

      call xtp_u(ub, u, v, vb, hord_mt)

      do j=js,je+1
         do i=is,ie+1
            ke(i,j) = 0.5*(ke(i,j) + ub(i,j)*vb(i,j))
         enddo
      enddo

!-----------------------------------------
! Fix KE at the 4 corners of the face:
!-----------------------------------------
!  if ( gnomonic_grid ) then
      dt6 = dt / 6.
      if ( sw_corner ) then
           ke(1,1) = dt6*( (ut(1,1) + ut(1,0)) * u(1,1) +  &
                           (vt(1,1) + vt(0,1)) * v(1,1) +  &
                           (ut(1,1) + vt(1,1)) * u(0,1) )
      endif
      if ( se_corner ) then
           i = npx
           ke(i,1) = dt6*( (ut(i,1) + ut(i,  0)) * u(i-1,1) + &
                           (vt(i,1) + vt(i-1,1)) * v(i,  1) + &
                           (ut(i,1) - vt(i-1,1)) * u(i,  1) )
      endif
      if ( ne_corner ) then
           i = npx;      j = npy
           ke(i,j) = dt6*( (ut(i,j  ) + ut(i,j-1)) * u(i-1,j) +  &
                           (vt(i,j  ) + vt(i-1,j)) * v(i,j-1) +  &
                           (ut(i,j-1) + vt(i-1,j)) * u(i,j  )  )
      endif
      if ( nw_corner ) then
           j = npy
           ke(1,j) = dt6*( (ut(1,  j) + ut(1,j-1)) * u(1,j  ) +  &
                           (vt(1,  j) + vt(0,  j)) * v(1,j-1) +  &
                           (ut(1,j-1) - vt(1,  j)) * u(0,j  )  )
      endif
!  elseif (grid_type < 3) then
!     call mp_corner_comm(ke, npx, npy) 
!     if (sw_corner) ke(1,    1) = r3*(ke(2,      1)+ke(1,      2)+ke(0,      1))
!     if (se_corner) ke(npx,  1) = r3*(ke(npx+1,  1)+ke(npx,    2)+ke(npx-1,  1))
!     if (ne_corner) ke(npx,npy) = r3*(ke(npx+1,npy)+ke(npx,npy-1)+ke(npx-1,npy))
!     if (nw_corner) ke(1,  npy) = r3*(ke(2,    npy)+ke(1,  npy-1)+ke(0,    npy))
!  endif

! Compute vorticity:
       do j=jsd,jed+1
          do i=isd,ied
             vt(i,j) = u(i,j)*dx(i,j)
          enddo
       enddo
       do j=jsd,jed
          do i=isd,ied+1
             ut(i,j) = v(i,j)*dy(i,j)
          enddo
       enddo

! wk is "volume-mean" relative vorticity
       do j=jsd,jed
          do i=isd,ied
             wk(i,j) = rarea(i,j)*(vt(i,j)-vt(i,j+1)-ut(i,j)+ut(i+1,j))
          enddo
       enddo

!-----------------------------
! Compute divergence damping
!-----------------------------
   damp = dddmp * da_min_c

   if ( nord==0 ) then
!         area ~ dxb*dyb*sin(alpha)
      do j=js,je+1
         if ( j==1 .or. j==npy ) then
            do i=is-1,ie+1
               if (vc(i,j) > 0) then
                  ptc(i,j) = u(i,j)*dyc(i,j)*sin_sg(i,j-1,4)
               else
                  ptc(i,j) = u(i,j)*dyc(i,j)*sin_sg(i,j,2)
               end if
            enddo
         else
            do i=is-1,ie+1
               ptc(i,j) = (u(i,j)-0.5*(va(i,j-1)+va(i,j))*cosa_v(i,j))   &
                        *dyc(i,j)*sina_v(i,j)
            enddo
         endif
      enddo

      do j=js-1,je+1
         do i=is2,ie1
            vort(i,j) = (v(i,j) - 0.5*(ua(i-1,j)+ua(i,j))*cosa_u(i,j))  &
                        *dxc(i,j)*sina_u(i,j)
         enddo
         if ( is ==  1 ) then
            if (uc(1,j) > 0) then
               vort(1,  j) = v(1,  j)*dxc(1,  j)*sin_sg(0,j,3)
            else
               vort(1,  j) = v(1,  j)*dxc(1,  j)*sin_sg(1,j,1)
            end if
         end if
         if ( (ie+1)==npx ) then 
            if (uc(npx,j) > 0) then
               vort(npx,j) = v(npx,j)*dxc(npx,j)* & 
                    sin_sg(npx-1,j,3)
            else
               vort(npx,j) = v(npx,j)*dxc(npx,j)* &
                    sin_sg(npx,j,1)
            end if
         end if
      enddo

      do j=js,je+1
         do i=is,ie+1
            delpc(i,j) = vort(i,j-1) - vort(i,j) + ptc(i-1,j) - ptc(i,j)
         enddo
      enddo

! Remove the extra term at the corners:
      if (sw_corner) delpc(1,    1) = delpc(1,    1) - vort(1,    0)
      if (se_corner) delpc(npx,  1) = delpc(npx,  1) - vort(npx,  0)
      if (ne_corner) delpc(npx,npy) = delpc(npx,npy) + vort(npx,npy)
      if (nw_corner) delpc(1,  npy) = delpc(1,  npy) + vort(1,  npy)

      do j=js,je+1
         do i=is,ie+1
            delpc(i,j) = rarea_c(i,j)*delpc(i,j)
                damp = da_min_c*max(d2_bg, min(0.20, dddmp*abs(delpc(i,j)*dt)))
                vort(i,j) = damp*delpc(i,j)
                ke(i,j) = ke(i,j) + vort(i,j)
         enddo
      enddo
   else
!--------------------------
! Higher order divg damping
!--------------------------
     do j=js,je+1
        do i=is,ie+1
! Save divergence for external mode filter
           delpc(i,j) = divg_d(i,j)
        enddo
     enddo

     n2 = nord + 1
     do n=1,nord
        nt = nord-n

        fill_c = (nt/=0) .and. (grid_type<3) .and.               &
                 ( sw_corner .or. se_corner .or. ne_corner .or. nw_corner )

        if ( fill_c ) call fill_corners(divg_d, npx, npy, FILL=XDir, BGRID=.true.)
        do j=js-nt,je+1+nt
           do i=is-1-nt,ie+1+nt
              vc(i,j) = (divg_d(i+1,j)-divg_d(i,j))*divg_u(i,j)
           enddo
        enddo

        if ( fill_c ) call fill_corners(divg_d, npx, npy, FILL=YDir, BGRID=.true.)
        do j=js-1-nt,je+1+nt
           do i=is-nt,ie+1+nt
              uc(i,j) = (divg_d(i,j+1)-divg_d(i,j))*divg_v(i,j)
           enddo
        enddo

        if ( fill_c ) call fill_corners(vc, uc, npx, npy, VECTOR=.true., DGRID=.true.)
        do j=js-nt,je+1+nt
           do i=is-nt,ie+1+nt
              divg_d(i,j) = uc(i,j-1) - uc(i,j) + vc(i-1,j) - vc(i,j)
           enddo
        enddo

! Remove the extra term at the corners:
        if (sw_corner) divg_d(1,    1) = divg_d(1,    1) - uc(1,    0)
        if (se_corner) divg_d(npx,  1) = divg_d(npx,  1) - uc(npx,  0)
        if (ne_corner) divg_d(npx,npy) = divg_d(npx,npy) + uc(npx,npy)
        if (nw_corner) divg_d(1,  npy) = divg_d(1,  npy) + uc(1,  npy)

     if ( .not. stretched_grid ) then
        do j=js-nt,je+1+nt
           do i=is-nt,ie+1+nt
              divg_d(i,j) = divg_d(i,j)*rarea_c(i,j)
           enddo
        enddo
     endif

     enddo

     if ( dddmp<1.E-5) then
          vort = 0.
     else
        do j=js,je+1
           do i=is,ie+1
              vort(i,j) = abs(dt*delpc(i,j))
           enddo
        enddo
     endif

     if (stretched_grid ) then
! Stretched grid with variable damping ~ area
         dd8 = da_min * d4_bg**n2
     else
         dd8 = ( da_min_c*d4_bg )**n2
     endif

     do j=js,je+1
        do i=is,ie+1
           damp2 =  da_min_c*max(d2_bg, min(0.20, dddmp*vort(i,j)))  ! del-2
           vort(i,j) = damp2*delpc(i,j) + dd8*divg_d(i,j)
             ke(i,j) = ke(i,j) + vort(i,j)
        enddo
     enddo

   endif

   if ( d_con > 1.e-5 ) then
      do j=js,je+1
         do i=is,ie
            ub(i,j) = vort(i,j) - vort(i+1,j)
         enddo
      enddo
      do j=js,je
         do i=is,ie+1
            vb(i,j) = vort(i,j) - vort(i,j+1)
         enddo
      enddo
   endif

! Vorticity transport
    do j=jsd,jed
       do i=isd,ied
          vort(i,j) = wk(i,j) + f0(i,j)
       enddo
    enddo

    call fv_tp_2d(vort, crx_adv, cry_adv, npx, npy, hord_vt, fx, fy, &
                  xfx_adv,yfx_adv, area, ra_x, ra_y)
    do j=js,je+1
       do i=is,ie
          u(i,j) = vt(i,j) + ke(i,j) - ke(i+1,j) + fy(i,j)
       enddo
    enddo
    do j=js,je
       do i=is,ie+1
          v(i,j) = ut(i,j) + ke(i,j) - ke(i,j+1) - fx(i,j)
       enddo
    enddo

!--------------------------------------------------------
! damping applied to relative vorticity (wk):
   if ( damp_v>1.E-5 ) then
        damp4 = (damp_v*da_min_c)**(nord_v+1)
        call del6_vt_flux(nord_v, npx, npy, damp4, wk, vort, ut, vt)
   endif

   if ( d_con > 1.e-5 ) then
      do j=js,je+1
         do i=is,ie
            ub(i,j) = (ub(i,j) + vt(i,j))*rdx(i,j)
            fy(i,j) =  u(i,j)*rdx(i,j)
            gy(i,j) = fy(i,j)*ub(i,j)
         enddo
      enddo
      do j=js,je
         do i=is,ie+1
            vb(i,j) = (vb(i,j) - ut(i,j))*rdy(i,j)
            fx(i,j) =  v(i,j)*rdy(i,j)
            gx(i,j) = fx(i,j)*vb(i,j)
         enddo
      enddo
!----------------------------------
! Heating due to damping:
!----------------------------------
      damp = 0.25*d_con
      do j=js,je
         do i=is,ie
            u2 = fy(i,j) + fy(i,j+1)
           du2 = ub(i,j) + ub(i,j+1)
            v2 = fx(i,j) + fx(i+1,j)
           dv2 = vb(i,j) + vb(i+1,j)
! Total energy conserving:
! Convert lost KE due to divergence damping to "heat"
         heat_source(i,j) = delp(i,j)*(heat_source(i,j) - damp*rsin2(i,j)*( &
                  (ub(i,j)**2 + ub(i,j+1)**2 + vb(i,j)**2 + vb(i+1,j)**2)  &
                              + 2.*(gy(i,j)+gy(i,j+1)+gx(i,j)+gx(i+1,j))   &
                              - cosa_s(i,j)*(u2*dv2 + v2*du2 + du2*dv2)) )
         enddo
      enddo
   endif

! Add diffusive fluxes to the momentum equation:
   if ( damp_v>1.E-5 ) then
      do j=js,je+1
         do i=is,ie
            u(i,j) = u(i,j) + vt(i,j)
         enddo
      enddo
      do j=js,je
         do i=is,ie+1
            v(i,j) = v(i,j) - ut(i,j)
         enddo
      enddo
   endif

#ifdef SW_DYNAMICS
      endif ! test_case
#endif

 end subroutine d_sw


 subroutine del6_vt_flux(nord, npx, npy, damp, q, d2, fx2, fy2)
! Del-nord damping for the relative vorticity
! nord must be <= 2
!------------------
! nord = 0:   del-2
! nord = 1:   del-4
! nord = 2:   del-6
!------------------
   integer, intent(in):: nord, npx, npy
   real, intent(in):: damp
   real, intent(inout):: q(isd:ied, jsd:jed)  ! rel. vorticity ghosted on input
! Work arrays:
   real, intent(out):: d2(isd:ied, jsd:jed)
   real, intent(out):: fx2(isd:ied+1,jsd:jed), fy2(isd:ied,jsd:jed+1)
   integer i,j, nt, n, i1, i2, j1, j2

   i1 = is-1-nord;    i2 = ie+1+nord
   j1 = js-1-nord;    j2 = je+1+nord

   do j=j1, j2
      do i=i1, i2
         d2(i,j) = damp*q(i,j)
      enddo
   enddo

   if( nord>0 ) call copy_corners(d2, npx, npy, 1)
   do j=js-nord,je+nord
      do i=is-nord,ie+nord+1
!        fx2(i,j) = dy(i,j)*sina_u(i,j)*(d2(i-1,j)-d2(i,j))*rdxc(i,j)
         fx2(i,j) = 0.5*(sin_sg(i-1,j,3)+sin_sg(i,j,1))*dy(i,j)*(d2(i-1,j)-d2(i,j))*rdxc(i,j)
      enddo
   enddo

   if( nord>0 ) call copy_corners(d2, npx, npy, 2)
   do j=js-nord,je+nord+1
      do i=is-nord,ie+nord
!        fy2(i,j) = dx(i,j)*sina_v(i,j)*(d2(i,j-1)-d2(i,j))*rdyc(i,j)
         fy2(i,j) = 0.5*(sin_sg(i,j-1,4)+sin_sg(i,j,2))*dx(i,j)*(d2(i,j-1)-d2(i,j))*rdyc(i,j)
      enddo
   enddo

   if ( nord>0 ) then
   do n=1, nord
      nt = nord-n
      do j=js-nt-1,je+nt+1
         do i=is-nt-1,ie+nt+1
            d2(i,j) = (fx2(i,j)-fx2(i+1,j)+fy2(i,j)-fy2(i,j+1))*rarea(i,j)
         enddo
      enddo

      call copy_corners(d2, npx, npy, 1)

      do j=js-nt,je+nt
         do i=is-nt,ie+nt+1
!           fx2(i,j) = dy(i,j)*sina_u(i,j)*(d2(i,j)-d2(i-1,j))*rdxc(i,j)
            fx2(i,j) = 0.5*(sin_sg(i-1,j,3)+sin_sg(i,j,1))*dy(i,j)*(d2(i,j)-d2(i-1,j))*rdxc(i,j)
         enddo
      enddo

      call copy_corners(d2, npx, npy, 2)

      do j=js-nt,je+nt+1
         do i=is-nt,ie+nt
!           fy2(i,j) = dx(i,j)*sina_v(i,j)*(d2(i,j)-d2(i,j-1))*rdyc(i,j)
            fy2(i,j) = 0.5*(sin_sg(i,j-1,4)+sin_sg(i,j,2))*dx(i,j)*(d2(i,j)-d2(i,j-1))*rdyc(i,j)
         enddo
      enddo
   enddo
   endif

 end subroutine del6_vt_flux



 subroutine divergence_corner(u, v, ua, va, divg_d)
 real, intent(in),  dimension(isd:ied,  jsd:jed+1):: u
 real, intent(in),  dimension(isd:ied+1,jsd:jed  ):: v
 real, intent(in),  dimension(isd:ied,jsd:jed):: ua, va
 real, intent(out), dimension(isd:ied+1,jsd:jed+1):: divg_d
! local
 real uf(is-2:ie+2,js-1:je+2)
 real vf(is-1:ie+2,js-2:je+2)
 integer i,j
 integer is2, ie1

  is2 = max(2,is); ie1 = min(npx-1,ie+1)

  if (grid_type==4) then
        do j=js-1,je+2
           do i=is-2,ie+2
              uf(i,j) = u(i,j)*dyc(i,j)
           enddo
        enddo
        do j=js-2,je+2
           do i=is-1,ie+2
              vf(i,j) = v(i,j)*dxc(i,j)
           enddo
        enddo
        do j=js-1,je+2
           do i=is-1,ie+2
              divg_d(i,j) = rarea_c(i,j)*(vf(i,j-1)-vf(i,j)+uf(i-1,j)-uf(i,j))
           enddo
        enddo
  else
#ifdef USE_UPWIND
    do j=js,je+1
       if ( j==1 .or. j==npy ) then
            do i=is-1,ie+1
               if ( va(i,j) + va(i,j-1) > 0. ) then 
                   uf(i,j) = u(i,j)*dyc(i,j) * sin_sg(i,j-1,4)
               else
                   uf(i,j) = u(i,j)*dyc(i,j) * sin_sg(i,j,2)
               end if
            enddo
       else
            do i=is-1,ie+1
               uf(i,j) = (u(i,j)-0.5*(va(i,j-1)+va(i,j))*cosa_v(i,j))   &
                        *dyc(i,j)*sina_v(i,j)
            enddo
       endif
    enddo

    do j=js-1,je+1
       do i=is2,ie1
          vf(i,j) = (v(i,j) - 0.5*(ua(i-1,j)+ua(i,j))*cosa_u(i,j))  &
                    *dxc(i,j)*sina_u(i,j)
       enddo
       if (  is == 1 ) then
          if (ua(1,j) + ua(0,j) > 0.) then
             vf(1,j) = v(1,j)*dxc(1,j)*sin_sg(0,j,3)
          else
             vf(1,j) = v(1,j)*dxc(1,j)*sin_sg(1,j,1)
          end if
       end if
       if ( (ie+1)==npx ) then
          if (ua(npx-1,j) + ua(npx,j) > 0.) then
             vf(npx,j) = v(npx,j)*dxc(npx,j)*sin_sg(npx-1,j,3)
          else
             vf(npx,j) = v(npx,j)*dxc(npx,j)*sin_sg(npx,j,1)
          end if
       end if
    enddo
#else
!     9---4---8
!     |       |
!     1   5   3
!     |       |
!     6---2---7
    do j=js,je+1
       if ( j==1 .or. j==npy ) then
         do i=is-1,ie+1
            uf(i,j) = u(i,j)*dyc(i,j)*0.5*(sin_sg(i,j-1,4)+sin_sg(i,j,2))
         enddo
       else
         do i=is-1,ie+1
            uf(i,j) = (u(i,j)-0.25*(va(i,j-1)+va(i,j))*(cos_sg(i,j-1,4)+cos_sg(i,j,2)))   &
                                        * dyc(i,j)*0.5*(sin_sg(i,j-1,4)+sin_sg(i,j,2))
         enddo
       endif
    enddo

    do j=js-1,je+1
       do i=is2,ie1
          vf(i,j) = (v(i,j) - 0.25*(ua(i-1,j)+ua(i,j))*(cos_sg(i-1,j,3)+cos_sg(i,j,1)))  &
                                         *dxc(i,j)*0.5*(sin_sg(i-1,j,3)+sin_sg(i,j,1))
       enddo
       if (  is   ==  1 ) vf(1,  j) = v(1,  j)*dxc(1,  j)*0.5*(sin_sg(0,j,3)+sin_sg(1,j,1))
       if ( (ie+1)==npx ) vf(npx,j) = v(npx,j)*dxc(npx,j)*0.5*(sin_sg(npx-1,j,3)+sin_sg(npx,j,1))
    enddo
#endif

    do j=js,je+1
       do i=is,ie+1
          divg_d(i,j) = vf(i,j-1) - vf(i,j) + uf(i-1,j) - uf(i,j)
       enddo
    enddo

! Remove the extra term at the corners:
    if (sw_corner) divg_d(1,    1) = divg_d(1,    1) - vf(1,    0)
    if (se_corner) divg_d(npx,  1) = divg_d(npx,  1) - vf(npx,  0)
    if (ne_corner) divg_d(npx,npy) = divg_d(npx,npy) + vf(npx,npy)
    if (nw_corner) divg_d(1,  npy) = divg_d(1,  npy) + vf(1,  npy)

    do j=js,je+1
       do i=is,ie+1
          divg_d(i,j) = rarea_c(i,j)*divg_d(i,j)
       enddo
    enddo

  endif

 end subroutine divergence_corner


 subroutine xtp_u(c, u, v, flux, iord)

 real, INTENT(IN)  ::   u(isd:ied,jsd:jed+1)
 real, INTENT(IN)  ::   v(isd:ied+1,jsd:jed)
 real, INTENT(IN)  ::   c(is:ie+1,js:je+1)
 real, INTENT(out):: flux(is:ie+1,js:je+1)
 integer, INTENT(IN) :: iord
! Local
 logical extm(is-2:ie+2)
 real al(is-1:ie+2), dm(is-2:ie+2)
 real bl(is-1:ie+1)
 real br(is-1:ie+1)
 real dq(is-3:ie+2)
 real dl, dr, xt, pmp, lac, cfl
 real pmp_1, lac_1, pmp_2, lac_2
 real x0, x1, x0L, x0R
 integer i, j

 select case ( iord )

 case (1)

     do j=js,je+1
        do i=is,ie+1
           if( c(i,j)>0. ) then
               flux(i,j) = u(i-1,j)
           else
               flux(i,j) = u(i,j)
           endif
        enddo
     enddo

 case (2)

     do j=js,je+1

        do i=is-2,ie+2
              xt = 0.25*(u(i+1,j) - u(i-1,j))
           dm(i) = sign(min(abs(xt), max(u(i-1,j), u(i,j), u(i+1,j)) - u(i,j),  &
                            u(i,j) - min(u(i-1,j), u(i,j), u(i+1,j))), xt)
        enddo

! Fix slopes near edges:
      if (grid_type < 3) then
        if ( is==1 ) then
           if ( j==1 .or. j==npy ) then
              dm(0) = 0.
              dm(1) = 0.
           else
              x0L = 0.5*((2.*dx(0,j)+dx(-1,j))*(u(0,j))   &
                 - dx(0,j)*(u(-1,j)))/(dx(0,j)+dx(-1,j))
              x0R = 0.5*((2.*dx(1,j)+dx(2,j))*(u(1,j))   &
                 - dx(1,j)*(u(2,j)))/(dx(1,j)+dx(2,j))
              x0 = x0L + x0R
              x1 = s15*u(0,j) + s11*u(-1,j) + s14*dm(-1)
!          dm(0) = u(0,j) - x1
           dm(0) = 0.5*(x0 - x1)
           dm(0) = sign(min(abs(dm(0)), max(u(0,j), x0, x1) - u(0,j),   &
                               u(0,j) - min(u(0,j), x0, x1)), dm(0))
              x1 = s15*u(1,j) + s11*u(2,j) - s14*dm(2)
!          dm(1) = x1 - u(1,j)
           dm(1) = 0.5*(x1 - x0)
           dm(1) = sign(min(abs(dm(1)), max(u(1,j), x0, x1) - u(1,j),   &
                               u(1,j) - min(u(1,j), x0, x1)), dm(1))
           endif
        endif

        if ( (ie+1)==npx ) then
           if ( j==1 .or. j==npy ) then
              dm(npx-1) = 0.
              dm(npx  ) = 0.
           else
              x0L = 0.5*( (2.*dx(npx-1,j)+dx(npx-2,j))*(u(npx-1,j))  &
                - dx(npx-1,j)*(u(npx-2,j)))/(dx(npx-1,j)+dx(npx-2,j))
              x0R = 0.5*( (2.*dx(npx,j)+dx(npx+1,j))*(u(npx,j))  &
                - dx(npx,j)*(u(npx+1,j)))/(dx(npx,j)+dx(npx+1,j))
              x0 = x0L + x0R
              x1 = s15*u(npx-1,j) + s11*u(npx-2,j) + s14*dm(npx-2)
!          dm(npx-1) = u(npx-1,j) - x1
           dm(npx-1) = 0.5*(x0 - x1)
           dm(npx-1) = sign(min(abs(dm(npx-1)), max(u(npx-1,j), x0, x1) - u(npx-1,j),  &
                                   u(npx-1,j) - min(u(npx-1,j), x0, x1)), dm(npx-1))
                x1 = s15*u(npx,j) + s11*u(npx+1,j) - s14*dm(npx+1)
!          dm(npx) = x1 - u(npx,j)
           dm(npx) = 0.5*(x1 - x0)
           dm(npx) = sign(min(abs(dm(npx)), max(u(npx,j), x0, x1) - u(npx,j),   &
                                 u(npx,j) - min(u(npx,j), x0, x1)), dm(npx))
           endif
        endif
      endif

       do i=is,ie+1
          if( c(i,j)>0. ) then
             flux(i,j) = u(i-1,j) + (1.-c(i,j)*rdx(i-1,j))*dm(i-1)
          else
             flux(i,j) = u(i,  j) - (1.+c(i,j)*rdx(i,  j))*dm(i)
          endif
       enddo
     enddo

 case (4)

     do j=js,je+1

        do i=is-2,ie+2
           xt = 0.25*(u(i+1,j) - u(i-1,j))
           dm(i) = sign(min(abs(xt), max(u(i-1,j), u(i,j), u(i+1,j)) - u(i,j),  &
                            u(i,j) - min(u(i-1,j), u(i,j), u(i+1,j))), xt)
        enddo

        do i=max(3,is-1),min(npx-2,ie+2)
           al(i) = 0.5*(u(i-1,j)+u(i,j)) + r3*(dm(i-1) - dm(i))
        enddo

! Fix slopes near edges:
      if (grid_type < 3) then
        if ( is==1 ) then
           x0L = 0.5*((2.*dx(0,j)+dx(-1,j))*(u(0,j))   &
                - dx(0,j)*(u(-1,j)))/(dx(0,j)+dx(-1,j))
           x0R = 0.5*((2.*dx(1,j)+dx(2,j))*(u(1,j))   &
                - dx(1,j)*(u(2,j)))/(dx(1,j)+dx(2,j))
           x0 = x0L + x0R
          if ( j==1 .or. j==npy ) then
              dm(0) = 0.
              dm(1) = 0.
           al(0) = 0.5*(u(-1,j)+u(0,j)) + r3*dm(-1)
           al(1) = x0
           al(2) = 0.5*(u(1,j)+u(2,j)) - r3*dm(2)
          else
           x1 = s15*u(1,j) + s11*u(2,j) - s14*dm(2)
           dm(1) = 0.5*(x1 - x0)
!          dm(1) = sign(min(abs(dm(1)), max(u(1,j), x0, x1) - u(1,j),   &
!                              u(1,j) - min(u(1,j), x0, x1)), dm(1))
              x1 = s15*u(0,j) + s11*u(-1,j) + s14*dm(-1)
           dm(0) = 0.5*(x0 - x1)
!          dm(0) = sign(min(abs(dm(0)), max(u(0,j), x0, x1) - u(0,j),   &
!                              u(0,j) - min(u(0,j), x0, x1)), dm(0))
           al(0) = 0.5*(u(-1,j)+u(0,j)) + r3*(dm(-1)-dm(0))
           al(1) = x0
           al(2) = 0.5*(u(1,j)+u(2,j)) + r3*(dm(1)-dm(2))
          endif
        endif

        if ( (ie+1)==npx ) then
           x0L = 0.5*( (2.*dx(npx-1,j)+dx(npx-2,j))*(u(npx-1,j))  &
                - dx(npx-1,j)*(u(npx-2,j)))/(dx(npx-1,j)+dx(npx-2,j))
           x0R = 0.5*( (2.*dx(npx,j)+dx(npx+1,j))*(u(npx,j))  &
                - dx(npx,j)*(u(npx+1,j)))/(dx(npx,j)+dx(npx+1,j))
           x0 = x0L + x0R
           if ( j==1 .or. j==npy ) then
              dm(npx-1) = 0.
              dm(npx  ) = 0.
              al(npx-1) = 0.5*(u(npx-2,j)+u(npx-1,j)) + r3*dm(npx-2)
              al(npx  ) = x0
              al(npx+1) = 0.5*(u(npx,j)+u(npx+1,j)) - r3*dm(npx+1)
           else
              x1 = s15*u(npx-1,j) + s11*u(npx-2,j) + s14*dm(npx-2)
           dm(npx-1) = 0.5*(x0 - x1)
!          dm(npx-1) = sign(min(abs(dm(npx-1)), max(u(npx-1,j), x0, x1) - u(npx-1,j),  &
!                                  u(npx-1,j) - min(u(npx-1,j), x0, x1)), dm(npx-1))
                x1 = s15*u(npx,j) + s11*u(npx+1,j) - s14*dm(npx+1)
           dm(npx) = 0.5*(x1 - x0)
!          dm(npx) = sign(min(abs(dm(npx)), max(u(npx,j), x0, x1) - u(npx,j),   &
!                                u(npx,j) - min(u(npx,j), x0, x1)), dm(npx))
           al(npx-1) = 0.5*(u(npx-2,j)+u(npx-1,j)) + r3*(dm(npx-2) - dm(npx-1))
           al(npx  ) = x0
           al(npx+1) = 0.5*(u(npx,j)+u(npx+1,j)) + r3*(dm(npx) - dm(npx+1))
           endif
        endif
      endif

        do i=is,ie+1
          if( c(i,j)>0. ) then
             xt = 2.*dm(i-1)
             dl = sign(min(abs(xt), abs(al(i-1)-u(i-1,j))), xt)
             dr = sign(min(abs(xt), abs(al(i  )-u(i-1,j))), xt)
             cfl = c(i,j) * rdx(i-1,j)
             flux(i,j) = u(i-1,j) + (1.-cfl)*(dr + cfl*(dl-dr))
          else
             xt = 2.*dm(i)
             dl = sign(min(abs(xt), abs(al(i  )-u(i,j))), xt)
             dr = sign(min(abs(xt), abs(al(i+1)-u(i,j))), xt)
             cfl = c(i,j) * rdx(i,j)
             flux(i,j) = u(i,j) - (1.+cfl)*(dl + cfl*(dl-dr))
          endif
        enddo
     enddo

 case (6)

     do j=js,je+1

        do i=is-3,ie+2
           dq(i) = u(i+1,j) - u(i,j)
        enddo

        do i=is-2, ie+2
           if ( dq(i-1)*dq(i) > 0. ) then
                extm(i) = .false.
           else
                extm(i) = .true.
           endif
        enddo
        do i=max(3,is-1),min(npx-3,ie+1)
!!!!          if ( extm(i) .and. (extm(i-1) .or. extm(i+1)) ) then
           if ( extm(i-1) .and. extm(i) .and. extm(i+1) ) then
                bl(i) = 0.
                br(i) = 0.
           else
                bl(i) = b5*u(i-2,j) + b4*u(i-1,j) + b3*u(i,j) + b2*u(i+1,j) + b1*u(i+2,j)
                br(i) = b1*u(i-2,j) + b2*u(i-1,j) + b3*u(i,j) + b4*u(i+1,j) + b5*u(i+2,j)
           endif
        enddo

        if (grid_type < 3) then
        if ( is==1 ) then
             br(2) = p1*(u(2,j)+u(3,j)) + p2*(u(1,j)+u(4,j)) - u(2,j)
                xt = c3*u(1,j) + c2*u(2,j) + c1*u(3,j)
             bl(2) = xt - u(2,j)
             if( j==1 .or. j==npy ) then
                 bl(0) = 0.   ! out
                 br(0) = 0.   ! edge
                 bl(1) = 0.   ! edge
                 br(1) = 0.   ! in
             else
             br(1) = xt - u(1,j)
              x0L = 0.5*((2.*dx(0,j)+dx(-1,j))*(u(0,j))   &
                 - dx(0,j)*(u(-1,j)))/(dx(0,j)+dx(-1,j))
              x0R = 0.5*((2.*dx(1,j)+dx(2,j))*(u(1,j))   &
                 - dx(1,j)*(u(2,j)))/(dx(1,j)+dx(2,j))
             xt = x0L + x0R
             bl(1) = xt - u(1,j)
             br(0) = xt - u(0,j)
                xt = c1*u(-2,j) + c2*u(-1,j) + c3*u(0,j)
             bl(0) = xt - u(0,j)
             endif
        endif

        if ( (ie+1)==npx ) then
             bl(npx-2) = p1*(u(npx-2,j)+u(npx-3,j)) + p2*(u(npx-4,j)+u(npx-1,j)) - u(npx-2,j)
             xt = c1*u(npx-3,j) + c2*u(npx-2,j) + c3*u(npx-1,j)
             br(npx-2) = xt - u(npx-2,j)
             if( j==1 .or. j==npy ) then
                 bl(npx-1) = 0.  ! in
                 br(npx-1) = 0.  ! edge
                 bl(npx  ) = 0.  ! edge
                 br(npx  ) = 0.  ! out
             else
                 bl(npx-1) = xt - u(npx-1,j)
                 x0L = 0.5*( (2.*dx(npx-1,j)+dx(npx-2,j))*(u(npx-1,j))  &
                      - dx(npx-1,j)*(u(npx-2,j)))/(dx(npx-1,j)+dx(npx-2,j))
                 x0R = 0.5*( (2.*dx(npx,j)+dx(npx+1,j))*(u(npx,j))  &
                      - dx(npx,j)*(u(npx+1,j)))/(dx(npx,j)+dx(npx+1,j))
                 xt = x0L + x0R
                 br(npx-1) = xt - u(npx-1,j)
                 bl(npx  ) = xt - u(npx  ,j)
                      xt = c3*u(npx,j) + c2*u(npx+1,j) + c1*u(npx+2,j)
                 br(npx) = xt - u(npx,j)
             endif
        endif
        endif

        do i=is,ie+1
           if( c(i,j)>0. ) then
                     cfl = c(i,j)*rdx(i-1,j)
               flux(i,j) = u(i-1,j) + (1.-cfl)*(br(i-1)-cfl*(bl(i-1)+br(i-1)))
           else
                     cfl = c(i,j)*rdx(i,j)
               flux(i,j) = u(i,  j) + (1.+cfl)*(bl(i  )+cfl*(bl(i  )+br(i  )))
           endif
        enddo
     enddo

 case default
 ! iord = 8, 9, 10

     do j=js,je+1
        do i=is-2,ie+2
           xt = 0.25*(u(i+1,j) - u(i-1,j))
           dm(i) = sign(min(abs(xt), max(u(i-1,j), u(i,j), u(i+1,j)) - u(i,j),  &
                            u(i,j) - min(u(i-1,j), u(i,j), u(i+1,j))), xt)
        enddo
        do i=is-3,ie+2
           dq(i) = u(i+1,j) - u(i,j)
        enddo

        if (grid_type < 3) then

           do i=max(3,is-1),min(npx-2,ie+2)
              al(i) = 0.5*(u(i-1,j)+u(i,j)) + r3*(dm(i-1) - dm(i))
           enddo

! Perturbation form:
           if( iord==8 ) then
             do i=max(3,is-1),min(npx-3,ie+1)
                xt = 2.*dm(i)
                bl(i) = -sign(min(abs(xt), abs(al(i  )-u(i,j))), xt)
                br(i) =  sign(min(abs(xt), abs(al(i+1)-u(i,j))), xt)
             enddo
           elseif( iord==9 ) then
             do i=max(3,is-1),min(npx-3,ie+1)
              pmp_1 = -2.*dq(i)
              lac_1 = pmp_1 + 1.5*dq(i+1)
              bl(i) = min(max(0., pmp_1, lac_1), max(al(i  )-u(i,j), min(0., pmp_1, lac_1)))
              pmp_2 = 2.*dq(i-1)
              lac_2 = pmp_2 - 1.5*dq(i-2)
              br(i) = min(max(0., pmp_2, lac_2), max(al(i+1)-u(i,j), min(0., pmp_2, lac_2)))
           enddo
         elseif( iord==10 ) then
           do i=max(3,is-1), min(npx-3,ie+1)
              bl(i) = al(i  ) - u(i,j)
              br(i) = al(i+1) - u(i,j)
!             if ( abs(dm(i-1))+abs(dm(i))+abs(dm(i+1)) < near_zero ) then
              if ( abs(dm(i)) < near_zero ) then
                if ( abs(dm(i-1))+abs(dm(i+1)) < near_zero ) then
! 2-delta-x structure detected within 3 cells
                   bl(i) = 0.
                   br(i) = 0.
                endif
              elseif( abs(3.*(bl(i)+br(i))) > abs(bl(i)-br(i)) ) then
                   pmp_1 = -2.*dq(i)
                   lac_1 = pmp_1 + 1.5*dq(i+1)
                   bl(i) = min(max(0., pmp_1, lac_1), max(bl(i), min(0., pmp_1, lac_1)))
                   pmp_2 = 2.*dq(i-1)
                   lac_2 = pmp_2 - 1.5*dq(i-2)
                   br(i) = min(max(0., pmp_2, lac_2), max(br(i), min(0., pmp_2, lac_2)))
              endif
             enddo
           else
! un-limited:
             do i=max(3,is-1),min(npx-3,ie+1)
                bl(i) = al(i  ) - u(i,j)
                br(i) = al(i+1) - u(i,j)
             enddo
           endif

!--------------
! fix the edges
!--------------
           if ( is==1 ) then
              br(2) = al(3) - u(2,j)
              xt = s15*u(1,j) + s11*u(2,j) - s14*dm(2)
              bl(2) = xt - u(2,j)
              br(1) = xt - u(1,j)
              if( j==1 .or. j==npy ) then
                 bl(0) = 0.   ! out
                 br(0) = 0.   ! edge
                 bl(1) = 0.   ! edge
                 br(1) = 0.   ! in
              else
                 bl(0) = s14*dm(-1) - s11*dq(-1)
!---------------------------------------------------------------
#ifdef ONE_SIDE
                 xt = t14*u(0,j) + t12*u(-1,j) + t15*u(-2,j)
                 br(0) = 2.*xt - u(0,j)
                 xt = t14*u(1,j) + t12*u(2,j) + t15*u(3,j)
                 bl(1) = 2.*xt - u(1,j)
#else
                 x0L = 0.5*((2.*dx(0,j)+dx(-1,j))*(u(0,j))   &
                      - dx(0,j)*(u(-1,j)))/(dx(0,j)+dx(-1,j))
                 x0R = 0.5*((2.*dx(1,j)+dx(2,j))*(u(1,j))   &
                      - dx(1,j)*(u(2,j)))/(dx(1,j)+dx(2,j))
                 xt = x0L + x0R
                 br(0) = xt - u(0,j)
                 bl(1) = xt - u(1,j)
!                 br(0) = xt - 0.5*(v(1,j-1)+v(1,j))*cosa(1,j) - u(0,j)
!                 bl(1) = xt + 0.5*(v(1,j-1)+v(1,j))*cosa(1,j) - u(1,j)
#endif
!---------------------------------------------------------------
                 endif
                 call pert_ppm(1, u(2,j), bl(2), br(2), -1)
              endif

           if ( (ie+1)==npx ) then
              bl(npx-2) = al(npx-2) - u(npx-2,j)
              xt = s15*u(npx-1,j) + s11*u(npx-2,j) + s14*dm(npx-2)
              br(npx-2) = xt - u(npx-2,j)
              bl(npx-1) = xt - u(npx-1,j)
              if( j==1 .or. j==npy ) then
                 bl(npx-1) = 0.   ! in
                 br(npx-1) = 0.   ! edge
                 bl(npx  ) = 0.   ! edge
                 br(npx  ) = 0.   ! out
              else
                 br(npx) = s11*dq(npx) - s14*dm(npx+1)
#ifdef ONE_SIDE
                 xt = t14*u(npx-1,j) + t12*u(npx-2,j) + t15*u(npx-3,j)
                 br(npx-1) = 2.*xt - u(npx-1,j)
                 xt = t14*u(npx,j) + t12*u(npx+1,j) + t15*u(npx+2,j)
                 bl(npx  ) = 2.*xt - u(npx  ,j)
#else
                 x0L = 0.5*( (2.*dx(npx-1,j)+dx(npx-2,j))*(u(npx-1,j))  &
                      - dx(npx-1,j)*(u(npx-2,j)))/(dx(npx-1,j)+dx(npx-2,j))
                 x0R = 0.5*( (2.*dx(npx,j)+dx(npx+1,j))*(u(npx,j))  &
                      - dx(npx,j)*(u(npx+1,j)))/(dx(npx,j)+dx(npx+1,j))
                 xt = x0L + x0R
                 br(npx-1) = xt - u(npx-1,j)
                 bl(npx  ) = xt - u(npx  ,j)
!                 br(npx-1) = xt + 0.5*(v(npx,j-1)+v(npx,j))*cosa(npx,j) - u(npx-1,j)
!                 bl(npx  ) = xt - 0.5*(v(npx,j-1)+v(npx,j))*cosa(npx,j) - u(npx  ,j)
#endif
                 endif
                 call pert_ppm(1, u(npx-2,j), bl(npx-2), br(npx-2), -1)
              endif
          else
          
              do i=is-1,ie+2
                 al(i) = 0.5*(u(i-1,j)+u(i,j)) + r3*(dm(i-1) - dm(i))
              enddo

              do i=is-1,ie+1
                 pmp = -2.*dq(i)
                 lac = pmp + 1.5*dq(i+1)
                 bl(i) = min(max(0., pmp, lac), max(al(i  )-u(i,j), min(0.,pmp, lac)))
                 pmp = 2.*dq(i-1)
                 lac = pmp - 1.5*dq(i-2)
                 br(i) = min(max(0., pmp, lac), max(al(i+1)-u(i,j), min(0.,pmp, lac)))
              enddo
          endif
       
          do i=is,ie+1
             if( c(i,j)>0. ) then
                cfl = c(i,j)*rdx(i-1,j)
                flux(i,j) = u(i-1,j) + (1.-cfl)*(br(i-1)-cfl*(bl(i-1)+br(i-1)))
             else
                cfl = c(i,j)*rdx(i,j)
                flux(i,j) = u(i,  j) + (1.+cfl)*(bl(i  )+cfl*(bl(i  )+br(i  )))
             endif
          enddo
     enddo

 end select

 end subroutine xtp_u


 subroutine ytp_v(c, u, v, flux, jord)
 integer, intent(IN):: jord
 real, INTENT(IN)  ::   u(isd:ied,jsd:jed+1)
 real, INTENT(IN)  ::   v(isd:ied+1,jsd:jed)
 real, INTENT(IN) ::    c(is:ie+1,js:je+1)   !  Courant   N (like FLUX)
 real, INTENT(OUT):: flux(is:ie+1,js:je+1)
! Local:
 logical extm(is:ie+1,js-2:je+2)
 real dm(is:ie+1,js-2:je+2)
 real al(is:ie+1,js-1:je+2)
 real bl(is:ie+1,js-1:je+1)
 real br(is:ie+1,js-1:je+1)
 real dq(is:ie+1,js-3:je+2)
 real xt, dl, dr, pmp, lac, cfl
 real pmp_1, lac_1, pmp_2, lac_2
 real x0, x1, x0R, x0L
 integer i, j

 select case ( jord )
 case (1)

      do j=js,je+1
         do i=is,ie+1
            if( c(i,j)>0. ) then
               flux(i,j) = v(i,j-1)
            else
               flux(i,j) = v(i,j)
            endif
         enddo
      enddo

 case (2)

   do j=js-2,je+2
      do i=is,ie+1
              xt = 0.25*(v(i,j+1) - v(i,j-1))
         dm(i,j) = sign(min(abs(xt), max(v(i,j-1), v(i,j), v(i,j+1)) - v(i,j),   &
                            v(i,j) - min(v(i,j-1), v(i,j), v(i,j+1))), xt)
      enddo
   enddo

   if (grid_type < 3) then
   if( js==1 ) then
         do i=is,ie+1
            x0L = 0.5*( (2.*dy(i,0)+dy(i,-1))*(v(i,0))   &
               - dy(i,0)*(v(i,-1)))/(dy(i,0)+dy(i,-1))
            x0R = 0.5*( (2.*dy(i,1)+dy(i,2))*(v(i,1))   &
               - dy(i,1)*(v(i,2)))/(dy(i,1)+dy(i,2))
            x0 = x0L + x0R
            x1 = s15*v(i,1) + s11*v(i,2) - s14*dm(i,2)
!           dm(i,1) = x1 - v(i,1)
            dm(i,1) = 0.5*(x1 - x0)
            dm(i,1) = sign(min(abs(dm(i,1)), max(v(i,1), x0, x1) - v(i,1),   &
                                    v(i,1) - min(v(i,1), x0, x1)), dm(i,1))
            x1 = s15*v(i,0) + s11*v(i,-1) + s14*dm(i,-1)
!           dm(i,0) = v(i,0) - x1
            dm(i,0) = 0.5*(x0 - x1)
            dm(i,0) = sign(min(abs(dm(i,0)), max(v(i,0), x0, x1) - v(i,0),   &
                                    v(i,0) - min(v(i,0), x0, x1)), dm(i,0))
         enddo
      if (     is == 1   ) then
           dm(1,0) = 0.
           dm(1,1) = 0.
      endif
      if ( (ie+1) == npx ) then
           dm(npx,0) = 0.
           dm(npx,1) = 0.
      endif
   endif

   if( (je+1)==npy ) then
         do i=is,ie+1
            x0L= 0.5*((2.*dy(i,npy-1)+dy(i,npy-2))*(v(i,npy-1)) -  &
                 dy(i,npy-1)*(v(i,npy-2)))/(dy(i,npy-1)+dy(i,npy-2))
            x0R= 0.5*((2.*dy(i,npy)+dy(i,npy+1))*(v(i,npy)) -  &
                 dy(i,npy)*(v(i,npy+1)))/(dy(i,npy)+dy(i,npy+1))
            x0 = x0L + x0R
            x1 = s15*v(i,npy-1) + s11*v(i,npy-2) + s14*dm(i,npy-2)
!           dm(i,npy-1) = v(i,npy-1) - x1
            dm(i,npy-1) = 0.5*(x0 - x1)
            dm(i,npy-1) = sign(min(abs(dm(i,npy-1)), max(v(i,npy-1), x0, x1) - v(i,npy-1),  &
                                        v(i,npy-1) - min(v(i,npy-1), x0, x1)), dm(i,npy-1))
            x1 = s15*v(i,npy) + s11*v(i,npy+1) - s14*dm(i,npy+1)
!           dm(i,npy) = x1 - v(i,npy)
            dm(i,npy) = 0.5*(x1 - x0)
            dm(i,npy) = sign(min(abs(dm(i,npy)), max(v(i,npy), x0, x1) - v(i,npy),   &
                                      v(i,npy) - min(v(i,npy), x0, x1)), dm(i,npy))
         enddo
      if (     is == 1   ) then
           dm(1,npy-1) = 0.
           dm(1,npy  ) = 0.
      endif
      if ( (ie+1) == npx ) then
           dm(npx,npy-1) = 0.
           dm(npx,npy  ) = 0.
      endif
   endif
   endif

   do j=js,je+1
      do i=is,ie+1
         if( c(i,j)>0. ) then
            flux(i,j) = v(i,j-1) + (1.-c(i,j)*rdy(i,j-1))*dm(i,j-1)
         else
            flux(i,j) = v(i,j  ) - (1.+c(i,j)*rdy(i,j  ))*dm(i,j)
         endif
      enddo
   enddo

 case (4)

   do j=js-2,je+2
      do i=is,ie+1
              xt = 0.25*(v(i,j+1) - v(i,j-1))
         dm(i,j) = sign(min(abs(xt), max(v(i,j-1), v(i,j), v(i,j+1)) - v(i,j),   &
                            v(i,j) - min(v(i,j-1), v(i,j), v(i,j+1))), xt)
      enddo
   enddo

   do j=js-1,je+2
      do i=is,ie+1
         al(i,j) = 0.5*(v(i,j-1)+v(i,j)) + r3*(dm(i,j-1) - dm(i,j))
      enddo
   enddo

   if (grid_type < 3) then
   if( js==1 ) then
         do i=is,ie+1
            x0L = 0.5*( (2.*dy(i,0)+dy(i,-1))*(v(i,0))   &
               - dy(i,0)*(v(i,-1)))/(dy(i,0)+dy(i,-1))
            x0R = 0.5*( (2.*dy(i,1)+dy(i,2))*(v(i,1))   &
               - dy(i,1)*(v(i,2)))/(dy(i,1)+dy(i,2))
            x0 = x0L + x0R
            x1 = s15*v(i,1) + s11*v(i,2) - s14*dm(i,2)
            dm(i,1) = 0.5*(x1 - x0)
!           dm(i,1) = sign(min(abs(dm(i,1)), max(v(i,1), x0, x1) - v(i,1),   &
!                                   v(i,1) - min(v(i,1), x0, x1)), dm(i,1))
            x1 = s15*v(i,0) + s11*v(i,-1) + s14*dm(i,-1)
            dm(i,0) = 0.5*(x0 - x1)
!           dm(i,0) = sign(min(abs(dm(i,0)), max(v(i,0), x0, x1) - v(i,0),   &
!                                   v(i,0) - min(v(i,0), x0, x1)), dm(i,0))
            al(i,0) = 0.5*(v(i,-1)+v(i,0)) + r3*(dm(i,-1) - dm(i,0))
            al(i,1) = x0
            al(i,2) = 0.5*(v(i,1)+v(i,2)) + r3*(dm(i,1) - dm(i,2))
         enddo

         if (     is == 1   ) then
             dm(1,0) = 0.
             dm(1,1) = 0.
            i = 1
            al(i,0) = 0.5*(v(i,-1)+v(i,0)) + r3*(dm(i,-1) - dm(i,0))
            al(i,2) = 0.5*(v(i, 1)+v(i,2)) + r3*(dm(i, 1) - dm(i,2))
         endif
         if ( (ie+1) == npx ) then
             dm(npx,0) = 0.
             dm(npx,1) = 0.
            i = npx
            al(i,0) = 0.5*(v(i,-1)+v(i,0)) + r3*dm(i,-1)
            al(i,2) = 0.5*(v(i, 1)+v(i,2)) - r3*dm(i,2)
         endif
   endif

   if( (je+1)==npy ) then
         do i=is,ie+1
            x0L= 0.5*((2.*dy(i,npy-1)+dy(i,npy-2))*(v(i,npy-1)) -  &
                 dy(i,npy-1)*(v(i,npy-2)))/(dy(i,npy-1)+dy(i,npy-2))
            x0R= 0.5*((2.*dy(i,npy)+dy(i,npy+1))*(v(i,npy)) -  &
                 dy(i,npy)*(v(i,npy+1)))/(dy(i,npy)+dy(i,npy+1))
            x0 = x0L + x0R
            x1 = s15*v(i,npy-1) + s11*v(i,npy-2) + s14*dm(i,npy-2)
            dm(i,npy-1) = 0.5*(x0 - x1)
!           dm(i,npy-1) = sign(min(abs(dm(i,npy-1)), max(v(i,npy-1), x0, x1) - v(i,npy-1),  &
!                                       v(i,npy-1) - min(v(i,npy-1), x0, x1)), dm(i,npy-1))
            x1 = s15*v(i,npy) + s11*v(i,npy+1) - s14*dm(i,npy+1)
            dm(i,npy) = 0.5*(x1 - x0)
!           dm(i,npy) = sign(min(abs(dm(i,npy)), max(v(i,npy), x0, x1) - v(i,npy),   &
!                                     v(i,npy) - min(v(i,npy), x0, x1)), dm(i,npy))
            al(i,npy-1) = 0.5*(v(i,npy-2)+v(i,npy-1)) + r3*(dm(i,npy-2) - dm(i,npy-1))
            al(i,npy  ) = x0
            al(i,npy+1) = 0.5*(v(i,npy)+v(i,npy+1)) + r3*(dm(i,npy) - dm(i,npy+1))
         enddo
         if (     is == 1   ) then
              dm(1,npy-1) = 0.
              dm(1,npy  ) = 0.
            i = 1
            al(i,npy-1) = 0.5*(v(i,npy-2)+v(i,npy-1)) + r3*dm(i,npy-2)
            al(i,npy+1) = 0.5*(v(i,npy  )+v(i,npy+1)) - r3*dm(i,npy+1)
         endif
         if ( (ie+1) == npx ) then
              dm(npx,npy-1) = 0.
              dm(npx,npy  ) = 0.
            i = npx
            al(i,npy-1) = 0.5*(v(i,npy-2)+v(i,npy-1)) + r3*dm(i,npy-2)
            al(i,npy+1) = 0.5*(v(i,npy  )+v(i,npy+1)) - r3*dm(i,npy+1)
        endif
   endif
   endif


   do j=js,je+1
      do i=is,ie+1
         if(c(i,j)>0.) then
            xt = 2.*dm(i,j-1)
            dl = sign(min(abs(xt), abs(al(i,j-1)-v(i,j-1))), xt)
            dr = sign(min(abs(xt), abs(al(i,j)-v(i,j-1))),   xt)
            cfl = c(i,j)*rdy(i,j-1)
            flux(i,j) = v(i,j-1) + (1.-cfl)*(dr + cfl*(dl-dr))
         else
            xt = 2.*dm(i,j)
            dl = sign(min(abs(xt), abs(al(i,j)-v(i,j))),   xt)
            dr = sign(min(abs(xt), abs(al(i,j+1)-v(i,j))), xt)
            cfl = c(i,j)*rdy(i,j)
            flux(i,j) = v(i,j) - (1.+cfl)*(dl + cfl*(dl-dr))
         endif
      enddo
   enddo

 case (6)

   do j=js-3,je+2
      do i=is,ie+1
         dq(i,j) = v(i,j+1) - v(i,j)
      enddo
   enddo
   do j=js-2,je+2
      do i=is,ie+1
         if ( dq(i,j-1)*dq(i,j) > 0. ) then
              extm(i,j) = .false.
         else
              extm(i,j) = .true.
         endif
      enddo
   enddo
   do j=max(3,js-1),min(npy-3,je+1)
      do i=is,ie+1
!!!!        if ( extm(i,j) .and. (extm(i,j-1) .or. extm(i,j+1)) ) then
         if ( extm(i,j-1) .and. extm(i,j) .and. extm(i,j+1) ) then
              bl(i,j) = 0.
              br(i,j) = 0.
         else
              bl(i,j) = b5*v(i,j-2) + b4*v(i,j-1) + b3*v(i,j) + b2*v(i,j+1) + b1*v(i,j+2)
              br(i,j) = b1*v(i,j-2) + b2*v(i,j-1) + b3*v(i,j) + b4*v(i,j+1) + b5*v(i,j+2)
         endif
      enddo
   enddo

   if (grid_type < 3) then
   if( js==1 ) then
       do i=is,ie+1
          br(i,2) = p1*(v(i,2)+v(i,3)) + p2*(v(i,1)+v(i,4)) - v(i,2)
               xt = c3*v(i,1) + c2*v(i,2) + c1*v(i,3)
          br(i,1) = xt - v(i,1)
          bl(i,2) = xt - v(i,2)

          bl(i,0) = c1*v(i,-2) + c2*v(i,-1) + c3*v(i,0) - v(i,0)

            x0L = 0.5*( (2.*dy(i,0)+dy(i,-1))*(v(i,0))   &
               - dy(i,0)*(v(i,-1)))/(dy(i,0)+dy(i,-1))
            x0R = 0.5*( (2.*dy(i,1)+dy(i,2))*(v(i,1))   &
               - dy(i,1)*(v(i,2)))/(dy(i,1)+dy(i,2))
            xt = x0L + x0R

          bl(i,1) = xt - v(i,1)
          br(i,0) = xt - v(i,0)
       enddo
       if ( is==1 ) then
            bl(1,0) = 0.  ! out
            br(1,0) = 0.  ! edge
            bl(1,1) = 0.  ! edge
            br(1,1) = 0.  ! in
       endif
       if ( (ie+1)==npx ) then
            bl(npx,0) = 0.   ! out
            br(npx,0) = 0.   ! edge
            bl(npx,1) = 0.   ! edge
            br(npx,1) = 0.   ! in
       endif
   endif

   if( (je+1)==npy ) then
       do i=is,ie+1
            bl(i,npy-2) = p1*(v(i,npy-3)+v(i,npy-2)) + p2*(v(i,npy-4)+v(i,npy-1)) - v(i,npy-2)
            xt = c1*v(i,npy-3) + c2*v(i,npy-2) + c3*v(i,npy-1)

            br(i,npy-2) = xt - v(i,npy-2)
            bl(i,npy-1) = xt - v(i,npy-1)
            br(i,npy) = c3*v(i,npy)+ c2*v(i,npy+1) + c1*v(i,npy+2) - v(i,npy)
            !xt = 0.5*((2.*dy(i,npy-1)+dy(i,npy-2))*(v(i,npy-1)+v(i,npy)) -  &
            !     dy(i,npy-1)*(v(i,npy-2)+v(i,npy+1)))/(dy(i,npy-1)+dy(i,npy-2))
            x0L= 0.5*((2.*dy(i,npy-1)+dy(i,npy-2))*(v(i,npy-1)) -  &
                 dy(i,npy-1)*(v(i,npy-2)))/(dy(i,npy-1)+dy(i,npy-2))
            x0R= 0.5*((2.*dy(i,npy)+dy(i,npy+1))*(v(i,npy)) -  &
                 dy(i,npy)*(v(i,npy+1)))/(dy(i,npy)+dy(i,npy+1))
            xt = x0L + x0R
            br(i,npy-1) = xt - v(i,npy-1)
            bl(i,npy  ) = xt - v(i,npy)
       enddo
       if ( is==1 ) then
            bl(1,npy-1) = 0.  ! in
            br(1,npy-1) = 0.  ! edge
            bl(1,npy  ) = 0.  ! edge
            br(1,npy  ) = 0.  ! out
       endif
       if ( (ie+1)==npx ) then
            bl(npx,npy-1) = 0.  ! in
            br(npx,npy-1) = 0.  ! edge
            bl(npx,npy  ) = 0.  ! edge
            br(npx,npy  ) = 0.  ! out
       endif
   endif
   endif

   do j=js,je+1
      do i=is,ie+1
         if( c(i,j)>0. ) then
                  cfl = c(i,j)*rdy(i,j-1)
            flux(i,j) = v(i,j-1) + (1.-cfl)*(br(i,j-1)-cfl*(bl(i,j-1)+br(i,j-1)))
         else
                  cfl = c(i,j)*rdy(i,j)
            flux(i,j) = v(i,j  ) + (1.+cfl)*(bl(i,j  )+cfl*(bl(i,j  )+br(i,j  )))
         endif
      enddo
   enddo

 case default
! jord= 8, 9, 10

   do j=js-2,je+2
      do i=is,ie+1
         xt = 0.25*(v(i,j+1) - v(i,j-1))
         dm(i,j) = sign(min(abs(xt), max(v(i,j-1), v(i,j), v(i,j+1)) - v(i,j),   &
                            v(i,j) - min(v(i,j-1), v(i,j), v(i,j+1))), xt)
      enddo
   enddo

   do j=js-3,je+2
      do i=is,ie+1
         dq(i,j) = v(i,j+1) - v(i,j)
      enddo
   enddo

   if (grid_type < 3) then
      do j=max(3,js-1),min(npy-2,je+2)
         do i=is,ie+1
            al(i,j) = 0.5*(v(i,j-1)+v(i,j)) + r3*(dm(i,j-1)-dm(i,j))
         enddo
      enddo
      
      if ( jord==8 ) then
        do j=max(3,js-1),min(npy-3,je+1)
           do i=is,ie+1
              xt =  2.*dm(i,j)
              bl(i,j) = -sign(min(abs(xt), abs(al(i,j)-v(i,j))),   xt)
              br(i,j) =  sign(min(abs(xt), abs(al(i,j+1)-v(i,j))), xt)
           enddo
        enddo
      elseif ( jord==9 ) then
        do j=max(3,js-1),min(npy-3,je+1)
           do i=is,ie+1
              pmp_1 = -2.*dq(i,j) 
              lac_1 = pmp_1 + 1.5*dq(i,j+1)
            bl(i,j) = min(max(0., pmp_1, lac_1), max(al(i,j)-v(i,j), min(0., pmp_1, lac_1)))
              pmp_2 = 2.*dq(i,j-1)
              lac_2 = pmp_2 - 1.5*dq(i,j-2)
            br(i,j) = min(max(0., pmp_2, lac_2), max(al(i,j+1)-v(i,j), min(0., pmp_2, lac_2)))
         enddo
      enddo
    elseif ( jord==10 ) then
      do j=max(3,js-1),min(npy-3,je+1)
         do i=is,ie+1
            bl(i,j) = al(i,j  ) - v(i,j)
            br(i,j) = al(i,j+1) - v(i,j)
!           if ( abs(dm(i,j-1))+abs(dm(i,j))+abs(dm(i,j+1)) < near_zero ) then
            if ( abs(dm(i,j)) < near_zero ) then
              if ( abs(dm(i,j-1))+abs(dm(i,j+1)) < near_zero ) then
                 bl(i,j) = 0.
                 br(i,j) = 0.
              endif
            elseif( abs(3.*(bl(i,j)+br(i,j))) > abs(bl(i,j)-br(i,j)) ) then
                  pmp_1 = -2.*dq(i,j) 
                  lac_1 = pmp_1 + 1.5*dq(i,j+1)
                bl(i,j) = min(max(0., pmp_1, lac_1), max(bl(i,j), min(0., pmp_1, lac_1)))
                  pmp_2 = 2.*dq(i,j-1)
                  lac_2 = pmp_2 - 1.5*dq(i,j-2)
                br(i,j) = min(max(0., pmp_2, lac_2), max(br(i,j), min(0., pmp_2, lac_2)))
            endif
           enddo
        enddo
      else
! Unlimited:
        do j=max(3,js-1),min(npy-3,je+1)
           do i=is,ie+1
              bl(i,j) = al(i,j  ) - v(i,j)
              br(i,j) = al(i,j+1) - v(i,j)
           enddo
        enddo
      endif
      
!--------------
! fix the edges
!--------------
      if( js==1 ) then
         do i=is,ie+1
            br(i,2) = al(i,3) - v(i,2)
            xt = s15*v(i,1) + s11*v(i,2) - s14*dm(i,2)
            br(i,1) = xt - v(i,1)
            bl(i,2) = xt - v(i,2)

            bl(i,0) = s14*dm(i,-1) - s11*dq(i,-1)

#ifdef ONE_SIDE
            xt =  t14*v(i,1) +  t12*v(i,2) + t15*v(i,3)
            bl(i,1) = 2.*xt - v(i,1)
            xt =  t14*v(i,0) +  t12*v(i,-1) + t15*v(i,-2)
            br(i,0) = 2.*xt - v(i,0)
#else
            x0L = 0.5*( (2.*dy(i,0)+dy(i,-1))*(v(i,0))   &
               - dy(i,0)*(v(i,-1)))/(dy(i,0)+dy(i,-1))
            x0R = 0.5*( (2.*dy(i,1)+dy(i,2))*(v(i,1))   &
               - dy(i,1)*(v(i,2)))/(dy(i,1)+dy(i,2))
            xt = x0L + x0R

             bl(i,1) = xt - v(i,1)
             br(i,0) = xt - v(i,0)
!             br(i,0) = xt - 0.5*(u(i-1,1)+u(i,1))*cosa(i,1) - v(i,0)
!             bl(i,1) = xt + 0.5*(u(i-1,1)+u(i,1))*cosa(i,1) - v(i,1)
#endif
         enddo
         if ( is==1 ) then
               bl(1,0) = 0.   ! out
               br(1,0) = 0.   ! edge
               bl(1,1) = 0.   ! edge
               br(1,1) = 0.   ! in
         endif
         if ( (ie+1)==npx ) then
               bl(npx,0) = 0.   ! out
               br(npx,0) = 0.   ! edge
               bl(npx,1) = 0.   ! edge
               br(npx,1) = 0.   ! in
         endif
         j=2
         call pert_ppm(ie-is+2, v(is,j), bl(is,j), br(is,j), -1)
      endif

      if( (je+1)==npy ) then
         do i=is,ie+1
            bl(i,npy-2) = al(i,npy-2) - v(i,npy-2)
            xt = s15*v(i,npy-1) + s11*v(i,npy-2) + s14*dm(i,npy-2)
            br(i,npy-2) = xt - v(i,npy-2)
            bl(i,npy-1) = xt - v(i,npy-1)
            br(i,npy) = s11*dq(i,npy) - s14*dm(i,npy+1)
#ifdef ONE_SIDE
            xt = t14*v(i,npy-1) + t12*v(i,npy-2) + t15*v(i,npy-3)
            br(i,npy-1) = 2.*xt - v(i,npy-1)
            xt = t14*v(i,npy) + t12*v(i,npy+1) + t15*v(i,npy+2)
            bl(i,npy  ) = 2.*xt - v(i,npy)
#else
            x0L= 0.5*((2.*dy(i,npy-1)+dy(i,npy-2))*(v(i,npy-1)) -  &
                 dy(i,npy-1)*(v(i,npy-2)))/(dy(i,npy-1)+dy(i,npy-2))
            x0R= 0.5*((2.*dy(i,npy)+dy(i,npy+1))*(v(i,npy)) -  &
                 dy(i,npy)*(v(i,npy+1)))/(dy(i,npy)+dy(i,npy+1))
            xt = x0L + x0R

            br(i,npy-1) = xt - v(i,npy-1)
            bl(i,npy  ) = xt - v(i,npy)
!            br(i,npy-1) = xt + 0.5*(u(i-1,npy)+u(i,npy))*cosa(i,npy) - v(i,npy-1)
!            bl(i,npy  ) = xt - 0.5*(u(i-1,npy)+u(i,npy))*cosa(i,npy) - v(i,npy)
#endif
         enddo
         if ( is==1 ) then
               bl(1,npy-1) = 0.   ! in
               br(1,npy-1) = 0.   ! edge
               bl(1,npy  ) = 0.   ! edge
               br(1,npy  ) = 0.   ! out
         endif
         if ( (ie+1)==npx ) then
               bl(npx,npy-1) = 0.   ! in
               br(npx,npy-1) = 0.   ! edge
               bl(npx,npy  ) = 0.   ! edge
               br(npx,npy  ) = 0.   ! out
         endif
         j=npy-2
         call pert_ppm(ie-is+2, v(is,j), bl(is,j), br(is,j), -1)
      endif

   else

      do j=js-1,je+2
         do i=is,ie+1
            al(i,j) = 0.5*(v(i,j-1)+v(i,j)) + r3*(dm(i,j-1)-dm(i,j))
         enddo
      enddo
      
      do j=js-1,je+1
         do i=is,ie+1
            pmp = 2.*dq(i,j-1)
            lac = pmp - 1.5*dq(i,j-2)
            br(i,j) = min(max(0.,pmp,lac), max(al(i,j+1)-v(i,j), min(0.,pmp,lac)))
            pmp = -2.*dq(i,j) 
            lac = pmp + 1.5*dq(i,j+1)
            bl(i,j) = min(max(0.,pmp,lac), max(al(i,j)-v(i,j), min(0.,pmp,lac)))
         enddo
      enddo
      
   endif

   do j=js,je+1
      do i=is,ie+1
         if(c(i,j)>0.) then
            cfl = c(i,j)*rdy(i,j-1)
            flux(i,j) = v(i,j-1) + (1.-cfl)*(br(i,j-1)-cfl*(bl(i,j-1)+br(i,j-1)))
         else
            cfl = c(i,j)*rdy(i,j)
            flux(i,j) = v(i,j  ) + (1.+cfl)*(bl(i,j  )+cfl*(bl(i,j  )+br(i,j  )))
         endif
      enddo
   enddo

 end select

end subroutine ytp_v



 subroutine d2a2c_vect(u, v, ua, va, uc, vc, ut, vt, dord4)
  logical, intent(in):: dord4
  real, intent(in) ::  u(isd:ied,jsd:jed+1)
  real, intent(in) ::  v(isd:ied+1,jsd:jed)
  real, intent(out), dimension(isd:ied+1,jsd:jed  ):: uc
  real, intent(out), dimension(isd:ied  ,jsd:jed+1):: vc
  real, intent(out), dimension(isd:ied  ,jsd:jed  ):: ua, va, ut, vt
! Local 
  real, dimension(isd:ied,jsd:jed):: utmp, vtmp
  integer npt, i, j, ifirst, ilast, id

  if ( dord4 ) then
       id = 1
  else
       id = 0
  endif


  if (grid_type < 3) then
     npt = 4
  else
     npt = -2
  endif

! Initialize the non-existing corner regions
  utmp=1.e25
  vtmp=1.e25

!----------
! Interior:
!----------
  do j=max(npt,js-1),min(npy-npt,je+1)
     do i=max(npt,isd),min(npx-npt,ied)
        utmp(i,j) = a2*(u(i,j-1)+u(i,j+2)) + a1*(u(i,j)+u(i,j+1))
     enddo
  enddo
  do j=max(npt,jsd),min(npy-npt,jed)
     do i=max(npt,is-1),min(npx-npt,ie+1)
        vtmp(i,j) = a2*(v(i-1,j)+v(i+2,j)) + a1*(v(i,j)+v(i+1,j))
     enddo
  enddo

!----------
! edges:
!----------
  if (grid_type < 3) then

  if ( js==1 .or. jsd<npt) then
      do j=jsd,npt-1
         do i=isd,ied
            utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
            vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
         enddo
      enddo
  endif

  if ( (je+1)==npy .or. jed>=(npy-npt)) then
      do j=npy-npt+1,jed
         do i=isd,ied
            utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
            vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
         enddo
      enddo
  endif

  if ( is==1 .or. isd<npt ) then
      do j=max(npt,jsd),min(npy-npt,jed)
         do i=isd,npt-1
            utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
            vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
         enddo
      enddo
  endif

  if ( (ie+1)==npx .or. ied>=(npx-npt)) then
      do j=max(npt,jsd),min(npy-npt,jed)
         do i=npx-npt+1,ied
            utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
            vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
         enddo
      enddo
  endif

  endif

  do j=js-1-id,je+1+id
     do i=is-1-id,ie+1+id
        ua(i,j) = (utmp(i,j)-vtmp(i,j)*cosa_s(i,j)) * rsin2(i,j)
        va(i,j) = (vtmp(i,j)-utmp(i,j)*cosa_s(i,j)) * rsin2(i,j)
     enddo
  enddo

! A -> C
!--------------
! Fix the edges
!--------------
! Xdir:
     if( sw_corner ) then
         do i=-2,0
            utmp(i,0) = -vtmp(0,1-i)
         enddo
     endif
     if( se_corner ) then
         do i=0,2
            utmp(npx+i,0) = vtmp(npx,i+1)
         enddo
     endif
     if( ne_corner ) then
         do i=0,2
            utmp(npx+i,npy) = -vtmp(npx,je-i)
         enddo
     endif
     if( nw_corner ) then
         do i=-2,0
            utmp(i,npy) = vtmp(0,je+i)
         enddo
     endif

  if (grid_type < 3) then
     ifirst = max(3,    is-1)
     ilast  = min(npx-2,ie+2)
  else
     ifirst = is-1
     ilast  = ie+2
  endif
!---------------------------------------------
! 4th order interpolation for interior points:
!---------------------------------------------
     do j=js-1,je+1
        do i=ifirst,ilast
           uc(i,j) = a1*(utmp(i-1,j)+utmp(i,j))+a2*(utmp(i-2,j)+utmp(i+1,j))
           ut(i,j) = (uc(i,j) - v(i,j)*cosa_u(i,j))*rsin_u(i,j)
        enddo
     enddo

     if (grid_type < 3) then
#ifndef TEST_NEW
! Xdir:
     if( sw_corner ) then
         ua(-1,0) = -va(0,2)
         ua( 0,0) = -va(0,1) 
     endif
     if( se_corner ) then
         ua(npx,  0) = va(npx,1)
         ua(npx+1,0) = va(npx,2) 
     endif
     if( ne_corner ) then
         ua(npx,  npy) = -va(npx,npy-1)
         ua(npx+1,npy) = -va(npx,npy-2) 
     endif
     if( nw_corner ) then
         ua(-1,npy) = va(0,npy-2)
         ua( 0,npy) = va(0,npy-1) 
     endif
#endif

     if( is==1 ) then
        do j=js-1,je+1
           uc(0,j) = c1*utmp(-2,j) + c2*utmp(-1,j) + c3*utmp(0,j) 
#ifndef TEST_NEW
           !ut(1,j) = 0.25*(-ua(-1,j) + 3.*(ua(0,j)+ua(1,j)) - ua(2,j))
!           uc(1,j) = 0.25*(-utmp(-1,j) + 3.*(utmp(0,j)+utmp(1,j)) - utmp(2,j))
           ut(1,j) = edge_interpolate4(ua(-1:2,j), dxa(-1:2,j))
           !Using the DOWNSTREAM sin value actually gives LOWER errors in stretched case 2, although only a bit (Order 10-6 improvement)
           if (ut(1,j) < 0.) then
              uc(1,j) = ut(1,j)*sin_sg(0,j,3)
           else
              uc(1,j) = ut(1,j)*sin_sg(1,j,1)
           end if
#else
! 3-pt extrapolation: grid symmetry assumed --------------------------------
           uc(1,j) = ( t14*(utmp( 0,j)+utmp(1,j))    &
                     + t12*(utmp(-1,j)+utmp(2,j))    &
                     + t15*(utmp(-2,j)+utmp(3,j)) )*rsin_u(1,j)
           ut(1,j) =  uc(1,j) * rsin_u(1,j)
! 3-pt extrapolation: grid symmetry assumed --------------------------------
#endif
           uc(2,j) = c1*utmp(3,j) + c2*utmp(2,j) + c3*utmp(1,j)
           ut(0,j) = (uc(0,j) - v(0,j)*cosa_u(0,j))*rsin_u(0,j)
           ut(2,j) = (uc(2,j) - v(2,j)*cosa_u(2,j))*rsin_u(2,j)
        enddo
     endif

     if( (ie+1)==npx ) then
        do j=js-1,je+1
           uc(npx-1,j) = c1*utmp(npx-3,j)+c2*utmp(npx-2,j)+c3*utmp(npx-1,j) 
#ifndef TEST_NEW
        i=npx
                   ut(i,j) = 0.25*(-ua(i-2,j) + 3.*(ua(i-1,j)+ua(i,j)) - ua(i+1,j))
        ut(i,j) = edge_interpolate4(ua(i-2:i+1,j), dxa(i-2:i+1,j))
        if (ut(i,j) < 0.) then
           uc(i,j) = ut(i,j)*sin_sg(i-1,j,3)
        else
           uc(i,j) = ut(i,j)*sin_sg(i,j,1)
        end if
#else
! 3-pt extrapolation --------------------------------------------------------
           uc(npx,j) = (t14*(utmp(npx-1,j)+utmp(npx,j))+      &
                        t12*(utmp(npx-2,j)+utmp(npx+1,j))     &
                      + t15*(utmp(npx-3,j)+utmp(npx+2,j)))*rsin_u(npx,j)
           ut(npx,  j) =  uc(npx,j) * rsin_u(npx,j)
! 3-pt extrapolation --------------------------------------------------------
#endif
           uc(npx+1,j) = c3*utmp(npx,j)+c2*utmp(npx+1,j)+c1*utmp(npx+2,j) 
           ut(npx-1,j) = (uc(npx-1,j)-v(npx-1,j)*cosa_u(npx-1,j))*rsin_u(npx-1,j)
           ut(npx+1,j) = (uc(npx+1,j)-v(npx+1,j)*cosa_u(npx+1,j))*rsin_u(npx+1,j)
        enddo
     endif

     endif

!------
! Ydir:
!------
     if( sw_corner ) then
         do j=-2,0
            vtmp(0,j) = -utmp(1-j,0)
         enddo
     endif
     if( nw_corner ) then
         do j=0,2
            vtmp(0,npy+j) = utmp(j+1,npy)
         enddo
     endif
     if( se_corner ) then
         do j=-2,0
            vtmp(npx,j) = utmp(ie+j,0)
         enddo
     endif
     if( ne_corner ) then
         do j=0,2
            vtmp(npx,npy+j) = -utmp(ie-j,npy)
         enddo
     endif
#ifndef TEST_NEW
     if( sw_corner ) then
         va(0,-1) = -ua(2,0)
         va(0, 0) = -ua(1,0)
     endif
     if( se_corner ) then
         va(npx, 0) = ua(npx-1,0)
         va(npx,-1) = ua(npx-2,0)
     endif
     if( ne_corner ) then
         va(npx,npy  ) = -ua(npx-1,npy)
         va(npx,npy+1) = -ua(npx-2,npy)
     endif
     if( nw_corner ) then
         va(0,npy)   = ua(1,npy)
         va(0,npy+1) = ua(2,npy)
     endif
#endif

     if (grid_type < 3) then

     do j=js-1,je+2
      if ( j==1 ) then
        do i=is-1,ie+1
#ifndef TEST_NEW
           !vt(i,j) = 0.25*(-va(i,j-2) + 3.*(va(i,j-1)+va(i,j)) - va(i,j+1))
!           vc(i,j) = 0.25*(-vtmp(i,j-2) + 3.*(vtmp(i,j-1)+vtmp(i,j)) - vtmp(i,j+1))
           vt(i,j) = edge_interpolate4(va(i,-1:2), dya(i,-1:2))
           if (vt(i,j) < 0.) then
              vc(i,j) = vt(i,j)*sin_sg(i,j-1,4)
           else
              vc(i,j) = vt(i,j)*sin_sg(i,j,2)
           end if
#else
! 3-pt extrapolation -----------------------------------------
           vc(i,1) = (t14*(vtmp(i, 0)+vtmp(i,1))    &
                    + t12*(vtmp(i,-1)+vtmp(i,2))    &
                    + t15*(vtmp(i,-2)+vtmp(i,3)))*rsin_v(i,1)
           vt(i,1) = vc(i,1) * rsin_v(i,1)
! 3-pt extrapolation -----------------------------------------
#endif
        enddo
      elseif ( j==0 .or. j==(npy-1) ) then
        do i=is-1,ie+1
           vc(i,j) = c1*vtmp(i,j-2) + c2*vtmp(i,j-1) + c3*vtmp(i,j)
           vt(i,j) = (vc(i,j) - u(i,j)*cosa_v(i,j))*rsin_v(i,j)
        enddo
      elseif ( j==2 .or. j==(npy+1) ) then
        do i=is-1,ie+1
           vc(i,j) = c1*vtmp(i,j+1) + c2*vtmp(i,j) + c3*vtmp(i,j-1)
           vt(i,j) = (vc(i,j) - u(i,j)*cosa_v(i,j))*rsin_v(i,j)
        enddo
      elseif ( j==npy ) then
        do i=is-1,ie+1
#ifndef TEST_NEW
           vt(i,j) = 0.25*(-va(i,j-2) + 3.*(va(i,j-1)+va(i,j)) - va(i,j+1))
!           vc(i,j) = 0.25*(-vtmp(i,j-2) + 3.*(vtmp(i,j-1)+vtmp(i,j)) - vtmp(i,j+1))
           vt(i,j) = edge_interpolate4(va(i,j-2:j+1), dya(i,j-2:j+1))
           if (vt(i,j) < 0.) then
              vc(i,j) = vt(i,j)*sin_sg(i,j-1,4)
           else
              vc(i,j) = vt(i,j)*sin_sg(i,j,2)
           end if
#else
! 3-pt extrapolation --------------------------------------------------------
           vc(i,npy) = (t14*(vtmp(i,npy-1)+vtmp(i,npy))    &
                      + t12*(vtmp(i,npy-2)+vtmp(i,npy+1))  &
                      + t15*(vtmp(i,npy-3)+vtmp(i,npy+2)))*rsin_v(i,npy)
           vt(i,npy) = vc(i,npy) * rsin_v(i,npy)
! 3-pt extrapolation -----------------------------------------
#endif
        enddo
      else
! 4th order interpolation for interior points:
        do i=is-1,ie+1
           vc(i,j) = a2*(vtmp(i,j-2)+vtmp(i,j+1))+a1*(vtmp(i,j-1)+vtmp(i,j))
           vt(i,j) = (vc(i,j) - u(i,j)*cosa_v(i,j))*rsin_v(i,j)
        enddo
      endif
     enddo
    else
! 4th order interpolation:
       do j=js-1,je+2
          do i=is-1,ie+1
             vc(i,j) = a2*(vtmp(i,j-2)+vtmp(i,j+1))+a1*(vtmp(i,j-1)+vtmp(i,j))
             vt(i,j) = vc(i,j)
          enddo
       enddo
    endif

 end subroutine d2a2c_vect
 

 subroutine d2a2c_vect_v2( u, v, ua, va, uc, vc, ut, vt )
  real, intent(in) ::  u(isd:ied,jsd:jed+1)
  real, intent(in) ::  v(isd:ied+1,jsd:jed)
  real, intent(out), dimension(isd:ied+1,jsd:jed  ):: uc
  real, intent(out), dimension(isd:ied  ,jsd:jed+1):: vc
  real, intent(out), dimension(isd:ied  ,jsd:jed  ):: ua, va, ut, vt
! Local 
    real, dimension(is-2:ie+2,js-2:je+2):: wk
    real :: utmp, vtmp
    integer i, j

! needs only ut[is-1:ie+2,js-1:je+1], vt[is-1:ie+1,js-1:je+2]

     do j=js-2,je+2
        do i=is-2,ie+3
           uc(i,j) = v(i,j)*dy(i,j)
        enddo
     enddo
     do j=js-2,je+3
        do i=is-2,ie+2
           vc(i,j) = u(i,j)*dx(i,j)
        enddo
     enddo

! D --> A
! Co-variant to Co-variant "vorticity-conserving" interpolation
     do j=js-2,je+2
        do i=is-2,ie+2
           utmp = 0.5*(vc(i,j) + vc(i,j+1)) * rdxa(i,j)
           vtmp = 0.5*(uc(i,j) + uc(i+1,j)) * rdya(i,j)
           ua(i,j) = (utmp-vtmp*cosa_s(i,j))*rsin2(i,j)
           va(i,j) = (vtmp-utmp*cosa_s(i,j))*rsin2(i,j)
        enddo
     enddo

! Xdir:
     if( sw_corner ) then
         ua(-1,0) = -va(0,2)
         ua( 0,0) = -va(0,1) 
     endif
     if( se_corner ) then
         ua(npx,  0) = va(npx,1)
         ua(npx+1,0) = va(npx,2) 
     endif
     if( ne_corner ) then
         ua(npx,  npy) = -va(npx,npy-1)
         ua(npx+1,npy) = -va(npx,npy-2) 
     endif
     if( nw_corner ) then
         ua(-1,npy) = va(0,npy-2)
         ua( 0,npy) = va(0,npy-1) 
     endif

! A -> C
!--------------------------------------------
! Divergence conserving interp to cell walls
!--------------------------------------------
     do j=js-1,je+1
        do i=is-2,ie+2
           wk(i,j) = ua(i,j)*dya(i,j)*sina_s(i,j)
        enddo
     enddo
     do j=js-1,je+1
        do i=is-1,ie+2
           ut(i,j) = 0.5*(wk(i-1,j)+wk(i,j)) / (dy(i,j)*sina_u(i,j))
           uc(i,j) = ut(i,j) + 0.5*(va(i-1,j)*cosa_s(i-1,j)+va(i,j)*cosa_s(i,j))
        enddo
     enddo

     if (grid_type < 3) then
     if ( is==1 ) then
        i=1
        do j=js-1,je+1
           ut(i,j) = 0.25*(-ua(-1,j) + 3.*(ua(0,j)+ua(1,j)) - ua(2,j))
           uc(i,j) = ut(i,j)*sina_u(i,j)
        enddo
     endif

     if ( (ie+1)==npx ) then
        i=npx
        do j=js-1,je+1
           ut(i,j) = 0.25*(-ua(i-2,j) + 3.*(ua(i-1,j)+ua(i,j)) - ua(i+1,j))
           uc(i,j) = ut(i,j)*sina_u(i,j)
        enddo
     endif
     endif

! Ydir:
     if( sw_corner ) then
         va(0,-1) = -ua(2,0)
         va(0, 0) = -ua(1,0)
     endif
     if( se_corner ) then
         va(npx, 0) = ua(npx-1,0)
         va(npx,-1) = ua(npx-2,0)
     endif
     if( ne_corner ) then
         va(npx,npy  ) = -ua(npx-1,npy)
         va(npx,npy+1) = -ua(npx-2,npy)
     endif
     if( nw_corner ) then
         va(0,npy)   = ua(1,npy)
         va(0,npy+1) = ua(2,npy)
     endif

     do j=js-2,je+2
        do i=is-1,ie+1
           wk(i,j) = va(i,j)*dxa(i,j)*sina_s(i,j)
        enddo
     enddo

     if (grid_type < 3) then
     do j=js-1,je+2
        if ( j==1 .or. j==npy ) then
          do i=is-1,ie+1
             vt(i,j) = 0.25*(-va(i,j-2) + 3.*(va(i,j-1)+va(i,j)) - va(i,j+1))
             vc(i,j) = vt(i,j)*sina_v(i,j)
          enddo
        else
          do i=is-1,ie+1
             vt(i,j) = 0.5*(wk(i,j-1)+wk(i,j)) / (dx(i,j)*sina_v(i,j))
             vc(i,j) = vt(i,j) + 0.5*(ua(i,j-1)*cosa_s(i,j-1)+ua(i,j)*cosa_s(i,j))
          enddo
        endif
     enddo
     else
        do j=js-1,je+2
           do i=is-1,ie+1
              vt(i,j) = 0.5*(wk(i,j-1)+wk(i,j)) / (dx(i,j)*sina_v(i,j))
              vc(i,j) = vt(i,j) + 0.5*(ua(i,j-1)*cosa_s(i,j-1)+ua(i,j)*cosa_s(i,j))
           enddo
        enddo
     endif

 end subroutine d2a2c_vect_v2

  real function edge_interpolate4(ua, dxa)

   real, intent(in) :: ua(4)
   real, intent(in) :: dxa(4)

   real u0L, u0R

   u0L = 0.5*((2.*dxa(2)+dxa(1))*ua(2) - dxa(2)*ua(1)) / ( dxa(1)+dxa(2) )
   u0R = 0.5*((2.*dxa(3)+dxa(4))*ua(3) - dxa(3)*ua(4)) / ( dxa(3)+dxa(4) )
   edge_interpolate4 = u0L + u0R

   !This is the original edge-interpolation code, which makes
   ! a relatively small increase in the error in unstretched case 2.

!   edge_interpolate4 = 0.25*( 3*(ua(2)+ua(3)) - (ua(1)+ua(4))  )

 end function edge_interpolate4

      
 subroutine d2a2c_vect_v1( u,v, ua,va, uc,vc, ut,vt )
  real, intent(in) ::  u(isd:ied,jsd:jed+1)
  real, intent(in) ::  v(isd:ied+1,jsd:jed)
  real, intent(out), dimension(isd:ied+1,jsd:jed  ):: uc
  real, intent(out), dimension(isd:ied  ,jsd:jed+1):: vc
  real, intent(out), dimension(isd:ied  ,jsd:jed  ):: ua, va, ut, vt
! Local 
  real, dimension(isd:ied,jsd:jed):: v1, v2, v3
  real, dimension(isd:ied,jsd:jed):: utmp, vtmp
    real vw1, vw2, vw3
    real vs1, vs2, vs3
    real up, vp
    integer i, j

! Needs ut[is-1:ie+2,js-1:je+1], vt[is-1:ie+1,js-1:je+2]

     do j=jsd,jed
        do i=isd,ied+1
           uc(i,j) = v(i,j)*dy(i,j)
        enddo
     enddo
     do j=jsd,jed+1
        do i=isd,ied
           vc(i,j) = u(i,j)*dx(i,j)
        enddo
     enddo

! D --> A
     do j=jsd,jed
        do i=isd,ied
           up = 0.5*(vc(i,j) + vc(i,j+1)) * rdxa(i,j)
           vp = 0.5*(uc(i,j) + uc(i+1,j)) * rdya(i,j)
           ua(i,j) = (up-vp*cosa_s(i,j)) * rsin2(i,j)
           va(i,j) = (vp-up*cosa_s(i,j)) * rsin2(i,j)
           v1(i,j) = ua(i,j)*ec1(1,i,j) + va(i,j)*ec2(1,i,j)
           v2(i,j) = ua(i,j)*ec1(2,i,j) + va(i,j)*ec2(2,i,j)
           v3(i,j) = ua(i,j)*ec1(3,i,j) + va(i,j)*ec2(3,i,j)
        enddo
     enddo

! A -> C (across face averaging taking place here):
! Xdir
     call fill3_4corners(v1, v2, v3, 1)
!    call copy_corners(v1, npx, npy, 1)
!    call copy_corners(v2, npx, npy, 1)
!    call copy_corners(v3, npx, npy, 1)

! 4th order interpolation:
     do j=js-1,je+1
        do i=max(3,is-1),min(npx-2,ie+2)
           vw1 = a2*(v1(i-2,j)+v1(i+1,j)) + a1*(v1(i-1,j)+v1(i,j))
           vw2 = a2*(v2(i-2,j)+v2(i+1,j)) + a1*(v2(i-1,j)+v2(i,j))
           vw3 = a2*(v3(i-2,j)+v3(i+1,j)) + a1*(v3(i-1,j)+v3(i,j))
           uc(i,j) = vw1*ew(1,i,j,1) + vw2*ew(2,i,j,1) + vw3*ew(3,i,j,1)
           ut(i,j) = (uc(i,j)-v(i,j)*cosa_u(i,j)) * rsin_u(i,j)
        enddo
     enddo

! Fix the edge:
     if ( is==1 ) then
        do j=js-1,je+1
        i=0
           vw1 = c1*v1(-2,j) + c2*v1(-1,j) + c3*v1(0,j) 
           vw2 = c1*v2(-2,j) + c2*v2(-1,j) + c3*v2(0,j) 
           vw3 = c1*v3(-2,j) + c2*v3(-1,j) + c3*v3(0,j) 
           uc(i,j) = vw1*ew(1,i,j,1) + vw2*ew(2,i,j,1) + vw3*ew(3,i,j,1)
           ut(i,j) = (uc(i,j)-v(i,j)*cosa_u(i,j)) * rsin_u(i,j)
        i=1
           vw1 = 3.*(v1(0,j)+v1(1,j)) - (v1(-1,j)+v1(2,j))
           vw2 = 3.*(v2(0,j)+v2(1,j)) - (v2(-1,j)+v2(2,j))
           vw3 = 3.*(v3(0,j)+v3(1,j)) - (v3(-1,j)+v3(2,j))
           uc(i,j) = 0.25*(vw1*ew(1,i,j,1) + vw2*ew(2,i,j,1) + vw3*ew(3,i,j,1))
           ut(i,j) = uc(i,j)*rsin_u(i,j)
        i=2
           vw1 = c3*v1(1,j) + c2*v1(2,j) + c1*v1(3,j)
           vw2 = c3*v2(1,j) + c2*v2(2,j) + c1*v2(3,j)
           vw3 = c3*v3(1,j) + c2*v3(2,j) + c1*v3(3,j)
           uc(i,j) = vw1*ew(1,i,j,1) + vw2*ew(2,i,j,1) + vw3*ew(3,i,j,1)
           ut(i,j) = (uc(i,j)-v(i,j)*cosa_u(i,j)) * rsin_u(i,j)
        enddo
     endif

     if ( (ie+1)==npx ) then
        do j=js-1,je+1
        i=npx-1
           vw1 = c1*v1(npx-3,j) + c2*v1(npx-2,j) + c3*v1(npx-1,j) 
           vw2 = c1*v2(npx-3,j) + c2*v2(npx-2,j) + c3*v2(npx-1,j) 
           vw3 = c1*v3(npx-3,j) + c2*v3(npx-2,j) + c3*v3(npx-1,j) 
           uc(i,j) = vw1*ew(1,i,j,1) + vw2*ew(2,i,j,1) + vw3*ew(3,i,j,1)
           ut(i,j) = (uc(i,j)-v(i,j)*cosa_u(i,j)) * rsin_u(i,j)
        i=npx
           vw1 = 3.*(v1(i-1,j)+v1(i,j)) - (v1(i-2,j)+v1(i+1,j))
           vw2 = 3.*(v2(i-1,j)+v2(i,j)) - (v2(i-2,j)+v2(i+1,j))
           vw3 = 3.*(v3(i-1,j)+v3(i,j)) - (v3(i-2,j)+v3(i+1,j))
           uc(i,j) = 0.25*(vw1*ew(1,i,j,1) + vw2*ew(2,i,j,1) + vw3*ew(3,i,j,1))
           ut(i,j) = uc(i,j)*rsin_u(i,j)
        i=npx+1
           vw1 = c3*v1(npx,j) + c2*v1(npx+1,j) + c1*v1(npx+2,j) 
           vw2 = c3*v2(npx,j) + c2*v2(npx+1,j) + c1*v2(npx+2,j) 
           vw3 = c3*v3(npx,j) + c2*v3(npx+1,j) + c1*v3(npx+2,j) 
           uc(i,j) = vw1*ew(1,i,j,1) + vw2*ew(2,i,j,1) + vw3*ew(3,i,j,1)
           ut(i,j) = (uc(i,j)-v(i,j)*cosa_u(i,j)) * rsin_u(i,j)
        enddo
     endif

! Ydir:
     call fill3_4corners(v1, v2, v3, 2)
!    call copy_corners(v1, npx, npy, 2)
!    call copy_corners(v2, npx, npy, 2)
!    call copy_corners(v3, npx, npy, 2)

     do j=js-1,je+2
        if( j==0 .or. j==(npy-1) ) then
          do i=is-1,ie+1
             vs1 = c1*v1(i,j-2) + c2*v1(i,j-1) + c3*v1(i,j)
             vs2 = c1*v2(i,j-2) + c2*v2(i,j-1) + c3*v2(i,j)
             vs3 = c1*v3(i,j-2) + c2*v3(i,j-1) + c3*v3(i,j)
             vc(i,j) = vs1*es(1,i,j,2) + vs2*es(2,i,j,2) + vs3*es(3,i,j,2)
             vt(i,j) = (vc(i,j)-u(i,j)*cosa_v(i,j)) * rsin_v(i,j)
          enddo
        elseif ( j==2 .or. j==(npy+1) ) then
          do i=is-1,ie+1
             vs1 = c3*v1(i,j-1) + c2*v1(i,j) + c1*v1(i,j+1)
             vs2 = c3*v2(i,j-1) + c2*v2(i,j) + c1*v2(i,j+1)
             vs3 = c3*v3(i,j-1) + c2*v3(i,j) + c1*v3(i,j+1)
             vc(i,j) = vs1*es(1,i,j,2) + vs2*es(2,i,j,2) + vs3*es(3,i,j,2)
             vt(i,j) = (vc(i,j)-u(i,j)*cosa_v(i,j)) * rsin_v(i,j)
          enddo
        elseif ( j==1 .or. j==npy ) then
          do i=is-1,ie+1
              vs1 = 3.*(v1(i,j-1)+v1(i,j)) - (v1(i,j-2)+v1(i,j+1))
              vs2 = 3.*(v2(i,j-1)+v2(i,j)) - (v2(i,j-2)+v2(i,j+1))
              vs3 = 3.*(v3(i,j-1)+v3(i,j)) - (v3(i,j-2)+v3(i,j+1))
              vc(i,j) = 0.25*(vs1*es(1,i,j,2) + vs2*es(2,i,j,2) + vs3*es(3,i,j,2))
              vt(i,j) = vc(i,j)*rsin_v(i,j)
          enddo
        else
! Interior: 4th order
          do i=is-1,ie+1
             vs1 = a2*(v1(i,j-2)+v1(i,j+1)) + a1*(v1(i,j-1)+v1(i,j))
             vs2 = a2*(v2(i,j-2)+v2(i,j+1)) + a1*(v2(i,j-1)+v2(i,j))
             vs3 = a2*(v3(i,j-2)+v3(i,j+1)) + a1*(v3(i,j-1)+v3(i,j))
             vc(i,j) = vs1*es(1,i,j,2) + vs2*es(2,i,j,2) + vs3*es(3,i,j,2)
             vt(i,j) = (vc(i,j)-u(i,j)*cosa_v(i,j)) * rsin_v(i,j)
          enddo
        endif
     enddo

 end subroutine d2a2c_vect_v1
      

 subroutine fill3_4corners(q1, q2, q3, dir)
! This routine fill the 4 corners of the scalar fileds only as needed by c_core
  integer, intent(in):: dir                ! 1: x-dir; 2: y-dir
  real, intent(inout):: q1(isd:ied,jsd:jed)
  real, intent(inout):: q2(isd:ied,jsd:jed)
  real, intent(inout):: q3(isd:ied,jsd:jed)
  integer i,j

  select case(dir)
  case(1)
      if ( sw_corner ) then
          q1(-1,0) = q1(0,2); q1(0,0) = q1(0,1); q1(0,-1) = q1(-1,1)
          q2(-1,0) = q2(0,2); q2(0,0) = q2(0,1); q2(0,-1) = q2(-1,1)
          q3(-1,0) = q3(0,2); q3(0,0) = q3(0,1); q3(0,-1) = q3(-1,1)
      endif
      if ( se_corner ) then
          q1(npx+1,0) = q1(npx,2); q1(npx,0) = q1(npx,1); q1(npx,-1) = q1(npx+1,1)
          q2(npx+1,0) = q2(npx,2); q2(npx,0) = q2(npx,1); q2(npx,-1) = q2(npx+1,1)
          q3(npx+1,0) = q3(npx,2); q3(npx,0) = q3(npx,1); q3(npx,-1) = q3(npx+1,1)
      endif
      if ( ne_corner ) then
          q1(npx,npy) = q1(npx,npy-1); q1(npx+1,npy) = q1(npx,npy-2); q1(npx,npy+1) = q1(npx+1,npy-1)
          q2(npx,npy) = q2(npx,npy-1); q2(npx+1,npy) = q2(npx,npy-2); q2(npx,npy+1) = q2(npx+1,npy-1)
          q3(npx,npy) = q3(npx,npy-1); q3(npx+1,npy) = q3(npx,npy-2); q3(npx,npy+1) = q3(npx+1,npy-1)
      endif
      if ( nw_corner ) then
          q1(0,npy) = q1(0,npy-1); q1(-1,npy) = q1(0,npy-2); q1(0,npy+1) = q1(-1,npy-1)
          q2(0,npy) = q2(0,npy-1); q2(-1,npy) = q2(0,npy-2); q2(0,npy+1) = q2(-1,npy-1)
          q3(0,npy) = q3(0,npy-1); q3(-1,npy) = q3(0,npy-2); q3(0,npy+1) = q3(-1,npy-1)
      endif

  case(2)
      if ( sw_corner ) then
          q1(0,0) = q1(1,0); q1(0,-1) = q1(2,0); q1(-1,0) = q1(1,-1)
          q2(0,0) = q2(1,0); q2(0,-1) = q2(2,0); q2(-1,0) = q2(1,-1)
          q3(0,0) = q3(1,0); q3(0,-1) = q3(2,0); q3(-1,0) = q3(1,-1)
      endif
      if ( se_corner ) then
          q1(npx,0) = q1(npx-1,0); q1(npx,-1) = q1(npx-2,0); q1(npx+1,0) = q1(npx-1,-1)
          q2(npx,0) = q2(npx-1,0); q2(npx,-1) = q2(npx-2,0); q2(npx+1,0) = q2(npx-1,-1)
          q3(npx,0) = q3(npx-1,0); q3(npx,-1) = q3(npx-2,0); q3(npx+1,0) = q3(npx-1,-1)
      endif
      if ( ne_corner ) then
          q1(npx,npy) = q1(npx-1,npy); q1(npx,npy+1) = q1(npx-2,npy); q1(npx+1,npy) = q1(npx-1,npy+1)
          q2(npx,npy) = q2(npx-1,npy); q2(npx,npy+1) = q2(npx-2,npy); q2(npx+1,npy) = q2(npx-1,npy+1)
          q3(npx,npy) = q3(npx-1,npy); q3(npx,npy+1) = q3(npx-2,npy); q3(npx+1,npy) = q3(npx-1,npy+1)
      endif
      if ( nw_corner ) then
          q1(0,npy) = q1(1,npy); q1(0,npy+1) = q1(2,npy); q1(-1,npy) = q1(1,npy+1)
          q2(0,npy) = q2(1,npy); q2(0,npy+1) = q2(2,npy); q2(-1,npy) = q2(1,npy+1)
          q3(0,npy) = q3(1,npy); q3(0,npy+1) = q3(2,npy); q3(-1,npy) = q3(1,npy+1)
      endif

  end select
 end subroutine fill3_4corners


 subroutine fill2_4corners(q1, q2, dir)
! This routine fill the 4 corners of the scalar fileds only as needed by c_core
  integer, intent(in):: dir                ! 1: x-dir; 2: y-dir
  real, intent(inout):: q1(isd:ied,jsd:jed)
  real, intent(inout):: q2(isd:ied,jsd:jed)

  select case(dir)
  case(1)
      if ( sw_corner ) then
          q1(-1,0) = q1(0,2);    q1(0,0) = q1(0,1)
          q2(-1,0) = q2(0,2);    q2(0,0) = q2(0,1)
      endif
      if ( se_corner ) then
          q1(npx+1,0) = q1(npx,2); q1(npx,0) = q1(npx,1)
          q2(npx+1,0) = q2(npx,2); q2(npx,0) = q2(npx,1)
      endif
      if ( nw_corner ) then
          q1(0,npy) = q1(0,npy-1); q1(-1,npy) = q1(0,npy-2)
          q2(0,npy) = q2(0,npy-1); q2(-1,npy) = q2(0,npy-2)
      endif
      if ( ne_corner ) then
          q1(npx,npy) = q1(npx,npy-1); q1(npx+1,npy) = q1(npx,npy-2)
          q2(npx,npy) = q2(npx,npy-1); q2(npx+1,npy) = q2(npx,npy-2)
      endif

  case(2)
      if ( sw_corner ) then
          q1(0,0) = q1(1,0); q1(0,-1) = q1(2,0)
          q2(0,0) = q2(1,0); q2(0,-1) = q2(2,0)
      endif
      if ( se_corner ) then
          q1(npx,0) = q1(npx-1,0); q1(npx,-1) = q1(npx-2,0)
          q2(npx,0) = q2(npx-1,0); q2(npx,-1) = q2(npx-2,0)
      endif
      if ( nw_corner ) then
          q1(0,npy) = q1(1,npy); q1(0,npy+1) = q1(2,npy)
          q2(0,npy) = q2(1,npy); q2(0,npy+1) = q2(2,npy)
      endif
      if ( ne_corner ) then
          q1(npx,npy) = q1(npx-1,npy); q1(npx,npy+1) = q1(npx-2,npy)
          q2(npx,npy) = q2(npx-1,npy); q2(npx,npy+1) = q2(npx-2,npy)
      endif

  end select

 end subroutine fill2_4corners

 subroutine fill_4corners(q, dir)
! This routine fill the 4 corners of the scalar fileds only as needed by c_core
  integer, intent(in):: dir                ! 1: x-dir; 2: y-dir
  real, intent(inout):: q(isd:ied,jsd:jed)

  select case(dir)
  case(1)
      if ( sw_corner ) then
          q(-1,0) = q(0,2)
          q( 0,0) = q(0,1)
      endif
      if ( se_corner ) then
          q(npx+1,0) = q(npx,2)
          q(npx,  0) = q(npx,1)
      endif
      if ( nw_corner ) then
          q( 0,npy) = q(0,npy-1)
          q(-1,npy) = q(0,npy-2)
      endif
      if ( ne_corner ) then
          q(npx,  npy) = q(npx,npy-1)
          q(npx+1,npy) = q(npx,npy-2)
      endif

  case(2)
      if ( sw_corner ) then
          q(0, 0) = q(1,0)
          q(0,-1) = q(2,0)
      endif
      if ( se_corner ) then
          q(npx, 0) = q(npx-1,0)
          q(npx,-1) = q(npx-2,0)
      endif
      if ( nw_corner ) then
          q(0,npy  ) = q(1,npy)
          q(0,npy+1) = q(2,npy)
      endif
      if ( ne_corner ) then
          q(npx,npy  ) = q(npx-1,npy)
          q(npx,npy+1) = q(npx-2,npy)
      endif

  end select

 end subroutine fill_4corners



 subroutine d2a2c(u, v, um, vm,  ua, va, uc, vc, dord4)
  real, intent(in), dimension(isd:ied,jsd:jed+1):: u, um
  real, intent(in), dimension(isd:ied+1,jsd:jed):: v, vm
  logical, intent(in):: dord4
  real, intent(out), dimension(isd:ied+1,jsd:jed  ):: uc
  real, intent(out), dimension(isd:ied  ,jsd:jed+1):: vc
  real, intent(out), dimension(isd:ied  ,jsd:jed  ):: ua, va
! Local 
  real, dimension(isd:ied,jsd:jed):: utmp, vtmp
  integer npt, i, j, ifirst, ilast, id

  if ( dord4 ) then
       id = 1
  else
       id = 0
  endif


  if (grid_type < 3) then
     npt = 4
  else
     npt = -2
  endif

!----------
! Interior:
!----------
  do j=max(npt,js-1),min(npy-npt,je+1)
     do i=max(npt,isd),min(npx-npt,ied)
        utmp(i,j) = a2*(u(i,j-1)+u(i,j+2)) + a1*(u(i,j)+u(i,j+1))
     enddo
  enddo
  do j=max(npt,jsd),min(npy-npt,jed)
     do i=max(npt,is-1),min(npx-npt,ie+1)
        vtmp(i,j) = a2*(v(i-1,j)+v(i+2,j)) + a1*(v(i,j)+v(i+1,j))
     enddo
  enddo

!----------
! edges:
!----------
  if (grid_type < 3) then

  if ( js==1 .or. jsd<npt) then
      do j=jsd,npt-1
         do i=isd,ied
            utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
            vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
         enddo
      enddo
  endif
  if ( (je+1)==npy .or. jed>=(npy-npt)) then
      do j=npy-npt+1,jed
         do i=isd,ied
            utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
            vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
         enddo
      enddo
  endif
  if ( is==1 .or. isd<npt ) then
      do j=max(npt,jsd),min(npy-npt,jed)
         do i=isd,npt-1
            utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
            vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
         enddo
      enddo
  endif
  if ( (ie+1)==npx .or. ied>=(npx-npt)) then
      do j=max(npt,jsd),min(npy-npt,jed)
         do i=npx-npt+1,ied
            utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
            vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
         enddo
      enddo
  endif

  endif

  do j=js-1-id,je+1+id
     do i=is-1-id,ie+1+id
        ua(i,j) = (utmp(i,j)-vtmp(i,j)*cosa_s(i,j)) * rsin2(i,j)
        va(i,j) = (vtmp(i,j)-utmp(i,j)*cosa_s(i,j)) * rsin2(i,j)
     enddo
  enddo

! Re-compute (utmp, vtmp) using (um,vm)
!----------
! Interior:
!----------
  do j=max(npt,js-1),min(npy-npt,je+1)
     do i=max(npt,isd),min(npx-npt,ied)
        utmp(i,j) = a2*(um(i,j-1)+um(i,j+2)) + a1*(um(i,j)+um(i,j+1))
     enddo
  enddo
  do j=max(npt,jsd),min(npy-npt,jed)
     do i=max(npt,is-1),min(npx-npt,ie+1)
        vtmp(i,j) = a2*(vm(i-1,j)+vm(i+2,j)) + a1*(vm(i,j)+vm(i+1,j))
     enddo
  enddo

!----------
! edges:
!----------
  if (grid_type < 3) then

  if ( js==1 .or. jsd<npt) then
      do j=jsd,npt-1
         do i=isd,ied
            utmp(i,j) = 0.5*(um(i,j) + um(i,j+1))
            vtmp(i,j) = 0.5*(vm(i,j) + vm(i+1,j))
         enddo
      enddo
  endif
  if ( (je+1)==npy .or. jed>=(npy-npt)) then
      do j=npy-npt+1,jed
         do i=isd,ied
            utmp(i,j) = 0.5*(um(i,j) + um(i,j+1))
            vtmp(i,j) = 0.5*(vm(i,j) + vm(i+1,j))
         enddo
      enddo
  endif
  if ( is==1 .or. isd<npt ) then
      do j=max(npt,jsd),min(npy-npt,jed)
         do i=isd,npt-1
            utmp(i,j) = 0.5*(um(i,j) + um(i,j+1))
            vtmp(i,j) = 0.5*(vm(i,j) + vm(i+1,j))
         enddo
      enddo
  endif
  if ( (ie+1)==npx .or. ied>=(npx-npt)) then
      do j=max(npt,jsd),min(npy-npt,jed)
         do i=npx-npt+1,ied
            utmp(i,j) = 0.5*(um(i,j) + um(i,j+1))
            vtmp(i,j) = 0.5*(vm(i,j) + vm(i+1,j))
         enddo
      enddo
  endif

  endif

! A -> C
!--------------
! Fix the edges
!--------------
! Xdir:
     if( sw_corner ) then
         do i=-2,0
            utmp(i,0) = -vtmp(0,1-i)
         enddo
     endif
     if( se_corner ) then
         do i=0,2
            utmp(npx+i,0) = vtmp(npx,i+1)
         enddo
     endif
     if( ne_corner ) then
         do i=0,2
            utmp(npx+i,npy) = -vtmp(npx,je-i)
         enddo
     endif
     if( nw_corner ) then
         do i=-2,0
            utmp(i,npy) = vtmp(0,je+i)
         enddo
     endif

  if (grid_type < 3) then
     ifirst = max(3,    is  )
     ilast  = min(npx-2,ie+1)
  else
     ifirst = is
     ilast  = ie+1
  endif
!---------------------------------------------
! 4th order interpolation for interior points:
!---------------------------------------------
     do j=js,je
        do i=ifirst,ilast
           uc(i,j) = a1*(utmp(i-1,j)+utmp(i,j))+a2*(utmp(i-2,j)+utmp(i+1,j))
        enddo
     enddo

     if (grid_type < 3) then

     if( is==1 ) then
        do j=js,je
! 3-pt extrapolation --------------------------------------------------
           uc(1,j) = ( t14*(utmp( 0,j)+utmp(1,j))    &
                     + t12*(utmp(-1,j)+utmp(2,j))    &
                     + t15*(utmp(-2,j)+utmp(3,j)) )*rsin_u(1,j)
           uc(2,j) = c1*utmp(3,j) + c2*utmp(2,j) + c3*utmp(1,j)
        enddo
     endif

     if( (ie+1)==npx ) then
        do j=js,je
           uc(npx-1,j) = c1*utmp(npx-3,j)+c2*utmp(npx-2,j)+c3*utmp(npx-1,j) 
! 3-pt extrapolation --------------------------------------------------------
           uc(npx,j) = (t14*(utmp(npx-1,j)+utmp(npx,j))+      &
                        t12*(utmp(npx-2,j)+utmp(npx+1,j))     &
                      + t15*(utmp(npx-3,j)+utmp(npx+2,j)))*rsin_u(npx,j)
        enddo
     endif

     endif

!------
! Ydir:
!------
     if( sw_corner ) then
         do j=-2,0
            vtmp(0,j) = -utmp(1-j,0)
         enddo
     endif
     if( nw_corner ) then
         do j=0,2
            vtmp(0,npy+j) = utmp(j+1,npy)
         enddo
     endif
     if( se_corner ) then
         do j=-2,0
            vtmp(npx,j) = utmp(ie+j,0)
         enddo
     endif
     if( ne_corner ) then
         do j=0,2
            vtmp(npx,npy+j) = -utmp(ie-j,npy)
         enddo
     endif

     if (grid_type < 3) then

     do j=js,je+1
      if ( j==1 ) then
        do i=is,ie
! 3-pt extrapolation -----------------------------------------
           vc(i,1) = (t14*(vtmp(i, 0)+vtmp(i,1))    &
                    + t12*(vtmp(i,-1)+vtmp(i,2))    &
                    + t15*(vtmp(i,-2)+vtmp(i,3)))*rsin_v(i,1)
        enddo
      elseif ( j==(npy-1) ) then
        do i=is,ie
           vc(i,j) = c1*vtmp(i,j-2) + c2*vtmp(i,j-1) + c3*vtmp(i,j)
        enddo
      elseif ( j==2 ) then
        do i=is,ie
!          vc(i,j) = c1*vtmp(i,j+1) + c2*vtmp(i,j) + c3*vtmp(i,j-1)
           vc(i,j) = c3*vtmp(i,j-1) + c2*vtmp(i,j) + c1*vtmp(i,j+1)
        enddo
      elseif ( j==npy ) then
        do i=is,ie
! 3-pt extrapolation --------------------------------------------------------
           vc(i,npy) = (t14*(vtmp(i,npy-1)+vtmp(i,npy))    &
                      + t12*(vtmp(i,npy-2)+vtmp(i,npy+1))  &
                      + t15*(vtmp(i,npy-3)+vtmp(i,npy+2)))*rsin_v(i,npy)
        enddo
      else
! 4th order interpolation for interior points:
        do i=is,ie
           vc(i,j) = a2*(vtmp(i,j-2)+vtmp(i,j+1))+a1*(vtmp(i,j-1)+vtmp(i,j))
        enddo
      endif
     enddo

    else
! 4th order interpolation:
       do j=js,je+1
          do i=is,ie
             vc(i,j) = a2*(vtmp(i,j-2)+vtmp(i,j+1))+a1*(vtmp(i,j-1)+vtmp(i,j))
          enddo
       enddo
    endif

 end subroutine d2a2c

 end module sw_core_mod
