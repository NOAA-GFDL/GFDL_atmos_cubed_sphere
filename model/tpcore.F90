module tpcore
!BOP
!
! !MODULE: tp_core --- A collection of routines to support FV transport
!
 use mp_mod,     only: is,js,ie,je, ng, isd,jsd,ied,jed
 use grid_utils, only: sw_corner, se_corner, ne_corner, nw_corner,   &
                       cx1, cx2, cy1, cy2
 use grid_tools, only: dxa, dya, grid_type

 implicit none

 private
 public fv_tp_2d, pert_ppm, copy_corners

 real, parameter:: r3 = 1./3.

#ifdef WAVE_FORM
! Suresh & Huynh scheme 2.2 (purtabation form)
! The wave-form is more diffusive than the default "5th" order
 real, parameter:: b1 =   0.0375
 real, parameter:: b2 =  -7./30.
 real, parameter:: b3 =  -23./120.
 real, parameter:: b4 =  13./30.
 real, parameter:: b5 = -11./240.
#else
! scheme 2.1: perturbation form
 real, parameter:: b1 =   1./30.
 real, parameter:: b2 = -13./60.
 real, parameter:: b3 = -13./60.
 real, parameter:: b4 =  0.45
 real, parameter:: b5 = -0.05
#endif
 real, parameter:: t11 = 27./28, t12=-13./28., t13=3./7.

!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------

CONTAINS

 subroutine fv_tp_2d(q, crx, cry, npx, npy, hord, fx, fy, &
                     xfx, yfx, area, ra_x, ra_y, uni_ppm, mfx, mfy)
   integer, intent(in):: npx, npy
   integer, intent(in)::hord
   logical, intent(IN) :: uni_ppm

   real, intent(in)::  crx(is:ie+1,jsd:jed)  !
   real, intent(in)::  xfx(is:ie+1,jsd:jed)  !
   real, intent(in)::  cry(isd:ied,js:je+1 )  !
   real, intent(in)::  yfx(isd:ied,js:je+1 )  !
   real, intent(in):: area(isd:ied,jsd:jed)
   real, intent(in):: ra_x(is:ie,jsd:jed)
   real, intent(in):: ra_y(isd:ied,js:je)
   real, intent(inout):: q(isd:ied,jsd:jed)  ! transported scalar
   real, intent(out)::fx(is:ie+1 ,js:je)    ! Flux in x ( E )
   real, intent(out)::fy(is:ie,   js:je+1 )    ! Flux in y ( N )
! optional Arguments:
   real, OPTIONAL, intent(in):: mfx(is:ie+1,js:je  )  ! Mass Flux X-Dir
   real, OPTIONAL, intent(in):: mfy(is:ie  ,js:je+1)  ! Mass Flux Y-Dir
! Local:
   real q_i(isd:ied,js:je)
   real q_j(is:ie,jsd:jed)
   real   fx2(is:ie+1,jsd:jed)
   real   fy2(isd:ied,js:je+1)
   real   fyy(isd:ied,js:je+1)
   real   fx1(is:ie+1)
   integer i, j

   call copy_corners(q, npx, npy, 2)
   call ytp(fy2, q, cry, hord, isd, ied, js, je, npx, npy, uni_ppm)

   do j=js,je+1
      do i=isd,ied
         fyy(i,j) = yfx(i,j) * fy2(i,j) 
      enddo
   enddo
   do j=js,je
      do i=isd,ied
         q_i(i,j) = (q(i,j)*area(i,j) + fyy(i,j)-fyy(i,j+1))/ra_y(i,j)
      enddo
  enddo
  call xtp(fx, q_i, crx(is,js), hord, is, ie, js, je, npx, npy, uni_ppm)

  call copy_corners(q, npx, npy, 1)
  call xtp(fx2, q, crx, hord, is, ie, jsd,jed, npx, npy, uni_ppm)

  do j=jsd,jed
     do i=is,ie+1
        fx1(i) =  xfx(i,j) * fx2(i,j)
     enddo
     do i=is,ie
        q_j(i,j) = (q(i,j)*area(i,j) + fx1(i)-fx1(i+1))/ra_x(i,j)
     enddo
  enddo
  call ytp(fy, q_j, cry, hord, is, ie, js, je, npx, npy, uni_ppm)

!----------------
! Flux averaging:
!----------------

   if ( present(mfx) .and. present(mfy) ) then
!---------------------------------
! For transport of pt and tracers
!---------------------------------
      do j=js,je
         do i=is,ie+1
            fx(i,j) = 0.5*(fx(i,j) + fx2(i,j)) * mfx(i,j)
         enddo
      enddo
      do j=js,je+1
         do i=is,ie
            fy(i,j) = 0.5*(fy(i,j) + fy2(i,j)) * mfy(i,j)
         enddo
      enddo
   else
!---------------------------------
! For transport of delp, vorticity
!---------------------------------
      do j=js,je
         do i=is,ie+1
            fx(i,j) = 0.5*(fx(i,j) + fx2(i,j)) * xfx(i,j)
         enddo
      enddo
      do j=js,je+1
         do i=is,ie
            fy(i,j) = 0.5*(fy(i,j) + fy2(i,j)) * yfx(i,j)
         enddo
      enddo
   endif

 end subroutine fv_tp_2d


 subroutine copy_corners(q, npx, npy, dir)
 integer, intent(in):: npx, npy, dir
 real, intent(inout):: q(isd:ied,jsd:jed)
 integer  i,j

 if ( dir == 1 ) then
! XDir:
    if ( sw_corner ) then
         do j=1-ng,0
            do i=1-ng,0
               q(i,j) = q(j,1-i)
            enddo
         enddo
    endif
    if ( se_corner ) then
         do j=1-ng,0
            do i=npx,npx+ng-1
               q(i,j) = q(npy-j,i-npx+1)
            enddo
         enddo
    endif
    if ( ne_corner ) then
         do j=npy,npy+ng-1
            do i=npx,npx+ng-1
               q(i,j) = q(j,2*npx-1-i)
            enddo
         enddo
    endif
    if ( nw_corner ) then
         do j=npy,npy+ng-1
            do i=1-ng,0
               q(i,j) = q(npy-j,i-1+npx)
            enddo
         enddo
    endif

 elseif ( dir == 2 ) then
! YDir:

    if ( sw_corner ) then
         do j=1-ng,0
            do i=1-ng,0
               q(i,j) = q(1-j,i)
            enddo
         enddo
    endif
    if ( se_corner ) then
         do j=1-ng,0
            do i=npx,npx+ng-1
               q(i,j) = q(npy+j-1,npx-i)
            enddo
         enddo
    endif
    if ( ne_corner ) then
         do j=npy,npy+ng-1
            do i=npx,npx+ng-1
               q(i,j) = q(2*npy-1-j,i)
            enddo
         enddo
    endif
    if ( nw_corner ) then
         do j=npy,npy+ng-1
            do i=1-ng,0
               q(i,j) = q(j+1-npx,npy-i)
            enddo
         enddo
    endif

 endif
      
 end subroutine copy_corners


 subroutine xtp(fx,  q,  c, iord, ifirst, ilast, jfirst, jlast,  &
                npx, npy, uniform_ppm)
   integer, intent(IN):: ifirst, ilast   !  X-Dir strip
   integer, intent(IN):: jfirst, jlast   !  Y-Dir strip
   integer, intent(IN):: npx, npy
   integer, intent(IN):: iord
   logical, intent(IN):: uniform_ppm
   real   , intent(in):: c(is :ie+1, jfirst:jlast)      ! Courant numbers
   real   , intent(in):: q(isd:ied,  jfirst:jlast)
   real,    intent(out):: fx(ifirst:ilast+1,jfirst:jlast)           
! Local:
   real   dm(ifirst-ng:ilast+ng,jfirst:jlast)
   integer i, j, iu

   if (iord==1) then
      do j=jfirst,jlast
         do i=ifirst,ilast+1
            iu = FLOOR(real(i) - c(i,j))
            fx(i,j) = q(iu,j)
         enddo
      enddo
   elseif (abs(iord)==2) then
      call xmist_2d(q, dm, ng, iord, ifirst, ilast, jfirst, jlast)

      do j=jfirst,jlast
         do i=ifirst,ilast+1
            iu = FLOOR(real(i) - c(i,j))
            fx(i,j) = q(iu,j) + (sign(1.,c(i,j)) - c(i,j))*dm(iu,j)
        enddo
      enddo
   else
      call fxppm(c, q, fx, iord, ifirst, ilast, jfirst, jlast, npx, npy, uniform_ppm)
   endif

 end subroutine xtp



 subroutine ytp(fy, q, c, jord, ifirst, ilast, jfirst, jlast,  &
                npx, npy, uniform_ppm)
 integer, intent(in) :: npx, npy
 integer, INTENT(IN) :: ifirst, ilast  !  X-Dir strip
 integer, INTENT(IN) :: jfirst, jlast  !  Y-Dir strip
 integer, intent(in):: jord
 logical, intent(IN) :: uniform_ppm
 real, intent(in)::   q(ifirst:ilast,jfirst-ng:jlast+ng) 
 real, intent(in)::   c(isd:ied,js:je+1 )  ! Courant number
 real, intent(out):: fy(ifirst:ilast,jfirst:jlast+1)     !  Flux
! !LOCAL VARIABLES:
 real   dm(ifirst:ilast,jfirst-ng:jlast+ng)
 integer i, j, n, jt

   if(jord==1) then
      do j=jfirst,jlast+1
         do i=ifirst,ilast
            jt = FLOOR(real(j) - c(i,j))
            fy(i,j) = q(i,jt)
         enddo
      enddo
   elseif(abs(jord)==2) then
      call ymist(q, dm, ng, jord, ifirst, ilast, jfirst, jlast)

      do j=jfirst,jlast+1
         do i=ifirst,ilast
            jt = FLOOR(real(j) - c(i,j))
            fy(i,j) = q(i,jt) + (sign(1.,c(i,j))-c(i,j))*dm(i,jt)
         enddo
      enddo
   else
      call fyppm(c, q, fy, jord, ifirst,ilast,jfirst,jlast,  &
                 npx, npy, uniform_ppm, dm)
   endif

 end subroutine ytp



 subroutine fxppm(c, q, flux, iord, ifirst, ilast, jfirst,   &
                  jlast,npx,npy, uniform_ppm)
! !INPUT PARAMETERS:
 integer, INTENT(IN) :: ifirst, ilast               !  X-Dir strip
 integer, INTENT(IN) :: jfirst, jlast               !  Y-Dir strip
 integer, INTENT(IN) :: iord
 integer, INTENT(IN) :: npx, npy
 logical, intent(IN) :: uniform_ppm
 real   , INTENT(IN) :: q(ifirst-ng:ilast+ng,jfirst:jlast)
 real   , INTENT(IN) :: c(ifirst   :ilast+1 ,jfirst:jlast) ! Courant   N (like FLUX)
! !OUTPUT PARAMETERS:
 real   , INTENT(OUT) :: flux(ifirst:ilast+1,jfirst:jlast) !  Flux
! Local
 real dm1(ifirst-2:ilast+2)
 real  al(ifirst-1:ilast+2)
 real  bl(ifirst-1:ilast+1)
 real  br(ifirst-1:ilast+1)
 real  dq(ifirst-3:ilast+2)
 real dl, dr, xt, pmp, lac, c1, qe
 integer i, j, is3, ie3, it

 if (iord<=4) then
     do j=jfirst,jlast

        if ( uniform_ppm ) then
             do i=ifirst-2,ilast+2
                xt = 0.25*(q(i+1,j) - q(i-1,j))
                dm1(i) = sign(min(abs(xt), max(q(i-1,j), q(i,j), q(i+1,j)) - q(i,j),  &
                                  q(i,j) - min(q(i-1,j), q(i,j), q(i+1,j))), xt)
             enddo
        else
! Variable grid size:
              do i=ifirst-3,ilast+2
                 dq(i) = q(i+1,j) - q(i,j)
             enddo
             do i=ifirst-2,ilast+2
                xt = cx1(i,j)*dq(i-1) + cx2(i,j)*dq(i)
                dm1(i) = sign(min(abs(xt), abs(dq(i-1)), abs(dq(i))), xt)
             enddo
        endif

        do i=ifirst-1,ilast+2
           al(i) = 0.5*(q(i-1,j)+q(i,j)) + r3*(dm1(i-1) - dm1(i))
        enddo

        do i=ifirst,ilast+1
          if(c(i,j)>0.) then
             xt = 2.*dm1(i-1)
             dl = sign(min(abs(xt), abs(al(i-1)-q(i-1,j))), xt)
             dr = sign(min(abs(xt), abs(al(i  )-q(i-1,j))), xt)
             flux(i,j) = q(i-1,j) + (1.-c(i,j))*(c(i,j)*(dl-dr) + dr)
          else
             xt = 2.*dm1(i)
             dl = sign(min(abs(xt), abs(al(i  )-q(i,j))), xt)
             dr = sign(min(abs(xt), abs(al(i+1)-q(i,j))), xt)
             flux(i,j) = q(i,j) - (1.+c(i,j))*(c(i,j)*(dl-dr) + dl)
          endif 
        enddo
     enddo
 elseif (iord==5) then
! PPM with Hunyh's 2nd constraint
     do j=jfirst,jlast
        do i=ifirst-3,ilast+2
           dq(i) = q(i+1,j) - q(i,j)
        enddo

        if ( uniform_ppm ) then
             do i=ifirst-2,ilast+2
                xt = 0.25*(q(i+1,j) - q(i-1,j))
                dm1(i) = sign(min(abs(xt), max(q(i-1,j), q(i,j), q(i+1,j)) - q(i,j),  &
                                  q(i,j) - min(q(i-1,j), q(i,j), q(i+1,j))), xt)
             enddo
        else
             do i=ifirst-2,ilast+2
                xt = cx1(i,j)*dq(i-1) + cx2(i,j)*dq(i)
                dm1(i) = sign(min(abs(xt), abs(dq(i-1)), abs(dq(i))), xt)
             enddo
        endif

        do i=ifirst-1,ilast+2
           al(i) = 0.5*(q(i-1,j)+q(i,j)) + r3*(dm1(i-1) - dm1(i))
        enddo

        do i=ifirst-1,ilast+1
           pmp = -2.*dq(i)
           lac = pmp + 1.5*dq(i+1)
           bl(i) = min(max(0., pmp, lac), max(al(i)-q(i,j), min(0.,pmp, lac)))
           pmp = 2.*dq(i-1)
           lac = pmp - 1.5*dq(i-2)
           br(i) = min(max(0., pmp, lac), max(al(i+1)-q(i,j), min(0.,pmp, lac)))
        enddo

        do i=ifirst,ilast+1
           if(c(i,j)>0.) then
              flux(i,j) = q(i-1,j) + (1.-c(i,j))*(br(i-1)-c(i,j)*(bl(i-1)+br(i-1)))
           else
              flux(i,j) = q(i,  j) + (1.+c(i,j))*(bl(i  )+c(i,j)*(bl(i  )+br(i  )))
           endif
        enddo
     enddo
 elseif ( iord==6 .or. iord==7 ) then
! Non-monotonic fast "5th order" scheme (not really 5th order)
     do j=jfirst,jlast
        do i=ifirst-1,ilast+1
              xt = b3*q(i,j)
           bl(i) = b5*q(i-2,j) + b4*q(i-1,j) + xt + b2*q(i+1,j) + b1*q(i+2,j)
           br(i) = b1*q(i-2,j) + b2*q(i-1,j) + xt + b4*q(i+1,j) + b5*q(i+2,j)
        enddo

        if ( iord == 7 ) then
!-- Huynh's 2nd constraint + simple mono limiter ---------------------------------
        do i=ifirst-2,ilast+3
           dq(i) = q(i,j) - q(i-1,j)
        enddo
        do i=ifirst-1,ilast+1
!Left:
           dl = -sign(min(abs(bl(i)), abs(dq(i))), dq(i))
           pmp = -2.*dq(i+1)
           lac = pmp + 1.5*dq(i+2)
           bl(i) = min(max(0., pmp, lac), max(dl, min(0.,pmp, lac)))
!Right:
           dr = sign(min(abs(br(i)), abs(dq(i+1))), dq(i+1))
           pmp = 2.*dq(i)
           lac = pmp - 1.5*dq(i-1)
           br(i) = min(max(0., pmp, lac), max(dr, min(0.,pmp, lac)))
        enddo
!---------------------------------------------------------------------------------
        endif

        do i=ifirst,ilast+1
           if(c(i,j)>0.) then
              flux(i,j) = q(i-1,j) + (1.-c(i,j))*(br(i-1)-c(i,j)*(bl(i-1)+br(i-1)))
           else
              flux(i,j) = q(i,  j) + (1.+c(i,j))*(bl(i  )+c(i,j)*(bl(i  )+br(i  )))
           endif
        enddo
     enddo
 elseif( iord==9 .or. iord==10 ) then
! Based on 5; with absolute upwind algorithm at edges
     is3 = max(3,is-1);   ie3 = min(npx-3,ie+1)

     do j=jfirst,jlast
        do i=is-3,ie+2
           dq(i) = q(i+1,j) - q(i,j)
        enddo

        if ( uniform_ppm ) then
             do i=is-2,ie+2
                xt = 0.25*(q(i+1,j) - q(i-1,j))
                dm1(i) = sign(min(abs(xt), max(q(i-1,j), q(i,j), q(i+1,j)) - q(i,j),  &
                                  q(i,j) - min(q(i-1,j), q(i,j), q(i+1,j))), xt)
             enddo
        else
             do i=is-2,ie+2
                xt = cx1(i,j)*dq(i-1) + cx2(i,j)*dq(i)
                dm1(i) = sign(min(abs(xt), abs(dq(i-1)), abs(dq(i))), xt)
             enddo
        endif

        if (grid_type < 3) then

           do i=is3,min(npx-2,ie+2)
              al(i) = 0.5*(q(i-1,j)+q(i,j)) + r3*(dm1(i-1)-dm1(i))
           enddo

           do i=is3, ie3
              pmp = -2.*dq(i)
              lac = pmp + 1.5*dq(i+1)
              bl(i) = min(max(0., pmp, lac), max(al(i  )-q(i,j), min(0.,pmp, lac)))
              pmp = 2.*dq(i-1)
              lac = pmp - 1.5*dq(i-2)
              br(i) = min(max(0., pmp, lac), max(al(i+1)-q(i,j), min(0.,pmp, lac)))
           enddo

!--------------
! fix the edges
!--------------
           if ( is==1 ) then
              br(2) = al(3) - q(2,j)
!             xt = 0.75*(q(0,j)+q(1,j)) - 0.25*(q(-1,j)+q(2,j))
              xt = 0.5*((2.*dxa(1,j)+dxa(2,j))*(q(0,j)+q(1,j))   &
                 - dxa(1,j)*(q(-1,j)+q(2,j)))/ ( dxa(1,j)+dxa(2,j))
!             xt = 6./7.*(q(0,j)+q(1,j)) - 13./28.*(q(-1,j)+q(2,j))  &
!                + 3./28.*(q(-2,j)+q(3,j))
#ifdef EDGE_MONO
              if ( (j>-2 .and. j<1) .or. (j>(npy-3) .and. j<npy) ) then
                 xt = min( xt, max(q(0,j), q(1,j), q(0,j-1), q(1,j-1)) )
                 xt = max( xt, min(q(0,j), q(1,j), q(0,j-1), q(1,j-1)) )
              elseif ( (j<3 .and.j>0) .or. (j>(npy-1).and.j<(npy+2)) ) then
                 xt = min( xt, max(q(0,j), q(1,j), q(0,j+1), q(1,j+1)) )
                 xt = max( xt, min(q(0,j), q(1,j), q(0,j+1), q(1,j+1)) )
              endif
#endif
              
              bl(1) = xt - q(1,j)
              br(0) = xt - q(0,j)
!
              xt = 4./7.*dm1(-1) - 11./14.*dq(-1) + q(0,j)
#ifdef S_MONO
              xt = max( xt, min(q(-1,j),q(0,j)) )
              xt = min( xt, max(q(-1,j),q(0,j)) )
#endif

              bl(0) = xt - q(0,j)
              xt = 3./14.*q(1,j) + 11./14.*q( 2,j) - 4./7.*dm1( 2)
#ifdef S_MONO
              xt = max( xt, min(q(1,j),q(2,j)) )
              xt = min( xt, max(q(1,j),q(2,j)) )
#endif

              br(1) = xt - q(1,j)
              bl(2) = xt - q(2,j)
              if(iord==9) call pert_ppm(3, q(0,j), bl(0), br(0), 1)
           endif

           if ( (ie+1)==npx ) then
              bl(npx-2) = al(npx-2) - q(npx-2,j)
!             xt = 0.75*(q(npx-1,j)+q(npx,j)) - 0.25*(q(npx-2,j)+q(npx+1,j))
              xt = 0.5*( (2.*dxa(npx-1,j)+dxa(npx-2,j))*(q(npx-1,j)+q(npx,j))   &
                 - dxa(npx-1,j)*(q(npx-2,j)+q(npx+1,j)))/( dxa(npx-1,j)+dxa(npx-2,j))
!             xt = 6./7.*(q(npx-1,j)+q(npx,j)) - 13./28.*(q(npx-2,j)+q(npx+1,j))  &
!                + 3./28.*(q(npx-3,j)+q(npx+2,j))
#ifdef EDGE_MONO
              if ( (j>-2 .and. j<1) .or. (j>(npy-3) .and. j<npy) ) then
                 xt = min( xt, max(q(npx-1,j), q(npx,j), q(npx-1,j-1), q(npx,j-1)))
                 xt = max( xt, min(q(npx-1,j), q(npx,j), q(npx-1,j-1), q(npx,j-1)))
              elseif ( (j<3 .and.j>0) .or. (j>(npy-1).and.j<(npy+2)) ) then
                 xt = min( xt, max(q(npx-1,j), q(npx,j), q(npx-1,j+1), q(npx,j+1)))
                 xt = max( xt, min(q(npx-1,j), q(npx,j), q(npx-1,j+1), q(npx,j+1)))
              endif
#endif

              br(npx-1) = xt - q(npx-1,j)
              bl(npx  ) = xt - q(npx  ,j)
!             br(npx) = 11./14.*dq(npx) - 4./7.*dm1(npx+1)
              xt = 11./14.*dq(npx) - 4./7.*dm1(npx+1) + q(npx,j)
#ifdef S_MONO
              xt = min( xt, max(q(npx,j), q(npx+1,j)) )
              xt = max( xt, min(q(npx,j), q(npx+1,j)) )
#endif

              br(npx) = xt - q(npx,j)
              xt = 3./14.*q(npx-1,j) + 11./14.*q(npx-2,j) + 4./7.*dm1(npx-2)
#ifdef S_MONO
              xt = min( xt, max(q(npx-2,j), q(npx-1,j)) )
              xt = max( xt, min(q(npx-2,j), q(npx-1,j)) )
#endif

              br(npx-2) = xt - q(npx-2,j)
              bl(npx-1) = xt - q(npx-1,j)
              if(iord==9) call pert_ppm(3, q(npx-2,j), bl(npx-2), br(npx-2), 1)
           endif
        else

           do i=is-1,ie+2
              al(i) = 0.5*(q(i-1,j)+q(i,j)) + r3*(dm1(i-1)-dm1(i))
           enddo
           
           do i=is-1, ie+1
              pmp = -2.*dq(i)
              lac = pmp + 1.5*dq(i+1)
              bl(i) = min(max(0., pmp, lac), max(al(i  )-q(i,j), min(0.,pmp, lac)))
              pmp = 2.*dq(i-1)
              lac = pmp - 1.5*dq(i-2)
              br(i) = min(max(0., pmp, lac), max(al(i+1)-q(i,j), min(0.,pmp, lac)))
           enddo
           
        endif

        do i=is,ie+1
#ifdef INTEL_OPT
             c1 = c(i,j)
             if( c1>0. ) then
                it = i-1
                qe = br(i-1) 
             else
                it = i
                qe = bl(i) 
             endif
             c1 = -abs(c1)
             flux(i,j) = q(it,j) + (1.+c1)*( qe + c1*(bl(it)+br(it)) )
#else
             if( c(i,j)>0. ) then
                flux(i,j) = q(i-1,j) + (1.-c(i,j))*(br(i-1)-c(i,j)*(bl(i-1)+br(i-1)))
             else
                flux(i,j) = q(i,  j) + (1.+c(i,j))*(bl(i  )+c(i,j)*(bl(i  )+br(i  )))
             endif
#endif
          enddo
       enddo
    else
!------------------------------
! For positive definite tracers:
!------------------------------
! iord=11: PPM mono constraint (Lin 2004)
! iord=12: Huynh 2nd constraint (Lin 2004) + positive definite (Lin & Rood 1996)
! iord>12: positive definite only (Lin & Rood 1996)

       do j=jfirst,jlast

          do i=ifirst-2,ilast+2
             xt = 0.25*(q(i+1,j) - q(i-1,j))
             dm1(i) = sign(min(abs(xt), max(q(i-1,j), q(i,j), q(i+1,j)) - q(i,j),  &
                               q(i,j) - min(q(i-1,j), q(i,j), q(i+1,j))), xt)
          enddo

          if (grid_type < 3) then

             is3 = max(3,is-1);   ie3 = min(npx-3,ie+1)

             do i=is3,min(npx-2,ie+2)
                al(i) = 0.5*(q(i-1,j)+q(i,j)) + r3*(dm1(i-1)-dm1(i))
             enddo

             if ( iord ==11 ) then
                do i=is3,ie3
                   xt = 2.*dm1(i)
                   bl(i) =-sign(min(abs(xt), abs(al(i)  -q(i,j))), xt)
                   br(i) = sign(min(abs(xt), abs(al(i+1)-q(i,j))), xt)
                enddo
             elseif( iord==12 ) then
                do i=is-3,ie+2
                   dq(i) = q(i+1,j) - q(i,j)
                enddo
                do i=is3,ie3
                   pmp = -2.*dq(i)
                   lac = pmp + 1.5*dq(i+1)
                   bl(i) = min(max(0., pmp, lac), max(al(i  )-q(i,j), min(0.,pmp, lac)))
                   pmp = 2.*dq(i-1)
                   lac = pmp - 1.5*dq(i-2)
                   br(i) = min(max(0., pmp, lac), max(al(i+1)-q(i,j), min(0.,pmp, lac)))
                enddo
             else
                do i=is3,ie3
                   bl(i) = al(i  ) - q(i,j)
                   br(i) = al(i+1) - q(i,j)
                enddo
             endif

! Positive definite constraint:
             if(iord/=11) call pert_ppm(ie3-is3+1, q(is3,j), bl(is3), br(is3), 0)

!--------------
! fix the edges
!--------------
             if ( is==1 ) then
                br(2) = al(3) - q(2,j)
!               xt = t11*(q(0,j)+q(1,j)) + t12*(q(-1,j)+q(2,j)) + t13*(dm1(2)-dm1(-1))
!!!             xt = 0.75*(q(0,j)+q(1,j)) - 0.25*(q(-1,j)+q(2,j))
                xt = 0.5*( (2.*dxa(1,j)+dxa(2,j))*(q(0,j)+q(1,j))  &
                   - dxa(1,j)*(q(-1,j)+q(2,j)) ) / ( dxa(1,j)+dxa(2,j) )
#ifdef EDGE_MONO
                if ( (j>-2 .and. j<1) .or. (j>(npy-3) .and. j<npy) ) then
                   xt = min( xt, max(q(0,j), q(1,j), q(0,j-1), q(1,j-1)) )
                   xt = max( xt, min(q(0,j), q(1,j), q(0,j-1), q(1,j-1)) )
                elseif ( (j<3 .and.j>0) .or. (j>(npy-1).and.j<(npy+2)) ) then
                   xt = min( xt, max(q(0,j), q(1,j), q(0,j+1), q(1,j+1)) )
                   xt = max( xt, min(q(0,j), q(1,j), q(0,j+1), q(1,j+1)) )
                else
                   xt = max(0., xt)
                endif
#else
                xt = max(0., xt)
#endif
                bl(1) = xt - q(1,j)
                br(0) = xt - q(0,j)
                xt = 4./7.*dm1(-1) + 11./14.*q(-1,j) + 3./14.*q(0,j)
#ifdef S_MONO
                xt = min( xt, max(q(-1,j), q(0,j)) )
                xt = max( xt, min(q(-1,j), q(0,j)) )
#else
                xt = max(0., xt)
#endif
                bl(0) =  xt - q(0,j)
                xt = 3./14.*q(1,j) + 11./14.*q(2,j) - 4./7.*dm1(2)
#ifdef S_MONO
                xt = min( xt, max(q(1,j), q(2,j)) )
                xt = max( xt, min(q(1,j), q(2,j)) )
#else
                xt = max(0., xt)
#endif
                br(1) = xt - q(1,j)
                bl(2) = xt - q(2,j)
                call pert_ppm(3, q(0,j), bl(0), br(0), 1)
             endif

             if ( (ie+1)==npx ) then
                bl(npx-2) = al(npx-2) - q(npx-2,j)
!               xt = t11*(q(npx-1,j)+q(npx,j)) + t12*(q(npx-2,j)+q(npx+1,j))   &
!                  + t13*(dm1(npx+1)-dm1(npx-2))
!!!             xt = 0.75*(q(npx-1,j)+q(npx,j)) - 0.25*(q(npx-2,j)+q(npx+1,j))
                xt = 0.5*((2.*dxa(npx-1,j)+dxa(npx-2,j))*(q(npx-1,j)+q(npx,j)) -   &
                     dxa(npx-1,j)*(q(npx-2,j)+q(npx+1,j)) )  &
                 / ( dxa(npx-1,j)+dxa(npx-2,j) )
#ifdef EDGE_MONO
                if ( (j>-2 .and. j<1) .or. (j>(npy-3) .and. j<npy) ) then
                   xt = min( xt, max(q(npx-1,j), q(npx,j), q(npx-1,j-1), q(npx,j-1)))
                   xt = max( xt, min(q(npx-1,j), q(npx,j), q(npx-1,j-1), q(npx,j-1)))
                elseif ( (j<3 .and.j>0) .or. (j>(npy-1).and.j<(npy+2)) ) then
                   xt = min( xt, max(q(npx-1,j), q(npx,j), q(npx-1,j+1), q(npx,j+1)))
                   xt = max( xt, min(q(npx-1,j), q(npx,j), q(npx-1,j+1), q(npx,j+1)))
                else
                   xt = max(0., xt)
                endif
#else
                xt = max(0., xt)
#endif
                br(npx-1) = xt - q(npx-1,j)
                bl(npx  ) = xt - q(npx  ,j)
!               br(npx) = 11./14.*q(npx+1,j) + 3./14.*q(npx,j) - 4./7.*dm1(npx+1)
                xt = 11./14.*q(npx+1,j) + 3./14.*q(npx,j) - 4./7.*dm1(npx+1)
#ifdef S_MONO
                xt = min( xt, max(q(npx,j), q(npx+1,j)) )
                xt = max( xt, min(q(npx,j), q(npx+1,j)) )
#else
                xt = max(0., xt)
#endif
                br(npx) = xt - q(npx,j)
                xt = 3./14.*q(npx-1,j) + 11./14.*q(npx-2,j) + 4./7.*dm1(npx-2)
#ifdef S_MONO
                xt = min( xt, max(q(npx-2,j), q(npx-1,j)) )
                xt = max( xt, min(q(npx-2,j), q(npx-1,j)) )
#else
                xt = max(0., xt)
#endif
                br(npx-2) = xt - q(npx-2,j)
                bl(npx-1) = xt - q(npx-1,j)
                call pert_ppm(3, q(npx-2,j), bl(npx-2), br(npx-2), 1)
             endif
          else

             do i=is-1,ie+2
                al(i) = 0.5*(q(i-1,j)+q(i,j)) + r3*(dm1(i-1)-dm1(i))
             enddo
             
             if ( iord ==11 ) then
                do i=is-1,ie+1
                   xt = 2.*dm1(i)
                   bl(i) =-sign(min(abs(xt), abs(al(i)  -q(i,j))), xt)
                   br(i) = sign(min(abs(xt), abs(al(i+1)-q(i,j))), xt)
                enddo
             elseif( iord==12 ) then
                do i=is-3,ie+2
                   dq(i) = q(i+1,j) - q(i,j)
                enddo
                do i=is-1,ie+1
                   pmp = -2.*dq(i)
                   lac = pmp + 1.5*dq(i+1)
                   bl(i) = min(max(0., pmp, lac), max(al(i  )-q(i,j), min(0.,pmp, lac)))
                   pmp = 2.*dq(i-1)
                   lac = pmp - 1.5*dq(i-2)
                   br(i) = min(max(0., pmp, lac), max(al(i+1)-q(i,j), min(0.,pmp, lac)))
                enddo
             else
                do i=is-1,ie+1
                   bl(i) = al(i  ) - q(i,j)
                   br(i) = al(i+1) - q(i,j)
                enddo
             endif

! Positive definite constraint:
             if(iord/=11) call pert_ppm(ie-is+1, q(is-1,j), bl(is-1), br(is-1), 0)

          endif
          
          do i=is,ie+1
             if( c(i,j)>0. ) then
                flux(i,j) = q(i-1,j) + (1.-c(i,j))*(br(i-1)-c(i,j)*(bl(i-1)+br(i-1)))
             else
                flux(i,j) = q(i,  j) + (1.+c(i,j))*(bl(i  )+c(i,j)*(bl(i  )+br(i  )))
             endif
          enddo
       enddo
    endif

 end subroutine fxppm



 subroutine fyppm(c,  q,  flux, jord, ifirst, ilast, jfirst, jlast,    &
                  npx, npy, uniform_ppm, dm)
 integer, INTENT(IN) :: ifirst, ilast               !  X-Dir strip
 integer, INTENT(IN) :: jfirst, jlast               !  Y-Dir strip
 integer, INTENT(IN) :: jord
 integer, INTENT(IN) :: npx, npy
 logical, intent(IN) :: uniform_ppm
 real   , INTENT(IN) :: q(ifirst:ilast,jfirst-ng:jlast+ng)
 real   , intent(in) :: c(isd:ied,js:je+1 )  ! Courant number
 real   , INTENT(OUT):: flux(ifirst:ilast,jfirst:jlast+1)   !  Flux
 real   , INTENT(OUT)::   dm(ifirst:ilast,jfirst-ng:jlast+ng)
! Local:
 real al(ifirst:ilast,jfirst-1:jlast+2)
 real bl(ifirst:ilast,jfirst-1:jlast+1)
 real br(ifirst:ilast,jfirst-1:jlast+1)
 real dq(ifirst:ilast,jfirst-3:jlast+2)
 real xt, dl, dr, pmp, lac, c1, qe
 integer i, j, js3, je3, jt

 if (jord<=4) then

     if ( uniform_ppm ) then
          do j=jfirst-2,jlast+2
             do i=ifirst,ilast
                xt = 0.25*(q(i,j+1) - q(i,j-1))
                dm(i,j) = sign(min(abs(xt), max(q(i,j-1), q(i,j), q(i,j+1)) - q(i,j),   &
                                   q(i,j) - min(q(i,j-1), q(i,j), q(i,j+1))), xt)
             enddo
          enddo
     else
! Variable grid size:
         do j=jfirst-3,jlast+2
            do i=ifirst,ilast
               dq(i,j) = q(i,j+1) - q(i,j)
            enddo
         enddo
         do j=jfirst-2,jlast+2
            do i=ifirst,ilast
               xt = cy1(i,j)*dq(i,j-1) + cy2(i,j)*dq(i,j)
               dm(i,j) = sign(min(abs(xt), abs(dq(i,j-1)), abs(dq(i,j))), xt)
            enddo
         enddo
     endif

   do j=jfirst-1,jlast+2
      do i=ifirst,ilast
         al(i,j) = 0.5*(q(i,j-1)+q(i,j)) + r3*(dm(i,j-1) - dm(i,j))
      enddo
   enddo

   do j=jfirst,jlast+1
      do i=ifirst,ilast
         if(c(i,j)>0.) then
            xt = 2.*dm(i,j-1)
            dl = sign(min(abs(xt), abs(al(i,j-1)-q(i,j-1))), xt)
            dr = sign(min(abs(xt), abs(al(i,j)-q(i,j-1))),   xt)
            flux(i,j) = q(i,j-1) + (1.-c(i,j))*(c(i,j)*(dl-dr)+dr)
         else
            xt = 2.*dm(i,j)
            dl = sign(min(abs(xt), abs(al(i,j)-q(i,j))),   xt)
            dr = sign(min(abs(xt), abs(al(i,j+1)-q(i,j))), xt)
            flux(i,j) = q(i,j) - (1.+c(i,j))*(c(i,j)*(dl-dr)+dl)
         endif
      enddo
   enddo
 elseif (jord==5) then
! PPM with Hunyh's 2nd constraint

   do j=jfirst-3, jlast+2
      do i=ifirst,ilast
         dq(i,j) = q(i,j+1) - q(i,j)
      enddo
   enddo

   if ( uniform_ppm ) then
        do j=jfirst-2,jlast+2
           do i=ifirst,ilast
              xt = 0.25*(q(i,j+1) - q(i,j-1))
              dm(i,j) = sign(min(abs(xt), max(q(i,j-1), q(i,j), q(i,j+1)) - q(i,j),   &
                                 q(i,j) - min(q(i,j-1), q(i,j), q(i,j+1))), xt)
           enddo
        enddo
   else
        do j=jfirst-2,jlast+2
           do i=ifirst,ilast
              xt = cy1(i,j)*dq(i,j-1) + cy2(i,j)*dq(i,j)
              dm(i,j) = sign(min(abs(xt), abs(dq(i,j-1)), abs(dq(i,j))), xt)
           enddo
        enddo
   endif

   do j=jfirst-1,jlast+2
      do i=ifirst,ilast
         al(i,j) = 0.5*(q(i,j-1)+q(i,j)) + r3*(dm(i,j-1) - dm(i,j))
      enddo
   enddo

   do j=jfirst-1,jlast+1
      do i=ifirst,ilast
            pmp = -2.*dq(i,j) 
            lac = pmp + 1.5*dq(i,j+1)
            bl(i,j) = min(max(0.,pmp,lac), max(al(i,j)-q(i,j), min(0.,pmp,lac)))
            pmp = 2.*dq(i,j-1)
            lac = pmp - 1.5*dq(i,j-2)
            br(i,j) = min(max(0.,pmp,lac), max(al(i,j+1)-q(i,j), min(0.,pmp,lac)))
      enddo
   enddo

   do j=jfirst,jlast+1
      do i=ifirst,ilast
         if(c(i,j)>0.) then
            flux(i,j) = q(i,j-1) + (1.-c(i,j))*(br(i,j-1)-c(i,j)*(bl(i,j-1)+br(i,j-1)))
         else
            flux(i,j) = q(i,j  ) + (1.+c(i,j))*(bl(i,j  )+c(i,j)*(bl(i,j  )+br(i,j  )))
         endif
      enddo
   enddo
 elseif( jord==6 .or. jord==7 ) then
! Non-monotonic "5th order" scheme (not really 5th order)
   do j=jfirst-1,jlast+1
      do i=ifirst,ilast
         xt = b3*q(i,j)
         bl(i,j) = b5*q(i,j-2) + b4*q(i,j-1) + xt + b2*q(i,j+1) + b1*q(i,j+2)
         br(i,j) = b1*q(i,j-2) + b2*q(i,j-1) + xt + b4*q(i,j+1) + b5*q(i,j+2)
      enddo
   enddo

!-- Huynh's 2nd constraint + simple mono limiter ---------------------------------
   if ( jord==7 ) then
   do j=jfirst-2,jlast+3
      do i=ifirst,ilast
         dq(i,j) = q(i,j) - q(i,j-1)
      enddo
   enddo
   do j=jfirst-1,jlast+1
      do i=ifirst,ilast
! 1st constraint:
          dl = -sign(min(abs(bl(i,j)), abs(dq(i,j))), dq(i,j))
         pmp = -2.*dq(i,j+1)
         lac = pmp + 1.5*dq(i,j+2)
         bl(i,j) = min(max(0.,pmp, lac), max(dl,  min(0.,pmp, lac)))
! 1st constraint:
          dr = sign(min(abs(br(i,j)), abs(dq(i,j+1))), dq(i,j+1))
         pmp = 2.*dq(i,j)
         lac = pmp - 1.5*dq(i,j-1)
         br(i,j) =  min(max(0.,pmp, lac), max(dr,  min(0.,pmp, lac)))
      enddo
   enddo
   endif
!---------------------------------------------------------------------------------

   do j=jfirst,jlast+1
      do i=ifirst,ilast
         if(c(i,j)>0.) then
            flux(i,j) = q(i,j-1) + (1.-c(i,j))*(br(i,j-1)-c(i,j)*(bl(i,j-1)+br(i,j-1)))
         else
            flux(i,j) = q(i,j  ) + (1.+c(i,j))*(bl(i,j  )+c(i,j)*(bl(i,j  )+br(i,j  )))
         endif
      enddo
   enddo

 elseif( jord==9 .or. jord==10 ) then
! Based on scheme-5

   do j=js-3,je+2
      do i=ifirst,ilast
         dq(i,j) = q(i,j+1) - q(i,j)
      enddo
   enddo

   if ( uniform_ppm ) then
        do j=js-2,je+2
           do i=ifirst,ilast
              xt = 0.25*(q(i,j+1) - q(i,j-1))
              dm(i,j) = sign(min(abs(xt), max(q(i,j-1), q(i,j), q(i,j+1)) - q(i,j),   &
                                 q(i,j) - min(q(i,j-1), q(i,j), q(i,j+1))), xt)
           enddo
        enddo
   else
        do j=js-2,je+2
           do i=ifirst,ilast
              xt = cy1(i,j)*dq(i,j-1) + cy2(i,j)*dq(i,j)
              dm(i,j) = sign(min(abs(xt), abs(dq(i,j-1)), abs(dq(i,j))), xt)
           enddo
        enddo
   endif

   if (grid_type < 3) then

      do j=max(3,js-1),min(npy-2,je+2)
         do i=ifirst,ilast
            al(i,j) = 0.5*(q(i,j-1)+q(i,j)) + r3*(dm(i,j-1) - dm(i,j))
         enddo
      enddo
      
      do j=max(3,js-1),min(npy-3,je+1)
         do i=ifirst,ilast
            pmp = -2.*dq(i,j) 
            lac = pmp + 1.5*dq(i,j+1)
            bl(i,j) = min(max(0.,pmp,lac), max(al(i,j  )-q(i,j), min(0.,pmp,lac)))
            pmp = 2.*dq(i,j-1)
            lac = pmp - 1.5*dq(i,j-2)
            br(i,j) = min(max(0.,pmp,lac), max(al(i,j+1)-q(i,j), min(0.,pmp,lac)))
         enddo
      enddo

!--------------
! Fix the edges:
!--------------
      if( js==1 ) then
         do i=ifirst,ilast
            br(i,2) = al(i,3) - q(i,2)
!           xt = 6./7.*(q(i,0)+q(i,1)) - 13./28.*(q(i,-1)+q(i,2))  &
!              + 3./28.*(q(i,-2)+q(i,3))
!           xt = 0.75*(q(i,0)+q(i,1)) - 0.25*(q(i,-1)+q(i,2))
            xt = 0.5*((2.*dya(i,1)+dya(i,2))*(q(i,0)+q(i,1))   &
               -dya(i,1)*(q(i,-1)+q(i,2))) / ( dya(i,1)+dya(i,2) )
#ifdef EDGE_MONO
            if ( (i>-2 .and. i<1) .or. (i>(npx-3) .and. i<npx) ) then
               xt = min( xt, max(q(i,0), q(i,1), q(i-1,0), q(i-1,1)) )
               xt = max( xt, min(q(i,0), q(i,1), q(i-1,0), q(i-1,1)) )
            elseif ( (i<3 .and.i>0) .or. (i>(npx-1).and.i<(npx+2)) ) then
               xt = min( xt, max(q(i,0), q(i,1), q(i+1,0), q(i+1,1)) )
               xt = max( xt, min(q(i,0), q(i,1), q(i+1,0), q(i+1,1)) )
            endif
#endif

            bl(i,1) = xt - q(i,1)
            br(i,0) = xt - q(i,0)

!           bl(i,0) = 4./7.*dm(i,-1) - 11./14.*dq(i,-1)
            xt = 4./7.*dm(i,-1) - 11./14.*dq(i,-1) + q(i,0)
#ifdef S_MONO
            xt = min( xt, max(q(i,-1), q(i,0)) )
            xt = max( xt, min(q(i,-1), q(i,0)) )
#endif
            bl(i,0) = xt - q(i,0)
            xt = 3./14.*q(i,1) + 11./14.*q(i,2) - 4./7.*dm(i,2)
#ifdef S_MONO
            xt = min( xt, max(q(i,1), q(i,2)) )
            xt = max( xt, min(q(i,1), q(i,2)) )
#endif
            br(i,1) = xt - q(i,1)
            bl(i,2) = xt - q(i,2)
         enddo
         if ( jord==9 ) then
            do j=0,2
               call pert_ppm(ilast-ifirst+1, q(ifirst,j), bl(ifirst,j), br(ifirst,j), 1)
            enddo
         endif
      endif

      if( (je+1)==npy ) then
         do i=ifirst,ilast
            bl(i,npy-2) = al(i,npy-2) - q(i,npy-2)
!           xt = 6./7.*(q(i,npy-1)+q(i,npy))-13./28.*(q(i,npy-2)+q(i,npy+1))  &
!              +3./28.*(q(i,npy-3)+q(i,npy+2))
!           xt = 0.75*(q(i,npy-1)+q(i,npy)) - 0.25*(q(i,npy-2)+q(i,npy+1))
            xt = 0.5*((2.*dya(i,npy-1)+dya(i,npy-2))*(q(i,npy-1)+q(i,npy))  &
               -dya(i,npy-1)*(q(i,npy-2)+q(i,npy+1)))/(dya(i,npy-1)+dya(i,npy-2))
#ifdef EDGE_MONO
            if ( (i>-2 .and. i<1) .or. (i>(npx-3) .and. i<npx) ) then
               xt = min( xt, max(q(i,npy-1), q(i,npy), q(i-1,npy-1), q(i-1,npy)) )
               xt = max( xt, min(q(i,npy-1), q(i,npy), q(i-1,npy-1), q(i-1,npy)) )
            elseif ( (i<3 .and.i>0) .or. (i>(npx-1).and.i<(npx+2)) ) then
               xt = min( xt, max(q(i,npy-1), q(i,npy), q(i+1,npy-1), q(i+1,npy)) )
               xt = max( xt, min(q(i,npy-1), q(i,npy), q(i+1,npy-1), q(i+1,npy)) )
            endif
#endif
            br(i,npy-1) = xt - q(i,npy-1)
            bl(i,npy  ) = xt - q(i,npy)
!           br(i,npy) = 11./14.*dq(i,npy) - 4./7.*dm(i,npy+1)
            xt = 11./14.*dq(i,npy) - 4./7.*dm(i,npy+1) + q(i,npy)
#ifdef S_MONO
            xt = min( xt, max( q(i,npy), q(i,npy+1)) )
            xt = max( xt, min( q(i,npy), q(i,npy+1)) )
#endif
            br(i,npy) = xt - q(i,npy)
            xt = 3./14.*q(i,npy-1) + 11./14.*q(i,npy-2) + 4./7.*dm(i,npy-2)
#ifdef S_MONO
            xt = min( xt, max( q(i,npy-2), q(i,npy-1)) )
            xt = max( xt, min( q(i,npy-2), q(i,npy-1)) )
#endif
            br(i,npy-2) = xt - q(i,npy-2)
            bl(i,npy-1) = xt - q(i,npy-1)
         enddo
         if ( jord==9 ) then
            do j=npy-2,npy
               call pert_ppm(ilast-ifirst+1, q(ifirst,j), bl(ifirst,j), br(ifirst,j), 1)
            enddo
         endif
      endif

   else

      do j=js-1,je+2
         do i=ifirst,ilast
            al(i,j) = 0.5*(q(i,j-1)+q(i,j)) + r3*(dm(i,j-1) - dm(i,j))
         enddo
      enddo
      
      do j=js-1,je+1
         do i=ifirst,ilast
            pmp = -2.*dq(i,j) 
            lac = pmp + 1.5*dq(i,j+1)
            bl(i,j) = min(max(0.,pmp,lac), max(al(i,j  )-q(i,j), min(0.,pmp,lac)))
            pmp = 2.*dq(i,j-1)
            lac = pmp - 1.5*dq(i,j-2)
            br(i,j) = min(max(0.,pmp,lac), max(al(i,j+1)-q(i,j), min(0.,pmp,lac)))
         enddo
      enddo

   endif

   do j=js,je+1
      do i=ifirst,ilast
#ifdef INTEL_OPT
         c1 = c(i,j)
         if( c1>0. ) then
             jt = j-1
             qe = br(i,j-1) 
         else
             jt = j
             qe = bl(i,j) 
         endif
         c1 = -abs(c1)
         flux(i,j) = q(i,jt) + (1.+c1)*( qe + c1*(bl(i,jt)+br(i,jt)) )
#else
         if( c(i,j)>0. ) then
            flux(i,j) = q(i,j-1) + (1.-c(i,j))*(br(i,j-1)-c(i,j)*(bl(i,j-1)+br(i,j-1)))
         else
            flux(i,j) = q(i,j  ) + (1.+c(i,j))*(bl(i,j  )+c(i,j)*(bl(i,j  )+br(i,j  )))
         endif
#endif
      enddo
   enddo
 else
!-------------------------------
! For positive definite tracers:
!-------------------------------
! jord=11: PPM mono constraint (Lin 2004)
! jord=12: Huynh 2nd constraint (Lin 2004) + positive definite (Lin & Rood 1996)
! jord>12: positive definite only (Lin & Rood 1996)

   js3 = max(3,js-1); je3 = min(npy-3,je+1)

   do j=js-2,je+2
      do i=ifirst,ilast
         xt = 0.25*(q(i,j+1) - q(i,j-1))
         dm(i,j) = sign(min(abs(xt), max(q(i,j-1), q(i,j), q(i,j+1)) - q(i,j),   &
                            q(i,j) - min(q(i,j-1), q(i,j), q(i,j+1))), xt)
      enddo
   enddo

   if (grid_type < 3) then

      do j=js3,min(npy-2,je+2)
         do i=ifirst,ilast
            al(i,j) = 0.5*(q(i,j-1)+q(i,j)) + r3*(dm(i,j-1) - dm(i,j))
         enddo
      enddo

      if ( jord==11 ) then
         do j=js3,je3
            do i=ifirst,ilast
               xt = 2.*dm(i,j)
               bl(i,j) = -sign(min(abs(xt), abs(al(i,j  )-q(i,j))), xt)
               br(i,j) =  sign(min(abs(xt), abs(al(i,j+1)-q(i,j))), xt)
            enddo
         enddo
      elseif( jord==12 ) then
         do j=js-3,je+2
            do i=ifirst,ilast
               dq(i,j) = q(i,j+1) - q(i,j)
            enddo
         enddo
         do j=js3,je3
            do i=ifirst,ilast
               pmp = -2.*dq(i,j) 
               lac = pmp + 1.5*dq(i,j+1)
               bl(i,j) = min(max(0.,pmp,lac), max(al(i,j  )-q(i,j), min(0.,pmp,lac)))
               pmp = 2.*dq(i,j-1)
               lac = pmp - 1.5*dq(i,j-2)
               br(i,j) = min(max(0.,pmp,lac), max(al(i,j+1)-q(i,j), min(0.,pmp,lac)))
            enddo
         enddo
      else
         do j=js3,je3
            do i=ifirst,ilast
               bl(i,j) = al(i,j  ) - q(i,j)
               br(i,j) = al(i,j+1) - q(i,j)
            enddo
         enddo
      endif
      
      if ( jord/=11 ) then
! Positive definite constraint:
         do j=js3,je3
            call pert_ppm(ilast-ifirst+1, q(ifirst,j), bl(ifirst,j), br(ifirst,j), 0)
         enddo
      endif

!--------------
! Fix the edges:
!--------------
      if( js==1 ) then
         do i=ifirst,ilast
            br(i,2) = al(i,3) - q(i,2)
!           xt = t11*(q(i,0)+q(i,1)) + t12*(q(i,-1)+q(i,2))   &
!              + t13*(dm(i,2)-dm(i,-1))
!!!         xt = 0.75*(q(i,0)+q(i,1)) - 0.25*(q(i,-1)+q(i,2))
            xt = 0.5*((2.*dya(i,1)+dya(i,2))*(q(i,0)+q(i,1))  &
               -dya(i,1)*(q(i,-1)+q(i,2))) / (dya(i,1)+dya(i,2))
#ifdef EDGE_MONO
            if ( (i>-2 .and. i<1) .or. (i>(npx-3) .and. i<npx) ) then
               xt = min( xt, max(q(i,0), q(i,1), q(i-1,0), q(i-1,1)) )
               xt = max( xt, min(q(i,0), q(i,1), q(i-1,0), q(i-1,1)) )
            elseif ( (i<3 .and.i>0) .or. (i>(npx-1).and.i<(npx+2)) ) then
               xt = min( xt, max(q(i,0), q(i,1), q(i+1,0), q(i+1,1)) )
               xt = max( xt, min(q(i,0), q(i,1), q(i+1,0), q(i+1,1)) )
            else
               xt = max(0., xt)
            endif
#else
            xt = max(0., xt)
#endif
            bl(i,1) = xt - q(i,1)
            br(i,0) = xt - q(i,0)
            xt = 4./7.*dm(i,-1) + 11./14.*q(i,-1) + 3./14.*q(i,0)
#ifdef S_MONO
            xt = min( xt, max(q(i,-1), q(i,0)) )
            xt = max( xt, min(q(i,-1), q(i,0)) )
#else
            xt = max(0., xt)
#endif
            bl(i,0) = xt - q(i,0)

            xt = 3./14.*q(i,1) + 11./14.*q(i,2) - 4./7.*dm(i,2)
#ifdef S_MONO
            xt = min( xt, max(q(i,1), q(i,2)) )
            xt = max( xt, min(q(i,1), q(i,2)) )
#else
            xt = max(0., xt)
#endif
            br(i,1) = xt - q(i,1)
            bl(i,2) = xt - q(i,2)
         enddo
         do j=0,2
            call pert_ppm(ilast-ifirst+1, q(ifirst,j), bl(ifirst,j), br(ifirst,j), 1)
         enddo
      endif

      if( (je+1)==npy ) then
         do i=ifirst,ilast
            bl(i,npy-2) = al(i,npy-2) - q(i,npy-2)
!           xt = t11*(q(i,npy-1)+q(i,npy)) + t12*(q(i,npy-2)+q(i,npy+1))   &
!               + t13*(dm(i,npy+1)-dm(i,npy-2))
!!!         xt = 0.75*(q(i,npy-1)+q(i,npy)) - 0.25*(q(i,npy-2)+q(i,npy+1))
            xt = 0.5*((2.*dya(i,npy-1)+dya(i,npy-2))*(q(i,npy-1)+q(i,npy)) &
               - dya(i,npy-1)*(q(i,npy-2)+q(i,npy+1)))  &
                / ( dya(i,npy-1)+dya(i,npy-2) )
#ifdef EDGE_MONO
            if ( (i>-2 .and. i<1) .or. (i>(npx-3) .and. i<npx) ) then
               xt = min( xt, max(q(i,npy-1), q(i,npy), q(i-1,npy-1), q(i-1,npy)) )
               xt = max( xt, min(q(i,npy-1), q(i,npy), q(i-1,npy-1), q(i-1,npy)) )
            elseif ( (i<3 .and.i>0) .or. (i>(npx-1).and.i<(npx+2)) ) then
               xt = min( xt, max(q(i,npy-1), q(i,npy), q(i+1,npy-1), q(i+1,npy)) )
               xt = max( xt, min(q(i,npy-1), q(i,npy), q(i+1,npy-1), q(i+1,npy)) )
            else
               xt = max(0., xt)
            endif
#else
            xt = max(0., xt)
#endif
            br(i,npy-1) = xt - q(i,npy-1)
            bl(i,npy  ) = xt - q(i,npy)
            xt = 3./14.*q(i,npy) + 11./14.*q(i,npy+1) - 4./7.*dm(i,npy+1)
#ifdef S_MONO
            xt = min( xt, max(q(i,npy), q(i,npy+1)) )
            xt = max( xt, min(q(i,npy), q(i,npy+1)) )
#else
            xt = max(0., xt)
#endif
            br(i,npy) = xt - q(i,npy)
            xt = 3./14.*q(i,npy-1) + 11./14.*q(i,npy-2) + 4./7.*dm(i,npy-2)
#ifdef S_MONO
            xt = min( xt, max(q(i,npy-2), q(i,npy-1)) )
            xt = max( xt, min(q(i,npy-2), q(i,npy-1)) )
#else
            xt = max(0., xt)
#endif
            br(i,npy-2) = xt - q(i,npy-2)
            bl(i,npy-1) = xt - q(i,npy-1)
         enddo
         do j=npy-2,npy
            call pert_ppm(ilast-ifirst+1, q(ifirst,j), bl(ifirst,j), br(ifirst,j), 1)
         enddo
      endif
   else

      do j=js-1,je+2
         do i=ifirst,ilast
            al(i,j) = 0.5*(q(i,j-1)+q(i,j)) + r3*(dm(i,j-1) - dm(i,j))
         enddo
      enddo

      if ( jord==11 ) then
         do j=js-1,je+1
            do i=ifirst,ilast
               xt = 2.*dm(i,j)
               bl(i,j) = -sign(min(abs(xt), abs(al(i,j  )-q(i,j))), xt)
               br(i,j) =  sign(min(abs(xt), abs(al(i,j+1)-q(i,j))), xt)
            enddo
         enddo
      elseif( jord==12 ) then
         do j=js-3,je+2
            do i=ifirst,ilast
               dq(i,j) = q(i,j+1) - q(i,j)
            enddo
         enddo
         do j=js-1,je+1
            do i=ifirst,ilast
               pmp = -2.*dq(i,j) 
               lac = pmp + 1.5*dq(i,j+1)
               bl(i,j) = min(max(0.,pmp,lac), max(al(i,j  )-q(i,j), min(0.,pmp,lac)))
               pmp = 2.*dq(i,j-1)
               lac = pmp - 1.5*dq(i,j-2)
               br(i,j) = min(max(0.,pmp,lac), max(al(i,j+1)-q(i,j), min(0.,pmp,lac)))
            enddo
         enddo
      else
         do j=js-1,je+1
            do i=ifirst,ilast
               bl(i,j) = al(i,j  ) - q(i,j)
               br(i,j) = al(i,j+1) - q(i,j)
            enddo
         enddo
      endif

      if ( jord/=11 ) then
! Positive definite constraint:
         do j=js-1,je+1
            call pert_ppm(ilast-ifirst+1, q(ifirst,j), bl(ifirst,j), br(ifirst,j), 0)
         enddo
      endif

   endif

   do j=js,je+1
      do i=ifirst,ilast
         if( c(i,j)>0. ) then
            flux(i,j) = q(i,j-1) + (1.-c(i,j))*(br(i,j-1)-c(i,j)*(bl(i,j-1)+br(i,j-1)))
         else
            flux(i,j) = q(i,j  ) + (1.+c(i,j))*(bl(i,j  )+c(i,j)*(bl(i,j  )+br(i,j  )))
         endif
      enddo
   enddo
 endif

 end subroutine fyppm

                       
 subroutine xmist_2d(q, dm, ng, iord, ifirst, ilast, jfirst, jlast)
 integer, INTENT(IN) :: ifirst, ilast             !  X-Dir strip
 integer, INTENT(IN) :: jfirst, jlast             !  Y-Dir strip
 integer, INTENT(IN) :: ng                        !  Ghost Points 
 integer, INTENT(IN) :: iord
 real   , INTENT(IN) :: q(ifirst-ng:ilast+ng,jfirst:jlast)
 real   , INTENT(INOUT) :: dm(ifirst-ng:ilast+ng,jfirst:jlast)
!---------------------------------------------------------------------
! Local variables
   integer i, j
   real   qmax, qmin, tmp

      do j=jfirst,jlast
        do i=ifirst-ng+1,ilast+ng-1
           dm(i,j) = 0.25*(q(i+1,j) - q(i-1,j))
        enddo
      enddo

   if( iord > 0 ) then
!
! Applies monotonic slope constraint (off if iord less than zero)
!
        do j=jfirst,jlast
          do i=ifirst-ng+1,ilast+ng-1
            qmax = max(q(i-1,j),q(i,j),q(i+1,j)) - q(i,j)
            qmin = q(i,j) - min(q(i-1,j),q(i,j),q(i+1,j))
            dm(i,j) = sign(min(abs(dm(i,j)),qmin,qmax),dm(i,j))
          enddo
        enddo
   endif
 end subroutine xmist_2d



 subroutine ymist(q, dm, ng, jord, ifirst, ilast, jfirst, jlast)
! !INPUT PARAMETERS:
 integer, INTENT(IN) :: ifirst, ilast             !  X-Dir strip
 integer, INTENT(IN) :: jfirst, jlast             !  Y-Dir strip
 integer, INTENT(IN) :: ng                        !  Ghost Points 
 integer, INTENT(IN) :: jord
 real   , INTENT(IN) :: q(ifirst:ilast,jfirst-ng:jlast+ng)
! !OUTPUT PARAMETERS:
 real   , INTENT(INOUT) :: dm(ifirst:ilast,jfirst-ng:jlast+ng)
! Local variables

   integer i, j
   real   qmax, qmin, tmp

      do j=jfirst-ng+1,jlast+ng-1
        do i=ifirst,ilast
           dm(i,j) = 0.25*(q(i,j+1) - q(i,j-1))
        enddo
      enddo

   if( jord > 0 ) then
!
! Applies monotonic slope constraint (off if jord less than zero)
!
        do j=jfirst-ng+1,jlast+ng-1
          do i=ifirst,ilast
            qmax = max(q(i,j-1),q(i,j),q(i,j+1)) - q(i,j)
            qmin = q(i,j) - min(q(i,j-1),q(i,j),q(i,j+1))
            dm(i,j) = sign(min(abs(dm(i,j)),qmin,qmax),dm(i,j))
          enddo
        enddo
   endif
!EOC
 end subroutine ymist


 subroutine mp_ghost_ew(im, jm, km, nq, ifirst, ilast, jfirst, jlast, &
                              kfirst, klast, ng_w, ng_e, ng_s, ng_n, q_ghst, q)
!
! !INPUT PARAMETERS:
      integer, intent(in):: im, jm, km, nq
      integer, intent(in):: ifirst, ilast
      integer, intent(in):: jfirst, jlast
      integer, intent(in):: kfirst, klast
      integer, intent(in):: ng_e      ! eastern  zones to ghost
      integer, intent(in):: ng_w      ! western  zones to ghost
      integer, intent(in):: ng_s      ! southern zones to ghost
      integer, intent(in):: ng_n      ! northern zones to ghost
      real, intent(inout):: q_ghst(ifirst-ng_w:ilast+ng_e,jfirst-ng_s:jlast+ng_n,kfirst:klast,nq)
      real, optional, intent(in):: q(ifirst:ilast,jfirst:jlast,kfirst:klast,nq)
!
! !DESCRIPTION:
!
!     Ghost 4d east/west 
!
! !REVISION HISTORY:
!    2005.08.22   Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: i,j,k,n

      if (present(q)) then
         q_ghst(ifirst:ilast,jfirst:jlast,kfirst:klast,1:nq) = &
              q(ifirst:ilast,jfirst:jlast,kfirst:klast,1:nq)
      endif

!      Assume Periodicity in X-dir and not overlapping
      do n=1,nq
         do k=kfirst,klast
            do j=jfirst-ng_s,jlast+ng_n
               do i=1, ng_w
                  q_ghst(ifirst-i,j,k,n) = q_ghst(ilast-i+1,j,k,n)
               enddo
               do i=1, ng_e
                  q_ghst(ilast+i,j,k,n) = q_ghst(ifirst+i-1,j,k,n)
               enddo
            enddo
         enddo
      enddo

 end subroutine mp_ghost_ew



 subroutine pert_ppm(im, a0, al, ar, iv)
 integer, intent(in):: im
 integer, intent(in):: iv
 real, intent(in)   :: a0(im)
 real, intent(inout):: al(im), ar(im)
! Local:
 real a4, da1, da2, a6da, fmin
 integer i
 real, parameter:: r12 = 1./12.

!-----------------------------------
! Optimized PPM in perturbation form:
!-----------------------------------

 if ( iv==0 ) then
! Positive definite constraint
    do i=1,im
        a4 = -3.*(ar(i) + al(i))
       da1 =      ar(i) - al(i)
      if( abs(da1) < -a4 ) then
         fmin = a0(i) + 0.25/a4*da1**2 + a4*r12
         if( fmin < 0. ) then
             if( ar(i)>0. .and. al(i)>0. ) then
                 ar(i) = 0.
                 al(i) = 0.
             elseif( da1 > 0. ) then
                 ar(i) = -2.*al(i)
             else
                 al(i) = -2.*ar(i)
             endif
         endif
      endif
    enddo
 else
! Standard PPM constraint
    do i=1,im
       if ( al(i)*ar(i) < 0. ) then
            da1 = al(i) - ar(i)
            da2 = da1**2
            a6da = 3.*(al(i)+ar(i))*da1
            if( a6da < -da2 ) then
                ar(i) = -2.*al(i)
            elseif( a6da > da2 ) then
                al(i) = -2.*ar(i)
            endif
       else
! effect of dm=0 included here
            al(i) = 0.
            ar(i) = 0.
       endif
  enddo
 endif

 end subroutine pert_ppm

end module tpcore
