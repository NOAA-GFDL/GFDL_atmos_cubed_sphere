module tpcore
!BOP
!
! !MODULE: tp_core --- A collection of routines to support FV transport
!
 use mp_mod,     only: is,js,ie,je, ng, isd,jsd,ied,jed
 use grid_utils, only: sw_corner, se_corner, ne_corner, nw_corner,   &
                       cx1, cx2, cy1, cy2

 implicit none

 private
 public fv_tp_2d, std_ppm

 real, parameter:: huge = 1.e+20
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

 integer:: iad=2

!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------

CONTAINS

 subroutine fv_tp_2d(q, crx, cry, im, jm, hord, fx, fy, &
                     xfx, yfx, area, uniform_ppm, mfx, mfy)
!-----------------------------------------------------------------------
! !INPUT PARAMETERS:
   integer, intent(in)::im, jm         ! Dimensions
   integer, intent(in)::hord

   real, intent(in):: crx(is:ie+1,jsd:jed)  !
   real, intent(in):: xfx(is:ie+1,jsd:jed)  !

   real, intent(in):: cry(isd:ied,js:je+1 )  !
   real, intent(in):: yfx(isd:ied,js:je+1 )  !
   logical, intent(IN) :: uniform_ppm
   real, intent(in):: area(isd:ied,jsd:jed)
   real, intent(inout):: q(isd:ied,jsd:jed)  ! transported scalar

! !OUTPUT PARAMETERS:
   real, intent(out)::fx(is:ie+1 ,js:je)    ! Flux in x ( E )
   real, intent(out)::fy(is:ie,   js:je+1 )    ! Flux in y ( N )

! optional Arguments:
   real, OPTIONAL, intent(in):: mfx(is:ie+1,js:je  )  ! Mass Flux X-Dir
   real, OPTIONAL, intent(in):: mfy(is:ie  ,js:je+1)  ! Mass Flux Y-Dir

! Local:
   real q2d_i(isd:ied,js:je)
   real q2d_j(is:ie,jsd:jed)
   real   fx2(is:ie+1,jsd:jed)
   real   fy2(isd:ied,js:je+1)
   integer i, j, npx, npy
   integer ord

   npx = im+1
   npy = jm+1

   if ( present(mfx) .and. present(mfy) ) then
!--------------------------------------------------
! Fully symmetrical revised Lin-Rood 1996 algorithm
!--------------------------------------------------
#ifdef QUICK_INNER
          if ( hord==6 ) then
               ord = -2       ! un-constrained linear scheme
          else
               ord = min(iad,hord)
          endif
#else
          ord = hord
#endif

! E-W flux
          call copy_corners(q, npx, npy, 2)
          call ytp(fy2, q, cry, yfx, ng, ord, isd, ied, js, je, npx, npy, uniform_ppm)

          do j=js,je
              do i=isd,ied
                 q2d_i(i,j) = 0.5 * (    q(i,j) +                    &
                             ( q(i,j)*area(i,j) + fy2(i,j)-fy2(i,j+1) ) /  &
                             (        area(i,j) + yfx(i,j)-yfx(i,j+1) ) )
              enddo
          enddo

          call xtp(ng, fx, q2d_i(isd,js), crx(is,js), hord, mfx,   &
                   is, ie, js,je, npx, npy, uniform_ppm)
! N-S flux

          call copy_corners(q, npx, npy, 1)
          call xtp(ng, fx2, q, crx, ord, xfx, is, ie, jsd,jed, npx, npy, uniform_ppm)

          do j=jsd,jed
             do i=is,ie
                q2d_j(i,j) = 0.5 * (    q(i,j) +                    &
                            ( q(i,j)*area(i,j) + fx2(i,j)-fx2(i+1,j) ) /  &
                                   ( area(i,j) + xfx(i,j)-xfx(i+1,j) ) )
             enddo
          enddo

          call ytp(fy, q2d_j, cry(is:ie,js:je+1),  &
                   mfy, ng, hord, is, ie, js, je, npx, npy, uniform_ppm)
   else
!------------------------------------------------
! Directionally symmetrical Flux Averaging Method
! Fully monotonic transport with ORD=1,2,3,4
!------------------------------------------------
! This is only good for (delp, vort) transport because mass 
! fluxes would need to be ghosted otherwise.
!
!----------
! 1st sweep
!----------
          call copy_corners(q, npx, npy, 2)
          call ytp(fy2, q, cry, yfx, ng, hord, isd, ied, js, je, npx, npy, uniform_ppm)

          do j=js,je
              do i=isd,ied
                 q2d_i(i,j) = ( q(i,j)*area(i,j) + fy2(i,j)-fy2(i,j+1) ) /  &
                              (        area(i,j) + yfx(i,j)-yfx(i,j+1) )  
              enddo
          enddo

          call xtp(ng, fx, q2d_i(isd,js), crx(is,js), hord, &
                   xfx(is,js), is, ie, js,je, npx, npy, uniform_ppm)
!----------
! 2nd sweep
!----------
          call copy_corners(q, npx, npy, 1)
          call xtp(ng, fx2, q, crx, hord, xfx, is, ie, jsd,jed, npx, npy, uniform_ppm)

          do j=jsd,jed
             do i=is,ie
                q2d_j(i,j) = ( q(i,j)*area(i,j) + fx2(i,j)-fx2(i+1,j) ) /  &
                             (        area(i,j) + xfx(i,j)-xfx(i+1,j) )
             enddo
          enddo

          call ytp(fy, q2d_j, cry(is:ie,js:je+1),   &
                   yfx(is:ie,js:je+1), ng, hord, is, ie, js, je, npx, npy, uniform_ppm)

!----------------
! Flux averaging:
!----------------
! SJL notes:
! By averaging two anti-symmetrical "monotonic" schemes
! we obtain a (twice expensive) symmetrical monotonic scheme
! The averaging is necessary because of the Cubed Sphere layout
          do j=js,je
             do i=is,ie+1
                fx(i,j) = 0.5*(fx(i,j) + fx2(i,j))  
             enddo
          enddo
          do j=js,je+1
             do i=is,ie
                fy(i,j) = 0.5*(fy(i,j) + fy2(i,j))  
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


 subroutine xtp(ng, fx,  q,  c, iord,  xfx, &
                ifirst, ilast, jfirst, jlast, npx, npy, uniform_ppm)
!---------------------------------------------------------------

! !INPUT PARAMETERS:
   integer, intent(IN):: ifirst, ilast   !  X-Dir strip
   integer, intent(IN):: jfirst, jlast   !  Y-Dir strip
   integer, intent(IN):: npx, npy
   integer, intent(IN):: ng              !  Ghost dependencies
   integer, intent(in):: iord
   real   , intent(in)::   c(ifirst   :ilast+1 ,jfirst:jlast)      ! Courant numbers
   real   , intent(in)::   q(ifirst-ng:ilast+ng,jfirst:jlast)
   real   , intent(in):: xfx(ifirst   :ilast+1 ,jfirst:jlast)  
   logical, intent(IN) :: uniform_ppm
! !OUTPUT PARAMETERS: 
   real, intent(out):: fx(ifirst:ilast+1,jfirst:jlast)           
! Local:
   integer iu
   integer i, j
   real   dm(ifirst-ng:ilast+ng,jfirst:jlast)

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
      call fxppm(c, q, fx, ng, iord, ifirst, ilast, jfirst, jlast, npx, npy, uniform_ppm)
   endif

   do j=jfirst,jlast
      do i=ifirst,ilast+1
         fx(i,j) = fx(i,j)*xfx(i,j)
      enddo
   enddo

 end subroutine xtp



 subroutine fxppm(c, q, flux, ng, iord, ifirst, ilast, jfirst,   &
                  jlast,npx,npy, uniform_ppm)
!-----------------------------------------------------------------------

! !USES:

! !INPUT PARAMETERS:
 integer, INTENT(IN)    :: ifirst, ilast               !  X-Dir strip
 integer, INTENT(IN)    :: jfirst, jlast               !  Y-Dir strip
 integer, INTENT(IN)    :: ng                          !  Ghost Points
 integer, INTENT(IN)    :: iord                        !  Approximation order
 integer, INTENT(IN)    :: npx, npy
 logical, intent(IN) :: uniform_ppm
 real   , INTENT(IN)    ::  q(ifirst-ng:ilast+ng,jfirst:jlast) !  mean value needed only N*2 S*2
 real   , INTENT(IN)    ::  c(ifirst   :ilast+1 ,jfirst:jlast) !  Courant   N (like FLUX)
 
! !OUTPUT PARAMETERS:
 real   , INTENT(OUT) :: flux(ifirst:ilast+1,jfirst:jlast)    !  Flux
!---------------------------------------------------------------------
! Local
 real al(ifirst-1:ilast+2), dm1(ifirst-2:ilast+2)
 real bl(ifirst-1:ilast+1)
 real br(ifirst-1:ilast+1)
 real dq(ifirst-3:ilast+2)
 real dl, dr, xt, pmp, lac, dqt, q0
 integer i, j

 if (iord<=4) then
     do j=jfirst,jlast

           do i=ifirst-3,ilast+2
              dq(i) = q(i+1,j) - q(i,j)
           enddo

           if ( uniform_ppm ) then
             do i=ifirst-2,ilast+2
                xt = 0.25*(dq(i-1) + dq(i))
                dm1(i) = sign(min(abs(xt), abs(dq(i-1)), abs(dq(i))), xt)
             enddo
           else
! Variable grid size:
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
                xt = 0.25*(dq(i-1) + dq(i))
                dm1(i) = sign(min(abs(xt), abs(dq(i-1)), abs(dq(i))), xt)
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
 elseif( iord==8 ) then
! iord==8sed on 4
     do j=jfirst,jlast
        do i=ifirst-2,ilast+2
               xt = 0.25*(q(i+1,j) - q(i-1,j))
           dm1(i) = sign(min(abs(xt), max(q(i-1,j), q(i,j), q(i+1,j)) - q(i,j),  &
                             q(i,j) - min(q(i-1,j), q(i,j), q(i+1,j))), xt)
        enddo

        do i=max(3,ifirst-1),min(npx-2,ilast+2)
           al(i) = 0.5*(q(i-1,j)+q(i,j)) + r3*(dm1(i-1) - dm1(i))
        enddo

        do i=max(3,ifirst-1),min(npx-3,ilast+1)
              xt = 2.*dm1(i)
           bl(i) = sign(min(abs(xt), abs(al(i)  -q(i,j))), xt)
           br(i) = sign(min(abs(xt), abs(al(i+1)-q(i,j))), xt)
        enddo

!--------------
! fix the edges
!--------------
             if ( is==1 ) then
                 q0 = 3./7.*(dm1(2)-dm1(-1)) + 27./28.*(q(0,j)+q(1,j)) -  &
                                               13./28.*(q(-1,j)+q(2,j))
! Inside:
                   dqt = q(1,j) - q(2,j)
!                  bl(1) =  6./7.*dm1(2) + 13./14.*dqt
                   bl(1) =  q0 - q(1,j)
                   br(1) = -4./7.*dm1(2) - 11./14.*dqt
                      xt = 0.5*(br(1) - bl(1))
                  dm1(1) = sign(min(abs(xt), abs(dqt)), xt)
                     xt = 2.*dm1(1)
                  bl(1) = sign( min(abs(xt), abs(bl(1))), xt )
                  br(1) = sign( min(abs(xt), abs(br(1))), xt )
!
                  bl(2) = 0.5*dqt + r3*(dm1(1) - dm1(2))
                     xt = 2.*dm1(2)
                  bl(2) = sign(min(abs(xt), abs(bl(2))), xt)
                  br(2) = sign(min(abs(xt), abs(al(3)-q(2,j))), xt)
! Outside: (mirrors what's done inside)
                   dqt = q(0,j) - q(-1,j)
!                  br(0) = -6./7.*dm1(-1) + 13./14.*dqt   ! face edge
                   br(0) = q0 - q(0,j)
                   bl(0) =  4./7.*dm1(-1) - 11./14.*dqt
                      xt = 0.5*(br(0) - bl(0))
                  dm1(0) = sign(min(abs(xt), abs(dqt)), xt)
                     xt = 2.*dm1(0)
                  bl(0) = sign( min(abs(xt), abs(bl(0))), xt )
                  br(0) = sign( min(abs(xt), abs(br(0))), xt )
             endif
             if ( (ie+1)==npx ) then
                 q0 = 3./7.*(dm1(npx+1)-dm1(npx-2)) + 27./28.*(q(npx-1,j)+q(npx,  j)) - &
                                                      13./28.*(q(npx-2,j)+q(npx+1,j))
! Inside:
                   dqt = q(npx-1,j) - q(npx-2,j)
!                  br(npx-1) = -6./7.*dm1(npx-2) + 13./14.*dqt
                   br(npx-1) = q0 - q(npx-1,j)
                   bl(npx-1) =  4./7.*dm1(npx-2) - 11./14.*dqt
                      xt = 0.5*(br(npx-1) - bl(npx-1))
                  dm1(npx-1) = sign(min(abs(xt), abs(dqt)), xt)
                     xt = 2.*dm1(npx-1)
                  bl(npx-1) = sign( min(abs(xt), abs(bl(npx-1))), xt )
                  br(npx-1) = sign( min(abs(xt), abs(br(npx-1))), xt )
!
                  br(npx-2) = 0.5*dqt + r3*(dm1(npx-2) - dm1(npx-1))
                         xt = 2.*dm1(npx-2)
                  br(npx-2) = sign(min(abs(xt), abs(br(npx-2))), xt)
                  bl(npx-2) = sign(min(abs(xt), abs(al(npx-2)-q(npx-2,j))), xt)
! Outside: (mirrors what's done inside)
                   dqt = q(npx,j) - q(npx+1,j)
!                  bl(npx) =  6./7.*dm1(npx+1) + 13./14.*dqt
                   bl(npx) =  q0 - q(npx,j)
                   br(npx) = -4./7.*dm1(npx+1) - 11./14.*dqt
                        xt = 0.5*(br(npx) - bl(npx))
                  dm1(npx) = sign(min(abs(xt), abs(dqt)), xt)
                     xt = 2.*dm1(npx)
                  bl(npx) = sign( min(abs(xt), abs(bl(npx))), xt )
                  br(npx) = sign( min(abs(xt), abs(br(npx))), xt )
        endif

        do i=ifirst,ilast+1
           if( c(i,j)>0. ) then
              flux(i,j) = q(i-1,j) + (1.-c(i,j))*(br(i-1)+c(i,j)*(bl(i-1)-br(i-1)))
           else
              flux(i,j) = q(i,  j) - (1.+c(i,j))*(bl(i  )+c(i,j)*(bl(i  )-br(i  )))
           endif
        enddo
     enddo
 elseif( iord==9 ) then
! Based on 5; with absolute upwind algorithm at edges
     do j=jfirst,jlast
        do i=is-3,ie+2
           dq(i) = q(i+1,j) - q(i,j)
        enddo

        if ( uniform_ppm ) then
             do i=is-2,ie+2
                xt = 0.25*(dq(i-1) + dq(i))
                dm1(i) = sign(min(abs(xt), abs(dq(i-1)), abs(dq(i))), xt)
             enddo
        else
             do i=is-2,ie+2
                xt = cx1(i,j)*dq(i-1) + cx2(i,j)*dq(i)
                dm1(i) = sign(min(abs(xt), abs(dq(i-1)), abs(dq(i))), xt)
             enddo
        endif

        do i=max(3,is-1),min(npx-2,ie+2)
           al(i) = 0.5*(q(i-1,j)+q(i,j)) + r3*(dm1(i-1)-dm1(i))
        enddo

        do i=max(3,is-1),min(npx-3,ie+1)
           pmp = -2.*dq(i)
           lac = pmp + 1.5*dq(i+1)
           bl(i) = min(max(0., pmp, lac), max(al(i)-q(i,j), min(0.,pmp, lac)))
           pmp = 2.*dq(i-1)
           lac = pmp - 1.5*dq(i-2)
           br(i) = min(max(0., pmp, lac), max(al(i+1)-q(i,j), min(0.,pmp, lac)))
        enddo

!--------------
! fix the edges
!--------------
        if ( is==1 ) then
             br(2) = al(3) - q(2,j)
                xt = 3./14.*q(1,j) + 11./14.*q( 2,j) - 4./7.*dm1( 2)
             bl(2) = xt - q(2,j)
             br(1) = xt - q(1,j)
             bl(0) = 11./14.*(q(-1,j) - q(0,j)) + 4./7.*dm1(-1)
             xt = 27./28.*(q(0,j)+q(1,j)) - 13./28.*(q(-1,j)+q(2,j))   &
                  + 3./7.*(dm1(2)-dm1(-1))
             bl(1) = xt - q(1,j)
             br(0) = xt - q(0,j)
           call std_ppm(3, q(0,j), bl(0), br(0))
        endif
        if ( (ie+1)==npx ) then
             bl(npx-2) = al(npx-2) - q(npx-2,j)
                    xt = 3./14.*q(npx-1,j) + 11./14.*q(npx-2,j) + 4./7.*dm1(npx-2)
             br(npx-2) = xt - q(npx-2,j)
             bl(npx-1) = xt - q(npx-1,j)
             br(npx) = 11./14.*(q(npx+1,j) - q(npx,j)) - 4./7.*dm1(npx+1)
             xt = 27./28.*(q(npx-1,j)+q(npx,j)) - 13./28.*(q(npx-2,j)+q(npx+1,j))   &
                  + 3./7.*(dm1(npx+1)-dm1(npx-2))
             br(npx-1) = xt - q(npx-1,j)
             bl(npx  ) = xt - q(npx  ,j)
           call std_ppm(3, q(npx-2,j), bl(npx-2), br(npx-2))
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



 subroutine fyppm(c,  q,  dm, flux, ng, jord, ifirst, ilast, jfirst, jlast,    &
                  npx, npy, uniform_ppm)
!-----------------------------------------------------------------------

! !USES:

! !INPUT PARAMETERS:
 integer, INTENT(IN) :: ifirst, ilast               !  X-Dir strip
 integer, INTENT(IN) :: jfirst, jlast               !  Y-Dir strip
 integer, INTENT(IN) :: ng                          !  Ghost Points
 integer, INTENT(IN) :: jord                        !  Approximation order
 integer, INTENT(IN) :: npx, npy
 logical, intent(IN) :: uniform_ppm
 real   , INTENT(IN) ::  q(ifirst:ilast,jfirst-ng:jlast+ng) !  mean value needed only N*2 S*2
 real   , INTENT(INOUT) :: dm(ifirst:ilast,jfirst-ng:jlast+ng) !  Slope     needed only N*2 S*2
 real   , INTENT(IN) ::  c(ifirst:ilast,jfirst   :jlast+1 ) !  Courant   N (like FLUX)

! !OUTPUT PARAMETERS:
 real   , INTENT(OUT) :: flux(ifirst:ilast,jfirst:jlast+1)    !  Flux
! Local:
 real al(ifirst:ilast,jfirst-1:jlast+2)
 real bl(ifirst:ilast,jfirst-1:jlast+1)
 real br(ifirst:ilast,jfirst-1:jlast+1)
 real dq(ifirst:ilast,jfirst-3:jlast+2)
 real xt, dl, dr, pmp, lac, dqt, q0
 integer i, j

 if (jord<=4) then

! Variable grid size:
       do j=jfirst-3,jlast+2
          do i=ifirst,ilast
             dq(i,j) = q(i,j+1) - q(i,j)
          enddo
       enddo

       if ( uniform_ppm ) then
           do j=jfirst-2,jlast+2
              do i=ifirst,ilast
                 xt = 0.25*(dq(i,j-1) + dq(i,j))
                 dm(i,j) = sign(min(abs(xt), abs(dq(i,j-1)), abs(dq(i,j))), xt)
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
              xt = 0.25*(dq(i,j-1) + dq(i,j))
              dm(i,j) = sign(min(abs(xt), abs(dq(i,j-1)), abs(dq(i,j))), xt)
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

 elseif( jord==8 ) then
! jord=8; based on 4
   do j=jfirst-2,jlast+2
      do i=ifirst,ilast
              xt = 0.25*(q(i,j+1) - q(i,j-1))
         dm(i,j) = sign(min(abs(xt), max(q(i,j-1), q(i,j), q(i,j+1)) - q(i,j),   &
                            q(i,j) - min(q(i,j-1), q(i,j), q(i,j+1))), xt)
      enddo
   enddo

   do j=max(3,jfirst-1),min(npy-2,jlast+2)
      do i=ifirst,ilast
         al(i,j) = 0.5*(q(i,j-1)+q(i,j)) + r3*(dm(i,j-1)-dm(i,j))
      enddo
   enddo

   do j=max(3,jfirst-1),min(npy-3,jlast+1)
      do i=ifirst,ilast
              xt = 2.*dm(i,j)
         bl(i,j) = sign(min(abs(xt), abs(al(i,j  )-q(i,j))), xt)
         br(i,j) = sign(min(abs(xt), abs(al(i,j+1)-q(i,j))), xt)
      enddo
   enddo

      if( js==1 ) then
        do i=ifirst,ilast
!       if( i==1 .or. i==npx ) then
          q0 = 3./7.*(dm(i,2)-dm(i,-1)) + 27./28.*(q(i, 0)+q(i,1)) -  &
                                          13./28.*(q(i,-1)+q(i,2))
! Inside:
          dqt = q(i,1) - q(i,2)
!         bl(i,1) =  6./7.*dm(i,2) + 13./14.*dqt
          bl(i,1) =  q0 - q(i,1)
          br(i,1) = -4./7.*dm(i,2) - 11./14.*dqt
               xt = 0.5*(br(i,1) - bl(i,1))
          dm(i,1) = sign(min(abs(xt), abs(dqt)), xt)
               xt = 2.*dm(i,1)
          bl(i,1) = sign( min(abs(xt), abs(bl(i,1))), xt )
          br(i,1) = sign( min(abs(xt), abs(br(i,1))), xt )
!
          bl(i,2) = 0.5*dqt + r3*(dm(i,1) - dm(i,2))
               xt = 2.*dm(i,2)
          bl(i,2) = sign(min(abs(xt), abs(bl(i,2))), xt)
          br(i,2) = sign(min(abs(xt), abs(al(i,3)-q(i,2))), xt)
! Outside: (mirrors what's done inside)
              dqt = q(i,0) - q(i,-1)
!         br(i,0) = -6./7.*dm(i,-1) + 13./14.*dqt
          br(i,0) = q0 - q(i,0)
          bl(i,0) =  4./7.*dm(i,-1) - 11./14.*dqt
               xt = 0.5*(br(i,0) - bl(i,0))
          dm(i,0) = sign(min(abs(xt), abs(dqt)), xt)
               xt = 2.*dm(i,0)
          bl(i,0) = sign( min(abs(xt), abs(bl(i,0))), xt )
          br(i,0) = sign( min(abs(xt), abs(br(i,0))), xt )
!       endif 
        enddo
      endif 
      if( (je+1)==npy ) then
        do i=ifirst,ilast
!       if( i==1 .or. i==npx ) then
           q0 = 3./7.*(dm(i,npy+1)-dm(i,npy-2)) + 27./28.*(q(i,npy-1)+q(i,npy  )) -  &
                                                  13./28.*(q(i,npy-2)+q(i,npy+1))
! Inside:
                   dqt = q(i,npy-1) - q(i,npy-2)
!          br(i,npy-1) = -6./7.*dm(i,npy-2) + 13./14.*dqt
           br(i,npy-1) = q0 - q(i,npy-1)
           bl(i,npy-1) =  4./7.*dm(i,npy-2) - 11./14.*dqt
                    xt = 0.5*(br(i,npy-1) - bl(i,npy-1))
           dm(i,npy-1) = sign(min(abs(xt), abs(dqt)), xt)
                    xt = 2.*dm(i,npy-1)
           bl(i,npy-1) = sign( min(abs(xt), abs(bl(i,npy-1))), xt )
           br(i,npy-1) = sign( min(abs(xt), abs(br(i,npy-1))), xt )
           br(i,npy-2) = 0.5*dqt + r3*(dm(i,npy-2) - dm(i,npy-1))
                    xt = 2.*dm(i,npy-2)
           br(i,npy-2) = sign(min(abs(xt), abs(br(i,npy-2))), xt)
           bl(i,npy-2) = sign(min(abs(xt), abs(al(i,npy-2)-q(i,npy-2))), xt)
! Outside: (mirrors what's done inside)
                 dqt = q(i,npy) - q(i,npy+1)
!          bl(i,npy) =  6./7.*dm(i,npy+1) + 13./14.*dqt
           bl(i,npy) =  q0 - q(i,npy)
           br(i,npy) = -4./7.*dm(i,npy+1) - 11./14.*dqt
                  xt = 0.5*(br(i,npy) - bl(i,npy))
           dm(i,npy) = sign(min(abs(xt), abs(dqt)), xt)
                  xt = 2.*dm(i,npy)
           bl(i,npy) = sign( min(abs(xt), abs(bl(i,npy))), xt )
           br(i,npy) = sign( min(abs(xt), abs(br(i,npy))), xt )
!          endif 
        enddo
      endif

   do j=jfirst,jlast+1
      do i=ifirst,ilast
         if(c(i,j)>0.) then
            flux(i,j) = q(i,j-1) + (1.-c(i,j))*(br(i,j-1)+c(i,j)*(bl(i,j-1)-br(i,j-1)))
         else
            flux(i,j) = q(i,j  ) - (1.+c(i,j))*(bl(i,j  )+c(i,j)*(bl(i,j  )-br(i,j  )))
         endif
      enddo
   enddo
 elseif( jord==9 ) then
! Based on 5
   do j=js-3,je+2
      do i=ifirst,ilast
         dq(i,j) = q(i,j+1) - q(i,j)
      enddo
   enddo

   if ( uniform_ppm ) then
        do j=js-2,je+2
           do i=ifirst,ilast
              xt = 0.25*(dq(i,j-1) + dq(i,j))
              dm(i,j) = sign(min(abs(xt), abs(dq(i,j-1)), abs(dq(i,j))), xt)
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

   do j=max(3,js-1),min(npy-2,je+2)
      do i=ifirst,ilast
         al(i,j) = 0.5*(q(i,j-1)+q(i,j)) + r3*(dm(i,j-1) - dm(i,j))
      enddo
   enddo

   do j=max(3,js-1),min(npy-3,je+1)
      do i=ifirst,ilast
            pmp = -2.*dq(i,j) 
            lac = pmp + 1.5*dq(i,j+1)
            bl(i,j) = min(max(0.,pmp,lac), max(al(i,j)-q(i,j), min(0.,pmp,lac)))
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
          xt = 27./28.*(q(i,0)+q(i,1)) - 13./28.*(q(i,-1)+q(i,2))   &
              + 3./7.*(dm(i,2)-dm(i,-1))
          bl(i,1) = xt - q(i,1)
          br(i,0) = xt - q(i,0)
          bl(i,0) = 11./14.*(q(i,-1)-q(i,0)) + 4./7.*dm(i,-1)
          xt = 3./14.*q(i,1) + 11./14.*q(i,2) - 4./7.*dm(i,2)
          bl(i,2) = xt - q(i,2)
          br(i,1) = xt - q(i,1)
          br(i,2) = al(i,3) - q(i,2)
       enddo
       do j=0,2
          call std_ppm(ilast-ifirst+1, q(ifirst,j), bl(ifirst,j), br(ifirst,j))
       enddo
   endif 

   if( (je+1)==npy ) then
       do i=ifirst,ilast
          bl(i,npy-2) = al(i,npy-2) - q(i,npy-2)
          xt = 3./14.*q(i,npy-1) + 11./14.*q(i,npy-2) + 4./7.*dm(i,npy-2)
          br(i,npy-2) = xt - q(i,npy-2)
          bl(i,npy-1) = xt - q(i,npy-1)
          br(i,npy) = 11./14.*(q(i,npy+1)-q(i,npy)) - 4./7.*dm(i,npy+1)
          xt = 27./28.*(q(i,npy-1)+q(i,npy)) - 13./28.*(q(i,npy-2)+q(i,npy+1))   &
              + 3./7.*(dm(i,npy+1)-dm(i,npy-2))
          br(i,npy-1) = xt - q(i,npy-1)
          bl(i,npy  ) = xt - q(i,npy)
       enddo
       do j=npy-2,npy
          call std_ppm(ilast-ifirst+1, q(ifirst,j), bl(ifirst,j), br(ifirst,j))
       enddo
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



 subroutine ytp(fy, q, c, yfx, ng, jord, &
                ifirst, ilast, jfirst, jlast, npx, npy, uniform_ppm)
! !INPUT PARAMETERS:
 integer, intent(in) :: npx, npy
 integer, INTENT(IN) :: ifirst, ilast  !  X-Dir strip
 integer, INTENT(IN) :: jfirst, jlast  !  Y-Dir strip
 integer, INTENT(IN) :: ng             !  Ghost Points
 integer, intent(in):: jord            !  order of subgrid dist
 logical, intent(IN) :: uniform_ppm
 real, intent(in)::   q(ifirst:ilast,jfirst-ng:jlast+ng) !  advected scalar
 real, intent(in)::   c(ifirst:ilast,jfirst:jlast+1)     !  Courant   N (like FY)
 real, intent(in):: yfx(ifirst:ilast,jfirst:jlast+1)     !  Backgrond flux

! !OUTPUT PARAMETERS:
 real, intent(out):: fy(ifirst:ilast,jfirst:jlast+1)     !  Flux
!
! !LOCAL VARIABLES:
 integer i, j, n, jt
! work arrays (should pass in eventually for performance enhancement):
 real   dm(ifirst:ilast,jfirst-ng:jlast+ng)

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
      call fyppm(c, q, dm, fy, ng, jord, ifirst,ilast,jfirst,jlast,  &
                 npx, npy, uniform_ppm)
   endif

   do j=jfirst,jlast+1
      do i=ifirst,ilast
         fy(i,j) = fy(i,j)*yfx(i,j)
      enddo
   enddo

 end subroutine ytp


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: xmist_2d
!                      
! !INTERFACE:
 subroutine xmist_2d(q, dm, ng, iord, ifirst, ilast, jfirst, jlast)
!-----------------------------------------------------------------------

! !USES:

! !INPUT PARAMETERS:
 integer, INTENT(IN) :: ifirst, ilast             !  X-Dir strip
 integer, INTENT(IN) :: jfirst, jlast             !  Y-Dir strip
 integer, INTENT(IN) :: ng                        !  Ghost Points 
 integer, INTENT(IN) :: iord                      !  order of subgrid distribution
 real   , INTENT(IN) :: q(ifirst-ng:ilast+ng,jfirst:jlast)  ! Ghosted transported scalar

! !OUTPUT PARAMETERS:
 real   , INTENT(INOUT) :: dm(ifirst-ng:ilast+ng,jfirst:jlast)  !  Ghosted Slope
!EOP
!---------------------------------------------------------------------
!BOC

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
!EOC
 end subroutine xmist_2d
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ymist
!                      
! !INTERFACE:
 subroutine ymist(q, dm, ng, jord, ifirst, ilast, jfirst, jlast)
!-----------------------------------------------------------------------
   
! !USES:
   
! !INPUT PARAMETERS:
 integer, INTENT(IN) :: ifirst, ilast             !  X-Dir strip
 integer, INTENT(IN) :: jfirst, jlast             !  Y-Dir strip
 integer, INTENT(IN) :: ng                        !  Ghost Points 
 integer, INTENT(IN) :: jord                      !  order of subgrid distribution
 real   , INTENT(IN) :: q(ifirst:ilast,jfirst-ng:jlast+ng)  ! Ghosted transported scalar

! !OUTPUT PARAMETERS:
 real   , INTENT(INOUT) :: dm(ifirst:ilast,jfirst-ng:jlast+ng)  !  Ghosted Slope
   
! !DESCRIPTION:
!     Calculate the slope of the pressure.  The number of ghost
!     latitudes (NG) depends on what method (JORD) will be used
!     subsequentally.    NG is equal to MIN(ABS(JORD),3).
!  
! !CALLED FROM:
!     ytp
!  
! !REVISION HISTORY:
!  
!  SJL 99.04.13:  Delivery
!  WS  99.04.13:  Added jfirst:jlast concept
!  WS  99.09.09:  Documentation; indentation; cleaning
!  SJL 00.01.06:  Documentation
!  WS  00.05.14:  Renamed ghost indices as per Kevin's definitions
!  BMP 05.08.08:  Removed lat/lon poles specifics to assume ghosting exists beyond poles
!  BMP 05.08.10:  Added ifirst:ilast concept
!
!EOP
!---------------------------------------------------------------------
!BOC

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
!-----------------------------------------------------------------------


!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_ghost_ew --- Ghost 4d east/west "lat/lon periodic
!
! !INTERFACE:
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

!EOC
 end subroutine mp_ghost_ew



 subroutine std_ppm(im, a0, al, ar)
 integer, intent(in):: im
 real, intent(in)   :: a0(im)
 real, intent(inout):: al(im), ar(im)
 real da1, da2, a6da
 integer i

! Optimized "Standard" PPM in perturbation form:
! effect of dm=0 not included
      do i=1,im
         da1 = al(i) - ar(i)
         da2 = da1**2
         a6da = 3.*(al(i)+ar(i))*da1
         if( a6da < -da2 ) then
             ar(i) = -2.*al(i)
         elseif( a6da > da2 ) then
             al(i) = -2.*ar(i)
         endif
      enddo
 end subroutine std_ppm


end module tpcore

