module a2b_edge_mod

  use grid_utils, only : edge_w, edge_e, edge_s, edge_n, sw_corner, se_corner,  &
                         nw_corner, ne_corner
  use grid_tools, only: dxa, dya, grid_type
  implicit none

!----------------------------
! 4-pt Lagrange interpolation
!----------------------------
  real, parameter:: a1 =  0.5625
  real, parameter:: a2 = -0.0625
!----------------------
! PPM volume mean form:
!----------------------
  real, parameter:: b1 =  7./12.     ! 0.58333333
  real, parameter:: b2 = -1./12.

  private
  public :: a2b_ord2, a2b_ord4

contains

  subroutine a2b_ord4(qin, qout, npx, npy, is, ie, js, je, ng, replace)
  integer, intent(IN):: npx, npy, is, ie, js, je, ng
  real, intent(INOUT)::  qin(is-ng:ie+ng,js-ng:je+ng)   ! A-grid field
  real, intent(INOUT):: qout(is-ng:ie+ng,js-ng:je+ng)   ! Output  B-grid field
  logical, optional, intent(IN):: replace
! local:
  real, parameter:: r3 = 1./3.
  real qx(is:ie+1,js-ng:je+ng)
  real qy(is-ng:ie+ng,js:je+1)
  real qxx(is-ng:ie+ng,js-ng:je+ng)
  real qyy(is-ng:ie+ng,js-ng:je+ng)
  real gratio
  integer :: i, j, is1, js1, is2, js2, ie1, je1


  if (grid_type < 3) then

    is1 = max(1,is-1)
    js1 = max(1,js-1)
    is2 = max(2,is)
    js2 = max(2,js)

    ie1 = min(npx-1,ie+1)
    je1 = min(npy-1,je+1)

    if ( sw_corner ) qout(1,    1) = r3*(qin(1,        1)+qin(1,      0)+qin(0,      1))
    if ( se_corner ) qout(npx,  1) = r3*(qin(npx-1,    1)+qin(npx-1,  0)+qin(npx,    1))
    if ( ne_corner ) qout(npx,npy) = r3*(qin(npx-1,npy-1)+qin(npx,npy-1)+qin(npx-1,npy))
    if ( nw_corner ) qout(1,  npy) = r3*(qin(1,    npy-1)+qin(0,  npy-1)+qin(1,    npy))

! West Edges:
    if ( is==1 ) then
       do j=max(1,js-2),min(npy-1,je+2)
           gratio = dxa(2,j) / dxa(1,j)
          qx(1,j) = 0.5*((2.+gratio)*(qin(0,j)+qin(1,j))    &
                  - (qin(-1,j)+qin(2,j))) / (1.+gratio)
          qx(2,j) = (2.*gratio*(gratio+1.)*qin(1,j)+qin(2,j) -     &
                     gratio*(gratio+0.5)*qx(1,j))/(1.+gratio*(gratio+1.5))
       enddo
!      do j=max(2,js),min(npy-1,je+1)
!         qout(1,j) = 0.5*(qx(1,j-1) + qx(1,j))
!      enddo
       if( js==1 ) qout(1,2) = 0.5*(qx(1,1) + qx(1,2))
       do j=max(3,js),min(npy-2,je+1)
          qout(1,j) = a1*(qx(1,j-1)+qx(1,j)) + a2*(qx(1,j-2)+qx(1,j+1))
       enddo
       if( (je+1)==npy) qout(1,npy-1) = 0.5*(qx(1,npy-2) + qx(1,npy-1))
    endif

! East Edges:
    if ( (ie+1)==npx ) then
       do j=max(1,js-2),min(npy-1,je+2)
               gratio = dxa(npx-2,j) / dxa(npx-1,j)
          qx(npx  ,j) = 0.5*((2.+gratio)*(qin(npx-1,j)+qin(npx,j))   &
                        - (qin(npx-2,j)+qin(npx+1,j))) / (1.+gratio )
          qx(npx-1,j) = (2.*gratio*(gratio+1.)*qin(npx-1,j)+qin(npx-2,j) -  &
                         gratio*(gratio+0.5)*qx(npx,j))/(1.+gratio*(gratio+1.5))
       enddo
!      do j=max(2,js),min(npy-1,je+1)
!         qout(npx,j) = 0.5*(qx(npx,j-1) + qx(npx,j))
!      enddo
       if( js==1 ) qout(npx,2) = 0.5*(qx(npx,1) + qx(npx,2))
       do j=max(3,js),min(npy-2,je+1)
          qout(npx,j) = a1*(qx(npx,j-1)+qx(npx,j)) + a2*(qx(npx,j-2)+qx(npx,j+1))
       enddo
       if( (je+1)==npy ) qout(npx,npy-1) = 0.5*(qx(npx,npy-2) + qx(npx,npy-1))
    endif

! South Edges:
    if ( js==1 ) then
       do i=max(1,is-2),min(npx-1,ie+2)
           gratio = dya(i,2) / dya(i,1)
          qy(i,1) = 0.5*((2.+gratio)*(qin(i,0)+qin(i,1))   &
                  - (qin(i,-1)+qin(i,2))) / (1.+gratio )
          qy(i,2) = (2.*gratio*(gratio+1.)*qin(i,1)+qin(i,2) -     &
                     gratio*(gratio+0.5)*qy(i,1))/(1.+gratio*(gratio+1.5))
       enddo
!      do i=max(2,is),min(npx-1,ie+1)
!         qout(i,1) = 0.5*(qy(i-1,1) + qy(i,1))
!      enddo
       if( is==1 ) qout(2,1) = 0.5*(qy(1,1) + qy(2,1))
       do i=max(3,is),min(npx-2,ie+1)
          qout(i,1) = a1*(qy(i-1,1)+qy(i,1)) + a2*(qy(i-2,1)+qy(i+1,1))
       enddo
       if( (ie+1)==npx ) qout(npx-1,1) = 0.5*(qy(npx-2,1) + qy(npx-1,1))
    endif

    ! North Edges:
    if ( (je+1)==npy ) then
       do i=max(1,is-2),min(npx-1,ie+2)
               gratio = dya(i,npy-2) / dya(i,npy-1)
          qy(i,npy  ) = 0.5*((2.+gratio)*(qin(i,npy-1)+qin(i,npy))  &
                      - (qin(i,npy-2)+qin(i,npy+1))) / (1.+gratio)
          qy(i,npy-1) = (2.*gratio*(gratio+1.)*qin(i,npy-1)+qin(i,npy-2) - &
                         gratio*(gratio+0.5)*qy(i,npy))/(1.+gratio*(gratio+1.5))
       enddo
!      do i=max(2,is),min(npx-1,ie+1)
!         qout(i,npy) = 0.5*(qy(i-1,npy) + qy(i,npy))
!      enddo
       if( is==1 ) qout(2,npy) = 0.5*(qy(1,npy) + qy(2,npy))
       do i=max(3,is),min(npx-2,ie+1)
          qout(i,npy) = a1*(qy(i-1,npy)+qy(i,npy)) + a2*(qy(i-2,npy)+qy(i+1,npy))
       enddo
       if( (ie+1)==npx ) qout(npx-1,npy) = 0.5*(qy(npx-2,npy) + qy(npx-1,npy))
    endif

!----------
! Interior:
!----------
! X-sweep:
    do j=max(1,js-2),min(npy-1,je+2)
       do i=max(3,is), min(npx-2,ie+1)
          qx(i,j) = a1*(qin(i-1,j)+qin(i,j)) + a2*(qin(i-2,j)+qin(i+1,j))
       enddo
    enddo
    
    do j=max(2,js),min(npy-1,je+1)
       if ( j==2 .or. j==(npy-1) ) then
            do i=max(2,is),min(npx-1,ie+1)
               qxx(i,j) = 0.5*(qx(i,j-1)+qx(i,j))
            enddo
       else
            do i=max(2,is),min(npx-1,ie+1)
               qxx(i,j) = a1*(qx(i,j-1)+qx(i,j)) + a2*(qx(i,j-2)+qx(i,j+1))
            enddo
       endif
    enddo

! Y-sweep:
    do j=max(3,js),min(npy-2,je+1)
       do i=max(1,is-2), min(npx-1,ie+2)
          qy(i,j) = a1*(qin(i,j-1)+qin(i,j)) + a2*(qin(i,j-2)+qin(i,j+1))
       enddo
    enddo
    
    do j=max(2,js),min(npy-1,je+1)
       if ( is==1 ) qyy(2,j) = 0.5*(qy(1,j)+qy(2,j))
       do i=max(3,is),min(npx-2,ie+1)
          qyy(i,j) = a1*(qy(i-1,j)+qy(i,j)) + a2*(qy(i-2,j)+qy(i+1,j))
       enddo
       if( (ie+1)==npx ) qyy(npx-1,j) = 0.5*(qy(npx-2,j)+qy(npx-1,j))
    enddo
 
!Avg:
    do j=max(2,js),min(npy-1,je+1)
    do i=max(2,is),min(npx-1,ie+1)
       qout(i,j) = 0.5*(qxx(i,j)+qyy(i,j))
    enddo
    enddo
    
 else  ! grid_type>=3

! X-sweep: PPM
    do j=js-2,je+2
       do i=is,ie+1
          qx(i,j) = b1*(qin(i-1,j)+qin(i,j)) + b2*(qin(i-2,j)+qin(i+1,j))
       enddo
    enddo
    
    do j=js,je+1
       do i=is,ie+1
          qxx(i,j) = a1*(qx(i,j-1)+qx(i,j)) + a2*(qx(i,j-2)+qx(i,j+1))
       enddo
    enddo

! Y-sweep: PPM
    do j=js,je+1
       do i=is-2,ie+2
          qy(i,j) = b1*(qin(i,j-1)+qin(i,j)) + b2*(qin(i,j-2)+qin(i,j+1))
       enddo
    enddo
    
    do j=js,je+1
       do i=is,ie+1
          qyy(i,j) = a1*(qy(i-1,j)+qy(i,j)) + a2*(qy(i-2,j)+qy(i+1,j))
       enddo
    enddo
 
!Avg:
    do j=js,je+1
       do i=is,ie+1
          qout(i,j) = 0.5*(qxx(i,j)+qyy(i,j))
       enddo
    enddo
    
 endif

    if ( present(replace) ) then
       if ( replace ) then
          do j=js,je+1
          do i=is,ie+1
             qin(i,j) = qout(i,j)
          enddo
          enddo
       endif
    endif
    
  end subroutine a2b_ord4


  subroutine a2b_ord2(qin, qout, npx, npy, is, ie, js, je, ng, replace)
    integer, intent(IN   ) :: npx, npy, is, ie, js, je, ng
    real   , intent(INOUT) ::  qin(is-ng:ie+ng,js-ng:je+ng)   ! A-grid field
    real   , intent(  OUT) :: qout(is-ng:ie+ng,js-ng:je+ng)   ! Output  B-grid field
    logical, optional, intent(IN) ::  replace
    ! local:
    real q1(npx), q2(npy)
    integer :: i,j
    integer :: is1, js1, is2, js2, ie1, je1
    real, parameter:: r3 = 1./3.

    if (grid_type < 3) then

    is1 = max(1,is-1)
    js1 = max(1,js-1)
    is2 = max(2,is)
    js2 = max(2,js)

    ie1 = min(npx-1,ie+1)
    je1 = min(npy-1,je+1)

    do j=js2,je1
       do i=is2,ie1
          qout(i,j) = 0.25*(qin(i-1,j-1)+qin(i,j-1)+qin(i-1,j)+qin(i,j))
       enddo
    enddo

! Fix the 4 Corners:
    if ( sw_corner ) qout(1,    1) = r3*(qin(1,        1)+qin(1,      0)+qin(0,      1))
    if ( se_corner ) qout(npx,  1) = r3*(qin(npx-1,    1)+qin(npx-1,  0)+qin(npx,    1))
    if ( ne_corner ) qout(npx,npy) = r3*(qin(npx-1,npy-1)+qin(npx,npy-1)+qin(npx-1,npy))
    if ( nw_corner ) qout(1,  npy) = r3*(qin(1,    npy-1)+qin(0,  npy-1)+qin(1,    npy))

    ! *** West Edges:
    if ( is==1 ) then
       do j=js1, je1
          q2(j) = 0.5*(qin(0,j) + qin(1,j))
       enddo
       do j=js2, je1
          qout(1,j) = edge_w(j)*q2(j-1) + (1.-edge_w(j))*q2(j)
       enddo
    endif

    ! East Edges:
    if ( (ie+1)==npx ) then
       do j=js1, je1
          q2(j) = 0.5*(qin(npx-1,j) + qin(npx,j))
       enddo
       do j=js2, je1
          qout(npx,j) = edge_e(j)*q2(j-1) + (1.-edge_e(j))*q2(j)
       enddo
    endif

    ! South Edges:
    if ( js==1 ) then
       do i=is1, ie1
          q1(i) = 0.5*(qin(i,0) + qin(i,1))
       enddo
       do i=is2, ie1
          qout(i,1) = edge_s(i)*q1(i-1) + (1.-edge_s(i))*q1(i)
       enddo
    endif

    ! North Edges:
    if ( (je+1)==npy ) then
       do i=is1, ie1
          q1(i) = 0.5*(qin(i,npy-1) + qin(i,npy))
       enddo
       do i=is2, ie1
          qout(i,npy) = edge_n(i)*q1(i-1) + (1.-edge_n(i))*q1(i)
       enddo
    endif

 else

    do j=js,je+1
       do i=is,ie+1
          qout(i,j) = 0.25*(qin(i-1,j-1)+qin(i,j-1)+qin(i-1,j)+qin(i,j))
       enddo
    enddo

 endif

    
    if ( present(replace) ) then
       if ( replace ) then
          do j=js,je+1
             do i=is,ie+1
                qin(i,j) = qout(i,j)
             enddo
          enddo
       endif
    endif
    
  end subroutine a2b_ord2
  
end module a2b_edge_mod
