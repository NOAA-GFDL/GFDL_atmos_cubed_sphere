module a2b_edge_mod

  use grid_utils, only : edge_w, edge_e, edge_s, edge_n, sw_corner, se_corner,  &
                         nw_corner, ne_corner
  implicit none

  private
  public :: a2b_edge

contains

  subroutine a2b_edge(qin, qout, npx, npy, is, ie, js, je, ng, replace)

    integer, intent(IN   ) :: npx, npy, is, ie, js, je, ng
    real   , intent(INOUT) ::  qin(is-ng:ie+ng,js-ng:je+ng)   ! A-grid field
    real   , intent(  OUT) :: qout(is-ng:ie+ng,js-ng:je+ng)   ! Output  B-grid field
    logical, optional, intent(IN) ::  replace
    ! local:
    real q1(npx), q2(npy)
    integer :: i,j
    integer :: is1, js1, is2, js2, ie1, je1
    real, parameter:: r3 = 1./3.

    is1 = max(1,is-1)
    js1 = max(1,js-1)

    is2 = max(2,is)
    js2 = max(2,js)

    ie1 = min(npx-1,ie+1)
    je1 = min(npy-1,je+1)

#ifdef TEST_A2B
    do j=js,je+1
       do i=is,ie+1
          qout(i,j) = 0.25*(qin(i-1,j-1)+qin(i,j-1)+qin(i-1,j)+qin(i,j))
       enddo
    enddo
#else
    do j=js2,je1
       do i=is2,ie1
          qout(i,j) = 0.25*(qin(i-1,j-1)+qin(i,j-1)+qin(i-1,j)+qin(i,j))
       enddo
    enddo

    ! *** West Edges:
    if ( is==1 ) then
       i=1
       do j=js1, je1
          q2(j) = 0.5*(qin(i-1,j) + qin(i,j))
       enddo
       do j=js2, je1
          qout(i,j) = edge_w(j)*q2(j-1) + (1.-edge_w(j))*q2(j)
       enddo
    endif

    ! East Edges:
    if ( (ie+1)==npx ) then
       i=npx
       do j=js1, je1
          q2(j) = 0.5*(qin(i-1,j) + qin(i,j))
       enddo
       do j=js2, je1
          qout(i,j) = edge_e(j)*q2(j-1) + (1.-edge_e(j))*q2(j)
       enddo
    endif

    ! South Edges:
    if ( js==1 ) then
       j=1
       do i=is1, ie1
          q1(i) = 0.5*(qin(i,j-1) + qin(i,j))
       enddo
       do i=is2, ie1
          qout(i,j) = edge_s(i)*q1(i-1) + (1.-edge_s(i))*q1(i)
       enddo
    endif

    ! North Edges:
    if ( (je+1)==npy ) then
       j=npy
       do i=is1, ie1
          q1(i) = 0.5*(qin(i,j-1) + qin(i,j))
       enddo
       do i=is2, ie1
          qout(i,j) = edge_n(i)*q1(i-1) + (1.-edge_n(i))*q1(i)
       enddo
    endif
#endif

! Fix the 4 Corners:
    if ( sw_corner ) qout(1,    1) = r3*(qin(1,        1)+qin(1,      0)+qin(0,      1))
    if ( se_corner ) qout(npx,  1) = r3*(qin(npx-1,    1)+qin(npx-1,  0)+qin(npx,    1))
    if ( ne_corner ) qout(npx,npy) = r3*(qin(npx-1,npy-1)+qin(npx,npy-1)+qin(npx-1,npy))
    if ( nw_corner ) qout(1,  npy) = r3*(qin(1,    npy-1)+qin(0,  npy-1)+qin(1,    npy))
    
    if ( present(replace) ) then
       if ( replace ) then
          do j=js,je+1
             do i=is,ie+1
                qin(i,j) = qout(i,j)
             enddo
          enddo
       endif
    endif
    
  end subroutine a2b_edge
  
end module a2b_edge_mod
