!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be
!* useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

! =======================================================================
! This module contains utilities for calculating flow deformation.
! Developers: Lucas Harris, Kun Gao
! =======================================================================

module deform_utils_mod

  use fv_arrays_mod,      only: fv_grid_bounds_type, fv_grid_type
  use constants_mod,      only: grav

  implicit none
  private
  public :: cal_3d_tke_budget

contains

  subroutine cal_3d_tke_budget(u, v, ua, va, w, &
                   delz, npz, ak, bk, gridstruct, bd, &
                   deform_1h, deform_1v)

! =======================================================================
! cal_3d_tke_budget :: calculate 3D TKE shear production
! - this code is created based on Ping Zhu's anisotropic
!   turbulence shear production formulation (Zhu et al. 2025; doi: 10.1038/s41612-025-01117-6)
!   and also his original code for 3D TKE shear production calculations
! - in this verion the native FV3 D-grid wind is used for
!   computing strain-rate tensor elements (Gao et al. 2026; doi: 10.1029/2026MS005916)
! - each element is strictly defined at FV3 D-grid cell center and layer center
! - the algorithms are updated to be consistent with finite-volume method
! - contact: Kun Gao (kun.gao@noaa.gov)
! =======================================================================

    integer, intent(in) :: npz
    type(fv_grid_bounds_type), intent(IN) :: bd
    type(fv_grid_type), intent(IN), target :: gridstruct
    real, intent(in) ::     ak(npz+1), bk(npz+1)
    real, intent(in) ::      u(bd%isd:bd%ied,  bd%jsd:bd%jed+1,npz)
    real, intent(in) ::      v(bd%isd:bd%ied+1,bd%jsd:bd%jed,  npz)
    real, intent(in) ::     ua(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    real, intent(in) ::     va(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    real, intent(in) ::      w(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    real, intent(in) ::   delz(bd%is:bd%ie,   bd%js:bd%je,   npz)

    real, intent(out) ::    deform_1h(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    real, intent(out) ::    deform_1v(bd%isd:bd%ied, bd%jsd:bd%jed, npz)

! local

    real, pointer, dimension(:,:) :: dx, dy, rarea, dxc, dyc, rarea_c, dxa, dya
    integer :: is, ie, js, je, isd, ied, jsd, jed, i, j, k
    real:: dudx(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: dudy(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: dvdx(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: dvdy(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: dwdx(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: dwdy(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: dudz(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: dvdz(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: dwdz(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: ut(bd%isd:bd%ied+1,bd%jsd:bd%jed)
    real:: vt(bd%isd:bd%ied,  bd%jsd:bd%jed+1)
    real:: wk1(bd%isd:bd%ied,bd%jsd:bd%jed) 
    real:: wk2(bd%isd:bd%ied,bd%jsd:bd%jed)
    real:: u_e(bd%is:bd%ie,npz+1)
    real:: v_e(bd%is:bd%ie,npz+1)
    real:: w_e(bd%is:bd%ie,npz+1)
    
    real :: dp_ref(npz)
    real :: tmp, dz

    real, parameter:: p1 =  7./12.
    real, parameter:: p2 = -1./12.

    dxc => gridstruct%dxc
    dyc => gridstruct%dyc
    dx  => gridstruct%dx
    dy  => gridstruct%dy
    dxa  => gridstruct%dxa
    dya  => gridstruct%dya
    rarea   => gridstruct%rarea
    rarea_c => gridstruct%rarea_c

    is  = bd%is
    ie  = bd%ie
    js  = bd%js
    je  = bd%je
    isd  = bd%isd
    ied  = bd%ied
    jsd  = bd%jsd
    jed  = bd%jed

!===========================================================
! Calculate deform_1h and deform_1v
!
! See Eq.11 in Zhu et al. 2025
!
!      shr_prod = kh * deform_1h + kv * deform_1v 
!
!     deform_1h = (2*(du/dx**2 + dv/dy**2)
!                  + (du/dy + dv/dx)**2
!                  + (dw/dx**2 + dw/dy**2)
!                  + (dw/dx*du/dz + dw/dy*dv/dz))
!
!     deform_1v = (2*dw/dz**2 
!                  + (du/dz**2 + dv/dz**2)
!                  + (dw/dx*du/dz + dw/dy*dv/dz))
!===========================================================

!$OMP parallel do default(none) shared(npz,is,ie,js,je,u,v,w, &
!$OMP                                  dx,dy,rarea,dxc,dyc,rarea_c,dxa,dya, &
!$OMP                                  dudx,dudy,dvdx,dvdy,dwdx,dwdy) &
!$OMP                          private(ut,vt,wk1,wk2)
   do k=1,npz

!---------------------------------------------
! Calculate du/dx and dv/dy:
!  get corner values first, and 
!  then average the 4 corners to get cell-mean
!---------------------------------------------

      do j=js-1,je+2
         do i=is-2,ie+2
            ut(i,j) = u(i,j,k)*dyc(i,j)
         enddo
      enddo
      do j=js-2,je+2
         do i=is-1,ie+2
            vt(i,j) = v(i,j,k)*dxc(i,j)
         enddo
      enddo
      do j=js,je+1
         do i=is,ie+1
            wk1(i,j) = rarea_c(i,j)*( -ut(i-1,j) + ut(i,j) ) ! du/dx at corners
            wk2(i,j) = rarea_c(i,j)*( -vt(i,j-1) + vt(i,j) ) ! dv/dy at corners
         enddo
      enddo
! Fix the corners if grid_type /= 4
      do j=js,je
         do i=is,ie
            ! mean of 4 corners
            dudx(i,j,k) = 0.25*(wk1(i,j) + wk1(i,j+1) + wk1(i+1,j) + wk1(i+1,j+1))
            dvdy(i,j,k) = 0.25*(wk2(i,j) + wk2(i,j+1) + wk2(i+1,j) + wk2(i+1,j+1))
         enddo
      enddo

!-------------------------------------
! Calculate du/dy and dv/dx
!-------------------------------------

      do j=js,je+1
         do i=is,ie
            vt(i,j) = u(i,j,k)*dx(i,j)
         enddo
      enddo
      do j=js,je
         do i=is,ie+1
            ut(i,j) = v(i,j,k)*dy(i,j)
         enddo
      enddo

      do j=js,je
         do i=is,ie
            dudy(i,j,k) = rarea(i,j)*( -vt(i,j) + vt(i,j+1) )
            dvdx(i,j,k) = rarea(i,j)*( -ut(i,j) + ut(i+1,j) )
         enddo
      enddo

!-------------------------------------------
! Calculate dw/dx and dw/dy:
!  get edge values of w first, and then 
!  use the two edge values to get cell-mean 
!  dw/dx and dw/dy
!-------------------------------------------

      do j=js-1,je
         do i=is-1,ie
           ! using the 1d ppm method to get w values at cell edges
           ! q(i+0.5) = p1*(q(i-1)+q(i)) + p2*(q(i-2)+q(i+1))
           wk1(i+1,j) = p1*(w(i-1,j,k)+w(i,j,k)) + p2*(w(i-2,j,k)+w(i+1,j,k))
           wk2(i,j+1) = p1*(w(i,j-1,k)+w(i,j,k)) + p2*(w(i,j-2,k)+w(i,j+1,k))
         enddo
      enddo

      do j=js,je
         do i=is,ie
            dwdx(i,j,k) = 1./dxa(i,j) * (-wk1(i,j) + wk1(i+1,j))
            dwdy(i,j,k) = 1./dya(i,j) * (-wk2(i,j) + wk2(i,j+1))
         enddo
      enddo

   enddo   !z loop

!-------------------------------------
! get du/dz, dv/dz and dw/dz
!-------------------------------------

   do k=1,npz
      dp_ref(k) = ak(k+1)-ak(k) + (bk(k+1)-bk(k))*1.E5
   enddo

   do j=js,je
      call edge_profile1(ua(is:ie,j,:), u_e, is, ie, npz, dp_ref, 1)
      call edge_profile1(va(is:ie,j,:), v_e, is, ie, npz, dp_ref, 1)
      call edge_profile1(w(is:ie,j,:),  w_e, is, ie, npz, dp_ref, 1)
      do k=1,npz
         do i=is,ie
            dz = delz(i,j,k) !zh(i,j,k) - zh(i,j,k+1) 
            dudz(i,j,k) = (u_e(i,k)-u_e(i,k+1))/dz
            dvdz(i,j,k) = (v_e(i,k)-v_e(i,k+1))/dz
            dwdz(i,j,k) = (w_e(i,k)-w_e(i,k+1))/dz
         enddo
      enddo
   enddo

!-----------------------------------------------
! get deform_1h and deform_1v based on all terms
!-----------------------------------------------

!$OMP parallel do default(none) shared(npz,is,ie,js,je, &
!$OMP                                  deform_1h, deform_1v, &
!$OMP                                  dudx,dudy,dudz,dvdx, &
!$OMP                                  dvdy,dvdz,dwdx,dwdy,dwdz) &
!$OMP                          private(tmp)
   do k=1,npz
       do j=js,je
          do i=is,ie

             tmp = dwdx(i,j,k)*dudz(i,j,k)+dwdy(i,j,k)*dvdz(i,j,k)

             deform_1h(i,j,k)=                            &
               (2*(dudx(i,j,k)**2 + dvdy(i,j,k)**2)       &
                + (dudy(i,j,k) + dvdx(i,j,k))**2+         &
                + dwdx(i,j,k)**2 + dwdy(i,j,k)**2 + tmp)


             deform_1v(i,j,k)=                            &
               (2*dwdz(i,j,k)**2                          &
                + dudz(i,j,k)**2 + dvdz(i,j,k)**2 + tmp)

             ! apply limiters
             deform_1h(i,j,k) = min( 1.e3, deform_1h(i,j,k) )
             !deform_1v(i,j,k) = min( 1.e3, deform_1v(i,j,k) )
          enddo
       enddo
   enddo

  end subroutine cal_3d_tke_budget

  subroutine compute_dudz(bd, npz, u, v, dudz, dvdz, gz, dp_ref)

 !routine to compute vertical gradients in winds
 ! for 2D smag damping
 ! Call AFTER updating gz
 !TODO needs cubed-sphere support (don't compute in corners)

   type(fv_grid_bounds_type), intent(IN) :: bd
   integer, intent(IN) :: npz
   real, intent(in) :: u(bd%isd:bd%ied,  bd%jsd:bd%jed+1,npz)
   real, intent(in) :: v(bd%isd:bd%ied+1,bd%jsd:bd%jed,  npz)
   real, intent(in) :: gz(bd%isd:bd%ied, bd%jsd:bd%jed,  npz+1)
   real, intent(IN) :: dp_ref(npz)
   real, intent(OUT) :: dudz(bd%isd:bd%ied,bd%jsd:bd%jed+1,npz)
   real, intent(OUT) :: dvdz(bd%isd:bd%ied+1,bd%jsd:bd%jed,npz)

   real :: dz
   real :: ue(bd%isd:bd%ied  ,npz+1)
   real :: ve(bd%isd:bd%ied+1,npz+1)
   integer :: i,j,k
   integer :: is,  ie,  js,  je
   integer :: isd, ied, jsd, jed

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
   isd  = bd%isd
   ied  = bd%ied
   jsd  = bd%jsd
   jed  = bd%jed

   dudz = -1.e50
   dvdz = -1.e50

   do j=jsd,jed
      !TODO: pass by reference and not copy
      call edge_profile1(v(isd:ied+1,j,:), ve, isd,  ied+1, npz, dp_ref, 0)
      do k=1,npz
         do i=isd+1,ied
            dz = gz(i,j,k) + gz(i-1,j,k)
            dz = dz - (gz(i,j,k+1) + gz(i-1,j,k+1))
            dz = 0.5*dz/grav
            dvdz(i,j,k) = (ve(i,k)-ve(i,k+1))/dz
         enddo
      enddo
   enddo

   do j=jsd+1,jed
      call edge_profile1(u(isd:ied,j,:), ue, isd, ied, npz, dp_ref, 0)
      do k=1,npz
         do i=isd,ied
            dz = gz(i,j,k) + gz(i,j-1,k)
            dz = dz - (gz(i,j,k+1) + gz(i,j-1,k+1))
            dz = 0.5*dz/grav
            dudz(i,j,k) = (ue(i,k)-ue(i,k+1))/dz
         enddo
      enddo
   enddo

  end subroutine compute_dudz

  subroutine edge_profile1(q1, q1e, i1, i2, km, dp0, limiter)
  ! Edge profiles for a single scalar quantity
  ! also see edge_profile in model/nh_utils.F90
  integer, intent(in):: i1, i2
  integer, intent(in):: km
  integer, intent(in):: limiter
  real, intent(in):: dp0(km)
  real, intent(in),  dimension(i1:i2,km):: q1
  real, intent(out), dimension(i1:i2,km+1):: q1e
  !-----------------------------------------------------------------------
  real, dimension(i1:i2,km+1):: qe1, gam  ! edge values
  real  gak(km)
  real  bet, r2o3, r4o3
  real  g0, gk, xt1, xt2, a_bot
  integer i, k

  ! Assuming grid varying in vertical only
   g0 = dp0(2) / dp0(1)
  xt1 = 2.*g0*(g0+1. )
  bet =    g0*(g0+0.5)
  do i=i1,i2
      qe1(i,1) = ( xt1*q1(i,1) + q1(i,2) ) / bet
      gam(i,1) = ( 1. + g0*(g0+1.5) ) / bet
  enddo

  do k=2,km
     gk = dp0(k-1) / dp0(k)
     do i=i1,i2
             bet =  2. + 2.*gk - gam(i,k-1)
        qe1(i,k) = ( 3.*(q1(i,k-1)+gk*q1(i,k)) - qe1(i,k-1) ) / bet
        gam(i,k) = gk / bet
     enddo
  enddo

  a_bot = 1. + gk*(gk+1.5)
    xt1 =   2.*gk*(gk+1.)
  do i=i1,i2
             xt2 = gk*(gk+0.5) - a_bot*gam(i,km)
     qe1(i,km+1) = ( xt1*q1(i,km) + q1(i,km-1) - a_bot*qe1(i,km) ) / xt2
  enddo

  do k=km,1,-1
     do i=i1,i2
        qe1(i,k) = qe1(i,k) - gam(i,k)*qe1(i,k+1)
     enddo
  enddo

  !------------------
  ! Apply constraints
  !------------------
    if ( limiter/=0 ) then   ! limit the top & bottom winds
         do i=i1,i2
  ! Top
            if ( q1(i,1)*qe1(i,1) < 0. ) qe1(i,1) = 0.
  ! Surface:
            if ( q1(i,km)*qe1(i,km+1) < 0. ) qe1(i,km+1) = 0.
         enddo
    endif

    do k=1,km+1
       do i=i1,i2
          q1e(i,k) = qe1(i,k)
       enddo
    enddo

  end subroutine edge_profile1

end module deform_utils_mod
