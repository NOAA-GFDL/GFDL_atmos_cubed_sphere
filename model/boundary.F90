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
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

!>@brief The module 'boundary' contains utility routines for grid nesting
!! and boundary conditions.

module boundary_mod

! Modules Included:
! <table>
! <tr>
!    <th>Module Name</th>
!     <th>Functions Included</th>
!  </tr>
!   <tr>
!     <td>constants_mod</td>
!     <td>grav</td>
!   <tr>
!     <td>fv_arrays_mod</td>
!     <td>fv_atmos_type, fv_nest_BC_type_3D, fv_grid_bounds_type</td>
!   </tr>
!   <tr>
!     <td>fv_mp_mod</td>
!     <td>ng, isc,jsc,iec,jec, isd,jsd,ied,jed, is,js,ie,je, is_master, mp_bcst</td>
!   </tr>
!   <tr>
!     <td>fv_timing_mod</td>
!     <td>timing_on, timing_off</td>
!   </tr>
!   <tr>
!     <td>mpp_mod/td>
!     <td>mpp_error, FATAL, mpp_sum, mpp_sync, mpp_npes, mpp_broadcast, WARNING, mpp_pe,
!         mpp_send, mpp_recv</td>
!   </tr>
!   <tr>
!     <td>mpp_domains_mod/td>
!     <td>mpp_get_compute_domain, mpp_get_data_domain, mpp_get_global_domain,
!         ENTER, CORNER, NORTH, EAST,nest_domain_type, WEST, SOUTH,
!         mpp_get_C2F_index, mpp_update_nest_fine,mpp_global_field, mpp_get_pelist
!         mpp_get_F2C_index, mpp_update_nest_coarse</td>
!   </tr>
! </table>

  use fv_mp_mod,         only: is_master
  use constants_mod,     only: grav

  use mpp_domains_mod,    only: mpp_get_compute_domain, mpp_get_data_domain, mpp_get_global_domain
  use mpp_domains_mod,    only: CENTER, CORNER, NORTH, EAST
  use mpp_domains_mod,    only: mpp_global_field, mpp_get_pelist
  use mpp_domains_mod,    only: AGRID, BGRID_NE, CGRID_NE, DGRID_NE
  use mpp_mod,            only: mpp_error, FATAL, mpp_sum, mpp_sync, mpp_npes, mpp_broadcast, WARNING, mpp_pe

  use fv_mp_mod,          only: mp_bcst
  use fv_arrays_mod,      only: fv_atmos_type, fv_nest_BC_type_3D, fv_grid_bounds_type
  use mpp_mod,            only: mpp_send, mpp_recv
  use fv_timing_mod,      only: timing_on, timing_off
  use mpp_domains_mod, only : nest_domain_type, WEST, SOUTH
  use mpp_domains_mod, only : mpp_get_C2F_index, mpp_update_nest_fine
  use mpp_domains_mod, only : mpp_get_F2C_index, mpp_update_nest_coarse
  !use mpp_domains_mod, only : mpp_get_domain_shift

  implicit none
  public extrapolation_BC
  public nested_grid_bc, update_coarse_grid
  public fill_nested_grid, nested_grid_BC_apply_intT
  public nested_grid_BC_send, nested_grid_BC_recv, nested_grid_BC_save_proc

!>@briefThe interface 'nested_grid_BC' includes subroutines 'nested_grid_BC_2d' and 'nested_grid_BC_3d'
!! that fetch coarse-grid data, interpolate it to nested-grid boundary cells,
!! apply the interpolated data directly to the boundary halo cells without saving the datatype.
  interface nested_grid_BC
     module procedure nested_grid_BC_2d
!     module procedure nested_grid_BC_mpp_2d
     module procedure nested_grid_BC_mpp_3d
     module procedure nested_grid_BC_mpp_send_2d
     module procedure nested_grid_BC_mpp_send_3d
     module procedure nested_grid_BC_2D_mpp
     module procedure nested_grid_BC_3d
     module procedure nested_grid_BC_mpp_3d_vector
  end interface

  interface nested_grid_BC_send
     module procedure nested_grid_BC_send_scalar
     module procedure nested_grid_BC_send_vector
  end interface

  interface nested_grid_BC_recv
     module procedure nested_grid_BC_recv_scalar
     module procedure nested_grid_BC_recv_vector
  end interface
!>@brief The interface 'fill_nested_grid' includes subroutines 'fill_nested_grid_2d' and 'fill_nested_grid_3d'
!! that fill nested-grid data with interpolated data from the coarse grid.
!>@details This is one method to create a new nested grid, and may be useful when cold-starting.

  interface fill_nested_grid
     module procedure fill_nested_grid_2d
     module procedure fill_nested_grid_3d
  end interface

!>@brief The interface'update_coarse_grid_mpp'contains subroutines that
!! fetch data from the nested grid and
!! interpolate it to the coarse grid using the method described by
!! \cite harris2013two.
  interface update_coarse_grid
     module procedure update_coarse_grid_mpp
     module procedure update_coarse_grid_mpp_2d
     module procedure update_coarse_grid_mpp_vector
  end interface

contains

!>@brief The subroutine 'extrapolation_BC' performs linear extrapolation into the halo region.
!Not to be confused with extrapolated-in-time nested BCs
  subroutine extrapolation_BC(q, istag, jstag, npx, npy, bd, pd_in, debug_in)

    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(in) :: istag, jstag, npx, npy
    real, intent(inout), dimension(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag) :: q
    logical, intent(in), OPTIONAL :: pd_in, debug_in

    integer :: i,j, istart, iend, jstart, jend
    logical :: pd, debug

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

    istart = max(isd, 1)
    iend = min(ied,npx-1)
    jstart = max(jsd, 1)
    jend = min(jed,npy-1)

    !Positive-definite extrapolation: shift from linear extrapolation to zero-gradient when the extrapolated value turns negative.
    if (present(pd_in)) then
       pd    = pd_in
    else
       pd = .false.
    end if

    if (present(debug_in)) then
       debug = debug_in
    else
       debug = .false.
    end if

    if (is == 1) then

       if (pd) then

          do j = jstart,jend+jstag
          do i = 0,isd,-1

             if  (real(i) <= 1. - q(1,j)/(q(2,j) - q(1,j) + 1.e-12) .and. q(1,j) < q(2,j)) then
                q(i,j) = q(i+1,j)
             else
                q(i,j) = real(2-i)*q(1,j) - real(1-i)*q(2,j)
             end if

          end do
          end do

       else

          do j = jstart,jend+jstag
          do i = 0,isd,-1

             q(i,j) = real(2-i)*q(1,j) - real(1-i)*q(2,j)

          end do
          end do

       end if

    end if

    if (js == 1) then

       if (pd) then

          do j = 0,jsd,-1
          do i = istart,iend+istag

             if  (real(j) <= 1. - q(i,1)/(q(i,2) - q(i,1) + 1.e-12) .and. q(i,1) < q(i,2)) then
                q(i,j) = q(i,j+1)
             else
                q(i,j) = real(2-j)*q(i,1) - real(1-j)*q(i,2)
             end if

          end do
          end do

       else

          do j = 0,jsd,-1
          do i = istart,iend+istag

             q(i,j) = real(2-j)*q(i,1) - real(1-j)*q(i,2)

          end do
          end do

       end if

    end if

    if (ie == npx - 1) then

       if (pd) then

          do j=jstart,jend+jstag
          do i=ie+1+istag,ied+istag

             if (real(i) >= ie+istag + q(ie+istag,j)/(q(ie+istag-1,j)-q(ie+istag,j)+1.e-12) .and. &
                  q(ie+istag,j) < q(ie+istag-1,j)) then
                q(i,j) = q(i-1,j)
             else
                q(i,j) = real(i - (ie+istag-1))*q(ie+istag,j) + real((ie+istag) - i)*q(ie+istag-1,j)
             end if

          end do
          end do

       else

          do j=jstart,jend+jstag
          do i=ie+1+istag,ied+istag

             q(i,j) = real(i - (ie+istag-1))*q(ie+istag,j) + real((ie+istag) - i)*q(ie+istag-1,j)

          end do
          end do

       end if

    end if

    if (je == npy - 1) then

       if (pd) then

          do j=je+1+jstag,jed+jstag
          do i=istart,iend+istag

             if (real(j) >= je+jstag + q(i,je+jstag)/(q(i,je+jstag-1)-q(i,je+jstag)+1.e-12) .and. &
                  q(i,je+jstag-1) > q(i,je+jstag)) then
                q(i,j) = q(i,j-1)
             else
                q(i,j) = real(j - (je+jstag-1))*q(i,je+jstag) + real((je+jstag) - j)*q(i,je+jstag-1)
             end if

          end do
          end do

       else

          do j=je+1+jstag,jed+jstag
          do i=istart,iend+istag

             q(i,j) = real(j - (je+jstag-1))*q(i,je+jstag) + real((je+jstag) - j)*q(i,je+jstag-1)

          end do
          end do

       end if

    end if


    !CORNERS: Average of extrapolations

    if (is == 1 .and. js == 1) then

       if (pd) then

          do j=0,jsd,-1
          do i=0,isd,-1

             if (real(i) <= 1. - q(1,j)/(q(2,j) - q(1,j) + 1.e-12) .and. q(2,j) > q(1,j)) then
                q(i,j) = 0.5*q(i+1,j)
             else
                q(i,j) = 0.5*( real(2-i)*q(1,j) - real(1-i)*q(2,j) )
             end if

             if  (real(j) <= 1. - q(i,1)/(q(i,2) - q(i,1) + 1.e-12) .and. q(i,2) > q(i,1)) then
                q(i,j) = q(i,j) + 0.5*q(i,j+1)

             else
                q(i,j) = q(i,j) + 0.5*(real(2-j)*q(i,1) - real(1-j)*q(i,2))
             end if

          end do
          end do

       else

          do j=jsd,0
          do i=isd,0

             q(i,j) = 0.5*( real(2-i)*q(1,j) - real(1-i)*q(2,j) ) + &
                  0.5*( real(2-j)*q(i,1) - real(1-j)*q(i,2) )

          end do
          end do

       end if

    end if

    if (is == 1 .and. je == npy-1) then

       if (pd) then

          do j=je+1+jstag,jed+jstag
          do i=0,isd,-1

             if (real(i) <= 1. - q(1,j)/(q(2,j) - q(1,j) + 1.e-12) .and. q(2,j) > q(1,j)) then
                q(i,j) = 0.5*q(i+1,j)
             else
                q(i,j) = 0.5*( real(2-i)*q(1,j) - real(1-i)*q(2,j) )
             end if

             !'Unary plus' removed to appease IBM compiler
             !if (real(j) >= je+jstag - q(i,je+jstag)/(q(i,je+jstag-1)-q(i,je+jstag)+1.e-12) .and. &
             if (real(j) >= je+jstag - q(i,je+jstag)/(q(i,je+jstag-1)-q(i,je+jstag)+1.e-12) .and. &
                  q(i,je+jstag-1) > q(i,je+jstag) ) then
                q(i,j) = q(i,j) + 0.5*q(i,j-1)
             else
                q(i,j) = q(i,j) + 0.5*( real(j - (je+jstag-1))*q(i,je+jstag) + real((je+jstag) - j)*q(i,je+jstag-1) )
             end if

          end do
          end do

       else

          do j=je+1+jstag,jed+jstag
          do i=isd,0

             q(i,j) = 0.5*( real(2-i)*q(1,j) - real(1-i)*q(2,j) ) + &
                      0.5*( real(j - (je+jstag-1))*q(i,je+jstag) + real((je+jstag) - j)*q(i,je+jstag-1) )

          end do
          end do

       end if

    end if

    if (ie == npx-1 .and. je == npy-1) then

       if (pd) then

          do j=je+1+jstag,jed+jstag
          do i=ie+1+istag,ied+istag


             if (real(i) >= ie+istag + q(ie+istag,j)/(q(ie+istag-1,j)-q(ie+istag,j)+1.e-12) .and. &
                  q(ie+istag-1,j) > q(ie+istag,j)) then
                q(i,j) = 0.5*q(i-1,j)
             else
                q(i,j) = 0.5*(real(i - (ie+istag-1))*q(ie+istag,j) + real((ie+istag) - i)*q(ie+istag-1,j))
             end if

             if (real(j) >= je+jstag + q(i,je+jstag)/(q(i,je+jstag-1)-q(i,je+jstag)+1.e-12) .and. &
                  q(i,je+jstag-1) > q(i,je+jstag)) then
                q(i,j) = q(i,j) + 0.5*q(i,j-1)
             else
                q(i,j) = q(i,j) + 0.5*( real(j - (je+jstag-1))*q(i,je+jstag) + real((je+jstag) - j)*q(i,je+jstag-1) )
             end if

          end do
          end do

       else

          do j=je+1+jstag,jed+jstag
          do i=ie+1+istag,ied+istag

             q(i,j) = 0.5*( real(i - (ie+istag-1))*q(ie+istag,j) + real((ie+istag) - i)*q(ie+istag-1,j) ) + &
                      0.5*( real(j - (je+jstag-1))*q(i,je+jstag) + real((je+jstag) - j)*q(i,je+jstag-1) )

          end do
          end do

       end if

    end if

    if (ie == npx-1 .and. js == 1) then

       if (pd) then

          do j=0,jsd,-1
          do i=ie+1+istag,ied+istag


             if (real(i) >= ie+istag + q(ie+istag,j)/(q(ie+istag-1,j)-q(ie+istag,j)+1.e-12) .and. &
                  q(ie+istag-1,j) > q(ie+istag,j)) then
                q(i,j) = 0.5*q(i-1,j)
             else
                q(i,j) = 0.5*(real(i - (ie+istag-1))*q(ie+istag,j) + real((ie+istag) - i)*q(ie+istag-1,j))
             end if

             if  (real(j) <= 1. - q(i,1)/(q(i,2) - q(i,1) + 1.e-12) .and. &
                  q(i,2) > q(i,1)) then
                q(i,j) = q(i,j) + 0.5*q(i,j+1)
             else
                q(i,j) = q(i,j) + 0.5*(real(2-j)*q(i,1) - real(1-j)*q(i,2))
             end if

          end do
          end do


       else

          do j=jsd,0
          do i=ie+1+istag,ied+istag

             q(i,j) = 0.5*( real(i - (ie+istag-1))*q(ie+istag,j) + real((ie+istag) - i)*q(ie+istag-1,j) ) + &
                      0.5*( real(2-j)*q(i,1) - real(1-j)*q(i,2) )

          end do
          end do

       end if

    end if


  end subroutine extrapolation_BC

  subroutine fill_nested_grid_2D(var_nest, var_coarse, ind, wt, istag, jstag,  &
      isg, ieg, jsg, jeg, bd, istart_in, iend_in, jstart_in, jend_in)

   type(fv_grid_bounds_type), intent(IN) :: bd
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag), intent(INOUT) :: var_nest
   real, dimension(isg:ieg+istag,jsg:jeg+jstag), intent(IN) :: var_coarse
   integer, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,2), intent(IN) :: ind
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,4), intent(IN) :: wt
   integer, intent(IN) :: istag, jstag, isg, ieg, jsg, jeg
   integer, intent(IN), OPTIONAL :: istart_in, iend_in, jstart_in, jend_in

   integer :: i,j, ic, jc
   integer :: istart, iend, jstart, jend

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

   if (present(istart_in)) then
      istart = istart_in
   else
      istart = isd
   end if
   if (present(iend_in)) then
      iend = iend_in+istag
   else
      iend = ied+istag
   end if

   if (present(jstart_in)) then
      jstart = jstart_in
   else
      jstart = jsd
   end if
   if (present(jend_in)) then
      jend = jend_in+jstag
   else
      jend = jed+jstag
   end if

   do j=jstart,jend
   do i=istart,iend

      ic = ind(i,j,1)
      jc = ind(i,j,2)

      var_nest(i,j) = &
           wt(i,j,1)*var_coarse(ic,  jc) +  &
           wt(i,j,2)*var_coarse(ic,  jc+1) +  &
           wt(i,j,3)*var_coarse(ic+1,jc+1) +  &
           wt(i,j,4)*var_coarse(ic+1,jc)

   end do
   end do

 end subroutine fill_nested_grid_2D

  subroutine fill_nested_grid_3D(var_nest, var_coarse, ind, wt, istag, jstag,  &
      isg, ieg, jsg, jeg, npz, bd, istart_in, iend_in, jstart_in, jend_in)

   type(fv_grid_bounds_type), intent(IN) :: bd
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,npz), intent(INOUT) :: var_nest
   real, dimension(isg:ieg+istag,jsg:jeg+jstag,npz), intent(IN) :: var_coarse
   integer, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,2), intent(IN) :: ind
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,4), intent(IN) :: wt
   integer, intent(IN) :: istag, jstag, isg, ieg, jsg, jeg, npz
   integer, intent(IN), OPTIONAL :: istart_in, iend_in, jstart_in, jend_in

   integer :: i,j, ic, jc, k
   integer :: istart, iend, jstart, jend

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

   if (present(istart_in)) then
      istart = istart_in
   else
      istart = isd
   end if
   if (present(iend_in)) then
      iend = iend_in+istag
   else
      iend = ied+istag
   end if

   if (present(jstart_in)) then
      jstart = jstart_in
   else
      jstart = jsd
   end if
   if (present(jend_in)) then
      jend = jend_in+jstag
   else
      jend = jed+jstag
   end if

   do k=1,npz

   do j=jstart,jend
   do i=istart,iend

      ic = ind(i,j,1)
      jc = ind(i,j,2)

      var_nest(i,j,k) = &
           wt(i,j,1)*var_coarse(ic,  jc,  k) +  &
           wt(i,j,2)*var_coarse(ic,  jc+1,k) +  &
           wt(i,j,3)*var_coarse(ic+1,jc+1,k) +  &
           wt(i,j,4)*var_coarse(ic+1,jc,  k)

   end do
   end do

   end do

 end subroutine fill_nested_grid_3D

!!$ subroutine nested_grid_BC_mpp_2d(var_nest, nest_domain, ind, wt, istag, jstag, &
!!$      npx, npy, bd, isg, ieg, jsg, jeg, nstep_in, nsplit_in, proc_in)
!!$
!!$   type(fv_grid_bounds_type), intent(IN) :: bd
!!$   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag), intent(INOUT) :: var_nest
!!$   real, dimension(isg:ieg+istag,jsg:jeg+jstag), intent(IN) :: var_coarse
!!$   type(nest_domain_type), intent(INOUT) :: nest_domain
!!$   integer, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,2), intent(IN) :: ind
!!$   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,4), intent(IN) :: wt
!!$   integer, intent(IN) :: istag, jstag, npx, npy, isg, ieg, jsg, jeg
!!$   integer, intent(IN), OPTIONAL :: nstep_in, nsplit_in
!!$   logical, intent(IN), OPTIONAL :: proc_in
!!$
!!$   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,1) :: var_nest_3d
!!$
!!$   integer :: i,j
!!$
!!$   do j=bd%jsd,bd%jed+jstag
!!$   do i=bd%isd,bd%ied+istag
!!$      var_nest_3d(i,j,1) = var_nest(i,j)
!!$   enddo
!!$   enddo
!!$
!!$   call nested_grid_BC_mpp_3d(var_nest_3d, nest_domain, ind, wt, istag, jstag, &
!!$      npx, npy, 1, bd, isg, ieg, jsg, jeg, nstep_in, nsplit_in, proc_in)
!!$
!!$
!!$ end subroutine nested_grid_BC_mpp_2d

 subroutine nested_grid_BC_mpp_3d(var_nest, var_coarse, nest_domain, ind, wt, istag, jstag, &
      npx, npy, npz, bd, isg, ieg, jsg, jeg, nest_level, nstep_in, nsplit_in, proc_in)

   type(fv_grid_bounds_type), intent(IN) :: bd
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,npz), intent(INOUT) :: var_nest
   real, dimension(isg:ieg+istag,jsg:jeg+jstag,npz), intent(IN) :: var_coarse
   type(nest_domain_type), intent(INOUT) :: nest_domain
   integer, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,2), intent(IN) :: ind
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,4), intent(IN) :: wt
   integer, intent(IN) :: istag, jstag, npx, npy, npz, isg, ieg, jsg, jeg
   integer, intent(IN) :: nest_level
   integer, intent(IN), OPTIONAL :: nstep_in, nsplit_in
   logical, intent(IN), OPTIONAL :: proc_in

   integer                      :: isw_f, iew_f, jsw_f, jew_f, isw_c, iew_c, jsw_c, jew_c
   integer                      :: ise_f, iee_f, jse_f, jee_f, ise_c, iee_c, jse_c, jee_c
   integer                      :: iss_f, ies_f, jss_f, jes_f, iss_c, ies_c, jss_c, jes_c
   integer                      :: isn_f, ien_f, jsn_f, jen_f, isn_c, ien_c, jsn_c, jen_c
   real,    allocatable         :: wbuffer(:,:,:)
   real,    allocatable         :: ebuffer(:,:,:)
   real,    allocatable         :: sbuffer(:,:,:)
   real,    allocatable         :: nbuffer(:,:,:)

   integer :: i,j, ic, jc, istart, iend, k

   integer :: position
   logical :: process

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

   if (PRESENT(proc_in)) then
      process = proc_in
   else
      process = .true.
   endif

   if (istag == 1 .and. jstag == 1) then
      position = CORNER
   else if (istag == 0 .and. jstag == 1) then
      position = NORTH
   else if (istag == 1 .and. jstag == 0) then
      position = EAST
   else
      position = CENTER
   end if

   call mpp_get_C2F_index(nest_domain, isw_f, iew_f, jsw_f, jew_f, isw_c, iew_c, jsw_c, jew_c, &
        WEST, nest_level=nest_level, position=position)
   call mpp_get_C2F_index(nest_domain, ise_f, iee_f, jse_f, jee_f, ise_c, iee_c, jse_c, jee_c, &
        EAST, nest_level=nest_level, position=position)
   call mpp_get_C2F_index(nest_domain, iss_f, ies_f, jss_f, jes_f, iss_c, ies_c, jss_c, jes_c, &
        SOUTH, nest_level=nest_level, position=position)
   call mpp_get_C2F_index(nest_domain, isn_f, ien_f, jsn_f, jen_f, isn_c, ien_c, jsn_c, jen_c, &
        NORTH, nest_level=nest_level, position=position)

   if( iew_c .GE. isw_c .AND. jew_c .GE. jsw_c ) then
      allocate(wbuffer(isw_c:iew_c, jsw_c:jew_c,npz))
   else
      allocate(wbuffer(1,1,1))
   endif
   wbuffer = 0

   if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
      allocate(ebuffer(ise_c:iee_c, jse_c:jee_c,npz))
   else
      allocate(ebuffer(1,1,1))
   endif
   ebuffer = 0

   if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
      allocate(sbuffer(iss_c:ies_c, jss_c:jes_c,npz))
   else
      allocate(sbuffer(1,1,1))
   endif
   sbuffer = 0

   if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
      allocate(nbuffer(isn_c:ien_c, jsn_c:jen_c,npz))
   else
      allocate(nbuffer(1,1,1))
   endif
   nbuffer = 0


       call timing_on ('COMM_TOTAL')
   call mpp_update_nest_fine(var_coarse, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, &
        nest_level=nest_level, position=position)
       call timing_off('COMM_TOTAL')

   if (process) then

   if (is == 1) then
!OMP parallel do default(none) shared(npz,jsd,jed,jstag,isd,ind,var_nest,wt,wbuffer) private(ic,jc)
      do k=1,npz
      do j=jsd,jed+jstag
         do i=isd,0

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            var_nest(i,j,k) = &
                 wt(i,j,1)*wbuffer(ic,  jc,  k) +  &
                 wt(i,j,2)*wbuffer(ic,  jc+1,k) +  &
                 wt(i,j,3)*wbuffer(ic+1,jc+1,k) +  &
                 wt(i,j,4)*wbuffer(ic+1,jc,  k)

         end do
      end do
      end do
   end if

   if (js == 1) then

      if (is == 1) then
         istart = is
      else
         istart = isd
      end if

      if (ie == npx-1) then
         iend = ie
      else
         iend = ied
      end if

!OMP parallel do default(none) shared(npz,jsd,istart,iend,istag,ind,var_nest,wt,sbuffer) private(ic,jc)
      do k=1,npz
      do j=jsd,0
         do i=istart,iend+istag

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            var_nest(i,j,k) = &
                 wt(i,j,1)*sbuffer(ic,  jc,  k) +  &
                 wt(i,j,2)*sbuffer(ic,  jc+1,k) +  &
                 wt(i,j,3)*sbuffer(ic+1,jc+1,k) +  &
                 wt(i,j,4)*sbuffer(ic+1,jc,  k)

         end do
      end do
      end do
   end if


   if (ie == npx-1) then
!OMP parallel do default(none) shared(npz,jsd,jed,jstag,npx,ied,istag,ind,var_nest,wt,ebuffer) private(ic,jc)
      do k=1,npz
      do j=jsd,jed+jstag
         do i=npx+istag,ied+istag

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            var_nest(i,j,k) = &
                 wt(i,j,1)*ebuffer(ic,  jc,  k) +  &
                 wt(i,j,2)*ebuffer(ic,  jc+1,k) +  &
                 wt(i,j,3)*ebuffer(ic+1,jc+1,k) +  &
                 wt(i,j,4)*ebuffer(ic+1,jc,  k)

         end do
      end do
      end do
   end if

   if (je == npy-1) then

      if (is == 1) then
         istart = is
      else
         istart = isd
      end if

      if (ie == npx-1) then
         iend = ie
      else
         iend = ied
      end if

!OMP parallel do default(none) shared(npz,jstag,npy,jed,istart,iend,istag,ind,var_nest,wt,nbuffer) private(ic,jc)
      do k=1,npz
      do j=npy+jstag,jed+jstag
         do i=istart,iend+istag

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            var_nest(i,j,k) = &
                 wt(i,j,1)*nbuffer(ic,  jc,  k) +  &
                 wt(i,j,2)*nbuffer(ic,  jc+1,k) +  &
                 wt(i,j,3)*nbuffer(ic+1,jc+1,k) +  &
                 wt(i,j,4)*nbuffer(ic+1,jc,  k)

         end do
      end do
      end do
   end if

   endif !process

   deallocate(wbuffer, ebuffer, sbuffer, nbuffer)

 end subroutine nested_grid_BC_mpp_3d

 subroutine get_vector_position(position_x, position_y, gridtype)
   integer,          intent(OUT) :: position_x, position_y
   integer, optional, intent(IN) :: gridtype

   integer :: grid_offset_type

   grid_offset_type = AGRID
   if(present(gridtype)) grid_offset_type = gridtype

   select case(grid_offset_type)
   case (AGRID)
      position_x = CENTER
      position_y = CENTER
   case (BGRID_NE)
      position_x = CORNER
      position_y = CORNER
   case (CGRID_NE)
      position_x = EAST
      position_y = NORTH
   case (DGRID_NE)
      position_y = EAST
      position_x = NORTH
   case default
      call mpp_error(FATAL, "get_vector_position: invalid value of gridtype")
   end select


 end subroutine get_vector_position

 subroutine init_buffer(nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, npz, nest_level, position)
   type(nest_domain_type),            intent(INOUT) :: nest_domain
   real, allocatable, dimension(:,:,:), intent(OUT) :: wbuffer, sbuffer, ebuffer, nbuffer
   integer,                             intent(IN)  :: npz, position, nest_level
   integer :: isw_f, iew_f, jsw_f, jew_f, isw_c, iew_c, jsw_c, jew_c
   integer :: ise_f, iee_f, jse_f, jee_f, ise_c, iee_c, jse_c, jee_c
   integer :: iss_f, ies_f, jss_f, jes_f, iss_c, ies_c, jss_c, jes_c
   integer :: isn_f, ien_f, jsn_f, jen_f, isn_c, ien_c, jsn_c, jen_c

   call mpp_get_C2F_index(nest_domain, isw_f, iew_f, jsw_f, jew_f, isw_c, iew_c, jsw_c, jew_c, &
        WEST, nest_level=nest_level, position=position)
   call mpp_get_C2F_index(nest_domain, ise_f, iee_f, jse_f, jee_f, ise_c, iee_c, jse_c, jee_c, &
        EAST, nest_level=nest_level, position=position)
   call mpp_get_C2F_index(nest_domain, iss_f, ies_f, jss_f, jes_f, iss_c, ies_c, jss_c, jes_c, &
        SOUTH, nest_level=nest_level, position=position)
   call mpp_get_C2F_index(nest_domain, isn_f, ien_f, jsn_f, jen_f, isn_c, ien_c, jsn_c, jen_c, &
        NORTH, nest_level=nest_level, position=position)

   if( iew_c .GE. isw_c .AND. jew_c .GE. jsw_c ) then
      allocate(wbuffer(isw_c:iew_c, jsw_c:jew_c,npz))
   else
      allocate(wbuffer(1,1,1))
   endif
   wbuffer = 0

   if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
      allocate(ebuffer(ise_c:iee_c, jse_c:jee_c,npz))
   else
      allocate(ebuffer(1,1,1))
   endif
   ebuffer = 0

   if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
      allocate(sbuffer(iss_c:ies_c, jss_c:jes_c,npz))
   else
      allocate(sbuffer(1,1,1))
   endif
   sbuffer = 0

   if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
      allocate(nbuffer(isn_c:ien_c, jsn_c:jen_c,npz))
   else
      allocate(nbuffer(1,1,1))
   endif
   nbuffer = 0

 end subroutine init_buffer


 subroutine nested_grid_BC_mpp_3d_vector(u_nest, v_nest, u_coarse, v_coarse, nest_domain, ind_u, ind_v, wt_u, wt_v, &
      istag_u, jstag_u, istag_v, jstag_v, npx, npy, npz, bd, isg, ieg, jsg, jeg, nest_level, nstep_in, nsplit_in, proc_in, &
      flags, gridtype)

   type(fv_grid_bounds_type), intent(IN) :: bd
   integer, intent(IN) :: istag_u, jstag_u, istag_v, jstag_v, npx, npy, npz, isg, ieg, jsg, jeg
   real, dimension(bd%isd:bd%ied+istag_u,bd%jsd:bd%jed+jstag_u,npz), intent(INOUT) :: u_nest
   real, dimension(bd%isd:bd%ied+istag_v,bd%jsd:bd%jed+jstag_v,npz), intent(INOUT) :: v_nest
   real, dimension(isg:ieg+istag_u,jsg:jeg+jstag_u,npz), intent(IN) :: u_coarse
   real, dimension(isg:ieg+istag_v,jsg:jeg+jstag_v,npz), intent(IN) :: v_coarse
   type(nest_domain_type), intent(INOUT) :: nest_domain
   integer, dimension(bd%isd:bd%ied+istag_u,bd%jsd:bd%jed+jstag_u,2), intent(IN) :: ind_u
   integer, dimension(bd%isd:bd%ied+istag_v,bd%jsd:bd%jed+jstag_v,2), intent(IN) :: ind_v
   real, dimension(bd%isd:bd%ied+istag_u,bd%jsd:bd%jed+jstag_u,4), intent(IN) :: wt_u
   real, dimension(bd%isd:bd%ied+istag_v,bd%jsd:bd%jed+jstag_v,4), intent(IN) :: wt_v
   integer, intent(IN)           :: nest_level
   integer, intent(IN), OPTIONAL :: nstep_in, nsplit_in
   logical, intent(IN), OPTIONAL :: proc_in
   integer, intent(IN), OPTIONAL :: flags, gridtype

   real,    allocatable         :: wbufferx(:,:,:), wbuffery(:,:,:)
   real,    allocatable         :: ebufferx(:,:,:), ebuffery(:,:,:)
   real,    allocatable         :: sbufferx(:,:,:), sbuffery(:,:,:)
   real,    allocatable         :: nbufferx(:,:,:), nbuffery(:,:,:)

   integer :: i,j, ic, jc, istart, iend, k

   integer :: position_x, position_y
   logical :: process

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

   if (PRESENT(proc_in)) then
      process = proc_in
   else
      process = .true.
   endif

   call get_vector_position(position_x, position_y, gridtype)
   call init_buffer(nest_domain, wbufferx, sbufferx, ebufferx, nbufferx, npz, nest_level, position_x)
   call init_buffer(nest_domain, wbuffery, sbuffery, ebuffery, nbuffery, npz, nest_level, position_x)

       call timing_on ('COMM_TOTAL')
   call mpp_update_nest_fine(u_coarse, v_coarse, nest_domain, wbufferx, wbuffery, sbufferx, sbuffery, &
                             ebufferx, ebuffery, nbufferx, nbuffery, flags=flags, nest_level=nest_level, gridtype=gridtype)
       call timing_off('COMM_TOTAL')

   if (process) then

   if (is == 1) then
!OMP parallel do default(none) shared(npz,jsd,jed,jstag,isd,ind,var_nest,wt,wbuffer) private(ic,jc)
      do k=1,npz
      do j=jsd,jed+jstag_u
         do i=isd,0

            ic = ind_u(i,j,1)
            jc = ind_u(i,j,2)

            u_nest(i,j,k) = &
                 wt_u(i,j,1)*wbufferx(ic,  jc,  k) +  &
                 wt_u(i,j,2)*wbufferx(ic,  jc+1,k) +  &
                 wt_u(i,j,3)*wbufferx(ic+1,jc+1,k) +  &
                 wt_u(i,j,4)*wbufferx(ic+1,jc,  k)

         end do
      end do
      do j=jsd,jed+jstag_v
         do i=isd,0

            ic = ind_v(i,j,1)
            jc = ind_v(i,j,2)

            v_nest(i,j,k) = &
                 wt_v(i,j,1)*wbuffery(ic,  jc,  k) +  &
                 wt_v(i,j,2)*wbuffery(ic,  jc+1,k) +  &
                 wt_v(i,j,3)*wbuffery(ic+1,jc+1,k) +  &
                 wt_v(i,j,4)*wbuffery(ic+1,jc,  k)

         end do
      end do
      end do

   end if

   if (js == 1) then

      if (is == 1) then
         istart = is
      else
         istart = isd
      end if

      if (ie == npx-1) then
         iend = ie
      else
         iend = ied
      end if

!OMP parallel do default(none) shared(npz,jsd,istart,iend,istag,ind,var_nest,wt,sbuffer) private(ic,jc)
      do k=1,npz
      do j=jsd,0
         do i=istart,iend+istag_u

            ic = ind_u(i,j,1)
            jc = ind_u(i,j,2)

            u_nest(i,j,k) = &
                 wt_u(i,j,1)*sbufferx(ic,  jc,  k) +  &
                 wt_u(i,j,2)*sbufferx(ic,  jc+1,k) +  &
                 wt_u(i,j,3)*sbufferx(ic+1,jc+1,k) +  &
                 wt_u(i,j,4)*sbufferx(ic+1,jc,  k)

         end do
      end do
      do j=jsd,0
         do i=istart,iend+istag_v

            ic = ind_v(i,j,1)
            jc = ind_v(i,j,2)

            v_nest(i,j,k) = &
                 wt_v(i,j,1)*sbuffery(ic,  jc,  k) +  &
                 wt_v(i,j,2)*sbuffery(ic,  jc+1,k) +  &
                 wt_v(i,j,3)*sbuffery(ic+1,jc+1,k) +  &
                 wt_v(i,j,4)*sbuffery(ic+1,jc,  k)

         end do
      end do
      end do
   end if


   if (ie == npx-1) then
!OMP parallel do default(none) shared(npz,jsd,jed,jstag,npx,ied,istag,ind,var_nest,wt,ebuffer) private(ic,jc)
      do k=1,npz
      do j=jsd,jed+jstag_u
         do i=npx+istag_u,ied+istag_u

            ic = ind_u(i,j,1)
            jc = ind_u(i,j,2)

            u_nest(i,j,k) = &
                 wt_u(i,j,1)*ebufferx(ic,  jc,  k) +  &
                 wt_u(i,j,2)*ebufferx(ic,  jc+1,k) +  &
                 wt_u(i,j,3)*ebufferx(ic+1,jc+1,k) +  &
                 wt_u(i,j,4)*ebufferx(ic+1,jc,  k)

         end do
      end do
      do j=jsd,jed+jstag_v
         do i=npx+istag_v,ied+istag_v

            ic = ind_v(i,j,1)
            jc = ind_v(i,j,2)

            v_nest(i,j,k) = &
                 wt_v(i,j,1)*ebuffery(ic,  jc,  k) +  &
                 wt_v(i,j,2)*ebuffery(ic,  jc+1,k) +  &
                 wt_v(i,j,3)*ebuffery(ic+1,jc+1,k) +  &
                 wt_v(i,j,4)*ebuffery(ic+1,jc,  k)

         end do
      end do
      end do
   end if

   if (je == npy-1) then

      if (is == 1) then
         istart = is
      else
         istart = isd
      end if

      if (ie == npx-1) then
         iend = ie
      else
         iend = ied
      end if

!OMP parallel do default(none) shared(npz,jstag,npy,jed,istart,iend,istag,ind,var_nest,wt,nbuffer) private(ic,jc)
      do k=1,npz
      do j=npy+jstag_u,jed+jstag_u
         do i=istart,iend+istag_u

            ic = ind_u(i,j,1)
            jc = ind_u(i,j,2)

            u_nest(i,j,k) = &
                 wt_u(i,j,1)*nbufferx(ic,  jc,  k) +  &
                 wt_u(i,j,2)*nbufferx(ic,  jc+1,k) +  &
                 wt_u(i,j,3)*nbufferx(ic+1,jc+1,k) +  &
                 wt_u(i,j,4)*nbufferx(ic+1,jc,  k)

         end do
      end do
      do j=npy+jstag_v,jed+jstag_v
         do i=istart,iend+istag_v

            ic = ind_v(i,j,1)
            jc = ind_v(i,j,2)

            v_nest(i,j,k) = &
                 wt_v(i,j,1)*nbuffery(ic,  jc,  k) +  &
                 wt_v(i,j,2)*nbuffery(ic,  jc+1,k) +  &
                 wt_v(i,j,3)*nbuffery(ic+1,jc+1,k) +  &
                 wt_v(i,j,4)*nbuffery(ic+1,jc,  k)

         end do
      end do
      end do
   end if

   endif !process

   deallocate(wbufferx, ebufferx, sbufferx, nbufferx)
   deallocate(wbuffery, ebuffery, sbuffery, nbuffery)

 end subroutine nested_grid_BC_mpp_3d_vector


 subroutine nested_grid_BC_mpp_send_3d(var_coarse, nest_domain, istag, jstag, nest_level)

   real, dimension(:,:,:), intent(IN) :: var_coarse
   type(nest_domain_type), intent(INOUT) :: nest_domain
   integer, intent(IN) :: istag, jstag
   integer, intent(IN) :: nest_level

   real,    allocatable         :: wbuffer(:,:,:)
   real,    allocatable         :: ebuffer(:,:,:)
   real,    allocatable         :: sbuffer(:,:,:)
   real,    allocatable         :: nbuffer(:,:,:)

   integer :: i,j, ic, jc, istart, iend, k

   integer :: position


   if (istag == 1 .and. jstag == 1) then
      position = CORNER
   else if (istag == 0 .and. jstag == 1) then
      position = NORTH
   else if (istag == 1 .and. jstag == 0) then
      position = EAST
   else
      position = CENTER
   end if


      allocate(wbuffer(1,1,1))

      allocate(ebuffer(1,1,1))

      allocate(sbuffer(1,1,1))

      allocate(nbuffer(1,1,1))


       call timing_on ('COMM_TOTAL')
   call mpp_update_nest_fine(var_coarse, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level=nest_level, position=position)
       call timing_off('COMM_TOTAL')


   deallocate(wbuffer, ebuffer, sbuffer, nbuffer)

 end subroutine nested_grid_BC_mpp_send_3d

 subroutine nested_grid_BC_mpp_send_2d(var_coarse, nest_domain, istag, jstag, nest_level)

   real, dimension(:,:), intent(IN) :: var_coarse
   type(nest_domain_type), intent(INOUT) :: nest_domain
   integer, intent(IN) :: istag, jstag
   integer, intent(IN) :: nest_level

   real,    allocatable         :: wbuffer(:,:)
   real,    allocatable         :: ebuffer(:,:)
   real,    allocatable         :: sbuffer(:,:)
   real,    allocatable         :: nbuffer(:,:)

   integer :: i,j, ic, jc, istart, iend, k

   integer :: position


   if (istag == 1 .and. jstag == 1) then
      position = CORNER
   else if (istag == 0 .and. jstag == 1) then
      position = NORTH
   else if (istag == 1 .and. jstag == 0) then
      position = EAST
   else
      position = CENTER
   end if


   allocate(wbuffer(1,1))

   allocate(ebuffer(1,1))

   allocate(sbuffer(1,1))

   allocate(nbuffer(1,1))


       call timing_on ('COMM_TOTAL')
   call mpp_update_nest_fine(var_coarse, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level=nest_level, position=position)
       call timing_off('COMM_TOTAL')


   deallocate(wbuffer, ebuffer, sbuffer, nbuffer)

 end subroutine nested_grid_BC_mpp_send_2d

 subroutine nested_grid_BC_2D_mpp(var_nest, var_coarse, nest_domain, ind, wt, istag, jstag, &
      npx, npy, bd, isg, ieg, jsg, jeg, nest_level, nstep_in, nsplit_in, proc_in)

   type(fv_grid_bounds_type), intent(IN) :: bd
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag), intent(INOUT) :: var_nest
   real, dimension(isg:ieg+istag,jsg:jeg+jstag), intent(IN) :: var_coarse
   type(nest_domain_type), intent(INOUT) :: nest_domain
   integer, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,2), intent(IN) :: ind
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,4), intent(IN) :: wt
   integer, intent(IN) :: istag, jstag, npx, npy, isg, ieg, jsg, jeg
   integer, intent(IN), OPTIONAL :: nest_level
   integer, intent(IN), OPTIONAL :: nstep_in, nsplit_in
   logical, intent(IN), OPTIONAL :: proc_in

   integer                      :: isw_f, iew_f, jsw_f, jew_f, isw_c, iew_c, jsw_c, jew_c
   integer                      :: ise_f, iee_f, jse_f, jee_f, ise_c, iee_c, jse_c, jee_c
   integer                      :: iss_f, ies_f, jss_f, jes_f, iss_c, ies_c, jss_c, jes_c
   integer                      :: isn_f, ien_f, jsn_f, jen_f, isn_c, ien_c, jsn_c, jen_c
   real,    allocatable         :: wbuffer(:,:)
   real,    allocatable         :: ebuffer(:,:)
   real,    allocatable         :: sbuffer(:,:)
   real,    allocatable         :: nbuffer(:,:)

   integer :: i,j, ic, jc, istart, iend, k
   integer :: nl = 1 !nest_level

   integer :: position
   logical :: process

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

   if (PRESENT(proc_in)) then
      process = proc_in
   else
      process = .true.
   endif

   if (PRESENT(nest_level)) then
      nl = nest_level
   endif

   if (istag == 1 .and. jstag == 1) then
      position = CORNER
   else if (istag == 0 .and. jstag == 1) then
      position = NORTH
   else if (istag == 1 .and. jstag == 0) then
      position = EAST
   else
      position = CENTER
   end if

   call mpp_get_C2F_index(nest_domain, isw_f, iew_f, jsw_f, jew_f, isw_c, iew_c, jsw_c, jew_c, &
        WEST, nest_level=nl, position=position)
   call mpp_get_C2F_index(nest_domain, ise_f, iee_f, jse_f, jee_f, ise_c, iee_c, jse_c, jee_c, &
        EAST, nest_level=nl, position=position)
   call mpp_get_C2F_index(nest_domain, iss_f, ies_f, jss_f, jes_f, iss_c, ies_c, jss_c, jes_c, &
        SOUTH, nest_level=nl, position=position)
   call mpp_get_C2F_index(nest_domain, isn_f, ien_f, jsn_f, jen_f, isn_c, ien_c, jsn_c, jen_c, &
        NORTH, nest_level=nl, position=position)

   if( iew_c .GE. isw_c .AND. jew_c .GE. jsw_c ) then
      allocate(wbuffer(isw_c:iew_c, jsw_c:jew_c))
   else
      allocate(wbuffer(1,1))
   endif
   wbuffer = 0

   if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
      allocate(ebuffer(ise_c:iee_c, jse_c:jee_c))
   else
      allocate(ebuffer(1,1))
   endif
   ebuffer = 0

   if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
      allocate(sbuffer(iss_c:ies_c, jss_c:jes_c))
   else
      allocate(sbuffer(1,1))
   endif
   sbuffer = 0

   if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
      allocate(nbuffer(isn_c:ien_c, jsn_c:jen_c))
   else
      allocate(nbuffer(1,1))
   endif
   nbuffer = 0

       call timing_on ('COMM_TOTAL')
   call mpp_update_nest_fine(var_coarse, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level=nl, position=position)
       call timing_off('COMM_TOTAL')

   if (process) then

   if (is == 1) then
      do j=jsd,jed+jstag
         do i=isd,0

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            var_nest(i,j) = &
                 wt(i,j,1)*wbuffer(ic,  jc) +  &
                 wt(i,j,2)*wbuffer(ic,  jc+1) +  &
                 wt(i,j,3)*wbuffer(ic+1,jc+1) +  &
                 wt(i,j,4)*wbuffer(ic+1,jc)

         end do
      end do
   end if

   if (js == 1) then

      if (is == 1) then
         istart = is
      else
         istart = isd
      end if

      if (ie == npx-1) then
         iend = ie
      else
         iend = ied
      end if

      do j=jsd,0
         do i=istart,iend+istag

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            var_nest(i,j) = &
                 wt(i,j,1)*sbuffer(ic,  jc) +  &
                 wt(i,j,2)*sbuffer(ic,  jc+1) +  &
                 wt(i,j,3)*sbuffer(ic+1,jc+1) +  &
                 wt(i,j,4)*sbuffer(ic+1,jc)

         end do
      end do
   end if


   if (ie == npx-1) then
      do j=jsd,jed+jstag
         do i=npx+istag,ied+istag

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            var_nest(i,j) = &
                 wt(i,j,1)*ebuffer(ic,  jc) +  &
                 wt(i,j,2)*ebuffer(ic,  jc+1) +  &
                 wt(i,j,3)*ebuffer(ic+1,jc+1) +  &
                 wt(i,j,4)*ebuffer(ic+1,jc)

         end do
      end do
   end if

   if (je == npy-1) then

      if (is == 1) then
         istart = is
      else
         istart = isd
      end if

      if (ie == npx-1) then
         iend = ie
      else
         iend = ied
      end if

      do j=npy+jstag,jed+jstag
         do i=istart,iend+istag

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            var_nest(i,j) = &
                 wt(i,j,1)*nbuffer(ic,  jc) +  &
                 wt(i,j,2)*nbuffer(ic,  jc+1) +  &
                 wt(i,j,3)*nbuffer(ic+1,jc+1) +  &
                 wt(i,j,4)*nbuffer(ic+1,jc)

         end do
      end do
   end if

   endif !process

   deallocate(wbuffer, ebuffer, sbuffer, nbuffer)

 end subroutine nested_grid_BC_2D_mpp

 subroutine nested_grid_BC_2D(var_nest, var_coarse, ind, wt, istag, jstag, &
      npx, npy, bd, isg, ieg, jsg, jeg, nstep_in, nsplit_in)

   type(fv_grid_bounds_type), intent(IN) :: bd
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag), intent(INOUT) :: var_nest
   real, dimension(isg:ieg+istag,jsg:jeg+jstag), intent(IN) :: var_coarse
   integer, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,2), intent(IN) :: ind
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,4), intent(IN) :: wt
   integer, intent(IN) :: istag, jstag, npx, npy, isg, ieg, jsg, jeg
   integer, intent(IN), OPTIONAL :: nstep_in, nsplit_in

   integer :: nstep, nsplit

   integer :: i,j, ic, jc, istart, iend

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

   if ( .not. present(nstep_in) .or. .not. present(nsplit_in) ) then
      nstep = 1
      nsplit = 2
   else
      nstep = nstep_in
      nsplit = nsplit_in
   end if

   if (is == 1) then
      do j=jsd,jed+jstag
         do i=isd,0

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            var_nest(i,j) = &
                 wt(i,j,1)*var_coarse(ic,  jc) +  &
                 wt(i,j,2)*var_coarse(ic,  jc+1) +  &
                 wt(i,j,3)*var_coarse(ic+1,jc+1) +  &
                 wt(i,j,4)*var_coarse(ic+1,jc)

         end do
      end do
   end if

   if (js == 1) then

      if (is == 1) then
         istart = is
      else
         istart = isd
      end if

      if (ie == npx-1) then
         iend = ie
      else
         iend = ied
      end if

      do j=jsd,0
         do i=istart,iend+istag

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            var_nest(i,j) = &
                 wt(i,j,1)*var_coarse(ic,  jc) +  &
                 wt(i,j,2)*var_coarse(ic,  jc+1) +  &
                 wt(i,j,3)*var_coarse(ic+1,jc+1) +  &
                 wt(i,j,4)*var_coarse(ic+1,jc)

         end do
      end do
   end if


   if (ie == npx-1) then
      do j=jsd,jed+jstag
         do i=npx+istag,ied+istag

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            var_nest(i,j) = &
                 wt(i,j,1)*var_coarse(ic,  jc) +  &
                 wt(i,j,2)*var_coarse(ic,  jc+1) +  &
                 wt(i,j,3)*var_coarse(ic+1,jc+1) +  &
                 wt(i,j,4)*var_coarse(ic+1,jc)

         end do
      end do
   end if

   if (je == npy-1) then

      if (is == 1) then
         istart = is
      else
         istart = isd
      end if

      if (ie == npx-1) then
         iend = ie
      else
         iend = ied
      end if


      do j=npy+jstag,jed+jstag
         do i=istart,iend+istag

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            var_nest(i,j) = &
                 wt(i,j,1)*var_coarse(ic,  jc) +  &
                 wt(i,j,2)*var_coarse(ic,  jc+1) +  &
                 wt(i,j,3)*var_coarse(ic+1,jc+1) +  &
                 wt(i,j,4)*var_coarse(ic+1,jc)

         end do
      end do
   end if



 end subroutine nested_grid_BC_2D

 subroutine nested_grid_BC_3D(var_nest, var_coarse, ind, wt, istag, jstag, &
      npx, npy, npz, bd, isg, ieg, jsg, jeg, nstep_in, nsplit_in)

   type(fv_grid_bounds_type), intent(IN) :: bd
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,npz), intent(INOUT) :: var_nest
   real, dimension(isg:ieg+istag,jsg:jeg+jstag,npz), intent(IN) :: var_coarse
   integer, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,2), intent(IN) :: ind
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,4), intent(IN) :: wt
   integer, intent(IN) :: istag, jstag, npx, npy, isg, ieg, jsg, jeg, npz
   integer, intent(IN), OPTIONAL :: nstep_in, nsplit_in

   integer :: nstep, nsplit

   integer :: i,j, ic, jc, istart, iend, k

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

   if ( .not. present(nstep_in) .or. .not. present(nsplit_in) ) then
      nstep = 1
      nsplit = 2
   else
      nstep = nstep_in
      nsplit = nsplit_in
   end if

   if (is == 1) then
!OMP parallel do default(none) shared(npz,jsd,jed,jstag,isd,ind,var_nest,wt,var_coarse) private(ic,jc)
      do k=1,npz
      do j=jsd,jed+jstag
         do i=isd,0

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            var_nest(i,j,k) = &
                 wt(i,j,1)*var_coarse(ic,  jc,  k) +  &
                 wt(i,j,2)*var_coarse(ic,  jc+1,k) +  &
                 wt(i,j,3)*var_coarse(ic+1,jc+1,k) +  &
                 wt(i,j,4)*var_coarse(ic+1,jc,  k)

         end do
      end do
      end do
   end if

   if (js == 1) then

      if (is == 1) then
         istart = is
      else
         istart = isd
      end if

      if (ie == npx-1) then
         iend = ie
      else
         iend = ied
      end if

!OMP parallel do default(none) shared(npz,jsd,istart,iend,istag,ind,var_nest,wt,var_coarse) private(ic,jc)
      do k=1,npz
      do j=jsd,0
         do i=istart,iend+istag

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            var_nest(i,j,k) = &
                 wt(i,j,1)*var_coarse(ic,  jc,  k) +  &
                 wt(i,j,2)*var_coarse(ic,  jc+1,k) +  &
                 wt(i,j,3)*var_coarse(ic+1,jc+1,k) +  &
                 wt(i,j,4)*var_coarse(ic+1,jc,  k)

         end do
      end do
      end do
   end if


   if (ie == npx-1) then
!OMP parallel do default(none) shared(npz,jsd,jed,jstag,npx,ied,istag,ind,var_nest,wt,var_coarse) private(ic,jc)
      do k=1,npz
      do j=jsd,jed+jstag
         do i=npx+istag,ied+istag

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            var_nest(i,j,k) = &
                 wt(i,j,1)*var_coarse(ic,  jc,  k) +  &
                 wt(i,j,2)*var_coarse(ic,  jc+1,k) +  &
                 wt(i,j,3)*var_coarse(ic+1,jc+1,k) +  &
                 wt(i,j,4)*var_coarse(ic+1,jc,  k)

         end do
      end do
      end do
   end if

   if (je == npy-1) then

      if (is == 1) then
         istart = is
      else
         istart = isd
      end if

      if (ie == npx-1) then
         iend = ie
      else
         iend = ied
      end if

!OMP parallel do default(none) shared(npz,npy,jed,jstag,istart,iend,istag,ind,var_nest,wt,var_coarse) private(ic,jc)
      do k=1,npz
      do j=npy+jstag,jed+jstag
         do i=istart,iend+istag

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            var_nest(i,j,k) = &
                 wt(i,j,1)*var_coarse(ic,  jc,  k) +  &
                 wt(i,j,2)*var_coarse(ic,  jc+1,k) +  &
                 wt(i,j,3)*var_coarse(ic+1,jc+1,k) +  &
                 wt(i,j,4)*var_coarse(ic+1,jc,  k)

         end do
      end do
      end do
   end if



 end subroutine nested_grid_BC_3D
!>@brief The subroutine 'nested_grid_BC_send' sends coarse-grid data to create boundary conditions.
 subroutine nested_grid_BC_send_scalar(var_coarse, nest_domain, istag, jstag, nest_level)

   real, dimension(:,:,:), intent(IN) :: var_coarse
   type(nest_domain_type), intent(INOUT) :: nest_domain
   integer, intent(IN) :: istag, jstag
   integer, intent(IN) :: nest_level

   integer                      :: position

   real :: wbuffer(1,1,1)
   real :: ebuffer(1,1,1)
   real :: sbuffer(1,1,1)
   real :: nbuffer(1,1,1)


   if (istag == 1 .and. jstag == 1) then
      position = CORNER
   else if (istag == 0 .and. jstag == 1) then
      position = NORTH
   else if (istag == 1 .and. jstag == 0) then
      position = EAST
   else
      position = CENTER
   end if

       call timing_on ('COMM_TOTAL')
   call mpp_update_nest_fine(var_coarse, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level=nest_level, position=position)
       call timing_off('COMM_TOTAL')

 end subroutine nested_grid_BC_send_scalar

 subroutine nested_grid_BC_recv_scalar(nest_domain, istag, jstag, npz, &
      bd, nest_BC_buffers, nest_level)

   type(fv_grid_bounds_type), intent(IN) :: bd
   type(nest_domain_type), intent(INOUT) :: nest_domain
   integer, intent(IN) :: istag, jstag, npz
   integer, intent(IN) :: nest_level

   type(fv_nest_BC_type_3d), intent(INOUT), target :: nest_BC_buffers

   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,npz) :: var_coarse_dummy

   integer                      :: position

!!$   integer                      :: isw_f, iew_f, jsw_f, jew_f, isw_c, iew_c, jsw_c, jew_c
!!$   integer                      :: ise_f, iee_f, jse_f, jee_f, ise_c, iee_c, jse_c, jee_c
!!$   integer                      :: iss_f, ies_f, jss_f, jes_f, iss_c, ies_c, jss_c, jes_c
!!$   integer                      :: isn_f, ien_f, jsn_f, jen_f, isn_c, ien_c, jsn_c, jen_c

   integer :: i,j, k

   if (istag == 1 .and. jstag == 1) then
      position = CORNER
   else if (istag == 0 .and. jstag == 1) then
      position = NORTH
   else if (istag == 1 .and. jstag == 0) then
      position = EAST
   else
      position = CENTER
   end if

   if (.not. allocated(nest_BC_buffers%west_t1) ) then
      call init_nest_bc_type(nest_domain, nest_BC_buffers, npz, nest_level, position)
   endif

       call timing_on ('COMM_TOTAL')
   call mpp_update_nest_fine(var_coarse_dummy, nest_domain, nest_BC_buffers%west_t1, nest_BC_buffers%south_t1, &
            nest_BC_buffers%east_t1, nest_BC_buffers%north_t1, nest_level=nest_level, position=position)
       call timing_off('COMM_TOTAL')

 end subroutine nested_grid_BC_recv_scalar

 subroutine nested_grid_BC_send_vector(u_coarse, v_coarse, nest_domain, nest_level, flags, gridtype)
   real, dimension(:,:,:), intent(IN)    :: u_coarse, v_coarse
   type(nest_domain_type), intent(INOUT) :: nest_domain
   integer,                intent(IN)    :: nest_level
   integer, optional,      intent(IN)    :: flags, gridtype

   real :: wbufferx(1,1,1), wbuffery(1,1,1)
   real :: ebufferx(1,1,1), ebuffery(1,1,1)
   real :: sbufferx(1,1,1), sbuffery(1,1,1)
   real :: nbufferx(1,1,1), nbuffery(1,1,1)

   integer :: nl = 1

       call timing_on ('COMM_TOTAL')
   call mpp_update_nest_fine(u_coarse, v_coarse, nest_domain, wbufferx,wbuffery, sbufferx, sbuffery,  &
                             ebufferx, ebuffery, nbufferx, nbuffery, nest_level=nest_level, flags=flags, gridtype=gridtype)
       call timing_off('COMM_TOTAL')

 end subroutine nested_grid_BC_send_vector

 subroutine init_nest_bc_type(nest_domain, nest_BC_buffers, npz, nest_level, position)
   type(nest_domain_type),   intent(INOUT) :: nest_domain
   type(fv_nest_BC_type_3d), intent(INOUT) :: nest_BC_buffers
   integer,                  intent(IN)    :: npz, position, nest_level

   integer                      :: isw_f, iew_f, jsw_f, jew_f, isw_c, iew_c, jsw_c, jew_c
   integer                      :: ise_f, iee_f, jse_f, jee_f, ise_c, iee_c, jse_c, jee_c
   integer                      :: iss_f, ies_f, jss_f, jes_f, iss_c, ies_c, jss_c, jes_c
   integer                      :: isn_f, ien_f, jsn_f, jen_f, isn_c, ien_c, jsn_c, jen_c
   integer                      :: i, j, k

      call mpp_get_C2F_index(nest_domain, isw_f, iew_f, jsw_f, jew_f, isw_c, iew_c, jsw_c, jew_c, &
           WEST, nest_level=nest_level, position=position)
      call mpp_get_C2F_index(nest_domain, ise_f, iee_f, jse_f, jee_f, ise_c, iee_c, jse_c, jee_c, &
           EAST, nest_level=nest_level, position=position)
      call mpp_get_C2F_index(nest_domain, iss_f, ies_f, jss_f, jes_f, iss_c, ies_c, jss_c, jes_c, &
           SOUTH, nest_level=nest_level, position=position)
      call mpp_get_C2F_index(nest_domain, isn_f, ien_f, jsn_f, jen_f, isn_c, ien_c, jsn_c, jen_c, &
           NORTH, nest_level=nest_level, position=position)

      if( iew_c .GE. isw_c .AND. jew_c .GE. jsw_c ) then
         If (.not. allocated(nest_BC_buffers%west_t1)) allocate(nest_BC_buffers%west_t1(isw_c:iew_c, jsw_c:jew_c,npz))
         !compatible with first touch principle
!OMP parallel do default(none) shared(npz,jsw_c,jew_c,isw_c,iew_c,nest_BC_buffers)
         do k=1,npz
         do j=jsw_c,jew_c
         do i=isw_c,iew_c
            nest_BC_buffers%west_t1(i,j,k) = 1.e24
         enddo
         enddo
         enddo
      else
         allocate(nest_BC_buffers%west_t1(1,1,1))
         nest_BC_buffers%west_t1(1,1,1) = 1.e24
      endif

      if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
         If (.not. allocated(nest_BC_buffers%east_t1)) allocate(nest_BC_buffers%east_t1(ise_c:iee_c, jse_c:jee_c,npz))
!OMP parallel do default(none) shared(npz,jse_c,jee_c,ise_c,iee_c,nest_BC_buffers)
         do k=1,npz
         do j=jse_c,jee_c
         do i=ise_c,iee_c
            nest_BC_buffers%east_t1(i,j,k) = 1.e24
         enddo
         enddo
         enddo
      else
         allocate(nest_BC_buffers%east_t1(1,1,1))
         nest_BC_buffers%east_t1(1,1,1) = 1.e24
      endif

      if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
         If (.not. allocated(nest_BC_buffers%south_t1)) allocate(nest_BC_buffers%south_t1(iss_c:ies_c, jss_c:jes_c,npz))
!OMP parallel do default(none) shared(npz,jss_c,jes_c,iss_c,ies_c,nest_BC_buffers)
         do k=1,npz
         do j=jss_c,jes_c
         do i=iss_c,ies_c
            nest_BC_buffers%south_t1(i,j,k) = 1.e24
         enddo
         enddo
         enddo
      else
         allocate(nest_BC_buffers%south_t1(1,1,1))
         nest_BC_buffers%south_t1(1,1,1) = 1.e24
      endif

      if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
         If (.not. allocated(nest_BC_buffers%north_t1)) allocate(nest_BC_buffers%north_t1(isn_c:ien_c, jsn_c:jen_c,npz))
!OMP parallel do default(none) shared(npz,jsn_c,jen_c,isn_c,ien_c,nest_BC_buffers)
         do k=1,npz
         do j=jsn_c,jen_c
         do i=isn_c,ien_c
            nest_BC_buffers%north_t1(i,j,k) = 1.e24
         enddo
         enddo
         enddo
      else
         allocate(nest_BC_buffers%north_t1(1,1,1))
         nest_BC_buffers%north_t1(1,1,1) = 1.e24
      endif


 end subroutine init_nest_bc_type

 subroutine nested_grid_BC_recv_vector(nest_domain, npz, bd, nest_BC_u_buffers, nest_BC_v_buffers, nest_level, flags, gridtype)

   type(fv_grid_bounds_type), intent(IN) :: bd
   type(nest_domain_type), intent(INOUT) :: nest_domain
   integer, intent(IN) :: npz
   type(fv_nest_BC_type_3d), intent(INOUT), target :: nest_BC_u_buffers, nest_BC_v_buffers
   integer,                  intent(IN)  :: nest_level
   integer, optional, intent(IN) :: flags, gridtype

   real, dimension(1,1,npz) :: u_coarse_dummy, v_coarse_dummy

   integer :: i,j, k
   integer :: position_x, position_y

   call get_vector_position(position_x, position_y, gridtype)

   if (.not. allocated(nest_BC_u_buffers%west_t1) ) then
      call init_nest_bc_type(nest_domain, nest_BC_u_buffers, npz, nest_level, position_x)
   endif
   if (.not. allocated(nest_BC_v_buffers%west_t1) ) then
      call init_nest_bc_type(nest_domain, nest_BC_v_buffers, npz, nest_level, position_y)
   endif

       call timing_on ('COMM_TOTAL')
   call mpp_update_nest_fine(u_coarse_dummy, v_coarse_dummy, nest_domain, &
        nest_BC_u_buffers%west_t1, nest_BC_v_buffers%west_t1, nest_BC_u_buffers%south_t1, nest_BC_v_buffers%south_t1,  &
        nest_BC_u_buffers%east_t1, nest_BC_v_buffers%east_t1, nest_BC_u_buffers%north_t1, nest_BC_v_buffers%north_t1,  &
        nest_level, flags, gridtype)
       call timing_off('COMM_TOTAL')

 end subroutine nested_grid_BC_recv_vector

!>@brief The subroutine 'nested_grid_BC_save_proc' saves data received by 'nested_grid_BC_recv'
!! into the datatype 'fv_nest_BC_type'.
 subroutine nested_grid_BC_save_proc(nest_domain, ind, wt, istag, jstag, &
      npx, npy, npz, bd, nest_BC, nest_BC_buffers, pd_in)

   type(fv_grid_bounds_type), intent(IN) :: bd
   type(nest_domain_type), intent(INOUT) :: nest_domain
   integer, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,2), intent(IN) :: ind
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,4), intent(IN) :: wt
   integer, intent(IN) :: istag, jstag, npx, npy, npz
   logical, intent(IN), OPTIONAL :: pd_in

   !!NOTE: if declaring an ALLOCATABLE array with intent(OUT), the resulting dummy array
   !!      will NOT be allocated! This goes for allocatable members of derived types as well.
   type(fv_nest_BC_type_3d), intent(INOUT), target :: nest_BC, nest_BC_buffers

   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,npz) :: var_coarse_dummy

   real, dimension(:,:,:), pointer :: var_east, var_west, var_south, var_north
   real, dimension(:,:,:), pointer :: buf_east, buf_west, buf_south, buf_north

   integer                      :: position


   integer :: i,j, k, ic, jc, istart, iend
   logical :: process, pd = .false.

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


   if (present(pd_in)) then
      pd = pd_in
   else
      pd = .false.
   endif


      var_east  => nest_BC%east_t1
      var_west  => nest_BC%west_t1
      var_north => nest_BC%north_t1
      var_south => nest_BC%south_t1

   buf_east  => nest_BC_buffers%east_t1
   buf_west  => nest_BC_buffers%west_t1
   buf_north => nest_BC_buffers%north_t1
   buf_south => nest_BC_buffers%south_t1
   ! ?buffer has uninterpolated coarse-grid data; need to perform interpolation ourselves
   !To do this more securely, instead of using is/etc we could use the fine-grid indices defined above
   if (is == 1  ) then

!$OMP parallel do default(none) shared(npz,isd,ied,jsd,jed,jstag,ind,var_west,wt,buf_west) private(ic,jc)
      do k=1,npz
      do j=jsd,jed+jstag
         do i=isd,0

            ic = ind(i,j,1)
            jc = ind(i,j,2)


            var_west(i,j,k) = &
                 wt(i,j,1)*buf_west(ic,  jc,k) +  &
                 wt(i,j,2)*buf_west(ic,  jc+1,k) +  &
                 wt(i,j,3)*buf_west(ic+1,jc+1,k) +  &
                 wt(i,j,4)*buf_west(ic+1,jc,k)

         end do
      end do
      end do

      if (pd) then
!$OMP parallel do default(none) shared(npz,jsd,jed,jstag,isd,var_west,nest_BC)
         do k=1,npz
         do j=jsd,jed+jstag
         do i=isd,0

            var_west(i,j,k) = max(var_west(i,j,k), 0.5*nest_BC%west_t0(i,j,k))
         end do
         end do
         end do
      endif

   end if

   if (js == 1  ) then

      if (is == 1) then
         istart = is
      else
         istart = isd
      end if

      if (ie == npx-1) then
         iend = ie
      else
         iend = ied
      end if

!$OMP parallel do default(none) shared(npz,istart,iend,jsd,jed,istag,ind,var_south,wt,buf_south) private(ic,jc)
      do k=1,npz
      do j=jsd,0
         do i=istart,iend+istag

            ic = ind(i,j,1)
            jc = ind(i,j,2)


            var_south(i,j,k) = &
                 wt(i,j,1)*buf_south(ic,  jc,k) +  &
                 wt(i,j,2)*buf_south(ic,  jc+1,k) +  &
                 wt(i,j,3)*buf_south(ic+1,jc+1,k) +  &
                 wt(i,j,4)*buf_south(ic+1,jc,k)

         end do
      end do
      end do

      if (pd) then
!$OMP parallel do default(none) shared(npz,jsd,jed,istart,iend,istag,var_south,nest_BC)
         do k=1,npz
         do j=jsd,0
         do i=istart,iend+istag

            var_south(i,j,k) = max(var_south(i,j,k), 0.5*nest_BC%south_t0(i,j,k))

         end do
         end do
         end do
      endif

   end if


   if (ie == npx-1 ) then

!$OMP parallel do default(none) shared(npx,npz,isd,ied,jsd,jed,istag,jstag,ind,var_east,wt,buf_east) private(ic,jc)
      do k=1,npz
      do j=jsd,jed+jstag
         do i=npx+istag,ied+istag

            ic = ind(i,j,1)
            jc = ind(i,j,2)


            var_east(i,j,k) = &
                 wt(i,j,1)*buf_east(ic,  jc,k) +  &
                 wt(i,j,2)*buf_east(ic,  jc+1,k) +  &
                 wt(i,j,3)*buf_east(ic+1,jc+1,k) +  &
                 wt(i,j,4)*buf_east(ic+1,jc,k)

         end do
      end do
      end do

      if (pd) then
!$OMP parallel do default(none) shared(npx,npz,jsd,jed,istag,jstag,ied,var_east,nest_BC)
         do k=1,npz
         do j=jsd,jed+jstag
         do i=npx+istag,ied+istag

            var_east(i,j,k) = max(var_east(i,j,k), 0.5*nest_BC%east_t0(i,j,k))

         end do
         end do
         end do
      endif

   end if

   if (je == npy-1 ) then

      if (is == 1) then
         istart = is
      else
         istart = isd
      end if

      if (ie == npx-1) then
         iend = ie
      else
         iend = ied
      end if

!$OMP parallel do default(none) shared(npy,npz,istart,iend,jsd,jed,istag,jstag,ind,var_north,wt,buf_north) private(ic,jc)
      do k=1,npz
      do j=npy+jstag,jed+jstag
         do i=istart,iend+istag

            ic = ind(i,j,1)
            jc = ind(i,j,2)


            var_north(i,j,k) = &
                 wt(i,j,1)*buf_north(ic,  jc,k) +  &
                 wt(i,j,2)*buf_north(ic,  jc+1,k) +  &
                 wt(i,j,3)*buf_north(ic+1,jc+1,k) +  &
                 wt(i,j,4)*buf_north(ic+1,jc,k)

         end do
      end do
      end do

      if (pd) then
!$OMP parallel do default(none) shared(npy,npz,jsd,jed,istart,iend,istag,jstag,ied,var_north,nest_BC)
         do k=1,npz
         do j=npy+jstag,jed+jstag
         do i=istart,iend+istag

            var_north(i,j,k) = max(var_north(i,j,k), 0.5*nest_BC%north_t0(i,j,k))

         end do
         end do
         end do
      endif

   end if

 end subroutine nested_grid_BC_save_proc


  ! A NOTE ON BCTYPE: currently only an interpolation BC is implemented,
  ! bctype >= 2 currently correspond
  ! to a flux BC on the tracers ONLY, which is implemented in fv_tracer.

!>@brief The subroutine 'nested_grid_BC_apply_intT' performs linear interpolation or
!! extrapolation in time for saved BC data, then applies the interlpolated
!! data to nested-grid boundary cells.
 subroutine nested_grid_BC_apply_intT(var_nest, istag, jstag, &
      npx, npy, npz, bd, step, split, &
      BC, bctype)

   type(fv_grid_bounds_type), intent(IN) :: bd
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag, npz), intent(INOUT) :: var_nest
   integer, intent(IN) :: istag, jstag, npx, npy, npz
   real, intent(IN) :: split, step
   integer, intent(IN) :: bctype

   type(fv_nest_BC_type_3D), intent(IN), target :: BC
   real, pointer, dimension(:,:,:) :: var_t0, var_t1

   integer :: i,j, istart, iend, k
   real :: denom

   logical, save :: printdiag = .true.

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

   denom = 1./split
   if (is == 1  ) then
      var_t0 => BC%west_t0
      var_t1 => BC%west_t1
!OMP parallel do default(none) shared(npz,jsd,jed,jstag,isd,var_nest,var_t0,var_t1,split,step,denom)
      do k=1,npz
      do j=jsd,jed+jstag
      do i=isd,0
            var_nest(i,j,k) = (var_t0(i,j,k)*(split-step) + step*var_t1(i,j,k))*denom
      end do
      end do
      end do
   end if

   if (js == 1  ) then

      if (is == 1) then
         istart = is
      else
         istart = isd
      end if

      if (ie == npx-1) then
         iend = ie
      else
         iend = ied
      end if

      var_t0 => BC%south_t0
      var_t1 => BC%south_t1
!OMP parallel do default(none) shared(npz,jsd,istart,iend,istag,var_nest,var_t0,var_t1,split,step,denom)
      do k=1,npz
      do j=jsd,0
      do i=istart,iend+istag
            var_nest(i,j,k) = (var_t0(i,j,k)*(split-step) + step*var_t1(i,j,k))*denom
      end do
      end do
      end do
   end if


   if (ie == npx-1 ) then
      var_t0 => BC%east_t0
      var_t1 => BC%east_t1
!OMP parallel do default(none) shared(npz,jsd,jed,jstag,npx,isd,istag,var_nest,var_t0,var_t1,split,step,denom)
      do k=1,npz
      do j=jsd,jed+jstag
      do i=npx+istag,ied+istag
         var_nest(i,j,k) = (var_t0(i,j,k)*(split-step) + step*var_t1(i,j,k))*denom
      end do
      end do
      end do
   end if

   if (je == npy-1 ) then

      if (is == 1) then
         istart = is
      else
         istart = isd
      end if

      if (ie == npx-1) then
         iend = ie
      else
         iend = ied
      end if

      var_t0 => BC%north_t0
      var_t1 => BC%north_t1
!OMP parallel do default(none) shared(npz,npy,jed,jstag,istart,iend,istag,var_nest,var_t0,var_t1,split,step,denom)
      do k=1,npz
      do j=npy+jstag,jed+jstag
      do i=istart,iend+istag
         var_nest(i,j,k) = (var_t0(i,j,k)*(split-step) + step*var_t1(i,j,k))*denom
      end do
      end do
      end do

   end if


 end subroutine nested_grid_BC_apply_intT

 subroutine update_coarse_grid_mpp_2d(var_coarse, var_nest, nest_domain, dx, dy, area, &
      bd, isd_p, ied_p, jsd_p, jed_p, is_n, ie_n, js_n, je_n, isu, ieu, jsu, jeu, npx, npy, &
      istag, jstag, r, nestupdate, upoff, nsponge, parent_proc, child_proc, parent_grid, nest_level)

   type(fv_grid_bounds_type), intent(IN) :: bd
   integer, intent(IN) :: isd_p, ied_p, jsd_p, jed_p, is_n, ie_n, js_n, je_n
   integer, intent(IN) :: isu, ieu, jsu, jeu
   integer, intent(IN) :: istag, jstag, r, nestupdate, upoff, nsponge
   integer, intent(IN) :: npx, npy
   real, intent(IN), target    :: var_nest(is_n:ie_n+istag,js_n:je_n+jstag)
   real, intent(INOUT), target :: var_coarse(isd_p:ied_p+istag,jsd_p:jed_p+jstag)
   real, intent(IN)    :: dx(bd%isd:bd%ied,  bd%jsd:bd%jed+1)
   real, intent(IN)    :: dy(bd%isd:bd%ied+1,bd%jsd:bd%jed)
   real, intent(IN)    :: area(bd%isd:bd%ied,bd%jsd:bd%jed)
   logical, intent(IN) :: parent_proc, child_proc
   type(fv_atmos_type), pointer, intent(IN) :: parent_grid
   type(nest_domain_type), intent(INOUT) :: nest_domain
   integer, intent(IN) :: nest_level

   real :: var_nest_3d(is_n:ie_n+istag,js_n:je_n+jstag,1)
   real :: var_coarse_3d(isd_p:ied_p+istag,jsd_p:jed_p+jstag,1)
   integer( KIND = 8) :: ptr_nest=0
   integer( KIND = 8) :: ptr_coarse=0
   pointer(ptr_nest, var_nest_3d)
   pointer(ptr_coarse, var_coarse_3d)

   if (child_proc .and. size(var_nest) > 1) ptr_nest = LOC(var_nest)
   if (parent_proc .and. size(var_coarse) > 1) ptr_coarse = LOC(var_coarse)

   call update_coarse_grid_mpp(var_coarse_3d, var_nest_3d, &
        nest_domain, dx, dy, area, &
        bd, isd_p, ied_p, jsd_p, jed_p, is_n, ie_n, js_n, je_n, &
        isu, ieu, jsu, jeu, npx, npy, 1, &
        istag, jstag, r, nestupdate, upoff, nsponge, &
        parent_proc, child_proc, parent_grid, nest_level )

 end subroutine update_coarse_grid_mpp_2d


  subroutine update_coarse_grid_mpp(var_coarse, var_nest, nest_domain, dx, dy, area, &
      bd, isd_p, ied_p, jsd_p, jed_p, is_n, ie_n, js_n, je_n, &
      isu, ieu, jsu, jeu, npx, npy, npz, &
      istag, jstag, r, nestupdate, upoff, nsponge, &
      parent_proc, child_proc, parent_grid, nest_level)

   !This routine assumes the coarse and nested grids are properly
   ! aligned, and that in particular for odd refinement ratios all
   ! coarse-grid cells (faces) coincide with nested-grid cells (faces)

   type(fv_grid_bounds_type), intent(IN) :: bd
   integer, intent(IN) :: isd_p, ied_p, jsd_p, jed_p, is_n, ie_n, js_n, je_n
   integer, intent(IN) :: isu, ieu, jsu, jeu
   integer, intent(IN) :: istag, jstag, npx, npy, npz, r, nestupdate, upoff, nsponge
   real, intent(IN)    :: var_nest(is_n:ie_n+istag,js_n:je_n+jstag,npz)
   real, intent(INOUT) :: var_coarse(isd_p:ied_p+istag,jsd_p:jed_p+jstag,npz)
   real, intent(IN)    :: area(bd%isd:bd%ied,bd%jsd:bd%jed)
   real, intent(IN)    :: dx(bd%isd:bd%ied,bd%jsd:bd%jed+1)
   real, intent(IN)    :: dy(bd%isd:bd%ied+1,bd%jsd:bd%jed)
   logical, intent(IN) :: parent_proc, child_proc
   type(fv_atmos_type), pointer, intent(IN) :: parent_grid
   type(nest_domain_type), intent(INOUT) :: nest_domain
   integer, intent(IN) :: nest_level

   integer :: in, jn, ini, jnj, s, qr
   integer :: is_c, ie_c, js_c, je_c, is_f, ie_f, js_f, je_f
   integer :: istart, istop, jstart, jstop, ishift, jshift, j, i, k
   real :: val
   real, allocatable, dimension(:,:,:) :: coarse_dat_send
   real, allocatable ::  coarse_dat_recv(:,:,:)
   integer :: position

   if (istag == 1 .and. jstag == 1) then
      position = CORNER
   else if (istag == 0 .and. jstag == 1) then
      position = NORTH
   else if (istag == 1 .and. jstag == 0) then
      position = EAST
   else
      position = CENTER
   end if

   !Note that *_c does not have values on the parent_proc.
   !Must use isu, etc. to get bounds of update region on parent.
   call mpp_get_F2C_index(nest_domain, is_c, ie_c, js_c, je_c, is_f, ie_f, js_f, je_f, nest_level=nest_level, position=position)
   if (child_proc) then
      allocate(coarse_dat_send(is_c:ie_c, js_c:je_c,npz))
      coarse_dat_send = -1200.
   endif
   allocate(coarse_dat_recv(isd_p:ied_p+istag, jsd_p:jed_p+jstag, npz))

   if (child_proc) then
      call fill_coarse_data_send(coarse_dat_send, var_nest, dx, dy, area, &
            bd, is_c, ie_c, js_c, je_c, is_f, js_f, is_n, ie_n, js_n, je_n, &
            npx, npy, npz, istag, jstag, r, nestupdate)
   endif

      call timing_on('COMM_TOTAL')
   call mpp_update_nest_coarse(field_in=coarse_dat_send, nest_domain=nest_domain, field_out=coarse_dat_recv, &
        nest_level=nest_level, position=position)

   if (allocated(coarse_dat_send)) then
      deallocate(coarse_dat_send)
   end if

      call timing_off('COMM_TOTAL')

   s = r/2 !rounds down (since r > 0)
   qr = r*upoff + nsponge - s

   if (parent_proc .and. .not. (ieu < isu .or. jeu < jsu)) then
      call fill_var_coarse(var_coarse, coarse_dat_recv, isd_p, ied_p, jsd_p, jed_p, &
           isu, ieu, jsu, jeu, npx, npy, npz, istag, jstag, nestupdate, parent_grid)
   endif

   if (allocated(coarse_dat_recv)) deallocate(coarse_dat_recv)

 end subroutine update_coarse_grid_mpp

 subroutine fill_coarse_data_send(coarse_dat_send, var_nest, dx, dy, area, &
      bd, is_c, ie_c, js_c, je_c, is_f, js_f, is_n, ie_n, js_n, je_n, &
      npx, npy, npz, istag, jstag, r, nestupdate)
   type(fv_grid_bounds_type), intent(IN) :: bd
   integer, intent(IN) :: is_c, ie_c, js_c, je_c, is_n, ie_n, js_n, je_n
   integer, intent(IN) :: is_f, js_f
   integer, intent(IN) :: istag, jstag
   integer, intent(IN) :: npx, npy, npz, r, nestupdate
   real, intent(INOUT) :: coarse_dat_send(is_c:ie_c,js_c:je_c,npz)
   real, intent(IN)    :: var_nest(is_n:ie_n+istag,js_n:je_n+jstag,npz)
   real, intent(IN)    :: area(bd%isd:bd%ied,bd%jsd:bd%jed)
   real, intent(IN)    :: dx(bd%isd:bd%ied,bd%jsd:bd%jed+1)
   real, intent(IN)    :: dy(bd%isd:bd%ied+1,bd%jsd:bd%jed)
   integer :: in, jn, ini, jnj, k, j, i
   real :: val


   if (istag == 0 .and. jstag == 0) then
      select case (nestupdate)
      case (1,2,6,7,8)

!$OMP parallel do default(none) shared(npz,js_c,je_c,is_c,ie_c,js_f,is_f,coarse_dat_send,var_nest,area,r) private(in,jn,val)
         do k=1,npz
            jn = js_f
         do j=js_c,je_c
            in = is_f
         do i=is_c,ie_c

            val = 0.
            do jnj=jn,jn+r-1
               do ini=in,in+r-1
                  val = val + var_nest(ini,jnj,k)*area(ini,jnj)
               end do
            end do
            coarse_dat_send(i,j,k) = val !divide area on coarse grid

            in = in + r
         end do
            jn = jn + r
         end do
         end do

      end select
   else if (istag == 0 .and. jstag > 0) then

      select case (nestupdate)
      case (1,6,7,8)

!$OMP parallel do default(none) shared(npz,js_c,je_c,is_c,ie_c,js_f,is_f,coarse_dat_send,var_nest,dx,r) private(in,jn,val)
         do k=1,npz
            jn = js_f
         do j=js_c,je_c!+1
            in = is_f
         do i=is_c,ie_c

            val = 0.
            do ini=in,in+r-1
               val = val + var_nest(ini,jn,k)*dx(ini,jn)
            end do
            coarse_dat_send(i,j,k) = val

            in = in + r
         end do
            jn = jn + r
         end do
         end do

      case default

         call mpp_error(FATAL, 'nestupdate type not implemented')

      end select

   else if (istag > 0 .and. jstag == 0) then
      select case (nestupdate)

      case (1,6,7,8)   !averaging update; in-line average for face-averaged values instead of areal average

!$OMP parallel do default(none) shared(npz,js_c,je_c,is_c,ie_c,js_f,is_f,coarse_dat_send,var_nest,dy,r) private(in,jn,val)
         do k=1,npz
            jn = js_f
         do j=js_c,je_c
            in = is_f
         do i=is_c,ie_c!+1

            val = 0.
            do jnj=jn,jn+r-1
                  val = val + var_nest(in,jnj,k)*dy(in,jnj)
            end do
            coarse_dat_send(i,j,k) = val

            in = in + r
         end do
            jn = jn + r
         end do
         end do

      case default

         call mpp_error(FATAL, 'nestupdate type not implemented')

      end select

   else

      call mpp_error(FATAL, "Cannot have both nonzero istag and jstag.")

   endif


 end subroutine fill_coarse_data_send

 subroutine fill_var_coarse(var_coarse, coarse_dat_recv, isd_p, ied_p, jsd_p, jed_p, &
      isu, ieu, jsu, jeu, npx, npy, npz, istag, jstag, nestupdate, parent_grid)

   !This routine assumes the coarse and nested grids are properly
   ! aligned, and that in particular for odd refinement ratios all
   ! coarse-grid cells (faces) coincide with nested-grid cells (faces)

   integer, intent(IN) :: isd_p, ied_p, jsd_p, jed_p
   integer, intent(IN) :: isu, ieu, jsu, jeu
   integer, intent(IN) :: istag, jstag
   integer, intent(IN) :: npx, npy, npz, nestupdate
   real, intent(INOUT) :: var_coarse(isd_p:ied_p+istag,jsd_p:jed_p+jstag,npz)
   real, intent(INOUT) :: coarse_dat_recv(isd_p:ied_p+istag,jsd_p:jed_p+jstag,npz)
   type(fv_atmos_type), intent(IN) :: parent_grid

   integer :: i, j, k

   if (istag == 0 .and. jstag == 0) then

      select case (nestupdate)
      case (1,2,6,7,8) ! 1 = Conserving update on all variables; 2 = conserving update for cell-centered values; 6 = conserving remap-update

!$OMP parallel do default(none) shared(npz,jsu,jeu,isu,ieu,coarse_dat_recv,parent_grid,var_coarse)
         do k=1,npz
         do j=jsu,jeu
         do i=isu,ieu
            var_coarse(i,j,k) = coarse_dat_recv(i,j,k)*parent_grid%gridstruct%rarea(i,j)
         end do
         end do
         end do


      case default

         call mpp_error(FATAL, 'nestupdate type not implemented')


      end select

   else if (istag == 0 .and. jstag > 0) then

      select case (nestupdate)
      case (1,6,7,8)

!$OMP parallel do default(none) shared(npz,jsu,jeu,isu,ieu,coarse_dat_recv,parent_grid,var_coarse)
         do k=1,npz
         do j=jsu,jeu+1
         do i=isu,ieu
            var_coarse(i,j,k) = coarse_dat_recv(i,j,k)*parent_grid%gridstruct%rdx(i,j)
         end do
         end do
         end do

      case default

         call mpp_error(FATAL, 'nestupdate type not implemented')

      end select

   else if (istag > 0 .and. jstag == 0) then

      select case (nestupdate)
      case (1,6,7,8)   !averaging update; in-line average for face-averaged values instead of areal average

!$OMP parallel do default(none) shared(npz,jsu,jeu,isu,ieu,coarse_dat_recv,parent_grid,var_coarse)
         do k=1,npz
         do j=jsu,jeu
         do i=isu,ieu+1
            var_coarse(i,j,k) = coarse_dat_recv(i,j,k)*parent_grid%gridstruct%rdy(i,j)
         end do
         end do
         end do

      case default

         call mpp_error(FATAL, 'nestupdate type not implemented')

      end select

   end if


 end subroutine fill_var_coarse

  subroutine update_coarse_grid_mpp_vector(u_coarse, v_coarse, u_nest, v_nest, nest_domain, dx, dy, area, &
      bd, isd_p, ied_p, jsd_p, jed_p, is_n, ie_n, js_n, je_n, &
      isu, ieu, jsu, jeu, npx, npy, npz, istag_u, jstag_u, istag_v, jstag_v, &
      r, nestupdate, upoff, nsponge, &
      parent_proc, child_proc, parent_grid, nest_level, flags, gridtype)

   !This routine assumes the coarse and nested grids are properly
   ! aligned, and that in particular for odd refinement ratios all
   ! coarse-grid cells (faces) coincide with nested-grid cells (faces)

   type(fv_grid_bounds_type), intent(IN) :: bd
   integer, intent(IN) :: isd_p, ied_p, jsd_p, jed_p, is_n, ie_n, js_n, je_n
   integer, intent(IN) :: isu, ieu, jsu, jeu
   integer, intent(IN) :: istag_u, jstag_u, istag_v, jstag_v
   integer, intent(IN) :: npx, npy, npz, r, nestupdate, upoff, nsponge
   real, intent(IN)    :: u_nest(is_n:ie_n+istag_u,js_n:je_n+jstag_u,npz)
   real, intent(INOUT) :: u_coarse(isd_p:ied_p+istag_u,jsd_p:jed_p+jstag_u,npz)
   real, intent(IN)    :: v_nest(is_n:ie_n+istag_v,js_n:je_n+jstag_v,npz)
   real, intent(INOUT) :: v_coarse(isd_p:ied_p+istag_v,jsd_p:jed_p+jstag_v,npz)
   real, intent(IN)    :: area(bd%isd:bd%ied,bd%jsd:bd%jed)
   real, intent(IN)    :: dx(bd%isd:bd%ied,bd%jsd:bd%jed+1)
   real, intent(IN)    :: dy(bd%isd:bd%ied+1,bd%jsd:bd%jed)
   logical, intent(IN) :: parent_proc, child_proc
   type(fv_atmos_type), intent(INOUT) :: parent_grid
   integer, intent(IN) :: nest_level
   type(nest_domain_type), intent(INOUT) :: nest_domain
   integer, optional,      intent(IN)    :: flags, gridtype

   integer :: s, qr
   integer :: is_cx, ie_cx, js_cx, je_cx, is_fx, ie_fx, js_fx, je_fx
   integer :: is_cy, ie_cy, js_cy, je_cy, is_fy, ie_fy, js_fy, je_fy
   integer :: istart, istop, jstart, jstop, ishift, jshift, j, i, k
   real :: val
   real, allocatable, dimension(:,:,:) :: coarse_dat_send_u, coarse_dat_send_v
   real, allocatable ::  coarse_dat_recv_u(:,:,:), coarse_dat_recv_v(:,:,:)
   integer :: position_x, position_y

   call get_vector_position(position_x, position_y, gridtype)

   call mpp_get_F2C_index(nest_domain, is_cx, ie_cx, js_cx, je_cx, is_fx, ie_fx, js_fx, je_fx, &
        nest_level=nest_level, position=position_x)
   call mpp_get_F2C_index(nest_domain, is_cy, ie_cy, js_cy, je_cy, is_fy, ie_fy, js_fy, je_fy, &
        nest_level=nest_level, position=position_y)
   if (child_proc) then
      allocate(coarse_dat_send_u(is_cx:ie_cx, js_cx:je_cx,npz))
      allocate(coarse_dat_send_v(is_cy:ie_cy, js_cy:je_cy,npz))
      coarse_dat_send_u = -1200.
      coarse_dat_send_v = -1200.
   endif
   allocate(coarse_dat_recv_u(isd_p:ied_p+istag_u, jsd_p:jed_p+jstag_u, npz))
   allocate(coarse_dat_recv_v(isd_p:ied_p+istag_v, jsd_p:jed_p+jstag_v, npz))

   if (child_proc) then
      call fill_coarse_data_send(coarse_dat_send_u, u_nest, dx, dy, area, &
           bd, is_cx, ie_cx, js_cx, je_cx, is_fx, js_fx, is_n, ie_n, js_n, je_n, &
           npx, npy, npz, istag_u, jstag_u, r, nestupdate)
      call fill_coarse_data_send(coarse_dat_send_v, v_nest, dx, dy, area, &
           bd, is_cy, ie_cy, js_cy, je_cy, is_fy, js_fy, is_n, ie_n, js_n, je_n, &
           npx, npy, npz, istag_v, jstag_v, r, nestupdate)
   endif

      call timing_on('COMM_TOTAL')
   call mpp_update_nest_coarse(coarse_dat_send_u, coarse_dat_send_v, nest_domain, coarse_dat_recv_u, &
                               coarse_dat_recv_v, nest_level, flags, gridtype)

   if (allocated(coarse_dat_send_u)) deallocate(coarse_dat_send_u)
   if (allocated(coarse_dat_send_v)) deallocate(coarse_dat_send_v)

      call timing_off('COMM_TOTAL')

   s = r/2 !rounds down (since r > 0)
   qr = r*upoff + nsponge - s

   if (parent_proc .and. .not. (ieu < isu .or. jeu < jsu)) then
      call fill_var_coarse(u_coarse, coarse_dat_recv_u, isd_p, ied_p, jsd_p, jed_p, &
           isu, ieu, jsu, jeu, npx, npy, npz, istag_u, jstag_u, nestupdate, parent_grid)
   endif
   if (parent_proc .and. .not. (ieu < isu .or. jeu < jsu)) then
      call fill_var_coarse(v_coarse, coarse_dat_recv_v, isd_p, ied_p, jsd_p, jed_p, &
           isu, ieu, jsu, jeu, npx, npy, npz, istag_v, jstag_v, nestupdate, parent_grid)
   endif

   if (allocated(coarse_dat_recv_u)) deallocate(coarse_dat_recv_u)
   if (allocated(coarse_dat_recv_v)) deallocate(coarse_dat_recv_v)

 end subroutine update_coarse_grid_mpp_vector

end module boundary_mod
