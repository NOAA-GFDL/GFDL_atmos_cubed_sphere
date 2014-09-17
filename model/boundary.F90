module boundary_mod

  use fv_mp_mod,         only: ng, isc,jsc,iec,jec, isd,jsd,ied,jed, is,js,ie,je, is_master
  use constants_mod,     only: grav

  use mpp_domains_mod,    only: mpp_get_compute_domain, mpp_get_data_domain, mpp_get_global_domain
  use mpp_domains_mod,    only: CENTER, CORNER, NORTH, EAST
  use mpp_domains_mod,    only: mpp_global_field, mpp_get_pelist
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
  public outflow_u, outflow_v, outflow_x, outflow_y, extrapolation_BC
  public nested_grid_bc, nested_grid_bc_save, gather_grid, update_coarse_grid
  public fill_nested_grid, nested_grid_bc_apply, nested_grid_bc_apply_intT

  real :: outflowspd_east  = 3.
  real :: outflowspd_west  = 3.
  real :: outflowspd_north = 3.
  real :: outflowspd_south = 3.

  interface gather_grid
     module procedure gather_grid_2d
     module procedure gather_grid_3d
  end interface

  interface nested_grid_BC_apply
     module procedure nested_grid_BC_apply_2d
     module procedure nested_grid_BC_apply_3d
  end interface

  interface nested_grid_BC_apply_intT
     module procedure nested_grid_BC_apply_intT_2d
     module procedure nested_grid_BC_apply_intT_3d
  end interface

  interface nested_grid_BC
     module procedure nested_grid_BC_2d
     module procedure nested_grid_BC_mpp
     module procedure nested_grid_BC_mpp_send
     module procedure nested_grid_BC_2D_mpp
     module procedure nested_grid_BC_3d
  end interface

  interface nested_grid_BC_save
     module procedure nested_grid_BC_save_2d
     module procedure nested_grid_BC_save_mpp
     module procedure nested_grid_BC_save_send
     module procedure nested_grid_BC_save_3d
  end interface

  interface fill_nested_grid
     module procedure fill_nested_grid_2d
     module procedure fill_nested_grid_3d
  end interface

  interface update_coarse_grid
     module procedure update_coarse_grid_mpp
     module procedure update_coarse_grid_mpp_2d
  end interface

  !---- version number -----
  character(len=128) :: version = '$Id: boundary.F90,v 20.0 2013/12/13 23:04:16 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal_201409 $'

contains

  !Outflow BC routines based off of the formulations of Klemp and
  ! Wilhelmson (1978), modified for d-grid or c-grid as necessary.

  subroutine outflow_u(u, rdx, dt, npx, cgrid, outflow_east, outflow_west, bd)

    implicit none
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(in) :: cgrid     ! 0 if d-grid, 1 if c-grid
    integer, intent(in) :: npx
    real, intent(inout), dimension(bd%isd:bd%ied+cgrid  ,bd%jsd:bd%jed+1-cgrid) :: u
    real, intent(in), dimension(bd%isd:bd%ied+cgrid  ,bd%jsd:bd%jed+1-cgrid) :: rdx
    real, intent(in) :: dt
    logical, intent(in) :: outflow_east, outflow_west

    integer :: j, ix
    real :: ul, ur

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

    !do we need to include any metric terms when computing these derivatives? 
    !Also these routines assume the coordinate line is perpendicular to the boundary

    if (outflow_west .AND. isc == 1) then

       do j = jsd,jed+1-cgrid
          ul = min( u(1,j) - outflowspd_west, 0.0 )
          u(1,j) = u(1,j) - dt * ul * rdx(1,j) * ( u(2,j) - u(1,j) ) 
       end do

    end if

    if (outflow_east .AND. iec+1 == npx) then

       ix = iec + cgrid

       do  j = jsd,jed+1-cgrid
          ur = max( u(ix,j) + outflowspd_east, 0.0)
          u(ix,j) = u(ix,j) - dt * ur * rdx(ix,j) * ( u(ix,j) - u(ix-1,j) )
       end do

    end if

  end subroutine outflow_u

  subroutine outflow_v(v, rdy, dt, npy, cgrid, outflow_north, outflow_south, bd)

    implicit none
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(in) :: cgrid     ! 0 if d-grid, 1 if c-grid
    integer, intent(in) :: npy
    real, intent(inout), dimension(bd%isd:bd%ied+1-cgrid  ,bd%jsd:bd%jed+cgrid) :: v
    real, intent(in), dimension(bd%isd:bd%ied+1-cgrid  ,bd%jsd:bd%jed+cgrid) :: rdy
    real, intent(in) :: dt
    logical, intent(in) :: outflow_north, outflow_south

    integer :: i, jx
    real :: vn, vs

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

    if (outflow_south .AND. jsc == 1) then

       do i = isd,ied+1-cgrid
          vs = min( v(i,1) - outflowspd_south , 0.0 )
          v(i,1) = v(i,1) - dt * vs * rdy(i,1) * ( v(i,2) - v(i,1) )
       end do

    end if

    if (outflow_north .AND. jec+1 == npy) then

       jx = jec + cgrid

       do i = isd,ied+1-cgrid
          vn = max( v(i,jx) + outflowspd_north, 0.0 )
          v(i,jx) = v(i,jx) - dt * vn * rdy(i,jx) * ( v(i,jx) - v(i,jx-1) )
       end do

    end if

  end subroutine outflow_v

  !For unstaggered scalars or boundary-parallel winds. u is assumed to be 
  ! co-located with the input variable, dx is assumed to be for the q points.
  ! The outflow wave speed is not used for non-normal flow components nor scalars.
  ! The flow speed used is that of the boundary point of the variable in question
  ! (hence the collocation); the differencing for the variable is first-order upstream
  subroutine outflow_x(u, q, rdx, dt, npx, istag, jstag, outflow_east, outflow_west, bd)

    implicit none
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(in) :: istag, jstag, npx
    real, intent(in), dimension(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag) :: u, rdx
    real, intent(inout), dimension(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag) :: q
    real, intent(in) :: dt
    logical, intent(in) :: outflow_east, outflow_west

    integer :: j, ix
    real :: ul, ur

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

    !For scalars, we need to have u interpolated to the a-grid
    !For c-grid v, we need to have u interpolated to the d-grid
    !For d-grid u, we need to have u interpolated to the c-grid

    if (outflow_west .AND. isc == 1) then 

       do j = jsd,jed+jstag
          ul = min( u(1,j), 0.0 )
          q(1,j) = q(1,j) - dt * ul * rdx(1,j) * ( q(2,j) - q(1,j) ) 
       end do

    end if

    if (outflow_east .AND. iec+1 == npx) then

       ix = iec + istag

       do  j = jsd,jed+jstag
          ur = max( u(ix,j), 0.0)
          q(ix,j) = q(ix,j) - dt * ur * rdx(ix,j) * ( q(ix,j) - q(ix-1,j) )
       end do

    end if


  end subroutine outflow_x

  subroutine outflow_y(v, q, rdy, dt, npy, istag, jstag, outflow_north, outflow_south, bd)

    implicit none
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(in) :: istag, jstag, npy
    real, intent(inout), dimension(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag) :: q
    real, intent(in), dimension(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag) :: v, rdy
    real, intent(in) :: dt
    logical, intent(in) :: outflow_north, outflow_south

    integer :: i, jx
    real :: vn, vs

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

    if (outflow_south .AND. jsc == 1) then

       do i = isd,ied+istag
          vs = min( v(i,1), 0.0 )
          q(i,1) = q(i,1) - dt * vs * rdy(i,1) * ( q(i,2) - q(i,1) )
       end do

    end if

    if (outflow_north .AND. jec+1 == npy) then

       jx = jec + jstag

       do i = isd,ied+istag
          vn = max( v(i,jx), 0.0 )
          q(i,jx) = q(i,jx) - dt * vn * rdy(i,jx) * ( q(i,jx) - q(i,jx-1) )
       end do

    end if    

  end subroutine outflow_y

  !Linear extrapolation into halo region
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
 
 subroutine nested_grid_BC_mpp(var_nest, var_coarse, nest_domain, ind, wt, istag, jstag, &
      npx, npy, npz, bd, isg, ieg, jsg, jeg, nstep_in, nsplit_in, proc_in)

   type(fv_grid_bounds_type), intent(IN) :: bd
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,npz), intent(INOUT) :: var_nest
   real, dimension(isg:ieg+istag,jsg:jeg+jstag,npz), intent(IN) :: var_coarse
   type(nest_domain_type), intent(INOUT) :: nest_domain
   integer, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,2), intent(IN) :: ind
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,4), intent(IN) :: wt
   integer, intent(IN) :: istag, jstag, npx, npy, npz, isg, ieg, jsg, jeg
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
        WEST,  position=position)
   call mpp_get_C2F_index(nest_domain, ise_f, iee_f, jse_f, jee_f, ise_c, iee_c, jse_c, jee_c, &
        EAST,  position=position)
   call mpp_get_C2F_index(nest_domain, iss_f, ies_f, jss_f, jes_f, iss_c, ies_c, jss_c, jes_c, &
        SOUTH,  position=position)
   call mpp_get_C2F_index(nest_domain, isn_f, ien_f, jsn_f, jen_f, isn_c, ien_c, jsn_c, jen_c, &
        NORTH,  position=position)

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
   call mpp_update_nest_fine(var_coarse, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer,  position=position)
       call timing_off('COMM_TOTAL')

   if (process) then

   if (is == 1) then
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
!!$   call nested_grid_BC_save_mpp(var_coarse, nest_domain, ind, wt, istag, jstag, &
!!$      npx, npy, npz, isg, ieg, jsg, jeg,var_east= ebuffer, var_west=wbuffer, &
!!$      var_north=nbuffer, var_south=sbuffer, ns=0)
!!$   
!!$   if( iew_c .GE. isw_c .AND. jew_c .GE. jsw_c ) then
!!$      var_nest(isw_c:iew_c, jsw_c:jew_c,:) = wbuffer(isw_c:iew_c, jsw_c:jew_c,:)
!!$   endif
!!$
!!$   if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
!!$      var_nest(ise_c:iee_c, jse_c:jee_c,:) = ebuffer(ise_c:iee_c, jse_c:jee_c,:)
!!$   endif
!!$
!!$   if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
!!$      var_nest(iss_c:ies_c, jss_c:jes_c,:) = sbuffer(iss_c:ies_c, jss_c:jes_c,:)
!!$   endif
!!$
!!$   if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
!!$      var_nest(isn_c:ien_c, jsn_c:jen_c,:) = nbuffer(isn_c:ien_c, jsn_c:jen_c,:)
!!$   endif

   deallocate(wbuffer, ebuffer, sbuffer, nbuffer)

 end subroutine nested_grid_BC_mpp

 subroutine nested_grid_BC_mpp_send(var_coarse, nest_domain, istag, jstag)

   real, dimension(:,:,:), intent(IN) :: var_coarse
   type(nest_domain_type), intent(INOUT) :: nest_domain
   integer, intent(IN) :: istag, jstag

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
   call mpp_update_nest_fine(var_coarse, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer,  position=position)
       call timing_off('COMM_TOTAL')


   deallocate(wbuffer, ebuffer, sbuffer, nbuffer)

 end subroutine nested_grid_BC_mpp_send

 subroutine nested_grid_BC_2D_mpp(var_nest, var_coarse, nest_domain, ind, wt, istag, jstag, &
      npx, npy, bd, isg, ieg, jsg, jeg, nstep_in, nsplit_in, proc_in)

   type(fv_grid_bounds_type), intent(IN) :: bd
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag), intent(INOUT) :: var_nest
   real, dimension(isg:ieg+istag,jsg:jeg+jstag), intent(IN) :: var_coarse
   type(nest_domain_type), intent(INOUT) :: nest_domain
   integer, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,2), intent(IN) :: ind
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,4), intent(IN) :: wt
   integer, intent(IN) :: istag, jstag, npx, npy, isg, ieg, jsg, jeg
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
        WEST,  position=position)
   call mpp_get_C2F_index(nest_domain, ise_f, iee_f, jse_f, jee_f, ise_c, iee_c, jse_c, jee_c, &
        EAST,  position=position)
   call mpp_get_C2F_index(nest_domain, iss_f, ies_f, jss_f, jes_f, iss_c, ies_c, jss_c, jes_c, &
        SOUTH,  position=position)
   call mpp_get_C2F_index(nest_domain, isn_f, ien_f, jsn_f, jen_f, isn_c, ien_c, jsn_c, jen_c, &
        NORTH,  position=position)

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
   call mpp_update_nest_fine(var_coarse, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer,  position=position)
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
!!$   call nested_grid_BC_save_mpp(var_coarse, nest_domain, ind, wt, istag, jstag, &
!!$      npx, npy, npz, isg, ieg, jsg, jeg,var_east= ebuffer, var_west=wbuffer, &
!!$      var_north=nbuffer, var_south=sbuffer, ns=0)
!!$   
!!$   if( iew_c .GE. isw_c .AND. jew_c .GE. jsw_c ) then
!!$      var_nest(isw_c:iew_c, jsw_c:jew_c,:) = wbuffer(isw_c:iew_c, jsw_c:jew_c,:)
!!$   endif
!!$
!!$   if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
!!$      var_nest(ise_c:iee_c, jse_c:jee_c,:) = ebuffer(ise_c:iee_c, jse_c:jee_c,:)
!!$   endif
!!$
!!$   if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
!!$      var_nest(iss_c:ies_c, jss_c:jes_c,:) = sbuffer(iss_c:ies_c, jss_c:jes_c,:)
!!$   endif
!!$
!!$   if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
!!$      var_nest(isn_c:ien_c, jsn_c:jen_c,:) = nbuffer(isn_c:ien_c, jsn_c:jen_c,:)
!!$   endif

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

 subroutine nested_grid_BC_save_mpp(nest_domain, ind, wt, istag, jstag, &
      npx, npy, npz, bd, nest_BC, &
      proc_in, ns_in)

   type(fv_grid_bounds_type), intent(IN) :: bd
   type(nest_domain_type), intent(INOUT) :: nest_domain
   integer, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,2), intent(IN) :: ind
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,4), intent(IN) :: wt
   integer, intent(IN) :: istag, jstag, npx, npy, npz
   integer, intent(IN), OPTIONAL :: ns_in
   logical, intent(IN), OPTIONAL :: proc_in

   !!NOTE: if declaring an ALLOCATABLE array with intent(OUT), the resulting dummy array
   !!      will NOT be allocated! This goes for allocatable members of derived types as well.
   type(fv_nest_BC_type_3d), intent(INOUT), target :: nest_BC
   
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,npz) :: var_coarse_dummy

   real, dimension(:,:,:), pointer :: var_east, var_west, var_south, var_north
   !real, dimension(bd%ie+1+istag-ns:bd%ied+istag,bd%jsd:bd%jed+jstag,npz), intent(out) :: var_east
   !real, dimension(bd%isd:bd%is-1+ns,bd%jsd:bd%jed+jstag,npz), intent(out) :: var_west
   !real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%js-1+ns,npz), intent(out) :: var_south
   !real, dimension(bd%isd:bd%ied+istag,bd%je+1+jstag-ns:bd%jed+jstag,npz), intent(out) :: var_north

   integer                      :: position

   integer                      :: isw_f, iew_f, jsw_f, jew_f, isw_c, iew_c, jsw_c, jew_c
   integer                      :: ise_f, iee_f, jse_f, jee_f, ise_c, iee_c, jse_c, jee_c
   integer                      :: iss_f, ies_f, jss_f, jes_f, iss_c, ies_c, jss_c, jes_c
   integer                      :: isn_f, ien_f, jsn_f, jen_f, isn_c, ien_c, jsn_c, jen_c
   real,    allocatable         :: wbuffer(:,:,:)
   real,    allocatable         :: ebuffer(:,:,:)
   real,    allocatable         :: sbuffer(:,:,:)
   real,    allocatable         :: nbuffer(:,:,:)

   integer :: i,j, k, ic, jc, istart, iend
   integer :: ns = 0
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

   if (present(proc_in)) then
      process = proc_in
   else
      process = .true.
   endif

   if (present(ns_in)) then
      ns = ns_in
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
        WEST,  position=position)
   call mpp_get_C2F_index(nest_domain, ise_f, iee_f, jse_f, jee_f, ise_c, iee_c, jse_c, jee_c, &
        EAST,  position=position)
   call mpp_get_C2F_index(nest_domain, iss_f, ies_f, jss_f, jes_f, iss_c, ies_c, jss_c, jes_c, &
        SOUTH,  position=position)
   call mpp_get_C2F_index(nest_domain, isn_f, ien_f, jsn_f, jen_f, isn_c, ien_c, jsn_c, jen_c, &
        NORTH,  position=position)

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
   call mpp_update_nest_fine(var_coarse_dummy, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer,  position=position)
       call timing_off('COMM_TOTAL')

   if (process) then

      var_east  => nest_BC%east_t1
      var_west  => nest_BC%west_t1
      var_north => nest_BC%north_t1
      var_south => nest_BC%south_t1

   ! ?buffer has uninterpolated coarse-grid data; need to perform interpolation ourselves
   !To do this more securely, instead of using is/etc we could use the fine-grid indices defined above
   if (is == 1  ) then
      do k=1,npz
      do j=jsd,jed+jstag
         do i=isd,ns

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            if (ic < isw_c .or. ic > iew_c) then
               print*, ic, isw_c, iew_c, is
               call mpp_error(FATAL, 'wbuffer i index out of bounds')
            end if

            if (jc < jsw_c .or. jc > jew_c) then
               print*, jc, jsw_c, jew_c
               call mpp_error(FATAL, 'wbuffer j index out of bounds')
            end if

            var_west(i,j,k) = &
                 wt(i,j,1)*wbuffer(ic,  jc,k) +  &
                 wt(i,j,2)*wbuffer(ic,  jc+1,k) +  &
                 wt(i,j,3)*wbuffer(ic+1,jc+1,k) +  &
                 wt(i,j,4)*wbuffer(ic+1,jc,k) 

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

      do k=1,npz
      do j=jsd,ns
         do i=istart,iend+istag

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            if (ic < iss_c .or. ic > ies_c) then
               print*, ic, iss_c, ies_c, istart, iend
               call mpp_error(FATAL, 'sbuffer i index out of bounds')
            end if

            if (jc < jss_c .or. jc > jes_c) then
               print*, jc, jss_c, jes_c
               call mpp_error(FATAL, 'sbuffer j index out of bounds')
            end if

            var_south(i,j,k) = &
                 wt(i,j,1)*sbuffer(ic,  jc,k) +  &
                 wt(i,j,2)*sbuffer(ic,  jc+1,k) +  &
                 wt(i,j,3)*sbuffer(ic+1,jc+1,k) +  &
                 wt(i,j,4)*sbuffer(ic+1,jc,k) 

         end do
      end do
      end do
   end if


   if (ie == npx-1 ) then

      do k=1,npz
      do j=jsd,jed+jstag
         do i=npx+istag-ns,ied+istag

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            if (ic < ise_c .or. ic > iee_c) then
               print*, ic, ise_c, iee_c, ie
               call mpp_error(FATAL, 'ebuffer i index out of bounds')
            end if

            if (jc < jse_c .or. jc > jee_c) then
               print*, jc, jse_c, jee_c
               call mpp_error(FATAL, 'ebuffer j index out of bounds')
            end if

            var_east(i,j,k) = &
                 wt(i,j,1)*ebuffer(ic,  jc,k) +  &
                 wt(i,j,2)*ebuffer(ic,  jc+1,k) +  &
                 wt(i,j,3)*ebuffer(ic+1,jc+1,k) +  &
                 wt(i,j,4)*ebuffer(ic+1,jc,k) 

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


      do k=1,npz
      do j=npy+jstag-ns,jed+jstag
         do i=istart,iend+istag

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            if (ic < isn_c .or. ic > ien_c) then
               print*, ic, isn_c, ien_c, istart, iend
               call mpp_error(FATAL, 'nbuffer i index out of bounds')
            end if

            if (jc < jsn_c .or. jc > jen_c) then
               print*, jc, jsn_c, jen_c
               call mpp_error(FATAL, 'nbuffer j index out of bounds')
            end if

            var_north(i,j,k) = &
                 wt(i,j,1)*nbuffer(ic,  jc,k) +  &
                 wt(i,j,2)*nbuffer(ic,  jc+1,k) +  &
                 wt(i,j,3)*nbuffer(ic+1,jc+1,k) +  &
                 wt(i,j,4)*nbuffer(ic+1,jc,k) 

         end do
      end do
      end do
   end if

  endif !process

   deallocate(wbuffer)
   deallocate(ebuffer)
   deallocate(sbuffer)
   deallocate(nbuffer)


 end subroutine nested_grid_BC_save_mpp

 subroutine nested_grid_BC_save_send(var_coarse, nest_domain, istag, jstag, isg, ieg, jsg, jeg, npz)

   real, dimension(isg:ieg+istag,jsg:jeg+jstag,npz), intent(IN) :: var_coarse
   type(nest_domain_type), intent(INOUT) :: nest_domain
   integer, intent(IN) :: istag, jstag, isg, ieg, jsg, jeg, npz

   integer                      :: position

   real,    allocatable         :: wbuffer(:,:,:)
   real,    allocatable         :: ebuffer(:,:,:)
   real,    allocatable         :: sbuffer(:,:,:)
   real,    allocatable         :: nbuffer(:,:,:)


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
   wbuffer = 0

   allocate(ebuffer(1,1,1))
   ebuffer = 0

   allocate(sbuffer(1,1,1))
   sbuffer = 0

   allocate(nbuffer(1,1,1))
   nbuffer = 0


       call timing_on ('COMM_TOTAL')
!!$       print*, mpp_pe, 'CALLING UPDATE_NEST_FINE'
   call mpp_update_nest_fine(var_coarse, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer,  position=position)
       call timing_off('COMM_TOTAL')

   deallocate(wbuffer)
   deallocate(ebuffer)
   deallocate(sbuffer)
   deallocate(nbuffer)


 end subroutine nested_grid_BC_save_send

 subroutine nested_grid_BC_save_2D(var_coarse, ind, wt, istag, jstag, &
      npx, npy, bd, isg, ieg, jsg, jeg, var_east, var_west, var_north, var_south, ns)

   type(fv_grid_bounds_type), intent(IN) :: bd
   real, dimension(isg:ieg+istag,jsg:jeg+jstag), intent(IN) :: var_coarse
   integer, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,2), intent(IN) :: ind
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,4), intent(IN) :: wt
   integer, intent(IN) :: istag, jstag, npx, npy, isg, ieg, jsg, jeg
   integer, intent(IN) :: ns
   
   real, dimension(bd%ie+1+istag-ns:bd%ied+istag,bd%jsd:bd%jed+jstag), intent(out) :: var_east
   real, dimension(bd%isd:bd%is-1+ns,bd%jsd:bd%jed+jstag), intent(out) :: var_west
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%js-1+ns), intent(out) :: var_south
   real, dimension(bd%isd:bd%ied+istag,bd%je+1+jstag-ns:bd%jed+jstag), intent(out) :: var_north

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

   if (is == 1  ) then
      do j=jsd,jed+jstag
         do i=isd,ns

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            var_west(i,j) = &
                 wt(i,j,1)*var_coarse(ic,  jc) +  &
                 wt(i,j,2)*var_coarse(ic,  jc+1) +  &
                 wt(i,j,3)*var_coarse(ic+1,jc+1) +  &
                 wt(i,j,4)*var_coarse(ic+1,jc) 

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

      do j=jsd,ns
         do i=istart,iend+istag

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            var_south(i,j) = &
                 wt(i,j,1)*var_coarse(ic,  jc) +  &
                 wt(i,j,2)*var_coarse(ic,  jc+1) +  &
                 wt(i,j,3)*var_coarse(ic+1,jc+1) +  &
                 wt(i,j,4)*var_coarse(ic+1,jc) 

         end do
      end do
   end if


   if (ie == npx-1 ) then
      do j=jsd,jed+jstag
         do i=npx+istag-ns,ied+istag

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            var_east(i,j) = &
                 wt(i,j,1)*var_coarse(ic,  jc) +  &
                 wt(i,j,2)*var_coarse(ic,  jc+1) +  &
                 wt(i,j,3)*var_coarse(ic+1,jc+1) +  &
                 wt(i,j,4)*var_coarse(ic+1,jc) 

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


      do j=npy+jstag-ns,jed+jstag
         do i=istart,iend+istag

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            var_north(i,j) = &
                 wt(i,j,1)*var_coarse(ic,  jc) +  &
                 wt(i,j,2)*var_coarse(ic,  jc+1) +  &
                 wt(i,j,3)*var_coarse(ic+1,jc+1) +  &
                 wt(i,j,4)*var_coarse(ic+1,jc) 

         end do
      end do
   end if



 end subroutine nested_grid_BC_save_2D

 subroutine nested_grid_BC_save_3D(var_coarse, ind, wt, istag, jstag, &
      npx, npy, npz, bd, isg, ieg, jsg, jeg, var_east, var_west, var_north, var_south, ns)

   type(fv_grid_bounds_type), intent(IN) :: bd
   real, dimension(isg:ieg+istag,jsg:jeg+jstag,npz), intent(IN) :: var_coarse
   integer, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,2), intent(IN) :: ind
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,4), intent(IN) :: wt
   integer, intent(IN) :: istag, jstag, npx, npy, isg, ieg, jsg, jeg, npz
   integer, intent(IN) :: ns
   
   real, dimension(bd%ie+1+istag-ns:bd%ied+istag,bd%jsd:bd%jed+jstag,npz), intent(out) :: var_east
   real, dimension(bd%isd:bd%is-1+ns,bd%jsd:bd%jed+jstag,npz), intent(out) :: var_west
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%js-1+ns,npz), intent(out) :: var_south
   real, dimension(bd%isd:bd%ied+istag,bd%je+1+jstag-ns:bd%jed+jstag,npz), intent(out) :: var_north

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

   if (is == 1  ) then
      do k=1,npz
      do j=jsd,jed+jstag
         do i=isd,ns

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            var_west(i,j,k) = &
                 wt(i,j,1)*var_coarse(ic,  jc,  k) +  &
                 wt(i,j,2)*var_coarse(ic,  jc+1,k) +  &
                 wt(i,j,3)*var_coarse(ic+1,jc+1,k) +  &
                 wt(i,j,4)*var_coarse(ic+1,jc,  k) 

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

      do k=1,npz
      do j=jsd,ns
         do i=istart,iend+istag

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            var_south(i,j,k) = &
                 wt(i,j,1)*var_coarse(ic,  jc,  k) +  &
                 wt(i,j,2)*var_coarse(ic,  jc+1,k) +  &
                 wt(i,j,3)*var_coarse(ic+1,jc+1,k) +  &
                 wt(i,j,4)*var_coarse(ic+1,jc,  k) 

         end do
      end do
      end do
   end if


   if (ie == npx-1 ) then
      do k=1,npz
      do j=jsd,jed+jstag
         do i=npx+istag-ns,ied+istag

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            var_east(i,j,k) = &
                 wt(i,j,1)*var_coarse(ic,  jc,  k) +  &
                 wt(i,j,2)*var_coarse(ic,  jc+1,k) +  &
                 wt(i,j,3)*var_coarse(ic+1,jc+1,k) +  &
                 wt(i,j,4)*var_coarse(ic+1,jc,  k) 

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

      do k=1,npz
      do j=npy+jstag-ns,jed+jstag
         do i=istart,iend+istag

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            var_north(i,j,k) = &
                 wt(i,j,1)*var_coarse(ic,  jc,  k) +  &
                 wt(i,j,2)*var_coarse(ic,  jc+1,k) +  &
                 wt(i,j,3)*var_coarse(ic+1,jc+1,k) +  &
                 wt(i,j,4)*var_coarse(ic+1,jc,  k) 

         end do
      end do
      end do
   end if



 end subroutine nested_grid_BC_save_3D


  ! A NOTE ON BCTYPE: currently only an interpolation BC is implemented, although there
  ! is an option to use sponge points in the interior. bctype >= 2 currently correspond
  ! to a flux BC on the tracers ONLY, which is implemented in fv_tracer.
 
 !This routine only works for integer interpolation, and applies an
 ! interative algorithm that does not allow skipped steps. If direct
 ! interpolation to a particular time (not necessarily an integer
 ! number of small timesteps) is desired, and two timelevels are
 ! available, then use nested_grid_BC_intT; if only one time level
 ! is available use nested_grid_BC.

 subroutine nested_grid_BC_apply_3D(var_nest, istag, jstag, &
      npx, npy, npz, bd, step, split, var_east, var_west, var_north, var_south, bctype, nsponge, s_weight)

   type(fv_grid_bounds_type), intent(IN) :: bd
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag, npz), intent(INOUT) :: var_nest
   integer, intent(IN) :: istag, jstag, npx, npy, npz
   integer, intent(IN) :: split, step
   integer, intent(IN) :: bctype, nsponge
   real, intent(IN) :: s_weight
   
   real, dimension(bd%ie+1+istag-nsponge:bd%ied+istag,bd%jsd:bd%jed+jstag,npz), intent(in) :: var_east
   real, dimension(bd%isd:bd%is-1+nsponge,bd%jsd:bd%jed+jstag,npz), intent(in) :: var_west
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%js-1+nsponge,npz), intent(in) :: var_south
   real, dimension(bd%isd:bd%ied+istag,bd%je+1+jstag-nsponge:bd%jed+jstag,npz), intent(in) :: var_north

   integer :: i,j, istart, iend, k
   real :: denom
   real :: weights(nsponge)

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

   denom = 1./real(split - step + 1)

   do i=1,nsponge
      weights(i) = s_weight*(1. + real(nsponge - i))/real(nsponge)
   end do

   if (printdiag .and. is_master() .and. nsponge > 0) then
      print*, 'SPONGE WEIGHTS'
      do i=1,nsponge
         write(*,'(I2, E14.8)') i, weights(i)
      end do
      printdiag = .false.
   end if

   if (is == 1  ) then
      do k=1,npz
      do j=jsd,jed+jstag
         do i=isd,0

            var_nest(i,j,k) = (var_nest(i,j,k)*real(split-step) + var_west(i,j,k))*denom

         end do
      end do
      end do

      if (nsponge > 0) then
         do k=1,npz
         do j=js,je+jstag
         do i=1,nsponge
            var_nest(i,j,k) = var_nest(i,j,k) - weights(i)*( var_nest(i,j,k) - var_west(i,j,k) )
         end do
         end do
         end do
      end if

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

      do k=1,npz
      do j=jsd,0
         do i=istart,iend+istag

            var_nest(i,j,k) = (var_nest(i,j,k)*real(split-step) + var_south(i,j,k))*denom

         end do
      end do
      end do

      if (nsponge > 0) then
         if (is == 1)     istart = is + nsponge
         if (ie == npx-1) iend   = ie - nsponge

         do k=1,npz
         do j=1,nsponge
         do i=istart,iend+istag
            var_nest(i,j,k) = var_nest(i,j,k) - weights(j)*( var_nest(i,j,k) - var_south(i,j,k) )
         end do
         end do
         end do
      end if

   end if


   if (ie == npx-1 ) then
      do k=1,npz
      do j=jsd,jed+jstag
         do i=npx+istag,ied+istag

            var_nest(i,j,k) = (var_nest(i,j,k)*real(split-step) + var_east(i,j,k))*denom

         end do
      end do
      end do

      if (nsponge > 0) then
         do k=1,npz
         do j=js,je+jstag
         do i=npx+istag-nsponge,npx-1+istag
            var_nest(i,j,k) = var_nest(i,j,k) - weights(npx+istag-i)*( var_nest(i,j,k) - var_east(i,j,k) )
         end do
         end do
         end do
      end if

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

      do k=1,npz
      do j=npy+jstag,jed+jstag
         do i=istart,iend+istag

            var_nest(i,j,k) = (var_nest(i,j,k)*real(split-step) + var_north(i,j,k))*denom

         end do
      end do
      end do

      if (nsponge > 0) then
         if (is == 1)     istart = is + nsponge
         if (ie == npx-1) iend   = ie - nsponge
         
         do k=1,npz
         do j=npy+jstag-nsponge,npy-1+jstag
         do i=istart,iend+istag
            var_nest(i,j,k) = var_nest(i,j,k) - weights(npy+jstag-j)*( var_nest(i,j,k) - var_north(i,j,k) )
         end do
         end do
         end do
      end if
   end if


 end subroutine nested_grid_BC_apply_3D

 subroutine nested_grid_BC_apply_2D(var_nest, istag, jstag, &
      npx, npy, bd, step, split, var_east, var_west, var_north, var_south, bctype, nsponge, s_weight)

   type(fv_grid_bounds_type), intent(IN) :: bd
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag), intent(INOUT) :: var_nest
   integer, intent(IN) :: istag, jstag, npx, npy
   integer, intent(IN) :: split, step
   integer, intent(IN) :: bctype, nsponge
   real, intent(IN) :: s_weight
   
   real, dimension(bd%ie+1+istag-nsponge:bd%ied+istag,bd%jsd:bd%jed+jstag), intent(in) :: var_east
   real, dimension(bd%isd:bd%is-1+nsponge,bd%jsd:bd%jed+jstag), intent(in) :: var_west
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%js-1+nsponge), intent(in) :: var_south
   real, dimension(bd%isd:bd%ied+istag,bd%je+1+jstag-nsponge:bd%jed+jstag), intent(in) :: var_north

   integer :: i,j, istart, iend
   real :: denom
   real :: weights(nsponge)

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

   denom = 1./real(split - step + 1)

   do i=1,nsponge
      weights(i) = s_weight*(1. + real(nsponge - i))/real(nsponge)
   end do

   if (printdiag .and. is_master() .and. nsponge > 0) then
      print*, 'SPONGE WEIGHTS'
      do i=1,nsponge
         write(*,'(I2, E12.8)') i, weights(i)
      end do
      printdiag = .false.
   end if

   if (is == 1  ) then
      do j=jsd,jed+jstag
         do i=isd,0

            var_nest(i,j) = (var_nest(i,j)*real(split-step) + var_west(i,j))*denom

         end do
      end do

      if (nsponge > 0) then
         do j=js,je+jstag
         do i=1,nsponge
            var_nest(i,j) = var_nest(i,j) - weights(i)*( var_nest(i,j) - var_west(i,j) )
         end do
         end do
      end if

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

      do j=jsd,0
         do i=istart,iend+istag

            var_nest(i,j) = (var_nest(i,j)*real(split-step) + var_south(i,j))*denom

         end do
      end do

      if (nsponge > 0) then
         if (is == 1)     istart = is + nsponge
         if (ie == npx-1) iend   = ie - nsponge

         do j=1,nsponge
         do i=istart,iend+istag
            var_nest(i,j) = var_nest(i,j) - weights(j)*( var_nest(i,j) - var_south(i,j) )
         end do
         end do
      end if

   end if


   if (ie == npx-1 ) then
      do j=jsd,jed+jstag
         do i=npx+istag,ied+istag

            var_nest(i,j) = (var_nest(i,j)*real(split-step) + var_east(i,j))*denom

         end do
      end do

      if (nsponge > 0) then
         do j=js,je+jstag
         do i=npx+istag-nsponge,npx-1+istag
            var_nest(i,j) = var_nest(i,j) - weights(npx+istag-i)*( var_nest(i,j) - var_east(i,j) )
         end do
         end do
      end if

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


      do j=npy+jstag,jed+jstag
         do i=istart,iend+istag

            var_nest(i,j) = (var_nest(i,j)*real(split-step) + var_north(i,j))*denom

         end do
      end do

      if (nsponge > 0) then
         if (is == 1)     istart = is + nsponge
         if (ie == npx-1) iend   = ie - nsponge
         
         do j=npy+jstag-nsponge,npy-1+jstag
         do i=istart,iend+istag
            var_nest(i,j) = var_nest(i,j) - weights(npy+jstag-j)*( var_nest(i,j) - var_north(i,j) )
         end do
         end do
      end if
   end if


 end subroutine nested_grid_BC_apply_2D


 subroutine nested_grid_BC_apply_intT_3D(var_nest, istag, jstag, &
      npx, npy, npz, bd, step, split, &
      var_east_t0, var_west_t0, var_north_t0, var_south_t0, &
      var_east_t1, var_west_t1, var_north_t1, var_south_t1, &
      bctype, nsponge, s_weight)

   type(fv_grid_bounds_type), intent(IN) :: bd
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag, npz), intent(INOUT) :: var_nest
   integer, intent(IN) :: istag, jstag, npx, npy, npz
   real, intent(IN) :: split, step
   integer, intent(IN) :: bctype, nsponge
   real, intent(IN) :: s_weight
   
   real, dimension(bd%ie+1+istag-nsponge:bd%ied+istag,bd%jsd:bd%jed+jstag,npz), intent(in) :: var_east_t0
   real, dimension(bd%isd:bd%is-1+nsponge,bd%jsd:bd%jed+jstag,npz), intent(in) :: var_west_t0
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%js-1+nsponge,npz), intent(in) :: var_south_t0
   real, dimension(bd%isd:bd%ied+istag,bd%je+1+jstag-nsponge:bd%jed+jstag,npz), intent(in) :: var_north_t0

   real, dimension(bd%ie+1+istag-nsponge:bd%ied+istag,bd%jsd:bd%jed+jstag,npz), intent(in) :: var_east_t1
   real, dimension(bd%isd:bd%is-1+nsponge,bd%jsd:bd%jed+jstag,npz), intent(in) :: var_west_t1
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%js-1+nsponge,npz), intent(in) :: var_south_t1
   real, dimension(bd%isd:bd%ied+istag,bd%je+1+jstag-nsponge:bd%jed+jstag,npz), intent(in) :: var_north_t1

   integer :: i,j, istart, iend, k
   real :: denom
   real :: weights(nsponge)

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

   do i=1,nsponge
      weights(i) = s_weight*(1. + real(nsponge - i))/real(nsponge)
   end do

   if (is == 1  ) then
      do k=1,npz
      do j=jsd,jed+jstag
         do i=isd,0

            var_nest(i,j,k) = (var_west_t0(i,j,k)*(split-step) + step*var_west_t1(i,j,k))*denom

         end do
      end do
      end do

      if (nsponge > 0) then
         do k=1,npz
         do j=js,je+jstag
         do i=1,nsponge
            var_nest(i,j,k) = var_nest(i,j,k) - weights(i)*( var_nest(i,j,k) - (var_west_t0(i,j,k)*(split-step) + step*var_west_t1(i,j,k))*denom )
         end do
         end do
         end do
      end if

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

      do k=1,npz
      do j=jsd,0
         do i=istart,iend+istag

            var_nest(i,j,k) = (var_south_t0(i,j,k)*(split-step) + step*var_south_t1(i,j,k))*denom

         end do
      end do
      end do

      if (nsponge > 0) then
         if (is == 1)     istart = is + nsponge
         if (ie == npx-1) iend   = ie - nsponge

         do k=1,npz
         do j=1,nsponge
         do i=istart,iend+istag
            var_nest(i,j,k) = var_nest(i,j,k) - weights(j)*( var_nest(i,j,k) - (var_south_t0(i,j,k)*(split-step) + step*var_south_t1(i,j,k))*denom )
         end do
         end do
         end do
      end if

   end if


   if (ie == npx-1 ) then
      do k=1,npz
      do j=jsd,jed+jstag
         do i=npx+istag,ied+istag

            var_nest(i,j,k) = (var_east_t0(i,j,k)*(split-step) + step*var_east_t1(i,j,k))*denom

         end do
      end do
      end do

      if (nsponge > 0) then
         do k=1,npz
         do j=js,je+jstag
         do i=npx+istag-nsponge,npx-1+istag
            var_nest(i,j,k) = var_nest(i,j,k) - weights(npx+istag-i)*( var_nest(i,j,k) - (var_east_t0(i,j,k)*(split-step) + step*var_east_t1(i,j,k))*denom )
         end do
         end do
         end do
      end if

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

      do k=1,npz
      do j=npy+jstag,jed+jstag
         do i=istart,iend+istag

            var_nest(i,j,k) = (var_north_t0(i,j,k)*(split-step) + step*var_north_t1(i,j,k))*denom

         end do
      end do
      end do

      if (nsponge > 0) then
         if (is == 1)     istart = is + nsponge
         if (ie == npx-1) iend   = ie - nsponge
         
         do k=1,npz
         do j=npy+jstag-nsponge,npy-1+jstag
         do i=istart,iend+istag
            var_nest(i,j,k) = var_nest(i,j,k) - weights(npy+jstag-j)*( var_nest(i,j,k) - (var_north_t0(i,j,k)*(split-step) + step*var_north_t1(i,j,k))*denom )
         end do
         end do
         end do
      end if
   end if


 end subroutine nested_grid_BC_apply_intT_3D
 
 subroutine nested_grid_BC_apply_intT_2D(var_nest, istag, jstag, &
      npx, npy, bd, step, split, &
      var_east_t0, var_west_t0, var_north_t0, var_south_t0, &
      var_east_t1, var_west_t1, var_north_t1, var_south_t1, &
      bctype, nsponge, s_weight)

   type(fv_grid_bounds_type), intent(IN) :: bd
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag), intent(INOUT) :: var_nest
   integer, intent(IN) :: istag, jstag, npx, npy
   real, intent(IN) :: split, step
   integer, intent(IN) :: bctype, nsponge
   real, intent(IN) :: s_weight
   
   real, dimension(bd%ie+1+istag-nsponge:bd%ied+istag,bd%jsd:bd%jed+jstag), intent(in) :: var_east_t0
   real, dimension(bd%isd:bd%is-1+nsponge,bd%jsd:bd%jed+jstag), intent(in) :: var_west_t0
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%js-1+nsponge), intent(in) :: var_south_t0
   real, dimension(bd%isd:bd%ied+istag,bd%je+1+jstag-nsponge:bd%jed+jstag), intent(in) :: var_north_t0

   real, dimension(bd%ie+1+istag-nsponge:bd%ied+istag,bd%jsd:bd%jed+jstag), intent(in) :: var_east_t1
   real, dimension(bd%isd:bd%is-1+nsponge,bd%jsd:bd%jed+jstag), intent(in) :: var_west_t1
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%js-1+nsponge), intent(in) :: var_south_t1
   real, dimension(bd%isd:bd%ied+istag,bd%je+1+jstag-nsponge:bd%jed+jstag), intent(in) :: var_north_t1

   integer :: i,j, istart, iend
   real :: denom
   real :: weights(nsponge)

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

   do i=1,nsponge
      weights(i) = s_weight*(1. + real(nsponge - i))/real(nsponge)
   end do

   if (is == 1  ) then
      do j=jsd,jed+jstag
         do i=isd,0

            var_nest(i,j) = (var_west_t0(i,j)*(split-step) + step*var_west_t1(i,j))*denom

         end do
      end do

      if (nsponge > 0) then
         do j=js,je+jstag
         do i=1,nsponge
            var_nest(i,j) = var_nest(i,j) - weights(i)*( var_nest(i,j) - (var_west_t0(i,j)*(split-step) + step*var_west_t1(i,j))*denom )
         end do
         end do
      end if

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

      do j=jsd,0
         do i=istart,iend+istag

            var_nest(i,j) = (var_south_t0(i,j)*(split-step) + step*var_south_t1(i,j))*denom

         end do
      end do

      if (nsponge > 0) then
         if (is == 1)     istart = is + nsponge
         if (ie == npx-1) iend   = ie - nsponge

         do j=1,nsponge
         do i=istart,iend+istag
            var_nest(i,j) = var_nest(i,j) - weights(j)*( var_nest(i,j) - (var_south_t0(i,j)*(split-step) + step*var_south_t1(i,j))*denom )
         end do
         end do
      end if

   end if


   if (ie == npx-1 ) then
      do j=jsd,jed+jstag
         do i=npx+istag,ied+istag

            var_nest(i,j) = (var_east_t0(i,j)*(split-step) + step*var_east_t1(i,j))*denom

         end do
      end do

      if (nsponge > 0) then
         do j=js,je+jstag
         do i=npx+istag-nsponge,npx-1+istag
            var_nest(i,j) = var_nest(i,j) - weights(npx+istag-i)*( var_nest(i,j) - (var_east_t0(i,j)*(split-step) + step*var_east_t1(i,j))*denom  )
         end do
         end do
      end if

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


      do j=npy+jstag,jed+jstag
         do i=istart,iend+istag

            var_nest(i,j) = (var_north_t0(i,j)*(split-step) + step*var_north_t1(i,j))*denom

         end do
      end do

      if (nsponge > 0) then
         if (is == 1)     istart = is + nsponge
         if (ie == npx-1) iend   = ie - nsponge
         
         do j=npy+jstag-nsponge,npy-1+jstag
         do i=istart,iend+istag
            var_nest(i,j) = var_nest(i,j) - weights(npy+jstag-j)*( var_nest(i,j) - (var_north_t0(i,j)*(split-step) + step*var_north_t1(i,j))*denom )
         end do
         end do
      end if
   end if


 end subroutine nested_grid_BC_apply_intT_2D


 subroutine nested_grid_sponge_apply(var_nest, istag, jstag, &
      npx, npy, bd, step, split, var_east, var_west, var_north, var_south, nsponge)
   !sponge BC: just relaxing to the advanced (not time-interpolated) coarse-grid solution.

   type(fv_grid_bounds_type), intent(IN) :: bd
   real, dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag), intent(INOUT) :: var_nest
   integer, intent(IN) :: istag, jstag, npx, npy
   integer, intent(IN) :: split, step
   integer, intent(in) :: nsponge
   
   real, dimension(npx-nsponge:bd%ied+istag,bd%jsd:bd%jed+jstag)        , intent(in) :: var_east
   real, dimension(bd%isd:nsponge          ,bd%jsd:bd%jed+jstag)        , intent(in) :: var_west
   real, dimension(bd%isd:bd%ied+istag        ,bd%jsd:nsponge)          , intent(in) :: var_south
   real, dimension(bd%isd:bd%ied+istag        ,npy-nsponge:bd%jed+jstag), intent(in) :: var_north

   integer :: i,j, istart, iend, bctype
   real, parameter :: s_weight = 0.1
   real :: denom
   real :: weights(nsponge)

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

   if (nsponge == 0) return

   if (nsponge <= 0) call mpp_error(FATAL, 'nsponge must be non-negative')

   do i=1,nsponge
      weights(i) = s_weight*(1. + real(nsponge - i))/real(nsponge)
   end do

   if (printdiag .and. is_master() .and. nsponge > 0) then
      print*, 'SPONGE WEIGHTS'
      do i=1,nsponge
         write(*,'(I2, E14.8)') i, weights(i)
      end do
   end if
   printdiag = .false.

   if (is == 1) then

      do j=js,je+jstag
         do i=1,nsponge
            var_nest(i,j) = var_nest(i,j) - weights(i)*( var_nest(i,j) - var_west(i,j) )
         end do
      end do
   end if

   if (js == 1) then

      if (is == 1) then
         istart = is+nsponge  
      else
         istart = is
      end if

      if (ie == npx-1) then
         iend = ie-nsponge
      else
         iend = ie
      end if

      do j=1,nsponge
         do i=istart,iend+istag
            var_nest(i,j) = var_nest(i,j) - weights(j)*( var_nest(i,j) - var_south(i,j) )
         end do
      end do

   end if


   if (ie == npx-1) then

      do j=js,je+jstag
         do i=npx+istag-nsponge,npx-1+istag
            var_nest(i,j) = var_nest(i,j) - weights(npx+istag-i)*( var_nest(i,j) - var_east(i,j) )
         end do
      end do

   end if

   if (je == npy-1 ) then

      if (is == 1) then
         istart = is+nsponge 
      else
         istart = is
      end if

      if (ie == npx-1) then
         iend = ie-nsponge
      else
         iend = ie
      end if

      do j=npy+jstag-nsponge,npy-1+jstag
         do i=istart,iend+istag
            var_nest(i,j) = var_nest(i,j) - weights(npy+jstag-j)*( var_nest(i,j) - var_north(i,j) )
         end do
      end do

   end if


 end subroutine nested_grid_sponge_apply
 
 !!! CLEANUP: remove concurrent argument?
 subroutine gather_grid_2d(parent_grid, var_coarse_this_grid, var_coarse, &
      isd_p, ied_p, jsd_p, jed_p, isg, ieg, jsg, jeg, istag, jstag, &
      parent_tile, concurrent, recv_pes, recv_npes)

   type(fv_atmos_type), intent(in) :: parent_grid
   real, dimension(isd_p:ied_p+istag,jsd_p:jed_p+jstag), intent(IN) :: var_coarse_this_grid
   real, dimension(isg:ieg+istag,jsg:jeg+jstag,1), intent(out) :: var_coarse
   integer, intent(IN) :: isd_p, ied_p, jsd_p, jed_p, isg, ieg, jsg, jeg, &
        istag, jstag, parent_tile, recv_npes, recv_pes(recv_npes)
   logical, intent(IN) :: concurrent

   integer :: position, npx, npy, n, p, sending_proc

   integer, allocatable, dimension(:) :: tilelist
   
   if (istag == 0 .and. jstag == 0) position=CENTER
   if (istag == 1 .and. jstag == 0) position=EAST
   if (istag == 0 .and. jstag == 1) position=NORTH
   if (istag == 1 .and. jstag == 1) position=CORNER


   if (concurrent) then

      sending_proc = parent_grid%pelist(1) + (parent_tile-1)*parent_grid%npes_per_tile

       !Call gather grid on the procs that have the required data.
       !Then broadcast from the head PE to the receiving PEs
       if (ANY(parent_grid%pelist == mpp_pe()) .and. parent_tile == parent_grid%tile) then
          call mpp_global_field( &
               parent_grid%domain, &
               var_coarse_this_grid, var_coarse(:,:,1), position=position)
          if (mpp_pe() == sending_proc) then 
             do p=1,recv_npes
                call mpp_send(var_coarse,size(var_coarse),recv_pes(p))
             enddo
          endif
       endif
       if (ANY(recv_pes == mpp_pe())) then
          call mpp_recv(var_coarse, size(var_coarse), sending_proc)
       endif

   else
      if (parent_grid%tile == parent_tile) then
         call mpp_global_field( &
              parent_grid%domain, &
              var_coarse_this_grid, var_coarse(:,:,1), position=position)
      end if

      npx = ieg-isg+1+istag
      npy = jeg-jsg+1+jstag

      allocate(tilelist(0:mpp_npes()-1))
      tilelist = 0
      tilelist(mpp_pe()) = parent_grid%tile
      call mpp_sum(tilelist,mpp_npes())
!!$   if (gid == 0) then
!!$      print*, 'tilelist:'
!!$      do n=0,mpp_npes()-1
!!$         print*, n, tilelist(n)
!!$      end do
!!$   end if
!!$   call mpp_sync

      !Don't need to do communication to the other tiles if there is only one tile
      if (maxval(tilelist) == minval(tilelist)) return

      do n=0,mpp_npes()-1
         if (tilelist(n) == parent_tile) then
            if (mpp_pe() == n) then
               !            print*, 'gather_grid: proc ', gid, ' broadcasting '
               !            print*, 'gather_grid: SIZE: ', size(var_coarse), npx, npy
               call mpp_broadcast(var_coarse, npx*npy, n, parent_grid%pelist)
            else 
               !            print*, 'gather_grid: proc ', gid, ' receiving from proc ', n
               call mpp_broadcast(var_coarse, npx*npy, n, parent_grid%pelist)                         
            end if
            exit
         end if
      end do

   endif

 end subroutine gather_grid_2d

 subroutine gather_grid_3d(parent_grid, var_coarse_this_grid, var_coarse, &
      isd_p, ied_p, jsd_p, jed_p, isg, ieg, jsg, jeg, npz, istag, jstag, &
      parent_tile, concurrent, recv_pes, recv_npes)

   type(fv_atmos_type), intent(in) :: parent_grid
   real, dimension(isg:ieg+istag,jsg:jeg+jstag,npz), intent(out) :: var_coarse
   real, dimension(isd_p:ied_p+istag,jsd_p:jed_p+jstag,npz), intent(IN) :: var_coarse_this_grid
   integer, intent(IN) :: isd_p, ied_p, jsd_p, jed_p, isg, ieg, jsg, jeg, &
        istag, jstag, npz, parent_tile, recv_npes, recv_pes(recv_npes)
   logical, intent(IN) :: concurrent

   integer :: position, npx, npy, n, k, sending_proc, p

   integer, allocatable, dimension(:) :: tilelist

   logical, parameter :: z_nest = .false.

   if (z_nest) then
      do k=1,npz
         call gather_grid_2d(parent_grid, var_coarse_this_grid(:,:,k), var_coarse(:,:,k), &
              isd_p, ied_p, jsd_p, jed_p, isg, ieg, jsg, jeg, istag, jstag, parent_tile, &
              concurrent, recv_pes, recv_npes)
      end do
      return
   end if
   
   if (istag == 0 .and. jstag == 0) position=CENTER
   if (istag == 1 .and. jstag == 0) position=EAST
   if (istag == 0 .and. jstag == 1) position=NORTH
   if (istag == 1 .and. jstag == 1) position=CORNER



   if (concurrent) then

      sending_proc = parent_grid%pelist(1) + (parent_tile-1)*parent_grid%npes_per_tile

       !Call gather grid on the procs that have the required data.
       !Then broadcast from the head PE to the receiving PEs
       if (ANY(parent_grid%pelist == mpp_pe()) .and. parent_tile == parent_grid%tile) then
          call mpp_global_field( &
               parent_grid%domain, &
               var_coarse_this_grid, var_coarse, position=position)
          if (mpp_pe() == sending_proc) then 
             do p=1,recv_npes
                call mpp_send(var_coarse,size(var_coarse),recv_pes(p))
             enddo
          endif
       endif
       if (ANY(recv_pes == mpp_pe())) then
          call mpp_recv(var_coarse, size(var_coarse), sending_proc)
       endif

      

   else

      if (parent_grid%tile == parent_tile) then
         call mpp_global_field( &
              parent_grid%domain, &
              var_coarse_this_grid, var_coarse, position=position)
      end if

      npx = ieg-isg+1+istag
      npy = jeg-jsg+1+jstag

      allocate(tilelist(0:mpp_npes()-1))
      tilelist = 0
      tilelist(mpp_pe()) = parent_grid%tile
      call mpp_sum(tilelist,mpp_npes())
!!$   if (gid == 0) then
!!$      print*, 'tilelist:'
!!$      do n=0,mpp_npes()-1
!!$         print*, n, tilelist(n)
!!$      end do
!!$   end if
!!$   call mpp_sync

      !Don't need to do communication to the other tiles if there is only one tile
      if (maxval(tilelist) == minval(tilelist)) return

      do n=0,mpp_npes()-1
         if (tilelist(n) == parent_tile) then
            if (mpp_pe() == n) then
               !            print*, 'gather_grid: proc ', gid, ' broadcasting '
               !            print*, 'gather_grid: SIZE: ', size(var_coarse), npx, npy
               call mpp_broadcast(var_coarse, npx*npy*npz, n, parent_grid%pelist)
            else 
               !            print*, 'gather_grid: proc ', gid, ' receiving from proc ', n
               call mpp_broadcast(var_coarse, npx*npy*npz, n, parent_grid%pelist)                         
            end if
            exit
         end if
      end do

   endif

 end subroutine gather_grid_3d

 subroutine update_coarse_grid_mpp_2d(var_coarse, var_nest, nest_domain, ind_update, dx, dy, area, &
      isd_p, ied_p, jsd_p, jed_p, is_n, ie_n, js_n, je_n, npx, npy, &
      istag, jstag, r, nestupdate, upoff, nsponge, parent_proc, child_proc, parent_grid)

   integer, intent(IN) :: isd_p, ied_p, jsd_p, jed_p, is_n, ie_n, js_n, je_n
   integer, intent(IN) :: istag, jstag, r, nestupdate, upoff, nsponge
   integer, intent(IN) :: ind_update(isd_p:ied_p+1,jsd_p:jed_p+1,2)
   integer, intent(IN) :: npx, npy
   real, intent(IN)    :: var_nest(is_n:ie_n+istag,js_n:je_n+jstag)
   real, intent(INOUT) :: var_coarse(isd_p:ied_p+istag,jsd_p:jed_p+jstag)
   real, intent(IN)    :: dx(isd:ied,jsd:jed+1)
   real, intent(IN)    :: dy(isd:ied+1,jsd:jed)
   real, intent(IN)    :: area(isd:ied,jsd:jed)
   logical, intent(IN) :: parent_proc, child_proc
   type(fv_atmos_type), intent(INOUT) :: parent_grid
   type(nest_domain_type), intent(INOUT) :: nest_domain

   real :: var_nest_3d(is_n:ie_n+istag,js_n:je_n+jstag,1)
   real :: var_coarse_3d(isd_p:ied_p+istag,jsd_p:jed_p+jstag,1)

   if (child_proc .and. size(var_nest) > 1) var_nest_3d(is_n:ie_n+istag,js_n:je_n+jstag,1) = var_nest(is_n:ie_n+istag,js_n:je_n+jstag)
   if (parent_proc .and. size(var_coarse) > 1) var_coarse_3d(isd_p:ied_p+istag,jsd_p:jed_p,1) = var_coarse(isd_p:ied_p+istag,jsd_p:jed_p+jstag)

   call update_coarse_grid_mpp(var_coarse_3d, var_nest_3d, nest_domain, ind_update, dx, dy, area, &
      isd_p, ied_p, jsd_p, jed_p, is_n, ie_n, js_n, je_n, npx, npy, 1, &
      istag, jstag, r, nestupdate, upoff, nsponge, parent_proc, child_proc, parent_grid)

   if (size(var_coarse) > 1 .and. parent_proc) var_coarse(isd_p:ied_p+istag,jsd_p:jed_p+jstag) = var_coarse_3d(isd_p:ied_p+istag,jsd_p:jed_p,1)

 end subroutine update_coarse_grid_mpp_2d


  subroutine update_coarse_grid_mpp(var_coarse, var_nest, nest_domain, ind_update, dx, dy, area, &
      isd_p, ied_p, jsd_p, jed_p, is_n, ie_n, js_n, je_n, npx, npy, npz, &
      istag, jstag, r, nestupdate, upoff, nsponge, parent_proc, child_proc, parent_grid)

   !This routine assumes the coarse and nested grids are properly
   ! aligned, and that in particular for odd refinement ratios all
   ! coarse-grid points coincide with nested-grid points

   integer, intent(IN) :: isd_p, ied_p, jsd_p, jed_p, is_n, ie_n, js_n, je_n
   integer, intent(IN) :: istag, jstag, npx, npy, npz, r, nestupdate, upoff, nsponge
   integer, intent(IN) :: ind_update(isd_p:ied_p+1,jsd_p:jed_p+1,2)
   real, intent(IN)    :: var_nest(is_n:ie_n+istag,js_n:je_n+jstag,npz)
   real, intent(INOUT) :: var_coarse(isd_p:ied_p+istag,jsd_p:jed_p+jstag,npz)
   real, intent(IN)    :: area(isd:ied,jsd:jed)
   real, intent(IN)    :: dx(isd:ied,jsd:jed+1)
   real, intent(IN)    :: dy(isd:ied+1,jsd:jed)
   logical, intent(IN) :: parent_proc, child_proc
   type(fv_atmos_type), intent(INOUT) :: parent_grid

   type(nest_domain_type), intent(INOUT) :: nest_domain

   integer :: in, jn, ini, jnj, s, qr
   integer :: is_c, ie_c, js_c, je_c, is_f, ie_f, js_f, je_f
   integer :: istart, istop, jstart, jstop, ishift, jshift, j, i, k
   real :: val
   real, allocatable, dimension(:,:,:) :: nest_dat
   real ::  var_nest_send(is_n:ie_n+istag,js_n:je_n+jstag,npz)
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

   call mpp_get_F2C_index(nest_domain, is_c, ie_c, js_c, je_c, is_f, ie_f, js_f, je_f, position=position)
   if (ie_f > is_f .and. je_f > js_f) then
      allocate(nest_dat (is_f:ie_f, js_f:je_f,npz))
   else
      allocate(nest_dat(1,1,1))
   endif
   nest_dat = -600

   if (child_proc) then
!! IF an area average (for istag == jstag == 0) or a linear average then multiply in the areas before sending data
   if (istag == 0 .and. jstag == 0) then
      select case (nestupdate)
      case (0,3,4,5)
         
         do k=1,npz
         do j=js_n,je_n
         do i=is_n,ie_n

            var_nest_send(i,j,k) = var_nest(i,j,k)

         end do
         end do
         end do

      case (1,2,6,7,8)
         
         do k=1,npz
         do j=js_n,je_n
         do i=is_n,ie_n

            var_nest_send(i,j,k) = var_nest(i,j,k)*area(i,j)

         end do
         end do
         end do

      end select
   else if (istag == 0 .and. jstag > 0) then

      select case (nestupdate) 
      case (0,2,3,4,5)

         do k=1,npz
         do j=js_n,je_n+1
         do i=is_n,ie_n

            var_nest_send(i,j,k) = var_nest(i,j,k)
            
         end do
         end do
         end do

      case (1,6,7,8)

         do k=1,npz
         do j=js_n,je_n+1
         do i=is_n,ie_n

            var_nest_send(i,j,k) = var_nest(i,j,k)*dx(i,j)
            
         end do
         end do
         end do

      case default

         call mpp_error(FATAL, 'nestupdate type not implemented')

      end select

   else if (istag > 0 .and. jstag == 0) then
      select case (nestupdate) 
      case (0,2,3,4,5)

         do k=1,npz
         do j=js_n,je_n
         do i=is_n,ie_n+1

            var_nest_send(i,j,k) = var_nest(i,j,k)

         end do
         end do
         end do

      case (1,6,7,8)   !averaging update; in-line average for face-averaged values instead of areal average

         do k=1,npz
         do j=js_n,je_n
         do i=is_n,ie_n+1

            var_nest_send(i,j,k) = var_nest(i,j,k)*dy(i,j)

         end do
         end do
         end do

      case default

         call mpp_error(FATAL, 'nestupdate type not implemented')

      end select

   else
      
      call mpp_error(FATAL, "Cannot have both nonzero istag and jstag.")

   endif
   endif

      call timing_on('COMM_TOTAL')
   call mpp_update_nest_coarse(var_nest_send, nest_domain, nest_dat, position=position)
      call timing_off('COMM_TOTAL')

   s = r/2 !rounds down (since r > 0)
   qr = r*upoff + nsponge - s

   if (parent_proc) then
   if (istag == 0 .and. jstag == 0) then

      select case (nestupdate) 
      case (0,3,4,5) !Interpolation update (0 = all variables, 4 = not on delp or q, 3 = only on winds, 5 = "remap-update")

         do k=1,npz
         do j=jsd_p,jed_p
         do i=isd_p,ied_p

            in = ind_update(i,j,1)
            jn = ind_update(i,j,2)

            if (in < max(1+qr,is_f) .or. in > min(npx-1-qr-s,ie_f) .or. &
                 jn < max(1+qr,js_f) .or. jn > min(npy-1-qr-s,je_f)) cycle

            if (mod(r,2) == 1) then
               val = nest_dat(in+s,jn+s,k)
            else
               !SUBTRACTING 1 because the first index LARGER THAN the median half-index is jn+s
               val = 0.25*( &
                    nest_dat(in+s,jn+s,k)   + nest_dat(in+s-1,jn+s,k) + &
                    nest_dat(in+s,jn+s-1,k) + nest_dat(in+s-1,jn+s-1,k) )
            end if
            var_coarse(i,j,k) = val

         end do
         end do
         end do

      case (1,2,6,7,8) ! 1 = Conserving update on all variables; 2 = conserving update for cell-centered values; 6 = conserving remap-update

         do k=1,npz
         do j=jsd_p,jed_p
         do i=isd_p,ied_p

            in = ind_update(i,j,1)
            jn = ind_update(i,j,2)

            if (in < max(1+qr,is_f) .or. in > min(npx-1-qr-r+1,ie_f) .or. &
                 jn < max(1+qr,js_f) .or. jn > min(npy-1-qr-r+1,je_f)) cycle

            val = 0.
            do jnj=jn,jn+r-1
               do ini=in,in+r-1
                  val = val + nest_dat(ini,jnj,k)
               end do
            end do            

            !var_coarse(i,j,k) = val/r**2.

            !!! CLEANUP: Couldn't rarea and rdx and rdy be built into the weight arrays?
            !!!    Two-way updates do not yet have weights, tho
            var_coarse(i,j,k) = val*parent_grid%gridstruct%rarea(i,j)

         end do
         end do
         end do


      case default

         call mpp_error(FATAL, 'nestupdate type not implemented')


      end select

   else if (istag == 0 .and. jstag > 0) then


      select case (nestupdate) 
      case (0,2,3,4,5)

         do k=1,npz
         do j=jsd_p,jed_p+1
         do i=isd_p,ied_p

            in = ind_update(i,j,1)
            jn = ind_update(i,j,2)

            if (in < max(1+qr,is_f) .or. in > min(npx-1-qr-s,ie_f) .or. &
                 jn < max(1+qr,js_f) .or. jn > min(npy-1-qr+1-s,je_f)) cycle

            if (mod(r,2) == 1) then
               val = nest_dat(in+s,jn,k)
            else
               val = 0.5*(nest_dat(in+s,jn,k) + nest_dat(in+s-1,jn,k) )
            end if
            var_coarse(i,j,k) = val
         end do
         end do
         end do

      case (1,6,7,8)

         do k=1,npz
         do j=jsd_p,jed_p+1
         do i=isd_p,ied_p

            in = ind_update(i,j,1)
            jn = ind_update(i,j,2)

            if (in < max(1+qr,is_f) .or. in > min(npx-1-qr-r+1,ie_f) .or. &
                 jn < max(1+qr+s,js_f) .or. jn > min(npy-1-qr-s+1,je_f)) cycle

            val = 0.
               do ini=in,in+r-1
                  val = val + nest_dat(ini,jn,k)
               end do

!            var_coarse(i,j,k) = val/r
            var_coarse(i,j,k) = val*parent_grid%gridstruct%rdx(i,j)

         end do
         end do
         end do

      case default

         call mpp_error(FATAL, 'nestupdate type not implemented')

      end select

   else if (istag > 0 .and. jstag == 0) then

      select case (nestupdate) 
      case (0,2,3,4,5)

         do k=1,npz
         do j=jsd_p,jed_p
         do i=isd_p,ied_p+1

            in = ind_update(i,j,1)
            jn = ind_update(i,j,2)

            if (in < max(1+qr,is_f) .or. in > min(npx-1-qr-s+1,ie_f) .or. &
                 jn < max(1+qr,js_f) .or. jn > min(npy-1-qr-s,je_f)) cycle

            if (mod(r,2) == 1) then
               val = nest_dat(in,jn+s,k)
            else
               val = 0.5*(nest_dat(in,jn+s,k)+nest_dat(in,jn+s-1,k))
            end if
            var_coarse(i,j,k) = val

         end do
         end do
         end do

      case (1,6,7,8)   !averaging update; in-line average for face-averaged values instead of areal average

         do k=1,npz
         do j=jsd_p,jed_p
         do i=isd_p,ied_p+1

            in = ind_update(i,j,1)
            jn = ind_update(i,j,2)

            if (in < max(1+qr+s,is_f) .or. in > min(npx-1-qr-s+1,ie_f) .or. &
                 jn < max(1+qr,js_f) .or. jn > min(npy-1-qr-r+1,je_f)) cycle

            val = 0.
            do jnj=jn,jn+r-1
                  val = val + nest_dat(in,jnj,k)
            end do

!            var_coarse(i,j,k) = val/r
            var_coarse(i,j,k) = val*parent_grid%gridstruct%rdy(i,j)

         end do
         end do
         end do

      case default

         call mpp_error(FATAL, 'nestupdate type not implemented')

      end select

   end if


   endif
   deallocate(nest_dat)
   
 end subroutine update_coarse_grid_mpp


   
end module boundary_mod
