!***********************************************************************
!>@brief!   Provides subroutines to debug and log moving nest functionality 
!!>@author Bill Ramstrom, AOML/HRD   01/15/2021
!
!   This code is in a separate module so that code review and optimization 
!     can focus on the algorithm code in fv_moving_nest.F90 that implements
!     the core functionality.  These routines will likely be disabled or 
!     removed before operational implementation.
! 
! =======================================================================!
!


! Notes


module fv_moving_nest_logging_mod

#ifdef MOVING_NEST

  use mpp_mod,                  only: mpp_pe
  use fv_arrays_mod,            only: R_GRID
  use fv_moving_nest_utils_mod, only: grid_geometry, fill_grid_from_supergrid
  use fv_arrays_mod,            only: fv_grid_type, fv_nest_type, fv_atmos_type


  interface check_array
     module procedure check_2d_array
     module procedure check_2d64_array
     module procedure check_3d_array
     module procedure check_4d_array
     module procedure check_4d64_array
  end interface check_array

  interface check_local_array
     module procedure check_2d_local_array
     module procedure check_3d_local_array
  end interface check_local_array


contains


  !==================================================================================================
  !
  !  Array Checking Section
  !
  !==================================================================================================


  subroutine check_2d_array(array, this_pe, var_name, min_range, max_range)
    real, intent(in), allocatable          :: array(:,:)
    integer, intent(in)                    :: this_pe
    character(len=*), intent(in)           :: var_name
    real, intent(in)                       :: min_range, max_range

    integer  :: i,j
    integer  :: num_invalid
    integer  :: num_valid
    real     :: eps = 0.0001
    real     :: invalid_last

    invalid_last = 0.0

    if (allocated(array)) then

       print '("[INFO] WDR 2Darray allocated  npe=",I0," ",A32,"(",I4,":",I4,",",I4,":",I4,")")', this_pe, var_name, lbound(array,1), ubound(array,1), lbound(array,2), ubound(array,2)

       num_invalid = 0
       num_valid = 0

       do i = lbound(array,1), ubound(array,1)
          do j =lbound(array,2), ubound(array,2)
             if (array(i,j) < min_range - eps) then
                num_invalid = num_invalid + 1
                invalid_last = array(i,j)
             elseif (array(i,j) > max_range + eps) then
                num_invalid = num_invalid + 1
                invalid_last = array(i,j)
             else 
                num_valid = num_valid + 1
             end if
          end do
       end do

       if (num_invalid > 0 ) then
          print '("[ERROR] WDR 2Darray invalid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0," last invalid=",E12.5)', this_pe, var_name, num_invalid, num_valid, invalid_last
       else
          print '("[INFO] WDR 2Darray all valid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
       end if

    else
       print '("[INFO] WDR 2Darray not allocated  npe=",I0," ",A32)', this_pe, var_name
    end if


  end subroutine check_2d_array

  subroutine check_2d64_array(array, this_pe, var_name, min_range, max_range)
    real(kind=R_GRID), intent(in), allocatable  :: array(:,:)
    integer, intent(in)                    :: this_pe
    character(len=*), intent(in)           :: var_name
    real(kind=R_GRID), intent(in)                       :: min_range, max_range

    integer  :: i,j
    integer  :: num_invalid
    integer  :: num_valid
    real     :: eps = 0.0001
    real     :: invalid_last

    invalid_last = 0.0

    if (allocated(array)) then

       print '("[INFO] WDR 2D64array allocated  npe=",I0," ",A32,"(",I4,":",I4,",",I4,":",I4,")")', this_pe, var_name, lbound(array,1), ubound(array,1), lbound(array,2), ubound(array,2)

       num_invalid = 0
       num_valid = 0

       do i = lbound(array,1), ubound(array,1)
          do j =lbound(array,2), ubound(array,2)
             if (array(i,j) < min_range - eps) then
                num_invalid = num_invalid + 1
                invalid_last = array(i,j)
             elseif (array(i,j) > max_range + eps) then
                num_invalid = num_invalid + 1
                invalid_last = array(i,j)
             else 
                num_valid = num_valid + 1
             end if
          end do
       end do

       if (num_invalid > 0 ) then
          print '("[ERROR] WDR 2D64array invalid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0," last invalid=",E12.5)', this_pe, var_name, num_invalid, num_valid, invalid_last
       else
          print '("[INFO] WDR 2D64array all valid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
       end if

    else
       print '("[INFO] WDR 2D64array not allocated  npe=",I0," ",A32)', this_pe, var_name
    end if


  end subroutine check_2d64_array


  subroutine check_2d_local_array(array, this_pe, var_name, min_range, max_range)
    real, intent(in)  :: array(:,:)
    integer, intent(in)                    :: this_pe
    character(len=*), intent(in)           :: var_name
    real, intent(in)                       :: min_range, max_range

    integer  :: i,j
    integer  :: num_invalid
    integer  :: num_valid
    real     :: eps = 0.0001

    print '("[INFO] WDR 2DLarray allocated  npe=",I0," ",A32,"(",I4,":",I4,",",I4,":",I4,")")', this_pe, var_name, lbound(array,1), ubound(array,1), lbound(array,2), ubound(array,2)

    num_invalid = 0
    num_valid = 0

    do i = lbound(array,1), ubound(array,1)
       do j =lbound(array,2), ubound(array,2)
          if (array(i,j) < min_range - eps) then
             num_invalid = num_invalid + 1
          elseif (array(i,j) > max_range + eps) then
             num_invalid = num_invalid + 1
          else 
             num_valid = num_valid + 1
          end if
       end do
    end do

    if (num_invalid > 0 ) then
       print '("[ERROR] WDR 2DLarray invalid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
    else
       print '("[INFO] WDR 2DLarray all valid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
    end if

  end subroutine check_2d_local_array


  subroutine check_3d_array(array, this_pe, var_name, min_range, max_range)
    real, intent(in), allocatable  :: array(:,:,:)
    integer, intent(in)                    :: this_pe
    character(len=*), intent(in)           :: var_name
    real, intent(in)                       :: min_range, max_range

    integer  :: i,j,k
    integer  :: num_invalid
    integer  :: num_valid
    real     :: eps = 0.0001

    if (allocated(array)) then

       print '("[INFO] WDR 3Darray allocated  npe=",I0," ",A32,"(",I4,":",I4,",",I4,":",I4,",",I4,":",I4,")")', this_pe, var_name, lbound(array,1), ubound(array,1), lbound(array,2), ubound(array,2), lbound(array,3), ubound(array,3)

       num_invalid = 0
       num_valid = 0

       do i = lbound(array,1), ubound(array,1)
          do j =lbound(array,2), ubound(array,2)
             do k =lbound(array,3), ubound(array,3)
                if (isnan(array(i,j,k))) then 
                   num_invalid = num_invalid + 1
                elseif (array(i,j,k) < min_range - eps) then
                   num_invalid = num_invalid + 1
                elseif (array(i,j,k) > max_range + eps) then
                   num_invalid = num_invalid + 1
                else 
                   num_valid = num_valid + 1
                end if
             end do
          end do
       end do

       if (num_invalid > 0 ) then
          print '("[ERROR] WDR 3Darray invalid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
       else
          print '("[INFO] WDR 3Darray all valid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
       end if

    else
       print '("[INFO] WDR 3Darray not allocated  npe=",I0," ",A32)', this_pe, var_name
    end if


  end subroutine check_3d_array

  subroutine check_3d_local_array(array, this_pe, var_name, min_range, max_range)
    real, intent(in) :: array(:,:,:)
    integer, intent(in)                    :: this_pe
    character(len=*), intent(in)           :: var_name
    real, intent(in)                       :: min_range, max_range

    integer  :: i,j,k
    integer  :: num_invalid
    integer  :: num_valid
    real     :: eps = 0.0001


    print '("[INFO] WDR 3DLarray bounds npe=",I0," ",A32,"(",I4,":",I4,",",I4,":",I4,",",I4,":",I4,")")', this_pe, var_name, lbound(array,1), ubound(array,1), lbound(array,2), ubound(array,2), lbound(array,3), ubound(array,3)

    num_invalid = 0
    num_valid = 0

    do i = lbound(array,1), ubound(array,1)
       do j =lbound(array,2), ubound(array,2)
          do k =lbound(array,3), ubound(array,3)
             if (isnan(array(i,j,k))) then 
                num_invalid = num_invalid + 1
             elseif (array(i,j,k) < min_range - eps) then
                num_invalid = num_invalid + 1
             elseif (array(i,j,k) > max_range + eps) then
                num_invalid = num_invalid + 1
             else 
                num_valid = num_valid + 1
             end if
          end do
       end do
    end do

    if (num_invalid > 0 ) then
       print '("[ERROR] WDR 3DLarray invalid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
    else
       print '("[INFO] WDR 3DLarray all valid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
    end if

  end subroutine check_3d_local_array


  subroutine check_4d_array(array, this_pe, var_name, min_range, max_range)
    real, intent(in), allocatable                      :: array(:,:,:,:)
    integer, intent(in)                    :: this_pe
    character(len=*), intent(in)           :: var_name
    real, intent(in)                       :: min_range, max_range

    integer  :: i,j,k,v
    integer  :: num_invalid
    integer  :: num_valid
    real     :: eps = 0.0001

    if (allocated(array)) then

       print '("[INFO] WDR 4Darray allocated  npe=",I0," ",A32,"(",I4,":",I4,",",I4,":",I4,",",I4,":",I4,",",I4,":",I4,")")', this_pe, var_name, lbound(array,1), ubound(array,1), lbound(array,2), ubound(array,2), lbound(array,3), ubound(array,3), lbound(array,4), ubound(array,4)

       num_invalid = 0
       num_valid = 0

       do i = lbound(array,1), ubound(array,1)
          do j =lbound(array,2), ubound(array,2)
             do k =lbound(array,3), ubound(array,3)
                do v =lbound(array,4), ubound(array,4)
                   if (isnan(array(i,j,k,v))) then 
                      num_invalid = num_invalid + 1
                   elseif (array(i,j,k,v) < min_range - eps) then
                      num_invalid = num_invalid + 1
                   elseif (array(i,j,k,v) > max_range + eps) then
                      num_invalid = num_invalid + 1
                   else 
                      num_valid = num_valid + 1
                   end if
                end do
             end do
          end do
       end do

       if (num_invalid > 0 ) then
          print '("[ERROR] WDR 4Darray invalid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
       else
          print '("[INFO] WDR 4Darray all valid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
       end if
    else
       print '("[INFO] WDR 4Darray not allocated  npe=",I0," ",A32)', this_pe, var_name
    end if
  end subroutine check_4d_array

  subroutine check_4d64_array(array, this_pe, var_name, min_range, max_range)
    real(kind=R_GRID), intent(in), allocatable                      :: array(:,:,:,:)
    integer, intent(in)                    :: this_pe
    character(len=*), intent(in)           :: var_name
    real, intent(in)                       :: min_range, max_range

    integer  :: i,j,k,v
    integer  :: num_invalid
    integer  :: num_valid
    real     :: eps = 0.0001

    if (allocated(array)) then

       print '("[INFO] WDR 4Darray allocated  npe=",I0," ",A32,"(",I4,":",I4,",",I4,":",I4,",",I4,":",I4,",",I4,":",I4,")")', this_pe, var_name, lbound(array,1), ubound(array,1), lbound(array,2), ubound(array,2), lbound(array,3), ubound(array,3), lbound(array,4), ubound(array,4)

       num_invalid = 0
       num_valid = 0

       do i = lbound(array,1), ubound(array,1)
          do j =lbound(array,2), ubound(array,2)
             do k =lbound(array,3), ubound(array,3)
                do v =lbound(array,4), ubound(array,4)
                   if (isnan(array(i,j,k,v))) then 
                      num_invalid = num_invalid + 1
                   elseif (array(i,j,k,v) < min_range - eps) then
                      num_invalid = num_invalid + 1
                   elseif (array(i,j,k,v) > max_range + eps) then
                      num_invalid = num_invalid + 1
                   else 
                      num_valid = num_valid + 1
                   end if
                end do
             end do
          end do
       end do

       if (num_invalid > 0 ) then
          print '("[ERROR] WDR 4Darray invalid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
       else
          print '("[INFO] WDR 4Darray all valid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
       end if
    else
       print '("[INFO] WDR 4Darray not allocated  npe=",I0," ",A32)', this_pe, var_name
    end if
  end subroutine check_4d64_array





  !==================================================================================================
  !
  !  Debugging Function Section
  !
  !==================================================================================================



  subroutine grid_equal(grid1, grid2, tag, this_pe, is_equal)
    real, allocatable, intent(in)    :: grid1(:,:,:)
    real, allocatable, intent(in)    :: grid2(:,:,:)
    character(len=*), intent(in)     :: tag
    integer, intent(in)              :: this_pe
    logical, intent(out)             :: is_equal

    integer :: x,y,z

    real                 :: pi = 4 * atan(1.0d0)
    real                 :: rad2deg

    rad2deg = 180.0 / pi

    is_equal = .true.

    do x=1,3
       if (lbound(grid1,x) /= lbound(grid2,x)) then
          print '("[ERROR] WDR grid_equal ",A16," npe=",I0," lbound mismatch ",I0, I0,I0)', tag, x, lbound(grid1,x), lbound(grid2,x)
          is_equal = .false.
       end if
       if (ubound(grid1,x) /= ubound(grid2,x)) then
          print '("[ERROR] WDR grid_equal ",A16," npe=",I0," ubound mismatch ",I0, I0,I0)', tag, x, ubound(grid1,x), ubound(grid2,x)
          is_equal = .false.
       end if
    end do

    if (is_equal) then
       do x=lbound(grid1,1), ubound(grid1,1)
          do y=lbound(grid1,2), ubound(grid1,2)
             do z=lbound(grid1,3), ubound(grid1,3)
                if ( abs(grid1(x,y,z) - grid2(x,y,z)) > 0.0001 ) then
                   print '("[ERROR] WDR grid_equal ",A16," npe=",I0," DEG value mismatch at (",I0,",",I0,",",I0,") ",F15.11, " ",F15.11, " ",F15.11)', tag, this_pe, x, y, z, grid1(x,y,z)*rad2deg, grid2(x,y,z)*rad2deg, grid1(x,y,z)*rad2deg - grid2(x,y,z)*rad2deg

                   print '("[ERROR] WDR grid_equal ",A16," npe=",I0," RAD value mismatch at (",I0,",",I0,",",I0,") ",F15.11, " ",F15.11, " ",F15.11)', tag, this_pe, x, y, z, grid1(x,y,z), grid2(x,y,z), grid1(x,y,z) - grid2(x,y,z)
                   is_equal = .false.
                else
                   print '("[INFO]  WDR grid_equal ",A16," npe=",I0," DEG value match    at (",I0,",",I0,",",I0,") ",F15.11, " ",F15.11, " ",F15.11)', tag, this_pe, x, y, z, grid1(x,y,z)*rad2deg, grid2(x,y,z)*rad2deg, grid1(x,y,z)*rad2deg - grid2(x,y,z)*rad2deg

                   print '("[INFO]  WDR grid_equal ",A16," npe=",I0," RAD value match    at (",I0,",",I0,",",I0,") ",F15.11, " ",F15.11, " ",F15.11)', tag, this_pe, x, y, z, grid1(x,y,z), grid2(x,y,z), grid1(x,y,z) - grid2(x,y,z)
                end if
             end do
          end do
       end do
    end if

    if (is_equal) then
       print '("[INFO] WDR grid_equal ",A16," npe=",I0," MATCH.")', tag, this_pe
    else
       print '("[ERROR] WDR grid_equal ",A16," npe=",I0," MISMATCH.")', tag, this_pe
    end if

  end subroutine grid_equal


  subroutine show_atm_grids(Atm, n)
    type(fv_atmos_type), allocatable, intent(in) :: Atm(:)
    integer, intent(in)                          :: n

    integer :: x,y
    real(kind=R_GRID) :: pi = 4 * atan(1.0d0)
    real                :: pi180
    real                :: rad2deg, deg2rad

    pi180 = pi / 180.0
    deg2rad = pi / 180.0
    rad2deg = 1.0 / pi180


    print *, "[INFO] WDR MV_NST2 shape(Atm(1)%grid_global)=", shape(Atm(1)%grid_global)
    print '("[INFO] WDR MV_NST2 bounds1 (Atm(1)%grid_global)=",I0,"-",I0)', lbound(Atm(1)%grid_global,1), ubound(Atm(1)%grid_global,1)
    print '("[INFO] WDR MV_NST2 bounds2 (Atm(1)%grid_global)=",I0,"-",I0)', lbound(Atm(1)%grid_global,2), ubound(Atm(1)%grid_global,2)
    print '("[INFO] WDR MV_NST2 bounds3 (Atm(1)%grid_global)=",I0,"-",I0)', lbound(Atm(1)%grid_global,3), ubound(Atm(1)%grid_global,3)
    print '("[INFO] WDR MV_NST2 bounds4 (Atm(1)%grid_global)=",I0,"-",I0)', lbound(Atm(1)%grid_global,4), ubound(Atm(1)%grid_global,4)


    print *, "[INFO] WDR MV_NST2 shape(Atm(n)%grid_global)=", shape(Atm(n)%grid_global)
    print '("[INFO] WDR MV_NST2 bounds1 (Atm(n)%grid_global)=",I0,"-",I0)', lbound(Atm(n)%grid_global,1), ubound(Atm(n)%grid_global,1)
    print '("[INFO] WDR MV_NST2 bounds2 (Atm(n)%grid_global)=",I0,"-",I0)', lbound(Atm(n)%grid_global,2), ubound(Atm(n)%grid_global,2)
    print '("[INFO] WDR MV_NST2 bounds3 (Atm(n)%grid_global)=",I0,"-",I0)', lbound(Atm(n)%grid_global,3), ubound(Atm(n)%grid_global,3)
    print '("[INFO] WDR MV_NST2 bounds4 (Atm(n)%grid_global)=",I0,"-",I0)', lbound(Atm(n)%grid_global,4), ubound(Atm(n)%grid_global,4)



    print *, "[INFO] WDR MV_NST2 shape(Atm(n)%gridstruct%grid)=", shape(Atm(n)%gridstruct%grid)
    print '("[INFO] WDR MV_NST2 bounds1 (Atm(n)%gridstruct%grid)=",I0,"-",I0)', lbound(Atm(n)%gridstruct%grid,1), ubound(Atm(n)%gridstruct%grid,1)
    print '("[INFO] WDR MV_NST2 bounds2 (Atm(n)%gridstruct%grid)=",I0,"-",I0)', lbound(Atm(n)%gridstruct%grid,2), ubound(Atm(n)%gridstruct%grid,2)
    print '("[INFO] WDR MV_NST2 bounds3 (Atm(n)%gridstruct%grid)=",I0,"-",I0)', lbound(Atm(n)%gridstruct%grid,3), ubound(Atm(n)%gridstruct%grid,3)


    print *, "[INFO] WDR MV_NST2 shape(Atm(n)%gridstruct%agrid)=", shape(Atm(n)%gridstruct%agrid)
    print '("[INFO] WDR MV_NST2 bounds1 (Atm(n)%gridstruct%agrid)=",I0,"-",I0)', lbound(Atm(n)%gridstruct%agrid,1), ubound(Atm(n)%gridstruct%agrid,1)
    print '("[INFO] WDR MV_NST2 bounds2 (Atm(n)%gridstruct%agrid)=",I0,"-",I0)', lbound(Atm(n)%gridstruct%agrid,2), ubound(Atm(n)%gridstruct%agrid,2)
    print '("[INFO] WDR MV_NST2 bounds3 (Atm(n)%gridstruct%agrid)=",I0,"-",I0)', lbound(Atm(n)%gridstruct%agrid,3), ubound(Atm(n)%gridstruct%agrid,3)


    x = lbound(Atm(n)%gridstruct%agrid,1)
    y = lbound(Atm(n)%gridstruct%agrid,2)
    print '("[INFO] WDR GRD_SHOa atmosphere.F90 Atm(n)%agrid(",I0,",",I0,")=",F10.5, F10.5)', x, y, Atm(n)%gridstruct%agrid(x,y,2)*rad2deg, Atm(n)%gridstruct%agrid(x,y,1)*rad2deg 

    x = ubound(Atm(n)%gridstruct%agrid,1)
    y = ubound(Atm(n)%gridstruct%agrid,2)
    print '("[INFO] WDR GRD_SHOb atmosphere.F90 Atm(n)%agrid(",I0,",",I0,")=",F10.5, F10.5)', x, y, Atm(n)%gridstruct%agrid(x,y,2)*rad2deg, Atm(n)%gridstruct%agrid(x,y,1)*rad2deg 


    print *, "[INFO] WDR MV_NST2 shape(Atm(n)%gridstruct%grid_64)=", shape(Atm(n)%gridstruct%grid_64)
    print '("[INFO] WDR MV_NST2 bounds1 (Atm(n)%gridstruct%grid_64)=",I0,"-",I0)', lbound(Atm(n)%gridstruct%grid_64,1), ubound(Atm(n)%gridstruct%grid_64,1)
    print '("[INFO] WDR MV_NST2 bounds2 (Atm(n)%gridstruct%grid_64)=",I0,"-",I0)', lbound(Atm(n)%gridstruct%grid_64,2), ubound(Atm(n)%gridstruct%grid_64,2)
    print '("[INFO] WDR MV_NST2 bounds3 (Atm(n)%gridstruct%grid_64)=",I0,"-",I0)', lbound(Atm(n)%gridstruct%grid_64,3), ubound(Atm(n)%gridstruct%grid_64,3)

    print *, "[INFO] WDR MV_NST2 shape(Atm(n)%gridstruct%agrid_64)=", shape(Atm(n)%gridstruct%agrid_64)
    print '("[INFO] WDR MV_NST2 bounds1 (Atm(n)%gridstruct%agrid_64)=",I0,"-",I0)', lbound(Atm(n)%gridstruct%agrid_64,1), ubound(Atm(n)%gridstruct%agrid_64,1)
    print '("[INFO] WDR MV_NST2 bounds2 (Atm(n)%gridstruct%agrid_64)=",I0,"-",I0)', lbound(Atm(n)%gridstruct%agrid_64,2), ubound(Atm(n)%gridstruct%agrid_64,2)
    print '("[INFO] WDR MV_NST2 bounds3 (Atm(n)%gridstruct%agrid_64)=",I0,"-",I0)', lbound(Atm(n)%gridstruct%agrid_64,3), ubound(Atm(n)%gridstruct%agrid_64,3)

  end subroutine show_atm_grids


  subroutine show_tile_geo(tile_geo, this_pe, var_name)
    type(grid_geometry)                    :: tile_geo
    integer, intent(in)                    :: this_pe
    character(len=*), intent(in)           :: var_name


    print '("[INFO] WDR 2Darray npe=",I0," ",A32, "nx=", I0," ny=", I0," nxp=",I0," nyp=",I0)', this_pe, var_name, tile_geo%nx, tile_geo%ny, tile_geo%nxp, tile_geo%nyp


    call check_2d_array(tile_geo%lats, this_pe, var_name // "%lats", -90.0, 90.0)
    call check_2d_array(tile_geo%lons, this_pe, var_name // "%lons", -360.0, 360.0)
    call check_2d_array(tile_geo%dx, this_pe, var_name // "%dx", 0.0, 1.0e9)
    call check_2d_array(tile_geo%dy, this_pe, var_name // "%dy", 0.0, 1.0e9)
    call check_2d_array(tile_geo%area, this_pe, var_name // "%area", 0.0, 1.0e9)


  end subroutine show_tile_geo


  subroutine show_atm_array4(tag, array, array_name, atm_n, this_pe)
    character(len=*), intent(in)                                    :: tag
    real(kind=R_GRID), allocatable, dimension(:,:,:,:), intent(in)  :: array
    character(len=*), intent(in)                                    :: array_name
    integer, intent(in)                                             :: atm_n, this_pe

    if (allocated(array)) then
       print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%",A12,"(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', tag, this_pe, atm_n, trim(array_name), lbound(array, 1), ubound(array, 1), lbound(array, 2), ubound(array, 2), lbound(array, 3), ubound(array, 3), lbound(array, 4), ubound(array, 4)

    else
       print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%",A12," is not allocated.")', tag, this_pe, trim(array_name)
    end if

  end subroutine show_atm_array4

  subroutine show_atm_neststruct(tag, neststruct, atm_n, this_pe)
    character(len=*), intent(in)                                    :: tag
    type(fv_nest_type), intent(in) :: neststruct
    integer, intent(in)  :: atm_n, this_pe

    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%parent_tile=",I0," %refinement=",I0)', tag, this_pe, atm_n, neststruct%parent_tile, neststruct%refinement
    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%nested=",L1," %ioffset=",I0," %joffset=",I0)', tag, this_pe, atm_n, neststruct%nested, neststruct%ioffset, neststruct%joffset

    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%nested=",L1," %isu=",I0," %ieu=",I0," %jsu=",I0," %jeu=",I0)', tag, this_pe, atm_n, neststruct%nested, neststruct%isu,  neststruct%ieu,  neststruct%jsu,  neststruct%jeu

    ! WDR ind_update_h seems to have been removed in recent version of the dycore
    !    if (allocated(neststruct%ind_update_h)) then
    !       print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%ind_update_h(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', tag, this_pe, atm_n, &
    !            lbound(neststruct%ind_update_h,1),  ubound(neststruct%ind_update_h,1), &
    !            lbound(neststruct%ind_update_h,2),  ubound(neststruct%ind_update_h,2), &
    !            lbound(neststruct%ind_update_h,3),  ubound(neststruct%ind_update_h,3)
    !
    !       if (ubound(neststruct%ind_update_h,1) > lbound(neststruct%ind_update_h,1)) then
    !          print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") val neststruct%ind_update_h(",I0,",",I0,",",I0,")=",I0)', tag, this_pe, atm_n, &
    !               lbound(neststruct%ind_update_h,1),  lbound(neststruct%ind_update_h,2),  lbound(neststruct%ind_update_h,3), & 
    !               neststruct%ind_update_h(lbound(neststruct%ind_update_h,1),  lbound(neststruct%ind_update_h,2),  lbound(neststruct%ind_update_h,3))
    !          print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") val neststruct%ind_update_h(",I0,",",I0,",",I0,")=",I0)', tag, this_pe, atm_n, &
    !               lbound(neststruct%ind_update_h,1),  lbound(neststruct%ind_update_h,2),  ubound(neststruct%ind_update_h,3), & 
    !               neststruct%ind_update_h(lbound(neststruct%ind_update_h,1),  lbound(neststruct%ind_update_h,2),  ubound(neststruct%ind_update_h,3))
    !
    !          print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") val neststruct%ind_update_h(",I0,",",I0,",",I0,")=",I0)', tag, this_pe, atm_n, &
    !               lbound(neststruct%ind_update_h,1)+4,  lbound(neststruct%ind_update_h,2)+4,  lbound(neststruct%ind_update_h,3), & 
    !               neststruct%ind_update_h(lbound(neststruct%ind_update_h,1)+4,  lbound(neststruct%ind_update_h,2)+4,  lbound(neststruct%ind_update_h,3))
    !          print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") val neststruct%ind_update_h(",I0,",",I0,",",I0,")=",I0)', tag, this_pe, atm_n, &
    !               lbound(neststruct%ind_update_h,1)+4,  lbound(neststruct%ind_update_h,2)+4,  ubound(neststruct%ind_update_h,3), & 
    !               neststruct%ind_update_h(lbound(neststruct%ind_update_h,1)+4,  lbound(neststruct%ind_update_h,2)+4,  ubound(neststruct%ind_update_h,3))
    !
    !
    !       end if
    !    else
    !       print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%ind_update_h is not allocated.")', tag, this_pe, atm_n
    !    end if



    ! WDR nest_domain_all appears to be obsolete in new dycore
    !if (allocated(neststruct%nest_domain_all)) then
    !   print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%nest_domain_all(",I0,"-",I0,")")', tag, this_pe, atm_n, lbound(neststruct%nest_domain_all), ubound(neststruct%nest_domain_all)
    !else
    !   print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%nest_domain_all is not allocated.")', tag, this_pe, atm_n
    !end if

    ! WDR nest_domain has moved to fv_mp_mod.F90 as a global
    !print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%nest_domain%tile_fine=",I0," %tile_coarse=",I0)', tag, this_pe, atm_n, neststruct%nest_domain%tile_fine, neststruct%nest_domain%tile_coarse

    !print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%nest_domain%istart_fine=",I0," %iend_fine=",I0)', tag, this_pe, atm_n, neststruct%nest_domain%istart_fine,  neststruct%nest_domain%iend_fine
    !print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%nest_domain%jstart_fine=",I0," %jend_fine=",I0)', tag, this_pe, atm_n, neststruct%nest_domain%jstart_fine,  neststruct%nest_domain%jend_fine

    !print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%nest_domain%istart_coarse=",I0," %iend_coarse=",I0)', tag, this_pe, atm_n, neststruct%nest_domain%istart_coarse,  neststruct%nest_domain%iend_coarse
    !print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%nest_domain%jstart_coarse=",I0," %jend_coarse=",I0)', tag, this_pe, atm_n, neststruct%nest_domain%jstart_coarse,  neststruct%nest_domain%jend_coarse


  end subroutine show_atm_neststruct


  subroutine show_atm_gridstruct(tag, gridstruct, atm_n, this_pe)
    character(len=*), intent(in)   :: tag
    type(fv_grid_type), intent(in) :: gridstruct
    integer, intent(in)            :: atm_n, this_pe

    ! nested is a pointer.  
    if (associated(gridstruct%nested)) then
       print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") gridstruct%nested=",L1)', tag, this_pe, atm_n, gridstruct%nested
    else
       print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") gridstruct%nested is not set.")', tag, this_pe, atm_n
    end if

    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") gridstruct%cubed_sphere=",L1)', tag, this_pe, atm_n, gridstruct%cubed_sphere
    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") gridstruct%have_north_pole=",L1)', tag, this_pe, atm_n, gridstruct%have_north_pole
    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") gridstruct%have_south_pole=",L1)', tag, this_pe, atm_n, gridstruct%have_south_pole

    if (allocated(gridstruct%agrid)) then
       print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") gridstruct%agrid(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', tag, this_pe, atm_n, lbound(gridstruct%agrid, 1), ubound(gridstruct%agrid, 1), lbound(gridstruct%agrid, 2), ubound(gridstruct%agrid, 2), lbound(gridstruct%agrid, 3), ubound(gridstruct%agrid, 3)
    else
       print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") gridstruct%agrid is not allocated.")', tag, this_pe, atm_n
    end if

    if (allocated(gridstruct%grid)) then
       print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") gridstruct%grid(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', tag, this_pe, atm_n, lbound(gridstruct%grid, 1), ubound(gridstruct%grid, 1), lbound(gridstruct%grid, 2), ubound(gridstruct%grid, 2), lbound(gridstruct%grid, 3), ubound(gridstruct%grid, 3)
    else
       print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") gridstruct%grid is not allocated.")', tag, this_pe, atm_n
    end if

  end subroutine show_atm_gridstruct


  subroutine show_atm(tag, Atm, atm_n, this_pe)
    implicit none
    character(len=*), intent(in)     :: tag
    type(fv_atmos_type), intent(in)  :: Atm
    integer, intent(in)              :: atm_n, this_pe

    integer is, ie, i


    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,")===============================================================")', tag, this_pe, atm_n
    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") allocated=",L1," dummy=",L1)', tag, this_pe, atm_n, Atm%allocated, Atm%dummy
    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") grid_number=",I0," ncnst=",I0," ng=",I0)', tag, this_pe, atm_n, Atm%grid_number, Atm%ncnst, Atm%ng
    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") npx=",I0," npy=",I0," npz=",I0)', tag, this_pe, atm_n, Atm%npx, Atm%npy, Atm%npz

    if (allocated(Atm%pelist)) then
       is = lbound(Atm%pelist, 1)
       ie = ubound(Atm%pelist, 1)
       print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") pelist(",I0,"-",I0,")=",I0,"...",I0)', tag, this_pe, atm_n, is, ie, Atm%pelist(is),  Atm%pelist(ie)
       !do i = is, ie
       !   print '("[INFO]    show_atm ",A8," npe=",I0," Atm(",I0,") pelist(",I0,")=",I0)', tag, this_pe, atm_n, i, Atm%pelist(i)
       !end do
    else
       print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") pelist is not allocated.")', tag, this_pe, atm_n
    end if

    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") bd%(is-ie)=",I0,"-",I0,") (js-je)=",I0,"-",I0,")"  )', tag, this_pe, atm_n, Atm%bd%is, Atm%bd%ie, Atm%bd%js, Atm%bd%je 
    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") bd%(isd-ied)=",I0,"-",I0,") (jsd-jed)=",I0,"-",I0,")"  )', tag, this_pe, atm_n, Atm%bd%isd, Atm%bd%ied, Atm%bd%jsd, Atm%bd%jed 
    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") bd%(isc-iec)=",I0,"-",I0,") (jsc-jec)=",I0,"-",I0,")"  )', tag, this_pe, atm_n, Atm%bd%isc, Atm%bd%iec, Atm%bd%jsc, Atm%bd%jec 


    call show_atm_neststruct(tag, Atm%neststruct, atm_n, this_pe)
    call show_atm_gridstruct(tag, Atm%gridstruct, atm_n, this_pe)
    call show_atm_array4(tag, Atm%grid_global, "grid_global", atm_n, this_pe)

  end subroutine show_atm




  subroutine show_gridstruct(gridstruct, this_pe)
    type(fv_grid_type), intent(in)     :: gridstruct
    integer, intent(in)                :: this_pe

    !real(kind=R_GRID)   :: pi = 4 * atan(1.0d0)
    real                 :: pi = 4 * atan(1.0d0)



    call check_3d_array(gridstruct%grid, this_pe, "SG gridstruct%grid", -2.0*pi, 2.0*pi)
    call check_3d_array(gridstruct%agrid, this_pe, "SG gridstruct%agrid", -2.0*pi, 2.0*pi)

    call check_2d_array(gridstruct%area, this_pe, "SG gridstruct%area", 0.0, 1.0e12)
    call check_2d_array(gridstruct%area_c, this_pe, "SG gridstruct%area_c", 0.0, 1.0e12)

    call check_2d_array(gridstruct%rarea, this_pe, "SG gridstruct%rarea", 0.0, 1.0e12)
    call check_2d_array(gridstruct%rarea_c, this_pe, "SG gridstruct%rarea_c", 0.0, 1.0e12)

    call check_2d_array(gridstruct%sina, this_pe, "SG gridstruct%sina",  -1.0, 1.0)
    call check_2d_array(gridstruct%cosa, this_pe, "SG gridstruct%cosa", -1.0, 1.0)

    call check_2d_array(gridstruct%dx, this_pe, "SG gridstruct%dx", 0.0, 1.0e12)
    call check_2d_array(gridstruct%dy, this_pe, "SG gridstruct%dy", 0.0, 1.0e12)

    call check_2d_array(gridstruct%dxc, this_pe, "SG gridstruct%dxc", 0.0, 1.0e12)
    call check_2d_array(gridstruct%dyc, this_pe, "SG gridstruct%dyc", 0.0, 1.0e12)

    call check_2d64_array(gridstruct%dxc_64, this_pe, "SG gridstruct%dxc_64", 0D0, 1.0D12)
    call check_2d64_array(gridstruct%dyc_64, this_pe, "SG gridstruct%dyc_64", 0D0, 1.0D12)


  end subroutine show_gridstruct




  subroutine show_nest_grid(Atm, this_pe, step_num)
    type(fv_atmos_type), intent(inout) :: Atm
    integer, intent(in)                :: this_pe, step_num

    integer             :: x,y
    integer             :: nhalo = 3  !! TODO get value from namelist
    real                :: crn_lat(4), crn_lon(4)
    real(kind=R_GRID)   :: pi = 4 * atan(1.0d0)
    real                :: pi180
    real                :: rad2deg, deg2rad

    pi180 = pi / 180.0
    deg2rad = pi / 180.0
    rad2deg = 1.0 / pi180


    print '("WDR NEST GRID bd, ",I0,",",I0," is,js=(",I0,":",I0,",",I0,":",I0,")"  )', &
         this_pe, step_num, Atm%bd%is, Atm%bd%ie, Atm%bd%js, Atm%bd%je

    print '("WDR NEST GRID bd, ",I0,",",I0," isd,jsd=(",I0,":",I0,",",I0,":",I0,")"  )', &
         this_pe, step_num, Atm%bd%isd, Atm%bd%ied, Atm%bd%jsd, Atm%bd%jed


    !do x = lbound(Atm%gridstruct%grid,1), ubound(Atm%gridstruct%grid,1)
    !   do y = lbound(Atm%gridstruct%grid,2), ubound(Atm%gridstruct%grid,2)
    !      print '("WDR NEST_GRID, ",I0,",",I0,",",I0,",",I0,",",F10.5,",",F10.5)', this_pe, step_num, x, y, &
    !           Atm%gridstruct%grid(x,y,2) * rad2deg, Atm%gridstruct%grid(x,y,1) * rad2deg - 360.0
    !   end do
    !end do

    ! Log the bounds of this PE's grid

    x = lbound(Atm%gridstruct%grid, 1)
    y = lbound(Atm%gridstruct%grid, 2)       
    crn_lon(1) = Atm%gridstruct%grid(x,y,1)
    crn_lat(1) = Atm%gridstruct%grid(x,y,2)

    x = ubound(Atm%gridstruct%grid, 1)
    y = lbound(Atm%gridstruct%grid, 2)       
    crn_lon(2) = Atm%gridstruct%grid(x,y,1)
    crn_lat(2) = Atm%gridstruct%grid(x,y,2)

    x = ubound(Atm%gridstruct%grid, 1)
    y = ubound(Atm%gridstruct%grid, 2)       
    crn_lon(3) = Atm%gridstruct%grid(x,y,1)
    crn_lat(3) = Atm%gridstruct%grid(x,y,2)

    x = lbound(Atm%gridstruct%grid, 1)
    y = ubound(Atm%gridstruct%grid, 2)       
    crn_lon(4) = Atm%gridstruct%grid(x,y,1)
    crn_lat(4) = Atm%gridstruct%grid(x,y,2)


    crn_lon(:) = crn_lon(:) * rad2deg
    crn_lat(:) = crn_lat(:) * rad2deg

    do x=1,4
       if (crn_lon(x) .gt. 180.0) then
          crn_lon(x) = crn_lon(x) - 360.0
       end if
    end do

    print '("PLOT",I0,"_data_corners,",I4.4 ,",",F10.5,",",F10.5,",",F10.5,",",F10.5,",",F10.5,",",F10.5,",",F10.5,",",F10.5)', &
         step_num, this_pe, crn_lat(1),  crn_lon(1),  crn_lat(2),  crn_lon(2),  crn_lat(3),  crn_lon(3),  crn_lat(4),  crn_lon(4)

    ! Assume that nhalo is the same as all the other halo values
    x = lbound(Atm%gridstruct%grid, 1) + nhalo
    y = lbound(Atm%gridstruct%grid, 2) + nhalo
    crn_lon(1) = Atm%gridstruct%grid(x,y,1)
    crn_lat(1) = Atm%gridstruct%grid(x,y,2)

    x = ubound(Atm%gridstruct%grid, 1) - nhalo
    y = lbound(Atm%gridstruct%grid, 2) + nhalo
    crn_lon(2) = Atm%gridstruct%grid(x,y,1)
    crn_lat(2) = Atm%gridstruct%grid(x,y,2)

    x = ubound(Atm%gridstruct%grid, 1) - nhalo
    y = ubound(Atm%gridstruct%grid, 2) - nhalo      
    crn_lon(3) = Atm%gridstruct%grid(x,y,1)
    crn_lat(3) = Atm%gridstruct%grid(x,y,2)

    x = lbound(Atm%gridstruct%grid, 1) + nhalo
    y = ubound(Atm%gridstruct%grid, 2) - nhalo
    crn_lon(4) = Atm%gridstruct%grid(x,y,1)
    crn_lat(4) = Atm%gridstruct%grid(x,y,2)


    crn_lon(:) = crn_lon(:) * rad2deg
    crn_lat(:) = crn_lat(:) * rad2deg

    do x=1,4
       if (crn_lon(x) .gt. 180.0) then
          crn_lon(x) = crn_lon(x) - 360.0
       end if
    end do

    print '("PLOT",I0,"_compute_corners,",I4.4 ,",",F10.5,",",F10.5,",",F10.5,",",F10.5,",",F10.5,",",F10.5,",",F10.5,",",F10.5)', &
         step_num, this_pe, crn_lat(1),  crn_lon(1),  crn_lat(2),  crn_lon(2),  crn_lat(3),  crn_lon(3),  crn_lat(4),  crn_lon(4)


  end subroutine show_nest_grid


  subroutine validate_hires_parent(fp_super_tile_geo, grid, agrid, x_refine, y_refine, ioffset, joffset)
    type(grid_geometry), intent(in)                  :: fp_super_tile_geo
    real, allocatable, intent(in), dimension(:,:,:)  :: grid, agrid
    integer, intent(in)                              :: x_refine, y_refine, ioffset, joffset


    real, allocatable               :: local_grid(:,:,:), local_agrid(:,:,:)
    real(kind=R_GRID), allocatable  :: local_agrid64(:,:,:)
    logical                         :: is_equal
    integer                         :: x, y, z, this_pe, stagger
    real(kind=R_GRID)               :: pi = 4 * atan(1.0d0)
    real                            :: rad2deg

    rad2deg = 180.0 / pi
    this_pe = mpp_pe()

    !! Begin test creating of grid and agrid aligned with initial nest
    !! This is for testing/validation, and will not be needed in operations

    ! Allocate grid/agrid to proper size/bounds

    allocate(local_grid(lbound(grid,1) : ubound(grid,1), &
         lbound(grid,2) : ubound(grid,2), &
         lbound(grid,3) : ubound(grid,3)))

    allocate(local_agrid(lbound(agrid,1) : ubound(agrid,1), &
         lbound(agrid,2) : ubound(agrid,2), &
         lbound(agrid,3) : ubound(agrid,3)))

    allocate(local_agrid64(lbound(agrid,1) : ubound(agrid,1), &
         lbound(agrid,2) : ubound(agrid,2), &
         lbound(agrid,3) : ubound(agrid,3)))


    ! Fill in values from high resolution, full panel, supergrid 

    stagger = CORNER
    call fill_grid_from_supergrid(local_grid, stagger, fp_super_tile_geo, ioffset, joffset, &
         x_refine, y_refine)
    stagger = CENTER
    call fill_grid_from_supergrid(local_agrid, stagger, fp_super_tile_geo, ioffset, joffset, &
         x_refine, y_refine)
    stagger = CENTER
    call fill_grid_from_supergrid(local_agrid64, stagger, fp_super_tile_geo, ioffset, joffset, &
         x_refine, y_refine)


    ! Verify that values are equivalent to the unmodified values in gridstruct

    call grid_equal(local_grid, grid, "GRID", this_pe, is_equal)
    call grid_equal(local_agrid, agrid, "AGRID", this_pe, is_equal)

    do x = lbound(grid,1), lbound(grid,1)+4
       do y = lbound(grid,2), lbound(grid,2)+4
          do z = lbound(grid,3), ubound(grid,3)
             print '("[INFO]  WDR grid_comp ",A16," npe=",I0," DEG value at (",I0,",",I0,",",I0,") ",F15.11, " ",F15.11, " ",F15.11)', "GRID", this_pe, x, y, z, local_grid(x,y,z)*rad2deg, grid(x,y,z)*rad2deg, local_grid(x,y,z)*rad2deg - grid(x,y,z)*rad2deg
             print '("[INFO]  WDR grid_comp ",A16," npe=",I0," RAD value at (",I0,",",I0,",",I0,") ",F15.11, " ",F15.11, " ",F15.11)', "GRID", this_pe, x, y, z, local_grid(x,y,z), grid(x,y,z), local_grid(x,y,z) - grid(x,y,z)
          end do
       end do
    end do

    do x = lbound(agrid,1), lbound(agrid,1)+4
       do y = lbound(agrid,2), lbound(agrid,2)+4
          do z = lbound(agrid,3), ubound(agrid,3)
             print '("[INFO]  WDR agrid_comp ",A16," npe=",I0," DEG value at (",I0,",",I0,",",I0,") ",F15.11, " ",F15.11, " ",F15.11)', "AGRID", this_pe, x, y, z, local_agrid(x,y,z)*rad2deg, agrid(x,y,z)*rad2deg, local_agrid(x,y,z)*rad2deg - agrid(x,y,z)*rad2deg
             print '("[INFO]  WDR agrid_comp ",A16," npe=",I0," RAD value at (",I0,",",I0,",",I0,") ",F15.11, " ",F15.11, " ",F15.11)', "AGRID", this_pe, x, y, z, local_agrid(x,y,z), agrid(x,y,z), local_agrid(x,y,z) - agrid(x,y,z)
          end do
       end do
    end do

    ! Validate at the end
    !call find_nest_alignment(tile_geo, fp_super_tile_geo, nest_x, nest_y, parent_x, parent_y)

  end subroutine validate_hires_parent

#endif ! MOVING_NEST

end module fv_moving_nest_logging_mod
