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

module statistics_mod

implicit none

interface mode
   module procedure mode_1d_real4
   module procedure mode_2d_real4
   module procedure masked_mode_2d_real4
   module procedure mode_1d_real8
   module procedure mode_2d_real8
   module procedure masked_mode_2d_real8
end interface mode

contains

  ! qksrt implementation copied and adapted for real arrays from implementation
  ! in FMS: FMS/drifters/quicksort.F90
  function qksrt_partition_real4(n, list, start, end) result(top)
    implicit none
    integer, intent(in) :: n
    real(kind=4), intent(inout) :: list(n)
    integer, intent(in) :: start, end

    real(kind=4) :: pivot
    integer :: bottom, top
    logical :: done

    pivot = list(end)                          ! Partition around the last value
    bottom = start-1                           ! Start outside the area to be partitioned
    top = end                                  ! Ditto

    done = .false.
    do while (.not. done)                      ! Until all elements are partitioned...

       do while (.not. done)                  ! Until we find an out of place element...
          bottom = bottom+1                  ! ... move the bottom up.

          if(bottom == top) then             ! If we hit the top...
             done = .true.                  ! ... we are done.
             exit
          endif

          if(list(bottom) > pivot) then           ! Is the bottom out of place?
             list(top) = list(bottom)       ! Then put it at the top...
              exit                          ! ... and start searching from the top.
           endif
        enddo

        do while (.not. done)                        ! Until we find an out of place element...
           top = top-1                        ! ... move the top down.

           if(top == bottom) then                  ! If we hit the bottom...
              done = .true.                      ! ... we are done.
              exit
           endif

           if(list(top) < pivot) then              ! Is the top out of place?
              list(bottom) = list(top)       ! Then put it at the bottom...
              exit                          ! ...and start searching from the bottom.
           endif
        enddo
    enddo

    list(top) = pivot                          ! Put the pivot in its place.
    ! Return the split point
  end function qksrt_partition_real4

  recursive subroutine qksrt_quicksort_real4(n, list, start, end)
    implicit none
    integer, intent(in) :: n
    real(kind=4), intent(inout) :: list(n)
    integer, intent(in) :: start, end
    integer :: split

    if(start < end) then                            ! If there are two or more elements...
        split = qksrt_partition_real4(n, list, start, end)    ! ... partition the sublist...
        call qksrt_quicksort_real4(n, list,  start, split-1)        ! ... and sort both halves.
        call qksrt_quicksort_real4(n, list, split+1, end)
    endif
  end subroutine qksrt_quicksort_real4

  ! This procedure produces the same results as scipy.stats.mode; if there is a
  ! tie in counts, the minimum mode value is returned.
  function mode_1d_real4(array)
    real(kind=4), dimension(:), intent(in) :: array

    real(kind=4) :: mode_1d_real4

    integer :: i, run, max_run
    real(kind=4), dimension(size(array)) :: sorted_array

    run = 1
    max_run = 0

    sorted_array = array
    call qksrt_quicksort_real4(size(sorted_array), sorted_array, 1, size(sorted_array))

    if (size(sorted_array) == 1) then
       mode_1d_real4 = sorted_array(1)
    else
       do i = 2, size(sorted_array)
          if (sorted_array(i) == sorted_array(i - 1)) then
             run = run + 1
          else
             run = 1
          endif
          if (run > max_run) then
             max_run = run
             mode_1d_real4 = sorted_array(i - 1)
          endif
       enddo
    endif
  end function mode_1d_real4

  function mode_2d_real4(array)
    real(kind=4), dimension(:,:), intent(in) :: array

    real(kind=4) :: mode_2d_real4

    mode_2d_real4 = mode_1d_real4(pack(array, .true.))
  end function mode_2d_real4

  function masked_mode_2d_real4(array, mask)
    real(kind=4), dimension(:,:), intent(in) :: array
    logical, dimension(:,:), intent(in) :: mask
    real(kind=4) :: masked_mode_2d_real4

    masked_mode_2d_real4 = mode_1d_real4(pack(array, mask))
  end function masked_mode_2d_real4

  ! qksrt implementation copied and adapted for real arrays from implementation
  ! in FMS: FMS/drifters/quicksort.F90
  function qksrt_partition_real8(n, list, start, end) result(top)
    implicit none
    integer, intent(in) :: n
    real(kind=8), intent(inout) :: list(n)
    integer, intent(in) :: start, end

    real(kind=8) :: pivot
    integer :: bottom, top
    logical :: done

    pivot = list(end)                          ! Partition around the last value
    bottom = start-1                           ! Start outside the area to be partitioned
    top = end                                  ! Ditto

    done = .false.
    do while (.not. done)                      ! Until all elements are partitioned...

       do while (.not. done)                  ! Until we find an out of place element...
          bottom = bottom+1                  ! ... move the bottom up.

          if(bottom == top) then             ! If we hit the top...
             done = .true.                  ! ... we are done.
             exit
          endif

          if(list(bottom) > pivot) then           ! Is the bottom out of place?
             list(top) = list(bottom)       ! Then put it at the top...
              exit                          ! ... and start searching from the top.
           endif
        enddo

        do while (.not. done)                        ! Until we find an out of place element...
           top = top-1                        ! ... move the top down.

           if(top == bottom) then                  ! If we hit the bottom...
              done = .true.                      ! ... we are done.
              exit
           endif

           if(list(top) < pivot) then              ! Is the top out of place?
              list(bottom) = list(top)       ! Then put it at the bottom...
              exit                          ! ...and start searching from the bottom.
           endif
        enddo
    enddo

    list(top) = pivot                          ! Put the pivot in its place.
    ! Return the split point
  end function qksrt_partition_real8

  recursive subroutine qksrt_quicksort_real8(n, list, start, end)
    implicit none
    integer, intent(in) :: n
    real(kind=8), intent(inout) :: list(n)
    integer, intent(in) :: start, end
    integer :: split

    if(start < end) then                            ! If there are two or more elements...
        split = qksrt_partition_real8(n, list, start, end)    ! ... partition the sublist...
        call qksrt_quicksort_real8(n, list,  start, split-1)        ! ... and sort both halves.
        call qksrt_quicksort_real8(n, list, split+1, end)
    endif
  end subroutine qksrt_quicksort_real8

  ! This procedure produces the same results as scipy.stats.mode; if there is a
  ! tie in counts, the minimum mode value is returned.
  function mode_1d_real8(array)
    real(kind=8), dimension(:), intent(in) :: array

    real(kind=8) :: mode_1d_real8

    integer :: i, run, max_run
    real(kind=8), dimension(size(array)) :: sorted_array

    run = 1
    max_run = 0

    sorted_array = array
    call qksrt_quicksort_real8(size(sorted_array), sorted_array, 1, size(sorted_array))

    if (size(sorted_array) == 1) then
       mode_1d_real8 = sorted_array(1)
    else
       do i = 2, size(sorted_array)
          if (sorted_array(i) == sorted_array(i - 1)) then
             run = run + 1
          else
             run = 1
          endif
          if (run > max_run) then
             max_run = run
             mode_1d_real8 = sorted_array(i - 1)
          endif
       enddo
    endif
  end function mode_1d_real8

  function mode_2d_real8(array)
    real(kind=8), dimension(:,:), intent(in) :: array

    real(kind=8) :: mode_2d_real8

    mode_2d_real8 = mode_1d_real8(pack(array, .true.))
  end function mode_2d_real8

  function masked_mode_2d_real8(array, mask)
    real(kind=8), dimension(:,:), intent(in) :: array
    logical, dimension(:,:), intent(in) :: mask
    real(kind=8) :: masked_mode_2d_real8

    masked_mode_2d_real8 = mode_1d_real8(pack(array, mask))
  end function masked_mode_2d_real8
end module statistics_mod
