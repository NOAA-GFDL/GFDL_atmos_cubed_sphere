
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
module read_climate_nudge_data_mod

use fms_mod, only: check_nml_error, &
                   stdlog, mpp_pe, mpp_root_pe, write_version_number, &
                   string, error_mesg, FATAL, NOTE
use fms2_io_mod,   only: open_file, close_file, get_num_dimensions, &
                         get_dimension_names, get_dimension_size, FmsNetcdfFile_t, &
                         get_num_variables, get_variable_names, get_variable_size, &
                         get_variable_num_dimensions, get_variable_units, &
                         get_time_calendar, read_data, variable_att_exists, &
                         is_dimension_unlimited
use mpp_mod,       only: input_nml_file, mpp_npes, mpp_get_current_pelist
use constants_mod, only: PI, GRAV, RDGAS, RVGAS

implicit none
private

public :: read_climate_nudge_data_init, read_time, read_grid,  &
          read_climate_nudge_data, read_climate_nudge_data_end
public :: read_sub_domain_init

interface read_climate_nudge_data
   module procedure read_climate_nudge_data_2d
   module procedure read_climate_nudge_data_3d
end interface

! version number of this module
! Include variable "version" to be written to log file.
#include<file_version.h>

  real, parameter :: P0 = 1.e5
  real, parameter :: D608 = RVGAS/RDGAS - 1.

  integer, parameter :: NUM_REQ_AXES = 3
  integer, parameter :: INDEX_LON = 1, INDEX_LAT = 2, INDEX_LEV = 3, INDEX_TIME = 4
  character(len=8), dimension(NUM_REQ_AXES) :: required_axis_names = &
                                     (/ 'lon', 'lat', 'lev' /)

  integer, parameter :: NUM_REQ_FLDS = 9
  integer, parameter :: INDEX_P0 = 1, INDEX_AK = 2, INDEX_BK = 3, &
                        INDEX_ZS = 4, INDEX_PS = 5,               &
                        INDEX_T  = 6, INDEX_Q  = 7,               &
                        INDEX_U  = 8, INDEX_V  = 9
  character(len=8), dimension(NUM_REQ_FLDS) :: required_field_names = &
       (/ 'P0  ', 'hyai', 'hybi', 'PHI ', 'PS  ', 'T   ', 'Q   ', 'U   ', 'V   ' /)

  integer, parameter :: MAXFILES = 53
  character(len=256) :: filenames(MAXFILES)
  character(len=256) :: filename_tails(MAXFILES)
  character(len=256) :: filename_head
  integer :: read_buffer_size
  integer :: nfiles = 0
  logical :: module_is_initialized = .false.

  namelist /read_climate_nudge_data_nml/ filename_tails, read_buffer_size, &
                                         filename_head

! dimensions for checking
  integer :: global_axis_size(NUM_REQ_AXES), numtime, sub_domain_latitude_size
  integer, allocatable :: file_index(:)

type filedata_type
  integer, pointer :: length_axes(:)
  integer :: ndim, nvar, ntim
  integer :: time_offset
  integer, dimension(NUM_REQ_FLDS) :: field_index   ! varid for variables
  integer, dimension(NUM_REQ_AXES) :: axis_index    ! varid for dimensions
  character(len=32), dimension(NUM_REQ_FLDS+1) :: axes
  character(len=32), dimension(NUM_REQ_FLDS) :: fields
end type

  type(filedata_type), allocatable :: Files(:)

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_climate_nudge_data_init (nlon, nlat, nlev, ntime)
integer, intent(out) :: nlon, nlat, nlev, ntime
! returns dimension lengths of input data set
!   nlon, nlat  lat/lon grid size
!   nlev        number of levels
!   ntime       number of time levels

  integer :: iunit, ierr, io
  integer :: istat, i, j, k, n, i1, i2
  character(len=32), allocatable :: fields(:)
  type(FmsNetcdfFile_t) :: fileobj
  integer, allocatable, dimension(:) :: pes !< Array of ther pes in the current pelist
  character(len=128), dimension(:), allocatable :: names

  if (module_is_initialized) return
  ! initial file names to blanks
  do n = 1, MAXFILES
     do i = 1, len(filename_tails(n))
        filename_tails(n)(i:i) = ' '
        filenames(n)(i:i) = ' '
     enddo
  enddo

!----- read namelist -----
  read (input_nml_file, nml=read_climate_nudge_data_nml, iostat=io)
  ierr = check_nml_error (io, 'read_climate_nudge_data_nml')

!----- write version and namelist to log file -----

  iunit = stdlog()
  call write_version_number ( 'READ_CLIMATE_NUDGE_DATA_MOD', version )
  if (mpp_pe() == mpp_root_pe()) write (iunit, nml=read_climate_nudge_data_nml)

  ! determine the number of files
  do n = 1, MAXFILES
     if (filename_tails(n)(1:1) .eq. ' ') exit
     nfiles = n
  enddo
  do n=1,nfiles
    filenames(n) = trim(filename_head)//trim(filename_tails(n))
  end do

  allocate(Files(nfiles))
  numtime = 0

! open input file(s)
  allocate(pes(mpp_npes()))
  call mpp_get_current_pelist(pes)

  do n = 1, nfiles
     if (open_file(fileobj, trim(filenames(n)), "read", pelist=pes)) then
         Files(n)%ndim=get_num_dimensions(fileobj)

         allocate(files(n)%length_axes(Files(n)%ndim))
         allocate(names(files(n)%ndim))
         call get_dimension_names(fileobj, names)
         files(n)%axes(:) = ""

         ! inquire dimension sizes
         do i = 1, Files(n)%ndim
            do j = 1, NUM_REQ_AXES
               if (trim(names(i)) .eq. trim(required_axis_names(j))) then
                  files(n)%axes(j) = trim(names(i))
                  call get_dimension_size(fileobj, Files(n)%axes(j), Files(n)%length_axes(j))
                  call check_axis_size (j,Files(n)%length_axes(j))
                  Files(n)%axis_index(j) = i
                  exit
               endif
            enddo
            if (j .gt. num_req_axes) then
              if (is_dimension_unlimited(fileobj, trim(names(i)))) then
                ! time axis indexing
                files(n)%axes(index_time) = trim(names(i))
                call get_dimension_size(fileobj, files(n)%axes(index_time), &
                                        files(n)%length_axes(index_time))
              endif
            endif
         enddo
         deallocate(names)
         Files(n)%ntim = Files(n)%length_axes(INDEX_TIME)
         Files(n)%time_offset = numtime
         numtime = numtime + Files(n)%ntim

         Files(n)%nvar = get_num_variables(fileobj)
         allocate(fields(Files(n)%nvar))
         call get_variable_names(fileobj, fields)
         Files(n)%field_index = 0
         do i = 1, Files(n)%nvar
            do j = 1, NUM_REQ_FLDS
               if (trim(fields(i)) .eq. trim(required_field_names(j))) then
                  Files(n)%field_index(j) = i
                  Files(n)%fields(j) = fields(i)
                  if (j .gt. 3) then
                     call check_resolution (fileobj, fields(i))
                  endif
                  exit
               endif
            enddo

            ! special case for surface geopotential (sometimes the name is PHIS)
            if (trim(fields(i)) .eq. 'PHIS') then
               Files(n)%field_index(INDEX_ZS) = i
               Files(n)%fields(INDEX_ZS) = fields(i)
               call check_resolution (fileobj, fields(i))
            endif
         enddo
         deallocate(fields)
         call close_file(fileobj)
     endif
  enddo ! "n" files loop
  deallocate(pes)

  ! setup file indexing
  allocate(file_index(numtime))
  i2 = 0
  do n = 1, nfiles
     i1 = i2+1
     i2 = i2+Files(n)%ntim
     file_index(i1:i2) = n
  enddo

      sub_domain_latitude_size = global_axis_size(INDEX_LAT)

    ! output arguments
      nlon = global_axis_size(INDEX_LON)
      nlat = global_axis_size(INDEX_LAT)
      nlev = global_axis_size(INDEX_LEV)
      ntime = numtime

      module_is_initialized = .true.

end subroutine read_climate_nudge_data_init

!###############################################################################

subroutine read_time ( times, units, calendar )
real(8),          intent(out) :: times(:)
character(len=*), intent(out) :: units, calendar
integer :: istat, i1, i2, n
real :: l_times(size(times))
character(len=32), save:: time_axis
type(FmsNetcdfFile_t) :: fileobj
integer, allocatable, dimension(:) :: pes !< Array of ther pes in the current pelist

   if (.not.module_is_initialized) then
     call error_mesg ('read_climate_nudge_data_mod/read_time',  &
                                        'module not initialized', FATAL)
   endif

   if (size(times(:)) < numtime) then
      call error_mesg ('read_climate_nudge_data_mod', 'argument times too small in read_time', FATAL)
   endif

 ! data
   i2 = 0
   allocate(pes(mpp_npes()))
   call mpp_get_current_pelist(pes)
   do n = 1, nfiles
      if (open_file(fileobj, trim(filenames(n)), "read", pelist=pes)) then
          i1 = i2+1
          i2 = i2+Files(n)%ntim
          if( n == 1) then
             time_axis = Files(n)%axes(INDEX_TIME)
             call get_variable_units(fileobj, time_axis, units)
             if( variable_att_exists(fileobj, time_axis, "calendar") ) then
                 call get_time_calendar(fileobj, time_axis, calendar)
             else
                 calendar = 'gregorian  '
             endif
          endif
          call read_data(fileobj,time_axis, l_times(i1:i2))
          times(i1:i2) = l_times(i1:i2)
          call close_file(fileobj)
      endif
   enddo
   deallocate(pes)

! NOTE: need to do the conversion to time_type in this routine
!       this will allow different units and calendars for each file

end subroutine read_time

!###############################################################################

subroutine read_grid ( lon, lat, ak, bk )
real, intent(out), dimension(:) :: lon, lat, ak, bk

 real :: pref
 integer :: istat
 type(FmsNetcdfFile_t) :: fileobj
 integer, allocatable, dimension(:) :: pes !< Array of ther pes in the current pelist

   if (.not.module_is_initialized) then
     call error_mesg ('read_climate_nudge_data_mod/read_grid',  &
                                        'module not initialized', FATAL)
   endif


    ! static fields from first file only
      allocate(pes(mpp_npes()))
      call mpp_get_current_pelist(pes)
      if (open_file(fileobj, trim(filenames(1)), "read", pelist=pes)) then
          call read_data(fileobj, Files(1)%axes(INDEX_LON), lon)
          call read_data(fileobj, Files(1)%axes(INDEX_LAT), lat)

        ! units are assumed to be degrees east and north
        ! convert to radians
          lon = lon * PI/180.
          lat = lat * PI/180.

        ! vertical coodinate
          if (Files(1)%field_index(INDEX_AK) .gt. 0) then
             call read_data(fileobj, Files(1)%fields(INDEX_AK), ak)
             if (Files(1)%field_index(INDEX_P0) .gt. 0) then
                call read_data(fileobj, Files(1)%fields(INDEX_P0), pref)
             else
                pref = P0
             endif
             ak = ak*pref
          else
             ak = 0.
          endif

          call read_data(fileobj, Files(1)%fields(INDEX_BK), bk)
          call close_file(fileobj)
      endif
      deallocate(pes)

end subroutine read_grid

!###############################################################################

subroutine read_sub_domain_init ( ylo, yhi, ydat, js, je )
 real,    intent(in)  :: ylo, yhi, ydat(:)
 integer, intent(out) :: js, je
 integer :: j

   if (.not.module_is_initialized) then
     call error_mesg ('read_climate_nudge_data_mod/read_sub_domain_init',  &
                                        'module not initialized', FATAL)
   endif
   ! increasing data
   if (ydat(1) < ydat(2)) then
      js = 1
      do j = 1, size(ydat(:))-1
         if (ylo >= ydat(j) .and. ylo <= ydat(j+1)) then
            js = j
            exit
         endif
      enddo

      if (ylo < -1.5) then
         print *, 'ylo=',ylo
         print *, 'js,ydat=',js,ydat(js)
         print *, 'ydat=',ydat(:js+2)
      endif

      je = size(ydat(:))
      do j = js, size(ydat(:))-1
         if (yhi >= ydat(j) .and. yhi <= ydat(j+1)) then
            je = j+1
            exit
         endif
      enddo

      if (yhi > 1.5) then
         print *, 'yhi=',yhi
         print *, 'je,ydat=',je,ydat(je)
         print *, 'ydat=',ydat(je-2:)
      endif

   ! decreasing data (may not work)
   else
      call error_mesg ('read_climate_nudge_data_mod', 'latitude values for observational data decrease with increasing index', NOTE)
      je = size(ydat(:))-1
      do j = 1, size(ydat(:))-1
         if (ylo >= ydat(j+1) .and. ylo <= ydat(j)) then
            je = j+1
            exit
         endif
      enddo

      js = 1
      do j = 1, je
         if (yhi >= ydat(j+1) .and. yhi <= ydat(j)) then
            js = j
            exit
         endif
      enddo

   endif

   sub_domain_latitude_size = je-js+1

 end subroutine read_sub_domain_init

!###############################################################################

subroutine read_climate_nudge_data_2d (itime, field, dat, is, js)
integer,          intent(in) :: itime
character(len=4), intent(in) :: field
real,             intent(out), dimension(:,:) :: dat
integer,          intent(in),  optional       :: is, js
integer :: istat, atime, n, this_index
integer :: nread(4), start(4)
type(FmsNetcdfFile_t) :: fileobj
integer, allocatable, dimension(:) :: pes !< Array of ther pes in the current pelist

   if (.not.module_is_initialized) then
     call error_mesg ('read_climate_nudge_data_mod',  &
                                        'module not initialized', FATAL)
   endif
     ! time index check
      if (itime < 1 .or. itime > numtime) then
         call error_mesg ('read_climate_nudge_data_mod', 'itime out of range', FATAL)
      endif

     ! check dimensions
     if (present(js)) then
        if (size(dat,1) .ne. global_axis_size(INDEX_LON) .or. &
            size(dat,2) .ne. sub_domain_latitude_size) then
            !write (*,'(a)') 'climate_nudge_data_mod: size dat2d = '//trim(string(size(dat,1)))//' x '//trim(string(size(dat,2)))// &
            !         '  <-vs->  '//trim(string(global_axis_size(INDEX_LON)))//' x '//trim(string(sub_domain_latitude_size))
            call error_mesg ('read_climate_nudge_data_mod', 'incorrect 2d array dimensions', FATAL)
        endif
     else
        if (size(dat,1) .ne. global_axis_size(INDEX_LON) .or. &
            size(dat,2) .ne. global_axis_size(INDEX_LAT))     &
            call error_mesg ('read_climate_nudge_data_mod', 'incorrect 2d array dimensions', FATAL)
     endif

     ! check field
     if (field .eq. 'phis') then
        this_index = INDEX_ZS
     else if (field .eq. 'psrf') then
        this_index = INDEX_PS
     else
         call error_mesg ('read_climate_nudge_data_mod', 'incorrect field requested in read_climate_nudge_data_2d', FATAL)
     endif

     ! file index and actual time index in file
     n = file_index(itime)
     atime = itime - Files(n)%time_offset

     start = 1
     if (present(is)) start(1) = is
     if (present(js)) start(2) = js

     nread(1) = size(dat,1)
     nread(2) = size(dat,2)

     allocate(pes(mpp_npes()))
     call mpp_get_current_pelist(pes)
     if (open_file(fileobj, trim(filenames(n)), "read", pelist=pes)) then
         call read_data(fileobj, Files(n)%fields(this_index), dat, unlim_dim_level=atime, &
                        corner=start(1:2), edge_lengths=nread(1:2))
         call close_file(fileobj)
     endif
     deallocate(pes)

      ! geopotential height (convert to m2/s2 if necessary)
     if (field .eq. 'phis') then
        if (maxval(dat) > 1000.*GRAV) then
          ! do nothing
        else
          dat = dat * GRAV
        endif
     endif

end subroutine read_climate_nudge_data_2d

!###############################################################################

subroutine read_climate_nudge_data_3d (itime, field, dat, is, js)
integer,          intent(in) :: itime
character(len=4), intent(in) :: field
real,             intent(out), dimension(:,:,:) :: dat
integer,          intent(in),  optional         :: is, js
integer :: istat, atime, n, this_index, start(4), nread(4)
type(FmsNetcdfFile_t) :: fileobj
integer, allocatable, dimension(:) :: pes !< Array of ther pes in the current pelist
!logical :: convert_virt_temp = .false.

   if (.not.module_is_initialized) then
     call error_mesg ('read_climate_nudge_data_mod',  &
                                        'module not initialized', FATAL)
   endif

     ! time index check
     if (itime < 1 .or. itime > numtime) then
        call error_mesg ('read_climate_nudge_data_mod', 'itime out of range', FATAL)
     endif

     ! check dimensions
     if (present(js)) then
        if (size(dat,1) .ne. global_axis_size(INDEX_LON) .or. &
            size(dat,2) .ne. sub_domain_latitude_size    .or. &
            size(dat,3) .ne. global_axis_size(INDEX_LEV)) then
            !write (*,'(a)') 'climate_nudge_data_mod: size dat3d = '//trim(string(size(dat,1)))//' x '//trim(string(size(dat,2)))// &
            !                                        ' x '//trim(string(size(dat,3)))
            call error_mesg ('read_climate_nudge_data_mod', 'incorrect 3d array dimensions', FATAL)
        endif
     else
        if (size(dat,1) .ne. global_axis_size(INDEX_LON) .or. &
            size(dat,2) .ne. global_axis_size(INDEX_LAT) .or. &
            size(dat,3) .ne. global_axis_size(INDEX_LEV))     &
            call error_mesg ('read_climate_nudge_mod', 'incorrect 3d array dimensions', FATAL)
     endif

     ! check field
     if (field .eq. 'temp') then
        this_index = INDEX_T
     else if (field .eq. 'qhum') then
        this_index = INDEX_Q
     else if (field .eq. 'uwnd') then
        this_index = INDEX_U
     else if (field .eq. 'vwnd') then
        this_index = INDEX_V
     else
        call error_mesg ('read_climate_nudge_data_mod', 'incorrect field requested in read_climate_nudge_data_3d', FATAL)
     endif


     ! file index and actual time index in file
     n = file_index(itime)
     atime = itime - Files(n)%time_offset

     start = 1
     if (present(is)) start(1) = is
     if (present(js)) start(2) = js

     nread(1) = size(dat,1)
     nread(2) = size(dat,2)
     nread(3) = size(dat,3)

     allocate(pes(mpp_npes()))
     call mpp_get_current_pelist(pes)
     if (open_file(fileobj, trim(filenames(n)), "read", pelist=pes)) then
         call read_data(fileobj, Files(n)%fields(this_index), dat, unlim_dim_level=atime, &
                        corner=start(1:3), edge_lengths=nread(1:3))
         call close_file(fileobj)
     endif
     deallocate(pes)

     ! convert virtual temp to temp
     ! necessary for some of the high resol AVN analyses
     !if (convert_virt_temp) then
     !   temp = temp/(1.+D608*qhum)
     !endif

end subroutine read_climate_nudge_data_3d

!###############################################################################

subroutine read_climate_nudge_data_end
integer :: istat, n

  if ( .not.module_is_initialized) return
  deallocate (Files)
  module_is_initialized = .false.

end subroutine read_climate_nudge_data_end

!###############################################################################

 subroutine check_axis_size (ind,lendim)
 integer, intent(in) :: ind,lendim

   ! once the axis size is set all subsuquent axes must be the same
   if (global_axis_size(ind) .gt. 0) then
      if (global_axis_size(ind) .ne. lendim) then
         call error_mesg ('read_climate_nudge_data_mod', 'incorrect axis size for axis = '//trim(required_axis_names(ind)), FATAL)
      endif
   else
      global_axis_size(ind) = lendim
   endif

 end subroutine check_axis_size

!------------------------------------

 subroutine check_resolution (fileobj, field_name)

 type(FmsNetcdfFile_t), intent(in) :: fileobj
 character(len=*), intent(in) :: field_name

 integer :: nd
 integer, allocatable :: siz(:)

   nd = get_variable_num_dimensions(fileobj, field_name)

   if (nd .lt. 2) then
      call error_mesg ('read_climate_nudge_data_mod', 'incorrect number of array dimensions', FATAL)
   endif

   allocate(siz(nd))
   call get_variable_size(fileobj, field_name, siz)

   if (siz(1) .ne. global_axis_size(INDEX_LON)) then
      call error_mesg ('read_climate_nudge_data_mod', 'incorrect array dimension one', FATAL)
   endif
   if (siz(2) .ne. global_axis_size(INDEX_LAT)) then
      call error_mesg ('read_climate_nudge_data_mod', 'incorrect array dimension two', FATAL)
   endif
   if (nd .gt. 3) then
      if (siz(3) .ne. global_axis_size(INDEX_LEV)) then
         call error_mesg ('read_climate_nudge_data_mod', 'incorrect array dimension three', FATAL)
      endif
   endif
   deallocate(siz)

 end subroutine check_resolution

!###############################################################################

end module read_climate_nudge_data_mod

