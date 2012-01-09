
module read_climate_nudge_data_mod

use fms_mod, only: open_namelist_file, check_nml_error, close_file, &
                   stdlog, mpp_pe, mpp_root_pe, write_version_number, &
                   string, error_mesg, FATAL, NOTE, file_exist
use mpp_mod, only: input_nml_file
use mpp_io_mod,    only: mpp_open, MPP_NETCDF, MPP_RDONLY,MPP_MULTI, MPP_SINGLE
use mpp_io_mod,    only: axistype, fieldtype, mpp_get_time_axis, mpp_get_atts
use mpp_io_mod,    only: mpp_get_fields, mpp_get_info, mpp_get_axes, mpp_get_times
use mpp_io_mod,    only: mpp_get_axis_data, mpp_read, mpp_close, mpp_get_default_calendar
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

  character(len=128) :: version = '$Id: read_climate_nudge_data.F90,v 19.0 2012/01/06 19:59:22 fms Exp $'
  character(len=128) :: tagname = '$Name: siena $'
  real, parameter :: P0 = 1.e5
  real, parameter :: D608 = RVGAS/RDGAS - 1.

  integer, parameter :: NUM_REQ_AXES = 3
  integer, parameter :: INDEX_LON = 1, INDEX_LAT = 2, INDEX_LEV = 3
  character(len=8), dimension(NUM_REQ_AXES) :: required_axis_names = &
                                     (/ 'lon', 'lat', 'lev' /)

  integer, parameter :: NUM_REQ_FLDS = 9
  integer, parameter :: INDEX_P0 = 1, INDEX_AK = 2, INDEX_BK = 3, &
                        INDEX_ZS = 4, INDEX_PS = 5,               &
                        INDEX_T  = 6, INDEX_Q  = 7,               &
                        INDEX_U  = 8, INDEX_V  = 9
  character(len=8), dimension(NUM_REQ_FLDS) :: required_field_names = &
                 (/ 'P0', 'hyai', 'hybi', 'PHI', 'PS', 'T', 'Q', 'U', 'V' /)
 
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
  integer :: ncid
  integer,  pointer :: length_axes(:)  ! length of all dimensions in file
  integer :: ndim, nvar, natt, ntim, varid_time
  integer :: time_offset
  integer, dimension(NUM_REQ_FLDS) :: field_index   ! varid for variables
  integer, dimension(NUM_REQ_AXES) :: axis_index    ! varid for dimensions
  type(axistype),  dimension(NUM_REQ_FLDS) :: axes
  type(fieldtype), dimension(NUM_REQ_FLDS) :: fields  
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
  character(len=128) :: name
  integer :: istat, i, j, k, n, nd, siz(4), i1, i2
  type(axistype), allocatable :: axes(:)
  type(fieldtype), allocatable :: fields(:)

  if (module_is_initialized) return
  ! initial file names to blanks
  do n = 1, MAXFILES
     do i = 1, len(filename_tails(n))
        filename_tails(n)(i:i) = ' '
        filenames(n)(i:i) = ' '
     enddo
  enddo

!----- read namelist -----
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=read_climate_nudge_data_nml, iostat=io)
  ierr = check_nml_error (io, 'read_climate_nudge_data_nml')
#else
  if (file_exist('input.nml') ) then
    iunit = open_namelist_file()
    ierr=1
    do while (ierr /= 0)
      read (iunit, nml=read_climate_nudge_data_nml, iostat=io, end=10)
      ierr = check_nml_error (io, 'read_climate_nudge_data_nml')
    enddo
10  call close_file (iunit)
  endif
#endif

!----- write version and namelist to log file -----

  iunit = stdlog()
  call write_version_number ( version, tagname )
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
  do n = 1, nfiles
     call mpp_open(Files(n)%ncid, trim(filenames(n)), form=MPP_NETCDF,action=MPP_RDONLY,threading=MPP_MULTI, &
             fileset=MPP_SINGLE )
!     call error_mesg ('read_climate_nudge_data_mod', 'Buffer size for reading = '//trim(string(read_buffer_size)), NOTE)

     call mpp_get_info(Files(n)%ncid, Files(n)%ndim, Files(n)%nvar, Files(n)%natt, Files(n)%ntim)

     allocate (Files(n)%length_axes(Files(n)%ndim))
     allocate (axes(Files(n)%ndim))
     call mpp_get_axes(Files(n)%ncid,axes)

     ! inquire dimension sizes
     do i = 1, Files(n)%ndim
        call mpp_get_atts(axes(i), name=name, len=Files(n)%length_axes(i))
        do j = 1, NUM_REQ_AXES
           if (trim(name) .eq. trim(required_axis_names(j))) then
              call check_axis_size (j,Files(n)%length_axes(i))
              Files(n)%axes(j) = axes(i)
              Files(n)%axis_index(j) = i
              exit
           endif
        enddo
     enddo
     deallocate(axes)
     ! time axis indexing
     Files(n)%time_offset = numtime
     numtime = numtime + Files(n)%ntim

     allocate(fields(Files(n)%nvar))
     call mpp_get_fields(Files(n)%ncid,fields)
     Files(n)%field_index = 0
     do i = 1, Files(n)%nvar
        call mpp_get_atts(fields(i), name=name, ndim=nd, siz=siz)
        do j = 1, NUM_REQ_FLDS
           if (trim(name) .eq. trim(required_field_names(j))) then
              Files(n)%field_index(j) = i
              Files(n)%fields(j) = fields(i)
              if (j .gt. 3) then
                 call check_resolution (siz(1:nd))
              endif
              exit
           endif
        enddo

        ! special case for surface geopotential (sometimes the name is PHIS)
        ! z1l: the following is not needed
!        if (trim(name) .eq. 'PHIS') then
!           Files(n)%field_index(INDEX_ZS) = i
!           call check_resolution (siz(1:nd))
!        endif
     enddo
     deallocate(fields)

  enddo ! "n" files loop

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
type(axistype), save:: time_axis
character(len=32) :: default_calendar

   if (.not.module_is_initialized) then
     call error_mesg ('read_climate_nudge_data_mod/read_time',  &
                                        'module not initialized', FATAL)
   endif

   if (size(times(:)) < numtime) then
      call error_mesg ('read_climate_nudge_data_mod', 'argument times too small in read_time', FATAL)
   endif

 ! data
   i2 = 0
   do n = 1, nfiles
      i1 = i2+1
      i2 = i2+Files(n)%ntim
      if( n == 1) then
         default_calendar = mpp_get_default_calendar()
         call mpp_get_time_axis( Files(n)%ncid, time_axis)
         call mpp_get_atts(time_axis, units=units, calendar=calendar)
         if( trim(calendar) == trim(default_calendar)) calendar = 'gregorian  '
      endif
      call mpp_get_times(Files(n)%ncid, times(i1:i2))
   enddo

! NOTE: need to do the conversion to time_type in this routine
!       this will allow different units and calendars for each file

end subroutine read_time

!###############################################################################

subroutine read_grid ( lon, lat, ak, bk )
real, intent(out), dimension(:) :: lon, lat, ak, bk

 real :: pref
 integer :: istat

   if (.not.module_is_initialized) then
     call error_mesg ('read_climate_nudge_data_mod/read_grid',  &
                                        'module not initialized', FATAL)
   endif


    ! static fields from first file only
      call mpp_get_axis_data(Files(1)%axes(INDEX_LON), lon)
      call mpp_get_axis_data(Files(1)%axes(INDEX_LAT), lat)

    ! units are assumed to be degrees east and north
    ! convert to radians
      lon = lon * PI/180.
      lat = lat * PI/180.

    ! vertical coodinate
      if (Files(1)%field_index(INDEX_AK) .gt. 0) then
         call mpp_read(Files(1)%ncid, Files(1)%fields(INDEX_AK), ak)
         if (Files(1)%field_index(INDEX_P0) .gt. 0) then
            call mpp_read(Files(1)%ncid, Files(1)%fields(INDEX_P0), ak)
         else
            pref = P0
         endif
         ak = ak*pref
      else
         ak = 0.
      endif
 
      call mpp_read(Files(1)%ncid, Files(1)%fields(INDEX_BK), bk)


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
     start(3) = atime

     nread = 1
     nread(1) = size(dat,1)
     nread(2) = size(dat,2)
     
     call mpp_read(Files(n)%ncid, Files(n)%fields(this_index), dat, start, nread)
  
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
     start(4) = atime

     nread = 1
     nread(1) = size(dat,1)
     nread(2) = size(dat,2)
     nread(3) = size(dat,3)

     call mpp_read(Files(n)%ncid, Files(n)%fields(this_index), dat, start, nread)

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
  do n = 1, nfiles
     call mpp_close(Files(n)%ncid)
  enddo
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

 subroutine check_resolution (axis_len)
 integer, intent(in) :: axis_len(:)

   if (size(axis_len(:)) .lt. 2) then
      call error_mesg ('read_climate_nudge_data_mod', 'incorrect number of array dimensions', FATAL)
   endif
   if (axis_len(1) .ne. global_axis_size(INDEX_LON)) then
      call error_mesg ('read_climate_nudge_data_mod', 'incorrect array dimension one', FATAL)
   endif
   if (axis_len(2) .ne. global_axis_size(INDEX_LAT)) then
      call error_mesg ('read_climate_nudge_data_mod', 'incorrect array dimension two', FATAL)
   endif
   if (size(axis_len(:)) .gt. 3) then
      if (axis_len(3) .ne. global_axis_size(INDEX_LEV)) then
         call error_mesg ('read_climate_nudge_data_mod', 'incorrect array dimension three', FATAL)
      endif
   endif

 end subroutine check_resolution

!###############################################################################

end module read_climate_nudge_data_mod

