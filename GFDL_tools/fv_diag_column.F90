module fv_diag_column_mod

  use fv_arrays_mod,      only: fv_atmos_type, fv_grid_type, fv_diag_type, fv_grid_bounds_type, &
                                R_GRID
  use fv_grid_utils_mod,  only: great_circle_dist
  use time_manager_mod,   only: time_type, get_date, get_time, month_name
  use constants_mod,      only: grav, rdgas, kappa, cp_air, TFREEZE, pi=>pi_8
  use fms_mod,            only: write_version_number, lowercase
  use mpp_mod,            only: mpp_error, FATAL, stdlog, mpp_pe, mpp_root_pe, mpp_sum, &
                                mpp_max, NOTE, input_nml_file, get_unit
  use mpp_io_mod,         only: mpp_flush
  use fv_sg_mod,          only: qsmith

  implicit none
  private

   integer, parameter :: MAX_DIAG_COLUMN = 100
 integer, parameter :: diag_name_len = 16
 integer, allocatable, dimension(:) :: diag_debug_units
 character(diag_name_len), dimension(MAX_DIAG_COLUMN) :: diag_debug_names
 real, dimension(MAX_DIAG_COLUMN) :: diag_debug_lon, diag_debug_lat

 integer, allocatable, dimension(:) :: diag_sonde_units
 character(diag_name_len), dimension(MAX_DIAG_COLUMN) :: diag_sonde_names
 real, dimension(MAX_DIAG_COLUMN) :: diag_sonde_lon, diag_sonde_lat
 integer, dimension(MAX_DIAG_COLUMN) :: diag_debug_i, diag_debug_j, diag_debug_tile
 integer, dimension(MAX_DIAG_COLUMN) :: diag_sonde_i, diag_sonde_j, diag_sonde_tile

 logical :: do_diag_debug = .false. !< Whether to enable "diagnostic" debug columns, read from column_table
 logical :: do_diag_debug_dyn = .false. !< Whether to write out debug columns every acoustic timestep instead of just every fv_diag timestep. Requires do_diag_debug to be .true.
 logical :: do_diag_sonde = .false. !< Whether to enable point (column) sounding output, in the University of Wyoming format, read from column_table
 integer :: sound_freq = 3 !< frequency (in hours) to write out diagnostic soundings
 integer :: num_diag_debug = 0
 integer :: num_diag_sonde = 0
 character(100) :: runname = 'test' !< Name for this run, used in sonde output
 integer :: diag_debug_kbottom !< bottom level (noting k=1 is the top) of diagnostic debug output. Used to limit the copious diagnostic sounding output to the layers of interest. Default is npz.
 integer :: diag_debug_nlevels !< number of levels, counting upwards (to smaller k) from diag_debug_kbottom of diagnostic debug output. Default is npz/3.

 character(10) :: init_str
 real, parameter    ::     rad2deg = 180./pi

 public :: do_diag_debug_dyn, debug_column, debug_column_dyn, fv_diag_column_init, sounding_column


 namelist /fv_diag_column_nml/ do_diag_debug, do_diag_debug_dyn, do_diag_sonde, &
      sound_freq, runname, diag_debug_kbottom, diag_debug_nlevels

! version number of this module
! Include variable "version" to be written to log file.
#include<file_version.h>

contains

 subroutine fv_diag_column_init(Atm, yr_init, mo_init, dy_init, hr_init, do_diag_debug_out, do_diag_sonde_out, sound_freq_out)

   type(fv_atmos_type), intent(inout), target :: Atm
   integer, intent(IN) :: yr_init, mo_init, dy_init, hr_init
   logical, intent(OUT) :: do_diag_debug_out, do_diag_sonde_out
   integer, intent(OUT) :: sound_freq_out

   integer :: ios, nlunit
   logical :: exists

   call write_version_number ( 'FV_DIAG_COLUMN_MOD', version )

    diag_debug_names(:) = ''
    diag_debug_lon(:) = -999.
    diag_debug_lat(:) = -999.
    diag_debug_i(:) = -999
    diag_debug_j(:) = -999
    diag_debug_tile(:) = -999
    diag_debug_kbottom = Atm%npz
    diag_debug_nlevels = Atm%npz/3

    diag_sonde_names(:) = ''
    diag_sonde_lon(:) = -999.
    diag_sonde_lat(:) = -999.
    diag_sonde_i(:) = -999
    diag_sonde_j(:) = -999
    diag_sonde_tile(:) = -99

#ifdef INTERNAL_FILE_NML
    read(input_nml_file, nml=fv_diag_column_nml,iostat=ios)
#else
    inquire (file=trim(Atm%nml_filename), exist=exists)
    if (.not. exists) then
      write(errmsg,*) 'fv_diag_column_nml: namelist file ',trim(Atm%nml_filename),' does not exist'
      call mpp_error(FATAL, errmsg)
    else
      open (unit=nlunit, file=Atm%nml_filename, READONLY, status='OLD', iostat=ios)
    endif
    rewind(nlunit)
    read (nlunit, nml=fv_diag_column_nml, iostat=ios)
    close (nlunit)
#endif

    if (do_diag_debug .or. do_diag_sonde) then
       call read_column_table
    endif

    if (do_diag_debug) then
       allocate(diag_debug_units(num_diag_debug))
       call find_diagnostic_column("DEBUG", diag_debug_names, diag_debug_i, diag_debug_j, diag_debug_tile, diag_debug_lat, diag_debug_lon, diag_debug_units, Atm%gridstruct%grid_64, Atm%gridstruct%agrid_64, num_diag_debug, Atm%gridstruct%ntiles_g, Atm%bd, Atm%global_tile, Atm%npx, Atm%npy)
    endif
    if (do_diag_sonde) then
       allocate(diag_sonde_units(num_diag_sonde))
       call find_diagnostic_column("Sonde ", diag_sonde_names, diag_sonde_i, diag_sonde_j, diag_sonde_tile, diag_sonde_lat, diag_sonde_lon, diag_sonde_units, Atm%gridstruct%grid_64, Atm%gridstruct%agrid_64, num_diag_sonde, Atm%gridstruct%ntiles_g, Atm%bd, Atm%global_tile, Atm%npx, Atm%npy)
    endif

    write(init_str,400) yr_init, mo_init, dy_init, hr_init
400 format(I4, I2.2, I2.2, I2.2 )

    do_diag_debug_out = do_diag_debug
    do_diag_sonde_out = do_diag_sonde
    sound_freq_out    = sound_freq


 end subroutine fv_diag_column_init


!-----------------------------------------------------------------------
!use diag_debug_[ij] for everything

  subroutine read_column_table
!< EXAMPLE COLUMN_TABLE file:
!< #Use space-delineated fields (no commas)
!< DEBUG index  ORD  2 30 5
!< DEBUG index  Princeton 2 37 5
!< DEBUG lonlat ORD2 272. 42.
!< DEBUG lonlat Princeton 285.33 40.36
!< DEBUG lonlat NP 0. 90.
!< DEBUG lonlat SP 0. -90.
!< sonde lonlat OUN          -97.47 35.22
!< sonde lonlat Amarillo    -101.70 35.22
!< sonde lonlat DelRio      -100.92 29.37
!< sonde lonlat Jackson      -90.08 32.32
!< sonde lonlat ILX          -89.34 40.15
!< sonde lonlat AtlanticCity -74.56 39.45
!< sonde lonlat DodgeCity    -99.97 37.77

    integer :: iunit, io, nline
    character(len=256)    :: record
    character(len=10)     :: dum1, dum2

    iunit = get_unit()
    open(iunit, file='column_table', action='READ', iostat=io)
    if(io/=0) call mpp_error(FATAL, ' find_diagnostic_column: Error in opening column_table')

    num_diag_debug=0
    num_diag_sonde=0
    nline=0
    do while (num_diag_debug < MAX_DIAG_COLUMN .and. num_diag_sonde < MAX_DIAG_COLUMN .and. nline < MAX_DIAG_COLUMN*4)
       nline = nline + 1
       read(iunit,'(a)',end=100) record
       if (record(1:1) == '#') cycle
       if (record(1:10) == '          ') cycle

       !Debug record with index point (index point not supported for sonde output)
       !if (is_master()) print*, index(lowercase(record), "debug"), index(lowercase(record), "index"), trim(record)
       if (index(lowercase(record), "debug") .ne. 0 .and. index(lowercase(record), "index") .ne. 0) then
          if (num_diag_debug >= MAX_DIAG_COLUMN) continue
          num_diag_debug = num_diag_debug + 1
          read(record,*,iostat=io) dum1, dum2, diag_debug_names(num_diag_debug), diag_debug_i(num_diag_debug), diag_debug_j(num_diag_debug), diag_debug_tile(num_diag_debug)
          if (io/=0) then
             print*, ' read_column_table: error on line ', nline
             call mpp_error(FATAL,'error in column_table format')
          endif
       else !debug or sonde record with specified lat-lon
          if (index(lowercase(record), "debug") .ne. 0 ) then
          if (num_diag_debug >= MAX_DIAG_COLUMN) continue
             num_diag_debug = num_diag_debug + 1
             read(record,*,iostat=io) dum1, dum2, diag_debug_names(num_diag_debug), diag_debug_lon(num_diag_debug), diag_debug_lat(num_diag_debug)
             if (io/=0) then
                print*, ' read_column_table: error on line ', nline
                call mpp_error(FATAL,'error in column_table format')
             endif
          else
          if (num_diag_sonde >= MAX_DIAG_COLUMN) continue
             num_diag_sonde = num_diag_sonde + 1
             read(record,*,iostat=io) dum1, dum2, diag_sonde_names(num_diag_sonde), diag_sonde_lon(num_diag_sonde), diag_sonde_lat(num_diag_sonde)
             if (io/=0) then
                print*, ' read_column_table: error on line ', nline
                call mpp_error(FATAL,'error in column_table format')
             endif
          endif

       endif

    enddo
100 continue

  end subroutine read_column_table

 !Note that output lat-lon are in degrees but input is in radians
  subroutine find_diagnostic_column(diag_class, diag_names, diag_i, diag_j, diag_tile, diag_lat, diag_lon, diag_units, grid, agrid, num_diag, ntiles, bd, tile, npx, npy)

    implicit none
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(IN) :: num_diag, tile, ntiles, npx, npy
    character(*), intent(IN) :: diag_class
    character(diag_name_len), intent(IN) :: diag_names(MAX_DIAG_COLUMN)
    integer, dimension(MAX_DIAG_COLUMN), intent(INOUT) :: diag_i, diag_j, diag_tile
    real, dimension(MAX_DIAG_COLUMN), intent(INOUT) :: diag_lat, diag_lon
    integer, dimension(num_diag), intent(OUT) :: diag_units
    real(kind=R_GRID), intent(IN) :: grid(bd%isd+1:bd%ied+1,bd%jsd+1:bd%jed+1,2)
    real(kind=R_GRID), intent(IN) :: agrid(bd%isd:bd%ied,bd%jsd:bd%jed,2)

    integer :: i,j,m,io
    character(80) :: filename
    real(kind=R_GRID), dimension(2):: pp
    real(kind=R_GRID), dimension(3):: vp_12, vp_23, vp_34, vp_14
    real :: dmin, dist
    integer :: isc, iec, jsc, jec
    integer :: isd, ied, jsd, jed
    logical :: point_found

    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    isc = bd%isc
    iec = bd%iec
    jsc = bd%jsc
    jec = bd%jec


    do m=1,num_diag

       point_found = .false.

       !Index specified
       if (diag_i(m) >= -10 .and. diag_j(m) >= -10) then

          if ((diag_tile(m) < 0 .or. diag_tile(m) > ntiles)) then
             if (ntiles > 1) then
                call mpp_error(FATAL, ' find_diagnostic_column: diag_tile must be specified for '//trim(diag_class)//' point '//trim(diag_names(m))//' since ntiles > 1')
             else
                diag_tile(m) = 1
             endif
          endif

          i=diag_i(m)
          j=diag_j(m)

          if (diag_tile(m) == tile .and. i >= isc .and. i <= iec .and. &
               j >= jsc .and. j <= jec) then
             diag_lon(m) = agrid(i,j,1)*rad2deg
             diag_lat(m) = agrid(i,j,2)*rad2deg
             point_found = .true.
          else
             diag_i(m) = -999
             diag_j(m) = -999
             diag_lon(m) = -999.
             diag_lat(m) = -999.
             diag_tile(m) = -1
             point_found = .false.
          endif

       else ! lat-lon specified: find nearest grid cell center

          !diag_lon and diag_lat are in degrees
          ! great_circle_dist wants radians
          pp = (/ diag_lon(m)/rad2deg, diag_lat(m)/rad2deg /)
          !find nearest grid cell: if it is in the halo skip
          dmin = 9.e20
          diag_i(m) = -999
          diag_j(m) = -999
          do j=jsd,jed
             do i=isd,ied
                !no corners
                if ( i < 1    .and. j < 1 ) cycle
                if ( i >= npx .and. j < 1 ) cycle
                if ( i < 1    .and. j >= npy ) cycle
                if ( i >= npx .and. j >= npy ) cycle
                dist =  great_circle_dist(pp, agrid(i,j,:))
                if (dmin >= dist) then
                   diag_i(m) = i
                   diag_j(m) = j
                   dmin = dist
                endif
             enddo
          enddo
          !print*, 'lat-lon point:', mpp_pe(), dmin, diag_i(m), diag_j(m), isc, iec, jsc, jec

          if ( diag_i(m) < isc .or. diag_i(m) > iec .or. diag_j(m) < jsc .or. diag_j(m) > jec ) then
             diag_i(m) = -999
             diag_j(m) = -999
             diag_lon(m) = -999.
             diag_lat(m) = -999.
             diag_tile(m) = -1
             point_found = .false.
          else
             diag_lon(m) = agrid(diag_i(m), diag_j(m), 1)*rad2deg
             diag_lat(m) = agrid(diag_i(m), diag_j(m), 2)*rad2deg
             diag_tile(m) = tile
             point_found = .true.
          endif

       endif

       if (point_found) then

          !Initialize output file
          diag_units(m) = get_unit()
          write(filename, 202) trim(diag_names(m)), trim(diag_class)
202       format(A, '.', A, '.out')
          open(diag_units(m), file=trim(filename), action='WRITE', position='rewind', iostat=io)
          if(io/=0) call mpp_error(FATAL, ' find_diagnostic_column: Error in opening file '//trim(filename))
          !Print debug message
          write(*,'(A, 1x, A, 1x, 1x, A, 2F8.3, 2I5, I3, I04)') trim(diag_class), 'point: ', diag_names(m), diag_lon(m), diag_lat(m), diag_i(m), diag_j(m), diag_tile(m), mpp_pe()

       endif

    enddo

  end subroutine find_diagnostic_column

  subroutine debug_column(pt, delp, delz, u, v, w, q, npz, ncnst, sphum, nwat, zvir, ptop, hydrostatic, bd, Time)

    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(IN) :: npz, ncnst, sphum, nwat
    real, intent(IN) :: zvir, ptop
    logical, intent(IN) :: hydrostatic
    real, dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz), intent(IN) :: pt, delp, w
    real, dimension(bd%is:, bd%js:,1:), intent(IN) :: delz
    real, dimension(bd%isd:bd%ied,bd%jsd:bd%jed+1,npz), intent(IN) :: u
    real, dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed,npz), intent(IN) :: v
    real, dimension(bd%isd:bd%ied, bd%jsd:bd%jed, npz, ncnst), intent(IN) :: q


    type(time_type), intent(IN) :: Time
    integer :: i,j,k,n,l, unit
    real cond, pres, rdg, preshyd(npz), pehyd(npz+1), presdry, preshyddry(npz), pehyddry(npz+1)
    integer :: yr, mon, dd, hr, mn, days, seconds

    rdg = -rdgas/grav

    do n=1,size(diag_debug_i)

       i=diag_debug_i(n)
       j=diag_debug_j(n)
       unit=diag_debug_units(n)

       !Sanity check
       if (i < bd%is .or. i > bd%ie) cycle
       if (j < bd%js .or. j > bd%je) cycle

!< EXAMPLE FORMAT FOR DIAG OUTPUT HEADER
!<               PRINTING ORD               DIAGNOSTICS
!<
!< time stamp:  2016  August   6   0   7  30
!< DIAGNOSTIC POINT COORDINATES, point #   1
!<
!< longitude =  271.354 latitude  =   42.063
!< on processor #    162 :   processor i =     2 ,   processor j =    30

       write(unit, *) "DEBUG POINT ",  diag_debug_names(n)
       write(unit, *)
       call get_date(Time, yr, mon, dd, hr, mn, seconds)
       write(unit, '(A, I6, A12, 4I4)') " Time: ", yr, month_name(mon), dd, hr, mn, seconds
       write(unit, *)
       write(unit, '(A, F8.3, A, F8.3)') ' longitude = ', diag_debug_lon(n), ' latitude = ', diag_debug_lat(n)
       write(unit, '(A, I8, A, I6, A, I6, A, I3)') ' on processor # ', mpp_pe(), ' :  local i = ', i, ',   local j = ', j, ' tile = ', diag_debug_tile(n)
       write(unit, *)

       write(unit,500) 'k', 'T', 'delp', 'delz',   'u',   'v',   'w', 'sphum', 'cond', 'pres', 'NHprime'!, 'pdry', 'NHpdry'
       write(unit,500) ' ', 'K',   'mb',    'm', 'm/s', 'm/s', 'm/s',  'g/kg', 'g/kg', 'mb',   'mb'!,    !  'mb',   'mb'
500    format(A4, A7, A8, A6, A8, A8, A8, A8, A9, A9, A9)
       if (hydrostatic) then
          call mpp_error(NOTE, 'Hydrostatic debug sounding not yet supported')
       else
          pehyd = ptop
          pehyddry = ptop
          do k=1,npz
             pehyd(k+1) = pehyd(k) + delp(i,j,k)
             preshyd(k) = (pehyd(k+1) - pehyd(k))/log(pehyd(k+1)/pehyd(k))
             !pehyddry(k+1) = pehyddry(k) + delp(i,j,k)*(1.-sum(q(i,j,k,1:nwat)))
             !preshyddry(k) = (pehyddry(k+1) - pehyddry(k))/log(pehyddry(k+1)/pehyddry(k))
          enddo

          !do k=2*npz/3,npz
          do k=max(diag_debug_kbottom-diag_debug_nlevels,1),min(diag_debug_kbottom,npz)
             cond = 0.
             do l=2,nwat
                cond = cond + q(i,j,k,l)
             enddo
             pres = rdg*delp(i,j,k)*(1.-cond)/delz(i,j,k)*pt(i,j,k)*(1.+zvir*q(i,j,k,sphum))
             !presdry = rdg*delp(i,j,k)*(1.-cond-q(i,j,k,sphum))/delz(i,j,k)*pt(i,j,k)
             write(unit,'(I4, F7.2, F8.3, I6, F8.3, F8.3, F8.3, F8.3, F9.5, F9.3, F9.3)') &
                  k, pt(i,j,k), delp(i,j,k)*0.01, -int(delz(i,j,k)), u(i,j,k), v(i,j,k), w(i,j,k), &
                  q(i,j,k,sphum)*1000., cond*1000., pres*1.e-2, (pres-preshyd(k))*1.e-2!, presdry*1.e-2, (presdry-preshyddry(k))*1.e-2
          enddo
       endif

       write(unit, *) '==================================================================='
       write(unit, *)

       call mpp_flush(unit)


    enddo

  end subroutine debug_column

  subroutine debug_column_dyn(pt, delp, delz, u, v, w, q, heat_source, cappa, akap, &
       use_heat_source, npz, ncnst, sphum, nwat, zvir, ptop, hydrostatic, bd, Time, k_step, n_step)

    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(IN) :: npz, ncnst, sphum, nwat, k_step, n_step
    real, intent(IN) :: akap, zvir, ptop
    logical, intent(IN) :: hydrostatic, use_heat_source
    real, dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz), intent(IN) :: pt, delp, w, heat_source
    real, dimension(bd%is:, bd%js:,1:), intent(IN) :: delz
    real, dimension(bd%isd:bd%ied,bd%jsd:bd%jed+1,npz), intent(IN) :: u
    real, dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed,npz), intent(IN) :: v
    real, dimension(bd%isd:bd%ied, bd%jsd:bd%jed, npz, ncnst), intent(IN) :: q
    real, dimension(bd%isd:,bd%jsd:,1:), intent(IN) :: cappa

    !Will need to convert variables from internal dyn_core values into logical external values
    ! esp. pt from theta_v to T

    type(time_type), intent(IN) :: Time
    integer :: i,j,k,n,l, unit
    real cond, pres, rdg, Tv, temp, heats, virt, pk, cv_air
    real preshyd(npz), pehyd(npz+1)
    integer yr, mon, dd, hr, mn, seconds

    rdg = -rdgas/grav
    cv_air = cp_air - rdgas

    do n=1,size(diag_debug_i)

       i=diag_debug_i(n)
       j=diag_debug_j(n)
       unit=diag_debug_units(n)

       if (i < bd%is .or. i > bd%ie) cycle
       if (j < bd%js .or. j > bd%je) cycle

       write(unit, *) "DEBUG POINT ",  diag_debug_names(n)
       write(unit, *)
       call get_date(Time, yr, mon, dd, hr, mn, seconds)
       write(unit, '(A, I6, A12, 4I4)') " Time: ", yr, month_name(mon), dd, hr, mn, seconds
       write(unit,*) 'k_split = ', k_step, ', n_split = ', n_step
       write(unit, *)
       write(unit, '(A, F8.3, A, F8.3)') ' longitude = ', diag_debug_lon(n), ' latitude = ', diag_debug_lat(n)
       write(unit, '(A, I8, A, I6, A, I6)') ' on processor # ', mpp_pe(), ' :  local i = ', i, ',   local j = ', j
       write(unit, *)

       write(unit,500) 'k', 'T', 'delp', 'delz',   'u',   'v',   'w', 'sphum', 'cond', 'pres', 'NHprime', 'heat'
       write(unit,500)  ' ', 'K',   'mb',    'm', 'm/s', 'm/s', 'm/s',  'g/kg', 'g/kg', 'mb', 'mb', 'K'
500    format(A4, A7, A8, A6, A8, A8, A8, A8, A9, A9, A9, A8)
          if (hydrostatic) then
             call mpp_error(NOTE, 'Hydrostatic debug sounding not yet supported')
          else
             pehyd = ptop
             do k=1,npz
                pehyd(k+1) = pehyd(k) + delp(i,j,k)
                preshyd(k) = (pehyd(k+1) - pehyd(k))/log(pehyd(k+1)/pehyd(k))
             enddo
             !do k=2*npz/3,npz
             do k=max(diag_debug_kbottom-diag_debug_nlevels,1),min(diag_debug_kbottom,npz)
                cond = 0.
                do l=2,nwat
                   cond = cond + q(i,j,k,l)
                enddo
                virt = (1.+zvir*q(i,j,k,sphum))
#ifdef MOIST_CAPPA
                pres = exp(1./(1.-cappa(i,j,k))*log(rdg*(delp(i,j,k)-cond)/delz(i,j,k)*pt(i,j,k)) )
                pk = exp(cappa(i,j,k)*log(pres))
#else
                pres = exp(1./(1.-akap)*log(rdg*(delp(i,j,k))/delz(i,j,k)*pt(i,j,k)) )
                pk = exp(akap*log(pres))
#endif
                temp = pt(i,j,k)*pk/virt
                if (use_heat_source) then
                   heats = heat_source(i,j,k) / (cv_air*delp(i,j,k))
                else
                   heats = 0.0
                endif
                write(unit,'(I4, F7.2, F8.3, I6, F8.3, F8.3, F8.3, F8.3, F9.5, F9.3, F9.3, G9.3 )') &
                     k, temp, delp(i,j,k)*0.01, -int(delz(i,j,k)), u(i,j,k), v(i,j,k), w(i,j,k), &
                     q(i,j,k,sphum)*1000., cond*1000., pres*1.e-2, (pres-preshyd(k))*1.e-2, heats
             enddo
          endif

       write(unit, *) '==================================================================='
       write(unit, *)

       call mpp_flush(unit)

    enddo

  end subroutine debug_column_dyn

  subroutine sounding_column( pt, delp, delz, u, v, q, peln, pkz, thetae, phis, &
       npz, ncnst, sphum, nwat, hydrostatic, zvir, ng, bd, Time )

    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(IN) :: npz, ncnst, sphum, nwat, ng
    real,    intent(IN) :: zvir
    logical, intent(IN) :: hydrostatic
    real, dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz), intent(IN) :: pt, delp
    real, dimension(bd%is:, bd%js:, 1:), intent(IN) :: delz
    real, dimension(bd%isd:bd%ied,bd%jsd:bd%jed+1,npz), intent(IN) :: u
    real, dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed,npz), intent(IN) :: v
    real, dimension(bd%isd:bd%ied, bd%jsd:bd%jed, npz, ncnst), intent(IN) :: q
    real, dimension(bd%is:bd%ie,npz+1,bd%js:bd%je), intent(in):: peln
    real, dimension(bd%is:bd%ie,bd%js:bd%je,npz),   intent(in):: pkz, thetae
    real, dimension(bd%isd:bd%ied,bd%jsd:bd%jed),   intent(IN) :: phis
    type(time_type), intent(IN) :: Time

    real :: Tv, pres, hght(npz), dewpt, rh, mixr, tmp, qs(1), wspd, wdir, rpk, theta, thetav

    real, PARAMETER :: rgrav = 1./grav
    real, PARAMETER :: rdg = -rdgas*rgrav
    real, PARAMETER :: sounding_top = 10.e2
    real, PARAMETER :: ms_to_knot = 1.9438445
    real, PARAMETER :: p0 = 1000.e2

    integer :: i, j, k, n, unit
    integer :: yr_v, mo_v, dy_v, hr_v, mn_v, sec_v ! need to get numbers for these

    call get_date(Time, yr_v, mo_v, dy_v, hr_v, mn_v, sec_v)

    do n=1,size(diag_sonde_i)

       i=diag_sonde_i(n)
       j=diag_sonde_j(n)
       unit=diag_sonde_units(n)

       if (i < bd%is .or. i > bd%ie) cycle
       if (j < bd%js .or. j > bd%je) cycle


          write(unit,600)        &
               trim(diag_sonde_names(n)), yr_v, mo_v, dy_v, hr_v, init_str, trim(runname)
600       format(A,'.v', I4, I2.2, I2.2, I2.2, '.i', A10, '.', A, '.dat########################################################')
          write(unit,601) trim(diag_sonde_names(n)), yr_v, mo_v, dy_v, hr_v, init_str(1:8),init_str(9:10)
601       format(3x, A16, ' Valid ', I4, I2.2, I2.2, '.', I2.2, 'Z  Init ', A8, '.', A2, 'Z')
          write(unit,'(5x, A, 2F8.3)') trim(runname), diag_sonde_lon(n), diag_sonde_lat(n)
          write(unit,*)
          write(unit,*)        '-------------------------------------------------------------------------------'
          write(unit,'(11A7)') 'PRES', 'HGHT', "TEMP", "DWPT", "RELH", "MIXR", "DRCT", "SKNT", "THTA", "THTE", "THTV"
          write(unit,'(11A7)') 'hPa', 'm', 'C', 'C', '%', 'g/kg', 'deg', 'knot', 'K', 'K', 'K'
          write(unit,*)        '-------------------------------------------------------------------------------'

          if (hydrostatic) then
             call mpp_error(NOTE, 'Hydrostatic diagnostic sounding not yet supported')
          else
             hght(npz) = phis(i,j)*rgrav - 0.5*delz(i,j,npz)
             do k=npz-1,1,-1
                hght(k) = hght(k+1) - 0.5*(delz(i,j,k)+delz(i,j,k+1))
             enddo

             do k=npz,1,-1

                Tv = pt(i,j,k)*(1.+zvir*q(i,j,k,sphum))
                pres = delp(i,j,k)/delz(i,j,k)*rdg*Tv
                !if (pres < sounding_top) cycle

                call qsmith(1, 1, 1, pt(i,j,k:k),   &
                     (/pres/), q(i,j,k:k,sphum), qs)

                mixr = q(i,j,k,sphum)/(1.-sum(q(i,j,k,1:nwat))) ! convert from sphum to mixing ratio
                rh = q(i,j,k,sphum)/qs(1)
                tmp = ( log(max(rh,1.e-2))/ 17.27  + ( pt(i,j,k) - 273.14 )/ ( -35.84 + pt(i,j,k)) )
                dewpt = 237.3* tmp/ ( 1. - tmp ) ! deg C
                wspd = 0.5*sqrt((u(i,j,k)+u(i,j+1,k))*(u(i,j,k)+u(i,j+1,k)) + (v(i,j,k)+v(i+1,j,k))*(v(i,j,k)+v(i+1,j,k)))*ms_to_knot ! convert to knots
                if (wspd > 0.01) then
                   !https://www.eol.ucar.edu/content/wind-direction-quick-reference
                   wdir = atan2(u(i,j,k)+u(i,j+1,k),v(i,j,k)+v(i+1,j,k)) * rad2deg
                else
                   wdir = 0.
                endif
                rpk = exp(-kappa*log(pres/p0))
                theta = pt(i,j,k)*rpk
                thetav = Tv*rpk

                write(unit,'(F7.1, I7, F7.1, F7.1, I7, F7.2, I7, F7.2, F7.1, F7.1, F7.1)') &
                     pres*1.e-2, int(hght(k)), pt(i,j,k)-TFREEZE, dewpt, int(rh*100.), mixr*1.e3, int(wdir), wspd, theta, thetae(i,j,k), thetav
             enddo
          endif

          call mpp_flush(unit)

    enddo


  end subroutine sounding_column



end module fv_diag_column_mod
