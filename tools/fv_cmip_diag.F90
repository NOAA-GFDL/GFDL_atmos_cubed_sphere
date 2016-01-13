module fv_cmip_diag_mod

use mpp_mod,            only: input_nml_file
use fms_mod,            only: open_namelist_file, check_nml_error, &
                              close_file, stdlog, mpp_pe, mpp_root_pe, &
                              write_version_number, file_exist, &
                              error_mesg, FATAL, WARNING
use fms_io_mod,         only: set_domain, nullify_domain, string
use time_manager_mod,   only: time_type
use mpp_domains_mod,    only: domain2d
use diag_manager_mod,   only: diag_axis_init, register_diag_field, &
                              send_data, get_diag_field_id
use tracer_manager_mod, only: get_tracer_index
use field_manager_mod,  only: MODEL_ATMOS

use fv_arrays_mod,      only: fv_atmos_type
use fv_eta_mod,         only: get_eta_level
use fv_diagnostics_mod, only: get_height_given_pressure, &
                              interpolate_vertical, rh_calc, &
                              get_height_field, prt_maxmin

!----------------------------------------------------------------------

implicit none
private

!----------------------------------------------------------------------

public :: fv_cmip_diag_init, fv_cmip_diag, fv_cmip_diag_end

integer :: sphum
integer :: istep = 0

! vertical pressure grids
!    plevs = Table Amon = standard 17 levels + 6 ext
!    plev8 = Table day
!    plev3 = Table 6hrPlev
!    plev7 = Table cfMon,cfDay

real, dimension(23) :: plevs = &
              (/ 100000., 92500., 85000., 70000., 60000., 50000., &
                  40000., 30000., 25000., 20000., 15000., 10000., &
                   7000.,  5000.,  3000.,  2000.,  1000.,   700., &
                    500.,   300.,   200.,   100.,    50. /)
real, dimension(8) :: plev8 = &
              (/ 100000., 85000., 70000., 50000., &
                  25000., 10000.,  5000.,  1000. /)
real, dimension(3) :: plev3 = &
              (/  85000., 50000., 25000. /)

! this set of levels has not been implemented
! MAY NOT BE NEEDED - ONLY USED FOR PHYSICS VARIABLES (cloud amounts)
! (needs conservative pressure-weighted interp)
real, dimension(7) :: plev7 = &
              (/  90000., 74000., 62000., 50000., 37500., 24500., 9000. /)
real, dimension(2,7) :: plev7_bnds = &
              (/ 100000., 80000., 80000., 68000., 68000., 56000., &
                  56000., 44000., 44000., 31000., 31000., 18000., &
                  18000.,     0. /)

!-----------------------------------------------------------------------
!--- namelist ---

logical :: use_extra_levels = .true.  ! use more than the standard
                                      ! 17 pressure levels when possible
integer :: print_freq = 9

namelist /fv_cmip_diag_nml/ use_extra_levels, print_freq

!-----------------------------------------------------------------------

integer, parameter :: MAXPLEVS = 3  ! max plev sets
integer, dimension(0:MAXPLEVS) :: id_ta, id_ua, id_va, id_hus, id_hur, id_wap, id_zg
integer                        :: id_ps
integer, dimension(MAXPLEVS) :: num_pres_levs
real,    dimension(MAXPLEVS,50) :: pressure_levels  ! max 50 levels per set

character(len=16) :: mod_name = 'cmip'
real :: missing_value = -1.e10

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

logical :: module_is_initialized=.false.

CONTAINS

!######################################################################

subroutine fv_cmip_diag_init ( Atm, axes, Time )
type(fv_atmos_type), intent(inout), target :: Atm(:)
integer,             intent(in) :: axes(4)
type(time_type),     intent(in) :: Time

integer :: axes3d(3), area_id, k, id, np, id_plevs, num_std_plevs
integer :: iunit, ierr, io
logical :: used
character(len=8) :: suffix
real :: ptop
character(len=16) :: axis_name, mod_name2

  if (module_is_initialized) then
    call error_mesg ('fv_cmip_diag_mod', &
                     'module has already been initialized', WARNING)
    return
  endif

!----- read namelist -----
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=fv_cmip_diag_nml, iostat=io)
  ierr = check_nml_error (io, 'fv_cmip_diag_nml')
#else
  if (file_exist('input.nml') ) then
    iunit = open_namelist_file()
    ierr=1
    do while (ierr /= 0)
      read (iunit, nml=fv_cmip_diag_nml, iostat=io, end=10)
      ierr = check_nml_error (io, 'fv_cmip_diag_nml')
    enddo
10  call close_file (iunit)
  endif
#endif

!----- write version and namelist to log file -----

  iunit = stdlog()
  call write_version_number ( version, tagname )
  if (mpp_pe() == mpp_root_pe()) write (iunit, nml=fv_cmip_diag_nml)


! Set domain so that diag_manager can access tile information
  call set_domain(Atm(1)%domain)

! axis identifiers
  axes3d(1:2) = axes(1:2)
  area_id = get_diag_field_id ('dynamics', 'area')

! initialize other id's
  id_hur = 0
  id_wap = 0
  id_zg = 0

  sphum = get_tracer_index (MODEL_ATMOS, 'sphum')

!-----------------------------------------------------------------------
! determine the maximum number of standard pressure levels
! first get the pressure (based on ps=1000hPa) at top model level

  if (use_extra_levels) then
    ptop = get_top_press_level (Atm, 1000.e2)
    do k = 23, 1, -1
      if (plevs(k) .gt. ptop) then
        num_std_plevs = k
        exit
      endif
    enddo
  else
    num_std_plevs = 17
  endif

!-----------------------------------------------------------------------
! register variables on model levels

    id_ta(0) = register_diag_field (mod_name, 'ta', axes(1:3), Time, &
                                    'Temperature', 'K',  &
                                    standard_name='air_temperature', &
                                    area=area_id, missing_value=missing_value )

    id_ua(0) = register_diag_field (mod_name, 'ua', axes(1:3), Time, &
                                    'Eastward Wind', 'm s-1',  &
                                    standard_name='eastward_wind', &
                                    area=area_id, missing_value=missing_value )

    id_va(0) = register_diag_field (mod_name, 'va', axes(1:3), Time, &
                                    'Northward Wind', 'm s-1',  &
                                    standard_name='northward_wind', &
                                    area=area_id, missing_value=missing_value )

    id_hus(0) = register_diag_field (mod_name, 'hus', axes(1:3), Time, &
                                     'Specific Humidity', '1',  &
                                     standard_name='specific_humidity', &
                                     area=area_id, missing_value=missing_value )

    id_ps = register_diag_field (mod_name, 'ps', axes(1:2), Time, &
                                 'Surface Air Pressure', 'Pa', &
                                 standard_name='surface_air_pressure', &
                                 area=area_id, missing_value=missing_value )

!-----------------------------------------------------------------------
! loop through all possible pressure level sets
! initialize the pressure axis
! define new 3d grid
! define all 3d state variable on this 3d grid

  do id = 1, MAXPLEVS
    if (id .eq. 1) then
      np = num_std_plevs
      pressure_levels(id,1:np) = plevs(1:np)
      axis_name = 'plevs'
    else if (id .eq. 2) then
      np = size(plev8,1)
      pressure_levels(id,1:np) = plev8
      axis_name = 'plev8'
    else if (id .eq. 3) then
      np = size(plev3,1)
      pressure_levels(id,1:np) = plev3
      axis_name = 'plev3'
    endif

    num_pres_levs(id) = np
    id_plevs = diag_axis_init(axis_name, pressure_levels(id,1:np), &
                              'Pa', 'z', 'Pressure Level', direction=1, set_name="dynamics")
    axes3d(3) = id_plevs
    mod_name2 = mod_name//'_'//trim(axis_name)

    id_ta(id) = register_diag_field (mod_name2, 'ta', axes3d, Time, &
                                     'Temperature', 'K',  &
                                     standard_name='air_temperature', &
                                     area=area_id, missing_value=missing_value )

    id_ua(id) = register_diag_field (mod_name2, 'ua', axes3d, Time, &
                                     'Eastward Wind', 'm s-1',  &
                                     standard_name='eastward_wind', &
                                     area=area_id, missing_value=missing_value )

    id_va(id) = register_diag_field (mod_name2, 'va', axes3d, Time, &
                                     'Northward Wind', 'm s-1',  &
                                     standard_name='northward_wind', &
                                     area=area_id, missing_value=missing_value )

    id_hus(id) = register_diag_field (mod_name2, 'hus', axes3d, Time, &
                                     'Specific Humidity', '1',  &
                                     standard_name='specific_humidity', &
                                     area=area_id, missing_value=missing_value )

    id_hur(id) = register_diag_field (mod_name2, 'hur', axes3d, Time, &
                                     'Relative Humidity', '%',  &
                                     standard_name='relative_humidity', &
                                     area=area_id, missing_value=missing_value )

    id_wap(id) = register_diag_field (mod_name2, 'wap', axes3d, Time, &
                                     'omega (=dp/dt)', 'Pa s-1',  &
                                     standard_name='lagrangian_tendency_of_air_pressure', &
                                     area=area_id, missing_value=missing_value )

    id_zg(id) = register_diag_field (mod_name2, 'zg', axes3d, Time, &
                                     'Geopotential Height', 'm',  &
                                     standard_name='geopotential_height', &
                                     area=area_id, missing_value=missing_value )

  enddo

!--- done ---
  call nullify_domain()
  module_is_initialized=.true.

end subroutine fv_cmip_diag_init

!######################################################################

subroutine fv_cmip_diag ( Atm, zvir, Time )
type(fv_atmos_type), intent(inout) :: Atm(:)
real,                intent(in) :: zvir
type(time_type),     intent(in) :: Time

integer :: isc, iec, jsc, jec, n, i, j, k, id, np
integer :: ngc, npz
logical :: used
logical :: do_print

real, dimension(Atm(1)%bd%isc:Atm(1)%bd%iec, &
                Atm(1)%bd%jsc:Atm(1)%bd%jec) :: pfull

real, dimension(Atm(1)%bd%isc:Atm(1)%bd%iec, &
                Atm(1)%bd%jsc:Atm(1)%bd%jec, &
                Atm(1)%npz) :: rhum

real, dimension(Atm(1)%bd%isc:Atm(1)%bd%iec, &
                Atm(1)%bd%jsc:Atm(1)%bd%jec, &
                Atm(1)%npz+1) :: wz

real, allocatable :: pdat(:,:,:)
real, allocatable :: lnplevs(:)

  if (.not.module_is_initialized) call error_mesg ('fv_cmip_diag_mod', &
                               'module has not been initialized', FATAL)

  n = 1
  isc = Atm(n)%bd%isc; iec = Atm(n)%bd%iec
  jsc = Atm(n)%bd%jsc; jec = Atm(n)%bd%jec
  ngc = Atm(n)%ng
  npz = Atm(n)%npz

  call set_domain(Atm(n)%domain)

  if( print_freq == 0 ) then
    do_print = .false.
  else
    istep = istep + 1
    do_print = mod(istep, print_freq) == 0
  endif

  ! compute relative humidity at model levels (if needed)
  if (count(id_hur(:)>0) > 0) then
    do k=1,npz
      do j=jsc,jec
      do i=isc,iec
        pfull(i,j) = Atm(n)%delp(i,j,k)/(Atm(n)%peln(i,k+1,j)-Atm(n)%peln(i,k,j))
      enddo
      enddo
      call rh_calc (pfull, Atm(n)%pt(isc:iec,jsc:jec,k), &
                    Atm(n)%q(isc:iec,jsc:jec,k,sphum), rhum(isc:iec,jsc:jec,k), do_cmip=.true.)
    enddo
  endif

  ! height field (wz) if needed
  if (count(id_zg(:)>0) > 0) then
    call get_height_field(isc, iec, jsc, jec, ngc, npz, wz, Atm(n)%pt, Atm(n)%q, Atm(n)%peln, zvir)
  endif

  ! process fields on model levels (could add hur,wap)
  if (id_ua(0)  > 0) used = send_data(id_ua(0),  Atm(n)%ua(isc:iec,jsc:jec,:), Time)
  if (id_va(0)  > 0) used = send_data(id_va(0),  Atm(n)%va(isc:iec,jsc:jec,:), Time)
  if (id_ta(0)  > 0) used = send_data(id_ta(0),  Atm(n)%pt(isc:iec,jsc:jec,:), Time)
  if (id_hus(0) > 0) used = send_data(id_hus(0), Atm(n)%q (isc:iec,jsc:jec,:,sphum), Time)
  if (id_ps     > 0) used = send_data(id_ps,     Atm(n)%ps(isc:iec,jsc:jec), Time)

  ! process fields on pressure levels
  do id = 1, MAXPLEVS
    np = num_pres_levs(id)
    allocate(pdat(isc:iec,jsc:jec,np))
    allocate(lnplevs(np))
    lnplevs = log(pressure_levels(id,1:np))

    ! relative humidity
    if (id_hur(id) > 0) then
      call interpolate_vertical_3d (isc, iec, jsc, jec, npz, np, pressure_levels(id,1:np), &
                                    Atm(n)%peln, rhum(isc:iec,jsc:jec,:), pdat(:,:,1:np))
      if (do_print .and. id.eq.1) call prt_maxmin('hur',pdat,isc,iec,jsc,jec,0,np,1.0)
      used = send_data(id_hur(id), pdat(:,:,1:np), Time)
    endif

    ! geopotential height
    if (id_zg(id) > 0) then
      call get_height_given_pressure ( isc, iec, jsc, jec, ngc, npz, wz, &
                                       np, lnplevs, Atm(n)%peln, pdat(:,:,1:np))
      if (do_print .and. id.eq.1) call prt_maxmin('zg',pdat,isc,iec,jsc,jec,0,np,1.0)
      used = send_data(id_zg(id), pdat(:,:,1:np), Time)
    endif

    ! wind components, temperature, spec hum, omega
    if (id_ua(id) > 0) then
      call interpolate_vertical_3d (isc, iec, jsc, jec, npz, np, pressure_levels(id,1:np), Atm(n)%peln, &
                                    Atm(n)%ua(isc:iec,jsc:jec,:), pdat(:,:,1:np))
      if (do_print .and. id.eq.1) call prt_maxmin('ua',pdat,isc,iec,jsc,jec,0,np,1.0)
      used = send_data(id_ua(id), pdat(:,:,1:np), Time)
    endif
    if (id_va(id) > 0) then
      call interpolate_vertical_3d (isc, iec, jsc, jec, npz, np, pressure_levels(id,1:np), Atm(n)%peln, &
                                    Atm(n)%va(isc:iec,jsc:jec,:), pdat(:,:,1:np))
      if (do_print .and. id.eq.1) call prt_maxmin('va',pdat,isc,iec,jsc,jec,0,np,1.0)
      used = send_data(id_va(id), pdat(:,:,1:np), Time)
    endif
    if (id_ta(id) > 0) then
      call interpolate_vertical_3d (isc, iec, jsc, jec, npz, np, pressure_levels(id,1:np), Atm(n)%peln, &
                                    Atm(n)%pt(isc:iec,jsc:jec,:), pdat(:,:,1:np))
      if (do_print .and. id.eq.1) call prt_maxmin('ta',pdat,isc,iec,jsc,jec,0,np,1.0)
      used = send_data(id_ta(id), pdat(:,:,1:np), Time)
    endif
    if (id_hus(id) > 0) then
      call interpolate_vertical_3d (isc, iec, jsc, jec, npz, np, pressure_levels(id,1:np), Atm(n)%peln, &
                                    Atm(n)%q(isc:iec,jsc:jec,:,sphum), pdat(:,:,1:np))
      if (do_print .and. id.eq.1) call prt_maxmin('hus',pdat,isc,iec,jsc,jec,0,np,1.0)
      used = send_data(id_hus(id), pdat(:,:,1:np), Time)
    endif
    if (id_wap(id) > 0) then
      call interpolate_vertical_3d (isc, iec, jsc, jec, npz, np, pressure_levels(id,1:np), Atm(n)%peln, &
                                    Atm(n)%omga(isc:iec,jsc:jec,:), pdat(:,:,1:np))
      if (do_print .and. id.eq.1) call prt_maxmin('wap',pdat,isc,iec,jsc,jec,0,np,1.0)
      used = send_data(id_wap(id), pdat(:,:,1:np), Time)
    endif

    deallocate(pdat)
    deallocate(lnplevs)
  enddo

  call nullify_domain()

!----------------------------------------------------------------------

end subroutine fv_cmip_diag

!######################################################################

subroutine fv_cmip_diag_end

! do nothing, no way to unregister diag fields

end subroutine fv_cmip_diag_end

!######################################################################

subroutine interpolate_vertical_3d (is, ie, js, je, npz, np, plev, peln, a3, ap)
 integer,  intent(in):: is, ie, js, je, npz, np
 real, intent(in):: peln(is:ie,npz+1,js:je)
 real, intent(in):: a3(is:ie,js:je,npz)
 real, intent(in):: plev(np)
 real, intent(out):: ap(is:ie,js:je,np)

 integer :: k

  do k = 1, np
    call interpolate_vertical (is, ie, js, je, npz, plev(k), peln, a3, ap(:,:,k))
  enddo

end subroutine interpolate_vertical_3d

!######################################################################

function get_top_press_level (Atm, psref)
 type(fv_atmos_type), intent(in), target :: Atm(:)
 real,                intent(in)         :: psref

 real :: get_top_press_level
 real, dimension(Atm(1)%npz)   :: pfull
 real, dimension(Atm(1)%npz+1) :: phalf

 ! first get the reference (full & half level) pressure profiles
  call get_eta_level (Atm(1)%npz, psref, pfull, phalf, Atm(1)%ak, Atm(1)%bk)

 ! return pressure at level 1
  get_top_press_level = pfull(1)

end function get_top_press_level

!######################################################################

end module fv_cmip_diag_mod

