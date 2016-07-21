module fv_cmip_diag_mod

use mpp_mod,            only: input_nml_file
use fms_mod,            only: open_namelist_file, check_nml_error, &
                              close_file, stdlog, mpp_pe, mpp_root_pe, &
                              write_version_number, file_exist, &
                              error_mesg, FATAL, WARNING, NOTE
use fms_io_mod,         only: set_domain, nullify_domain, string
use time_manager_mod,   only: time_type
use mpp_domains_mod,    only: domain2d
use diag_manager_mod,   only: register_diag_field, &
                              send_data, get_diag_field_id, &
                              register_static_field, &
                              diag_field_add_attribute, &
                              DIAG_FIELD_NOT_FOUND
use diag_data_mod,      only: CMOR_MISSING_VALUE, null_axis_id
use tracer_manager_mod, only: get_tracer_index
use field_manager_mod,  only: MODEL_ATMOS
use constants_mod,      only: GRAV

use fv_arrays_mod,      only: fv_atmos_type
use fv_diagnostics_mod, only: interpolate_vertical, &
                              get_height_given_pressure, &
                              rh_calc, get_height_field

use atmos_cmip_diag_mod, only: register_cmip_diag_field_2d, &
                               register_cmip_diag_field_3d, &
                               send_cmip_data_3d, cmip_diag_id_type, &
                               query_cmip_diag_id

!----------------------------------------------------------------------

implicit none
private

!----------------------------------------------------------------------

public :: fv_cmip_diag_init, fv_cmip_diag, fv_cmip_diag_end

integer :: sphum

!-----------------------------------------------------------------------
!--- namelist ---

integer :: dummy = 0

namelist /fv_cmip_diag_nml/ dummy

!-----------------------------------------------------------------------

type(cmip_diag_id_type) :: ID_ta, ID_ua, ID_va, ID_hus, ID_hur, ID_wap, ID_zg
integer              :: id_ps, id_orog
integer              :: id_plev200, id_plev850
integer              :: id_ua200, id_va200, id_ua850, id_va850, id_ta850, id_zg500

character(len=5) :: mod_name = 'atmos'

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

logical :: module_is_initialized=.false.

CONTAINS

!######################################################################

subroutine fv_cmip_diag_init ( Atm, axes, Time )
type(fv_atmos_type), intent(inout), target :: Atm(:)
integer,             intent(in) :: axes(4)
type(time_type),     intent(in) :: Time

!-----------------------------------------------------------------------
! local data

integer :: area_id, k, id, np
integer :: n, isc, iec, jsc, jec
integer :: iunit, ierr, io
logical :: used

!-----------------------------------------------------------------------

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
  area_id = get_diag_field_id ('dynamics', 'area')
  if (area_id .eq. DIAG_FIELD_NOT_FOUND) call error_mesg &
        ('fv_cmip_diag_init', 'diagnostic field "dynamics", '// &
         '"area" is not in the diag_table', NOTE)


!-----------------------------------------------------------------------

  sphum = get_tracer_index (MODEL_ATMOS, 'sphum')

!-----------------------------------------------------------------------
! register cmip 3D variables (on model levels and pressure levels)

    ID_ta = register_cmip_diag_field_3d (mod_name, 'ta', Time, &
                      'Air Temperature', 'K', standard_name='air_temperature')

    ID_ua = register_cmip_diag_field_3d (mod_name, 'ua', Time, &
                       'Eastward Wind', 'm s-1', standard_name='eastward_wind')

    ID_va = register_cmip_diag_field_3d (mod_name, 'va', Time, &
                       'Northward Wind', 'm s-1', standard_name='northward_wind')

    ID_hus = register_cmip_diag_field_3d (mod_name, 'hus', Time, &
                       'Specific Humidity', '1.0', standard_name='specific_humidity')

    ID_wap = register_cmip_diag_field_3d (mod_name, 'wap', Time, &
                       'omega (=dp/dt)', 'Pa s-1', standard_name='lagrangian_tendency_of_air_pressure')

    ID_hur = register_cmip_diag_field_3d (mod_name, 'hur', Time, &
                       'Relative Humidity', '%', standard_name='relative_humidity')

    ID_zg = register_cmip_diag_field_3d (mod_name, 'zg', Time, &
                       'Geopotential Height', 'm', standard_name='geopotential_height')


    ! 2D fields

    id_ps = register_diag_field (mod_name, 'ps', axes(1:2), Time, &
                                 'Surface Air Pressure', 'Pa', &
                                 standard_name='surface_air_pressure', &
                                 area=area_id, missing_value=CMOR_MISSING_VALUE )

! surface altitude

  n = 1
  isc = Atm(n)%bd%isc; iec = Atm(n)%bd%iec
  jsc = Atm(n)%bd%jsc; jec = Atm(n)%bd%jec

#ifndef DYNAMICS_ZS
    id_orog = register_static_field (mod_name, 'orog', axes(1:2), &
                                    'Surface Altitude', 'm', &
                                    standard_name='surface_altitude', &
                                    area=area_id)
    if (id_orog > 0) used = send_data (id_orog, Atm(n)%phis(isc:iec,jsc:jec)/GRAV, Time)
#else
!--- for now output this as 'zsurf' from fv_diagnostics ---
!   id_orog = register_diag_field (mod_name, 'orog', axes(1:2), Time, &
!                                  'Surface Altitude', 'm', &
!                                  standard_name='surface_altitude', &
!                                  area=area_id)
#endif

!-----------------------------------------------------------------------
!---- register fields on specific pressure levels ----

  !---- first register pressure levels as scalar variables ----
    id_plev200 = register_static_field (mod_name, 'plev200', (/null_axis_id/), &
                          'Pressure Level', 'Pa', standard_name = 'pressure')
    if(id_plev200 > 0) then
      call diag_field_add_attribute (id_plev200, 'axis', 'Z')
      call diag_field_add_attribute (id_plev200, 'positive', 'up' )
      used = send_data (id_plev200, 200.e2, Time)
    endif

    id_plev850 = register_static_field (mod_name, 'plev850', (/null_axis_id/), &
                          'Pressure Level', 'Pa', standard_name = 'pressure')
    if(id_plev850 > 0) then
      call diag_field_add_attribute (id_plev850, 'axis', 'Z')
      call diag_field_add_attribute (id_plev850, 'positive', 'up' )
      used = send_data (id_plev850, 850.e2, Time)
    endif
  !----

    id_ua200 = register_cmip_diag_field_2d (mod_name, 'ua200', Time, &
                       'Eastward Wind', 'm s-1', standard_name='eastward_wind')
    if (id_ua200 > 0 .and. id_plev200 > 0) &
        call diag_field_add_attribute (id_ua200, 'coordinates', 'plev200')

    id_va200 = register_cmip_diag_field_2d (mod_name, 'va200', Time, &
                       'Northward Wind', 'm s-1', standard_name='northward_wind')
    if (id_va200 > 0 .and. id_plev200 > 0) &
        call diag_field_add_attribute (id_va200, 'coordinates', 'plev200')

  !---- 850 hPa ----
    id_ua850 = register_cmip_diag_field_2d (mod_name, 'ua850', Time, &
                       'Eastward Wind at 850hPa', 'm s-1', standard_name='eastward_wind')
    if (id_ua850 > 0 .and. id_plev850 > 0) &
        call diag_field_add_attribute (id_ua850, 'coordinates', 'plev850')

    id_va850 = register_cmip_diag_field_2d (mod_name, 'va850', Time, &
                       'Northward Wind', 'm s-1', standard_name='northward_wind')
    if (id_va850 > 0 .and. id_plev850 > 0) &
        call diag_field_add_attribute (id_va850, 'coordinates', 'plev850')

    id_ta850 = register_cmip_diag_field_2d (mod_name, 'ta850', Time, &
                       'Air Temperature', 'K', standard_name='air_temperature')
    if (id_ta850 > 0 .and. id_plev850 > 0) &
        call diag_field_add_attribute (id_ta850, 'coordinates', 'plev850')

!   id_zg500 = register_cmip_diag_field_2d (mod_name, 'zg500', Time, &
!                      'Geopotential Height at 500 hPa', 'm', standard_name='geopotential_height')

!--- done ---
  call nullify_domain()
  module_is_initialized=.true.

!-----------------------------------------------------------------------

end subroutine fv_cmip_diag_init

!######################################################################

subroutine fv_cmip_diag ( Atm, zvir, Time )
type(fv_atmos_type), intent(inout) :: Atm(:)
real,                intent(in) :: zvir
type(time_type),     intent(in) :: Time

integer :: isc, iec, jsc, jec, n, i, j, k, id
integer :: ngc, npz
logical :: used

real, dimension(Atm(1)%bd%isc:Atm(1)%bd%iec, &
                Atm(1)%bd%jsc:Atm(1)%bd%jec) :: pfull, dat2

real, dimension(Atm(1)%bd%isc:Atm(1)%bd%iec, &
                Atm(1)%bd%jsc:Atm(1)%bd%jec, &
                Atm(1)%npz) :: rhum

real, dimension(Atm(1)%bd%isc:Atm(1)%bd%iec, &
                Atm(1)%bd%jsc:Atm(1)%bd%jec, &
                Atm(1)%npz+1) :: wz

!-----------------------------------------------------------------------

  if (.not.module_is_initialized) call error_mesg ('fv_cmip_diag_mod', &
                               'module has not been initialized', FATAL)

  n = 1
  isc = Atm(n)%bd%isc; iec = Atm(n)%bd%iec
  jsc = Atm(n)%bd%jsc; jec = Atm(n)%bd%jec
  ngc = Atm(n)%ng
  npz = Atm(n)%npz

  call set_domain(Atm(n)%domain)

  ! compute relative humidity at model levels (if needed)
  if (count(ID_hur%field_id(:)>0) > 0) then
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
  if (count(ID_zg%field_id(:)>0) > 0) then
    call get_height_field(isc, iec, jsc, jec, ngc, npz, wz, Atm(n)%pt, Atm(n)%q, Atm(n)%peln, zvir)
  endif
  
!----------------------------------------------------------------------
! process 2D fields

  if (id_ps     > 0) used = send_data  (id_ps,     Atm(n)%ps  (isc:iec,jsc:jec), Time)

!----------------------------------------------------------------------
  ! process fields on model levels and pressure levels

  if (query_cmip_diag_id(ID_ua)) &
          used = send_cmip_data_3d (ID_ua,  Atm(n)%ua  (isc:iec,jsc:jec,:), Time, phalf=Atm(n)%peln, opt=1)

  if (query_cmip_diag_id(ID_va)) &
          used = send_cmip_data_3d (ID_va,  Atm(n)%va  (isc:iec,jsc:jec,:), Time, phalf=Atm(n)%peln, opt=1)

  if (query_cmip_diag_id(ID_ta)) &
          used = send_cmip_data_3d (ID_ta,  Atm(n)%pt  (isc:iec,jsc:jec,:), Time, phalf=Atm(n)%peln, opt=1)

  if (query_cmip_diag_id(ID_wap)) &
          used = send_cmip_data_3d (ID_wap, Atm(n)%omga(isc:iec,jsc:jec,:), Time, phalf=Atm(n)%peln, opt=1)

  if (query_cmip_diag_id(ID_hus)) &
          used = send_cmip_data_3d (ID_hus, Atm(n)%q   (isc:iec,jsc:jec,:,sphum), Time, phalf=Atm(n)%peln, opt=1)

  ! relative humidity
  if (query_cmip_diag_id(ID_hur)) &
          used = send_cmip_data_3d (ID_hur, rhum(isc:iec,jsc:jec,:), Time, phalf=Atm(n)%peln, opt=1)

  ! geopotential height
  if (query_cmip_diag_id(ID_zg)) &
          used = send_cmip_data_3d (ID_zg, wz, Time, phalf=Atm(n)%peln, opt=1, ext=.true.)

!----------------------------------------------------------------------
! process 2D fields on specific pressure levels
! 
  if (id_ua200 > 0) then
    call interpolate_vertical (isc, iec, jsc, jec, npz, 200.e2, Atm(n)%peln, &
                               Atm(n)%ua(isc:iec,jsc:jec,:), dat2)          
    used = send_data (id_ua200, dat2, Time)
  endif

  if (id_va200 > 0) then
    call interpolate_vertical (isc, iec, jsc, jec, npz, 200.e2, Atm(n)%peln, &
                               Atm(n)%va(isc:iec,jsc:jec,:), dat2)          
    used = send_data (id_va200, dat2, Time)
  endif

  if (id_ua850 > 0) then
    call interpolate_vertical (isc, iec, jsc, jec, npz, 850.e2, Atm(n)%peln, &
                               Atm(n)%ua(isc:iec,jsc:jec,:), dat2)          
    used = send_data (id_ua850, dat2, Time)
  endif

  if (id_va850 > 0) then
    call interpolate_vertical (isc, iec, jsc, jec, npz, 850.e2, Atm(n)%peln, &
                               Atm(n)%va(isc:iec,jsc:jec,:), dat2)          
    used = send_data (id_va850, dat2, Time)
  endif

  if (id_ta850 > 0) then
    call interpolate_vertical (isc, iec, jsc, jec, npz, 850.e2, Atm(n)%peln, &
                               Atm(n)%pt(isc:iec,jsc:jec,:), dat2)          
    used = send_data (id_ta850, dat2, Time)
  endif

 !if (id_zg500 > 0) then
 !  call get_height_given_pressure (isc, iec, jsc, jec, ngc, npz, wz, 1, (/500.e2/), Atm(n)%peln, dat3)
 !  used = send_data (id_ta850, dat2, Time)
 !endif

!----------------------------------------------------------------------

  call nullify_domain()

!----------------------------------------------------------------------

end subroutine fv_cmip_diag

!######################################################################

subroutine fv_cmip_diag_end

  if (.not. module_is_initialized) then
    call error_mesg ('fv_cmip_diag_mod', &
                     'module is not initialized, nothing to deallocate', WARNING)
    return
  endif

! deallocate module arrays
! ID_ta%dealloc; ID_ua%dealloc; ID_ta%dealloc; ID_hus%dealloc
! ID_wap%dealloc; ID_hur%dealloc; ID_zg%dealloc

  module_is_initialized = .false.

end subroutine fv_cmip_diag_end

!######################################################################

end module fv_cmip_diag_mod

