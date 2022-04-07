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
module fv_cmip_diag_mod

use mpp_mod,            only: input_nml_file
use fms_mod,            only: check_nml_error, string, &
                              stdlog, mpp_pe, mpp_root_pe, &
                              write_version_number, &
                              error_mesg, FATAL, WARNING, NOTE
use time_manager_mod,   only: time_type
use mpp_domains_mod,    only: domain2d
use diag_manager_mod,   only: register_diag_field, diag_axis_init, &
                              send_data, get_diag_field_id, &
                              register_static_field, &
                              diag_field_add_attribute, &
                              DIAG_FIELD_NOT_FOUND
use diag_data_mod,      only: CMOR_MISSING_VALUE, null_axis_id
use tracer_manager_mod, only: get_tracer_index
use field_manager_mod,  only: MODEL_ATMOS
use constants_mod,      only: GRAV, RDGAS

use fv_mapz_mod,        only: E_Flux
use fv_arrays_mod,      only: fv_atmos_type
use fv_diagnostics_mod, only: interpolate_vertical, &
                              get_height_given_pressure, &
                              rh_calc, get_height_field, get_vorticity

use atmos_cmip_diag_mod, only: register_cmip_diag_field_2d, &
                               register_cmip_diag_field_3d, &
                               send_cmip_data_3d, cmip_diag_id_type, &
                               query_cmip_diag_id

!----------------------------------------------------------------------

implicit none
private

!----------------------------------------------------------------------

public :: fv_cmip_diag_init, fv_cmip_diag, fv_cmip_diag_end

integer :: sphum, nql, nqi, nqa

!-----------------------------------------------------------------------
!--- namelist ---

integer :: dummy = 0

namelist /fv_cmip_diag_nml/ dummy

!-----------------------------------------------------------------------

type(cmip_diag_id_type) :: ID_ta, ID_ua, ID_va, ID_hus, ID_hur, ID_wap, ID_zg, &
                           ID_u2, ID_v2, ID_t2, ID_wap2, ID_uv, ID_ut, ID_vt,  &
                           ID_uwap, ID_vwap, ID_twap, ID_wa, &
                           ID_cls, ID_clws, ID_clis

integer              :: id_ps, id_orog
integer              :: id_ua200, id_va200, id_ua850, id_va850, &
                        id_ta500, id_ta700, id_ta850, id_zg500, &
                        id_zg100, id_zg10,  id_zg1000,          &
                        id_hus850, id_wap500, id_ua10
integer              :: id_rv200, id_rv500, id_rv850, id_vortmean

character(len=5) :: mod_name = 'atmos'

! version number of this module
! Include variable "version" to be written to log file.
#include<file_version.h>

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

! data for single pressure levels
integer, dimension(7) :: plevels = (/10,100,200,500,700,850,1000/)
integer, dimension(7) :: id_plevels
integer, parameter    :: id_p10=1, id_p100=2, id_p200=3, id_p500=4, id_p700=5, id_p850=6, id_p1000=7
character(len=4)      :: plabel
integer               :: id_pl700, id_pl700_bnds, id_nv
!-----------------------------------------------------------------------

  if (module_is_initialized) then
    call error_mesg ('fv_cmip_diag_mod', &
                     'module has already been initialized', WARNING)
    return
  endif

!----- read namelist -----
  read (input_nml_file, nml=fv_cmip_diag_nml, iostat=io)
  ierr = check_nml_error (io, 'fv_cmip_diag_nml')

!----- write version and namelist to log file -----

  iunit = stdlog()
  call write_version_number ( 'FV_CMIP_DIAG_MOD', version )
  if (mpp_pe() == mpp_root_pe()) write (iunit, nml=fv_cmip_diag_nml)

! axis identifiers
  area_id = get_diag_field_id ('dynamics', 'area')
  if (area_id .eq. DIAG_FIELD_NOT_FOUND) call error_mesg &
        ('fv_cmip_diag_init', 'diagnostic field "dynamics", '// &
         '"area" is not in the diag_table', NOTE)


!-----------------------------------------------------------------------

  sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
  nql   = get_tracer_index (MODEL_ATMOS, 'liq_wat')
  nqi   = get_tracer_index (MODEL_ATMOS, 'ice_wat')
  nqa   = get_tracer_index (MODEL_ATMOS, 'cld_amt')

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

    ID_wa = register_cmip_diag_field_3d (mod_name, 'wa', Time, &
                       'Upward Air Velocity', 'm s-1', standard_name='upward_air_velocity')

    ID_zg = register_cmip_diag_field_3d (mod_name, 'zg', Time, &
                       'Geopotential Height', 'm', standard_name='geopotential_height', axis='half')

!-----------------------------------------------------------------------
! register products of 3D variables (on model levels and pressure levels)

    ID_u2 = register_cmip_diag_field_3d (mod_name, 'u2', Time, &
                      'Square of Eastward Wind', 'm2 s-2', standard_name='square_of_eastward_wind')

    ID_v2 = register_cmip_diag_field_3d (mod_name, 'v2', Time, &
                       'Square of Northward Wind', 'm2 s-2', standard_name='square_of_northward_wind')

    ID_t2 = register_cmip_diag_field_3d (mod_name, 't2', Time, &
                       'Square of Air Temperature', 'K2', standard_name='square_of_air_temperature')

    ID_wap2 = register_cmip_diag_field_3d (mod_name, 'wap2', Time, &
                       'Square of Omega', 'Pa2 s-2', standard_name='square_of_omega')

    ID_uv = register_cmip_diag_field_3d (mod_name, 'uv', Time, &
                        'Eastward Wind times Northward Wind', 'm2 s-2', &
                         standard_name='product_of_eastward_wind_and_northward_wind')

    ID_ut = register_cmip_diag_field_3d (mod_name, 'ut', Time, &
                        'Air Temperature times Eastward Wind', 'K m s-1', &
                         standard_name='product_of_eastward_wind_and_air_temperature')

    ID_vt = register_cmip_diag_field_3d (mod_name, 'vt', Time, &
                        'Air Temperature times Northward Wind', 'K m s-1', &
                         standard_name='product_of_northward_wind_and_air_temperature')

    ID_uwap = register_cmip_diag_field_3d (mod_name, 'uwap', Time, &
                           'Eastward Wind times Omega', 'K m s-1', &
                           standard_name='product_of_eastward_wind_and_omega')

    ID_vwap = register_cmip_diag_field_3d (mod_name, 'vwap', Time, &
                           'Northward Wind times Omega', 'K m s-1', &
                           standard_name='product_of_northward_wind_and_omega')

    ID_twap = register_cmip_diag_field_3d (mod_name, 'twap', Time, &
                           'Air Temperature times Omega', 'K m s-1', &
                           standard_name='product_of_omega_and_air_temperature')

!-----------------------------------------------------------------------
! stratiform cloud tracers

    if (nql > 0) then
      ID_clws = register_cmip_diag_field_3d (mod_name, 'clws', Time, &
                        'Mass Fraction of Stratiform Cloud Liquid Water', '1.0', &
                 standard_name='mass_fraction_of_stratiform_cloud_liquid_water_in_air')
    endif
    if (nqi > 0) then
      ID_clis = register_cmip_diag_field_3d (mod_name, 'clis', Time, &
                        'Mass Fraction of Stratiform Cloud Ice', '1.0', &
                 standard_name='mass_fraction_of_convective_cloud_ice_in_air')
    endif
    if (nqa > 0) then
      ID_cls = register_cmip_diag_field_3d (mod_name, 'cls', Time, &
                        'Percentage Cover of Stratiform Cloud', '%', &
                 standard_name='stratiform_cloud_area_fraction_in_atmosphere_layer')
    endif
!-----------------------------------------------------------------------
    ! 2D fields

    id_ps = register_cmip_diag_field_2d (mod_name, 'ps', Time, &
                                 'Surface Air Pressure', 'Pa', &
                                 standard_name='surface_air_pressure' )

! surface altitude

  n = 1
  isc = Atm(n)%bd%isc; iec = Atm(n)%bd%iec
  jsc = Atm(n)%bd%jsc; jec = Atm(n)%bd%jec

#ifndef DYNAMICS_ZS
    id_orog = register_static_field (mod_name, 'orog', axes(1:2), &
                                    'Surface Altitude', 'm', &
                                    standard_name='surface_altitude', &
                                    area=area_id, interp_method='conserve_order1')
    if (id_orog > 0) used = send_data (id_orog, Atm(n)%phis(isc:iec,jsc:jec)/GRAV, Time)
#else
!--- for now output this as 'zsurf' from fv_diagnostics ---
!   id_orog = register_diag_field (mod_name, 'orog', axes(1:2), Time, &
!                                  'Surface Altitude', 'm', &
!                                  standard_name='surface_altitude', &
!                                  area=area_id, interp_method='conserve_order1')
#endif

!-----------------------------------------------------------------------
!---- register fields on specific pressure levels ----
!-----------------------------------------------------------------------

  !---- first register pressure levels as scalar variables ----
    do k = 1, size(plevels,1)
      write(plabel,'(i4)') plevels(k)
      plabel = adjustl(plabel)
      id_plevels(k) = register_static_field (mod_name, 'p'//trim(plabel), (/null_axis_id/), &
                             trim(plabel)//' hPa', 'Pa', standard_name='air_pressure')
      if (id_plevels(k) > 0) then
        call diag_field_add_attribute (id_plevels(k), 'axis', 'Z')
        call diag_field_add_attribute (id_plevels(k), 'positive', 'down' )
        used = send_data (id_plevels(k), real(plevels(k))*100., Time)
      endif
    enddo

    id_pl700 = register_static_field (mod_name, 'pl700', (/null_axis_id/), &
                        '700 hPa Average', 'Pa', standard_name='air_pressure')
    if (id_pl700 > 0) then
      call diag_field_add_attribute (id_pl700, 'axis', 'Z')
      call diag_field_add_attribute (id_pl700, 'positive', 'down' )
      call diag_field_add_attribute (id_pl700, 'comment', 'average at levels 600,700,850 hPa' )
      ! add bounds
      id_nv = diag_axis_init('nv', (/1.,2./), 'none', 'N', 'vertex number', set_name='nv')
      id_pl700_bnds = register_static_field (mod_name, 'pl700_bnds', (/id_nv,null_axis_id/), &
                                     '700 hPa boundaries', 'Pa', standard_name='air_pressure')
      if (id_pl700_bnds > 0) then
        call diag_field_add_attribute (id_pl700, 'bounds', 'pl700_bnds' )
        used = send_data (id_pl700_bnds, (/850.e2,600.e2/), Time)
      endif
      used = send_data (id_pl700, 700.e2, Time)
    endif


  !---- register field on single pressure levels ----

    id_ua10 = register_cmip_diag_field_2d (mod_name, 'ua10', Time, &
                       'Eastward Wind at 10hPa', 'm s-1', standard_name='eastward_wind')
    if (id_ua10 > 0 .and. id_plevels(id_p10) > 0) &
        call diag_field_add_attribute (id_ua10, 'coordinates', 'p10')

    id_ua200 = register_cmip_diag_field_2d (mod_name, 'ua200', Time, &
                       'Eastward Wind', 'm s-1', standard_name='eastward_wind')
    if (id_ua200 > 0 .and. id_plevels(id_p200) > 0) &
        call diag_field_add_attribute (id_ua200, 'coordinates', 'p200')

    id_va200 = register_cmip_diag_field_2d (mod_name, 'va200', Time, &
                       'Northward Wind', 'm s-1', standard_name='northward_wind')
    if (id_va200 > 0 .and. id_plevels(id_p200) > 0) &
        call diag_field_add_attribute (id_va200, 'coordinates', 'p200')

  !---- wind components at 850 hPa ----
    id_ua850 = register_cmip_diag_field_2d (mod_name, 'ua850', Time, &
                       'Eastward Wind at 850hPa', 'm s-1', standard_name='eastward_wind')
    if (id_ua850 > 0 .and. id_plevels(id_p850) > 0) &
        call diag_field_add_attribute (id_ua850, 'coordinates', 'p850')

    id_va850 = register_cmip_diag_field_2d (mod_name, 'va850', Time, &
                       'Northward Wind', 'm s-1', standard_name='northward_wind')
    if (id_va850 > 0 .and. id_plevels(id_p850) > 0) &
        call diag_field_add_attribute (id_va850, 'coordinates', 'p850')

  !---- temperature ----

    id_ta500 = register_cmip_diag_field_2d (mod_name, 'ta500', Time, &
                       'Air Temperature', 'K', standard_name='air_temperature')
    if (id_ta500 > 0 .and. id_plevels(id_p500) > 0) &
        call diag_field_add_attribute (id_ta500, 'coordinates', 'p500')

    id_ta700 = register_cmip_diag_field_2d (mod_name, 'ta700', Time, &
                       'Air Temperature', 'K', standard_name='air_temperature')
    if (id_ta700 > 0 .and. id_plevels(id_p700) > 0) &
        call diag_field_add_attribute (id_ta700, 'coordinates', 'p700')

    id_ta850 = register_cmip_diag_field_2d (mod_name, 'ta850', Time, &
                       'Air Temperature', 'K', standard_name='air_temperature')
    if (id_ta850 > 0 .and. id_plevels(id_p850) > 0) &
        call diag_field_add_attribute (id_ta850, 'coordinates', 'p850')

  !---- spec humidity at 850 hPa ----

    id_hus850 = register_cmip_diag_field_2d (mod_name, 'hus850', Time, &
                       'Specific Humidity', '1.0', standard_name='specific_humidity')
    if (id_hus850 > 0 .and. id_plevels(id_p850) > 0) &
        call diag_field_add_attribute (id_hus850, 'coordinates', 'p850')

  !---- relative vorticity at 200, 500, 850 hPa ----
    id_rv200 = register_cmip_diag_field_2d (mod_name, 'rv200', Time, &
                  'Relative Vorticity at 200 hPa', 's-1', standard_name='atmosphere_relative_vorticity')
    if (id_rv200 > 0 .and. id_plevels(id_p200) > 0) &
        call diag_field_add_attribute (id_rv200, 'coordinates', 'p200')

    id_rv500 = register_cmip_diag_field_2d (mod_name, 'rv500', Time, &
                  'Relative Vorticity at 500 hPa', 's-1', standard_name='atmosphere_relative_vorticity')
    if (id_rv500 > 0 .and. id_plevels(id_p500) > 0) &
        call diag_field_add_attribute (id_rv500, 'coordinates', 'p500')

    id_rv850 = register_cmip_diag_field_2d (mod_name, 'rv850', Time, &
                  'Relative Vorticity at 850 hPa', 's-1', standard_name='atmosphere_relative_vorticity')
    if (id_rv850 > 0 .and. id_plevels(id_p850) > 0) &
        call diag_field_add_attribute (id_rv850, 'coordinates', 'p850')

  !---- mean relative vorticity 600, 700, 850 hPa ----

    id_vortmean = register_cmip_diag_field_2d (mod_name, 'vortmean', Time, &
                 'Mean Relative Vorticity over 600-850 hPa', 's-1',  &
                        standard_name='atmosphere_relative_vorticity')
    if (id_vortmean > 0 .and. id_pl700 > 0) &
        call diag_field_add_attribute (id_vortmean, 'coordinates', 'pl700')

  !---- omega at 500 hPa ----

    id_wap500 = register_cmip_diag_field_2d (mod_name, 'wap500', Time, &
                       'omega (=dp/dt)', 'Pa s-1', standard_name='lagrangian_tendency_of_air_pressure')
    if (id_wap500 > 0 .and. id_plevels(id_p500) > 0) &
        call diag_field_add_attribute (id_wap500, 'coordinates', 'p500')

  !---- geopotential height ----

    id_zg10 = register_cmip_diag_field_2d (mod_name, 'zg10', Time, &
                       'Geopotential Height at 10 hPa', 'm', standard_name='geopotential_height')
    if (id_zg10 > 0 .and. id_plevels(id_p10) > 0) &
        call diag_field_add_attribute (id_zg10, 'coordinates', 'p10')

    id_zg100 = register_cmip_diag_field_2d (mod_name, 'zg100', Time, &
                       'Geopotential Height at 100 hPa', 'm', standard_name='geopotential_height')
    if (id_zg100 > 0 .and. id_plevels(id_p1000) > 0) &
        call diag_field_add_attribute (id_zg100, 'coordinates', 'p100')

    id_zg500 = register_cmip_diag_field_2d (mod_name, 'zg500', Time, &
                       'Geopotential Height at 500 hPa', 'm', standard_name='geopotential_height')
    if (id_zg500 > 0 .and. id_plevels(id_p500) > 0) &
        call diag_field_add_attribute (id_zg500, 'coordinates', 'p500')

    id_zg1000 = register_cmip_diag_field_2d (mod_name, 'zg1000', Time, &
                       'Geopotential Height at 1000 hPa', 'm', standard_name='geopotential_height')
    if (id_zg1000 > 0 .and. id_plevels(id_p1000) > 0) &
        call diag_field_add_attribute (id_zg1000, 'coordinates', 'p1000')

!--- done ---
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
logical :: compute_wa , compute_rh

real, dimension(Atm(1)%bd%isc:Atm(1)%bd%iec, &
                Atm(1)%bd%jsc:Atm(1)%bd%jec) :: pfull, dat2, &
                                                rv850, rv700, rv600

real, dimension(Atm(1)%bd%isc:Atm(1)%bd%iec, &
                Atm(1)%bd%jsc:Atm(1)%bd%jec,1) :: dat3

real, dimension(Atm(1)%bd%isc:Atm(1)%bd%iec, &
                Atm(1)%bd%jsc:Atm(1)%bd%jec, &
                Atm(1)%npz) :: rhum, wa, rv

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

  ! set flags for computing quantities
  compute_rh = .false.
  compute_wa = .false.
  if (count(ID_hur%field_id(:)>0) > 0) compute_rh = .true.
  if (count(ID_wa%field_id(:)>0)  > 0) compute_wa = .true.

  ! compute relative humidity at model levels (if needed)
  if (compute_rh .or. compute_wa) then
    do k=1,npz
      do j=jsc,jec
      do i=isc,iec
        pfull(i,j) = Atm(n)%delp(i,j,k)/(Atm(n)%peln(i,k+1,j)-Atm(n)%peln(i,k,j))
      enddo
      enddo
      ! compute relative humidity
      if (compute_rh) then
        call rh_calc (pfull, Atm(n)%pt(isc:iec,jsc:jec,k), &
                    Atm(n)%q(isc:iec,jsc:jec,k,sphum), rhum(isc:iec,jsc:jec,k), do_cmip=.true.)
      endif
      ! compute vertical velocity
      if (compute_wa) then
        wa(isc:iec,jsc:jec,k) = -(Atm(n)%omga(isc:iec,jsc:jec,k)*Atm(n)%pt(isc:iec,jsc:jec,k)/ &
                                        pfull(isc:iec,jsc:jec))*(RDGAS/GRAV)
      endif
    enddo
  endif

  ! height field (wz) if needed
  if (count(ID_zg%field_id(:)>0) > 0 .or. any((/id_zg10,id_zg100,id_zg500,id_zg1000/) > 0)) then
    call get_height_field(isc, iec, jsc, jec, ngc, npz, Atm(n)%flagstruct%hydrostatic, Atm(n)%delz, &
                          wz, Atm(n)%pt, Atm(n)%q, Atm(n)%peln, zvir)
  endif

  ! relative vorticity
  if (any((/id_rv200,id_rv500,id_rv850,id_vortmean/) > 0)) then
    call get_vorticity(isc, iec, jsc, jec, Atm(n)%bd%isd, Atm(n)%bd%ied, Atm(n)%bd%jsd, Atm(n)%bd%jed, npz, &
                       Atm(n)%u, Atm(n)%v, rv, Atm(n)%gridstruct%dx, Atm(n)%gridstruct%dy, Atm(n)%gridstruct%rarea)
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

    ! vertical velocity
  if (query_cmip_diag_id(ID_wa)) &
          used = send_cmip_data_3d (ID_wa, wa(isc:iec,jsc:jec,:), Time, phalf=Atm(n)%peln, opt=1)

    ! geopotential height
  if (query_cmip_diag_id(ID_zg)) &
          used = send_cmip_data_3d (ID_zg, wz, Time, phalf=Atm(n)%peln, opt=1, ext=.true.)

!----------------------------------------------------------------------
  ! process product of fields on model levels and/or pressure levels

  if (query_cmip_diag_id(ID_u2)) &
          used = send_cmip_data_3d (ID_u2,  Atm(n)%ua  (isc:iec,jsc:jec,:)*Atm(n)%ua  (isc:iec,jsc:jec,:), &
                                    Time, phalf=Atm(n)%peln, opt=1)

  if (query_cmip_diag_id(ID_v2)) &
          used = send_cmip_data_3d (ID_v2,  Atm(n)%va  (isc:iec,jsc:jec,:)*Atm(n)%va  (isc:iec,jsc:jec,:), &
                                    Time, phalf=Atm(n)%peln, opt=1)

  if (query_cmip_diag_id(ID_t2)) &
          used = send_cmip_data_3d (ID_t2,  Atm(n)%pt  (isc:iec,jsc:jec,:)*Atm(n)%pt  (isc:iec,jsc:jec,:), &
                                    Time, phalf=Atm(n)%peln, opt=1)

  if (query_cmip_diag_id(ID_wap2)) &
          used = send_cmip_data_3d (ID_wap2, Atm(n)%omga(isc:iec,jsc:jec,:)*Atm(n)%omga(isc:iec,jsc:jec,:), &
                                    Time, phalf=Atm(n)%peln, opt=1)

  if (query_cmip_diag_id(ID_uv)) &
          used = send_cmip_data_3d (ID_uv,  Atm(n)%ua  (isc:iec,jsc:jec,:)*Atm(n)%va  (isc:iec,jsc:jec,:), &
                                    Time, phalf=Atm(n)%peln, opt=1)

  if (query_cmip_diag_id(ID_ut)) &
          used = send_cmip_data_3d (ID_ut,  Atm(n)%ua  (isc:iec,jsc:jec,:)*Atm(n)%pt  (isc:iec,jsc:jec,:), &
                                    Time, phalf=Atm(n)%peln, opt=1)

  if (query_cmip_diag_id(ID_vt)) &
          used = send_cmip_data_3d (ID_vt,  Atm(n)%va  (isc:iec,jsc:jec,:)*Atm(n)%pt  (isc:iec,jsc:jec,:), &
                                    Time, phalf=Atm(n)%peln, opt=1)

  if (query_cmip_diag_id(ID_uwap)) &
          used = send_cmip_data_3d (ID_uwap, Atm(n)%ua  (isc:iec,jsc:jec,:)*Atm(n)%omga(isc:iec,jsc:jec,:), &
                                    Time, phalf=Atm(n)%peln, opt=1)

  if (query_cmip_diag_id(ID_vwap)) &
          used = send_cmip_data_3d (ID_vwap, Atm(n)%va  (isc:iec,jsc:jec,:)*Atm(n)%omga(isc:iec,jsc:jec,:), &
                                    Time, phalf=Atm(n)%peln, opt=1)

  if (query_cmip_diag_id(ID_twap)) &
          used = send_cmip_data_3d (ID_twap, Atm(n)%pt  (isc:iec,jsc:jec,:)*Atm(n)%omga(isc:iec,jsc:jec,:), &
                                    Time, phalf=Atm(n)%peln, opt=1)

!----------------------------------------------------------------------
! stratiform cloud tracers (only on model levels)

  if (query_cmip_diag_id(ID_cls))  used = send_cmip_data_3d (ID_cls,  Atm(n)%q(isc:iec,jsc:jec,:,nqa)*100., Time)
  if (query_cmip_diag_id(ID_clws)) used = send_cmip_data_3d (ID_clws, Atm(n)%q(isc:iec,jsc:jec,:,nql), Time)
  if (query_cmip_diag_id(ID_clis)) used = send_cmip_data_3d (ID_clis, Atm(n)%q(isc:iec,jsc:jec,:,nqi), Time)

!----------------------------------------------------------------------
! process 2D fields on specific pressure levels
!
  if (id_ua10 > 0) then
    call interpolate_vertical (isc, iec, jsc, jec, npz, 10.e2, Atm(n)%peln, &
                               Atm(n)%ua(isc:iec,jsc:jec,:), dat2)
    used = send_data (id_ua10, dat2, Time)
  endif

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

  if (id_ta500 > 0) then
    call interpolate_vertical (isc, iec, jsc, jec, npz, 500.e2, Atm(n)%peln, &
                               Atm(n)%pt(isc:iec,jsc:jec,:), dat2)
    used = send_data (id_ta500, dat2, Time)
  endif

  if (id_ta700 > 0) then
    call interpolate_vertical (isc, iec, jsc, jec, npz, 700.e2, Atm(n)%peln, &
                               Atm(n)%pt(isc:iec,jsc:jec,:), dat2)
    used = send_data (id_ta700, dat2, Time)
  endif

  if (id_ta850 > 0) then
    call interpolate_vertical (isc, iec, jsc, jec, npz, 850.e2, Atm(n)%peln, &
                               Atm(n)%pt(isc:iec,jsc:jec,:), dat2)
    used = send_data (id_ta850, dat2, Time)
    endif

  if (id_hus850 > 0) then
    call interpolate_vertical (isc, iec, jsc, jec, npz, 850.e2, Atm(n)%peln, &
                               Atm(n)%q(isc:iec,jsc:jec,:,sphum), dat2)
    used = send_data (id_hus850, dat2, Time)
  endif

  if (id_wap500 > 0) then
    call interpolate_vertical (isc, iec, jsc, jec, npz, 500.e2, Atm(n)%peln, &
                               Atm(n)%omga(isc:iec,jsc:jec,:), dat2)
    used = send_data (id_wap500, dat2, Time)
  endif

  if (id_rv200 > 0) then
    call interpolate_vertical (isc, iec, jsc, jec, npz, 200.e2, Atm(n)%peln, rv, dat2)
    used = send_data (id_rv200, dat2, Time)
  endif

  if (id_rv500 > 0) then
    call interpolate_vertical (isc, iec, jsc, jec, npz, 500.e2, Atm(n)%peln, rv, dat2)
    used = send_data (id_rv500, dat2, Time)
  endif

  if (id_rv850 > 0 .or. id_vortmean > 0 ) then
    call interpolate_vertical (isc, iec, jsc, jec, npz, 850.e2, Atm(n)%peln, rv, rv850)
    if (id_rv850 > 0) used = send_data (id_rv850, rv850, Time)
    if (id_vortmean > 0) then
      call interpolate_vertical (isc, iec, jsc, jec, npz, 600.e2, Atm(n)%peln, rv, rv600)
      call interpolate_vertical (isc, iec, jsc, jec, npz, 700.e2, Atm(n)%peln, rv, rv700)
      used = send_data (id_vortmean, (rv600+rv700+rv850)/3., Time)
    endif
  endif

  if (id_zg10 > 0) then
    call get_height_given_pressure (isc, iec, jsc, jec, npz, wz, 1, (/id_zg10/), &
                                    (/log(10.e2)/), Atm(n)%peln, dat3)
    used = send_data (id_zg10, dat3(:,:,1), Time)
  endif

  if (id_zg100 > 0) then
    call get_height_given_pressure (isc, iec, jsc, jec, npz, wz, 1, (/id_zg100/), &
                                    (/log(100.e2)/), Atm(n)%peln, dat3)
    used = send_data (id_zg100, dat3(:,:,1), Time)
  endif

  if (id_zg500 > 0) then
    call get_height_given_pressure (isc, iec, jsc, jec, npz, wz, 1, (/id_zg500/), &
                                    (/log(500.e2)/), Atm(n)%peln, dat3)
    used = send_data (id_zg500, dat3(:,:,1), Time)
  endif

  if (id_zg1000 > 0) then
    call get_height_given_pressure (isc, iec, jsc, jec, npz, wz, 1, (/id_zg1000/), &
                                    (/log(1000.e2)/), Atm(n)%peln, dat3)
    used = send_data (id_zg1000, dat3(:,:,1), Time)
  endif

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

end module fv_cmip_diag_mod

