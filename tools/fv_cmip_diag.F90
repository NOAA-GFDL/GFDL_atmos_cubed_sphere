module fv_cmip_diag_mod

use mpp_mod,            only: input_nml_file
use fms_mod,            only: open_namelist_file, check_nml_error, &
                              close_file, stdlog, mpp_pe, mpp_root_pe, &
                              write_version_number, file_exist, &
                              error_mesg, FATAL, WARNING, NOTE
use fms_io_mod,         only: set_domain, nullify_domain, string
use time_manager_mod,   only: time_type
use mpp_domains_mod,    only: domain2d
use diag_manager_mod,   only: diag_axis_init, register_diag_field, &
                              send_data, get_diag_field_id, &
                              register_static_field, &
                              diag_axis_add_attribute, &
                              diag_field_add_attribute, &
                              DIAG_FIELD_NOT_FOUND
use diag_data_mod,      only: CMOR_MISSING_VALUE
use tracer_manager_mod, only: get_tracer_index
use field_manager_mod,  only: MODEL_ATMOS
use constants_mod,      only: GRAV

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
!    plev  = Table Amon = standard 17 levels + 6 ext
!    plev8 = Table day
!    plev3 = Table 6hrPlev
!    plev7 = Table cfMon,cfDay

real, dimension(23) :: plev = &
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
integer                        :: id_ps, id_orog
integer                        :: id_lev, id_nv, id_ap, id_b, &
                                  id_ap_bnds, id_b_bnds, id_lev_bnds
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

!-----------------------------------------------------------------------
! local data

integer :: axes3d(3), area_id, k, id, np, id_plev, num_std_plevs
integer :: n, isc, iec, jsc, jec
integer :: iunit, ierr, io
logical :: used
character(len=8) :: suffix
real :: ptop
character(len=16) :: axis_name, mod_name2

real                          :: p0
real, dimension(Atm(1)%npz)   :: ap, b, lev
real, dimension(2,Atm(1)%npz) :: ap_bnds, b_bnds, lev_bnds

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
  axes3d(1:2) = axes(1:2)
  area_id = get_diag_field_id ('dynamics', 'area')
  if (area_id .eq. DIAG_FIELD_NOT_FOUND) call error_mesg &
        ('fv_cmip_diag_init', 'diagnostic field "dynamics", '// &
         '"area" is not in the diag_table', NOTE)

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
      if (plev(k) .gt. ptop) then
        num_std_plevs = k
        exit
      endif
    enddo
  else
    num_std_plevs = 17
  endif

!-----------------------------------------------------------------------
! vertical coordinate variables

    p0 = 100000.
    do k = 1, Atm(1)%npz
      ap_bnds(1,k) = Atm(1)%ak(k)
      ap_bnds(2,k) = Atm(1)%ak(k+1)
      b_bnds(1,k) = Atm(1)%bk(k)
      b_bnds(2,k) = Atm(1)%bk(k+1)
      ap(k) = (ap_bnds(1,k)+ap_bnds(2,k))*0.5
      b(k) = (b_bnds(1,k)+b_bnds(2,k))*0.5
    enddo
    lev = ap/p0 + b                 ! definition for CMIP purposes
    lev_bnds = ap_bnds/p0 + b_bnds  

    ! register new axes

    id_lev = diag_axis_init('lev', lev, '1', 'Z', &
                            'hybrid sigma pressure coordinate', &
                             direction=-1, set_name='cmip')
    call diag_axis_add_attribute (id_lev, 'formula', 'p(n,k,j,i) = ap(k) + b(k)*ps(n,j,i)')
    call diag_axis_add_attribute (id_lev, 'formula_terms', 'ap: ap b: b ps: ps')
    call diag_axis_add_attribute (id_lev, 'bounds', 'lev_bnds')

    ! vertex number for bounds dimension
    id_nv = diag_axis_init('nv', (/1.,2./), 'none', 'N', 'vertex number', set_name='nv')

    ! register new static variables

    id_ap = register_static_field (mod_name, 'ap', (/id_lev/), &
                                 'vertical coordinate formula term: ap(k)', 'Pa')

    id_b = register_static_field (mod_name, 'b', (/id_lev/), &
                                 'vertical coordinate formula term: b(k)', '1')

    id_ap_bnds = register_static_field (mod_name, 'ap_bnds', (/id_nv,id_lev/), &
                                 'vertical coordinate formula term: ap(k+1/2)', 'Pa')

    id_b_bnds = register_static_field (mod_name, 'b_bnds', (/id_nv,id_lev/), &
                                 'vertical coordinate formula term: b(k+1/2)', '1')

    id_lev_bnds = register_static_field (mod_name, 'lev_bnds', (/id_nv,id_lev/), &
                                        'hybrid sigma pressure coordinate', '1', &
                          standard_name='atmosphere_hybrid sigma pressure coordinate')
    if (id_lev_bnds > 0) then
      call diag_field_add_attribute ( id_lev_bnds, 'formula', 'p(n,k+1/2,j,i) = ap(k+1/2) + b(k+1/2)*ps(n,j,i)')
      call diag_field_add_attribute ( id_lev_bnds, 'formula_terms', 'ap: ap_bnds b: b_bnds ps: ps')
    endif

   ! save static data
   if (id_ap > 0) used = send_data ( id_ap, ap, Time )
   if (id_b  > 0) used = send_data ( id_b , b , Time )
   if (id_ap_bnds  > 0) used = send_data ( id_ap_bnds, ap_bnds, Time )
   if (id_b_bnds   > 0) used = send_data ( id_b_bnds, b_bnds, Time )
   if (id_lev_bnds > 0) used = send_data ( id_lev_bnds, lev_bnds, Time )

!-----------------------------------------------------------------------
! register variables on model levels

    axes3d(3) = id_lev

    id_ta(0) = register_diag_field (mod_name, 'ta', axes3d, Time, &
                                    'Air Temperature', 'K',  &
                                    standard_name='air_temperature', &
                                    area=area_id, missing_value=CMOR_MISSING_VALUE )

    id_ua(0) = register_diag_field (mod_name, 'ua', axes3d, Time, &
                                    'Eastward Wind', 'm s-1',  &
                                    standard_name='eastward_wind', &
                                    area=area_id, missing_value=CMOR_MISSING_VALUE )

    id_va(0) = register_diag_field (mod_name, 'va', axes3d, Time, &
                                    'Northward Wind', 'm s-1',  &
                                    standard_name='northward_wind', &
                                    area=area_id, missing_value=CMOR_MISSING_VALUE )

    id_hus(0) = register_diag_field (mod_name, 'hus', axes3d, Time, &
                                     'Specific Humidity', '1',  &
                                     standard_name='specific_humidity', &
                                     area=area_id, missing_value=CMOR_MISSING_VALUE )


    id_wap(0) = register_diag_field (mod_name, 'wap', axes3d, Time, &
                                     'omega (=dp/dt)', 'Pa s-1',  &
                                     standard_name='lagrangian_tendency_of_air_pressure', &
                                     area=area_id, missing_value=CMOR_MISSING_VALUE )

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
! loop through all possible pressure level sets
! initialize the pressure axis
! define new 3d grid
! define all 3d state variable on this 3d grid

  do id = 1, MAXPLEVS
    if (id .eq. 1) then
      np = num_std_plevs
      pressure_levels(id,1:np) = plev(1:np)
      axis_name = 'plev'
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
    id_plev = diag_axis_init(axis_name, pressure_levels(id,1:np), &
                             'Pa', 'z', 'pressure', direction=-1, set_name="dynamics")
    axes3d(3) = id_plev
    mod_name2 = trim(mod_name)//'_'//trim(axis_name)

    id_ta(id) = register_diag_field (mod_name2, 'ta', axes3d, Time, &
                                     'Air Temperature', 'K',  &
                                     standard_name='air_temperature', &
                                     area=area_id, missing_value=CMOR_MISSING_VALUE )

    id_ua(id) = register_diag_field (mod_name2, 'ua', axes3d, Time, &
                                     'Eastward Wind', 'm s-1',  &
                                     standard_name='eastward_wind', &
                                     area=area_id, missing_value=CMOR_MISSING_VALUE )

    id_va(id) = register_diag_field (mod_name2, 'va', axes3d, Time, &
                                     'Northward Wind', 'm s-1',  &
                                     standard_name='northward_wind', &
                                     area=area_id, missing_value=CMOR_MISSING_VALUE )

    id_hus(id) = register_diag_field (mod_name2, 'hus', axes3d, Time, &
                                     'Specific Humidity', '1',  &
                                     standard_name='specific_humidity', &
                                     area=area_id, missing_value=CMOR_MISSING_VALUE )

    id_hur(id) = register_diag_field (mod_name2, 'hur', axes3d, Time, &
                                     'Relative Humidity', '%',  &
                                     standard_name='relative_humidity', &
                                     area=area_id, missing_value=CMOR_MISSING_VALUE )

    id_wap(id) = register_diag_field (mod_name2, 'wap', axes3d, Time, &
                                     'omega (=dp/dt)', 'Pa s-1',  &
                                     standard_name='lagrangian_tendency_of_air_pressure', &
                                     area=area_id, missing_value=CMOR_MISSING_VALUE )

    id_zg(id) = register_diag_field (mod_name2, 'zg', axes3d, Time, &
                                     'Geopotential Height', 'm',  &
                                     standard_name='geopotential_height', &
                                     area=area_id, missing_value=CMOR_MISSING_VALUE )

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
  if (id_ua(0)  > 0) used = send_data(id_ua(0),  Atm(n)%ua  (isc:iec,jsc:jec,:), Time)
  if (id_va(0)  > 0) used = send_data(id_va(0),  Atm(n)%va  (isc:iec,jsc:jec,:), Time)
  if (id_ta(0)  > 0) used = send_data(id_ta(0),  Atm(n)%pt  (isc:iec,jsc:jec,:), Time)
  if (id_hus(0) > 0) used = send_data(id_hus(0), Atm(n)%q   (isc:iec,jsc:jec,:,sphum), Time)
  if (id_wap(0) > 0) used = send_data(id_wap(0), Atm(n)%omga(isc:iec,jsc:jec,:), Time)
  if (id_ps     > 0) used = send_data(id_ps,     Atm(n)%ps  (isc:iec,jsc:jec), Time)

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

