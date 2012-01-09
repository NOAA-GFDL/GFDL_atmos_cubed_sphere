module fv_climate_nudge_mod

use fms_mod,          only: open_namelist_file, check_nml_error,  &
                            close_file, stdlog, mpp_pe, mpp_root_pe, &
                            write_version_number, string, error_mesg, &
                            FATAL, WARNING, NOTE, file_exist
use mpp_mod,          only: input_nml_file
use diag_manager_mod, only: register_diag_field, send_data,   &
                            register_static_field
use read_climate_nudge_data_mod, only:   &
                            read_climate_nudge_data_init, read_time,  &
                            read_grid, read_climate_nudge_data,  &
                            read_climate_nudge_data_end, &
                            read_sub_domain_init
use time_manager_mod, only: time_type, set_time, print_date,  &
                            operator(<), operator(+)
use  time_interp_mod, only: time_interp
use get_cal_time_mod, only: get_cal_time
use          mpp_mod, only: mpp_min, mpp_max
use    constants_mod, only: RDGAS, RVGAS, PI, KAPPA, CP_AIR
use fv_mapz_mod,      only: mappm
implicit none
private

public :: fv_climate_nudge_init, fv_climate_nudge,  &
          fv_climate_nudge_end, do_ps

character(len=128), parameter :: version = '$Id: fv_climate_nudge.F90,v 19.0 2012/01/06 19:59:01 fms Exp $'
character(len=128), parameter :: tagname = '$Name: siena $'

type var_state_type
   integer :: is, ie, js, je, npz
   integer :: time_level
   real, pointer, dimension(:,:)   :: ps
   real, pointer, dimension(:,:,:) :: u, v, t
   real, pointer, dimension(:,:,:,:) :: q
end type

interface assignment(=)
  module procedure var_state_assignment
end interface

interface remap_xy
   module procedure remap_xy_2d
   module procedure remap_xy_3d
end interface


integer, allocatable :: id1(:,:), id2(:,:), jdc(:,:)
real,    allocatable :: s2c(:,:,:)

integer :: nlon_obs, nlat_obs, nlev_obs, ntime_obs
integer :: jsd, jed
real, allocatable, dimension(:) :: lon_obs, lat_obs, ak_obs, bk_obs
type(time_type), allocatable :: Timelist(:)
type(var_state_type) :: State(2)
logical :: do_state_alloc = .true.
logical :: module_is_initialized = .false.

integer :: freq = 0   ! frequency in seconds
real ::  u_tau = -1.  ! relaxation time in seconds (no insertion if < 0)
real ::  v_tau = -1.
real ::  t_tau = -1.
real ::  q_tau = -1.
real :: ps_tau = -1.
integer :: skip_top_v = 2            ! momentum
integer :: skip_bot_v = 0
integer :: skip_top_t = 0            ! temperature
integer :: skip_bot_t = 21
integer :: skip_top_q = 8            ! specific humidity
integer :: skip_bot_q = 0
logical :: use_pdep_nudge = .false.  ! impose nudging strength that varies with pressure
logical :: use_sub_domain = .false.  ! only read data needed
integer :: verbose = 0               ! 0 .le. verbose .ge. 2

namelist /fv_climate_nudge_nml/ freq, u_tau, v_tau, t_tau, q_tau, ps_tau, &
                                  skip_top_v, skip_bot_v,               &
                                  skip_top_t, skip_bot_t,               &
                                  skip_top_q, skip_bot_q,               &
                                  use_pdep_nudge,                       &
                                  use_sub_domain, verbose

type(time_type) :: Time_next
integer :: id_udt, id_vdt, id_tdt, id_qdt, id_psdt
integer :: id_uerr, id_verr, id_terr, id_qerr, id_pserr
integer :: id_uobs, id_vobs, id_tobs, id_qobs, id_psobs
logical :: do_u, do_v, do_t, do_q, do_ps
logical :: get_wind, get_temp, get_qhum

integer :: id_index, id_coeff

real, parameter :: ZVIR = RVGAS/RDGAS-1.

CONTAINS

!###################################################################################

subroutine fv_climate_nudge_init ( Time, axes, flag )
type (time_type),      intent(in)  :: Time
integer, dimension(3), intent(in)  :: axes
logical, optional,     intent(out) :: flag
integer :: n, ierr, io, unit
character(len=128) :: units, calendar, desc
real(8), allocatable :: tlevels(:)
real :: eps = 1.e-10
real :: missing_value = -1.e10

   if (module_is_initialized) return

 ! read namelist
#ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=fv_climate_nudge_nml, iostat=io)
   ierr = check_nml_error (io, 'fv_climate_nudge_nml')
#else
   if (file_exist('input.nml') ) then
     unit = open_namelist_file()
     ierr=1  
     do while (ierr /= 0)
       read (unit, nml=fv_climate_nudge_nml, iostat=io, end=10) 
       ierr = check_nml_error (io, 'fv_climate_nudge_nml')
     enddo   
10   call close_file (unit)
   endif
#endif

!----- write version and namelist to log file -----

   unit = stdlog()
   call write_version_number (version, tagname)
   if (mpp_pe() == mpp_root_pe()) write (unit, nml=fv_climate_nudge_nml)

 ! initialize flags
   do_u  = .false.; if ( u_tau > -eps) do_u  = .true.
   do_v  = .false.; if ( v_tau > -eps) do_v  = .true.
   do_t  = .false.; if ( t_tau > -eps) do_t  = .true.
   do_q  = .false.; if ( q_tau > -eps) do_q  = .true.
   do_ps = .false.; if (ps_tau > -eps) do_ps = .true.

 ! namelist dummy checks
 ! if no overrides turned on then set freq = 0
   if (freq > 0) then
       if ( .not.do_u .and. .not.do_v .and. .not.do_t .and. &
            .not.do_q .and. .not.do_ps ) then
            call error_mesg ('fv_climate_nudge_mod', 'no variables specified '//&
                             'for override', WARNING)
       endif
   else
       if ( do_u .or. do_v .or. do_t .or.  do_q .or. do_ps ) then
            call error_mesg ('fv_climate_nudge_mod', 'variables specified '//&
                             'for override when freq = 0', FATAL)
       endif
       freq = 0
   endif

 ! return flag = true when override is needed
   if (present(flag)) then
       flag = freq .gt. 0
   endif

 ! what is the next time for data insertion

   Time_next = Time + set_time(freq)

 ! initialize diagnostics

   desc = ' tendency due to data insertion'

   id_udt = register_diag_field ('atmos_nudge', 'udt_nudge', axes, Time, &
                                 'zonal wind'//trim(desc), 'm/s2', missing_value=missing_value)
   id_vdt = register_diag_field ('atmos_nudge', 'vdt_nudge', axes, Time, &
                                 'meridional wind'//trim(desc), 'm/s2',missing_value=missing_value)
   id_tdt = register_diag_field ('atmos_nudge', 'tdt_nudge', axes, Time, &
                                 'temperature'//trim(desc), 'degK/s',missing_value=missing_value)
   id_qdt = register_diag_field ('atmos_nudge', 'qdt_nudge', axes, Time, &
                                 'specific humidity'//trim(desc), 'kg/kg/s',missing_value=missing_value)
   id_psdt = register_diag_field ('atmos_nudge', 'psdt_nudge', axes(1:2), Time, &
                                 'surface pressure'//trim(desc), 'Pa/s',missing_value=missing_value)

   desc = ' bias'

   id_uerr = register_diag_field ('atmos_nudge', 'u_bias', axes, Time, &
                                 'zonal wind'//trim(desc), 'm/s2', missing_value=missing_value)
   id_verr = register_diag_field ('atmos_nudge', 'v_bias', axes, Time, &
                                 'meridional wind'//trim(desc), 'm/s2',missing_value=missing_value)
   id_terr = register_diag_field ('atmos_nudge', 't_bias', axes, Time, &
                                 'temperature'//trim(desc), 'degK/s',missing_value=missing_value)
   id_qerr = register_diag_field ('atmos_nudge', 'q_bias', axes, Time, &
                                 'specific humidity'//trim(desc), 'kg/kg/s',missing_value=missing_value)
   id_pserr = register_diag_field ('atmos_nudge', 'ps_bias', axes(1:2), Time, &
                                 'surface pressure'//trim(desc), 'Pa/s',missing_value=missing_value)

   desc = ' observation used for nudging'

   id_uobs = register_diag_field ('atmos_nudge', 'u_obs', axes, Time, &
                                 'zonal wind'//trim(desc), 'm/s2', missing_value=missing_value)
   id_vobs = register_diag_field ('atmos_nudge', 'v_obs', axes, Time, &
                                 'meridional wind'//trim(desc), 'm/s2',missing_value=missing_value)
   id_tobs = register_diag_field ('atmos_nudge', 't_obs', axes, Time, &
                                 'temperature'//trim(desc), 'degK/s',missing_value=missing_value)
   id_qobs = register_diag_field ('atmos_nudge', 'q_obs', axes, Time, &
                                 'specific humidity'//trim(desc), 'kg/kg/s',missing_value=missing_value)
   id_psobs = register_diag_field ('atmos_nudge', 'ps_obs', axes(1:2), Time, &
                                 'surface pressure'//trim(desc), 'Pa/s',missing_value=missing_value)

   id_index = register_static_field ('atmos_nudge', 'jdc', axes(1:2), 'interp index', 'none' )
   id_coeff = register_static_field ('atmos_nudge', 's2c', axes(1:2), 'interp coeff', 'none' )

 ! set flags
   get_wind = do_u .or. id_udt>0 .or. id_uerr>0 .or. id_uobs>0 .or. do_v .or. id_vdt>0 .or. id_verr>0 .or. id_vobs>0
   get_temp = do_t .or. id_tdt>0 .or. id_terr>0 .or. id_tobs>0
   get_qhum = do_q .or. id_qdt>0 .or. id_qerr>0 .or. id_qobs>0

   if (.not.get_wind .and. .not.get_temp .and. .not.get_qhum .and. freq.eq.0) return
!--------------------------------------------
! initialize input data from file
! get the size of the global data
  call error_mesg ('fv_climate_nudge_mod', 'initializing nudging', NOTE)
  call read_climate_nudge_data_init (nlon_obs, nlat_obs, nlev_obs, ntime_obs)
  if (verbose .gt. 1 .and. mpp_pe() .eq. mpp_root_pe()) then
     print '(a,4i10)', 'fv_climate_nudge: nlon_obs, nlat_obs, nlev_obs, ntime_obs = ', &
                        nlon_obs, nlat_obs, nlev_obs, ntime_obs
  endif

! read the time level information
  allocate ( tlevels (ntime_obs) )
  allocate ( TimeList(ntime_obs) )
  call read_time ( tlevels, units, calendar )
  do n = 1, ntime_obs
     Timelist(n) = get_cal_time( tlevels(n), trim(units), trim(calendar), &
                                 permit_calendar_conversion=.true. )
  enddo
  deallocate ( tlevels )

! read the grid information
  allocate ( lon_obs(nlon_obs), lat_obs(nlat_obs), ak_obs(nlev_obs+1), bk_obs(nlev_obs+1) )
  call read_grid ( lon_obs, lat_obs, ak_obs, bk_obs )

  if (mpp_pe() .eq. mpp_root_pe()) then
     if (verbose .gt. 1) then
        print *, 'fv_climate_nudge: ak_obs=',ak_obs
        print *, 'fv_climate_nudge: bk_obs=',bk_obs
     endif
     if (verbose .gt. 0) then
        call print_date (Timelist(1),'nudging data, start date: ')
        call print_date (Timelist(ntime_obs),'nudging data, final date: ')
     endif
  endif

  module_is_initialized = .true.

end subroutine fv_climate_nudge_init

!###################################################################################

subroutine fv_climate_nudge (Time, dt, is, ie, js, je, npz, pfull, &
                                   lon, lat, phis, ak, bk,        &
                                   ps, u, v, t, q, psdt, udt, vdt, tdt, qdt )
type(time_type), intent(in) :: Time
real,            intent(in) :: dt
integer,         intent(in) :: is, ie, js, je, npz

real, intent(in)    :: phis(is:ie,js:je)
real, intent(in)    :: lon (is:ie,js:je)
real, intent(in)    :: lat (is:ie,js:je)
real, intent(in)    :: ak  (npz+1)
real, intent(in)    :: bk  (npz+1)
real, intent(in)    :: pfull (npz)

real, intent(inout) :: ps  (is:ie,js:je)
real, intent(inout) :: psdt(is:ie,js:je)

real, intent(inout) :: u(is:ie,js:je,npz)
real, intent(inout) :: v(is:ie,js:je,npz)
real, intent(inout) :: t(is:ie,js:je,npz)
real, intent(inout) :: q(is:ie,js:je,npz,1)

real, intent(inout) :: udt(is:ie,js:je,npz)
real, intent(inout) :: vdt(is:ie,js:je,npz)
real, intent(inout) :: tdt(is:ie,js:je,npz)
real, intent(inout) :: qdt(is:ie,js:je,npz,1)

! local arrays
real, dimension(is:ie,js:je,nlev_obs) :: u_obs, v_obs, t_obs
real, dimension(is:ie,js:je,nlev_obs) :: q_obs
real, dimension(is:ie,js:je) :: ps_obs, phis_obs

real, dimension(is:ie,js:je,nlev_obs+1) :: phaf_obs, lphaf_obs
real, dimension(is:ie,js:je,npz+1)      :: phaf, lphaf

!real, dimension(nlon_obs,jsd:jed,nlev_obs) :: dat3
real, allocatable :: dat3(:,:,:)

real, dimension(is:ie,js:je,npz) :: obs, tend
real    :: factor(npz,3), wght(2)
logical :: sent
integer :: itime(2), k, n
character(len=128) :: str
character(len=8)   :: rstr(2)
! this variable addresses a bug with high resolution GFS analyses
! before a certain date preprocessing scripts did not convert virtual temperature to temperature
logical :: virtual_temp_obs = .false.

   if (.not.module_is_initialized) then
       call error_mesg ('fv_climate_nudge_mod', 'module not initialized', FATAL)
   endif

 ! is it time for data forcing
   if (Time < Time_next) then
       return
   endif
   Time_next = Time_next + set_time(freq)

 ! only one tracer (specific humidity) can be interpolated
 ! if (size(q,4).ne.1 .or. size(qdt,4).ne.1) then
 !     call error_mesg ('fv_climate_nudge_mod', 'only one tracer (sphum) allowed', FATAL)
 ! endif

 ! vertically dependent factor
   call get_factor (npz,pfull, factor)
 ! first time allocate state 
   if (do_state_alloc) then
      call var_state_init ( is, ie, js, je, npz, State(1) )
      call var_state_init ( is, ie, js, je, npz, State(2) )

      ! sub-domain initialization
      if (use_sub_domain) then
         call read_sub_domain_init (minval(lat), maxval(lat), lat_obs, jsd, jed)
      else
         jsd=1; jed=nlat_obs
      endif

      ! also initalize coefficient for horizontal interpolation
      allocate ( id1(is:ie,js:je), id2(is:ie,js:je), jdc(is:ie,js:je), s2c(is:ie,js:je,4) )
      call remap_coef ( 1, nlon_obs, jsd, jed, lon_obs, lat_obs(jsd:jed), &
                        is, ie, js, je, lon, lat, id1, id2, jdc, s2c )
      do_state_alloc = .false.
      if (id_index > 0) sent = send_data (id_index, real(jdc), Time)
      if (id_coeff > 0) sent = send_data (id_coeff, real(s2c(:,:,1)), Time)
   endif

!//////////////////////////////////////////////////////

  ! get the time indices needed
    call time_interp (Time, Timelist, wght(2), itime(1), itime(2))
    wght(1) = 1. - wght(2)
    if (verbose .gt. 1 .and. mpp_pe() .eq. mpp_root_pe()) then
       write (rstr(1),'(f8.6)') wght(1)
       write (rstr(2),'(f8.6)') wght(2)
       str = 'Data Nudging: itime='//trim(string(itime(1)))//','//trim(string(itime(2)))// &
             '; wght='//rstr(1)//','//rstr(2)//'; Date: '
       call print_date (Time,trim(str))
    endif

    do n = 1, 2

       if (itime(n) .ne. State(n)%time_level) then
          if (n .eq. 1) then
             if (itime(1) .eq. State(2)%time_level) then
                 State(1) = State(2)
                 cycle
             endif
          endif
          if (.not.allocated(dat3)) allocate (dat3(nlon_obs,jsd:jed,nlev_obs))

         ! ---- horizontal interpolation ----
         ! geop hght
           call read_climate_nudge_data (itime(n), 'phis', dat3(:,:,1), 1, jsd)
           if (verbose .gt. 1) call prt_minmax_2d('phis_obs(mn,mx)=',dat3(:,:,1))
           call remap_xy (1, nlon_obs, jsd, jed, dat3(:,:,1), is, ie, js, je, id1, id2, jdc, s2c, phis_obs)
         ! surf pres
           call read_climate_nudge_data (itime(n), 'psrf', dat3(:,:,1), 1, jsd)
           if (verbose .gt. 1) call prt_minmax_2d('ps_obs(mn,mx)=',dat3(:,:,1))
           call remap_xy (1, nlon_obs, jsd, jed, dat3(:,:,1), is, ie, js, je, id1, id2, jdc, s2c, ps_obs)
          !if (verbose .gt. 1) call prt_minmax_2d('ps_obs_adj(mn,mx)=',ps_obs)

         ! compute pressure and ln pres at layer interfaces
           do k = 1, nlev_obs+1
              phaf_obs(:,:,k) = ak_obs(k) + bk_obs(k)*ps_obs
           enddo
           do k = 2, nlev_obs+1
              lphaf_obs(:,:,k) = log(phaf_obs(:,:,k))
           enddo
           where (phaf_obs(:,:,1) .gt. 0.0)
              lphaf_obs(:,:,1) = log(phaf_obs(:,:,1))
           elsewhere
              lphaf_obs(:,:,1) = lphaf_obs(:,:,2)-2.
           endwhere
           lphaf_obs(:,:,1) = max(lphaf_obs(:,:,1),0.0)

         ! wind interpolation
           if (get_wind) then
             call read_climate_nudge_data (itime(n), 'uwnd', dat3, 1, jsd)
             if (verbose .gt. 1) call prt_minmax_3d('u_obs(mn,mx)=',dat3)
             call remap_xy (1, nlon_obs, jsd, jed, nlev_obs, dat3, is, ie, js, je, id1, id2, jdc, s2c, u_obs)
             call read_climate_nudge_data (itime(n), 'vwnd', dat3, 1, jsd)
             if (verbose .gt. 1) call prt_minmax_3d('v_obs(mn,mx)=',dat3)
             call remap_xy (1, nlon_obs, jsd, jed, nlev_obs, dat3, is, ie, js, je, id1, id2, jdc, s2c, v_obs)
           endif
         ! spec hum (always need for virt temp)
         ! if (get_qhum) then
             call read_climate_nudge_data (itime(n), 'qhum', dat3, 1, jsd)
             if (verbose .gt. 1) call prt_minmax_3d('q_obs(mn,mx)=',dat3)
             call remap_xy (1, nlon_obs, jsd, jed, nlev_obs, dat3, is, ie, js, je, id1, id2, jdc, s2c, q_obs)
         ! endif
         ! temperature (always need for surf pres remapping)
         ! if (get_temp) then
             call read_climate_nudge_data (itime(n), 'temp', dat3, 1, jsd)
             if (verbose .gt. 1) call prt_minmax_3d('t_obs(mn,mx)=',dat3)
             call remap_xy (1, nlon_obs, jsd, jed, nlev_obs, dat3, is, ie, js, je, id1, id2, jdc, s2c, t_obs)
             if (.not.virtual_temp_obs) then
                t_obs = t_obs*(1.+ZVIR*q_obs) ! virtual effect
             endif
         ! endif

         ! VERTICAL REMAPPING
         ! remap surface pressure for different surface height
           call remap_ps (is, ie, js, je, nlev_obs, phis_obs, phaf_obs, lphaf_obs, t_obs, phis, State(n)%ps)

         ! compute pressure and ln pres at layer interfaces
          do k = 1, npz+1
              phaf(:,:,k) = ak(k) + bk(k)*State(n)%ps
             lphaf(:,:,k) = log(phaf(:,:,k))
          enddo

          if (get_wind) then
             call remap_3d (is, ie, js, je, nlev_obs, npz, phaf_obs, u_obs, phaf, State(n)%u, -1)
             call remap_3d (is, ie, js, je, nlev_obs, npz, phaf_obs, v_obs, phaf, State(n)%v, -1)
          endif
          if (get_qhum .or. get_temp) then
             call remap_3d (is, ie, js, je, nlev_obs, npz, phaf_obs, q_obs, phaf, State(n)%q(:,:,:,1), 0)
          endif
          if (get_temp) then
             ! use logp
             call remap_3d (is, ie, js, je, nlev_obs, npz, lphaf_obs, t_obs, lphaf, State(n)%t, 1)
             State(n)%t = State(n)%t/(1.+ZVIR*State(n)%q(:,:,:,1)) ! virtual effect
          endif

          State(n)%time_level = itime(n)
          if (verbose .gt. 1) then
             !call prt_minmax_2d('phis(mn,mx)=',phis)
             call prt_minmax_2d('ps(mn,mx)=',State(n)%ps)
             if (get_wind) call prt_minmax_3d('u(mn,mx)=',State(n)%u)
             if (get_wind) call prt_minmax_3d('v(mn,mx)=',State(n)%v)
             if (get_temp) call prt_minmax_3d('t(mn,mx)=',State(n)%t)
             if (get_qhum) call prt_minmax_3d('q(mn,mx)=',State(n)%q(:,:,:,1))
          endif
       endif
    enddo
    if (allocated(dat3)) deallocate (dat3)

!//////////////////////////////////////////////////////

! zonal wind component
  if (do_u .or. id_udt>0 .or. id_uerr>0 .or. id_uobs>0) then
     obs = State(1)%u*wght(1) + State(2)%u*wght(2)
     if (do_u) then
        do k = 1, npz
           tend(:,:,k) = (obs(:,:,k) - u(:,:,k)) / (u_tau + dt) * factor(k,1)
           u(:,:,k) = u(:,:,k) + dt*tend(:,:,k)
           udt(:,:,k) = udt(:,:,k) + tend(:,:,k)
        enddo
     else if (id_udt > 0) then
        tend = 0.0
     endif
     if (id_udt > 0) sent = send_data (id_udt, tend, Time)
     if (id_uerr > 0) then
        tend = u - obs
        sent = send_data (id_uerr, tend, Time)
     endif
     if (id_uobs > 0) sent = send_data (id_uobs, obs, Time)
  endif

! meridional wind component
  if (do_v .or. id_vdt>0 .or. id_verr>0 .or. id_vobs>0) then
     obs = State(1)%v*wght(1) + State(2)%v*wght(2)
     if (do_v) then
        do k = 1, npz
           tend(:,:,k) = (obs(:,:,k) - v(:,:,k)) / (v_tau + dt) * factor(k,1)
           v(:,:,k) = v(:,:,k) + dt*tend(:,:,k)
           vdt(:,:,k) = vdt(:,:,k) + tend(:,:,k)
        enddo
     else if (id_vdt > 0) then
        tend = 0.0
     endif
     if (id_vdt > 0) sent = send_data (id_vdt, tend, Time)
     if (id_verr > 0) then
        tend = v - obs
        sent = send_data (id_verr, tend, Time)
     endif
     if (id_vobs > 0) sent = send_data (id_vobs, obs, Time)
  endif

! temperature
  if (do_t .or. id_tdt>0 .or. id_terr>0 .or. id_tobs>0) then
     obs = State(1)%t*wght(1) + State(2)%t*wght(2)
     if (do_t) then
        do k = 1, npz
           tend(:,:,k) = (obs(:,:,k) - t(:,:,k)) / (t_tau + dt) * factor(k,2)
           t(:,:,k) = t(:,:,k) + dt*tend(:,:,k)
           tdt(:,:,k) = tdt(:,:,k) + tend(:,:,k)
        enddo
     else if (id_tdt > 0) then
        tend = 0.0
     endif
     if (id_tdt > 0) sent = send_data (id_tdt, tend, Time)
     if (id_terr > 0) then
        tend = t - obs
        sent = send_data (id_terr, tend, Time)
     endif
     if (id_tobs > 0) sent = send_data (id_tobs, obs, Time)
  endif

! specific humidity
  if (do_q .or. id_qdt>0 .or. id_qerr>0 .or. id_qobs>0) then
     obs = State(1)%q(:,:,:,1)*wght(1) + State(2)%q(:,:,:,1)*wght(2)
     if (do_q) then
        do k = 1, npz
           tend(:,:,k) = (obs(:,:,k) - q(:,:,k,1)) / (q_tau + dt) * factor(k,2)
           q(:,:,k,1) = q(:,:,k,1) + dt*tend(:,:,k)
           qdt(:,:,k,1) = qdt(:,:,k,1) + tend(:,:,k)
        enddo
     else if (id_qdt > 0) then
        tend = 0.0
     endif
     if (id_qdt > 0) sent = send_data (id_qdt, tend, Time)
     if (id_qerr > 0) then
        tend = q(:,:,:,1) - obs
        sent = send_data (id_qerr, tend, Time)
     endif
     if (id_qobs > 0) sent = send_data (id_qobs, obs, Time)
  endif

! surface pressure
  if (do_ps .or. id_psdt>0 .or. id_pserr>0 .or. id_psobs>0) then
     obs(:,:,1) = State(1)%ps*wght(1) + State(2)%ps*wght(2)
     if (do_ps) then
        tend(:,:,1) = (obs(:,:,1) - ps) / (ps_tau + dt)
        ps = ps + dt*tend(:,:,1)
        psdt = psdt + tend(:,:,1)
     else if (id_psdt > 0) then
        tend(:,:,1) = 0.0
     endif
     if (id_psdt > 0) sent = send_data (id_psdt, tend(:,:,1), Time)
     if (id_pserr > 0) then
        tend(:,:,1) = ps - obs(:,:,1)
        sent = send_data (id_pserr, tend(:,:,1), Time)
     endif
     if (id_psobs > 0) sent = send_data (id_psobs, obs(:,:,1), Time)
  endif

end subroutine fv_climate_nudge

!###################################################################################

subroutine get_factor (nlev,pfull,factor)
integer, intent(in)  :: nlev
real,    intent(in)  :: pfull(nlev)
real,    intent(out) :: factor(nlev,3)
integer :: k
real    :: psurf

! vertically dependent Tau
! very crude - zero at top/bottom + linear increase with level downward/upward

   factor = 1.

!------------------------------------------------------------------
! Momentum:
   if (skip_top_v > 0) then
!++amf
!      factor(1,1) = 0.
!      do k = 2, skip_top_v
!         factor(k,1) = factor(k-1,1) + 1./real(skip_top_v)
!      enddo
    factor(skip_top_v+1,1) = 0.25
    factor(skip_top_v+2,1) = 0.5
    do k = 1, skip_top_v
       factor(k,1) = 0.
    enddo
   endif
   if (skip_bot_v > 0) then
      factor(nlev,1) = 0.
      do k = nlev-1, nlev-skip_bot_v+1, -1
         factor(k,1) = factor(k+1,1) + 1./real(skip_bot_v)
      enddo
   endif

! temperature
   if (skip_top_t > 0) then
!++amf
!      factor(1,2) = 0.
!      do k = 2, skip_top_t
!         factor(k,2) = factor(k-1,2) + 1./real(skip_top_t)
!      enddo
    factor(skip_top_t+1,2) = 0.25
    factor(skip_top_t+2,2) = 0.5
    do k = 1, skip_top_t
       factor(k,2) = 0.
    enddo
!--amf
   endif
   if (skip_bot_t > 0) then
         factor(nlev-skip_bot_t-1,2) = 0.5
         factor(nlev-skip_bot_t,  2) = 0.25
      do k=nlev-skip_bot_t+1,nlev
         factor(k,2) = 0.
      enddo
   endif
   
! Specific humidity
   if (skip_top_q > 0) then
      do k = 1, skip_top_q
         factor(k,3) = 0.
      enddo
         factor(skip_top_q+1,3) = 0.25
         factor(skip_top_q+2,3) = 0.5
   endif
   if (skip_bot_q > 0) then
      factor(nlev,3) = 0.
      do k = nlev-1, nlev-skip_bot_q+1, -1
         factor(k,3) = factor(k+1,3) + 1./real(skip_bot_q)
      enddo
   endif

!++amf Apply scaling such that nudging strength falls off with pressure.

   if (use_pdep_nudge) then
      psurf = pfull(nlev)
      do k = 1, nlev
        factor(k,:) = factor(k,:)*(pfull(k)/psurf)
 !       print*, 'AMF: k,pfull(k),psurf, factor = ', k, pfull(k), psurf, pfull(k)/psurf
      enddo
   endif

!--amf

!------------------------------------------------------------------

end subroutine get_factor

!###################################################################################

subroutine var_state_init ( is, ie, js, je, npz, State )
type(var_state_type), intent(inout) :: State
integer, intent(in) :: is, ie, js, je, npz

   State%is = is
   State%ie = ie
   State%js = js
   State%je = je
   State%npz = npz
   State%time_level = -1
   allocate (State%ps(is:ie,js:je))
   if (get_wind) then
      allocate (State%u(is:ie,js:je,1:npz))
      allocate (State%v(is:ie,js:je,1:npz))
   endif
   if (get_temp) then
      allocate (State%t(is:ie,js:je,1:npz))
   endif
   if (get_qhum .or. get_temp) then
      allocate (State%q(is:ie,js:je,1:npz,1:1))
   endif

end subroutine var_state_init

!-----------------------------------------------------

subroutine var_state_assignment (State1,State2)
type(var_state_type), intent(out) :: State1
type(var_state_type), intent(in)  :: State2

 ! resolution must match
   if (State1%is /= State2%is .or. State1%ie /= State2%ie .or. &
       State1%js /= State2%js .or. State1%je /= State2%je .or. &
       State1%npz /= State2%npz) then
            call error_mesg ('fv_climate_nudge_mod', 'invalid var_state assignment: '// &
                             'dimensions must match', FATAL)
   endif

 ! copy data - must be allocated
   State1%time_level = State2%time_level
   State1%ps = State2%ps
   if (associated(State1%u))  State1%u = State2%u
   if (associated(State1%v))  State1%v = State2%v
   if (associated(State1%t))  State1%t = State2%t
   if (associated(State1%q))  State1%q = State2%q

end subroutine var_state_assignment

!-----------------------------------------------------

subroutine var_state_del ( State )
type(var_state_type), intent(inout) :: State

   State%is = 1
   State%ie = 0
   State%js = 1
   State%je = 0
   State%npz = 0
   State%time_level = -1
   deallocate (State%ps)
   if (get_wind) then
      if (associated(State%u)) deallocate (State%u)
      if (associated(State%v)) deallocate (State%v)
   endif
   if (get_temp) then
      if (associated(State%t)) deallocate (State%t)
   endif
   if (get_qhum .or. get_temp) then
      if (associated(State%q)) deallocate (State%q)
   endif

end subroutine var_state_del

!###################################################################################

subroutine fv_climate_nudge_end

  if (.not.module_is_initialized) return

  deallocate ( TimeList )
  deallocate ( lon_obs, lat_obs, ak_obs, bk_obs )

  if (.not.do_state_alloc) then
     call var_state_del ( State(1) )
     call var_state_del ( State(2) )
     deallocate ( id1, id2, jdc, s2c )
     do_state_alloc = .true.
  endif

  call read_climate_nudge_data_end

  module_is_initialized = .false.
end subroutine fv_climate_nudge_end

!###################################################################################

subroutine prt_minmax_2d (str,a)
character(len=*), intent(in) :: str
real, intent(in) :: a(:,:)
real :: local_min, local_max

   local_min = minval(a)
   local_max = maxval(a)
   call mpp_min(local_min)
   call mpp_max(local_max)
   if (mpp_pe()==mpp_root_pe()) then
      print *, trim(str),local_min,local_max
   endif

end subroutine prt_minmax_2d

!------------------------------------------

subroutine prt_minmax_3d (str,a)
character(len=*), intent(in) :: str
real, intent(in) :: a(:,:,:)
real :: local_min, local_max

   local_min = minval(a)
   local_max = maxval(a)
   call mpp_min(local_min)
   call mpp_max(local_max)
   if (mpp_pe()==mpp_root_pe()) then
      print *, trim(str),local_min,local_max
   endif

end subroutine prt_minmax_3d

!###################################################################################
!###################################################################################
! For data nudging with the FV cubed sphere dynamical core:
! Contact S.-J. Lin, NOAA/GFDL, for more information

  subroutine remap_coef( isd, ied, jsd, jed, lon_in, lat_in, &
                         is, ie, js, je, lon_out, lat_out,   &
                         id1, id2, jdc, s2c )
!--------
! Input:
!--------
! Data Input: data must be global
  integer, intent(in):: isd, ied                  ! Data x-dimension (W->E; must be periodic)
  integer, intent(in):: jsd, jed                  ! Data y-dimension (S->N)
  real,    intent(in):: lon_in(isd:ied)           ! Data longitude (Radian; periodic)
  real,    intent(in):: lat_in(jsd:jed)           ! Data latitude (increases from SP to NP)

! Model input:
  integer, intent(in):: is, ie, js, je        ! model horizontal dimensions (un-ghosted sub-domian)
  real,    intent(in):: lon_out(is:ie,js:je)  ! model longitude (Radian)
  real,    intent(in):: lat_out(is:ie,js:je)  ! model latitude (Radian)

!--------
! Output:
!--------
!
  integer, intent(out), dimension(is:ie,js:je  ):: id1, id2, jdc
  real,    intent(out), dimension(is:ie,js:je,4):: s2c
 
!===============================================================================================

! local:
  real:: rdlon(isd:ied)
  real:: rdlat(jsd:jed)
  real:: a1, b1
  integer i, j, i1, i2, jc, i0, j0

 !pk0(1) = ak_in(1)**KAPPA 
 !pn_top = log(ak_in(1))

  do i=isd,ied-1
     rdlon(i) = 1. / (lon_in(i+1) - lon_in(i))
  enddo
     rdlon(ied) = 1. / (lon_in(isd) + 2.*PI - lon_in(ied))  ! periodic assumption

  do j=jsd,jed-1
     rdlat(j) = 1. / (lat_in(j+1) - lat_in(j))
  enddo

! * Interpolate to cubed sphere cell center
  do 5000 j=js,je

     do i=is,ie

       if ( lon_out(i,j) .gt. lon_in(ied) ) then
            i1 = ied;     i2 = isd
            a1 = (lon_out(i,j)-lon_in(ied)) * rdlon(ied)
       elseif ( lon_out(i,j) .lt. lon_in(1) ) then
            i1 = ied;     i2 = isd
            a1 = (lon_out(i,j)+2.*PI-lon_in(ied)) * rdlon(ied)
       else
            do i0=isd,ied-1
            if ( lon_out(i,j) .ge. lon_in(i0) .and. lon_out(i,j) .le. lon_in(i0+1) ) then
               i1 = i0;  i2 = i0+1
               a1 = (lon_out(i,j)-lon_in(i1)) * rdlon(i0)
               go to 111
            endif
            enddo
       endif

111    continue

       if ( lat_out(i,j) .lt. lat_in(jsd) ) then
            jc = jsd
            b1 = 0.
       elseif ( lat_out(i,j) .gt. lat_in(jed) ) then
            jc = jed-1
            b1 = 1.
       else
          do j0=jsd,jed-1
          if ( lat_out(i,j) .ge. lat_in(j0) .and. lat_out(i,j) .le. lat_in(j0+1) ) then
               jc = j0
               b1 = (lat_out(i,j)-lat_in(jc)) * rdlat(jc)
               go to 222
          endif
          enddo
       endif
222    continue

! Debug codes:
!      if ( a1<0.0 .or. a1>1.0 .or.  b1<0.0 .or. b1>1.0 ) then
!           write(*,*) i,j,a1, b1
!      endif

       s2c(i,j,1) = (1.-a1) * (1.-b1)
       s2c(i,j,2) =     a1  * (1.-b1)
       s2c(i,j,3) =     a1  *     b1
       s2c(i,j,4) = (1.-a1) *     b1
       id1(i,j) = i1
       id2(i,j) = i2
       jdc(i,j) = jc

     enddo
5000 continue

  end subroutine remap_coef

!---------------------------------------------------

  subroutine remap_xy_3d( isd, ied, jsd, jed, km, q_in, &
                          is, ie, js, je, id1, id2, jdc, s2c, &
                          q_out )
!--------
! Input:
!--------
! Data Input: data must be global
  integer, intent(in):: isd, ied                  ! Data x-dimension (W->E; must be periodic)
  integer, intent(in):: jsd, jed                  ! Data y-dimension (S->N)
  integer, intent(in):: km                        ! Data vertical-dimension
  real,    intent(in), dimension(isd:ied,jsd:jed,km) :: q_in    ! Data on input grid

! Model input:
  integer, intent(in):: is, ie, js, je        ! model horizontal dimensions (un-ghosted sub-domian)

! Indices and coefficients for bilinear interpolation
  integer, intent(in), dimension(is:ie,js:je  ):: id1, id2, jdc
  real,    intent(in), dimension(is:ie,js:je,4):: s2c

!--------
! Output:
!--------
! Data Output: on local model horizontal grid and input data vertical grid
  real,    intent(out), dimension(is:ie,js:je,km):: q_out

  integer :: i, j, k

    do k=1,km
    do j=js,je
    do i=is,ie
        q_out(i,j,k) = s2c(i,j,1)*q_in(id1(i,j),jdc(i,j),  k) + s2c(i,j,2)*q_in(id2(i,j),jdc(i,j),  k) +  &
                       s2c(i,j,3)*q_in(id2(i,j),jdc(i,j)+1,k) + s2c(i,j,4)*q_in(id1(i,j),jdc(i,j)+1,k)
    enddo
    enddo
    enddo

  end subroutine remap_xy_3d

!---------------------------------------------------

  subroutine remap_xy_2d( isd, ied, jsd, jed, q_in, &
                          is, ie, js, je, id1, id2, jdc, s2c, &
                          q_out )
!--------
! Input:
!--------
  integer, intent(in):: isd, ied                            ! Data x-dimension (W->E; must be periodic)
  integer, intent(in):: jsd, jed                            ! Data y-dimension (S->N)
  real,    intent(in), dimension(isd:ied,jsd:jed) :: q_in   ! Data on input grid
  integer, intent(in):: is, ie, js, je                      ! model horizontal dimensions (un-ghosted sub-domian)
  integer, intent(in), dimension(is:ie,js:je  ):: id1, id2, jdc
  real,    intent(in), dimension(is:ie,js:je,4):: s2c
!--------
! Output:
!--------
  real,    intent(out), dimension(is:ie,js:je):: q_out

! Local:
  real, dimension(isd:ied,jsd:jed,1) :: q3d_in
  real, dimension(is:ie,js:je,1)     :: q3d_out

     q3d_in(:,:,1) = q_in
     call remap_xy_3d( isd, ied, jsd, jed, 1, q3d_in, &
                          is, ie, js, je, id1, id2, jdc, s2c, &
                          q3d_out )
     q_out = q3d_out(:,:,1)

  end subroutine remap_xy_2d

!---------------------------------------------------

  subroutine remap_ps( is, ie, js, je, km, &
                       gz_dat, ph_dat, pn_dat, tp_dat, phis, ps )

!--------
! Input:
!--------
  integer, intent(in):: is, ie, js, je           ! model horizontal dimensions (un-ghosted sub-domian)
  integer, intent(in):: km                       ! model horizontal dimensions (un-ghosted sub-domian)
    real,  intent(in):: gz_dat(is:ie,js:je)      ! Data surface geop height (m2/s2)
    real,  intent(in):: ph_dat(is:ie,js:je,km+1) ! Data pressure at layer interfaces (Pa)
    real,  intent(in):: pn_dat(is:ie,js:je,km+1) ! Data natural log pressure at layer interfaces (Pa)
    real,  intent(in):: tp_dat(is:ie,js:je,km)   ! Data temperature in layers (K)
    real,  intent(in):: phis  (is:ie,js:je)      ! Model surface geop height (m2/s2)

!--------
! Output:
!--------
  real,    intent(out):: ps (is:ie,js:je)      ! Model surface pressure (pa)

! local
  integer :: i, j, k
  real    :: pk0(km+1), gz(km+1), pt0, pst
! needed:  real, parameter :: kappa, rdgas, cp_air

  do j = js,je
  do i = is,ie

! Adjust interpolated ps to model terrain
       gz(km+1) = gz_dat(i,j)
       pk0(km+1) = ph_dat(i,j,km+1)**KAPPA
       do k=km,1,-1
           gz(k) = gz(k+1) + RDGAS*tp_dat(i,j,k)*(pn_dat(i,j,k+1)-pn_dat(i,j,k)) 
           pk0(k) = ph_dat(i,j,k)**KAPPA
       enddo
       if ( phis(i,j) .gt. gz_dat(i,j) ) then
           do k=km,1,-1
              if( phis(i,j) <  gz(k)  .and.    &
                  phis(i,j) >= gz(k+1) ) then
                  pst = pk0(k) + (pk0(k+1)-pk0(k))*(gz(k)-phis(i,j))/(gz(k)-gz(k+1))
                  go to 123
              endif
           enddo
       else
! Extrapolation into the ground
! lowest layer potential temp is needed
           pt0= tp_dat(i,j,km)/(pk0(km+1)-pk0(km))*(KAPPA*(pn_dat(i,j,km+1)-pn_dat(i,j,km)))
           pst = pk0(km+1) + (gz_dat(i,j)-phis(i,j))/(CP_AIR*pt0)
       endif
  123  ps(i,j) = pst**(1./KAPPA)

  enddo
  enddo


  end subroutine remap_ps

!---------------------------------------------------

  subroutine remap_3d( is, ie, js, je, km, npz, &
                       pe0, qn0, pe1, qn1, n )

!--------
! Input:
!--------
  integer, intent(in):: is, ie, js, je         ! model horizontal dimensions (un-ghosted sub-domian)
  integer, intent(in):: km                     ! vertical dimension for input data
  integer, intent(in):: npz                    ! vertical dimension for model data
    real,  intent(in):: pe0(is:ie,js:je,km+1)  ! pressure at layer interfaces for input data
    real,  intent(in):: qn0(is:ie,js:je,km)    ! scalar quantity on input data levels
    real,  intent(in):: pe1(is:ie,js:je,npz+1) ! pressure at layer interfaces for model data
  integer, intent(in):: n                      ! -1 wind; 0 sphum; +1 ptemp

!--------
! Output:
!--------
  real,    intent(out):: qn1(is:ie,js:je,npz)  ! Scalar quantity on 3d model grid

! local
  integer :: i, j, k

    do j = js,je
       call mappm(km, pe0(is:ie,j,:), qn0(is:ie,j,:), npz, pe1(is:ie,j,:), qn1(is:ie,j,:), is,ie, n, 8)
    enddo

  end subroutine remap_3d

!###################################################################################
!###################################################################################

end module fv_climate_nudge_mod

