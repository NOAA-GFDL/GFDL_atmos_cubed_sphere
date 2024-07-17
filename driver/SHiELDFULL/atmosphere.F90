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

module atmosphere_mod
#include <fms_platform.h>

!-----------------------------------------------------------------------
!
! Interface for Cubed_Sphere fv dynamical core
!
!-----------------------------------------------------------------------

!-----------------
! FMS modules:
!-----------------
use block_control_mod,      only: block_control_type
use constants_mod,          only: cp_air, rdgas, grav, rvgas, kappa, pstd_mks, pi
use time_manager_mod,       only: time_type, get_time, set_time, operator(+), &
                                  operator(-), operator(/), time_type_to_real
use fms_mod,                only: error_mesg, FATAL,                 &
                                  check_nml_error, stdlog,           &
                                  write_version_number,              &
                                  mpp_clock_id, mpp_clock_begin,     &
                                  mpp_clock_end, CLOCK_SUBCOMPONENT, &
                                  clock_flag_default
use fms2_io_mod,            only: file_exists
use mpp_mod,                only: mpp_error, stdout, FATAL, WARNING, NOTE, &
                                  input_nml_file, mpp_root_pe,    &
                                  mpp_npes, mpp_pe, mpp_chksum,   &
                                  mpp_get_current_pelist,         &
                                  mpp_set_current_pelist, mpp_sync
use mpp_parameter_mod,      only: EUPDATE, WUPDATE, SUPDATE, NUPDATE
use mpp_domains_mod,        only: domain2d, mpp_update_domains
use xgrid_mod,              only: grid_box_type
use field_manager_mod,      only: MODEL_ATMOS
use tracer_manager_mod,     only: get_tracer_index, get_number_tracers, &
                                  NO_TRACER, get_tracer_names
use IPD_typedefs,           only: IPD_data_type, kind_phys
use data_override_mod,      only: data_override_init
use fv_iau_mod,             only: IAU_external_data_type
use atmos_cmip_diag_mod,   only: atmos_cmip_diag_init
use atmos_global_diag_mod, only: atmos_global_diag_init, &
                                 atmos_global_diag_end
!-----------------
! FV core modules:
!-----------------
use fv_arrays_mod,      only: fv_atmos_type, R_GRID, fv_grid_bounds_type, phys_diag_type
use fv_control_mod,     only: fv_control_init, fv_end, ngrids
use fv_eta_mod,         only: get_eta_level
use fv_fill_mod,        only: fill_gfs
use dyn_core_mod,       only: del2_cubed
use fv_dynamics_mod,    only: fv_dynamics
use fv_nesting_mod,     only: twoway_nesting
use fv_diagnostics_mod, only: fv_diag_init, fv_diag, fv_time, prt_maxmin, prt_height
use fv_nggps_diags_mod, only: fv_nggps_diag_init, fv_nggps_diag
use fv_restart_mod,     only: fv_restart, fv_write_restart
use fv_timing_mod,      only: timing_on, timing_off, timing_init, timing_prt
use fv_mp_mod,          only: is_master
use fv_sg_mod,          only: fv_subgrid_z
use fv_update_phys_mod, only: fv_update_phys
use fv_io_mod,          only: fv_io_register_nudge_restart
use fv_nwp_nudge_mod,   only: fv_nwp_nudge_init, fv_nwp_nudge_end, do_adiabatic_init
use fv_regional_mod,    only: start_regional_restart, read_new_bc_data
use fv_regional_mod,    only: a_step, p_step
use fv_regional_mod,    only: current_time_in_seconds
use fv_grid_utils_mod,  only: g_sum

use mpp_domains_mod, only:  mpp_get_data_domain, mpp_get_compute_domain
use gfdl_mp_mod,        only: gfdl_mp_init, gfdl_mp_end
use diag_manager_mod,   only: send_data
use external_aero_mod,  only: load_aero, read_aero, clean_aero
use coarse_graining_mod, only: coarse_graining_init
use coarse_grained_diagnostics_mod, only: fv_coarse_diag_init, fv_coarse_diag
use coarse_grained_restart_files_mod, only: fv_coarse_restart_init

implicit none
private

!--- driver routines
public :: atmosphere_init, atmosphere_end, atmosphere_restart, &
          atmosphere_dynamics, atmosphere_state_update

!--- utility routines
public :: atmosphere_resolution, atmosphere_grid_bdry, &
          atmosphere_grid_ctr, atmosphere_domain, &
          atmosphere_control_data, atmosphere_pref, &
          atmosphere_diag_axes, atmosphere_etalvls, &
          atmosphere_hgt, atmosphere_scalar_field_halo, &
!rab          atmosphere_tracer_postinit, &
          atmosphere_diss_est, & ! dissipation estimate for SKEB
          atmosphere_nggps_diag, &
          get_bottom_mass, get_bottom_wind,   &
          get_stock_pe, set_atmosphere_pelist, &
          atmosphere_coarse_diag_axes

!--- physics/radiation data exchange routines
public :: atmos_phys_driver_statein
public :: atmosphere_coarse_graining_parameters
public :: atmosphere_coarsening_strategy

!-----------------------------------------------------------------------
! version number of this module
! Include variable "version" to be written to log file.
#include<file_version.h>
character(len=20)   :: mod_name = 'SHiELD/atmosphere_mod'

!---- private data ----
  type (time_type) :: Time_step_atmos
  public Atm, mygrid

  !These are convenience variables for local use only, and are set to values in Atm%
  real    :: dt_atmos
  real    :: zvir
  integer :: npx, npy, npz, ncnst, pnats
  integer :: isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed
  integer :: nq                       ! transported tracers
  integer :: sec, seconds, days
  integer :: id_dycore, id_fv_diag, id_update, id_dynam, id_subgrid
  logical :: cold_start = .false.       ! read in initial condition

  integer, dimension(:), allocatable :: id_tracerdt_dyn
  integer :: sphum, liq_wat, rainwat, ice_wat, snowwat, graupel  !condensate species

  integer :: mygrid = 1
  integer :: p_split = 1
  integer, allocatable :: pelist(:)
  logical, allocatable :: grids_on_this_pe(:)
  type(fv_atmos_type), allocatable, target :: Atm(:)

  integer :: id_udt_dyn, id_vdt_dyn

  real, parameter:: w0_big = 60.  ! to prevent negative w-tracer diffusion

!---dynamics tendencies for use in fv_subgrid_z and during fv_update_phys
  real, allocatable, dimension(:,:,:)   :: u_dt, v_dt, t_dt, qv_dt
  real, allocatable :: pref(:,:), dum1d(:), ps_dt(:,:)

  logical :: first_diag = .true.

contains



 subroutine atmosphere_init (Time_init, Time, Time_step, Grid_box, area, IAU_Data)
   type (time_type),    intent(in)    :: Time_init, Time, Time_step
   type(grid_box_type), intent(inout) :: Grid_box
   type(iau_external_data_type), intent(out) :: IAU_Data
   real(kind=kind_phys), pointer, dimension(:,:), intent(inout) :: area
!--- local variables ---
   integer :: i, n
   integer :: itrac
   logical :: do_atmos_nudge
   character(len=32) :: tracer_name, tracer_units
   real :: ps1, ps2

   integer :: nlunit = 9999
   character (len = 64) :: fn_nml = 'input.nml'

   call timing_init
   call timing_on('ATMOS_TOTAL')
   call timing_on('ATMOS_INIT')

   !For regional
   a_step = 0
   current_time_in_seconds = time_type_to_real( Time - Time_init )
   if (mpp_pe() == 0) write(0,"('atmosphere_init: current_time_seconds = ',f9.1)")current_time_in_seconds

   allocate(pelist(mpp_npes()))
   call mpp_get_current_pelist(pelist)

   zvir = rvgas/rdgas - 1.

!---- compute physics/atmos time step in seconds ----

   Time_step_atmos = Time_step
   call get_time (Time_step_atmos, sec)
   dt_atmos = real(sec)

!----- initialize FV dynamical core -----
   !NOTE do we still need the second file_exist call?
   cold_start = (.not.file_exists('INPUT/fv_core.res.nc') .and. .not.file_exists('INPUT/fv_core.res.tile1.nc'))

   call fv_control_init( Atm, dt_atmos, mygrid, grids_on_this_pe, p_split )  ! allocates Atm components; sets mygrid

   if (Atm(mygrid)%coarse_graining%write_coarse_restart_files .or. &
        Atm(mygrid)%coarse_graining%write_coarse_diagnostics) then
      call coarse_graining_init(Atm(mygrid)%flagstruct%npx, Atm(mygrid)%npz, &
           Atm(mygrid)%layout, Atm(mygrid)%bd%is, Atm(mygrid)%bd%ie, &
           Atm(mygrid)%bd%js, Atm(mygrid)%bd%je, Atm(mygrid)%coarse_graining%factor, &
           Atm(mygrid)%coarse_graining%nx_coarse, &
           Atm(mygrid)%coarse_graining%strategy, &
           Atm(mygrid)%coarse_graining%domain)
   endif

   Atm(mygrid)%Time_init = Time_init

!----- write version and namelist to log file -----
   call write_version_number ( mod_name, version )

!-----------------------------------

   npx   = Atm(mygrid)%npx
   npy   = Atm(mygrid)%npy
   npz   = Atm(mygrid)%npz
   ncnst = Atm(mygrid)%ncnst
   pnats = Atm(mygrid)%flagstruct%pnats

   isc = Atm(mygrid)%bd%isc
   iec = Atm(mygrid)%bd%iec
   jsc = Atm(mygrid)%bd%jsc
   jec = Atm(mygrid)%bd%jec

   isd = isc - Atm(mygrid)%bd%ng
   ied = iec + Atm(mygrid)%bd%ng
   jsd = jsc - Atm(mygrid)%bd%ng
   jed = jec + Atm(mygrid)%bd%ng

   nq = ncnst-pnats
   sphum   = get_tracer_index (MODEL_ATMOS, 'sphum' )
   liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat' )
   ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat' )
   rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat' )
   snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat' )
   graupel = get_tracer_index (MODEL_ATMOS, 'graupel' )

   if (max(sphum,liq_wat,ice_wat,rainwat,snowwat,graupel) > Atm(mygrid)%flagstruct%nwat) then
      call mpp_error (FATAL,' atmosphere_init: condensate species are not first in the list of &
                            &tracers defined in the field_table')
   endif

   ! Allocate grid variables to be used to calculate gradient in 2nd order flux exchange
   ! This data is only needed for the COARSEST grid.
   !call switch_current_Atm(Atm(mygrid))

   allocate(Grid_box%dx    (   isc:iec  , jsc:jec+1))
   allocate(Grid_box%dy    (   isc:iec+1, jsc:jec  ))
   allocate(Grid_box%area  (   isc:iec  , jsc:jec  ))
   allocate(Grid_box%edge_w(              jsc:jec+1))
   allocate(Grid_box%edge_e(              jsc:jec+1))
   allocate(Grid_box%edge_s(   isc:iec+1           ))
   allocate(Grid_box%edge_n(   isc:iec+1           ))
   allocate(Grid_box%en1   (3, isc:iec  , jsc:jec+1))
   allocate(Grid_box%en2   (3, isc:iec+1, jsc:jec  ))
   allocate(Grid_box%vlon  (3, isc:iec  , jsc:jec  ))
   allocate(Grid_box%vlat  (3, isc:iec  , jsc:jec  ))
   Grid_box%dx    (   isc:iec  , jsc:jec+1) = Atm(mygrid)%gridstruct%dx    (   isc:iec,   jsc:jec+1)
   Grid_box%dy    (   isc:iec+1, jsc:jec  ) = Atm(mygrid)%gridstruct%dy    (   isc:iec+1, jsc:jec  )
   Grid_box%area  (   isc:iec  , jsc:jec  ) = Atm(mygrid)%gridstruct%area  (   isc:iec  , jsc:jec  )
   Grid_box%edge_w(              jsc:jec+1) = Atm(mygrid)%gridstruct%edge_w(              jsc:jec+1)
   Grid_box%edge_e(              jsc:jec+1) = Atm(mygrid)%gridstruct%edge_e(              jsc:jec+1)
   Grid_box%edge_s(   isc:iec+1           ) = Atm(mygrid)%gridstruct%edge_s(   isc:iec+1)
   Grid_box%edge_n(   isc:iec+1           ) = Atm(mygrid)%gridstruct%edge_n(   isc:iec+1)
   Grid_box%en1   (:, isc:iec  , jsc:jec+1) = Atm(mygrid)%gridstruct%en1   (:, isc:iec  , jsc:jec+1)
   Grid_box%en2   (:, isc:iec+1, jsc:jec  ) = Atm(mygrid)%gridstruct%en2   (:, isc:iec+1, jsc:jec  )
   do i = 1,3
     Grid_box%vlon  (i, isc:iec  , jsc:jec  ) = Atm(mygrid)%gridstruct%vlon  (isc:iec ,  jsc:jec, i )
     Grid_box%vlat  (i, isc:iec  , jsc:jec  ) = Atm(mygrid)%gridstruct%vlat  (isc:iec ,  jsc:jec, i )
   enddo
   allocate (area(isc:iec  , jsc:jec  ))
   area(isc:iec,jsc:jec) = Atm(mygrid)%gridstruct%area_64(isc:iec,jsc:jec)

!----- allocate and zero out the dynamics (and accumulated) tendencies
   allocate( u_dt(isd:ied,jsd:jed,npz), &
             v_dt(isd:ied,jsd:jed,npz), &
             t_dt(isc:iec,jsc:jec,npz), &
             qv_dt(isc:iec,jsc:jec,npz), &
             ps_dt(isd:ied,jsd:jed) )

!--- allocate pref
   allocate(pref(npz+1,2), dum1d(npz+1))

   call gfdl_mp_init(input_nml_file, stdlog(), Atm(mygrid)%flagstruct%hydrostatic)

   call timing_on('FV_RESTART')
   call fv_restart(Atm(mygrid)%domain, Atm, dt_atmos, seconds, days, cold_start, &
                   Atm(mygrid)%gridstruct%grid_type, mygrid)
   call timing_off('FV_RESTART')

   fv_time = Time

!----- initialize atmos_axes and fv_dynamics diagnostics
       !I've had trouble getting this to work with multiple grids at a time; worth revisiting?
   call fv_diag_init(Atm(mygrid:mygrid), Atm(mygrid)%atmos_axes, Time, npx, npy, npz, Atm(mygrid)%flagstruct%p_ref)

   if (Atm(mygrid)%flagstruct%do_aerosol) then
     call load_aero(Atm(mygrid), Time)
     call read_aero(isc, iec, jsc, jec, npz, nq, Time, Atm(mygrid)%pe(isc:iec,:,jsc:jec), &
       Atm(mygrid)%peln(isc:iec,:,jsc:jec), Atm(mygrid)%q(isc:iec,jsc:jec,:,:), &
       Atm(mygrid)%flagstruct%kord_tr, Atm(mygrid)%flagstruct%fill)
   endif

   if (Atm(mygrid)%coarse_graining%write_coarse_diagnostics) then
      call fv_coarse_diag_init(Atm, Time, Atm(mygrid)%atmos_axes(3), &
           Atm(mygrid)%atmos_axes(4), Atm(mygrid)%coarse_graining)
   endif
   if (Atm(mygrid)%coarse_graining%write_coarse_restart_files) then
      call fv_coarse_restart_init(Atm(mygrid)%npz, Atm(mygrid)%flagstruct%nt_prog, &
           Atm(mygrid)%flagstruct%nt_phys, Atm(mygrid)%flagstruct%hydrostatic, &
           Atm(mygrid)%flagstruct%hybrid_z, Atm(mygrid)%flagstruct%fv_land, &
           Atm(mygrid)%coarse_graining%write_coarse_dgrid_vel_rst, &
           Atm(mygrid)%coarse_graining%write_coarse_agrid_vel_rst, &
           Atm(mygrid)%coarse_graining%restart)
   endif

!---------- reference profile -----------
    ps1 = 101325.
    ps2 =  81060.
    pref(npz+1,1) = ps1
    pref(npz+1,2) = ps2
    call get_eta_level ( npz, ps1, pref(1,1), dum1d, Atm(mygrid)%ak, Atm(mygrid)%bk )
    call get_eta_level ( npz, ps2, pref(1,2), dum1d, Atm(mygrid)%ak, Atm(mygrid)%bk )

!  --- initialize clocks for dynamics, physics_down and physics_up
   id_dycore  = mpp_clock_id ('---FV Dycore',flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
   id_update  = mpp_clock_id ('---FV Update',flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
   id_fv_diag = mpp_clock_id ('---FV Diag',  flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
   id_dynam   = mpp_clock_id ('----FV Dynamics',flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
   id_subgrid = mpp_clock_id ('----FV Subgrid_z',flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )

!---- initialize cmip diagnostic output ----
   call atmos_cmip_diag_init   ( Atm(mygrid)%ak, Atm(mygrid)%bk, pref(1,1), Atm(mygrid)%atmos_axes, Time )
if ( is_master() ) write(*,*) 'CALL atmos_global_diag_init'
   call atmos_global_diag_init ( Atm(mygrid)%atmos_axes, Atm(mygrid)%gridstruct%area(isc:iec,jsc:jec) )

!  --- initiate the start for a restarted regional forecast
   if ( Atm(mygrid)%gridstruct%regional .and. Atm(mygrid)%flagstruct%warm_start ) then

     call start_regional_restart(Atm(1),       &
                                 isc, iec, jsc, jec, &
                                 isd, ied, jsd, jed )
   endif


   if ( Atm(mygrid)%flagstruct%nudge ) then
      call fv_nwp_nudge_init( Time, Atm(mygrid)%atmos_axes, npz, zvir, Atm(mygrid)%ak, Atm(mygrid)%bk, Atm(mygrid)%ts, &
           Atm(mygrid)%phis, Atm(mygrid)%gridstruct, Atm(mygrid)%ks, Atm(mygrid)%npx, Atm(mygrid)%neststruct, Atm(mygrid)%bd)
      call mpp_error(NOTE, 'NWP nudging is active')
   endif
   call fv_io_register_nudge_restart ( Atm )

   if ( Atm(mygrid)%flagstruct%na_init>0 ) then
      if ( .not. Atm(mygrid)%flagstruct%hydrostatic ) then
           call prt_maxmin('Before adi: W', Atm(mygrid)%w, isc, iec, jsc, jec, Atm(mygrid)%ng, npz, 1.)
      endif
      call adiabatic_init(zvir,Atm(mygrid)%flagstruct%nudge_dz)
      if ( .not. Atm(mygrid)%flagstruct%hydrostatic ) then
           call prt_maxmin('After adi: W', Atm(mygrid)%w, isc, iec, jsc, jec, Atm(mygrid)%ng, npz, 1.)
! Not nested?
           call prt_height('na_ini Z500', isc,iec, jsc,jec, 3, npz, 500.E2, Atm(mygrid)%phis, Atm(mygrid)%delz,    &
                Atm(mygrid)%peln, Atm(mygrid)%gridstruct%area_64(isc:iec,jsc:jec), Atm(mygrid)%gridstruct%agrid_64(isc:iec,jsc:jec,2))
      endif
   else
      call mpp_error(NOTE,'No adiabatic initialization correction in use')
   endif

#ifdef DEBUG
   call fv_diag(Atm(mygrid:mygrid), zvir, Time, -1)
   if (Atm(mygrid)%coarse_graining%write_coarse_diagnostics) then
      call fv_coarse_diag(Atm(mygrid:mygrid), fv_time, zvir)
   endif
#endif

   if ( trim(Atm(mygrid)%flagstruct%grid_file) .NE. "Inline" .and. trim(Atm(mygrid)%flagstruct%grid_file) .NE. "" &
      & .and. .NOT.Atm(mygrid)%gridstruct%bounded_domain ) then
      call data_override_init(Atm_domain_in = Atm(mygrid)%domain)
   endif

   call timing_off('ATMOS_INIT')

 end subroutine atmosphere_init


 subroutine p_adi(km, ng, ifirst, ilast, jfirst, jlast, ptop,   &
                  delp, pt, ps, pe, peln, pk, pkz, hydrostatic)
! Given (ptop, delp) computes (ps, pk, pe, peln, pkz)
! Input:
   integer,  intent(in):: km, ng
   integer,  intent(in):: ifirst, ilast            ! Longitude strip
   integer,  intent(in):: jfirst, jlast            ! Latitude strip
   logical, intent(in)::  hydrostatic
   real, intent(in):: ptop
   real, intent(in)::   pt(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng, km)
   real, intent(in):: delp(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng, km)
! Output:
   real, intent(out) ::   ps(ifirst-ng:ilast+ng, jfirst-ng:jlast+ng)
   real, intent(out) ::   pk(ifirst:ilast, jfirst:jlast, km+1)
   real, intent(out) ::   pe(ifirst-1:ilast+1,km+1,jfirst-1:jlast+1) ! Ghosted Edge pressure
   real, intent(out) :: peln(ifirst:ilast, km+1, jfirst:jlast)    ! Edge pressure
   real, intent(out) ::  pkz(ifirst:ilast, jfirst:jlast, km)
! Local
   real pek
   integer i, j, k

   pek = ptop ** kappa
!$OMP parallel do default (none) &
!$OMP              shared (ifirst,ilast,jfirst,jlast,km,ptop,pek,pe,pk, &
!$OMP                      ps,delp,peln,hydrostatic,pkz) &
!$OMP             private (j, i, k)
   do j=jfirst,jlast
      do i=ifirst,ilast
         pe(i,1,j) = ptop
         pk(i,j,1) = pek
      enddo

      do k=2,km+1
         do i=ifirst,ilast
            pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
            peln(i,k,j) = log(pe(i,k,j))
            pk(i,j,k) = exp( kappa*peln(i,k,j) )
         enddo
      enddo

      do i=ifirst,ilast
         ps(i,j) = pe(i,km+1,j)
      enddo

      if ( hydrostatic ) then
         do k=1,km
            do i=ifirst,ilast
               pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k))/(kappa*(peln(i,k+1,j)-peln(i,k,j)))
            enddo
         enddo
      endif
   enddo

 end subroutine p_adi


 subroutine atmosphere_dynamics ( Time )
   type(time_type),intent(in) :: Time
   integer :: itrac, n, psc
   integer :: k, w_diff, nt_dyn
   logical :: used

   type(time_type) :: atmos_time
   integer :: atmos_time_step
   real :: rdt
   real :: time_total
!---- Call FV dynamics -----

   call timing_on('ATMOS_DYNAMICS')

   call mpp_clock_begin (id_dycore)
   call mpp_clock_begin (id_dynam)

   n = mygrid
   a_step = a_step + 1
!
!*** If this is a regional run then read in the next boundary data when it is time.
!
   if(Atm(n)%flagstruct%regional)then

     call read_new_bc_data(Atm(n), Time, Time_step_atmos, p_split, &
                           isd, ied, jsd, jed )
   endif

   if (Atm(mygrid)%flagstruct%do_aerosol) then
     call read_aero(isc, iec, jsc, jec, npz, nq, Time, Atm(mygrid)%pe(isc:iec,:,jsc:jec), &
	     Atm(mygrid)%peln(isc:iec,:,jsc:jec), Atm(mygrid)%q(isc:iec,jsc:jec,:,:), &
         Atm(mygrid)%flagstruct%kord_tr, Atm(mygrid)%flagstruct%fill)
   endif

!save ps to ps_dt before dynamics update
   ps_dt(:,:)=Atm(n)%ps(:,:)

   do psc=1,abs(p_split)
      p_step = psc

     call timing_on('FV_DYNAMICS')
!uc/vc only need be same on coarse grid? However BCs do need to be the same
     call fv_dynamics(npx, npy, npz, nq, Atm(n)%ng, dt_atmos/real(abs(p_split)),&
                      Atm(n)%flagstruct%consv_te, Atm(n)%flagstruct%fill,  &
                      Atm(n)%flagstruct%reproduce_sum, kappa, cp_air, zvir,&
                      Atm(n)%ptop, Atm(n)%ks, nq,                          &
                      Atm(n)%flagstruct%n_split, Atm(n)%flagstruct%q_split,&
                      Atm(n)%u0, Atm(n)%v0, Atm(n)%u, Atm(n)%v, Atm(n)%w,  &
                      Atm(n)%delz, Atm(n)%flagstruct%hydrostatic,          &
                      Atm(n)%pt, Atm(n)%delp, Atm(n)%q, Atm(n)%ps,         &
                      Atm(n)%pe, Atm(n)%pk, Atm(n)%peln,                   &
                      Atm(n)%pkz, Atm(n)%phis, Atm(n)%q_con,               &
                      Atm(n)%omga, Atm(n)%ua, Atm(n)%va, Atm(n)%uc,        &
                      Atm(n)%vc, Atm(n)%ak, Atm(n)%bk, Atm(n)%mfx,         &
                      Atm(n)%mfy, Atm(n)%cx, Atm(n)%cy, Atm(n)%ze0,        &
                      Atm(n)%flagstruct%hybrid_z,                          &
                      Atm(n)%gridstruct, Atm(n)%flagstruct,                &
                      Atm(n)%neststruct, Atm(n)%idiag, Atm(n)%bd,          &
                      Atm(n)%parent_grid, Atm(n)%domain, Atm(n)%inline_mp, &
                      Atm(n)%diss_est,time_total=time_total)

     call timing_off('FV_DYNAMICS')

    if (ngrids > 1 .and. (psc < p_split .or. p_split < 0)) then
       call mpp_sync()
       call timing_on('TWOWAY_UPDATE')
       call twoway_nesting(Atm, ngrids, grids_on_this_pe, zvir, fv_time, mygrid)
       call timing_off('TWOWAY_UPDATE')
    endif

    end do !p_split

    if (.not. Atm(n)%flagstruct%hydrostatic .and. .not. Atm(n)%flagstruct%pass_full_omega_to_physics_in_non_hydrostatic_mode) then
        Atm(n)%local_omga(isc:iec,jsc:jec,1:npz) = Atm(n)%delp(isc:iec,jsc:jec,1:npz) / Atm(n)%delz(isc:iec,jsc:jec,1:npz) * Atm(n)%w(isc:iec,jsc:jec,1:npz)
        if(Atm(n)%flagstruct%nf_omega>0)   then
           call del2_cubed(&
                Atm(n)%local_omga, &
                0.18*Atm(n)%gridstruct%da_min, &
                Atm(n)%gridstruct, &
                Atm(n)%domain, &
                Atm(n)%npx, &
                Atm(n)%npy, &
                Atm(n)%npz, &
                Atm(n)%flagstruct%nf_omega, &
                Atm(n)%bd)
        endif
    endif

   call mpp_clock_end (id_dynam)
   call mpp_clock_begin (id_subgrid)

!-----------------------------------------------------
!--- COMPUTE SUBGRID Z
!-----------------------------------------------------
!--- zero out tendencies

    call timing_on('FV_SUBGRID_Z')

    u_dt(:,:,:)   = 0. ! These are updated by fv_subgrid_z
    v_dt(:,:,:)   = 0.
! t_dt is used for two different purposes:
!    1 - to calculate the diagnostic temperature tendency from fv_subgrid_z
!    2 - as an accumulator for the IAU increment and physics tendency
! because of this, it will need to be zeroed out after the diagnostic is calculated
    t_dt(:,:,:)   = Atm(n)%pt(isc:iec,jsc:jec,:)
    qv_dt(:,:,:)  = Atm(n)%q (isc:iec,jsc:jec,:,sphum)

    rdt = 1./dt_atmos

    w_diff = get_tracer_index (MODEL_ATMOS, 'w_diff' )
    if ( Atm(n)%flagstruct%fv_sg_adj > 0 ) then
      nt_dyn = nq
      if ( w_diff /= NO_TRACER ) then
        nt_dyn = nq - 1
      endif
      call fv_subgrid_z(isd, ied, jsd, jed, isc, iec, jsc, jec, Atm(n)%npz, &
                        nt_dyn, dt_atmos, Atm(n)%flagstruct%fv_sg_adj,      &
                        Atm(n)%flagstruct%nwat, Atm(n)%delp, Atm(n)%pe,     &
                        Atm(n)%peln, Atm(n)%pkz, Atm(n)%pt, Atm(n)%q,       &
                        Atm(n)%ua, Atm(n)%va, Atm(n)%flagstruct%hydrostatic,&
                        Atm(n)%w, Atm(n)%delz, u_dt, v_dt, t_dt, Atm(n)%flagstruct%n_sponge)
    endif

#ifdef USE_Q_DT
    if ( .not. Atm(n)%flagstruct%hydrostatic .and. w_diff /= NO_TRACER ) then
!$OMP parallel do default (none) &
!$OMP              shared (isc, iec, jsc, jec, w_diff, n, Atm, q_dt) &
!$OMP             private (k)
       do k=1, Atm(n)%npz
          Atm(n)%q(isc:iec,jsc:jec,k,w_diff) = Atm(n)%w(isc:iec,jsc:jec,k) + w0_big
          q_dt(:,:,k,w_diff) = 0.
        enddo
    endif
#endif

    if (allocated(Atm(n)%sg_diag%u_dt)) then
       Atm(n)%sg_diag%u_dt = u_dt(isc:iec,jsc:jec,:)
    endif
    if (allocated(Atm(n)%sg_diag%v_dt)) then
       Atm(n)%sg_diag%v_dt = v_dt(isc:iec,jsc:jec,:)
    endif
    if (allocated(Atm(n)%sg_diag%t_dt)) then
       t_dt(:,:,:) = rdt*(Atm(1)%pt(isc:iec,jsc:jec,:) - t_dt(:,:,:))
       Atm(n)%sg_diag%t_dt = t_dt(isc:iec,jsc:jec,:)
    endif
    if (allocated(Atm(n)%sg_diag%qv_dt)) then
       qv_dt(:,:,:) = rdt*(Atm(1)%q(isc:iec,jsc:jec,:,sphum) - qv_dt(:,:,:))
       Atm(n)%sg_diag%qv_dt = qv_dt(isc:iec,jsc:jec,:)
    endif

! zero out t_dt for use as an accumulator
   t_dt = 0.

   call timing_off('FV_SUBGRID_Z')

   call mpp_clock_end (id_subgrid)
   call mpp_clock_end (id_dycore)

   call timing_off('ATMOS_DYNAMICS')

 end subroutine atmosphere_dynamics


 subroutine atmosphere_end (Time, Grid_box )!rab, Radiation, Physics)
   type (time_type),      intent(in)    :: Time
   type(grid_box_type),   intent(inout) :: Grid_box
!rab   type (radiation_type), intent(inout) :: Radiation
!rab   type (physics_type),   intent(inout) :: Physics

   call timing_on('ATMOS_END')

  ! initialize domains for writing global physics data
  if ( Atm(mygrid)%flagstruct%nudge ) call fv_nwp_nudge_end

   if (Atm(mygrid)%flagstruct%do_aerosol) then
     call clean_aero()
   endif

   call gfdl_mp_end ( )

   if (first_diag) then
      call timing_on('FV_DIAG')
      call fv_diag(Atm(mygrid:mygrid), zvir, fv_time, Atm(mygrid)%flagstruct%print_freq)
      call fv_nggps_diag(Atm(mygrid:mygrid), zvir, fv_time)
      if (Atm(mygrid)%coarse_graining%write_coarse_diagnostics) then
         call fv_coarse_diag(Atm(mygrid:mygrid), fv_time, zvir)
      endif
      first_diag = .false.
      call timing_off('FV_DIAG')
   endif

   call timing_on('FV_DIAG')
   call atmos_global_diag_end
   call fv_end(Atm, mygrid)
   call timing_off('FV_DIAG')
   deallocate (Atm)

   deallocate( u_dt, v_dt, t_dt, qv_dt, ps_dt, pref, dum1d )

   call timing_off('ATMOS_END')
   call timing_off('ATMOS_TOTAL')
   call timing_prt( mpp_pe() )

 end subroutine atmosphere_end



  !#######################################################################
  ! <SUBROUTINE NAME="atmosphere_restart">
  ! <DESCRIPTION>
  !  Write out restart files registered through register_restart_file
  ! </DESCRIPTION>
  subroutine atmosphere_restart(timestamp)
    character(len=*),  intent(in) :: timestamp

    call fv_write_restart(Atm(mygrid), timestamp)

  end subroutine atmosphere_restart
  ! </SUBROUTINE>


 subroutine atmosphere_resolution (i_size, j_size, global)
   integer, intent(out)          :: i_size, j_size
   logical, intent(in), optional :: global
   logical :: local

   local = .true.
   if( PRESENT(global) ) local = .NOT.global

   if( local ) then
       i_size = iec - isc + 1
       j_size = jec - jsc + 1
   else
       i_size = npx - 1
       j_size = npy - 1
   end if

 end subroutine atmosphere_resolution


 subroutine atmosphere_pref (p_ref)
   real, dimension(:,:), intent(inout) :: p_ref

   p_ref = pref

 end subroutine atmosphere_pref


 subroutine atmosphere_control_data (i1, i2, j1, j2, kt, p_hydro, hydro, tile_num, &
                                     do_inline_mp, do_cosp)
   integer, intent(out)           :: i1, i2, j1, j2, kt
   logical, intent(out), optional :: p_hydro, hydro
   integer, intent(out), optional :: tile_num
   logical, intent(out), optional :: do_inline_mp, do_cosp
   i1 = Atm(mygrid)%bd%isc
   i2 = Atm(mygrid)%bd%iec
   j1 = Atm(mygrid)%bd%jsc
   j2 = Atm(mygrid)%bd%jec
   kt = Atm(mygrid)%npz

   if (present(p_hydro)) p_hydro = Atm(mygrid)%flagstruct%phys_hydrostatic
   if (present(  hydro))   hydro = Atm(mygrid)%flagstruct%hydrostatic
   if (present(tile_num)) tile_num = Atm(mygrid)%global_tile
   if (present(do_inline_mp)) do_inline_mp = Atm(mygrid)%flagstruct%do_inline_mp
   if (present(do_cosp)) do_cosp = Atm(mygrid)%flagstruct%do_cosp

 end subroutine atmosphere_control_data


 subroutine atmosphere_grid_ctr (lon, lat)
!---------------------------------------------------------------
!    returns the longitude and latitude cell centers
!---------------------------------------------------------------
    real(kind=kind_phys), intent(out) :: lon(:,:), lat(:,:)   ! Unit: radian
! Local data:
    integer i,j

    do j=jsc,jec
       do i=isc,iec
          lon(i-isc+1,j-jsc+1) = Atm(mygrid)%gridstruct%agrid_64(i,j,1)
          lat(i-isc+1,j-jsc+1) = Atm(mygrid)%gridstruct%agrid_64(i,j,2)
       enddo
    end do

 end subroutine atmosphere_grid_ctr


 subroutine atmosphere_grid_bdry (blon, blat, global)
!---------------------------------------------------------------
!    returns the longitude and latitude grid box edges
!    for either the local PEs grid (default) or the global grid
!---------------------------------------------------------------
    real,    intent(out) :: blon(:,:), blat(:,:)   ! Unit: radian
    logical, intent(in), optional :: global
! Local data:
    integer i,j

    if( PRESENT(global) ) then
      if (global) call mpp_error(FATAL, '==> global grid is no longer available &
                               & in the Cubed Sphere')
    endif

    do j=jsc,jec+1
       do i=isc,iec+1
          blon(i-isc+1,j-jsc+1) = Atm(mygrid)%gridstruct%grid(i,j,1)
          blat(i-isc+1,j-jsc+1) = Atm(mygrid)%gridstruct%grid(i,j,2)
       enddo
    end do

 end subroutine atmosphere_grid_bdry


 subroutine set_atmosphere_pelist ()
   call mpp_set_current_pelist(Atm(mygrid)%pelist, no_sync=.TRUE.)
 end subroutine set_atmosphere_pelist


 subroutine atmosphere_domain ( fv_domain, rd_domain, layout, regional, bounded_domain )
   type(domain2d), intent(out) :: fv_domain, rd_domain
   integer, intent(out) :: layout(2)
   logical, intent(out) :: regional
   logical, intent(out) :: bounded_domain
!  returns the domain2d variable associated with the coupling grid
!  note: coupling is done using the mass/temperature grid with no halos

   fv_domain = Atm(mygrid)%domain_for_coupler
   rd_domain = Atm(mygrid)%domain_for_read
   layout(1:2) =  Atm(mygrid)%layout(1:2)

   regional = Atm(mygrid)%flagstruct%regional
   bounded_domain = Atm(mygrid)%gridstruct%bounded_domain

 end subroutine atmosphere_domain


 subroutine atmosphere_diag_axes ( axes )
   integer, intent(out) :: axes (:)

!----- returns the axis indices for the atmospheric (mass) grid -----
   if ( size(axes(:)) < 0 .or. size(axes(:)) > 4 ) call error_mesg (    &
                               'get_atmosphere_axes in atmosphere_mod', &
                               'size of argument is incorrect', FATAL   )

   axes (1:size(axes(:))) = Atm(mygrid)%atmos_axes (1:size(axes(:)))

 end subroutine atmosphere_diag_axes

 !>@brief The subroutine 'atmosphere_coarse_diag_axes' is an API to return the axis indices
 !! for the coarse atmospheric (mass) grid.
  subroutine atmosphere_coarse_diag_axes(coarse_axes)
    integer, intent(out) :: coarse_axes(4)

    coarse_axes = (/ &
         Atm(mygrid)%coarse_graining%id_xt_coarse, &
         Atm(mygrid)%coarse_graining%id_yt_coarse, &
         Atm(mygrid)%coarse_graining%id_pfull, &
         Atm(mygrid)%coarse_graining%id_phalf /)
  end subroutine atmosphere_coarse_diag_axes

 subroutine atmosphere_etalvls (ak, bk, flip)
   real(kind=kind_phys), pointer, dimension(:), intent(inout) :: ak, bk
   logical, intent(in) :: flip

   allocate(ak(npz+1))
   allocate(bk(npz+1))

   if (flip) then
     ak(1:npz+1) = Atm(mygrid)%ak(npz+1:1:-1)
     bk(1:npz+1) = Atm(mygrid)%bk(npz+1:1:-1)
   else
     ak(1:npz+1) = Atm(mygrid)%ak(1:npz+1)
     bk(1:npz+1) = Atm(mygrid)%bk(1:npz+1)
   endif
 end subroutine atmosphere_etalvls


 subroutine atmosphere_hgt (hgt, position, relative, flip)
   real(kind=kind_phys), pointer, dimension(:,:,:), intent(inout) :: hgt
   character(len=5), intent(in) :: position
   logical, intent(in) :: relative
   logical, intent(in) :: flip
   !--- local variables ---
   integer:: lev, k, j, i
   real(kind=kind_phys), allocatable, dimension(:,:,:) :: z, dz

   if ((position .ne. "layer") .and. (position .ne. "level")) then
     call mpp_error (FATAL, 'atmosphere_hgt:: incorrect position specification')
   endif

   allocate(z(iec-isc+1,jec-jsc+1,npz+1))
   allocate(dz(iec-isc+1,jec-jsc+1,npz))
   z  = 0
   dz = 0

   if (Atm(mygrid)%flagstruct%hydrostatic) then
     !--- generate dz using hydrostatic assumption
     do j = jsc, jec
       do i = isc, iec
         dz(i-isc+1,j-jsc+1,1:npz) = (rdgas/grav)*Atm(mygrid)%pt(i,j,1:npz)  &
                         * (Atm(mygrid)%peln(i,1:npz,j) - Atm(mygrid)%peln(i,2:npz+1,j))
       enddo
     enddo
   else
     !--- use non-hydrostatic delz directly
     do j = jsc, jec
       do i = isc, iec
         dz(i-isc+1,j-jsc+1,1:npz) = Atm(mygrid)%delz(i,j,1:npz)
       enddo
     enddo
   endif

   !--- calculate geometric heights at the interfaces (levels)
   !--- if needed, flip the indexing during this step
   if (flip) then
     if (.not. relative) then
       do j = jsc, jec
       do i = isc, iec
         z(i-isc+1,j-jsc+1,1) = Atm(mygrid)%phis(i,j)/grav
       enddo
       enddo
     endif
     do k = 2,npz+1
       z(:,:,k) = z(:,:,k-1) - dz(:,:,npz+2-k)
     enddo
   else
     if (.not. relative) then
       do j = jsc, jec
       do i = isc, iec
         z(i-isc+1,j-jsc+1,npz+1) = Atm(mygrid)%phis(i,j)/grav
       enddo
       enddo
     endif
     do k = npz,1,-1
       z(:,:,k) = z(:,:,k+1) - dz(:,:,k)
     enddo
   endif

   !--- allocate and set either the level or layer height for return
   if (position == "level") then
     allocate (hgt(iec-isc+1,jec-jsc+1,npz+1))
     hgt = z
   elseif (position == "layer") then
     allocate (hgt(iec-isc+1,jec-jsc+1,npz))
     hgt(:,:,1:npz) = 0.5d0 * (z(:,:,1:npz) + z(:,:,2:npz+1))
   endif

   deallocate (z)
   deallocate (dz)

 end subroutine atmosphere_hgt


 subroutine atmosphere_scalar_field_halo (data, halo, isize, jsize, ksize, data_p)
   !--------------------------------------------------------------------
   ! data   - output array to return the field with halo (i,j,k)
   !          optionally input for field already in (i,j,k) form
   !          sized to include the halo of the field (+ 2*halo)
   ! halo   - size of the halo (must be less than 3)
   ! ied    - horizontal resolution in i-dir with haloes
   ! jed    - horizontal resolution in j-dir with haloes
   ! ksize  - vertical resolution
   ! data_p - optional input field in packed format (ix,k)
   !--------------------------------------------------------------------
   !--- interface variables ---
   real(kind=kind_phys), dimension(1:isize,1:jsize,ksize), intent(inout) :: data
   integer, intent(in) :: halo
   integer, intent(in) :: isize
   integer, intent(in) :: jsize
   integer, intent(in) :: ksize
   real(kind=kind_phys), dimension(:,:), optional, intent(in) :: data_p
   !--- local variables ---
   integer :: i, j, k
   integer :: ic, jc
   character(len=44) :: modname = 'atmosphere_mod::atmosphere_scalar_field_halo'
   integer :: mpp_flags

   !--- perform error checking
   if (halo .gt. 3) call mpp_error(FATAL, modname//' - halo.gt.3 requires extending the MPP domain to support')
   ic = isize - 2 * halo
   jc = jsize - 2 * halo

   !--- if packed data is present, unpack it into the two-dimensional data array
   if (present(data_p)) then
     if (ic*jc .ne. size(data_p,1)) call mpp_error(FATAL, modname//' - incorrect sizes for incoming &
                                                  &variables data and data_p')
     data = 0.
!$OMP parallel do default (none) &
!$OMP              shared (data, data_p, halo, ic, jc, ksize) &
!$OMP             private (i, j, k)
     do k = 1, ksize
       do j = 1, jc
         do i = 1, ic
           data(i+halo, j+halo, k) = data_p(i + (j-1)*ic, k)
         enddo
       enddo
     enddo
   endif

   mpp_flags = EUPDATE + WUPDATE + SUPDATE + NUPDATE
   if (halo == 1) then
     call mpp_update_domains(data, Atm(mygrid)%domain_for_coupler, flags=mpp_flags, complete=.true.)
   elseif (halo == 3) then
     call mpp_update_domains(data, Atm(mygrid)%domain, flags=mpp_flags, complete=.true.)
   else
     call mpp_error(FATAL, modname//' - unsupported halo size')
   endif

   !--- fill the halo points when at a corner of the cubed-sphere tile
   !--- interior domain corners are handled correctly
   if ( (isc==1) .or. (jsc==1) .or. (iec==npx-1) .or. (jec==npy-1) ) then
     do k = 1, ksize
       do j=1,halo
         do i=1,halo
           if ((isc==    1) .and. (jsc==    1)) data(halo+1-j ,halo+1-i ,k) = data(halo+i     ,halo+1-j ,k)  !SW Corner
           if ((isc==    1) .and. (jec==npy-1)) data(halo+1-j ,halo+jc+i,k) = data(halo+i     ,halo+jc+j,k)  !NW Corner
           if ((iec==npx-1) .and. (jsc==    1)) data(halo+ic+j,halo+1-i ,k) = data(halo+ic-i+1,halo+1-j ,k)  !SE Corner
           if ((iec==npx-1) .and. (jec==npy-1)) data(halo+ic+j,halo+jc+i,k) = data(halo+ic-i+1,halo+jc+j,k)  !NE Corner
         enddo
       enddo
     enddo
   endif

   return
 end subroutine atmosphere_scalar_field_halo

 subroutine atmosphere_diss_est (npass)
   use dyn_core_mod, only: del2_cubed
   !--- interface variables ---
   integer, intent(in) :: npass
   !--- local variables
   integer:: k

   !horizontally smooth dissiapation estimate for SKEB
   ! 3 passes before taking absolute value
   do k = 1,min(3,npass)
     call del2_cubed(Atm(mygrid)%diss_est, 0.25*Atm(mygrid)%gridstruct%da_min, Atm(mygrid)%gridstruct, &
                     Atm(mygrid)%domain, npx, npy, npz, 3, Atm(mygrid)%bd)
   enddo

   Atm(mygrid)%diss_est=abs(Atm(mygrid)%diss_est)

   do k = 4,npass
     call del2_cubed(Atm(mygrid)%diss_est, 0.25*Atm(mygrid)%gridstruct%da_min, Atm(mygrid)%gridstruct, &
                     Atm(mygrid)%domain, npx, npy, npz, 3, Atm(mygrid)%bd)
   enddo
   ! provide back sqrt of dissipation estimate
   Atm(mygrid)%diss_est=sqrt(abs(Atm(mygrid)%diss_est))

 end subroutine atmosphere_diss_est

 subroutine atmosphere_nggps_diag (Time, init)
   !----------------------------------------------
   ! api for output of NCEP/EMC diagnostics
   !
   ! if register is present and set to .true.
   ! will make the initialization call
   !
   ! outputs 3D state fields via either
   ! NCEP write_component or GFDL/FMS diag_manager
   !----------------------------------------------
   type(time_type),   intent(in) :: Time
   logical, optional, intent(in) :: init

   if (PRESENT(init)) then
     if (init .eqv. .true.) then
       call fv_nggps_diag_init(Atm(mygrid:mygrid), Atm(mygrid)%atmos_axes, Time)
       return
     else
       call mpp_error(FATAL, 'atmosphere_nggps_diag - calling with init present, but set to .false.')
     endif
   endif

   call fv_nggps_diag(Atm(mygrid:mygrid), zvir, Time)

 end subroutine atmosphere_nggps_diag


!--- Need to know the formulation of the mixing ratio being imported into FV3
!--- in order to adjust it in a consistent manner for advection
!rab subroutine atmosphere_tracer_postinit (IPD_Data, Atm_block)
!rab   !--- interface variables ---
!rab   type(IPD_data_type),       intent(in) :: IPD_Data(:)
!rab   type(block_control_type),  intent(in) :: Atm_block
!rab   !--- local variables ---
!rab   integer :: i, j, ix, k, k1, n, nwat, nb, blen
!rab   real(kind=kind_phys) :: qwat
!rab
!rab   if( nq<3 ) call mpp_error(FATAL, 'GFS phys must have 3 interactive tracers')
!rab
!rab   n = mygrid
!rab   nwat = Atm(n)%flagstruct%nwat
!rab
!rab!$OMP parallel do default (none) &
!rab!$OMP              shared (Atm_block, Atm, IPD_Data, npz, nq, ncnst, n, nwat) &
!rab!$OMP             private (nb, k, k1, ix, i, j, qwat)
!rab   do nb = 1,Atm_block%nblks
!rab     do k = 1, npz
!rab       k1 = npz+1-k !reverse the k direction
!rab       do ix = 1, Atm_block%blksz(nb)
!rab         i = Atm_block%index(nb)%ii(ix)
!rab         j = Atm_block%index(nb)%jj(ix)
!rab         qwat = sum(Atm(n)%q(i,j,k1,1:nwat))
!rab         Atm(n)%q(i,j,k1,1:nq) = Atm(n)%q(i,j,k1,1:nq) + IPD_Data(nb)%Stateout%gq0(ix,k,1:nq) * (1.0 - qwat)
!rab         if (nq .gt. ncnst) then
!rab           Atm(n)%qdiag(i,j,k1,nq+1:ncnst) = Atm(n)%qdiag(i,j,k1,nq+1:ncnst) + IPD_Data(nb)%Stateout%gq0(ix,k,nq+1:ncnst)
!rab         endif
!rab       enddo
!rab     enddo
!rab   enddo
!rab
!rab   call mpp_update_domains (Atm(n)%q, Atm(n)%domain, complete=.true.)
!rab
!rab   return
!rab end subroutine atmosphere_tracer_postinit


 subroutine get_bottom_mass ( t_bot, tr_bot, p_bot, z_bot, p_surf, slp )
!--------------------------------------------------------------
! returns temp, sphum, pres, height at the lowest model level
! and surface pressure
!--------------------------------------------------------------
   real, intent(out), dimension(isc:iec,jsc:jec):: t_bot, p_bot, z_bot, p_surf
   real, intent(out), optional, dimension(isc:iec,jsc:jec):: slp
   real, intent(out), dimension(isc:iec,jsc:jec,nq):: tr_bot
   integer :: i, j, m, k, kr
   real    :: rrg, sigtop, sigbot
   real, dimension(isc:iec,jsc:jec) :: tref
   real, parameter :: tlaps = 6.5e-3

   rrg  = rdgas / grav

   do j=jsc,jec
      do i=isc,iec
         p_surf(i,j) = Atm(mygrid)%ps(i,j)
         t_bot(i,j) = Atm(mygrid)%pt(i,j,npz)
         p_bot(i,j) = Atm(mygrid)%delp(i,j,npz)/(Atm(mygrid)%peln(i,npz+1,j)-Atm(mygrid)%peln(i,npz,j))
         z_bot(i,j) = rrg*t_bot(i,j)*(1.+zvir*Atm(mygrid)%q(i,j,npz,1)) *  &
                      (1. - Atm(mygrid)%pe(i,npz,j)/p_bot(i,j))
      enddo
   enddo

   if ( present(slp) ) then
     ! determine 0.8 sigma reference level
     sigtop = Atm(mygrid)%ak(1)/pstd_mks+Atm(mygrid)%bk(1)
     do k = 1, npz
        sigbot = Atm(mygrid)%ak(k+1)/pstd_mks+Atm(mygrid)%bk(k+1)
        if (sigbot+sigtop > 1.6) then
           kr = k
           exit
        endif
        sigtop = sigbot
     enddo
     do j=jsc,jec
        do i=isc,iec
           ! sea level pressure
           tref(i,j) = Atm(mygrid)%pt(i,j,kr) * (Atm(mygrid)%delp(i,j,kr)/ &
                            ((Atm(mygrid)%peln(i,kr+1,j)-Atm(mygrid)%peln(i,kr,j))*Atm(mygrid)%ps(i,j)))**(-rrg*tlaps)
           slp(i,j) = Atm(mygrid)%ps(i,j)*(1.+tlaps*Atm(mygrid)%phis(i,j)/(tref(i,j)*grav))**(1./(rrg*tlaps))
        enddo
     enddo
   endif

! Copy tracers
   do m=1,nq
      do j=jsc,jec
         do i=isc,iec
            tr_bot(i,j,m) = Atm(mygrid)%q(i,j,npz,m)
         enddo
      enddo
   enddo

 end subroutine get_bottom_mass


 subroutine get_bottom_wind ( u_bot, v_bot )
!-----------------------------------------------------------
! returns u and v on the mass grid at the lowest model level
!-----------------------------------------------------------
   real, intent(out), dimension(isc:iec,jsc:jec):: u_bot, v_bot
   integer i, j

   do j=jsc,jec
      do i=isc,iec
         u_bot(i,j) = Atm(mygrid)%u_srf(i,j)
         v_bot(i,j) = Atm(mygrid)%v_srf(i,j)
      enddo
   enddo

 end subroutine get_bottom_wind



 subroutine get_stock_pe(index, value)
   integer, intent(in) :: index
   real,   intent(out) :: value

#ifdef USE_STOCK
   include 'stock.inc'
#endif

   real wm(isc:iec,jsc:jec)
   integer i,j,k
   real, pointer :: area(:,:)

   area => Atm(mygrid)%gridstruct%area

   select case (index)

#ifdef USE_STOCK
   case (ISTOCK_WATER)
#else
   case (1)
#endif

!----------------------
! Perform vertical sum:
!----------------------
     wm = 0.
     do j=jsc,jec
        do k=1,npz
           do i=isc,iec
! Warning: the following works only with AM2 physics: water vapor; cloud water, cloud ice.
              wm(i,j) = wm(i,j) + Atm(mygrid)%delp(i,j,k) * ( Atm(mygrid)%q(i,j,k,1) +    &
                                                         Atm(mygrid)%q(i,j,k,2) +    &
                                                         Atm(mygrid)%q(i,j,k,3) )
           enddo
        enddo
     enddo

!----------------------
! Horizontal sum:
!----------------------
     value = 0.
     do j=jsc,jec
        do i=isc,iec
           value = value + wm(i,j)*area(i,j)
        enddo
     enddo
     value = value/grav

   case default
     value = 0.0
   end select

 end subroutine get_stock_pe


 subroutine atmosphere_state_update (Time, IPD_Data, IAU_Data, Atm_block)
   !--- interface variables ---
   type(time_type),           intent(in) :: Time
   type(IPD_data_type),       intent(in) :: IPD_Data(:)
   type(IAU_external_data_type), intent(in) :: IAU_Data
   type(block_control_type),  intent(in) :: Atm_block
   !--- local variables ---
   type(time_type) :: Time_prev, Time_next
   integer :: i, j, ix, k, k1, n, w_diff, nt_dyn, iq
   integer :: nb, blen, nwat, dnats, nq_adv
   real(kind=kind_phys):: rcp, q0, qwat(nq), qt, rdt
   real :: psum, qsum, psumb, qsumb, betad, psdt_mean
   real :: tracer_clock, lat_thresh, fhr
   character(len=32) :: tracer_name

   call timing_on('ATMOS_UPDATE')

   call mpp_clock_begin (id_update)

   Time_prev = Time
   Time_next = Time + Time_step_atmos
   rdt = 1.d0 / dt_atmos

   n = mygrid
   nwat = Atm(n)%flagstruct%nwat
   dnats = Atm(mygrid)%flagstruct%dnats
   nq_adv = nq - dnats
   nt_dyn = ncnst-pnats   !nothing more than nq

   if( nq<3 ) call mpp_error(FATAL, 'GFS phys must have 3 interactive tracers')

   if (IAU_Data%in_interval) then
      if (IAU_Data%drymassfixer) then
         ! global mean total pressure and water before IAU
         psumb = g_sum(Atm(n)%domain,sum(Atm(n)%delp(isc:iec,jsc:jec,1:npz),dim=3),&
                 isc,iec,jsc,jec,Atm(n)%ng,Atm(n)%gridstruct%area_64,1,reproduce=.true.)
         qsumb = g_sum(Atm(n)%domain,&
                 sum(Atm(n)%delp(isc:iec,jsc:jec,1:npz)*sum(Atm(n)%q(isc:iec,jsc:jec,1:npz,1:nwat),4),dim=3),&
                 isc,iec,jsc,jec,Atm(n)%ng,Atm(n)%gridstruct%area_64,1,reproduce=.true.)
         if (is_master()) then
           print *,'dry ps before IAU/physics',psumb+Atm(n)%ptop-qsumb
         endif
      endif

!     IAU increments are in units of 1/sec

!     add analysis increment to u,v,t tendencies
!     directly update delp with analysis increment
      do k = 1, npz
         do j = jsc,jec
            do i = isc,iec
               u_dt(i,j,k) = u_dt(i,j,k) + IAU_Data%ua_inc(i,j,k)
               v_dt(i,j,k) = v_dt(i,j,k) + IAU_Data%va_inc(i,j,k)
               t_dt(i,j,k) = t_dt(i,j,k) + IAU_Data%temp_inc(i,j,k)
               Atm(n)%delp(i,j,k) = Atm(n)%delp(i,j,k)  + IAU_Data%delp_inc(i,j,k)*dt_atmos
          enddo
        enddo
      enddo
      if (.not. Atm(n)%flagstruct%hydrostatic) then
         do k = 1, npz
            do j = jsc,jec
               do i = isc,iec
                  Atm(n)%delz(i,j,k) = Atm(n)%delz(i,j,k)  + IAU_Data%delz_inc(i,j,k)*dt_atmos
             enddo
           enddo
         enddo
      endif
!     add analysis increment to tracers to output from physics
      do nb = 1,Atm_block%nblks
         !if (nb.EQ.1) print*,'in block_update',IAU_Data%in_interval,IAU_Data%temp_inc(isc,jsc,30)
         blen = Atm_block%blksz(nb)
         do k = 1, npz
            k1 = npz+1-k !reverse the k direction
            do ix = 1, blen
               i = Atm_block%index(nb)%ii(ix)
               j = Atm_block%index(nb)%jj(ix)
               IPD_Data(nb)%Stateout%gq0(ix,k,:) = IPD_Data(nb)%Stateout%gq0(ix,k,:) + IAU_Data%tracer_inc(i,j,k1,:)*dt_atmos
            enddo
         enddo
      enddo
   endif

   call atmos_phys_qdt_diag(Atm(n)%q, Atm(n)%phys_diag, nt_dyn, dt_atmos, .true.)

!--- put u/v tendencies into haloed arrays u_dt and v_dt
!$OMP parallel do default (none) &
!$OMP              shared (rdt, n, nq, dnats, npz, ncnst, nwat, mygrid, u_dt, v_dt, t_dt,&
!$OMP                      Atm, IPD_Data, Atm_block, sphum, liq_wat, rainwat, ice_wat,   &
!$OMP                      snowwat, graupel, nq_adv)   &
!$OMP             private (nb, blen, i, j, k, k1, ix, q0, qwat, qt,tracer_name)
   do nb = 1,Atm_block%nblks

!SJL: perform vertical filling to fix the negative humidity if the SAS convection scheme is used
!     This call may be commented out if RAS or other positivity-preserving CPS is used.
     blen = Atm_block%blksz(nb)
     if (Atm(n)%flagstruct%fill_gfs) call fill_gfs(blen, npz, IPD_Data(nb)%Statein%prsi, IPD_Data(nb)%Stateout%gq0, 1.e-9_kind_phys)

!LMH 28sep18: If the name of a tracer ends in 'nopbl' then do NOT update it;
     !override this by setting Stateout%gq0(:,:,iq) to the input value
     do iq = 1, nq
        call get_tracer_names (MODEL_ATMOS, iq, tracer_name)
        if (index(tracer_name, 'nopbl') > 0) then
           IPD_Data(nb)%Stateout%gq0(:,:,iq) = IPD_Data(nb)%Statein%qgrs(:,:,iq)
        endif
     enddo

     do k = 1, npz
       k1 = npz+1-k !reverse the k direction
       do ix = 1, blen
         i = Atm_block%index(nb)%ii(ix)
         j = Atm_block%index(nb)%jj(ix)
         u_dt(i,j,k1) = u_dt(i,j,k1) + (IPD_Data(nb)%Stateout%gu0(ix,k) - IPD_Data(nb)%Statein%ugrs(ix,k)) * rdt
         v_dt(i,j,k1) = v_dt(i,j,k1) + (IPD_Data(nb)%Stateout%gv0(ix,k) - IPD_Data(nb)%Statein%vgrs(ix,k)) * rdt
         t_dt(i,j,k1) = t_dt(i,j,k1) + (IPD_Data(nb)%Stateout%gt0(ix,k) - IPD_Data(nb)%Statein%tgrs(ix,k)) * rdt
         !t_dt(i,j,k1) = (IPD_Data(nb)%Stateout%gt0(ix,k) - IPD_Data(nb)%Statein%tgrs(ix,k)) * rdt
! SJL notes:
! ---- DO not touch the code below; dry mass conservation may change due to 64bit <-> 32bit conversion
! GFS total air mass = dry_mass + water_vapor (condensate excluded)
! GFS mixing ratios  = tracer_mass / (dry_mass + vapor_mass)
! FV3 total air mass = dry_mass + [water_vapor + condensate ]
! FV3 mixing ratios  = tracer_mass / (dry_mass+vapor_mass+cond_mass)
         q0 = IPD_Data(nb)%Statein%prsi(ix,k) - IPD_Data(nb)%Statein%prsi(ix,k+1)
         qwat(1:nq_adv) = q0*IPD_Data(nb)%Stateout%gq0(ix,k,1:nq_adv)
! **********************************************************************************************************
! Dry mass: the following way of updating delp is key to mass conservation with hybrid 32-64 bit computation
! **********************************************************************************************************
! The following example is for 2 water species.
!        q0 = Atm(n)%delp(i,j,k1)*(1.-(Atm(n)%q(i,j,k1,1)+Atm(n)%q(i,j,k1,2))) + q1 + q2
         qt = sum(qwat(1:nwat))
         q0 = Atm(n)%delp(i,j,k1)*(1.-sum(Atm(n)%q(i,j,k1,1:nwat))) + qt
         Atm(n)%delp(i,j,k1) = q0
         Atm(n)%q(i,j,k1,1:nq_adv) = qwat(1:nq_adv) / q0
         if (dnats .gt. 0) Atm(n)%q(i,j,k1,nq_adv+1:nq) = IPD_Data(nb)%Stateout%gq0(ix,k,nq_adv+1:nq)
       enddo
     enddo

!GFDL     if ( dnats > 0 ) then
!GFDL       do iq = nq-dnats+1, nq
!GFDL         do k = 1, npz
!GFDL           k1 = npz+1-k !reverse the k direction
!GFDL           do ix = 1,blen
!GFDL             i = Atm_block%index(nb)%ii(ix)
!GFDL             j = Atm_block%index(nb)%jj(ix)
!GFDL             Atm(n)%q(i,j,k1,iq) = Atm(n)%q(i,j,k1,iq) + (IPD_Data(nb)%Stateout%gq0(ix,k,iq)-IPD_Data(nb)%Statein%qgrs(ix,k,iq))*rdt
!GFDL!                                * (IPD_Data%Statein%prsi(ix,k)-IPD_Data(nb)%Statein%prsi(ix,k+1))/Atm(n)%delp(i,j,k1)
!GFDL           enddo
!GFDL         enddo
!GFDL       enddo
!GFDL     endif

     !--- diagnostic tracers are assumed to be updated in-place
     !--- SHOULD THESE DIAGNOSTIC TRACERS BE MASS ADJUSTED???
     !--- See Note in statein...
     do iq = nq+1, ncnst
       do k = 1, npz
         k1 = npz+1-k !reverse the k direction
         do ix = 1, blen
           i = Atm_block%index(nb)%ii(ix)
           j = Atm_block%index(nb)%jj(ix)
           Atm(mygrid)%qdiag(i,j,k1,iq) = IPD_Data(nb)%Stateout%gq0(ix,k,iq)
         enddo
       enddo
     enddo

   enddo  ! nb-loop

! dry mass fixer in IAU interval following
! https://onlinelibrary.wiley.com/doi/full/10.1111/j.1600-0870.2007.00299.x
   if (IAU_Data%in_interval .and. IAU_data%drymassfixer) then
      ! global mean total pressure
      psum = g_sum(Atm(n)%domain,sum(Atm(n)%delp(isc:iec,jsc:jec,1:npz),dim=3),&
             isc,iec,jsc,jec,Atm(n)%ng,Atm(n)%gridstruct%area_64,1,reproduce=.true.)
      ! global mean total water (before adjustment)
      qsum = g_sum(Atm(n)%domain,&
             sum(Atm(n)%delp(isc:iec,jsc:jec,1:npz)*sum(Atm(n)%q(isc:iec,jsc:jec,1:npz,1:nwat),4),dim=3),&
             isc,iec,jsc,jec,Atm(n)%ng,Atm(n)%gridstruct%area_64,1,reproduce=.true.)
      betad = (psum - (psumb - qsumb))/qsum
      !if (is_master()) print *,'dry ps after IAU/physics',psum+Atm(n)%ptop-qsum, betad
      Atm(n)%q(:,:,:,1:nwat) = betad*Atm(n)%q(:,:,:,1:nwat)
      qsum = g_sum(Atm(n)%domain,&
             sum(Atm(n)%delp(isc:iec,jsc:jec,1:npz)*sum(Atm(n)%q(isc:iec,jsc:jec,1:npz,1:nwat),4),dim=3),&
             isc,iec,jsc,jec,Atm(n)%ng,Atm(n)%gridstruct%area_64,1)
      !if (is_master()) print *,'dry ps after iau_drymassfixer',psum+Atm(n)%ptop-qsum
   endif

   call atmos_phys_qdt_diag(Atm(n)%q, Atm(n)%phys_diag, nt_dyn, dt_atmos, .false.)

   w_diff = get_tracer_index (MODEL_ATMOS, 'w_diff' )
   if ( w_diff /= NO_TRACER ) then
      nt_dyn = nt_dyn - 1
   endif

!--- adjust w and heat tendency for non-hydrostatic case
#ifdef USE_Q_DT
    if ( .not.Atm(n)%flagstruct%hydrostatic .and. w_diff /= NO_TRACER ) then
      rcp = 1. / cp_air
!$OMP parallel do default (none) &
!$OMP              shared (jsc, jec, isc, iec, n, w_diff, Atm, q_dt, t_dt, rcp, dt_atmos) &
!$OMP             private (i, j, k)
       do k=1, Atm(n)%npz
         do j=jsc, jec
           do i=isc, iec
             Atm(n)%q(i,j,k,w_diff) = q_dt(i,j,k,w_diff) ! w tendency due to phys
! Heating due to loss of KE (vertical diffusion of w)
             t_dt(i,j,k) = t_dt(i,j,k) - q_dt(i,j,k,w_diff)*rcp*&
                                     (Atm(n)%w(i,j,k)+0.5*dt_atmos*q_dt(i,j,k,w_diff))
             Atm(n)%w(i,j,k) = Atm(n)%w(i,j,k) + dt_atmos*Atm(n)%q(i,j,k,w_diff)
           enddo
         enddo
       enddo
    endif
#endif

    call timing_on('FV_UPDATE_PHYS')
    call fv_update_phys( dt_atmos, isc, iec, jsc, jec, isd, ied, jsd, jed, Atm(n)%ng, nt_dyn, &
                         Atm(n)%u,  Atm(n)%v,   Atm(n)%w,  Atm(n)%delp, Atm(n)%pt,         &
                         Atm(n)%q,  Atm(n)%qdiag,                                          &
                         Atm(n)%ua, Atm(n)%va,  Atm(n)%ps, Atm(n)%pe,   Atm(n)%peln,       &
                         Atm(n)%pk, Atm(n)%pkz, Atm(n)%ak, Atm(n)%bk,   Atm(n)%phis,       &
                         Atm(n)%u_srf, Atm(n)%v_srf, Atm(n)%ts, Atm(n)%delz,               &
                         Atm(n)%flagstruct%hydrostatic, u_dt, v_dt, t_dt,                  &
                         .true., Time_next, Atm(n)%flagstruct%nudge, Atm(n)%gridstruct,    &
                         Atm(n)%gridstruct%agrid(:,:,1), Atm(n)%gridstruct%agrid(:,:,2),   &
                         Atm(n)%npx, Atm(n)%npy, Atm(n)%npz, Atm(n)%flagstruct,            &
                         Atm(n)%neststruct, Atm(n)%bd, Atm(n)%domain, &
                         Atm(n)%ptop, Atm(n)%phys_diag, Atm(n)%nudge_diag)
    call timing_off('FV_UPDATE_PHYS')

!MT surface pressure tendency (hPa/3hr)
   ps_dt(:,:)=(Atm(n)%ps(:,:)-ps_dt(:,:))*rdt*108.
   psdt_mean = g_sum(Atm(n)%domain,ABS(ps_dt(isc:iec,jsc:jec)),isc,iec,jsc,jec, &
                     Atm(n)%ng,Atm(n)%gridstruct%area_64,1,reproduce=.true.)
   do nb = 1,Atm_block%nblks
      blen = Atm_block%blksz(nb)
      do ix = 1, blen
         i = Atm_block%index(nb)%ii(ix)
         j = Atm_block%index(nb)%jj(ix)
         IPD_Data(nb)%intdiag%ps_dt(ix)=ps_dt(i,j)
      enddo
   enddo

!LMH 7jan2020: Update PBL and other clock tracers, if present
   tracer_clock = time_type_to_real(Time_next - Atm(n)%Time_init)*1.e-6
   lat_thresh = 15.*pi/180.
   do iq = 1, nq
      call get_tracer_names (MODEL_ATMOS, iq, tracer_name)
      if (trim(tracer_name) == 'pbl_clock' .or. trim(tracer_name) == 'tro_pbl_clock') then
         do nb = 1,Atm_block%nblks
            blen = Atm_block%blksz(nb)
            do ix = 1, blen
               i = Atm_block%index(nb)%ii(ix)
               j = Atm_block%index(nb)%jj(ix)
               if (trim(tracer_name) == 'tro_pbl_clock' .and. abs(Atm(n)%gridstruct%agrid(i,j,2)) > lat_thresh) cycle
               do k=1,npz
                  k1 = npz+1-k !reverse the k direction
                  Atm(n)%q(i,j,k1,iq) = tracer_clock
                  if (IPD_Data(nb)%Statein%phii(ix,k) > IPD_Data(nb)%intdiag%hpbl(ix)*grav) exit
               enddo
            enddo
         enddo
      else if (trim(tracer_name) == 'sfc_clock') then
         do j=jsc,jec
         do i=isc,iec
            Atm(n)%q(i,j,npz,iq) = tracer_clock
         enddo
         enddo
      else if (trim(tracer_name) == 'itcz_clock' ) then
         do k=1,npz
         do j=jsc,jec
         do i=isc,iec
            if (abs(Atm(n)%gridstruct%agrid(i,j,2)) < lat_thresh .and. Atm(n)%w(i,j,k) > 1.5) then
               Atm(n)%q(i,j,npz,iq) = tracer_clock
            endif
         enddo
         enddo
         enddo
      endif
  enddo

!--- nesting update after updating atmospheric variables with
!--- physics tendencies
    if (ngrids > 1 .and. p_split > 0) then
       call mpp_sync()
       call timing_on('TWOWAY_UPDATE')
       call twoway_nesting(Atm, ngrids, grids_on_this_pe, zvir, fv_time, mygrid)
       call timing_off('TWOWAY_UPDATE')
    endif

  call mpp_clock_end (id_update)

  call mpp_clock_begin(id_fv_diag)

  !---- diagnostics for FV dynamics -----
   if (Atm(mygrid)%flagstruct%print_freq /= -99) then

     fv_time = Time_next
     call get_time (fv_time, seconds,  days)

     call timing_on('FV_DIAG')
     call fv_diag(Atm(mygrid:mygrid), zvir, fv_time, Atm(mygrid)%flagstruct%print_freq)
     if (Atm(mygrid)%coarse_graining%write_coarse_diagnostics) then
        call fv_coarse_diag(Atm(mygrid:mygrid), fv_time, zvir)
     endif
     first_diag = .false.
     call timing_off('FV_DIAG')

   endif

   call mpp_clock_end(id_fv_diag)

   call timing_off('ATMOS_UPDATE')

 end subroutine atmosphere_state_update


 subroutine adiabatic_init(zvir,nudge_dz)
   real, allocatable, dimension(:,:,:):: u0, v0, t0, dz0, dp0
   real, intent(in):: zvir
   logical, intent(inout):: nudge_dz
!  real, parameter:: wt = 1.  ! was 2.
   real, parameter:: wt = 2.
!***********
! Haloe Data
!***********
   real, parameter::    q1_h2o = 2.2E-6
   real, parameter::    q7_h2o = 3.8E-6
   real, parameter::  q100_h2o = 3.8E-6
   real, parameter:: q1000_h2o = 3.1E-6
   real, parameter:: q2000_h2o = 2.8E-6
   real, parameter:: q3000_h2o = 3.0E-6
   real:: xt, p00, q00
   integer:: isc, iec, jsc, jec, npz
   integer:: m, n, i,j,k, ngc

   character(len=80) :: errstr

   xt = 1./(1.+wt)

   write(errstr,'(A, I4, A)') 'Performing adiabatic init',  Atm(mygrid)%flagstruct%na_init, ' times'
   call mpp_error(NOTE, errstr)
   sphum = get_tracer_index (MODEL_ATMOS, 'sphum' )

    npz = Atm(mygrid)%npz

    isc = Atm(mygrid)%bd%isc
    iec = Atm(mygrid)%bd%iec
    jsc = Atm(mygrid)%bd%jsc
    jec = Atm(mygrid)%bd%jec

    ngc = Atm(mygrid)%ng
    isd = isc - ngc
    ied = iec + ngc
    jsd = jsc - ngc
    jed = jec + ngc

     do_adiabatic_init = .true.

     allocate ( u0(isc:iec,  jsc:jec+1, npz) )
     allocate ( v0(isc:iec+1,jsc:jec,   npz) )
     allocate (dp0(isc:iec,jsc:jec, npz) )

     if ( Atm(mygrid)%flagstruct%hydrostatic ) nudge_dz = .false.

     if ( nudge_dz ) then
          allocate (dz0(isc:iec,jsc:jec, npz) )
     else
          allocate ( t0(isc:iec,jsc:jec, npz) )
     endif

!$omp parallel do default (none) &
!$omp              shared (nudge_dz, npz, jsc, jec, isc, iec, n, sphum, u0, v0, t0, dz0, dp0, Atm, zvir, mygrid) &
!$omp             private (k, j, i)
       do k=1,npz
          do j=jsc,jec+1
             do i=isc,iec
                u0(i,j,k) = Atm(mygrid)%u(i,j,k)
             enddo
          enddo
          do j=jsc,jec
             do i=isc,iec+1
                v0(i,j,k) = Atm(mygrid)%v(i,j,k)
             enddo
          enddo
          if ( nudge_dz ) then
             do j=jsc,jec
                do i=isc,iec
                   dp0(i,j,k) = Atm(mygrid)%delp(i,j,k)
                   dz0(i,j,k) = Atm(mygrid)%delz(i,j,k)
                enddo
             enddo
          else
             do j=jsc,jec
                do i=isc,iec
                   t0(i,j,k) = Atm(mygrid)%pt(i,j,k)*(1.+zvir*Atm(mygrid)%q(i,j,k,sphum))  ! virt T
                   dp0(i,j,k) = Atm(mygrid)%delp(i,j,k)
                enddo
             enddo
          endif
       enddo

     do m=1,Atm(mygrid)%flagstruct%na_init
! Forward call
    call fv_dynamics(Atm(mygrid)%npx, Atm(mygrid)%npy, npz,  nq, Atm(mygrid)%ng, dt_atmos, 0.,      &
                     Atm(mygrid)%flagstruct%fill, Atm(mygrid)%flagstruct%reproduce_sum, kappa, cp_air, zvir,  &
                     Atm(mygrid)%ptop, Atm(mygrid)%ks, nq, Atm(mygrid)%flagstruct%n_split,        &
                     Atm(mygrid)%flagstruct%q_split, Atm(mygrid)%u0, Atm(mygrid)%v0, Atm(mygrid)%u, &
                     Atm(mygrid)%v, Atm(mygrid)%w, Atm(mygrid)%delz, Atm(mygrid)%flagstruct%hydrostatic, &
                     Atm(mygrid)%pt, Atm(mygrid)%delp, Atm(mygrid)%q, Atm(mygrid)%ps,                     &
                     Atm(mygrid)%pe, Atm(mygrid)%pk, Atm(mygrid)%peln, Atm(mygrid)%pkz, Atm(mygrid)%phis,      &
                     Atm(mygrid)%q_con, Atm(mygrid)%omga, Atm(mygrid)%ua, Atm(mygrid)%va, Atm(mygrid)%uc, Atm(mygrid)%vc, &
                     Atm(mygrid)%ak, Atm(mygrid)%bk, Atm(mygrid)%mfx, Atm(mygrid)%mfy,                    &
                     Atm(mygrid)%cx, Atm(mygrid)%cy, Atm(mygrid)%ze0, Atm(mygrid)%flagstruct%hybrid_z,    &
                     Atm(mygrid)%gridstruct, Atm(mygrid)%flagstruct,                            &
                     Atm(mygrid)%neststruct, Atm(mygrid)%idiag, Atm(mygrid)%bd, Atm(mygrid)%parent_grid,  &
                     Atm(mygrid)%domain, Atm(mygrid)%inline_mp, Atm(mygrid)%diss_est)
! Backward
    call fv_dynamics(Atm(mygrid)%npx, Atm(mygrid)%npy, npz,  nq, Atm(mygrid)%ng, -dt_atmos, 0.,      &
                     Atm(mygrid)%flagstruct%fill, Atm(mygrid)%flagstruct%reproduce_sum, kappa, cp_air, zvir,  &
                     Atm(mygrid)%ptop, Atm(mygrid)%ks, nq, Atm(mygrid)%flagstruct%n_split,        &
                     Atm(mygrid)%flagstruct%q_split, Atm(mygrid)%u0, Atm(mygrid)%v0, Atm(mygrid)%u, &
                     Atm(mygrid)%v, Atm(mygrid)%w, Atm(mygrid)%delz, Atm(mygrid)%flagstruct%hydrostatic, &
                     Atm(mygrid)%pt, Atm(mygrid)%delp, Atm(mygrid)%q, Atm(mygrid)%ps,                     &
                     Atm(mygrid)%pe, Atm(mygrid)%pk, Atm(mygrid)%peln, Atm(mygrid)%pkz, Atm(mygrid)%phis,      &
                     Atm(mygrid)%q_con, Atm(mygrid)%omga, Atm(mygrid)%ua, Atm(mygrid)%va, Atm(mygrid)%uc, Atm(mygrid)%vc, &
                     Atm(mygrid)%ak, Atm(mygrid)%bk, Atm(mygrid)%mfx, Atm(mygrid)%mfy,                    &
                     Atm(mygrid)%cx, Atm(mygrid)%cy, Atm(mygrid)%ze0, Atm(mygrid)%flagstruct%hybrid_z,    &
                     Atm(mygrid)%gridstruct, Atm(mygrid)%flagstruct,                            &
                     Atm(mygrid)%neststruct, Atm(mygrid)%idiag, Atm(mygrid)%bd, Atm(mygrid)%parent_grid,  &
                     Atm(mygrid)%domain, Atm(mygrid)%inline_mp, Atm(mygrid)%diss_est)
! Nudging back to IC
!$omp parallel do default (none) &
!$omp              shared (pref, npz, jsc, jec, isc, iec, n, sphum, Atm, u0, v0, t0, dp0, xt, zvir, mygrid, nudge_dz, dz0) &
!$omp             private (i, j, k, p00, q00)
       do k=1,npz
          do j=jsc,jec+1
             do i=isc,iec
                Atm(mygrid)%u(i,j,k) = xt*(Atm(mygrid)%u(i,j,k) + wt*u0(i,j,k))
             enddo
          enddo
          do j=jsc,jec
             do i=isc,iec+1
                Atm(mygrid)%v(i,j,k) = xt*(Atm(mygrid)%v(i,j,k) + wt*v0(i,j,k))
             enddo
          enddo
          if( Atm(mygrid)%flagstruct%nudge_qv ) then
! SJL note: Nudging water vaport towards HALOE climatology:
! In case of better IC (IFS) this step may not be necessary
             p00 = Atm(mygrid)%pe(isc,k,jsc)
             if ( p00 < 30.E2 ) then
                if ( p00 < 1. ) then
                     q00 = q1_h2o
                elseif ( p00 <= 7. .and. p00 >= 1. ) then
                     q00 = q1_h2o + (q7_h2o-q1_h2o)*log(pref(k,1)/1.)/log(7.)
                elseif ( p00 < 100. .and. p00 >= 7. ) then
                     q00 = q7_h2o + (q100_h2o-q7_h2o)*log(pref(k,1)/7.)/log(100./7.)
                elseif ( p00 < 1000. .and. p00 >= 100. ) then
                     q00 = q100_h2o + (q1000_h2o-q100_h2o)*log(pref(k,1)/1.E2)/log(10.)
                elseif ( p00 < 2000. .and. p00 >= 1000. ) then
                     q00 = q1000_h2o + (q2000_h2o-q1000_h2o)*log(pref(k,1)/1.E3)/log(2.)
                else
                     q00 = q2000_h2o + (q3000_h2o-q2000_h2o)*log(pref(k,1)/2.E3)/log(1.5)
                endif
                do j=jsc,jec
                   do i=isc,iec
                      Atm(mygrid)%q(i,j,k,sphum) = xt*(Atm(mygrid)%q(i,j,k,sphum) + wt*q00)
                   enddo
                enddo
             endif
          endif
          if ( nudge_dz ) then
             do j=jsc,jec
                do i=isc,iec
                   Atm(mygrid)%delp(i,j,k) = xt*(Atm(mygrid)%delp(i,j,k) + wt*dp0(i,j,k))
                   Atm(mygrid)%delz(i,j,k) = xt*(Atm(mygrid)%delz(i,j,k) + wt*dz0(i,j,k))
                enddo
             enddo
          else
             do j=jsc,jec
                do i=isc,iec
                   Atm(mygrid)%pt(i,j,k) = xt*(Atm(mygrid)%pt(i,j,k) + wt*t0(i,j,k)/(1.+zvir*Atm(mygrid)%q(i,j,k,sphum)))
                   Atm(mygrid)%delp(i,j,k) = xt*(Atm(mygrid)%delp(i,j,k) + wt*dp0(i,j,k))
                enddo
             enddo
          endif

       enddo

! Backward
    call fv_dynamics(Atm(mygrid)%npx, Atm(mygrid)%npy, npz,  nq, Atm(mygrid)%ng, -dt_atmos, 0.,      &
                     Atm(mygrid)%flagstruct%fill, Atm(mygrid)%flagstruct%reproduce_sum, kappa, cp_air, zvir,  &
                     Atm(mygrid)%ptop, Atm(mygrid)%ks, nq, Atm(mygrid)%flagstruct%n_split,        &
                     Atm(mygrid)%flagstruct%q_split, Atm(mygrid)%u0, Atm(mygrid)%v0, Atm(mygrid)%u, &
                     Atm(mygrid)%v, Atm(mygrid)%w, Atm(mygrid)%delz, Atm(mygrid)%flagstruct%hydrostatic, &
                     Atm(mygrid)%pt, Atm(mygrid)%delp, Atm(mygrid)%q, Atm(mygrid)%ps,                     &
                     Atm(mygrid)%pe, Atm(mygrid)%pk, Atm(mygrid)%peln, Atm(mygrid)%pkz, Atm(mygrid)%phis,      &
                     Atm(mygrid)%q_con, Atm(mygrid)%omga, Atm(mygrid)%ua, Atm(mygrid)%va, Atm(mygrid)%uc, Atm(mygrid)%vc, &
                     Atm(mygrid)%ak, Atm(mygrid)%bk, Atm(mygrid)%mfx, Atm(mygrid)%mfy,                    &
                     Atm(mygrid)%cx, Atm(mygrid)%cy, Atm(mygrid)%ze0, Atm(mygrid)%flagstruct%hybrid_z,    &
                     Atm(mygrid)%gridstruct, Atm(mygrid)%flagstruct,                            &
                     Atm(mygrid)%neststruct, Atm(mygrid)%idiag, Atm(mygrid)%bd, Atm(mygrid)%parent_grid,  &
                     Atm(mygrid)%domain, Atm(mygrid)%inline_mp, Atm(mygrid)%diss_est)
! Forward call
    call fv_dynamics(Atm(mygrid)%npx, Atm(mygrid)%npy, npz,  nq, Atm(mygrid)%ng, dt_atmos, 0.,      &
                     Atm(mygrid)%flagstruct%fill, Atm(mygrid)%flagstruct%reproduce_sum, kappa, cp_air, zvir,  &
                     Atm(mygrid)%ptop, Atm(mygrid)%ks, nq, Atm(mygrid)%flagstruct%n_split,        &
                     Atm(mygrid)%flagstruct%q_split, Atm(mygrid)%u0, Atm(mygrid)%v0, Atm(mygrid)%u, &
                     Atm(mygrid)%v, Atm(mygrid)%w, Atm(mygrid)%delz, Atm(mygrid)%flagstruct%hydrostatic, &
                     Atm(mygrid)%pt, Atm(mygrid)%delp, Atm(mygrid)%q, Atm(mygrid)%ps,                     &
                     Atm(mygrid)%pe, Atm(mygrid)%pk, Atm(mygrid)%peln, Atm(mygrid)%pkz, Atm(mygrid)%phis,      &
                     Atm(mygrid)%q_con, Atm(mygrid)%omga, Atm(mygrid)%ua, Atm(mygrid)%va, Atm(mygrid)%uc, Atm(mygrid)%vc, &
                     Atm(mygrid)%ak, Atm(mygrid)%bk, Atm(mygrid)%mfx, Atm(mygrid)%mfy,                    &
                     Atm(mygrid)%cx, Atm(mygrid)%cy, Atm(mygrid)%ze0, Atm(mygrid)%flagstruct%hybrid_z,    &
                     Atm(mygrid)%gridstruct, Atm(mygrid)%flagstruct,                            &
                     Atm(mygrid)%neststruct, Atm(mygrid)%idiag, Atm(mygrid)%bd, Atm(mygrid)%parent_grid,  &
                     Atm(mygrid)%domain, Atm(mygrid)%inline_mp, Atm(mygrid)%diss_est)
! Nudging back to IC
!$omp parallel do default (none) &
!$omp              shared (nudge_dz,npz, jsc, jec, isc, iec, n, sphum, Atm, u0, v0, t0, dz0, dp0, xt, zvir, mygrid) &
!$omp             private (i, j, k)
       do k=1,npz
          do j=jsc,jec+1
             do i=isc,iec
                Atm(mygrid)%u(i,j,k) = xt*(Atm(mygrid)%u(i,j,k) + wt*u0(i,j,k))
             enddo
          enddo
          do j=jsc,jec
             do i=isc,iec+1
                Atm(mygrid)%v(i,j,k) = xt*(Atm(mygrid)%v(i,j,k) + wt*v0(i,j,k))
             enddo
          enddo
          if ( nudge_dz ) then
             do j=jsc,jec
             do i=isc,iec
                Atm(mygrid)%delp(i,j,k) = xt*(Atm(mygrid)%delp(i,j,k) + wt*dp0(i,j,k))
                Atm(mygrid)%delz(i,j,k) = xt*(Atm(mygrid)%delz(i,j,k) + wt*dz0(i,j,k))
             enddo
             enddo
          else
             do j=jsc,jec
             do i=isc,iec
                Atm(mygrid)%pt(i,j,k) = xt*(Atm(mygrid)%pt(i,j,k) + wt*t0(i,j,k)/(1.+zvir*Atm(mygrid)%q(i,j,k,sphum)))
                Atm(mygrid)%delp(i,j,k) = xt*(Atm(mygrid)%delp(i,j,k) + wt*dp0(i,j,k))
             enddo
             enddo
          endif
       enddo

     enddo

     deallocate ( u0 )
     deallocate ( v0 )
     deallocate (dp0 )
     if ( allocated(t0) )  deallocate ( t0 )
     if ( allocated(dz0) ) deallocate ( dz0 )

     do_adiabatic_init = .false.

 end subroutine adiabatic_init



#if defined(OVERLOAD_R4)
#define _DBL_(X) DBLE(X)
#define _RL_(X) REAL(X,KIND=4)
#else
#define _DBL_(X) X
#define _RL_(X) X
#endif
 subroutine atmos_phys_driver_statein (IPD_Data, Atm_block)
   type (IPD_data_type),      intent(inout) :: IPD_Data(:)
   type (block_control_type), intent(in)    :: Atm_block
!--------------------------------------
! Local GFS-phys consistent parameters:
!--------------------------------------
   real(kind=kind_phys), parameter:: p00 = 1.e5
   real(kind=kind_phys), parameter:: qmin = 1.0e-10
   real(kind=kind_phys):: pk0inv, ptop, pktop
   real(kind=kind_phys) :: rTv, dm, qgrs_rad
   integer :: nb, blen, npz, i, j, k, ix, k1, dnats, nq_adv, isd, ied, jsd, jed

   real, pointer :: omega_for_physics(:,:,:)
!!! NOTES: lmh 6nov15
!!! - "Layer" means "layer mean", ie. the average value in a layer
!!! - "Level" means "level interface", ie the point values at the top or bottom of a layer

   ptop =  _DBL_(_RL_(Atm(mygrid)%ak(1)))
   pktop  = (ptop/p00)**kappa
   pk0inv = (1.0_kind_phys/p00)**kappa

   npz = Atm_block%npz
   dnats = Atm(mygrid)%flagstruct%dnats
   nq_adv = nq - dnats

   if (.not. Atm(mygrid)%flagstruct%hydrostatic .and. .not. Atm(mygrid)%flagstruct%pass_full_omega_to_physics_in_non_hydrostatic_mode) then
      omega_for_physics => Atm(mygrid)%local_omga
   else
      omega_for_physics => Atm(mygrid)%omga
   endif

!---------------------------------------------------------------------
! use most up to date atmospheric properties when running serially
!---------------------------------------------------------------------
!$OMP parallel do default (none) &
!$OMP             shared  (Atm_block, Atm, IPD_Data, npz, nq, ncnst, sphum, liq_wat, &
!$OMP                      ice_wat, rainwat, snowwat, graupel, pk0inv, ptop,   &
!$OMP                      pktop, zvir, mygrid, dnats, nq_adv, omega_for_physics) &
!$OMP             private (dm, nb, blen, i, j, ix, k1, rTv, qgrs_rad)

   do nb = 1,Atm_block%nblks
! gas_phase_mass <-- prsl
! log(pe) <-- prsik

     blen = Atm_block%blksz(nb)

     !-- level interface geopotential height (relative to the surface)
     IPD_Data(nb)%Statein%phii(:,1) = 0.0_kind_phys
     IPD_Data(nb)%Statein%prsik(:,:) = 1.e25_kind_phys

     if (Atm(mygrid)%flagstruct%do_inline_mp) then
         do ix = 1, blen
           i = Atm_block%index(nb)%ii(ix)
           j = Atm_block%index(nb)%jj(ix)
           IPD_Data(nb)%Statein%prew(ix) = _DBL_(_RL_(Atm(mygrid)%inline_mp%prew(i,j)))
           IPD_Data(nb)%Statein%prer(ix) = _DBL_(_RL_(Atm(mygrid)%inline_mp%prer(i,j)))
           IPD_Data(nb)%Statein%prei(ix) = _DBL_(_RL_(Atm(mygrid)%inline_mp%prei(i,j)))
           IPD_Data(nb)%Statein%pres(ix) = _DBL_(_RL_(Atm(mygrid)%inline_mp%pres(i,j)))
           IPD_Data(nb)%Statein%preg(ix) = _DBL_(_RL_(Atm(mygrid)%inline_mp%preg(i,j)))
           if (Atm(mygrid)%flagstruct%do_cosp) then
             do k = 1, npz
               k1 = npz+1-k ! flipping the index
               IPD_Data(nb)%Statein%prefluxw(ix,k) = _DBL_(_RL_(Atm(mygrid)%inline_mp%prefluxw(i,j,k1)))
               IPD_Data(nb)%Statein%prefluxr(ix,k) = _DBL_(_RL_(Atm(mygrid)%inline_mp%prefluxr(i,j,k1)))
               IPD_Data(nb)%Statein%prefluxi(ix,k) = _DBL_(_RL_(Atm(mygrid)%inline_mp%prefluxi(i,j,k1)))
               IPD_Data(nb)%Statein%prefluxs(ix,k) = _DBL_(_RL_(Atm(mygrid)%inline_mp%prefluxs(i,j,k1)))
               IPD_Data(nb)%Statein%prefluxg(ix,k) = _DBL_(_RL_(Atm(mygrid)%inline_mp%prefluxg(i,j,k1)))
             enddo
           endif
         enddo
     endif

     do k = 1, npz
       do ix = 1, blen
         i = Atm_block%index(nb)%ii(ix)
         j = Atm_block%index(nb)%jj(ix)

            !Indices for FV's vertical coordinate, for which 1 = top
            !here, k is the index for GFS's vertical coordinate, for which 1 = bottom
         k1 = npz+1-k ! flipping the index
         IPD_Data(nb)%Statein%tgrs(ix,k) = _DBL_(_RL_(Atm(mygrid)%pt(i,j,k1)))
         IPD_Data(nb)%Statein%ugrs(ix,k) = _DBL_(_RL_(Atm(mygrid)%ua(i,j,k1)))
         IPD_Data(nb)%Statein%vgrs(ix,k) = _DBL_(_RL_(Atm(mygrid)%va(i,j,k1)))
          IPD_Data(nb)%Statein%vvl(ix,k) = _DBL_(_RL_(omega_for_physics(i,j,k1)))
         IPD_Data(nb)%Statein%prsl(ix,k) = _DBL_(_RL_(Atm(mygrid)%delp(i,j,k1)))   ! Total mass
         if (Atm(mygrid)%flagstruct%do_diss_est)IPD_Data(nb)%Statein%diss_est(ix,k) = _DBL_(_RL_(Atm(mygrid)%diss_est(i,j,k1)))

         if (.not.Atm(mygrid)%flagstruct%hydrostatic .and. (.not.Atm(mygrid)%flagstruct%use_hydro_pressure))  &
           IPD_Data(nb)%Statein%phii(ix,k+1) = IPD_Data(nb)%Statein%phii(ix,k) - _DBL_(_RL_(Atm(mygrid)%delz(i,j,k1)*grav))

! Convert to tracer mass:
         IPD_Data(nb)%Statein%qgrs(ix,k,1:nq_adv) =  _DBL_(_RL_(Atm(mygrid)%q(i,j,k1,1:nq_adv))) &
                                                          * IPD_Data(nb)%Statein%prsl(ix,k)

         if (dnats .gt. 0) &
             IPD_Data(nb)%Statein%qgrs(ix,k,nq_adv+1:nq) =  _DBL_(_RL_(Atm(mygrid)%q(i,j,k1,nq_adv+1:nq)))
         !--- SHOULD THESE BE CONVERTED TO MASS SINCE THE DYCORE DOES NOT TOUCH THEM IN ANY WAY???
         !--- See Note in state update...
         if ( ncnst > nq) &
             IPD_Data(nb)%Statein%qgrs(ix,k,nq+1:ncnst) = _DBL_(_RL_(Atm(mygrid)%qdiag(i,j,k1,nq+1:ncnst)))
! Remove the contribution of condensates to delp (mass):
         if ( Atm(mygrid)%flagstruct%nwat .eq. 6 ) then
            IPD_Data(nb)%Statein%prsl(ix,k) = IPD_Data(nb)%Statein%prsl(ix,k) &
                                            - IPD_Data(nb)%Statein%qgrs(ix,k,liq_wat)   &
                                            - IPD_Data(nb)%Statein%qgrs(ix,k,ice_wat)   &
                                            - IPD_Data(nb)%Statein%qgrs(ix,k,rainwat)   &
                                            - IPD_Data(nb)%Statein%qgrs(ix,k,snowwat)   &
                                            - IPD_Data(nb)%Statein%qgrs(ix,k,graupel)
         else !variable condensate numbers
            IPD_Data(nb)%Statein%prsl(ix,k) = IPD_Data(nb)%Statein%prsl(ix,k) &
                                            - sum(IPD_Data(nb)%Statein%qgrs(ix,k,2:Atm(mygrid)%flagstruct%nwat))
         endif
       enddo
     enddo

! Re-compute pressure (dry_mass + water_vapor) derived fields:
     do i=1,blen
        IPD_Data(nb)%Statein%prsi(i,npz+1) = ptop
     enddo
     do k=npz,1,-1
        do i=1,blen
           IPD_Data(nb)%Statein%prsi(i,k)  = IPD_Data(nb)%Statein%prsi(i,k+1)  &
                                           + IPD_Data(nb)%Statein%prsl(i,k)
           IPD_Data(nb)%Statein%prsik(i,k) = log( IPD_Data(nb)%Statein%prsi(i,k) )
! Redefine mixing ratios for GFS == tracer_mass / (dry_air_mass + water_vapor_mass)
           IPD_Data(nb)%Statein%qgrs(i,k,1:nq_adv) = IPD_Data(nb)%Statein%qgrs(i,k,1:nq_adv) &
                                                   / IPD_Data(nb)%Statein%prsl(i,k)
        enddo
     enddo
     do i=1,blen
        IPD_Data(nb)%Statein%pgr(i)         = IPD_Data(nb)%Statein%prsi(i,1)    ! surface pressure for GFS
        IPD_Data(nb)%Statein%prsik(i,npz+1) = log(ptop)
     enddo

     do k=1,npz
        do i=1,blen
! Geo-potential at interfaces:
           qgrs_rad = max(qmin,IPD_Data(nb)%Statein%qgrs(i,k,sphum))
           rTv = rdgas*IPD_Data(nb)%Statein%tgrs(i,k)*(1.+zvir*qgrs_rad)
           if ( Atm(mygrid)%flagstruct%hydrostatic .or. Atm(mygrid)%flagstruct%use_hydro_pressure )   &
                IPD_Data(nb)%Statein%phii(i,k+1) = IPD_Data(nb)%Statein%phii(i,k) &
                                                     + rTv*(IPD_Data(nb)%Statein%prsik(i,k) &
                                                          - IPD_Data(nb)%Statein%prsik(i,k+1))
! Layer mean pressure by perfect gas law:
           dm = IPD_Data(nb)%Statein%prsl(i,k)
           IPD_Data(nb)%Statein%prsl(i,k) = dm*rTv/(IPD_Data(nb)%Statein%phii(i,k+1) &
                                                  - IPD_Data(nb)%Statein%phii(i,k))

!!! Ensure subgrid MONOTONICITY of Pressure: SJL 09/11/2016
           if ( .not.Atm(mygrid)%flagstruct%hydrostatic ) then
! If violated, replaces it with hydrostatic pressure
              IPD_Data(nb)%Statein%prsl(i,k) = min(IPD_Data(nb)%Statein%prsl(i,k), &
                                                   IPD_Data(nb)%Statein%prsi(i,k)   - 0.01*dm)
              IPD_Data(nb)%Statein%prsl(i,k) = max(IPD_Data(nb)%Statein%prsl(i,k), &
                                                   IPD_Data(nb)%Statein%prsi(i,k+1) + 0.01*dm)
           endif
        enddo
     enddo

     do k = 1,npz
        do i=1,blen
! Exner function layer center: large sensitivity to non-hydro runs with moist kappa
           IPD_Data(nb)%Statein%prslk(i,k) = exp( kappa*log(IPD_Data(nb)%Statein%prsl(i,k)/p00) )
!--  layer center geopotential; geometric midpoint
           IPD_Data(nb)%Statein%phil(i,k) = 0.5_kind_phys*(IPD_Data(nb)%Statein%phii(i,k) &
                                                         + IPD_Data(nb)%Statein%phii(i,k+1))
        enddo
     enddo

! Compute Exner function at layer "interfaces"
    do i=1,blen
!      Top & Bottom edges computed hydrostaticlly
       IPD_Data(nb)%Statein%prsik(i,    1) = exp( kappa*IPD_Data(nb)%Statein%prsik(i,1) )*pk0inv  ! bottom
       IPD_Data(nb)%Statein%prsik(i,npz+1) = pktop                           ! TOA
    enddo

    if ( Atm(mygrid)%flagstruct%hydrostatic .or. Atm(mygrid)%flagstruct%use_hydro_pressure ) then
        do k=2,npz
           do i=1,blen
              IPD_Data(nb)%Statein%prsik(i,k) = exp( kappa*IPD_Data(nb)%Statein%prsik(i,k) )*pk0inv
           enddo
        enddo
    endif
    IPD_Data(nb)%Statein%nwat = Atm(mygrid)%flagstruct%nwat
  enddo

 end subroutine atmos_phys_driver_statein

 subroutine atmos_phys_qdt_diag(q, phys_diag, nq, dt, begin)

   integer, intent(IN) :: nq
   real, intent(IN) :: dt
   logical, intent(IN) :: begin
   real, intent(IN) :: q(isd:ied,jsd:jed,npz,nq)
   type(phys_diag_type), intent(INOUT) :: phys_diag

   integer  sphum, liq_wat, ice_wat        ! GFDL AM physics
   integer  rainwat, snowwat, graupel      ! GFDL Cloud Microphysics

   sphum   = get_tracer_index (MODEL_ATMOS, 'sphum')
   liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
   ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
   rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat')
   snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat')
   graupel = get_tracer_index (MODEL_ATMOS, 'graupel')

   if (begin) then
      if (allocated(phys_diag%phys_qv_dt)) phys_diag%phys_qv_dt = q(isc:iec,jsc:jec,:,sphum)
      if (allocated(phys_diag%phys_ql_dt)) then
         if (liq_wat < 0) call mpp_error(FATAL, " phys_ql_dt needs at least one liquid water tracer defined")
         phys_diag%phys_ql_dt = q(isc:iec,jsc:jec,:,liq_wat)
      endif
      if (allocated(phys_diag%phys_qi_dt)) then
         if (ice_wat < 0) then
            call mpp_error(WARNING, " phys_qi_dt needs at least one ice water tracer defined")
            phys_diag%phys_qi_dt = 0.
         endif
         phys_diag%phys_qi_dt = q(isc:iec,jsc:jec,:,ice_wat)
      endif

      if (liq_wat > 0) then
         if (allocated(phys_diag%phys_liq_wat_dt)) phys_diag%phys_liq_wat_dt = q(isc:iec,jsc:jec,:,liq_wat)
      endif

      if (rainwat > 0) then
         if (allocated(phys_diag%phys_qr_dt)) phys_diag%phys_qr_dt = q(isc:iec,jsc:jec,:,rainwat)
      endif

      if (ice_wat > 0) then
         if (allocated(phys_diag%phys_ice_wat_dt)) phys_diag%phys_ice_wat_dt = q(isc:iec,jsc:jec,:,ice_wat)
      endif

      if (graupel > 0) then
         if (allocated(phys_diag%phys_qg_dt)) phys_diag%phys_qg_dt = q(isc:iec,jsc:jec,:,graupel)
      endif

      if (snowwat > 0) then
         if (allocated(phys_diag%phys_qs_dt)) phys_diag%phys_qs_dt = q(isc:iec,jsc:jec,:,snowwat)
      endif
   else
      if (allocated(phys_diag%phys_qv_dt)) phys_diag%phys_qv_dt = q(isc:iec,jsc:jec,:,sphum) - phys_diag%phys_qv_dt
      if (allocated(phys_diag%phys_ql_dt)) then
         phys_diag%phys_ql_dt = q(isc:iec,jsc:jec,:,liq_wat) - phys_diag%phys_ql_dt
      endif
      if (allocated(phys_diag%phys_qi_dt)) then
         phys_diag%phys_qi_dt = q(isc:iec,jsc:jec,:,ice_wat) - phys_diag%phys_qi_dt
      endif
   endif

   if (allocated(phys_diag%phys_ql_dt)) then
      if (rainwat > 0) phys_diag%phys_ql_dt = q(isc:iec,jsc:jec,:,rainwat) + phys_diag%phys_ql_dt
   endif
   if (allocated(phys_diag%phys_qi_dt)) then
      if (snowwat > 0) phys_diag%phys_qi_dt = q(isc:iec,jsc:jec,:,snowwat) + phys_diag%phys_qi_dt
      if (graupel > 0) phys_diag%phys_qi_dt = q(isc:iec,jsc:jec,:,graupel) + phys_diag%phys_qi_dt
   endif

   if (.not. begin) then
      if (allocated(phys_diag%phys_qv_dt)) phys_diag%phys_qv_dt = phys_diag%phys_qv_dt / dt
      if (allocated(phys_diag%phys_ql_dt)) phys_diag%phys_ql_dt = phys_diag%phys_ql_dt / dt
      if (allocated(phys_diag%phys_qi_dt)) phys_diag%phys_qi_dt = phys_diag%phys_qi_dt / dt

      if (liq_wat > 0) then
         if (allocated(phys_diag%phys_liq_wat_dt)) phys_diag%phys_liq_wat_dt = (q(isc:iec,jsc:jec,:,liq_wat) - phys_diag%phys_liq_wat_dt) / dt
      endif

      if (rainwat > 0) then
         if (allocated(phys_diag%phys_qr_dt)) phys_diag%phys_qr_dt = (q(isc:iec,jsc:jec,:,rainwat) - phys_diag%phys_qr_dt) / dt
      endif

      if (ice_wat > 0) then
         if (allocated(phys_diag%phys_ice_wat_dt)) phys_diag%phys_ice_wat_dt = (q(isc:iec,jsc:jec,:,ice_wat) - phys_diag%phys_ice_wat_dt) / dt
      endif

      if (graupel > 0) then
         if (allocated(phys_diag%phys_qg_dt)) phys_diag%phys_qg_dt = (q(isc:iec,jsc:jec,:,graupel) - phys_diag%phys_qg_dt) / dt
      endif

      if (snowwat > 0) then
         if (allocated(phys_diag%phys_qs_dt)) phys_diag%phys_qs_dt = (q(isc:iec,jsc:jec,:,snowwat) - phys_diag%phys_qs_dt) / dt
      endif
   endif

 end subroutine atmos_phys_qdt_diag

 subroutine atmosphere_coarse_graining_parameters(coarse_domain, write_coarse_restart_files, write_only_coarse_intermediate_restarts)
   type(domain2d), intent(out) :: coarse_domain
   logical, intent(out) :: write_coarse_restart_files, write_only_coarse_intermediate_restarts

   coarse_domain = Atm(mygrid)%coarse_graining%domain
   write_coarse_restart_files = Atm(mygrid)%coarse_graining%write_coarse_restart_files
   write_only_coarse_intermediate_restarts = Atm(mygrid)%coarse_graining%write_only_coarse_intermediate_restarts
 end subroutine atmosphere_coarse_graining_parameters

 subroutine atmosphere_coarsening_strategy(coarsening_strategy)
   character(len=64), intent(out) :: coarsening_strategy

   coarsening_strategy = Atm(mygrid)%coarse_graining%strategy
 end subroutine atmosphere_coarsening_strategy

end module atmosphere_mod
