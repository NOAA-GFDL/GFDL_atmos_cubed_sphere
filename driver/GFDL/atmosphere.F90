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
use atmos_co2_mod,         only: atmos_co2_rad, co2_radiation_override
use block_control_mod,     only: block_control_type
use constants_mod,         only: cp_air, rdgas, grav, rvgas, kappa, pstd_mks
use time_manager_mod,      only: time_type, get_time, set_time, operator(+), &
                                 operator(-), operator(/), time_type_to_real
use fms_mod,               only: error_mesg, FATAL,                 &
                                 check_nml_error, stdlog,           &
                                 write_version_number,              &
                                 mpp_clock_id, mpp_clock_begin,     &
                                 mpp_clock_end, CLOCK_SUBCOMPONENT, &
                                 clock_flag_default
use fms2_io_mod,           only: file_exists
use mpp_mod,               only: mpp_error, FATAL, NOTE, input_nml_file, &
                                 mpp_npes, mpp_get_current_pelist, &
                                 mpp_set_current_pelist, stdout, &
                                 mpp_pe, mpp_root_pe, mpp_chksum, mpp_sync
use mpp_domains_mod,       only: domain2d
use xgrid_mod,             only: grid_box_type
!miz
use diag_manager_mod,      only: register_diag_field, send_data
use field_manager_mod,     only: MODEL_ATMOS
use tracer_manager_mod,    only: get_tracer_index, get_number_tracers, &
                                 NO_TRACER, get_tracer_names
use physics_driver_mod,    only: surf_diff_type
use physics_types_mod,     only: physics_type, &
                                 physics_tendency_type
use radiation_types_mod,   only: radiation_type, compute_g_avg
use atmos_cmip_diag_mod,   only: atmos_cmip_diag_init, &
                                 register_cmip_diag_field_3d, &
                                 send_cmip_data_3d, cmip_diag_id_type, &
                                 query_cmip_diag_id
use atmos_global_diag_mod, only: atmos_global_diag_init, &
                                 atmos_global_diag_end

!-----------------
! FV core modules:
!-----------------
use fv_arrays_mod,      only: fv_atmos_type
use fv_control_mod,     only: fv_control_init, fv_end, ngrids
use fv_eta_mod,         only: get_eta_level
use fv_dynamics_mod,    only: fv_dynamics
use fv_nesting_mod,     only: twoway_nesting
use fv_diagnostics_mod, only: fv_diag_init, fv_diag, fv_time, prt_maxmin, prt_height, Mw_air!_3d
use fv_cmip_diag_mod,   only: fv_cmip_diag_init, fv_cmip_diag, fv_cmip_diag_end
use fv_restart_mod,     only: fv_restart, fv_write_restart
use fv_timing_mod,      only: timing_on, timing_off
use fv_mp_mod,          only: is_master
use fv_sg_mod,          only: fv_subgrid_z
use fv_update_phys_mod, only: fv_update_phys
use fv_io_mod,          only: fv_io_register_nudge_restart
use fv_regional_mod,    only: start_regional_restart, read_new_bc_data
use fv_regional_mod,    only: a_step, p_step
use fv_regional_mod,    only: current_time_in_seconds
#if defined (ATMOS_NUDGE)
use atmos_nudge_mod,      only: atmos_nudge_init, atmos_nudge_end
#elif defined (CLIMATE_NUDGE)
use fv_climate_nudge_mod, only: fv_climate_nudge_init,fv_climate_nudge_end
#elif defined (ADA_NUDGE)
use fv_ada_nudge_mod,     only: fv_ada_nudge_init, fv_ada_nudge_end
#else
use fv_nwp_nudge_mod,     only: fv_nwp_nudge_init, fv_nwp_nudge_end, do_adiabatic_init
use amip_interp_mod,      only: forecast_mode
#endif

use mpp_domains_mod, only:  mpp_get_data_domain, mpp_get_compute_domain
use gfdl_mp_mod,        only: gfdl_mp_init, gfdl_mp_end
use coarse_graining_mod, only: coarse_graining_init
use coarse_grained_diagnostics_mod, only: fv_coarse_diag_init, fv_coarse_diag
use coarse_grained_restart_files_mod, only: fv_coarse_restart_init

implicit none
private

!--- driver routines
public :: atmosphere_init, atmosphere_end, atmosphere_restart, &
          atmosphere_dynamics, atmosphere_state_update

!--- utility routines
public :: atmosphere_resolution, atmosphere_boundary, &
          atmosphere_grid_center, atmosphere_domain, &
          atmosphere_cell_area, atmosphere_control_data, &
          atmosphere_pref, &
          get_atmosphere_axes, get_bottom_mass, &
          get_bottom_wind, get_stock_pe, &
          set_atmosphere_pelist, reset_atmos_tracers

!--- physics/radiation data exchange routines
public :: atmos_radiation_driver_inputs, atmos_physics_driver_inputs

!-----------------------------------------------------------------------
! version number of this module
! Include variable "version" to be written to log file.
#include<file_version.h>
character(len=20)   :: mod_name = 'GFDL/atmosphere_mod'

!---- private data ----
  type (time_type) :: Time_step_atmos
  public Atm

  !These are convenience variables for local use only, and are set to values in Atm%
  real    :: dt_atmos
  real    :: zvir
  integer :: npx, npy, npz, ncnst, pnats
  integer :: isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed
  integer :: nq                       ! transported tracers
  integer :: sec, seconds, days
  integer :: id_dynam, id_fv_diag, id_subgridz
  logical :: cold_start = .false.       ! read in initial condition

  integer, dimension(:), allocatable :: id_tracerdt_dyn, id_tracer, id_tracerdt_dyn_dp
  integer :: num_tracers = 0

!miz
  !Diagnostics
  type(cmip_diag_id_type) :: ID_tnta, ID_tnhusa, ID_tnt, ID_tnhus
  integer :: id_udt_dyn, id_vdt_dyn, id_tdt_dyn, id_qdt_dyn
  integer :: id_qldt_dyn, id_qidt_dyn, id_qadt_dyn
  logical :: used
  character(len=64) :: field
  real, allocatable :: ttend(:,:,:)
  real, allocatable :: qtendyyf(:,:,:,:), qtendyyf_dp(:,:,:,:), mw_air_store(:,:,:)
  real, allocatable :: qtend(:,:,:,:)
  real              :: mv = -1.e10   ! missing value for diagnostics
  integer :: sphum, liq_wat, rainwat, ice_wat, snowwat, graupel  !condensate species
  integer :: cld_amt
!miz

  integer :: mygrid = 1
  integer :: p_split = 1
  integer, allocatable :: pelist(:)
  logical, allocatable :: grids_on_this_pe(:)
  type(fv_atmos_type), allocatable, target :: Atm(:)

  real, parameter:: w0_big = 60.  ! to prevent negative w-tracer diffusion

!---dynamics tendencies for use in fv_subgrid_z and during fv_update_phys
  real, allocatable, dimension(:,:,:)   :: u_dt, v_dt, t_dt, qv_dt
  real, allocatable, dimension(:,:,:,:) :: q_dt
  real, allocatable :: pref(:,:), dum1d(:)

!---need to define do_adiabatic_init to satisfy a reference when nwp_nudge is not the default
#if defined(ATMOS_NUDGE) || defined(CLIMATE_NUDGE) || defined(ADA_NUDGE)
   logical :: do_adiabatic_init
#endif

   integer, parameter :: kind_phys=8

   logical, allocatable :: is_vmr(:)

contains



 subroutine atmosphere_init (Time_init, Time, Time_step, Surf_diff, Grid_box)
   type (time_type),      intent(in)    :: Time_init, Time, Time_step
   type(surf_diff_type),  intent(inout) :: Surf_diff
   type(grid_box_type),   intent(inout) :: Grid_box

!--- local variables ---
   integer :: i, n
   integer :: itrac
   logical :: do_atmos_nudge
   character(len=32) :: tracer_name, tracer_units
   real :: ps1, ps2

   integer :: nlunit = 9999
   character (len = 64) :: fn_nml = 'input.nml'

   !For regional
   a_step = 0
   current_time_in_seconds = time_type_to_real( Time - Time_init )
   if (mpp_pe() == 0) write(0,"('atmosphere_init: current_time_seconds = ',f9.1)")current_time_in_seconds

                    call timing_on('ATMOS_INIT')
   allocate(pelist(mpp_npes()))
   call mpp_get_current_pelist(pelist)

   call get_number_tracers(MODEL_ATMOS, num_prog= num_tracers)

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
   cld_amt = get_tracer_index (MODEL_ATMOS, 'cld_amt' )

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

!----- allocate and zero out the dynamics (and accumulated) tendencies
   allocate( u_dt(isd:ied,jsd:jed,npz), &
             v_dt(isd:ied,jsd:jed,npz), &
             t_dt(isc:iec,jsc:jec,npz), &
             qv_dt(isc:iec,jsc:jec,npz), &
             q_dt(isc:iec,jsc:jec,npz,nq) )
!--- allocate pref
   allocate(pref(npz+1,2), dum1d(npz+1))

   if (Atm(mygrid)%flagstruct%do_inline_mp) then
     call gfdl_mp_init(input_nml_file, stdlog(), Atm(mygrid)%flagstruct%hydrostatic)
   endif

   call fv_restart(Atm(mygrid)%domain, Atm, dt_atmos, seconds, days, cold_start, &
                   Atm(mygrid)%gridstruct%grid_type, mygrid)

   fv_time = Time

!----- initialize atmos_axes and fv_dynamics diagnostics
       !I've had trouble getting this to work with multiple grids at a time; worth revisiting?
   call fv_diag_init(Atm(mygrid:mygrid), Atm(mygrid)%atmos_axes, Time, npx, npy, npz, Atm(mygrid)%flagstruct%p_ref)

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
   id_dynam     = mpp_clock_id ('FV dy-core',  flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
   id_subgridz  = mpp_clock_id ('FV subgrid_z',flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
   id_fv_diag   = mpp_clock_id ('FV Diag',     flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
!---- initialize cmip diagnostic output ----
   call atmos_cmip_diag_init   ( Atm(mygrid)%ak, Atm(mygrid)%bk, pref(1,1), Atm(mygrid)%atmos_axes, Time )
   call atmos_global_diag_init ( Atm(mygrid)%atmos_axes, Atm(mygrid)%gridstruct%area(isc:iec,jsc:jec) )
   call fv_cmip_diag_init      ( Atm(mygrid:mygrid), Atm(mygrid)%atmos_axes, Time )

!--- initialize nudging module ---
#if defined (ATMOS_NUDGE)
    call atmos_nudge_init ( Time, Atm(mygrid)%atmos_axes(1:3), flag=do_atmos_nudge )
    if ( do_atmos_nudge .and. Atm(mygrid)%flagstruct%nudge ) then
         call mpp_error(NOTE, 'Code compiled with atmospheric nudging, but fv_core_nml nudge is also set to .true.')
    elseif ( do_atmos_nudge) then
         call mpp_error(NOTE, 'Code compiled with and using atmospheric nudging')
    endif
    Atm(mygrid)%flagstruct%nudge = do_atmos_nudge
#elif defined (CLIMATE_NUDGE)
    call fv_climate_nudge_init ( Time, Atm(mygrid)%atmos_axes(1:3), flag=do_atmos_nudge )
    if ( do_atmos_nudge .and. Atm(mygrid)%flagstruct%nudge ) then
         call mpp_error(NOTE, 'Code compiled with climate nudging, but fv_core_nml nudge is also set to .true.')
    elseif ( do_atmos_nudge ) then
         call mpp_error(NOTE, 'Code compiled with and using climate nudging')
    endif
    Atm(mygrid)%flagstruct%nudge = do_atmos_nudge
#elif defined (ADA_NUDGE)
    if ( Atm(mygrid)%flagstruct%nudge ) then
        call fv_ada_nudge_init( Time, Atm(mygrid)%atmos_axes, npz, zvir, Atm(mygrid)%ak, Atm(mygrid)%bk, &
           Atm(mygrid)%ts, Atm(mygrid)%phis, Atm(mygrid)%gridstruct, Atm(mygrid)%ks, Atm(mygrid)%npx,    &
           Atm(mygrid)%neststruct, Atm(mygrid)%bd, Atm(mygrid)%domain_for_read)
        call mpp_error(NOTE, 'ADA nudging is active')
     endif
#else
   !Only do nudging on coarse grid for now
   if ( Atm(mygrid)%flagstruct%nudge ) then
      call fv_nwp_nudge_init( Time, Atm(mygrid)%atmos_axes, npz, zvir, Atm(mygrid)%ak, Atm(mygrid)%bk, &
           Atm(mygrid)%ts, Atm(mygrid)%phis, Atm(mygrid)%gridstruct, Atm(mygrid)%ks, Atm(mygrid)%npx,  &
           Atm(mygrid)%neststruct, Atm(mygrid)%bd)
        call mpp_error(NOTE, 'NWP nudging is active')
   endif
#endif

!  --- initiate the start for a restarted regional forecast
   if ( Atm(mygrid)%gridstruct%regional .and. Atm(mygrid)%flagstruct%warm_start ) then
     call start_regional_restart(Atm(1),       &
                                 isc, iec, jsc, jec, &
                                 isd, ied, jsd, jed )
   endif

! This call needs to be separate from the register nudging restarts after initialization
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

   !This appears to all be diagnostics through the end of this routine,
   !and so for now we will only define for the coarsest grid

!miz
!---allocate id_tracer_*
   allocate (id_tracerdt_dyn    (num_tracers))
   allocate (id_tracerdt_dyn_dp (num_tracers))
   allocate (is_vmr             (num_tracers))
   is_vmr(:) = .false.
   if ( Atm(mygrid)%flagstruct%write_3d_diags) then
      id_udt_dyn    =register_diag_field(mod_name,'udt_dyn', Atm(mygrid)%atmos_axes(1:3),  &
           Time,'udt_dyn',    'm/s/s', missing_value=mv)
      id_vdt_dyn    =register_diag_field(mod_name,'vdt_dyn', Atm(mygrid)%atmos_axes(1:3),  &
           Time,'vdt_dyn',    'm/s/s', missing_value=mv)
      id_tdt_dyn    =register_diag_field(mod_name,'tdt_dyn', Atm(mygrid)%atmos_axes(1:3),  &
           Time,'tdt_dyn',    'K/s', missing_value=mv)
      id_qdt_dyn    =register_diag_field(mod_name,'qdt_dyn', Atm(mygrid)%atmos_axes(1:3),  &
           Time,'qdt_dyn',    'kg/kg/s', missing_value=mv)
      id_qldt_dyn   =register_diag_field(mod_name,'qldt_dyn', Atm(mygrid)%atmos_axes(1:3),  &
           Time,'qldt_dyn',   'kg/kg/s', missing_value=mv)
      id_qidt_dyn   =register_diag_field(mod_name,'qidt_dyn', Atm(mygrid)%atmos_axes(1:3),  &
           Time,'qidt_dyn',   'kg/kg/s', missing_value=mv)
      id_qadt_dyn   =register_diag_field(mod_name,'qadt_dyn', Atm(mygrid)%atmos_axes(1:3),  &
           Time,'qadt_dyn',   '1/s', missing_value=mv)
      !--- register cmip tendency fields ---
      ID_tnta = register_cmip_diag_field_3d (mod_name, 'tnta', Time, &
           'Tendency of Air Temperature due to Advection', 'K s-1', &
           standard_name='tendency_of_air_temperature_due_to_advection')
      ID_tnhusa = register_cmip_diag_field_3d (mod_name, 'tnhusa', Time, &
           'Tendency of Specific Humidity due to Advection', 's-1', &
           standard_name='tendency_of_specific_humidity_due_to_advection')
      ID_tnt = register_cmip_diag_field_3d (mod_name, 'tnt', Time, &
           'Tendency of Air Temperature', 'K s-1', &
           standard_name='tendency_of_air_temperature')
      ID_tnhus = register_cmip_diag_field_3d (mod_name, 'tnhus', Time, &
           'Tendency of Specific Humidity', 's-1', &
           standard_name='tendency_of_specific_humidity')

      !---loop for tracers
      do itrac = 1, num_tracers
         call get_tracer_names (MODEL_ATMOS, itrac, name = tracer_name, units = tracer_units)
         if (get_tracer_index(MODEL_ATMOS,tracer_name)>0) then
            id_tracerdt_dyn(itrac) = register_diag_field(mod_name, TRIM(tracer_name)//'dt_dyn',  &
                 Atm(mygrid)%atmos_axes(1:3),Time,                       &
                 TRIM(tracer_name)//' total tendency from advection',    &
                 TRIM(tracer_units)//'/s', missing_value = mv)

            if ( trim(tracer_units) .eq. 'vmr' ) then
               is_vmr(itrac) = .true.
               id_tracerdt_dyn_dp(itrac) = register_diag_field(mod_name, TRIM(tracer_name)//'dt_dyn_dp',  &
                    Atm(mygrid)%atmos_axes(1:3),Time,                       &
                    TRIM(tracer_name)//' total tendency from advection',    &
                    'mol/m2/s', missing_value = mv)
            else
               id_tracerdt_dyn_dp(itrac) = register_diag_field(mod_name, TRIM(tracer_name)//'dt_dyn_dp',  &
                    Atm(mygrid)%atmos_axes(1:3),Time,                       &
                    TRIM(tracer_name)//' total tendency from advection',    &
                    'kg/m2/s', missing_value = mv)
            end if
         endif
      enddo
   endif
   if (any(id_tracerdt_dyn(:)>0)) allocate(qtendyyf(isc:iec, jsc:jec,1:npz,num_tracers))
   if (any(id_tracerdt_dyn_dp(:)>0)) allocate(qtendyyf_dp(isc:iec, jsc:jec,1:npz,num_tracers))
   allocate(mw_air_store(isc:iec, jsc:jec,1:npz))
   if ( id_tdt_dyn>0 .or. query_cmip_diag_id(ID_tnta) .or. query_cmip_diag_id(ID_tnt) ) &
                                                      allocate(ttend(isc:iec, jsc:jec, 1:npz))
   if ( any((/ id_qdt_dyn, id_qldt_dyn, id_qidt_dyn, id_qadt_dyn /) > 0) .or. &
        query_cmip_diag_id(ID_tnhusa) .or. query_cmip_diag_id(ID_tnhus) )  allocate(qtend(isc:iec, jsc:jec, 1:npz, 4))

! could zero out diagnostics if tracer field not defined
   if (sphum > size(qtend,4)) id_qdt_dyn = 0
   if (liq_wat > size(qtend,4)) id_qldt_dyn = 0
   if (ice_wat > size(qtend,4)) id_qidt_dyn = 0
   if (cld_amt > size(qtend,4)) id_qadt_dyn = 0
!miz

!  --- initialize clocks for dynamics, physics_down and physics_up
   id_dynam     = mpp_clock_id ('FV dy-core',  flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
   id_subgridz  = mpp_clock_id ('FV subgrid_z',flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
   id_fv_diag   = mpp_clock_id ('FV Diag',     flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )

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


 subroutine atmosphere_dynamics ( Time, Surf_diff )
   type(time_type),intent(in) :: Time
   integer :: itrac, n, psc
   integer :: k, w_diff, nt_dyn
   type(surf_diff_type), intent(inout):: Surf_diff
   logical :: used
   type(time_type) :: atmos_time
   integer :: atmos_time_step
   real :: rdt
!---- Call FV dynamics -----

   call mpp_clock_begin (id_dynam)

   Surf_diff%ddp_dyn(:,:,:) = Atm(mygrid)%delp(isc:iec, jsc:jec, :)
   Surf_diff%tdt_dyn(:,:,:) = Atm(mygrid)%pt(isc:iec, jsc:jec, :)
   Surf_diff%qdt_dyn(:,:,:) = Atm(mygrid)%q (isc:iec, jsc:jec, :, sphum) + &
                              Atm(mygrid)%q (isc:iec, jsc:jec, :, liq_wat) + &
                              Atm(mygrid)%q (isc:iec, jsc:jec, :, ice_wat)

!miz
   if ( id_tdt_dyn>0 .or. query_cmip_diag_id(ID_tnta) ) ttend(:, :, :) = Atm(mygrid)%pt(isc:iec, jsc:jec, :)
   if ( any((/ id_qdt_dyn, id_qldt_dyn, id_qidt_dyn, id_qadt_dyn /) > 0) .or. &
        query_cmip_diag_id(ID_tnhusa) ) qtend(:, :, :, 1:4) = Atm(mygrid)%q (isc:iec, jsc:jec, :, 1:4)

   mw_air_store = Mw_air(Atm(mygrid)%q (isc:iec, jsc:jec, :, sphum)   + &
                         Atm(mygrid)%q (isc:iec, jsc:jec, :, liq_wat) + &
                         Atm(mygrid)%q (isc:iec, jsc:jec, :, ice_wat)) !g/mol

   do itrac = 1, num_tracers
     if (id_tracerdt_dyn (itrac) >0 ) &
          qtendyyf(:,:,:,itrac) = Atm(mygrid)%q(isc:iec,jsc:jec,:,itrac)
   end do

   !convert vmr tracers to mmr tracers before dynamics
   if (Atm(mygrid)%flagstruct%adj_mass_vmr .eq. 2) then
      do itrac = 1, num_tracers
         if (is_vmr(itrac)) then
            Atm(mygrid)%q(isc:iec,jsc:jec,:,itrac) = Atm(mygrid)%q(isc:iec,jsc:jec,:,itrac)/mw_air_store
         end if
      end do
   end if

   do itrac = 1, num_tracers
     if (id_tracerdt_dyn_dp (itrac) >0 ) then
        qtendyyf_dp(:,:,:,itrac) = Atm(mygrid)%q(isc:iec,jsc:jec,:,itrac)*Atm(mygrid)%delp(isc:iec,jsc:jec,:)/grav
     end if
  enddo

   n = mygrid
   a_step = a_step + 1
!
!*** If this is a regional run then read in the next boundary data when it is time.
!
   if(Atm(n)%flagstruct%regional)then

     call read_new_bc_data(Atm(n), Time, Time_step_atmos, p_split, &
                           isd, ied, jsd, jed )
   endif

   do psc=1,abs(p_split)
      p_step = psc
                    call timing_on('fv_dynamics')
!uc/vc only need be same on coarse grid? However BCs do need to be the same
     call fv_dynamics(npx, npy, npz, nq, Atm(n)%ng, dt_atmos/real(abs(p_split)),&
                      Atm(n)%flagstruct%consv_te, Atm(n)%flagstruct%fill,  &
                      Atm(n)%flagstruct%reproduce_sum, kappa, cp_air, zvir,&
                      Atm(n)%ptop, Atm(n)%ks, nq,                          &
                      Atm(n)%flagstruct%n_split, Atm(n)%flagstruct%q_split,&
                      Atm(n)%u, Atm(n)%v, Atm(n)%w, Atm(n)%delz,           &
                      Atm(n)%flagstruct%hydrostatic,                       &
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
                      Atm(n)%diss_est)
     call timing_off('fv_dynamics')

    if (ngrids > 1 .and. (psc < p_split .or. p_split < 0)) then
       call mpp_sync()
       call timing_on('TWOWAY_UPDATE')
       call twoway_nesting(Atm, ngrids, grids_on_this_pe, zvir, fv_time, mygrid)
       call timing_off('TWOWAY_UPDATE')
    endif

    end do !p_split
    call mpp_clock_end (id_dynam)

   Surf_diff%ddp_dyn(:,:,:) =(Atm(mygrid)%delp(isc:iec,jsc:jec,:)-Surf_diff%ddp_dyn(:,:,:))/dt_atmos
   Surf_diff%tdt_dyn(:,:,:) =(Atm(mygrid)%pt(isc:iec,jsc:jec,:)  -Surf_diff%tdt_dyn(:,:,:))/dt_atmos
   Surf_diff%qdt_dyn(:,:,:) =(Atm(mygrid)%q (isc:iec,jsc:jec,:,sphum) + &
                              Atm(mygrid)%q (isc:iec,jsc:jec,:,liq_wat) + &
                              Atm(mygrid)%q (isc:iec,jsc:jec,:,ice_wat) - Surf_diff%qdt_dyn(:,:,:))/dt_atmos

!miz
   if (id_udt_dyn>0)  used = send_data( id_udt_dyn, 2.0/dt_atmos*Atm(mygrid)%ua(isc:iec,jsc:jec,:), Time)
   if (id_vdt_dyn>0)  used = send_data( id_vdt_dyn, 2.0/dt_atmos*Atm(mygrid)%va(isc:iec,jsc:jec,:), Time)
   if (id_tdt_dyn > 0) used = send_data( id_tdt_dyn, (Atm(mygrid)%pt(isc:iec,jsc:jec,:)-ttend(:,:,:))/dt_atmos, Time)
   if (query_cmip_diag_id(ID_tnta)) &
                 used = send_cmip_data_3d ( ID_tnta, (Atm(mygrid)%pt(isc:iec,jsc:jec,:)-ttend(:,:,:))/dt_atmos, Time)

   if (id_qdt_dyn  > 0) used = send_data( id_qdt_dyn , (Atm(mygrid)%q(isc:iec,jsc:jec,:,sphum)-qtend(:,:,:,sphum))/dt_atmos, Time)
   if (id_qldt_dyn > 0) used = send_data( id_qldt_dyn, (Atm(mygrid)%q(isc:iec,jsc:jec,:,liq_wat)-qtend(:,:,:,liq_wat))/dt_atmos, Time)
   if (id_qidt_dyn > 0) used = send_data( id_qidt_dyn, (Atm(mygrid)%q(isc:iec,jsc:jec,:,ice_wat)-qtend(:,:,:,ice_wat))/dt_atmos, Time)
   if (id_qadt_dyn > 0) used = send_data( id_qadt_dyn, (Atm(mygrid)%q(isc:iec,jsc:jec,:,cld_amt)-qtend(:,:,:,cld_amt))/dt_atmos, Time)
   if (query_cmip_diag_id(ID_tnhusa)) &
                  used = send_cmip_data_3d (ID_tnhusa, (Atm(mygrid)%q(isc:iec,jsc:jec,:,sphum)-qtend(:,:,:,sphum))/dt_atmos, Time)
!miz

   mw_air_store = Mw_air(Atm(mygrid)%q (isc:iec, jsc:jec, :, sphum)   + &
                         Atm(mygrid)%q (isc:iec, jsc:jec, :, liq_wat) + &
                         Atm(mygrid)%q (isc:iec, jsc:jec, :, ice_wat)) !g/mol

   do itrac = 1, num_tracers
      if(id_tracerdt_dyn(itrac)>0) then
         if (is_vmr(itrac) .and. Atm(mygrid)%flagstruct%adj_mass_vmr.eq. 2 ) then         !if adj_mass_vmr == 1, no conversion to mmr is done, so do not convert
            qtendyyf(:,:,:,itrac) = (Atm(mygrid)%q (isc:iec, jsc:jec, :,itrac)*mw_air_store -  &
                 qtendyyf(:,:,:,itrac))/dt_atmos
         else
            qtendyyf(:,:,:,itrac) = (Atm(mygrid)%q (isc:iec, jsc:jec, :,itrac) -  &
                 qtendyyf(:,:,:,itrac))/dt_atmos
         end if
       used = send_data(id_tracerdt_dyn(itrac), qtendyyf(:,:,:,itrac), Time)
     endif

     if(id_tracerdt_dyn_dp(itrac)>0) then
        qtendyyf_dp(:,:,:,itrac) = (Atm(mygrid)%q (isc:iec, jsc:jec, :,itrac)*Atm(mygrid)%delp(isc:iec,jsc:jec,:)/grav -  qtendyyf_dp(:,:,:,itrac))/dt_atmos

        !this is in kg/m2/s. For vmr we need to convert to mol/m2/s [we used 1g/mol for conversion purposes]
        if (is_vmr(itrac)) then
           qtendyyf_dp(:,:,:,itrac) = qtendyyf_dp(:,:,:,itrac) * 1e3
        end if
       used = send_data(id_tracerdt_dyn_dp(itrac), qtendyyf_dp(:,:,:,itrac), Time)
     endif

  enddo

!-----------------------------------------------------
!--- COMPUTE SUBGRID Z
!-----------------------------------------------------
!--- zero out tendencies
    call mpp_clock_begin (id_subgridz)
    u_dt(:,:,:)   = 0. ! These are updated by fv_subgrid_z
    v_dt(:,:,:)   = 0.
! t_dt is used for two different purposes:
!    1 - to calculate the diagnostic temperature tendency from fv_subgrid_z
!    2 - as an accumulator for the IAU increment and physics tendency
! because of this, it will need to be zeroed out after the diagnostic is calculated
    t_dt(:,:,:)   = Atm(n)%pt(isc:iec,jsc:jec,:)
    qv_dt(:,:,:)  = Atm(n)%q (isc:iec,jsc:jec,:,sphum)
    q_dt(:,:,:,:) = 0.

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
                        Atm(n)%w, Atm(n)%delz, u_dt, v_dt, t_dt, q_dt,      &
                        Atm(n)%flagstruct%n_sponge)
   endif

   !convert back to vmr if needed
   if (Atm(mygrid)%flagstruct%adj_mass_vmr.eq. 2) then
      mw_air_store = Mw_air(Atm(mygrid)%q (isc:iec, jsc:jec, :, sphum)   + &
                            Atm(mygrid)%q (isc:iec, jsc:jec, :, liq_wat) + &
                            Atm(mygrid)%q (isc:iec, jsc:jec, :, ice_wat)) !g/mol

      do itrac = 1, num_tracers
         if (is_vmr(itrac)) then
            Atm(mygrid)%q(isc:iec,jsc:jec,:,itrac) = Atm(mygrid)%q(isc:iec,jsc:jec,:,itrac)*mw_air_store
         end if
      end do
   end if


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

    if (Atm(1)%idiag%id_u_dt_sg > 0) then
       used = send_data(Atm(1)%idiag%id_u_dt_sg, u_dt(isc:iec,jsc:jec,:), fv_time)
    end if
    if (Atm(1)%idiag%id_v_dt_sg > 0) then
       used = send_data(Atm(1)%idiag%id_v_dt_sg, v_dt(isc:iec,jsc:jec,:), fv_time)
    end if
    if (Atm(1)%idiag%id_t_dt_sg > 0) then
       t_dt(:,:,:) = rdt*(Atm(1)%pt(isc:iec,jsc:jec,:) - t_dt(:,:,:))
       used = send_data(Atm(1)%idiag%id_t_dt_sg, t_dt, fv_time)
    end if
    if (Atm(1)%idiag%id_qv_dt_sg > 0) then
       qv_dt(:,:,:) = rdt*(Atm(1)%q(isc:iec,jsc:jec,:,sphum) - qv_dt(:,:,:))
       used = send_data(Atm(1)%idiag%id_qv_dt_sg, qv_dt, fv_time)
    end if

! zero out t_dt for use as an accumulator
    t_dt = 0.

   call mpp_clock_end (id_subgridz)

 end subroutine atmosphere_dynamics


 subroutine atmosphere_end (Time, Grid_box )
   type (time_type),      intent(in)    :: Time
   type(grid_box_type),   intent(inout) :: Grid_box

!--- end nudging module ---
#if defined (ATMOS_NUDGE)
   if ( Atm(mygrid)%flagstruct%nudge ) call atmos_nudge_end
#elif defined (CLIMATE_NUDGE)
   if ( Atm(mygrid)%flagstruct%nudge ) call fv_climate_nudge_end
#elif defined (ADA_NUDGE)
   if ( Atm(mygrid)%flagstruct%nudge ) call fv_ada_nudge_end
#else
   if ( Atm(mygrid)%flagstruct%nudge ) call fv_nwp_nudge_end
#endif

   if (Atm(mygrid)%flagstruct%do_inline_mp) then
     call gfdl_mp_end ( )
   endif

      call timing_on('FV_DIAG')
   call atmos_global_diag_end
   call fv_cmip_diag_end
   call fv_end(Atm, mygrid)
      call timing_off('FV_DIAG')

   deallocate ( Atm )
   deallocate ( u_dt, v_dt, t_dt, qv_dt, q_dt, pref, dum1d )
   deallocate ( is_vmr )

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

 subroutine atmosphere_control_data (i1, i2, j1, j2, kt, p_hydro, hydro, do_uni_zfull) !miz
   integer, intent(out)           :: i1, i2, j1, j2, kt
   logical, intent(out), optional :: p_hydro, hydro, do_uni_zfull !miz
   i1 = Atm(mygrid)%bd%isc
   i2 = Atm(mygrid)%bd%iec
   j1 = Atm(mygrid)%bd%jsc
   j2 = Atm(mygrid)%bd%jec
   kt = Atm(mygrid)%npz

   if (present(p_hydro)) p_hydro = Atm(mygrid)%flagstruct%phys_hydrostatic
   if (present(  hydro))   hydro = Atm(mygrid)%flagstruct%hydrostatic
   if (present(do_uni_zfull)) do_uni_zfull = Atm(mygrid)%flagstruct%do_uni_zfull

 end subroutine atmosphere_control_data


 subroutine atmosphere_cell_area  (area_out)
   real, dimension(:,:),  intent(out)          :: area_out

   area_out(1:iec-isc+1, 1:jec-jsc+1) =  Atm(mygrid)%gridstruct%area (isc:iec,jsc:jec)

 end subroutine atmosphere_cell_area



 subroutine atmosphere_grid_center (lon, lat)
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

 end subroutine atmosphere_grid_center



 subroutine atmosphere_boundary (blon, blat, global)
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

 end subroutine atmosphere_boundary


 subroutine set_atmosphere_pelist ()
   call mpp_set_current_pelist(Atm(mygrid)%pelist, no_sync=.TRUE.)
 end subroutine set_atmosphere_pelist


 subroutine atmosphere_domain ( fv_domain )
   type(domain2d), intent(out) :: fv_domain
!  returns the domain2d variable associated with the coupling grid
!  note: coupling is done using the mass/temperature grid with no halos

   fv_domain = Atm(mygrid)%domain_for_coupler

 end subroutine atmosphere_domain



 subroutine get_atmosphere_axes ( axes )
   integer, intent(out) :: axes (:)

!----- returns the axis indices for the atmospheric (mass) grid -----
   if ( size(axes(:)) < 0 .or. size(axes(:)) > 4 ) call error_mesg (    &
                               'get_atmosphere_axes in atmosphere_mod', &
                               'size of argument is incorrect', FATAL   )

   axes (1:size(axes(:))) = Atm(mygrid)%atmos_axes (1:size(axes(:)))

 end subroutine get_atmosphere_axes



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
         z_bot(i,j) = rrg*t_bot(i,j)*(1.+zvir*Atm(mygrid)%q(i,j,npz,sphum)) *  &
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
              wm(i,j) = wm(i,j) + Atm(mygrid)%delp(i,j,k) * ( Atm(mygrid)%q(i,j,k,sphum)   +  &
                                                              Atm(mygrid)%q(i,j,k,liq_wat) +  &
                                                              Atm(mygrid)%q(i,j,k,ice_wat) )
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


 subroutine atmosphere_state_update (Time, Physics_tendency, Physics, Atm_block)
   type(time_type),intent(in)      :: Time
   type (physics_tendency_type),   intent(in) :: Physics_tendency
   type (physics_type),    intent(in) :: Physics
   type (block_control_type), intent(in) :: Atm_block
   type(time_type) :: Time_prev, Time_next
!--- local variables ---
   integer :: i, j, k, n, w_diff, nt_dyn
   integer :: nb, ibs, ibe, jbs, jbe
   real ::  rcp

   Time_prev = Time
   Time_next = Time + Time_step_atmos

   n = mygrid

!--- put u/v tendencies into haloed arrays u_dt and v_dt
!$OMP parallel do default(shared) private(nb, ibs, ibe, jbs, jbe)
   do nb = 1,Atm_block%nblks
     ibs = Atm_block%ibs(nb)
     ibe = Atm_block%ibe(nb)
     jbs = Atm_block%jbs(nb)
     jbe = Atm_block%jbe(nb)

     u_dt(ibs:ibe,jbs:jbe,:)   = Physics_tendency%block(nb)%u_dt
     v_dt(ibs:ibe,jbs:jbe,:)   = Physics_tendency%block(nb)%v_dt
     t_dt(ibs:ibe,jbs:jbe,:)   = Physics_tendency%block(nb)%t_dt
     q_dt(ibs:ibe,jbs:jbe,:,:) = Physics_tendency%block(nb)%q_dt

!--- diagnostic tracers are being updated in-place
!--- tracer fields must be returned to the Atm structure
     Atm(mygrid)%qdiag(ibs:ibe,jbs:jbe,:,:) = Physics_tendency%block(nb)%qdiag

   enddo

   w_diff = get_tracer_index (MODEL_ATMOS, 'w_diff' )
   nt_dyn = ncnst-pnats   !nothing more than nq
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

   call mpp_clock_begin (id_dynam)
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
                         Atm(n)%ptop, Atm(n)%phys_diag, Atm(n)%nudge_diag, q_dt)
       call timing_off('FV_UPDATE_PHYS')
   call mpp_clock_end (id_dynam)

!--- nesting update after updating atmospheric variables with
!--- physics tendencies
    if (ngrids > 1 .and. p_split > 0) then
       call timing_on('TWOWAY_UPDATE')
       call twoway_nesting(Atm, ngrids, grids_on_this_pe, zvir, fv_time, mygrid)
       call timing_off('TWOWAY_UPDATE')
    endif

!--- cmip6 total tendencies of temperature and specific humidity
   if (query_cmip_diag_id(ID_tnt)) &
                 used = send_cmip_data_3d ( ID_tnt, (Atm(mygrid)%pt(isc:iec,jsc:jec,:)-ttend(:,:,:))/dt_atmos, Time)
   if (query_cmip_diag_id(ID_tnhus)) &
                  used = send_cmip_data_3d (ID_tnhus, (Atm(mygrid)%q(isc:iec,jsc:jec,:,sphum)-qtend(:,:,:,sphum))/dt_atmos, Time)

  !---- diagnostics for FV dynamics -----
   if (Atm(mygrid)%flagstruct%print_freq /= -99) then
     call mpp_clock_begin(id_fv_diag)
     call timing_on('FV_DIAG')

     fv_time = Time_next
     call get_time (fv_time, seconds,  days)

     call fv_diag(Atm(mygrid:mygrid), zvir, fv_time, Atm(mygrid)%flagstruct%print_freq)
      if (Atm(mygrid)%coarse_graining%write_coarse_diagnostics) then
         call fv_coarse_diag(Atm(mygrid:mygrid), fv_time, zvir)
      endif
     call fv_cmip_diag(Atm(mygrid:mygrid), zvir, fv_time)

     call timing_off('FV_DIAG')
     call mpp_clock_end(id_fv_diag)
   endif

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

     call timing_on('adiabatic_init')
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
                     Atm(mygrid)%flagstruct%q_split, Atm(mygrid)%u, Atm(mygrid)%v, Atm(mygrid)%w,         &
                     Atm(mygrid)%delz, Atm(mygrid)%flagstruct%hydrostatic,                      &
                     Atm(mygrid)%pt, Atm(mygrid)%delp, Atm(mygrid)%q, Atm(mygrid)%ps,                     &
                     Atm(mygrid)%pe, Atm(mygrid)%pk, Atm(mygrid)%peln, Atm(mygrid)%pkz, Atm(mygrid)%phis,      &
                     Atm(mygrid)%q_con, Atm(mygrid)%omga, Atm(mygrid)%ua, Atm(mygrid)%va, Atm(mygrid)%uc, Atm(mygrid)%vc, &
                     Atm(mygrid)%ak, Atm(mygrid)%bk, Atm(mygrid)%mfx, Atm(mygrid)%mfy,                    &
                     Atm(mygrid)%cx, Atm(mygrid)%cy, Atm(mygrid)%ze0, Atm(mygrid)%flagstruct%hybrid_z,    &
                     Atm(mygrid)%gridstruct, Atm(mygrid)%flagstruct,                            &
                     Atm(mygrid)%neststruct, Atm(mygrid)%idiag, Atm(mygrid)%bd, Atm(mygrid)%parent_grid,  &
                     Atm(mygrid)%domain, Atm(mygrid)%inline_mp, Atm(n)%diss_est)
! Backward
    call fv_dynamics(Atm(mygrid)%npx, Atm(mygrid)%npy, npz,  nq, Atm(mygrid)%ng, -dt_atmos, 0.,      &
                     Atm(mygrid)%flagstruct%fill, Atm(mygrid)%flagstruct%reproduce_sum, kappa, cp_air, zvir,  &
                     Atm(mygrid)%ptop, Atm(mygrid)%ks, nq, Atm(mygrid)%flagstruct%n_split,        &
                     Atm(mygrid)%flagstruct%q_split, Atm(mygrid)%u, Atm(mygrid)%v, Atm(mygrid)%w,         &
                     Atm(mygrid)%delz, Atm(mygrid)%flagstruct%hydrostatic,                      &
                     Atm(mygrid)%pt, Atm(mygrid)%delp, Atm(mygrid)%q, Atm(mygrid)%ps,                     &
                     Atm(mygrid)%pe, Atm(mygrid)%pk, Atm(mygrid)%peln, Atm(mygrid)%pkz, Atm(mygrid)%phis,      &
                     Atm(mygrid)%q_con, Atm(mygrid)%omga, Atm(mygrid)%ua, Atm(mygrid)%va, Atm(mygrid)%uc, Atm(mygrid)%vc, &
                     Atm(mygrid)%ak, Atm(mygrid)%bk, Atm(mygrid)%mfx, Atm(mygrid)%mfy,                    &
                     Atm(mygrid)%cx, Atm(mygrid)%cy, Atm(mygrid)%ze0, Atm(mygrid)%flagstruct%hybrid_z,    &
                     Atm(mygrid)%gridstruct, Atm(mygrid)%flagstruct,                            &
                     Atm(mygrid)%neststruct, Atm(mygrid)%idiag, Atm(mygrid)%bd, Atm(mygrid)%parent_grid,  &
                     Atm(mygrid)%domain, Atm(mygrid)%inline_mp, Atm(n)%diss_est)
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
                     Atm(mygrid)%flagstruct%q_split, Atm(mygrid)%u, Atm(mygrid)%v, Atm(mygrid)%w,         &
                     Atm(mygrid)%delz, Atm(mygrid)%flagstruct%hydrostatic,                      &
                     Atm(mygrid)%pt, Atm(mygrid)%delp, Atm(mygrid)%q, Atm(mygrid)%ps,                     &
                     Atm(mygrid)%pe, Atm(mygrid)%pk, Atm(mygrid)%peln, Atm(mygrid)%pkz, Atm(mygrid)%phis,      &
                     Atm(mygrid)%q_con, Atm(mygrid)%omga, Atm(mygrid)%ua, Atm(mygrid)%va, Atm(mygrid)%uc, Atm(mygrid)%vc, &
                     Atm(mygrid)%ak, Atm(mygrid)%bk, Atm(mygrid)%mfx, Atm(mygrid)%mfy,                    &
                     Atm(mygrid)%cx, Atm(mygrid)%cy, Atm(mygrid)%ze0, Atm(mygrid)%flagstruct%hybrid_z,    &
                     Atm(mygrid)%gridstruct, Atm(mygrid)%flagstruct,                            &
                     Atm(mygrid)%neststruct, Atm(mygrid)%idiag, Atm(mygrid)%bd, Atm(mygrid)%parent_grid,  &
                     Atm(mygrid)%domain, Atm(mygrid)%inline_mp, Atm(n)%diss_est)
! Forward call
    call fv_dynamics(Atm(mygrid)%npx, Atm(mygrid)%npy, npz,  nq, Atm(mygrid)%ng, dt_atmos, 0.,      &
                     Atm(mygrid)%flagstruct%fill, Atm(mygrid)%flagstruct%reproduce_sum, kappa, cp_air, zvir,  &
                     Atm(mygrid)%ptop, Atm(mygrid)%ks, nq, Atm(mygrid)%flagstruct%n_split,        &
                     Atm(mygrid)%flagstruct%q_split, Atm(mygrid)%u, Atm(mygrid)%v, Atm(mygrid)%w,         &
                     Atm(mygrid)%delz, Atm(mygrid)%flagstruct%hydrostatic,                      &
                     Atm(mygrid)%pt, Atm(mygrid)%delp, Atm(mygrid)%q, Atm(mygrid)%ps,                     &
                     Atm(mygrid)%pe, Atm(mygrid)%pk, Atm(mygrid)%peln, Atm(mygrid)%pkz, Atm(mygrid)%phis,      &
                     Atm(mygrid)%q_con, Atm(mygrid)%omga, Atm(mygrid)%ua, Atm(mygrid)%va, Atm(mygrid)%uc, Atm(mygrid)%vc, &
                     Atm(mygrid)%ak, Atm(mygrid)%bk, Atm(mygrid)%mfx, Atm(mygrid)%mfy,                    &
                     Atm(mygrid)%cx, Atm(mygrid)%cy, Atm(mygrid)%ze0, Atm(mygrid)%flagstruct%hybrid_z,    &
                     Atm(mygrid)%gridstruct, Atm(mygrid)%flagstruct,                            &
                     Atm(mygrid)%neststruct, Atm(mygrid)%idiag, Atm(mygrid)%bd, Atm(mygrid)%parent_grid,  &
                     Atm(mygrid)%domain, Atm(mygrid)%inline_mp, Atm(n)%diss_est)
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
     call timing_off('adiabatic_init')

 end subroutine adiabatic_init


 subroutine atmos_physics_driver_inputs (Physics, Atm_block, Physics_tendency)
   type (physics_type),  intent(inout) :: Physics
   type (block_control_type), intent(in) :: Atm_block
   type (physics_tendency_type), intent(inout), optional :: Physics_tendency
!--- local variabls
   integer :: nb, ibs, ibe, jbs, jbe

!---------------------------------------------------------------------
! use most up to date atmospheric properties when running serially
!---------------------------------------------------------------------
!$OMP parallel do default(shared) private(nb, ibs, ibe, jbs, jbe)
   do nb = 1, Atm_block%nblks
     ibs = Atm_block%ibs(nb)
     ibe = Atm_block%ibe(nb)
     jbs = Atm_block%jbs(nb)
     jbe = Atm_block%jbe(nb)

     Physics%block(nb)%phis = Atm(mygrid)%phis(ibs:ibe,jbs:jbe)
     Physics%block(nb)%u    = Atm(mygrid)%ua(ibs:ibe,jbs:jbe,:)
     Physics%block(nb)%v    = Atm(mygrid)%va(ibs:ibe,jbs:jbe,:)
     Physics%block(nb)%t    = Atm(mygrid)%pt(ibs:ibe,jbs:jbe,:)
     Physics%block(nb)%q    = Atm(mygrid)%q(ibs:ibe,jbs:jbe,:,:)
     Physics%block(nb)%omega= Atm(mygrid)%omga(ibs:ibe,jbs:jbe,:)
     Physics%block(nb)%pe   = Atm(mygrid)%pe(ibs:ibe,:,jbs:jbe)
     Physics%block(nb)%peln = Atm(mygrid)%peln(ibs:ibe,:,jbs:jbe)
     Physics%block(nb)%delp = Atm(mygrid)%delp(ibs:ibe,jbs:jbe,:)
     if (.not.Physics%control%phys_hydrostatic) then
        Physics%block(nb)%delz = Atm(mygrid)%delz(ibs:ibe,jbs:jbe,:)
        Physics%block(nb)%w    = Atm(mygrid)%w(ibs:ibe,jbs:jbe,:)
     endif
     if (_ALLOCATED(Physics%block(nb)%tmp_4d)) &
        Physics%block(nb)%tmp_4d = Atm(mygrid)%qdiag(ibs:ibe,jbs:jbe,:,:)

     call fv_compute_p_z (Atm_block%npz, Physics%block(nb)%phis, Physics%block(nb)%pe, &
                          Physics%block(nb)%peln, Physics%block(nb)%delp, Physics%block(nb)%delz, &
                          Physics%block(nb)%t, Physics%block(nb)%q(:,:,:,Physics%control%sphum), &
                          Physics%block(nb)%p_full, Physics%block(nb)%p_half, &
                          Physics%block(nb)%z_full, Physics%block(nb)%z_half, &
#ifdef USE_COND
                          Atm(mygrid)%q_con(ibs:ibe,jbs:jbe,:), &
#else
                          Atm(mygrid)%q_con, &
#endif
                          Physics%control%phys_hydrostatic, Physics%control%do_uni_zfull) !miz

     if (PRESENT(Physics_tendency)) then
!--- copy the dynamics tendencies into the physics tendencies
!--- if one wants to run physics concurrent with dynamics,
!--- these values would be zeroed out and accumulated
!--- in the atmosphere_state_update

       Physics_tendency%block(nb)%u_dt = u_dt(ibs:ibe,jbs:jbe,:)
       Physics_tendency%block(nb)%v_dt = v_dt(ibs:ibe,jbs:jbe,:)
       Physics_tendency%block(nb)%t_dt = t_dt(ibs:ibe,jbs:jbe,:)
       Physics_tendency%block(nb)%q_dt = q_dt(ibs:ibe,jbs:jbe,:,:)
       Physics_tendency%block(nb)%qdiag = Atm(mygrid)%qdiag(ibs:ibe,jbs:jbe,:,:)
     endif
   enddo

 end subroutine atmos_physics_driver_inputs


 subroutine atmos_radiation_driver_inputs (Time, Radiation, Atm_block)
   type (time_type),          intent(in)    :: Time
   type (radiation_type),     intent(inout) :: Radiation
   type (block_control_type), intent(in)    :: Atm_block
!--- local variables
   integer :: nb, ibs, ibe, jbs, jbe

!---------------------------------------------------------------------
! use most up to date atmospheric properties when running serially
!---------------------------------------------------------------------
!$OMP parallel do default(shared) private(nb, ibs, ibe, jbs, jbe)
   do nb = 1,Atm_block%nblks
     ibs = Atm_block%ibs(nb)
     ibe = Atm_block%ibe(nb)
     jbs = Atm_block%jbs(nb)
     jbe = Atm_block%jbe(nb)

     Radiation%block(nb)%phis = Atm(mygrid)%phis(ibs:ibe,jbs:jbe)
     Radiation%block(nb)%t    = Atm(mygrid)%pt(ibs:ibe,jbs:jbe,:)
     Radiation%block(nb)%q    = Atm(mygrid)%q(ibs:ibe,jbs:jbe,:,:)
     Radiation%block(nb)%pe   = Atm(mygrid)%pe(ibs:ibe,:,jbs:jbe)
     Radiation%block(nb)%peln = Atm(mygrid)%peln(ibs:ibe,:,jbs:jbe)
     Radiation%block(nb)%delp = Atm(mygrid)%delp(ibs:ibe,jbs:jbe,:)
     if (.not.Radiation%control%phys_hydrostatic) &
        Radiation%block(nb)%delz = Atm(mygrid)%delz(ibs:ibe,jbs:jbe,:)

     call fv_compute_p_z (Atm_block%npz, Radiation%block(nb)%phis, Radiation%block(nb)%pe, &
                          Radiation%block(nb)%peln, Radiation%block(nb)%delp, Radiation%block(nb)%delz, &
                          Radiation%block(nb)%t, Radiation%block(nb)%q(:,:,:,Radiation%control%sphum), &
                          Radiation%block(nb)%p_full, Radiation%block(nb)%p_half, &
                          Radiation%block(nb)%z_full, Radiation%block(nb)%z_half, &
#ifdef USE_COND
                          Atm(mygrid)%q_con(ibs:ibe,jbs:jbe,:), &
#else
                          Atm(mygrid)%q_con, &
#endif
                          Radiation%control%phys_hydrostatic, Radiation%control%do_uni_zfull) !miz
   enddo

!----------------------------------------------------------------------
! obtain pressure-weighted global mean co2 dry volume mixing ratio for
! use by radiation package.
!----------------------------------------------------------------------
! compute_g_avg must be called here because it contains
! mpp_sums that cannot be called during the concurrent radiation
! phase due to the way in which MPI interacts with nested OpenMP
!----------------------------------------------------------------------
   call compute_g_avg(Time, 'co2', Radiation, Atm_block)
   call compute_g_avg(Time, 'ch4', Radiation, Atm_block)

 end subroutine atmos_radiation_driver_inputs



 subroutine fv_compute_p_z (npz, phis, pe, peln, delp, delz, pt, q_sph, &
                            p_full, p_half, z_full, z_half, q_con, hydrostatic, do_uni_zfull) !miz
    integer, intent(in)  :: npz
    real, dimension(:,:),   intent(in)  :: phis
    real, dimension(:,:,:), intent(in)  :: pe, peln, delp, delz, q_con, pt, q_sph
    real, dimension(:,:,:), intent(out) :: p_full, p_half, z_full, z_half
    logical, intent(in)  :: hydrostatic, do_uni_zfull !miz
!--- local variables
    integer i,j,k,isiz,jsiz
    real    tvm
    real    :: zvir, rrg, ginv
#ifdef USE_COND
    real, dimension(size(pe,1),size(pe,3),size(pe,2)):: peg, pelng
    real:: dlg
#endif

    isiz=size(phis,1)
    jsiz=size(phis,2)
    zvir = rvgas/rdgas - 1.
    ginv = 1./ grav
    rrg  = rdgas / grav

!----------------------------------------------------
! Compute pressure and height at full and half levels
!----------------------------------------------------
    z_half(:,:,npz+1) = phis(:,:) * ginv

    do k=1,npz+1
      do j=1,jsiz
         do i=1,isiz
           p_half(i,j,k) = pe(i,k,j)
        enddo
      enddo
    enddo

!--------- Hydrostatic option ----------------------------------------------
    if (hydrostatic ) then
#ifdef USE_COND
    do j=1,jsiz
       do i=1,isiz
          peg(i,j,1) = pe(i,1,j)
       enddo
    end do
    do k=2,npz+1
       do j=1,jsiz
          do i=1,isiz
             peg(i,j,k) = peg(i,j,k-1) + delp(i,j,k-1)*(1.-q_con(i,j,k-1))
          enddo
       enddo
    enddo
#endif
      do k=npz,1,-1
        do j=1,jsiz
          do i=1,isiz
            tvm = rrg*pt(i,j,k)*(1.+zvir*q_sph(i,j,k))
            p_full(i,j,k) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
#ifdef USE_COND
            dlg = log(peg(i,j,k+1)/peg(i,j,k))
            z_full(i,j,k) = z_half(i,j,k+1) + tvm*(1.-peg(i,j,k)*dlg/(peg(i,j,k+1)-peg(i,j,k)))
            z_half(i,j,k) = z_half(i,j,k+1) + tvm*dlg
#else
            z_full(i,j,k) = z_half(i,j,k+1) + tvm*(1.-p_half(i,j,k)/p_full(i,j,k))
            z_half(i,j,k) = z_half(i,j,k+1) + tvm*(peln(i,k+1,j)-peln(i,k,j))
#endif
          enddo
        enddo
      enddo
    else
!--------- Non-Hydrostatic option ------------------------------------------
      do k=npz,1,-1
        do j=1,jsiz
          do i=1,isiz
            p_full(i,j,k) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
            z_half(i,j,k) = z_half(i,j,k+1) - delz(i,j,k)
            z_full(i,j,k) = 0.5*(z_half(i,j,k) + z_half(i,j,k+1))
          enddo
        enddo
      enddo
    endif
    if (do_uni_zfull) then
       do k=1,npz
         z_full(:,:,k)=0.5*(z_half(:,:,k)+z_half(:,:,k+1))
       enddo
    endif
  end subroutine fv_compute_p_z


 subroutine reset_atmos_tracers (Physics, Physics_tendency, Atm_block)
   type (physics_type), intent(in) :: Physics
   type (physics_tendency_type), intent(in) :: Physics_tendency
   type (block_control_type), intent(in) :: Atm_block
!--- local variables
   integer :: nb, ibs, ibe, jbs, jbe

!--- After initialization by the physics, tracer fields must be
!--- returned to the Atm structure.  This is because tracer driver
!--- can reset the initial values
!$OMP parallel do default(shared) private(nb, ibs, ibe, jbs, jbe)
    do nb = 1, Atm_block%nblks
      ibs = Atm_block%ibs(nb)
      ibe = Atm_block%ibe(nb)
      jbs = Atm_block%jbs(nb)
      jbe = Atm_block%jbe(nb)

      Atm(mygrid)%q(ibs:ibe,jbs:jbe,:,:) = Physics%block(nb)%q
      Atm(mygrid)%qdiag(ibs:ibe,jbs:jbe,:,:) = Physics_tendency%block(nb)%qdiag
    enddo

 end subroutine reset_atmos_tracers

end module atmosphere_mod
