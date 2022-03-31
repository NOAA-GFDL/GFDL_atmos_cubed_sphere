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
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

!***********************************************************************
!> @file
!! @brief Provides top-level interface for moving nest functionality
!! @author W. Ramstrom, AOML/HRD   05/27/2021
!! @email William.Ramstrom@noaa.gov
! =======================================================================!

module fv_moving_nest_main_mod
#ifdef MOVING_NEST

#include <fms_platform.h>

  !-----------------
  ! FMS modules:
  !-----------------
  use block_control_mod,      only: block_control_type
  use constants_mod,          only: cp_air, rdgas, grav, rvgas, kappa, pstd_mks
  use time_manager_mod,       only: time_type, get_time, get_date, set_time, operator(+), &
      operator(-), operator(/), time_type_to_real
  use fms_mod,                only: file_exist, open_namelist_file,    &
      close_file, error_mesg, FATAL,     &
      check_nml_error, stdlog,           &
      write_version_number,              &
      mpp_clock_id, mpp_clock_begin,     &
      mpp_clock_end, CLOCK_SUBCOMPONENT, &
      clock_flag_default
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
  use DYCORE_typedefs,        only: DYCORE_data_type
#ifdef GFS_TYPES
  use GFS_typedefs,           only: IPD_data_type => GFS_data_type, &
      IPD_control_type => GFS_control_type, kind_phys
#else
  use IPD_typedefs,           only: IPD_data_type, IPD_control_type, kind_phys => IPD_kind_phys
#endif

  use fv_iau_mod,             only: IAU_external_data_type
#ifdef MULTI_GASES
  use multi_gases_mod,  only: virq, virq_max, num_gas, ri, cpi
#endif

  !-----------------
  ! FV core modules:
  !-----------------
  use atmosphere_mod,     only: Atm, mygrid, p_split, dt_atmos
  use fv_arrays_mod,      only: fv_atmos_type, R_GRID, fv_grid_bounds_type, phys_diag_type
  use fv_moving_nest_types_mod, only: allocate_fv_moving_nest_prog_type, allocate_fv_moving_nest_physics_type
  use fv_moving_nest_types_mod, only: Moving_nest
  use fv_diagnostics_mod, only: fv_diag_init, fv_diag_reinit, fv_diag, fv_time, prt_maxmin, prt_height
  use fv_restart_mod,     only: fv_restart, fv_write_restart
  use fv_timing_mod,      only: timing_on, timing_off
  use fv_mp_mod,          only: is_master
  use fv_regional_mod,    only: start_regional_restart, read_new_bc_data, a_step, p_step, current_time_in_seconds

  !-----------------------------------------
  !  External routines
  !-----------------------------------------
  use fms_io_mod,         only: fms_io_exit
  use mpp_domains_mod,    only: NORTH, NORTH_EAST, EAST, SOUTH_EAST, CORNER, CENTER
  use mpp_domains_mod,    only: nest_domain_type
  use mpp_mod,            only: mpp_sync, mpp_exit
  use mpp_domains_mod,    only: mpp_get_global_domain
  use mpp_mod,            only: mpp_send, mpp_sync_self, mpp_broadcast

  use fv_mp_mod,          only: global_nest_domain

  use tracer_manager_mod, only: get_tracer_names
  use field_manager_mod,  only: MODEL_ATMOS
  use fv_io_mod,          only: fv_io_exit
  !!use fv_restart_mod,     only: d2c_setup

  !------------------------------------
  !  Moving Nest Routines
  !------------------------------------

  !      Prognostic variable routines
  use fv_moving_nest_mod,         only: mn_prog_fill_intern_nest_halos, mn_prog_fill_nest_halos_from_parent, &
      mn_prog_dump_to_netcdf, mn_prog_shift_data
  !      Physics variable routines
  use fv_moving_nest_physics_mod, only: mn_phys_fill_intern_nest_halos, mn_phys_fill_nest_halos_from_parent, &
      mn_phys_dump_to_netcdf, mn_phys_shift_data, mn_phys_reset_sfc_props, move_nsst

  !      Metadata routines
  use fv_moving_nest_mod,         only: mn_meta_move_nest, mn_meta_recalc, mn_meta_reset_gridstruct, mn_shift_index

  !      Temporary variable routines (delz)
  use fv_moving_nest_mod,         only: mn_prog_fill_temp_variables, mn_prog_apply_temp_variables
  use fv_moving_nest_physics_mod, only: mn_phys_fill_temp_variables, mn_phys_apply_temp_variables

  !      Load static datasets
  use fv_moving_nest_mod,         only: mn_latlon_read_hires_parent, mn_latlon_load_parent
  use fv_moving_nest_mod,         only: mn_orog_read_hires_parent, mn_static_read_hires
  use fv_moving_nest_utils_mod,   only: load_nest_latlons_from_nc, compare_terrain, set_smooth_nest_terrain, set_blended_terrain

  use fv_moving_nest_physics_mod, only: mn_reset_phys_latlon, mn_surface_grids

  !      Grid reset routines
  use fv_moving_nest_mod,         only: grid_geometry, assign_n_p_grids, move_nest_geo
  use fv_moving_nest_utils_mod,   only: fill_grid_from_supergrid, fill_weight_grid

  !      Physics moving logical variables
  use fv_moving_nest_physics_mod, only: move_physics, move_nsst

  !      Recalculation routines
  use fv_moving_nest_mod,         only: reallocate_BC_buffers, recalc_aux_pressures

  !      Logging and debugging information
  use fv_moving_nest_mod,         only: check_array
  use fv_moving_nest_utils_mod,   only: show_atm, show_atm_grids, show_tile_geo, show_nest_grid, show_gridstruct, grid_equal
  use fv_moving_nest_utils_mod,   only: validate_hires_parent

  use fv_tracker_mod,             only: Tracker, allocate_tracker

  implicit none

  !-----------------------------------------------------------------------
  ! version number of this module
  ! Include variable "version" to be written to log file.
#include<file_version.h>
  character(len=20)   :: mod_name = 'fvGFS/fv_moving_nest_main_mod'

#ifdef OVERLOAD_R4
  real, parameter:: real_snan=x'FFBFFFFF'
#else
  real, parameter:: real_snan=x'FFF7FFFFFFFFFFFF'
#endif

  ! Enable these for more debugging outputs
  logical :: debug_log = .false.    ! Produces logging to out.* file
  logical :: tsvar_out = .false.    ! Produces netCDF outputs; be careful to not exceed file number limits set in namelist

  !  --- Clock ids for moving_nest performance metering
  integer :: id_movnest1, id_movnest1_9, id_movnest2, id_movnest3, id_movnest4, id_movnest5
  integer :: id_movnest6, id_movnest7_0, id_movnest7_1, id_movnest7_2, id_movnest7_3, id_movnest8, id_movnest9
  integer :: id_movnestTot
  logical :: use_timers = .False. ! Set this to true for detailed performance profiling.  False only profiles total moving nest time.
  integer, save :: output_step = 0

contains

  !>@brief The subroutine 'update_moving_nest' decides whether the nest should be moved, and if so, performs the move.
  !>@details This subroutine evaluates the automatic storm tracker (or prescribed motion configuration), then decides
  !!  if the nest should be moved.  If it should be moved, it calls fv_moving_nest_exec() to perform the nest move.
  subroutine update_moving_nest(Atm_block, IPD_control, IPD_data, time_step)
    type(block_control_type), intent(in) :: Atm_block     !< Physics block layout
    type(IPD_control_type), intent(in)   :: IPD_control   !< Physics metadata
    type(IPD_data_type), intent(inout)   :: IPD_data(:)   !< Physics variable data
    type(time_type), intent(in)          :: time_step     !< Current timestep

    logical :: do_move
    integer :: delta_i_c, delta_j_c
    integer :: parent_grid_num, child_grid_num, nest_num
    integer, allocatable :: global_pelist(:)
    integer :: n
    integer :: this_pe

    this_pe = mpp_pe()

    do_move = .false.

    ! dt_atmos was initialized in atmosphere.F90::atmosphere_init()

    n = mygrid   ! Public variable from atmosphere.F90

    ! Hard-coded for now - these will need to be looked up on each PE when multiple and telescoped nests are enabled.
    parent_grid_num = 1
    child_grid_num = 2
    nest_num = 1

    call eval_move_nest(Atm, a_step, parent_grid_num, child_grid_num, do_move, delta_i_c, delta_j_c, dt_atmos)

    allocate(global_pelist(Atm(parent_grid_num)%npes_this_grid+Atm(child_grid_num)%npes_this_grid))
    global_pelist=(/Atm(parent_grid_num)%pelist, Atm(child_grid_num)%pelist/)

    call mpp_set_current_pelist(global_pelist)
    call mpp_broadcast( delta_i_c, Atm(child_grid_num)%pelist(1), global_pelist )
    call mpp_broadcast( delta_j_c, Atm(child_grid_num)%pelist(1), global_pelist )
    call mpp_broadcast( do_move, Atm(child_grid_num)%pelist(1), global_pelist )
    call mpp_set_current_pelist(Atm(n)%pelist)

    if (do_move) then
      call fv_moving_nest_exec(Atm, Atm_block, IPD_control, IPD_data, delta_i_c, delta_j_c, n, nest_num, parent_grid_num, child_grid_num, dt_atmos)
    endif

  end subroutine update_moving_nest

  !>@brief The subroutine 'dump_moving_nest' outputs native grid format data to netCDF files
  !>@details This subroutine exports model variables using FMS IO to netCDF files if tsvar_out is set to .True.
  subroutine dump_moving_nest(Atm_block, IPD_control, IPD_data, time_step)
    type(block_control_type), intent(in) :: Atm_block     !< Physics block layout
    type(IPD_control_type), intent(in)   :: IPD_control   !< Physics metadata
    type(IPD_data_type), intent(in)      :: IPD_data(:)   !< Physics variable data
    type(time_type), intent(in)          :: time_step     !< Current timestep

    type(domain2d), pointer              :: domain_coarse, domain_fine
    logical                              :: is_fine_pe
    integer :: parent_grid_num, child_grid_num, nz, this_pe, n

    this_pe = mpp_pe()
    n = mygrid

    parent_grid_num = 1
    child_grid_num = 2

    domain_fine => Atm(child_grid_num)%domain
    domain_coarse => Atm(parent_grid_num)%domain
    is_fine_pe = Atm(n)%neststruct%nested .and. ANY(Atm(n)%pelist(:) == this_pe)
    nz = Atm(n)%npz

    if (debug_log) print '("[INFO] WDR ptbounds 3 atmosphere.F90  npe=",I0," pt(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', this_pe, lbound(Atm(n)%pt,1), ubound(Atm(n)%pt,1), lbound(Atm(n)%pt,2), ubound(Atm(n)%pt,2), lbound(Atm(n)%pt,3), ubound(Atm(n)%pt,3)

    ! Enable this to dump debug netCDF files.  Make sure to enable fms_io_exit() in fv_control.F90 so that files are written and closed.
    !if (mod(a_step, 20) .eq. 0 ) then
    !  if (tsvar_out) call mn_prog_dump_to_netcdf(Atm(n), a_step, "tsavar", is_fine_pe, domain_coarse, domain_fine, nz)
    !  if (tsvar_out) call mn_phys_dump_to_netcdf(Atm(n), Atm_block, IPD_control, IPD_data, a_step, "tsavar", is_fine_pe, domain_coarse, domain_fine, nz)
    !endif

  end subroutine dump_moving_nest

  !>@brief The subroutine 'fv_moving_nest_init_clocks' intializes performance profiling timers of sections of the moving nest code.
  !>@details Starts timers for subcomponents of moving nest code to determine performance.  mpp routines group them into separate
  !! sections for parent and nest PEs.
  subroutine fv_moving_nest_init_clocks()

    !  --- initialize clocks for moving_nest
    if (use_timers) then
      id_movnest1     = mpp_clock_id ('MN Part 1 Init',  flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
      id_movnest1_9   = mpp_clock_id ('MN Part 1.9 Copy delz',  flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
      id_movnest2     = mpp_clock_id ('MN Part 2 Fill Halos from Parent',  flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
      id_movnest3     = mpp_clock_id ('MN Part 3 Meta Move Nest',  flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
      id_movnest4     = mpp_clock_id ('MN Part 4 Fill Intern Nest Halos',  flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
      id_movnest5     = mpp_clock_id ('MN Part 5 Recalc Weights',  flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
      id_movnest6     = mpp_clock_id ('MN Part 6 EOSHIFT',  flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )

      id_movnest7_0   = mpp_clock_id ('MN Part 7.0 Recalc gridstruct',  flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
      id_movnest7_1   = mpp_clock_id ('MN Part 7.1 Refill halos from Parent',  flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
      id_movnest7_2   = mpp_clock_id ('MN Part 7.2 Refill Intern Nest Halos',  flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
      id_movnest7_3   = mpp_clock_id ('MN Part 7.3 Fill delz',  flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )

      id_movnest8     = mpp_clock_id ('MN Part 8 Dump to netCDF',  flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
      id_movnest9     = mpp_clock_id ('MN Part 9 Aux Pressure',  flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
    endif

    id_movnestTot     = mpp_clock_id ('Moving Nest Total',  flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  end subroutine fv_moving_nest_init_clocks

  !>@brief The subroutine 'eval_move_nest' determines whether the nest should be moved and in which direction.
  !>@details  This subroutine can execute prescribed motion or automated storm tracking based on namelist options.
  subroutine eval_move_nest(Atm, a_step, parent_grid_num, child_grid_num, do_move, delta_i_c, delta_j_c, dt_atmos)
    type(fv_atmos_type), intent(inout)   :: Atm(:)       !< Input atmospheric data
    integer, intent(in)                  :: a_step       !< Timestep
    integer, intent(in)                  :: parent_grid_num, child_grid_num  !< Grid numbers of parent and child
    logical, intent(out)                 :: do_move      !< Logical for whether to move nest
    integer, intent(out)                 :: delta_i_c, delta_j_c  !< Each can be -1, 0, or +1
    real, intent(in)                     :: dt_atmos     !< only needed for the simple version of this subroutine

    integer       :: n
    integer       :: cx, cy
    real          :: xdiff, ydiff
    integer       :: nest_i_c, nest_j_c
    integer       :: nis, nie, njs, nje
    integer       :: this_pe
    character*255 :: message

    ! On the tropical channel configuration, tile 6 numbering starts at 0,0 off the coast of Spain
    !  delta_i_c = +1 is westward
    !  delta_i_c = -1 is eastward
    !
    !  delta_j_c = +1 is southward
    !  delta_j_c = -1 is northward

    this_pe = mpp_pe()
    n = mygrid   ! Public variable from atmosphere.F90
    do_move = .false.
    delta_i_c = 0
    delta_j_c = 0

    if ( Moving_nest(n)%mn_flag%vortex_tracker .eq. 0  .or. Atm(n)%grid_number .eq. 1) then
      ! No need to move
      do_move = .false.
      delta_i_c = 0
      delta_j_c = 0
    else if ( Moving_nest(n)%mn_flag%vortex_tracker .eq. 1 ) then
      ! Prescribed move according to ntrack, move_cd_x and move_cd_y
      ! Move every ntrack of dt_atmos time step
      if ( mod(a_step,Moving_nest(n)%mn_flag%ntrack) .eq. 0) then
        do_move = .true.
        delta_i_c = Moving_nest(n)%mn_flag%move_cd_x
        delta_j_c = Moving_nest(n)%mn_flag%move_cd_y
      endif
    else if ( Moving_nest(n)%mn_flag%vortex_tracker .eq. 2 .or. &
        Moving_nest(n)%mn_flag%vortex_tracker .eq. 6 .or. &
        Moving_nest(n)%mn_flag%vortex_tracker .eq. 7 ) then
      ! Automatic moving following the internal storm tracker
      if ( mod(a_step,Moving_nest(n)%mn_flag%ntrack) .eq. 0) then
        if(Tracker(n)%tracker_gave_up) then
          call mpp_error(NOTE,'Not moving: tracker decided the storm dissapated')
          return
        endif
        if(.not.Tracker(n)%tracker_havefix) then
          call mpp_error(NOTE,'Not moving: tracker did not find a storm')
          return
        endif
        ! Calcuate domain center indexes
        cx=(Atm(n)%npx-1)/2+1
        cy=(Atm(n)%npy-1)/2+1
        ! Calculate distance in parent grid index space between storm
        ! center and domain center
        ! Consider using xydiff as integers in the future?
        xdiff=(Tracker(n)%tracker_ifix-real(cx))/Atm(n)%neststruct%refinement
        ydiff=(Tracker(n)%tracker_jfix-real(cy))/Atm(n)%neststruct%refinement
        if(xdiff .ge. 1.0) then
          Moving_nest(n)%mn_flag%move_cd_x=1
        else if(xdiff .le. -1.0) then
          Moving_nest(n)%mn_flag%move_cd_x=-1
        else
          Moving_nest(n)%mn_flag%move_cd_x=0
        endif
        if(ydiff .ge. 1.0) then
          Moving_nest(n)%mn_flag%move_cd_y=1
        else if(ydiff .le. -1.0) then
          Moving_nest(n)%mn_flag%move_cd_y=-1
        else
          Moving_nest(n)%mn_flag%move_cd_y=0
        endif
        if(abs(Moving_nest(n)%mn_flag%move_cd_x)>0 .or. abs(Moving_nest(n)%mn_flag%move_cd_y)>0) then
          call mpp_error(NOTE,'Moving: tracker center shifted from nest center')
          do_move = .true.
          delta_i_c = Moving_nest(n)%mn_flag%move_cd_x
          delta_j_c = Moving_nest(n)%mn_flag%move_cd_y
        else
          call mpp_error(NOTE,'Not moving: tracker center is near nest center')
          do_move = .false.
          delta_i_c = 0
          delta_j_c = 0
        endif
      endif
    else
      write(message,*) 'Wrong vortex_tracker option: ', Moving_nest(n)%mn_flag%vortex_tracker
      call mpp_error(FATAL,message)
    endif

    ! Override to prevent move on first timestep
    if (a_step .eq. 0) then
      do_move = .false.
      delta_i_c = 0
      delta_j_c = 0
    endif

    ! Check whether or not the nest move is permitted
    if (n==child_grid_num) then
      !  Figure out the bounds of the cube face

      ! x parent bounds: 1 to Atm(parent_grid_num)%flagstruct%npx
      ! y parent bounds: 1 to Atm(parent_grid_num)%flagstruct%npy

      !  Figure out the bounds of the nest

      ! x nest bounds: 1 to Atm(child_grid_num)%flagstruct%npx
      ! y nest bounds: 1 to Atm(child_grid_num)%flagstruct%npy

      ! Nest refinement: Atm(child_grid_num)%neststruct%refinement
      ! Nest starting cell in x direction:  Atm(child_grid_num)%neststruct%ioffset
      ! Nest starting cell in y direction:  Atm(child_grid_num)%neststruct%joffset

      nest_i_c = ( Atm(child_grid_num)%flagstruct%npx - 1 ) / Atm(child_grid_num)%neststruct%refinement
      nest_j_c = ( Atm(child_grid_num)%flagstruct%npy - 1 ) / Atm(child_grid_num)%neststruct%refinement

      nis = Atm(child_grid_num)%neststruct%ioffset + delta_i_c
      nie = Atm(child_grid_num)%neststruct%ioffset + nest_i_c + delta_i_c

      njs = Atm(child_grid_num)%neststruct%joffset + delta_j_c
      nje = Atm(child_grid_num)%neststruct%joffset + nest_j_c + delta_j_c

      !  Will the nest motion push the nest over one of the edges?
      !  Handle each direction individually, so that nest could slide along edge

      ! Causes a crash if we use .le. 1
      if (nis .le. Moving_nest(child_grid_num)%mn_flag%corral_x) then
        delta_i_c = 0
        !      block_moves = .true.
        write(message,*) 'eval_move_nest motion in x direction blocked.  small nis: ', nis
        call mpp_error(WARNING,message)
      endif
      if (njs .le. Moving_nest(child_grid_num)%mn_flag%corral_y) then
        delta_j_c = 0
        !      block_moves = .true.
        write(message,*) 'eval_move_nest motion in y direction blocked.  small njs: ', njs
        call mpp_error(WARNING,message)
      endif

      if (nie .ge. Atm(parent_grid_num)%flagstruct%npx - Moving_nest(child_grid_num)%mn_flag%corral_x) then
        delta_i_c = 0
        !      block_moves = .true.
        write(message,*) 'eval_move_nest motion in x direction blocked.  large nie: ', nie
        call mpp_error(WARNING,message)
      endif
      if (nje .ge. Atm(parent_grid_num)%flagstruct%npy - Moving_nest(child_grid_num)%mn_flag%corral_y) then
        delta_j_c = 0
        !      block_moves = .true.
        write(message,*) 'eval_move_nest motion in y direction blocked.  large nje: ', nje
        call mpp_error(WARNING,message)
      endif

      if (delta_i_c .eq. 0 .and. delta_j_c .eq. 0) then
        do_move = .false.
      endif

    endif

    write(message, *) 'eval_move_nest: move_cd_x=', delta_i_c, 'move_cd_y=', delta_j_c, 'do_move=', do_move
    call mpp_error(NOTE,message)

  end subroutine eval_move_nest

  !>@brief The subroutine 'fv_moving_nest_exec' performs the nest move - most work occurs on nest PEs but some on parent PEs.
  !>@details This subroutine shifts the prognostic and physics/surface variables.
  !!  It also updates metadata and interpolation weights.
  subroutine fv_moving_nest_exec(Atm, Atm_block, IPD_control, IPD_data, delta_i_c, delta_j_c, n, nest_num, parent_grid_num, child_grid_num, dt_atmos)
    implicit none
    type(fv_atmos_type), allocatable, target, intent(inout) :: Atm(:)                !< Atmospheric variables
    type(block_control_type), intent(in)                    :: Atm_block             !< Physics block
    type(IPD_control_type), intent(in)                      :: IPD_control           !< Physics metadata
    type(IPD_data_type), intent(inout)                      :: IPD_data(:)           !< Physics variable data
    integer, intent(in)                                     :: delta_i_c, delta_j_c  !< Nest motion increments
    integer, intent(in)                                     :: n, nest_num           !< Nest indices
    integer, intent(in)                                     :: parent_grid_num, child_grid_num  !< Grid numbers
    real, intent(in)                                        :: dt_atmos              !< Timestep in seconds

    !---- Moving Nest local variables  -----
    integer                                        :: this_pe
    integer, pointer                               :: ioffset, joffset
    real, pointer, dimension(:,:,:)                :: grid, agrid
    type(domain2d), pointer                        :: domain_coarse, domain_fine
    real(kind=R_GRID), pointer, dimension(:,:,:,:) :: grid_global

    ! Constants for mpp calls
    integer  :: position      = CENTER
    integer  :: position_u    = NORTH
    integer  :: position_v    = EAST
    logical  :: do_move = .True.
    integer  :: x_refine, y_refine  ! Currently equal, but allows for future flexibility
    logical  :: is_fine_pe

    ! TODO read halo size from the namelist instead to allow nest refinement > 3
    integer  :: ehalo = 3
    integer  :: whalo = 3
    integer  :: nhalo = 3
    integer  :: shalo = 3
    integer  :: extra_halo = 0   ! Extra halo for moving nest routines

    integer  :: istart_fine, iend_fine, jstart_fine, jend_fine
    integer  :: istart_coarse, iend_coarse, jstart_coarse, jend_coarse
    integer  :: nx, ny, nz, nx_cubic, ny_cubic
    integer  :: p_istart_fine, p_iend_fine, p_jstart_fine, p_jend_fine

    ! Parent tile data, saved between timesteps
    logical, save                          :: first_nest_move = .true.
    type(grid_geometry), save              :: parent_geo
    type(grid_geometry), save              :: fp_super_tile_geo
    type(mn_surface_grids), save           :: mn_static

    type(grid_geometry)              :: tile_geo, tile_geo_u, tile_geo_v
    real(kind=R_GRID), allocatable   :: p_grid(:,:,:), n_grid(:,:,:)
    real(kind=R_GRID), allocatable   :: p_grid_u(:,:,:), n_grid_u(:,:,:)
    real(kind=R_GRID), allocatable   :: p_grid_v(:,:,:), n_grid_v(:,:,:)
    real, allocatable  :: wt_h(:,:,:)
    real, allocatable  :: wt_u(:,:,:)
    real, allocatable  :: wt_v(:,:,:)
    !real :: ua(isd:ied,jsd:jed)
    !real :: va(isd:ied,jsd:jed)

    logical :: filtered_terrain = .True.   ! TODO set this from namelist
    integer :: i, j, x, y, z, p, nn, n_moist
    integer :: parent_tile
    logical :: found_nest_domain = .false.

    ! Variables to enable debugging use of mpp_sync
    logical              :: debug_sync = .false.
    integer, allocatable :: full_pelist(:)
    integer              :: pp, p1, p2

    ! Variables for parent side of setup_aligned_nest()
    integer :: isg, ieg, jsg, jeg, gid
    integer :: isc_p, iec_p, jsc_p, jec_p
    integer :: upoff, jind
    integer :: ng, refinement
    integer :: npx, npy, npz, ncnst, pnats
    integer :: isc, iec, jsc, jec
    integer :: isd, ied, jsd, jed
    integer :: nq                       !  number of transported tracers
    integer :: is, ie, js, je, k        ! For recalculation of omga
    integer, save :: output_step = 0
    integer, allocatable :: pelist(:)
    character(len=16) :: errstring
    logical :: is_moving_nest  !! TODO Refine this per Atm(n) structure to allow some static and some moving nests in same run
    integer             :: year, month, day, hour, minute, second
    real(kind=R_GRID)   :: pi = 4 * atan(1.0d0)
    real                :: rad2deg

    rad2deg = 180.0 / pi

    gid = mpp_pe()
    this_pe = mpp_pe()

    allocate(pelist(mpp_npes()))
    call mpp_get_current_pelist(pelist)

    ! Get month to use for reading static datasets
    call get_date(Atm(n)%Time_init, year, month, day, hour, minute, second)

    ! mygrid and n are the same in atmosphere.F90
    npx   = Atm(n)%npx
    npy   = Atm(n)%npy
    npz   = Atm(n)%npz
    ncnst = Atm(n)%ncnst
    pnats = Atm(n)%flagstruct%pnats

    isc = Atm(n)%bd%isc
    iec = Atm(n)%bd%iec
    jsc = Atm(n)%bd%jsc
    jec = Atm(n)%bd%jec

    isd = isc - Atm(n)%bd%ng
    ied = iec + Atm(n)%bd%ng
    jsd = jsc - Atm(n)%bd%ng
    jed = jec + Atm(n)%bd%ng

    is = Atm(n)%bd%is
    ie = Atm(n)%bd%ie
    js = Atm(n)%bd%js
    je = Atm(n)%bd%je

    nq = ncnst-pnats

    is_fine_pe = Atm(n)%neststruct%nested .and. ANY(Atm(n)%pelist(:) == this_pe)


    if (first_nest_move) then
      if (debug_log) print '("[INFO] WDR Start Clocks npe=",I0," n=",I0)', this_pe, n
      call fv_moving_nest_init_clocks()

      ! If NSST is turned off, do not move the NSST variables.
      !  Namelist switches are confusing; this should be the correct way to distinguish, not using nst_anl
      if (IPD_Control%nstf_name(1) == 0) then
        move_nsst=.false.
      else
        move_nsst=.true.
      endif

      ! This will only allocate the mn_prog and mn_phys for the active Atm(n), not all of them
      !  The others can safely remain unallocated.
      if (debug_log) print '("[INFO] WDR call allocate_fv_moving_nest_prog npe=",I0," n=",I0)', this_pe, n
      call allocate_fv_moving_nest_prog_type(isd, ied, jsd, jed, npz, Moving_nest(n)%mn_prog)
      call allocate_fv_moving_nest_physics_type(isd, ied, jsd, jed, npz, move_physics, move_nsst, &
          IPD_Control%lsoil, IPD_Control%nmtvr, IPD_Control%levs, IPD_Control%ntot2d, IPD_Control%ntot3d, &
          Moving_nest(n)%mn_phys)

    endif

    !==================================================================================================
    !
    !  Begin moving nest code
    !      W. Ramstrom - AOML/HRD/CIMAS 01/15/2021
    !
    !==================================================================================================

    !!================================================================
    !! Step 1 -- Initialization
    !!================================================================

    if (debug_log) print '("WDR_NEST_HALO_RECV,",I0,"===STEP 1====")', this_pe
    if (debug_log) print '("[INFO] WDR MV_NST1 run step 1 fv_moving_nest_main.F90 npe=",I0)', this_pe

    domain_fine => Atm(child_grid_num)%domain
    parent_tile = Atm(child_grid_num)%neststruct%parent_tile
    domain_coarse => Atm(parent_grid_num)%domain
    is_moving_nest = Moving_nest(child_grid_num)%mn_flag%is_moving_nest
    nz = Atm(n)%npz

    if (debug_log) then
      if (is_fine_pe) then
        print '("[INFO] WDR move_nest FINE. npe=",I0, " ", I2.2," do_move=",L1," delta_i_c=",I0," delta_j_c=",I0)', this_pe, n, do_move, delta_i_c, delta_j_c
      else
        print '("[INFO] WDR move_nest COARSE. npe=",I0, " ", I2.2)', this_pe, n
      endif

      do nn = 1, size(Atm)
        call show_atm("1", Atm(nn), nn, this_pe)
      enddo
      print '("[INFO] WDR diag Atm DONE npe=",I0," Atm(",I0,")")', this_pe, n
    endif

    if (is_moving_nest .and. do_move) then
      call mpp_clock_begin (id_movnestTot)
      if (use_timers) call mpp_clock_begin (id_movnest1)

      !!================================================================
      !! Step 1.1 -- Show the nest grids
      !!================================================================

      if (debug_log .and. this_pe .eq. 0) then
        !call show_nest_grid(Atm(n), this_pe, 0)
        print '("[INFO] WDR BD init fv_moving_nest_main.F90 npe=",I0," is=",I0," ie=",I0," js=",I0," je=",I0)', this_pe, Atm(n)%bd%is,  Atm(n)%bd%ie,  Atm(n)%bd%js,  Atm(n)%bd%je
        print '("[INFO] WDR BD init fv_moving_nest_main.F90 npe=",I0," isd=",I0," ied=",I0," jsd=",I0," jed=",I0)', this_pe, Atm(n)%bd%isd,  Atm(n)%bd%ied,  Atm(n)%bd%jsd,  Atm(n)%bd%jed
        print '("[INFO] WDR BD init fv_moving_nest_main.F90 npe=",I0," isc=",I0," iec=",I0," jsc=",I0," jec=",I0)', this_pe, Atm(n)%bd%isc,  Atm(n)%bd%iec,  Atm(n)%bd%jsc, Atm(n)%bd%jec
      endif

      !!================================================================
      !! Step 1.2 -- Configure local variables
      !!================================================================

      x_refine = Atm(child_grid_num)%neststruct%refinement
      y_refine = x_refine
      ioffset => Atm(child_grid_num)%neststruct%ioffset
      joffset => Atm(child_grid_num)%neststruct%joffset

      if (debug_log) print '("[INFO] WDR MV_NST0 fv_moving_nest_main.F90 processing Atm(n) npe=",I0," n=",I0," ioffset=",I0," joffset=",I0)', this_pe, n, ioffset, joffset

      istart_fine = global_nest_domain%istart_fine(nest_num)
      iend_fine = global_nest_domain%iend_fine(nest_num)
      jstart_fine = global_nest_domain%jstart_fine(nest_num)
      jend_fine = global_nest_domain%jend_fine(nest_num)

      istart_coarse = global_nest_domain%istart_coarse(nest_num)
      iend_coarse = global_nest_domain%iend_coarse(nest_num)
      jstart_coarse = global_nest_domain%jstart_coarse(nest_num)
      jend_coarse = global_nest_domain%jend_coarse(nest_num)

      ! Allocate the local weight arrays.  TODO OPTIMIZE change to use the ones from the gridstruct
      if (is_fine_pe) then
        allocate(wt_h(Atm(child_grid_num)%bd%isd:Atm(child_grid_num)%bd%ied, Atm(child_grid_num)%bd%jsd:Atm(child_grid_num)%bd%jed, 4))
        wt_h = real_snan

        allocate(wt_u(Atm(child_grid_num)%bd%isd:Atm(child_grid_num)%bd%ied, Atm(child_grid_num)%bd%jsd:Atm(child_grid_num)%bd%jed+1, 4))
        wt_u = real_snan

        allocate(wt_v(Atm(child_grid_num)%bd%isd:Atm(child_grid_num)%bd%ied+1, Atm(child_grid_num)%bd%jsd:Atm(child_grid_num)%bd%jed, 4))
        wt_v = real_snan
      else
        allocate(wt_h(1,1,4))
        wt_h = 0.0

        allocate(wt_u(1,1,4))
        wt_u = 0.0

        allocate(wt_v(1,1,4))
        wt_v = 0.0
      endif

      ! This full list of PEs is used for the mpp_sync for debugging.  Can later be removed.
      p1 = size(Atm(1)%pelist)   ! Parent PEs
      p2 = size(Atm(2)%pelist)   ! Nest PEs

      allocate(full_pelist(p1 + p2))
      do pp=1,p1
        full_pelist(pp) = Atm(1)%pelist(pp)
      enddo
      do pp=1,p2
        full_pelist(p1+pp) = Atm(2)%pelist(pp)
      enddo

      !!============================================================================
      !! Step 1.3 -- Dump the prognostic variables before we do the nest motion.
      !!============================================================================

      if (debug_log) print '("[INFO] WDR MV_NST0 run step 0 fv_moving_nest_main.F90 npe=",I0)', this_pe
      output_step = output_step + 1

      !!============================================================================
      !! Step 1.4 -- Read in the full panel grid definition
      !!============================================================================

      if (debug_log) then
        print '("[INFO] WDR check grid_global fv_moving_nest_main.F90 npe=",I0," n=",I0)', this_pe, 1
        call check_array(Atm(1)%grid_global, this_pe, "grid_global", -2.0*3.1415926536, 2.0*3.1415926536)
        print '("[INFO] WDR check grid_global fv_moving_nest_main.F90 npe=",I0," n=",I0)', this_pe, 2
        call check_array(Atm(2)%grid_global, this_pe, "grid_global", -2.0*3.1415926536, 2.0*3.1415926536)
      endif

      if (is_fine_pe) then

        nx_cubic = Atm(1)%npx - 1
        ny_cubic = Atm(1)%npy - 1

        nx = Atm(n)%npx - 1
        ny = Atm(n)%npy - 1

        grid => Atm(n)%gridstruct%grid
        agrid => Atm(n)%gridstruct%agrid

        if (debug_log) print '("[INFO] WDR MV_NST0 fv_moving_nest_main.F90 processing Atm(n) npe=",I0," nx_cubic=",I0," ny_cubic=",I0," nx=",I0," ny=",I0)', this_pe, nx_cubic, ny_cubic, nx ,ny

        ! Read in static lat/lon data for parent at nest resolution; returns fp_ full panel variables
        ! Also read in other static variables from the orography and surface files

        if (first_nest_move) then
          if (debug_log) print '("[INFO] WDR mn_latlon_read_hires_parent READING static fine file on npe=",I0)', this_pe

          call mn_latlon_read_hires_parent(Atm(1)%npx, Atm(1)%npy, x_refine, fp_super_tile_geo, &
              Moving_nest(child_grid_num)%mn_flag%surface_dir,  parent_tile)

          call mn_orog_read_hires_parent(Atm(1)%npx, Atm(1)%npy, x_refine, &
              Moving_nest(child_grid_num)%mn_flag%surface_dir, filtered_terrain, &
              mn_static%orog_grid, mn_static%orog_std_grid, mn_static%ls_mask_grid, mn_static%land_frac_grid,  parent_tile)

          ! If terrain_smoother method 1 is chosen, we need the parent coarse terrain
          if (Moving_nest(n)%mn_flag%terrain_smoother .eq. 1) then
            if (filtered_terrain) then
              call mn_static_read_hires(Atm(1)%npx, Atm(1)%npy, 1, Moving_nest(child_grid_num)%mn_flag%surface_dir, "oro_data", "orog_filt", mn_static%parent_orog_grid,  parent_tile)
            else
              call mn_static_read_hires(Atm(1)%npx, Atm(1)%npy, 1, Moving_nest(child_grid_num)%mn_flag%surface_dir, "oro_data", "orog_raw", mn_static%parent_orog_grid,  parent_tile)
            endif
          endif

          call mn_static_read_hires(Atm(1)%npx, Atm(1)%npy, x_refine, trim(Moving_nest(child_grid_num)%mn_flag%surface_dir), "substrate_temperature", "substrate_temperature", mn_static%deep_soil_temp_grid,  parent_tile)
          ! set any -999s to +4C
          call mn_replace_low_values(mn_static%deep_soil_temp_grid, -100.0, 277.0)

          call mn_static_read_hires(Atm(1)%npx, Atm(1)%npy, x_refine, trim(Moving_nest(child_grid_num)%mn_flag%surface_dir), "soil_type", "soil_type", mn_static%soil_type_grid,  parent_tile)
          ! To match initialization behavior, set any -999s to 0 in soil_type
          call mn_replace_low_values(mn_static%soil_type_grid, -100.0, 0.0)


          !! TODO investigate reading high-resolution veg_frac and veg_greenness
          !call mn_static_read_hires(Atm(1)%npx, Atm(1)%npy, x_refine, trim(Moving_nest(child_grid_num)%mn_flag%surface_dir), "", mn_static%veg_frac_grid)

          call mn_static_read_hires(Atm(1)%npx, Atm(1)%npy, x_refine, trim(Moving_nest(child_grid_num)%mn_flag%surface_dir), "vegetation_type", "vegetation_type", mn_static%veg_type_grid,  parent_tile)
          ! To match initialization behavior, set any -999s to 0 in veg_type
          call mn_replace_low_values(mn_static%veg_type_grid, -100.0, 0.0)


          call mn_static_read_hires(Atm(1)%npx, Atm(1)%npy, x_refine, trim(Moving_nest(child_grid_num)%mn_flag%surface_dir), "slope_type", "slope_type", mn_static%slope_type_grid,  parent_tile)
          ! To match initialization behavior, set any -999s to 0 in slope_type
          call mn_replace_low_values(mn_static%slope_type_grid, -100.0, 0.0)


          call mn_static_read_hires(Atm(1)%npx, Atm(1)%npy, x_refine, trim(Moving_nest(child_grid_num)%mn_flag%surface_dir), "maximum_snow_albedo", "maximum_snow_albedo", mn_static%max_snow_alb_grid,  parent_tile)
          ! Set any -999s to 0.5
          call mn_replace_low_values(mn_static%max_snow_alb_grid, -100.0, 0.5)

          ! Albedo fraction -- read and calculate
          call mn_static_read_hires(Atm(1)%npx, Atm(1)%npy, x_refine, trim(Moving_nest(child_grid_num)%mn_flag%surface_dir), "facsf", "facsf", mn_static%facsf_grid,  parent_tile)

          allocate(mn_static%facwf_grid(lbound(mn_static%facsf_grid,1):ubound(mn_static%facsf_grid,1),lbound(mn_static%facsf_grid,2):ubound(mn_static%facsf_grid,2)))

          ! For land points, set facwf = 1.0 - facsf
          ! To match initialization behavior, set any -999s to 0
          do i=lbound(mn_static%facsf_grid,1),ubound(mn_static%facsf_grid,1)
            do j=lbound(mn_static%facsf_grid,2),ubound(mn_static%facsf_grid,2)
              if (mn_static%facsf_grid(i,j) .lt. -100) then
                mn_static%facsf_grid(i,j) = 0
                mn_static%facwf_grid(i,j) = 0
              else
                mn_static%facwf_grid(i,j) = 1.0 - mn_static%facsf_grid(i,j)
              endif
            enddo
          enddo

          ! Additional albedo variables
          !  black sky = strong cosz -- direct sunlight
          !  white sky = weak cosz -- diffuse light

          ! alvsf = visible strong cosz = visible_black_sky_albedo
          ! alvwf = visible weak cosz = visible_white_sky_albedo
          ! alnsf = near IR strong cosz = near_IR_black_sky_albedo
          ! alnwf = near IR weak cosz = near_IR_white_sky_albedo

          call mn_static_read_hires(Atm(1)%npx, Atm(1)%npy, x_refine, trim(Moving_nest(child_grid_num)%mn_flag%surface_dir), "snowfree_albedo", "visible_black_sky_albedo", mn_static%alvsf_grid,  parent_tile, time=month)
          call mn_static_read_hires(Atm(1)%npx, Atm(1)%npy, x_refine, trim(Moving_nest(child_grid_num)%mn_flag%surface_dir), "snowfree_albedo", "visible_white_sky_albedo", mn_static%alvwf_grid,  parent_tile, time=month)

          call mn_static_read_hires(Atm(1)%npx, Atm(1)%npy, x_refine, trim(Moving_nest(child_grid_num)%mn_flag%surface_dir), "snowfree_albedo", "near_IR_black_sky_albedo", mn_static%alnsf_grid,  parent_tile, time=month)
          call mn_static_read_hires(Atm(1)%npx, Atm(1)%npy, x_refine, trim(Moving_nest(child_grid_num)%mn_flag%surface_dir), "snowfree_albedo", "near_IR_white_sky_albedo", mn_static%alnwf_grid,  parent_tile, time=month)

          ! Set the -999s to small value of 0.06, matching initialization code in chgres

          call mn_replace_low_values(mn_static%alvsf_grid, -100.0, 0.06)
          call mn_replace_low_values(mn_static%alvwf_grid, -100.0, 0.06)
          call mn_replace_low_values(mn_static%alnsf_grid, -100.0, 0.06)
          call mn_replace_low_values(mn_static%alnwf_grid, -100.0, 0.06)

        endif

        !  Validation/logging calls that can be disabled
        if (debug_log) then
          call show_tile_geo(fp_super_tile_geo, this_pe, "fp_super_tile_geo")
          call show_gridstruct(Atm(n)%gridstruct, this_pe)
          !call validate_hires_parent(fp_super_tile_geo, Atm(n)%gridstruct%grid, Atm(n)%gridstruct%agrid, x_refine, y_refine, ioffset, joffset)
        endif
      endif

      if (first_nest_move) first_nest_move = .false.

      if (use_timers) call mpp_clock_end (id_movnest1)
      if (use_timers) call mpp_clock_begin (id_movnest1_9)

      !!=====================================================================================
      !! Step 1.9 -- Allocate and fill the temporary variable(s)
      !!=====================================================================================

      call mn_prog_fill_temp_variables(Atm, n, child_grid_num, is_fine_pe, npz)
      call mn_phys_fill_temp_variables(Atm, Atm_block, IPD_control, IPD_data, n, child_grid_num, is_fine_pe, npz)

      if (use_timers) call mpp_clock_end (id_movnest1_9)
      if (use_timers) call mpp_clock_begin (id_movnest2)

      !!============================================================================
      !! Step 2 -- Fill in the halos from the coarse grids
      !!============================================================================
      if (debug_log) print '("WDR_NEST_HALO_RECV,",I0,"===STEP 2====")', this_pe
      if (debug_log) print '("[INFO] WDR MV_NST2 run step 2 fv_moving_nest_main.F90 npe=",I0)', this_pe

      !  The halos seem to be empty at least on the first model timestep.
      !  These calls need to be executed by the parent and nest PEs in order to do the communication
      !  This is before any nest motion has occurred

      call mn_prog_fill_nest_halos_from_parent(Atm, n, child_grid_num, is_fine_pe, global_nest_domain, nz)
      call mn_phys_fill_nest_halos_from_parent(Atm, IPD_control, IPD_data, mn_static, n, child_grid_num, is_fine_pe, global_nest_domain, nz)

      if (use_timers) call mpp_clock_end (id_movnest2)
      if (use_timers) call mpp_clock_begin (id_movnest3)

      !!============================================================================
      !! Step 3 -- Redefine the nest domain to new location
      !!   This calls mpp_define_nest_domains.  Following the code in fv_control.F90, only should
      !!   be executed on the nest PEs. Operates only on indices.
      !!  --  Similar to med_nest_configure() from HWRF
      !!============================================================================

      if (debug_log) print '("WDR_NEST_HALO_RECV,",I0,"===STEP 3====")', this_pe
      if (debug_log) print '("[INFO] WDR MV_NST3 run step 3 fv_moving_nest_main.F90 npe=",I0)', this_pe

      if (debug_log) print '("[INFO] WDR MV_NST3 run step 3 fv_moving_nest_main.F90 processing Atm(n) npe=",I0," n=",I0," ioffset=",I0," joffset=",I0)', this_pe, n, ioffset, joffset

      call mn_meta_move_nest(delta_i_c, delta_j_c, pelist, is_fine_pe, extra_halo, &
          global_nest_domain, domain_fine, domain_coarse, &
          istart_coarse, iend_coarse, jstart_coarse, jend_coarse,  &
          istart_fine, iend_fine, jstart_fine, jend_fine)

      ! This code updates the values in neststruct; ioffset/joffset are pointers:  ioffset => Atm(child_grid_num)%neststruct%ioffset
      ioffset = ioffset + delta_i_c
      joffset = joffset + delta_j_c

      if (debug_sync) call mpp_sync(full_pelist)   ! Used to make debugging easier.  Can be removed.

      if (use_timers) call mpp_clock_end (id_movnest3)
      if (use_timers) call mpp_clock_begin (id_movnest4)

      !!============================================================================
      !! Step 4  -- Fill the internal nest halos for the prognostic variables,
      !!           then physics variables
      !!    Only acts on the nest PEs
      !!    --  similar to med_nest_initial
      !!============================================================================

      if (debug_log) print '("WDR_NEST_HALO_RECV,",I0,"===STEP 4====")', this_pe
      if (debug_log) print '("[INFO] WDR MV_NST4 run step 4 fv_moving_nest_main.F90 npe=",I0)', this_pe

      ! TODO should/can this run before the mn_meta_move_nest?
      if (is_fine_pe) then
        call mn_prog_fill_intern_nest_halos(Atm(n), domain_fine, is_fine_pe)
        call mn_phys_fill_intern_nest_halos(Moving_nest(n), IPD_control, IPD_data, domain_fine, is_fine_pe)
      else
        if (debug_log) print '("[INFO] WDR MV_NST4 skip step 4 fv_moving_nest_main.F90 npe=",I0)', this_pe
      endif

      if (debug_sync) call mpp_sync(full_pelist)   ! Used to make debugging easier.  Can be removed.

      if (use_timers) call mpp_clock_end (id_movnest4)
      if (use_timers) call mpp_clock_begin (id_movnest5)

      !!============================================================================
      !! Step 5  --  Recalculate nest halo weights (for fine PEs only) and indices
      !!   -- Similiar to med_nest_weights
      !!============================================================================

      if (debug_log) print '("WDR_NEST_HALO_RECV,",I0,"===STEP 5====")', this_pe
      if (debug_log) print '("[INFO] WDR MV_NST5 run step 5 fv_moving_nest_main.F90 npe=",I0)', this_pe

      if (is_fine_pe) then
        !!============================================================================
        !! Step 5.1 -- Fill the p_grid* and n_grid* variables
        !!============================================================================
        if (debug_log) print '("[INFO] WDR MV_NST5 run step 5 fv_moving_nest_main.F90 npe=",I0, " tile_geo%lats allocated:",L1)', this_pe, allocated(tile_geo%lats)
        if (debug_log) print '("[INFO] WDR MV_NST5 run step 5 fv_moving_nest_main.F90 npe=",I0, " parent_geo%lats allocated:",L1)', this_pe, allocated(parent_geo%lats)

        ! parent_geo is only loaded first time; afterwards it is reused.
        ! This is the coarse resolution data for the parent
        call mn_latlon_load_parent(Moving_nest(child_grid_num)%mn_flag%surface_dir, Atm, n, parent_tile, &
            delta_i_c, delta_j_c, child_grid_num, &
            parent_geo, tile_geo, tile_geo_u, tile_geo_v, fp_super_tile_geo, &
            p_grid, n_grid, p_grid_u, n_grid_u, p_grid_v, n_grid_v)

        ! tile_geo holds the center lat/lons for the entire nest (all PEs).
        call mn_reset_phys_latlon(Atm, n, tile_geo, fp_super_tile_geo, Atm_block, IPD_control, IPD_data)

        !!============================================================================
        !! Step 5.2 -- Fill the wt* variables for each stagger
        !!============================================================================

        call mn_meta_recalc( delta_i_c, delta_j_c, x_refine, y_refine, tile_geo, parent_geo, fp_super_tile_geo, &
            is_fine_pe, global_nest_domain, position, p_grid, n_grid, wt_h, istart_coarse, jstart_coarse)

        call mn_meta_recalc( delta_i_c, delta_j_c, x_refine, y_refine, tile_geo_u, parent_geo, fp_super_tile_geo, &
            is_fine_pe, global_nest_domain, position_u, p_grid_u, n_grid_u, wt_u, istart_coarse, jstart_coarse)

        call mn_meta_recalc( delta_i_c, delta_j_c, x_refine, y_refine, tile_geo_v, parent_geo, fp_super_tile_geo, &
            is_fine_pe, global_nest_domain, position_v, p_grid_v, n_grid_v, wt_v, istart_coarse, jstart_coarse)

      endif

      !!============================================================================
      !! Step 5.3 -- Adjust the indices by the values of delta_i_c, delta_j_c
      !!============================================================================

      call mn_shift_index(delta_i_c, delta_j_c, Atm(child_grid_num)%neststruct%ind_h)
      call mn_shift_index(delta_i_c, delta_j_c, Atm(child_grid_num)%neststruct%ind_u)
      call mn_shift_index(delta_i_c, delta_j_c, Atm(child_grid_num)%neststruct%ind_v)
      call mn_shift_index(delta_i_c, delta_j_c, Atm(child_grid_num)%neststruct%ind_b)

      if (debug_sync) call mpp_sync(full_pelist)   ! Used to make debugging easier.  Can be removed.

      if (use_timers) call mpp_clock_end (id_movnest5)
      if (use_timers) call mpp_clock_begin (id_movnest6)

      !!============================================================================
      !! Step 6   Shift the data on each nest PE
      !!            -- similar to med_nest_move in HWRF
      !!============================================================================

      if (debug_log) print '("WDR_NEST_HALO_RECV,",I0,"===STEP 6====")', this_pe
      if (debug_log) print '("[INFO] WDR MV_NST6 run step 6 fv_moving_nest_main.F90 npe=",I0," n=",I0)', this_pe, n

      call mn_prog_shift_data(Atm, n, child_grid_num, wt_h, wt_u, wt_v, &
          delta_i_c, delta_j_c, x_refine, y_refine, &
          is_fine_pe, global_nest_domain, nz)

      call mn_phys_shift_data(Atm, IPD_control, IPD_data, n, child_grid_num, wt_h, wt_u, wt_v, &
          delta_i_c, delta_j_c, x_refine, y_refine, &
          is_fine_pe, global_nest_domain, nz)

      if (debug_log) print '("[INFO] WDR MV_NST6 complete step 6 fv_moving_nest_main.F90 npe=",I0)', this_pe

      if (debug_sync) call mpp_sync(full_pelist)   ! Used to make debugging easier.  Can be removed.

      if (use_timers) call mpp_clock_end (id_movnest6)
      if (use_timers) call mpp_clock_begin (id_movnest7_0)

      !!=====================================================================================
      !! Step 7 --  Reset the grid definition data and buffer sizes and weights after the nest motion
      !!             Mostly needed when dynamics is executed
      !!=====================================================================================

      if (debug_log) print '("WDR_NEST_HALO_RECV,",I0,"===STEP 7====")', this_pe
      if (debug_log) print '("[INFO] WDR MV_NST7 run step 7 fv_moving_nest_main.F90 npe=",I0," n=",I0)', this_pe, n

      call mn_meta_reset_gridstruct(Atm, n, child_grid_num, global_nest_domain, fp_super_tile_geo, x_refine, y_refine, is_fine_pe, wt_h, wt_u, wt_v, a_step, dt_atmos)

      if (use_timers) call mpp_clock_end (id_movnest7_0)
      if (use_timers) call mpp_clock_begin (id_movnest7_1)

      !!=====================================================================================
      !! Step 7.01 --  Reset the orography data that was read from the hires static file
      !!
      !!=====================================================================================

      if (is_fine_pe) then
        ! phis is allocated in fv_arrays.F90 as:  allocate ( Atm%phis(isd:ied  ,jsd:jed  ) )
        ! 0 -- all high-resolution data, 1 - static nest smoothing algorithm, 5 - 5 point smoother, 9 - 9 point smoother
        ! Defaults to 1 - static nest smoothing algorithm; this seems to produce the most stable solutions
        !print '("[INFO] WDR Moving Nest terrain_smoother=",I0," High-resolution terrain. npe=",I0)', Atm(n)%neststruct%terrain_smoother, this_pe

        select case(Moving_nest(n)%mn_flag%terrain_smoother)
        case (0)
          ! High-resolution terrain for entire nest
          if (debug_log) print '("[INFO] WDR Moving Nest terrain_smoother=0 High-resolution terrain. npe=",I0)', this_pe
          Atm(n)%phis(isd:ied, jsd:jed) = mn_static%orog_grid((ioffset-1)*x_refine+isd:(ioffset-1)*x_refine+ied, (joffset-1)*y_refine+jsd:(joffset-1)*y_refine+jed) * grav
        case (1)
          ! Static nest smoothing algorithm - interpolation of coarse terrain in halo zone and 5 point blending zone of coarse and fine data
          if (debug_log) print '("[INFO] WDR Moving Nest terrain_smoother=1 Blending5 algorithm. npe=",I0)', this_pe
          call set_blended_terrain(Atm(n), mn_static%parent_orog_grid, mn_static%orog_grid, x_refine, Atm(n)%bd%ng, 5, a_step)
        case (2)
          ! Static nest smoothing algorithm - interpolation of coarse terrain in halo zone and 5 point blending zone of coarse and fine data
          if (debug_log) print '("[INFO] WDR Moving Nest terrain_smoother=1 Blending10 algorithm. npe=",I0)', this_pe
          call set_blended_terrain(Atm(n), mn_static%parent_orog_grid, mn_static%orog_grid, x_refine, Atm(n)%bd%ng, 10, a_step)
        case (5)
          ! 5 pt smoother.  blend zone of 5 to match static nest
          if (debug_log) print '("[INFO] WDR Moving Nest terrain_smoother=5  5-point smoother. npe=",I0)', this_pe
          call set_smooth_nest_terrain(Atm(n), mn_static%orog_grid, x_refine, 5, Atm(n)%bd%ng, 5)
        case (9)
          ! 9 pt smoother.  blend zone of 5 to match static nest
          if (debug_log) print '("[INFO] WDR Moving Nest terrain_smoother=9  9-point smoother. npe=",I0)', this_pe
          call set_smooth_nest_terrain(Atm(n), mn_static%orog_grid, x_refine, 9, Atm(n)%bd%ng, 5)
        case default
          write (errstring, "(I0)") Moving_nest(n)%mn_flag%terrain_smoother
          call mpp_error(FATAL,'Invalid terrain_smoother in fv_moving_nest_main '//errstring)
        end select

        ! Reinitialize diagnostics -- zsurf which is g * Atm%phis
        call fv_diag_reinit(Atm(n:n))

        ! sgh and oro were only fully allocated if fv_land is True
        !      if false, oro is (1,1), and sgh is not allocated
        if ( Atm(n)%flagstruct%fv_land ) then
          if (debug_log) print '("[INFO] WDR shift orography data fv_land TRUE npe=",I0)', this_pe
          ! oro and sgh are allocated only for the compute domain -- they do not have halos

          !fv_arrays.F90 oro() !< land fraction (1: all land; 0: all water)
          !real, _ALLOCATABLE :: oro(:,:)      _NULL  !< land fraction (1: all land; 0: all water)
          Atm(n)%oro(isc:iec, jsc:jec) = mn_static%land_frac_grid((ioffset-1)*x_refine+isc:(ioffset-1)*x_refine+iec, (joffset-1)*y_refine+jsc:(joffset-1)*y_refine+jec)

          !real, _ALLOCATABLE :: sgh(:,:)      _NULL  !< Terrain standard deviation
          Atm(n)%sgh(isc:iec, jsc:jec) = mn_static%orog_std_grid((ioffset-1)*x_refine+isc:(ioffset-1)*x_refine+iec, (joffset-1)*y_refine+jsc:(joffset-1)*y_refine+jec)
        else
          if (debug_log) print '("[INFO] WDR shift orography data fv_land FALSE npe=",I0)', this_pe
        endif

        call mn_phys_reset_sfc_props(Atm, n, mn_static, Atm_block, IPD_data, ioffset, joffset, x_refine)
      endif

      !!=====================================================================================
      !! Step 7.1   Refill the nest edge halos from parent grid after nest motion
      !!            Parent and nest PEs need to execute these subroutines
      !!=====================================================================================

      ! Refill the halos around the edge of the nest from the parent
      call mn_prog_fill_nest_halos_from_parent(Atm, n, child_grid_num, is_fine_pe, global_nest_domain, nz)
      call mn_phys_fill_nest_halos_from_parent(Atm, IPD_control, IPD_data, mn_static, n, child_grid_num, is_fine_pe, global_nest_domain, nz)

      if (use_timers) call mpp_clock_end (id_movnest7_1)

      if (is_fine_pe) then
        if (use_timers) call mpp_clock_begin (id_movnest7_2)

        ! Refill the internal halos after nest motion
        call mn_prog_fill_intern_nest_halos(Atm(n), domain_fine, is_fine_pe)
        call mn_phys_fill_intern_nest_halos(Moving_nest(n), IPD_control, IPD_data, domain_fine, is_fine_pe)

        if (use_timers) call mpp_clock_end (id_movnest7_2)
      endif

      if (debug_sync) call mpp_sync(full_pelist)   ! Used to make debugging easier.  Can be removed.

      !!=====================================================================================
      !! Step 7.3 -- Apply the temporary variable to the prognostics and physics structures
      !!=====================================================================================
      if (use_timers) call mpp_clock_begin (id_movnest7_3)

      call mn_prog_apply_temp_variables(Atm, n, child_grid_num, is_fine_pe, npz)
      call mn_phys_apply_temp_variables(Atm, Atm_block, IPD_control, IPD_data, n, child_grid_num, is_fine_pe, npz)

      if (use_timers) call mpp_clock_end (id_movnest7_3)
      if (use_timers) call mpp_clock_begin (id_movnest8)

      !!============================================================================
      !!  Step 8 -- Dump to netCDF
      !!============================================================================

      if (debug_log) print '("WDR_NEST_HALO_RECV,",I0,"===STEP 8====")', this_pe
      if (debug_log) print '("[INFO] WDR MV_NST8 run step 8 fv_moving_nest_main.F90 npe=",I0)', this_pe

      if (is_fine_pe) then
        do i=isc,iec
          do j=jsc,jec
            ! WDR EMIS PATCH - Force to positive at all locations matching the landmask
            !if (Moving_nest(n)%mn_phys%slmsk(i,j) .eq. 1 .and. Moving_nest(n)%mn_phys%emis_lnd(i,j) .lt. 0.0) Moving_nest(n)%mn_phys%emis_lnd(i,j) = 0.5
            !if (Moving_nest(n)%mn_phys%slmsk(i,j) .eq. 2 .and. Moving_nest(n)%mn_phys%emis_ice(i,j) .lt. 0.0) Moving_nest(n)%mn_phys%emis_ice(i,j) = 0.5
            !if (Moving_nest(n)%mn_phys%slmsk(i,j) .eq. 0 .and. Moving_nest(n)%mn_phys%emis_wat(i,j) .lt. 0.0) Moving_nest(n)%mn_phys%emis_wat(i,j) = 0.5
            !if (Moving_nest(n)%mn_phys%slmsk(i,j) .eq. 1 .and. Moving_nest(n)%mn_phys%albdirvis_lnd(i,j) .lt. 0.0) Moving_nest(n)%mn_phys%albdirvis_lnd(i,j) = 0.5
            !if (Moving_nest(n)%mn_phys%slmsk(i,j) .eq. 1 .and. Moving_nest(n)%mn_phys%albdirnir_lnd(i,j) .lt. 0.0) Moving_nest(n)%mn_phys%albdirvis_lnd(i,j) = 0.5
            !if (Moving_nest(n)%mn_phys%slmsk(i,j) .eq. 1 .and. Moving_nest(n)%mn_phys%albdifvis_lnd(i,j) .lt. 0.0) Moving_nest(n)%mn_phys%albdifvis_lnd(i,j) = 0.5
            !if (Moving_nest(n)%mn_phys%slmsk(i,j) .eq. 1 .and. Moving_nest(n)%mn_phys%albdifnir_lnd(i,j) .lt. 0.0) Moving_nest(n)%mn_phys%albdifnir_lnd(i,j) = 0.5

            ! WDR EMIS PATCH - Force to positive at all locations.
            if (Moving_nest(n)%mn_phys%emis_lnd(i,j) .lt. 0.0) Moving_nest(n)%mn_phys%emis_lnd(i,j) = 0.5
            if (Moving_nest(n)%mn_phys%emis_ice(i,j) .lt. 0.0) Moving_nest(n)%mn_phys%emis_ice(i,j) = 0.5
            if (Moving_nest(n)%mn_phys%emis_wat(i,j) .lt. 0.0) Moving_nest(n)%mn_phys%emis_wat(i,j) = 0.5
            if (Moving_nest(n)%mn_phys%albdirvis_lnd(i,j) .lt. 0.0) Moving_nest(n)%mn_phys%albdirvis_lnd(i,j) = 0.5
            if (Moving_nest(n)%mn_phys%albdirnir_lnd(i,j) .lt. 0.0) Moving_nest(n)%mn_phys%albdirvis_lnd(i,j) = 0.5
            if (Moving_nest(n)%mn_phys%albdifvis_lnd(i,j) .lt. 0.0) Moving_nest(n)%mn_phys%albdifvis_lnd(i,j) = 0.5
            if (Moving_nest(n)%mn_phys%albdifnir_lnd(i,j) .lt. 0.0) Moving_nest(n)%mn_phys%albdifnir_lnd(i,j) = 0.5

            !if (Moving_nest(n)%mn_phys%semis(i,j) .lt. 0.0) then
            !   print '("[INFO] WDR SEMIS fv_moving_nest_main.F90 npe=",I0," semis(",I0,",",I0,")=",F15.5)', this_pe, i, j, Moving_nest(n)%mn_phys%semis(i,j)
            !endif
            !if (Moving_nest(n)%mn_phys%semisbase(i,j) .lt. 0.0) then
            !   print '("[INFO] WDR SEMISBASE fv_moving_nest_main.F90 npe=",I0," semisbase(",I0,",",I0,")=",F15.5)', this_pe, i, j, Moving_nest(n)%mn_phys%semisbase(i,j)
            !endif

            if ( Moving_nest(n)%mn_phys%slmsk(i,j) .eq. 1 .and.  Moving_nest(n)%mn_phys%emis_lnd(i,j) .lt. 0.0) then
              print '("[INFO] WDR SEMISLND fv_moving_nest_main.F90 npe=",I0," emis_lnd(",I0,",",I0,")=",F15.5)', this_pe, i, j, Moving_nest(n)%mn_phys%emis_lnd(i,j)
            endif
            if ( Moving_nest(n)%mn_phys%slmsk(i,j) .eq. 2 .and. Moving_nest(n)%mn_phys%emis_ice(i,j) .lt. 0.0) then
              print '("[INFO] WDR SEMISLND fv_moving_nest_main.F90 npe=",I0," emis_ice(",I0,",",I0,")=",F15.5)', this_pe, i, j, Moving_nest(n)%mn_phys%emis_ice(i,j)
            endif
            if ( Moving_nest(n)%mn_phys%slmsk(i,j) .eq. 0 .and. Moving_nest(n)%mn_phys%emis_wat(i,j) .lt. 0.0) then
              print '("[INFO] WDR SEMISLND fv_moving_nest_main.F90 npe=",I0," emis_wat(",I0,",",I0,")=",F15.5)', this_pe, i, j, Moving_nest(n)%mn_phys%emis_wat(i,j)
            endif
            if ( Moving_nest(n)%mn_phys%slmsk(i,j) .eq. 1 .and. Moving_nest(n)%mn_phys%albdirvis_lnd(i,j) .lt. 0.0) then
              print '("[INFO] WDR ALBLND fv_moving_nest_main.F90 npe=",I0," albdirvis_lnd(",I0,",",I0,")=",F15.5)', this_pe, i, j, Moving_nest(n)%mn_phys%albdirvis_lnd(i,j)
            endif
          enddo
        enddo
      endif

      output_step = output_step + 1

      if (debug_sync) call mpp_sync(full_pelist)   ! Used to make debugging easier.  Can be removed.

      if (use_timers) call mpp_clock_end (id_movnest8)
      if (use_timers) call mpp_clock_begin (id_movnest9)

      !!=========================================================================================
      !! Step 9 -- Recalculate auxiliary pressures
      !!           Should help stabilize the fields before dynamics runs
      !! TODO Consider whether vertical remapping, recalculation of omega, interpolation of winds
      !!  to A or C grids, and/or divergence recalculation are needed here.
      !!=========================================================================================

      if (is_fine_pe) then
        if (debug_log) print '("[INFO] WDR MV_NST L2E before recalc auxiliary pressures fv_moving_nest_main.F90 npe=",I0)', this_pe
        call recalc_aux_pressures(Atm(n))
        if (debug_log) print '("[INFO] WDR MV_NST L2E after recalc auxiliary pressures fv_moving_nest_main.F90 npe=",I0)', this_pe
      endif

      if (debug_log) print '("[INFO] WDR PTVAL fv_dynamics.F90 npe=",I0," AfterNestMove ================================================")', this_pe
      output_step = output_step + 1
    endif

    if (use_timers) call mpp_clock_end (id_movnest9)
    call mpp_clock_end (id_movnestTot)

    if (debug_sync) call mpp_sync(full_pelist)   ! Used to make debugging easier.  Can be removed.

    !! Uncomment to exit and force file IO after single nest move, without dynamics
    !    call fms_io_exit()   !! Force the output of the buffered NC files
    !    if (debug_log) print '("[INFO] WDR calling mpp_exit after moving nest fv_moving_nest_main.F90 npe=",I0)', this_pe
    !    call mpp_exit()
    !    if (debug_log) print '("[INFO] WDR calling STOP after moving nest fv_moving_nest_main.F90 npe=",I0)', this_pe
    !    stop
    !!  else
    !!    if (debug_log) print '("[INFO] WDR move_nest not nested PE  npe=",I0)', this_pe
    !!  endif

    !call compare_terrain("phis", Atm(n)%phis, 1, Atm(n)%neststruct%ind_h, x_refine, y_refine, is_fine_pe, global_nest_domain)

    if (debug_log) call show_nest_grid(Atm(n), this_pe, 99)

  end subroutine fv_moving_nest_exec

  !>@brief The subroutine 'mn_replace_low_values' replaces low values with a default value.
  subroutine mn_replace_low_values(data_grid, low_value, new_value)
    real, _ALLOCATABLE, intent(inout)   :: data_grid(:,:)  !< 2D grid of data
    real, intent(in)                    :: low_value       !< Low value to check for; e.g. negative or fill value
    real, intent(in)                    :: new_value       !< Value to replace low value with

    integer :: i, j

    do i=lbound(data_grid,1),ubound(data_grid,1)
      do j=lbound(data_grid,2),ubound(data_grid,2)
        if (data_grid(i,j) .le. low_value) data_grid(i,j) = new_value
      enddo
    enddo
  end subroutine mn_replace_low_values

#endif ! MOVING_NEST

end module fv_moving_nest_main_mod

