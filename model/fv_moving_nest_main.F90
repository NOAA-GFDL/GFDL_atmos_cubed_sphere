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

!>@brief The module 'fv_moving_nest_main' provides the top-level 
!! interface for moving nests in the Cubed-Sphere FV dynamical core

module fv_moving_nest_main_mod
#ifdef MOVING_NEST
!-----------------------------------------
! Moving Nest Top Level Functionality
! W. Ramstrom - AOML/HRD/CIMAS 05/27/2021
!-----------------------------------------


#include <fms_platform.h>

!-----------------
! FMS modules:
!-----------------
use block_control_mod,      only: block_control_type
use constants_mod,          only: cp_air, rdgas, grav, rvgas, kappa, pstd_mks
use time_manager_mod,       only: time_type, get_time, set_time, operator(+), &
                                  operator(-), operator(/), time_type_to_real
use fms_mod,                only: file_exist, open_namelist_file,    &
                                  close_file, error_mesg, FATAL,     &
                                  check_nml_error, stdlog,           &
                                  write_version_number,              &
                                  set_domain,   &
                                  mpp_clock_id, mpp_clock_begin,     &
                                  mpp_clock_end, CLOCK_SUBCOMPONENT, &
                                  clock_flag_default, nullify_domain
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
use atmosphere_mod,     only: Atm, mygrid, p_split
use fv_arrays_mod,      only: fv_atmos_type, R_GRID, fv_grid_bounds_type, phys_diag_type
use fv_control_mod,     only: fv_control_init, fv_end, ngrids
use fv_eta_mod,         only: get_eta_level
use fv_fill_mod,        only: fill_gfs
use fv_dynamics_mod,    only: fv_dynamics
use fv_nesting_mod,     only: twoway_nesting
use fv_diagnostics_mod, only: fv_diag_init, fv_diag_reinit, fv_diag, fv_time, prt_maxmin, prt_height
use fv_nggps_diags_mod, only: fv_nggps_diag_init, fv_nggps_diag, fv_nggps_tavg
use fv_restart_mod,     only: fv_restart, fv_write_restart
use fv_timing_mod,      only: timing_on, timing_off
use fv_mp_mod,          only: is_master
use fv_sg_mod,          only: fv_subgrid_z
use fv_update_phys_mod, only: fv_update_phys
use fv_io_mod,          only: fv_io_register_nudge_restart
use fv_nwp_nudge_mod,   only: fv_nwp_nudge_init, fv_nwp_nudge_end, do_adiabatic_init
use fv_regional_mod,    only: start_regional_restart, read_new_bc_data, &
                              a_step, p_step, current_time_in_seconds
use fv_grid_utils_mod,  only: g_sum
use mpp_domains_mod, only:  mpp_get_data_domain, mpp_get_compute_domain
use coarse_graining_mod, only: coarse_graining_init
use coarse_grained_diagnostics_mod, only: fv_coarse_diag_init, fv_coarse_diag
use coarse_grained_restart_files_mod, only: fv_coarse_restart_init
use diag_manager_mod,   only: send_data


!-----------------------------------------
!  External routines
!-----------------------------------------
use fms_io_mod,         only: fms_io_exit
use mpp_domains_mod,    only: NORTH, NORTH_EAST, EAST, SOUTH_EAST, CORNER, CENTER
use mpp_domains_mod,    only: nest_domain_type
use mpp_mod,            only: mpp_sync, mpp_exit
use mpp_domains_mod,    only: mpp_get_global_domain
use mpp_mod,            only: mpp_send, mpp_sync_self

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
use fv_moving_nest_mod,         only: mn_phys_fill_intern_nest_halos, mn_phys_fill_nest_halos_from_parent, &
     mn_phys_dump_to_netcdf, mn_phys_shift_data

!      Metadata routines
use fv_moving_nest_mod,         only: mn_meta_move_nest, mn_meta_recalc, mn_meta_reset_gridstruct, mn_shift_index

!      Temporary variable routines (delz)
use fv_moving_nest_mod,         only: mn_prog_fill_temp_variables, mn_prog_apply_temp_variables
use fv_moving_nest_mod,         only: mn_phys_fill_temp_variables, mn_phys_apply_temp_variables

!      Load static datasets
use fv_moving_nest_mod,         only: mn_latlon_read_hires_parent, mn_latlon_load_parent, mn_reset_phys_latlon
use fv_moving_nest_mod,         only: mn_orog_read_hires_parent, mn_static_read_hires
use fv_moving_nest_utils_mod,   only: load_nest_latlons_from_nc, compare_terrain

!      Bounds checking routines
use fv_moving_nest_mod,         only: permit_move_nest

!      Grid reset routines
use fv_moving_nest_mod,         only: grid_geometry, assign_n_p_grids, move_nest_geo
use fv_moving_nest_utils_mod,   only: fill_grid_from_supergrid, fill_weight_grid

!      Recalculation routines
use fv_moving_nest_mod,         only: reallocate_BC_buffers, recalc_aux_pressures, vertical_remap_nest !, reinit_parent_indices

!      Logging and debugging information
use fv_moving_nest_mod,         only: check_array
use fv_moving_nest_logging_mod, only: show_atm, show_atm_grids, show_tile_geo, show_nest_grid, show_gridstruct, grid_equal
use fv_moving_nest_logging_mod, only: validate_hires_parent



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

logical :: debug_log = .false.
logical :: tsvar_out = .true.
logical :: wxvar_out = .false.


!  --- Clock ids for moving_nest performance metering
integer :: id_movnest1, id_movnest1_9, id_movnest2, id_movnest3, id_movnest4, id_movnest5
integer :: id_movnest6, id_movnest7_0, id_movnest7_1, id_movnest7_2, id_movnest7_3, id_movnest8, id_movnest9
integer :: id_movnestTot
!integer, save :: a_step = 0
integer, save :: output_step = 0

type mn_surface_grids
   real, allocatable  :: orog_grid(:,:)                _NULL  ! orography -- raw or filtered depending on namelist option, in meters
   real, allocatable  :: orog_std_grid(:,:)            _NULL  ! terrain standard deviation for gravity wave drag, in meters (?)
   real, allocatable  :: ls_mask_grid(:,:)             _NULL  ! land sea mask -- 0 for ocean/lakes, 1, for land.  Perhaps 2 for sea ice. 
   real, allocatable  :: land_frac_grid(:,:)           _NULL  ! Continuous land fraction - 0.0 ocean, 0.5 half of each, 1.0 all land

   ! Soil variables
   real, allocatable  :: deep_soil_temp_grid(:,:)      _NULL  ! deep soil temperature at 5m, in degrees K
   real, allocatable  :: soil_type_grid(:,:)           _NULL  ! STATSGO soil type

   ! Vegetation variables
   real, allocatable  :: veg_frac_grid(:,:)           _NULL  ! vegetation fraction 
   real, allocatable  :: veg_type_grid(:,:)           _NULL  ! IGBP vegetation type
   real, allocatable  :: veg_greenness_grid(:,:)      _NULL  ! NESDIS vegetation greenness; netCDF file has monthly values

   ! Orography variables
   real, allocatable  :: slope_type_grid(:,:)         _NULL  ! legacy 1 degree GFS slope type 

   ! Albedo variables
   real, allocatable  :: facsf_grid(:,:)              _NULL  ! legacy 1 degree GFS fractional coverage for strong/weak zenith angle dependent albedo
   real, allocatable  :: max_snow_alb_grid(:,:)       _NULL  ! max snow albedo
   ! Snow free albedo
   real, allocatable  :: vis_black_alb_grid(:,:)      _NULL  ! Visible black sky albeo; netCDF file has monthly values
   real, allocatable  :: vis_white_alb_grid(:,:)      _NULL  ! Visible white sky albeo; netCDF file has monthly values
   real, allocatable  :: ir_black_alb_grid(:,:)       _NULL  ! Near IR black sky albeo; netCDF file has monthly values
   real, allocatable  :: ir_white_alb_grid(:,:)       _NULL  ! Near IR white sky albeo; netCDF file has monthly values

end type mn_surface_grids



contains

  subroutine update_moving_nest(Atm_block, IPD_control, IPD_data, time_step)
    type(block_control_type), intent(in) :: Atm_block
    type(IPD_control_type), intent(in) :: IPD_control
    type(IPD_data_type), intent(inout) :: IPD_data(:)
    type(time_type), intent(in)     :: time_step

    logical :: is_moving_nest = .True.  !! TODO connect to namelist for each nest

    real :: dt_atmos = 90.0  !! TODO connect to timestep

    logical :: do_move
    integer :: delta_i_c, delta_j_c
    integer :: parent_grid_num, child_grid_num, nest_num
    integer :: n
    
    do_move = .false.
    
    ! Get dt_atmos from Atm()%Time_step_atmos -  seems like some transformations might be needed
    
    n = mygrid   ! Public variable from atmosphere.F90
    ! These will need to be looked up on each PE when multiple and telescoped nests are enabled.
    parent_grid_num = 1 
    child_grid_num = 2
    nest_num = 1
    
    if (is_moving_nest) then
       call eval_move_nest(Atm, a_step, do_move, delta_i_c, delta_j_c, dt_atmos)
       if (do_move) then
          ! Verifies if nest motion is permitted
          ! If nest would cross the cube face edge, do_move is reset to .false. and delta_i_c, delta_j_c are set to 0
          call permit_move_nest(Atm, a_step, parent_grid_num, child_grid_num, delta_i_c, delta_j_c, do_move)
          
       end if
       
       if (do_move) then
          call fv_moving_nest_exec(Atm, Atm_block, IPD_control, IPD_data, delta_i_c, delta_j_c, n, nest_num, parent_grid_num, child_grid_num, dt_atmos)
       end if
       
    end if
    
    !a_step = a_step + 1
    
  end subroutine update_moving_nest


  subroutine dump_moving_nest(Atm_block, IPD_control, IPD_data, time_step)
    type(block_control_type), intent(in) :: Atm_block
    type(IPD_control_type), intent(in) :: IPD_control
    type(IPD_data_type), intent(inout) :: IPD_data(:)
    type(time_type), intent(in)     :: time_step

    !logical :: tsvar_out = .true.
    !logical :: wxvar_out = .false.

    type(domain2d), pointer           :: domain_coarse, domain_fine
    logical :: is_fine_pe
    integer :: parent_grid_num, child_grid_num, nz
    integer :: this_pe, n

    this_pe = mpp_pe()
    n = mygrid

    parent_grid_num = 1 
    child_grid_num = 2

    domain_fine => Atm(child_grid_num)%domain
    domain_coarse => Atm(parent_grid_num)%domain
    is_fine_pe = Atm(n)%neststruct%nested .and. ANY(Atm(n)%pelist(:) == this_pe)
    nz = Atm(n)%npz


    if (debug_log) print '("[INFO] WDR ptbounds 3 atmosphere.F90  npe=",I0," pt(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', this_pe, lbound(Atm(n)%pt,1), ubound(Atm(n)%pt,1), lbound(Atm(n)%pt,2), ubound(Atm(n)%pt,2), lbound(Atm(n)%pt,3), ubound(Atm(n)%pt,3)
    
    !if (is_fine_pe) then
    !   if (wxvar_out) call mn_prog_dump_to_netcdf(Atm(n), output_step, "wxavar", is_fine_pe, domain_coarse, domain_fine, nz)
    !   if (wxvar_out) call mn_phys_dump_to_netcdf(Atm(n), Atm_block, IPD_control, IPD_data, output_step, "wxavar", is_fine_pe, domain_coarse, domain_fine, nz)
    !   output_step = output_step + 1
    !   if (debug_log) print '("[INFO] WDR after outputting to netCDF fv_dynamics atmosphere.F90 npe=",I0, " psc=",I0)', this_pe, psc
    !end if
    
    if (a_step .lt. 10 .or. mod(a_step, 5) .eq. 0) then
    !if (mod(a_step, 20) .eq. 0) then
       if (tsvar_out) call mn_prog_dump_to_netcdf(Atm(n), a_step, "tsavar", is_fine_pe, domain_coarse, domain_fine, nz)
       if (tsvar_out) call mn_phys_dump_to_netcdf(Atm(n), Atm_block, IPD_control, IPD_data, a_step, "tsavar", is_fine_pe, domain_coarse, domain_fine, nz)
    endif
    

  end subroutine dump_moving_nest


  subroutine fv_moving_nest_init_clocks()
  
    !  --- initialize clocks for moving_nest
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
    
    
    id_movnestTot     = mpp_clock_id ('Moving Nest Total',  flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  end subroutine fv_moving_nest_init_clocks
  
  !>@brief The subroutine 'eval_move_nest' determines whether the nest should be moved
  !  and in which direction.  
  !  This is a simple prescribed motion routine; and will be replaced by code that performs 
  !  the storm tracking algorithm
  subroutine eval_move_nest(Atm, a_step, do_move, delta_i_c, delta_j_c, dt_atmos)
    type(fv_atmos_type), intent(in)   :: Atm(:)
    integer, intent(in)               :: a_step
    logical, intent(out)              :: do_move
    integer, intent(out)              :: delta_i_c, delta_j_c
    real, intent(out)              :: dt_atmos  !! only needed for the simple version of this subroutine
    
    logical, save :: first_time = .true.
    integer       :: move_incr
    logical       :: move_diag
    
    ! On the tropical channel configuration, tile 6 numbering starts at 0,0 off the coast of Spain
    !  delta_i_c = +1 is westward
    !  delta_i_c = -1 is eastward
    !
    !  delta_j_c = +1 is southward
    !  delta_j_c = -1 is northward
    
    if (dt_atmos > 100) then
       ! move once every 40 timesteps for 300s = 3 hours 20 minutes; appropriate for C96
       move_incr = 40
    else
       ! move once every 20 timesteps for 90s = 30 minutes; appropriate for C768
       move_incr = 20
    end if
    
    move_incr = 5
    
    move_diag = .false.
    
    if (move_diag) then
       ! If moving diagonal, only have to shift half as often.
       !if (a_step .eq. 1 .or. mod(a_step,2*move_incr) .eq. 0) then
       if ( mod(a_step,2*move_incr) .eq. 0) then
          do_move = .true.
          delta_i_c = 1
          delta_j_c = -1
          first_time = .false.
       else
          do_move = .false.
          delta_i_c = 0
          delta_j_c = 0
       end if
       
    else
       
       !if (a_step .eq. 1 .or. mod(a_step,2*move_incr) .eq. 0) then
       if ( mod(a_step,2*move_incr) .eq. 0) then
          do_move = .true.
          delta_i_c = 1
          delta_j_c = 0
          first_time = .false.
       else if (mod(a_step,move_incr) .eq. 0) then
          do_move = .true.
          !delta_i_c = 0
          !delta_j_c = -1
          delta_i_c = 1
          delta_j_c = 0
          first_time = .false.
       else
          do_move = .false.
          delta_i_c = 0
          delta_j_c = 0
       end if
    end if

    ! Override to prevent move on first timestep
    if (a_step .eq. 0) then
       do_move = .false.
       delta_i_c = 0
       delta_j_c = 0
    end if
    
    
  end subroutine eval_move_nest
  
  ! TODO  clarify the naming of all the grid indices;  are some of these repeated?  
  subroutine fv_moving_nest_exec(Atm, Atm_block, IPD_control, IPD_data, delta_i_c, delta_j_c, n, nest_num, parent_grid_num, child_grid_num, dt_atmos)
    type(fv_atmos_type), allocatable, target, intent(inout) :: Atm(:)
    type(block_control_type), intent(in) :: Atm_block
    type(IPD_control_type), intent(in) :: IPD_control
    type(IPD_data_type), intent(inout) :: IPD_data(:)
    integer, intent(in)  :: delta_i_c, delta_j_c
    integer, intent(in)  :: n, nest_num, parent_grid_num, child_grid_num 
    real, intent(in)     :: dt_atmos
   !  Variables from atmos_dynamics()

   integer  :: this_pe

   
   

  !---- Moving Nest local variables  -----
   integer, pointer                  :: ioffset, joffset
   real, pointer, dimension(:,:,:)   :: grid, agrid
   !type(nest_domain_type), pointer   :: nest_domain
   type(domain2d), pointer           :: domain_coarse, domain_fine
   real(kind=R_GRID), pointer, dimension(:,:,:,:) :: grid_global

   ! Constants for mpp calls
   integer  :: position      = CENTER ! CENTER, NORTH, EAST
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

   integer  :: tile_fine, tile_coarse, this_tile
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

   logical            :: filtered_terrain = .True.   ! TODO set this from namelist

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

   ! For iterating through physics/surface vector data
   integer :: nb, blen, ix, i_pe, j_pe, i_idx, j_idx


   integer, save :: output_step = 0
   
   integer, allocatable :: pelist(:)
   
   !! TODO Refine this per Atm(n) structure to allow some static and some moving nests in same run
   logical  :: is_moving_nest = .true.  ! Attach this to the namelist reading 

   real(kind=R_GRID)   :: pi = 4 * atan(1.0d0)
   real                :: rad2deg

   rad2deg = 180.0 / pi
  
   !---- Call FV dynamics -----

   gid = mpp_pe()
   this_pe = mpp_pe()


   allocate(pelist(mpp_npes()))
   call mpp_get_current_pelist(pelist)

!-----------------------------------


   !print '("[INFO] WDR NESTIDX fv_moving_nest_main.F90 npe=",I0, " n=",I0," n=",I0," nest_num=",I0," parent_grid_num=",I0," child_grid_num=",I0)', this_pe, n, n, nest_num, parent_grid_num, child_grid_num    
   
   if (first_nest_move) then
      !print '("[INFO] WDR Start Clocks npe=",I0)', this_pe
      call fv_moving_nest_init_clocks()
   end if



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

   nq = ncnst-pnats


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


    !! This tile is 1-6 for the global tiles, 7 for nest 1, 8 for nest 2, etc.
    this_tile = Atm(n)%global_tile
    
    domain_fine => Atm(child_grid_num)%domain
    parent_tile = Atm(child_grid_num)%neststruct%parent_tile
    domain_coarse => Atm(parent_grid_num)%domain

    is_fine_pe = Atm(n)%neststruct%nested .and. ANY(Atm(n)%pelist(:) == this_pe)
    nz = Atm(n)%npz
    
    if (debug_log) then
       if (is_fine_pe) then
          print '("[INFO] WDR move_nest FINE. npe=",I0, " ", I2.2, " ", I2.2," do_move=",L1," delta_i_c=",I0," delta_j_c=",I0)', this_pe, n, this_tile, do_move, delta_i_c, delta_j_c
       else
          print '("[INFO] WDR move_nest COARSE. npe=",I0, " ", I2.2, " ", I2.2)', this_pe, n, this_tile
       end if

       do nn = 1, size(Atm)
          call show_atm("1", Atm(nn), nn, this_pe)
       end do
       print '("[INFO] WDR diag Atm DONE npe=",I0," Atm(",I0,") this_tile=",I0)', this_pe, n, this_tile
    end if
       

    !if (a_step .eq. 1) then
    !   if (tsvar_out) call mn_prog_dump_to_netcdf(Atm(n), a_step-1, "tsvar", is_fine_pe, domain_coarse, domain_fine, nz)
    !   if (tsvar_out) call mn_phys_dump_to_netcdf(Atm(n), a_step-1, "tsvar", is_fine_pe, domain_coarse, domain_fine, nz)
    !end if


    !if (Atm(child_grid_num)%gridstruct%nested .and. is_moving_nest .and. do_move) then
    if (is_moving_nest .and. do_move) then
       call mpp_clock_begin (id_movnestTot)
       call mpp_clock_begin (id_movnest1)

       !!================================================================
       !! Step 1.1 -- Show the nest grids
       !!================================================================

       !if (debug_log) then
       if (debug_log .and. this_pe .eq. 0) then
          !call show_nest_grid(Atm(n), this_pe, 0)
          print '("[INFO] WDR BD init fv_moving_nest_main.F90 npe=",I0," is=",I0," ie=",I0," js=",I0," je=",I0)', this_pe, Atm(n)%bd%is,  Atm(n)%bd%ie,  Atm(n)%bd%js,  Atm(n)%bd%je
          print '("[INFO] WDR BD init fv_moving_nest_main.F90 npe=",I0," isd=",I0," ied=",I0," jsd=",I0," jed=",I0)', this_pe, Atm(n)%bd%isd,  Atm(n)%bd%ied,  Atm(n)%bd%jsd,  Atm(n)%bd%jed
          print '("[INFO] WDR BD init fv_moving_nest_main.F90 npe=",I0," isc=",I0," iec=",I0," jsc=",I0," jec=",I0)', this_pe, Atm(n)%bd%isc,  Atm(n)%bd%iec,  Atm(n)%bd%jsc, Atm(n)%bd%jec
       end if


       !!================================================================
       !! Step 1.2 -- Configure local variables
       !!================================================================

       x_refine = Atm(child_grid_num)%neststruct%refinement
       y_refine = x_refine


       if (debug_log) print '("[INFO] WDR global_nest_domain npe=",I0," tile_file=",I0," tile_coarse=",I0," istart_coarse=",I0)', this_pe, global_nest_domain%tile_fine, global_nest_domain%tile_coarse, global_nest_domain%istart_coarse


       ioffset => Atm(child_grid_num)%neststruct%ioffset
       joffset => Atm(child_grid_num)%neststruct%joffset


       if (debug_log) print '("[INFO] WDR MV_NST0 fv_moving_nest_main.F90 processing Atm(n) npe=",I0," n=",I0," ioffset=",I0," joffset=",I0)', this_pe, n, ioffset, joffset

       ! TODO update tile_fine to support multiple nests
       !tile_fine = n + 6 - 1  ! TODO: Not clear this is correct for both parent and child values of n
       tile_fine = 7
       tile_coarse = parent_tile

       istart_fine = global_nest_domain%istart_fine(nest_num)
       iend_fine = global_nest_domain%iend_fine(nest_num)
       jstart_fine = global_nest_domain%jstart_fine(nest_num)
       jend_fine = global_nest_domain%jend_fine(nest_num)

       istart_coarse = global_nest_domain%istart_coarse(nest_num)
       iend_coarse = global_nest_domain%iend_coarse(nest_num)
       jstart_coarse = global_nest_domain%jstart_coarse(nest_num)
       jend_coarse = global_nest_domain%jend_coarse(nest_num)

       ! Allocate the local weight arrays.  TODO change to use the ones from the gridstruct
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
       end if

       ! This full list of PEs is used for the mpp_sync for debugging.  Can later be removed.
       p1 = size(Atm(1)%pelist)   ! Parent PEs
       p2 = size(Atm(2)%pelist)   ! Nest PEs

       allocate(full_pelist(p1 + p2))
       do pp=1,p1
          full_pelist(pp) = Atm(1)%pelist(pp)
       end do
       do pp=1,p2
          full_pelist(p1+pp) = Atm(2)%pelist(pp)
       end do

       !!============================================================================
       !! Step 1.3 -- Dump the prognostic variables before we do the nest motion.
       !!============================================================================

       if (debug_log) print '("[INFO] WDR MV_NST0 run step 0 fv_moving_nest_main.F90 npe=",I0)', this_pe
       !if (wxvar_out) call mn_prog_dump_to_netcdf(Atm(n), output_step, "wxvar", is_fine_pe, domain_coarse, domain_fine, nz)
       !if (wxvar_out) call mn_phys_dump_to_netcdf(Atm(n), IPD_control, IPD_data, output_step, "wxvar", is_fine_pe, domain_coarse, domain_fine, nz)
       !if (wxvar_out) call mn_phys_dump_to_netcdf(Atm(n), output_step, "wxvar", is_fine_pe, domain_coarse, domain_fine, nz)
       output_step = output_step + 1

       !!============================================================================
       !! Step 1.4 -- Read in the full panel grid definition
       !!============================================================================

       if (debug_log) then
          print '("[INFO] WDR check grid_global fv_moving_nest_main.F90 npe=",I0," n=",I0)', this_pe, 1
          call check_array(Atm(1)%grid_global, this_pe, "grid_global", -2.0*3.1415926536, 2.0*3.1415926536)
          print '("[INFO] WDR check grid_global fv_moving_nest_main.F90 npe=",I0," n=",I0)', this_pe, 2
          call check_array(Atm(2)%grid_global, this_pe, "grid_global", -2.0*3.1415926536, 2.0*3.1415926536)
       end if

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
             print '("[INFO] WDR mn_latlon_read_hires_parent READING static fine file on npe=",I0)', this_pe

             call mn_latlon_read_hires_parent(Atm(1)%npx, Atm(1)%npy, x_refine, fp_super_tile_geo, &
                  Atm(child_grid_num)%neststruct%surface_dir)

             !print '("[INFO] WDR mn_orog_read_hires_parent BEFORE READING static orog fine file on npe=",I0)', this_pe
             call mn_orog_read_hires_parent(Atm(1)%npx, Atm(1)%npy, x_refine, Atm(child_grid_num)%neststruct%surface_dir, filtered_terrain, &
                  mn_static%orog_grid, mn_static%orog_std_grid, mn_static%ls_mask_grid, mn_static%land_frac_grid)
             !print '("[INFO] WDR mn_orog_read_hires_parent COMPLETED READING static orog fine file on npe=",I0)', this_pe




             call mn_static_read_hires(Atm(1)%npx, Atm(1)%npy, x_refine, trim(Atm(child_grid_num)%neststruct%surface_dir) // "/fix_sfc", "substrate_temperature", "substrate_temperature", mn_static%deep_soil_temp_grid)
             call mn_static_read_hires(Atm(1)%npx, Atm(1)%npy, x_refine, trim(Atm(child_grid_num)%neststruct%surface_dir) // "/fix_sfc", "soil_type", "soil_type", mn_static%soil_type_grid)

             !call mn_static_read_hires(Atm(1)%npx, Atm(1)%npy, x_refine, trim(Atm(child_grid_num)%neststruct%surface_dir) // "/fix_sfc", "", mn_static%veg_frac_grid)
             call mn_static_read_hires(Atm(1)%npx, Atm(1)%npy, x_refine, trim(Atm(child_grid_num)%neststruct%surface_dir) // "/fix_sfc", "vegetation_type", "vegetation_type", mn_static%veg_type_grid)

             call mn_static_read_hires(Atm(1)%npx, Atm(1)%npy, x_refine, trim(Atm(child_grid_num)%neststruct%surface_dir) // "/fix_sfc", "slope_type", "slope_type", mn_static%slope_type_grid)

             call mn_static_read_hires(Atm(1)%npx, Atm(1)%npy, x_refine, trim(Atm(child_grid_num)%neststruct%surface_dir) // "/fix_sfc", "facsf", "facsf", mn_static%facsf_grid)
             call mn_static_read_hires(Atm(1)%npx, Atm(1)%npy, x_refine, trim(Atm(child_grid_num)%neststruct%surface_dir) // "/fix_sfc", "maximum_snow_albedo", "maximum_snow_albedo", mn_static%max_snow_alb_grid)

             ! Monthly static data
             !call mn_static_read_hires_parent(Atm(1)%npx, Atm(1)%npy, x_refine, Atm(child_grid_num)%neststruct%surface_dir, "vegetation_greenness", mn_static%veg_greenness_grid)
             ! Add snowfree albedo variables here

             first_nest_move = .false.
          !else
           !  print '("[INFO] WDR mn_latlon_read_hires_parent SKIPPING static fine file on npe=",I0)', this_pe

             ! Debug outputs of hires terrain 
             !call mn_var_dump_to_netcdf(Atm(n)%phis(isd:ied, jsd:jed), is_fine_pe, domain_coarse, domain_fine, position, 1, &
             !     time_val, Atm%global_tile, "terrain", "PHIS")
             !call mn_var_dump_to_netcdf(mn_static%orog_grid(ioffset*x_refine+isd:ioffset*x_refine+ied, joffset*y_refine+jsd:joffset*y_refine+jed) * grav, &
             !     is_fine_pe, domain_coarse, domain_fine, position, 1, &
             !     time_val, Atm%global_tile, "terrain", "ORG")

          end if

          !  Validation/logging calls that can be disabled
          if (debug_log) then
             call show_tile_geo(fp_super_tile_geo, this_pe, "fp_super_tile_geo")
             call show_gridstruct(Atm(n)%gridstruct, this_pe)
             !call validate_hires_parent(fp_super_tile_geo, Atm(n)%gridstruct%grid, Atm(n)%gridstruct%agrid, &
             !     x_refine, y_refine, ioffset, joffset)
          end if
       end if

       call mpp_clock_end (id_movnest1)
       call mpp_clock_begin (id_movnest1_9)

       !!=====================================================================================           
       !! Step 1.9 -- Allocate and fill the temporary variable(s)                                              
       !!=====================================================================================           

       call mn_prog_fill_temp_variables(Atm, n, child_grid_num, is_fine_pe, npz)
       call mn_phys_fill_temp_variables(Atm, Atm_block, IPD_control, IPD_data, n, child_grid_num, is_fine_pe, npz)

       call mpp_clock_end (id_movnest1_9)
       call mpp_clock_begin (id_movnest2)

       !!============================================================================
       !! Step 2 -- Fill in the halos from the coarse grids
       !!============================================================================
       if (debug_log) print '("WDR_NEST_HALO_RECV,",I0,"===STEP 2====")', this_pe
       if (debug_log) print '("[INFO] WDR MV_NST2 run step 2 fv_moving_nest_main.F90 npe=",I0)', this_pe


       !  The halos seem to be empty at least on the first model timestep.
       !  These calls need to be executed by the parent and nest PEs in order to do the communication
       !  This is before any nest motion has occurred

       call mn_prog_fill_nest_halos_from_parent(Atm, n, child_grid_num, is_fine_pe, global_nest_domain, nz)
       call mn_phys_fill_nest_halos_from_parent(Atm, IPD_control, IPD_data, n, child_grid_num, is_fine_pe, global_nest_domain, nz)
       
       call mpp_clock_end (id_movnest2)
       call mpp_clock_begin (id_movnest3)

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
            global_nest_domain, domain_fine, domain_coarse, tile_fine, tile_coarse, &
            istart_coarse, iend_coarse, jstart_coarse, jend_coarse,  &
            istart_fine, iend_fine, jstart_fine, jend_fine)


       ! This code updates the values in neststruct; ioffset/joffset are pointers:  ioffset => Atm(child_grid_num)%neststruct%ioffset
       ioffset = ioffset + delta_i_c
       joffset = joffset + delta_j_c

       if (debug_sync) call mpp_sync(full_pelist)   ! Used to make debugging easier.  Can be removed.

       call mpp_clock_end (id_movnest3)
       call mpp_clock_begin (id_movnest4)

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
          call mn_phys_fill_intern_nest_halos(Atm(n), IPD_control, IPD_data, domain_fine, is_fine_pe)
       else
          if (debug_log) print '("[INFO] WDR MV_NST4 skip step 4 fv_moving_nest_main.F90 npe=",I0)', this_pe
       end if

       if (debug_sync) call mpp_sync(full_pelist)   ! Used to make debugging easier.  Can be removed.

       call mpp_clock_end (id_movnest4)
       call mpp_clock_begin (id_movnest5)

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
          call mn_latlon_load_parent(Atm, n, delta_i_c, delta_j_c, child_grid_num, &
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

       end if
       
       !!============================================================================
       !! Step 5.3 -- Adjust the indices by the values of delta_i_c, delta_j_c
       !!============================================================================

       call mn_shift_index(delta_i_c, delta_j_c, Atm(child_grid_num)%neststruct%ind_h)
       call mn_shift_index(delta_i_c, delta_j_c, Atm(child_grid_num)%neststruct%ind_u)
       call mn_shift_index(delta_i_c, delta_j_c, Atm(child_grid_num)%neststruct%ind_v)
       call mn_shift_index(delta_i_c, delta_j_c, Atm(child_grid_num)%neststruct%ind_b)

       if (debug_sync) call mpp_sync(full_pelist)   ! Used to make debugging easier.  Can be removed.

       call mpp_clock_end (id_movnest5)
       call mpp_clock_begin (id_movnest6)


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

       call mpp_clock_end (id_movnest6)
       call mpp_clock_begin (id_movnest7_0)

       !!=====================================================================================
       !! Step 7 --  Reset the grid definition data and buffer sizes and weights after the nest motion
       !!             Mostly needed when dynamics is executed
       !!=====================================================================================

       if (debug_log) print '("WDR_NEST_HALO_RECV,",I0,"===STEP 7====")', this_pe
       if (debug_log) print '("[INFO] WDR MV_NST7 run step 7 fv_moving_nest_main.F90 npe=",I0," n=",I0)', this_pe, n

       call mn_meta_reset_gridstruct(Atm, n, child_grid_num, global_nest_domain, fp_super_tile_geo, x_refine, y_refine, is_fine_pe, wt_h, wt_u, wt_v, a_step, dt_atmos)

       call mpp_clock_end (id_movnest7_0)
       call mpp_clock_begin (id_movnest7_1)
       
       !!=====================================================================================
       !! Step 7.01 --  Reset the orography data that was read from the hires static file
       !!
       !!=====================================================================================
       
       if (is_fine_pe) then
          ! phis is allocated in fv_arrays.F90 as:  allocate ( Atm%phis(isd:ied  ,jsd:jed  ) )
          Atm(n)%phis(isd:ied, jsd:jed) = mn_static%orog_grid((ioffset-1)*x_refine+isd:(ioffset-1)*x_refine+ied, (joffset-1)*y_refine+jsd:(joffset-1)*y_refine+jed) * grav

          ! Reinitialize diagnostics -- zsurf which is g * Atm%phis
          call fv_diag_reinit(Atm(n:n))

          ! sgh and oro were only fully allocated if fv_land is True
          !      if false, oro is (1,1), and sgh is not allocated
          if ( Atm(n)%flagstruct%fv_land ) then
             !print '("[INFO] WDR shift orography data fv_land TRUE npe=",I0)', this_pe
             ! oro and sgh are allocated only for the compute domain -- they do not have halos
          
             !fv_arrays.F90 oro() !< land fraction (1: all land; 0: all water)
             !real, _ALLOCATABLE :: oro(:,:)      _NULL  !< land fraction (1: all land; 0: all water)
             Atm(n)%oro(isc:iec, jsc:jec) = mn_static%land_frac_grid(ioffset*x_refine+isc:ioffset*x_refine+iec, joffset*y_refine+jsc:joffset*y_refine+jec)
          
             !real, _ALLOCATABLE :: sgh(:,:)      _NULL  !< Terrain standard deviation
             Atm(n)%sgh(isc:iec, jsc:jec) = mn_static%orog_std_grid(ioffset*x_refine+isc:ioffset*x_refine+iec, joffset*y_refine+jsc:joffset*y_refine+jec)
          !else
          !   print '("[INFO] WDR shift orography data fv_land FALSE npe=",I0)', this_pe
          end if

          ! Reset the land sea mask from the hires parent data
          !  Reset the variables from the fix_sfc files
          do nb = 1,Atm_block%nblks
             blen = Atm_block%blksz(nb)
             do ix = 1, blen
                i_pe = Atm_block%index(nb)%ii(ix)  
                j_pe = Atm_block%index(nb)%jj(ix)
                
                i_idx = ioffset*x_refine + i_pe
                j_idx = joffset*y_refine + j_pe
                
                IPD_data(nb)%Sfcprop%slmsk(ix) = mn_static%ls_mask_grid(i_idx, j_idx)
                
                !  IFD values are 0 for land, and 1 for oceans/lakes -- reverse of the land sea mask
                !  Land Sea Mask has values of 0 for oceans/lakes, 1 for land, 2 for sea ice
                !  TODO figure out what ifd should be for sea ice
                if (mn_static%ls_mask_grid(i_idx, j_idx) .eq. 1 ) then
                   IPD_data(nb)%Sfcprop%ifd(ix) = 0   ! Land
                else
                   IPD_data(nb)%Sfcprop%ifd(ix) = 1   ! Ocean
                end if

                IPD_data(nb)%Sfcprop%tg3(ix) = mn_static%deep_soil_temp_grid(i_idx, j_idx)
                IPD_data(nb)%Sfcprop%stype(ix) = mn_static%soil_type_grid(i_idx, j_idx)
                
                !IPD_data(nb)%Sfcprop%vfrac(ix) = mn_static%veg_frac_grid(i_idx, j_idx) 
                IPD_data(nb)%Sfcprop%vtype(ix) = mn_static%veg_type_grid(i_idx, j_idx)
                ! Add veg_greenness_grid here, monthly
                
                IPD_data(nb)%Sfcprop%slope(ix) = mn_static%slope_type_grid(i_idx, j_idx)
                
                IPD_data(nb)%Sfcprop%facsf(ix) = mn_static%facsf_grid(i_idx, j_idx)      ! fractional coverage for strong zenith angle albedo
                !IPD_data(nb)%Sfcprop%facwf(ix) = mn_static%facsf_grid(i_idx, j_idx)     ! fractional coverage for weak zenith angle albedo
                
                IPD_data(nb)%Sfcprop%snoalb(ix) = mn_static%max_snow_alb_grid(i_idx, j_idx)
                ! Add Vis/Near IR black/white sky albedo, monthly
                
                
                
             end do
          end do
       end if



       !!=====================================================================================
       !! Step 7.1   Refill the nest edge halos from parent grid after nest motion
       !!            Parent and nest PEs need to execute these subroutines
       !!=====================================================================================

       ! Refill the halos around the edge of the nest from the parent
       call mn_prog_fill_nest_halos_from_parent(Atm, n, child_grid_num, is_fine_pe, global_nest_domain, nz)
       call mn_phys_fill_nest_halos_from_parent(Atm, IPD_control, IPD_data, n, child_grid_num, is_fine_pe, global_nest_domain, nz)

       call mpp_clock_end (id_movnest7_1)

       if (is_fine_pe) then
          call mpp_clock_begin (id_movnest7_2)

          ! Refill the internal halos after nest motion
          call mn_prog_fill_intern_nest_halos(Atm(n), domain_fine, is_fine_pe)
          call mn_phys_fill_intern_nest_halos(Atm(n), IPD_control, IPD_data, domain_fine, is_fine_pe)

          call mpp_clock_end (id_movnest7_2)
       end if

       if (debug_sync) call mpp_sync(full_pelist)   ! Used to make debugging easier.  Can be removed.


       !!=====================================================================================           
       !! Step 7.3 -- Allocate and fill the temporary variable(s)                                              
       !!=====================================================================================           
       call mpp_clock_begin (id_movnest7_3)

       call mn_prog_apply_temp_variables(Atm, n, child_grid_num, is_fine_pe, npz)
       call mn_phys_apply_temp_variables(Atm, Atm_block, IPD_control, IPD_data, n, child_grid_num, is_fine_pe, npz)

       call mpp_clock_end (id_movnest7_3)
       call mpp_clock_begin (id_movnest8)


       !!============================================================================
       !!  Step 8 -- Dump to netCDF 
       !!============================================================================

       if (debug_log) print '("WDR_NEST_HALO_RECV,",I0,"===STEP 8====")', this_pe
       if (debug_log) print '("[INFO] WDR MV_NST8 run step 8 fv_moving_nest_main.F90 npe=",I0)', this_pe

       !if (wxvar_out) call mn_prog_dump_to_netcdf(Atm(n), output_step, "wxvar", is_fine_pe, domain_coarse, domain_fine, nz)
       !if (wxvar_out) call mn_phys_dump_to_netcdf(Atm(n), output_step, "wxvar", is_fine_pe, domain_coarse, domain_fine, nz)
       output_step = output_step + 1

       if (debug_sync) call mpp_sync(full_pelist)   ! Used to make debugging easier.  Can be removed.

       call mpp_clock_end (id_movnest8)
       call mpp_clock_begin (id_movnest9)


       !!=========================================================================================
       !! Step 9 -- Perform vertical remapping on nest(s) and recalculate auxiliary pressures
       !!           Should help stabilize the fields before dynamics runs
       !!=========================================================================================
       
       if (is_fine_pe) then


          if (debug_log) print '("[INFO] WDR MV_NST L2E before recalc auxiliary pressures fv_moving_nest_main.F90 npe=",I0)', this_pe
          call recalc_aux_pressures(Atm(n))
          if (debug_log) print '("[INFO] WDR MV_NST L2E after recalc auxiliary pressures fv_moving_nest_main.F90 npe=",I0)', this_pe


          !! TODO Enable vertical remapping.  Likely needed when nests move over land.
          if (debug_log) print '("[INFO] WDR MV_NST L2E before vertical remapping fv_moving_nest_main.F90 npe=",I0)', this_pe
          call vertical_remap_nest(Atm(n), dt_atmos, p_split)          
          if (debug_log) print '("[INFO] WDR MV_NST L2E after vertical remapping fv_moving_nest_main.F90 npe=",I0)', this_pe


          !! Recalculate omga; not sure if Lagrangian_to_Eulerian did this earlier
          !if( .not. Atm(n)%flagstruct%hydrostatic ) then
          !   ! !$O M P  parallel do default(none) shared(is,ie,js,je,npz,omga,delp,delz,w)
          !   do k=1,npz
          !      do j=js,je
          !         do i=is,ie
          !            Atm(n)%omga(i,j,k) = Atm(n)%delp(i,j,k)/Atm(n)%delz(i,j,k)*Atm(n)%w(i,j,k)
          !         enddo
          !      enddo
          !   enddo
          !end if

          
#ifdef REPROJ_WINDS
          ! TODO d2c_setup was removed in recent version of dycore.  Find appropriate subroutine to 
          !   interpolate A and C grid winds from the prognostic D grid winds.

          ! TODO reenable the parallel setup here
!! !$ O M P parallel do default(none) shared(isd,jsd,ied,jed,is,ie,js,je,npx,npy,npz, &
!! !$ O M P       gridstruct,flagstruct,bd,u,v,uc,vc,nested,divg) &
!! !$ O M P       private(ua,va)
          do k=1,npz
             ! Runs on individual vertical levels
             call d2c_setup(Atm(n)%u(isd,jsd,k),  Atm(n)%v(isd,jsd,k),   &
                  ua, va, &
                  Atm(n)%uc(isd,jsd,k), Atm(n)%vc(isd,jsd,k), Atm(n)%flagstruct%nord>0, &
                  isd,ied,jsd,jed, is,ie,js,je, npx,npy, &
                  Atm(n)%gridstruct%grid_type, Atm(n)%gridstruct%nested, &
                  Atm(n)%gridstruct%se_corner, Atm(n)%gridstruct%sw_corner, &
                  Atm(n)%gridstruct%ne_corner, Atm(n)%gridstruct%nw_corner, &
                  Atm(n)%gridstruct%rsin_u, Atm(n)%gridstruct%rsin_v, &
                  Atm(n)%gridstruct%cosa_s, Atm(n)%gridstruct%rsin2, Atm(n)%flagstruct%regional )

             ! These calculate divergence -- do we need?
             !if (nested) then
             !   call divergence_corner_nest(Atm(n)%u(isd,jsd,k), Atm(n)%v(isd,jsd,k), ua, va, divg(isd,jsd,k), Atm(n)%gridstruct, Atm(n)%flagstruct, Atm(n)%bd)
             !else
             !   call divergence_corner(Atm(n)%u(isd,jsd,k), Atm(n)%v(isd,jsd,k), ua, va, divg(isd,jsd,k), Atm(n)%gridstruct, Atm(n)%flagstruct, Atm(n)%bd)
             !endif
          end do



          end if

#endif REPROJ_WINDS


          if (debug_log) print '("[INFO] WDR PTVAL fv_dynamics.F90 npe=",I0," AfterNestMove ================================================")', this_pe


          !if (wxvar_out) call mn_prog_dump_to_netcdf(Atm(n), output_step, "wxvar", is_fine_pe, domain_coarse, domain_fine, nz)
          !if (wxvar_out) call mn_phys_dump_to_netcdf(Atm(n), IPD_control, IPD_data, output_step, "wxvar", is_fine_pe, domain_coarse, domain_fine, nz)
          !if (wxvar_out) call mn_phys_dump_to_netcdf(Atm(n), output_step, "wxvar", is_fine_pe, domain_coarse, domain_fine, nz)
          output_step = output_step + 1

       end if

       call mpp_clock_end (id_movnest9)
       call mpp_clock_end (id_movnestTot)


       if (debug_sync) call mpp_sync(full_pelist)   ! Used to make debugging easier.  Can be removed.

       !! Uncomment to exit and force file IO after single nest move, without dynamics
       !    call fms_io_exit()   !! Force the output of the buffered NC files
       !    if (debug_log) print '("[INFO] WDR calling mpp_exit after moving nest fv_moving_nest_main.F90 npe=",I0)', this_pe
       !    call mpp_exit()
       !    if (debug_log) print '("[INFO] WDR calling STOP after moving nest fv_moving_nest_main.F90 npe=",I0)', this_pe
       !    stop
    else
       if (debug_log) print '("[INFO] WDR move_nest not nested PE  npe=",I0)', this_pe
    end if

    !call compare_terrain("phis", Atm(n)%phis, 1, Atm(n)%neststruct%ind_h, x_refine, y_refine, is_fine_pe, global_nest_domain)
    
    if (debug_log) call show_nest_grid(Atm(n), this_pe, 99)


end subroutine fv_moving_nest_exec

#endif ! MOVING_NEST

end module fv_moving_nest_main_mod
