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
!! @brief Provides Moving Nest functionality in FV3 dynamic core
!! @author W. Ramstrom, AOML/HRD  01/15/2021
!! @email William.Ramstrom@noaa.gov
!=======================================================================!


!=======================================================================!
!
! Notes
!
!------------------------------------------------------------------------
! Moving Nest Subroutine Naming Convention
!-----------------------------------------------------------------------
!
! mn_meta_* subroutines perform moving nest operations for FV3 metadata.
!               These routines will run only once per nest move.
!
! mn_var_*  subroutines perform moving nest operations for an individual FV3 variable.
!               These routines will run many times per nest move.
!
! mn_prog_* subroutines perform moving nest operations for the list of prognostic fields.
!               These routines will run only once per nest move.
!
! mn_phys_* subroutines perform moving nest operations for the list of physics fields.
!               These routines will run only once per nest move.
!
! =======================================================================!

#define REMAP 1

module fv_moving_nest_mod
#ifdef MOVING_NEST

  use block_control_mod,      only : block_control_type
  use fms_mod,                only : mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE, clock_flag_default
  use mpp_mod,                only : mpp_pe, mpp_sync, mpp_sync_self, mpp_send, mpp_error, NOTE, FATAL
  use mpp_domains_mod,        only : mpp_update_domains, mpp_get_data_domain, mpp_get_global_domain
  use mpp_domains_mod,        only : mpp_define_nest_domains, mpp_shift_nest_domains, nest_domain_type, domain2d
  use mpp_domains_mod,        only : mpp_get_C2F_index, mpp_update_nest_fine
  use mpp_domains_mod,        only : mpp_get_F2C_index, mpp_update_nest_coarse
  use mpp_domains_mod,        only : NORTH, SOUTH, EAST, WEST, CORNER, CENTER
  use mpp_domains_mod,        only : NUPDATE, SUPDATE, EUPDATE, WUPDATE, DGRID_NE

#ifdef GFS_TYPES
  use GFS_typedefs,           only: IPD_data_type => GFS_data_type, &
      IPD_control_type => GFS_control_type, kind_phys
#else
  use IPD_typedefs,           only: IPD_data_type, IPD_control_type, kind_phys => IPD_kind_phys
#endif
  use GFS_init,               only: GFS_grid_populate

  use boundary_mod,           only: update_coarse_grid, update_coarse_grid_mpp
  use bounding_box_mod,       only: bbox, bbox_get_C2F_index, fill_bbox, show_bbox
  use constants_mod,          only: cp_air, omega, rdgas, grav, rvgas, kappa, pstd_mks, hlv
  use field_manager_mod,      only: MODEL_ATMOS
  use fms_io_mod,             only: read_data, write_data, get_global_att_value, fms_io_init, fms_io_exit
  use fv_arrays_mod,          only: fv_atmos_type, fv_nest_type, fv_grid_type, R_GRID
  use fv_arrays_mod,          only: allocate_fv_nest_bc_type, deallocate_fv_nest_bc_type
  use fv_grid_tools_mod,      only: init_grid
  use fv_grid_utils_mod,      only: grid_utils_init, ptop_min, dist2side_latlon
  use fv_mapz_mod,            only: Lagrangian_to_Eulerian, moist_cv, compute_total_energy
  use fv_moving_nest_utils_mod, only: check_array, check_local_array, show_atm, show_atm_grids, show_nest_grid, show_tile_geo, grid_equal
  use fv_nesting_mod,         only: dealloc_nested_buffers
  use fv_nwp_nudge_mod,       only: do_adiabatic_init
  use init_hydro_mod,         only: p_var
  use tracer_manager_mod,     only: get_tracer_index, get_tracer_names
  use fv_moving_nest_types_mod, only: fv_moving_nest_prog_type, fv_moving_nest_physics_type, Moving_nest
  use fv_moving_nest_utils_mod,  only: alloc_halo_buffer, load_nest_latlons_from_nc, grid_geometry, output_grid_to_nc, find_nest_alignment
  use fv_moving_nest_utils_mod,  only: fill_nest_from_buffer, fill_nest_from_buffer_cell_center, fill_nest_from_buffer_nearest_neighbor
  use fv_moving_nest_utils_mod,  only: fill_nest_halos_from_parent, fill_grid_from_supergrid, fill_weight_grid
  use fv_moving_nest_utils_mod,  only: alloc_read_data

  implicit none

#ifdef NO_QUAD_PRECISION
  ! 64-bit precision (kind=8)
  integer, parameter:: f_p = selected_real_kind(15)
#else
  ! Higher precision (kind=16) for grid geometrical factors:
  integer, parameter:: f_p = selected_real_kind(20)
#endif

#ifdef OVERLOAD_R4
  real, parameter:: real_snan=x'FFBFFFFF'
#else
  real, parameter:: real_snan=x'FFF7FFFFFFFFFFFF'
#endif

  logical :: debug_log = .false.

#include <fms_platform.h>

  !! Step 2
  interface mn_var_fill_intern_nest_halos
    module procedure mn_var_fill_intern_nest_halos_r4_2d
    module procedure mn_var_fill_intern_nest_halos_r4_3d
    module procedure mn_var_fill_intern_nest_halos_r4_4d

    module procedure mn_var_fill_intern_nest_halos_r8_2d
    module procedure mn_var_fill_intern_nest_halos_r8_3d
    module procedure mn_var_fill_intern_nest_halos_r8_4d

    module procedure mn_var_fill_intern_nest_halos_wind
  end interface mn_var_fill_intern_nest_halos


  !! Step 6
  interface mn_var_shift_data
    module procedure mn_var_shift_data_r4_2d
    module procedure mn_var_shift_data_r4_3d
    module procedure mn_var_shift_data_r4_4d

    module procedure mn_var_shift_data_r8_2d
    module procedure mn_var_shift_data_r8_3d
    module procedure mn_var_shift_data_r8_4d
  end interface mn_var_shift_data

  !! Step 8
  interface mn_var_dump_to_netcdf
    module procedure mn_var_dump_2d_to_netcdf
    module procedure mn_var_dump_3d_to_netcdf
  end interface mn_var_dump_to_netcdf

  interface mn_static_read_hires
    module procedure  mn_static_read_hires_r4
    module procedure  mn_static_read_hires_r8
  end interface mn_static_read_hires

contains

  !!=====================================================================================
  !! Step 1.9 -- Allocate and fill the temporary variable(s)
  !!            This is to manage variables that are not allocated with a halo
  !!            on the Atm structure
  !!=====================================================================================

  !>@brief The subroutine 'mn_prog_fill_temp_variables' fills the temporary variable for delz
  !>@details The delz variable does not have haloes so we need a temporary variable to move it.
  subroutine mn_prog_fill_temp_variables(Atm, n, child_grid_num, is_fine_pe, npz)
    type(fv_atmos_type), allocatable, target, intent(in)     :: Atm(:)     !< Array of atmospheric data
    integer, intent(in)                              :: n, child_grid_num  !< This level and nest level
    logical, intent(in)                              :: is_fine_pe         !< Is this the nest PE?
    integer, intent(in)                              :: npz                !< Number of vertical levels

    integer :: isd, ied, jsd, jed
    integer :: is, ie, js, je
    integer :: this_pe
    type(fv_moving_nest_prog_type), pointer :: mn_prog

    mn_prog => Moving_nest(n)%mn_prog

    this_pe = mpp_pe()

    if (debug_log) print '("[INFO] WDR start mn_prog_fill_temp_variables. npe=",I0," n=",I0)', this_pe, n

    isd = Atm(n)%bd%isd
    ied = Atm(n)%bd%ied
    jsd = Atm(n)%bd%jsd
    jed = Atm(n)%bd%jed

    if (debug_log) print '("[INFO] WDR mn_prog_fill_temp_variables. npe=",I0," isd=",I0," ied=",I0," jsd=",I0," jed=",I0)', this_pe, isd, ied, jsd, jed

    is = Atm(n)%bd%is
    ie = Atm(n)%bd%ie
    js = Atm(n)%bd%js
    je = Atm(n)%bd%je

    if (debug_log) print '("[INFO] WDR mn_prog_fill_temp_variables. npe=",I0," is=",I0," ie=",I0," js=",I0," je=",I0)', this_pe, is, ie, js, je

    ! Reset this to a dummy value, to help flag if the halos don't get updated later.
    mn_prog%delz = +99999.9
    mn_prog%delz(is:ie, js:je, 1:npz) =  Atm(n)%delz(is:ie, js:je, 1:npz)

    if (debug_log) print '("[INFO] WDR Z mn_prog_fill_temp_variables. npe=",I0," npz=",I0," ",I0," ",I0)', this_pe, npz, lbound(Atm(n)%delz,3), ubound(Atm(n)%delz,3)
    if (debug_log) print '("[INFO] WDR end mn_prog_fill_temp_variables. npe=",I0," n=",I0)', this_pe, n

  end subroutine mn_prog_fill_temp_variables

  !>@brief The subroutine 'mn_prog_apply_temp_variables' fills the Atm%delz value from the temporary variable after nest move
  !>@details The delz variable does not have haloes so we need a temporary variable to move it.
  subroutine mn_prog_apply_temp_variables(Atm, n, child_grid_num, is_fine_pe, npz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)             !< Array of atmospheric data
    integer, intent(in)                                      :: n, child_grid_num  !< This level and nest level
    logical, intent(in)                                      :: is_fine_pe         !< Is this the nest PE?
    integer, intent(in)                                      :: npz                !< Number of vertical levels

    integer :: is, ie, js, je
    integer :: this_pe
    integer :: i,j,k
    integer :: bad_values, good_values
    type(fv_moving_nest_prog_type), pointer :: mn_prog

    mn_prog => Moving_nest(n)%mn_prog

    this_pe = mpp_pe()

    if (debug_log) print '("[INFO] WDR start mn_prog_apply_temp_variables. npe=",I0," n=",I0)', this_pe, n

    ! Check if the variables were filled in properly.

    if (debug_log) then
      good_values = 0
      bad_values = 0

      if (is_fine_pe) then
        do i = Atm(n)%bd%isd, Atm(n)%bd%ied
          do j = Atm(n)%bd%jsd, Atm(n)%bd%jed
            do k = 1, npz
              if (mn_prog%delz(i,j,k) .gt. 20000.0) then
                print '("[WARN] WDR BAD NEST mn_prog%delz value. npe=",I0," mn_prog%delz(",I0,",",I0,",",I0,")=",F12.3)', this_pe, i, j, k, mn_prog%delz(i,j,k)
                bad_values = bad_values + 1
              else
                good_values = good_values + 1
              endif
            enddo
          enddo
        enddo
      else
        do i = Atm(n)%bd%is, Atm(n)%bd%ie
          do j = Atm(n)%bd%js, Atm(n)%bd%je
            do k = 1, npz
              if (mn_prog%delz(i,j,k) .gt. 20000.0) then
                print '("[WARN] WDR BAD GLOBAL mn_prog%delz value. npe=",I0," mn_prog%delz(",I0,",",I0,",",I0,")=",F12.3)', this_pe, i, j, k, mn_prog%delz(i,j,k)
                bad_values = bad_values + 1
              else
                good_values = good_values + 1
              endif
            enddo
          enddo
        enddo
      endif

      i = Atm(n)%bd%is
      j = Atm(n)%bd%js
      k = npz

      print '("[WARN] WDR Surface mn_prog%delz value. npe=",I0," mn_prog%delz(",I0,",",I0,",",I0,")=",F18.3)', this_pe, i, j, k, mn_prog%delz(i,j,k)

      print '("INFO] WDR mn_prog%delz values. npe=",I0," good_values=",I0," bad_values=",I0)', this_pe, good_values, bad_values
    endif

    if (is_fine_pe) then
      is = Atm(n)%bd%is
      ie = Atm(n)%bd%ie
      js = Atm(n)%bd%js
      je = Atm(n)%bd%je

      if (debug_log) print '("[INFO] WDR mn_prog_apply_temp_variables. npe=",I0," is=",I0," ie=",I0," js=",I0," je=",I0)', this_pe, is, ie, js, je

      Atm(n)%delz(is:ie, js:je, 1:npz) =  mn_prog%delz(is:ie, js:je, 1:npz)
    endif

    if (debug_log) print '("[INFO] WDR end mn_prog_apply_temp_variables. npe=",I0," n=",I0)', this_pe, n

  end subroutine mn_prog_apply_temp_variables


  !!=====================================================================================
  !! Step 2 -- Fill the nest edge halos from parent grid before nest motion
  !!            OR Refill the nest edge halos from parent grid after nest motion
  !!            Parent and nest PEs need to execute these subroutines
  !!=====================================================================================

  !>@brief The subroutine 'mn_prog_fill_nest_halos_from_parent' fills the nest edge halos from the parent
  !>@details Parent and nest PEs must run this subroutine.  It transfers data and interpolates onto fine nest.
  subroutine mn_prog_fill_nest_halos_from_parent(Atm, n, child_grid_num, is_fine_pe, nest_domain, nz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)             !< Array of atmospheric data
    integer, intent(in)                                      :: n, child_grid_num  !< This level and nest level
    logical, intent(in)                                      :: is_fine_pe         !< Is this the nest PE?
    type(nest_domain_type), intent(inout)                    :: nest_domain        !< Domain structure for nest
    integer, intent(in)                                      :: nz                 !< Number of vertical levels

    integer  :: position, position_u, position_v
    integer  :: interp_type, interp_type_u, interp_type_v
    integer  :: x_refine, y_refine
    type(fv_moving_nest_prog_type), pointer :: mn_prog

    mn_prog => Moving_nest(n)%mn_prog

    !  TODO Rename this from interp_type to stagger_type
    interp_type = 1    ! cell-centered A-grid
    interp_type_u = 4  ! D-grid
    interp_type_v = 4  ! D-grid

    position = CENTER
    position_u = NORTH
    position_v = EAST

    x_refine = Atm(child_grid_num)%neststruct%refinement
    y_refine = x_refine

    !  Fill centered-grid variables
    call fill_nest_halos_from_parent("q_con", Atm(n)%q_con, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
        Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)
    call fill_nest_halos_from_parent("pt", Atm(n)%pt, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
        Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)
    call fill_nest_halos_from_parent("w", Atm(n)%w, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
        Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)
    !call fill_nest_halos_from_parent("omga", Atm(n)%omga, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
    !     Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)
    call fill_nest_halos_from_parent("delp", Atm(n)%delp, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
        Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)
    call fill_nest_halos_from_parent("delz", mn_prog%delz, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
        Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)
    call fill_nest_halos_from_parent("q", Atm(n)%q, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
        Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)

    !  Move the A-grid winds.  TODO consider recomputing them from D grid instead
    call fill_nest_halos_from_parent("ua", Atm(n)%ua, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
        Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)
    call fill_nest_halos_from_parent("va", Atm(n)%va, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
        Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)

    !  Fill staggered D-grid variables
    call fill_nest_halos_from_parent("u", Atm(n)%u, interp_type_u, Atm(child_grid_num)%neststruct%wt_u, &
        Atm(child_grid_num)%neststruct%ind_u, x_refine, y_refine, is_fine_pe, nest_domain, position_u, nz)
    call fill_nest_halos_from_parent("v", Atm(n)%v, interp_type_v, Atm(child_grid_num)%neststruct%wt_v, &
        Atm(child_grid_num)%neststruct%ind_v, x_refine, y_refine, is_fine_pe, nest_domain, position_v, nz)

  end subroutine mn_prog_fill_nest_halos_from_parent

  !!============================================================================
  !! Step 3 -- Redefine the nest domain to new location
  !!   This calls mpp_shift_nest_domains.
  !!  --  Similar to med_nest_configure() from HWRF
  !!============================================================================

  !>@brief The subroutine 'mn_meta_move_nest' resets the metadata for the nest
  !>@details Parent and  nest PEs run this subroutine.
  subroutine mn_meta_move_nest(delta_i_c, delta_j_c, pelist, is_fine_pe, extra_halo, nest_domain, domain_fine, domain_coarse, &
      istart_coarse, iend_coarse, jstart_coarse, jend_coarse,  istart_fine, iend_fine, jstart_fine, jend_fine)

    implicit none

    integer, intent(in)                   :: delta_i_c, delta_j_c                                    !< Coarse grid delta i,j for nest move
    integer, allocatable, intent(in)      :: pelist(:)                                               !< List of involved PEs
    logical, intent(in)                   :: is_fine_pe                                              !< Is this a nest PE?
    integer, intent(in)                   :: extra_halo                                              !< Extra halo points (not fully implemented)
    type(nest_domain_type), intent(inout) :: nest_domain                                             !< Nest domain structure
    type(domain2d), intent(inout)         :: domain_coarse, domain_fine                              !< Coarse and fine domain structures
    integer, intent(inout)                :: istart_coarse, iend_coarse, jstart_coarse, jend_coarse  !< Bounds of coarse grid
    integer, intent(in)                   :: istart_fine, iend_fine, jstart_fine, jend_fine          !< Bounds of fine grid

    !  Local variables
    integer   :: num_nest
    integer   :: this_pe

    integer   :: delta_i_coarse(1), delta_j_coarse(1)

    this_pe = mpp_pe()

    if (debug_log) print '("[INFO] WDR start mn_meta_move_nest. npe=",I0)', this_pe

    !  Initial implementation only supports single moving nest.  Update this later.
    !  mpp_shift_nest_domains has a call signature to support multiple moving nests, though has not been tested for correctness.
    delta_i_coarse(1) = delta_i_c
    delta_j_coarse(1) = delta_j_c

    !!===========================================================
    !!
    !! Relocate where the nest is aligned on the parent
    !!
    !!===========================================================

    if (debug_log) print '("[INFO] WDR NRD0. npe=",I0," ",I0," ",I0," ",I0," ",I0," num_nest=",I0," delta_i_c=",I0," delta_j_c=",I0)', this_pe, istart_coarse, iend_coarse, jstart_coarse, jend_coarse, num_nest, delta_i_c, delta_j_c

    istart_coarse = istart_coarse + delta_i_c
    iend_coarse = iend_coarse + delta_i_c

    jstart_coarse = jstart_coarse + delta_j_c
    jend_coarse = jend_coarse + delta_j_c

    ! The fine nest will maintain the same indices

    num_nest = nest_domain%num_nest

    if (debug_log) print '("[INFO] WDR NRD1 about to call mpp_shift_nest_domains. npe=",I0," ",I0," ",I0," ",I0," ",I0," num_nest=",I0," delta_i_c=",I0," delta_j_c=",I0)', this_pe, istart_coarse, iend_coarse, jstart_coarse, jend_coarse, num_nest, delta_i_c, delta_j_c


    ! WDR TODO Verify whether rerunning this will cause (small) memory leaks.
    if (is_fine_pe) then
      call mpp_shift_nest_domains(nest_domain, domain_fine, delta_i_coarse, delta_j_coarse, extra_halo)
    else
      call mpp_shift_nest_domains(nest_domain, domain_coarse, delta_i_coarse, delta_j_coarse, extra_halo)
    endif

    if (debug_log) print '("[INFO] WDR NRD2 after call to mpp_define_nest_domains. npe=",I0)', this_pe

  end subroutine mn_meta_move_nest


  !================================================================================
  !! Step 4 --  Updates the internal nest tile halos
  !================================================================================

  !>@brief The subroutine 'mn_prog_fill_intern_nest_halos' fill internal nest halos for prognostic variables
  !>@details Only nest PEs call this subroutine.
  subroutine mn_prog_fill_intern_nest_halos(Atm, domain_fine, is_fine_pe)
    type(fv_atmos_type), target, intent(inout)  :: Atm           !< Single instance of atmospheric data
    type(domain2d), intent(inout)               :: domain_fine   !< Domain structure for nest
    logical, intent(in)                         :: is_fine_pe    !< Is this a nest PE?

    integer :: this_pe
    type(fv_moving_nest_prog_type), pointer :: mn_prog

    mn_prog => Moving_nest(2)%mn_prog  ! TODO allow nest number to vary
    this_pe = mpp_pe()

    call mn_var_fill_intern_nest_halos(Atm%q_con, domain_fine, is_fine_pe)
    call mn_var_fill_intern_nest_halos(Atm%pt, domain_fine, is_fine_pe)
    call mn_var_fill_intern_nest_halos(Atm%w, domain_fine, is_fine_pe)
    !call mn_var_fill_intern_nest_halos(Atm%omga, domain_fine, is_fine_pe)
    call mn_var_fill_intern_nest_halos(Atm%delp, domain_fine, is_fine_pe)
    call mn_var_fill_intern_nest_halos(mn_prog%delz, domain_fine, is_fine_pe)

    call mn_var_fill_intern_nest_halos(Atm%ua, domain_fine, is_fine_pe)
    call mn_var_fill_intern_nest_halos(Atm%va, domain_fine, is_fine_pe)

    if (debug_log) then
      call check_array(Atm%u, this_pe, "Atm%u", -300.0, 300.0)
      call check_array(Atm%v, this_pe, "Atm%v", -300.0, 300.0)
    endif

    ! The vector form of the subroutine takes care of the staggering of the wind variables internally.
    call mn_var_fill_intern_nest_halos(Atm%u, Atm%v, domain_fine, is_fine_pe)

    call mn_var_fill_intern_nest_halos(Atm%q, domain_fine, is_fine_pe)

  end subroutine mn_prog_fill_intern_nest_halos


  !================================================================================
  !
  !   Step 4 -- Per variable fill internal nest halos
  !
  !================================================================================

  !>@brief The subroutine 'mn_var_fill_intern_nest_halos_r4_2d' fills internal nest halos
  !>@details This version of the subroutine is for 2D arrays of single precision reals.
  subroutine mn_var_fill_intern_nest_halos_r4_2d(data_var, domain_fine, is_fine_pe)
    real*4, allocatable, intent(inout)          :: data_var(:,:)  !< Model variable data
    type(domain2d), intent(inout)               :: domain_fine    !< Nest domain structure
    logical, intent(in)                         :: is_fine_pe     !< Is this the nest PE?

    integer                      :: this_pe
    this_pe = mpp_pe()

    if (is_fine_pe) then
      if (debug_log) print '("[INFO] WDR INH2 before call to mpp_update_domains. npe=",I0)', this_pe
      ! mpp_update_domains fills the halo region of the fine grids for the interior of the nest.
      ! The fine nest boundary with the coarse grid remains unchanged.
      ! seems that this only performs communication between fine nest PEs
      ! Just transfers halo data between tiles of same resolution -- doesn't perform any interpolation!
      call mpp_update_domains(data_var, domain_fine,  flags=NUPDATE + EUPDATE + SUPDATE + WUPDATE)
      if (debug_log) print '("[INFO] WDR INH2 after call to mpp_update_domains. npe=",I0)', this_pe
    endif

  end subroutine mn_var_fill_intern_nest_halos_r4_2d

  !>@brief The subroutine 'mn_var_fill_intern_nest_halos_r8_2d' fills internal nest halos
  !>@details This version of the subroutine is for 2D arrays of double precision reals.
  subroutine mn_var_fill_intern_nest_halos_r8_2d(data_var, domain_fine, is_fine_pe)
    real*8, allocatable, intent(inout)          :: data_var(:,:)  !< Double precision model variable
    type(domain2d), intent(inout)               :: domain_fine    !< Nest domain structure
    logical, intent(in)                         :: is_fine_pe     !< Is this a nest PE?

    integer                      :: this_pe
    this_pe = mpp_pe()

    if (is_fine_pe) then
      if (debug_log) print '("[INFO] WDR INH2p before call to mpp_update_domains. npe=",I0)', this_pe
      call mpp_update_domains(data_var, domain_fine,  flags=NUPDATE + EUPDATE + SUPDATE + WUPDATE)
      if (debug_log) print '("[INFO] WDR INH2p after call to mpp_update_domains. npe=",I0)', this_pe
    endif

  end subroutine mn_var_fill_intern_nest_halos_r8_2d

  !>@brief The subroutine 'mn_var_fill_intern_nest_halos_r4_3d' fills internal nest halos
  !>@details This version of the subroutine is for 3D arrays of single precision reals.
  subroutine mn_var_fill_intern_nest_halos_r4_3d(data_var, domain_fine, is_fine_pe)
    real*4, allocatable, intent(inout)          :: data_var(:,:,:)  !< Single precision model variable
    type(domain2d), intent(inout)               :: domain_fine      !< Nest domain structure
    logical, intent(in)                         :: is_fine_pe       !< Is this a nest PE?

    integer                      :: this_pe
    this_pe = mpp_pe()

    if (is_fine_pe) then
      if (debug_log) print '("[INFO] WDR INH3 before call to mpp_update_domains. npe=",I0)', this_pe
      call mpp_update_domains(data_var, domain_fine,  flags=NUPDATE + EUPDATE + SUPDATE + WUPDATE)
      if (debug_log) print '("[INFO] WDR INH3 after call to mpp_update_domains. npe=",I0)', this_pe
    endif

  end subroutine mn_var_fill_intern_nest_halos_r4_3d

  !>@brief The subroutine 'mn_var_fill_intern_nest_halos_r8_3d' fills internal nest halos
  !>@details This version of the subroutine is for 3D arrays of double precision reals.
  subroutine mn_var_fill_intern_nest_halos_r8_3d(data_var, domain_fine, is_fine_pe)
    real*8, allocatable, intent(inout)          :: data_var(:,:,:)  !< Double precision model variable
    type(domain2d), intent(inout)               :: domain_fine      !< Nest domain structure
    logical, intent(in)                         :: is_fine_pe       !< Is this a nest PE?

    integer                      :: this_pe
    this_pe = mpp_pe()

    if (is_fine_pe) then
      if (debug_log) print '("[INFO] WDR INH3p before call to mpp_update_domains. npe=",I0)', this_pe
      call mpp_update_domains(data_var, domain_fine,  flags=NUPDATE + EUPDATE + SUPDATE + WUPDATE)
      if (debug_log) print '("[INFO] WDR INH3p after call to mpp_update_domains. npe=",I0)', this_pe
    endif

  end subroutine mn_var_fill_intern_nest_halos_r8_3d

  !>@brief The subroutine 'mn_var_fill_intern_nest_halos_wind' fills internal nest halos for u and v wind
  !>@details This version of the subroutine is for 3D arrays of single precision reals for each wind component
  subroutine mn_var_fill_intern_nest_halos_wind(u_var, v_var, domain_fine, is_fine_pe)
    real, allocatable, intent(inout)            :: u_var(:,:,:) !< Staggered u wind
    real, allocatable, intent(inout)            :: v_var(:,:,:) !< Staggered v wind
    type(domain2d), intent(inout)               :: domain_fine  !< Nest domain structure
    logical, intent(in)                         :: is_fine_pe   !< Is this a nest PE?

    integer                      :: this_pe
    this_pe = mpp_pe()

    if (is_fine_pe) then
      if (debug_log) print '("[INFO] WDR INH3W before call to mpp_update_domains. npe=",I0)', this_pe
      call mpp_update_domains(u_var, v_var, domain_fine,  flags=NUPDATE + EUPDATE + SUPDATE + WUPDATE, gridtype=DGRID_NE)
      if (debug_log) print '("[INFO] WDR INH3W after call to mpp_update_domains. npe=",I0)', this_pe
    endif

  end subroutine mn_var_fill_intern_nest_halos_wind


  !>@brief The subroutine 'mn_var_fill_intern_nest_halos_r4_4d' fills internal nest halos
  !>@details This version of the subroutine is for 4D arrays of single precision reals.
  subroutine mn_var_fill_intern_nest_halos_r4_4d(data_var, domain_fine, is_fine_pe)
    real*4, allocatable, intent(inout)          :: data_var(:,:,:,:)  !< Single prevision variable
    type(domain2d), intent(inout)               :: domain_fine        !< Nest domain structure
    logical, intent(in)                         :: is_fine_pe         !< Is this a nest PE?

    integer                      :: this_pe
    this_pe = mpp_pe()

    if (is_fine_pe) then
      if (debug_log) print '("[INFO] WDR INH4 before call to mpp_update_domains. npe=",I0)', this_pe
      call mpp_update_domains(data_var, domain_fine,  flags=NUPDATE + EUPDATE + SUPDATE + WUPDATE)
      if (debug_log) print '("[INFO] WDR INH4 after call to mpp_update_domains. npe=",I0)', this_pe
    endif

  end subroutine mn_var_fill_intern_nest_halos_r4_4d

  !>@brief The subroutine 'mn_var_fill_intern_nest_halos_r8_4d' fills internal nest halos
  !>@details This version of the subroutine is for 4D arrays of double precision reals.
  subroutine mn_var_fill_intern_nest_halos_r8_4d(data_var, domain_fine, is_fine_pe)
    real*8, allocatable, intent(inout)          :: data_var(:,:,:,:)  !< Double precision variable
    type(domain2d), intent(inout)               :: domain_fine        !< Nest domain structure
    logical, intent(in)                         :: is_fine_pe         !< Is this a nest PE?

    integer                      :: this_pe
    this_pe = mpp_pe()

    if (is_fine_pe) then
      if (debug_log) print '("[INFO] WDR INH4 before call to mpp_update_domains. npe=",I0)', this_pe
      call mpp_update_domains(data_var, domain_fine,  flags=NUPDATE + EUPDATE + SUPDATE + WUPDATE)
      if (debug_log) print '("[INFO] WDR INH4 after call to mpp_update_domains. npe=",I0)', this_pe
    endif

  end subroutine mn_var_fill_intern_nest_halos_r8_4d


  !!============================================================================
  !! Step 5.1 -- Load the latlon data from NetCDF
  !!             update parent_geo, tile_geo*, p_grid*, n_grid*
  !!============================================================================

  !>@brief The subroutine 'mn_latlon_load_parent' loads parent latlon data from netCDF
  !>@details Updates parent_geo, tile_geo*, p_grid*, n_grid*
  subroutine mn_latlon_load_parent(surface_dir, Atm, n, parent_tile, delta_i_c, delta_j_c, child_grid_num, parent_geo, tile_geo, tile_geo_u, tile_geo_v, fp_super_tile_geo, p_grid, n_grid, p_grid_u, n_grid_u, p_grid_v, n_grid_v)
    character(len=*), intent(in)                 :: surface_dir                                   !< Directory for static files
    type(fv_atmos_type), allocatable, intent(in) :: Atm(:)                                        !< Atm data array
    integer, intent(in)                          :: n, parent_tile, child_grid_num                !< Grid numbers
    integer, intent(in)                          :: delta_i_c, delta_j_c                          !< Nest motion in delta i,j
    type(grid_geometry), intent(inout)           :: parent_geo, tile_geo, tile_geo_u, tile_geo_v  !< Tile geometries
    type(grid_geometry), intent(in)              :: fp_super_tile_geo                             !< Parent grid at high-resolution geometry
    real(kind=R_GRID), allocatable, intent(out)  :: p_grid(:,:,:), n_grid(:,:,:)                  !< A-stagger lat/lon grids
    real(kind=R_GRID), allocatable, intent(out)  :: p_grid_u(:,:,:), n_grid_u(:,:,:)              !< u-wind staggered lat/lon grids
    real(kind=R_GRID), allocatable, intent(out)  :: p_grid_v(:,:,:), n_grid_v(:,:,:)              !< v-wind staggered lat/lon grids

    character(len=256) :: grid_filename
    logical, save  :: first_nest_move = .true.
    integer, save  :: p_istart_fine, p_iend_fine, p_jstart_fine, p_jend_fine
    integer :: x, y, fp_i, fp_j
    integer :: position, position_u, position_v
    integer :: x_refine, y_refine
    integer :: nest_x, nest_y, parent_x, parent_y
    integer :: this_pe

    this_pe = mpp_pe()

    position = CENTER
    position_u = NORTH
    position_v = EAST

    x_refine = Atm(child_grid_num)%neststruct%refinement
    y_refine = x_refine

    !  Setup parent_geo with the values for the parent tile
    !  Note that lat/lon are stored in the model in RADIANS
    !  Only the netCDF files use degrees

    if (first_nest_move) then
      if (debug_log) print '("[INFO] WDR mn_latlon_load_parent READING static coarse file on npe=",I0)', this_pe

      call mn_static_filename(surface_dir, parent_tile, 'grid', 1, grid_filename)
      call load_nest_latlons_from_nc(grid_filename, Atm(1)%npx, Atm(1)%npy, 1, &
          parent_geo, p_istart_fine, p_iend_fine, p_jstart_fine, p_jend_fine)

      first_nest_move = .false.
    endif

    parent_geo%nxp = Atm(1)%npx
    parent_geo%nyp = Atm(1)%npy

    parent_geo%nx = Atm(1)%npx - 1
    parent_geo%ny = Atm(1)%npy - 1

    if (debug_log) then
      call show_tile_geo(parent_geo, this_pe, "parent_geo")
      call show_atm_grids(Atm, n)
    endif

    !===========================================================
    !  Begin tile_geo per PE.
    !===========================================================

    !------------------------
    ! Grid Definitions
    !------------------------
    !
    ! tile_geo - lat/lons on A-grid (cell centers) for nest, on data domain (includes halo) for each PE
    ! parent_geo - lat/lons of supergrid for parent
    ! n_grid - lat/lons of cell centers for nest
    ! p_grid - lat/lons of cell centers for parent
    !
    ! gridstruct%agrid - cell centers for each PE
    ! gridstruct%grid - cell corners for each PE

    ! Allocate tile_geo just for this PE, copied from Atm(n)%gridstruct%agrid
    tile_geo%nx = ubound(Atm(n)%gridstruct%agrid, 1) - lbound(Atm(n)%gridstruct%agrid, 1)
    tile_geo%ny = ubound(Atm(n)%gridstruct%agrid, 2) - lbound(Atm(n)%gridstruct%agrid, 2)
    tile_geo%nxp = tile_geo%nx + 1
    tile_geo%nyp = tile_geo%ny + 1

    allocate(tile_geo%lons(lbound(Atm(n)%gridstruct%agrid, 1):ubound(Atm(n)%gridstruct%agrid, 1), lbound(Atm(n)%gridstruct%agrid, 2):ubound(Atm(n)%gridstruct%agrid, 2)))
    allocate(tile_geo%lats(lbound(Atm(n)%gridstruct%agrid, 1):ubound(Atm(n)%gridstruct%agrid, 1), lbound(Atm(n)%gridstruct%agrid, 2):ubound(Atm(n)%gridstruct%agrid, 2)))

    tile_geo%lats = -999.9
    tile_geo%lons = -999.9

    do x = lbound(Atm(n)%gridstruct%agrid, 1), ubound(Atm(n)%gridstruct%agrid, 1)
      do y = lbound(Atm(n)%gridstruct%agrid, 2), ubound(Atm(n)%gridstruct%agrid, 2)
        tile_geo%lons(x,y) = Atm(n)%gridstruct%agrid(x,y,1)
        tile_geo%lats(x,y) = Atm(n)%gridstruct%agrid(x,y,2)
      enddo
    enddo

    if (debug_log) call show_tile_geo(tile_geo, this_pe, "tile_geo")
    call find_nest_alignment(tile_geo, fp_super_tile_geo, nest_x, nest_y, parent_x, parent_y)

    if (parent_x .eq. -999) then
      print '("[ERROR] WDR mn_latlon_load_parent on npe=",I0," parent and nest grids are not aligned!")', this_pe
      call mpp_error(FATAL, "mn_latlon_load_parent parent and nest grids are not aligned.")
    endif

    ! Allocate tile_geo_u just for this PE, copied from Atm(n)%gridstruct%grid
    ! grid is 1 larger than agrid
    ! u(npx, npy+1)
    tile_geo_u%nx = ubound(Atm(n)%gridstruct%agrid, 1) - lbound(Atm(n)%gridstruct%agrid, 1)
    tile_geo_u%ny = ubound(Atm(n)%gridstruct%grid, 2) - lbound(Atm(n)%gridstruct%grid, 2)
    tile_geo_u%nxp = tile_geo_u%nx + 1
    tile_geo_u%nyp = tile_geo_u%ny + 1

    allocate(tile_geo_u%lons(lbound(Atm(n)%gridstruct%agrid, 1):ubound(Atm(n)%gridstruct%agrid, 1), lbound(Atm(n)%gridstruct%grid, 2):ubound(Atm(n)%gridstruct%grid, 2)))
    allocate(tile_geo_u%lats(lbound(Atm(n)%gridstruct%agrid, 1):ubound(Atm(n)%gridstruct%agrid, 1), lbound(Atm(n)%gridstruct%grid, 2):ubound(Atm(n)%gridstruct%grid, 2)))

    tile_geo_u%lons = -999.9
    tile_geo_u%lats = -999.9

    do x = lbound(tile_geo_u%lats, 1), ubound(tile_geo_u%lats, 1)
      do y = lbound(tile_geo_u%lats, 2), ubound(tile_geo_u%lats, 2)
        fp_i = (x - nest_x) * 2 + parent_x - 1
        fp_j = (y - nest_y) * 2 + parent_y

        !print '("[INFO] WDR mn_latlon_load_parent on npe=",I0," fp_i=",I0," fp_j=",I0,4I6)', this_pe, fp_i, fp_j, nest_x, nest_y, parent_x, parent_y

        tile_geo_u%lons(x,y) = fp_super_tile_geo%lons(fp_i, fp_j)
        tile_geo_u%lats(x,y) = fp_super_tile_geo%lats(fp_i, fp_j)
      enddo
    enddo

    if (debug_log) call show_tile_geo(tile_geo_u, this_pe, "tile_geo_u")

    ! Allocate tile_geo_v just for this PE, copied from Atm(n)%gridstruct%grid
    ! grid is 1 larger than agrid
    ! u(npx, npy+1)
    tile_geo_v%nx = ubound(Atm(n)%gridstruct%grid, 1) - lbound(Atm(n)%gridstruct%grid, 1)
    tile_geo_v%ny = ubound(Atm(n)%gridstruct%agrid, 2) - lbound(Atm(n)%gridstruct%agrid, 2)
    tile_geo_v%nxp = tile_geo_v%nx + 1
    tile_geo_v%nyp = tile_geo_v%ny + 1

    allocate(tile_geo_v%lons(lbound(Atm(n)%gridstruct%grid, 1):ubound(Atm(n)%gridstruct%grid, 1), lbound(Atm(n)%gridstruct%agrid, 2):ubound(Atm(n)%gridstruct%agrid, 2)))
    allocate(tile_geo_v%lats(lbound(Atm(n)%gridstruct%grid, 1):ubound(Atm(n)%gridstruct%grid, 1), lbound(Atm(n)%gridstruct%agrid, 2):ubound(Atm(n)%gridstruct%agrid, 2)))

    tile_geo_v%lons = -999.9
    tile_geo_v%lats = -999.9

    do x = lbound(tile_geo_v%lats, 1), ubound(tile_geo_v%lats, 1)
      do y = lbound(tile_geo_v%lats, 2), ubound(tile_geo_v%lats, 2)
        fp_i = (x - nest_x) * 2 + parent_x
        fp_j = (y - nest_y) * 2 + parent_y - 1

        tile_geo_v%lons(x,y) = fp_super_tile_geo%lons(fp_i, fp_j)
        tile_geo_v%lats(x,y) = fp_super_tile_geo%lats(fp_i, fp_j)
      enddo
    enddo

    if (debug_log) call show_tile_geo(tile_geo_v, this_pe, "tile_geo_v")

    !===========================================================
    !  End tile_geo per PE.
    !===========================================================

    allocate(p_grid(1:parent_geo%nxp, 1:parent_geo%nyp,2))
    allocate(n_grid(Atm(child_grid_num)%bd%isd:Atm(child_grid_num)%bd%ied, Atm(child_grid_num)%bd%jsd:Atm(child_grid_num)%bd%jed, 2))
    n_grid = real_snan

    allocate(p_grid_u(1:parent_geo%nxp, 1:parent_geo%nyp+1,2))
    allocate(n_grid_u(Atm(child_grid_num)%bd%isd:Atm(child_grid_num)%bd%ied, Atm(child_grid_num)%bd%jsd:Atm(child_grid_num)%bd%jed+1, 2))
    n_grid_u = real_snan

    allocate(p_grid_v(1:parent_geo%nxp+1, 1:parent_geo%nyp,2))
    allocate(n_grid_v(Atm(child_grid_num)%bd%isd:Atm(child_grid_num)%bd%ied+1, Atm(child_grid_num)%bd%jsd:Atm(child_grid_num)%bd%jed, 2))
    n_grid_v = real_snan

    ! TODO - propagate tile_geo information back to Atm structure
    ! TODO - deallocate tile_geo lat/lons
    ! TODO - ensure the allocation of tile_geo lat/lons is only performed once - outside the loop

    if (debug_log) print '("[INFO] WDR MV_NST2 run step 2 atmosphere.F90 npe=",I0, " tile_geo: nxp=",I0," nyp=",I0," nx=",I0," ny=", I0)', this_pe, tile_geo%nxp, tile_geo%nyp, tile_geo%nx, tile_geo%ny
    if (debug_log) print *, "[INFO] WDR MV_NST2 run step 2 atmosphere.F90 shape(tile_geo%lats)=", shape(tile_geo%lats)
    if (debug_log) print '("[INFO] WDR MV_NST2 bounds1 (tile_geo%lats)=",I0,"-",I0)', lbound(tile_geo%lats,1), ubound(tile_geo%lats,1)
    if (debug_log) print '("[INFO] WDR MV_NST2 bounds2 (tile_geo%lats)=",I0,"-",I0)', lbound(tile_geo%lats,2), ubound(tile_geo%lats,2)

    call move_nest_geo(tile_geo, tile_geo_u, tile_geo_v, fp_super_tile_geo, delta_i_c, delta_j_c, x_refine, y_refine)

    call assign_n_p_grids(parent_geo, tile_geo, p_grid, n_grid, position)
    call assign_n_p_grids(parent_geo, tile_geo_u, p_grid_u, n_grid_u, position_u)
    call assign_n_p_grids(parent_geo, tile_geo_v, p_grid_v, n_grid_v, position_v)

  end subroutine mn_latlon_load_parent

  !>@brief The subroutine 'mn_static_filename' generates the full pathname for a static file for each run
  !>@details Constructs the full pathname for a variable and refinement level and tests whether it exists
  subroutine mn_static_filename(surface_dir, tile_num, tag, refine, grid_filename)
    character(len=*), intent(in)       :: surface_dir     !< Directory
    character(len=*), intent(in)       :: tag             !< Variable name
    integer, intent(in)                :: tile_num        !< Tile number
    integer, intent(in)                :: refine          !< Nest refinement
    character(len=*), intent(out)      :: grid_filename   !< Output pathname to netCDF file

    character(len=256) :: refine_str, parent_str
    character(len=1)   :: divider
    logical            :: file_exists

    write(parent_str, '(I0)'), tile_num

    if (refine .eq. 1 .and. (tag .eq. 'grid' .or. tag .eq. 'oro_data')) then
      ! For 1x files in INPUT directory; go at the symbolic link
      grid_filename = trim(trim(surface_dir) // '/' // trim(tag) // '.tile' // trim(parent_str) // '.nc')
    else
      if (refine .eq. 1) then
        grid_filename = trim(trim(surface_dir) // '/' // trim(tag) // '.tile' // trim(parent_str) // '.nc')
      else
        write(refine_str, '(I0,A1)'), refine, 'x'
        grid_filename = trim(trim(surface_dir) // '/' // trim(tag) // '.tile' // trim(parent_str) // '.' // trim(refine_str) // '.nc')
      endif
    endif

    grid_filename = trim(grid_filename)

    inquire(FILE=grid_filename, EXIST=file_exists)
    if (.not. file_exists) then
      print '("[ERROR] WDR mn_static_filename DOES NOT EXIST npe=",I0," exists="L1," ",A256)', mpp_pe(), file_exists, grid_filename
    endif

  end subroutine mn_static_filename

  !>@brief The subroutine 'mn_latlon_read_hires_parent' reads in static data from a netCDF file
  subroutine mn_latlon_read_hires_parent(npx, npy, refine, fp_super_tile_geo, surface_dir, parent_tile)
    integer, intent(in)                :: npx, npy, refine     !< Number of points in x,y, and refinement
    type(grid_geometry), intent(inout) :: fp_super_tile_geo    !< Geometry of supergrid for parent tile at high resolution
    character(len=*), intent(in)       :: surface_dir          !< Surface directory to read netCDF file from
    integer, intent(in)                :: parent_tile          !< Parent tile number

    integer                            :: fp_super_istart_fine, fp_super_jstart_fine,fp_super_iend_fine, fp_super_jend_fine
    character(len=256)                 :: grid_filename

    call mn_static_filename(surface_dir, parent_tile, 'grid',  refine, grid_filename)

    call load_nest_latlons_from_nc(trim(grid_filename), npx, npy, refine, fp_super_tile_geo, &
        fp_super_istart_fine, fp_super_iend_fine, fp_super_jstart_fine, fp_super_jend_fine)

  end subroutine mn_latlon_read_hires_parent

  !>@brief The subroutine 'mn_orog_read_hires_parent' loads parent orography data from netCDF
  !>@details Gathers a number of terrain-related variables from the netCDF file
  subroutine mn_orog_read_hires_parent(npx, npy, refine, surface_dir, filtered_terrain, orog_grid, orog_std_grid, ls_mask_grid, land_frac_grid, parent_tile)
    integer, intent(in)                :: npx, npy, refine   !< Number of points in x,y, and refinement
    character(len=*), intent(in)       :: surface_dir        !< Surface directory to read netCDF file from
    logical, intent(in)                :: filtered_terrain   !< Whether to use filtered terrain
    real, allocatable, intent(out)     :: orog_grid(:,:)     !< Output orography grid
    real, allocatable, intent(out)     :: orog_std_grid(:,:) !< Output orography standard deviation grid
    real, allocatable, intent(out)     :: ls_mask_grid(:,:)  !< Output land sea mask grid
    real, allocatable, intent(out)     :: land_frac_grid(:,:)!< Output land fraction grid
    integer, intent(in)                :: parent_tile        !< Parent tile number

    integer :: nx_cubic, nx, ny, fp_nx, fp_ny, mid_nx, mid_ny
    integer :: fp_istart_fine, fp_iend_fine, fp_jstart_fine, fp_jend_fine
    character(len=512) :: nc_filename
    character(len=16)  :: orog_var_name
    integer :: this_pe

    this_pe = mpp_pe()

    nx_cubic = npx - 1
    nx = npx - 1
    ny = npy - 1

    fp_istart_fine = 0
    fp_iend_fine = nx * refine
    fp_jstart_fine = 0
    fp_jend_fine = ny * refine

    fp_nx = fp_iend_fine - fp_istart_fine
    fp_ny = fp_jend_fine - fp_jstart_fine

    mid_nx = (fp_iend_fine - fp_istart_fine) / 2
    mid_ny = (fp_jend_fine - fp_jstart_fine) / 2

    call mn_static_filename(surface_dir, parent_tile, 'oro_data', refine, nc_filename)

    if (filtered_terrain) then
      orog_var_name = 'orog_filt'
    else
      orog_var_name = 'orog_raw'
    endif

    if (debug_log) print '("[INFO] WDR NCREAD LOFC mn_orog_read_hires_parent npe=",I0,I4,I4,I4,I4," ",A12," ",A128)', this_pe, fp_nx, fp_ny, mid_nx,mid_ny, orog_var_name, nc_filename

    call alloc_read_data(nc_filename, orog_var_name, fp_nx, fp_ny, orog_grid)
    !call check_array(orog_grid, this_pe, "parent coarse" // orog_var_name, -1000.0, 5000.0)
    call alloc_read_data(nc_filename, 'slmsk', fp_nx, fp_ny, ls_mask_grid)
    !call check_array(ls_mask_grid, this_pe, 'slmsk', 0.0, 3.0)

    call alloc_read_data(nc_filename, 'stddev', fp_nx, fp_ny, orog_std_grid)      ! TODO validate if this is needed
    call alloc_read_data(nc_filename, 'land_frac', fp_nx, fp_ny, land_frac_grid)  ! TODO validate if this is needed

  end subroutine mn_orog_read_hires_parent

  !>@brief The subroutine 'mn_static_read_hires_r4' loads high resolution data from netCDF
  !>@details Gathers a single variable from the netCDF file
  subroutine mn_static_read_hires_r4(npx, npy, refine, surface_dir, file_prefix, var_name, data_grid, parent_tile, time)
    integer, intent(in)                :: npx, npy, refine           !< Number of x,y points and nest refinement
    character(len=*), intent(in)       :: surface_dir, file_prefix   !< Surface directory and file tag
    character(len=*), intent(in)       :: var_name                   !< Variable name in netCDF file
    real*4, allocatable, intent(out)   :: data_grid(:,:)             !< Output data grid
    integer, intent(in)                :: parent_tile                !< Parent tile number
    integer, intent(in), optional      :: time                       !< Optional month number for time-varying parameters

    character(len=256) :: res_str, parent_str
    character(len=16)  :: halo
    character(len=512) :: nc_filename
    integer :: nx_cubic, nx, ny, fp_nx, fp_ny
    integer :: fp_istart_fine, fp_iend_fine, fp_jstart_fine, fp_jend_fine
    integer :: this_pe

    this_pe = mpp_pe()

    nx_cubic = npx - 1
    nx = npx - 1
    ny = npy - 1

    fp_istart_fine = 0
    fp_iend_fine = nx * refine
    fp_jstart_fine = 0
    fp_jend_fine = ny * refine

    fp_nx = fp_iend_fine - fp_istart_fine
    fp_ny = fp_jend_fine - fp_jstart_fine

    if (debug_log) print '("[INFO] WDR NCREAD LOFC mn_static_read_hires npe=",I0,I4,I4," ",A128," ",A128)', this_pe, fp_nx, fp_ny, var_name, nc_filename

    call mn_static_filename(surface_dir, parent_tile, file_prefix, refine, nc_filename)

    if (present(time)) then
      call alloc_read_data(nc_filename, var_name, fp_nx, fp_ny, data_grid, time)
    else
      call alloc_read_data(nc_filename, var_name, fp_nx, fp_ny, data_grid)
    endif

  end subroutine mn_static_read_hires_r4

  !>@brief The subroutine 'mn_static_read_hires_r8' loads high resolution data from netCDF
  !>@details Gathers a single variable from the netCDF file
  subroutine mn_static_read_hires_r8(npx, npy, refine, surface_dir, file_prefix, var_name, data_grid, parent_tile)
    integer, intent(in)                :: npx, npy, refine           !< Number of x,y points and nest refinement
    character(len=*), intent(in)       :: surface_dir, file_prefix   !< Surface directory and file tag
    character(len=*), intent(in)       :: var_name                   !< Variable name in netCDF file
    real*8, allocatable, intent(out)   :: data_grid(:,:)             !< Output data grid
    integer, intent(in)                :: parent_tile                !< Parent tile number

    character(len=256) :: res_str, parent_str
    character(len=16)  :: halo
    character(len=512) :: nc_filename

    integer :: nx_cubic, nx, ny, fp_nx, fp_ny
    integer :: fp_istart_fine, fp_iend_fine, fp_jstart_fine, fp_jend_fine
    integer :: this_pe

    this_pe = mpp_pe()

    nx_cubic = npx - 1
    nx = npx - 1
    ny = npy - 1

    fp_istart_fine = 0
    fp_iend_fine = nx * refine
    fp_jstart_fine = 0
    fp_jend_fine = ny * refine

    fp_nx = fp_iend_fine - fp_istart_fine
    fp_ny = fp_jend_fine - fp_jstart_fine

    if (debug_log) print '("[INFO] WDR NCREAD LOFC mn_static_read_hires npe=",I0,I4,I4," ",A128," ",A128)', this_pe, fp_nx, fp_ny, var_name, nc_filename

    call mn_static_filename(surface_dir, parent_tile, file_prefix, refine, nc_filename)

    call alloc_read_data(nc_filename, var_name, fp_nx, fp_ny, data_grid)

  end subroutine mn_static_read_hires_r8


  !!============================================================================
  !! Step 5.2 -- Recalculate nest halo weights
  !!============================================================================

  !>@brief The subroutine 'mn_meta_recalc' recalculates nest halo weights
  subroutine mn_meta_recalc( delta_i_c, delta_j_c, x_refine, y_refine, tile_geo, parent_geo, fp_super_tile_geo, &
      is_fine_pe, nest_domain, position, p_grid, n_grid, wt, istart_coarse, jstart_coarse)
    integer, intent(in)                           :: delta_i_c, delta_j_c                     !< Nest motion in delta i,j
    integer, intent(in)                           :: x_refine, y_refine                       !< Nest refinement
    type(grid_geometry), intent(inout)            :: tile_geo, parent_geo, fp_super_tile_geo  !< tile geometries
    logical, intent(in)                           :: is_fine_pe                               !< Is this a nest PE?
    type(nest_domain_type), intent(in)            :: nest_domain                              !< Nest domain structure
    real(kind=R_GRID), allocatable, intent(inout) :: p_grid(:,:,:)                            !< Parent lat/lon grid
    real(kind=R_GRID), allocatable, intent(inout) :: n_grid(:,:,:)                            !< Nest lat/lon grid
    real, allocatable, intent(inout)              :: wt(:,:,:)                                !< Interpolation weights
    integer, intent(inout)                        :: position                                 !< Stagger
    integer, intent(in)                           :: istart_coarse, jstart_coarse             !< Initian nest offsets

    type(bbox) :: wt_fine, wt_coarse
    integer    :: this_pe

    this_pe = mpp_pe()

    ! Update the coarse and fine indices after shifting the nest
    if (is_fine_pe) then

      if (debug_log) print '("[INFO] WDR NRD4 is_fine_pe=TRUE about to call bbox_get_C2F_index. npe=",I0, " position=",I0)', this_pe, position

      !!===========================================================
      !!
      !!  Recalculate halo weights
      !!
      !!===========================================================

      call bbox_get_C2F_index(nest_domain, wt_fine, wt_coarse, EAST,  position)
      call calc_nest_halo_weights(wt_fine, wt_coarse, p_grid, n_grid, wt, istart_coarse, jstart_coarse, x_refine, y_refine)

      call bbox_get_C2F_index(nest_domain, wt_fine, wt_coarse, WEST,  position)
      call calc_nest_halo_weights(wt_fine, wt_coarse, p_grid, n_grid, wt, istart_coarse, jstart_coarse, x_refine, y_refine)

      call bbox_get_C2F_index(nest_domain, wt_fine, wt_coarse, NORTH,  position)
      call calc_nest_halo_weights(wt_fine, wt_coarse, p_grid, n_grid, wt, istart_coarse, jstart_coarse, x_refine, y_refine)

      call bbox_get_C2F_index(nest_domain, wt_fine, wt_coarse, SOUTH,  position)
      call calc_nest_halo_weights(wt_fine, wt_coarse, p_grid, n_grid, wt, istart_coarse, jstart_coarse, x_refine, y_refine)

    endif

  end subroutine mn_meta_recalc


  !!============================================================================
  !! Step 5.3 -- Adjust index by delta_i_c, delta_j_c
  !!============================================================================

  !>@brief The subroutine 'mn_shift_index' adjusts the index array for a nest move
  !>@details Fast routine to increment indices by the delta in i,j direction
  subroutine mn_shift_index(delta_i_c, delta_j_c, ind)
    integer, intent(in)                    :: delta_i_c, delta_j_c    !< Nest move deltas in i,j
    integer, allocatable, intent(inout)    :: ind(:,:,:)              !< Nest to parent index

    ! Shift the index by the delta of this nest move.
    ! TODO -- validate that we are not moving off the edge of the parent grid.
    integer  :: i, j

    do i = lbound(ind,1), ubound(ind,1)
      do j = lbound(ind,2), ubound(ind,2)
        ind(i,j,1) = ind(i,j,1) + delta_i_c
        ind(i,j,2) = ind(i,j,2) + delta_j_c
      enddo
    enddo

  end subroutine mn_shift_index


  !================================================================================
  !
  !  Prognostic and Physics Variable Nest Motion
  !
  !================================================================================

  !!============================================================================
  !! Step 6   Shift the data on each nest PE
  !!            -- similar to med_nest_move in HWRF
  !!============================================================================

  !>@brief The subroutine 'mn_prog_shift_data' shifts the data on each nest PE
  !>@details Iterates through the prognostic variables
  subroutine mn_prog_shift_data(Atm, n, child_grid_num, wt_h, wt_u, wt_v, &
      delta_i_c, delta_j_c, x_refine, y_refine, &
      is_fine_pe, nest_domain, nz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)                                      !< Atm data array
    integer, intent(in)                                      :: n, child_grid_num                           !< Grid numbers
    real, allocatable, intent(in)                            :: wt_h(:,:,:), wt_u(:,:,:), wt_v(:,:,:)       !< Interpolation weights
    integer, intent(in)                                      :: delta_i_c, delta_j_c, x_refine, y_refine    !< Delta i,j, nest refinement
    logical, intent(in)                                      :: is_fine_pe                                  !< Is this is a nest PE?
    type(nest_domain_type), intent(inout)                    :: nest_domain                                 !< Nest domain structure
    integer, intent(in)                                      :: nz                                          !< Number of vertical levels

    ! Constants for mpp calls
    integer  :: interp_type   = 1    ! cell-centered A-grid
    integer  :: interp_type_u = 4    ! D-grid
    integer  :: interp_type_v = 4    ! D-grid
    integer  :: position      = CENTER ! CENTER, NORTH, EAST
    integer  :: position_u    = NORTH
    integer  :: position_v    = EAST

    type(fv_moving_nest_prog_type), pointer :: mn_prog

    mn_prog => Moving_nest(n)%mn_prog

    call mn_var_shift_data(Atm(n)%q_con, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
        delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)

    call mn_var_shift_data(Atm(n)%pt, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
        delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)

    call mn_var_shift_data(Atm(n)%w, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
        delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)

    !call mn_var_shift_data(Atm(n)%omga, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
    !     delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)

    call mn_var_shift_data(Atm(n)%delp, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
        delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)

    !call mn_var_shift_data(Atm(n)%delz, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
    !     delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)

    call mn_var_shift_data(mn_prog%delz, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
        delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)

    call mn_var_shift_data(Atm(n)%ua, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
        delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)

    call mn_var_shift_data(Atm(n)%va, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
        delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)

    call mn_var_shift_data(Atm(n)%q, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
        delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)

    !if (debug_log) print '("[INFO] WDR MV_NST6 show wt_u run step 6 atmosphere.F90 npe=",I0," n=",I0)', this_pe, n
    !call check_array(Atm(n)%neststruct%wt_u, this_pe, "Atm(n)%neststruct%wt_u", 0.0, 1.0)
    !call check_array(wt_u, this_pe, "wt_u", 0.0, 1.0)
    !if (debug_log) print '("[INFO] WDR MV_NST6 stagger run step 6 atmosphere.F90 npe=",I0," n=",I0)', this_pe, n

    call mn_var_shift_data(Atm(n)%u, interp_type_u, wt_u, Atm(child_grid_num)%neststruct%ind_u, &
        delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position_u, nz)

    call mn_var_shift_data(Atm(n)%v, interp_type_v, wt_v, Atm(child_grid_num)%neststruct%ind_v, &
        delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position_v, nz)

  end subroutine mn_prog_shift_data


  !!============================================================================
  !! Step 6 - per variable
  !!============================================================================

  !>@brief The subroutine 'mn_prog_shift_data_r4_2d' shifts the data for a variable on each nest PE
  !>@details For single variable
  subroutine mn_var_shift_data_r4_2d(data_var, interp_type, wt, ind, delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
    real*4, allocatable, intent(inout)          :: data_var(:,:)                                !< Data variable
    integer, intent(in)                         :: interp_type                                  !< Interpolation stagger type
    real, allocatable, intent(in)               :: wt(:,:,:)                                    !< Interpolation weight array
    integer, allocatable, intent(in)            :: ind(:,:,:)                                   !< Fine to coarse index array
    integer, intent(in)                         :: delta_i_c, delta_j_c, x_refine, y_refine     !< delta i,j for nest move.  Nest refinement.
    logical, intent(in)                         :: is_fine_pe                                   !< Is nest PE?
    type(nest_domain_type), intent(inout)       :: nest_domain                                  !< Nest domain structure
    integer, intent(in)                         :: position                                     !< Grid offset

    real*4, dimension(:,:), allocatable :: nbuffer, sbuffer, ebuffer, wbuffer
    logical         :: parent_proc, child_proc
    type(bbox)      :: north_fine, north_coarse ! step 4
    type(bbox)      :: south_fine, south_coarse
    type(bbox)      :: east_fine, east_coarse
    type(bbox)      :: west_fine, west_coarse
    integer         :: my_stat
    character(256)  :: my_errmsg
    integer         :: is, ie, js, je
    integer         :: this_pe

    integer         :: nest_level = 1  ! WDR TODO allow to vary

    this_pe = mpp_pe()

    !!===========================================================
    !!
    !! Fill halo buffers
    !!
    !!===========================================================

    if (debug_log) print '("[INFO] WDR NRD5. npe=",I0)', this_pe

    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%tile_fine=",I0," %tile_coarse=",I0)', this_pe, nest_domain%tile_fine, nest_domain%tile_coarse

    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%istart_fine=",I0," %iend_fine=",I0)', this_pe, nest_domain%istart_fine,  nest_domain%iend_fine
    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%jstart_fine=",I0," %jend_fine=",I0)', this_pe, nest_domain%jstart_fine,  nest_domain%jend_fine

    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%istart_coarse=",I0," %iend_coarse=",I0)', this_pe, nest_domain%istart_coarse,  nest_domain%iend_coarse
    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%jstart_coarse=",I0," %jend_coarse=",I0)', this_pe, nest_domain%jstart_coarse,  nest_domain%jend_coarse

    if (debug_log) print '("[INFO] data_var npe=",I0," data_var(",I0,"-",I0,",",I0,"-",I0,")")', this_pe,  lbound(data_var, 1), ubound(data_var, 1), lbound(data_var, 2), ubound(data_var, 2)

    if (debug_log) print '("[INFO] wt npe=",I0," wt(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', this_pe,  lbound(wt, 1), ubound(wt, 1), lbound(wt, 2), ubound(wt, 2), lbound(wt, 3), ubound(wt, 3)

    !====================================================
    if (debug_log) print '("[INFO] WDR ALL1. npe=",I0," position=",I0)', this_pe, position
    call alloc_halo_buffer(nbuffer, north_fine, north_coarse, nest_domain, NORTH,  position)
    call alloc_halo_buffer(sbuffer, south_fine, south_coarse, nest_domain, SOUTH,  position)
    call alloc_halo_buffer(ebuffer, east_fine,  east_coarse,  nest_domain, EAST,   position)
    call alloc_halo_buffer(wbuffer, west_fine,  west_coarse,  nest_domain, WEST,   position)

    if (debug_log) print '("[INFO] WDR allocate_halo_buffers DONE. npe=",I0)', this_pe

    !====================================================

    if (debug_log) print '("[INFO] WDR NRF0.d mn_var_shift_data npe=",I0," data_var(",I0,",",I0,")")', this_pe, size(data_var,1), size(data_var,2)
    if (debug_log) print '("[INFO] WDR NRF1 mn_var_shift_data start. npe=",I0)', this_pe

    ! Passes data from coarse grid to fine grid's halo buffers; requires nest_domain to be intent(inout)
    call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

    if (debug_log) print '("[INFO] WDR NRF2 mn_var_shift_data start. npe=",I0)', this_pe

    if (is_fine_pe) then
      if (debug_log) print '("[INFO] WDR NRF3 mn_var_shift_data start. npe=",I0)', this_pe

      !!===========================================================
      !!
      !! Shift grids internal to each nest PE
      !!
      !!===========================================================

      if ( delta_i_c .ne. 0 ) then
        if (debug_log) print '("[INFO] WDR NREX mn_var_shift_data start. npe=",I0)', this_pe
        data_var = eoshift(data_var, x_refine * delta_i_c, DIM=1)
      endif

      if (delta_j_c .ne.  0) then
        if (debug_log) print '("[INFO] WDR NREY mn_var_shift_data start. npe=",I0)', this_pe
        data_var = eoshift(data_var, y_refine * delta_j_c, DIM=2)
      endif

      !!===========================================================
      !!
      !! Apply halo data
      !!
      !!===========================================================

      if (debug_log) print '("[INFO] WDR NRFI mn_var_shift_data start. npe=",I0)', this_pe

      call fill_nest_from_buffer(interp_type, data_var, nbuffer, north_fine, north_coarse, NORTH, x_refine, y_refine, wt, ind)
      if (debug_log) print '("[INFO] WDR NRF N mn_var_shift_data start. npe=",I0)', this_pe

      call fill_nest_from_buffer(interp_type, data_var, sbuffer, south_fine, south_coarse, SOUTH, x_refine, y_refine, wt, ind)
      if (debug_log) print '("[INFO] WDR NRF S mn_var_shift_data start. npe=",I0)', this_pe

      call fill_nest_from_buffer(interp_type, data_var, ebuffer, east_fine, east_coarse, EAST, x_refine, y_refine, wt, ind)
      if (debug_log) print '("[INFO] WDR NRF E mn_var_shift_data start. npe=",I0)', this_pe

      call fill_nest_from_buffer(interp_type, data_var, wbuffer, west_fine, west_coarse, WEST, x_refine, y_refine, wt, ind)
      if (debug_log) print '("[INFO] WDR NRF W mn_var_shift_data start. npe=",I0)', this_pe

    endif

    deallocate(nbuffer)
    deallocate(sbuffer)
    deallocate(ebuffer)
    deallocate(wbuffer)

  end subroutine mn_var_shift_data_r4_2d

  !>@brief The subroutine 'mn_prog_shift_data_r8_2d' shifts the data for a variable on each nest PE
  !>@details For one double precision 2D variable
  subroutine mn_var_shift_data_r8_2d(data_var, interp_type, wt, ind, delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)

    real*8, allocatable, intent(inout)          :: data_var(:,:)                            !< Data variable
    integer, intent(in)                         :: interp_type                              !< Interpolation stagger type
    real, allocatable, intent(in)               :: wt(:,:,:)                                !< Interpolation weight array
    integer, allocatable, intent(in)            :: ind(:,:,:)                               !< Fine to coarse index array
    integer, intent(in)                         :: delta_i_c, delta_j_c, x_refine, y_refine !< delta i,j for nest move.  Nest refinement.
    logical, intent(in)                         :: is_fine_pe                                   !< Is nest PE?
    type(nest_domain_type), intent(inout)       :: nest_domain                                  !< Nest domain structure
    integer, intent(in)                         :: position                                     !< Grid offset

    real*8, dimension(:,:), allocatable :: nbuffer, sbuffer, ebuffer, wbuffer
    logical         :: parent_proc, child_proc
    type(bbox)      :: north_fine, north_coarse ! step 4
    type(bbox)      :: south_fine, south_coarse
    type(bbox)      :: east_fine, east_coarse
    type(bbox)      :: west_fine, west_coarse
    integer         :: my_stat
    character(256)  :: my_errmsg
    integer         :: is, ie, js, je
    integer         :: this_pe

    integer         :: nest_level = 1  ! WDR TODO allow to vary

    this_pe = mpp_pe()

    !!===========================================================
    !!
    !! Fill halo buffers
    !!
    !!===========================================================

    if (debug_log) print '("[INFO] WDR NRD5. npe=",I0)', this_pe


    call alloc_halo_buffer(nbuffer, north_fine, north_coarse, nest_domain, NORTH,  position)
    call alloc_halo_buffer(sbuffer, south_fine, south_coarse, nest_domain, SOUTH,  position)
    call alloc_halo_buffer(ebuffer, east_fine,  east_coarse,  nest_domain, EAST,   position)
    call alloc_halo_buffer(wbuffer, west_fine,  west_coarse,  nest_domain, WEST,   position)

    ! Passes data from coarse grid to fine grid's halo buffers; requires nest_domain to be intent(inout)
    call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

    if (is_fine_pe) then

      !!===========================================================
      !!
      !! Shift grids internal to each nest PE
      !!
      !!===========================================================

      if ( delta_i_c .ne. 0 ) then
        data_var = eoshift(data_var, x_refine * delta_i_c, DIM=1)
      endif

      if (delta_j_c .ne.  0) then
        data_var = eoshift(data_var, y_refine * delta_j_c, DIM=2)
      endif

      !!===========================================================
      !!
      !! Apply halo data
      !!
      !!===========================================================

      call fill_nest_from_buffer(interp_type, data_var, nbuffer, north_fine, north_coarse, NORTH, x_refine, y_refine, wt, ind)
      call fill_nest_from_buffer(interp_type, data_var, sbuffer, south_fine, south_coarse, SOUTH, x_refine, y_refine, wt, ind)
      call fill_nest_from_buffer(interp_type, data_var, ebuffer, east_fine, east_coarse, EAST, x_refine, y_refine, wt, ind)
      call fill_nest_from_buffer(interp_type, data_var, wbuffer, west_fine, west_coarse, WEST, x_refine, y_refine, wt, ind)
    endif

    deallocate(nbuffer)
    deallocate(sbuffer)
    deallocate(ebuffer)
    deallocate(wbuffer)

  end subroutine mn_var_shift_data_r8_2d

  !>@brief The subroutine 'mn_prog_shift_data_r4_3d' shifts the data for a variable on each nest PE
  !>@details For one single precision 3D variable
  subroutine mn_var_shift_data_r4_3d(data_var, interp_type, wt, ind, delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)

    real*4, allocatable, intent(inout)          :: data_var(:,:,:)                          !< Data variable
    integer, intent(in)                         :: interp_type                              !< Interpolation stagger type
    real, allocatable, intent(in)               :: wt(:,:,:)                                !< Interpolation weight array
    integer, allocatable, intent(in)            :: ind(:,:,:)                               !< Fine to coarse index array
    integer, intent(in)                         :: delta_i_c, delta_j_c, x_refine, y_refine !< delta i,j for nest move.  Nest refinement.
    logical, intent(in)                         :: is_fine_pe                               !< Is nest PE?
    type(nest_domain_type), intent(inout)       :: nest_domain                              !< Nest domain structure
    integer, intent(in)                         :: position, nz                             !< Grid offset, number of vertical levels

    real*4, dimension(:,:,:), allocatable :: nbuffer, sbuffer, ebuffer, wbuffer
    logical         :: parent_proc, child_proc
    type(bbox)      :: north_fine, north_coarse ! step 4
    type(bbox)      :: south_fine, south_coarse
    type(bbox)      :: east_fine, east_coarse
    type(bbox)      :: west_fine, west_coarse
    integer         :: my_stat
    character(256)  :: my_errmsg
    integer         :: is, ie, js, je
    integer         :: this_pe

    integer         :: nest_level = 1  ! WDR TODO allow to vary

    this_pe = mpp_pe()

    !!===========================================================
    !!
    !! Fill halo buffers
    !!
    !!===========================================================

    call alloc_halo_buffer(nbuffer, north_fine, north_coarse, nest_domain, NORTH,  position, nz)
    call alloc_halo_buffer(sbuffer, south_fine, south_coarse, nest_domain, SOUTH,  position, nz)
    call alloc_halo_buffer(ebuffer, east_fine,  east_coarse,  nest_domain, EAST,   position, nz)
    call alloc_halo_buffer(wbuffer, west_fine,  west_coarse,  nest_domain, WEST,   position, nz)


    !====================================================
    ! Passes data from coarse grid to fine grid's halo buffers; requires nest_domain to be intent(inout)
    call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

    if (is_fine_pe) then

      !!===========================================================
      !!
      !! Shift grids internal to each nest PE
      !!
      !!===========================================================

      if ( delta_i_c .ne. 0 ) then
        data_var = eoshift(data_var, x_refine * delta_i_c, DIM=1)
      endif

      if (delta_j_c .ne.  0) then
        data_var = eoshift(data_var, y_refine * delta_j_c, DIM=2)
      endif

      !!===========================================================
      !!
      !! Apply halo data
      !!
      !!===========================================================

      call fill_nest_from_buffer(interp_type, data_var, nbuffer, north_fine, north_coarse, nz, NORTH, x_refine, y_refine, wt, ind)
      call fill_nest_from_buffer(interp_type, data_var, sbuffer, south_fine, south_coarse, nz, SOUTH, x_refine, y_refine, wt, ind)
      call fill_nest_from_buffer(interp_type, data_var, ebuffer, east_fine, east_coarse, nz, EAST, x_refine, y_refine, wt, ind)
      call fill_nest_from_buffer(interp_type, data_var, wbuffer, west_fine, west_coarse, nz, WEST, x_refine, y_refine, wt, ind)
    endif

    deallocate(nbuffer)
    deallocate(sbuffer)
    deallocate(ebuffer)
    deallocate(wbuffer)

  end subroutine mn_var_shift_data_r4_3d


  !>@brief The subroutine 'mn_prog_shift_data_r8_3d' shifts the data for a variable on each nest PE
  !>@details For one double precision 3D variable
  subroutine mn_var_shift_data_r8_3d(data_var, interp_type, wt, ind, delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)

    real*8, allocatable, intent(inout)          :: data_var(:,:,:)                              !< Data variable
    integer, intent(in)                         :: interp_type                                  !< Interpolation stagger type
    real, allocatable, intent(in)               :: wt(:,:,:)                                    !< Interpolation weight array
    integer, allocatable, intent(in)            :: ind(:,:,:)                                   !< Fine to coarse index array
    integer, intent(in)                         :: delta_i_c, delta_j_c, x_refine, y_refine     !< delta i,j for nest move.  Nest refinement.
    logical, intent(in)                         :: is_fine_pe                                   !< Is nest PE?
    type(nest_domain_type), intent(inout)       :: nest_domain                                  !< Nest domain structure
    integer, intent(in)                         :: position, nz                                 !< Grid offset, number vertical levels

    real*8, dimension(:,:,:), allocatable :: nbuffer, sbuffer, ebuffer, wbuffer
    logical         :: parent_proc, child_proc
    type(bbox)      :: north_fine, north_coarse ! step 4
    type(bbox)      :: south_fine, south_coarse
    type(bbox)      :: east_fine, east_coarse
    type(bbox)      :: west_fine, west_coarse
    integer         :: my_stat
    character(256)  :: my_errmsg
    integer         :: is, ie, js, je
    integer         :: this_pe

    integer         :: nest_level = 1  ! WDR TODO allow to vary

    this_pe = mpp_pe()

    !!===========================================================
    !!
    !! Fill halo buffers
    !!
    !!===========================================================

    call alloc_halo_buffer(nbuffer, north_fine, north_coarse, nest_domain, NORTH,  position, nz)
    call alloc_halo_buffer(sbuffer, south_fine, south_coarse, nest_domain, SOUTH,  position, nz)
    call alloc_halo_buffer(ebuffer, east_fine,  east_coarse,  nest_domain, EAST,   position, nz)
    call alloc_halo_buffer(wbuffer, west_fine,  west_coarse,  nest_domain, WEST,   position, nz)

    !====================================================
    ! Passes data from coarse grid to fine grid's halo buffers; requires nest_domain to be intent(inout)
    call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

    if (is_fine_pe) then
      !!===========================================================
      !!
      !! Shift grids internal to each nest PE
      !!
      !!===========================================================

      if ( delta_i_c .ne. 0 ) then
        data_var = eoshift(data_var, x_refine * delta_i_c, DIM=1)
      endif

      if (delta_j_c .ne.  0) then
        data_var = eoshift(data_var, y_refine * delta_j_c, DIM=2)
      endif

      !!===========================================================
      !!
      !! Apply halo data
      !!
      !!===========================================================

      call fill_nest_from_buffer(interp_type, data_var, nbuffer, north_fine, north_coarse, nz, NORTH, x_refine, y_refine, wt, ind)
      call fill_nest_from_buffer(interp_type, data_var, sbuffer, south_fine, south_coarse, nz, SOUTH, x_refine, y_refine, wt, ind)
      call fill_nest_from_buffer(interp_type, data_var, ebuffer, east_fine, east_coarse, nz, EAST, x_refine, y_refine, wt, ind)
      call fill_nest_from_buffer(interp_type, data_var, wbuffer, west_fine, west_coarse, nz, WEST, x_refine, y_refine, wt, ind)
    endif

    deallocate(nbuffer)
    deallocate(sbuffer)
    deallocate(ebuffer)
    deallocate(wbuffer)

  end subroutine mn_var_shift_data_r8_3d


  !>@brief The subroutine 'mn_prog_shift_data_r4_4d' shifts the data for a variable on each nest PE
  !>@details For one single precision 4D variable
  subroutine mn_var_shift_data_r4_4d(data_var, interp_type, wt, ind, delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)
    real*4, allocatable, intent(inout)          :: data_var(:,:,:,:)                            !< Data variable
    integer, intent(in)                         :: interp_type                                  !< Interpolation stagger type
    real, allocatable, intent(in)               :: wt(:,:,:)                                    !< Interpolation weight array
    integer, allocatable, intent(in)            :: ind(:,:,:)                                   !< Fine to coarse index array
    integer, intent(in)                         :: delta_i_c, delta_j_c, x_refine, y_refine     !< delta i,j for nest move.  Nest refinement.
    logical, intent(in)                         :: is_fine_pe                                   !< Is nest PE?
    type(nest_domain_type), intent(inout)       :: nest_domain                                  !< Nest domain structure
    integer, intent(in)                         :: position, nz                                 !< Grid offset, number of vertical levels

    real*4, dimension(:,:,:,:), allocatable     :: nbuffer, sbuffer, ebuffer, wbuffer
    logical         :: parent_proc, child_proc
    type(bbox)      :: north_fine, north_coarse ! step 4
    type(bbox)      :: south_fine, south_coarse
    type(bbox)      :: east_fine, east_coarse
    type(bbox)      :: west_fine, west_coarse
    integer         :: my_stat
    character(256)  :: my_errmsg
    integer         :: n4d
    integer         :: this_pe
    integer         :: is, ie, js, je
    integer         :: nest_level = 1  ! WDR TODO allow to vary

    this_pe = mpp_pe()

    n4d = ubound(data_var, 4)

    !!===========================================================
    !!
    !! Fill halo buffers
    !!
    !!===========================================================

    call alloc_halo_buffer(nbuffer, north_fine, north_coarse, nest_domain, NORTH,  position, nz, n4d)
    call alloc_halo_buffer(sbuffer, south_fine, south_coarse, nest_domain, SOUTH,  position, nz, n4d)
    call alloc_halo_buffer(ebuffer, east_fine,  east_coarse,  nest_domain, EAST,   position, nz, n4d)
    call alloc_halo_buffer(wbuffer, west_fine,  west_coarse,  nest_domain, WEST,   position, nz, n4d)

    !====================================================

    ! Passes data from coarse grid to fine grid's halo
    call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

    if (is_fine_pe) then
      !!===========================================================
      !!
      !! Shift grids internal to each nest PE
      !!
      !!===========================================================

      if ( delta_i_c .ne. 0 ) then
        data_var = eoshift(data_var, x_refine * delta_i_c, DIM=1)
      endif

      if (delta_j_c .ne.  0) then
        data_var = eoshift(data_var, y_refine * delta_j_c, DIM=2)
      endif

      !!===========================================================
      !!
      !! Apply halo data
      !!
      !!===========================================================

      call fill_nest_from_buffer(interp_type, data_var, nbuffer, north_fine, north_coarse, nz, NORTH, x_refine, y_refine, wt, ind)
      call fill_nest_from_buffer(interp_type, data_var, sbuffer, south_fine, south_coarse, nz, SOUTH, x_refine, y_refine, wt, ind)
      call fill_nest_from_buffer(interp_type, data_var, ebuffer, east_fine, east_coarse, nz, EAST, x_refine, y_refine, wt, ind)
      call fill_nest_from_buffer(interp_type, data_var, wbuffer, west_fine, west_coarse, nz, WEST, x_refine, y_refine, wt, ind)
    endif

    deallocate(nbuffer)
    deallocate(sbuffer)
    deallocate(ebuffer)
    deallocate(wbuffer)

  end subroutine mn_var_shift_data_r4_4d


  !>@brief The subroutine 'mn_prog_shift_data_r8_4d' shifts the data for a variable on each nest PE
  !>@details For one double precision 4D variable
  subroutine mn_var_shift_data_r8_4d(data_var, interp_type, wt, ind, delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)
    real*8, allocatable, intent(inout)          :: data_var(:,:,:,:)                            !< Data variable
    integer, intent(in)                         :: interp_type                                  !< Interpolation stagger type
    real, allocatable, intent(in)               :: wt(:,:,:)                                    !< Interpolation weight array
    integer, allocatable, intent(in)            :: ind(:,:,:)                                   !< Fine to coarse index array
    integer, intent(in)                         :: delta_i_c, delta_j_c, x_refine, y_refine     !< delta i,j for nest move.  Nest refinement.
    logical, intent(in)                         :: is_fine_pe                                   !< Is nest PE?
    type(nest_domain_type), intent(inout)       :: nest_domain                                  !< Nest domain structure
    integer, intent(in)                         :: position, nz                                 !< Grid offset, number of vertical levels

    real*8, dimension(:,:,:,:), allocatable     :: nbuffer, sbuffer, ebuffer, wbuffer
    logical         :: parent_proc, child_proc
    type(bbox)      :: north_fine, north_coarse ! step 4
    type(bbox)      :: south_fine, south_coarse
    type(bbox)      :: east_fine, east_coarse
    type(bbox)      :: west_fine, west_coarse
    integer         :: my_stat
    character(256)  :: my_errmsg
    integer         :: n4d
    integer         :: this_pe
    integer         :: is, ie, js, je
    integer         :: nest_level = 1  ! WDR TODO allow to vary

    this_pe = mpp_pe()

    n4d = ubound(data_var, 4)

    !!===========================================================
    !!
    !! Fill halo buffers
    !!
    !!===========================================================

    call alloc_halo_buffer(nbuffer, north_fine, north_coarse, nest_domain, NORTH,  position, nz, n4d)
    call alloc_halo_buffer(sbuffer, south_fine, south_coarse, nest_domain, SOUTH,  position, nz, n4d)
    call alloc_halo_buffer(ebuffer, east_fine,  east_coarse,  nest_domain, EAST,   position, nz, n4d)
    call alloc_halo_buffer(wbuffer, west_fine,  west_coarse,  nest_domain, WEST,   position, nz, n4d)

    !====================================================
    ! Passes data from coarse grid to fine grid's halo
    call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

    if (is_fine_pe) then
      !!===========================================================
      !!
      !! Shift grids internal to each nest PE
      !!
      !!===========================================================

      if ( delta_i_c .ne. 0 ) then
        data_var = eoshift(data_var, x_refine * delta_i_c, DIM=1)
      endif

      if (delta_j_c .ne.  0) then
        data_var = eoshift(data_var, y_refine * delta_j_c, DIM=2)
      endif

      !!===========================================================
      !!
      !! Apply halo data
      !!
      !!===========================================================

      call fill_nest_from_buffer(interp_type, data_var, nbuffer, north_fine, north_coarse, nz, NORTH, x_refine, y_refine, wt, ind)
      call fill_nest_from_buffer(interp_type, data_var, sbuffer, south_fine, south_coarse, nz, SOUTH, x_refine, y_refine, wt, ind)
      call fill_nest_from_buffer(interp_type, data_var, ebuffer, east_fine, east_coarse, nz, EAST, x_refine, y_refine, wt, ind)
      call fill_nest_from_buffer(interp_type, data_var, wbuffer, west_fine, west_coarse, nz, WEST, x_refine, y_refine, wt, ind)
    endif

    deallocate(nbuffer)
    deallocate(sbuffer)
    deallocate(ebuffer)
    deallocate(wbuffer)

  end subroutine mn_var_shift_data_r8_4d


  !================================================================================
  !
  !  Step 7 -- Gridstruct resetting and reallocation of static buffers
  !      init_grid() also updates the wt arrays
  !================================================================================

  !>@brief The subroutine 'mn_meta_reset_gridstruct' resets navigation data and reallocates needed data in the gridstruct after nest move
  !>@details This routine is computationally demanding and is a target for later optimization.
  subroutine mn_meta_reset_gridstruct(Atm, n, child_grid_num, nest_domain, fp_super_tile_geo, x_refine, y_refine, is_fine_pe, wt_h, wt_u, wt_v, a_step, dt_atmos)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)                               !< Atm data array
    integer, intent(in)                              :: n, child_grid_num                            !< This level and nest level
    type(nest_domain_type),     intent(in)           :: nest_domain                                  !< Nest domain structure
    type(grid_geometry), intent(in)                  :: fp_super_tile_geo                            !< Parent high-resolution geometry
    integer, intent(in)                              :: x_refine, y_refine                           !< Nest refinement
    logical, intent(in)                              :: is_fine_pe                                   !< Is nest PE?
    real, allocatable, intent(in)                    :: wt_h(:,:,:), wt_u(:,:,:), wt_v(:,:,:)        !< Interpolation weights
    integer, intent(in)                              :: a_step                                       !< Which timestep
    real, intent(in)                                 :: dt_atmos                                     !< Timestep duration in seconds

    integer :: isg, ieg, jsg, jeg
    integer :: ng, pp, nn, parent_tile, refinement, ioffset, joffset
    integer :: this_pe, gid
    integer :: tile_coarse(2)
    integer :: half_x, half_y

    real(kind=R_GRID)   :: pi = 4 * atan(1.0d0)
    real                :: rad2deg, half_lat, half_lon

    ! Coriolis parameter variables
    real                :: alpha = 0.
    real, pointer, dimension(:,:,:) :: grid, agrid
    real, pointer, dimension(:,:) :: fC, f0
    integer             :: isd, ied, jsd, jed
    integer             :: i, j

    logical, save       :: first_time = .true.
    integer, save       :: id_reset1, id_reset2, id_reset3, id_reset4, id_reset5, id_reset6, id_reset7

    logical             :: use_timers = .false. !  Set this to true to generate performance profiling information in out.* file

    if (first_time .and. use_timers) then
      id_reset1     = mpp_clock_id ('MN 7 Reset 1',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
      id_reset2     = mpp_clock_id ('MN 7 Reset 2',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
      id_reset3     = mpp_clock_id ('MN 7 Reset 3',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
      id_reset4     = mpp_clock_id ('MN 7 Reset 4',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
      id_reset5     = mpp_clock_id ('MN 7 Reset 5',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
      id_reset6     = mpp_clock_id ('MN 7 Reset 6',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
      id_reset7     = mpp_clock_id ('MN 7 Reset 7',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
    endif

    rad2deg = 180.0 / pi

    this_pe = mpp_pe()
    gid = this_pe

    parent_tile = Atm(child_grid_num)%neststruct%parent_tile
    ioffset = Atm(child_grid_num)%neststruct%ioffset
    joffset = Atm(child_grid_num)%neststruct%joffset

    ! Log the bounds of this PE's grid after nest motion.  TODO replace step 4 with timestep
    if (is_fine_pe .and. debug_log) then
      call show_nest_grid(Atm(n), this_pe, 4)
    endif

    !  Reset the gridstruct values for the nest
    if (is_fine_pe) then
      ! Fill in values from high resolution, full panel, supergrid
      if (use_timers) call mpp_clock_begin (id_reset1)

      call fill_grid_from_supergrid(Atm(n)%gridstruct%grid, CORNER, fp_super_tile_geo, ioffset, joffset, &
          x_refine, y_refine)
      call fill_grid_from_supergrid(Atm(n)%gridstruct%agrid, CENTER, fp_super_tile_geo, ioffset, joffset, &
          x_refine, y_refine)
      call fill_grid_from_supergrid(Atm(n)%gridstruct%grid_64, CORNER, fp_super_tile_geo, &
          ioffset, joffset, x_refine, y_refine)
      call fill_grid_from_supergrid(Atm(n)%gridstruct%agrid_64, CENTER, fp_super_tile_geo, &
          ioffset, joffset, x_refine, y_refine)

      ! What's the status of Atm(n)%grid_global?
      if (debug_log) print '("[INFO] WDR Atm(1) GLOBAL npe=",I0," grid_global(",I0,"-",I0",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', this_pe, &
          lbound(Atm(1)%grid_global,1), ubound(Atm(1)%grid_global,1), &
          lbound(Atm(1)%grid_global,2), ubound(Atm(1)%grid_global,2), &
          lbound(Atm(1)%grid_global,3), ubound(Atm(1)%grid_global,3), &
          lbound(Atm(1)%grid_global,4), ubound(Atm(1)%grid_global,4)

      if (debug_log) print '("[INFO] WDR Atm(n) GLOBAL npe=",I0," grid_global(",I0,"-",I0",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', this_pe, &
          lbound(Atm(n)%grid_global,1), ubound(Atm(n)%grid_global,1), &
          lbound(Atm(n)%grid_global,2), ubound(Atm(n)%grid_global,2), &
          lbound(Atm(n)%grid_global,3), ubound(Atm(n)%grid_global,3), &
          lbound(Atm(n)%grid_global,4), ubound(Atm(n)%grid_global,4)


      ! Reset the coriolis parameters, using code from external_ic.F90::get_external_ic()

      isd = Atm(n)%bd%isd
      ied = Atm(n)%bd%ied
      jsd = Atm(n)%bd%jsd
      jed = Atm(n)%bd%jed

      grid  => Atm(n)%gridstruct%grid
      agrid => Atm(n)%gridstruct%agrid
      fC    => Atm(n)%gridstruct%fC
      f0    => Atm(n)%gridstruct%f0

      ! * Initialize coriolis param:

      do j=jsd,jed+1
        do i=isd,ied+1
          fC(i,j) = 2.*omega*( -1.*cos(grid(i,j,1))*cos(grid(i,j,2))*sin(alpha) + &
              sin(grid(i,j,2))*cos(alpha) )
        enddo
      enddo

      do j=jsd,jed
        do i=isd,ied
          f0(i,j) = 2.*omega*( -1.*cos(agrid(i,j,1))*cos(agrid(i,j,2))*sin(alpha) + &
              sin(agrid(i,j,2))*cos(alpha) )
        enddo
      enddo





      !! Let this get reset in init_grid()/setup_aligned_nest()
      !call fill_grid_from_supergrid(Atm(n)%grid_global, CORNER, fp_super_tile_geo, &
      !     ioffset, joffset, x_refine, y_refine)

      if (use_timers) call mpp_clock_end (id_reset1)
      if (use_timers) call mpp_clock_begin (id_reset2)

      ! TODO should these get reset by init_grid instead??
      call fill_weight_grid(Atm(n)%neststruct%wt_h, wt_h)
      call fill_weight_grid(Atm(n)%neststruct%wt_u, wt_u)
      call fill_weight_grid(Atm(n)%neststruct%wt_v, wt_v)
      ! WDR TODO -- Seems like this is not used anywhere, other than being allocated, filled, deallocated
      !call fill_weight_grid(Atm(n)%neststruct%wt_b, wt_b)

      if (use_timers) call mpp_clock_end (id_reset2)

    endif

    if (debug_log) print '("[INFO] WDR INIT_GRID AP1 fv_moving_nest.F90 npe=",I0," n=",I0)', this_pe, n

    if (use_timers) call mpp_clock_begin (id_reset3)

    ! TODO Write clearer comments on what is happening here.

    ! This code runs several communications steps:
    !  1.  As npe=0, it gets the global_grid domain setup
    !  2.  sends the global_grid to the other parent PEs
    !  3.  global_grid is received in call to setup_aligned_nest() in fv_grid_tools.F90::init_grid()
    !  Other communication is contained full within setup_aligned_nest().

    ! Sends around data from the parent grids, and recomputes the update indices
    ! This code copied from fv_control.F90
    ! Need to SEND grid_global to any child grids; this is received in setup_aligned_nest in fv_grid_tools
    ! if (Atm(pp)%neststruct%nested) then

    ! TODO phrase this more carefully to choose the parent master PE grid if we are operating in a nested setup.
    ! Unlike in fv_control.F90, this will be running on Atm(1) when it's on pe=0, so we don't need to navigate to parent_grid.

    first_time = .false.

    ! Seems like we do not need to resend this -- setup_aligned_nest now saves the parent tile information during model initialization,
    !  which happens before we enter the moving nest code.
    if (this_pe .eq. 0 .and. first_time) then

      ! This is the Atm index for the nest values.
      pp = child_grid_num

      if (debug_log) print '("[INFO] WDR INIT_GRID AP2 fv_moving_nest.F90 npe=",I0," n=",I0," pp=",I0)', this_pe, n, pp

      refinement = x_refine
      ng = Atm(n)%ng

      call mpp_get_global_domain( Atm(n)%domain, isg, ieg, jsg, jeg)

      !if (debug_log) print '("[INFO] WDR INIT_GRID AP3.1 fv_moving_nest.F90 npe=",I0," gid=",I0," associated(parent_grid)=",L1)', this_pe, gid, associated(Atm(pp)%parent_grid)
      if (debug_log) print '("[INFO] WDR INIT_GRID AP3.1 fv_moving_nest.F90 npe=",I0," gid=",I0," parent_tile=",I0)', this_pe, gid, parent_tile
      if (debug_log) print '("[INFO] WDR INIT_GRID AP3.2 fv_moving_nest.F90 npe=",I0," gid=",I0," size(pelist)=",I0)', this_pe, gid, size(Atm(pp)%pelist)
      if (debug_log) print '("[INFO] WDR INIT_GRID AP3.3 fv_moving_nest.F90 npe=",I0," gid=",I0," pelist1=",I0)', this_pe, gid, Atm(pp)%pelist(1)
      !FIXME: Should replace this by generating the global grid (or at least one face thereof) on the
      ! nested PEs instead of sending it around.
      !if (gid == Atm(pp)%parent_grid%pelist(1)) then
      if (debug_log) print '("[INFO] WDR INIT_GRID XFER AP4 fv_moving_nest.F90 npe=",I0," send to pe=",I0," size=",I0)', this_pe, Atm(pp)%pelist(1), size(Atm(n)%grid_global(isg-ng:ieg+1+ng,jsg-ng:jeg+1+ng,1:2,parent_tile))

      call mpp_send(Atm(n)%grid_global(isg-ng:ieg+1+ng,jsg-ng:jeg+1+ng,1:2,parent_tile), &
          size(Atm(n)%grid_global(isg-ng:ieg+1+ng,jsg-ng:jeg+1+ng,1:2,parent_tile)), &
          Atm(pp)%pelist(1)) !send to p_ind in setup_aligned_nest
      if (debug_log) print '("[INFO] WDR INIT_GRID AP5 fv_moving_nest.F90 npe=",I0)', this_pe
      call mpp_sync_self()
      if (debug_log) print '("[INFO] WDR INIT_GRID AP6 fv_moving_nest.F90 npe=",I0)', this_pe
      !endif
    endif

    if (debug_log) print '("[INFO] WDR INIT_GRID AP9 fv_moving_nest.F90 npe=",I0)', this_pe

    !if (ngrids > 1) call setup_update_regions   ! Originally from fv_control.F90
    call mn_setup_update_regions(Atm, n, nest_domain)

    if (use_timers) call mpp_clock_end (id_reset3)
    if (use_timers) call mpp_clock_begin (id_reset4)

    if (Atm(n)%neststruct%nested) then
      if (debug_log) print '("[INFO] WDR INIT_GRID setup_aligned_nestA fv_moving_nest.F90 npe=",I0," n=",I0)', this_pe, n

      ! New code from fv_control.F90
      ! call init_grid(Atm(this_grid), Atm(this_grid)%flagstruct%grid_name, Atm(this_grid)%flagstruct%grid_file, &
      !    Atm(this_grid)%flagstruct%npx, Atm(this_grid)%flagstruct%npy, Atm(this_grid)%flagstruct%npz, Atm(this_grid)%flagstruct%ndims, Atm(this_grid)%flagstruct%ntiles, Atm(this_grid)%ng, tile_coarse)

      ! Atm(n)%neststruct%parent_tile            = tile_coarse(n)

      ! Old Code
      !call init_grid(Atm(n), Atm(n)%flagstruct%grid_name, Atm(n)%flagstruct%grid_file, &
      !     Atm(n)%npx, Atm(n)%npy, Atm(n)%npz, Atm(n)%flagstruct%ndims, Atm(n)%flagstruct%ntiles, Atm(n)%ng)

      !tile_coarse(1) = Atm(n)%neststruct%parent_tile
      tile_coarse(1) = parent_tile
      tile_coarse(2) = parent_tile

      call init_grid(Atm(n), Atm(n)%flagstruct%grid_name, Atm(n)%flagstruct%grid_file, &
          Atm(n)%flagstruct%npx, Atm(n)%flagstruct%npy, Atm(n)%flagstruct%npz, &
          Atm(n)%flagstruct%ndims, Atm(n)%flagstruct%ntiles, Atm(n)%ng, tile_coarse)
      if (debug_log) print '("[INFO] WDR INIT_GRID setup_aligned_nestB fv_moving_nest.F90 npe=",I0)', this_pe
    endif

    if (use_timers) call mpp_clock_end (id_reset4)
    if (use_timers) call mpp_clock_begin (id_reset5)

    !  Reset the gridstruct values for the nest
    if (is_fine_pe) then
      if (debug_log) print '("[INFO] WDR INIT_GRID AA fv_moving_nest.F90 npe=",I0)', this_pe
      if (debug_log) print '("[INFO] WDR INIT_GRID BB fv_moving_nest.F90 npe=",I0)', this_pe

      call grid_utils_init(Atm(n), Atm(n)%npx, Atm(n)%npy, Atm(n)%npz, &
          Atm(n)%flagstruct%non_ortho, Atm(n)%flagstruct%grid_type, Atm(n)%flagstruct%c2l_ord)

      if (debug_log) print '("[INFO] WDR INIT_GRID CC fv_moving_nest.F90 npe=",I0)', this_pe
    endif

    if (use_timers) call mpp_clock_end (id_reset5)
    if (use_timers) call mpp_clock_begin (id_reset6)

    if (debug_log) print '("[INFO] WDR NEST_DOMAIN ZZ fv_moving_nest.F90 npe=",I0)', this_pe

    if (debug_log) print '("[INFO] WDR REINIT1 CT fv_moving_nest.F90. npe=",I0," twowaynest=",L1" Atm(1)%neststruct%parent_tile=",I0)', &
        this_pe, Atm(1)%neststruct%twowaynest, Atm(1)%neststruct%parent_tile

    if (debug_log) print '("[INFO] WDR REINIT2 CT fv_moving_nest.F90. npe=",I0," twowaynest=",L1," Atm(2)%neststruct%parent_tile=",I0," n=",I0)', &
        this_pe, Atm(2)%neststruct%twowaynest, Atm(2)%neststruct%parent_tile, n

    !call mpp_sync(full_pelist)   ! Used to make debugging easier.  Can be removed.

    ! Needs to run for parent and nest Atm(2)
    !    Nest PEs update ind_update_h  -- this now seems obsolete
    !    Parent tile PEs update isu, ieu, jsu, jeu
    !    Global tiles that are not parent have no changes
    if (debug_log) print '("[INFO] WDR REINIT CV fv_moving_nest.F90. npe=",I0, " n=",I0)', this_pe, n

    ! WDR  This is now accomplished with the earlier call to setup_update_regions()
    !call reinit_parent_indices(Atm(2))
    !!call reinit_parent_indices(Atm(n))
    !if (debug_log) print '("[INFO] WDR REINIT CW fv_moving_nest.F90. npe=",I0)', this_pe

    do nn = 1, size(Atm)
      if (debug_log) call show_atm("3", Atm(nn), nn, this_pe)
    enddo


    ! Output the center lat/lon of the nest
    !   only the PE that holds the center point will output this information to the logfile
    !   lat = agrid(:,:,2) and lon = agrid(:,:,1), in radians
    if (is_fine_pe) then
      half_x = Atm(child_grid_num)%npx / 2
      half_y = Atm(child_grid_num)%npy / 2

      if (half_x .ge. Atm(child_grid_num)%bd%is .and. half_x .le. Atm(child_grid_num)%bd%ie .and. half_y .ge. Atm(child_grid_num)%bd%js .and. half_y .le. Atm(child_grid_num)%bd%je) then

        half_lat = Atm(child_grid_num)%gridstruct%agrid(half_x, half_y,2) * rad2deg
        half_lon = Atm(child_grid_num)%gridstruct%agrid(half_x, half_y,1) * rad2deg
        if (half_lon .gt. 180.0) half_lon = half_lon - 360.0

        print '("[INFO] fv_moving_nest.F90 NEST MOVED to npe=",I0," x=",I0," y=",I0," lat=",F6.2," lon=",F7.2," a_step=",I8," fcst_hr=",F12.3)', this_pe, \
        half_x, half_y, half_lat, half_lon, a_step,  a_step * dt_atmos / 3600.0
      endif

    endif

    ! Reallocate the halo buffers in the neststruct, as some are now the wrong size
    !   Optimization would be to only deallocate the edges that have changed.

    ! TODO Write comments on the t0 and t1 buffers
    if (use_timers) call mpp_clock_end (id_reset6)
    if (use_timers) call mpp_clock_begin (id_reset7)

    if (is_fine_pe) then
      !call reallocate_BC_buffers(Atm(child_grid_num))
      call reallocate_BC_buffers(Atm(1))
      if (debug_log) print '("[INFO] WDR INIT_GRID DD fv_moving_nest.F90 npe=",I0)', this_pe

      ! Reallocate buffers that are declared in fv_nesting.F90
      call dealloc_nested_buffers(Atm(1))

      if (debug_log) print '("[INFO] WDR INIT_GRID EE fv_moving_nest.F90 npe=",I0)', this_pe

      ! Set both to true so the call to setup_nested_grid_BCs() (at the beginning of fv_dynamics()) will reset t0 buffers
      ! They will be returned to false by setup_nested_grid_BCs()

      if (debug_log) print '("[INFO] WDR RESET_BCs first_step=.true.  fv_moving_nest.F90 npe=",I0)', this_pe
      Atm(n)%neststruct%first_step = .true.
      !Atm(n)%flagstruct%make_nh= .true.

      !! Fill in the BC time1 buffers
      !call setup_nested_grid_BCs(npx, npy, npz, zvir, ncnst, &
      !     u, v, w, pt, delp, delz, q, uc, vc, pkz, &
      !     neststruct%nested, flagstruct%inline_q, flagstruct%make_nh, ng, &
      !     gridstruct, flagstruct, neststruct, &
      !     neststruct%nest_timestep, neststruct%tracer_nest_timestep, &
      !     domain, bd, nwat)

      ! Transfer the BC time1 buffers to time0

      !call set_NH_BCs_t0(neststruct)
      !call set_BCs_t0(ncnst, flagstruct%hydrostatic, neststruct)

    endif
    if (use_timers) call mpp_clock_end (id_reset7)

  end subroutine mn_meta_reset_gridstruct


  ! Copied and adapted from fv_control.F90::setup_update_regions(); where it is an internal subroutine
  ! Modifications only to pass necessary variables as arguments

  !>@brief The subroutine 'mn_setup_update_regions' performs some of the tasks of fv_control.F90::setup_update_regions() for nest motion
  !>@details This routine only updates indices, so is computationally efficient
  subroutine mn_setup_update_regions(Atm, this_grid, nest_domain)
    type(fv_atmos_type), allocatable, intent(INOUT) :: Atm(:)                 !< Array of atmospheric data
    integer, intent(IN)                             :: this_grid              !< Parent or child grid number
    type(nest_domain_type),     intent(in)          :: nest_domain            !< Nest domain structure

    integer :: isu, ieu, jsu, jeu ! update regions
    integer :: isc, jsc, iec, jec
    integer :: upoff
    integer :: ngrids, n, nn
    integer :: isu_stag, isv_stag, jsu_stag, jsv_stag
    integer :: ieu_stag, iev_stag, jeu_stag, jev_stag
    integer :: this_pe

    this_pe = mpp_pe()

    ! Need to get the following variables from nest_domain
    !   tile_coarse()
    !   icount_coarse()
    !         from mpp_define_nest_domains.inc:  iend_coarse(n) = istart_coarse(n) + icount_coarse(n) - 1
    !         rearrange to: iend_coarse(n) - istart_coarse(n) + 1 = icount_coarse(n)
    !   jcount_coarse()
    !   nest_ioffsets()
    !      in fv_control.F90. pass nest_ioffsets as istart_coarse
    !   nest_joffsets()

    isc = Atm(this_grid)%bd%isc
    jsc = Atm(this_grid)%bd%jsc
    iec = Atm(this_grid)%bd%iec
    jec = Atm(this_grid)%bd%jec

    upoff = Atm(this_grid)%neststruct%upoff

    ngrids = size(Atm)

    if (debug_log) print '("[INFO] WDR SUR fv_moving_nest.F90. npe=",I0," ngrids=",I0," nest_domain%tile_coarse(",I0,"-",I0,")")', this_pe, ngrids, lbound(nest_domain%tile_coarse), ubound(nest_domain%tile_coarse)

    if (debug_log) print '("[INFO] WDR tile_coarse fv_moving_nest.F90 npe=",I0," tile_coarse(",I0,"-",I0") ngrids=",I0," tile_coarse(1)=",I0)', this_pe, &
        lbound(nest_domain%tile_coarse,1), ubound(nest_domain%tile_coarse,1), ngrids, nest_domain%tile_coarse(1)

    if (debug_log) print '("[INFO] WDR tile_coarse fv_moving_nest.F90 npe=",I0," istart_coarse(",I0,"-",I0")")', this_pe, &
        lbound(nest_domain%istart_coarse,1), ubound(nest_domain%istart_coarse,1)

    do n=2,ngrids
      nn = n - 1  !  WDR TODO revise this to handle multiple nests.  This adjusts to match fv_control.F90 where these
      !  arrays are passed in to mpp_define_nest_domains with bounds (2:ngrids)

      ! Updated code from new fv_control.F90   November 8. 2021 Ramstrom

      if (nest_domain%tile_coarse(nn) == Atm(this_grid)%global_tile) then

        !isu = nest_ioffsets(n)
        isu = nest_domain%istart_coarse(nn)
        !ieu = isu + icount_coarse(n) - 1
        ieu = isu + (nest_domain%iend_coarse(nn) - nest_domain%istart_coarse(nn) + 1) - 1

        !jsu = nest_joffsets(n)
        jsu = nest_domain%jstart_coarse(nn)
        !jeu = jsu + jcount_coarse(n) - 1
        jeu = jsu + (nest_domain%jend_coarse(nn) - nest_domain%jstart_coarse(nn) + 1) - 1

!!! Begin new
        isu_stag = isu
        jsu_stag = jsu
        ieu_stag = ieu
        jeu_stag = jeu+1

        isv_stag = isu
        jsv_stag = jsu
        iev_stag = ieu+1
        jev_stag = jeu
!!! End new


        !update offset adjustment
        isu = isu + upoff
        ieu = ieu - upoff
        jsu = jsu + upoff
        jeu = jeu - upoff

!!! Begin new
        isu_stag = isu_stag + upoff
        ieu_stag = ieu_stag - upoff
        jsu_stag = jsu_stag + upoff
        jeu_stag = jeu_stag - upoff

        isv_stag = isv_stag + upoff
        iev_stag = iev_stag - upoff
        jsv_stag = jsv_stag + upoff
        jev_stag = jev_stag - upoff

        ! Absolute boundary for the staggered point update region on the parent.
        ! This is used in remap_uv to control the update of the last staggered point
        ! when the the update region coincides with a pe domain to avoid cross-restart repro issues

        Atm(n)%neststruct%jeu_stag_boundary = jeu_stag
        Atm(n)%neststruct%iev_stag_boundary = iev_stag

        if (isu > iec .or. ieu < isc .or. &
            jsu > jec .or. jeu < jsc ) then
          isu = -999 ; jsu = -999 ; ieu = -1000 ; jeu = -1000
        else
          isu = max(isu,isc) ; jsu = max(jsu,jsc)
          ieu = min(ieu,iec) ; jeu = min(jeu,jec)
        endif

        ! Update region for staggered quantity to avoid cross repro issues when the pe domain boundary
        ! coincide with the nest. Basically write the staggered update on compute domains

        if (isu_stag > iec .or. ieu_stag < isc .or. &
            jsu_stag > jec .or. jeu_stag < jsc ) then
          isu_stag = -999 ; jsu_stag = -999 ; ieu_stag = -1000 ; jeu_stag = -1000
        else
          isu_stag = max(isu_stag,isc) ; jsu_stag = max(jsu_stag,jsc)
          ieu_stag = min(ieu_stag,iec) ; jeu_stag = min(jeu_stag,jec)
        endif

        if (isv_stag > iec .or. iev_stag < isc .or. &
            jsv_stag > jec .or. jev_stag < jsc ) then
          isv_stag = -999 ; jsv_stag = -999 ; iev_stag = -1000 ; jev_stag = -1000
        else
          isv_stag = max(isv_stag,isc) ; jsv_stag = max(jsv_stag,jsc)
          iev_stag = min(iev_stag,iec) ; jev_stag = min(jev_stag,jec)
        endif
!!! End new

        if (isu > iec .or. ieu < isc .or. &
            jsu > jec .or. jeu < jsc ) then
          isu = -999 ; jsu = -999 ; ieu = -1000 ; jeu = -1000
        else
          isu = max(isu,isc) ; jsu = max(jsu,jsc)
          ieu = min(ieu,iec) ; jeu = min(jeu,jec)
        endif

        ! lump indices
        isu=max(isu, isu_stag, isv_stag)
        jsu=max(jsu, jsu_stag, jsv_stag)
        jeu_stag=max(jeu, jeu_stag)
        jev_stag=max(jeu, jev_stag)
        ieu_stag=max(ieu ,ieu_stag)
        iev_stag=max(ieu ,iev_stag)

        Atm(n)%neststruct%isu = isu
        Atm(n)%neststruct%ieu = ieu_stag
        Atm(n)%neststruct%jsu = jsu
        Atm(n)%neststruct%jeu = jev_stag

        Atm(n)%neststruct%jeu_stag = jeu_stag
        Atm(n)%neststruct%iev_stag = iev_stag
      endif
    enddo

  end subroutine mn_setup_update_regions


  !==================================================================================================
  !
  !  Recalculation Section -- Buffers that have to change size after nest motion
  !
  !==================================================================================================

  !>@brief The subroutine 'reallocate_BC_buffers' reallocates boundary condition buffers - some need to change size after a nest move.
  !>@details Thought they would be reallocated in boundary.F90 nested_grid_BC_recv() when needed, but seem not to.
  subroutine reallocate_BC_buffers(Atm)
    type(fv_atmos_type), intent(inout)  :: Atm     !< Single instance of atmospheric data

    integer :: n, ns
    logical :: dummy = .false. ! same as grids_on_this_pe(n)

    call deallocate_fv_nest_BC_type(Atm%neststruct%delp_BC)
    call deallocate_fv_nest_BC_type(Atm%neststruct%u_BC)
    call deallocate_fv_nest_BC_type(Atm%neststruct%v_BC)
    call deallocate_fv_nest_BC_type(Atm%neststruct%uc_BC)
    call deallocate_fv_nest_BC_type(Atm%neststruct%vc_BC)
    call deallocate_fv_nest_BC_type(Atm%neststruct%divg_BC)

    if (allocated(Atm%neststruct%q_BC)) then
      do n=1,size(Atm%neststruct%q_BC)
        call deallocate_fv_nest_BC_type(Atm%neststruct%q_BC(n))
      enddo
    endif

#ifndef SW_DYNAMICS
    call deallocate_fv_nest_BC_type(Atm%neststruct%pt_BC)
#ifdef USE_COND
    call deallocate_fv_nest_BC_type(Atm%neststruct%q_con_BC)
#ifdef MOIST_CAPPA
    call deallocate_fv_nest_BC_type(Atm%neststruct%cappa_BC)
#endif
#endif
    if (.not.Atm%flagstruct%hydrostatic) then
      call deallocate_fv_nest_BC_type(Atm%neststruct%w_BC)
      call deallocate_fv_nest_BC_type(Atm%neststruct%delz_BC)
    endif
#endif

    ! Reallocate the buffers

    ns = Atm%neststruct%nsponge

    call allocate_fv_nest_BC_type(Atm%neststruct%delp_BC,Atm,ns,0,0,dummy)
    call allocate_fv_nest_BC_type(Atm%neststruct%u_BC,Atm,ns,0,1,dummy)
    call allocate_fv_nest_BC_type(Atm%neststruct%v_BC,Atm,ns,1,0,dummy)
    call allocate_fv_nest_BC_type(Atm%neststruct%uc_BC,Atm,ns,1,0,dummy)
    call allocate_fv_nest_BC_type(Atm%neststruct%vc_BC,Atm,ns,0,1,dummy)
    call allocate_fv_nest_BC_type(Atm%neststruct%divg_BC,Atm,ns,1,1,dummy)

    !  if (ncnst > 0) then
    !     allocate(Atm%neststruct%q_BC(ncnst))
    !     do n=1,ncnst
    !        call allocate_fv_nest_BC_type(Atm%neststruct%q_BC(n),Atm,ns,0,0,dummy)
    !     enddo
    !  endif

    if (allocated(Atm%neststruct%q_BC)) then
      do n=1,size(Atm%neststruct%q_BC)
        call allocate_fv_nest_BC_type(Atm%neststruct%q_BC(n),Atm,ns,0,0,dummy)
      enddo
    endif

#ifndef SW_DYNAMICS
    call allocate_fv_nest_BC_type(Atm%neststruct%pt_BC,Atm,ns,0,0,dummy)
#ifdef USE_COND
    call allocate_fv_nest_BC_type(Atm%neststruct%q_con_BC,Atm,ns,0,0,dummy)
#ifdef MOIST_CAPPA
    call allocate_fv_nest_BC_type(Atm%neststruct%cappa_BC,Atm,ns,0,0,dummy)
#endif
#endif
    if (.not.Atm%flagstruct%hydrostatic) then
      call allocate_fv_nest_BC_type(Atm%neststruct%w_BC,Atm,ns,0,0,dummy)
      call allocate_fv_nest_BC_type(Atm%neststruct%delz_BC,Atm,ns,0,0,dummy)
    endif
#endif

  end subroutine reallocate_BC_buffers


  !!============================================================================
  !!  Step 8 -- Moving Nest Output to NetCDF
  !!============================================================================

  !>@brief The subroutine 'mn_prog_dump_to_netcdf' dumps selected prognostic variables to netCDF file.
  !>@details Can be modified to output more of the prognostic variables if wanted.  Certain 3D variables were commented out for performance.
  subroutine mn_prog_dump_to_netcdf(Atm, time_val, file_prefix, is_fine_pe, domain_coarse, domain_fine, nz)
    type(fv_atmos_type), intent(in)            :: Atm                                  !< Single instance of atmospheric data
    integer, intent(in)                        :: time_val                             !< Timestep number
    character(len=*), intent(in)               :: file_prefix                          !< Filename prefix
    logical, intent(in)                        :: is_fine_pe                           !< Is nest PE?
    type(domain2d), intent(in)                 :: domain_coarse, domain_fine           !< Domain structures
    integer, intent(in)                        :: nz                                   !< Number of vertical levels

    integer            :: n_moist
    character(len=16)  :: out_var_name
    integer            :: position = CENTER
    !integer            :: position_u = NORTH
    !integer            :: position_v = EAST

    call mn_var_dump_to_netcdf(Atm%pt   , is_fine_pe, domain_coarse, domain_fine, position, nz, &
        time_val, Atm%global_tile, file_prefix, "tempK")
    call mn_var_dump_to_netcdf(Atm%pt(:,:,64)   , is_fine_pe, domain_coarse, domain_fine, position, nz, &
        time_val, Atm%global_tile, file_prefix, "T64")
    !call mn_var_dump_to_netcdf(Atm%delp , is_fine_pe, domain_coarse, domain_fine, position, nz, &
    !     time_val, Atm%global_tile, file_prefix, "DELP")
    call mn_var_dump_to_netcdf(Atm%delz , is_fine_pe, domain_coarse, domain_fine, position, nz, &
        time_val, Atm%global_tile, file_prefix, "DELZ")
    call mn_var_dump_to_netcdf(Atm%q_con, is_fine_pe, domain_coarse, domain_fine, position, nz, &
        time_val, Atm%global_tile, file_prefix, "qcon")

    !call mn_var_dump_to_netcdf(Atm%w    , is_fine_pe, domain_coarse, domain_fine, position, nz, &
    !     time_val, Atm%global_tile, file_prefix, "WWND")
    !call mn_var_dump_to_netcdf(Atm%ua   , is_fine_pe, domain_coarse, domain_fine, position, nz, &
    !     time_val, Atm%global_tile, file_prefix, "UA")
    !call mn_var_dump_to_netcdf(Atm%va   , is_fine_pe, domain_coarse, domain_fine, position, nz, &
    !     time_val, Atm%global_tile, file_prefix, "VA")

    call mn_var_dump_to_netcdf(Atm%ps   , is_fine_pe, domain_coarse, domain_fine, position, 1 , &
        time_val, Atm%global_tile, file_prefix, "PS")

    !! TODO figure out what to do with ze0;  different bounds - only compute domain

    !! TODO Wind worked fine when in its own file.  Can it merge in with the regular file??
    !!call mn_var_dump_to_netcdf(Atm%u, is_fine_pe, domain_coarse, domain_fine, position_u, nz, &
    !!     time_val, Atm%global_tile, "wxvarU", "UWND")
    !!call mn_var_dump_to_netcdf(Atm%v, is_fine_pe, domain_coarse, domain_fine, position_v, nz, &
    !!     time_val, Atm%global_tile, "wxvarU", "VWND")

    ! Latitude and longitude in radians
    call mn_var_dump_to_netcdf( Atm%gridstruct%agrid(:,:,2), is_fine_pe, domain_coarse, domain_fine, position, nz, &
        time_val, Atm%global_tile, file_prefix, "latrad")
    call mn_var_dump_to_netcdf( Atm%gridstruct%agrid(:,:,1), is_fine_pe, domain_coarse, domain_fine, position, nz, &
        time_val, Atm%global_tile, file_prefix, "lonrad")

    !do n_moist = lbound(Atm%q, 4), ubound(Atm%q, 4)
    !   call get_tracer_names(MODEL_ATMOS, n_moist, out_var_name)
    !   call mn_var_dump_to_netcdf( Atm%q(:,:,:,n_moist), is_fine_pe, domain_coarse, domain_fine, position, nz, &
    !        time_val, Atm%global_tile, file_prefix, trim(out_var_name))
    !enddo

  end subroutine mn_prog_dump_to_netcdf


  !!  Step 8 -- Moving Nest Output Individual Variables

  !>@brief The subroutine 'mn_var_dump_3d_to_netcdf' dumps a 3D single precision variable to netCDF file.
  subroutine mn_var_dump_3d_to_netcdf( data_var, is_fine_pe, domain_coarse, domain_fine, position, nz, time_step, this_tile, file_prefix, var_name)
    real, intent(in)                            :: data_var(:,:,:)                      !< Single precision model variable
    logical, intent(in)                         :: is_fine_pe                           !< Is nest PE?
    type(domain2d), intent(in)                  :: domain_coarse, domain_fine           !< Domain structures
    integer, intent(in)                         :: position, nz, time_step, this_tile   !< Stagger, number vertical levels, timestep, tile number
    character(len=*)                            :: file_prefix, var_name                !< Filename prefix, and netCDF variable name

    integer                      :: isd_coarse, ied_coarse, jsd_coarse, jed_coarse
    integer                      :: isd_fine, ied_fine, jsd_fine, jed_fine
    integer                      :: this_pe
    character(len=64)            :: prefix_fine, prefix_coarse

    this_pe = mpp_pe()

    prefix_fine = trim(file_prefix) // "_fine"
    prefix_coarse = trim(file_prefix) // "_coarse"

    !!===========================================================
    !!
    !! Output the grid data from both nest grids and parent grids to netCDF
    !!
    !!===========================================================

    if (is_fine_pe) then
      call mpp_get_data_domain(domain_fine, isd_fine, ied_fine, jsd_fine, jed_fine, position=position)

      if (debug_log) print '("[INFO] WDR NRF FG mn_var_dump_to_netcdf start. npe=",I0," size of x=(",I0,",",I0,",",I0")")', this_pe, size(data_var,1), size(data_var,2), size(data_var,3)
      if (debug_log) print '("[INFO] WDR NRF FG mn_var_dump_to_netcdf start. npe=",I0," Data domain i=",I0,"-",I0," j=",I0,"-",I0," (",I0,",",I0,")")', &
          this_pe, isd_fine, ied_fine, jsd_fine, jed_fine, ied_fine - isd_fine + 1,  jed_fine - jsd_fine + 1

      call output_grid_to_nc("GH", isd_fine, ied_fine, jsd_fine, jed_fine, nz, data_var, prefix_fine, var_name, time_step, domain_fine, position)

    else
      if (this_tile == 6) then
        !call mpp_get_compute_domain(domain_coarse, isc_coarse, iec_coarse, jsc_coarse, jec_coarse, position=position)
        call mpp_get_data_domain(domain_coarse, isd_coarse, ied_coarse, jsd_coarse, jed_coarse, position=position)
        !call mpp_get_memory_domain(domain_coarse, ism_coarse, iem_coarse, jsm_coarse, jem_coarse, position=position)

        if (debug_log) print '("[INFO] WDR NRF CG mn_var_dump_to_netcdf start. npe=",I0," size of x=(",I0,",",I0,",",I0")")', this_pe, size(data_var,1), size(data_var,2), size(data_var,3)
        if (debug_log) print '("[INFO] WDR NRF CG mn_var_dump_to_netcdf start. npe=",I0," Data domain i=",I0,"-",I0," j=",I0,"-",I0," (",I0,",",I0,")")', &
            this_pe, isd_coarse, ied_coarse, jsd_coarse, jed_coarse, ied_coarse - isd_coarse + 1,  jed_coarse - jsd_coarse + 1
        !if (debug_log) print '("[INFO] WDR NRF CG mn_var_dump_to_netcdf start. npe=",I0," Compute domain i=",I0,"-",I0," j=",I0,"-",I0," (",I0,",",I0,")")', &
        !     this_pe, isc_coarse, iec_coarse, jsc_coarse, jec_coarse, iec_coarse - isc_coarse + 1,  jec_coarse - jsc_coarse + 1
        !if (debug_log) print '("[INFO] WDR NRF CG mn_var_dump_to_netcdf start. npe=",I0," Memory domain i=",I0,"-",I0," j=",I0,"-",I0," (",I0,",",I0,")")', &
        !     this_pe, ism_coarse, iem_coarse, jsm_coarse, jem_coarse, iem_coarse - ism_coarse + 1,  jem_coarse - jsm_coarse + 1

        call output_grid_to_nc("GH", isd_coarse, ied_coarse, jsd_coarse, jed_coarse, nz, data_var, prefix_coarse, var_name, time_step, domain_coarse, position)

      endif
    endif

  end subroutine mn_var_dump_3d_to_netcdf

  !>@brief The subroutine 'mn_var_dump_2d_to_netcdf' dumps a 3D single precision variable to netCDF file.
  subroutine mn_var_dump_2d_to_netcdf( data_var, is_fine_pe, domain_coarse, domain_fine, position, nz, time_step, this_tile, file_prefix, var_name)
    implicit none
    real, intent(in)                            :: data_var(:,:)                         !< Data variable
    logical, intent(in)                         :: is_fine_pe                            !< Is nest PE?
    type(domain2d), intent(in)                  :: domain_coarse, domain_fine            !< Domain structures
    integer, intent(in)                         :: position, nz, time_step, this_tile    !< Stagger, number vertical levels, timestep, tile number
    character(len=*)                            :: file_prefix, var_name                 !< Filename prefix, and netCDF variable name

    integer                      :: isc_coarse, iec_coarse, jsc_coarse, jec_coarse
    integer                      :: isd_coarse, ied_coarse, jsd_coarse, jed_coarse
    integer                      :: isd_fine, ied_fine, jsd_fine, jed_fine
    integer                      :: isc_fine, iec_fine, jsc_fine, jec_fine

    integer                      :: ism_coarse, iem_coarse, jsm_coarse, jem_coarse
    integer                      :: ism_fine, iem_fine, jsm_fine, jem_fine

    integer                      :: this_pe

    character(len=64)            :: prefix_fine, prefix_coarse

    this_pe = mpp_pe()

    prefix_fine = trim(file_prefix) // "_fine"
    prefix_coarse = trim(file_prefix) // "_coarse"

    !!===========================================================
    !!
    !! Output the grid data from both nest grids and parent grids to netCDF
    !!
    !!===========================================================

    if (is_fine_pe) then
      ! Maybe don't need to call mpp_get_compute_domain here?
      !call mpp_get_compute_domain(domain_fine, isc_fine, iec_fine, jsc_fine, jec_fine, position=position)
      call mpp_get_data_domain(domain_fine, isd_fine, ied_fine, jsd_fine, jed_fine, position=position)
      !call mpp_get_memory_domain(domain_fine, ism_fine, iem_fine, jsm_fine, jem_fine, position=position)

      if (debug_log) print '("[INFO] WDR NRF FG mn_var_dump_to_netcdf start. npe=",I0," size of x=(",I0,",",I0,")")', this_pe, size(data_var,1), size(data_var,2)
      if (debug_log) print '("[INFO] WDR NRF FG mn_var_dump_to_netcdf start. npe=",I0," Data domain i=",I0,"-",I0," j=",I0,"-",I0," (",I0,",",I0,")")', &
          this_pe, isd_fine, ied_fine, jsd_fine, jed_fine, ied_fine - isd_fine + 1,  jed_fine - jsd_fine + 1
      !if (debug_log) print '("[INFO] WDR NRF FG mn_var_dump_to_netcdf start. npe=",I0," Compute domain i=",I0,"-",I0," j=",I0,"-",I0," (",I0,",",I0,")")', &
      !     this_pe, isc_fine, iec_fine, jsc_fine, jec_fine, iec_fine - isc_fine + 1,  jec_fine - jsc_fine + 1
      !if (debug_log) print '("[INFO] WDR NRF FG mn_var_dump_to_netcdf start. npe=",I0," Memory domain i=",I0,"-",I0," j=",I0,"-",I0," (",I0,",",I0,")")', &
      !     this_pe, ism_fine, iem_fine, jsm_fine, jem_fine, iem_fine - ism_fine + 1,  jem_fine - jsm_fine + 1

      call output_grid_to_nc("GH", isd_fine, ied_fine, jsd_fine, jed_fine, nz, data_var, prefix_fine, var_name, time_step, domain_fine, position)
    else

      if (this_tile == 6) then
        !call mpp_get_compute_domain(domain_coarse, isc_coarse, iec_coarse, jsc_coarse, jec_coarse, position=position)
        call mpp_get_data_domain(domain_coarse, isd_coarse, ied_coarse, jsd_coarse, jed_coarse, position=position)
        !call mpp_get_memory_domain(domain_coarse, ism_coarse, iem_coarse, jsm_coarse, jem_coarse, position=position)

        if (debug_log) print '("[INFO] WDR NRF CG mn_var_dump_to_netcdf start. npe=",I0," size of x=(",I0,",",I0,")")', this_pe, size(data_var,1), size(data_var,2)
        if (debug_log) print '("[INFO] WDR NRF CG mn_var_dump_to_netcdf start. npe=",I0," Data domain i=",I0,"-",I0," j=",I0,"-",I0," (",I0,",",I0,")")', &
            this_pe, isd_coarse, ied_coarse, jsd_coarse, jed_coarse, ied_coarse - isd_coarse + 1,  jed_coarse - jsd_coarse + 1
        !if (debug_log) print '("[INFO] WDR NRF CG mn_var_dump_to_netcdf start. npe=",I0," Compute domain i=",I0,"-",I0," j=",I0,"-",I0," (",I0,",",I0,")")', &
        !     this_pe, isc_coarse, iec_coarse, jsc_coarse, jec_coarse, iec_coarse - isc_coarse + 1,  jec_coarse - jsc_coarse + 1
        !if (debug_log) print '("[INFO] WDR NRF CG mn_var_dump_to_netcdf start. npe=",I0," Memory domain i=",I0,"-",I0," j=",I0,"-",I0," (",I0,",",I0,")")', &
        !     this_pe, ism_coarse, iem_coarse, jsm_coarse, jem_coarse, iem_coarse - ism_coarse + 1,  jem_coarse - jsm_coarse + 1

        call output_grid_to_nc("GH", isd_coarse, ied_coarse, jsd_coarse, jed_coarse, nz, data_var, prefix_coarse, var_name, time_step, domain_coarse, position)

      endif
    endif

  end subroutine mn_var_dump_2d_to_netcdf


  !!=========================================================================================
  !! Step 9 -- Perform vertical remapping on nest(s) and recalculate auxiliary pressures
  !!           Should help stabilize the fields before dynamics runs
  !!=========================================================================================

  !>@brief The subroutine 'recalc_aux_pressures' updates auxiliary pressures after a nest move.
  subroutine recalc_aux_pressures(Atm)
    type(fv_atmos_type), intent(inout) :: Atm      !< Single Atm structure

    !  Update the auxiliary pressure variables
    !  In nest moving code, we moved delp and delz; this will update ps, pk, pe, peln, and pkz
    !  Note this routine makes hydrostatic calculations (but has non-hydrostatic branches)
    !  Perhaps not appropriate for a non-hydrostatic run.
    !  May need to find or write a non-hydrostatic version of this routine

    ! TODO determine if this is the correct way to recalculate the auxiliary pressure variables

    call p_var(Atm%npz, Atm%bd%is, Atm%bd%ie, Atm%bd%js, Atm%bd%je, Atm%ptop, ptop_min,  &
        Atm%delp, Atm%delz, &
        Atm%pt, Atm%ps, &
        Atm%pe, Atm%peln,   &
        Atm%pk,   Atm%pkz, kappa, &
        Atm%q, Atm%ng, Atm%flagstruct%ncnst, Atm%gridstruct%area_64, 0.,  &
        .false.,  .false., & !mountain argument not used
        Atm%flagstruct%moist_phys,  Atm%flagstruct%hydrostatic, &
        Atm%flagstruct%nwat, Atm%domain, .false.)

  end subroutine recalc_aux_pressures


  !==================================================================================================
  !
  !  Utility Section  -- After Step 9
  !
  !==================================================================================================

  !>@brief The subroutine 'init_ijk_mem' was copied from dyn_core.F90 to avoid circular dependencies
  subroutine init_ijk_mem(i1, i2, j1, j2, km, array, var)
    integer, intent(in):: i1, i2, j1, j2, km
    real, intent(inout):: array(i1:i2,j1:j2,km)
    real, intent(in):: var
    integer:: i, j, k

    !$OMP parallel do default(none) shared(i1,i2,j1,j2,km,array,var)
    do k=1,km
      do j=j1,j2
        do i=i1,i2
          array(i,j,k) = var
        enddo
      enddo
    enddo

  end subroutine init_ijk_mem

  !>@brief The function 'almost_equal' tests whether real values are within a tolerance of one another.
  function almost_equal(a, b)
    logical :: almost_equal
    real, intent(in):: a,b

    real :: tolerance = 0.00001

    if ( abs(a - b) < tolerance ) then
      almost_equal = .true.
    else
      almost_equal = .false.
    endif
  end function almost_equal



  !>@brief The subroutine 'move_nest_geo' shifts tile_geo values using the data from fp_super_tile_geo
  subroutine move_nest_geo(tile_geo, tile_geo_u, tile_geo_v, fp_super_tile_geo, delta_i_c, delta_j_c, x_refine, y_refine)
    implicit none
    type(grid_geometry), intent(inout)  :: tile_geo                                     !< A-grid tile geometry
    type(grid_geometry), intent(inout)  :: tile_geo_u                                   !< u-wind tile geometry
    type(grid_geometry), intent(inout)  :: tile_geo_v                                   !< v-wind tile geometry
    type(grid_geometry), intent(in)     :: fp_super_tile_geo                            !< Parent high-resolution supergrid tile geometry
    integer, intent(in)                 :: delta_i_c, delta_j_c, x_refine, y_refine     !< delta i,j for nest move.  Nest refinement.

    integer :: nest_x, nest_y, parent_x, parent_y

    type(bbox)  :: tile_bbox, fp_tile_bbox, tile_bbox_u, tile_bbox_v
    integer   :: i, j, fp_i, fp_j

    ! tile_geo is cell-centered, at nest refinement
    ! fp_super_tile_geo is a supergrid, at nest refinement

    call find_nest_alignment(tile_geo, fp_super_tile_geo, nest_x, nest_y, parent_x, parent_y)

    call fill_bbox(tile_bbox, tile_geo%lats)
    call fill_bbox(tile_bbox_u, tile_geo_u%lats)
    call fill_bbox(tile_bbox_v, tile_geo_v%lats)
    call fill_bbox(fp_tile_bbox, fp_super_tile_geo%lats)

    ! Calculate new parent alignment -- supergrid at the refine ratio
    !  delta_{i,j}_c are at the coarse center grid resolution
    parent_x = parent_x + delta_i_c * 2 * x_refine
    parent_y = parent_y + delta_j_c * 2 * y_refine

    ! Brute force repopulation of full tile_geo grids.
    ! Optimization would be to use EOSHIFT and bring in just leading edge
    do i = tile_bbox%is, tile_bbox%ie
      do j = tile_bbox%js, tile_bbox%je
        fp_i = (i - nest_x) * 2 + parent_x
        fp_j = (j - nest_y) * 2 + parent_y

        if (fp_i < fp_tile_bbox%is .or. fp_i > fp_tile_bbox%ie) then
          if (debug_log) print '("[ERROR] WDR move_nest_geo invalid fp_i=",I0," is=",I0," ie=",I0)', fp_i, fp_tile_bbox%is, fp_tile_bbox%ie
          stop  ! replace with a fatal error
        endif
        if (fp_j < fp_tile_bbox%js .or. fp_j > fp_tile_bbox%je) then
          if (debug_log) print '("[ERROR] WDR move_nest_geo invalid fp_j=",I0," js=",I0," je=",I0)', fp_j, fp_tile_bbox%js, fp_tile_bbox%je
          stop  ! replace with a fatal error
        endif

        tile_geo%lats(i,j) = fp_super_tile_geo%lats(fp_i, fp_j)
        tile_geo%lons(i,j) = fp_super_tile_geo%lons(fp_i, fp_j)
      enddo
    enddo

    do i = tile_bbox_u%is, tile_bbox_u%ie
      do j = tile_bbox_u%js, tile_bbox_u%je
        fp_i = (i - nest_x) * 2 + parent_x
        fp_j = (j - nest_y) * 2 + parent_y - 1

        if (fp_i < fp_tile_bbox%is .or. fp_i > fp_tile_bbox%ie) then
          if (debug_log) print '("[ERROR] WDR move_nest_geo invalid fp_i=",I0," is=",I0," ie=",I0)', fp_i, fp_tile_bbox%is, fp_tile_bbox%ie
          stop  ! replace with a fatal error
        endif
        if (fp_j < fp_tile_bbox%js .or. fp_j > fp_tile_bbox%je) then
          if (debug_log) print '("[ERROR] WDR move_nest_geo invalid fp_j=",I0," js=",I0," je=",I0)', fp_j, fp_tile_bbox%js, fp_tile_bbox%je
          stop  ! replace with a fatal error
        endif

        tile_geo_u%lats(i,j) = fp_super_tile_geo%lats(fp_i, fp_j)
        tile_geo_u%lons(i,j) = fp_super_tile_geo%lons(fp_i, fp_j)
      enddo
    enddo

    do i = tile_bbox_v%is, tile_bbox_v%ie
      do j = tile_bbox_v%js, tile_bbox_v%je
        fp_i = (i - nest_x) * 2 + parent_x - 1
        fp_j = (j - nest_y) * 2 + parent_y

        if (fp_i < fp_tile_bbox%is .or. fp_i > fp_tile_bbox%ie) then
          if (debug_log) print '("[ERROR] WDR move_nest_geo invalid fp_i=",I0," is=",I0," ie=",I0)', fp_i, fp_tile_bbox%is, fp_tile_bbox%ie
          stop  ! replace with a fatal error
        endif
        if (fp_j < fp_tile_bbox%js .or. fp_j > fp_tile_bbox%je) then
          if (debug_log) print '("[ERROR] WDR move_nest_geo invalid fp_j=",I0," js=",I0," je=",I0)', fp_j, fp_tile_bbox%js, fp_tile_bbox%je
          stop  ! replace with a fatal error
        endif

        tile_geo_v%lats(i,j) = fp_super_tile_geo%lats(fp_i, fp_j)
        tile_geo_v%lons(i,j) = fp_super_tile_geo%lons(fp_i, fp_j)
      enddo
    enddo

    ! Validate at the end
    call find_nest_alignment(tile_geo, fp_super_tile_geo, nest_x, nest_y, parent_x, parent_y)

  end subroutine move_nest_geo

  !>@brief The subroutine 'assign_n_p_grids' sets values for parent and nest grid arrays from the grid_geometry structures.
  subroutine assign_n_p_grids(parent_geo, tile_geo, p_grid, n_grid, position)
    type(grid_geometry), intent(in)          ::  parent_geo, tile_geo                        !< Parent geometry, nest geometry
    real(kind=R_GRID), allocatable, intent(inout)            :: p_grid(:,:,:)                !< Parent grid
    real(kind=R_GRID), allocatable, intent(inout)            :: n_grid(:,:,:)                !< Nest grid
    integer, intent(in)                      :: position                                     !< Grid offset

    integer :: i,j

    if (position == CENTER) then
      do j = lbound(tile_geo%lats,2), ubound(tile_geo%lats,2)
        do i = lbound(tile_geo%lats,1), ubound(tile_geo%lats,1)
          ! centered grid version
          n_grid(i, j, 1) = tile_geo%lons(i, j)
          n_grid(i, j, 2) = tile_geo%lats(i, j)
          !if (debug_log) print '("[INFO] WDR populate ngrid npe=",I0, I4,I4, F12.4, F12.4)', this_pe, i, j, n_grid(i,j,1), n_grid(i,j,2)
        enddo
      enddo

      do j = 1, parent_geo%ny
        do i = 1, parent_geo%nx
          ! centered grid version
          p_grid(i, j, 1) = parent_geo%lons(2*i, 2*j)
          p_grid(i, j, 2) = parent_geo%lats(2*i, 2*j)
        enddo
      enddo

      ! u(npx, npy+1)
    elseif (position == NORTH) then  ! u wind on D-stagger
      do j = lbound(tile_geo%lats,2), ubound(tile_geo%lats,2)
        do i = lbound(tile_geo%lats,1), ubound(tile_geo%lats,1)
          ! centered grid version
          n_grid(i, j, 1) = tile_geo%lons(i, j)
          n_grid(i, j, 2) = tile_geo%lats(i, j)
          !if (debug_log) print '("[INFO] WDR populate ngrid_u npe=",I0, I4,I4, F12.4, F12.4)', this_pe, i, j, n_grid(i,j,1), n_grid(i,j,2)
        enddo
      enddo

      do j = 1, parent_geo%ny
        do i = 1, parent_geo%nx
          ! centered grid version
          p_grid(i, j, 1) = parent_geo%lons(2*i, 2*j-1)
          p_grid(i, j, 2) = parent_geo%lats(2*i, 2*j-1)
        enddo
      enddo

      ! v(npx+1, npy)
    elseif (position == EAST) then  ! v wind on D-stagger
      do j = lbound(tile_geo%lats,2), ubound(tile_geo%lats,2)
        do i = lbound(tile_geo%lats,1), ubound(tile_geo%lats,1)
          ! centered grid version
          n_grid(i, j, 1) = tile_geo%lons(i, j)
          n_grid(i, j, 2) = tile_geo%lats(i, j)
          !if (debug_log) print '("[INFO] WDR populate ngrid_v npe=",I0, I4,I4, F12.4, F12.4)', this_pe, i, j, n_grid(i,j,1), n_grid(i,j,2)
        enddo
      enddo

      do j = 1, parent_geo%ny
        do i = 1, parent_geo%nx
          ! centered grid version
          p_grid(i, j, 1) = parent_geo%lons(2*i-1, 2*j)
          p_grid(i, j, 2) = parent_geo%lats(2*i-1, 2*j)
        enddo
      enddo

    endif

  end subroutine assign_n_p_grids


  !>@brief The subroutine 'calc_nest_halo_weights' calculates the interpolation weights
  !>@details Computationally demanding; target for optimization after nest moves
  subroutine calc_nest_halo_weights(bbox_fine, bbox_coarse, p_grid, n_grid, wt, istart_coarse, jstart_coarse, x_refine, y_refine)
    implicit none

    type(bbox), intent(in)                       :: bbox_coarse, bbox_fine                            !< Bounding boxes of parent and nest
    real(kind=R_GRID), allocatable, intent(in)   :: p_grid(:,:,:), n_grid(:,:,:)                      !< Latlon rids of parent and nest in radians
    real, allocatable, intent(inout)             :: wt(:,:,:)                                         !< Interpolation weight array
    integer, intent(in)                          :: istart_coarse, jstart_coarse, x_refine, y_refine  !< Offsets and nest refinements

    integer       :: i,j, ic, jc
    real          :: dist1, dist2, dist3, dist4, sum
    logical       :: verbose = .false.
    !logical       :: verbose = .true.

    integer       :: this_pe

    real(kind=R_GRID)  :: pi = 4 * atan(1.0d0)
    real               :: pi180
    real               :: rad2deg, deg2rad

    pi180 = pi / 180.0
    deg2rad = pi / 180.0
    rad2deg = 1.0 / pi180

    this_pe = mpp_pe()

    if ( bbox_coarse%is == 0 .and. bbox_coarse%ie == -1 ) then
      ! Skip this one
      if (debug_log) print '("[INFO] WDR skip calc weights npe=",I0)', this_pe

    else
      if (debug_log) print '("[INFO] WDR run calc weights npe=",I0)', this_pe

      ! Calculate the bounding parent grid points for the nest grid point
      ! Rely on the nest being aligned
      ! code is from $CUBE/tools/fv_grid_tools.F90
      !

      do j = bbox_fine%js, bbox_fine%je
        ! F90 integer division truncates
        jc = jstart_coarse  + (j + y_refine/2 + 1) / y_refine
        do i = bbox_fine%is, bbox_fine%ie
          ic = istart_coarse  + (i + x_refine/2 + 1) / x_refine

          if (verbose) then
            if (debug_log) print '("[INFO] WDR MAP npe=",I0," istart_coarse, jstart_coarse,   ic,if,jc,jf",I3,I3," ",I3,I3,I3,I3)', this_pe, istart_coarse, jstart_coarse,ic,i,jc,j

            if (debug_log) print '("[INFO] WDR LATLON npe=",I0,"  p_grid(",I3,I3,")",F8.2,F8.2, F8.2)', this_pe, ic, jc, rad2deg*p_grid(ic,jc,1)-360.0 , rad2deg*p_grid(ic,jc,2), rad2deg*p_grid(ic,jc,1)
            if (debug_log) print '("[INFO] WDR LATLON npe=",I0,"  nest n_grid(",I3,I3,") ",F8.2,F8.2, F8.2)', this_pe, i, j, rad2deg*n_grid(i,j,1)-360.0, rad2deg*n_grid(i,j,2), rad2deg*n_grid(i,j,1)

            if (debug_log) print '("[INFO] WDR LOC npe=",I0,"  -------------------")', this_pe
            if (debug_log) print '("[INFO] WDR LOC npe=",I0,"  A p_grid(",I3,I3,")",F12.6,F12.6, F12.6)', this_pe, ic, jc, rad2deg*p_grid(ic,jc,1)-360.0, rad2deg*p_grid(ic,jc,2), rad2deg*p_grid(ic,jc,1)
            if (debug_log) print '("[INFO] WDR LOC npe=",I0,"  B p_grid(",I3,I3,")",F12.6,F12.6, F12.6)', this_pe, ic, jc+1, rad2deg*p_grid(ic,jc+1,1)-360.0, rad2deg*p_grid(ic,jc+1,2), rad2deg*p_grid(ic,jc+1,1)
            if (debug_log) print '("[INFO] WDR LOC npe=",I0,"  C p_grid(",I3,I3,")",F12.6,F12.6, F12.6)', this_pe, ic+1, jc+1, rad2deg*p_grid(ic+1,jc+1,1)-360.0, rad2deg*p_grid(ic+1,jc+1,2), rad2deg*p_grid(ic+1,jc+1,1)
            if (debug_log) print '("[INFO] WDR LOC npe=",I0,"  D p_grid(",I3,I3,")",F12.6,F12.6, F12.6)', this_pe, ic+1, jc, rad2deg*p_grid(ic+1,jc,1)-360.0, rad2deg*p_grid(ic+1,jc,2), rad2deg*p_grid(ic+1,jc,1)
            if (debug_log) print '("[INFO] WDR LOC npe=",I0,"  nest n_grid(",I3,I3,") ",F12.6,F12.6, F12.6)', this_pe, i, j, rad2deg*n_grid(i,j,1)-360.0, rad2deg*n_grid(i,j,2), rad2deg*n_grid(i,j,1)
          endif

          ! dist2side_latlon takes points in longitude-latitude coordinates.
          dist1 = dist2side_latlon(p_grid(ic,jc,:),     p_grid(ic,jc+1,:),   n_grid(i,j,:))

          if (verbose) then
            if (debug_log) print '("[INFO] WDR LATLON npe=",I0," dist1=",F9.4," p_grid(",I3,I3,")=",F9.4,F9.4," p_grid(",I3,I3,")=",F9.4,F9.4," n_grid(",I3,I3,")=",F9.4,F9.4)', this_pe, dist1, ic, jc, rad2deg*p_grid(ic,jc,1)-360.0,  rad2deg*p_grid(ic,jc,2),  ic, jc+1, rad2deg*p_grid(ic,jc+1,1)-360.0,  rad2deg*p_grid(ic,jc+1,2), i, j, rad2deg*n_grid(i,j,1)-360.0,  rad2deg*n_grid(i,j,2)
          endif
          dist2 = dist2side_latlon(p_grid(ic,jc+1,:),   p_grid(ic+1,jc+1,:), n_grid(i,j,:))
          dist3 = dist2side_latlon(p_grid(ic+1,jc+1,:), p_grid(ic+1,jc,:),   n_grid(i,j,:))
          dist4 = dist2side_latlon(p_grid(ic,jc,:),     p_grid(ic+1,jc,:),   n_grid(i,j,:))

          !if (debug_log) print '("[INFO] WDR LATLON npe=",I0,"  dists at (",I3,I3,"): dist: ",F12.4, F12.4, F12.4, F12.4)', this_pe, i, j, dist1*RADIUS, dist2*RADIUS, dist3*RADIUS, dist4*RADIUS
          if (verbose) then
            if (debug_log) print '("[INFO] WDR LATLON npe=",I0,"  dists at (",I3,I3,"): dist: ",F12.4, F12.4, F12.4, F12.4)', this_pe, i, j, dist1, dist2, dist3, dist4
          endif

          wt(i,j,1)=dist2*dist3      ! ic,   jc    weight
          wt(i,j,2)=dist3*dist4      ! ic,   jc+1  weight
          wt(i,j,3)=dist4*dist1      ! ic+1, jc+1  weight
          wt(i,j,4)=dist1*dist2      ! ic+1, jc    weight

          sum=wt(i,j,1)+wt(i,j,2)+wt(i,j,3)+wt(i,j,4)
          wt(i,j,:)=wt(i,j,:)/sum

          if (verbose) then
            if (debug_log) print '("[INFO] WDR LATLON npe=",I0,"  sum (",I3,I3,"): ",F12.2,"  wt: ",F12.6, F12.6, F12.6, F12.6)', this_pe, i, j, sum, wt(i,j,1), wt(i,j,2), wt(i,j,3), wt(i,j,4)
          endif

        enddo
      enddo
    endif

    if (debug_log) print '("[INFO] WDR DONE calc weights npe=",I0)', this_pe

  end subroutine calc_nest_halo_weights

#endif ! MOVING_NEST

end module fv_moving_nest_mod

