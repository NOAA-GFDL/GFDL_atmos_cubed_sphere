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
!! @brief Provides Moving Nest functionality for physics and surface variables
!! @author W. Ramstrom.  Collaboration with Bin Liu and Chunxi Zhang, EMC
!! @email William.Ramstrom@noaa.gov
! =======================================================================!


! =======================================================================!
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

module fv_moving_nest_physics_mod
#ifdef MOVING_NEST

  use block_control_mod,      only: block_control_type
  use fms_mod,                only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE, clock_flag_default
  use mpp_mod,                only: mpp_pe, mpp_sync, mpp_sync_self, mpp_send, mpp_error, NOTE, FATAL
  use mpp_domains_mod,        only: mpp_update_domains, mpp_get_data_domain, mpp_get_global_domain
  use mpp_domains_mod,        only: mpp_define_nest_domains, mpp_shift_nest_domains, nest_domain_type, domain2d
  use mpp_domains_mod,        only: mpp_get_C2F_index, mpp_update_nest_fine
  use mpp_domains_mod,        only: mpp_get_F2C_index, mpp_update_nest_coarse
  use mpp_domains_mod,        only: NORTH, SOUTH, EAST, WEST, CORNER, CENTER
  use mpp_domains_mod,        only: NUPDATE, SUPDATE, EUPDATE, WUPDATE, DGRID_NE

#ifdef GFS_TYPES
  use GFS_typedefs,           only: IPD_data_type => GFS_data_type, &
      IPD_control_type => GFS_control_type, kind_phys
#else
  use IPD_typedefs,           only: IPD_data_type, IPD_control_type, kind_phys => IPD_kind_phys
#endif
  use GFS_init,               only: GFS_grid_populate

  use boundary_mod,           only: update_coarse_grid, update_coarse_grid_mpp
  use bounding_box_mod,       only: bbox, bbox_get_C2F_index, fill_bbox, show_bbox
  use constants_mod,          only: cp_air, rdgas, grav, rvgas, kappa, pstd_mks, hlv
  use field_manager_mod,      only: MODEL_ATMOS
  use fms_io_mod,             only: read_data, write_data, get_global_att_value, fms_io_init, fms_io_exit
  use fv_arrays_mod,          only: fv_atmos_type, fv_nest_type, fv_grid_type, R_GRID
  use fv_moving_nest_types_mod,   only: fv_moving_nest_prog_type, fv_moving_nest_physics_type, mn_surface_grids, fv_moving_nest_type
  use fv_arrays_mod,          only: allocate_fv_nest_bc_type, deallocate_fv_nest_bc_type
  use fv_grid_tools_mod,      only: init_grid
  use fv_grid_utils_mod,      only: grid_utils_init, ptop_min, dist2side_latlon
  use fv_mapz_mod,            only: Lagrangian_to_Eulerian, moist_cv, compute_total_energy
  use fv_moving_nest_utils_mod, only: check_array, check_local_array, show_atm, show_atm_grids, show_nest_grid, show_tile_geo, grid_equal
  use fv_nesting_mod,         only: dealloc_nested_buffers
  use fv_nwp_nudge_mod,       only: do_adiabatic_init
  use init_hydro_mod,         only: p_var
  use tracer_manager_mod,     only: get_tracer_index, get_tracer_names
  use fv_moving_nest_utils_mod,  only: alloc_halo_buffer, load_nest_latlons_from_nc, grid_geometry, output_grid_to_nc, find_nest_alignment
  use fv_moving_nest_utils_mod,  only: fill_nest_from_buffer, fill_nest_from_buffer_cell_center, fill_nest_from_buffer_nearest_neighbor
  use fv_moving_nest_utils_mod,  only: fill_nest_halos_from_parent, fill_grid_from_supergrid, fill_weight_grid
  use fv_moving_nest_utils_mod,  only: alloc_read_data
  use fv_moving_nest_utils_mod,  only: fill_nest_from_buffer_cell_center_masked
  use fv_moving_nest_utils_mod,  only: fill_nest_halos_from_parent_masked

  use fv_moving_nest_mod,     only: mn_var_fill_intern_nest_halos, mn_var_dump_to_netcdf, mn_var_shift_data
  use fv_moving_nest_types_mod, only: Moving_nest
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
  logical :: move_physics = .true.       ! Always true, unless developer sets move_physics to .False. here for debugging.
  logical :: move_nsst = .true.          ! Value is reset in fv_moving_nest_main.F90 from namelist options

  ! Persistent variables to enable debug printing after range warnings.
  type (fv_atmos_type), pointer                 :: save_Atm_n
  type (block_control_type), pointer            :: save_Atm_block
  type(IPD_control_type), pointer               :: save_IPD_Control
  type(IPD_data_type), pointer                  :: save_IPD_Data(:)

#include <fms_platform.h>

contains

  !>@brief The subroutine 'mn_phys_reset_sfc_props' sets the static surface parameters from the high-resolution input file data
  !>@details This subroutine relies on earlier code reading the data from files into the mn_static data structure
  !!  This subroutine does not yet handle ice points or frac_grid - fractional landfrac/oceanfrac values
  subroutine mn_phys_reset_sfc_props(Atm, n, mn_static, Atm_block, IPD_data, ioffset, joffset, refine)
    type(fv_atmos_type), intent(inout),allocatable   :: Atm(:)              !< Array of atmospheric data
    integer, intent(in)                              :: n                   !< Current grid number
    type(mn_surface_grids), intent(in)               :: mn_static           !< Static surface data
    type(block_control_type), intent(in)             :: Atm_block           !< Physics block layout
    type(IPD_data_type), intent(inout)               :: IPD_data(:)         !< Physics variable data
    integer, intent(in)                              :: ioffset, joffset    !< Current nest offset in i,j direction
    integer, intent(in)                              :: refine              !< Nest refinement ratio

    ! For iterating through physics/surface vector data
    integer                 :: nb, blen, ix, i_pe, j_pe, i_idx, j_idx
    real(kind=kind_phys)    :: phys_oro

    ! Setup local land sea mask grid for masked interpolations
    do i_pe = Atm(n)%bd%isd, Atm(n)%bd%ied
      do j_pe = Atm(n)%bd%jsd, Atm(n)%bd%jed
        i_idx = (ioffset-1)*refine + i_pe
        j_idx = (joffset-1)*refine + j_pe

        Moving_nest(n)%mn_phys%slmsk(i_pe, j_pe) = mn_static%ls_mask_grid(i_idx, j_idx)
      enddo
    enddo

    !  Reset the variables from the fix_sfc files
    do nb = 1,Atm_block%nblks
      blen = Atm_block%blksz(nb)
      do ix = 1, blen
        i_pe = Atm_block%index(nb)%ii(ix)
        j_pe = Atm_block%index(nb)%jj(ix)

        i_idx = (ioffset-1)*refine + i_pe
        j_idx = (joffset-1)*refine + j_pe

        ! Reset the land sea mask from the hires parent data
        IPD_data(nb)%Sfcprop%slmsk(ix) = mn_static%ls_mask_grid(i_idx, j_idx)

        !  IFD values are 0 for land, and 1 for oceans/lakes -- reverse of the land sea mask
        !  Land Sea Mask has values of 0 for oceans/lakes, 1 for land, 2 for sea ice
        !  TODO figure out what ifd should be for sea ice
        if (mn_static%ls_mask_grid(i_idx, j_idx) .eq. 1 ) then
          if (move_nsst) IPD_data(nb)%Sfcprop%ifd(ix) = 0         ! Land
          IPD_data(nb)%Sfcprop%oceanfrac(ix) = 0   ! Land -- TODO permit fractions
          IPD_data(nb)%Sfcprop%landfrac(ix) = 1    ! Land -- TODO permit fractions
        else
          if (move_nsst) IPD_data(nb)%Sfcprop%ifd(ix) = 1         ! Ocean
          IPD_data(nb)%Sfcprop%oceanfrac(ix) = 1   ! Ocean -- TODO permit fractions
          IPD_data(nb)%Sfcprop%landfrac(ix) = 0    ! Ocean -- TODO permit fractions
        endif

        IPD_data(nb)%Sfcprop%tg3(ix) = mn_static%deep_soil_temp_grid(i_idx, j_idx)

        ! Follow logic from FV3/io/FV3GFS_io.F90 line 1187
        ! TODO this will need to be more complicated if we support frac_grid
        !if (nint(mn_static%soil_type_grid(i_idx, j_idx)) == 14 .or. int(mn_static%soil_type_grid(i_idx, j_idx)+0.5) <= 0) then
        !if (nint(mn_static%soil_type_grid(i_idx, j_idx)) == 14 .or.

        !if ( (mn_static%ls_mask_grid(i_idx, j_idx) .eq. 1 .and. nint(mn_static%land_frac_grid(i_idx, j_idx)) == 0) .or. &
        !    mn_static%soil_type_grid(i_idx, j_idx) < 0.5) then
        if (mn_static%ls_mask_grid(i_idx, j_idx) .eq. 1 .and. nint(mn_static%land_frac_grid(i_idx, j_idx)) == 0 ) then
          ! Water soil type == lake, etc. -- override the other variables and make this water
          print '("WDR mn_phys_reset_sfc_props LAKE SOIL npe=",I0," x,y=",I0,",",I0," lat=",F10.3," lon=",F10.3)', mpp_pe(), i_idx, j_idx, IPD_data(nb)%Grid%xlat_d(ix), IPD_data(nb)%Grid%xlon_d(ix)-360.0

          if (move_nsst) IPD_data(nb)%Sfcprop%ifd(ix) = 1         ! Ocean
          IPD_data(nb)%Sfcprop%oceanfrac(ix) = 1   ! Ocean -- TODO permit fractions
          IPD_data(nb)%Sfcprop%landfrac(ix) = 0    ! Ocean -- TODO permit fractions

          IPD_data(nb)%Sfcprop%stype(ix) = 0
          IPD_data(nb)%Sfcprop%slmsk(ix) = 0
        else
          IPD_data(nb)%Sfcprop%stype(ix) = nint(mn_static%soil_type_grid(i_idx, j_idx))
        endif

        !IPD_data(nb)%Sfcprop%vfrac(ix) = mn_static%veg_frac_grid(i_idx, j_idx)
        IPD_data(nb)%Sfcprop%vtype(ix) = nint(mn_static%veg_type_grid(i_idx, j_idx))
        IPD_data(nb)%Sfcprop%slope(ix) = nint(mn_static%slope_type_grid(i_idx, j_idx))
        IPD_data(nb)%Sfcprop%snoalb(ix) = mn_static%max_snow_alb_grid(i_idx, j_idx)

        IPD_data(nb)%Sfcprop%facsf(ix) = mn_static%facsf_grid(i_idx, j_idx)
        IPD_data(nb)%Sfcprop%facwf(ix) = mn_static%facwf_grid(i_idx, j_idx)

        IPD_data(nb)%Sfcprop%alvsf(ix) = mn_static%alvsf_grid(i_idx, j_idx)
        IPD_data(nb)%Sfcprop%alvwf(ix) = mn_static%alvwf_grid(i_idx, j_idx)
        IPD_data(nb)%Sfcprop%alnsf(ix) = mn_static%alnsf_grid(i_idx, j_idx)
        IPD_data(nb)%Sfcprop%alnwf(ix) = mn_static%alnwf_grid(i_idx, j_idx)

        ! Reset the orography in the physics arrays, using the smoothed values from above
        phys_oro =  Atm(n)%phis(i_pe, j_pe) / grav
        IPD_data(nb)%Sfcprop%oro(ix) = phys_oro
        IPD_data(nb)%Sfcprop%oro_uf(ix) = phys_oro

      enddo
    enddo

  end subroutine mn_phys_reset_sfc_props

  !>@brief The subroutine 'mn_phys_reset_phys_latlon' sets the lat/lons from the high-resolution input file data
  !>@details This subroutine sets lat/lons of the moved nest, then recalculates all the derived quantities (dx,dy,etc.)
  subroutine mn_reset_phys_latlon(Atm, n, tile_geo, fp_super_tile_geo, Atm_block, IPD_control, IPD_data)
    type(fv_atmos_type), intent(in)      :: Atm(:)               !< Array of atmospheric data
    integer, intent(in)                  :: n                    !< Current grid number
    type(grid_geometry), intent(in)      :: tile_geo             !< Bounds of this grid
    type(grid_geometry), intent(in)      :: fp_super_tile_geo    !< Bounds of high-resolution parent grid
    type(block_control_type), intent(in) :: Atm_block            !< Physics block layout
    type(IPD_control_type), intent(in)   :: IPD_control          !< Physics metadata
    type(IPD_data_type), intent(inout)   :: IPD_data(:)          !< Physics variable data

    integer :: isc, jsc, iec, jec
    integer :: x, y, fp_i, fp_j
    integer :: nest_x, nest_y, parent_x, parent_y
    integer :: this_pe

    real(kind=kind_phys), allocatable :: lats(:,:), lons(:,:), area(:,:)

    this_pe = mpp_pe()

    isc = Atm(n)%bd%isc
    jsc = Atm(n)%bd%jsc
    iec = Atm(n)%bd%iec
    jec = Atm(n)%bd%jec

    allocate(lats(isc:iec, jsc:jec))
    allocate(lons(isc:iec, jsc:jec))
    allocate(area(isc:iec, jsc:jec))

    ! This is going to be slow -- replace with better way to calculate index offsets, or pass them from earlier calculations
    ! TODO optimization here
    call find_nest_alignment(tile_geo, fp_super_tile_geo, nest_x, nest_y, parent_x, parent_y)
    !print '("WDR mn_reset_phys_latlon AB npe=",I0)', this_pe

    do x = isc, iec
      do y = jsc, jec
        fp_i = (x - nest_x) * 2 + parent_x
        fp_j = (y - nest_y) * 2 + parent_y

        lons(x,y) = fp_super_tile_geo%lons(fp_i, fp_j)
        lats(x,y) = fp_super_tile_geo%lats(fp_i, fp_j)

        ! Need to add the areas from 4 squares, because the netCDF file has areas calculated for the supergrid cells
        !  We need the area of the whole center of the cell.
        !  Example dimensions for C288_grid.tile6.nc
        !   longitude -- x(577,577)
        !   latitude  -- y(577,577)
        !   area      -- x(576,576)

        !  Extracting lat/lon/area from Supergrid
        !
        !   1,1----2,1----3,1
        !    |      |      |
        !    | a1,1 | a2,1 |
        !    |      |      |
        !   1,2----2,2----3,2
        !    |      |      |
        !    | a1,2 | a2,2 |
        !    |      |      |
        !   1,3----2,3----3,3
        !
        !  The model A-grid cell 1,1 is centered at supergrid location 2,2
        !    The area of the A-grid cell is the sum of the 4 supergrid areas   A = a(1,1) + a(1,2) + a(2,1) + a(2,2)

        area(x,y) = fp_super_tile_geo%area(fp_i - 1, fp_j - 1) + fp_super_tile_geo%area(fp_i - 1, fp_j) + &
            fp_super_tile_geo%area(fp_i, fp_j - 1) + fp_super_tile_geo%area(fp_i, fp_j)   ! TODO make sure these offsets are correct.
      enddo
    enddo

    call GFS_grid_populate(IPD_data%Grid, lons, lats, area)

    deallocate(lats)
    deallocate(lons)
    deallocate(area)

  end subroutine mn_reset_phys_latlon

  !>@brief The subroutine 'dump_surface_physics' outputs surface physics data for a given point and its neighbors to stdout
  !>@details This subroutine is appropriate to be called for debugging when range warnings are detected, in tools/fv_diagnostics.F90.
  subroutine dump_surface_physics(i_out, j_out, k_out)
    integer, intent(in)    :: i_out, j_out, k_out     !< i,j,k values of point to output

    integer :: nb, blen, ix, i, j, k, kk
    integer :: this_pe

    this_pe = mpp_pe()

    if (associated(save_Atm_block)) then
      print '("WDR dump_surface_physics npe=",I0)', this_pe
    else
      print '("WDR dump_surface_physics RANGE RETURN npe=",I0)', this_pe
      return
    end if

    k = k_out

    do nb = 1,save_Atm_block%nblks
      blen = save_Atm_block%blksz(nb)
      do ix = 1, blen
        ! Get the indices only once, before iterating through vertical levels or number of variables
        !  Was there a different efficiency from having the k loop outside?
        i = save_Atm_block%index(nb)%ii(ix)
        j = save_Atm_block%index(nb)%jj(ix)

        if (i .ge. i_out-2 .and. i .le. i_out+2 .and. j .ge. j_out-2 .and. j .le. j_out+2) then

          print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") slmsk=",F8.5, " lakefrac=",F10.5, " lakedepth=",F14.5, " landfrac=",F10.5, " oro=",F10.5, " oro_uf=",F10.5)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%slmsk(ix), save_IPD_Data(nb)%Sfcprop%lakefrac(ix), save_IPD_Data(nb)%Sfcprop%lakedepth(ix), save_IPD_Data(nb)%Sfcprop%landfrac(ix), save_IPD_Data(nb)%Sfcprop%oro(ix), save_IPD_Data(nb)%Sfcprop%oro_uf(ix)
          print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") oro=",F10.5, " oro_uf=",F10.5, " phis/g=",F10.5, " slope=",I0)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%oro(ix), save_IPD_Data(nb)%Sfcprop%oro_uf(ix), save_Atm_n%phis(i,j)/grav, save_IPD_Data(nb)%Sfcprop%slope(ix)
          print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") emis_lnd=",F10.4," emis_ice=",F10.4," emis_wat=",F10.4," hflx=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%emis_lnd(ix), save_IPD_Data(nb)%Sfcprop%emis_ice(ix), save_IPD_Data(nb)%Sfcprop%emis_wat(ix), save_IPD_Data(nb)%Sfcprop%hflx(ix)
          print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") albdirvis_lnd=",F10.4," albdirnir_lnd=",F10.4," albdifvis_lnd=",F10.4," albdifnir_lnd=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%albdirvis_lnd(ix), save_IPD_Data(nb)%Sfcprop%albdirnir_lnd(ix), save_IPD_Data(nb)%Sfcprop%albdifvis_lnd(ix), save_IPD_Data(nb)%Sfcprop%albdifnir_lnd(ix)
          !print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") albdirvis_ice=",F10.4," albdirnir_ice=",F10.4," albdifvis_ice=",F10.4," albdifnir_ice=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%albdirvis_ice(ix), save_IPD_Data(nb)%Sfcprop%albdirnir_ice(ix), save_IPD_Data(nb)%Sfcprop%albdifvis_ice(ix), save_IPD_Data(nb)%Sfcprop%albdifnir_ice(ix)
          if (associated(save_IPD_Data(nb)%Sfcprop%qss)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") qss=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%qss(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%evap)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") evap=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%evap(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%sncovr)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") sncovr",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%sncovr(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%sncovr_ice)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") sncovr_ice",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%sncovr_ice(ix)
          endif
          if (associated(save_IPD_Data(nb)%Intdiag%total_albedo)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") Intdiag%total_albedo=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Intdiag%total_albedo(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%ifd)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") ifd=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%ifd(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%semisbase)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") semisbase=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%semisbase(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%sfalb_lnd)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") sfalb_lnd=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%sfalb_lnd(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%sfalb_ice)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") sfalb_ice=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%sfalb_ice(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%sfalb_lnd_bck)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") sfalb_lnd_bck=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%sfalb_lnd_bck(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%emis_lnd)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") emis_lnd=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%emis_lnd(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%emis_ice)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") emis_ice=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%emis_ice(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%emis_wat)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") emis_wat=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%emis_wat(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%tvxy)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") veg temp tvxy=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%tvxy(ix)
          endif

          if (associated(save_IPD_Data(nb)%Sfcprop%tgxy)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") ground temp tgxy=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%tgxy(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%tg3)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") deep soil temp tg3=",F10.4," slmsk=",F8.3)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%tg3(ix), save_IPD_Data(nb)%Sfcprop%slmsk(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%alboldxy)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") alboldxy=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%alboldxy(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%shdmin)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") shdmin=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%shdmin(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%shdmax)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") shdmax=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%shdmax(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%stype)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") stype=",I0)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%stype(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%vtype)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") vtype=",I0)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%vtype(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%stype_save)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") stype_save=",I0)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%stype_save(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%vtype_save)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") vtype_save=",I0)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%vtype_save(ix)
          endif
          do kk = 1, save_IPD_Control%nmtvr
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") hprime(",I0,")=",F10.4)', this_pe, i, j, kk, save_IPD_Data(nb)%Sfcprop%hprime(ix,kk)
          enddo
          if (associated(save_IPD_Data(nb)%Sfcprop%snowd)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") snowd=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%snowd(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%weasd)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") weasd=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%weasd(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%ffmm)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") ffmm=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%ffmm(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%ffhh)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") ffhh=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%ffhh(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%f10m)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") f10m=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%f10m(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%uustar)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") uustar=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%uustar(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%z0base)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") z0base=",F18.6)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%z0base(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%zorl)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") zorl=",F18.6," ",E15.6)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%zorl(ix), save_IPD_Data(nb)%Sfcprop%zorl(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%zorlw)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") zorlw=",F18.6," ",E15.6)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%zorlw(ix), save_IPD_Data(nb)%Sfcprop%zorlw(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%zorll)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") zorll=",F18.6," ",E15.6)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%zorll(ix), save_IPD_Data(nb)%Sfcprop%zorll(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%zorli)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") zorli=",F15.6," ",E15.6)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%zorli(ix), save_IPD_Data(nb)%Sfcprop%zorli(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%zorlwav)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") zorlwav=",F15.6," ",E15.6)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%zorlwav(ix), save_IPD_Data(nb)%Sfcprop%zorlwav(ix)
          endif
          if (associated(save_IPD_Data(nb)%Coupling%tsfc_radtime)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") Coupling%tsfc_radtime=",F15.6)', this_pe, i, j, save_IPD_Data(nb)%Coupling%tsfc_radtime(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%canopy)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") canopy=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%canopy(ix)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%vfrac)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") vfrac=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%vfrac(ix)
          endif
          if (associated(save_IPD_Data(nb)%Radtend%sfalb)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") Radtend%sfalb=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Radtend%sfalb(ix)
          endif
          if (associated(save_IPD_Data(nb)%Radtend%semis)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") Radtend%semis=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Radtend%semis(ix)
          endif
          if (associated(save_IPD_Data(nb)%Radtend%sfcfsw)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") Radtend%sfcfsw upfxc=",F10.4,"  upfx0=",F10.4," dnfxc=",F10.4," dnfx0=",F10.4)', &
                this_pe, i, j, save_IPD_Data(nb)%Radtend%sfcfsw(ix)%upfxc, save_IPD_Data(nb)%Radtend%sfcfsw(ix)%upfx0, save_IPD_Data(nb)%Radtend%sfcfsw(ix)%dnfxc, save_IPD_Data(nb)%Radtend%sfcfsw(ix)%dnfx0
          endif
          if (associated(save_IPD_Data(nb)%Radtend%sfcflw)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") Radtend%sfcflw  upfxc=",F10.4," upfx0=",F10.4," dnfxc=",F10.4," dnfx0=",F10.4)', &
                this_pe, i, j, save_IPD_Data(nb)%Radtend%sfcflw(ix)%upfxc, save_IPD_Data(nb)%Radtend%sfcflw(ix)%upfx0, save_IPD_Data(nb)%Radtend%sfcflw(ix)%dnfxc, save_IPD_Data(nb)%Radtend%sfcflw(ix)%dnfx0
          endif
          if (associated(save_IPD_Data(nb)%Radtend%coszen)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") Radtend%coszen=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Radtend%coszen(ix)
          endif
          if (associated(save_IPD_Data(nb)%Radtend%coszdg)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") Radtend%coszdg=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Radtend%coszdg(ix)
          endif

          !if (associated(save_IPD_Data(nb)%Sfcprop%semisbase)) then
          !  print '("[INFO] WDR RANGEP AA this_pe= ",I0)', this_pe
          !  !if (associated (save_IPD_Data(nb)%Sfcprop%sfalb_lnd)) then
          !  print '("[INFO] WDR RANGEP AB this_pe= ",I0)', this_pe
          !  !if (associated(save_IPD_Data(nb)%Sfcprop%sfalb_ice)) then
          !  print '("[INFO] WDR RANGEP AC this_pe= ",I0)', this_pe
          !  if (associated(save_IPD_Data(nb)%Sfcprop%sfalb_lnd_bck)) then
          !    !print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") semisbase=",F10.4," sfalb_lnd=",F10.4," sfalb_ice=",F10.4," sfalb_lnd_bck=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%semisbase(ix), save_IPD_Data(nb)%Sfcprop%sfalb_lnd(ix), save_IPD_Data(nb)%Sfcprop%sfalb_ice(ix), save_IPD_Data(nb)%Sfcprop%sfalb_lnd_bck(ix)
          !    print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") semisbase=",F10.4," sfalb_lnd_bck=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%semisbase(ix), save_IPD_Data(nb)%Sfcprop%sfalb_lnd_bck(ix)
          !  endif
          !  !endif
          !  !endif
          !endif

          if (associated(save_IPD_Data(nb)%Sfcprop%alvsf)) then
            !print '("[INFO] WDR RANGEP BA this_pe= ",I0)', this_pe
            if (associated(save_IPD_Data(nb)%Sfcprop%alnsf)) then
              !print '("[INFO] WDR RANGEP BB this_pe= ",I0)', this_pe
              if (associated(save_IPD_Data(nb)%Sfcprop%alvwf)) then
                !print '("[INFO] WDR RANGEP BC this_pe= ",I0)', this_pe
                if (associated(save_IPD_Data(nb)%Sfcprop%alnwf)) then
                  print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") alvsf=",F10.4," alnsf=",F10.4," alvwf=",F10.4," alnwf=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%alvsf(ix), save_IPD_Data(nb)%Sfcprop%alnsf(ix), save_IPD_Data(nb)%Sfcprop%alvwf(ix), save_IPD_Data(nb)%Sfcprop%alnwf(ix)
                endif
              endif
            endif
          endif

          print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") sncovr=",F10.4," snoalb=",F10.4," facsf=",F10.4," facwf=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%sncovr(ix), save_IPD_Data(nb)%Sfcprop%snoalb(ix), save_IPD_Data(nb)%Sfcprop%facsf(ix), save_IPD_Data(nb)%Sfcprop%facwf(ix)

          if (associated(save_IPD_Data(nb)%Sfcprop%t2m)) then
            if (associated(save_IPD_Data(nb)%Sfcprop%th2m)) then
              print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") t2m=",F10.4," th2m=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%t2m(ix), save_IPD_Data(nb)%Sfcprop%th2m(ix)
            else
              print '("[INFO] WDR RANGEP CB this_pe= ",I0)', this_pe
            endif
          else
            print '("[INFO] WDR RANGEP CA this_pe= ",I0)', this_pe
          endif

          print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") tsfc=",F10.4," tsfco=",F10.4," tsfcl=",F10.4," tisfc=",F10.4," stc1=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%tsfc(ix), save_IPD_Data(nb)%Sfcprop%tsfco(ix), save_IPD_Data(nb)%Sfcprop%tsfcl(ix), save_IPD_Data(nb)%Sfcprop%tisfc(ix), save_IPD_Data(nb)%Sfcprop%stc(ix,1)
          print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") psurf=",F10.4," t2m=",F10.4," th2m=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%IntDiag%psurf(ix), save_IPD_Data(nb)%Sfcprop%t2m(ix), save_IPD_Data(nb)%Sfcprop%th2m(ix)

          if (associated(save_IPD_Data(nb)%Sfcprop%slc)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") soil moist slc1=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%slc(ix,1)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%smc)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") tot soil moist smc1=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%smc(ix,1)
          endif
          if (associated(save_IPD_Data(nb)%Sfcprop%stc)) then
            print '("[INFO] WDR RANGE this_pe= ",I0," i,j=(",I0,",",I0,") soil temp stc1=",F10.4)', this_pe, i, j, save_IPD_Data(nb)%Sfcprop%stc(ix,1)
          endif

        endif
      enddo
    enddo
  end subroutine dump_surface_physics

  !>@brief The subroutine 'mn_phys_fill_temp_variables' extracts 1D physics data into a 2D array for nest motion
  !>@details This subroutine fills in the mn_phys structure on the Atm object with 2D arrays of physics/surface variables.
  !!  Note that ice variables are not yet handled.
  subroutine mn_phys_fill_temp_variables(Atm, Atm_block, IPD_Control, IPD_Data, n, child_grid_num, is_fine_pe, npz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)            !< Array of atmospheric data
    type (block_control_type), target, intent(in)            :: Atm_block         !< Physics block layout
    type(IPD_control_type), target, intent(in)               :: IPD_Control       !< Physics metadata
    type(IPD_data_type), target, intent(inout)               :: IPD_Data(:)       !< Physics variable data
    integer, intent(in)                                      :: n, child_grid_num !< Current grid number, child grid number
    logical, intent(in)                                      :: is_fine_pe        !< Is this a nest PE?
    integer, intent(in)                                      :: npz               !< Number of vertical levels

    integer :: isd, ied, jsd, jed
    integer :: is, ie, js, je
    integer :: this_pe

    integer :: nb, blen, i, j, k, ix, nv
    type(fv_moving_nest_physics_type), pointer       :: mn_phys

    this_pe = mpp_pe()

    save_Atm_n => Atm(n)
    save_Atm_block => Atm_block
    save_IPD_Control => IPD_Control
    save_IPD_Data => IPD_Data

    if (debug_log) print '("[INFO] WDR start mn_phys_fill_temp_variables. npe=",I0," n=",I0)', this_pe, n

    isd = Atm(n)%bd%isd
    ied = Atm(n)%bd%ied
    jsd = Atm(n)%bd%jsd
    jed = Atm(n)%bd%jed

    !if (is_fine_pe) call dump_surface_physics(isd+8, jsd+8, npz-1)

    if (debug_log) print '("[INFO] WDR mn_phys_fill_temp_variables. npe=",I0," isd=",I0," ied=",I0," jsd=",I0," jed=",I0)', this_pe, isd, ied, jsd, jed

    is = Atm(n)%bd%is
    ie = Atm(n)%bd%ie
    js = Atm(n)%bd%js
    je = Atm(n)%bd%je

    if (debug_log) print '("[INFO] WDR mn_phys_fill_temp_variables. npe=",I0," is=",I0," ie=",I0," js=",I0," je=",I0)', this_pe, is, ie, js, je

    mn_phys => Moving_nest(n)%mn_phys

    mn_phys%ts(is:ie, js:je) =  Atm(n)%ts(is:ie, js:je)

    do nb = 1,Atm_block%nblks
      blen = Atm_block%blksz(nb)
      do ix = 1, blen
        ! Get the indices only once, before iterating through vertical levels or number of variables
        i = Atm_block%index(nb)%ii(ix)
        j = Atm_block%index(nb)%jj(ix)

        if (move_physics) then
          do k = 1, IPD_Control%lsoil
            mn_phys%smc(i,j,k) = IPD_Data(nb)%Sfcprop%smc(ix,k)
            mn_phys%stc(i,j,k) = IPD_Data(nb)%Sfcprop%stc(ix,k)
            mn_phys%slc(i,j,k) = IPD_Data(nb)%Sfcprop%slc(ix,k)
          enddo

          mn_phys%emis_lnd(i,j)      = IPD_Data(nb)%Sfcprop%emis_lnd(ix)
          mn_phys%emis_ice(i,j)      = IPD_Data(nb)%Sfcprop%emis_ice(ix)
          mn_phys%emis_wat(i,j)      = IPD_Data(nb)%Sfcprop%emis_wat(ix)

          !mn_phys%sfalb_lnd(i,j)     = IPD_Data(nb)%Sfcprop%sfalb_lnd(ix)
          !mn_phys%sfalb_lnd_bck(i,j) = IPD_Data(nb)%Sfcprop%sfalb_lnd_bck(ix)
          !mn_phys%semis(i,j)      = IPD_Data(nb)%Radtend%semis(ix)
          !mn_phys%semisbase(i,j)      = IPD_Data(nb)%Sfcprop%semisbase(ix)
          !mn_phys%sfalb(i,j)      = IPD_Data(nb)%Radtend%sfalb(ix)

          mn_phys%albdirvis_lnd(i,j) = IPD_Data(nb)%Sfcprop%albdirvis_lnd(ix)
          mn_phys%albdirnir_lnd(i,j) = IPD_Data(nb)%Sfcprop%albdirnir_lnd(ix)
          mn_phys%albdifvis_lnd(i,j) = IPD_Data(nb)%Sfcprop%albdifvis_lnd(ix)
          mn_phys%albdifnir_lnd(i,j) = IPD_Data(nb)%Sfcprop%albdifnir_lnd(ix)

          mn_phys%u10m(i,j)  = IPD_Data(nb)%IntDiag%u10m(ix)
          mn_phys%v10m(i,j)  = IPD_Data(nb)%IntDiag%v10m(ix)
          mn_phys%tprcp(i,j)  = IPD_Data(nb)%Sfcprop%tprcp(ix)

          do k = 1, IPD_Control%nmtvr
            mn_phys%hprime(i,j,k)  = IPD_Data(nb)%Sfcprop%hprime(ix,k)
          enddo

          mn_phys%lakefrac(i,j) = IPD_Data(nb)%Sfcprop%lakefrac(ix)
          mn_phys%lakedepth(i,j) = IPD_Data(nb)%Sfcprop%lakedepth(ix)

          mn_phys%canopy(i,j) = IPD_Data(nb)%Sfcprop%canopy(ix)
          mn_phys%vegfrac(i,j)= IPD_Data(nb)%Sfcprop%vfrac(ix)
          mn_phys%uustar(i,j) = IPD_Data(nb)%Sfcprop%uustar(ix)
          mn_phys%shdmin(i,j) = IPD_Data(nb)%Sfcprop%shdmin(ix)
          mn_phys%shdmax(i,j) = IPD_Data(nb)%Sfcprop%shdmax(ix)
          mn_phys%zorl(i,j)   = IPD_Data(nb)%Sfcprop%zorl(ix)
          mn_phys%zorll(i,j)  = IPD_Data(nb)%Sfcprop%zorll(ix)
          mn_phys%zorlwav(i,j)= IPD_Data(nb)%Sfcprop%zorlwav(ix)
          mn_phys%zorlw(i,j)  = IPD_Data(nb)%Sfcprop%zorlw(ix)
          mn_phys%tsfco(i,j)  = IPD_Data(nb)%Sfcprop%tsfco(ix)
          mn_phys%tsfcl(i,j)  = IPD_Data(nb)%Sfcprop%tsfcl(ix)
          mn_phys%tsfc(i,j)   = IPD_Data(nb)%Sfcprop%tsfc(ix)

          mn_phys%albdirvis_lnd(i,j)   = IPD_Data(nb)%Sfcprop%albdirvis_lnd(ix)
          mn_phys%albdirnir_lnd(i,j)   = IPD_Data(nb)%Sfcprop%albdirnir_lnd(ix)
          mn_phys%albdifvis_lnd(i,j)   = IPD_Data(nb)%Sfcprop%albdifvis_lnd(ix)
          mn_phys%albdifnir_lnd(i,j)   = IPD_Data(nb)%Sfcprop%albdifnir_lnd(ix)

          do nv = 1, IPD_Control%ntot2d
            mn_phys%phy_f2d(i,j,nv) = IPD_Data(nb)%Tbd%phy_f2d(ix, nv)
          enddo

          do k = 1, IPD_Control%levs
            do nv = 1, IPD_Control%ntot3d
              mn_phys%phy_f3d(i,j,k,nv) = IPD_Data(nb)%Tbd%phy_f3d(ix, k, nv)
            enddo
          enddo

          ! Cloud prop data has x,y dimensions
          mn_phys%cv(i,j)  = IPD_Data(nb)%Cldprop%cv(ix)
          mn_phys%cvt(i,j) = IPD_Data(nb)%Cldprop%cvt(ix)
          mn_phys%cvb(i,j) = IPD_Data(nb)%Cldprop%cvb(ix)
        endif

        if (move_nsst) then
          mn_phys%tref(i,j)   = IPD_Data(nb)%Sfcprop%tref(ix)
          mn_phys%z_c(i,j)    = IPD_Data(nb)%Sfcprop%z_c(ix)
          mn_phys%c_0(i,j)    = IPD_Data(nb)%Sfcprop%c_0(ix)
          mn_phys%c_d(i,j)    = IPD_Data(nb)%Sfcprop%c_d(ix)
          mn_phys%w_0(i,j)    = IPD_Data(nb)%Sfcprop%w_0(ix)
          mn_phys%w_d(i,j)    = IPD_Data(nb)%Sfcprop%w_d(ix)
          mn_phys%xt(i,j)     = IPD_Data(nb)%Sfcprop%xt(ix)
          mn_phys%xs(i,j)     = IPD_Data(nb)%Sfcprop%xs(ix)
          mn_phys%xu(i,j)     = IPD_Data(nb)%Sfcprop%xu(ix)
          mn_phys%xv(i,j)     = IPD_Data(nb)%Sfcprop%xv(ix)
          mn_phys%xz(i,j)     = IPD_Data(nb)%Sfcprop%xz(ix)
          mn_phys%zm(i,j)     = IPD_Data(nb)%Sfcprop%zm(ix)
          mn_phys%xtts(i,j)   = IPD_Data(nb)%Sfcprop%xtts(ix)
          mn_phys%xzts(i,j)   = IPD_Data(nb)%Sfcprop%xzts(ix)
          mn_phys%d_conv(i,j) = IPD_Data(nb)%Sfcprop%d_conv(ix)
          mn_phys%dt_cool(i,j)= IPD_Data(nb)%Sfcprop%dt_cool(ix)
          mn_phys%qrain(i,j)  = IPD_Data(nb)%Sfcprop%qrain(ix)
        endif
      enddo
    enddo

    if (debug_log) print '("[INFO] WDR end mn_phys_fill_temp_variables. npe=",I0," n=",I0)', this_pe, n

  end subroutine mn_phys_fill_temp_variables

  !>@brief The subroutine 'mn_phys_apply_temp_variables' copies moved 2D data back into 1D physics arryas for nest motion
  !>@details This subroutine fills the 1D physics arrays from the mn_phys structure on the Atm object
  !!  Note that ice variables are not yet handled.
  subroutine mn_phys_apply_temp_variables(Atm, Atm_block, IPD_Control, IPD_Data, n, child_grid_num, is_fine_pe, npz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)            !< Array of atmospheric data
    type (block_control_type), intent(in)                    :: Atm_block         !< Physics block layout
    type(IPD_control_type), intent(in)                       :: IPD_Control       !< Physics metadata
    type(IPD_data_type), intent(inout)                       :: IPD_Data(:)       !< Physics variable data
    integer, intent(in)                                      :: n, child_grid_num !< Current grid number, child grid number
    logical, intent(in)                                      :: is_fine_pe        !< Is this a nest PE?
    integer, intent(in)                                      :: npz               !< Number of vertical levels

    integer :: is, ie, js, je
    integer :: this_pe
    integer :: nb, blen, i, j ,k, ix, nv
    integer :: bad_values, good_values
    type(fv_moving_nest_physics_type), pointer       :: mn_phys

    this_pe = mpp_pe()
    mn_phys => Moving_nest(n)%mn_phys

    if (debug_log) print '("[INFO] WDR start mn_phys_apply_temp_variables. npe=",I0," n=",I0)', this_pe, n

    ! Check if the variables were filled in properly.

    if (debug_log) then
      good_values = 0
      bad_values = 0

      if (is_fine_pe) then
        do i = Atm(n)%bd%isd, Atm(n)%bd%ied
          do j = Atm(n)%bd%jsd, Atm(n)%bd%jed
            if (mn_phys%ts(i,j) .gt. 20000.0) then
              print '("[WARN] WDR BAD NEST ts value. npe=",I0," ts(",I0,",",I0,")=",F12.3)', this_pe, i, j, mn_phys%ts(i,j)
              bad_values = bad_values + 1
            else
              good_values = good_values + 1
            endif
          enddo
        enddo
      else
        do i = Atm(n)%bd%is, Atm(n)%bd%ie
          do j = Atm(n)%bd%js, Atm(n)%bd%je
            if (mn_phys%ts(i,j) .gt. 20000.0) then
              print '("[WARN] WDR BAD GLOBAL ts value. npe=",I0," ts(",I0,",",I0")=",F12.3)', this_pe, i, j, mn_phys%ts(i,j)
              bad_values = bad_values + 1
            else
              good_values = good_values + 1
            endif
          enddo
        enddo
      endif

      i = Atm(n)%bd%is
      j = Atm(n)%bd%js

      print '("[WARN] WDR Surface ts value. npe=",I0," ts(",I0,",",I0,")=",F18.3)', this_pe, i, j, mn_phys%ts(i,j)

      print '("INFO] WDR ts values. npe=",I0," good_values=",I0," bad_values=",I0)', this_pe, good_values, bad_values
    endif

    !  Needed to fill the local grids for parent and nest PEs in order to transmit/interpolate data from parent to nest
    !  But only the nest PE's have changed the values with nest motion, so they are the only ones that need to update the original arrays
    if (is_fine_pe) then
      is = Atm(n)%bd%is
      ie = Atm(n)%bd%ie
      js = Atm(n)%bd%js
      je = Atm(n)%bd%je

      if (debug_log) print '("[INFO] WDR mn_phys_apply_temp_variables. npe=",I0," is=",I0," ie=",I0," js=",I0," je=",I0)', this_pe, is, ie, js, je

      ! SST directly in Atm structure
      Atm(n)%ts(is:ie, js:je) =  mn_phys%ts(is:ie, js:je)

      do nb = 1,Atm_block%nblks
        blen = Atm_block%blksz(nb)
        do ix = 1, blen
          i = Atm_block%index(nb)%ii(ix)
          j = Atm_block%index(nb)%jj(ix)

          if (move_physics) then
            ! Surface properties
            do k = 1, IPD_Control%lsoil
              IPD_Data(nb)%Sfcprop%smc(ix,k) = mn_phys%smc(i,j,k)
              IPD_Data(nb)%Sfcprop%stc(ix,k) = mn_phys%stc(i,j,k)
              IPD_Data(nb)%Sfcprop%slc(ix,k) = mn_phys%slc(i,j,k)
            enddo

            ! WDR EMIS PATCH - Force to positive at all locations.
            if (mn_phys%emis_lnd(i,j) .ge. 0.0) then
              IPD_Data(nb)%Sfcprop%emis_lnd(ix) = mn_phys%emis_lnd(i,j)
            else
              IPD_Data(nb)%Sfcprop%emis_lnd(ix) = 0.5
            endif
            if (mn_phys%emis_ice(i,j) .ge. 0.0) then
              IPD_Data(nb)%Sfcprop%emis_ice(ix) = mn_phys%emis_ice(i,j)
            else
              IPD_Data(nb)%Sfcprop%emis_ice(ix) = 0.5
            endif
            if (mn_phys%emis_wat(i,j) .ge. 0.0) then
              IPD_Data(nb)%Sfcprop%emis_wat(ix) = mn_phys%emis_wat(i,j)
            else
              IPD_Data(nb)%Sfcprop%emis_wat(ix) = 0.5
            endif

            !IPD_Data(nb)%Sfcprop%sfalb_lnd(ix) = mn_phys%sfalb_lnd(i,j)
            !IPD_Data(nb)%Sfcprop%sfalb_lnd_bck(ix) = mn_phys%sfalb_lnd_bck(i,j)
            !IPD_Data(nb)%Radtend%semis(ix) = mn_phys%semis(i,j)
            !IPD_Data(nb)%Sfcprop%semisbase(ix) = mn_phys%semisbase(i,j)
            !IPD_Data(nb)%Radtend%sfalb(ix) = mn_phys%sfalb(i,j)

            IPD_Data(nb)%IntDiag%u10m(ix) = mn_phys%u10m(i,j)
            IPD_Data(nb)%IntDiag%v10m(ix) = mn_phys%v10m(i,j)
            IPD_Data(nb)%Sfcprop%tprcp(ix) = mn_phys%tprcp(i,j)

            do k = 1, IPD_Control%nmtvr
              IPD_Data(nb)%Sfcprop%hprime(ix,k) = mn_phys%hprime(i,j,k)
            enddo

            IPD_Data(nb)%Sfcprop%lakefrac(ix) = mn_phys%lakefrac(i,j)
            IPD_Data(nb)%Sfcprop%lakedepth(ix) = mn_phys%lakedepth(i,j)

            IPD_Data(nb)%Sfcprop%canopy(ix) = mn_phys%canopy(i,j)
            IPD_Data(nb)%Sfcprop%vfrac(ix)  = mn_phys%vegfrac(i,j)
            IPD_Data(nb)%Sfcprop%uustar(ix) = mn_phys%uustar(i,j)
            IPD_Data(nb)%Sfcprop%shdmin(ix) = mn_phys%shdmin(i,j)
            IPD_Data(nb)%Sfcprop%shdmax(ix) = mn_phys%shdmax(i,j)

            ! Set roughness lengths to physically reasonable values if they have fill value (possible at coastline)
            ! sea/land mask array (sea:0,land:1,sea-ice:2)
            if (nint(IPD_data(nb)%Sfcprop%slmsk(ix)) .eq. 1 .and. mn_phys%zorll(i,j) .gt. 1e6) then
              IPD_Data(nb)%Sfcprop%zorll(ix)  = 82.0   !
            else
              IPD_Data(nb)%Sfcprop%zorll(ix)  = mn_phys%zorll(i,j)
            endif

            if (nint(IPD_data(nb)%Sfcprop%slmsk(ix)) .eq. 0 .and. mn_phys%zorlw(i,j) .gt. 1e6) then
              IPD_Data(nb)%Sfcprop%zorlw(ix)  = 83.0   !
            else
              IPD_Data(nb)%Sfcprop%zorlw(ix)  = mn_phys%zorlw(i,j)
            endif

            if (nint(IPD_data(nb)%Sfcprop%slmsk(ix)) .eq. 0 .and. mn_phys%zorlwav(i,j) .gt. 1e6) then
              IPD_Data(nb)%Sfcprop%zorlwav(ix)  = 84.0   !
            else
              IPD_Data(nb)%Sfcprop%zorlwav(ix)  = mn_phys%zorlwav(i,j)
            endif

            if (mn_phys%zorl(i,j) .gt. 1e6) then
              IPD_Data(nb)%Sfcprop%zorl(ix)   = 85.0
            else
              IPD_Data(nb)%Sfcprop%zorl(ix)   = mn_phys%zorl(i,j)
            endif

            IPD_Data(nb)%Sfcprop%tsfco(ix)  = mn_phys%tsfco(i,j)
            IPD_Data(nb)%Sfcprop%tsfcl(ix)  = mn_phys%tsfcl(i,j)
            IPD_Data(nb)%Sfcprop%tsfc(ix)   = mn_phys%tsfc(i,j)

            ! Set albedo values to physically reasonable values if they have negative fill values.
            if (mn_phys%albdirvis_lnd (i,j) .ge. 0.0) then
              IPD_Data(nb)%Sfcprop%albdirvis_lnd (ix)   = mn_phys%albdirvis_lnd (i,j)
            else
              IPD_Data(nb)%Sfcprop%albdirvis_lnd (ix)   = 0.5
            endif

            if (mn_phys%albdirnir_lnd (i,j) .ge. 0.0) then
              IPD_Data(nb)%Sfcprop%albdirnir_lnd (ix)   = mn_phys%albdirnir_lnd (i,j)
            else
              IPD_Data(nb)%Sfcprop%albdirnir_lnd (ix)   = 0.5
            endif

            if (mn_phys%albdifvis_lnd (i,j) .ge. 0.0) then
              IPD_Data(nb)%Sfcprop%albdifvis_lnd (ix)   = mn_phys%albdifvis_lnd (i,j)
            else
              IPD_Data(nb)%Sfcprop%albdifvis_lnd (ix)   = 0.5
            endif

            if (mn_phys%albdifnir_lnd (i,j) .ge. 0.0) then
              IPD_Data(nb)%Sfcprop%albdifnir_lnd (ix)   = mn_phys%albdifnir_lnd (i,j)
            else
              IPD_Data(nb)%Sfcprop%albdifnir_lnd (ix)   = 0.5
            endif

            ! Cloud properties
            IPD_Data(nb)%Cldprop%cv(ix) = mn_phys%cv(i,j)
            IPD_Data(nb)%Cldprop%cvt(ix) = mn_phys%cvt(i,j)
            IPD_Data(nb)%Cldprop%cvb(ix) = mn_phys%cvb(i,j)

            do nv = 1, IPD_Control%ntot2d
              IPD_Data(nb)%Tbd%phy_f2d(ix, nv) = mn_phys%phy_f2d(i,j,nv)
            enddo

            do k = 1, IPD_Control%levs
              do nv = 1, IPD_Control%ntot3d
                IPD_Data(nb)%Tbd%phy_f3d(ix, k, nv) = mn_phys%phy_f3d(i,j,k,nv)
              enddo
            enddo
          endif

          if (move_nsst) then
            IPD_Data(nb)%Sfcprop%tref(ix)    = mn_phys%tref(i,j)
            IPD_Data(nb)%Sfcprop%z_c(ix)     = mn_phys%z_c(i,j)
            IPD_Data(nb)%Sfcprop%c_0(ix)     = mn_phys%c_0(i,j)
            IPD_Data(nb)%Sfcprop%c_d(ix)     = mn_phys%c_d(i,j)
            IPD_Data(nb)%Sfcprop%w_0(ix)     = mn_phys%w_0(i,j)
            IPD_Data(nb)%Sfcprop%w_d(ix)     = mn_phys%w_d(i,j)
            IPD_Data(nb)%Sfcprop%xt(ix)      = mn_phys%xt(i,j)
            IPD_Data(nb)%Sfcprop%xs(ix)      = mn_phys%xs(i,j)
            IPD_Data(nb)%Sfcprop%xu(ix)      = mn_phys%xu(i,j)
            IPD_Data(nb)%Sfcprop%xv(ix)      = mn_phys%xv(i,j)
            IPD_Data(nb)%Sfcprop%xz(ix)      = mn_phys%xz(i,j)
            IPD_Data(nb)%Sfcprop%zm(ix)      = mn_phys%zm(i,j)
            IPD_Data(nb)%Sfcprop%xtts(ix)    = mn_phys%xtts(i,j)
            IPD_Data(nb)%Sfcprop%xzts(ix)    = mn_phys%xzts(i,j)
            IPD_Data(nb)%Sfcprop%d_conv(ix)  = mn_phys%d_conv(i,j)
            IPD_Data(nb)%Sfcprop%dt_cool(ix) = mn_phys%dt_cool(i,j)
            IPD_Data(nb)%Sfcprop%qrain(ix)   = mn_phys%qrain(i,j)
          endif

          ! Check if stype and vtype are properly set for land points.  Set to reasonable values if they have fill values.
          if ( (int(IPD_data(nb)%Sfcprop%slmsk(ix)) .eq. 1) )  then

            if (IPD_data(nb)%Sfcprop%vtype(ix) .lt. 0.5) then
              print '("[INFO] WDR FIXPHYS resetting vtype from 0. npe=",I0," i,j=",I0,",",I0," lat=",F10.3," lon=",F10.3)', this_pe, i,j, IPD_data(nb)%Grid%xlat_d(ix), IPD_data(nb)%Grid%xlon_d(ix)-360.0
              IPD_data(nb)%Sfcprop%vtype(ix) = 7    ! Force to grassland
            endif

            if (IPD_data(nb)%Sfcprop%stype(ix) .lt. 0.5) then
              print '("[INFO] WDR FIXPHYS resetting stype from 0. npe=",I0," i,j=",I0,",",I0," lat=",F10.3," lon=",F10.3)', this_pe, i,j, IPD_data(nb)%Grid%xlat_d(ix), IPD_data(nb)%Grid%xlon_d(ix)-360.0
              IPD_data(nb)%Sfcprop%stype(ix) = 3    ! Force to sandy loam
            endif

            if (IPD_data(nb)%Sfcprop%vtype_save(ix) .lt. 0.5) then
              IPD_data(nb)%Sfcprop%vtype_save(ix) = 7    ! Force to grassland
            endif
            if (IPD_data(nb)%Sfcprop%stype_save(ix) .lt. 0.5) then
              IPD_data(nb)%Sfcprop%stype_save(ix) = 3    ! Force to sandy loam
            endif

          endif
        enddo
      enddo
    endif

    if (debug_log) print '("[INFO] WDR end mn_phys_apply_temp_variables. npe=",I0," n=",I0)', this_pe, n

  end subroutine mn_phys_apply_temp_variables


  !>@brief The subroutine 'mn_physfill_nest_halos_from_parent' transfers data from the coarse grid to the nest edge
  !>@details This subroutine must run on parent and nest PEs to complete the data transfers
  subroutine mn_phys_fill_nest_halos_from_parent(Atm, IPD_Control, IPD_Data, mn_static, n, child_grid_num, is_fine_pe, nest_domain, nz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)            !< Array of atmospheric data
    type(IPD_control_type), intent(in)                       :: IPD_Control       !< Physics metadata
    type(IPD_data_type), intent(inout)                       :: IPD_Data(:)       !< Physics variable data
    type(mn_surface_grids), intent(in)                       :: mn_static         !< Static data
    integer, intent(in)                                      :: n, child_grid_num !< Current grid number, child grid number
    logical, intent(in)                                      :: is_fine_pe        !< Is this a nest PE?
    type(nest_domain_type), intent(inout)                    :: nest_domain       !< Nest domain for FMS
    integer, intent(in)                                      :: nz                !< Number of vertical levels

    integer  :: position, position_u, position_v
    integer  :: interp_type, interp_type_u, interp_type_v, interp_type_lmask
    integer  :: x_refine, y_refine
    type(fv_moving_nest_physics_type), pointer :: mn_phys

    interp_type = 1        ! cell-centered A-grid
    interp_type_u = 4      ! D-grid
    interp_type_v = 4      ! D-grid
    interp_type_lmask = 7  ! land mask, cell-centered A-grid

    position = CENTER
    position_u = NORTH
    position_v = EAST

    x_refine = Atm(child_grid_num)%neststruct%refinement
    y_refine = x_refine

    mn_phys => Moving_nest(n)%mn_phys

    !  Fill centered-grid variables

    call fill_nest_halos_from_parent("ts", mn_phys%ts, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
        Atm(child_grid_num)%neststruct%ind_h, &
        x_refine, y_refine, &
        is_fine_pe, nest_domain, position)

    if (move_physics) then
      call fill_nest_halos_from_parent("smc", mn_phys%smc, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, IPD_Control%lsoil)
      call fill_nest_halos_from_parent("stc", mn_phys%stc, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, IPD_Control%lsoil)
      call fill_nest_halos_from_parent("slc", mn_phys%slc, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, IPD_Control%lsoil)

      call fill_nest_halos_from_parent("phy_f2d", mn_phys%phy_f2d, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, IPD_Control%ntot2d)

      call fill_nest_halos_from_parent("phy_f3d", mn_phys%phy_f3d, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, IPD_Control%levs)

      !!  Surface variables

      !call fill_nest_halos_from_parent("sfalb_lnd", mn_phys%sfalb_lnd, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
      !     Atm(child_grid_num)%neststruct%ind_h, &
      !     x_refine, y_refine, &
      !     is_fine_pe, nest_domain, position)

      ! sea/land mask array (sea:0,land:1,sea-ice:2)

      call fill_nest_halos_from_parent_masked("emis_lnd", mn_phys%emis_lnd, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, 1, 0.5D0)

      call fill_nest_halos_from_parent_masked("emis_ice", mn_phys%emis_ice, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, 2, 0.5D0)

      call fill_nest_halos_from_parent_masked("emis_wat", mn_phys%emis_wat, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, 0, 0.5D0)

      !call fill_nest_halos_from_parent("sfalb_lnd_bck", mn_phys%sfalb_lnd_bck, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
      !     Atm(child_grid_num)%neststruct%ind_h, &
      !     x_refine, y_refine, &
      !     is_fine_pe, nest_domain, position)


      !call fill_nest_halos_from_parent("semis", mn_phys%semis, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
      !     Atm(child_grid_num)%neststruct%ind_h, &
      !     x_refine, y_refine, &
      !     is_fine_pe, nest_domain, position)
      !call fill_nest_halos_from_parent("semisbase", mn_phys%semisbase, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
      !     Atm(child_grid_num)%neststruct%ind_h, &
      !     x_refine, y_refine, &
      !     is_fine_pe, nest_domain, position)
      !call fill_nest_halos_from_parent("sfalb", mn_phys%sfalb, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
      !     Atm(child_grid_num)%neststruct%ind_h, &
      !     x_refine, y_refine, &
      !     is_fine_pe, nest_domain, position)


      call fill_nest_halos_from_parent("u10m", mn_phys%u10m, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("v10m", mn_phys%v10m, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("tprcp", mn_phys%tprcp, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)

      call fill_nest_halos_from_parent("hprime", mn_phys%hprime, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, IPD_Control%nmtvr)

      call fill_nest_halos_from_parent("lakefrac", mn_phys%lakefrac, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("lakedepth", mn_phys%lakedepth, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)

      call fill_nest_halos_from_parent("canopy", mn_phys%canopy, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("vegfrac", mn_phys%vegfrac, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("uustar", mn_phys%uustar, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("shdmin", mn_phys%shdmin, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("shdmax", mn_phys%shdmax, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("zorl", mn_phys%zorl, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)

      call fill_nest_halos_from_parent_masked("zorll", mn_phys%zorll, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, 1, 86.0D0)
      call fill_nest_halos_from_parent_masked("zorlwav", mn_phys%zorlwav, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, 0, 77.0D0)
      call fill_nest_halos_from_parent_masked("zorlw", mn_phys%zorlw, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, 0, 78.0D0)

      call fill_nest_halos_from_parent("tsfco", mn_phys%tsfco, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("tsfcl", mn_phys%tsfcl, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("tsfc", mn_phys%tsfc, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)

      call fill_nest_halos_from_parent_masked("albdirvis_lnd", mn_phys%albdirvis_lnd, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, 1, 0.5D0)
      call fill_nest_halos_from_parent_masked("albdirnir_lnd", mn_phys%albdirnir_lnd, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, 1, 0.5D0)
      call fill_nest_halos_from_parent_masked("albdifvis_lnd", mn_phys%albdifvis_lnd, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, 1, 0.5D0)
      call fill_nest_halos_from_parent_masked("albdifnir_lnd", mn_phys%albdifnir_lnd, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, 1, 0.5D0)



      call fill_nest_halos_from_parent("cv", mn_phys%cv, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("cvt", mn_phys%cvt, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("cvb", mn_phys%cvb, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
    endif

    if (move_nsst) then

      call fill_nest_halos_from_parent("tref", mn_phys%tref, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("z_c", mn_phys%z_c, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("c_0", mn_phys%c_0, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("c_d", mn_phys%c_d, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("w_0", mn_phys%w_0, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("w_d", mn_phys%w_d, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("xt", mn_phys%xt, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("xs", mn_phys%xs, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("xu", mn_phys%xu, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("xv", mn_phys%xv, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("xz", mn_phys%xz, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("zm", mn_phys%zm, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("xtts", mn_phys%xtts, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("xzts", mn_phys%xzts, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("d_conv", mn_phys%d_conv, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("dt_cool", mn_phys%dt_cool, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("qrain", mn_phys%qrain, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)

    endif

  end subroutine mn_phys_fill_nest_halos_from_parent

  !>@brief The subroutine 'mn_phys_fill_intern_nest_halos' fills the intenal nest halos for the physics variables
  !>@details This subroutine is only called for the nest PEs.
  subroutine mn_phys_fill_intern_nest_halos(moving_nest, IPD_Control, IPD_Data, domain_fine, is_fine_pe)
    type(fv_moving_nest_type), target, intent(inout) :: moving_nest         !< Single instance of moving nest data
    type(IPD_control_type), intent(in)               :: IPD_Control         !< Physics metadata
    type(IPD_data_type), intent(inout)               :: IPD_Data(:)         !< Physics variable data
    type(domain2d), intent(inout)                    :: domain_fine         !< Domain structure for this nest
    logical, intent(in)                              :: is_fine_pe          !< Is nest PE - should be True.  Argument is redundant.

    type(fv_moving_nest_physics_type), pointer :: mn_phys

    mn_phys => moving_nest%mn_phys

    call mn_var_fill_intern_nest_halos(mn_phys%ts, domain_fine, is_fine_pe)   !! Skin Temp/SST
    if (move_physics) then
      call mn_var_fill_intern_nest_halos(mn_phys%smc, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%stc, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%slc, domain_fine, is_fine_pe)

      call mn_var_fill_intern_nest_halos(mn_phys%phy_f2d, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%phy_f3d, domain_fine, is_fine_pe)

      call mn_var_fill_intern_nest_halos(mn_phys%emis_lnd, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%emis_ice, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%emis_wat, domain_fine, is_fine_pe)

      !call mn_var_fill_intern_nest_halos(mn_phys%sfalb_lnd, domain_fine, is_fine_pe)
      !call mn_var_fill_intern_nest_halos(mn_phys%sfalb_lnd_bck, domain_fine, is_fine_pe)
      !call mn_var_fill_intern_nest_halos(mn_phys%semis, domain_fine, is_fine_pe)
      !call mn_var_fill_intern_nest_halos(mn_phys%semisbase, domain_fine, is_fine_pe)
      !call mn_var_fill_intern_nest_halos(mn_phys%sfalb, domain_fine, is_fine_pe)

      call mn_var_fill_intern_nest_halos(mn_phys%u10m, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%v10m, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%tprcp, domain_fine, is_fine_pe)

      call mn_var_fill_intern_nest_halos(mn_phys%hprime, domain_fine, is_fine_pe)

      call mn_var_fill_intern_nest_halos(mn_phys%lakefrac, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%lakedepth, domain_fine, is_fine_pe)

      call mn_var_fill_intern_nest_halos(mn_phys%canopy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%vegfrac, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%uustar, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%shdmin, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%shdmax, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%zorl, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%zorll, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%zorlwav, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%zorlw, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%tsfco, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%tsfcl, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%tsfc, domain_fine, is_fine_pe)

      call mn_var_fill_intern_nest_halos(mn_phys%albdirvis_lnd, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%albdirnir_lnd, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%albdifvis_lnd, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%albdifnir_lnd, domain_fine, is_fine_pe)

      call mn_var_fill_intern_nest_halos(mn_phys%cv, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%cvt, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%cvb, domain_fine, is_fine_pe)
    endif

    if (move_nsst) then
      call mn_var_fill_intern_nest_halos(mn_phys%tref, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%z_c, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%c_0, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%c_d, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%w_0, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%w_d, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%xt, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%xs, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%xu, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%xv, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%xz, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%zm, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%xtts, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%xzts, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%d_conv, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%dt_cool, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%qrain, domain_fine, is_fine_pe)
    endif

  end subroutine mn_phys_fill_intern_nest_halos

  !>@brief The subroutine 'mn_phys_shift_data' shifts the variable in the nest, including interpolating at the leading edge
  !>@details This subroutine is called for the nest and parent PEs.
  subroutine mn_phys_shift_data(Atm, IPD_Control, IPD_Data, n, child_grid_num, wt_h, wt_u, wt_v, &
      delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, nz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)                  !< Array of atmospheric data
    type(IPD_control_type), intent(in)                       :: IPD_Control             !< Physics metadata
    type(IPD_data_type), intent(inout)                       :: IPD_Data(:)             !< Physics variable data
    integer, intent(in)                                      :: n, child_grid_num       !< Current grid number, child grid number
    real, allocatable, intent(in)                            :: wt_h(:,:,:), wt_u(:,:,:), wt_v(:,:,:) !< Interpolation weights
    integer, intent(in)                                      :: delta_i_c, delta_j_c    !< Nest motion in i,j direction
    integer, intent(in)                                      :: x_refine, y_refine      !< Nest refinement
    logical, intent(in)                                      :: is_fine_pe              !< Is this the nest PE?
    type(nest_domain_type), intent(inout)                    :: nest_domain             !< Nest domain structure
    integer, intent(in)                                      :: nz                      !< Number of vertical levels

    ! Constants for mpp calls
    integer  :: interp_type   = 1    ! cell-centered A-grid
    integer  :: interp_type_u = 4    ! D-grid
    integer  :: interp_type_v = 4    ! D-grid
    integer  :: position      = CENTER
    integer  :: position_u    = NORTH
    integer  :: position_v    = EAST
    type(fv_moving_nest_physics_type), pointer :: mn_phys

    mn_phys => Moving_nest(n)%mn_phys

    !! Skin temp/SST
    call mn_var_shift_data(mn_phys%ts, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
        delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)

    if (move_physics) then
      !! Soil variables
      call mn_var_shift_data(mn_phys%smc, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, IPD_Control%lsoil)
      call mn_var_shift_data(mn_phys%stc, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, IPD_Control%lsoil)
      call mn_var_shift_data(mn_phys%slc, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, IPD_Control%lsoil)

      !! Physics arrays
      call mn_var_shift_data(mn_phys%phy_f2d, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, IPD_control%ntot2d)

      call mn_var_shift_data(mn_phys%phy_f3d, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, IPD_Control%levs)

      ! Surface variables

      call mn_var_shift_data(mn_phys%emis_lnd, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%emis_ice, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%emis_wat, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)


      !call mn_var_shift_data(mn_phys%sfalb_lnd, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
      !  delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      !call mn_var_shift_data(mn_phys%sfalb_lnd_bck, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
      !  delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      !call mn_var_shift_data(mn_phys%semis, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
      !  delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      !call mn_var_shift_data(mn_phys%semisbase, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
      !  delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      !call mn_var_shift_data(mn_phys%sfalb, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
      !  delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)

      call mn_var_shift_data(mn_phys%u10m, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%v10m, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%tprcp, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%hprime, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, IPD_Control%nmtvr)
      call mn_var_shift_data(mn_phys%lakefrac, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%lakedepth, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%canopy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%vegfrac, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%uustar, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%shdmin, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%shdmax, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%zorl, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%zorll, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%zorlwav, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%zorlw, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%tsfco, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%tsfcl, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%tsfc, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%albdirvis_lnd, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%albdirnir_lnd, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%albdifvis_lnd, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%albdifnir_lnd, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%cv, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%cvt, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%cvb, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
    endif

    if (move_nsst) then
      call mn_var_shift_data(mn_phys%tref, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%z_c, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%c_0, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%c_d, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%w_0, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%w_d, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%xt, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%xs, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%xu, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%xv, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%xz, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%zm, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%xtts, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%xzts, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%d_conv, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%dt_cool, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%qrain, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
    endif

  end subroutine mn_phys_shift_data

  !>@brief The subroutine 'mn_phys_dump_to_netcdf' dumps physics variables to debugging netCDF files
  !>@details This subroutine is called for the nest and parent PEs.
  subroutine mn_phys_dump_to_netcdf(Atm, Atm_block, IPD_Control, IPD_Data, time_val, file_prefix, is_fine_pe, domain_coarse, domain_fine, nz)
    type(fv_atmos_type), intent(in)            :: Atm                           !< Single instance of atmospheric data
    type (block_control_type), intent(in)      :: Atm_block                     !< Physics block layout
    type(IPD_control_type), intent(in)         :: IPD_Control                   !< Physics metadata
    type(IPD_data_type), intent(in)            :: IPD_Data(:)                   !< Physics variable data
    integer, intent(in)                        :: time_val                      !< Timestep number for filename
    character(len=*), intent(in)               :: file_prefix                   !< Prefix for output netCDF filenames
    logical, intent(in)                        :: is_fine_pe                    !< Is this the nest PE?
    type(domain2d), intent(in)                 :: domain_coarse, domain_fine    !< Domain structures for parent and nest
    integer, intent(in)                        :: nz                            !< Number of vertical levels

    integer :: is, ie, js, je
    integer :: nb, blen, i, j, k, ix, nv
    integer :: this_pe

    integer            :: n_moist
    character(len=16)  :: out_var_name, phys_var_name
    integer            :: position = CENTER

    ! Coerce the double precision variables from physics into single precision for debugging netCDF output
    ! Does not affect values used in calculations.
    ! TODO do we want to dump these as double precision??
    real, allocatable :: smc_pr_local (:,:,:)  !< soil moisture content
    real, allocatable :: stc_pr_local (:,:,:)  !< soil temperature
    real, allocatable :: slc_pr_local (:,:,:)  !< soil liquid water content
    real, allocatable, dimension(:,:) :: sealand_pr_local, deep_soil_t_pr_local, soil_type_pr_local, veg_type_pr_local, slope_type_pr_local, max_snow_alb_pr_local
    real, allocatable, dimension(:,:) :: tsfco_pr_local, tsfcl_pr_local, tsfc_pr_local, vegfrac_pr_local
    real, allocatable, dimension(:,:) :: tref_pr_local, c_0_pr_local, xt_pr_local,  xu_pr_local,  xv_pr_local, ifd_pr_local
    real, allocatable, dimension(:,:) :: facsf_pr_local, facwf_pr_local
    real, allocatable, dimension(:,:) :: alvsf_pr_local, alvwf_pr_local, alnsf_pr_local, alnwf_pr_local
    real, allocatable, dimension(:,:) :: zorl_pr_local, zorll_pr_local, zorlw_pr_local, zorli_pr_local
    real, allocatable :: phy_f2d_pr_local (:,:,:)
    real, allocatable :: phy_f3d_pr_local (:,:,:,:)
    real, allocatable :: lakefrac_pr_local (:,:)  !< lake fraction
    real, allocatable :: landfrac_pr_local (:,:)  !< land fraction
    real, allocatable :: emis_lnd_pr_local (:,:)  !< emissivity land

    this_pe = mpp_pe()

    !  Skin temp/SST
    call mn_var_dump_to_netcdf(Atm%ts, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "SSTK")
    !  Terrain height == phis / grav
    call mn_var_dump_to_netcdf(Atm%phis / grav, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "orog")

    ! sgh and oro were only fully allocated if fv_land is True
    !      if false, oro is (1,1), and sgh is not allocated
    if ( Atm%flagstruct%fv_land ) then
      ! land frac --  called oro in fv_array.F90
      call mn_var_dump_to_netcdf(Atm%oro, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "LFRAC")
      ! terrain standard deviation --  called sgh in fv_array.F90
      call mn_var_dump_to_netcdf(Atm%sgh, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "STDDEV")
    endif

    is = Atm%bd%is
    ie = Atm%bd%ie
    js = Atm%bd%js
    je = Atm%bd%je

    if (debug_log) print '("[INFO] WDR mn_phys_dump_to_netcdf. npe=",I0," is=",I0," ie=",I0," js=",I0," je=",I0)', this_pe, is, ie, js, je

    ! Just allocate compute domain size here for outputs;  the nest moving code also has halos added, but we don't need them here.
    if (move_physics) then
      allocate ( smc_pr_local(is:ie, js:je, IPD_Control%lsoil) )
      allocate ( stc_pr_local(is:ie, js:je, IPD_Control%lsoil) )
      allocate ( slc_pr_local(is:ie, js:je, IPD_Control%lsoil) )
      allocate ( sealand_pr_local(is:ie, js:je) )
      allocate ( lakefrac_pr_local(is:ie, js:je) )
      allocate ( landfrac_pr_local(is:ie, js:je) )
      allocate ( emis_lnd_pr_local(is:ie, js:je) )
      allocate ( phy_f2d_pr_local(is:ie, js:je, IPD_Control%ntot2d) )
      allocate ( phy_f3d_pr_local(is:ie, js:je, IPD_Control%levs, IPD_Control%ntot3d) )
      allocate ( tsfco_pr_local(is:ie, js:je) )
      allocate ( tsfcl_pr_local(is:ie, js:je) )
      allocate ( tsfc_pr_local(is:ie, js:je) )
      allocate ( vegfrac_pr_local(is:ie, js:je) )
      allocate ( alvsf_pr_local(is:ie, js:je) )
      allocate ( alvwf_pr_local(is:ie, js:je) )
      allocate ( alnsf_pr_local(is:ie, js:je) )
      allocate ( alnwf_pr_local(is:ie, js:je) )
      allocate ( deep_soil_t_pr_local(is:ie, js:je) )
      allocate ( soil_type_pr_local(is:ie, js:je) )
      !allocate ( veg_frac_pr_local(is:ie, js:je) )
      allocate ( veg_type_pr_local(is:ie, js:je) )
      allocate ( slope_type_pr_local(is:ie, js:je) )
      allocate ( max_snow_alb_pr_local(is:ie, js:je) )
      allocate ( facsf_pr_local(is:ie, js:je) )
      allocate ( facwf_pr_local(is:ie, js:je) )
      allocate ( zorl_pr_local(is:ie, js:je) )
      allocate ( zorll_pr_local(is:ie, js:je) )
      allocate ( zorlw_pr_local(is:ie, js:je) )
      allocate ( zorli_pr_local(is:ie, js:je) )
    endif

    if (move_nsst) then
      allocate ( tref_pr_local(is:ie, js:je) )
      allocate ( c_0_pr_local(is:ie, js:je) )
      allocate ( xt_pr_local(is:ie, js:je) )
      allocate ( xu_pr_local(is:ie, js:je) )
      allocate ( xv_pr_local(is:ie, js:je) )
      allocate ( ifd_pr_local(is:ie, js:je) )
    endif

    if (move_physics) then
      smc_pr_local = +99999.9
      stc_pr_local = +99999.9
      slc_pr_local = +99999.9
      sealand_pr_local = +99999.9
      lakefrac_pr_local = +99999.9
      landfrac_pr_local = +99999.9
      emis_lnd_pr_local = +99999.9
      phy_f2d_pr_local = +99999.9
      phy_f3d_pr_local = +99999.9
      tsfco_pr_local = +99999.9
      tsfcl_pr_local = +99999.9
      tsfc_pr_local = +99999.9
      vegfrac_pr_local = +99999.9
      alvsf_pr_local = +99999.9
      alvwf_pr_local = +99999.9
      alnsf_pr_local = +99999.9
      alnwf_pr_local = +99999.9
    endif
    if (move_nsst) then
      tref_pr_local = +99999.9
      c_0_pr_local = +99999.9
      xt_pr_local = +99999.9
      xu_pr_local = +99999.9
      xv_pr_local = +99999.9
      ifd_pr_local = +99999.9
    endif

    do nb = 1,Atm_block%nblks
      blen = Atm_block%blksz(nb)
      do ix = 1, blen
        i = Atm_block%index(nb)%ii(ix)
        j = Atm_block%index(nb)%jj(ix)

        if (move_physics) then
          do k = 1, IPD_Control%lsoil
            ! Use real() to lower the precision
            smc_pr_local(i,j,k) = real(IPD_Data(nb)%Sfcprop%smc(ix,k))
            stc_pr_local(i,j,k) = real(IPD_Data(nb)%Sfcprop%stc(ix,k))
            slc_pr_local(i,j,k) = real(IPD_Data(nb)%Sfcprop%slc(ix,k))
          enddo

          sealand_pr_local(i,j) = real(IPD_Data(nb)%Sfcprop%slmsk(ix))
          lakefrac_pr_local(i,j) = real(IPD_Data(nb)%Sfcprop%lakefrac(ix))
          landfrac_pr_local(i,j) = real(IPD_Data(nb)%Sfcprop%landfrac(ix))
          emis_lnd_pr_local(i,j) = real(IPD_Data(nb)%Sfcprop%emis_lnd(ix))
          deep_soil_t_pr_local(i, j) = IPD_data(nb)%Sfcprop%tg3(ix)
          soil_type_pr_local(i, j) = IPD_data(nb)%Sfcprop%stype(ix)
          !veg_frac_pr_local(i, j) = IPD_data(nb)%Sfcprop%vfrac(ix)
          veg_type_pr_local(i, j) = IPD_data(nb)%Sfcprop%vtype(ix)
          slope_type_pr_local(i, j) = IPD_data(nb)%Sfcprop%slope(ix)
          facsf_pr_local(i, j) = IPD_data(nb)%Sfcprop%facsf(ix)
          facwf_pr_local(i, j) = IPD_data(nb)%Sfcprop%facwf(ix)
          zorl_pr_local(i, j) = IPD_data(nb)%Sfcprop%zorl(ix)
          zorlw_pr_local(i, j) = IPD_data(nb)%Sfcprop%zorlw(ix)
          zorll_pr_local(i, j) = IPD_data(nb)%Sfcprop%zorll(ix)
          zorli_pr_local(i, j) = IPD_data(nb)%Sfcprop%zorli(ix)
          max_snow_alb_pr_local(i, j) = IPD_data(nb)%Sfcprop%snoalb(ix)
          tsfco_pr_local(i, j) = IPD_data(nb)%Sfcprop%tsfco(ix)
          tsfcl_pr_local(i, j) = IPD_data(nb)%Sfcprop%tsfcl(ix)
          tsfc_pr_local(i, j)  = IPD_data(nb)%Sfcprop%tsfc(ix)
          vegfrac_pr_local(i, j) = IPD_data(nb)%Sfcprop%vfrac(ix)
          alvsf_pr_local(i, j) = IPD_data(nb)%Sfcprop%alvsf(ix)
          alvwf_pr_local(i, j) = IPD_data(nb)%Sfcprop%alvwf(ix)
          alnsf_pr_local(i, j) = IPD_data(nb)%Sfcprop%alnsf(ix)
          alnwf_pr_local(i, j) = IPD_data(nb)%Sfcprop%alnwf(ix)

          do nv = 1, IPD_Control%ntot2d
            ! Use real() to lower the precision
            phy_f2d_pr_local(i,j,nv) = real(IPD_Data(nb)%Tbd%phy_f2d(ix, nv))
          enddo

          do k = 1, IPD_Control%levs
            do nv = 1, IPD_Control%ntot3d
              ! Use real() to lower the precision
              phy_f3d_pr_local(i,j,k,nv) = real(IPD_Data(nb)%Tbd%phy_f3d(ix, k, nv))
            enddo
          enddo
        endif

        if (move_nsst) then
          tref_pr_local(i,j) = IPD_data(nb)%Sfcprop%tref(ix)
          c_0_pr_local(i,j) = IPD_data(nb)%Sfcprop%c_0(ix)
          xt_pr_local(i,j) = IPD_data(nb)%Sfcprop%xt(ix)
          xu_pr_local(i,j) = IPD_data(nb)%Sfcprop%xu(ix)
          xv_pr_local(i,j) = IPD_data(nb)%Sfcprop%xv(ix)
          ifd_pr_local(i,j) = IPD_data(nb)%Sfcprop%ifd(ix)
        endif
      enddo
    enddo

    if (move_physics) then
      call mn_var_dump_to_netcdf(stc_pr_local, is_fine_pe, domain_coarse, domain_fine, position, IPD_Control%lsoil, time_val, Atm%global_tile, file_prefix, "SOILT")
      call mn_var_dump_to_netcdf(smc_pr_local, is_fine_pe, domain_coarse, domain_fine, position, IPD_Control%lsoil, time_val, Atm%global_tile, file_prefix, "SOILM")
      call mn_var_dump_to_netcdf(slc_pr_local, is_fine_pe, domain_coarse, domain_fine, position, IPD_Control%lsoil, time_val, Atm%global_tile, file_prefix, "SOILL")
      call mn_var_dump_to_netcdf(sealand_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "LMASK")
      call mn_var_dump_to_netcdf(lakefrac_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "LAKEFRAC")
      call mn_var_dump_to_netcdf(landfrac_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "LANDFRAC")
      call mn_var_dump_to_netcdf(emis_lnd_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "EMISLAND")
      call mn_var_dump_to_netcdf(deep_soil_t_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "DEEPSOIL")
      call mn_var_dump_to_netcdf(soil_type_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "SOILTP")
      !call mn_var_dump_to_netcdf(veg_frac_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "VEGFRAC")
      call mn_var_dump_to_netcdf(veg_type_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "VEGTYPE")
      call mn_var_dump_to_netcdf(slope_type_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "SLOPE")
      call mn_var_dump_to_netcdf(max_snow_alb_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "SNOWALB")
      call mn_var_dump_to_netcdf(tsfco_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "TSFCO")
      call mn_var_dump_to_netcdf(tsfcl_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "TSFCL")
      call mn_var_dump_to_netcdf(tsfc_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "TSFC")
      call mn_var_dump_to_netcdf(vegfrac_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "VEGFRAC")
      call mn_var_dump_to_netcdf(alvsf_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "ALVSF")
      call mn_var_dump_to_netcdf(alvwf_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "ALVWF")
      call mn_var_dump_to_netcdf(alnsf_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "ALNSF")
      call mn_var_dump_to_netcdf(alnwf_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "ALNWF")
      call mn_var_dump_to_netcdf(facsf_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "FACSF")
      call mn_var_dump_to_netcdf(facwf_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "FACWF")
      call mn_var_dump_to_netcdf(zorl_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "ZORL")
      call mn_var_dump_to_netcdf(zorlw_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "ZORLW")
      call mn_var_dump_to_netcdf(zorll_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "ZORLL")
      call mn_var_dump_to_netcdf(zorli_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "ZORLI")

      do nv = 1, IPD_Control%ntot2d
        write (phys_var_name, "(A4,I0.3)")  'PH2D', nv
        call mn_var_dump_to_netcdf(phy_f2d_pr_local(:,:,nv), is_fine_pe, domain_coarse, domain_fine, position, 1, &
            time_val, Atm%global_tile, file_prefix, phys_var_name)
      enddo

      do nv = 1, IPD_Control%ntot3d
        write (phys_var_name, "(A4,I0.3)")  'PH3D', nv
        call mn_var_dump_to_netcdf(phy_f3d_pr_local(:,:,:,nv), is_fine_pe, domain_coarse, domain_fine, position, IPD_Control%levs, &
            time_val, Atm%global_tile, file_prefix, phys_var_name)
      enddo
    endif

    if (move_nsst) then
      call mn_var_dump_to_netcdf(tref_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "TREF")
      call mn_var_dump_to_netcdf(c_0_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "C_0")
      call mn_var_dump_to_netcdf(xt_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "XT")
      call mn_var_dump_to_netcdf(xu_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "XU")
      call mn_var_dump_to_netcdf(xv_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "XV")
      call mn_var_dump_to_netcdf(ifd_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "IFD")
    endif

    if (move_physics) then
      deallocate(smc_pr_local)
      deallocate(stc_pr_local)
      deallocate(slc_pr_local)
      deallocate(lakefrac_pr_local)
      deallocate(landfrac_pr_local)
      deallocate(emis_lnd_pr_local)
      deallocate(sealand_pr_local, deep_soil_t_pr_local, soil_type_pr_local, veg_type_pr_local, max_snow_alb_pr_local)
      deallocate(tsfco_pr_local, tsfcl_pr_local, tsfc_pr_local, vegfrac_pr_local)
      deallocate(alvsf_pr_local, alvwf_pr_local, alnsf_pr_local, alnwf_pr_local)
      deallocate(facsf_pr_local, facwf_pr_local)
      deallocate(zorl_pr_local, zorlw_pr_local, zorll_pr_local, zorli_pr_local)
      deallocate(phy_f2d_pr_local)
      deallocate(phy_f3d_pr_local)
    endif

    if (move_nsst) deallocate(tref_pr_local, c_0_pr_local, xt_pr_local,  xu_pr_local,  xv_pr_local, ifd_pr_local)

    if (debug_log) print '("[INFO] WDR end mn_phys_dump_tp_netcdf npe=",I0)', this_pe

  end subroutine mn_phys_dump_to_netcdf

#endif MOVING_NEST

end module fv_moving_nest_physics_mod
