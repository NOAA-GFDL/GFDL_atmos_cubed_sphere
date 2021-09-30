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

!----------------------------------------------------------
! Moving Nest Physics Variables    W. Ramstrom - 09/29/2021
!----------------------------------------------------------


!*************************************************************************
!>@brief!   Provides Moving Nest functionality for physics variables
!!>@author W. Ramstrom, AOML/HRD  09/29/2021
!
! 
! =======================================================================!
!

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

#define REMAP 1 

module fv_moving_nest_physics_mod
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
  use constants_mod,          only: cp_air, rdgas, grav, rvgas, kappa, pstd_mks, hlv
  use field_manager_mod,      only: MODEL_ATMOS
  use fms_io_mod,             only: read_data, write_data, get_global_att_value, fms_io_init, fms_io_exit
  use fv_arrays_mod,          only: fv_atmos_type, fv_nest_type, fv_grid_type, R_GRID
  use fv_arrays_mod,          only: fv_moving_nest_prog_type, fv_moving_nest_physics_type
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
  use fv_moving_nest_mod,     only: mn_var_fill_intern_nest_halos, mn_var_dump_to_netcdf, mn_var_shift_data

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
  logical :: move_physics = .true.
  logical :: move_nsst = .true.

#include <fms_platform.h>

  type mn_surface_grids
     real, allocatable  :: orog_grid(:,:)                _NULL  ! orography -- raw or filtered depending on namelist option, in meters
     real, allocatable  :: orog_std_grid(:,:)            _NULL  ! terrain standard deviation for gravity wave drag, in meters (?)
     real, allocatable  :: ls_mask_grid(:,:)             _NULL  ! land sea mask -- 0 for ocean/lakes, 1, for land.  Perhaps 2 for sea ice. 
     real, allocatable  :: land_frac_grid(:,:)           _NULL  ! Continuous land fraction - 0.0 ocean, 0.5 half of each, 1.0 all land
     
     real, allocatable  :: parent_orog_grid(:,:)         _NULL  ! parent orography -- only used for terrain_smoother=1.  
     !     raw or filtered depending on namelist option, in meters
     
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
     real, allocatable  :: max_snow_alb_grid(:,:)       _NULL  ! max snow albedo
     ! Snow free albedo
     real, allocatable  :: vis_black_alb_grid(:,:)      _NULL  ! Visible black sky albeo; netCDF file has monthly values
     real, allocatable  :: vis_white_alb_grid(:,:)      _NULL  ! Visible white sky albeo; netCDF file has monthly values
     real, allocatable  :: ir_black_alb_grid(:,:)       _NULL  ! Near IR black sky albeo; netCDF file has monthly values
     real, allocatable  :: ir_white_alb_grid(:,:)       _NULL  ! Near IR white sky albeo; netCDF file has monthly values
     
  end type mn_surface_grids


contains

  subroutine mn_phys_reset_sfc_props(Atm, n, mn_static, Atm_block, IPD_data, ioffset, joffset, refine)
    type(fv_atmos_type), intent(in),allocatable      :: Atm(:)
    integer, intent(in)                              :: n
    type(mn_surface_grids), intent(in)               :: mn_static
    type(block_control_type), intent(in)             :: Atm_block
    type(IPD_data_type), intent(inout)               :: IPD_data(:)
    integer, intent(in)                              :: ioffset, joffset, refine

    ! For iterating through physics/surface vector data
    integer                 :: nb, blen, ix, i_pe, j_pe, i_idx, j_idx
    real(kind=kind_phys)    :: phys_oro
    
    ! Reset the land sea mask from the hires parent data
    !  Reset the variables from the fix_sfc files
    !  Reset the oro and oro_uf from the smoothed values calculated above
    do nb = 1,Atm_block%nblks
       blen = Atm_block%blksz(nb)
       do ix = 1, blen
          i_pe = Atm_block%index(nb)%ii(ix)  
          j_pe = Atm_block%index(nb)%jj(ix)
          
          i_idx = (ioffset-1)*refine + i_pe
          j_idx = (joffset-1)*refine + j_pe
          
          IPD_data(nb)%Sfcprop%slmsk(ix) = mn_static%ls_mask_grid(i_idx, j_idx)
          
          !  IFD values are 0 for land, and 1 for oceans/lakes -- reverse of the land sea mask
          !  Land Sea Mask has values of 0 for oceans/lakes, 1 for land, 2 for sea ice
          !  TODO figure out what ifd should be for sea ice
          if (mn_static%ls_mask_grid(i_idx, j_idx) .eq. 1 ) then
             IPD_data(nb)%Sfcprop%ifd(ix) = 0         ! Land
             IPD_data(nb)%Sfcprop%oceanfrac(ix) = 0   ! Land -- TODO permit fractions
             IPD_data(nb)%Sfcprop%landfrac(ix) = 1    ! Land -- TODO permit fractions
          else
             IPD_data(nb)%Sfcprop%ifd(ix) = 1         ! Ocean
             IPD_data(nb)%Sfcprop%oceanfrac(ix) = 1   ! Ocean -- TODO permit fractions
             IPD_data(nb)%Sfcprop%landfrac(ix) = 0    ! Ocean -- TODO permit fractions
          end if
          
          IPD_data(nb)%Sfcprop%tg3(ix) = mn_static%deep_soil_temp_grid(i_idx, j_idx)
          IPD_data(nb)%Sfcprop%stype(ix) = mn_static%soil_type_grid(i_idx, j_idx)
          
          !IPD_data(nb)%Sfcprop%vfrac(ix) = mn_static%veg_frac_grid(i_idx, j_idx) 
          IPD_data(nb)%Sfcprop%vtype(ix) = mn_static%veg_type_grid(i_idx, j_idx)
          ! Add veg_greenness_grid here, monthly
          
          IPD_data(nb)%Sfcprop%slope(ix) = mn_static%slope_type_grid(i_idx, j_idx)
          
          IPD_data(nb)%Sfcprop%snoalb(ix) = mn_static%max_snow_alb_grid(i_idx, j_idx)
          ! Add Vis/Near IR black/white sky albedo, monthly
          
          ! Reset the orography in the physics arrays, using the smoothed values from above
          phys_oro =  Atm(n)%phis(i_pe, j_pe) / grav
          !if (phys_oro .gt. 15000.0) then
          !   phys_oro = 0.0
          !end if
          IPD_data(nb)%Sfcprop%oro(ix) = phys_oro
          IPD_data(nb)%Sfcprop%oro_uf(ix) = phys_oro
          
       end do
    end do
    
  end subroutine mn_phys_reset_sfc_props

  subroutine mn_reset_phys_latlon(Atm, n, tile_geo, fp_super_tile_geo, Atm_block, IPD_control, IPD_data)
    type(fv_atmos_type), intent(in)      :: Atm(:)
    integer, intent(in)                  :: n
    type(grid_geometry), intent(in)      :: tile_geo
    type(grid_geometry), intent(in)      :: fp_super_tile_geo
    type(block_control_type), intent(in) :: Atm_block
    type(IPD_control_type), intent(in) :: IPD_control
    type(IPD_data_type), intent(inout) :: IPD_data(:)
    !type(time_type), intent(in)     :: time_step
    
    integer :: isc, jsc, iec, jec
    integer :: x, y, fp_i, fp_j !, fpa_i, fpa_j
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

    !print '("WDR mn_reset_phys_latlon AA npe=",I0)', this_pe

    ! This is going to be slow -- replace with better way to calculate index offsets, or pass them from earlier calculations
    call find_nest_alignment(tile_geo, fp_super_tile_geo, nest_x, nest_y, parent_x, parent_y)
    !print '("WDR mn_reset_phys_latlon AB npe=",I0)', this_pe

    do x = isc, iec
       do y = jsc, jec
          fp_i = (x - nest_x) * 2 + parent_x
          fp_j = (y - nest_y) * 2 + parent_y

          !fpa_i = (x - nest_x) * 2 + parent_x - 1
          !fpa_j = (y - nest_y) * 2 + parent_y - 1

          !print '("WDR mn_reset_phys_latlon BB npe=",I0," ix=",I0," x,y=",I0,",",I0)', this_pe, x,y
          lons(x,y) = fp_super_tile_geo%lons(fp_i, fp_j)
          !print '("WDR mn_reset_phys_latlon CC npe=",I0," ix=",I0," x,y=",I0,",",I0)', this_pe, x,y
          lats(x,y) = fp_super_tile_geo%lats(fp_i, fp_j)
          !print '("WDR mn_reset_phys_latlon DD npe=",I0," ix=",I0," x,y=",I0,",",I0)', this_pe, x,y

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
          !print '("WDR mn_reset_phys_latlon EE npe=",I0," ix=",I0," x,y=",I0,",",I0)', this_pe, x,y
       end do
    end do

    !print '("WDR mn_reset_phys_latlon FF npe=",I0," xlon_d=",F15.10," xlat_d=",F15.10," area=",F20.5)', this_pe, IPD_data(1)%Grid%xlon_d(1), IPD_data(1)%Grid%xlat_d(1), IPD_data(1)%Grid%area(1)

    call GFS_grid_populate(IPD_data%Grid, lons, lats, area)
    
    !print '("WDR mn_reset_phys_latlon GG npe=",I0," xlon_d=",F15.10," xlat_d=",F15.10," area=",F20.5)', this_pe, IPD_data(1)%Grid%xlon_d(1), IPD_data(1)%Grid%xlat_d(1), IPD_data(1)%Grid%area(1)

    deallocate(lats)
    deallocate(lons)
    deallocate(area)

  end subroutine mn_reset_phys_latlon




  subroutine mn_phys_fill_temp_variables(Atm, Atm_block, IPD_Control, IPD_Data, n, child_grid_num, is_fine_pe, npz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)
    type (block_control_type), intent(in)            :: Atm_block
    type(IPD_control_type), intent(in)               :: IPD_Control
    type(IPD_data_type), intent(inout)               :: IPD_Data(:)
    integer, intent(in)                              :: n, child_grid_num
    logical, intent(in)                              :: is_fine_pe
    integer, intent(in)                              :: npz

    integer :: isd, ied, jsd, jed
    integer :: is, ie, js, je
    integer :: this_pe

    integer :: nb, blen, i, j, k, ix, nv
    type(fv_moving_nest_physics_type), pointer       :: mn_phys

    this_pe = mpp_pe()


    if (debug_log) print '("[INFO] WDR start mn_phys_fill_temp_variables. npe=",I0," n=",I0)', this_pe, n

    isd = Atm(n)%bd%isd
    ied = Atm(n)%bd%ied
    jsd = Atm(n)%bd%jsd
    jed = Atm(n)%bd%jed
    
    if (debug_log) print '("[INFO] WDR mn_phys_fill_temp_variables. npe=",I0," isd=",I0," ied=",I0," jsd=",I0," jed=",I0)', this_pe, isd, ied, jsd, jed

    is = Atm(n)%bd%is
    ie = Atm(n)%bd%ie
    js = Atm(n)%bd%js
    je = Atm(n)%bd%je
    
    if (debug_log) print '("[INFO] WDR mn_phys_fill_temp_variables. npe=",I0," is=",I0," ie=",I0," js=",I0," je=",I0)', this_pe, is, ie, js, je
    
    mn_phys => Atm(n)%mn_phys

    mn_phys%ts(is:ie, js:je) =  Atm(n)%ts(is:ie, js:je) 

    do nb = 1,Atm_block%nblks
      blen = Atm_block%blksz(nb)
      do ix = 1, blen
         ! Get the indices only once, before iterating through vertical levels or number of variables
         !  Was there a different efficiency from having the k loop outside?
         i = Atm_block%index(nb)%ii(ix)  
         j = Atm_block%index(nb)%jj(ix)

         if (move_physics) then
            do k = 1, IPD_Control%lsoil            
               mn_phys%smc(i,j,k) = IPD_Data(nb)%Sfcprop%smc(ix,k)
               mn_phys%stc(i,j,k) = IPD_Data(nb)%Sfcprop%stc(ix,k)
               mn_phys%slc(i,j,k) = IPD_Data(nb)%Sfcprop%slc(ix,k)
            enddo

            mn_phys%u10m(i,j)  = IPD_Data(nb)%IntDiag%u10m(ix)
            mn_phys%v10m(i,j)  = IPD_Data(nb)%IntDiag%v10m(ix)
            mn_phys%tprcp(i,j)  = IPD_Data(nb)%Sfcprop%tprcp(ix)

            do k = 1, IPD_Control%nmtvr
               mn_phys%hprime(i,j,k)  = IPD_Data(nb)%Sfcprop%hprime(ix,k)
            end do
            
            mn_phys%zorl(i,j)  = IPD_Data(nb)%Sfcprop%zorl(ix)
            mn_phys%alvsf(i,j) = IPD_Data(nb)%Sfcprop%alvsf(ix)
            mn_phys%alvwf(i,j) = IPD_Data(nb)%Sfcprop%alvwf(ix)
            mn_phys%alnsf(i,j) = IPD_Data(nb)%Sfcprop%alnsf(ix)
            mn_phys%alnwf(i,j) = IPD_Data(nb)%Sfcprop%alnwf(ix)
            
            mn_phys%facsf(i,j) = IPD_data(nb)%Sfcprop%facsf(ix)   ! fractional coverage for strong zenith angle albedo
            mn_phys%facwf(i,j) = IPD_data(nb)%Sfcprop%facwf(ix)   ! fractional coverage for weak zenith angle albedo

            mn_phys%canopy(i,j) = IPD_Data(nb)%Sfcprop%canopy(ix)
            mn_phys%vegfrac(i,j)= IPD_Data(nb)%Sfcprop%vfrac(ix)
            mn_phys%uustar(i,j) = IPD_Data(nb)%Sfcprop%uustar(ix)
            mn_phys%shdmin(i,j) = IPD_Data(nb)%Sfcprop%shdmin(ix)
            mn_phys%shdmax(i,j) = IPD_Data(nb)%Sfcprop%shdmax(ix)
            mn_phys%zorll(i,j)  = IPD_Data(nb)%Sfcprop%zorll(ix)
            mn_phys%zorlo(i,j)  = IPD_Data(nb)%Sfcprop%zorlo(ix)
            mn_phys%zorlw(i,j)  = IPD_Data(nb)%Sfcprop%zorlw(ix)
            mn_phys%tsfco(i,j)  = IPD_Data(nb)%Sfcprop%tsfco(ix)
            mn_phys%tsfcl(i,j)  = IPD_Data(nb)%Sfcprop%tsfcl(ix)
            mn_phys%tsfc(i,j)   = IPD_Data(nb)%Sfcprop%tsfc(ix)

            do nv = 1, IPD_Control%ntot2d
               mn_phys%phy_f2d(i,j,nv) = IPD_Data(nb)%Tbd%phy_f2d(ix, nv)
            end do
            
            do k = 1, IPD_Control%levs
               do nv = 1, IPD_Control%ntot3d
                  mn_phys%phy_f3d(i,j,k,nv) = IPD_Data(nb)%Tbd%phy_f3d(ix, k, nv)
               end do
            end do
            
            ! Cloud prop data has x,y dimensions
            mn_phys%cv(i,j)  = IPD_Data(nb)%Cldprop%cv(ix)
            mn_phys%cvt(i,j) = IPD_Data(nb)%Cldprop%cvt(ix)
            mn_phys%cvb(i,j) = IPD_Data(nb)%Cldprop%cvb(ix)
         end if
         
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
            !mn_phys%ifd(i,j)    = IPD_Data(nb)%Sfcprop%ifd(ix)
            mn_phys%dt_cool(i,j)= IPD_Data(nb)%Sfcprop%dt_cool(ix)
            mn_phys%qrain(i,j)  = IPD_Data(nb)%Sfcprop%qrain(ix)
         end if
      enddo
    enddo

    if (debug_log) print '("[INFO] WDR end mn_phys_fill_temp_variables. npe=",I0," n=",I0)', this_pe, n
    
  end subroutine mn_phys_fill_temp_variables

  subroutine mn_phys_apply_temp_variables(Atm, Atm_block, IPD_Control, IPD_Data, n, child_grid_num, is_fine_pe, npz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)
    type (block_control_type), intent(in)            :: Atm_block
    type(IPD_control_type), intent(in)               :: IPD_Control
    type(IPD_data_type), intent(inout)               :: IPD_Data(:)
    integer, intent(in)                              :: n, child_grid_num
    logical, intent(in)                              :: is_fine_pe
    integer, intent(in)                              :: npz

    integer :: is, ie, js, je
    integer :: this_pe
    integer :: nb, blen, i, j ,k, ix, nv
    integer :: bad_values, good_values
    type(fv_moving_nest_physics_type), pointer       :: mn_phys

    this_pe = mpp_pe()
    mn_phys => Atm(n)%mn_phys

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
                end if
             end do
          end do
       else
          do i = Atm(n)%bd%is, Atm(n)%bd%ie
             do j = Atm(n)%bd%js, Atm(n)%bd%je
                if (mn_phys%ts(i,j) .gt. 20000.0) then
                   print '("[WARN] WDR BAD GLOBAL ts value. npe=",I0," ts(",I0,",",I0")=",F12.3)', this_pe, i, j, mn_phys%ts(i,j)
                   bad_values = bad_values + 1
                else
                   good_values = good_values + 1
                end if
             end do
          end do
       end if
       

       i = Atm(n)%bd%is
       j = Atm(n)%bd%js
       
       print '("[WARN] WDR Surface ts value. npe=",I0," ts(",I0,",",I0,")=",F18.3)', this_pe, i, j, mn_phys%ts(i,j)

       print '("INFO] WDR ts values. npe=",I0," good_values=",I0," bad_values=",I0)', this_pe, good_values, bad_values
    end if
       
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

                IPD_Data(nb)%IntDiag%u10m(ix) = mn_phys%u10m(i,j)
                IPD_Data(nb)%IntDiag%v10m(ix) = mn_phys%v10m(i,j)
                IPD_Data(nb)%Sfcprop%tprcp(ix) = mn_phys%tprcp(i,j)

                do k = 1, IPD_Control%nmtvr
                   IPD_Data(nb)%Sfcprop%hprime(ix,k) = mn_phys%hprime(i,j,k)
                end do
                
                IPD_Data(nb)%Sfcprop%zorl(ix) = mn_phys%zorl(i,j)
                IPD_Data(nb)%Sfcprop%alvsf(ix) = mn_phys%alvsf(i,j)
                IPD_Data(nb)%Sfcprop%alvwf(ix) = mn_phys%alvwf(i,j)
                IPD_Data(nb)%Sfcprop%alnsf(ix) = mn_phys%alnsf(i,j)
                IPD_Data(nb)%Sfcprop%alnwf(ix) = mn_phys%alnwf(i,j)

                IPD_Data(nb)%Sfcprop%facsf(ix) = mn_phys%facsf(i,j)
                IPD_Data(nb)%Sfcprop%facwf(ix) = mn_phys%facwf(i,j)
                
                IPD_Data(nb)%Sfcprop%canopy(ix) = mn_phys%canopy(i,j)
                IPD_Data(nb)%Sfcprop%vfrac(ix)  = mn_phys%vegfrac(i,j)
                IPD_Data(nb)%Sfcprop%uustar(ix) = mn_phys%uustar(i,j)
                IPD_Data(nb)%Sfcprop%shdmin(ix) = mn_phys%shdmin(i,j)
                IPD_Data(nb)%Sfcprop%shdmax(ix) = mn_phys%shdmax(i,j)
                IPD_Data(nb)%Sfcprop%zorll(ix)  = mn_phys%zorll(i,j)
                IPD_Data(nb)%Sfcprop%zorlo(ix)  = mn_phys%zorlo(i,j)
                IPD_Data(nb)%Sfcprop%zorlw(ix)  = mn_phys%zorlw(i,j)
                IPD_Data(nb)%Sfcprop%tsfco(ix)  = mn_phys%tsfco(i,j)
                IPD_Data(nb)%Sfcprop%tsfcl(ix)  = mn_phys%tsfcl(i,j)
                IPD_Data(nb)%Sfcprop%tsfc(ix)   = mn_phys%tsfc(i,j)
                
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
             end if

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
                !IPD_Data(nb)%Sfcprop%ifd(ix)    = mn_phys%ifd(i,j)
                IPD_Data(nb)%Sfcprop%dt_cool(ix) = mn_phys%dt_cool(i,j)
                IPD_Data(nb)%Sfcprop%qrain(ix)   = mn_phys%qrain(i,j)
             end if
          enddo
       enddo       
   end if
      
   if (debug_log) print '("[INFO] WDR end mn_phys_apply_temp_variables. npe=",I0," n=",I0)', this_pe, n
   
  end subroutine mn_phys_apply_temp_variables


  subroutine mn_phys_fill_nest_halos_from_parent(Atm, IPD_Control, IPD_Data, n, child_grid_num, is_fine_pe, nest_domain, nz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)
    type(IPD_control_type), intent(in)               :: IPD_Control
    type(IPD_data_type), intent(inout)               :: IPD_Data(:)
    integer, intent(in)                              :: n, child_grid_num
    logical, intent(in)                              :: is_fine_pe
    type(nest_domain_type), intent(inout)            :: nest_domain
    integer, intent(in)                              :: nz


    integer  :: position, position_u, position_v
    integer  :: interp_type, interp_type_u, interp_type_v
    integer  :: x_refine, y_refine
    type(fv_moving_nest_physics_type), pointer :: mn_phys

    !  TODO: examine how the static nesting code handles this
    !  TODO move terrain and surface parameters, including phis

    interp_type = 1    ! cell-centered A-grid
    interp_type_u = 4  ! D-grid
    interp_type_v = 4  ! D-grid

    position = CENTER ! CENTER, NORTH, EAST
    position_u = NORTH
    position_v = EAST

    x_refine = Atm(child_grid_num)%neststruct%refinement
    y_refine = x_refine

    mn_phys => Atm(n)%mn_phys

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
       
       call fill_nest_halos_from_parent("zorl", mn_phys%zorl, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       
       call fill_nest_halos_from_parent("alvsf", mn_phys%alvsf, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       
       call fill_nest_halos_from_parent("alvwf", mn_phys%alvwf, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       
       call fill_nest_halos_from_parent("alnsf", mn_phys%alnsf, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       
       call fill_nest_halos_from_parent("alnwf", mn_phys%alnwf, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)


       call fill_nest_halos_from_parent("facsf", mn_phys%facsf, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)

       call fill_nest_halos_from_parent("facwf", mn_phys%facwf, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       
       !!
       
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
       call fill_nest_halos_from_parent("zorll", mn_phys%zorll, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("zorlo", mn_phys%zorlo, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("zorlw", mn_phys%zorlw, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)


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
    end if

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
       !call fill_nest_halos_from_parent("ifd", mn_phys%ifd, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
       !     Atm(child_grid_num)%neststruct%ind_h, &
       !     x_refine, y_refine, &
       !     is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("dt_cool", mn_phys%dt_cool, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("qrain", mn_phys%qrain, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)

       end if

  end subroutine mn_phys_fill_nest_halos_from_parent


  ! Fill internal nest halos for physics variables
  ! 
  subroutine mn_phys_fill_intern_nest_halos(Atm, IPD_Control, IPD_Data, domain_fine, is_fine_pe)
    type(fv_atmos_type), target, intent(inout)       :: Atm  
    type(IPD_control_type), intent(in)               :: IPD_Control
    type(IPD_data_type), intent(inout)               :: IPD_Data(:)
    type(domain2d), intent(inout)                    :: domain_fine
    logical, intent(in)                              :: is_fine_pe

    type(fv_moving_nest_physics_type), pointer :: mn_phys

    mn_phys => Atm%mn_phys

    call mn_var_fill_intern_nest_halos(mn_phys%ts, domain_fine, is_fine_pe)   !! Skin Temp/SST
    if (move_physics) then
       call mn_var_fill_intern_nest_halos(mn_phys%smc, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%stc, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%slc, domain_fine, is_fine_pe)
       
       call mn_var_fill_intern_nest_halos(mn_phys%phy_f2d, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%phy_f3d, domain_fine, is_fine_pe)

       call mn_var_fill_intern_nest_halos(mn_phys%u10m, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%v10m, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%tprcp, domain_fine, is_fine_pe)

       call mn_var_fill_intern_nest_halos(mn_phys%hprime, domain_fine, is_fine_pe)
       
       call mn_var_fill_intern_nest_halos(mn_phys%zorl, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%alvsf, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%alvwf, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%alnsf, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%alnwf, domain_fine, is_fine_pe)

       call mn_var_fill_intern_nest_halos(mn_phys%facsf, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%facwf, domain_fine, is_fine_pe)
       
       call mn_var_fill_intern_nest_halos(mn_phys%canopy, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%vegfrac, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%uustar, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%shdmin, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%shdmax, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%zorll, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%zorlo, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%zorlw, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%tsfco, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%tsfcl, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%tsfc, domain_fine, is_fine_pe)
       
       call mn_var_fill_intern_nest_halos(mn_phys%cv, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%cvt, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%cvb, domain_fine, is_fine_pe)
    end if

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
       !call mn_var_fill_intern_nest_halos(mn_phys%ifd, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%dt_cool, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%qrain, domain_fine, is_fine_pe)
    end if
    !call check_array(Atm%u, this_pe, "Atm%pt", 100.0, 400.0)

  end subroutine mn_phys_fill_intern_nest_halos


  subroutine mn_phys_shift_data(Atm, IPD_Control, IPD_Data, n, child_grid_num, wt_h, wt_u, wt_v, &
       delta_i_c, delta_j_c, x_refine, y_refine, &
       is_fine_pe, nest_domain, nz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)
    type(IPD_control_type), intent(in)               :: IPD_Control
    type(IPD_data_type), intent(inout)               :: IPD_Data(:)
    integer, intent(in)                              :: n, child_grid_num
    real, allocatable, intent(in)                    :: wt_h(:,:,:), wt_u(:,:,:), wt_v(:,:,:)
    integer, intent(in)                              :: delta_i_c, delta_j_c, x_refine, y_refine   
    logical, intent(in)                              :: is_fine_pe
    type(nest_domain_type), intent(inout)            :: nest_domain
    integer, intent(in)                              :: nz

    ! Constants for mpp calls
    integer  :: interp_type   = 1    ! cell-centered A-grid
    integer  :: interp_type_u = 4    ! D-grid
    integer  :: interp_type_v = 4    ! D-grid
    integer  :: position      = CENTER ! CENTER, NORTH, EAST
    integer  :: position_u    = NORTH
    integer  :: position_v    = EAST
    type(fv_moving_nest_physics_type), pointer :: mn_phys

    mn_phys => Atm(n)%mn_phys

    !! Skin temp/SST
    call mn_var_shift_data(mn_phys%ts, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
         delta_i_c, delta_j_c, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position)

    if (move_physics) then
       !! Soil variables
       call mn_var_shift_data(mn_phys%smc, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position, IPD_Control%lsoil)
       call mn_var_shift_data(mn_phys%stc, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position, IPD_Control%lsoil)
       call mn_var_shift_data(mn_phys%slc, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position, IPD_Control%lsoil)
       
       !! Physics arrays
       call mn_var_shift_data(mn_phys%phy_f2d, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position, IPD_control%ntot2d)
       
       call mn_var_shift_data(mn_phys%phy_f3d, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position, IPD_Control%levs)
       
       
       ! Surface variables
       call mn_var_shift_data(mn_phys%u10m, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
         delta_i_c, delta_j_c, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%v10m, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
         delta_i_c, delta_j_c, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%tprcp, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
         delta_i_c, delta_j_c, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position)

       call mn_var_shift_data(mn_phys%hprime, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
         delta_i_c, delta_j_c, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position, IPD_Control%nmtvr)
       

       call mn_var_shift_data(mn_phys%zorl, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
         delta_i_c, delta_j_c, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%alvsf, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%alvwf, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%alnsf, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%alnwf, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)

       call mn_var_shift_data(mn_phys%facsf, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%facwf, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       
       
       call mn_var_shift_data(mn_phys%canopy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%vegfrac, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%uustar, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%shdmin, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%shdmax, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%zorll, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%zorlo, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%zorlw, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%tsfco, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%tsfcl, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%tsfc, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       
       call mn_var_shift_data(mn_phys%cv, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%cvt, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%cvb, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       
       end if
       
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
          !call mn_var_shift_data(mn_phys%ifd, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          !     delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
          call mn_var_shift_data(mn_phys%dt_cool, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
               delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
          call mn_var_shift_data(mn_phys%qrain, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
               delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)          
       end if

  end subroutine mn_phys_shift_data



  subroutine mn_phys_dump_to_netcdf(Atm, Atm_block, IPD_Control, IPD_Data, time_val, file_prefix, is_fine_pe, domain_coarse, domain_fine, nz)
    type(fv_atmos_type), intent(in)            :: Atm
    type (block_control_type), intent(in)      :: Atm_block
    type(IPD_control_type), intent(in)         :: IPD_Control
    type(IPD_data_type), intent(in)            :: IPD_Data(:)
    integer, intent(in)                        :: time_val
    character(len=*), intent(in)               :: file_prefix
    logical, intent(in)                        :: is_fine_pe
    type(domain2d), intent(in)                 :: domain_coarse, domain_fine
    integer, intent(in)                        :: nz

    integer :: is, ie, js, je
    integer :: nb, blen, i, j, k, ix, nv 
    integer :: this_pe

    integer            :: n_moist
    character(len=16)  :: out_var_name, phys_var_name
    integer            :: position = CENTER 
    !integer            :: position_u = NORTH
    !integer            :: position_v = EAST


    ! TODO PHYSICS Add processing of physics variables, following example in mn_prog_dump_to_netcdf
    !real (kind=kind_phys), allocatable :: smc_pr_local (:,:,:)  !< soil moisture content
    !real (kind=kind_phys), allocatable :: stc_pr_local (:,:,:)  !< soil temperature
    !real (kind=kind_phys), allocatable :: slc_pr_local (:,:,:)  !< soil liquid water content
    

    ! Coerce the high precision variables from physics into regular precision for debugging netCDF output
    ! Does not affect values used in calculations.
    ! TODO do we want to dump these as double precision??
    real, allocatable :: smc_pr_local (:,:,:)  !< soil moisture content
    real, allocatable :: stc_pr_local (:,:,:)  !< soil temperature
    real, allocatable :: slc_pr_local (:,:,:)  !< soil liquid water content

    real, allocatable, dimension(:,:) :: sealand_pr_local, deep_soil_t_pr_local, soil_type_pr_local, veg_type_pr_local, slope_type_pr_local, max_snow_alb_pr_local

    real, allocatable, dimension(:,:) ::  tsfco_pr_local, tsfcl_pr_local, tsfc_pr_local, vegfrac_pr_local, alvsf_pr_local, alvwf_pr_local, facsf_pr_local, facwf_pr_local
    real, allocatable, dimension(:,:) :: tref_pr_local, c_0_pr_local, xt_pr_local,  xu_pr_local,  xv_pr_local, ifd_pr_local

    real, allocatable :: phy_f2d_pr_local (:,:,:)
    real, allocatable :: phy_f3d_pr_local (:,:,:,:)

    this_pe = mpp_pe()

    !  Skin temp/SST
    call mn_var_dump_to_netcdf(Atm%ts   , is_fine_pe, domain_coarse, domain_fine, position, 1, &
         time_val, Atm%global_tile, file_prefix, "SSTK")

    !  Terrain height == phis / grav
    call mn_var_dump_to_netcdf(Atm%phis / grav   , is_fine_pe, domain_coarse, domain_fine, position, 1, &
         time_val, Atm%global_tile, file_prefix, "orog")
    
    ! sgh and oro were only fully allocated if fv_land is True
    !      if false, oro is (1,1), and sgh is not allocated
    if ( Atm%flagstruct%fv_land ) then
       !print '("[INFO] WDR mn_phys_dump_to_netcdf fv_land TRUE npe=",I0," size(oro)=(",I0,",",I0,")")', this_pe, size(Atm%oro, 1), size(Atm%oro, 1)
       !print '("[INFO] WDR mn_phys_dump_to_netcdf fv_land TRUE npe=",I0," size(sgh)=(",I0,",",I0,")")', this_pe, size(Atm%sgh, 1), size(Atm%sgh, 1)
       ! land frac --  called oro in fv_array.F90
       call mn_var_dump_to_netcdf(Atm%oro   , is_fine_pe, domain_coarse, domain_fine, position, 1, &
            time_val, Atm%global_tile, file_prefix, "LFRAC")
       
       ! terrain standard deviation --  called sgh in fv_array.F90
       call mn_var_dump_to_netcdf(Atm%sgh   , is_fine_pe, domain_coarse, domain_fine, position, 1, &
            time_val, Atm%global_tile, file_prefix, "STDDEV")
    !else 
       !print '("[INFO] WDR mn_phys_dump_to_netcdf fv_land FALSE npe=",I0, " size(oro)=(",I0,",",I0,")")', this_pe, size(Atm%oro, 1), size(Atm%oro, 1)
    end if

    ! Latitude and longitude in radians
    !call mn_var_dump_to_netcdf( Atm%gridstruct%agrid(:,:,2), is_fine_pe, domain_coarse, domain_fine, position, nz, &
    !     time_val, Atm%global_tile, file_prefix, "latrad")
    !call mn_var_dump_to_netcdf( Atm%gridstruct%agrid(:,:,1), is_fine_pe, domain_coarse, domain_fine, position, nz, &
    !     time_val, Atm%global_tile, file_prefix, "lonrad")



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
       
       allocate ( phy_f2d_pr_local(is:ie, js:je, IPD_Control%ntot2d) )
       allocate ( phy_f3d_pr_local(is:ie, js:je, IPD_Control%levs, IPD_Control%ntot3d) )
       
       
       allocate ( tsfco_pr_local(is:ie, js:je) )
       allocate ( tsfcl_pr_local(is:ie, js:je) )
       allocate ( tsfc_pr_local(is:ie, js:je) )
       allocate ( vegfrac_pr_local(is:ie, js:je) )
       allocate ( alvsf_pr_local(is:ie, js:je) )
       allocate ( alvwf_pr_local(is:ie, js:je) )
       
       ! Static file variables
       allocate ( deep_soil_t_pr_local(is:ie, js:je) )
       allocate ( soil_type_pr_local(is:ie, js:je) )
       
       !allocate ( veg_frac_pr_local(is:ie, js:je) )
       allocate ( veg_type_pr_local(is:ie, js:je) )
       
       allocate ( slope_type_pr_local(is:ie, js:je) )
       
       allocate ( max_snow_alb_pr_local(is:ie, js:je) )

       allocate ( facsf_pr_local(is:ie, js:je) )
       allocate ( facwf_pr_local(is:ie, js:je) )

    end if

    if (move_nsst) then
       allocate ( tref_pr_local(is:ie, js:je) )
       allocate ( c_0_pr_local(is:ie, js:je) )
       allocate ( xt_pr_local(is:ie, js:je) )
       allocate ( xu_pr_local(is:ie, js:je) )
       allocate ( xv_pr_local(is:ie, js:je) )
       allocate ( ifd_pr_local(is:ie, js:je) )
    end if

    
    if (move_physics) then
       smc_pr_local = +99999.9
       stc_pr_local = +99999.9
       slc_pr_local = +99999.9
       
       sealand_pr_local = +99999.9
       
       phy_f2d_pr_local = +99999.9
       phy_f3d_pr_local = +99999.9
       
       tsfco_pr_local = +99999.9
       tsfcl_pr_local = +99999.9
       tsfc_pr_local = +99999.9
       vegfrac_pr_local = +99999.9
       alvsf_pr_local = +99999.9
       alvwf_pr_local = +99999.9



    end if
    if (move_nsst) then
       tref_pr_local = +99999.9
       c_0_pr_local = +99999.9
       xt_pr_local = +99999.9
       xu_pr_local = +99999.9
       xv_pr_local = +99999.9
       ifd_pr_local = +99999.9
    end if 

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
             
             deep_soil_t_pr_local(i, j) = IPD_data(nb)%Sfcprop%tg3(ix)
             soil_type_pr_local(i, j) = IPD_data(nb)%Sfcprop%stype(ix)
             
             !veg_frac_pr_local(i, j) = IPD_data(nb)%Sfcprop%vfrac(ix)
             veg_type_pr_local(i, j) = IPD_data(nb)%Sfcprop%vtype(ix)
             
             slope_type_pr_local(i, j) = IPD_data(nb)%Sfcprop%slope(ix)
             
             facsf_pr_local(i, j) = IPD_data(nb)%Sfcprop%facsf(ix)
             facwf_pr_local(i, j) = IPD_data(nb)%Sfcprop%facwf(ix)

             max_snow_alb_pr_local(i, j) = IPD_data(nb)%Sfcprop%snoalb(ix)
             
             
             tsfco_pr_local(i, j) = IPD_data(nb)%Sfcprop%tsfco(ix)
             tsfcl_pr_local(i, j) = IPD_data(nb)%Sfcprop%tsfcl(ix)
             tsfc_pr_local(i, j)  = IPD_data(nb)%Sfcprop%tsfc(ix)
             vegfrac_pr_local(i, j) = IPD_data(nb)%Sfcprop%vfrac(ix)
             alvsf_pr_local(i, j) = IPD_data(nb)%Sfcprop%alvsf(ix)
             alvwf_pr_local(i, j) = IPD_data(nb)%Sfcprop%alvwf(ix)
             
             do nv = 1, IPD_Control%ntot2d
                ! Use real() to lower the precision
                phy_f2d_pr_local(i,j,nv) = real(IPD_Data(nb)%Tbd%phy_f2d(ix, nv))
             end do
             
             do k = 1, IPD_Control%levs
                do nv = 1, IPD_Control%ntot3d
                   ! Use real() to lower the precision
                   phy_f3d_pr_local(i,j,k,nv) = real(IPD_Data(nb)%Tbd%phy_f3d(ix, k, nv))
                end do
             end do
          end if
          
          if (move_nsst) then
             tref_pr_local(i, j) = IPD_data(nb)%Sfcprop%tref(ix)
             c_0_pr_local(i, j) = IPD_data(nb)%Sfcprop%c_0(ix)
             xt_pr_local(i, j) = IPD_data(nb)%Sfcprop%xt(ix)
             xu_pr_local(i, j) = IPD_data(nb)%Sfcprop%xu(ix)
             xv_pr_local(i, j) = IPD_data(nb)%Sfcprop%xv(ix)

             ifd_pr_local(i, j) = IPD_data(nb)%Sfcprop%ifd(ix)
             
          end if
          
 
       enddo
    enddo


    call mn_var_dump_to_netcdf(stc_pr_local   , is_fine_pe, domain_coarse, domain_fine, position, IPD_Control%lsoil, &
         time_val, Atm%global_tile, file_prefix, "SOILT")

    if (move_physics) then
       call mn_var_dump_to_netcdf(smc_pr_local   , is_fine_pe, domain_coarse, domain_fine, position, IPD_Control%lsoil, &
            time_val, Atm%global_tile, file_prefix, "SOILM")
       
       call mn_var_dump_to_netcdf(slc_pr_local   , is_fine_pe, domain_coarse, domain_fine, position, IPD_Control%lsoil, &
            time_val, Atm%global_tile, file_prefix, "SOILL")
       
       call mn_var_dump_to_netcdf(sealand_pr_local   , is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "LMASK")
       
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
       
       call mn_var_dump_to_netcdf(facsf_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "FACSF")
       call mn_var_dump_to_netcdf(facwf_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "FACWF")
       
       
       do nv = 1, IPD_Control%ntot2d
          write (phys_var_name, "(A4,I0.3)")  'PH2D', nv
          call mn_var_dump_to_netcdf(phy_f2d_pr_local(:,:,nv), is_fine_pe, domain_coarse, domain_fine, position, 1, &
               time_val, Atm%global_tile, file_prefix, phys_var_name)       
       end do
       
       do nv = 1, IPD_Control%ntot3d
          write (phys_var_name, "(A4,I0.3)")  'PH3D', nv
          call mn_var_dump_to_netcdf(phy_f3d_pr_local(:,:,:,nv), is_fine_pe, domain_coarse, domain_fine, position, IPD_Control%levs, &
               time_val, Atm%global_tile, file_prefix, phys_var_name)       
       end do
    endif
       
    if (move_nsst) then
       call mn_var_dump_to_netcdf(tref_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "TREF")
       call mn_var_dump_to_netcdf(c_0_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "C_0")
       call mn_var_dump_to_netcdf(xt_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "XT")
       call mn_var_dump_to_netcdf(xu_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "XU")
       call mn_var_dump_to_netcdf(xv_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "XV")

       call mn_var_dump_to_netcdf(ifd_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "IFD")

    end if

    ! TODO does skin temp need to be deallocated?
    if (move_physics) then
       deallocate(smc_pr_local)
       deallocate(stc_pr_local)
       deallocate(slc_pr_local)
       
       deallocate(sealand_pr_local, deep_soil_t_pr_local, soil_type_pr_local, veg_type_pr_local, max_snow_alb_pr_local)
       
       deallocate(tsfco_pr_local, tsfcl_pr_local, tsfc_pr_local, vegfrac_pr_local, alvsf_pr_local, alvwf_pr_local)
       deallocate(facsf_pr_local, facwf_pr_local)

       deallocate(phy_f2d_pr_local)
       deallocate(phy_f3d_pr_local)
    end if
    if (move_nsst) deallocate(tref_pr_local, c_0_pr_local, xt_pr_local,  xu_pr_local,  xv_pr_local, ifd_pr_local)

    if (debug_log) print '("[INFO] WDR end mn_phys_dump_tp_netcdf npe=",I0)', this_pe


  end subroutine mn_phys_dump_to_netcdf

#endif MOVING_NEST

end module fv_moving_nest_physics_mod
