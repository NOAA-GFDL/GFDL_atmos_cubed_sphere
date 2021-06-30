!***********************************************************************
!>@brief!   Provides subroutines to enable moving nest functionality in FV3 dynamic core.  
!!>@author Bill Ramstrom, AOML/HRD   01/15/2021
!
! =======================================================================!
!


! Notes


module fv_moving_nest_utils_mod

#ifdef MOVING_NEST

  use mpp_mod,         only : FATAL, WARNING, MPP_DEBUG, NOTE, MPP_CLOCK_SYNC,MPP_CLOCK_DETAILED
  use mpp_mod,         only : mpp_pe, mpp_npes, mpp_root_pe, mpp_error, mpp_set_warn_level
  use mpp_mod,         only : mpp_declare_pelist, mpp_set_current_pelist, mpp_sync, mpp_sync_self
  use mpp_mod,         only : mpp_clock_begin, mpp_clock_end, mpp_clock_id
  use mpp_mod,         only : mpp_init, mpp_exit, mpp_chksum, stdout, stderr
  use mpp_mod,         only : input_nml_file
  use mpp_mod,         only : mpp_get_current_pelist, mpp_broadcast
  use mpp_domains_mod, only : GLOBAL_DATA_DOMAIN, BITWISE_EXACT_SUM, BGRID_NE, CGRID_NE, DGRID_NE, AGRID
  use mpp_parameter_mod, only: AGRID_PARAM=>AGRID,CGRID_NE_PARAM=>CGRID_NE,SCALAR_PAIR
  use mpp_domains_mod, only : FOLD_SOUTH_EDGE, FOLD_NORTH_EDGE, FOLD_WEST_EDGE, FOLD_EAST_EDGE
  use mpp_domains_mod, only : MPP_DOMAIN_TIME, CYCLIC_GLOBAL_DOMAIN, NUPDATE,EUPDATE, XUPDATE, YUPDATE, SCALAR_PAIR
  use mpp_domains_mod, only : domain1D, domain2D, DomainCommunicator2D, BITWISE_EFP_SUM
  use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain, mpp_domains_set_stack_size
  use mpp_domains_mod, only : mpp_global_field, mpp_global_sum, mpp_global_max, mpp_global_min
  use mpp_domains_mod, only : mpp_domains_init, mpp_domains_exit, mpp_broadcast_domain
  use mpp_domains_mod, only : mpp_update_domains, mpp_check_field, mpp_redistribute, mpp_get_memory_domain
  use mpp_domains_mod, only : mpp_define_layout, mpp_define_domains, mpp_modify_domain
  use mpp_domains_mod, only : mpp_define_io_domain
  use mpp_domains_mod, only : mpp_get_neighbor_pe, mpp_define_mosaic, mpp_nullify_domain_list
  use mpp_domains_mod, only : NORTH, NORTH_EAST, EAST, SOUTH_EAST, CORNER, CENTER
  use mpp_domains_mod, only : SOUTH, SOUTH_WEST, WEST, NORTH_WEST, mpp_define_mosaic_pelist
  use mpp_domains_mod, only : mpp_get_global_domain, ZERO, NINETY, MINUS_NINETY
  use mpp_domains_mod, only : mpp_get_boundary, mpp_start_update_domains, mpp_complete_update_domains
  use mpp_domains_mod, only : mpp_define_nest_domains, nest_domain_type
  use mpp_domains_mod, only : mpp_get_C2F_index, mpp_update_nest_fine
  use mpp_domains_mod, only : mpp_get_F2C_index, mpp_update_nest_coarse
  use mpp_domains_mod, only : mpp_get_domain_shift, EDGEUPDATE, mpp_deallocate_domain
  use mpp_domains_mod, only : mpp_group_update_type, mpp_create_group_update
  use mpp_domains_mod, only : mpp_do_group_update, mpp_clear_group_update
  use mpp_domains_mod, only : mpp_start_group_update, mpp_complete_group_update
  use mpp_domains_mod, only : WUPDATE, SUPDATE, mpp_get_compute_domains, NONSYMEDGEUPDATE
  use mpp_domains_mod, only : domainUG, mpp_define_unstruct_domain, mpp_get_UG_domain_tile_id
  use mpp_domains_mod, only : mpp_get_UG_compute_domain, mpp_pass_SG_to_UG, mpp_pass_UG_to_SG
  use mpp_domains_mod, only : mpp_get_ug_global_domain, mpp_global_field_ug
  use mpp_memutils_mod, only : mpp_memuse_begin, mpp_memuse_end

#ifdef GFS_TYPES
use GFS_typedefs,           only: kind_phys
#else
use IPD_typedefs,           only: kind_phys => IPD_kind_phys
#endif




  ! Added WDR
  use boundary_mod,      only: update_coarse_grid, update_coarse_grid_mpp
  use bounding_box_mod,  only: bbox, bbox_get_C2F_index, fill_bbox, show_bbox
  use fms_io_mod,        only: read_data, write_data, get_global_att_value, fms_io_init, fms_io_exit
  use fv_arrays_mod,     only: R_GRID
  use fv_arrays_mod,     only: fv_atmos_type

  !use bounding_box,              only: bbox
  !use fv_moving_nest,     only: bbox, grid_geometry, fill_nest_from_buffer
  !use fv_moving_nest,     only: grid_geometry
  implicit none

#ifdef NO_QUAD_PRECISION
  ! 64-bit precision (kind=8)
  integer, parameter:: f_p = selected_real_kind(15)
#else
  ! Higher precision (kind=16) for grid geometrical factors:
  integer, parameter:: f_p = selected_real_kind(20)
#endif

  integer, parameter :: UWIND = 1
  integer, parameter :: VWIND = 2

  logical :: debug_log = .false.
  

#include <fms_platform.h>

  ! Encapsulates the grid definition data, such as read from the netCDF files
  type grid_geometry
     integer   :: nx, ny, nxp, nyp

     real(kind=kind_phys), allocatable  :: lats(:,:)
     real(kind=kind_phys), allocatable  :: lons(:,:)

     real, allocatable  :: dx(:,:)
     real, allocatable  :: dy(:,:)
     real(kind=kind_phys), allocatable  :: area(:,:)
  end type grid_geometry


  interface alloc_read_data
     module procedure alloc_read_data_dyn
     module procedure alloc_read_data_kind_phys
  end interface alloc_read_data

  interface fill_nest_halos_from_parent
     module procedure fill_nest_halos_from_parent2D
     module procedure fill_nest_halos_from_parent3D
     module procedure fill_nest_halos_from_parent3D_kindphys
     module procedure fill_nest_halos_from_parent4D
     module procedure fill_nest_halos_from_parent4D_kindphys
  end interface fill_nest_halos_from_parent
  
  interface alloc_halo_buffer
     module procedure alloc_2D_halo_buffer
     module procedure alloc_3D_halo_buffer
     module procedure alloc_3D_halo_buffer_kindphys
     module procedure alloc_4D_halo_buffer
     module procedure alloc_4D_halo_buffer_kindphys
  end interface alloc_halo_buffer

  interface fill_nest_from_buffer
     module procedure fill_nest_from_buffer2D
     module procedure fill_nest_from_buffer3D
     module procedure fill_nest_from_buffer3D_kindphys
     module procedure fill_nest_from_buffer4D
     module procedure fill_nest_from_buffer4D_kindphys
  end interface fill_nest_from_buffer

  interface fill_nest_from_buffer_cell_center
     module procedure fill_nest_from_buffer_cell_center2D
     module procedure fill_nest_from_buffer_cell_center3D
     module procedure fill_nest_from_buffer_cell_center3D_kindphys
     module procedure fill_nest_from_buffer_cell_center4D
     module procedure fill_nest_from_buffer_cell_center4D_kindphys
  end interface fill_nest_from_buffer_cell_center

  interface output_grid_to_nc
     module procedure output_grid_2d_to_nc
     module procedure output_grid_3d_to_nc
  end interface output_grid_to_nc


  interface fill_grid_from_supergrid
     module procedure fill_grid32_from_supergrid
     module procedure fill_grid64_from_supergrid
     module procedure fill_grid64_4D_from_supergrid
  end interface fill_grid_from_supergrid

contains


  !==================================================================================================
  !
  !  Fill Nest Halos from Parent
  !
  !==================================================================================================



subroutine fill_nest_halos_from_parent2D(var_name, data_var, interp_type, wt, ind, x_refine, y_refine, is_fine_pe, nest_domain, position)
  character(len=*), intent(in)                :: var_name
  real, allocatable, intent(inout)            :: data_var(:,:)
  integer, intent(in)                         :: interp_type
  real, allocatable, intent(in)               :: wt(:,:,:)
  integer, allocatable, intent(in)            :: ind(:,:,:)
  integer, intent(in)                         :: x_refine, y_refine
  logical, intent(in)                         :: is_fine_pe
  type(nest_domain_type), intent(inout)       :: nest_domain
  integer, intent(in)                         :: position


  real, dimension(:,:), allocatable   :: nbuffer, sbuffer, ebuffer, wbuffer
  type(bbox)                          :: north_fine, north_coarse
  type(bbox)                          :: south_fine, south_coarse
  type(bbox)                          :: east_fine, east_coarse
  type(bbox)                          :: west_fine, west_coarse
  integer                             :: this_pe
  integer                             :: nest_level = 1  ! WDR TODO allow to vary

  this_pe = mpp_pe()

  !!===========================================================
  !!
  !! Fill halo buffers
  !!
  !!===========================================================

  if (debug_log) then

     print '("[INFO] WDR Start fill_nest_halos_from_parent2D. npe=",I0," var_name=",A16)', this_pe, var_name
     print '("[INFO] fill_nest_halos npe=",I0," nest_domain%tile_fine=",I0," %tile_coarse=",I0)', this_pe, nest_domain%tile_fine, nest_domain%tile_coarse
     print '("[INFO] fill_nest_halos npe=",I0," nest_domain%istart_fine=",I0," %iend_fine=",I0)', this_pe, nest_domain%istart_fine,  nest_domain%iend_fine
     print '("[INFO] fill_nest_halos npe=",I0," nest_domain%jstart_fine=",I0," %jend_fine=",I0)', this_pe, nest_domain%jstart_fine,  nest_domain%jend_fine
     print '("[INFO] fill_nest_halos npe=",I0," nest_domain%istart_coarse=",I0," %iend_coarse=",I0)', this_pe, nest_domain%istart_coarse,  nest_domain%iend_coarse
     print '("[INFO] fill_nest_halos npe=",I0," nest_domain%jstart_coarse=",I0," %jend_coarse=",I0)', this_pe, nest_domain%jstart_coarse,  nest_domain%jend_coarse
     print '("[INFO] data_var npe=",I0," var_name=",A16," data_var(",I0,"-",I0,",",I0,"-",I0,")")', &
          this_pe,  var_name, lbound(data_var, 1), ubound(data_var, 1), lbound(data_var, 2), ubound(data_var, 2)
     print '("[INFO] wt npe=",I0," var_name=",A16," wt(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', &
          this_pe,  var_name, lbound(wt, 1), ubound(wt, 1), lbound(wt, 2), ubound(wt, 2), lbound(wt, 3), ubound(wt, 3)
  end if

  !====================================================

  if (debug_log) print '("[INFO] WDR ALL1. npe=",I0," position=",I0)', this_pe, position

  call alloc_halo_buffer(nbuffer, north_fine, north_coarse, nest_domain, NORTH,  position)
  call alloc_halo_buffer(sbuffer, south_fine, south_coarse, nest_domain, SOUTH,  position)
  call alloc_halo_buffer(ebuffer, east_fine,  east_coarse,  nest_domain, EAST,   position)
  call alloc_halo_buffer(wbuffer, west_fine,  west_coarse,  nest_domain, WEST,   position)

  if (debug_log) print '("[INFO] WDR NRF0.d mn_var_shift_data npe=",I0," data_var(",I0,",",I0,")")', this_pe, size(data_var,1), size(data_var,2)

  ! Passes data from coarse grid to fine grid's halo
  call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

  if (is_fine_pe) then

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

  end if


  deallocate(nbuffer)
  deallocate(sbuffer)
  deallocate(ebuffer)
  deallocate(wbuffer)

  if (debug_log) print '("[INFO] WDR End fill_nest_halos_from_parent2D. npe=",I0," var_name=",A16)', this_pe, var_name

end subroutine fill_nest_halos_from_parent2D


subroutine fill_nest_halos_from_parent3D(var_name, data_var, interp_type, wt, ind, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)
  character(len=*), intent(in)                :: var_name
  real, allocatable, intent(inout)            :: data_var(:,:,:)
  integer, intent(in)                         :: interp_type
  real, allocatable, intent(in)               :: wt(:,:,:)
  integer, allocatable, intent(in)            :: ind(:,:,:)
  integer, intent(in)                         :: x_refine, y_refine
  logical, intent(in)                         :: is_fine_pe
  type(nest_domain_type), intent(inout)       :: nest_domain
  integer, intent(in)                         :: position, nz


  real, dimension(:,:,:), allocatable :: nbuffer, sbuffer, ebuffer, wbuffer
  type(bbox)                          :: north_fine, north_coarse
  type(bbox)                          :: south_fine, south_coarse
  type(bbox)                          :: east_fine, east_coarse
  type(bbox)                          :: west_fine, west_coarse
  integer                             :: this_pe
  integer                             :: nest_level = 1  ! WDR TODO allow to vary

  this_pe = mpp_pe()

  !!===========================================================
  !!
  !! Fill halo buffers
  !!
  !!===========================================================

  if (debug_log) then

     print '("[INFO] WDR Start fill_nest_halos_from_parent3D. npe=",I0," var_name=",A16)', this_pe, var_name
     print '("[INFO] fill_nest_halos npe=",I0," nest_domain%tile_fine=",I0," %tile_coarse=",I0)', this_pe, nest_domain%tile_fine, nest_domain%tile_coarse
     print '("[INFO] fill_nest_halos npe=",I0," nest_domain%istart_fine=",I0," %iend_fine=",I0)', this_pe, nest_domain%istart_fine,  nest_domain%iend_fine
     print '("[INFO] fill_nest_halos npe=",I0," nest_domain%jstart_fine=",I0," %jend_fine=",I0)', this_pe, nest_domain%jstart_fine,  nest_domain%jend_fine
     print '("[INFO] fill_nest_halos npe=",I0," nest_domain%istart_coarse=",I0," %iend_coarse=",I0)', this_pe, nest_domain%istart_coarse,  nest_domain%iend_coarse
     print '("[INFO] fill_nest_halos npe=",I0," nest_domain%jstart_coarse=",I0," %jend_coarse=",I0)', this_pe, nest_domain%jstart_coarse,  nest_domain%jend_coarse
     print '("[INFO] data_var npe=",I0," var_name=",A16," data_var(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', &
          this_pe,  var_name, lbound(data_var, 1), ubound(data_var, 1), lbound(data_var, 2), ubound(data_var, 2), lbound(data_var, 3), ubound(data_var, 3)
     print '("[INFO] wt npe=",I0," var_name=",A16," wt(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', &
          this_pe,  var_name, lbound(wt, 1), ubound(wt, 1), lbound(wt, 2), ubound(wt, 2), lbound(wt, 3), ubound(wt, 3)
  end if

  !====================================================

  if (debug_log) print '("[INFO] WDR ALL1. npe=",I0," position=",I0," nz=",I0)', this_pe, position, nz

  call alloc_halo_buffer(nbuffer, north_fine, north_coarse, nest_domain, NORTH,  position, nz)
  call alloc_halo_buffer(sbuffer, south_fine, south_coarse, nest_domain, SOUTH,  position, nz)
  call alloc_halo_buffer(ebuffer, east_fine,  east_coarse,  nest_domain, EAST,   position, nz)
  call alloc_halo_buffer(wbuffer, west_fine,  west_coarse,  nest_domain, WEST,   position, nz)

  if (debug_log) print '("[INFO] WDR NRF0.d mn_var_shift_data npe=",I0," data_var(",I0,",",I0,",",I0,")")', this_pe, size(data_var,1), size(data_var,2), size(data_var,3)

  ! Passes data from coarse grid to fine grid's halo
  call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

  if (is_fine_pe) then

     !!===========================================================
     !!
     !! Apply halo data
     !!
     !!===========================================================

     if (debug_log) print '("[INFO] WDR NRFI mn_var_shift_data start. npe=",I0)', this_pe

     call fill_nest_from_buffer(interp_type, data_var, nbuffer, north_fine, north_coarse, nz, NORTH, x_refine, y_refine, wt, ind)
     if (debug_log) print '("[INFO] WDR NRF N mn_var_shift_data start. npe=",I0)', this_pe

     call fill_nest_from_buffer(interp_type, data_var, sbuffer, south_fine, south_coarse, nz, SOUTH, x_refine, y_refine, wt, ind)
     if (debug_log) print '("[INFO] WDR NRF S mn_var_shift_data start. npe=",I0)', this_pe

     call fill_nest_from_buffer(interp_type, data_var, ebuffer, east_fine, east_coarse, nz, EAST, x_refine, y_refine, wt, ind)
     if (debug_log) print '("[INFO] WDR NRF E mn_var_shift_data start. npe=",I0)', this_pe

     call fill_nest_from_buffer(interp_type, data_var, wbuffer, west_fine, west_coarse, nz, WEST, x_refine, y_refine, wt, ind)
     if (debug_log) print '("[INFO] WDR NRF W mn_var_shift_data start. npe=",I0)', this_pe

  end if


  deallocate(nbuffer)
  deallocate(sbuffer)
  deallocate(ebuffer)
  deallocate(wbuffer)

  if (debug_log) print '("[INFO] WDR End fill_nest_halos_from_parent3D. npe=",I0," var_name=",A16)', this_pe, var_name

end subroutine fill_nest_halos_from_parent3D


subroutine fill_nest_halos_from_parent3D_kindphys(var_name, data_var, interp_type, wt, ind, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)
  character(len=*), intent(in)                :: var_name
  real(kind=kind_phys), allocatable, intent(inout)            :: data_var(:,:,:)
  integer, intent(in)                         :: interp_type
  real, allocatable, intent(in)               :: wt(:,:,:)
  integer, allocatable, intent(in)            :: ind(:,:,:)
  integer, intent(in)                         :: x_refine, y_refine
  logical, intent(in)                         :: is_fine_pe
  type(nest_domain_type), intent(inout)       :: nest_domain
  integer, intent(in)                         :: position, nz


  real(kind=kind_phys), dimension(:,:,:), allocatable :: nbuffer, sbuffer, ebuffer, wbuffer
  type(bbox)                          :: north_fine, north_coarse
  type(bbox)                          :: south_fine, south_coarse
  type(bbox)                          :: east_fine, east_coarse
  type(bbox)                          :: west_fine, west_coarse
  integer                             :: this_pe
  integer                             :: nest_level = 1  ! WDR TODO allow to vary

  this_pe = mpp_pe()

  !!===========================================================
  !!
  !! Fill halo buffers
  !!
  !!===========================================================

  if (debug_log) then

     print '("[INFO] WDR Start fill_nest_halos_from_parent3D. npe=",I0," var_name=",A16)', this_pe, var_name
     print '("[INFO] fill_nest_halos npe=",I0," nest_domain%tile_fine=",I0," %tile_coarse=",I0)', this_pe, nest_domain%tile_fine, nest_domain%tile_coarse
     print '("[INFO] fill_nest_halos npe=",I0," nest_domain%istart_fine=",I0," %iend_fine=",I0)', this_pe, nest_domain%istart_fine,  nest_domain%iend_fine
     print '("[INFO] fill_nest_halos npe=",I0," nest_domain%jstart_fine=",I0," %jend_fine=",I0)', this_pe, nest_domain%jstart_fine,  nest_domain%jend_fine
     print '("[INFO] fill_nest_halos npe=",I0," nest_domain%istart_coarse=",I0," %iend_coarse=",I0)', this_pe, nest_domain%istart_coarse,  nest_domain%iend_coarse
     print '("[INFO] fill_nest_halos npe=",I0," nest_domain%jstart_coarse=",I0," %jend_coarse=",I0)', this_pe, nest_domain%jstart_coarse,  nest_domain%jend_coarse
     print '("[INFO] data_var npe=",I0," var_name=",A16," data_var(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', &
          this_pe,  var_name, lbound(data_var, 1), ubound(data_var, 1), lbound(data_var, 2), ubound(data_var, 2), lbound(data_var, 3), ubound(data_var, 3)
     print '("[INFO] wt npe=",I0," var_name=",A16," wt(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', &
          this_pe,  var_name, lbound(wt, 1), ubound(wt, 1), lbound(wt, 2), ubound(wt, 2), lbound(wt, 3), ubound(wt, 3)
  end if

  !====================================================

  if (debug_log) print '("[INFO] WDR ALL1. npe=",I0," position=",I0," nz=",I0)', this_pe, position, nz

  call alloc_halo_buffer(nbuffer, north_fine, north_coarse, nest_domain, NORTH,  position, nz)
  call alloc_halo_buffer(sbuffer, south_fine, south_coarse, nest_domain, SOUTH,  position, nz)
  call alloc_halo_buffer(ebuffer, east_fine,  east_coarse,  nest_domain, EAST,   position, nz)
  call alloc_halo_buffer(wbuffer, west_fine,  west_coarse,  nest_domain, WEST,   position, nz)

  if (debug_log) print '("[INFO] WDR NRF0.d mn_var_shift_data npe=",I0," data_var(",I0,",",I0,",",I0,")")', this_pe, size(data_var,1), size(data_var,2), size(data_var,3)

  ! Passes data from coarse grid to fine grid's halo
  call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

  if (is_fine_pe) then

     !!===========================================================
     !!
     !! Apply halo data
     !!
     !!===========================================================

     if (debug_log) print '("[INFO] WDR NRFI mn_var_shift_data start. npe=",I0)', this_pe

     call fill_nest_from_buffer(interp_type, data_var, nbuffer, north_fine, north_coarse, nz, NORTH, x_refine, y_refine, wt, ind)
     if (debug_log) print '("[INFO] WDR NRF N mn_var_shift_data start. npe=",I0)', this_pe

     call fill_nest_from_buffer(interp_type, data_var, sbuffer, south_fine, south_coarse, nz, SOUTH, x_refine, y_refine, wt, ind)
     if (debug_log) print '("[INFO] WDR NRF S mn_var_shift_data start. npe=",I0)', this_pe

     call fill_nest_from_buffer(interp_type, data_var, ebuffer, east_fine, east_coarse, nz, EAST, x_refine, y_refine, wt, ind)
     if (debug_log) print '("[INFO] WDR NRF E mn_var_shift_data start. npe=",I0)', this_pe

     call fill_nest_from_buffer(interp_type, data_var, wbuffer, west_fine, west_coarse, nz, WEST, x_refine, y_refine, wt, ind)
     if (debug_log) print '("[INFO] WDR NRF W mn_var_shift_data start. npe=",I0)', this_pe

  end if


  deallocate(nbuffer)
  deallocate(sbuffer)
  deallocate(ebuffer)
  deallocate(wbuffer)

  if (debug_log) print '("[INFO] WDR End fill_nest_halos_from_parent3D. npe=",I0," var_name=",A16)', this_pe, var_name

end subroutine fill_nest_halos_from_parent3D_kindphys


subroutine fill_nest_halos_from_parent4D(var_name, data_var, interp_type, wt, ind, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)
  character(len=*), intent(in)                :: var_name
  real, allocatable, intent(inout)            :: data_var(:,:,:,:)
  integer, intent(in)                         :: interp_type
  real, allocatable, intent(in)               :: wt(:,:,:)
  integer, allocatable, intent(in)            :: ind(:,:,:)
  integer, intent(in)                         :: x_refine, y_refine
  logical, intent(in)                         :: is_fine_pe
  type(nest_domain_type), intent(inout)       :: nest_domain
  integer, intent(in)                         :: position, nz

  real, dimension(:,:,:,:), allocatable :: nbuffer, sbuffer, ebuffer, wbuffer
  type(bbox)                            :: north_fine, north_coarse 
  type(bbox)                            :: south_fine, south_coarse
  type(bbox)                            :: east_fine, east_coarse
  type(bbox)                            :: west_fine, west_coarse
  integer                               :: n4d, this_pe
  integer                               :: nest_level = 1  ! WDR TODO allow to vary

  this_pe = mpp_pe()

  !!===========================================================
  !!
  !! Fill halo buffers
  !!
  !!===========================================================

  if (debug_log) print '("[INFO] WDR Start fill_nest_halos_from_parent4D. npe=",I0," var_name=",A16)', this_pe, var_name

  n4d = ubound(data_var, 4)

  if (debug_log) print '("[INFO] fill_nest_halos npe=",I0," nest_domain%tile_fine=",I0," %tile_coarse=",I0)', this_pe, nest_domain%tile_fine, nest_domain%tile_coarse

  if (debug_log) print '("[INFO] fill_nest_halos npe=",I0," nest_domain%istart_fine=",I0," %iend_fine=",I0)', this_pe, nest_domain%istart_fine,  nest_domain%iend_fine
  if (debug_log) print '("[INFO] fill_nest_halos npe=",I0," nest_domain%jstart_fine=",I0," %jend_fine=",I0)', this_pe, nest_domain%jstart_fine,  nest_domain%jend_fine

  if (debug_log) print '("[INFO] fill_nest_halos npe=",I0," nest_domain%istart_coarse=",I0," %iend_coarse=",I0)', this_pe, nest_domain%istart_coarse,  nest_domain%iend_coarse
  if (debug_log) print '("[INFO] fill_nest_halos npe=",I0," nest_domain%jstart_coarse=",I0," %jend_coarse=",I0)', this_pe, nest_domain%jstart_coarse,  nest_domain%jend_coarse

  if (debug_log) print '("[INFO] data_var 4D npe=",I0," var_name=",A16," data_var(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', this_pe,  var_name, lbound(data_var, 1), ubound(data_var, 1), lbound(data_var, 2), ubound(data_var, 2), lbound(data_var, 3), ubound(data_var, 3), lbound(data_var, 4), ubound(data_var, 4)


  if (debug_log) print '("[INFO] wt npe=",I0," var_name=",A16," wt(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', this_pe, var_name, lbound(wt, 1), ubound(wt, 1), lbound(wt, 2), ubound(wt, 2), lbound(wt, 3), ubound(wt, 3)


  !====================================================

  if (debug_log) print '("[INFO] WDR ALL1. npe=",I0," position=",I0," nz=",I0)', this_pe, position, nz

  call alloc_halo_buffer(nbuffer, north_fine, north_coarse, nest_domain, NORTH,  position, nz, n4d)
  call alloc_halo_buffer(sbuffer, south_fine, south_coarse, nest_domain, SOUTH,  position, nz, n4d)
  call alloc_halo_buffer(ebuffer, east_fine,  east_coarse,  nest_domain, EAST,   position, nz, n4d)
  call alloc_halo_buffer(wbuffer, west_fine,  west_coarse,  nest_domain, WEST,   position, nz, n4d)

  if (debug_log) print '("[INFO] WDR allocate_halo_buffers DONE. npe=",I0)', this_pe

  !====================================================

  if (debug_log) print '("[INFO] WDR NRF0.d mn_var_shift_data npe=",I0," data_var(",I0,",",I0,",",I0,")")', this_pe, size(data_var,1), size(data_var,2), size(data_var,3)

  ! Passes data from coarse grid to fine grid's halo
  ! Coarse parent PEs send data from data_var
  ! Fine halo PEs receive data into one or more of the halo buffers
  call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

  if (debug_log) print '("[INFO] WDR NRF2 mn_var_shift_data start. npe=",I0)', this_pe

  if (is_fine_pe) then

     !!===========================================================
     !!
     !! Apply halo data
     !!
     !!===========================================================

     if (debug_log) print '("[INFO] WDR NRFI mn_var_shift_data start. npe=",I0)', this_pe

     call fill_nest_from_buffer(interp_type, data_var, nbuffer, north_fine, north_coarse, nz, NORTH, x_refine, y_refine, wt, ind)
     if (debug_log) print '("[INFO] WDR NRF N mn_var_shift_data start. npe=",I0)', this_pe

     call fill_nest_from_buffer(interp_type, data_var, sbuffer, south_fine, south_coarse, nz, SOUTH, x_refine, y_refine, wt, ind)
     if (debug_log) print '("[INFO] WDR NRF S mn_var_shift_data start. npe=",I0)', this_pe

     call fill_nest_from_buffer(interp_type, data_var, ebuffer, east_fine, east_coarse, nz, EAST, x_refine, y_refine, wt, ind)
     if (debug_log) print '("[INFO] WDR NRF E mn_var_shift_data start. npe=",I0)', this_pe

     call fill_nest_from_buffer(interp_type, data_var, wbuffer, west_fine, west_coarse, nz, WEST, x_refine, y_refine, wt, ind)
     if (debug_log) print '("[INFO] WDR NRF W mn_var_shift_data start. npe=",I0)', this_pe

  end if

  deallocate(nbuffer)
  deallocate(sbuffer)
  deallocate(ebuffer)
  deallocate(wbuffer)

  if (debug_log) print '("[INFO] WDR End fill_nest_halos_from_parent4D. npe=",I0," var_name=",A16)', this_pe, var_name

end subroutine fill_nest_halos_from_parent4D


subroutine fill_nest_halos_from_parent4D_kindphys(var_name, data_var, interp_type, wt, ind, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)
  character(len=*), intent(in)                :: var_name
  real(kind=kind_phys), allocatable, intent(inout)            :: data_var(:,:,:,:)
  integer, intent(in)                         :: interp_type
  real, allocatable, intent(in)               :: wt(:,:,:)
  integer, allocatable, intent(in)            :: ind(:,:,:)
  integer, intent(in)                         :: x_refine, y_refine
  logical, intent(in)                         :: is_fine_pe
  type(nest_domain_type), intent(inout)       :: nest_domain
  integer, intent(in)                         :: position, nz

  real(kind=kind_phys), dimension(:,:,:,:), allocatable :: nbuffer, sbuffer, ebuffer, wbuffer
  type(bbox)                            :: north_fine, north_coarse 
  type(bbox)                            :: south_fine, south_coarse
  type(bbox)                            :: east_fine, east_coarse
  type(bbox)                            :: west_fine, west_coarse
  integer                               :: n4d, this_pe
  integer                               :: nest_level = 1  ! WDR TODO allow to vary

  this_pe = mpp_pe()

  !!===========================================================
  !!
  !! Fill halo buffers
  !!
  !!===========================================================

  if (debug_log) print '("[INFO] WDR Start fill_nest_halos_from_parent4D_kindphys. npe=",I0," var_name=",A16)', this_pe, var_name

  n4d = ubound(data_var, 4)

  if (debug_log) print '("[INFO] fill_nest_halos npe=",I0," nest_domain%tile_fine=",I0," %tile_coarse=",I0)', this_pe, nest_domain%tile_fine, nest_domain%tile_coarse

  if (debug_log) print '("[INFO] fill_nest_halos npe=",I0," nest_domain%istart_fine=",I0," %iend_fine=",I0)', this_pe, nest_domain%istart_fine,  nest_domain%iend_fine
  if (debug_log) print '("[INFO] fill_nest_halos npe=",I0," nest_domain%jstart_fine=",I0," %jend_fine=",I0)', this_pe, nest_domain%jstart_fine,  nest_domain%jend_fine

  if (debug_log) print '("[INFO] fill_nest_halos npe=",I0," nest_domain%istart_coarse=",I0," %iend_coarse=",I0)', this_pe, nest_domain%istart_coarse,  nest_domain%iend_coarse
  if (debug_log) print '("[INFO] fill_nest_halos npe=",I0," nest_domain%jstart_coarse=",I0," %jend_coarse=",I0)', this_pe, nest_domain%jstart_coarse,  nest_domain%jend_coarse

  if (debug_log) print '("[INFO] data_var 4D npe=",I0," var_name=",A16," data_var(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', this_pe,  var_name, lbound(data_var, 1), ubound(data_var, 1), lbound(data_var, 2), ubound(data_var, 2), lbound(data_var, 3), ubound(data_var, 3), lbound(data_var, 4), ubound(data_var, 4)


  if (debug_log) print '("[INFO] wt npe=",I0," var_name=",A16," wt(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', this_pe, var_name, lbound(wt, 1), ubound(wt, 1), lbound(wt, 2), ubound(wt, 2), lbound(wt, 3), ubound(wt, 3)


  !====================================================

  if (debug_log) print '("[INFO] WDR ALL1. npe=",I0," position=",I0," nz=",I0)', this_pe, position, nz

  call alloc_halo_buffer(nbuffer, north_fine, north_coarse, nest_domain, NORTH,  position, nz, n4d)
  call alloc_halo_buffer(sbuffer, south_fine, south_coarse, nest_domain, SOUTH,  position, nz, n4d)
  call alloc_halo_buffer(ebuffer, east_fine,  east_coarse,  nest_domain, EAST,   position, nz, n4d)
  call alloc_halo_buffer(wbuffer, west_fine,  west_coarse,  nest_domain, WEST,   position, nz, n4d)

  if (debug_log) print '("[INFO] WDR allocate_halo_buffers DONE. npe=",I0)', this_pe

  !====================================================

  if (debug_log) print '("[INFO] WDR NRF0.d mn_var_shift_data npe=",I0," data_var(",I0,",",I0,",",I0,")")', this_pe, size(data_var,1), size(data_var,2), size(data_var,3)

  ! Passes data from coarse grid to fine grid's halo
  ! Coarse parent PEs send data from data_var
  ! Fine halo PEs receive data into one or more of the halo buffers
  call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

  if (debug_log) print '("[INFO] WDR NRF2 mn_var_shift_data start. npe=",I0)', this_pe

  if (is_fine_pe) then

     !!===========================================================
     !!
     !! Apply halo data
     !!
     !!===========================================================

     if (debug_log) print '("[INFO] WDR NRFI mn_var_shift_data start. npe=",I0)', this_pe

     call fill_nest_from_buffer(interp_type, data_var, nbuffer, north_fine, north_coarse, nz, NORTH, x_refine, y_refine, wt, ind)
     if (debug_log) print '("[INFO] WDR NRF N mn_var_shift_data start. npe=",I0)', this_pe

     call fill_nest_from_buffer(interp_type, data_var, sbuffer, south_fine, south_coarse, nz, SOUTH, x_refine, y_refine, wt, ind)
     if (debug_log) print '("[INFO] WDR NRF S mn_var_shift_data start. npe=",I0)', this_pe

     call fill_nest_from_buffer(interp_type, data_var, ebuffer, east_fine, east_coarse, nz, EAST, x_refine, y_refine, wt, ind)
     if (debug_log) print '("[INFO] WDR NRF E mn_var_shift_data start. npe=",I0)', this_pe

     call fill_nest_from_buffer(interp_type, data_var, wbuffer, west_fine, west_coarse, nz, WEST, x_refine, y_refine, wt, ind)
     if (debug_log) print '("[INFO] WDR NRF W mn_var_shift_data start. npe=",I0)', this_pe

  end if

  deallocate(nbuffer)
  deallocate(sbuffer)
  deallocate(ebuffer)
  deallocate(wbuffer)

  if (debug_log) print '("[INFO] WDR End fill_nest_halos_from_parent4D_kindphys. npe=",I0," var_name=",A16)', this_pe, var_name

end subroutine fill_nest_halos_from_parent4D_kindphys




  !==================================================================================================
  !
  !  Allocate halo buffers
  !
  !==================================================================================================


  subroutine alloc_2D_halo_buffer(buffer, bbox_fine, bbox_coarse, nest_domain, direction, position)
    real, dimension(:,:), allocatable, intent(out)   :: buffer
    type(bbox), intent(out)                          :: bbox_fine, bbox_coarse 
    type(nest_domain_type), intent(inout)            :: nest_domain
    integer, intent(in)                              :: direction, position


    integer                             :: my_stat
    character(256)                      :: my_errmsg
    integer                             :: this_pe

    this_pe = mpp_pe()

    call bbox_get_C2F_index(nest_domain, bbox_fine, bbox_coarse, direction,  position)
    if (debug_log) print '("[INFO] WDR FNHC npe=",I0," direction=",I0," bbox_coarse(",I0,"-",I0,",",I0,"-",I0,")")', this_pe, direction, bbox_coarse.is, bbox_coarse.ie, bbox_coarse.js, bbox_coarse.je
    if (debug_log) print '("[INFO] WDR FNHF npe=",I0," direction=",I0,"   bbox_fine(",I0,"-",I0,",",I0,"-",I0,")")', this_pe, direction, bbox_fine.is, bbox_fine.ie, bbox_fine.js, bbox_fine.je


    if( bbox_coarse.ie .GE. bbox_coarse.is .AND. bbox_coarse.je .GE. bbox_coarse.js ) then
       if (debug_log) print '("[INFO] WDR BUFR Allocating large buffer. npe=",I0," i=",I0,"-",I0," j=",I0,"-",I0)', this_pe, bbox_coarse.is, bbox_coarse.ie, bbox_coarse.js, bbox_coarse.je
       allocate(buffer(bbox_coarse.is:bbox_coarse.ie, bbox_coarse.js:bbox_coarse.je), stat=my_stat, errmsg=my_errmsg)
       if (my_stat .ne. 0) print '("[ERROR] WDR NBFR error allocating buffer. npe=",I0,I0,A80)', this_pe, my_stat, my_errmsg

    else
       ! The buffer must have some storage allocated, whether it's a useful buffer or just a dummy.
       if (debug_log) print '("[INFO] WDR NBFR only allocating single entry buffer. npe=",I0," direction=",I0," i=",I0,"-",I0," j=",I0,"-",I0)', this_pe, direction, bbox_coarse.is, bbox_coarse.ie, bbox_coarse.js, bbox_coarse.je
       allocate(buffer(1,1))
    endif

    buffer = 0

  end subroutine alloc_2D_halo_buffer

  subroutine alloc_3D_halo_buffer(buffer, bbox_fine, bbox_coarse, nest_domain, direction, position, nz)
    real, dimension(:,:,:), allocatable, intent(out) :: buffer
    type(bbox), intent(out)                          :: bbox_fine, bbox_coarse 
    type(nest_domain_type), intent(inout)            :: nest_domain
    integer, intent(in)                              :: direction, position, nz


    integer                             :: my_stat
    character(256)                      :: my_errmsg
    integer                             :: this_pe

    this_pe = mpp_pe()

    call bbox_get_C2F_index(nest_domain, bbox_fine, bbox_coarse, direction,  position)
    if (debug_log) print '("[INFO] WDR FNHC npe=",I0," direction=",I0," bbox_coarse(",I0,"-",I0,",",I0,"-",I0,")")', this_pe, direction, bbox_coarse.is, bbox_coarse.ie, bbox_coarse.js, bbox_coarse.je
    if (debug_log) print '("[INFO] WDR FNHF npe=",I0," direction=",I0,"   bbox_fine(",I0,"-",I0,",",I0,"-",I0,")")', this_pe, direction, bbox_fine.is, bbox_fine.ie, bbox_fine.js, bbox_fine.je


    if( bbox_coarse.ie .GE. bbox_coarse.is .AND. bbox_coarse.je .GE. bbox_coarse.js ) then
       if (debug_log) print '("[INFO] WDR BUFR Allocating large buffer. npe=",I0," i=",I0,"-",I0," j=",I0,"-",I0," k=",I0,"-",I0)', this_pe, bbox_coarse.is, bbox_coarse.ie, bbox_coarse.js, bbox_coarse.je, 1, nz
       allocate(buffer(bbox_coarse.is:bbox_coarse.ie, bbox_coarse.js:bbox_coarse.je,1:nz), stat=my_stat, errmsg=my_errmsg)
       if (my_stat .ne. 0) print '("[ERROR] WDR NBFR error allocating buffer. npe=",I0,I0,A80)', this_pe, my_stat, my_errmsg

    else
       ! The buffer must have some storage allocated, whether it's a useful buffer or just a dummy.
       if (debug_log) print '("[INFO] WDR NBFR only allocating single entry buffer. npe=",I0," direction=",I0," i=",I0,"-",I0," j=",I0,"-",I0)', this_pe, direction, bbox_coarse.is, bbox_coarse.ie, bbox_coarse.js, bbox_coarse.je
       allocate(buffer(1,1,1))
    endif

    buffer = 0

  end subroutine alloc_3D_halo_buffer

  subroutine alloc_3D_halo_buffer_kindphys(buffer, bbox_fine, bbox_coarse, nest_domain, direction, position, nz)
    real(kind=kind_phys), dimension(:,:,:), allocatable, intent(out) :: buffer
    type(bbox), intent(out)                          :: bbox_fine, bbox_coarse 
    type(nest_domain_type), intent(inout)            :: nest_domain
    integer, intent(in)                              :: direction, position, nz


    integer                             :: my_stat
    character(256)                      :: my_errmsg
    integer                             :: this_pe

    this_pe = mpp_pe()

    call bbox_get_C2F_index(nest_domain, bbox_fine, bbox_coarse, direction,  position)
    if (debug_log) print '("[INFO] WDR FNHC npe=",I0," direction=",I0," bbox_coarse(",I0,"-",I0,",",I0,"-",I0,")")', this_pe, direction, bbox_coarse.is, bbox_coarse.ie, bbox_coarse.js, bbox_coarse.je
    if (debug_log) print '("[INFO] WDR FNHF npe=",I0," direction=",I0,"   bbox_fine(",I0,"-",I0,",",I0,"-",I0,")")', this_pe, direction, bbox_fine.is, bbox_fine.ie, bbox_fine.js, bbox_fine.je


    if( bbox_coarse.ie .GE. bbox_coarse.is .AND. bbox_coarse.je .GE. bbox_coarse.js ) then
       if (debug_log) print '("[INFO] WDR BUFR Allocating large buffer. npe=",I0," i=",I0,"-",I0," j=",I0,"-",I0," k=",I0,"-",I0)', this_pe, bbox_coarse.is, bbox_coarse.ie, bbox_coarse.js, bbox_coarse.je, 1, nz
       allocate(buffer(bbox_coarse.is:bbox_coarse.ie, bbox_coarse.js:bbox_coarse.je,1:nz), stat=my_stat, errmsg=my_errmsg)
       if (my_stat .ne. 0) print '("[ERROR] WDR NBFR error allocating buffer. npe=",I0,I0,A80)', this_pe, my_stat, my_errmsg

    else
       ! The buffer must have some storage allocated, whether it's a useful buffer or just a dummy.
       if (debug_log) print '("[INFO] WDR NBFR only allocating single entry buffer. npe=",I0," direction=",I0," i=",I0,"-",I0," j=",I0,"-",I0)', this_pe, direction, bbox_coarse.is, bbox_coarse.ie, bbox_coarse.js, bbox_coarse.je
       allocate(buffer(1,1,1))
    endif

    buffer = 0

  end subroutine alloc_3D_halo_buffer_kindphys



  subroutine alloc_4D_halo_buffer(buffer, bbox_fine, bbox_coarse, nest_domain, direction, position, nz, n4d)
    real, dimension(:,:,:,:), allocatable, intent(out) :: buffer
    type(bbox), intent(out)                          :: bbox_fine, bbox_coarse 
    type(nest_domain_type), intent(inout)            :: nest_domain
    integer, intent(in)                              :: direction, position, nz, n4d


    integer                             :: my_stat
    character(256)                      :: my_errmsg
    integer                             :: this_pe

    this_pe = mpp_pe()

    call bbox_get_C2F_index(nest_domain, bbox_fine, bbox_coarse, direction,  position)
    if (debug_log) print '("[INFO] WDR FNHC4 npe=",I0," direction=",I0," bbox_coarse(",I0,"-",I0,",",I0,"-",I0,")")', this_pe, direction, bbox_coarse.is, bbox_coarse.ie, bbox_coarse.js, bbox_coarse.je
    if (debug_log) print '("[INFO] WDR FNHF4 npe=",I0," direction=",I0,"   bbox_fine(",I0,"-",I0,",",I0,"-",I0,")")', this_pe, direction, bbox_fine.is, bbox_fine.ie, bbox_fine.js, bbox_fine.je


    if( bbox_coarse.ie .GE. bbox_coarse.is .AND. bbox_coarse.je .GE. bbox_coarse.js ) then
       if (debug_log) print '("[INFO] WDR BUFR4 Allocating large buffer. npe=",I0," i=",I0,"-",I0," j=",I0,"-",I0," k=",I0," n4d=",I0)', this_pe, bbox_coarse.is, bbox_coarse.ie, bbox_coarse.js, bbox_coarse.je, nz, n4d
       allocate(buffer(bbox_coarse.is:bbox_coarse.ie, bbox_coarse.js:bbox_coarse.je, 1:nz, 1:n4d), stat=my_stat, errmsg=my_errmsg)
       if (my_stat .ne. 0)  print '("[ERROR] WDR NBFR4 error allocating buffer. npe=",I0,I0,A80)', this_pe, my_stat, my_errmsg

    else
       ! The buffer must have some storage allocated, whether it's a useful buffer or just a dummy.
       if (debug_log) print '("[INFO] WDR NBFR4 only allocating single entry buffer. npe=",I0," direction=",I0," i=",I0,"-",I0," j=",I0,"-",I0)', this_pe, direction, bbox_coarse.is, bbox_coarse.ie, bbox_coarse.js, bbox_coarse.je
       allocate(buffer(1,1,1,1))
    endif

    buffer = 0

  end subroutine alloc_4D_halo_buffer


  subroutine alloc_4D_halo_buffer_kindphys(buffer, bbox_fine, bbox_coarse, nest_domain, direction, position, nz, n4d)
    real(kind=kind_phys), dimension(:,:,:,:), allocatable, intent(out) :: buffer
    type(bbox), intent(out)                          :: bbox_fine, bbox_coarse 
    type(nest_domain_type), intent(inout)            :: nest_domain
    integer, intent(in)                              :: direction, position, nz, n4d


    integer                             :: my_stat
    character(256)                      :: my_errmsg
    integer                             :: this_pe

    this_pe = mpp_pe()

    call bbox_get_C2F_index(nest_domain, bbox_fine, bbox_coarse, direction,  position)
    if (debug_log) print '("[INFO] WDR FNHC4 npe=",I0," direction=",I0," bbox_coarse(",I0,"-",I0,",",I0,"-",I0,")")', this_pe, direction, bbox_coarse.is, bbox_coarse.ie, bbox_coarse.js, bbox_coarse.je
    if (debug_log) print '("[INFO] WDR FNHF4 npe=",I0," direction=",I0,"   bbox_fine(",I0,"-",I0,",",I0,"-",I0,")")', this_pe, direction, bbox_fine.is, bbox_fine.ie, bbox_fine.js, bbox_fine.je


    if( bbox_coarse.ie .GE. bbox_coarse.is .AND. bbox_coarse.je .GE. bbox_coarse.js ) then
       if (debug_log) print '("[INFO] WDR BUFR4 Allocating large buffer. npe=",I0," i=",I0,"-",I0," j=",I0,"-",I0," k=",I0," n4d=",I0)', this_pe, bbox_coarse.is, bbox_coarse.ie, bbox_coarse.js, bbox_coarse.je, nz, n4d
       allocate(buffer(bbox_coarse.is:bbox_coarse.ie, bbox_coarse.js:bbox_coarse.je, 1:nz, 1:n4d), stat=my_stat, errmsg=my_errmsg)
       if (my_stat .ne. 0)  print '("[ERROR] WDR NBFR4 error allocating buffer. npe=",I0,I0,A80)', this_pe, my_stat, my_errmsg

    else
       ! The buffer must have some storage allocated, whether it's a useful buffer or just a dummy.
       if (debug_log) print '("[INFO] WDR NBFR4 only allocating single entry buffer. npe=",I0," direction=",I0," i=",I0,"-",I0," j=",I0,"-",I0)', this_pe, direction, bbox_coarse.is, bbox_coarse.ie, bbox_coarse.js, bbox_coarse.je
       allocate(buffer(1,1,1,1))
    endif

    buffer = 0

  end subroutine alloc_4D_halo_buffer_kindphys




  !==================================================================================================
  !
  !  Load static data from netCDF files
  !
  !==================================================================================================


  ! Load the full panel nest latlons from netCDF file
  ! character(*), parameter      :: nc_filename = '/scratch2/NAGAPE/aoml-hafs1/William.Ramstrom/static_grids/C384_grid.tile6.nc' 
  ! Read in the lat/lon in degrees, convert to radians

  subroutine load_nest_latlons_from_nc(nc_filename, nxp, nyp, refine, &
       fp_tile_geo, &
       fp_istart_fine, fp_iend_fine, fp_jstart_fine, fp_jend_fine)
    implicit none

    character(*), intent(in)              :: nc_filename
    integer, intent(in)                   :: nxp, nyp, refine
    type(grid_geometry), intent(out)      :: fp_tile_geo
    integer, intent(out)                  :: fp_istart_fine, fp_iend_fine, fp_jstart_fine, fp_jend_fine

    !========================================================================================
    !
    !  Determine which tile this PE is operating on
    !  Load the lat/lon data from netCDF file
    !  If fine nest, also determine the parent tile
    !  load the lat/lon data from that tile
    !  This code will only operate for nest motion within a single tile
    !
    !========================================================================================

    !  read lat/lon for this tile
    !  lat is y from grid file
    !  lon is x from grid file


    integer                      :: nx, ny

    integer                      :: nn
    integer                      :: super_nxp, super_nyp, mid_nx, mid_ny
    integer                      :: super_nx, super_ny
    type(grid_geometry)          :: temp_tile_geo
    ! Full panel nest data
    integer                      :: i, j, fi, fj
    integer                      :: this_pe


    real(kind=kind_phys) :: pi = 4d0 * atan(1.0d0)
    !real(kind=kind_phys) :: pi180, rad2deg, deg2rad
    real(kind=kind_phys) :: deg2rad

    !pi180 = pi / 180.0d0
    deg2rad = pi / 180.0d0
    !rad2deg = 1.0d0 / pi180


    this_pe = mpp_pe()

    if (debug_log) print '("[INFO] WDR NCREAD LLFE load_nest_latlons_from_nc fp interp_single_nest start, nread npe=",I0," nxp=",I0," nyp=",I0," refine=",I0)', this_pe, nxp, nyp, refine


    nx = nxp - 1
    ny = nyp - 1


    ! Global tiles don't have a halo in lat/lon data
    ! Nests have a halo in the lat/lon data
    !start = 1
    !nread = 1


    ! single fine nest
    ! full panel variables
    !fp_istart_fine = 12
    !fp_iend_fine = 269
    !fp_jstart_fine = 12
    !fp_jend_fine = 269
    !super_nx = 2*(fp_iend_fine - fp_istart_fine  + 1) + ( ehalo + whalo )
    !super_ny = 2*(fp_jend_fine - fp_jstart_fine  + 1) + ( nhalo + shalo )


    fp_istart_fine = 1
    fp_iend_fine = nx * refine
    fp_jstart_fine = 1
    fp_jend_fine = ny * refine
    super_nx = 2*(fp_iend_fine - fp_istart_fine  + 1)
    super_ny = 2*(fp_jend_fine - fp_jstart_fine  + 1)


    super_nxp = super_nx + 1
    super_nyp = super_ny + 1

    mid_nx = (fp_iend_fine - fp_istart_fine) 
    mid_ny = (fp_jend_fine - fp_jstart_fine) 



    if (debug_log) print '("[INFO] WDR LLFB load_nest_latlons_from_nc allocate fp fine temp_tile_geo%lats npe=",I0," dims: ",I4,":",I4,I4,":",I4,I4)', this_pe, 1, super_nxp, 1, super_nyp


    if (debug_log) print '("[INFO] WDR NCREAD LLFC load_nest_latlons_from_nc fp interp_single_nest. npe=",I0,I4,I4,I4,I4," ",A128)', this_pe, super_nxp, super_nyp, mid_nx,mid_ny, nc_filename


    call alloc_read_data(nc_filename, 'x', super_nxp, super_nyp, fp_tile_geo%lons)
    call alloc_read_data(nc_filename, 'y', super_nxp, super_nyp, fp_tile_geo%lats)
    call alloc_read_data(nc_filename, 'area', super_nx, super_ny, fp_tile_geo%area)

    !  double dx(nyp, nx) 
    !call alloc_read_data(nc_filename, 'dx', super_nx, super_nyp, fp_tile_geo%dx)
    ! double dy(ny, nxp)
    !call alloc_read_data(nc_filename, 'dy', super_nxp, super_ny, fp_tile_geo%dy)
    ! double area(ny, nx)
    !call alloc_read_data(nc_filename, 'area', super_nx, super_ny, fp_tile_geo%area)



    if (debug_log) print '("[INFO] WDR NCREAD LLFE load_nest_latlons_from_nc fp interp_single_nest start, nread npe=",I0)', this_pe


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! 
    !! Setup the lat/lons of the actual nest, read from the larger array
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !super_nxp = 2*(iend_fine - istart_fine  + 1) + 2 * ( ehalo + whalo ) + 1
    !super_nyp = 2*(jend_fine - jstart_fine  + 1) + 2 * ( nhalo + shalo ) + 1
    !mid_nx = (iend_fine - istart_fine) 
    !mid_ny = (jend_fine - jstart_fine) 

    ! end reading in nest

    if (debug_log) print '("[INFO] WDR CTR load_nest_latlons_from_nc center lat/lon. npe=",I0, " ", I0," ", I0)', this_pe, mid_nx, mid_ny
    if (debug_log) print '("[INFO] WDR CTR load_nest_latlons_from_nc center lat/lon. npe=",I0, " ", I0," ", I0)', this_pe, size(fp_tile_geo%lats, 1), size(fp_tile_geo%lats, 2)
    if (debug_log) print '("[INFO] WDR CTR load_nest_latlons_from_nc DEGS center lat/lon. npe=",I0,F8.2,F8.2," ",A128)', this_pe, fp_tile_geo%lats(mid_nx,mid_ny), fp_tile_geo%lons(mid_nx,mid_ny), nc_filename


    fp_tile_geo%lats = fp_tile_geo%lats * deg2rad
    fp_tile_geo%lons = fp_tile_geo%lons * deg2rad

    if (debug_log) print '("[INFO] WDR CTR load_nest_latlons_from_nc RADS center lat/lon. npe=",I0,F8.2,F8.2," ",A128)', this_pe, fp_tile_geo%lats(mid_nx,mid_ny), fp_tile_geo%lons(mid_nx,mid_ny), nc_filename


  end subroutine load_nest_latlons_from_nc



  subroutine load_nest_orog_from_nc(nc_filename, nxp, nyp, refine, orog_var_name, orog_grid, ls_mask_grid, land_frac_grid)
    implicit none

    character(*), intent(in)              :: nc_filename
    integer, intent(in)                   :: nxp, nyp, refine
    character(*), intent(in)              :: orog_var_name
    real, allocatable, intent(inout)      :: orog_grid(:,:)
    real, allocatable, intent(inout)      :: ls_mask_grid(:,:)
    real, allocatable, intent(inout)      :: land_frac_grid(:,:)

!    type(grid_geometry), intent(out)      :: fp_tile_geo


    integer                      :: fp_istart_fine, fp_iend_fine, fp_jstart_fine, fp_jend_fine


    integer                      :: nx, ny

    !integer                      :: nn
    integer                      :: fp_nx, fp_ny, mid_nx, mid_ny
    !integer                      :: super_nx, super_ny
    !type(grid_geometry)          :: temp_tile_geo
    ! Full panel nest data
    !integer                      :: i, j, fi, fj
    integer                      :: this_pe

    this_pe = mpp_pe()


    if (debug_log) print '("[INFO] WDR NCSTATIC load_nest_orog_from_nc start npe=",I0," nxp=",I0," nyp=",I0)', this_pe, nxp, nyp


    nx = nxp - 1
    ny = nyp - 1

    fp_istart_fine = 0
    fp_iend_fine = nx * refine
    fp_jstart_fine = 0
    fp_jend_fine = ny * refine

    fp_nx = fp_iend_fine - fp_istart_fine
    fp_ny = fp_jend_fine - fp_jstart_fine

    mid_nx = (fp_iend_fine - fp_istart_fine) / 2 
    mid_ny = (fp_jend_fine - fp_jstart_fine) / 2



    !if (debug_log) print '("[INFO] WDR LOFB load_nest_orog_from_nc temp_tile_geo%lats npe=",I0," dims: ",I4,":",I4,I4,":",I4,I4)', this_pe, 1, super_nxp, 1, super_nyp


    if (debug_log) print '("[INFO] WDR NCREAD LOFC load_nest_orog_from_nc npe=",I0,I4,I4,I4,I4," ",A128)', this_pe, fp_nx, fp_ny, mid_nx,mid_ny, nc_filename


    call alloc_read_data(nc_filename, orog_var_name, fp_nx, fp_ny, orog_grid)

    call alloc_read_data(nc_filename, orog_var_name, fp_nx, fp_ny, ls_mask_grid)
    call alloc_read_data(nc_filename, orog_var_name, fp_nx, fp_ny, land_frac_grid)


  end subroutine load_nest_orog_from_nc



  subroutine alloc_read_data_dyn(nc_filename, var_name, x_size, y_size, data_array)
    character(len=*), intent(in)           :: nc_filename, var_name
    integer, intent(in)                    :: x_size, y_size
    real, allocatable, intent(inout)       :: data_array(:,:)


    integer                      :: start(4), nread(4)    
    integer                      :: this_pe

    ! Allocate data_array to match the expected data size, then read in the data
    ! This subroutine consolidates the allocation and reading of data to ensure consistency of data sizing and simplify code 
    ! Could later extend this function to determine data size based on netCDF file metadata

    this_pe = mpp_pe()

    allocate(data_array(x_size, y_size))
    data_array = -9999.9

    if (debug_log) print '("[INFO] WDR alloc_read_data allocate  npe=",I0," ",A16," dims: ",I4,":",I4,I4,":",I4,I4)', this_pe, var_name, 1, x_size, 1, y_size

    start = 1
    nread = 1

    start(1) = 1
    start(2) = 1
    nread(1) = x_size
    nread(2) = y_size

    if (debug_log) print '("[INFO] WDR NCREAD NCRA alloc_read_data. npe=",I0," ",A96," ", A16)', this_pe, trim(nc_filename), var_name
    if (debug_log) print '("[INFO] WDR NCREAD NCRB alloc_read_data, nread npe=",I0, " ", A16,I4,I4,I4,I4)', this_pe, var_name, start(1), start(2), nread(1), nread(2)

    call read_data(nc_filename, var_name, data_array, start, nread, no_domain=.TRUE.)

    if (debug_log) print '("[INFO] WDR NCREAD NCRC alloc_read_data, nread npe=",I0, " ", A16,I4,I4,I4,I4)', this_pe, var_name, start(1), start(2), nread(1), nread(2)


  end subroutine alloc_read_data_dyn





  subroutine alloc_read_data_kind_phys(nc_filename, var_name, x_size, y_size, data_array)
    character(len=*), intent(in)           :: nc_filename, var_name
    integer, intent(in)                    :: x_size, y_size
    real(kind=kind_phys), allocatable, intent(inout)       :: data_array(:,:)


    integer                      :: start(4), nread(4)    
    integer                      :: this_pe

    ! Allocate data_array to match the expected data size, then read in the data
    ! This subroutine consolidates the allocation and reading of data to ensure consistency of data sizing and simplify code 
    ! Could later extend this function to determine data size based on netCDF file metadata

    this_pe = mpp_pe()

    allocate(data_array(x_size, y_size))
    data_array = -9999.9

    if (debug_log) print '("[INFO] WDR alloc_read_data allocate  npe=",I0," ",A16," dims: ",I4,":",I4,I4,":",I4,I4)', this_pe, var_name, 1, x_size, 1, y_size

    start = 1
    nread = 1

    start(1) = 1
    start(2) = 1
    nread(1) = x_size
    nread(2) = y_size

    if (debug_log) print '("[INFO] WDR NCREAD NCRA alloc_read_data. npe=",I0," ",A96," ", A16)', this_pe, trim(nc_filename), var_name
    if (debug_log) print '("[INFO] WDR NCREAD NCRB alloc_read_data, nread npe=",I0, " ", A16,I4,I4,I4,I4)', this_pe, var_name, start(1), start(2), nread(1), nread(2)

    call read_data(nc_filename, var_name, data_array, start, nread, no_domain=.TRUE.)

    if (debug_log) print '("[INFO] WDR NCREAD NCRC alloc_read_data, nread npe=",I0, " ", A16,I4,I4,I4,I4)', this_pe, var_name, start(1), start(2), nread(1), nread(2)


  end subroutine alloc_read_data_kind_phys



!  nest_geo and parent_geo can be centered or supergrids. 
!  Assumes and validates that nest_geo is smaller, and inside parent_geo
subroutine find_nest_alignment(nest_geo, parent_geo, nest_x, nest_y, parent_x, parent_y)
  implicit none
  type(grid_geometry), intent(in)     :: nest_geo, parent_geo
  integer, intent(out)                :: nest_x, nest_y, parent_x, parent_y

  type(bbox)  :: nest_bbox, parent_bbox
  integer     :: x,y
  logical     :: found = .false.

  real(kind=R_GRID) :: pi = 4 * atan(1.0d0)
  real                :: rad2deg

  rad2deg =  180.0 / pi

  parent_x = -999
  parent_y = -999
  nest_x = -999
  nest_y = -999


  if (debug_log) print '("[INFO] WDR start find_nest_alignment")'


  call fill_bbox(nest_bbox, nest_geo%lats)
  call show_bbox('nest', nest_bbox, nest_geo%lats, nest_geo%lons)
  call fill_bbox(parent_bbox, parent_geo%lats)
  call show_bbox('parent', parent_bbox, parent_geo%lats, parent_geo%lons)

  !parent_bbox%is = lbound(parent_geo%lats, 1)
  !parent_bbox%ie = ubound(parent_geo%lats, 1)
  !parent_bbox%js = lbound(parent_geo%lats, 2)
  !parent_bbox%je = ubound(parent_geo%lats, 2)

  do x = parent_bbox.is, parent_bbox.ie
     do y = parent_bbox.js, parent_bbox.je

        if (abs(parent_geo%lats(x,y) - nest_geo%lats(nest_bbox.is, nest_bbox.js)) .lt. 0.0001) then
           if (abs(parent_geo%lons(x,y) - nest_geo%lons(nest_bbox.is, nest_bbox.js)) .lt. 0.0001) then
              found = .true.

              parent_x = x
              parent_y = y
              nest_x = nest_bbox.is
              nest_y = nest_bbox.js

              if (debug_log) print '("[INFO] WDR find_nest_alignment parent(",I0,",",I0,") nest(",I0,",",I0,")")', x,y,nest_bbox.is, nest_bbox.js
              if (debug_log) print '("[INFO] WDR find_nest_alignment ",F10.5, F10.5)', parent_geo%lats(x,y)*rad2deg, parent_geo%lons(x,y)*rad2deg
           end if
        end if
     end do
  end do

  if (found) then
     if (debug_log) print '("[INFO] WDR find_nest_alignment MATCH FOUND",F10.5, F10.5)', nest_geo%lats(nest_bbox.is, nest_bbox.js)*rad2deg, nest_geo%lons(nest_bbox.is, nest_bbox.js)*rad2deg
  else
     if (debug_log) print '("[INFO] WDR find_nest_alignment NO MATCH FOUND",F10.5, F10.5)', nest_geo%lats(nest_bbox.is, nest_bbox.js)*rad2deg, nest_geo%lons(nest_bbox.is, nest_bbox.js)*rad2deg
  end if

end subroutine find_nest_alignment




  !==================================================================================================
  !
  !  NetCDF Function Section
  !
  !==================================================================================================



  subroutine output_grid_3d_to_nc(flag, istart, iend, jstart, jend, k, grid, file_str, var_name, time_step, dom, position)
    implicit none

    character(len=*), intent(in)  :: flag
    integer, intent(in)           :: istart, iend, jstart, jend, k
    !real, intent(in)                 :: grid(istart:iend, jstart:jend, k)
    !real, allocatable, intent(in)    :: grid(:, :, :)
    !real, allocatable, dimension(:,:,:),   intent(in)   :: grid
    real, dimension(:,:,:),   intent(in)   :: grid


    character(len=*), intent(in)  :: file_str, var_name
    integer, intent(in)           :: time_step
    type(domain2d), intent(in)    :: dom
    integer, intent(in)           :: position

    integer              :: this_pe
    character(len=256)   :: filename

    this_pe = mpp_pe()

    if (debug_log) print '("[INFO] WDR output_grid_3d_to_nc calling write_data. ",A8," npe=",I0, " i=",I0,"-",I0, " j=",I0,"-",I0," grid(",I0,",",I0,",",I0,")")', &
         flag, this_pe, istart, iend, jstart, jend, size(grid,1), size(grid,2), size(grid,3)


    write (filename, "(A,A1,I0.3,A)") trim(file_str), "_", time_step, ".nc"

    ! Resolves to:
    !subroutine write_data_3d_new(filename, fieldname, data, domain, no_domain, scalar_or_1d, &
    !                         position, tile_count, data_default)
    !character(len=*),         intent(in)         :: filename, fieldname
    !real, dimension(:,:,:),   intent(in)         :: data
    !type(domain2d), optional, intent(in), target :: domain
    !real,           optional, intent(in)         :: data_default
    !logical,        optional, intent(in)         :: no_domain
    !logical,        optional, intent(in)         :: scalar_or_1d
    !integer,        optional, intent(in)         :: position, tile_count

    call write_data(filename, var_name, grid, dom, position=position)

  end subroutine output_grid_3d_to_nc

  subroutine output_grid_2d_to_nc(flag, istart, iend, jstart, jend, k, grid, file_str, var_name, time_step, dom, position)
    implicit none

    character(len=*), intent(in)  :: flag
    integer, intent(in)           :: istart, iend, jstart, jend, k
    real, dimension(:,:),   intent(in)   :: grid


    character(len=*), intent(in)  :: file_str, var_name
    integer, intent(in)           :: time_step
    type(domain2d), intent(in)    :: dom
    integer, intent(in)           :: position

    integer              :: this_pe
    character(len=256) :: filename

    this_pe = mpp_pe()

    if (debug_log) print '("[INFO] WDR output_grid_2d_to_nc calling write_data. ",A8," npe=",I0, " i=",I0,"-",I0, " j=",I0,"-",I0," grid(",I0,")")', &
         flag, this_pe, istart, iend, jstart, jend, size(grid,1), size(grid,2)


    write (filename, "(A,A1,I0.3,A)") trim(file_str), "_", time_step, ".nc"

    call write_data(filename, var_name, grid, dom, position=position)

  end subroutine output_grid_2d_to_nc





  !==================================================================================================
  !
  !  Fill Section
  !
  !==================================================================================================


  subroutine fill_grid32_from_supergrid(in_grid, stagger_type, fp_super_tile_geo, ioffset, joffset, x_refine, y_refine)
    implicit none
    real, allocatable, intent(inout)    :: in_grid(:,:,:)
    integer, intent(in)                 :: stagger_type   ! CENTER, CORNER
    type(grid_geometry), intent(in)     :: fp_super_tile_geo
    integer, intent(in)                 :: ioffset, joffset, x_refine, y_refine


    integer :: nest_x, nest_y, parent_x, parent_y

    type(bbox)  :: tile_bbox, fp_tile_bbox
    integer   :: i, j, fp_i, fp_j

    ! tile_geo is cell-centered, at nest refinement
    ! fp_super_tile_geo is a supergrid, at nest refinement

    !call find_nest_alignment(tile_geo, fp_super_tile_geo, nest_x, nest_y, parent_x, parent_y)

    !  There are a few different offsets operating here:
    !  1. ioffset,joffset is how far the start of the (centered/corner?) grid is from the start of the parent grid
    !      i.e. the index of the parent center cell (not supergrid!) where the nest compute domain begins 
    !  2. nest_x, nest_y are the initial indices of this tile of the nest (the patch running on the PE)
    !  2. parent_x, parent_y are the initial indices of this tile of the parent supergrid (the patch running on the PE)
    !  3. parent_x = ((ioffset -1) * x_refine + nest_x) * 2
    !  


    call fill_bbox(tile_bbox, in_grid)
    call fill_bbox(fp_tile_bbox, fp_super_tile_geo%lats)

    ! Calculate new parent alignment -- supergrid at the refine ratio
    nest_x = tile_bbox%is
    nest_y = tile_bbox%js

    parent_x = ((ioffset - 1) * x_refine + nest_x) * 2
    parent_y = ((joffset - 1) * y_refine + nest_y) * 2


    do i = tile_bbox%is, tile_bbox%ie
       do j = tile_bbox%js, tile_bbox%je
          if (stagger_type == CENTER) then
             fp_i = (i - nest_x) * 2 + parent_x
             fp_j = (j - nest_y) * 2 + parent_y
          elseif (stagger_type == CORNER) then
             fp_i = (i - nest_x) * 2 + parent_x - 1
             fp_j = (j - nest_y) * 2 + parent_y - 1
          end if

          ! Make sure we don't run off the edge of the parent supergrid
          if (fp_i < fp_tile_bbox%is .or. fp_i > fp_tile_bbox%ie) then
             print '("[ERROR] WDR move_nest_geo invalid fp_i=",I0," is=",I0," ie=",I0)', fp_i, fp_tile_bbox%is, fp_tile_bbox%ie
             stop  ! TODO replace with a fatal error
          end if
          if (fp_j < fp_tile_bbox%js .or. fp_j > fp_tile_bbox%je) then
             print '("[ERROR] WDR move_nest_geo invalid fp_j=",I0," js=",I0," je=",I0)', fp_j, fp_tile_bbox%js, fp_tile_bbox%je
             stop  ! TODO replace with a fatal error
          end if

          in_grid(i,j,2) = fp_super_tile_geo%lats(fp_i, fp_j)
          in_grid(i,j,1) = fp_super_tile_geo%lons(fp_i, fp_j)
       end do
    end do

    ! Validate at the end
    !call find_nest_alignment(tile_geo, fp_super_tile_geo, nest_x, nest_y, parent_x, parent_y)


  end subroutine fill_grid32_from_supergrid


  subroutine fill_grid64_from_supergrid(in_grid, stagger_type, fp_super_tile_geo, ioffset, joffset, x_refine, y_refine)
    implicit none
    real(kind=R_GRID), allocatable, intent(inout)    :: in_grid(:,:,:)
    integer, intent(in)                 :: stagger_type   ! CENTER, CORNER
    type(grid_geometry), intent(in)     :: fp_super_tile_geo
    integer, intent(in)                 :: ioffset, joffset, x_refine, y_refine


    integer :: nest_x, nest_y, parent_x, parent_y

    type(bbox)  :: tile_bbox, fp_tile_bbox
    integer   :: i, j, fp_i, fp_j

    ! tile_geo is cell-centered, at nest refinement
    ! fp_super_tile_geo is a supergrid, at nest refinement

    !call find_nest_alignment(tile_geo, fp_super_tile_geo, nest_x, nest_y, parent_x, parent_y)

    !  There are a few different offsets operating here:
    !  1. ioffset,joffset is how far the start of the (centered/corner?) grid is from the start of the parent grid
    !      i.e. the index of the parent center cell (not supergrid!) where the nest compute domain begins 
    !  2. nest_x, nest_y are the initial indices of this tile of the nest (the patch running on the PE)
    !  2. parent_x, parent_y are the initial indices of this tile of the parent supergrid (the patch running on the PE)
    !  3. parent_x = ((ioffset -1) * x_refine + nest_x) * 2
    !  


    call fill_bbox(tile_bbox, in_grid)
    call fill_bbox(fp_tile_bbox, fp_super_tile_geo%lats)

    ! Calculate new parent alignment -- supergrid at the refine ratio
    nest_x = tile_bbox%is
    nest_y = tile_bbox%js

    parent_x = ((ioffset - 1) * x_refine + nest_x) * 2
    parent_y = ((joffset - 1) * y_refine + nest_y) * 2


    do i = tile_bbox%is, tile_bbox%ie
       do j = tile_bbox%js, tile_bbox%je
          if (stagger_type == CENTER) then
             fp_i = (i - nest_x) * 2 + parent_x
             fp_j = (j - nest_y) * 2 + parent_y
          elseif (stagger_type == CORNER) then
             fp_i = (i - nest_x) * 2 + parent_x - 1
             fp_j = (j - nest_y) * 2 + parent_y - 1
          end if

          ! Make sure we don't run off the edge of the parent supergrid
          if (fp_i < fp_tile_bbox%is .or. fp_i > fp_tile_bbox%ie) then
             print '("[ERROR] WDR move_nest_geo invalid fp_i=",I0," is=",I0," ie=",I0)', fp_i, fp_tile_bbox%is, fp_tile_bbox%ie
             stop  ! TODO replace with a fatal error
          end if
          if (fp_j < fp_tile_bbox%js .or. fp_j > fp_tile_bbox%je) then
             print '("[ERROR] WDR move_nest_geo invalid fp_j=",I0," js=",I0," je=",I0)', fp_j, fp_tile_bbox%js, fp_tile_bbox%je
             stop  ! TODO replace with a fatal error
          end if

          in_grid(i,j,2) = fp_super_tile_geo%lats(fp_i, fp_j)
          in_grid(i,j,1) = fp_super_tile_geo%lons(fp_i, fp_j)
       end do
    end do

    ! Validate at the end
    !call find_nest_alignment(tile_geo, fp_super_tile_geo, nest_x, nest_y, parent_x, parent_y)


  end subroutine fill_grid64_from_supergrid



  subroutine fill_grid64_4D_from_supergrid(in_grid, stagger_type, fp_super_tile_geo, ioffset, joffset, x_refine, y_refine)
    implicit none
    real(kind=R_GRID), allocatable, intent(inout)    :: in_grid(:,:,:,:)
    integer, intent(in)                 :: stagger_type   ! CENTER, CORNER
    type(grid_geometry), intent(in)     :: fp_super_tile_geo
    integer, intent(in)                 :: ioffset, joffset, x_refine, y_refine


    integer :: nest_x, nest_y, parent_x, parent_y

    type(bbox)  :: tile_bbox, fp_tile_bbox
    integer   :: i, j, fp_i, fp_j

    ! tile_geo is cell-centered, at nest refinement
    ! fp_super_tile_geo is a supergrid, at nest refinement

    !call find_nest_alignment(tile_geo, fp_super_tile_geo, nest_x, nest_y, parent_x, parent_y)

    !  There are a few different offsets operating here:
    !  1. ioffset,joffset is how far the start of the (centered/corner?) grid is from the start of the parent grid
    !      i.e. the index of the parent center cell (not supergrid!) where the nest compute domain begins 
    !  2. nest_x, nest_y are the initial indices of this tile of the nest (the patch running on the PE)
    !  2. parent_x, parent_y are the initial indices of this tile of the parent supergrid (the patch running on the PE)
    !  3. parent_x = ((ioffset -1) * x_refine + nest_x) * 2
    !  


    call fill_bbox(tile_bbox, in_grid)
    call fill_bbox(fp_tile_bbox, fp_super_tile_geo%lats)

    ! Calculate new parent alignment -- supergrid at the refine ratio
    nest_x = tile_bbox%is
    nest_y = tile_bbox%js

    parent_x = ((ioffset - 1) * x_refine + nest_x) * 2
    parent_y = ((joffset - 1) * y_refine + nest_y) * 2


    do i = tile_bbox%is, tile_bbox%ie
       do j = tile_bbox%js, tile_bbox%je
          if (stagger_type == CENTER) then
             fp_i = (i - nest_x) * 2 + parent_x
             fp_j = (j - nest_y) * 2 + parent_y
          elseif (stagger_type == CORNER) then
             fp_i = (i - nest_x) * 2 + parent_x - 1
             fp_j = (j - nest_y) * 2 + parent_y - 1
          end if

          ! Make sure we don't run off the edge of the parent supergrid
          if (fp_i < fp_tile_bbox%is .or. fp_i > fp_tile_bbox%ie) then
             print '("[ERROR] WDR move_nest_geo invalid fp_i=",I0," is=",I0," ie=",I0)', fp_i, fp_tile_bbox%is, fp_tile_bbox%ie
             stop  ! TODO replace with a fatal error
          end if
          if (fp_j < fp_tile_bbox%js .or. fp_j > fp_tile_bbox%je) then
             print '("[ERROR] WDR move_nest_geo invalid fp_j=",I0," js=",I0," je=",I0)', fp_j, fp_tile_bbox%js, fp_tile_bbox%je
             stop  ! TODO replace with a fatal error
          end if

          in_grid(i,j,2,1) = fp_super_tile_geo%lats(fp_i, fp_j)
          in_grid(i,j,1,1) = fp_super_tile_geo%lons(fp_i, fp_j)
       end do
    end do

    ! Validate at the end
    !call find_nest_alignment(tile_geo, fp_super_tile_geo, nest_x, nest_y, parent_x, parent_y)


  end subroutine fill_grid64_4D_from_supergrid


  !>@brief  This subroutine fills the nest halo data from the coarse grid data by downscaling.
  !>@details  Applicable to any interpolation type

  subroutine fill_nest_from_buffer2D(interp_type, x, buffer, bbox_fine, bbox_coarse, dir, x_refine, y_refine, wt, ind)
    implicit none

    integer, intent(in)                         :: interp_type
    real,    allocatable, intent(inout)         :: x(:,:)
    real,    allocatable, intent(in)            :: buffer(:,:)
    type(bbox), intent(in)                      :: bbox_fine, bbox_coarse
    integer, intent(in)                         :: dir, x_refine, y_refine
    real, allocatable, intent(in)               :: wt(:,:,:)    ! The final dimension is always 4                          
    integer, allocatable, intent(in)            :: ind(:,:,:) 

    integer   :: this_pe
    this_pe = mpp_pe()


    ! Output the interpolation type                                                                                          
    select case (interp_type)
    case (1)
       if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= cell centered")', this_pe, interp_type
       call fill_nest_from_buffer_cell_center("A", x, buffer, bbox_fine, bbox_coarse, dir, x_refine, y_refine, wt, ind)
       !     case (3)                                                                                                             
       !        if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= C grid staggered")', this_pe, interp_type     
    case (4)                                                                                                             
       if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= D grid staggered")', this_pe, interp_type
       call fill_nest_from_buffer_cell_center("D", x, buffer, bbox_fine, bbox_coarse, dir, x_refine, y_refine, wt, ind)     
    case (9)
       if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= nearest neighbor cell centered")', this_pe, interp_type
       !call fill_nest_from_buffer_nearest_neighbor(x, buffer, bbox_fine, bbox_coarse, dir, wt)
       call mpp_error(FATAL, '2D fill_nest_from_buffer_nearest_neighbor not yet implemented.')
    case default
       if (debug_log) print '("[ERROR] WDR FNB this_tile. npe=",I0," UNDEFINED interp_type=",I0)', this_pe, interp_type
       call mpp_error(FATAL, 'interp_single_nest got invalid value for interp_type from namelist.')
    end select

  end subroutine fill_nest_from_buffer2D

  subroutine fill_nest_from_buffer3D(interp_type, x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
    implicit none

    integer, intent(in)                         :: interp_type
    real,    allocatable, intent(inout)         :: x(:,:,:)
    real,    allocatable, intent(in)            :: buffer(:,:,:)
    type(bbox), intent(in)                      :: bbox_fine, bbox_coarse
    integer, intent(in)                         :: nz
    integer, intent(in)                         :: dir, x_refine, y_refine
    real, allocatable, intent(in)               :: wt(:,:,:)    ! The final dimension is always 4                          
    integer, allocatable, intent(in)               :: ind(:,:,:) 

    integer   :: this_pe
    this_pe = mpp_pe()


    ! Output the interpolation type                                                                                          
    select case (interp_type)
    case (1)
       if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= cell centered")', this_pe, interp_type
       call fill_nest_from_buffer_cell_center("A", x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
       !     case (3)                                                                                                             
       !        if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= C grid staggered")', this_pe, interp_type     
    case (4)                                                                                                             
       if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= D grid staggered")', this_pe, interp_type
       call fill_nest_from_buffer_cell_center("D", x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)     
    case (9)
       if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= nearest neighbor cell centered")', this_pe, interp_type
       call fill_nest_from_buffer_nearest_neighbor(x, buffer, bbox_fine, bbox_coarse, nz, dir, wt)
    case default
       if (debug_log) print '("[ERROR] WDR FNB this_tile. npe=",I0," UNDEFINED interp_type=",I0)', this_pe, interp_type
       call mpp_error(FATAL, 'interp_single_nest got invalid value for interp_type from namelist.')
    end select

  end subroutine fill_nest_from_buffer3D



  subroutine fill_nest_from_buffer3D_kindphys(interp_type, x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
    implicit none

    integer, intent(in)                         :: interp_type
    real(kind=kind_phys),    allocatable, intent(inout)         :: x(:,:,:)
    real(kind=kind_phys),    allocatable, intent(in)            :: buffer(:,:,:)
    type(bbox), intent(in)                      :: bbox_fine, bbox_coarse
    integer, intent(in)                         :: nz
    integer, intent(in)                         :: dir, x_refine, y_refine
    real, allocatable, intent(in)               :: wt(:,:,:)    ! The final dimension is always 4                          
    integer, allocatable, intent(in)               :: ind(:,:,:) 

    integer   :: this_pe
    this_pe = mpp_pe()


    ! Output the interpolation type                                                                                          
    select case (interp_type)
    case (1)
       if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= cell centered")', this_pe, interp_type
       call fill_nest_from_buffer_cell_center("A", x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
       !     case (3)                                                                                                             
       !        if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= C grid staggered")', this_pe, interp_type     
    case (4)                                                                                                             
       if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= D grid staggered")', this_pe, interp_type
       call fill_nest_from_buffer_cell_center("D", x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)     
    case (9)
       if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= nearest neighbor cell centered")', this_pe, interp_type
       call mpp_error(FATAL, 'nearest_neighbor is not yet implemented for fv_moving_nest_utils.F90::fill_nest_from_buffer_3D_kindphys')
       !call fill_nest_from_buffer_nearest_neighbor(x, buffer, bbox_fine, bbox_coarse, nz, dir, wt)
    case default
       if (debug_log) print '("[ERROR] WDR FNB this_tile. npe=",I0," UNDEFINED interp_type=",I0)', this_pe, interp_type
       call mpp_error(FATAL, 'interp_single_nest got invalid value for interp_type from namelist.')
    end select

  end subroutine fill_nest_from_buffer3D_kindphys

  !>@brief  This subroutine fills the nest halo data from the coarse grid data by downscaling.
  !>@details  Applicable to any interpolation type

  subroutine fill_nest_from_buffer4D(interp_type, x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
    implicit none

    integer, intent(in)                         :: interp_type
    real,    allocatable, intent(inout)         :: x(:,:,:,:)
    real,    allocatable, intent(in)            :: buffer(:,:,:,:)
    type(bbox), intent(in)                      :: bbox_fine, bbox_coarse
    integer, intent(in)                         :: nz
    integer, intent(in)                         :: dir, x_refine, y_refine
    real, allocatable, intent(in)               :: wt(:,:,:)    ! The final dimension is always 4                          
    integer, allocatable, intent(in)            :: ind(:,:,:) 

    integer   :: this_pe
    this_pe = mpp_pe()


    ! Output the interpolation type                                                                                          
    select case (interp_type)
    case (1)
       if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= cell centered")', this_pe, interp_type
       call fill_nest_from_buffer_cell_center("A", x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
       !     case (3)                                                                                                             
       !        if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= C grid staggered")', this_pe, interp_type     
    case (4)                                                                                                             
       if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= D grid staggered")', this_pe, interp_type
       call fill_nest_from_buffer_cell_center("D", x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)     
    case (9)
       if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= nearest neighbor cell centered")', this_pe, interp_type
       !call fill_nest_from_buffer_nearest_neighbor(x, buffer, bbox_fine, bbox_coarse, nz, dir, wt)
       call mpp_error(FATAL, '4D fill_nest_from_buffer_nearest_neighbor not yet implemented.')
    case default
       if (debug_log) print '("[ERROR] WDR FNB this_tile. npe=",I0," UNDEFINED interp_type=",I0)', this_pe, interp_type
       call mpp_error(FATAL, 'interp_single_nest got invalid value for interp_type from namelist.')
    end select

  end subroutine fill_nest_from_buffer4D


  subroutine fill_nest_from_buffer4D_kindphys(interp_type, x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
    implicit none

    integer, intent(in)                         :: interp_type
    real(kind=kind_phys), allocatable, intent(inout)         :: x(:,:,:,:)
    real(kind=kind_phys), allocatable, intent(in)            :: buffer(:,:,:,:)
    type(bbox), intent(in)                      :: bbox_fine, bbox_coarse
    integer, intent(in)                         :: nz
    integer, intent(in)                         :: dir, x_refine, y_refine
    real, allocatable, intent(in)               :: wt(:,:,:)    ! The final dimension is always 4                          
    integer, allocatable, intent(in)            :: ind(:,:,:) 

    integer   :: this_pe
    this_pe = mpp_pe()


    ! Output the interpolation type                                                                                          
    select case (interp_type)
    case (1)
       if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= cell centered")', this_pe, interp_type
       call fill_nest_from_buffer_cell_center("A", x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
       !     case (3)                                                                                                             
       !        if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= C grid staggered")', this_pe, interp_type     
    case (4)                                                                                                             
       if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= D grid staggered")', this_pe, interp_type
       call fill_nest_from_buffer_cell_center("D", x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)     
    case (9)
       if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= nearest neighbor cell centered")', this_pe, interp_type
       !call fill_nest_from_buffer_nearest_neighbor(x, buffer, bbox_fine, bbox_coarse, nz, dir, wt)
       call mpp_error(FATAL, '4D fill_nest_from_buffer_nearest_neighbor not yet implemented.')
    case default
       if (debug_log) print '("[ERROR] WDR FNB this_tile. npe=",I0," UNDEFINED interp_type=",I0)', this_pe, interp_type
       call mpp_error(FATAL, 'interp_single_nest got invalid value for interp_type from namelist.')
    end select

  end subroutine fill_nest_from_buffer4D_kindphys




  !>@brief  This subroutine fills the nest halo data from the coarse grid data by downscaling.  It can accommodate all grid staggers, using the stagger variable.  [The routine needs to be renamed since "_from_cell_center" has become incorrect.)
  !>@details  Applicable to any interpolation type

  subroutine fill_nest_from_buffer_cell_center2D(stagger, x, buffer, bbox_fine, bbox_coarse, dir, x_refine, y_refine, wt, ind)
    implicit none
    character ( len = 1 ), intent(in)             :: stagger
    real,    allocatable, intent(inout)           :: x(:,:)
    real,    allocatable, intent(in)              :: buffer(:,:)
    type(bbox), intent(in)                        :: bbox_fine, bbox_coarse
    integer, intent(in)                           :: dir, x_refine, y_refine
    real, allocatable, intent(in)                 :: wt(:,:,:)    ! The final dimension is always 4
    integer, allocatable, intent(in)              :: ind(:,:,:) 

    character(len=8)       :: dir_str
    integer                :: i, j, k, ic, jc
    integer                :: focus_i = 1
    integer                :: focus_j = 1
    integer                :: this_pe

    this_pe = mpp_pe()


    select case(dir)
    case (NORTH)
       dir_str = "NORTH"
    case (SOUTH)
       dir_str = "SOUTH"
    case (EAST)
       dir_str = "EAST"
    case (WEST)
       dir_str = "WEST"
    case default
       dir_str = "ERR DIR"
    end select

    if (debug_log) print '("[INFO] WDR FNFBCC start if (debug_log) print ",A1," ",A8,"  buffer. npe=",I0," buffer(",I0,"-",I0,",",I0,"-",I0,")")', stagger, dir_str, this_pe, lbound(buffer,1), ubound(buffer,1), lbound(buffer,2), ubound(buffer,2)

    if (debug_log) print '("[INFO] WDR FNFBCCX start print ",A1," ",A8,"  x. npe=",I0," x(",I0,"-",I0,",",I0,"-",I0,")")', stagger, dir_str, this_pe, lbound(x,1), ubound(x,1), lbound(x,2), ubound(x,2)

    if( bbox_coarse%ie .GE. bbox_coarse%is .AND. bbox_coarse%je .GE. bbox_coarse%js ) then

       if (debug_log) print '("[INFO] WDR BUFR print ",A8," large buffer. npe=",I0," buffer(is_c, js_c)=",F12.5," buffer(ie_c-1, je_c-1)=",F12.5)', dir_str, this_pe, buffer(bbox_coarse%is, bbox_coarse%js),  buffer(bbox_coarse%ie-1, bbox_coarse%je-1)

       if (debug_log) print '("[INFO] WDR ",A8," BOUNDS i npe=",I0," is_f=",I0," ie_f=",I0," is_c=",I0," ie_c=",I0)', dir_str, this_pe, bbox_fine%is, bbox_fine%ie, bbox_coarse%is, bbox_coarse%ie
       if (debug_log) print '("[INFO] WDR ",A8," BOUNDS j npe=",I0," js_f=",I0," je_f=",I0," js_c=",I0," je_c=",I0)', dir_str, this_pe, bbox_fine%js, bbox_fine%je, bbox_coarse%js, bbox_coarse%je

       do j=bbox_fine%js, bbox_fine%je
          do i=bbox_fine%is, bbox_fine%ie
             !if (stagger == "A") then
             !else if (stagger == "C") then
             !else if (stagger == "D") then
             !end if
             
             ic = ind(i,j,1)
             jc = ind(i,j,2)
             
             x(i,j) = &
                  wt(i,j,1)*buffer(ic,  jc  ) +  &
                  wt(i,j,2)*buffer(ic,  jc+1) +  &
                  wt(i,j,3)*buffer(ic+1,jc+1) +  &
                  wt(i,j,4)*buffer(ic+1,jc  )
             
             !call check_array(buffer, this_pe, "buffer"//dir_str, -300.0, 300.0)
             !call check_array(wt, this_pe, "wt"//dir_str, 0.0, 1.0)
             if (debug_log) print '("[INFO] WDR FILL WEIGHTS ",A8,"  npe=",I0," (",I0,",",I0,") ic,jc=(",I0,",",I0,"): wt:",F12.5,F12.5,F12.5,F12.5)', dir_str, this_pe, i, j, ic, jc, wt(i,j,1), wt(i,j,2), wt(i,j,3), wt(i,j,4)
          end do
       end do
    else
       if (debug_log) print '("[INFO] WDR NIL BUFR print ",A8,"  buffer. npe=",I0)', dir_str, this_pe       
    endif

    if (debug_log) print '("[INFO] WDR FILLNEST DONE print ",A8,"  buffer. npe=",I0)', dir_str, this_pe       
    
  end subroutine fill_nest_from_buffer_cell_center2D


  subroutine fill_nest_from_buffer_cell_center3D(stagger, x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
    implicit none
    character ( len = 1 ), intent(in)             :: stagger
    real,    allocatable, intent(inout)           :: x(:,:,:)
    real,    allocatable, intent(in)              :: buffer(:,:,:)
    !real, intent(inout)                          :: x(:,:,:)
    !real, intent(in)                             :: buffer(:,:,:)
    type(bbox), intent(in)                        :: bbox_fine, bbox_coarse
    integer, intent(in)                           :: nz
    integer, intent(in)                           :: dir, x_refine, y_refine
    !real, intent(in)                             :: wt(:,:,:)    ! The final dimension is always 4
    real, allocatable, intent(in)                 :: wt(:,:,:)    ! The final dimension is always 4
    integer, allocatable, intent(in)              :: ind(:,:,:) 



    character(len=8)       :: dir_str
    integer                :: i, j, k, ic, jc
    integer                :: focus_i = 1
    integer                :: focus_j = 1
    integer                :: this_pe

    this_pe = mpp_pe()


    select case(dir)
    case (NORTH)
       dir_str = "NORTH"
    case (SOUTH)
       dir_str = "SOUTH"
    case (EAST)
       dir_str = "EAST"
    case (WEST)
       dir_str = "WEST"
    case default
       dir_str = "ERR DIR"
    end select

    if (debug_log) print '("[INFO] WDR FNFBCC start if (debug_log) print ",A1," ",A8,"  buffer. npe=",I0," buffer(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', stagger, dir_str, this_pe, lbound(buffer,1), ubound(buffer,1), lbound(buffer,2), ubound(buffer,2), lbound(buffer,3), ubound(buffer,3)

    if (debug_log) print '("[INFO] WDR FNFBCCX start print ",A1," ",A8,"  x. npe=",I0," x(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', stagger, dir_str, this_pe, lbound(x,1), ubound(x,1), lbound(x,2), ubound(x,2), lbound(x,3), ubound(x,3)

    if( bbox_coarse%ie .GE. bbox_coarse%is .AND. bbox_coarse%je .GE. bbox_coarse%js ) then

       if (debug_log) print '("[INFO] WDR BUFR print ",A8," large buffer. npe=",I0," buffer(is_c, js_c, nz)=",F12.5," buffer(ie_c-1, je_c-1, nz)=",F12.5)', dir_str, this_pe, buffer(bbox_coarse%is, bbox_coarse%js, nz),  buffer(bbox_coarse%ie-1, bbox_coarse%je-1, nz)

       if (debug_log) print '("[INFO] WDR ",A8," BOUNDS i npe=",I0," is_f=",I0," ie_f=",I0," is_c=",I0," ie_c=",I0)', dir_str, this_pe, bbox_fine%is, bbox_fine%ie, bbox_coarse%is, bbox_coarse%ie
       if (debug_log) print '("[INFO] WDR ",A8," BOUNDS j npe=",I0," js_f=",I0," je_f=",I0," js_c=",I0," je_c=",I0)', dir_str, this_pe, bbox_fine%js, bbox_fine%je, bbox_coarse%js, bbox_coarse%je

       do k=1,nz
          do j=bbox_fine%js, bbox_fine%je
             do i=bbox_fine%is, bbox_fine%ie
                !if (stagger == "A") then
                !else if (stagger == "C") then
                !else if (stagger == "D") then
                !end if

                ic = ind(i,j,1)
                jc = ind(i,j,2)

                x(i,j,k) = &
                     wt(i,j,1)*buffer(ic,  jc,  k) +  &
                     wt(i,j,2)*buffer(ic,  jc+1,k) +  &
                     wt(i,j,3)*buffer(ic+1,jc+1,k) +  &
                     wt(i,j,4)*buffer(ic+1,jc,  k)

                !call check_array(buffer, this_pe, "buffer"//dir_str, -300.0, 300.0)
                !call check_array(wt, this_pe, "wt"//dir_str, 0.0, 1.0)
                if (debug_log) print '("[INFO] WDR FILL WEIGHTS ",A8,"  npe=",I0," (",I0,",",I0,",",I0,") ic,jc=(",I0,",",I0,"): wt:",F12.5,F12.5,F12.5,F12.5)', dir_str, this_pe, i, j, k, ic, jc, wt(i,j,1), wt(i,j,2), wt(i,j,3), wt(i,j,4)
                if (debug_log) print '("[INFO] WDR FILL WEIGHTS ",A8,"  npe=",I0," (",I0,",",I0,",",I0,") : buffer:",F12.5,F12.5,F12.5,F12.5)', dir_str, this_pe, i, j, k, buffer(ic,jc,k), buffer(ic,jc+1,k), buffer(ic+1,jc+1,k), buffer(ic+1,jc,k)
                if (debug_log) print '("[INFO] WDR after FILL nest from ",A8," buffer. npe=",I0," x(",I0,",",I0,",",I0,")=",F12.5)', dir_str, this_pe, i, j, k, x(i,j,k)
                if (debug_log) print '("[INFO] WDR FILLNEST from ",A8," buffer. npe=",I0," i,j=(",I0,",",I0,") ic,jc=(",I0,",",I0") x=",F12.5)', dir_str, this_pe, i, j, ic, jc, x(i,j,k) 
                !! Debugging printing
                !if ( ( i == focus_i ) .and. ( j == focus_j ) ) then
                !  if (debug_log) print '("[INFO] WDR FOCUS FILL WEIGHTS ",A8,"  npe=",I0," (",I0,",",I0,",",I0,") ic,jc=(",I0,",",I0,"): wt:",F12.5,F12.5,F12.5,F12.5)', dir_str, this_pe, i, j, k, ic, jc, wt(i,j,1), wt(i,j,2), wt(i,j,3), wt(i,j,4)
                !   if (debug_log) print '("[INFO] WDR FOCUS FILL WEIGHTS ",A8,"  npe=",I0," (",I0,",",I0,",",I0,") : buffer:",F12.5,F12.5,F12.5,F12.5)', dir_str, this_pe, i, j, k, buffer(ic,jc,k), buffer(ic,jc+1,k), buffer(ic+1,jc+1,k), buffer(ic+1,jc,k)
                !   if (debug_log) print '("[INFO] WDR FOCUS after FILL nest from ",A8," buffer. npe=",I0," x(",I0,",",I0,",",I0,")=",F12.5)', dir_str, this_pe, i, j, k, x(i,j,k)
                !   if (debug_log) print '("[INFO] WDR FOCUS FILLNEST from ",A8," buffer. npe=",I0," i,j=(",I0,",",I0,") ic,jc=(",I0,",",I0") x=",F12.5)', dir_str, this_pe, i, j, ic, jc, x(i,j,k) 
                !end if

             end do
          end do
       end do
    else
       if (debug_log) print '("[INFO] WDR NIL BUFR print ",A8,"  buffer. npe=",I0)', dir_str, this_pe       
       !if (debug_log) print '("[INFO WDR NIL BUFR ",A8," BOUNDS i npe=",I0,"is_f=",I0," ie_f=",I0,"is_c=",I0," ie_c=",I0)', dir_str, this_pe, is_f, ie_f, is_c, ie_c
       !if (debug_log) print '("[INFO WDR NIL BUFR ",A8," BOUNDS j npe=",I0,"js_f=",I0," je_f=",I0,"js_c=",I0," je_c=",I0)', dir_str, this_pe, js_f, je_f, js_c, je_c

    endif

    if (debug_log) print '("[INFO] WDR FILLNEST DONE print ",A8,"  buffer. npe=",I0)', dir_str, this_pe       

  end subroutine fill_nest_from_buffer_cell_center3D

  subroutine fill_nest_from_buffer_cell_center3D_kindphys(stagger, x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
    implicit none
    character ( len = 1 ), intent(in)             :: stagger
    real(kind=kind_phys), allocatable, intent(inout)           :: x(:,:,:)
    real(kind=kind_phys), allocatable, intent(in)              :: buffer(:,:,:)
    type(bbox), intent(in)                        :: bbox_fine, bbox_coarse
    integer, intent(in)                           :: nz
    integer, intent(in)                           :: dir, x_refine, y_refine
    real, allocatable, intent(in)                 :: wt(:,:,:)    ! The final dimension is always 4
    integer, allocatable, intent(in)              :: ind(:,:,:) 

    character(len=8)       :: dir_str
    integer                :: i, j, k, ic, jc
    integer                :: focus_i = 1
    integer                :: focus_j = 1
    integer                :: this_pe

    this_pe = mpp_pe()


    select case(dir)
    case (NORTH)
       dir_str = "NORTH"
    case (SOUTH)
       dir_str = "SOUTH"
    case (EAST)
       dir_str = "EAST"
    case (WEST)
       dir_str = "WEST"
    case default
       dir_str = "ERR DIR"
    end select

    if (debug_log) print '("[INFO] WDR FNFBCC start if (debug_log) print ",A1," ",A8,"  buffer. npe=",I0," buffer(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', stagger, dir_str, this_pe, lbound(buffer,1), ubound(buffer,1), lbound(buffer,2), ubound(buffer,2), lbound(buffer,3), ubound(buffer,3)

    if (debug_log) print '("[INFO] WDR FNFBCCX start print ",A1," ",A8,"  x. npe=",I0," x(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', stagger, dir_str, this_pe, lbound(x,1), ubound(x,1), lbound(x,2), ubound(x,2), lbound(x,3), ubound(x,3)

    if( bbox_coarse%ie .GE. bbox_coarse%is .AND. bbox_coarse%je .GE. bbox_coarse%js ) then

       if (debug_log) print '("[INFO] WDR BUFR print ",A8," large buffer. npe=",I0," buffer(is_c, js_c, nz)=",F12.5," buffer(ie_c-1, je_c-1, nz)=",F12.5)', dir_str, this_pe, buffer(bbox_coarse%is, bbox_coarse%js, nz),  buffer(bbox_coarse%ie-1, bbox_coarse%je-1, nz)

       if (debug_log) print '("[INFO] WDR ",A8," BOUNDS i npe=",I0," is_f=",I0," ie_f=",I0," is_c=",I0," ie_c=",I0)', dir_str, this_pe, bbox_fine%is, bbox_fine%ie, bbox_coarse%is, bbox_coarse%ie
       if (debug_log) print '("[INFO] WDR ",A8," BOUNDS j npe=",I0," js_f=",I0," je_f=",I0," js_c=",I0," je_c=",I0)', dir_str, this_pe, bbox_fine%js, bbox_fine%je, bbox_coarse%js, bbox_coarse%je

       do k=1,nz
          do j=bbox_fine%js, bbox_fine%je
             do i=bbox_fine%is, bbox_fine%ie
                !if (stagger == "A") then
                !else if (stagger == "C") then
                !else if (stagger == "D") then
                !end if

                ic = ind(i,j,1)
                jc = ind(i,j,2)

                x(i,j,k) = &
                     wt(i,j,1)*buffer(ic,  jc,  k) +  &
                     wt(i,j,2)*buffer(ic,  jc+1,k) +  &
                     wt(i,j,3)*buffer(ic+1,jc+1,k) +  &
                     wt(i,j,4)*buffer(ic+1,jc,  k)

                !call check_array(buffer, this_pe, "buffer"//dir_str, -300.0, 300.0)
                !call check_array(wt, this_pe, "wt"//dir_str, 0.0, 1.0)
                if (debug_log) print '("[INFO] WDR FILL WEIGHTS ",A8,"  npe=",I0," (",I0,",",I0,",",I0,") ic,jc=(",I0,",",I0,"): wt:",F12.5,F12.5,F12.5,F12.5)', dir_str, this_pe, i, j, k, ic, jc, wt(i,j,1), wt(i,j,2), wt(i,j,3), wt(i,j,4)
                if (debug_log) print '("[INFO] WDR FILL WEIGHTS ",A8,"  npe=",I0," (",I0,",",I0,",",I0,") : buffer:",F12.5,F12.5,F12.5,F12.5)', dir_str, this_pe, i, j, k, buffer(ic,jc,k), buffer(ic,jc+1,k), buffer(ic+1,jc+1,k), buffer(ic+1,jc,k)
                if (debug_log) print '("[INFO] WDR after FILL nest from ",A8," buffer. npe=",I0," x(",I0,",",I0,",",I0,")=",F12.5)', dir_str, this_pe, i, j, k, x(i,j,k)
                if (debug_log) print '("[INFO] WDR FILLNEST from ",A8," buffer. npe=",I0," i,j=(",I0,",",I0,") ic,jc=(",I0,",",I0") x=",F12.5)', dir_str, this_pe, i, j, ic, jc, x(i,j,k) 
                !! Debugging printing
                !if ( ( i == focus_i ) .and. ( j == focus_j ) ) then
                !  if (debug_log) print '("[INFO] WDR FOCUS FILL WEIGHTS ",A8,"  npe=",I0," (",I0,",",I0,",",I0,") ic,jc=(",I0,",",I0,"): wt:",F12.5,F12.5,F12.5,F12.5)', dir_str, this_pe, i, j, k, ic, jc, wt(i,j,1), wt(i,j,2), wt(i,j,3), wt(i,j,4)
                !   if (debug_log) print '("[INFO] WDR FOCUS FILL WEIGHTS ",A8,"  npe=",I0," (",I0,",",I0,",",I0,") : buffer:",F12.5,F12.5,F12.5,F12.5)', dir_str, this_pe, i, j, k, buffer(ic,jc,k), buffer(ic,jc+1,k), buffer(ic+1,jc+1,k), buffer(ic+1,jc,k)
                !   if (debug_log) print '("[INFO] WDR FOCUS after FILL nest from ",A8," buffer. npe=",I0," x(",I0,",",I0,",",I0,")=",F12.5)', dir_str, this_pe, i, j, k, x(i,j,k)
                !   if (debug_log) print '("[INFO] WDR FOCUS FILLNEST from ",A8," buffer. npe=",I0," i,j=(",I0,",",I0,") ic,jc=(",I0,",",I0") x=",F12.5)', dir_str, this_pe, i, j, ic, jc, x(i,j,k) 
                !end if

             end do
          end do
       end do
    else
       if (debug_log) print '("[INFO] WDR NIL BUFR print ",A8,"  buffer. npe=",I0)', dir_str, this_pe       
       !if (debug_log) print '("[INFO WDR NIL BUFR ",A8," BOUNDS i npe=",I0,"is_f=",I0," ie_f=",I0,"is_c=",I0," ie_c=",I0)', dir_str, this_pe, is_f, ie_f, is_c, ie_c
       !if (debug_log) print '("[INFO WDR NIL BUFR ",A8," BOUNDS j npe=",I0,"js_f=",I0," je_f=",I0,"js_c=",I0," je_c=",I0)', dir_str, this_pe, js_f, je_f, js_c, je_c

    endif

    if (debug_log) print '("[INFO] WDR FILLNEST DONE print ",A8,"  buffer. npe=",I0)', dir_str, this_pe       

  end subroutine fill_nest_from_buffer_cell_center3D_kindphys


  subroutine fill_nest_from_buffer_cell_center4D(stagger, x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
    implicit none
    character ( len = 1 ), intent(in)             :: stagger
    real,    allocatable, intent(inout)           :: x(:,:,:,:)
    real,    allocatable, intent(in)              :: buffer(:,:,:,:)
    type(bbox), intent(in)                        :: bbox_fine, bbox_coarse
    integer, intent(in)                           :: nz
    integer, intent(in)                           :: dir, x_refine, y_refine
    real, allocatable, intent(in)                 :: wt(:,:,:)    ! The final dimension is always 4
    integer, allocatable, intent(in)              :: ind(:,:,:) 

    character(len=8)       :: dir_str
    integer                :: i, j, k, v, ic, jc
    integer                :: focus_i = 1
    integer                :: focus_j = 1
    integer                :: this_pe

    this_pe = mpp_pe()


    select case(dir)
    case (NORTH)
       dir_str = "NORTH"
    case (SOUTH)
       dir_str = "SOUTH"
    case (EAST)
       dir_str = "EAST"
    case (WEST)
       dir_str = "WEST"
    case default
       dir_str = "ERR DIR"
    end select

    if (debug_log) print '("[INFO] WDR FNFBCC4D start print ",A1," ",A8,"  buffer. npe=",I0," buffer(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', stagger, dir_str, this_pe, lbound(buffer,1), ubound(buffer,1), lbound(buffer,2), ubound(buffer,2), lbound(buffer,3), ubound(buffer,3), lbound(buffer,4), ubound(buffer,4)

    if( bbox_coarse%ie .GE. bbox_coarse%is .AND. bbox_coarse%je .GE. bbox_coarse%js ) then

       if (debug_log) print '("[INFO] WDR BUFR print ",A8," large buffer. npe=",I0," buffer(is_c, js_c, nz, 1)=",F12.5," buffer(ie_c-1, je_c-1, nz, 1)=",F12.5)', dir_str, this_pe, buffer(bbox_coarse%is, bbox_coarse%js, nz, 1),  buffer(bbox_coarse%ie-1, bbox_coarse%je-1, nz, 1)

       if (debug_log) print '("[INFO] WDR ",A8," BOUNDS i npe=",I0,"is_f=",I0," ie_f=",I0," is_c=",I0," ie_c=",I0)', dir_str, this_pe, bbox_fine%is, bbox_fine%ie, bbox_coarse%is, bbox_coarse%ie
       if (debug_log) print '("[INFO] WDR ",A8," BOUNDS j npe=",I0,"js_f=",I0," je_f=",I0," js_c=",I0," je_c=",I0)', dir_str, this_pe, bbox_fine%js, bbox_fine%je, bbox_coarse%js, bbox_coarse%je

       do v=1,ubound(buffer,4)
          do k=1,nz
             do j=bbox_fine%js, bbox_fine%je
                do i=bbox_fine%is, bbox_fine%ie
                   ic = ind(i,j,1)
                   jc = ind(i,j,2)


                   !if (debug_log) print '("[INFO] WDR fill_nest from ",A8," buffer. npe=",I0," i,j=(",I0,",",I0,") ic,jc=(",I0,",",I0")")', dir_str, this_pe, i, j, ic, jc


                   !if (debug_log) print '("[INFO] WDR before FILL nest from ",A8," buffer. npe=",I0," x(",I0,",",I0,",",I0,")=",F12.5)', dir_str, this_pe, i, j, k, x(i,j,k,v)


                   !  Fill in with weighted interpolation
                   !                x(i,j,k) = &
                   !                     wt(i,j,1)*buffer(ic,  jc,  k) +  &
                   !                     wt(i,j,2)*buffer(ic,  jc+1,k) +  &
                   !                     wt(i,j,3)*buffer(ic+1,jc+1,k) +  &
                   !                     wt(i,j,4)*buffer(ic+1,jc,  k)

                   !        wt(iw,jw,1)=dist2*dist3      ! ic,   jc    weight
                   !        wt(iw,jw,2)=dist3*dist4      ! ic,   jc+2  weight
                   !        wt(iw,jw,3)=dist4*dist1      ! ic+2, jc+2  weight
                   !        wt(iw,jw,4)=dist1*dist2      ! ic+2, jc    weight



                   x(i,j,k,v) = &
                        wt(i,j,1)*buffer(ic,  jc,  k, v) +  &
                        wt(i,j,2)*buffer(ic,  jc+1,k, v) +  &
                        wt(i,j,3)*buffer(ic+1,jc+1,k, v) +  &
                        wt(i,j,4)*buffer(ic+1,jc,  k, v)

                   !call check_array(buffer, this_pe, "buffer"//dir_str, -300.0, 300.0)
                   !call check_array(wt, this_pe, "wt"//dir_str, 0.0, 1.0)


                   !if (debug_log) print '("[INFO] WDR FILL WEIGHTS ",A8,"  npe=",I0," (",I0,",",I0,",",I0,",",I0,") ic,jc=(",I0,",",I0,"): wt:",F12.5,F12.5,F12.5,F12.5)', dir_str, this_pe, i, j, k, v, ic, jc, wt(i,j,1), wt(i,j,2), wt(i,j,3), wt(i,j,4)
                   !if (debug_log) print '("[INFO] WDR FILL WEIGHTS ",A8,"  npe=",I0," (",I0,",",I0,",",I0,",",I0,") : buffer:",F12.5,F12.5,F12.5,F12.5)', dir_str, this_pe, i, j, k, v, buffer(ic,jc,k,v), buffer(ic,jc+1,k,v), buffer(ic+1,jc+1,k,v), buffer(ic+1,jc,k,v)


                   !if (debug_log) print '("[INFO] WDR after FILL nest from ",A8," buffer. npe=",I0," x(",I0,",",I0,",",I0,",",I0,")=",F12.5)', dir_str, this_pe, i, j, k, v,x(i,j,k,v)

                   !if (debug_log) print '("[INFO] WDR FILLNEST4D from ",A8," buffer. npe=",I0," i,j=(",I0,",",I0,") ic,jc=(",I0,",",I0") x=",F12.5)', dir_str, this_pe, i, j, ic, jc, x(i,j,k,v) 

                end do
             end do
          end do
       end do
    else
       if (debug_log) print '("[INFO] WDR NIL BUFR print ",A8,"  buffer. npe=",I0)', dir_str, this_pe       
    endif

    if (debug_log) print '("[INFO] WDR FILLNEST4D DONE print ",A8,"  buffer. npe=",I0)', dir_str, this_pe       

  end subroutine fill_nest_from_buffer_cell_center4D


  subroutine fill_nest_from_buffer_cell_center4D_kindphys(stagger, x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
    implicit none
    character ( len = 1 ), intent(in)             :: stagger
    real(kind=kind_phys), allocatable, intent(inout)           :: x(:,:,:,:)
    real(kind=kind_phys), allocatable, intent(in)              :: buffer(:,:,:,:)
    type(bbox), intent(in)                        :: bbox_fine, bbox_coarse
    integer, intent(in)                           :: nz
    integer, intent(in)                           :: dir, x_refine, y_refine
    real, allocatable, intent(in)                 :: wt(:,:,:)    ! The final dimension is always 4
    integer, allocatable, intent(in)              :: ind(:,:,:) 

    character(len=8)       :: dir_str
    integer                :: i, j, k, v, ic, jc
    integer                :: focus_i = 1
    integer                :: focus_j = 1
    integer                :: this_pe

    this_pe = mpp_pe()


    select case(dir)
    case (NORTH)
       dir_str = "NORTH"
    case (SOUTH)
       dir_str = "SOUTH"
    case (EAST)
       dir_str = "EAST"
    case (WEST)
       dir_str = "WEST"
    case default
       dir_str = "ERR DIR"
    end select

    if (debug_log) print '("[INFO] WDR FNFBCC4D start print ",A1," ",A8,"  buffer. npe=",I0," buffer(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', stagger, dir_str, this_pe, lbound(buffer,1), ubound(buffer,1), lbound(buffer,2), ubound(buffer,2), lbound(buffer,3), ubound(buffer,3), lbound(buffer,4), ubound(buffer,4)

    if( bbox_coarse%ie .GE. bbox_coarse%is .AND. bbox_coarse%je .GE. bbox_coarse%js ) then

       if (debug_log) print '("[INFO] WDR BUFR print ",A8," large buffer. npe=",I0," buffer(is_c, js_c, nz, 1)=",F12.5," buffer(ie_c-1, je_c-1, nz, 1)=",F12.5)', dir_str, this_pe, buffer(bbox_coarse%is, bbox_coarse%js, nz, 1),  buffer(bbox_coarse%ie-1, bbox_coarse%je-1, nz, 1)

       if (debug_log) print '("[INFO] WDR ",A8," BOUNDS i npe=",I0,"is_f=",I0," ie_f=",I0," is_c=",I0," ie_c=",I0)', dir_str, this_pe, bbox_fine%is, bbox_fine%ie, bbox_coarse%is, bbox_coarse%ie
       if (debug_log) print '("[INFO] WDR ",A8," BOUNDS j npe=",I0,"js_f=",I0," je_f=",I0," js_c=",I0," je_c=",I0)', dir_str, this_pe, bbox_fine%js, bbox_fine%je, bbox_coarse%js, bbox_coarse%je

       do v=1,ubound(buffer,4)
          do k=1,nz
             do j=bbox_fine%js, bbox_fine%je
                do i=bbox_fine%is, bbox_fine%ie
                   ic = ind(i,j,1)
                   jc = ind(i,j,2)


                   !if (debug_log) print '("[INFO] WDR fill_nest from ",A8," buffer. npe=",I0," i,j=(",I0,",",I0,") ic,jc=(",I0,",",I0")")', dir_str, this_pe, i, j, ic, jc


                   !if (debug_log) print '("[INFO] WDR before FILL nest from ",A8," buffer. npe=",I0," x(",I0,",",I0,",",I0,")=",F12.5)', dir_str, this_pe, i, j, k, x(i,j,k,v)


                   !  Fill in with weighted interpolation
                   !                x(i,j,k) = &
                   !                     wt(i,j,1)*buffer(ic,  jc,  k) +  &
                   !                     wt(i,j,2)*buffer(ic,  jc+1,k) +  &
                   !                     wt(i,j,3)*buffer(ic+1,jc+1,k) +  &
                   !                     wt(i,j,4)*buffer(ic+1,jc,  k)

                   !        wt(iw,jw,1)=dist2*dist3      ! ic,   jc    weight
                   !        wt(iw,jw,2)=dist3*dist4      ! ic,   jc+2  weight
                   !        wt(iw,jw,3)=dist4*dist1      ! ic+2, jc+2  weight
                   !        wt(iw,jw,4)=dist1*dist2      ! ic+2, jc    weight



                   x(i,j,k,v) = &
                        wt(i,j,1)*buffer(ic,  jc,  k, v) +  &
                        wt(i,j,2)*buffer(ic,  jc+1,k, v) +  &
                        wt(i,j,3)*buffer(ic+1,jc+1,k, v) +  &
                        wt(i,j,4)*buffer(ic+1,jc,  k, v)

                   !call check_array(buffer, this_pe, "buffer"//dir_str, -300.0, 300.0)
                   !call check_array(wt, this_pe, "wt"//dir_str, 0.0, 1.0)


                   !if (debug_log) print '("[INFO] WDR FILL WEIGHTS ",A8,"  npe=",I0," (",I0,",",I0,",",I0,",",I0,") ic,jc=(",I0,",",I0,"): wt:",F12.5,F12.5,F12.5,F12.5)', dir_str, this_pe, i, j, k, v, ic, jc, wt(i,j,1), wt(i,j,2), wt(i,j,3), wt(i,j,4)
                   !if (debug_log) print '("[INFO] WDR FILL WEIGHTS ",A8,"  npe=",I0," (",I0,",",I0,",",I0,",",I0,") : buffer:",F12.5,F12.5,F12.5,F12.5)', dir_str, this_pe, i, j, k, v, buffer(ic,jc,k,v), buffer(ic,jc+1,k,v), buffer(ic+1,jc+1,k,v), buffer(ic+1,jc,k,v)


                   !if (debug_log) print '("[INFO] WDR after FILL nest from ",A8," buffer. npe=",I0," x(",I0,",",I0,",",I0,",",I0,")=",F12.5)', dir_str, this_pe, i, j, k, v,x(i,j,k,v)

                   !if (debug_log) print '("[INFO] WDR FILLNEST4D from ",A8," buffer. npe=",I0," i,j=(",I0,",",I0,") ic,jc=(",I0,",",I0") x=",F12.5)', dir_str, this_pe, i, j, ic, jc, x(i,j,k,v) 

                end do
             end do
          end do
       end do
    else
       if (debug_log) print '("[INFO] WDR NIL BUFR print ",A8,"  buffer. npe=",I0)', dir_str, this_pe       
    endif

    if (debug_log) print '("[INFO] WDR FILLNEST4D DONE print ",A8,"  buffer. npe=",I0)', dir_str, this_pe       

  end subroutine fill_nest_from_buffer_cell_center4D_kindphys


  subroutine fill_nest_from_buffer_nearest_neighbor(x, buffer, bbox_fine, bbox_coarse, nz, dir, wt)
    implicit none

    real,    allocatable, intent(inout)         :: x(:,:,:)
    real,    allocatable, intent(in)            :: buffer(:,:,:)
    type(bbox), intent(in)                      :: bbox_fine, bbox_coarse
    integer, intent(in)                         :: dir
    real, allocatable, intent(in)               :: wt(:,:,:)    ! The final dimension is always 4
    integer, intent(in)                         :: nz



    character(len=8)       :: dir_str
    integer                :: i, j, k, ic, jc
    integer                :: nearest_idx

    integer   :: this_pe

    this_pe = mpp_pe()


    select case(dir)
    case (NORTH)
       dir_str = "NORTH"
    case (SOUTH)
       dir_str = "SOUTH"
    case (EAST)
       dir_str = "EAST"
    case (WEST)
       dir_str = "WEST"
    case default
       dir_str = "ERR DIR"
    end select

    if( bbox_coarse%ie .GE. bbox_coarse%is .AND. bbox_coarse%je .GE. bbox_coarse%js ) then
       if (debug_log) print '("[INFO] WDR BUFR print ",A8," large buffer. npe=",I0," buffer(is_c, js_c, nz)=",F12.5," buffer(ie_c-1, je_c-1, nz)=",F12.5)', dir_str, this_pe, buffer(bbox_coarse%is, bbox_coarse%js, nz),  buffer(bbox_coarse%ie-1, bbox_coarse%je-1, nz)

       if (debug_log) print '("[INFO WDR ",A8," BOUNDS i npe=",I0,"is_f=",I0," ie_f=",I0,"is_c=",I0," ie_c=",I0)', dir_str, this_pe, bbox_fine%is, bbox_fine%ie, bbox_coarse%is, bbox_coarse%ie
       if (debug_log) print '("[INFO WDR ",A8," BOUNDS j npe=",I0,"js_f=",I0," je_f=",I0,"js_c=",I0," je_c=",I0)', dir_str, this_pe, bbox_fine%js, bbox_fine%je, bbox_coarse%js, bbox_coarse%je


       do j=bbox_fine%js, bbox_fine%je
          do i=bbox_fine%is, bbox_fine%ie

             ! ic = (ie_c - is_c) / (ie_f - is_c)
             ic = bbox_coarse%is + 1
             jc = bbox_coarse%js + 1

             do k=1,nz

                if (debug_log) print '("[INFO] WDR before FILL nest from ",A8," buffer. npe=",I0," x(",I0,",",I0,",",I0,")=",F12.5)', dir_str, this_pe, i, j, k, x(i,j,k)

                ! Pick the maximum weight of the 4
                !   If two are tied for the max weight, use whichever one maxloc returns first
                !   TODO Might need a more deterministic algorithm here for reproducibility;  e.g. take the lowest index, etc. 
                nearest_idx = maxloc(wt(i, j, :), 1)
                if (debug_log) print '("[INFO] WDR Nearest Neighbor algorithm index ",I0," buffer. npe=",I0)', nearest_idx, this_pe


                !!  Fill in with weighted interpolation
                !x(i,j,k) = &
                !     wt(i,j,1)*buffer(ic,  jc,  k) +  &
                !     wt(i,j,2)*buffer(ic,  jc+1,k) +  &
                !     wt(i,j,3)*buffer(ic+1,jc+1,k) +  &
                !     wt(i,j,4)*buffer(ic+1,jc,  k)

                select case (nearest_idx)
                case (1)
                   x(i,j,k) = buffer(ic,  jc,  k)
                case (2)
                   x(i,j,k) = buffer(ic,  jc+1,k) 
                case (3)
                   x(i,j,k) = buffer(ic+1,jc+1,k) 
                case (4)
                   x(i,j,k) = buffer(ic+1,jc,  k)
                case default
                   ! Fill in with first value and warn
                   x(i,j,k) = buffer(ic,  jc,  k)
                   if (debug_log) print '("[WARN] WDR Nearest Neighbor algorithm mismatch index ",I0," buffer. npe=",I0," x(",I0,",",I0,",",I0,")=",F12.5)', nearest_idx, this_pe, i, j, k, x(i,j,k)
                end select

                if (debug_log) print '("[INFO] WDR FILL WEIGHTS ",A8,"  npe=",I0," (",I0,",",I0,",",I0,") : wt:",F12.5,F12.5,F12.5,F12.5)', dir_str, this_pe, i, j, k, wt(i,j,1), wt(i,j,2), wt(i,j,3), wt(i,j,4)


                if (debug_log) print '("[INFO] WDR after FILL nest from ",A8," buffer. npe=",I0," x(",I0,",",I0,",",I0,")=",F12.5)', dir_str, this_pe, i, j, k, x(i,j,k)
             end do
          end do
       end do
    else
       if (debug_log) print '("[INFO] WDR NIL BUFR print ",A8,"  buffer. npe=",I0)', dir_str, this_pe       
    endif

  end subroutine fill_nest_from_buffer_nearest_neighbor



  subroutine fill_weight_grid(atm_wt, new_wt)
    real, allocatable, intent(inout) :: atm_wt(:,:,:)
    real, allocatable, intent(in) :: new_wt(:,:,:)

    integer :: x,y,z,n
    integer :: this_pe

    this_pe = mpp_pe()

    do n=1,3
       if (lbound(atm_wt, n) .ne. lbound(new_wt, n)) then
          print '("[ERROR] WDR fill_weight_grid lbound mismatch fv_moving_nest.F90 npe=",I0," n=",I0, I0, I0)', this_pe, n, lbound(atm_wt, n), lbound(new_wt, n)
          stop
       end if
       if (ubound(atm_wt, n) .ne. ubound(new_wt, n)) then
          print '("[ERROR] WDR fill_weight_grid ubound mismatch fv_moving_nest.F90 npe=",I0," n=",I0, I0, I0)', this_pe, n, ubound(atm_wt, n), ubound(new_wt, n)
          stop
       end if
    end do

    if (debug_log) print '("[INFO] WDR running fill_weight_grid fv_moving_nest.F90 npe=",I0)', this_pe
    do x = lbound(atm_wt,1),ubound(atm_wt,1)
       do y = lbound(atm_wt,2),ubound(atm_wt,2)
          do z = 1,4
             atm_wt(x,y,z) = new_wt(x,y,z)
          end do
       end do
    end do

  end subroutine fill_weight_grid

#endif ! MOVING_NEST

end module fv_moving_nest_utils_mod
