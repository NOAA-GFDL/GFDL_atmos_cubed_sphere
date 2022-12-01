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
!! @brief   Provides subroutines to enable moving nest functionality in FV3 dynamic core.
!! @author W. Ramstrom, AOML/HRD   01/15/2021
!! @email William.Ramstrom@noaa.gov
! =======================================================================!

module fv_moving_nest_utils_mod

#ifdef MOVING_NEST
  use fms_mod,           only : mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE, clock_flag_default
  use mpp_mod,           only: FATAL, WARNING, MPP_DEBUG, NOTE, MPP_CLOCK_SYNC,MPP_CLOCK_DETAILED
  use mpp_mod,           only: mpp_pe, mpp_npes, mpp_root_pe, mpp_error, mpp_set_warn_level
  use mpp_mod,           only: mpp_declare_pelist, mpp_set_current_pelist, mpp_sync, mpp_sync_self
  use mpp_mod,           only: mpp_clock_begin, mpp_clock_end, mpp_clock_id
  use mpp_mod,           only: mpp_init, mpp_exit, mpp_chksum, stdout, stderr
  use mpp_mod,           only: input_nml_file
  use mpp_mod,           only: mpp_get_current_pelist, mpp_broadcast
  use mpp_domains_mod,   only: GLOBAL_DATA_DOMAIN, BITWISE_EXACT_SUM, BGRID_NE, CGRID_NE, DGRID_NE, AGRID
  use mpp_parameter_mod, only: AGRID_PARAM=>AGRID,CGRID_NE_PARAM=>CGRID_NE,SCALAR_PAIR
  use mpp_domains_mod,   only: FOLD_SOUTH_EDGE, FOLD_NORTH_EDGE, FOLD_WEST_EDGE, FOLD_EAST_EDGE
  use mpp_domains_mod,   only: MPP_DOMAIN_TIME, CYCLIC_GLOBAL_DOMAIN, NUPDATE,EUPDATE, XUPDATE, YUPDATE, SCALAR_PAIR
  use mpp_domains_mod,   only: domain1D, domain2D, DomainCommunicator2D, BITWISE_EFP_SUM
  use mpp_domains_mod,   only: mpp_get_compute_domain, mpp_get_data_domain, mpp_domains_set_stack_size
  use mpp_domains_mod,   only: mpp_global_field, mpp_global_sum, mpp_global_max, mpp_global_min
  use mpp_domains_mod,   only: mpp_domains_init, mpp_domains_exit, mpp_broadcast_domain
  use mpp_domains_mod,   only: mpp_update_domains, mpp_check_field, mpp_redistribute, mpp_get_memory_domain
  use mpp_domains_mod,   only: mpp_define_layout, mpp_define_domains, mpp_modify_domain
  use mpp_domains_mod,   only: mpp_define_io_domain
  use mpp_domains_mod,   only: mpp_get_neighbor_pe, mpp_define_mosaic, mpp_nullify_domain_list
  use mpp_domains_mod,   only: NORTH, NORTH_EAST, EAST, SOUTH_EAST, CORNER, CENTER
  use mpp_domains_mod,   only: SOUTH, SOUTH_WEST, WEST, NORTH_WEST, mpp_define_mosaic_pelist
  use mpp_domains_mod,   only: mpp_get_global_domain, ZERO, NINETY, MINUS_NINETY
  use mpp_domains_mod,   only: mpp_get_boundary, mpp_start_update_domains, mpp_complete_update_domains
  use mpp_domains_mod,   only: mpp_define_nest_domains, nest_domain_type
  use mpp_domains_mod,   only: mpp_get_C2F_index, mpp_update_nest_fine
  use mpp_domains_mod,   only: mpp_get_F2C_index, mpp_update_nest_coarse
  use mpp_domains_mod,   only: mpp_get_domain_shift, EDGEUPDATE, mpp_deallocate_domain
  use mpp_domains_mod,   only: mpp_group_update_type, mpp_create_group_update
  use mpp_domains_mod,   only: mpp_do_group_update, mpp_clear_group_update
  use mpp_domains_mod,   only: mpp_start_group_update, mpp_complete_group_update
  use mpp_domains_mod,   only: WUPDATE, SUPDATE, mpp_get_compute_domains, NONSYMEDGEUPDATE
  use mpp_domains_mod,   only: domainUG, mpp_define_unstruct_domain, mpp_get_UG_domain_tile_id
  use mpp_domains_mod,   only: mpp_get_UG_compute_domain, mpp_pass_SG_to_UG, mpp_pass_UG_to_SG
  use mpp_domains_mod,   only: mpp_get_ug_global_domain, mpp_global_field_ug
  use mpp_memutils_mod,  only: mpp_memuse_begin, mpp_memuse_end

#ifdef GFS_TYPES
  use GFS_typedefs,      only: kind_phys
#else
  use IPD_typedefs,      only: kind_phys => IPD_kind_phys
#endif

  use constants_mod,     only: grav

  use boundary_mod,      only: update_coarse_grid, update_coarse_grid_mpp
  use bounding_box_mod,  only: bbox, bbox_get_C2F_index, fill_bbox
  use fms2_io_mod,       only: read_data, write_data, open_file, close_file, register_axis, register_field
  use fms2_io_mod,       only: FmsNetcdfDomainFile_t, FmsNetcdfFile_t, is_dimension_registered

  use fv_arrays_mod,     only: R_GRID
  use fv_arrays_mod,     only: fv_grid_type, fv_nest_type, fv_atmos_type
  use fv_surf_map_mod,   only: FV3_zs_filter
  use fv_moving_nest_types_mod, only: grid_geometry
  use ifport,            only: getcwd

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


  interface alloc_read_data
#ifdef OVERLOAD_R8
    module procedure alloc_read_data_r4_2d
#endif
    module procedure alloc_read_data_r8_2d
  end interface alloc_read_data

  interface fill_nest_halos_from_parent
    module procedure fill_nest_halos_from_parent_r4_2d
    module procedure fill_nest_halos_from_parent_r4_3d
    module procedure fill_nest_halos_from_parent_r4_4d

    module procedure fill_nest_halos_from_parent_r8_2d
    module procedure fill_nest_halos_from_parent_r8_3d
    module procedure fill_nest_halos_from_parent_r8_4d
  end interface fill_nest_halos_from_parent

  interface alloc_halo_buffer
    module procedure alloc_halo_buffer_r4_2d
    module procedure alloc_halo_buffer_r4_3d
    module procedure alloc_halo_buffer_r4_4d

    module procedure alloc_halo_buffer_r8_2d
    module procedure alloc_halo_buffer_r8_3d
    module procedure alloc_halo_buffer_r8_4d
  end interface alloc_halo_buffer

  interface fill_nest_from_buffer
    module procedure fill_nest_from_buffer_r4_2d
    module procedure fill_nest_from_buffer_r4_3d
    module procedure fill_nest_from_buffer_r4_4d

    module procedure fill_nest_from_buffer_r8_2d
    module procedure fill_nest_from_buffer_r8_3d
    module procedure fill_nest_from_buffer_r8_4d
  end interface fill_nest_from_buffer

  interface fill_nest_from_buffer_cell_center
    module procedure fill_nest_from_buffer_cell_center_r4_2d
    module procedure fill_nest_from_buffer_cell_center_r4_3d
    module procedure fill_nest_from_buffer_cell_center_r4_4d

    module procedure fill_nest_from_buffer_cell_center_r8_2d
    module procedure fill_nest_from_buffer_cell_center_r8_3d
    module procedure fill_nest_from_buffer_cell_center_r8_4d
  end interface fill_nest_from_buffer_cell_center

  interface output_grid_to_nc
    module procedure output_grid_to_nc_2d
    module procedure output_grid_to_nc_3d
  end interface output_grid_to_nc

  interface fill_grid_from_supergrid
    module procedure fill_grid_from_supergrid_r4_3d
    module procedure fill_grid_from_supergrid_r8_3d
    module procedure fill_grid_from_supergrid_r8_4d
  end interface fill_grid_from_supergrid


contains

  ! GEMPAK 5-point smoother
  !SM5S  Smooth scalar grid using a 5-point smoother
  !      SM5S ( S ) = .5 * S (i,j) + .125 * ( S (i+1,j) + S (i,j+1) +
  !                                           S (i-1,j) + S (i,j-1) )
  ! GEMPAK 9-point smoother
  !SM9S  Smooth scalar grid using a 9-point smoother
  !      SM5S ( S ) = .25 * S (i,j) + .125 * ( S (i+1,j) + S (i,j+1) +
  !                                            S (i-1,j) + S (i,j-1) )
  !                                 + .0625 * ( S (i+1,j+1) +
  !                                             S (i+1,j-1) +
  !                                             S (i-1,j+1) +
  !                                             S (i-1,j-1) )


  subroutine smooth_5_point(data_var, i, j, val)
    real, allocatable, intent(in)               :: data_var(:,:)
    integer                                     :: i,j
    real, intent(out)                           :: val

    ! Stay in bounds of the array
    if ( (i-1) .ge. lbound(data_var,1) .and. i .le. ubound(data_var,1) .and. (j-1) .ge. lbound(data_var,2) .and. j .le. ubound(data_var,2) ) then
      val = .5 * data_var(i,j) + .125 * ( data_var(i+1,j) + data_var(i,j+1) + data_var(i-1,j) + data_var(i,j-1) )
    else
      ! Don't smooth if at the edge.  Could do partial smoothing here also, but don't expect moving nest to reach the edge.
      val = data_var(i,j)
    endif

  end subroutine smooth_5_point


  subroutine smooth_9_point(data_var, i, j, val)
    real, allocatable, intent(in)               :: data_var(:,:)
    integer                                     :: i,j
    real, intent(out)                           :: val

    ! Stay in bounds of the array
    if ( (i-1) .ge. lbound(data_var,1) .and. i .le. ubound(data_var,1) .and. (j-1) .ge. lbound(data_var,2) .and. j .le. ubound(data_var,2) ) then
      val = .25 * data_var(i,j) + .125 * ( data_var(i+1,j) + data_var(i,j+1) + data_var(i-1,j) + data_var(i,j-1) ) &
          + .0625 * ( data_var(i+1,j+1) + data_var(i+1,j-1) + data_var(i-1,j+1) + data_var(i-1,j-1) )
    else
      ! Don't smooth if at the edge.  Could do partial smoothing here also, but don't expect moving nest to reach the edge.
      val = data_var(i,j)
    endif

  end subroutine smooth_9_point

  ! blend_size is 5 for static nests.  We may increase it for moving nests.
  !  This is only called for fine PEs.
  !  Blends a few points into the nest.  Calls zs filtering if enabled in namelist.
  subroutine set_blended_terrain(Atm, parent_orog_grid, nest_orog_grid, refine, halo_size, blend_size, a_step)
    type(fv_atmos_type), intent(inout), target :: Atm
    real, allocatable, intent(in)              :: parent_orog_grid(:,:)   ! Coarse grid orography
    real, allocatable, intent(in)              :: nest_orog_grid(:,:)     ! orography for the full panel of the parent, at high-resolution
    integer, intent(in)                        :: refine, halo_size, blend_size, a_step

    integer            :: i, j, ic, jc
    integer            :: ioffset, joffset
    integer            :: npx, npy, isd, ied, jsd, jed
    real               :: smoothed_orog, hires_orog, blend_wt, blend_orog

    real, pointer, dimension(:,:,:) :: wt
    integer, pointer, dimension(:,:,:) :: ind
    integer :: this_pe

    this_pe = mpp_pe()

    npx   = Atm%npx
    npy   = Atm%npy

    isd = Atm%bd%isc - halo_size
    ied = Atm%bd%iec + halo_size
    jsd = Atm%bd%jsc - halo_size
    jed = Atm%bd%jec + halo_size

    ioffset = Atm%neststruct%ioffset
    joffset = Atm%neststruct%joffset

    wt => Atm%neststruct%wt_h
    ind => Atm%neststruct%ind_h

    do j=jsd, jed
      do i=isd, ied
        ic = ind(i,j,1)
        jc = ind(i,j,2)

        smoothed_orog = &
            wt(i,j,1)*parent_orog_grid(ic,  jc  ) +  &
            wt(i,j,2)*parent_orog_grid(ic,  jc+1) +  &
            wt(i,j,3)*parent_orog_grid(ic+1,jc+1) +  &
            wt(i,j,4)*parent_orog_grid(ic+1,jc  )

        hires_orog = nest_orog_grid((ioffset-1)*refine+i, (joffset-1)*refine+j)

        ! From tools/external_ic.F90
        if (blend_size .eq. 10) then
          blend_wt = max(0.,min(1.,real(10 - min(i,j,npx-i,npy-j,10))/10. ))
        else
          blend_wt = max(0.,min(1.,real(5 - min(i,j,npx-i,npy-j,5))/5. ))
        end if

        !blend_wt = max(0.,min(1.,real(blend_size - min(i,j,npx-i,npy-j,blend_size))/real(blend_size) ))
        blend_orog = (1.-blend_wt)*hires_orog + blend_wt*smoothed_orog

        Atm%phis(i,j) = blend_orog * grav

      enddo
    enddo


    ! From tools/fv_surf_map.F90::surfdrv()
    if ( Atm%flagstruct%full_zs_filter ) then
      !if(is_master()) then
      !  write(*,*) 'Applying terrain filters. zero_ocean is', zero_ocean
      !endif
      !call FV3_zs_filter (bd, isd, ied, jsd, jed, npx, npy, npx_global,  &
      !    stretch_fac, bounded_domain, domain, area, dxa, dya, dx, dy, dxc, dyc, grid,  &
      !    agrid, sin_sg, phis, oro_g)

      call FV3_zs_filter (Atm%bd, isd, ied, jsd, jed, Atm%npx, Atm%npy, Atm%neststruct%npx_global,  &
          Atm%flagstruct%stretch_fac, Atm%gridstruct%bounded_domain, Atm%domain, &
          Atm%gridstruct%area_64, Atm%gridstruct%dxa, Atm%gridstruct%dya, &
          Atm%gridstruct%dx, Atm%gridstruct%dy, &
          Atm%gridstruct%dxc, Atm%gridstruct%dyc, &
          Atm%gridstruct%grid_64,  &
          Atm%gridstruct%agrid_64, Atm%gridstruct%sin_sg, Atm%phis, parent_orog_grid)

      call mpp_update_domains(Atm%phis, Atm%domain)
    endif          ! end terrain filter

  end subroutine set_blended_terrain

  subroutine set_smooth_nest_terrain(Atm, fp_orog, refine, num_points, halo_size, blend_size)
    type(fv_atmos_type), intent(inout) :: Atm
    real, allocatable, intent(in)      :: fp_orog(:,:)   ! orography for the full panel of the parent, at high-resolution
    integer, intent(in)                :: refine, num_points, halo_size, blend_size

    integer            :: i,j
    integer            :: ioffset, joffset
    integer            :: npx, npy, isd, ied, jsd, jed
    integer            :: smooth_i_lo, smooth_i_hi, smooth_j_lo, smooth_j_hi
    real               :: smoothed_orog
    character(len=16)  :: errstring

    npx   = Atm%npx
    npy   = Atm%npy

    isd = Atm%bd%isc - halo_size
    ied = Atm%bd%iec + halo_size
    jsd = Atm%bd%jsc - halo_size
    jed = Atm%bd%jec + halo_size

    ioffset = Atm%neststruct%ioffset
    joffset = Atm%neststruct%joffset

    smooth_i_lo = 1 + blend_size
    smooth_i_hi = npx - blend_size - halo_size

    smooth_j_lo = 1 + blend_size
    smooth_j_hi = npy - blend_size - halo_size

    !Atm(n)%phis(isd:ied, jsd:jed) = mn_static%orog_grid((ioffset-1)*x_refine+isd:(ioffset-1)*x_refine+ied, (joffset-1)*y_refine+jsd:(joffset-1)*y_refine+jed) * grav

    select case(num_points)
    case (5)

      do j=jsd, jed
        do i=isd, ied
          if (i .lt. smooth_i_lo .or. i .gt. smooth_i_hi .or. j .lt. smooth_j_lo .or. j .gt. smooth_j_hi) then
            call smooth_5_point(fp_orog, (ioffset-1)*refine + i, (joffset-1)*refine + j, smoothed_orog)
            Atm%phis(i,j) = smoothed_orog * grav
          else
            Atm%phis(i,j) = fp_orog((ioffset-1)*refine + i, (joffset-1)*refine + j) * grav
          endif
        enddo
      enddo

    case (9)

      do j=jsd, jed
        do i=isd, ied
          if (i .lt. smooth_i_lo .or. i .gt. smooth_i_hi .or. j .lt. smooth_j_lo .or. j .gt. smooth_j_hi) then
            call smooth_9_point(fp_orog, (ioffset-1)*refine + i, (joffset-1)*refine + j, smoothed_orog)
            Atm%phis(i,j) = smoothed_orog * grav
          else
            Atm%phis(i,j) = fp_orog((ioffset-1)*refine + i, (joffset-1)*refine + j) * grav
          endif
        enddo
      enddo

    case default
      write (errstring, "(I0)") num_points
      call mpp_error(FATAL,'Invalid terrain_smoother in set_smooth_nest_terrain '//errstring)
    end select

  end subroutine set_smooth_nest_terrain

  !==================================================================================================
  !
  !  Fill Nest Halos from Parent
  !
  !==================================================================================================

  subroutine fill_nest_halos_from_parent_r4_2d(var_name, data_var, interp_type, wt, ind, x_refine, y_refine, is_fine_pe, nest_domain, position)
    character(len=*), intent(in)                :: var_name
    real*4, allocatable, intent(inout)          :: data_var(:,:)
    integer, intent(in)                         :: interp_type
    real, allocatable, intent(in)               :: wt(:,:,:)
    integer, allocatable, intent(in)            :: ind(:,:,:)
    integer, intent(in)                         :: x_refine, y_refine
    logical, intent(in)                         :: is_fine_pe
    type(nest_domain_type), intent(inout)       :: nest_domain
    integer, intent(in)                         :: position

    real*4, dimension(:,:), allocatable :: nbuffer, sbuffer, ebuffer, wbuffer
    type(bbox)                          :: north_fine, north_coarse
    type(bbox)                          :: south_fine, south_coarse
    type(bbox)                          :: east_fine, east_coarse
    type(bbox)                          :: west_fine, west_coarse
    integer                             :: this_pe
    integer                             :: nest_level = 1  ! TODO allow to vary

    this_pe = mpp_pe()

    !!===========================================================
    !!
    !! Fill halo buffers
    !!
    !!===========================================================

    call alloc_halo_buffer(nbuffer, north_fine, north_coarse, nest_domain, NORTH,  position)
    call alloc_halo_buffer(sbuffer, south_fine, south_coarse, nest_domain, SOUTH,  position)
    call alloc_halo_buffer(ebuffer, east_fine,  east_coarse,  nest_domain, EAST,   position)
    call alloc_halo_buffer(wbuffer, west_fine,  west_coarse,  nest_domain, WEST,   position)

    ! Passes data from coarse grid to fine grid's halo
    call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

    if (is_fine_pe) then

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

  end subroutine fill_nest_halos_from_parent_r4_2d


  subroutine fill_nest_halos_from_parent_r8_2d(var_name, data_var, interp_type, wt, ind, x_refine, y_refine, is_fine_pe, nest_domain, position)
    character(len=*), intent(in)                :: var_name
    real*8, allocatable, intent(inout)          :: data_var(:,:)
    integer, intent(in)                         :: interp_type
    real, allocatable, intent(in)               :: wt(:,:,:)   ! TODO should this also be real*8?
    integer, allocatable, intent(in)            :: ind(:,:,:)
    integer, intent(in)                         :: x_refine, y_refine
    logical, intent(in)                         :: is_fine_pe
    type(nest_domain_type), intent(inout)       :: nest_domain
    integer, intent(in)                         :: position


    real*8, dimension(:,:), allocatable :: nbuffer, sbuffer, ebuffer, wbuffer
    type(bbox)                          :: north_fine, north_coarse
    type(bbox)                          :: south_fine, south_coarse
    type(bbox)                          :: east_fine, east_coarse
    type(bbox)                          :: west_fine, west_coarse
    integer                             :: this_pe
    integer                             :: nest_level = 1  ! TODO allow to vary

    this_pe = mpp_pe()

    !!===========================================================
    !!
    !! Fill halo buffers
    !!
    !!===========================================================

    call alloc_halo_buffer(nbuffer, north_fine, north_coarse, nest_domain, NORTH,  position)
    call alloc_halo_buffer(sbuffer, south_fine, south_coarse, nest_domain, SOUTH,  position)
    call alloc_halo_buffer(ebuffer, east_fine,  east_coarse,  nest_domain, EAST,   position)
    call alloc_halo_buffer(wbuffer, west_fine,  west_coarse,  nest_domain, WEST,   position)

    ! Passes data from coarse grid to fine grid's halo
    call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

    if (is_fine_pe) then

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

  end subroutine fill_nest_halos_from_parent_r8_2d


  subroutine fill_nest_halos_from_parent_masked(var_name, data_var, interp_type, wt, ind, x_refine, y_refine, is_fine_pe, nest_domain, position, mask_var, mask_val, default_val)
    character(len=*), intent(in)                :: var_name
    real*8, allocatable, intent(inout)          :: data_var(:,:)
    integer, intent(in)                         :: interp_type
    real, allocatable, intent(in)               :: wt(:,:,:)   ! TODO should this also be real*8?
    integer, allocatable, intent(in)            :: ind(:,:,:)
    integer, intent(in)                         :: x_refine, y_refine
    logical, intent(in)                         :: is_fine_pe
    type(nest_domain_type), intent(inout)       :: nest_domain
    integer, intent(in)                         :: position
    real*4, allocatable, intent(in)             :: mask_var(:,:)
    integer, intent(in)                         :: mask_val
    real*8, intent(in)                          :: default_val

    real*8, dimension(:,:), allocatable :: nbuffer, sbuffer, ebuffer, wbuffer
    type(bbox)                          :: north_fine, north_coarse
    type(bbox)                          :: south_fine, south_coarse
    type(bbox)                          :: east_fine, east_coarse
    type(bbox)                          :: west_fine, west_coarse
    integer                             :: this_pe
    integer                             :: nest_level = 1  ! TODO allow to vary

    this_pe = mpp_pe()

    !!===========================================================
    !!
    !! Fill halo buffers
    !!
    !!===========================================================

    call alloc_halo_buffer(nbuffer, north_fine, north_coarse, nest_domain, NORTH,  position)
    call alloc_halo_buffer(sbuffer, south_fine, south_coarse, nest_domain, SOUTH,  position)
    call alloc_halo_buffer(ebuffer, east_fine,  east_coarse,  nest_domain, EAST,   position)
    call alloc_halo_buffer(wbuffer, west_fine,  west_coarse,  nest_domain, WEST,   position)

    ! Passes data from coarse grid to fine grid's halo
    call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

    if (is_fine_pe) then

      !!===========================================================
      !!
      !! Apply halo data
      !!
      !!===========================================================

      call fill_nest_from_buffer_masked(interp_type, data_var, nbuffer, north_fine, north_coarse, NORTH, x_refine, y_refine, wt, ind, mask_var, mask_val, default_val)
      call fill_nest_from_buffer_masked(interp_type, data_var, sbuffer, south_fine, south_coarse, SOUTH, x_refine, y_refine, wt, ind, mask_var, mask_val, default_val)
      call fill_nest_from_buffer_masked(interp_type, data_var, ebuffer, east_fine, east_coarse, EAST, x_refine, y_refine, wt, ind, mask_var, mask_val, default_val)
      call fill_nest_from_buffer_masked(interp_type, data_var, wbuffer, west_fine, west_coarse, WEST, x_refine, y_refine, wt, ind, mask_var, mask_val, default_val)

    endif

    deallocate(nbuffer)
    deallocate(sbuffer)
    deallocate(ebuffer)
    deallocate(wbuffer)

  end subroutine fill_nest_halos_from_parent_masked


  subroutine fill_nest_halos_from_parent_r4_3d(var_name, data_var, interp_type, wt, ind, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)
    character(len=*), intent(in)                :: var_name
    real*4, allocatable, intent(inout)          :: data_var(:,:,:)
    integer, intent(in)                         :: interp_type
    real, allocatable, intent(in)               :: wt(:,:,:)
    integer, allocatable, intent(in)            :: ind(:,:,:)
    integer, intent(in)                         :: x_refine, y_refine
    logical, intent(in)                         :: is_fine_pe
    type(nest_domain_type), intent(inout)       :: nest_domain
    integer, intent(in)                         :: position, nz

    real*4, dimension(:,:,:), allocatable :: nbuffer, sbuffer, ebuffer, wbuffer
    type(bbox)                          :: north_fine, north_coarse
    type(bbox)                          :: south_fine, south_coarse
    type(bbox)                          :: east_fine, east_coarse
    type(bbox)                          :: west_fine, west_coarse
    integer                             :: this_pe
    integer                             :: nest_level = 1  ! TODO allow to vary

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

    ! Passes data from coarse grid to fine grid's halo
    call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

    if (is_fine_pe) then

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

  end subroutine fill_nest_halos_from_parent_r4_3d


  subroutine fill_nest_halos_from_parent_r8_3d(var_name, data_var, interp_type, wt, ind, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)
    character(len=*), intent(in)                :: var_name
    real*8, allocatable, intent(inout)          :: data_var(:,:,:)
    integer, intent(in)                         :: interp_type
    real, allocatable, intent(in)               :: wt(:,:,:)   ! TODO should this be real*8?
    integer, allocatable, intent(in)            :: ind(:,:,:)
    integer, intent(in)                         :: x_refine, y_refine
    logical, intent(in)                         :: is_fine_pe
    type(nest_domain_type), intent(inout)       :: nest_domain
    integer, intent(in)                         :: position, nz

    real*8, dimension(:,:,:), allocatable :: nbuffer, sbuffer, ebuffer, wbuffer
    type(bbox)                          :: north_fine, north_coarse
    type(bbox)                          :: south_fine, south_coarse
    type(bbox)                          :: east_fine, east_coarse
    type(bbox)                          :: west_fine, west_coarse
    integer                             :: this_pe
    integer                             :: nest_level = 1  ! TODO allow to vary

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

    ! Passes data from coarse grid to fine grid's halo
    call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

    if (is_fine_pe) then

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

  end subroutine fill_nest_halos_from_parent_r8_3d


  subroutine fill_nest_halos_from_parent_r4_4d(var_name, data_var, interp_type, wt, ind, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)
    character(len=*), intent(in)                :: var_name
    real*4, allocatable, intent(inout)          :: data_var(:,:,:,:)
    integer, intent(in)                         :: interp_type
    real, allocatable, intent(in)               :: wt(:,:,:)
    integer, allocatable, intent(in)            :: ind(:,:,:)
    integer, intent(in)                         :: x_refine, y_refine
    logical, intent(in)                         :: is_fine_pe
    type(nest_domain_type), intent(inout)       :: nest_domain
    integer, intent(in)                         :: position, nz

    real*4, dimension(:,:,:,:), allocatable :: nbuffer, sbuffer, ebuffer, wbuffer
    type(bbox)                            :: north_fine, north_coarse
    type(bbox)                            :: south_fine, south_coarse
    type(bbox)                            :: east_fine, east_coarse
    type(bbox)                            :: west_fine, west_coarse
    integer                               :: n4d, this_pe
    integer                               :: nest_level = 1  ! TODO allow to vary

    this_pe = mpp_pe()

    !!===========================================================
    !!
    !! Fill halo buffers
    !!
    !!===========================================================

    n4d = ubound(data_var, 4)

    call alloc_halo_buffer(nbuffer, north_fine, north_coarse, nest_domain, NORTH,  position, nz, n4d)
    call alloc_halo_buffer(sbuffer, south_fine, south_coarse, nest_domain, SOUTH,  position, nz, n4d)
    call alloc_halo_buffer(ebuffer, east_fine,  east_coarse,  nest_domain, EAST,   position, nz, n4d)
    call alloc_halo_buffer(wbuffer, west_fine,  west_coarse,  nest_domain, WEST,   position, nz, n4d)

    !====================================================
    ! Passes data from coarse grid to fine grid's halo
    ! Coarse parent PEs send data from data_var
    ! Fine halo PEs receive data into one or more of the halo buffers
    !====================================================

    call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

    if (is_fine_pe) then

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

  end subroutine fill_nest_halos_from_parent_r4_4d


  subroutine fill_nest_halos_from_parent_r8_4d(var_name, data_var, interp_type, wt, ind, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)
    character(len=*), intent(in)                :: var_name
    real*8, allocatable, intent(inout)          :: data_var(:,:,:,:)
    integer, intent(in)                         :: interp_type
    real, allocatable, intent(in)               :: wt(:,:,:) ! TODO should this be real*8?
    integer, allocatable, intent(in)            :: ind(:,:,:)
    integer, intent(in)                         :: x_refine, y_refine
    logical, intent(in)                         :: is_fine_pe
    type(nest_domain_type), intent(inout)       :: nest_domain
    integer, intent(in)                         :: position, nz

    real*8, dimension(:,:,:,:), allocatable :: nbuffer, sbuffer, ebuffer, wbuffer
    type(bbox)                            :: north_fine, north_coarse
    type(bbox)                            :: south_fine, south_coarse
    type(bbox)                            :: east_fine, east_coarse
    type(bbox)                            :: west_fine, west_coarse
    integer                               :: n4d, this_pe
    integer                               :: nest_level = 1  ! TODO allow to vary

    this_pe = mpp_pe()

    !!===========================================================
    !!
    !! Fill halo buffers
    !!
    !!===========================================================

    n4d = ubound(data_var, 4)

    call alloc_halo_buffer(nbuffer, north_fine, north_coarse, nest_domain, NORTH,  position, nz, n4d)
    call alloc_halo_buffer(sbuffer, south_fine, south_coarse, nest_domain, SOUTH,  position, nz, n4d)
    call alloc_halo_buffer(ebuffer, east_fine,  east_coarse,  nest_domain, EAST,   position, nz, n4d)
    call alloc_halo_buffer(wbuffer, west_fine,  west_coarse,  nest_domain, WEST,   position, nz, n4d)

    !====================================================
    ! Passes data from coarse grid to fine grid's halo
    ! Coarse parent PEs send data from data_var
    ! Fine halo PEs receive data into one or more of the halo buffers
    !====================================================

    call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

    if (is_fine_pe) then

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

  end subroutine fill_nest_halos_from_parent_r8_4d


  !==================================================================================================
  !
  !  Allocate halo buffers
  !
  !==================================================================================================

  subroutine alloc_halo_buffer_r8_2d(buffer, bbox_fine, bbox_coarse, nest_domain, direction, position)
    real*8, dimension(:,:), allocatable, intent(out) :: buffer
    type(bbox), intent(out)                          :: bbox_fine, bbox_coarse
    type(nest_domain_type), intent(in)               :: nest_domain
    integer, intent(in)                              :: direction, position

    call bbox_get_C2F_index(nest_domain, bbox_fine, bbox_coarse, direction,  position)

    if( bbox_coarse.ie .GE. bbox_coarse.is .AND. bbox_coarse.je .GE. bbox_coarse.js ) then
      allocate(buffer(bbox_coarse.is:bbox_coarse.ie, bbox_coarse.js:bbox_coarse.je))
    else
      ! The buffer must have some storage allocated, whether it's a useful buffer or just a dummy.
      allocate(buffer(1,1))
    endif

    buffer = 0

  end subroutine alloc_halo_buffer_r8_2d


  subroutine alloc_halo_buffer_r4_2d(buffer, bbox_fine, bbox_coarse, nest_domain, direction, position)
    real*4, dimension(:,:), allocatable, intent(out) :: buffer
    type(bbox), intent(out)                          :: bbox_fine, bbox_coarse
    type(nest_domain_type), intent(in)               :: nest_domain
    integer, intent(in)                              :: direction, position

    call bbox_get_C2F_index(nest_domain, bbox_fine, bbox_coarse, direction,  position)

    if( bbox_coarse.ie .GE. bbox_coarse.is .AND. bbox_coarse.je .GE. bbox_coarse.js ) then
      allocate(buffer(bbox_coarse.is:bbox_coarse.ie, bbox_coarse.js:bbox_coarse.je))
    else
      ! The buffer must have some storage allocated, whether it's a useful buffer or just a dummy.
      allocate(buffer(1,1))
    endif

    buffer = 0

  end subroutine alloc_halo_buffer_r4_2d


  subroutine alloc_halo_buffer_r4_3d(buffer, bbox_fine, bbox_coarse, nest_domain, direction, position, nz)
    real*4, dimension(:,:,:), allocatable, intent(out) :: buffer
    type(bbox), intent(out)                            :: bbox_fine, bbox_coarse
    type(nest_domain_type), intent(in)                 :: nest_domain
    integer, intent(in)                                :: direction, position, nz


    call bbox_get_C2F_index(nest_domain, bbox_fine, bbox_coarse, direction,  position)

    if( bbox_coarse.ie .GE. bbox_coarse.is .AND. bbox_coarse.je .GE. bbox_coarse.js ) then
      allocate(buffer(bbox_coarse.is:bbox_coarse.ie, bbox_coarse.js:bbox_coarse.je,1:nz))
    else
      ! The buffer must have some storage allocated, whether it's a useful buffer or just a dummy.
      allocate(buffer(1,1,1))
    endif

    buffer = 0

  end subroutine alloc_halo_buffer_r4_3d


  subroutine alloc_halo_buffer_r8_3d(buffer, bbox_fine, bbox_coarse, nest_domain, direction, position, nz)
    real*8, dimension(:,:,:), allocatable, intent(out) :: buffer
    type(bbox), intent(out)                            :: bbox_fine, bbox_coarse
    type(nest_domain_type), intent(in)                 :: nest_domain
    integer, intent(in)                                :: direction, position, nz

    call bbox_get_C2F_index(nest_domain, bbox_fine, bbox_coarse, direction,  position)

    if( bbox_coarse.ie .GE. bbox_coarse.is .AND. bbox_coarse.je .GE. bbox_coarse.js ) then
      allocate(buffer(bbox_coarse.is:bbox_coarse.ie, bbox_coarse.js:bbox_coarse.je,1:nz))
    else
      ! The buffer must have some storage allocated, whether it's a useful buffer or just a dummy.
      allocate(buffer(1,1,1))
    endif

    buffer = 0

  end subroutine alloc_halo_buffer_r8_3d


  subroutine alloc_halo_buffer_r4_4d(buffer, bbox_fine, bbox_coarse, nest_domain, direction, position, nz, n4d)
    real*4, dimension(:,:,:,:), allocatable, intent(out) :: buffer
    type(bbox), intent(out)                              :: bbox_fine, bbox_coarse
    type(nest_domain_type), intent(in)                   :: nest_domain
    integer, intent(in)                                  :: direction, position, nz, n4d

    call bbox_get_C2F_index(nest_domain, bbox_fine, bbox_coarse, direction,  position)

    if( bbox_coarse.ie .GE. bbox_coarse.is .AND. bbox_coarse.je .GE. bbox_coarse.js ) then
      allocate(buffer(bbox_coarse.is:bbox_coarse.ie, bbox_coarse.js:bbox_coarse.je, 1:nz, 1:n4d))
    else
      ! The buffer must have some storage allocated, whether it's a useful buffer or just a dummy.
      allocate(buffer(1,1,1,1))
    endif

    buffer = 0

  end subroutine alloc_halo_buffer_r4_4d


  subroutine alloc_halo_buffer_r8_4d(buffer, bbox_fine, bbox_coarse, nest_domain, direction, position, nz, n4d)
    real*8, dimension(:,:,:,:), allocatable, intent(out) :: buffer
    type(bbox), intent(out)                              :: bbox_fine, bbox_coarse
    type(nest_domain_type), intent(in)                   :: nest_domain
    integer, intent(in)                                  :: direction, position, nz, n4d

    call bbox_get_C2F_index(nest_domain, bbox_fine, bbox_coarse, direction,  position)

    if( bbox_coarse.ie .GE. bbox_coarse.is .AND. bbox_coarse.je .GE. bbox_coarse.js ) then
      allocate(buffer(bbox_coarse.is:bbox_coarse.ie, bbox_coarse.js:bbox_coarse.je, 1:nz, 1:n4d))
    else
      ! The buffer must have some storage allocated, whether it's a useful buffer or just a dummy.
      allocate(buffer(1,1,1,1))
    endif

    buffer = 0

  end subroutine alloc_halo_buffer_r8_4d


  !==================================================================================================
  !
  !  Load static data from netCDF files
  !
  !==================================================================================================

  ! Load the full panel nest latlons from netCDF file
  ! character(*), parameter      :: nc_filename = '/scratch2/NAGAPE/aoml-hafs1/William.Ramstrom/static_grids/C384_grid.tile6.nc'
  ! Read in the lat/lon in degrees, convert to radians

  subroutine load_nest_latlons_from_nc(nc_filename, nxp, nyp, refine, pelist, &
      fp_tile_geo, fp_istart_fine, fp_iend_fine, fp_jstart_fine, fp_jend_fine)
    implicit none

    character(*), intent(in)              :: nc_filename
    integer, intent(in)                   :: nxp, nyp, refine
    integer, allocatable, intent(in)      :: pelist(:)
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
    real(kind=kind_phys) :: deg2rad

    deg2rad = pi / 180.0d0

    this_pe = mpp_pe()

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

    call alloc_read_data(nc_filename, 'x', super_nxp, super_nyp, fp_tile_geo%lons, pelist)
    call alloc_read_data(nc_filename, 'y', super_nxp, super_nyp, fp_tile_geo%lats, pelist)
    call alloc_read_data(nc_filename, 'area', super_nx, super_ny, fp_tile_geo%area, pelist)

    !  double dx(nyp, nx)
    !call alloc_read_data(nc_filename, 'dx', super_nx, super_nyp, fp_tile_geo%dx)
    ! double dy(ny, nxp)
    !call alloc_read_data(nc_filename, 'dy', super_nxp, super_ny, fp_tile_geo%dy)
    ! double area(ny, nx)
    !call alloc_read_data(nc_filename, 'area', super_nx, super_ny, fp_tile_geo%area)


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

    fp_tile_geo%lats = fp_tile_geo%lats * deg2rad
    fp_tile_geo%lons = fp_tile_geo%lons * deg2rad

  end subroutine load_nest_latlons_from_nc

#ifdef OVERLOAD_R8
  subroutine alloc_read_data_r4_2d(nc_filename, var_name, x_size, y_size, data_array, pes, time)
    character(len=*), intent(in)           :: nc_filename, var_name
    integer, intent(in)                    :: x_size, y_size
    real*4, allocatable, intent(inout)     :: data_array(:,:)
    integer, allocatable, intent(in)       :: pes(:)
    integer, intent(in),optional           :: time

    type(FmsNetcdfFile_t)        :: fileobj        !< Fms2_io fileobj
    real*4, allocatable          :: time_array(:,:,:)
    integer                      :: this_pe

    ! Allocate data_array to match the expected data size, then read in the data
    ! This subroutine consolidates the allocation and reading of data to ensure consistency of data sizing and simplify code
    ! Could later extend this function to determine data size based on netCDF file metadata

    this_pe = mpp_pe()

    allocate(data_array(x_size, y_size))
    data_array = -9999.9

    if (present(time)) then
      allocate(time_array(x_size, y_size, 12)) ! assume monthly data; allocate 12 slots
      if (open_file(fileobj, nc_filename, "read", pelist=pes, is_restart=.false.)) then
        call read_data(fileobj, var_name, time_array)
        call close_file(fileobj)
      endif
      
      data_array = time_array(:,:,time)
      deallocate(time_array)
    else
      ! Following transition documents at https://github.com/NOAA-GFDL/FMS/tree/2021.03.01/fms2_io
      if (open_file(fileobj, nc_filename, "read", pelist=pes, is_restart=.false.)) then
        call read_data(fileobj, var_name, data_array)
        call close_file(fileobj)
      endif
    endif
    
  end subroutine alloc_read_data_r4_2d
#endif

  subroutine alloc_read_data_r8_2d(nc_filename, var_name, x_size, y_size, data_array, pes, time)
    character(len=*), intent(in)           :: nc_filename, var_name
    integer, intent(in)                    :: x_size, y_size
    real*8, allocatable, intent(inout)     :: data_array(:,:)
    integer, allocatable, intent(in)       :: pes(:)
    integer, intent(in),optional           :: time

    real*8, allocatable          :: time_array(:,:,:)
    type(FmsNetcdfFile_t)        :: fileobj        !< Fms2_io fileobj
    integer                      :: this_pe

    ! Allocate data_array to match the expected data size, then read in the data
    ! This subroutine consolidates the allocation and reading of data to ensure consistency of data sizing and simplify code
    ! Could later extend this function to determine data size based on netCDF file metadata

    this_pe = mpp_pe()

    allocate(data_array(x_size, y_size))
    data_array = -9999.9

    ! Following transition documents at https://github.com/NOAA-GFDL/FMS/tree/2021.03.01/fms2_io
    if (present(time)) then
      allocate(time_array(x_size, y_size, 12)) ! assume monthly data; allocate 12 slots
      if (open_file(fileobj, nc_filename, "read", pelist=pes, is_restart=.false.)) then
        call read_data(fileobj, var_name, time_array)
        call close_file(fileobj)
      endif
      
      data_array = time_array(:,:,time)
      deallocate(time_array)
    else
      if (open_file(fileobj, nc_filename, "read", pelist=pes, is_restart=.false.)) then
        call read_data(fileobj, var_name, data_array)
        call close_file(fileobj)
      endif
    endif
    
  end subroutine alloc_read_data_r8_2d


  !==================================================================================================
  !
  !  NetCDF Function Section
  !
  !==================================================================================================

  subroutine output_grid_to_nc_3d(flag, istart, iend, jstart, jend, k, grid, file_str, var_name, time_step, dom, pos)
    implicit none

    character(len=*), intent(in)  :: flag
    integer, intent(in)           :: istart, iend, jstart, jend, k
    real, dimension(:,:,:), intent(in)   :: grid
    character(len=*), intent(in)  :: file_str, var_name
    integer, intent(in)           :: time_step
    type(domain2d), intent(in)    :: dom
    integer, intent(in)           :: pos

    logical                     :: new_file
    integer                     :: this_pe
    character(len=512)          :: dirname
    character(len=512)          :: filename
    type(FmsNetcdfDomainFile_t) :: fileobj        !< Fms2_io domain decomposed fileobj
    character(len=10)           :: dim_names(3)   !< Array of dimension names
    integer                     :: istat
    logical                     :: file_exists
    character(len=12)           :: mode

    istat = getcwd(dirname)
    write (filename, "(A,A1,A,A1,A,A1,I0.3,A)") trim(dirname), "/", trim(file_str), "_", trim(var_name), "_", time_step, ".nc"

   if (pos .eq. CENTER) then
      dim_names(1) = "xaxis_1"
      dim_names(2) = "yaxis_1"
    elseif (pos .eq. NORTH) then
      dim_names(1) = "xaxis_2"
      dim_names(2) = "yaxis_2"
    elseif (pos .eq. EAST) then
      dim_names(1) = "xaxis_3"
      dim_names(2) = "yaxis_3"
    endif

    !dim_names(3) = "zaxis_1"
    write (dim_names(3),'(A,I0)') "zaxis_", k
    
    !inquire(FILE=filename, EXIST=file_exists)
    !if (file_exists) then
    !  mode = "append"
    !else
    !  mode = "overwrite"
    !endif

    new_file = .true.

    if (new_file) then
      mode = "write"
    else
      mode = "append"
    endif

    mode = "write"

    if (open_file(fileobj, filename, mode, dom)) then
      
      if (new_file) then
        call register_axis(fileobj, dim_names(1), "x", CENTER)  ! TODO investigate handling of non-centered position
        call register_axis(fileobj, dim_names(2), "y", CENTER)  ! TODO investigate handling of non-centered position
        call register_axis(fileobj, trim(dim_names(3)), k)
      endif

      call register_field(fileobj, trim(var_name), 'float', dim_names)
      call write_data(fileobj, trim(var_name), grid)
      call close_file(fileobj)
    endif

!      if (.not. is_dimension_registered(fileobj, dim_names(1))) then
!        call register_axis(fileobj, dim_names(1), "x")  ! TODO investigate handling of non-centered position
!      endif
!     if (.not. is_dimension_registered(fileobj, dim_names(2))) call register_axis(fileobj, dim_names(2), "y")  ! TODO investigate handling of non-centered position
!      if (.not. is_dimension_registered(fileobj, trim(dim_names(3)))) then
!        call register_axis(fileobj, trim(dim_names(3)), k)
!      endif
 
  end subroutine output_grid_to_nc_3d


  subroutine output_grid_to_nc_2d(flag, istart, iend, jstart, jend, grid, file_str, var_name, time_step, dom, pos)
    implicit none

    character(len=*), intent(in)  :: flag
    integer, intent(in)           :: istart, iend, jstart, jend
    real, dimension(:,:), intent(in)   :: grid
    character(len=*), intent(in)  :: file_str, var_name
    integer, intent(in)           :: time_step
    type(domain2d), intent(in)    :: dom
    integer, intent(in)           :: pos

    logical                     :: new_file
    integer                     :: istat
    character(len=512)          :: dirname
    character(len=512)          :: filename
    type(FmsNetcdfDomainFile_t) :: fileobj        !< Fms2_io domain decomposed fileobj                                                               
    character(len=8)            :: dim_names(2)   !< Array of dimension names                                                                        
    character(len=12)           :: mode

    istat = getcwd(dirname)
    write (filename, "(A,A1,A,A1,A,A1,I0.3,A)") trim(dirname), "/", trim(file_str), "_", trim(var_name), "_", time_step, ".nc"

    if (pos .eq. CENTER) then
      dim_names(1) = "xaxis_1"
      dim_names(2) = "yaxis_1"
    elseif (pos .eq. NORTH) then
      dim_names(1) = "xaxis_2"
      dim_names(2) = "yaxis_2"
    elseif (pos .eq. EAST) then
      dim_names(1) = "xaxis_3"
      dim_names(2) = "yaxis_3"
    endif

    new_file = .true.

    if (new_file) then
      mode = "write"
    else
      mode = "append"
    endif

    if (open_file(fileobj, filename, mode, dom)) then
      if (new_file) then
        call register_axis(fileobj, dim_names(1), "x", CENTER)  ! TODO investigate handling of non-centered position                                 
        call register_axis(fileobj, dim_names(2), "y", CENTER)  ! TODO investigate handling of non-centered position                                 
      endif

      call register_field(fileobj, trim(var_name), 'float', dim_names)
      call write_data(fileobj, trim(var_name), grid)
      call close_file(fileobj)
    endif

  end subroutine output_grid_to_nc_2d



  !==================================================================================================
  !
  !  Fill Section
  !
  !==================================================================================================

  subroutine fill_grid_from_supergrid_r4_3d(in_grid, stagger_type, fp_super_tile_geo, ioffset, joffset, x_refine, y_refine)
    implicit none
    real*4, allocatable, intent(inout)  :: in_grid(:,:,:)
    integer, intent(in)                 :: stagger_type   ! CENTER, CORNER
    type(grid_geometry), intent(in)     :: fp_super_tile_geo
    integer, intent(in)                 :: ioffset, joffset, x_refine, y_refine

    integer           :: nest_x, nest_y, parent_x, parent_y
    type(bbox)        :: tile_bbox, fp_tile_bbox
    integer           :: i, j, fp_i, fp_j
    character(len=64) :: errstring

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
        endif

        ! Make sure we don't run off the edge of the parent supergrid
        if (fp_i < fp_tile_bbox%is .or. fp_i > fp_tile_bbox%ie) then
          write(errstring, "(A,I0,A,I0,A,I0)") "fp_i=", fp_i," is=",fp_tile_bbox%is," ie=",fp_tile_bbox%ie
          call mpp_error(FATAL, "fill_grid_from_supergrid_r4_3d invalid bounds i " // errstring)
        endif
        if (fp_j < fp_tile_bbox%js .or. fp_j > fp_tile_bbox%je) then
          write(errstring, "(A,I0,A,I0,A,I0)") "fp_j=", fp_j," js=",fp_tile_bbox%js," je=",fp_tile_bbox%je
          call mpp_error(FATAL, "fill_grid_from_supergrid_r4_3d invalid bounds j " // errstring)
        endif

        in_grid(i,j,2) = fp_super_tile_geo%lats(fp_i, fp_j)
        in_grid(i,j,1) = fp_super_tile_geo%lons(fp_i, fp_j)
      enddo
    enddo

    ! Validate at the end
    !call find_nest_alignment(tile_geo, fp_super_tile_geo, nest_x, nest_y, parent_x, parent_y)

  end subroutine fill_grid_from_supergrid_r4_3d


  subroutine fill_grid_from_supergrid_r8_3d(in_grid, stagger_type, fp_super_tile_geo, ioffset, joffset, x_refine, y_refine)
    implicit none
    real*8, allocatable, intent(inout)  :: in_grid(:,:,:)
    integer, intent(in)                 :: stagger_type   ! CENTER, CORNER
    type(grid_geometry), intent(in)     :: fp_super_tile_geo
    integer, intent(in)                 :: ioffset, joffset, x_refine, y_refine

    integer           :: nest_x, nest_y, parent_x, parent_y
    type(bbox)        :: tile_bbox, fp_tile_bbox
    integer           :: i, j, fp_i, fp_j
    character(len=64) :: errstring

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
        endif

        ! Make sure we don't run off the edge of the parent supergrid
        if (fp_i < fp_tile_bbox%is .or. fp_i > fp_tile_bbox%ie) then
          write(errstring, "(A,I0,A,I0,A,I0)") "fp_i=", fp_i," is=",fp_tile_bbox%is," ie=",fp_tile_bbox%ie
          call mpp_error(FATAL, "fill_grid_from_supergrid_r8_3d invalid bounds i " // errstring)
        endif
        if (fp_j < fp_tile_bbox%js .or. fp_j > fp_tile_bbox%je) then
          write(errstring, "(A,I0,A,I0,A,I0)") "fp_j=", fp_j," js=",fp_tile_bbox%js," je=",fp_tile_bbox%je
          call mpp_error(FATAL, "fill_grid_from_supergrid_r8_3d invalid bounds j " // errstring)
        endif

        in_grid(i,j,2) = fp_super_tile_geo%lats(fp_i, fp_j)
        in_grid(i,j,1) = fp_super_tile_geo%lons(fp_i, fp_j)
      enddo
    enddo

    ! Validate at the end
    !call find_nest_alignment(tile_geo, fp_super_tile_geo, nest_x, nest_y, parent_x, parent_y)

  end subroutine fill_grid_from_supergrid_r8_3d


  subroutine fill_grid_from_supergrid_r8_4d(in_grid, stagger_type, fp_super_tile_geo, ioffset, joffset, x_refine, y_refine)
    implicit none
    real*8, allocatable, intent(inout)  :: in_grid(:,:,:,:)
    integer, intent(in)                 :: stagger_type   ! CENTER, CORNER
    type(grid_geometry), intent(in)     :: fp_super_tile_geo
    integer, intent(in)                 :: ioffset, joffset, x_refine, y_refine

    integer           :: nest_x, nest_y, parent_x, parent_y
    type(bbox)        :: tile_bbox, fp_tile_bbox
    integer           :: i, j, fp_i, fp_j
    character(len=64) :: errstring

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
        endif

        ! Make sure we don't run off the edge of the parent supergrid
        if (fp_i < fp_tile_bbox%is .or. fp_i > fp_tile_bbox%ie) then
          write(errstring, "(A,I0,A,I0,A,I0)") "fp_i=", fp_i," is=",fp_tile_bbox%is," ie=",fp_tile_bbox%ie
          call mpp_error(FATAL, "fill_grid_from_supergrid_r8_4d invalid bounds i " // errstring)
        endif
        if (fp_j < fp_tile_bbox%js .or. fp_j > fp_tile_bbox%je) then
          write(errstring, "(A,I0,A,I0,A,I0)") "fp_j=", fp_j," js=",fp_tile_bbox%js," je=",fp_tile_bbox%je
          call mpp_error(FATAL, "fill_grid_from_supergrid_r8_4d invalid bounds j " // errstring)
        endif

        in_grid(i,j,2,1) = fp_super_tile_geo%lats(fp_i, fp_j)
        in_grid(i,j,1,1) = fp_super_tile_geo%lons(fp_i, fp_j)
      enddo
    enddo

    ! Validate at the end
    !call find_nest_alignment(tile_geo, fp_super_tile_geo, nest_x, nest_y, parent_x, parent_y)

  end subroutine fill_grid_from_supergrid_r8_4d


  !>@brief  This subroutine fills the nest halo data from the coarse grid data by downscaling.
  !>@details  Applicable to any interpolation type

  subroutine fill_nest_from_buffer_r4_2d(interp_type, x, buffer, bbox_fine, bbox_coarse, dir, x_refine, y_refine, wt, ind)
    implicit none

    integer, intent(in)                         :: interp_type
    real*4,  allocatable, intent(inout)         :: x(:,:)
    real*4,  allocatable, intent(in)            :: buffer(:,:)
    type(bbox), intent(in)                      :: bbox_fine, bbox_coarse
    integer, intent(in)                         :: dir, x_refine, y_refine
    real, allocatable, intent(in)               :: wt(:,:,:)    ! The final dimension is always 4
    integer, allocatable, intent(in)            :: ind(:,:,:)

    integer   :: this_pe
    this_pe = mpp_pe()

    ! Output the interpolation type
    select case (interp_type)
    case (1)
      call fill_nest_from_buffer_cell_center("A", x, buffer, bbox_fine, bbox_coarse, dir, x_refine, y_refine, wt, ind)
      !     case (3)    ! C grid staggered
    case (4)
      call fill_nest_from_buffer_cell_center("D", x, buffer, bbox_fine, bbox_coarse, dir, x_refine, y_refine, wt, ind)
    case (9)
      !call fill_nest_from_buffer_nearest_neighbor(x, buffer, bbox_fine, bbox_coarse, dir, wt)
      call mpp_error(FATAL, '2D fill_nest_from_buffer_nearest_neighbor not yet implemented.')
    case default
      call mpp_error(FATAL, 'interp_single_nest got invalid value for interp_type from namelist.')
    end select

  end subroutine fill_nest_from_buffer_r4_2d


  !>@brief  This subroutine fills the nest halo data from the coarse grid data by downscaling.
  !>@details  Applicable to any interpolation type

  subroutine fill_nest_from_buffer_r8_2d(interp_type, x, buffer, bbox_fine, bbox_coarse, dir, x_refine, y_refine, wt, ind)
    implicit none

    integer, intent(in)                         :: interp_type
    real*8, allocatable, intent(inout)          :: x(:,:)
    real*8, allocatable, intent(in)             :: buffer(:,:)
    type(bbox), intent(in)                      :: bbox_fine, bbox_coarse
    integer, intent(in)                         :: dir, x_refine, y_refine
    real, allocatable, intent(in)               :: wt(:,:,:)    ! The final dimension is always 4
    integer, allocatable, intent(in)            :: ind(:,:,:)

    integer   :: this_pe
    this_pe = mpp_pe()

    ! Output the interpolation type
    select case (interp_type)
    case (1)
      call fill_nest_from_buffer_cell_center("A", x, buffer, bbox_fine, bbox_coarse, dir, x_refine, y_refine, wt, ind)
      !     case (3)   ! C grid staggered
    case (4)
      call fill_nest_from_buffer_cell_center("D", x, buffer, bbox_fine, bbox_coarse, dir, x_refine, y_refine, wt, ind)
    case (9)
      !call fill_nest_from_buffer_nearest_neighbor(x, buffer, bbox_fine, bbox_coarse, dir, wt)
      call mpp_error(FATAL, '2D fill_nest_from_buffer_nearest_neighbor not yet implemented.')
    case default
      call mpp_error(FATAL, 'interp_single_nest got invalid value for interp_type from namelist.')
    end select

  end subroutine fill_nest_from_buffer_r8_2d


  subroutine fill_nest_from_buffer_masked(interp_type, x, buffer, bbox_fine, bbox_coarse, dir, x_refine, y_refine, wt, ind, mask_var, mask_val, default_val)
    implicit none

    integer, intent(in)                         :: interp_type
    real*8, allocatable, intent(inout)          :: x(:,:)
    real*8, allocatable, intent(in)             :: buffer(:,:)
    type(bbox), intent(in)                      :: bbox_fine, bbox_coarse
    integer, intent(in)                         :: dir, x_refine, y_refine
    real, allocatable, intent(in)               :: wt(:,:,:)    ! The final dimension is always 4
    integer, allocatable, intent(in)            :: ind(:,:,:)
    real, allocatable, intent(in)               :: mask_var(:,:)
    integer, intent(in)                         :: mask_val
    real*8, intent(in)                          :: default_val

    integer   :: this_pe
    this_pe = mpp_pe()

    ! Output the interpolation type
    select case (interp_type)
    case (1)
      call fill_nest_from_buffer_cell_center("A", x, buffer, bbox_fine, bbox_coarse, dir, x_refine, y_refine, wt, ind)
      !     case (3)  ! C grid staggered
    case (4)
      call fill_nest_from_buffer_cell_center("D", x, buffer, bbox_fine, bbox_coarse, dir, x_refine, y_refine, wt, ind)
    case (7)
      call fill_nest_from_buffer_cell_center_masked("A", x, buffer, bbox_fine, bbox_coarse, dir, x_refine, y_refine, wt, ind, mask_var, mask_val, default_val)
    case (9)
      !call fill_nest_from_buffer_nearest_neighbor(x, buffer, bbox_fine, bbox_coarse, dir, wt)
      call mpp_error(FATAL, '2D fill_nest_from_buffer_nearest_neighbor not yet implemented.')
    case default
      call mpp_error(FATAL, 'interp_single_nest got invalid value for interp_type from namelist.')
    end select

  end subroutine fill_nest_from_buffer_masked



  subroutine fill_nest_from_buffer_r4_3d(interp_type, x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
    implicit none

    integer, intent(in)                         :: interp_type
    real*4, allocatable, intent(inout)          :: x(:,:,:)
    real*4, allocatable, intent(in)             :: buffer(:,:,:)
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
      call fill_nest_from_buffer_cell_center("A", x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
      !     case (3) ! C grid staggered
    case (4)
      call fill_nest_from_buffer_cell_center("D", x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
    case (9)
      !call fill_nest_from_buffer_nearest_neighbor(x, buffer, bbox_fine, bbox_coarse, nz, dir, wt)
      call mpp_error(FATAL, 'fill_nest_from_buffer_nearest_neighbor is not yet implemented.')
    case default
      call mpp_error(FATAL, 'interp_single_nest got invalid value for interp_type from namelist.')
    end select

  end subroutine fill_nest_from_buffer_r4_3d


  subroutine fill_nest_from_buffer_r8_3d(interp_type, x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
    implicit none

    integer, intent(in)                         :: interp_type
    real*8, allocatable, intent(inout)          :: x(:,:,:)
    real*8, allocatable, intent(in)             :: buffer(:,:,:)
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
      call fill_nest_from_buffer_cell_center("A", x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
      !     case (3)  ! C grid staggered
    case (4)
      call fill_nest_from_buffer_cell_center("D", x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
    case (9)
      call mpp_error(FATAL, 'nearest_neighbor is not yet implemented for fv_moving_nest_utils.F90::fill_nest_from_buffer_3D_kindphys')
      !call fill_nest_from_buffer_nearest_neighbor(x, buffer, bbox_fine, bbox_coarse, nz, dir, wt)
    case default
      call mpp_error(FATAL, 'interp_single_nest got invalid value for interp_type from namelist.')
    end select

  end subroutine fill_nest_from_buffer_r8_3d


  !>@brief  This subroutine fills the nest halo data from the coarse grid data by downscaling.
  !>@details  Applicable to any interpolation type

  subroutine fill_nest_from_buffer_r4_4d(interp_type, x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
    implicit none

    integer, intent(in)                         :: interp_type
    real*4, allocatable, intent(inout)          :: x(:,:,:,:)
    real*4, allocatable, intent(in)             :: buffer(:,:,:,:)
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
      call fill_nest_from_buffer_cell_center("A", x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
      !     case (3)  ! C grid staggered
    case (4)
      call fill_nest_from_buffer_cell_center("D", x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
    case (9)
      !call fill_nest_from_buffer_nearest_neighbor(x, buffer, bbox_fine, bbox_coarse, nz, dir, wt)
      call mpp_error(FATAL, '4D fill_nest_from_buffer_nearest_neighbor not yet implemented.')
    case default
      call mpp_error(FATAL, 'interp_single_nest got invalid value for interp_type from namelist.')
    end select

  end subroutine fill_nest_from_buffer_r4_4d


  subroutine fill_nest_from_buffer_r8_4d(interp_type, x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
    implicit none

    integer, intent(in)                         :: interp_type
    real*8, allocatable, intent(inout)          :: x(:,:,:,:)
    real*8, allocatable, intent(in)             :: buffer(:,:,:,:)
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
      call fill_nest_from_buffer_cell_center("A", x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
      !     case (3) ! C grid staggered
    case (4)
      call fill_nest_from_buffer_cell_center("D", x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
    case (9)
      !call fill_nest_from_buffer_nearest_neighbor(x, buffer, bbox_fine, bbox_coarse, nz, dir, wt)
      call mpp_error(FATAL, '4D fill_nest_from_buffer_nearest_neighbor not yet implemented.')
    case default
      call mpp_error(FATAL, 'interp_single_nest got invalid value for interp_type from namelist.')
    end select

  end subroutine fill_nest_from_buffer_r8_4d


  !>@brief  This subroutine fills the nest halo data from the coarse grid data by downscaling.  It can accommodate all grid staggers, using the stagger variable.  [The routine needs to be renamed since "_from_cell_center" has become incorrect.)
  !>@details  Applicable to any interpolation type

  subroutine fill_nest_from_buffer_cell_center_r4_2d(stagger, x, buffer, bbox_fine, bbox_coarse, dir, x_refine, y_refine, wt, ind)
    implicit none
    character ( len = 1 ), intent(in)             :: stagger
    real*4, allocatable, intent(inout)            :: x(:,:)
    real*4, allocatable, intent(in)               :: buffer(:,:)
    type(bbox), intent(in)                        :: bbox_fine, bbox_coarse
    integer, intent(in)                           :: dir, x_refine, y_refine
    real, allocatable, intent(in)                 :: wt(:,:,:)    ! The final dimension is always 4
    integer, allocatable, intent(in)              :: ind(:,:,:)

    character(len=8)       :: dir_str
    integer                :: i, j, k, ic, jc

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
      do j=bbox_fine%js, bbox_fine%je
        do i=bbox_fine%is, bbox_fine%ie
          !if (stagger == "A") then
          !else if (stagger == "C") then
          !else if (stagger == "D") then
          !endif

          ic = ind(i,j,1)
          jc = ind(i,j,2)

          x(i,j) = &
              wt(i,j,1)*buffer(ic,  jc  ) +  &
              wt(i,j,2)*buffer(ic,  jc+1) +  &
              wt(i,j,3)*buffer(ic+1,jc+1) +  &
              wt(i,j,4)*buffer(ic+1,jc  )

        enddo
      enddo
    endif

  end subroutine fill_nest_from_buffer_cell_center_r4_2d


  subroutine fill_nest_from_buffer_cell_center_r8_2d(stagger, x, buffer, bbox_fine, bbox_coarse, dir, x_refine, y_refine, wt, ind)
    implicit none
    character ( len = 1 ), intent(in)             :: stagger
    real*8, allocatable, intent(inout)            :: x(:,:)
    real*8, allocatable, intent(in)               :: buffer(:,:)
    type(bbox), intent(in)                        :: bbox_fine, bbox_coarse
    integer, intent(in)                           :: dir, x_refine, y_refine
    real, allocatable, intent(in)                 :: wt(:,:,:)    ! The final dimension is always 4
    integer, allocatable, intent(in)              :: ind(:,:,:)

    character(len=8)       :: dir_str
    integer                :: i, j, k, ic, jc

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
      do j=bbox_fine%js, bbox_fine%je
        do i=bbox_fine%is, bbox_fine%ie
          !if (stagger == "A") then
          !else if (stagger == "C") then
          !else if (stagger == "D") then
          !endif

          ic = ind(i,j,1)
          jc = ind(i,j,2)

          x(i,j) = &
              wt(i,j,1)*buffer(ic,  jc  ) +  &
              wt(i,j,2)*buffer(ic,  jc+1) +  &
              wt(i,j,3)*buffer(ic+1,jc+1) +  &
              wt(i,j,4)*buffer(ic+1,jc  )

        enddo
      enddo
    endif

  end subroutine fill_nest_from_buffer_cell_center_r8_2d


  subroutine fill_nest_from_buffer_cell_center_masked(stagger, x, buffer, bbox_fine, bbox_coarse, dir, x_refine, y_refine, wt, ind, mask_var, mask_val, default_val)
    implicit none
    character ( len = 1 ), intent(in)             :: stagger
    real*8, allocatable, intent(inout)            :: x(:,:)
    real*8, allocatable, intent(in)               :: buffer(:,:)
    type(bbox), intent(in)                        :: bbox_fine, bbox_coarse
    integer, intent(in)                           :: dir, x_refine, y_refine
    real, allocatable, intent(in)                 :: wt(:,:,:)    ! The final dimension is always 4
    integer, allocatable, intent(in)              :: ind(:,:,:)
    real, allocatable, intent(in)                 :: mask_var(:,:)
    integer, intent(in)                           :: mask_val
    real*8, intent(in)                            :: default_val

    character(len=8)       :: dir_str
    integer                :: i, j, k, ic, jc
    real                   :: tw

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
      do j=bbox_fine%js, bbox_fine%je
        do i=bbox_fine%is, bbox_fine%ie

          ic = ind(i,j,1)
          jc = ind(i,j,2)

          !x(i,j) = &
          !     wt(i,j,1)*buffer(ic,  jc  ) +  &
          !     wt(i,j,2)*buffer(ic,  jc+1) +  &
          !     wt(i,j,3)*buffer(ic+1,jc+1) +  &
          !     wt(i,j,4)*buffer(ic+1,jc  )

          ! Land type
          !if (mask_var(i,j) .eq. mask_val) then
          x(i,j) = 0.0
          tw = 0.0
          if (buffer(ic,jc) .gt. -1.0)     x(i,j) = x(i,j) + wt(i,j,1)*buffer(ic,  jc  )
          if (buffer(ic,jc+1) .gt. -1.0)   x(i,j) = x(i,j) + wt(i,j,1)*buffer(ic,  jc+1)
          if (buffer(ic+1,jc+1) .gt. -1.0) x(i,j) = x(i,j) + wt(i,j,1)*buffer(ic+1,jc+1)
          if (buffer(ic+1,jc) .gt. -1.0)   x(i,j) = x(i,j) + wt(i,j,1)*buffer(ic+1,jc  )

          if (buffer(ic,jc) .gt. -1.0)     tw = tw + wt(i,j,1)
          if (buffer(ic,jc+1) .gt. -1.0)   tw = tw + wt(i,j,1)
          if (buffer(ic+1,jc+1) .gt. -1.0) tw = tw + wt(i,j,1)
          if (buffer(ic+1,jc) .gt. -1.0)   tw = tw + wt(i,j,1)

          if (tw .gt. 0.0) then
            x(i,j) = x(i,j) / tw
          else
            x(i,j) = default_val
          endif

        enddo
      enddo
    endif

  end subroutine fill_nest_from_buffer_cell_center_masked


  subroutine fill_nest_from_buffer_cell_center_r4_3d(stagger, x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
    implicit none
    character ( len = 1 ), intent(in)             :: stagger
    real*4,  allocatable, intent(inout)           :: x(:,:,:)
    real*4,  allocatable, intent(in)              :: buffer(:,:,:)
    type(bbox), intent(in)                        :: bbox_fine, bbox_coarse
    integer, intent(in)                           :: nz
    integer, intent(in)                           :: dir, x_refine, y_refine
    real, allocatable, intent(in)                 :: wt(:,:,:)    ! The final dimension is always 4
    integer, allocatable, intent(in)              :: ind(:,:,:)

    character(len=8)       :: dir_str
    integer                :: i, j, k, ic, jc

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
      do k=1,nz
        do j=bbox_fine%js, bbox_fine%je
          do i=bbox_fine%is, bbox_fine%ie
            !if (stagger == "A") then
            !else if (stagger == "C") then
            !else if (stagger == "D") then
            !endif

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            x(i,j,k) = &
                wt(i,j,1)*buffer(ic,  jc,  k) +  &
                wt(i,j,2)*buffer(ic,  jc+1,k) +  &
                wt(i,j,3)*buffer(ic+1,jc+1,k) +  &
                wt(i,j,4)*buffer(ic+1,jc,  k)

          enddo
        enddo
      enddo
    endif

  end subroutine fill_nest_from_buffer_cell_center_r4_3d

  subroutine fill_nest_from_buffer_cell_center_r8_3d(stagger, x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
    implicit none
    character ( len = 1 ), intent(in)             :: stagger
    real*8, allocatable, intent(inout)            :: x(:,:,:)
    real*8, allocatable, intent(in)               :: buffer(:,:,:)
    type(bbox), intent(in)                        :: bbox_fine, bbox_coarse
    integer, intent(in)                           :: nz
    integer, intent(in)                           :: dir, x_refine, y_refine
    real, allocatable, intent(in)                 :: wt(:,:,:)    ! The final dimension is always 4
    integer, allocatable, intent(in)              :: ind(:,:,:)

    character(len=8)       :: dir_str
    integer                :: i, j, k, ic, jc

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
      do k=1,nz
        do j=bbox_fine%js, bbox_fine%je
          do i=bbox_fine%is, bbox_fine%ie
            !if (stagger == "A") then
            !else if (stagger == "C") then
            !else if (stagger == "D") then
            !endif

            ic = ind(i,j,1)
            jc = ind(i,j,2)

            x(i,j,k) = &
                wt(i,j,1)*buffer(ic,  jc,  k) +  &
                wt(i,j,2)*buffer(ic,  jc+1,k) +  &
                wt(i,j,3)*buffer(ic+1,jc+1,k) +  &
                wt(i,j,4)*buffer(ic+1,jc,  k)
          enddo
        enddo
      enddo
    endif

  end subroutine fill_nest_from_buffer_cell_center_r8_3d


  subroutine fill_nest_from_buffer_cell_center_r4_4d(stagger, x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
    implicit none
    character ( len = 1 ), intent(in)             :: stagger
    real*4, allocatable, intent(inout)            :: x(:,:,:,:)
    real*4, allocatable, intent(in)               :: buffer(:,:,:,:)
    type(bbox), intent(in)                        :: bbox_fine, bbox_coarse
    integer, intent(in)                           :: nz
    integer, intent(in)                           :: dir, x_refine, y_refine
    real, allocatable, intent(in)                 :: wt(:,:,:)    ! The final dimension is always 4
    integer, allocatable, intent(in)              :: ind(:,:,:)

    character(len=8)       :: dir_str
    integer                :: i, j, k, v, ic, jc

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
      do v=1,ubound(buffer,4)
        do k=1,nz
          do j=bbox_fine%js, bbox_fine%je
            do i=bbox_fine%is, bbox_fine%ie
              ic = ind(i,j,1)
              jc = ind(i,j,2)

              x(i,j,k,v) = &
                  wt(i,j,1)*buffer(ic,  jc,  k, v) +  &
                  wt(i,j,2)*buffer(ic,  jc+1,k, v) +  &
                  wt(i,j,3)*buffer(ic+1,jc+1,k, v) +  &
                  wt(i,j,4)*buffer(ic+1,jc,  k, v)
            enddo
          enddo
        enddo
      enddo
    endif

  end subroutine fill_nest_from_buffer_cell_center_r4_4d


  subroutine fill_nest_from_buffer_cell_center_r8_4d(stagger, x, buffer, bbox_fine, bbox_coarse, nz, dir, x_refine, y_refine, wt, ind)
    implicit none
    character ( len = 1 ), intent(in)             :: stagger
    real*8, allocatable, intent(inout)            :: x(:,:,:,:)
    real*8, allocatable, intent(in)               :: buffer(:,:,:,:)
    type(bbox), intent(in)                        :: bbox_fine, bbox_coarse
    integer, intent(in)                           :: nz
    integer, intent(in)                           :: dir, x_refine, y_refine
    real, allocatable, intent(in)                 :: wt(:,:,:)    ! The final dimension is always 4
    integer, allocatable, intent(in)              :: ind(:,:,:)

    character(len=8)       :: dir_str
    integer                :: i, j, k, v, ic, jc

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
      do v=1,ubound(buffer,4)
        do k=1,nz
          do j=bbox_fine%js, bbox_fine%je
            do i=bbox_fine%is, bbox_fine%ie
              ic = ind(i,j,1)
              jc = ind(i,j,2)

              x(i,j,k,v) = &
                  wt(i,j,1)*buffer(ic,  jc,  k, v) +  &
                  wt(i,j,2)*buffer(ic,  jc+1,k, v) +  &
                  wt(i,j,3)*buffer(ic+1,jc+1,k, v) +  &
                  wt(i,j,4)*buffer(ic+1,jc,  k, v)
            enddo
          enddo
        enddo
      enddo
    endif

  end subroutine fill_nest_from_buffer_cell_center_r8_4d


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
      do j=bbox_fine%js, bbox_fine%je
        do i=bbox_fine%is, bbox_fine%ie

          ic = bbox_coarse%is + 1
          jc = bbox_coarse%js + 1

          do k=1,nz

            ! Pick the maximum weight of the 4
            !   If two are tied for the max weight, use whichever one maxloc returns first
            !   TODO Might need a more deterministic algorithm here for reproducibility;  e.g. take the lowest index, etc.
            nearest_idx = maxloc(wt(i, j, :), 1)

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
              !if (debug_log) print '("[WARN] Nearest Neighbor algorithm mismatch index ",I0," buffer. npe=",I0," x(",I0,",",I0,",",I0,")=",F12.5)', nearest_idx, this_pe, i, j, k, x(i,j,k)
            end select
          enddo
        enddo
      enddo
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
        call mpp_error(FATAL, "fill_weight_grid invalid lower bounds")
      endif
      if (ubound(atm_wt, n) .ne. ubound(new_wt, n)) then
        call mpp_error(FATAL, "fill_weight_grid invalid upper bounds")
      endif
    enddo

    do x = lbound(atm_wt,1),ubound(atm_wt,1)
      do y = lbound(atm_wt,2),ubound(atm_wt,2)
        do z = 1,4
          atm_wt(x,y,z) = new_wt(x,y,z)
        enddo
      enddo
    enddo

  end subroutine fill_weight_grid


#endif ! MOVING_NEST

end module fv_moving_nest_utils_mod
