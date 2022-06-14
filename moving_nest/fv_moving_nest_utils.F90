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

  ! Added WDR
  use boundary_mod,      only: update_coarse_grid, update_coarse_grid_mpp
  use bounding_box_mod,  only: bbox, bbox_get_C2F_index, fill_bbox, show_bbox
  use fms_io_mod,        only: read_data, write_data, get_global_att_value, fms_io_init, fms_io_exit
  use fv_arrays_mod,     only: R_GRID
  use fv_arrays_mod,     only: fv_grid_type, fv_nest_type, fv_atmos_type
  use fv_surf_map_mod,   only: FV3_zs_filter
  use fv_moving_nest_types_mod, only: grid_geometry

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

  interface check_array
    module procedure check_array_r4_2d
    module procedure check_array_r4_3d
    module procedure check_array_r4_4d

    module procedure check_array_r8_2d
    module procedure check_array_r8_3d
    module procedure check_array_r8_4d
  end interface check_array

  interface check_local_array
    module procedure check_local_array_r4_2d
    module procedure check_local_array_r4_3d
    module procedure check_local_array_r8_2d
    module procedure check_local_array_r8_3d
  end interface check_local_array

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

        !if (this_pe .ge. 96) then
        !   print '("[INFO] WDR BLEND npe=",I0," a_step=",I0," i,j=",I0,",",I0," smoothed_orog=",F10.5," hires_orog=",F10.5," blend_wt=",F6.4," blend_orog=",F10.5)', this_pe, a_step, i, j, smoothed_orog, hires_orog, blend_wt, blend_orog
        !endif

      enddo
    enddo


    ! From tools/fv_surf_map.F90::surfdrv()
    !print '("[INFO] WDR BLEND npe=",I0," full_zs_filter=",L1," blend_size=",I0)', this_pe, Atm%flagstruct%full_zs_filter, blend_size
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

  ! Compare terrain for parent and nest cells - for debugging

  subroutine compare_terrain(var_name, data_var, interp_type, ind, x_refine, y_refine, is_fine_pe, nest_domain)
    character(len=*), intent(in)                :: var_name
    real, allocatable, intent(in)               :: data_var(:,:)
    integer, intent(in)                         :: interp_type
    integer, allocatable, intent(in)            :: ind(:,:,:)
    integer, intent(in)                         :: x_refine, y_refine
    logical, intent(in)                         :: is_fine_pe
    type(nest_domain_type), intent(inout)          :: nest_domain

    integer                             :: position = CENTER

    real, dimension(:,:), allocatable   :: nbuffer, sbuffer, ebuffer, wbuffer
    type(bbox)                          :: north_fine, north_coarse
    type(bbox)                          :: south_fine, south_coarse
    type(bbox)                          :: east_fine, east_coarse
    type(bbox)                          :: west_fine, west_coarse
    integer                             :: this_pe
    integer                             :: nest_level = 1  ! WDR TODO allow to vary

    this_pe = mpp_pe()

    !  Get the parent terrain through halo mechanism
    !print '("[INFO] WDR compare_terrain AA. npe=",I0)', this_pe
    call alloc_halo_buffer(nbuffer, north_fine, north_coarse, nest_domain, NORTH,  position)
    !print '("[INFO] WDR compare_terrain BB. npe=",I0)', this_pe
    call alloc_halo_buffer(sbuffer, south_fine, south_coarse, nest_domain, SOUTH,  position)
    !print '("[INFO] WDR compare_terrain CC. npe=",I0)', this_pe
    call alloc_halo_buffer(ebuffer, east_fine,  east_coarse,  nest_domain, EAST,   position)
    !print '("[INFO] WDR compare_terrain DD. npe=",I0)', this_pe
    call alloc_halo_buffer(wbuffer, west_fine,  west_coarse,  nest_domain, WEST,   position)
    !print '("[INFO] WDR compare_terrain EE. npe=",I0)', this_pe

    if (debug_log) print '("[INFO] WDR NRF0.d mn_var_shift_data npe=",I0," data_var(",I0,",",I0,")")', this_pe, size(data_var,1), size(data_var,2)

    ! Passes data from coarse grid to fine grid's halo
    call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

    !print '("[INFO] WDR compare_terrain FF. npe=",I0)', this_pe

    ! Figure out alignment of parent and child data and compare
    ! At most one of the buffers will have any data in it from the parent

    if (is_fine_pe) then
      call compare_buffer(north_coarse, north_fine, ind, nbuffer, data_var)
      !print '("[INFO] WDR compare_terrain GG. npe=",I0)', this_pe
      call compare_buffer(south_coarse, south_fine, ind, sbuffer, data_var)
      !print '("[INFO] WDR compare_terrain HH. npe=",I0)', this_pe
      call compare_buffer(east_coarse, east_fine, ind, ebuffer, data_var)
      !print '("[INFO] WDR compare_terrain II. npe=",I0)', this_pe
      call compare_buffer(west_coarse, west_fine, ind, wbuffer, data_var)
      !print '("[INFO] WDR compare_terrain JJ. npe=",I0)', this_pe
    endif

    print '("[INFO] WDR compare_terrain ZZ. npe=",I0)', this_pe

    deallocate(nbuffer)
    deallocate(sbuffer)
    deallocate(ebuffer)
    deallocate(wbuffer)

  end subroutine compare_terrain


  subroutine compare_buffer(bbox_coarse, bbox_fine, ind, buffer, fine_var)
    type(bbox), intent(in)                      :: bbox_coarse, bbox_fine
    integer, allocatable, intent(in)            :: ind(:,:,:)
    real, allocatable, intent(in)               :: buffer(:,:)
    real, allocatable, intent(in)               :: fine_var(:,:)


    integer :: i, j, ic, jc
    integer :: this_pe

    this_pe = mpp_pe()

    if ( bbox_coarse%ie .GE. bbox_coarse%is .AND. bbox_coarse%je .GE. bbox_coarse%js ) then
      !debug_log = .true.

      !if (debug_log) print '("[INFO] WDR BUFR print large buffer. npe=",I0," buffer(is_c, js_c)=",F12.5," buffer(ie_c-1, je_c-1)=",F12.5)', this_pe, buffer(bbox_coarse%is, bbox_coarse%js),  buffer(bbox_coarse%ie-1, bbox_coarse%je-1)

      if (debug_log) print '("[INFO] WDR BOUNDS i npe=",I0," is_f=",I0," ie_f=",I0," is_c=",I0," ie_c=",I0)', this_pe, bbox_fine%is, bbox_fine%ie, bbox_coarse%is, bbox_coarse%ie
      if (debug_log) print '("[INFO] WDR BOUNDS j npe=",I0," js_f=",I0," je_f=",I0," js_c=",I0," je_c=",I0)', this_pe, bbox_fine%js, bbox_fine%je, bbox_coarse%js, bbox_coarse%je

      if (debug_log) print '("[INFO] WDR BOUNDS fine_var npe=",I0," fine_var(",I0,"-",I0,",",I0,"-",I0,")")', this_pe, lbound(fine_var,1), ubound(fine_var,1), lbound(fine_var,2), ubound(fine_var,2)
      if (debug_log) print '("[INFO] WDR BOUNDS buffer npe=",I0," buffer(",I0,"-",I0,",",I0,"-",I0,")")', this_pe, lbound(buffer,1), ubound(buffer,1), lbound(buffer,2), ubound(buffer,2)

      do i=bbox_fine%is, bbox_fine%ie
        do j=bbox_fine%js, bbox_fine%je

          ic = ind(i,j,1)
          jc = ind(i,j,2)

          !print '("[INFO] WDR BOUNDS_ITER  npe=",I0," i=",I0," j=",I0," ic=",I0," jc=",I0)', this_pe, i, j, ic, jc
          !print '("[INFO] WDR BOUNDS_FINE npe=",I0," i=",I0," j=",I0," fine_var=",F12.5)', this_pe, i, j, fine_var(i,j)
          !print '("[INFO] WDR BOUNDS_BUFFER1 npe=",I0," ic=",I0," jc=",I0," buffer=",F12.5)', this_pe, ic, jc, buffer(ic,jc)
          !print '("[INFO] WDR BOUNDS_BUFFER2 npe=",I0," ic=",I0," jc=",I0," buffer=",F12.5)', this_pe, ic, jc+1, buffer(ic,jc+1)
          !print '("[INFO] WDR BOUNDS_BUFFER3 npe=",I0," ic=",I0," jc=",I0," buffer=",F12.5)', this_pe, ic+1, jc+1, buffer(ic+1,jc+1)
          !print '("[INFO] WDR BOUNDS_BUFFER4 npe=",I0," ic=",I0," jc=",I0," buffer=",F12.5)', this_pe, ic+1, jc, buffer(ic+1,jc)

          if ( (fine_var(i,j) .gt. 0.01) .or. &
              (buffer(ic,jc) .gt. 0.01) .or. &
              (buffer(ic,jc+1) .gt. 0.01) .or. &
              (buffer(ic+1,jc+1) .gt. 0.01) .or. &
              (buffer(ic+1,jc) .gt. 0.01)) then
            print '("[INFO] WDR COMP_TERR npe=",I0," i=",I0," j=",I0," ic=",I0," jc=",I0,F10.3," ",F10.3," ",F10.3," ",F10.3," ",F10.3)', this_pe, i, j, ic, jc, fine_var(i,j), buffer(ic,  jc  ), buffer(ic,  jc+1), buffer(ic+1,jc+1), buffer(ic+1,jc  )
          endif

          !wt(i,j,1)*buffer(ic,  jc  ) +  &
          !wt(i,j,2)*buffer(ic,  jc+1) +  &
          !wt(i,j,3)*buffer(ic+1,jc+1) +  &
          !wt(i,j,4)*buffer(ic+1,jc  )

        enddo
      enddo
      !print '("[INFO] WDR BOUNDS_DONE npe=",I0," i=",I0," j=",I0)', this_pe, i, j

      debug_log = .false.
      !else
      !   print '("[INFO] WDR NIL BUFR. npe=",I0)', this_pe
    endif
  end subroutine compare_buffer

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
    endif

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

    endif

    if (debug_log) print '("[INFO] WDR End fill_nest_halos_from_parent2D. npe=",I0," var_name=",A16)', this_pe, var_name

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
    integer                             :: nest_level = 1  ! WDR TODO allow to vary

    this_pe = mpp_pe()

    !!===========================================================
    !!
    !! Fill halo buffers
    !!
    !!===========================================================

    if (debug_log) then

      print '("[INFO] WDR Start fill_nest_halos_from_parent2D_kindphys. npe=",I0," var_name=",A16)', this_pe, var_name
      print '("[INFO] fill_nest_halos npe=",I0," nest_domain%tile_fine=",I0," %tile_coarse=",I0)', this_pe, nest_domain%tile_fine, nest_domain%tile_coarse
      print '("[INFO] fill_nest_halos npe=",I0," nest_domain%istart_fine=",I0," %iend_fine=",I0)', this_pe, nest_domain%istart_fine,  nest_domain%iend_fine
      print '("[INFO] fill_nest_halos npe=",I0," nest_domain%jstart_fine=",I0," %jend_fine=",I0)', this_pe, nest_domain%jstart_fine,  nest_domain%jend_fine
      print '("[INFO] fill_nest_halos npe=",I0," nest_domain%istart_coarse=",I0," %iend_coarse=",I0)', this_pe, nest_domain%istart_coarse,  nest_domain%iend_coarse
      print '("[INFO] fill_nest_halos npe=",I0," nest_domain%jstart_coarse=",I0," %jend_coarse=",I0)', this_pe, nest_domain%jstart_coarse,  nest_domain%jend_coarse
      print '("[INFO] data_var npe=",I0," var_name=",A16," data_var(",I0,"-",I0,",",I0,"-",I0,")")', &
          this_pe,  var_name, lbound(data_var, 1), ubound(data_var, 1), lbound(data_var, 2), ubound(data_var, 2)
      print '("[INFO] wt npe=",I0," var_name=",A16," wt(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', &
          this_pe,  var_name, lbound(wt, 1), ubound(wt, 1), lbound(wt, 2), ubound(wt, 2), lbound(wt, 3), ubound(wt, 3)
    endif

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

    endif

    deallocate(nbuffer)
    deallocate(sbuffer)
    deallocate(ebuffer)
    deallocate(wbuffer)

    if (debug_log) print '("[INFO] WDR End fill_nest_halos_from_parent2D_kindphys. npe=",I0," var_name=",A16)', this_pe, var_name

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
    integer                             :: nest_level = 1  ! WDR TODO allow to vary

    this_pe = mpp_pe()

    !!===========================================================
    !!
    !! Fill halo buffers
    !!
    !!===========================================================

    if (debug_log) then

      print '("[INFO] WDR Start fill_nest_halos_from_parent2D_kindphys. npe=",I0," var_name=",A16)', this_pe, var_name
      print '("[INFO] fill_nest_halos npe=",I0," nest_domain%tile_fine=",I0," %tile_coarse=",I0)', this_pe, nest_domain%tile_fine, nest_domain%tile_coarse
      print '("[INFO] fill_nest_halos npe=",I0," nest_domain%istart_fine=",I0," %iend_fine=",I0)', this_pe, nest_domain%istart_fine,  nest_domain%iend_fine
      print '("[INFO] fill_nest_halos npe=",I0," nest_domain%jstart_fine=",I0," %jend_fine=",I0)', this_pe, nest_domain%jstart_fine,  nest_domain%jend_fine
      print '("[INFO] fill_nest_halos npe=",I0," nest_domain%istart_coarse=",I0," %iend_coarse=",I0)', this_pe, nest_domain%istart_coarse,  nest_domain%iend_coarse
      print '("[INFO] fill_nest_halos npe=",I0," nest_domain%jstart_coarse=",I0," %jend_coarse=",I0)', this_pe, nest_domain%jstart_coarse,  nest_domain%jend_coarse
      print '("[INFO] data_var npe=",I0," var_name=",A16," data_var(",I0,"-",I0,",",I0,"-",I0,")")', &
          this_pe,  var_name, lbound(data_var, 1), ubound(data_var, 1), lbound(data_var, 2), ubound(data_var, 2)
      print '("[INFO] wt npe=",I0," var_name=",A16," wt(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', &
          this_pe,  var_name, lbound(wt, 1), ubound(wt, 1), lbound(wt, 2), ubound(wt, 2), lbound(wt, 3), ubound(wt, 3)
    endif

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

      call fill_nest_from_buffer_masked(interp_type, data_var, nbuffer, north_fine, north_coarse, NORTH, x_refine, y_refine, wt, ind, mask_var, mask_val, default_val)
      if (debug_log) print '("[INFO] WDR NRF N mn_var_shift_data start. npe=",I0)', this_pe

      call fill_nest_from_buffer_masked(interp_type, data_var, sbuffer, south_fine, south_coarse, SOUTH, x_refine, y_refine, wt, ind, mask_var, mask_val, default_val)
      if (debug_log) print '("[INFO] WDR NRF S mn_var_shift_data start. npe=",I0)', this_pe

      call fill_nest_from_buffer_masked(interp_type, data_var, ebuffer, east_fine, east_coarse, EAST, x_refine, y_refine, wt, ind, mask_var, mask_val, default_val)
      if (debug_log) print '("[INFO] WDR NRF E mn_var_shift_data start. npe=",I0)', this_pe

      call fill_nest_from_buffer_masked(interp_type, data_var, wbuffer, west_fine, west_coarse, WEST, x_refine, y_refine, wt, ind, mask_var, mask_val, default_val)
      if (debug_log) print '("[INFO] WDR NRF W mn_var_shift_data start. npe=",I0)', this_pe

    endif

    deallocate(nbuffer)
    deallocate(sbuffer)
    deallocate(ebuffer)
    deallocate(wbuffer)

    if (debug_log) print '("[INFO] WDR End fill_nest_halos_from_parent2D_kindphys. npe=",I0," var_name=",A16)', this_pe, var_name

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
    endif

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

    endif

    deallocate(nbuffer)
    deallocate(sbuffer)
    deallocate(ebuffer)
    deallocate(wbuffer)

    if (debug_log) print '("[INFO] WDR End fill_nest_halos_from_parent3D. npe=",I0," var_name=",A16)', this_pe, var_name

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
    endif

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

    endif

    deallocate(nbuffer)
    deallocate(sbuffer)
    deallocate(ebuffer)
    deallocate(wbuffer)

    if (debug_log) print '("[INFO] WDR End fill_nest_halos_from_parent3D. npe=",I0," var_name=",A16)', this_pe, var_name

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

    endif

    deallocate(nbuffer)
    deallocate(sbuffer)
    deallocate(ebuffer)
    deallocate(wbuffer)

    if (debug_log) print '("[INFO] WDR End fill_nest_halos_from_parent4D. npe=",I0," var_name=",A16)', this_pe, var_name

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

    endif

    deallocate(nbuffer)
    deallocate(sbuffer)
    deallocate(ebuffer)
    deallocate(wbuffer)

    if (debug_log) print '("[INFO] WDR End fill_nest_halos_from_parent4D_kindphys. npe=",I0," var_name=",A16)', this_pe, var_name

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

  end subroutine alloc_halo_buffer_r8_2d


  subroutine alloc_halo_buffer_r4_2d(buffer, bbox_fine, bbox_coarse, nest_domain, direction, position)
    real*4, dimension(:,:), allocatable, intent(out) :: buffer
    type(bbox), intent(out)                          :: bbox_fine, bbox_coarse
    type(nest_domain_type), intent(in)               :: nest_domain
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

  end subroutine alloc_halo_buffer_r4_2d


  subroutine alloc_halo_buffer_r4_3d(buffer, bbox_fine, bbox_coarse, nest_domain, direction, position, nz)
    real*4, dimension(:,:,:), allocatable, intent(out) :: buffer
    type(bbox), intent(out)                            :: bbox_fine, bbox_coarse
    type(nest_domain_type), intent(in)                 :: nest_domain
    integer, intent(in)                                :: direction, position, nz

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

  end subroutine alloc_halo_buffer_r4_3d


  subroutine alloc_halo_buffer_r8_3d(buffer, bbox_fine, bbox_coarse, nest_domain, direction, position, nz)
    real*8, dimension(:,:,:), allocatable, intent(out) :: buffer
    type(bbox), intent(out)                            :: bbox_fine, bbox_coarse
    type(nest_domain_type), intent(in)                 :: nest_domain
    integer, intent(in)                                :: direction, position, nz

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

  end subroutine alloc_halo_buffer_r8_3d


  subroutine alloc_halo_buffer_r4_4d(buffer, bbox_fine, bbox_coarse, nest_domain, direction, position, nz, n4d)
    real*4, dimension(:,:,:,:), allocatable, intent(out) :: buffer
    type(bbox), intent(out)                              :: bbox_fine, bbox_coarse
    type(nest_domain_type), intent(in)                   :: nest_domain
    integer, intent(in)                                  :: direction, position, nz, n4d

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

  end subroutine alloc_halo_buffer_r4_4d


  subroutine alloc_halo_buffer_r8_4d(buffer, bbox_fine, bbox_coarse, nest_domain, direction, position, nz, n4d)
    real*8, dimension(:,:,:,:), allocatable, intent(out) :: buffer
    type(bbox), intent(out)                              :: bbox_fine, bbox_coarse
    type(nest_domain_type), intent(in)                   :: nest_domain
    integer, intent(in)                                  :: direction, position, nz, n4d

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

  end subroutine alloc_halo_buffer_r8_4d


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

#ifdef OVERLOAD_R8
  subroutine alloc_read_data_r4_2d(nc_filename, var_name, x_size, y_size, data_array, time)
    character(len=*), intent(in)           :: nc_filename, var_name
    integer, intent(in)                    :: x_size, y_size
    real*4, allocatable, intent(inout)     :: data_array(:,:)
    integer, intent(in),optional           :: time

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

    if (present(time)) then
      start(3) = time
      nread(3) = 1
    endif

    if (debug_log) print '("[INFO] WDR NCREAD NCRA alloc_read_data. npe=",I0," ",A96," ", A16)', this_pe, trim(nc_filename), var_name
    if (debug_log) print '("[INFO] WDR NCREAD NCRB alloc_read_data, nread npe=",I0, " ", A16,I4,I4,I4,I4)', this_pe, var_name, start(1), start(2), nread(1), nread(2)

    call read_data(nc_filename, var_name, data_array, start, nread, no_domain=.TRUE.)  ! r4_2d

    if (debug_log) print '("[INFO] WDR NCREAD NCRC alloc_read_data, nread npe=",I0, " ", A16,I4,I4,I4,I4)', this_pe, var_name, start(1), start(2), nread(1), nread(2)

  end subroutine alloc_read_data_r4_2d
#endif

  subroutine alloc_read_data_r8_2d(nc_filename, var_name, x_size, y_size, data_array)
    character(len=*), intent(in)           :: nc_filename, var_name
    integer, intent(in)                    :: x_size, y_size
    real*8, allocatable, intent(inout)     :: data_array(:,:)

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

    call read_data(nc_filename, var_name, data_array, start, nread, no_domain=.TRUE.)  ! r8_2d

    if (debug_log) print '("[INFO] WDR NCREAD NCRC alloc_read_data, nread npe=",I0, " ", A16,I4,I4,I4,I4)', this_pe, var_name, start(1), start(2), nread(1), nread(2)

  end subroutine alloc_read_data_r8_2d


  !  nest_geo and parent_geo can be centered or supergrids.
  !  Assumes and validates that nest_geo is smaller, and inside parent_geo
  subroutine find_nest_alignment(nest_geo, parent_geo, nest_x, nest_y, parent_x, parent_y)
    implicit none
    type(grid_geometry), intent(in)     :: nest_geo, parent_geo
    integer, intent(out)                :: nest_x, nest_y, parent_x, parent_y

    type(bbox)  :: nest_bbox, parent_bbox
    integer     :: x,y
    logical     :: found

    real(kind=R_GRID) :: pi = 4 * atan(1.0d0)
    real                :: rad2deg
    integer     :: this_pe

    this_pe = mpp_pe()

    rad2deg =  180.0 / pi

    found = .false.
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
          endif

          if (   abs(abs(parent_geo%lons(x,y) - nest_geo%lons(nest_bbox.is, nest_bbox.js)) - 2*pi)  .lt. 0.0001) then
            found = .true.
            if (debug_log) print '("[INFO] WDR find_nest_alignment nest WRAP MATCH FOUND npe=",I0,F10.5, F10.5)', this_pe, nest_geo%lats(nest_bbox.is, nest_bbox.js)*rad2deg, nest_geo%lons(nest_bbox.is, nest_bbox.js)*rad2deg
            if (debug_log) print '("[INFO] WDR find_nest_alignment WRAP MATCH ",F10.5, F10.5)', parent_geo%lats(x,y)*rad2deg, parent_geo%lons(x,y)*rad2deg

            parent_x = x
            parent_y = y
            nest_x = nest_bbox.is
            nest_y = nest_bbox.js

            if (debug_log) print '("[INFO] WDR find_nest_alignment parent(",I0,",",I0,") nest(",I0,",",I0,")")', x,y,nest_bbox.is, nest_bbox.js
            if (debug_log) print '("[INFO] WDR find_nest_alignment ",F10.5, F10.5)', parent_geo%lats(x,y)*rad2deg, parent_geo%lons(x,y)*rad2deg
          endif
        endif
      enddo
    enddo

    if (found) then
      if (debug_log) print '("[INFO] WDR find_nest_alignment MATCH FOUND",F10.5, F10.5)', nest_geo%lats(nest_bbox.is, nest_bbox.js)*rad2deg, nest_geo%lons(nest_bbox.is, nest_bbox.js)*rad2deg
    endif

    if (.not. found .and. debug_log) then
      print '("[INFO] WDR find_nest_alignment nest NO MATCH FOUND npe=",I0,F10.5, F10.5)', this_pe, nest_geo%lats(nest_bbox.is, nest_bbox.js)*rad2deg, nest_geo%lons(nest_bbox.is, nest_bbox.js)*rad2deg
      print '("[INFO] WDR find_nest_alignment nest NO MATCH FOUND npe=",I0,F10.5, F10.5)', this_pe, nest_geo%lats(nest_bbox.is, nest_bbox.je)*rad2deg, nest_geo%lons(nest_bbox.is, nest_bbox.je)*rad2deg
      print '("[INFO] WDR find_nest_alignment nest NO MATCH FOUND npe=",I0,F10.5, F10.5)', this_pe, nest_geo%lats(nest_bbox.ie, nest_bbox.je)*rad2deg, nest_geo%lons(nest_bbox.ie, nest_bbox.je)*rad2deg
      print '("[INFO] WDR find_nest_alignment nest NO MATCH FOUND npe=",I0,F10.5, F10.5)', this_pe, nest_geo%lats(nest_bbox.ie, nest_bbox.js)*rad2deg, nest_geo%lons(nest_bbox.ie, nest_bbox.js)*rad2deg

      do x = parent_bbox.is, parent_bbox.ie
        do y = parent_bbox.js, parent_bbox.je          
          print '("[INFO] WDR find_nest_alignment parent NO MATCH FOUND npe="I0," ",I0," ",I0," ",F10.5, F10.5)', this_pe, x, y, parent_geo%lats(x,y)*rad2deg, parent_geo%lons(x,y)*rad2deg
        enddo
      enddo
    endif

  end subroutine find_nest_alignment


  !==================================================================================================
  !
  !  NetCDF Function Section
  !
  !==================================================================================================

  subroutine output_grid_to_nc_3d(flag, istart, iend, jstart, jend, k, grid, file_str, var_name, time_step, dom, position)
    implicit none

    character(len=*), intent(in)  :: flag
    integer, intent(in)           :: istart, iend, jstart, jend, k
    real, dimension(:,:,:), intent(in)   :: grid

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

    call write_data(filename, var_name, grid, dom, position=position)   ! r4_3d

  end subroutine output_grid_to_nc_3d


  subroutine output_grid_to_nc_2d(flag, istart, iend, jstart, jend, k, grid, file_str, var_name, time_step, dom, position)
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

    call write_data(filename, var_name, grid, dom, position=position)  ! r4_2d

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
        endif

        ! Make sure we don't run off the edge of the parent supergrid
        if (fp_i < fp_tile_bbox%is .or. fp_i > fp_tile_bbox%ie) then
          print '("[ERROR] WDR move_nest_geo invalid fp_i=",I0," is=",I0," ie=",I0)', fp_i, fp_tile_bbox%is, fp_tile_bbox%ie
          stop  ! TODO replace with a fatal error
        endif
        if (fp_j < fp_tile_bbox%js .or. fp_j > fp_tile_bbox%je) then
          print '("[ERROR] WDR move_nest_geo invalid fp_j=",I0," js=",I0," je=",I0)', fp_j, fp_tile_bbox%js, fp_tile_bbox%je
          stop  ! TODO replace with a fatal error
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
        endif

        ! Make sure we don't run off the edge of the parent supergrid
        if (fp_i < fp_tile_bbox%is .or. fp_i > fp_tile_bbox%ie) then
          print '("[ERROR] WDR move_nest_geo invalid fp_i=",I0," is=",I0," ie=",I0)', fp_i, fp_tile_bbox%is, fp_tile_bbox%ie
          stop  ! TODO replace with a fatal error
        endif
        if (fp_j < fp_tile_bbox%js .or. fp_j > fp_tile_bbox%je) then
          print '("[ERROR] WDR move_nest_geo invalid fp_j=",I0," js=",I0," je=",I0)', fp_j, fp_tile_bbox%js, fp_tile_bbox%je
          stop  ! TODO replace with a fatal error
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
        endif

        ! Make sure we don't run off the edge of the parent supergrid
        if (fp_i < fp_tile_bbox%is .or. fp_i > fp_tile_bbox%ie) then
          print '("[ERROR] WDR move_nest_geo invalid fp_i=",I0," is=",I0," ie=",I0)', fp_i, fp_tile_bbox%is, fp_tile_bbox%ie
          stop  ! TODO replace with a fatal error
        endif
        if (fp_j < fp_tile_bbox%js .or. fp_j > fp_tile_bbox%je) then
          print '("[ERROR] WDR move_nest_geo invalid fp_j=",I0," js=",I0," je=",I0)', fp_j, fp_tile_bbox%js, fp_tile_bbox%je
          stop  ! TODO replace with a fatal error
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
      if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= cell centered")', this_pe, interp_type
      call fill_nest_from_buffer_cell_center("A", x, buffer, bbox_fine, bbox_coarse, dir, x_refine, y_refine, wt, ind)
      !     case (3)
      !        if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= C grid staggered")', this_pe, interp_type
    case (4)
      if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= D grid staggered")', this_pe, interp_type
      call fill_nest_from_buffer_cell_center("D", x, buffer, bbox_fine, bbox_coarse, dir, x_refine, y_refine, wt, ind)
    case (7)
      if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= MASKED")', this_pe, interp_type
      call fill_nest_from_buffer_cell_center_masked("A", x, buffer, bbox_fine, bbox_coarse, dir, x_refine, y_refine, wt, ind, mask_var, mask_val, default_val)
    case (9)
      if (debug_log) print '("[INFO] WDR FNB this_tile. npe=",I0," interp_type=",I0,"= nearest neighbor cell centered")', this_pe, interp_type
      !call fill_nest_from_buffer_nearest_neighbor(x, buffer, bbox_fine, bbox_coarse, dir, wt)
      call mpp_error(FATAL, '2D fill_nest_from_buffer_nearest_neighbor not yet implemented.')
    case default
      if (debug_log) print '("[ERROR] WDR FNB this_tile. npe=",I0," UNDEFINED interp_type=",I0)', this_pe, interp_type
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
      call mpp_error(FATAL, 'fill_nest_from_buffer_nearest_neighbor is not yet implemented.')
    case default
      if (debug_log) print '("[ERROR] WDR FNB this_tile. npe=",I0," UNDEFINED interp_type=",I0)', this_pe, interp_type
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
          !endif

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
        enddo
      enddo
    else
      if (debug_log) print '("[INFO] WDR NIL BUFR print ",A8,"  buffer. npe=",I0)', dir_str, this_pe
    endif

    if (debug_log) print '("[INFO] WDR FILLNEST DONE print ",A8,"  buffer. npe=",I0)', dir_str, this_pe

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
          !endif

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
        enddo
      enddo
    else
      if (debug_log) print '("[INFO] WDR NIL BUFR print ",A8,"  buffer. npe=",I0)', dir_str, this_pe
    endif

    if (debug_log) print '("[INFO] WDR FILLNEST DONE print ",A8,"  buffer. npe=",I0)', dir_str, this_pe

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
    integer                :: focus_i = 1
    integer                :: focus_j = 1
    integer                :: this_pe
    real                   :: tw

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


          if (x(i,j) .lt. 0.0) print '("[WARN] WDR MASK npe=",I0," i,j=",I5,I5," x()=",F15.5," tw=",F10.5)', this_pe, i, j, x(i,j), tw

          !else
          !   x(i,j) = &
          !        wt(i,j,1)*buffer(ic,  jc  ) +  &
          !        wt(i,j,2)*buffer(ic,  jc+1) +  &
          !        wt(i,j,3)*buffer(ic+1,jc+1) +  &
          !        wt(i,j,4)*buffer(ic+1,jc  )
          !endif

          !call check_array(buffer, this_pe, "buffer"//dir_str, -300.0, 300.0)
          !call check_array(wt, this_pe, "wt"//dir_str, 0.0, 1.0)
          if (debug_log) print '("[INFO] WDR FILL WEIGHTS ",A8,"  npe=",I0," (",I0,",",I0,") ic,jc=(",I0,",",I0,"): wt:",F12.5,F12.5,F12.5,F12.5)', dir_str, this_pe, i, j, ic, jc, wt(i,j,1), wt(i,j,2), wt(i,j,3), wt(i,j,4)
        enddo
      enddo
    else
      if (debug_log) print '("[INFO] WDR NIL BUFR print ",A8,"  buffer. npe=",I0)', dir_str, this_pe
    endif

    if (debug_log) print '("[INFO] WDR FILLNEST DONE print ",A8,"  buffer. npe=",I0)', dir_str, this_pe

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
            !endif

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
            !endif

          enddo
        enddo
      enddo
    else
      if (debug_log) print '("[INFO] WDR NIL BUFR print ",A8,"  buffer. npe=",I0)', dir_str, this_pe
      !if (debug_log) print '("[INFO WDR NIL BUFR ",A8," BOUNDS i npe=",I0,"is_f=",I0," ie_f=",I0,"is_c=",I0," ie_c=",I0)', dir_str, this_pe, is_f, ie_f, is_c, ie_c
      !if (debug_log) print '("[INFO WDR NIL BUFR ",A8," BOUNDS j npe=",I0,"js_f=",I0," je_f=",I0,"js_c=",I0," je_c=",I0)', dir_str, this_pe, js_f, je_f, js_c, je_c

    endif

    if (debug_log) print '("[INFO] WDR FILLNEST DONE print ",A8,"  buffer. npe=",I0)', dir_str, this_pe

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
            !endif

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
            !endif

          enddo
        enddo
      enddo
    else
      if (debug_log) print '("[INFO] WDR NIL BUFR print ",A8,"  buffer. npe=",I0)', dir_str, this_pe
      !if (debug_log) print '("[INFO WDR NIL BUFR ",A8," BOUNDS i npe=",I0,"is_f=",I0," ie_f=",I0,"is_c=",I0," ie_c=",I0)', dir_str, this_pe, is_f, ie_f, is_c, ie_c
      !if (debug_log) print '("[INFO WDR NIL BUFR ",A8," BOUNDS j npe=",I0,"js_f=",I0," je_f=",I0,"js_c=",I0," je_c=",I0)', dir_str, this_pe, js_f, je_f, js_c, je_c

    endif

    if (debug_log) print '("[INFO] WDR FILLNEST DONE print ",A8,"  buffer. npe=",I0)', dir_str, this_pe

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

            enddo
          enddo
        enddo
      enddo
    else
      if (debug_log) print '("[INFO] WDR NIL BUFR print ",A8,"  buffer. npe=",I0)', dir_str, this_pe
    endif

    if (debug_log) print '("[INFO] WDR FILLNEST4D DONE print ",A8,"  buffer. npe=",I0)', dir_str, this_pe

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

            enddo
          enddo
        enddo
      enddo
    else
      if (debug_log) print '("[INFO] WDR NIL BUFR print ",A8,"  buffer. npe=",I0)', dir_str, this_pe
    endif

    if (debug_log) print '("[INFO] WDR FILLNEST4D DONE print ",A8,"  buffer. npe=",I0)', dir_str, this_pe

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
          enddo
        enddo
      enddo
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
      endif
      if (ubound(atm_wt, n) .ne. ubound(new_wt, n)) then
        print '("[ERROR] WDR fill_weight_grid ubound mismatch fv_moving_nest.F90 npe=",I0," n=",I0, I0, I0)', this_pe, n, ubound(atm_wt, n), ubound(new_wt, n)
        stop
      endif
    enddo

    if (debug_log) print '("[INFO] WDR running fill_weight_grid fv_moving_nest.F90 npe=",I0)', this_pe
    do x = lbound(atm_wt,1),ubound(atm_wt,1)
      do y = lbound(atm_wt,2),ubound(atm_wt,2)
        do z = 1,4
          atm_wt(x,y,z) = new_wt(x,y,z)
        enddo
      enddo
    enddo

  end subroutine fill_weight_grid


  !==================================================================================================
  !
  !  Array Checking Section
  !
  !==================================================================================================

  subroutine check_array_r4_2d(array, this_pe, var_name, min_range, max_range)
    real*4, intent(in), allocatable        :: array(:,:)
    integer, intent(in)                    :: this_pe
    character(len=*), intent(in)           :: var_name
    real, intent(in)                       :: min_range, max_range

    integer  :: i,j
    integer  :: num_invalid
    integer  :: num_valid
    real     :: eps = 0.0001
    real     :: invalid_last

    invalid_last = 0.0

    if (allocated(array)) then

      print '("[INFO] WDR 2Darray allocated  npe=",I0," ",A32,"(",I4,":",I4,",",I4,":",I4,")")', this_pe, var_name, lbound(array,1), ubound(array,1), lbound(array,2), ubound(array,2)

      num_invalid = 0
      num_valid = 0

      do i = lbound(array,1), ubound(array,1)
        do j =lbound(array,2), ubound(array,2)
          if (array(i,j) < min_range - eps) then
            num_invalid = num_invalid + 1
            invalid_last = array(i,j)
          elseif (array(i,j) > max_range + eps) then
            num_invalid = num_invalid + 1
            invalid_last = array(i,j)
          else
            num_valid = num_valid + 1
          endif
        enddo
      enddo

      if (num_invalid > 0 ) then
        print '("[ERROR] WDR 2Darray invalid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0," last invalid=",E12.5)', this_pe, var_name, num_invalid, num_valid, invalid_last
      else
        print '("[INFO] WDR 2Darray all valid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
      endif

    else
      print '("[INFO] WDR 2Darray not allocated  npe=",I0," ",A32)', this_pe, var_name
    endif

  end subroutine check_array_r4_2d


  subroutine check_array_r8_2d(array, this_pe, var_name, min_range, max_range)
    real*8, intent(in), allocatable        :: array(:,:)
    integer, intent(in)                    :: this_pe
    character(len=*), intent(in)           :: var_name
    real(kind=R_GRID), intent(in)          :: min_range, max_range

    integer  :: i,j
    integer  :: num_invalid
    integer  :: num_valid
    real     :: eps = 0.0001
    real     :: invalid_last

    invalid_last = 0.0

    if (allocated(array)) then

      print '("[INFO] WDR 2D64array allocated  npe=",I0," ",A32,"(",I4,":",I4,",",I4,":",I4,")")', this_pe, var_name, lbound(array,1), ubound(array,1), lbound(array,2), ubound(array,2)

      num_invalid = 0
      num_valid = 0

      do i = lbound(array,1), ubound(array,1)
        do j =lbound(array,2), ubound(array,2)
          if (array(i,j) < min_range - eps) then
            num_invalid = num_invalid + 1
            invalid_last = array(i,j)
          elseif (array(i,j) > max_range + eps) then
            num_invalid = num_invalid + 1
            invalid_last = array(i,j)
          else
            num_valid = num_valid + 1
          endif
        enddo
      enddo

      if (num_invalid > 0 ) then
        print '("[ERROR] WDR 2D64array invalid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0," last invalid=",E12.5)', this_pe, var_name, num_invalid, num_valid, invalid_last
      else
        print '("[INFO] WDR 2D64array all valid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
      endif

    else
      print '("[INFO] WDR 2D64array not allocated  npe=",I0," ",A32)', this_pe, var_name
    endif

  end subroutine check_array_r8_2d


  subroutine check_local_array_r4_2d(array, this_pe, var_name, min_range, max_range)
    real*4, intent(in)                     :: array(:,:)
    integer, intent(in)                    :: this_pe
    character(len=*), intent(in)           :: var_name
    real, intent(in)                       :: min_range, max_range

    integer  :: i,j
    integer  :: num_invalid
    integer  :: num_valid
    real     :: eps = 0.0001

    print '("[INFO] WDR 2DLarray allocated  npe=",I0," ",A32,"(",I4,":",I4,",",I4,":",I4,")")', this_pe, var_name, lbound(array,1), ubound(array,1), lbound(array,2), ubound(array,2)

    num_invalid = 0
    num_valid = 0

    do i = lbound(array,1), ubound(array,1)
      do j =lbound(array,2), ubound(array,2)
        if (array(i,j) < min_range - eps) then
          num_invalid = num_invalid + 1
        elseif (array(i,j) > max_range + eps) then
          num_invalid = num_invalid + 1
        else
          num_valid = num_valid + 1
        endif
      enddo
    enddo

    if (num_invalid > 0 ) then
      print '("[ERROR] WDR 2DLarray invalid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
    else
      print '("[INFO] WDR 2DLarray all valid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
    endif

  end subroutine check_local_array_r4_2d

  subroutine check_local_array_r8_2d(array, this_pe, var_name, min_range, max_range)
    real*8, intent(in)                     :: array(:,:)
    integer, intent(in)                    :: this_pe
    character(len=*), intent(in)           :: var_name
    real, intent(in)                       :: min_range, max_range

    integer  :: i,j
    integer  :: num_invalid
    integer  :: num_valid
    real     :: eps = 0.0001

    print '("[INFO] WDR 2DLarray allocated  npe=",I0," ",A32,"(",I4,":",I4,",",I4,":",I4,")")', this_pe, var_name, lbound(array,1), ubound(array,1), lbound(array,2), ubound(array,2)

    num_invalid = 0
    num_valid = 0

    do i = lbound(array,1), ubound(array,1)
      do j =lbound(array,2), ubound(array,2)
        if (array(i,j) < min_range - eps) then
          num_invalid = num_invalid + 1
        elseif (array(i,j) > max_range + eps) then
          num_invalid = num_invalid + 1
        else
          num_valid = num_valid + 1
        endif
      enddo
    enddo

    if (num_invalid > 0 ) then
      print '("[ERROR] WDR 2DLarray invalid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
    else
      print '("[INFO] WDR 2DLarray all valid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
    endif

  end subroutine check_local_array_r8_2d

  subroutine check_array_r4_3d(array, this_pe, var_name, min_range, max_range)
    real*4, intent(in), allocatable        :: array(:,:,:)
    integer, intent(in)                    :: this_pe
    character(len=*), intent(in)           :: var_name
    real, intent(in)                       :: min_range, max_range

    integer  :: i,j,k
    integer  :: num_invalid
    integer  :: num_valid
    real     :: eps = 0.0001

    if (allocated(array)) then

      print '("[INFO] WDR 3Darray allocated  npe=",I0," ",A32,"(",I4,":",I4,",",I4,":",I4,",",I4,":",I4,")")', this_pe, var_name, lbound(array,1), ubound(array,1), lbound(array,2), ubound(array,2), lbound(array,3), ubound(array,3)

      num_invalid = 0
      num_valid = 0

      do i = lbound(array,1), ubound(array,1)
        do j =lbound(array,2), ubound(array,2)
          do k =lbound(array,3), ubound(array,3)
            if (isnan(array(i,j,k))) then
              num_invalid = num_invalid + 1
            elseif (array(i,j,k) < min_range - eps) then
              num_invalid = num_invalid + 1
            elseif (array(i,j,k) > max_range + eps) then
              num_invalid = num_invalid + 1
            else
              num_valid = num_valid + 1
            endif
          enddo
        enddo
      enddo

      if (num_invalid > 0 ) then
        print '("[ERROR] WDR 3Darray invalid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
      else
        print '("[INFO] WDR 3Darray all valid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
      endif

    else
      print '("[INFO] WDR 3Darray not allocated  npe=",I0," ",A32)', this_pe, var_name
    endif

  end subroutine check_array_r4_3d

  subroutine check_array_r8_3d(array, this_pe, var_name, min_range, max_range)
    real*8, intent(in), allocatable        :: array(:,:,:)
    integer, intent(in)                    :: this_pe
    character(len=*), intent(in)           :: var_name
    real, intent(in)                       :: min_range, max_range

    integer  :: i,j,k
    integer  :: num_invalid
    integer  :: num_valid
    real     :: eps = 0.0001

    if (allocated(array)) then

      print '("[INFO] WDR 3Darray allocated  npe=",I0," ",A32,"(",I4,":",I4,",",I4,":",I4,",",I4,":",I4,")")', this_pe, var_name, lbound(array,1), ubound(array,1), lbound(array,2), ubound(array,2), lbound(array,3), ubound(array,3)

      num_invalid = 0
      num_valid = 0

      do i = lbound(array,1), ubound(array,1)
        do j =lbound(array,2), ubound(array,2)
          do k =lbound(array,3), ubound(array,3)
            if (isnan(array(i,j,k))) then
              num_invalid = num_invalid + 1
            elseif (array(i,j,k) < min_range - eps) then
              num_invalid = num_invalid + 1
            elseif (array(i,j,k) > max_range + eps) then
              num_invalid = num_invalid + 1
            else
              num_valid = num_valid + 1
            endif
          enddo
        enddo
      enddo

      if (num_invalid > 0 ) then
        print '("[ERROR] WDR 3Darray invalid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
      else
        print '("[INFO] WDR 3Darray all valid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
      endif

    else
      print '("[INFO] WDR 3Darray not allocated  npe=",I0," ",A32)', this_pe, var_name
    endif

  end subroutine check_array_r8_3d

  subroutine check_local_array_r4_3d(array, this_pe, var_name, min_range, max_range)
    real*4, intent(in)                     :: array(:,:,:)
    integer, intent(in)                    :: this_pe
    character(len=*), intent(in)           :: var_name
    real, intent(in)                       :: min_range, max_range

    integer  :: i,j,k
    integer  :: num_invalid
    integer  :: num_valid
    real     :: eps = 0.0001

    print '("[INFO] WDR 3DLarray bounds npe=",I0," ",A32,"(",I4,":",I4,",",I4,":",I4,",",I4,":",I4,")")', this_pe, var_name, lbound(array,1), ubound(array,1), lbound(array,2), ubound(array,2), lbound(array,3), ubound(array,3)

    num_invalid = 0
    num_valid = 0

    do i = lbound(array,1), ubound(array,1)
      do j =lbound(array,2), ubound(array,2)
        do k =lbound(array,3), ubound(array,3)
          if (isnan(array(i,j,k))) then
            num_invalid = num_invalid + 1
          elseif (array(i,j,k) < min_range - eps) then
            num_invalid = num_invalid + 1
          elseif (array(i,j,k) > max_range + eps) then
            num_invalid = num_invalid + 1
          else
            num_valid = num_valid + 1
          endif
        enddo
      enddo
    enddo

    if (num_invalid > 0 ) then
      print '("[ERROR] WDR 3DLarray invalid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
    else
      print '("[INFO] WDR 3DLarray all valid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
    endif

  end subroutine check_local_array_r4_3d

  subroutine check_local_array_r8_3d(array, this_pe, var_name, min_range, max_range)
    real*8, intent(in)                     :: array(:,:,:)
    integer, intent(in)                    :: this_pe
    character(len=*), intent(in)           :: var_name
    real, intent(in)                       :: min_range, max_range

    integer  :: i,j,k
    integer  :: num_invalid
    integer  :: num_valid
    real     :: eps = 0.0001

    print '("[INFO] WDR 3DLarray bounds npe=",I0," ",A32,"(",I4,":",I4,",",I4,":",I4,",",I4,":",I4,")")', this_pe, var_name, lbound(array,1), ubound(array,1), lbound(array,2), ubound(array,2), lbound(array,3), ubound(array,3)

    num_invalid = 0
    num_valid = 0

    do i = lbound(array,1), ubound(array,1)
      do j =lbound(array,2), ubound(array,2)
        do k =lbound(array,3), ubound(array,3)
          if (isnan(array(i,j,k))) then
            num_invalid = num_invalid + 1
          elseif (array(i,j,k) < min_range - eps) then
            num_invalid = num_invalid + 1
          elseif (array(i,j,k) > max_range + eps) then
            num_invalid = num_invalid + 1
          else
            num_valid = num_valid + 1
          endif
        enddo
      enddo
    enddo

    if (num_invalid > 0 ) then
      print '("[ERROR] WDR 3DLarray invalid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
    else
      print '("[INFO] WDR 3DLarray all valid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
    endif

  end subroutine check_local_array_r8_3d


  subroutine check_array_r4_4d(array, this_pe, var_name, min_range, max_range)
    real*4, intent(in), allocatable        :: array(:,:,:,:)
    integer, intent(in)                    :: this_pe
    character(len=*), intent(in)           :: var_name
    real, intent(in)                       :: min_range, max_range

    integer  :: i,j,k,v
    integer  :: num_invalid
    integer  :: num_valid
    real     :: eps = 0.0001

    if (allocated(array)) then

      print '("[INFO] WDR 4Darray allocated  npe=",I0," ",A32,"(",I4,":",I4,",",I4,":",I4,",",I4,":",I4,",",I4,":",I4,")")', this_pe, var_name, lbound(array,1), ubound(array,1), lbound(array,2), ubound(array,2), lbound(array,3), ubound(array,3), lbound(array,4), ubound(array,4)

      num_invalid = 0
      num_valid = 0

      do i = lbound(array,1), ubound(array,1)
        do j =lbound(array,2), ubound(array,2)
          do k =lbound(array,3), ubound(array,3)
            do v =lbound(array,4), ubound(array,4)
              if (isnan(array(i,j,k,v))) then
                num_invalid = num_invalid + 1
              elseif (array(i,j,k,v) < min_range - eps) then
                num_invalid = num_invalid + 1
              elseif (array(i,j,k,v) > max_range + eps) then
                num_invalid = num_invalid + 1
              else
                num_valid = num_valid + 1
              endif
            enddo
          enddo
        enddo
      enddo

      if (num_invalid > 0 ) then
        print '("[ERROR] WDR 4Darray invalid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
      else
        print '("[INFO] WDR 4Darray all valid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
      endif
    else
      print '("[INFO] WDR 4Darray not allocated  npe=",I0," ",A32)', this_pe, var_name
    endif
  end subroutine check_array_r4_4d


  subroutine check_array_r8_4d(array, this_pe, var_name, min_range, max_range)
    real*8, intent(in), allocatable        :: array(:,:,:,:)
    integer, intent(in)                    :: this_pe
    character(len=*), intent(in)           :: var_name
    real, intent(in)                       :: min_range, max_range

    integer  :: i,j,k,v
    integer  :: num_invalid
    integer  :: num_valid
    real     :: eps = 0.0001

    if (allocated(array)) then

      print '("[INFO] WDR 4Darray allocated  npe=",I0," ",A32,"(",I4,":",I4,",",I4,":",I4,",",I4,":",I4,",",I4,":",I4,")")', this_pe, var_name, lbound(array,1), ubound(array,1), lbound(array,2), ubound(array,2), lbound(array,3), ubound(array,3), lbound(array,4), ubound(array,4)

      num_invalid = 0
      num_valid = 0

      do i = lbound(array,1), ubound(array,1)
        do j =lbound(array,2), ubound(array,2)
          do k =lbound(array,3), ubound(array,3)
            do v =lbound(array,4), ubound(array,4)
              if (isnan(array(i,j,k,v))) then
                num_invalid = num_invalid + 1
              elseif (array(i,j,k,v) < min_range - eps) then
                num_invalid = num_invalid + 1
              elseif (array(i,j,k,v) > max_range + eps) then
                num_invalid = num_invalid + 1
              else
                num_valid = num_valid + 1
              endif
            enddo
          enddo
        enddo
      enddo

      if (num_invalid > 0 ) then
        print '("[ERROR] WDR 4Darray invalid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
      else
        print '("[INFO] WDR 4Darray all valid entries npe=",I0," ",A32," num_invalid=",I0," num_valid=",I0)', this_pe, var_name, num_invalid, num_valid
      endif
    else
      print '("[INFO] WDR 4Darray not allocated  npe=",I0," ",A32)', this_pe, var_name
    endif
  end subroutine check_array_r8_4d


  !==================================================================================================
  !
  !  Debugging Function Section
  !
  !==================================================================================================

  subroutine grid_equal(grid1, grid2, tag, this_pe, is_equal)
    real, allocatable, intent(in)    :: grid1(:,:,:)
    real, allocatable, intent(in)    :: grid2(:,:,:)
    character(len=*), intent(in)     :: tag
    integer, intent(in)              :: this_pe
    logical, intent(out)             :: is_equal

    integer :: x,y,z

    real                 :: pi = 4 * atan(1.0d0)
    real                 :: rad2deg

    rad2deg = 180.0 / pi

    is_equal = .true.

    do x=1,3
      if (lbound(grid1,x) /= lbound(grid2,x)) then
        print '("[ERROR] WDR grid_equal ",A16," npe=",I0," lbound mismatch ",I0, I0,I0)', tag, x, lbound(grid1,x), lbound(grid2,x)
        is_equal = .false.
      endif
      if (ubound(grid1,x) /= ubound(grid2,x)) then
        print '("[ERROR] WDR grid_equal ",A16," npe=",I0," ubound mismatch ",I0, I0,I0)', tag, x, ubound(grid1,x), ubound(grid2,x)
        is_equal = .false.
      endif
    enddo

    if (is_equal) then
      do x=lbound(grid1,1), ubound(grid1,1)
        do y=lbound(grid1,2), ubound(grid1,2)
          do z=lbound(grid1,3), ubound(grid1,3)
            if ( abs(grid1(x,y,z) - grid2(x,y,z)) > 0.0001 ) then
              print '("[ERROR] WDR grid_equal ",A16," npe=",I0," DEG value mismatch at (",I0,",",I0,",",I0,") ",F15.11, " ",F15.11, " ",F15.11)', tag, this_pe, x, y, z, grid1(x,y,z)*rad2deg, grid2(x,y,z)*rad2deg, grid1(x,y,z)*rad2deg - grid2(x,y,z)*rad2deg

              print '("[ERROR] WDR grid_equal ",A16," npe=",I0," RAD value mismatch at (",I0,",",I0,",",I0,") ",F15.11, " ",F15.11, " ",F15.11)', tag, this_pe, x, y, z, grid1(x,y,z), grid2(x,y,z), grid1(x,y,z) - grid2(x,y,z)
              is_equal = .false.
            else
              print '("[INFO]  WDR grid_equal ",A16," npe=",I0," DEG value match    at (",I0,",",I0,",",I0,") ",F15.11, " ",F15.11, " ",F15.11)', tag, this_pe, x, y, z, grid1(x,y,z)*rad2deg, grid2(x,y,z)*rad2deg, grid1(x,y,z)*rad2deg - grid2(x,y,z)*rad2deg

              print '("[INFO]  WDR grid_equal ",A16," npe=",I0," RAD value match    at (",I0,",",I0,",",I0,") ",F15.11, " ",F15.11, " ",F15.11)', tag, this_pe, x, y, z, grid1(x,y,z), grid2(x,y,z), grid1(x,y,z) - grid2(x,y,z)
            endif
          enddo
        enddo
      enddo
    endif

    if (is_equal) then
      print '("[INFO] WDR grid_equal ",A16," npe=",I0," MATCH.")', tag, this_pe
    else
      print '("[ERROR] WDR grid_equal ",A16," npe=",I0," MISMATCH.")', tag, this_pe
    endif

  end subroutine grid_equal


  subroutine show_atm_grids(Atm, n)
    type(fv_atmos_type), allocatable, intent(in) :: Atm(:)
    integer, intent(in)                          :: n

    integer :: x,y
    real(kind=R_GRID) :: pi = 4 * atan(1.0d0)
    real                :: pi180
    real                :: rad2deg, deg2rad

    pi180 = pi / 180.0
    deg2rad = pi / 180.0
    rad2deg = 1.0 / pi180

    print *, "[INFO] WDR MV_NST2 shape(Atm(1)%grid_global)=", shape(Atm(1)%grid_global)
    print '("[INFO] WDR MV_NST2 bounds1 (Atm(1)%grid_global)=",I0,"-",I0)', lbound(Atm(1)%grid_global,1), ubound(Atm(1)%grid_global,1)
    print '("[INFO] WDR MV_NST2 bounds2 (Atm(1)%grid_global)=",I0,"-",I0)', lbound(Atm(1)%grid_global,2), ubound(Atm(1)%grid_global,2)
    print '("[INFO] WDR MV_NST2 bounds3 (Atm(1)%grid_global)=",I0,"-",I0)', lbound(Atm(1)%grid_global,3), ubound(Atm(1)%grid_global,3)
    print '("[INFO] WDR MV_NST2 bounds4 (Atm(1)%grid_global)=",I0,"-",I0)', lbound(Atm(1)%grid_global,4), ubound(Atm(1)%grid_global,4)

    print *, "[INFO] WDR MV_NST2 shape(Atm(n)%grid_global)=", shape(Atm(n)%grid_global)
    print '("[INFO] WDR MV_NST2 bounds1 (Atm(n)%grid_global)=",I0,"-",I0)', lbound(Atm(n)%grid_global,1), ubound(Atm(n)%grid_global,1)
    print '("[INFO] WDR MV_NST2 bounds2 (Atm(n)%grid_global)=",I0,"-",I0)', lbound(Atm(n)%grid_global,2), ubound(Atm(n)%grid_global,2)
    print '("[INFO] WDR MV_NST2 bounds3 (Atm(n)%grid_global)=",I0,"-",I0)', lbound(Atm(n)%grid_global,3), ubound(Atm(n)%grid_global,3)
    print '("[INFO] WDR MV_NST2 bounds4 (Atm(n)%grid_global)=",I0,"-",I0)', lbound(Atm(n)%grid_global,4), ubound(Atm(n)%grid_global,4)

    print *, "[INFO] WDR MV_NST2 shape(Atm(n)%gridstruct%grid)=", shape(Atm(n)%gridstruct%grid)
    print '("[INFO] WDR MV_NST2 bounds1 (Atm(n)%gridstruct%grid)=",I0,"-",I0)', lbound(Atm(n)%gridstruct%grid,1), ubound(Atm(n)%gridstruct%grid,1)
    print '("[INFO] WDR MV_NST2 bounds2 (Atm(n)%gridstruct%grid)=",I0,"-",I0)', lbound(Atm(n)%gridstruct%grid,2), ubound(Atm(n)%gridstruct%grid,2)
    print '("[INFO] WDR MV_NST2 bounds3 (Atm(n)%gridstruct%grid)=",I0,"-",I0)', lbound(Atm(n)%gridstruct%grid,3), ubound(Atm(n)%gridstruct%grid,3)

    print *, "[INFO] WDR MV_NST2 shape(Atm(n)%gridstruct%agrid)=", shape(Atm(n)%gridstruct%agrid)
    print '("[INFO] WDR MV_NST2 bounds1 (Atm(n)%gridstruct%agrid)=",I0,"-",I0)', lbound(Atm(n)%gridstruct%agrid,1), ubound(Atm(n)%gridstruct%agrid,1)
    print '("[INFO] WDR MV_NST2 bounds2 (Atm(n)%gridstruct%agrid)=",I0,"-",I0)', lbound(Atm(n)%gridstruct%agrid,2), ubound(Atm(n)%gridstruct%agrid,2)
    print '("[INFO] WDR MV_NST2 bounds3 (Atm(n)%gridstruct%agrid)=",I0,"-",I0)', lbound(Atm(n)%gridstruct%agrid,3), ubound(Atm(n)%gridstruct%agrid,3)

    x = lbound(Atm(n)%gridstruct%agrid,1)
    y = lbound(Atm(n)%gridstruct%agrid,2)
    print '("[INFO] WDR GRD_SHOa atmosphere.F90 Atm(n)%agrid(",I0,",",I0,")=",F10.5, F10.5)', x, y, Atm(n)%gridstruct%agrid(x,y,2)*rad2deg, Atm(n)%gridstruct%agrid(x,y,1)*rad2deg

    x = ubound(Atm(n)%gridstruct%agrid,1)
    y = ubound(Atm(n)%gridstruct%agrid,2)
    print '("[INFO] WDR GRD_SHOb atmosphere.F90 Atm(n)%agrid(",I0,",",I0,")=",F10.5, F10.5)', x, y, Atm(n)%gridstruct%agrid(x,y,2)*rad2deg, Atm(n)%gridstruct%agrid(x,y,1)*rad2deg

    print *, "[INFO] WDR MV_NST2 shape(Atm(n)%gridstruct%grid_64)=", shape(Atm(n)%gridstruct%grid_64)
    print '("[INFO] WDR MV_NST2 bounds1 (Atm(n)%gridstruct%grid_64)=",I0,"-",I0)', lbound(Atm(n)%gridstruct%grid_64,1), ubound(Atm(n)%gridstruct%grid_64,1)
    print '("[INFO] WDR MV_NST2 bounds2 (Atm(n)%gridstruct%grid_64)=",I0,"-",I0)', lbound(Atm(n)%gridstruct%grid_64,2), ubound(Atm(n)%gridstruct%grid_64,2)
    print '("[INFO] WDR MV_NST2 bounds3 (Atm(n)%gridstruct%grid_64)=",I0,"-",I0)', lbound(Atm(n)%gridstruct%grid_64,3), ubound(Atm(n)%gridstruct%grid_64,3)

    print *, "[INFO] WDR MV_NST2 shape(Atm(n)%gridstruct%agrid_64)=", shape(Atm(n)%gridstruct%agrid_64)
    print '("[INFO] WDR MV_NST2 bounds1 (Atm(n)%gridstruct%agrid_64)=",I0,"-",I0)', lbound(Atm(n)%gridstruct%agrid_64,1), ubound(Atm(n)%gridstruct%agrid_64,1)
    print '("[INFO] WDR MV_NST2 bounds2 (Atm(n)%gridstruct%agrid_64)=",I0,"-",I0)', lbound(Atm(n)%gridstruct%agrid_64,2), ubound(Atm(n)%gridstruct%agrid_64,2)
    print '("[INFO] WDR MV_NST2 bounds3 (Atm(n)%gridstruct%agrid_64)=",I0,"-",I0)', lbound(Atm(n)%gridstruct%agrid_64,3), ubound(Atm(n)%gridstruct%agrid_64,3)

  end subroutine show_atm_grids


  subroutine show_tile_geo(tile_geo, this_pe, var_name)
    type(grid_geometry)                    :: tile_geo
    integer, intent(in)                    :: this_pe
    character(len=*), intent(in)           :: var_name

    print '("[INFO] WDR 2Darray npe=",I0," ",A32, "nx=", I0," ny=", I0," nxp=",I0," nyp=",I0)', this_pe, var_name, tile_geo%nx, tile_geo%ny, tile_geo%nxp, tile_geo%nyp

    call check_array(tile_geo%lats, this_pe, var_name // "%lats", -90.0D0, 90.0D0)
    call check_array(tile_geo%lons, this_pe, var_name // "%lons", -360.0D0, 360.0D0)
    !call check_array(tile_geo%dx, this_pe, var_name // "%dx", 0.0, 1.0e9)
    !call check_array(tile_geo%dy, this_pe, var_name // "%dy", 0.0, 1.0e9)
    call check_array(tile_geo%area, this_pe, var_name // "%area", 0.0D0, 1.0D9)

  end subroutine show_tile_geo


  subroutine show_atm_array4(tag, array, array_name, atm_n, this_pe)
    character(len=*), intent(in)                                    :: tag
    real(kind=R_GRID), allocatable, dimension(:,:,:,:), intent(in)  :: array
    character(len=*), intent(in)                                    :: array_name
    integer, intent(in)                                             :: atm_n, this_pe

    if (allocated(array)) then
      print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%",A12,"(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', tag, this_pe, atm_n, trim(array_name), lbound(array, 1), ubound(array, 1), lbound(array, 2), ubound(array, 2), lbound(array, 3), ubound(array, 3), lbound(array, 4), ubound(array, 4)

    else
      print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%",A12," is not allocated.")', tag, this_pe, trim(array_name)
    endif

  end subroutine show_atm_array4


  subroutine show_atm_neststruct(tag, neststruct, atm_n, this_pe)
    character(len=*), intent(in)                                    :: tag
    type(fv_nest_type), intent(in) :: neststruct
    integer, intent(in)  :: atm_n, this_pe

    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%parent_tile=",I0," %refinement=",I0)', tag, this_pe, atm_n, neststruct%parent_tile, neststruct%refinement
    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%nested=",L1," %ioffset=",I0," %joffset=",I0)', tag, this_pe, atm_n, neststruct%nested, neststruct%ioffset, neststruct%joffset

    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%nested=",L1," %isu=",I0," %ieu=",I0," %jsu=",I0," %jeu=",I0)', tag, this_pe, atm_n, neststruct%nested, neststruct%isu,  neststruct%ieu,  neststruct%jsu,  neststruct%jeu

    ! WDR ind_update_h seems to have been removed in recent version of the dycore
    !    if (allocated(neststruct%ind_update_h)) then
    !       print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%ind_update_h(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', tag, this_pe, atm_n, &
    !            lbound(neststruct%ind_update_h,1),  ubound(neststruct%ind_update_h,1), &
    !            lbound(neststruct%ind_update_h,2),  ubound(neststruct%ind_update_h,2), &
    !            lbound(neststruct%ind_update_h,3),  ubound(neststruct%ind_update_h,3)
    !
    !       if (ubound(neststruct%ind_update_h,1) > lbound(neststruct%ind_update_h,1)) then
    !          print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") val neststruct%ind_update_h(",I0,",",I0,",",I0,")=",I0)', tag, this_pe, atm_n, &
    !               lbound(neststruct%ind_update_h,1),  lbound(neststruct%ind_update_h,2),  lbound(neststruct%ind_update_h,3), &
    !               neststruct%ind_update_h(lbound(neststruct%ind_update_h,1),  lbound(neststruct%ind_update_h,2),  lbound(neststruct%ind_update_h,3))
    !          print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") val neststruct%ind_update_h(",I0,",",I0,",",I0,")=",I0)', tag, this_pe, atm_n, &
    !               lbound(neststruct%ind_update_h,1),  lbound(neststruct%ind_update_h,2),  ubound(neststruct%ind_update_h,3), &
    !               neststruct%ind_update_h(lbound(neststruct%ind_update_h,1),  lbound(neststruct%ind_update_h,2),  ubound(neststruct%ind_update_h,3))
    !
    !          print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") val neststruct%ind_update_h(",I0,",",I0,",",I0,")=",I0)', tag, this_pe, atm_n, &
    !               lbound(neststruct%ind_update_h,1)+4,  lbound(neststruct%ind_update_h,2)+4,  lbound(neststruct%ind_update_h,3), &
    !               neststruct%ind_update_h(lbound(neststruct%ind_update_h,1)+4,  lbound(neststruct%ind_update_h,2)+4,  lbound(neststruct%ind_update_h,3))
    !          print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") val neststruct%ind_update_h(",I0,",",I0,",",I0,")=",I0)', tag, this_pe, atm_n, &
    !               lbound(neststruct%ind_update_h,1)+4,  lbound(neststruct%ind_update_h,2)+4,  ubound(neststruct%ind_update_h,3), &
    !               neststruct%ind_update_h(lbound(neststruct%ind_update_h,1)+4,  lbound(neststruct%ind_update_h,2)+4,  ubound(neststruct%ind_update_h,3))
    !
    !
    !       endif
    !    else
    !       print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%ind_update_h is not allocated.")', tag, this_pe, atm_n
    !    endif

    ! WDR nest_domain_all appears to be obsolete in new dycore
    !if (allocated(neststruct%nest_domain_all)) then
    !   print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%nest_domain_all(",I0,"-",I0,")")', tag, this_pe, atm_n, lbound(neststruct%nest_domain_all), ubound(neststruct%nest_domain_all)
    !else
    !   print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%nest_domain_all is not allocated.")', tag, this_pe, atm_n
    !endif

    ! WDR nest_domain has moved to fv_mp_mod.F90 as a global
    !print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%nest_domain%tile_fine=",I0," %tile_coarse=",I0)', tag, this_pe, atm_n, neststruct%nest_domain%tile_fine, neststruct%nest_domain%tile_coarse

    !print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%nest_domain%istart_fine=",I0," %iend_fine=",I0)', tag, this_pe, atm_n, neststruct%nest_domain%istart_fine,  neststruct%nest_domain%iend_fine
    !print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%nest_domain%jstart_fine=",I0," %jend_fine=",I0)', tag, this_pe, atm_n, neststruct%nest_domain%jstart_fine,  neststruct%nest_domain%jend_fine

    !print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%nest_domain%istart_coarse=",I0," %iend_coarse=",I0)', tag, this_pe, atm_n, neststruct%nest_domain%istart_coarse,  neststruct%nest_domain%iend_coarse
    !print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") neststruct%nest_domain%jstart_coarse=",I0," %jend_coarse=",I0)', tag, this_pe, atm_n, neststruct%nest_domain%jstart_coarse,  neststruct%nest_domain%jend_coarse

  end subroutine show_atm_neststruct


  subroutine show_atm_gridstruct(tag, gridstruct, atm_n, this_pe)
    character(len=*), intent(in)   :: tag
    type(fv_grid_type), intent(in) :: gridstruct
    integer, intent(in)            :: atm_n, this_pe

    ! nested is a pointer.
    if (associated(gridstruct%nested)) then
      print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") gridstruct%nested=",L1)', tag, this_pe, atm_n, gridstruct%nested
    else
      print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") gridstruct%nested is not set.")', tag, this_pe, atm_n
    endif

    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") gridstruct%cubed_sphere=",L1)', tag, this_pe, atm_n, gridstruct%cubed_sphere
    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") gridstruct%have_north_pole=",L1)', tag, this_pe, atm_n, gridstruct%have_north_pole
    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") gridstruct%have_south_pole=",L1)', tag, this_pe, atm_n, gridstruct%have_south_pole
    if (allocated(gridstruct%agrid)) then
      print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") gridstruct%agrid(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', tag, this_pe, atm_n, lbound(gridstruct%agrid, 1), ubound(gridstruct%agrid, 1), lbound(gridstruct%agrid, 2), ubound(gridstruct%agrid, 2), lbound(gridstruct%agrid, 3), ubound(gridstruct%agrid, 3)
    else
      print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") gridstruct%agrid is not allocated.")', tag, this_pe, atm_n
    endif

    if (allocated(gridstruct%grid)) then
      print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") gridstruct%grid(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', tag, this_pe, atm_n, lbound(gridstruct%grid, 1), ubound(gridstruct%grid, 1), lbound(gridstruct%grid, 2), ubound(gridstruct%grid, 2), lbound(gridstruct%grid, 3), ubound(gridstruct%grid, 3)
    else
      print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") gridstruct%grid is not allocated.")', tag, this_pe, atm_n
    endif

  end subroutine show_atm_gridstruct


  subroutine show_atm(tag, Atm, atm_n, this_pe)
    implicit none
    character(len=*), intent(in)     :: tag
    type(fv_atmos_type), intent(in)  :: Atm
    integer, intent(in)              :: atm_n, this_pe

    integer is, ie, i

    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,")===============================================================")', tag, this_pe, atm_n
    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") allocated=",L1," dummy=",L1)', tag, this_pe, atm_n, Atm%allocated, Atm%dummy
    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") grid_number=",I0," ncnst=",I0," ng=",I0)', tag, this_pe, atm_n, Atm%grid_number, Atm%ncnst, Atm%ng
    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") npx=",I0," npy=",I0," npz=",I0)', tag, this_pe, atm_n, Atm%npx, Atm%npy, Atm%npz

    if (allocated(Atm%pelist)) then
      is = lbound(Atm%pelist, 1)
      ie = ubound(Atm%pelist, 1)
      print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") pelist(",I0,"-",I0,")=",I0,"...",I0)', tag, this_pe, atm_n, is, ie, Atm%pelist(is),  Atm%pelist(ie)
      !do i = is, ie
      !   print '("[INFO]    show_atm ",A8," npe=",I0," Atm(",I0,") pelist(",I0,")=",I0)', tag, this_pe, atm_n, i, Atm%pelist(i)
      !enddo
    else
      print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") pelist is not allocated.")', tag, this_pe, atm_n
    endif

    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") bd%(is-ie)=",I0,"-",I0,") (js-je)=",I0,"-",I0,")"  )', tag, this_pe, atm_n, Atm%bd%is, Atm%bd%ie, Atm%bd%js, Atm%bd%je
    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") bd%(isd-ied)=",I0,"-",I0,") (jsd-jed)=",I0,"-",I0,")"  )', tag, this_pe, atm_n, Atm%bd%isd, Atm%bd%ied, Atm%bd%jsd, Atm%bd%jed
    print '("[INFO] show_atm ",A8," npe=",I0," Atm(",I0,") bd%(isc-iec)=",I0,"-",I0,") (jsc-jec)=",I0,"-",I0,")"  )', tag, this_pe, atm_n, Atm%bd%isc, Atm%bd%iec, Atm%bd%jsc, Atm%bd%jec

    call show_atm_neststruct(tag, Atm%neststruct, atm_n, this_pe)
    call show_atm_gridstruct(tag, Atm%gridstruct, atm_n, this_pe)
    call show_atm_array4(tag, Atm%grid_global, "grid_global", atm_n, this_pe)

  end subroutine show_atm


  subroutine show_gridstruct(gridstruct, this_pe)
    type(fv_grid_type), intent(in)     :: gridstruct
    integer, intent(in)                :: this_pe

    !real(kind=R_GRID)   :: pi = 4 * atan(1.0d0)
    real                 :: pi = 4 * atan(1.0d0)

    call check_array(gridstruct%grid, this_pe, "SG gridstruct%grid", -2.0*pi, 2.0*pi)
    call check_array(gridstruct%agrid, this_pe, "SG gridstruct%agrid", -2.0*pi, 2.0*pi)

    call check_array(gridstruct%area, this_pe, "SG gridstruct%area", 0.0, 1.0e12)
    call check_array(gridstruct%area_c, this_pe, "SG gridstruct%area_c", 0.0, 1.0e12)

    call check_array(gridstruct%rarea, this_pe, "SG gridstruct%rarea", 0.0, 1.0e12)
    call check_array(gridstruct%rarea_c, this_pe, "SG gridstruct%rarea_c", 0.0, 1.0e12)

    call check_array(gridstruct%sina, this_pe, "SG gridstruct%sina",  -1.0, 1.0)
    call check_array(gridstruct%cosa, this_pe, "SG gridstruct%cosa", -1.0, 1.0)

    call check_array(gridstruct%dx, this_pe, "SG gridstruct%dx", 0.0, 1.0e12)
    call check_array(gridstruct%dy, this_pe, "SG gridstruct%dy", 0.0, 1.0e12)

    call check_array(gridstruct%dxc, this_pe, "SG gridstruct%dxc", 0.0, 1.0e12)
    call check_array(gridstruct%dyc, this_pe, "SG gridstruct%dyc", 0.0, 1.0e12)

    call check_array(gridstruct%dxc_64, this_pe, "SG gridstruct%dxc_64", 0D0, 1.0D12)
    call check_array(gridstruct%dyc_64, this_pe, "SG gridstruct%dyc_64", 0D0, 1.0D12)

  end subroutine show_gridstruct


  subroutine show_nest_grid(Atm, this_pe, step_num)
    type(fv_atmos_type), intent(inout) :: Atm
    integer, intent(in)                :: this_pe, step_num

    integer             :: x,y
    integer             :: nhalo = 3  !! TODO get value from namelist
    real                :: crn_lat(4), crn_lon(4)
    real(kind=R_GRID)   :: pi = 4 * atan(1.0d0)
    real                :: pi180
    real                :: rad2deg, deg2rad

    pi180 = pi / 180.0
    deg2rad = pi / 180.0
    rad2deg = 1.0 / pi180

    print '("WDR NEST GRID bd, ",I0,",",I0," is,js=(",I0,":",I0,",",I0,":",I0,")"  )', &
        this_pe, step_num, Atm%bd%is, Atm%bd%ie, Atm%bd%js, Atm%bd%je

    print '("WDR NEST GRID bd, ",I0,",",I0," isd,jsd=(",I0,":",I0,",",I0,":",I0,")"  )', &
        this_pe, step_num, Atm%bd%isd, Atm%bd%ied, Atm%bd%jsd, Atm%bd%jed

    !do x = lbound(Atm%gridstruct%grid,1), ubound(Atm%gridstruct%grid,1)
    !   do y = lbound(Atm%gridstruct%grid,2), ubound(Atm%gridstruct%grid,2)
    !      print '("WDR NEST_GRID, ",I0,",",I0,",",I0,",",I0,",",F10.5,",",F10.5)', this_pe, step_num, x, y, &
    !           Atm%gridstruct%grid(x,y,2) * rad2deg, Atm%gridstruct%grid(x,y,1) * rad2deg - 360.0
    !   enddo
    !enddo

    ! Log the bounds of this PE's grid

    x = lbound(Atm%gridstruct%grid, 1)
    y = lbound(Atm%gridstruct%grid, 2)
    crn_lon(1) = Atm%gridstruct%grid(x,y,1)
    crn_lat(1) = Atm%gridstruct%grid(x,y,2)

    x = ubound(Atm%gridstruct%grid, 1)
    y = lbound(Atm%gridstruct%grid, 2)
    crn_lon(2) = Atm%gridstruct%grid(x,y,1)
    crn_lat(2) = Atm%gridstruct%grid(x,y,2)

    x = ubound(Atm%gridstruct%grid, 1)
    y = ubound(Atm%gridstruct%grid, 2)
    crn_lon(3) = Atm%gridstruct%grid(x,y,1)
    crn_lat(3) = Atm%gridstruct%grid(x,y,2)

    x = lbound(Atm%gridstruct%grid, 1)
    y = ubound(Atm%gridstruct%grid, 2)
    crn_lon(4) = Atm%gridstruct%grid(x,y,1)
    crn_lat(4) = Atm%gridstruct%grid(x,y,2)

    crn_lon(:) = crn_lon(:) * rad2deg
    crn_lat(:) = crn_lat(:) * rad2deg

    do x=1,4
      if (crn_lon(x) .gt. 180.0) then
        crn_lon(x) = crn_lon(x) - 360.0
      endif
    enddo

    print '("PLOT",I0,"_data_corners,",I4.4 ,",",F10.5,",",F10.5,",",F10.5,",",F10.5,",",F10.5,",",F10.5,",",F10.5,",",F10.5)', &
        step_num, this_pe, crn_lat(1),  crn_lon(1),  crn_lat(2),  crn_lon(2),  crn_lat(3),  crn_lon(3),  crn_lat(4),  crn_lon(4)

    ! Assume that nhalo is the same as all the other halo values
    x = lbound(Atm%gridstruct%grid, 1) + nhalo
    y = lbound(Atm%gridstruct%grid, 2) + nhalo
    crn_lon(1) = Atm%gridstruct%grid(x,y,1)
    crn_lat(1) = Atm%gridstruct%grid(x,y,2)

    x = ubound(Atm%gridstruct%grid, 1) - nhalo
    y = lbound(Atm%gridstruct%grid, 2) + nhalo
    crn_lon(2) = Atm%gridstruct%grid(x,y,1)
    crn_lat(2) = Atm%gridstruct%grid(x,y,2)

    x = ubound(Atm%gridstruct%grid, 1) - nhalo
    y = ubound(Atm%gridstruct%grid, 2) - nhalo
    crn_lon(3) = Atm%gridstruct%grid(x,y,1)
    crn_lat(3) = Atm%gridstruct%grid(x,y,2)

    x = lbound(Atm%gridstruct%grid, 1) + nhalo
    y = ubound(Atm%gridstruct%grid, 2) - nhalo
    crn_lon(4) = Atm%gridstruct%grid(x,y,1)
    crn_lat(4) = Atm%gridstruct%grid(x,y,2)

    crn_lon(:) = crn_lon(:) * rad2deg
    crn_lat(:) = crn_lat(:) * rad2deg

    do x=1,4
      if (crn_lon(x) .gt. 180.0) then
        crn_lon(x) = crn_lon(x) - 360.0
      endif
    enddo

    print '("PLOT",I0,"_compute_corners,",I4.4 ,",",F10.5,",",F10.5,",",F10.5,",",F10.5,",",F10.5,",",F10.5,",",F10.5,",",F10.5)', &
        step_num, this_pe, crn_lat(1),  crn_lon(1),  crn_lat(2),  crn_lon(2),  crn_lat(3),  crn_lon(3),  crn_lat(4),  crn_lon(4)

  end subroutine show_nest_grid


  subroutine validate_hires_parent(fp_super_tile_geo, grid, agrid, x_refine, y_refine, ioffset, joffset)
    type(grid_geometry), intent(in)                  :: fp_super_tile_geo
    real, allocatable, intent(in), dimension(:,:,:)  :: grid, agrid
    integer, intent(in)                              :: x_refine, y_refine, ioffset, joffset

    real, allocatable               :: local_grid(:,:,:), local_agrid(:,:,:)
    real(kind=R_GRID), allocatable  :: local_agrid64(:,:,:)
    logical                         :: is_equal
    integer                         :: x, y, z, this_pe, stagger
    real(kind=R_GRID)               :: pi = 4 * atan(1.0d0)
    real                            :: rad2deg

    rad2deg = 180.0 / pi
    this_pe = mpp_pe()

    !! Begin test creating of grid and agrid aligned with initial nest
    !! This is for testing/validation, and will not be needed in operations

    ! Allocate grid/agrid to proper size/bounds

    allocate(local_grid(lbound(grid,1) : ubound(grid,1), &
        lbound(grid,2) : ubound(grid,2), &
        lbound(grid,3) : ubound(grid,3)))

    allocate(local_agrid(lbound(agrid,1) : ubound(agrid,1), &
        lbound(agrid,2) : ubound(agrid,2), &
        lbound(agrid,3) : ubound(agrid,3)))

    allocate(local_agrid64(lbound(agrid,1) : ubound(agrid,1), &
        lbound(agrid,2) : ubound(agrid,2), &
        lbound(agrid,3) : ubound(agrid,3)))

    ! Fill in values from high resolution, full panel, supergrid

    stagger = CORNER
    call fill_grid_from_supergrid(local_grid, stagger, fp_super_tile_geo, ioffset, joffset, &
        x_refine, y_refine)
    stagger = CENTER
    call fill_grid_from_supergrid(local_agrid, stagger, fp_super_tile_geo, ioffset, joffset, &
        x_refine, y_refine)
    stagger = CENTER
    call fill_grid_from_supergrid(local_agrid64, stagger, fp_super_tile_geo, ioffset, joffset, &
        x_refine, y_refine)

    ! Verify that values are equivalent to the unmodified values in gridstruct

    call grid_equal(local_grid, grid, "GRID", this_pe, is_equal)
    call grid_equal(local_agrid, agrid, "AGRID", this_pe, is_equal)

    do x = lbound(grid,1), lbound(grid,1)+4
      do y = lbound(grid,2), lbound(grid,2)+4
        do z = lbound(grid,3), ubound(grid,3)
          print '("[INFO]  WDR grid_comp ",A16," npe=",I0," DEG value at (",I0,",",I0,",",I0,") ",F15.11, " ",F15.11, " ",F15.11)', "GRID", this_pe, x, y, z, local_grid(x,y,z)*rad2deg, grid(x,y,z)*rad2deg, local_grid(x,y,z)*rad2deg - grid(x,y,z)*rad2deg
          print '("[INFO]  WDR grid_comp ",A16," npe=",I0," RAD value at (",I0,",",I0,",",I0,") ",F15.11, " ",F15.11, " ",F15.11)', "GRID", this_pe, x, y, z, local_grid(x,y,z), grid(x,y,z), local_grid(x,y,z) - grid(x,y,z)
        enddo
      enddo
    enddo

    do x = lbound(agrid,1), lbound(agrid,1)+4
      do y = lbound(agrid,2), lbound(agrid,2)+4
        do z = lbound(agrid,3), ubound(agrid,3)
          print '("[INFO]  WDR agrid_comp ",A16," npe=",I0," DEG value at (",I0,",",I0,",",I0,") ",F15.11, " ",F15.11, " ",F15.11)', "AGRID", this_pe, x, y, z, local_agrid(x,y,z)*rad2deg, agrid(x,y,z)*rad2deg, local_agrid(x,y,z)*rad2deg - agrid(x,y,z)*rad2deg
          print '("[INFO]  WDR agrid_comp ",A16," npe=",I0," RAD value at (",I0,",",I0,",",I0,") ",F15.11, " ",F15.11, " ",F15.11)', "AGRID", this_pe, x, y, z, local_agrid(x,y,z), agrid(x,y,z), local_agrid(x,y,z) - agrid(x,y,z)
        enddo
      enddo
    enddo

    ! Validate at the end
    !call find_nest_alignment(tile_geo, fp_super_tile_geo, nest_x, nest_y, parent_x, parent_y)

  end subroutine validate_hires_parent

#endif ! MOVING_NEST

end module fv_moving_nest_utils_mod
