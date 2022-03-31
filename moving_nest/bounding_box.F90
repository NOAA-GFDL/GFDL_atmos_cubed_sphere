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
!! @brief Provides subroutines for grid bounding boxes for moving nest
!! @author W. Ramstrom, AOML/HRD  07/28/2021
!! @email William.Ramstrom@noaa.gov
!=======================================================================!


module bounding_box_mod
  use mpp_domains_mod, only : mpp_get_C2F_index, nest_domain_type
  use mpp_mod,         only : mpp_pe
  use fv_arrays_mod,   only : R_GRID

#ifdef GFS_TYPES
  use GFS_typedefs,      only : kind_phys
#else
  use IPD_typedefs,      only : kind_phys => IPD_kind_phys
#endif

  ! Simple aggregation of the start and end indices of a 2D grid
  ! Makes argument lists clearer to read
  type bbox
    integer :: is, ie, js, je
  end type bbox

  interface fill_bbox
    module procedure fill_bbox_r4_2d
    module procedure fill_bbox_r4_3d
    module procedure fill_bbox_r4_4d
    module procedure fill_bbox_r8_2d
    module procedure fill_bbox_r8_3d
    module procedure fill_bbox_r8_4d
  end interface fill_bbox

contains

  subroutine fill_bbox_r4_2d(out_bbox, in_grid)
    type(bbox), intent(out)         :: out_bbox
    real*4, allocatable, intent(in) :: in_grid(:,:)

    out_bbox%is = lbound(in_grid, 1)
    out_bbox%ie = ubound(in_grid, 1)
    out_bbox%js = lbound(in_grid, 2)
    out_bbox%je = ubound(in_grid, 2)
  end subroutine fill_bbox_r4_2d


  subroutine fill_bbox_r4_3d(out_bbox, in_grid)
    type(bbox), intent(out)         :: out_bbox
    real*4, allocatable, intent(in) :: in_grid(:,:,:)

    out_bbox%is = lbound(in_grid, 1)
    out_bbox%ie = ubound(in_grid, 1)
    out_bbox%js = lbound(in_grid, 2)
    out_bbox%je = ubound(in_grid, 2)
  end subroutine fill_bbox_r4_3d

  subroutine fill_bbox_r4_4d(out_bbox, in_grid)
    type(bbox), intent(out)         :: out_bbox
    real*4, allocatable, intent(in) :: in_grid(:,:,:,:)

    out_bbox%is = lbound(in_grid, 1)
    out_bbox%ie = ubound(in_grid, 1)
    out_bbox%js = lbound(in_grid, 2)
    out_bbox%je = ubound(in_grid, 2)
  end subroutine fill_bbox_r4_4d


  subroutine fill_bbox_r8_2d(out_bbox, in_grid)
    type(bbox), intent(out)         :: out_bbox
    real*8, allocatable, intent(in) :: in_grid(:,:)

    out_bbox%is = lbound(in_grid, 1)
    out_bbox%ie = ubound(in_grid, 1)
    out_bbox%js = lbound(in_grid, 2)
    out_bbox%je = ubound(in_grid, 2)
  end subroutine fill_bbox_r8_2d

  subroutine fill_bbox_r8_3d(out_bbox, in_grid)
    type(bbox), intent(out)         :: out_bbox
    real*8, allocatable, intent(in) :: in_grid(:,:,:)

    out_bbox%is = lbound(in_grid, 1)
    out_bbox%ie = ubound(in_grid, 1)
    out_bbox%js = lbound(in_grid, 2)
    out_bbox%je = ubound(in_grid, 2)
  end subroutine fill_bbox_r8_3d


  subroutine fill_bbox_r8_4d(out_bbox, in_grid)
    type(bbox), intent(out)          :: out_bbox
    real*8, allocatable, intent(in)  :: in_grid(:,:,:,:)

    out_bbox%is = lbound(in_grid, 1)
    out_bbox%ie = ubound(in_grid, 1)
    out_bbox%js = lbound(in_grid, 2)
    out_bbox%je = ubound(in_grid, 2)
  end subroutine fill_bbox_r8_4d

  subroutine show_bbox(tag, in_bbox, lats, lons)
    character(len=*) :: tag
    type(bbox), intent(out) :: in_bbox
    real(kind=kind_phys), allocatable, intent(in)        :: lats(:,:), lons(:,:)

    integer  :: x,y
    real(kind=R_GRID) :: pi = 4 * atan(1.0d0)
    real                :: pi180
    real                :: rad2deg, deg2rad

    pi180 = pi / 180.0
    deg2rad = pi / 180.0
    rad2deg = 1.0 / pi180

    x = in_bbox%is
    y = in_bbox%js
    !print '("[INFO] WDR show_bbox ",A8," lats(",I0,",",I0,")=",F10.5," lons(",I0,",",I0,")="F10.5)', tag, x, y, lats(x,y)*rad2deg,  x, y, lons(x,y)*rad2deg
    x = in_bbox%ie
    y = in_bbox%js
    !print '("[INFO] WDR show_bbox ",A8," lats(",I0,",",I0,")=",F10.5," lons(",I0,",",I0,")="F10.5)', tag, x, y, lats(x,y)*rad2deg,  x, y, lons(x,y)*rad2deg
    x = in_bbox%is
    y = in_bbox%je
    !print '("[INFO] WDR show_bbox ",A8," lats(",I0,",",I0,")=",F10.5," lons(",I0,",",I0,")="F10.5)', tag, x, y, lats(x,y)*rad2deg,  x, y, lons(x,y)*rad2deg
    x = in_bbox%ie
    y = in_bbox%je
    !print '("[INFO] WDR show_bbox ",A8," lats(",I0,",",I0,")=",F10.5," lons(",I0,",",I0,")="F10.5)', tag, x, y, lats(x,y)*rad2deg,  x, y, lons(x,y)*rad2deg

  end subroutine show_bbox


  !>@brief This subroutine returns the nest grid indices that correspond to the input nest domain, direction, and position
  !>@details  Simplifies the call signature with the bbox type rather than 4 separate integers
  subroutine bbox_get_C2F_index(nest_domain, bbox_fine, bbox_coarse, direction,  position)
    implicit none
    type(nest_domain_type), intent(in)     :: nest_domain
    type(bbox), intent(out) :: bbox_fine, bbox_coarse
    integer, intent(in) :: direction, position

    integer        :: this_pe
    integer        :: nest_level = 1   ! WDR TODO allow to vary
    this_pe = mpp_pe()

    !print '("[INFO] WDR enter bbox_get_C2F_index npe=",I0)', this_pe

    call mpp_get_C2F_index(nest_domain, bbox_fine%is, bbox_fine%ie, bbox_fine%js, bbox_fine%je, &
        bbox_coarse%is, bbox_coarse%ie, bbox_coarse%js, bbox_coarse%je, direction,  nest_level, position=position)

    !print '("[INFO] WDR bbox_get_C2F_index npe=",I0," dir=",I0," pos=",I0," fine (",I0,"-",I0,",",I0,"-",I0,")")', this_pe, direction, position, bbox_fine%is, bbox_fine%ie, bbox_fine%js, bbox_fine%je

    !print '("[INFO] WDR bbox_get_C2F_index npe=",I0," dir=",I0," pos=",I0," coarse (",I0,"-",I0,",",I0,"-",I0,")")', this_pe, direction, position, bbox_coarse%is, bbox_coarse%ie, bbox_coarse%js, bbox_coarse%je

  end subroutine bbox_get_C2F_index

end module bounding_box_mod
