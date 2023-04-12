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

module coarse_graining_mod

  use fms_mod, only: check_nml_error
  use mpp_domains_mod, only: domain2d, mpp_define_io_domain, mpp_define_mosaic, mpp_get_compute_domain
  use mpp_mod, only: FATAL, input_nml_file, mpp_error, mpp_npes
  use statistics_mod, only: mode

  implicit none
  private

  public :: block_sum, compute_mass_weights, get_fine_array_bounds, &
       get_coarse_array_bounds, coarse_graining_init, weighted_block_average, &
       weighted_block_edge_average_x, weighted_block_edge_average_y, MODEL_LEVEL, &
       block_upsample, mask_area_weights, PRESSURE_LEVEL, vertical_remapping_requirements, &
       vertically_remap_field, remap_edges_along_x, remap_edges_along_y, &
       block_edge_sum_x, block_edge_sum_y, block_mode, block_min, block_max, &
       eddy_covariance, eddy_covariance_2d_weights, eddy_covariance_3d_weights

  interface block_sum
     module procedure block_sum_2d_real4
     module procedure block_sum_2d_real8
     module procedure masked_block_sum_2d_real4
     module procedure masked_block_sum_2d_real8
  end interface block_sum

  interface block_edge_sum_x
     module procedure block_edge_sum_x_2d_full_input
  end interface block_edge_sum_x

  interface block_edge_sum_y
     module procedure block_edge_sum_y_2d_full_input
  end interface block_edge_sum_y

  interface weighted_block_average
     module procedure weighted_block_average_2d_real4
     module procedure weighted_block_average_2d_real8
     module procedure masked_weighted_block_average_2d_real4
     module procedure masked_weighted_block_average_2d_real8
     module procedure masked_weighted_block_average_3d_field_2d_weights_real4
     module procedure masked_weighted_block_average_3d_field_2d_weights_real8
     module procedure weighted_block_average_3d_field_2d_weights_real4
     module procedure weighted_block_average_3d_field_2d_weights_real8
     module procedure weighted_block_average_3d_field_3d_weights_real4
     module procedure weighted_block_average_3d_field_3d_weights_real8
  end interface weighted_block_average

  interface weighted_block_edge_average_x
     module procedure weighted_block_edge_average_x_2d
     module procedure weighted_block_edge_average_x_3d_field_2d_weights
  end interface weighted_block_edge_average_x

  interface weighted_block_edge_average_y
     module procedure weighted_block_edge_average_y_2d
     module procedure weighted_block_edge_average_y_3d_field_2d_weights
  end interface weighted_block_edge_average_y

  interface block_upsample
     module procedure block_upsample_2d_real4
     module procedure block_upsample_3d_real4
     module procedure block_upsample_2d_real8
     module procedure block_upsample_3d_real8
  end interface block_upsample

  interface weighted_block_edge_average_x_pre_downsampled
     module procedure weighted_block_edge_average_x_pre_downsampled_unmasked
     module procedure weighted_block_edge_average_x_pre_downsampled_masked
  end interface weighted_block_edge_average_x_pre_downsampled

  interface weighted_block_edge_average_y_pre_downsampled
     module procedure weighted_block_edge_average_y_pre_downsampled_unmasked
     module procedure weighted_block_edge_average_y_pre_downsampled_masked
  end interface weighted_block_edge_average_y_pre_downsampled

  interface eddy_covariance
     module procedure eddy_covariance_2d_weights
     module procedure eddy_covariance_3d_weights
  end interface eddy_covariance

  interface block_mode
     module procedure block_mode_2d_real4
     module procedure masked_block_mode_2d_real4
     module procedure block_mode_2d_real8
     module procedure masked_block_mode_2d_real8
  end interface block_mode

  interface block_min
     module procedure masked_block_min_2d_real4
     module procedure masked_block_min_2d_real8
  end interface block_min

  interface block_max
     module procedure masked_block_max_2d_real4
     module procedure masked_block_max_2d_real8
  end interface block_max

  interface vertical_remapping_requirements
     module procedure vertical_remapping_requirements_real4
     module procedure vertical_remapping_requirements_real8
  end interface vertical_remapping_requirements

  interface compute_phalf_from_delp
     module procedure compute_phalf_from_delp_real4
     module procedure compute_phalf_from_delp_real8
  end interface compute_phalf_from_delp

  interface mask_area_weights
     module procedure mask_area_weights_real4
     module procedure mask_area_weights_real8
  end interface mask_area_weights

  interface vertically_remap_field
     module procedure vertically_remap_field_real4
     module procedure vertically_remap_field_real8
  end interface vertically_remap_field

  interface mappm
     module procedure mappm_real4
     module procedure mappm_real8
  end interface mappm

  interface cs_profile
     module procedure cs_profile_real4
     module procedure cs_profile_real8
  end interface cs_profile

  interface cs_limiters
     module procedure cs_limiters_real4
     module procedure cs_limiters_real8
  end interface cs_limiters

  interface ppm_profile
     module procedure ppm_profile_real4
     module procedure ppm_profile_real8
  end interface ppm_profile

  interface ppm_limiters
     module procedure ppm_limiters_real4
     module procedure ppm_limiters_real8
  end interface ppm_limiters

  ! Global variables for the module, initialized in coarse_graining_init
  integer :: is, ie, js, je, npz
  integer :: is_coarse, ie_coarse, js_coarse, je_coarse
  character(len=11) :: MODEL_LEVEL = 'model_level'
  character(len=14) :: PRESSURE_LEVEL = 'pressure_level'

  ! GLobal variables for mappm
  real(kind=4), parameter:: r3_real4 = 1./3., r23_real4 = 2./3., r12_real4 = 1./12.
  real(kind=8), parameter:: r3_real8 = 1./3., r23_real8 = 2./3., r12_real8 = 1./12.

  ! Namelist parameters initialized with default values
  integer :: coarsening_factor = 8  !< factor the coarse grid is downsampled by (e.g. 8 if coarsening from C384 to C48 resolution)
  integer :: coarse_io_layout(2) = (/1, 1/)  !< I/O layout for coarse-grid fields
  character(len=64) :: strategy = 'model_level'  !< Valid values are 'model_level' and 'pressure_level'

  namelist /coarse_graining_nml/ coarsening_factor, coarse_io_layout, strategy

contains

  subroutine coarse_graining_init(npx, atm_npz, layout, is_fine, ie_fine, &
       js_fine, je_fine, factor, nx_coarse, coarse_graining_strategy, coarse_domain)
    integer, intent(in) :: npx
    integer, intent(in) :: atm_npz
    integer, intent(in) :: layout(2)
    integer, intent(in) :: is_fine, ie_fine, js_fine, je_fine
    integer, intent(out) :: factor
    integer, intent(out) :: nx_coarse
    character(len=64), intent(out) :: coarse_graining_strategy
    type(domain2d), intent(out) :: coarse_domain

    character(len=256) :: error_message
    logical :: exists
    integer :: error_code, iostat

    read(input_nml_file, coarse_graining_nml, iostat=iostat)
    error_code = check_nml_error(iostat, 'coarse_graining_nml')

    call assert_valid_strategy(strategy)
    call compute_nx_coarse(npx, coarsening_factor, nx_coarse)
    call assert_valid_domain_layout(nx_coarse, layout)
    call define_cubic_mosaic(coarse_domain, nx_coarse, nx_coarse, layout)
    call mpp_define_io_domain(coarse_domain, coarse_io_layout)
    call mpp_get_compute_domain(coarse_domain, is_coarse, ie_coarse, js_coarse, je_coarse)
    call set_fine_array_bounds(is_fine, ie_fine, js_fine, je_fine)
    npz = atm_npz
    factor = coarsening_factor
    coarse_graining_strategy = strategy
  end subroutine coarse_graining_init

  subroutine compute_nx_coarse(npx, coarsening_factor, nx_coarse)
    integer, intent(in) :: npx
    integer, intent(in) :: coarsening_factor
    integer, intent(out) :: nx_coarse

    character(len=256) :: error_message
    integer :: nx

    nx = npx - 1
    if (mod(nx, coarsening_factor) > 0) then
       write(error_message, *) 'coarse_graining_init: coarsening_factor does not evenly divide the native resolution.'
       call mpp_error(FATAL, error_message)
    endif
    nx_coarse = nx / coarsening_factor
  end subroutine compute_nx_coarse

  subroutine assert_valid_domain_layout(nx_coarse, layout)
    integer, intent(in) :: nx_coarse
    integer, intent(in) :: layout(2)

    character(len=256) :: error_message
    integer :: layout_x, layout_y
    layout_x = layout(1)
    layout_y = layout(2)

    if (mod(nx_coarse, layout_x) > 0 .or. mod(nx_coarse, layout_y) > 0) then
       write(error_message, *) 'coarse_graining_init: domain decomposition layout does not evenly divide the coarse grid.'
       call mpp_error(FATAL, error_message)
    endif
  end subroutine assert_valid_domain_layout

  subroutine assert_valid_strategy(strategy)
    character(len=64), intent(in) :: strategy

    character(len=256) :: error_message

    if (trim(strategy) .ne. MODEL_LEVEL .and. trim(strategy) .ne. PRESSURE_LEVEL) then
       write(error_message, *) 'Invalid coarse graining strategy provided.'
       call mpp_error(FATAL, error_message)
    endif
  end subroutine assert_valid_strategy

  subroutine set_fine_array_bounds(is_in, ie_in, js_in, je_in)
    integer, intent(in) :: is_in, ie_in, js_in, je_in

    is = is_in
    ie = ie_in
    js = js_in
    je = je_in
  end subroutine set_fine_array_bounds

  subroutine get_fine_array_bounds(is_out, ie_out, js_out, je_out)
    integer, intent(out) :: is_out, ie_out, js_out, je_out

    is_out = is
    ie_out = ie
    js_out = js
    je_out = je
  end subroutine get_fine_array_bounds

  subroutine get_coarse_array_bounds(is_out, ie_out, js_out, je_out)
    integer, intent(out) :: is_out, ie_out, js_out, je_out

    is_out = is_coarse
    ie_out = ie_coarse
    js_out = js_coarse
    je_out = je_coarse
  end subroutine get_coarse_array_bounds

  subroutine compute_mass_weights(area, delp, mass)
    real, intent(in) :: area(is:ie,js:je)
    real, intent(in) :: delp(is:ie,js:je,1:npz)
    real, intent(out) :: mass(is:ie,js:je,1:npz)

    integer :: k

    do k = 1, npz
       mass(:,:,k) = area * delp(:,:,k)
    enddo
  end subroutine compute_mass_weights

  subroutine block_sum_2d_real4(fine, coarse)
    real(kind=4), intent(in) :: fine(is:ie,js:je)
    real(kind=4), intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse)

    integer :: i, j, i_coarse, j_coarse, offset

    offset = coarsening_factor - 1
    do i = is, ie, coarsening_factor
       i_coarse = (i - 1) / coarsening_factor + 1
       do j = js, je, coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          coarse(i_coarse,j_coarse) = sum(fine(i:i+offset,j:j+offset))
       enddo
    enddo
  end subroutine block_sum_2d_real4

  subroutine block_sum_2d_real8(fine, coarse)
    real(kind=8), intent(in) :: fine(is:ie,js:je)
    real(kind=8), intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse)

    integer :: i, j, i_coarse, j_coarse, offset

    offset = coarsening_factor - 1
    do i = is, ie, coarsening_factor
       i_coarse = (i - 1) / coarsening_factor + 1
       do j = js, je, coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          coarse(i_coarse,j_coarse) = sum(fine(i:i+offset,j:j+offset))
       enddo
    enddo
  end subroutine block_sum_2d_real8

  subroutine masked_block_sum_2d_real4(fine, mask, coarse)
    real(kind=4), intent(in) :: fine(is:ie,js:je)
    logical, intent(in) :: mask(is:ie,js:je)
    real(kind=4), intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse)

    integer :: i, j, i_coarse, j_coarse, offset

    offset = coarsening_factor - 1
    do i = is, ie, coarsening_factor
       i_coarse = (i - 1) / coarsening_factor + 1
       do j = js, je, coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          coarse(i_coarse,j_coarse) = sum(fine(i:i+offset,j:j+offset), mask=mask(i:i+offset,j:j+offset))
       enddo
    enddo
  end subroutine masked_block_sum_2d_real4

  subroutine masked_block_sum_2d_real8(fine, mask, coarse)
    real(kind=8), intent(in) :: fine(is:ie,js:je)
    logical, intent(in) :: mask(is:ie,js:je)
    real(kind=8), intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse)

    integer :: i, j, i_coarse, j_coarse, offset

    offset = coarsening_factor - 1
    do i = is, ie, coarsening_factor
       i_coarse = (i - 1) / coarsening_factor + 1
       do j = js, je, coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          coarse(i_coarse,j_coarse) = sum(fine(i:i+offset,j:j+offset), mask=mask(i:i+offset,j:j+offset))
       enddo
    enddo
  end subroutine masked_block_sum_2d_real8

  subroutine weighted_block_average_2d_real4(weights, fine, coarse)
    real(kind=4), intent(in) :: weights(is:ie,js:je), fine(is:ie,js:je)
    real(kind=4), intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse)

    real(kind=4), allocatable :: weighted_fine(:,:), weighted_block_sum(:,:), block_sum_weights(:,:)

    allocate(weighted_fine(is:ie,js:je))
    allocate(weighted_block_sum(is_coarse:ie_coarse,js_coarse:je_coarse))
    allocate(block_sum_weights(is_coarse:ie_coarse,js_coarse:je_coarse))

    weighted_fine = weights * fine
    call block_sum_2d_real4(weighted_fine, weighted_block_sum)
    call block_sum_2d_real4(weights, block_sum_weights)
    coarse = weighted_block_sum / block_sum_weights
  end subroutine weighted_block_average_2d_real4

  subroutine weighted_block_average_2d_real8(weights, fine, coarse)
    real(kind=8), intent(in) :: weights(is:ie,js:je), fine(is:ie,js:je)
    real(kind=8), intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse)

    real(kind=8), allocatable :: weighted_fine(:,:), weighted_block_sum(:,:), block_sum_weights(:,:)

    allocate(weighted_fine(is:ie,js:je))
    allocate(weighted_block_sum(is_coarse:ie_coarse,js_coarse:je_coarse))
    allocate(block_sum_weights(is_coarse:ie_coarse,js_coarse:je_coarse))

    weighted_fine = weights * fine
    call block_sum_2d_real8(weighted_fine, weighted_block_sum)
    call block_sum_2d_real8(weights, block_sum_weights)
    coarse = weighted_block_sum / block_sum_weights
  end subroutine weighted_block_average_2d_real8

  subroutine masked_weighted_block_average_2d_real4(weights, fine, mask, coarse)
    real(kind=4), intent(in) :: weights(is:ie,js:je), fine(is:ie,js:je)
    logical, intent(in) :: mask(is:ie,js:je)
    real(kind=4), intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse)

    real(kind=4), allocatable :: weighted_fine(:,:), weighted_block_sum(:,:), block_sum_weights(:,:)

    allocate(weighted_fine(is:ie,js:je))
    allocate(weighted_block_sum(is_coarse:ie_coarse,js_coarse:je_coarse))
    allocate(block_sum_weights(is_coarse:ie_coarse,js_coarse:je_coarse))

    weighted_fine = weights * fine
    call masked_block_sum_2d_real4(weighted_fine, mask, weighted_block_sum)
    call masked_block_sum_2d_real4(weights, mask, block_sum_weights)
    coarse = weighted_block_sum / block_sum_weights
  end subroutine masked_weighted_block_average_2d_real4

  subroutine masked_weighted_block_average_2d_real8(weights, fine, mask, coarse)
    real(kind=8), intent(in) :: weights(is:ie,js:je), fine(is:ie,js:je)
    logical, intent(in) :: mask(is:ie,js:je)
    real(kind=8), intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse)

    real(kind=8), allocatable :: weighted_fine(:,:), weighted_block_sum(:,:), block_sum_weights(:,:)

    allocate(weighted_fine(is:ie,js:je))
    allocate(weighted_block_sum(is_coarse:ie_coarse,js_coarse:je_coarse))
    allocate(block_sum_weights(is_coarse:ie_coarse,js_coarse:je_coarse))

    weighted_fine = weights * fine
    call masked_block_sum_2d_real8(weighted_fine, mask, weighted_block_sum)
    call masked_block_sum_2d_real8(weights, mask, block_sum_weights)
    coarse = weighted_block_sum / block_sum_weights
  end subroutine masked_weighted_block_average_2d_real8

  subroutine masked_weighted_block_average_3d_field_2d_weights_real4(weights, fine, mask, nz, coarse)
    real(kind=4), intent(in) :: weights(is:ie,js:je), fine(is:ie,js:je,1:nz)
    logical, intent(in) :: mask(is:ie,js:je)
    integer, intent(in) :: nz
    real(kind=4), intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:nz)

    integer :: k

    do k = 1, nz
       call masked_weighted_block_average_2d_real4(weights, fine(is:ie,js:je,k), mask, coarse(is_coarse:ie_coarse,js_coarse:je_coarse,k))
    enddo
  end subroutine masked_weighted_block_average_3d_field_2d_weights_real4

  subroutine masked_weighted_block_average_3d_field_2d_weights_real8(weights, fine, mask, nz, coarse)
    real(kind=8), intent(in) :: weights(is:ie,js:je), fine(is:ie,js:je,1:nz)
    logical, intent(in) :: mask(is:ie,js:je)
    integer, intent(in) :: nz
    real(kind=8), intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:nz)

    integer :: k

    do k = 1, nz
       call masked_weighted_block_average_2d_real8(weights, fine(is:ie,js:je,k), mask, coarse(is_coarse:ie_coarse,js_coarse:je_coarse,k))
    enddo
  end subroutine masked_weighted_block_average_3d_field_2d_weights_real8

  subroutine weighted_block_average_3d_field_2d_weights_real4(weights, fine, coarse)
    real(kind=4), intent(in) :: weights(is:ie,js:je), fine(is:ie,js:je,1:npz)
    real(kind=4), intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz)

    integer :: k

    do k = 1, npz
       call weighted_block_average_2d_real4(weights, fine(is:ie,js:je,k), coarse(is_coarse:ie_coarse,js_coarse:je_coarse,k))
    enddo
  end subroutine weighted_block_average_3d_field_2d_weights_real4

  subroutine weighted_block_average_3d_field_2d_weights_real8(weights, fine, coarse)
    real(kind=8), intent(in) :: weights(is:ie,js:je), fine(is:ie,js:je,1:npz)
    real(kind=8), intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz)

    integer :: k

    do k = 1, npz
       call weighted_block_average_2d_real8(weights, fine(is:ie,js:je,k), coarse(is_coarse:ie_coarse,js_coarse:je_coarse,k))
    enddo
  end subroutine weighted_block_average_3d_field_2d_weights_real8

  subroutine weighted_block_average_3d_field_3d_weights_real4(weights, fine, coarse)
    real(kind=4), intent(in) :: weights(is:ie,js:je,1:npz), fine(is:ie,js:je,1:npz)
    real(kind=4), intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz)

    integer :: k

    do k = 1, npz
       call weighted_block_average_2d_real4(weights(is:ie,js:je,k), fine(is:ie,js:je,k), coarse(is_coarse:ie_coarse,js_coarse:je_coarse,k))
    enddo
  end subroutine weighted_block_average_3d_field_3d_weights_real4

  subroutine weighted_block_average_3d_field_3d_weights_real8(weights, fine, coarse)
    real(kind=8), intent(in) :: weights(is:ie,js:je,1:npz), fine(is:ie,js:je,1:npz)
    real(kind=8), intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz)

    integer :: k

    do k = 1, npz
       call weighted_block_average_2d_real8(weights(is:ie,js:je,k), fine(is:ie,js:je,k), coarse(is_coarse:ie_coarse,js_coarse:je_coarse,k))
    enddo
  end subroutine weighted_block_average_3d_field_3d_weights_real8


  subroutine block_edge_sum_x_2d(fine, coarse)
    real, intent(in) :: fine(is:ie,js_coarse:je_coarse+1)
    real, intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse+1)

    integer :: i, i_coarse, j_coarse, offset

    offset = coarsening_factor - 1
    do i = is, ie, coarsening_factor
       i_coarse = (i - 1) / coarsening_factor + 1
       do j_coarse = js_coarse, je_coarse + 1
          coarse(i_coarse,j_coarse) = sum(fine(i:i+offset,j_coarse))
       enddo
    enddo
  end subroutine block_edge_sum_x_2d

  subroutine weighted_block_edge_average_x_2d(weights, fine, coarse)
    real, intent(in) :: weights(is:ie,js:je+1)
    real, intent(in) :: fine(is:ie,js:je+1)
    real, intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse+1)

    real, allocatable :: weighted_fine(:,:), weighted_block_sum(:,:), block_sum_weights(:,:)

    allocate(weighted_fine(is:ie,js_coarse:je_coarse+1))
    allocate(weighted_block_sum(is_coarse:ie_coarse,js_coarse:je_coarse+1))
    allocate(block_sum_weights(is_coarse:ie_coarse,js_coarse:je_coarse+1))

    weighted_fine = weights(is:ie,js:je+1:coarsening_factor) * fine(is:ie,js:je+1:coarsening_factor)
    call block_edge_sum_x_2d(weighted_fine, weighted_block_sum)
    call block_edge_sum_x_2d(weights(is:ie,js:je+1:coarsening_factor), block_sum_weights)
    coarse = weighted_block_sum / block_sum_weights
  end subroutine weighted_block_edge_average_x_2d

  subroutine weighted_block_edge_average_x_3d_field_2d_weights(weights, fine, coarse)
    real, intent(in) :: weights(is:ie,js:je+1)
    real, intent(in) :: fine(is:ie,js:je+1,1:npz)
    real, intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse+1,1:npz)

    integer :: k

    do k = 1, npz
       call weighted_block_edge_average_x_2d(weights, fine(is:ie,js:je+1,k), &
            coarse(is_coarse:ie_coarse,js_coarse:je_coarse+1,k))
    enddo
  end subroutine weighted_block_edge_average_x_3d_field_2d_weights

  subroutine block_edge_sum_y_2d(fine, coarse)
    real, intent(in) :: fine(is_coarse:ie_coarse+1,js:je)
    real, intent(out) :: coarse(is_coarse:ie_coarse+1,js_coarse:je_coarse)

    integer :: j, i_coarse, j_coarse, offset

    offset = coarsening_factor - 1
    do i_coarse = is_coarse, ie_coarse + 1
       do j = js, je, coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          coarse(i_coarse,j_coarse) = sum(fine(i_coarse,j:j+offset))
       enddo
    enddo
  end subroutine block_edge_sum_y_2d

  subroutine weighted_block_edge_average_y_2d(weights, fine, coarse)
    real, intent(in) :: weights(is:ie+1,js:je)
    real, intent(in) :: fine(is:ie+1,js:je)
    real, intent(out) :: coarse(is_coarse:ie_coarse+1,js_coarse:je_coarse)

    real, allocatable :: weighted_fine(:,:), weighted_block_sum(:,:), block_sum_weights(:,:)

    allocate(weighted_fine(is_coarse:ie_coarse+1,js:je))
    allocate(weighted_block_sum(is_coarse:ie_coarse+1,js_coarse:je_coarse))
    allocate(block_sum_weights(is_coarse:ie_coarse+1,js_coarse:je_coarse))

    weighted_fine = weights(is:ie+1:coarsening_factor,js:je) * fine(is:ie+1:coarsening_factor,js:je)
    call block_edge_sum_y_2d(weighted_fine, weighted_block_sum)
    call block_edge_sum_y_2d(weights(is:ie+1:coarsening_factor,js:je), block_sum_weights)
    coarse = weighted_block_sum / block_sum_weights
  end subroutine weighted_block_edge_average_y_2d

  subroutine weighted_block_edge_average_y_3d_field_2d_weights(weights, fine, coarse)
    real, intent(in) :: weights(is:ie+1,js:je)
    real, intent(in) :: fine(is:ie+1,js:je,1:npz)
    real, intent(out) :: coarse(is_coarse:ie_coarse+1,js_coarse:je_coarse,1:npz)

    integer :: k

    do k = 1, npz
       call weighted_block_edge_average_y_2d(weights, fine(is:ie+1,js:je,k), &
            coarse(is_coarse:ie_coarse+1,js_coarse:je_coarse,k))
    enddo
  end subroutine weighted_block_edge_average_y_3d_field_2d_weights

  subroutine block_mode_2d_real8(fine, coarse)
    real(kind=8), intent(in) :: fine(is:ie,js:je)
    real(kind=8), intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse)

    integer :: i, j, i_coarse, j_coarse, offset

    offset = coarsening_factor - 1
    do i = is, ie, coarsening_factor
       i_coarse = (i - 1) / coarsening_factor + 1
       do j = js, je, coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          coarse(i_coarse, j_coarse) = mode(fine(i:i+offset,j:j+offset))
       enddo
    enddo
  end subroutine block_mode_2d_real8

  subroutine masked_block_mode_2d_real8(fine, mask, coarse)
    real(kind=8), intent(in) :: fine(is:ie,js:je)
    logical, intent(in) :: mask(is:ie,js:je)
    real(kind=8), intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse)

    integer :: i, j, i_coarse, j_coarse, offset

    offset = coarsening_factor - 1
    do i = is, ie, coarsening_factor
       i_coarse = (i - 1) / coarsening_factor + 1
       do j = js, je, coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          coarse(i_coarse, j_coarse) = mode(fine(i:i+offset,j:j+offset), mask(i:i+offset,j:j+offset))
       enddo
    enddo
  end subroutine masked_block_mode_2d_real8

  subroutine block_mode_2d_real4(fine, coarse)
    real(kind=4), intent(in) :: fine(is:ie,js:je)
    real(kind=4), intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse)

    integer :: i, j, i_coarse, j_coarse, offset

    offset = coarsening_factor - 1
    do i = is, ie, coarsening_factor
       i_coarse = (i - 1) / coarsening_factor + 1
       do j = js, je, coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          coarse(i_coarse, j_coarse) = mode(fine(i:i+offset,j:j+offset))
       enddo
    enddo
  end subroutine block_mode_2d_real4

  subroutine masked_block_mode_2d_real4(fine, mask, coarse)
    real(kind=4), intent(in) :: fine(is:ie,js:je)
    logical, intent(in) :: mask(is:ie,js:je)
    real(kind=4), intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse)

    integer :: i, j, i_coarse, j_coarse, offset

    offset = coarsening_factor - 1
    do i = is, ie, coarsening_factor
       i_coarse = (i - 1) / coarsening_factor + 1
       do j = js, je, coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          coarse(i_coarse, j_coarse) = mode(fine(i:i+offset,j:j+offset), mask(i:i+offset,j:j+offset))
       enddo
    enddo
  end subroutine masked_block_mode_2d_real4

  subroutine masked_block_min_2d_real8(fine, mask, coarse)
    real(kind=8), intent(in) :: fine(is:ie,js:je)
    logical, intent(in) :: mask(is:ie,js:je)
    real(kind=8), intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse)

    integer :: i, j, i_coarse, j_coarse, offset

    offset = coarsening_factor - 1
    do i = is, ie, coarsening_factor
       i_coarse = (i - 1) / coarsening_factor + 1
       do j = js, je, coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          coarse(i_coarse, j_coarse) = minval(fine(i:i+offset,j:j+offset), mask=mask(i:i + offset,j:j + offset))
       enddo
    enddo
  end subroutine masked_block_min_2d_real8

  subroutine masked_block_max_2d_real8(fine, mask, coarse)
    real(kind=8), intent(in) :: fine(is:ie,js:je)
    logical, intent(in) :: mask(is:ie,js:je)
    real(kind=8), intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse)

    integer :: i, j, i_coarse, j_coarse, offset

    offset = coarsening_factor - 1
    do i = is, ie, coarsening_factor
       i_coarse = (i - 1) / coarsening_factor + 1
       do j = js, je, coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          coarse(i_coarse, j_coarse) = maxval(fine(i:i+offset,j:j+offset), mask=mask(i:i + offset,j:j + offset))
       enddo
    enddo
  end subroutine masked_block_max_2d_real8

  subroutine masked_block_min_2d_real4(fine, mask, coarse)
    real(kind=4), intent(in) :: fine(is:ie,js:je)
    logical, intent(in) :: mask(is:ie,js:je)
    real(kind=4), intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse)

    integer :: i, j, i_coarse, j_coarse, offset

    offset = coarsening_factor - 1
    do i = is, ie, coarsening_factor
       i_coarse = (i - 1) / coarsening_factor + 1
       do j = js, je, coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          coarse(i_coarse, j_coarse) = minval(fine(i:i+offset,j:j+offset), mask=mask(i:i + offset,j:j + offset))
       enddo
    enddo
  end subroutine masked_block_min_2d_real4

  subroutine masked_block_max_2d_real4(fine, mask, coarse)
    real(kind=4), intent(in) :: fine(is:ie,js:je)
    logical, intent(in) :: mask(is:ie,js:je)
    real(kind=4), intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse)

    integer :: i, j, i_coarse, j_coarse, offset

    offset = coarsening_factor - 1
    do i = is, ie, coarsening_factor
       i_coarse = (i - 1) / coarsening_factor + 1
       do j = js, je, coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          coarse(i_coarse, j_coarse) = maxval(fine(i:i+offset,j:j+offset), mask=mask(i:i + offset,j:j + offset))
       enddo
    enddo
  end subroutine masked_block_max_2d_real4

  subroutine vertically_remap_field_real4(phalf_in, field, phalf_out, ptop, field_out)
    real(kind=4), intent(in) :: phalf_in(is:ie,js:je,1:npz+1), phalf_out(is:ie,js:je,1:npz+1)
    real(kind=4), intent(in) :: field(is:ie,js:je,1:npz)
    real(kind=4), intent(in) :: ptop
    real(kind=4), intent(out) :: field_out(is:ie,js:je,1:npz)

    integer :: kn, km, kord, iv, j, q2

    kn = npz
    km = npz

    ! Hard code values of kord and iv for now
    kord = 1
    iv = 1
    q2 = 1

    do j = js, je
       call mappm(km, phalf_in(is:ie,j,:), field(is:ie,j,:), kn, &
            phalf_out(is:ie,j,:), field_out(is:ie,j,:), is, ie, iv, kord, ptop)
    enddo
  end subroutine vertically_remap_field_real4

  subroutine vertically_remap_field_real8(phalf_in, field, phalf_out, ptop, field_out)
    real(kind=8), intent(in) :: phalf_in(is:ie,js:je,1:npz+1), phalf_out(is:ie,js:je,1:npz+1)
    real(kind=8), intent(in) :: field(is:ie,js:je,1:npz)
    real(kind=8), intent(in) :: ptop
    real(kind=8), intent(out) :: field_out(is:ie,js:je,1:npz)

    integer :: kn, km, kord, iv, j, q2

    kn = npz
    km = npz

    ! Hard code values of kord and iv for now
    kord = 1
    iv = 1
    q2 = 1

    do j = js, je
       call mappm(km, phalf_in(is:ie,j,:), field(is:ie,j,:), kn, &
            phalf_out(is:ie,j,:), field_out(is:ie,j,:), is, ie, iv, kord, ptop)
    enddo
  end subroutine vertically_remap_field_real8

  subroutine block_upsample_2d_real4(coarse, fine)
    real(kind=4), intent(in) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse)
    real(kind=4), intent(out) :: fine(is:ie,js:je)

    integer :: i, j, i_coarse, j_coarse, offset

    offset = coarsening_factor - 1
    do i = is, ie, coarsening_factor
      i_coarse = (i - 1) / coarsening_factor + 1
      do j = js, je, coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          fine(i:i+offset,j:j+offset) = coarse(i_coarse, j_coarse)
      enddo
    enddo
  end subroutine block_upsample_2d_real4

  subroutine block_upsample_3d_real4(coarse, fine, nz)
    integer, intent(in) :: nz
    real(kind=4), intent(in) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:nz)
    real(kind=4), intent(out) :: fine(is:ie,js:je,1:nz)

    integer :: k

    do k = 1, nz
      call block_upsample_2d_real4(coarse(is_coarse:ie_coarse,js_coarse:je_coarse,k), fine(is:ie,js:je,k))
    enddo
  end subroutine block_upsample_3d_real4

  subroutine block_upsample_2d_real8(coarse, fine)
    real(kind=8), intent(in) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse)
    real(kind=8), intent(out) :: fine(is:ie,js:je)

    integer :: i, j, i_coarse, j_coarse, offset

    offset = coarsening_factor - 1
    do i = is, ie, coarsening_factor
      i_coarse = (i - 1) / coarsening_factor + 1
      do j = js, je, coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          fine(i:i+offset,j:j+offset) = coarse(i_coarse, j_coarse)
      enddo
    enddo
  end subroutine block_upsample_2d_real8

  subroutine block_upsample_3d_real8(coarse, fine, nz)
    integer, intent(in) :: nz
    real(kind=8), intent(in) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:nz)
    real(kind=8), intent(out) :: fine(is:ie,js:je,1:nz)

    integer :: k

    do k = 1, nz
      call block_upsample_2d_real8(coarse(is_coarse:ie_coarse,js_coarse:je_coarse,k), fine(is:ie,js:je,k))
    enddo
  end subroutine block_upsample_3d_real8

  ! This subroutine is copied from FMS/test_fms/horiz_interp/test2_horiz_interp.F90.
  ! domain_decomp in fv_mp_mod.F90 does something similar, but it does a
  ! few other unnecessary things (and requires more arguments).
  subroutine define_cubic_mosaic(domain, ni, nj, layout)
    type(domain2d), intent(inout) :: domain
    integer,        intent(in)    :: layout(2)
    integer,        intent(in)    :: ni, nj
    integer   :: pe_start(6), pe_end(6)
    integer   :: global_indices(4,6), layout2d(2,6)
    integer, dimension(12)        :: istart1, iend1, jstart1, jend1, tile1
    integer, dimension(12)        :: istart2, iend2, jstart2, jend2, tile2
    integer                       :: ntiles, num_contact
    integer :: p, npes_per_tile, i

    ntiles = 6
    num_contact = 12
    p = 0
    npes_per_tile = mpp_npes()/ntiles
    do i = 1, 6
       layout2d(:,i) = layout(:)
       global_indices(1,i) = 1
       global_indices(2,i) = ni
       global_indices(3,i) = 1
       global_indices(4,i) = nj
       pe_start(i) = p
       p = p + npes_per_tile
       pe_end(i) = p-1
    enddo

    !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
    tile1(1) = 1;     tile2(1) = 2
    istart1(1) = ni;  iend1(1) = ni;  jstart1(1) = 1;      jend1(1) = nj
    istart2(1) = 1;   iend2(1) = 1;   jstart2(1) = 1;      jend2(1) = nj

    !--- Contact line 2, between tile 1 (NORTH) and tile 3 (WEST)
    tile1(2) = 1; tile2(2) = 3
    istart1(2) = 1;      iend1(2) = ni;  jstart1(2) = nj;  jend1(2) = nj
    istart2(2) = 1;      iend2(2) = 1;   jstart2(2) = nj;  jend2(2) = 1

    !--- Contact line 3, between tile 1 (WEST) and tile 5 (NORTH)
    tile1(3) = 1;     tile2(3) = 5
    istart1(3) = 1;   iend1(3) = 1;      jstart1(3) = 1;   jend1(3) = nj
    istart2(3) = ni;  iend2(3) = 1;      jstart2(3) = nj;  jend2(3) = nj

    !--- Contact line 4, between tile 1 (SOUTH) and tile 6 (NORTH)
    tile1(4) = 1; tile2(4) = 6
    istart1(4) = 1;      iend1(4) = ni;  jstart1(4) = 1;   jend1(4) = 1
    istart2(4) = 1;      iend2(4) = ni;  jstart2(4) = nj;  jend2(4) = nj

    !--- Contact line 5, between tile 2 (NORTH) and tile 3 (SOUTH)
    tile1(5) = 2;        tile2(5) = 3
    istart1(5) = 1;      iend1(5) = ni;  jstart1(5) = nj;  jend1(5) = nj
    istart2(5) = 1;      iend2(5) = ni;  jstart2(5) = 1;   jend2(5) = 1

    !--- Contact line 6, between tile 2 (EAST) and tile 4 (SOUTH)
    tile1(6) = 2; tile2(6) = 4
    istart1(6) = ni;  iend1(6) = ni;  jstart1(6) = 1;      jend1(6) = nj
    istart2(6) = ni;  iend2(6) = 1;   jstart2(6) = 1;      jend2(6) = 1

    !--- Contact line 7, between tile 2 (SOUTH) and tile 6 (EAST)
    tile1(7) = 2; tile2(7) = 6
    istart1(7) = 1;   iend1(7) = ni;  jstart1(7) = 1;   jend1(7) = 1
    istart2(7) = ni;  iend2(7) = ni;  jstart2(7) = nj;  jend2(7) = 1

    !--- Contact line 8, between tile 3 (EAST) and tile 4 (WEST)
    tile1(8) = 3; tile2(8) = 4
    istart1(8) = ni;  iend1(8) = ni;  jstart1(8) = 1;      jend1(8) = nj
    istart2(8) = 1;   iend2(8) = 1;   jstart2(8) = 1;      jend2(8) = nj

    !--- Contact line 9, between tile 3 (NORTH) and tile 5 (WEST)
    tile1(9) = 3; tile2(9) = 5
    istart1(9) = 1;      iend1(9) = ni;  jstart1(9) = nj;  jend1(9) = nj
    istart2(9) = 1;      iend2(9) = 1;   jstart2(9) = nj;  jend2(9) = 1

    !--- Contact line 10, between tile 4 (NORTH) and tile 5 (SOUTH)
    tile1(10) = 4; tile2(10) = 5
    istart1(10) = 1;     iend1(10) = ni; jstart1(10) = nj; jend1(10) = nj
    istart2(10) = 1;     iend2(10) = ni; jstart2(10) = 1;  jend2(10) = 1

    !--- Contact line 11, between tile 4 (EAST) and tile 6 (SOUTH)
    tile1(11) = 4; tile2(11) = 6
    istart1(11) = ni; iend1(11) = ni; jstart1(11) = 1;     jend1(11) = nj
    istart2(11) = ni; iend2(11) = 1;  jstart2(11) = 1;     jend2(11) = 1

    !--- Contact line 12, between tile 5 (EAST) and tile 6 (WEST)
    tile1(12) = 5; tile2(12) = 6
    istart1(12) = ni; iend1(12) = ni; jstart1(12) = 1;     jend1(12) = nj
    istart2(12) = 1;  iend2(12) = 1;  jstart2(12) = 1;     jend2(12) = nj

    call mpp_define_mosaic(global_indices, layout2d, domain, ntiles, &
         num_contact, tile1, tile2, istart1, iend1, jstart1, jend1, &
         istart2, iend2, jstart2, jend2, pe_start, pe_end, &
         symmetry=.true., name='coarse cubic mosaic')
  end subroutine define_cubic_mosaic

  subroutine compute_phalf_from_delp_real4(delp, ptop, i_start, i_end, j_start, j_end, phalf)
    integer, intent(in) :: i_start, i_end, j_start, j_end
    real(kind=4), intent(in) :: delp(i_start:i_end,j_start:j_end,1:npz)
    real(kind=4), intent(in) :: ptop
    real(kind=4), intent(out) :: phalf(i_start:i_end,j_start:j_end,1:npz+1)

    integer :: i, j, k

    phalf(:,:,1) = ptop  ! Top level interface pressure is the model top

    ! Integrate delp from top of model to the surface.
    do i = i_start, i_end
       do j = j_start, j_end
          do k = 2, npz + 1
             phalf(i,j,k) = phalf(i,j,k-1) + delp(i,j,k-1)
          enddo
       enddo
    enddo
  end subroutine compute_phalf_from_delp_real4

  subroutine compute_phalf_from_delp_real8(delp, ptop, i_start, i_end, j_start, j_end, phalf)
    integer, intent(in) :: i_start, i_end, j_start, j_end
    real(kind=8), intent(in) :: delp(i_start:i_end,j_start:j_end,1:npz)
    real(kind=8), intent(in) :: ptop
    real(kind=8), intent(out) :: phalf(i_start:i_end,j_start:j_end,1:npz+1)

    integer :: i, j, k

    phalf(:,:,1) = ptop  ! Top level interface pressure is the model top

    ! Integrate delp from top of model to the surface.
    do i = i_start, i_end
       do j = j_start, j_end
          do k = 2, npz + 1
             phalf(i,j,k) = phalf(i,j,k-1) + delp(i,j,k-1)
          enddo
       enddo
    enddo
  end subroutine compute_phalf_from_delp_real8

 ! Routine for computing the common requirements for pressure-level coarse-graining.
  subroutine vertical_remapping_requirements_real4(delp, area, ptop, phalf, upsampled_coarse_phalf)
    real(kind=4), intent(in) :: delp(is:ie,js:je,1:npz)
    real(kind=4), intent(in) :: area(is:ie,js:je)
    real(kind=4), intent(in) :: ptop
    real(kind=4), intent(out) :: phalf(is:ie,js:je,1:npz+1)
    real(kind=4), intent(out) :: upsampled_coarse_phalf(is:ie,js:je,1:npz+1)

    real(kind=4), allocatable :: coarse_delp(:,:,:), coarse_phalf(:,:,:)

    allocate(coarse_delp(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz))
    allocate(coarse_phalf(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz+1))

    call compute_phalf_from_delp(delp(is:ie,js:je,1:npz), ptop, is, ie, js, je, phalf)
    call weighted_block_average(area(is:ie,js:je), delp(is:ie,js:je,1:npz), coarse_delp)
    call compute_phalf_from_delp(coarse_delp, ptop, is_coarse, ie_coarse, js_coarse, je_coarse, coarse_phalf)
    call block_upsample(coarse_phalf, upsampled_coarse_phalf, npz+1)

    deallocate(coarse_delp)
    deallocate(coarse_phalf)
   end subroutine vertical_remapping_requirements_real4

  subroutine vertical_remapping_requirements_real8(delp, area, ptop, phalf, upsampled_coarse_phalf)
    real(kind=8), intent(in) :: delp(is:ie,js:je,1:npz)
    real(kind=8), intent(in) :: area(is:ie,js:je)
    real(kind=8), intent(in) :: ptop
    real(kind=8), intent(out) :: phalf(is:ie,js:je,1:npz+1)
    real(kind=8), intent(out) :: upsampled_coarse_phalf(is:ie,js:je,1:npz+1)

    real(kind=8), allocatable :: coarse_delp(:,:,:), coarse_phalf(:,:,:)

    allocate(coarse_delp(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz))
    allocate(coarse_phalf(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz+1))

    call compute_phalf_from_delp(delp(is:ie,js:je,1:npz), ptop, is, ie, js, je, phalf)
    call weighted_block_average(area(is:ie,js:je), delp(is:ie,js:je,1:npz), coarse_delp)
    call compute_phalf_from_delp(coarse_delp, ptop, is_coarse, ie_coarse, js_coarse, je_coarse, coarse_phalf)
    call block_upsample(coarse_phalf, upsampled_coarse_phalf, npz+1)

    deallocate(coarse_delp)
    deallocate(coarse_phalf)
   end subroutine vertical_remapping_requirements_real8

   subroutine mask_area_weights_real4(area, phalf, upsampled_coarse_phalf, masked_area_weights)
    real(kind=4), intent(in) :: area(is:ie,js:je)
    real(kind=4), intent(in) :: phalf(is:ie,js:je,1:npz+1)
    real(kind=4), intent(in) :: upsampled_coarse_phalf(is:ie,js:je,1:npz+1)
    real(kind=4), intent(out) :: masked_area_weights(is:ie,js:je,1:npz)

    integer :: k

    do k = 1, npz
      where (upsampled_coarse_phalf(is:ie,js:je,k+1) .lt. phalf(is:ie,js:je,npz+1))
        masked_area_weights(is:ie,js:je,k) = area(is:ie,js:je)
      elsewhere
        masked_area_weights(is:ie,js:je,k) = 0.0
      endwhere
    enddo
   end subroutine mask_area_weights_real4

   subroutine mask_area_weights_real8(area, phalf, upsampled_coarse_phalf, masked_area_weights)
    real(kind=8), intent(in) :: area(is:ie,js:je)
    real(kind=8), intent(in) :: phalf(is:ie,js:je,1:npz+1)
    real(kind=8), intent(in) :: upsampled_coarse_phalf(is:ie,js:je,1:npz+1)
    real(kind=8), intent(out) :: masked_area_weights(is:ie,js:je,1:npz)

    integer :: k

    do k = 1, npz
      where (upsampled_coarse_phalf(is:ie,js:je,k+1) .lt. phalf(is:ie,js:je,npz+1))
        masked_area_weights(is:ie,js:je,k) = area(is:ie,js:je)
      elsewhere
        masked_area_weights(is:ie,js:je,k) = 0.0
      endwhere
    enddo
   end subroutine mask_area_weights_real8

   ! A naive routine for interpolating a field from the A-grid to the y-boundary
   ! of the D-grid; this is a specialized function that automatically
   ! downsamples to the coarse-grid on the downsampling dimension.
   subroutine interpolate_to_d_grid_and_downsample_along_y(field_in, field_out, nz)
    integer, intent(in) :: nz
    real, intent(in) :: field_in(is-1:ie+1,js-1:je+1,1:nz)
    real, intent(out) :: field_out(is:ie,js_coarse:je_coarse+1,1:nz)

    integer :: i, j, k, j_coarse

    do i = is,ie
       do j = js,je+1,coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          do k = 1,nz
             field_out(i,j_coarse,k) = 0.5 * (field_in(i,j,k) + field_in(i,j-1,k))
          enddo
       enddo
    enddo
   end subroutine interpolate_to_d_grid_and_downsample_along_y

   subroutine weighted_block_edge_average_x_pre_downsampled_unmasked(fine, dx, coarse, nz)
    integer, intent(in) :: nz
    real, intent(in) :: fine(is:ie,js_coarse:je_coarse+1,1:nz)
    real, intent(in) :: dx(is:ie,js:je+1)
    real, intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse+1,1:nz)

    integer :: i, j, k, a, i_coarse, j_coarse

    a = coarsening_factor - 1
    do k = 1, nz
        do i = is, ie, coarsening_factor
          i_coarse = (i - 1) / coarsening_factor + 1
          do j = js, je + 1, coarsening_factor
              j_coarse = (j - 1) / coarsening_factor + 1
              coarse(i_coarse,j_coarse,k) = sum(dx(i:i+a,j) * fine(i:i+a,j_coarse,k)) / sum(dx(i:i+a,j))
          enddo
        enddo
    enddo
   end subroutine weighted_block_edge_average_x_pre_downsampled_unmasked

   subroutine weighted_block_edge_average_x_pre_downsampled_masked(fine, dx,&
    coarse, mask, nz)
    integer, intent(in) :: nz
    real, intent(in) :: fine(is:ie,js_coarse:je_coarse+1,1:nz)
    real, intent(in) :: dx(is:ie,js:je+1)
    logical, intent(in) :: mask(is:ie,js_coarse:je_coarse+1,1:nz)
    real, intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse+1,1:nz)

    real, allocatable :: weights(:,:), downsampled_dx(:,:)

    integer :: i, j, k, a, i_coarse, j_coarse

    allocate(weights(is:ie,js_coarse:je_coarse+1))
    allocate(downsampled_dx(is:ie,js_coarse:je_coarse+1))

    downsampled_dx = dx(:,js:je+1:coarsening_factor)

    a = coarsening_factor - 1
    do k = 1, nz
        where (mask(:,:,k))
          weights = downsampled_dx
        elsewhere
          weights = 0.0
        endwhere
        do i = is, ie, coarsening_factor
          i_coarse = (i - 1) / coarsening_factor + 1
          do j = js, je + 1, coarsening_factor
              j_coarse = (j - 1) / coarsening_factor + 1
              coarse(i_coarse,j_coarse,k) = sum(weights(i:i+a,j_coarse) * fine(i:i+a,j_coarse,k)) / sum(weights(i:i+a,j_coarse))
          enddo
        enddo
    enddo
   end subroutine weighted_block_edge_average_x_pre_downsampled_masked

   subroutine upsample_d_grid_x(field_in, field_out, nz)
    integer, intent(in) :: nz
    real, intent(in) :: field_in(is_coarse:ie_coarse,js_coarse:je_coarse+1,1:nz)
    real, intent(out) :: field_out(is:ie,js_coarse:je_coarse+1,1:nz)

    integer :: i, j, k, a, i_coarse
    a = coarsening_factor - 1
    do i = is, ie, coarsening_factor
       i_coarse = (i - 1) / coarsening_factor + 1
       do j = js_coarse, je_coarse + 1
          do k = 1, nz
             field_out(i:i+a,j,k) = field_in(i_coarse,j,k)
          enddo
       enddo
    enddo
   end subroutine upsample_d_grid_x

   subroutine remap_edges_along_x(field, phalf, dx, ptop, result)
    real, intent(in) :: field(is:ie,js:je+1,1:npz)
    real, intent(in) :: phalf(is-1,ie+1,js-1,je+1,1:npz+1)
    real, intent(in) :: dx(is:ie,js:je+1)
    real, intent(in) :: ptop
    real, intent(out) :: result(is_coarse:ie_coarse,js_coarse:je_coarse+1,1:npz)

    real, allocatable, dimension(:,:,:) :: phalf_d_grid, coarse_phalf_d_grid, coarse_phalf_d_grid_on_fine, remapped
    logical, allocatable :: mask(:,:,:)

    integer :: i, i_coarse, j, j_coarse, k, kn, km, kord, iv

    allocate(phalf_d_grid(is:ie,js_coarse:je_coarse+1,1:npz+1))
    allocate(coarse_phalf_d_grid(is_coarse:ie_coarse,js_coarse:je_coarse+1,1:npz+1))
    allocate(coarse_phalf_d_grid_on_fine(is:ie,js_coarse:je_coarse+1,1:npz+1))
    allocate(remapped(is:ie,js_coarse:je_coarse+1,1:npz))
    allocate(mask(is:ie,js_coarse:je_coarse+1,1:npz))

    ! Hard-code parameters related to mappm.
    kn = npz
    km = npz
    kord = 1
    iv = 1

    ! 1. Interpolate and downsample phalf
    call interpolate_to_d_grid_and_downsample_along_y(phalf, phalf_d_grid, npz+1)

    ! 2. Coarsen phalf on the D-grid
    call weighted_block_edge_average_x_pre_downsampled(phalf_d_grid, dx, coarse_phalf_d_grid, npz+1)

    ! 3. Upsample coarsened phalf back to the original resolution
    call upsample_d_grid_x(coarse_phalf_d_grid, coarse_phalf_d_grid_on_fine, npz+1)

    do j = js, je + 1, coarsening_factor
      j_coarse = (j - 1) / coarsening_factor + 1
      call mappm(km, phalf_d_grid(is:ie,j_coarse,:), field(is:ie,j,:), kn, &
        coarse_phalf_d_grid_on_fine(is:ie,j_coarse,:), &
        remapped(is:ie,j_coarse,:), is, ie, iv, kord, ptop)
    enddo

    ! 5. Create mask
    do k = 1, npz
      where (coarse_phalf_d_grid_on_fine(:,:,k+1) .lt. phalf_d_grid(:,:,npz+1))
        mask(:,:,k) = .true.
      elsewhere
        mask(:,:,k) = .false.
      endwhere
    enddo

    ! 6. Coarsen the remapped field
    call weighted_block_edge_average_x_pre_downsampled(remapped, dx, result, mask, npz)
  end subroutine remap_edges_along_x

  ! A naive routine for interpolating a field from the A-grid to the x-boundary
  ! of the D-grid; this is a specialized function that automatically
  ! downsamples to the coarse-grid on the downsampling dimension.
  subroutine interpolate_to_d_grid_and_downsample_along_x(field_in, field_out, nz)
    integer, intent(in) :: nz
    real, intent(in) :: field_in(is-1:ie+1,js-1:je+1,1:nz)
    real, intent(out) :: field_out(is_coarse:ie_coarse+1,js:je,1:nz)

    integer :: i, j, k, i_coarse

    do i = is,ie+1,coarsening_factor
       i_coarse = (i - 1) / coarsening_factor + 1
       do j = js,je
          do k = 1,nz
             field_out(i_coarse,j,k) = 0.5 * (field_in(i,j,k) + field_in(i-1,j,k))
          enddo
       enddo
    enddo
  end subroutine interpolate_to_d_grid_and_downsample_along_x

  subroutine weighted_block_edge_average_y_pre_downsampled_unmasked(fine, dy, coarse, nz)
    integer, intent(in) :: nz
    real, intent(in) :: fine(is_coarse:ie_coarse+1,js:je,1:nz)
    real, intent(in) :: dy(is:ie+1,js:je)
    real, intent(out) :: coarse(is_coarse:ie_coarse+1,js_coarse:je_coarse,1:nz)

    integer :: i, j, k, a, i_coarse, j_coarse

    a = coarsening_factor - 1
    do k = 1, nz
       do i = is, ie + 1, coarsening_factor
          i_coarse = (i - 1) / coarsening_factor + 1
          do j = js, je, coarsening_factor
             j_coarse = (j - 1) / coarsening_factor + 1
             coarse(i_coarse,j_coarse,k) = sum(dy(i,j:j+a) * fine(i_coarse,j:j+a,k)) / sum(dy(i,j:j+a))
          enddo
       enddo
    enddo
  end subroutine weighted_block_edge_average_y_pre_downsampled_unmasked

  subroutine weighted_block_edge_average_y_pre_downsampled_masked(fine, dy,&
    coarse, mask, nz)
    integer, intent(in) :: nz
    real, intent(in) :: fine(is_coarse:ie_coarse+1,js:je,1:nz)
    real, intent(in) :: dy(is:ie+1,js:je)
    logical, intent(in) :: mask(is_coarse:ie_coarse+1,js:je,1:nz)
    real, intent(out) :: coarse(is_coarse:ie_coarse+1,js_coarse:je_coarse,1:nz)

    real, allocatable :: weights(:,:), downsampled_dy(:,:)

    integer :: i, j, k, a, i_coarse, j_coarse


    allocate(weights(is_coarse:ie_coarse+1,js:je))
    allocate(downsampled_dy(is_coarse:ie_coarse+1,js:je))

    downsampled_dy = dy(is:ie+1:coarsening_factor,:)

    a = coarsening_factor - 1
    do k = 1, nz
        where (mask(:,:,k))
          weights = downsampled_dy
        elsewhere
          weights = 0.0
        endwhere
        do i = is, ie + 1, coarsening_factor
          i_coarse = (i - 1) / coarsening_factor + 1
          do j = js, je, coarsening_factor
              j_coarse = (j - 1) / coarsening_factor + 1
              coarse(i_coarse,j_coarse,k) = sum(weights(i_coarse,j:j+a) * fine(i_coarse,j:j+a,k)) / sum(weights(i_coarse,j:j+a))
          enddo
        enddo
    enddo
  end subroutine weighted_block_edge_average_y_pre_downsampled_masked

  subroutine upsample_d_grid_y(field_in, field_out, nz)
    integer, intent(in) :: nz
    real, intent(in) :: field_in(is_coarse:ie_coarse+1,js_coarse:je_coarse,1:nz)
    real, intent(out) :: field_out(is_coarse:ie_coarse+1,js:je,1:nz)

    integer :: i, j, k, a, j_coarse
    a = coarsening_factor - 1
    do i = is_coarse, ie_coarse + 1
       do j = js,je,coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          do k = 1, nz
             field_out(i,j:j+a,k) = field_in(i,j_coarse,k)
          enddo
       enddo
    enddo
  end subroutine upsample_d_grid_y

  subroutine remap_edges_along_y(field, phalf, dy, ptop, result)
    real, intent(in) :: field(is:ie+1,js:je,1:npz)
    real, intent(in) :: phalf(is-1:ie+1,js-1:je+1,1:npz+1)
    real, intent(in) :: dy(is:ie+1,js:je)
    real, intent(in) :: ptop
    real, intent(out) :: result(is_coarse:ie_coarse+1,js_coarse:je_coarse,1:npz)

    real, allocatable, dimension(:,:,:) :: phalf_d_grid, coarse_phalf_d_grid, coarse_phalf_d_grid_on_fine, remapped
    logical, allocatable :: mask(:,:,:)

    integer :: i, i_coarse, j, j_coarse, k, kn, km, kord, iv

    allocate(phalf_d_grid(is_coarse:ie_coarse+1,js:je,1:npz+1))
    allocate(coarse_phalf_d_grid(is_coarse:ie_coarse+1,js_coarse:je_coarse,1:npz+1))
    allocate(coarse_phalf_d_grid_on_fine(is_coarse:ie_coarse+1,js:je,1:npz+1))
    allocate(remapped(is_coarse:ie_coarse+1,js:je,1:npz))
    allocate(mask(is_coarse:ie_coarse+1,js:je,1:npz))

    ! Hard-code parameters related to mappm.
    kn = npz
    km = npz
    kord = 1
    iv = 1

    ! 1. Interpolate and downsample phalf
    call interpolate_to_d_grid_and_downsample_along_x(phalf, phalf_d_grid, npz+1)

    ! 2. Coarsen phalf on the D-grid
    call weighted_block_edge_average_y_pre_downsampled(phalf_d_grid, dy, coarse_phalf_d_grid, npz+1)

    ! 3. Upsample coarsened phalf back to the original resolution
    call upsample_d_grid_y(coarse_phalf_d_grid, coarse_phalf_d_grid_on_fine, npz+1)

    do i = is, ie + 1, coarsening_factor
      i_coarse = (i - 1) / coarsening_factor + 1
      call mappm(km, phalf_d_grid(i_coarse,js:je,:), field(i,js:je,:), kn, &
        coarse_phalf_d_grid_on_fine(i_coarse,js:je,:), &
        remapped(i_coarse,js:je,:), js, je, iv, kord, ptop)
    enddo

    ! 5. Create mask
    do k = 1, npz
      where (coarse_phalf_d_grid_on_fine(:,:,k+1) .lt. phalf_d_grid(:,:,npz+1))
        mask(:,:,k) = .true.
      elsewhere
        mask(:,:,k) = .false.
      endwhere
    enddo

    ! 6. Coarsen the remapped field
    call weighted_block_edge_average_y_pre_downsampled(remapped, dy, result, mask, npz)
  end subroutine remap_edges_along_y

  subroutine block_edge_sum_x_2d_full_input(fine, coarse)
    real, intent(in) :: fine(is:ie,js:je+1)
    real, intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse+1)

    integer :: i, j, i_coarse, j_coarse, offset

    offset = coarsening_factor - 1
    do i = is, ie, coarsening_factor
       i_coarse = (i - 1) / coarsening_factor + 1
       do j = js, je + 1, coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          coarse(i_coarse,j_coarse) = sum(fine(i:i+offset,j))
       enddo
    enddo
  end subroutine block_edge_sum_x_2d_full_input

  subroutine block_edge_sum_y_2d_full_input(fine, coarse)
    real, intent(in) :: fine(is:ie+1,js:je)
    real, intent(out) :: coarse(is_coarse:ie_coarse+1,js_coarse:je_coarse)

    integer :: i, j, i_coarse, j_coarse, offset

    offset = coarsening_factor - 1
    do i = is, ie + 1, coarsening_factor
       i_coarse = (i - 1) / coarsening_factor + 1
       do j = js, je, coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          coarse(i_coarse,j_coarse) = sum(fine(i,j:j+offset))
       enddo
    enddo
  end subroutine block_edge_sum_y_2d_full_input

  ! Needed for computing variances over coarse grid cells.
  subroutine anomaly_2d(weights, fine, anom)
    real, intent(in) :: weights(is:ie,js:je), fine(is:ie,js:je)
    real, intent(out) :: anom(is:ie,js:je)

    integer :: i, j, i_coarse, j_coarse, offset
    real, allocatable :: coarse(:,:)

    allocate(coarse(is_coarse:ie_coarse,js_coarse:je_coarse))

    ! First compute the coarse-grained field
    call weighted_block_average(weights, fine, coarse)

    ! Then subtract it off
    offset = coarsening_factor - 1
    do i = is, ie, coarsening_factor
       i_coarse = (i - 1) / coarsening_factor + 1
       do j = js, je, coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          anom(i:i+offset,j:j+offset) = fine(i:i+offset,j:j+offset) - coarse(i_coarse,j_coarse)
       enddo
    enddo
  end subroutine anomaly_2d

  subroutine anomaly_3d_weights_3d_array(weights, fine, anom)
    real, intent(in) :: weights(is:ie,js:je,1:npz), fine(is:ie,js:je,1:npz)
    real, intent(out) :: anom(is:ie,js:je,1:npz)

    integer :: k
    do k = 1,npz
       call anomaly_2d(weights(is:ie,js:je,k), fine(is:ie,js:je,k), anom(is:ie,js:je,k))
    enddo
  end subroutine anomaly_3d_weights_3d_array

  subroutine anomaly_2d_weights_3d_array(weights, fine, anom)
    real, intent(in) :: weights(is:ie,js:je)
    real, intent(in) :: fine(is:ie,js:je,1:npz)
    real, intent(out) :: anom(is:ie,js:je,1:npz)

    integer :: k
    do k = 1,npz
       call anomaly_2d(weights(is:ie,js:je), fine(is:ie,js:je,k), anom(is:ie,js:je,k))
    enddo
  end subroutine anomaly_2d_weights_3d_array

  subroutine eddy_covariance_2d_weights(weights, field_a, field_b, coarse)
    real, intent(in) :: weights(is:ie,js:je)
    real, intent(in) :: field_a(is:ie,js:je,1:npz)
    real, intent(in) :: field_b(is:ie,js:je,1:npz)
    real, intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz)

    real, allocatable :: anom_a(:,:,:)
    real, allocatable :: anom_b(:,:,:)

    allocate(anom_a(is:ie,js:je,1:npz))
    allocate(anom_b(is:ie,js:je,1:npz))

    call anomaly_2d_weights_3d_array(weights, field_a, anom_a)
    call anomaly_2d_weights_3d_array(weights, field_b, anom_b)
    call weighted_block_average(weights, anom_a * anom_b, coarse)
  end subroutine eddy_covariance_2d_weights

  subroutine eddy_covariance_3d_weights(weights, field_a, field_b, coarse)
    real, intent(in) :: weights(is:ie,js:je,1:npz)
    real, intent(in) :: field_a(is:ie,js:je,1:npz)
    real, intent(in) :: field_b(is:ie,js:je,1:npz)
    real, intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz)

    real, allocatable :: anom_a(:,:,:)
    real, allocatable :: anom_b(:,:,:)

    allocate(anom_a(is:ie,js:je,1:npz))
    allocate(anom_b(is:ie,js:je,1:npz))

    call anomaly_3d_weights_3d_array(weights, field_a, anom_a)
    call anomaly_3d_weights_3d_array(weights, field_b, anom_b)
    call weighted_block_average(weights, anom_a * anom_b, coarse)
  end subroutine eddy_covariance_3d_weights

! Port mappm for single and double precision
 subroutine mappm_real4(km, pe1, q1, kn, pe2, q2, i1, i2, iv, kord, ptop)

! IV = 0: constituents
! IV = 1: potential temp
! IV =-1: winds

! Mass flux preserving mapping: q1(im,km) -> q2(im,kn)

! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate

 integer, intent(in):: i1, i2, km, kn, kord, iv
 real(kind=4), intent(in ):: pe1(i1:i2,km+1), pe2(i1:i2,kn+1)
 real(kind=4), intent(in )::  q1(i1:i2,km)
 real(kind=4), intent(out)::  q2(i1:i2,kn)
 real(kind=4), intent(IN) :: ptop
! local
      real(kind=4)  qs(i1:i2)
      real(kind=4) dp1(i1:i2,km)
      real(kind=4) a4(4,i1:i2,km)
      integer i, k, l
      integer k0, k1
      real(kind=4) pl, pr, tt, delp, qsum, dpsum, esl

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
            a4(1,i,k) = q1(i,k)
         enddo
      enddo

      if ( kord >7 ) then
           call  cs_profile( qs, a4, dp1, km, i1, i2, iv, kord )
      else
           call ppm_profile( a4, dp1, km, i1, i2, iv, kord )
      endif

!------------------------------------
! Lowest layer: constant distribution
!------------------------------------
#ifdef NGGPS_SUBMITTED
      do i=i1,i2
         a4(2,i,km) = q1(i,km)
         a4(3,i,km) = q1(i,km)
         a4(4,i,km) = 0.
      enddo
#endif

      do 5555 i=i1,i2
         k0 = 1
      do 555 k=1,kn

         if(pe2(i,k) .le. pe1(i,1)) then
! above old ptop
            q2(i,k) = q1(i,1)
         elseif(pe2(i,k) .ge. pe1(i,km+1)) then
! Entire grid below old ps
#ifdef NGGPS_SUBMITTED
            q2(i,k) = a4(3,i,km)   ! this is not good.
#else
            q2(i,k) = q1(i,km)
#endif
         else

         do 45 L=k0,km
! locate the top edge at pe2(i,k)
         if( pe2(i,k) .ge. pe1(i,L) .and.        &
             pe2(i,k) .le. pe1(i,L+1)    ) then
             k0 = L
             PL = (pe2(i,k)-pe1(i,L)) / dp1(i,L)
             if(pe2(i,k+1) .le. pe1(i,L+1)) then

! entire new grid is within the original grid
               PR = (pe2(i,k+1)-pe1(i,L)) / dp1(i,L)
               TT = r3_real4*(PR*(PR+PL)+PL**2)
               q2(i,k) = a4(2,i,L) + 0.5*(a4(4,i,L)+a4(3,i,L)  &
                       - a4(2,i,L))*(PR+PL) - a4(4,i,L)*TT
              goto 555
             else
! Fractional area...
              delp = pe1(i,L+1) - pe2(i,k)
              TT   = r3_real4*(1.+PL*(1.+PL))
              qsum = delp*(a4(2,i,L)+0.5*(a4(4,i,L)+            &
                     a4(3,i,L)-a4(2,i,L))*(1.+PL)-a4(4,i,L)*TT)
              dpsum = delp
              k1 = L + 1
             goto 111
             endif
         endif
45       continue

111      continue
         do 55 L=k1,km
         if( pe2(i,k+1) .gt. pe1(i,L+1) ) then

! Whole layer..

            qsum  =  qsum + dp1(i,L)*q1(i,L)
            dpsum = dpsum + dp1(i,L)
         else
           delp = pe2(i,k+1)-pe1(i,L)
           esl  = delp / dp1(i,L)
           qsum = qsum + delp * (a4(2,i,L)+0.5*esl*            &
                 (a4(3,i,L)-a4(2,i,L)+a4(4,i,L)*(1.-r23_real4*esl)) )
          dpsum = dpsum + delp
           k0 = L
           goto 123
         endif
55       continue
        delp = pe2(i,k+1) - pe1(i,km+1)
        if(delp > 0.) then
! Extended below old ps
#ifdef NGGPS_SUBMITTED
           qsum = qsum + delp * a4(3,i,km)    ! not good.
#else
           qsum = qsum + delp * q1(i,km)
#endif
          dpsum = dpsum + delp
        endif
123     q2(i,k) = qsum / dpsum
      endif
555   continue
5555  continue

 end subroutine mappm_real4

 subroutine mappm_real8(km, pe1, q1, kn, pe2, q2, i1, i2, iv, kord, ptop)

! IV = 0: constituents
! IV = 1: potential temp
! IV =-1: winds

! Mass flux preserving mapping: q1(im,km) -> q2(im,kn)

! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate

 integer, intent(in):: i1, i2, km, kn, kord, iv
 real(kind=8), intent(in ):: pe1(i1:i2,km+1), pe2(i1:i2,kn+1)
 real(kind=8), intent(in )::  q1(i1:i2,km)
 real(kind=8), intent(out)::  q2(i1:i2,kn)
 real(kind=8), intent(IN) :: ptop
! local
      real(kind=8)  qs(i1:i2)
      real(kind=8) dp1(i1:i2,km)
      real(kind=8) a4(4,i1:i2,km)
      integer i, k, l
      integer k0, k1
      real(kind=8) pl, pr, tt, delp, qsum, dpsum, esl

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
            a4(1,i,k) = q1(i,k)
         enddo
      enddo

      if ( kord >7 ) then
           call  cs_profile( qs, a4, dp1, km, i1, i2, iv, kord )
      else
           call ppm_profile( a4, dp1, km, i1, i2, iv, kord )
      endif

!------------------------------------
! Lowest layer: constant distribution
!------------------------------------
#ifdef NGGPS_SUBMITTED
      do i=i1,i2
         a4(2,i,km) = q1(i,km)
         a4(3,i,km) = q1(i,km)
         a4(4,i,km) = 0.
      enddo
#endif

      do 5555 i=i1,i2
         k0 = 1
      do 555 k=1,kn

         if(pe2(i,k) .le. pe1(i,1)) then
! above old ptop
            q2(i,k) = q1(i,1)
         elseif(pe2(i,k) .ge. pe1(i,km+1)) then
! Entire grid below old ps
#ifdef NGGPS_SUBMITTED
            q2(i,k) = a4(3,i,km)   ! this is not good.
#else
            q2(i,k) = q1(i,km)
#endif
         else

         do 45 L=k0,km
! locate the top edge at pe2(i,k)
         if( pe2(i,k) .ge. pe1(i,L) .and.        &
             pe2(i,k) .le. pe1(i,L+1)    ) then
             k0 = L
             PL = (pe2(i,k)-pe1(i,L)) / dp1(i,L)
             if(pe2(i,k+1) .le. pe1(i,L+1)) then

! entire new grid is within the original grid
               PR = (pe2(i,k+1)-pe1(i,L)) / dp1(i,L)
               TT = r3_real8*(PR*(PR+PL)+PL**2)
               q2(i,k) = a4(2,i,L) + 0.5*(a4(4,i,L)+a4(3,i,L)  &
                       - a4(2,i,L))*(PR+PL) - a4(4,i,L)*TT
              goto 555
             else
! Fractional area...
              delp = pe1(i,L+1) - pe2(i,k)
              TT   = r3_real8*(1.+PL*(1.+PL))
              qsum = delp*(a4(2,i,L)+0.5*(a4(4,i,L)+            &
                     a4(3,i,L)-a4(2,i,L))*(1.+PL)-a4(4,i,L)*TT)
              dpsum = delp
              k1 = L + 1
             goto 111
             endif
         endif
45       continue

111      continue
         do 55 L=k1,km
         if( pe2(i,k+1) .gt. pe1(i,L+1) ) then

! Whole layer..

            qsum  =  qsum + dp1(i,L)*q1(i,L)
            dpsum = dpsum + dp1(i,L)
         else
           delp = pe2(i,k+1)-pe1(i,L)
           esl  = delp / dp1(i,L)
           qsum = qsum + delp * (a4(2,i,L)+0.5*esl*            &
                 (a4(3,i,L)-a4(2,i,L)+a4(4,i,L)*(1.-r23_real8*esl)) )
          dpsum = dpsum + delp
           k0 = L
           goto 123
         endif
55       continue
        delp = pe2(i,k+1) - pe1(i,km+1)
        if(delp > 0.) then
! Extended below old ps
#ifdef NGGPS_SUBMITTED
           qsum = qsum + delp * a4(3,i,km)    ! not good.
#else
           qsum = qsum + delp * q1(i,km)
#endif
          dpsum = dpsum + delp
        endif
123     q2(i,k) = qsum / dpsum
      endif
555   continue
5555  continue

 end subroutine mappm_real8

 subroutine cs_profile_real4(qs, a4, delp, km, i1, i2, iv, kord)
! Optimized vertical profile reconstruction:
! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
 integer, intent(in):: i1, i2
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: iv      ! iv =-2: vertical velocity
                               ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
 integer, intent(in):: kord
 real(kind=4), intent(in)   ::   qs(i1:i2)
 real(kind=4), intent(in)   :: delp(i1:i2,km)     ! layer pressure thickness
 real(kind=4), intent(inout):: a4(4,i1:i2,km)     ! Interpolated values
!-----------------------------------------------------------------------
 logical, dimension(i1:i2,km):: extm, ext5, ext6
 real(kind=4)  gam(i1:i2,km)
 real(kind=4)    q(i1:i2,km+1) ! interface values
 real(kind=4)   d4(i1:i2)
 real(kind=4)   bet, a_bot, grat
 real(kind=4)   pmp_1, lac_1, pmp_2, lac_2, x0, x1
 integer i, k, im

 if ( iv .eq. -2 ) then
      do i=i1,i2
         gam(i,2) = 0.5
           q(i,1) = 1.5*a4(1,i,1)
      enddo
      do k=2,km-1
         do i=i1, i2
                  grat = delp(i,k-1) / delp(i,k)
                   bet =  2. + grat + grat - gam(i,k)
                q(i,k) = (3.*(a4(1,i,k-1)+a4(1,i,k)) - q(i,k-1))/bet
            gam(i,k+1) = grat / bet
         enddo
      enddo
      do i=i1,i2
            grat = delp(i,km-1) / delp(i,km)
         q(i,km) = (3.*(a4(1,i,km-1)+a4(1,i,km)) - grat*qs(i) - q(i,km-1)) /  &
                   (2. + grat + grat - gam(i,km))
         q(i,km+1) = qs(i)
      enddo
      do k=km-1,1,-1
        do i=i1,i2
           q(i,k) = q(i,k) - gam(i,k+1)*q(i,k+1)
        enddo
      enddo

else ! all others
  do i=i1,i2
         grat = delp(i,2) / delp(i,1)   ! grid ratio
          bet = grat*(grat+0.5)
       q(i,1) = ( (grat+grat)*(grat+1.)*a4(1,i,1) + a4(1,i,2) ) / bet
     gam(i,1) = ( 1. + grat*(grat+1.5) ) / bet
  enddo

  if (iv.eq.-3) then !LBC for vertical velocities
    do k=2,km-1
       do i=i1,i2
             d4(i) = delp(i,k-1) / delp(i,k)
               bet =  2. + d4(i) + d4(i) - gam(i,k-1)
            q(i,k) = ( 3.*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - q(i,k-1) )/bet
          gam(i,k) = d4(i) / bet
       enddo
    enddo

    do i=i1,i2
       !    a_bot = 1. + d4(i)*(d4(i)+1.5)
       !q(i,km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1,i,km)+a4(1,i,km-1)-a_bot*q(i,km))  &
       !          / ( d4(i)*(d4(i)+0.5) - a_bot*gam(i,km) )
       d4(i) = delp(i,km-1) / delp(i,km)
       bet =  2. + d4(i) + d4(i) - gam(i,km-1)
       grat = delp(i,km-1) / delp(i,km)
       q(i,km) = ( 3.*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - grat*qs(i) - q(i,k-1) )/bet
       q(i,km+1) = qs(i)
    enddo

  else ! all others
    do k=2,km
       do i=i1,i2
             d4(i) = delp(i,k-1) / delp(i,k)
               bet =  2. + d4(i) + d4(i) - gam(i,k-1)
            q(i,k) = ( 3.*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - q(i,k-1) )/bet
          gam(i,k) = d4(i) / bet
       enddo
    enddo

    do i=i1,i2
           a_bot = 1. + d4(i)*(d4(i)+1.5)
       q(i,km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1,i,km)+a4(1,i,km-1)-a_bot*q(i,km))  &
                 / ( d4(i)*(d4(i)+0.5) - a_bot*gam(i,km) )
    enddo
  endif

  do k=km,1,-1
     do i=i1,i2
        q(i,k) = q(i,k) - gam(i,k)*q(i,k+1)
     enddo
  enddo
 endif
!----- Perfectly linear scheme --------------------------------
 if ( abs(kord) > 16 ) then
  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
  enddo
  return
 endif
!----- Perfectly linear scheme --------------------------------

!------------------
! Apply constraints
!------------------
  im = i2 - i1 + 1

! Apply *large-scale* constraints
  do i=i1,i2
     q(i,2) = min( q(i,2), max(a4(1,i,1), a4(1,i,2)) )
     q(i,2) = max( q(i,2), min(a4(1,i,1), a4(1,i,2)) )
  enddo

  do k=2,km
     do i=i1,i2
        gam(i,k) = a4(1,i,k) - a4(1,i,k-1) ! now dq
     enddo
  enddo

! Interior:
  do k=3,km-1
     do i=i1,i2
        if ( gam(i,k-1)*gam(i,k+1)>0. ) then
! Apply large-scale constraint to ALL fields if not local max/min
             q(i,k) = min( q(i,k), max(a4(1,i,k-1),a4(1,i,k)) )
             q(i,k) = max( q(i,k), min(a4(1,i,k-1),a4(1,i,k)) )
        else
          if ( gam(i,k-1) > 0. ) then
! There exists a local max
               q(i,k) = max(q(i,k), min(a4(1,i,k-1),a4(1,i,k)))
          else
! There exists a local min
                 q(i,k) = min(q(i,k), max(a4(1,i,k-1),a4(1,i,k)))
               if ( iv==0 ) q(i,k) = max(0., q(i,k)) ! positive-definite
          endif
        endif
     enddo
  enddo

! Bottom:
  do i=i1,i2
     q(i,km) = min( q(i,km), max(a4(1,i,km-1), a4(1,i,km)) )
     q(i,km) = max( q(i,km), min(a4(1,i,km-1), a4(1,i,km)) )
  enddo

  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
     enddo
  enddo

  do k=1,km
     if ( k==1 .or. k==km ) then
       do i=i1,i2
          extm(i,k) = (a4(2,i,k)-a4(1,i,k)) * (a4(3,i,k)-a4(1,i,k)) > 0.
       enddo
     else
       do i=i1,i2
          extm(i,k) = gam(i,k)*gam(i,k+1) < 0.
       enddo
     endif
     if ( abs(kord) > 9 ) then
       do i=i1,i2
          x0 = 2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k))
          x1 = abs(a4(2,i,k)-a4(3,i,k))
          a4(4,i,k) = 3.*x0
          ext5(i,k) = abs(x0) > x1
          ext6(i,k) = abs(a4(4,i,k)) > x1
       enddo
     endif
  enddo

!---------------------------
! Apply subgrid constraints:
!---------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
! Top 2 and bottom 2 layers always use monotonic mapping

  if ( iv==0 ) then
     do i=i1,i2
        a4(2,i,1) = max(0., a4(2,i,1))
     enddo
  elseif ( iv==-1 ) then
      do i=i1,i2
         if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
      enddo
  elseif ( iv==2 ) then
     do i=i1,i2
        a4(2,i,1) = a4(1,i,1)
        a4(3,i,1) = a4(1,i,1)
        a4(4,i,1) = 0.
     enddo
  endif

  if ( iv/=2 ) then
     do i=i1,i2
        a4(4,i,1) = 3.*(2.*a4(1,i,1) - (a4(2,i,1)+a4(3,i,1)))
     enddo
     call cs_limiters_real4(im, extm(i1,1), a4(1,i1,1), 1)
  endif

! k=2
   do i=i1,i2
      a4(4,i,2) = 3.*(2.*a4(1,i,2) - (a4(2,i,2)+a4(3,i,2)))
   enddo
   call cs_limiters_real4(im, extm(i1,2), a4(1,i1,2), 2)

!-------------------------------------
! Huynh's 2nd constraint for interior:
!-------------------------------------
  do k=3,km-2
     if ( abs(kord)<9 ) then
       do i=i1,i2
! Left  edges
          pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
          lac_1 = pmp_1 + 1.5*gam(i,k+2)
          a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                         max(a4(1,i,k), pmp_1, lac_1) )
! Right edges
          pmp_2 = a4(1,i,k) + 2.*gam(i,k)
          lac_2 = pmp_2 - 1.5*gam(i,k-1)
          a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                         max(a4(1,i,k), pmp_2, lac_2) )

          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo

     elseif ( abs(kord)==9 ) then
       do i=i1,i2
          if ( extm(i,k) .and. extm(i,k-1) ) then  ! c90_mp122
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else if ( extm(i,k) .and. extm(i,k+1) ) then  ! c90_mp122
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==10 ) then
       do i=i1,i2
          if( ext5(i,k) ) then
              if( ext5(i,k-1) .or. ext5(i,k+1) ) then
                   a4(2,i,k) = a4(1,i,k)
                   a4(3,i,k) = a4(1,i,k)
              elseif ( ext6(i,k-1) .or. ext6(i,k+1) ) then
                   pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                   lac_1 = pmp_1 + 1.5*gam(i,k+2)
                   a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                                  max(a4(1,i,k), pmp_1, lac_1) )
                   pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                   lac_2 = pmp_2 - 1.5*gam(i,k-1)
                   a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                                  max(a4(1,i,k), pmp_2, lac_2) )
              endif
          elseif( ext6(i,k) ) then
              if( ext5(i,k-1) .or. ext5(i,k+1) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
                  a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                                 max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
                  a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                                 max(a4(1,i,k), pmp_2, lac_2) )
              endif
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     elseif ( abs(kord)==12 ) then
       do i=i1,i2
          if( extm(i,k) ) then
! grid-scale 2-delta-z wave detected
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
          else        ! not a local extremum
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==13 ) then
       do i=i1,i2
          if( ext6(i,k) ) then
             if ( ext6(i,k-1) .and. ext6(i,k+1) ) then
! grid-scale 2-delta-z wave detected
                 a4(2,i,k) = a4(1,i,k)
                 a4(3,i,k) = a4(1,i,k)
             endif
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     elseif ( abs(kord)==14 ) then

       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo

     elseif ( abs(kord)==15 ) then   ! revised kord=9 scehem
       do i=i1,i2
          if ( ext5(i,k) ) then  ! c90_mp122
             if ( ext5(i,k-1) .or. ext5(i,k+1) ) then  ! c90_mp122
! grid-scale 2-delta-z wave detected
                  a4(2,i,k) = a4(1,i,k)
                  a4(3,i,k) = a4(1,i,k)
             endif
          elseif( ext6(i,k) ) then
! Check within the smooth region if subgrid profile is non-monotonic
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     elseif ( abs(kord)==16 ) then
       do i=i1,i2
          if( ext5(i,k) ) then
             if ( ext5(i,k-1) .or. ext5(i,k+1) ) then
                 a4(2,i,k) = a4(1,i,k)
                 a4(3,i,k) = a4(1,i,k)
             elseif ( ext6(i,k-1) .or. ext6(i,k+1) ) then
                 ! Left  edges
                 pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                 lac_1 = pmp_1 + 1.5*gam(i,k+2)
                 a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                     max(a4(1,i,k), pmp_1, lac_1) )
                 ! Right edges
                 pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                 lac_2 = pmp_2 - 1.5*gam(i,k-1)
                 a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                     max(a4(1,i,k), pmp_2, lac_2) )
             endif
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     else      ! kord = 11
       do i=i1,i2
         if ( ext5(i,k) .and. (ext5(i,k-1) .or. ext5(i,k+1)) ) then
! Noisy region:
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
         else
              a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         endif
       enddo
     endif

! Additional constraint to ensure positivity
     if ( iv==0 ) call cs_limiters_real4(im, extm(i1,k), a4(1,i1,k), 0)

  enddo      ! k-loop

!----------------------------------
! Bottom layer subgrid constraints:
!----------------------------------
  if ( iv==0 ) then
     do i=i1,i2
        a4(3,i,km) = max(0., a4(3,i,km))
     enddo
  elseif ( iv .eq. -1 ) then
      do i=i1,i2
         if ( a4(3,i,km)*a4(1,i,km) <= 0. )  a4(3,i,km) = 0.
      enddo
  endif

  do k=km-1,km
     do i=i1,i2
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
     if(k==(km-1)) call cs_limiters_real4(im, extm(i1,k), a4(1,i1,k), 2)
     if(k== km   ) call cs_limiters_real4(im, extm(i1,k), a4(1,i1,k), 1)
  enddo

 end subroutine cs_profile_real4

 subroutine cs_limiters_real4(im, extm, a4, iv)
 integer, intent(in) :: im
 integer, intent(in) :: iv
 logical, intent(in) :: extm(im)
 real(kind=4) , intent(inout) :: a4(4,im)   ! PPM array
! !LOCAL VARIABLES:
 real(kind=4)  da1, da2, a6da
 integer i

 if ( iv==0 ) then
! Positive definite constraint
    do i=1,im
    if( a4(1,i)<=0.) then
        a4(2,i) = a4(1,i)
        a4(3,i) = a4(1,i)
        a4(4,i) = 0.
    else
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
         if( (a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12_real4) < 0. ) then
! local minimum is negative
             if( a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i) ) then
                 a4(3,i) = a4(1,i)
                 a4(2,i) = a4(1,i)
                 a4(4,i) = 0.
             elseif( a4(3,i) > a4(2,i) ) then
                 a4(4,i) = 3.*(a4(2,i)-a4(1,i))
                 a4(3,i) = a4(2,i) - a4(4,i)
             else
                 a4(4,i) = 3.*(a4(3,i)-a4(1,i))
                 a4(2,i) = a4(3,i) - a4(4,i)
             endif
         endif
      endif
    endif
    enddo
 elseif ( iv==1 ) then
    do i=1,im
      if( (a4(1,i)-a4(2,i))*(a4(1,i)-a4(3,i))>=0. ) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
    enddo
 else
! Standard PPM constraint
    do i=1,im
      if( extm(i) ) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
    enddo
 endif
 end subroutine cs_limiters_real4

 subroutine cs_profile_real8(qs, a4, delp, km, i1, i2, iv, kord)
! Optimized vertical profile reconstruction:
! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
 integer, intent(in):: i1, i2
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: iv      ! iv =-2: vertical velocity
                               ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
 integer, intent(in):: kord
 real(kind=8), intent(in)   ::   qs(i1:i2)
 real(kind=8), intent(in)   :: delp(i1:i2,km)     ! layer pressure thickness
 real(kind=8), intent(inout):: a4(4,i1:i2,km)     ! Interpolated values
!-----------------------------------------------------------------------
 logical, dimension(i1:i2,km):: extm, ext5, ext6
 real(kind=8)  gam(i1:i2,km)
 real(kind=8)    q(i1:i2,km+1) ! interface values
 real(kind=8)   d4(i1:i2)
 real(kind=8)   bet, a_bot, grat
 real(kind=8)   pmp_1, lac_1, pmp_2, lac_2, x0, x1
 integer i, k, im

 if ( iv .eq. -2 ) then
      do i=i1,i2
         gam(i,2) = 0.5
           q(i,1) = 1.5*a4(1,i,1)
      enddo
      do k=2,km-1
         do i=i1, i2
                  grat = delp(i,k-1) / delp(i,k)
                   bet =  2. + grat + grat - gam(i,k)
                q(i,k) = (3.*(a4(1,i,k-1)+a4(1,i,k)) - q(i,k-1))/bet
            gam(i,k+1) = grat / bet
         enddo
      enddo
      do i=i1,i2
            grat = delp(i,km-1) / delp(i,km)
         q(i,km) = (3.*(a4(1,i,km-1)+a4(1,i,km)) - grat*qs(i) - q(i,km-1)) /  &
                   (2. + grat + grat - gam(i,km))
         q(i,km+1) = qs(i)
      enddo
      do k=km-1,1,-1
        do i=i1,i2
           q(i,k) = q(i,k) - gam(i,k+1)*q(i,k+1)
        enddo
      enddo

else ! all others
  do i=i1,i2
         grat = delp(i,2) / delp(i,1)   ! grid ratio
          bet = grat*(grat+0.5)
       q(i,1) = ( (grat+grat)*(grat+1.)*a4(1,i,1) + a4(1,i,2) ) / bet
     gam(i,1) = ( 1. + grat*(grat+1.5) ) / bet
  enddo

  if (iv.eq.-3) then !LBC for vertical velocities
    do k=2,km-1
       do i=i1,i2
             d4(i) = delp(i,k-1) / delp(i,k)
               bet =  2. + d4(i) + d4(i) - gam(i,k-1)
            q(i,k) = ( 3.*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - q(i,k-1) )/bet
          gam(i,k) = d4(i) / bet
       enddo
    enddo

    do i=i1,i2
       !    a_bot = 1. + d4(i)*(d4(i)+1.5)
       !q(i,km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1,i,km)+a4(1,i,km-1)-a_bot*q(i,km))  &
       !          / ( d4(i)*(d4(i)+0.5) - a_bot*gam(i,km) )
       d4(i) = delp(i,km-1) / delp(i,km)
       bet =  2. + d4(i) + d4(i) - gam(i,km-1)
       grat = delp(i,km-1) / delp(i,km)
       q(i,km) = ( 3.*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - grat*qs(i) - q(i,k-1) )/bet
       q(i,km+1) = qs(i)
    enddo

  else ! all others
    do k=2,km
       do i=i1,i2
             d4(i) = delp(i,k-1) / delp(i,k)
               bet =  2. + d4(i) + d4(i) - gam(i,k-1)
            q(i,k) = ( 3.*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - q(i,k-1) )/bet
          gam(i,k) = d4(i) / bet
       enddo
    enddo

    do i=i1,i2
           a_bot = 1. + d4(i)*(d4(i)+1.5)
       q(i,km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1,i,km)+a4(1,i,km-1)-a_bot*q(i,km))  &
                 / ( d4(i)*(d4(i)+0.5) - a_bot*gam(i,km) )
    enddo
  endif

  do k=km,1,-1
     do i=i1,i2
        q(i,k) = q(i,k) - gam(i,k)*q(i,k+1)
     enddo
  enddo
 endif
!----- Perfectly linear scheme --------------------------------
 if ( abs(kord) > 16 ) then
  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
  enddo
  return
 endif
!----- Perfectly linear scheme --------------------------------

!------------------
! Apply constraints
!------------------
  im = i2 - i1 + 1

! Apply *large-scale* constraints
  do i=i1,i2
     q(i,2) = min( q(i,2), max(a4(1,i,1), a4(1,i,2)) )
     q(i,2) = max( q(i,2), min(a4(1,i,1), a4(1,i,2)) )
  enddo

  do k=2,km
     do i=i1,i2
        gam(i,k) = a4(1,i,k) - a4(1,i,k-1) ! now dq
     enddo
  enddo

! Interior:
  do k=3,km-1
     do i=i1,i2
        if ( gam(i,k-1)*gam(i,k+1)>0. ) then
! Apply large-scale constraint to ALL fields if not local max/min
             q(i,k) = min( q(i,k), max(a4(1,i,k-1),a4(1,i,k)) )
             q(i,k) = max( q(i,k), min(a4(1,i,k-1),a4(1,i,k)) )
        else
          if ( gam(i,k-1) > 0. ) then
! There exists a local max
               q(i,k) = max(q(i,k), min(a4(1,i,k-1),a4(1,i,k)))
          else
! There exists a local min
                 q(i,k) = min(q(i,k), max(a4(1,i,k-1),a4(1,i,k)))
               if ( iv==0 ) q(i,k) = max(0., q(i,k)) ! positive-definite
          endif
        endif
     enddo
  enddo

! Bottom:
  do i=i1,i2
     q(i,km) = min( q(i,km), max(a4(1,i,km-1), a4(1,i,km)) )
     q(i,km) = max( q(i,km), min(a4(1,i,km-1), a4(1,i,km)) )
  enddo

  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
     enddo
  enddo

  do k=1,km
     if ( k==1 .or. k==km ) then
       do i=i1,i2
          extm(i,k) = (a4(2,i,k)-a4(1,i,k)) * (a4(3,i,k)-a4(1,i,k)) > 0.
       enddo
     else
       do i=i1,i2
          extm(i,k) = gam(i,k)*gam(i,k+1) < 0.
       enddo
     endif
     if ( abs(kord) > 9 ) then
       do i=i1,i2
          x0 = 2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k))
          x1 = abs(a4(2,i,k)-a4(3,i,k))
          a4(4,i,k) = 3.*x0
          ext5(i,k) = abs(x0) > x1
          ext6(i,k) = abs(a4(4,i,k)) > x1
       enddo
     endif
  enddo

!---------------------------
! Apply subgrid constraints:
!---------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
! Top 2 and bottom 2 layers always use monotonic mapping

  if ( iv==0 ) then
     do i=i1,i2
        a4(2,i,1) = max(0., a4(2,i,1))
     enddo
  elseif ( iv==-1 ) then
      do i=i1,i2
         if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
      enddo
  elseif ( iv==2 ) then
     do i=i1,i2
        a4(2,i,1) = a4(1,i,1)
        a4(3,i,1) = a4(1,i,1)
        a4(4,i,1) = 0.
     enddo
  endif

  if ( iv/=2 ) then
     do i=i1,i2
        a4(4,i,1) = 3.*(2.*a4(1,i,1) - (a4(2,i,1)+a4(3,i,1)))
     enddo
     call cs_limiters_real8(im, extm(i1,1), a4(1,i1,1), 1)
  endif

! k=2
   do i=i1,i2
      a4(4,i,2) = 3.*(2.*a4(1,i,2) - (a4(2,i,2)+a4(3,i,2)))
   enddo
   call cs_limiters_real8(im, extm(i1,2), a4(1,i1,2), 2)

!-------------------------------------
! Huynh's 2nd constraint for interior:
!-------------------------------------
  do k=3,km-2
     if ( abs(kord)<9 ) then
       do i=i1,i2
! Left  edges
          pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
          lac_1 = pmp_1 + 1.5*gam(i,k+2)
          a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                         max(a4(1,i,k), pmp_1, lac_1) )
! Right edges
          pmp_2 = a4(1,i,k) + 2.*gam(i,k)
          lac_2 = pmp_2 - 1.5*gam(i,k-1)
          a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                         max(a4(1,i,k), pmp_2, lac_2) )

          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo

     elseif ( abs(kord)==9 ) then
       do i=i1,i2
          if ( extm(i,k) .and. extm(i,k-1) ) then  ! c90_mp122
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else if ( extm(i,k) .and. extm(i,k+1) ) then  ! c90_mp122
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==10 ) then
       do i=i1,i2
          if( ext5(i,k) ) then
              if( ext5(i,k-1) .or. ext5(i,k+1) ) then
                   a4(2,i,k) = a4(1,i,k)
                   a4(3,i,k) = a4(1,i,k)
              elseif ( ext6(i,k-1) .or. ext6(i,k+1) ) then
                   pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                   lac_1 = pmp_1 + 1.5*gam(i,k+2)
                   a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                                  max(a4(1,i,k), pmp_1, lac_1) )
                   pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                   lac_2 = pmp_2 - 1.5*gam(i,k-1)
                   a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                                  max(a4(1,i,k), pmp_2, lac_2) )
              endif
          elseif( ext6(i,k) ) then
              if( ext5(i,k-1) .or. ext5(i,k+1) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
                  a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                                 max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
                  a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                                 max(a4(1,i,k), pmp_2, lac_2) )
              endif
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     elseif ( abs(kord)==12 ) then
       do i=i1,i2
          if( extm(i,k) ) then
! grid-scale 2-delta-z wave detected
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
          else        ! not a local extremum
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==13 ) then
       do i=i1,i2
          if( ext6(i,k) ) then
             if ( ext6(i,k-1) .and. ext6(i,k+1) ) then
! grid-scale 2-delta-z wave detected
                 a4(2,i,k) = a4(1,i,k)
                 a4(3,i,k) = a4(1,i,k)
             endif
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     elseif ( abs(kord)==14 ) then

       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo

     elseif ( abs(kord)==15 ) then   ! revised kord=9 scehem
       do i=i1,i2
          if ( ext5(i,k) ) then  ! c90_mp122
             if ( ext5(i,k-1) .or. ext5(i,k+1) ) then  ! c90_mp122
! grid-scale 2-delta-z wave detected
                  a4(2,i,k) = a4(1,i,k)
                  a4(3,i,k) = a4(1,i,k)
             endif
          elseif( ext6(i,k) ) then
! Check within the smooth region if subgrid profile is non-monotonic
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     elseif ( abs(kord)==16 ) then
       do i=i1,i2
          if( ext5(i,k) ) then
             if ( ext5(i,k-1) .or. ext5(i,k+1) ) then
                 a4(2,i,k) = a4(1,i,k)
                 a4(3,i,k) = a4(1,i,k)
             elseif ( ext6(i,k-1) .or. ext6(i,k+1) ) then
                 ! Left  edges
                 pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                 lac_1 = pmp_1 + 1.5*gam(i,k+2)
                 a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                     max(a4(1,i,k), pmp_1, lac_1) )
                 ! Right edges
                 pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                 lac_2 = pmp_2 - 1.5*gam(i,k-1)
                 a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                     max(a4(1,i,k), pmp_2, lac_2) )
             endif
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     else      ! kord = 11
       do i=i1,i2
         if ( ext5(i,k) .and. (ext5(i,k-1) .or. ext5(i,k+1)) ) then
! Noisy region:
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
         else
              a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         endif
       enddo
     endif

! Additional constraint to ensure positivity
     if ( iv==0 ) call cs_limiters_real8(im, extm(i1,k), a4(1,i1,k), 0)

  enddo      ! k-loop

!----------------------------------
! Bottom layer subgrid constraints:
!----------------------------------
  if ( iv==0 ) then
     do i=i1,i2
        a4(3,i,km) = max(0., a4(3,i,km))
     enddo
  elseif ( iv .eq. -1 ) then
      do i=i1,i2
         if ( a4(3,i,km)*a4(1,i,km) <= 0. )  a4(3,i,km) = 0.
      enddo
  endif

  do k=km-1,km
     do i=i1,i2
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
     if(k==(km-1)) call cs_limiters_real8(im, extm(i1,k), a4(1,i1,k), 2)
     if(k== km   ) call cs_limiters_real8(im, extm(i1,k), a4(1,i1,k), 1)
  enddo

 end subroutine cs_profile_real8

 subroutine cs_limiters_real8(im, extm, a4, iv)
 integer, intent(in) :: im
 integer, intent(in) :: iv
 logical, intent(in) :: extm(im)
 real(kind=8) , intent(inout) :: a4(4,im)   ! PPM array
! !LOCAL VARIABLES:
 real(kind=8)  da1, da2, a6da
 integer i

 if ( iv==0 ) then
! Positive definite constraint
    do i=1,im
    if( a4(1,i)<=0.) then
        a4(2,i) = a4(1,i)
        a4(3,i) = a4(1,i)
        a4(4,i) = 0.
    else
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
         if( (a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12_real8) < 0. ) then
! local minimum is negative
             if( a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i) ) then
                 a4(3,i) = a4(1,i)
                 a4(2,i) = a4(1,i)
                 a4(4,i) = 0.
             elseif( a4(3,i) > a4(2,i) ) then
                 a4(4,i) = 3.*(a4(2,i)-a4(1,i))
                 a4(3,i) = a4(2,i) - a4(4,i)
             else
                 a4(4,i) = 3.*(a4(3,i)-a4(1,i))
                 a4(2,i) = a4(3,i) - a4(4,i)
             endif
         endif
      endif
    endif
    enddo
 elseif ( iv==1 ) then
    do i=1,im
      if( (a4(1,i)-a4(2,i))*(a4(1,i)-a4(3,i))>=0. ) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
    enddo
 else
! Standard PPM constraint
    do i=1,im
      if( extm(i) ) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
    enddo
 endif
 end subroutine cs_limiters_real8

 subroutine ppm_profile_real4(a4, delp, km, i1, i2, iv, kord)

! !INPUT PARAMETERS:
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
                               ! iv = 2: temp (if remap_t) and w (iv=-2)
 integer, intent(in):: i1      ! Starting longitude
 integer, intent(in):: i2      ! Finishing longitude
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: kord    ! Order (or more accurately method no.):
                               !
 real(kind=4) , intent(in):: delp(i1:i2,km)     ! layer pressure thickness

! !INPUT/OUTPUT PARAMETERS:
 real(kind=4) , intent(inout):: a4(4,i1:i2,km)  ! Interpolated values

! DESCRIPTION:
!
!   Perform the piecewise parabolic reconstruction
!
! !REVISION HISTORY:
! S.-J. Lin   revised at GFDL 2007
!-----------------------------------------------------------------------
! local arrays:
      real(kind=4)    dc(i1:i2,km)
      real(kind=4)    h2(i1:i2,km)
      real(kind=4)  delq(i1:i2,km)
      real(kind=4)   df2(i1:i2,km)
      real(kind=4)    d4(i1:i2,km)

! local scalars:
      integer i, k, km1, lmt, it
      real(kind=4)  fac
      real(kind=4)  a1, a2, c1, c2, c3, d1, d2
      real(kind=4)  qm, dq, lac, qmp, pmp

      km1 = km - 1
       it = i2 - i1 + 1

      do k=2,km
         do i=i1,i2
            delq(i,k-1) =   a4(1,i,k) - a4(1,i,k-1)
              d4(i,k  ) = delp(i,k-1) + delp(i,k)
         enddo
      enddo

      do k=2,km1
         do i=i1,i2
                 c1  = (delp(i,k-1)+0.5*delp(i,k))/d4(i,k+1)
                 c2  = (delp(i,k+1)+0.5*delp(i,k))/d4(i,k)
            df2(i,k) = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /      &
                                    (d4(i,k)+delp(i,k+1))
            dc(i,k) = sign( min(abs(df2(i,k)),              &
                            max(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))-a4(1,i,k),  &
                  a4(1,i,k)-min(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))), df2(i,k) )
         enddo
      enddo

!-----------------------------------------------------------
! 4th order interpolation of the provisional cell edge value
!-----------------------------------------------------------

      do k=3,km1
         do i=i1,i2
            c1 = delq(i,k-1)*delp(i,k-1) / d4(i,k)
            a1 = d4(i,k-1) / (d4(i,k) + delp(i,k-1))
            a2 = d4(i,k+1) / (d4(i,k) + delp(i,k))
            a4(2,i,k) = a4(1,i,k-1) + c1 + 2./(d4(i,k-1)+d4(i,k+1)) *    &
                      ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -          &
                        delp(i,k-1)*a1*dc(i,k  ) )
         enddo
      enddo

!     if(km>8 .and. kord>4) call steepz(i1, i2, km, a4, df2, dc, delq, delp, d4)

! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top
      do i=i1,i2
         d1 = delp(i,1)
         d2 = delp(i,2)
         qm = (d2*a4(1,i,1)+d1*a4(1,i,2)) / (d1+d2)
         dq = 2.*(a4(1,i,2)-a4(1,i,1)) / (d1+d2)
         c1 = 4.*(a4(2,i,3)-qm-d2*dq) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) )
         c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1*d1)
         a4(2,i,2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
! Top edge:
!-------------------------------------------------------
         a4(2,i,1) = d1*(2.*c1*d1**2-c3) + a4(2,i,2)
!-------------------------------------------------------
!        a4(2,i,1) = (12./7.)*a4(1,i,1)-(13./14.)*a4(1,i,2)+(3./14.)*a4(1,i,3)
!-------------------------------------------------------
! No over- and undershoot condition
         a4(2,i,2) = max( a4(2,i,2), min(a4(1,i,1), a4(1,i,2)) )
         a4(2,i,2) = min( a4(2,i,2), max(a4(1,i,1), a4(1,i,2)) )
         dc(i,1) =  0.5*(a4(2,i,2) - a4(1,i,1))
      enddo

! Enforce monotonicity  within the top layer

      if( iv==0 ) then
         do i=i1,i2
            a4(2,i,1) = max(0., a4(2,i,1))
            a4(2,i,2) = max(0., a4(2,i,2))
         enddo
      elseif( iv==-1 ) then
         do i=i1,i2
            if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
         enddo
      elseif( abs(iv)==2 ) then
         do i=i1,i2
            a4(2,i,1) = a4(1,i,1)
            a4(3,i,1) = a4(1,i,1)
         enddo
      endif

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
      do i=i1,i2
         d1 = delp(i,km)
         d2 = delp(i,km1)
         qm = (d2*a4(1,i,km)+d1*a4(1,i,km1)) / (d1+d2)
         dq = 2.*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2)
         c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
         c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1*d1)
         a4(2,i,km) = qm - c1*d1*d2*(d2+3.*d1)
! Bottom edge:
!-----------------------------------------------------
         a4(3,i,km) = d1*(8.*c1*d1**2-c3) + a4(2,i,km)
!        dc(i,km) = 0.5*(a4(3,i,km) - a4(1,i,km))
!-----------------------------------------------------
!        a4(3,i,km) = (12./7.)*a4(1,i,km)-(13./14.)*a4(1,i,km-1)+(3./14.)*a4(1,i,km-2)
! No over- and under-shoot condition
         a4(2,i,km) = max( a4(2,i,km), min(a4(1,i,km), a4(1,i,km1)) )
         a4(2,i,km) = min( a4(2,i,km), max(a4(1,i,km), a4(1,i,km1)) )
         dc(i,km) = 0.5*(a4(1,i,km) - a4(2,i,km))
      enddo


! Enforce constraint on the "slope" at the surface

#ifdef BOT_MONO
      do i=i1,i2
         a4(4,i,km) = 0
         if( a4(3,i,km) * a4(1,i,km) <= 0. ) a4(3,i,km) = 0.
         d1 = a4(1,i,km) - a4(2,i,km)
         d2 = a4(3,i,km) - a4(1,i,km)
         if ( d1*d2 < 0. ) then
              a4(2,i,km) = a4(1,i,km)
              a4(3,i,km) = a4(1,i,km)
         else
              dq = sign(min(abs(d1),abs(d2),0.5*abs(delq(i,km-1))), d1)
              a4(2,i,km) = a4(1,i,km) - dq
              a4(3,i,km) = a4(1,i,km) + dq
         endif
      enddo
#else
      if( iv==0 ) then
          do i=i1,i2
             a4(2,i,km) = max(0.,a4(2,i,km))
             a4(3,i,km) = max(0.,a4(3,i,km))
          enddo
      elseif( iv<0 ) then
          do i=i1,i2
             if( a4(1,i,km)*a4(3,i,km) <= 0. )  a4(3,i,km) = 0.
          enddo
      endif
#endif

   do k=1,km1
      do i=i1,i2
         a4(3,i,k) = a4(2,i,k+1)
      enddo
   enddo

!-----------------------------------------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
!-----------------------------------------------------------
! Top 2 and bottom 2 layers always use monotonic mapping
      do k=1,2
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call ppm_limiters_real4(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

      if(kord >= 7) then
!-----------------------
! Huynh's 2nd constraint
!-----------------------
      do k=2,km1
         do i=i1,i2
! Method#1
!           h2(i,k) = delq(i,k) - delq(i,k-1)
! Method#2 - better
            h2(i,k) = 2.*(dc(i,k+1)/delp(i,k+1) - dc(i,k-1)/delp(i,k-1))  &
                     / ( delp(i,k)+0.5*(delp(i,k-1)+delp(i,k+1)) )        &
                     * delp(i,k)**2
! Method#3
!!!            h2(i,k) = dc(i,k+1) - dc(i,k-1)
         enddo
      enddo

      fac = 1.5           ! original quasi-monotone

      do k=3,km-2
        do i=i1,i2
! Right edges
!        qmp   = a4(1,i,k) + 2.0*delq(i,k-1)
!        lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
!
         pmp   = 2.*dc(i,k)
         qmp   = a4(1,i,k) + pmp
         lac   = a4(1,i,k) + fac*h2(i,k-1) + dc(i,k)
         a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), qmp, lac)),    &
                                        max(a4(1,i,k), qmp, lac) )
! Left  edges
!        qmp   = a4(1,i,k) - 2.0*delq(i,k)
!        lac   = a4(1,i,k) + fac*h2(i,k+1) - 0.5*delq(i,k)
!
         qmp   = a4(1,i,k) - pmp
         lac   = a4(1,i,k) + fac*h2(i,k+1) - dc(i,k)
         a4(2,i,k) = min(max(a4(2,i,k),  min(a4(1,i,k), qmp, lac)),   &
                     max(a4(1,i,k), qmp, lac))
!-------------
! Recompute A6
!-------------
         a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
! Additional constraint to ensure positivity when kord=7
         if (iv == 0 .and. kord >= 6 )                      &
             call ppm_limiters_real4(dc(i1,k), a4(1,i1,k), it, 2)
      enddo

      else

         lmt = kord - 3
         lmt = max(0, lmt)
         if (iv == 0) lmt = min(2, lmt)

         do k=3,km-2
            if( kord /= 4) then
              do i=i1,i2
                 a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
              enddo
            endif
            if(kord/=6) call ppm_limiters_real4(dc(i1,k), a4(1,i1,k), it, lmt)
         enddo
      endif

      do k=km1,km
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call ppm_limiters_real4(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

 end subroutine ppm_profile_real4

 subroutine ppm_limiters_real4(dm, a4, itot, lmt)

! !INPUT PARAMETERS:
      real(kind=4) , intent(in):: dm(*)     ! the linear slope
      integer, intent(in) :: itot      ! Total Longitudes
      integer, intent(in) :: lmt       ! 0: Standard PPM constraint
                                       ! 1: Improved full monotonicity constraint (Lin)
                                       ! 2: Positive definite constraint
                                       ! 3: do nothing (return immediately)
! !INPUT/OUTPUT PARAMETERS:
      real(kind=4) , intent(inout) :: a4(4,*)   ! PPM array
                                           ! AA <-- a4(1,i)
                                           ! AL <-- a4(2,i)
                                           ! AR <-- a4(3,i)
                                           ! A6 <-- a4(4,i)
! !LOCAL VARIABLES:
      real(kind=4)  qmp
      real(kind=4)  da1, da2, a6da
      real(kind=4)  fmin
      integer i

! Developer: S.-J. Lin

      if ( lmt == 3 ) return

      if(lmt == 0) then
! Standard PPM constraint
      do i=1,itot
      if(dm(i) == 0.) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
      enddo

      elseif (lmt == 1) then

! Improved full monotonicity constraint (Lin 2004)
! Note: no need to provide first guess of A6 <-- a4(4,i)
      do i=1, itot
           qmp = 2.*dm(i)
         a4(2,i) = a4(1,i)-sign(min(abs(qmp),abs(a4(2,i)-a4(1,i))), qmp)
         a4(3,i) = a4(1,i)+sign(min(abs(qmp),abs(a4(3,i)-a4(1,i))), qmp)
         a4(4,i) = 3.*( 2.*a4(1,i) - (a4(2,i)+a4(3,i)) )
      enddo

      elseif (lmt == 2) then

! Positive definite constraint
      do i=1,itot
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
      fmin = a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12_real4
         if( fmin < 0. ) then
         if(a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i)) then
            a4(3,i) = a4(1,i)
            a4(2,i) = a4(1,i)
            a4(4,i) = 0.
         elseif(a4(3,i) > a4(2,i)) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         else
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
         endif
      endif
      enddo

      endif

 end subroutine ppm_limiters_real4

 subroutine ppm_profile_real8(a4, delp, km, i1, i2, iv, kord)

! !INPUT PARAMETERS:
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
                               ! iv = 2: temp (if remap_t) and w (iv=-2)
 integer, intent(in):: i1      ! Starting longitude
 integer, intent(in):: i2      ! Finishing longitude
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: kord    ! Order (or more accurately method no.):
                               !
 real(kind=8) , intent(in):: delp(i1:i2,km)     ! layer pressure thickness

! !INPUT/OUTPUT PARAMETERS:
 real(kind=8) , intent(inout):: a4(4,i1:i2,km)  ! Interpolated values

! DESCRIPTION:
!
!   Perform the piecewise parabolic reconstruction
!
! !REVISION HISTORY:
! S.-J. Lin   revised at GFDL 2007
!-----------------------------------------------------------------------
! local arrays:
      real(kind=8)    dc(i1:i2,km)
      real(kind=8)    h2(i1:i2,km)
      real(kind=8)  delq(i1:i2,km)
      real(kind=8)   df2(i1:i2,km)
      real(kind=8)    d4(i1:i2,km)

! local scalars:
      integer i, k, km1, lmt, it
      real(kind=8)  fac
      real(kind=8)  a1, a2, c1, c2, c3, d1, d2
      real(kind=8)  qm, dq, lac, qmp, pmp

      km1 = km - 1
       it = i2 - i1 + 1

      do k=2,km
         do i=i1,i2
            delq(i,k-1) =   a4(1,i,k) - a4(1,i,k-1)
              d4(i,k  ) = delp(i,k-1) + delp(i,k)
         enddo
      enddo

      do k=2,km1
         do i=i1,i2
                 c1  = (delp(i,k-1)+0.5*delp(i,k))/d4(i,k+1)
                 c2  = (delp(i,k+1)+0.5*delp(i,k))/d4(i,k)
            df2(i,k) = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /      &
                                    (d4(i,k)+delp(i,k+1))
            dc(i,k) = sign( min(abs(df2(i,k)),              &
                            max(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))-a4(1,i,k),  &
                  a4(1,i,k)-min(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))), df2(i,k) )
         enddo
      enddo

!-----------------------------------------------------------
! 4th order interpolation of the provisional cell edge value
!-----------------------------------------------------------

      do k=3,km1
         do i=i1,i2
            c1 = delq(i,k-1)*delp(i,k-1) / d4(i,k)
            a1 = d4(i,k-1) / (d4(i,k) + delp(i,k-1))
            a2 = d4(i,k+1) / (d4(i,k) + delp(i,k))
            a4(2,i,k) = a4(1,i,k-1) + c1 + 2./(d4(i,k-1)+d4(i,k+1)) *    &
                      ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -          &
                        delp(i,k-1)*a1*dc(i,k  ) )
         enddo
      enddo

!     if(km>8 .and. kord>4) call steepz(i1, i2, km, a4, df2, dc, delq, delp, d4)

! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top
      do i=i1,i2
         d1 = delp(i,1)
         d2 = delp(i,2)
         qm = (d2*a4(1,i,1)+d1*a4(1,i,2)) / (d1+d2)
         dq = 2.*(a4(1,i,2)-a4(1,i,1)) / (d1+d2)
         c1 = 4.*(a4(2,i,3)-qm-d2*dq) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) )
         c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1*d1)
         a4(2,i,2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
! Top edge:
!-------------------------------------------------------
         a4(2,i,1) = d1*(2.*c1*d1**2-c3) + a4(2,i,2)
!-------------------------------------------------------
!        a4(2,i,1) = (12./7.)*a4(1,i,1)-(13./14.)*a4(1,i,2)+(3./14.)*a4(1,i,3)
!-------------------------------------------------------
! No over- and undershoot condition
         a4(2,i,2) = max( a4(2,i,2), min(a4(1,i,1), a4(1,i,2)) )
         a4(2,i,2) = min( a4(2,i,2), max(a4(1,i,1), a4(1,i,2)) )
         dc(i,1) =  0.5*(a4(2,i,2) - a4(1,i,1))
      enddo

! Enforce monotonicity  within the top layer

      if( iv==0 ) then
         do i=i1,i2
            a4(2,i,1) = max(0., a4(2,i,1))
            a4(2,i,2) = max(0., a4(2,i,2))
         enddo
      elseif( iv==-1 ) then
         do i=i1,i2
            if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
         enddo
      elseif( abs(iv)==2 ) then
         do i=i1,i2
            a4(2,i,1) = a4(1,i,1)
            a4(3,i,1) = a4(1,i,1)
         enddo
      endif

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
      do i=i1,i2
         d1 = delp(i,km)
         d2 = delp(i,km1)
         qm = (d2*a4(1,i,km)+d1*a4(1,i,km1)) / (d1+d2)
         dq = 2.*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2)
         c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
         c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1*d1)
         a4(2,i,km) = qm - c1*d1*d2*(d2+3.*d1)
! Bottom edge:
!-----------------------------------------------------
         a4(3,i,km) = d1*(8.*c1*d1**2-c3) + a4(2,i,km)
!        dc(i,km) = 0.5*(a4(3,i,km) - a4(1,i,km))
!-----------------------------------------------------
!        a4(3,i,km) = (12./7.)*a4(1,i,km)-(13./14.)*a4(1,i,km-1)+(3./14.)*a4(1,i,km-2)
! No over- and under-shoot condition
         a4(2,i,km) = max( a4(2,i,km), min(a4(1,i,km), a4(1,i,km1)) )
         a4(2,i,km) = min( a4(2,i,km), max(a4(1,i,km), a4(1,i,km1)) )
         dc(i,km) = 0.5*(a4(1,i,km) - a4(2,i,km))
      enddo


! Enforce constraint on the "slope" at the surface

#ifdef BOT_MONO
      do i=i1,i2
         a4(4,i,km) = 0
         if( a4(3,i,km) * a4(1,i,km) <= 0. ) a4(3,i,km) = 0.
         d1 = a4(1,i,km) - a4(2,i,km)
         d2 = a4(3,i,km) - a4(1,i,km)
         if ( d1*d2 < 0. ) then
              a4(2,i,km) = a4(1,i,km)
              a4(3,i,km) = a4(1,i,km)
         else
              dq = sign(min(abs(d1),abs(d2),0.5*abs(delq(i,km-1))), d1)
              a4(2,i,km) = a4(1,i,km) - dq
              a4(3,i,km) = a4(1,i,km) + dq
         endif
      enddo
#else
      if( iv==0 ) then
          do i=i1,i2
             a4(2,i,km) = max(0.,a4(2,i,km))
             a4(3,i,km) = max(0.,a4(3,i,km))
          enddo
      elseif( iv<0 ) then
          do i=i1,i2
             if( a4(1,i,km)*a4(3,i,km) <= 0. )  a4(3,i,km) = 0.
          enddo
      endif
#endif

   do k=1,km1
      do i=i1,i2
         a4(3,i,k) = a4(2,i,k+1)
      enddo
   enddo

!-----------------------------------------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
!-----------------------------------------------------------
! Top 2 and bottom 2 layers always use monotonic mapping
      do k=1,2
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call ppm_limiters_real8(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

      if(kord >= 7) then
!-----------------------
! Huynh's 2nd constraint
!-----------------------
      do k=2,km1
         do i=i1,i2
! Method#1
!           h2(i,k) = delq(i,k) - delq(i,k-1)
! Method#2 - better
            h2(i,k) = 2.*(dc(i,k+1)/delp(i,k+1) - dc(i,k-1)/delp(i,k-1))  &
                     / ( delp(i,k)+0.5*(delp(i,k-1)+delp(i,k+1)) )        &
                     * delp(i,k)**2
! Method#3
!!!            h2(i,k) = dc(i,k+1) - dc(i,k-1)
         enddo
      enddo

      fac = 1.5           ! original quasi-monotone

      do k=3,km-2
        do i=i1,i2
! Right edges
!        qmp   = a4(1,i,k) + 2.0*delq(i,k-1)
!        lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
!
         pmp   = 2.*dc(i,k)
         qmp   = a4(1,i,k) + pmp
         lac   = a4(1,i,k) + fac*h2(i,k-1) + dc(i,k)
         a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), qmp, lac)),    &
                                        max(a4(1,i,k), qmp, lac) )
! Left  edges
!        qmp   = a4(1,i,k) - 2.0*delq(i,k)
!        lac   = a4(1,i,k) + fac*h2(i,k+1) - 0.5*delq(i,k)
!
         qmp   = a4(1,i,k) - pmp
         lac   = a4(1,i,k) + fac*h2(i,k+1) - dc(i,k)
         a4(2,i,k) = min(max(a4(2,i,k),  min(a4(1,i,k), qmp, lac)),   &
                     max(a4(1,i,k), qmp, lac))
!-------------
! Recompute A6
!-------------
         a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
! Additional constraint to ensure positivity when kord=7
         if (iv == 0 .and. kord >= 6 )                      &
             call ppm_limiters_real8(dc(i1,k), a4(1,i1,k), it, 2)
      enddo

      else

         lmt = kord - 3
         lmt = max(0, lmt)
         if (iv == 0) lmt = min(2, lmt)

         do k=3,km-2
            if( kord /= 4) then
              do i=i1,i2
                 a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
              enddo
            endif
            if(kord/=6) call ppm_limiters_real8(dc(i1,k), a4(1,i1,k), it, lmt)
         enddo
      endif

      do k=km1,km
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call ppm_limiters_real8(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

 end subroutine ppm_profile_real8


 subroutine ppm_limiters_real8(dm, a4, itot, lmt)

! !INPUT PARAMETERS:
      real(kind=8) , intent(in):: dm(*)     ! the linear slope
      integer, intent(in) :: itot      ! Total Longitudes
      integer, intent(in) :: lmt       ! 0: Standard PPM constraint
                                       ! 1: Improved full monotonicity constraint (Lin)
                                       ! 2: Positive definite constraint
                                       ! 3: do nothing (return immediately)
! !INPUT/OUTPUT PARAMETERS:
      real(kind=8) , intent(inout) :: a4(4,*)   ! PPM array
                                           ! AA <-- a4(1,i)
                                           ! AL <-- a4(2,i)
                                           ! AR <-- a4(3,i)
                                           ! A6 <-- a4(4,i)
! !LOCAL VARIABLES:
      real(kind=8)  qmp
      real(kind=8)  da1, da2, a6da
      real(kind=8)  fmin
      integer i

! Developer: S.-J. Lin

      if ( lmt == 3 ) return

      if(lmt == 0) then
! Standard PPM constraint
      do i=1,itot
      if(dm(i) == 0.) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
      enddo

      elseif (lmt == 1) then

! Improved full monotonicity constraint (Lin 2004)
! Note: no need to provide first guess of A6 <-- a4(4,i)
      do i=1, itot
           qmp = 2.*dm(i)
         a4(2,i) = a4(1,i)-sign(min(abs(qmp),abs(a4(2,i)-a4(1,i))), qmp)
         a4(3,i) = a4(1,i)+sign(min(abs(qmp),abs(a4(3,i)-a4(1,i))), qmp)
         a4(4,i) = 3.*( 2.*a4(1,i) - (a4(2,i)+a4(3,i)) )
      enddo

      elseif (lmt == 2) then

! Positive definite constraint
      do i=1,itot
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
      fmin = a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12_real8
         if( fmin < 0. ) then
         if(a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i)) then
            a4(3,i) = a4(1,i)
            a4(2,i) = a4(1,i)
            a4(4,i) = 0.
         elseif(a4(3,i) > a4(2,i)) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         else
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
         endif
      endif
      enddo

      endif

 end subroutine ppm_limiters_real8
end module coarse_graining_mod
