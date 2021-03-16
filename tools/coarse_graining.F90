module coarse_graining_mod

  use fms_mod, only: check_nml_error, close_file, open_namelist_file
  use mpp_domains_mod, only: domain2d, mpp_define_io_domain, mpp_define_mosaic, mpp_get_compute_domain
  use fv_mapz_mod, only: mappm
  use mpp_mod, only: FATAL, input_nml_file, mpp_error, mpp_npes

  implicit none
  private

  public :: block_sum, compute_mass_weights, get_fine_array_bounds, &
       get_coarse_array_bounds, coarse_graining_init, weighted_block_average, &
       weighted_block_edge_average_x, weighted_block_edge_average_y, MODEL_LEVEL, &
       block_upsample, mask_area_weights, PRESSURE_LEVEL, vertical_remapping_requirements, &
       vertically_remap_field, mask_mass_weights, remap_edges_along_x, remap_edges_along_y, &
       block_edge_sum_x, block_edge_sum_y

  interface block_sum
     module procedure block_sum_2d
  end interface block_sum

  interface block_edge_sum_x
     module procedure block_edge_sum_x_2d_full_input
  end interface block_edge_sum_x

  interface block_edge_sum_y
     module procedure block_edge_sum_y_2d_full_input
  end interface block_edge_sum_y

  interface weighted_block_average
     module procedure weighted_block_average_2d
     module procedure weighted_block_average_3d_field_2d_weights
     module procedure weighted_block_average_3d_field_3d_weights
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
     module procedure block_upsample_2d
     module procedure block_upsample_3d
  end interface block_upsample

  interface weighted_block_edge_average_x_pre_downsampled
     module procedure weighted_block_edge_average_x_pre_downsampled_unmasked
     module procedure weighted_block_edge_average_x_pre_downsampled_masked
  end interface weighted_block_edge_average_x_pre_downsampled

  interface weighted_block_edge_average_y_pre_downsampled
     module procedure weighted_block_edge_average_y_pre_downsampled_unmasked
     module procedure weighted_block_edge_average_y_pre_downsampled_masked
  end interface weighted_block_edge_average_y_pre_downsampled

  ! Global variables for the module, initialized in coarse_graining_init
  integer :: is, ie, js, je, npz
  integer :: is_coarse, ie_coarse, js_coarse, je_coarse
  character(len=11) :: MODEL_LEVEL = 'model_level'
  character(len=14) :: PRESSURE_LEVEL = 'pressure_level'

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

  subroutine block_sum_2d(fine, coarse)
    real, intent(in) :: fine(is:ie,js:je)
    real, intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse)

    integer :: i, j, i_coarse, j_coarse, offset

    offset = coarsening_factor - 1
    do i = is, ie, coarsening_factor
       i_coarse = (i - 1) / coarsening_factor + 1
       do j = js, je, coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          coarse(i_coarse,j_coarse) = sum(fine(i:i+offset,j:j+offset))
       enddo
    enddo
  end subroutine

  subroutine weighted_block_average_2d(weights, fine, coarse)
    real, intent(in) :: weights(is:ie,js:je), fine(is:ie,js:je)
    real, intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse)

    real, allocatable :: weighted_fine(:,:), weighted_block_sum(:,:), block_sum_weights(:,:)

    allocate(weighted_fine(is:ie,js:je))
    allocate(weighted_block_sum(is_coarse:ie_coarse,js_coarse:je_coarse))
    allocate(block_sum_weights(is_coarse:ie_coarse,js_coarse:je_coarse))

    weighted_fine = weights * fine
    call block_sum_2d(weighted_fine, weighted_block_sum)
    call block_sum_2d(weights, block_sum_weights)
    coarse = weighted_block_sum / block_sum_weights
  end subroutine weighted_block_average_2d

  subroutine weighted_block_average_3d_field_2d_weights(weights, fine, coarse)
    real, intent(in) :: weights(is:ie,js:je), fine(is:ie,js:je,1:npz)
    real, intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz)

    integer :: k

    do k = 1, npz
       call weighted_block_average_2d(weights, fine(is:ie,js:je,k), coarse(is_coarse:ie_coarse,js_coarse:je_coarse,k))
    enddo
  end subroutine weighted_block_average_3d_field_2d_weights

  subroutine weighted_block_average_3d_field_3d_weights(weights, fine, coarse)
    real, intent(in) :: weights(is:ie,js:je,1:npz), fine(is:ie,js:je,1:npz)
    real, intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz)

    integer :: k

    do k = 1, npz
       call weighted_block_average_2d(weights(is:ie,js:je,k), fine(is:ie,js:je,k), coarse(is_coarse:ie_coarse,js_coarse:je_coarse,k))
    enddo
  end subroutine weighted_block_average_3d_field_3d_weights

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

  subroutine vertically_remap_field(phalf_in, field, phalf_out, ptop, field_out)
    real, intent(in) :: phalf_in(is:ie,js:je,1:npz+1), phalf_out(is:ie,js:je,1:npz+1)
    real, intent(in) :: field(is:ie,js:je,1:npz)
    real, intent(in) :: ptop
    real, intent(out) :: field_out(is:ie,js:je,1:npz)

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
  end subroutine vertically_remap_field

  subroutine block_upsample_2d(coarse, fine)
    real, intent(in) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse)
    real, intent(out) :: fine(is:ie,js:je)

    integer :: i, j, i_coarse, j_coarse, offset

    offset = coarsening_factor - 1
    do i = is, ie, coarsening_factor
      i_coarse = (i - 1) / coarsening_factor + 1
      do j = js, je, coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          fine(i:i+offset,j:j+offset) = coarse(i_coarse, j_coarse)
      enddo
    enddo
  end subroutine block_upsample_2d

  subroutine block_upsample_3d(coarse, fine, nz)
    integer, intent(in) :: nz
    real, intent(in) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:nz)
    real, intent(out) :: fine(is:ie,js:je,1:nz)

    integer :: k

    do k = 1, nz
      call block_upsample_2d(coarse(is_coarse:ie_coarse,js_coarse:je_coarse,k), fine(is:ie,js:je,k))
    enddo
  end subroutine block_upsample_3d

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

  subroutine compute_phalf_from_delp(delp, ptop, i_start, i_end, j_start, j_end, phalf)
    integer, intent(in) :: i_start, i_end, j_start, j_end
    real, intent(in) :: delp(i_start:i_end,j_start:j_end,1:npz)
    real, intent(in) :: ptop
    real, intent(out) :: phalf(i_start:i_end,j_start:j_end,1:npz+1)

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
  end subroutine compute_phalf_from_delp

 ! Routine for computing the common requirements for pressure-level coarse-graining.
  subroutine vertical_remapping_requirements(delp, area, ptop, phalf, upsampled_coarse_phalf)
    real, intent(in) :: delp(is:ie,js:je,1:npz)
    real, intent(in) :: area(is:ie,js:je)
    real, intent(in) :: ptop
    real, intent(out) :: phalf(is:ie,js:je,1:npz+1)
    real, intent(out) :: upsampled_coarse_phalf(is:ie,js:je,1:npz+1)

    real, allocatable :: coarse_delp(:,:,:), coarse_phalf(:,:,:)

    allocate(coarse_delp(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz))
    allocate(coarse_phalf(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz+1))

    call compute_phalf_from_delp(delp(is:ie,js:je,1:npz), ptop, is, ie, js, je, phalf)
    call weighted_block_average(area(is:ie,js:je), delp(is:ie,js:je,1:npz), coarse_delp)
    call compute_phalf_from_delp(coarse_delp, ptop, is_coarse, ie_coarse, js_coarse, je_coarse, coarse_phalf)
    call block_upsample(coarse_phalf, upsampled_coarse_phalf, npz+1)

    deallocate(coarse_delp)
    deallocate(coarse_phalf)
   end subroutine vertical_remapping_requirements

   subroutine mask_area_weights(area, phalf, upsampled_coarse_phalf, masked_area_weights)
    real, intent(in) :: area(is:ie,js:je)
    real, intent(in) :: phalf(is:ie,js:je,1:npz+1)
    real, intent(in) :: upsampled_coarse_phalf(is:ie,js:je,1:npz+1)
    real, intent(out) :: masked_area_weights(is:ie,js:je,1:npz)

    integer :: k

    do k = 1, npz
      where (upsampled_coarse_phalf(is:ie,js:je,k+1) .lt. phalf(is:ie,js:je,npz+1))
        masked_area_weights(is:ie,js:je,k) = area(is:ie,js:je)
      elsewhere
        masked_area_weights(is:ie,js:je,k) = 0.0
      endwhere
    enddo
   end subroutine mask_area_weights

   subroutine mask_mass_weights(area, delp, phalf, upsampled_coarse_phalf, &
    masked_mass_weights)
    real, intent(in) :: area(is:ie,js:je)
    real, intent(in) :: delp(is:ie,js:je,1:npz)
    real, intent(in) :: phalf(is:ie,js:je,1:npz+1)
    real, intent(in) :: upsampled_coarse_phalf(is:ie,js:je,1:npz+1)
    real, intent(out) :: masked_mass_weights(is:ie,js:je,1:npz)

    integer :: k

    do k = 1, npz
      where (upsampled_coarse_phalf(:,:,k+1) .lt. phalf(is:ie,js:je,npz+1))
        masked_mass_weights(:,:,k) = delp(:,:,k) * area(:,:)
      elsewhere
        masked_mass_weights(:,:,k) = 0.0
      endwhere
    enddo
   end subroutine mask_mass_weights

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

end module coarse_graining_mod
