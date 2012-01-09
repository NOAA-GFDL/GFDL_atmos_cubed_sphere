module fv_grid_tools_mod

  use constants_mod, only: radius, pi, omega, grav
  use fv_arrays_mod, only: fv_atmos_type
  use fv_grid_utils_mod, only: gnomonic_grids, great_circle_dist,  &
                           mid_pt_sphere, spherical_angle,     &
                           project_sphere_v,  cell_center2,    &
                           get_area, inner_prod, deglat,       &
                           sw_corner, se_corner, ne_corner, nw_corner, fill_ghost, &
                           Gnomonic_grid, direct_transform
  use fv_timing_mod,  only: timing_on, timing_off
  use fv_mp_mod,      only: gid, masterproc, domain, tile, &
                            is,js,ie,je,isd,jsd,ied,jed, ng, &
                            fill_corners, XDir, YDir, &
                            mp_gather, mp_bcst, mp_reduce_max, mp_stop, &
                            npes_x, npes_y
  use sorted_index_mod,  only: sorted_inta, sorted_intb
  use mpp_mod,           only: mpp_error, FATAL, get_unit, mpp_chksum, mpp_pe, stdout, &
                               mpp_send, mpp_recv, mpp_sync_self, EVENT_RECV, mpp_npes, &
                               mpp_sum, mpp_max, mpp_min, mpp_root_pe, mpp_broadcast
  use mpp_domains_mod,   only: mpp_update_domains, mpp_get_boundary, &
                               mpp_get_ntile_count, mpp_get_pelist, &
                               mpp_get_compute_domains, mpp_global_field
  use mpp_io_mod,        only: mpp_get_att_value     

  use mpp_parameter_mod, only: AGRID_PARAM=>AGRID,       & 
                               CGRID_NE_PARAM=>CGRID_NE, &
                               CGRID_SW_PARAM=>CGRID_SW, &
                               BGRID_NE_PARAM=>BGRID_NE, &
                               BGRID_SW_PARAM=>BGRID_SW, & 
                               SCALAR_PAIR,              &
                               CORNER, CENTER, XUPDATE
  use fms_mod,           only: get_mosaic_tile_grid
  use fms_io_mod,        only: file_exist, field_exist, read_data, &
                               get_global_att_value, get_var_att_value
  use mosaic_mod,       only : get_mosaic_ntiles
  implicit none
  private
#include <netcdf.inc>

  real :: csFac = -999
  real            :: zeta = 1.0                ! non-linear flag 
  real , parameter:: todeg = 180.0/pi          ! convert to degrees
  real , parameter:: torad = pi/180.0          ! convert to radians
  real , parameter:: missing = 1.e25
  real    :: stretch               ! Optional stretching factor for the grid 
  logical :: dxdy_area = .false.   ! define area using dx*dy else spherical excess formula
  logical :: latlon = .false.
  logical :: cubed_sphere = .false.
  logical :: double_periodic = .false.
  logical :: latlon_patch = .false.
  logical :: latlon_strip = .false.
  logical :: channel = .false.
  logical :: have_south_pole = .false.
  logical :: have_north_pole = .false.
  logical :: uniform_ppm = .true.     ! assume uniform grid spacing for PPM calcs, else variable dx,dy
  integer :: interpOrder = 1
  logical :: debug_message_size = .false.
  logical :: write_grid_char_file = .false.
  logical :: stretched_grid = .false.

  ! grid descriptors

  ! Horizontal
  integer :: npx_g, npy_g, npz_g, ntiles_g ! global domain
#ifndef NO_GRID_G
  real, allocatable, target, dimension(:,:,:) :: grid_g
#endif
  real, allocatable, target, dimension(:,:,:) :: grid, agrid
  real, allocatable, dimension(:,:) :: area, area_c
  real, allocatable, dimension(:,:) :: sina, cosa
  real, allocatable, dimension(:,:,:) :: e1,e2
  real, allocatable, dimension(:,:) :: dx, dy
  real, allocatable, dimension(:,:) :: dxc, dyc
  real, allocatable, dimension(:,:) :: dxa, dya
  real, allocatable, dimension(:,:) :: rarea, rarea_c
  real, allocatable, dimension(:,:) :: rdx, rdy
  real, allocatable, dimension(:,:) :: rdxc, rdyc
  real, allocatable, dimension(:,:) :: rdxa, rdya
  real  :: acapN, acapS
  real  :: globalarea  ! total Global Area
  real, allocatable :: cose(:,:)
  real, allocatable :: cosp(:,:)
  real, allocatable :: acosp(:,:)
  
  integer, dimension(:,:,:), allocatable :: iinta, jinta, iintb, jintb
  
  integer :: grid_type = 0    ! -1: read from file; 0: ED Gnomonic
                              !  0: the "true" equal-distance Gnomonic grid
                              !  1: the traditional equal-distance Gnomonic grid
                              !  2: the equal-angular Gnomonic grid
                              !  3: the lat-lon grid -- to be implemented
                              !  4: double periodic boundary condition on Cartesian grid
                              !  5: latlon patch
                              !  6: latlon strip (cyclic in longitude)
                              !  7: channel flow on Cartesian grid

  real :: dx_const = 1000.    ! spatial resolution for double periodic boundary configuration [m]
  real :: dy_const = 1000.
  real :: deglon_start = -30., deglon_stop = 30., &  ! boundaries of latlon patch
          deglat_start = -30., deglat_stop = 30.

#ifndef NO_GRID_G
  public :: grid_g
#endif
  public :: npx_g, npy_g, npz_g, grid, agrid, stretch, todeg, &
            interpOrder, uniform_ppm, zeta, missing, &
            cubed_sphere, latlon, have_south_pole, have_north_pole, &
            double_periodic, channel, &
            dx,dy, dxa,dya, dxc,dyc, rdx,rdy, rdxc,rdyc,  &
            sina, cosa, area, rarea, area_c, rarea_c,  &
            acapN, acapS, cosp, cose, acosp, init_grid, read_grid, &
            rdxa, rdya, d2a2c, ctoa, atod, dtoa, atoc, atob_s,   &
            mp_update_dwinds, rotate_winds, &
            spherical_to_cartesian, globalsum, &
            get_unit_vector, unit_vect2
  public :: grid_type, dx_const, dy_const
  public :: deglon_start, deglon_stop, deglat_start, deglat_stop
  public :: debug_message_size, write_grid_char_file

  INTERFACE get_unit_vector
     MODULE PROCEDURE get_unit_vector_3pts
     MODULE PROCEDURE get_unit_vector_2pts
  END INTERFACE

  INTERFACE mp_update_dwinds
     MODULE PROCEDURE mp_update_dwinds_2d
     MODULE PROCEDURE mp_update_dwinds_3d
  END INTERFACE

!---- version number -----
  character(len=128) :: version = '$Id: fv_grid_tools.F90,v 19.0 2012/01/06 19:59:07 fms Exp $'
  character(len=128) :: tagname = '$Name: siena $'

contains

  subroutine read_grid(Atm, grid_name, grid_file, npx, npy, npz, ndims, nregions, ng)
    !     read_grid :: read grid from mosaic grid file.
    type(fv_atmos_type), intent(inout) :: Atm
    character(len=*),    intent(IN)    :: grid_name
    character(len=*),    intent(IN)    :: grid_file
    integer,             intent(IN)    :: npx, npy, npz
    integer,             intent(IN)    :: ndims
    integer,             intent(IN)    :: nregions
    integer,             intent(IN)    :: ng

    real, allocatable, dimension(:,:)  :: tmpx, tmpy
    real, allocatable, dimension(:)    :: ebuffer, wbuffer, sbuffer, nbuffer
    character(len=128)                 :: units = ""
    character(len=256)                 :: atm_mosaic, atm_hgrid, grid_form
    character(len=1024)                :: attvalue
    integer                            :: ntiles, i, j
    integer                            :: isg, ieg, jsg, jeg
    integer                            :: isc2, iec2, jsc2, jec2
    real                               :: p1(3), p2(3), p3(3), p4(3)
    integer                            :: start(4), nread(4)
    real                               :: angN,angM,angAV,ang
    real                               :: aspN,aspM,aspAV,asp
    real                               ::  dxN, dxM, dxAV
    real                               :: dx_local, dy_local
    real, allocatable, dimension(:,:)  :: tmp, g_tmp, angs, asps, dxs
    character(len=80)                  :: gcharFile
    integer                            :: fileLun, n, stdunit
    real                               :: p_lL(ndims) ! lower Left
    real                               :: p_uL(ndims) ! upper Left
    real                               :: p_lR(ndims) ! lower Right
    real                               :: p_uR(ndims) ! upper Right
    real                               :: d1, d2, mydx, mydy

    Gnomonic_grid = .true.   
    cubed_sphere = .true.
    uniform_ppm = .true.
    npx_g = npx
    npy_g = npy
    npz_g = npz
    ntiles_g = nregions

    if ( Atm%do_schmidt .and. abs(atm%stretch_fac-1.) > 1.E-5 ) stretched_grid = .true.

    if(.not. file_exist(grid_file)) call mpp_error(FATAL, 'fv_grid_tools(read_grid): file '// &
         trim(grid_file)//' does not exist')

    !--- make sure the grid file is mosaic file.
    if( field_exist(grid_file, 'atm_mosaic_file') .OR. field_exist(grid_file, 'gridfiles') ) then
       stdunit = stdout()
       write(stdunit,*) '==>Note from fv_grid_tools_mod(read_grid): read atmosphere grid from mosaic version grid'
    else
       call mpp_error(FATAL, 'fv_grid_tools(read_grid): neither atm_mosaic_file nor gridfiles exists in file ' &
            //trim(grid_file))
    endif

    if(field_exist(grid_file, 'atm_mosaic_file')) then
       call read_data(grid_file, "atm_mosaic_file", atm_mosaic)
       atm_mosaic = "INPUT/"//trim(atm_mosaic)
    else 
       atm_mosaic = trim(grid_file)
    endif

    call get_mosaic_tile_grid(atm_hgrid, atm_mosaic, domain)

    grid_form = "none"    
    if( get_global_att_value(atm_hgrid, "history", attvalue) ) then
       if( index(attvalue, "gnomonic_ed") > 0) grid_form = "gnomonic_ed"
    endif
    if(grid_form .NE. "gnomonic_ed") call mpp_error(FATAL, &
         "fv_grid_tools(read_grid): the grid should be 'gnomonic_ed' when reading from grid file, contact developer")

    ntiles = get_mosaic_ntiles(atm_mosaic)
    if(ntiles .NE. 6) call mpp_error(FATAL, &
       'fv_grid_tools(read_grid): ntiles should be 6 in mosaic file '//trim(atm_mosaic) )
    if(nregions .NE. 6) call mpp_error(FATAL, &
       'fv_grid_tools(read_grid): nregions should be 6 when reading from mosaic file '//trim(grid_file) )

    !-------------------------------------------------------------------
    !   memory allocation for module variable or public variable
    !------------------------------------------------------------------
    allocate (  area(isd:ied  ,jsd:jed  ) )   ! Cell Centered
    allocate ( rarea(isd:ied  ,jsd:jed  ) )   ! Cell Centered
    
    allocate (  area_c(isd:ied+1,jsd:jed+1) )  ! Cell Corners
    allocate ( rarea_c(isd:ied+1,jsd:jed+1) )  ! Cell Corners
    
    allocate (  dx(isd:ied  ,jsd:jed+1) )
    allocate ( rdx(isd:ied  ,jsd:jed+1) )
    allocate (  dy(isd:ied+1,jsd:jed  ) )
    allocate ( rdy(isd:ied+1,jsd:jed  ) )
    
    allocate (  dxc(isd:ied+1,jsd:jed  ) )
    allocate ( rdxc(isd:ied+1,jsd:jed  ) )
    allocate (  dyc(isd:ied  ,jsd:jed+1) )
    allocate ( rdyc(isd:ied  ,jsd:jed+1) )
    
    allocate (  dxa(isd:ied  ,jsd:jed  ) )
    allocate ( rdxa(isd:ied  ,jsd:jed  ) )
    allocate (  dya(isd:ied  ,jsd:jed  ) )
    allocate ( rdya(isd:ied  ,jsd:jed  ) )
    
    allocate ( grid (isd:ied+1,jsd:jed+1,1:ndims) )
    allocate ( agrid(isd:ied  ,jsd:jed  ,1:ndims) )
    
    Atm%grid  =>grid
    Atm%agrid =>agrid
    
    allocate ( sina(isd:ied+1,jsd:jed+1) )   ! SIN(angle of intersection)
    allocate ( cosa(isd:ied+1,jsd:jed+1) )   ! COS(angle of intersection)
    
    allocate (   e1(3,isd:ied+1,jsd:jed+1) )
    allocate (   e2(3,isd:ied+1,jsd:jed+1) )

    call get_var_att_value(atm_hgrid, 'x', 'units', units)

    !--- get the geographical coordinates of super-grid.
    isc2 = 2*is-1; iec2 = 2*ie+1
    jsc2 = 2*js-1; jec2 = 2*je+1  
    allocate(tmpx(isc2:iec2, jsc2:jec2) )
    allocate(tmpy(isc2:iec2, jsc2:jec2) )
    start = 1; nread = 1
    start(1) = isc2; nread(1) = iec2 - isc2 + 1
    start(2) = jsc2; nread(2) = jec2 - jsc2 + 1
    call read_data(atm_hgrid, 'x', tmpx, start, nread, no_domain=.TRUE.)
    call read_data(atm_hgrid, 'y', tmpy, start, nread, no_domain=.TRUE.)

    !--- geographic grid at cell corner
    grid(isd: is-1, jsd:js-1,1:ndims)=0.
    grid(isd: is-1, je+2:jed+1,1:ndims)=0.
    grid(ie+2:ied+1,jsd:js-1,1:ndims)=0.
    grid(ie+2:ied+1,je+2:jed+1,1:ndims)=0.
    if(len_trim(units) < 6) call mpp_error(FATAL, &
          "fv_grid_tools_mod(read_grid): the length of units must be no less than 6")
    if(units(1:6) == 'degree') then
       do j = js, je+1
          do i = is, ie+1
             grid(i,j,1) = tmpx(2*i-1,2*j-1)*pi/180.
             grid(i,j,2) = tmpy(2*i-1,2*j-1)*pi/180.
          enddo
       enddo
    else if(units(1:6) == 'radian') then
       do j = js, je+1
          do i = is, ie+1
             grid(i,j,1) = tmpx(2*i-1,2*j-1)
             grid(i,j,2) = tmpy(2*i-1,2*j-1)
          enddo
       enddo
    else
       print*, 'units is ' , trim(units), len_trim(units), mpp_pe()
       call mpp_error(FATAL, 'fv_grid_tools_mod(read_grid): units must start with degree or radian')
    endif

    call mpp_update_domains( grid, domain, position=CORNER)    

    !--- geographic grid at cell center
    agrid(:,:,:) = -1.e25
    if(units(1:6) == 'degree') then
       do j = js, je
          do i = is, ie
             agrid(i,j,1) = tmpx(2*i,2*j)*pi/180.
             agrid(i,j,2) = tmpy(2*i,2*j)*pi/180.
          enddo
       enddo
    else if(units(1:6) == 'radian') then
       do j = js, je
          do i = is, ie
             agrid(i,j,1) = tmpx(2*i,2*j)
             agrid(i,j,2) = tmpy(2*i,2*j)
          enddo
       enddo
    endif
    
    call mpp_update_domains( agrid, domain)       
    call fill_corners(agrid(:,:,1), npx, npy, XDir, AGRID=.true.)
    call fill_corners(agrid(:,:,2), npx, npy, YDir, AGRID=.true.)
    deallocate(tmpx, tmpy)

    !--- dx and dy         
    do j = js, je+1
       do i = is, ie
          p1(1) = grid(i  ,j,1)
          p1(2) = grid(i  ,j,2)
          p2(1) = grid(i+1,j,1)
          p2(2) = grid(i+1,j,2)
          dx(i,j) = great_circle_dist( p2, p1, radius )
       enddo
    enddo
    call get_symmetry(dx(is:ie,js:je+1), dy(is:ie+1,js:je), 0, 1 )     
    allocate(ebuffer(js:je), wbuffer(js:je), sbuffer(is:ie), nbuffer(is:ie))
    call mpp_get_boundary( dy, dx, domain, ebufferx=ebuffer, wbufferx=wbuffer, sbuffery=sbuffer, nbuffery=nbuffer,&
         flags=SCALAR_PAIR+XUPDATE, gridtype=CGRID_NE_PARAM)
    if(is == 1 .AND. mod(tile,2) .NE. 0) then ! on the west boundary
       dy(is, js:je) = wbuffer(js:je)
    endif
    if(ie == npx-1) then  ! on the east boundary
       dy(ie+1, js:je) = ebuffer(js:je)
    endif
    deallocate(wbuffer, ebuffer, sbuffer, nbuffer)

    call mpp_update_domains( dy, dx, domain, flags=SCALAR_PAIR,      &
         gridtype=CGRID_NE_PARAM, complete=.true.)

    call fill_corners(dx, dy, npx, npy, DGRID=.true.)

    !--- dxa and dya

    do j=jsd,jed
       do i=isd,ied
          !        do j=js,je
          !           do i=is,ie
          call mid_pt_sphere(grid(i,  j,1:2), grid(i,  j+1,1:2), p1)
          call mid_pt_sphere(grid(i+1,j,1:2), grid(i+1,j+1,1:2), p2)
          dxa(i,j) = great_circle_dist( p2, p1, radius )
          !
          call mid_pt_sphere(grid(i,j  ,1:2), grid(i+1,j  ,1:2), p1)
          call mid_pt_sphere(grid(i,j+1,1:2), grid(i+1,j+1,1:2), p2)
          dya(i,j) = great_circle_dist( p2, p1, radius )
       enddo
    enddo
    !      call mpp_update_domains( dxa, dya, domain, flags=SCALAR_PAIR, gridtype=AGRID_PARAM)
    call fill_corners(dxa, dya, npx, npy, AGRID=.true.)

    !--- dxc and dyc
    do j=js,je
       do i=is,ie+1
          p1(1) = agrid(i-1,j,1)
          p1(2) = agrid(i-1,j,2)
          p2(1) = agrid(i  ,j,1)
          p2(2) = agrid(i  ,j,2)
          dxc(i,j) = great_circle_dist( p2, p1, radius )
       enddo
    enddo
    do j=js,je+1
       do i=is,ie
          p1(1) = agrid(i,j-1,1)
          p1(2) = agrid(i,j-1,2)
          p2(1) = agrid(i,j  ,1)
          p2(2) = agrid(i,j  ,2)
          dyc(i,j) = great_circle_dist( p2, p1, radius )
       enddo
    enddo

    !--- area and area_c
    allocate (iinta(4, isd:ied ,jsd:jed), jinta(4, isd:ied ,jsd:jed),  &
              iintb(4, is:ie+1 ,js:je+1), jintb(4, is:ie+1 ,js:je+1))
    call sorted_inta(isd, ied, jsd, jed, cubed_sphere, grid, iinta, jinta)
    call sorted_intb(isd, ied, jsd, jed, is, ie, js, je, npx, npy, &
         cubed_sphere, agrid, iintb, jintb)
    call grid_area( npx, npy, ndims, nregions )
    deallocate(iintb, jintb)

#ifndef ORIG_AREA_C
    ! Compute area_c, rarea_c, dxc, dyc
    if ( is==1 ) then
       i = 1
       do j=js,je+1
          call mid_pt_sphere(grid(i,j-1,1:2), grid(i,j,  1:2), p1)
          call mid_pt_sphere(grid(i,j  ,1:2), grid(i,j+1,1:2), p4)
          p2(1:2) = agrid(i,j-1,1:2)
          p3(1:2) = agrid(i,j,  1:2)
          area_c(i,j) = 2.*get_area(p1, p4, p2, p3, radius)
       enddo
       do j=js,je
          call mid_pt_sphere(grid(i,j,1:2), grid(i,j+1,1:2), p1)
          p2(1:2) = agrid(i,j,1:2)
          dxc(i,j) = 2.*great_circle_dist( p1, p2, radius )
       enddo
    endif
    if ( (ie+1)==npx ) then
       i = npx
       do j=js,je+1
          p1(1:2) = agrid(i-1,j-1,1:2)
          call mid_pt_sphere(grid(i,j-1,1:2), grid(i,j,  1:2), p2)
          call mid_pt_sphere(grid(i,j  ,1:2), grid(i,j+1,1:2), p3)
          p4(1:2) = agrid(i-1,j,1:2)
          area_c(i,j) = 2.*get_area(p1, p4, p2, p3, radius)
       enddo
       do j=js,je
          p1(1:2) = agrid(i-1,j,1:2)
          call mid_pt_sphere(grid(i,j,1:2), grid(i,j+1,1:2), p2)
          dxc(i,j) = 2.*great_circle_dist( p1, p2, radius )
       enddo
    endif
    if ( js==1 ) then
       j = 1
       do i=is,ie+1
          call mid_pt_sphere(grid(i-1,j,1:2), grid(i,  j,1:2), p1)
          call mid_pt_sphere(grid(i,  j,1:2), grid(i+1,j,1:2), p2)
          p3(1:2) = agrid(i,  j,1:2)
          p4(1:2) = agrid(i-1,j,1:2)
          area_c(i,j) = 2.*get_area(p1, p4, p2, p3, radius)
       enddo
       do i=is,ie
          call mid_pt_sphere(grid(i,j,1:2), grid(i+1,j,1:2), p1)
          p2(1:2) = agrid(i,j,1:2)
          dyc(i,j) = 2.*great_circle_dist( p1, p2, radius )
       enddo
    endif
    if ( (je+1)==npy ) then
       j = npy
       do i=is,ie+1
          p1(1:2) = agrid(i-1,j-1,1:2)
          p2(1:2) = agrid(i  ,j-1,1:2)
          call mid_pt_sphere(grid(i,  j,1:2), grid(i+1,j,1:2), p3)
          call mid_pt_sphere(grid(i-1,j,1:2), grid(i,  j,1:2), p4)
          area_c(i,j) = 2.*get_area(p1, p4, p2, p3, radius)
       enddo
       do i=is,ie
          p1(1:2) = agrid(i,j-1,1:2)
          call mid_pt_sphere(grid(i,j,1:2), grid(i+1,j,1:2), p2)
          dyc(i,j) = 2.*great_circle_dist( p1, p2, radius )
       enddo
    endif
    if ( sw_corner ) then
       i=1; j=1
       p1(1:2) = grid(i,j,1:2)
       call mid_pt_sphere(grid(i,j,1:2), grid(i+1,j,1:2), p2)
       p3(1:2) = agrid(i,j,1:2)
       call mid_pt_sphere(grid(i,j,1:2), grid(i,j+1,1:2), p4)
       area_c(i,j) = 3.*get_area(p1, p4, p2, p3, radius)
    endif
    if ( se_corner ) then
       i=npx; j=1
       call mid_pt_sphere(grid(i-1,j,1:2), grid(i,j,1:2), p1)
       p2(1:2) = grid(i,j,1:2)
       call mid_pt_sphere(grid(i,j,1:2), grid(i,j+1,1:2), p3)
       p4(1:2) = agrid(i,j,1:2)
       area_c(i,j) = 3.*get_area(p1, p4, p2, p3, radius)
    endif
    if ( ne_corner ) then
       i=npx; j=npy
       p1(1:2) = agrid(i-1,j-1,1:2)
       call mid_pt_sphere(grid(i,j-1,1:2), grid(i,j,1:2), p2)
       p3(1:2) = grid(i,j,1:2)
       call mid_pt_sphere(grid(i-1,j,1:2), grid(i,j,1:2), p4)
       area_c(i,j) = 3.*get_area(p1, p4, p2, p3, radius)
    endif
    if ( nw_corner ) then
       i=1; j=npy
       call mid_pt_sphere(grid(i,j-1,1:2), grid(i,j,1:2), p1)
       p2(1:2) = agrid(i,j-1,1:2)
       call mid_pt_sphere(grid(i,j,1:2), grid(i+1,j,1:2), p3)
       p4(1:2) = grid(i,j,1:2)
       area_c(i,j) = 3.*get_area(p1, p4, p2, p3, radius)
    endif
#endif

    call mpp_update_domains( dxc, dyc, domain, flags=SCALAR_PAIR,   &
         gridtype=CGRID_NE_PARAM, complete=.true.)
    call fill_corners(dxc, dyc, npx, npy, CGRID=.true.)

    call mpp_update_domains( area,   domain, complete=.true. )
    call mpp_update_domains( area_c, domain, position=CORNER, complete=.true.)

    ! Handle corner Area ghosting
    call fill_ghost(area, npx, npy, -1.E35)  ! fill in garbage values
    call fill_corners(area_c, npx, npy, FILL=XDir, BGRID=.true.)

    do j=jsd,jed+1
       do i=isd,ied
          rdx(i,j) = 1.0/dx(i,j)
       enddo
    enddo
    do j=jsd,jed
       do i=isd,ied+1
          rdy(i,j) = 1.0/dy(i,j)
       enddo
    enddo
    do j=jsd,jed
       do i=isd,ied+1
          rdxc(i,j) = 1.0/dxc(i,j)
       enddo
    enddo
    do j=jsd,jed+1
       do i=isd,ied
          rdyc(i,j) = 1.0/dyc(i,j)
       enddo
    enddo
    do j=jsd,jed
       do i=isd,ied
          rarea(i,j) = 1.0/area(i,j)
          rdxa(i,j) = 1./dxa(i,j)
          rdya(i,j) = 1./dya(i,j)
       enddo
    enddo
    do j=jsd,jed+1
       do i=isd,ied+1
          rarea_c(i,j) = 1.0/area_c(i,j)
       enddo
    enddo

200    format(A,f9.2,A,f9.2,A,f9.2)
201    format(A,f9.2,A,f9.2,A,f9.2,A,f9.2)
202    format(A,A,i4.4,A,i4.4,A)

    ! Get and print Grid Statistics, Only from tile 1
    dxAV =0.0
    angAV=0.0
    aspAV=0.0
    dxN  =  missing
    dxM  = -missing
    angN =  missing
    angM = -missing
    aspN =  missing
    aspM = -missing
    allocate(angs(is:ie,js:je), asps(is:ie,js:je), dxs(is:ie,js:je) )
    if (tile == 1) then
       do j=js, je
          do i=is, ie
             if(i>ceiling(npx/2.) .OR. j>ceiling(npy/2.)) cycle
             ang  = get_angle(2, grid(i,j+1,1:2), grid(i,j,1:2), grid(i+1,j,1:2))
             ang  = ABS(90.0 - ang)
             angs(i,j) = ang

             if ( (i==1) .and. (j==1) ) then
             else 
                angAV = angAV + ang
                angM  = MAX(angM,ang)
                angN  = MIN(angN,ang)
             endif

             dx_local = dx(i,j)
             dy_local = dy(i,j)

             dxAV  = dxAV + 0.5 * (dx_local + dy_local)
             dxM   = MAX(dxM,dx_local)
             dxM   = MAX(dxM,dy_local)
             dxN   = MIN(dxN,dx_local)
             dxN   = MIN(dxN,dy_local)
             dxs(i,j) = dy_local !0.5 * (dx_local + dy_local)

             asp   = ABS(dx_local/dy_local)
             if (asp < 1.0) asp = 1.0/asp
             asps(i,j) = asp 
             aspAV = aspAV + asp
             aspM  = MAX(aspM,asp)
             aspN  = MIN(aspN,asp)
          enddo
       enddo
    else
       angs = 0
       asps = 0
       dxs  = 0 
    endif
    call mpp_sum(angAv)
    call mpp_sum(dxAV)
    call mpp_sum(aspAV)
    call mpp_max(angM)
    call mpp_min(angN)
    call mpp_max(dxM)
    call mpp_min(dxN)
    call mpp_max(aspM)
    call mpp_min(aspN)

    if( gid==masterproc ) then
       angAV = angAV / ( (ceiling(npy/2.0))*(ceiling(npx/2.0)) - 1 )
       dxAV  = dxAV  / ( (ceiling(npy/2.0))*(ceiling(npx/2.0)) )
       aspAV = aspAV / ( (ceiling(npy/2.0))*(ceiling(npx/2.0)) )
       write(*,*  ) ''
       write(*,*  ) ' Cubed-Sphere Grid Stats : ', npx,'x',npy,'x',nregions
       write(*,201) '      Grid Length               : min: ', dxN,' max: ', dxM,' avg: ', dxAV, ' min/max: ',dxN/dxM
       write(*,200) '      Deviation from Orthogonal : min: ',angN,' max: ',angM,' avg: ',angAV
       write(*,200) '      Aspect Ratio              : min: ',aspN,' max: ',aspM,' avg: ',aspAV
       write(*,*  ) ''
    endif

    if(write_grid_char_file) then
       allocate(g_tmp(npx-1,npy-1))
       call mpp_global_field(domain, angs, g_tmp)
       if( gid==masterproc ) then
          write(gcharFile,202) TRIM(grid_name),'_chars_',npx,'x',npy,'.dat'
          fileLun=get_unit()
          open(unit=fileLun,file=gcharFile, form='unformatted', access='direct',  &
               recl=((npx/2)+1)*((npy/2)+1)*8, status='unknown')
          allocate(tmp(1:(npx/2)+1, 1:(npy/2)+1))
          do j = 1,ceiling(npy/2.)
             do i=1,ceiling(npx/2.)
                tmp(i,j) = g_tmp(i,j)
             enddo
          enddo
          write(fileLun,rec=1) tmp
       endif

       call mpp_global_field(domain, asps, g_tmp)
       if( gid==masterproc ) then
          do j = 1,ceiling(npy/2.)
             do i=1,ceiling(npx/2.)
                tmp(i,j) = g_tmp(i,j)
             enddo
          enddo
          write(fileLun,rec=2) tmp
       endif

       call mpp_global_field(domain, dxs,  g_tmp)
       if( gid==masterproc ) then
          do j = 1,ceiling(npy/2.)
             do i=1,ceiling(npx/2.)
                tmp(i,j) = g_tmp(i,j)
             enddo
          enddo
          write(fileLun,rec=3) tmp
       endif

       if(tile == 1) then
          do j=js, je
             do i=is, ie
                if(i>(npx/2.0)+1 .OR. j>(npy/2.0)+1) cycle
                do n=1,ndims
                   p_lL(n) = grid(i  ,j  ,n)
                   p_uL(n) = grid(i  ,j+1,n)
                   p_lR(n) = grid(i+1,j  ,n)
                   p_uR(n) = grid(i+1,j+1,n)
                enddo
                if ((latlon) .or. (dxdy_area)) then
                   ! DX_*DY_
                   d1 = dx(i  ,j  )
                   d2 = dx(i  ,j+1)
                   mydx = 0.5 * ( d1+d2 )
                   d1 = dy(i  ,j)
                   d2 = dy(i+1,j)
                   mydy = 0.5 * ( d1+d2 )
                   angs(i,j) = (mydx*mydy)
                else
                   ! Spherical Excess Formula
                   angs(i,j) = get_area(p_lL, p_uL, p_lR, p_uR, radius)
                endif
             enddo
          enddo
       else
          angs = 0
       endif
       call mpp_global_field(domain, angs,  g_tmp)
       if( gid==masterproc ) then
          do j = 1,npy/2+1
             do i=1,npx/2+1
                tmp(i,j) = g_tmp(i,j)
             enddo
          enddo
          write(fileLun,rec=4) tmp
          close(unit=fileLun)
          deallocate(tmp ) 
       endif
       deallocate(angs, asps, dxs, g_tmp)
    endif

#ifdef GLOBAL_TRIG
    call mpp_error(FATAL, 'fv_grid_tools(read_grid): when reading from '// &
         trim(grid_file)//', -DGLOBAL_TRIG should not be present when compiling')
#endif

  end subroutine read_grid



  !#################################################################################
  subroutine get_symmetry(data_in, data_out, ishift, jshift)
    integer,                                            intent(in)  :: ishift, jshift
    real, dimension(is:ie+ishift, js:je+jshift ), intent(in)  :: data_in
    real, dimension(is:ie+jshift,js:je+ishift  ), intent(out) :: data_out      
    real,    dimension(:), allocatable :: send_buffer
    real,    dimension(:), allocatable :: recv_buffer
    integer, dimension(:), allocatable :: is_recv, ie_recv, js_recv, je_recv, pe_recv
    integer, dimension(:), allocatable :: is_send, ie_send, js_send, je_send, pe_send
    integer, dimension(:), allocatable :: isl, iel, jsl, jel, pelist, msg1, msg2
    integer                            :: msgsize, pos, ntiles, npes_per_tile, npes
    integer                            :: send_buf_size, recv_buf_size, buffer_pos
    integer                            :: is0, ie0, js0, je0
    integer                            :: is1, ie1, js1, je1
    integer                            :: is2, ie2, js2, je2
    integer                            :: i, j, p, nrecv, nsend, tile_you, is3, ie3, nlist
    integer                            :: start_pe, ipos, jpos, from_pe, to_pe
    
    !--- This routine will be called only for cubic sphere grid. so 6 tiles will be assumed
    !--- also number of processors on each tile will be the same.
    ntiles = mpp_get_ntile_count(domain)
    npes = mpp_npes()

    if(ntiles .NE. 6 ) call mpp_error(FATAL, 'fv_grid_tools(get_symmetry): ntiles should be 6 ')
    if(mod(npes,ntiles) /= 0) call mpp_error(FATAL, 'fv_grid_tools(get_symmetry): npes should be divided by ntiles')
    npes_per_tile = npes/ntiles

!   if(npes_x == npes_y) then ! even, simple communication
    if(npes_x == npes_y .AND. mod(npx_g-1,npes_x) == 0 ) then ! even, 
       msgsize = (ie-is+1+jshift)*(je-js+1+ishift)

       pos = mod((mpp_pe()-mpp_root_pe()), npes_x*npes_y)       
       start_pe = mpp_pe() - pos
       ipos = mod(pos, npes_x)
       jpos = pos/npes_x
       from_pe = start_pe + ipos*npes_x + jpos
       to_pe   = from_pe
       allocate(recv_buffer(msgsize))
       call mpp_recv(recv_buffer(1), glen=msgsize, from_pe=from_pe, block=.FALSE. )

       pos = 0
       allocate(send_buffer(msgsize))
       do j = js, je+jshift
          do i = is, ie+ishift
             pos = pos + 1
             send_buffer(pos) = data_in(i,j)
          enddo
       enddo

       call mpp_send(send_buffer(1), plen=msgsize, to_pe=to_pe)
       call mpp_sync_self(check=EVENT_RECV) ! To ensure recv is completed.

       !--unpack buffer
       pos = 0
       do i = is, ie+jshift
          do j = js, je+ishift
             pos = pos + 1
             data_out(i,j) = recv_buffer(pos)
          enddo
       enddo

       call mpp_sync_self()     
       deallocate(send_buffer, recv_buffer)
    else

       allocate(is_recv(0:npes_per_tile-1), ie_recv(0:npes_per_tile-1))
       allocate(js_recv(0:npes_per_tile-1), je_recv(0:npes_per_tile-1))
       allocate(is_send(0:npes_per_tile-1), ie_send(0:npes_per_tile-1))
       allocate(js_send(0:npes_per_tile-1), je_send(0:npes_per_tile-1))
       allocate(pe_send(0:npes_per_tile-1), pe_recv(0:npes_per_tile-1))
       if(debug_message_size) then
          allocate(msg1   (0:npes_per_tile-1), msg2   (0:npes_per_tile-1))
          msg1 = 0
          msg2 = 0
       endif

       allocate(pelist(0:npes-1))
       call mpp_get_pelist(domain, pelist)
       allocate(isl(0:npes-1), iel(0:npes-1), jsl(0:npes-1), jel(0:npes-1) )
       call mpp_get_compute_domains(domain, xbegin=isl, xend=iel, ybegin=jsl, yend=jel)
       !--- pre-post receiving 
       buffer_pos = 0  
       nrecv = 0
       nsend = 0
       recv_buf_size = 0

       !--- first set up the receiving index
       nlist = 0
       do p = 0, npes-1
          tile_you = p/(npes_x*npes_y) + 1
          if(tile_you .NE. tile) cycle

          !--- my index for data_out after rotation
          is1 = js; ie1 = je + ishift;
          js1 = is; je1 = ie + jshift;
          !--- your index for data_out
          is2 = isl(p); ie2 = iel(p) + ishift;
          js2 = jsl(p); je2 = jel(p) + jshift;
          is0 = max(is1,is2); ie0 = min(ie1,ie2)
          js0 = max(js1,js2); je0 = min(je1,je2)             
          msgsize = 0             
          if(ie0 .GE. is0 .AND. je0 .GE. js0) then
             msgsize = (ie0-is0+1)*(je0-js0+1)
             recv_buf_size = recv_buf_size + msgsize
             pe_recv(nrecv) = pelist(p)
             !--- need to rotate back the index
             is_recv(nrecv) = js0; ie_recv(nrecv) = je0
             js_recv(nrecv) = is0; je_recv(nrecv) = ie0
             nrecv = nrecv+1
          endif
          if(debug_message_size) then
             msg1(nlist) = msgsize
             call mpp_recv(msg2(nlist), glen=1, from_pe=pelist(p), block=.FALSE. )
             nlist = nlist + 1
          endif
       enddo

       !--- Then setup the sending index.
       send_buf_size = 0
       do p = 0, npes-1
          tile_you = p/(npes_x*npes_y) + 1
          if(tile_you .NE. tile) cycle
          !--- my index on data_in
          is1 = is; ie1 = ie + ishift;
          js1 = js; je1 = je + jshift;
          !--- your index on data_out after rotate
          is2 = jsl(p); ie2 = jel(p) + ishift;
          js2 = isl(p); je2 = iel(p) + jshift;
          is0 = max(is1,is2); ie0 = min(ie1,ie2)
          js0 = max(js1,js2); je0 = min(je1,je2)
          msgsize = 0
          if(ie0 .GE. is0 .AND. je0 .GE. js0 )then
             msgsize = (ie0-is0+1)*(je0-js0+1)
             send_buf_size = send_buf_size + msgsize
             pe_send(nsend) = pelist(p)
             is_send(nsend) = is0; ie_send(nsend) = ie0
             js_send(nsend) = js0; je_send(nsend) = je0
             nsend = nsend+1
          endif
          IF(debug_message_size) call mpp_send(msgsize, plen=1, to_pe=pelist(p) )
       enddo

       !--- check to make sure send and recv size match.
       if(debug_message_size) then
          call mpp_sync_self(check=EVENT_RECV) ! To ensure recv is completed.
          do p = 0, nlist-1
             if(msg1(p) .NE. msg2(p)) then
                call mpp_error(FATAL, "fv_grid_tools_mod(get_symmetry): mismatch on send and recv size")
             endif
          enddo
          call mpp_sync_self()
          deallocate(msg1, msg2)
       endif

       !--- pre-post data
       allocate(recv_buffer(recv_buf_size))
       buffer_pos = 0
       do p = 0, nrecv-1
          is0 = is_recv(p); ie0 = ie_recv(p)
          js0 = js_recv(p); je0 = je_recv(p)
          msgsize = (ie0-is0+1)*(je0-js0+1)
          call mpp_recv(recv_buffer(buffer_pos+1), glen=msgsize, from_pe=pe_recv(p), block=.FALSE. )
          buffer_pos = buffer_pos + msgsize       
       enddo

       !--- send the data
       buffer_pos = 0
       allocate(send_buffer(send_buf_size))
       do p = 0, nsend-1
          is0 = is_send(p); ie0 = ie_send(p)
          js0 = js_send(p); je0 = je_send(p)
          msgsize = (ie0-is0+1)*(je0-js0+1)
          pos = buffer_pos
          do j = js0, je0
             do i = is0, ie0
                pos = pos+1
                send_buffer(pos) = data_in(i,j)
             enddo
          enddo
          call mpp_send(send_buffer(buffer_pos+1), plen=msgsize, to_pe=pe_send(p) )
          buffer_pos = buffer_pos + msgsize       
       enddo

       call mpp_sync_self(check=EVENT_RECV) ! To ensure recv is completed.

       !--- unpack buffer
       pos = 0
       do p = 0, nrecv-1
          is0 = is_recv(p); ie0 = ie_recv(p)       
          js0 = js_recv(p); je0 = je_recv(p)

          do i = is0, ie0
             do j = js0, je0
                pos = pos + 1
                data_out(i,j) = recv_buffer(pos)
             enddo
          enddo
       enddo

       call mpp_sync_self()
       deallocate(isl, iel, jsl, jel, pelist)
       deallocate(is_recv, ie_recv, js_recv, je_recv, pe_recv)
       deallocate(is_send, ie_send, js_send, je_send, pe_send)
       deallocate(recv_buffer, send_buffer)
     endif

  end subroutine get_symmetry

  subroutine init_grid(Atm, grid_name, grid_file, npx, npy, npz, ndims, nregions, ng)
 
!     init_grid :: read grid from input file and setup grid descriptors
 
!--------------------------------------------------------
    type(fv_atmos_type), intent(inout) :: Atm
    character(len=80), intent(IN) :: grid_name
    character(len=120),intent(IN) :: grid_file
    integer,      intent(IN) :: npx, npy, npz
    integer,      intent(IN) :: ndims
    integer,      intent(IN) :: nregions
    integer,      intent(IN) :: ng
!--------------------------------------------------------
    real   ::  xs(npx,npy)
    real   ::  ys(npx,npy)
    real(kind=8) ::  grid_R8(npx,npy)

    real  :: dp, dl
    real  :: x1,x2,y1,y2,z1,z2
    integer :: i,j,k,n,nreg
    integer :: fileLun

    real  :: p_lL(ndims) ! lower Left
    real  :: p_uL(ndims) ! upper Left
    real  :: p_lR(ndims) ! lower Right
    real  :: p_uR(ndims) ! upper Right
    real  :: d1, d2, mydx, mydy, tmp

    real  :: p1(3), p2(3), p3(3), p4(3)
    real  :: dist,dist1,dist2, pa(2), pa1(2), pa2(2), pb(2)
    real  :: pt(3), pt1(3), pt2(3), pt3(3)
    real :: ee1(3), ee2(3)

    real  :: angN,angM,angAV,ang
    real  :: aspN,aspM,aspAV,asp
    real  ::  dxN, dxM, dxAV
    real  :: dx_local, dy_local

    real  :: vec1(3), vec2(3), vec3(3), vec4(3)
    real  :: vecAvg(3), vec3a(3), vec3b(3), vec4a(3), vec4b(3)
    real  :: xyz1(3), xyz2(3)

    real  :: angs(1:(npx/2)+1, 1:(npy/2)+1)
    real  :: asps(1:(npx/2)+1, 1:(npy/2)+1)
    real  ::  dxs(1:(npx/2)+1, 1:(npy/2)+1)
    character(len=80) :: gcharFile

    real :: grid_global(1-ng:npx  +ng,1-ng:npy  +ng,ndims,1:nregions)
    real ::   dx_global(1:npx-1,1:npy  ,1:nregions)
    real ::   dy_global(1:npx  ,1:npy-1,1:nregions)

    character(len=80) :: evalue
    integer :: ios, ip, jp
    
    integer :: igrid
    
    integer :: tmplun
    character(len=80) :: tmpFile   

    npx_g = npx
    npy_g = npy
    npz_g = npz
    ntiles_g = nregions
    latlon = .false.
    cubed_sphere = .false.
    if ( grid_type < 0 ) then
       Gnomonic_grid = .false.
    else
       Gnomonic_grid = .true.
    endif

    if ( Atm%do_schmidt .and. abs(atm%stretch_fac-1.) > 1.E-5 ) stretched_grid = .true.
    
    allocate (  area(isd:ied  ,jsd:jed  ) )   ! Cell Centered
    allocate ( rarea(isd:ied  ,jsd:jed  ) )   ! Cell Centered
    
    allocate (  area_c(isd:ied+1,jsd:jed+1) )  ! Cell Corners
    allocate ( rarea_c(isd:ied+1,jsd:jed+1) )  ! Cell Corners
    
    allocate (  dx(isd:ied  ,jsd:jed+1) )
    allocate ( rdx(isd:ied  ,jsd:jed+1) )
    allocate (  dy(isd:ied+1,jsd:jed  ) )
    allocate ( rdy(isd:ied+1,jsd:jed  ) )
    
    allocate (  dxc(isd:ied+1,jsd:jed  ) )
    allocate ( rdxc(isd:ied+1,jsd:jed  ) )
    allocate (  dyc(isd:ied  ,jsd:jed+1) )
    allocate ( rdyc(isd:ied  ,jsd:jed+1) )
    
    allocate (  dxa(isd:ied  ,jsd:jed  ) )
    allocate ( rdxa(isd:ied  ,jsd:jed  ) )
    allocate (  dya(isd:ied  ,jsd:jed  ) )
    allocate ( rdya(isd:ied  ,jsd:jed  ) )
    
    allocate ( grid (isd:ied+1,jsd:jed+1,1:ndims) )
    allocate ( agrid(isd:ied  ,jsd:jed  ,1:ndims) )
    
#ifndef NO_GRID_G
    allocate ( grid_g(1:npx,1:npy,1:ndims) )
    Atm%grid_g =>grid_g
#endif
    Atm%grid  =>grid
    Atm%agrid =>agrid
    
    allocate ( sina(isd:ied+1,jsd:jed+1) )   ! SIN(angle of intersection)
    allocate ( cosa(isd:ied+1,jsd:jed+1) )   ! COS(angle of intersection)
    
    allocate (   e1(3,isd:ied+1,jsd:jed+1) )
    allocate (   e2(3,isd:ied+1,jsd:jed+1) )
    
    if( .not. stretched_grid )      &
    allocate (iinta(4, isd:ied ,jsd:jed), jinta(4, isd:ied ,jsd:jed),  &
              iintb(4, is:ie+1 ,js:je+1), jintb(4, is:ie+1 ,js:je+1))

    if (grid_type>3) then
       uniform_ppm = .true.
       if (grid_type == 4) then
          double_periodic = .true.
          call setup_cartesian(npx, npy)
       else
          call mpp_error(FATAL, 'init_grid: unsupported grid type')
       endif
    else

          cubed_sphere = .true.
          uniform_ppm = .true.

          if (grid_type>=0) call gnomonic_grids(grid_type, npx-1, xs, ys)

          if (gid == masterproc) then

             if (grid_type>=0) then
                do j=1,npy
                   do i=1,npx
                      grid_global(i,j,1,1) = xs(i,j)
                      grid_global(i,j,2,1) = ys(i,j)
                   enddo
                enddo
! mirror_grid assumes that the tile=1 is centered on equator and greenwich meridian Lon[-pi,pi] 
                call mirror_grid(grid_global, ng, npx, npy, 2, 6)
                do n=1,nregions
                   do j=1,npy
                      do i=1,npx
!---------------------------------
! Shift the corner away from Japan
!---------------------------------
!--------------------- This will result in the corner close to east coast of China ------------------
                         if ( .not.Atm%do_schmidt .and. (Atm%shift_fac)>1.E-4 )   &
                              grid_global(i,j,1,n) = grid_global(i,j,1,n) - pi/Atm%shift_fac
!----------------------------------------------------------------------------------------------------
                         if ( grid_global(i,j,1,n) < 0. )              &
                              grid_global(i,j,1,n) = grid_global(i,j,1,n) + 2.*pi
                         if (ABS(grid_global(i,j,1,1)) < 1.e-10) grid_global(i,j,1,1) = 0.0
                         if (ABS(grid_global(i,j,2,1)) < 1.e-10) grid_global(i,j,2,1) = 0.0
                      enddo
                   enddo
                enddo
             else
                call mpp_error(FATAL, "fv_grid_tools: reading of ASCII grid files no longer supported")
             endif

             grid_global(  1,1:npy,:,2)=grid_global(npx,1:npy,:,1)
             grid_global(  1,1:npy,:,3)=grid_global(npx:1:-1,npy,:,1)
             grid_global(1:npx,npy,:,5)=grid_global(1,npy:1:-1,:,1)
             grid_global(1:npx,npy,:,6)=grid_global(1:npx,1,:,1)
             
             grid_global(1:npx,  1,:,3)=grid_global(1:npx,npy,:,2)
             grid_global(1:npx,  1,:,4)=grid_global(npx,npy:1:-1,:,2)
             grid_global(npx,1:npy,:,6)=grid_global(npx:1:-1,1,:,2)
             
             grid_global(  1,1:npy,:,4)=grid_global(npx,1:npy,:,3)
             grid_global(  1,1:npy,:,5)=grid_global(npx:1:-1,npy,:,3)
             
             grid_global(npx,1:npy,:,3)=grid_global(1,1:npy,:,4)
             grid_global(1:npx,  1,:,5)=grid_global(1:npx,npy,:,4)
             grid_global(1:npx,  1,:,6)=grid_global(npx,npy:1:-1,:,4)
             
             grid_global(  1,1:npy,:,6)=grid_global(npx,1:npy,:,5)

!------------------------
! Schmidt transformation:
!------------------------
             if ( Atm%do_schmidt ) then
             do n=1,nregions
                call direct_transform(Atm%stretch_fac, 1, npx, 1, npy, Atm%target_lon, Atm%target_lat, &
                                      n, grid_global(1:npx,1:npy,1,n), grid_global(1:npx,1:npy,2,n))
             enddo
             endif

! Compute dx_global:
             do n=1,nregions
                do j=1,npy
                   do i=1,npx-1
                      dx_global(i,j,n) = great_circle_dist(grid_global(i,j,:,n), grid_global(i+1,j,:,n), radius)
                   enddo
                enddo
             enddo
             dx_global(1:npx-1,  1,1) = dx_global(1:npx-1,npy,6)
             dx_global(1:npx-1,npy,2) = dx_global(1:npx-1,1,  3)
             dx_global(1:npx-1,npy,4) = dx_global(1:npx-1,1,  5)


! Compute dy_global:
        if( stretched_grid ) then
           do n=1,nregions
                do j=1,npy-1
                   do i=1,npx
                      dy_global(i,j,n) = great_circle_dist(grid_global(i,j,:,n), grid_global(i,j+1,:,n), radius)
                   enddo
                enddo
             enddo
        else
! Grid symmetry is assumed here: for non-stretched grids
           do n=1,nregions
              do j=1,npy-1
                 do i=1,npx-1
                    dy_global(j,i,n) = dx_global(i,j,n)
                 enddo
              enddo
           enddo
        endif
           dy_global(npx,1:npy-1,1)=dy_global(1,1:npy-1,2)
           dy_global(npx,1:npy-1,3)=dy_global(1,1:npy-1,4)
           dy_global(npx,1:npy-1,5)=dy_global(1,1:npy-1,6)

           dy_global(  1,1:npy-1,3)=dx_global(npx-1:1:-1, npy,1)
           dy_global(npx,1:npy-1,6)=dx_global(npx-1:1:-1, 1,  2)
           dy_global(  1,1:npy-1,5)=dx_global(npx-1:1:-1, npy,3)
           dy_global(npx,1:npy-1,2)=dx_global(npx-1:1:-1, 1,  4)
           dy_global(  1,1:npy-1,1)=dx_global(npx-1:1:-1, npy,5)
           dy_global(npx,1:npy-1,4)=dx_global(npx-1:1:-1, 1,  6)
             
          endif ! masterproc

       call mpp_broadcast(grid_global, size(grid_global), masterproc)
       call mpp_broadcast(dx_global, size(dx_global), masterproc)
       call mpp_broadcast(dy_global, size(dy_global), masterproc)
      
       do n=1,ndims
          do j=js,je+1
             do i=is,ie+1
                grid(i,j,n) = grid_global(i,j,n,tile)
             enddo
          enddo
       enddo
!
! SJL: For phys/exchange grid, etc
!
#ifndef NO_GRID_G
       do j=1,npy
          do i=1,npx
             grid_g(i,j,1) = grid_global(i,j,1,tile)
             grid_g(i,j,2) = grid_global(i,j,2,tile)
          enddo
       enddo
#endif
       call mpp_update_domains( grid, domain, position=CORNER)
       call fill_corners(grid(:,:,1), npx, npy, FILL=XDir, BGRID=.true.)
       call fill_corners(grid(:,:,2), npx, npy, FILL=XDir, BGRID=.true.)

       if( .not. stretched_grid )         &
       call sorted_inta(isd, ied, jsd, jed, cubed_sphere, grid, iinta, jinta)

       agrid(:,:,:) = -1.e25
 
       do j=js,je
          do i=is,ie
             if ( stretched_grid ) then
                  call cell_center2(grid(i,j,  1:2), grid(i+1,j,  1:2),   &
                                    grid(i,j+1,1:2), grid(i+1,j+1,1:2),   &
                                    agrid(i,j,1:2) )
             else
                  call cell_center2(grid(iinta(1,i,j),jinta(1,i,j),1:2),  &
                                    grid(iinta(2,i,j),jinta(2,i,j),1:2),  &
                                    grid(iinta(3,i,j),jinta(3,i,j),1:2),  &
                                    grid(iinta(4,i,j),jinta(4,i,j),1:2),  &
                                    agrid(i,j,1:2) )
             endif
          enddo
       enddo

       call mpp_update_domains( agrid, domain, position=CENTER, complete=.true. )
       call fill_corners(agrid(:,:,1), npx, npy, XDir, AGRID=.true.)
       call fill_corners(agrid(:,:,2), npx, npy, YDir, AGRID=.true.)

       do j=js,je+1
          do i=is,ie
             dx(i,j) = dx_global(i,j,tile)
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             dy(i,j) = dy_global(i,j,tile)
          enddo
       enddo

       call mpp_update_domains( dy, dx, domain, flags=SCALAR_PAIR,      &
                                gridtype=CGRID_NE_PARAM, complete=.true.)
       if (cubed_sphere) call fill_corners(dx, dy, npx, npy, DGRID=.true.)

       do j=jsd,jed
          do i=isd,ied
             call mid_pt_sphere(grid(i,  j,1:2), grid(i,  j+1,1:2), p1)
             call mid_pt_sphere(grid(i+1,j,1:2), grid(i+1,j+1,1:2), p2)
             dxa(i,j) = great_circle_dist( p2, p1, radius )
!
             call mid_pt_sphere(grid(i,j  ,1:2), grid(i+1,j  ,1:2), p1)
             call mid_pt_sphere(grid(i,j+1,1:2), grid(i+1,j+1,1:2), p2)
             dya(i,j) = great_circle_dist( p2, p1, radius )
          enddo
       enddo
!      call mpp_update_domains( dxa, dya, domain, flags=SCALAR_PAIR, gridtype=AGRID_PARAM)
       if (cubed_sphere) call fill_corners(dxa, dya, npx, npy, AGRID=.true.)

       do j=js,je
          do i=is,ie+1
             dxc(i,j) = great_circle_dist(agrid(i,j,:), agrid(i-1,j,:), radius)
          enddo
       enddo
       do j=js,je+1
          do i=is,ie
             dyc(i,j) = great_circle_dist(agrid(i,j,:), agrid(i,j-1,:), radius)
          enddo
       enddo

       if( .not. stretched_grid )      &
       call sorted_intb(isd, ied, jsd, jed, is, ie, js, je, npx, npy, &
                        cubed_sphere, agrid, iintb, jintb)

       call grid_area( npx, npy, ndims, nregions )
!      stretched_grid = .false.

!----------------------------------
! Compute area_c, rarea_c, dxc, dyc
!----------------------------------
  if ( .not. stretched_grid ) then
! For symmetrical grids:
       if ( is==1 ) then
          i = 1
          do j=js,je+1
             call mid_pt_sphere(grid(i,j-1,1:2), grid(i,j,  1:2), p1)
             call mid_pt_sphere(grid(i,j  ,1:2), grid(i,j+1,1:2), p4)
             p2(1:2) = agrid(i,j-1,1:2)
             p3(1:2) = agrid(i,j,  1:2)
             area_c(i,j) = 2.*get_area(p1, p4, p2, p3, radius)
          enddo
          do j=js,je
             call mid_pt_sphere(grid(i,j,1:2), grid(i,j+1,1:2), p1)
             p2(1:2) = agrid(i,j,1:2)
             dxc(i,j) = 2.*great_circle_dist( p1, p2, radius )
          enddo
       endif
       if ( (ie+1)==npx ) then
          i = npx
          do j=js,je+1
             p1(1:2) = agrid(i-1,j-1,1:2)
             call mid_pt_sphere(grid(i,j-1,1:2), grid(i,j,  1:2), p2)
             call mid_pt_sphere(grid(i,j  ,1:2), grid(i,j+1,1:2), p3)
             p4(1:2) = agrid(i-1,j,1:2)
             area_c(i,j) = 2.*get_area(p1, p4, p2, p3, radius)
          enddo
          do j=js,je
             p1(1:2) = agrid(i-1,j,1:2)
             call mid_pt_sphere(grid(i,j,1:2), grid(i,j+1,1:2), p2)
             dxc(i,j) = 2.*great_circle_dist( p1, p2, radius )
          enddo
       endif
       if ( js==1 ) then
          j = 1
          do i=is,ie+1
             call mid_pt_sphere(grid(i-1,j,1:2), grid(i,  j,1:2), p1)
             call mid_pt_sphere(grid(i,  j,1:2), grid(i+1,j,1:2), p2)
             p3(1:2) = agrid(i,  j,1:2)
             p4(1:2) = agrid(i-1,j,1:2)
             area_c(i,j) = 2.*get_area(p1, p4, p2, p3, radius)
          enddo
          do i=is,ie
             call mid_pt_sphere(grid(i,j,1:2), grid(i+1,j,1:2), p1)
             p2(1:2) = agrid(i,j,1:2)
             dyc(i,j) = 2.*great_circle_dist( p1, p2, radius )
          enddo
       endif
       if ( (je+1)==npy ) then
          j = npy
          do i=is,ie+1
             p1(1:2) = agrid(i-1,j-1,1:2)
             p2(1:2) = agrid(i  ,j-1,1:2)
             call mid_pt_sphere(grid(i,  j,1:2), grid(i+1,j,1:2), p3)
             call mid_pt_sphere(grid(i-1,j,1:2), grid(i,  j,1:2), p4)
             area_c(i,j) = 2.*get_area(p1, p4, p2, p3, radius)
          enddo
          do i=is,ie
             p1(1:2) = agrid(i,j-1,1:2)
             call mid_pt_sphere(grid(i,j,1:2), grid(i+1,j,1:2), p2)
             dyc(i,j) = 2.*great_circle_dist( p1, p2, radius )
          enddo
       endif

       if ( sw_corner ) then
             i=1; j=1
             p1(1:2) = grid(i,j,1:2)
             call mid_pt_sphere(grid(i,j,1:2), grid(i+1,j,1:2), p2)
             p3(1:2) = agrid(i,j,1:2)
             call mid_pt_sphere(grid(i,j,1:2), grid(i,j+1,1:2), p4)
             area_c(i,j) = 3.*get_area(p1, p4, p2, p3, radius)
       endif
       if ( se_corner ) then
             i=npx; j=1
             call mid_pt_sphere(grid(i-1,j,1:2), grid(i,j,1:2), p1)
             p2(1:2) = grid(i,j,1:2)
             call mid_pt_sphere(grid(i,j,1:2), grid(i,j+1,1:2), p3)
             p4(1:2) = agrid(i,j,1:2)
             area_c(i,j) = 3.*get_area(p1, p4, p2, p3, radius)
       endif
       if ( ne_corner ) then
             i=npx; j=npy
             p1(1:2) = agrid(i-1,j-1,1:2)
             call mid_pt_sphere(grid(i,j-1,1:2), grid(i,j,1:2), p2)
             p3(1:2) = grid(i,j,1:2)
             call mid_pt_sphere(grid(i-1,j,1:2), grid(i,j,1:2), p4)
             area_c(i,j) = 3.*get_area(p1, p4, p2, p3, radius)
       endif
       if ( nw_corner ) then
             i=1; j=npy
             call mid_pt_sphere(grid(i,j-1,1:2), grid(i,j,1:2), p1)
             p2(1:2) = agrid(i,j-1,1:2)
             call mid_pt_sphere(grid(i,j,1:2), grid(i+1,j,1:2), p3)
             p4(1:2) = grid(i,j,1:2)
             area_c(i,j) = 3.*get_area(p1, p4, p2, p3, radius)
       endif
   endif
!-----------------

       call mpp_update_domains( dxc, dyc, domain, flags=SCALAR_PAIR,   &
                                gridtype=CGRID_NE_PARAM, complete=.true.)
       if (cubed_sphere) call fill_corners(dxc, dyc, npx, npy, CGRID=.true.)
       
       call mpp_update_domains( area,   domain, complete=.true. )
       call mpp_update_domains( area_c, domain, position=CORNER, complete=.true.)
       
       ! Handle corner Area ghosting
       if (cubed_sphere) then
          call fill_ghost(area, npx, npy, -1.E35)  ! fill in garbage values
          call fill_corners(area_c, npx, npy, FILL=XDir, BGRID=.true.)
       endif
       
       do j=jsd,jed+1
          do i=isd,ied
             rdx(i,j) = 1.0/dx(i,j)
          enddo
       enddo
       do j=jsd,jed
          do i=isd,ied+1
             rdy(i,j) = 1.0/dy(i,j)
          enddo
       enddo
       do j=jsd,jed
          do i=isd,ied+1
             rdxc(i,j) = 1.0/dxc(i,j)
          enddo
       enddo
       do j=jsd,jed+1
          do i=isd,ied
             rdyc(i,j) = 1.0/dyc(i,j)
          enddo
       enddo
       do j=jsd,jed
          do i=isd,ied
             rarea(i,j) = 1.0/area(i,j)
             rdxa(i,j) = 1./dxa(i,j)
             rdya(i,j) = 1./dya(i,j)
          enddo
       enddo
       do j=jsd,jed+1
          do i=isd,ied+1
             rarea_c(i,j) = 1.0/area_c(i,j)
          enddo
       enddo

200    format(A,f9.2,A,f9.2,A,f9.2)
201    format(A,f9.2,A,f9.2,A,f9.2,A,f9.2)
202    format(A,A,i4.4,A,i4.4,A)
       
! Get and print Grid Statistics
       if ((gid==masterproc) .and. (cubed_sphere)) then
          dxAV =0.0
          angAV=0.0
          aspAV=0.0
          dxN  =  missing
          dxM  = -missing
          angN =  missing
          angM = -missing
          aspN =  missing
          aspM = -missing
          angs(1,1) = get_angle(2, grid_global(1,2,1:2,1), grid_global(1,1,1:2,1), grid_global(2,1,1:2,1))
          angs(1,1) = ABS(90.0 - angs(1,1))
          do j=1,ceiling(npy/2.)
             do i=1,ceiling(npx/2.)
                ang  = get_angle(2, grid_global(i,j+1,1:2,1), grid_global(i,j,1:2,1), grid_global(i+1,j,1:2,1))
                ang  = ABS(90.0 - ang)
                angs(i,j) = ang

                if ( (i==1) .and. (j==1) ) then
                else 
                   angAV = angAV + ang
                   angM  = MAX(angM,ang)
                   angN  = MIN(angN,ang)
                endif

                dx_local = dx_global(i,j,1)
                dy_local = dy_global(i,j,1)

                dxAV  = dxAV + 0.5 * (dx_local + dy_local)
                dxM   = MAX(dxM,dx_local)
                dxM   = MAX(dxM,dy_local)
                dxN   = MIN(dxN,dx_local)
                dxN   = MIN(dxN,dy_local)
                dxs(i,j) = dy_local !0.5 * (dx_local + dy_local)

                asp   = ABS(dx_local/dy_local)
                if (asp < 1.0) asp = 1.0/asp
                asps(i,j) = asp 
                aspAV = aspAV + asp
                aspM  = MAX(aspM,asp)
                aspN  = MIN(aspN,asp)
             enddo
          enddo
          angAV = angAV / ( (ceiling(npy/2.0))*(ceiling(npx/2.0)) - 1 )
          dxAV  = dxAV  / ( (ceiling(npy/2.0))*(ceiling(npx/2.0)) )
          aspAV = aspAV / ( (ceiling(npy/2.0))*(ceiling(npx/2.0)) )
          write(*,*  ) ''
          write(*,*  ) ' Cubed-Sphere Grid Stats : ', npx,'x',npy,'x',nregions
          write(*,201) '      Grid Length               : min: ', dxN,' max: ', dxM,' avg: ', dxAV, ' min/max: ',dxN/dxM
          write(*,200) '      Deviation from Orthogonal : min: ',angN,' max: ',angM,' avg: ',angAV
          write(*,200) '      Aspect Ratio              : min: ',aspN,' max: ',aspM,' avg: ',aspAV
          write(*,*  ) ''
          write(gcharFile,202) TRIM(grid_name),'_chars_',npx,'x',npy,'.dat'
          fileLun=get_unit()
          open(unit=fileLun,file=gcharFile, form='unformatted', access='direct',  &
               recl=((npx/2)+1)*((npy/2)+1)*8, status='unknown')
          write(fileLun,rec=1) angs
          write(fileLun,rec=2) asps
          write(fileLun,rec=3)  dxs
          do j=1,(npy/2.0)+1
             do i=1,(npx/2.0)+1
                do n=1,ndims
                   p_lL(n) = grid_global(i  ,j  ,n,1)
                   p_uL(n) = grid_global(i  ,j+1,n,1)
                   p_lR(n) = grid_global(i+1,j  ,n,1)
                   p_uR(n) = grid_global(i+1,j+1,n,1)
                enddo
                   ! Spherical Excess Formula
                   angs(i,j) = get_area(p_lL, p_uL, p_lR, p_uR, radius)
             enddo
          enddo
          write(fileLun,rec=4) angs
          close(unit=fileLun)

       endif
    endif

  contains

    subroutine setup_cartesian(npx, npy)
       integer, intent(in):: npx, npy
       real lat_rad, lon_rad, domain_rad
       integer i,j

       domain_rad = pi/16.   ! arbitrary
       lat_rad = deglat * pi/180.
       lon_rad = 0.          ! arbitrary

       dx(:,:)  = dx_const
       rdx(:,:) = 1./dx_const
       dy(:,:)  = dy_const
       rdy(:,:) = 1./dy_const
       
       dxc(:,:)  = dx_const
       rdxc(:,:) = 1./dx_const
       dyc(:,:)  = dy_const
       rdyc(:,:) = 1./dy_const
       
       dxa(:,:)  = dx_const
       rdxa(:,:) = 1./dx_const
       dya(:,:)  = dy_const
       rdya(:,:) = 1./dy_const
       
       area(:,:)  = dx_const*dy_const
       rarea(:,:) = 1./(dx_const*dy_const)
       
       area_c(:,:)  = dx_const*dy_const
       rarea_c(:,:) = 1./(dx_const*dy_const)
       
! The following is a hack to get pass the am2 phys init:
       do j=max(1,jsd),min(jed,npy)
          do i=max(1,isd),min(ied,npx)
             grid(i,j,1) = lon_rad - 0.5*domain_rad + real(i-1)/real(npx-1)*domain_rad
             grid(i,j,2) = lat_rad - 0.5*domain_rad + real(j-1)/real(npy-1)*domain_rad
          enddo
       enddo

       agrid(:,:,1)  = lon_rad
       agrid(:,:,2)  = lat_rad
       
       sina(:,:) = 1.
       cosa(:,:) = 0.

       e1(1,:,:) = 1.
       e1(2,:,:) = 0.
       e1(3,:,:) = 0.

       e2(1,:,:) = 0.
       e2(2,:,:) = 1.
       e2(3,:,:) = 0.

    end subroutine setup_cartesian


    subroutine setup_latlon()
      real, parameter :: big_number = 1.e30
      real :: lon_start, lat_start, area_j

      dl = (deglon_stop-deglon_start)*pi/(180.*(npx-1))
      dp = (deglat_stop-deglat_start)*pi/(180.*(npy-1))

      lon_start = deglon_start*pi/180.
      lat_start = deglat_start*pi/180.
      
      do j=jsd,jed+1
         do i=isd,ied+1
            grid(i,j,1) = lon_start + real(i-1)*dl
            grid(i,j,2) = lat_start + real(j-1)*dp
         enddo
      enddo

      do j=jsd,jed
         do i=isd,ied
            agrid(i,j,1) = (grid(i,j,1) + grid(i+1,j,1))/2.
            agrid(i,j,2) = (grid(i,j,2) + grid(i,j+1,2))/2.
         enddo
      enddo


      do j=jsd,jed
         do i=isd,ied+1
            dxc(i,j) = dl*radius*cos(agrid(is,j,2))
            rdxc(i,j) = 1./dxc(i,j)
         enddo
      enddo
      do j=jsd,jed+1
         do i=isd,ied
            dyc(i,j) = dp*radius
            rdyc(i,j) = 1./dyc(i,j)
         enddo
      enddo

      do j=jsd,jed
         do i=isd,ied
            dxa(i,j) = dl*radius*cos(agrid(i,j,2))
            dya(i,j) = dp*radius
            rdxa(i,j) = 1./dxa(i,j)
            rdya(i,j) = 1./dya(i,j)
         enddo
      enddo
          
      do j=jsd,jed+1
         do i=isd,ied
            dx(i,j) = dl*radius*cos(grid(i,j,2))
            rdx(i,j) = 1./dx(i,j)
         enddo
      enddo
      do j=jsd,jed
         do i=isd,ied+1
            dy(i,j) = dp*radius
            rdy(i,j) = 1./dy(i,j)
         enddo
      enddo

      do j=jsd,jed
         area_j = radius*radius*dl*(sin(grid(is,j+1,2))-sin(grid(is,j,2)))
         do i=isd,ied
            area(i,j) = area_j
            rarea(i,j) = 1./area_j
         enddo
      enddo

      do j=jsd+1,jed
         area_j = radius*radius*dl*(sin(agrid(is,j,2))-sin(agrid(is,j-1,2)))
         do i=isd,ied+1
            area_c(i,j) = area_j
            rarea_c(i,j) = 1./area_j
         enddo
      enddo
      if (jsd==1) then
         j=1
         area_j = radius*radius*dl*(sin(agrid(is,j,2))-sin(agrid(is,j,2)-dp))
         do i=isd,ied+1
            area_c(i,j) = area_j
            rarea_c(i,j) = 1./area_j
         enddo
      endif
      if (jed+1==npy) then
         j=npy
         area_j = radius*radius*dl*(sin(agrid(is,j-1,2)+dp)-sin(agrid(is,j-1,2)))
         do i=isd,ied+1
            area_c(i,j) = area_j
            rarea_c(i,j) = 1./area_j
         enddo
      endif
      call mpp_update_domains( area_c, domain, position=CORNER, complete=.true.)

      sina(:,:) = 1.
      cosa(:,:) = 0.
      
      e1(1,:,:) = 1.
      e1(2,:,:) = 0.
      e1(3,:,:) = 0.
      
      e2(1,:,:) = 0.
      e2(2,:,:) = 1.
      e2(3,:,:) = 0.

    end subroutine setup_latlon
  
   end subroutine init_grid


      subroutine cartesian_to_spherical(x, y, z, lon, lat, r) 
      real , intent(IN)  :: x, y, z
      real , intent(OUT) :: lon, lat, r

      r = SQRT(x*x + y*y + z*z)
      if ( (abs(x) + abs(y)) < 1.E-10 ) then       ! poles:
           lon = 0.
      else
           lon = ATAN2(y,x)    ! range: [-pi,pi]
      endif 

#ifdef RIGHT_HAND
      lat = asin(z/r)
#else
      lat = ACOS(z/r) - pi/2.
#endif

      end subroutine cartesian_to_spherical
 

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     spherical_to_cartesian :: convert from spheircal coordinates to xyz coords
! 
      subroutine spherical_to_cartesian(lon, lat, r, x, y, z)

         real , intent(IN)  :: lon, lat, r
         real , intent(OUT) :: x, y, z

         x = r * COS(lon) * cos(lat)
         y = r * SIN(lon) * cos(lat)

#ifdef RIGHT_HAND
         z =  r * SIN(lat)
#else
         z = -r * sin(lat)
#endif

      end subroutine spherical_to_cartesian



!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     rot_3d :: rotate points on a sphere in xyz coords (convert angle from
!               degrees to radians if necessary)
!
      subroutine rot_3d(axis, x1in, y1in, z1in, angle, x2out, y2out, z2out, degrees, convert)

         integer, intent(IN) :: axis         ! axis of rotation 1=x, 2=y, 3=z
         real , intent(IN)    :: x1in, y1in, z1in
         real , intent(INOUT) :: angle        ! angle to rotate in radians
         real , intent(OUT)   :: x2out, y2out, z2out
         integer, intent(IN), optional :: degrees ! if present convert angle 
                                                  ! from degrees to radians
         integer, intent(IN), optional :: convert ! if present convert input point
                                                  ! from spherical to cartesian, rotate, 
                                                  ! and convert back

         real  :: c, s
         real  :: x1,y1,z1, x2,y2,z2

         if ( present(convert) ) then
           call spherical_to_cartesian(x1in, y1in, z1in, x1, y1, z1)
         else
           x1=x1in
           y1=y1in
           z1=z1in
         endif

         if ( present(degrees) ) then
            angle = angle*torad
         endif

         c = COS(angle)
         s = SIN(angle)

         SELECT CASE(axis)
             
            CASE(1)
               x2 =  x1
               y2 =  c*y1 + s*z1
               z2 = -s*y1 + c*z1
            CASE(2)
               x2 = c*x1 - s*z1
               y2 = y1
               z2 = s*x1 + c*z1
            CASE(3)
               x2 =  c*x1 + s*y1
               y2 = -s*x1 + c*y1
               z2 = z1
            CASE DEFAULT
              write(*,*) "Invalid axis: must be 1 for X, 2 for Y, 3 for Z."
 
         END SELECT

         if ( present(convert) ) then
           call cartesian_to_spherical(x2, y2, z2, x2out, y2out, z2out)
         else
           x2out=x2
           y2out=y2
           z2out=z2
         endif

      end subroutine rot_3d





      real  function get_area_tri(ndims, p_1, p_2, p_3) &
                        result (myarea)
 
!     get_area_tri :: get the surface area of a cell defined as a triangle
!                  on the sphere. Area is computed as the spherical excess
!                  [area units are based on the units of radius]
 

      integer, intent(IN)    :: ndims          ! 2=lat/lon, 3=xyz
      real , intent(IN)    :: p_1(ndims) ! 
      real , intent(IN)    :: p_2(ndims) ! 
      real , intent(IN)    :: p_3(ndims) ! 

      real  :: angA, angB, angC

        if ( ndims==3 ) then
            angA = spherical_angle(p_1, p_2, p_3)
            angB = spherical_angle(p_2, p_3, p_1)
            angC = spherical_angle(p_3, p_1, p_2)
        else
            angA = get_angle(ndims, p_1, p_2, p_3, 1)
            angB = get_angle(ndims, p_2, p_3, p_1, 1)
            angC = get_angle(ndims, p_3, p_1, p_2, 1)
        endif

        myarea = (angA+angB+angC - pi) * radius**2

      end function get_area_tri
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     grid_area :: get surface area on grid in lat/lon coords or xyz coords
!                    (determined by ndims argument 2=lat/lon, 3=xyz)
!                    [area is returned in m^2 on Unit sphere]
!
      subroutine grid_area(nx, ny, ndims, nregions )

         integer, intent(IN) :: nx, ny, ndims, nregions

         real  :: p_lL(ndims) ! lower Left
         real  :: p_uL(ndims) ! upper Left
         real  :: p_lR(ndims) ! lower Right
         real  :: p_uR(ndims) ! upper Right
         real  :: a1, d1, d2, mydx, mydy

         real  :: p1(ndims), p2(ndims), p3(ndims), pi1(ndims), pi2(ndims)

         real  :: maxarea, minarea

         integer :: i,j,n, nreg

         real, allocatable :: p_R8(:,:,:) 

         maxarea = -1.e25
         minarea =  1.e25

         globalarea = 0.0
         do j=js,je
            do i=is,ie
               do n=1,ndims
               if ( stretched_grid ) then
                  p_lL(n) = grid(i  ,j  ,n)
                  p_uL(n) = grid(i  ,j+1,n)
                  p_lR(n) = grid(i+1,j  ,n)
                  p_uR(n) = grid(i+1,j+1,n)
               else
                  p_lL(n) = grid(iinta(1,i,j), jinta(1,i,j), n)
                  p_uL(n) = grid(iinta(2,i,j), jinta(2,i,j), n)
                  p_lR(n) = grid(iinta(4,i,j), jinta(4,i,j), n)
                  p_uR(n) = grid(iinta(3,i,j), jinta(3,i,j), n)
               endif
               enddo

           ! Spherical Excess Formula
              area(i,j) = get_area(p_lL, p_uL, p_lR, p_uR, radius)
              maxarea=MAX(area(i,j),maxarea)
              minarea=MIN(area(i,j),minarea)
              globalarea = globalarea + area(i,j)
            enddo
         enddo

         allocate( p_R8(nx-1,ny-1,ntiles_g) )   ! this is a "global" array
         do j=js,je
            do i=is,ie
               p_R8(i,j,tile) = area(i,j)
            enddo
         enddo
         call mp_gather(p_R8, is,ie, js,je, nx-1, ny-1, ntiles_g)
         if (gid == masterproc) then
            globalarea = 0.0
            do n=1,ntiles_g
               do j=1,ny-1
                  do i=1,nx-1
                     globalarea = globalarea + p_R8(i,j,n)
                  enddo
               enddo
            enddo
         endif

         call mpp_broadcast(globalarea, masterproc)

         deallocate( p_R8 )

         call mp_reduce_max(maxarea)
         minarea = -minarea                  
         call mp_reduce_max(minarea)
         minarea = -minarea

        if (gid == masterproc) write(*,209) 'MAX    AREA (m*m):', maxarea,            '          MIN AREA (m*m):', minarea
        if (gid == masterproc) write(*,209) 'GLOBAL AREA (m*m):', globalarea, ' IDEAL GLOBAL AREA (m*m):', 4.0*pi*radius**2
 209  format(A,e21.14,A,e21.14)

         do j=js,je+1
            do i=is,ie+1
               do n=1,ndims
               if ( stretched_grid ) then
                  p_lL(n) = agrid(i-1,j-1,n)
                  p_lR(n) = agrid(i  ,j-1,n)
                  p_uL(n) = agrid(i-1,j  ,n)
                  p_uR(n) = agrid(i  ,j  ,n)
               else
                  p_lL(n) = agrid(iintb(1,i,j), jintb(1,i,j), n)
                  p_lR(n) = agrid(iintb(2,i,j), jintb(2,i,j), n)
                  p_uL(n) = agrid(iintb(4,i,j), jintb(4,i,j), n)
                  p_uR(n) = agrid(iintb(3,i,j), jintb(3,i,j), n)
               endif
               enddo
              ! Spherical Excess Formula
                area_c(i,j) = get_area(p_lL, p_uL, p_lR, p_uR, radius)
            enddo
         enddo

! Corners: assuming triangular cells
         if (cubed_sphere) then
! SW:
            i=1
            j=1
            if ( (is==1) .and. (js==1) ) then
              do n=1,ndims
               if ( stretched_grid ) then
                    p1(n) = agrid(i-1,j  ,n)
                    p2(n) = agrid(i  ,j  ,n)
                    p3(n) = agrid(i  ,j-1,n)
               else
                    p1(n) = agrid(iintb(1,i,j), jintb(1,i,j), n)
                    p2(n) = agrid(iintb(2,i,j), jintb(2,i,j), n)
                    p3(n) = agrid(iintb(3,i,j), jintb(3,i,j), n)
               endif
              enddo
              area_c(i,j) = get_area_tri(ndims, p1, p2, p3)
            endif

            i=nx
            j=1
            if ( (ie+1==nx) .and. (js==1) ) then
              do n=1,ndims
               if ( stretched_grid ) then
               p1(n) = agrid(i  ,j  ,n)
               p2(n) = agrid(i-1,j  ,n)
               p3(n) = agrid(i-1,j-1,n)
               else
               p1(n) = agrid(iintb(1,i,j), jintb(1,i,j), n)
               p2(n) = agrid(iintb(2,i,j), jintb(2,i,j), n)
               p3(n) = agrid(iintb(3,i,j), jintb(3,i,j), n)
               endif
              enddo
              area_c(i,j) = get_area_tri(ndims, p1, p2, p3)
            endif

            i=nx
            j=ny
            if ( (ie+1==nx) .and. (je+1==ny) ) then
              do n=1,ndims
               if ( stretched_grid ) then
               p1(n) = agrid(i-1,j  ,n)
               p2(n) = agrid(i-1,j-1,n)
               p3(n) = agrid(i  ,j-1,n)
               else
               p1(n) = agrid(iintb(1,i,j), jintb(1,i,j), n)
               p2(n) = agrid(iintb(2,i,j), jintb(2,i,j), n)
               p3(n) = agrid(iintb(3,i,j), jintb(3,i,j), n)
               endif
              enddo
              area_c(i,j) = get_area_tri(ndims, p1, p2, p3)
            endif

            i=1
            j=ny
            if ( (is==1) .and. (je+1==ny) ) then
              do n=1,ndims
               if ( stretched_grid ) then
               p1(n) = agrid(i  ,j  ,n)
               p2(n) = agrid(i  ,j-1,n)
               p3(n) = agrid(i-1,j-1,n)
               else
               p1(n) = agrid(iintb(1,i,j), jintb(1,i,j), n)
               p2(n) = agrid(iintb(2,i,j), jintb(2,i,j), n)
               p3(n) = agrid(iintb(3,i,j), jintb(3,i,j), n)
              endif
              enddo
              area_c(i,j) = get_area_tri(ndims, p1, p2, p3)
            endif
         endif

      end subroutine grid_area



      real  function get_angle(ndims, p1, p2, p3, rad) result (angle)
!     get_angle :: get angle between 3 points on a sphere in lat/lon coords or
!                  xyz coords (determined by ndims argument 2=lat/lon, 3=xyz)
!                  [angle is returned in degrees]

         integer, intent(IN) :: ndims         ! 2=lat/lon, 3=xyz
         real , intent(IN)   :: p1(ndims)
         real , intent(IN)   :: p2(ndims)
         real , intent(IN)   :: p3(ndims)
         integer, intent(in), optional:: rad

         real  :: e1(3), e2(3), e3(3)

         if (ndims == 2) then
            call spherical_to_cartesian(p2(1), p2(2), 1., e1(1), e1(2), e1(3))
            call spherical_to_cartesian(p1(1), p1(2), 1., e2(1), e2(2), e2(3))
            call spherical_to_cartesian(p3(1), p3(2), 1., e3(1), e3(2), e3(3))
         else
            e1 = p2; e2 = p1; e3 = p3
         endif

! High precision version:
         if ( present(rad) ) then
           angle = spherical_angle(e1, e2, e3)
         else
           angle = todeg * spherical_angle(e1, e2, e3)
         endif

      end function get_angle
 

 
      subroutine mp_update_dwinds_2d(u, v, npx, npy)
        use mpp_parameter_mod, only: DGRID_NE
         real  , intent(INOUT)   :: u(isd:ied  ,jsd:jed+1) ! D-grid u-wind field
         real  , intent(INOUT)   :: v(isd:ied+1,jsd:jed  ) ! D-grid v-wind field
         integer,      intent(IN) :: npx, npy

         call mpp_update_domains( u, v, domain, gridtype=DGRID_NE, complete=.true.)
!        call fill_corners(u , v , npx, npy, VECTOR=.true., DGRID=.true.)

      end subroutine mp_update_dwinds_2d
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
      subroutine mp_update_dwinds_3d(u, v, npx, npy, npz)
        use mpp_parameter_mod, only: DGRID_NE
         real  , intent(INOUT)   :: u(isd:ied  ,jsd:jed+1,npz) ! D-grid u-wind field
         real  , intent(INOUT)   :: v(isd:ied+1,jsd:jed  ,npz) ! D-grid v-wind field
         integer,      intent(IN) :: npx, npy, npz
         integer k

      call mpp_update_domains( u, v, domain, gridtype=DGRID_NE, complete=.true.)
!     do k=1,npz
!        call fill_corners(u(isd:,jsd:,k) , v(isd:,jsd:,k) , npx, npy, VECTOR=.true., DGRID=.true.)
!     enddo

      end subroutine mp_update_dwinds_3d



      subroutine atob_s(qin, qout, npx, npy, altInterp)

!     atob_s :: interpolate scalar from the A-Grid to the B-grid
!
         integer,      intent(IN) :: npx, npy
         real  , intent(IN)    ::  qin(isd:ied  ,jsd:jed  )    ! A-grid field
         real  , intent(OUT)   :: qout(isd:ied+1,jsd:jed+1)    ! Output  B-grid field
         integer, OPTIONAL, intent(IN) :: altInterp 

         integer :: i,j,n

         real :: tmp1j(jsd:jed+1)
         real :: tmp2j(jsd:jed+1)
         real :: tmp3j(jsd:jed+1)
         real :: tmp1i(isd:ied+1)
         real :: tmp2i(isd:ied+1)
         real :: tmp3i(isd:ied+1)
         real :: tmpq(isd:ied  ,jsd:jed  )
         real :: tmpq1(isd:ied+1,jsd:jed+1)
         real :: tmpq2(isd:ied+1,jsd:jed+1)

         if (present(altInterp)) then

         tmpq(:,:) = qin(:,:)

         call fill_corners(tmpq  , npx, npy, FILL=XDir, AGRID=.true.)
! ATOC
         do j=jsd,jed
            call interp_left_edge_1d(tmpq1(:,j), tmpq(:,j), dxa(:,j), isd, ied, altInterp) 
         enddo

         call fill_corners(tmpq  , npx, npy, FILL=YDir, AGRID=.true.)
! ATOD
         do i=isd,ied
            tmp1j(jsd:jed) = 0.0 
            tmp2j(jsd:jed) = tmpq(i,jsd:jed)
            tmp3j(jsd:jed) = dya(i,jsd:jed)
            call interp_left_edge_1d(tmp1j, tmp2j, tmp3j, jsd, jed, altInterp)
            tmpq2(i,jsd:jed) = tmp1j(jsd:jed)
         enddo

! CTOB
         do i=isd,ied
            tmp1j(:) = tmpq1(i,:)
            tmp2j(:) = tmpq1(i,:)
            tmp3j(:) = 1.0  ! Uniform Weighting missing first value so will not reproduce
            call interp_left_edge_1d(tmp1j, tmp2j, tmp3j, jsd, jed+1, altInterp) 
            tmpq1(i,:) = tmp1j(:)
         enddo

! DTOB
         do j=jsd,jed
            tmp1i(:) = tmpq2(:,j)
            tmp2i(:) = tmpq2(:,j)
            tmp3i(:) = 1.0  ! Uniform Weighting missing first value so will not reproduce
            call interp_left_edge_1d(tmp1i, tmp2i, tmp3i, isd, ied+1, altInterp)
            tmpq2(:,j) = tmp1i(:)
         enddo

! Average 
         do j=jsd,jed+1
            do i=isd,ied+1
               qout(i,j) = 0.5 * (tmpq1(i,j) + tmpq2(i,j))
            enddo
         enddo

! Fix Corners
         if (cubed_sphere) then
            i=1
            j=1
            if ( (is==i) .and. (js==j) ) then
               qout(i,j) = (1./3.) * (qin(i,j) + qin(i-1,j) + qin(i,j-1))
            endif

            i=npx
            j=1
            if ( (ie+1==i) .and. (js==j) ) then
               qout(i,j) = (1./3.) * (qin(i-1,j) + qin(i-1,j-1) + qin(i,j))
            endif

            i=1
            j=npy
            if ( (is==i) .and. (je+1==j) ) then
               qout(i,j) = (1./3.) * (qin(i,j-1) + qin(i-1,j-1) + qin(i,j))
            endif

            i=npx
            j=npy
            if ( (ie+1==i) .and. (je+1==j) ) then
               qout(i,j) = (1./3.) * (qin(i-1,j-1) + qin(i,j-1) + qin(i-1,j))
            endif
        endif

        else ! altInterp

            do j=js,je+1
               do i=is,ie+1
                  qout(i,j) = 0.25 * (qin(i-1,j) + qin(i-1,j-1) + &
                                      qin(i  ,j) + qin(i  ,j-1))
               enddo
            enddo
            i=1
            j=1
            if ( (is==i) .and. (js==j) ) then
               qout(i,j) = (1./3.) * (qin(i,j) + qin(i-1,j) + qin(i,j-1))
            endif

            i=npx
            j=1
            if ( (ie+1==i) .and. (js==j) ) then
               qout(i,j) = (1./3.) * (qin(i-1,j) + qin(i-1,j-1) + qin(i,j))
            endif

            i=1
            j=npy
            if ( (is==i) .and. (je+1==j) ) then
               qout(i,j) = (1./3.) * (qin(i,j-1) + qin(i-1,j-1) + qin(i,j))
            endif

            i=npx
            j=npy
            if ( (ie+1==i) .and. (je+1==j) ) then
               qout(i,j) = (1./3.) * (qin(i-1,j-1) + qin(i,j-1) + qin(i-1,j))
            endif

        endif ! altInterp

      end subroutine atob_s
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     atod :: interpolate from the A-Grid to the D-grid
!
      subroutine atod(uin, vin, uout, vout, npx, npy, ng)


         integer,      intent(IN) :: npx, npy, ng
         real  , intent(IN)    ::  uin(isd:ied  ,jsd:jed  ) ! A-grid u-wind field
         real  , intent(IN)    ::  vin(isd:ied  ,jsd:jed  ) ! A-grid v-wind field
         real  , intent(OUT)   :: uout(isd:ied  ,jsd:jed+1) ! D-grid u-wind field
         real  , intent(OUT)   :: vout(isd:ied+1,jsd:jed  ) ! D-grid v-wind field

         integer :: i,j
         real :: tmp1i(isd:ied+1)
         real :: tmp2i(isd:ied)
         real :: tmp3i(isd:ied)
         real :: tmp1j(jsd:jed+1)
         real :: tmp2j(jsd:jed)
         real :: tmp3j(jsd:jed)

         do j=jsd+1,jed
            tmp1i(:) = 0.0
            tmp2i(:) = vin(:,j)*dxa(:,j)
            tmp3i(:) = dxa(:,j)
            call interp_left_edge_1d(tmp1i, tmp2i, tmp3i, isd, ied, interpOrder)
            vout(:,j) = tmp1i(:)/dxc(:,j)
         enddo
         do i=isd+1,ied
            tmp1j(:) = 0.0
            tmp2j(:) = uin(i,:)*dya(i,:)
            tmp3j(:) = dya(i,:)
            call interp_left_edge_1d(tmp1j, tmp2j, tmp3j, jsd, jed, interpOrder)
            uout(i,:) = tmp1j(:)/dyc(i,:)
         enddo
         call mp_update_dwinds(uout, vout, npx, npy)
         call fill_corners(uout, vout, npx, npy, VECTOR=.true., DGRID=.true.)
      end subroutine atod
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     dtoa :: interpolate from the D-Grid to the A-grid
!
      subroutine dtoa(uin, vin, uout, vout, npx, npy, ng)

         integer,      intent(IN) :: npx, npy, ng
         real  , intent(IN)    ::  uin(isd:ied  ,jsd:jed+1)    ! D-grid u-wind field
         real  , intent(IN)    ::  vin(isd:ied+1,jsd:jed  )    ! D-grid v-wind field
         real  , intent(OUT)   :: uout(isd:ied  ,jsd:jed  )    ! A-grid u-wind field
         real  , intent(OUT)   :: vout(isd:ied  ,jsd:jed  )    ! A-grid v-wind field

         integer :: i,j,n

         real :: tmp1i(isd:ied+1)
         real :: tmp2i(isd:ied+1)
         real :: tmp3i(isd:ied+1)
         real :: tmp1j(jsd:jed+1)
         real :: tmp2j(jsd:jed+1)
         real :: tmp3j(jsd:jed+1)

#ifdef VORT_ON
! circulation (therefore, vort) conserving:
         do j=jsd,jed
            do i=isd,ied
                uout(i,j) = 0.5*(uin(i,j)*dx(i,j)+uin(i,j+1)*dx(i,j+1))/dxa(i,j)
                vout(i,j) = 0.5*(vin(i,j)*dy(i,j)+vin(i+1,j)*dy(i+1,j))/dya(i,j)
            enddo
         enddo
#else
         do i=isd,ied
            tmp1j(:) = 0.0
            tmp2j(:) = uin(i,:)*dyc(i,:)
            tmp3j(:) = dyc(i,:)
            call interp_left_edge_1d(tmp1j, tmp2j, tmp3j, jsd, jed+1, interpOrder) 
            uout(i,jsd:jed) = tmp1j(jsd+1:jed+1)/dya(i,jsd:jed)
         enddo
         do j=jsd,jed
            tmp1i(:) = 0.0
            tmp2i(:) = vin(:,j)*dxc(:,j)
            tmp3i(:) = dxc(:,j)
            call interp_left_edge_1d(tmp1i, tmp2i, tmp3i, isd, ied+1, interpOrder) 
            vout(isd:ied,j) = tmp1i(isd+1:ied+1)/dxa(isd:ied,j)
         enddo
#endif

      end subroutine dtoa
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     atoc :: interpolate from the A-Grid to the C-grid
!
      subroutine atoc(uin, vin, uout, vout, npx, npy, ng, noComm)


         integer,      intent(IN) :: npx, npy, ng
         real  , intent(IN)    ::  uin(isd:ied  ,jsd:jed  ) ! A-grid u-wind field
         real  , intent(IN)    ::  vin(isd:ied  ,jsd:jed  ) ! A-grid v-wind field
         real  , intent(OUT)   :: uout(isd:ied+1,jsd:jed  ) ! C-grid u-wind field
         real  , intent(OUT)   :: vout(isd:ied  ,jsd:jed+1) ! C-grid v-wind field
         logical, OPTIONAL, intent(IN)   :: noComm

         real :: ang1
         integer :: i,j,n

         real :: tmp1i(isd:ied+1)
         real :: tmp2i(isd:ied)
         real :: tmp3i(isd:ied)
         real :: tmp1j(jsd:jed+1)
         real :: tmp2j(jsd:jed)
         real :: tmp3j(jsd:jed)

#if !defined(ALT_INTERP)
#ifdef VORT_ON
! Circulation conserving
         do j=jsd,jed
            do i=isd+1,ied
               uout(i,j) = ( uin(i,j)*dxa(i,j) + uin(i-1,j)*dxa(i-1,j) )    &
                           /        ( dxa(i,j) +            dxa(i-1,j) )
            enddo
         enddo
         do j=jsd+1,jed
            do i=isd,ied
               vout(i,j) = ( vin(i,j)*dya(i,j) + vin(i,j-1)*dya(i,j-1) )    &
                           /        ( dya(i,j) +            dya(i,j-1) )
            enddo
         enddo
#else
         do j=jsd,jed
            call interp_left_edge_1d(uout(:,j), uin(:,j), dxa(:,j), isd, ied, interpOrder)
         enddo
         do i=isd,ied
!!$            tmp1j(:) = vout(i,:)
            tmp2j(:) = vin(i,:)
            tmp3j(:) = dya(i,:)
            call interp_left_edge_1d(tmp1j, tmp2j, tmp3j, jsd, jed, interpOrder)
            vout(i,:) = tmp1j(:)
         enddo 
#endif
#else

         do j=jsd,jed
!!$            tmp1i(:) = uout(:,j)
            tmp2i(:) = uin(:,j)*dya(:,j)
            tmp3i(:) = dxa(:,j)
            call interp_left_edge_1d(tmp1i, tmp2i, tmp3i, isd, ied, interpOrder)
            uout(:,j) = tmp1i(:)/dy(:,j)
         enddo
         do i=isd,ied
!!$            tmp1j(:) = vout(i,:)
            tmp2j(:) = vin(i,:)*dxa(i,:)
            tmp3j(:) = dya(i,:)
            call interp_left_edge_1d(tmp1j, tmp2j, tmp3j, jsd, jed, interpOrder)
            vout(i,:) = tmp1j(:)/dx(i,:)
         enddo

       if (cubed_sphere) then
         csFac = COS(30.0*PI/180.0)
      ! apply Corner scale factor for interp on Cubed-Sphere
         if ( (is==1) .and. (js==1) ) then
            i=1
            j=1
            uout(i,j)=uout(i,j)*csFac
            uout(i,j-1)=uout(i,j-1)*csFac
            vout(i,j)=vout(i,j)*csFac
            vout(i-1,j)=vout(i-1,j)*csFac
         endif
         if ( (is==1) .and. (je==npy-1) ) then
            i=1
            j=npy-1
            uout(i,j)=uout(i,j)*csFac
            uout(i,j+1)=uout(i,j+1)*csFac
            vout(i,j+1)=vout(i,j+1)*csFac
            vout(i-1,j+1)=vout(i-1,j+1)*csFac
         endif
         if ( (ie==npx-1) .and. (je==npy-1) ) then
            i=npx-1
            j=npy-1
            uout(i+1,j)=uout(i+1,j)*csFac
            uout(i+1,j+1)=uout(i+1,j+1)*csFac
            vout(i,j+1)=vout(i,j+1)*csFac
            vout(i+1,j+1)=vout(i+1,j+1)*csFac
         endif
         if ( (ie==npx-1) .and. (js==1) ) then
            i=npx-1
            j=1
            uout(i+1,j)=uout(i+1,j)*csFac
            uout(i+1,j-1)=uout(i+1,j-1)*csFac
            vout(i,j)=vout(i,j)*csFac
            vout(i+1,j)=vout(i+1,j)*csFac
         endif
       endif

#endif

         if (present(noComm)) then
            if (.not. noComm) call mpp_update_domains( uout,vout, domain, gridtype=CGRID_NE_PARAM, complete=.true.)
         else
            call mpp_update_domains( uout,vout, domain, gridtype=CGRID_NE_PARAM, complete=.true.)
         endif
         call fill_corners(uout, vout, npx, npy, VECTOR=.true., CGRID=.true.)

      end subroutine atoc
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     ctoa :: interpolate from the C-Grid to the A-grid
!
      subroutine ctoa(uin, vin, uout, vout, npx, npy, ng)


         integer,      intent(IN) :: npx, npy, ng 
         real  , intent(IN)    ::  uin(isd:ied+1,jsd:jed  )    ! C-grid u-wind field
         real  , intent(IN)    ::  vin(isd:ied  ,jsd:jed+1)    ! C-grid v-wind field
         real  , intent(OUT)   :: uout(isd:ied  ,jsd:jed  )    ! A-grid u-wind field
         real  , intent(OUT)   :: vout(isd:ied  ,jsd:jed  )    ! A-grid v-wind field

         integer :: i,j

         real :: tmp1i(isd:ied+1)
         real :: tmp2i(isd:ied+1)
         real :: tmp3i(isd:ied+1)
         real :: tmp1j(jsd:jed+1)
         real :: tmp2j(jsd:jed+1)
         real :: tmp3j(jsd:jed+1)

        ! do j=jsd,jed
        !    do i=isd,ied
        !       uout(i,j) = 0.5 * (uin(i,j)*dy(i,j) + uin(i+1,j)*dy(i+1,j))/dya(i,j)
        !    enddo
        !  enddo
        ! do j=jsd,jed
        !    do i=isd,ied
        !       vout(i,j) = 0.5 * (vin(i,j)*dx(i,j) + vin(i,j+1)*dx(i,j+1))/dxa(i,j)
        !    enddo
        ! enddo
         do i=isd,ied
            tmp1j(:) = 0.0
            tmp2j(:) = vin(i,:)*dx(i,:)
            tmp3j(:) = dyc(i,:)
            call interp_left_edge_1d(tmp1j, tmp2j, tmp3j, jsd, jed+1, interpOrder)
            vout(i,jsd:jed) = tmp1j(jsd+1:jed+1)/dxa(i,jsd:jed)
         enddo
         do j=jsd,jed
            tmp1i(:) = 0.0
            tmp2i(:) = uin(:,j)*dy(:,j)
            tmp3i(:) = dxc(:,j)
            call interp_left_edge_1d(tmp1i, tmp2i, tmp3i, isd, ied+1, interpOrder)
            uout(isd:ied,j) = tmp1i(isd+1:ied+1)/dya(isd:ied,j)
         enddo

      end subroutine ctoa
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     rotate_winds :: rotate winds from the sphere-to-cube || cube-to-sphere
!
      subroutine rotate_winds(myU, myV, p1, p2, p3, p4, t1, ndims, dir)


         integer,      intent(IN) :: ndims
         real  , intent(INOUT) :: myU    ! u-wind field
         real  , intent(INOUT) :: myV    ! v-wind field
         real  , intent(IN)    :: p1(ndims)    !             p4     
         real  , intent(IN)    :: p2(ndims)    !                    
         real  , intent(IN)    :: p3(ndims)    !        p1   t1   p3
         real  , intent(IN)    :: p4(ndims)    !                    
         real  , intent(IN)    :: t1(ndims)    !             p2     
         integer,   intent(IN)    :: dir   ! Direction ; 1=>sphere-to-cube  2=> cube-to-sphere

         real :: ee1(3), ee2(3), ee3(3), elon(3), elat(3)

         real :: g11, g12, g21, g22

         real :: newu, newv

         call get_unit_vector(p3, t1, p1, ee1)
         call get_unit_vector(p4, t1, p2, ee2)
         elon(1) = -SIN(t1(1) - pi)
         elon(2) =  COS(t1(1) - pi)
         elon(3) = 0.0
         elat(1) = -SIN(t1(2))*COS(t1(1) - pi)
         elat(2) = -SIN(t1(2))*SIN(t1(1) - pi)
         elat(3) =  COS(t1(2))

         g11 = inner_prod(ee1,elon)
         g12 = inner_prod(ee1,elat)
         g21 = inner_prod(ee2,elon)
         g22 = inner_prod(ee2,elat)

         if (dir == 1) then    ! Sphere to Cube Rotation
            newu = myU*g11 + myV*g12
            newv = myU*g21 + myV*g22
         else
            newu = ( myU*g22 - myV*g12)/(g11*g22 - g21*g12) 
            newv = (-myU*g21 + myV*g11)/(g11*g22 - g21*g12)
         endif
         myU = newu
         myV = newv

      end subroutine rotate_winds




      subroutine mirror_grid(grid_global,ng,npx,npy,ndims,nregions)
         integer, intent(IN)    :: ng,npx,npy,ndims,nregions
         real   , intent(INOUT) :: grid_global(1-ng:npx  +ng,1-ng:npy  +ng,ndims,1:nregions)
         integer :: i,j,n,n1,n2,nreg
         real :: x1,y1,z1, x2,y2,z2, ang
!
!    Mirror Across the 0-longitude
!
         nreg = 1
         do j=1,ceiling(npy/2.)
            do i=1,ceiling(npx/2.)

            x1 = 0.25 * (ABS(grid_global(i        ,j        ,1,nreg)) + &
                         ABS(grid_global(npx-(i-1),j        ,1,nreg)) + &
                         ABS(grid_global(i        ,npy-(j-1),1,nreg)) + &
                         ABS(grid_global(npx-(i-1),npy-(j-1),1,nreg)))
            grid_global(i        ,j        ,1,nreg) = SIGN(x1,grid_global(i        ,j        ,1,nreg))
            grid_global(npx-(i-1),j        ,1,nreg) = SIGN(x1,grid_global(npx-(i-1),j        ,1,nreg))
            grid_global(i        ,npy-(j-1),1,nreg) = SIGN(x1,grid_global(i        ,npy-(j-1),1,nreg))
            grid_global(npx-(i-1),npy-(j-1),1,nreg) = SIGN(x1,grid_global(npx-(i-1),npy-(j-1),1,nreg))

            y1 = 0.25 * (ABS(grid_global(i        ,j        ,2,nreg)) + &   
                         ABS(grid_global(npx-(i-1),j        ,2,nreg)) + &
                         ABS(grid_global(i        ,npy-(j-1),2,nreg)) + &
                         ABS(grid_global(npx-(i-1),npy-(j-1),2,nreg)))
            grid_global(i        ,j        ,2,nreg) = SIGN(y1,grid_global(i        ,j        ,2,nreg))
            grid_global(npx-(i-1),j        ,2,nreg) = SIGN(y1,grid_global(npx-(i-1),j        ,2,nreg))
            grid_global(i        ,npy-(j-1),2,nreg) = SIGN(y1,grid_global(i        ,npy-(j-1),2,nreg))
            grid_global(npx-(i-1),npy-(j-1),2,nreg) = SIGN(y1,grid_global(npx-(i-1),npy-(j-1),2,nreg))
             
           ! force dateline/greenwich-meridion consitency
            if (mod(npx,2) /= 0) then
              if ( (i==1+(npx-1)/2.0) ) then
                 grid_global(i,j        ,1,nreg) = 0.0
                 grid_global(i,npy-(j-1),1,nreg) = 0.0
              endif
            endif

            enddo
         enddo

         do nreg=2,nregions
           do j=1,npy
             do i=1,npx

               x1 = grid_global(i,j,1,1)
               y1 = grid_global(i,j,2,1)
               z1 = radius

               if (nreg == 2) then
                  ang = -90.
                  call rot_3d( 3, x1, y1, z1, ang, x2, y2, z2, 1, 1)  ! rotate about the z-axis
               elseif (nreg == 3) then
                  ang = -90.
                  call rot_3d( 3, x1, y1, z1, ang, x2, y2, z2, 1, 1)  ! rotate about the z-axis
                  ang = 90.
                  call rot_3d( 1, x2, y2, z2, ang, x1, y1, z1, 1, 1)  ! rotate about the x-axis
                  x2=x1
                  y2=y1
                  z2=z1

           ! force North Pole and dateline/greenwich-meridion consitency
                  if (mod(npx,2) /= 0) then
                     if ( (i==1+(npx-1)/2.0) .and. (i==j) ) then
                        x2 = 0.0
                        y2 = pi/2.0
                     endif
                     if ( (j==1+(npy-1)/2.0) .and. (i < 1+(npx-1)/2.0) ) then
                        x2 = 0.0
                     endif
                     if ( (j==1+(npy-1)/2.0) .and. (i > 1+(npx-1)/2.0) ) then
                        x2 = pi
                     endif
                  endif

               elseif (nreg == 4) then
                  ang = -180.
                  call rot_3d( 3, x1, y1, z1, ang, x2, y2, z2, 1, 1)  ! rotate about the z-axis
                  ang = 90.
                  call rot_3d( 1, x2, y2, z2, ang, x1, y1, z1, 1, 1)  ! rotate about the x-axis
                  x2=x1
                  y2=y1
                  z2=z1

               ! force dateline/greenwich-meridion consitency
                  if (mod(npx,2) /= 0) then
                    if ( (j==1+(npy-1)/2.0) ) then
                       x2 = pi
                    endif
                  endif

               elseif (nreg == 5) then
                  ang = 90.
                  call rot_3d( 3, x1, y1, z1, ang, x2, y2, z2, 1, 1)  ! rotate about the z-axis
                  ang = 90.
                  call rot_3d( 2, x2, y2, z2, ang, x1, y1, z1, 1, 1)  ! rotate about the y-axis
                  x2=x1
                  y2=y1
                  z2=z1
               elseif (nreg == 6) then
                  ang = 90.
                  call rot_3d( 2, x1, y1, z1, ang, x2, y2, z2, 1, 1)  ! rotate about the y-axis
                  ang = 0.
                  call rot_3d( 3, x2, y2, z2, ang, x1, y1, z1, 1, 1)  ! rotate about the z-axis
                  x2=x1
                  y2=y1
                  z2=z1

           ! force South Pole and dateline/greenwich-meridion consitency
                  if (mod(npx,2) /= 0) then
                     if ( (i==1+(npx-1)/2.0) .and. (i==j) ) then
                        x2 = 0.0
                        y2 = -pi/2.0
                     endif
                     if ( (i==1+(npx-1)/2.0) .and. (j > 1+(npy-1)/2.0) ) then
                        x2 = 0.0
                     endif
                     if ( (i==1+(npx-1)/2.0) .and. (j < 1+(npy-1)/2.0) ) then
                        x2 = pi
                     endif
                  endif

               endif

               grid_global(i,j,1,nreg) = x2
               grid_global(i,j,2,nreg) = y2

              enddo
            enddo
          enddo

  end subroutine mirror_grid




 subroutine d2a2c(im,jm,km, ifirst,ilast, jfirst,jlast, ng, &
                  u,v, ua,va, uc,vc)

! Input
  integer, intent(IN) :: im,jm,km
  integer, intent(IN) :: ifirst,ilast
  integer, intent(IN) :: jfirst,jlast
  integer, intent(IN) :: ng
  !real   , intent(in) :: sinlon(im,jm)
  !real   , intent(in) :: coslon(im,jm)
  !real   , intent(in) :: sinl5(im,jm)
  !real   , intent(in) :: cosl5(im,jm)

! Output
 ! real   , intent(inout) ::  u(ifirst-ng:ilast+ng,jfirst-ng:jlast+1+ng)
 ! real   , intent(inout) ::  v(ifirst-ng:ilast+1+ng,jfirst-ng:jlast+ng)
 ! real   , intent(inout) :: ua(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng)
 ! real   , intent(inout) :: va(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng)
 ! real   , intent(inout) :: uc(ifirst-ng:ilast+1+ng,jfirst-ng:jlast+ng)
 ! real   , intent(inout) :: vc(ifirst-ng:ilast+ng,jfirst-ng:jlast+1+ng)

  real   , intent(inout) ::  u(isd:ied,jsd:jed+1) !ifirst-ng:ilast+ng,jfirst-ng:jlast+1+ng)
  real   , intent(inout) ::  v(isd:ied+1,jsd:jed) !ifirst-ng:ilast+1+ng,jfirst-ng:jlast+ng)
  real   , intent(inout) :: ua(isd:ied,jsd:jed)   !ifirst-ng:ilast+ng,jfirst-ng:jlast+ng)
  real   , intent(inout) :: va(isd:ied,jsd:jed)   !(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng)
  real   , intent(inout) :: uc(isd:ied+1,jsd:jed) !(ifirst-ng:ilast+1+ng,jfirst-ng:jlast+ng)
  real   , intent(inout) :: vc(isd:ied,jsd:jed+1) !(ifirst-ng:ilast+ng,jfirst-ng:jlast+1+ng)

!--------------------------------------------------------------
! Local 

  real   :: sinlon(im,jm)
  real   :: coslon(im,jm)
  real   :: sinl5(im,jm)
  real   :: cosl5(im,jm)

    real :: tmp1(jsd:jed+1)
    real :: tmp2(jsd:jed)
    real :: tmp3(jsd:jed)

    real  mag,mag1,mag2, ang,ang1,ang2 
    real  us, vs, un, vn
    integer i, j, k, im2
    integer js1g1
    integer js2g1
    integer js2g2
    integer js2gc
    integer js2gc1
    integer js2gcp1
    integer js2gd
    integer jn2gc
    integer jn1g1
    integer jn1g2
    integer jn2gd
    integer jn2gsp1

 if (cubed_sphere) then

    call dtoa( u, v,ua,va,im,jm,ng)
    call fill_corners(ua, va, im, jm, VECTOR=.true., AGRID=.true.)
    call atoc(ua,va,uc,vc,im,jm,ng, noComm=.true.)
    call fill_corners(uc, vc, im, jm, VECTOR=.true., CGRID=.true.)

 else  ! Lat-Lon

    im2 = im/2

! Set loop limits

    js1g1   = jfirst-1
    js2g1   = jfirst-1
    js2g2   = jfirst-2
    js2gc   = jfirst-ng
    js2gcp1 = jfirst-ng-1
    js2gd   = jfirst-ng
    jn1g1   = jlast+1
    jn1g2   = jlast+2
    jn2gc   = jlast+ng
    jn2gd   = jlast+ng-1
    jn2gsp1 = jlast+ng-1

    if (have_south_pole) then
       js1g1   = 1
       js2g1   = 2
       js2g2   = 2
       js2gc   = 2
       js2gcp1 = 2   ! NG-1 latitudes on S (starting at 2)
       js2gd   = 2
    endif
    if (have_north_pole) then
       jn1g1   = jm
       jn1g2   = jm
       jn2gc   = jm-1  ! NG latitudes on N (ending at jm-1)
       jn2gd   = jm-1
       jn2gsp1 = jm-1
    endif
!
! Treat the special case of ng = 1
!
    if ( ng == 1 .AND. ng > 1 ) THEN
        js2gc1 = js2gc
    else
        js2gc1 = jfirst-ng+1
        if (have_south_pole) js2gc1 = 2  ! NG-1 latitudes on S (starting at 2)
    endif

  do k=1,km

       if ((have_south_pole) .or. (have_north_pole)) then
! Get D-grid V-wind at the poles.
          call vpol5(u(1:im,:), v(1:im,:), im, jm,            &
                     coslon, sinlon, cosl5, sinl5, ng, ng, jfirst, jlast )
          call mp_ghost_ew(im,jm,1,1, ifirst,ilast, jfirst,jlast, 1,1, ng,ng, ng,ng, v(:,:))
       endif

       call dtoa(u, v, ua, va, im, jm, ng)
       call fill_corners(ua, va, im, jm, VECTOR=.true., AGRID=.true.)

       if ( have_south_pole ) then
! Projection at SP
          us = 0.
          vs = 0.
          do i=1,im2
            us = us + (ua(i+im2,2)-ua(i,2))*sinlon(i,2)         &
                    + (va(i,2)-va(i+im2,2))*coslon(i,2)
            vs = vs + (ua(i+im2,2)-ua(i,2))*coslon(i,2)         &
                    + (va(i+im2,2)-va(i,2))*sinlon(i,2)
          enddo
          us = us/im
          vs = vs/im
! SP
          do i=1,im2
            ua(i,1)  = -us*sinlon(i,1) - vs*coslon(i,1)
            va(i,1)  =  us*coslon(i,1) - vs*sinlon(i,1)
            ua(i+im2,1)  = -ua(i,1)
            va(i+im2,1)  = -va(i,1)
          enddo
          ua(0   ,1) = ua(im,1)
          ua(im+1,1) = ua(1 ,1)
          va(im+1,1) = va(1 ,1)
        endif

        if ( have_north_pole ) then
! Projection at NP
          un = 0.
          vn = 0.
          j = jm-1
          do i=1,im2
            un = un + (ua(i+im2,j)-ua(i,j))*sinlon(i,j)        &
                    + (va(i+im2,j)-va(i,j))*coslon(i,j)
            vn = vn + (ua(i,j)-ua(i+im2,j))*coslon(i,j)        &
                    + (va(i+im2,j)-va(i,j))*sinlon(i,j)
          enddo
          un = un/im
          vn = vn/im
! NP
          do i=1,im2
            ua(i,jm) = -un*sinlon(i,jm) + vn*coslon(i,jm)
            va(i,jm) = -un*coslon(i,jm) - vn*sinlon(i,jm)
            ua(i+im2,jm) = -ua(i,jm)
            va(i+im2,jm) = -va(i,jm)
          enddo
          ua(0   ,jm) = ua(im,jm)
          ua(im+1,jm) = ua(1 ,jm)
          va(im+1,jm) = va(1 ,jm)
        endif

        if (latlon) call mp_ghost_ew(im,jm,1,1, ifirst,ilast, jfirst,jlast, 1,1, ng,ng, ng,ng, ua(:,:))
        if (latlon) call mp_ghost_ew(im,jm,1,1, ifirst,ilast, jfirst,jlast, 1,1, ng,ng, ng,ng, va(:,:))

! A -> C
        call atoc(ua, va, uc, vc, im, jm, ng, noComm=.true.)

     enddo ! km loop

     call fill_corners(uc, vc, im, jm, VECTOR=.true., CGRID=.true.)
   endif


 end subroutine d2a2c

!----------------------------------------------------------------------- 
!----------------------------------------------------------------------- 
!BOP
!
 subroutine vpol5(u, v, im, jm, coslon, sinlon, cosl5, sinl5,    &
                  ng_d,  ng_s,  jfirst, jlast)

! !INPUT PARAMETERS:
      integer im                       ! Total longitudes
      integer jm                       ! Total latitudes
      integer jfirst                   ! First PE latitude (no ghosting)
      integer jlast                    ! Last  PE latitude (no ghosting)
      integer, intent(in):: ng_s, ng_d
      real, intent(in):: coslon(im,jm), sinlon(im,jm)
      real, intent(in):: cosl5(im,jm),sinl5(im,jm)
      real, intent(in):: u(im,jfirst-ng_d:jlast+ng_s)

! !INPUT/OUTPUT PARAMETERS:
      real, intent(inout):: v(im,jfirst-ng_d:jlast+ng_d)

! !DESCRIPTION:
!
!   Treat the V winds at the poles.  This requires an average 
!   of the U- and V-winds, weighted by their angles of incidence
!   at the pole points.     
!
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      integer i, imh
      real  uanp(im), uasp(im), vanp(im), vasp(im)
      real  un, vn, us, vs, r2im

! WS 99.05.25 :  Replaced conversions of IMR with IM
      r2im = 0.5d0/dble(im)
      imh  = im / 2

! WS 990726 :  Added condition to decide if poles are on this processor

   if ( jfirst-ng_d <= 1 ) then
         do i=1,im
            uasp(i) = u(i,  2) + u(i,3)
         enddo

         do i=1,im-1
            vasp(i)  = v(i,  2) + v(i+1,2)
         enddo
            vasp(im) = v(im,2) + v(1,2)

! Projection at SP
      us = 0.; vs = 0.

      do i=1,imh
         us = us + (uasp(i+imh)-uasp(i))*sinlon(i,1)    &
                 + (vasp(i)-vasp(i+imh))*coslon(i,1)
         vs = vs + (uasp(i+imh)-uasp(i))*coslon(i,1)    &
                 + (vasp(i+imh)-vasp(i))*sinlon(i,1)
      enddo
      us = us*r2im
      vs = vs*r2im

! get V-wind at SP

      do i=1,imh
         v(i,    1) =  us*cosl5(i,1) - vs*sinl5(i,1)
         v(i+imh,1) = -v(i,1)
      enddo

   endif

   if ( jlast+ng_d >= jm ) then

      do i=1,im
         uanp(i) = u(i,jm-1) + u(i,jm)
      enddo

      do i=1,im-1
         vanp(i) = v(i,jm-1) + v(i+1,jm-1)
      enddo
         vanp(im) = v(im,jm-1) + v(1,jm-1)

! Projection at NP

      un = 0.
      vn = 0.
      do i=1,imh
         un = un + (uanp(i+imh)-uanp(i))*sinlon(i,jm)   &
                 + (vanp(i+imh)-vanp(i))*coslon(i,jm)
         vn = vn + (uanp(i)-uanp(i+imh))*coslon(i,jm)   &
                 + (vanp(i+imh)-vanp(i))*sinlon(i,jm)
      enddo
      un = un*r2im
      vn = vn*r2im

! get V-wind at NP

      do i=1,imh
         v(i,    jm) = -un*cosl5(i,jm) - vn*sinl5(i,jm)
         v(i+imh,jm) = -v(i,jm)
      enddo

   endif

 end subroutine vpol5

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_ghost_ew --- Ghost 4d east/west "lat/lon periodic
!
! !INTERFACE:
      subroutine mp_ghost_ew(im, jm, km, nq, ifirst, ilast, jfirst, jlast, &
                              kfirst, klast, ng_w, ng_e, ng_s, ng_n, q_ghst, q)
!
! !INPUT PARAMETERS:
      integer, intent(in):: im, jm, km, nq
      integer, intent(in):: ifirst, ilast
      integer, intent(in):: jfirst, jlast
      integer, intent(in):: kfirst, klast
      integer, intent(in):: ng_e      ! eastern  zones to ghost
      integer, intent(in):: ng_w      ! western  zones to ghost
      integer, intent(in):: ng_s      ! southern zones to ghost
      integer, intent(in):: ng_n      ! northern zones to ghost
      real, intent(inout):: q_ghst(ifirst-ng_w:ilast+ng_e,jfirst-ng_s:jlast+ng_n,kfirst:klast,nq)
      real, optional, intent(in):: q(ifirst:ilast,jfirst:jlast,kfirst:klast,nq)
!
! !DESCRIPTION:
!
!     Ghost 4d east/west 
!
! !REVISION HISTORY:
!    2005.08.22   Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: i,j,k,n

      if (present(q)) then
         q_ghst(ifirst:ilast,jfirst:jlast,kfirst:klast,1:nq) = &
              q(ifirst:ilast,jfirst:jlast,kfirst:klast,1:nq)
      endif

!      Assume Periodicity in X-dir and not overlapping
      do n=1,nq
         do k=kfirst,klast
            do j=jfirst-ng_s,jlast+ng_n
               do i=1, ng_w
                  q_ghst(ifirst-i,j,k,n) = q_ghst(ilast-i+1,j,k,n)
               enddo
               do i=1, ng_e
                  q_ghst(ilast+i,j,k,n) = q_ghst(ifirst+i-1,j,k,n)
               enddo
            enddo
         enddo
      enddo

!EOC
      end subroutine mp_ghost_ew





      subroutine unit_vect2( p1, p2, uvect )
! No normal projection version
      real, intent(in):: p1(2), p2(2)        ! input position unit vectors (spherical coordinates)
      real, intent(out):: uvect(3)           ! output unit vspherical cartesian
! local        
      integer :: n
      real :: xyz1(3), xyz2(3)

      call spherical_to_cartesian(p1(1), p1(2), 1.0, xyz1(1), xyz1(2), xyz1(3))
      call spherical_to_cartesian(p2(1), p2(2), 1.0, xyz2(1), xyz2(2), xyz2(3))
      do n=1,3
         uvect(n) = xyz2(n)-xyz1(n)
      enddo
      call normalize_vect(1, uvect)

      end subroutine unit_vect2


 subroutine get_unit_vector_3pts( p1, p2, p3, uvect )
 real, intent(in):: p1(2), p2(2), p3(2) ! input position unit vectors (spherical coordinates)
 real, intent(out):: uvect(3)           ! output unit vspherical cartesian
! local
 integer :: n 
 real :: xyz1(3), xyz2(3), xyz3(3)
 real :: dp(3) 
 real :: dp_dot_p2

  call spherical_to_cartesian(p1(1), p1(2), 1.0, xyz1(1), xyz1(2), xyz1(3))
  call spherical_to_cartesian(p2(1), p2(2), 1.0, xyz2(1), xyz2(2), xyz2(3))
  call spherical_to_cartesian(p3(1), p3(2), 1.0, xyz3(1), xyz3(2), xyz3(3))
  do n=1,3
     uvect(n) = xyz3(n)-xyz1(n)
  enddo
  call project_sphere_v(1, uvect,xyz2)
  call normalize_vect(1, uvect)

 end subroutine get_unit_vector_3pts


 subroutine get_unit_vector_2pts( p1, p2, uvect )
 real, intent(in):: p1(2), p2(2)        ! input position unit vectors (spherical coordinates)
 real, intent(out):: uvect(3)           ! output unit vspherical cartesian
! local        
 integer :: n 
 real :: xyz1(3), xyz2(3)         
 real :: dp_dot_xyz1
                  
  call spherical_to_cartesian(p1(1), p1(2), 1.0, xyz1(1), xyz1(2), xyz1(3))
  call spherical_to_cartesian(p2(1), p2(2), 1.0, xyz2(1), xyz2(2), xyz2(3))
  do n=1,3                 
     uvect(n) = xyz2(n)-xyz1(n)   
  enddo                 
  call project_sphere_v(1, uvect,xyz1)
  call normalize_vect(1, uvect)

 end subroutine get_unit_vector_2pts




 subroutine normalize_vect(np, e)
!
! Make e an unit vector
!
 implicit none
 integer, intent(in):: np
 real, intent(inout):: e(3,np)
! local:
 integer k, n
 real pdot

 do n=1,np
    pdot = sqrt(e(1,n)**2+e(2,n)**2+e(3,n)**2)
    do k=1,3
       e(k,n) = e(k,n) / pdot
    enddo
 enddo

 end subroutine normalize_vect

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!     
!     interp_left_edge_1d :: interpolate to left edge of a cell either
!               order = 1 -> Linear average
!               order = 2 -> Uniform PPM
!               order = 3 -> Non-Uniform PPM  
!
 subroutine interp_left_edge_1d(qout, qin, dx, ifirst, ilast, order)
 integer, intent(in):: ifirst,ilast
 real, intent(out)  :: qout(ifirst:)
 real, intent(in)   ::  qin(ifirst:)
 real, intent(in)   ::   dx(ifirst:)
 integer, intent(in):: order
 integer :: i

 real :: dm(ifirst:ilast),qmax,qmin
 real :: r3, da1, da2, a6da, a6, al, ar  
 real :: qLa, qLb1, qLb2
 real :: x

 r3 = 1./3.

 qout(:) = 0.0 
 if (order==1) then 
! 1st order Uniform linear averaging
    do i=ifirst+1,ilast
       qout(i) = 0.5 * (qin(i-1) + qin(i))
    enddo
 elseif (order==2) then
! Non-Uniform 1st order average 
    do i=ifirst+1,ilast
       qout(i) = (dx(i-1)*qin(i-1) + dx(i)*qin(i))/(dx(i-1)+dx(i))
    enddo
 elseif (order==3) then 

! PPM - Uniform 
    do i=ifirst+1,ilast-1
       dm(i) = 0.25*(qin(i+1) - qin(i-1))
    enddo
!
! Applies monotonic slope constraint
!
     do i=ifirst+1,ilast-1
        qmax = max(qin(i-1),qin(i),qin(i+1)) - qin(i)
        qmin = qin(i) - min(qin(i-1),qin(i),qin(i+1))
        dm(i) = sign(min(abs(dm(i)),qmin,qmax),dm(i))
     enddo

     do i=ifirst+1,ilast-1
         qout(i) = 0.5*(qin(i-1)+qin(i)) + r3*(dm(i-1) - dm(i))
       ! al = 0.5*(qin(i-1)+qin(i)) + r3*(dm(i-1) - dm(i))
       ! da1 = dm(i) + dm(i)
       ! qout(i) = qin(i) - sign(min(abs(da1),abs(al-qin(i))), da1)
     enddo

! First order average to fill in end points
     qout(ifirst+1) = 0.5 * (qin(ifirst) + qin(ifirst+1))
     qout(ilast) = 0.5 * (qin(ilast-1) + qin(ilast))

 elseif (order==4) then

  ! Non-Uniform PPM
     do i=ifirst+1,ilast-1
        dm(i) = ( (2.*dx(i-1) + dx(i) ) /                         &
                  (   dx(i+1) + dx(i) )  )  * ( qin(i+1) - qin(i) ) + &
                ( (dx(i)   + 2.*dx(i+1)) /                        &
                  (dx(i-1) +    dx(i)  )  ) * ( qin(i) - qin(i-1) )
        dm(i) = ( dx(i) / ( dx(i-1) + dx(i) + dx(i+1) ) ) * dm(i)
        if ( (qin(i+1)-qin(i))*(qin(i)-qin(i-1)) > 0.) then
           dm(i) = SIGN( MIN( ABS(dm(i)), 2.*ABS(qin(i)-qin(i-1)), 2.*ABS(qin(i+1)-qin(i)) ) , dm(i) )
        else
           dm(i) = 0.
        endif
     enddo

     do i=ifirst+2,ilast-1
        qLa = ( (dx(i-2) + dx(i-1)) / (2.*dx(i-1) +  dx(i)) ) - &
              ( (dx(i+1) + dx(i)) / (2.*dx(i) +  dx(i-1)) )
        qLa = ( (2.*dx(i) * dx(i-1))  / (dx(i-1) + dx(i)) ) * qLa * &
                (qin(i) - qin(i-1))
        qLb1 = dx(i-1) * ( (dx(i-2) + dx(i-1)) / (2.*dx(i-1) + dx(i)) ) * &
              dm(i)
        qLb2 = dx(i) * ( (dx(i) + dx(i+1)) / (dx(i-1) + 2.*dx(i)) ) * &
              dm(i-1)

        qout(i) = 1. / ( dx(i-2) + dx(i-1) + dx(i) + dx(i+1) )
        qout(i) = qout(i) * ( qLa - qLb1 + qLb2 )
        qout(i) = qin(i-1) + ( dx(i-1) / ( dx(i-1) + dx(i) ) ) * (qin(i) - qin(i-1)) + qout(i)
     enddo

 elseif (order==5) then
  
     ! Linear Spline
    do i=ifirst+1,ilast-1
       x = FLOAT(i-(ifirst+1))*FLOAT(ilast-ifirst+1-1)/FLOAT(ilast-ifirst-1) 
       qout(i) = qin(ifirst+NINT(x)) + (x - NINT(x)) * (qin(ifirst+NINT(x+1)) - qin(ifirst+NINT(x)))
      ! if (tile==1) print*, ifirst+NINT(x+1), ifirst+NINT(x), (x - NINT(x)) 
      ! if (tile==1) print*, 0.5*(qin(i-1)+qin(i)), qout(i)
    enddo

   if (tile==1) print*,'x=fltarr(28)'
    do i=ifirst,ilast
       if (tile==1) print*, 'x(',i-ifirst,')=',qin(i)
    enddo


	call mp_stop
	stop

 endif

 end subroutine interp_left_edge_1d
!------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     gsum :: get global sum
!
      real  function globalsum(p, npx, npy, ifirst, ilast, jfirst, jlast) result (gsum)
             
         integer,   intent(IN)    :: npx, npy
         integer,   intent(IN)    :: ifirst, ilast
         integer,   intent(IN)    :: jfirst, jlast
         real  , intent(IN)    :: p(ifirst:ilast,jfirst:jlast)      ! field to be summed
         
         integer :: i,j,k,n
         integer :: j1, j2
         real  :: gsum0
         real :: p_R8(npx-1,npy-1,ntiles_g)
         
         gsum = 0.
            
         if (latlon) then          
            j1 = 2                          
            j2 = npy-2
            gsum = gsum + p(1,1)*acapS
            gsum = gsum + p(1,npy-1)*acapN
            do j=j1,j2
               do i=1,npx-1
                  gsum = gsum + p(i,j)*cos(agrid(i,j,2))
               enddo
            enddo
         else

            do n=tile,tile            
               do j=jfirst,jlast
                  do i=ifirst,ilast
                     p_R8(i,j,n) = p(i,j)*area(i,j)
                  enddo
               enddo
            enddo
            call mp_gather(p_R8, ifirst,ilast, jfirst,jlast, npx-1, npy-1, ntiles_g)
            if (gid == masterproc) then
               do n=1,ntiles_g
                  do j=1,npy-1
                     do i=1,npx-1
                        gsum = gsum + p_R8(i,j,n)
                     enddo
                  enddo
               enddo
               gsum = gsum/globalarea
            endif
            call mpp_broadcast(gsum, masterproc)

         endif

      end function globalsum
 

      end module fv_grid_tools_mod

