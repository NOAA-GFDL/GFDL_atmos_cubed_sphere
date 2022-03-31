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
module fv_grid_tools_mod

! <table>
! <tr>
!     <th>Module Name</th>
!     <th>Functions Included</th>
!   </tr>
!   <tr>
!     <td>constants_mod</td>
!     <td>grav, omega, pi=>pi_8, cnst_radius=>radius</td>
!   </tr>
!   <tr>
!     <td>fms_mod</td>
!     <td>get_mosaic_tile_grid</td>
!   </tr>
!   <tr>
!     <td>fms_io_mod</td>
!     <td>file_exist, field_exist, read_data, get_global_att_value, get_var_att_value</td>
!   </tr>
!   <tr>
!     <td>fv_arrays_mod</td>
!     <td>fv_atmos_type, fv_grid_type, fv_grid_bounds_type, R_GRID</td>
!   </tr>
!   <tr>
!     <td>fv_control_mod</td>
!     <td>fv_init, fv_end, ngrids</td>
!   </tr>
!   <tr>
!     <td>fv_diagnostics_mod</td>
!     <td>prt_maxmin, prt_gb_nh_sh, prt_height</td>
!   </tr>
!   <tr>
!     <td>fv_eta_mod</td>
!     <td>set_eta, set_external_eta</td>
!   </tr>
!   <tr>
!     <td>fv_fill_mod</td>
!     <td>fillz</td>
!   </tr>
!   <tr>
!     <td>fv_grid_utils_mod</td>
!     <td>gnomonic_grids, great_circle_dist,mid_pt_sphere, spherical_angle,
!        cell_center2, get_area, inner_prod, fill_ghost, direct_transform,
!        dist2side_latlon,spherical_linear_interpolation, big_number</td>
!   </tr>
!   <tr>
!     <td>fv_io_mod</td>
!     <td>fv_io_read_tracers</td>
!   </tr>
!   <tr>
!     <td>fv_mp_mod</td>
!     <td>ng, is_master, fill_corners, XDir, YDir,mp_gather,
!         mp_bcst, mp_reduce_max, mp_stop </td>
!   </tr>
!   <tr>
!     <td>fv_timing_mod</td>
!     <td>timing_on, timing_off</td>
!   </tr>
!   <tr>
!     <td>mosaic_mod</td>
!     <td>get_mosaic_ntiles</td>
!   </tr>
!   <tr>
!     <td>mpp_mod</td>
!     <td>mpp_error, FATAL, get_unit, mpp_chksum, mpp_pe, stdout,
!         mpp_send, mpp_recv, mpp_sync_self, EVENT_RECV, mpp_npes,
!         mpp_sum, mpp_max, mpp_min, mpp_root_pe, mpp_broadcast, mpp_transmit</td>
!   </tr>
!   <tr>
!     <td>mpp_domains_mod</td>
!     <td>domain2d,mpp_update_domains, mpp_get_boundary,mpp_get_ntile_count,
!         mpp_get_pelist, mpp_get_compute_domains, mpp_global_field,
!         mpp_get_data_domain, mpp_get_compute_domain,mpp_get_global_domain,
!         mpp_global_sum, mpp_global_max, mpp_global_min</td>
!   </tr>
!   <tr>
!     <td>mpp_io_mod</td>
!     <td>mpp_get_att_value</td>
!   </tr>
!   <tr>
!     <td>mpp_parameter_mod</td>
!     <td>AGRID_PARAM=>AGRID,DGRID_NE_PARAM=>DGRID_NE,
!         CGRID_NE_PARAM=>CGRID_NE,CGRID_SW_PARAM=>CGRID_SW,
!         BGRID_NE_PARAM=>BGRID_NE,BGRID_SW_PARAM=>BGRID_SW,
!         SCALAR_PAIR,CORNER, CENTER, XUPDATE</td>
!   </tr>
!   <tr>
!     <td>sorted_index_mod</td>
!     <td>sorted_inta, sorted_intb</td>
!   </tr>
!   <tr>
!     <td>tracer_manager_mod</td>
!     <td>get_tracer_names, get_number_tracers, get_tracer_index, set_tracer_profile</td>
!   </tr>
! </table>


  use constants_mod,     only: grav, omega, pi=>pi_8, cnst_radius=>radius, small_fac
  use fms_mod,           only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE, clock_flag_default
  use fv_arrays_mod,     only: fv_atmos_type, fv_grid_type, fv_grid_bounds_type, R_GRID
  use fv_grid_utils_mod, only: gnomonic_grids, great_circle_dist, &
                               mid_pt_sphere, spherical_angle, &
                               cell_center2, get_area, inner_prod, fill_ghost, &
                               direct_transform, cube_transform, dist2side_latlon, &
                               spherical_linear_interpolation, big_number
  use fv_timing_mod,     only: timing_on, timing_off
  use fv_mp_mod,         only: is_master, fill_corners, XDir, YDir
  use fv_mp_mod,         only: mp_bcst, mp_reduce_max, mp_stop, grids_master_procs
  use sorted_index_mod,  only: sorted_inta, sorted_intb
  use mpp_mod,           only: mpp_error, FATAL, get_unit, mpp_chksum, mpp_pe, stdout, &
                               mpp_send, mpp_recv, mpp_sync_self, EVENT_RECV, mpp_npes, &
                               mpp_sum, mpp_max, mpp_min, mpp_root_pe, mpp_broadcast, mpp_gather
  use mpp_domains_mod,   only: mpp_update_domains, mpp_get_boundary, &
                               mpp_get_ntile_count, mpp_get_pelist, &
                               mpp_get_compute_domains, mpp_global_field, &
                               mpp_get_data_domain, mpp_get_compute_domain, &
                               mpp_get_global_domain, mpp_global_sum, mpp_global_max, mpp_global_min
  use mpp_domains_mod,   only: domain2d

  use mpp_parameter_mod, only: AGRID_PARAM=>AGRID,       &
                               DGRID_NE_PARAM=>DGRID_NE, &
                               CGRID_NE_PARAM=>CGRID_NE, &
                               CGRID_SW_PARAM=>CGRID_SW, &
                               BGRID_NE_PARAM=>BGRID_NE, &
                               BGRID_SW_PARAM=>BGRID_SW, &
                               SCALAR_PAIR,              &
                               CORNER, CENTER, XUPDATE
  use fms2_io_mod,       only: file_exists, variable_exists, open_file, read_data, &
                               get_global_attribute, get_variable_attribute, &
                               close_file, get_mosaic_tile_grid, FmsNetcdfFile_t
  use mosaic_mod,        only: get_mosaic_ntiles

  implicit none
  private
#include <netcdf.inc>

  real(kind=R_GRID), parameter:: radius = cnst_radius

  real(kind=R_GRID), parameter:: todeg = 180.0d0/pi          !< convert to degrees
  real(kind=R_GRID), parameter:: torad = pi/180.0d0          !< convert to radians
  real(kind=R_GRID), parameter:: missing = 1.d25

  real(kind=R_GRID) :: csFac

  logical, parameter :: debug_message_size = .false.
  logical :: write_grid_char_file = .false.


  public :: todeg, missing, init_grid, spherical_to_cartesian

contains

!>@brief The subroutine 'read_grid' reads the grid from the mosaic grid file.
  subroutine read_grid(Atm, grid_file, ndims, nregions, ng)
    !                  only reads in the grid CORNERS; other metrics (agrid, dx, dy, etc.)
    !                  still need to be computed
    type(fv_atmos_type), intent(inout), target :: Atm
    character(len=*),    intent(IN)    :: grid_file
    integer,             intent(IN)    :: ndims
    integer,             intent(IN)    :: nregions
    integer,             intent(IN)    :: ng

    type(FmsNetcdfFile_t) :: Grid_input
    real, allocatable, dimension(:,:)  :: tmpx, tmpy
    real(kind=R_GRID), pointer, dimension(:,:,:)    :: grid
    character(len=128)                 :: units = ""
    character(len=256)                 :: atm_mosaic, atm_hgrid, grid_form
    character(len=1024)                :: attvalue
    integer                            :: ntiles, i, j, stdunit
    integer                            :: isc2, iec2, jsc2, jec2
    integer                            :: start(4), nread(4)
    integer                            :: is,  ie,  js,  je
    integer                            :: isd, ied, jsd, jed
      integer,save :: halo=3 ! for regional domain external tools

    is  = Atm%bd%is
    ie  = Atm%bd%ie
    js  = Atm%bd%js
    je  = Atm%bd%je
    isd = Atm%bd%isd
    ied = Atm%bd%ied
    jsd = Atm%bd%jsd
    jed = Atm%bd%jed
    grid  => Atm%gridstruct%grid_64

    if(.not. file_exists(grid_file)) call mpp_error(FATAL, 'fv_grid_tools(read_grid): file '// &
         trim(grid_file)//' does not exist')

    !--- make sure the grid file is mosaic file.
    if( open_file(Grid_input, grid_file, "read") ) then
       if( variable_exists(Grid_input, 'atm_mosaic_file') .OR. variable_exists(Grid_input, 'gridfiles') ) then
          stdunit = stdout()
          write(stdunit,*) '==>Note from fv_grid_tools_mod(read_grid): read atmosphere grid from mosaic version grid'
       else
          call mpp_error(FATAL, 'fv_grid_tools(read_grid): neither atm_mosaic_file nor gridfiles exists in file ' &
               //trim(grid_file))
       endif

       if(variable_exists(Grid_input, 'atm_mosaic_file') ) then
          call read_data(Grid_input, "atm_mosaic_file", atm_mosaic)
          atm_mosaic = "INPUT/"//trim(atm_mosaic)
       else
          atm_mosaic = trim(grid_file)
       endif
       call close_file(Grid_input)
    endif

    call get_mosaic_tile_grid(atm_hgrid, atm_mosaic, Atm%domain)

    grid_form = "none"
    if (open_file(Grid_input, atm_hgrid, "read")) then
       call get_global_attribute(Grid_input, "history", attvalue)
       if( index(attvalue, "gnomonic_ed") > 0) grid_form = "gnomonic_ed"
    if(grid_form .NE. "gnomonic_ed") call mpp_error(FATAL, &
         "fv_grid_tools(read_grid): the grid should be 'gnomonic_ed' when reading from grid file, contact developer")

    ntiles = get_mosaic_ntiles(atm_mosaic)
    if( .not. Atm%gridstruct%bounded_domain) then  !<-- The regional setup has only 1 tile so do not shutdown in that case.
       if(ntiles .NE. 6) call mpp_error(FATAL, &
            'fv_grid_tools(read_grid): ntiles should be 6 in mosaic file '//trim(atm_mosaic) )
       if(nregions .NE. 6) call mpp_error(FATAL, &
            'fv_grid_tools(read_grid): nregions should be 6 when reading from mosaic file '//trim(grid_file) )
    endif

       call get_variable_attribute(Grid_input, 'x', 'units', units)

       !--- get the geographical coordinates of super-grid.
       isc2 = 2*is-1; iec2 = 2*ie+1
       jsc2 = 2*js-1; jec2 = 2*je+1
       if( Atm%gridstruct%bounded_domain ) then
         isc2 = 2*(isd+halo)-1; iec2 = 2*(ied+1+halo)-1   ! For the regional domain the cell corner locations must be transferred
         jsc2 = 2*(jsd+halo)-1; jec2 = 2*(jed+1+halo)-1   ! from the entire supergrid to the compute grid, including the halo region.
       endif
       allocate(tmpx(isc2:iec2, jsc2:jec2) )
       allocate(tmpy(isc2:iec2, jsc2:jec2) )
       start = 1; nread = 1
       start(1) = isc2; nread(1) = iec2 - isc2 + 1
       start(2) = jsc2; nread(2) = jec2 - jsc2 + 1
       call read_data(Grid_input, 'x', tmpx, corner=start, edge_lengths=nread)  !<-- tmpx (lon, deg east) is on the supergrid
       call read_data(Grid_input, 'y', tmpy, corner=start, edge_lengths=nread)  !<-- tmpy (lat, deg) is on the supergrid
       call close_file(Grid_input)
    endif

    !--- geographic grid at cell corner
    grid(isd: is-1, jsd:js-1,1:ndims)=0.
    grid(isd: is-1, je+2:jed+1,1:ndims)=0.
    grid(ie+2:ied+1,jsd:js-1,1:ndims)=0.
    grid(ie+2:ied+1,je+2:jed+1,1:ndims)=0.
    if(len_trim(units) < 6) call mpp_error(FATAL, &
          "fv_grid_tools_mod(read_grid): the length of units must be no less than 6")
    if(units(1:6) == 'degree') then
    if( .not. Atm%gridstruct%bounded_domain) then
       do j = js, je+1
          do i = is, ie+1
             grid(i,j,1) = tmpx(2*i-1,2*j-1)*pi/180.
             grid(i,j,2) = tmpy(2*i-1,2*j-1)*pi/180.
          enddo
       enddo
    else
!
!***  In the regional case the halo surrounding the domain was included in the read.
!***  Transfer the compute and halo regions to the compute grid.
!
       do j = jsd, jed+1
          do i = isd, ied+1
             grid(i,j,1) = tmpx(2*i+halo+2,2*j+halo+2)*pi/180.
             grid(i,j,2) = tmpy(2*i+halo+2,2*j+halo+2)*pi/180.
          enddo
       enddo
    endif
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

    deallocate(tmpx, tmpy)
    nullify(grid)

  end subroutine read_grid



  !#################################################################################
  subroutine get_symmetry(data_in, data_out, ishift, jshift, npes_x, npes_y, domain, tile, npx_g, bd)
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer,                                      intent(in)  :: ishift, jshift, npes_x, npes_y
    real(kind=R_GRID), dimension(bd%is:bd%ie+ishift, bd%js:bd%je+jshift ), intent(in)  :: data_in
    real(kind=R_GRID), dimension(bd%is:bd%ie+jshift, bd%js:bd%je+ishift ), intent(out) :: data_out
    real(kind=R_GRID),    dimension(:), allocatable :: send_buffer
    real(kind=R_GRID),    dimension(:), allocatable :: recv_buffer
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
    type(domain2d)                     :: domain
    integer                            :: tile, npx_g
    integer :: is,  ie,  js,  je
    integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

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

!>@brief The subroutine 'init_grid' reads the grid from the input file
!! and sets up grid descriptors.
  subroutine init_grid(Atm, grid_name, grid_file, npx, npy, npz, ndims, nregions, ng, tile_coarse)
!--------------------------------------------------------
    type(fv_atmos_type), intent(inout), target :: Atm
    character(len=80), intent(IN) :: grid_name
    character(len=120),intent(IN) :: grid_file
    integer,      intent(IN) :: npx, npy, npz
    integer,      intent(IN) :: ndims
    integer,      intent(IN) :: nregions
    integer,      intent(IN) :: ng
    integer,      intent(IN) :: tile_coarse(:)
!--------------------------------------------------------
    real(kind=R_GRID)   ::  xs(npx,npy)
    real(kind=R_GRID)   ::  ys(npx,npy)

    real(kind=R_GRID)  :: dp, dl
    real(kind=R_GRID)  :: x1,x2,y1,y2,z1,z2
    integer :: i,j,k,n,nreg
    integer :: fileLun

    real(kind=R_GRID)  :: p1(3), p2(3), p3(3), p4(3)
    real(kind=R_GRID)  :: dist,dist1,dist2, pa(2), pa1(2), pa2(2), pb(2)
    real(kind=R_GRID)  :: pt(3), pt1(3), pt2(3), pt3(3)
    real(kind=R_GRID) :: ee1(3), ee2(3)

    real(kind=R_GRID)  :: angN,angM,angAV,ang
    real(kind=R_GRID)  :: aspN,aspM,aspAV,asp
    real(kind=R_GRID)  ::  dxN, dxM, dxAV
    real(kind=R_GRID)  :: dx_local, dy_local

    real(kind=R_GRID)  :: vec1(3), vec2(3), vec3(3), vec4(3)
    real(kind=R_GRID)  :: vecAvg(3), vec3a(3), vec3b(3), vec4a(3), vec4b(3)
    real(kind=R_GRID)  :: xyz1(3), xyz2(3)

!    real(kind=R_GRID) :: grid_global(1-ng:npx  +ng,1-ng:npy  +ng,ndims,1:nregions)
    integer :: ios, ip, jp

    integer :: igrid

    integer :: tmplun
    character(len=80) :: tmpFile

    real(kind=R_GRID), dimension(Atm%bd%is:Atm%bd%ie) :: sbuffer, nbuffer
    real(kind=R_GRID), dimension(Atm%bd%js:Atm%bd%je) :: wbuffer, ebuffer

    real(kind=R_GRID), pointer, dimension(:,:,:) :: agrid, grid
    real(kind=R_GRID), pointer, dimension(:,:) :: area, area_c
    real(kind=R_GRID), pointer, dimension(:,:) ::  area_u,  area_v
    real(kind=R_GRID), pointer, dimension(:,:) ::  dx6,  dy6
    real, pointer, dimension(:,:) :: rarea_u, rarea_v
    real, pointer, dimension(:,:) :: rdx6, rdy6

    real(kind=R_GRID), pointer, dimension(:,:) :: sina, cosa, dx, dy, dxc, dyc, dxa, dya
    real, pointer, dimension(:,:,:) :: e1, e2

    real, pointer, dimension(:,:) :: rarea, rarea_c
    real, pointer, dimension(:,:) :: rdx, rdy, rdxc, rdyc, rdxa, rdya

    integer, pointer, dimension(:,:,:) ::  iinta, jinta, iintb, jintb

    real(kind=R_GRID), pointer, dimension(:,:,:,:) :: grid_global

    integer, pointer :: npx_g, npy_g, ntiles_g, tile
    logical, pointer :: sw_corner, se_corner, ne_corner, nw_corner
    logical, pointer :: latlon, cubed_sphere, have_south_pole, have_north_pole, stretched_grid

    type(domain2d), pointer :: domain
    integer :: is,  ie,  js,  je
    integer :: isd, ied, jsd, jed
    integer :: istart, iend, jstart, jend
    integer :: isection_s, isection_e, jsection_s, jsection_e

    !  Setup timing variables

    logical, save       :: first_time = .true.
    integer, save       :: id_timer1, id_timer2, id_timer3, id_timer3a, id_timer4, id_timer5, id_timer6, id_timer7, id_timer8
    logical             :: use_timer = .false.  ! Set to True for detailed performance profiling
    logical             :: debug_log = .false.
    integer             :: this_pe

    this_pe = mpp_pe()

    if (first_time) then
       if (use_timer) then
          id_timer1     = mpp_clock_id ('init_grid Step 1',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
          id_timer2     = mpp_clock_id ('init_grid Step 2',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
          id_timer3     = mpp_clock_id ('init_grid Step 3',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
          id_timer3a    = mpp_clock_id ('init_grid Step 3a',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
          id_timer4     = mpp_clock_id ('init_grid Step 4',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
          id_timer5     = mpp_clock_id ('init_grid Step 5',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
          id_timer6     = mpp_clock_id ('init_grid Step 6',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
          id_timer7     = mpp_clock_id ('init_grid Step 7',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
          id_timer8     = mpp_clock_id ('init_grid Step 8',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
       end if
       first_time = .false.
    end if

    if (use_timer) call mpp_clock_begin (id_timer1)

    is  = Atm%bd%is
    ie  = Atm%bd%ie
    js  = Atm%bd%js
    je  = Atm%bd%je
    isd = Atm%bd%isd
    ied = Atm%bd%ied
    jsd = Atm%bd%jsd
    jed = Atm%bd%jed

    !!! Associate pointers
    agrid => Atm%gridstruct%agrid_64
    grid  => Atm%gridstruct%grid_64

    area   => Atm%gridstruct%area_64
    area_c  => Atm%gridstruct%area_c_64
    rarea   => Atm%gridstruct%rarea
    rarea_c => Atm%gridstruct%rarea_c

!   For MOLECULAR_DIFFUSION
    if ( Atm%flagstruct%molecular_diffusion ) then
       area_u  => Atm%gridstruct%area_u_64
       area_v  => Atm%gridstruct%area_v_64
       dx6     => Atm%gridstruct%dx6_64
       dy6     => Atm%gridstruct%dy6_64
       rarea_u => Atm%gridstruct%rarea_u
       rarea_v => Atm%gridstruct%rarea_v
       rdx6    => Atm%gridstruct%rdx6
       rdy6    => Atm%gridstruct%rdy6
    endif

    sina   => Atm%gridstruct%sina_64
    cosa   => Atm%gridstruct%cosa_64
    dx     => Atm%gridstruct%dx_64
    dy     => Atm%gridstruct%dy_64
    dxc    => Atm%gridstruct%dxc_64
    dyc    => Atm%gridstruct%dyc_64
    dxa    => Atm%gridstruct%dxa_64
    dya    => Atm%gridstruct%dya_64
    rdx    => Atm%gridstruct%rdx
    rdy    => Atm%gridstruct%rdy
    rdxc   => Atm%gridstruct%rdxc
    rdyc   => Atm%gridstruct%rdyc
    rdxa   => Atm%gridstruct%rdxa
    rdya   => Atm%gridstruct%rdya
    e1     => Atm%gridstruct%e1
    e2     => Atm%gridstruct%e2

    if (Atm%neststruct%nested .or. ANY(Atm%neststruct%child_grids)) then
       if (debug_log) print '("[INFO] WDR grid_global => Atm%grid_global in init_grid fv_grid_tools.F90. npe=",I0)', this_pe
        grid_global => Atm%grid_global
    else if( trim(grid_file) .EQ. 'Inline') then
       if (debug_log) print '("[INFO] WDR inline, allocating grid_global in init_grid fv_grid_tools.F90. npe=",I0)', this_pe
       allocate(grid_global(1-ng:npx  +ng,1-ng:npy  +ng,ndims,1:nregions))
    endif

    iinta                         => Atm%gridstruct%iinta
    jinta                         => Atm%gridstruct%jinta
    iintb                         => Atm%gridstruct%iintb
    jintb                         => Atm%gridstruct%jintb
    npx_g                         => Atm%gridstruct%npx_g
    npy_g                         => Atm%gridstruct%npy_g
    ntiles_g                      => Atm%gridstruct%ntiles_g
    sw_corner                     => Atm%gridstruct%sw_corner
    se_corner                     => Atm%gridstruct%se_corner
    ne_corner                     => Atm%gridstruct%ne_corner
    nw_corner                     => Atm%gridstruct%nw_corner
    latlon                        => Atm%gridstruct%latlon
    cubed_sphere                  => Atm%gridstruct%cubed_sphere
    have_south_pole               => Atm%gridstruct%have_south_pole
    have_north_pole               => Atm%gridstruct%have_north_pole
    stretched_grid                => Atm%gridstruct%stretched_grid

    tile                          => Atm%tile_of_mosaic

    domain                        => Atm%domain

    npx_g = npx
    npy_g = npy
    ntiles_g = nregions
    latlon = .false.
    cubed_sphere = .false.

    if ( (Atm%flagstruct%do_schmidt .or. Atm%flagstruct%do_cube_transform) .and. abs(atm%flagstruct%stretch_fac-1.) > 1.E-5 ) then
       stretched_grid = .true.
       if (Atm%flagstruct%do_schmidt .and. Atm%flagstruct%do_cube_transform) then
          call mpp_error(FATAL, ' Cannot set both do_schmidt and do_cube_transform to .true.')
       endif
    endif

    if (use_timer) call mpp_clock_end (id_timer1)

    if (Atm%flagstruct%grid_type>3) then
       if (use_timer) call mpp_clock_begin (id_timer2)

       if (Atm%flagstruct%grid_type == 4) then
          call setup_cartesian(npx, npy, Atm%flagstruct%dx_const, Atm%flagstruct%dy_const, &
               Atm%flagstruct%deglat, Atm%bd)
       elseif (Atm%flagstruct%grid_type == 5) then
          call setup_orthogonal_grid(npx, npy, Atm%bd, grid_file)
       else
          call mpp_error(FATAL, 'init_grid: unsupported grid type')
       endif
       if (use_timer) call mpp_clock_end (id_timer2)

    else
       if (use_timer) call mpp_clock_begin (id_timer3)

          cubed_sphere = .true.

          if (Atm%neststruct%nested) then
             !Read grid if it exists

             if (use_timer) call mpp_clock_begin (id_timer3a)
             if (Atm%flagstruct%grid_type < 0) then
                !Note that read_grid only reads in grid corners. Will still need to compute all other grid metrics.
                !NOTE: cannot currently read in mosaic for both coarse and nested grids simultaneously
                call read_grid(Atm, grid_file, ndims, 1, ng)
             endif
             ! still need to set up weights
             call setup_aligned_nest(Atm)
             if (use_timer) call mpp_clock_end (id_timer3a)

          else
             if(trim(grid_file) .NE. 'Inline' .or. Atm%flagstruct%grid_type < 0) then
                call read_grid(Atm, grid_file, ndims, nregions, ng)

             ! Here if we are reading from grid_spec and the grid has a nest we need to assemble
             ! the global grid array 'grid_global' to be sent at the end of this routine to the nest
                if (ANY(Atm%neststruct%child_grids)) then
                   grid_global(:,:,:,1)=-99999
                   isection_s = is
                   isection_e = ie
                   jsection_s = js
                   jsection_e = je

                   if ( isd < 0 )     isection_s = isd
                   if ( ied > npx-1 ) isection_e = ied
                   if ( jsd < 0 )     jsection_s = jsd
                   if ( jed > npy-1 ) jsection_e = jed
                   ! if there is a nest, we need to setup grid_global on pe master
                   ! to send it to the nest at the end of init_grid
                   call mpp_gather(isection_s,isection_e,jsection_s,jsection_e,atm%pelist, &
                                   grid(isection_s:isection_e,jsection_s:jsection_e,1),grid_global(1-ng:npx+ng,1-ng:npy+ng,1,1),is_master(),ng,ng)
                   call mpp_gather(isection_s,isection_e,jsection_s,jsection_e,atm%pelist, &
                                   grid(isection_s:isection_e,jsection_s:jsection_e,2),grid_global(1-ng:npx+ng,1-ng:npy+ng,2,1),is_master(),ng,ng)
                   !do we need the haloes?!
                   !do j=jsd,jed
                   !do i=isd,ied
                     !grid_global(i,j,1,1)=grid(i,j,1)
                     !grid_global(i,j,2,1)=grid(i,j,2)
                   !enddo
                   !enddo
                   !do j=1,npy
                   !do i=1,npx
                     !call mpp_max(grid_global(i,j,1,1),atm%pelist)
                     !call mpp_max(grid_global(i,j,2,1),atm%pelist)
                   !enddo
                   !enddo
                endif

             else

                if (Atm%flagstruct%grid_type>=0) call gnomonic_grids(Atm%flagstruct%grid_type, npx-1, xs, ys)

                if (is_master()) then

                if (Atm%flagstruct%grid_type>=0) then
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
                      if ( .not. ( Atm%flagstruct%do_schmidt .or. Atm%flagstruct%do_cube_transform) .and. (Atm%flagstruct%shift_fac)>1.E-4 )   &
                           grid_global(i,j,1,n) = grid_global(i,j,1,n) - pi/Atm%flagstruct%shift_fac
!----------------------------------------------------------------------------------------------------
                      if ( grid_global(i,j,1,n) < 0. )              &
                           grid_global(i,j,1,n) = grid_global(i,j,1,n) + 2.*pi
                      if (ABS(grid_global(i,j,1,1)) < 1.d-10) grid_global(i,j,1,1) = 0.0
                      if (ABS(grid_global(i,j,2,1)) < 1.d-10) grid_global(i,j,2,1) = 0.0
                      !Change from Github PR #39 - this changes answers
                      !if (ABS(grid_global(i,j,1,n)) < 1.d-10) grid_global(i,j,1,n) = 0.0
                      !if (ABS(grid_global(i,j,2,n)) < 1.d-10) grid_global(i,j,2,n) = 0.0
                   enddo
                   enddo
                   enddo
                else
                   call mpp_error(FATAL, "fv_grid_tools: reading of ASCII grid files no longer supported")
                endif ! Atm%flagstruct%grid_type>=0

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
                if ( Atm%flagstruct%do_schmidt ) then
                   do n=1,nregions
                      call direct_transform(Atm%flagstruct%stretch_fac, 1, npx, 1, npy, &
                           Atm%flagstruct%target_lon, Atm%flagstruct%target_lat, &
                           n, grid_global(1:npx,1:npy,1,n), grid_global(1:npx,1:npy,2,n))
                   enddo
                elseif ( Atm%flagstruct%do_cube_transform) then
                   do n=1,nregions
                      call cube_transform(Atm%flagstruct%stretch_fac, 1, npx, 1, npy, &
                           Atm%flagstruct%target_lon, Atm%flagstruct%target_lat, &
                           n, grid_global(1:npx,1:npy,1,n), grid_global(1:npx,1:npy,2,n))
                   enddo
                endif
                endif !is master
                call mpp_broadcast(grid_global, size(grid_global), mpp_root_pe())
                !--- copy grid to compute domain
                do n=1,ndims
                do j=js,je+1
                do i=is,ie+1
                   grid(i,j,n) = grid_global(i,j,n,tile)
                enddo
                enddo
                enddo
             endif !(trim(grid_file) == 'INPUT/grid_spec.nc')
!
! SJL: For phys/exchange grid, etc
!
             call mpp_update_domains( grid, Atm%domain, position=CORNER)
             if (.not. (Atm%gridstruct%bounded_domain)) then
                call fill_corners(grid(:,:,1), npx, npy, FILL=XDir, BGRID=.true.)
                call fill_corners(grid(:,:,2), npx, npy, FILL=XDir, BGRID=.true.)
             endif

             !--- dx and dy
             if( .not. Atm%gridstruct%bounded_domain) then
                istart=is
                iend=ie
                jstart=js
                jend=je
             else
                istart=isd
                iend=ied
                jstart=jsd
                jend=jed
             endif

             do j = jstart, jend+1
             do i = istart, iend
                p1(1) = grid(i  ,j,1)
                p1(2) = grid(i  ,j,2)
                p2(1) = grid(i+1,j,1)
                p2(2) = grid(i+1,j,2)
                dx(i,j) = great_circle_dist( p2, p1, radius )
             enddo
             enddo
             if( stretched_grid .or. Atm%gridstruct%bounded_domain ) then
                do j = jstart, jend
                do i = istart, iend+1
                   p1(1) = grid(i,j,  1)
                   p1(2) = grid(i,j,  2)
                   p2(1) = grid(i,j+1,1)
                   p2(2) = grid(i,j+1,2)
                   dy(i,j) = great_circle_dist( p2, p1, radius )
                enddo
                enddo
             else
                call get_symmetry(dx(is:ie,js:je+1), dy(is:ie+1,js:je), 0, 1, Atm%layout(1), Atm%layout(2), &
                     Atm%domain, Atm%tile_of_mosaic, Atm%gridstruct%npx_g, Atm%bd)
             endif

             call mpp_get_boundary( dy, dx, Atm%domain, ebufferx=ebuffer, wbufferx=wbuffer, sbuffery=sbuffer, nbuffery=nbuffer,&
                  flags=SCALAR_PAIR+XUPDATE, gridtype=CGRID_NE_PARAM)
             if( .not. Atm%gridstruct%bounded_domain ) then
                if(is == 1 .AND. mod(tile,2) .NE. 0) then ! on the west boundary
                   dy(is, js:je) = wbuffer(js:je)
                endif
                if(ie == npx-1) then  ! on the east boundary
                   dy(ie+1, js:je) = ebuffer(js:je)
                endif
             endif

             call mpp_update_domains( dy, dx, Atm%domain, flags=SCALAR_PAIR,      &
                  gridtype=CGRID_NE_PARAM, complete=.true.)
             if (cubed_sphere .and. (.not. (Atm%gridstruct%bounded_domain))) then
                call fill_corners(dx, dy, npx, npy, DGRID=.true.)
             endif

       if( .not. stretched_grid )         &
           call sorted_inta(isd, ied, jsd, jed, cubed_sphere, grid, iinta, jinta)

       agrid(:,:,:) = -1.e25

          !--- compute agrid (use same indices as for dx/dy above)

       do j=jstart,jend
          do i=istart,iend
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

             call mpp_update_domains( agrid, Atm%domain, position=CENTER, complete=.true. )
             if (.not. (Atm%gridstruct%bounded_domain)) then
                call fill_corners(agrid(:,:,1), npx, npy, XDir, AGRID=.true.)
                call fill_corners(agrid(:,:,2), npx, npy, YDir, AGRID=.true.)
             endif

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
!      call mpp_update_domains( dxa, dya, Atm%domain, flags=SCALAR_PAIR, gridtype=AGRID_PARAM)
             if (cubed_sphere  .and. (.not. (Atm%gridstruct%bounded_domain))) then
                call fill_corners(dxa, dya, npx, npy, AGRID=.true.)
             endif


          end if !if nested

    if (use_timer) call mpp_clock_end (id_timer3)
    if (use_timer) call mpp_clock_begin (id_timer4)

       do j=jsd,jed
          do i=isd+1,ied
             dxc(i,j) = great_circle_dist(agrid(i,j,:), agrid(i-1,j,:), radius)
          enddo
!xxxxxx
      !Are the following 2 lines appropriate for the regional domain?
!xxxxxx
          dxc(isd,j)   = dxc(isd+1,j)
          dxc(ied+1,j) = dxc(ied,j)
       enddo

!       do j=js,je+1
!          do i=is,ie
       do j=jsd+1,jed
          do i=isd,ied
             dyc(i,j) = great_circle_dist(agrid(i,j,:), agrid(i,j-1,:), radius)
          enddo
       enddo
!xxxxxx
      !Are the following 2 lines appropriate for the regional domain?
!xxxxxx
       do i=isd,ied
          dyc(i,jsd)   = dyc(i,jsd+1)
          dyc(i,jed+1) = dyc(i,jed)
       end do

       if (use_timer) call mpp_clock_end (id_timer4)
       if (use_timer) call mpp_clock_begin (id_timer5)

       if ( Atm%flagstruct%molecular_diffusion ) then
! dx6, dy6
          do j=jsd,jed+1
             do i=isd+1,ied
                call mid_pt_sphere(grid(i-1,j,1:2), grid(i  ,j,1:2), p1)
                call mid_pt_sphere(grid(i  ,j,1:2), grid(i+1,j,1:2), p2)
                dx6(i,j) = great_circle_dist( p2, p1, radius )
             enddo
!xxxxxx
      !Are the following 2 lines appropriate for the regional domain?
!xxxxxx
            dx6(isd  ,j) = dx6(isd+1,j)
            dx6(ied+1,j) = dx6(ied  ,j)
         enddo

         do j=jsd+1,jed
            do i=isd,ied+1
               call mid_pt_sphere(grid(i,j-1,1:2), grid(i,j  ,1:2), p1)
               call mid_pt_sphere(grid(i,j  ,1:2), grid(i,j+1,1:2), p2)
               dy6(i,j) = great_circle_dist( p2, p1, radius )
            enddo
         enddo
!xxxxxx
      !Are the following 2 lines appropriate for the regional domain?
!xxxxxx
         do i=isd,ied+1
            dy6(i,jsd)   = dy6(i,jsd+1)
            dy6(i,jed+1) = dy6(i,jed)
         enddo
! area_u, area_v
         do j=jsd,jed
            do i=isd+1,ied
               call mid_pt_sphere(grid(i-1,j  ,1:2), grid(i  ,j  ,1:2), p1)
               call mid_pt_sphere(grid(i-1,j+1,1:2), grid(i  ,j+1,1:2), p4)
               call mid_pt_sphere(grid(i  ,j  ,1:2), grid(i+1,j  ,1:2), p2)
               call mid_pt_sphere(grid(i  ,j+1,1:2), grid(i+1,j+1,1:2), p3)
               area_v(i,j) = get_area(p1, p4, p2, p3, radius)
            enddo
!xxxxxx
      !Are the following 2 lines appropriate for the regional domain?
!xxxxxx
            area_v(isd  ,j) = area_v(isd+1,j)
            area_v(ied+1,j) = area_v(ied  ,j)
         enddo

         do j=jsd+1,jed
            do i=isd,ied
               call mid_pt_sphere(grid(i  ,j-1,1:2), grid(i  ,j  ,1:2), p1)
               call mid_pt_sphere(grid(i  ,j  ,1:2), grid(i  ,j+1,1:2), p4)
               call mid_pt_sphere(grid(i+1,j-1,1:2), grid(i+1,j  ,1:2), p2)
               call mid_pt_sphere(grid(i+1,j  ,1:2), grid(i+1,j+1,1:2), p3)
               area_u(i,j) = get_area(p1, p4, p2, p3, radius)
            enddo
         enddo
!xxxxxx
      !Are the following 2 lines appropriate for the regional domain?
!xxxxxx
         do i=isd,ied
            area_u(i,jsd)   = area_u(i,jsd+1)
            area_u(i,jed+1) = area_u(i,jed)
         enddo

! we are not really using the outmost, so no need to tune this
! so update and fill corners should be enough

         call mpp_update_domains( dx6, Atm%domain, position=CORNER, complete=.true.)
         call mpp_update_domains( dy6, Atm%domain, position=CORNER, complete=.true.)
         call mpp_update_domains( area_v, area_u, Atm%domain, flags=SCALAR_PAIR, &
                                gridtype=CGRID_NE_PARAM, complete=.true.)

         if (cubed_sphere  .and. (.not. (Atm%neststruct%nested .or. Atm%flagstruct%regional))) then
            call fill_corners( dx6, npx, npy, FILL=XDir, BGRID=.true.)
            call fill_corners( dy6, npx, npy, FILL=XDir, BGRID=.true.)
            call fill_corners( area_v, area_u, npx, npy, CGRID=.true.)
         endif

       endif ! MOLECULAR_DIFFUSION

       if( .not. stretched_grid )      &
           call sorted_intb(isd, ied, jsd, jed, is, ie, js, je, npx, npy, &
                            cubed_sphere, agrid, iintb, jintb)

       call grid_area( npx, npy, ndims, nregions, Atm%gridstruct%bounded_domain, Atm%gridstruct, Atm%domain, Atm%bd )
!      stretched_grid = .false.

       if (use_timer) call mpp_clock_end (id_timer5)

!----------------------------------
! Compute area_c, rarea_c, dxc, dyc
!----------------------------------
       if (use_timer) call mpp_clock_begin (id_timer6)
  if ( .not. stretched_grid .and. (.not. (Atm%gridstruct%bounded_domain))) then
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
       if (use_timer) call mpp_clock_end (id_timer6)
       if (use_timer) call mpp_clock_begin (id_timer7)

       call mpp_update_domains( dxc, dyc, Atm%domain, flags=SCALAR_PAIR,   &
                                gridtype=CGRID_NE_PARAM, complete=.true.)
       if (cubed_sphere  .and. (.not. (Atm%gridstruct%bounded_domain))) then
         call fill_corners(dxc, dyc, npx, npy, CGRID=.true.)
       endif

       call mpp_update_domains( area,   Atm%domain, complete=.true. )


       !Handling outermost ends for area_c
       if (Atm%gridstruct%bounded_domain) then
          if (is == 1) then
             do j=jsd,jed
                area_c(isd,j) = area_c(isd+1,j)
             end do
             if (js == 1)     area_c(isd,jsd) = area_c(isd+1,jsd+1)
             if (js == npy-1) area_c(isd,jed+1) = area_c(isd+1,jed)
          end if
          if (ie == npx-1) then
             do j=jsd,jed
                area_c(ied+1,j) = area_c(ied,j)
             end do
             if (js == 1)     area_c(ied+1,jsd) = area_c(ied,jsd+1)
             if (js == npy-1) area_c(ied+1,jed+1) = area_c(ied,jed)
          end if
          if (js == 1) then
             do i=isd,ied
                area_c(i,jsd) = area_c(i,jsd+1)
             end do
          end if
          if (je == npy-1) then
             do i=isd,ied
                area_c(i,jed+1) = area_c(i,jed)
             end do
          end if
       end if

       call mpp_update_domains( area_c, Atm%domain, position=CORNER, complete=.true.)

       ! Handle corner Area ghosting
       if (cubed_sphere .and. (.not. (Atm%gridstruct%bounded_domain))) then
          call fill_ghost(area, npx, npy, -big_number, Atm%bd)  ! fill in garbage values
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

       if (use_timer) call mpp_clock_end (id_timer7)
       if (use_timer) call mpp_clock_begin (id_timer8)

       if ( Atm%flagstruct%molecular_diffusion ) then
          do j=jsd,jed+1
             do i=isd,ied
                rarea_u(i,j) = 1.0/area_u(i,j)
             enddo
          enddo
          do j=jsd,jed
             do i=isd,ied+1
                rarea_v(i,j) = 1.0/area_v(i,j)
             enddo
          enddo
          do j=jsd,jed+1
             do i=isd,ied+1
                rdx6(i,j) = 1.0/dx6(i,j)
             enddo
          enddo
          do j=jsd,jed+1
             do i=isd,ied+1
                rdy6(i,j) = 1.0/dy6(i,j)
             enddo
          enddo
       endif

200    format(A,f9.2,A,f9.2,A,f9.2)
201    format(A,f9.2,A,f9.2,A,f9.2,A,f9.2)
202    format(A,A,i4.4,A,i4.4,A)

       ! Get and print Grid Statistics
       dxAV =0.0
       angAV=0.0
       aspAV=0.0
       dxN  =  missing
       dxM  = -missing
       angN =  missing
       angM = -missing
       aspN =  missing
       aspM = -missing
       !if (tile == 1) then ! doing a GLOBAL domain search on each grid
          do j=js, je
             do i=is, ie
                if(i>ceiling(npx/2.) .OR. j>ceiling(npy/2.)) cycle
                ang  = get_angle(2, grid(i,j+1,1:2), grid(i,j,1:2), grid(i+1,j,1:2))
                ang  = ABS(90.0 - ang)

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

                asp   = ABS(dx_local/dy_local)
                if (asp < 1.0) asp = 1.0/asp
                aspAV = aspAV + asp
                aspM  = MAX(aspM,asp)
                aspN  = MIN(aspN,asp)
             enddo
          enddo
       !endif
       call mpp_sum(angAv)
       call mpp_sum(dxAV)
       call mpp_sum(aspAV)
       call mpp_max(angM)
       call mpp_min(angN)
       call mpp_max(dxM)
       call mpp_min(dxN)
       call mpp_max(aspM)
       call mpp_min(aspN)

       if( is_master() ) then

          angAV = angAV / ( (ceiling(npy/2.0))*(ceiling(npx/2.0)) - 1 )
          dxAV  = dxAV  / ( (ceiling(npy/2.0))*(ceiling(npx/2.0)) )
          aspAV = aspAV / ( (ceiling(npy/2.0))*(ceiling(npx/2.0)) )
          write(*,*  ) ''
          write(*,*) ' Radius is ', radius, ', omega is ', omega, ' small_fac = ', small_fac
          write(*,*  ) ' Cubed-Sphere Grid Stats : ', npx,'x',npy,'x',nregions
          print*, dxN, dxM, dxAV, dxN, dxM
          write(*,'(A,f11.2,A,f11.2,A,f11.2,A,f11.2)') '      Grid Length               : min: ', dxN,' max: ', dxM,' avg: ', dxAV, ' min/max: ',dxN/dxM
          write(*,'(A,e21.14,A,e21.14,A,e21.14)') '      Deviation from Orthogonal : min: ',angN,' max: ',angM,' avg: ',angAV
          write(*,'(A,e21.14,A,e21.14,A,e21.14)') '      Aspect Ratio              : min: ',aspN,' max: ',aspM,' avg: ',aspAV
          write(*,*  ) ''

       endif
    endif!if gridtype > 3

    !SEND grid global if any child nests
    !Matching receive in setup_aligned_nest
    do n=1,size(Atm%neststruct%child_grids)
       if (Atm%neststruct%child_grids(n) .and. is_master()) then
          !need to get tile_coarse AND determine local number for tile
          if (ntiles_g > 1) then ! coarse grid only!!
             call mpp_send(grid_global(:,:,:,tile_coarse(n)), &
                  size(grid_global)/Atm%flagstruct%ntiles,grids_master_procs(n))
          else
             call mpp_send(grid_global(:,:,:,1),size(grid_global),grids_master_procs(n))
          endif
          call mpp_sync_self()
       endif
    enddo

    if (Atm%neststruct%nested .or. ANY(Atm%neststruct%child_grids)) then
       nullify(grid_global)
    else if( trim(grid_file) .EQ. 'Inline') then
       deallocate(grid_global)
    endif

    nullify(agrid)
    nullify(grid)

    nullify( area)
    nullify(rarea)
    nullify( area_c)
    nullify(rarea_c)

    if ( Atm%flagstruct%molecular_diffusion ) then
       nullify( area_u)
       nullify( area_v)
       nullify(rarea_u)
       nullify(rarea_v)
       nullify( dx6)
       nullify( dy6)
       nullify(rdx6)
       nullify(rdy6)
    endif

    nullify(sina)
    nullify(cosa)
    nullify(dx)
    nullify(dy)
    nullify(dxc)
    nullify(dyc)
    nullify(dxa)
    nullify(dya)
    nullify(rdx)
    nullify(rdy)
    nullify(rdxc)
    nullify(rdyc)
    nullify(rdxa)
    nullify(rdya)
    nullify(e1)
    nullify(e2)

    nullify(iinta)
    nullify(jinta)
    nullify(iintb)
    nullify(jintb)
    nullify(npx_g)
    nullify(npy_g)
    nullify(ntiles_g)
    nullify(sw_corner)
    nullify(se_corner)
    nullify(ne_corner)
    nullify(nw_corner)
    nullify(latlon)
    nullify(cubed_sphere)
    nullify(have_south_pole)
    nullify(have_north_pole)
    nullify(stretched_grid)

    nullify(tile)

    nullify(domain)

    if (use_timer) call mpp_clock_end (id_timer8)

  contains

    subroutine setup_cartesian(npx, npy, dx_const, dy_const, deglat, bd)

      type(fv_grid_bounds_type), intent(IN) :: bd
       integer, intent(in):: npx, npy
       real(kind=R_GRID), intent(IN) :: dx_const, dy_const, deglat
       real(kind=R_GRID) lat_rad, lon_rad, domain_rad
       integer i,j
       integer :: is,  ie,  js,  je
       integer :: isd, ied, jsd, jed

       is  = bd%is
       ie  = bd%ie
       js  = bd%js
       je  = bd%je
       isd = bd%isd
       ied = bd%ied
       jsd = bd%jsd
       jed = bd%jed

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


    ! Subroutine to be used by setup_aligned_nest to configure the nest grid -- either the entire grid, or just the leading edge
    !  based on the input dimensions in range_x and range_y.  Algorithm copied from setup_aligned_nest.
    subroutine compute_nest_points(p_grid, p_ind, out_grid, refinement, ioffset, joffset, range_x, range_y, isg, ieg, jsg, jeg)
      real(kind=R_GRID), allocatable, intent(in)      :: p_grid(:,:,:)
      !integer, intent(inout)                          :: p_ind(:,:,:)
      integer, intent(inout)                          :: p_ind(1-ng:npx  +ng,1-ng:npy  +ng,4)
      real(kind=R_GRID), allocatable, intent(inout)   :: out_grid(:,:,:,:)
      integer, intent(in)                             :: refinement, ioffset, joffset
      integer, intent(in)                             :: range_x(2), range_y(2)
      integer, intent(in)                             :: isg, ieg, jsg, jeg

      real(kind=R_GRID), dimension(2) :: q1, q2
      integer  :: i, j, ic, jc, imod, jmod
      integer  :: this_pe

      ! Need isg, ieg, jsg, jeg

      this_pe = mpp_pe()

      if (debug_log) print '("[INFO] Filling out_grid(",I0,"-",I0,",",I0,"-",I0,",1-2,1) in compute_nest_points fv_grid_tools.F90. npe=",I0)', range_x(1), range_x(2), range_y(1), range_y(2), this_pe

      do j=range_y(1), range_y(2)
         jc = joffset + (j-1)/refinement !int( real(j-1) / real(refinement) )
         jmod = mod(j-1,refinement)
         if (j-1 < 0 .and. jmod /= 0) jc = jc - 1
         if (jmod < 0) jmod = jmod + refinement

         do i=range_x(1), range_x(2)
            ic = ioffset + (i-1)/refinement !int( real(i-1) / real(refinement) )
            imod = mod(i-1,refinement)
            if (i-1 < 0 .and. imod /= 0) ic = ic - 1
            if (imod < 0) imod = imod + refinement

            if (ic+1 > ieg+1 .or. ic < isg .or. jc+1 > jeg+1 .or. jc < jsg) then
               print*, 'p_grid:',  i, j,  ' OUT OF BOUNDS'
               print*, ic, jc
               print*, isg, ieg, jsg, jeg
               print*, imod, jmod
            end if

            if (jmod == 0) then
               q1 = p_grid(ic, jc, 1:2)
               q2 = p_grid(ic+1,jc,1:2)
            else
               call spherical_linear_interpolation( real(jmod,kind=R_GRID)/real(refinement,kind=R_GRID),  &
                    p_grid(ic, jc, 1:2), p_grid(ic, jc+1, 1:2), q1)
               call spherical_linear_interpolation( real(jmod,kind=R_GRID)/real(refinement,kind=R_GRID),  &
                    p_grid(ic+1, jc, 1:2), p_grid(ic+1, jc+1, 1:2), q2)
            end if

            if (imod == 0) then
               out_grid(i,j,:,1) = q1
            else
               call spherical_linear_interpolation( real(imod,kind=R_GRID)/real(refinement,kind=R_GRID),  &
                    q1,q2,out_grid(i,j,:,1))
            end if

            !SW coarse-grid index; assumes grid does
            !not overlie other cube panels. (These indices
            !are also for the corners and thus need modification
            !to be used for cell-centered and edge-
            !centered variables; see below)
            p_ind(i,j,1) = ic
            p_ind(i,j,2) = jc
            p_ind(i,j,3) = imod
            p_ind(i,j,4) = jmod

            if (out_grid(i,j,1,1) > 2.*pi) out_grid(i,j,1,1) = out_grid(i,j,1,1) - 2.*pi
            if (out_grid(i,j,1,1) < 0.) out_grid(i,j,1,1) = out_grid(i,j,1,1) + 2.*pi

         end do
      end do
    end subroutine compute_nest_points

    subroutine setup_orthogonal_grid(npx, npy, bd, grid_file)
      type(fv_grid_bounds_type), intent(IN) :: bd
      character(len=*),    intent(IN)    :: grid_file
      integer,      intent(IN) :: npx, npy

      ! real(kind=R_GRID), pointer, dimension(:,:,:) :: agrid, grid
      ! real(kind=R_GRID), pointer, dimension(:,:) :: area, area_c
      ! real(kind=R_GRID), pointer, dimension(:,:) :: dx, dy, dxc, dyc, dxa, dya

      ! real, pointer, dimension(:,:) :: rarea, rarea_c
      ! real, pointer, dimension(:,:) :: rdx, rdy, rdxc, rdyc, rdxa, rdya
      ! real, pointer, dimension(:,:,:) :: e1, e2


      type(FmsNetcdfFile_t)              :: Grid_input
      character(len=256)                 :: atm_mosaic, atm_hgrid
      real, allocatable, dimension(:,:)  :: tmpx, tmpy, tmpu, tmpv, tmpa

      integer i, j, stdunit
      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed
      integer :: isc2, iec2, jsc2, jec2
      integer :: start(4), nread(4)
      integer,save :: halo=3

      real(kind=R_GRID)  :: dxN, dxM, dxAV
      real(kind=R_GRID)  :: dx_local, dy_local
      real(kind=R_GRID)  :: maxarea, minarea, globalarea


      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed


      if(.not. file_exists(grid_file)) call mpp_error(FATAL, 'fv_grid_tools(read_grid): file '// &
      trim(grid_file)//' does not exist')

      !--- make sure the grid file is mosaic file.
      if( open_file(Grid_input, grid_file, "read") ) then
        if( variable_exists(Grid_input, 'atm_mosaic_file') .OR. variable_exists(Grid_input, 'gridfiles') ) then
          stdunit = stdout()
          write(stdunit,*) '==>Note from fv_grid_tools_mod(read_grid): read atmosphere grid from mosaic version grid'
        else
          call mpp_error(FATAL, 'fv_grid_tools(read_grid): neither atm_mosaic_file nor gridfiles exists in file ' &
               //trim(grid_file))
        endif

        if(variable_exists(Grid_input, 'atm_mosaic_file') ) then
          call read_data(Grid_input, "atm_mosaic_file", atm_mosaic)
          atm_mosaic = "INPUT/"//trim(atm_mosaic)
        else
          atm_mosaic = trim(grid_file)
        endif
        call close_file(Grid_input)
      endif

      call get_mosaic_tile_grid(atm_hgrid, atm_mosaic, Atm%domain)


      !--- get the geographical coordinates of super-grid.

      isc2 = 2*(isd+halo)-1; iec2 = 2*(ied+1+halo)-1   ! For the regional domain the cell corner locations must be transferred
      jsc2 = 2*(jsd+halo)-1; jec2 = 2*(jed+1+halo)-1   ! from the entire supergrid to the compute grid, including the halo region.


      allocate(tmpx(isc2:iec2, jsc2:jec2) )
      allocate(tmpy(isc2:iec2, jsc2:jec2) )
      start = 1; nread = 1
      start(1) = isc2; nread(1) = iec2 - isc2 + 1
      start(2) = jsc2; nread(2) = jec2 - jsc2 + 1
      if (open_file(Grid_input, atm_hgrid, "read")) then
        call read_data(Grid_input, 'x', tmpx, corner=start, edge_lengths=nread)  !<-- tmpx (lon, deg east) is on the supergrid
        call read_data(Grid_input, 'y', tmpy, corner=start, edge_lengths=nread)  !<-- tmpy (lat, deg) is on the supergrid

        !--- geographic grid at cell corner
        grid(isd: is-1, jsd:js-1,1:ndims)=0.
        grid(isd: is-1, je+2:jed+1,1:ndims)=0.
        grid(ie+2:ied+1,jsd:js-1,1:ndims)=0.
        grid(ie+2:ied+1,je+2:jed+1,1:ndims)=0.


        do j = jsd, jed+1
          do i = isd, ied+1
            grid(i,j,1) = tmpx(2*i+halo+2,2*j+halo+2)*pi/180.
            grid(i,j,2) = tmpy(2*i+halo+2,2*j+halo+2)*pi/180.
          enddo
        enddo

        call mpp_update_domains( grid, Atm%domain, position=CORNER)

        iec2 = 2*(ied+1+halo)-2   ! For the regional domain the cell corner locations must be transferred
        jec2 = 2*(jed+1+halo)-1   ! from the entire supergrid to the compute grid, including the halo region.

        allocate(tmpu(isc2:iec2, jsc2:jec2) )

        nread(1) = iec2 - isc2 + 1
        nread(2) = jec2 - jsc2 + 1
        call read_data(Grid_input, 'dx', tmpu, corner=start, edge_lengths=nread)


        do j = jsd, jed+1
          do i = isd, ied
            dx(i,j) = tmpu(2*i+halo+2,2*j+halo+2) + tmpu(2*i+halo+3,2*j+halo+2)
          enddo
        enddo

        iec2 = 2*(ied+1+halo)-1   ! For the regional domain the cell corner locations must be transferred
        jec2 = 2*(jed+1+halo)-2   ! from the entire supergrid to the compute grid, including the halo region.

        allocate(tmpv(isc2:iec2, jsc2:jec2) )

        nread(1) = iec2 - isc2 + 1
        nread(2) = jec2 - jsc2 + 1
        call read_data(Grid_input, 'dy', tmpv, corner=start, edge_lengths=nread)


        do j = jsd, jed
          do i = isd, ied+1
            dy(i,j) = tmpv(2*i+halo+2,2*j+halo+2) + tmpv(2*i+halo+2,2*j+halo+3)
          enddo
        enddo


        call mpp_update_domains( dy, dx, Atm%domain, flags=SCALAR_PAIR,      &
        gridtype=CGRID_NE_PARAM, complete=.true.)

        iec2 = 2*(ied+1+halo)-2   ! For the regional domain the cell corner locations must be transferred
        jec2 = 2*(jed+1+halo)-2   ! from the entire supergrid to the compute grid, including the halo region.

        allocate(tmpa(isc2:iec2, jsc2:jec2) )

        nread(1) = iec2 - isc2 + 1
        nread(2) = jec2 - jsc2 + 1
        call read_data(Grid_input, 'area', tmpa, corner=start, edge_lengths=nread) !<-- tmpx (lon, deg east) is on the supergrid
        call close_file(Grid_input)
      endif


      !agrid(:,:,:) = -1.e25
      area_c(:,:) = -missing ! To prevent divide by zero error


      do j = jsd, jed
        do i = isd, ied
          agrid(i,j,1) = tmpx(2*i+halo+3,2*j+halo+3)*pi/180.
          agrid(i,j,2) = tmpy(2*i+halo+3,2*j+halo+3)*pi/180.

              dxa(i,j) = tmpu(2*i+halo+2,2*j+halo+3) + tmpu(2*i+halo+3,2*j+halo+3)
              dya(i,j) = tmpv(2*i+halo+3,2*j+halo+2) + tmpv(2*i+halo+3,2*j+halo+3)

             area(i,j) = tmpa(2*i+halo+2,2*j+halo+2) + tmpa(2*i+halo+3,2*j+halo+2) + tmpa(2*i+halo+2,2*j+halo+3) + tmpa(2*i+halo+3,2*j+halo+3)

        enddo
      enddo

      call mpp_update_domains( agrid, Atm%domain, position=CENTER, complete=.true. )
      call mpp_update_domains( area,   Atm%domain, complete=.true. )
      call mpp_update_domains( dxa, dya, Atm%domain, flags=SCALAR_PAIR, gridtype=AGRID_PARAM)

      do j = jsd+1, jed
        do i = isd+1, ied
          area_c(i,j) = tmpa(2*i+halo+2,2*j+halo+2) + tmpa(2*i+halo+1,2*j+halo+2) + tmpa(2*i+halo+2,2*j+halo+1) + tmpa(2*i+halo+1,2*j+halo+1)
        enddo
      enddo

      if (is == 1) then
        do j=jsd,jed
          area_c(isd,j) = area_c(isd+1,j)
        end do
        if (js == 1)     area_c(isd,jsd) = area_c(isd+1,jsd+1)
        if (js == npy-1) area_c(isd,jed+1) = area_c(isd+1,jed)
      end if
      if (ie == npx-1) then
        do j=jsd,jed
          area_c(ied+1,j) = area_c(ied,j)
        end do
        if (js == 1)     area_c(ied+1,jsd) = area_c(ied,jsd+1)
        if (js == npy-1) area_c(ied+1,jed+1) = area_c(ied,jed)
      end if
      if (js == 1) then
        do i=isd,ied
          area_c(i,jsd) = area_c(i,jsd+1)
        end do
      end if
      if (je == npy-1) then
        do i=isd,ied
          area_c(i,jed+1) = area_c(i,jed)
        end do
      end if


      do j=jsd,jed
        do i=isd+1,ied
          dxc(i,j) = tmpu(2*i+halo+1,2*j+halo+3) + tmpu(2*i+halo+2,2*j+halo+3)

        enddo
        !xxxxxx
        !Are the following 2 lines appropriate for the regional domain?
        !xxxxxx
        dxc(isd,j)   = dxc(isd+1,j)
        dxc(ied+1,j) = dxc(ied,j)
      enddo


      do j=jsd+1,jed
        do i=isd,ied
          dyc(i,j) = tmpv(2*i+halo+3,2*j+halo+1) + tmpv(2*i+halo+3,2*j+halo+2)
        enddo
      enddo
      !xxxxxx
      !Are the following 2 lines appropriate for the regional domain?
      !xxxxxx
      do i=isd,ied
        dyc(i,jsd)   = dyc(i,jsd+1)
        dyc(i,jed+1) = dyc(i,jed)
      end do

      call mpp_update_domains( dxc, dyc, Atm%domain, flags=SCALAR_PAIR,   &
      gridtype=CGRID_NE_PARAM, complete=.true.)

      call mpp_update_domains( area_c, Atm%domain, position=CORNER, complete=.true.)


      do j=jsd,jed+1
        do i=isd,ied
          rdx(i,j) = 1.0/dx(i,j)
          rdyc(i,j) = 1.0/dyc(i,j)
        enddo
      enddo
      do j=jsd,jed
        do i=isd,ied+1
          rdy(i,j) = 1.0/dy(i,j)
          rdxc(i,j) = 1.0/dxc(i,j)
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



      ! Get and print Grid Statistics
      dxAV =0.0

      dxN  =  missing
      dxM  = -missing

      do j=js, je
        do i=is, ie
          if(i>ceiling(npx/2.) .OR. j>ceiling(npy/2.)) cycle

          dx_local = dx(i,j)
          dy_local = dy(i,j)

          dxAV  = dxAV + 0.5 * (dx_local + dy_local)
          dxM   = MAX(dxM,dx_local)
          dxM   = MAX(dxM,dy_local)
          dxN   = MIN(dxN,dx_local)
          dxN   = MIN(dxN,dy_local)

        enddo
      enddo


      call mpp_sum(dxAV)
      call mpp_max(dxM)
      call mpp_min(dxN)

      globalarea = mpp_global_sum(domain, area)
      maxarea = mpp_global_max(domain, area)
      minarea = mpp_global_min(domain, area)

      if( is_master() ) then

        dxAV  = dxAV  / ( (ceiling(npy/2.0))*(ceiling(npx/2.0)) )

        write(*,*  ) ''
        write(*,*  ) ' Lambert Grid Stats : ', npx,'x',npy,'x 1'
        write(*,201) '      Grid Length   : min: ', dxN,' max: ', dxM,' avg: ', dxAV, ' min/max: ',dxN/dxM
        write(*,*  ) ''
        write(*,209) '   MAX    AREA (m*m):', maxarea,            '          MIN AREA (m*m):', minarea
        write(*,210) '   GLOBAL AREA (m*m):', globalarea
        write(*,*  ) ''

201  format(A,f11.2,A,f11.2,A,f11.2,A,f11.2)
209  format(A,e21.14,A,e21.14)
210  format(A,e21.14)

      endif

!      sina(:,:) = 1.
!      cosa(:,:) = 0.

      e1(1,:,:) = 1.
      e1(2,:,:) = 0.
      e1(3,:,:) = 0.

      e2(1,:,:) = 0.
      e2(2,:,:) = 1.
      e2(3,:,:) = 0.


      deallocate(tmpx, tmpy, tmpu, tmpv, tmpa)

    end subroutine setup_orthogonal_grid


    !This routine currently does two things:
    ! 1) Create the nested grid on-the-fly from the parent
    ! 2) Compute the weights and indices for the boundary conditions
    ! We should split these into two routines in case we can
    !   read the nest from the input mosaic. Then we only need
    !   to set up the weights.
    ! When creating the nest on-the-fly we need the global parent grid,
    !   as we are doing now. For nests crossing a cube edge
    !   new code is needed.
    ! Creating the indices should be relatvely straightforward procedure
    !   since we will always know ioffset and joffset, which are needed
    !   to initialize the mpp nesting structure
    ! Computing the weights can be simplified by simply retreiving the
    !   BC agrid/grid structures?
    subroutine setup_aligned_nest(Atm)

      type(fv_atmos_type), intent(INOUT), target :: Atm

      integer :: isd_p, ied_p, jsd_p, jed_p
      integer :: isg, ieg, jsg, jeg
      integer :: ic, jc, imod, jmod

      !  Hold these between executions if moving nest
      real(kind=R_GRID), allocatable, dimension(:,:,:), save  :: p_grid_u, p_grid_v, pa_grid, p_grid
      real(kind=R_GRID), allocatable, dimension(:,:,:) :: c_grid_u, c_grid_v
      integer ::    p_ind(1-ng:npx  +ng,1-ng:npy  +ng,4) !< First two entries along dim 3 are
                                                         !! for the corner source indices;
                                                         !! the last two are for the remainders

      integer, allocatable, save  ::   shift_p_ind(:,:,:)

      integer i,j,k, p
      real(kind=R_GRID) sum
      real(kind=R_GRID) :: dist1, dist2, dist3, dist4
      real(kind=R_GRID), dimension(2) :: q1, q2

      integer, pointer :: parent_tile, refinement, ioffset, joffset
      integer, pointer, dimension(:,:,:) :: ind_h, ind_u, ind_v
      real,    pointer, dimension(:,:,:) :: wt_h, wt_u, wt_v

      integer, pointer, dimension(:,:,:) :: ind_b
      real,    pointer, dimension(:,:,:) :: wt_b

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

    !  Setup timing variables

    logical, save       :: first_time = .true.
    integer, save       :: id_timer1, id_timer2, id_timer3a,  id_timer3b,  id_timer3c,  id_timer3d, id_timer4, id_timer5, id_timer6, id_timer7, id_timer8
    integer, save       :: prev_ioffset, prev_joffset    ! not pointers, because we want to save them between runs of this subroutine
    integer, save       :: move_step
    integer             :: delta_i_c, delta_j_c
    integer             :: range_x(2), range_y(2)

    real(kind=R_GRID), allocatable, dimension(:,:,:,:) :: out_grid

    logical             :: moving_nest = .true.  ! TODO set this from the Atm structure

    if (first_time .and. use_timer) then
       id_timer1     = mpp_clock_id ('setup_aligned_nest Step 1',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
       id_timer2     = mpp_clock_id ('setup_aligned_nest Step 2 sph_lin_interp',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
       id_timer3a    = mpp_clock_id ('setup_aligned_nest Step 3a mid_pt_sphere',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
       id_timer3b    = mpp_clock_id ('setup_aligned_nest Step 3b mid_pt_sphere',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
       id_timer3c    = mpp_clock_id ('setup_aligned_nest Step 3c cell_ctr',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
       id_timer3d    = mpp_clock_id ('setup_aligned_nest Step 3d',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
       id_timer4     = mpp_clock_id ('setup_aligned_nest Step 4',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
       id_timer5     = mpp_clock_id ('setup_aligned_nest Step 5',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
       id_timer6     = mpp_clock_id ('setup_aligned_nest Step 6',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
       id_timer7     = mpp_clock_id ('setup_aligned_nest Step 7',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
       id_timer8     = mpp_clock_id ('setup_aligned_nest Step 8',  flags = clock_flag_default, grain=CLOCK_ROUTINE )

       prev_ioffset = Atm%neststruct%ioffset
       prev_joffset = Atm%neststruct%joffset

       !first_time = .false.
    end if

    if (use_timer) call mpp_clock_begin (id_timer1)

      is  = Atm%bd%is
      ie  = Atm%bd%ie
      js  = Atm%bd%js
      je  = Atm%bd%je
      isd = Atm%bd%isd
      ied = Atm%bd%ied
      jsd = Atm%bd%jsd
      jed = Atm%bd%jed



      parent_tile => Atm%neststruct%parent_tile
      refinement => Atm%neststruct%refinement
      ioffset => Atm%neststruct%ioffset
      joffset => Atm%neststruct%joffset

      ind_h => Atm%neststruct%ind_h
      ind_u => Atm%neststruct%ind_u
      ind_v => Atm%neststruct%ind_v

      wt_h => Atm%neststruct%wt_h
      wt_u => Atm%neststruct%wt_u
      wt_v => Atm%neststruct%wt_v

      ind_b => Atm%neststruct%ind_b
      wt_b => Atm%neststruct%wt_b

      ! For moving nest
      if (first_time) then
         delta_i_c = 0
         delta_j_c = 0
         prev_ioffset = ioffset
         prev_joffset = joffset
      else
         delta_i_c = ioffset - prev_ioffset
         delta_j_c = joffset - prev_joffset
      end if

      if (debug_log) print '("[INFO] WDR setup_aligned_nest fv_grid_tools.F90. npe=",I0," delta_i_c=",I0," delta_j_c=",I0," ioffset=",I0," joffset=",I0)', this_pe, delta_i_c, delta_j_c, ioffset, joffset

      call mpp_get_data_domain( Atm%parent_grid%domain, &
           isd_p,  ied_p,  jsd_p,  jed_p  )
      call mpp_get_global_domain( Atm%parent_grid%domain, &
           isg, ieg, jsg, jeg)

      p_ind = -1000000000

      if (first_time) then
         !! Initial allocation of p_grid_u, pgrid_v, pa_grid, and p_grid
         !! Save these parent grids between executions if using moving nest.

         allocate(p_grid_u(isg:ieg  ,jsg:jeg+1,1:2))
         allocate(p_grid_v(isg:ieg+1,jsg:jeg  ,1:2))
         allocate(pa_grid(isg:ieg,jsg:jeg  ,1:2))

         allocate(p_grid( isg-ng:ieg+1+ng, jsg-ng:jeg+1+ng,1:2) )
         p_grid = 1.e25

      end if

      ! Note this will be called during model initialization, then not repeated once moving nest functionality is used
      !  Moving nest will rely on the saved data in p_grid, which does not change (as long as nest remains on same parent tile).
      if (first_time) then

         !Need to RECEIVE parent grid_global;
         !matching mpp_send of grid_global from parent grid is in init_grid()
         if( is_master() ) then

            call mpp_recv(p_grid( isg-ng:ieg+1+ng, jsg-ng:jeg+1+ng,1:2), size(p_grid( isg-ng:ieg+1+ng, jsg-ng:jeg+1+ng,1:2)), &
                 Atm%parent_grid%pelist(1))

         endif

         call mpp_broadcast( p_grid(isg-ng:ieg+ng+1, jsg-ng:jeg+ng+1, :), &
              (ieg-isg+2+2*ng)*(jeg-jsg+2+2*ng)*ndims, mpp_root_pe() )

       end if

       if (use_timer) call mpp_clock_end (id_timer1)
       if (use_timer) call mpp_clock_begin (id_timer2)

    !!  Setup full grid for nest; confusingly called grid_global.  Each nest PE is computing the grid lat/lons for the entire nest
    !!    not just its section.
    !!  INPUTS:  ioffset, joffset, p_grid
    !!  OUTPUTS:  grid_global

      ! Begin calculate shifted version of global_grid

      if (first_time) allocate(shift_p_ind(1-ng:npx  +ng,1-ng:npy  +ng,4))   ! TODO need to deallocate this somewhere

      if (.not. first_time) then

         ! Make copies of grid_global and p_ind to validate that code is correct
         allocate( out_grid( lbound(grid_global,1):ubound(grid_global,1), &
              lbound(grid_global,2):ubound(grid_global,2), &
              lbound(grid_global,3):ubound(grid_global,3), &
              lbound(grid_global,4):ubound(grid_global,4) ) )

         if (debug_log) print '("[INFO] WDR bounds grid_global setup_nest_grid npe=",I0," grid_global(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', this_pe, lbound(grid_global,1), ubound(grid_global,1), &
              lbound(grid_global,2), ubound(grid_global,2), &
              lbound(grid_global,3), ubound(grid_global,3), &
              lbound(grid_global,4), ubound(grid_global,4)

         if (debug_log) print '("[INFO] WDR bounds out_grid setup_nest_grid npe=",I0," out_grid(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', this_pe, lbound(out_grid,1), ubound(out_grid,1), &
              lbound(out_grid,2), ubound(out_grid,2), &
              lbound(out_grid,3), ubound(out_grid,3), &
              lbound(out_grid,4), ubound(out_grid,4)

         out_grid = grid_global

         if ( delta_i_c .ne. 0 ) then
            if (debug_log) print '("[INFO] setup_nest_grid EOSHIFT delta_i_c=",I0," start. npe=",I0)', delta_i_c, this_pe
            out_grid = eoshift(out_grid, refinement * delta_i_c, DIM=1)
         end if

         if (delta_j_c .ne.  0) then
            if (debug_log) print '("[INFO] setup_nest_grid EOSHIFT delta_j_c=",I0," start. npe=",I0)', delta_j_c, this_pe
            out_grid = eoshift(out_grid, refinement * delta_j_c, DIM=2)
         end if

         shift_p_ind(:,:,1) = shift_p_ind(:,:,1) + delta_i_c
         shift_p_ind(:,:,2) = shift_p_ind(:,:,2) + delta_j_c

         !  Compute nest points on any of the halo edges that are empty.  This could be 1 leading edge for N,S,E, or W motion
         !    or two leading edges for NW, NE, SW, or SE motion.
         range_y(1) = 1-ng
         range_y(2) = npy+ng
         if (delta_i_c .lt. 0) then
            range_x(1) = 1-ng
            range_x(2) = 0
            call compute_nest_points(p_grid, shift_p_ind, out_grid, refinement, ioffset, joffset, range_x, range_y, isg, ieg, jsg, jeg)
         elseif  (delta_i_c .gt. 0) then
            range_x(1) = npx
            range_x(2) = npx+ng
            call compute_nest_points(p_grid, shift_p_ind, out_grid, refinement, ioffset, joffset, range_x, range_y, isg, ieg, jsg, jeg)
         end if

         range_x(1) = 1-ng
         range_x(2) = npx+ng
         if (delta_j_c .lt. 0) then
            range_y(1) = 1-ng
            range_y(2) = 0
            call compute_nest_points(p_grid, shift_p_ind, out_grid, refinement, ioffset, joffset, range_x, range_y, isg, ieg, jsg, jeg)
         elseif (delta_j_c .gt. 0) then
            range_y(1) = npy
            range_y(2) = npy+ng
            call compute_nest_points(p_grid, shift_p_ind, out_grid, refinement, ioffset, joffset, range_x, range_y, isg, ieg, jsg, jeg)
         end if

      end if

         ! End calculate shifted version of global_grid
         !  Validate that they match

         if (debug_log) print '("[INFO] Filling grid_global(",I0,"-",I0,",",I0,"-",I0,",1-2,1) in setup_aligned_grid fv_grid_tools.F90. npe=",I0)', 1-ng, npx+ng, 1-ng, npy+ng, this_pe

         if (first_time) then
      ! Generate grid global and parent_grid indices
      ! Grid global only needed in case we create a new child nest on-the-fly?
      !TODO If reading in grid from disk then simply mpp_GATHER grid global from local grid arrays
      !     in fact for nest we should ONLY gather it WHEN NECESSARY.

         do j=1-ng,npy+ng
            jc = joffset + (j-1)/refinement !int( real(j-1) / real(refinement) )
            jmod = mod(j-1,refinement)
            if (j-1 < 0 .and. jmod /= 0) jc = jc - 1
            if (jmod < 0) jmod = jmod + refinement

            do i=1-ng,npx+ng
               ic = ioffset + (i-1)/refinement !int( real(i-1) / real(refinement) )
               imod = mod(i-1,refinement)
               if (i-1 < 0 .and. imod /= 0) ic = ic - 1
               if (imod < 0) imod = imod + refinement

               if (ic+1 > ieg+1 .or. ic < isg .or. jc+1 > jeg+1 .or. jc < jsg) then
                  print*, 'p_grid:',  i, j,  ' OUT OF BOUNDS'
                  print*, ic, jc
                  print*, isg, ieg, jsg, jeg
                  print*, imod, jmod
               end if

               if (jmod == 0) then
                  q1 = p_grid(ic, jc, 1:2)
                  q2 = p_grid(ic+1,jc,1:2)
               else
                  call spherical_linear_interpolation( real(jmod,kind=R_GRID)/real(refinement,kind=R_GRID),  &
                       p_grid(ic, jc, 1:2), p_grid(ic, jc+1, 1:2), q1)
                  call spherical_linear_interpolation( real(jmod,kind=R_GRID)/real(refinement,kind=R_GRID),  &
                       p_grid(ic+1, jc, 1:2), p_grid(ic+1, jc+1, 1:2), q2)
               end if

               if (imod == 0) then
                  grid_global(i,j,:,1) = q1
               else
                  call spherical_linear_interpolation( real(imod,kind=R_GRID)/real(refinement,kind=R_GRID),  &
                       q1,q2,grid_global(i,j,:,1))
               end if

               !SW coarse-grid index; assumes grid does
               !not overlie other cube panels. (These indices
               !are also for the corners and thus need modification
               !to be used for cell-centered and edge-
               !centered variables; see below)
               p_ind(i,j,1) = ic
               p_ind(i,j,2) = jc
               p_ind(i,j,3) = imod
               p_ind(i,j,4) = jmod

               if (grid_global(i,j,1,1) > 2.*pi) grid_global(i,j,1,1) = grid_global(i,j,1,1) - 2.*pi
               if (grid_global(i,j,1,1) < 0.) grid_global(i,j,1,1) = grid_global(i,j,1,1) + 2.*pi

            end do
         end do
         else
            p_ind = shift_p_ind
            grid_global = out_grid
         end if

         if (use_timer) call mpp_clock_end (id_timer2)

         ! Move this elsewhere later.
         if (.not. first_time) deallocate(out_grid)

         if (first_time) shift_p_ind = p_ind

         !if (.not. first_time) then
         if (.false.) then
            !  Do fully recomputed p_ind and grid_global match with shifted grids?
            do i=1-ng,npx+ng
               do j=1-ng,npy+ng
                  do k=1,4
                     if (p_ind(i,j,k) .ne. shift_p_ind(i,j,k)) then
                        print '("[ERROR] WDR setup_nest_grid MISMATCH p_ind(",I0,",",I0,",",I0,")=",I0," shift_p_ind(",I0,",",I0,",",I0,")=",I0," npe=",I0, " move_step=",I0," ")', &
                             i, j, k, p_ind(i,j,k), i, j, k, shift_p_ind(i,j,k), this_pe, move_step
                     end if
                  end do
               end do
            end do

            do i=1-ng,npx+ng
               do j=1-ng,npy+ng
                  if (abs(grid_global(i,j,1,1) -  out_grid(i,j,1,1)) .gt. 0.01) then
                     print '("[ERROR] WDR setup_nest_grid MISMATCH grid_global(",I0,",",I0,",",I0,",1)=",F18.12," out_grid(",I0,",",I0,",",I0,",1)=",F18.12," npe=",I0," move_step=",I0," ")', &
                          i, j, 1, grid_global(i,j,1,1)*180.0/pi, i, j, 1, out_grid(i,j,1,1)*180.0/pi, this_pe, move_step
                  end if
                  if (abs(grid_global(i,j,2,1) -  out_grid(i,j,2,1)) .gt. 0.01) then
                     print '("[ERROR] WDR setup_nest_grid MISMATCH grid_global(",I0,",",I0,",",I0,",1)=",F18.12," out_grid(",I0,",",I0,",",I0,",1)=",F18.12," npe=",I0, " move_step=",I0," ")', &
                          i, j, 2, grid_global(i,j,2,1)*180.0/pi, i, j, 2, out_grid(i,j,2,1)*180.0/pi, this_pe, move_step
                  end if
               end do
            end do

            ! Move this elsewhere later.
            deallocate(out_grid)

         end if

         if (first_time) then
            ! These are the various staggers of the parent grid
            !  They do not vary if the nest moves.  Safe to preserve them between
            !  calls to this routine to save processing time.

            if (use_timer) call mpp_clock_begin (id_timer3a)

            !!  Setup parent staggered grids
            !!  INPUTS:  p_grid
            !!  OUTPUTS:  p_grid_u, p_grid_v, pa_grid

            ! Set up parent grids for interpolation purposes
            do j=jsg,jeg+1
               do i=isg,ieg
                  call mid_pt_sphere(p_grid(i,  j,1:2), p_grid(i+1,  j,1:2), p_grid_u(i,j,:))
                  !call mid_pt_sphere(p_grid(i,  j,1:2), p_grid(i,  j+1,1:2), p_grid_u(i,j,:))
               end do
            end do

            if (use_timer) call mpp_clock_end (id_timer3a)
            if (use_timer) call mpp_clock_begin (id_timer3b)

            do j=jsg,jeg
               do i=isg,ieg+1
                  call mid_pt_sphere(p_grid(i,  j,1:2), p_grid(i,  j+1,1:2), p_grid_v(i,j,:))
                  !call mid_pt_sphere(p_grid(i,  j,1:2), p_grid(i+1,  j,1:2), p_grid_v(i,j,:))
               end do
            end do
            if (use_timer) call mpp_clock_end (id_timer3b)
            if (use_timer) call mpp_clock_begin (id_timer3c)

            do j=jsg,jeg
               do i=isg,ieg
                  call cell_center2(p_grid(i,j,  1:2), p_grid(i+1,j,  1:2),   &
                       p_grid(i,j+1,1:2), p_grid(i+1,j+1,1:2),   &
                       pa_grid(i,j,1:2) )
               end do
            end do

!!$      !TODO: can we just send around ONE grid and re-calculate
!!$      ! staggered grids from that??
!!$      call mpp_broadcast(grid_global(1-ng:npx+ng,  1-ng:npy+ng  ,:,1), &
!!$           ((npx+ng)-(1-ng)+1)*((npy+ng)-(1-ng)+1)*ndims, mpp_root_pe() )
!!$      call mpp_broadcast(      p_ind(1-ng:npx+ng,  1-ng:npy+ng  ,1:4),   &
!!$           ((npx+ng)-(1-ng)+1)*((npy+ng)-(1-ng)+1)*4, mpp_root_pe() )
!!$      call mpp_broadcast(    pa_grid( isg:ieg  , jsg:jeg  , :), &
!!$           ((ieg-isg+1))*(jeg-jsg+1)*ndims, mpp_root_pe())
!!$      call mpp_broadcast(  p_grid_u( isg:ieg  , jsg:jeg+1, :), &
!!$           (ieg-isg+1)*(jeg-jsg+2)*ndims, mpp_root_pe())
!!$      call mpp_broadcast(  p_grid_v( isg:ieg+1, jsg:jeg  , :), &
!!$           (ieg-isg+2)*(jeg-jsg+1)*ndims, mpp_root_pe())

            if (use_timer) call mpp_clock_end (id_timer3c)

         end if

         if (use_timer) call mpp_clock_begin (id_timer3d)

         !!  Setup "grid" -- what is this doing??
         !!  INPUTS:  grid_global
         !!  OUTPUTS:  grid

         if (Atm%flagstruct%grid_type >= 0) then
            do n=1,ndims
            do j=jsd,jed+1
            do i=isd,ied+1
               grid(i,j,n) = grid_global(i,j,n,1)
            enddo
            enddo
            enddo
         endif

    if (use_timer) call mpp_clock_end (id_timer3d)
    if (use_timer) call mpp_clock_begin (id_timer4)

         ind_h = -999999999
         do j=jsd,jed
         do i=isd,ied
            ic = p_ind(i,j,1)
            jc = p_ind(i,j,2)
            imod = p_ind(i,j,3)
            jmod = p_ind(i,j,4)


            if (imod < refinement/2) then
               ind_h(i,j,1) = ic - 1
            else
               ind_h(i,j,1) = ic
            end if

            if (jmod < refinement/2) then
               ind_h(i,j,2) = jc - 1
            else
               ind_h(i,j,2) = jc
            end if
            ind_h(i,j,3) = imod
            ind_h(i,j,4) = jmod

         end do
         end do

         ind_b = -999999999
         do j=jsd,jed+1
         do i=isd,ied+1
            ic = p_ind(i,j,1)
            jc = p_ind(i,j,2)
            imod = p_ind(i,j,3)
            jmod = p_ind(i,j,4)

            ind_b(i,j,1) = ic
            ind_b(i,j,2) = jc

            ind_b(i,j,3) = imod
            ind_b(i,j,4) = jmod
         enddo
         enddo

         ind_u = -99999999
         !New BCs for wind components:
         ! For aligned grid segments (mod(j-1,R) == 0) set
         !     identically equal to the coarse-grid value
         ! Do linear interpolation in the y-dir elsewhere

         do j=jsd,jed+1
         do i=isd,ied
            ic = p_ind(i,j,1)
            jc = p_ind(i,j,2)
            imod = p_ind(i,j,3)

#ifdef NEW_BC
            ind_u(i,j,1) = ic
#else
            if (imod < refinement/2) then
               ind_u(i,j,1) = ic - 1
            else
               ind_u(i,j,1) = ic
            end if
#endif

            ind_u(i,j,2) = jc
            ind_u(i,j,3) = imod
            ind_u(i,j,4) = p_ind(i,j,4)

         end do
         end do

         ind_v = -999999999

         do j=jsd,jed
         do i=isd,ied+1
            ic = p_ind(i,j,1)
            jc = p_ind(i,j,2)
            jmod = p_ind(i,j,4)

            ind_v(i,j,1) = ic

#ifdef NEW_BC
            ind_v(i,j,2) = jc
#else
            if (jmod < refinement/2) then
               ind_v(i,j,2) = jc - 1
            else
               ind_v(i,j,2) = jc
            end if
#endif

            ind_v(i,j,4) = jmod
            ind_v(i,j,3) = p_ind(i,j,3)
         end do
         end do



         agrid(:,:,:) = -1.e25

         do j=jsd,jed
         do i=isd,ied
            call cell_center2(grid(i,j,  1:2), grid(i+1,j,  1:2),   &
                 grid(i,j+1,1:2), grid(i+1,j+1,1:2),   &
                 agrid(i,j,1:2) )
         enddo
         enddo

         call mpp_update_domains( agrid, Atm%domain, position=CENTER, complete=.true. )

    if (use_timer) call mpp_clock_end (id_timer4)
    if (use_timer) call mpp_clock_begin (id_timer5)

      ! Compute dx
      do j=jsd,jed+1
         do i=isd,ied
            dx(i,j) = great_circle_dist(grid_global(i,j,:,1), grid_global(i+1,j,:,1), radius)
         enddo
      enddo

      ! Compute dy
      do j=jsd,jed
         do i=isd,ied+1
            dy(i,j) = great_circle_dist(grid_global(i,j,:,1), grid_global(i,j+1,:,1), radius)
         enddo
      enddo

      !We will use Michael Herzog's algorithm for computing the weights.

      !h points - need center (a-grid) points of parent grid
      do j=jsd,jed
         do i=isd,ied

            ic = ind_h(i,j,1)
            jc = ind_h(i,j,2)

            dist1 = dist2side_latlon(pa_grid(ic,jc,:)    ,pa_grid(ic,jc+1,:),agrid(i,j,:))
            dist2 = dist2side_latlon(pa_grid(ic,jc+1,:)  ,pa_grid(ic+1,jc+1,:),agrid(i,j,:))
            dist3 = dist2side_latlon(pa_grid(ic+1,jc+1,:),pa_grid(ic+1,jc,:),agrid(i,j,:))
            dist4 = dist2side_latlon(pa_grid(ic,jc,:)  ,pa_grid(ic+1,jc,:),agrid(i,j,:))

            wt_h(i,j,1)=dist2*dist3      ! ic,   jc    weight
            wt_h(i,j,2)=dist3*dist4      ! ic,   jc+1  weight
            wt_h(i,j,3)=dist4*dist1      ! ic+1, jc+1  weight
            wt_h(i,j,4)=dist1*dist2      ! ic+1, jc    weight

            sum=wt_h(i,j,1)+wt_h(i,j,2)+wt_h(i,j,3)+wt_h(i,j,4)
            wt_h(i,j,:)=wt_h(i,j,:)/sum

         end do
      end do

      if (.not. moving_nest) deallocate(pa_grid)

    if (use_timer) call mpp_clock_end (id_timer5)
    if (use_timer) call mpp_clock_begin (id_timer6)

      do j=jsd,jed+1
      do i=isd,ied+1

         ic = ind_b(i,j,1)
         jc = ind_b(i,j,2)

         dist1 = dist2side_latlon(p_grid(ic,jc,:),     p_grid(ic,jc+1,:),   grid(i,j,:))
         dist2 = dist2side_latlon(p_grid(ic,jc+1,:),   p_grid(ic+1,jc+1,:), grid(i,j,:))
         dist3 = dist2side_latlon(p_grid(ic+1,jc+1,:), p_grid(ic+1,jc,:),   grid(i,j,:))
         dist4 = dist2side_latlon(p_grid(ic,jc,:),     p_grid(ic+1,jc,:),   grid(i,j,:))

         wt_b(i,j,1)=dist2*dist3      ! ic,   jc    weight
         wt_b(i,j,2)=dist3*dist4      ! ic,   jc+1  weight
         wt_b(i,j,3)=dist4*dist1      ! ic+1, jc+1  weight
         wt_b(i,j,4)=dist1*dist2      ! ic+1, jc    weight

         sum=wt_b(i,j,1)+wt_b(i,j,2)+wt_b(i,j,3)+wt_b(i,j,4)
         wt_b(i,j,:)=wt_b(i,j,:)/sum

      enddo
      enddo

      if (.not. moving_nest) deallocate(p_grid)


      allocate(c_grid_u(isd:ied+1,jsd:jed,2))
      allocate(c_grid_v(isd:ied,jsd:jed+1,2))

      do j=jsd,jed
         do i=isd,ied+1
            call mid_pt_sphere(grid(i,  j,1:2), grid(i,  j+1,1:2), c_grid_u(i,j,:))
         end do
      end do

      do j=jsd,jed+1
         do i=isd,ied
            call mid_pt_sphere(grid(i,j  ,1:2), grid(i+1,j  ,1:2), c_grid_v(i,j,:))
         end do
      end do


      do j=jsd,jed
         do i=isd,ied
            dxa(i,j) = great_circle_dist(c_grid_u(i,j,:), c_grid_u(i+1,j,:), radius)
         end do
      end do

      do j=jsd,jed
         do i=isd,ied
            dya(i,j) = great_circle_dist(c_grid_v(i,j,:), c_grid_v(i,j+1,:), radius)
         end do
      end do

      if (use_timer) call mpp_clock_end (id_timer6)
      if (use_timer) call mpp_clock_begin (id_timer7)

      !Compute interpolation weights. (Recall that the weights are defined with respect to a d-grid)

      !U weights

      do j=jsd,jed+1
         do i=isd,ied

            ic = ind_u(i,j,1)
            jc = ind_u(i,j,2)

            if (ic+1 > ieg .or. ic < isg .or. jc+1 > jeg+1 .or. jc < jsg) then
               print*, 'IND_U ', i, j, ' OUT OF BOUNDS'
               print*, ic, jc
               print*, isg, ieg, jsg, jeg
            end if


#ifdef NEW_BC
            !New vorticity-conserving weights. Note that for the C-grid winds these
            ! become divergence-conserving weights!!
            jmod = p_ind(i,j,4)
            if (jmod == 0) then
               wt_u(i,j,1) = 1. ; wt_u(i,j,2) = 0.
            else
               dist1 = dist2side_latlon(p_grid_u(ic,jc,:), p_grid_u(ic+1,jc,:), c_grid_v(i,j,:))
               dist2 = dist2side_latlon(p_grid_u(ic,jc+1,:), p_grid_u(ic+1,jc+1,:), c_grid_v(i,j,:))
               sum = dist1+dist2
               wt_u(i,j,1) = dist2/sum
               wt_u(i,j,2) = dist1/sum
            endif
            wt_u(i,j,3:4) = 0.
#else
            dist1 = dist2side_latlon(p_grid_u(ic,jc,:)    ,p_grid_u(ic,jc+1,:),  c_grid_v(i,j,:))
            dist2 = dist2side_latlon(p_grid_u(ic,jc+1,:)  ,p_grid_u(ic+1,jc+1,:),c_grid_v(i,j,:))
            dist3 = dist2side_latlon(p_grid_u(ic+1,jc+1,:),p_grid_u(ic+1,jc,:),  c_grid_v(i,j,:))
            !dist4 = dist2side_latlon(p_grid_u(ic+1,jc,:)  ,p_grid_u(ic,jc,:),    c_grid_v(i,j,:))
            dist4 = dist2side_latlon(p_grid_u(ic,jc,:)  ,p_grid_u(ic+1,jc,:),    c_grid_v(i,j,:))

            wt_u(i,j,1)=dist2*dist3      ! ic,   jc    weight
            wt_u(i,j,2)=dist3*dist4      ! ic,   jc+1  weight
            wt_u(i,j,3)=dist4*dist1      ! ic+1, jc+1  weight
            wt_u(i,j,4)=dist1*dist2      ! ic+1, jc    weight

            sum=wt_u(i,j,1)+wt_u(i,j,2)+wt_u(i,j,3)+wt_u(i,j,4)
            wt_u(i,j,:)=wt_u(i,j,:)/sum
#endif

         end do
      end do
      !v weights

      if (use_timer) call mpp_clock_end (id_timer7)
      if (use_timer) call mpp_clock_begin (id_timer8)

      do j=jsd,jed
         do i=isd,ied+1

            ic = ind_v(i,j,1)
            jc = ind_v(i,j,2)

            if (ic+1 > ieg .or. ic < isg .or. jc+1 > jeg+1 .or. jc < jsg) then
               print*, 'IND_V ', i, j, ' OUT OF BOUNDS'
               print*, ic, jc
               print*, isg, ieg, jsg, jeg
            end if

#ifdef NEW_BC
            imod = p_ind(i,j,3)
            if (imod == 0) then
               wt_v(i,j,1) = 1. ; wt_v(i,j,4) = 0.
            else
               dist1 = dist2side_latlon(p_grid_v(ic,jc,:), p_grid_v(ic,jc+1,:), c_grid_u(i,j,:))
               dist2 = dist2side_latlon(p_grid_v(ic+1,jc,:), p_grid_v(ic+1,jc+1,:), c_grid_u(i,j,:))
               sum = dist1+dist2
               wt_v(i,j,1) = dist2/sum
               wt_v(i,j,4) = dist1/sum
            endif
            wt_v(i,j,2) = 0. ; wt_v(i,j,3) = 0.
#else
            dist1 = dist2side_latlon(p_grid_v(ic,jc,:)    ,p_grid_v(ic,jc+1,:),  c_grid_u(i,j,:))
            dist2 = dist2side_latlon(p_grid_v(ic,jc+1,:)  ,p_grid_v(ic+1,jc+1,:),c_grid_u(i,j,:))
            dist3 = dist2side_latlon(p_grid_v(ic+1,jc+1,:),p_grid_v(ic+1,jc,:),  c_grid_u(i,j,:))
            dist4 = dist2side_latlon(p_grid_v(ic,jc,:)  ,p_grid_v(ic+1,jc,:),    c_grid_u(i,j,:))

            wt_v(i,j,1)=dist2*dist3      ! ic,   jc    weight
            wt_v(i,j,2)=dist3*dist4      ! ic,   jc+1  weight
            wt_v(i,j,3)=dist4*dist1      ! ic+1, jc+1  weight
            wt_v(i,j,4)=dist1*dist2      ! ic+1, jc    weight

            sum=wt_v(i,j,1)+wt_v(i,j,2)+wt_v(i,j,3)+wt_v(i,j,4)
            wt_v(i,j,:)=wt_v(i,j,:)/sum
#endif

         end do
      end do

      deallocate(c_grid_u)
      deallocate(c_grid_v)

      if (.not. moving_nest) deallocate(p_grid_u)
      if (.not. moving_nest) deallocate(p_grid_v)

      if (use_timer) call mpp_clock_end (id_timer8)

      if (is_master()) then
         if (Atm%neststruct%nested) then
            !Nesting position information
            !BUG multiply by 180 not 90....
            write(*,*) 'NESTED GRID ', Atm%grid_number
            ic = p_ind(1,1,1) ; jc = p_ind(1,1,1)
            write(*,'(A, 2I5, 4F10.4)') 'SW CORNER: ', ic, jc, grid_global(1,1,:,1)*90./pi
            ic = p_ind(1,npy,1) ; jc = p_ind(1,npy,1)
            write(*,'(A, 2I5, 4F10.4)') 'NW CORNER: ', ic, jc, grid_global(1,npy,:,1)*90./pi
            ic = p_ind(npx,npy,1) ; jc = p_ind(npx,npy,1)
            write(*,'(A, 2I5, 4F10.4)') 'NE CORNER: ', ic, jc, grid_global(npx,npy,:,1)*90./pi
            ic = p_ind(npx,1,1) ; jc = p_ind(npx,1,1)
            write(*,'(A, 2I5, 4F10.4)') 'SE CORNER: ', ic, jc, grid_global(npx,1,:,1)*90./pi
         else
            write(*,*) 'PARENT GRID ', Atm%parent_grid%grid_number, Atm%parent_grid%global_tile
            ic = p_ind(1,1,1) ; jc = p_ind(1,1,1)
            write(*,'(A, 2I5, 4F10.4)') 'SW CORNER: ', ic, jc, Atm%parent_grid%grid_global(ic,jc,:,parent_tile)*90./pi
            ic = p_ind(1,npy,1) ; jc = p_ind(1,npy,1)
            write(*,'(A, 2I5, 4F10.4)') 'NW CORNER: ', ic, jc, Atm%parent_grid%grid_global(ic,jc,:,parent_tile)*90./pi
            ic = p_ind(npx,npy,1) ; jc = p_ind(npx,npy,1)
            write(*,'(A, 2I5, 4F10.4)') 'NE CORNER: ', ic, jc, Atm%parent_grid%grid_global(ic,jc,:,parent_tile)*90./pi
            ic = p_ind(npx,1,1) ; jc = p_ind(npx,1,1)
            write(*,'(A, 2I5, 4F10.4)') 'SE CORNER: ', ic, jc, Atm%parent_grid%grid_global(ic,jc,:,parent_tile)*90./pi
         endif
      end if

      !  Finalize variables in case moving nest calls this again
      first_time = .false.
      move_step = move_step + 1
      prev_ioffset = ioffset
      prev_joffset = joffset

    end subroutine setup_aligned_nest

    subroutine setup_latlon(deglon_start,deglon_stop, deglat_start, deglat_stop, bd )

      type(fv_grid_bounds_type), intent(IN) :: bd
      real(kind=R_GRID), intent(IN) :: deglon_start,deglon_stop, deglat_start, deglat_stop
      real(kind=R_GRID) :: lon_start, lat_start, area_j
      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed


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
      call mpp_update_domains( area_c, Atm%domain, position=CORNER, complete=.true.)

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
      real(kind=R_GRID) , intent(IN)  :: x, y, z
      real(kind=R_GRID) , intent(OUT) :: lon, lat, r

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
 subroutine spherical_to_cartesian(lon, lat, r, x, y, z)
         real(kind=R_GRID) , intent(IN)  :: lon, lat, r
         real(kind=R_GRID) , intent(OUT) :: x, y, z

         x = r * COS(lon) * cos(lat)
         y = r * SIN(lon) * cos(lat)
#ifdef RIGHT_HAND
         z =  r * SIN(lat)
#else
         z = -r * sin(lat)
#endif
 end subroutine spherical_to_cartesian

!>@brief The subroutine 'rot_3d' rotates points on a sphere in xyz coordinates.
!>@details Converts angle from degrees to radians if necessary
      subroutine rot_3d(axis, x1in, y1in, z1in, angle, x2out, y2out, z2out, degrees, convert)
         integer, intent(IN) :: axis         !< axis of rotation 1=x, 2=y, 3=z
         real(kind=R_GRID) , intent(IN)    :: x1in, y1in, z1in
         real(kind=R_GRID) , intent(INOUT) :: angle        !< angle to rotate in radians
         real(kind=R_GRID) , intent(OUT)   :: x2out, y2out, z2out
         integer, intent(IN), optional :: degrees !< if present convert angle
                                                  !! from degrees to radians
         integer, intent(IN), optional :: convert !< if present convert input point
                                                  !! from spherical to cartesian, rotate,
                                                  !! and convert back

         real(kind=R_GRID)  :: c, s
         real(kind=R_GRID)  :: x1,y1,z1, x2,y2,z2

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



!>brief The function 'get_area_tri' gets the surface area of a cell defined as a triangle
!!on the sphere.
!>@details The area is computed as the spherical excess [area units are based on the units of radius]
      real(kind=R_GRID)  function get_area_tri(ndims, p_1, p_2, p_3) &
                        result (myarea)

      integer, intent(IN)    :: ndims          !< 2=lat/lon, 3=xyz
      real(kind=R_GRID) , intent(IN)    :: p_1(ndims) !
      real(kind=R_GRID) , intent(IN)    :: p_2(ndims) !
      real(kind=R_GRID) , intent(IN)    :: p_3(ndims) !

      real(kind=R_GRID)  :: angA, angB, angC

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

!>@brief The subroutine 'grid_area' gets the surface area on a grid in lat/lon or xyz coordinates.
!>@details Determined by 'ndims' argument: 2=lat/lon, 3=xyz)
!! The area is returned in m^2 on a unit sphere.
      subroutine grid_area(nx, ny, ndims, nregions, bounded_domain, gridstruct, domain, bd )
        type(fv_grid_bounds_type), intent(IN) :: bd
        integer, intent(IN) :: nx, ny, ndims, nregions
        logical, intent(IN) :: bounded_domain
        type(fv_grid_type), intent(IN), target :: gridstruct
        type(domain2d), intent(INOUT) :: domain

         real(kind=R_GRID)  :: p_lL(ndims) !< lower Left
         real(kind=R_GRID)  :: p_uL(ndims) !< upper Left
         real(kind=R_GRID)  :: p_lR(ndims) !< lower Right
         real(kind=R_GRID)  :: p_uR(ndims) !< upper Right
         real(kind=R_GRID)  :: a1, d1, d2, mydx, mydy, globalarea

         real(kind=R_GRID)  :: p1(ndims), p2(ndims), p3(ndims), pi1(ndims), pi2(ndims)

         real(kind=R_GRID)  :: maxarea, minarea

         integer :: i,j,n, nreg
         integer :: nh = 0

         real(kind=R_GRID), allocatable :: p_R8(:,:,:)

         real(kind=R_GRID),    pointer, dimension(:,:,:) :: grid, agrid
         integer, pointer, dimension(:,:,:) :: iinta, jinta, iintb, jintb
         real(kind=R_GRID),    pointer, dimension(:,:)   :: area, area_c

         integer :: is,  ie,  js,  je
         integer :: isd, ied, jsd, jed, ng

         is  = bd%is
         ie  = bd%ie
         js  = bd%js
         je  = bd%je
         isd = bd%isd
         ied = bd%ied
         jsd = bd%jsd
         jed = bd%jed
         ng  = bd%ng

         grid  => gridstruct%grid_64
         agrid => gridstruct%agrid_64
         iinta => gridstruct%iinta
         jinta => gridstruct%jinta
         iintb => gridstruct%iintb
         jintb => gridstruct%jintb

         area   => gridstruct%area_64
         area_c => gridstruct%area_c_64

         if (bounded_domain) nh = ng

         maxarea = -1.e25
         minarea =  1.e25

         globalarea = 0.0
         do j=js-nh,je+nh
            do i=is-nh,ie+nh
               do n=1,ndims
               if ( gridstruct%stretched_grid .or. bounded_domain ) then
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
              !maxarea=MAX(area(i,j),maxarea)
              !minarea=MIN(area(i,j),minarea)
              !globalarea = globalarea + area(i,j)
            enddo
         enddo

         globalarea = mpp_global_sum(domain, area)
         maxarea = mpp_global_max(domain, area)
         minarea = mpp_global_min(domain, area)

        if (is_master()) write(*,209) 'MAX    AREA (m*m):', maxarea,            '          MIN AREA (m*m):', minarea
        if (is_master()) write(*,209) 'GLOBAL AREA (m*m):', globalarea, ' IDEAL GLOBAL AREA (m*m):', 4.0*pi*radius**2
 209  format(A,e21.14,A,e21.14)

        if (bounded_domain) then
           nh = ng-1 !cannot get rarea_c on boundary directly
           area_c = 1.e30
        end if

         do j=js-nh,je+nh+1
            do i=is-nh,ie+nh+1
               do n=1,ndims
               if ( gridstruct%stretched_grid .or. bounded_domain ) then
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
         if (gridstruct%cubed_sphere .and. .not. bounded_domain) then
! SW:
            i=1
            j=1
            if ( (is==1) .and. (js==1) ) then
              do n=1,ndims
               if ( gridstruct%stretched_grid ) then
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
               if ( gridstruct%stretched_grid ) then
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
               if ( gridstruct%stretched_grid ) then
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
               if ( gridstruct%stretched_grid ) then
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

!>@brief The function 'get_angle' gets the angle between 3 points on a sphere in lat/lon or
!! xyz coordinates.
!>@details Determined by the 'ndims' argument: 2=lat/lon, 3=xyz
!! The angle is returned in degrees.
      real(kind=R_GRID)  function get_angle(ndims, p1, p2, p3, rad) result (angle)

         integer, intent(IN) :: ndims         !< 2=lat/lon, 3=xyz
         real(kind=R_GRID) , intent(IN)   :: p1(ndims)
         real(kind=R_GRID) , intent(IN)   :: p2(ndims)
         real(kind=R_GRID) , intent(IN)   :: p3(ndims)
         integer, intent(in), optional:: rad

         real(kind=R_GRID)  :: e1(3), e2(3), e3(3)

         if (ndims == 2) then
            call spherical_to_cartesian(p2(1), p2(2), real(1.,kind=R_GRID), e1(1), e1(2), e1(3))
            call spherical_to_cartesian(p1(1), p1(2), real(1.,kind=R_GRID), e2(1), e2(2), e2(3))
            call spherical_to_cartesian(p3(1), p3(2), real(1.,kind=R_GRID), e3(1), e3(2), e3(3))
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

!>@brief The subroutine 'mirror_grid' mirrors the grid across the 0-longitude line
      subroutine mirror_grid(grid_global,ng,npx,npy,ndims,nregions)
         integer, intent(IN)    :: ng,npx,npy,ndims,nregions
         real(kind=R_GRID)   , intent(INOUT) :: grid_global(1-ng:npx  +ng,1-ng:npy  +ng,ndims,1:nregions)
         integer :: i,j,n,n1,n2,nreg
         real(kind=R_GRID) :: x1,y1,z1, x2,y2,z2, ang

         nreg = 1
         do j=1,ceiling(npy/2.)
            do i=1,ceiling(npx/2.)

            x1 = 0.25d0 * (ABS(grid_global(i        ,j        ,1,nreg)) + &
                           ABS(grid_global(npx-(i-1),j        ,1,nreg)) + &
                           ABS(grid_global(i        ,npy-(j-1),1,nreg)) + &
                           ABS(grid_global(npx-(i-1),npy-(j-1),1,nreg)))
            grid_global(i        ,j        ,1,nreg) = SIGN(x1,grid_global(i        ,j        ,1,nreg))
            grid_global(npx-(i-1),j        ,1,nreg) = SIGN(x1,grid_global(npx-(i-1),j        ,1,nreg))
            grid_global(i        ,npy-(j-1),1,nreg) = SIGN(x1,grid_global(i        ,npy-(j-1),1,nreg))
            grid_global(npx-(i-1),npy-(j-1),1,nreg) = SIGN(x1,grid_global(npx-(i-1),npy-(j-1),1,nreg))

            y1 = 0.25d0 * (ABS(grid_global(i        ,j        ,2,nreg)) + &
                           ABS(grid_global(npx-(i-1),j        ,2,nreg)) + &
                           ABS(grid_global(i        ,npy-(j-1),2,nreg)) + &
                           ABS(grid_global(npx-(i-1),npy-(j-1),2,nreg)))
            grid_global(i        ,j        ,2,nreg) = SIGN(y1,grid_global(i        ,j        ,2,nreg))
            grid_global(npx-(i-1),j        ,2,nreg) = SIGN(y1,grid_global(npx-(i-1),j        ,2,nreg))
            grid_global(i        ,npy-(j-1),2,nreg) = SIGN(y1,grid_global(i        ,npy-(j-1),2,nreg))
            grid_global(npx-(i-1),npy-(j-1),2,nreg) = SIGN(y1,grid_global(npx-(i-1),npy-(j-1),2,nreg))

           ! force dateline/greenwich-meridion consitency
            if (mod(npx,2) /= 0) then
              if ( (i==1+(npx-1)/2.0d0) ) then
                 grid_global(i,j        ,1,nreg) = 0.0d0
                 grid_global(i,npy-(j-1),1,nreg) = 0.0d0
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
                  ang = -90.d0
                  call rot_3d( 3, x1, y1, z1, ang, x2, y2, z2, 1, 1)  ! rotate about the z-axis
               elseif (nreg == 3) then
                  ang = -90.d0
                  call rot_3d( 3, x1, y1, z1, ang, x2, y2, z2, 1, 1)  ! rotate about the z-axis
                  ang = 90.d0
                  call rot_3d( 1, x2, y2, z2, ang, x1, y1, z1, 1, 1)  ! rotate about the x-axis
                  x2=x1
                  y2=y1
                  z2=z1

           ! force North Pole and dateline/greenwich-meridion consitency
                  if (mod(npx,2) /= 0) then
                     if ( (i==1+(npx-1)/2.0d0) .and. (i==j) ) then
                        x2 = 0.0d0
                        y2 = pi/2.0d0
                     endif
                     if ( (j==1+(npy-1)/2.0d0) .and. (i < 1+(npx-1)/2.0d0) ) then
                        x2 = 0.0d0
                     endif
                     if ( (j==1+(npy-1)/2.0d0) .and. (i > 1+(npx-1)/2.0d0) ) then
                        x2 = pi
                     endif
                  endif

               elseif (nreg == 4) then
                  ang = -180.d0
                  call rot_3d( 3, x1, y1, z1, ang, x2, y2, z2, 1, 1)  ! rotate about the z-axis
                  ang = 90.d0
                  call rot_3d( 1, x2, y2, z2, ang, x1, y1, z1, 1, 1)  ! rotate about the x-axis
                  x2=x1
                  y2=y1
                  z2=z1

               ! force dateline/greenwich-meridion consitency
                  if (mod(npx,2) /= 0) then
                    if ( (j==1+(npy-1)/2.0d0) ) then
                       x2 = pi
                    endif
                  endif

               elseif (nreg == 5) then
                  ang = 90.d0
                  call rot_3d( 3, x1, y1, z1, ang, x2, y2, z2, 1, 1)  ! rotate about the z-axis
                  ang = 90.d0
                  call rot_3d( 2, x2, y2, z2, ang, x1, y1, z1, 1, 1)  ! rotate about the y-axis
                  x2=x1
                  y2=y1
                  z2=z1
               elseif (nreg == 6) then
                  ang = 90.d0
                  call rot_3d( 2, x1, y1, z1, ang, x2, y2, z2, 1, 1)  ! rotate about the y-axis
                  ang = 0.d0
                  call rot_3d( 3, x2, y2, z2, ang, x1, y1, z1, 1, 1)  ! rotate about the z-axis
                  x2=x1
                  y2=y1
                  z2=z1

           ! force South Pole and dateline/greenwich-meridion consitency
                  if (mod(npx,2) /= 0) then
                     if ( (i==1+(npx-1)/2.0d0) .and. (i==j) ) then
                        x2 = 0.0d0
                        y2 = -pi/2.0d0
                     endif
                     if ( (i==1+(npx-1)/2.0d0) .and. (j > 1+(npy-1)/2.0d0) ) then
                        x2 = 0.0d0
                     endif
                     if ( (i==1+(npx-1)/2.0d0) .and. (j < 1+(npy-1)/2.0d0) ) then
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

      end module fv_grid_tools_mod

