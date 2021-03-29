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

!>@brief 'The module 'tread_da_increment' contains routines for treating the increments
!! of the prognostic variables that are calculated by the DA process
!>@details This module includes functions to read in the externally calculated increments
!! and applies the increments to the restart variables. Specifically, if the increments are
!! zero, FV3 should reproduce directly from the restart files.
!>@note Please treat the following subroutines as API interfaces, and consult the FV3 team
!! code modification proposal.
!>@warning Expanding the list of increments without the proper knowledge of the FV3 dynamical
!! core is EXTREMELY RISKY, especially for the non-hydrostatic scenario. Such a modification
!! could break the scientific consistency, possibly leading to a simulation crash.
!> @author Xi.Chen <xi.chen@noaa.gov>
!> @date 02/12/2016
!
!  REVISION HISTORY:
!  02/12/2016 - Initial Version
!-------------------------------------------------------------------------------

#ifdef OVERLOAD_R4
#define _GET_VAR1 get_var1_real
#else
#define _GET_VAR1 get_var1_double
#endif

module fv_treat_da_inc_mod

! <table>
! <tr>
!     <th>Module Name</th>
!     <th>Functions Included</th>
!   </tr>
!   <tr>
!     <td>constants_mod</td>
!     <td>pi=>pi_8, omega, grav, kappa, rdgas, rvgas, cp_air</td>
!   </tr>
!   <tr>
!     <td>field_manager_mod</td>
!     <td>MODEL_ATMOS</td>
!   </tr>
!   <tr>
!     <td>fms_mod</td>
!     <td>file_exist, open_namelist_file,close_file, error_mesg, FATAL,
!         check_nml_error, stdlog,write_version_number,set_domain,
!         mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_SUBCOMPONENT,
!         clock_flag_default, nullify_domain</td>
!   </tr>
!   <tr>
!     <td>fv_arrays_mod</td>
!     <td>fv_atmos_type, fv_grid_type, fv_grid_bounds_type,
!          R_GRID</td>
!   </tr>
!   <tr>
!     <td>fv_control_mod</td>
!     <td>fv_init, fv_end, ngrids</td>
!   </tr>
!   <tr>
!     <td>fv_grid_utils_mod</td>
!     <td>ptop_min, g_sum, mid_pt_sphere, get_unit_vect2,
!         get_latlon_vector, inner_prod, cubed_to_latlon</td>
!   </tr>
!   <tr>
!     <td>fms_mod</td>
!     <td>file_exist, read_data, field_exist, write_version_number</td>
!   </tr>
!   <tr>
!     <td>fv_mp_mod</td>
!     <td>ng,is_master,fill_corners,YDir,mp_reduce_min, mp_reduce_max</td>
!   </tr>
!   <tr>
!     <td>mpp_mod</td>
!     <td>mpp_error, FATAL, NOTE, mpp_pe</td>
!   </tr>
!   <tr>
!     <td>mpp_domains_mod</td>
!     <td>mpp_get_tile_id, domain2d, mpp_update_domains,NORTH, EAST</td>
!   </tr>
!   <tr>
!     <td>sim_nc_mod</td>
!     <td>open_ncfile, close_ncfile, get_ncdim1, get_var1_double,
!         get_var2_real, get_var3_r4, get_var1_real</td>
!   </tr>
!   <tr>
!     <td>tracer_manager_mod</td>
!     <td>get_tracer_names, get_number_tracers, get_tracer_index</td>
!   </tr>
! </table>

  use fms_mod,           only: file_exist, read_data, &
                               field_exist, write_version_number
  use mpp_mod,           only: mpp_error, FATAL, NOTE, mpp_pe
  use mpp_domains_mod,   only: mpp_get_tile_id, &
                               domain2d, &
                               mpp_update_domains, &
                               NORTH, &
                               EAST
  use tracer_manager_mod,only: get_tracer_names, &
                               get_number_tracers, &
                               get_tracer_index
  use field_manager_mod, only: MODEL_ATMOS

  use constants_mod,     only: pi=>pi_8, omega, grav, kappa, &
                               rdgas, rvgas, cp_air
  use fv_arrays_mod,     only: fv_atmos_type, &
                               fv_grid_type, &
                               fv_grid_bounds_type, &
                               R_GRID
  use fv_grid_utils_mod, only: ptop_min, g_sum, &
                               mid_pt_sphere, get_unit_vect2, &
                               get_latlon_vector, inner_prod, &
                               cubed_to_latlon
  use fv_mp_mod,         only: is_master, &
                               fill_corners, &
                               YDir, &
                               mp_reduce_min, &
                               mp_reduce_max
  use sim_nc_mod,        only: open_ncfile, &
                               close_ncfile, &
                               get_ncdim1, &
                               get_var1_double, &
                               get_var2_real,   &
                               get_var3_r4, &
                               get_var1_real, &
                               check_var_exists
  implicit none
  private

  public :: read_da_inc,remap_coef

contains
  !=============================================================================
  !>@brief The subroutine 'read_da_inc' reads the increments of the diagnostic variables
  !! from the DA-generated files.
  !>@details Additional support of prognostic variables such as tracers can be assessed
  !! and added upon request.
  !>@author Xi.Chen <xi.chen@noaa.gov>
  !>@date 02/12/2016
   subroutine read_da_inc(Atm, fv_domain, bd, npz_in, nq, u, v, q, delp, pt, delz, is_in, js_in, ie_in, je_in,  &
                                                                                   isc_in, jsc_in, iec_in, jec_in )
    type(fv_atmos_type),       intent(inout) :: Atm
    type(domain2d),            intent(inout) :: fv_domain
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer,                   intent(IN) :: npz_in, nq, is_in, js_in, ie_in, je_in, isc_in, jsc_in, iec_in, jec_in
    real, intent(inout), dimension(is_in:ie_in,  js_in:je_in+1,npz_in):: u  ! D grid zonal wind (m/s)
    real, intent(inout), dimension(is_in:ie_in+1,js_in:je_in  ,npz_in):: v  ! D grid meridional wind (m/s)
    real, intent(inout) :: delp(is_in:ie_in  ,js_in:je_in  ,npz_in)  ! pressure thickness (pascal)
    real, intent(inout) :: pt(  is_in:ie_in  ,js_in:je_in  ,npz_in)  ! temperature (K)
    real, intent(inout) :: q(   is_in:ie_in  ,js_in:je_in  ,npz_in, nq)  !
    real, intent(inout) :: delz(isc_in:iec_in  ,jsc_in:jec_in  ,npz_in)  !

    ! local
    real :: deg2rad
    character(len=128) :: fname
    real(kind=4), allocatable:: wk1(:), wk2(:,:), wk3(:,:,:)
    real(kind=4), allocatable:: wk3_u(:,:,:), wk3_v(:,:,:)
    real, allocatable:: tp(:,:,:), qp(:,:,:)
    real, dimension(:,:,:), allocatable:: u_inc, v_inc, ud_inc, vd_inc
    real, allocatable:: lat(:), lon(:)
    real, allocatable:: pt_c(:,:,:), pt_d(:,:,:)
    real:: s2c(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je,4)
    real:: s2c_c(Atm%bd%is:Atm%bd%ie+1,Atm%bd%js:Atm%bd%je,4)
    real:: s2c_d(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je+1,4)
    integer, dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je):: &
        id1, id2, jdc
    integer, dimension(Atm%bd%is:Atm%bd%ie+1,Atm%bd%js:Atm%bd%je)::&
        id1_c, id2_c, jdc_c
    integer, dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je+1)::&
        id1_d, id2_d, jdc_d

    integer:: i, j, k, im, jm, km, npt
    integer:: i1, i2, j1, ncid
    integer:: jbeg, jend
    integer tsize(3)
    real(kind=R_GRID), dimension(2):: p1, p2, p3
    real(kind=R_GRID), dimension(3):: e1, e2, ex, ey

    logical:: found
    integer :: is,  ie,  js,  je
    integer :: isd, ied, jsd, jed
    integer :: sphum, liq_wat
#ifdef MULTI_GASES
    integer :: spfo, spfo2, spfo3
#else
    integer :: o3mr
#endif

    is  = Atm%bd%is
    ie  = Atm%bd%ie
    js  = Atm%bd%js
    je  = Atm%bd%je
    isd = Atm%bd%isd
    ied = Atm%bd%ied
    jsd = Atm%bd%jsd
    jed = Atm%bd%jed

    deg2rad = pi/180.

    fname = 'INPUT/'//Atm%flagstruct%res_latlon_dynamics

    if( file_exist(fname) ) then
      call open_ncfile( fname, ncid )        ! open the file
      call get_ncdim1( ncid, 'lon',   tsize(1) )
      call get_ncdim1( ncid, 'lat',   tsize(2) )
      call get_ncdim1( ncid, 'lev', tsize(3) )

      im = tsize(1); jm = tsize(2); km = tsize(3)

      if (km.ne.npz_in) then
        if (is_master()) print *, 'km = ', km
        call mpp_error(FATAL, &
            '==> Error in read_da_inc: km is not equal to npz_in')
      endif

      if(is_master())  write(*,*) fname, ' DA increment dimensions:', tsize

      allocate (  lon(im) )
      allocate (  lat(jm) )

      call _GET_VAR1 (ncid, 'lon', im, lon )
      call _GET_VAR1 (ncid, 'lat', jm, lat )

      ! Convert to radian
      do i=1,im
        lon(i) = lon(i) * deg2rad  ! lon(1) = 0.
      enddo
      do j=1,jm
        lat(j) = lat(j) * deg2rad
      enddo

    else
      call mpp_error(FATAL,'==> Error in read_da_inc: Expected file '&
          //trim(fname)//' for DA increment does not exist')
    endif

    ! Initialize lat-lon to Cubed bi-linear interpolation coeff:
    call remap_coef( is, ie, js, je, isd, ied, jsd, jed, &
        im, jm, lon, lat, id1, id2, jdc, s2c, &
        Atm%gridstruct%agrid)

    ! Find bounding latitudes:
    jbeg = jm-1;         jend = 2
    do j=js,je
      do i=is,ie
          j1 = jdc(i,j)
        jbeg = min(jbeg, j1)
        jend = max(jend, j1+1)
      enddo
    enddo

    sphum   = get_tracer_index(MODEL_ATMOS, 'sphum')
#ifdef MULTI_GASES
    spfo3   = get_tracer_index(MODEL_ATMOS, 'spfo3')
    spfo    = get_tracer_index(MODEL_ATMOS, 'spfo')
    spfo2   = get_tracer_index(MODEL_ATMOS, 'spfo2')
#else
    o3mr    = get_tracer_index(MODEL_ATMOS, 'o3mr')
#endif
    liq_wat = get_tracer_index(MODEL_ATMOS, 'liq_wat')

    ! perform increments on scalars
    allocate ( wk3(1:im,jbeg:jend, 1:km) )
    allocate (  tp(is:ie,js:je,km) )

    call apply_inc_on_3d_scalar('T_inc',pt, is_in, js_in, ie_in, je_in)
    call apply_inc_on_3d_scalar('delp_inc',delp, is_in, js_in, ie_in, je_in)
    if ( .not. Atm%flagstruct%hydrostatic ) then
        call apply_inc_on_3d_scalar('delz_inc',delz, isc_in, jsc_in, iec_in, jec_in)
    endif
    call apply_inc_on_3d_scalar('sphum_inc',q(:,:,:,sphum), is_in, js_in, ie_in, je_in)
    call apply_inc_on_3d_scalar('liq_wat_inc',q(:,:,:,liq_wat), is_in, js_in, ie_in, je_in)
#ifdef MULTI_GASES
    call apply_inc_on_3d_scalar('spfo3_inc',q(:,:,:,spfo3), is_in, js_in, ie_in, je_in)
    call apply_inc_on_3d_scalar('spfo_inc',q(:,:,:,spfo), is_in, js_in, ie_in, je_in)
    call apply_inc_on_3d_scalar('spfo2_inc',q(:,:,:,spfo2), is_in, js_in, ie_in, je_in)
#else
    call apply_inc_on_3d_scalar('o3mr_inc',q(:,:,:,o3mr), is_in, js_in, ie_in, je_in)
#endif

    deallocate ( tp )
    deallocate ( wk3 )

    ! perform increments on winds
    allocate (pt_c(isd:ied+1,jsd:jed  ,2))
    allocate (pt_d(isd:ied  ,jsd:jed+1,2))
    allocate (ud_inc(is:ie  , js:je+1, km))
    allocate (vd_inc(is:ie+1, js:je  , km))

    call get_staggered_grid( &
        is, ie, js, je, &
        isd, ied, jsd, jed, &
        Atm%gridstruct%grid, pt_c, pt_d)

    !------ pt_c part ------
    ! Initialize lat-lon to Cubed bi-linear interpolation coeff:
    call remap_coef( is, ie+1, js, je, isd, ied+1, jsd, jed, &
        im, jm, lon, lat, id1_c, id2_c, jdc_c, s2c_c, &
        pt_c)

    ! Find bounding latitudes:
    jbeg = jm-1;         jend = 2
    do j=js,je
      do i=is,ie+1
          j1 = jdc_c(i,j)
        jbeg = min(jbeg, j1)
        jend = max(jend, j1+1)
      enddo
    enddo

    allocate ( wk3_u(1:im,jbeg:jend, 1:km) )
    allocate ( wk3_v(1:im,jbeg:jend, 1:km) )
    allocate (  u_inc(is:ie+1,js:je,km) )
    allocate (  v_inc(is:ie+1,js:je,km) )

    call get_var3_r4( ncid, 'u_inc', 1,im, jbeg,jend, 1,km, wk3_u )
    call get_var3_r4( ncid, 'v_inc', 1,im, jbeg,jend, 1,km, wk3_v )

    do k=1,km
      do j=js,je
        do i=is,ie+1
          i1 = id1_c(i,j)
          i2 = id2_c(i,j)
          j1 = jdc_c(i,j)
          u_inc(i,j,k) = s2c_c(i,j,1)*wk3_u(i1,j1  ,k) + &
                         s2c_c(i,j,2)*wk3_u(i2,j1  ,k) + &
                         s2c_c(i,j,3)*wk3_u(i2,j1+1,k) + &
                         s2c_c(i,j,4)*wk3_u(i1,j1+1,k)
          v_inc(i,j,k) = s2c_c(i,j,1)*wk3_v(i1,j1  ,k) + &
                         s2c_c(i,j,2)*wk3_v(i2,j1  ,k) + &
                         s2c_c(i,j,3)*wk3_v(i2,j1+1,k) + &
                         s2c_c(i,j,4)*wk3_v(i1,j1+1,k)
          p1(:) = Atm%gridstruct%grid(i,j  ,1:2)
          p2(:) = Atm%gridstruct%grid(i,j+1,1:2)
          call  mid_pt_sphere(p1, p2, p3)
          call get_unit_vect2(p1, p2, e2)
          call get_latlon_vector(p3, ex, ey)
          vd_inc(i,j,k) = u_inc(i,j,k)*inner_prod(e2,ex) + &
                          v_inc(i,j,k)*inner_prod(e2,ey)
          v(i,j,k) = v(i,j,k) + vd_inc(i,j,k)
        enddo
      enddo
    enddo

    deallocate ( u_inc, v_inc )
    deallocate ( wk3_u, wk3_v )

    !------ pt_d part ------
    ! Initialize lat-lon to Cubed bi-linear interpolation coeff:
    call remap_coef( is, ie, js, je+1, isd, ied, jsd, jed+1, &
        im, jm, lon, lat, id1_d, id2_d, jdc_d, s2c_d, &
        pt_d)

    ! Find bounding latitudes:
    jbeg = jm-1;         jend = 2
    do j=js,je+1
      do i=is,ie
          j1 = jdc_d(i,j)
        jbeg = min(jbeg, j1)
        jend = max(jend, j1+1)
      enddo
    enddo

    allocate ( wk3_u(1:im,jbeg:jend, 1:km) )
    allocate ( wk3_v(1:im,jbeg:jend, 1:km) )
    allocate (  u_inc(is:ie,js:je+1,km) )
    allocate (  v_inc(is:ie,js:je+1,km) )

    call get_var3_r4( ncid, 'u_inc', 1,im, jbeg,jend, 1,km, wk3_u )
    call get_var3_r4( ncid, 'v_inc', 1,im, jbeg,jend, 1,km, wk3_v )

    do k=1,km
      do j=js,je+1
        do i=is,ie
          i1 = id1_d(i,j)
          i2 = id2_d(i,j)
          j1 = jdc_d(i,j)
          u_inc(i,j,k) = s2c_d(i,j,1)*wk3_u(i1,j1  ,k) + &
                         s2c_d(i,j,2)*wk3_u(i2,j1  ,k) + &
                         s2c_d(i,j,3)*wk3_u(i2,j1+1,k) + &
                         s2c_d(i,j,4)*wk3_u(i1,j1+1,k)
          v_inc(i,j,k) = s2c_d(i,j,1)*wk3_v(i1,j1  ,k) + &
                         s2c_d(i,j,2)*wk3_v(i2,j1  ,k) + &
                         s2c_d(i,j,3)*wk3_v(i2,j1+1,k) + &
                         s2c_d(i,j,4)*wk3_v(i1,j1+1,k)
          p1(:) = Atm%gridstruct%grid(i,  j,1:2)
          p2(:) = Atm%gridstruct%grid(i+1,j,1:2)
          call  mid_pt_sphere(p1, p2, p3)
          call get_unit_vect2(p1, p2, e1)
          call get_latlon_vector(p3, ex, ey)
          ud_inc(i,j,k) = u_inc(i,j,k)*inner_prod(e1,ex) + &
                          v_inc(i,j,k)*inner_prod(e1,ey)
          u(i,j,k) = u(i,j,k) + ud_inc(i,j,k)
        enddo
      enddo
    enddo

    deallocate ( u_inc, v_inc )
    deallocate ( wk3_u, wk3_v )

!rab The following is not necessary as ua/va will be re-calculated during model startup
!rab    call cubed_to_latlon(Atm%u, Atm%v, Atm%ua, Atm%va, &
!rab        Atm%gridstruct, Atm%flagstruct%npx, Atm%flagstruct%npy, &
!rab        Atm%flagstruct%npz, 1, Atm%gridstruct%grid_type, &
!rab        fv_domain, Atm%gridstruct%nested, &
!rab        Atm%flagstruct%c2l_ord, Atm%bd)

    !------ winds clean up ------
    deallocate ( pt_c, pt_d, ud_inc, vd_inc )
    !------ all clean up ------
    deallocate ( lon, lat )

  contains
   !---------------------------------------------------------------------------
   !> @brief The subroutine 'apply_inc_on3d_scalar' applies the input increments
   !! to the prognostic variables.
    subroutine apply_inc_on_3d_scalar(field_name,var, is_in, js_in, ie_in, je_in)
      character(len=*), intent(in) :: field_name
      integer, intent(IN) :: is_in, js_in, ie_in, je_in
      real, dimension(is_in:ie_in,js_in:je_in,1:km), intent(inout) :: var
      integer :: ierr

      call check_var_exists(ncid, field_name, ierr)
      if (ierr == 0) then
         call get_var3_r4( ncid, field_name, 1,im, jbeg,jend, 1,km, wk3 )
      else
         if (is_master()) print *,'warning: no increment for ',trim(field_name),' found, assuming zero'
         wk3 = 0.
      endif

      do k=1,km
        do j=js,je
          do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            tp(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k)+&
                        s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
            var(i,j,k) = var(i,j,k)+tp(i,j,k)
          enddo
        enddo
      enddo

    end subroutine apply_inc_on_3d_scalar
    !---------------------------------------------------------------------------
  end subroutine read_da_inc

  !=============================================================================
  !>@brief The subroutine 'remap_coef' calculates the coefficients for horizonal regridding.

  subroutine remap_coef( is, ie, js, je, isd, ied, jsd, jed, &
      im, jm, lon, lat, id1, id2, jdc, s2c, agrid )

    integer, intent(in):: is, ie, js, je, isd, ied, jsd, jed
    integer, intent(in):: im, jm
    real,    intent(in):: lon(im), lat(jm)
    real,    intent(out):: s2c(is:ie,js:je,4)
    integer, intent(out), dimension(is:ie,js:je):: id1, id2, jdc
    real,    intent(in):: agrid(isd:ied,jsd:jed,2)
    ! local:
    real :: rdlon(im)
    real :: rdlat(jm)
    real:: a1, b1
    integer i,j, i1, i2, jc, i0, j0
    do i=1,im-1
      rdlon(i) = 1. / (lon(i+1) - lon(i))
    enddo
    rdlon(im) = 1. / (lon(1) + 2.*pi - lon(im))

    do j=1,jm-1
      rdlat(j) = 1. / (lat(j+1) - lat(j))
    enddo

    ! * Interpolate to cubed sphere cell center
    do 5000 j=js,je

      do i=is,ie

        if ( agrid(i,j,1)>lon(im) ) then
          i1 = im;     i2 = 1
          a1 = (agrid(i,j,1)-lon(im)) * rdlon(im)
        elseif ( agrid(i,j,1)<lon(1) ) then
          i1 = im;     i2 = 1
          a1 = (agrid(i,j,1)+2.*pi-lon(im)) * rdlon(im)
        else
          do i0=1,im-1
            if ( agrid(i,j,1)>=lon(i0) .and. agrid(i,j,1)<=lon(i0+1) ) then
              i1 = i0;  i2 = i0+1
              a1 = (agrid(i,j,1)-lon(i1)) * rdlon(i0)
              go to 111
            endif
          enddo
        endif
111     continue

        if ( agrid(i,j,2)<lat(1) ) then
          jc = 1
          b1 = 0.
        elseif ( agrid(i,j,2)>lat(jm) ) then
          jc = jm-1
          b1 = 1.
        else
          do j0=1,jm-1
            if ( agrid(i,j,2)>=lat(j0) .and. agrid(i,j,2)<=lat(j0+1) ) then
              jc = j0
              b1 = (agrid(i,j,2)-lat(jc)) * rdlat(jc)
              go to 222
            endif
          enddo
        endif
222     continue

        if ( a1<0.0 .or. a1>1.0 .or.  b1<0.0 .or. b1>1.0 ) then
             write(*,*) 'gid=', mpp_pe(), i,j,a1, b1
        endif

        s2c(i,j,1) = (1.-a1) * (1.-b1)
        s2c(i,j,2) =     a1  * (1.-b1)
        s2c(i,j,3) =     a1  *     b1
        s2c(i,j,4) = (1.-a1) *     b1
        id1(i,j) = i1
        id2(i,j) = i2
        jdc(i,j) = jc
      enddo   !i-loop
5000 continue   ! j-loop

  end subroutine remap_coef

  !=============================================================================
  !>@brief The subroutine 'get_staggered_grid' gets the lat-lon locations of the
  !! staggered points.
  subroutine get_staggered_grid( &
      is, ie, js, je, &
      isd, ied, jsd, jed, &
      pt_b, pt_c, pt_d)
    integer, intent(in) :: is, ie, js, je, isd, ied, jsd, jed
    real, dimension(isd:ied+1,jsd:jed+1,2), intent(in) :: pt_b
    real, dimension(isd:ied+1,jsd:jed  ,2), intent(out) :: pt_c
    real, dimension(isd:ied  ,jsd:jed+1,2), intent(out) :: pt_d
    ! local
    real(kind=R_GRID), dimension(2):: p1, p2, p3
    integer :: i, j

    do j = js,je+1
      do i = is,ie
        p1(:) = pt_b(i,  j,1:2)
        p2(:) = pt_b(i+1,j,1:2)
        call  mid_pt_sphere(p1, p2, p3)
        pt_d(i,j,1:2) = p3(:)
      enddo
    enddo
    do j = js,je
      do i = is,ie+1
        p1(:) = pt_b(i,j  ,1:2)
        p2(:) = pt_b(i,j+1,1:2)
        call  mid_pt_sphere(p1, p2, p3)
        pt_c(i,j,1:2) = p3(:)
      enddo
    enddo

  end subroutine get_staggered_grid
  !=============================================================================
end module fv_treat_da_inc_mod

