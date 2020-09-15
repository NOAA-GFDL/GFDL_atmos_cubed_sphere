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

module fv_regional_mod

   use netcdf
   use mpp_domains_mod,   only: domain2d
   use mpp_domains_mod,   only: domain1D, mpp_get_domain_components,    &
                                mpp_get_global_domain,                  &
                                mpp_get_data_domain,                    &
                                mpp_get_compute_domain,                 &
                                NORTH, SOUTH, EAST, WEST,               &
                                CENTER, CORNER,                         &
                                mpp_domains_set_stack_size,             &
                                mpp_update_domains, mpp_get_neighbor_pe
   use mpp_mod,           only: FATAL, input_nml_file,                  &
                                mpp_error ,mpp_pe, mpp_sync,            &
                                mpp_npes, mpp_root_pe, mpp_gather,      &
                                mpp_get_current_pelist, NOTE, NULL_PE
   use mpp_io_mod
   use tracer_manager_mod,only: get_tracer_index,get_tracer_names
   use field_manager_mod, only: MODEL_ATMOS
   use time_manager_mod,  only: get_time                                &
                               ,operator(-),operator(/)                 &
                               ,time_type,time_type_to_real 
   use constants_mod,     only: cp_air, cp_vapor, grav, kappa           &
                               ,pi=>pi_8,rdgas, rvgas
   use fv_arrays_mod,     only: fv_atmos_type                           &
                               ,fv_grid_bounds_type                     &
                               ,fv_regional_bc_bounds_type              &
                               ,R_GRID                                  &
                               ,fv_nest_BC_type_3D                      &
                               ,allocate_fv_nest_BC_type

   use fv_diagnostics_mod,only: prt_gb_nh_sh, prt_height
   use fv_grid_utils_mod, only: g_sum,mid_pt_sphere,get_unit_vect2      &
                               ,get_latlon_vector,inner_prod            &
                               ,cell_center2
   use fv_mapz_mod,       only: mappm, moist_cp, moist_cv
   use fv_mp_mod,         only: is_master, mp_reduce_min, mp_reduce_max
   use fv_fill_mod,       only: fillz
   use fv_eta_mod,        only: get_eta_level
   use fms_mod,           only: check_nml_error,file_exist
   use fms_io_mod,        only: read_data,get_global_att_value
   use boundary_mod,      only: fv_nest_BC_type_3D

   implicit none 

      private

      public ak_in, bk_in                                               &
            ,bc_hour                                                    &
            ,bc_time_interval                                           &
            ,BC_t0,BC_t1                                                &
            ,begin_regional_restart,exch_uv                             &
            ,ntimesteps_per_bc_update                                   &
            ,read_new_bc_data                                           &
            ,regional_bc_data                                           &
            ,regional_bc_t1_to_t0                                       &
            ,regional_boundary_update                                   &
            ,next_time_to_read_bcs                                      &
            ,set_regional_BCs                                            &
            ,setup_regional_BC                                          &
            ,start_regional_cold_start                                  &
            ,start_regional_restart                                     &
            ,dump_field                                                 &
            ,current_time_in_seconds                                    &
            ,a_step, p_step, k_step, n_step, get_data_source            &
            ,write_full_fields
      integer,parameter :: bc_time_interval=3                           &
                          ,nhalo_data =4                                &
                          ,nhalo_model=3
!
      integer, public, parameter :: H_STAGGER = 1
      integer, public, parameter :: U_STAGGER = 2
      integer, public, parameter :: V_STAGGER = 3

      !These parameters are ONLY used for the dump_field debugging routines
      real, parameter :: stretch_factor = 1.5
      real, parameter :: target_lon = -97.5
      real, parameter :: target_lat = 35.5
      integer, parameter :: parent_tile = 6
      integer, parameter :: refine_ratio = 3

      integer, parameter :: cube_res = 96
      integer, parameter :: istart_nest = 26
      integer, parameter :: jstart_nest = 36
      integer, parameter :: iend_nest = 167
      integer, parameter :: jend_nest = 165

!     integer, parameter :: cube_res = 768
!     integer, parameter :: istart_nest = 191
!     integer, parameter :: jstart_nest = 327
!     integer, parameter :: iend_nest = 1346
!     integer, parameter :: jend_nest = 1290

      integer,parameter :: nvars_core=7                                 &  !<-- # of prognostic variables in core restart file
                          ,ndims_core=6                                 &  !<-- # of core restart dimensions
                          ,ndims_tracers=4                                 !<-- # of tracer restart dimensions
!
      real,parameter :: blend_exp1=0.5,blend_exp2=10.                      !<-- Define the exponential dropoff of weights
                                                                           !    for prescribed external values in the
                                                                           !    blending rows inside the domain boundary.
      real :: current_time_in_seconds
!
      integer,save :: isd_mod,ied_mod,jsd_mod,jed_mod
!
      integer,save :: ncid,next_time_to_read_bcs,nfields_tracers       &
                     ,npz,ntracers
!
      integer,save :: k_split,n_split
!
      integer,save :: bc_hour, ntimesteps_per_bc_update
!
      integer,save :: cld_amt_index                                    &   !<--
                     ,graupel_index                                    &   !
                     ,ice_water_index                                  &   !    Locations of 
                     ,liq_water_index                                  &   !    tracer vbls 
                     ,o3mr_index                                       &   !    in the tracers
                     ,rain_water_index                                 &   !    array.
                     ,snow_water_index                                 &   !
                     ,sphum_index                                          !<--
!
      integer,save :: lbnd_x_tracers,lbnd_y_tracers                    &   !<-- Local lower bounds of x,y for tracer arrays
                     ,ubnd_x_tracers,ubnd_y_tracers                        !<-- Local upper bounds of x,y for tracer arrays
!
      integer,save :: nrows_blend                                          !<-- # of blending rows in the BC data files.
!
      real,save :: dt_atmos                                            &   !<-- The physics (large) timestep (sec)
                  ,dyn_timestep                                            !<-- The dynamics timestep (sec)
!
      real(kind=R_GRID),dimension(:,:,:),allocatable :: agrid_reg      &   !<-- Lon/lat of cell centers
                                                       ,grid_reg           !<-- Lon/lat of cell corners

      real,dimension(:,:),allocatable :: phis_reg                          !<-- Filtered sfc geopotential

      real,dimension(:),allocatable :: ak_in, bk_in

      logical,save :: north_bc,south_bc,east_bc,west_bc                &
                     ,begin_regional_restart=.true.

      logical,dimension(:),allocatable,save :: blend_this_tracer

      character(len=50) :: filename_core='INPUT/fv_core.res.temp.nc'
      character(len=50) :: filename_core_new='RESTART/fv_core.res.tile1_new.nc'
      character(len=50) :: filename_tracers='INPUT/fv_tracer.res.temp.nc'
      character(len=50) :: filename_tracers_new='RESTART/fv_tracer.res.tile1_new.nc'

      type fv_regional_BC_variables
        real,dimension(:,:,:),allocatable :: delp_BC, divgd_BC, u_BC, v_BC, uc_BC, vc_BC
        real,dimension(:,:,:,:),allocatable :: q_BC
#ifndef SW_DYNAMICS
        real,dimension(:,:,:),allocatable :: pt_BC, w_BC, delz_BC
#ifdef USE_COND
        real,dimension(:,:,:),allocatable :: q_con_BC
#ifdef MOIST_CAPPA
        real,dimension(:,:,:),allocatable :: cappa_BC
#endif
#endif
#endif
      end type fv_regional_BC_variables

      type fv_domain_sides
        type(fv_regional_BC_variables) :: north, south, east, west
      end type fv_domain_sides

      type single_vbl3D_sides
        real,dimension(:,:,:),pointer :: north, south, east, west
      end type single_vbl3D_sides

      type vars_2d
        real,dimension(:,:),pointer :: ptr
        character(len=10) :: name
      end type vars_2d

      type vars_3d
        real,dimension(:,:,:),pointer :: ptr
        character(len=10) :: name
      end type vars_3d

      type(fv_domain_sides),target,save :: BC_t0, BC_t1  !<-- Boundary values for all BC variables at successive times from the regional BC file

      type(fv_regional_BC_variables),pointer,save :: bc_north_t0        &
                                                    ,bc_south_t0        &
                                                    ,bc_west_t0         &
                                                    ,bc_east_t0         &
                                                    ,bc_north_t1        &
                                                    ,bc_south_t1        &
                                                    ,bc_west_t1         &
                                                    ,bc_east_t1

      type(fv_regional_BC_variables),pointer :: bc_side_t0,bc_side_t1

      type(fv_regional_bc_bounds_type),pointer,save :: regional_bounds

      type(vars_3d),dimension(:),allocatable :: fields_core             &
                                               ,fields_tracers
      type(fv_nest_BC_type_3D), public :: delz_regBC ! lmh

      type(single_vbl3D_sides) :: delz_auxiliary    !<-- Boundary delz that follows integration through forecast time.

      integer :: ns = 0 ! lmh

      real,parameter :: tice=273.16                                     &
                       ,t_i0=15.
      real, parameter :: c_liq = 4185.5 ! gfdl: heat capacity of liquid at 15 deg c
      real, parameter :: c_ice = 1972.0 ! gfdl: heat capacity of ice at - 15 deg c
      real, parameter :: zvir = rvgas/rdgas - 1.                        &
                        ,cv_air = cp_air - rdgas                        &
                        ,cv_vap = cp_vapor - rvgas

      real,dimension(:),allocatable :: dum1d, pref

      character(len=100) :: grid_data='grid.tile7.halo4.nc'             &
                           ,oro_data ='oro_data.tile7.halo4.nc'

#ifdef OVERLOAD_R4
      real, parameter:: real_snan=x'FFBFFFFF'
#else
      real, parameter:: real_snan=x'FFF7FFFFFFFFFFFF'
#endif
      real(kind=R_GRID), parameter:: dbl_snan=x'FFF7FFFFFFFFFFFF'

      interface dump_field
        module procedure dump_field_3d
        module procedure dump_field_2d
      end interface dump_field

      integer,save :: bc_update_interval, nrows_blend_user

      integer :: a_step, p_step, k_step, n_step
!
      character(len=80) :: data_source
contains

!-----------------------------------------------------------------------
!
      subroutine setup_regional_BC(Atm, dt_atmos                        &
                                  ,isd,ied,jsd,jed                      &
                                  ,npx,npy )
!
!-----------------------------------------------------------------------
!***  Regional boundary data is obtained from the external BC file.
!-----------------------------------------------------------------------
      use netcdf
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!
!---------------------
!***  Input variables
!---------------------
!
      integer,intent(in) :: isd,ied,jsd,jed,npx,npy
!
      real,intent(in) :: dt_atmos                                          !<-- The large (physics) timestep (sec)
!
      type(fv_atmos_type),target,intent(inout) :: Atm                      !<-- Atm object for the current domain
!
!--------------------
!*** Local variables
!--------------------
!
      integer :: dimid,i,i_start,i_end,j,j_start,j_end,klev_out         &
                ,nrows_bc_data,nrows_blend_in_data,sec
!
      real :: ps1
!
      character(len=2) :: char2_1,char2_2
      character(len=3) :: int_to_char
      character(len=6) :: fmt='(i3.3)'
      character(len=50) :: file_name
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  The boundary data is laid out so that the pieces for the north
!***  and south sides span the entire distance from the east side of
!***  of the east halo to the west side of the west halo.  Therefore
!***  there the # of cells in the x direction in the north/south BC
!***  data is nx+2*nhalo where nx is the # of cells in the x direction
!***  on the compute domain.  This means the # of cells spanned in the
!***  west/east side BC data is just ny (the # of cells in the y
!***  direction on the compute domain) and not ny+2*nhalo since the
!***  halo values on the south and north ends of the east/west sides
!***  are already part of the BC data on the north/south sides.
!-----------------------------------------------------------------------
!
!                             nhalo_model=3
!
!                    |----------- nxp-1 -----------| <-- east/west compute points
!                 |---------- north BC data ----------|
!                 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!                 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!                 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!       ---       ooo           ---j=1---           ooo     ---         ---
!        |        ooo                               ooo      |           |
!        |        ooo                              |ooo      |           |
!                 ooo                        i=1-->|ooo
!   west BC data  ooo|                             |ooo east BC data    nyp-1 <-- north/south compute points
!                 ooo|<--i=isd-nhalo_model          ooo
!        |        ooo|                              ooo      |           |
!        |        ooo                               ooo      |           |
!       ---       ooo    ---j=jsd-nhalo_model---    ooo     ---         ---
!                 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!                 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!                 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!                 |---------- south BC data ----------|
!
!-----------------------------------------------------------------------
!
      north_bc=.false.
      south_bc=.false.
      east_bc =.false.
      west_bc =.false.
!
!  write(0,*)' enter setup_regional_BC isd=',isd,' ied=',ied,' jsd=',jsd,' jed=',jed    
      isd_mod=isd
      ied_mod=ied
      jsd_mod=jsd
      jed_mod=jed
!
!-----------------------------------------------------------------------
!***  Which side(s) of the domain does this task lie on if any?
!-----------------------------------------------------------------------
!
      if(jsd<0)then
        north_bc=.true.
      endif

      if(jed>npy-1)then
        south_bc=.true.
      endif

      if(isd<0)then
        east_bc=.true.
      endif

      if(ied>npx-1)then
        west_bc=.true.
      endif
!
      bc_update_interval=Atm%flagstruct%bc_update_interval
!
      k_split=Atm%flagstruct%k_split
      n_split=Atm%flagstruct%n_split
!
      dyn_timestep=dt_atmos/real(k_split*n_split)
!
!-----------------------------------------------------------------------
!***  Is blending row data present in the BC file and if so how many
!***  rows of data are there?  All blending data that is present will
!***  be read even if the user requests fewer rows be applied.
!***  Construct the name of the regional BC file to be read.
!***  We must know whether this is a standard BC file from chgres
!***  or a new BC file generated from DA-updated data from enlarged
!***  restart files that include the boundary rows.
!-----------------------------------------------------------------------
!
      write(int_to_char,fmt) bc_hour
      if(.not.Atm%flagstruct%regional_bcs_from_gsi)then
        file_name='INPUT/gfs_bndy.tile7.'//int_to_char//'.nc'              !<-- The standard BC file from chgres.
      else
        file_name='INPUT/gfs_bndy.tile7.'//int_to_char//'_gsi.nc'          !<-- The DA-updated BC file.
      endif
!
      if (is_master()) then
        write(*,20011)trim(file_name)
20011   format(' regional_bc_data file_name=',a)
      endif
!-----------------------------------------------------------------------
!***  Open the regional BC file.
!-----------------------------------------------------------------------
!
      call check(nf90_open(file_name,nf90_nowrite,ncid))                   !<-- Open the BC file; get the file ID.
      if (is_master()) then
        write(0,*)' opened BC file ',trim(file_name)
      endif
!
!-----------------------------------------------------------------------
!***  Check if the desired number of blending rows are present in
!***  the boundary files.
!-----------------------------------------------------------------------
!
      nrows_blend_user=Atm%flagstruct%nrows_blend                          !<-- # of blending rows the user wnats to apply.
!
      call check(nf90_inq_dimid(ncid,'halo',dimid))                        !<-- ID of the halo dimension.
      call check(nf90_inquire_dimension(ncid,dimid,len=nrows_bc_data))     !<-- Total # of rows of BC data (bndry + blending)
!
      nrows_blend_in_data=nrows_bc_data-nhalo_data                         !<-- # of blending rows in the BC files.
!
      if(nrows_blend_user>nrows_blend_in_data)then                         !<-- User wants more blending rows than are in the BC file.
        write(char2_1,'(I2.2)')nrows_blend_user
        write(char2_2,'(I2.2)')nrows_blend_in_data
        call mpp_error(FATAL,'User wants to use '//char2_1//' blending rows but only '//char2_2//' blending rows are in the BC file!')
      else
        nrows_blend=nrows_blend_in_data                                    !<-- # of blending rows in the BC files.
      endif
!
      call check(nf90_close(ncid))                                         !<-- Close the BC file for now.
!
!-----------------------------------------------------------------------
!***  Compute the index limits within the boundary region on each
!***  side of the domain for both scalars and winds.  Since the
!***  domain does not move then the computations need to be done
!***  only once.
!-----------------------------------------------------------------------
!
      call compute_regional_bc_indices(Atm%regional_bc_bounds)
!
!-----------------------------------------------------------------------
!
      if(.not.(north_bc.or.south_bc.or.east_bc.or.west_bc))then
        return                                                             !<-- This task is not on the domain boundary so exit.
      endif
!
!-----------------------------------------------------------------------
!
!     ntracers=Atm%ncnst - Atm%flagstruct%dnats                            !<-- # of advected tracers
      ntracers=Atm%ncnst                                                   !<-- Total # of tracers
      npz=Atm%npz                                                          !<-- # of layers in vertical configuration of integration
      klev_out=npz
!
      regional_bounds=>Atm%regional_bc_bounds
!
!-----------------------------------------------------------------------
!***  Allocate the objects that will hold the boundary variables
!***  at the two time levels surrounding each piece of the regional
!***  domain's integration.  Data is read from the BC files into
!***  time level t1 while t0 holds the data from the preceding 
!***  BC file.
!-----------------------------------------------------------------------
!***  Point pointers at each side's boundary data for both time levels.
!***  Those are needed when the actual update of boundary points is
!***  executed.
!-----------------------------------------------------------------------
!
      if(north_bc)then
        call allocate_regional_BC_arrays('north'                             &
                                        ,north_bc,south_bc                   &
                                        ,east_bc,west_bc                     &
                                        ,Atm%regional_bc_bounds%is_north     &
                                        ,Atm%regional_bc_bounds%ie_north     &
                                        ,Atm%regional_bc_bounds%js_north     &
                                        ,Atm%regional_bc_bounds%je_north     &
                                        ,Atm%regional_bc_bounds%is_north_uvs &
                                        ,Atm%regional_bc_bounds%ie_north_uvs &
                                        ,Atm%regional_bc_bounds%js_north_uvs &
                                        ,Atm%regional_bc_bounds%je_north_uvs &
                                        ,Atm%regional_bc_bounds%is_north_uvw &
                                        ,Atm%regional_bc_bounds%ie_north_uvw &
                                        ,Atm%regional_bc_bounds%js_north_uvw &
                                        ,Atm%regional_bc_bounds%je_north_uvw &
                                        ,klev_out                            &
                                        ,ntracers                            &
                                        ,BC_t1%north                         &
                                        ,delz_auxiliary%north )
!
        call allocate_regional_BC_arrays('north'                             &
                                        ,north_bc,south_bc                   &
                                        ,east_bc,west_bc                     &
                                        ,Atm%regional_bc_bounds%is_north     &
                                        ,Atm%regional_bc_bounds%ie_north     &
                                        ,Atm%regional_bc_bounds%js_north     &
                                        ,Atm%regional_bc_bounds%je_north     &
                                        ,Atm%regional_bc_bounds%is_north_uvs &
                                        ,Atm%regional_bc_bounds%ie_north_uvs &
                                        ,Atm%regional_bc_bounds%js_north_uvs &
                                        ,Atm%regional_bc_bounds%je_north_uvs &
                                        ,Atm%regional_bc_bounds%is_north_uvw &
                                        ,Atm%regional_bc_bounds%ie_north_uvw &
                                        ,Atm%regional_bc_bounds%js_north_uvw &
                                        ,Atm%regional_bc_bounds%je_north_uvw &
                                        ,klev_out                            &
                                        ,ntracers                            &
                                        ,BC_t0%north )
!
        bc_north_t0=>BC_t0%north
        bc_north_t1=>BC_t1%north
!
      endif

      if(south_bc)then
        call allocate_regional_BC_arrays('south'                             &
                                        ,north_bc,south_bc                   &
                                        ,east_bc,west_bc                     &
                                        ,Atm%regional_bc_bounds%is_south     &
                                        ,Atm%regional_bc_bounds%ie_south     &
                                        ,Atm%regional_bc_bounds%js_south     &
                                        ,Atm%regional_bc_bounds%je_south     &
                                        ,Atm%regional_bc_bounds%is_south_uvs &
                                        ,Atm%regional_bc_bounds%ie_south_uvs &
                                        ,Atm%regional_bc_bounds%js_south_uvs &
                                        ,Atm%regional_bc_bounds%je_south_uvs &
                                        ,Atm%regional_bc_bounds%is_south_uvw &
                                        ,Atm%regional_bc_bounds%ie_south_uvw &
                                        ,Atm%regional_bc_bounds%js_south_uvw &
                                        ,Atm%regional_bc_bounds%je_south_uvw &
                                        ,klev_out                            &
                                        ,ntracers                            &
                                        ,BC_t1%south                         &
                                        ,delz_auxiliary%south )
!
        call allocate_regional_BC_arrays('south'                             &
                                        ,north_bc,south_bc                   &
                                        ,east_bc,west_bc                     &
                                        ,Atm%regional_bc_bounds%is_south     &
                                        ,Atm%regional_bc_bounds%ie_south     &
                                        ,Atm%regional_bc_bounds%js_south     &
                                        ,Atm%regional_bc_bounds%je_south     &
                                        ,Atm%regional_bc_bounds%is_south_uvs &
                                        ,Atm%regional_bc_bounds%ie_south_uvs &
                                        ,Atm%regional_bc_bounds%js_south_uvs &
                                        ,Atm%regional_bc_bounds%je_south_uvs &
                                        ,Atm%regional_bc_bounds%is_south_uvw &
                                        ,Atm%regional_bc_bounds%ie_south_uvw &
                                        ,Atm%regional_bc_bounds%js_south_uvw &
                                        ,Atm%regional_bc_bounds%je_south_uvw &
                                        ,klev_out                            &
                                        ,ntracers                            &
                                        ,BC_t0%south )
!
        bc_south_t0=>BC_t0%south
        bc_south_t1=>BC_t1%south
!
      endif
!
      if(east_bc)then
        call allocate_regional_BC_arrays('east '                             &
                                        ,north_bc,south_bc                   &
                                        ,east_bc,west_bc                     &
                                        ,Atm%regional_bc_bounds%is_east      &
                                        ,Atm%regional_bc_bounds%ie_east      &
                                        ,Atm%regional_bc_bounds%js_east      &
                                        ,Atm%regional_bc_bounds%je_east      &
                                        ,Atm%regional_bc_bounds%is_east_uvs  &
                                        ,Atm%regional_bc_bounds%ie_east_uvs  &
                                        ,Atm%regional_bc_bounds%js_east_uvs  &
                                        ,Atm%regional_bc_bounds%je_east_uvs  &
                                        ,Atm%regional_bc_bounds%is_east_uvw  &
                                        ,Atm%regional_bc_bounds%ie_east_uvw  &
                                        ,Atm%regional_bc_bounds%js_east_uvw  &
                                        ,Atm%regional_bc_bounds%je_east_uvw  &
                                        ,klev_out                            &
                                        ,ntracers                            &
                                        ,BC_t1%east                          &
                                        ,delz_auxiliary%east )
!
        call allocate_regional_BC_arrays('east '                             &
                                        ,north_bc,south_bc                   &
                                        ,east_bc,west_bc                     &
                                        ,Atm%regional_bc_bounds%is_east      &
                                        ,Atm%regional_bc_bounds%ie_east      &
                                        ,Atm%regional_bc_bounds%js_east      &
                                        ,Atm%regional_bc_bounds%je_east      &
                                        ,Atm%regional_bc_bounds%is_east_uvs  &
                                        ,Atm%regional_bc_bounds%ie_east_uvs  &
                                        ,Atm%regional_bc_bounds%js_east_uvs  &
                                        ,Atm%regional_bc_bounds%je_east_uvs  &
                                        ,Atm%regional_bc_bounds%is_east_uvw  &
                                        ,Atm%regional_bc_bounds%ie_east_uvw  &
                                        ,Atm%regional_bc_bounds%js_east_uvw  &
                                        ,Atm%regional_bc_bounds%je_east_uvw  &
                                        ,klev_out                            &
                                        ,ntracers                            &
                                        ,BC_t0%east )
!
        bc_east_t0=>BC_t0%east
        bc_east_t1=>BC_t1%east
!
      endif
!
      if(west_bc)then
        call allocate_regional_BC_arrays('west '                             &
                                        ,north_bc,south_bc                   &
                                        ,east_bc,west_bc                     &
                                        ,Atm%regional_bc_bounds%is_west      &
                                        ,Atm%regional_bc_bounds%ie_west      &
                                        ,Atm%regional_bc_bounds%js_west      &
                                        ,Atm%regional_bc_bounds%je_west      &
                                        ,Atm%regional_bc_bounds%is_west_uvs  &
                                        ,Atm%regional_bc_bounds%ie_west_uvs  &
                                        ,Atm%regional_bc_bounds%js_west_uvs  &
                                        ,Atm%regional_bc_bounds%je_west_uvs  &
                                        ,Atm%regional_bc_bounds%is_west_uvw  &
                                        ,Atm%regional_bc_bounds%ie_west_uvw  &
                                        ,Atm%regional_bc_bounds%js_west_uvw  &
                                        ,Atm%regional_bc_bounds%je_west_uvw  &
                                        ,klev_out                            &
                                        ,ntracers                            &
                                        ,BC_t1%west                          &
                                        ,delz_auxiliary%west )
!
        call allocate_regional_BC_arrays('west '                             &
                                        ,north_bc,south_bc                   &
                                        ,east_bc,west_bc                     &
                                        ,Atm%regional_bc_bounds%is_west      &
                                        ,Atm%regional_bc_bounds%ie_west      &
                                        ,Atm%regional_bc_bounds%js_west      &
                                        ,Atm%regional_bc_bounds%je_west      &
                                        ,Atm%regional_bc_bounds%is_west_uvs  &
                                        ,Atm%regional_bc_bounds%ie_west_uvs  &
                                        ,Atm%regional_bc_bounds%js_west_uvs  &
                                        ,Atm%regional_bc_bounds%je_west_uvs  &
                                        ,Atm%regional_bc_bounds%is_west_uvw  &
                                        ,Atm%regional_bc_bounds%ie_west_uvw  &
                                        ,Atm%regional_bc_bounds%js_west_uvw  &
                                        ,Atm%regional_bc_bounds%je_west_uvw  &
                                        ,klev_out                            &
                                        ,ntracers                            &
                                        ,BC_t0%west )
!
        bc_west_t0=>BC_t0%west
        bc_west_t1=>BC_t1%west
!
      endif

      call allocate_fv_nest_BC_type(delz_regBC,Atm,ns,0,0,.false.)
!
!-----------------------------------------------------------------------
!***  We need regional versions of the arrays for surface elevation,
!***  latitude/longitude of grid cell corners, and lat/lon of the
!***  cell centers because those variables are needed an extra row
!***  beyond FV3's normal bounday region width of nhalo_model rows.
!-----------------------------------------------------------------------
!
      allocate(phis_reg(isd-1:ied+1,jsd-1:jed+1))   ; phis_reg=real_snan  !<-- Sfc elevation of filtered topography.
!
      allocate(agrid_reg(isd-1:ied+1,jsd-1:jed+1,2)); agrid_reg=dbl_snan  !<-- Center lat/lon of grid cells.
      allocate(grid_reg(isd-1:ied+2,jsd-1:jed+2,2)) ; grid_reg=dbl_snan   !<-- Lon/lat of grid cell corners.
!
!-----------------------------------------------------------------------
!***  From the data holding nhalo_model rows of boundary values 
!***  read in the lat/lon of the grid cell corners and fill in 
!***  the values of the grid cell centers.  The regional mode needs 
!***  the extra row of data.
!-----------------------------------------------------------------------
!
      call read_regional_lon_lat
!
!-----------------------------------------------------------------------
!***  From the data holding nhalo_model rows of filtered topography
!***  read in those values.  The regional mode needs the extra row
!***  of data.
!-----------------------------------------------------------------------
!
      call read_regional_filtered_topo
!
!-----------------------------------------------------------------------
!***  In the init step Atm%phis is given values only in the integration
!***  domain but in a regional run values are also needed in the
!***  boundary rows.  Since the same data is read in the preceding
!***  subroutine call as when Atm%phis was first filled, fill its
!***  boundary rows now.
!-----------------------------------------------------------------------
!
      if(north_bc)then
        i_start=isd
        i_end  =ied
        j_start=jsd
        if(.not.Atm%flagstruct%warm_start)then                             !<-- NOT a restarted run.
          j_end=jsd+nhalo_model-1
        else                                                               !<-- A restarted run.
          j_end=jsd+nhalo_model+1         
        endif
        do j=j_start,j_end
        do i=i_start,i_end
          Atm%phis(i,j)=phis_reg(i,j)
        enddo
        enddo
      endif
!
      if(south_bc)then
        i_start=isd
        i_end  =ied
        j_end  =jed
        if(.not.Atm%flagstruct%warm_start)then                             !<-- NOT a restarted run.
          j_start=jed-nhalo_model+1
        else                                                               !<-- A restarted run.
          j_start=jed-nhalo_model-1
        endif
        do j=j_start,j_end
        do i=i_start,i_end
          Atm%phis(i,j)=phis_reg(i,j)
        enddo
        enddo
      endif
      if(east_bc)then
        i_start=isd
        j_start=jsd
        j_end  =jed
        if(.not.Atm%flagstruct%warm_start)then                             !<-- NOT a restarted run.
          i_end=isd+nhalo_model-1
        else                                                               !<-- A restarted run.
          i_end=isd+nhalo_model+1
        endif
        do j=j_start,j_end
        do i=i_start,i_end
          Atm%phis(i,j)=phis_reg(i,j)
        enddo
        enddo
      endif
      if(west_bc)then
        i_end  =ied
        j_start=jsd
        j_end  =jed
        if(.not.Atm%flagstruct%warm_start)then                             !<-- NOT a restarted run.
          i_start=ied-nhalo_model+1
        else                                                               !<-- A restarted run.
          i_start=ied-nhalo_model-1
        endif
        do j=j_start,j_end
        do i=i_start,i_end
          Atm%phis(i,j)=phis_reg(i,j)
        enddo
        enddo
      endif
!
      sphum_index      = get_tracer_index(MODEL_ATMOS, 'sphum')
      liq_water_index  = get_tracer_index(MODEL_ATMOS, 'liq_wat')
      ice_water_index  = get_tracer_index(MODEL_ATMOS, 'ice_wat')
      rain_water_index = get_tracer_index(MODEL_ATMOS, 'rainwat')
      snow_water_index = get_tracer_index(MODEL_ATMOS, 'snowwat')
      graupel_index    = get_tracer_index(MODEL_ATMOS, 'graupel')
      cld_amt_index    = get_tracer_index(MODEL_ATMOS, 'cld_amt')
      o3mr_index       = get_tracer_index(MODEL_ATMOS, 'o3mr')
!  write(0,*)' setup_regional_bc'
!  write(0,*)' sphum_index=',sphum_index
!  write(0,*)' liq_water_index=',liq_water_index
!  write(0,*)' ice_water_index=',ice_water_index
!  write(0,*)' rain_water_index=',rain_water_index
!  write(0,*)' snow_water_index=',snow_water_index
!  write(0,*)' graupel_index=',graupel_index
!  write(0,*)' cld_amt_index=',cld_amt_index
!  write(0,*)' o3mr_index=',o3mr_index
!
!-----------------------------------------------------------------------
!***  When nudging of specific humidity is selected then we need a 
!***  reference pressure profile.  Compute it now.
!-----------------------------------------------------------------------
!
      allocate(pref(npz+1))
      allocate(dum1d(npz+1))
!
      ps1=101325.
      pref(npz+1)=ps1
      call get_eta_level(npz,ps1,pref(1),dum1d,Atm%ak,Atm%bk )
!
!-----------------------------------------------------------------------

      contains

!-----------------------------------------------------------------------
!
      subroutine compute_regional_bc_indices(regional_bc_bounds)
!
!-----------------------------------------------------------------------
!***  This routine computes the starting and ending indices for
!***  working arrays of task subdomains that lie on the edges
!***  of the FV3 regional domain.  These arrays will hold boundary
!***  region values of scalar variables located at the grid cell
!***  centers and wind components lying on the east/west sides
!***  and north/south sides of each cell.  Note that the width
!***  of the domain's boundary region (4 rows) is currently 
!***  greater than the fundamental width of the task subdomains'
!***  halo regions (3 rows).  The variables isd,ied,jsd,jed are
!***  the task subdomain index limits including their halos.
!***  The diagram in subroutine regional_bc_data will help to
!***  understand these index limits have the values they do.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      type(fv_regional_bc_bounds_type),intent(out) :: regional_bc_bounds
!
!---------------------
!***  Local variables
!---------------------
!
      integer, parameter :: invalid_index = -99
      integer :: halo_diff
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      regional_bc_bounds%is_north = invalid_index
      regional_bc_bounds%ie_north = invalid_index
      regional_bc_bounds%js_north = invalid_index
      regional_bc_bounds%je_north = invalid_index
      regional_bc_bounds%is_north_uvs = invalid_index
      regional_bc_bounds%ie_north_uvs = invalid_index
      regional_bc_bounds%js_north_uvs = invalid_index
      regional_bc_bounds%je_north_uvs = invalid_index
      regional_bc_bounds%is_north_uvw = invalid_index
      regional_bc_bounds%ie_north_uvw = invalid_index
      regional_bc_bounds%js_north_uvw = invalid_index
      regional_bc_bounds%je_north_uvw = invalid_index

      regional_bc_bounds%is_south = invalid_index
      regional_bc_bounds%ie_south = invalid_index
      regional_bc_bounds%js_south = invalid_index
      regional_bc_bounds%je_south = invalid_index
      regional_bc_bounds%is_south_uvs = invalid_index
      regional_bc_bounds%ie_south_uvs = invalid_index
      regional_bc_bounds%js_south_uvs = invalid_index
      regional_bc_bounds%je_south_uvs = invalid_index
      regional_bc_bounds%is_south_uvw = invalid_index
      regional_bc_bounds%ie_south_uvw = invalid_index
      regional_bc_bounds%js_south_uvw = invalid_index
      regional_bc_bounds%je_south_uvw = invalid_index

      regional_bc_bounds%is_east = invalid_index
      regional_bc_bounds%ie_east = invalid_index
      regional_bc_bounds%js_east = invalid_index
      regional_bc_bounds%je_east = invalid_index
      regional_bc_bounds%is_east_uvs = invalid_index
      regional_bc_bounds%ie_east_uvs = invalid_index
      regional_bc_bounds%js_east_uvs = invalid_index
      regional_bc_bounds%je_east_uvs = invalid_index
      regional_bc_bounds%is_east_uvw = invalid_index
      regional_bc_bounds%ie_east_uvw = invalid_index
      regional_bc_bounds%js_east_uvw = invalid_index
      regional_bc_bounds%je_east_uvw = invalid_index

      regional_bc_bounds%is_west = invalid_index
      regional_bc_bounds%ie_west = invalid_index
      regional_bc_bounds%js_west = invalid_index
      regional_bc_bounds%je_west = invalid_index
      regional_bc_bounds%is_west_uvs = invalid_index
      regional_bc_bounds%ie_west_uvs = invalid_index
      regional_bc_bounds%js_west_uvs = invalid_index
      regional_bc_bounds%je_west_uvs = invalid_index
      regional_bc_bounds%is_west_uvw = invalid_index
      regional_bc_bounds%ie_west_uvw = invalid_index
      regional_bc_bounds%js_west_uvw = invalid_index
      regional_bc_bounds%je_west_uvw = invalid_index
!
!-----------------------------------------------------------------------
!***  Scalar BC indices
!-----------------------------------------------------------------------
!***  These must reach one row beyond nhalo_model since we must
!***  surround the wind points on the cell edges with mass points.
!
!***  NOTE:  The value of nrows_blend is the total number of
!***          blending rows in the BC files.
!-----------------------------------------------------------------------
!
      halo_diff=nhalo_data-nhalo_model
!
!-----------
!***  North
!-----------
!
      if (north_bc) then
        regional_bc_bounds%is_north=isd-1
        regional_bc_bounds%ie_north=ied+1
!
        regional_bc_bounds%js_north=jsd-1
        regional_bc_bounds%je_north=nrows_blend
      endif
!
!-----------
!***  South
!-----------
!
      if (south_bc) then
        regional_bc_bounds%is_south=isd-1
        regional_bc_bounds%ie_south=ied+1
!
        regional_bc_bounds%js_south=jed-nhalo_model-nrows_blend+1
        regional_bc_bounds%je_south=jed+1
      endif
!
!----------
!***  East
!----------
!
      if (east_bc) then
        regional_bc_bounds%is_east=isd-1
        regional_bc_bounds%ie_east=nrows_blend
!
        regional_bc_bounds%js_east=jsd-1
        if(north_bc)then
          regional_bc_bounds%js_east=1
        endif
!
        regional_bc_bounds%je_east=jed+1
        if(south_bc)then
          regional_bc_bounds%je_east=jed-nhalo_model
        endif
      endif
!
!----------
!***  West
!----------
!
      if (west_bc) then
        regional_bc_bounds%is_west=ied-nhalo_model-nrows_blend+1
        regional_bc_bounds%ie_west=ied+1
!
        regional_bc_bounds%js_west=jsd-1
        if(north_bc)then
          regional_bc_bounds%js_west=1
        endif
!
        regional_bc_bounds%je_west=jed+1
        if(south_bc)then
          regional_bc_bounds%je_west=jed-nhalo_model
        endif
      endif
!
!-----------------------------------------------------------------------
!*** Wind component BC indices
!-----------------------------------------------------------------------
!
!-----------
!***  North
!-----------
!
      if (north_bc) then
        regional_bc_bounds%is_north_uvs=isd
        regional_bc_bounds%ie_north_uvs=ied
!
        regional_bc_bounds%js_north_uvs=jsd
        regional_bc_bounds%je_north_uvs=nrows_blend+1
!
        regional_bc_bounds%is_north_uvw=isd
        regional_bc_bounds%ie_north_uvw=ied+1
!
        regional_bc_bounds%js_north_uvw=jsd
        regional_bc_bounds%je_north_uvw=nrows_blend
      endif
!
!-----------
!***  South
!-----------
!
      if (south_bc) then
        regional_bc_bounds%is_south_uvs=isd
        regional_bc_bounds%ie_south_uvs=ied
!
        regional_bc_bounds%js_south_uvs=jed-nhalo_model-nrows_blend+1
        regional_bc_bounds%je_south_uvs=jed+1
!
        regional_bc_bounds%is_south_uvw=isd
        regional_bc_bounds%ie_south_uvw=ied+1
!
        regional_bc_bounds%js_south_uvw=jed-nhalo_model-nrows_blend+1
        regional_bc_bounds%je_south_uvw=jed
      endif
!
!----------
!***  East
!----------
!
      if (east_bc) then
        regional_bc_bounds%is_east_uvs=isd
        regional_bc_bounds%ie_east_uvs=nrows_blend
!
        regional_bc_bounds%js_east_uvs=jsd
        if(north_bc)then
          regional_bc_bounds%js_east_uvs=1  !<-- north side of cell at j=1 (north bdry contains north side of j=1)
        endif
!
        regional_bc_bounds%je_east_uvs=jed+1
        if(south_bc)then
          regional_bc_bounds%je_east_uvs=jed-nhalo_model+1
        endif
!
!       regional_bc_bounds%is_east_uvw=isd-1
        regional_bc_bounds%is_east_uvw=isd
        regional_bc_bounds%ie_east_uvw=nrows_blend     !<-- east side of cell at i=nrows_blend
!
        regional_bc_bounds%js_east_uvw=jsd
        if(north_bc)then
          regional_bc_bounds%js_east_uvw=1
        endif
        regional_bc_bounds%je_east_uvw=jed
        if(south_bc)then
          regional_bc_bounds%je_east_uvw=jed-nhalo_model
        endif
      endif
!
!----------
!***  West
!----------
!
      if (west_bc) then
        regional_bc_bounds%is_west_uvs=ied-nhalo_model-nrows_blend+1
        regional_bc_bounds%ie_west_uvs=ied
!
        regional_bc_bounds%js_west_uvs=jsd
        if(north_bc)then
          regional_bc_bounds%js_west_uvs=1
        endif
!
        regional_bc_bounds%je_west_uvs=jed+1
        if(south_bc)then
          regional_bc_bounds%je_west_uvs=jed-nhalo_model+1
        endif
!
        regional_bc_bounds%is_west_uvw=ied-nhalo_model-nrows_blend+1
        regional_bc_bounds%ie_west_uvw=ied+1
!
        regional_bc_bounds%js_west_uvw=jsd
        if(north_bc)then
          regional_bc_bounds%js_west_uvw=1
        endif
!
        regional_bc_bounds%je_west_uvw=jed
        if(south_bc)then
          regional_bc_bounds%je_west_uvw=jed-nhalo_model
        endif
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine compute_regional_bc_indices
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      subroutine read_regional_lon_lat
!
!-----------------------------------------------------------------------
!***  Read the longitude/latitude of the grid cell corners from
!***  the external file holding the additional row of data required
!***  by the regional domain.
!-----------------------------------------------------------------------
      use netcdf
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!
!--------------------
!***  Local variables
!--------------------
!
      integer :: i_start_data,istat,j_start_data,n,ncid_grid,var_id
!
      character(len=150) :: filename,vname
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Open the grid data file.
!-----------------------------------------------------------------------
!
      filename='INPUT/'//trim(grid_data)
!
      call check(nf90_open(filename,nf90_nowrite,ncid_grid))               !<-- Open the grid data netcdf file; get the file ID.
!
      call  mpp_error(NOTE,' opened grid file '//trim(filename))
!
!-----------------------------------------------------------------------
!***  The longitude and latitude are on the super grid.  We need only
!***  the points on each corner of the grid cells which is every other
!***  point on the super grid.
!-----------------------------------------------------------------------
!
      i_start_data=2*(isd+nhalo_model)-1
      j_start_data=2*(jsd+nhalo_model)-1
!
!---------------
!***  Longitude
!---------------
!
      vname='x'                                                            !<-- Geographic_longitude (degrees east) in netcdf file
      call check(nf90_inq_varid(ncid_grid,vname,var_id))                   !<-- Get the variable ID.
      call check(nf90_get_var(ncid_grid,var_id                          &
                             ,grid_reg(isd-1:ied+2,jsd-1:jed+2,1)       &  !<-- Longitude of grid cell corners
                             ,start=(/i_start_data,j_start_data/)       &
                             ,stride=(/2,2/) ) )
!
!--------------
!***  Latitude
!--------------
!
      vname='y'                                                            !<-- Geographic_latitude (degrees north) in netcdf file
      call check(nf90_inq_varid(ncid_grid,vname,var_id))                   !<-- Get the variable ID.
      call check(nf90_get_var(ncid_grid,var_id                          &
                             ,grid_reg(isd-1:ied+2,jsd-1:jed+2,2)       &  !<-- Latitude of grid cell corners
                             ,start=(/i_start_data,j_start_data/)       &
                             ,stride=(/2,2/) ) )
!
      call check(nf90_close(ncid_grid))
!
!-----------------------------------------------------------------------
!***  Convert from degrees to radians.
!-----------------------------------------------------------------------
!
      do n=1,2
        do j=jsd-1,jed+2
        do i=isd-1,ied+2
          grid_reg(i,j,n)=grid_reg(i,j,n)*pi/180.
        enddo
        enddo
      enddo
!
!-----------------------------------------------------------------------
!***  Compute the longitude/latitude in the cell centers.
!-----------------------------------------------------------------------
!
      do j=jsd-1,jed+1
      do i=isd-1,ied+1
        call cell_center2(grid_reg(i,j,  1:2), grid_reg(i+1,j,  1:2),   &
                          grid_reg(i,j+1,1:2), grid_reg(i+1,j+1,1:2),   &
                          agrid_reg(i,j,1:2) )
      enddo
      enddo
!
!-----------------------------------------------------------------------
!
      end subroutine read_regional_lon_lat
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      subroutine read_regional_filtered_topo
!
!-----------------------------------------------------------------------
!***  Read the filtered topography including the extra outer row.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!
!---------------------
!***  Local variables
!---------------------
!
      integer :: i,i_start_data,istat,j,j_start_data,ncid_oro,var_id
!
      character(len=150) :: filename,vname
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Get the name of the working directory.  Open the topography data
!***  file.
!-----------------------------------------------------------------------
!
      filename='INPUT/'//trim(oro_data)

      if (is_master()) then
        write(*,23421)trim(filename)
23421   format(' topo filename=',a)
      endif
!
      call check(nf90_open(filename,nf90_nowrite,ncid_oro))                !<-- Open the netcdf file; get the file ID.
!
!-----------------------------------------------------------------------
!***  Read in the data including the extra outer row.
!-----------------------------------------------------------------------
!
      i_start_data=isd+nhalo_model
      j_start_data=jsd+nhalo_model
!
      vname='orog_filt'                                                    !<-- Filtered topography (m) in netcdf file
      call check(nf90_inq_varid(ncid_oro,vname,var_id))                    !<-- Get the variable ID.
      call check(nf90_get_var(ncid_oro,var_id                           &
                             ,phis_reg(isd-1:ied+1,jsd-1:jed+1)         &  !<-- Extracted filtered topography (m)
                             ,start=(/i_start_data,j_start_data/)))
!
      call check(nf90_close(ncid_oro))
!
!-----------------------------------------------------------------------
!***  We want the geopotential.
!-----------------------------------------------------------------------
!
      do j=jsd-1,jed+1
      do i=isd-1,ied+1
        phis_reg(i,j)=phis_reg(i,j)*grav
      enddo
      enddo
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      end subroutine read_regional_filtered_topo
!
!-----------------------------------------------------------------------
!
      end subroutine setup_regional_BC
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      subroutine start_regional_cold_start(Atm, dt_atmos, ak, bk, levp  &
                                          ,is ,ie ,js ,je               &
                                          ,isd,ied,jsd,jed )
!
!-----------------------------------------------------------------------
!***  Prepare the regional run for a cold start.
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      type(fv_atmos_type),intent(inout) :: Atm                             !<-- Atm object for the current domain
!
      integer ,intent(in) :: is ,ie ,js ,je                             &  !<-- Integration limits of task subdomain
                            ,isd,ied,jsd,jed                            &  !<-- Memory limits of task subdomain
                            ,levp 
!
      real,intent(in) :: dt_atmos                                          !<-- The large (physics) timestep (sec)
      real,intent(in) :: ak(1:levp+1), bk(1:levp+1)
!
!---------------------
!***  Local variables
!---------------------
!
      integer :: k
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Get the source of the input data
!-----------------------------------------------------------------------
!
      call get_data_source(data_source,Atm%flagstruct%regional)
!
      call setup_regional_BC(Atm, dt_atmos                              &
                            ,isd, ied, jsd, jed                         &
                            ,Atm%npx, Atm%npy )
!
      bc_hour=0
      call regional_bc_data(Atm, bc_hour                                &  !<-- Fill time level t1 from BC file at 0 hours.
                           ,is, ie, js, je                              &
                           ,isd, ied, jsd, jed                          &
                           ,ak, bk )
      call regional_bc_t1_to_t0(BC_t1, BC_t0                            &  !
                               ,Atm%npz                                 &  !<-- Move BC t1 data to t0.
                               ,ntracers                                &
                               ,Atm%regional_bc_bounds )                   !
!
      bc_hour=bc_hour+bc_update_interval
!
!-----------------------------------------------------------------------
!***  If this is a DA run and the first BC file was updated by
!***  the GSI then reset the gsi flag so that all subsequent
!***  BC files are read normally.
!-----------------------------------------------------------------------
!
      if(Atm%flagstruct%regional_bcs_from_gsi)then
        Atm%flagstruct%regional_bcs_from_gsi=.false.
      endif
!
      call regional_bc_data(Atm, bc_hour                                &  !<-- Fill time level t1 
                           ,is, ie, js, je                              &  !    from the 2nd time level
                           ,isd, ied, jsd, jed                          &  !    in the BC file.
                           ,ak, bk )                                       !
!
      allocate (ak_in(1:levp+1))                                           !<-- Save the input vertical structure for
      allocate (bk_in(1:levp+1))                                           !    remapping BC updates during the forecast.
      do k=1,levp+1
        ak_in(k)=ak(k)
        bk_in(k)=bk(k)
      enddo
!
!-----------------------------------------------------------------------
!***  If the GSI will need restart files that includes the
!***  fields' boundary rows.  Those files were already created.
!***  Prepare the objects that hold their variables' names and
!***  values.
!-----------------------------------------------------------------------
!
      if(Atm%flagstruct%write_restart_with_bcs)then
        call prepare_full_fields(Atm)
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine start_regional_cold_start
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      subroutine start_regional_restart(Atm, dt_atmos                   &
                                       ,isc,iec,jsc,jec                 &
                                       ,isd,ied,jsd,jed )
!
!-----------------------------------------------------------------------
!***  Prepare the regional forecast for a restart.
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      type(fv_atmos_type),intent(inout) :: Atm                             !<-- Atm object for the current domain
!
      real,intent(in) :: dt_atmos                                          !<-- The large (physics) timestep (sec)
!
      integer ,intent(in) :: isc,iec,jsc,jec                            &  !<-- Integration limits of task subdomain
                            ,isd,ied,jsd,jed                               !<-- Memory limits of task subdomain
!
!---------------------
!***  Local variables
!---------------------
!
      integer :: ierr, ios
      real, allocatable :: wk2(:,:)
!
      logical :: filtered_terrain = .true.
      logical :: gfs_dwinds       = .true.
      integer :: levp             = 64
      logical :: checker_tr       = .false.
      integer :: nt_checker       = 0
      namelist /external_ic_nml/ filtered_terrain, levp, gfs_dwinds     &
                                ,checker_tr, nt_checker
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!*** Read the number of model layers in the external forecast (=levp).
!-----------------------------------------------------------------------
!
      read (input_nml_file,external_ic_nml,iostat=ios)
      ierr = check_nml_error(ios,'external_ic_nml')
      if(ierr/=0)then
        write(0,11011)ierr
11011   format(' start_regional_restart failed to read external_ic_nml ierr=',i3)
      endif
!
!-----------------------------------------------------------------------
!***  Get the source of the input data.
!-----------------------------------------------------------------------
!
      call get_data_source(data_source,Atm%flagstruct%regional)
!
!-----------------------------------------------------------------------
!***  Preliminary setup for the forecast.
!-----------------------------------------------------------------------
!
      call setup_regional_BC(Atm, dt_atmos                              &
                            ,isd, ied, jsd, jed                         &
                            ,Atm%npx, Atm%npy )
!
      allocate (wk2(levp+1,2))
      allocate (ak_in(levp+1))                                             !<-- Save the input vertical structure for
      allocate (bk_in(levp+1))                                             !    remapping BC updates during the forecast.
      call read_data('INPUT/gfs_ctrl.nc','vcoord',wk2, no_domain=.TRUE.)
      ak_in(1:levp+1) = wk2(1:levp+1,1)
      ak_in(1) = 1.e-9
      bk_in(1:levp+1) = wk2(1:levp+1,2)
      deallocate(wk2)
      bc_hour=nint(current_time_in_seconds/3600.)
!
!-----------------------------------------------------------------------
!***  Fill time level t1 from the BC file at the restart time.
!-----------------------------------------------------------------------
!
      call regional_bc_data(Atm, bc_hour                                & 
                           ,isc, iec, jsc, jec                          & 
                           ,isd, ied, jsd, jed                          & 
                           ,ak_in, bk_in )
!
!-----------------------------------------------------------------------
!***  If this is a DA run and the first BC file was updated by
!***  the GSI then that file was read differently in the preceding
!***  call to subroutine regional_bc_data.  Now reset the gsi
!***  flag so that all subsequent BC files are read normally.
!-----------------------------------------------------------------------
!
      if(Atm%flagstruct%regional_bcs_from_gsi)then
        Atm%flagstruct%regional_bcs_from_gsi=.false.
      endif
!
!-----------------------------------------------------------------------
!***  If the GSI will need restart files that include the
!***  fields' boundary rows after this forecast or forecast
!***  segment completes then prepare the objects that will
!***  hold their variables' names and values.
!-----------------------------------------------------------------------
!
      if(Atm%flagstruct%write_restart_with_bcs)then
        call prepare_full_fields(Atm)
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine start_regional_restart
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      subroutine read_new_bc_data(Atm, Time, Time_step_atmos, p_split   &
                                 ,isd,ied,jsd,jed )
!
!-----------------------------------------------------------------------
!***  When it is time to read new boundary data from the external files
!***  move time level t1 to t0 and then read the data into t1.
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      type(fv_atmos_type),intent(inout) :: Atm                             !<-- Atm object for the current domain
      type(time_type),intent(in) :: Time                                   !<-- Current forecast time
      type (time_type),intent(in) :: Time_step_atmos                       !<-- Large (physics) timestep

      integer,intent(in) :: isd,ied,jsd,jed                             &  !<-- Memory limits of task subdomain
                           ,p_split
!
!---------------------
!***  Local variables
!---------------------
!
      integer :: atmos_time_step, sec
      real :: dt_atmos
      type(time_type) :: atmos_time
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      atmos_time = Time - Atm%Time_init
      atmos_time_step = atmos_time / Time_step_atmos
      current_time_in_seconds = time_type_to_real( atmos_time )
      if (mpp_pe() == 0 .and. Atm%flagstruct%fv_debug) write(*,"('current_time_seconds = ',f9.1)")current_time_in_seconds
!
      call get_time (Time_step_atmos, sec)
      dt_atmos = real(sec)
!
      if(atmos_time_step==0.or.Atm%flagstruct%warm_start)then
        ntimesteps_per_bc_update=nint(bc_update_interval*3600./(dt_atmos/real(abs(p_split))))
      endif
!
      if(atmos_time_step+1>=ntimesteps_per_bc_update.and.mod(atmos_time_step,ntimesteps_per_bc_update)==0 &
                                                    .or.                                                  &
         Atm%flagstruct%warm_start.and.begin_regional_restart)then
!
        begin_regional_restart=.false.
        bc_hour=bc_hour+bc_update_interval
!
!-----------------------------------------------------------------------
!***  Transfer the time level t1 data to t0.
!-----------------------------------------------------------------------
!
        call regional_bc_t1_to_t0(BC_t1, BC_t0                          &  
                                 ,Atm%npz                               & 
                                 ,ntracers                              &
                                 ,Atm%regional_bc_bounds )
!
!-----------------------------------------------------------------------
!***  Fill time level t1 from the BC file containing data from
!***  the next time level.
!-----------------------------------------------------------------------
!
        call regional_bc_data(Atm, bc_hour                              & 
                             ,Atm%bd%is, Atm%bd%ie                      &
                             ,Atm%bd%js, Atm%bd%je                      & 
                             ,isd, ied, jsd, jed                        &
                             ,ak_in, bk_in )
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine read_new_bc_data
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      subroutine regional_bc_data(Atm,bc_hour                           &
                                 ,is,ie,js,je                           &
                                 ,isd,ied,jsd,jed                       &
                                 ,ak,bk )
!
!-----------------------------------------------------------------------
!***  Regional boundary data is obtained from the external BC file.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
!-----------
!***  Input
!-----------
!
      integer,intent(in) :: bc_hour                                        !<-- The forecast hour of the BC file to be read.
!
      integer,intent(in) :: is,ie,js,je                                 &  !<-- Compute limits of task subdomain
                           ,isd,ied,jsd,jed                                !<-- Halo limits of task subdomain
!
      real,dimension(:),intent(in) :: ak,bk
!
!-----------------
!*** Input/output
!-----------------
!
      type(fv_atmos_type),target,intent(inout) :: Atm                      !<-- Atm object for the current domain
!
!---------------------
!***  Local variables
!---------------------
!
      integer :: dimid,i,j,k,klev_in,klev_out,n,nlev
!
      integer :: is_north,is_south,is_east,is_west                      &
                ,ie_north,ie_south,ie_east,ie_west                      &
                ,js_north,js_south,js_east,js_west                      &
                ,je_north,je_south,je_east,je_west
!
      integer :: is_u,ie_u,js_u,je_u                                    &
                ,is_v,ie_v,js_v,je_v
!
      integer :: is_input,ie_input,js_input,je_input
!
      integer :: i_start,i_end,j_start,j_end
!
      integer :: nside,nt,index
!
      real,dimension(:,:,:),allocatable :: ud,vd,uc,vc
!
      real,dimension(:,:),allocatable :: ps_reg
      real,dimension(:,:,:),allocatable :: delp_input,delz_input        &
                                          ,ps_input,t_input             &
                                          ,w_input,zh_input
      real,dimension(:,:,:),allocatable :: u_s_input,v_s_input          &
                                          ,u_w_input,v_w_input
      real,dimension(:,:,:,:),allocatable :: tracers_input
!
      real(kind=R_GRID), dimension(2):: p1, p2, p3, p4
      real(kind=R_GRID), dimension(3):: e1, e2, ex, ey

#undef USE_FMS_READ
#ifdef USE_FMS_READ
      integer :: isc2, iec2, jsc2, jec2
      real(kind=R_GRID), allocatable, dimension(:,:)  :: tmpx, tmpy
      integer :: start(4), nread(4)
      real(kind=R_GRID), allocatable, dimension(:,:,:) :: reg_grid
      real(kind=R_GRID), allocatable, dimension(:,:,:) :: reg_agrid
#endif
!
      logical,save :: computed_regional_bc_indices=.false.
!
      character(len=3) :: int_to_char
      character(len=5) :: side
      character(len=6) :: fmt='(i3.3)'
!
      character(len=50) :: file_name
!
      integer,save :: kount1=0,kount2=0
      integer :: istart, iend, jstart, jend
      integer :: npx, npy
!
      character(len=60) :: var_name_root
      logical :: required
!
      logical :: call_remap
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Only boundary tasks are needed.
!-----------------------------------------------------------------------
!
      if(.not.(north_bc.or.south_bc.or.east_bc.or.west_bc))then
        return
      endif
!
!-----------------------------------------------------------------------
!
      klev_out=Atm%npz                                                     !<-- # of layers in vertical configuration of integration
!
!-----------------------------------------------------------------------
!***  Construct the name of the regional BC file to be read.
!***  We must know whether this is a standard BC file from chgres
!***  or a new BC file generated from DA-updated data from enlarged
!***  restart files that include the boundary rows.
!-----------------------------------------------------------------------
!
      write(int_to_char,fmt) bc_hour
      if(.not.Atm%flagstruct%regional_bcs_from_gsi)then
        file_name='INPUT/gfs_bndy.tile7.'//int_to_char//'.nc'              !<-- The standard BC file from chgres.
      else
        file_name='INPUT/gfs_bndy.tile7.'//int_to_char//'_gsi.nc'          !<-- The DA-updated BC file.
      endif
!
      if (is_master()) then
        write(*,22211)trim(file_name)
22211   format(' regional_bc_data file_name=',a)
      endif
!-----------------------------------------------------------------------
!***  Open the regional BC file.
!***  Find the # of layers (klev_in) in the BC input.
!-----------------------------------------------------------------------
!
      call check(nf90_open(file_name,nf90_nowrite,ncid))                   !<-- Open the BC file; get the file ID.
      if (is_master()) then
        write(0,*)' opened BC file ',trim(file_name)
      endif
!
      call check(nf90_inq_dimid(ncid,'lev',dimid))                         !<-- Get the vertical dimension's NetCDF ID.
      call check(nf90_inquire_dimension(ncid,dimid,len=klev_in))           !<-- Get the vertical dimension's value (klev_in).
!
!-----------------------------------------------------------------------
!***  Allocate the boundary variables and initialize them to garbage.
!-----------------------------------------------------------------------
!
      is_input=is-nhalo_data
      ie_input=ie+nhalo_data
      js_input=js-nhalo_data
      je_input=je+nhalo_data
      npx = Atm%npx
      npy = Atm%npy
!
      allocate( ps_input(is_input:ie_input,js_input:je_input,1)) ; ps_input=real_snan                 !<-- Sfc pressure
      allocate(  t_input(is_input:ie_input,js_input:je_input,1:klev_in)) ; t_input=real_snan          !<-- Sensible temperature
      allocate(  w_input(is_input:ie_input,js_input:je_input,1:klev_in)) ; w_input=real_snan          !<-- Vertical velocity
      allocate(u_s_input(is_input:ie_input,js_input:je_input,1:klev_in)) ; u_s_input=real_snan        !<-- D-grid u component
      allocate(v_s_input(is_input:ie_input,js_input:je_input,1:klev_in)) ; v_s_input=real_snan        !<-- C-grid v component
      allocate(u_w_input(is_input:ie_input,js_input:je_input,1:klev_in)) ; u_w_input=real_snan        !<-- C-grid u component
      allocate(v_w_input(is_input:ie_input,js_input:je_input,1:klev_in)) ; v_w_input=real_snan        !<-- D-grid v component
!
      if(Atm%flagstruct%regional_bcs_from_gsi)then
        allocate(delp_input(is_input:ie_input,js_input:je_input,1:klev_in)) ; delp_input=real_snan    !<-- Lyr pressure depth (Pa)
        allocate(delz_input(is_input:ie_input,js_input:je_input,1:klev_in)) ; delz_input=real_snan    !<-- Lyr geometric depth (m)
      else
        allocate( zh_input(is_input:ie_input,js_input:je_input,1:klev_in+1)) ; zh_input=real_snan     !<-- Lyr interface heights (m)
      endif
!
      allocate(tracers_input(is_input:ie_input,js_input:je_input,klev_in,ntracers)) ; tracers_input=real_snan
!
!-----------------------------------------------------------------------
!***  Extract each variable from the regional BC file.  The final
!***  argument is the object being filled.
!-----------------------------------------------------------------------
!
!------------------
!***  Sfc pressure
!------------------
!
      nlev=1
      var_name_root='ps'
      call read_regional_bc_file(is_input,ie_input,js_input,je_input    &
                                ,nlev                                   &
                                ,ntracers                               &
!                               ,Atm%regional_bc_bounds                 &
                                ,var_name_root                          &
                                ,array_3d=ps_input )                       !<-- ps is 2D but for simplicity here use a 3rd dim of 1
!
!-----------------------
!***  Vertical velocity
!-----------------------
!
      nlev=klev_in
      var_name_root='w'
      call read_regional_bc_file(is_input,ie_input,js_input,je_input    &
                                ,nlev                                   &
                                ,ntracers                               &
!                               ,Atm%regional_bc_bounds                 &
                                ,var_name_root                          &
                                ,array_3d=w_input)
!
!-----------------------
!***  Interface heights
!-----------------------
!
      if(.not.Atm%flagstruct%regional_bcs_from_gsi)then
        nlev=klev_in+1
        var_name_root='zh'
        call read_regional_bc_file(is_input,ie_input,js_input,je_input  &
                                  ,nlev                                 &
                                  ,ntracers                             &
                                  ,var_name_root                        &
                                  ,array_3d=zh_input)
      endif
!
!--------------------------
!***  Sensible temperature
!--------------------------
!
      if (data_source == 'FV3GFS GAUSSIAN NEMSIO FILE') then
        nlev=klev_in
        var_name_root='t'
        call read_regional_bc_file(is_input,ie_input,js_input,je_input  &
                                  ,nlev                                 &
                                  ,ntracers                             &
!                                 ,Atm%regional_bc_bounds               &
                                  ,var_name_root                        &
                                  ,array_3d=t_input)
      endif
!
!-----------------------------
!***  U component south/north
!-----------------------------
!
      nlev=klev_in
      var_name_root='u_s'
      call read_regional_bc_file(is_input,ie_input,js_input,je_input    &
                                ,nlev                                   &
                                ,ntracers                               &
!                               ,Atm%regional_bc_bounds                 &
                                ,var_name_root                          &
                                ,array_3d=u_s_input)
!
!-----------------------------
!***  V component south/north
!-----------------------------
!
      nlev=klev_in
      var_name_root='v_s'
      call read_regional_bc_file(is_input,ie_input,js_input,je_input    &
                                ,nlev                                   &
                                ,ntracers                               &
!                               ,Atm%regional_bc_bounds                 &
                                ,var_name_root                          &
                                ,array_3d=v_s_input)
!
!---------------------------
!***  U component east/west
!---------------------------
!
      nlev=klev_in
      var_name_root='u_w'
      call read_regional_bc_file(is_input,ie_input,js_input,je_input    &
                                ,nlev                                   &
                                ,ntracers                               &
!                               ,Atm%regional_bc_bounds                 &
                                ,var_name_root                          &
                                ,array_3d=u_w_input)
!
!---------------------------
!***  V component east/west
!---------------------------
!
      nlev=klev_in
      var_name_root='v_w'
      call read_regional_bc_file(is_input,ie_input,js_input,je_input    &
                                ,nlev                                   &
                                ,ntracers                               &
!                               ,Atm%regional_bc_bounds                 &
                                ,var_name_root                          &
                                ,array_3d=v_w_input)
!
!-----------------------------------------------------------------------
!***  If this is a DA-updated BC file then also read in the layer
!***  pressure depths.
!-----------------------------------------------------------------------
!
      if(Atm%flagstruct%regional_bcs_from_gsi)then
        nlev=klev_in
        var_name_root='delp'
        call read_regional_bc_file(is_input,ie_input,js_input,je_input  &
                                  ,nlev                                 &
                                  ,ntracers                             &
                                  ,var_name_root                        &
                                  ,array_3d=delp_input)
        var_name_root='delz'
        call read_regional_bc_file(is_input,ie_input,js_input,je_input  &
                                  ,nlev                                 &
                                  ,ntracers                             &
                                  ,var_name_root                        &
                                  ,array_3d=delz_input)
      endif
!
!-------------
!***  Tracers
!-------------

      nlev=klev_in
!
!-----------------------------------------------------------------------
!***  Read the tracers specified in the field_table.  If they are not
!***  in the input data then print a warning and set them to 0 in the
!***  boundary. Some tracers are mandatory to have, because they are 
!***  used later for calculating virtual potential temperature etc.
!-----------------------------------------------------------------------
!
      do nt = 1, ntracers 
        call get_tracer_names(MODEL_ATMOS, nt, var_name_root)
        index= get_tracer_index(MODEL_ATMOS,trim(var_name_root))
        if (index==liq_water_index .or. index==sphum_index) then
          required = .true.
        else
          required = .false.
        endif
        call read_regional_bc_file(is_input,ie_input,js_input,je_input  &
                                  ,nlev                                 &
                                  ,ntracers                             &
                                  ,var_name_root                        &
                                  ,array_4d=tracers_input               &
                                  ,tlev=index                           &
                                  ,required=required )
      enddo
!
!-----------------------------------------------------------------------
!***  For a DA-updated BC file we can simply transfer the data 
!***  from the *_input arrays into the model's boundary arrays
!***  since they came out of restart files.  Otherwise proceed
!***  with vertical remapping from input layers to model forecast
!***  layers and rotate the winds from geographic lat/lon to the
!***  integration grid.
!-----------------------------------------------------------------------
!
      data_to_BC: if(Atm%flagstruct%regional_bcs_from_gsi)then             !<-- Fill BC arrays directly from the BC file data
!
!-----------------------------------------------------------------------
!
        call fill_BC_for_DA
!      
!-----------------------------------------------------------------------
!
      else                                                                 !<-- Rotate winds and vertically remap BC file data
!
!-----------------------------------------------------------------------
!***  One final array needs to be allocated.  It is the sfc pressure
!***  in the domain's boundary region that is derived from the input
!***  sfc pressure from the BC files.  The derived sfc pressure will
!***  be needed in the vertical remapping of the wind components to
!***  the integration levels.
!-----------------------------------------------------------------------
!
        allocate(ps_reg(is_input:ie_input,js_input:je_input)) ; ps_reg=-9999999 ! for now don't set to snan until remap dwinds is changed
!
!-----------------------------------------------------------------------
!***  We have the boundary variables from the BC file on the levels
!***  of the input data.  Remap the scalars (tracers, vertical 
!***  velocity, ozone) to the FV3 domain levels.  Scalar remapping
!***  must be done on all four sides before remapping of the winds
!***  since pressures are needed on each side of wind points and so
!***  for a given wind component those pressures could include values
!***  from two different boundary side regions.
!-----------------------------------------------------------------------
!
! Definitions in this module greatly differ from those in existing nesting
!  code or elsewhere in FMS. North <--> South, East <--> West, and 
!  North and South always span  [isd-1 , ied+1] while East and West do not
!  go into the outermost corners (so the they span [1, je], always.)
!-----------------------------------------------------------------------
        sides_scalars: do nside=1,4
!-----------------------------------------------------------------------
!-----------
!***  North
!-----------
!
          call_remap=.false.
!
        if(nside==1)then
          if(north_bc)then
            call_remap=.true.
            side='north'
            bc_side_t1=>BC_t1%north
            bc_side_t0=>BC_t0%north
          endif
        endif
!
        if(nside==2)then
          if(south_bc)then
            call_remap=.true.
            side='south'
            bc_side_t1=>BC_t1%south
            bc_side_t0=>BC_t0%south
          endif
        endif
!
          if(nside==3)then
            if(east_bc)then
              call_remap=.true.
              side='east '
              bc_side_t1=>BC_t1%east
              bc_side_t0=>BC_t0%east  
            endif
          endif
!
        if(nside==4)then
          if(west_bc)then
            call_remap=.true.
            side='west '
            bc_side_t1=>BC_t1%west
            bc_side_t0=>BC_t0%west
          endif
        endif
!
          if(call_remap)then
            call remap_scalar_nggps_regional_bc(Atm                     &
                                               ,side                    &

                                               ,isd,ied,jsd,jed         &  !<-- Atm array indices w/halo

                                               ,is_input                &  !<--
                                               ,ie_input                &  !  Input array
                                               ,js_input                &  !  index limits.
                                               ,je_input                &  !<--

                                               ,klev_in, klev_out       &
                                               ,ntracers                &
                                               ,ak, bk                  &

                                               ,ps_input                &  !<--
                                               ,t_input                 &  !  BC vbls
                                               ,tracers_input           &  !  on input
                                               ,w_input                 &  !  model levels
                                               ,zh_input                &  !<--

                                               ,phis_reg                &  !<-- Filtered topography

                                               ,ps_reg                  &  !<-- Derived FV3 psfc in regional domain boundary region

                                               ,bc_side_t1 )               !<-- BC vbls on final integration levels
!
            call set_delp_and_tracers(bc_side_t1,Atm%npz,Atm%flagstruct%nwat)
!
           if(nside==1)then
              if(north_bc)then
                if (is == 1) then
                   istart = 1
                else
                   istart = isd
                endif
                if (ie == npx-1) then
                     iend = npx-1
                else
                   iend = ied
                endif

                do k=1,npz
                do j=jsd,0
                do i=istart,iend 
                      delz_regBC%south_t1(i,j,k) = bc_side_t1%delz_BC(i,j,k)
                      delz_regBC%south_t0(i,j,k) = bc_side_t0%delz_BC(i,j,k)
                enddo
                enddo
                enddo

          ! North, south include all corners
                if (is == 1) then
                  do k=1,npz
                  do j=jsd,0
                  do i=isd,0
                     delz_regBC%west_t1(i,j,k) = bc_side_t1%delz_BC(i,j,k)
             	     delz_regBC%west_t0(i,j,k) = bc_side_t0%delz_BC(i,j,k)
                  enddo
          	  enddo
          	  enddo
                endif

                if (ie == npx-1) then
                  do k=1,npz
                  do j=jsd,0
                  do i=npx,ied
                     delz_regBC%east_t1(i,j,k) = bc_side_t1%delz_BC(i,j,k)
                     delz_regBC%east_t0(i,j,k) = bc_side_t0%delz_BC(i,j,k)
                  enddo
                  enddo
                  enddo
                endif
              endif 
           endif

           if(nside==2)then
              if(south_bc)then
                if (is == 1) then
                   istart = 1
                else
                   istart = isd
                endif
                if (ie == npx-1) then
                     iend = npx-1
                else
                   iend = ied
                endif

                do k=1,npz
                do j=npy,jed
                do i=istart,iend 
                      delz_regBC%north_t1(i,j,k) = bc_side_t1%delz_BC(i,j,k)
                      delz_regBC%north_t0(i,j,k) = bc_side_t0%delz_BC(i,j,k)
                enddo
                enddo
                enddo

          ! North, south include all corners
                if (is == 1) then
                  do k=1,npz
                  do j=npy,jed
                  do i=isd,0
                     delz_regBC%west_t1(i,j,k) = bc_side_t1%delz_BC(i,j,k)
             	     delz_regBC%west_t0(i,j,k) = bc_side_t0%delz_BC(i,j,k)
                  enddo
          	  enddo
          	  enddo
                endif

 
                if (ie == npx-1) then
                  do k=1,npz
                  do j=npy,jed
                  do i=npx,ied
                     delz_regBC%east_t1(i,j,k) = bc_side_t1%delz_BC(i,j,k)
                     delz_regBC%east_t0(i,j,k) = bc_side_t0%delz_BC(i,j,k)
                  enddo
                  enddo
                  enddo
                endif
              endif 
           endif
      
!          
            
           if(nside==3)then
              if(east_bc)then
                if (js == 1) then
                    jstart = 1
                else
                    jstart = jsd
                endif
                if (je == npy-1) then
                   jend = je
                else
                   jend = jed
                endif


                do k=1,npz
                do j=jstart,jend
                do i=isd,0
                     delz_regBC%west_t1(i,j,k) = bc_side_t1%delz_BC(i,j,k)
                     delz_regBC%west_t0(i,j,k) = bc_side_t0%delz_BC(i,j,k)
                enddo
                enddo
                enddo
              endif
            endif

          if(nside==4)then
              if(west_bc)then
                if (js == 1) then
                    jstart = 1
                else
                    jstart = jsd
                endif
                if (je == npy-1) then
                   jend = je
                else
                   jend = jed
                endif


                do k=1,npz
                do j=jstart,jend
                do i=npx,ied
                     delz_regBC%east_t1(i,j,k) = bc_side_t1%delz_BC(i,j,k)
                     delz_regBC%east_t0(i,j,k) = bc_side_t0%delz_BC(i,j,k)
                enddo
                enddo
                enddo
              endif
            endif

        endif
!
!-----------------------------------------------------------------------
        enddo sides_scalars
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Now that we have the pressure throughout the boundary region
!***  including a row beyond the boundary winds we are ready to
!***  finalize those winds.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Transform the D-grid wind components on the four sides of
!***  the regional domain then remap them from the input levels
!***  to the integration levels.
!-----------------------------------------------------------------------
!
#ifdef USE_FMS_READ
        isc2 = 2*(isd-1+nhalo_data)-1
        iec2 = 2*(ied+2+nhalo_data)-1
        jsc2 = 2*(jsd-1+nhalo_data)-1
        jec2 = 2*(jed+2+nhalo_data)-1
        allocate(tmpx(isc2:iec2, jsc2:jec2)) ; tmpx=dbl_snan
        allocate(tmpy(isc2:iec2, jsc2:jec2)) ; tmpy=dbl_snan
        start = 1; nread = 1
        start(1) = isc2; nread(1) = iec2 - isc2 + 1
        start(2) = jsc2; nread(2) = jec2 - jsc2 + 1
        call read_data("INPUT/grid.tile7.halo4.nc", 'x', tmpx, start, nread, no_domain=.TRUE.)
        call read_data("INPUT/grid.tile7.halo4.nc", 'y', tmpy, start, nread, no_domain=.TRUE.)

        allocate(reg_grid(isd-1:ied+2,jsd-1:jed+2,1:2)) ; reg_grid=dbl_snan
        do j = jsd-1, jed+2
        do i = isd-1, ied+2
            reg_grid(i,j,1) = tmpx(2*(i+nhalo_data)-1, 2*(j+nhalo_data)-1)*pi/180.
            reg_grid(i,j,2) = tmpy(2*(i+nhalo_data)-1, 2*(j+nhalo_data)-1)*pi/180.
            if ( reg_grid(i,j,1) /= grid_reg(i,j,1) ) then
               write(0,*)' reg_grid(i,j,1) /= grid_reg(i,j,1) ',i,j, reg_grid(i,j,1),grid_reg(i,j,1)
            endif
        enddo
        enddo

        allocate(reg_agrid(isd-1:ied+1,jsd-1:jed+1,1:2)) ; reg_agrid=dbl_snan
        do j=jsd-1,jed+1
        do i=isd-1,ied+1
            call cell_center2(reg_grid(i,j,  1:2), reg_grid(i+1,j,  1:2),   &
                              reg_grid(i,j+1,1:2), reg_grid(i+1,j+1,1:2),   &
                              reg_agrid(i,j,1:2) )
        enddo
        enddo
#endif
!
!-----------------------------------------------------------------------
!***  Loop through the four sides of the domain.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
        sides_winds: do nside=1,4
!-----------------------------------------------------------------------
!
          call_remap=.false.

          if(nside==1)then
            if(north_bc)then
              call_remap=.true.
              bc_side_t1=>BC_t1%north
!
              is_u=Atm%regional_bc_bounds%is_north_uvs
              ie_u=Atm%regional_bc_bounds%ie_north_uvs
              js_u=Atm%regional_bc_bounds%js_north_uvs
              je_u=Atm%regional_bc_bounds%je_north_uvs
!
              is_v=Atm%regional_bc_bounds%is_north_uvw
              ie_v=Atm%regional_bc_bounds%ie_north_uvw
              js_v=Atm%regional_bc_bounds%js_north_uvw
              je_v=Atm%regional_bc_bounds%je_north_uvw
            endif
          endif
!
          if(nside==2)then
            if(south_bc)then
              call_remap=.true.
              bc_side_t1=>BC_t1%south
!
              is_u=Atm%regional_bc_bounds%is_south_uvs
              ie_u=Atm%regional_bc_bounds%ie_south_uvs
              js_u=Atm%regional_bc_bounds%js_south_uvs
              je_u=Atm%regional_bc_bounds%je_south_uvs
!
              is_v=Atm%regional_bc_bounds%is_south_uvw
              ie_v=Atm%regional_bc_bounds%ie_south_uvw
              js_v=Atm%regional_bc_bounds%js_south_uvw
              je_v=Atm%regional_bc_bounds%je_south_uvw
            endif
          endif
!
          if(nside==3)then
            if(east_bc)then
              call_remap=.true.
              bc_side_t1=>BC_t1%east
!
              is_u=Atm%regional_bc_bounds%is_east_uvs
              ie_u=Atm%regional_bc_bounds%ie_east_uvs
              js_u=Atm%regional_bc_bounds%js_east_uvs
              je_u=Atm%regional_bc_bounds%je_east_uvs
!
              is_v=Atm%regional_bc_bounds%is_east_uvw
              ie_v=Atm%regional_bc_bounds%ie_east_uvw
              js_v=Atm%regional_bc_bounds%js_east_uvw
              je_v=Atm%regional_bc_bounds%je_east_uvw
            endif
          endif
!
          if(nside==4)then
            if(west_bc)then
              call_remap=.true.
              bc_side_t1=>BC_t1%west
!
              is_u=Atm%regional_bc_bounds%is_west_uvs
              ie_u=Atm%regional_bc_bounds%ie_west_uvs
              js_u=Atm%regional_bc_bounds%js_west_uvs
              je_u=Atm%regional_bc_bounds%je_west_uvs
!
              is_v=Atm%regional_bc_bounds%is_west_uvw
              ie_v=Atm%regional_bc_bounds%ie_west_uvw
              js_v=Atm%regional_bc_bounds%js_west_uvw
              je_v=Atm%regional_bc_bounds%je_west_uvw
            endif
          endif
!
          if(call_remap)then
!
            allocate(ud(is_u:ie_u,js_u:je_u,1:nlev)) ; ud=real_snan
            allocate(vd(is_v:ie_v,js_v:je_v,1:nlev)) ; vd=real_snan
            allocate(vc(is_u:ie_u,js_u:je_u,1:nlev)) ; vc=real_snan
            allocate(uc(is_v:ie_v,js_v:je_v,1:nlev)) ; uc=real_snan
!
            do k=1,nlev
              do j=js_u,je_u
              do i=is_u,ie_u
                p1(:) = grid_reg(i,  j,1:2)
                p2(:) = grid_reg(i+1,j,1:2)
                call  mid_pt_sphere(p1, p2, p3)
                call get_unit_vect2(p1, p2, e1)
                call get_latlon_vector(p3, ex, ey)
                ud(i,j,k) = u_s_input(i,j,k)*inner_prod(e1,ex)+v_s_input(i,j,k)*inner_prod(e1,ey)
                p4(:) = agrid_reg(i,j,1:2) ! cell centroid
                call get_unit_vect2(p3, p4, e2) !C-grid V-wind unit vector
                vc(i,j,k) = u_s_input(i,j,k)*inner_prod(e2,ex)+v_s_input(i,j,k)*inner_prod(e2,ey)
              enddo
              enddo
!
              do j=js_v,je_v
                do i=is_v,ie_v
                  p1(:) = grid_reg(i,j  ,1:2)
                  p2(:) = grid_reg(i,j+1,1:2)
                  call  mid_pt_sphere(p1, p2, p3)
                  call get_unit_vect2(p1, p2, e2)
                  call get_latlon_vector(p3, ex, ey)
                  vd(i,j,k) = u_w_input(i,j,k)*inner_prod(e2,ex)+v_w_input(i,j,k)*inner_prod(e2,ey)
                  p4(:) = agrid_reg(i,j,1:2) ! cell centroid
                  call get_unit_vect2(p3, p4, e1) !C-grid U-wind unit vector
                  uc(i,j,k) = u_w_input(i,j,k)*inner_prod(e1,ex)+v_w_input(i,j,k)*inner_prod(e1,ey)
                enddo
              enddo
            enddo
!
            call remap_dwinds_regional_bc(Atm                           &

                                         ,is_input                      &  !<--
                                         ,ie_input                      &  !  Index limits for scalars
                                         ,js_input                      &  !  at center of north BC region grid cells.
                                         ,je_input                      &  !<--

                                         ,is_u                          &  !<--
                                         ,ie_u                          &  !  Index limits for u component
                                         ,js_u                          &  !  on north edge of BC region grid cells.
                                         ,je_u                          &  !<--

                                         ,is_v                          &  !<--
                                         ,ie_v                          &  !  Index limits for v component
                                         ,js_v                          &  !  on north edge of BC region grid cells.
                                         ,je_v                          &  !<--

                                         ,klev_in, klev_out             &  !<-- data / model levels
                                         ,ak, bk                        &

                                         ,ps_reg                        &  !<-- BC values of sfc pressure
                                         ,ud ,vd                        &  !<-- BC values of D-grid u and v
                                         ,uc ,vc                        &  !<-- BC values of C-grid u and v
                                         ,bc_side_t1 )                     !<-- North BC vbls on final integration levels
!
            deallocate(ud,vd,uc,vc)
!
          endif
!
!-----------------------------------------------------------------------
        enddo sides_winds
!-----------------------------------------------------------------------
!
      endif data_to_BC
!
!-----------------------------------------------------------------------
!***  Close the boundary file.
!-----------------------------------------------------------------------
!
      call check(nf90_close(ncid))
!
!-----------------------------------------------------------------------
!***  Deallocate working arrays.
!-----------------------------------------------------------------------
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(allocated(ps_input))then
        deallocate(ps_input)
      endif
      if(allocated(t_input))then
        deallocate(t_input)
      endif
      if(allocated(zh_input))then
        deallocate(zh_input)
      endif
      if(allocated(w_input))then
        deallocate(w_input)
      endif
      if(allocated(tracers_input))then
        deallocate(tracers_input)
      endif
      if(allocated(u_s_input))then
        deallocate(u_s_input)
      endif
      if(allocated(u_w_input))then
        deallocate(u_w_input)
      endif
      if(allocated(v_s_input))then
        deallocate(v_s_input)
      endif
      if(allocated(v_w_input))then
        deallocate(v_w_input)
      endif
      if(allocated(delp_input))then
        deallocate(delp_input)
      endif
      if(allocated(delz_input))then
        deallocate(delz_input)
      endif
!
!-----------------------------------------------------------------------
!***  Fill the remaining boundary arrays starting with the divergence.
!-----------------------------------------------------------------------
!
      call fill_divgd_BC
!
!-----------------------------------------------------------------------
!***  Fill the total condensate in the regional boundary array.
!-----------------------------------------------------------------------
!
#ifdef USE_COND
      call fill_q_con_BC
#endif
!
!-----------------------------------------------------------------------
!***  Fill moist kappa in the regional domain boundary array.
!-----------------------------------------------------------------------
!
#ifdef MOIST_CAPPA
      call fill_cappa_BC
#endif
!
!-----------------------------------------------------------------------
!***  Convert the boundary region sensible temperature array to 
!***  FV3's modified virtual potential temperature.
!-----------------------------------------------------------------------
!
      call convert_to_virt_pot_temp(isd,ied,jsd,jed,npz)
!
!-----------------------------------------------------------------------
!***  If nudging of the specific humidity has been selected then
!***  nudge the boundary values in the same way as is done for the
!***  interior.
!-----------------------------------------------------------------------
!
      if(Atm%flagstruct%nudge_qv)then
        call nudge_qv_bc(Atm,isd,ied,jsd,jed)
      endif
!
!-----------------------------------------------------------------------

      contains

!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      subroutine fill_BC_for_DA
!
!-----------------------------------------------------------------------
!***  Transfer the input boundary data directly into the BC object.
!-----------------------------------------------------------------------
!
!---------------------
!***  Local variables
!---------------------
!
      integer :: i,j,k,n
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Since corner tasks are on more than one side we cannot 
!***  generalize the transfer of data into a given side's
!***  arrays.  Do each side separately.
!
!***  Simply obtain the loop limits from the bounds of one of the
!***  BC arrays to be filled.
!-----------------------------------------------------------------------
!
!-----------
!***  North
!-----------
!
      if(north_bc)then
!
        is_input=lbound(BC_t1%north%delp_BC,1)                             !<-- 
        ie_input=ubound(BC_t1%north%delp_BC,1)                             !  Index limits for
        js_input=lbound(BC_t1%north%delp_BC,2)                             !  mass variables.
        je_input=ubound(BC_t1%north%delp_BC,2)                             !<--
!
        do k=1,klev_in
        do j=js_input,je_input
        do i=is_input,ie_input
          BC_t1%north%delp_BC(i,j,k)=delp_input(i,j,k)
          BC_t1%north%pt_BC(i,j,k)=t_input(i,j,k)
          BC_t1%north%w_BC(i,j,k)=w_input(i,j,k)
          BC_t1%north%delz_BC(i,j,k)=delz_input(i,j,k)
        enddo
        enddo
        enddo
!
        do n=1,ntracers
        do k=1,klev_in
        do j=js_input,je_input
        do i=is_input,ie_input
          BC_t1%north%q_BC(i,j,k,n)=tracers_input(i,j,k,n)
        enddo
        enddo
        enddo
        enddo
!
        is_input=lbound(BC_t1%north%u_BC,1)                                !<-- 
        ie_input=ubound(BC_t1%north%u_BC,1)                                !  Index limits for
        js_input=lbound(BC_t1%north%u_BC,2)                                !  D-grid u and C-grid v.
        je_input=ubound(BC_t1%north%u_BC,2)                                !<--
!
        do k=1,klev_in
        do j=js_input,je_input
        do i=is_input,ie_input
          BC_t1%north%u_BC(i,j,k)=u_s_input(i,j,k)
          BC_t1%north%vc_BC(i,j,k)=v_s_input(i,j,k)
        enddo
        enddo
        enddo
!
        is_input=lbound(BC_t1%north%v_BC,1)                                !<-- 
        ie_input=ubound(BC_t1%north%v_BC,1)                                !  Index limits for
        js_input=lbound(BC_t1%north%v_BC,2)                                !  D-grid v and C-grid u.
        je_input=ubound(BC_t1%north%v_BC,2)                                !<--
!
        do k=1,klev_in
        do j=js_input,je_input
        do i=is_input,ie_input
          BC_t1%north%v_BC(i,j,k)=v_w_input(i,j,k)
          BC_t1%north%uc_BC(i,j,k)=u_w_input(i,j,k)
        enddo
        enddo
        enddo
!
      endif
!
!-----------
!***  South
!-----------
!
      if(south_bc)then
        is_input=lbound(BC_t1%south%delp_BC,1)                             !<---
        ie_input=ubound(BC_t1%south%delp_BC,1)                             !   Index limits for
        js_input=lbound(BC_t1%south%delp_BC,2)                             !   mass variables.
        je_input=ubound(BC_t1%south%delp_BC,2)                             !<--
!
        do k=1,klev_in
        do j=js_input,je_input
        do i=is_input,ie_input
          BC_t1%south%delp_BC(i,j,k)=delp_input(i,j,k)
          BC_t1%south%pt_BC(i,j,k)=t_input(i,j,k)
          BC_t1%south%w_BC(i,j,k)=w_input(i,j,k)
          BC_t1%south%delz_BC(i,j,k)=delz_input(i,j,k)
        enddo
        enddo
        enddo
!
        do n=1,ntracers
        do k=1,klev_in
        do j=js_input,je_input
        do i=is_input,ie_input
          BC_t1%south%q_BC(i,j,k,n)=tracers_input(i,j,k,n)
        enddo
        enddo
        enddo
        enddo
!
        is_input=lbound(BC_t1%south%u_BC,1)                                !<-- 
        ie_input=ubound(BC_t1%south%u_BC,1)                                !  Index limits for
        js_input=lbound(BC_t1%south%u_BC,2)                                !  D-grid u and C-grid v.
        je_input=ubound(BC_t1%south%u_BC,2)                                !<--
!
        do k=1,klev_in
        do j=js_input,je_input
        do i=is_input,ie_input
          BC_t1%south%u_BC(i,j,k)=u_s_input(i,j,k)
          BC_t1%south%vc_BC(i,j,k)=v_s_input(i,j,k)
        enddo
        enddo
        enddo
!
        is_input=lbound(BC_t1%south%v_BC,1)                                !<-- 
        ie_input=ubound(BC_t1%south%v_BC,1)                                !  Index limits for
        js_input=lbound(BC_t1%south%v_BC,2)                                !  D-grid v and C-grid u.
        je_input=ubound(BC_t1%south%v_BC,2)                                !<--
!
        do k=1,klev_in
        do j=js_input,je_input
        do i=is_input,ie_input
          BC_t1%south%v_BC(i,j,k)=v_w_input(i,j,k)
          BC_t1%south%uc_BC(i,j,k)=u_w_input(i,j,k)
        enddo
        enddo
        enddo
!
      endif
!
!----------
!***  East
!----------
!
      if(east_bc)then
        is_input=lbound(BC_t1%east%delp_BC,1)                               !<--
        ie_input=ubound(BC_t1%east%delp_BC,1)                               !  Index limits
        js_input=lbound(BC_t1%east%delp_BC,2)                               !  for mass variables.
        je_input=ubound(BC_t1%east%delp_BC,2)                               !<--
!
        do k=1,klev_in
        do j=js_input,je_input
        do i=is_input,ie_input
          BC_t1%east%delp_BC(i,j,k)=delp_input(i,j,k)
          BC_t1%east%pt_BC(i,j,k)=t_input(i,j,k)
          BC_t1%east%w_BC(i,j,k)=w_input(i,j,k)
          BC_t1%east%delz_BC(i,j,k)=delz_input(i,j,k)
        enddo
        enddo
        enddo
!
        do n=1,ntracers
        do k=1,klev_in
        do j=js_input,je_input
        do i=is_input,ie_input
          BC_t1%east%q_BC(i,j,k,n)=tracers_input(i,j,k,n)
        enddo
        enddo
        enddo
        enddo
!
        is_input=lbound(BC_t1%east%u_BC,1)                                 !<-- 
        ie_input=ubound(BC_t1%east%u_BC,1)                                 !  Index limits for
        js_input=lbound(BC_t1%east%u_BC,2)                                 !  D-grid u and C-grid v.
        je_input=ubound(BC_t1%east%u_BC,2)                                 !<--
!
        do k=1,klev_in
        do j=js_input,je_input
        do i=is_input,ie_input
          BC_t1%east%u_BC(i,j,k)=u_s_input(i,j,k)
          BC_t1%east%vc_BC(i,j,k)=v_s_input(i,j,k)
        enddo
        enddo
        enddo
!
        is_input=lbound(BC_t1%east%v_BC,1)                                 !<-- 
        ie_input=ubound(BC_t1%east%v_BC,1)                                 !  Index limits for
        js_input=lbound(BC_t1%east%v_BC,2)                                 !  D-grid v and C-grid u.
        je_input=ubound(BC_t1%east%v_BC,2)                                 !<--
!
        do k=1,klev_in
        do j=js_input,je_input
        do i=is_input,ie_input
          BC_t1%east%v_BC(i,j,k)=v_w_input(i,j,k)
          BC_t1%east%uc_BC(i,j,k)=u_w_input(i,j,k)
        enddo
        enddo
        enddo
!
      endif
!
!----------
!***  West
!----------
!
      if(west_bc)then
        is_input=lbound(BC_t1%west%delp_BC,1)                              !<--
        ie_input=ubound(BC_t1%west%delp_BC,1)                              !  Index limits for
        js_input=lbound(BC_t1%west%delp_BC,2)                              !  mass variables.
        je_input=ubound(BC_t1%west%delp_BC,2)                              !<--
!
        do k=1,klev_in
        do j=js_input,je_input
        do i=is_input,ie_input
          BC_t1%west%delp_BC(i,j,k)=delp_input(i,j,k)
          BC_t1%west%pt_BC(i,j,k)=t_input(i,j,k)
          BC_t1%west%w_BC(i,j,k)=w_input(i,j,k)
          BC_t1%west%delz_BC(i,j,k)=delz_input(i,j,k)
        enddo
        enddo
        enddo
!
        do n=1,ntracers
        do k=1,klev_in
        do j=js_input,je_input
        do i=is_input,ie_input
          BC_t1%west%q_BC(i,j,k,n)=tracers_input(i,j,k,n)
        enddo
        enddo
        enddo
        enddo
!
        is_input=lbound(BC_t1%west%u_BC,1)                                 !<-- 
        ie_input=ubound(BC_t1%west%u_BC,1)                                 !  Index limits for
        js_input=lbound(BC_t1%west%u_BC,2)                                 !  D-grid u and C-grid v.
        je_input=ubound(BC_t1%west%u_BC,2)                                 !<--
!
        do k=1,klev_in
        do j=js_input,je_input
        do i=is_input,ie_input
          BC_t1%west%u_BC(i,j,k)=u_s_input(i,j,k)
          BC_t1%west%vc_BC(i,j,k)=v_s_input(i,j,k)
        enddo
        enddo
        enddo
!
        is_input=lbound(BC_t1%west%v_BC,1)                                 !<-- 
        ie_input=ubound(BC_t1%west%v_BC,1)                                 !  Index limits for
        js_input=lbound(BC_t1%west%v_BC,2)                                 !  D-grid v and C-grid u.
        je_input=ubound(BC_t1%west%v_BC,2)                                 !<--
!
        do k=1,klev_in
        do j=js_input,je_input
        do i=is_input,ie_input
          BC_t1%west%v_BC(i,j,k)=v_w_input(i,j,k)
          BC_t1%west%uc_BC(i,j,k)=u_w_input(i,j,k)
        enddo
        enddo
        enddo
!
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine fill_BC_for_DA
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      subroutine fill_divgd_BC
!
!-----------------------------------------------------------------------
!***  For now fill the boundary divergence with zero.
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!
!--------------------
!***  Local variables
!--------------------
!
      integer :: i,ie0,is0,j,je0,js0,k,nside
!
      logical :: call_set
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Loop through the four sides.
!-----------------------------------------------------------------------
!
      do nside=1,4
!
        call_set=.false.
!
        if(nside==1)then
          if(north_bc)then
            call_set=.true.
            bc_side_t1=>BC_t1%north
            is0=lbound(BC_t1%north%divgd_BC,1)
            ie0=ubound(BC_t1%north%divgd_BC,1)
            js0=lbound(BC_t1%north%divgd_BC,2)
            je0=ubound(BC_t1%north%divgd_BC,2)
          endif
        endif
!
        if(nside==2)then
          if(south_bc)then
            call_set=.true.
            bc_side_t1=>BC_t1%south
            is0=lbound(BC_t1%south%divgd_BC,1)
            ie0=ubound(BC_t1%south%divgd_BC,1)
            js0=lbound(BC_t1%south%divgd_BC,2)
            je0=ubound(BC_t1%south%divgd_BC,2)
          endif
        endif
!
        if(nside==3)then
          if(east_bc)then
            call_set=.true.
            bc_side_t1=>BC_t1%east
            is0=lbound(BC_t1%east%divgd_BC,1)
            ie0=ubound(BC_t1%east%divgd_BC,1)
            js0=lbound(BC_t1%east%divgd_BC,2)
            je0=ubound(BC_t1%east%divgd_BC,2)
          endif
        endif
!
        if(nside==4)then
          if(west_bc)then
            call_set=.true.
            bc_side_t1=>BC_t1%west
            is0=lbound(BC_t1%west%divgd_BC,1)
            ie0=ubound(BC_t1%west%divgd_BC,1)
            js0=lbound(BC_t1%west%divgd_BC,2)
            je0=ubound(BC_t1%west%divgd_BC,2)
          endif
        endif
!
        if(call_set)then
          do k=1,klev_out
            do j=js0,je0
            do i=is0,ie0
              bc_side_t1%divgd_BC(i,j,k)=0.
            enddo
            enddo
          enddo
        endif
!
      enddo
!
!-----------------------------------------------------------------------
!
      end subroutine fill_divgd_BC
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
#ifdef USE_COND
      subroutine fill_q_con_BC
!
!-----------------------------------------------------------------------
!***  For now fill the total condensate in the boundary regiona
!***  with only the liquid water content.
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!
!--------------------
!***  Local variables
!--------------------
!
      integer :: i,ie0,is0,j,je0,js0,k,nside
!
      logical :: call_set
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Loop through the four sides.
!-----------------------------------------------------------------------
!
      do nside=1,4
        call_set=.false.
!
        if(nside==1)then
          if(north_bc)then
            call_set=.true.
            bc_side_t1=>BC_t1%north
            is0=lbound(BC_t1%north%q_con_BC,1)
            ie0=ubound(BC_t1%north%q_con_BC,1)
            js0=lbound(BC_t1%north%q_con_BC,2)
            je0=ubound(BC_t1%north%q_con_BC,2)
          endif
        endif
!
        if(nside==2)then
          if(south_bc)then
            call_set=.true.
            bc_side_t1=>BC_t1%south
            is0=lbound(BC_t1%south%q_con_BC,1)
            ie0=ubound(BC_t1%south%q_con_BC,1)
            js0=lbound(BC_t1%south%q_con_BC,2)
            je0=ubound(BC_t1%south%q_con_BC,2)
          endif
        endif
!
        if(nside==3)then
          if(east_bc)then
            call_set=.true.
            bc_side_t1=>BC_t1%east
            is0=lbound(BC_t1%east%q_con_BC,1)
            ie0=ubound(BC_t1%east%q_con_BC,1)
            js0=lbound(BC_t1%east%q_con_BC,2)
            je0=ubound(BC_t1%east%q_con_BC,2)
          endif
        endif
!
        if(nside==4)then
          if(west_bc)then
            call_set=.true.
            bc_side_t1=>BC_t1%west
            is0=lbound(BC_t1%west%q_con_BC,1)
            ie0=ubound(BC_t1%west%q_con_BC,1)
            js0=lbound(BC_t1%west%q_con_BC,2)
            je0=ubound(BC_t1%west%q_con_BC,2)
          endif
        endif
!
        if(call_set)then
          do k=1,klev_out
            do j=js0,je0
            do i=is0,ie0
              bc_side_t1%q_con_BC(i,j,k)=bc_side_t1%q_BC(i,j,k,liq_water_index)
            enddo
            enddo
          enddo
        endif
!
      enddo
!
!-----------------------------------------------------------------------
!
      end subroutine fill_q_con_BC
#endif
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
#ifdef MOIST_CAPPA
      subroutine fill_cappa_BC
!
!-----------------------------------------------------------------------
!***  Compute cappa in the regional domain boundary area following
!***  Zhao-Carr microphysics.
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!
!---------------------
!***  Local variables
!---------------------
!
      integer :: i1,i2,j1,j2,nside
!
      real,dimension(:,:,:),pointer :: cappa,temp,liq_wat,sphum
!
      logical :: call_compute
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      do nside=1,4
        call_compute=.false.
!
        if(nside==1)then
          if(north_bc)then
            call_compute=.true.
            bc_side_t1=>BC_t1%north
          endif
        endif
!
        if(nside==2)then
          if(south_bc)then
            call_compute=.true.
            bc_side_t1=>BC_t1%south
          endif
        endif
!
        if(nside==3)then
          if(east_bc)then
            call_compute=.true.
            bc_side_t1=>BC_t1%east
          endif
        endif
!
        if(nside==4)then
          if(west_bc)then
            call_compute=.true.
            bc_side_t1=>BC_t1%west
          endif
        endif
!
        if(call_compute)then
          i1=lbound(bc_side_t1%cappa_BC,1)
          i2=ubound(bc_side_t1%cappa_BC,1)
          j1=lbound(bc_side_t1%cappa_BC,2)
          j2=ubound(bc_side_t1%cappa_BC,2)
          cappa  =>bc_side_t1%cappa_BC
          temp   =>bc_side_t1%pt_BC
          liq_wat=>bc_side_t1%q_BC(:,:,:,liq_water_index)
          sphum  =>bc_side_t1%q_BC(:,:,:,sphum_index)
          call compute_cappa(i1,i2,j1,j2,cappa,temp,liq_wat,sphum)
        endif
!
      enddo
!
!-----------------------------------------------------------------------
!
      end subroutine fill_cappa_BC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      subroutine compute_cappa(i1,i2,j1,j2,cappa,temp,liq_wat,sphum)
!
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      integer,intent(in) :: i1,i2,j1,j2
!
      real,dimension(i1:i2,j1:j2,1:npz),intent(in) :: temp,liq_wat,sphum
      real,dimension(i1:i2,j1:j2,1:npz),intent(inout) :: cappa
!
!---------------------
!***  Local variables
!---------------------
!
      integer :: i,ie,is,j,je,js,k
!
      real :: cvm,qd,ql,qs,qv
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      is=lbound(cappa,1)
      ie=ubound(cappa,1)
      js=lbound(cappa,2)
      je=ubound(cappa,2)
!
      do k=1,klev_out
        do j=js,je
        do i=is,ie
          qd=max(0.,liq_wat(i,j,k))
          if( temp(i,j,k) > tice )then
            qs=0.
          elseif( temp(i,j,k) < tice-t_i0 )then
            qs=qd
          else
            qs=qd*(tice-temp(i,j,k))/t_i0
          endif
          ql=qd-qs
          qv=max(0.,sphum(i,j,k))
          cvm=(1.-(qv+qd))*cv_air + qv*cv_vap + ql*c_liq + qs*c_ice
!
          cappa(i,j,k)=rdgas/(rdgas+cvm/(1.+zvir*sphum(i,j,k)))
!
        enddo
        enddo
      enddo
!
!-----------------------------------------------------------------------
!
      end subroutine compute_cappa
#endif
!
!-----------------------------------------------------------------------
!
      end subroutine regional_bc_data

!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      subroutine read_regional_bc_file(is_input,ie_input                &
                                      ,js_input,je_input                &
                                      ,nlev                             &
                                      ,ntracers                         &
                                      ,var_name_root                    &
                                      ,array_3d                         &
                                      ,array_4d                         &
                                      ,tlev                             &
                                      ,required )
!-----------------------------------------------------------------------
!***  Read the boundary data from the external file generated by
!***  chgres.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
!----------
!*** Input
!----------
!
      integer,intent(in) :: is_input,ie_input,js_input,je_input,nlev
      integer,intent(in) :: ntracers
!
      integer,intent(in),optional :: tlev                                  !<-- Position of current tracer among all of them
!
      character(len=*),intent(in) :: var_name_root                         !<-- Root of variable name in the boundary file
      logical,intent(in),optional :: required
!
!------------
!***  Output
!------------
!
      real,dimension(is_input:ie_input,js_input:je_input,1:nlev),intent(out),optional :: array_3d  !<-- The input 3-D variable's coverage of task subdomain
!
      real,dimension(is_input:ie_input,js_input:je_input,1:nlev,1:ntracers),intent(out),optional :: array_4d  !<-- The input 4-D variable's coverage of subdomain
!
!---------------------
!***  Local variables
!---------------------
!
      integer :: halo,lat,lev,lon
!
      integer :: i_count,i_start_array,i_start_data,i_end_array         &
                ,j_count,j_start_array,j_start_data,j_end_array
!
      integer :: dim_id,nctype,ndims,var_id
      integer :: nside,status
!
      character(len=5) :: dim_name_x                                    &  !<-- Dimension names in
                         ,dim_name_y                                       !    the BC file
!
      character(len=80) :: var_name                                        !<-- Variable name in the boundary NetCDF file
!
      logical :: call_get_var,is_root_pe
      logical :: required_local
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Process optional argument required, default value is .true.
!-----------------------------------------------------------------------
!
      if(present(required)) then
        required_local=required
      else
        required_local=.true.
      endif
!
      is_root_pe=(mpp_pe()==mpp_root_pe())
!
!-----------------------------------------------------------------------
!***  Loop through the four sides of the domain.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      sides: do nside=1,4
!-----------------------------------------------------------------------
!
        call_get_var=.false.
!
!-----------------------------------------------------------------------
!***  Construct the variable's name in the NetCDF file and set
!***  the start locations and point counts for the data file and
!***  for the BC arrays being filled.  The input array begins
!***  receiving data at (i_start_array,j_start_array), etc.
!***  The read of the data for the given input array begins at
!***  (i_start_data,j_start_data) and encompasses i_count by
!***  j_count datapoints in each direction.
!-----------------------------------------------------------------------
!
!-----------
!***  North
!-----------
!
        if(nside==1)then
          if(north_bc)then
            call_get_var=.true.
!
            var_name=trim(var_name_root)//"_bottom"
!
            i_start_array=is_input
            i_end_array  =ie_input
            j_start_array=js_input
            if(trim(var_name_root)=='u_s'.or.trim(var_name_root)=='v_s')then
              j_end_array=js_input+nhalo_data+nrows_blend
            else
              j_end_array=js_input+nhalo_data+nrows_blend-1
            endif
!
            i_start_data=i_start_array+nhalo_data
            i_count=i_end_array-i_start_array+1
            j_start_data=1
            j_count=j_end_array-j_start_array+1
!
          endif
        endif
!
!-----------
!***  South
!-----------
!
        if(nside==2)then
          if(south_bc)then
            call_get_var=.true.
!
            var_name=trim(var_name_root)//"_top"
!
            i_start_array=is_input
            i_end_array  =ie_input
            j_start_array=je_input-nhalo_data-nrows_blend+1
            j_end_array  =je_input
!
            i_start_data=i_start_array+nhalo_data
            i_count=i_end_array-i_start_array+1
            j_start_data=1
            j_count=j_end_array-j_start_array+1
!
          endif
        endif
!
!----------
!***  East
!----------
!
        if(nside==3)then
          if(east_bc)then
            call_get_var=.true.
!
            var_name=trim(var_name_root)//"_left"
!
            j_start_array=js_input
            j_end_array  =je_input
!
            i_start_array=is_input
!
            if(trim(var_name_root)=='u_w'.or.trim(var_name_root)=='v_w')then
              i_end_array=is_input+nhalo_data+nrows_blend
            else
              i_end_array=is_input+nhalo_data+nrows_blend-1
            endif
!
            if(north_bc)then
              if(trim(var_name_root)=='u_s'.or.trim(var_name_root)=='v_s')then
                j_start_array=js_input+nhalo_data+1
              else
                j_start_array=js_input+nhalo_data
              endif
            endif
            if(south_bc)then
              j_end_array  =je_input-nhalo_data
            endif
!
            i_start_data=1
            i_count=i_end_array-i_start_array+1
            if(trim(var_name_root)=='u_s'.or.trim(var_name_root)=='v_s')then
              j_start_data=j_start_array-1
            else
              j_start_data=j_start_array
            endif
            j_count=j_end_array-j_start_array+1
!
          endif
        endif
!
!----------
!***  West
!----------
!
        if(nside==4)then
          if(west_bc)then
            call_get_var=.true.
!
            var_name=trim(var_name_root)//"_right"
!
            j_start_array=js_input
            j_end_array  =je_input
!
            i_start_array=ie_input-nhalo_data-nrows_blend+1
            i_end_array=ie_input
!
            if(north_bc)then
              if(trim(var_name_root)=='u_s'.or.trim(var_name_root)=='v_s')then
                j_start_array=js_input+nhalo_data+1
              else
                j_start_array=js_input+nhalo_data
              endif
            endif
!
            if(south_bc)then
              j_end_array  =je_input-nhalo_data
            endif
!
            i_start_data=1
            i_count=i_end_array-i_start_array+1
            if(trim(var_name_root)=='u_s'.or.trim(var_name_root)=='v_s')then
              j_start_data=j_start_array-1
            else
              j_start_data=j_start_array
            endif
            j_count=j_end_array-j_start_array+1
!
          endif
        endif
!
!-----------------------------------------------------------------------
!***  Fill this task's subset of boundary data for this 3-D
!***  or 4-D variable.  This includes the data in the domain's
!***  halo region as well as the blending region that overlaps
!***  the outer nhalo_blend rows of the integration domain.
!***  If the variable is a tracer then check if it is present
!***  in the input data.  If it is not then print a warning
!***  and set it to zero.
!-----------------------------------------------------------------------
!
        if(call_get_var)then
          if (present(array_4d)) then   !<-- 4-D variable
            status=nf90_inq_varid(ncid,trim(var_name),var_id)                !<-- Get this variable's ID.
            if (required_local) then
              call check(status)
            endif
            if (status /= nf90_noerr) then
              if (east_bc.and.is_master()) write(0,*)' WARNING: Tracer ',trim(var_name),' not in input file'
              array_4d(:,:,:,tlev)=0.                                        !<-- Tracer not in input so set to zero in boundary.
!
              blend_this_tracer(tlev)=.false.                                !<-- Tracer not in input so do not apply blending.
!
            else
              call check(nf90_get_var(ncid,var_id                                         &
                                     ,array_4d(i_start_array:i_end_array                  &  !<-- Fill this task's domain boundary halo.
                                              ,j_start_array:j_end_array                  &
                                              ,1:nlev, tlev)                              &
                                              ,start=(/i_start_data,j_start_data,1,tlev/) &  !<-- Start reading the data array here.
                                              ,count=(/i_count,j_count,nlev,1/)))            !<-- Extent of data to read in each dimension.
!
            endif
!
          else                         !<-- 3-D variable
            call check(nf90_inq_varid(ncid,trim(var_name),var_id))                    !<-- Get this variable's ID.
            call check(nf90_get_var(ncid,var_id                                    &
                                   ,array_3d(i_start_array:i_end_array             &  !<-- Fill this task's domain boundary halo.
                                            ,j_start_array:j_end_array             &
                                            ,1:nlev)                               &
                                            ,start=(/i_start_data,j_start_data,1/) &  !<-- Start reading the data array here.
                                            ,count=(/i_count,j_count,nlev/)))         !<-- Extent of data to read in each dimension.
!
          endif
        endif
!
      enddo sides
!
!-----------------------------------------------------------------------
!
      end subroutine read_regional_bc_file
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      subroutine check(status)
      integer,intent(in) :: status
!
      if(status /= nf90_noerr) then
        write(0,*)' check netcdf status=',status
        call mpp_error(FATAL, ' NetCDF error ' // trim(nf90_strerror(status)))
      endif
!
      end subroutine check
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      subroutine allocate_regional_BC_arrays(side                       &
                                            ,north_bc,south_bc          &
                                            ,east_bc,west_bc            &
                                            ,is_0,ie_0,js_0,je_0        &
                                            ,is_sn,ie_sn,js_sn,je_sn    &
                                            ,is_we,ie_we,js_we,je_we    &
                                            ,klev                       &
                                            ,ntracers                   &
                                            ,BC_side                    &
                                            ,delz_side )
!
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      integer,intent(in) :: klev,ntracers
!
      integer,intent(in) :: is_0,ie_0,js_0,je_0                          !<-- Start/end BC indices for cell centers
      integer,intent(in) :: is_sn,ie_sn,js_sn,je_sn                      !<-- Start/end BC indices for south/north cell edges
      integer,intent(in) :: is_we,ie_we,js_we,je_we                      !<-- Start/end BC indices for west/east cell edges
!
      character(len=5),intent(in) :: side                                !<-- Which side are we allocating?
!
      logical,intent(in) :: north_bc,south_bc,east_bc,west_bc            !<-- Which sides is this task on?
!
      type(fv_regional_BC_variables),intent(out) :: BC_side
!
      real,dimension(:,:,:),pointer,intent(inout),optional :: delz_side  !<-- Boundary delz that follows integration through time.
!
!---------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------
!
      if(allocated(BC_side%delp_BC))then
        return                                                           !<-- The BC arrays are already allocated so exit.
      endif
!
      allocate(BC_side%delp_BC (is_0:ie_0,js_0:je_0,klev)) ; BC_side%delp_BC=real_snan
      allocate(BC_side%divgd_BC(is_0:ie_0,js_0:je_0,klev)) ; BC_side%divgd_BC=real_snan
!
      allocate(BC_side%q_BC    (is_0:ie_0,js_0:je_0,1:klev,1:ntracers)) ; BC_side%q_BC=real_snan
!
      if(.not.allocated(blend_this_tracer))then
        allocate(blend_this_tracer(1:ntracers))
        blend_this_tracer=.true.                                        !<-- Start with blending all tracers.
      endif
!
#ifndef SW_DYNAMICS
      allocate(BC_side%pt_BC   (is_0:ie_0,js_0:je_0,klev)) ; BC_side%pt_BC=real_snan
      allocate(BC_side%w_BC    (is_0:ie_0,js_0:je_0,klev)) ; BC_side%w_BC=real_snan
      allocate(BC_side%delz_BC (is_0:ie_0,js_0:je_0,klev)) ; BC_side%delz_BC=real_snan
      if(present(delz_side))then
        if(.not.associated(delz_side))then
          allocate(delz_side (is_0:ie_0,js_0:je_0,klev)) ; delz_side=real_snan
        endif
      endif
#ifdef USE_COND
      allocate(BC_side%q_con_BC(is_0:ie_0,js_0:je_0,klev)) ; BC_side%q_con_BC=real_snan
#ifdef MOIST_CAPPA
      allocate(BC_side%cappa_BC(is_0:ie_0,js_0:je_0,klev)) ; BC_side%cappa_BC=real_snan
#endif
#endif
#endif
!
!--------------------
!*** Wind components
!--------------------
!
!** D-grid u, C-grid v
!
      allocate(BC_side%u_BC (is_sn:ie_sn, js_sn:je_sn, klev)) ; BC_side%u_BC=real_snan
      allocate(BC_side%vc_BC(is_sn:ie_sn, js_sn:je_sn, klev)) ; BC_side%vc_BC=real_snan
!
!** C-grid u, D-grid v
!
      allocate(BC_side%uc_BC(is_we:ie_we, js_we:je_we, klev)) ; BC_side%uc_BC=real_snan
      allocate(BC_side%v_BC (is_we:ie_we, js_we:je_we, klev)) ; BC_side%v_BC=real_snan
!
!---------------------------------------------------------------------
!
      end subroutine allocate_regional_BC_arrays
!
!---------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------------------------------------------------------------------

subroutine remap_scalar_nggps_regional_bc(Atm                         &
                                         ,side                        &
                                         ,isd,ied,jsd,jed             &
                                         ,is_bc,ie_bc,js_bc,je_bc     &
                                         ,km, npz, ncnst, ak0, bk0    &
                                         ,psc, t_in, qa, omga, zh     &
                                         ,phis_reg                    &
                                         ,ps                          &
                                         ,BC_side )

  type(fv_atmos_type), intent(inout) :: Atm
  integer, intent(in):: isd,ied,jsd,jed          !<-- index limits of the Atm arrays w/halo=nhalo_model
  integer, intent(in):: is_bc,ie_bc,js_bc,je_bc  !<-- index limits of working arrays on boundary task subdomains (halo=nhalo_data)
  integer, intent(in):: km    &                  !<-- # of levels in 3-D input variables
                       ,npz   &                  !<-- # of levels in final 3-D integration variables
                       ,ncnst                    !<-- # of tracer variables
  real,    intent(in):: ak0(km+1), bk0(km+1)
  real,    intent(in), dimension(is_bc:ie_bc,js_bc:je_bc):: psc
  real,    intent(in), dimension(is_bc:ie_bc,js_bc:je_bc,km):: t_in
  real,    intent(in), dimension(is_bc:ie_bc,js_bc:je_bc,km):: omga
  real,    intent(in), dimension(is_bc:ie_bc,js_bc:je_bc,km,ncnst):: qa
  real,    intent(in), dimension(is_bc:ie_bc,js_bc:je_bc,km+1):: zh
  real,    intent(inout), dimension(isd-1:ied+1,jsd-1:jed+1):: phis_reg   !<-- Filtered sfc geopotential from preprocessing.
  real,    intent(out),dimension(is_bc:ie_bc,js_bc:je_bc) :: ps  !<-- sfc p in regional domain boundary region
  character(len=5),intent(in) :: side
  type(fv_regional_BC_variables),intent(inout) :: BC_side   !<-- The BC variables on a domain side at the final integration levels.

! local:
!
  real, dimension(:,:),allocatable :: pe0
  real, dimension(:,:),allocatable :: qn1
  real, dimension(:,:),allocatable :: dp2
  real, dimension(:,:),allocatable :: pe1
  real, dimension(:,:),allocatable :: qp
!
  real wk(is_bc:ie_bc,js_bc:je_bc)
  real, dimension(is_bc:ie_bc,js_bc:je_bc):: phis

!!! High-precision
  real(kind=R_GRID), dimension(is_bc:ie_bc,npz+1):: pn1
  real(kind=R_GRID):: gz_fv(npz+1)
  real(kind=R_GRID), dimension(2*km+1):: gz, pn
  real(kind=R_GRID), dimension(is_bc:ie_bc,km+1):: pn0
  real(kind=R_GRID):: pst
!!! High-precision
  integer i,ie,is,j,je,js,k,l,m, k2,iq
  integer  sphum, o3mr, liq_wat, ice_wat, rainwat, snowwat, graupel, cld_amt
!
!---------------------------------------------------------------------------------
!
  sphum   = sphum_index
  liq_wat = liq_water_index
  ice_wat = ice_water_index
  rainwat = rain_water_index
  snowwat = snow_water_index
  graupel = graupel_index
  cld_amt = cld_amt_index
  o3mr    = o3mr_index

  k2 = max(10, km/2)

  if (mpp_pe()==1) then
    print *, 'sphum = ', sphum
    print *, 'clwmr = ', liq_wat
    print *, ' o3mr = ', o3mr
    print *, 'ncnst = ', ncnst
    print *, 'ntracers = ', ntracers
  endif

  if ( sphum/=1 ) then
       call mpp_error(FATAL,'SPHUM must be 1st tracer')
  endif
!
!---------------------------------------------------------------------------------
!***  First compute over the extended boundary regions with halo=nhalo_data.
!***  This is needed to obtain pressures that will surround the wind points.
!---------------------------------------------------------------------------------
!
      is=is_bc
      if(side=='west')then
        is=ie_bc-nhalo_data-nrows_blend+1
      endif
!
      ie=ie_bc
      if(side=='east')then
        ie=is_bc+nhalo_data+nrows_blend-1
      endif
!
      js=js_bc
      if(side=='south')then
        js=je_bc-nhalo_data-nrows_blend+1
      endif
!
      je=je_bc
      if(side=='north')then
        je=js_bc+nhalo_data+nrows_blend-1
      endif
!
      allocate(pe0(is:ie,km+1)) ; pe0=real_snan
      allocate(qn1(is:ie,npz)) ; qn1=real_snan
      allocate(dp2(is:ie,npz)) ; dp2=real_snan
      allocate(pe1(is:ie,npz+1)) ; pe1=real_snan
      allocate(qp (is:ie,km)) ; qp=real_snan
!
!---------------------------------------------------------------------------------
      jloop1: do j=js,je
!---------------------------------------------------------------------------------
!
     do k=1,km+1
        do i=is,ie
           pe0(i,k) = ak0(k) + bk0(k)*psc(i,j)
           pn0(i,k) = log(pe0(i,k))
        enddo
     enddo

     do i=is,ie
        do k=1,km+1
           pn(k) = pn0(i,k)
           gz(k) = zh(i,j,k)*grav
        enddo
! Use log-p for interpolation/extrapolation
! mirror image method:
        do k=km+2, km+k2
               l = 2*(km+1) - k
           gz(k) = 2.*gz(km+1) - gz(l)
           pn(k) = 2.*pn(km+1) - pn(l)
        enddo

        do k=km+k2-1, 2, -1
          if( phis_reg(i,j).le.gz(k) .and. phis_reg(i,j).ge.gz(k+1) ) then
            pst = pn(k) + (pn(k+1)-pn(k))*(gz(k)-phis_reg(i,j))/(gz(k)-gz(k+1))
            go to 123
          endif
        enddo
  123   ps(i,j) = exp(pst)

     enddo   ! i-loop

!---------------------------------------------------------------------------------
     enddo jloop1
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
!***  Transfer values from the expanded boundary array for sfc pressure into
!***  the Atm object.
!---------------------------------------------------------------------------------
!
      is=lbound(Atm%ps,1)
      ie=ubound(Atm%ps,1)
      js=lbound(Atm%ps,2)
      je=ubound(Atm%ps,2)
!
      do j=js,je
      do i=is,ie
        Atm%ps(i,j)=ps(i,j)
      enddo
      enddo
!
!---------------------------------------------------------------------------------
!***  Now compute over the normal boundary regions with halo=nhalo_model
!***  extended through nrows_blend rows into the integration domain.
!***  Use the dimensions of one of the permanent BC variables in Atm
!***  as the loop limits so any side of the domain can be addressed.
!---------------------------------------------------------------------------------
!
      is=lbound(BC_side%delp_BC,1)
      ie=ubound(BC_side%delp_BC,1)
      js=lbound(BC_side%delp_BC,2)
      je=ubound(BC_side%delp_BC,2)
!
!---------------------------------------------------------------------------------
    jloop2: do j=js,je
!---------------------------------------------------------------------------------
     do k=1,km+1
        do i=is,ie
           pe0(i,k) = ak0(k) + bk0(k)*psc(i,j)
           pn0(i,k) = log(pe0(i,k))
        enddo
      enddo
!
     do i=is,ie
        pe1(i,1) = Atm%ak(1)
        pn1(i,1) = log(pe1(i,1))
     enddo
     do k=2,npz+1
       do i=is,ie
          pe1(i,k) = Atm%ak(k) + Atm%bk(k)*ps(i,j)
          pn1(i,k) = log(pe1(i,k))
       enddo
     enddo

! * Compute delp
      do k=1,npz
        do i=is,ie
          dp2(i,k) = pe1(i,k+1) - pe1(i,k)
          BC_side%delp_BC(i,j,k) = dp2(i,k)
        enddo
      enddo

! map shpum, o3mr, liq_wat tracers
      do iq=1,ncnst
!       if (iq == sphum .or. iq == liq_wat .or.  iq == o3mr) then ! only remap if the data is already set
       if (iq /= cld_amt) then ! don't remap cld_amt
         do k=1,km
            do i=is,ie
               qp(i,k) = qa(i,j,k,iq)
            enddo
         enddo

         call mappm(km, pe0, qp, npz, pe1,  qn1, is,ie, 0, 8, Atm%ptop)

         if ( iq==sphum ) then
            call fillq(ie-is+1, npz, 1, qn1, dp2)
         else
            call fillz(ie-is+1, npz, 1, qn1, dp2)
         endif
! The HiRam step of blending model sphum with NCEP data is obsolete because nggps is always cold starting...
         do k=1,npz
           do i=is,ie
             BC_side%q_BC(i,j,k,iq) = qn1(i,k)
           enddo
         enddo
       endif
      enddo

!---------------------------------------------------
! Retrieve temperature using GFS geopotential height
!---------------------------------------------------
!
      i_loop: do i=is,ie
!
! Make sure FV3 top is lower than GFS; can not do extrapolation above the top at this point
        if ( pn1(i,1) .lt. pn0(i,1) ) then
          call mpp_error(FATAL,'FV3 top higher than NCEP/GFS')
        endif

        do k=1,km+1
           pn(k) = pn0(i,k)
           gz(k) = zh(i,j,k)*grav
        enddo
!-------------------------------------------------
        do k=km+2, km+k2
           l = 2*(km+1) - k
           gz(k) = 2.*gz(km+1) - gz(l)
           pn(k) = 2.*pn(km+1) - pn(l)
        enddo
!-------------------------------------------------

        gz_fv(npz+1) = phis_reg(i,j)

        m = 1

        do k=1,npz
! Searching using FV3 log(pe): pn1
#ifdef USE_ISOTHERMO
           do l=m,km
              if ( (pn1(i,k).le.pn(l+1)) .and. (pn1(i,k).ge.pn(l)) ) then
                  gz_fv(k) = gz(l) + (gz(l+1)-gz(l))*(pn1(i,k)-pn(l))/(pn(l+1)-pn(l))
                  goto 555
              elseif ( pn1(i,k) .gt. pn(km+1) ) then
! Isothermal under ground; linear in log-p extra-polation
                  gz_fv(k) = gz(km+1) + (gz_fv(npz+1)-gz(km+1))*(pn1(i,k)-pn(km+1))/(pn1(i,npz+1)-pn(km+1))
                  goto 555
              endif
           enddo
#else
           do l=m,km+k2-1
              if ( (pn1(i,k).le.pn(l+1)) .and. (pn1(i,k).ge.pn(l)) ) then
                  gz_fv(k) = gz(l) + (gz(l+1)-gz(l))*(pn1(i,k)-pn(l))/(pn(l+1)-pn(l))
                  goto 555
              endif
           enddo
#endif
555     m = l
        enddo

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxx  DO WE NEED Atm%peln to have values in the boundary region?
!xxx  FOR NOW COMMENT IT OUT.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxx  do k=1,npz+1
!xxx     Atm%peln(i,k,j) = pn1(i,k)
!xxx  enddo

! Compute true temperature using hydrostatic balance if not read from input.

        if (data_source /= 'FV3GFS GAUSSIAN NEMSIO FILE') then
          do k=1,npz
            BC_side%pt_BC(i,j,k) = (gz_fv(k)-gz_fv(k+1))/( rdgas*(pn1(i,k+1)-pn1(i,k))*(1.+zvir*BC_side%q_BC(i,j,k,sphum)) )
          enddo
        endif


        if ( .not. Atm%flagstruct%hydrostatic ) then
          do k=1,npz
            BC_side%delz_BC(i,j,k) = (gz_fv(k+1) - gz_fv(k)) / grav
          enddo
        endif

      enddo i_loop

!-----------------------------------------------------------------------
! separate cloud water and cloud ice
! From Jan-Huey Chen's HiRAM code
!-----------------------------------------------------------------------
!
! If the source is FV3GFS GAUSSIAN NEMSIO FILE then all the tracers are in the boundary files
! and will be read in.
! If the source is from old GFS or operational GSM then the tracers will be fixed in the boundaries
! and may not provide a very good result
! 
!  if (cld_amt .gt. 0) BC_side%q_BC(:,:,:,cld_amt) = 0.
  if (trim(data_source) /= 'FV3GFS GAUSSIAN NEMSIO FILE') then
   if ( Atm%flagstruct%nwat .eq. 6 ) then
      do k=1,npz
         do i=is,ie
            qn1(i,k) = BC_side%q_BC(i,j,k,liq_wat)
            BC_side%q_BC(i,j,k,rainwat) = 0.
            BC_side%q_BC(i,j,k,snowwat) = 0.
            BC_side%q_BC(i,j,k,graupel) = 0.
            if ( BC_side%pt_BC(i,j,k) > 273.16 ) then       ! > 0C all liq_wat
               BC_side%q_BC(i,j,k,liq_wat) = qn1(i,k)
               BC_side%q_BC(i,j,k,ice_wat) = 0.
#ifdef ORIG_CLOUDS_PART
            else if ( BC_side%pt_BC(i,j,k) < 258.16 ) then  ! < -15C all ice_wat
               BC_side%q_BC(i,j,k,liq_wat) = 0.
               BC_side%q_BC(i,j,k,ice_wat) = qn1(i,k)
            else                                     ! between -15~0C: linear interpolation
               BC_side%q_BC(i,j,k,liq_wat) = qn1(i,k)*((BC_side%pt_BC(i,j,k)-258.16)/15.)
               BC_side%q_BC(i,j,k,ice_wat) = qn1(i,k) - BC_side%q_BC(i,j,k,liq_wat)
            endif
#else
            else if ( BC_side%pt_BC(i,j,k) < 233.16 ) then  ! < -40C all ice_wat
               BC_side%q_BC(i,j,k,liq_wat) = 0.
               BC_side%q_BC(i,j,k,ice_wat) = qn1(i,k)
            else
               if ( k.eq.1 ) then  ! between [-40,0]: linear interpolation
                  BC_side%q_BC(i,j,k,liq_wat) = qn1(i,k)*((BC_side%pt_BC(i,j,k)-233.16)/40.)
                  BC_side%q_BC(i,j,k,ice_wat) = qn1(i,k) - BC_side%q_BC(i,j,k,liq_wat)
               else
                 if (BC_side%pt_BC(i,j,k)<258.16 .and. BC_side%q_BC(i,j,k-1,ice_wat)>1.e-5 ) then
                    BC_side%q_BC(i,j,k,liq_wat) = 0.
                    BC_side%q_BC(i,j,k,ice_wat) = qn1(i,k)
                 else  ! between [-40,0]: linear interpolation
                    BC_side%q_BC(i,j,k,liq_wat) = qn1(i,k)*((BC_side%pt_BC(i,j,k)-233.16)/40.)
                    BC_side%q_BC(i,j,k,ice_wat) = qn1(i,k) - BC_side%q_BC(i,j,k,liq_wat)
                 endif
               endif
            endif
#endif
            call mp_auto_conversion(BC_side%q_BC(i,j,k,liq_wat), BC_side%q_BC(i,j,k,rainwat),  &
                                    BC_side%q_BC(i,j,k,ice_wat), BC_side%q_BC(i,j,k,snowwat) )
         enddo
      enddo
   endif
  endif ! data source /= FV3GFS GAUSSIAN NEMSIO FILE
!
! For GFS spectral input, omega in pa/sec is stored as w in the input data so actual w(m/s) is calculated
! For GFS nemsio input, omega is 0, so best not to use for input since boundary data will not exist for w
! For FV3GFS NEMSIO input, w is already in m/s (but the code reads in as omga) and just needs to be remapped
!-------------------------------------------------------------
! map omega
!------- ------------------------------------------------------
   if ( .not. Atm%flagstruct%hydrostatic ) then
      do k=1,km
         do i=is,ie
            qp(i,k) = omga(i,j,k)
         enddo
      enddo

      call mappm(km, pe0, qp, npz, pe1, qn1, is,ie, -1, 4, Atm%ptop)

      if (data_source == 'FV3GFS GAUSSIAN NEMSIO FILE') then
        do k=1,npz
          do i=is,ie
            BC_side%w_BC(i,j,k) = qn1(i,k)
          enddo
        enddo
!------------------------------
! Remap input T linearly in p.
!------------------------------
        do k=1,km
          do i=is,ie
            qp(i,k) = t_in(i,j,k)
          enddo
        enddo

        call mappm(km, pe0, qp, npz, pe1, qn1, is,ie, 2, 4, Atm%ptop)

        do k=1,npz
          do i=is,ie
            BC_side%pt_BC(i,j,k) = qn1(i,k)
          enddo
        enddo

      else          !<-- datasource /= 'FV3GFS GAUSSIAN NEMSIO FILE'
        do k=1,npz
          do i=is,ie
            BC_side%w_BC(i,j,k) = qn1(i,k)/BC_side%delp_BC(i,j,k)*BC_side%delz_BC(i,j,k)
          enddo
        enddo
      endif

   endif   !.not. Atm%flagstruct%hydrostatic

   enddo jloop2

! Add some diagnostics:
!   call p_maxmin('PS_model (mb)', Atm%ps(is:ie,js:je), is, ie, js, je, 1, 0.01)
!   call p_maxmin('PT_model', Atm%pt(is:ie,js:je,1:npz), is, ie, js, je, npz, 1.)
  do j=js,je
     do i=is,ie
        wk(i,j) = phis_reg(i,j)/grav - zh(i,j,km+1)
     enddo
  enddo
!   call pmaxmn('ZS_diff (m)', wk, is, ie, js, je, 1, 1., Atm%gridstruct%area_64, Atm%domain)

  do j=js,je
     do i=is,ie
        wk(i,j) = ps(i,j) - psc(i,j)
     enddo
  enddo
!   call pmaxmn('PS_diff (mb)', wk, is, ie, js, je, 1, 0.01, Atm%gridstruct%area_64, Atm%domain)
  deallocate (pe0,qn1,dp2,pe1,qp)
  if (is_master()) write(*,*) 'done remap_scalar_nggps_regional_bc'
!---------------------------------------------------------------------

 end subroutine remap_scalar_nggps_regional_bc

!---------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------------------------------------------------------------------

 subroutine remap_dwinds_regional_bc(Atm                              &
                                    ,is_input,ie_input                &
                                    ,js_input,je_input                &
                                    ,is_u,ie_u,js_u,je_u              &
                                    ,is_v,ie_v,js_v,je_v              &
                                    ,km, npz                          &
                                    ,ak0, bk0                         &
                                    ,psc, ud, vd, uc, vc              &
                                    ,BC_side )
  type(fv_atmos_type), intent(inout) :: Atm
  integer, intent(in):: is_input, ie_input, js_input, je_input   !<-- index limits of the boundary arrays with nahlo=nhalo_data
  integer, intent(in):: is_u,ie_u,js_u,je_u          !<-- index limits of D-grid u in this boundary region
  integer, intent(in):: is_v,ie_v,js_v,je_v          !<-- index limits of D-grid v in this boundary region
  integer, intent(in):: km    &                      !<-- # of levels in 3-D input variables
                       ,npz                          !<-- # of levels in final 3-D integration variables
  real,    intent(in):: ak0(km+1), bk0(km+1)

  real, intent(in) :: psc(is_input:ie_input,js_input:je_input)

  real,    intent(in)::  ud(is_u:ie_u,js_u:je_u,km)
  real,    intent(in)::  vc(is_u:ie_u,js_u:je_u,km)
  real,    intent(in)::  vd(is_v:ie_v,js_v:je_v,km)
  real,    intent(in)::  uc(is_v:ie_v,js_v:je_v,km)
  type(fv_regional_BC_variables),intent(inout) :: BC_side   !<-- The BC variables on a domain side at the final integration levels.
! local:
      real, dimension(:,:),allocatable :: pe0
      real, dimension(:,:),allocatable :: pe1
      real, dimension(:,:),allocatable :: qn1_d,qn1_c
  integer i,j,k

      allocate(pe0  (is_u:ie_u, km+1)) ; pe0=real_snan
      allocate(pe1  (is_u:ie_u, npz+1)) ; pe1=real_snan
      allocate(qn1_d(is_u:ie_u, npz)) ; qn1_d=real_snan
      allocate(qn1_c(is_u:ie_u, npz)) ; qn1_c=real_snan

!----------------------------------------------------------------------------------------------
    j_loopu: do j=js_u,je_u
!----------------------------------------------------------------------------------------------

!------
! map u
!------
     do k=1,km+1
        do i=is_u,ie_u
           pe0(i,k) = ak0(k) + bk0(k)*0.5*(psc(i,j-1)+psc(i,j))
        enddo
     enddo
     do k=1,npz+1
        do i=is_u,ie_u
           pe1(i,k) = Atm%ak(k) + Atm%bk(k)*0.5*(psc(i,j-1)+psc(i,j))
        enddo
     enddo
     call mappm(km, pe0(is_u:ie_u,1:km+1), ud(is_u:ie_u,j,1:km), npz, pe1(is_u:ie_u,1:npz+1),   &
                qn1_d(is_u:ie_u,1:npz), is_u,ie_u, -1, 8, Atm%ptop )
     call mappm(km, pe0(is_u:ie_u,1:km+1), vc(is_u:ie_u,j,1:km), npz, pe1(is_u:ie_u,1:npz+1),   &
                qn1_c(is_u:ie_u,1:npz), is_u,ie_u, -1, 8, Atm%ptop )
     do k=1,npz
        do i=is_u,ie_u
           BC_side%u_BC(i,j,k) = qn1_d(i,k)
           BC_side%vc_BC(i,j,k) = qn1_c(i,k)
        enddo
     enddo

     enddo j_loopu

      deallocate(pe0)
      deallocate(pe1)
      deallocate(qn1_d)
      deallocate(qn1_c)

      allocate(pe0  (is_v:ie_v, km+1)) ; pe0=real_snan
      allocate(pe1  (is_v:ie_v, npz+1)) ; pe1=real_snan
      allocate(qn1_d(is_v:ie_v, npz)) ; qn1_d=real_snan
      allocate(qn1_c(is_v:ie_v, npz)) ; qn1_c=real_snan

!----------------------------------------------------------------------------------------------
  j_loopv: do j=js_v,je_v
!----------------------------------------------------------------------------------------------
!
!------
! map v
!------

     do k=1,km+1
        do i=is_v,ie_v
           pe0(i,k) = ak0(k) + bk0(k)*0.5*(psc(i-1,j)+psc(i,j))
        enddo
     enddo
     do k=1,npz+1
        do i=is_v,ie_v
           pe1(i,k) = Atm%ak(k) + Atm%bk(k)*0.5*(psc(i-1,j)+psc(i,j))
        enddo
     enddo
     call mappm(km, pe0(is_v:ie_v,1:km+1), vd(is_v:ie_v,j,1:km), npz, pe1(is_v:ie_v,1:npz+1),  &
                qn1_d(is_v:ie_v,1:npz), is_v,ie_v, -1, 8, Atm%ptop)
     call mappm(km, pe0(is_v:ie_v,1:km+1), uc(is_v:ie_v,j,1:km), npz, pe1(is_v:ie_v,1:npz+1),  &
                qn1_c(is_v:ie_v,1:npz), is_v,ie_v, -1, 8, Atm%ptop)
     do k=1,npz
        do i=is_v,ie_v
           BC_side%v_BC(i,j,k) = qn1_d(i,k)
           BC_side%uc_BC(i,j,k) = qn1_c(i,k)
        enddo
     enddo

      enddo j_loopv

      deallocate(pe0)
      deallocate(pe1)
      deallocate(qn1_d)
      deallocate(qn1_c)

  if (is_master()) write(*,*) 'done remap_dwinds'

 end subroutine remap_dwinds_regional_bc

!---------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------------------------------------------------------------------

      subroutine set_regional_BCs(delp,delz,w,pt                      &
#ifdef USE_COND
                                 ,q_con                               &
#endif
#ifdef MOIST_CAPPA
                                 ,cappa                               &
#endif
                                 ,q                                   &
                                 ,u,v,uc,vc                           &
                                 ,bd, nlayers                        &
                                 ,fcst_time )
!
!---------------------------------------------------------------------
!***  Select the boundary variables' boundary data at the two
!***  bracketing time levels and apply them to the updating 
!***  of the variables' boundary regions at the appropriate
!***  forecast time.  This is done at the beginning of every 
!***  large timestep in fv_dynamics.
!---------------------------------------------------------------------
      implicit none
!---------------------------------------------------------------------
!
!--------------------
!***  Input variables
!--------------------
!
      integer,intent(in) :: nlayers
!
      real,intent(in) :: fcst_time                                       !<-- Current forecast time (sec)
!
      type(fv_grid_bounds_type),intent(in) :: bd                         !<-- Task subdomain indices
!
!----------------------
!***  Output variables
!----------------------
!
      real,dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz),intent(out) ::  &
                                                                delp  &
                                                               ,pt 
! 
      real,dimension(bd%isd:,bd%jsd:,1:),intent(out) :: w
      real,dimension(bd%is:,bd%js:,1:),intent(out) :: delz
#ifdef USE_COND
      real,dimension(bd%isd:,bd%jsd:,1:),intent(out) :: q_con
#endif

!
      real,dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz,ntracers),intent(out) :: q
!
#ifdef MOIST_CAPPA
      real,dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz),intent(out) :: cappa
!#else
!      real,dimension(isd:isd,jsd:jsd,1),intent(out) :: cappa
#endif
!
      real,dimension(bd%isd:bd%ied,bd%jsd:bd%jed+1,npz),intent(out) :: u,vc
!
      real,dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed,npz),intent(out) :: uc,v
!
!---------------------
!***  Local variables
!---------------------
!
      real :: fraction_interval
!
!---------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
!***  The current forecast time is this fraction of the way from
!***  time level 0 to time level 1.
!---------------------------------------------------------------------
!
      fraction_interval=mod(fcst_time,(bc_update_interval*3600.))     &
                       /(bc_update_interval*3600.)
!
!---------------------------------------------------------------------
!
      if(north_bc)then
        call bc_values_into_arrays(BC_t0%north,BC_t1%north            &
                                  ,'north'                            &
                                  ,bd%isd                             &
                                  ,bd%ied                             &
                                  ,bd%jsd                             &
                                  ,bd%js-1                            &
                                  ,bd%isd                             &
                                  ,bd%ied                             &
                                  ,bd%jsd                             &
                                  ,bd%js-1                            &
                                  ,bd%isd                             &
                                  ,bd%ied+1                           &
                                  ,bd%jsd                             &
                                  ,bd%js-1)                       
      endif
!
      if(south_bc)then
        call bc_values_into_arrays(BC_t0%south,BC_t1%south            &
                                  ,'south'                            &
                                  ,bd%isd                             &
                                  ,bd%ied                             &
                                  ,bd%je+1                            &
                                  ,bd%jed                             &
                                  ,bd%isd                             &
                                  ,bd%ied                             &
                                  ,bd%je+2                            &
                                  ,bd%jed+1                           &
                                  ,bd%isd                             &
                                  ,bd%ied+1                           &
                                  ,bd%je+1                            &
                                  ,bd%jed )                       
      endif
!
      if(east_bc)then
        call bc_values_into_arrays(BC_t0%east,BC_t1%east              &
                                  ,'east '                            &
                                  ,bd%isd                             &
                                  ,bd%is-1                            &
                                  ,bd%js                              &
                                  ,bd%je                              &
                                  ,bd%isd                             &
                                  ,bd%is-1                            &
                                  ,bd%js                              &
                                  ,bd%je+1                            &
                                  ,bd%isd                             &
                                  ,bd%is-1                            &
                                  ,bd%js                              &
                                  ,bd%je  )                       
      endif
!
      if(west_bc)then
        call bc_values_into_arrays(BC_t0%west,BC_t1%west             &
                                  ,'west '                            &
                                  ,bd%ie+1                            &
                                  ,bd%ied                             &
                                  ,bd%js                              &
                                  ,bd%je                              &
                                  ,bd%ie+1                            &
                                  ,bd%ied                             &
                                  ,bd%js                              &
                                  ,bd%je+1                            &
                                  ,bd%ie+2                            &
                                  ,bd%ied+1                           &
                                  ,bd%js                              &
                                  ,bd%je  )                       
      endif
!
!---------------------------------------------------------------------

      contains

!---------------------------------------------------------------------
!
      subroutine bc_values_into_arrays(side_t0,side_t1                &
                                      ,side                           &
                                      ,i1,i2,j1,j2                    &
                                      ,i1_uvs,i2_uvs,j1_uvs,j2_uvs    &
                                      ,i1_uvw,i2_uvw,j1_uvw,j2_uvw )
!
!---------------------------------------------------------------------
!***  Apply boundary values to the prognostic arrays at the 
!***  desired time.
!---------------------------------------------------------------------
      implicit none
!---------------------------------------------------------------------
!
!---------------------
!***  Input arguments
!---------------------
!
      type(fv_regional_BC_variables),intent(in) :: side_t0            &
                                                  ,side_t1
!
      character(len=*),intent(in) :: side
!
      integer,intent(in) :: i1,i2,j1,j2                               &
                           ,i1_uvs,i2_uvs,j1_uvs,j2_uvs               &
                           ,i1_uvw,i2_uvw,j1_uvw,j2_uvw
!
!---------------------
!***  Local arguments
!---------------------
!
      integer :: i,ie,j,je,jend,jend_uvs,jend_uvw                     &
                ,jstart,jstart_uvs,jstart_uvw,k,nt,nz
!
      real,dimension(:,:,:),pointer :: delz_ptr
!
!---------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------
!
      jstart=j1
      jend  =j2
      jstart_uvs=j1_uvs
      jend_uvs  =j2_uvs
      jstart_uvw=j1_uvw
      jend_uvw  =j2_uvw
      if((trim(side)=='east'.or.trim(side)=='west').and..not.north_bc)then
        jstart=j1-nhalo_model
        jstart_uvs=j1_uvs-nhalo_model
        jstart_uvw=j1_uvw-nhalo_model
      endif
      if((trim(side)=='east'.or.trim(side)=='west').and..not.south_bc)then
        jend=j2+nhalo_model
        jend_uvs=j2_uvs+nhalo_model
        jend_uvw=j2_uvw+nhalo_model
      endif
!
      select case (trim(side))
        case ('north')
          delz_ptr=>delz_auxiliary%north
        case ('south')
          delz_ptr=>delz_auxiliary%south
        case ('east') 
          delz_ptr=>delz_auxiliary%east
        case ('west') 
          delz_ptr=>delz_auxiliary%west
      end select
!
      do k=1,nlayers
        do j=jstart,jend
        do i=i1,i2
          delp(i,j,k)=side_t0%delp_BC(i,j,k)                          &
                     +(side_t1%delp_BC(i,j,k)-side_t0%delp_BC(i,j,k)) &
                      *fraction_interval
          pt(i,j,k)=side_t0%pt_BC(i,j,k)                              &
                     +(side_t1%pt_BC(i,j,k)-side_t0%pt_BC(i,j,k))     &
                      *fraction_interval
!          delz(i,j,k)=side_t0%delz_BC(i,j,k)                            &
!                     +(side_t1%delz_BC(i,j,k)-side_t0%delz_BC(i,j,k))   &
!                      *fraction_interval
           delz_ptr(i,j,k)=side_t0%delz_BC(i,j,k)                            &
                          +(side_t1%delz_BC(i,j,k)-side_t0%delz_BC(i,j,k))   &
                           *fraction_interval
#ifdef MOIST_CAPPA
          cappa(i,j,k)=side_t0%cappa_BC(i,j,k)                          &
                     +(side_t1%cappa_BC(i,j,k)-side_t0%cappa_BC(i,j,k)) &
                      *fraction_interval
#endif
#ifdef USE_COND
          q_con(i,j,k)=side_t0%q_con_BC(i,j,k)                          &
                     +(side_t1%q_con_BC(i,j,k)-side_t0%q_con_BC(i,j,k)) &
                      *fraction_interval
#endif
          w(i,j,k)=side_t0%w_BC(i,j,k)                                  &
                     +(side_t1%w_BC(i,j,k)-side_t0%w_BC(i,j,k))         &
                      *fraction_interval
        enddo
        enddo
!
        do j=jstart_uvs,jend_uvs
        do i=i1_uvs,i2_uvs
          u(i,j,k)=side_t0%u_BC(i,j,k)                                &
                     +(side_t1%u_BC(i,j,k)-side_t0%u_BC(i,j,k))       &
                      *fraction_interval
          vc(i,j,k)=side_t0%vc_BC(i,j,k)                              &
                     +(side_t1%vc_BC(i,j,k)-side_t0%vc_BC(i,j,k))     &
                      *fraction_interval
        enddo
        enddo
!
        do j=jstart_uvw,jend_uvw
        do i=i1_uvw,i2_uvw
          v(i,j,k)=side_t0%v_BC(i,j,k)                                &
                     +(side_t1%v_BC(i,j,k)-side_t0%v_BC(i,j,k))       &
                      *fraction_interval
          uc(i,j,k)=side_t0%uc_BC(i,j,k)                              &
                     +(side_t1%uc_BC(i,j,k)-side_t0%uc_BC(i,j,k))     &
                      *fraction_interval
        enddo
        enddo
      enddo
!
      ie=min(ubound(side_t0%w_BC,1),ubound(w,1))
      je=min(ubound(side_t0%w_BC,2),ubound(w,2))
      nz=ubound(w,3)
!
      do nt=1,ntracers
        do k=1,nz
          do j=jstart,jend
          do i=i1,i2
            q(i,j,k,nt)=side_t0%q_BC(i,j,k,nt)                            &
                       +(side_t1%q_BC(i,j,k,nt)-side_t0%q_BC(i,j,k,nt))   &
                        *fraction_interval
          enddo
          enddo
        enddo
      enddo
!
!---------------------------------------------------------------------
!
      end subroutine bc_values_into_arrays
!
!---------------------------------------------------------------------
!
      end subroutine set_regional_BCs
!
!---------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------------------------------------------------------------------
      subroutine regional_boundary_update(array                       &
                                         ,bc_vbl_name                 &
                                         ,lbnd_x,ubnd_x               &
                                         ,lbnd_y,ubnd_y               &
                                         ,ubnd_z                      &
                                         ,is,ie,js,je                 &
                                         ,isd,ied,jsd,jed             &
                                         ,fcst_time                   &
                                         ,index4 )
!
!---------------------------------------------------------------------
!***  Select the given variable's boundary data at the two
!***  bracketing time levels and apply them to the updating
!***  of the variable's boundary region at the appropriate
!***  forecast time.
!---------------------------------------------------------------------
      implicit none
!---------------------------------------------------------------------
!
!--------------------
!***  Input variables
!--------------------
!
      integer,intent(in) :: lbnd_x,ubnd_x,lbnd_y,ubnd_y,ubnd_z           !<-- Dimensions of full prognostic array to be updated.
!
      integer,intent(in) :: is,ie,js,je                               &  !<-- Compute limits
                           ,isd,ied,jsd,jed                              !<-- Memory limits
!
      integer,intent(in),optional :: index4                              !<-- Index for the 4-D tracer array.
!
      real,intent(in) :: fcst_time                                       !<-- Forecast time (sec) at which BC update is applied.
!
      character(len=*),intent(in) :: bc_vbl_name                         !<-- Name of the variable to be updated.
!
!----------------------
!***  Output variables
!----------------------
!
      real,dimension(lbnd_x:ubnd_x,lbnd_y:ubnd_y,1:ubnd_z)            &
                                              ,intent(out) :: array      !<-- Update this full array's boundary region.
!
!---------------------
!***  Local variables
!---------------------
!
      integer :: i1,i2,j1,j2                                             !<-- Horizontal limits of region updated.
      integer :: i_bc,j_bc                                               !<-- Innermost bndry index (anchor point for blending)
      integer :: i1_blend,i2_blend,j1_blend,j2_blend                     !<-- Limits of updated blending region.
      integer :: lbnd1,ubnd1,lbnd2,ubnd2                                 !<-- Horizontal limits of BC update arrays.
      integer :: iq                                                      !<-- Tracer index
      integer :: nside
!
      real,dimension(:,:,:),pointer :: bc_t0,bc_t1                       !<-- Boundary data at the two bracketing times.
!
      logical :: blend,call_interp
!
!---------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------
!
      if(.not.(north_bc.or.south_bc.or.east_bc.or.west_bc))then
        return
      endif
!
      blend=.true.
      iq=0
      if(present(index4))then
        iq=index4
        blend=blend_this_tracer(iq)
      endif
!
!---------------------------------------------------------------------
!***  Loop through the sides of the domain and find the limits
!***  of the region to update in the boundary.
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      sides: do nside=1,4
!---------------------------------------------------------------------
!
        call_interp=.false.
!
!-----------
!***  North
!-----------
!
        if(nside==1)then
          if(north_bc)then
            call_interp=.true.
            bc_side_t0=>bc_north_t0
            bc_side_t1=>bc_north_t1
!
            i1=isd
            i2=ied
            if(trim(bc_vbl_name)=='uc'.or.trim(bc_vbl_name)=='v')then
              i2=ied+1
            endif
!
            j1=jsd
            j2=js-1
!
            i1_blend=is
            i2_blend=ie
            if(trim(bc_vbl_name)=='uc'.or.trim(bc_vbl_name)=='v')then
              i2_blend=ie+1
            endif
            j1_blend=js
            j2_blend=js+nrows_blend_user-1
            i_bc=-9e9
            j_bc=j2
!
          endif
        endif
!
!-----------
!***  South
!-----------
!
        if(nside==2)then
          if(south_bc)then
            call_interp=.true.
            bc_side_t0=>bc_south_t0
            bc_side_t1=>bc_south_t1
!
            i1=isd
            i2=ied
            if(trim(bc_vbl_name)=='uc'.or.trim(bc_vbl_name)=='v')then
              i2=ied+1
            endif
!
            j1=je+1
            j2=jed
            if(trim(bc_vbl_name)=='u'.or.trim(bc_vbl_name)=='vc')then
              j1=je+2
              j2=jed+1
            endif
!
            i1_blend=is
            i2_blend=ie
            if(trim(bc_vbl_name)=='uc'.or.trim(bc_vbl_name)=='v')then
              i2_blend=ie+1
            endif
            j2_blend=je
            if(trim(bc_vbl_name)=='u'.or.trim(bc_vbl_name)=='vc')then
              j2_blend=je+1
            endif
            j1_blend=j2_blend-nrows_blend_user+1
            i_bc=-9e9
            j_bc=j1
!
          endif
        endif
!
!----------
!***  East
!----------
!
        if(nside==3)then
          if(east_bc)then
            call_interp=.true.
            bc_side_t0=>bc_east_t0
            bc_side_t1=>bc_east_t1
!
            j1=jsd
            j2=jed
!
            i1=isd
            i2=is-1
!
            if(north_bc)then
              j1=js
            endif
            if(south_bc)then
              j2=je
              if(trim(bc_vbl_name)=='u'.or.trim(bc_vbl_name)=='vc')then
                j2=je+1
              endif
            endif
!
            i1_blend=is
            i2_blend=is+nrows_blend_user-1
            j1_blend=js
            j2_blend=je
            if(north_bc)then
              j1_blend=js+nrows_blend_user     !<-- North BC already handles nrows_blend_user blending rows
            endif
            if(south_bc)then
              j2_blend=je-nrows_blend_user     !<-- South BC already handles nrows_blend_user blending rows
            endif
            if(trim(bc_vbl_name)=='u'.or.trim(bc_vbl_name)=='vc')then
              j2_blend=j2_blend+1
            endif
            i_bc=i2
            j_bc=-9e9
!
          endif
        endif
!
!----------
!***  West
!----------
!
        if(nside==4)then
          if(west_bc)then
            call_interp=.true.
            bc_side_t0=>bc_west_t0
            bc_side_t1=>bc_west_t1
!
            j1=jsd
            j2=jed
!
            i1=ie+1
            i2=ied
            if(trim(bc_vbl_name)=='uc'.or.trim(bc_vbl_name)=='v')then
              i1=ie+2
              i2=ied+1
            endif
!
            if(north_bc)then
              j1=js
            endif
            if(south_bc)then
              j2=je
              if(trim(bc_vbl_name)=='u'.or.trim(bc_vbl_name)=='vc')then
                j2=je+1
              endif
            endif
!
            i1_blend=i1-nrows_blend_user     
            i2_blend=i1-1
            j1_blend=js
            j2_blend=je
            if(north_bc)then
              j1_blend=js+nrows_blend_user   !<-- North BC already handled nrows_blend_user blending rows.
            endif
            if(south_bc)then
              j2_blend=je-nrows_blend_user   !<-- South BC already handled nrows_blend_user blending rows.
            endif
            if(trim(bc_vbl_name)=='u'.or.trim(bc_vbl_name)=='vc')then
              j2_blend=j2_blend+1
            endif
            i_bc=i1
            j_bc=-9e9
!
          endif
        endif
!
!---------------------------------------------------------------------
!***  Get the pointers pointing at the boundary arrays holding the
!***  two time levels of the given prognostic array's boundary region
!***  then update the boundary points.
!---------------------------------------------------------------------
!
        if(call_interp)then
!
          call retrieve_bc_variable_data(bc_vbl_name                  &
                                        ,bc_side_t0,bc_side_t1        &  !<-- Boundary data objects
                                        ,bc_t0,bc_t1                  &  !<-- Pointer to boundary arrays
                                        ,lbnd1,ubnd1,lbnd2,ubnd2      &  !<-- Bounds of the boundary data objects
                                        ,iq )
!
          call bc_time_interpolation(array                               &
                                    ,lbnd_x,ubnd_x,lbnd_y,ubnd_y,ubnd_z  &
                                    ,bc_t0,bc_t1                         &
                                    ,lbnd1,ubnd1,lbnd2,ubnd2             &
                                    ,i1,i2,j1,j2                         &
                                    ,is,ie,js,je                         &
                                    ,fcst_time                           &
                                    ,bc_update_interval                  &
                                    ,i1_blend,i2_blend,j1_blend,j2_blend &
                                    ,i_bc,j_bc,nside,bc_vbl_name,blend )
        endif
!
!---------------------------------------------------------------------
!
      enddo sides
!
!---------------------------------------------------------------------
!
      end subroutine regional_boundary_update

!---------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------------------------------------------------------------------

      subroutine retrieve_bc_variable_data(bc_vbl_name                &
                                          ,bc_side_t0,bc_side_t1      &
                                          ,bc_t0,bc_t1                &
                                          ,lbnd1,ubnd1,lbnd2,ubnd2    &
                                          ,iq )
                                      
!---------------------------------------------------------------------
!***  Select the boundary variable associated with the prognostic
!***  array that needs its boundary region to be updated.
!---------------------------------------------------------------------
      implicit none
!---------------------------------------------------------------------
!
!---------------------
!***  Input variables
!---------------------
!
      integer,intent(in) :: iq                                           !<-- Index used by 4-D tracer array.
!
      character(len=*),intent(in) :: bc_vbl_name
!
      type(fv_regional_BC_variables),pointer :: bc_side_t0,bc_side_t1    !<-- Boundary states for the given domain side.
!
!
!----------------------
!***  Output variables
!----------------------
!
      integer,intent(out) :: lbnd1,ubnd1,lbnd2,ubnd2                     !<-- Horizontal dimensions of boundary array
!
      real,dimension(:,:,:),pointer :: bc_t0,bc_t1                       !<-- Boundary state values for the desired variable.
!
!---------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------
!
      select case (bc_vbl_name)
!
        case ('delp')
          bc_t0=>bc_side_t0%delp_BC
          bc_t1=>bc_side_t1%delp_BC
        case ('delz')
          bc_t0=>bc_side_t0%delz_BC
          bc_t1=>bc_side_t1%delz_BC
        case ('pt')
          bc_t0=>bc_side_t0%pt_BC
          bc_t1=>bc_side_t1%pt_BC
        case ('w')
          bc_t0=>bc_side_t0%w_BC
          bc_t1=>bc_side_t1%w_BC
        case ('divgd')
          bc_t0=>bc_side_t0%divgd_BC
          bc_t1=>bc_side_t1%divgd_BC
#ifdef MOIST_CAPPA
        case ('cappa')
          bc_t0=>bc_side_t0%cappa_BC
          bc_t1=>bc_side_t1%cappa_BC
#endif
#ifdef USE_COND
        case ('q_con')
          bc_t0=>bc_side_t0%q_con_BC
          bc_t1=>bc_side_t1%q_con_BC
#endif
        case ('q')
          if(iq<1)then
            call mpp_error(FATAL,' iq<1 is not a valid index for q_BC array in retrieve_bc_variable_data')
          endif
          lbnd1=lbound(bc_side_t0%q_BC,1)
          lbnd2=lbound(bc_side_t0%q_BC,2)
          ubnd1=ubound(bc_side_t0%q_BC,1)
          ubnd2=ubound(bc_side_t0%q_BC,2)
          bc_t0=>bc_side_t0%q_BC(:,:,:,iq)
          bc_t1=>bc_side_t1%q_BC(:,:,:,iq)
        case ('u')
          bc_t0=>bc_side_t0%u_BC
          bc_t1=>bc_side_t1%u_BC
        case ('v')
          bc_t0=>bc_side_t0%v_BC
          bc_t1=>bc_side_t1%v_BC
        case ('uc')
          bc_t0=>bc_side_t0%uc_BC
          bc_t1=>bc_side_t1%uc_BC
        case ('vc')
          bc_t0=>bc_side_t0%vc_BC
          bc_t1=>bc_side_t1%vc_BC
!
      end select
!
      if(trim(bc_vbl_name)/='q')then
        lbnd1=lbound(bc_t0,1)
        lbnd2=lbound(bc_t0,2)
        ubnd1=ubound(bc_t0,1)
        ubnd2=ubound(bc_t0,2)
      endif
!
!---------------------------------------------------------------------
!
      end subroutine retrieve_bc_variable_data
!
!---------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------------------------------------------------------------------
!
      subroutine bc_time_interpolation(array                               &
                                      ,lbnd_x, ubnd_x                      &
                                      ,lbnd_y, ubnd_y                      &
                                      ,ubnd_z                              &
                                      ,bc_t0, bc_t1                        &
                                      ,lbnd1, ubnd1                        &
                                      ,lbnd2, ubnd2                        &
                                      ,i1,i2,j1,j2                         &
                                      ,is,ie,js,je                         &
                                      ,fcst_time                           &
                                      ,bc_update_interval                  &
                                      ,i1_blend,i2_blend,j1_blend,j2_blend &
				      ,i_bc,j_bc,nside,bc_vbl_name,blend )

!---------------------------------------------------------------------
!***  Update the boundary region of the input array at the given
!***  forecast time that is within the interval bracketed by the
!***  two current boundary region states.
!---------------------------------------------------------------------
      implicit none
!---------------------------------------------------------------------
!
!---------------------
!***  Input variables
!---------------------
!
      integer,intent(in) :: lbnd_x,ubnd_x,lbnd_y,ubnd_y,ubnd_z           !<-- Dimensions of the array to be updated.
!
      integer,intent(in) :: lbnd1,ubnd1,lbnd2,ubnd2                      !<-- Index limits of the BC arrays.
!
      integer,intent(in) :: i1,i2,j1,j2                               &  !<-- Index limits of the updated boundary region.
                           ,i_bc,j_bc                                 &  !<-- Innermost bndry indices (anchor pts for blending)
                           ,i1_blend,i2_blend,j1_blend,j2_blend       &  !<-- Index limits of the updated blending region.
                           ,nside
!
      integer,intent(in) :: is,ie,js,je                                  !<-- Min/Max index limits on task's computational subdomain
!
      integer,intent(in) :: bc_update_interval                           !<-- Time (hours) between BC data states
!
      real,intent(in) :: fcst_time                                       !<-- Current forecast time (sec)
!
      real,dimension(lbnd1:ubnd1,lbnd2:ubnd2,1:ubnd_z),intent(in) :: bc_t0  & !<-- Interpolate between these
                                                                    ,bc_t1    !    two boundary region states.
!
      character(len=*),intent(in) :: bc_vbl_name
!
      logical,intent(in) :: blend                                        !<-- Can blending be applied to this variable?
!
!---------------------
!*** Output variables
!---------------------
!
      real,dimension(lbnd_x:ubnd_x,lbnd_y:ubnd_y,1:ubnd_z)            &
                                               ,intent(out) :: array     !<-- Update boundary points in this array.
!
!---------------------
!***  Local variables
!---------------------
!
      integer :: i,j,k
!
      real :: blend_value,factor_dist,fraction_interval,rdenom
!
!---------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
!***  The current forecast time is this fraction of the way from
!***  time level 0 to time level 1.
!---------------------------------------------------------------------
!
      fraction_interval=mod(fcst_time,(bc_update_interval*3600.))     &
                       /(bc_update_interval*3600.)
!
!---------------------------------------------------------------------
!
      do k=1,ubnd_z
        do j=j1,j2
        do i=i1,i2
          array(i,j,k)=bc_t0(i,j,k)                                   &
                      +(bc_t1(i,j,k)-bc_t0(i,j,k))*fraction_interval
        enddo
        enddo
      enddo
!
!---------------------------------------------------------------------
!***  If this tracer is not in the external BC file then it will not
!***  be blended.
!---------------------------------------------------------------------
!
      if(.not.blend)then
        return
      endif
!
!---------------------------------------------------------------------
!***  Use specified external data to blend with integration values
!***  across nrows_blend rows immediately within the domain's
!***  boundary rows.  The weighting of the external data drops
!***  off exponentially.
!---------------------------------------------------------------------
!
!-----------
!***  North
!-----------
!
      if(nside==1.and.north_bc)then
        rdenom=1./real(j2_blend-j_bc-1)
        do k=1,ubnd_z
          do j=j1_blend,j2_blend
            factor_dist=exp(-(blend_exp1+blend_exp2*(j-j_bc-1)*rdenom)) !<-- Exponential falloff of blending weights.
            do i=i1_blend,i2_blend
              blend_value=bc_t0(i,j,k)                                &  !<-- Blend data interpolated 
                         +(bc_t1(i,j,k)-bc_t0(i,j,k))*fraction_interval  !    between t0 and t1.
!
              array(i,j,k)=(1.-factor_dist)*array(i,j,k)+factor_dist*blend_value
            enddo
          enddo
        enddo
      endif
!
!-----------
!***  South
!-----------
!
      if(nside==2.and.south_bc)then
        rdenom=1./real(j_bc-j1_blend-1)
        do k=1,ubnd_z
          do j=j1_blend,j2_blend
            factor_dist=exp(-(blend_exp1+blend_exp2*(j_bc-j-1)*rdenom)) !<-- Exponential falloff of blending weights.
            do i=i1_blend,i2_blend
              blend_value=bc_t0(i,j,k)                                &  !<-- Blend data interpolated 
                         +(bc_t1(i,j,k)-bc_t0(i,j,k))*fraction_interval  !    between t0 and t1.
              array(i,j,k)=(1.-factor_dist)*array(i,j,k)+factor_dist*blend_value
            enddo
          enddo
        enddo
      endif
!
!----------
!***  East
!----------
!
      if(nside==3.and.east_bc)then
        rdenom=1./real(i2_blend-i_bc-1)
        do k=1,ubnd_z
          do j=j1_blend,j2_blend
            do i=i1_blend,i2_blend
!
              blend_value=bc_t0(i,j,k)                                  &  !<-- Blend data interpolated 
                         +(bc_t1(i,j,k)-bc_t0(i,j,k))*fraction_interval    !    between t0 and t1.
!
              factor_dist=exp(-(blend_exp1+blend_exp2*(i-i_bc-1)*rdenom))  !<-- Exponential falloff of blending weights.
!
              array(i,j,k)=(1.-factor_dist)*array(i,j,k)+factor_dist*blend_value
            enddo
          enddo
        enddo
      endif
!
!----------
!***  West
!----------
!
      if(nside==4.and.west_bc)then
        rdenom=1./real(i_bc-i1_blend-1)
        do k=1,ubnd_z
          do j=j1_blend,j2_blend
            do i=i1_blend,i2_blend
!
              blend_value=bc_t0(i,j,k)                                  &  !<-- Blend data interpolated 
                         +(bc_t1(i,j,k)-bc_t0(i,j,k))*fraction_interval    !    between t0 and t1.
!
              factor_dist=exp(-(blend_exp1+blend_exp2*(i_bc-i-1)*rdenom))  !<-- Exponential falloff of blending weights.
!
              array(i,j,k)=(1.-factor_dist)*array(i,j,k)+factor_dist*blend_value
            enddo
          enddo
        enddo
      endif
!
!---------------------------------------------------------------------
!
      end subroutine bc_time_interpolation
!
!---------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------------------------------------------------------------------
!
      subroutine regional_bc_t1_to_t0(BC_t1,BC_t0                     &
                                     ,nlev,ntracers,bnds )
!
!---------------------------------------------------------------------
!***  BC data has been read into the time level t1 object.  Now
!***  move the t1 data into the t1 object before refilling t1
!***  with the next data from the BC file.
!---------------------------------------------------------------------
      implicit none
!---------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      integer,intent(in) :: nlev                                      &  !<-- # of model layers.
                           ,ntracers                                     !<-- # of advected tracers
!
      type(fv_regional_bc_bounds_type),intent(in) :: bnds                !<-- Index limits for all types of vbls in boundary region
!
      type(fv_domain_sides),target,intent(in) :: BC_t1
!
      type(fv_domain_sides),target,intent(inout) :: BC_t0
!
!---------------------
!***  Local variables
!---------------------
!
      integer :: i,ie_c,ie_s,ie_w,is_c,is_s,is_w                      &
                ,j,je_c,je_s,je_w,js_c,js_s,js_w                      &
                ,k,n,nside
!
      logical :: move
!
!---------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
!***  Loop through the four sides of the domain.
!---------------------------------------------------------------------
!
      sides: do nside=1,4
!
        move=.false.
!
        if(nside==1)then
          if(north_bc)then
            move=.true.
            bc_side_t0=>BC_t0%north
            bc_side_t1=>BC_t1%north
!
            is_c=bnds%is_north                                               !<--
            ie_c=bnds%ie_north                                               !  North BC index limits
            js_c=bnds%js_north                                               !  for centers of grid cells.
            je_c=bnds%je_north                                               !<--
!
            is_s=bnds%is_north_uvs                                           !<--
            ie_s=bnds%ie_north_uvs                                           !  North BC index limits
            js_s=bnds%js_north_uvs                                           !  for winds on N/S sides of grid cells.
            je_s=bnds%je_north_uvs                                           !<--
!
            is_w=bnds%is_north_uvw                                           !<--
            ie_w=bnds%ie_north_uvw                                           !  North BC index limits
            js_w=bnds%js_north_uvw                                           !  for winds on E/W sides of grid cells.
            je_w=bnds%je_north_uvw                                           !<--
          endif
        endif
!
        if(nside==2)then
          if(south_bc)then
            move=.true.
            bc_side_t0=>BC_t0%south
            bc_side_t1=>BC_t1%south
!
            is_c=bnds%is_south                                               !<--
            ie_c=bnds%ie_south                                               !  South BC index limits
            js_c=bnds%js_south                                               !  for centers of grid cells.
            je_c=bnds%je_south                                               !<--
!
            is_s=bnds%is_south_uvs                                           !<--
            ie_s=bnds%ie_south_uvs                                           !  South BC index limits
            js_s=bnds%js_south_uvs                                           !  for winds on N/S sides of grid cells.
            je_s=bnds%je_south_uvs                                           !<--
!
            is_w=bnds%is_south_uvw                                           !<--
            ie_w=bnds%ie_south_uvw                                           !  South BC index limits
            js_w=bnds%js_south_uvw                                           !  for winds on E/W sides of grid cells.
            je_w=bnds%je_south_uvw                                           !<--
          endif
        endif
!
        if(nside==3)then
          if(east_bc)then
            move=.true.
            bc_side_t0=>BC_t0%east
            bc_side_t1=>BC_t1%east
!
            is_c=bnds%is_east                                                !<--
            ie_c=bnds%ie_east                                                !  East BC index limits
            js_c=bnds%js_east                                                !  for centers of grid cells.
            je_c=bnds%je_east                                                !<--
!
            is_s=bnds%is_east_uvs                                            !<--
            ie_s=bnds%ie_east_uvs                                            !  East BC index limits
            js_s=bnds%js_east_uvs                                            !  for winds on N/S sides of grid cells.
            je_s=bnds%je_east_uvs                                            !<--
!
            is_w=bnds%is_east_uvw                                            !<--
            ie_w=bnds%ie_east_uvw                                            !  East BC index limits
            js_w=bnds%js_east_uvw                                            !  for winds on E/W sides of grid cells.
            je_w=bnds%je_east_uvw                                            !<--
          endif
        endif
!
        if(nside==4)then
          if(west_bc)then
            move=.true.
            bc_side_t0=>BC_t0%west
            bc_side_t1=>BC_t1%west
!
            is_c=bnds%is_west                                                !<--
            ie_c=bnds%ie_west                                                !  West BC index limits
            js_c=bnds%js_west                                                !  for centers of grid cells.
            je_c=bnds%je_west                                                !<--
!
            is_s=bnds%is_west_uvs                                            !<--
            ie_s=bnds%ie_west_uvs                                            !  West BC index limits
            js_s=bnds%js_west_uvs                                            !  for winds on N/S sides of grid cells.
            je_s=bnds%je_west_uvs                                            !<--
!
            is_w=bnds%is_west_uvw                                            !<--
            ie_w=bnds%ie_west_uvw                                            !  West BC index limits
            js_w=bnds%js_west_uvw                                            !  for winds on E/W sides of grid cells.
            je_w=bnds%je_west_uvw                                            !<--
          endif
        endif
!
        if(move)then
          do k=1,nlev
            do j=js_c,je_c
            do i=is_c,ie_c
              bc_side_t0%delp_BC(i,j,k) =bc_side_t1%delp_BC(i,j,k)
              bc_side_t0%divgd_BC(i,j,k)=bc_side_t1%divgd_BC(i,j,k)
            enddo
            enddo
          enddo
!
          do n=1,ntracers
            do k=1,nlev
              do j=js_c,je_c
              do i=is_c,ie_c
                bc_side_t0%q_BC(i,j,k,n)=bc_side_t1%q_BC(i,j,k,n)
              enddo
              enddo
            enddo
          enddo
!
          do k=1,nlev
            do j=js_c,je_c
            do i=is_c,ie_c
#ifndef SW_DYNAMICS
              bc_side_t0%w_BC(i,j,k)    =bc_side_t1%w_BC(i,j,k)
              bc_side_t0%pt_BC(i,j,k)   =bc_side_t1%pt_BC(i,j,k)
              bc_side_t0%delz_BC(i,j,k) =bc_side_t1%delz_BC(i,j,k)
#ifdef USE_COND
              bc_side_t0%q_con_BC(i,j,k)=bc_side_t1%q_con_BC(i,j,k)
#ifdef MOIST_CAPPA
              bc_side_t0%cappa_BC(i,j,k)=bc_side_t1%cappa_BC(i,j,k)
#endif
#endif
#endif
            enddo
            enddo
          enddo
!
          do k=1,nlev
            do j=js_s,je_s
            do i=is_s,ie_s
              bc_side_t0%u_BC(i,j,k) =bc_side_t1%u_BC(i,j,k)
              bc_side_t0%vc_BC(i,j,k)=bc_side_t1%vc_BC(i,j,k)
            enddo
            enddo
          enddo
!
          do k=1,nlev
            do j=js_w,je_w
            do i=is_w,ie_w
              bc_side_t0%v_BC(i,j,k) =bc_side_t1%v_BC(i,j,k)
              bc_side_t0%uc_BC(i,j,k)=bc_side_t1%uc_BC(i,j,k)
            enddo
            enddo
          enddo
!
        endif
!
      enddo sides
!
!---------------------------------------------------------------------
!
      end subroutine regional_bc_t1_to_t0
!
!---------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------------------------------------------------------------------
!
      subroutine convert_to_virt_pot_temp(isd,ied,jsd,jed,npz)
!
!-----------------------------------------------------------------------
!***  Convert the incoming sensible temperature to virtual potential
!***  temperature.
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      integer,intent(in) :: isd,ied,jsd,jed,npz
!
!---------------------
!***  Local variables
!---------------------
!
      integer :: i1,i2,j1,j2
!
      real :: rdg
!
      real,dimension(:,:,:),pointer :: delp,delz,pt
#ifdef USE_COND
      real,dimension(:,:,:),pointer :: q_con
#endif
#ifdef MOIST_CAPPA
      real,dimension(:,:,:),pointer ::cappa
#endif
!
      real,dimension(:,:,:,:),pointer :: q
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      if(.not.(north_bc.or.south_bc.or.east_bc.or.west_bc))then
        return
      endif
!
      rdg=-rdgas/grav
!
      if(north_bc)then
        i1=regional_bounds%is_north
        i2=regional_bounds%ie_north
        j1=regional_bounds%js_north
        j2=regional_bounds%je_north
        q    =>BC_t1%north%q_BC
#ifdef USE_COND
        q_con=>BC_t1%north%q_con_BC
#endif
        delp =>BC_t1%north%delp_BC
        delz =>BC_t1%north%delz_BC
#ifdef MOIST_CAPPA
        cappa=>BC_t1%north%cappa_BC
#endif
        pt   =>BC_t1%north%pt_BC
        call compute_vpt             !<-- Compute the virtual potential temperature.
      endif
!
      if(south_bc)then
        i1=regional_bounds%is_south
        i2=regional_bounds%ie_south
        j1=regional_bounds%js_south
        j2=regional_bounds%je_south
        q    =>BC_t1%south%q_BC
#ifdef USE_COND
        q_con=>BC_t1%south%q_con_BC
#endif
        delp =>BC_t1%south%delp_BC
        delz =>BC_t1%south%delz_BC
#ifdef MOIST_CAPPA
        cappa=>BC_t1%south%cappa_BC
#endif
        pt   =>BC_t1%south%pt_BC
        call compute_vpt             !<-- Compute the virtual potential temperature.
      endif
!
      if(east_bc)then
        i1=regional_bounds%is_east
        i2=regional_bounds%ie_east
        j1=regional_bounds%js_east
        j2=regional_bounds%je_east
        q    =>BC_t1%east%q_BC
#ifdef USE_COND
        q_con=>BC_t1%east%q_con_BC
#endif
        delp =>BC_t1%east%delp_BC
        delz =>BC_t1%east%delz_BC
#ifdef MOIST_CAPPA
        cappa=>BC_t1%east%cappa_BC
#endif
        pt   =>BC_t1%east%pt_BC
        call compute_vpt             !<-- Compute the virtual potential temperature.
      endif
!
      if(west_bc)then
        i1=regional_bounds%is_west
        i2=regional_bounds%ie_west
        j1=regional_bounds%js_west
        j2=regional_bounds%je_west
        q    =>BC_t1%west%q_BC
#ifdef USE_COND
        q_con=>BC_t1%west%q_con_BC
#endif
        delp =>BC_t1%west%delp_BC
        delz =>BC_t1%west%delz_BC
#ifdef MOIST_CAPPA
        cappa=>BC_t1%west%cappa_BC
#endif
        pt   =>BC_t1%west%pt_BC
        call compute_vpt             !<-- Compute the virtual potential temperature.
      endif
!
!-----------------------------------------------------------------------

      contains

!-----------------------------------------------------------------------
!
      subroutine compute_vpt
!
!-----------------------------------------------------------------------
!***  Compute the virtual potential temperature as done in fv_dynamics.
!-----------------------------------------------------------------------
!
!---------------------
!***  Local variables
!---------------------
!
      integer :: i,j,k
!
      real :: cvm,dp1,pkz
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      do k=1,npz
!
        do j=j1,j2
        do i=i1,i2
          dp1 = zvir*q(i,j,k,sphum_index)
#ifdef USE_COND
#ifdef MOIST_CAPPA
          cvm=(1.-q(i,j,k,sphum_index)+q_con(i,j,k))*cv_air             &
             +q(i,j,k,sphum_index)*cv_vap+q(i,j,k,liq_water_index)*c_liq
          pkz=exp(cappa(i,j,k)*log(rdg*delp(i,j,k)*pt(i,j,k)            &
              *(1.+dp1)*(1.-q_con(i,j,k))/delz(i,j,k)))
#else
          pkz=exp(kappa*log(rdg*delp(i,j,k)*pt(i,j,k)                   &
              *(1.+dp1)*(1.-q_con(i,j,k))/delz(i,j,k)))
#endif
          pt(i,j,k)=pt(i,j,k)*(1.+dp1)*(1.-q_con(i,j,k))/pkz
#else
          pkz=exp(kappa*log(rdg*delp(i,j,k)*pt(i,j,k)                   &
              *(1.+dp1)/delz(i,j,k)))
          pt(i,j,k)=pt(i,j,k)*(1.+dp1)/pkz          
#endif
        enddo
        enddo
!
      enddo
!
!-----------------------------------------------------------------------
!
      end subroutine compute_vpt
!
!-----------------------------------------------------------------------
!
      end subroutine convert_to_virt_pot_temp
!
!-----------------------------------------------------------------------
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------------------------------------------------------------------
!***  The following four subroutines are exact copies from
!***  external_ic_mod.  That module must USE this module therefore
!***  this module cannout USE external_IC_mod to get at those
!***  subroutines.  The routines may be moved to their own module.
!---------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------------------------------------------------------------------
 subroutine p_maxmin(qname, q, is, ie, js, je, km, fac)
      character(len=*), intent(in)::  qname
      integer, intent(in):: is, ie, js, je, km
      real, intent(in)::    q(is:ie, js:je, km)
      real, intent(in)::    fac
      real qmin, qmax
      integer i,j,k

      qmin = q(is,js,1)
      qmax = qmin
      do k=1,km
      do j=js,je
         do i=is,ie
            if( q(i,j,k) < qmin ) then
                qmin = q(i,j,k)
            elseif( q(i,j,k) > qmax ) then
                qmax = q(i,j,k)
            endif
          enddo
      enddo
      enddo
      call mp_reduce_min(qmin)
      call mp_reduce_max(qmax)
      if(is_master()) write(6,*) qname, qmax*fac, qmin*fac

 end subroutine p_maxmin


 subroutine pmaxmn(qname, q, is, ie, js, je, km, fac, area, domain)
      character(len=*), intent(in)::  qname
      integer, intent(in):: is, ie, js, je
      integer, intent(in):: km
      real, intent(in)::    q(is:ie, js:je, km)
      real, intent(in)::    fac
      real(kind=R_GRID), intent(IN)::  area(is-3:ie+3, js-3:je+3)
      type(domain2d), intent(INOUT) :: domain
!---local variables
      real qmin, qmax, gmean
      integer i,j,k

      qmin = q(is,js,1)
      qmax = qmin
      gmean = 0.

      do k=1,km
      do j=js,je
         do i=is,ie
            if( q(i,j,k) < qmin ) then
                qmin = q(i,j,k)
            elseif( q(i,j,k) > qmax ) then
                qmax = q(i,j,k)
            endif
          enddo
      enddo
      enddo

      call mp_reduce_min(qmin)
      call mp_reduce_max(qmax)

      gmean = g_sum(domain, q(is,js,km), is, ie, js, je, 3, area, 1, reproduce=.true.)
      if(is_master()) write(6,*) qname, qmax*fac, qmin*fac, gmean*fac

 end subroutine pmaxmn


 subroutine fillq(im, km, nq, q, dp)
   integer,  intent(in):: im                ! No. of longitudes
   integer,  intent(in):: km                ! No. of levels
   integer,  intent(in):: nq                ! Total number of tracers
   real , intent(in)::  dp(im,km)       ! pressure thickness
   real , intent(inout) :: q(im,km,nq)   ! tracer mixing ratio
! !LOCAL VARIABLES:
   integer i, k, ic, k1

   do ic=1,nq
! Bottom up:
      do k=km,2,-1
         k1 = k-1
         do i=1,im
           if( q(i,k,ic) < 0. ) then
               q(i,k1,ic) = q(i,k1,ic) + q(i,k,ic)*dp(i,k)/dp(i,k1)
               q(i,k ,ic) = 0.
           endif
         enddo
      enddo
! Top down:
      do k=1,km-1
         k1 = k+1
         do i=1,im
            if( q(i,k,ic) < 0. ) then
                q(i,k1,ic) = q(i,k1,ic) + q(i,k,ic)*dp(i,k)/dp(i,k1)
                q(i,k ,ic) = 0.
            endif
         enddo
      enddo

   enddo

 end subroutine fillq

 subroutine mp_auto_conversion(ql, qr, qi, qs)
 real, intent(inout):: ql, qr, qi, qs
 real, parameter:: qi0_max = 2.0e-3
 real, parameter:: ql0_max = 2.5e-3

! Convert excess cloud water into rain:
  if ( ql > ql0_max ) then
       qr = ql - ql0_max
       ql = ql0_max
  endif
! Convert excess cloud ice into snow:
  if ( qi > qi0_max ) then
       qs = qi - qi0_max
       qi = qi0_max
  endif

 end subroutine mp_auto_conversion

!-----------------------------------------------------------------------
!
      subroutine nudge_qv_bc(Atm,isd,ied,jsd,jed)
!
!-----------------------------------------------------------------------
!***  When nudging of specific humidity is selected then we must also
!***  nudge the values in the regional boundary.
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      integer,intent(in) :: isd,ied,jsd,jed                                !<-- Memory limits of task subdomain
!
      type(fv_atmos_type),target,intent(inout) :: Atm                      !<-- Atm object for the current domain
!
!---------------------
!***  Local variables
!---------------------
!
      integer :: i,i_x,ie,is,j,j_x,je,js,k
!
      real, parameter::    q1_h2o = 2.2E-6
      real, parameter::    q7_h2o = 3.8E-6
      real, parameter::  q100_h2o = 3.8E-6
      real, parameter:: q1000_h2o = 3.1E-6
      real, parameter:: q2000_h2o = 2.8E-6
      real, parameter:: q3000_h2o = 3.0E-6
      real, parameter:: wt=2., xt=1./(1.+wt)
!
      real :: p00,q00
!
      type(fv_regional_bc_bounds_type),pointer :: bnds
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      bnds=>Atm%regional_bc_bounds
!
!-----------
!***  North
!-----------
!
      if(north_bc)then
        is=lbound(BC_t1%north%q_BC,1)
        ie=ubound(BC_t1%north%q_BC,1)
        js=lbound(BC_t1%north%q_BC,2)
        je=ubound(BC_t1%north%q_BC,2)
!
        i_x=isd                                                          !<--  Use column at
        j_x=jsd                                                          !     this location.
!
        p00=Atm%ptop                                                     !<-- Use layer interface pressures.
!
        n_loopk: do k=1,npz
          if(p00<3000.)then                                              !<-- Apply nudging only if pressure < 30 mb. 
            call get_q00
            do j=js,je
            do i=is,ie
              BC_t1%north%q_BC(i,j,k,sphum_index)=                    &  !<-- Nudge the north boundary sphum at time t1.
                        xt*(BC_t1%north%q_BC(i,j,k,sphum_index)+wt*q00)
            enddo
            enddo
          else
            exit n_loopk
          endif
          p00=p00+BC_t1%north%delp_BC(i_x,j_x,k)
        enddo n_loopk
      endif
!
!-----------
!***  South
!-----------
!
      if(south_bc)then
        is=lbound(BC_t1%south%q_BC,1)
        ie=ubound(BC_t1%south%q_BC,1)
        js=lbound(BC_t1%south%q_BC,2)
        je=ubound(BC_t1%south%q_BC,2)
!
        i_x=isd                                                          !<--  Use column at
        j_x=jed                                                          !     this location.
!
        p00=Atm%ptop                                                     !<-- Use layer interface pressures.
!
        s_loopk: do k=1,npz      
          if(p00<3000.)then                                              !<-- Apply nudging only if pressure < 30 mb. 
            call get_q00
            do j=js,je
            do i=is,ie
              BC_t1%south%q_BC(i,j,k,sphum_index)=                    &  !<-- Nudge the south boundary sphum at time t1.
                        xt*(BC_t1%south%q_BC(i,j,k,sphum_index)+wt*q00)
            enddo
            enddo
          else
            exit s_loopk
          endif
          p00=p00+BC_t1%south%delp_BC(i_x,j_x,k)
        enddo s_loopk
      endif
!
!----------
!***  East
!----------
!
      if(east_bc)then
        is=lbound(BC_t1%east%q_BC,1)
        ie=ubound(BC_t1%east%q_BC,1)
        js=lbound(BC_t1%east%q_BC,2)
        je=ubound(BC_t1%east%q_BC,2)
!
        i_x=isd                                                          !<--  Use column at
        j_x=jsd+nhalo_model                                              !     this location.
!
        p00=Atm%ptop                                                     !<-- Use layer interface pressures.
!
        e_loopk: do k=1,npz      
          if(p00<3000.)then                                              !<-- Apply nudging only if pressure < 30 mb. 
            call get_q00
            do j=js,je
            do i=is,ie
              BC_t1%east%q_BC(i,j,k,sphum_index)=                     &  !<-- Nudge the east boundary sphum at time t1.
                        xt*(BC_t1%east%q_BC(i,j,k,sphum_index)+wt*q00)
            enddo
            enddo
          else
            exit e_loopk
          endif
          p00=p00+BC_t1%east%delp_BC(i_x,j_x,k)
        enddo e_loopk
      endif
!
!----------
!***  West
!----------
!
      if(west_bc)then
        is=lbound(BC_t1%west%q_BC,1)
        ie=ubound(BC_t1%west%q_BC,1)
        js=lbound(BC_t1%west%q_BC,2)
        je=ubound(BC_t1%west%q_BC,2)
!
        i_x=ied                                                          !<--  Use column at
        j_x=jsd+nhalo_model                                              !     this location.
!
        p00=Atm%ptop                                                     !<-- Use layer interface pressures.
!
        w_loopk: do k=1,npz      
          if(p00<3000.)then                                              !<-- Apply nudging only if pressure < 30 mb. 
            call get_q00
            do j=js,je
            do i=is,ie
              BC_t1%west%q_BC(i,j,k,sphum_index)=                     &  !<-- Nudge the west boundary sphum at time t1.
                        xt*(BC_t1%west%q_BC(i,j,k,sphum_index)+wt*q00)
            enddo
            enddo
          else
            exit w_loopk
          endif
          p00=p00+BC_t1%west%delp_BC(i_x,j_x,k)
        enddo w_loopk
      endif
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine get_q00
!
!-----------------------------------------------------------------------
!***  This is an internal subroutine to subroutine nudge_qv_bc that
!***  computes the climatological contribution to the nudging ot the
!***  input specific humidity.
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      if ( p00 < 30.E2 ) then
        if ( p00 < 1. ) then
          q00 = q1_h2o
        elseif ( p00 <= 7. .and. p00 >= 1. ) then
          q00 = q1_h2o + (q7_h2o-q1_h2o)*log(pref(k)/1.)/log(7.)
        elseif ( p00 < 100. .and. p00 >= 7. ) then
          q00 = q7_h2o + (q100_h2o-q7_h2o)*log(pref(k)/7.)/log(100./7.)
        elseif ( p00 < 1000. .and. p00 >= 100. ) then
          q00 = q100_h2o + (q1000_h2o-q100_h2o)*log(pref(k)/1.E2)/log(10.)
        elseif ( p00 < 2000. .and. p00 >= 1000. ) then
          q00 = q1000_h2o + (q2000_h2o-q1000_h2o)*log(pref(k)/1.E3)/log(2.)
        else
          q00 = q2000_h2o + (q3000_h2o-q2000_h2o)*log(pref(k)/2.E3)/log(1.5)
        endif
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine get_q00
!
!-----------------------------------------------------------------------
!
      end subroutine nudge_qv_bc
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

  subroutine dump_field_3d (domain, name, field, isd, ied, jsd, jed, nlev, stag)

!-----------------------------------------------------------------------
!***  Subroutines dump_field_2d and dump_field_3d are module 
!***  procedures with the generic interface 'dump_field'.
!***  Use these routines to write out NetCDF files containing 
!***  FULL fields that include the variables' boundary region.
!***  See the following four examples for guidance on how to 
!***  call the routines.
!-----------------------------------------------------------------------
!   call dump_field(Atm(1)%domain,"atm_pt",   Atm(1)%pt,   isd, ied,   jsd, jed,   Atm(1)%npz, stag=H_STAGGER)
!   call dump_field(Atm(1)%domain,"atm_u",    Atm(1)%u,    isd, ied,   jsd, jed+1, Atm(1)%npz, stag=U_STAGGER)
!   call dump_field(Atm(1)%domain,"atm_v",    Atm(1)%v,    isd, ied+1, jsd, jed,   Atm(1)%npz, stag=V_STAGGER)
!   call dump_field(Atm(1)%domain,"atm_phis", Atm(1)%phis, isd, ied,   jsd, jed,               stag=H_STAGGER)

    type(domain2d),         intent(INOUT) :: domain
    character(len=*),       intent(IN)    :: name
    real, dimension(isd:ied,jsd:jed,1:nlev), intent(INOUT) :: field
    integer,                intent(IN)    :: isd, ied, jsd, jed, nlev
    integer,                intent(IN)    :: stag

    integer                             :: unit
    character(len=128)                  :: fname
    type(axistype)                      :: x, y, z
    type(fieldtype)                     :: f
    type(domain1D)                      :: xdom, ydom
    integer                             :: nz
    integer                             :: is, ie, js, je
    integer                             :: isg, ieg, jsg, jeg, nxg, nyg, npx, npy
    integer                             :: i, j, halo, iext, jext
    logical                             :: is_root_pe
    real, allocatable, dimension(:,:,:) :: glob_field
    integer, allocatable, dimension(:)  :: pelist
    character(len=1)                    :: stagname
    integer                             :: isection_s, isection_e, jsection_s, jsection_e

    write(fname,"(A,A,A,I1.1,A)") "regional_",name,".tile", 7 , ".nc"
    write(0,*)'dump_field_3d: file name = |', trim(fname) , '|'

    call mpp_get_domain_components( domain, xdom, ydom )
    call mpp_get_compute_domain( domain, is,  ie,  js,  je )
    call mpp_get_global_domain ( domain, isg, ieg, jsg, jeg, xsize=npx, ysize=npy, position=CENTER )

    halo = is - isd
    if ( halo /= 3 ) then
       write(0,*) 'dusan- halo should be 3 ', halo
    endif

    iext = 0
    jext = 0
    stagname = "h";
    if (stag == U_STAGGER) then
      jext = 1
      stagname = "u";
    endif
    if (stag == V_STAGGER) then
      iext = 1
      stagname = "v";
    endif

    nxg = npx + 2*halo + iext
    nyg = npy + 2*halo + jext
    nz = size(field,dim=3)

    allocate( glob_field(isg-halo:ieg+halo+iext, jsg-halo:jeg+halo+jext, 1:nz) )

    isection_s = is
    isection_e = ie
    jsection_s = js
    jsection_e = je

    if ( isd < 0 )     isection_s = isd
    if ( ied > npx-1 ) isection_e = ied
    if ( jsd < 0 )     jsection_s = jsd
    if ( jed > npy-1 ) jsection_e = jed

    allocate( pelist(mpp_npes()) )
    call mpp_get_current_pelist(pelist)

    is_root_pe = (mpp_pe()==mpp_root_pe())

    call mpp_gather(isection_s,isection_e,jsection_s,jsection_e, nz, &
                    pelist, field(isection_s:isection_e,jsection_s:jsection_e,:), glob_field, is_root_pe, halo, halo)

    call mpp_open( unit, trim(fname), action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_SINGLE)

    call mpp_write_meta( unit, x, 'grid_xt', 'km', 'X distance', 'X', domain=xdom, data=(/(i*1.0,i=1,nxg)/) )
    call mpp_write_meta( unit, y, 'grid_yt', 'km', 'Y distance', 'Y', domain=ydom, data=(/(j*1.0,j=1,nyg)/) )
    call mpp_write_meta( unit, z, 'lev',     'km', 'Z distance',                   data=(/(i*1.0,i=1,nz)/) )

    call mpp_write_meta( unit, f, (/x,y,z/), name, 'unit', name)
    call mpp_write_meta( unit, "stretch_factor", rval=stretch_factor )
    call mpp_write_meta( unit, "target_lon", rval=target_lon )
    call mpp_write_meta( unit, "target_lat", rval=target_lat )
    call mpp_write_meta( unit, "cube_res", ival= cube_res)
    call mpp_write_meta( unit, "parent_tile", ival=parent_tile )
    call mpp_write_meta( unit, "refine_ratio", ival=refine_ratio )
    call mpp_write_meta( unit, "istart_nest", ival=istart_nest )
    call mpp_write_meta( unit, "jstart_nest", ival=jstart_nest )
    call mpp_write_meta( unit, "iend_nest", ival=iend_nest )
    call mpp_write_meta( unit, "jend_nest", ival=jend_nest )
    call mpp_write_meta( unit, "ihalo_shift", ival=halo )
    call mpp_write_meta( unit, "jhalo_shift", ival=halo )
    call mpp_write_meta( unit, mpp_get_id(f), "hstagger", cval=stagname )
    call mpp_write( unit, x )
    call mpp_write( unit, y )
    call mpp_write( unit, z )
    call mpp_write( unit, f, glob_field )

    call mpp_close( unit )

  end subroutine dump_field_3d

  subroutine dump_field_2d (domain, name, field, isd, ied, jsd, jed, stag)

    type(domain2d),         intent(INOUT) :: domain
    character(len=*),       intent(IN)    :: name
    real, dimension(isd:ied,jsd:jed), intent(INOUT) :: field
    integer,                intent(IN)    :: isd, ied, jsd, jed
    integer,                intent(IN)    :: stag

    integer                             :: unit
    character(len=128)                  :: fname
    type(axistype)                      :: x, y
    type(fieldtype)                     :: f
    type(domain1D)                      :: xdom, ydom
    integer                             :: is, ie, js, je
    integer                             :: isg, ieg, jsg, jeg, nxg, nyg, npx, npy
    integer                             :: i, j, halo, iext, jext
    logical                             :: is_root_pe
    real, allocatable, dimension(:,:)   :: glob_field
    integer, allocatable, dimension(:)  :: pelist
    character(len=1)                    :: stagname
    integer                             :: isection_s, isection_e, jsection_s, jsection_e

    write(fname,"(A,A,A,I1.1,A)") "regional_",name,".tile", 7 , ".nc"
    write(0,*)'dump_field_3d: file name = |', trim(fname) , '|'

    call mpp_get_domain_components( domain, xdom, ydom )
    call mpp_get_compute_domain( domain, is,  ie,  js,  je )
    call mpp_get_global_domain ( domain, isg, ieg, jsg, jeg, xsize=npx, ysize=npy, position=CENTER )

    halo = is - isd
    if ( halo /= 3 ) then
       write(0,*) 'dusan- halo should be 3 ', halo
    endif

    iext = 0
    jext = 0
    stagname = "h";
    if (stag == U_STAGGER) then
      jext = 1
      stagname = "u";
    endif
    if (stag == V_STAGGER) then
      iext = 1
      stagname = "v";
    endif

    nxg = npx + 2*halo + iext
    nyg = npy + 2*halo + jext

    allocate( glob_field(isg-halo:ieg+halo+iext, jsg-halo:jeg+halo+jext) )

    isection_s = is
    isection_e = ie
    jsection_s = js
    jsection_e = je

    if ( isd < 0 )     isection_s = isd
    if ( ied > npx-1 ) isection_e = ied
    if ( jsd < 0 )     jsection_s = jsd
    if ( jed > npy-1 ) jsection_e = jed

    allocate( pelist(mpp_npes()) )
    call mpp_get_current_pelist(pelist)

    is_root_pe = (mpp_pe()==mpp_root_pe())

    call mpp_gather(isection_s,isection_e,jsection_s,jsection_e, &
                    pelist, field(isection_s:isection_e,jsection_s:jsection_e), glob_field, is_root_pe, halo, halo)

    call mpp_open( unit, trim(fname), action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_SINGLE)

    call mpp_write_meta( unit, x, 'grid_xt', 'km', 'X distance', 'X', domain=xdom, data=(/(i*1.0,i=1,nxg)/) )
    call mpp_write_meta( unit, y, 'grid_yt', 'km', 'Y distance', 'Y', domain=ydom, data=(/(j*1.0,j=1,nyg)/) )

    call mpp_write_meta( unit, f, (/x,y/), name, 'unit', name)
    call mpp_write_meta( unit, "stretch_factor", rval=stretch_factor )
    call mpp_write_meta( unit, "target_lon", rval=target_lon )
    call mpp_write_meta( unit, "target_lat", rval=target_lat )
    call mpp_write_meta( unit, "cube_res", ival= cube_res)
    call mpp_write_meta( unit, "parent_tile", ival=parent_tile )
    call mpp_write_meta( unit, "refine_ratio", ival=refine_ratio )
    call mpp_write_meta( unit, "istart_nest", ival=istart_nest )
    call mpp_write_meta( unit, "jstart_nest", ival=jstart_nest )
    call mpp_write_meta( unit, "iend_nest", ival=iend_nest )
    call mpp_write_meta( unit, "jend_nest", ival=jend_nest )
    call mpp_write_meta( unit, "ihalo_shift", ival=halo )
    call mpp_write_meta( unit, "jhalo_shift", ival=halo )
    call mpp_write_meta( unit, mpp_get_id(f), "hstagger", cval=stagname )
    call mpp_write( unit, x )
    call mpp_write( unit, y )
    call mpp_write( unit, f, glob_field )

    call mpp_close( unit )

  end subroutine dump_field_2d

!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      subroutine prepare_full_fields(Atm)
!
!-----------------------------------------------------------------------
!***  Prepare the objects that will hold the names and values of
!***  the core and tracer fields to be written into the expanded
!***  restart files that include the boundary rows so the GSI
!***  can update both the interior and BCs.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      type(fv_atmos_type),target,intent(inout) :: Atm                      !<-- Atm object for the current domain
!
!---------------------
!***  Local variables
!---------------------
!
      integer :: index,istat,n                                          &
                ,ncid_core_new                                          &
                ,ncid_tracers_new                                       &
                ,ndims,nkount,nv_core,nv_tracers                        &
                ,var_id
!
      integer :: lbnd1,lbnd2,lbnd3,ubnd1,ubnd2,ubnd3
!
      integer,dimension(ndims_core) :: dim_lengths_core
!
      integer,dimension(ndims_tracers) :: dim_lengths_tracers
!
      integer,dimension(1:4) :: dimids=(/0,0,0,0/)
!
      real,dimension(:),allocatable :: dim_values
!
      character(len=50) :: att_name,var_name
!
      character(len=9),dimension(ndims_core) :: dim_names_core=(/           &
                                                                 'xaxis_1'  &
                                                                ,'xaxis_2'  &
                                                                ,'yaxis_1'  &
                                                                ,'yaxis_2'  &
                                                                ,'zaxis_1'  &
                                                                ,'Time   '  &
                                                                /)
!
      character(len=9),dimension(ndims_tracers) :: dim_names_tracers=(/           &
                                                                       'xaxis_1'  &
                                                                      ,'yaxis_1'  &
                                                                      ,'zaxis_1'  &
                                                                      ,'Time   '  &
                                                                      /)
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  The first file to be handled is the core restart file.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  All tasks are given pointers into the model data that will
!***  be written to the new restart file.  The following are the
!***  prognostic variables in the core restart file.  Note that
!***  we must add the halo region back into DZ since we need the
!***  domain boundary points for all the fields.
!-----------------------------------------------------------------------
!
      allocate(fields_core(1:nvars_core))
!
      fields_core(1)%ptr=>Atm%u
      fields_core(1)%name='u'
!
      fields_core(2)%ptr=>Atm%v
      fields_core(2)%name='v'
!
      fields_core(3)%ptr=>Atm%w
      fields_core(3)%name='W'
!
      lbnd1=lbound(Atm%delz,1)
      ubnd1=ubound(Atm%delz,1)
      lbnd2=lbound(Atm%delz,2)
      ubnd2=ubound(Atm%delz,2)
      lbnd3=1
      ubnd3=ubound(Atm%delz,3)
      allocate(fields_core(4)%ptr(lbnd1-nhalo_model:ubnd1+nhalo_model   &
                                 ,lbnd2-nhalo_model:ubnd2+nhalo_model   &
                                 ,lbnd3:ubnd3))
      fields_core(4)%name='DZ'
!
      fields_core(5)%ptr=>Atm%pt
      fields_core(5)%name='T'
!
      fields_core(6)%ptr=>Atm%delp
      fields_core(6)%name='delp'
!
      allocate(fields_core(7)%ptr(lbound(Atm%phis,1):ubound(Atm%phis,1) &
                            ,lbound(Atm%phis,2):ubound(Atm%phis,2)      &
                            ,1:1))
      fields_core(7)%ptr(:,:,1)=Atm%phis(:,:)                              !<-- For generality treat the 2-D phis as 3-D
      fields_core(7)%name='phis'
!
!-----------------------------------------------------------------------
!***  We need to point at the tracers in the model's tracer array.
!***  Those tracers depend on the physics that was selected so they
!***  cannot be pre-specified like the variables in the core restart
!***  file were.  Read them from the expanded tracer restart file
!***  that was created prior to the start for the forecast.
!-----------------------------------------------------------------------
!
      call check(nf90_open(path=filename_tracers_new                    &  !<-- The expanded tracer restart file.
                          ,mode=nf90_nowrite                            &  !<-- File access.
                          ,ncid=ncid_tracers_new ))                        !<-- The expanded tracer restart file's ID
!
      call check(nf90_inquire(ncid      =ncid_tracers_new               &  !<-- The expanded tracer restart file's ID.
                             ,nvariables=nv_tracers       ))               !<-- The TOTAL number of tracer restart file variables.
!
      nfields_tracers=nv_tracers-ndims_tracers                             !<-- # of 3-D tracer fields
      allocate(fields_tracers(1:nfields_tracers),stat=istat)
      if(istat/=0)then
        call mpp_error(FATAL,' Failed to allocate fields_tracers.')
      else
        if(is_master())then
          write(0,33012)nfields_tracers
33012     format(' Allocated fields_tracers(1:',i3,')')
        endif
      endif
      nkount=0
!
      do n=1,nv_tracers
        var_id=n
        call check(nf90_inquire_variable(ncid =ncid_tracers_new         &  !<-- The file's ID.
                                        ,varid=var_id                   &  !<-- The variable's ID.
                                        ,name =var_name ))                 !<-- The variable's name.
!
        if(n>ndims_tracers)then
          nkount=nkount+1
          fields_tracers(nkount)%name=trim(var_name)
          index=get_tracer_index(MODEL_ATMOS, trim(var_name))
          fields_tracers(nkount)%ptr=>Atm%q(:,:,:, index)
        endif
!
      enddo
!
!-----------------------------------------------------------------------
!
      call check(nf90_close(ncid_tracers_new))
!
!-----------------------------------------------------------------------
!
      end subroutine prepare_full_fields
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      subroutine write_full_fields(Atm)
!
!-----------------------------------------------------------------------
!***  Write out full fields of the primary restart variables
!***  INCLUDING BOUNDARY ROWS so the GSI can include BCs in its
!***  update.  This is done in a restart look-alike file.
!-----------------------------------------------------------------------
!
      type(fv_atmos_type), intent(inout), target :: Atm
!
      integer :: count_i,count_j
      integer :: iend,istart,jend,jstart,kend,kstart,nz
      integer :: iend_ptr,istart_ptr,jend_ptr,jstart_ptr
      integer :: iend_g,istart_g,jend_g,jstart_g
      integer :: ieg,iext,isg,jeg,jext,jsg,k
      integer :: n,ncid_core_new,ncid_tracers_new,nv,var_id
      integer :: halo
!
      integer,dimension(:),allocatable :: pelist
!
      real,dimension(:,:,:),allocatable :: global_field
      real,dimension(:,:,:),pointer :: field_3d
!
      character(len=10) :: var_name
!
      logical :: is_root_pe
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      allocate( pelist(mpp_npes()) )
      call mpp_get_current_pelist(pelist)
!     write(0,*)' pelist=',pelist
!
      halo=nhalo_model
!
      is_root_pe = (mpp_pe()==mpp_root_pe())
      if(is_root_pe)then
        call check(nf90_open(filename_core_new,nf90_write,ncid_core_new))  !<-- Open the new netcdf file
        write(0,*)' Opened core restart with BCs: ',trim(filename_core_new)
      endif
!
!-----------------------------------------------------------------------
!***  Save the global limits of the domain and its vertical extent.
!-----------------------------------------------------------------------
!
      call mpp_get_global_domain (Atm%domain, isg, ieg, jsg, jeg, position=CENTER )
!
!-----------------------------------------------------------------------
!***  Begin with the core restart file.
!***  Loop through that file's prognostic variables.
!-----------------------------------------------------------------------
!
      vbls_core: do nv=1,nvars_core
!
        var_name=fields_core(nv)%name
        if(is_root_pe)then
          call check(nf90_inq_varid(ncid_core_new,var_name,var_id))        !<-- Get this variable's ID
        endif
!
!-----------------------------------------------------------------------
!***  What is the full domain extent of this variable including 
!***  boundary rows?
!-----------------------------------------------------------------------
!
        iext=0
        jext=0
        if(var_name=='u'.or.var_name=='vc')then
          jext=1
        endif
        if(var_name=='v'.or.var_name=='uc')then
          iext=1
        endif
!
        call mpp_get_global_domain (atm%domain, isg, ieg, jsg, jeg, position=CENTER )
        istart_g=isg-halo
        iend_g  =ieg+halo+iext
        jstart_g=jsg-halo
        jend_g  =jeg+halo+jext
!
        count_i=iend_g-istart_g+1
        count_j=jend_g-jstart_g+1
!
        nz=size(fields_core(nv)%ptr,3)
!
        allocate( global_field(istart_g:iend_g, jstart_g:jend_g, 1:nz) )
!
!-----------------------------------------------------------------------
!***  What is the local extent of the variable on the task subdomain? 
!***  We must exclude inner halo data since the data is not updated
!***  there in some of the variables.  Of course the outer halo data 
!***  around the domain boundary is included.
!-----------------------------------------------------------------------
!
        istart=lbound(fields_core(nv)%ptr,1)   
        if(istart>1)then
          istart=istart+halo
        endif
!
        iend  =ubound(fields_core(nv)%ptr,1)  
        if(iend<ieg-halo)then
          iend=iend-halo
        endif
!
        jstart=lbound(fields_core(nv)%ptr,2)   
        if(jstart>1)then
          jstart=jstart+halo
        endif
!
        jend  =ubound(fields_core(nv)%ptr,2)  
        if(jend<jeg-halo)then
          jend=jend-halo
        endif
!
!-----------------------------------------------------------------------
!***  The interior values of the pt array are the sensible
!***  temperature.  The halo points though remain as the 
!***  special potential temperature used inside the dynamics 
!***  since those halo values never needed to be converted
!***  back to sensible.  We are now writing out the full
!***  field including boundary rows for the GSI so the domain
!***  halo points need to be converted to sensible temperature.
!***  Each boundary task will now do the conversion before the
!***  values are gathered onto the root task for writing out.
!***  Also since the DZ array no longer has halo points we must
!***  insert values back into the domain boundary rows of the
!***  object holding DZ.
!-----------------------------------------------------------------------
!
        if(trim(fields_core(nv)%name)=='T'                              &
                        .or.                                            &
           trim(fields_core(nv)%name)=='DZ') then
!
          call apply_delz_boundary(istart,iend,jstart,jend,nz           &
                                  ,Atm                                  &
                                  ,fields_core(nv)%name                 &
                                  ,fields_core(nv)%ptr(istart:iend,jstart:jend,:))
        endif
!
!-----------------------------------------------------------------------
!***  Gather onto a single task one layer at a time.  That task
!***  writes the full data to the new larger restart file.
!-----------------------------------------------------------------------
!
        do k=1,nz
          call mpp_gather(istart,iend,jstart,jend                                &
                         ,pelist, fields_core(nv)%ptr(istart:iend,jstart:jend,k) &
                         ,global_field(:,:,k), is_root_pe, halo, halo)
!
          if(is_root_pe)then
            call check(nf90_put_var(ncid_core_new,var_id                         &
                                   ,global_field(:,:,k)                          &
                                   ,start=(/1,1,k/)                              &
                                   ,count=(/count_i,count_j,1/)))
          endif 
        enddo
!
        deallocate(global_field)
!
      enddo vbls_core
!
      if(is_root_pe)then
        call check(nf90_close(ncid_core_new))
      endif
!
!-----------------------------------------------------------------------
!***  Now open the new tracer restart file.
!-----------------------------------------------------------------------
!
      if(is_root_pe)then
        call check(nf90_open(filename_tracers_new,nf90_write,ncid_tracers_new))  !<-- Open the new netcdf file
        write(0,*)' Opened tracer restart with BCs: ',trim(filename_tracers_new)
      endif
!
!-----------------------------------------------------------------------
!***  What is the full domain extent of this variable including 
!***  boundary rows?
!-----------------------------------------------------------------------
!
      call mpp_get_global_domain (Atm%domain, isg, ieg, jsg, jeg, position=CENTER )
      istart_g=isg-halo
      iend_g  =ieg+halo
      jstart_g=jsg-halo
      jend_g  =jeg+halo
!
      count_i=iend_g-istart_g+1
      count_j=jend_g-jstart_g+1
      nz=size(fields_tracers(1)%ptr,3)
!
      allocate( global_field(istart_g:iend_g, jstart_g:jend_g, 1:nz) )
!
!-----------------------------------------------------------------------
!***  What is the local extent of the variable on the task subdomain? 
!***  We must exclude inner halo data since the data is not updated
!***  there in some of the variables.  Of course the outer halo data 
!***  around the domain boundary is included.  These values are the
!***  same for all the tracers.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  The following are local bounds based on the global indexing.
!-----------------------------------------------------------------------
!
      lbnd_x_tracers=lbound(Atm%q,1)
      ubnd_x_tracers=ubound(Atm%q,1)
      lbnd_y_tracers=lbound(Atm%q,2)
      ubnd_y_tracers=ubound(Atm%q,2)
!
      istart=lbnd_x_tracers
      if(istart>1)then
        istart=istart+halo
      endif
!
      iend  =ubnd_x_tracers
      if(iend<ieg-halo)then
        iend=iend-halo
      endif
!
      jstart=lbnd_y_tracers
      if(jstart>1)then
        jstart=jstart+halo
      endif
!
      jend  =ubnd_y_tracers
      if(jend<jeg-halo)then
        jend=jend-halo
      endif
!
!-----------------------------------------------------------------------
!***  The following are local bounds based on the indexing of
!***  the pointers to the tracer arrays in memory.  They are all
!***  relative to 1 since each task pointed at the arrays as
!***  ptr => Atm%q(:,:,:,sphum_index) rather than as was done
!***  for the core arrays which was ptr => Atm%u .
!-----------------------------------------------------------------------
!
      istart_ptr=halo+1
      iend_ptr  =ubnd_x_tracers-lbnd_x_tracers+1-halo
      jstart_ptr=halo+1
      jend_ptr  =ubnd_y_tracers-lbnd_y_tracers+1-halo
!
      if(north_bc)then
        jstart_ptr=1
      endif
      if(south_bc)then
        jend_ptr=ubnd_y_tracers-lbnd_y_tracers+1
      endif
      if(east_bc)then
        istart_ptr=1
      endif
      if(west_bc)then
        iend_ptr=ubnd_x_tracers-lbnd_x_tracers+1
      endif
!
!-----------------------------------------------------------------------
!***  Loop through that file's prognostic tracers.
!-----------------------------------------------------------------------
!
      vbls_tracers: do nv=1,nfields_tracers
!
        var_name=fields_tracers(nv)%name
        if(is_root_pe)then
          call check(nf90_inq_varid(ncid_tracers_new,var_name,var_id))     !<-- Get this variable's ID
        endif
!
!-----------------------------------------------------------------------
!***  Gather onto a single task one layer at a time.  That task
!***  writes the full data to the new larger restart file.
!-----------------------------------------------------------------------
!
        do k=1,nz
          call mpp_gather(istart,iend,jstart,jend                                    &
                         ,pelist, fields_tracers(nv)%ptr(istart_ptr:iend_ptr,jstart_ptr:jend_ptr,k)  &
                         ,global_field(:,:,k), is_root_pe, halo, halo)
!
          if(is_root_pe)then
            call check(nf90_put_var(ncid_tracers_new,var_id             &
                                   ,global_field(:,:,k)                 &
                                   ,start=(/1,1,k/)                     &
                                   ,count=(/count_i,count_j,1/)))
          endif 
        enddo
!
      enddo vbls_tracers
!
      deallocate(global_field)
!
      if(is_root_pe)then
        call check(nf90_close(ncid_tracers_new))
      endif
!
!---------------------------------------------------------------------
      end subroutine write_full_fields
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!
      subroutine apply_delz_boundary(istart,iend,jstart,jend,nz       &
                                    ,Atm                              &
                                    ,name                             &
                                    ,field)
!
!---------------------------------------------------------------------
!***  Use the current boundary values of delz to convert the
!***  boundary potential temperature to sensible temperature
!***  and to fill in the boundary rows of the 3D delz array.
!---------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      integer,intent(in) :: istart,iend,jstart,jend
!
      character(len=*),intent(in) :: name
!
      type(fv_atmos_type),intent(inout) :: Atm
!
      real,dimension(istart:iend,jstart:jend,1:nz),intent(inout) :: field
!
!---------------------
!***  Local variables
!---------------------
!
      integer :: i1,i2,j1,j2,nz
      integer :: lbnd1,lbnd2,ubnd1,ubnd2,i,j,k
!
      real :: rdg
!
      real,dimension(:,:,:),pointer :: delz_ptr
!
!---------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------
!
!fill interior delz points before dealing with boundaries
!
!---------------------------------------------------------------------
!***  Fill the interior of the full delz array using Atm%delz
!***  which does not have a boundary.
!---------------------------------------------------------------------
!
      if (trim(name)=='DZ') then
        lbnd1=lbound(Atm%delz,1)
        ubnd1=ubound(Atm%delz,1)
        lbnd2=lbound(Atm%delz,2)
        ubnd2=ubound(Atm%delz,2)
!
        do k=1,nz
        do j=lbnd2,ubnd2
        do i=lbnd1,ubnd1
          field(i,j,k)=Atm%delz(i,j,k)
        enddo
        enddo
        enddo
      endif
!
      if(.not.(north_bc.or.south_bc.or.east_bc.or.west_bc))then
        return                                                           !<-- Tasks not on the boundary may exit.
      endif
!
      rdg=-rdgas/grav
!
!---------------------------------------------------------------------
!
      if(north_bc)then
        i1=istart
        i2=iend
        j1=jstart
        j2=jstart+nhalo_model-1
        delz_ptr=>delz_auxiliary%north
!
        if(trim(name)=='T')then
          call compute_halo_t
        elseif(trim(name)=='DZ')then
          call fill_delz
        endif
!
      endif
!
      if(south_bc)then
        i1=istart
        i2=iend
        j1=jend-nhalo_model+1
        j2=jend
        delz_ptr=>delz_auxiliary%south
!
        if(trim(name)=='T')then
          call compute_halo_t
        elseif(trim(name)=='DZ')then
          call fill_delz
        endif
!
      endif
!
      if(east_bc)then
        i1=istart
        i2=istart+nhalo_model-1
        j1=jstart
        j2=jend
        if(north_bc)then
          j1=jstart+nhalo_model
        elseif(south_bc)then
          j2=jend-nhalo_model
        endif
        delz_ptr=>delz_auxiliary%east
!
        if(trim(name)=='T')then
          call compute_halo_t
        elseif(trim(name)=='DZ')then
          call fill_delz
        endif
!
      endif
!
      if(west_bc)then
        i1=iend-nhalo_model+1
        i2=iend
        j1=jstart
        j2=jend
        if(north_bc)then
          j1=jstart+nhalo_model
        elseif(south_bc)then
          j2=jend-nhalo_model
        endif
        delz_ptr=>delz_auxiliary%west
!
        if(trim(name)=='T')then
          call compute_halo_t
        elseif(trim(name)=='DZ')then
          call fill_delz
        endif
!
      endif
!
!---------------------------------------------------------------------
      contains
!---------------------------------------------------------------------
!
      subroutine compute_halo_t
!
!---------------------------------------------------------------------
!
      integer :: i,j,k
!
      real :: cappa,cvm,dp1,part1,part2
!
!---------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------
!
      do k=1,nz
      do j=j1,j2
      do i=i1,i2
        dp1 = zvir*Atm%q(i,j,k,sphum_index)
        cvm=(1.-Atm%q(i,j,k,sphum_index)+Atm%q_con(i,j,k))*cv_air     &
            +Atm%q(i,j,k,sphum_index)*cv_vap                          &
            +Atm%q(i,j,k,liq_water_index)*c_liq
        cappa=rdgas/(rdgas+cvm/(1.+dp1))
!
        part1=(1.+dp1)*(1.-Atm%q_con(i,j,k))
        part2=rdg*Atm%delp(i,j,k)*(1.+dp1)*(1.-Atm%q_con(i,j,k))      &
              /delz_ptr(i,j,k)
        field(i,j,k)=exp((log(field(i,j,k))-log(part1)+cappa*log(part2)) &
                        /(1.-cappa))
      enddo
      enddo
      enddo
!
!---------------------------------------------------------------------
      end subroutine compute_halo_t
!---------------------------------------------------------------------
!
      subroutine fill_delz
!
!---------------------------------------------------------------------
!
      integer :: i,j,k
      integer :: lbnd1,lbnd2,ubnd1,ubnd2
!
!---------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
!***  Now fill the boundary rows using data from the BC files.
!---------------------------------------------------------------------
!
      do k=1,nz
      do j=j1,j2
      do i=i1,i2
        field(i,j,k)=delz_ptr(i,j,k)
      enddo
      enddo
      enddo
!
!---------------------------------------------------------------------
      end subroutine fill_delz
!---------------------------------------------------------------------
!
      end subroutine apply_delz_boundary
!
!---------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------------------------------------------------------------------

  subroutine exch_uv(domain, bd, npz, u, v)
    use mpi

    implicit none

    type(domain2d), intent(inout) :: domain
    type(fv_grid_bounds_type), intent(in) :: bd
    integer, intent(in) :: npz
    real, intent(inout) :: u   (bd%isd:bd%ied  ,bd%jsd:bd%jed+1,1:npz)
    real, intent(inout) :: v   (bd%isd:bd%ied+1,bd%jsd:bd%jed  ,1:npz)

    real, dimension(:), allocatable :: buf1,buf2,buf3,buf4
    integer :: ihandle1,ihandle2,ihandle3,ihandle4
    integer,dimension(MPI_STATUS_SIZE) :: istat
    integer :: ic, i, j, k, is, ie, js, je
    integer :: irecv, isend, ierr

    integer :: mype
    integer :: north_pe, south_pe, east_pe, west_pe

    mype = mpp_pe()
    call mpp_get_neighbor_pe( domain, NORTH, north_pe)
    call mpp_get_neighbor_pe( domain, SOUTH, south_pe)
    call mpp_get_neighbor_pe( domain, WEST,  west_pe)
    call mpp_get_neighbor_pe( domain, EAST,  east_pe)

    ! write(0,*) ' north_pe = ', north_pe
    ! write(0,*) ' south_pe = ', south_pe
    ! write(0,*) ' west_pe  = ', west_pe
    ! write(0,*) ' east_pe  = ', east_pe

    is=bd%is
    ie=bd%ie
    js=bd%js
    je=bd%je

    ! The size of these buffers must match the number of indices
    ! required below to send/receive the data. In particular,
    ! buf1 and buf4 must be of the same size (sim. for buf2 and buf3).
    ! Changes to the code below should be tested with debug flags
    ! enabled (out-of-bounds reads/writes).
    allocate(buf1(1:24*npz))
    allocate(buf2(1:36*npz))
    allocate(buf3(1:36*npz))
    allocate(buf4(1:24*npz))

! FIXME: MPI_COMM_WORLD

#ifdef OVERLOAD_R4
#define _DYN_MPI_REAL MPI_REAL
#else
#define _DYN_MPI_REAL MPI_DOUBLE_PRECISION
#endif

! Receive from north
    if( north_pe /= NULL_PE )then
       call MPI_Irecv(buf1,size(buf1),_DYN_MPI_REAL,north_pe,north_pe &
                     ,MPI_COMM_WORLD,ihandle1,irecv)
    endif

! Receive from south
    if( south_pe /= NULL_PE )then
       call MPI_Irecv(buf2,size(buf2),_DYN_MPI_REAL,south_pe,south_pe &
                     ,MPI_COMM_WORLD,ihandle2,irecv)
    endif

! Send to north
    if( north_pe /= NULL_PE )then
       ic=0
       do k=1,npz

         do j=je-3+1,je-1+1
         do i=is-3,is-1
           ic=ic+1
           buf3(ic)=u(i,j,k)
         enddo
         do i=ie+1,ie+3
           ic=ic+1
           buf3(ic)=u(i,j,k)
         enddo
         enddo

         do j=je-2,je
         do i=is-3,is-1
           ic=ic+1
           buf3(ic)=v(i,j,k)
         enddo
         do i=ie+1,ie+3
           ic=ic+1
           buf3(ic)=v(i,j,k)
         enddo
         enddo
       enddo
       if (ic/=size(buf2).or.ic/=size(buf3)) &
         call mpp_error(FATAL,'Buffer sizes buf2 and buf3 in routine exch_uv do not match actual message size')
       call MPI_Issend(buf3,size(buf3),_DYN_MPI_REAL,north_pe,mype &
                      ,MPI_COMM_WORLD,ihandle3,isend)
    endif

! Send to south
    if( south_pe /= NULL_PE )then
       ic=0
       do k=1,npz

         do j=js+2,js+3
         do i=is-3,is-1
           ic=ic+1
           buf4(ic)=u(i,j,k)
         enddo
         do i=ie+1,ie+3
           ic=ic+1
           buf4(ic)=u(i,j,k)
         enddo
         enddo

         do j=js+1,js+2
         do i=is-3,is-1
           ic=ic+1
           buf4(ic)=v(i,j,k)
         enddo
         do i=ie+1,ie+3
           ic=ic+1
           buf4(ic)=v(i,j,k)
         enddo
         enddo

       enddo
       if (ic/=size(buf1).or.ic/=size(buf4)) &
         call mpp_error(FATAL,'Buffer sizes buf1 and buf4 in routine exch_uv do not match actual message size')
       call MPI_Issend(buf4,size(buf4),_DYN_MPI_REAL,south_pe,mype &
                      ,MPI_COMM_WORLD,ihandle4,isend)
    endif

! Store from south
    if( south_pe /= NULL_PE )then
       ic=0
       call MPI_Wait(ihandle2,istat,ierr)
       do k=1,npz

         do j=js-3,js-1
         do i=is-3,is-1
           ic=ic+1
           u(i,j,k)=buf2(ic)
         enddo
         do i=ie+1,ie+3
           ic=ic+1
           u(i,j,k)=buf2(ic)
         enddo
         enddo

         do j=js-3,js-1
         do i=is-3,is-1
           ic=ic+1
           v(i,j,k)=buf2(ic)
         enddo
         do i=ie+1,ie+3
           ic=ic+1
           v(i,j,k)=buf2(ic)
         enddo
         enddo

       enddo
    endif

! Store from north
    if( north_pe /= NULL_PE )then
       ic=0
       call MPI_Wait(ihandle1,istat,ierr)
       do k=1,npz

         do j=je+2+1,je+3+1
         do i=is-3,is-1
           ic=ic+1
           u(i,j,k)=buf1(ic)
         enddo
         do i=ie+1,ie+3
           ic=ic+1
           u(i,j,k)=buf1(ic)
         enddo
         enddo

         do j=je+2,je+3
         do i=is-3,is-1
           ic=ic+1
           v(i,j,k)=buf1(ic)
         enddo
         do i=ie+1,ie+3
           ic=ic+1
           v(i,j,k)=buf1(ic)
         enddo
         enddo

       enddo
    endif

    deallocate(buf1)
    deallocate(buf2)
    deallocate(buf3)
    deallocate(buf4)

  end subroutine exch_uv

!---------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------------------------------------------------------------------

  subroutine get_data_source(source,regional)
!
! This routine extracts the data source information if it is present in the datafile.
!
      character (len = 80) :: source
      integer              :: ncids,sourceLength
      logical :: lstatus,regional
!
! Use the fms call here so we can actually get the return code value.
!
      if (regional) then
       lstatus = get_global_att_value('INPUT/gfs_data.nc',"source", source)
      else
       lstatus = get_global_att_value('INPUT/gfs_data.tile1.nc',"source", source)
      endif
      if (.not. lstatus) then
       if (mpp_pe() == 0) write(0,*) 'INPUT source not found ',lstatus,' set source=No Source Attribute'
       source='No Source Attribute'
      endif
  end subroutine get_data_source

!---------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------------------------------------------------------------------

  subroutine set_delp_and_tracers(BC_side,npz,nwat)
!
! This routine mimics what is done in external_ic to add mass back to delp
! and remove it from the tracers
!
  integer :: npz,nwat
  type(fv_regional_BC_variables),intent(inout) :: BC_side   !<-- The BC variables on a domain side at the final integration levels.
!
! local variables
!
  integer :: k, j, i, iq, is, ie, js, je
  integer :: liq_wat, ice_wat, rainwat, snowwat, graupel, cld_amt
  real    :: qt, wt, m_fac

   is=lbound(BC_side%delp_BC,1)
   ie=ubound(BC_side%delp_BC,1)
   js=lbound(BC_side%delp_BC,2)
   je=ubound(BC_side%delp_BC,2)
!
   liq_wat = get_tracer_index(MODEL_ATMOS, 'liq_wat')
   ice_wat = get_tracer_index(MODEL_ATMOS, 'ice_wat')
   rainwat = get_tracer_index(MODEL_ATMOS, 'rainwat')
   snowwat = get_tracer_index(MODEL_ATMOS, 'snowwat')
   graupel = get_tracer_index(MODEL_ATMOS, 'graupel')
   cld_amt = get_tracer_index(MODEL_ATMOS, 'cld_amt')
!
   source: if (trim(data_source) == 'FV3GFS GAUSSIAN NEMSIO FILE') then
!
!    if (cld_amt > 0) BC_side%q_BC(:,:,:,cld_amt) = 0.0    ! Moorthi
     do k=1,npz
     do j=js,je
     do i=is,ie
       wt = BC_side%delp_BC(i,j,k)
       if ( nwat == 6 ) then
         qt = wt*(1. + BC_side%q_BC(i,j,k,liq_wat) + &
                       BC_side%q_BC(i,j,k,ice_wat) + &
                       BC_side%q_BC(i,j,k,rainwat) + &
                       BC_side%q_BC(i,j,k,snowwat) + &
                       BC_side%q_BC(i,j,k,graupel))
       else   ! all other values of nwat
         qt = wt*(1. + sum(BC_side%q_BC(i,j,k,2:nwat)))
       endif
!--- Adjust delp with tracer mass.
       BC_side%delp_BC(i,j,k) = qt
     enddo
     enddo
     enddo
!
   else source   ! This else block is for all sources other than FV3GFS GAUSSIAN NEMSIO FILE
!
! 20160928: Adjust the mixing ratios consistently...
     do k=1,npz
     do j=js,je
     do i=is,ie
       wt = BC_side%delp_BC(i,j,k)
       if ( nwat == 6 ) then
         qt = wt*(1. + BC_side%q_BC(i,j,k,liq_wat) + &
                       BC_side%q_BC(i,j,k,ice_wat) + &
                       BC_side%q_BC(i,j,k,rainwat) + &
                       BC_side%q_BC(i,j,k,snowwat) + &
                       BC_side%q_BC(i,j,k,graupel))
       else   ! all other values of nwat
         qt = wt*(1. + sum(BC_side%q_BC(i,j,k,2:nwat)))
       endif
       m_fac = wt / qt
       do iq=1,ntracers
         BC_side%q_BC(i,j,k,iq) = m_fac * BC_side%q_BC(i,j,k,iq)
       enddo
       BC_side%delp_BC(i,j,k) = qt
     enddo
     enddo
     enddo
!
   endif source
!
  end subroutine set_delp_and_tracers

!---------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------------------------------------------------------------------

end module fv_regional_mod

!---------------------------------------------------------------------
