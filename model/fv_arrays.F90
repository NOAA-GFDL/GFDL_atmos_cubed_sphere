module fv_arrays_mod
#include <fms_platform.h>
 use mpp_domains_mod,  only: domain2d
  use fms_io_mod,       only: restart_file_type
  use time_manager_mod, only: time_type
  use horiz_interp_type_mod, only:  horiz_interp_type
  use mpp_domains_mod, only : nest_domain_type
  public

  !Several 'auxiliary' structures are introduced here. These are for
  ! the internal use by certain modules, and although fv_atmos_type
  !  contains one of each of these structures all memory management
  !   is performed by the module in question.

  integer, parameter:: max_step = 1000
  type fv_diag_type


 integer ::id_ps, id_slp, id_ua, id_va, id_pt, id_omga, id_vort,  &
           id_tm, id_pv, id_zsurf, id_oro, id_sgh, id_divg, id_w, &
           id_te, id_zs, id_ze, id_mq, id_vorts, id_us, id_vs,    &
           id_tq, id_rh, id_c15, id_c25, id_c35, id_c45,          &
                         id_f15, id_f25, id_f35, id_f45,          &
           id_ppt, id_ts, id_pmask, id_pmaskv2,                   &
           id_delp, id_delz, id_zratio, id_ws

! Selected p-level fields from 3D variables:
 integer :: id_vort850, id_w850,  &
            id_w200, id_s200, id_sl12, id_sl13
! IPCC diag
 integer :: id_u10,  id_v10,  id_t10,  id_q10,  id_rh10,  id_omg10,  id_h10,  &
            id_u50,  id_v50,  id_t50,  id_q50,  id_rh50,  id_omg50,  id_h50,  &
            id_u100, id_v100, id_t100, id_q100, id_rh100, id_omg100, id_h100, &
            id_u200, id_v200, id_t200, id_q200, id_rh200, id_omg200, id_h200, &
            id_u250, id_v250, id_t250, id_q250, id_rh250, id_omg250, id_h250, &
            id_u300, id_v300, id_t300, id_q300, id_rh300, id_omg300, id_h300, &
            id_u500, id_v500, id_t500, id_q500, id_rh500, id_omg500, id_h500, &
            id_u700, id_v700, id_t700, id_q700, id_rh700, id_omg700, id_h700, &
            id_u850, id_v850, id_t850, id_q850, id_rh850, id_omg850, id_h850, &
            id_u1000,id_v1000,id_t1000,id_q1000,id_rh1000,id_omg1000,id_h1000
 integer :: id_rh1000_cmip, id_rh850_cmip, id_rh700_cmip, id_rh500_cmip, &
            id_rh300_cmip, id_rh250_cmip, id_rh100_cmip, id_rh50_cmip, id_rh10_cmip
 integer :: id_hght

#ifdef MARS_GCM
 integer ::  id_t05
 integer ::  id_tdust, id_sfc_dust
#endif MARS_GCM

     ! For initial conditions:
     integer ic_ps, ic_ua, ic_va, ic_ppt
#ifdef LASPRAT
     integer ic_sphum
#endif
     integer, allocatable :: id_tracer(:)
! ESM requested diagnostics  -  dry mass/volume mixing ratios
 integer, allocatable :: id_tracer_dmmr(:)
 integer, allocatable :: id_tracer_dvmr(:)
 real,    allocatable :: w_mr(:)

     real, allocatable :: phalf(:)
     real, allocatable :: zsurf(:,:)
     real, allocatable :: zxg(:,:)
     real, allocatable :: pt1(:)


     logical :: initialized = .false.
     real  sphum, liq_wat, ice_wat       ! GFDL physics
     real  rainwat, snowwat, graupel

     real :: efx(max_step), efx_sum, efx_nest(max_step), efx_sum_nest, mtq(max_step), mtq_sum
     integer :: steps

  end type fv_diag_type


  type fv_atmos_type

     logical :: allocated = .false.
     logical :: dummy = .false.
     integer :: grid_number = 1

     !Timestep-related variables.
     !Each grid should have its own set of timing utilities
     real :: dt_grid
     integer :: refinement ! number of this grid's timesteps per model-wide dt_atmos
     type(time_type) :: Time_init, Time, Run_length, Time_end, Time_step_atmos

     logical :: grid_active = .true. !Always active for now

!-----------------------------------------------------------------------
! Five prognostic state variables for the f-v dynamics
!-----------------------------------------------------------------------
! dyn_state:
! D-grid prognostatic variables: u, v, and delp (and other scalars)
!
!     o--------u(i,j+1)----------o
!     |           |              |
!     |           |              |
!  v(i,j)------scalar(i,j)----v(i+1,j)
!     |           |              |
!     |           |              |
!     o--------u(i,j)------------o
!
! The C grid component is "diagnostic" in that it is predicted every time step
! from the D grid variables.
    real, _ALLOCATABLE :: u(:,:,:)    _NULL  ! D grid zonal wind (m/s)
    real, _ALLOCATABLE :: v(:,:,:)    _NULL  ! D grid meridional wind (m/s)
    real, _ALLOCATABLE :: pt(:,:,:)   _NULL  ! temperature (K)
    real, _ALLOCATABLE :: delp(:,:,:) _NULL  ! pressure thickness (pascal)
    real, _ALLOCATABLE :: q(:,:,:,:)  _NULL  ! specific humidity and constituents

!----------------------
! non-hydrostatic state:
!----------------------------------------------------------------------
    real, _ALLOCATABLE ::     w(:,:,:)  _NULL  ! cell center vertical wind (m/s)
    real, _ALLOCATABLE ::  delz(:,:,:)  _NULL  ! layer thickness (meters)
    real, _ALLOCATABLE ::   ze0(:,:,:)  _NULL  ! height at layer edges for remapping

!-----------------------------------------------------------------------
! Auxilliary pressure arrays:
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
    real, _ALLOCATABLE :: ps (:,:)      _NULL  ! Surface pressure (pascal)
    real, _ALLOCATABLE :: pe (:,:,: )   _NULL  ! edge pressure (pascal)
    real, _ALLOCATABLE :: pk  (:,:,:)   _NULL  ! pe**cappa
    real, _ALLOCATABLE :: peln(:,:,:)   _NULL  ! ln(pe)
    real, _ALLOCATABLE :: pkz (:,:,:)   _NULL  ! finite-volume mean pk

! For phys coupling:
    real, _ALLOCATABLE :: u_srf(:,:)    _NULL  ! Surface u-wind
    real, _ALLOCATABLE :: v_srf(:,:)    _NULL  ! Surface v-wind
    real, _ALLOCATABLE :: sgh(:,:)      _NULL  ! Terrain standard deviation
    real, _ALLOCATABLE :: oro(:,:)      _NULL  ! land fraction (1: all land; 0: all water)
    real, _ALLOCATABLE :: ts(:,:)       _NULL  ! skin temperature (sst) from NCEP/GFS (K) -- tile
 
!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
    real, _ALLOCATABLE :: phis(:,:)     _NULL  ! Surface geopotential (g*Z_surf)
    real, _ALLOCATABLE :: omga(:,:,:)   _NULL  ! Vertical pressure velocity (pa/s)
    real, _ALLOCATABLE :: ua(:,:,:)     _NULL  ! (ua, va) are mostly used as the A grid winds
    real, _ALLOCATABLE :: va(:,:,:)     _NULL
    real, _ALLOCATABLE :: uc(:,:,:)     _NULL  ! (uc, vc) are mostly used as the C grid winds
    real, _ALLOCATABLE :: vc(:,:,:)     _NULL

    real, _ALLOCATABLE :: ak(:)  _NULL
    real, _ALLOCATABLE :: bk(:)  _NULL

! Accumulated Mass flux arrays
    real, _ALLOCATABLE ::  mfx(:,:,:)  _NULL
    real, _ALLOCATABLE ::  mfy(:,:,:)  _NULL
! Accumulated Courant number arrays
    real, _ALLOCATABLE ::  cx(:,:,:)  _NULL
    real, _ALLOCATABLE ::  cy(:,:,:)  _NULL

    integer :: isd, ied, jsd, jed
    integer :: is, ie, js, je
    integer :: isc, iec, jsc, jec 

!!!!! fv_control

!-----------------------------------------------------------------------
! Grid descriptor file setup
!-----------------------------------------------------------------------
   character(len=80) :: grid_name = 'Gnomonic'
   character(len=120):: grid_file = 'Inline'
  integer      :: grid_type = 0     ! -1: read from file; 0: ED Gnomonic
!                                    !  0: the "true" equal-distance Gnomonic grid
!                                    !  1: the traditional equal-distance Gnomonic grid
!                                    !  2: the equal-angular Gnomonic grid
!                                    !  3: the lat-lon grid -- to be implemented
!                                    !  4: double periodic boundary condition on Cartesian grid
!                                    !  5: channel flow on Cartesian grid
!  -> moved to grid_tools

! Momentum (or KE) options:
   integer :: hord_mt = 9    ! the best option for Gnomonic grids  
   integer :: kord_mt = 8    ! vertical mapping option for (u,v)
   integer :: kord_wz = 8    ! vertical mapping option for w

! Vorticity & w transport options:
   integer :: hord_vt = 9    ! 10 not recommended (noisy case-5) 

! Heat & air mass (delp) transport options:
   integer :: hord_tm = 9    ! virtual potential temperature
   integer :: hord_dp = 9    ! delp (positive definite)
   integer :: kord_tm =-8    !

! Tracer transport options:
   integer :: hord_tr = 12   !11: PPM mono constraint (Lin 2004); fast 
                             !12: Huynh 2nd constraint (Lin 2004) +
                             !    positive definite (Lin & Rood 1996); slower
                             !>12: positive definite only (Lin & Rood 1996); fastest
   integer :: kord_tr = 8    ! 
   real    :: scale_z = 0.  ! diff_z = scale_z**2 * 0.25
   real    :: w_max = 75.    ! max w (m/s) threshold for hydostatiic adjustment 
   real    :: z_min = 0.05   ! min ration of dz_nonhydrostatic/dz_hydrostatic

   integer :: nord=1         ! 0: del-2, 1: del-4, 2: del-6, 3: del-8 divergence damping
                             ! Alternative setting for high-res: nord=1; d4_bg = 0.075
   real    :: dddmp = 0.0    ! coefficient for del-2 divergence damping (0.2)
                             ! for C90 or lower: 0.2
   real    :: d2_bg = 0.0    ! coefficient for background del-2 divergence damping
   real    :: d4_bg = 0.16   ! coefficient for background del-4(6) divergence damping
                             ! for stability, d4_bg must be <=0.16 if nord=3
   real    :: vtdm4 = 0.0    ! coefficient for del-4 vorticity damping
   real    :: d2_bg_k1 = 4.         ! factor for d2_bg (k=1)
   real    :: d2_bg_k2 = 2.         ! factor for d2_bg (k=2)
   real    :: d2_divg_max_k1 = 0.15 ! d2_divg max value (k=1)
   real    :: d2_divg_max_k2 = 0.08 ! d2_divg max value (k=2)
   real    :: damp_k_k1 = 0.2       ! damp_k value (k=1)
   real    :: damp_k_k2 = 0.12      ! damp_k value (k=2)

! Additional (after the fact) terrain filter (to further smooth the terrain after cold start)
   integer ::    n_zs_filter=0      !  number of application of the terrain filter
   integer :: nord_zs_filter=4      !  use del-2 (2) OR del-4 (4)

   logical :: no_dycore = .false.   ! skip the dycore
   logical :: replace_w = .false.   ! replace w (m/sec) with omega (pa/sec) 
                                    ! this is useful for getting a good initial estimate of w
                                    ! suugest usage: from a hydrostatic IC, set hydrostatic = .F., make_nh = .T.
                                    !                replace_w = .T., then run the model adiabatically for 1 sec
   logical :: convert_ke = .false. 
   logical :: do_vort_damp = .false. 
   logical :: use_old_omega = .true. 
! PG off centering:
   real    :: beta  = 0.0    ! 0.5 is "neutral" but it may not be stable
#ifdef SW_DYNAMICS
   integer :: n_sponge = 0   ! Number of sponge layers at the top of the atmosphere
   real    :: d_ext = 0.    
   integer :: nwat  = 0      ! Number of water species
   logical :: warm_start = .false. 
   logical :: inline_q = .true.
#else
   integer :: n_sponge = 1   ! Number of sponge layers at the top of the atmosphere
   real    :: d_ext = 0.02   ! External model damping (was 0.02)
   integer :: nwat  = 3      ! Number of water species
   logical :: warm_start = .true. 
                             ! Set to .F. if cold_start is desired (including terrain generation)
   logical :: inline_q = .false.
#endif
!-----------------------------------------------------------
! Grid shifting, rotation, and the Schmidt transformation:
!-----------------------------------------------------------
   real :: shift_fac   =  18.   ! shift west by 180/shift_fac = 10 degrees
! Defaults for Schmidt transformation:
   logical :: do_schmidt = .false. 
   real :: stretch_fac =   1.   ! No stretching
   real :: target_lat  = -90.   ! -90: no grid rotation 
   real :: target_lon  =   0.   ! 

!-----------------------------------------------------------------------------------------------
! Example #1a: US regional climate simulation, center located over Oklahoma city: (262.4, 35.4)
!              stretching factor: 2.5
! Example #1b: US Hurricane model, center over Miami: (279.7, 25.8)
!              stretching factor: 3-5
! Example #2a: South-East Asia Regional Climate H*** (SERACH), Central Taiwan: (121.0, 23.5)
! Example #2b: Typhoon Model: (124.0, 22.0)
!              stretching factor: 5-10
!-----------------------------------------------------------------------------------------------

   logical :: reset_eta = .false. 
   real    :: p_fac  = 0.05
   real    :: a_imp  = 0.75  ! Off center parameter for the implicit solver [0.5,1.0]
   integer :: n_split = 0    ! Number of time splits for the lagrangian dynamics
                             ! Default = 0 (automatic computation of best value)
   integer :: m_split = 0    ! Number of time splits for Riemann solver
   integer :: k_split = 1    ! Number of time splits for Remapping

!            For doubly periodic domain with sim_phys
!                     5km        150         20 (7.5 s)  2
!
!                     Estimates for Gnomonic grids:
            !===================================================
            !        dx (km)    dt (sc)    n_split    m_split
            !===================================================
            ! C1000:  ~10        150         16          3
            ! C2000:   ~5         90         18 (5 s)    2
            !===================================================
! The nonhydrostatic algorithm is described in Lin 2006, QJ, (submitted)
! C2000 should easily scale to at least 6 * 100 * 100 = 60,000 CPUs  
! For a 1024 system: try 6 x 13 * 13 = 1014 CPUs
  
   integer :: q_split = 0    ! Number of time splits for tracer transport
   integer :: print_freq = 0 ! Print max/min of selected fields
                             ! 0: off
                             ! positive n: every n hours
                             ! negative n: every time step

!------------------------------------------
! Model Domain parameters
!------------------------------------------
   integer :: npx                     ! Number of Grid Points in X- dir
   integer :: npy                     ! Number of Grid Points in Y- dir
   integer :: npz                     ! Number of Vertical Levels
   integer :: npz_rst = 0             ! Original Vertical Levels (in the restart)
                                      ! 0: no change (default)
   integer :: ncnst = 0               ! Number of advected consituents
   integer :: pnats = 0               ! Number of non-advected consituents
   integer :: ntiles = 1                 ! Number or tiles that make up the Grid 
   integer :: ndims = 2     ! Lat-Lon Dims for Grid in Radians
   integer :: nf_omega  = 1           ! Filter omega "nf_omega" times
   integer :: fv_sg_adj = -1          ! Perform grid-scale dry adjustment if > 0
                                      ! Relaxzation time  scale (sec) if positive
   integer :: na_init = 0             ! Perform adiabatic initialization
#ifdef MARS_GCM
   real    :: p_ref = 600.
   real    :: reference_sfc_pres = 7.7E2
   real    :: sponge_damp=   1.0
   real    :: dry_mass = 7.7E2
#else
   real    :: p_ref = 1.E5
   real    :: dry_mass = 98290.
#endif
   integer :: nt_prog = 0
   integer :: nt_phys = 0
   real    :: tau_h2o = 0.            ! Time scale (days) for ch4_chem

   real    :: d_con = 0.
   real    :: consv_te = 0.
   real    :: tau = 0.                ! Time scale (days) for Rayleigh friction
   real    :: rf_center = 0.          ! Center position of the hyper-tan profile
                                      ! 0: use the top layer center
                                      ! > 0, [Pascal]
   logical :: tq_filter = .false.
   logical :: filter_phys = .false.
   logical :: dwind_2d = .false.
   logical :: breed_vortex_inline = .false.
   logical :: range_warn = .false.
   logical :: fill = .false.
   logical :: fill_dp = .false.
   logical :: fill_wz = .false.
   logical :: non_ortho = .true.
   logical :: adiabatic = .false.     ! Run without physics (full or idealized).
   logical :: moist_phys = .true.     ! Run with moist physics
   logical :: do_Held_Suarez = .false.
   logical :: reproduce_sum = .true.  ! Make global sum for consv_te reproduce
   logical :: adjust_dry_mass = .false.
   logical :: fv_debug  = .false.
   logical :: srf_init  = .false.
   logical :: mountain  = .true.
   logical :: remap_t  = .true.
   logical :: z_tracer = .false.      ! transport tracers layer by layer with independent
                                      ! time split; use this if tracer number is huge and/or
                                      ! high resolution (nsplt > 1)

   logical :: old_divg_damp = .false. ! parameter to revert damping parameters back to values
                                      ! defined in a previous revision
                                      ! old_values:
                                      !    d2_bg_k1 = 6.           d2_bg_k2 = 4.
                                      !    d2_divg_max_k1 = 0.02   d2_divg_max_k2 = 0.01
                                      !    damp_k_k1 = 0.          damp_k_k2 = 0.
                                      ! current_values:
                                      !    d2_bg_k1 = 4.           d2_bg_k2 = 2.
                                      !    d2_divg_max_k1 = 0.15   d2_divg_max_k2 = 0.08
                                      !    damp_k_k1 = 0.2         damp_k_k2 = 0.12
   logical :: master

   logical :: fv_land = .false.       ! To cold starting the model with USGS terrain
!--------------------------------------------------------------------------------------
! The following options are useful for NWP experiments using datasets on the lat-lon grid
!--------------------------------------------------------------------------------------
   logical :: nudge = .false.         ! Perform nudging
   logical :: ncep_ic = .false.       ! use NCEP ICs 
   logical :: fv_diag_ic = .false.    ! reconstruct IC from fv_diagnostics on lat-lon grid
   logical :: external_ic = .false.   ! use ICs from external sources; e.g. lat-lon FV core
                                      ! or NCEP re-analysis; both vertical remapping & horizontal
                                      ! (lat-lon to cubed sphere) interpolation will be done
! Default restart files from the "Memphis" latlon FV core:
   character(len=128) :: res_latlon_dynamics = 'INPUT/fv_rst.res.nc'
   character(len=128) :: res_latlon_tracers  = 'INPUT/atmos_tracers.res.nc'
! The user also needs to copy the "cold start" cubed sphere restart files (fv_core.res.tile1-6)
! to the INPUT dir during runtime
!------------------------------------------------
! Parameters related to non-hydrostatic dynamics:
!------------------------------------------------
   logical :: hydrostatic = .true.
   logical :: phys_hydrostatic = .true.    ! heating/cooling term from the physics is hydrostatic
   logical :: hybrid_z    = .false. ! use hybrid_z for remapping
   logical :: Make_NH     = .false. ! Initialize (w, delz) from hydro restart file 
   logical :: make_hybrid_z  = .false. ! transform hydrostatic eta-coord IC into non-hydrostatic hybrid_z

   integer :: a2b_ord = 4    ! order for interpolation from A to B Grid (corners)
   integer :: c2l_ord = 4    ! order for interpolation from D to lat-lon A winds for phys & output

   integer :: ks

!!!!!!!!!!!!!!!!!!
     ! From fv_mp_mod !
!!!!!!!!!!!!!!!!!!

     integer :: ng = 3 !this SHOULD be a constant, but structure elements are not allowed to be constants
     type(domain2D) :: domain
#if defined(SPMD)

     type(domain2D) :: domain_for_coupler ! domain used in coupled model with halo = 1.

     integer :: num_contact, npes_per_tile, tile, npes_this_grid
     logical :: square_domain = .false.
     integer :: npes_x, npes_y
     integer :: layout(2), io_layout(2) = (/ 1,1 /)

#endif
     !These do not actually belong to the grid, but to the process
     !integer :: masterproc
     !integer :: gid 

     type(nest_domain_type) :: nest_domain !Structure holding link from this grid to its parent
     logical :: this_proc_sends, this_proc_recvs

!!!!!!!!!!!!!!!!
! From fv_grid_tools
!!!!!!!!!!!!!!!!

     
  real :: csFac = -999
  real            :: zeta = 1.0                ! non-linear flag 
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
  integer :: interpOrder = 1
  logical :: debug_message_size = .false.
  logical :: write_grid_char_file = .false.
  logical :: stretched_grid = .false.
  ! grid descriptors

  ! Horizontal
  integer :: npx_g, npy_g, npz_g, ntiles_g ! global domain
#ifndef NO_GRID_G
  real, allocatable, dimension(:,:,:) :: grid_g
#endif
  real, allocatable, dimension(:,:,:) :: grid, agrid
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
  
  integer, dimension(:,:,:), allocatable :: iinta, jinta, iintb, jintb

  real :: dx_const = 1000.    ! spatial resolution for double periodic boundary configuration [m]
  real :: dy_const = 1000.
  real :: deglon_start = -30., deglon_stop = 30., &  ! boundaries of latlon patch
          deglat_start = -30., deglat_stop = 30.

!!!!!!!!!!!!!!!!
!fv_diagnostics!
!!!!!!!!!!!!!!!!

     type(fv_diag_type) :: idiag
     integer steps
     real(kind=4):: efx(max_step), efx_nest(max_step), mtq(max_step)
     real(kind=4):: efx_sum,       efx_sum_nest, mtq_sum


!!!!!!!!!!!!!!!!!!!!!!
     ! From fv_grid_utils !
!!!!!!!!!!!!!!!!!!!!!!

     ! Scalars:
     real, allocatable :: edge_s(:)
     real, allocatable :: edge_n(:)
     real, allocatable :: edge_w(:)
     real, allocatable :: edge_e(:)
     ! Vector:
     real, allocatable :: edge_vect_s(:)
     real, allocatable :: edge_vect_n(:)
     real, allocatable :: edge_vect_w(:)
     real, allocatable :: edge_vect_e(:)
     ! scalar:
     real, allocatable :: ex_s(:)
     real, allocatable :: ex_n(:)
     real, allocatable :: ex_w(:)
     real, allocatable :: ex_e(:)

     ! divergence Damping:
     real, allocatable :: divg_u(:,:), divg_v(:,:)    !
     ! Cubed_2_latlon:
     real, allocatable :: a11(:,:)
     real, allocatable :: a12(:,:)
     real, allocatable :: a21(:,:)
     real, allocatable :: a22(:,:)
     ! latlon_2_cubed:
     real, allocatable :: z11(:,:)
     real, allocatable :: z12(:,:)
     real, allocatable :: z21(:,:)
     real, allocatable :: z22(:,:)

     real :: global_area
     real:: stretch_factor = 1.
     logical :: g_sum_initialized = .false. !Not currently used but can be useful
     logical:: gnomonic_grid
     logical:: sw_corner, se_corner, ne_corner, nw_corner
     real, allocatable :: cosa_u(:,:)
     real, allocatable :: cosa_v(:,:)
     real, allocatable :: cosa_s(:,:)
     real, allocatable :: sina_s(:,:)
     real, allocatable :: sina_u(:,:)
     real, allocatable :: sina_v(:,:)
     real, allocatable :: rsin_u(:,:)
     real, allocatable :: rsin_v(:,:)
     real, allocatable ::  rsina(:,:)
     real, allocatable ::  rsin2(:,:)
     real, allocatable :: ee1(:,:,:)
     real, allocatable :: ee2(:,:,:)
     real, allocatable :: ec1(:,:,:)
     real, allocatable :: ec2(:,:,:)
     real, allocatable :: ew(:,:,:,:)
     real, allocatable :: es(:,:,:,:)


     !- 3D Super grid to contain all geometrical factors --
     ! the 3rd dimension is 9
     real, allocatable :: sin_sg(:,:,:)
     real, allocatable :: cos_sg(:,:,:)
     !--------------------------------------------------

     ! Unit Normal vectors at cell edges:
     real, allocatable :: en1(:,:,:)
     real, allocatable :: en2(:,:,:)

     ! Extended Cubed cross-edge winds
     real, allocatable :: eww(:,:)
     real, allocatable :: ess(:,:)

     ! Unit vectors for lat-lon grid
     real, allocatable :: vlon(:,:,:), vlat(:,:,:)
     real, allocatable :: fC(:,:), f0(:,:)

     real :: deglat=15.

     real    :: ptop
     real :: da_min, da_max, da_min_c, da_max_c


!!!!!!!!!!!!!!
! From fv_sg !
!!!!!!!!!!!!!!

     real, allocatable :: fv_olr(:,:), fv_abs_sw(:,:)
     integer:: irad = 0

!!!!!!!!!!!!!!!!!!!!
! From fv_surf_map !
!!!!!!!!!!!!!!!!!!!!

     real, allocatable:: zs_g(:,:)
     real, allocatable:: sgh_g(:,:), oro_g(:,:)

!!!!!!!!!!!!!!!!!!!
! From test_cases !
!!!!!!!!!!!!!!!!!!!

     integer :: test_case = -999 !If not specified, this at least tells us
                                 !that the nested-grid halo topography needs
                                 !to be filled by interpolating from the coarse grid.
     ! alpha = angle of axis rotation about the poles
     real   :: alpha = 0.
     ! Ubar = initial wind speed parameter
     real   :: Ubar
     ! gh0 = initial surface height parameter
     real   :: gh0

     !  case 9 parameters
     real  , allocatable :: case9_B(:,:)
     real   :: AofT(2)


     !  Validating fields used in statistics
     real  , allocatable :: phi0(:,:,:) ! Validating Field
     real  , allocatable :: ua0(:,:,:)  ! Validating U-Wind
     real  , allocatable :: va0(:,:,:)  ! Validating V-Windfms_io_exit, get_tile_string, &

     real  , allocatable :: gh_table(:), lats_table(:)
     logical :: gh_initialized = .false.

     !  Initial Conservation statistics ; total mass ; enstrophy ; energy
     real   :: tmass_orig
     real   :: tvort_orig
     real   :: tener_orig



!!!!!!!!!!!!!!!!
! From fv_phys !
!!!!!!!!!!!!!!!!
     !This is exclusively for nudging

     logical :: nudge_initialized
     real, allocatable:: u0(:,:,:), v0(:,:,:), t0(:,:,:), dp(:,:,:)

!!!!!!!!!!!!!
! From hswf !
!!!!!!!!!!!!!

#ifdef MARS_GCM
     logical :: tmars_initialized = .false.
     real,  allocatable, dimension(:,:,:):: tmars
#endif

!!!!!!!!!!!!!!
! From fv_io !
!!!!!!!!!!!!!!
     type(restart_file_type) :: Fv_restart, SST_restart, Fv_tile_restart, &
          Rsf_restart, Mg_restart, Lnd_restart, Tra_restart

!!!!!!!!!!!!!!!!!!!!!!!
! Nesting information !
!!!!!!!!!!!!!!!!!!!!!!!

     type(fv_atmos_type), pointer :: parent_grid _NULL
     integer :: parent_tile = 1     !Tile (of cubed sphere) in which nested grid lies 
     integer :: nest_refinement = 3 !Refinement wrt parent
     logical :: nested = .false.
     integer :: nestbctype = 1
     integer :: nsponge = 0
     integer :: nestupdate = 0       
     logical :: twowaynest = .false. 
     integer :: ioffset, joffset !Position of nest within parent grid
     integer :: nest_timestep = 0 !Counter for nested-grid timesteps
     integer :: tracer_nest_timestep = 0 !Counter for nested-grid timesteps
     real    :: s_weight = 1.e-6 !sponge weight
     logical :: first_step = .true.
     integer :: refinement_of_global = 1
     integer :: npx_global

     !Interpolation arrays for grid nesting
     integer, allocatable, dimension(:,:,:) :: ind_h, ind_u, ind_v
     real, allocatable, dimension(:,:,:) :: wt_h, wt_u, wt_v
     integer, allocatable, dimension(:,:,:) :: ind_update_h

     !These arrays are not allocated by allocate_fv_atmos_type; but instead
     !allocated for all grids, regardless of whether the grid is
     !on a PE of a concurrent run.
     logical, allocatable, dimension(:) :: child_grids
     integer, allocatable, dimension(:) :: pelist

     logical :: parent_proc, child_proc

!!$     !Indices for coarse-to-nested interpolation
!!$   integer                      :: isw_f, iew_f, jsw_f, jew_f, isw_c, iew_c, jsw_c, jew_c
!!$   integer                      :: ise_f, iee_f, jse_f, jee_f, ise_c, iee_c, jse_c, jee_c
!!$   integer                      :: iss_f, ies_f, jss_f, jes_f, iss_c, ies_c, jss_c, jes_c
!!$   integer                      :: isn_f, ien_f, jsn_f, jen_f, isn_c, ien_c, jsn_c, jen_c

     !Arrays to hold data for time interpolation
     real, allocatable, dimension(:,:,:) :: h_west, h_east, h_north, h_south
     real, allocatable, dimension(:,:,:) :: u_west, u_east, u_north, u_south
     real, allocatable, dimension(:,:,:) :: v_west, v_east, v_north, v_south

     !These are for time-extrapolated BCs when doing concurrent nesting
     real, allocatable, dimension(:,:,:) :: h_west_t0, h_east_t0, h_north_t0, h_south_t0
     real, allocatable, dimension(:,:,:) :: u_west_t0, u_east_t0, u_north_t0, u_south_t0
     real, allocatable, dimension(:,:,:) :: v_west_t0, v_east_t0, v_north_t0, v_south_t0
     !uc and vc always need two time levels because they are interpolated to half-timesteps, not full ones
     real, allocatable, dimension(:,:,:) :: uc_west_t0, uc_east_t0, uc_north_t0, uc_south_t0
     real, allocatable, dimension(:,:,:) :: vc_west_t0, vc_east_t0, vc_north_t0, vc_south_t0
     real, allocatable, dimension(:,:,:) :: uc_west_t1, uc_east_t1, uc_north_t1, uc_south_t1
     real, allocatable, dimension(:,:,:) :: vc_west_t1, vc_east_t1, vc_north_t1, vc_south_t1
     real, allocatable, dimension(:,:,:,:) :: q_west, q_east, q_north, q_south
     real, allocatable, dimension(:,:,:,:) :: q_west_t0, q_east_t0, q_north_t0, q_south_t0
#ifndef SW_DYNAMICS
     real, allocatable, dimension(:,:,:) :: pt_west, pt_east, pt_north, pt_south
     real, allocatable, dimension(:,:,:) :: w_west, w_east, w_north, w_south
     real, allocatable, dimension(:,:,:) :: pt_west_t0, pt_east_t0, pt_north_t0, pt_south_t0
     real, allocatable, dimension(:,:,:) :: w_west_t0, w_east_t0, w_north_t0, w_south_t0
#endif

     !These are for tracer flux BCs
     real, allocatable, dimension(:,:,:) :: nest_fx_west_accum, nest_fx_east_accum, nest_fx_south_accum, nest_fx_north_accum
     logical :: do_flux_BCs, do_2way_flux_BCs !For a parent grid; determine whether there is a need to send BCs

     integer, allocatable, dimension(:,:) :: process_bounds

     !Hold on to coarse-grid global grid, so we don't have to waste processor time getting it again when starting to do grid nesting
     real, allocatable, dimension(:,:,:,:) :: grid_global

!!!!!!!!!!!!!!!!!!!
! From fv_physics !
!!!!!!!!!!!!!!!!!!!
     
   real, allocatable, dimension(:,:,:)   :: t_phys
   real, allocatable, dimension(:,:,:,:) :: q_phys
   real, allocatable, dimension(:,:,:)   :: u_dt, v_dt, t_dt
   real, allocatable, dimension(:,:,:,:) :: q_dt  ! potentially a huge array
   real, allocatable, dimension(:,:,:)   :: p_full, z_full, p_half, z_half
   integer :: nx_win, ny_win      ! iew-isw+1, jew-jsw+1 (window sizes)
   integer :: nx_dom, ny_dom      ! ie-is+1, je-js+1 (compute domain sizes)
   integer, allocatable, dimension(:)  :: physics_window_x, physics_window_y !Allocated in fv_physics
   integer :: ny_per_thread, num_phys_windows

  integer :: atmos_axes(4)


  end type fv_atmos_type

  type(fv_atmos_type), allocatable :: Atm(:)

contains

  subroutine allocate_fv_atmos_type(Atm, isd_in, ied_in, jsd_in, jed_in, is_in, ie_in, js_in, je_in, &
       npx_in, npy_in, npz_in, ndims_in, ncnst_in, ng_in, dummy, concurrent, alloc_2d)

    !WARNING: Before calling this routine, be sure to have set up the
    ! proper domain parameters from the namelists (as is done in
    ! fv_control.F90)

    implicit none
    type(fv_atmos_type), intent(INOUT) :: Atm
    integer, intent(IN) :: isd_in, ied_in, jsd_in, jed_in, is_in, ie_in, js_in, je_in
    integer, intent(IN) :: npx_in, npy_in, npz_in, ndims_in, ncnst_in, ng_in
    logical, intent(IN) :: dummy, concurrent, alloc_2d
    integer:: isd, ied, jsd, jed, is, ie, js, je
    integer:: npx, npy, npz, ndims, ncnst, ng

    !For 2D utility arrays
    integer:: isd_2d, ied_2d, jsd_2d, jed_2d, is_2d, ie_2d, js_2d, je_2d
    integer:: npx_2d, npy_2d, npz_2d, ndims_2d, ncnst_2d, ng_2d

    integer :: ns

    if (Atm%allocated) return

    if (dummy) then
       isd     =  0   
       ied=   -1   
       jsd=   0   
       jed=   -1   
       is=   0    
       ie=   -1    
       js=   0    
       je=   -1    
       npx=   1   
       npy=   1   
       npz=   1   
       ndims=   1 
       ncnst=   1 
       ng     =   1   
    else
       isd     =  isd_in   
       ied=   ied_in   
       jsd=   jsd_in   
       jed=   jed_in   
       is=   is_in    
       ie=   ie_in    
       js=   js_in    
       je=   je_in    
       npx=   npx_in   
       npy=   npy_in   
       npz=   npz_in   
       ndims=   ndims_in 
       ncnst=   ncnst_in 
       ng     =   ng_in    
    endif

    if ((.not. dummy) .or. alloc_2d) then
       isd_2d     =  isd_in   
       ied_2d=   ied_in   
       jsd_2d=   jsd_in   
       jed_2d=   jed_in   
       is_2d=   is_in    
       ie_2d=   ie_in    
       js_2d=   js_in    
       je_2d=   je_in    
       npx_2d=   npx_in   
       npy_2d=   npy_in   
       npz_2d=   npz_in   
       ndims_2d=   ndims_in 
       ncnst_2d=   ncnst_in 
       ng_2d     =   ng_in 
    else
       isd_2d     =  0   
       ied_2d=   -1   
       jsd_2d=   0   
       jed_2d=   -1   
       is_2d=   0    
       ie_2d=   -1    
       js_2d=   0    
       je_2d=   -1    
       npx_2d=   1   
       npy_2d=   1   
       npz_2d=   0 !for ak, bk   
       ndims_2d=   1 
       ncnst_2d=   1 
       ng_2d     =   1        
    endif

    allocate (    Atm%u(isd:ied  ,jsd:jed+1,npz) )
    allocate (    Atm%v(isd:ied+1,jsd:jed  ,npz) )

    allocate (   Atm%pt(isd:ied  ,jsd:jed  ,npz) )
    allocate ( Atm%delp(isd:ied  ,jsd:jed  ,npz) )
    allocate (    Atm%q(isd:ied  ,jsd:jed  ,npz, ncnst) )

    ! Allocate Auxilliary pressure arrays
    allocate (   Atm%ps(isd:ied  ,jsd:jed) )
    allocate (   Atm%pe(is-1:ie+1, npz+1,js-1:je+1) )
    allocate (   Atm%pk(is:ie    ,js:je  , npz+1) )
    allocate ( Atm%peln(is:ie,npz+1,js:je) )
    allocate (  Atm%pkz(is:ie,js:je,npz) )

    allocate ( Atm%u_srf(is:ie,js:je) )
    allocate ( Atm%v_srf(is:ie,js:je) )

    if ( Atm%fv_land ) then
       allocate ( Atm%sgh(is:ie,js:je) )
       allocate ( Atm%oro(is:ie,js:je) )
    else
       allocate ( Atm%oro(1,1) )
    endif

    ! Allocate others
    allocate ( Atm%ts(is:ie,js:je) )
    allocate ( Atm%phis(isd:ied  ,jsd:jed  ) )
    allocate ( Atm%omga(isd:ied  ,jsd:jed  ,npz) ); Atm%omga=0.
    allocate (   Atm%ua(isd:ied  ,jsd:jed  ,npz) )
    allocate (   Atm%va(isd:ied  ,jsd:jed  ,npz) )
    allocate (   Atm%uc(isd:ied+1,jsd:jed  ,npz) )
    allocate (   Atm%vc(isd:ied  ,jsd:jed+1,npz) )
    ! For tracer transport:
    allocate ( Atm%mfx(is:ie+1, js:je,  npz) )
    allocate ( Atm%mfy(is:ie  , js:je+1,npz) )
    allocate (  Atm%cx(is:ie+1, jsd:jed, npz) )
    allocate (  Atm%cy(isd:ied ,js:je+1, npz) )

    allocate (  Atm%ak(npz_2d+1) )
    allocate (  Atm%bk(npz_2d+1) )

    !--------------------------
    ! Non-hydrostatic dynamics:
    !--------------------------
    if ( Atm%hydrostatic ) then
       allocate (    Atm%w(1, 1  ,1) )
       allocate ( Atm%delz(1, 1  ,1) )
       allocate (  Atm%ze0(1, 1  ,1) )
    else
       allocate (    Atm%w(isd:ied, jsd:jed  ,npz  ) )
       allocate ( Atm%delz(is:ie, js:je  ,npz) )
       if( Atm%hybrid_z ) then
          allocate (  Atm%ze0(is:ie, js:je ,npz+1) )
       else
          allocate (  Atm%ze0(1, 1  ,1) )
       endif
       !         allocate ( mono(isd:ied, jsd:jed, npz))
    endif

    allocate ( Atm% area(isd_2d:ied_2d  ,jsd_2d:jed_2d  ) )   ! Cell Centered
    allocate ( Atm%rarea(isd_2d:ied_2d  ,jsd_2d:jed_2d  ) )   ! Cell Centered
    
    allocate ( Atm% area_c(isd_2d:ied_2d+1,jsd_2d:jed_2d+1) )  ! Cell Corners
    allocate ( Atm%rarea_c(isd_2d:ied_2d+1,jsd_2d:jed_2d+1) )  ! Cell Corners
    
    allocate ( Atm% dx(isd_2d:ied_2d  ,jsd_2d:jed_2d+1) )
    allocate ( Atm%rdx(isd_2d:ied_2d  ,jsd_2d:jed_2d+1) )
    allocate ( Atm% dy(isd_2d:ied_2d+1,jsd_2d:jed_2d  ) )
    allocate ( Atm%rdy(isd_2d:ied_2d+1,jsd_2d:jed_2d  ) )
    
    allocate ( Atm% dxc(isd_2d:ied_2d+1,jsd_2d:jed_2d  ) )
    allocate ( Atm%rdxc(isd_2d:ied_2d+1,jsd_2d:jed_2d  ) )
    allocate ( Atm% dyc(isd_2d:ied_2d  ,jsd_2d:jed_2d+1) )
    allocate ( Atm%rdyc(isd_2d:ied_2d  ,jsd_2d:jed_2d+1) )
    
    allocate ( Atm% dxa(isd_2d:ied_2d  ,jsd_2d:jed_2d  ) )
    allocate ( Atm%rdxa(isd_2d:ied_2d  ,jsd_2d:jed_2d  ) )
    allocate ( Atm% dya(isd_2d:ied_2d  ,jsd_2d:jed_2d  ) )
    allocate ( Atm%rdya(isd_2d:ied_2d  ,jsd_2d:jed_2d  ) )
    
    allocate ( Atm%grid (isd_2d:ied_2d+1,jsd_2d:jed_2d+1,1:ndims_2d) )
    allocate ( Atm%agrid(isd_2d:ied_2d  ,jsd_2d:jed_2d  ,1:ndims_2d) )
    allocate ( Atm%sina(isd_2d:ied_2d+1,jsd_2d:jed_2d+1) )   ! SIN(angle of intersection)
    allocate ( Atm%rsina(is_2d:ie_2d+1,js_2d:je_2d+1) )      ! Why is the size different?
    allocate ( Atm%cosa(isd_2d:ied_2d+1,jsd_2d:jed_2d+1) )   ! COS(angle of intersection)
    
    allocate ( Atm%  e1(3,isd_2d:ied_2d+1,jsd_2d:jed_2d+1) )
    allocate ( Atm%  e2(3,isd_2d:ied_2d+1,jsd_2d:jed_2d+1) )

    allocate (Atm%iinta(4, isd_2d:ied_2d ,jsd_2d:jed_2d), &
         Atm%jinta(4, isd_2d:ied_2d ,jsd_2d:jed_2d),  &
         Atm%iintb(4, is_2d:ie_2d+1 ,js_2d:je_2d+1), &
         Atm%jintb(4, is_2d:ie_2d+1 ,js_2d:je_2d+1) )

#ifndef NO_GRID_G
    allocate ( Atm%grid_g(1:npx_2d,1:npy_2d,1:ndims_2d) )
#endif

    allocate ( Atm%edge_s(npx_2d) )
    allocate ( Atm%edge_n(npx_2d) )
    allocate ( Atm%edge_w(npy_2d) )
    allocate ( Atm%edge_e(npy_2d) )

    allocate ( Atm%edge_vect_s(isd_2d:ied_2d) )
    allocate ( Atm%edge_vect_n(isd_2d:ied_2d) )
    allocate ( Atm%edge_vect_w(jsd_2d:jed_2d) )
    allocate ( Atm%edge_vect_e(jsd_2d:jed_2d) )

    allocate ( Atm%ex_s(npx_2d) )
    allocate ( Atm%ex_n(npx_2d) )
    allocate ( Atm%ex_w(npy_2d) )
    allocate ( Atm%ex_e(npy_2d) )


    ! For diveregnce damping:
    allocate (  Atm%divg_u(isd_2d:ied_2d,  jsd_2d:jed_2d+1) )
    allocate (  Atm%divg_v(isd_2d:ied_2d+1,jsd_2d:jed_2d) )

    allocate (  Atm%z11(is_2d-1:ie_2d+1,js_2d-1:je_2d+1) )
    allocate (  Atm%z12(is_2d-1:ie_2d+1,js_2d-1:je_2d+1) )
    allocate (  Atm%z21(is_2d-1:ie_2d+1,js_2d-1:je_2d+1) )
    allocate (  Atm%z22(is_2d-1:ie_2d+1,js_2d-1:je_2d+1) )

    allocate (  Atm%a11(is_2d-1:ie_2d+1,js_2d-1:je_2d+1) )
    allocate (  Atm%a12(is_2d-1:ie_2d+1,js_2d-1:je_2d+1) )
    allocate (  Atm%a21(is_2d-1:ie_2d+1,js_2d-1:je_2d+1) )
    allocate (  Atm%a22(is_2d-1:ie_2d+1,js_2d-1:je_2d+1) )
    allocate ( Atm%vlon(is_2d-2:ie_2d+2,js_2d-2:je_2d+2,3) )
    allocate ( Atm%vlat(is_2d-2:ie_2d+2,js_2d-2:je_2d+2,3) )
    ! Coriolis parameters:
    allocate ( Atm%f0(isd_2d:ied_2d  ,jsd_2d:jed_2d  ) )
    allocate ( Atm%fC(isd_2d:ied_2d+1,jsd_2d:jed_2d+1) )

    ! Corner unit vectors:
    allocate( Atm%ee1(3,isd_2d:ied_2d+1,jsd_2d:jed_2d+1) )
    allocate( Atm%ee2(3,isd_2d:ied_2d+1,jsd_2d:jed_2d+1) )

    ! Center unit vectors:
    allocate( Atm%ec1(3,isd_2d:ied_2d,jsd_2d:jed_2d) )
    allocate( Atm%ec2(3,isd_2d:ied_2d,jsd_2d:jed_2d) )

    ! Edge unit vectors:
    allocate( Atm%ew(3,isd_2d:ied_2d+1,jsd_2d:jed_2d,  2) )
    allocate( Atm%es(3,isd_2d:ied_2d  ,jsd_2d:jed_2d+1,2) )

    ! Edge unit "Normal" vectors: (for omega computation)
    allocate( Atm%en1(3,is_2d:ie_2d,  js_2d:je_2d+1) )   ! E-W edges
    allocate( Atm%en2(3,is_2d:ie_2d+1,js_2d:je_2d  ) )   ! N-S egdes

    allocate ( Atm%cosa_u(isd_2d:ied_2d+1,jsd_2d:jed_2d) )
    allocate ( Atm%sina_u(isd_2d:ied_2d+1,jsd_2d:jed_2d) )
    allocate ( Atm%rsin_u(isd_2d:ied_2d+1,jsd_2d:jed_2d) )

    allocate ( Atm%cosa_v(isd_2d:ied_2d,jsd_2d:jed_2d+1) )
    allocate ( Atm%sina_v(isd_2d:ied_2d,jsd_2d:jed_2d+1) )
    allocate ( Atm%rsin_v(isd_2d:ied_2d,jsd_2d:jed_2d+1) )

    allocate ( Atm%cosa_s(isd_2d:ied_2d,jsd_2d:jed_2d) )    ! cell center
    allocate ( Atm%sina_s(isd_2d:ied_2d,jsd_2d:jed_2d) )    ! cell center

    allocate (  Atm%rsin2(isd_2d:ied_2d,jsd_2d:jed_2d) )    ! cell center


    ! Super (composite) grid:

    !     9---4---8
    !     |       |
    !     1   5   3
    !     |       |
    !     6---2---7

    allocate ( Atm%cos_sg(isd_2d:ied_2d,jsd_2d:jed_2d,9) )
    allocate ( Atm%sin_sg(isd_2d:ied_2d,jsd_2d:jed_2d,9) )

    allocate( Atm%eww(3,4) )
    allocate( Atm%ess(3,4) )

    allocate ( Atm%   fv_olr(is:ie,js:je) )
    allocate ( Atm%fv_abs_sw(is:ie,js:je) )

    allocate ( Atm% zs_g(is:ie, js:je) )
    allocate ( Atm%oro_g(isd:ied, jsd:jed) )
    allocate ( Atm%sgh_g(isd:ied, jsd:jed) )

    allocate ( Atm%phi0(isd:ied  ,jsd:jed  ,npz) )
    allocate (  Atm%ua0(isd:ied  ,jsd:jed  ,npz) )
    allocate (  Atm%va0(isd:ied  ,jsd:jed  ,npz) )

    if (Atm%test_case == 9) then
       allocate( Atm%case9_B(isd:ied,jsd:jed) )
    end if

    if (Atm%test_case == 7) then
       allocate (   Atm%gh_table(4*npy) )
       allocate ( Atm%lats_table(4*npy) )
    end if

    allocate ( Atm%u0(is:ie,  js:je+1,npz) )
    allocate ( Atm%v0(is:ie+1,js:je  ,npz) )
    allocate ( Atm%dp(is:ie,js:je,npz) )
    allocate ( Atm%t0(is:ie,js:je,npz) )

    if (Atm%nested) then

       allocate(Atm%ind_h(isd:ied,jsd:jed,2))
       allocate(Atm%ind_u(isd:ied,jsd:jed+1,2))
       allocate(Atm%ind_v(isd:ied+1,jsd:jed,2))

       allocate(Atm%wt_h(isd:ied,   jsd:jed,  4))
       allocate(Atm%wt_u(isd:ied,   jsd:jed+1,4))
       allocate(Atm%wt_v(isd:ied+1, jsd:jed,  4))

       ns = Atm%nsponge

       if (ie == npx-1) then
          allocate(Atm%h_east(ie+1-ns:ied,jsd:jed,npz))
          allocate(Atm%u_east(ie+1-ns:ied,jsd:jed+1,npz))
          allocate(Atm%v_east(ie+2-ns:ied+1,jsd:jed,npz))
          if (concurrent) then
             allocate(Atm%h_east_t0(ie+1-ns:ied,jsd:jed,npz))
             allocate(Atm%u_east_t0(ie+1-ns:ied,jsd:jed+1,npz))
             allocate(Atm%v_east_t0(ie+2-ns:ied+1,jsd:jed,npz))
          endif
          allocate(Atm%vc_east_t0(ie+1-ns:ied,jsd:jed+1,npz))
          allocate(Atm%uc_east_t0(ie+2-ns:ied+1,jsd:jed,npz))
          allocate(Atm%vc_east_t1(ie+1-ns:ied,jsd:jed+1,npz))
          allocate(Atm%uc_east_t1(ie+2-ns:ied+1,jsd:jed,npz))
          if (ncnst > 0) allocate(Atm%q_east(ie+1-ns:ied,jsd:jed,npz,ncnst))
          if (ncnst > 0) allocate(Atm%q_east_t0(ie+1-ns:ied,jsd:jed,npz,ncnst))
       else
          allocate(Atm%h_east(1,1,npz))
          allocate(Atm%u_east(1,1,npz))
          allocate(Atm%v_east(1,1,npz))
          if (concurrent) then
             allocate(Atm%h_east_t0(1,1,npz))
             allocate(Atm%u_east_t0(1,1,npz))
             allocate(Atm%v_east_t0(1,1,npz))
          endif
          allocate(Atm%uc_east_t0(1,1,npz))
          allocate(Atm%vc_east_t0(1,1,npz))
          allocate(Atm%uc_east_t1(1,1,npz))
          allocate(Atm%vc_east_t1(1,1,npz))
          if (ncnst > 0) allocate(Atm%q_east(1,1,npz,ncnst))          
          if (ncnst > 0) allocate(Atm%q_east_t0(1,1,npz,ncnst))          
       end if

       if (js == 1) then
          allocate(Atm%h_south(isd:ied,jsd:js-1+ns,npz))
          allocate(Atm%u_south(isd:ied,jsd:js-1+ns,npz))
          allocate(Atm%v_south(isd:ied+1,jsd:js-1+ns,npz))
          if (concurrent) then
             allocate(Atm%h_south_t0(isd:ied,jsd:js-1+ns,npz))
             allocate(Atm%u_south_t0(isd:ied,jsd:js-1+ns,npz))
             allocate(Atm%v_south_t0(isd:ied+1,jsd:js-1+ns,npz))
          endif
          allocate(Atm%vc_south_t0(isd:ied,jsd:js-1+ns,npz))
          allocate(Atm%uc_south_t0(isd:ied+1,jsd:js-1+ns,npz))
          allocate(Atm%vc_south_t1(isd:ied,jsd:js-1+ns,npz))
          allocate(Atm%uc_south_t1(isd:ied+1,jsd:js-1+ns,npz))
          if (ncnst > 0) allocate(Atm%q_south(isd:ied,jsd:js-1+ns,npz,ncnst))
          if (ncnst > 0) allocate(Atm%q_south_t0(isd:ied,jsd:js-1+ns,npz,ncnst))
       else
          allocate(Atm%h_south(1,1,npz))
          allocate(Atm%u_south(1,1,npz))
          allocate(Atm%v_south(1,1,npz))
          if (concurrent) then
             allocate(Atm%h_south_t0(1,1,npz))
             allocate(Atm%u_south_t0(1,1,npz))
             allocate(Atm%v_south_t0(1,1,npz))
          endif
          allocate(Atm%uc_south_t0(1,1,npz))
          allocate(Atm%vc_south_t0(1,1,npz))
          allocate(Atm%uc_south_t1(1,1,npz))
          allocate(Atm%vc_south_t1(1,1,npz))
          if (ncnst > 0) allocate(Atm%q_south(1,1,npz,ncnst))
          if (ncnst > 0) allocate(Atm%q_south_t0(1,1,npz,ncnst))
       end if

       if (is == 1) then
          allocate(Atm%h_west(isd:is-1+ns,jsd:jed,npz))
          allocate(Atm%u_west(isd:is-1+ns,jsd:jed+1,npz))
          allocate(Atm%v_west(isd:is-1+ns,jsd:jed,npz))
          if (concurrent) then
             allocate(Atm%h_west_t0(isd:is-1+ns,jsd:jed,npz))
             allocate(Atm%u_west_t0(isd:is-1+ns,jsd:jed+1,npz))
             allocate(Atm%v_west_t0(isd:is-1+ns,jsd:jed,npz))
          endif
          allocate(Atm%vc_west_t0(isd:is-1+ns,jsd:jed+1,npz))
          allocate(Atm%uc_west_t0(isd:is-1+ns,jsd:jed,npz))
          allocate(Atm%vc_west_t1(isd:is-1+ns,jsd:jed+1,npz))
          allocate(Atm%uc_west_t1(isd:is-1+ns,jsd:jed,npz))
          if (ncnst > 0) allocate(Atm%q_west(isd:is-1+ns,jsd:jed,npz,ncnst))
          if (ncnst > 0) allocate(Atm%q_west_t0(isd:is-1+ns,jsd:jed,npz,ncnst))
       else
          allocate(Atm%h_west(1,1,npz))
          allocate(Atm%u_west(1,1,npz))
          allocate(Atm%v_west(1,1,npz))
          if (concurrent) then
             allocate(Atm%h_west_t0(1,1,npz))
             allocate(Atm%u_west_t0(1,1,npz))
             allocate(Atm%v_west_t0(1,1,npz))
          endif
          allocate(Atm%uc_west_t0(1,1,npz))
          allocate(Atm%vc_west_t0(1,1,npz))
          allocate(Atm%uc_west_t1(1,1,npz))
          allocate(Atm%vc_west_t1(1,1,npz))
          if (ncnst > 0) allocate(Atm%q_west(1,1,npz,ncnst))
          if (ncnst > 0) allocate(Atm%q_west_t0(1,1,npz,ncnst))
       end if

       if (je == npy-1) then
          allocate(Atm%h_north(isd:ied,je+1-ns:jed,npz))
          allocate(Atm%u_north(isd:ied,je+2-ns:jed+1,npz))
          allocate(Atm%v_north(isd:ied+1,je+1-ns:jed,npz))
          if (concurrent) then
             allocate(Atm%h_north_t0(isd:ied,je+1-ns:jed,npz))
             allocate(Atm%u_north_t0(isd:ied,je+2-ns:jed+1,npz))
             allocate(Atm%v_north_t0(isd:ied+1,je+1-ns:jed,npz))
          endif
          allocate(Atm%vc_north_t0(isd:ied,je+2-ns:jed+1,npz))
          allocate(Atm%uc_north_t0(isd:ied+1,je+1-ns:jed,npz))
          allocate(Atm%vc_north_t1(isd:ied,je+2-ns:jed+1,npz))
          allocate(Atm%uc_north_t1(isd:ied+1,je+1-ns:jed,npz))
          if (ncnst > 0) allocate(Atm%q_north(isd:ied,je+1-ns:jed,npz,ncnst))
          if (ncnst > 0) allocate(Atm%q_north_t0(isd:ied,je+1-ns:jed,npz,ncnst))
       else
          allocate(Atm%h_north(1,1,npz))
          allocate(Atm%u_north(1,1,npz))
          allocate(Atm%v_north(1,1,npz))
          if (concurrent) then
             allocate(Atm%h_north_t0(1,1,npz))
             allocate(Atm%u_north_t0(1,1,npz))
             allocate(Atm%v_north_t0(1,1,npz))
          endif
          allocate(Atm%uc_north_t0(1,1,npz))
          allocate(Atm%vc_north_t0(1,1,npz))
          allocate(Atm%uc_north_t1(1,1,npz))
          allocate(Atm%vc_north_t1(1,1,npz))
          if (ncnst > 0) allocate(Atm%q_north(1,1,npz,ncnst))
          if (ncnst > 0) allocate(Atm%q_north_t0(1,1,npz,ncnst))
       end if

       if (Atm%twowaynest) allocate(Atm%ind_update_h(Atm%parent_grid%isd:Atm%parent_grid%ied+1,Atm%parent_grid%jsd:Atm%parent_grid%jed+1,2))

#ifndef SW_DYNAMICS

       if (ie == npx-1) then
          allocate(Atm%pt_east(ie+1-ns:ied,jsd:jed,npz))
          if (.not. Atm%hydrostatic) allocate(Atm%w_east(ie+1-ns:ied,jsd:jed+1,npz))
          if (concurrent) then
             allocate(Atm%pt_east_t0(ie+1-ns:ied,jsd:jed,npz))
             if (.not. Atm%hydrostatic) allocate(Atm%w_east_t0(ie+1-ns:ied,jsd:jed+1,npz))
          endif
       else
          allocate(Atm%pt_east(1,1,npz))
          if (.not. Atm%hydrostatic) allocate(Atm%w_east(1,1,npz))
          if (concurrent) then
             allocate(Atm%pt_east_t0(1,1,npz))
             if (.not. Atm%hydrostatic) allocate(Atm%w_east_t0(1,1,npz))
          endif
       end if

       if (js == 1) then
          allocate(Atm%pt_south(isd:ied,jsd:js-1+ns,npz))
          if (.not. Atm%hydrostatic) allocate(Atm%w_south(isd:ied,jsd:js-1+ns,npz))
          if (concurrent) then
             allocate(Atm%pt_south_t0(isd:ied,jsd:js-1+ns,npz))
             if (.not. Atm%hydrostatic) allocate(Atm%w_south_t0(isd:ied,jsd:js-1+ns,npz))
          endif
       else
          allocate(Atm%pt_south(1,1,npz))
          if (.not. Atm%hydrostatic) allocate(Atm%w_south(1,1,npz))
          if (concurrent) then
             allocate(Atm%pt_south_t0(1,1,npz))
             if (.not. Atm%hydrostatic) allocate(Atm%w_south_t0(1,1,npz))
          endif
       end if

       if (is == 1) then
          allocate(Atm%pt_west(isd:is-1+ns,jsd:jed,npz))
          if (.not. Atm%hydrostatic) allocate(Atm%w_west(isd:is-1+ns,jsd:jed+1,npz))
          if (concurrent) then
             allocate(Atm%pt_west_t0(isd:is-1+ns,jsd:jed,npz))
             if (.not. Atm%hydrostatic) allocate(Atm%w_west_t0(isd:is-1+ns,jsd:jed+1,npz))
          endif
       else
          allocate(Atm%pt_west(1,1,npz))
          if (.not. Atm%hydrostatic) allocate(Atm%w_west(1,1,npz))
          if (concurrent) then
             allocate(Atm%pt_west_t0(1,1,npz))
             if (.not. Atm%hydrostatic) allocate(Atm%w_west_t0(1,1,npz))
          endif
       end if

       if (je == npy-1) then
          allocate(Atm%pt_north(isd:ied,je+1-ns:jed,npz))
          if (.not. Atm%hydrostatic) allocate(Atm%w_north(isd:ied,je+2-ns:jed+1,npz))
          if (concurrent) then
             allocate(Atm%pt_north_t0(isd:ied,je+1-ns:jed,npz))
             if (.not. Atm%hydrostatic) allocate(Atm%w_north_t0(isd:ied,je+2-ns:jed+1,npz))
          endif
       else
          allocate(Atm%pt_north(1,1,npz))
          if (.not. Atm%hydrostatic) allocate(Atm%w_north(1,1,npz))
          if (concurrent) then
             allocate(Atm%pt_north_t0(1,1,npz))
             if (.not. Atm%hydrostatic) allocate(Atm%w_north_t0(1,1,npz))
          endif
       end if

#endif

       if (Atm%nestbctype > 1 .and. ncnst > 0) then
          if (is == 1) then
             allocate(Atm%nest_fx_west_accum(js:je,npz,ncnst))
          else
             allocate(Atm%nest_fx_west_accum(1,1,1))
          endif
          if (ie == npx-1) then
             allocate(Atm%nest_fx_east_accum(js:je,npz,ncnst))
          else
             allocate(Atm%nest_fx_east_accum(1,1,1))
          endif

          if (js == 1) then
             allocate(Atm%nest_fx_south_accum(is:ie,npz,ncnst))
          else
             allocate(Atm%nest_fx_south_accum(1,1,1))
          endif
          if (je == npy-1) then
             allocate(Atm%nest_fx_north_accum(is:ie,npz,ncnst))
          else
             allocate(Atm%nest_fx_north_accum(1,1,1))
          endif
       endif

    end if

       if (Atm%grid_type < 4) then
          if (Atm%nested) then
             allocate(Atm%grid_global(1-ng_2d:npx_2d  +ng_2d,1-ng_2d:npy_2d  +ng_2d,2,1))
          else
             allocate(Atm%grid_global(1-ng_2d:npx_2d  +ng_2d,1-ng_2d:npy_2d  +ng_2d,2,1:6))
          endif
       end if


#ifdef MARS_GCM
    allocate( Atm%tmars(is:ie,js:je,npz) )
#endif

#ifndef SW_DYNAMICS
    allocate( Atm%u_dt(isd:ied,jsd:jed, npz) )
    allocate( Atm%v_dt(isd:ied,jsd:jed, npz) )
    allocate( Atm%t_dt(is:ie,js:je, npz) )
    allocate( Atm%q_dt(is:ie,js:je, npz, Atm%nt_prog) )
    allocate( Atm%p_full(is:ie,js:je,npz) )
    allocate( Atm%z_full(is:ie,js:je,npz) )
    allocate( Atm%p_half(is:ie,js:je,npz+1) )
    allocate( Atm%z_half(is:ie,js:je,npz+1) )
! For phys_filter:
    if ( Atm%tq_filter ) then
         allocate (  Atm%t_phys(isd:ied,jsd:jed,npz) )
         allocate (  Atm%q_phys(isd:ied,jsd:jed,npz,Atm%nt_prog) )
    endif

#endif


    Atm%isd = isd_in
    Atm%ied = ied_in
    Atm%jsd = jsd_in
    Atm%jed = jed_in

    Atm%is = is_in
    Atm%ie = ie_in
    Atm%js = js_in
    Atm%je = je_in

    Atm%isc = Atm%is
    Atm%iec = Atm%ie
    Atm%jsc = Atm%js
    Atm%jec = Atm%je

    Atm%npx = npx_in
    Atm%npy = npy_in
    Atm%npz = npz_in
    Atm%ndims = ndims_in

    Atm%allocated = .true.
    if (dummy) Atm%dummy = .true.
    
  end subroutine allocate_fv_atmos_type

  subroutine deallocate_fv_atmos_type(Atm)

    implicit none
    type(fv_atmos_type), intent(INOUT) :: Atm

    if (.not.Atm%allocated) return

    deallocate (    Atm%u )
    deallocate (    Atm%v )
    deallocate (   Atm%pt )
    deallocate ( Atm%delp )
    deallocate (    Atm%q )
    deallocate (   Atm%ps )
    deallocate (   Atm%pe )
    deallocate (   Atm%pk )
    deallocate ( Atm%peln )
    deallocate (  Atm%pkz )
    deallocate ( Atm%phis )
    deallocate ( Atm%omga )
    deallocate (   Atm%ua )
    deallocate (   Atm%va )
    deallocate (   Atm%uc )
    deallocate (   Atm%vc )
    deallocate ( Atm%mfx )
    deallocate ( Atm%mfy )
    deallocate (  Atm%cx )
    deallocate (  Atm%cy )
    deallocate (  Atm%ak )
    deallocate (  Atm%bk )

    deallocate ( Atm%u_srf )
    deallocate ( Atm%v_srf )
    if( Atm%fv_land ) deallocate ( Atm%sgh )
    deallocate ( Atm%oro )

    deallocate ( Atm%w )
    deallocate ( Atm%delz  )
    deallocate ( Atm%ze0   )

    deallocate ( Atm% area )   ! Cell Centered
    deallocate ( Atm%rarea )   ! Cell Centered
    
    deallocate ( Atm% area_c )  ! Cell Corners
    deallocate ( Atm%rarea_c )  ! Cell Corners
    
    deallocate ( Atm% dx )
    deallocate ( Atm%rdx )
    deallocate ( Atm% dy )
    deallocate ( Atm%rdy )
    
    deallocate ( Atm% dxc )
    deallocate ( Atm%rdxc )
    deallocate ( Atm% dyc )
    deallocate ( Atm%rdyc )
    
    deallocate ( Atm% dxa )
    deallocate ( Atm%rdxa )
    deallocate ( Atm% dya )
    deallocate ( Atm%rdya )
    
    deallocate ( Atm%grid  )
    deallocate ( Atm%agrid )
    deallocate ( Atm%sina )   ! SIN(angle of intersection)
    deallocate ( Atm%cosa )   ! COS(angle of intersection)
    
    deallocate ( Atm%  e1 )
    deallocate ( Atm%  e2 )




    deallocate (Atm%iinta, &
         Atm%jinta,  &
         Atm%iintb, &
         Atm%jintb )

#ifndef NO_GRID_G
    deallocate ( Atm%grid_g )
#endif

    deallocate ( Atm%edge_s )
    deallocate ( Atm%edge_n )
    deallocate ( Atm%edge_w )
    deallocate ( Atm%edge_e )

    deallocate ( Atm%edge_vect_s )
    deallocate ( Atm%edge_vect_n )
    deallocate ( Atm%edge_vect_w )
    deallocate ( Atm%edge_vect_e )

    deallocate ( Atm%ex_s )
    deallocate ( Atm%ex_n )
    deallocate ( Atm%ex_w )
    deallocate ( Atm%ex_e )


    ! For diveregnce damping:
    deallocate (  Atm%divg_u )
    deallocate (  Atm%divg_v )

    deallocate (  Atm%z11 )
    deallocate (  Atm%z12 )
    deallocate (  Atm%z21 )
    deallocate (  Atm%z22 )

    deallocate (  Atm%a11 )
    deallocate (  Atm%a12 )
    deallocate (  Atm%a21 )
    deallocate (  Atm%a22 )
    deallocate ( Atm%vlon )
    deallocate ( Atm%vlat )
    ! Coriolis parameters:
    deallocate ( Atm%f0 )
    deallocate ( Atm%fC )

    ! Corner unit vectors:
    deallocate( Atm%ee1 )
    deallocate( Atm%ee2 )

    ! Center unit vectors:
    deallocate( Atm%ec1 )
    deallocate( Atm%ec2 )

    ! Edge unit vectors:
    deallocate( Atm%ew )
    deallocate( Atm%es )

    ! Edge unit "Normal" vectors: (for omega computation)
    deallocate( Atm%en1 )   ! E-W edges
    deallocate( Atm%en2 )   ! N-S egdes

    deallocate ( Atm%cosa_u )
    deallocate ( Atm%sina_u )
    deallocate ( Atm%rsin_u )

    deallocate ( Atm%cosa_v )
    deallocate ( Atm%sina_v )
    deallocate ( Atm%rsin_v )

    deallocate ( Atm%cosa_s )    ! cell center
    deallocate ( Atm%sina_s )    ! cell center

    deallocate (  Atm%rsin2 )    ! cell center


    ! Super (composite) grid:

    !     9---4---8
    !     |       |
    !     1   5   3
    !     |       |
    !     6---2---7

    deallocate ( Atm%cos_sg )
    deallocate ( Atm%sin_sg )

    deallocate( Atm%eww )
    deallocate( Atm%ess )

    deallocate ( Atm%   fv_olr )
    deallocate ( Atm%fv_abs_sw )

    deallocate ( Atm% zs_g )
    deallocate ( Atm%oro_g )
    deallocate ( Atm%sgh_g )

    deallocate ( Atm%phi0 )
    deallocate (  Atm%ua0 )
    deallocate (  Atm%va0 )

    if (Atm%test_case == 9) then
       deallocate( Atm%case9_B )
    end if

    if (Atm%test_case == 7) then
       deallocate (   Atm%gh_table )
       deallocate ( Atm%lats_table )
    end if

    deallocate(Atm%u0)
    deallocate ( Atm%v0)
    deallocate ( Atm%dp)
    deallocate ( Atm%t0)

    if (Atm%nested) then
       deallocate(Atm%ind_h)
       deallocate(Atm%ind_u)
       deallocate(Atm%ind_v)

       deallocate(Atm%wt_h)
       deallocate(Atm%wt_u)
       deallocate(Atm%wt_v)

          deallocate(Atm%h_east)
          deallocate(Atm%u_east)
          deallocate(Atm%v_east)
          deallocate(Atm%uc_east_t0)
          deallocate(Atm%vc_east_t0)
          deallocate(Atm%uc_east_t1)
          deallocate(Atm%vc_east_t1)
          if (Atm%ncnst > 0) deallocate(Atm%q_east)
          if (Atm%ncnst > 0) deallocate(Atm%q_east_t0)

          deallocate(Atm%h_south)
          deallocate(Atm%u_south)
          deallocate(Atm%v_south)
          deallocate(Atm%uc_south_t0)
          deallocate(Atm%vc_south_t0)
          deallocate(Atm%uc_south_t1)
          deallocate(Atm%vc_south_t1)
          if (Atm%ncnst > 0) deallocate(Atm%q_south)
          if (Atm%ncnst > 0) deallocate(Atm%q_south_t0)

          deallocate(Atm%h_west)
          deallocate(Atm%u_west)
          deallocate(Atm%v_west)
          deallocate(Atm%uc_west_t0)
          deallocate(Atm%vc_west_t0)
          deallocate(Atm%uc_west_t1)
          deallocate(Atm%vc_west_t1)
          if (Atm%ncnst > 0) deallocate(Atm%q_west)
          if (Atm%ncnst > 0) deallocate(Atm%q_west_t0)

          deallocate(Atm%h_north)
          deallocate(Atm%u_north)
          deallocate(Atm%v_north)
          deallocate(Atm%uc_north_t0)
          deallocate(Atm%vc_north_t0)
          deallocate(Atm%uc_north_t1)
          deallocate(Atm%vc_north_t1)
          if (Atm%ncnst > 0) deallocate(Atm%q_north)
          if (Atm%ncnst > 0) deallocate(Atm%q_north_t0)

          if (Atm%ncnst > 0 .and. Atm%nestbctype > 1) deallocate( Atm%nest_fx_west_accum, Atm%nest_fx_east_accum, Atm%nest_fx_south_accum, Atm%nest_fx_north_accum) 

       if (Atm%twowaynest) deallocate(Atm%ind_update_h)

#ifndef SW_DYNAMICS

          deallocate(Atm%pt_east    )
          if (.not. Atm%hydrostatic) deallocate(Atm%w_east      )

          deallocate(Atm%pt_south    )
          if (.not. Atm%hydrostatic) deallocate(Atm%w_south      )

          deallocate(Atm%pt_west  )
          if (.not. Atm%hydrostatic) deallocate(Atm%w_west   )

          deallocate(Atm%pt_north    )
          if (.not. Atm%hydrostatic) deallocate(Atm%w_north      )

#endif
    end if

#ifdef MARS_GCM
    deallocate( Atm%tmars)
#endif

#ifndef SW_DYNAMICS
    deallocate( Atm%u_dt)
    deallocate( Atm%v_dt)
    deallocate( Atm%t_dt)
    deallocate( Atm%q_dt)
    deallocate( Atm%p_full)
    deallocate( Atm%z_full)
    deallocate( Atm%p_half)
    deallocate( Atm%z_half)
    if ( Atm%tq_filter ) then
         deallocate (  Atm%t_phys)
         deallocate (  Atm%q_phys)
    endif
#endif

       if (Atm%grid_type < 4) then
          deallocate(Atm%grid_global)
       end if

       Atm%allocated = .false.

  end subroutine deallocate_fv_atmos_type

end module fv_arrays_mod
