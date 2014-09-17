module fv_arrays_mod
#include <fms_platform.h>
 use mpp_domains_mod,  only: domain2d
  use fms_io_mod,       only: restart_file_type
  use time_manager_mod, only: time_type
  use horiz_interp_type_mod, only:  horiz_interp_type
  use mpp_domains_mod, only : nest_domain_type
  use mpp_mod, only: mpp_broadcast
  public

  !Several 'auxiliary' structures are introduced here. These are for
  ! the internal use by certain modules, and although fv_atmos_type
  !  contains one of each of these structures all memory management
  !   is performed by the module in question.

  integer, parameter:: max_step = 1000
  type fv_diag_type


 integer ::id_ps, id_slp, id_ua, id_va, id_pt, id_omga, id_vort,  &
           id_tm, id_pv, id_zsurf, id_oro, id_sgh, id_divg, id_w, &
           id_ke, id_te, id_zs, id_ze, id_mq, id_vorts, id_us, id_vs,    &
           id_tq, id_rh, id_c15, id_c25, id_c35, id_c45,          &
                         id_f15, id_f25, id_f35, id_f45,          &
           id_ppt, id_ts, id_pmask, id_pmaskv2,                   &
           id_delp, id_delz, id_zratio, id_ws, id_iw, id_lw

! Selected p-level fields from 3D variables:
 integer :: id_vort850, id_w850, id_x850, id_srh,  &
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


  !fv_grid_type is made up of grid-dependent information from fv_grid_tools and fv_grid_utils.
  ! It should not contain any user options (that goes in a different structure) nor data which
  ! is altered outside of those two modules.
  type fv_grid_type
     real, allocatable, dimension(:,:,:) :: grid, agrid
     real, allocatable, dimension(:,:) :: area, area_c
     real, allocatable, dimension(:,:) :: rarea, rarea_c     

  real, allocatable, dimension(:,:) :: sina, cosa
  real, allocatable, dimension(:,:,:) :: e1,e2
  real, allocatable, dimension(:,:) :: dx, dy
  real, allocatable, dimension(:,:) :: dxc, dyc
  real, allocatable, dimension(:,:) :: dxa, dya
  real, allocatable, dimension(:,:) :: rdx, rdy
  real, allocatable, dimension(:,:) :: rdxc, rdyc
  real, allocatable, dimension(:,:) :: rdxa, rdya

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

     integer, dimension(:,:,:), allocatable :: iinta, jinta, iintb, jintb
  
     !Scalar data
     
     integer :: npx_g, npy_g, ntiles_g ! global domain

     real :: global_area
     logical :: g_sum_initialized = .false. !Not currently used but can be useful
     logical:: sw_corner, se_corner, ne_corner, nw_corner

     real :: da_min, da_max, da_min_c, da_max_c

     real  :: acapN, acapS
     real  :: globalarea  ! total Global Area
     
     logical :: latlon = .false.
     logical :: cubed_sphere = .false.
     logical :: have_south_pole = .false.
     logical :: have_north_pole = .false.
     logical :: stretched_grid = .false.

     logical :: square_domain = .false.


     !! Convenience pointers

     integer, pointer :: grid_type
     logical, pointer :: nested

  end type fv_grid_type

  type fv_flags_type

     !! FOR EACH VARIABLE IN FV_FLAGS:
     !! 1. Must be defined here:
     !! 2. Must be broadcast in fv_atmos_data
     !! 3. If a namelist entry, a pointer must
     !!    be defined and associated in fv_control
     !! 4. Must NOT appear in fv_current_grid_mod.
     !!    (this module will soon be removed)
     !! 5. Must be referenced through Atm%flagstruct,
     !!    not Atm%, unless a convenience
     !!    pointer is defined

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
   real    :: scale_z = 0.   ! diff_z = scale_z**2 * 0.25
   real    :: w_max = 75.    ! max w (m/s) threshold for hydostatiic adjustment 
   real    :: z_min = 0.05   ! min ratio of dz_nonhydrostatic/dz_hydrostatic

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

   logical :: do_sat_adj= .false.   ! 
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

   logical :: use_logp = .false.

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
   real    :: rf_cutoff = 30.E2       ! cutoff pressure level for RF
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
   logical :: do_reed_physics = .false.
   logical :: reed_cond_only = .false.
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
   logical :: phys_hydrostatic = .true.  ! heating/cooling term from the physics is hydrostatic
   logical :: hybrid_z    = .false.      ! use hybrid_z for remapping
   logical :: Make_NH     = .false.      ! Initialize (w, delz) from hydro restart file 
   logical :: make_hybrid_z  = .false.   ! transform hydrostatic eta-coord IC into non-hydrostatic hybrid_z
   real    :: add_noise = -1.            !Amplitude of random noise added upon model startup; <=0 means no noise added

   integer :: a2b_ord = 4    ! order for interpolation from A to B Grid (corners)
   integer :: c2l_ord = 4    ! order for interpolation from D to lat-lon A winds for phys & output

  real :: dx_const = 1000.    ! spatial resolution for double periodic boundary configuration [m]
  real :: dy_const = 1000.
  real :: deglat=15.
  !The following deglat_*, deglon_* options are not used.
  real :: deglon_start = -30., deglon_stop = 30., &  ! boundaries of latlon patch
          deglat_start = -30., deglat_stop = 30.

  !Convenience pointers
  integer, pointer :: grid_number

  
  !integer, pointer :: test_case
  !real,    pointer :: alpha

  end type fv_flags_type

  type fv_nest_BC_type_3D

     !!! CLEANUP: could we have pointers to np[xyz], nest_domain, and the index/weight arrays?

     real, allocatable, dimension(:,:,:) :: west_t1, east_t1, south_t1, north_t1
     real, allocatable, dimension(:,:,:) :: west_t0, east_t0, south_t0, north_t0

     integer :: istag, jstag

     logical :: allocated = .false.

  end type fv_nest_BC_type_3D

  type fv_nest_BC_type_4D

     real, allocatable, dimension(:,:,:,:) :: west_t1, east_t1, south_t1, north_t1
     real, allocatable, dimension(:,:,:,:) :: west_t0, east_t0, south_t0, north_t0

     integer :: istag, jstag

     logical :: allocated = .false.

  end type fv_nest_BC_type_4D

  type fv_nest_type

!nested grid flags:

     integer :: refinement = 3  !Refinement wrt parent

     integer :: parent_tile = 1     !Tile (of cubed sphere) in which nested grid lies 
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
     

     type(nest_domain_type) :: nest_domain !Structure holding link from this grid to its parent
     type(nest_domain_type), allocatable :: nest_domain_all(:)

     !Interpolation arrays for grid nesting
     integer, allocatable, dimension(:,:,:) :: ind_h, ind_u, ind_v
     real, allocatable, dimension(:,:,:) :: wt_h, wt_u, wt_v
     integer, allocatable, dimension(:,:,:) :: ind_update_h

     !These arrays are not allocated by allocate_fv_atmos_type; but instead
     !allocated for all grids, regardless of whether the grid is
     !on a PE of a concurrent run.
     logical, allocatable, dimension(:) :: child_grids

     logical :: parent_proc, child_proc
     logical :: parent_of_twoway = .false.
   
     !These are for time-extrapolated BCs
     type(fv_nest_BC_type_3D) :: delp_BC, u_BC, v_BC, uc_BC, vc_BC
     type(fv_nest_BC_type_3D), allocatable, dimension(:) :: q_BC
#ifndef SW_DYNAMICS
     type(fv_nest_BC_type_3D) :: pt_BC, w_BC, delz_BC
#endif

     !These are for tracer flux BCs
     logical :: do_flux_BCs, do_2way_flux_BCs !For a parent grid; determine whether there is a need to send BCs
     type(restart_file_type) :: BCfile_ne, BCfile_sw

  end type fv_nest_type

  interface allocate_fv_nest_BC_type
     module procedure allocate_fv_nest_BC_type_3D
     module procedure allocate_fv_nest_BC_type_3D_Atm
  end interface

  interface deallocate_fv_nest_BC_type
     module procedure deallocate_fv_nest_BC_type_3D
  end interface

  type fv_grid_bounds_type

     integer :: is,  ie,  js,  je
     integer :: isd, ied, jsd, jed
     integer :: isc, iec, jsc, jec

     integer :: ng

  end type fv_grid_bounds_type

  type fv_atmos_type

     logical :: allocated = .false.
     logical :: dummy = .false.
     integer :: grid_number = 1

     !Timestep-related variables.

     type(time_type) :: Time_init, Time, Run_length, Time_end, Time_step_atmos

     logical :: grid_active = .true. !Always active for now

     !This is kept here instead of in neststruct% simply for convenience
     type(fv_atmos_type), pointer :: parent_grid _NULL

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
#ifdef PKC
    real, _ALLOCATABLE :: pkc (:,:,:)   _NULL  ! finite-volume edge pk
#endif

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

   integer :: ks

! Accumulated Mass flux arrays
    real, _ALLOCATABLE ::  mfx(:,:,:)  _NULL
    real, _ALLOCATABLE ::  mfy(:,:,:)  _NULL
! Accumulated Courant number arrays
    real, _ALLOCATABLE ::  cx(:,:,:)  _NULL
    real, _ALLOCATABLE ::  cy(:,:,:)  _NULL

    type(fv_flags_type) :: flagstruct
    
    !! Convenience pointers
    integer, pointer :: npx, npy, npz, ncnst, ng

     integer, allocatable, dimension(:) :: pelist

     type(fv_grid_bounds_type) :: bd

     type(domain2D) :: domain
#if defined(SPMD)

     type(domain2D) :: domain_for_coupler ! domain used in coupled model with halo = 1.

     integer :: num_contact, npes_per_tile, tile, npes_this_grid
     integer :: layout(2), io_layout(2) = (/ 1,1 /)

#endif
     !These do not actually belong to the grid, but to the process
     !integer :: masterproc
     !integer :: gid 

!!!!!!!!!!!!!!!!
! From fv_grid_tools
!!!!!!!!!!!!!!!!


     real    :: ptop

  type(fv_grid_type) :: gridstruct
  

!!!!!!!!!!!!!!!!
!fv_diagnostics!
!!!!!!!!!!!!!!!!

     type(fv_diag_type) :: idiag

!!!!!!!!!!!!!!
! From fv_io !
!!!!!!!!!!!!!!
     type(restart_file_type) :: Fv_restart, SST_restart, Fv_tile_restart, &
          Rsf_restart, Mg_restart, Lnd_restart, Tra_restart

     type(fv_nest_type) :: neststruct

     !Hold on to coarse-grid global grid, so we don't have to waste processor time getting it again when starting to do grid nesting
     real, allocatable, dimension(:,:,:,:) :: grid_global

  integer :: atmos_axes(4)


  end type fv_atmos_type

!---- version number -----
  character(len=128) :: version = '$Id: fv_arrays.F90,v 20.0.2.1 2014/02/07 00:17:12 Rusty.Benson Exp $'
  character(len=128) :: tagname = '$Name: tikal_201409 $'

contains

  subroutine allocate_fv_atmos_type(Atm, isd_in, ied_in, jsd_in, jed_in, is_in, ie_in, js_in, je_in, &
       npx_in, npy_in, npz_in, ndims_in, ncnst_in, ng_in, dummy, alloc_2d)

    !WARNING: Before calling this routine, be sure to have set up the
    ! proper domain parameters from the namelists (as is done in
    ! fv_control.F90)

    implicit none
    type(fv_atmos_type), intent(INOUT), target :: Atm
    integer, intent(IN) :: isd_in, ied_in, jsd_in, jed_in, is_in, ie_in, js_in, je_in
    integer, intent(IN) :: npx_in, npy_in, npz_in, ndims_in, ncnst_in, ng_in
    logical, intent(IN) :: dummy, alloc_2d
    integer:: isd, ied, jsd, jed, is, ie, js, je
    integer:: npx, npy, npz, ndims, ncnst, ng

    !For 2D utility arrays
    integer:: isd_2d, ied_2d, jsd_2d, jed_2d, is_2d, ie_2d, js_2d, je_2d
    integer:: npx_2d, npy_2d, npz_2d, ndims_2d, ncnst_2d, ng_2d

    integer :: ns, n

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

!This should be set up in fv_mp_mod
!!$    Atm%bd%isd = isd_in
!!$    Atm%bd%ied = ied_in
!!$    Atm%bd%jsd = jsd_in
!!$    Atm%bd%jed = jed_in
!!$
!!$    Atm%bd%is = is_in
!!$    Atm%bd%ie = ie_in
!!$    Atm%bd%js = js_in
!!$    Atm%bd%je = je_in
!!$
!!$    Atm%bd%isc = Atm%bd%is
!!$    Atm%bd%iec = Atm%bd%ie
!!$    Atm%bd%jsc = Atm%bd%js
!!$    Atm%bd%jec = Atm%bd%je

    Atm%bd%ng  = ng

    !Convenience pointers
    Atm%npx => Atm%flagstruct%npx
    Atm%npy => Atm%flagstruct%npy
    Atm%npz => Atm%flagstruct%npz
    Atm%ncnst => Atm%flagstruct%ncnst

    Atm%ng => Atm%bd%ng

!!$    Atm%npx = npx_in
!!$    Atm%npy = npy_in
!!$    Atm%npz = npz_in
    Atm%flagstruct%ndims = ndims_in

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
#ifdef PKC
    allocate (  Atm%pkc(isd:ied,jsd:jed,npz+1) )
#endif

    allocate ( Atm%u_srf(is:ie,js:je) )
    allocate ( Atm%v_srf(is:ie,js:je) )

    if ( Atm%flagstruct%fv_land ) then
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
    if ( Atm%flagstruct%hydrostatic ) then
       allocate (    Atm%w(1, 1  ,1) )
       allocate ( Atm%delz(1, 1  ,1) )
       allocate (  Atm%ze0(1, 1  ,1) )
    else
       allocate (    Atm%w(isd:ied, jsd:jed  ,npz  ) )
       allocate ( Atm%delz(isd:ied, jsd:jed  ,npz) )
       if( Atm%flagstruct%hybrid_z ) then
          allocate (  Atm%ze0(is:ie, js:je ,npz+1) )
       else
          allocate (  Atm%ze0(1, 1  ,1) )
       endif
       !         allocate ( mono(isd:ied, jsd:jed, npz))
    endif

    allocate ( Atm%gridstruct% area(isd_2d:ied_2d  ,jsd_2d:jed_2d  ) )   ! Cell Centered
    allocate ( Atm%gridstruct%rarea(isd_2d:ied_2d  ,jsd_2d:jed_2d  ) )   ! Cell Centered
    
    allocate ( Atm%gridstruct% area_c(isd_2d:ied_2d+1,jsd_2d:jed_2d+1) )  ! Cell Corners
    allocate ( Atm%gridstruct%rarea_c(isd_2d:ied_2d+1,jsd_2d:jed_2d+1) )  ! Cell Corners
    
    allocate ( Atm%gridstruct% dx(isd_2d:ied_2d  ,jsd_2d:jed_2d+1) )
    allocate ( Atm%gridstruct%rdx(isd_2d:ied_2d  ,jsd_2d:jed_2d+1) )
    allocate ( Atm%gridstruct% dy(isd_2d:ied_2d+1,jsd_2d:jed_2d  ) )
    allocate ( Atm%gridstruct%rdy(isd_2d:ied_2d+1,jsd_2d:jed_2d  ) )
    
    allocate ( Atm%gridstruct% dxc(isd_2d:ied_2d+1,jsd_2d:jed_2d  ) )
    allocate ( Atm%gridstruct%rdxc(isd_2d:ied_2d+1,jsd_2d:jed_2d  ) )
    allocate ( Atm%gridstruct% dyc(isd_2d:ied_2d  ,jsd_2d:jed_2d+1) )
    allocate ( Atm%gridstruct%rdyc(isd_2d:ied_2d  ,jsd_2d:jed_2d+1) )
    
    allocate ( Atm% gridstruct%dxa(isd_2d:ied_2d  ,jsd_2d:jed_2d  ) )
    allocate ( Atm%gridstruct%rdxa(isd_2d:ied_2d  ,jsd_2d:jed_2d  ) )
    allocate ( Atm% gridstruct%dya(isd_2d:ied_2d  ,jsd_2d:jed_2d  ) )
    allocate ( Atm%gridstruct%rdya(isd_2d:ied_2d  ,jsd_2d:jed_2d  ) )
    
    allocate ( Atm%gridstruct%grid (isd_2d:ied_2d+1,jsd_2d:jed_2d+1,1:ndims_2d) )
    allocate ( Atm%gridstruct%agrid(isd_2d:ied_2d  ,jsd_2d:jed_2d  ,1:ndims_2d) )
    allocate ( Atm%gridstruct%sina(isd_2d:ied_2d+1,jsd_2d:jed_2d+1) )   ! SIN(angle of intersection)
    allocate ( Atm%gridstruct%rsina(is_2d:ie_2d+1,js_2d:je_2d+1) )      ! Why is the size different?
    allocate ( Atm%gridstruct%cosa(isd_2d:ied_2d+1,jsd_2d:jed_2d+1) )   ! COS(angle of intersection)
    
    allocate ( Atm%gridstruct%  e1(3,isd_2d:ied_2d+1,jsd_2d:jed_2d+1) )
    allocate ( Atm%gridstruct%  e2(3,isd_2d:ied_2d+1,jsd_2d:jed_2d+1) )

    allocate (Atm%gridstruct%iinta(4, isd_2d:ied_2d ,jsd_2d:jed_2d), &
         Atm%gridstruct%jinta(4, isd_2d:ied_2d ,jsd_2d:jed_2d),  &
         Atm%gridstruct%iintb(4, is_2d:ie_2d+1 ,js_2d:je_2d+1), &
         Atm%gridstruct%jintb(4, is_2d:ie_2d+1 ,js_2d:je_2d+1) )

    allocate ( Atm%gridstruct%edge_s(npx_2d) )
    allocate ( Atm%gridstruct%edge_n(npx_2d) )
    allocate ( Atm%gridstruct%edge_w(npy_2d) )
    allocate ( Atm%gridstruct%edge_e(npy_2d) )

    allocate ( Atm%gridstruct%edge_vect_s(isd_2d:ied_2d) )
    allocate ( Atm%gridstruct%edge_vect_n(isd_2d:ied_2d) )
    allocate ( Atm%gridstruct%edge_vect_w(jsd_2d:jed_2d) )
    allocate ( Atm%gridstruct%edge_vect_e(jsd_2d:jed_2d) )

    allocate ( Atm%gridstruct%ex_s(npx_2d) )
    allocate ( Atm%gridstruct%ex_n(npx_2d) )
    allocate ( Atm%gridstruct%ex_w(npy_2d) )
    allocate ( Atm%gridstruct%ex_e(npy_2d) )


    ! For diveregnce damping:
    allocate (  Atm%gridstruct%divg_u(isd_2d:ied_2d,  jsd_2d:jed_2d+1) )
    allocate (  Atm%gridstruct%divg_v(isd_2d:ied_2d+1,jsd_2d:jed_2d) )

    allocate (  Atm%gridstruct%z11(is_2d-1:ie_2d+1,js_2d-1:je_2d+1) )
    allocate (  Atm%gridstruct%z12(is_2d-1:ie_2d+1,js_2d-1:je_2d+1) )
    allocate (  Atm%gridstruct%z21(is_2d-1:ie_2d+1,js_2d-1:je_2d+1) )
    allocate (  Atm%gridstruct%z22(is_2d-1:ie_2d+1,js_2d-1:je_2d+1) )

    allocate (  Atm%gridstruct%a11(is_2d-1:ie_2d+1,js_2d-1:je_2d+1) )
    allocate (  Atm%gridstruct%a12(is_2d-1:ie_2d+1,js_2d-1:je_2d+1) )
    allocate (  Atm%gridstruct%a21(is_2d-1:ie_2d+1,js_2d-1:je_2d+1) )
    allocate (  Atm%gridstruct%a22(is_2d-1:ie_2d+1,js_2d-1:je_2d+1) )
    allocate ( Atm%gridstruct%vlon(is_2d-2:ie_2d+2,js_2d-2:je_2d+2,3) )
    allocate ( Atm%gridstruct%vlat(is_2d-2:ie_2d+2,js_2d-2:je_2d+2,3) )
    ! Coriolis parameters:
    allocate ( Atm%gridstruct%f0(isd_2d:ied_2d  ,jsd_2d:jed_2d  ) )
    allocate ( Atm%gridstruct%fC(isd_2d:ied_2d+1,jsd_2d:jed_2d+1) )

    ! Corner unit vectors:
    allocate( Atm%gridstruct%ee1(3,isd_2d:ied_2d+1,jsd_2d:jed_2d+1) )
    allocate( Atm%gridstruct%ee2(3,isd_2d:ied_2d+1,jsd_2d:jed_2d+1) )

    ! Center unit vectors:
    allocate( Atm%gridstruct%ec1(3,isd_2d:ied_2d,jsd_2d:jed_2d) )
    allocate( Atm%gridstruct%ec2(3,isd_2d:ied_2d,jsd_2d:jed_2d) )

    ! Edge unit vectors:
    allocate( Atm%gridstruct%ew(3,isd_2d:ied_2d+1,jsd_2d:jed_2d,  2) )
    allocate( Atm%gridstruct%es(3,isd_2d:ied_2d  ,jsd_2d:jed_2d+1,2) )

    ! Edge unit "Normal" vectors: (for omega computation)
    allocate( Atm%gridstruct%en1(3,is_2d:ie_2d,  js_2d:je_2d+1) )   ! E-W edges
    allocate( Atm%gridstruct%en2(3,is_2d:ie_2d+1,js_2d:je_2d  ) )   ! N-S egdes

    allocate ( Atm%gridstruct%cosa_u(isd_2d:ied_2d+1,jsd_2d:jed_2d) )
    allocate ( Atm%gridstruct%sina_u(isd_2d:ied_2d+1,jsd_2d:jed_2d) )
    allocate ( Atm%gridstruct%rsin_u(isd_2d:ied_2d+1,jsd_2d:jed_2d) )

    allocate ( Atm%gridstruct%cosa_v(isd_2d:ied_2d,jsd_2d:jed_2d+1) )
    allocate ( Atm%gridstruct%sina_v(isd_2d:ied_2d,jsd_2d:jed_2d+1) )
    allocate ( Atm%gridstruct%rsin_v(isd_2d:ied_2d,jsd_2d:jed_2d+1) )

    allocate ( Atm%gridstruct%cosa_s(isd_2d:ied_2d,jsd_2d:jed_2d) )    ! cell center
    allocate ( Atm%gridstruct%sina_s(isd_2d:ied_2d,jsd_2d:jed_2d) )    ! cell center

    allocate (  Atm%gridstruct%rsin2(isd_2d:ied_2d,jsd_2d:jed_2d) )    ! cell center


    ! Super (composite) grid:

    !     9---4---8
    !     |       |
    !     1   5   3
    !     |       |
    !     6---2---7

    allocate ( Atm%gridstruct%cos_sg(isd_2d:ied_2d,jsd_2d:jed_2d,9) )
    allocate ( Atm%gridstruct%sin_sg(isd_2d:ied_2d,jsd_2d:jed_2d,9) )

    allocate( Atm%gridstruct%eww(3,4) )
    allocate( Atm%gridstruct%ess(3,4) )

    if (Atm%neststruct%nested) then

       allocate(Atm%neststruct%ind_h(isd:ied,jsd:jed,2))
       allocate(Atm%neststruct%ind_u(isd:ied,jsd:jed+1,2))
       allocate(Atm%neststruct%ind_v(isd:ied+1,jsd:jed,2))

       allocate(Atm%neststruct%wt_h(isd:ied,   jsd:jed,  4))
       allocate(Atm%neststruct%wt_u(isd:ied,   jsd:jed+1,4))
       allocate(Atm%neststruct%wt_v(isd:ied+1, jsd:jed,  4))

       ns = Atm%neststruct%nsponge

       call allocate_fv_nest_BC_type(Atm%neststruct%delp_BC,Atm,ns,0,0,dummy)
       call allocate_fv_nest_BC_type(Atm%neststruct%u_BC,Atm,ns,0,1,dummy)
       call allocate_fv_nest_BC_type(Atm%neststruct%v_BC,Atm,ns,1,0,dummy)
       call allocate_fv_nest_BC_type(Atm%neststruct%uc_BC,Atm,ns,1,0,dummy)
       call allocate_fv_nest_BC_type(Atm%neststruct%vc_BC,Atm,ns,0,1,dummy)

       if (ncnst > 0) then
          allocate(Atm%neststruct%q_BC(ncnst))
          do n=1,ncnst
             call allocate_fv_nest_BC_type(Atm%neststruct%q_BC(n),Atm,ns,0,0,dummy)
          enddo
       endif

#ifndef SW_DYNAMICS

       call allocate_fv_nest_BC_type(Atm%neststruct%pt_BC,Atm,ns,0,0,dummy)
       if (.not.Atm%flagstruct%hydrostatic) then
          call allocate_fv_nest_BC_type(Atm%neststruct%w_BC,Atm,ns,0,0,dummy)
          call allocate_fv_nest_BC_type(Atm%neststruct%delz_BC,Atm,ns,0,0,dummy)
       endif

#endif

       if (Atm%neststruct%twowaynest) allocate(&
            Atm%neststruct%ind_update_h( &
              Atm%parent_grid%bd%isd:Atm%parent_grid%bd%ied+1, &
              Atm%parent_grid%bd%jsd:Atm%parent_grid%bd%jed+1,2))

    end if

       if (Atm%flagstruct%grid_type < 4) then
          if (Atm%neststruct%nested) then
             allocate(Atm%grid_global(1-ng_2d:npx_2d  +ng_2d,1-ng_2d:npy_2d  +ng_2d,2,1))
          else
             allocate(Atm%grid_global(1-ng_2d:npx_2d  +ng_2d,1-ng_2d:npy_2d  +ng_2d,2,1:6))
          endif
       end if


    Atm%allocated = .true.
    if (dummy) Atm%dummy = .true.
    
  end subroutine allocate_fv_atmos_type

  subroutine deallocate_fv_atmos_type(Atm)

    implicit none
    type(fv_atmos_type), intent(INOUT) :: Atm

    integer :: n

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
#ifdef PKC
    deallocate (  Atm%pkc )
#endif
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
    if( Atm%flagstruct%fv_land ) deallocate ( Atm%sgh )
    deallocate ( Atm%oro )

    deallocate ( Atm%w )
    deallocate ( Atm%delz  )
    deallocate ( Atm%ze0   )

    deallocate ( Atm%gridstruct% area )   ! Cell Centered
    deallocate ( Atm%gridstruct%rarea )   ! Cell Centered
    
    deallocate ( Atm%gridstruct% area_c )  ! Cell Corners
    deallocate ( Atm%gridstruct%rarea_c )  ! Cell Corners
    
    deallocate ( Atm%gridstruct% dx )
    deallocate ( Atm%gridstruct%rdx )
    deallocate ( Atm%gridstruct% dy )
    deallocate ( Atm%gridstruct%rdy )
    
    deallocate ( Atm%gridstruct% dxc )
    deallocate ( Atm%gridstruct%rdxc )
    deallocate ( Atm%gridstruct% dyc )
    deallocate ( Atm%gridstruct%rdyc )
    
    deallocate ( Atm%gridstruct% dxa )
    deallocate ( Atm%gridstruct%rdxa )
    deallocate ( Atm%gridstruct% dya )
    deallocate ( Atm%gridstruct%rdya )
    
    deallocate ( Atm%gridstruct%grid  )
    deallocate ( Atm%gridstruct%agrid )
    deallocate ( Atm%gridstruct%sina )   ! SIN(angle of intersection)
    deallocate ( Atm%gridstruct%cosa )   ! COS(angle of intersection)
    
    deallocate ( Atm%gridstruct%  e1 )
    deallocate ( Atm%gridstruct%  e2 )




    deallocate (Atm%gridstruct%iinta, &
         Atm%gridstruct%jinta,  &
         Atm%gridstruct%iintb, &
         Atm%gridstruct%jintb )

    deallocate ( Atm%gridstruct%edge_s )
    deallocate ( Atm%gridstruct%edge_n )
    deallocate ( Atm%gridstruct%edge_w )
    deallocate ( Atm%gridstruct%edge_e )

    deallocate ( Atm%gridstruct%edge_vect_s )
    deallocate ( Atm%gridstruct%edge_vect_n )
    deallocate ( Atm%gridstruct%edge_vect_w )
    deallocate ( Atm%gridstruct%edge_vect_e )

    deallocate ( Atm%gridstruct%ex_s )
    deallocate ( Atm%gridstruct%ex_n )
    deallocate ( Atm%gridstruct%ex_w )
    deallocate ( Atm%gridstruct%ex_e )


    ! For diveregnce damping:
    deallocate (  Atm%gridstruct%divg_u )
    deallocate (  Atm%gridstruct%divg_v )

    deallocate (  Atm%gridstruct%z11 )
    deallocate (  Atm%gridstruct%z12 )
    deallocate (  Atm%gridstruct%z21 )
    deallocate (  Atm%gridstruct%z22 )

    deallocate (  Atm%gridstruct%a11 )
    deallocate (  Atm%gridstruct%a12 )
    deallocate (  Atm%gridstruct%a21 )
    deallocate (  Atm%gridstruct%a22 )
    deallocate ( Atm%gridstruct%vlon )
    deallocate ( Atm%gridstruct%vlat )
    ! Coriolis parameters:
    deallocate ( Atm%gridstruct%f0 )
    deallocate ( Atm%gridstruct%fC )

    ! Corner unit vectors:
    deallocate( Atm%gridstruct%ee1 )
    deallocate( Atm%gridstruct%ee2 )

    ! Center unit vectors:
    deallocate( Atm%gridstruct%ec1 )
    deallocate( Atm%gridstruct%ec2 )

    ! Edge unit vectors:
    deallocate( Atm%gridstruct%ew )
    deallocate( Atm%gridstruct%es )

    ! Edge unit "Normal" vectors: (for omega computation)
    deallocate( Atm%gridstruct%en1 )   ! E-W edges
    deallocate( Atm%gridstruct%en2 )   ! N-S egdes

    deallocate ( Atm%gridstruct%cosa_u )
    deallocate ( Atm%gridstruct%sina_u )
    deallocate ( Atm%gridstruct%rsin_u )

    deallocate ( Atm%gridstruct%cosa_v )
    deallocate ( Atm%gridstruct%sina_v )
    deallocate ( Atm%gridstruct%rsin_v )

    deallocate ( Atm%gridstruct%cosa_s )    ! cell center
    deallocate ( Atm%gridstruct%sina_s )    ! cell center

    deallocate (  Atm%gridstruct%rsin2 )    ! cell center


    ! Super (composite) grid:

    !     9---4---8
    !     |       |
    !     1   5   3
    !     |       |
    !     6---2---7

    deallocate ( Atm%gridstruct%cos_sg )
    deallocate ( Atm%gridstruct%sin_sg )

    deallocate( Atm%gridstruct%eww )
    deallocate( Atm%gridstruct%ess )

    if (Atm%neststruct%nested) then
       deallocate(Atm%neststruct%ind_h)
       deallocate(Atm%neststruct%ind_u)
       deallocate(Atm%neststruct%ind_v)

       deallocate(Atm%neststruct%wt_h)
       deallocate(Atm%neststruct%wt_u)
       deallocate(Atm%neststruct%wt_v)

       call deallocate_fv_nest_BC_type(Atm%neststruct%delp_BC)
       call deallocate_fv_nest_BC_type(Atm%neststruct%u_BC)
       call deallocate_fv_nest_BC_type(Atm%neststruct%v_BC)
       call deallocate_fv_nest_BC_type(Atm%neststruct%uc_BC)
       call deallocate_fv_nest_BC_type(Atm%neststruct%vc_BC)
       if (allocated(Atm%neststruct%q_BC)) then
          do n=1,size(Atm%neststruct%q_BC)
             call deallocate_fv_nest_BC_type(Atm%neststruct%q_BC(n))
          enddo
       endif

#ifndef SW_DYNAMICS
       call deallocate_fv_nest_BC_type(Atm%neststruct%pt_BC)
       if (.not.Atm%flagstruct%hydrostatic) then
          call deallocate_fv_nest_BC_type(Atm%neststruct%w_BC)
          call deallocate_fv_nest_BC_type(Atm%neststruct%delz_BC)
       endif
#endif


       if (Atm%neststruct%twowaynest) deallocate(Atm%neststruct%ind_update_h)

    end if

    if (Atm%flagstruct%grid_type < 4) then
       deallocate(Atm%grid_global)
    end if
    
    Atm%allocated = .false.

  end subroutine deallocate_fv_atmos_type


subroutine allocate_fv_nest_BC_type_3D_Atm(BC,Atm,ns,istag,jstag,dummy)

  type(fv_nest_BC_type_3D), intent(INOUT) :: BC
  type(fv_atmos_type), intent(IN) :: Atm
  integer, intent(IN) :: ns, istag, jstag
  logical, intent(IN) :: dummy

  integer :: is, ie, js, je, isd, ied, jsd, jed, npx, npy, npz, ng

  if (BC%allocated) return

  is = Atm%bd%is
  ie = Atm%bd%ie
  js = Atm%bd%js
  je = Atm%bd%je

  isd = Atm%bd%isd
  ied = Atm%bd%ied
  jsd = Atm%bd%jsd
  jed = Atm%bd%jed

  npx = Atm%npx
  npy = Atm%npy
  npz = Atm%npz

  ng = Atm%ng

  call allocate_fv_nest_BC_type_3D(BC,is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz,ng,ns,istag,jstag,dummy)


end subroutine allocate_fv_nest_BC_type_3D_Atm

subroutine allocate_fv_nest_BC_type_3D(BC,is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz,ng,ns,istag,jstag,dummy)

  type(fv_nest_BC_type_3D), intent(INOUT) :: BC
  integer, intent(IN) :: ns, istag, jstag
  logical, intent(IN) :: dummy

  integer, intent(IN) :: is, ie, js, je, isd, ied, jsd, jed, npx, npy, npz, ng

  if (BC%allocated) return


  if (ie == npx-1 .and. .not. dummy) then
     allocate(BC%east_t1(ie+1-ns+istag:ied+istag,jsd:jed+jstag,npz))
     allocate(BC%east_t0(ie+1-ns+istag:ied+istag,jsd:jed+jstag,npz))
  else
     allocate(BC%east_t1(1,1,npz))
     allocate(BC%east_t0(1,1,npz))
  end if

  if (js == 1 .and. .not. dummy) then
     allocate(BC%south_t1(isd:ied+istag,jsd:js-1+ns,npz))
     allocate(BC%south_t0(isd:ied+istag,jsd:js-1+ns,npz))
  else
     allocate(BC%south_t1(1,1,npz))
     allocate(BC%south_t0(1,1,npz))
  end if

  if (is == 1 .and. .not. dummy) then
     allocate(BC%west_t1(isd:is-1+ns,jsd:jed+jstag,npz))
     allocate(BC%west_t0(isd:is-1+ns,jsd:jed+jstag,npz))
  else
     allocate(BC%west_t1(1,1,npz))
     allocate(BC%west_t0(1,1,npz))
  end if

  if (je == npy-1 .and. .not. dummy) then
     allocate(BC%north_t1(isd:ied+istag,je+1-ns+jstag:jed+jstag,npz))
     allocate(BC%north_t0(isd:ied+istag,je+1-ns+jstag:jed+jstag,npz))
  else
     allocate(BC%north_t1(1,1,npz))
     allocate(BC%north_t0(1,1,npz))
  end if

  BC%allocated = .true.

end subroutine allocate_fv_nest_BC_type_3D

subroutine deallocate_fv_nest_BC_type_3d(BC)

  type(fv_nest_BC_type_3d) :: BC

  if (.not. BC%allocated) return

     deallocate(BC%north_t1)
     deallocate(BC%south_t1)
     deallocate(BC%west_t1)
     deallocate(BC%east_t1)

  if (allocated(BC%north_t0)) then
     deallocate(BC%north_t0)
     deallocate(BC%south_t0)
     deallocate(BC%west_t0)
     deallocate(BC%east_t0)
  endif

  BC%allocated = .false.

end subroutine deallocate_fv_nest_BC_type_3d

  subroutine broadcast_fv_atmos_data(Atm, fromproc, to_pelist)

    type(fv_atmos_type) :: Atm
    integer, intent(IN) :: fromproc, to_pelist(:)

    character(len=128) :: sendchar(4)

     call mpp_broadcast( Atm%grid_number, fromproc, to_pelist)
     call mpp_broadcast( Atm%grid_active, fromproc, to_pelist)

     !NOTE: mpp_broadcast can only send ARRAYS of character strings, and not individual
     !character strings. The code here is a work-around.
     sendchar(1) = Atm%flagstruct%grid_name
     sendchar(2) = Atm%flagstruct%grid_file
     sendchar(3) = Atm%flagstruct%res_latlon_dynamics
     sendchar(4) = Atm%flagstruct%res_latlon_tracers
     call mpp_broadcast( sendchar, size(sendchar), fromproc, to_pelist)
     Atm%flagstruct%grid_name = sendchar(1)
     Atm%flagstruct%grid_file = sendchar(2)
     Atm%flagstruct%res_latlon_dynamics = sendchar(3)
     Atm%flagstruct%res_latlon_tracers = sendchar(4)
     call mpp_broadcast( Atm%flagstruct%grid_type, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%hord_mt, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%kord_mt, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%kord_wz, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%hord_vt, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%hord_tm, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%hord_dp, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%kord_tm, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%hord_tr, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%kord_tr, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%scale_z, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%w_max, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%z_min, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%nord, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%dddmp, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%d2_bg, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%d4_bg, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%vtdm4, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%d2_bg_k1, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%d2_bg_k2, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%d2_divg_max_k1, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%d2_divg_max_k2, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%damp_k_k1, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%damp_k_k2, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%n_zs_filter, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%nord_zs_filter, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%do_sat_adj, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%no_dycore, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%replace_w, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%convert_ke, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%do_vort_damp, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%use_old_omega, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%beta, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%n_sponge, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%d_ext, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%nwat, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%warm_start, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%inline_q, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%n_sponge, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%d_ext, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%nwat, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%warm_start, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%inline_q, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%shift_fac, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%do_schmidt, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%stretch_fac, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%target_lat, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%target_lon, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%reset_eta, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%p_fac, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%a_imp, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%n_split, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%m_split, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%k_split, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%q_split, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%print_freq, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%npx, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%npy, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%npz, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%npz_rst, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%ncnst, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%pnats, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%ntiles, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%ndims, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%nf_omega, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%fv_sg_adj, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%na_init, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%p_ref, fromproc, to_pelist)
#ifdef MARS_GCM
     call mpp_broadcast( Atm%flagstruct%reference_sfc_pres, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%sponge_damp, fromproc, to_pelist)
#endif
     call mpp_broadcast( Atm%flagstruct%dry_mass, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%p_ref, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%nt_prog, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%nt_phys, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%tau_h2o, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%d_con, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%consv_te, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%tau, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%rf_center, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%rf_cutoff, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%filter_phys, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%dwind_2d, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%breed_vortex_inline, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%range_warn, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%fill, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%fill_dp, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%fill_wz, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%non_ortho, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%adiabatic, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%moist_phys, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%do_Held_Suarez, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%reproduce_sum, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%adjust_dry_mass, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%fv_debug, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%srf_init, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%mountain, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%remap_t, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%z_tracer, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%old_divg_damp, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%fv_land, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%nudge, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%ncep_ic, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%fv_diag_ic, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%external_ic, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%hydrostatic, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%phys_hydrostatic, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%hybrid_z, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%Make_NH, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%make_hybrid_z, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%add_noise, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%a2b_ord, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%c2l_ord, fromproc, to_pelist)
     call mpp_broadcast( Atm%ks, fromproc, to_pelist)
     call mpp_broadcast( Atm%bd%ng, fromproc, to_pelist)
     call mpp_broadcast( Atm%num_contact, fromproc, to_pelist)
     call mpp_broadcast( Atm%npes_per_tile, fromproc, to_pelist)
     call mpp_broadcast( Atm%tile, fromproc, to_pelist)
     call mpp_broadcast( Atm%npes_this_grid, fromproc, to_pelist)
     call mpp_broadcast( Atm%gridstruct%square_domain, fromproc, to_pelist)
     call mpp_broadcast( Atm%layout, 2, fromproc, to_pelist)
     call mpp_broadcast( Atm%io_layout, 2, fromproc, to_pelist)
     call mpp_broadcast( Atm%gridstruct%latlon, fromproc, to_pelist)
     call mpp_broadcast( Atm%gridstruct%cubed_sphere, fromproc, to_pelist)
     call mpp_broadcast( Atm%gridstruct%have_south_pole, fromproc, to_pelist)
     call mpp_broadcast( Atm%gridstruct%have_north_pole, fromproc, to_pelist)
     call mpp_broadcast( Atm%gridstruct%stretched_grid, fromproc, to_pelist)
     call mpp_broadcast( Atm%gridstruct%npx_g, fromproc, to_pelist)
     call mpp_broadcast( Atm%gridstruct%npy_g, fromproc, to_pelist)
     call mpp_broadcast( Atm%gridstruct%ntiles_g, fromproc, to_pelist)
     call mpp_broadcast( Atm%gridstruct%acapN, fromproc, to_pelist)
     call mpp_broadcast( Atm%gridstruct%acapS, fromproc, to_pelist)
     call mpp_broadcast( Atm%gridstruct%globalarea, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%dx_const, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%dy_const, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%deglon_start, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%deglon_stop, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%deglat_start, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%deglat_stop, fromproc, to_pelist)
     call mpp_broadcast( Atm%flagstruct%deglat, fromproc, to_pelist)
     call mpp_broadcast( Atm%ptop, fromproc, to_pelist)
     call mpp_broadcast( Atm%gridstruct%da_min, fromproc, to_pelist)
     call mpp_broadcast( Atm%gridstruct%da_max, fromproc, to_pelist)
     call mpp_broadcast( Atm%gridstruct%da_min_c, fromproc, to_pelist)
     call mpp_broadcast( Atm%gridstruct%da_max_c, fromproc, to_pelist)
     call mpp_broadcast( Atm%neststruct%parent_tile, fromproc, to_pelist)
     call mpp_broadcast( Atm%neststruct%refinement, fromproc, to_pelist)
     call mpp_broadcast( Atm%neststruct%nested, fromproc, to_pelist)
     call mpp_broadcast( Atm%neststruct%nestbctype, fromproc, to_pelist)
     call mpp_broadcast( Atm%neststruct%nsponge, fromproc, to_pelist)
     call mpp_broadcast( Atm%neststruct%nestupdate, fromproc, to_pelist)
     call mpp_broadcast( Atm%neststruct%twowaynest, fromproc, to_pelist)
     call mpp_broadcast( Atm%neststruct%ioffset, fromproc, to_pelist)
     call mpp_broadcast( Atm%neststruct%joffset, fromproc, to_pelist)
     call mpp_broadcast( Atm%neststruct%s_weight, fromproc, to_pelist)

   end subroutine broadcast_fv_atmos_data


end module fv_arrays_mod
