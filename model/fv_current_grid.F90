module fv_current_grid_mod

  use fv_arrays_mod

  implicit none
  public

     type(fv_atmos_type), pointer :: current_Atm
     integer, pointer :: grid_number

     !Timestep-related variables.
     !Each grid should have its own set of timing utilities
     real , pointer :: dt_grid
     integer , pointer :: refinement ! number of this grid's timesteps per model-wide dt_atmos
     type(time_type) , pointer :: Time_init, Time, Run_length, Time_end, Time_step_atmos

     logical , pointer :: grid_active 

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
     real,  pointer :: u(:,:,:)      ! D grid zonal wind (m/s)
     real,  pointer :: v(:,:,:)      ! D grid meridional wind (m/s)
     real,  pointer :: pt(:,:,:)     ! temperature (K)
     real,  pointer :: delp(:,:,:)   ! pressure thickness (pascal)
     real,  pointer :: q(:,:,:,:)    ! specific humidity and constituents

     !----------------------
     ! non-hydrostatic state:
     !----------------------------------------------------------------------
     real,  pointer ::     w(:,:,:)    ! cell center vertical wind (m/s)
     real,  pointer ::  delz(:,:,:)    ! layer thickness (meters)
     real,  pointer ::   ze0(:,:,:)    ! height at layer edges for remapping

     !-----------------------------------------------------------------------
     ! Auxilliary pressure arrays:
     ! The 5 vars below can be re-computed from delp and ptop.
     !-----------------------------------------------------------------------
     ! dyn_aux:
     real,  pointer :: ps (:,:)        ! Surface pressure (pascal)
     real,  pointer :: pe (:,:,: )     ! edge pressure (pascal)
     real,  pointer :: pk  (:,:,:)     ! pe**cappa
     real,  pointer :: peln(:,:,:)     ! ln(pe)
     real,  pointer :: pkz (:,:,:)     ! finite-volume mean pk

     ! For phys coupling:
     real,  pointer :: u_srf(:,:)      ! Surface u-wind
     real,  pointer :: v_srf(:,:)      ! Surface v-wind
     real,  pointer :: sgh(:,:)        ! Terrain standard deviation
     real,  pointer :: oro(:,:)        ! land fraction (1: all land; 0: all water)
     real,  pointer :: ts(:,:)         ! skin temperature (sst) from NCEP/GFS (K) -- tile

     !-----------------------------------------------------------------------
     ! Others:
     !-----------------------------------------------------------------------
     real,  pointer :: phis(:,:)       ! Surface geopotential (g*Z_surf)
     real,  pointer :: omga(:,:,:)     ! Vertical pressure velocity (pa/s)
     real,  pointer :: ua(:,:,:)       ! (ua, va) are mostly used as the A grid winds
     real,  pointer :: va(:,:,:)     
     real,  pointer :: uc(:,:,:)       ! (uc, vc) are mostly used as the C grid winds
     real,  pointer :: vc(:,:,:)     

     real,  pointer :: ak(:)  
     real,  pointer :: bk(:)  

     ! Accumulated Mass flux arrays
     real,  pointer ::  mfx(:,:,:)  
     real,  pointer ::  mfy(:,:,:)  
     ! Accumulated Courant number arrays
     real,  pointer ::  cx(:,:,:)  
     real,  pointer ::  cy(:,:,:)  


!!!!!!!!!!!!!!!!!!
     ! From fv_control!
!!!!!!!!!!!!!!!!!!
   character(len=80) , pointer :: grid_name
   character(len=120), pointer :: grid_file
   integer, pointer :: grid_type
   integer , pointer :: hord_mt 
   integer , pointer :: kord_mt 
   integer , pointer :: kord_wz 
   integer , pointer :: hord_vt 
   integer , pointer :: hord_tm 
   integer , pointer :: hord_dp 
   integer , pointer :: kord_tm 
   integer , pointer :: hord_tr 
                             !12: Huynh 2nd constraint (Lin 2004) +
                             !    positive definite (Lin & Rood 1996); slower
                             !>12: positive definite only (Lin & Rood 1996); fastest
   integer , pointer :: kord_tr 
   real    , pointer :: scale_z 
   real    , pointer :: w_max 
   real    , pointer :: z_min 

   integer , pointer :: nord
                             ! Alternative setting for high-res: nord
   real    , pointer :: dddmp 
                             ! for C90 or lower: 0.2
   real    , pointer :: d2_bg 
   real    , pointer :: d4_bg 
                             ! for stability, d4_bg must be <
   real    , pointer :: vtdm4 
   real    , pointer :: d2_bg_k1 
   real    , pointer :: d2_bg_k2 
   real    , pointer :: d2_divg_max_k1 
   real    , pointer :: d2_divg_max_k2 
   real    , pointer :: damp_k_k1 
   real    , pointer :: damp_k_k2 
   integer , pointer ::    n_zs_filter
   integer , pointer :: nord_zs_filter

   logical , pointer :: no_dycore 
   logical , pointer :: replace_w 
                                    ! this is useful for getting a good initial estimate of w
                                    ! suugest usage: from a hydrostatic IC, set hydrostatic 
                                    !                replace_w 
   logical , pointer :: convert_ke 
   logical , pointer :: do_vort_damp 
   logical , pointer :: use_old_omega 
! PG off centering:
   real    , pointer :: beta  
#ifdef SW_DYNAMICS
   integer , pointer :: n_sponge 
   real    , pointer :: d_ext 
   integer , pointer :: nwat  
   logical , pointer :: warm_start 
   logical , pointer :: inline_q 
#else
   integer , pointer :: n_sponge 
   real    , pointer :: d_ext 
   integer , pointer :: nwat  
   logical , pointer :: warm_start 
                             ! Set to .F. if cold_start is desired (including terrain generation)
   logical , pointer :: inline_q 
#endif
   real , pointer :: shift_fac   
   logical , pointer :: do_schmidt 
   real , pointer :: stretch_fac 
   real , pointer :: target_lat  
   real , pointer :: target_lon  

   logical , pointer :: reset_eta 
   real    , pointer :: p_fac
   real    , pointer :: a_imp
   integer , pointer :: n_split 
                             ! Default 
   integer , pointer :: m_split 
   integer , pointer :: k_split 

   integer , pointer :: q_split 
   integer , pointer :: print_freq 
                             ! 0: off
                             ! positive n: every n hours
                             ! negative n: every time step

   integer , pointer :: npx                     ! Number of Grid Points in X- dir
   integer , pointer :: npy                     ! Number of Grid Points in Y- dir
   integer , pointer :: npz                     ! Number of Vertical Levels
   integer , pointer :: npz_rst 
                                      ! 0: no change (default)
   integer , pointer :: ncnst 
   integer , pointer :: pnats 
   integer , pointer :: ntiles                  ! Number or tiles that make up the Grid 
   integer , pointer :: nf_omega  
   integer , pointer :: fv_sg_adj 
                                      ! Relaxzation time  scale (sec) if positive
   integer , pointer :: na_init 
#ifdef MARS_GCM
   real    , pointer :: p_ref 
   real    , pointer :: reference_sfc_pres 
   real    , pointer :: sponge_damp
   real    , pointer :: dry_mass 
#else
   real    , pointer :: p_ref 
   real    , pointer :: dry_mass 
#endif
   integer , pointer :: nt_prog 
   integer , pointer :: nt_phys 
   real    , pointer :: tau_h2o 

   real    , pointer :: d_con 
   real    , pointer :: consv_te 
   real    , pointer :: tau 
   real    , pointer :: rf_center 
                                      ! 0: use the top layer center
                                      ! > 0, [Pascal]
   logical , pointer :: tq_filter 
   logical , pointer :: filter_phys 
   logical , pointer :: dwind_2d 
   logical , pointer :: breed_vortex_inline 
   logical , pointer :: range_warn 
   logical , pointer :: fill 
   logical , pointer :: fill_dp 
   logical , pointer :: fill_wz 
   logical , pointer :: non_ortho 
   logical , pointer :: adiabatic 
   logical , pointer :: moist_phys 
   logical , pointer :: do_Held_Suarez 
   logical , pointer :: reproduce_sum 
   logical , pointer :: adjust_dry_mass 
   logical , pointer :: fv_debug  
   logical , pointer :: srf_init  
   logical , pointer :: mountain  
   logical , pointer :: remap_t  
   logical , pointer :: z_tracer 
                                      ! time split; use this if tracer number is huge and/or
                                      ! high resolution (nsplt > 1)

   logical , pointer :: old_divg_damp 
                                      ! defined in a previous revision
                                      ! old_values:
                                      !    d2_bg_k1 
                                      !    d2_divg_max_k1 
                                      !    damp_k_k1 
                                      ! current_values:
                                      !    d2_bg_k1 
                                      !    d2_divg_max_k1 
                                      !    damp_k_k1 
   logical , pointer :: master

   logical , pointer :: fv_land 
   logical , pointer :: nudge 
   logical , pointer :: ncep_ic 
   logical , pointer :: fv_diag_ic 
   logical , pointer :: external_ic 
   character(len=128) , pointer :: res_latlon_dynamics
   character(len=128) , pointer :: res_latlon_tracers 
   logical , pointer :: hydrostatic 
   logical , pointer :: phys_hydrostatic 
   logical , pointer :: hybrid_z    
   logical , pointer :: Make_NH     
   logical , pointer :: make_hybrid_z  

   integer , pointer :: a2b_ord 
   integer , pointer :: c2l_ord 

   integer, pointer :: ks
   integer, pointer :: ndims

!!!!!!!!!!!!!!!!!!
! From fv_mp_mod !
!!!!!!!!!!!!!!!!!!

     integer , pointer :: ng !this SHOULD be a constant, but structure elements are not allowed to be constants
     type(domain2D) , pointer :: domain
#if defined(SPMD)

     type(domain2D) , pointer :: domain_for_coupler ! domain used in coupled model with halo = 1.

     integer , pointer :: num_contact, npes_per_tile, tile, npes_this_grid
     logical , pointer :: square_domain
     integer , pointer :: npes_x, npes_y
     integer, pointer :: layout(:), io_layout(:)

#endif
     integer, pointer :: is, ie, js, je
     integer, pointer :: isd, ied, jsd, jed
     integer, pointer :: isc, iec, jsc, jec

     type(nest_domain_type) , pointer :: nest_domain !Structure holding link from this grid to its parent
     logical , pointer :: this_proc_sends, this_proc_recvs

!!!!!!!!!!!!!!!!
! From fv_grid_tools
!!!!!!!!!!!!!!!!


  real, pointer :: csFac 
  real, pointer            :: zeta 
  real, pointer    :: stretch               ! Optional stretching factor for the grid
  logical, pointer :: dxdy_area 
  logical, pointer :: latlon 
  logical, pointer :: cubed_sphere 
  logical, pointer :: double_periodic 
  logical, pointer :: latlon_patch 
  logical, pointer :: latlon_strip 
  logical, pointer :: channel 
  logical, pointer :: have_south_pole 
  logical, pointer :: have_north_pole 
  integer, pointer :: interpOrder 
  logical, pointer :: debug_message_size 
  logical, pointer :: write_grid_char_file 
  logical, pointer :: stretched_grid 
  ! grid descriptors

  ! Horizontal
  integer, pointer :: npx_g, npy_g, npz_g, ntiles_g ! global domain
#ifndef NO_GRID_G
  real, pointer, dimension(:,:,:) :: grid_g
#endif
  real, pointer, dimension(:,:,:) :: grid, agrid
  real, pointer, dimension(:,:) :: area, area_c
  real, pointer, dimension(:,:) :: sina, cosa
  real, pointer, dimension(:,:,:) :: e1,e2
  real, pointer, dimension(:,:) :: dx, dy
  real, pointer, dimension(:,:) :: dxc, dyc
  real, pointer, dimension(:,:) :: dxa, dya
  real, pointer, dimension(:,:) :: rarea, rarea_c
  real, pointer, dimension(:,:) :: rdx, rdy
  real, pointer, dimension(:,:) :: rdxc, rdyc
  real, pointer, dimension(:,:) :: rdxa, rdya
  real, pointer  :: acapN, acapS
  real, pointer  :: globalarea  ! total Global Area
  real, pointer :: cose(:,:)
  real, pointer :: cosp(:,:)
  real, pointer :: acosp(:,:)

  integer, pointer, dimension(:,:,:) :: iinta, jinta, iintb, jintb

  real, pointer :: dx_const
  real, pointer :: dy_const
  real, pointer :: deglon_start, deglon_stop, &  ! boundaries of latlon patch
          deglat_start, deglat_stop

!!!!!!!!!!!!!!!!
!fv_diagnostics!
!!!!!!!!!!!!!!!!

     type(fv_diag_type), pointer :: idiag
     integer, pointer :: steps
     real(kind=4), pointer :: efx(:), mtq(:), efx_nest(:)
     real(kind=4), pointer :: efx_sum, efx_sum_nest,      mtq_sum


!!!!!!!!!!!!!!!!!!!!!!
     ! From fv_grid_utils !
!!!!!!!!!!!!!!!!!!!!!!

     ! Scalars:
     real, pointer :: edge_s(:)
     real, pointer :: edge_n(:)
     real, pointer :: edge_w(:)
     real, pointer :: edge_e(:)
     ! Vector:
     real, pointer :: edge_vect_s(:)
     real, pointer :: edge_vect_n(:)
     real, pointer :: edge_vect_w(:)
     real, pointer :: edge_vect_e(:)
     ! scalar:
     real, pointer :: ex_s(:)
     real, pointer :: ex_n(:)
     real, pointer :: ex_w(:)
     real, pointer :: ex_e(:)

     ! divergence Damping:
     real, pointer :: divg_u(:,:), divg_v(:,:)    !
     ! Cubed_2_latlon:
     real, pointer :: a11(:,:)
     real, pointer :: a12(:,:)
     real, pointer :: a21(:,:)
     real, pointer :: a22(:,:)
     ! latlon_2_cubed:
     real, pointer :: z11(:,:)
     real, pointer :: z12(:,:)
     real, pointer :: z21(:,:)
     real, pointer :: z22(:,:)

     real, pointer :: global_area
     real, pointer:: stretch_factor
     logical, pointer :: g_sum_initialized !Not currently used but can be useful
     logical, pointer :: gnomonic_grid
     logical, pointer :: sw_corner, se_corner, ne_corner, nw_corner
     real, pointer :: cosa_u(:,:)
     real, pointer :: cosa_v(:,:)
     real, pointer :: cosa_s(:,:)
     real, pointer :: sina_s(:,:)
     real, pointer :: sina_u(:,:)
     real, pointer :: sina_v(:,:)
     real, pointer :: rsin_u(:,:)
     real, pointer :: rsin_v(:,:)
     real, pointer ::  rsina(:,:)
     real, pointer ::  rsin2(:,:)
     real, pointer :: ee1(:,:,:)
     real, pointer :: ee2(:,:,:)
     real, pointer :: ec1(:,:,:)
     real, pointer :: ec2(:,:,:)
     real, pointer :: ew(:,:,:,:)
     real, pointer :: es(:,:,:,:)


     !- 3D Super grid to contain all geometrical factors --
     ! the 3rd dimension is 9
     real, pointer :: sin_sg(:,:,:)
     real, pointer :: cos_sg(:,:,:)
     !--------------------------------------------------

     ! Unit Normal vectors at cell edges:
     real, pointer :: en1(:,:,:)
     real, pointer :: en2(:,:,:)

     ! Extended Cubed cross-edge winds
     real, pointer :: eww(:,:)
     real, pointer :: ess(:,:)

     ! Unit vectors for lat-lon grid
     real, pointer :: vlon(:,:,:), vlat(:,:,:)
     real, pointer :: fC(:,:), f0(:,:)

     real, pointer :: deglat

     real   , pointer :: ptop
     real, pointer :: da_min, da_max, da_min_c, da_max_c


!!!!!!!!!!!!!!
! From fv_sg !
!!!!!!!!!!!!!!

     real, pointer :: fv_olr(:,:), fv_abs_sw(:,:)
     integer, pointer :: irad

!!!!!!!!!!!!!!!!!!!!
! From fv_surf_map !
!!!!!!!!!!!!!!!!!!!!

     real, pointer :: zs_g(:,:)
     real, pointer :: sgh_g(:,:), oro_g(:,:)

!!!!!!!!!!!!!!!!!!!
! From test_cases !
!!!!!!!!!!!!!!!!!!!

     integer, pointer :: test_case !If not specified, this at least tells us
                                 !that the nested-grid halo topography needs
                                 !to be filled by interpolating from the coarse grid.
     ! alpha = angle of axis rotation about the poles
     real  , pointer :: alpha
     ! Ubar = initial wind speed parameter
     real  , pointer :: Ubar
     ! gh0 = initial surface height parameter
     real  , pointer :: gh0

     !  case 9 parameters
     real  , pointer :: case9_B(:,:)
     real  , pointer :: AofT(:)


     !  Validating fields used in statistics
     real  , pointer :: phi0(:,:,:) ! Validating Field
     real  , pointer :: ua0(:,:,:)  ! Validating U-Wind
     real  , pointer :: va0(:,:,:)  ! Validating V-Windfms_io_exit, get_tile_string, &

     real  , pointer :: gh_table(:), lats_table(:)
     logical, pointer :: gh_initialized

     !  Initial Conservation statistics ; total mass ; enstrophy ; energy
     real  , pointer :: tmass_orig
     real  , pointer :: tvort_orig
     real  , pointer :: tener_orig

!!!!!!!!!!!!!!!!
! From fv_phys !
!!!!!!!!!!!!!!!!
     !This is exclusively for nudging

     logical , pointer :: nudge_initialized
     real, pointer :: u0(:,:,:), v0(:,:,:), t0(:,:,:), dp(:,:,:)

!!!!!!!!!!!!!
! From hswf !
!!!!!!!!!!!!!

#ifdef MARS_GCM
     logical , pointer :: tmars_initialized
     real,   dimension(:,:,:), pointer :: tmars
#endif

!!!!!!!!!!!!!!
! From fv_io !
!!!!!!!!!!!!!!
     type(restart_file_type) , pointer :: Fv_restart, SST_restart, Fv_tile_restart, &
          Rsf_restart, Mg_restart, Lnd_restart, Tra_restart

!!!!!!!!!!!!!!!!!!!!!!!
! Nesting information !
!!!!!!!!!!!!!!!!!!!!!!!

     type(fv_atmos_type), pointer :: parent_grid 
     integer , pointer :: parent_tile     !Tile (of cubed sphere) in which nested grid lies 
     integer , pointer :: nest_refinement  !Refinement wrt parent
     logical , pointer :: nested
     integer , pointer :: nestbctype
     integer , pointer :: nsponge 
     integer , pointer :: nestupdate       
     logical , pointer :: twowaynest 
     integer , pointer :: ioffset, joffset !Position of nest within parent grid
     integer , pointer :: nest_timestep !Counter for nested-grid timesteps
     integer , pointer :: tracer_nest_timestep !Counter for nested-grid timesteps
     real    , pointer :: s_weight !sponge weight
     logical , pointer :: first_step
     integer , pointer :: refinement_of_global
     integer , pointer :: npx_global
     logical, dimension(:), pointer :: child_grids
     integer, dimension(:), pointer :: pelist

     logical, pointer :: parent_proc, child_proc

     !Interpolation arrays for grid nesting
     integer,  dimension(:,:,:) , pointer :: ind_h, ind_u, ind_v
     real,  dimension(:,:,:) , pointer :: wt_h, wt_u, wt_v
     integer,  dimension(:,:,:) , pointer :: ind_update_h

     !Arrays to hold data for time interpolation
     real,  dimension(:,:,:) , pointer :: h_west, h_east, h_north, h_south
     real,  dimension(:,:,:) , pointer :: u_west, u_east, u_north, u_south
     real,  dimension(:,:,:) , pointer :: v_west, v_east, v_north, v_south
     real,  dimension(:,:,:) , pointer :: h_west_t0, h_east_t0, h_north_t0, h_south_t0
     real,  dimension(:,:,:) , pointer :: u_west_t0, u_east_t0, u_north_t0, u_south_t0
     real,  dimension(:,:,:) , pointer :: v_west_t0, v_east_t0, v_north_t0, v_south_t0
     !uc and vc need two time levels because they are interpolated to half-timesteps, not full ones
     real,  dimension(:,:,:) , pointer :: uc_west_t0, uc_east_t0, uc_north_t0, uc_south_t0
     real,  dimension(:,:,:) , pointer :: vc_west_t0, vc_east_t0, vc_north_t0, vc_south_t0
     real,  dimension(:,:,:) , pointer :: uc_west_t1, uc_east_t1, uc_north_t1, uc_south_t1
     real,  dimension(:,:,:) , pointer :: vc_west_t1, vc_east_t1, vc_north_t1, vc_south_t1
     real,  dimension(:,:,:,:) , pointer :: q_west, q_east, q_north, q_south
     real,  dimension(:,:,:,:) , pointer :: q_west_t0, q_east_t0, q_north_t0, q_south_t0
#ifndef SW_DYNAMICS
     real,  dimension(:,:,:) , pointer :: pt_west, pt_east, pt_north, pt_south
     real,  dimension(:,:,:) , pointer :: w_west, w_east, w_north, w_south
     real,  dimension(:,:,:) , pointer :: pt_west_t0, pt_east_t0, pt_north_t0, pt_south_t0
     real,  dimension(:,:,:) , pointer :: w_west_t0, w_east_t0, w_north_t0, w_south_t0
#endif
     real, dimension(:,:,:), pointer ::  nest_fx_west_accum, nest_fx_east_accum, nest_fx_south_accum, nest_fx_north_accum
     logical, pointer :: do_flux_BCs, do_2way_flux_BCs !For a parent grid; determine whether there is a need to send BCs

     integer,  dimension(:,:) , pointer :: process_bounds

     !Hold on to coarse-grid global grid, so we don't have to waste processor time getting it again when starting to do grid nesting
     real,  dimension(:,:,:,:) , pointer :: grid_global

!!!!!!!!!!!!!!!!!!!
! From fv_physics !
!!!!!!!!!!!!!!!!!!!
     
   real, pointer, dimension(:,:,:)   :: t_phys
   real, pointer, dimension(:,:,:,:) :: q_phys
   real, pointer, dimension(:,:,:)   :: u_dt, v_dt, t_dt
   real, pointer, dimension(:,:,:,:) :: q_dt  ! potentially a huge array
   real, pointer, dimension(:,:,:)   :: p_full, z_full, p_half, z_half
   integer, pointer :: nx_win, ny_win      ! iew-isw+1, jew-jsw+1 (window sizes)
   integer, pointer :: nx_dom, ny_dom      ! ie-is+1, je-js+1 (compute domain sizes)
   integer, pointer, dimension(:)  :: physics_window_x, physics_window_y !Allocated in fv_physics
   integer, pointer :: ny_per_thread, num_phys_windows

  integer, pointer :: atmos_axes(:)


contains

  subroutine switch_current_grid_pointers(Atm)
    type(fv_atmos_type), intent(IN), target :: Atm

     grid_number                   => Atm%grid_number
     dt_grid                       => Atm%dt_grid
     refinement                    => Atm%refinement
     Time_init                     => Atm%Time_init
     Time                          => Atm%Time
     Run_length                    => Atm%Run_length
     Time_end                      => Atm%Time_end
     Time_step_atmos               => Atm%Time_step_atmos
     grid_active                   => Atm%grid_active
     u                             => Atm%u
     v                             => Atm%v
     pt                            => Atm%pt
     delp                          => Atm%delp
     q                             => Atm%q
     w                             => Atm%w
     delz                          => Atm%delz
     ze0                           => Atm%ze0
     ps                            => Atm%ps
     pe                            => Atm%pe
     pk                            => Atm%pk
     peln                          => Atm%peln
     pkz                           => Atm%pkz
     u_srf                         => Atm%u_srf
     v_srf                         => Atm%v_srf
     sgh                           => Atm%sgh
     oro                           => Atm%oro
     ts                            => Atm%ts
     phis                          => Atm%phis
     omga                          => Atm%omga
     ua                            => Atm%ua
     va                            => Atm%va
     uc                            => Atm%uc
     vc                            => Atm%vc
     ak                            => Atm%ak
     bk                            => Atm%bk
     mfx                           => Atm%mfx
     mfy                           => Atm%mfy
     cx                            => Atm%cx
     cy                            => Atm%cy
#ifndef NO_GRID_G
     grid_g                        => Atm%grid_g
#endif
     isc                           => Atm%isc
     iec                           => Atm%iec
     jsc                           => Atm%jsc
     jec                           => Atm%jec
     res_latlon_dynamics           => Atm%res_latlon_dynamics
     res_latlon_tracers            => Atm%res_latlon_tracers

     grid_type                     => Atm%grid_type
     grid_name                     => Atm%grid_name
     grid_file                     => Atm%grid_file
     hord_mt                       => Atm%hord_mt
     kord_mt                       => Atm%kord_mt
     kord_wz                       => Atm%kord_wz
     hord_vt                       => Atm%hord_vt
     hord_tm                       => Atm%hord_tm
     hord_dp                       => Atm%hord_dp
     kord_tm                       => Atm%kord_tm
     hord_tr                       => Atm%hord_tr
     kord_tr                       => Atm%kord_tr
     scale_z                       => Atm%scale_z
     w_max                         => Atm%w_max
     z_min                         => Atm%z_min
     nord                          => Atm%nord
     dddmp                         => Atm%dddmp
     d2_bg                         => Atm%d2_bg
     d4_bg                         => Atm%d4_bg
     vtdm4                         => Atm%vtdm4
     d2_bg_k1                      => Atm%d2_bg_k1
     d2_bg_k2                      => Atm%d2_bg_k2
     d2_divg_max_k1                => Atm%d2_divg_max_k1
     d2_divg_max_k2                => Atm%d2_divg_max_k2
     damp_k_k1                     => Atm%damp_k_k1
     damp_k_k2                     => Atm%damp_k_k2
     n_zs_filter                   => Atm%n_zs_filter
     nord_zs_filter                => Atm%nord_zs_filter
     no_dycore                     => Atm%no_dycore
     replace_w                     => Atm%replace_w
     convert_ke                    => Atm%convert_ke
     do_vort_damp                  => Atm%do_vort_damp
     use_old_omega                 => Atm%use_old_omega
     beta                          => Atm%beta
     n_sponge                      => Atm%n_sponge
     d_ext                         => Atm%d_ext
     nwat                          => Atm%nwat
     warm_start                    => Atm%warm_start
     inline_q                      => Atm%inline_q
     n_sponge                      => Atm%n_sponge
     d_ext                         => Atm%d_ext
     nwat                          => Atm%nwat
     warm_start                    => Atm%warm_start
     inline_q                      => Atm%inline_q
     shift_fac                     => Atm%shift_fac
     do_schmidt                    => Atm%do_schmidt
     stretch_fac                   => Atm%stretch_fac
     target_lat                    => Atm%target_lat
     target_lon                    => Atm%target_lon
     reset_eta                     => Atm%reset_eta
     p_fac                         => Atm%p_fac
     a_imp                         => Atm%a_imp
     n_split                       => Atm%n_split
     m_split                       => Atm%m_split
     k_split                       => Atm%k_split
     q_split                       => Atm%q_split
     print_freq                    => Atm%print_freq
     npx                           => Atm%npx
     npy                           => Atm%npy
     npz                           => Atm%npz
     npz_rst                       => Atm%npz_rst
     ncnst                         => Atm%ncnst
     pnats                         => Atm%pnats
     ntiles                        => Atm%ntiles
     nf_omega                      => Atm%nf_omega
     fv_sg_adj                     => Atm%fv_sg_adj
     na_init                       => Atm%na_init
     p_ref                         => Atm%p_ref
#ifdef MARS_GCM
     reference_sfc_pres            => Atm%reference_sfc_pres
     sponge_damp                   => Atm%sponge_damp
#endif
     dry_mass                      => Atm%dry_mass
     p_ref                         => Atm%p_ref
     dry_mass                      => Atm%dry_mass
     nt_prog                       => Atm%nt_prog
     nt_phys                       => Atm%nt_phys
     tau_h2o                       => Atm%tau_h2o
     d_con                         => Atm%d_con
     consv_te                      => Atm%consv_te
     tau                           => Atm%tau
     rf_center                     => Atm%rf_center
     tq_filter                     => Atm%tq_filter
     filter_phys                   => Atm%filter_phys
     dwind_2d                      => Atm%dwind_2d
     breed_vortex_inline           => Atm%breed_vortex_inline
     range_warn                    => Atm%range_warn
     fill                          => Atm%fill
     fill_dp                       => Atm%fill_dp
     fill_wz                       => Atm%fill_wz
     non_ortho                     => Atm%non_ortho
     adiabatic                     => Atm%adiabatic
     moist_phys                    => Atm%moist_phys
     do_Held_Suarez                => Atm%do_Held_Suarez
     reproduce_sum                 => Atm%reproduce_sum
     adjust_dry_mass               => Atm%adjust_dry_mass
     fv_debug                      => Atm%fv_debug
     srf_init                      => Atm%srf_init
     mountain                      => Atm%mountain
     remap_t                       => Atm%remap_t
     z_tracer                      => Atm%z_tracer
     old_divg_damp                 => Atm%old_divg_damp
     master                        => Atm%master
     fv_land                       => Atm%fv_land
     nudge                         => Atm%nudge
     ncep_ic                       => Atm%ncep_ic
     fv_diag_ic                    => Atm%fv_diag_ic
     external_ic                   => Atm%external_ic
     res_latlon_dynamics           => Atm%res_latlon_dynamics
     res_latlon_tracers            => Atm%res_latlon_tracers
     hydrostatic                   => Atm%hydrostatic
     phys_hydrostatic              => Atm%phys_hydrostatic
     hybrid_z                      => Atm%hybrid_z
     Make_NH                       => Atm%Make_NH
     make_hybrid_z                 => Atm%make_hybrid_z
     a2b_ord                       => Atm%a2b_ord
     c2l_ord                       => Atm%c2l_ord
     ks                            => Atm%ks
     ndims                         => Atm%ndims

     ng                            => Atm%ng
     domain                        => Atm%domain
     domain_for_coupler            => Atm%domain_for_coupler
     num_contact                   => Atm%num_contact
     npes_per_tile                 => Atm%npes_per_tile
     tile                          => Atm%tile
     npes_this_grid                => Atm%npes_this_grid
     is                            => Atm%is
     ie                            => Atm%ie
     js                            => Atm%js
     je                            => Atm%je
     isd                           => Atm%isd
     ied                           => Atm%ied
     jsd                           => Atm%jsd
     jed                           => Atm%jed
     isc                           => Atm%isc
     iec                           => Atm%iec
     jsc                           => Atm%jsc
     jec                           => Atm%jec
     square_domain                 => Atm%square_domain
     npes_x                        => Atm%npes_x
     npes_y                        => Atm%npes_y
     layout                        => Atm%layout
     io_layout                      => Atm%io_layout
     nest_domain                   => Atm%nest_domain
     this_proc_sends               => Atm%this_proc_sends
     this_proc_recvs               => Atm%this_proc_recvs

     csFac                         => Atm%csFac
     zeta                          => Atm%zeta
     stretch                       => Atm%stretch
     dxdy_area                     => Atm%dxdy_area
     latlon                        => Atm%latlon
     cubed_sphere                  => Atm%cubed_sphere
     double_periodic               => Atm%double_periodic
     latlon_patch                  => Atm%latlon_patch
     latlon_strip                  => Atm%latlon_strip
     channel                       => Atm%channel
     have_south_pole               => Atm%have_south_pole
     have_north_pole               => Atm%have_north_pole
     interpOrder                   => Atm%interpOrder
     debug_message_size            => Atm%debug_message_size
     write_grid_char_file          => Atm%write_grid_char_file
     stretched_grid                => Atm%stretched_grid


     npx_g                         => Atm%npx_g
     npy_g                         => Atm%npy_g
     npz_g                         => Atm%npz_g
     ntiles_g                      => Atm%ntiles_g
     grid                          => Atm%grid
     agrid                         => Atm%agrid
     area                          => Atm%area
     area_c                        => Atm%area_c
     sina                          => Atm%sina
     cosa                          => Atm%cosa
     e1                            => Atm%e1
     e2                            => Atm%e2
     dx                            => Atm%dx
     dy                            => Atm%dy
     dxc                           => Atm%dxc
     dyc                           => Atm%dyc
     dxa                           => Atm%dxa
     dya                           => Atm%dya
     rarea                         => Atm%rarea
     rarea_c                       => Atm%rarea_c
     rdx                           => Atm%rdx
     rdy                           => Atm%rdy
     rdxc                          => Atm%rdxc
     rdyc                          => Atm%rdyc
     rdxa                          => Atm%rdxa
     rdya                          => Atm%rdya
     acapN                         => Atm%acapN
     acapS                         => Atm%acapS
     globalarea                    => Atm%globalarea
     iinta                         => Atm%iinta
     jinta                         => Atm%jinta
     iintb                         => Atm%iintb
     jintb                         => Atm%jintb
     dx_const                      => Atm%dx_const
     dy_const                      => Atm%dy_const
     deglon_start                  => Atm%deglon_start
     deglon_stop                   => Atm%deglon_stop
     deglat_start                  => Atm%deglat_start
     deglat_stop                   => Atm%deglat_stop

     idiag                         => Atm%idiag
     steps                         => Atm%steps
     efx                           => Atm%efx
     efx_nest                      => Atm%efx_nest
     mtq                           => Atm%mtq
     efx_sum                       => Atm%efx_sum
     efx_sum_nest                  => Atm%efx_sum_nest
     mtq_sum                       => Atm%mtq_sum
     edge_s                        => Atm%edge_s
     edge_n                        => Atm%edge_n
     edge_w                        => Atm%edge_w
     edge_e                        => Atm%edge_e
     edge_vect_s                   => Atm%edge_vect_s
     edge_vect_n                   => Atm%edge_vect_n
     edge_vect_w                   => Atm%edge_vect_w
     edge_vect_e                   => Atm%edge_vect_e
     ex_s                          => Atm%ex_s
     ex_n                          => Atm%ex_n
     ex_w                          => Atm%ex_w
     ex_e                          => Atm%ex_e
     divg_u                        => Atm%divg_u
     divg_v                        => Atm%divg_v
     a11                           => Atm%a11
     a12                           => Atm%a12
     a21                           => Atm%a21
     a22                           => Atm%a22
     z11                           => Atm%z11
     z12                           => Atm%z12
     z21                           => Atm%z21
     z22                           => Atm%z22
     global_area                   => Atm%global_area
     stretch_factor                => Atm%stretch_factor
     g_sum_initialized             => Atm%g_sum_initialized
     gnomonic_grid                 => Atm%gnomonic_grid
     sw_corner                     => Atm%sw_corner
     se_corner                     => Atm%se_corner
     ne_corner                     => Atm%ne_corner
     nw_corner                     => Atm%nw_corner
     cosa_u                        => Atm%cosa_u
     cosa_v                        => Atm%cosa_v
     cosa_s                        => Atm%cosa_s
     sina_s                        => Atm%sina_s
     sina_u                        => Atm%sina_u
     sina_v                        => Atm%sina_v
     rsin_u                        => Atm%rsin_u
     rsin_v                        => Atm%rsin_v
     rsina                         => Atm%rsina
     rsin2                         => Atm%rsin2
     ee1                           => Atm%ee1
     ee2                           => Atm%ee2
     ec1                           => Atm%ec1
     ec2                           => Atm%ec2
     ew                            => Atm%ew
     es                            => Atm%es
     sin_sg                        => Atm%sin_sg
     cos_sg                        => Atm%cos_sg
     en1                           => Atm%en1
     en2                           => Atm%en2
     eww                           => Atm%eww
     ess                           => Atm%ess
     vlon                          => Atm%vlon
     vlat                          => Atm%vlat
     fC                            => Atm%fC
     f0                            => Atm%f0
     deglat                        => Atm%deglat
     ptop                          => Atm%ptop
     da_min                        => Atm%da_min
     da_max                        => Atm%da_max
     da_min_c                      => Atm%da_min_c
     da_max_c                      => Atm%da_max_c
     fv_olr                        => Atm%fv_olr
     fv_abs_sw                     => Atm%fv_abs_sw
     irad                          => Atm%irad
     zs_g                          => Atm%zs_g
     sgh_g                         => Atm%sgh_g
     oro_g                         => Atm%oro_g
     test_case                     => Atm%test_case
     alpha                         => Atm%alpha
     Ubar                          => Atm%Ubar
     gh0                           => Atm%gh0
     case9_B                       => Atm%case9_B
     AofT                          => Atm%AofT
     phi0                          => Atm%phi0
     ua0                           => Atm%ua0
     va0                           => Atm%va0
     gh_table                      => Atm%gh_table
     lats_table                    => Atm%lats_table
     gh_initialized                => Atm%gh_initialized
     tmass_orig                    => Atm%tmass_orig
     tvort_orig                    => Atm%tvort_orig
     tener_orig                    => Atm%tener_orig


     nudge_initialized             => Atm%nudge_initialized
     u0                            => Atm%u0
     v0                            => Atm%v0
     t0                            => Atm%t0
     dp                            => Atm%dp
#ifdef MARS_GCM
     tmars_initialized             => Atm%tmars_initialized
     tmars                         => Atm%tmars
#endif
     Fv_restart                    => Atm%Fv_restart
     SST_restart                   => Atm%SST_restart
     Fv_tile_restart               => Atm%Fv_tile_restart
     Rsf_restart                   => Atm%Rsf_restart
     Mg_restart                    => Atm%Mg_restart
     Lnd_restart                   => Atm%Lnd_restart
     Tra_restart                   => Atm%Tra_restart
     parent_grid                   => Atm%parent_grid
     parent_tile                   => Atm%parent_tile
     nest_refinement               => Atm%nest_refinement
     nested                        => Atm%nested
     nestbctype                    => Atm%nestbctype
     nsponge                       => Atm%nsponge
     nestupdate                    => Atm%nestupdate
     twowaynest                    => Atm%twowaynest
     ioffset                       => Atm%ioffset
     joffset                       => Atm%joffset
     nest_timestep                 => Atm%nest_timestep
     tracer_nest_timestep          => Atm%tracer_nest_timestep
     s_weight                      => Atm%s_weight
     first_step                    => Atm%first_step
     refinement_of_global          => Atm%refinement_of_global
     npx_global                    => Atm%npx_global
     child_grids                   => Atm%child_grids
     pelist                        => Atm%pelist
     parent_proc                   => Atm%parent_proc
     child_proc                    => Atm%child_proc
     ind_h                         => Atm%ind_h
     ind_u                         => Atm%ind_u
     ind_v                         => Atm%ind_v
     wt_h                          => Atm%wt_h
     wt_u                          => Atm%wt_u
     wt_v                          => Atm%wt_v
     ind_update_h                  => Atm%ind_update_h
     h_west                        => Atm%h_west
     h_east                        => Atm%h_east
     h_north                       => Atm%h_north
     h_south                       => Atm%h_south
     u_west                        => Atm%u_west
     u_east                        => Atm%u_east
     u_north                       => Atm%u_north
     u_south                       => Atm%u_south
     v_west                        => Atm%v_west
     v_east                        => Atm%v_east
     v_north                       => Atm%v_north
     v_south                       => Atm%v_south
     h_west_t0                        => Atm%h_west_t0
     h_east_t0                        => Atm%h_east_t0
     h_north_t0                       => Atm%h_north_t0
     h_south_t0                       => Atm%h_south_t0
     u_west_t0                        => Atm%u_west_t0
     u_east_t0                        => Atm%u_east_t0
     u_north_t0                       => Atm%u_north_t0
     u_south_t0                       => Atm%u_south_t0
     v_west_t0                        => Atm%v_west_t0
     v_east_t0                        => Atm%v_east_t0
     v_north_t0                       => Atm%v_north_t0
     v_south_t0                       => Atm%v_south_t0

     uc_west_t0                    => Atm%uc_west_t0
     uc_east_t0                    => Atm%uc_east_t0
     uc_north_t0                   => Atm%uc_north_t0
     uc_south_t0                   => Atm%uc_south_t0
     vc_west_t0                    => Atm%vc_west_t0
     vc_east_t0                    => Atm%vc_east_t0
     vc_north_t0                   => Atm%vc_north_t0
     vc_south_t0                   => Atm%vc_south_t0
     uc_west_t1                    => Atm%uc_west_t1
     uc_east_t1                    => Atm%uc_east_t1
     uc_north_t1                   => Atm%uc_north_t1
     uc_south_t1                   => Atm%uc_south_t1
     vc_west_t1                    => Atm%vc_west_t1
     vc_east_t1                    => Atm%vc_east_t1
     vc_north_t1                   => Atm%vc_north_t1
     vc_south_t1                   => Atm%vc_south_t1
     q_west                        => Atm%q_west
     q_east                        => Atm%q_east
     q_north                       => Atm%q_north
     q_south                       => Atm%q_south
     q_west_t0                        => Atm%q_west_t0
     q_east_t0                        => Atm%q_east_t0
     q_north_t0                       => Atm%q_north_t0
     q_south_t0                       => Atm%q_south_t0
#ifndef SW_DYNAMICS
     pt_west                       => Atm%pt_west
     pt_east                       => Atm%pt_east
     pt_north                      => Atm%pt_north
     pt_south                      => Atm%pt_south
     w_west                        => Atm%w_west
     w_east                        => Atm%w_east
     w_north                       => Atm%w_north
     w_south                       => Atm%w_south
     pt_west_t0                       => Atm%pt_west_t0
     pt_east_t0                       => Atm%pt_east_t0
     pt_north_t0                      => Atm%pt_north_t0
     pt_south_t0                      => Atm%pt_south_t0
     w_west_t0                        => Atm%w_west_t0
     w_east_t0                        => Atm%w_east_t0
     w_north_t0                       => Atm%w_north_t0
     w_south_t0                       => Atm%w_south_t0
#endif

     if (Atm%nestbctype > 1 .and. ncnst > 0) then
        nest_fx_west_accum              => Atm%nest_fx_west_accum
        nest_fx_east_accum              => Atm%nest_fx_east_accum
        nest_fx_south_accum             => Atm%nest_fx_south_accum
        nest_fx_north_accum             => Atm%nest_fx_north_accum
     endif
     do_flux_BCs                   => Atm%do_flux_BCs
     do_2way_flux_BCs              => Atm%do_2way_flux_BCs

     process_bounds                => Atm%process_bounds
     grid_global                   => Atm%grid_global
     t_phys                        => Atm%t_phys
     q_phys                        => Atm%q_phys
     u_dt                          => Atm%u_dt
     v_dt                          => Atm%v_dt
     t_dt                          => Atm%t_dt
     q_dt                          => Atm%q_dt
     p_full                         => Atm%p_full
     z_full                         => Atm%z_full
     p_half                         => Atm%p_half
     z_half                         => Atm%z_half
     nx_win                         => Atm%nx_win          
     ny_win                         => Atm%ny_win          
     nx_dom                         => Atm%nx_dom          
     ny_dom                         => Atm%ny_dom          
     physics_window_x               => Atm%physics_window_x
     physics_window_y               => Atm%physics_window_y
     ny_per_thread                  => Atm%ny_per_thread   
     num_phys_windows               => Atm%num_phys_windows
     atmos_axes                     => Atm%atmos_axes
   end subroutine switch_current_grid_pointers

 end module fv_current_grid_mod
