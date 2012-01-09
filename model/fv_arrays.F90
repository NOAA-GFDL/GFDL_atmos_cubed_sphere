module fv_arrays_mod
#include <fms_platform.h>
 use mpp_domains_mod,  only: domain2d
public
  type fv_atmos_type
     type(domain2d), pointer :: domain =>NULL()
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

! Horizontal Grid descriptors
    real, pointer :: grid(:,:,:)  _NULL  ! Leave as a pointer for now
    real, pointer :: agrid(:,:,:)  _NULL  ! Leave as a pointer for now
    real, pointer :: grid_g(:,:,:) _NULL  ! "global" grid (one face of a cube)

    real :: consv_te
    real :: shift_fac, stretch_fac, target_lat, target_lon
    logical :: do_schmidt

    integer :: isc, iec, jsc, jec
    integer :: isd, ied, jsd, jed
    integer :: ks, npx, npy, npz, npz_rst, ng, ntiles
    integer :: n_sponge    ! Number of sponge layers at the top of the atmosphere
    integer :: k_top       ! Starting layer for non-hydrostatic dynamics
    integer :: ncnst, pnats, ndims, k_split, n_split, m_split, q_split, print_freq
    integer :: nwat        ! water substance
    integer :: fv_sg_adj
    integer :: na_init     ! number of iteraation for the adiabatic initialization
    integer :: n_zs_filter, nord_zs_filter

! Namelist control values
    logical :: fill
    logical :: range_warn
    logical :: z_tracer
    logical :: do_Held_Suarez
    logical :: reproduce_sum
    logical :: moist_phys
    logical :: srf_init
    logical :: mountain
    logical :: non_ortho
    logical :: adjust_dry_mass
    logical :: hydrostatic, phys_hydrostatic
    logical :: hybrid_z, Make_NH, make_hybrid_z
    logical :: external_ic
    logical :: ncep_ic
    logical :: fv_diag_ic
    logical :: fv_land
    logical :: nudge
    logical :: tq_filter
    logical :: warm_start

    character(len=128) :: res_latlon_dynamics  ! restart file from the latlon FV core
    character(len=128) :: res_latlon_tracers   ! tracer restart file from the latlon core

    real    :: dry_mass

  end type fv_atmos_type

!---- version number -----
  character(len=128) :: version = '$Id: fv_arrays.F90,v 19.0 2012/01/06 19:57:34 fms Exp $'
  character(len=128) :: tagname = '$Name: siena $'

end module fv_arrays_mod
