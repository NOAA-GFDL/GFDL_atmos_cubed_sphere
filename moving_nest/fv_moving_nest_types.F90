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

!***********************************************************************
!> @file
!! @brief Provides data structures for moving nest functionality
!! @author W. Ramstrom, AOML/HRD   03/24/2022
!! @email William.Ramstrom@noaa.gov
! =======================================================================!

module fv_moving_nest_types_mod

#ifdef MOVING_NEST
#include <fms_platform.h>

#ifdef GFS_TYPES
  use GFS_typedefs,           only: kind_phys
#else
  use IPD_typedefs,           only: kind_phys => IPD_kind_phys
#endif

  use fms_mod,         only: check_nml_error
  use fv_arrays_mod,   only: fv_atmos_type
  use fv_mp_mod,       only: MAX_NNEST
  use mpp_mod,         only: input_nml_file, mpp_pe

  implicit none

  type fv_moving_nest_flag_type
    ! Moving Nest Namelist Variables
    logical               :: is_moving_nest = .false.
    character(len=120)    :: surface_dir = "INPUT/moving_nest"
    integer               :: terrain_smoother = 1
    integer               :: vortex_tracker = 0
    integer               :: ntrack = 1
    integer               :: corral_x = 5
    integer               :: corral_y = 5

    integer               :: outatcf_lun = 600

    ! Moving nest related variables
    integer               :: move_cd_x = 0
    integer               :: move_cd_y = 0
    logical               :: do_move = .false.
  end type fv_moving_nest_flag_type

  ! Encapsulates the grid definition data, such as read from the netCDF files
  type grid_geometry
    integer   :: nx, ny, nxp, nyp

    real(kind=kind_phys), allocatable  :: lats(:,:)
    real(kind=kind_phys), allocatable  :: lons(:,:)

    !real, allocatable  :: dx(:,:)
    !real, allocatable  :: dy(:,:)
    real(kind=kind_phys), allocatable  :: area(:,:)
  end type grid_geometry

  type fv_moving_nest_prog_type
    real, _ALLOCATABLE                  :: delz(:,:,:)      _NULL   !< layer thickness (meters)
  end type fv_moving_nest_prog_type

  ! TODO deallocate these at end of model run.  They are only allocated once, at first nest move, inside mn_static_read_hires().
  !  Note these are only 32 bits for now; matching the precision of the input netCDF files
  !  though the model generally handles physics variables with 64 bit precision
  type mn_surface_grids
    real, allocatable  :: orog_grid(:,:)               _NULL  ! orography -- raw or filtered depending on namelist option, in meters
    real, allocatable  :: orog_std_grid(:,:)           _NULL  ! terrain standard deviation for gravity wave drag, in meters (?)
    real, allocatable  :: ls_mask_grid(:,:)            _NULL  ! land sea mask -- 0 for ocean/lakes, 1, for land.  Perhaps 2 for sea ice.
    real, allocatable  :: land_frac_grid(:,:)          _NULL  ! Continuous land fraction - 0.0 ocean, 0.5 half of each, 1.0 all land

    real, allocatable  :: parent_orog_grid(:,:)        _NULL  ! parent orography -- only used for terrain_smoother=1.
    !     raw or filtered depending on namelist option,in meters

    ! Soil variables
    real, allocatable  :: deep_soil_temp_grid(:,:)     _NULL  ! deep soil temperature at 5m, in degrees K
    real, allocatable  :: soil_type_grid(:,:)          _NULL  ! STATSGO soil type

    ! Vegetation variables
    real, allocatable  :: veg_frac_grid(:,:)           _NULL  ! vegetation fraction
    real, allocatable  :: veg_type_grid(:,:)           _NULL  ! IGBP vegetation type
    real, allocatable  :: veg_greenness_grid(:,:)      _NULL  ! NESDIS vegetation greenness; netCDF file has monthly values

    ! Orography variables
    real, allocatable  :: slope_type_grid(:,:)         _NULL  ! legacy 1 degree GFS slope type

    ! Albedo variables
    real, allocatable  :: max_snow_alb_grid(:,:)       _NULL  ! max snow albedo
    real, allocatable  :: facsf_grid(:,:)              _NULL  ! fractional coverage with strong cosz dependency
    real, allocatable  :: facwf_grid(:,:)              _NULL  ! fractional coverage with weak cosz dependency

    ! Snow free albedo
    !   strong cosz angle dependence = black sky
    !   weak cosz angle dependence = white sky
    !  From the chgres code in static_data.F90, we see the linkage of variable names:
    !   type(esmf_field), public           :: alvsf_target_grid !< visible black sky albedo
    !   type(esmf_field), public           :: alvwf_target_grid !< visible white sky albedo
    !   type(esmf_field), public           :: alnsf_target_grid !< near ir black sky albedo
    !   type(esmf_field), public           :: alnwf_target_grid !< near ir white sky albedo

    real, allocatable  :: alvsf_grid(:,:)              _NULL  ! Visible black sky albedo; netCDF file has monthly values
    real, allocatable  :: alvwf_grid(:,:)              _NULL  ! Visible white sky albedo; netCDF file has monthly values
    real, allocatable  :: alnsf_grid(:,:)              _NULL  ! Near IR black sky albedo; netCDF file has monthly values
    real, allocatable  :: alnwf_grid(:,:)              _NULL  ! Near IR white sky albedo; netCDF file has monthly values

  end type mn_surface_grids

  type fv_moving_nest_physics_type
    real, _ALLOCATABLE                  :: ts(:,:)          _NULL   !< 2D skin temperature/SST
    real, _ALLOCATABLE                  :: slmsk(:,:)       _NULL   !< land sea mask -- 0 for ocean/lakes, 1, for land.  Perhaps 2 for sea ice.
    real (kind=kind_phys), _ALLOCATABLE :: smc (:,:,:)      _NULL   !< soil moisture content
    real (kind=kind_phys), _ALLOCATABLE :: stc (:,:,:)      _NULL   !< soil temperature
    real (kind=kind_phys), _ALLOCATABLE :: slc (:,:,:)      _NULL   !< soil liquid water content

    real (kind=kind_phys), _ALLOCATABLE :: u10m (:,:)       _NULL   !< 10m u wind (a-grid?)
    real (kind=kind_phys), _ALLOCATABLE :: v10m (:,:)       _NULL   !< 10m v wind (a-grid?)
    real (kind=kind_phys), _ALLOCATABLE :: hprime (:,:,:)   _NULL   !< orographic metrics (maybe standard deviation?)

    real (kind=kind_phys), _ALLOCATABLE :: tprcp (:,:)      _NULL   !< total (of all precip types) precipitation rate

    real (kind=kind_phys), _ALLOCATABLE :: zorl (:,:)       _NULL   !< roughness length
    real (kind=kind_phys), _ALLOCATABLE :: zorll (:,:)      _NULL   !< land roughness length
    !real (kind=kind_phys), _ALLOCATABLE :: zorli (:,:)     _NULL   !< ice surface roughness length ! TODO do we need this?
    real (kind=kind_phys), _ALLOCATABLE :: zorlw (:,:)      _NULL   !< wave surface roughness length
    real (kind=kind_phys), _ALLOCATABLE :: zorlwav (:,:)    _NULL   !< wave surface roughness in cm derived from wave model

    real (kind=kind_phys), _ALLOCATABLE :: sfalb_lnd(:,:)   _NULL   !< surface albedo over land for LSM
    real (kind=kind_phys), _ALLOCATABLE :: emis_lnd(:,:)    _NULL   !< surface emissivity over land for LSM
    real (kind=kind_phys), _ALLOCATABLE :: emis_ice(:,:)    _NULL   !< surface emissivity over ice for LSM
    real (kind=kind_phys), _ALLOCATABLE :: emis_wat(:,:)    _NULL   !< surface emissivity over water for LSM
    real (kind=kind_phys), _ALLOCATABLE :: sfalb_lnd_bck(:,:) _NULL !< snow-free albedo over land

    !real (kind=kind_phys), _ALLOCATABLE :: semis(:,:)       _NULL   !< surface lw emissivity in fraction 
    !real (kind=kind_phys), _ALLOCATABLE :: semisbase(:,:)   _NULL   !< background surface emissivity 
    !real (kind=kind_phys), _ALLOCATABLE :: sfalb(:,:)       _NULL   !< mean surface diffused sw albedo 

    real (kind=kind_phys), _ALLOCATABLE :: alvsf(:,:)       _NULL   !< visible black sky albedo
    real (kind=kind_phys), _ALLOCATABLE :: alvwf(:,:)       _NULL   !< visible white sky albedo
    real (kind=kind_phys), _ALLOCATABLE :: alnsf(:,:)       _NULL   !< near IR black sky albedo
    real (kind=kind_phys), _ALLOCATABLE :: alnwf(:,:)       _NULL   !< near IR white sky albedo

    real (kind=kind_phys), _ALLOCATABLE :: albdirvis_lnd(:,:)       _NULL   !<
    real (kind=kind_phys), _ALLOCATABLE :: albdirnir_lnd(:,:)       _NULL   !<
    real (kind=kind_phys), _ALLOCATABLE :: albdifvis_lnd(:,:)       _NULL   !<
    real (kind=kind_phys), _ALLOCATABLE :: albdifnir_lnd(:,:)       _NULL   !<

    real (kind=kind_phys), _ALLOCATABLE :: facsf(:,:)       _NULL   !< fractional coverage for strong zenith angle albedo
    real (kind=kind_phys), _ALLOCATABLE :: facwf(:,:)       _NULL   !< fractional coverage for strong zenith angle albedo

    real (kind=kind_phys), _ALLOCATABLE :: lakefrac (:,:)   _NULL   !< lake  fraction [0:1] 
    real (kind=kind_phys), _ALLOCATABLE :: lakedepth (:,:)  _NULL   !< lake  depth [ m ]    

    real (kind=kind_phys), _ALLOCATABLE :: canopy (:,:)     _NULL   !< canopy water content
    real (kind=kind_phys), _ALLOCATABLE :: vegfrac (:,:)    _NULL   !< vegetation fraction
    real (kind=kind_phys), _ALLOCATABLE :: uustar (:,:)     _NULL   !< u* wind in similarity theory
    real (kind=kind_phys), _ALLOCATABLE :: shdmin (:,:)     _NULL   !< min fractional coverage of green vegetation
    real (kind=kind_phys), _ALLOCATABLE :: shdmax (:,:)     _NULL   !< max fractional coverage of green vegetation
    real (kind=kind_phys), _ALLOCATABLE :: tsfco (:,:)      _NULL   !< surface temperature ocean
    real (kind=kind_phys), _ALLOCATABLE :: tsfcl (:,:)      _NULL   !< surface temperature land
    real (kind=kind_phys), _ALLOCATABLE :: tsfc (:,:)       _NULL   !< surface temperature
    !real (kind=kind_phys), _ALLOCATABLE :: tsfc_radtime (:,:) _NULL !< surface temperature on radiative timestep

    real (kind=kind_phys), _ALLOCATABLE :: cv  (:,:)        _NULL   !< fraction of convective cloud
    real (kind=kind_phys), _ALLOCATABLE :: cvt (:,:)        _NULL   !< convective cloud top pressure
    real (kind=kind_phys), _ALLOCATABLE :: cvb (:,:)        _NULL   !< convective cloud bottom pressure

    real (kind=kind_phys), _ALLOCATABLE :: phy_f2d (:,:,:)  _NULL   !< 2D physics variables
    real (kind=kind_phys), _ALLOCATABLE :: phy_f3d(:,:,:,:) _NULL   !< 3D physics variables

    ! NSST Variables

    real (kind=kind_phys), _ALLOCATABLE :: tref (:,:)       _NULL   !< reference temperature for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: z_c (:,:)        _NULL   !< coefficient for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: c_0 (:,:)        _NULL   !< coefficient for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: c_d (:,:)        _NULL   !< coefficient for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: w_0 (:,:)        _NULL   !< coefficient for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: w_d (:,:)        _NULL   !< coefficient for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: xt (:,:)         _NULL   !< heat content for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: xs (:,:)         _NULL   !< salinity for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: xu (:,:)         _NULL   !< u current constant for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: xv (:,:)         _NULL   !< v current constant for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: xz (:,:)         _NULL   !< DTL thickness for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: zm (:,:)         _NULL   !< MXL for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: xtts (:,:)       _NULL   !< d(xt)/d(ts) for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: xzts (:,:)       _NULL   !< d(xz)/d(ts) for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: d_conv (:,:)     _NULL   !< think of free convection layer for NSSTM
    ! real (kind=kind_phys), _ALLOCATABLE :: ifd (:,:)      _NULL   !< index to start DTM run  for NSSTM   ! TODO Probably can't interpolate an index.
    !  IFD values are 0 for land, and 1 for oceans/lakes -- reverse of the land sea mask
    !  Land Sea Mask has values of 0 for oceans/lakes, 1 for land, 2 for sea ice
    real (kind=kind_phys), _ALLOCATABLE :: dt_cool (:,:)    _NULL   !< sub-layer cooling amount for NSSTM
    real (kind=kind_phys), _ALLOCATABLE :: qrain (:,:)      _NULL   !< sensible heat flux due to rainfall for NSSTM

  end type fv_moving_nest_physics_type

  type fv_moving_nest_type
    type(fv_moving_nest_flag_type)    :: mn_flag   ! Mostly namelist variables    
    type(mn_surface_grids)            :: mn_static
    type(fv_moving_nest_prog_type)    :: mn_prog
    type(fv_moving_nest_physics_type) :: mn_phys

    type(grid_geometry)               :: parent_geo
    type(grid_geometry)               :: fp_super_tile_geo
  end type fv_moving_nest_type

  ! Moving Nest Namelist Variables
  logical, dimension(MAX_NNEST) :: is_moving_nest = .False.
  character(len=120)            :: surface_dir = "INPUT/moving_nest"
  integer, dimension(MAX_NNEST) :: terrain_smoother = 1  ! 0 -- all high-resolution data, 1 - static nest smoothing algorithm with blending zone of 5 points, 2 - blending zone of 10 points, 5 - 5 point smoother, 9 - 9 point smoother
  integer, dimension(MAX_NNEST) :: vortex_tracker = 0 ! 0 - not a moving nest, tracker not needed
  ! 1 - prescribed nest moving
  ! 2 - following child domain center
  ! 3 - tracking Min MSLP
  ! 6 - simplified version of GFDL tracker, adopted from HWRF's internal vortex tracker.
  ! 7 - nearly the full storm tracking algorithm from GFDL vortex tracker. The only part that is missing is the part that gives up when the storm dissipates, which is left out intentionally. Adopted from HWRF's internal vortex tracker.
  integer, dimension(MAX_NNEST) :: ntrack = 1 ! number of dt_atmos steps to call the vortex tracker, tracker time step = ntrack*dt_atmos
  integer, dimension(MAX_NNEST) :: move_cd_x = 0 ! the number of parent domain grid cells to move in i direction
  integer, dimension(MAX_NNEST) :: move_cd_y = 0 ! the number of parent domain grid cells to move in j direction
  ! used to control prescribed nest moving, when vortex_tracker=1
  ! the move happens every ntrack*dt_atmos seconds
  ! positive is to move in increasing i and j direction, and
  ! negative is to move in decreasing i and j direction.
  ! 0 means no move. The limitation is to move only 1 grid cell at each move.
  integer, dimension(MAX_NNEST) :: corral_x = 5 ! Minimum parent gridpoints on each side of nest in i direction
  integer, dimension(MAX_NNEST) :: corral_y = 5 ! Minimum parent gridpoints on each side of nest in j direction
  
  integer, dimension(MAX_NNEST) :: outatcf_lun = 600  ! base fortran unit number to write out the partial atcfunix file from the internal tracker

  type(fv_moving_nest_type), _ALLOCATABLE, target    :: Moving_nest(:)

contains

  subroutine fv_moving_nest_init(Atm)
    type(fv_atmos_type), allocatable, intent(in) :: Atm(:)

    integer :: n, ngrids

    ! Allocate the array of fv_moving_nest_type structures of the proper length
    allocate(Moving_nest(size(Atm)))

    ! Configure namelist variables

    ngrids = size(Atm)

    ! Read in namelist

    call read_namelist_moving_nest_nml

    do n=1,ngrids
      if (Atm(n)%neststruct%nested) then
        Moving_nest(n)%mn_flag%is_moving_nest         = is_moving_nest(n)
        Moving_nest(n)%mn_flag%surface_dir            = trim(surface_dir)
        Moving_nest(n)%mn_flag%terrain_smoother       = terrain_smoother(n)
        Moving_nest(n)%mn_flag%vortex_tracker         = vortex_tracker(n)
        Moving_nest(n)%mn_flag%ntrack                 = ntrack(n)
        Moving_nest(n)%mn_flag%move_cd_x              = move_cd_x(n)
        Moving_nest(n)%mn_flag%move_cd_y              = move_cd_y(n)
        Moving_nest(n)%mn_flag%corral_x               = corral_x(n)
        Moving_nest(n)%mn_flag%corral_y               = corral_y(n)
        Moving_nest(n)%mn_flag%outatcf_lun            = outatcf_lun(n)
      else
        Moving_nest(n)%mn_flag%is_moving_nest         = .false.
        Moving_nest(n)%mn_flag%vortex_tracker         = 0
        Moving_nest(n)%mn_flag%ntrack                 = 1
        Moving_nest(n)%mn_flag%move_cd_x              = 0
        Moving_nest(n)%mn_flag%move_cd_y              = 0
        Moving_nest(n)%mn_flag%corral_x               = 5
        Moving_nest(n)%mn_flag%corral_y               = 5
        Moving_nest(n)%mn_flag%outatcf_lun            = 600
      endif
    enddo
  end subroutine fv_moving_nest_init

  subroutine read_namelist_moving_nest_nml
    integer :: f_unit, ios, ierr
    namelist /fv_moving_nest_nml/ surface_dir, is_moving_nest, terrain_smoother, &
        vortex_tracker, ntrack, move_cd_x, move_cd_y, corral_x, corral_y, outatcf_lun

#ifdef INTERNAL_FILE_NML
    read (input_nml_file,fv_moving_nest_nml,iostat=ios)
    ierr = check_nml_error(ios,'fv_moving_nest_nml')
#else
    f_unit=open_namelist_file()
    rewind (f_unit)
    read (f_unit,fv_moving_nest_nml,iostat=ios)
    ierr = check_nml_error(ios,'fv_moving_nest_nml')
    call close_file(f_unit)
#endif

  end subroutine read_namelist_moving_nest_nml

  subroutine deallocate_fv_moving_nests(n)
    integer, intent(in)   :: n

    integer :: i

    do i=1,n
      call deallocate_fv_moving_nest(i)
    enddo
    deallocate(Moving_nest)
  end subroutine deallocate_fv_moving_nests

  subroutine deallocate_fv_moving_nest(n)
    integer, intent(in)   :: n

    call deallocate_fv_moving_nest_prog_type(Moving_nest(n)%mn_prog)
    call deallocate_fv_moving_nest_physics_type(Moving_nest(n)%mn_phys)

  end subroutine deallocate_fv_moving_nest


  subroutine  allocate_fv_moving_nest_prog_type(isd, ied, jsd, jed, npz, mn_prog)
    integer, intent(in)                           :: isd, ied, jsd, jed, npz
    type(fv_moving_nest_prog_type), intent(inout) :: mn_prog

    allocate ( mn_prog%delz(isd:ied, jsd:jed, 1:npz) )
    mn_prog%delz = +99999.9

  end subroutine allocate_fv_moving_nest_prog_type

  subroutine  deallocate_fv_moving_nest_prog_type(mn_prog)
    type(fv_moving_nest_prog_type), intent(inout) :: mn_prog

    if (allocated(mn_prog%delz)) deallocate(mn_prog%delz)

  end subroutine deallocate_fv_moving_nest_prog_type

  subroutine  allocate_fv_moving_nest_physics_type(isd, ied, jsd, jed, npz, move_physics, move_nsst, lsoil, nmtvr, levs, ntot2d, ntot3d, mn_phys)
    integer, intent(in)                           :: isd, ied, jsd, jed, npz
    logical, intent(in)                           :: move_physics, move_nsst
    integer, intent(in)                           :: lsoil, nmtvr, levs, ntot2d, ntot3d    ! From IPD_Control
    type(fv_moving_nest_physics_type), intent(inout) :: mn_phys

    ! The local/temporary variables need to be allocated to the larger data (compute + halos) domain so that the nest motion code has halos to use
    allocate ( mn_phys%ts(isd:ied, jsd:jed) )

    if (move_physics) then
      allocate ( mn_phys%slmsk(isd:ied, jsd:jed) )
      allocate ( mn_phys%smc(isd:ied, jsd:jed, lsoil) )
      allocate ( mn_phys%stc(isd:ied, jsd:jed, lsoil) )
      allocate ( mn_phys%slc(isd:ied, jsd:jed, lsoil) )

      allocate ( mn_phys%sfalb_lnd(isd:ied, jsd:jed) )
      allocate ( mn_phys%emis_lnd(isd:ied, jsd:jed) )
      allocate ( mn_phys%emis_ice(isd:ied, jsd:jed) )
      allocate ( mn_phys%emis_wat(isd:ied, jsd:jed) )
      allocate ( mn_phys%sfalb_lnd_bck(isd:ied, jsd:jed) )

      !allocate ( mn_phys%semis(isd:ied, jsd:jed) )
      !allocate ( mn_phys%semisbase(isd:ied, jsd:jed) )
      !allocate ( mn_phys%sfalb(isd:ied, jsd:jed) )

      allocate ( mn_phys%u10m(isd:ied, jsd:jed) )
      allocate ( mn_phys%v10m(isd:ied, jsd:jed) )
      allocate ( mn_phys%tprcp(isd:ied, jsd:jed) )

      allocate ( mn_phys%hprime(isd:ied, jsd:jed, nmtvr) )

      allocate ( mn_phys%zorl(isd:ied, jsd:jed) )
      allocate ( mn_phys%zorll(isd:ied, jsd:jed) )
      allocate ( mn_phys%zorlwav(isd:ied, jsd:jed) )
      allocate ( mn_phys%zorlw(isd:ied, jsd:jed) )

      allocate ( mn_phys%alvsf(isd:ied, jsd:jed) )
      allocate ( mn_phys%alvwf(isd:ied, jsd:jed) )
      allocate ( mn_phys%alnsf(isd:ied, jsd:jed) )
      allocate ( mn_phys%alnwf(isd:ied, jsd:jed) )

      allocate ( mn_phys%facsf(isd:ied, jsd:jed) )
      allocate ( mn_phys%facwf(isd:ied, jsd:jed) )

      allocate ( mn_phys%lakefrac(isd:ied, jsd:jed) )
      allocate ( mn_phys%lakedepth(isd:ied, jsd:jed) )

      allocate ( mn_phys%canopy(isd:ied, jsd:jed) )
      allocate ( mn_phys%vegfrac(isd:ied, jsd:jed) )
      allocate ( mn_phys%uustar(isd:ied, jsd:jed) )
      allocate ( mn_phys%shdmin(isd:ied, jsd:jed) )
      allocate ( mn_phys%shdmax(isd:ied, jsd:jed) )
      allocate ( mn_phys%tsfco(isd:ied, jsd:jed) )
      allocate ( mn_phys%tsfcl(isd:ied, jsd:jed) )
      allocate ( mn_phys%tsfc(isd:ied, jsd:jed) )
      !allocate ( mn_phys%tsfc_radtime(isd:ied, jsd:jed) )


      allocate ( mn_phys%albdirvis_lnd (isd:ied, jsd:jed) )
      allocate ( mn_phys%albdirnir_lnd (isd:ied, jsd:jed) )
      allocate ( mn_phys%albdifvis_lnd (isd:ied, jsd:jed) )
      allocate ( mn_phys%albdifnir_lnd (isd:ied, jsd:jed) )

      allocate ( mn_phys%cv(isd:ied, jsd:jed) )
      allocate ( mn_phys%cvt(isd:ied, jsd:jed) )
      allocate ( mn_phys%cvb(isd:ied, jsd:jed) )

      allocate ( mn_phys%phy_f2d(isd:ied, jsd:jed, ntot2d) )
      allocate ( mn_phys%phy_f3d(isd:ied, jsd:jed, levs, ntot3d) )
    end if

    if (move_nsst) then
      allocate ( mn_phys%tref(isd:ied, jsd:jed) )
      allocate ( mn_phys%z_c(isd:ied, jsd:jed) )
      allocate ( mn_phys%c_0(isd:ied, jsd:jed) )
      allocate ( mn_phys%c_d(isd:ied, jsd:jed) )
      allocate ( mn_phys%w_0(isd:ied, jsd:jed) )
      allocate ( mn_phys%w_d(isd:ied, jsd:jed) )
      allocate ( mn_phys%xt(isd:ied, jsd:jed) )
      allocate ( mn_phys%xs(isd:ied, jsd:jed) )
      allocate ( mn_phys%xu(isd:ied, jsd:jed) )
      allocate ( mn_phys%xv(isd:ied, jsd:jed) )
      allocate ( mn_phys%xz(isd:ied, jsd:jed) )
      allocate ( mn_phys%zm(isd:ied, jsd:jed) )
      allocate ( mn_phys%xtts(isd:ied, jsd:jed) )
      allocate ( mn_phys%xzts(isd:ied, jsd:jed) )
      allocate ( mn_phys%d_conv(isd:ied, jsd:jed) )
      !allocate ( mn_phys%ifd(isd:ied, jsd:jed) )
      allocate ( mn_phys%dt_cool(isd:ied, jsd:jed) )
      allocate ( mn_phys%qrain(isd:ied, jsd:jed) )
    end if

    mn_phys%ts = +99999.9
    if (move_physics) then
      mn_phys%slmsk = +99999.9
      mn_phys%smc = +99999.9
      mn_phys%stc = +99999.9
      mn_phys%slc = +99999.9


      mn_phys%sfalb_lnd = +99999.9
      mn_phys%emis_lnd = +99999.9
      mn_phys%emis_ice = +99999.9
      mn_phys%emis_wat = +99999.9
      mn_phys%sfalb_lnd_bck = +99999.9

      !mn_phys%semis = +99999.9
      !mn_phys%semisbase = +99999.9
      !mn_phys%sfalb = +99999.9

      mn_phys%u10m = +99999.9
      mn_phys%v10m = +99999.9
      mn_phys%tprcp = +99999.9

      mn_phys%hprime = +99999.9

      mn_phys%zorl = +99999.9
      mn_phys%zorll = +99999.9
      mn_phys%zorlwav = +99999.9
      mn_phys%zorlw = +99999.9

      mn_phys%alvsf = +99999.9
      mn_phys%alvwf = +99999.9
      mn_phys%alnsf = +99999.9
      mn_phys%alnwf = +99999.9

      mn_phys%facsf = +99999.9
      mn_phys%facwf = +99999.9

      mn_phys%lakefrac = +99999.9
      mn_phys%lakedepth = +99999.9

      mn_phys%canopy = +99999.9
      mn_phys%vegfrac = +99999.9
      mn_phys%uustar = +99999.9
      mn_phys%shdmin = +99999.9
      mn_phys%shdmax = +99999.9
      mn_phys%tsfco = +99999.9
      mn_phys%tsfcl = +99999.9
      mn_phys%tsfc = +99999.9
      !mn_phys%tsfc_radtime = +99999.9

      mn_phys%albdirvis_lnd = +99999.9
      mn_phys%albdirnir_lnd = +99999.9
      mn_phys%albdifvis_lnd = +99999.9
      mn_phys%albdifnir_lnd = +99999.9

      mn_phys%cv = +99999.9
      mn_phys%cvt = +99999.9
      mn_phys%cvb = +99999.9

      mn_phys%phy_f2d = +99999.9
      mn_phys%phy_f3d = +99999.9
    end if

    if (move_nsst) then
      mn_phys%tref = +99999.9
      mn_phys%z_c = +99999.9
      mn_phys%c_0 = +99999.9
      mn_phys%c_d = +99999.9
      mn_phys%w_0 = +99999.9
      mn_phys%w_d = +99999.9
      mn_phys%xt = +99999.9
      mn_phys%xs = +99999.9
      mn_phys%xu = +99999.9
      mn_phys%xv = +99999.9
      mn_phys%xz = +99999.9
      mn_phys%zm = +99999.9
      mn_phys%xtts = +99999.9
      mn_phys%xzts = +99999.9
      mn_phys%d_conv = +99999.9
      !mn_phys%ifd = +99999.9
      mn_phys%dt_cool = +99999.9
      mn_phys%qrain = +99999.9
    end if

  end subroutine allocate_fv_moving_nest_physics_type


  subroutine  deallocate_fv_moving_nest_physics_type(mn_phys)
    type(fv_moving_nest_physics_type), intent(inout) :: mn_phys

    if (allocated(mn_phys%ts)) then
      deallocate ( mn_phys%ts )
    else
      ! If ts was not allocated, then none of this structure was allocated.
      return
    end if

    !  if move_phys
    if (allocated(mn_phys%smc)) then
      deallocate( mn_phys%slmsk )
      deallocate( mn_phys%smc )
      deallocate( mn_phys%stc )
      deallocate( mn_phys%slc )

      deallocate( mn_phys%sfalb_lnd )
      deallocate( mn_phys%emis_lnd )
      deallocate( mn_phys%emis_ice )
      deallocate( mn_phys%emis_wat )
      deallocate( mn_phys%sfalb_lnd_bck )

      !deallocate( mn_phys%semis )
      !deallocate( mn_phys%semisbase )
      !deallocate( mn_phys%sfalb )

      deallocate( mn_phys%u10m )
      deallocate( mn_phys%v10m )
      deallocate( mn_phys%tprcp )

      deallocate( mn_phys%hprime )

      deallocate( mn_phys%zorl )
      deallocate( mn_phys%zorll )
      deallocate( mn_phys%zorlwav )
      deallocate( mn_phys%zorlw )

      deallocate( mn_phys%alvsf )
      deallocate( mn_phys%alvwf )
      deallocate( mn_phys%alnsf )
      deallocate( mn_phys%alnwf )

      deallocate( mn_phys%facsf )
      deallocate( mn_phys%facwf )

      deallocate( mn_phys%lakefrac )
      deallocate( mn_phys%lakedepth )

      deallocate( mn_phys%canopy )
      deallocate( mn_phys%vegfrac )
      deallocate( mn_phys%uustar )
      deallocate( mn_phys%shdmin )
      deallocate( mn_phys%shdmax )
      deallocate( mn_phys%tsfco )
      deallocate( mn_phys%tsfcl )
      deallocate( mn_phys%tsfc )
      !deallocate( mn_phys%tsfc_radtime )

      deallocate( mn_phys%albdirvis_lnd )
      deallocate( mn_phys%albdirnir_lnd )
      deallocate( mn_phys%albdifvis_lnd )
      deallocate( mn_phys%albdifnir_lnd )

      deallocate( mn_phys%cv )
      deallocate( mn_phys%cvt )
      deallocate( mn_phys%cvb )

      deallocate( mn_phys%phy_f2d )
      deallocate( mn_phys%phy_f3d )
    end if

    ! if move_nsst
    if (allocated( mn_phys%tref )) then
      deallocate( mn_phys%tref )
      deallocate( mn_phys%z_c )
      deallocate( mn_phys%c_0 )
      deallocate( mn_phys%c_d )
      deallocate( mn_phys%w_0 )
      deallocate( mn_phys%w_d )
      deallocate( mn_phys%xt )
      deallocate( mn_phys%xs )
      deallocate( mn_phys%xu )
      deallocate( mn_phys%xv )
      deallocate( mn_phys%xz )
      deallocate( mn_phys%zm )
      deallocate( mn_phys%xtts )
      deallocate( mn_phys%xzts )
      deallocate( mn_phys%d_conv )
      !deallocate( mn_phys%ifd )
      deallocate( mn_phys%dt_cool )
      deallocate( mn_phys%qrain )
    end if

  end subroutine deallocate_fv_moving_nest_physics_type

#endif
end module fv_moving_nest_types_mod
