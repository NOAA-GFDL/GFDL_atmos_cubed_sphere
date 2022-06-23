module coarse_grained_diagnostics_mod

  use constants_mod, only: rdgas, grav, pi=>pi_8
  use diag_manager_mod, only: diag_axis_init, register_diag_field, register_static_field, send_data
  use field_manager_mod,  only: MODEL_ATMOS
  use fv_arrays_mod, only: fv_atmos_type, fv_coarse_graining_type
  use fv_diagnostics_mod, only: cs3_interpolator, get_height_given_pressure
  use fv_mapz_mod, only: moist_cp, moist_cv
  use mpp_domains_mod, only: domain2d, EAST, NORTH
  use mpp_mod, only: FATAL, mpp_error
  use coarse_graining_mod, only: block_sum, get_fine_array_bounds, get_coarse_array_bounds, MODEL_LEVEL, &
                                 weighted_block_average, PRESSURE_LEVEL, vertically_remap_field, &
                                 vertical_remapping_requirements, mask_area_weights, mask_mass_weights, &
                                 block_edge_sum_x, block_edge_sum_y
  use time_manager_mod, only: time_type
  use tracer_manager_mod, only: get_tracer_index, get_tracer_names

  implicit none
  private

  type data_subtype
    real, dimension(:,:),   pointer :: var2 => null()
    real, dimension(:,:,:), pointer :: var3 => null()
  end type data_subtype

  type coarse_diag_type
    integer :: id = -99
    integer :: axes  ! 2 or 3, depending on whether the variable is 2D or 3D
    character(len=64) :: module_name
    character(len=128) :: name
    character(len=128) :: description
    character(len=64) :: units
    character(len=64) :: reduction_method
    logical :: vertically_integrated = .false.
    logical :: scaled_by_specific_heat_and_vertically_integrated = .false.
    logical :: always_model_level_coarse_grain = .false.
    integer :: pressure_level = -1  ! If greater than 0, interpolate to this pressure level (in hPa)
    integer :: iv = 0  ! Controls type of pressure-level interpolation performed (-1, 0, or 1)
    character(len=64) :: special_case  ! E.g. height is computed differently on pressure surfaces
    type(data_subtype) :: data
  end type coarse_diag_type

  public :: fv_coarse_diag_init, fv_coarse_diag

  integer :: tile_count = 1  ! Following fv_diagnostics.F90
  integer :: DIAG_SIZE = 512
  type(coarse_diag_type), dimension(512) :: coarse_diagnostics

  ! Reduction methods
  character(len=11) :: AREA_WEIGHTED = 'area_weighted'
  character(len=11) :: MASS_WEIGHTED = 'mass_weighted'
  character(len=5) :: pressure_level_label

contains

  subroutine populate_coarse_diag_type(Atm, coarse_diagnostics)
    type(fv_atmos_type), intent(in), target :: Atm(:)
    type(coarse_diag_type), intent(out) :: coarse_diagnostics(:)

    integer :: is, ie, js, je, npz, n_tracers, n_prognostic, t, p, n_pressure_levels
    integer :: index = 1
    character(len=128) :: tracer_name
    character(len=256) :: tracer_long_name, tracer_units
    character(len=8) :: DYNAMICS = 'dynamics'
    integer :: pressure_levels(31)

    n_pressure_levels = 31
    pressure_levels = (/1,2,3,5,7,10,20,30,50,70,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,925,950,975,1000/)
    npz = Atm(tile_count)%npz
    n_prognostic = size(Atm(tile_count)%q, 4)
    n_tracers = Atm(tile_count)%ncnst
    call get_fine_array_bounds(is, ie, js, je)

    coarse_diagnostics(index)%axes = 3
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'omega_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained pressure velocity'
    coarse_diagnostics(index)%units = 'Pa/s'
    coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED
    coarse_diagnostics(index)%data%var3 => Atm(tile_count)%omga(is:ie,js:je,1:npz)

    index = index + 1
    coarse_diagnostics(index)%axes = 3
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'ucomp_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained zonal wind'
    coarse_diagnostics(index)%units = 'm/s'
    coarse_diagnostics(index)%reduction_method = MASS_WEIGHTED
    coarse_diagnostics(index)%data%var3 => Atm(tile_count)%ua(is:ie,js:je,1:npz)

    index = index + 1
    coarse_diagnostics(index)%axes = 3
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'vcomp_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained meridional wind'
    coarse_diagnostics(index)%units = 'm/s'
    coarse_diagnostics(index)%reduction_method = MASS_WEIGHTED
    coarse_diagnostics(index)%data%var3 => Atm(tile_count)%va(is:ie,js:je,1:npz)

    index = index + 1
    coarse_diagnostics(index)%axes = 3
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'temp_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained temperature'
    coarse_diagnostics(index)%units = 'K'
    coarse_diagnostics(index)%reduction_method = MASS_WEIGHTED
    coarse_diagnostics(index)%data%var3 => Atm(tile_count)%pt(is:ie,js:je,1:npz)

    index = index + 1
    coarse_diagnostics(index)%axes = 3
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'delp_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained pressure thickness'
    coarse_diagnostics(index)%units = 'Pa'
    coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED
    coarse_diagnostics(index)%always_model_level_coarse_grain = .true.
    coarse_diagnostics(index)%data%var3 => Atm(tile_count)%delp(is:ie,js:je,1:npz)

    index = index + 1
    coarse_diagnostics(index)%axes = 2
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'ps_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained surface pressure'
    coarse_diagnostics(index)%units = 'Pa'
    coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED
    coarse_diagnostics(index)%data%var2 => Atm(tile_count)%ps(is:ie,js:je)

    if (.not. Atm(tile_count)%flagstruct%hydrostatic) then
       index = index + 1
       coarse_diagnostics(index)%axes = 3
       coarse_diagnostics(index)%module_name = DYNAMICS
       coarse_diagnostics(index)%name = 'delz_coarse'
       coarse_diagnostics(index)%description = 'coarse-grained height thickness'
       coarse_diagnostics(index)%units = 'm'
       coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED
       coarse_diagnostics(index)%always_model_level_coarse_grain = .true.
       coarse_diagnostics(index)%data%var3 => Atm(tile_count)%delz(is:ie,js:je,1:npz)

       index = index + 1
       coarse_diagnostics(index)%axes = 3
       coarse_diagnostics(index)%module_name = DYNAMICS
       coarse_diagnostics(index)%name = 'w_coarse'
       coarse_diagnostics(index)%description = 'coarse-grained vertical wind'
       coarse_diagnostics(index)%units = 'm/s'
       coarse_diagnostics(index)%reduction_method = MASS_WEIGHTED
       coarse_diagnostics(index)%data%var3 => Atm(tile_count)%w(is:ie,js:je,1:npz)
    endif

    do t = 1, n_tracers
      call get_tracer_names(MODEL_ATMOS, t, tracer_name, tracer_long_name, tracer_units)
      index = index + 1
      coarse_diagnostics(index)%axes = 3
      coarse_diagnostics(index)%module_name = DYNAMICS
      coarse_diagnostics(index)%name = trim(tracer_name) // '_coarse'
      coarse_diagnostics(index)%description = 'coarse-grained ' // trim(tracer_long_name)
      coarse_diagnostics(index)%units = tracer_units
      coarse_diagnostics(index)%reduction_method = MASS_WEIGHTED
      if (t .gt. n_prognostic) then
        coarse_diagnostics(index)%data%var3 => Atm(tile_count)%qdiag(is:ie,js:je,1:npz,t)
      else
        coarse_diagnostics(index)%data%var3 => Atm(tile_count)%q(is:ie,js:je,1:npz,t)
      endif
    enddo

    ! Defer pointer association for these diagnostics in case their arrays have
    ! not been allocated yet.
    index = index + 1
    coarse_diagnostics(index)%axes = 3
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'qv_dt_phys_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained water vapor specific humidity tendency from physics'
    coarse_diagnostics(index)%units = 'kg/kg/s'
    coarse_diagnostics(index)%reduction_method = MASS_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 3
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'ql_dt_phys_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained total liquid water tendency from physics'
    coarse_diagnostics(index)%units = 'kg/kg/s'
    coarse_diagnostics(index)%reduction_method = MASS_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 3
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'qi_dt_phys_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained total ice water tendency from physics'
    coarse_diagnostics(index)%units = 'kg/kg/s'
    coarse_diagnostics(index)%reduction_method = MASS_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 3
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'liq_wat_dt_phys_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained liquid water tracer tendency from physics'
    coarse_diagnostics(index)%units = 'kg/kg/s'
    coarse_diagnostics(index)%reduction_method = MASS_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 3
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'ice_wat_dt_phys_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained ice water tracer tendency from physics'
    coarse_diagnostics(index)%units = 'kg/kg/s'
    coarse_diagnostics(index)%reduction_method = MASS_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 3
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'qr_dt_phys_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained rain water tracer tendency from physics'
    coarse_diagnostics(index)%units = 'kg/kg/s'
    coarse_diagnostics(index)%reduction_method = MASS_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 3
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'qs_dt_phys_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained snow water tracer tendency from physics'
    coarse_diagnostics(index)%units = 'kg/kg/s'
    coarse_diagnostics(index)%reduction_method = MASS_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 3
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'qg_dt_phys_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained graupel tracer tendency from physics'
    coarse_diagnostics(index)%units = 'kg/kg/s'
    coarse_diagnostics(index)%reduction_method = MASS_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 3
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 't_dt_phys_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained temperature tendency from physics'
    coarse_diagnostics(index)%units = 'K/s'
    coarse_diagnostics(index)%reduction_method = MASS_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 3
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'u_dt_phys_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained zonal wind tendency from physics'
    coarse_diagnostics(index)%units = 'm/s/s'
    coarse_diagnostics(index)%reduction_method = MASS_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 3
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'v_dt_phys_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained meridional wind tendency from physics'
    coarse_diagnostics(index)%units = 'm/s/s'
    coarse_diagnostics(index)%reduction_method = MASS_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 3
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 't_dt_nudge_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained temperature tendency from nudging'
    coarse_diagnostics(index)%units = 'K/s'
    coarse_diagnostics(index)%reduction_method = MASS_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 2
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'ps_dt_nudge_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained surface pressure tendency from nudging'
    coarse_diagnostics(index)%units = 'Pa/s'
    coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 3
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'delp_dt_nudge_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained pressure thickness tendency from nudging'
    coarse_diagnostics(index)%units = 'Pa/s'
    coarse_diagnostics(index)%always_model_level_coarse_grain = .true.
    coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 3
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'u_dt_nudge_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained zonal wind tendency from nudging'
    coarse_diagnostics(index)%units = 'm/s/s'
    coarse_diagnostics(index)%reduction_method = MASS_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 3
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'v_dt_nudge_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained meridional wind tendency from nudging'
    coarse_diagnostics(index)%units = 'm/s/s'
    coarse_diagnostics(index)%reduction_method = MASS_WEIGHTED

    ! Vertically integrated diagnostics
    index = index + 1
    coarse_diagnostics(index)%axes = 2
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'int_qv_dt_phys_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained vertically integrated water vapor specific humidity tendency from physics'
    coarse_diagnostics(index)%units = 'kg/m**2/s'
    coarse_diagnostics(index)%vertically_integrated = .true.
    coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 2
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'int_ql_dt_phys_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained vertically integrated total liquid water tendency from physics'
    coarse_diagnostics(index)%units = 'kg/m**2/s'
    coarse_diagnostics(index)%vertically_integrated = .true.
    coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 2
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'int_qi_dt_phys_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained vertically integrated total ice water tendency from physics'
    coarse_diagnostics(index)%units = 'kg/m**2/s'
    coarse_diagnostics(index)%vertically_integrated = .true.
    coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 2
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'int_liq_wat_dt_phys_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained vertically integrated liquid water tracer tendency from physics'
    coarse_diagnostics(index)%units = 'kg/m**2/s'
    coarse_diagnostics(index)%vertically_integrated = .true.
    coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 2
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'int_ice_wat_dt_phys_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained vertically integrated ice water tracer tendency from physics'
    coarse_diagnostics(index)%units = 'kg/m**2/s'
    coarse_diagnostics(index)%vertically_integrated = .true.
    coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 2
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'int_qr_dt_phys_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained vertically integrated rain water tracer tendency from physics'
    coarse_diagnostics(index)%units = 'kg/m**2/s'
    coarse_diagnostics(index)%vertically_integrated = .true.
    coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 2
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'int_qs_dt_phys_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained vertically integrated snow water tracer tendency from physics'
    coarse_diagnostics(index)%units = 'kg/m**2/s'
    coarse_diagnostics(index)%vertically_integrated = .true.
    coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 2
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'int_qg_dt_phys_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained vertically integrated graupel tracer tendency from physics'
    coarse_diagnostics(index)%units = 'kg/m**2/s'
    coarse_diagnostics(index)%vertically_integrated = .true.
    coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 2
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'int_t_dt_phys_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained vertically integrated temperature tendency from physics'
    coarse_diagnostics(index)%units = 'W/m**2'
    coarse_diagnostics(index)%scaled_by_specific_heat_and_vertically_integrated = .true.
    coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 2
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'int_u_dt_phys_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained vertically integrated zonal wind tendency from physics'
    coarse_diagnostics(index)%units = 'kg/m s/s'
    coarse_diagnostics(index)%vertically_integrated = .true.
    coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 2
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'int_v_dt_phys_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained vertically integrated meridional wind tendency from physics'
    coarse_diagnostics(index)%units = 'kg/m s/s'
    coarse_diagnostics(index)%vertically_integrated = .true.
    coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 2
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'int_t_dt_nudge_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained vertically integrated temperature tendency from nudging'
    coarse_diagnostics(index)%units = 'W/m**2'
    coarse_diagnostics(index)%scaled_by_specific_heat_and_vertically_integrated = .true.
    coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 2
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'int_u_dt_nudge_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained vertically integrated zonal wind tendency from nudging'
    coarse_diagnostics(index)%units = 'kg/m s/s'
    coarse_diagnostics(index)%vertically_integrated = .true.
    coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED

    index = index + 1
    coarse_diagnostics(index)%axes = 2
    coarse_diagnostics(index)%module_name = DYNAMICS
    coarse_diagnostics(index)%name = 'int_v_dt_nudge_coarse'
    coarse_diagnostics(index)%description = 'coarse-grained vertically integrated meridional wind tendency from nudging'
    coarse_diagnostics(index)%units = 'kg/m s/s'
    coarse_diagnostics(index)%vertically_integrated = .true.
    coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED

    ! iv =-1: winds
    ! iv = 0: positive definite scalars
    ! iv = 1: temperature
    do p = 1, n_pressure_levels
      ! Note all reference data for pressure-level variables is 3D, but the diagnostics
      ! themselves are 2D.
      write(pressure_level_label, '(I5)') pressure_levels(p)

      index = index + 1
      coarse_diagnostics(index)%pressure_level = pressure_levels(p)
      coarse_diagnostics(index)%axes = 2
      coarse_diagnostics(index)%module_name = DYNAMICS
      coarse_diagnostics(index)%name = 'u' // trim(adjustl(pressure_level_label)) // '_coarse'
      coarse_diagnostics(index)%description = 'coarse-grained ' // trim(adjustl(pressure_level_label)) // '-mb u'
      coarse_diagnostics(index)%units = 'm/s'
      coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED
      coarse_diagnostics(index)%data%var3 => Atm(tile_count)%ua(is:ie,js:je,1:npz)
      coarse_diagnostics(index)%iv = -1

      index = index + 1
      coarse_diagnostics(index)%pressure_level = pressure_levels(p)
      coarse_diagnostics(index)%axes = 2
      coarse_diagnostics(index)%module_name = DYNAMICS
      coarse_diagnostics(index)%name = 'v' // trim(adjustl(pressure_level_label)) // '_coarse'
      coarse_diagnostics(index)%description = 'coarse-grained ' // trim(adjustl(pressure_level_label)) // '-mb v'
      coarse_diagnostics(index)%units = 'm/s'
      coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED
      coarse_diagnostics(index)%data%var3 => Atm(tile_count)%va(is:ie,js:je,1:npz)
      coarse_diagnostics(index)%iv = -1

      index = index + 1
      coarse_diagnostics(index)%pressure_level = pressure_levels(p)
      coarse_diagnostics(index)%axes = 2
      coarse_diagnostics(index)%module_name = DYNAMICS
      coarse_diagnostics(index)%name = 't' // trim(adjustl(pressure_level_label)) // '_coarse'
      coarse_diagnostics(index)%description = 'coarse-grained ' // trim(adjustl(pressure_level_label)) // '-mb temperature'
      coarse_diagnostics(index)%units = 'K'
      coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED
      coarse_diagnostics(index)%data%var3 => Atm(tile_count)%pt(is:ie,js:je,1:npz)
      coarse_diagnostics(index)%iv = 1

      index = index + 1
      coarse_diagnostics(index)%pressure_level = pressure_levels(p)
      coarse_diagnostics(index)%axes = 2
      coarse_diagnostics(index)%module_name = DYNAMICS
      coarse_diagnostics(index)%name = 'omg' // trim(adjustl(pressure_level_label)) // '_coarse'
      coarse_diagnostics(index)%description = 'coarse-grained ' // trim(adjustl(pressure_level_label)) // '-mb omega'
      coarse_diagnostics(index)%units = 'Pa/s'
      coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED
      coarse_diagnostics(index)%data%var3 => Atm(tile_count)%omga(is:ie,js:je,1:npz)
      coarse_diagnostics(index)%iv = -1

      index = index + 1
      coarse_diagnostics(index)%pressure_level = pressure_levels(p)
      coarse_diagnostics(index)%axes = 2
      coarse_diagnostics(index)%module_name = DYNAMICS
      coarse_diagnostics(index)%name = 'z' // trim(adjustl(pressure_level_label)) // '_coarse'
      coarse_diagnostics(index)%description = 'coarse-grained ' // trim(adjustl(pressure_level_label)) // '-mb height'
      coarse_diagnostics(index)%units = 'm'
      coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED
      coarse_diagnostics(index)%special_case = 'height'

      do t = 1, n_tracers
        call get_tracer_names(MODEL_ATMOS, t, tracer_name, tracer_long_name, tracer_units)
        index = index + 1
        coarse_diagnostics(index)%pressure_level = pressure_levels(p)
        coarse_diagnostics(index)%axes = 2
        coarse_diagnostics(index)%module_name = DYNAMICS
        if (trim(tracer_name) .eq. 'sphum') then
           coarse_diagnostics(index)%name = 'q' // trim(adjustl(pressure_level_label)) // '_coarse'
        else
           coarse_diagnostics(index)%name = trim(tracer_name) // trim(adjustl(pressure_level_label)) // '_coarse'
        endif
        coarse_diagnostics(index)%description = 'coarse-grained ' // trim(adjustl(pressure_level_label)) // '-mb ' // trim(tracer_long_name)
        coarse_diagnostics(index)%units = tracer_units
        coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED
        if (t .gt. n_prognostic) then
          coarse_diagnostics(index)%data%var3 => Atm(tile_count)%qdiag(is:ie,js:je,1:npz,t)
        else
          coarse_diagnostics(index)%data%var3 => Atm(tile_count)%q(is:ie,js:je,1:npz,t)
        endif
        coarse_diagnostics(index)%iv = 0
      enddo

      if (.not. Atm(tile_count)%flagstruct%hydrostatic) then
         index = index + 1
         coarse_diagnostics(index)%pressure_level = pressure_levels(p)
         coarse_diagnostics(index)%axes = 2
         coarse_diagnostics(index)%module_name = DYNAMICS
         coarse_diagnostics(index)%name = 'w' // trim(adjustl(pressure_level_label)) // '_coarse'
         coarse_diagnostics(index)%description = 'coarse-grained ' // trim(adjustl(pressure_level_label)) // '-mb vertical wind'
         coarse_diagnostics(index)%units = 'm/s'
         coarse_diagnostics(index)%reduction_method = AREA_WEIGHTED
         coarse_diagnostics(index)%data%var3 => Atm(tile_count)%w(is:ie,js:je,1:npz)
         coarse_diagnostics(index)%iv = -1
      endif
    enddo
  end subroutine populate_coarse_diag_type

  subroutine register_coarse_diagnostics(Atm, coarse_diagnostics, Time, &
       id_xt_coarse, id_yt_coarse, id_pfull_coarse, id_x_coarse, id_y_coarse)
    type(fv_atmos_type), intent(inout) :: Atm(:)
    type(coarse_diag_type), intent(inout) :: coarse_diagnostics(:)
    type(time_type), intent(in) :: Time
    integer, intent(in) :: id_xt_coarse, id_yt_coarse, id_pfull_coarse
    integer, intent(in) :: id_x_coarse, id_y_coarse

    integer :: index, n_valid_diagnostics
    integer :: axes_t(3), axes(3)
    real :: missing_value = -1.0e10  ! Following fv_diagnostics.F90

    axes_t = (/  id_xt_coarse, id_yt_coarse, id_pfull_coarse /)
    axes = (/  id_x_coarse, id_y_coarse, id_pfull_coarse /)
    do index = 1, DIAG_SIZE
      if (trim(coarse_diagnostics(index)%name) == '') exit
      n_valid_diagnostics = index
    enddo

    do index = 1, n_valid_diagnostics
      coarse_diagnostics(index)%id = register_diag_field( &
        trim(coarse_diagnostics(index)%module_name), &
        trim(coarse_diagnostics(index)%name), &
        axes_t(1:coarse_diagnostics(index)%axes), &
        Time, &
        trim(coarse_diagnostics(index)%description), &
        trim(coarse_diagnostics(index)%units), &
        missing_value=missing_value &
      )
      call maybe_allocate_reference_array(Atm, coarse_diagnostics(index))
    enddo

    call register_coarse_static_diagnostics(Atm, Time, axes_t, axes)
  end subroutine register_coarse_diagnostics

  ! Some diagnostics may only have memory allocated for them if they are requested
  subroutine maybe_allocate_reference_array(Atm, coarse_diagnostic)
    type(fv_atmos_type), target, intent(inout) :: Atm(:)
    type(coarse_diag_type), intent(inout) :: coarse_diagnostic

    integer :: is, ie, js, je, npz

    call get_fine_array_bounds(is, ie, js, je)
    npz = Atm(tile_count)%npz

    ! It would be really nice if there were a cleaner way to do this;
    ! unfortunately it is not possible to check if an array associated with a
    ! pointer is allocated.
    if (coarse_diagnostic%id .gt. 0) then
       if (ends_with(coarse_diagnostic%name, 'qv_dt_phys_coarse')) then
          if (.not. allocated(Atm(tile_count)%phys_diag%phys_qv_dt)) then
             allocate(Atm(tile_count)%phys_diag%phys_qv_dt(is:ie,js:je,1:npz))
          endif
          coarse_diagnostic%data%var3 => Atm(tile_count)%phys_diag%phys_qv_dt(is:ie,js:je,1:npz)
       elseif (ends_with(coarse_diagnostic%name, 'ql_dt_phys_coarse')) then
          if (.not. allocated(Atm(tile_count)%phys_diag%phys_ql_dt)) then
             allocate(Atm(tile_count)%phys_diag%phys_ql_dt(is:ie,js:je,1:npz))
          endif
          coarse_diagnostic%data%var3 => Atm(tile_count)%phys_diag%phys_ql_dt(is:ie,js:je,1:npz)
       elseif (ends_with(coarse_diagnostic%name, 'qi_dt_phys_coarse')) then
          if (.not. allocated(Atm(tile_count)%phys_diag%phys_qi_dt)) then
             allocate(Atm(tile_count)%phys_diag%phys_qi_dt(is:ie,js:je,1:npz))
          endif
          coarse_diagnostic%data%var3 => Atm(tile_count)%phys_diag%phys_qi_dt(is:ie,js:je,1:npz)
       elseif (ends_with(coarse_diagnostic%name, 'liq_wat_dt_phys_coarse')) then
          if (.not. allocated(Atm(tile_count)%phys_diag%phys_liq_wat_dt)) then
             allocate(Atm(tile_count)%phys_diag%phys_liq_wat_dt(is:ie,js:je,1:npz))
          endif
          coarse_diagnostic%data%var3 => Atm(tile_count)%phys_diag%phys_liq_wat_dt(is:ie,js:je,1:npz)
       elseif (ends_with(coarse_diagnostic%name, 'ice_wat_dt_phys_coarse')) then
          if (.not. allocated(Atm(tile_count)%phys_diag%phys_ice_wat_dt)) then
             allocate(Atm(tile_count)%phys_diag%phys_ice_wat_dt(is:ie,js:je,1:npz))
          endif
          coarse_diagnostic%data%var3 => Atm(tile_count)%phys_diag%phys_ice_wat_dt(is:ie,js:je,1:npz)
       elseif (ends_with(coarse_diagnostic%name, 'qr_dt_phys_coarse')) then
          if (.not. allocated(Atm(tile_count)%phys_diag%phys_qr_dt)) then
             allocate(Atm(tile_count)%phys_diag%phys_qr_dt(is:ie,js:je,1:npz))
          endif
          coarse_diagnostic%data%var3 => Atm(tile_count)%phys_diag%phys_qr_dt(is:ie,js:je,1:npz)
       elseif (ends_with(coarse_diagnostic%name, 'qg_dt_phys_coarse')) then
          if (.not. allocated(Atm(tile_count)%phys_diag%phys_qg_dt)) then
             allocate(Atm(tile_count)%phys_diag%phys_qg_dt(is:ie,js:je,1:npz))
          endif
          coarse_diagnostic%data%var3 => Atm(tile_count)%phys_diag%phys_qg_dt(is:ie,js:je,1:npz)
       elseif (ends_with(coarse_diagnostic%name, 'qs_dt_phys_coarse')) then
          if (.not. allocated(Atm(tile_count)%phys_diag%phys_qs_dt)) then
             allocate(Atm(tile_count)%phys_diag%phys_qs_dt(is:ie,js:je,1:npz))
          endif
          coarse_diagnostic%data%var3 => Atm(tile_count)%phys_diag%phys_qs_dt(is:ie,js:je,1:npz)
       elseif (ends_with(coarse_diagnostic%name, 't_dt_phys_coarse')) then
          if (.not. allocated(Atm(tile_count)%phys_diag%phys_t_dt)) then
             allocate(Atm(tile_count)%phys_diag%phys_t_dt(is:ie,js:je,1:npz))
          endif
          coarse_diagnostic%data%var3 => Atm(tile_count)%phys_diag%phys_t_dt(is:ie,js:je,1:npz)
       elseif (ends_with(coarse_diagnostic%name, 'u_dt_phys_coarse')) then
          if (.not. allocated(Atm(tile_count)%phys_diag%phys_u_dt)) then
             allocate(Atm(tile_count)%phys_diag%phys_u_dt(is:ie,js:je,1:npz))
          endif
          coarse_diagnostic%data%var3 => Atm(tile_count)%phys_diag%phys_u_dt(is:ie,js:je,1:npz)
       elseif (ends_with(coarse_diagnostic%name, 'v_dt_phys_coarse')) then
          if (.not. allocated(Atm(tile_count)%phys_diag%phys_v_dt)) then
             allocate(Atm(tile_count)%phys_diag%phys_v_dt(is:ie,js:je,1:npz))
          endif
          coarse_diagnostic%data%var3 => Atm(tile_count)%phys_diag%phys_v_dt(is:ie,js:je,1:npz)
       elseif (ends_with(coarse_diagnostic%name, 't_dt_nudge_coarse')) then
          if (.not. allocated(Atm(tile_count)%nudge_diag%nudge_t_dt)) then
             allocate(Atm(tile_count)%nudge_diag%nudge_t_dt(is:ie,js:je,1:npz))
             Atm(tile_count)%nudge_diag%nudge_t_dt(is:ie,js:je,1:npz) = 0.0
          endif
          coarse_diagnostic%data%var3 => Atm(tile_count)%nudge_diag%nudge_t_dt(is:ie,js:je,1:npz)
       elseif (trim(coarse_diagnostic%name) .eq. 'ps_dt_nudge_coarse') then
          if (.not. allocated(Atm(tile_count)%nudge_diag%nudge_ps_dt)) then
             allocate(Atm(tile_count)%nudge_diag%nudge_ps_dt(is:ie,js:je))
             Atm(tile_count)%nudge_diag%nudge_ps_dt(is:ie,js:je) = 0.0
          endif
          coarse_diagnostic%data%var2 => Atm(tile_count)%nudge_diag%nudge_ps_dt(is:ie,js:je)
       elseif (trim(coarse_diagnostic%name) .eq. 'delp_dt_nudge_coarse') then
          if (.not. allocated(Atm(tile_count)%nudge_diag%nudge_delp_dt)) then
             allocate(Atm(tile_count)%nudge_diag%nudge_delp_dt(is:ie,js:je,1:npz))
             Atm(tile_count)%nudge_diag%nudge_delp_dt(is:ie,js:je,1:npz) = 0.0
          endif
          coarse_diagnostic%data%var3 => Atm(tile_count)%nudge_diag%nudge_delp_dt(is:ie,js:je,1:npz)
       elseif (ends_with(coarse_diagnostic%name, 'u_dt_nudge_coarse')) then
          if (.not. allocated(Atm(tile_count)%nudge_diag%nudge_u_dt)) then
             allocate(Atm(tile_count)%nudge_diag%nudge_u_dt(is:ie,js:je,1:npz))
             Atm(tile_count)%nudge_diag%nudge_u_dt(is:ie,js:je,1:npz) = 0.0
          endif
          coarse_diagnostic%data%var3 => Atm(tile_count)%nudge_diag%nudge_u_dt(is:ie,js:je,1:npz)
       elseif (ends_with(coarse_diagnostic%name, 'v_dt_nudge_coarse')) then
          if (.not. allocated(Atm(tile_count)%nudge_diag%nudge_v_dt)) then
             allocate(Atm(tile_count)%nudge_diag%nudge_v_dt(is:ie,js:je,1:npz))
             Atm(tile_count)%nudge_diag%nudge_v_dt(is:ie,js:je,1:npz) = 0.0
          endif
          coarse_diagnostic%data%var3 => Atm(tile_count)%nudge_diag%nudge_v_dt(is:ie,js:je,1:npz)
       endif
    endif
  end subroutine maybe_allocate_reference_array

  subroutine fv_coarse_diag_init(Atm, Time, id_pfull, id_phalf, coarse_graining)
    type(fv_atmos_type), intent(inout) :: Atm(:)
    type(time_type), intent(in) :: Time
    integer, intent(in) :: id_pfull, id_phalf
    type(fv_coarse_graining_type), intent(inout) :: coarse_graining

    integer :: is, ie, js, je, is_coarse, ie_coarse, js_coarse, je_coarse

    call get_fine_array_bounds(is, ie, js, je)
    call get_coarse_array_bounds(is_coarse, ie_coarse, js_coarse, je_coarse)
    call initialize_coarse_diagnostic_axes(coarse_graining%domain, coarse_graining%nx_coarse, &
         coarse_graining%id_x_coarse, coarse_graining%id_y_coarse, coarse_graining%id_xt_coarse, &
         coarse_graining%id_yt_coarse)

    coarse_graining%id_pfull = id_pfull
    coarse_graining%id_phalf = id_phalf

    call populate_coarse_diag_type(Atm, coarse_diagnostics)
    call register_coarse_diagnostics(Atm, coarse_diagnostics, Time, &
         coarse_graining%id_xt_coarse, coarse_graining%id_yt_coarse, id_pfull, &
         coarse_graining%id_x_coarse, coarse_graining%id_y_coarse)
  end subroutine fv_coarse_diag_init

  subroutine initialize_coarse_diagnostic_axes(coarse_domain, &
    nx_coarse, id_x_coarse, id_y_coarse, id_xt_coarse, id_yt_coarse)
    type(domain2d), intent(in) :: coarse_domain
    integer, intent(in) :: nx_coarse
    integer, intent(inout) :: id_x_coarse, id_y_coarse, id_xt_coarse, id_yt_coarse

    integer :: i, j
    real, allocatable :: grid_x_coarse(:), grid_y_coarse(:), grid_xt_coarse(:), grid_yt_coarse(:)

    allocate(grid_x_coarse(nx_coarse + 1))
    allocate(grid_y_coarse(nx_coarse + 1))
    allocate(grid_xt_coarse(nx_coarse))
    allocate(grid_yt_coarse(nx_coarse))

    grid_x_coarse = (/ (i, i=1, nx_coarse + 1) /)
    grid_y_coarse = (/ (j, j=1, nx_coarse + 1) /)
    grid_xt_coarse = (/ (i, i=1, nx_coarse) /)
    grid_yt_coarse = (/ (j, j=1, nx_coarse) /)

    id_xt_coarse = diag_axis_init('grid_xt_coarse', grid_xt_coarse, &
         'index', 'x', 'x-index of cell center points', set_name='coarse_grid', &
         Domain2=coarse_domain, tile_count=tile_count)
    id_yt_coarse = diag_axis_init('grid_yt_coarse', grid_yt_coarse, &
         'index', 'y', 'y-index of cell center points', set_name='coarse_grid', &
         Domain2=coarse_domain, tile_count=tile_count)

    id_x_coarse = diag_axis_init('grid_x_coarse', grid_x_coarse, &
         'index', 'x', 'x-index of cell corner points', set_name='coarse_grid', &
         Domain2=coarse_domain, tile_count=tile_count, domain_position=EAST)
    id_y_coarse = diag_axis_init('grid_y_coarse', grid_y_coarse, &
         'index', 'y', 'y-index of cell corner points', set_name='coarse_grid', &
         Domain2=coarse_domain, tile_count=tile_count, domain_position=NORTH)
  end subroutine initialize_coarse_diagnostic_axes

  subroutine fv_coarse_diag(Atm, Time)
    type(fv_atmos_type), intent(in), target :: Atm(:)
    type(time_type), intent(in) :: Time

    real, allocatable :: work_2d(:,:), work_2d_coarse(:,:), work_3d_coarse(:,:,:)
    real, allocatable :: mass(:,:,:), height_on_interfaces(:,:,:), masked_area(:,:,:)
    real, allocatable :: phalf(:,:,:), upsampled_coarse_phalf(:,:,:)
    integer :: is, ie, js, je, is_coarse, ie_coarse, js_coarse, je_coarse, npz
    logical :: used
    logical :: need_2d_work_array, need_3d_work_array, need_mass_array, need_height_array, need_masked_area_array
    integer :: index, i, j
    character(len=256) :: error_message

    call get_need_nd_work_array(2, need_2d_work_array)
    call get_need_nd_work_array(3, need_3d_work_array)
    call get_need_mass_array(need_mass_array)
    call get_need_height_array(need_height_array)

    if (trim(Atm(tile_count)%coarse_graining%strategy) .eq. PRESSURE_LEVEL) then
      call get_need_masked_area_array(need_masked_area_array)
    else
      need_masked_area_array = .false.
    endif

    call get_fine_array_bounds(is, ie, js, je)
    call get_coarse_array_bounds(is_coarse, ie_coarse, js_coarse, je_coarse)
    npz = Atm(tile_count)%npz

    if (need_2d_work_array) then
      allocate(work_2d_coarse(is_coarse:ie_coarse,js_coarse:je_coarse))
    endif

    if (need_3d_work_array) then
       allocate(work_3d_coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz))
       if (trim(Atm(tile_count)%coarse_graining%strategy) .eq. PRESSURE_LEVEL) then
        allocate(phalf(is:ie,js:je,1:npz+1))
        allocate(upsampled_coarse_phalf(is:ie,js:je,1:npz+1))

        call vertical_remapping_requirements( &
              Atm(tile_count)%delp(is:ie,js:je,1:npz), &
              Atm(tile_count)%gridstruct%area(is:ie,js:je), &
              Atm(tile_count)%ptop, &
              phalf, &
              upsampled_coarse_phalf)
       endif
    endif

    if (need_mass_array) then
      allocate(mass(is:ie,js:je,1:npz))
      if (trim(Atm(tile_count)%coarse_graining%strategy) .eq. MODEL_LEVEL) then
        call compute_mass(Atm(tile_count), is, ie, js, je, npz, mass)
      else if (trim(Atm(tile_count)%coarse_graining%strategy) .eq. PRESSURE_LEVEL) then
        call mask_mass_weights( &
             Atm(tile_count)%gridstruct%area(is:ie,js:je), &
             Atm(tile_count)%delp(is:ie,js:je,1:npz), &
             phalf, &
             upsampled_coarse_phalf, &
             mass)
      endif
    endif

    if (need_masked_area_array) then
      allocate(masked_area(is:ie,js:je,1:npz))
      call mask_area_weights( &
           Atm(tile_count)%gridstruct%area(is:ie,js:je), &
           phalf, &
           upsampled_coarse_phalf, &
           masked_area)
    endif

    if (need_height_array) then
      allocate(height_on_interfaces(is:ie,js:je,1:npz+1))
      if(Atm(tile_count)%flagstruct%hydrostatic) then
        call compute_height_on_interfaces_hydrostatic( &
          is, &
          ie, &
          js, &
          je, &
          npz, &
          Atm(tile_count)%pt(is:ie,js:je,1:npz), &
          Atm(tile_count)%peln(is:ie,1:npz+1,js:je), &
          height_on_interfaces(is:ie,js:je,1:npz) &
        )
      else
        call compute_height_on_interfaces_nonhydrostatic( &
          is, &
          ie, &
          js, &
          je, &
          npz, &
          Atm(tile_count)%delz(is:ie,js:je,1:npz), &
          height_on_interfaces(is:ie,js:je,1:npz) &
        )
      endif
      if (.not. allocated(work_2d_coarse)) allocate(work_2d_coarse(is_coarse:ie_coarse,js_coarse:je_coarse))
      allocate(work_2d(is:ie,js:je))
    endif

    do index = 1, DIAG_SIZE
      if (coarse_diagnostics(index)%id .gt. 0) then
        if (coarse_diagnostics(index)%axes .eq. 2) then
          call coarse_grain_2D_field(is, ie, js, je, npz, is_coarse, ie_coarse, js_coarse, je_coarse, &
                                     Atm(tile_count), coarse_diagnostics(index), height_on_interfaces, work_2d_coarse)
          used = send_data(coarse_diagnostics(index)%id, work_2d_coarse, Time)
       elseif (coarse_diagnostics(index)%axes .eq. 3) then
          if (trim(Atm(tile_count)%coarse_graining%strategy) .eq. MODEL_LEVEL .or. coarse_diagnostics(index)%always_model_level_coarse_grain) then
            call coarse_grain_3D_field_on_model_levels(is, ie, js, je, is_coarse, ie_coarse, js_coarse, je_coarse, npz, &
                                                       coarse_diagnostics(index), Atm(tile_count)%gridstruct%area(is:ie,js:je),&
                                                       mass, work_3d_coarse)
          else if (trim(Atm(tile_count)%coarse_graining%strategy) .eq. PRESSURE_LEVEL) then
             call coarse_grain_3D_field_on_pressure_levels(is, ie, js, je, is_coarse, ie_coarse, js_coarse, je_coarse, npz, &
                                                          coarse_diagnostics(index), masked_area, mass, phalf, &
                                                          upsampled_coarse_phalf, Atm(tile_count)%ptop, work_3d_coarse)
          else
            write(error_message, *) 'fv_coarse_diag: invalid coarse-graining strategy provided for 3D variables, ' // &
            trim(Atm(tile_count)%coarse_graining%strategy)
            call mpp_error(FATAL, error_message)
          endif
          used = send_data(coarse_diagnostics(index)%id, work_3d_coarse, Time)
        endif
      endif
    enddo
  end subroutine fv_coarse_diag

   subroutine coarse_grain_3D_field_on_model_levels(is, ie, js, je, is_coarse, ie_coarse, js_coarse, je_coarse, &
                                                    npz, coarse_diag, area, mass, result)
    integer, intent(in) :: is, ie, js, je, is_coarse, ie_coarse, js_coarse, je_coarse, npz
    type(coarse_diag_type) :: coarse_diag
    real, intent(in) :: mass(is:ie,js:je,1:npz), area(is:ie,js:je)
    real, intent(out) :: result(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz)

    character(len=256) :: error_message

    if (trim(coarse_diag%reduction_method) .eq. AREA_WEIGHTED) then
      call weighted_block_average( &
        area(is:ie,js:je), &
        coarse_diag%data%var3, &
        result &
      )
    elseif (trim(coarse_diag%reduction_method) .eq. MASS_WEIGHTED) then
      call weighted_block_average( &
        mass(is:ie,js:je,1:npz), &
        coarse_diag%data%var3, &
        result &
      )
    else
      write(error_message, *) 'coarse_grain_3D_field_on_model_levels: invalid reduction_method, ' // &
        trim(coarse_diag%reduction_method) // ', provided for 3D variable, ' // &
        trim(coarse_diag%name)
      call mpp_error(FATAL, error_message)
    endif
   end subroutine coarse_grain_3D_field_on_model_levels

   subroutine coarse_grain_3D_field_on_pressure_levels(is, ie, js, je, is_coarse, ie_coarse, js_coarse, je_coarse, &
                                                       npz, coarse_diag, masked_area, masked_mass, phalf, upsampled_coarse_phalf, &
                                                       ptop, result)
    integer, intent(in) :: is, ie, js, je, is_coarse, ie_coarse, js_coarse, je_coarse, npz
    type(coarse_diag_type) :: coarse_diag
    real, intent(in) :: masked_mass(is:ie,js:je,1:npz), masked_area(is:ie,js:je,1:npz)
    real, intent(in) :: phalf(is:ie,js:je,1:npz+1), upsampled_coarse_phalf(is:ie,js:je,1:npz+1)
    real, intent(in) :: ptop
    real, intent(out) :: result(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz)

    real, allocatable :: remapped_field(:,:,:)
    character(len=256) :: error_message

    allocate(remapped_field(is:ie,js:je,1:npz))

    call vertically_remap_field( &
      phalf, &
      coarse_diag%data%var3, &
      upsampled_coarse_phalf, &
      ptop, &
      remapped_field)
    if (trim(coarse_diag%reduction_method) .eq. AREA_WEIGHTED) then
      call weighted_block_average( &
        masked_area(is:ie,js:je,1:npz), &
        remapped_field(is:ie,js:je,1:npz), &
        result &
      )
    elseif (trim(coarse_diag%reduction_method) .eq. MASS_WEIGHTED) then
      call weighted_block_average( &
        masked_mass(is:ie,js:je,1:npz), &
        remapped_field(is:ie,js:je,1:npz), &
        result &
      )
    else
      write(error_message, *) 'coarse_grain_3D_field_on_pressure_levels: invalid reduction_method, ' // &
        trim(coarse_diag%reduction_method) // ', provided for 3D variable, ' // &
        trim(coarse_diag%name)
      call mpp_error(FATAL, error_message)
    endif
   end subroutine coarse_grain_3D_field_on_pressure_levels

   subroutine coarse_grain_2D_field(is, ie, js, je, npz, is_coarse, ie_coarse, js_coarse, je_coarse, &
                                    Atm, coarse_diag, height_on_interfaces, result)
    integer, intent(in) :: is, ie, js, je, npz, is_coarse, ie_coarse, js_coarse, je_coarse
    type(fv_atmos_type), intent(in) :: Atm
    type(coarse_diag_type), intent(in) :: coarse_diag
    real, intent(in) :: height_on_interfaces(is:ie,js:je,1:npz+1)
    real, intent(out) :: result(is_coarse:ie_coarse,js_coarse:je_coarse)

    character(len=256) :: error_message
    real, allocatable :: work_2d(:,:)

    if (coarse_diag%pressure_level > 0 .or. coarse_diag%vertically_integrated &
         .or. coarse_diag%scaled_by_specific_heat_and_vertically_integrated) then
      allocate(work_2d(is:ie,js:je))
    endif

    if (trim(coarse_diag%reduction_method) .eq. AREA_WEIGHTED) then
      if (coarse_diag%pressure_level < 0 &
         .and. .not. coarse_diag%vertically_integrated &
         .and. .not. coarse_diag%scaled_by_specific_heat_and_vertically_integrated) then
        call weighted_block_average( &
          Atm%gridstruct%area(is:ie,js:je), &
          coarse_diag%data%var2, &
          result &
        )
      elseif (trim(coarse_diag%special_case) .eq. 'height') then
        call height_given_pressure_level( &
          is, &
          ie, &
          js, &
          je, &
          npz, &
          height_on_interfaces(is:ie,js:je,1:npz+1), &
          Atm%peln(is:ie,1:npz+1,js:je), &
          coarse_diag%pressure_level, &
          work_2d(is:ie,js:je) &
        )
        call weighted_block_average( &
          Atm%gridstruct%area(is:ie,js:je), &
          work_2d, &
          result &
          )
      elseif (coarse_diag%vertically_integrated) then
        call vertically_integrate( &
             is, &
             ie, &
             js, &
             je, &
             npz, &
             Atm%delp(is:ie,js:je,1:npz), &
             coarse_diag%data%var3, &
             work_2d(is:ie,js:je))
        call weighted_block_average( &
          Atm%gridstruct%area(is:ie,js:je), &
          work_2d, &
          result &
          )
      elseif (coarse_diag%scaled_by_specific_heat_and_vertically_integrated) then
        call scale_by_specific_heat_and_vertically_integrate( &
             is, &
             ie, &
             js, &
             je, &
             npz, &
             Atm, &
             coarse_diag%data%var3, &
             work_2d(is:ie,js:je))
        call weighted_block_average( &
          Atm%gridstruct%area(is:ie,js:je), &
          work_2d, &
          result &
          )
      else
        call interpolate_to_pressure_level( &
          is, &
          ie, &
          js, &
          je, &
          npz, &
          coarse_diag%data%var3, &
          height_on_interfaces(is:ie,js:je,1:npz+1), &
          Atm%peln(is:ie,1:npz+1,js:je), &
          coarse_diag%pressure_level, &
          coarse_diag%iv, &
          work_2d &
        )
        call weighted_block_average( &
          Atm%gridstruct%area(is:ie,js:je), &
          work_2d, &
          result &
        )
      endif
    else
      write(error_message, *) 'coarse_grain_2D_field: invalid reduction_method, ' // &
        trim(coarse_diag%reduction_method) // ', provided for 2D variable, ' // &
        trim(coarse_diag%name)
      call mpp_error(FATAL, error_message)
    endif
   end subroutine coarse_grain_2D_field

   subroutine get_need_nd_work_array(dimension, need_nd_work_array)
     integer, intent(in) :: dimension
     logical, intent(out) :: need_nd_work_array

     integer :: index

     need_nd_work_array = .false.
     do index = 1, DIAG_SIZE
       if ((coarse_diagnostics(index)%axes == dimension) .and. (coarse_diagnostics(index)%id > 0)) then
         need_nd_work_array = .true.
         exit
       endif
     enddo
   end subroutine get_need_nd_work_array

   subroutine get_need_mass_array(need_mass_array)
     logical, intent(out) :: need_mass_array

     integer :: index

     need_mass_array = .false.
     do index = 1, DIAG_SIZE
       if ((coarse_diagnostics(index)%axes == 3) .and. &
           (trim(coarse_diagnostics(index)%reduction_method) .eq. MASS_WEIGHTED) .and. &
           (coarse_diagnostics(index)%id > 0)) then
           need_mass_array = .true.
           exit
       endif
     enddo
  end subroutine get_need_mass_array

  ! If we are interpolating the surfaces of constant pressure, we need
  ! to compute the height on model level interfaces.
  subroutine get_need_height_array(need_height_array)
    logical, intent(out) :: need_height_array

    integer :: index

    need_height_array = .false.
    do index = 1, DIAG_SIZE
      if ((coarse_diagnostics(index)%axes == 2) .and. &
          (coarse_diagnostics(index)%pressure_level > 0) .and. &
          (coarse_diagnostics(index)%id > 0)) then
          need_height_array = .true.
          exit
      endif
    enddo
 end subroutine get_need_height_array

  subroutine get_need_masked_area_array(need_masked_area_array)
    logical, intent(out) :: need_masked_area_array

    integer :: index

    need_masked_area_array = .false.
    do index = 1, DIAG_SIZE
      if ((coarse_diagnostics(index)%axes == 3) .and. &
          (trim(coarse_diagnostics(index)%reduction_method) .eq. AREA_WEIGHTED) .and. &
          (coarse_diagnostics(index)%id > 0)) then
         need_masked_area_array = .true.
         exit
      endif
   enddo
 end subroutine get_need_masked_area_array

  subroutine compute_mass(Atm, is, ie, js, je, npz, mass)
    type(fv_atmos_type), intent(in) :: Atm
    integer, intent(in) :: is, ie, js, je, npz
    real, intent(out) :: mass(is:ie,js:je,1:npz)

    integer :: k

    do k = 1, npz
      mass(is:ie,js:je,k) = Atm%delp(is:ie,js:je,k) * Atm%gridstruct%area(is:ie,js:je)
    enddo
  end subroutine compute_mass

  subroutine interpolate_to_pressure_level(is, ie, js, je, npz, field, height, phalf, pressure_level, iv, result)
    integer, intent(in) :: is, ie, js, je, npz, iv
    real, intent(in) :: field(is:ie,js:je,1:npz), height(is:ie,js:je,1:npz+1), phalf(is:ie,1:npz+1,js:je)
    integer, intent(in) :: pressure_level
    real, intent(out) :: result(is:ie,js:je)

    real, allocatable :: work(:,:,:)
    integer :: n_pressure_levels = 1
    real :: output_pressures(1)
    integer :: ids(1) = 1  ! Set > 0

    output_pressures = log(100.0 * real(pressure_level))  ! convert to Pa then take log to match expectation of cs3_interpolator
    allocate(work(is:ie,js:je,n_pressure_levels))

    call cs3_interpolator(is, ie, js, je, npz, field, n_pressure_levels, output_pressures, height, phalf, ids, work, iv)
    result = work(is:ie,js:je,1)
  end subroutine interpolate_to_pressure_level

  subroutine compute_height_on_interfaces_hydrostatic(is, ie, js, je, npz, temperature, phalf, height)
    integer, intent(in) :: is, ie, js, je, npz
    real, intent(in) :: temperature(is:ie,js:je,1:npz), phalf(is:ie,1:npz+1,js:je)
    real, intent(out) :: height(is:ie,js:je,1:npz+1)

    integer :: i, j, k
    real :: rgrav

    rgrav = 1.0 / grav

    do j = js, je
      do i = is, ie
        height(i,j,npz+1) = 0.0
        do k = npz, 1, -1
          height(i,j,k) = height(i,j,k+1) - (rdgas / grav) * temperature(i,j,k) * (phalf(i,k,j) - phalf(i,k+1,j))
        enddo
      enddo
    enddo
  end subroutine compute_height_on_interfaces_hydrostatic

  subroutine compute_height_on_interfaces_nonhydrostatic(is, ie, js, je, npz, delz, height)
    integer, intent(in) :: is, ie, js, je, npz
    real, intent(in) :: delz(is:ie,js:je,1:npz)
    real, intent(out) :: height(is:ie,js:je,1:npz+1)

    integer :: i, j, k

    do j = js, je
      do i = is, ie
        height(i,j,npz+1) = 0.0
        do k = npz, 1, -1
          height(i,j,k) = height(i,j,k+1) - delz(i,j,k)
        enddo
      enddo
    enddo
  end subroutine compute_height_on_interfaces_nonhydrostatic

  subroutine height_given_pressure_level(is, ie, js, je, npz, height, phalf, pressure_level, result)
    integer, intent(in) :: is, ie, js, je, npz, pressure_level
    real, intent(in) :: height(is:ie,js:je,1:npz+1), phalf(is:ie,1:npz+1,js:je)
    real, intent(out) :: result(is:ie,js:je)

    real, allocatable :: work(:,:,:)
    integer :: n_pressure_levels = 1
    real :: output_pressures(1)
    integer :: ids(1) = 1  ! Set > 0

    output_pressures = log(100 * real(pressure_level))
    allocate(work(is:ie,js:je,n_pressure_levels))

    call get_height_given_pressure(is, ie, js, je, npz, height, n_pressure_levels, ids, output_pressures, phalf, work)
    result(is:ie,js:je) = work(is:ie,js:je,1)
  end subroutine height_given_pressure_level

  function starts_with(string, prefix)
    character(len=64), intent(in) :: string, prefix
    logical :: starts_with

    starts_with = string(1:len(trim(prefix))) .eq. trim(prefix)
    return
  end function starts_with

  function ends_with(string, suffix)
    character(len=64), intent(in) :: string
    character(len=*), intent(in) :: suffix
    logical :: ends_with

    integer :: string_length, suffix_length

    string_length = len(trim(string))
    suffix_length = len(trim(suffix))
    if (string_length .lt. suffix_length) then
       ends_with = .false.
    else
       ends_with = string(string_length - suffix_length + 1:string_length) .eq. trim(suffix)
    endif
    return
  end function ends_with

  subroutine vertically_integrate(is, ie, js, je, npz, delp, field, integrated_field)
    integer, intent(in) :: is, ie, js, je, npz
    real, intent(in) :: delp(is:ie,js:je,1:npz), field(is:ie,js:je,1:npz)
    real, intent(out) :: integrated_field(is:ie,js:je)

    integrated_field = sum(delp * field, dim=3) / grav
  end subroutine vertically_integrate

  subroutine scale_by_specific_heat_and_vertically_integrate(is, ie, js, je, npz, Atm, field, integrated_field)
    integer, intent(in) :: is, ie, js, je, npz
    type(fv_atmos_type) :: Atm
    real, intent(in) :: field(is:ie,js:je,1:npz)
    real, intent(out) :: integrated_field(is:ie,js:je)

    real, allocatable, dimension(:,:,:) :: specific_heat

    allocate(specific_heat(is:ie,js:je,1:npz))

    if (.not. Atm%flagstruct%hydrostatic) then
       call compute_cvm(Atm%q, Atm%pt, is, ie, js, je, npz, Atm%bd%isd, Atm%bd%ied, Atm%bd%jsd, Atm%bd%jed, Atm%flagstruct%nwat, specific_heat)
    else
       call compute_cpm(Atm%q, Atm%pt, is, ie, js, je, npz, Atm%bd%isd, Atm%bd%ied, Atm%bd%jsd, Atm%bd%jed, Atm%flagstruct%nwat, specific_heat)
    endif
    integrated_field = sum(specific_heat * Atm%delp(is:ie,js:je,1:npz) * field, dim=3) / grav
  end subroutine scale_by_specific_heat_and_vertically_integrate

  subroutine compute_cvm(q, pt, isc, iec, jsc, jec, npz, isd, ied, jsd, jed, nwat, cvm)
    integer :: isc, iec, jsc, jec, npz, isd, ied, jsd, jed, nwat
    real, dimension(isd:ied,jsd:jed,1:npz,1:nwat), intent(in) :: q
    real, dimension(isd:ied,jsd:jed,1:npz), intent(in) :: pt
    real, dimension(isc:iec,jsc:jec,1:npz), intent(out) :: cvm
    real, dimension(isc:iec) :: qc, cvm_tmp
    integer :: j, k, sphum, liq_wat, ice_wat, rainwat, snowwat, graupel
    sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
    liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
    ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
    rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat')
    snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat')
    graupel = get_tracer_index (MODEL_ATMOS, 'graupel')
    do j = jsc, jec
       do k = 1, npz
          call moist_cv(isc, iec, isd, ied, jsd, jed, npz, j, k, nwat, sphum, &
               liq_wat, rainwat, ice_wat, snowwat, graupel, &
               q, qc, cvm_tmp, pt(isc:iec,j,k))
          cvm(isc:iec,j,k) = cvm_tmp
       enddo
    enddo
  end subroutine compute_cvm

 subroutine compute_cpm(q, pt, isc, iec, jsc, jec, npz, isd, ied, jsd, jed, nwat, cpm)
    integer :: isc, iec, jsc, jec, npz, isd, ied, jsd, jed, nwat
    real, dimension(isd:ied,jsd:jed,1:npz,1:nwat), intent(in) :: q
    real, dimension(isd:ied,jsd:jed,1:npz), intent(in) :: pt
    real, dimension(isc:iec,jsc:jec,1:npz), intent(out) :: cpm
    real, dimension(isc:iec) :: qc, cpm_tmp
    integer :: j, k, sphum, liq_wat, ice_wat, rainwat, snowwat, graupel
    sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
    liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
    ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
    rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat')
    snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat')
    graupel = get_tracer_index (MODEL_ATMOS, 'graupel')
    do j = jsc, jec
       do k = 1, npz
          call moist_cp(isc, iec, isd, ied, jsd, jed, npz, j, k, nwat, sphum, &
               liq_wat, rainwat, ice_wat, snowwat, graupel, &
               q, qc, cpm_tmp, pt(isc:iec,j,k))
          cpm(isc:iec,j,k) = cpm_tmp
       enddo
    enddo
  end subroutine compute_cpm

  subroutine register_coarse_static_diagnostics(Atm, Time, axes_t, axes)
    type(fv_atmos_type), intent(in) :: Atm(:)
    type(time_type), intent(in) :: Time
    integer, intent(in) :: axes_t(3), axes(3)

    integer :: id_area_coarse, id_dx_coarse, id_dy_coarse, id_grid_lon_coarse
    integer :: id_grid_lat_coarse, id_grid_lont_coarse, id_grid_latt_coarse
    integer :: id_zsurf_coarse
    integer :: is, ie, js, je
    integer :: is_coarse, ie_coarse, js_coarse, je_coarse
    logical :: used
    integer :: tile_count = 1
    character(len=8) :: DYNAMICS = 'dynamics'
    real :: rad2deg = 180. / pi


    real, allocatable, dimension(:,:,:) :: grid_coarse, gridt_coarse
    real, allocatable, dimension(:,:) :: work_2d_coarse

    call get_fine_array_bounds(is, ie, js, je)
    call get_coarse_array_bounds(is_coarse, ie_coarse, js_coarse, je_coarse)

    id_grid_lon_coarse = register_static_field(DYNAMICS, 'grid_lon_coarse', &
         axes(1:2), 'longitude', 'degrees_E')
    id_grid_lat_coarse = register_static_field(DYNAMICS, 'grid_lat_coarse', &
         axes(1:2), 'latitude', 'degrees_N')
    id_grid_lont_coarse = register_static_field(DYNAMICS, 'grid_lont_coarse', &
         axes_t(1:2), 'longitude', 'degrees_E')
    id_grid_latt_coarse = register_static_field(DYNAMICS, 'grid_latt_coarse', &
         axes_t(1:2), 'latitude', 'degrees_N')
    id_dx_coarse = register_static_field(DYNAMICS, 'dx_coarse', &
         (/ axes_t(1), axes(2) /), 'dx', 'm')
    id_dy_coarse = register_static_field(DYNAMICS, 'dy_coarse', &
         (/ axes(1), axes_t(2) /), 'dy', 'm')
    id_area_coarse = register_static_field(DYNAMICS, 'area_coarse', &
         axes_t(1:2), 'cell area', 'm**2')
    id_zsurf_coarse = register_static_field(DYNAMICS, 'zsurf_coarse', &
         axes_t(1:2), 'surface height', 'm')

    if (id_grid_lont_coarse .gt. 0 .or. id_grid_latt_coarse .gt. 0) then
       allocate(gridt_coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:2))
       call compute_gridt_coarse(is, ie, js, je, is_coarse, ie_coarse, js_coarse, je_coarse, Atm, gridt_coarse)
       if (id_grid_lont_coarse .gt. 0) used = send_data(id_grid_lont_coarse, rad2deg * gridt_coarse(:,:,1), Time)
       if (id_grid_latt_coarse .gt. 0) used = send_data(id_grid_latt_coarse, rad2deg * gridt_coarse(:,:,2), Time)
    endif
    if (id_grid_lon_coarse .gt. 0 .or. id_grid_lat_coarse .gt. 0) then
       allocate(grid_coarse(is_coarse:ie_coarse+1,js_coarse:je_coarse+1,1:2))
       call compute_grid_coarse(is, ie, js, je, is_coarse, ie_coarse, js_coarse, je_coarse, Atm, grid_coarse)
       if (id_grid_lon_coarse .gt. 0) used = send_data(id_grid_lon_coarse, rad2deg * grid_coarse(:,:,1), Time)
       if (id_grid_lat_coarse .gt. 0) used = send_data(id_grid_lat_coarse, rad2deg * grid_coarse(:,:,2), Time)
    endif
    if (id_area_coarse .gt. 0) then
       allocate(work_2d_coarse(is_coarse:ie_coarse,js_coarse:je_coarse))
       call block_sum(Atm(tile_count)%gridstruct%area(is:ie,js:je), work_2d_coarse)
       used = send_data(id_area_coarse, work_2d_coarse, Time)
       deallocate(work_2d_coarse)
    endif
    if (id_zsurf_coarse .gt. 0) then
       allocate(work_2d_coarse(is_coarse:ie_coarse,js_coarse:je_coarse))
       call weighted_block_average(Atm(tile_count)%gridstruct%area(is:ie,js:je), &
            Atm(tile_count)%phis(is:ie,js:je) / GRAV, work_2d_coarse)
       used = send_data(id_zsurf_coarse, work_2d_coarse, Time)
       deallocate(work_2d_coarse)
    endif
    if (id_dx_coarse .gt. 0) then
       allocate(work_2d_coarse(is_coarse:ie_coarse,js_coarse:je_coarse+1))
       call block_edge_sum_x(Atm(tile_count)%gridstruct%dx(is:ie,js:je+1), work_2d_coarse)
       used = send_data(id_dx_coarse, work_2d_coarse, Time)
       deallocate(work_2d_coarse)
    endif
    if (id_dy_coarse .gt. 0) then
       allocate(work_2d_coarse(is_coarse:ie_coarse+1,js_coarse:je_coarse))
       call block_edge_sum_y(Atm(tile_count)%gridstruct%dy(is:ie+1,js:je), work_2d_coarse)
       used = send_data(id_dy_coarse, work_2d_coarse, Time)
       deallocate(work_2d_coarse)
    endif
  end subroutine register_coarse_static_diagnostics

  subroutine compute_gridt_coarse(is, ie, js, je, is_coarse, ie_coarse, &
       js_coarse, je_coarse, Atm, gridt_coarse)
    integer, intent(in) :: is, ie, js, je, is_coarse, ie_coarse, js_coarse, je_coarse
    type(fv_atmos_type), intent(in) :: Atm(:)
    real, intent(out) :: gridt_coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:2)

    integer :: factor, offset
    integer :: tile_count = 1

    factor = Atm(tile_count)%coarse_graining%factor
    offset = factor / 2
    if (mod(factor, 2) .eq. 0) then
       gridt_coarse = Atm(tile_count)%gridstruct%grid(is+offset:ie+1:factor,js+offset:je+1:factor,:)
    else
       gridt_coarse = Atm(tile_count)%gridstruct%grid(is+offset:ie:factor,js+offset:je:factor,:)
    endif
  end subroutine compute_gridt_coarse

  subroutine compute_grid_coarse(is, ie, js, je, is_coarse, ie_coarse, &
       js_coarse, je_coarse, Atm, grid_coarse)
    integer, intent(in) :: is, ie, js, je, is_coarse, ie_coarse, js_coarse, je_coarse
    type(fv_atmos_type), intent(in) :: Atm(:)
    real, intent(out) :: grid_coarse(is_coarse:ie_coarse+1,js_coarse:je_coarse+1,1:2)

    integer :: factor, offset
    integer :: tile_count = 1

    factor = Atm(tile_count)%coarse_graining%factor
    grid_coarse = Atm(tile_count)%gridstruct%grid(is:ie+1:factor,js:je+1:factor,:)
  end subroutine compute_grid_coarse
end module coarse_grained_diagnostics_mod
