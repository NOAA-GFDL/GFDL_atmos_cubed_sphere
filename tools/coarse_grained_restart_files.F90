module coarse_grained_restart_files_mod

  use coarse_graining_mod, only: compute_mass_weights, get_coarse_array_bounds,&
       get_fine_array_bounds, MODEL_LEVEL, PRESSURE_LEVEL, weighted_block_average, &
       weighted_block_edge_average_x, weighted_block_edge_average_y, &
       mask_area_weights, mask_mass_weights, block_upsample, remap_edges_along_x, &
       remap_edges_along_y, vertically_remap_field
  use constants_mod, only: GRAV, RDGAS, RVGAS
  use field_manager_mod, only: MODEL_ATMOS
  use fms_io_mod,      only: register_restart_field, save_restart
  use fv_arrays_mod, only: coarse_restart_type, fv_atmos_type
  use mpp_domains_mod, only: domain2d, EAST, NORTH, mpp_update_domains
  use mpp_mod, only: FATAL, mpp_error
  use tracer_manager_mod, only: get_tracer_names, get_tracer_index, set_tracer_profile

  implicit none
  private

  public :: deallocate_coarse_restart_type, fv_coarse_restart_init, fv_io_write_restart_coarse

  ! Global variables for this module, initialized in fv_coarse_restart_init
  integer :: is, ie, js, je, npz
  integer :: is_coarse, ie_coarse, js_coarse, je_coarse
  integer :: n_prognostic_tracers, n_diagnostic_tracers, n_tracers

contains

  subroutine fv_coarse_restart_init(tile_count, nz, nt_prog, &
       nt_phys, hydrostatic, hybrid_z, fv_land, &
       write_coarse_dgrid_vel_rst, write_coarse_agrid_vel_rst, &
       coarse_domain, restart)
    integer, intent(in) :: tile_count, nz, nt_prog, nt_phys
    logical, intent(in) :: hydrostatic, hybrid_z, fv_land
    logical, intent(in) :: write_coarse_dgrid_vel_rst, write_coarse_agrid_vel_rst
    type(domain2d), intent(inout) :: coarse_domain
    type(coarse_restart_type), intent(inout) :: restart

    call get_fine_array_bounds(is, ie, js, je)
    call get_coarse_array_bounds(is_coarse, ie_coarse, js_coarse, je_coarse)
    n_prognostic_tracers = nt_prog
    n_diagnostic_tracers = nt_phys
    n_tracers = nt_prog + nt_phys
    npz = nz

    call allocate_coarse_restart_type(hydrostatic, hybrid_z, &
         fv_land, write_coarse_dgrid_vel_rst, write_coarse_agrid_vel_rst, &
         restart)
    call register_coarse_restart_files(tile_count, hydrostatic, &
         hybrid_z, fv_land, write_coarse_dgrid_vel_rst, &
         write_coarse_agrid_vel_rst, coarse_domain, restart)
  end subroutine fv_coarse_restart_init

  subroutine fv_io_write_restart_coarse(Atm, timestamp)
    type(fv_atmos_type), intent(inout) :: Atm
    character(len=*), optional, intent(in) :: timestamp

    integer :: tile_count, n_tiles

    call coarse_grain_restart_data(Atm)
    call save_restart(Atm%coarse_graining%restart%fv_core_coarse, timestamp)
    call save_restart(Atm%coarse_graining%restart%fv_tracer_coarse, timestamp)
    call save_restart(Atm%coarse_graining%restart%fv_srf_wnd_coarse, timestamp)
    if (Atm%flagstruct%fv_land) then
       call save_restart(Atm%coarse_graining%restart%mg_drag_coarse, timestamp)
       call save_restart(Atm%coarse_graining%restart%fv_land_coarse, timestamp)
    endif
  end subroutine fv_io_write_restart_coarse

  subroutine allocate_coarse_restart_type(hydrostatic, hybrid_z, &
       fv_land, write_coarse_dgrid_vel_rst, write_coarse_agrid_vel_rst, restart)
    logical, intent(in) :: hydrostatic, hybrid_z, fv_land
    logical, intent(in) :: write_coarse_dgrid_vel_rst, write_coarse_agrid_vel_rst
    type(coarse_restart_type), intent(inout) :: restart

    if (write_coarse_dgrid_vel_rst) then
       allocate(restart%u(is_coarse:ie_coarse,js_coarse:je_coarse+1,npz))
       allocate(restart%v(is_coarse:ie_coarse+1,js_coarse:je_coarse,npz))
    endif

    allocate(restart%u_srf(is_coarse:ie_coarse,js_coarse:je_coarse))
    allocate(restart%v_srf(is_coarse:ie_coarse,js_coarse:je_coarse))
    allocate(restart%delp(is_coarse:ie_coarse,js_coarse:je_coarse,npz))
    allocate(restart%pt(is_coarse:ie_coarse,js_coarse:je_coarse,npz))
    allocate(restart%q(is_coarse:ie_coarse,js_coarse:je_coarse,npz,n_prognostic_tracers))
    allocate(restart%qdiag(is_coarse:ie_coarse,js_coarse:je_coarse,npz,n_prognostic_tracers+1:n_tracers))
    allocate(restart%phis(is_coarse:ie_coarse,js_coarse:je_coarse))

    if (write_coarse_agrid_vel_rst) then
       allocate(restart%ua(is_coarse:ie_coarse,js_coarse:je_coarse,npz))
       allocate(restart%va(is_coarse:ie_coarse,js_coarse:je_coarse,npz))
    endif

    if (.not. hydrostatic) then
       allocate(restart%w(is_coarse:ie_coarse,js_coarse:je_coarse,npz))
       allocate(restart%delz(is_coarse:ie_coarse,js_coarse:je_coarse,npz))
       if (hybrid_z) allocate(restart%ze0(is_coarse:ie_coarse,js_coarse:je_coarse,npz))
    endif

    if (fv_land) then
       allocate(restart%sgh(is_coarse:ie_coarse,js_coarse:je_coarse))
       allocate(restart%oro(is_coarse:ie_coarse,js_coarse:je_coarse))
    endif
  end subroutine allocate_coarse_restart_type

  subroutine deallocate_coarse_restart_type(restart)
    type(coarse_restart_type), intent(inout) :: restart

    if (allocated(restart%u)) deallocate(restart%u)
    if (allocated(restart%v)) deallocate(restart%v)
    if (allocated(restart%w)) deallocate(restart%w)
    if (allocated(restart%pt)) deallocate(restart%pt)
    if (allocated(restart%delp)) deallocate(restart%delp)
    if (allocated(restart%delz)) deallocate(restart%delz)
    if (allocated(restart%ua)) deallocate(restart%ua)
    if (allocated(restart%va)) deallocate(restart%va)
    if (allocated(restart%phis)) deallocate(restart%phis)
    if (allocated(restart%q)) deallocate(restart%q)
    if (allocated(restart%qdiag)) deallocate(restart%qdiag)
    if (allocated(restart%u_srf)) deallocate(restart%u_srf)
    if (allocated(restart%v_srf)) deallocate(restart%v_srf)
    if (allocated(restart%sgh)) deallocate(restart%sgh)
    if (allocated(restart%oro)) deallocate(restart%oro)
    if (allocated(restart%ze0)) deallocate(restart%ze0)
  end subroutine deallocate_coarse_restart_type

  subroutine register_coarse_restart_files(tile_count, hydrostatic, &
       hybrid_z, fv_land, write_coarse_dgrid_vel_rst, write_coarse_agrid_vel_rst, &
       coarse_domain, restart)
    integer, intent(in) :: tile_count
    logical, intent(in) :: hydrostatic, hybrid_z, fv_land
    logical, intent(in) :: write_coarse_dgrid_vel_rst, write_coarse_agrid_vel_rst
    type(domain2d), intent(in) :: coarse_domain
    type(coarse_restart_type), intent(inout) :: restart

    call register_fv_core_coarse(tile_count, hydrostatic, hybrid_z, &
         write_coarse_dgrid_vel_rst, write_coarse_agrid_vel_rst, &
         coarse_domain, restart)
    call register_fv_tracer_coarse(tile_count, coarse_domain, restart)
    call register_fv_srf_wnd_coarse(tile_count, coarse_domain, restart)
    if (fv_land) then
       call register_mg_drag_coarse(tile_count, coarse_domain, restart)
       call register_fv_land_coarse(tile_count, coarse_domain, restart)
    endif
  end subroutine register_coarse_restart_files

  subroutine register_fv_core_coarse(tile_count, hydrostatic, hybrid_z, &
       write_coarse_dgrid_vel_rst, write_coarse_agrid_vel_rst, coarse_domain, &
       restart)
    integer, intent(in) :: tile_count
    logical, intent(in) :: hydrostatic, hybrid_z
    logical, intent(in) :: write_coarse_dgrid_vel_rst, write_coarse_agrid_vel_rst
    type(domain2d), intent(in) :: coarse_domain
    type(coarse_restart_type), intent(inout) :: restart

    character(len=64) :: filename
    integer :: id_restart

    filename = 'fv_core_coarse.res.nc'

    if (write_coarse_dgrid_vel_rst) then
       id_restart = register_restart_field(restart%fv_core_coarse, &
            filename, 'u', restart%u, domain=coarse_domain, position=NORTH, &
            tile_count=tile_count)
       id_restart = register_restart_field(restart%fv_core_coarse, &
            filename, 'v', restart%v, domain=coarse_domain, position=EAST, &
            tile_count=tile_count)
    endif

    if (write_coarse_agrid_vel_rst) then
       id_restart = register_restart_field(restart%fv_core_coarse, &
            filename, 'ua', restart%ua, domain=coarse_domain, tile_count=tile_count)
       id_restart = register_restart_field(restart%fv_core_coarse, &
            filename, 'va', restart%va, domain=coarse_domain, tile_count=tile_count)
    endif

    if (.not. hydrostatic) then
       id_restart = register_restart_field(restart%fv_core_coarse, &
            filename, 'W', restart%w, domain=coarse_domain, mandatory=.false., tile_count=tile_count)
       id_restart = register_restart_field(restart%fv_core_coarse, &
            filename, 'DZ', restart%delz, domain=coarse_domain, mandatory=.false., tile_count=tile_count)
       if (hybrid_z) then
          id_restart = register_restart_field(restart%fv_core_coarse, &
            filename, 'ZE0', restart%ze0, domain=coarse_domain, mandatory=.false., tile_count=tile_count)
       endif
    endif

    id_restart = register_restart_field(restart%fv_core_coarse, &
         filename, 'T', restart%pt, domain=coarse_domain, tile_count=tile_count)
    id_restart = register_restart_field(restart%fv_core_coarse, &
            filename, 'delp', restart%delp, domain=coarse_domain, tile_count=tile_count)
    id_restart = register_restart_field(restart%fv_core_coarse, &
            filename, 'phis', restart%phis, domain=coarse_domain, tile_count=tile_count)
  end subroutine register_fv_core_coarse

  subroutine register_fv_tracer_coarse(tile_count, coarse_domain, restart)
    integer, intent(in) :: tile_count
    type(domain2d), intent(in) :: coarse_domain
    type(coarse_restart_type), intent(inout) :: restart

    character(len=64) :: filename, tracer_name
    integer :: id_restart, n_tracer

    filename = 'fv_tracer_coarse.res.nc'

    do n_tracer = 1, n_prognostic_tracers
       call get_tracer_names(MODEL_ATMOS, n_tracer, tracer_name)
       call set_tracer_profile(MODEL_ATMOS, n_tracer, restart%q(:,:,:,n_tracer))
       id_restart = register_restart_field(restart%fv_tracer_coarse, &
            filename, tracer_name, restart%q(:,:,:,n_tracer), domain=coarse_domain, &
            mandatory=.false., tile_count=tile_count)
    enddo

    do n_tracer = n_prognostic_tracers + 1, n_tracers
       call get_tracer_names(MODEL_ATMOS, n_tracer, tracer_name)
       call set_tracer_profile(MODEL_ATMOS, n_tracer, restart%qdiag(:,:,:,n_tracer))
       id_restart = register_restart_field(restart%fv_tracer_coarse, &
            filename, tracer_name, restart%qdiag(:,:,:,n_tracer), domain=coarse_domain, &
            mandatory=.false., tile_count=tile_count)
    enddo
  end subroutine register_fv_tracer_coarse

  subroutine register_fv_srf_wnd_coarse(tile_count, coarse_domain, restart)
    integer, intent(in) :: tile_count
    type(domain2d), intent(in) :: coarse_domain
    type(coarse_restart_type), intent(inout) :: restart

    character(len=64) :: filename
    integer :: id_restart

    filename = 'fv_srf_wnd_coarse.res.nc'

    id_restart = register_restart_field(restart%fv_srf_wnd_coarse, &
         filename, 'u_srf', restart%u_srf, domain=coarse_domain, &
         tile_count=tile_count)
    id_restart = register_restart_field(restart%fv_srf_wnd_coarse, &
         filename, 'v_srf', restart%v_srf, domain=coarse_domain, &
         tile_count=tile_count)
  end subroutine register_fv_srf_wnd_coarse

  subroutine register_mg_drag_coarse(tile_count, coarse_domain, restart)
    integer, intent(in) :: tile_count
    type(domain2d), intent(in) :: coarse_domain
    type(coarse_restart_type), intent(out) :: restart

    character(len=64) :: filename
    integer :: id_restart

    filename = 'mg_drag_coarse.res.nc'

    id_restart = register_restart_field(restart%mg_drag_coarse, &
         filename, 'ghprime', restart%sgh, domain=coarse_domain, &
         tile_count=tile_count)
  end subroutine register_mg_drag_coarse

  subroutine register_fv_land_coarse(tile_count, coarse_domain, restart)
    integer, intent(in) :: tile_count
    type(domain2d), intent(in) :: coarse_domain
    type(coarse_restart_type), intent(inout) :: restart

    character(len=64) :: filename
    integer :: id_restart

    filename = 'fv_land_coarse.res.nc'

    id_restart = register_restart_field(restart%fv_land_coarse, &
         filename, 'oro', restart%oro, domain=coarse_domain, &
         tile_count=tile_count)
  end subroutine register_fv_land_coarse

  subroutine coarse_grain_restart_data(Atm)
    type(fv_atmos_type), intent(inout) :: Atm

    character(len=256) :: error_message

    if (trim(Atm%coarse_graining%strategy) .eq. MODEL_LEVEL) then
       call coarse_grain_restart_data_on_model_levels(Atm)
    elseif (trim(Atm%coarse_graining%strategy) .eq. PRESSURE_LEVEL) then
       call coarse_grain_restart_data_on_pressure_levels(Atm)
    else
       write(error_message, *) 'Currently only model_level and pressure_level coarse-graining are supported for restart files.'
       call mpp_error(FATAL, error_message)
    endif
  end subroutine coarse_grain_restart_data

  subroutine coarse_grain_restart_data_on_model_levels(Atm)
    type(fv_atmos_type), intent(inout) :: Atm

    real, allocatable :: mass(:,:,:)

    allocate(mass(is:ie,js:je,1:npz))
    call compute_mass_weights(Atm%gridstruct%area(is:ie,js:je), Atm%delp(is:ie,js:je,1:npz), mass)

    call coarse_grain_fv_core_restart_data_on_model_levels(Atm, mass)
    call coarse_grain_fv_tracer_restart_data_on_model_levels(Atm, mass)
    call coarse_grain_fv_srf_wnd_restart_data(Atm)
    if (Atm%flagstruct%fv_land) then
       call coarse_grain_mg_drag_restart_data(Atm)
       call coarse_grain_fv_land_restart_data(Atm)
    endif
  end subroutine coarse_grain_restart_data_on_model_levels

  subroutine coarse_grain_restart_data_on_pressure_levels(Atm)
     type(fv_atmos_type), intent(inout) :: Atm

     real, allocatable, dimension(:,:,:):: phalf, coarse_phalf, coarse_phalf_on_fine
     real, allocatable, dimension(:,:,:) :: masked_mass_weights, masked_area_weights

     allocate(phalf(is-1:ie+1,js-1:je+1,1:npz+1))  ! Require the halo here for the winds
     allocate(coarse_phalf(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz+1))
     allocate(coarse_phalf_on_fine(is:ie,js:je,1:npz+1))
     allocate(masked_mass_weights(is:ie,js:je,1:npz))
     allocate(masked_area_weights(is:ie,js:je,1:npz))

     ! delp and delz are coarse-grained on model levels; u, v, W, T, and all the tracers
     ! are all remapped to surfaces of constant pressure within coarse grid cells before
     ! coarse graining.  At the end, delz and phis are corrected to impose hydrostatic balance.
     call compute_pressure_level_coarse_graining_requirements( &
       Atm, phalf, coarse_phalf, coarse_phalf_on_fine, masked_mass_weights, masked_area_weights)
     call coarse_grain_fv_core_restart_data_on_pressure_levels( &
       Atm, phalf, coarse_phalf, coarse_phalf_on_fine, masked_mass_weights, masked_area_weights)
     call coarse_grain_fv_tracer_restart_data_on_pressure_levels( &
       Atm, phalf, coarse_phalf_on_fine, masked_mass_weights, masked_area_weights)
     call coarse_grain_fv_srf_wnd_restart_data(Atm)
     if (Atm%flagstruct%fv_land) then
       call coarse_grain_mg_drag_restart_data(Atm)
       call coarse_grain_fv_land_restart_data(Atm)
     endif
     call impose_hydrostatic_balance(Atm, coarse_phalf)
  end subroutine coarse_grain_restart_data_on_pressure_levels

  subroutine coarse_grain_fv_core_restart_data_on_model_levels(Atm, mass)
    type(fv_atmos_type), intent(inout) :: Atm
    real, intent(in) :: mass(is:ie,js:je,1:npz)

    if (Atm%coarse_graining%write_coarse_dgrid_vel_rst) then
       call weighted_block_edge_average_x(Atm%gridstruct%dx(is:ie,js:je+1), &
            Atm%u(is:ie,js:je+1,1:npz), Atm%coarse_graining%restart%u)
       call weighted_block_edge_average_y(Atm%gridstruct%dy(is:ie+1,js:je), &
            Atm%v(is:ie+1,js:je,1:npz), Atm%coarse_graining%restart%v)
    endif

    if (.not. Atm%flagstruct%hydrostatic) then
       call weighted_block_average(mass(is:ie,js:je,1:npz), &
            Atm%w(is:ie,js:je,1:npz), Atm%coarse_graining%restart%w)
       call weighted_block_average(Atm%gridstruct%area(is:ie,js:je), &
            Atm%delz(is:ie,js:je,1:npz), Atm%coarse_graining%restart%delz)
       if (Atm%flagstruct%hybrid_z) then
          call weighted_block_average(Atm%gridstruct%area(is:ie,js:je), &
            Atm%ze0(is:ie,js:je,1:npz), Atm%coarse_graining%restart%ze0)
       endif
    endif

    call weighted_block_average(mass(is:ie,js:je,1:npz), &
         Atm%pt(is:ie,js:je,1:npz), Atm%coarse_graining%restart%pt)
    call weighted_block_average(Atm%gridstruct%area(is:ie,js:je), &
         Atm%delp(is:ie,js:je,1:npz), Atm%coarse_graining%restart%delp)
    call weighted_block_average(Atm%gridstruct%area(is:ie,js:je), &
         Atm%phis(is:ie,js:je), Atm%coarse_graining%restart%phis)

    if (Atm%coarse_graining%write_coarse_agrid_vel_rst) then
       call weighted_block_average(mass(is:ie,js:je,1:npz), &
            Atm%ua(is:ie,js:je,1:npz), Atm%coarse_graining%restart%ua)
       call weighted_block_average(mass(is:ie,js:je,1:npz), &
            Atm%va(is:ie,js:je,1:npz), Atm%coarse_graining%restart%va)
    endif
  end subroutine coarse_grain_fv_core_restart_data_on_model_levels

  subroutine coarse_grain_fv_tracer_restart_data_on_model_levels(Atm, mass)
    type(fv_atmos_type), intent(inout) :: Atm
    real, intent(in) :: mass(is:ie,js:je,1:npz)

    character(len=64) :: tracer_name
    integer :: n_tracer

    do n_tracer = 1, n_prognostic_tracers
       call get_tracer_names(MODEL_ATMOS, n_tracer, tracer_name)
       if (trim(tracer_name) .eq. 'cld_amt') then
          call weighted_block_average(Atm%gridstruct%area(is:ie,js:je), &
               Atm%q(is:ie,js:je,1:npz,n_tracer), &
               Atm%coarse_graining%restart%q(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz,n_tracer))
       else
          call weighted_block_average(mass(is:ie,js:je,1:npz), &
               Atm%q(is:ie,js:je,1:npz,n_tracer), &
               Atm%coarse_graining%restart%q(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz,n_tracer))
       endif
    enddo

    do n_tracer = n_prognostic_tracers + 1, n_tracers
       call weighted_block_average(mass(is:ie,js:je,1:npz), &
               Atm%qdiag(is:ie,js:je,1:npz,n_tracer), &
               Atm%coarse_graining%restart%qdiag(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz,n_tracer))
    enddo
  end subroutine coarse_grain_fv_tracer_restart_data_on_model_levels

  subroutine coarse_grain_fv_srf_wnd_restart_data(Atm)
    type(fv_atmos_type), intent(inout) :: Atm

    call weighted_block_average(Atm%gridstruct%area(is:ie,js:je), &
         Atm%u_srf(is:ie,js:je), Atm%coarse_graining%restart%u_srf)
    call weighted_block_average(Atm%gridstruct%area(is:ie,js:je), &
         Atm%v_srf(is:ie,js:je), Atm%coarse_graining%restart%v_srf)
  end subroutine coarse_grain_fv_srf_wnd_restart_data

  subroutine coarse_grain_mg_drag_restart_data(Atm)
    type(fv_atmos_type), intent(inout) :: Atm

    call weighted_block_average(Atm%gridstruct%area(is:ie,js:je), &
         Atm%sgh(is:ie,js:je), Atm%coarse_graining%restart%sgh)
  end subroutine coarse_grain_mg_drag_restart_data

  subroutine coarse_grain_fv_land_restart_data(Atm)
    type(fv_atmos_type), intent(inout) :: Atm

    call weighted_block_average(Atm%gridstruct%area(is:ie,js:je), &
         Atm%oro(is:ie,js:je), Atm%coarse_graining%restart%oro)
  end subroutine coarse_grain_fv_land_restart_data

  subroutine coarse_grain_fv_core_restart_data_on_pressure_levels(&
     Atm, phalf, coarse_phalf, coarse_phalf_on_fine, masked_mass_weights, masked_area_weights)
     type(fv_atmos_type), intent(inout) :: Atm
     real, intent(in) :: phalf(is-1:ie+1,js-1:je+1,1:npz+1)
     real, intent(in) :: coarse_phalf(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz+1)
     real, intent(in) :: coarse_phalf_on_fine(is:ie,js:je,1:npz+1)
     real, intent(in), dimension(is:ie,js:je,1:npz) :: masked_mass_weights, masked_area_weights

     real, allocatable :: remapped(:,:,:)  ! Will re-use this to save memory

     allocate(remapped(is:ie,js:je,1:npz))

     if (Atm%coarse_graining%write_coarse_dgrid_vel_rst) then
        call remap_edges_along_x(Atm%u(is:ie,js:je+1,1:npz), &
             phalf(is-1:ie+1,js-1:je+1,1:npz+1), &
             Atm%gridstruct%dx(is:ie,js:je+1), &
             Atm%ptop, &
             Atm%coarse_graining%restart%u)
        call remap_edges_along_y(Atm%v(is:ie+1,js:je,1:npz), &
             phalf(is-1:ie+1,js-1:je+1,1:npz+1), &
             Atm%gridstruct%dy(is:ie+1,js:je), &
             Atm%ptop, &
             Atm%coarse_graining%restart%v)
     endif

     call vertically_remap_field(phalf(is:ie,js:je,1:npz+1), Atm%pt(is:ie,js:je,1:npz), coarse_phalf_on_fine, Atm%ptop, remapped)
     call weighted_block_average(masked_mass_weights, remapped, Atm%coarse_graining%restart%pt)

     if (.not. Atm%flagstruct%hydrostatic) then
       call vertically_remap_field(phalf(is:ie,js:je,1:npz+1), Atm%w(is:ie,js:je,1:npz), coarse_phalf_on_fine, Atm%ptop, remapped)
       call weighted_block_average(masked_mass_weights, remapped, Atm%coarse_graining%restart%w)
       call weighted_block_average(Atm%gridstruct%area(is:ie,js:je), Atm%delz(is:ie,js:je,1:npz), Atm%coarse_graining%restart%delz)
       if (Atm%flagstruct%hybrid_z) then
          call weighted_block_average(Atm%gridstruct%area(is:ie,js:je), Atm%ze0(is:ie,js:je,1:npz), Atm%coarse_graining%restart%ze0)
       endif
     endif

     call weighted_block_average(Atm%gridstruct%area(is:ie,js:je), Atm%phis(is:ie,js:je), Atm%coarse_graining%restart%phis)

     if (Atm%coarse_graining%write_coarse_agrid_vel_rst) then
       call vertically_remap_field(phalf(is:ie,js:je,1:npz+1), Atm%ua(is:ie,js:je,1:npz), coarse_phalf_on_fine, Atm%ptop, remapped)
       call weighted_block_average(masked_mass_weights, remapped, Atm%coarse_graining%restart%ua)
       call vertically_remap_field(phalf(is:ie,js:je,1:npz+1), Atm%va(is:ie,js:je,1:npz), coarse_phalf_on_fine, Atm%ptop, remapped)
       call weighted_block_average(masked_mass_weights, remapped, Atm%coarse_graining%restart%va)
     endif
  end subroutine coarse_grain_fv_core_restart_data_on_pressure_levels

  subroutine coarse_grain_fv_tracer_restart_data_on_pressure_levels( &
     Atm, phalf, coarse_phalf_on_fine, masked_mass_weights, masked_area_weights)
     type(fv_atmos_type), intent(inout) :: Atm
     real, intent(in) :: phalf(is-1:ie+1,js-1:je+1,1:npz+1)
     real, intent(in) :: coarse_phalf_on_fine(is:ie,js:je,1:npz+1)
     real, intent(in), dimension(is:ie,js:je,1:npz) :: masked_mass_weights, masked_area_weights

     real, allocatable :: remapped(:,:,:)
     character(len=64) :: tracer_name
     integer :: n_tracer

     allocate(remapped(is:ie,js:je,1:npz))

     do n_tracer = 1, n_prognostic_tracers
       call get_tracer_names(MODEL_ATMOS, n_tracer, tracer_name)
       call vertically_remap_field(phalf(is:ie,js:je,1:npz+1), &
         Atm%q(is:ie,js:je,1:npz,n_tracer), coarse_phalf_on_fine, Atm%ptop, remapped)
       if (trim(tracer_name) .eq. 'cld_amt') then
         call weighted_block_average(masked_area_weights, &
           remapped, &
           Atm%coarse_graining%restart%q(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz,n_tracer))
       else
         call weighted_block_average(masked_mass_weights, &
           remapped, &
           Atm%coarse_graining%restart%q(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz,n_tracer))
       endif
     enddo

     do n_tracer = n_prognostic_tracers + 1, n_tracers
       call vertically_remap_field(phalf(is:ie,js:je,1:npz+1), &
         Atm%qdiag(is:ie,js:je,1:npz,n_tracer), coarse_phalf_on_fine, Atm%ptop, remapped)
       call weighted_block_average(masked_mass_weights, &
         remapped, &
         Atm%coarse_graining%restart%qdiag(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz,n_tracer))
     enddo
   end subroutine coarse_grain_fv_tracer_restart_data_on_pressure_levels

   subroutine compute_top_height(delz, phis, top_height)
     real, intent(in) :: delz(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz)
     real, intent(in) :: phis(is_coarse:ie_coarse,js_coarse:je_coarse)
     real, intent(out) :: top_height(is_coarse:ie_coarse,js_coarse:je_coarse)

     top_height = (phis / GRAV) - sum(delz, dim=3)
   end subroutine compute_top_height

   subroutine hydrostatic_delz(phalf, temp, sphum, delz)
     real, intent(in) :: phalf(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz+1)
     real, intent(in) :: temp(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz)
     real, intent(in) :: sphum(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz)
     real, intent(out) :: delz(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz)

     real, allocatable :: virtual_temp(:,:,:), dlogp(:,:,:)
     integer :: k

     allocate(virtual_temp(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz))
     allocate(dlogp(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz))

     virtual_temp = temp * (1.0 + (RVGAS / RDGAS - 1.0) * sphum)
     do k = 1, npz
        dlogp(:,:,k) = log(phalf(:,:,k+1)) - log(phalf(:,:,k))
     enddo
     delz = -dlogp * RDGAS * virtual_temp / GRAV
   end subroutine hydrostatic_delz

   subroutine delz_and_top_height_to_phis(top_height, delz, phis)
     real, intent(in) :: top_height(is_coarse:ie_coarse,js_coarse:je_coarse)
     real, intent(in) :: delz(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz)
     real, intent(out) :: phis(is_coarse:ie_coarse,js_coarse:je_coarse)

     phis = GRAV * (top_height + sum(delz, dim=3))
   end subroutine delz_and_top_height_to_phis

   subroutine impose_hydrostatic_balance(Atm, coarse_phalf)
     type(fv_atmos_type), intent(inout) :: Atm
     real, intent(in) :: coarse_phalf(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz+1)

     integer :: sphum
     real, allocatable :: top_height(:,:)
     allocate(top_height(is_coarse:ie_coarse,js_coarse:je_coarse))

     sphum = get_tracer_index(MODEL_ATMOS, 'sphum')

     call compute_top_height(Atm%coarse_graining%restart%delz, Atm%coarse_graining%restart%phis, top_height)
     call hydrostatic_delz(coarse_phalf, Atm%coarse_graining%restart%pt, Atm%coarse_graining%restart%q(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz,sphum), Atm%coarse_graining%restart%delz)
     call delz_and_top_height_to_phis(top_height, Atm%coarse_graining%restart%delz, Atm%coarse_graining%restart%phis)
   end subroutine impose_hydrostatic_balance

  subroutine compute_pressure_level_coarse_graining_requirements( &
     Atm, phalf, coarse_phalf, coarse_phalf_on_fine, masked_mass_weights, masked_area_weights)
     type(fv_atmos_type), intent(inout) :: Atm
     real, intent(out) :: phalf(is-1:ie+1,js-1:je+1,1:npz+1)
     real, intent(out) :: coarse_phalf(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz+1)
     real, intent(out) :: coarse_phalf_on_fine(is:ie,js:je,1:npz+1)
     real, intent(out), dimension(is:ie,js:je,1:npz) :: masked_mass_weights, masked_area_weights

     ! Do a halo update on delp before proceeding here, because the remapping procedure
     ! for the winds requires interpolating across tile edges.
     call mpp_update_domains(Atm%delp, Atm%domain, complete=.true.)
     call compute_phalf(is-1, ie+1, js-1, je+1, Atm%delp(is-1:ie+1,js-1:je+1,1:npz), Atm%ptop, phalf)
     call weighted_block_average(Atm%gridstruct%area(is:ie,js:je), Atm%delp(is:ie,js:je,1:npz), Atm%coarse_graining%restart%delp)
     call compute_phalf(is_coarse, ie_coarse, js_coarse, je_coarse, Atm%coarse_graining%restart%delp, Atm%ptop, coarse_phalf)
     call block_upsample(coarse_phalf, coarse_phalf_on_fine, npz+1)
     call mask_mass_weights(Atm%gridstruct%area(is:ie,js:je), Atm%delp(is:ie,js:je,1:npz), phalf(is:ie,js:je,1:npz+1), coarse_phalf_on_fine, masked_mass_weights)
     call mask_area_weights(Atm%gridstruct%area(is:ie,js:je), phalf(is:ie,js:je,1:npz+1), coarse_phalf_on_fine, masked_area_weights)
  end subroutine compute_pressure_level_coarse_graining_requirements

  subroutine compute_phalf(i_start, i_end, j_start, j_end, delp, ptop, phalf)
     integer, intent(in) :: i_start, i_end, j_start, j_end
     real, intent(in) :: delp(i_start:i_end,j_start:j_end,1:npz)
     real, intent(in) :: ptop
     real, intent(out) :: phalf(i_start:i_end,j_start:j_end,1:npz+1)

     integer :: k

     phalf(:,:,1) = ptop
     do k = 2, npz + 1
       phalf(:,:,k) = phalf(:,:,k-1) + delp(:,:,k-1)
     enddo
  end subroutine compute_phalf
end module coarse_grained_restart_files_mod
