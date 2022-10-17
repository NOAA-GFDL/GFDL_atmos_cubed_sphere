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
!* useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

module coarse_grained_restart_files_mod

  use coarse_graining_mod, only: compute_mass_weights, get_coarse_array_bounds,&
       get_fine_array_bounds, MODEL_LEVEL, PRESSURE_LEVEL, weighted_block_average, &
       weighted_block_edge_average_x, weighted_block_edge_average_y, &
       mask_area_weights, block_upsample, remap_edges_along_x, &
       remap_edges_along_y, vertically_remap_field
  use constants_mod, only: GRAV, RDGAS, RVGAS
  use field_manager_mod, only: MODEL_ATMOS
  use fms2_io_mod,      only: register_restart_field, write_restart, open_file, close_file, register_variable_attribute, variable_exists
  use fv_arrays_mod, only: coarse_restart_type, fv_atmos_type
  use mpp_domains_mod, only: domain2d, EAST, NORTH, CENTER, mpp_update_domains
  use mpp_mod, only: FATAL, mpp_error
  use tracer_manager_mod, only: get_tracer_names, get_tracer_index, set_tracer_profile
  use fv_io_mod, only: fv_io_register_axis

  implicit none
  private

  public :: deallocate_coarse_restart_type, fv_coarse_restart_init, fv_io_write_restart_coarse

  ! Global variables for this module, initialized in fv_coarse_restart_init
  integer :: is, ie, js, je, npz
  integer :: is_coarse, ie_coarse, js_coarse, je_coarse
  integer :: n_prognostic_tracers, n_diagnostic_tracers, n_tracers

contains

  subroutine fv_coarse_restart_init(nz, nt_prog, &
       nt_phys, hydrostatic, hybrid_z, fv_land, &
       write_coarse_dgrid_vel_rst, write_coarse_agrid_vel_rst, &
       restart)
    integer, intent(in) :: nz, nt_prog, nt_phys
    logical, intent(in) :: hydrostatic, hybrid_z, fv_land
    logical, intent(in) :: write_coarse_dgrid_vel_rst, write_coarse_agrid_vel_rst
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
  end subroutine fv_coarse_restart_init

  subroutine fv_io_write_restart_coarse(Atm, timestamp)
    type(fv_atmos_type), intent(inout) :: Atm
    character(len=*), optional, intent(in) :: timestamp

    integer :: tile_count, n_tiles

    call register_coarse_restart_files(Atm%flagstruct%hydrostatic, &
           Atm%flagstruct%hybrid_z, Atm%flagstruct%fv_land, &
           Atm%coarse_graining%write_coarse_dgrid_vel_rst, &
           Atm%coarse_graining%write_coarse_agrid_vel_rst, &
           Atm%coarse_graining%domain, &
           Atm%coarse_graining%restart, timestamp)

    call coarse_grain_restart_data(Atm)

    if (Atm%coarse_graining%restart%fv_core_coarse_is_open) then
      call write_restart(Atm%coarse_graining%restart%fv_core_coarse)
      call close_file(Atm%coarse_graining%restart%fv_core_coarse)
      Atm%coarse_graining%restart%fv_core_coarse_is_open=.false.
    endif
    if (Atm%coarse_graining%restart%fv_tracer_coarse_is_open) then
      call write_restart(Atm%coarse_graining%restart%fv_tracer_coarse)
      call close_file(Atm%coarse_graining%restart%fv_tracer_coarse)
      Atm%coarse_graining%restart%fv_tracer_coarse_is_open=.false.
    endif
    if (Atm%coarse_graining%restart%fv_srf_wnd_coarse_is_open) then
      call write_restart(Atm%coarse_graining%restart%fv_srf_wnd_coarse)
      call close_file(Atm%coarse_graining%restart%fv_srf_wnd_coarse)
      Atm%coarse_graining%restart%fv_srf_wnd_coarse_is_open=.false.
    endif
    if (Atm%flagstruct%fv_land) then
       if (Atm%coarse_graining%restart%mg_drag_coarse_is_open) then
         call write_restart(Atm%coarse_graining%restart%mg_drag_coarse)
         call close_file(Atm%coarse_graining%restart%mg_drag_coarse)
         Atm%coarse_graining%restart%mg_drag_coarse_is_open=.false.
       endif
       if (Atm%coarse_graining%restart%fv_land_coarse_is_open) then
         call write_restart(Atm%coarse_graining%restart%fv_land_coarse)
         call close_file(Atm%coarse_graining%restart%fv_land_coarse)
         Atm%coarse_graining%restart%fv_land_coarse_is_open=.false.
       endif
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

  subroutine register_coarse_restart_files(hydrostatic, &
       hybrid_z, fv_land, write_coarse_dgrid_vel_rst, write_coarse_agrid_vel_rst, &
       coarse_domain, restart, timestamp)
    logical, intent(in) :: hydrostatic, hybrid_z, fv_land
    logical, intent(in) :: write_coarse_dgrid_vel_rst, write_coarse_agrid_vel_rst
    type(domain2d), intent(in) :: coarse_domain
    type(coarse_restart_type), intent(inout) :: restart
    character(len=*), optional, intent(in) :: timestamp

    call register_fv_core_coarse(hydrostatic, hybrid_z, &
         write_coarse_dgrid_vel_rst, write_coarse_agrid_vel_rst, &
         coarse_domain, restart, timestamp)
    call register_fv_tracer_coarse(coarse_domain, restart, timestamp)
    call register_fv_srf_wnd_coarse(coarse_domain, restart, timestamp)
    if (fv_land) then
       call register_mg_drag_coarse(coarse_domain, restart, timestamp)
       call register_fv_land_coarse(coarse_domain, restart, timestamp)
    endif
  end subroutine register_coarse_restart_files

  subroutine register_fv_core_coarse(hydrostatic, hybrid_z, &
       write_coarse_dgrid_vel_rst, write_coarse_agrid_vel_rst, coarse_domain, &
       restart, timestamp)
    logical, intent(in) :: hydrostatic, hybrid_z
    logical, intent(in) :: write_coarse_dgrid_vel_rst, write_coarse_agrid_vel_rst
    type(domain2d), intent(in) :: coarse_domain
    type(coarse_restart_type), intent(inout) :: restart
    character(len=*), optional, intent(in) :: timestamp

    character(len=8), dimension(4) :: dim_names_4d, dim_names_4d2, dim_names_4d3
    character(len=8), dimension(3) :: dim_names_3d
    character(len=64) :: filename
    integer, dimension(1) :: zsize

    dim_names_4d(1) = "xaxis_1"
    dim_names_4d(2) = "yaxis_1"
    dim_names_4d(3) = "zaxis_1"
    dim_names_4d(4) = "Time"
    dim_names_4d2 = dim_names_4d
    dim_names_4d2(1) = "xaxis_2"
    dim_names_4d2(2) = "yaxis_2"
    dim_names_4d3 = dim_names_4d
    dim_names_4d3(2) = "yaxis_2"
    dim_names_3d(1) = "xaxis_1"
    dim_names_3d(2) = "yaxis_2"
    dim_names_3d(3) = "Time"

    if (present(timestamp)) then
      filename = 'RESTART/'//trim(timestamp)//'.fv_core_coarse.res.nc'
    else
      filename = 'RESTART/fv_core_coarse.res.nc'
    endif

    restart%fv_core_coarse_is_open = open_file(restart%fv_core_coarse, filename, &
            "overwrite", coarse_domain, is_restart=.true.)
    if (restart%fv_core_coarse_is_open) then
      zsize = (/size(restart%u,3)/)
      call fv_io_register_axis(restart%fv_core_coarse, numx=2 ,numy=2, xpos=(/CENTER, EAST/) ,ypos=(/NORTH, CENTER/) ,numz=1, zsize=zsize)

      if (write_coarse_dgrid_vel_rst) then
         call register_restart_field(restart%fv_core_coarse, &
              'u', restart%u, dim_names_4d)
         call register_variable_attribute(restart%fv_core_coarse, &
              'u', "long_name", "u", str_len=len("u"))
         call register_variable_attribute(restart%fv_core_coarse, &
              'u', "units", "none", str_len=len("none"))
         call register_restart_field(restart%fv_core_coarse, &
              'v', restart%v, dim_names_4d2)
         call register_variable_attribute(restart%fv_core_coarse, &
              'v', "long_name", "v", str_len=len("v"))
         call register_variable_attribute(restart%fv_core_coarse, &
              'v', "units", "none", str_len=len("none"))
      endif

      if (write_coarse_agrid_vel_rst) then
         call register_restart_field(restart%fv_core_coarse, &
              'ua', restart%ua, dim_names_4d3)
         call register_variable_attribute(restart%fv_core_coarse, &
              'ua', "long_name", "ua", str_len=len("ua"))
         call register_variable_attribute(restart%fv_core_coarse, &
              'ua', "units", "none", str_len=len("none"))
         call register_restart_field(restart%fv_core_coarse, &
              'va', restart%va, dim_names_4d3)
         call register_variable_attribute(restart%fv_core_coarse, &
              'va', "long_name", "va", str_len=len("va"))
         call register_variable_attribute(restart%fv_core_coarse, &
              'va', "units", "none", str_len=len("none"))
      endif

      if (.not. hydrostatic) then
         call register_restart_field(restart%fv_core_coarse, &
              'W', restart%w, dim_names_4d3, is_optional=.true.)
         if (variable_exists(restart%fv_core_coarse, 'W')) then
            call register_variable_attribute(restart%fv_core_coarse, &
                 'W', "long_name", "W", str_len=len("W"))
            call register_variable_attribute(restart%fv_core_coarse, &
                 'W', "units", "none", str_len=len("none"))
         endif
         call register_restart_field(restart%fv_core_coarse, &
              'DZ', restart%delz, dim_names_4d3, is_optional=.true.)
         if (variable_exists(restart%fv_core_coarse, 'DZ')) then
            call register_variable_attribute(restart%fv_core_coarse, &
                 'DZ', "long_name", "DZ", str_len=len("DZ"))
            call register_variable_attribute(restart%fv_core_coarse, &
                 'DZ', "units", "none", str_len=len("none"))
         endif
         if (hybrid_z) then
            call register_restart_field(restart%fv_core_coarse, &
                 'ZE0', restart%ze0, dim_names_4d3, is_optional=.false.)
            call register_variable_attribute(restart%fv_core_coarse, &
                 'ZE0', "long_name", "ZE0", str_len=len("ZE0"))
            call register_variable_attribute(restart%fv_core_coarse, &
                 'ZE0', "units", "none", str_len=len("none"))
         endif
      endif

      call register_restart_field(restart%fv_core_coarse, &
           'T', restart%pt, dim_names_4d3)
      call register_variable_attribute(restart%fv_core_coarse, &
           'T', "long_name", "T", str_len=len("T"))
      call register_variable_attribute(restart%fv_core_coarse, &
           'T', "units", "none", str_len=len("none"))
      call register_restart_field(restart%fv_core_coarse, &
           'delp', restart%delp, dim_names_4d3)
      call register_variable_attribute(restart%fv_core_coarse, &
           'delp', "long_name", "delp", str_len=len("delp"))
      call register_variable_attribute(restart%fv_core_coarse, &
           'delp', "units", "none", str_len=len("none"))
      call register_restart_field(restart%fv_core_coarse, &
           'phis', restart%phis, dim_names_3d)
      call register_variable_attribute(restart%fv_core_coarse, &
           'phis', "long_name", "phis", str_len=len("phis"))
      call register_variable_attribute(restart%fv_core_coarse, &
           'phis', "units", "none", str_len=len("none"))
    endif
  end subroutine register_fv_core_coarse

  subroutine register_fv_tracer_coarse(coarse_domain, restart, timestamp)
    type(domain2d), intent(in) :: coarse_domain
    type(coarse_restart_type), intent(inout) :: restart
    character(len=*), optional, intent(in) :: timestamp

    character(len=8), dimension(4) :: dim_names_4d
    character(len=64) :: filename, tracer_name
    integer :: n_tracer
    integer, dimension(1) :: zsize

    dim_names_4d(1) = "xaxis_1"
    dim_names_4d(2) = "yaxis_1"
    dim_names_4d(3) = "zaxis_1"
    dim_names_4d(4) = "Time"

    if (present(timestamp)) then
      filename = 'RESTART/'//trim(timestamp)//'.fv_tracer_coarse.res.nc'
    else
      filename = 'RESTART/fv_tracer_coarse.res.nc'
    endif

    restart%fv_tracer_coarse_is_open = open_file(restart%fv_tracer_coarse, filename, &
            "overwrite", coarse_domain, is_restart=.true.)
    if (restart%fv_tracer_coarse_is_open) then
      zsize=(/size(restart%q,3)/)
      call fv_io_register_axis(restart%fv_tracer_coarse, numx=1 ,numy=1, xpos=(/CENTER/) ,ypos=(/CENTER/) ,numz=1, zsize=zsize)

      do n_tracer = 1, n_prognostic_tracers
         call get_tracer_names(MODEL_ATMOS, n_tracer, tracer_name)
         call set_tracer_profile(MODEL_ATMOS, n_tracer, restart%q(:,:,:,n_tracer))
         call register_restart_field(restart%fv_tracer_coarse, &
              tracer_name, restart%q(:,:,:,n_tracer), dim_names_4d, &
              is_optional=.true.)
         if (variable_exists(restart%fv_tracer_coarse, tracer_name)) then
            call register_variable_attribute(restart%fv_tracer_coarse, &
                 tracer_name, "long_name", tracer_name, str_len=len(tracer_name))
            call register_variable_attribute(restart%fv_tracer_coarse, &
                 tracer_name, "units", "none", str_len=len("none"))
         endif
      enddo

      do n_tracer = n_prognostic_tracers + 1, n_tracers
         call get_tracer_names(MODEL_ATMOS, n_tracer, tracer_name)
         call set_tracer_profile(MODEL_ATMOS, n_tracer, restart%qdiag(:,:,:,n_tracer))
         call register_restart_field(restart%fv_tracer_coarse, &
              tracer_name, restart%qdiag(:,:,:,n_tracer), dim_names_4d, &
              is_optional=.true.)
         if (variable_exists(restart%fv_tracer_coarse, tracer_name)) then
            call register_variable_attribute(restart%fv_tracer_coarse, &
                 tracer_name, "long_name", tracer_name, str_len=len(tracer_name))
            call register_variable_attribute(restart%fv_tracer_coarse, &
                 tracer_name, "units", "none", str_len=len("none"))
         endif
      enddo
    endif
  end subroutine register_fv_tracer_coarse

  subroutine register_fv_srf_wnd_coarse(coarse_domain, restart, timestamp)
    type(domain2d), intent(in) :: coarse_domain
    type(coarse_restart_type), intent(inout) :: restart
    character(len=*), optional, intent(in) :: timestamp

    character(len=8), dimension(3) :: dim_names_3d
    character(len=64) :: filename

    dim_names_3d(1) = "xaxis_1"
    dim_names_3d(2) = "yaxis_1"
    dim_names_3d(3) = "Time"

    if (present(timestamp)) then
      filename = 'RESTART/'//trim(timestamp)//'.fv_srf_wnd_coarse.res.nc'
    else
      filename = 'RESTART/fv_srf_wnd_coarse.res.nc'
    endif

    restart%fv_srf_wnd_coarse_is_open = open_file(restart%fv_srf_wnd_coarse, filename, &
            "overwrite", coarse_domain, is_restart=.true.)
    if (restart%fv_srf_wnd_coarse_is_open) then
      call fv_io_register_axis(restart%fv_srf_wnd_coarse, numx=1 ,numy=1, xpos=(/CENTER/) ,ypos=(/CENTER/))

      call register_restart_field(restart%fv_srf_wnd_coarse, &
           'u_srf', restart%u_srf, dim_names_3d)
      call register_variable_attribute(restart%fv_srf_wnd_coarse, &
           'u_srf', "long_name", "u_srf", str_len=len("u_srf"))
      call register_variable_attribute(restart%fv_srf_wnd_coarse, &
           'u_srf', "units", "none", str_len=len("none"))
      call register_restart_field(restart%fv_srf_wnd_coarse, &
           'v_srf', restart%v_srf, dim_names_3d)
      call register_variable_attribute(restart%fv_srf_wnd_coarse, &
           'v_srf', "long_name", "v_srf", str_len=len("v_srf"))
      call register_variable_attribute(restart%fv_srf_wnd_coarse, &
           'v_srf', "units", "none", str_len=len("none"))
    endif
  end subroutine register_fv_srf_wnd_coarse

  subroutine register_mg_drag_coarse(coarse_domain, restart, timestamp)
    type(domain2d), intent(in) :: coarse_domain
    type(coarse_restart_type), intent(out) :: restart
    character(len=*), optional, intent(in) :: timestamp

    character(len=8), dimension(3) :: dim_names_3d
    character(len=64) :: filename

    dim_names_3d(1) = "xaxis_1"
    dim_names_3d(2) = "yaxis_1"
    dim_names_3d(3) = "Time"

    if (present(timestamp)) then
      filename = 'RESTART/'//trim(timestamp)//'.mg_drag_coarse.res.nc'
    else
      filename = 'RESTART/mg_drag_coarse.res.nc'
    endif

    restart%mg_drag_coarse_is_open = open_file(restart%mg_drag_coarse, filename, &
            "overwrite", coarse_domain, is_restart=.true.)
    if (restart%mg_drag_coarse_is_open) then
      call fv_io_register_axis(restart%mg_drag_coarse, numx=1, numy=1, xpos=(/CENTER/), ypos=(/CENTER/))

      call register_restart_field(restart%mg_drag_coarse, &
          'ghprime', restart%sgh, dim_names_3d)
      call register_variable_attribute(restart%mg_drag_coarse, &
          'ghprime', "long_name", "ghprime", str_len=len("ghprime"))
      call register_variable_attribute(restart%mg_drag_coarse, &
          'ghprime', "units", "none", str_len=len("none"))
    endif
  end subroutine register_mg_drag_coarse

  subroutine register_fv_land_coarse(coarse_domain, restart, timestamp)
    type(domain2d), intent(in) :: coarse_domain
    type(coarse_restart_type), intent(inout) :: restart
    character(len=*), optional, intent(in) :: timestamp

    character(len=8), dimension(3) :: dim_names_3d
    character(len=64) :: filename

    dim_names_3d(1) = "xaxis_1"
    dim_names_3d(2) = "yaxis_1"
    dim_names_3d(3) = "Time"

    if (present(timestamp)) then
      filename = 'RESTART/'//trim(timestamp)//'.fv_land_coarse.res.nc'
    else
      filename = 'RESTART/fv_land_coarse.res.nc'
    endif

    restart%fv_land_coarse_is_open = open_file(restart%fv_land_coarse, filename, &
            "overwrite", coarse_domain, is_restart=.true.)
    if (restart%fv_land_coarse_is_open) then
      call fv_io_register_axis(restart%fv_land_coarse, numx=1, numy=1, xpos=(/CENTER/), ypos=(/CENTER/))

      call register_restart_field(restart%fv_land_coarse, &
          'oro', restart%oro, dim_names_3d)
      call register_variable_attribute(restart%fv_land_coarse, &
           'oro', "long_name", "oro", str_len=len("oro"))
      call register_variable_attribute(restart%fv_core_coarse, &
           'oro', "units", "none", str_len=len("none"))
    endif
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
     real, allocatable, dimension(:,:,:) :: masked_area_weights

     allocate(phalf(is-1:ie+1,js-1:je+1,1:npz+1))  ! Require the halo here for the winds
     allocate(coarse_phalf(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz+1))
     allocate(coarse_phalf_on_fine(is:ie,js:je,1:npz+1))
     allocate(masked_area_weights(is:ie,js:je,1:npz))

     ! delp and delz are coarse-grained on model levels; u, v, W, T, and all the tracers
     ! are all remapped to surfaces of constant pressure within coarse grid cells before
     ! coarse graining.  At the end, delz and phis are corrected to impose hydrostatic balance.
     call compute_pressure_level_coarse_graining_requirements( &
       Atm, phalf, coarse_phalf, coarse_phalf_on_fine, masked_area_weights)
     call coarse_grain_fv_core_restart_data_on_pressure_levels( &
       Atm, phalf, coarse_phalf, coarse_phalf_on_fine, masked_area_weights)
     call coarse_grain_fv_tracer_restart_data_on_pressure_levels( &
       Atm, phalf, coarse_phalf_on_fine, masked_area_weights)
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
     Atm, phalf, coarse_phalf, coarse_phalf_on_fine, masked_area_weights)
     type(fv_atmos_type), intent(inout) :: Atm
     real, intent(in) :: phalf(is-1:ie+1,js-1:je+1,1:npz+1)
     real, intent(in) :: coarse_phalf(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz+1)
     real, intent(in) :: coarse_phalf_on_fine(is:ie,js:je,1:npz+1)
     real, intent(in), dimension(is:ie,js:je,1:npz) :: masked_area_weights

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
     call weighted_block_average(masked_area_weights, remapped, Atm%coarse_graining%restart%pt)

     if (.not. Atm%flagstruct%hydrostatic) then
       call vertically_remap_field(phalf(is:ie,js:je,1:npz+1), Atm%w(is:ie,js:je,1:npz), coarse_phalf_on_fine, Atm%ptop, remapped)
       call weighted_block_average(masked_area_weights, remapped, Atm%coarse_graining%restart%w)
       call weighted_block_average(Atm%gridstruct%area(is:ie,js:je), Atm%delz(is:ie,js:je,1:npz), Atm%coarse_graining%restart%delz)
       if (Atm%flagstruct%hybrid_z) then
          call weighted_block_average(Atm%gridstruct%area(is:ie,js:je), Atm%ze0(is:ie,js:je,1:npz), Atm%coarse_graining%restart%ze0)
       endif
     endif

     call weighted_block_average(Atm%gridstruct%area(is:ie,js:je), Atm%phis(is:ie,js:je), Atm%coarse_graining%restart%phis)

     if (Atm%coarse_graining%write_coarse_agrid_vel_rst) then
       call vertically_remap_field(phalf(is:ie,js:je,1:npz+1), Atm%ua(is:ie,js:je,1:npz), coarse_phalf_on_fine, Atm%ptop, remapped)
       call weighted_block_average(masked_area_weights, remapped, Atm%coarse_graining%restart%ua)
       call vertically_remap_field(phalf(is:ie,js:je,1:npz+1), Atm%va(is:ie,js:je,1:npz), coarse_phalf_on_fine, Atm%ptop, remapped)
       call weighted_block_average(masked_area_weights, remapped, Atm%coarse_graining%restart%va)
     endif
  end subroutine coarse_grain_fv_core_restart_data_on_pressure_levels

  subroutine coarse_grain_fv_tracer_restart_data_on_pressure_levels( &
     Atm, phalf, coarse_phalf_on_fine, masked_area_weights)
     type(fv_atmos_type), intent(inout) :: Atm
     real, intent(in) :: phalf(is-1:ie+1,js-1:je+1,1:npz+1)
     real, intent(in) :: coarse_phalf_on_fine(is:ie,js:je,1:npz+1)
     real, intent(in), dimension(is:ie,js:je,1:npz) :: masked_area_weights

     real, allocatable :: remapped(:,:,:)
     integer :: n_tracer

     allocate(remapped(is:ie,js:je,1:npz))

     do n_tracer = 1, n_prognostic_tracers
       call vertically_remap_field(phalf(is:ie,js:je,1:npz+1), &
         Atm%q(is:ie,js:je,1:npz,n_tracer), coarse_phalf_on_fine, Atm%ptop, remapped)
       call weighted_block_average(masked_area_weights, &
         remapped, &
         Atm%coarse_graining%restart%q(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz,n_tracer))
     enddo

     do n_tracer = n_prognostic_tracers + 1, n_tracers
       call vertically_remap_field(phalf(is:ie,js:je,1:npz+1), &
         Atm%qdiag(is:ie,js:je,1:npz,n_tracer), coarse_phalf_on_fine, Atm%ptop, remapped)
       call weighted_block_average(masked_area_weights, &
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
     Atm, phalf, coarse_phalf, coarse_phalf_on_fine, masked_area_weights)
     type(fv_atmos_type), intent(inout) :: Atm
     real, intent(out) :: phalf(is-1:ie+1,js-1:je+1,1:npz+1)
     real, intent(out) :: coarse_phalf(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz+1)
     real, intent(out) :: coarse_phalf_on_fine(is:ie,js:je,1:npz+1)
     real, intent(out), dimension(is:ie,js:je,1:npz) :: masked_area_weights

     ! Do a halo update on delp before proceeding here, because the remapping procedure
     ! for the winds requires interpolating across tile edges.
     call mpp_update_domains(Atm%delp, Atm%domain, complete=.true.)
     call compute_phalf(is-1, ie+1, js-1, je+1, Atm%delp(is-1:ie+1,js-1:je+1,1:npz), Atm%ptop, phalf)
     call weighted_block_average(Atm%gridstruct%area(is:ie,js:je), Atm%delp(is:ie,js:je,1:npz), Atm%coarse_graining%restart%delp)
     call compute_phalf(is_coarse, ie_coarse, js_coarse, je_coarse, Atm%coarse_graining%restart%delp, Atm%ptop, coarse_phalf)
     call block_upsample(coarse_phalf, coarse_phalf_on_fine, npz+1)
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
