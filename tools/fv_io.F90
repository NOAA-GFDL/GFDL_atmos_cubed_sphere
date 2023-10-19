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

!!!NOTE: Merging in the seasonal forecast initialization code
!!!!     has proven problematic in the past, since many conflicts
!!!!     occur. Leaving this for now --- lmh 10aug15

module fv_io_mod

  !<OVERVIEW>
  ! Restart facilities for FV core
  !</OVERVIEW>
  !<DESCRIPTION>
  ! This module writes and reads restart files for the FV core. Additionally
  ! it provides setup and calls routines necessary to provide a complete restart
  ! for the model.
  !</DESCRIPTION>

  use fms2_io_mod,             only: FmsNetcdfFile_t, FmsNetcdfDomainFile_t, &
                                     register_restart_field, register_axis, unlimited, &
                                     open_file, read_restart, read_restart_bc, write_restart, &
                                     write_restart_bc, close_file, register_field, write_data, &
                                     get_global_io_domain_indices, register_variable_attribute, &
                                     variable_exists, read_data, set_filename_appendix
  use mpp_mod,                 only: mpp_error, FATAL, NOTE, WARNING, mpp_root_pe, &
                                     mpp_sync, mpp_pe, mpp_declare_pelist, mpp_get_current_pelist, &
                                     mpp_npes
  use mpp_domains_mod,         only: domain2d, EAST, WEST, NORTH, CENTER, SOUTH, CORNER, &
                                     mpp_get_compute_domain, mpp_get_data_domain, &
                                     mpp_get_layout, mpp_get_ntile_count, &
                                     mpp_get_global_domain, mpp_update_domains
  use tracer_manager_mod,      only: tr_get_tracer_names=>get_tracer_names, &
                                     get_tracer_names, get_number_tracers, &
                                     set_tracer_profile, &
                                     get_tracer_index
  use field_manager_mod,       only: MODEL_ATMOS
  use external_sst_mod,        only: sst_ncep, sst_anom, use_ncep_sst
  use fv_arrays_mod,           only: fv_atmos_type, fv_nest_BC_type_3D, R_GRID, &
                                     fv_grid_bounds_type, fv_grid_type
  use fv_eta_mod,              only: set_external_eta

  use fv_mp_mod,               only: mp_gather, is_master
  use fv_treat_da_inc_mod,     only: read_da_inc
  use mpp_parameter_mod,       only: DGRID_NE
  use fv_grid_utils_mod,       only: cubed_a2d
  
  implicit none
  private

  public :: fv_io_init, fv_io_exit, fv_io_read_restart, remap_restart, fv_io_write_restart
  public :: fv_io_read_tracers, fv_io_register_restart, fv_io_register_nudge_restart
  public :: fv_io_register_restart_BCs
  public :: fv_io_write_BCs, fv_io_read_BCs
  public :: fv_io_register_axis

  logical                       :: module_is_initialized = .FALSE.


  integer ::grid_xtdimid, grid_ytdimid, haloid, pfullid !For writing BCs
  integer ::grid_xtstagdimid, grid_ytstagdimid, oneid

#ifdef OVERLOAD_R4
  character(len=5), parameter :: axis_type = 'float'
#else
  character(len=6), parameter :: axis_type = 'double'
#endif

contains

  !#####################################################################
  ! <SUBROUTINE NAME="fv_io_init">
  !
  ! <DESCRIPTION>
  ! Initialize the fv core restart facilities
  ! </DESCRIPTION>
  !
  subroutine fv_io_init()
    module_is_initialized = .TRUE.
  end subroutine fv_io_init
  ! </SUBROUTINE> NAME="fv_io_init"


  !#####################################################################
  ! <SUBROUTINE NAME="fv_io_exit">
  !
  ! <DESCRIPTION>
  ! Close the fv core restart facilities
  ! </DESCRIPTION>
  !
  subroutine fv_io_exit
    module_is_initialized = .FALSE.
  end subroutine fv_io_exit
  ! </SUBROUTINE> NAME="fv_io_exit"


  !#####################################################################
  ! <SUBROUTINE NAME="fv_io_register_axis">
  !
  ! <DESCRIPTION>
  ! Register the fv axis for new fms2 io
  ! </DESCRIPTION>
  !
  subroutine fv_io_register_axis(file_obj, numx, xpos, numy, ypos, numz, zsize)
    type(FmsNetcdfDomainFile_t), intent(inout) ::  file_obj
    integer, intent(in), optional :: numx, numy, numz
    integer, dimension(:), intent(in), optional :: xpos, ypos, zsize

    integer :: i, ie, is, j
    integer, dimension(:), allocatable :: buffer
    character(len=1) :: suffix
    character(len=7) :: axisname

    if (present(numx)) then
    do i=1, numx
      write(suffix,'(I1)') i
      axisname = 'xaxis_'//suffix
      call register_axis(file_obj, axisname, 'X', domain_position=xpos(i))
      if (.not. file_obj%is_readonly) then !if writing file
        call register_field(file_obj, axisname, axis_type, (/axisname/))
        call register_variable_attribute(file_obj,axisname, "long_name", axisname, str_len=len(axisname))
        call register_variable_attribute(file_obj,axisname, "units", "none", str_len=len("none"))
        call register_variable_attribute(file_obj,axisname, "cartesian_axis", "X", str_len=1)
        call get_global_io_domain_indices(file_obj, axisname, is, ie, buffer)
        call write_data(file_obj, axisname, buffer)
      endif
    end do
    endif

    if (present(numy)) then
    do i=1, numy
      write(suffix,'(I1)') i
      axisname = 'yaxis_'//suffix
      call register_axis(file_obj, axisname, 'Y', domain_position=ypos(i))
      if (.not. file_obj%is_readonly) then !if writing file
        call register_field(file_obj, axisname, axis_type, (/axisname/))
        call register_variable_attribute(file_obj,axisname, "long_name", axisname, str_len=len(axisname))
        call register_variable_attribute(file_obj,axisname, "units", "none", str_len=len("none"))
        call register_variable_attribute(file_obj,axisname, "cartesian_axis", "Y", str_len=1)
        call get_global_io_domain_indices(file_obj, axisname, is, ie, buffer)
        call write_data(file_obj, axisname, buffer)
      endif
    end do
    endif

    if (present(numz)) then
      do i= 1, numz
        write(suffix,'(I1)') i
        axisname = 'zaxis_'//suffix
        call register_axis(file_obj, axisname, zsize(i))
        if (.not. file_obj%is_readonly) then !if writing file
          call register_field(file_obj, axisname, axis_type, (/axisname/))
          call register_variable_attribute(file_obj,axisname, "long_name", axisname, str_len=len(axisname))
          call register_variable_attribute(file_obj,axisname, "units", "none", str_len=len("none"))
          call register_variable_attribute(file_obj,axisname, "cartesian_axis", "Z", str_len=1)
          if (allocated(buffer)) deallocate(buffer)
          allocate(buffer(zsize(i)))
          do j = 1, zsize(i)
            buffer(j) = j
          end do
          call write_data(file_obj, axisname, buffer)
          deallocate(buffer)
        endif
      end do
    endif

    call register_axis(file_obj, "Time", unlimited)
    if (.not. file_obj%is_readonly) then !if writing file
       call register_field(file_obj, "Time", axis_type, (/"Time"/))
       call register_variable_attribute(file_obj, "Time", "long_name", "Time", &
                                        str_len=len("Time"))
       call register_variable_attribute(file_obj, "Time", "units", "time level", &
                                        str_len=len("time level"))
       call register_variable_attribute(file_obj, "Time", "cartesian_axis", "T", str_len=1)
       call write_data(file_obj, "Time", 1)
    endif

  end subroutine fv_io_register_axis
  ! </SUBROUTINE> NAME="fv_io_register_axis"

!#####################################################################
  ! <SUBROUTINE NAME="fv_io_register_restart">
  !
  ! <DESCRIPTION>
  !   register restart field to be written out to restart file.
  ! </DESCRIPTION>
  subroutine  fv_io_register_restart(Atm)

    type(fv_atmos_type), intent(inout) :: Atm
    character(len=64) :: tracer_name
    character(len=8), dimension(1)  :: dim_names
    character(len=8), dimension(2)  :: dim_names_2d
    character(len=8), dimension(4)  :: dim_names_4d, dim_names_4d2, dim_names_4d3
    character(len=8), dimension(3)  :: dim_names_3d, dim_names_3d2
    integer           :: i, j
    integer           :: nt, ntracers, ntprog, ntdiag
    integer, dimension(:), allocatable :: buffer
    integer, parameter :: numx=1, numx_2d=2, numy=1, numy_2d=2, numz=1
    integer, dimension(1) :: xpos
    integer, dimension(2) :: xpos_2d
    integer, dimension(1) :: ypos
    integer, dimension(2) :: ypos_2d
    integer, dimension(numz) :: zsize

    dim_names_2d(1) = "xaxis_1"
    dim_names_2d(2) = "Time"
    dim_names_3d(1) = "xaxis_1"
    dim_names_3d(2) = "yaxis_2"
    dim_names_3d(3) = "Time"
    dim_names_3d2 = dim_names_3d
    dim_names_3d2(2) = "yaxis_1"
    dim_names_4d(1) = "xaxis_1"
    dim_names_4d(2) = "yaxis_1"
    dim_names_4d(3) = "zaxis_1"
    dim_names_4d(4) = "Time"
    dim_names_4d2 = dim_names_4d
    dim_names_4d2(1) = "xaxis_2"
    dim_names_4d2(2) = "yaxis_2"
    dim_names_4d3 = dim_names_4d
    dim_names_4d3(2) = "yaxis_2"

    ntprog = size(Atm%q,4)
    ntdiag = size(Atm%qdiag,4)
    ntracers = ntprog+ntdiag

    xpos = (/CENTER/)
    xpos_2d = (/CENTER, EAST/)
    ypos = (/CENTER/)
    ypos_2d = (/NORTH, CENTER/)

    ! fname = 'fv_core.res.nc'
    if (Atm%Fv_restart_is_open) then
       call register_axis(Atm%Fv_restart, "xaxis_1", size(Atm%ak(:), 1))
       call register_axis(Atm%Fv_restart, "Time", unlimited)
       if (.not. Atm%Fv_restart%is_readonly) then !if writing file
          call register_field(Atm%Fv_restart, dim_names_2d(1), axis_type, (/dim_names_2d(1)/))
          call register_variable_attribute(Atm%Fv_restart, dim_names_2d(1), "long_name", dim_names_2d(1), str_len=len(dim_names_2d(1)))
          call register_variable_attribute(Atm%Fv_restart, dim_names_2d(1), "units", "none", str_len=len("none"))
          call register_variable_attribute(Atm%Fv_restart, dim_names_2d(1), "cartesian_axis", "X", str_len=1)
          if (allocated(buffer)) deallocate(buffer)
          allocate(buffer(size(Atm%ak(:), 1)))
          do j = 1, size(Atm%ak(:), 1)
             buffer(j) = j
          end do
          call write_data(Atm%Fv_restart, dim_names_2d(1), buffer)
          deallocate(buffer)
          call register_field(Atm%Fv_restart, dim_names_2d(2), axis_type, (/dim_names_2d(2)/))
          call register_variable_attribute(Atm%Fv_restart, dim_names_2d(2), "long_name", dim_names_2d(2), str_len=len(dim_names_2d(2)))
          call register_variable_attribute(Atm%Fv_restart, dim_names_2d(2), "units", "time level", str_len=len("time level"))
          call register_variable_attribute(Atm%Fv_restart, dim_names_2d(2), "cartesian_axis", "T", str_len=1)
          call write_data(Atm%Fv_restart, dim_names_2d(2), 1)
       endif
       call register_restart_field (Atm%Fv_restart, 'ak', Atm%ak(:), dim_names_2d)
       call register_restart_field (Atm%Fv_restart, 'bk', Atm%bk(:), dim_names_2d)
       if (.not. Atm%Fv_restart%is_readonly) then !if writing file
         call register_variable_attribute(Atm%Fv_restart, 'ak', "long_name", "ak", str_len=len("ak"))
         call register_variable_attribute(Atm%Fv_restart, 'ak', "units", "none", str_len=len("none"))
         call register_variable_attribute(Atm%Fv_restart, 'bk', "long_name", "bk", str_len=len("bk"))
         call register_variable_attribute(Atm%Fv_restart, 'bk', "units", "none", str_len=len("none"))
       endif

    ! fname= 'fv_core.res'//trim(stile_name)//'.nc'
    elseif (Atm%Fv_restart_tile_is_open) then
       zsize = (/size(Atm%u,3)/)
       call fv_io_register_axis(Atm%Fv_restart_tile, numx=numx_2d, numy=numy_2d, xpos=xpos_2d, ypos=ypos_2d, numz=numz, zsize=zsize)

       !--- optionally include D-grid winds even if restarting from A-grid winds
       if (Atm%flagstruct%is_ideal_case) then
          call register_restart_field(Atm%Fv_restart_tile, 'u0', Atm%u0, &
               dim_names_4d, is_optional=.true.)
          call register_restart_field(Atm%Fv_restart_tile, 'v0', Atm%v0, &
               dim_names_4d2, is_optional=.true.)
       endif
       if (Atm%flagstruct%write_optional_dgrid_vel_rst .and. Atm%flagstruct%restart_from_agrid_winds) then
          call register_restart_field(Atm%Fv_restart_tile, 'u', Atm%u, &
               dim_names_4d, is_optional=.true.)
          call register_restart_field(Atm%Fv_restart_tile, 'v', Atm%v, &
               dim_names_4d2, is_optional=.true.)
       endif
       
       !--- include agrid winds in restarts for use in data assimilation or for restarting
       if (Atm%flagstruct%agrid_vel_rst .or. Atm%flagstruct%restart_from_agrid_winds) then
          call register_restart_field(Atm%Fv_restart_tile, 'ua', Atm%ua, &
               dim_names_4d3)
          call register_restart_field(Atm%Fv_restart_tile, 'va', Atm%va, &
               dim_names_4d3)
       endif
       
       if (.not. Atm%flagstruct%restart_from_agrid_winds) then
          call register_restart_field(Atm%Fv_restart_tile, 'u', Atm%u, &
               dim_names_4d)
          call register_restart_field(Atm%Fv_restart_tile, 'v', Atm%v, &
               dim_names_4d2)
       endif

       if (.not.Atm%flagstruct%hydrostatic) then
          if (Atm%flagstruct%make_nh) then ! Hydrostatic restarts dont have these variables
               call register_restart_field(Atm%Fv_restart_tile,  'W', Atm%w, dim_names_4d3, is_optional=.true.)
               call register_restart_field(Atm%Fv_restart_tile,  'DZ', Atm%delz, dim_names_4d3, is_optional=.true.)
               if ( Atm%flagstruct%hybrid_z ) then
                   call register_restart_field(Atm%Fv_restart_tile,  'ZE0', Atm%ze0, dim_names_4d3, is_optional=.true.)
               endif
          else !The restart file has the non-hydrostatic variables
               call register_restart_field(Atm%Fv_restart_tile,  'W', Atm%w, dim_names_4d3)
               call register_restart_field(Atm%Fv_restart_tile,  'DZ', Atm%delz, dim_names_4d3)
               if ( Atm%flagstruct%hybrid_z ) then
                   call register_restart_field(Atm%Fv_restart_tile,  'ZE0', Atm%ze0, dim_names_4d3)
               endif
          endif
       endif
       call register_restart_field(Atm%Fv_restart_tile,  'T', Atm%pt, dim_names_4d3)
       call register_restart_field(Atm%Fv_restart_tile,  'delp', Atm%delp, dim_names_4d3)
       call register_restart_field(Atm%Fv_restart_tile,  'phis', Atm%phis, dim_names_3d)

       if (.not. Atm%Fv_restart_tile%is_readonly) then !if writing file
         if (Atm%flagstruct%is_ideal_case) then
           if (variable_exists(Atm%Fv_restart_tile, 'u0')) then
             call register_variable_attribute(Atm%Fv_restart_tile, 'u0', "long_name", "u0", str_len=len("u0"))
             call register_variable_attribute(Atm%Fv_restart_tile, 'u0', "units", "none", str_len=len("none"))
           endif
           if (variable_exists(Atm%Fv_restart_tile, 'v0')) then
             call register_variable_attribute(Atm%Fv_restart_tile, 'v0', "long_name", "v0", str_len=len("v0"))
             call register_variable_attribute(Atm%Fv_restart_tile, 'v0', "units", "none", str_len=len("none"))
           endif
         endif
         if (variable_exists(Atm%Fv_restart_tile, 'u')) then
           call register_variable_attribute(Atm%Fv_restart_tile, 'u', "long_name", "u", str_len=len("u"))
           call register_variable_attribute(Atm%Fv_restart_tile, 'u', "units", "none", str_len=len("none"))
         endif
         if (variable_exists(Atm%Fv_restart_tile, 'v')) then
           call register_variable_attribute(Atm%Fv_restart_tile, 'v', "long_name", "v", str_len=len("v"))
           call register_variable_attribute(Atm%Fv_restart_tile, 'v', "units", "none", str_len=len("none"))
         endif
         if (variable_exists(Atm%Fv_restart_tile, 'W')) then
           call register_variable_attribute(Atm%Fv_restart_tile, 'W', "long_name", "W", str_len=len("W"))
           call register_variable_attribute(Atm%Fv_restart_tile, 'W', "units", "none", str_len=len("none"))
         endif
         if (variable_exists(Atm%Fv_restart_tile, 'DZ')) then
           call register_variable_attribute(Atm%Fv_restart_tile, 'DZ', "long_name", "DZ", str_len=len("DZ"))
           call register_variable_attribute(Atm%Fv_restart_tile, 'DZ', "units", "none", str_len=len("none"))
         endif
         if ( Atm%flagstruct%hybrid_z .and. variable_exists(Atm%Fv_restart_tile, 'ZEO')) then
           call register_variable_attribute(Atm%Fv_restart_tile, 'ZE0', "long_name", "ZE0", str_len=len("ZEO"))
           call register_variable_attribute(Atm%Fv_restart_tile, 'ZEO', "units", "none", str_len=len("none"))
         endif
         call register_variable_attribute(Atm%Fv_restart_tile, 'T', "long_name", "T", str_len=len("T"))
         call register_variable_attribute(Atm%Fv_restart_tile, 'T', "units", "none", str_len=len("none"))
         call register_variable_attribute(Atm%Fv_restart_tile, 'delp', "long_name", "delp", str_len=len("delp"))
         call register_variable_attribute(Atm%Fv_restart_tile, 'delp', "units", "none", str_len=len("none"))
         call register_variable_attribute(Atm%Fv_restart_tile, 'phis', "long_name", "phis", str_len=len("phis"))
         call register_variable_attribute(Atm%Fv_restart_tile, 'phis', "units", "none", str_len=len("none"))
         if (variable_exists(Atm%Fv_restart_tile, 'ua')) then
           call register_variable_attribute(Atm%Fv_restart_tile, 'ua', "long_name", "ua", str_len=len("ua"))
           call register_variable_attribute(Atm%Fv_restart_tile, 'ua', "units", "none", str_len=len("none"))
         endif
         if (variable_exists(Atm%Fv_restart_tile, 'va')) then
           call register_variable_attribute(Atm%Fv_restart_tile, 'va', "long_name", "va", str_len=len("va"))
           call register_variable_attribute(Atm%Fv_restart_tile, 'va', "units", "none", str_len=len("none"))
         endif
       endif

    ! fname = 'fv_srf_wnd.res'//trim(stile_name)//'.nc
    elseif (Atm%Rsf_restart_is_open) then
       call fv_io_register_axis(Atm%Rsf_restart, numx=numx, numy=numy, xpos=xpos, ypos=ypos)
       call register_restart_field(Atm%Rsf_restart, 'u_srf', Atm%u_srf, dim_names_3d2)
       call register_restart_field(Atm%Rsf_restart, 'v_srf', Atm%v_srf, dim_names_3d2)
#ifdef SIM_PHYS
       call register_restart_field(Atm%Rsf_restart, 'ts', Atm%ts, dim_names_3d2)
#endif
       if (.not. Atm%Rsf_restart%is_readonly) then !if writing file
         call register_variable_attribute(Atm%Rsf_restart, 'u_srf', "long_name", "u_srf", str_len=len("u_srf"))
         call register_variable_attribute(Atm%Rsf_restart, 'u_srf', "units", "none", str_len=len("none"))
         call register_variable_attribute(Atm%Rsf_restart, 'v_srf', "long_name", "v_srf", str_len=len("v_srf"))
         call register_variable_attribute(Atm%Rsf_restart, 'v_srf', "units", "none", str_len=len("none"))
#ifdef SIM_PHYS
         call register_variable_attribute(Atm%Rsf_restart, 'ts', "long_name", "ts", str_len=len("ts"))
         call register_variable_attribute(Atm%Rsf_restart, 'ts', "units", "none", str_len=len("none"))
#endif
       endif


    ! fname = 'mg_drag.res'//trim(stile_name)//'.nc'
    elseif (Atm%Mg_restart_is_open) then
       call fv_io_register_axis(Atm%Mg_restart, numx=numx, numy=numy, xpos=xpos, ypos=ypos)
       call register_restart_field (Atm%Mg_restart, 'ghprime', Atm%sgh, dim_names_3d2)
       if (.not. Atm%Mg_restart%is_readonly) then !if writing file
         call register_variable_attribute(Atm%Mg_restart, 'ghprime', "long_name", "ghprime", str_len=len("ghprime"))
         call register_variable_attribute(Atm%Mg_restart, 'ghprime', "units", "none", str_len=len("none"))
       endif

    ! fname = 'fv_land.res'//trim(stile_name)//'.nc'
    elseif (Atm%Lnd_restart_is_open) then
       call fv_io_register_axis(Atm%Lnd_restart, numx=numx, numy=numy, xpos=xpos, ypos=ypos)
       call register_restart_field (Atm%Lnd_restart, 'oro', Atm%oro, dim_names_3d2)
       if (.not. Atm%Lnd_restart%is_readonly) then !if writing file
         call register_variable_attribute(Atm%Lnd_restart, 'oro', "long_name", "oro", str_len=len("oro"))
         call register_variable_attribute(Atm%Lnd_restart, 'oro', "units", "none", str_len=len("none"))
       endif

    ! fname = 'fv_tracer.res'//trim(stile_name)//'.nc'
    elseif (Atm%Tra_restart_is_open) then
       zsize = (/size(Atm%q,3)/)
       call fv_io_register_axis(Atm%Tra_restart, numx=numx, numy=numy, xpos=xpos, ypos=ypos, numz=numz, zsize=zsize)
       do nt = 1, ntprog
          call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
          call register_restart_field(Atm%Tra_restart, tracer_name, Atm%q(:,:,:,nt), &
                       dim_names_4d, is_optional=.true.)
          if (variable_exists(Atm%Tra_restart, tracer_name) .and. .not. Atm%Tra_restart%is_readonly) then
             call register_variable_attribute(Atm%Tra_restart, tracer_name, "long_name", tracer_name, str_len=len(tracer_name))
             call register_variable_attribute(Atm%Tra_restart, tracer_name, "units", "none", str_len=len("none"))
          endif
       enddo
       do nt = ntprog+1, ntracers
          call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
          call register_restart_field(Atm%Tra_restart, tracer_name, Atm%qdiag(:,:,:,nt), &
                       dim_names_4d, is_optional=.true.)
          if (variable_exists(Atm%Tra_restart, tracer_name) .and. .not. Atm%Tra_restart%is_readonly) then
             call register_variable_attribute(Atm%Tra_restart, tracer_name, "long_name", tracer_name, str_len=len(tracer_name))
             call register_variable_attribute(Atm%Tra_restart, tracer_name, "units", "none", str_len=len("none"))
          endif
       enddo
    endif
  end subroutine  fv_io_register_restart
  ! </SUBROUTINE> NAME="fv_io_register_restart"


  !#####################################################################
  ! <SUBROUTINE NAME="fv_io_read_restart">
  !
  ! <DESCRIPTION>
  ! Write the fv core restart quantities
  ! </DESCRIPTION>
  subroutine  fv_io_read_restart(fv_domain,Atm,prefix,directory)
    type(domain2d),      intent(inout) :: fv_domain
    type(fv_atmos_type), intent(inout) :: Atm(:)
    character(len=*), optional, intent(in) :: prefix
    character(len=*), optional, intent(in) :: directory

    character(len=64)    :: tracer_name
    integer              :: isc, iec, jsc, jec, n, nt, nk, ntracers
    integer              :: ntileMe
    integer              :: ks, ntiles
    real                 :: ptop

    character (len=:), allocatable :: dir, pre, suffix, fname
    character(len=128) :: tracer_longname, tracer_units
    character(len=1) :: tile_num
    integer, allocatable, dimension(:) :: pes !< Array of the pes in the current pelist

    pre = ''
    if (present(prefix)) pre = ''//trim(prefix)//'.'
    dir = 'INPUT'
    if (present(directory)) dir = trim(directory)

    allocate(pes(mpp_npes()))
    call mpp_get_current_pelist(pes)

    suffix = ''
    fname = ''//trim(dir)//'/'//trim(pre)//'fv_core.res.nc'
    Atm(1)%Fv_restart_is_open = open_file(Atm(1)%Fv_restart,fname,"read", is_restart=.true., pelist=pes)
    if (Atm(1)%Fv_restart_is_open) then
      call fv_io_register_restart(Atm(1))
      call read_restart(Atm(1)%Fv_restart)
      call close_file(Atm(1)%Fv_restart)
      Atm(1)%Fv_restart_is_open = .false.
    endif
    deallocate(pes)

    if (Atm(1)%flagstruct%external_eta) then
       call set_external_eta(Atm(1)%ak, Atm(1)%bk, Atm(1)%ptop, Atm(1)%ks)
    endif

    if ( use_ncep_sst .or. Atm(1)%flagstruct%nudge .or. Atm(1)%flagstruct%ncep_ic ) then
       call mpp_error(NOTE, 'READING FROM SST_RESTART DISABLED')
       !call restore_state(Atm(1)%SST_restart)
    endif

    ntiles = mpp_get_ntile_count(fv_domain)
    !If the number of tiles is equal to 1, and it is not a nested case add the ".tile1" suffix to the filename
    if (ntiles == 1 .and. .not. Atm(1)%neststruct%nested) then
       suffix = ''//trim(suffix)//'.tile1'
    endif

    fname = ''//trim(dir)//'/'//trim(pre)//'fv_core.res'//trim(suffix)//'.nc'
    Atm(1)%Fv_restart_tile_is_open = open_file(Atm(1)%Fv_restart_tile, fname, "read", fv_domain, is_restart=.true.)
    if (Atm(1)%Fv_restart_tile_is_open) then
      call fv_io_register_restart(Atm(1))
      call read_restart(Atm(1)%Fv_restart_tile, ignore_checksum=Atm(1)%flagstruct%ignore_rst_cksum)
      call close_file(Atm(1)%Fv_restart_tile)
      Atm(1)%Fv_restart_tile_is_open = .false.
      if (Atm(1)%flagstruct%restart_from_agrid_winds) then
         call cubed_a2d(Atm(1)%npx, Atm(1)%npy, Atm(1)%npz, &
              Atm(1)%ua, Atm(1)%va, Atm(1)%u, Atm(1)%v, &
              Atm(1)%gridstruct, Atm(1)%domain, Atm(1)%bd)
         call mpp_update_domains(Atm(1)%u, Atm(1)%v, Atm(1)%domain, gridtype=DGRID_NE, complete=.true.)
      endif
    endif

!--- restore data for fv_tracer - if it exists
    fname = ''//trim(dir)//'/'//trim(pre)//'fv_tracer.res'//trim(suffix)//'.nc'
    Atm(1)%Tra_restart_is_open = open_file(Atm(1)%Tra_restart, fname, "read", fv_domain, is_restart=.true.)
    if (Atm(1)%Tra_restart_is_open) then
      call fv_io_register_restart(Atm(1))
      call read_restart(Atm(1)%Tra_restart, ignore_checksum=Atm(1)%flagstruct%ignore_rst_cksum)
      call close_file(Atm(1)%Tra_restart)
      Atm(1)%Tra_restart_is_open = .false.
    else
      call mpp_error(NOTE,'==> Warning from fv_read_restart: Expected file '//trim(fname)//' does not exist')
    endif

!--- restore data for surface winds - if it exists
    fname = ''//trim(dir)//'/'//trim(pre)//'fv_srf_wnd.res'//trim(suffix)//'.nc'
    Atm(1)%Rsf_restart_is_open = open_file(Atm(1)%Rsf_restart, fname, "read", fv_domain, is_restart=.true.)
    if (Atm(1)%Rsf_restart_is_open) then
      Atm(1)%flagstruct%srf_init = .true.
      call fv_io_register_restart(Atm(1))
      call read_restart(Atm(1)%Rsf_restart, ignore_checksum=Atm(1)%flagstruct%ignore_rst_cksum)
      call close_file(Atm(1)%Rsf_restart)
      Atm(1)%Rsf_restart_is_open = .false.
    else
      call mpp_error(NOTE,'==> Warning from fv_read_restart: Expected file '//trim(fname)//' does not exist')
      Atm(1)%flagstruct%srf_init = .false.
    endif

    if ( Atm(1)%flagstruct%fv_land ) then
!--- restore data for mg_drag - if it exists
         fname = ''//trim(dir)//'/'//trim(pre)//'mg_drag.res'//trim(suffix)//'.nc'
         Atm(1)%Mg_restart_is_open = open_file(Atm(1)%Mg_restart, fname, "read", fv_domain, is_restart=.true.)
         if (Atm(1)%Mg_restart_is_open) then
           call fv_io_register_restart(Atm(1))
           call read_restart(Atm(1)%Mg_restart, ignore_checksum=Atm(1)%flagstruct%ignore_rst_cksum)
           call close_file(Atm(1)%Mg_restart)
           Atm(1)%Mg_restart_is_open = .false.
         else
           call mpp_error(NOTE,'==> Warning from fv_read_restart: Expected file '//trim(fname)//' does not exist')
         endif
!--- restore data for fv_land - if it exists
         fname = ''//trim(dir)//'/'//trim(pre)//'/fv_land.res'//trim(suffix)//'.nc'
         Atm(1)%Lnd_restart_is_open = open_file(Atm(1)%Lnd_restart, fname, "read", fv_domain, is_restart=.true.)
         if (Atm(1)%Lnd_restart_is_open) then
           call fv_io_register_restart(Atm(1))
           call read_restart(Atm(1)%Lnd_restart, ignore_checksum=Atm(1)%flagstruct%ignore_rst_cksum)
           call close_file(Atm(1)%Lnd_restart)
           Atm(1)%Lnd_restart_is_open = .false.
         else
           call mpp_error(NOTE,'==> Warning from fv_read_restart: Expected file '//trim(fname)//' does not exist')
         endif
    endif

    return

  end subroutine  fv_io_read_restart
  ! </SUBROUTINE> NAME="fv_io_read_restart"
  !#####################################################################


  subroutine fv_io_read_tracers(Atm)
    type(fv_atmos_type), intent(inout) :: Atm(:)
    integer :: ntracers, ntprog, nt, isc, iec, jsc, jec
    character(len=6) :: stile_name
    character(len=64):: fname, tracer_name
    type(FmsNetcdfDomainFile_t) :: Tra_restart_r
    integer :: ntiles

    isc = Atm(1)%bd%isc
    iec = Atm(1)%bd%iec
    jsc = Atm(1)%bd%jsc
    jec = Atm(1)%bd%jec
    call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers, num_prog=ntprog)

! fix for single tile runs where you need fv_core.res.nc and fv_core.res.tile1.nc
    ntiles = mpp_get_ntile_count(Atm(1)%domain_for_read)
    if(ntiles == 1 .and. .not. Atm(1)%neststruct%nested) then
       stile_name = '.tile1'
    else
       stile_name = ''
    endif

    fname = 'INPUT/fv_tracer.res'//trim(stile_name)//'.nc'

    if (open_file(Tra_restart_r,fname,"read",Atm(1)%domain_for_read, is_restart=.true.)) then
    do nt = 2, ntprog
       call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
       call set_tracer_profile (MODEL_ATMOS, nt, Atm(1)%q(isc:iec,jsc:jec,:,nt)  )
       if (variable_exists(Tra_restart_r, tracer_name)) then
          call read_data(Tra_restart_r, tracer_name, Atm(1)%q(:,:,:,nt))
       endif
    enddo
    do nt = ntprog+1, ntracers
       call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
       call set_tracer_profile (MODEL_ATMOS, nt, Atm(1)%qdiag(isc:iec,jsc:jec,:,nt)  )
       if (variable_exists(Tra_restart_r, tracer_name)) then
          call read_data (Tra_restart_r, tracer_name, Atm(1)%qdiag(:,:,:,nt))
       endif
    enddo
    call close_file(Tra_restart_r)
    else
      call mpp_error(NOTE,'==> Warning from fv_io_read_tracers: Expected file '//trim(fname)//' does not exist')
    endif
    return

  end subroutine  fv_io_read_tracers


  subroutine  remap_restart(Atm)
  use fv_mapz_mod,       only: rst_remap

    type(fv_atmos_type), intent(inout) :: Atm(:)

    character(len=64)    :: fname, tracer_name
    character(len=6)     :: stile_name
    integer              :: isc, iec, jsc, jec, nt, nk, ntracers, ntprog, ntdiag
    integer              :: isd, ied, jsd, jed
    integer              :: ntiles

    type(domain2d) :: fv_domain
    type(FmsNetcdfDomainFile_t) :: FV_tile_restart_r, Tra_restart_r
    type(FmsNetcdfFile_t)       :: Fv_restart_r
    integer, allocatable, dimension(:) :: pes !< Array of the pes in the current pelist

!
!-------------------------------------------------------------------------
    real, allocatable:: ak_r(:), bk_r(:), u0_r(:,:,:), v0_r(:,:,:)
    real, allocatable:: u_r(:,:,:), v_r(:,:,:), pt_r(:,:,:), delp_r(:,:,:)
    real, allocatable:: w_r(:,:,:), delz_r(:,:,:), ze0_r(:,:,:)
    real, allocatable:: q_r(:,:,:,:), qdiag_r(:,:,:,:)
!-------------------------------------------------------------------------
    integer npz, npz_rst, ng
    integer i,j,k

    fv_domain = Atm(1)%domain_for_read
    npz     = Atm(1)%npz       ! run time z dimension
    npz_rst = Atm(1)%flagstruct%npz_rst   ! restart z dimension
    isc = Atm(1)%bd%isc; iec = Atm(1)%bd%iec; jsc = Atm(1)%bd%jsc; jec = Atm(1)%bd%jec
    ng = Atm(1)%ng

    isd = isc - ng;  ied = iec + ng
    jsd = jsc - ng;  jed = jec + ng


!   call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)
    ntprog = size(Atm(1)%q,4)  ! Temporary until we get tracer manager integrated
    ntdiag = size(Atm(1)%qdiag,4)
    ntracers = ntprog+ntdiag

!    ntileMe = size(Atm(:))  ! This will have to be modified for mult tiles per PE


! Allocate arrays for reading old restart file:
    allocate ( ak_r(npz_rst+1) )
    allocate ( bk_r(npz_rst+1) )

    allocate ( u0_r(isc:iec,  jsc:jec+1,npz_rst) )
    allocate ( v0_r(isc:iec+1,jsc:jec  ,npz_rst) )
    allocate ( u_r(isc:iec,  jsc:jec+1,npz_rst) )
    allocate ( v_r(isc:iec+1,jsc:jec  ,npz_rst) )

    allocate (   pt_r(isc:iec, jsc:jec,  npz_rst) )
    allocate ( delp_r(isc:iec, jsc:jec,  npz_rst) )
    allocate (    q_r(isc:iec, jsc:jec,  npz_rst, ntprog) )
    allocate (qdiag_r(isc:iec, jsc:jec,  npz_rst, ntprog+1:ntracers) )

    if ( (.not.Atm(1)%flagstruct%hydrostatic) .and. (.not.Atm(1)%flagstruct%make_nh) ) then
           allocate (    w_r(isc:iec, jsc:jec,  npz_rst) )
           allocate ( delz_r(isc:iec, jsc:jec,  npz_rst) )
           if ( Atm(1)%flagstruct%hybrid_z )   &
           allocate ( ze0_r(isc:iec, jsc:jec,  npz_rst+1) )
    endif

    fname = 'INPUT/fv_core.res.nc'
    allocate(pes(mpp_npes()))
    call mpp_get_current_pelist(pes)
    if (open_file(Fv_restart_r,fname,"read", is_restart=.true., pelist=pes)) then
       call read_data(Fv_restart_r, 'ak', ak_r(:))
       call read_data(Fv_restart_r, 'bk', bk_r(:))
       call close_file(Fv_restart_r)
    endif
    deallocate(pes)

! fix for single tile runs where you need fv_core.res.nc and fv_core.res.tile1.nc
    ntiles = mpp_get_ntile_count(fv_domain)
    if(ntiles == 1 .and. .not. Atm(1)%neststruct%nested) then
       stile_name = '.tile1'
    else
       stile_name = ''
    endif

!!!! A NOTE about file names
!!! file_exist() needs the full relative path, including INPUT/
!!! But register_restart_field ONLY looks in INPUT/ and so JUST needs the file name!!

       fname = 'INPUT/fv_core.res'//trim(stile_name)//'.nc'
       if (open_file(Fv_tile_restart_r, fname, "read", fv_domain, is_restart=.true.)) then
          if (Atm(1)%flagstruct%is_ideal_case) then
             call read_data(Fv_tile_restart_r, 'u0', u0_r)
             call read_data(Fv_tile_restart_r, 'v0', v0_r)
          endif
          call read_data(Fv_tile_restart_r, 'u', u_r)
          call read_data(Fv_tile_restart_r, 'v', v_r)
          if (.not.Atm(1)%flagstruct%hydrostatic) then
             call read_data(Fv_tile_restart_r, 'W', w_r)
             call read_data(Fv_tile_restart_r, 'DZ', delz_r)
             if ( Atm(1)%flagstruct%hybrid_z ) then
                call read_data(Fv_tile_restart_r, 'ZE0', ze0_r)
             endif
          endif
          call read_data(Fv_tile_restart_r, 'T', pt_r)
          call read_data(Fv_tile_restart_r, 'delp', delp_r)
          call read_data(Fv_tile_restart_r, 'phis', Atm(1)%phis)
          call close_file(FV_tile_restart_r)
       endif

       fname = 'INPUT/fv_srf_wnd.res'//trim(stile_name)//'.nc'
       Atm(1)%Rsf_restart_is_open = open_file(Atm(1)%Rsf_restart, fname, "read", fv_domain, is_restart=.true.)
       if (Atm(1)%Rsf_restart_is_open) then
          Atm%flagstruct%srf_init = .true.
          call fv_io_register_restart(Atm(1))
          call read_restart(Atm(1)%Rsf_restart, ignore_checksum=Atm(1)%flagstruct%ignore_rst_cksum)
          call close_file(Atm(1)%Rsf_restart)
          Atm(1)%Rsf_restart_is_open = .false.
       else
          call mpp_error(NOTE,'==> Warning from remap_restart: Expected file '//trim(fname)//' does not exist')
          Atm%flagstruct%srf_init = .false.
       endif

       if ( Atm(1)%flagstruct%fv_land ) then
!--- restore data for mg_drag - if it exists
         fname = 'INPUT/mg_drag.res'//trim(stile_name)//'.nc'
         Atm(1)%Mg_restart_is_open = open_file(Atm(1)%Mg_restart, fname, "read", fv_domain, is_restart=.true.)
         if (Atm(1)%Mg_restart_is_open) then
            call read_data(Atm(1)%Mg_restart, 'ghprime', Atm(1)%sgh)
            call close_file(Atm(1)%Mg_restart)
            Atm(1)%Mg_restart_is_open = .false.
         else
           call mpp_error(NOTE,'==> Warning from remap_restart: Expected file '//trim(fname)//' does not exist')
         endif
!--- restore data for fv_land - if it exists
         fname = 'INPUT/fv_land.res'//trim(stile_name)//'.nc'
         Atm(1)%Lnd_restart_is_open = open_file(Atm(1)%Lnd_restart, fname, "read", fv_domain, is_restart=.true.)
         if (Atm(1)%Lnd_restart_is_open) then
           call read_data(Atm(1)%Lnd_restart, 'oro', Atm(1)%oro)
           call close_file(Atm(1)%Lnd_restart)
           Atm(1)%Lnd_restart_is_open = .false.
         else
           call mpp_error(NOTE,'==> Warning from remap_restart: Expected file '//trim(fname)//' does not exist')
         endif
       endif

       fname = 'INPUT/fv_tracer.res'//trim(stile_name)//'.nc'
       if (open_file(Tra_restart_r, fname, "read", fv_domain, is_restart=.true.)) then
         do nt = 1, ntprog
            call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
            call set_tracer_profile (MODEL_ATMOS, nt, q_r(isc:iec,jsc:jec,:,nt)  )
            if (variable_exists(Tra_restart_r, tracer_name)) then
               call read_data(Tra_restart_r, tracer_name, q_r(:,:,:,nt))
            endif
         enddo
         do nt = ntprog+1, ntracers
            call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
            call set_tracer_profile (MODEL_ATMOS, nt, qdiag_r(isc:iec,jsc:jec,:,nt)  )
            if (variable_exists(Tra_restart_r, tracer_name)) then
               call read_data (Tra_restart_r, tracer_name, qdiag_r(:,:,:,nt))
            endif
         enddo
         call close_file(Tra_restart_r)
       else
         call mpp_error(NOTE,'==> Warning from remap_restart: Expected file '//trim(fname)//' does not exist')
       endif

!      ====== PJP added DA functionality ======
       if (Atm(1)%flagstruct%read_increment) then
          ! print point in middle of domain for a sanity check
          i = (isc + iec)/2
          j = (jsc + jec)/2
          k = npz_rst/2
          if( is_master() ) write(*,*) 'Calling read_da_inc',pt_r(i,j,k)
          call read_da_inc(Atm(1), Atm(1)%domain)
          if( is_master() ) write(*,*) 'Back from read_da_inc',pt_r(i,j,k)
       endif
!      ====== end PJP added DA functionailty======

       call rst_remap(npz_rst, npz, isc, iec, jsc, jec, isd, ied, jsd, jed, ntracers, ntprog,      &
                      delp_r, u0_r, v0_r, u_r, v_r, w_r, delz_r, pt_r, q_r, qdiag_r, Atm(1)%delp,  &
                      Atm(1)%u0, Atm(1)%v0, Atm(1)%u, Atm(1)%v, Atm(1)%w, Atm(1)%delz, Atm(1)%pt,  &
                      Atm(1)%q, Atm(1)%qdiag, ak_r,  bk_r, Atm(1)%ptop, Atm(1)%ak, Atm(1)%bk,      &
                      Atm(1)%flagstruct%hydrostatic, Atm(1)%flagstruct%make_nh, Atm(1)%domain,     &
                      Atm(1)%gridstruct%square_domain, Atm(1)%flagstruct%is_ideal_case)
    !end do

    deallocate( ak_r )
    deallocate( bk_r )
    deallocate( u0_r )
    deallocate( v0_r )
    deallocate( u_r )
    deallocate( v_r )
    deallocate( pt_r )
    deallocate( delp_r )
    deallocate( q_r )
    deallocate( qdiag_r )

    if ( (.not.Atm(1)%flagstruct%hydrostatic) .and. (.not.Atm(1)%flagstruct%make_nh) ) then
         deallocate ( w_r )
         deallocate ( delz_r )
         if ( Atm(1)%flagstruct%hybrid_z ) deallocate ( ze0_r )
    endif

  end subroutine  remap_restart


  !#####################################################################
  ! <SUBROUTINE NAME="fv_io_register_nudge_restart">
  !
  ! <DESCRIPTION>
  !   register restart nudge field to be written out to restart file.
  ! </DESCRIPTION>
  subroutine  fv_io_register_nudge_restart(Atm)
    type(fv_atmos_type), intent(inout) :: Atm(:)
    character(len=64) :: fname
    integer           :: id_restart

! use_ncep_sst may not be initialized at this point?
    call mpp_error(NOTE, 'READING FROM SST_restart DISABLED')

  end subroutine  fv_io_register_nudge_restart

  ! </SUBROUTINE> NAME="fv_io_register_nudge_restart"

  !#####################################################################
  ! <SUBROUTINE NAME="fv_io_write_restart">
  !
  ! <DESCRIPTION>
  ! Write the fv core restart quantities
  ! </DESCRIPTION>
  subroutine  fv_io_write_restart(Atm, prefix, directory)

    type(fv_atmos_type),        intent(inout) :: Atm
    character(len=*), optional, intent(in) :: prefix
    character(len=*), optional, intent(in) :: directory

    character (len=:), allocatable :: dir, pre, fname, suffix
    integer :: ntiles
    logical :: tile_file_exists
    type(domain2d)     :: fv_domain
    character(len=1) :: tile_num
    integer, allocatable, dimension(:) :: pes !< Array of the pes in the current pelist
    fv_domain = Atm%domain

    if ( (use_ncep_sst .or. Atm%flagstruct%nudge) .and. .not. Atm%gridstruct%nested ) then
       !call save_restart(Atm%SST_restart, timestamp)
    endif

    pre = ''
    dir = 'RESTART'
    if (present(prefix)) pre = ''//trim(prefix)//'.'
    if (present(directory)) dir = trim(directory)

    suffix = ''

    fname = ''//trim(dir)//'/'//trim(pre)//'fv_core.res.nc'
    allocate(pes(mpp_npes()))
    call mpp_get_current_pelist(pes)
    Atm%Fv_restart_is_open = open_file(Atm%Fv_restart, fname, "overwrite", is_restart=.true., pelist=pes)
    if (Atm%Fv_restart_is_open) then
       call fv_io_register_restart(Atm)
       call write_restart(Atm%Fv_restart)
       call close_file(Atm%Fv_restart)
       Atm%Fv_restart_is_open = .false.
    endif
    deallocate(pes)

    ntiles = mpp_get_ntile_count(fv_domain)
    !If the number of tiles is equal to 1, and it is not a nested case add the ".tile1" suffix to the filename
    if (ntiles == 1 .and. .not. Atm%neststruct%nested) then
       suffix = ''//trim(suffix)//'.tile1'
    endif

    fname = ''//trim(dir)//'/'//trim(pre)//'fv_core.res'//trim(suffix)//'.nc'
    Atm%Fv_restart_tile_is_open = open_file(Atm%Fv_restart_tile, fname, "overwrite", fv_domain, is_restart=.true.)
    if (Atm%Fv_restart_tile_is_open) then
       call fv_io_register_restart(Atm)
       call write_restart (Atm%Fv_restart_tile)
       call close_file (Atm%Fv_restart_tile)
       Atm%Fv_restart_tile_is_open = .false.
    endif

    fname = ''//trim(dir)//'/'//trim(pre)//'fv_srf_wnd.res'//trim(suffix)//'.nc'
    Atm%Rsf_restart_is_open = open_file(Atm%Rsf_restart, fname, "overwrite", fv_domain, is_restart=.true.)
    if (Atm%Rsf_restart_is_open) then
       call fv_io_register_restart(Atm)
       call write_restart (Atm%Rsf_restart)
       call close_file (Atm%Rsf_restart)
       Atm%Rsf_restart_is_open = .false.
    endif

    if ( Atm%flagstruct%fv_land ) then
       fname = ''//trim(dir)//'/'//trim(pre)//'mg_drag.res'//trim(suffix)//'.nc'
       Atm%Mg_restart_is_open = open_file(Atm%Mg_restart, fname, "overwrite", fv_domain, is_restart=.true.)
       if (Atm%Mg_restart_is_open) then
          call fv_io_register_restart(Atm)
          call write_restart(Atm%Mg_restart)
          call close_file(Atm%Mg_restart)
          Atm%Mg_restart_is_open = .false.
       endif

       fname = ''//trim(dir)//'/'//trim(pre)//'/fv_land.res'//trim(suffix)//'.nc'
       Atm%Lnd_restart_is_open = open_file(Atm%Lnd_restart, fname, "overwrite", fv_domain, is_restart=.true.)
       if (Atm%Lnd_restart_is_open) then
          call fv_io_register_restart(Atm)
          call write_restart(Atm%Lnd_restart)
          call close_file(Atm%Lnd_restart)
          Atm%Lnd_restart_is_open = .false.
       endif
    endif

    fname = ''//trim(dir)//'/'//trim(pre)//'fv_tracer.res'//trim(suffix)//'.nc'
    Atm%Tra_restart_is_open = open_file(Atm%Tra_restart, fname, "overwrite", fv_domain, is_restart=.true.)
    if (Atm%Tra_restart_is_open) then
       call fv_io_register_restart(Atm)
       call write_restart(Atm%Tra_restart)
       call close_file(Atm%Tra_restart)
       Atm%Tra_restart_is_open = .false.
    endif

  end subroutine fv_io_write_restart


  subroutine register_bcs_2d(Atm, BCfile_ne, BCfile_sw, fname_ne, fname_sw, &
                             var_name, var, var_bc, istag, jstag)
    type(fv_atmos_type),      intent(in)    :: Atm
    type(FmsNetcdfFile_t),  intent(inout) :: BCfile_ne, BCfile_sw
    character(len=120),       intent(in)    :: fname_ne, fname_sw
    character(len=*),         intent(in)    :: var_name
    real, dimension(:,:),     intent(in), optional :: var
    type(fv_nest_BC_type_3D), intent(in), optional :: var_bc
    integer,                  intent(in), optional :: istag, jstag

    integer :: npx, npy, i_stag, j_stag
    integer :: is, ie, js, je, isd, ied, jsd, jed, n
    integer :: x_halo, y_halo, x_halo_ns, id_restart
    integer :: layout(2), global_size(2), indices(4)
    integer, allocatable, dimension(:) :: x1_pelist, y1_pelist
    integer, allocatable, dimension(:) :: x2_pelist, y2_pelist
    logical :: is_root_pe

    i_stag = 0
    j_stag = 0
    if (present(istag)) i_stag = i_stag
    if (present(jstag)) j_stag = j_stag
    call mpp_get_global_domain(Atm%domain, xsize = npx, ysize = npy, position=CORNER )
    call mpp_get_data_domain(Atm%domain, isd, ied, jsd, jed )
    call mpp_get_compute_domain(Atm%domain, is, ie, js, je )
    call mpp_get_layout(Atm%domain, layout)
    allocate (x1_pelist(layout(1)))
    allocate (y1_pelist(layout(2)))
    allocate (x2_pelist(layout(1)))
    allocate (y2_pelist(layout(2)))
    x_halo = is-isd
    y_halo = js-jsd
! define west and east pelist
    do n = 1,layout(2)
      y1_pelist(n)=mpp_root_pe()+layout(1)*n-1
      y2_pelist(n)=mpp_root_pe()+layout(1)*(n-1)
    enddo
! define south and north pelist
    do n = 1,layout(1)
      x1_pelist(n)=mpp_root_pe()+layout(1)*(layout(2)-1)+(n-1)
      x2_pelist(n)=mpp_root_pe()+(n-1)
    enddo
! declare the pelists inside of mpp (creates the MPI communicator)
    call mpp_declare_pelist(x1_pelist)
    call mpp_declare_pelist(x2_pelist)
    call mpp_declare_pelist(y1_pelist)
    call mpp_declare_pelist(y2_pelist)

!EAST & WEST
!set defaults for west/east halo regions
    indices(1) = 1
    indices(2) = x_halo
    indices(3) = jsd
    indices(4) = jed+j_stag
    global_size(1) = x_halo
    global_size(2) = npy-1+2*y_halo+j_stag

!define west root_pe
    is_root_pe = .FALSE.
    if (is.eq.1 .and. js.eq.1) is_root_pe = .TRUE.
!register west halo data in t1
    if (present(var_bc) .and. Atm%neststruct%BCfile_sw_is_open) then
       call register_restart_field(BCfile_sw, trim(var_name)//'_west_t1', &
                                   var_bc%west_t1, indices, global_size, &
                                   y2_pelist, is_root_pe, jshift=y_halo)
       if (.not. BCfile_sw%is_readonly) then !if writing file
         call register_variable_attribute(BCfile_sw, trim(var_name)//'_west_t1', &
                                          "long_name", trim(var_name)//'_west_t1', &
                                          str_len=len(trim(var_name)//'_west_t1'))
         call register_variable_attribute(BCfile_sw, trim(var_name)//'_west_t1', &
                                          "units", "none", str_len=len("none"))
       endif
    endif
!register west prognostic halo data
    if (present(var) .and. Atm%neststruct%BCfile_sw_is_open) then
       call register_restart_field(BCfile_sw, trim(var_name)//'_west', &
                                   var, indices, global_size, &
                                   y2_pelist, is_root_pe, jshift=y_halo)
       if (.not. BCfile_sw%is_readonly) then !if writing file
         call register_variable_attribute(BCfile_sw, trim(var_name)//'_west', &
                                          "long_name", trim(var_name)//'_west', &
                                          str_len=len(trim(var_name)//'_west'))
         call register_variable_attribute(BCfile_sw, trim(var_name)//'_west', &
                                          "units", "none", str_len=len("none"))
       endif
    endif

!define east root_pe
    is_root_pe = .FALSE.
    if (ie.eq.npx-1 .and. je.eq.npy-1) is_root_pe = .TRUE.
!register east halo data in t1
    if (present(var_bc) .and. Atm%neststruct%BCfile_ne_is_open) then
       call register_restart_field(BCfile_ne, trim(var_name)//'_east_t1', &
                                   var_bc%east_t1, indices, global_size, &
                                   y1_pelist, is_root_pe, jshift=y_halo)
       if (.not. BCfile_ne%is_readonly) then !if writing file
         call register_variable_attribute(BCfile_ne, trim(var_name)//'_east_t1', &
                                          "long_name", trim(var_name)//'_east_t1', &
                                          str_len=len(trim(var_name)//'_east_t1'))
         call register_variable_attribute(BCfile_ne, trim(var_name)//'_east_t1', &
                                          "units", "none", str_len=len("none"))
       endif
    endif

!reset indices for prognostic variables in the east halo
    indices(1) = ied-x_halo+1+i_stag
    indices(2) = ied+i_stag
!register east prognostic halo data
    if (present(var) .and. Atm%neststruct%BCfile_ne_is_open) then
       call register_restart_field(BCfile_ne, trim(var_name)//'_east', &
                                   var, indices, global_size, &
                                   y1_pelist, is_root_pe, jshift=y_halo, &
                                   x_halo=(size(var,1)-x_halo), ishift=-(ie+i_stag))
       if (.not. BCfile_ne%is_readonly) then !if writing file
         call register_variable_attribute(BCfile_ne, trim(var_name)//'_east', &
                                          "long_name", trim(var_name)//'_east', &
                                          str_len=len(trim(var_name)//'_east'))
         call register_variable_attribute(BCfile_ne, trim(var_name)//'_east', &
                                          "units", "none", str_len=len("none"))
       endif
    endif

!NORTH & SOUTH
!set defaults for north/south halo regions
    indices(1) = isd
    indices(2) = ied+i_stag
    indices(3) = 1
    indices(4) = y_halo
    global_size(1) = npx-1+i_stag
    global_size(2) = y_halo
!modify starts and ends for certain pes
    if (is.eq.1)     indices(1) = is
    if (ie.eq.npx-1) indices(2) = ie+i_stag
    x_halo_ns = 0
    if (is.eq.1) x_halo_ns=x_halo

!define south root_pe
    is_root_pe = .FALSE.
    if (is.eq.1 .and. js.eq.1) is_root_pe = .TRUE.
!register south halo data in t1
    if (present(var_bc) .and. Atm%neststruct%BCfile_sw_is_open) then
       call register_restart_field(BCfile_sw, trim(var_name)//'_south_t1', &
                                   var_bc%south_t1, indices, global_size, &
                                   x2_pelist, is_root_pe, x_halo=x_halo_ns)
       if (.not. BCfile_sw%is_readonly) then !if writing file
         call register_variable_attribute(BCfile_sw, trim(var_name)//'_south_t1', &
                                          "long_name", trim(var_name)//'_south_t1', &
                                          str_len=len(trim(var_name)//'_south_t1'))
         call register_variable_attribute(BCfile_sw, trim(var_name)//'_south_t1', &
                                          "units", "none", str_len=len("none"))
       endif
    endif
!register south prognostic halo data
    if (present(var) .and. Atm%neststruct%BCfile_sw_is_open) then
       call register_restart_field(BCfile_sw, trim(var_name)//'_south', &
                                   var, indices, global_size, &
                                   x2_pelist, is_root_pe, x_halo=x_halo_ns)
       if (.not. BCfile_sw%is_readonly) then !if writing file
         call register_variable_attribute(BCfile_sw, trim(var_name)//'_south', &
                                          "long_name", trim(var_name)//'_south', &
                                          str_len=len(trim(var_name)//'_south'))
         call register_variable_attribute(BCfile_sw, trim(var_name)//'_south', &
                                          "units", "none", str_len=len("none"))
       endif
    endif

!define north root_pe
    is_root_pe = .FALSE.
    if (ie.eq.npx-1 .and. je.eq.npy-1) is_root_pe = .TRUE.
!register north halo data in t1
    if (present(var_bc) .and. Atm%neststruct%BCfile_ne_is_open) then
       call register_restart_field(BCfile_ne, trim(var_name)//'_north_t1', &
                                   var_bc%north_t1, indices, global_size, &
                                   x1_pelist, is_root_pe, x_halo=x_halo_ns)
       if (.not. BCfile_ne%is_readonly) then !if writing file
         call register_variable_attribute(BCfile_ne, trim(var_name)//'_north_t1', &
                                          "long_name", trim(var_name)//'_north_t1', &
                                          str_len=len(trim(var_name)//'_north_t1'))
         call register_variable_attribute(BCfile_ne, trim(var_name)//'_north_t1', &
                                          "units", "none", str_len=len("none"))
       endif
    endif

!reset indices for prognostic variables in the north halo
    indices(3) = jed-y_halo+1+j_stag
    indices(4) = jed+j_stag
!register north prognostic halo data
    if (present(var) .and. Atm%neststruct%BCfile_ne_is_open) then
       call register_restart_field(BCfile_ne, trim(var_name)//'_north', &
                                   var, indices, global_size, &
                                   x1_pelist, is_root_pe, x_halo=x_halo_ns, &
                                   y_halo=(size(var,2)-y_halo), jshift=-(je+j_stag))
       if (.not. BCfile_ne%is_readonly) then !if writing file
         call register_variable_attribute(BCfile_ne, trim(var_name)//'_north', &
                                          "long_name", trim(var_name)//'_north', &
                                          str_len=len(trim(var_name)//'_north'))
         call register_variable_attribute(BCfile_ne, trim(var_name)//'_north', &
                                          "units", "none", str_len=len("none"))
       endif
    endif

    deallocate (x1_pelist)
    deallocate (y1_pelist)
    deallocate (x2_pelist)
    deallocate (y2_pelist)

  end subroutine register_bcs_2d


  subroutine register_bcs_3d(Atm, BCfile_ne, BCfile_sw, fname_ne, fname_sw, &
                             var_name, var, var_bc, istag, jstag, mandatory)
    type(fv_atmos_type),      intent(in)    :: Atm
    type(FmsNetcdfFile_t),  intent(inout) :: BCfile_ne, BCfile_sw
    character(len=120),       intent(in)    :: fname_ne, fname_sw
    character(len=*),         intent(in)    :: var_name
    real, dimension(:,:,:),   intent(in), optional :: var
    type(fv_nest_BC_type_3D), intent(in), optional :: var_bc
    integer,                  intent(in), optional :: istag, jstag
    logical,                  intent(IN), optional :: mandatory

    integer :: npx, npy, i_stag, j_stag
    integer :: is, ie, js, je, isd, ied, jsd, jed, n
    integer :: x_halo, y_halo, x_halo_ns, id_restart
    integer :: layout(2), global_size(3), indices(4)
    integer, allocatable, dimension(:) :: x1_pelist, y1_pelist
    integer, allocatable, dimension(:) :: x2_pelist, y2_pelist
    logical :: is_root_pe
    logical :: mandatory_flag

    mandatory_flag = .true.
    if (present(mandatory)) mandatory_flag = mandatory

    i_stag = 0
    j_stag = 0
    if (present(istag)) i_stag = istag
    if (present(jstag)) j_stag = jstag
    call mpp_get_global_domain(Atm%domain, xsize = npx, ysize = npy, position=CORNER )
    call mpp_get_data_domain(Atm%domain, isd, ied, jsd, jed )
    call mpp_get_compute_domain(Atm%domain, is, ie, js, je )
    call mpp_get_layout(Atm%domain, layout)
    allocate (x1_pelist(layout(1)))
    allocate (y1_pelist(layout(2)))
    allocate (x2_pelist(layout(1)))
    allocate (y2_pelist(layout(2)))
    x_halo = is-isd
    y_halo = js-jsd
! define west and east pelist
    do n = 1,layout(2)
      y1_pelist(n)=mpp_root_pe()+layout(1)*n-1
      y2_pelist(n)=mpp_root_pe()+layout(1)*(n-1)
    enddo
! define south and north pelist
    do n = 1,layout(1)
      x1_pelist(n)=mpp_root_pe()+layout(1)*(layout(2)-1)+(n-1)
      x2_pelist(n)=mpp_root_pe()+(n-1)
    enddo
! declare the pelists inside of mpp (creates the MPI communicator)
    call mpp_declare_pelist(x1_pelist)
    call mpp_declare_pelist(x2_pelist)
    call mpp_declare_pelist(y1_pelist)
    call mpp_declare_pelist(y2_pelist)

!EAST & WEST
!set defaults for west/east halo regions
    indices(1) = 1
    indices(2) = x_halo
    indices(3) = jsd
    indices(4) = jed + j_stag
    global_size(1) = x_halo
    global_size(2) = npy-1+2*y_halo + j_stag
    global_size(3) = Atm%npz

!define west root_pe
    is_root_pe = .FALSE.
    if (is.eq.1 .and. js.eq.1) is_root_pe = .TRUE.
!register west halo data in t1
    if (present(var_bc) .and. Atm%neststruct%BCfile_sw_is_open) then
       call register_restart_field(BCfile_sw, trim(var_name)//'_west_t1', &
                                   var_bc%west_t1, indices, global_size, &
                                   y2_pelist, is_root_pe, jshift=y_halo, &
                                   is_optional=.not.mandatory_flag)
       if (.not. BCfile_sw%is_readonly) then !if writing file
         call register_variable_attribute(BCfile_sw, trim(var_name)//'_west_t1', &
                                          "long_name", trim(var_name)//'_west_t1', &
                                          str_len=len(trim(var_name)//'_west_t1'))
         call register_variable_attribute(BCfile_sw, trim(var_name)//'_west_t1', &
                                          "units", "none", str_len=len("none"))
       endif
    endif

!register west prognostic halo data
    if (present(var) .and. Atm%neststruct%BCfile_sw_is_open) then
       call register_restart_field(BCfile_sw, trim(var_name)//'_west', &
                                   var, indices, global_size, &
                                   y2_pelist, is_root_pe, jshift=y_halo, &
                                   is_optional=.not.mandatory_flag)
       if (.not. BCfile_sw%is_readonly) then !if writing file
         call register_variable_attribute(BCfile_sw, trim(var_name)//'_west', &
                                          "long_name", trim(var_name)//'_west', &
                                          str_len=len(trim(var_name)//'_west'))
         call register_variable_attribute(BCfile_sw, trim(var_name)//'_west', &
                                          "units", "none", str_len=len("none"))
       endif
    endif

!define east root_pe
    is_root_pe = .FALSE.
    if (ie.eq.npx-1 .and. je.eq.npy-1) is_root_pe = .TRUE.
!register east halo data in t1
    if (present(var_bc) .and. Atm%neststruct%BCfile_ne_is_open) then
       call register_restart_field(BCfile_ne, trim(var_name)//'_east_t1', &
                                   var_bc%east_t1, indices, global_size, &
                                   y1_pelist, is_root_pe, jshift=y_halo, &
                                   is_optional=.not.mandatory_flag)
       if (.not. BCfile_ne%is_readonly) then !if writing file
         call register_variable_attribute(BCfile_ne, trim(var_name)//'_east_t1', &
                                          "long_name", trim(var_name)//'_east_t1', &
                                          str_len=len(trim(var_name)//'_east_t1'))
         call register_variable_attribute(BCfile_ne, trim(var_name)//'_east_t1', &
                                          "units", "none", str_len=len("none"))
       endif
    endif

!reset indices for prognostic variables in the east halo
    indices(1) = ied-x_halo+1+i_stag
    indices(2) = ied+i_stag
!register east prognostic halo data
    if (present(var) .and. Atm%neststruct%BCfile_ne_is_open) then
       call register_restart_field(BCfile_ne, trim(var_name)//'_east', &
                                   var, indices, global_size, &
                                   y1_pelist, is_root_pe, jshift=y_halo, &
                                   x_halo=(size(var,1)-x_halo), ishift=-(ie+i_stag), &
                                   is_optional=.not.mandatory_flag)
       if (.not. BCfile_ne%is_readonly) then !if writing file
         call register_variable_attribute(BCfile_ne, trim(var_name)//'_east', &
                                          "long_name", trim(var_name)//'_east', &
                                          str_len=len(trim(var_name)//'_east'))
         call register_variable_attribute(BCfile_ne, trim(var_name)//'_east', &
                                          "units", "none", str_len=len("none"))
       endif
    endif

!NORTH & SOUTH
!set defaults for north/south halo regions
    indices(1) = isd
    indices(2) = ied+i_stag
    indices(3) = 1
    indices(4) = y_halo
    global_size(1) = npx-1+i_stag
    global_size(2) = y_halo
    global_size(3) = Atm%npz
!modify starts and ends for certain pes
    if (is.eq.1)     indices(1) = is
    if (ie.eq.npx-1) indices(2) = ie+i_stag
    x_halo_ns = 0
    if (is.eq.1) x_halo_ns=x_halo

!define south root_pe
    is_root_pe = .FALSE.
    if (is.eq.1 .and. js.eq.1) is_root_pe = .TRUE.
!register south halo data in t1
    if (present(var_bc) .and. Atm%neststruct%BCfile_sw_is_open) then
       call register_restart_field(BCfile_sw, trim(var_name)//'_south_t1', &
                                   var_bc%south_t1, indices, global_size, &
                                   x2_pelist, is_root_pe, x_halo=x_halo_ns, &
                                   is_optional=.not.mandatory_flag)
       if (.not. BCfile_sw%is_readonly) then !if writing file
         call register_variable_attribute(BCfile_sw, trim(var_name)//'_south_t1', &
                                          "long_name", trim(var_name)//'_south_t1', &
                                          str_len=len(trim(var_name)//'_south_t1'))
         call register_variable_attribute(BCfile_sw, trim(var_name)//'_south_t1', &
                                          "units", "none", str_len=len("none"))
       endif
    endif
!register south prognostic halo data
    if (present(var) .and. Atm%neststruct%BCfile_sw_is_open) then
       call register_restart_field(BCfile_sw, trim(var_name)//'_south', &
                                   var, indices, global_size, &
                                   x2_pelist, is_root_pe, x_halo=x_halo_ns, &
                                   is_optional=.not.mandatory_flag)
       if (.not. BCfile_sw%is_readonly) then !if writing file
         call register_variable_attribute(BCfile_sw, trim(var_name)//'_south', &
                                          "long_name", trim(var_name)//'_south', &
                                          str_len=len(trim(var_name)//'_south'))
         call register_variable_attribute(BCfile_sw, trim(var_name)//'_south', &
                                          "units", "none", str_len=len("none"))
       endif
    endif

!define north root_pe
    is_root_pe = .FALSE.
    if (ie.eq.npx-1 .and. je.eq.npy-1) is_root_pe = .TRUE.
!register north halo data in t1
    if (present(var_bc) .and. Atm%neststruct%BCfile_ne_is_open) then
       call register_restart_field(BCfile_ne, trim(var_name)//'_north_t1', &
                                   var_bc%north_t1, indices, global_size, &
                                   x1_pelist, is_root_pe, x_halo=x_halo_ns, &
                                   is_optional=.not.mandatory_flag)
       if (.not. BCfile_ne%is_readonly) then !if writing file
         call register_variable_attribute(BCfile_ne, trim(var_name)//'_north_t1', &
                                          "long_name", trim(var_name)//'_north_t1', &
                                          str_len=len(trim(var_name)//'_north_t1'))
         call register_variable_attribute(BCfile_ne, trim(var_name)//'_north_t1', &
                                          "units", "none", str_len=len("none"))
       endif
    endif

!reset indices for prognostic variables in the north halo
    indices(3) = jed-y_halo+1+j_stag
    indices(4) = jed+j_stag
!register north prognostic halo data
    if (present(var) .and. Atm%neststruct%BCfile_ne_is_open) then
       call register_restart_field(BCfile_ne, trim(var_name)//'_north', &
                                   var, indices, global_size, &
                                   x1_pelist, is_root_pe, x_halo=x_halo_ns, &
                                   y_halo=(size(var,2)-y_halo), jshift=-(je+j_stag), &
                                   is_optional=.not.mandatory_flag)
       if (.not. BCfile_ne%is_readonly) then !if writing file
         call register_variable_attribute(BCfile_ne, trim(var_name)//'_north', &
                                          "long_name", trim(var_name)//'_north', &
                                          str_len=len(trim(var_name)//'_north'))
         call register_variable_attribute(BCfile_ne, trim(var_name)//'_north', &
                                          "units", "none", str_len=len("none"))
       endif
    endif
    deallocate (x1_pelist)
    deallocate (y1_pelist)
    deallocate (x2_pelist)
    deallocate (y2_pelist)

  end subroutine register_bcs_3d


  ! </SUBROUTINE> NAME="fv_io_register_restart_BCs"
  !#####################################################################

  subroutine fv_io_register_restart_BCs(Atm)
    type(fv_atmos_type),        intent(inout) :: Atm

    integer :: n, ntracers, ntprog, ntdiag
    character(len=120) :: tname, fname_ne, fname_sw

    fname_ne = 'fv_BC_ne.res.nc'
    fname_sw = 'fv_BC_sw.res.nc'

    ntprog=size(Atm%q,4)
    ntdiag=size(Atm%qdiag,4)
    ntracers=ntprog+ntdiag

    call register_bcs_2d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                         fname_ne, fname_sw, 'phis', var=Atm%phis)
    call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                         fname_ne, fname_sw, 'delp', Atm%delp, Atm%neststruct%delp_BC)
    do n=1,ntprog
       call get_tracer_names(MODEL_ATMOS, n, tname)
       call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                            fname_ne, fname_sw, trim(tname), Atm%q(:,:,:,n), Atm%neststruct%q_BC(n), mandatory=.false.)
    enddo
    do n=ntprog+1,ntracers
       call get_tracer_names(MODEL_ATMOS, n, tname)
       call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                            fname_ne, fname_sw, trim(tname), var=Atm%qdiag(:,:,:,n), mandatory=.false.)
    enddo
#ifndef SW_DYNAMICS
    call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                         fname_ne, fname_sw, 'pt', Atm%pt, Atm%neststruct%pt_BC)
    if ((.not.Atm%flagstruct%hydrostatic)) then
       if (is_master()) print*, 'fv_io_register_restart_BCs: REGISTERING NH BCs'
      call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                           fname_ne, fname_sw, 'w', Atm%w, Atm%neststruct%w_BC, mandatory=.false.)
      call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                           fname_ne, fname_sw, 'delz', var_bc=Atm%neststruct%delz_BC, mandatory=.false.)
!                           fname_ne, fname_sw, 'delz', Atm%delz, Atm%neststruct%delz_BC, mandatory=.false.)
    endif
#ifdef USE_COND
       call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                            fname_ne, fname_sw,'q_con', var_bc=Atm%neststruct%q_con_BC, mandatory=.false.)
#ifdef MOIST_CAPPA
       call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
            fname_ne, fname_sw, 'cappa', var_bc=Atm%neststruct%cappa_BC, mandatory=.false.)
#endif
#endif
#endif
    if (Atm%flagstruct%is_ideal_case) then
       call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                            fname_ne, fname_sw, 'u0', Atm%u0, Atm%neststruct%u_BC, jstag=1)
       call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                            fname_ne, fname_sw, 'v0', Atm%v0, Atm%neststruct%v_BC, istag=1)
    endif
    call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                         fname_ne, fname_sw, 'u', Atm%u, Atm%neststruct%u_BC, jstag=1)
    call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                         fname_ne, fname_sw, 'v', Atm%v, Atm%neststruct%v_BC, istag=1)
    call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                         fname_ne, fname_sw, 'uc', var_bc=Atm%neststruct%uc_BC, istag=1)
    call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                         fname_ne, fname_sw, 'vc', var_bc=Atm%neststruct%vc_BC, jstag=1)
    call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                         fname_ne, fname_sw, 'divg', var_bc=Atm%neststruct%divg_BC, istag=1,jstag=1, mandatory=.false.)

    return
  end subroutine fv_io_register_restart_BCs


  subroutine fv_io_write_BCs(Atm, timestamp)
    type(fv_atmos_type), intent(inout)        :: Atm
    character(len=*),    intent(in), optional :: timestamp
    integer, allocatable, dimension(:)        :: all_pelist
    integer, dimension(2)                     :: layout
    integer                                   :: n
    character(len=1)                   :: tile_num
    character(len=120)                 :: fname_ne, fname_sw

    if (present(timestamp)) then
      fname_ne = 'RESTART/'//trim(timestamp)//'fv_BC_ne.res.nc'
      fname_sw = 'RESTART/'//trim(timestamp)//'fv_BC_sw.res.nc'
    else
      fname_ne = 'RESTART/fv_BC_ne.res.nc'
      fname_sw = 'RESTART/fv_BC_sw.res.nc'
    endif

    allocate(all_pelist(mpp_npes()))
    call mpp_get_current_pelist(all_pelist)

    Atm%neststruct%BCfile_sw_is_open = open_file(Atm%neststruct%BCfile_sw, fname_sw, "overwrite", is_restart=.true., pelist=all_pelist)
    Atm%neststruct%BCfile_ne_is_open = open_file(Atm%neststruct%BCfile_ne, fname_ne, "overwrite", is_restart=.true., pelist=all_pelist)
    call fv_io_register_restart_BCs(Atm)

    if (Atm%neststruct%BCfile_sw_is_open) then
      call write_restart_bc(Atm%neststruct%BCfile_sw)
      call close_file(Atm%neststruct%BCfile_sw)
      Atm%neststruct%BCfile_sw_is_open = .false.
    endif

    if (Atm%neststruct%BCfile_ne_is_open) then
     call write_restart_bc(Atm%neststruct%BCfile_ne)
     call close_file(Atm%neststruct%BCfile_ne)
     Atm%neststruct%BCfile_ne_is_open = .false.
    endif

    deallocate(all_pelist)

    return
  end subroutine fv_io_write_BCs


  subroutine fv_io_read_BCs(Atm)
    type(fv_atmos_type), intent(inout) :: Atm
    integer, allocatable, dimension(:) :: all_pelist
    integer, dimension(2)              :: layout
    integer                            :: n
    character(len=1)                   :: tile_num
    character(len=120)                 :: fname_ne, fname_sw

    fname_ne = 'INPUT/fv_BC_ne.res.nc'
    fname_sw = 'INPUT/fv_BC_sw.res.nc'

    allocate(all_pelist(mpp_npes()))
    call mpp_get_current_pelist(all_pelist)

    Atm%neststruct%BCfile_sw_is_open = open_file(Atm%neststruct%BCfile_sw, fname_sw, "read", is_restart=.true., pelist=all_pelist)
    Atm%neststruct%BCfile_ne_is_open = open_file(Atm%neststruct%BCfile_ne, fname_ne, "read", is_restart=.true., pelist=all_pelist)
    call fv_io_register_restart_BCs(Atm)

    if (Atm%neststruct%BCfile_sw_is_open) then
      call read_restart_bc(Atm%neststruct%BCfile_sw, ignore_checksum=Atm%flagstruct%ignore_rst_cksum)
      call close_file(Atm%neststruct%BCfile_sw)
      Atm%neststruct%BCfile_sw_is_open = .false.
    endif

    if (Atm%neststruct%BCfile_ne_is_open) then
     call read_restart_bc(Atm%neststruct%BCfile_ne, ignore_checksum=Atm%flagstruct%ignore_rst_cksum)
     call close_file(Atm%neststruct%BCfile_ne)
     Atm%neststruct%BCfile_ne_is_open = .false.
    endif


    !These do not work yet
    !need to modify register_bcs_?d to get ids for registered variables, and then use query_initialized_id
    !Atm%neststruct%divg_BC%initialized = field_exist(fname_ne, 'divg_north_t1', Atm%domain)
    !Atm%neststruct%w_BC%initialized    = field_exist(fname_ne, 'w_north_t1', Atm%domain)
    !Atm%neststruct%delz_BC%initialized = field_exist(fname_ne, 'delz_north_t1', Atm%domain)
    !if (is_master()) print*, ' BCs: ', Atm%neststruct%divg_BC%initialized, Atm%neststruct%w_BC%initialized, Atm%neststruct%delz_BC%initialized

    deallocate(all_pelist)

    return
  end subroutine fv_io_read_BCs

end module fv_io_mod
