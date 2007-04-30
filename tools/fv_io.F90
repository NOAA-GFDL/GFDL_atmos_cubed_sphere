module fv_io_mod
  !-----------------------------------------------------------------------
  !                   GNU General Public License                        
  !                                                                      
  ! This program is free software; you can redistribute it and/or modify it and  
  ! are expected to follow the terms of the GNU General Public License  
  ! as published by the Free Software Foundation; either version 2 of   
  ! the License, or (at your option) any later version.                 
  !                                                                      
  ! MOM is distributed in the hope that it will be useful, but WITHOUT    
  ! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  
  ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
  ! License for more details.                                           
  !                                                                      
  ! For the full text of the GNU General Public License,                
  ! write to: Free Software Foundation, Inc.,                           
  !           675 Mass Ave, Cambridge, MA 02139, USA.                   
  ! or see:   http://www.gnu.org/licenses/gpl.html                      
  !-----------------------------------------------------------------------
  ! 
  ! <CONTACT EMAIL= "Jeffrey.Durachta@noaa.gov">Jeffrey Durachta </CONTACT>

  ! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

  !<OVERVIEW>
  ! Restart facilities for FV core
  !</OVERVIEW>
  !<DESCRIPTION>
  ! This module writes and reads restart files for the FV core. Additionally
  ! it provides setup and calls routines necessary to provide a complete restart
  ! for the model.
  !</DESCRIPTION>

  use fms_mod,                 only: file_exist
  use fms_mod,                 only: read_data, write_data
  use fms_io_mod,              only: fms_io_exit, get_tile_string
  use fv_arrays_mod,           only: fv_atmos_type
  use mpp_mod,                 only: mpp_error, FATAL, NOTE
  use mpp_domains_mod,         only: domain2d
  use mpp_domains_mod,         only: EAST, NORTH, mpp_get_tile_id
  use mpp_domains_mod,         only: mpp_get_compute_domain, mpp_get_data_domain
! use tracer_manager_mod,      only: get_tracer_names

  implicit none
  private

  public :: fv_io_init, fv_io_exit, fv_io_read_restart, fv_io_write_restart

  logical                       :: module_is_initialized = .FALSE.

  !--- version information variables ----
  character(len=128) :: version = '$Id: fv_io.F90,v 14.0 2007/03/15 21:58:41 fms Exp $'
  character(len=128) :: tagname = '$Name: nalanda_2007_04 $'

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
  ! <SUBROUTINE NAME="fv_io_read_restart">
  !
  ! <DESCRIPTION>
  ! Write the fv core restart quantities 
  ! </DESCRIPTION>
  subroutine  fv_io_read_restart(fv_domain,Atm)
    type(domain2d),      intent(inout) :: fv_domain
    type(fv_atmos_type), intent(inout) :: Atm(:)

    character(len=64)    :: fname, fname_nd, tracer_name
    integer              :: isc, iec, jsc, jec, n, nt, nk, ntracers
    integer              :: isd, ied, jsd, jed
    integer              :: ntileMe
    integer, allocatable :: tile_id(:)

    ntileMe = size(Atm(:))  ! This will have to be modified for mult tiles per PE
    allocate(tile_id(ntileMe))
    tile_id = mpp_get_tile_id(fv_domain)
 
!   call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)
    ntracers = size(Atm(1)%q,4)  ! Temporary until we get tracer manager integrated
 
    fname_nd = 'INPUT/fv_core.res.nc'
 
  ! write_data does not (yet?) support vector data and tiles
    call read_data(fname_nd, 'ak', Atm(1)%ak(:))
    call read_data(fname_nd, 'bk', Atm(1)%bk(:))
 
    do n = 1, ntileMe
       isc = Atm(n)%isc; iec = Atm(n)%iec; jsc = Atm(n)%jsc; jec = Atm(n)%jec
       call get_tile_string(fname, 'INPUT/fv_core.res.tile', tile_id(n), '.nc' )
       if(file_exist(fname))then
         call read_data(fname, 'u', Atm(n)%u(isc:iec,jsc:jec+1,:), domain=fv_domain, position=NORTH,tile_count=n)
         call read_data(fname, 'v', Atm(n)%v(isc:iec+1,jsc:jec,:), domain=fv_domain, position=EAST,tile_count=n)

         if ( (.not.Atm(n)%hydrostatic) .and. (.not.Atm(n)%Make_NH) ) then
              call read_data(fname, 'W',     Atm(n)%w(isc:iec,jsc:jec,:), domain=fv_domain, tile_count=n)
              call read_data(fname, 'DZ', Atm(n)%delz(isc:iec,jsc:jec,:), domain=fv_domain, tile_count=n)
              if ( Atm(n)%hybrid_z )   &
              call read_data(fname, 'ZE0', Atm(n)%ze0(isc:iec,jsc:jec,:), domain=fv_domain, tile_count=n)
         endif

         call read_data(fname, 'T', Atm(n)%pt(isc:iec,jsc:jec,:), domain=fv_domain, tile_count=n)
         call read_data(fname, 'delp', Atm(n)%delp(isc:iec,jsc:jec,:), domain=fv_domain, tile_count=n)
         call read_data(fname, 'phis', Atm(n)%phis(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
       else
         call mpp_error(FATAL,'==> Error from fv_read_restart: Expected file '//trim(fname)//' does not exist')
       endif

       call get_tile_string(fname, 'INPUT/fv_srf_wnd.res.tile', tile_id(n), '.nc' )
       if(file_exist(fname))then
         call read_data(fname, 'u_srf', Atm(n)%u_srf(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
         call read_data(fname, 'v_srf', Atm(n)%v_srf(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
         Atm(n)%srf_init = .true.
       else
         call mpp_error(NOTE,'==> Warning from fv_read_restart: Expected file '//trim(fname)//' does not exist')
         Atm(n)%srf_init = .false.
       endif

#ifdef FV_LAND
!----------------------------------------------------------------------------------------------------------------
! Optional terrain deviation (sgh) and land fraction (oro)
       call get_tile_string(fname, 'INPUT/fv_land.res.tile', tile_id(n), '.nc' )
       if(file_exist(fname))then
         call read_data(fname, 'sgh', Atm(n)%sgh(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
         call read_data(fname, 'oro', Atm(n)%oro(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
       else
         call mpp_error(NOTE,'==> Warning from fv_read_restart: Expected file '//trim(fname)//' does not exist')
       endif
!----------------------------------------------------------------------------------------------------------------
#endif
 
       call get_tile_string(fname, 'INPUT/fv_tracer.res.tile', tile_id(n), '.nc' )
       if(file_exist(fname))then
         do nt = 1, ntracers
!          call get_tracer_names(MODEL_ATMOS, n, tracer_name)
           ! Temporary until we get tracer manager (or at least a tracer list) integrated
           call get_tile_string(tracer_name, 'atm_T', nt)
           call read_data(fname, tracer_name, Atm(n)%q(isc:iec,jsc:jec,:,nt), domain=fv_domain, tile_count=n)
         end do
       else
         if (ntracers>0) call mpp_error(FATAL,'==> Error from fv_read_restart: Expected file '//trim(fname)//' does not exist')
       endif
 
    end do
 
    deallocate(tile_id)

  end subroutine  fv_io_read_restart
  ! </SUBROUTINE> NAME="fv_io_read_restart"
  !#####################################################################


  !#####################################################################
  ! <SUBROUTINE NAME="fv_io_write_restart">
  !
  ! <DESCRIPTION>
  ! Write the fv core restart quantities 
  ! </DESCRIPTION>
  subroutine  fv_io_write_restart(fv_domain,Atm)
    type(domain2d),      intent(in) :: fv_domain
    type(fv_atmos_type), intent(in) :: Atm(:)

    character(len=64)    :: fname, fname_nd, tracer_name
    integer              :: isc, iec, jsc, jec, n, nt, ntracers
    integer              :: ntileMe
    integer, allocatable :: tile_id(:)

    ntileMe = size(Atm(:))  ! This will need mods for more than 1 tile per pe
    allocate(tile_id(ntileMe))
    tile_id = mpp_get_tile_id(fv_domain)
 
!   call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)
    ntracers = size(Atm(1)%q,4)  ! Temporary until we get tracer manager integrated
 
    fname_nd = 'RESTART/fv_core.res.nc'
 
  ! write_data does not (yet?) support vector data and tiles
    call write_data(fname_nd, 'ak', Atm(1)%ak(:))
    call write_data(fname_nd, 'bk', Atm(1)%bk(:))
 
    do n = 1, ntileMe
       isc = Atm(n)%isc; iec = Atm(n)%iec; jsc = Atm(n)%jsc; jec = Atm(n)%jec  
       call get_tile_string(fname, 'RESTART/fv_core.res.tile', tile_id(n), '.nc' )
       call write_data(fname, 'u', Atm(n)%u(isc:iec,jsc:jec+1,:), domain=fv_domain, position=NORTH,tile_count=n)
       call write_data(fname, 'v', Atm(n)%v(isc:iec+1,jsc:jec,:), domain=fv_domain, position=EAST,tile_count=n)

       if ( .not.Atm(n)%hydrostatic ) then
            call write_data(fname, 'W',     Atm(n)%w(isc:iec,jsc:jec,:), domain=fv_domain, tile_count=n)
            call write_data(fname, 'DZ', Atm(n)%delz(isc:iec,jsc:jec,:), domain=fv_domain, tile_count=n)
            if ( Atm(n)%hybrid_z )  &
            call write_data(fname, 'ZE0', Atm(n)%ze0(isc:iec,jsc:jec,:), domain=fv_domain, tile_count=n)
       endif

       call write_data(fname, 'T', Atm(n)%pt(isc:iec,jsc:jec,:), domain=fv_domain, tile_count=n)
       call write_data(fname, 'delp', Atm(n)%delp(isc:iec,jsc:jec,:), domain=fv_domain, tile_count=n)
       call write_data(fname, 'phis', Atm(n)%phis(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
 
       call get_tile_string(fname, 'RESTART/fv_srf_wnd.res.tile', tile_id(n), '.nc' )
       call write_data(fname, 'u_srf', Atm(n)%u_srf(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
       call write_data(fname, 'v_srf', Atm(n)%v_srf(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)

#ifdef FV_LAND
       call get_tile_string(fname, 'RESTART/fv_land.res.tile', tile_id(n), '.nc' )
       call write_data(fname, 'sgh', Atm(n)%sgh(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
       call write_data(fname, 'oro', Atm(n)%oro(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
#endif

       call get_tile_string(fname, 'RESTART/fv_tracer.res.tile', tile_id(n), '.nc' )
       do nt = 1, ntracers
!        call get_tracer_names(MODEL_ATMOS, n, tracer_name)
         ! Temporary until we get tracer manager (or at least a tracer list) integrated
         call get_tile_string(tracer_name, 'atm_T', nt)
         call write_data(fname, tracer_name, Atm(n)%q(isc:iec,jsc:jec,:,nt), domain=fv_domain, tile_count=n)
       end do
    end do
 
    module_is_initialized = .false.

  end subroutine  fv_io_write_restart
  ! </SUBROUTINE> NAME="fv_io_write_restart"
  !#####################################################################

end module fv_io_mod
