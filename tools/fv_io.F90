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
  use tracer_manager_mod,      only: tr_get_tracer_names=>get_tracer_names, &
                                     get_tracer_names, get_number_tracers, &
                                     set_tracer_profile, &
                                     get_tracer_index, NO_TRACER
  use field_manager_mod,       only: MODEL_ATMOS  
  use fms_mod,                 only: field_exist    


  implicit none
  private

  public :: fv_io_init, fv_io_exit, fv_io_read_restart, remap_restart, fv_io_write_restart

  logical                       :: module_is_initialized = .FALSE.

  !--- version information variables ----
  character(len=128) :: version = '$Id: fv_io.F90,v 1.1.4.3.2.2.2.2.2.10.2.2.2.1.2.2 2007/11/01 18:35:21 sjl Exp $'
  character(len=128) :: tagname = '$Name: omsk_2008_03 $'

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
    integer              :: ntileMe
    integer, allocatable :: tile_id(:)

    character(len=128)           :: tracer_longname, tracer_units

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

  if ( Atm(n)%fv_land ) then
!----------------------------------------------------------------------------------------------------------------
! Optional terrain deviation (sgh) and land fraction (oro)
       call get_tile_string(fname, 'INPUT/mg_drag.res.tile', tile_id(n), '.nc' )
       if(file_exist(fname))then
         call read_data(fname, 'ghprime', Atm(n)%sgh(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
       else
         call mpp_error(NOTE,'==> Warning from fv_read_restart: Expected file '//trim(fname)//' does not exist')
       endif
! Land-water mask:
       call get_tile_string(fname, 'INPUT/fv_land.res.tile', tile_id(n), '.nc' )
       if(file_exist(fname))then
         call read_data(fname, 'oro', Atm(n)%oro(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
       else
         call mpp_error(NOTE,'==> Warning from fv_read_restart: Expected file '//trim(fname)//' does not exist')
       endif
!----------------------------------------------------------------------------------------------------------------
  endif
       call get_tile_string(fname, 'INPUT/fv_tracer.res.tile', tile_id(n), '.nc' )

         DO nt = 1, ntracers
           call get_tracer_names(MODEL_ATMOS, nt, tracer_name)

           if (file_exist(fname)) then
              if (field_exist(fname,tracer_name)) then
                 call read_data(fname, tracer_name, Atm(n)%q(isc:iec,jsc:jec,:,nt), &
                                                      domain=fv_domain, tile_count=n)
                 call mpp_error(NOTE,'==>  Have read tracer '//trim(tracer_name)//' from fv_tracer.res')
                 cycle
              endif
           endif

           call set_tracer_profile (MODEL_ATMOS, nt, Atm(n)%q(isc:iec,jsc:jec,:,nt)  )
           call mpp_error(NOTE,'==>  Setting tracer '//trim(tracer_name)//' from set_tracer')
         ENDDO

    end do
 
    deallocate(tile_id)

  end subroutine  fv_io_read_restart
  ! </SUBROUTINE> NAME="fv_io_read_restart"
  !#####################################################################


  subroutine  remap_restart(fv_domain,Atm)
  use fv_mapz_mod,       only: rst_remap

    type(domain2d),      intent(inout) :: fv_domain
    type(fv_atmos_type), intent(inout) :: Atm(:)

    character(len=64)    :: fname, fname_nd, tracer_name
    integer              :: isc, iec, jsc, jec, n, nt, nk, ntracers
    integer              :: isd, ied, jsd, jed
    integer              :: ntileMe
    integer, allocatable :: tile_id(:)
!
!-------------------------------------------------------------------------
    real, allocatable:: ak_r(:), bk_r(:)
    real, allocatable:: u_r(:,:,:), v_r(:,:,:), pt_r(:,:,:), delp_r(:,:,:)
    real, allocatable:: w_r(:,:,:), delz_r(:,:,:), ze0_r(:,:,:)
    real, allocatable:: q_r(:,:,:,:)
!-------------------------------------------------------------------------
    integer npz, npz_rst, ng

    npz     = Atm(1)%npz       ! run time z dimension
    npz_rst = Atm(1)%npz_rst   ! restart z dimension
    isc = Atm(1)%isc; iec = Atm(1)%iec; jsc = Atm(1)%jsc; jec = Atm(1)%jec
    ng = Atm(1)%ng

    isd = isc - ng;  ied = iec + ng
    jsd = jsc - ng;  jed = jec + ng


!   call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)
    ntracers = size(Atm(1)%q,4)  ! Temporary until we get tracer manager integrated

    ntileMe = size(Atm(:))  ! This will have to be modified for mult tiles per PE
    allocate(tile_id(ntileMe))
    tile_id = mpp_get_tile_id(fv_domain)


! Allocate arrays for reading old restart file:
    allocate ( ak_r(npz_rst+1) )
    allocate ( bk_r(npz_rst+1) )

    allocate ( u_r(isc:iec,  jsc:jec+1,npz_rst) )
    allocate ( v_r(isc:iec+1,jsc:jec  ,npz_rst) )

    allocate (   pt_r(isc:iec, jsc:jec,  npz_rst) )
    allocate ( delp_r(isc:iec, jsc:jec,  npz_rst) )
    allocate (    q_r(isc:iec, jsc:jec,  npz_rst, ntracers) )

    if ( (.not.Atm(1)%hydrostatic) .and. (.not.Atm(1)%Make_NH) ) then
           allocate (    w_r(isc:iec, jsc:jec,  npz_rst) )
           allocate ( delz_r(isc:iec, jsc:jec,  npz_rst) )
           if ( Atm(1)%hybrid_z )   &
           allocate ( ze0_r(isc:iec, jsc:jec,  npz_rst+1) )
    endif

    fname_nd = 'INPUT/fv_core.res.nc'

  ! write_data does not (yet?) support vector data and tiles
    call read_data(fname_nd, 'ak', ak_r(1:npz_rst+1))
    call read_data(fname_nd, 'bk', bk_r(1:npz_rst+1))

    do n = 1, ntileMe
       call get_tile_string(fname, 'INPUT/fv_core.res.tile', tile_id(n), '.nc' )
       if(file_exist(fname))then
         call read_data(fname, 'u', u_r(isc:iec,jsc:jec+1,:), domain=fv_domain, position=NORTH,tile_count=n)
         call read_data(fname, 'v', v_r(isc:iec+1,jsc:jec,:), domain=fv_domain, position=EAST,tile_count=n)

         if ( (.not.Atm(n)%hydrostatic) .and. (.not.Atm(n)%Make_NH) ) then
              call read_data(fname, 'W',     w_r(isc:iec,jsc:jec,:), domain=fv_domain, tile_count=n)
              call read_data(fname, 'DZ', delz_r(isc:iec,jsc:jec,:), domain=fv_domain, tile_count=n)
              if ( Atm(n)%hybrid_z )   &
              call read_data(fname, 'ZE0', ze0_r(isc:iec,jsc:jec,:), domain=fv_domain, tile_count=n)
         endif

         call read_data(fname, 'T', pt_r(isc:iec,jsc:jec,:), domain=fv_domain, tile_count=n)
         call read_data(fname, 'delp', delp_r(isc:iec,jsc:jec,:), domain=fv_domain, tile_count=n)
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

     if ( Atm(n)%fv_land ) then
! Optional terrain deviation (sgh)
       call get_tile_string(fname, 'INPUT/mg_drag.res.tile', tile_id(n), '.nc' )
       if(file_exist(fname))then
         call read_data(fname, 'ghprime', Atm(n)%sgh(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
       else
         call mpp_error(NOTE,'==> Warning from fv_read_restart: Expected file '//trim(fname)//' does not exist')
       endif
! Land-water mask
       call get_tile_string(fname, 'INPUT/fv_land.res.tile', tile_id(n), '.nc' )
       if(file_exist(fname))then
         call read_data(fname, 'oro', Atm(n)%oro(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
       else
         call mpp_error(NOTE,'==> Warning from fv_read_restart: Expected file '//trim(fname)//' does not exist')
       endif
     endif

       call get_tile_string(fname, 'INPUT/fv_tracer.res.tile', tile_id(n), '.nc' )

       do nt = 1, ntracers
           call get_tracer_names(MODEL_ATMOS, n, tracer_name)

          if(file_exist(fname))then
              if (field_exist(fname,tracer_name)) then
                 call read_data(fname, tracer_name, q_r(isc:iec,jsc:jec,:,nt), domain=fv_domain, tile_count=n)
                 call mpp_error(NOTE,'==>  Have read tracer '//trim(tracer_name)//' from fv_tracer.res')
                 cycle
              endif
          endif

          call set_tracer_profile (MODEL_ATMOS, nt, q_r(isc:iec,jsc:jec,:,nt)  )
          call mpp_error(NOTE,'==>  Setting tracer '//trim(tracer_name)//' from set_tracer')
       enddo

       call rst_remap(npz_rst, npz, isc, iec, jsc, jec, isd, ied, jsd, jed, ntracers,              &
                      delp_r,      u_r,      v_r,      w_r,      delz_r,      pt_r,      q_r,      &
                      Atm(n)%delp, Atm(n)%u, Atm(n)%v, Atm(n)%w, Atm(n)%delz, Atm(n)%pt, Atm(n)%q, &
                      ak_r,  bk_r, Atm(n)%ak, Atm(n)%bk, Atm(n)%hydrostatic)
    end do

    deallocate(tile_id)
    deallocate( ak_r )
    deallocate( bk_r )
    deallocate( u_r )
    deallocate( v_r )
    deallocate( pt_r )
    deallocate( delp_r )
    deallocate( q_r )

    if ( (.not.Atm(1)%hydrostatic) .and. (.not.Atm(1)%Make_NH) ) then
         deallocate ( w_r )
         deallocate ( delz_r )
         if ( Atm(1)%hybrid_z ) deallocate ( ze0_r )
    endif

  end subroutine  remap_restart


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
    character(len=128)   :: tracer_longname, tracer_units

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

     if ( Atm(n)%fv_land ) then
       call get_tile_string(fname, 'RESTART/mg_drag.res.tile', tile_id(n), '.nc' )
       call write_data(fname, 'ghprime', Atm(n)%sgh(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
       call get_tile_string(fname, 'RESTART/fv_land.res.tile', tile_id(n), '.nc' )
       call write_data(fname, 'oro', Atm(n)%oro(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
     endif

       call get_tile_string(fname, 'RESTART/fv_tracer.res.tile', tile_id(n), '.nc' )

       do nt = 1, ntracers
         call tr_get_tracer_names(MODEL_ATMOS, nt, &
                tracer_name, tracer_longname, tracer_units)
         call write_data(fname, tracer_name, Atm(n)%q(isc:iec,jsc:jec,:,nt), domain=fv_domain, tile_count=n)
       end do

    end do

    module_is_initialized = .false.

  end subroutine  fv_io_write_restart
  ! </SUBROUTINE> NAME="fv_io_write_restart"
  !#####################################################################

end module fv_io_mod
