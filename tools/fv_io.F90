! These routines may need some thought before they work with multiple grids. For now, just outputting one grid

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
  use fms_io_mod,              only: fms_io_exit, get_tile_string, &
                                     restart_file_type, register_restart_field, &
                                     save_restart, restore_state, &
                                     set_domain, nullify_domain, &
                                     get_mosaic_tile_file, get_instance_filename, & 
                                     save_restart_border, restore_state_border
  use mpp_mod,                 only: mpp_error, FATAL, NOTE, WARNING, mpp_root_pe, &
                                     mpp_sync, mpp_pe, mpp_declare_pelist
  use mpp_domains_mod,         only: domain2d, EAST, WEST, NORTH, CENTER, SOUTH, CORNER, &
                                     mpp_get_compute_domain, mpp_get_data_domain, & 
                                     mpp_get_layout, mpp_get_ntile_count, &
                                     mpp_get_global_domain
  use tracer_manager_mod,      only: tr_get_tracer_names=>get_tracer_names, &
                                     get_tracer_names, get_number_tracers, &
                                     set_tracer_profile, &
                                     get_tracer_index
  use field_manager_mod,       only: MODEL_ATMOS  
  use external_sst_mod,        only: sst_ncep, sst_anom, use_ncep_sst
  use fv_arrays_mod,           only: fv_atmos_type, fv_nest_BC_type_3D
  use fv_eta_mod,              only: set_eta

  use fv_mp_mod,               only: ng, mp_gather, is_master
  use fms_io_mod,              only: set_domain

  implicit none
  private

  public :: fv_io_init, fv_io_exit, fv_io_read_restart, remap_restart, fv_io_write_restart
  public :: fv_io_read_tracers, fv_io_register_restart, fv_io_register_nudge_restart
  public :: fv_io_register_restart_BCs, fv_io_write_BCs, fv_io_read_BCs

  logical                       :: module_is_initialized = .FALSE.


!---- version number -----
  character(len=128) :: version = '$Id: fv_io.F90,v 20.0.2.3 2014/03/06 17:53:38 Rusty.Benson Exp $'
  character(len=128) :: tagname = '$Name: tikal_201403 $'

  integer ::grid_xtdimid, grid_ytdimid, haloid, pfullid !For writing BCs
  integer ::grid_xtstagdimid, grid_ytstagdimid, oneid

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
    character(len=3)  :: gn
    integer              :: isc, iec, jsc, jec, n, nt, nk, ntracers
    integer              :: ntileMe
    integer              :: ks, ntiles
    real                 :: ptop

    character(len=128)           :: tracer_longname, tracer_units

    if (Atm(1)%grid_number > 1) then
       write(gn,'(A2, I1)') "_g", Atm(1)%grid_number
    else
       gn = ''
    end if

    ntileMe = size(Atm(:))  ! This will need mods for more than 1 tile per pe
    
    call set_domain(Atm(1)%domain)

    if (.not. Atm(1)%neststruct%nested) then
       call restore_state(Atm(1)%Fv_restart)
    endif

    if ( use_ncep_sst .or. Atm(1)%flagstruct%nudge .or. Atm(1)%flagstruct%ncep_ic ) then
       call mpp_error(NOTE, 'READING FROM SST_RESTART DISABLED')
       !call restore_state(Atm(1)%SST_restart)
    endif

    call nullify_domain
 
    do n = 1, ntileMe
!    n = 1
       call set_domain(Atm(n)%domain)
       call restore_state(Atm(n)%Fv_tile_restart)
       call restore_state(Atm(n)%Rsf_restart)
       fname = 'INPUT/fv_srf_wnd'//trim(gn)//'.res.nc'
       if (.not.Atm(n)%neststruct%nested) then
         call get_instance_filename(fname, fname_nd)
         call get_mosaic_tile_file(fname_nd, fname, .FALSE., fv_domain, n)
       end if
       if (file_exist(fname)) then
         Atm(n)%flagstruct%srf_init = .true.
       else
         call mpp_error(NOTE,'==> Warning from fv_read_restart: Expected file '//trim(fname)//' does not exist')
         Atm(n)%flagstruct%srf_init = .false.
       endif

       if ( Atm(n)%flagstruct%fv_land ) then
          call restore_state(Atm(n)%Mg_restart)
          call restore_state(Atm(n)%Lnd_restart)
       endif

       call restore_state(Atm(n)%Tra_restart)

       call nullify_domain

    end do

    return

  end subroutine  fv_io_read_restart
  ! </SUBROUTINE> NAME="fv_io_read_restart"
  !#####################################################################


  subroutine fv_io_read_tracers(fv_domain,Atm)
    type(domain2d),      intent(inout) :: fv_domain
    type(fv_atmos_type), intent(inout) :: Atm(:)
    integer :: n, ntracers, nt, isc, iec, jsc, jec, id_restart
    character(len=3) :: gn
    character(len=64):: fname, fname_nd, tracer_name
    type(restart_file_type) :: Tra_restart_r

    n = 1
    isc = Atm(n)%bd%isc
    iec = Atm(n)%bd%iec
    jsc = Atm(n)%bd%jsc
    jec = Atm(n)%bd%jec
    call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)

    if (Atm(n)%grid_number > 1) then
       write(gn,'(A2, I1)') "_g", Atm(1)%grid_number
    else
       gn = ''
    end if

    call set_domain(fv_domain)
    fname = 'fv_tracer'//trim(gn)//'.res.nc'
    if (.not.Atm(n)%neststruct%nested) then
       call get_instance_filename(fname, fname_nd)
       call get_mosaic_tile_file(fname_nd, fname, .FALSE., fv_domain, n)
    end if
    do nt = 2, ntracers
       call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
       call set_tracer_profile (MODEL_ATMOS, nt, Atm(n)%q(isc:iec,jsc:jec,:,nt)  )
       id_restart = register_restart_field(Tra_restart_r, fname_nd, tracer_name, Atm(n)%q(:,:,:,nt), &
                    domain=fv_domain, mandatory=.false., tile_count=n)
    enddo
    if (file_exist(fname_nd)) call restore_state(Tra_restart_r)
    call nullify_domain()

    return

  end subroutine  fv_io_read_tracers


  subroutine  remap_restart(fv_domain,Atm)
  use fv_mapz_mod,       only: rst_remap

    type(domain2d),      intent(inout) :: fv_domain
    type(fv_atmos_type), intent(inout) :: Atm(:)

    character(len=64)    :: fname, fname_nd, tracer_name
    character(len=3)  :: gn
    integer              :: isc, iec, jsc, jec, n, nt, nk, ntracers
    integer              :: isd, ied, jsd, jed
!    integer              :: ntileMe
    type(restart_file_type) :: FV_restart_r, FV_tile_restart_r, Tra_restart_r
    integer :: id_restart


!
!-------------------------------------------------------------------------
    real, allocatable:: ak_r(:), bk_r(:)
    real, allocatable:: u_r(:,:,:), v_r(:,:,:), pt_r(:,:,:), delp_r(:,:,:)
    real, allocatable:: w_r(:,:,:), delz_r(:,:,:), ze0_r(:,:,:)
    real, allocatable:: q_r(:,:,:,:)
!-------------------------------------------------------------------------
    integer npz, npz_rst, ng

    if (Atm(1)%grid_number > 1) then
       write(gn,'(A2, I1)') "_g", Atm(1)%grid_number
    else
       gn = ''
    end if

    npz     = Atm(1)%npz       ! run time z dimension
    npz_rst = Atm(1)%flagstruct%npz_rst   ! restart z dimension
    isc = Atm(1)%bd%isc; iec = Atm(1)%bd%iec; jsc = Atm(1)%bd%jsc; jec = Atm(1)%bd%jec
    ng = Atm(1)%ng

    isd = isc - ng;  ied = iec + ng
    jsd = jsc - ng;  jed = jec + ng


!   call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)
    ntracers = size(Atm(1)%q,4)  ! Temporary until we get tracer manager integrated

!    ntileMe = size(Atm(:))  ! This will have to be modified for mult tiles per PE


! Allocate arrays for reading old restart file:
    allocate ( ak_r(npz_rst+1) )
    allocate ( bk_r(npz_rst+1) )

    allocate ( u_r(isc:iec,  jsc:jec+1,npz_rst) )
    allocate ( v_r(isc:iec+1,jsc:jec  ,npz_rst) )

    allocate (   pt_r(isc:iec, jsc:jec,  npz_rst) )
    allocate ( delp_r(isc:iec, jsc:jec,  npz_rst) )
    allocate (    q_r(isc:iec, jsc:jec,  npz_rst, ntracers) )

    if ( (.not.Atm(1)%flagstruct%hydrostatic) .and. (.not.Atm(1)%flagstruct%make_nh) ) then
           allocate (    w_r(isc:iec, jsc:jec,  npz_rst) )
           allocate ( delz_r(isc:iec, jsc:jec,  npz_rst) )
           if ( Atm(1)%flagstruct%hybrid_z )   &
           allocate ( ze0_r(isc:iec, jsc:jec,  npz_rst+1) )
    endif

    call set_domain(Atm(1)%domain)

    fname_nd = 'fv_core.res.nc'
    id_restart = register_restart_field(Fv_restart_r, fname_nd, 'ak', ak_r(:), no_domain=.true.)
    id_restart = register_restart_field(Fv_restart_r, fname_nd, 'bk', bk_r(:), no_domain=.true.)
    call restore_state(Fv_restart_r)

!    do n = 1, ntileMe
    n = 1
       fname = 'fv_core'//trim(gn)//'.res.nc'
       if (.not.Atm(n)%neststruct%nested) then
          call get_instance_filename(fname, fname_nd)
          call get_mosaic_tile_file(fname_nd, fname, .FALSE., fv_domain, n)
       end if
       id_restart =  register_restart_field(Fv_tile_restart_r, fname_nd, 'u', u_r, &
                     domain=fv_domain, position=NORTH,tile_count=n)
       id_restart =  register_restart_field(Fv_tile_restart_r, fname_nd, 'v', v_r, &
                     domain=fv_domain, position=EAST,tile_count=n)
       if (.not.Atm(n)%flagstruct%hydrostatic) then
          id_restart =  register_restart_field(Fv_tile_restart_r, fname_nd, 'W', w_r, &
                        domain=fv_domain, mandatory=.false., tile_count=n)
          id_restart =  register_restart_field(Fv_tile_restart_r, fname_nd, 'DZ', delz_r, &
                        domain=fv_domain, mandatory=.false., tile_count=n)
          if ( Atm(n)%flagstruct%hybrid_z ) then
             id_restart =  register_restart_field(Fv_tile_restart_r, fname_nd, 'ZE0', ze0_r, &
                           domain=fv_domain, mandatory=.false., tile_count=n)
          endif
       endif
       id_restart =  register_restart_field(Fv_tile_restart_r, fname_nd, 'T', pt_r, &
                     domain=fv_domain, tile_count=n)
       id_restart =  register_restart_field(Fv_tile_restart_r, fname_nd, 'delp', delp_r, &
                     domain=fv_domain, tile_count=n)
       id_restart =  register_restart_field(Fv_tile_restart_r, fname_nd, 'phis', Atm(n)%phis, &
                     domain=fv_domain, tile_count=n)
       call restore_state(FV_tile_restart_r)

       fname = 'fv_srf_wnd'//trim(gn)//'.res.nc'
       if (.not.Atm(n)%neststruct%nested) then
         call get_instance_filename(fname, fname_nd)
         call get_mosaic_tile_file(fname_nd, fname, .FALSE., fv_domain, n)
       end if
       if (file_exist(fname)) then
         Atm(n)%flagstruct%srf_init = .true.
       else
         call mpp_error(NOTE,'==> Warning from fv_read_restart: Expected file '//trim(fname)//' does not exist')
         Atm(n)%flagstruct%srf_init = .false.
       endif

       if ( Atm(n)%flagstruct%fv_land ) then
! Optional terrain deviation (sgh)
          call restore_state(Atm(n)%Mg_restart)
          call restore_state(Atm(n)%Lnd_restart)
       endif

       fname = 'fv_tracer'//trim(gn)//'.res.nc'
       if (.not.Atm(n)%neststruct%nested) then
          call get_instance_filename(fname, fname_nd)
          call get_mosaic_tile_file(fname_nd, fname, .FALSE., fv_domain, n)
       end if
       do nt = 1, ntracers
          call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
          call set_tracer_profile (MODEL_ATMOS, nt, q_r(isc:iec,jsc:jec,:,nt)  )
          id_restart = register_restart_field(Tra_restart_r, fname_nd, tracer_name, q_r(:,:,:,nt), &
                       domain=fv_domain, mandatory=.false., tile_count=n)
       enddo
       call restore_state(Tra_restart_r)

       call rst_remap(npz_rst, npz, isc, iec, jsc, jec, isd, ied, jsd, jed, ntracers,              &
                      delp_r,      u_r,      v_r,      w_r,      delz_r,      pt_r,      q_r,      &
                      Atm(n)%delp, Atm(n)%u, Atm(n)%v, Atm(n)%w, Atm(n)%delz, Atm(n)%pt, Atm(n)%q, &
                      ak_r,  bk_r, Atm(n)%ptop, Atm(n)%ak, Atm(n)%bk,                              &
                      Atm(n)%flagstruct%hydrostatic, Atm(n)%flagstruct%make_nh, Atm(n)%domain,     &
                      Atm(n)%gridstruct%square_domain)
    !end do

    call nullify_domain()

    deallocate( ak_r )
    deallocate( bk_r )
    deallocate( u_r )
    deallocate( v_r )
    deallocate( pt_r )
    deallocate( delp_r )
    deallocate( q_r )

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
    character(len=64) :: fname_nd
    integer           :: id_restart

! use_ncep_sst may not be initialized at this point?
    call mpp_error(NOTE, 'READING FROM SST_restart DISABLED')
!!$    if ( use_ncep_sst .or. Atm(1)%nudge .or. Atm(1)%ncep_ic ) then
!!$       fname_nd = 'sst_ncep.res.nc'
!!$       id_restart = register_restart_field(Atm(1)%SST_restart, fname_nd, 'sst_ncep', sst_ncep)
!!$       id_restart = register_restart_field(Atm(1)%SST_restart, fname_nd, 'sst_anom', sst_anom)
!!$    endif

  end subroutine  fv_io_register_nudge_restart
  ! </SUBROUTINE> NAME="fv_io_register_nudge_restart"


  !#####################################################################
  ! <SUBROUTINE NAME="fv_io_register_restart">
  !
  ! <DESCRIPTION>
  !   register restart field to be written out to restart file. 
  ! </DESCRIPTION>
  subroutine  fv_io_register_restart(fv_domain,Atm)
    type(domain2d),      intent(inout) :: fv_domain
    type(fv_atmos_type), intent(inout) :: Atm(:)

    character(len=64) :: fname_nd, tracer_name
    character(len=3)  :: gn
    character(len=6)  :: stile_name
    integer           :: id_restart
    integer           :: n, nt, ntracers, ntileMe, ntiles

    ntileMe = size(Atm(:)) 
    ntracers = size(Atm(1)%q,4) 

    if (Atm(1)%grid_number > 1) then
       write(gn,'(A2, I1)') "_g", Atm(1)%grid_number
    else
       gn = ''
    end if

    call set_domain(Atm(1)%domain)

    fname_nd = 'fv_core.res.nc'
    id_restart = register_restart_field(Atm(1)%Fv_restart, fname_nd, 'ak', Atm(1)%ak(:), no_domain=.true.)
    id_restart = register_restart_field(Atm(1)%Fv_restart, fname_nd, 'bk', Atm(1)%bk(:), no_domain=.true.) 

! use_ncep_sst may not be initialized at this point?
#ifndef DYCORE_SOLO
    call mpp_error(NOTE, 'READING FROM SST_RESTART DISABLED')
!!$   if ( use_ncep_sst .or. Atm(1)%flagstruct%nudge .or. Atm(1)%flagstruct%ncep_ic ) then
!!$       fname_nd = 'sst_ncep'//trim(gn)//'.res.nc'
!!$       id_restart = register_restart_field(Atm(1)%SST_restart, fname_nd, 'sst_ncep', sst_ncep)
!!$       id_restart = register_restart_field(Atm(1)%SST_restart, fname_nd, 'sst_anom', sst_anom)
!!$   endif
#endif

! fix for single tile runs where you need fv_core.res.nc and fv_core.res.tile1.nc
    ntiles = mpp_get_ntile_count(fv_domain)
    if(ntiles == 1 .and. .not. Atm(1)%neststruct%nested) then
       stile_name = '.tile1'
    else
       stile_name = ''
    endif

    fname_nd = 'fv_core'//trim(gn)//'.res'//trim(stile_name)//'.nc'
    do n = 1, ntileMe
!    n = 1
       id_restart =  register_restart_field(Atm(n)%Fv_tile_restart, fname_nd, 'u', Atm(n)%u, &
                     domain=fv_domain, position=NORTH,tile_count=n)
       id_restart =  register_restart_field(Atm(n)%Fv_tile_restart, fname_nd, 'v', Atm(n)%v, &
                     domain=fv_domain, position=EAST,tile_count=n)
       if (.not.Atm(n)%flagstruct%hydrostatic) then
          id_restart =  register_restart_field(Atm(n)%Fv_tile_restart, fname_nd, 'W', Atm(n)%w, &
                        domain=fv_domain, mandatory=.false., tile_count=n)
          id_restart =  register_restart_field(Atm(n)%Fv_tile_restart, fname_nd, 'DZ', Atm(n)%delz, &
                        domain=fv_domain, mandatory=.false., tile_count=n)
          if ( Atm(n)%flagstruct%hybrid_z ) then
             id_restart =  register_restart_field(Atm(n)%Fv_tile_restart, fname_nd, 'ZE0', Atm(n)%ze0, &
                           domain=fv_domain, mandatory=.false., tile_count=n)
          endif
       endif
       id_restart =  register_restart_field(Atm(n)%Fv_tile_restart, fname_nd, 'T', Atm(n)%pt, &
                     domain=fv_domain, tile_count=n)
       id_restart =  register_restart_field(Atm(n)%Fv_tile_restart, fname_nd, 'delp', Atm(n)%delp, &
                     domain=fv_domain, tile_count=n)
       id_restart =  register_restart_field(Atm(n)%Fv_tile_restart, fname_nd, 'phis', Atm(n)%phis, &
                     domain=fv_domain, tile_count=n)

       fname_nd = 'fv_srf_wnd'//trim(gn)//'.res.nc'
       !if(ntiles == 1) fname_nd = 'fv_srf_wnd'//trim(gn)//'.res.tile1.nc'
       id_restart =  register_restart_field(Atm(n)%Rsf_restart, fname_nd, 'u_srf', Atm(n)%u_srf, &
                     domain=fv_domain, tile_count=n)
       id_restart =  register_restart_field(Atm(n)%Rsf_restart, fname_nd, 'v_srf', Atm(n)%v_srf, &
                     domain=fv_domain, tile_count=n)
#ifdef SIM_PHYS
       id_restart =  register_restart_field(Rsf_restart(n), fname_nd, 'ts', Atm(n)%ts, &
                     domain=fv_domain, tile_count=n)
#endif

       if ( Atm(n)%flagstruct%fv_land ) then
          !-------------------------------------------------------------------------------------------------
          ! Optional terrain deviation (sgh) and land fraction (oro)
          fname_nd = 'mg_drag'//trim(gn)//'.res.nc'
          id_restart =  register_restart_field(Atm(n)%Mg_restart, fname_nd, 'ghprime', Atm(n)%sgh, &
                        domain=fv_domain, tile_count=n)  

          fname_nd = 'fv_land'//trim(gn)//'.res.nc'
          id_restart = register_restart_field(Atm(n)%Lnd_restart, fname_nd, 'oro', Atm(n)%oro, &
                        domain=fv_domain, tile_count=n)
       endif
       fname_nd = 'fv_tracer'//trim(gn)//'.res.nc'
       !if(ntiles == 1) fname_nd = 'fv_tracer'//trim(gn)//'.res.tile1.nc'
       do nt = 1, ntracers
          call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
          ! set all tracers to an initial profile value
          call set_tracer_profile (MODEL_ATMOS, nt, Atm(n)%q(:,:,:,nt)  )
          id_restart = register_restart_field(Atm(n)%Tra_restart, fname_nd, tracer_name, Atm(n)%q(:,:,:,nt), &
                       domain=fv_domain, mandatory=.false., tile_count=n)
       enddo

    enddo

  end subroutine  fv_io_register_restart
  ! </SUBROUTINE> NAME="fv_io_register_restart"



  !#####################################################################
  ! <SUBROUTINE NAME="fv_io_write_restart">
  !
  ! <DESCRIPTION>
  ! Write the fv core restart quantities 
  ! </DESCRIPTION>
  subroutine  fv_io_write_restart(Atm, timestamp)

    type(fv_atmos_type),        intent(inout) :: Atm(:)
    character(len=*), optional, intent(in) :: timestamp
    integer                                :: n, ntileMe

    ntileMe = size(Atm(:))  ! This will need mods for more than 1 tile per pe
    
    call set_domain(Atm(1)%domain)

    if (.not. Atm(1)%neststruct%nested) then
       call save_restart(Atm(1)%Fv_restart, timestamp)
    endif

    if ( use_ncep_sst .or. Atm(1)%flagstruct%nudge .or. Atm(1)%flagstruct%ncep_ic ) then
       call mpp_error(NOTE, 'READING FROM SST_RESTART DISABLED')
       !call save_restart(Atm(1)%SST_restart, timestamp)
    endif

    call nullify_domain
 
    do n = 1, ntileMe
!    n = 1
       call set_domain(Atm(n)%domain)
       call save_restart(Atm(n)%Fv_tile_restart, timestamp)
       call save_restart(Atm(n)%Rsf_restart, timestamp)

       if ( Atm(n)%flagstruct%fv_land ) then
          call save_restart(Atm(n)%Mg_restart, timestamp)
          call save_restart(Atm(n)%Lnd_restart, timestamp)
       endif

       call save_restart(Atm(n)%Tra_restart, timestamp)

       call nullify_domain

    end do

  end subroutine  fv_io_write_restart

  subroutine register_bcs_2d(Atm, BCfile_ne, BCfile_sw, fname_ne, fname_sw, &
                             var_name, var, var_bc, istag, jstag)
    type(fv_atmos_type),      intent(in)    :: Atm
    type(restart_file_type),  intent(inout) :: BCfile_ne, BCfile_sw
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
    if (present(var_bc)) id_restart = register_restart_field(BCfile_sw, trim(fname_sw), &
                                        trim(var_name)//'_west_t1', &
                                        var_bc%west_t1, & 
                                        indices, global_size, y2_pelist, &
                                        is_root_pe, jshift=y_halo)
!register west prognostic halo data
    if (present(var)) id_restart = register_restart_field(BCfile_sw, trim(fname_sw), &
                                        trim(var_name)//'_west', &
                                        var, indices, global_size, &
                                        y2_pelist, is_root_pe, jshift=y_halo)

!define east root_pe
    is_root_pe = .FALSE.
    if (ie.eq.npx-1 .and. je.eq.npy-1) is_root_pe = .TRUE.
!register east halo data in t1
    if (present(var_bc)) id_restart = register_restart_field(BCfile_ne, trim(fname_ne), &
                                        trim(var_name)//'_east_t1', &
                                        var_bc%east_t1, & 
                                        indices, global_size, y1_pelist, &
                                        is_root_pe, jshift=y_halo)

!reset indices for prognostic variables in the east halo
    indices(1) = ied-x_halo+1+i_stag
    indices(2) = ied+i_stag
!register east prognostic halo data
    if (present(var)) id_restart = register_restart_field(BCfile_ne, trim(fname_ne), &
                                        trim(var_name)//'_east', &
                                        var, indices, global_size, &
                                        y1_pelist, is_root_pe, jshift=y_halo, &
                                        x_halo=(size(var,1)-x_halo), ishift=-(ie+i_stag))

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
    if (present(var_bc)) id_restart = register_restart_field(BCfile_sw, trim(fname_sw), &
                                        trim(var_name)//'_south_t1', &
                                        var_bc%south_t1, & 
                                        indices, global_size, x2_pelist, &
                                        is_root_pe, x_halo=x_halo_ns)
!register south prognostic halo data
    if (present(var)) id_restart = register_restart_field(BCfile_sw, trim(fname_sw), &
                                        trim(var_name)//'_south', &
                                        var, indices, global_size, &
                                        x2_pelist, is_root_pe, x_halo=x_halo_ns)

!define north root_pe
    is_root_pe = .FALSE.
    if (ie.eq.npx-1 .and. je.eq.npy-1) is_root_pe = .TRUE.
!register north halo data in t1
    if (present(var_bc)) id_restart = register_restart_field(BCfile_ne, trim(fname_ne), &
                                        trim(var_name)//'_north_t1', &
                                        var_bc%north_t1, & 
                                        indices, global_size, x1_pelist, &
                                        is_root_pe, x_halo=x_halo_ns)

!reset indices for prognostic variables in the north halo
    indices(3) = jed-y_halo+1+j_stag
    indices(4) = jed+j_stag
!register north prognostic halo data
    if (present(var)) id_restart = register_restart_field(BCfile_ne, trim(fname_ne), &
                                        trim(var_name)//'_north', &
                                        var, indices, global_size, &
                                        x1_pelist, is_root_pe, x_halo=x_halo_ns, &
                                        y_halo=(size(var,2)-y_halo), jshift=-(je+j_stag))

  end subroutine register_bcs_2d


  subroutine register_bcs_3d(Atm, BCfile_ne, BCfile_sw, fname_ne, fname_sw, &
                             var_name, var, var_bc, istag, jstag)
    type(fv_atmos_type),      intent(in)    :: Atm
    type(restart_file_type),  intent(inout) :: BCfile_ne, BCfile_sw
    character(len=120),       intent(in)    :: fname_ne, fname_sw
    character(len=*),         intent(in)    :: var_name
    real, dimension(:,:,:),   intent(in), optional :: var
    type(fv_nest_BC_type_3D), intent(in), optional :: var_bc
    integer,                  intent(in), optional :: istag, jstag

    integer :: npx, npy, i_stag, j_stag
    integer :: is, ie, js, je, isd, ied, jsd, jed, n
    integer :: x_halo, y_halo, x_halo_ns, id_restart
    integer :: layout(2), global_size(3), indices(4)
    integer, allocatable, dimension(:) :: x1_pelist, y1_pelist
    integer, allocatable, dimension(:) :: x2_pelist, y2_pelist
    logical :: is_root_pe

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
    if (present(var_bc)) id_restart = register_restart_field(BCfile_sw, trim(fname_sw), &
                                        trim(var_name)//'_west_t1', &
                                        var_bc%west_t1, & 
                                        indices, global_size, y2_pelist, &
                                        is_root_pe, jshift=y_halo)
!register west prognostic halo data
    if (present(var)) id_restart = register_restart_field(BCfile_sw, trim(fname_sw), &
                                        trim(var_name)//'_west', &
                                        var, indices, global_size, &
                                        y2_pelist, is_root_pe, jshift=y_halo)

!define east root_pe
    is_root_pe = .FALSE.
    if (ie.eq.npx-1 .and. je.eq.npy-1) is_root_pe = .TRUE.
!register east halo data in t1
    if (present(var_bc)) id_restart = register_restart_field(BCfile_ne, trim(fname_ne), &
                                        trim(var_name)//'_east_t1', &
                                        var_bc%east_t1, & 
                                        indices, global_size, y1_pelist, &
                                        is_root_pe, jshift=y_halo)

!reset indices for prognostic variables in the east halo
    indices(1) = ied-x_halo+1+i_stag
    indices(2) = ied+i_stag
!register east prognostic halo data
    if (present(var)) id_restart = register_restart_field(BCfile_ne, trim(fname_ne), &
                                        trim(var_name)//'_east', &
                                        var, indices, global_size, &
                                        y1_pelist, is_root_pe, jshift=y_halo, &
                                        x_halo=(size(var,1)-x_halo), ishift=-(ie+i_stag))

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
    if (present(var_bc)) id_restart = register_restart_field(BCfile_sw, trim(fname_sw), &
                                        trim(var_name)//'_south_t1', &
                                        var_bc%south_t1, & 
                                        indices, global_size, x2_pelist, &
                                        is_root_pe, x_halo=x_halo_ns)
!register south prognostic halo data
    if (present(var)) id_restart = register_restart_field(BCfile_sw, trim(fname_sw), &
                                        trim(var_name)//'_south', &
                                        var, indices, global_size, &
                                        x2_pelist, is_root_pe, x_halo=x_halo_ns)

!define north root_pe
    is_root_pe = .FALSE.
    if (ie.eq.npx-1 .and. je.eq.npy-1) is_root_pe = .TRUE.
!register north halo data in t1
    if (present(var_bc)) id_restart = register_restart_field(BCfile_ne, trim(fname_ne), &
                                        trim(var_name)//'_north_t1', &
                                        var_bc%north_t1, & 
                                        indices, global_size, x1_pelist, &
                                        is_root_pe, x_halo=x_halo_ns)

!reset indices for prognostic variables in the north halo
    indices(3) = jed-y_halo+1+j_stag
    indices(4) = jed+j_stag
!register north prognostic halo data
    if (present(var)) id_restart = register_restart_field(BCfile_ne, trim(fname_ne), &
                                        trim(var_name)//'_north', &
                                        var, indices, global_size, &
                                        x1_pelist, is_root_pe, x_halo=x_halo_ns, &
                                        y_halo=(size(var,2)-y_halo), jshift=-(je+j_stag))

  end subroutine register_bcs_3d


  ! </SUBROUTINE> NAME="fv_io_regsiter_restart_BCs"
  !#####################################################################

  subroutine fv_io_register_restart_BCs(Atm)
!!! CLEANUP: The make_nh option does not yet work with nesting since we don't have a routine yet to fill the halo with just w = 0 and setting up delz.
    type(fv_atmos_type),        intent(inout) :: Atm

    integer :: n
    character(len=120) :: tname, fname_ne, fname_sw
    character(len=2) :: gn

    if (Atm%grid_number > 1) then
      write(gn,'(i2.2)') Atm%grid_number
    else
      gn = ''
    end if
    fname_ne = 'fv_BCnest'//gn//'_ne.res.nc'
    fname_sw = 'fv_BCnest'//gn//'_sw.res.nc'

    call set_domain(Atm%domain)

    call register_bcs_2d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                         fname_ne, fname_sw, 'phis', var=Atm%phis)
    call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                         fname_ne, fname_sw, 'delp', Atm%delp, Atm%neststruct%delp_BC)
    do n=1,Atm%ncnst
       call get_tracer_names(MODEL_ATMOS, n, tname)
       call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                            fname_ne, fname_sw, trim(tname), Atm%q(:,:,:,n), Atm%neststruct%q_BC(n))
    enddo
#ifndef SW_DYNAMICS
    call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                         fname_ne, fname_sw, 'pt', Atm%pt, Atm%neststruct%pt_BC)
    if ((.not.Atm%flagstruct%hydrostatic) .and. (.not.Atm%flagstruct%make_nh)) then
      call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                           fname_ne, fname_sw, 'w', Atm%w, Atm%neststruct%w_BC)
      call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                           fname_ne, fname_sw, 'delz', Atm%delz, Atm%neststruct%delz_BC)
    endif
#endif
    call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                         fname_ne, fname_sw, 'u', Atm%u, Atm%neststruct%u_BC, jstag=1)
    call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                         fname_ne, fname_sw, 'v', Atm%v, Atm%neststruct%v_BC, istag=1)
    call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                         fname_ne, fname_sw, 'uc', var_bc=Atm%neststruct%uc_BC, istag=1)
    call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                         fname_ne, fname_sw, 'vc', var_bc=Atm%neststruct%vc_BC, jstag=1)

    call nullify_domain

    return
  end subroutine fv_io_register_restart_BCs


  subroutine fv_io_write_BCs(Atm, timestamp)
    type(fv_atmos_type), intent(inout) :: Atm
    character(len=*),    intent(in), optional :: timestamp

    call save_restart_border(Atm%neststruct%BCfile_ne, timestamp)
    call save_restart_border(Atm%neststruct%BCfile_sw, timestamp)

    return
  end subroutine fv_io_write_BCs


  subroutine fv_io_read_BCs(Atm)
    type(fv_atmos_type), intent(inout) :: Atm

    call restore_state_border(Atm%neststruct%BCfile_ne)
    call restore_state_border(Atm%neststruct%BCfile_sw)

    return
  end subroutine fv_io_read_BCs
    
end module fv_io_mod
