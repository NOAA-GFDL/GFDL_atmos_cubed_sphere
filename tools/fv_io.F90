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

  use fms_mod,                 only: file_exist, read_data, write_data, field_exist    
  use fms_io_mod,              only: fms_io_exit, get_tile_string, &
                                     restart_file_type, register_restart_field, &
                                     save_restart, restore_state, &
                                     set_domain, nullify_domain, &
                                     get_mosaic_tile_file, get_instance_filename
  use mpp_mod,                 only: mpp_error, FATAL, NOTE
  use mpp_domains_mod,         only: domain2d, EAST, NORTH, mpp_get_tile_id, &
                                     mpp_get_compute_domain, mpp_get_data_domain, &
                                     mpp_get_ntile_count
  use tracer_manager_mod,      only: tr_get_tracer_names=>get_tracer_names, &
                                     get_tracer_names, get_number_tracers, &
                                     set_tracer_profile, &
                                     get_tracer_index
  use field_manager_mod,       only: MODEL_ATMOS  
  use external_sst_mod,        only: sst_ncep, sst_anom, use_ncep_sst
  use fv_arrays_mod,           only: fv_atmos_type
  use fv_current_grid_mod,     only: Fv_restart, SST_restart, Fv_tile_restart, &
                                     Rsf_restart, Mg_restart, Lnd_restart, Tra_restart
  use fv_eta_mod,              only: set_eta

   use fv_mp_mod, only: concurrent
   use fms_io_mod, only: set_domain

  implicit none
  private

  public :: fv_io_init, fv_io_exit, fv_io_read_restart, remap_restart, fv_io_write_restart
  public :: fv_io_read_tracers, fv_io_register_restart, fv_io_register_nudge_restart
  public :: fv_io_write_BCs, fv_io_read_BCs

  logical                       :: module_is_initialized = .FALSE.


!---- version number -----
  character(len=128) :: version = '$Id: fv_io.F90,v 17.0.2.6.2.1.4.4.2.14.4.1 2013/02/27 17:07:03 Seth.Underwood Exp $'
  character(len=128) :: tagname = '$Name: siena_201303 $'

  integer ::grid_xtdimid, grid_ytdimid, haloid, pfullid !For writing BCs
  integer ::grid_xtstagdimid, grid_ytstagdimid, oneid

contains 

  !#####################################################################
  ! <SUBROUTINE NAME="fv_io_initFV_IO_REGISTER_NUDGE">
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
!    integer              :: ntileMe
    integer              :: ks
    real                 :: ptop

    character(len=128)           :: tracer_longname, tracer_units

    if (Atm(1)%grid_number > 1) then
       write(gn,'(A2, I1)') "_g", Atm(1)%grid_number
    else
       gn = ''
    end if

!    ntileMe = size(Atm(:))  ! This will have to be modified for mult tiles per PE
 
!   call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)
    ntracers = size(Atm(1)%q,4)  ! Temporary until we get tracer manager integrated

#ifndef DYCORE_SOLO
    call mpp_error(NOTE, 'READING FROM INPUT/sst_ncep DISABLED')
!    fname_nd = 'INPUT/sst_ncep'//trim(gn)//'.res.nc'
!    if(file_exist(fname_nd))then
! sst_ncep may be used in free-running forecast mode
!       call read_data(fname_nd, 'sst_ncep', sst_ncep)
!       if (field_exist(fname_nd,'sst_anom')) then
!           call read_data(fname_nd, 'sst_anom', sst_anom)
!       else
!           sst_anom(:,:) = 1.E35   ! make it big enough to cause blowup if used
!       endif
!    elseif ( use_ncep_sst .and. (.not. Atm(1)%nudge) ) then
!       call mpp_error(FATAL,'==> Error:'//trim(fname_nd)//' does not exist')
!    endif
#endif 
    call nullify_domain !Needed for the following read_data calls; else they
                        !will try to read from .res.tileN.nc, and crash when
                        !it can't find ak and bk in those files
  ! write_data does not (yet?) support vector data and tiles
    call get_instance_filename('INPUT/fv_core.res.nc', fname_nd)
    call read_data(fname_nd, 'ak', Atm(1)%ak(:))
    call read_data(fname_nd, 'bk', Atm(1)%bk(:))
    if ( Atm(1)%reset_eta ) call set_eta(Atm(1)%npz, ks, ptop, Atm(1)%ak, Atm(1)%bk)
    call set_domain(fv_domain)
   
!    do n = 1, ntileMe
    n = 1
       isc = Atm(n)%isc; iec = Atm(n)%iec; jsc = Atm(n)%jsc; jec = Atm(n)%jec

       fname = 'INPUT/fv_core'//trim(gn)//'.res.nc'
       if (.not.Atm(n)%nested) then
          call get_instance_filename(fname, fname_nd)
          call get_mosaic_tile_file(fname_nd, fname, .FALSE., fv_domain, n)
       end if
       if(file_exist(fname))then
         call read_data(fname, 'u', Atm(n)%u(isc:iec,jsc:jec+1,:), domain=fv_domain, position=NORTH,tile_count=n)
         call read_data(fname, 'v', Atm(n)%v(isc:iec+1,jsc:jec,:), domain=fv_domain, position=EAST,tile_count=n)

         if ( (.not.Atm(n)%hydrostatic) .and. (.not.Atm(n)%make_nh) ) then
              call read_data(fname, 'W',     Atm(n)%w(isc:iec,jsc:jec,:), domain=fv_domain, tile_count=n)
              call read_data(fname, 'DZ', Atm(n)%delz(isc:iec,jsc:jec,:), domain=fv_domain, tile_count=n)
              if ( Atm(n)%hybrid_z .and. (.not. Atm(n)%make_hybrid_z) )   &
              call read_data(fname, 'ZE0', Atm(n)%ze0(isc:iec,jsc:jec,:), domain=fv_domain, tile_count=n)
         endif

         call read_data(fname, 'T', Atm(n)%pt(isc:iec,jsc:jec,:), domain=fv_domain, tile_count=n)
         call read_data(fname, 'delp', Atm(n)%delp(isc:iec,jsc:jec,:), domain=fv_domain, tile_count=n)
         call read_data(fname, 'phis', Atm(n)%phis(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
       else
         call mpp_error(FATAL,'==> Error from fv_read_restart: Expected file '//trim(fname)//' does not exist')
       endif

       fname = 'INPUT/fv_srf_wnd'//trim(gn)//'.res.nc'
       if (.not.Atm(n)%nested) then
          call get_instance_filename(fname, fname_nd)
          call get_mosaic_tile_file(fname_nd, fname, .FALSE., fv_domain, n)
       end if
       if(file_exist(fname))then
         call read_data(fname, 'u_srf', Atm(n)%u_srf(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
         call read_data(fname, 'v_srf', Atm(n)%v_srf(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
         if (field_exist(fname,'ts'))   &
               call read_data(fname, 'ts', Atm(n)%ts(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
         Atm(n)%srf_init = .true.
       else
         call mpp_error(NOTE,'==> Warning from fv_read_restart: Expected file '//trim(fname)//' does not exist')
         Atm(n)%srf_init = .false.
       endif

  if ( Atm(n)%fv_land ) then
!----------------------------------------------------------------------------------------------------------------
! Optional terrain deviation (sgh) and land fraction (oro)
       fname = 'INPUT/mg_drag'//trim(gn)//'.res.nc'
       if (.not.Atm(n)%nested) then
          call get_instance_filename(fname, fname_nd)
          call get_mosaic_tile_file(fname_nd, fname, .FALSE., fv_domain, n)
       end if
       if(file_exist(fname))then
         call read_data(fname, 'ghprime', Atm(n)%sgh(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
       else
         call mpp_error(NOTE,'==> Warning from fv_read_restart: Expected file '//trim(fname)//' does not exist')
       endif
! Land-water mask:
       fname = 'INPUT/fv_land'//trim(gn)//'.res.nc'
       if (.not.Atm(n)%nested) then
          call get_instance_filename(fname, fname_nd)
          call get_mosaic_tile_file(fname_nd, fname, .FALSE., fv_domain, n)
       end if
       if(file_exist(fname))then
         call read_data(fname, 'oro', Atm(n)%oro(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
       else
         call mpp_error(NOTE,'==> Warning from fv_read_restart: Expected file '//trim(fname)//' does not exist')
       endif
!----------------------------------------------------------------------------------------------------------------
  endif
       fname = 'INPUT/fv_tracer'//trim(gn)//'.res.nc'
       if (.not.Atm(n)%nested) then
          call get_instance_filename(fname, fname_nd)
          call get_mosaic_tile_file(fname_nd, fname, .FALSE., fv_domain, n)
       end if

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

!    end do
 
    call nullify_domain()

  end subroutine  fv_io_read_restart
  ! </SUBROUTINE> NAME="fv_io_read_restart"
  !#####################################################################


  subroutine fv_io_read_tracers(fv_domain,Atm)
    type(domain2d),      intent(inout) :: fv_domain
    type(fv_atmos_type), intent(inout) :: Atm(:)

    character(len=64)    :: fname, fname_nd, tracer_name
    character(len=3)  :: gn
    integer              :: isc, iec, jsc, jec, n, nt, nk, ntracers
!    integer              :: ntileMe

    character(len=128)           :: tracer_longname, tracer_units

    if (Atm(1)%grid_number > 1) then
       write(gn,'(A2, I1)') "_g", Atm(1)%grid_number
    else
       gn = ''
    end if

!    ntileMe = size(Atm(:))  ! This will have to be modified for mult tiles per PE

!   call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)
    ntracers = size(Atm(1)%q,4)  ! Temporary until we get tracer manager integrated
    call set_domain(fv_domain)

! skip the first tracer, which is sphum
!    do n = 1, ntileMe
    n = 1
       isc = Atm(n)%isc; iec = Atm(n)%iec; jsc = Atm(n)%jsc; jec = Atm(n)%jec
       fname = 'INPUT/fv_tracer'//trim(gn)//'.res.nc'
       if (.not.Atm(n)%nested) then
          call get_instance_filename(fname, fname_nd)
          call get_mosaic_tile_file(fname_nd, fname, .FALSE., fv_domain, n)
       end if
         DO nt = 2, ntracers
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
!    end do

    call nullify_domain()

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
    npz_rst = Atm(1)%npz_rst   ! restart z dimension
    isc = Atm(1)%isc; iec = Atm(1)%iec; jsc = Atm(1)%jsc; jec = Atm(1)%jec
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

    if ( (.not.Atm(1)%hydrostatic) .and. (.not.Atm(1)%make_nh) ) then
           allocate (    w_r(isc:iec, jsc:jec,  npz_rst) )
           allocate ( delz_r(isc:iec, jsc:jec,  npz_rst) )
           if ( Atm(1)%hybrid_z )   &
           allocate ( ze0_r(isc:iec, jsc:jec,  npz_rst+1) )
    endif

    call get_instance_filename('INPUT/fv_core'//trim(gn)//'.res.nc', fname_nd)
  ! write_data does not (yet?) support vector data and tiles
    call read_data(fname_nd, 'ak', ak_r(1:npz_rst+1))
    call read_data(fname_nd, 'bk', bk_r(1:npz_rst+1))

    call set_domain(fv_domain)

!    do n = 1, ntileMe
    n = 1
      fname = 'INPUT/fv_core'//trim(gn)//'.res.nc'
       if (.not.Atm(n)%nested) then
          call get_instance_filename(fname, fname_nd)
          call get_mosaic_tile_file(fname_nd, fname, .FALSE., fv_domain, n)
       end if
       if(file_exist(fname))then
         call read_data(fname, 'u', u_r(isc:iec,jsc:jec+1,:), domain=fv_domain, position=NORTH,tile_count=n)
         call read_data(fname, 'v', v_r(isc:iec+1,jsc:jec,:), domain=fv_domain, position=EAST,tile_count=n)

         if ( (.not.Atm(n)%hydrostatic) .and. (.not.Atm(n)%make_nh) ) then
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


       fname = 'INPUT/fv_srf_wnd'//trim(gn)//'.res.nc'
       if (.not.Atm(n)%nested) then
          call get_instance_filename(fname, fname_nd)
          call get_mosaic_tile_file(fname_nd, fname, .FALSE., fv_domain, n)
       end if
       if(file_exist(fname))then
         call read_data(fname, 'u_srf', Atm(n)%u_srf(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
         call read_data(fname, 'v_srf', Atm(n)%v_srf(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
         if (field_exist(fname,'ts'))   &
               call read_data(fname, 'ts', Atm(n)%ts(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
         Atm(n)%srf_init = .true.
       else
         call mpp_error(NOTE,'==> Warning from fv_read_restart: Expected file '//trim(fname)//' does not exist')
         Atm(n)%srf_init = .false.
       endif

     if ( Atm(n)%fv_land ) then
! Optional terrain deviation (sgh)
       fname = 'INPUT/mg_drag'//trim(gn)//'.res.nc'
       if (.not.Atm(n)%nested) then
          call get_instance_filename(fname, fname_nd)
          call get_mosaic_tile_file(fname_nd, fname, .FALSE., fv_domain, n)
       end if
       if(file_exist(fname))then
         call read_data(fname, 'ghprime', Atm(n)%sgh(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
       else
         call mpp_error(NOTE,'==> Warning from fv_read_restart: Expected file '//trim(fname)//' does not exist')
       endif
! Land-water mask
       fname = 'INPUT/fv_land'//trim(gn)//'.res.nc'
       if (.not.Atm(n)%nested) then
          call get_instance_filename(fname, fname_nd)
          call get_mosaic_tile_file(fname_nd, fname, .FALSE., fv_domain, n)
       end if
       if(file_exist(fname))then
         call read_data(fname, 'oro', Atm(n)%oro(isc:iec,jsc:jec), domain=fv_domain, tile_count=n)
       else
         call mpp_error(NOTE,'==> Warning from fv_read_restart: Expected file '//trim(fname)//' does not exist')
       endif
     endif


       fname = 'INPUT/fv_tracer'//trim(gn)//'.res.nc'
       if (.not.Atm(n)%nested) then
          call get_instance_filename(fname, fname_nd)
          call get_mosaic_tile_file(fname_nd, fname, .FALSE., fv_domain, n)
       end if

       do nt = 1, ntracers
           call get_tracer_names(MODEL_ATMOS, nt, tracer_name)

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
    !end do

    call nullify_domain()

    deallocate( ak_r )
    deallocate( bk_r )
    deallocate( u_r )
    deallocate( v_r )
    deallocate( pt_r )
    deallocate( delp_r )
    deallocate( q_r )

    if ( (.not.Atm(1)%hydrostatic) .and. (.not.Atm(1)%make_nh) ) then
         deallocate ( w_r )
         deallocate ( delz_r )
         if ( Atm(1)%hybrid_z ) deallocate ( ze0_r )
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
    integer           :: id_restart
    integer           :: n, nt, ntracers, ntileMe, ntiles

    ntileMe = size(Atm(:)) 
    ntracers = size(Atm(1)%q,4) 

    if (Atm(1)%grid_number > 1) then
       write(gn,'(A2, I1)') "_g", Atm(1)%grid_number
    else
       gn = ''
    end if

    call nullify_domain !should this go here or before the previous 2 register_restart_field calls?

    fname_nd = 'fv_core.res.nc'

    id_restart = register_restart_field(Atm(1)%Fv_restart, fname_nd, 'ak', Atm(1)%ak(:))
    id_restart = register_restart_field(Atm(1)%Fv_restart, fname_nd, 'bk', Atm(1)%bk(:)) 

! use_ncep_sst may not be initialized at this point?
#ifndef DYCORE_SOLO
    call mpp_error(NOTE, 'READING FROM SST_RESTART DISABLED')
!!$   if ( use_ncep_sst .or. Atm(1)%nudge .or. Atm(1)%ncep_ic ) then
!!$       fname_nd = 'sst_ncep'//trim(gn)//'.res.nc'
!!$       id_restart = register_restart_field(Atm(1)%SST_restart, fname_nd, 'sst_ncep', sst_ncep)
!!$       id_restart = register_restart_field(Atm(1)%SST_restart, fname_nd, 'sst_anom', sst_anom)
!!$   endif
#endif

       call set_domain(Atm(1)%domain)

! fix for single tile runs where you need fv_core.res.nc and fv_core.res.tile1.nc
    fname_nd = 'fv_core'//trim(gn)//'.res.nc'
    ntiles = mpp_get_ntile_count(fv_domain)
    !if(ntiles == 1) fname_nd =  'fv_core'//trim(gn)//'.res.tile1.nc'

    do n = 1, ntileMe
!    n = 1
       id_restart =  register_restart_field(Atm(n)%Fv_tile_restart, fname_nd, 'u', Atm(n)%u, &
                     domain=fv_domain, position=NORTH,tile_count=n)
       id_restart =  register_restart_field(Atm(n)%Fv_tile_restart, fname_nd, 'v', Atm(n)%v, &
                     domain=fv_domain, position=EAST,tile_count=n)
       if (.not.Atm(n)%hydrostatic) then
          id_restart =  register_restart_field(Atm(n)%Fv_tile_restart, fname_nd, 'W', Atm(n)%w, &
                        domain=fv_domain, mandatory=.false., tile_count=n)
          id_restart =  register_restart_field(Atm(n)%Fv_tile_restart, fname_nd, 'DZ', Atm(n)%delz, &
                        domain=fv_domain, mandatory=.false., tile_count=n)
          if ( Atm(n)%hybrid_z ) then
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

       if ( Atm(n)%fv_land ) then
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
    use fv_mp_mod, only: masterproc, gid, ng !DEBUG

    type(fv_atmos_type),        intent(inout) :: Atm(:)
    character(len=*), optional, intent(in) :: timestamp
    integer                                :: n, ntileMe

    ntileMe = size(Atm(:))  ! This will need mods for more than 1 tile per pe
    
    if (Atm(1)%ntiles > 1) then
       call save_restart(Atm(1)%Fv_restart, timestamp)
    endif

    if ( use_ncep_sst .or. Atm(1)%nudge .or. Atm(1)%ncep_ic ) then
       call mpp_error(NOTE, 'READING FROM SST_RESTART DISABLED')
       !call save_restart(Atm(1)%SST_restart, timestamp)
    endif
 
    do n = 1, ntileMe
!    n = 1
       call save_restart(Atm(n)%Fv_tile_restart, timestamp)
       call save_restart(Atm(n)%Rsf_restart, timestamp)

       if ( Atm(n)%fv_land ) then
          call save_restart(Atm(n)%Mg_restart, timestamp)
          call save_restart(Atm(n)%Lnd_restart, timestamp)
       endif

       call save_restart(Atm(n)%Tra_restart, timestamp)

    end do

  end subroutine  fv_io_write_restart
  ! </SUBROUTINE> NAME="fv_io_write_restart"
  !#####################################################################

  subroutine fv_io_write_BCs(Atm)

    use sim_nc_mod, only: handle_err
    use fv_mp_mod, only: masterproc, gid, ng
    use mpp_domains_mod, only: mpp_get_global_domain, CORNER
    use tracer_manager_mod, only: get_tracer_names
    use field_manager_mod,       only: MODEL_ATMOS
#include <netcdf.inc>

    type(fv_atmos_type),        intent(inout) :: Atm
    character(len=120) :: fname, tname
    integer :: BCfileid, status
    character(len=3) :: gn
    integer :: npx, npy
    integer :: is, ie, js, je, isd, ied, jsd, jed, n


    !Easiest way to handle this: send all BC data to one PE.
    !Remember we need to save TWO times for concurrent simulations

    call mpp_get_global_domain(Atm%domain, xsize = npx, ysize = npy, position=CORNER )
    call mpp_get_data_domain( Atm%domain, isd, ied, jsd, jed )
    call mpp_get_compute_domain( Atm%domain, is, ie, js, je )

    if (gid == masterproc) then

       !Now save the data
    
       if (Atm%grid_number > 1) then
          write(gn,'(A2, I1)') "_g", Atm%grid_number
       else
          gn = ''
       end if
       fname = 'RESTART/fv_BCFile'//gn//'.nc'

       status = nf_create(trim(fname), NF_CLOBBER, BCfileid)
       if (status .ne. NF_NOERR) then
          print*, 'Could not create nested BC file ', trim(fname)
          call handle_err(status)
       endif

       !print*, 'Creating ', trim(fname), Atm%npz

       status = nf_def_dim(BCfileid, 'grid_xt'//gn, npx-1, grid_xtdimid) 
       if (status .ne. NF_NOERR) call handle_err(status)
       status = nf_def_dim(BCfileid, 'grid_yt'//gn, npy+2*ng-1, grid_ytdimid)
       if (status .ne. NF_NOERR) call handle_err(status)
       status = nf_def_dim(BCfileid, 'grid_xtstag'//gn, npx, grid_xtstagdimid) 
       if (status .ne. NF_NOERR) call handle_err(status)
       status = nf_def_dim(BCfileid, 'grid_ytstag'//gn, npy+2*ng, grid_ytstagdimid)
       if (status .ne. NF_NOERR) call handle_err(status)
       status = nf_def_dim(BCfileid, 'halo'//gn, ng, haloid)
       if (status .ne. NF_NOERR) call handle_err(status)
       status = nf_def_dim(BCfileid, 'pfull'//gn, Atm%npz, pfullid)
       if (status .ne. NF_NOERR) call handle_err(status)

    endif

    call fv_io_write_a_BC_2D(BCfileid, Atm%grid_number, 'phis', Atm%phis, &
         is, ie, js, je, isd, ied, jsd, jed, npx, npy, 0, 0, gid == masterproc)

    !!! These are to save the nested-grid boundary haloes, to ensure cross-restart consistency. Do not remove!!
    do n=1,Atm%ncnst
       call get_tracer_names(MODEL_ATMOS, n, tname)
       call fv_io_write_a_BC_3D(BCfileid, Atm%grid_number, trim(tname), &
            Atm%q(:,:,:,n), &
            is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 0, gid == masterproc)       
    enddo
#ifndef SW_DYNAMICS
    call fv_io_write_a_BC_3D(BCfileid, Atm%grid_number, 'pt', &
         Atm%pt, &
         is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 0, gid == masterproc)    
    if (.not.Atm%hydrostatic) then
       call fv_io_write_a_BC_3D(BCfileid, Atm%grid_number, 'w', &
            Atm%w, &
            is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 0, gid == masterproc)    
    endif
#endif


    call fv_io_write_a_BC_3D(BCfileid, Atm%grid_number, 'delp', &
         Atm%delp, &
         is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 0, gid == masterproc) 
    call fv_io_write_a_BC_3D(BCfileid, Atm%grid_number, 'v', &
         Atm%v, &
         is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 1, 0, gid == masterproc)
    call fv_io_write_a_BC_3D(BCfileid, Atm%grid_number, 'u', &
         Atm%u, &
         is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 1, gid == masterproc)


    call fv_io_write_a_BC_direct(BCfileid, Atm%grid_number, 'vc', &
         Atm%vc_west_t1, Atm%vc_east_t1, Atm%vc_south_t1, Atm%vc_north_t1, &
         is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 1, gid == masterproc)
    call fv_io_write_a_BC_direct(BCfileid, Atm%grid_number, 'uc', &
         Atm%uc_west_t1, Atm%uc_east_t1, Atm%uc_south_t1, Atm%uc_north_t1, &
         is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 1, 0, gid == masterproc)
    do n=1,Atm%ncnst
       call get_tracer_names(MODEL_ATMOS, n, tname)
       call fv_io_write_a_BC_direct(BCfileid, Atm%grid_number, trim(tname), &
            Atm%q_west(:,:,:,n), Atm%q_east(:,:,:,n), Atm%q_south(:,:,:,n), Atm%q_north(:,:,:,n), &
            is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 0, gid == masterproc)       
    enddo
#ifndef SW_DYNAMICS
    call fv_io_write_a_BC_direct(BCfileid, Atm%grid_number, 'pt', &
         Atm%pt_west, Atm%pt_east, Atm%pt_south, Atm%pt_north, &
         is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 0, gid == masterproc)    
    if (.not.Atm%hydrostatic) then
       call fv_io_write_a_BC_direct(BCfileid, Atm%grid_number, 'w', &
            Atm%w_west, Atm%w_east, Atm%w_south, Atm%w_north, &
            is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 0, gid == masterproc)    
    endif
#endif


    call fv_io_write_a_BC_direct(BCfileid, Atm%grid_number, 'delp', &
         Atm%h_west, Atm%h_east, Atm%h_south, Atm%h_north, &
         is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 0, gid == masterproc) 
    call fv_io_write_a_BC_direct(BCfileid, Atm%grid_number, 'v', &
         Atm%v_west, Atm%v_east, Atm%v_south, Atm%v_north, &
         is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 1, 0, gid == masterproc)
    call fv_io_write_a_BC_direct(BCfileid, Atm%grid_number, 'u', &
         Atm%u_west, Atm%u_east, Atm%u_south, Atm%u_north, &
         is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 1, gid == masterproc)
    call fv_io_write_a_BC_direct(BCfileid, Atm%grid_number, 'vc0', &
         Atm%vc_west_t0, Atm%vc_east_t0, Atm%vc_south_t0, Atm%vc_north_t0, &
         is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 1, gid == masterproc)
    call fv_io_write_a_BC_direct(BCfileid, Atm%grid_number, 'uc0', &
         Atm%uc_west_t0, Atm%uc_east_t0, Atm%uc_south_t0, Atm%uc_north_t0, &
         is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 1, 0, gid == masterproc)

    if (concurrent) then

       call fv_io_write_a_BC_direct(BCfileid, Atm%grid_number, 'delp0', &
            Atm%h_west_t0, Atm%h_east_t0, Atm%h_south_t0, Atm%h_north_t0, &
            is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 0, gid == masterproc) 
       call fv_io_write_a_BC_direct(BCfileid, Atm%grid_number, 'v0', &
            Atm%v_west_t0, Atm%v_east_t0, Atm%v_south_t0, Atm%v_north_t0, &
            is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 1, 0, gid == masterproc)
       call fv_io_write_a_BC_direct(BCfileid, Atm%grid_number, 'u0', &
            Atm%u_west_t0, Atm%u_east_t0, Atm%u_south_t0, Atm%u_north_t0, &
            is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 1, gid == masterproc)
       do n=1,Atm%ncnst
          call get_tracer_names(MODEL_ATMOS, n, tname)
          call fv_io_write_a_BC_direct(BCfileid, Atm%grid_number, trim(tname)//'0', &
               Atm%q_west_t0(:,:,:,n), Atm%q_east_t0(:,:,:,n), Atm%q_south_t0(:,:,:,n), Atm%q_north_t0(:,:,:,n), &
               is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 0, gid == masterproc)       
       enddo
#ifndef SW_DYNAMICS
       call fv_io_write_a_BC_direct(BCfileid, Atm%grid_number, 'pt0', &
            Atm%pt_west_t0, Atm%pt_east_t0, Atm%pt_south_t0, Atm%pt_north_t0, &
            is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 0, gid == masterproc)    
       if (.not.Atm%hydrostatic) then
          call fv_io_write_a_BC_direct(BCfileid, Atm%grid_number, 'w0', &
               Atm%w_west_t0, Atm%w_east_t0, Atm%w_south_t0, Atm%w_north_t0, &
               is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 0, gid == masterproc)    
       endif
#endif

    endif


    if (gid == Atm%pelist(1)) then

       status = nf_close(BCfileid)
       if (status .ne. NF_NOERR) call handle_err(status)

    endif



  end subroutine fv_io_write_BCs

  subroutine fv_io_write_a_BC_3D(BCfileid, gridnum, varname, var, &
       is, ie, js, je, isd, ied, jsd, jed, npx, npy, npz, istag, jstag, master)

    use sim_nc_mod, only: handle_err
    use mpp_domains_mod, only: mpp_get_layout, mpp_get_data_domain, mpp_get_compute_domain, mpp_get_global_domain
    use mpp_mod, only: mpp_npes, mpp_recv, mpp_send, mpp_sync_self, mpp_sync
    use fv_mp_mod, only: ng, mp_gather, gid

#include <netcdf.inc>

    integer, intent(IN) :: BCfileid, gridnum, npz, istag, jstag
    character(len=*), intent(IN) :: varname
    logical, intent(IN) :: master
    integer, intent(IN) :: is, ie, js, je, isd, ied, jsd, jed
    integer, intent(IN) :: npx, npy
    real, allocatable, dimension(:,:,:) :: alldat
    real, intent(IN) :: var(isd:ied+istag,jsd:jed+jstag,npz)
    character(len=3) :: gn
    integer :: status, westid, eastid, southid, northid
    integer :: isd2, ied2, jsd2, jed2
    integer :: k
    
    if (gridnum > 1) then
       write(gn,'(A2, I1)') "_g", gridnum
    else
       gn = ''
    end if

    allocate(alldat(npx+2*ng-1+istag,npy+2*ng-1+jstag,npz))
    alldat = -1000.

    isd2 = is
    jsd2 = js
    ied2 = ie
    jed2 = je
    if (is == 1) isd2 = isd
    if (js == 1) jsd2 = jsd
    if (ie == npx-1) ied2 = ied + istag
    if (je == npy-1) jed2 = jed + jstag

    if (master) then

!!$       print*, 'WRITING ', trim(varname) !! DEBUG

       if (jstag > 0) then
          status = nf_def_var(BCfileid, trim(varname)//'west'//gn, NF_DOUBLE, 3, (/ haloid, grid_ytstagdimid, pfullid /), westid)
          if (status .ne. NF_NOERR) call handle_err(status)
          status = nf_def_var(BCfileid, trim(varname)//'east'//gn, NF_DOUBLE, 3, (/ haloid, grid_ytstagdimid, pfullid /), eastid)
       else
          status = nf_def_var(BCfileid, trim(varname)//'west'//gn, NF_DOUBLE, 3, (/ haloid, grid_ytdimid, pfullid /), westid)
          if (status .ne. NF_NOERR) call handle_err(status)
          status = nf_def_var(BCfileid, trim(varname)//'east'//gn, NF_DOUBLE, 3, (/ haloid, grid_ytdimid, pfullid /), eastid)
       endif
       if (status .ne. NF_NOERR) call handle_err(status)

       if (istag > 0) then
          status = nf_def_var(BCfileid, trim(varname)//'south'//gn, NF_DOUBLE, 3, (/ grid_xtstagdimid, haloid, pfullid /), southid)
          if (status .ne. NF_NOERR) call handle_err(status)
          status = nf_def_var(BCfileid, trim(varname)//'north'//gn, NF_DOUBLE, 3, (/ grid_xtstagdimid, haloid, pfullid /), northid)
       else
          status = nf_def_var(BCfileid, trim(varname)//'south'//gn, NF_DOUBLE, 3, (/ grid_xtdimid, haloid, pfullid /), southid)
          if (status .ne. NF_NOERR) call handle_err(status)
          status = nf_def_var(BCfileid, trim(varname)//'north'//gn, NF_DOUBLE, 3, (/ grid_xtdimid, haloid, pfullid /), northid)
       endif
       if (status .ne. NF_NOERR) call handle_err(status)

       status = nf_enddef(BCfileid)

    endif

       do k=1,npz

          alldat(isd2+ng:ied2+ng,jsd2+ng:jed2+ng,k) = var(isd2:ied2,jsd2:jed2,k)
          call mp_gather(alldat(:,:,k:k), isd2+ng, ied2+ng, jsd2+ng, jed2+ng, npx+2*ng-1+istag, npy+2*ng-1+jstag,1)

          if (master) then

          status = nf_put_vara_double(BCfileid, westid, (/1,1,k/), (/ng, npy+2*ng-1+jstag, 1/), &
               alldat(1:ng,1:npy+2*ng-1+jstag,k))
          if (status .ne. NF_NOERR) call handle_err(status)
          status = nf_put_vara_double(BCfileid, eastid, (/1,1,k/), (/ng, npy+2*ng-1+jstag, 1/), &
               alldat(npx+ng+istag:npx+2*ng-1+istag,1:npy+2*ng-1+jstag,k))
          if (status .ne. NF_NOERR) call handle_err(status)
          status = nf_put_vara_double(BCfileid, southid, (/1,1,k/), (/npx-1+istag, ng, 1/), &
               alldat(ng+1:npx+ng-1+istag,1:ng,k))
          if (status .ne. NF_NOERR) call handle_err(status)
          status = nf_put_vara_double(BCfileid, northid, (/1,1,k/), (/npx-1+istag, ng, 1/), &
               alldat(ng+1:npx+ng-1+istag,npy+ng+jstag:npy+2*ng-1+jstag,k))
          if (status .ne. NF_NOERR) call handle_err(status)

          endif

       enddo

       status = nf_sync(BCfileid)
       status = nf_redef(BCfileid)

!!$       if (master) print*, 'DONE WRITING ', trim(varname) !! DEBUG

    deallocate(alldat)

  end subroutine fv_io_write_a_BC_3D

  subroutine fv_io_write_a_BC_direct(BCfileid, gridnum, varname, BCwest, BCeast, BCsouth, BCnorth, &
       is, ie, js, je, isd, ied, jsd, jed, npx, npy, npz, istag, jstag, master)

    use sim_nc_mod, only: handle_err
    use mpp_domains_mod, only: mpp_get_layout, mpp_get_data_domain, mpp_get_compute_domain, mpp_get_global_domain
    use mpp_mod, only: mpp_npes, mpp_recv, mpp_send, mpp_sync_self, mpp_sync
    use fv_mp_mod, only: ng, mp_gather, gid

#include <netcdf.inc>

    integer, intent(IN) :: BCfileid, gridnum, npz, istag, jstag
    character(len=*), intent(IN) :: varname
    logical, intent(IN) :: master
    integer, intent(IN) :: is, ie, js, je, isd, ied, jsd, jed
    integer, intent(IN) :: npx, npy
    real, dimension(ie+1+istag:ied+istag,jsd:jed+jstag,npz), intent(in) :: BCeast
    real, dimension(isd:is-1,jsd:jed+jstag,npz), intent(in) :: BCwest
    real, dimension(isd:ied+istag,jsd:js-1,npz), intent(in) :: BCsouth
    real, dimension(isd:ied+istag,je+1+jstag:jed+jstag,npz), intent(in) :: BCnorth
    real, allocatable, dimension(:,:,:) :: alldat
    character(len=3) :: gn
    integer :: status, westid, eastid, southid, northid
    integer :: isd2, ied2, jsd2, jed2
    integer :: k
    
    if (gridnum > 1) then
       write(gn,'(A2, I1)') "_g", gridnum
    else
       gn = ''
    end if

    allocate(alldat(npx+2*ng-1+istag,npy+2*ng-1+jstag,npz))
    alldat = -1000.

    isd2 = is
    jsd2 = js
    ied2 = ie
    jed2 = je
    if (is == 1) isd2 = isd
    if (js == 1) jsd2 = jsd
    if (ie == npx-1) ied2 = ied + istag
    if (je == npy-1) jed2 = jed + jstag

    if (master) then

       if (jstag > 0) then
          status = nf_def_var(BCfileid, trim(varname)//'west_BC'//gn, NF_DOUBLE, 3, (/ haloid, grid_ytstagdimid, pfullid /), westid)
          if (status .ne. NF_NOERR) call handle_err(status)
          status = nf_def_var(BCfileid, trim(varname)//'east_BC'//gn, NF_DOUBLE, 3, (/ haloid, grid_ytstagdimid, pfullid /), eastid)
       else
          status = nf_def_var(BCfileid, trim(varname)//'west_BC'//gn, NF_DOUBLE, 3, (/ haloid, grid_ytdimid, pfullid /), westid)
          if (status .ne. NF_NOERR) call handle_err(status)
          status = nf_def_var(BCfileid, trim(varname)//'east_BC'//gn, NF_DOUBLE, 3, (/ haloid, grid_ytdimid, pfullid /), eastid)
       endif
       if (status .ne. NF_NOERR) call handle_err(status)

       if (istag > 0) then
          status = nf_def_var(BCfileid, trim(varname)//'south_BC'//gn, NF_DOUBLE, 3, (/ grid_xtstagdimid, haloid, pfullid /), southid)
          if (status .ne. NF_NOERR) call handle_err(status)
          status = nf_def_var(BCfileid, trim(varname)//'north_BC'//gn, NF_DOUBLE, 3, (/ grid_xtstagdimid, haloid, pfullid /), northid)
       else
          status = nf_def_var(BCfileid, trim(varname)//'south_BC'//gn, NF_DOUBLE, 3, (/ grid_xtdimid, haloid, pfullid /), southid)
          if (status .ne. NF_NOERR) call handle_err(status)
          status = nf_def_var(BCfileid, trim(varname)//'north_BC'//gn, NF_DOUBLE, 3, (/ grid_xtdimid, haloid, pfullid /), northid)
       endif
       if (status .ne. NF_NOERR) call handle_err(status)

       status = nf_enddef(BCfileid)

    endif

       do k=1,npz

          if (is == 1)     alldat(isd+ng:is-1+ng,jsd2+ng:jed2+ng,k)             = BCwest(isd:is-1,jsd2:jed2,k)
          if (ie == npx-1) alldat(ie+1+istag+ng:ied+istag+ng,jsd2+ng:jed2+ng,k) = BCeast(ie+1+istag:ied+istag,jsd2:jed2,k)
          if (js == 1)     alldat(is+ng:ie+istag+ng,jsd+ng:js-1+ng,k)             = BCsouth(is:ie+istag,jsd:js-1,k)
          if (je == npy-1) alldat(is+ng:ie+istag+ng,je+1+jstag+ng:jed+jstag+ng,k) = BCnorth(is:ie+istag,je+1+jstag:jed+jstag,k)

          call mp_gather(alldat(:,:,k:k), isd2+ng, ied2+ng, jsd2+ng, jed2+ng, npx+2*ng-1+istag, npy+2*ng-1+jstag,1)

          if (master) then
          status = nf_put_vara_double(BCfileid, westid, (/1,1,k/), (/ng, npy+2*ng-1+jstag, 1/), &
               alldat(1:ng,1:npy+2*ng-1+jstag,k))
          if (status .ne. NF_NOERR) call handle_err(status)
          status = nf_put_vara_double(BCfileid, eastid, (/1,1,k/), (/ng, npy+2*ng-1+jstag, 1/), &
               alldat(npx+ng+istag:npx+2*ng-1+istag,1:npy+2*ng-1+jstag,k))
          if (status .ne. NF_NOERR) call handle_err(status)
          status = nf_put_vara_double(BCfileid, southid, (/1,1,k/), (/npx-1+istag, ng, 1/), &
               alldat(ng+1:npx+ng-1+istag,1:ng,k))
          if (status .ne. NF_NOERR) call handle_err(status)
          status = nf_put_vara_double(BCfileid, northid, (/1,1,k/), (/npx-1+istag, ng, 1/), &
               alldat(ng+1:npx+ng-1+istag,npy+ng+jstag:npy+2*ng-1+jstag,k))
          if (status .ne. NF_NOERR) call handle_err(status)

          endif

       enddo

       status = nf_sync(BCfileid)
       status = nf_redef(BCfileid)

    deallocate(alldat)

  end subroutine fv_io_write_a_BC_direct

  subroutine fv_io_write_a_BC_2D(BCfileid, gridnum, varname, var, &
       is, ie, js, je, isd, ied, jsd, jed, npx, npy, istag, jstag, master)

    use sim_nc_mod, only: handle_err
    use mpp_domains_mod, only: mpp_get_layout, mpp_get_data_domain, mpp_get_compute_domain, mpp_get_global_domain
    use mpp_mod, only: mpp_npes, mpp_recv, mpp_send, mpp_sync_self, mpp_sync
    use fv_mp_mod, only: ng, mp_gather, gid

#include <netcdf.inc>

    integer, intent(IN) :: BCfileid, gridnum, istag, jstag
    character(len=*), intent(IN) :: varname
    logical, intent(IN) :: master
    integer, intent(IN) :: is, ie, js, je, isd, ied, jsd, jed
    integer, intent(IN) :: npx, npy
    real, allocatable, dimension(:,:,:) :: alldat
    real, intent(IN) :: var(isd:ied+istag,jsd:jed+jstag)
    character(len=3) :: gn
    integer :: status, westid, eastid, southid, northid
    integer :: isd2, ied2, jsd2, jed2
    
    if (gridnum > 1) then
       write(gn,'(A2, I1)') "_g", gridnum
    else
       gn = ''
    end if

    allocate(alldat(npx+2*ng-1+istag,npy+2*ng-1+jstag,1))
    alldat = 0.

    isd2 = is
    jsd2 = js
    ied2 = ie
    jed2 = je
    if (is == 1) isd2 = isd
    if (js == 1) jsd2 = jsd
    if (ie == npx-1) ied2 = ied + istag
    if (je == npy-1) jed2 = jed + jstag

    alldat(isd2+ng:ied2+ng,jsd2+ng:jed2+ng,1) = var(isd2:ied2,jsd2:jed2)

    call mp_gather(alldat, isd2+ng, ied2+ng, jsd2+ng, jed2+ng, npx+2*ng-1+istag, npy+2*ng-1+jstag, 1)

    if (master) then

!!$       print*, 'WRITING ', trim(varname) !! DEBUG

       status = nf_def_var(BCfileid, trim(varname)//'west'//gn, NF_DOUBLE, 2, (/ haloid, grid_ytdimid /), westid)
       if (status .ne. NF_NOERR) call handle_err(status)
       if (jstag > 0) then
          status = nf_def_var(BCfileid, trim(varname)//'east'//gn, NF_DOUBLE, 2, (/ haloid, grid_ytstagdimid /), eastid)
       else
          status = nf_def_var(BCfileid, trim(varname)//'east'//gn, NF_DOUBLE, 2, (/ haloid, grid_ytdimid /), eastid)
       endif
       if (status .ne. NF_NOERR) call handle_err(status)
       status = nf_def_var(BCfileid, trim(varname)//'south'//gn, NF_DOUBLE, 2, (/ grid_xtdimid, haloid /), southid)
       if (status .ne. NF_NOERR) call handle_err(status)
       if (istag > 0) then
          status = nf_def_var(BCfileid, trim(varname)//'north'//gn, NF_DOUBLE, 2, (/ grid_xtstagdimid, haloid /), northid)
       else
          status = nf_def_var(BCfileid, trim(varname)//'north'//gn, NF_DOUBLE, 2, (/ grid_xtdimid, haloid /), northid)
       endif
       if (status .ne. NF_NOERR) call handle_err(status)

       status = nf_enddef(BCfileid)

       status = nf_put_var_double(BCfileid, westid, alldat(1:ng,1:npy+2*ng-1+jstag,1))
       if (status .ne. NF_NOERR) call handle_err(status)
       status = nf_put_var_double(BCfileid, eastid, alldat(npx+ng+istag:npx+2*ng-1+istag,1:npy+2*ng-1+jstag,1))
       if (status .ne. NF_NOERR) call handle_err(status)
       status = nf_put_var_double(BCfileid, southid, alldat(ng+1:npx+ng-1+istag,1:ng,1))
       if (status .ne. NF_NOERR) call handle_err(status)
       status = nf_put_var_double(BCfileid, northid, alldat(ng+1:npx+ng-1+istag,npy+ng+jstag:npy+2*ng-1+jstag,1))
       if (status .ne. NF_NOERR) call handle_err(status)

       status = nf_redef(BCfileid)

!!$       print*, 'DONE WRITING ', trim(varname) !! DEBUG

    end if

    deallocate(alldat)


  end subroutine fv_io_write_a_BC_2D

  subroutine fv_io_read_BCs(Atm)

    use sim_nc_mod, only: handle_err
    use fv_mp_mod, only: masterproc, gid, ng
    use mpp_domains_mod, only: mpp_get_global_domain, CORNER
    use tracer_manager_mod, only: get_tracer_names
    use field_manager_mod,       only: MODEL_ATMOS
    use boundary_mod, only : nested_grid_BC_apply

#include <netcdf.inc>

    type(fv_atmos_type),        intent(inout) :: Atm
    character(len=120) :: fname, tname
    integer :: BCfileid, status
    character(len=3) :: gn
    integer :: npx, npy
    integer :: is, ie, js, je, isd, ied, jsd, jed
    integer :: nxf, nyf, nzf, ngf
    integer :: n

    call mpp_get_global_domain(Atm%domain, xsize = npx, ysize = npy, position=CORNER )
    call mpp_get_data_domain( Atm%domain, isd, ied, jsd, jed )
    call mpp_get_compute_domain( Atm%domain, is, ie, js, je )

    if (gid == masterproc) then


       !Now save the data
    
       if (Atm%grid_number > 1) then
          write(gn,'(A2, I1)') "_g", Atm%grid_number
       else
          gn = ''
       end if
       fname = 'INPUT/fv_BCFile'//gn//'.nc'
       print*, 'Opening ', trim(fname)

       status = nf_open(trim(fname), NF_CLOBBER, BCfileid)
       if (status .ne. NF_NOERR) then
          print*, 'Could not open nested BC file ', trim(fname)
          call handle_err(status)
       endif

       !Check file
  
       status = nf_inq_dimid (BCfileid, 'grid_xt'//gn, grid_xtdimid)
       if (status .ne. NF_NOERR) call handle_err(status)
       status = nf_inq_dimlen (BCfileid, grid_xtdimid, nxf)
       if (status .ne. NF_NOERR) call handle_err(status)
       if (nxf /= npx-1) then
          print*, 'NX = ', nxf, ' NPX = ', npx
          call mpp_error(FATAL, 'NX in BC file does not match npx for the simulation')
       endif
       
       status = nf_inq_dimid (BCfileid, 'grid_yt'//gn, grid_ytdimid)
       if (status .ne. NF_NOERR) call handle_err(status)
       status = nf_inq_dimlen (BCfileid, grid_ytdimid, nyf)
       if (status .ne. NF_NOERR) call handle_err(status)
       if (nyf /= npy-1+2*ng) then
          print*, 'NY = ', nyf, ' NPY = ', npy
          call mpp_error(FATAL, 'NY in BC file does not match npy for the simulation')
       endif

       status = nf_inq_dimid (BCfileid, 'pfull'//gn, pfullid)
       if (status .ne. NF_NOERR) call handle_err(status)
       status = nf_inq_dimlen (BCfileid, pfullid, nzf)
       if (status .ne. NF_NOERR) call handle_err(status)
       if (nzf /= Atm%npz) then
          print*, 'NZ = ', nzf, ' NPZ = ', Atm%npz
          call mpp_error(FATAL, 'NZ in BC file does not match npz for the simulation')
       endif

       status = nf_inq_dimid (BCfileid, 'halo'//gn, haloid)
       if (status .ne. NF_NOERR) call handle_err(status)
       status = nf_inq_dimlen (BCfileid, haloid, ngf)
       if (status .ne. NF_NOERR) call handle_err(status)
       if (ngf /= ng) then
          print*, 'NG = ', ngf, ' ng = ', ng
          call mpp_error(FATAL, 'NG in BC file does not match ng for the simulation')
       endif
       
    end if

    !These read in the actual halo values, for cross-restart consistency. Do not remove!!
    call fv_io_read_a_BC_3D(BCfileid, Atm%grid_number, 'delp', Atm%delp, &
          is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 0, gid == masterproc)
    call fv_io_read_a_BC_3D(BCfileid, Atm%grid_number, 'u', Atm%u, &
          is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 1, gid == masterproc)
    call fv_io_read_a_BC_3D(BCfileid, Atm%grid_number, 'v', Atm%v, &
          is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 1, 0, gid == masterproc)
    do n=1,Atm%ncnst
       call get_tracer_names(MODEL_ATMOS, n, tname)
       call fv_io_read_a_BC_3D(BCfileid, Atm%grid_number, trim(tname), Atm%q(:,:,:,n), &
            is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 0, gid == masterproc)       
    enddo
#ifndef SW_DYNAMICS
    call fv_io_read_a_BC_3D(BCfileid, Atm%grid_number, 'pt', Atm%pt, &
          is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 0, gid == masterproc)    
    if (.not. Atm%hydrostatic ) then
        call fv_io_read_a_BC_3D(BCfileid, Atm%grid_number, 'w', Atm%w, &
             is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 0, gid == masterproc) 
    endif
#endif
    call fv_io_read_a_BC_2D(BCfileid, Atm%grid_number, 'phis', Atm%phis, &
         is, ie, js, je, isd, ied, jsd, jed, npx, npy, 0, 0, gid == masterproc)



    call fv_io_read_a_BC_direct(BCfileid, Atm%grid_number, 'delp', &
         Atm%h_west, Atm%h_east, Atm%h_south, Atm%h_north, &
         is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 0, gid == masterproc) 
    call fv_io_read_a_BC_direct(BCfileid, Atm%grid_number, 'v', &
         Atm%v_west, Atm%v_east, Atm%v_south, Atm%v_north, &
         is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 1, 0, gid == masterproc)
    call fv_io_read_a_BC_direct(BCfileid, Atm%grid_number, 'u', &
         Atm%u_west, Atm%u_east, Atm%u_south, Atm%u_north, &
         is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 1, gid == masterproc)
    call fv_io_read_a_BC_direct(BCfileid, Atm%grid_number, 'vc0', &
         Atm%vc_west_t0, Atm%vc_east_t0, Atm%vc_south_t0, Atm%vc_north_t0, &
         is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 1, gid == masterproc)
    call fv_io_read_a_BC_direct(BCfileid, Atm%grid_number, 'uc0', &
         Atm%uc_west_t0, Atm%uc_east_t0, Atm%uc_south_t0, Atm%uc_north_t0, &
         is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 1, 0, gid == masterproc)
    call fv_io_read_a_BC_direct(BCfileid, Atm%grid_number, 'vc', &
         Atm%vc_west_t1, Atm%vc_east_t1, Atm%vc_south_t1, Atm%vc_north_t1, &
         is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 1, gid == masterproc)
    call fv_io_read_a_BC_direct(BCfileid, Atm%grid_number, 'uc', &
         Atm%uc_west_t1, Atm%uc_east_t1, Atm%uc_south_t1, Atm%uc_north_t1, &
         is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 1, 0, gid == masterproc)
    do n=1,Atm%ncnst
       call get_tracer_names(MODEL_ATMOS, n, tname)
       call fv_io_read_a_BC_direct(BCfileid, Atm%grid_number, trim(tname), &
            Atm%q_west(:,:,:,n), Atm%q_east(:,:,:,n), Atm%q_south(:,:,:,n), Atm%q_north(:,:,:,n), &
            is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 0, gid == masterproc)       
    enddo
#ifndef SW_DYNAMICS
    call fv_io_read_a_BC_direct(BCfileid, Atm%grid_number, 'pt', &
         Atm%pt_west, Atm%pt_east, Atm%pt_south, Atm%pt_north, &
         is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 0, gid == masterproc)    
    if (.not.Atm%hydrostatic) then
       call fv_io_read_a_BC_direct(BCfileid, Atm%grid_number, 'w', &
            Atm%w_west, Atm%w_east, Atm%w_south, Atm%w_north, &
            is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 0, gid == masterproc)    
    endif
#endif

    if (concurrent) then
       call fv_io_read_a_BC_direct(BCfileid, Atm%grid_number, 'delp0', &
            Atm%h_west_t0, Atm%h_east_t0, Atm%h_south_t0, Atm%h_north_t0, &
            is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 0, gid == masterproc) 
       call fv_io_read_a_BC_direct(BCfileid, Atm%grid_number, 'v0', &
            Atm%v_west_t0, Atm%v_east_t0, Atm%v_south_t0, Atm%v_north_t0, &
            is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 1, 0, gid == masterproc)
       call fv_io_read_a_BC_direct(BCfileid, Atm%grid_number, 'u0', &
            Atm%u_west_t0, Atm%u_east_t0, Atm%u_south_t0, Atm%u_north_t0, &
            is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 1, gid == masterproc)
       do n=1,Atm%ncnst
          call get_tracer_names(MODEL_ATMOS, n, tname)
          call fv_io_read_a_BC_direct(BCfileid, Atm%grid_number, trim(tname)//'0', &
               Atm%q_west_t0(:,:,:,n), Atm%q_east_t0(:,:,:,n), Atm%q_south_t0(:,:,:,n), Atm%q_north_t0(:,:,:,n), &
               is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 0, gid == masterproc)       
       enddo
#ifndef SW_DYNAMICS
       call fv_io_read_a_BC_direct(BCfileid, Atm%grid_number, 'pt0', &
            Atm%pt_west_t0, Atm%pt_east_t0, Atm%pt_south_t0, Atm%pt_north_t0, &
            is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 0, gid == masterproc)    
       if (.not.Atm%hydrostatic) then
          call fv_io_read_a_BC_direct(BCfileid, Atm%grid_number, 'w0', &
               Atm%w_west_t0, Atm%w_east_t0, Atm%w_south_t0, Atm%w_north_t0, &
               is, ie, js, je, isd, ied, jsd, jed, npx, npy, Atm%npz, 0, 0, gid == masterproc)    
       endif
#endif


    endif

  end subroutine fv_io_read_BCs
    
  subroutine fv_io_read_a_BC_3D(BCfileid, gridnum, varname, var, &
       is, ie, js, je, isd, ied, jsd, jed, npx, npy, npz, istag, jstag, master)

    use sim_nc_mod, only: get_var3_double, get_var2_double
    use mpp_mod, only: mpp_broadcast
    use fv_mp_mod, only: masterproc, ng

#include <netcdf.inc>

    integer, intent(IN) :: BCfileid, gridnum, npz, istag, jstag
    character(len=*), intent(IN) :: varname
    logical, intent(IN) :: master
    integer, intent(IN) :: is, ie, js, je, isd, ied, jsd, jed
    integer, intent(IN) :: npx, npy
    real, intent(INOUT) :: var(isd:ied+istag,jsd:jed+jstag,npz)
    real :: varwest(-ng+1:0,-ng+1:npy+ng-1+jstag,npz)
    real :: vareast(npx+istag:npx+istag+ng-1,-ng+1:npy+ng-1+jstag,npz)
    real :: varsouth(1:npx-1+istag,-ng+1:0,npz)
    real :: varnorth(1:npx-1+istag,npy+jstag:npy+jstag+ng-1,npz)
    character(len=3) :: gn
    integer :: status, i, j, k
    integer is2, ie2, js2, je2
    
    if (gridnum > 1) then
       write(gn,'(A2, I1)') "_g", gridnum
    else
       gn = ''
    end if

    varwest = -1.e25
    vareast = -1.e25
    varsouth = -1.e25
    varnorth = -1.e25

    !This is rather simpler than collecting the BCs. What we can do here is just
    !load the BCs on the head PE, and then broadcast the data to every PE.

    if (master) then

       call get_var3_double( BCfileid, varname//'west'//gn, ng, npy+2*ng-1+jstag, npz, varwest)
       call get_var3_double( BCfileid, varname//'east'//gn, ng, npy+2*ng-1+jstag, npz, vareast)
       call get_var3_double( BCfileid, varname//'north'//gn, npx-1+istag, ng, npz, varnorth)
       call get_var3_double( BCfileid, varname//'south'//gn, npx-1+istag, ng, npz, varsouth)

    endif

    call mpp_broadcast(varwest,  size(varwest),  masterproc)
    call mpp_broadcast(vareast,  size(vareast),  masterproc)
    call mpp_broadcast(varnorth, size(varnorth), masterproc)
    call mpp_broadcast(varsouth, size(varsouth), masterproc)

    is2 = max(1,isd)
    ie2 = min(npx-1,ied)+istag

    if (is == 1) then
       do k=1,npz
       do j=jsd,jed+jstag
       do i=isd,0
          var(i,j,k) = varwest(i,j,k)
       enddo
       enddo
       enddo
    endif
    if (ie == npx-1) then
       do k=1,npz
       do j=jsd,jed+jstag
       do i=npx+istag,ied+istag
          var(i,j,k) = vareast(i,j,k)
       enddo
       enddo
       enddo
    endif

    if (js == 1) then
       do k=1,npz
       do j=jsd,0
       do i=is2,ie2
          var(i,j,k) = varsouth(i,j,k)
       enddo
       enddo
       enddo
    endif
    if (je == npy-1) then
       do k=1,npz
       do j=npy+jstag,jed+jstag
       do i=is2,ie2
          var(i,j,k) = varnorth(i,j,k)
       enddo
       enddo
       enddo
    endif

  end subroutine fv_io_read_a_BC_3D
    
  subroutine fv_io_read_a_BC_direct(BCfileid, gridnum, varname, BCwest, BCeast, BCsouth, BCnorth, &
       is, ie, js, je, isd, ied, jsd, jed, npx, npy, npz, istag, jstag, master)

    use sim_nc_mod, only: get_var3_double, get_var2_double
    use mpp_mod, only: mpp_broadcast
    use fv_mp_mod, only: masterproc, ng

#include <netcdf.inc>

    integer, intent(IN) :: BCfileid, gridnum, npz, istag, jstag
    character(len=*), intent(IN) :: varname
    logical, intent(IN) :: master
    integer, intent(IN) :: is, ie, js, je, isd, ied, jsd, jed
    integer, intent(IN) :: npx, npy
    real, dimension(ie+1+istag:ied+istag,jsd:jed+jstag,npz), intent(out) :: BCeast
    real, dimension(isd:is-1,jsd:jed+jstag,npz), intent(out) :: BCwest
    real, dimension(isd:ied+istag,jsd:js-1,npz), intent(out) :: BCsouth
    real, dimension(isd:ied+istag,je+1+jstag:jed+jstag,npz), intent(out) :: BCnorth

    real :: varwest(-ng+1:0,-ng+1:npy+ng-1+jstag,npz)
    real :: vareast(npx+istag:npx+istag+ng-1,-ng+1:npy+ng-1+jstag,npz)
    real :: varsouth(1:npx-1+istag,-ng+1:0,npz)
    real :: varnorth(1:npx-1+istag,npy+jstag:npy+jstag+ng-1,npz)
    character(len=3) :: gn
    integer :: status, i, j, k
    integer is2, ie2, js2, je2
    
    if (gridnum > 1) then
       write(gn,'(A2, I1)') "_g", gridnum
    else
       gn = ''
    end if

    varwest = -1.e25
    vareast = -1.e25
    varsouth = -1.e25
    varnorth = -1.e25

    !This is rather simpler than collecting the BCs. What we can do here is just
    !load the BCs on the head PE, and then broadcast the data to every PE.

    if (master) then

       call get_var3_double( BCfileid, varname//'west_BC'//gn, ng, npy+2*ng-1+jstag, npz, varwest)
       call get_var3_double( BCfileid, varname//'east_BC'//gn, ng, npy+2*ng-1+jstag, npz, vareast)
       call get_var3_double( BCfileid, varname//'north_BC'//gn, npx-1+istag, ng, npz, varnorth)
       call get_var3_double( BCfileid, varname//'south_BC'//gn, npx-1+istag, ng, npz, varsouth)

    endif

    call mpp_broadcast(varwest,  size(varwest),  masterproc)
    call mpp_broadcast(vareast,  size(vareast),  masterproc)
    call mpp_broadcast(varnorth, size(varnorth), masterproc)
    call mpp_broadcast(varsouth, size(varsouth), masterproc)

    is2 = max(1,isd)
    ie2 = min(npx-1,ied)+istag

    if (is == 1) then
       do k=1,npz
       do j=jsd,jed+jstag
       do i=isd,0
          BCwest(i,j,k) = varwest(i,j,k)
       enddo
       enddo
       enddo
    endif
    if (ie == npx-1) then
       do k=1,npz
       do j=jsd,jed+jstag
       do i=npx+istag,ied+istag
          BCeast(i,j,k) = vareast(i,j,k)
       enddo
       enddo
       enddo
    endif

    if (js == 1) then
       do k=1,npz
       do j=jsd,0
       do i=is2,ie2
          BCsouth(i,j,k) = varsouth(i,j,k)
       enddo
       enddo
       enddo
    endif
    if (je == npy-1) then
       do k=1,npz
       do j=npy+jstag,jed+jstag
       do i=is2,ie2
          BCnorth(i,j,k) = varnorth(i,j,k)
       enddo
       enddo
       enddo
    endif

  end subroutine fv_io_read_a_BC_direct


  subroutine fv_io_read_a_BC_2D(BCfileid, gridnum, varname, var, &
       is, ie, js, je, isd, ied, jsd, jed, npx, npy, istag, jstag, master)

    use sim_nc_mod, only: get_var3_double, get_var2_double
    use mpp_mod, only: mpp_broadcast
    use fv_mp_mod, only: masterproc, ng

#include <netcdf.inc>

    integer, intent(IN) :: BCfileid, gridnum, istag, jstag
    character(len=*), intent(IN) :: varname
    logical, intent(IN) :: master
    integer, intent(IN) :: is, ie, js, je, isd, ied, jsd, jed
    integer, intent(IN) :: npx, npy
    real, intent(INOUT) :: var(isd:ied+istag,jsd:jed+jstag)
    real :: varwest(-ng+1:0,-ng+1:npy+ng-1+jstag)
    real :: vareast(npx+istag:npx+istag+ng-1,-ng+1:npy+ng-1+jstag)
    real :: varsouth(1:npx-1+istag,-ng+1:0)
    real :: varnorth(1:npx-1+istag,npy+jstag:npy+jstag+ng-1)
    character(len=3) :: gn
    integer :: status, i, j, k
    integer is2, ie2, js2, je2
    
    if (gridnum > 1) then
       write(gn,'(A2, I1)') "_g", gridnum
    else
       gn = ''
    end if


    !This is rather simpler than collecting the BCs. What we can do here is just
    !load the BCs on the head PE, and then broadcast the data to every PE.

    if (master) then

       call get_var2_double( BCfileid, varname//'west'//gn, ng, npy+2*ng-1+jstag, varwest)
       call get_var2_double( BCfileid, varname//'east'//gn, ng, npy+2*ng-1+jstag, vareast)
       call get_var2_double( BCfileid, varname//'north'//gn, npx-1+istag, ng, varnorth)
       call get_var2_double( BCfileid, varname//'south'//gn, npx-1+istag, ng, varsouth)

    endif

    call mpp_broadcast(varwest,  size(varwest),  masterproc)
    call mpp_broadcast(vareast,  size(vareast),  masterproc)
    call mpp_broadcast(varnorth, size(varnorth), masterproc)
    call mpp_broadcast(varsouth, size(varsouth), masterproc)

    is2 = max(1,isd)
    ie2 = min(npx-1,ied)+istag

    if (is == 1) then
       do j=jsd,jed+jstag
       do i=isd,0
          var(i,j) = varwest(i,j)
       enddo
       enddo
    endif
    if (ie == npx-1) then
       do j=jsd,jed+jstag
       do i=npx+istag,ied+istag
          var(i,j) = vareast(i,j)
       enddo
       enddo
    endif

    if (js == 1) then
       do j=jsd,0
       do i=is2,ie2
          var(i,j) = varsouth(i,j)
       enddo
       enddo
    endif
    if (je == npy-1) then
       do j=npy+jstag,jed+jstag
       do i=is2,ie2
          var(i,j) = varnorth(i,j)
       enddo
       enddo
    endif

  end subroutine fv_io_read_a_BC_2D

end module fv_io_mod
