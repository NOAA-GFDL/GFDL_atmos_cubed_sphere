module fv_restart_mod
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

  use constants_mod,       only: kappa, pi, omega, rdgas, grav
  use fv_arrays_mod,       only: fv_atmos_type
  use fv_io_mod,           only: fv_io_init, fv_io_read_restart, fv_io_write_restart, &
                                 remap_restart, fv_io_register_restart
  use fv_grid_tools_mod,   only: area, dx, dy, rdxa, rdya
  use fv_grid_utils_mod,   only: fc, f0, ptop, ptop_min, fill_ghost, big_number,   &
                                 make_eta_level, deglat, cubed_to_latlon
  use fv_diagnostics_mod,  only: prt_maxmin
  use init_hydro_mod,      only: p_var
  use mpp_domains_mod,     only: mpp_update_domains, domain2d, DGRID_NE
  use mpp_mod,             only: mpp_chksum, stdout, mpp_error, FATAL, get_unit
  use test_cases_mod,      only: alpha, init_case, init_double_periodic, init_latlon
  use fv_mp_mod,           only: gid, masterproc
  use fv_surf_map_mod,     only: sgh_g, oro_g
  use fv_diagnostics_mod,  only: steps, efx, efx_sum, mtq, mtq_sum
  use tracer_manager_mod,  only: get_tracer_names
  use field_manager_mod,   only: MODEL_ATMOS
  use external_ic_mod,     only: get_external_ic
  use fv_eta_mod,          only: compute_dz_L32, set_hybrid_z


  implicit none
  private

  public :: fv_restart_init, fv_restart_end, fv_restart, fv_write_restart

  !--- private data type
  logical                       :: module_is_initialized = .FALSE.

  !--- version information variables ----
  character(len=128) :: version = '$Id: fv_restart.F90,v 18.0 2010/03/02 23:27:40 fms Exp $'
  character(len=128) :: tagname = '$Name: riga $'

contains 

  !#####################################################################
  ! <SUBROUTINE NAME="fv_restart_init">
  !
  ! <DESCRIPTION>
  ! Initialize the fv core restart facilities
  ! </DESCRIPTION>
  !
  subroutine fv_restart_init()
    call fv_io_init()
    module_is_initialized = .TRUE.
  end subroutine fv_restart_init
  ! </SUBROUTINE> NAME="fv_restart_init"


    !#####################################################################
  ! <SUBROUTINE NAME="fv_restart">
  !
  ! <DESCRIPTION>
  ! The fv core restart facility
  ! </DESCRIPTION>
  !
  subroutine fv_restart(fv_domain, Atm, dt_atmos, seconds, days, cold_start, grid_type)
    type(domain2d),      intent(inout) :: fv_domain
    type(fv_atmos_type), intent(inout) :: Atm(:)
    real,                intent(in)    :: dt_atmos
    integer,             intent(out)   :: seconds
    integer,             intent(out)   :: days
    logical,             intent(in)    :: cold_start
    integer,             intent(in)    :: grid_type


    integer :: i, j, k, n, ntileMe
    integer :: isc, iec, jsc, jec, npz, npz_rst, ncnst
    integer :: isd, ied, jsd, jed

    integer :: unit
    real, allocatable :: dz1(:)
    real rgrav, f00, ztop
    logical :: hybrid
    character(len=128):: tname

    rgrav = 1. / grav

    if(.not.module_is_initialized) call mpp_error(FATAL, 'You must call fv_restart_init.')

    ntileMe = size(Atm(:))
    npz     = Atm(1)%npz
    npz_rst = Atm(1)%npz_rst

    !--- call fv_io_register_restart to register restart field to be written out in fv_io_write_restart
    call fv_io_register_restart(fv_domain,Atm)

    if( .not.cold_start .and. (.not. Atm(1)%external_ic) ) then
        if ( npz_rst /= 0 .and. npz_rst /= npz ) then
!            Remap vertically the prognostic variables for the chosen vertical resolution
             if( gid==masterproc ) then
                 write(*,*) ' '
                 write(*,*) '***** Important Note from FV core ********************'
                 write(*,*) 'Remapping dynamic IC from', npz_rst, 'levels to ', npz,'levels'
                 write(*,*) '***** End Note from FV core **************************'
                 write(*,*) ' '
             endif
             call remap_restart( fv_domain, Atm )
             if( gid==masterproc ) write(*,*) 'Done remapping dynamical IC'
        else
             call fv_io_read_restart(fv_domain,Atm)
        endif
    endif

!---------------------------------------------------------------------------------------------
! Read, interpolate (latlon to cubed), then remap vertically with terrain adjustment if needed
!---------------------------------------------------------------------------------------------
    if ( Atm(1)%external_ic ) then
         call get_external_ic(Atm, fv_domain) 
         if( gid==masterproc ) write(*,*) 'IC generated from the specified external source'
    endif

    seconds = 0; days = 0   ! Restart needs to be modified to record seconds and days.

! Notes by Jeff D.
  ! This logic doesn't work very well.
  ! Shouldn't have read for all tiles then loop over tiles

    do n = 1, ntileMe

       isd = Atm(n)%isd
       ied = Atm(n)%ied
       jsd = Atm(n)%jsd
       jed = Atm(n)%jed
       ncnst = Atm(n)%ncnst
       isc = Atm(n)%isc; iec = Atm(n)%iec; jsc = Atm(n)%jsc; jec = Atm(n)%jec

      ! Init model data
      if(.not.cold_start)then  ! This is not efficient stacking if there are really more tiles than 1.
        if ( Atm(n)%mountain ) then
             call mpp_update_domains( Atm(n)%phis, fv_domain, complete=.true. )
        else
             Atm(n)%phis = 0.
            if( gid==masterproc ) write(*,*) 'phis set to zero'
        endif

#ifdef SW_DYNAMICS
        Atm(n)%pt(:,:,:)=1.
#else
        if (ptop/=Atm(n)%ak(1)) call mpp_error(FATAL,'FV restart: ptop not equal Atm(n)%ak(1)')
        call p_var(npz,         isc,         iec,       jsc,     jec,   ptop,     ptop_min,  &
                   Atm(n)%delp, Atm(n)%delz, Atm(n)%pt, Atm(n)%ps, Atm(n)%pe, Atm(n)%peln,   &
                   Atm(n)%pk,   Atm(n)%pkz, kappa, Atm(n)%q, Atm(n)%ng, ncnst,  Atm(n)%dry_mass,  &
                   Atm(n)%adjust_dry_mass,  Atm(n)%mountain, Atm(n)%moist_phys,  Atm(n)%hydrostatic, &
                   Atm(n)%k_top, Atm(n)%nwat, Atm(n)%make_nh)
#endif
        if ( grid_type < 7 .and. grid_type /= 4 ) then
! Fill big values in the non-existinng corner regions:
!          call fill_ghost(Atm(n)%phis, Atm(n)%npx, Atm(n)%npy, big_number)
           do j=jsd,jed+1
           do i=isd,ied+1
              fc(i,j) = 2.*omega*( -cos(Atm(n)%grid(i,j,1))*cos(Atm(n)%grid(i,j,2))*sin(alpha) + &
                                    sin(Atm(n)%grid(i,j,2))*cos(alpha) )
           enddo
           enddo
           do j=jsd,jed
           do i=isd,ied
             f0(i,j) = 2.*omega*( -cos(Atm(n)%agrid(i,j,1))*cos(Atm(n)%agrid(i,j,2))*sin(alpha) + &
                                    sin(Atm(n)%agrid(i,j,2))*cos(alpha) )
           enddo
           enddo
        else
           f00 = 2.*omega*sin(deglat/180.*pi)
           do j=jsd,jed+1
              do i=isd,ied+1
                 fc(i,j) = f00
              enddo
           enddo
           do j=jsd,jed
              do i=isd,ied
                 f0(i,j) = f00
              enddo
           enddo
        endif
      else
       if ( Atm(n)%warm_start ) then
         call mpp_error(FATAL, 'FV restart files not found; set warm_start = .F. if cold_start is desired.')
       endif
! Cold start
       if ( Atm(n)%make_hybrid_z ) then
         hybrid = .false.
       else
         hybrid = Atm(n)%hybrid_z
       endif
         if (grid_type < 4) then
            if ( .not. Atm(n)%external_ic ) then
            call init_case(Atm(n)%u,Atm(n)%v,Atm(n)%w,Atm(n)%pt,Atm(n)%delp,Atm(n)%q,Atm(n)%phis, Atm(n)%ps,Atm(n)%pe, &
                           Atm(n)%peln,Atm(n)%pk,Atm(n)%pkz, Atm(n)%uc,Atm(n)%vc, Atm(n)%ua,Atm(n)%va,        & 
                           Atm(n)%ak, Atm(n)%bk, Atm(n)%npx, Atm(n)%npy, npz, Atm(n)%ng, ncnst, Atm(n)%nwat,  &
                           Atm(n)%k_top, Atm(n)%ndims, Atm(n)%ntiles, Atm(n)%dry_mass, Atm(n)%mountain,       &
                           Atm(n)%moist_phys, Atm(n)%hydrostatic, hybrid, Atm(n)%delz, Atm(n)%ze0)
            endif
         elseif (grid_type == 4) then
            call init_double_periodic(Atm(n)%u,Atm(n)%v,Atm(n)%w,Atm(n)%pt,Atm(n)%delp,Atm(n)%q,Atm(n)%phis, Atm(n)%ps,Atm(n)%pe, &
                                      Atm(n)%peln,Atm(n)%pk,Atm(n)%pkz, Atm(n)%uc,Atm(n)%vc, Atm(n)%ua,Atm(n)%va,        & 
                                      Atm(n)%ak, Atm(n)%bk, Atm(n)%npx, Atm(n)%npy, npz, Atm(n)%ng, ncnst, Atm(n)%nwat,  &
                                      Atm(n)%k_top, Atm(n)%ndims, Atm(n)%ntiles, Atm(n)%dry_mass, Atm(n)%mountain,       &
                                      Atm(n)%moist_phys, Atm(n)%hydrostatic, hybrid, Atm(n)%delz, Atm(n)%ze0)
            if( gid==masterproc ) write(*,*) 'Doubly Periodic IC generated'
         elseif (grid_type == 5 .or. grid_type == 6) then
            call init_latlon(Atm(n)%u,Atm(n)%v,Atm(n)%pt,Atm(n)%delp,Atm(n)%q,Atm(n)%phis, Atm(n)%ps,Atm(n)%pe, &
                             Atm(n)%peln,Atm(n)%pk,Atm(n)%pkz, Atm(n)%uc,Atm(n)%vc, Atm(n)%ua,Atm(n)%va,        &
                             Atm(n)%ak, Atm(n)%bk, Atm(n)%npx, Atm(n)%npy, npz, Atm(n)%ng, ncnst, &
                             Atm(n)%k_top, Atm(n)%ndims, Atm(n)%ntiles, Atm(n)%dry_mass, Atm(n)%mountain,       &
                             Atm(n)%moist_phys, hybrid, Atm(n)%delz, Atm(n)%ze0)
         endif

        if ( Atm(n)%fv_land ) then
             do j=jsc,jec
                do i=isc,iec
                   Atm(n)%sgh(i,j) = sgh_g(i,j)
                   Atm(n)%oro(i,j) = oro_g(i,j)
                enddo
             enddo
        endif

        if ( Atm(n)%no_cgrid ) then
             Atm(n)%um(:,:,:) = Atm(n)%u(:,:,:)
             Atm(n)%vm(:,:,:) = Atm(n)%v(:,:,:)
        endif

      endif  !end cold_start check

!---------------------------------------------------------------------------------------------
! Transform the (starting) Eulerian vertical coordinate from sigma-p to hybrid_z
     if ( Atm(n)%hybrid_z ) then
       if ( Atm(n)%make_hybrid_z ) then
          allocate ( dz1(npz) )
          if( npz==32 ) then
              call compute_dz_L32(npz, ztop, dz1)
          else
              call mpp_error(FATAL, 'You must provide a specific routine for hybrid_z')
          endif
          call set_hybrid_z(isc, iec, jsc, jec, Atm(n)%ng, npz, ztop, dz1, rgrav,  &
                            Atm(n)%phis, Atm(n)%ze0)
          deallocate ( dz1 )
!         call prt_maxmin('ZE0', Atm(n)%ze0,  isc, iec, jsc, jec, 0, npz, 1.E-3, gid==masterproc)
!         call prt_maxmin('DZ0', Atm(n)%delz, isc, iec, jsc, jec, 0, npz, 1.   , gid==masterproc)
       endif
       call make_eta_level(npz, Atm(n)%pe, area, Atm(n)%ks, Atm(n)%ak, Atm(n)%bk)
      endif
!---------------------------------------------------------------------------------------------

      unit = stdout()
      write(unit,*)
      write(unit,*) 'fv_restart u    = ', mpp_chksum(Atm(n)%u(isc:iec,jsc:jec,:))
      write(unit,*) 'fv_restart v    = ', mpp_chksum(Atm(n)%v(isc:iec,jsc:jec,:))
      write(unit,*) 'fv_restart delp = ', mpp_chksum(Atm(n)%delp(isc:iec,jsc:jec,:))
      write(unit,*) 'fv_restart phis = ', mpp_chksum(Atm(n)%phis(isc:iec,jsc:jec))

#ifdef SW_DYNAMICS
      call prt_maxmin('H ', Atm(n)%delp, isc, iec, jsc, jec, Atm(n)%ng, 1, rgrav, gid==masterproc)
#else
      write(unit,*) 'fv_restart pt   = ', mpp_chksum(Atm(n)%pt(isc:iec,jsc:jec,:))
      if (ncnst>0) write(unit,*) 'fv_init nq =',ncnst, mpp_chksum(Atm(n)%q(isc:iec,jsc:jec,:,:))
!---------------
! Check Min/Max:
!---------------
      call prt_maxmin('ZS', Atm(n)%phis, isc, iec, jsc, jec, Atm(n)%ng, 1, rgrav, gid==masterproc)
!     call prt_maxmin('ORO',Atm(n)%oro, isc, iec, jsc, jec,          0, 1, 1., gid==masterproc)

      if ( (.not.Atm(n)%hydrostatic) .and. (.not.Atm(n)%make_nh) ) then
            call prt_maxmin('DZ', Atm(n)%delz, isc, iec, jsc, jec, 0, npz, 1., gid==masterproc)
            if ( Atm(n)%hybrid_z ) then
            call prt_maxmin('ZTOP(km)', Atm(n)%ze0, isc, iec, jsc, jec, 0, 1, 1.E-3, gid==masterproc)
            call prt_maxmin('DZ_top', Atm(n)%delz, isc, iec, jsc, jec, 0, 1, 1.E-3, gid==masterproc)
            endif
      endif

      call prt_maxmin('PS', Atm(n)%ps, isc, iec, jsc, jec, Atm(n)%ng, 1,    0.01, gid==masterproc)
      call prt_maxmin('T ', Atm(n)%pt, isc, iec, jsc, jec, Atm(n)%ng, npz, 1., gid==masterproc)

! Check tracers:
      do i=1, ncnst
          call get_tracer_names ( MODEL_ATMOS, i, tname )
          call prt_maxmin(trim(tname), Atm(n)%q(isd:ied,jsd:jed,1:npz,i), isc, iec, jsc, jec, Atm(n)%ng, npz, 1.,gid==masterproc)
      enddo
#endif
      call prt_maxmin('U ', Atm(n)%u(isc:iec,jsc:jec,1:npz), isc, iec, jsc, jec, 0, npz, 1., gid==masterproc)
      call prt_maxmin('V ', Atm(n)%v(isc:iec,jsc:jec,1:npz), isc, iec, jsc, jec, 0, npz, 1., gid==masterproc)
      if ( .not.Atm(n)%hydrostatic )   &
      call prt_maxmin('W ', Atm(n)%w, isc, iec, jsc, jec, Atm(n)%ng, npz, 1.,gid==masterproc)

      if ( (.not.Atm(n)%hydrostatic) .and. Atm(n)%make_nh ) then
         Atm(n)%w = 0.
         if ( .not.Atm(n)%hybrid_z ) then
             do k=1,npz
                do j=jsc,jec
                   do i=isc,iec
                      Atm(n)%delz(i,j,k) = (rdgas*rgrav)*Atm(n)%pt(i,j,k)*(Atm(n)%peln(i,k,j)-Atm(n)%peln(i,k+1,j))
                   enddo
                enddo
             enddo
         endif
      endif

      if (gid==masterproc) write(unit,*)

!--------------------------------------------
! Initialize surface winds for flux coupler:
!--------------------------------------------
    if ( .not. Atm(n)%srf_init ) then
         call cubed_to_latlon(Atm(n)%u, Atm(n)%v, Atm(n)%ua, Atm(n)%va, dx, dy, rdxa, rdya, npz, 1)
         do j=jsc,jec
            do i=isc,iec
               Atm(n)%u_srf(i,j) = Atm(n)%ua(i,j,npz)
               Atm(n)%v_srf(i,j) = Atm(n)%va(i,j,npz)
            enddo
         enddo
         Atm(n)%srf_init = .true.
    endif

    end do   ! n_tile

  end subroutine fv_restart
  ! </SUBROUTINE> NAME="fv_restart"



  !#######################################################################
  ! <SUBROUTINE NAME="fv_write_restart">
  ! <DESCRIPTION>
  !  Write out restart files registered through register_restart_file
  ! </DESCRIPTION>
  subroutine fv_write_restart(Atm, timestamp)
    type(fv_atmos_type), intent(inout) :: Atm(:)
    character(len=*),    intent(in)    :: timestamp

    call fv_io_write_restart(Atm, timestamp)

  end subroutine fv_write_restart
  ! </SUBROUTINE>



  !#####################################################################
  ! <SUBROUTINE NAME="fv_restart_end">
  !
  ! <DESCRIPTION>
  ! Initialize the fv core restart facilities
  ! </DESCRIPTION>
  !
  subroutine fv_restart_end(Atm)
    type(fv_atmos_type), intent(in) :: Atm(:)

    integer :: isc, iec, jsc, jec
    integer :: iq, n, ntileMe
    integer :: isd, ied, jsd, jed, npz
    integer :: unit
    integer :: file_unit

    ntileMe = size(Atm(:))

    do n = 1, ntileMe
      isc = Atm(n)%isc; iec = Atm(n)%iec; jsc = Atm(n)%jsc; jec = Atm(n)%jec

      isd = Atm(n)%isd
      ied = Atm(n)%ied
      jsd = Atm(n)%jsd
      jed = Atm(n)%jed
      npz = Atm(n)%npz
 
      unit = stdout()
      write(unit,*)
      write(unit,*) 'fv_restart_end u    = ', mpp_chksum(Atm(n)%u(isc:iec,jsc:jec,:))
      write(unit,*) 'fv_restart_end v    = ', mpp_chksum(Atm(n)%v(isc:iec,jsc:jec,:))
      if ( .not. Atm(n)%hydrostatic )    &
         write(unit,*) 'fv_restart_end w    = ', mpp_chksum(Atm(n)%w(isc:iec,jsc:jec,:))
      write(unit,*) 'fv_restart_end delp = ', mpp_chksum(Atm(n)%delp(isc:iec,jsc:jec,:))
      write(unit,*) 'fv_restart_end phis = ', mpp_chksum(Atm(n)%phis(isc:iec,jsc:jec))
#ifndef SW_DYNAMICS
      write(unit,*) 'fv_restart_end pt   = ', mpp_chksum(Atm(n)%pt(isc:iec,jsc:jec,:))
      do iq=1,min(7, Atm(n)%ncnst)     ! Check up to 7 tracers
        write(unit,*) 'fv_restart_end q    = ', mpp_chksum(Atm(n)%q(isc:iec,jsc:jec,:,iq))
      enddo

!---------------
! Check Min/Max:
!---------------
      call prt_maxmin('ZS', Atm(n)%phis, isc, iec, jsc, jec, Atm(n)%ng, 1, 1./grav, gid==masterproc)
      call prt_maxmin('PS', Atm(n)%ps, isc, iec, jsc, jec, Atm(n)%ng, 1, 0.01, gid==masterproc)
      call prt_maxmin('U ', Atm(n)%u(isd:ied,jsd:jed,1:npz), isc, iec, jsc, jec, Atm(n)%ng, npz, 1., gid==masterproc)
      call prt_maxmin('V ', Atm(n)%v(isd:ied,jsd:jed,1:npz), isc, iec, jsc, jec, Atm(n)%ng, npz, 1., gid==masterproc)
      if ( .not. Atm(n)%hydrostatic )    &
      call prt_maxmin('W ', Atm(n)%w , isc, iec, jsc, jec, Atm(n)%ng, npz, 1., gid==masterproc)
      call prt_maxmin('T ', Atm(n)%pt, isc, iec, jsc, jec, Atm(n)%ng, npz, 1., gid==masterproc)
! Write4 energy correction term
#endif
    end do

    call fv_io_write_restart(Atm)
    module_is_initialized = .FALSE.

#ifdef EFLUX_OUT
    if( gid==masterproc ) then
        write(*,*) steps, 'Mean equivalent Heat flux for this integration period=',efx_sum/real(max(1,steps)), &
                          'Mean mountain torque=',mtq_sum/real(max(1,steps))
        file_unit = get_unit()
        open (unit=file_unit, file='e_flux.data', form='unformatted',status='unknown', access='sequential')
        do n=1,steps
           write(file_unit) efx(n)
           write(file_unit) mtq(n)    ! time series global mountain torque
        enddo
        close(unit=file_unit)
    endif
#endif

  end subroutine fv_restart_end
  ! </SUBROUTINE> NAME="fv_restart_end"

end module fv_restart_mod
