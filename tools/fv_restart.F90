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

  use constants_mod,       only: kappa, omega, rdgas, grav
  use fv_arrays_mod,       only: fv_atmos_type
  use fv_io_mod,           only: fv_io_init, fv_io_read_restart, fv_io_write_restart
  use grid_tools,          only: area
  use grid_utils,          only: fC, f0, ptop, ptop_min, fill_ghost, big_number,   &
                                 make_eta_level
  use fv_diagnostics_mod,  only: prt_maxmin
  use init_hydro,          only: p_var
  use mpp_domains_mod,     only: mpp_update_domains, domain2d, DGRID_NE
  use mpp_mod,             only: mpp_chksum, stdout, mpp_error, FATAL
  use test_cases,          only: alpha, init_case
  use mp_mod,              only: gid, masterproc
#ifdef FV_LAND
  use surf_map,            only: sgh_g, oro_g
#endif


  implicit none
  private

  public :: fv_restart_init, fv_restart_end, fv_restart

  !--- private data type
  logical                       :: module_is_initialized = .FALSE.

  !--- version information variables ----
  character(len=128) :: version = '$Id: fv_restart.F90,v 14.0 2007/03/15 21:58:44 fms Exp $'
  character(len=128) :: tagname = '$Name: nalanda_2007_04 $'

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
  subroutine fv_restart(fv_domain, Atm, dt_atmos, seconds, days, cold_start)
    type(domain2d),      intent(inout) :: fv_domain
    type(fv_atmos_type), intent(inout) :: Atm(:)
    real,                intent(in)    :: dt_atmos
    integer,             intent(out)   :: seconds
    integer,             intent(out)   :: days
    logical,             intent(in)    :: cold_start


    integer :: i, j, k, n, ntileMe
    integer :: isc, iec, jsc, jec, npz, ncnst
    integer :: isd, ied, jsd, jed
    real rgrav

    rgrav = 1. / grav

    if(.not.module_is_initialized) call mpp_error(FATAL, 'You must call fv_restart_init.')

    ntileMe = size(Atm(:))

  ! This logic doesn't work very well.
  ! Shouldn't have read for all tiles then loop over tiles
    if(.not.cold_start) call fv_io_read_restart(fv_domain,Atm)

    seconds = 0; days = 0   ! Restart needs to be modified to record seconds and days.

    do n = 1, ntileMe

       isd = Atm(n)%isd
       ied = Atm(n)%ied
       jsd = Atm(n)%jsd
       jed = Atm(n)%jed
       npz = Atm(n)%npz
       ncnst = Atm(n)%ncnst
       isc = Atm(n)%isc; iec = Atm(n)%iec; jsc = Atm(n)%jsc; jec = Atm(n)%jec

      ! Init model data
      if(.not.cold_start)then  ! This is not efficient stacking if there are really more tiles than 1.

        call mpp_update_domains( Atm(n)%phis, fv_domain, complete=.true. )
! Fill big values in  the non-existinng corner regions:
        call fill_ghost(Atm(n)%phis, Atm(n)%npx, Atm(n)%npy, big_number)

#ifdef SW_DYNAMICS
        Atm(n)%pt(:,:,:)=1.
#else
        if (ptop/=Atm(n)%ak(1)) call mpp_error(FATAL,'FV restart: ptop not equal Atm(n)%ak(1)')
        call p_var(npz,         isc,         iec,       jsc,     jec,   ptop,     ptop_min,  &
                   Atm(n)%delp, Atm(n)%delz, Atm(n)%pt, Atm(n)%ps, Atm(n)%pe, Atm(n)%peln,   &
                   Atm(n)%pk,   Atm(n)%pkz, kappa, Atm(n)%q, Atm(n)%ng, ncnst,  Atm(n)%dry_mass,  &
                   Atm(n)%adjust_dry_mass,  Atm(n)%mountain, Atm(n)%full_phys,  Atm(n)%hydrostatic, &
                   Atm(n)%k_top)
#endif
        do j=jsd,jed+1
           do i=isd,ied+1
              fC(i,j) = 2.*Omega*( -cos(Atm(n)%grid(i,j,1))*cos(Atm(n)%grid(i,j,2))*sin(alpha) + &
                                    sin(Atm(n)%grid(i,j,2))*cos(alpha) )
           enddo
        enddo
        do j=jsd,jed
           do i=isd,ied
             f0(i,j) = 2.*Omega*( -cos(Atm(n)%agrid(i,j,1))*cos(Atm(n)%agrid(i,j,2))*sin(alpha) + &
                                    sin(Atm(n)%agrid(i,j,2))*cos(alpha) )
           enddo
        enddo
      else
! Setup Case to Run
          call init_case(Atm(n)%u,Atm(n)%v,Atm(n)%pt,Atm(n)%delp,Atm(n)%q,Atm(n)%phis, Atm(n)%ps,Atm(n)%pe, &
                         Atm(n)%peln,Atm(n)%pk,Atm(n)%pkz, Atm(n)%uc,Atm(n)%vc, Atm(n)%ua,Atm(n)%va,        & 
                         Atm(n)%ak, Atm(n)%bk, Atm(n)%npx, Atm(n)%npy, npz, Atm(n)%ng, ncnst, &
                         Atm(n)%k_top, Atm(n)%ndims, Atm(n)%ntiles, Atm(n)%dry_mass, Atm(n)%mountain,       &
#ifdef MAKE_HYBRID_Z
                         Atm(n)%full_phys,  .false.       , Atm(n)%delz, Atm(n)%ze0)
#else
                         Atm(n)%full_phys, Atm(n)%hybrid_z, Atm(n)%delz, Atm(n)%ze0)
#endif


#ifdef FV_LAND
        do j=jsc,jec
           do i=isc,iec
              Atm(n)%sgh(i,j) = sgh_g(i,j)
              Atm(n)%oro(i,j) = oro_g(i,j)
           enddo
        enddo
#endif
      endif

     if ( Atm(n)%hybrid_z ) then
#ifdef MAKE_HYBRID_Z
          do k=1,npz
             do j=jsc,jec
                do i=isc,iec
                   Atm(n)%delz(i,j,k) = (rdgas*rgrav)*Atm(n)%pt(i,j,k)*(Atm(n)%peln(i,k,j)-Atm(n)%peln(i,k+1,j))
                enddo
             enddo
          enddo

          do j=jsc,jec
             do i=isc,iec
                Atm(n)%ze0(i,j,npz+1) = Atm(n)%phis(i,j) * rgrav
             enddo
          enddo
          do k=npz,1,-1
             do j=jsc,jec
                do i=isc,iec
                   Atm(n)%ze0(i,j,k) = Atm(n)%ze0(i,j,k+1) - Atm(n)%delz(i,j,k)
                enddo
             enddo
          enddo
#endif
!         call prt_maxmin('ZE0', Atm(n)%ze0,  isc, iec, jsc, jec, 0, npz, 1.E-3, gid==masterproc)
!         call prt_maxmin('DZ0', Atm(n)%delz, isc, iec, jsc, jec, 0, npz, 1.   , gid==masterproc)
          call make_eta_level(npz, Atm(n)%pe, area, Atm(n)%ks, Atm(n)%ak, Atm(n)%bk)
      endif

      write(stdout(),*)
      write(stdout(),*) 'fv_restart u    = ', mpp_chksum(Atm(n)%u(isc:iec,jsc:jec,:))
      write(stdout(),*) 'fv_restart v    = ', mpp_chksum(Atm(n)%v(isc:iec,jsc:jec,:))
      write(stdout(),*) 'fv_restart delp = ', mpp_chksum(Atm(n)%delp(isc:iec,jsc:jec,:))
!     write(stdout(),*) 'fv_restart phis = ', mpp_chksum(Atm(n)%phis(isc:iec,jsc:jec))
#ifdef SW_DYNAMICS
      call prt_maxmin('H ', Atm(n)%delp, isc, iec, jsc, jec, Atm(n)%ng, 1, rgrav, gid==masterproc)
#else
      write(stdout(),*) 'fv_restart pt   = ', mpp_chksum(Atm(n)%pt(isc:iec,jsc:jec,:))
!     write(stdout(),*) 'fv_restart u_srf = ', mpp_chksum(Atm(n)%u_srf(isc:iec,jsc:jec))
!     write(stdout(),*) 'fv_restart v_srf = ', mpp_chksum(Atm(n)%v_srf(isc:iec,jsc:jec))
      if (ncnst>0) write(stdout(),*) 'fv_init q    = ', mpp_chksum(Atm(n)%q(isc:iec,jsc:jec,:,:))
!---------------
! Check Min/Max:
!---------------
      call prt_maxmin('ZS', Atm(n)%phis, isc, iec, jsc, jec, Atm(n)%ng, 1, rgrav, gid==masterproc)

      if ( (.not.Atm(n)%hydrostatic) .and. (.not.Atm(n)%Make_NH) ) then
            call prt_maxmin('DZ', Atm(n)%delz, isc, iec, jsc, jec, 0, npz, 1., gid==masterproc)
            if ( Atm(n)%hybrid_z ) then
            call prt_maxmin('ZTOP(km)', Atm(n)%ze0, isc, iec, jsc, jec, 0, 1, 1.E-3, gid==masterproc)
            call prt_maxmin('DZ_top', Atm(n)%delz, isc, iec, jsc, jec, 0, 1, 1.E-3, gid==masterproc)
            endif
      endif

      call prt_maxmin('PS', Atm(n)%ps, isc, iec, jsc, jec, Atm(n)%ng, 1,    0.01, gid==masterproc)
      call prt_maxmin('T ', Atm(n)%pt, isc, iec, jsc, jec, Atm(n)%ng, npz, 1., gid==masterproc)
      if (ncnst>0) call prt_maxmin('Q1', Atm(n)%q(isd,jsd,1,1), isc, iec, jsc, jec, Atm(n)%ng, npz, 1.,gid==masterproc)
      if (ncnst>1) call prt_maxmin('Q2', Atm(n)%q(isd,jsd,1,2), isc, iec, jsc, jec, Atm(n)%ng, npz, 1.,gid==masterproc)
      if (ncnst>2) call prt_maxmin('Q3', Atm(n)%q(isd,jsd,1,3), isc, iec, jsc, jec, Atm(n)%ng, npz, 1.,gid==masterproc)
      if (ncnst>3) call prt_maxmin('Q4', Atm(n)%q(isd,jsd,1,4), isc, iec, jsc, jec, Atm(n)%ng, npz, 1.,gid==masterproc)
#endif
      call prt_maxmin('U ', Atm(n)%u(isc:iec,jsc:jec,1:npz), isc, iec, jsc, jec, 0, npz, 1., gid==masterproc)
      call prt_maxmin('V ', Atm(n)%v(isc:iec,jsc:jec,1:npz), isc, iec, jsc, jec, 0, npz, 1., gid==masterproc)
      if ( .not.Atm(n)%hydrostatic )   &
      call prt_maxmin('W ', Atm(n)%w, isc, iec, jsc, jec, Atm(n)%ng, npz, 1.,gid==masterproc)

      if ( (.not.Atm(n)%hydrostatic) .and. Atm(n)%Make_NH ) then
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

      write(stdout(),*)
    end do

  end subroutine fv_restart
  ! </SUBROUTINE> NAME="fv_restart"


  !#####################################################################
  ! <SUBROUTINE NAME="fv_restart_end">
  !
  ! <DESCRIPTION>
  ! Initialize the fv core restart facilities
  ! </DESCRIPTION>
  !
  subroutine fv_restart_end(fv_domain,Atm)
    type(domain2d),      intent(in) :: fv_domain
    type(fv_atmos_type), intent(in) :: Atm(:)

    integer :: isc, iec, jsc, jec
    integer :: n, ntileMe
    integer :: isd, ied, jsd, jed, npz

    ntileMe = size(Atm(:))

    do n = 1, ntileMe
      isc = Atm(n)%isc; iec = Atm(n)%iec; jsc = Atm(n)%jsc; jec = Atm(n)%jec

       isd = Atm(n)%isd
       ied = Atm(n)%ied
       jsd = Atm(n)%jsd
       jed = Atm(n)%jed
       npz = Atm(n)%npz
 
      write(stdout(),*)
      write(stdout(),*) 'fv_restart_end u    = ', mpp_chksum(Atm(n)%u(isc:iec,jsc:jec,:))
      write(stdout(),*) 'fv_restart_end v    = ', mpp_chksum(Atm(n)%v(isc:iec,jsc:jec,:))
      if ( .not. Atm(n)%hydrostatic )    &
      write(stdout(),*) 'fv_restart_end w    = ', mpp_chksum(Atm(n)%w(isc:iec,jsc:jec,:))

      write(stdout(),*) 'fv_restart_end delp = ', mpp_chksum(Atm(n)%delp(isc:iec,jsc:jec,:))
!     write(stdout(),*) 'fv_restart_end phis = ', mpp_chksum(Atm(n)%phis(isc:iec,jsc:jec))
#ifndef SW_DYNAMICS
      write(stdout(),*) 'fv_restart_end pt   = ', mpp_chksum(Atm(n)%pt(isc:iec,jsc:jec,:))
!     write(stdout(),*) 'fv_restart_end u_srf = ', mpp_chksum(Atm(n)%u_srf(isc:iec,jsc:jec))
!     write(stdout(),*) 'fv_restart_end v_srf = ', mpp_chksum(Atm(n)%v_srf(isc:iec,jsc:jec))
      write(stdout(),*) 'fv_restart_end q    = ', mpp_chksum(Atm(n)%q(isc:iec,jsc:jec,:,:))

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
#endif
    end do

    call fv_io_write_restart(fv_domain,Atm)
    module_is_initialized = .FALSE.
  end subroutine fv_restart_end
  ! </SUBROUTINE> NAME="fv_restart_end"

end module fv_restart_mod
