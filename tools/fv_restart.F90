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

  use constants_mod,       only: kappa, pi, omega, rdgas, grav, rvgas, cp_air, radius
  use fv_arrays_mod,       only: fv_atmos_type
  use fv_io_mod,           only: fv_io_init, fv_io_read_restart, fv_io_write_restart, &
                                 remap_restart, fv_io_register_restart, fv_io_write_BCs, fv_io_read_BCs
  use fv_grid_tools_mod,   only: area, dx, dy, rdxa, rdya, dxc, dyc
  use fv_grid_utils_mod,   only: fc, f0, ptop, ptop_min, fill_ghost, big_number,   &
                                 make_eta_level, deglat, cubed_to_latlon, da_min, great_circle_dist
  use fv_diagnostics_mod,  only: prt_maxmin
  use init_hydro_mod,      only: p_var
  use mpp_domains_mod,     only: mpp_update_domains, domain2d, DGRID_NE
  use mpp_mod,             only: mpp_chksum, stdout, mpp_error, FATAL, NOTE, get_unit
  use test_cases_mod,      only: alpha, init_case, init_double_periodic, init_latlon
  use fv_mp_mod,           only: gid, masterproc, switch_current_Atm
  use fv_surf_map_mod,     only: sgh_g, oro_g
  use fv_diagnostics_mod,  only: steps, efx, efx_sum, mtq, mtq_sum
  use tracer_manager_mod,  only: get_tracer_names
  use field_manager_mod,   only: MODEL_ATMOS
  use external_ic_mod,     only: get_external_ic, get_cubed_sphere_terrain
  use fv_eta_mod,          only: compute_dz_var, compute_dz_L32, set_hybrid_z
  use fv_surf_map_mod,     only: del2_cubed_sphere, del4_cubed_sphere
  use boundary_mod,        only: gather_grid, fill_nested_grid, nested_grid_BC, update_coarse_grid, nested_grid_BC_save
   use tracer_manager_mod, only: get_tracer_index
   use field_manager_mod,  only: MODEL_ATMOS
   use fv_timing_mod,      only: timing_on, timing_off
  use mpp_domains_mod,     only: mpp_get_compute_domain, mpp_get_data_domain, mpp_get_global_domain
  use fv_current_grid_mod, only: efx_nest, efx_sum_nest, nest_domain, npes_per_tile, test_case, first_step
  use boundary_mod,        only: nested_grid_BC_apply
  use fv_mp_mod,           only: grids_on_this_pe, concurrent
  use mpp_mod,             only: mpp_send, mpp_recv, mpp_sync_self, mpp_set_current_pelist, mpp_get_current_pelist, mpp_npes
  use mpp_domains_mod,     only: CENTER, CORNER, NORTH, EAST,  mpp_get_C2F_index, WEST, SOUTH
  use mpp_domains_mod,     only: mpp_global_field
  use fms_mod,             only: file_exist

  implicit none
  private

  public :: fv_restart_init, fv_restart_end, fv_restart, fv_write_restart, setup_nested_boundary_halo
  public :: d2c_setup, d2a_setup, setup_nested_BC_ucvc

  !--- private data type
  logical                       :: module_is_initialized = .FALSE.

!---- version number -----
  character(len=128) :: version = '$Id: fv_restart.F90,v 17.0.2.2.2.4.2.20.2.1.4.1 2012/09/28 16:05:26 Rusty.Benson Exp $'
  character(len=128) :: tagname = '$Name: siena_201303 $'

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
    logical,             intent(inout)    :: cold_start
    integer,             intent(in)    :: grid_type


    integer :: i, j, k, n, ntileMe, nt
    integer :: isc, iec, jsc, jec, npz, npz_rst, ncnst
    integer :: isd, ied, jsd, jed
    integer isd_p, ied_p, jsd_p, jed_p, isc_p, iec_p, jsc_p, jec_p, isg, ieg, jsg,jeg, npx_p, npy_p
    real, allocatable :: g_dat(:,:,:)

    integer :: unit
    real, allocatable :: dz1(:)
    real rgrav, f00, ztop
    logical :: hybrid
    character(len=128):: tname, errstring, fname

    rgrav = 1. / grav

    if(.not.module_is_initialized) call mpp_error(FATAL, 'You must call fv_restart_init.')

    ntileMe = size(Atm(:))

    do n = 1, ntileMe

!       if (Atm(n)%nested) then
!          write(fname,'(A, I1, A)') 'INPUT/fv_core_g', Atm(n)%grid_number, '.res.nc'
!          cold_start = .not. file_exist(fname)
!          Atm(n)%warm_start = .not. cold_start
!       endif

       if (.not. grids_on_this_pe(n)) then
!!$          write(errstring,'(A, I4)') 'NOT RESTARTING GRID ', n
!!$          call mpp_error(NOTE, errstring)
          
          !Even if this grid is not on this PE, if it has child grids we must send
          !along the data that is needed. 

          !This block of code is only called for concurrent runs.

          if (Atm(n)%nested) then
             if (Atm(n)%external_ic) then
                !Need to set up topography halo
                call nested_grid_BC(Atm(n)%phis, Atm(n)%parent_grid%phis, Atm(n)%nest_domain, &
                     Atm(n)%ind_h, Atm(n)%wt_h, 0, 0, &
                     Atm(n)%npx, Atm(n)%npy, isg, ieg, jsg, jeg, proc_in=.false.)
             endif
             if (cold_start) then
                call fill_nested_grid_data(Atm(n:n), .false.)
                if (.not. concurrent) call setup_nested_BC_ucvc(Atm(n), npz, .false.)
             else
                if (.not. concurrent) call setup_nested_BC_ucvc(Atm(n), npz, .false., saveBCs_in=.false.)
             end if
          endif

          cycle

       endif
       call switch_current_Atm(Atm(n))



    npz     = Atm(1)%npz
    npz_rst = Atm(1)%npz_rst

    !--- call fv_io_register_restart to register restart field to be written out in fv_io_write_restart
    call fv_io_register_restart(Atm(n)%domain,Atm(n:n))
    if( .not.cold_start .and. (.not. Atm(n)%external_ic) ) then

       !This is needed for the idealized test cases
       if (Atm(n)%nested) call setup_nested_grid_topography(Atm(n))

        if ( npz_rst /= 0 .and. npz_rst /= npz ) then
!            Remap vertically the prognostic variables for the chosen vertical resolution
             if( gid==masterproc ) then
                 write(*,*) ' '
                 write(*,*) '***** Important Note from FV core ********************'
                 write(*,*) 'Remapping dynamic IC from', npz_rst, 'levels to ', npz,'levels'
                 write(*,*) '***** End Note from FV core **************************'
                 write(*,*) ' '
             endif
             call remap_restart( Atm(n)%domain, Atm(n:n) )
             if( gid==masterproc ) write(*,*) 'Done remapping dynamical IC'
        else
             call fv_io_read_restart(Atm(n)%domain,Atm(n:n))
        endif
    endif

!---------------------------------------------------------------------------------------------
! Read, interpolate (latlon to cubed), then remap vertically with terrain adjustment if needed
!---------------------------------------------------------------------------------------------
    if ( Atm(n)%external_ic ) then
         call get_external_ic(Atm(n:n), Atm(n)%domain) 
         if( gid==masterproc ) write(*,*) 'IC generated from the specified external source'
    endif

    seconds = 0; days = 0   ! Restart needs to be modified to record seconds and days.

! Notes by Jeff D.
  ! This logic doesn't work very well.
  ! Shouldn't have read for all tiles then loop over tiles

       isd = Atm(n)%isd
       ied = Atm(n)%ied
       jsd = Atm(n)%jsd
       jed = Atm(n)%jed
       ncnst = Atm(n)%ncnst
       isc = Atm(n)%isc; iec = Atm(n)%iec; jsc = Atm(n)%jsc; jec = Atm(n)%jec

      ! Init model data
      if(.not.cold_start)then  ! This is not efficient stacking if there are really more tiles than 1.
         first_step = .false.
        if ( Atm(n)%mountain ) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!! Additional terrain filter -- should not be called repeatedly !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if ( Atm(n)%n_zs_filter > 0 ) then
              if ( Atm(n)%nord_zs_filter == 2 ) then
                   call del2_cubed_sphere(Atm(n)%npx, Atm(n)%npy, Atm(n)%phis, area, dx, dy,   &
                                          dxc, dyc, Atm(n)%n_zs_filter, 0.20*da_min, .false., oro_g)
                   if ( gid==masterproc ) write(*,*) 'Warning !!! del-2 terrain filter has been applied ', Atm(n)%n_zs_filter, ' times'
              else if( Atm(n)%nord_zs_filter == 4 ) then
                   call del4_cubed_sphere(Atm(n)%npx, Atm(n)%npy, Atm(n)%phis, area, dx, dy,   &
                                        dxc, dyc, Atm(n)%n_zs_filter, .false., oro_g)
                 if ( gid==masterproc ) write(*,*) 'Warning !!! del-4 terrain filter has been applied ', Atm(n)%n_zs_filter, ' times'
              endif
            endif

            call mpp_update_domains( Atm(n)%phis, Atm(n)%domain, complete=.true. )
        else
             Atm(n)%phis = 0.
            if( gid==masterproc ) write(*,*) 'phis set to zero'
         endif !mountain

         if (Atm(n)%nested) then
            if ( npz_rst /= 0 .and. npz_rst /= npz ) then
               call setup_nested_boundary_halo(Atm(n)) 
            else
               call fv_io_read_BCs(Atm(n))
               !Following line to make sure u and v are consistent across processor subdomains
               call mpp_update_domains(Atm(n)%u, Atm(n)%v, Atm(n)%domain, gridtype=DGRID_NE, complete=.true.)
            endif
         endif

#ifdef SW_DYNAMICS
        Atm(n)%pt(:,:,:)=1.
#else
        if ( .not.Atm(n)%hybrid_z ) then
           if(ptop/=Atm(n)%ak(1)) call mpp_error(FATAL,'FV restart: ptop not equal Atm(n)%ak(1)')
        else
           ptop = Atm(n)%ak(1);  Atm(n)%ks = 0
        endif
        call p_var(npz,         isc,         iec,       jsc,     jec,   ptop,     ptop_min,  &
                   Atm(n)%delp, Atm(n)%delz, Atm(n)%pt, Atm(n)%ps, Atm(n)%pe, Atm(n)%peln,   &
                   Atm(n)%pk,   Atm(n)%pkz, kappa, Atm(n)%q, Atm(n)%ng, ncnst,  Atm(n)%dry_mass,  &
                   Atm(n)%adjust_dry_mass,  Atm(n)%mountain, Atm(n)%moist_phys,  Atm(n)%hydrostatic, &
                   Atm(n)%nwat, Atm(n)%make_nh)
#endif
        if ( grid_type < 7 .and. grid_type /= 4 ) then
! Fill big values in the non-existing corner regions:
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
                           Atm(n)%ndims, Atm(n)%ntiles, Atm(n)%dry_mass, Atm(n)%mountain,       &
                           Atm(n)%moist_phys, Atm(n)%hydrostatic, hybrid, Atm(n)%delz, Atm(n)%ze0)
            endif
         elseif (grid_type == 4) then
            call init_double_periodic(Atm(n)%u,Atm(n)%v,Atm(n)%w,Atm(n)%pt,Atm(n)%delp,Atm(n)%q,Atm(n)%phis, Atm(n)%ps,Atm(n)%pe, &
                                      Atm(n)%peln,Atm(n)%pk,Atm(n)%pkz, Atm(n)%uc,Atm(n)%vc, Atm(n)%ua,Atm(n)%va,        & 
                                      Atm(n)%ak, Atm(n)%bk, Atm(n)%npx, Atm(n)%npy, npz, Atm(n)%ng, ncnst, Atm(n)%nwat,  &
                                      Atm(n)%ndims, Atm(n)%ntiles, Atm(n)%dry_mass, Atm(n)%mountain,       &
                                      Atm(n)%moist_phys, Atm(n)%hydrostatic, hybrid, Atm(n)%delz, Atm(n)%ze0)
            if( gid==masterproc ) write(*,*) 'Doubly Periodic IC generated'
         elseif (grid_type == 5 .or. grid_type == 6) then
            call init_latlon(Atm(n)%u,Atm(n)%v,Atm(n)%pt,Atm(n)%delp,Atm(n)%q,Atm(n)%phis, Atm(n)%ps,Atm(n)%pe, &
                             Atm(n)%peln,Atm(n)%pk,Atm(n)%pkz, Atm(n)%uc,Atm(n)%vc, Atm(n)%ua,Atm(n)%va,        &
                             Atm(n)%ak, Atm(n)%bk, Atm(n)%npx, Atm(n)%npy, npz, Atm(n)%ng, ncnst, &
                             Atm(n)%ndims, Atm(n)%ntiles, Atm(n)%dry_mass, Atm(n)%mountain,       &
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


        !Set up nested grids
        If (Atm(n)%nested) call fill_nested_grid_data(Atm(n:n))

     endif  !end cold_start check

  end do


    do n = ntileMe,1,-1
       if (.not. grids_on_this_pe(n)) then
          if (Atm(n)%nested .and. cold_start) call fill_nested_grid_data_end(Atm(n), .false.)
          cycle
       endif
       if (Atm(n)%nested .and. cold_start) call fill_nested_grid_data_end(Atm(n))
    end do

    do n = 1, ntileMe
       if (.not. grids_on_this_pe(n)) cycle
       call switch_current_Atm(Atm(n))

       isd = Atm(n)%isd
       ied = Atm(n)%ied
       jsd = Atm(n)%jsd
       jed = Atm(n)%jed
       ncnst = Atm(n)%ncnst
       isc = Atm(n)%isc; iec = Atm(n)%iec; jsc = Atm(n)%jsc; jec = Atm(n)%jec

       if (.not. concurrent) then
          !If this grid has any grids nested inside of it we need to initialize uc, vc
          if (ANY(Atm(n)%child_grids)) then

             call mpp_update_domains(Atm(n)%u, Atm(n)%v, &
                  Atm(n)%domain, gridtype=DGRID_NE, complete=.true.)
             call mpp_sync_self
             do k=1,npz
                call d2c_setup(Atm(n)%u(isd,jsd,k),  Atm(n)%v(isd,jsd,k), &
                     Atm(n)%uc(isd,jsd,k), Atm(n)%vc(isd,jsd,k), &
                     Atm(n)%nord>0, &
                     isd,ied,jsd,jed, isc,iec,jsc,jec, &
                     Atm(n)%npx,Atm(n)%npy, &
                     Atm(n)%grid_type, Atm(n)%nested, &
                     Atm(n)%se_corner, Atm(n)%sw_corner, &
                     Atm(n)%ne_corner, Atm(n)%nw_corner, &
                     Atm(n)%rsin_u, Atm(n)%rsin_v, &
                     Atm(n)%cosa_s, Atm(n)%rsin2)
             end do
          endif


          if (Atm(n)%nested) then
             if (cold_start) then
                !This routine sets up the uc and vc time t1 BCs
                call setup_nested_BC_ucvc(Atm(n),npz) 
             else
                !This just sets up coarse-grid uc and vc so that
                !they are available for nested-grid BCs to be created
                if (concurrent) call setup_nested_BC_ucvc(Atm(n),npz, .false., .false.)
             endif
          endif
       endif

!---------------------------------------------------------------------------------------------
! Transform the (starting) Eulerian vertical coordinate from sigma-p to hybrid_z
     if ( Atm(n)%hybrid_z ) then
       if ( Atm(n)%make_hybrid_z ) then
          allocate ( dz1(npz) )
          if( npz==32 ) then
              call compute_dz_L32(npz, ztop, dz1)
          else
              ztop = 45.E3
              call compute_dz_var(npz, ztop, dz1)
          endif
          call set_hybrid_z(isc, iec, jsc, jec, Atm(n)%ng, npz, ztop, dz1, rgrav,  &
                            Atm(n)%phis, Atm(n)%ze0)
          deallocate ( dz1 )
!         call prt_maxmin('ZE0', Atm(n)%ze0,  isc, iec, jsc, jec, 0, npz, 1.E-3, gid==masterproc)
!         call prt_maxmin('DZ0', Atm(n)%delz, isc, iec, jsc, jec, 0, npz, 1.   , gid==masterproc)
       endif
!      call make_eta_level(npz, Atm(n)%pe, area, Atm(n)%ks, Atm(n)%ak, Atm(n)%bk)
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

  subroutine setup_nested_boundary_halo(Atm, proc_in)

    type(fv_atmos_type), intent(INOUT) :: Atm
    logical, INTENT(IN), OPTIONAL :: proc_in
    real, allocatable :: g_dat(:,:,:), g_dat2(:,:,:)
    real, allocatable :: pt_coarse(:,:,:)
    integer i,j,k,nq, sphum, ncnst, istart, iend, npz
    integer isc, iec, jsc, jec, isd, ied, jsd, jed
    integer isd_p, ied_p, jsd_p, jed_p, isc_p, iec_p, jsc_p, jec_p, isg, ieg, jsg,jeg, npx_p, npy_p
    real zvir
    logical process

    if (PRESENT(proc_in)) then
       process = proc_in
    else
       process = .true.
    endif

    isd = Atm%isd
    ied = Atm%ied
    jsd = Atm%jsd
    jed = Atm%jed
    ncnst = Atm%ncnst
    isc = Atm%isc; iec = Atm%iec; jsc = Atm%jsc; jec = Atm%jec
    npz     = Atm%npz    

    call mpp_get_data_domain( Atm%parent_grid%domain, &
         isd_p,  ied_p,  jsd_p,  jed_p  )
    call mpp_get_compute_domain( Atm%parent_grid%domain, &
         isc_p,  iec_p,  jsc_p,  jec_p  )
    call mpp_get_global_domain( Atm%parent_grid%domain, &
         isg, ieg, jsg, jeg, xsize=npx_p, ysize=npy_p)

    call nested_grid_BC(Atm%delp, Atm%parent_grid%delp, Atm%nest_domain, &
         Atm%ind_h, Atm%wt_h, 0, 0, &
         Atm%npx, Atm%npy, npz, isg, ieg, jsg, jeg, proc_in=process)
    do nq=1,ncnst
       call nested_grid_BC(Atm%q(:,:,:,nq), &
            Atm%parent_grid%q(:,:,:,nq), Atm%nest_domain, &
            Atm%ind_h, Atm%wt_h, 0, 0, &
            Atm%npx, Atm%npy, npz, isg, ieg, jsg, jeg, proc_in=process)
    end do


    !Filling phis?
    !In idealized test cases, where the topography is EXACTLY known (ex case 13),
    !interpolating the topography yields a much worse result. In comparison in
    !real topography cases little difference is seen.

    !This is probably because the halo phis, which is used to compute
    !geopotential height (gz, gh), only affects the interior by being
    !used to compute corner gz in a2b_ord[24]. We might suppose this
    !computation would be more accurate when using values of phis which
    !are more consistent with those on the interior (ie the exactly-known
    !values) than the crude values given through linear interpolation.

    !HOWEVER: Since restarts do not save data in the halo, we INTERPOLATE if test_case == 11 or -999 (not specified)


#ifndef SW_DYNAMICS
    !pt --- actually temperature

    call nested_grid_BC(Atm%pt, Atm%parent_grid%pt, Atm%nest_domain, &
         Atm%ind_h, Atm%wt_h, 0, 0, &
         Atm%npx, Atm%npy, npz, isg, ieg, jsg, jeg, proc_in=process)    

    !Need to fill boundaries with INITIAL coarse-grid potential temperature,
    !because the nested grid cannot do it by itself (values of pkz in the halo are not saved)
    !We will want to INTERPOLATE pkz from the coarse grid and then use it to compute
    !theta on the haloes of the nested grid. Computing theta on the coarse grid and then
    !interpolating yields a different (and inconsistent) answer with the way it is computed
    !in the interior

    if ( Atm%nwat > 0 ) then
       sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
    else
       sphum = 1
    endif
    if ( Atm%parent_grid%adiabatic .or. Atm%parent_grid%do_Held_Suarez ) then
       zvir = 0.         ! no virtual effect
    else
       zvir = rvgas/rdgas - 1.
    endif

    allocate(pt_coarse(isd:ied,jsd:jed,npz))

!! FIXME: nested_grid_BC does not work very well with pkz, since pkz does not have
    !! a halo; hence the having to allocate a new variable and such.

    pt_coarse = 0.
    allocate(g_dat(isd_p:ied_p,  jsd_p:jed_p  , npz))
    g_dat = 0.
    if (ANY(gid==Atm%parent_grid%pelist)) g_dat(isc_p:iec_p,  jsc_p:jec_p  , :) = Atm%parent_grid%pkz
!    call nested_grid_BC(pt_coarse, Atm%parent_grid%pkz, Atm%nest_domain, &
    call nested_grid_BC(pt_coarse, g_dat, Atm%nest_domain, &
         Atm%ind_h, Atm%wt_h, 0, 0, &
         Atm%npx, Atm%npy, npz, isg, ieg, jsg, jeg, proc_in=process)
    deallocate(g_dat)
    

    if (process) then

    if (Atm%is == 1) then
       do k=1,npz
          do j=Atm%jsd,Atm%jed
             do i=Atm%isd,0
                Atm%pt(i,j,k) = cp_air*Atm%pt(i,j,k)/pt_coarse(i,j,k)*(1.+zvir*Atm%q(i,j,k,sphum))
             end do
          end do
       end do
    end if

    if (Atm%js == 1) then
       if (Atm%is == 1) then
          istart = Atm%is
       else
          istart = Atm%isd
       end if
       if (Atm%ie == Atm%npx-1) then
          iend = Atm%ie
       else
          iend = Atm%ied
       end if

       do k=1,npz
          do j=Atm%jsd,0
             do i=istart,iend
                Atm%pt(i,j,k) = cp_air*Atm%pt(i,j,k)/pt_coarse(i,j,k)*(1.+zvir*Atm%q(i,j,k,sphum))
             end do
          end do
       end do
    end if

    if (Atm%ie == Atm%npx-1) then
       do k=1,npz
          do j=Atm%jsd,Atm%jed
             do i=Atm%npx,Atm%ied
                Atm%pt(i,j,k) = cp_air*Atm%pt(i,j,k)/pt_coarse(i,j,k)*(1.+zvir*Atm%q(i,j,k,sphum))
             end do
          end do
       end do
    end if

    if (Atm%je == Atm%npy-1) then
       if (Atm%is == 1) then
          istart = Atm%is
       else
          istart = Atm%isd
       end if
       if (Atm%ie == Atm%npx-1) then
          iend = Atm%ie
       else
          iend = Atm%ied
       end if

       do k=1,npz
          do j=Atm%npy,Atm%jed
             do i=istart,iend
                Atm%pt(i,j,k) = cp_air*Atm%pt(i,j,k)/pt_coarse(i,j,k)*(1.+zvir*Atm%q(i,j,k,sphum))
             end do
          end do
       end do
    end if

    end if !process

    deallocate(pt_coarse)

    if (.not. Atm%hydrostatic) then


       !w
       call nested_grid_BC(Atm%w(:,:,:), &
            Atm%parent_grid%w(:,:,:), &
            Atm%nest_domain, Atm%ind_h, Atm%wt_h, 0, 0, &
            Atm%npx, Atm%npy, npz, isg, ieg, jsg, jeg, proc_in=process)

    end if


#endif

    if (ANY(gid == Atm%pelist)) then
       call nested_grid_BC(Atm%u, Atm%parent_grid%u(:,:,:), &
            Atm%nest_domain, Atm%ind_u, Atm%wt_u, 0, 1, &
            Atm%npx, Atm%npy, npz, isg, ieg, jsg, jeg, proc_in=process)
       call nested_grid_BC(Atm%v, Atm%parent_grid%v(:,:,:), &
            Atm%nest_domain, Atm%ind_v, Atm%wt_v, 1, 0, &
            Atm%npx, Atm%npy, npz, isg, ieg, jsg, jeg, proc_in=process)
    else
       call nested_grid_BC(Atm%parent_grid%u(:,:,:), &
            Atm%nest_domain, 0, 1, &
            isg, ieg, jsg, jeg, npz)
       call nested_grid_BC(Atm%parent_grid%v(:,:,:), &
            Atm%nest_domain, 1, 0, &
            isg, ieg, jsg, jeg, npz)
    endif


    if (process) then
!!$#ifdef SW_DYNAMICS
!!$    !ps: first level only
!!$    !This is only valid for shallow-water simulations
!!$    do j=jsd,jed
!!$       do i=isd,ied
!!$
!!$          Atm%ps(i,j) = Atm%delp(i,j,1)/grav
!!$
!!$       end do
!!$    end do
!!$#endif
         call mpp_update_domains(Atm%u, Atm%v, Atm%domain, gridtype=DGRID_NE, complete=.true.)
      endif

  end subroutine setup_nested_boundary_halo


  subroutine setup_nested_BC_ucvc(Atm, npz, setBCs, saveBCs_in)

    use boundary_mod, only: gather_grid, nested_grid_BC_save, nested_grid_BC
   use fv_current_grid_mod, only: uc_east_t1, uc_west_t1, uc_south_t1, uc_north_t1
   use fv_current_grid_mod, only: vc_east_t1, vc_west_t1, vc_south_t1, vc_north_t1
   use fv_current_grid_mod, only: uc_east_t0, uc_west_t0, uc_south_t0, uc_north_t0
   use fv_current_grid_mod, only: vc_east_t0, vc_west_t0, vc_south_t0, vc_north_t0
   
    type(fv_atmos_type), intent(INOUT) :: Atm
    integer, intent(IN) :: npz
    logical, intent(IN), OPTIONAL :: setBCs, saveBCs_in

    real, allocatable :: g_dat(:,:,:)
    integer i,j,k,nq, sphum, ncnst, istart, iend
    integer isc, iec, jsc, jec, isd, ied, jsd, jed
    integer isd_p, ied_p, jsd_p, jed_p, isc_p, iec_p, jsc_p, jec_p, isg, ieg, jsg,jeg, npx_p, npy_p
    integer isg_n, ieg_n, jsg_n, jeg_n, npx_n, npy_n
    logical :: saveBCs

    saveBCs = .true.
    if (present(saveBCs_in)) then
       saveBCs = saveBCs_in
       if (.not. saveBCs .and. present(setBCs)) then
          if (setBCs) call mpp_error(FATAL, "setup_nested_BC_ucvc: cannot SET uc/vc BCs if they are not computed.")
       endif
    endif

      !The 3D test cases do not set up uc and vc, and rely on dyn_core to initialize them. However, 
      !    uc and vc are needed when doing grid nesting, because BCs for uc and vc are taken from
      !    a parent grid before dyn_core can be called. We call d2c here to set up a consistent
      !    uc and vc. 
            call mpp_get_data_domain( Atm%parent_grid%domain, &
                 isd_p,  ied_p,  jsd_p,  jed_p  )
            call mpp_get_compute_domain( Atm%parent_grid%domain, &
                 isc_p,  iec_p,  jsc_p,  jec_p  )
            call mpp_get_global_domain( Atm%parent_grid%domain, &
                 isg, ieg, jsg, jeg, xsize=npx_p, ysize=npy_p)
            call mpp_get_global_domain( Atm%domain, &
                 isg_n, ieg_n, jsg_n, jeg_n, xsize=npx_n, ysize=npy_n)


            !!! NOTE: we now do this elsewhere
         !Should only need to do this for a PARENT grid.
!         if (Atm%test_case > 9 .or. (.not. cold_start)) then
!         end if

         if (saveBCs) then
         call nested_grid_BC_save(Atm%parent_grid%vc, Atm%nest_domain, &
              Atm%ind_u, Atm%wt_u, 0, 1, Atm%npx, Atm%npy, npz, isg, ieg, jsg, jeg, &
              var_east=Atm%vc_east_t1, &
              var_west=Atm%vc_west_t1, &
              var_north=Atm%vc_north_t1, &
              var_south=Atm%vc_south_t1, ns=Atm%nsponge, proc_in=ANY(Atm%pelist == gid)  )

         call nested_grid_BC_save(Atm%parent_grid%uc, Atm%nest_domain, &
              Atm%ind_v, Atm%wt_v, 1, 0, Atm%npx, Atm%npy, npz, isg, ieg, jsg, jeg, &
              var_east=Atm%uc_east_t1, &
              var_west=Atm%uc_west_t1, &
              var_north=Atm%uc_north_t1, &
              var_south=Atm%uc_south_t1, ns=Atm%nsponge, proc_in=ANY(Atm%pelist == gid)  )

         if (present(setBCs)) then
            if (setBCs) then
               call nested_grid_BC_apply(Atm%uc, &
                 1, 0, Atm%npx, Atm%npy, npz, 1, 1, &
                 var_east=Atm%uc_east_t1, &
                 var_west=Atm%uc_west_t1, &
                 var_north=Atm%uc_north_t1, &
                 var_south=Atm%uc_south_t1, bctype=1, nsponge=0, s_weight=0.)
               call nested_grid_BC_apply(Atm%vc, &
                 0, 1, Atm%npx, Atm%npy, npz, 1, 1, &
                 var_east=Atm%vc_east_t1, &
                 var_west=Atm%vc_west_t1, &
                 var_north=Atm%vc_north_t1, &
                 var_south=Atm%vc_south_t1, bctype=1, nsponge=0, s_weight=0.)
            end if
         end if
         endif

  end subroutine setup_nested_BC_ucvc

  subroutine fill_nested_grid_data(Atm, proc_in)

    type(fv_atmos_type), intent(INOUT) :: Atm(:) !Only intended to be one element; needed for cubed_sphere_terrain
    logical, intent(IN), OPTIONAL :: proc_in
    real, allocatable :: g_dat(:,:,:), pt_coarse(:,:,:)
    integer :: i,j,k,nq, sphum, ncnst, istart, iend, npz
    integer :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer :: isd_p, ied_p, jsd_p, jed_p, isc_p, iec_p, jsc_p, jec_p
    integer :: isg, ieg, jsg,jeg, npx_p, npy_p
    integer :: isg_n, ieg_n, jsg_n, jeg_n, npx_n, npy_n
    real zvir, gh0, p1(2), p2(2), r, r0

    integer :: p, sending_proc
    logical process


    if (present(proc_in)) then
       process = proc_in
    else
       process = .true.
    endif

    isd = Atm(1)%isd
    ied = Atm(1)%ied
    jsd = Atm(1)%jsd
    jed = Atm(1)%jed
    ncnst = Atm(1)%ncnst
    isc = Atm(1)%isc; iec = Atm(1)%iec; jsc = Atm(1)%jsc; jec = Atm(1)%jec
    npz     = Atm(1)%npz    
    
    if (concurrent) then
       sending_proc = Atm(1)%parent_grid%pelist(1) + (Atm(1)%parent_tile-1)*Atm(1)%parent_grid%npes_per_tile
    else
       sending_proc = 0
    endif
       call mpp_get_data_domain( Atm(1)%parent_grid%domain, &
            isd_p,  ied_p,  jsd_p,  jed_p  )
       call mpp_get_compute_domain( Atm(1)%parent_grid%domain, &
            isc_p,  iec_p,  jsc_p,  jec_p  )
    call mpp_get_global_domain( Atm(1)%parent_grid%domain, &
         isg, ieg, jsg, jeg, xsize=npx_p, ysize=npy_p)


    if (process) then 
       
       ! * Initialize coriolis param:
       
       do j=jsd,jed+1
          do i=isd,ied+1
             Atm(1)%fc(i,j) = 2.*omega*( -1.*cos(Atm(1)%grid(i,j,1))*cos(Atm(1)%grid(i,j,2))*sin(alpha) + &
                  sin(Atm(1)%grid(i,j,2))*cos(alpha) )
          enddo
       enddo

       do j=jsd,jed
          do i=isd,ied
             Atm(1)%f0(i,j) = 2.*omega*( -1.*cos(Atm(1)%agrid(i,j,1))*cos(Atm(1)%agrid(i,j,2))*sin(alpha) + &
                  sin(Atm(1)%agrid(i,j,2))*cos(alpha) )
          enddo
       enddo

       call mpp_update_domains( Atm(1)%f0, Atm(1)%domain )

    endif
    if (test_case == 11 .or. test_case == -999) then
       if (gid == masterproc) print*, '  FILLING NESTED GRID HALO WITH INTERPOLATED TERRAIN'
       allocate(g_dat( isg:ieg, jsg:jeg, 1) )
       call timing_on('COMM_TOTAL')
       if (concurrent) then
          !Call gather grid on the procs that have the required data.
          !Then broadcast from the head PE to the receiving PEs
          if (ANY(Atm(1)%parent_grid%pelist == gid) .and. Atm(1)%parent_tile == Atm(1)%parent_grid%tile) then
             call mpp_global_field( &
                  Atm(1)%parent_grid%domain, &
                  Atm(1)%parent_grid%phis(isd_p:ied_p,jsd_p:jed_p), g_dat(:,:,1), position=CENTER)
             if (gid == sending_proc) then 
                do p=1,size(Atm(1)%pelist)
                   call mpp_send(g_dat,size(g_dat),Atm(1)%pelist(p))
                enddo
             endif
          endif
          if (ANY(Atm(1)%pelist == gid)) then
             call mpp_recv(g_dat, size(g_dat), sending_proc)
          endif
       else
          call mpp_error(FATAL, 'TOPOGRAPHY FOR SERIAL NESTING NOT IMPLEMENTED')
!!$          call gather_grid(Atm(1)%parent_grid, Atm(1)%parent_grid%phis(isd_p:ied_p,jsd_p:jed_p), g_dat(:,:,1), &
!!$               isd_p, ied_p, jsd_p, jed_p, isg, ieg, jsg, jeg, 0, 0, Atm(1)%parent_tile, &
!!$               concurrent, Atm(1)%pelist, size(Atm(1)%pelist))
       endif
       call timing_off('COMM_TOTAL')
       if (process) call fill_nested_grid(Atm(1)%phis, g_dat(:,:,1), &
            Atm(1)%ind_h, Atm(1)%wt_h, &
            0, 0,  isg, ieg, jsg, jeg)

       call mpp_sync_self
       deallocate(g_dat)
    endif

    !delp

    allocate(g_dat( isg:ieg, jsg:jeg, npz) )

    call timing_on('COMM_TOTAL')
    if (concurrent) then
       !Call gather grid on the procs that have the required data.
       !Then broadcast from the head PE to the receiving PEs
       if (ANY(Atm(1)%parent_grid%pelist == gid) .and. Atm(1)%parent_tile == Atm(1)%parent_grid%tile) then
          call mpp_global_field( &
               Atm(1)%parent_grid%domain, &
               Atm(1)%parent_grid%delp(isd_p:ied_p,jsd_p:jed_p,:), g_dat, position=CENTER)
          if (gid == sending_proc) then !crazy logic but what we have for now
             do p=1,size(Atm(1)%pelist)
                call mpp_send(g_dat,size(g_dat),Atm(1)%pelist(p))
             enddo
          endif
       endif
       if (ANY(Atm(1)%pelist == gid)) then
          call mpp_recv(g_dat, size(g_dat), sending_proc)
       endif
!!$!       call gather_grid(Atm(1)%parent_grid, Atm(1)%parent_grid%delp(isd_p:ied_p,jsd_p:jed_p,:), g_dat, &
!!$       call gather_grid(Atm(1)%parent_grid, Atm(1)%parent_grid%delp, g_dat, &
!!$            isd_p, ied_p, jsd_p, jed_p, isg, ieg, jsg, jeg, npz, 0, 0, Atm(1)%parent_tile, &
!!$            concurrent, Atm(1)%pelist, size(Atm(1)%pelist))
    else
       call gather_grid(Atm(1)%parent_grid, Atm(1)%parent_grid%delp(isd_p:ied_p,jsd_p:jed_p,:), g_dat, &
            isd_p, ied_p, jsd_p, jed_p, isg, ieg, jsg, jeg, npz, 0, 0, Atm(1)%parent_tile, &
            concurrent, Atm(1)%pelist, size(Atm(1)%pelist))
    endif
    call timing_off('COMM_TOTAL')
    if (process) call fill_nested_grid(Atm(1)%delp, g_dat, &
         Atm(1)%ind_h, Atm(1)%wt_h, &
         0, 0,  isg, ieg, jsg, jeg, npz)

    call mpp_sync_self

    !tracers
    do nq=1,ncnst

       call timing_on('COMM_TOTAL')
       if (concurrent) then
          if (ANY(Atm(1)%parent_grid%pelist == gid) .and. Atm(1)%parent_tile == Atm(1)%parent_grid%tile) then
             call mpp_global_field( &
               Atm(1)%parent_grid%domain, &
               Atm(1)%parent_grid%q(isd_p:ied_p,jsd_p:jed_p,:,nq), g_dat, position=CENTER)
             if (gid == sending_proc) then
                do p=1,size(Atm(1)%pelist)
                   call mpp_send(g_dat,size(g_dat),Atm(1)%pelist(p))
                enddo
             endif
          endif
          if (ANY(Atm(1)%pelist == gid)) then
             call mpp_recv(g_dat, size(g_dat), sending_proc)
          endif
       else
          call gather_grid(Atm(1)%parent_grid, Atm(1)%parent_grid%q(isd_p:ied_p,jsd_p:jed_p,:,nq), g_dat, &
               isd_p, ied_p, jsd_p, jed_p, isg, ieg, jsg, jeg, npz, 0, 0, Atm(1)%parent_tile, &
            concurrent, Atm(1)%pelist, size(Atm(1)%pelist))
       endif
       call timing_off('COMM_TOTAL')
       if (process) call fill_nested_grid(Atm(1)%q(isd:ied,jsd:jed,:,nq), g_dat, &
            Atm(1)%ind_h, Atm(1)%wt_h, &
            0, 0,  isg, ieg, jsg, jeg, npz)

    call mpp_sync_self

    end do

    !Note that we do NOT fill in phis (surface geopotential), which should 
    !be computed exactly instead of being interpolated.


#ifndef SW_DYNAMICS
    !pt --- actually temperature

    call timing_on('COMM_TOTAL')
    if (concurrent) then
       if (ANY(Atm(1)%parent_grid%pelist == gid) .and. Atm(1)%parent_tile == Atm(1)%parent_grid%tile) then
             call mpp_global_field( &
               Atm(1)%parent_grid%domain, &
               Atm(1)%parent_grid%pt(isd_p:ied_p,jsd_p:jed_p,:), g_dat, position=CENTER)
          if (gid == sending_proc) then
             do p=1,size(Atm(1)%pelist)
                call mpp_send(g_dat,size(g_dat),Atm(1)%pelist(p))
             enddo
          endif
       endif
       if (ANY(Atm(1)%pelist == gid)) then
          call mpp_recv(g_dat, size(g_dat), sending_proc)
       endif
    else
       call gather_grid(Atm(1)%parent_grid, Atm(1)%parent_grid%pt(isd_p:ied_p,jsd_p:jed_p,:), g_dat, &
            isd_p, ied_p, jsd_p, jed_p, isg, ieg, jsg, jeg, npz, 0, 0, Atm(1)%parent_tile, &
            concurrent, Atm(1)%pelist, size(Atm(1)%pelist))
    endif
    call mpp_sync_self

    call timing_off('COMM_TOTAL')
    if (process) call fill_nested_grid(Atm(1)%pt, g_dat, &
         Atm(1)%ind_h, Atm(1)%wt_h, &
         0, 0,  isg, ieg, jsg, jeg, npz)


    !Need to fill boundaries with INITIAL coarse-grid potential temperature,
    !because the nested grid cannot do it by itself (values of pkz in the halo are not saved)
    !We will want to INTERPOLATE pkz from the coarse grid and then use it to compute
    !theta on the haloes of the nested grid. Computing theta on the coarse grid and then
    !interpolating yields a different (and inconsistent) answer with the way it is computed
    !in the interior

    if ( Atm(1)%nwat > 0 ) then
       sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
    else
       sphum = 1
    endif
    if ( Atm(1)%parent_grid%adiabatic .or. Atm(1)%parent_grid%do_Held_Suarez ) then
       zvir = 0.         ! no virtual effect
    else
       zvir = rvgas/rdgas - 1.
    endif

    call timing_on('COMM_TOTAL')
    if (concurrent) then
       if (ANY(Atm(1)%parent_grid%pelist == gid) .and. Atm(1)%parent_tile == Atm(1)%parent_grid%tile) then
             call mpp_global_field( &
               Atm(1)%parent_grid%domain, &
               Atm(1)%parent_grid%pkz(isc_p:iec_p,jsc_p:jec_p,:), g_dat, position=CENTER)
          if (gid == sending_proc) then
             do p=1,size(Atm(1)%pelist)
                call mpp_send(g_dat,size(g_dat),Atm(1)%pelist(p))
             enddo
          endif
       endif
       if (ANY(Atm(1)%pelist == gid)) then
          call mpp_recv(g_dat, size(g_dat), sending_proc)
       endif
    else
       call gather_grid(Atm(1)%parent_grid, Atm(1)%parent_grid%pkz(isc_p:iec_p,jsc_p:jec_p,:), g_dat, &
            isc_p, iec_p, jsc_p, jec_p, isg, ieg, jsg, jeg, npz, 0, 0, Atm(1)%parent_tile, &
            concurrent, Atm(1)%pelist, size(Atm(1)%pelist))
    endif
    call mpp_sync_self

    call timing_off('COMM_TOTAL')
    if (process) then 
       allocate(pt_coarse(isd:ied,jsd:jed,npz))
       call fill_nested_grid(pt_coarse, g_dat, &
            Atm(1)%ind_h, Atm(1)%wt_h, &
            0, 0,  isg, ieg, jsg, jeg, npz)

       if (Atm(1)%is == 1) then
          do k=1,npz
             do j=Atm(1)%jsd,Atm(1)%jed
                do i=Atm(1)%isd,0
                   Atm(1)%pt(i,j,k) = cp_air*Atm(1)%pt(i,j,k)/pt_coarse(i,j,k)*(1.+zvir*Atm(1)%q(i,j,k,sphum))
                end do
             end do
          end do
       end if

       if (Atm(1)%js == 1) then
          if (Atm(1)%is == 1) then
             istart = Atm(1)%is
          else
             istart = Atm(1)%isd
          end if
          if (Atm(1)%ie == Atm(1)%npx-1) then
             iend = Atm(1)%ie
          else
             iend = Atm(1)%ied
          end if

          do k=1,npz
             do j=Atm(1)%jsd,0
                do i=istart,iend
                   Atm(1)%pt(i,j,k) = cp_air*Atm(1)%pt(i,j,k)/pt_coarse(i,j,k)*(1.+zvir*Atm(1)%q(i,j,k,sphum))
                end do
             end do
          end do
       end if

       if (Atm(1)%ie == Atm(1)%npx-1) then
          do k=1,npz
             do j=Atm(1)%jsd,Atm(1)%jed
                do i=Atm(1)%npx,Atm(1)%ied
                   Atm(1)%pt(i,j,k) = cp_air*Atm(1)%pt(i,j,k)/pt_coarse(i,j,k)*(1.+zvir*Atm(1)%q(i,j,k,sphum))
                end do
             end do
          end do
       end if

       if (Atm(1)%je == Atm(1)%npy-1) then
          if (Atm(1)%is == 1) then
             istart = Atm(1)%is
          else
             istart = Atm(1)%isd
          end if
          if (Atm(1)%ie == Atm(1)%npx-1) then
             iend = Atm(1)%ie
          else
             iend = Atm(1)%ied
          end if

          do k=1,npz
             do j=Atm(1)%npy,Atm(1)%jed
                do i=istart,iend
                   Atm(1)%pt(i,j,k) = cp_air*Atm(1)%pt(i,j,k)/pt_coarse(i,j,k)*(1.+zvir*Atm(1)%q(i,j,k,sphum))
                end do
             end do
          end do
       end if

       deallocate(pt_coarse)

    end if

    if (.not. Atm(1)%hydrostatic) then

       !delz
       !delz DOES NOT go into the halo!!
       call timing_on('COMM_TOTAL')
       if (concurrent) then
          if (ANY(Atm(1)%parent_grid%pelist == gid) .and. Atm(1)%parent_tile == Atm(1)%parent_grid%tile) then
             call mpp_global_field( &
               Atm(1)%parent_grid%domain, &
               Atm(1)%parent_grid%delz(isd_p:ied_p,jsd_p:jed_p,:), g_dat, position=CENTER)
             if (gid == sending_proc) then
                do p=1,size(Atm(1)%pelist)
                   call mpp_send(g_dat,size(g_dat),Atm(1)%pelist(p))
                enddo
             endif
          endif
          if (ANY(Atm(1)%pelist == gid)) then
             call mpp_recv(g_dat, size(g_dat), sending_proc)
          endif
       else
          call gather_grid(Atm(1)%parent_grid, Atm(1)%parent_grid%delz(isd_p:ied_p,jsd_p:jed_p,:), g_dat, &
               isd_p, ied_p, jsd_p, jed_p, isg, ieg, jsg, jeg, npz, 0, 0, Atm(1)%parent_tile, &
            concurrent, Atm(1)%pelist, size(Atm(1)%pelist))
       endif
    call mpp_sync_self

       call timing_off('COMM_TOTAL')
       if (process) call fill_nested_grid(Atm(1)%delz, g_dat, &
            Atm(1)%ind_h, Atm(1)%wt_h, &
            0, 0,  isg, ieg, jsg, jeg, npz, &
            istart_in=Atm(1)%is, jstart_in=Atm(1)%js, iend_in=Atm(1)%ie, jend_in=Atm(1)%je)

       !w

       call timing_on('COMM_TOTAL')
       if (concurrent) then
          if (ANY(Atm(1)%parent_grid%pelist == gid) .and. Atm(1)%parent_tile == Atm(1)%parent_grid%tile) then
             call mpp_global_field( &
               Atm(1)%parent_grid%domain, &
               Atm(1)%parent_grid%w(isd_p:ied_p,jsd_p:jed_p,:), g_dat, position=CENTER)
             if (gid == sending_proc) then
                do p=1,size(Atm(1)%pelist)
                   call mpp_send(g_dat,size(g_dat),Atm(1)%pelist(p))
                enddo
             endif
          endif
          if (ANY(Atm(1)%pelist == gid)) then
             call mpp_recv(g_dat, size(g_dat), sending_proc)
          endif
       else
          call gather_grid(Atm(1)%parent_grid, Atm(1)%parent_grid%w(isd_p:ied_p,jsd_p:jed_p,:), g_dat, &
               isd_p, ied_p, jsd_p, jed_p, isg, ieg, jsg, jeg, npz, 0, 0, Atm(1)%parent_tile, &
            concurrent, Atm(1)%pelist, size(Atm(1)%pelist))
       endif
    call mpp_sync_self

       call timing_off('COMM_TOTAL')
       if (process) call fill_nested_grid(Atm(1)%w, g_dat, &
            Atm(1)%ind_h, Atm(1)%wt_h, &
            0, 0,  isg, ieg, jsg, jeg, npz)

    end if

#endif
    deallocate(g_dat) 

    !u

    allocate(g_dat( isg:ieg, jsg:jeg+1, npz) )
    g_dat = 1.e25

    call timing_on('COMM_TOTAL')
    if (concurrent) then
       if (ANY(Atm(1)%parent_grid%pelist == gid) .and. Atm(1)%parent_tile == Atm(1)%parent_grid%tile) then
             call mpp_global_field( &
               Atm(1)%parent_grid%domain, &
               Atm(1)%parent_grid%u(isd_p:ied_p,jsd_p:jed_p+1,:), g_dat, position=NORTH)
          if (gid == sending_proc) then
             do p=1,size(Atm(1)%pelist)
                call mpp_send(g_dat,size(g_dat),Atm(1)%pelist(p))
             enddo
          endif
       endif
       if (ANY(Atm(1)%pelist == gid)) then
          call mpp_recv(g_dat, size(g_dat), sending_proc)
       endif
    else
       call gather_grid(Atm(1)%parent_grid, Atm(1)%parent_grid%u(isd_p:ied_p,jsd_p:jed_p+1,:), g_dat, &
            isd_p, ied_p, jsd_p, jed_p, isg, ieg, jsg, jeg, npz, 0, 1, Atm(1)%parent_tile, &
            concurrent, Atm(1)%pelist, size(Atm(1)%pelist))
    endif
    call mpp_sync_self

    call timing_off('COMM_TOTAL')
    call mpp_sync_self
    if (process) call fill_nested_grid(Atm(1)%u, g_dat, &
         Atm(1)%ind_u, Atm(1)%wt_u, &
         0, 1,  isg, ieg, jsg, jeg, npz)
    deallocate(g_dat)

    !v

    allocate(g_dat( isg:ieg+1, jsg:jeg, npz) )
    g_dat = 1.e25

    call timing_on('COMM_TOTAL')
    if (concurrent) then
       if (ANY(Atm(1)%parent_grid%pelist == gid) .and. Atm(1)%parent_tile == Atm(1)%parent_grid%tile) then
             call mpp_global_field( &
               Atm(1)%parent_grid%domain, &
               Atm(1)%parent_grid%v(isd_p:ied_p+1,jsd_p:jed_p,:), g_dat, position=EAST)
          if (gid == sending_proc) then
             do p=1,size(Atm(1)%pelist)
                call mpp_send(g_dat,size(g_dat),Atm(1)%pelist(p))
             enddo
          endif
       endif
       if (ANY(Atm(1)%pelist == gid)) then
          call mpp_recv(g_dat, size(g_dat), sending_proc)
       endif
    else
       call gather_grid(Atm(1)%parent_grid, Atm(1)%parent_grid%v(isd_p:ied_p+1,jsd_p:jed_p,:), g_dat, &
            isd_p, ied_p, jsd_p, jed_p, isg, ieg, jsg, jeg, npz, 1, 0, Atm(1)%parent_tile, &
            concurrent, Atm(1)%pelist, size(Atm(1)%pelist))
    endif
    call mpp_sync_self

    call timing_off('COMM_TOTAL')
    call mpp_sync_self
    if (process) call fill_nested_grid(Atm(1)%v, g_dat, &
         Atm(1)%ind_v, Atm(1)%wt_v, &
         1, 0,  isg, ieg, jsg, jeg, npz)

    deallocate(g_dat)

  end subroutine fill_nested_grid_data

  subroutine fill_nested_grid_data_end(Atm, proc_in)

    type(fv_atmos_type), intent(INOUT) :: Atm  
    logical, intent(IN), OPTIONAL :: proc_in
    real, allocatable :: g_dat(:,:,:), pt_coarse(:,:,:)
    integer :: i,j,k,nq, sphum, ncnst, istart, iend, npz
    integer :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer :: isd_p, ied_p, jsd_p, jed_p, isc_p, iec_p, jsc_p, jec_p
    integer :: isg, ieg, jsg,jeg, npx_p, npy_p
    integer :: isg_n, ieg_n, jsg_n, jeg_n, npx_n, npy_n
    real zvir
    logical :: process
    integer :: p , sending_proc

    if (present(proc_in)) then
       process = proc_in
    else
       process = .true.
    endif

    isd = Atm%isd
    ied = Atm%ied
    jsd = Atm%jsd
    jed = Atm%jed
    ncnst = Atm%ncnst
    isc = Atm%isc; iec = Atm%iec; jsc = Atm%jsc; jec = Atm%jec
    npz     = Atm%npz    
    
    if (concurrent) then
          isd_p = Atm%parent_grid%isd
          ied_p = Atm%parent_grid%ied
          jsd_p = Atm%parent_grid%jsd
          jed_p = Atm%parent_grid%jed
          isc_p = Atm%parent_grid%isc
          iec_p = Atm%parent_grid%iec
          jsc_p = Atm%parent_grid%jsc
          jec_p = Atm%parent_grid%jec
       sending_proc = Atm%parent_grid%pelist(1) + (Atm%parent_tile-1)*Atm%parent_grid%npes_per_tile
    else
       call mpp_get_data_domain( Atm%parent_grid%domain, &
            isd_p,  ied_p,  jsd_p,  jed_p  )
       call mpp_get_compute_domain( Atm%parent_grid%domain, &
            isc_p,  iec_p,  jsc_p,  jec_p  )
       sending_proc = 0
    endif
    call mpp_get_global_domain( Atm%parent_grid%domain, &
         isg, ieg, jsg, jeg, xsize=npx_p, ysize=npy_p)


    !NOW: what we do is to update the nested-grid terrain to the coarse grid,
    !to ensure consistency between the two grids.
    if (Atm%twowaynest) then
       if (ANY(Atm%parent_grid%pelist == gid) .or. ANY(Atm%pelist == gid) .or. .not. concurrent) then
          call update_coarse_grid(Atm%parent_grid%phis, &
               Atm%phis, Atm%nest_domain, &
               Atm%ind_update_h(isd_p:ied_p+1,jsd_p:jed_p+1,:), &
               isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, 0, 0, &
               Atm%refinement, Atm%nestupdate, 0, 0)
       end if

    end if




#ifdef SW_DYNAMICS
!!$    !ps: first level only
!!$    !This is only valid for shallow-water simulations
!!$    if (process) then
!!$    do j=jsd,jed
!!$       do i=isd,ied
!!$
!!$          Atm%ps(i,j) = Atm%delp(i,j,1)/grav
!!$
!!$       end do
!!$    end do
!!$    endif
#else
    !Sets up flow to be initially hydrostatic (shouldn't be the case for all ICs?)
    if (process) call p_var(npz, isc, iec, jsc, jec, Atm%ptop, ptop_min, Atm%delp, &
         Atm%delz, Atm%pt, Atm%ps,   &
         Atm%pe, Atm%peln, Atm%pk, Atm%pkz, kappa, Atm%q, &
         Atm%ng, ncnst, Atm%dry_mass, .false., Atm%mountain, &
         Atm%moist_phys, .true., Atm%nwat)
#endif

 

  end subroutine fill_nested_grid_data_end

  subroutine setup_nested_grid_topography(Atm)

    type(fv_atmos_type), intent(INOUT) :: Atm
    logical hybrid

       if ( Atm%make_hybrid_z ) then
          hybrid = .false.
       else
          hybrid = Atm%hybrid_z
       endif
       if (Atm%grid_type < 4) then
          if ( .not. Atm%external_ic ) then
             call init_case(Atm%u,Atm%v,Atm%w,Atm%pt,Atm%delp, &
                  Atm%q,Atm%phis, Atm%ps,Atm%pe, &
                  Atm%peln,Atm%pk,Atm%pkz, Atm%uc,Atm%vc, Atm%ua,Atm%va,        &
                  Atm%ak, Atm%bk, Atm%npx, Atm%npy, Atm%npz, Atm%ng, Atm%ncnst, Atm%nwat,  &
                  Atm%ndims, Atm%ntiles, Atm%dry_mass, Atm%mountain,       &
                  Atm%moist_phys, Atm%hydrostatic, hybrid, Atm%delz, Atm%ze0)
          endif

       elseif (Atm%grid_type == 4) then
          if ( .not. Atm%external_ic ) then
             call init_double_periodic(Atm%u,Atm%v,Atm%w,Atm%pt,Atm%delp,Atm%q,Atm%phis, &
                  Atm%ps,Atm%pe, &
                  Atm%peln,Atm%pk,Atm%pkz, Atm%uc,Atm%vc, Atm%ua,Atm%va,        &
                  Atm%ak, Atm%bk, Atm%npx, Atm%npy, Atm%npz, Atm%ng, Atm%ncnst, Atm%nwat,  &
                  Atm%ndims, Atm%ntiles, Atm%dry_mass, Atm%mountain,       &
                  Atm%moist_phys, Atm%hydrostatic, hybrid, Atm%delz, Atm%ze0)
          end if
       elseif (Atm%grid_type == 5 .or. Atm%grid_type == 6) then
          call init_latlon(Atm%u,Atm%v,Atm%pt,Atm%delp,Atm%q,Atm%phis, Atm%ps,Atm%pe, &
               Atm%peln,Atm%pk,Atm%pkz, Atm%uc,Atm%vc, Atm%ua,Atm%va,        &
               Atm%ak, Atm%bk, Atm%npx, Atm%npy, Atm%npz, Atm%ng, Atm%ncnst, &
               Atm%ndims, Atm%ntiles, Atm%dry_mass, Atm%mountain,       &
               Atm%moist_phys, hybrid, Atm%delz, Atm%ze0)
       endif

  end subroutine setup_nested_grid_topography

  !#######################################################################
  ! <SUBROUTINE NAME="fv_write_restart">
  ! <DESCRIPTION>
  !  Write out restart files registered through register_restart_file
  ! </DESCRIPTION>
  subroutine fv_write_restart(Atm, timestamp)
    type(fv_atmos_type), intent(inout) :: Atm(:)
    character(len=*),    intent(in)    :: timestamp
    integer :: n

    do n=1,size(Atm)
       call fv_io_write_restart(Atm(n:n), timestamp)
       if (Atm(n)%nested) then
          !FV_IO_WRITE_BCs needs to have the current grid set,
          !so that mp_gather can get the correct tile number.
          call switch_current_Atm(Atm(n))
          call fv_io_write_BCs(Atm(n))
       endif
    enddo
    if (size(Atm) > 1) call switch_current_Atm(Atm(1))

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
    type(fv_atmos_type), intent(inout) :: Atm(:)

    integer :: isc, iec, jsc, jec
    integer :: iq, n, ntileMe
    integer :: isd, ied, jsd, jed, npz
    integer :: unit
    integer :: file_unit
    integer, allocatable :: pelist(:)
    character(len=128):: tracer_name

    ntileMe = size(Atm(:))

    allocate(pelist(mpp_npes()))
    call mpp_get_current_pelist(pelist)
    do n = 1, ntileMe

       if (.not. grids_on_this_pe(n)) then
          cycle
       endif

       !FV_IO_WRITE_BCs needs to have the current grid set,
       !so that mp_gather can get the correct tile number.
       call switch_current_Atm(Atm(n))
       call mpp_set_current_pelist(Atm(n)%pelist)
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
        call get_tracer_names(MODEL_ATMOS, iq, tracer_name)
        write(unit,*) 'fv_restart_end '//trim(tracer_name)//' = ', mpp_chksum(Atm(n)%q(isc:iec,jsc:jec,:,iq))
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

      call fv_io_write_restart(Atm(n:n))
       if (Atm(n)%nested) then
          call mpp_set_current_pelist(Atm(n)%pelist)
          call fv_io_write_BCs(Atm(n))
      endif
    end do

    call mpp_set_current_pelist(pelist)   
    call switch_current_Atm(Atm(1))

    module_is_initialized = .FALSE.

#ifdef EFLUX_OUT
    if( gid==masterproc ) then
        write(*,*) steps, 'Mean equivalent Heat flux for this integration period=',efx_sum/real(max(1,steps)), &
                          'Mean nesting-related flux for this integration period=',efx_sum_nest/real(max(1,steps)), &
                          'Mean mountain torque=',mtq_sum/real(max(1,steps))
        file_unit = get_unit()
        open (unit=file_unit, file='e_flux.data', form='unformatted',status='unknown', access='sequential')
        do n=1,steps
           write(file_unit) efx(n)
           write(file_unit) mtq(n)    ! time series global mountain torque
           !write(file_unit) efx_nest(n)  
        enddo
        close(unit=file_unit)
    endif
#endif

  end subroutine fv_restart_end
  ! </SUBROUTINE> NAME="fv_restart_end"

 subroutine d2c_setup(u, v, uc, vc, dord4, &
      isd,ied,jsd,jed, is,ie,js,je, npx,npy, &
      grid_type, nested, &
      se_corner, sw_corner, ne_corner, nw_corner, &
      rsin_u,rsin_v,cosa_s,rsin2 )

  logical, intent(in):: dord4
  real, intent(in) ::  u(isd:ied,jsd:jed+1)
  real, intent(in) ::  v(isd:ied+1,jsd:jed)
  real, intent(out), dimension(isd:ied+1,jsd:jed  ):: uc
  real, intent(out), dimension(isd:ied  ,jsd:jed+1):: vc
  integer, intent(in) :: isd,ied,jsd,jed, is,ie,js,je, npx,npy,grid_type
  logical, intent(in) :: nested, se_corner, sw_corner, ne_corner, nw_corner
  real, intent(in) :: rsin_u(isd:ied+1,jsd:jed)
  real, intent(in) :: rsin_v(isd:ied,jsd:jed+1)
  real, intent(in) :: cosa_s(isd:ied,jsd:jed)
  real, intent(in) :: rsin2(isd:ied,jsd:jed)

! Local 
  real, dimension(isd:ied,jsd:jed):: utmp, vtmp
  real, parameter:: t11=27./28., t12=-13./28., t13=3./7., t14=6./7., t15=3./28.
  real, parameter:: a1 =  0.5625
  real, parameter:: a2 = -0.0625
  real, parameter:: c1 = -2./14.
  real, parameter:: c2 = 11./14.
  real, parameter:: c3 =  5./14.
  integer npt, i, j, ifirst, ilast, id

  if ( dord4) then
       id = 1
  else
       id = 0
  endif


  if (grid_type < 3 .and. .not. nested) then
     npt = 4
  else
     npt = -2
  endif

  if ( nested) then  

     do j=jsd+1,jed-1
        do i=isd,ied
           utmp(i,j) = a2*(u(i,j-1)+u(i,j+2)) + a1*(u(i,j)+u(i,j+1))
        enddo
     enddo
     do i=isd,ied
        j = jsd
        utmp(i,j) = 0.5*(u(i,j)+u(i,j+1))
        j = jed
        utmp(i,j) = 0.5*(u(i,j)+u(i,j+1))
     end do

     do j=jsd,jed
        do i=isd+1,ied-1
           vtmp(i,j) = a2*(v(i-1,j)+v(i+2,j)) + a1*(v(i,j)+v(i+1,j))
        enddo
        i = isd
        vtmp(i,j) = 0.5*(v(i,j)+v(i+1,j)) 
        i = ied
        vtmp(i,j) = 0.5*(v(i,j)+v(i+1,j))
     enddo

  else

     !----------
     ! Interior:
     !----------

     do j=max(npt,js-1),min(npy-npt,je+1)
        do i=max(npt,isd),min(npx-npt,ied)
           utmp(i,j) = a2*(u(i,j-1)+u(i,j+2)) + a1*(u(i,j)+u(i,j+1))
        enddo
     enddo
     do j=max(npt,jsd),min(npy-npt,jed)
        do i=max(npt,is-1),min(npx-npt,ie+1)
           vtmp(i,j) = a2*(v(i-1,j)+v(i+2,j)) + a1*(v(i,j)+v(i+1,j))
        enddo
     enddo

     !----------
     ! edges:
     !----------
     if (grid_type < 3) then

        if ( js==1 .or. jsd<npt) then
           do j=jsd,npt-1
              do i=isd,ied
                 utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
                 vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
              enddo
           enddo
        endif

        if ( (je+1)==npy .or. jed>=(npy-npt)) then
           do j=npy-npt+1,jed
              do i=isd,ied
                 utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
                 vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
              enddo
           enddo
        endif

        if ( is==1 .or. isd<npt ) then
           do j=max(npt,jsd),min(npy-npt,jed)
              do i=isd,npt-1
                 utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
                 vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
              enddo
           enddo
        endif

        if ( (ie+1)==npx .or. ied>=(npx-npt)) then
           do j=max(npt,jsd),min(npy-npt,jed)
              do i=npx-npt+1,ied
                 utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
                 vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
              enddo
           enddo
        endif

     endif

  end if

! A -> C
!--------------
! Fix the edges
!--------------
! Xdir:
     if( sw_corner ) then
         do i=-2,0
            utmp(i,0) = -vtmp(0,1-i)
         enddo
     endif
     if( se_corner ) then
         do i=0,2
            utmp(npx+i,0) = vtmp(npx,i+1)
         enddo
     endif
     if( ne_corner ) then
         do i=0,2
            utmp(npx+i,npy) = -vtmp(npx,je-i)
         enddo
     endif
     if( nw_corner ) then
         do i=-2,0
            utmp(i,npy) = vtmp(0,je+i)
         enddo
     endif

  if (grid_type < 3 .and. .not. nested) then
     ifirst = max(3,    is-1)
     ilast  = min(npx-2,ie+2)
  else
     ifirst = is-1
     ilast  = ie+2
  endif
!---------------------------------------------
! 4th order interpolation for interior points:
!---------------------------------------------
     do j=js-1,je+1
        do i=ifirst,ilast
           uc(i,j) = a1*(utmp(i-1,j)+utmp(i,j))+a2*(utmp(i-2,j)+utmp(i+1,j))
        enddo
     enddo

     if (grid_type < 3) then
! Xdir:
     if( is==1 .and. .not. nested ) then
        do j=js-1,je+1
           uc(0,j) = c1*utmp(-2,j) + c2*utmp(-1,j) + c3*utmp(0,j) 
           uc(1,j) = ( t14*(utmp( 0,j)+utmp(1,j))    &
                     + t12*(utmp(-1,j)+utmp(2,j))    &
                     + t15*(utmp(-2,j)+utmp(3,j)) )*rsin_u(1,j)
           uc(2,j) = c1*utmp(3,j) + c2*utmp(2,j) + c3*utmp(1,j)
        enddo
     endif

     if( (ie+1)==npx .and. .not. nested ) then
        do j=js-1,je+1
           uc(npx-1,j) = c1*utmp(npx-3,j)+c2*utmp(npx-2,j)+c3*utmp(npx-1,j) 
           uc(npx,j) = (t14*(utmp(npx-1,j)+utmp(npx,j))+      &
                        t12*(utmp(npx-2,j)+utmp(npx+1,j))     &
                      + t15*(utmp(npx-3,j)+utmp(npx+2,j)))*rsin_u(npx,j)
           uc(npx+1,j) = c3*utmp(npx,j)+c2*utmp(npx+1,j)+c1*utmp(npx+2,j) 
        enddo
     endif

     endif

!------
! Ydir:
!------
     if( sw_corner ) then
         do j=-2,0
            vtmp(0,j) = -utmp(1-j,0)
         enddo
     endif
     if( nw_corner ) then
         do j=0,2
            vtmp(0,npy+j) = utmp(j+1,npy)
         enddo
     endif
     if( se_corner ) then
         do j=-2,0
            vtmp(npx,j) = utmp(ie+j,0)
         enddo
     endif
     if( ne_corner ) then
         do j=0,2
            vtmp(npx,npy+j) = -utmp(ie-j,npy)
         enddo
     endif

     if (grid_type < 3) then

     do j=js-1,je+2
      if ( j==1  .and. .not. nested) then
        do i=is-1,ie+1
           vc(i,1) = (t14*(vtmp(i, 0)+vtmp(i,1))    &
                    + t12*(vtmp(i,-1)+vtmp(i,2))    &
                    + t15*(vtmp(i,-2)+vtmp(i,3)))*rsin_v(i,1)
        enddo
      elseif ( (j==0 .or. j==(npy-1))  .and. .not. nested) then
        do i=is-1,ie+1
           vc(i,j) = c1*vtmp(i,j-2) + c2*vtmp(i,j-1) + c3*vtmp(i,j)
        enddo
      elseif ( (j==2 .or. j==(npy+1))  .and. .not. nested) then
        do i=is-1,ie+1
           vc(i,j) = c1*vtmp(i,j+1) + c2*vtmp(i,j) + c3*vtmp(i,j-1)
        enddo
      elseif ( j==npy  .and. .not. nested) then
        do i=is-1,ie+1
           vc(i,npy) = (t14*(vtmp(i,npy-1)+vtmp(i,npy))    &
                      + t12*(vtmp(i,npy-2)+vtmp(i,npy+1))  &
                      + t15*(vtmp(i,npy-3)+vtmp(i,npy+2)))*rsin_v(i,npy)
        enddo
      else
! 4th order interpolation for interior points:
        do i=is-1,ie+1
           vc(i,j) = a2*(vtmp(i,j-2)+vtmp(i,j+1))+a1*(vtmp(i,j-1)+vtmp(i,j))
        enddo
     endif
     enddo
    else
! 4th order interpolation:
       do j=js-1,je+2
          do i=is-1,ie+1
             vc(i,j) = a2*(vtmp(i,j-2)+vtmp(i,j+1))+a1*(vtmp(i,j-1)+vtmp(i,j))
          enddo
       enddo
    endif

  end subroutine d2c_setup

 subroutine d2a_setup(u, v, ua, va, dord4, &
      isd,ied,jsd,jed, is,ie,js,je, npx,npy, &
      grid_type, nested, &
      cosa_s,rsin2 )

  logical, intent(in):: dord4
  real, intent(in) ::  u(isd:ied,jsd:jed+1)
  real, intent(in) ::  v(isd:ied+1,jsd:jed)
  real, intent(out), dimension(isd:ied  ,jsd:jed  ):: ua
  real, intent(out), dimension(isd:ied  ,jsd:jed  ):: va
  integer, intent(in) :: isd,ied,jsd,jed, is,ie,js,je, npx,npy,grid_type
  real, intent(in) :: cosa_s(isd:ied,jsd:jed)
  real, intent(in) :: rsin2(isd:ied,jsd:jed)
  logical, intent(in) :: nested

! Local 
  real, dimension(isd:ied,jsd:jed):: utmp, vtmp
  real, parameter:: t11=27./28., t12=-13./28., t13=3./7., t14=6./7., t15=3./28.
  real, parameter:: a1 =  0.5625
  real, parameter:: a2 = -0.0625
  real, parameter:: c1 = -2./14.
  real, parameter:: c2 = 11./14.
  real, parameter:: c3 =  5./14.
  integer npt, i, j, ifirst, ilast, id

  if ( dord4) then
       id = 1
  else
       id = 0
  endif


  if (grid_type < 3 .and. .not. nested) then
     npt = 4
  else
     npt = -2
  endif

  if ( nested) then  

     do j=jsd+1,jed-1
        do i=isd,ied
           utmp(i,j) = a2*(u(i,j-1)+u(i,j+2)) + a1*(u(i,j)+u(i,j+1))
        enddo
     enddo
     do i=isd,ied
        j = jsd
        utmp(i,j) = 0.5*(u(i,j)+u(i,j+1))
        j = jed
        utmp(i,j) = 0.5*(u(i,j)+u(i,j+1))
     end do

     do j=jsd,jed
        do i=isd+1,ied-1
           vtmp(i,j) = a2*(v(i-1,j)+v(i+2,j)) + a1*(v(i,j)+v(i+1,j))
        enddo
        i = isd
        vtmp(i,j) = 0.5*(v(i,j)+v(i+1,j)) 
        i = ied
        vtmp(i,j) = 0.5*(v(i,j)+v(i+1,j))
     enddo

  else

     !----------
     ! Interior:
     !----------

     do j=max(npt,js-1),min(npy-npt,je+1)
        do i=max(npt,isd),min(npx-npt,ied)
           utmp(i,j) = a2*(u(i,j-1)+u(i,j+2)) + a1*(u(i,j)+u(i,j+1))
        enddo
     enddo
     do j=max(npt,jsd),min(npy-npt,jed)
        do i=max(npt,is-1),min(npx-npt,ie+1)
           vtmp(i,j) = a2*(v(i-1,j)+v(i+2,j)) + a1*(v(i,j)+v(i+1,j))
        enddo
     enddo

     !----------
     ! edges:
     !----------
     if (grid_type < 3) then

        if ( js==1 .or. jsd<npt) then
           do j=jsd,npt-1
              do i=isd,ied
                 utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
                 vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
              enddo
           enddo
        endif

        if ( (je+1)==npy .or. jed>=(npy-npt)) then
           do j=npy-npt+1,jed
              do i=isd,ied
                 utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
                 vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
              enddo
           enddo
        endif

        if ( is==1 .or. isd<npt ) then
           do j=max(npt,jsd),min(npy-npt,jed)
              do i=isd,npt-1
                 utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
                 vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
              enddo
           enddo
        endif

        if ( (ie+1)==npx .or. ied>=(npx-npt)) then
           do j=max(npt,jsd),min(npy-npt,jed)
              do i=npx-npt+1,ied
                 utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
                 vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
              enddo
           enddo
        endif

     endif

  end if



  do j=js-1-id,je+1+id
     do i=is-1-id,ie+1+id
        ua(i,j) = (utmp(i,j)-vtmp(i,j)*cosa_s(i,j)) * rsin2(i,j)
        va(i,j) = (vtmp(i,j)-utmp(i,j)*cosa_s(i,j)) * rsin2(i,j)
     enddo
  enddo

end subroutine d2a_setup

end module fv_restart_mod
