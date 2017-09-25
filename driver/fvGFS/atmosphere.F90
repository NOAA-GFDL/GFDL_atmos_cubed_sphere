!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of fvGFS.                                       *
!*                                                                     *
!* fvGFS is free software; you can redistribute it and/or modify it    *
!* and are expected to follow the terms of the GNU General Public      *
!* License as published by the Free Software Foundation; either        *
!* version 2 of the License, or (at your option) any later version.    *
!*                                                                     *
!* fvGFS is distributed in the hope that it will be useful, but        *
!* WITHOUT ANY WARRANTY; without even the implied warranty of          *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU   *
!* General Public License for more details.                            *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************
module atmosphere_mod
#include <fms_platform.h>

!-----------------------------------------------------------------------
!
! Interface for Cubed_Sphere fv dynamical core
!
!-----------------------------------------------------------------------

!-----------------
! FMS modules:
!-----------------
use block_control_mod,      only: block_control_type
use constants_mod,          only: cp_air, rdgas, grav, rvgas, kappa, pstd_mks
use time_manager_mod,       only: time_type, get_time, set_time, operator(+), &
                                  operator(-)
use fms_mod,                only: file_exist, open_namelist_file,    &
                                  close_file, error_mesg, FATAL,     &
                                  check_nml_error, stdlog,           &
                                  write_version_number,              &
                                  set_domain,   &
                                  mpp_clock_id, mpp_clock_begin,     &
                                  mpp_clock_end, CLOCK_SUBCOMPONENT, &
                                  clock_flag_default, nullify_domain
use mpp_mod,                only: mpp_error, stdout, FATAL, NOTE, &
                                  input_nml_file, mpp_root_pe,    &
                                  mpp_npes, mpp_pe, mpp_chksum,   &
                                  mpp_get_current_pelist,         &
                                  mpp_set_current_pelist
use mpp_parameter_mod,      only: EUPDATE, WUPDATE, SUPDATE, NUPDATE
use mpp_domains_mod,        only: domain2d, mpp_update_domains
use xgrid_mod,              only: grid_box_type
use field_manager_mod,      only: MODEL_ATMOS
use tracer_manager_mod,     only: get_tracer_index, get_number_tracers, &
                                  NO_TRACER
use IPD_typedefs,           only: IPD_data_type, kind_phys
use fv_iau_mod,             only: IAU_external_data_type

!-----------------
! FV core modules:
!-----------------
use fv_arrays_mod,      only: fv_atmos_type, R_GRID
use fv_control_mod,     only: fv_init, fv_end, ngrids
use fv_eta_mod,         only: get_eta_level
use fv_fill_mod,        only: fill_gfs
use fv_dynamics_mod,    only: fv_dynamics
use fv_nesting_mod,     only: twoway_nesting
use fv_diagnostics_mod, only: fv_diag_init, fv_diag, fv_time, prt_maxmin, prt_height
use fv_nggps_diags_mod, only: fv_nggps_diag_init, fv_nggps_diag
use fv_restart_mod,     only: fv_restart, fv_write_restart
use fv_timing_mod,      only: timing_on, timing_off
use fv_mp_mod,          only: switch_current_Atm
use fv_sg_mod,          only: fv_subgrid_z
use fv_update_phys_mod, only: fv_update_phys
use fv_nwp_nudge_mod,   only: fv_nwp_nudge_init, fv_nwp_nudge_end, do_adiabatic_init

use mpp_domains_mod, only:  mpp_get_data_domain, mpp_get_compute_domain
use boundary_mod, only: update_coarse_grid

implicit none
private

!--- driver routines
public :: atmosphere_init, atmosphere_end, atmosphere_restart, &
          atmosphere_dynamics, atmosphere_state_update

!--- utility routines
public :: atmosphere_resolution, atmosphere_grid_bdry, &
          atmosphere_grid_ctr, atmosphere_domain, &
          atmosphere_control_data, atmosphere_pref, &
          atmosphere_diag_axes, atmosphere_etalvls, &
          atmosphere_hgt, atmosphere_scalar_field_halo, &
!rab          atmosphere_tracer_postinit, &
          atmosphere_diss_est, &
          get_bottom_mass, get_bottom_wind,   &
          get_stock_pe, set_atmosphere_pelist

!--- physics/radiation data exchange routines
public :: atmos_phys_driver_statein

!-----------------------------------------------------------------------
! version number of this module
! Include variable "version" to be written to log file.
#include<file_version.h>
character(len=20)   :: mod_name = 'fvGFS/atmosphere_mod'

!---- private data ----
  type (time_type) :: Time_step_atmos
  public Atm, mytile

  !These are convenience variables for local use only, and are set to values in Atm%
  real    :: dt_atmos
  real    :: zvir
  integer :: npx, npy, npz, ncnst, pnats
  integer :: isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed
  integer :: nq                       ! transported tracers
  integer :: sec, seconds, days
  integer :: id_dynam, id_fv_diag, id_subgridz
  logical :: cold_start = .false.       ! read in initial condition

  integer, dimension(:), allocatable :: id_tracerdt_dyn
  integer :: sphum, liq_wat, rainwat, ice_wat, snowwat, graupel  !condensate species

  integer :: mytile = 1
  integer :: p_split = 1
  integer, allocatable :: pelist(:)
  logical, allocatable :: grids_on_this_pe(:)
  type(fv_atmos_type), allocatable, target :: Atm(:)

  integer :: id_udt_dyn, id_vdt_dyn

  real, parameter:: w0_big = 60.  ! to prevent negative w-tracer diffusion

!---dynamics tendencies for use in fv_subgrid_z and during fv_update_phys
  real, allocatable, dimension(:,:,:)   :: u_dt, v_dt, t_dt
  real, allocatable :: pref(:,:), dum1d(:)

  logical :: first_diag = .true.

contains



 subroutine atmosphere_init (Time_init, Time, Time_step, Grid_box, area)
   type (time_type),    intent(in)    :: Time_init, Time, Time_step
   type(grid_box_type), intent(inout) :: Grid_box
   real(kind=kind_phys), pointer, dimension(:,:), intent(inout) :: area
!--- local variables ---
   integer :: i, n
   integer :: itrac
   logical :: do_atmos_nudge
   character(len=32) :: tracer_name, tracer_units
   real :: ps1, ps2

                    call timing_on('ATMOS_INIT')
   allocate(pelist(mpp_npes()))
   call mpp_get_current_pelist(pelist)

   zvir = rvgas/rdgas - 1.

!---- compute physics/atmos time step in seconds ----

   Time_step_atmos = Time_step
   call get_time (Time_step_atmos, sec)
   dt_atmos = real(sec)

!----- initialize FV dynamical core -----
   !NOTE do we still need the second file_exist call?
   cold_start = (.not.file_exist('INPUT/fv_core.res.nc') .and. .not.file_exist('INPUT/fv_core.res.tile1.nc'))

   call fv_init( Atm, dt_atmos, grids_on_this_pe, p_split )  ! allocates Atm components

   do n=1,ngrids
      if (grids_on_this_pe(n)) mytile = n
   enddo

   Atm(mytile)%Time_init = Time_init

!----- write version and namelist to log file -----
   call write_version_number ( 'fvGFS/ATMOSPHERE_MOD', version )

!-----------------------------------

   npx   = Atm(mytile)%npx
   npy   = Atm(mytile)%npy
   npz   = Atm(mytile)%npz
   ncnst = Atm(mytile)%ncnst
   pnats = Atm(mytile)%flagstruct%pnats

   isc = Atm(mytile)%bd%isc
   iec = Atm(mytile)%bd%iec
   jsc = Atm(mytile)%bd%jsc
   jec = Atm(mytile)%bd%jec

   isd = isc - Atm(mytile)%bd%ng
   ied = iec + Atm(mytile)%bd%ng
   jsd = jsc - Atm(mytile)%bd%ng
   jed = jec + Atm(mytile)%bd%ng

   nq = ncnst-pnats
   sphum   = get_tracer_index (MODEL_ATMOS, 'sphum' )
   liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat' )
   ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat' )
   rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat' )
   snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat' )
   graupel = get_tracer_index (MODEL_ATMOS, 'graupel' )

   if (max(sphum,liq_wat,ice_wat,rainwat,snowwat,graupel) > Atm(mytile)%flagstruct%nwat) then
      call mpp_error (FATAL,' atmosphere_init: condensate species are not first in the list of &
                            &tracers defined in the field_table')
   endif

   ! Allocate grid variables to be used to calculate gradient in 2nd order flux exchange
   ! This data is only needed for the COARSEST grid.
   call switch_current_Atm(Atm(mytile))

   allocate(Grid_box%dx    (   isc:iec  , jsc:jec+1))
   allocate(Grid_box%dy    (   isc:iec+1, jsc:jec  ))
   allocate(Grid_box%area  (   isc:iec  , jsc:jec  ))
   allocate(Grid_box%edge_w(              jsc:jec+1))
   allocate(Grid_box%edge_e(              jsc:jec+1))
   allocate(Grid_box%edge_s(   isc:iec+1           ))
   allocate(Grid_box%edge_n(   isc:iec+1           ))
   allocate(Grid_box%en1   (3, isc:iec  , jsc:jec+1))
   allocate(Grid_box%en2   (3, isc:iec+1, jsc:jec  ))
   allocate(Grid_box%vlon  (3, isc:iec  , jsc:jec  ))
   allocate(Grid_box%vlat  (3, isc:iec  , jsc:jec  ))
   Grid_box%dx    (   isc:iec  , jsc:jec+1) = Atm(mytile)%gridstruct%dx    (   isc:iec,   jsc:jec+1)
   Grid_box%dy    (   isc:iec+1, jsc:jec  ) = Atm(mytile)%gridstruct%dy    (   isc:iec+1, jsc:jec  )
   Grid_box%area  (   isc:iec  , jsc:jec  ) = Atm(mytile)%gridstruct%area  (   isc:iec  , jsc:jec  )
   Grid_box%edge_w(              jsc:jec+1) = Atm(mytile)%gridstruct%edge_w(              jsc:jec+1)
   Grid_box%edge_e(              jsc:jec+1) = Atm(mytile)%gridstruct%edge_e(              jsc:jec+1)
   Grid_box%edge_s(   isc:iec+1           ) = Atm(mytile)%gridstruct%edge_s(   isc:iec+1)
   Grid_box%edge_n(   isc:iec+1           ) = Atm(mytile)%gridstruct%edge_n(   isc:iec+1)
   Grid_box%en1   (:, isc:iec  , jsc:jec+1) = Atm(mytile)%gridstruct%en1   (:, isc:iec  , jsc:jec+1)
   Grid_box%en2   (:, isc:iec+1, jsc:jec  ) = Atm(mytile)%gridstruct%en2   (:, isc:iec+1, jsc:jec  )
   do i = 1,3
     Grid_box%vlon  (i, isc:iec  , jsc:jec  ) = Atm(mytile)%gridstruct%vlon  (isc:iec ,  jsc:jec, i )
     Grid_box%vlat  (i, isc:iec  , jsc:jec  ) = Atm(mytile)%gridstruct%vlat  (isc:iec ,  jsc:jec, i )
   enddo
   allocate (area(isc:iec  , jsc:jec  ))
   area(isc:iec,jsc:jec) = Atm(mytile)%gridstruct%area_64(isc:iec,jsc:jec)

!----- allocate and zero out the dynamics (and accumulated) tendencies
   allocate( u_dt(isd:ied,jsd:jed,npz), &
             v_dt(isd:ied,jsd:jed,npz), &
             t_dt(isc:iec,jsc:jec,npz) )
!--- allocate pref
   allocate(pref(npz+1,2), dum1d(npz+1))

   call set_domain ( Atm(mytile)%domain )
   call fv_restart(Atm(mytile)%domain, Atm, dt_atmos, seconds, days, cold_start, Atm(mytile)%gridstruct%grid_type, grids_on_this_pe)

   fv_time = Time

!----- initialize atmos_axes and fv_dynamics diagnostics
       !I've had trouble getting this to work with multiple grids at a time; worth revisiting?
   call fv_diag_init(Atm(mytile:mytile), Atm(mytile)%atmos_axes, Time, npx, npy, npz, Atm(mytile)%flagstruct%p_ref)
   call fv_nggps_diag_init(Atm(mytile:mytile), Atm(mytile)%atmos_axes, Time)

!---------- reference profile -----------
    ps1 = 101325.
    ps2 =  81060.
    pref(npz+1,1) = ps1
    pref(npz+1,2) = ps2
    call get_eta_level ( npz, ps1, pref(1,1), dum1d, Atm(mytile)%ak, Atm(mytile)%bk )
    call get_eta_level ( npz, ps2, pref(1,2), dum1d, Atm(mytile)%ak, Atm(mytile)%bk )

!  --- initialize clocks for dynamics, physics_down and physics_up
   id_dynam     = mpp_clock_id ('FV dy-core',  flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
   id_subgridz  = mpp_clock_id ('FV subgrid_z',flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
   id_fv_diag   = mpp_clock_id ('FV Diag',     flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )

                    call timing_off('ATMOS_INIT')

   if ( Atm(mytile)%flagstruct%na_init>0 ) then
      call nullify_domain ( )
      if ( .not. Atm(mytile)%flagstruct%hydrostatic ) then
           call prt_maxmin('Before adi: W', Atm(mytile)%w, isc, iec, jsc, jec, Atm(mytile)%ng, npz, 1.)
      endif
      call adiabatic_init(zvir,Atm(mytile)%flagstruct%nudge_dz)
      if ( .not. Atm(mytile)%flagstruct%hydrostatic ) then
           call prt_maxmin('After adi: W', Atm(mytile)%w, isc, iec, jsc, jec, Atm(mytile)%ng, npz, 1.)
! Not nested?
           call prt_height('na_ini Z500', isc,iec, jsc,jec, 3, npz, 500.E2, Atm(mytile)%phis, Atm(mytile)%delz,    &
                Atm(mytile)%peln, Atm(mytile)%gridstruct%area_64(isc:iec,jsc:jec), Atm(mytile)%gridstruct%agrid_64(isc:iec,jsc:jec,2))
      endif
   else
      call mpp_error(NOTE,'No adiabatic initialization correction in use')
   endif

#ifdef DEBUG
   call nullify_domain()
   call fv_diag(Atm(mytile:mytile), zvir, Time, -1)
#endif

   n = mytile
   call switch_current_Atm(Atm(n)) 
      
 end subroutine atmosphere_init


 subroutine p_adi(km, ng, ifirst, ilast, jfirst, jlast, ptop,   &
                  delp, pt, ps, pe, peln, pk, pkz, hydrostatic)
! Given (ptop, delp) computes (ps, pk, pe, peln, pkz)
! Input:
   integer,  intent(in):: km, ng
   integer,  intent(in):: ifirst, ilast            ! Longitude strip
   integer,  intent(in):: jfirst, jlast            ! Latitude strip
   logical, intent(in)::  hydrostatic
   real, intent(in):: ptop
   real, intent(in)::   pt(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng, km)
   real, intent(in):: delp(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng, km)
! Output:
   real, intent(out) ::   ps(ifirst-ng:ilast+ng, jfirst-ng:jlast+ng)
   real, intent(out) ::   pk(ifirst:ilast, jfirst:jlast, km+1)
   real, intent(out) ::   pe(ifirst-1:ilast+1,km+1,jfirst-1:jlast+1) ! Ghosted Edge pressure
   real, intent(out) :: peln(ifirst:ilast, km+1, jfirst:jlast)    ! Edge pressure
   real, intent(out) ::  pkz(ifirst:ilast, jfirst:jlast, km)
! Local
   real pek
   integer i, j, k

   pek = ptop ** kappa
!$OMP parallel do default (none) &
!$OMP              shared (ifirst,ilast,jfirst,jlast,km,ptop,pek,pe,pk, &
!$OMP                      ps,delp,peln,hydrostatic,pkz) &
!$OMP             private (j, i, k)
   do j=jfirst,jlast
      do i=ifirst,ilast
         pe(i,1,j) = ptop
         pk(i,j,1) = pek
      enddo

      do k=2,km+1
         do i=ifirst,ilast
            pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
            peln(i,k,j) = log(pe(i,k,j))
            pk(i,j,k) = exp( kappa*peln(i,k,j) )
         enddo
      enddo

      do i=ifirst,ilast
         ps(i,j) = pe(i,km+1,j)
      enddo

      if ( hydrostatic ) then
         do k=1,km
            do i=ifirst,ilast
               pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k))/(kappa*(peln(i,k+1,j)-peln(i,k,j)))
            enddo
         enddo
      endif
   enddo

 end subroutine p_adi


 subroutine atmosphere_dynamics ( Time )
   type(time_type),intent(in) :: Time
   integer :: itrac, n, psc
   integer :: k, w_diff, nt_dyn

!---- Call FV dynamics -----

   call mpp_clock_begin (id_dynam)

   n = mytile
   do psc=1,abs(p_split)
                    call timing_on('fv_dynamics')
!uc/vc only need be same on coarse grid? However BCs do need to be the same
     call fv_dynamics(npx, npy, npz, nq, Atm(n)%ng, dt_atmos/real(abs(p_split)),&
                      Atm(n)%flagstruct%consv_te, Atm(n)%flagstruct%fill,  &
                      Atm(n)%flagstruct%reproduce_sum, kappa, cp_air, zvir,&
                      Atm(n)%ptop, Atm(n)%ks, nq,                          &
                      Atm(n)%flagstruct%n_split, Atm(n)%flagstruct%q_split,&
                      Atm(n)%u, Atm(n)%v, Atm(n)%w, Atm(n)%delz,           &
                      Atm(n)%flagstruct%hydrostatic,                       &
                      Atm(n)%pt, Atm(n)%delp, Atm(n)%q, Atm(n)%ps,         &
                      Atm(n)%pe, Atm(n)%pk, Atm(n)%peln,                   &
                      Atm(n)%pkz, Atm(n)%phis, Atm(n)%q_con,               &
                      Atm(n)%omga, Atm(n)%ua, Atm(n)%va, Atm(n)%uc,        &
                      Atm(n)%vc, Atm(n)%ak, Atm(n)%bk, Atm(n)%mfx,         &
                      Atm(n)%mfy, Atm(n)%cx, Atm(n)%cy, Atm(n)%ze0,        &
                      Atm(n)%flagstruct%hybrid_z,                          &
                      Atm(n)%gridstruct, Atm(n)%flagstruct,                &
                      Atm(n)%neststruct, Atm(n)%idiag, Atm(n)%bd,          &
                      Atm(n)%parent_grid, Atm(n)%domain,Atm(n)%diss_est)

     call timing_off('fv_dynamics')

    if (ngrids > 1 .and. (psc < p_split .or. p_split < 0)) then
       call timing_on('TWOWAY_UPDATE')
       call twoway_nesting(Atm, ngrids, grids_on_this_pe, zvir)
       call timing_off('TWOWAY_UPDATE')
    endif

    end do !p_split
    call mpp_clock_end (id_dynam)

!-----------------------------------------------------
!--- COMPUTE SUBGRID Z
!-----------------------------------------------------
!--- zero out tendencies
    call mpp_clock_begin (id_subgridz)
    u_dt(:,:,:)   = 0.
    v_dt(:,:,:)   = 0.
    t_dt(:,:,:)   = 0.

    w_diff = get_tracer_index (MODEL_ATMOS, 'w_diff' )
    if ( Atm(n)%flagstruct%fv_sg_adj > 0 ) then
      nt_dyn = nq
      if ( w_diff /= NO_TRACER ) then
        nt_dyn = nq - 1
      endif
      call fv_subgrid_z(isd, ied, jsd, jed, isc, iec, jsc, jec, Atm(n)%npz, &
                        nt_dyn, dt_atmos, Atm(n)%flagstruct%fv_sg_adj,      &
                        Atm(n)%flagstruct%nwat, Atm(n)%delp, Atm(n)%pe,     &
                        Atm(n)%peln, Atm(n)%pkz, Atm(n)%pt, Atm(n)%q,       &
                        Atm(n)%ua, Atm(n)%va, Atm(n)%flagstruct%hydrostatic,&
                        Atm(n)%w, Atm(n)%delz, u_dt, v_dt, t_dt, Atm(n)%flagstruct%n_sponge)
    endif

#ifdef USE_Q_DT
    if ( .not. Atm(n)%flagstruct%hydrostatic .and. w_diff /= NO_TRACER ) then
!$OMP parallel do default (none) &
!$OMP              shared (isc, iec, jsc, jec, w_diff, n, Atm, q_dt) &
!$OMP             private (k)
       do k=1, Atm(n)%npz
          Atm(n)%q(isc:iec,jsc:jec,k,w_diff) = Atm(n)%w(isc:iec,jsc:jec,k) + w0_big
          q_dt(:,:,k,w_diff) = 0.
        enddo
    endif
#endif

   call mpp_clock_end (id_subgridz)

 end subroutine atmosphere_dynamics


 subroutine atmosphere_end (Time, Grid_box )!rab, Radiation, Physics)
   type (time_type),      intent(in)    :: Time
   type(grid_box_type),   intent(inout) :: Grid_box
!rab   type (radiation_type), intent(inout) :: Radiation
!rab   type (physics_type),   intent(inout) :: Physics

  ! initialize domains for writing global physics data
   call set_domain ( Atm(mytile)%domain )

   call nullify_domain ( )
   if (first_diag) then
      call timing_on('FV_DIAG')
      call fv_diag(Atm(mytile:mytile), zvir, fv_time, Atm(mytile)%flagstruct%print_freq)
      call fv_nggps_diag(Atm(mytile:mytile), zvir, fv_time)
      first_diag = .false.
      call timing_off('FV_DIAG')
   endif

   call fv_end(Atm, grids_on_this_pe)
   deallocate (Atm)

   deallocate( u_dt, v_dt, t_dt, pref, dum1d )

 end subroutine atmosphere_end



  !#######################################################################
  ! <SUBROUTINE NAME="atmosphere_restart">
  ! <DESCRIPTION>
  !  Write out restart files registered through register_restart_file
  ! </DESCRIPTION>
  subroutine atmosphere_restart(timestamp)
    character(len=*),  intent(in) :: timestamp

    call fv_write_restart(Atm, grids_on_this_pe, timestamp)

  end subroutine atmosphere_restart
  ! </SUBROUTINE>


 subroutine atmosphere_resolution (i_size, j_size, global)
   integer, intent(out)          :: i_size, j_size
   logical, intent(in), optional :: global
   logical :: local

   local = .true.
   if( PRESENT(global) ) local = .NOT.global

   if( local ) then
       i_size = iec - isc + 1
       j_size = jec - jsc + 1
   else
       i_size = npx - 1
       j_size = npy - 1
   end if

 end subroutine atmosphere_resolution


 subroutine atmosphere_pref (p_ref)
   real, dimension(:,:), intent(inout) :: p_ref

   p_ref = pref

 end subroutine atmosphere_pref


 subroutine atmosphere_control_data (i1, i2, j1, j2, kt, p_hydro, hydro)
   integer, intent(out)           :: i1, i2, j1, j2, kt
   logical, intent(out), optional :: p_hydro, hydro
   i1 = Atm(mytile)%bd%isc
   i2 = Atm(mytile)%bd%iec
   j1 = Atm(mytile)%bd%jsc
   j2 = Atm(mytile)%bd%jec
   kt = Atm(mytile)%npz

   if (present(p_hydro)) p_hydro = Atm(mytile)%flagstruct%phys_hydrostatic
   if (present(  hydro))   hydro = Atm(mytile)%flagstruct%hydrostatic

 end subroutine atmosphere_control_data


 subroutine atmosphere_grid_ctr (lon, lat)
!---------------------------------------------------------------
!    returns the longitude and latitude cell centers
!---------------------------------------------------------------
    real(kind=kind_phys), intent(out) :: lon(:,:), lat(:,:)   ! Unit: radian
! Local data:
    integer i,j

    do j=jsc,jec
       do i=isc,iec
          lon(i-isc+1,j-jsc+1) = Atm(mytile)%gridstruct%agrid_64(i,j,1)
          lat(i-isc+1,j-jsc+1) = Atm(mytile)%gridstruct%agrid_64(i,j,2)
       enddo
    end do

 end subroutine atmosphere_grid_ctr


 subroutine atmosphere_grid_bdry (blon, blat, global)
!---------------------------------------------------------------
!    returns the longitude and latitude grid box edges
!    for either the local PEs grid (default) or the global grid
!---------------------------------------------------------------
    real,    intent(out) :: blon(:,:), blat(:,:)   ! Unit: radian
    logical, intent(in), optional :: global
! Local data:
    integer i,j

    if( PRESENT(global) ) then
      if (global) call mpp_error(FATAL, '==> global grid is no longer available &
                               & in the Cubed Sphere')
    endif

    do j=jsc,jec+1
       do i=isc,iec+1
          blon(i-isc+1,j-jsc+1) = Atm(mytile)%gridstruct%grid(i,j,1)
          blat(i-isc+1,j-jsc+1) = Atm(mytile)%gridstruct%grid(i,j,2)
       enddo
    end do

 end subroutine atmosphere_grid_bdry


 subroutine set_atmosphere_pelist ()
   call mpp_set_current_pelist(Atm(mytile)%pelist, no_sync=.TRUE.)
 end subroutine set_atmosphere_pelist


 subroutine atmosphere_domain ( fv_domain, layout )
   type(domain2d), intent(out) :: fv_domain
   integer, intent(out) :: layout(2)
!  returns the domain2d variable associated with the coupling grid
!  note: coupling is done using the mass/temperature grid with no halos

   fv_domain = Atm(mytile)%domain_for_coupler
   layout(1:2) =  Atm(mytile)%layout(1:2)

 end subroutine atmosphere_domain


 subroutine atmosphere_diag_axes ( axes )
   integer, intent(out) :: axes (:)

!----- returns the axis indices for the atmospheric (mass) grid -----
   if ( size(axes(:)) < 0 .or. size(axes(:)) > 4 ) call error_mesg (    &
                               'get_atmosphere_axes in atmosphere_mod', &
                               'size of argument is incorrect', FATAL   )

   axes (1:size(axes(:))) = Atm(mytile)%atmos_axes (1:size(axes(:)))

 end subroutine atmosphere_diag_axes


 subroutine atmosphere_etalvls (ak, bk, flip)
   real(kind=kind_phys), pointer, dimension(:), intent(inout) :: ak, bk
   logical, intent(in) :: flip

   allocate(ak(npz+1))
   allocate(bk(npz+1))

   if (flip) then
     ak(1:npz+1) = Atm(mytile)%ak(npz+1:1:-1)
     bk(1:npz+1) = Atm(mytile)%bk(npz+1:1:-1)
   else
     ak(1:npz+1) = Atm(mytile)%ak(1:npz+1)
     bk(1:npz+1) = Atm(mytile)%bk(1:npz+1)
   endif
 end subroutine atmosphere_etalvls


 subroutine atmosphere_hgt (hgt, position, relative, flip)
   real(kind=kind_phys), pointer, dimension(:,:,:), intent(inout) :: hgt
   character(len=5), intent(in) :: position
   logical, intent(in) :: relative
   logical, intent(in) :: flip
   !--- local variables ---
   integer:: lev, k, j, i
   real(kind=kind_phys), allocatable, dimension(:,:,:) :: z, dz

   if ((position .ne. "layer") .and. (position .ne. "level")) then
     call mpp_error (FATAL, 'atmosphere_hgt:: incorrect position specification')
   endif

   allocate(z(iec-isc+1,jec-jsc+1,npz+1))
   allocate(dz(iec-isc+1,jec-jsc+1,npz))
   z  = 0
   dz = 0 

   if (Atm(mytile)%flagstruct%hydrostatic) then
     !--- generate dz using hydrostatic assumption
     do j = jsc, jec
       do i = isc, iec
         dz(i-isc+1,j-jsc+1,1:npz) = (rdgas/grav)*Atm(mytile)%pt(i,j,1:npz)  &
                         * (Atm(mytile)%peln(i,1:npz,j) - Atm(mytile)%peln(i,2:npz+1,j))
       enddo
     enddo
   else
     !--- use non-hydrostatic delz directly
     do j = jsc, jec
       do i = isc, iec
         dz(i-isc+1,j-jsc+1,1:npz) = Atm(mytile)%delz(i,j,1:npz)
       enddo
     enddo
   endif

   !--- calculate geometric heights at the interfaces (levels)
   !--- if needed, flip the indexing during this step
   if (flip) then
     if (.not. relative) then
       z(:,:,1) = Atm(mytile)%phis(:,:)/grav
     endif
     do k = 2,npz+1
       z(:,:,k) = z(:,:,k-1) - dz(:,:,npz+2-k)
     enddo
   else
     if (.not. relative) then
       z(:,:,npz+1) = Atm(mytile)%phis(:,:)/grav
     endif
     do k = npz,1,-1
       z(:,:,k) = z(:,:,k+1) - dz(:,:,k)
     enddo
   endif

   !--- allocate and set either the level or layer height for return
   if (position == "level") then
     allocate (hgt(iec-isc+1,jec-jsc+1,npz+1))
     hgt = z
   elseif (position == "layer") then
     allocate (hgt(iec-isc+1,jec-jsc+1,npz))
     hgt(:,:,1:npz) = 0.5d0 * (z(:,:,1:npz) + z(:,:,2:npz+1))
   endif

   deallocate (z)
   deallocate (dz)

 end subroutine atmosphere_hgt


 subroutine atmosphere_scalar_field_halo (data, halo, isize, jsize, ksize, data_p)
   !--------------------------------------------------------------------
   ! data   - output array to return the field with halo (i,j,k)
   !          optionally input for field already in (i,j,k) form
   !          sized to include the halo of the field (+ 2*halo)
   ! halo   - size of the halo (must be less than 3)
   ! ied    - horizontal resolution in i-dir with haloes
   ! jed    - horizontal resolution in j-dir with haloes
   ! ksize  - vertical resolution
   ! data_p - optional input field in packed format (ix,k)  
   !--------------------------------------------------------------------
   !--- interface variables ---
   real(kind=kind_phys), dimension(1:isize,1:jsize,ksize), intent(inout) :: data
   integer, intent(in) :: halo
   integer, intent(in) :: isize
   integer, intent(in) :: jsize
   integer, intent(in) :: ksize
   real(kind=kind_phys), dimension(:,:), optional, intent(in) :: data_p
   !--- local variables ---
   integer :: i, j, k
   integer :: ic, jc
   character(len=44) :: modname = 'atmosphere_mod::atmosphere_scalar_field_halo'
   integer :: mpp_flags

   !--- perform error checking
   if (halo .gt. 3) call mpp_error(FATAL, modname//' - halo.gt.3 requires extending the MPP domain to support')
   ic = isize - 2 * halo
   jc = jsize - 2 * halo

   !--- if packed data is present, unpack it into the two-dimensional data array
   if (present(data_p)) then
     if (ic*jc .ne. size(data_p,1)) call mpp_error(FATAL, modname//' - incorrect sizes for incoming &
                                                  &variables data and data_p')
     data = 0.
!$OMP parallel do default (none) &
!$OMP              shared (data, data_p, halo, ic, jc, ksize) &
!$OMP             private (i, j, k)
     do k = 1, ksize
       do j = 1, jc
         do i = 1, ic
           data(i+halo, j+halo, k) = data_p(i + (j-1)*ic, k)
         enddo
       enddo
     enddo
   endif

   mpp_flags = EUPDATE + WUPDATE + SUPDATE + NUPDATE
   if (halo == 1) then
     call mpp_update_domains(data, Atm(mytile)%domain_for_coupler, flags=mpp_flags, complete=.true.)
   elseif (halo == 3) then
     call mpp_update_domains(data, Atm(mytile)%domain, flags=mpp_flags, complete=.true.)
   else
     call mpp_error(FATAL, modname//' - unsupported halo size')
   endif

   !--- fill the halo points when at a corner of the cubed-sphere tile 
   !--- interior domain corners are handled correctly
   if ( (isc==1) .or. (jsc==1) .or. (iec==npx-1) .or. (jec==npy-1) ) then
     do k = 1, ksize
       do j=1,halo
         do i=1,halo
           if ((isc==    1) .and. (jsc==    1)) data(halo+1-j ,halo+1-i ,k) = data(halo+i     ,halo+1-j ,k)  !SW Corner
           if ((isc==    1) .and. (jec==npy-1)) data(halo+1-j ,halo+jc+i,k) = data(halo+i     ,halo+jc+j,k)  !NW Corner
           if ((iec==npx-1) .and. (jsc==    1)) data(halo+ic+j,halo+1-i ,k) = data(halo+ic-i+1,halo+1-j ,k)  !SE Corner
           if ((iec==npx-1) .and. (jec==npy-1)) data(halo+ic+j,halo+jc+i,k) = data(halo+ic-i+1,halo+jc+j,k)  !NE Corner
         enddo
       enddo
     enddo
   endif

   return
 end subroutine atmosphere_scalar_field_halo


 subroutine atmosphere_diss_est (npass)
   use dyn_core_mod, only: del2_cubed
   !--- interface variables ---
   integer, intent(in) :: npass
   !--- local variables
   integer:: k

   !horizontally smooth dissiapation estimate for SKEB
   ! 3 passes before taking absolute value
   do k = 1,min(3,npass)
     call del2_cubed(Atm(mytile)%diss_est, 0.25*Atm(mytile)%gridstruct%da_min, Atm(mytile)%gridstruct, &
                     Atm(mytile)%domain, npx, npy, npz, 3, Atm(mytile)%bd)
   enddo

   Atm(mytile)%diss_est=abs(Atm(mytile)%diss_est)

   do k = 4,npass
     call del2_cubed(Atm(mytile)%diss_est, 0.25*Atm(mytile)%gridstruct%da_min, Atm(mytile)%gridstruct, &
                     Atm(mytile)%domain, npx, npy, npz, 3, Atm(mytile)%bd)
   enddo
   ! provide back sqrt of dissipation estimate
   Atm(mytile)%diss_est=sqrt(Atm(mytile)%diss_est)

 end subroutine atmosphere_diss_est


!--- Need to know the formulation of the mixing ratio being imported into FV3
!--- in order to adjust it in a consistent manner for advection
!rab subroutine atmosphere_tracer_postinit (IPD_Data, Atm_block)
!rab   !--- interface variables ---
!rab   type(IPD_data_type),       intent(in) :: IPD_Data(:)
!rab   type(block_control_type),  intent(in) :: Atm_block
!rab   !--- local variables ---
!rab   integer :: i, j, ix, k, k1, n, nwat, nb, blen
!rab   real(kind=kind_phys) :: qwat
!rab
!rab   if( nq<3 ) call mpp_error(FATAL, 'GFS phys must have 3 interactive tracers')
!rab
!rab   n = mytile
!rab   nwat = Atm(n)%flagstruct%nwat
!rab
!rab!$OMP parallel do default (none) &
!rab!$OMP              shared (Atm_block, Atm, IPD_Data, npz, nq, ncnst, n, nwat) &
!rab!$OMP             private (nb, k, k1, ix, i, j, qwat)
!rab   do nb = 1,Atm_block%nblks
!rab     do k = 1, npz
!rab       k1 = npz+1-k !reverse the k direction
!rab       do ix = 1, Atm_block%blksz(nb)
!rab         i = Atm_block%index(nb)%ii(ix)
!rab         j = Atm_block%index(nb)%jj(ix)
!rab         qwat = sum(Atm(n)%q(i,j,k1,1:nwat))
!rab         Atm(n)%q(i,j,k1,1:nq) = Atm(n)%q(i,j,k1,1:nq) + IPD_Data(nb)%Stateout%gq0(ix,k,1:nq) * (1.0 - qwat)
!rab         if (nq .gt. ncnst) then
!rab           Atm(n)%qdiag(i,j,k1,nq+1:ncnst) = Atm(n)%qdiag(i,j,k1,nq+1:ncnst) + IPD_Data(nb)%Stateout%gq0(ix,k,nq+1:ncnst)
!rab         endif
!rab       enddo
!rab     enddo
!rab   enddo
!rab
!rab   call mpp_update_domains (Atm(n)%q, Atm(n)%domain, complete=.true.)
!rab
!rab   return
!rab end subroutine atmosphere_tracer_postinit


 subroutine get_bottom_mass ( t_bot, tr_bot, p_bot, z_bot, p_surf, slp )
!--------------------------------------------------------------
! returns temp, sphum, pres, height at the lowest model level
! and surface pressure
!--------------------------------------------------------------
   real, intent(out), dimension(isc:iec,jsc:jec):: t_bot, p_bot, z_bot, p_surf
   real, intent(out), optional, dimension(isc:iec,jsc:jec):: slp
   real, intent(out), dimension(isc:iec,jsc:jec,nq):: tr_bot
   integer :: i, j, m, k, kr
   real    :: rrg, sigtop, sigbot
   real, dimension(isc:iec,jsc:jec) :: tref
   real, parameter :: tlaps = 6.5e-3

   rrg  = rdgas / grav

   do j=jsc,jec
      do i=isc,iec
         p_surf(i,j) = Atm(mytile)%ps(i,j)
         t_bot(i,j) = Atm(mytile)%pt(i,j,npz)
         p_bot(i,j) = Atm(mytile)%delp(i,j,npz)/(Atm(mytile)%peln(i,npz+1,j)-Atm(mytile)%peln(i,npz,j))
         z_bot(i,j) = rrg*t_bot(i,j)*(1.+zvir*Atm(mytile)%q(i,j,npz,1)) *  &
                      (1. - Atm(mytile)%pe(i,npz,j)/p_bot(i,j))
      enddo
   enddo

   if ( present(slp) ) then
     ! determine 0.8 sigma reference level
     sigtop = Atm(mytile)%ak(1)/pstd_mks+Atm(mytile)%bk(1)
     do k = 1, npz
        sigbot = Atm(mytile)%ak(k+1)/pstd_mks+Atm(mytile)%bk(k+1)
        if (sigbot+sigtop > 1.6) then
           kr = k
           exit
        endif
        sigtop = sigbot
     enddo
     do j=jsc,jec
        do i=isc,iec
           ! sea level pressure
           tref(i,j) = Atm(mytile)%pt(i,j,kr) * (Atm(mytile)%delp(i,j,kr)/ &
                            ((Atm(mytile)%peln(i,kr+1,j)-Atm(mytile)%peln(i,kr,j))*Atm(mytile)%ps(i,j)))**(-rrg*tlaps)
           slp(i,j) = Atm(mytile)%ps(i,j)*(1.+tlaps*Atm(mytile)%phis(i,j)/(tref(i,j)*grav))**(1./(rrg*tlaps))
        enddo
     enddo
   endif

! Copy tracers
   do m=1,nq
      do j=jsc,jec
         do i=isc,iec
            tr_bot(i,j,m) = Atm(mytile)%q(i,j,npz,m)
         enddo
      enddo
   enddo

 end subroutine get_bottom_mass


 subroutine get_bottom_wind ( u_bot, v_bot )
!-----------------------------------------------------------
! returns u and v on the mass grid at the lowest model level
!-----------------------------------------------------------
   real, intent(out), dimension(isc:iec,jsc:jec):: u_bot, v_bot
   integer i, j

   do j=jsc,jec
      do i=isc,iec
         u_bot(i,j) = Atm(mytile)%u_srf(i,j)
         v_bot(i,j) = Atm(mytile)%v_srf(i,j)
      enddo
   enddo

 end subroutine get_bottom_wind



 subroutine get_stock_pe(index, value)
   integer, intent(in) :: index
   real,   intent(out) :: value

#ifdef USE_STOCK
   include 'stock.inc'
#endif

   real wm(isc:iec,jsc:jec)
   integer i,j,k
   real, pointer :: area(:,:)

   area => Atm(mytile)%gridstruct%area

   select case (index)

#ifdef USE_STOCK
   case (ISTOCK_WATER)
#else
   case (1)
#endif

!----------------------
! Perform vertical sum:
!----------------------
     wm = 0.
     do j=jsc,jec
        do k=1,npz
           do i=isc,iec
! Warning: the following works only with AM2 physics: water vapor; cloud water, cloud ice.
              wm(i,j) = wm(i,j) + Atm(mytile)%delp(i,j,k) * ( Atm(mytile)%q(i,j,k,1) +    &
                                                         Atm(mytile)%q(i,j,k,2) +    &
                                                         Atm(mytile)%q(i,j,k,3) )
           enddo
        enddo
     enddo

!----------------------
! Horizontal sum:
!----------------------
     value = 0.
     do j=jsc,jec
        do i=isc,iec
           value = value + wm(i,j)*area(i,j)
        enddo
     enddo
     value = value/grav

   case default
     value = 0.0
   end select

 end subroutine get_stock_pe


 subroutine atmosphere_state_update (Time, IPD_Data, IAU_Data, Atm_block)
   type(time_type),              intent(in) :: Time
   type(IPD_data_type),          intent(in) :: IPD_Data(:)
   type(IAU_external_data_type), intent(in) :: IAU_Data
   type(block_control_type),     intent(in) :: Atm_block
   !--- local variables ---
   type(time_type) :: Time_prev, Time_next
   integer :: i, j, ix, k, k1, n, w_diff, nt_dyn, iq
   integer :: nb, blen, nwat, dnats, nq_adv
   real(kind=kind_phys):: rcp, q0, qwat(nq), qt, rdt

   Time_prev = Time
   Time_next = Time + Time_step_atmos
   rdt = 1.d0 / dt_atmos

   n = mytile
   nwat = Atm(n)%flagstruct%nwat
   dnats = Atm(mytile)%flagstruct%dnats
   nq_adv = nq - dnats

   if( nq<3 ) call mpp_error(FATAL, 'GFS phys must have 3 interactive tracers')

   if (IAU_Data%in_interval) then
!     IAU increments are in units of 1/sec

!     add analysis increment to u,v,t tendencies
!     directly update delp with analysis increment
      do k = 1, npz
         do j = jsc,jec
            do i = isc,iec
               u_dt(i,j,k) = u_dt(i,j,k) + IAU_Data%ua_inc(i,j,k)
               v_dt(i,j,k) = v_dt(i,j,k) + IAU_Data%va_inc(i,j,k)
               t_dt(i,j,k) = t_dt(i,j,k) + IAU_Data%temp_inc(i,j,k)
               Atm(n)%delp(i,j,k) = Atm(n)%delp(i,j,k)  + IAU_Data%delp_inc(i,j,k)*dt_atmos
          enddo
        enddo
      enddo
!     add analysis increment to tracers to output from physics
      do nb = 1,Atm_block%nblks
         !if (nb.EQ.1) print*,'in block_update',IAU_Data%in_interval,IAU_Data%temp_inc(isc,jsc,30)
         blen = Atm_block%blksz(nb)
         do k = 1, npz
            k1 = npz+1-k !reverse the k direction 
            do ix = 1, blen
               i = Atm_block%index(nb)%ii(ix)
               j = Atm_block%index(nb)%jj(ix)
               IPD_Data(nb)%Stateout%gq0(ix,k,:) = IPD_Data(nb)%Stateout%gq0(ix,k,:) + IAU_Data%tracer_inc(i,j,k1,:)*dt_atmos
            enddo
         enddo
      enddo
   endif

   call set_domain ( Atm(mytile)%domain )

   call timing_on('GFS_TENDENCIES')
!--- put u/v tendencies into haloed arrays u_dt and v_dt
!$OMP parallel do default (none) & 
!$OMP              shared (rdt, n, nq, dnats, npz, ncnst, nwat, mytile, u_dt, v_dt, t_dt,&
!$OMP                      Atm, IPD_Data, Atm_block, sphum, liq_wat, rainwat, ice_wat,   &
!$OMP                      snowwat, graupel, nq_adv)   &
!$OMP             private (nb, blen, i, j, k, k1, ix, q0, qwat, qt)
   do nb = 1,Atm_block%nblks

!SJL: perform vertical filling to fix the negative humidity if the SAS convection scheme is used
!     This call may be commented out if RAS or other positivity-preserving CPS is used.
     blen = Atm_block%blksz(nb)
     call fill_gfs(blen, npz, IPD_Data(nb)%Statein%prsi, IPD_Data(nb)%Stateout%gq0, 1.e-9_kind_phys)

     do k = 1, npz
       k1 = npz+1-k !reverse the k direction 
       do ix = 1, blen
         i = Atm_block%index(nb)%ii(ix)
         j = Atm_block%index(nb)%jj(ix)
         u_dt(i,j,k1) = u_dt(i,j,k1) + (IPD_Data(nb)%Stateout%gu0(ix,k) - IPD_Data(nb)%Statein%ugrs(ix,k)) * rdt
         v_dt(i,j,k1) = v_dt(i,j,k1) + (IPD_Data(nb)%Stateout%gv0(ix,k) - IPD_Data(nb)%Statein%vgrs(ix,k)) * rdt
!         t_dt(i,j,k1) = (IPD_Data(nb)%Stateout%gt0(ix,k) - IPD_Data(nb)%Statein%tgrs(ix,k)) * rdt
         t_dt(i,j,k1) = t_dt(i,j,k1) + (IPD_Data(nb)%Stateout%gt0(ix,k) - IPD_Data(nb)%Statein%tgrs(ix,k)) * rdt
! SJL notes:
! ---- DO not touch the code below; dry mass conservation may change due to 64bit <-> 32bit conversion
! GFS total air mass = dry_mass + water_vapor (condensate excluded)
! GFS mixing ratios  = tracer_mass / (dry_mass + vapor_mass)
! FV3 total air mass = dry_mass + [water_vapor + condensate ]
! FV3 mixing ratios  = tracer_mass / (dry_mass+vapor_mass+cond_mass)
         q0 = IPD_Data(nb)%Statein%prsi(ix,k) - IPD_Data(nb)%Statein%prsi(ix,k+1)
         qwat(1:nq_adv) = q0*IPD_Data(nb)%Stateout%gq0(ix,k,1:nq_adv)
! **********************************************************************************************************
! Dry mass: the following way of updating delp is key to mass conservation with hybrid 32-64 bit computation
! **********************************************************************************************************
! The following example is for 2 water species. 
!        q0 = Atm(n)%delp(i,j,k1)*(1.-(Atm(n)%q(i,j,k1,1)+Atm(n)%q(i,j,k1,2))) + q1 + q2
         qt = sum(qwat(1:nwat))
         q0 = Atm(n)%delp(i,j,k1)*(1.-sum(Atm(n)%q(i,j,k1,1:nwat))) + qt 
         Atm(n)%delp(i,j,k1) = q0
         Atm(n)%q(i,j,k1,1:nq_adv) = qwat(1:nq_adv) / q0
!        if (dnats .gt. 0) Atm(n)%q(i,j,k1,nq_adv+1:nq) = IPD_Data(nb)%Stateout%gq0(ix,k,nq_adv+1:nq)
       enddo
     enddo

     !--- diagnostic tracers are assumed to be updated in-place
     !--- SHOULD THESE DIAGNOSTIC TRACERS BE MASS ADJUSTED???
     !--- See Note in statein...
     do iq = nq+1, ncnst
       do k = 1, npz
         k1 = npz+1-k !reverse the k direction 
         do ix = 1, blen
           i = Atm_block%index(nb)%ii(ix)
           j = Atm_block%index(nb)%jj(ix)
           Atm(mytile)%qdiag(i,j,k1,iq) = IPD_Data(nb)%Stateout%gq0(ix,k,iq)
         enddo
       enddo
     enddo

   enddo  ! nb-loop

   call timing_off('GFS_TENDENCIES')

   w_diff = get_tracer_index (MODEL_ATMOS, 'w_diff' )
   nt_dyn = ncnst-pnats   !nothing more than nq
   if ( w_diff /= NO_TRACER ) then
      nt_dyn = nt_dyn - 1
   endif

!--- adjust w and heat tendency for non-hydrostatic case
#ifdef USE_Q_DT
    if ( .not.Atm(n)%flagstruct%hydrostatic .and. w_diff /= NO_TRACER ) then
      rcp = 1. / cp_air
!$OMP parallel do default (none) &
!$OMP              shared (jsc, jec, isc, iec, n, w_diff, Atm, q_dt, t_dt, rcp, dt_atmos) &
!$OMP             private (i, j, k)
       do k=1, Atm(n)%npz
         do j=jsc, jec
           do i=isc, iec
             Atm(n)%q(i,j,k,w_diff) = q_dt(i,j,k,w_diff) ! w tendency due to phys
! Heating due to loss of KE (vertical diffusion of w)
             t_dt(i,j,k) = t_dt(i,j,k) - q_dt(i,j,k,w_diff)*rcp*&
                                     (Atm(n)%w(i,j,k)+0.5*dt_atmos*q_dt(i,j,k,w_diff))
             Atm(n)%w(i,j,k) = Atm(n)%w(i,j,k) + dt_atmos*Atm(n)%q(i,j,k,w_diff)
           enddo
         enddo
       enddo
    endif
#endif

   call mpp_clock_begin (id_dynam)
       call timing_on('FV_UPDATE_PHYS')
    call fv_update_phys( dt_atmos, isc, iec, jsc, jec, isd, ied, jsd, jed, Atm(n)%ng, nt_dyn, &
                         Atm(n)%u,  Atm(n)%v,   Atm(n)%w,  Atm(n)%delp, Atm(n)%pt,         &
                         Atm(n)%q,  Atm(n)%qdiag,                                          &
                         Atm(n)%ua, Atm(n)%va,  Atm(n)%ps, Atm(n)%pe,   Atm(n)%peln,       &
                         Atm(n)%pk, Atm(n)%pkz, Atm(n)%ak, Atm(n)%bk,   Atm(n)%phis,       &
                         Atm(n)%u_srf, Atm(n)%v_srf, Atm(n)%ts, Atm(n)%delz,               &
                         Atm(n)%flagstruct%hydrostatic, u_dt, v_dt, t_dt,                  &
                         .true., Time_next, Atm(n)%flagstruct%nudge, Atm(n)%gridstruct,    &
                         Atm(n)%gridstruct%agrid(:,:,1), Atm(n)%gridstruct%agrid(:,:,2),   &
                         Atm(n)%npx, Atm(n)%npy, Atm(n)%npz, Atm(n)%flagstruct,            &
                         Atm(n)%neststruct, Atm(n)%bd, Atm(n)%domain, Atm(n)%ptop)
       call timing_off('FV_UPDATE_PHYS')
   call mpp_clock_end (id_dynam)

!--- nesting update after updating atmospheric variables with
!--- physics tendencies
    if (ngrids > 1 .and. p_split > 0) then
       call timing_on('TWOWAY_UPDATE')
       call twoway_nesting(Atm, ngrids, grids_on_this_pe, zvir)
       call timing_off('TWOWAY_UPDATE')
    endif   
   call nullify_domain()

  !---- diagnostics for FV dynamics -----
   if (Atm(mytile)%flagstruct%print_freq /= -99) then
     call mpp_clock_begin(id_fv_diag)

     fv_time = Time_next
     call get_time (fv_time, seconds,  days)

     call nullify_domain()
     call timing_on('FV_DIAG')
     call fv_diag(Atm(mytile:mytile), zvir, fv_time, Atm(mytile)%flagstruct%print_freq)
     first_diag = .false.

     fv_time = Time_next - Atm(n)%Time_init
     call get_time (fv_time, seconds,  days)
    !--- perform diagnostics on GFS fdiag schedule
     if (ANY(nint(Atm(mytile)%fdiag(:)*3600.) == (days*24*3600+seconds)) .or. (days*24*3600+seconds) == dt_atmos ) then
       if (mpp_pe() == mpp_root_pe()) write(6,*) 'NGGPS:FV3 DIAG STEP', (real(days)*24. + real(seconds)/3600.)
       call fv_nggps_diag(Atm(mytile:mytile), zvir, Time_next)
     endif
     call timing_off('FV_DIAG')

     call mpp_clock_end(id_fv_diag)
   endif

 end subroutine atmosphere_state_update


 subroutine adiabatic_init(zvir,nudge_dz)
   real, allocatable, dimension(:,:,:):: u0, v0, t0, dz0, dp0
   real, intent(in):: zvir
   logical, intent(inout):: nudge_dz
!  real, parameter:: wt = 1.  ! was 2.
   real, parameter:: wt = 2.
!***********
! Haloe Data
!***********
   real, parameter::    q1_h2o = 2.2E-6
   real, parameter::    q7_h2o = 3.8E-6
   real, parameter::  q100_h2o = 3.8E-6
   real, parameter:: q1000_h2o = 3.1E-6
   real, parameter:: q2000_h2o = 2.8E-6
   real, parameter:: q3000_h2o = 3.0E-6
   real:: xt, p00, q00
   integer:: isc, iec, jsc, jec, npz
   integer:: m, n, i,j,k, ngc

   character(len=80) :: errstr

   xt = 1./(1.+wt)

   write(errstr,'(A, I4, A)') 'Performing adiabatic init',  Atm(mytile)%flagstruct%na_init, ' times'
   call mpp_error(NOTE, errstr)
   sphum = get_tracer_index (MODEL_ATMOS, 'sphum' )

    npz = Atm(mytile)%npz

    isc = Atm(mytile)%bd%isc
    iec = Atm(mytile)%bd%iec
    jsc = Atm(mytile)%bd%jsc
    jec = Atm(mytile)%bd%jec

    ngc = Atm(mytile)%ng
    isd = isc - ngc
    ied = iec + ngc
    jsd = jsc - ngc
    jed = jec + ngc

     call timing_on('adiabatic_init')
     do_adiabatic_init = .true.

     allocate ( u0(isc:iec,  jsc:jec+1, npz) )
     allocate ( v0(isc:iec+1,jsc:jec,   npz) )
     allocate (dp0(isc:iec,jsc:jec, npz) )

     if ( Atm(mytile)%flagstruct%hydrostatic ) nudge_dz = .false.

     if ( nudge_dz ) then
          allocate (dz0(isc:iec,jsc:jec, npz) )
     else
          allocate ( t0(isc:iec,jsc:jec, npz) )
     endif

!$omp parallel do default (none) & 
!$omp              shared (nudge_dz, npz, jsc, jec, isc, iec, n, sphum, u0, v0, t0, dz0, dp0, Atm, zvir, mytile) &
!$omp             private (k, j, i) 
       do k=1,npz
          do j=jsc,jec+1
             do i=isc,iec
                u0(i,j,k) = Atm(mytile)%u(i,j,k)
             enddo
          enddo
          do j=jsc,jec
             do i=isc,iec+1
                v0(i,j,k) = Atm(mytile)%v(i,j,k)
             enddo
          enddo
          if ( nudge_dz ) then
             do j=jsc,jec
                do i=isc,iec
                   dp0(i,j,k) = Atm(mytile)%delp(i,j,k)
                   dz0(i,j,k) = Atm(mytile)%delz(i,j,k)
                enddo
             enddo
          else
             do j=jsc,jec
                do i=isc,iec
                   t0(i,j,k) = Atm(mytile)%pt(i,j,k)*(1.+zvir*Atm(mytile)%q(i,j,k,sphum))  ! virt T
                   dp0(i,j,k) = Atm(mytile)%delp(i,j,k)
                enddo
             enddo
          endif
       enddo

     do m=1,Atm(mytile)%flagstruct%na_init
! Forward call
    call fv_dynamics(Atm(mytile)%npx, Atm(mytile)%npy, npz,  nq, Atm(mytile)%ng, dt_atmos, 0.,      &
                     Atm(mytile)%flagstruct%fill, Atm(mytile)%flagstruct%reproduce_sum, kappa, cp_air, zvir,  &
                     Atm(mytile)%ptop, Atm(mytile)%ks, nq, Atm(mytile)%flagstruct%n_split,        &
                     Atm(mytile)%flagstruct%q_split, Atm(mytile)%u, Atm(mytile)%v, Atm(mytile)%w,         &
                     Atm(mytile)%delz, Atm(mytile)%flagstruct%hydrostatic,                      & 
                     Atm(mytile)%pt, Atm(mytile)%delp, Atm(mytile)%q, Atm(mytile)%ps,                     &
                     Atm(mytile)%pe, Atm(mytile)%pk, Atm(mytile)%peln, Atm(mytile)%pkz, Atm(mytile)%phis,      &
                     Atm(mytile)%q_con, Atm(mytile)%omga, Atm(mytile)%ua, Atm(mytile)%va, Atm(mytile)%uc, Atm(mytile)%vc, &
                     Atm(mytile)%ak, Atm(mytile)%bk, Atm(mytile)%mfx, Atm(mytile)%mfy,                    &
                     Atm(mytile)%cx, Atm(mytile)%cy, Atm(mytile)%ze0, Atm(mytile)%flagstruct%hybrid_z,    &
                     Atm(mytile)%gridstruct, Atm(mytile)%flagstruct,                            &
                     Atm(mytile)%neststruct, Atm(mytile)%idiag, Atm(mytile)%bd, Atm(mytile)%parent_grid,  &
                     Atm(mytile)%domain,Atm(mytile)%diss_est)
! Backward
    call fv_dynamics(Atm(mytile)%npx, Atm(mytile)%npy, npz,  nq, Atm(mytile)%ng, -dt_atmos, 0.,      &
                     Atm(mytile)%flagstruct%fill, Atm(mytile)%flagstruct%reproduce_sum, kappa, cp_air, zvir,  &
                     Atm(mytile)%ptop, Atm(mytile)%ks, nq, Atm(mytile)%flagstruct%n_split,        &
                     Atm(mytile)%flagstruct%q_split, Atm(mytile)%u, Atm(mytile)%v, Atm(mytile)%w,         &
                     Atm(mytile)%delz, Atm(mytile)%flagstruct%hydrostatic,                      & 
                     Atm(mytile)%pt, Atm(mytile)%delp, Atm(mytile)%q, Atm(mytile)%ps,                     &
                     Atm(mytile)%pe, Atm(mytile)%pk, Atm(mytile)%peln, Atm(mytile)%pkz, Atm(mytile)%phis,      &
                     Atm(mytile)%q_con, Atm(mytile)%omga, Atm(mytile)%ua, Atm(mytile)%va, Atm(mytile)%uc, Atm(mytile)%vc, &
                     Atm(mytile)%ak, Atm(mytile)%bk, Atm(mytile)%mfx, Atm(mytile)%mfy,                    &
                     Atm(mytile)%cx, Atm(mytile)%cy, Atm(mytile)%ze0, Atm(mytile)%flagstruct%hybrid_z,    &
                     Atm(mytile)%gridstruct, Atm(mytile)%flagstruct,                            &
                     Atm(mytile)%neststruct, Atm(mytile)%idiag, Atm(mytile)%bd, Atm(mytile)%parent_grid,  &
                     Atm(mytile)%domain,Atm(mytile)%diss_est)
!Nudging back to IC
!$omp parallel do default (none) &
!$omp              shared (pref, npz, jsc, jec, isc, iec, n, sphum, Atm, u0, v0, t0, dp0, xt, zvir, mytile, nudge_dz, dz0) &
!$omp             private (i, j, k, p00, q00)
       do k=1,npz
          do j=jsc,jec+1
             do i=isc,iec
                Atm(mytile)%u(i,j,k) = xt*(Atm(mytile)%u(i,j,k) + wt*u0(i,j,k))
             enddo
          enddo
          do j=jsc,jec
             do i=isc,iec+1
                Atm(mytile)%v(i,j,k) = xt*(Atm(mytile)%v(i,j,k) + wt*v0(i,j,k))
             enddo
          enddo
          if( Atm(mytile)%flagstruct%nudge_qv ) then
! SJL note: Nudging water vaport towards HALOE climatology:
! In case of better IC (IFS) this step may not be necessary
             p00 = Atm(mytile)%pe(isc,k,jsc)
             if ( p00 < 30.E2 ) then
                if ( p00 < 1. ) then
                     q00 = q1_h2o
                elseif ( p00 <= 7. .and. p00 >= 1. ) then
                     q00 = q1_h2o + (q7_h2o-q1_h2o)*log(pref(k,1)/1.)/log(7.)
                elseif ( p00 < 100. .and. p00 >= 7. ) then
                     q00 = q7_h2o + (q100_h2o-q7_h2o)*log(pref(k,1)/7.)/log(100./7.)
                elseif ( p00 < 1000. .and. p00 >= 100. ) then
                     q00 = q100_h2o + (q1000_h2o-q100_h2o)*log(pref(k,1)/1.E2)/log(10.)
                elseif ( p00 < 2000. .and. p00 >= 1000. ) then
                     q00 = q1000_h2o + (q2000_h2o-q1000_h2o)*log(pref(k,1)/1.E3)/log(2.)
                else
                     q00 = q2000_h2o + (q3000_h2o-q2000_h2o)*log(pref(k,1)/2.E3)/log(1.5)
                endif
                do j=jsc,jec
                   do i=isc,iec
                      Atm(mytile)%q(i,j,k,sphum) = xt*(Atm(mytile)%q(i,j,k,sphum) + wt*q00)
                   enddo
                enddo
             endif
          endif
          if ( nudge_dz ) then
             do j=jsc,jec
                do i=isc,iec
                   Atm(mytile)%delp(i,j,k) = xt*(Atm(mytile)%delp(i,j,k) + wt*dp0(i,j,k))
                   Atm(mytile)%delz(i,j,k) = xt*(Atm(mytile)%delz(i,j,k) + wt*dz0(i,j,k))
                enddo
             enddo
          else
             do j=jsc,jec
                do i=isc,iec
                   Atm(mytile)%pt(i,j,k) = xt*(Atm(mytile)%pt(i,j,k) + wt*t0(i,j,k)/(1.+zvir*Atm(mytile)%q(i,j,k,sphum)))
                   Atm(mytile)%delp(i,j,k) = xt*(Atm(mytile)%delp(i,j,k) + wt*dp0(i,j,k))
                enddo
             enddo
          endif

       enddo

! Backward
    call fv_dynamics(Atm(mytile)%npx, Atm(mytile)%npy, npz,  nq, Atm(mytile)%ng, -dt_atmos, 0.,      &
                     Atm(mytile)%flagstruct%fill, Atm(mytile)%flagstruct%reproduce_sum, kappa, cp_air, zvir,  &
                     Atm(mytile)%ptop, Atm(mytile)%ks, nq, Atm(mytile)%flagstruct%n_split,        &
                     Atm(mytile)%flagstruct%q_split, Atm(mytile)%u, Atm(mytile)%v, Atm(mytile)%w,         &
                     Atm(mytile)%delz, Atm(mytile)%flagstruct%hydrostatic,                      & 
                     Atm(mytile)%pt, Atm(mytile)%delp, Atm(mytile)%q, Atm(mytile)%ps,                     &
                     Atm(mytile)%pe, Atm(mytile)%pk, Atm(mytile)%peln, Atm(mytile)%pkz, Atm(mytile)%phis,      &
                     Atm(mytile)%q_con, Atm(mytile)%omga, Atm(mytile)%ua, Atm(mytile)%va, Atm(mytile)%uc, Atm(mytile)%vc, &
                     Atm(mytile)%ak, Atm(mytile)%bk, Atm(mytile)%mfx, Atm(mytile)%mfy,                    &
                     Atm(mytile)%cx, Atm(mytile)%cy, Atm(mytile)%ze0, Atm(mytile)%flagstruct%hybrid_z,    &
                     Atm(mytile)%gridstruct, Atm(mytile)%flagstruct,                            &
                     Atm(mytile)%neststruct, Atm(mytile)%idiag, Atm(mytile)%bd, Atm(mytile)%parent_grid,  &
                     Atm(mytile)%domain,Atm(mytile)%diss_est)
! Forward call
    call fv_dynamics(Atm(mytile)%npx, Atm(mytile)%npy, npz,  nq, Atm(mytile)%ng, dt_atmos, 0.,      &
                     Atm(mytile)%flagstruct%fill, Atm(mytile)%flagstruct%reproduce_sum, kappa, cp_air, zvir,  &
                     Atm(mytile)%ptop, Atm(mytile)%ks, nq, Atm(mytile)%flagstruct%n_split,        &
                     Atm(mytile)%flagstruct%q_split, Atm(mytile)%u, Atm(mytile)%v, Atm(mytile)%w,         &
                     Atm(mytile)%delz, Atm(mytile)%flagstruct%hydrostatic,                      & 
                     Atm(mytile)%pt, Atm(mytile)%delp, Atm(mytile)%q, Atm(mytile)%ps,                     &
                     Atm(mytile)%pe, Atm(mytile)%pk, Atm(mytile)%peln, Atm(mytile)%pkz, Atm(mytile)%phis,      &
                     Atm(mytile)%q_con, Atm(mytile)%omga, Atm(mytile)%ua, Atm(mytile)%va, Atm(mytile)%uc, Atm(mytile)%vc, &
                     Atm(mytile)%ak, Atm(mytile)%bk, Atm(mytile)%mfx, Atm(mytile)%mfy,                    &
                     Atm(mytile)%cx, Atm(mytile)%cy, Atm(mytile)%ze0, Atm(mytile)%flagstruct%hybrid_z,    &
                     Atm(mytile)%gridstruct, Atm(mytile)%flagstruct,                            &
                     Atm(mytile)%neststruct, Atm(mytile)%idiag, Atm(mytile)%bd, Atm(mytile)%parent_grid,  &
                     Atm(mytile)%domain,Atm(mytile)%diss_est)
! Nudging back to IC
!$omp parallel do default (none) &
!$omp              shared (nudge_dz,npz, jsc, jec, isc, iec, n, sphum, Atm, u0, v0, t0, dz0, dp0, xt, zvir, mytile) &
!$omp             private (i, j, k)
       do k=1,npz
          do j=jsc,jec+1
             do i=isc,iec
                Atm(mytile)%u(i,j,k) = xt*(Atm(mytile)%u(i,j,k) + wt*u0(i,j,k))
             enddo
          enddo
          do j=jsc,jec
             do i=isc,iec+1
                Atm(mytile)%v(i,j,k) = xt*(Atm(mytile)%v(i,j,k) + wt*v0(i,j,k))
             enddo
          enddo
          if ( nudge_dz ) then
             do j=jsc,jec
             do i=isc,iec
                Atm(mytile)%delp(i,j,k) = xt*(Atm(mytile)%delp(i,j,k) + wt*dp0(i,j,k))
                Atm(mytile)%delz(i,j,k) = xt*(Atm(mytile)%delz(i,j,k) + wt*dz0(i,j,k))
             enddo
             enddo
          else
             do j=jsc,jec
             do i=isc,iec
                Atm(mytile)%pt(i,j,k) = xt*(Atm(mytile)%pt(i,j,k) + wt*t0(i,j,k)/(1.+zvir*Atm(mytile)%q(i,j,k,sphum)))
                Atm(mytile)%delp(i,j,k) = xt*(Atm(mytile)%delp(i,j,k) + wt*dp0(i,j,k))
             enddo
             enddo
          endif
       enddo

     enddo

     deallocate ( u0 )
     deallocate ( v0 )
     deallocate (dp0 )
     if ( allocated(t0) )  deallocate ( t0 )
     if ( allocated(dz0) ) deallocate ( dz0 )

     do_adiabatic_init = .false.
     call timing_off('adiabatic_init')

 end subroutine adiabatic_init



#if defined(OVERLOAD_R4)
#define _DBL_(X) DBLE(X)
#define _RL_(X) REAL(X,KIND=4)
#else
#define _DBL_(X) X
#define _RL_(X) X
#endif
 subroutine atmos_phys_driver_statein (IPD_Data, Atm_block)
   type (IPD_data_type),      intent(inout) :: IPD_Data(:)
   type (block_control_type), intent(in)    :: Atm_block
!--------------------------------------
! Local GFS-phys consistent parameters:
!--------------------------------------
   real(kind=kind_phys), parameter:: p00 = 1.e5
   real(kind=kind_phys), parameter:: qmin = 1.0e-10   
   real(kind=kind_phys):: pk0inv, ptop, pktop
   real(kind=kind_phys) :: rTv, dm, qgrs_rad
   integer :: nb, blen, npz, i, j, k, ix, k1, dnats, nq_adv

!!! NOTES: lmh 6nov15
!!! - "Layer" means "layer mean", ie. the average value in a layer
!!! - "Level" means "level interface", ie the point values at the top or bottom of a layer

   ptop =  _DBL_(_RL_(Atm(mytile)%ak(1)))
   pktop  = (ptop/p00)**kappa
   pk0inv = (1.0_kind_phys/p00)**kappa

   npz = Atm_block%npz
   dnats = Atm(mytile)%flagstruct%dnats
   nq_adv = nq - dnats

!---------------------------------------------------------------------
! use most up to date atmospheric properties when running serially
!---------------------------------------------------------------------
!$OMP parallel do default (none) & 
!$OMP             shared  (Atm_block, Atm, IPD_Data, npz, nq, ncnst, sphum, liq_wat, &
!$OMP                      ice_wat, rainwat, snowwat, graupel, pk0inv, ptop,   &
!$OMP                      pktop, zvir, mytile, dnats, nq_adv) &
!$OMP             private (dm, nb, blen, i, j, ix, k1, rTv, qgrs_rad)

   do nb = 1,Atm_block%nblks
! gas_phase_mass <-- prsl
! log(pe) <-- prsik

     blen = Atm_block%blksz(nb)

     !-- level interface geopotential height (relative to the surface)
     IPD_Data(nb)%Statein%phii(:,1) = 0.0_kind_phys
     IPD_Data(nb)%Statein%prsik(:,:) = 1.e25_kind_phys

     do k = 1, npz
       do ix = 1, blen
         i = Atm_block%index(nb)%ii(ix)
         j = Atm_block%index(nb)%jj(ix)

            !Indices for FV's vertical coordinate, for which 1 = top
            !here, k is the index for GFS's vertical coordinate, for which 1 = bottom
         k1 = npz+1-k ! flipping the index
         IPD_Data(nb)%Statein%tgrs(ix,k) = _DBL_(_RL_(Atm(mytile)%pt(i,j,k1)))
         IPD_Data(nb)%Statein%ugrs(ix,k) = _DBL_(_RL_(Atm(mytile)%ua(i,j,k1)))
         IPD_Data(nb)%Statein%vgrs(ix,k) = _DBL_(_RL_(Atm(mytile)%va(i,j,k1)))
         IPD_Data(nb)%Statein%vvl(ix,k)  = _DBL_(_RL_(Atm(mytile)%omga(i,j,k1)))
         IPD_Data(nb)%Statein%prsl(ix,k) = _DBL_(_RL_(Atm(mytile)%delp(i,j,k1)))   ! Total mass
         if (Atm(mytile)%flagstruct%do_skeb)IPD_Data(nb)%Statein%diss_est(ix,k) = _DBL_(_RL_(Atm(mytile)%diss_est(i,j,k1)))

         if (.not.Atm(mytile)%flagstruct%hydrostatic .and. (.not.Atm(mytile)%flagstruct%use_hydro_pressure))  &
           IPD_Data(nb)%Statein%phii(ix,k+1) = IPD_Data(nb)%Statein%phii(ix,k) - _DBL_(_RL_(Atm(mytile)%delz(i,j,k1)*grav))

! Convert to tracer mass:
         IPD_Data(nb)%Statein%qgrs(ix,k,1:nq_adv) =  _DBL_(_RL_(Atm(mytile)%q(i,j,k1,1:nq_adv))) &
                                                          * IPD_Data(nb)%Statein%prsl(ix,k)
         if (dnats .gt. 0) &
             IPD_Data(nb)%Statein%qgrs(ix,k,nq_adv+1:nq) =  _DBL_(_RL_(Atm(mytile)%q(i,j,k1,nq_adv+1:nq)))
         !--- SHOULD THESE BE CONVERTED TO MASS SINCE THE DYCORE DOES NOT TOUCH THEM IN ANY WAY???
         !--- See Note in state update...
         if ( ncnst > nq) &
             IPD_Data(nb)%Statein%qgrs(ix,k,nq+1:ncnst) = _DBL_(_RL_(Atm(mytile)%qdiag(i,j,k1,nq+1:ncnst)))
! Remove the contribution of condensates to delp (mass):
         if ( Atm(mytile)%flagstruct%nwat .eq. 6 ) then
            IPD_Data(nb)%Statein%prsl(ix,k) = IPD_Data(nb)%Statein%prsl(ix,k) &
                                            - IPD_Data(nb)%Statein%qgrs(ix,k,liq_wat)   &
                                            - IPD_Data(nb)%Statein%qgrs(ix,k,ice_wat)   &
                                            - IPD_Data(nb)%Statein%qgrs(ix,k,rainwat)   &
                                            - IPD_Data(nb)%Statein%qgrs(ix,k,snowwat)   &
                                            - IPD_Data(nb)%Statein%qgrs(ix,k,graupel)
         else !variable condensate numbers
            IPD_Data(nb)%Statein%prsl(ix,k) = IPD_Data(nb)%Statein%prsl(ix,k) &
                                            - sum(IPD_Data(nb)%Statein%qgrs(ix,k,2:Atm(mytile)%flagstruct%nwat))   
         endif
       enddo
     enddo

! Re-compute pressure (dry_mass + water_vapor) derived fields:
     do i=1,blen
        IPD_Data(nb)%Statein%prsi(i,npz+1) = ptop 
     enddo
     do k=npz,1,-1
        do i=1,blen
           IPD_Data(nb)%Statein%prsi(i,k)  = IPD_Data(nb)%Statein%prsi(i,k+1)  &
                                           + IPD_Data(nb)%Statein%prsl(i,k)
           IPD_Data(nb)%Statein%prsik(i,k) = log( IPD_Data(nb)%Statein%prsi(i,k) )
! Redefine mixing ratios for GFS == tracer_mass / (dry_air_mass + water_vapor_mass)
           IPD_Data(nb)%Statein%qgrs(i,k,1:nq_adv) = IPD_Data(nb)%Statein%qgrs(i,k,1:nq_adv) &
                                                   / IPD_Data(nb)%Statein%prsl(i,k)
        enddo
     enddo
     do i=1,blen
        IPD_Data(nb)%Statein%pgr(i)         = IPD_Data(nb)%Statein%prsi(i,1)    ! surface pressure for GFS
        IPD_Data(nb)%Statein%prsik(i,npz+1) = log(ptop)
     enddo

     do k=1,npz
        do i=1,blen
! Geo-potential at interfaces:
           qgrs_rad = max(qmin,IPD_Data(nb)%Statein%qgrs(i,k,sphum))
           rTv = rdgas*IPD_Data(nb)%Statein%tgrs(i,k)*(1.+zvir*qgrs_rad)
           if ( Atm(mytile)%flagstruct%hydrostatic .or. Atm(mytile)%flagstruct%use_hydro_pressure )   &
                IPD_Data(nb)%Statein%phii(i,k+1) = IPD_Data(nb)%Statein%phii(i,k) &
                                                     + rTv*(IPD_Data(nb)%Statein%prsik(i,k) &
                                                          - IPD_Data(nb)%Statein%prsik(i,k+1))
! Layer mean pressure by perfect gas law:
           dm = IPD_Data(nb)%Statein%prsl(i,k)
           IPD_Data(nb)%Statein%prsl(i,k) = dm*rTv/(IPD_Data(nb)%Statein%phii(i,k+1) &
                                                  - IPD_Data(nb)%Statein%phii(i,k))

!!! Ensure subgrid MONOTONICITY of Pressure: SJL 09/11/2016
           if ( .not.Atm(mytile)%flagstruct%hydrostatic ) then
! If violated, replaces it with hydrostatic pressure
              IPD_Data(nb)%Statein%prsl(i,k) = min(IPD_Data(nb)%Statein%prsl(i,k), &
                                                   IPD_Data(nb)%Statein%prsi(i,k)   - 0.01*dm)
              IPD_Data(nb)%Statein%prsl(i,k) = max(IPD_Data(nb)%Statein%prsl(i,k), &
                                                   IPD_Data(nb)%Statein%prsi(i,k+1) + 0.01*dm)
           endif
        enddo
     enddo

     do k = 1,npz
        do i=1,blen
! Exner function layer center: large sensitivity to non-hydro runs with moist kappa
           IPD_Data(nb)%Statein%prslk(i,k) = exp( kappa*log(IPD_Data(nb)%Statein%prsl(i,k)/p00) )
!--  layer center geopotential; geometric midpoint
           IPD_Data(nb)%Statein%phil(i,k) = 0.5_kind_phys*(IPD_Data(nb)%Statein%phii(i,k) &
                                                         + IPD_Data(nb)%Statein%phii(i,k+1))
        enddo
     enddo

! Compute Exner function at layer "interfaces"
    do i=1,blen
!      Top & Bottom edges computed hydrostaticlly
       IPD_Data(nb)%Statein%prsik(i,    1) = exp( kappa*IPD_Data(nb)%Statein%prsik(i,1) )*pk0inv  ! bottom
       IPD_Data(nb)%Statein%prsik(i,npz+1) = pktop                           ! TOA
    enddo

    if ( Atm(mytile)%flagstruct%hydrostatic .or. Atm(mytile)%flagstruct%use_hydro_pressure ) then
        do k=2,npz
           do i=1,blen
              IPD_Data(nb)%Statein%prsik(i,k) = exp( kappa*IPD_Data(nb)%Statein%prsik(i,k) )*pk0inv 
           enddo
        enddo
    endif
  enddo

 end subroutine atmos_phys_driver_statein

end module atmosphere_mod
