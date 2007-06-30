module atmosphere_mod

!-----------------------------------------------------------------------
!
! Interface for Cubed_Sphere fv dynamical core and physics
!
!-----------------------------------------------------------------------

!-----------------
! FMS modules:
!-----------------
use constants_mod,    only: cp_air, rdgas, grav, rvgas, kappa
use time_manager_mod, only: time_type, get_time, set_time, operator(+)
use fms_mod,          only: file_exist, open_namelist_file,    &
                            close_file, error_mesg, FATAL,     &
                            check_nml_error, stdlog,           &
                            write_version_number,              &
                            mpp_pe, mpp_root_pe, set_domain,   &
                            mpp_clock_id, mpp_clock_begin,     &
                            mpp_clock_end, CLOCK_SUBCOMPONENT, &
                            clock_flag_default, nullify_domain
use mpp_domains_mod,  only: domain2d
!-----------------
! FV core modules:
!-----------------
use grid_tools,         only: area
use fv_arrays_mod,      only: fv_atmos_type
use fv_pack_mod,        only: fv_init, domain, fv_end
use fv_dynamics_mod,    only: fv_dynamics
use fv_diagnostics_mod, only: fv_diag_init, fv_diag, fv_time
use fv_restart_mod,     only: fv_restart
use timingModule,       only: timing_on, timing_off
use fv_physics_mod,     only: fv_physics_down, fv_physics_up,  &
                           fv_physics_init, fv_physics_end, &
                           surf_diff_type

implicit none
private

public  atmosphere_down,       atmosphere_up,       &
        atmosphere_init,       atmosphere_end,      &
        atmosphere_resolution, atmosphere_boundary, &
        get_atmosphere_axes,   atmosphere_domain,   &
        get_bottom_mass,       get_bottom_wind,     &
        get_stock_pe,           surf_diff_type

!-----------------------------------------------------------------------

character(len=128) :: version = '$Id: atmosphere.F90,v 14.0.2.1 2007/05/23 15:46:07 z1l Exp $'
character(len=128) :: tagname = '$Name: nalanda_2007_06 $'

!---- namelist (saved in file input.nml) ----
!
! physics_window  The number of "i" by "j" rows processed each time
!                 the modular physics is called. To process the entire
!                 domain use physics_window = (/0,0/).
!                   [integer, default: physics_window = 0,0]

   integer, dimension(2) :: physics_window = (/0,0/)
   namelist /atmosphere_nml/ physics_window

!---- private data ----
  type (time_type) :: Time_step_atmos
  type (fv_atmos_type), allocatable :: Atm(:)
  real    :: dt_atmos
  real    :: zvir
  integer :: npx, npy, npz, ncnst
  integer :: isc, iec, jsc, jec
  integer :: nq                       ! transported tracers
  integer :: sec, seconds, days
  integer :: atmos_axes(4)
  integer :: ntiles=1
  integer :: id_dynam, id_phys_down, id_phys_up, id_fv_diag
  logical :: cold_start = .false.       ! read in initial condition

contains

 subroutine atmosphere_init (Time_init, Time, Time_step, Surf_diff)

 type (time_type),     intent(in)    :: Time_init, Time, Time_step
 type(surf_diff_type), intent(inout) :: Surf_diff

  integer :: unit, ierr, io
  integer :: ss, ds

  zvir = rvgas/rdgas - 1.

!----- read namelist -----
    if ( file_exist('input.nml') ) then
        unit = open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
           read (unit, nml=atmosphere_nml, iostat=io, end=10)
           ierr = check_nml_error (io, 'atmosphere_nml')
        enddo
 10     call close_file (unit)
    endif

!----- write version and namelist to log file -----

    call write_version_number ( version, tagname )
    if ( mpp_pe() == mpp_root_pe() ) write (stdlog(), nml=atmosphere_nml)

!---- compute physics/atmos time step in seconds ----

   Time_step_atmos = Time_step
   call get_time (Time_step_atmos, sec)
   dt_atmos = real(sec)

!----- initialize FV dynamical core -----
   cold_start = (.not.file_exist('INPUT/fv_core.res.nc'))

   allocate( Atm(ntiles) )

   call fv_init( Atm(:), dt_atmos )  ! allocates Atm components

   npx   = Atm(1)%npx
   npy   = Atm(1)%npy
   npz   = Atm(1)%npz
   ncnst = Atm(1)%ncnst

   isc = Atm(1)%isc
   iec = Atm(1)%iec
   jsc = Atm(1)%jsc
   jec = Atm(1)%jec

   nq = ncnst

   call fv_restart(domain, Atm, dt_atmos, seconds, days, cold_start)

   fv_time = Time

!----- initialize atmos_axes and fv_dynamics diagnostics

   call fv_diag_init(Atm, atmos_axes, Time, npx, npy, npz)

!----- initialize physics interface -----
!----- initialize domains for reading global physics data -----

   call set_domain ( domain )

   call fv_physics_init (Atm, atmos_axes, Time, physics_window, Surf_diff)

   call nullify_domain ( )

!  --- initialize clocks for dynamics, physics_down and physics_up
   id_dynam     = mpp_clock_id ('FV dynamical core',   &
          flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
   id_phys_down = mpp_clock_id ('Physics_down',   &
          flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
   id_phys_up   = mpp_clock_id ('Physics_up',   &
          flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
   id_fv_diag   = mpp_clock_id ('FV Diag',   &
          flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )

 end subroutine atmosphere_init




 subroutine atmosphere_down ( Time,    frac_land,          &
                           t_surf,  albedo,                &
                           albedo_vis_dir, albedo_nir_dir, &
                           albedo_vis_dif, albedo_nir_dif, &
                           rough_mom,                      &
                           u_star,  b_star, q_star,        &
                           dtau_du, dtau_dv, tau_x, tau_y, &
                           gust, coszen, flux_sw,          &
                           flux_sw_dir, flux_sw_dif,       &
                           flux_sw_down_vis_dir,           &          
                           flux_sw_down_vis_dif,           &
                           flux_sw_down_total_dir,         &
                           flux_sw_down_total_dif,         &
                           flux_sw_vis,                    &
                           flux_sw_vis_dir,                &
                           flux_sw_vis_dif,                &
                           flux_lw,                        &
                           Surf_diff                       )
!
!        Time = time at the current time level
!
   type(time_type),intent(in)      :: Time
   real, intent(in), dimension(:,:):: frac_land, t_surf, albedo,      &
                                      albedo_vis_dir, albedo_nir_dir, &
                                      albedo_vis_dif, albedo_nir_dif, &
                                      rough_mom, u_star, b_star,      &
                                      q_star, dtau_du, dtau_dv
   real, intent(inout), dimension(:,:):: tau_x,  tau_y
   real, intent(out),   dimension(:,:):: gust, coszen, flux_sw,    &
                                         flux_sw_dir, flux_sw_dif, &
                                         flux_sw_down_vis_dir,     &
                                         flux_sw_down_total_dir,   &
                                         flux_sw_down_vis_dif,     &
                                         flux_sw_down_total_dif,   &
                                         flux_sw_vis,              &
                                         flux_sw_vis_dir,          &
                                         flux_sw_vis_dif, flux_lw
   type(surf_diff_type), intent(inout):: Surf_diff
   type(time_type) :: Time_prev, Time_next

   Time_prev = Time                       ! two time-level scheme
   Time_next = Time + Time_step_atmos

!---- Call FV dynamics -----

   call mpp_clock_begin (id_dynam)
                    call timing_on('fv_dynamics')

   call fv_dynamics(npx, npy, npz, nq, Atm(1)%ng, dt_atmos, Atm(1)%consv_te,         &
                    Atm(1)%fill,  Atm(1)%reproduce_sum, kappa, cp_air, zvir,         &
                    Atm(1)%ks,    ncnst,        Atm(1)%n_split,   Atm(1)%q_split,    &
                    Atm(1)%u, Atm(1)%v, Atm(1)%w, Atm(1)%delz, Atm(1)%hydrostatic,   & 
                    Atm(1)%pt, Atm(1)%delp, Atm(1)%q, Atm(1)%ps, &
                    Atm(1)%pe, Atm(1)%pk, Atm(1)%peln, Atm(1)%pkz, Atm(1)%phis,      &
                    Atm(1)%omga, Atm(1)%ua, Atm(1)%va, Atm(1)%uc, Atm(1)%vc,         &
                    Atm(1)%ak, Atm(1)%bk, Atm(1)%mfx, Atm(1)%mfy,                    &
                    Atm(1)%cx, Atm(1)%cy, Atm(1)%u_srf, Atm(1)%v_srf,                &
                    Atm(1)%srf_init, Atm(1)%ze0, Atm(1)%hybrid_z)

                    call timing_off('fv_dynamics')
   call mpp_clock_end (id_dynam)


   call set_domain ( domain )
   call mpp_clock_begin (id_phys_down)
                         call timing_on('fv_physics_down')
   call fv_physics_down (Atm, dt_atmos, Time_prev, Time, Time_next,     &
                         frac_land, albedo,              &
                         albedo_vis_dir, albedo_nir_dir, &
                         albedo_vis_dif, albedo_nir_dif, &
                         rough_mom,  t_surf,             &
                         u_star,  b_star, q_star,        &
                         dtau_du, dtau_dv, tau_x, tau_y, &
                         flux_sw, flux_sw_dir,           &
                         flux_sw_dif,                    &
                         flux_sw_down_vis_dir,           &
                         flux_sw_down_vis_dif,           &
                         flux_sw_down_total_dir,         &
                         flux_sw_down_total_dif,         &
                         flux_sw_vis, flux_sw_vis_dir,   &
                         flux_sw_vis_dif, flux_lw,       &
                         coszen, gust, Surf_diff )
                         call timing_off('fv_physics_down')
   call mpp_clock_end (id_phys_down)
   call nullify_domain ( )

 end subroutine atmosphere_down



 subroutine atmosphere_up ( Time,  frac_land, Surf_diff, lprec, fprec, gust, &
                            u_star, b_star, q_star )

   type(time_type),intent(in)         :: Time
   type(surf_diff_type), intent(inout):: Surf_diff
   real, intent(in),  dimension(:,:)  :: frac_land
   real, intent(inout), dimension(:,:):: gust
   real, intent(out), dimension(:,:)  :: lprec,   fprec
   real, intent(in), dimension(:,:)   :: u_star, b_star, q_star

   type(time_type) :: Time_prev, Time_next

   Time_prev = Time                       ! two time-level scheme
   Time_next = Time + Time_step_atmos

   call set_domain ( domain )
   call mpp_clock_begin (id_phys_up)
!-----------------------------------------------------------------------
                       call timing_on('fv_physics_up')
   call fv_physics_up( Atm, dt_atmos, Time_prev, Time, Time_next,      &
                       frac_land, Surf_diff, lprec, fprec, gust ,      &
                       u_star, b_star, q_star   )
                       call timing_off('fv_physics_up')
!-----------------------------------------------------------------------
   call mpp_clock_end (id_phys_up)

   call mpp_clock_begin(id_fv_diag)

   fv_time = Time_next
   call get_time (fv_time, seconds,  days)
!-----------------------------------------------------------------------
                call timing_on('FV_DIAG')
   call fv_diag( Atm, zvir, fv_time, Atm(1)%print_freq )
                call timing_off('FV_DIAG')     
!-----------------------------------------------------------------------
   call mpp_clock_end(id_fv_diag)
   call nullify_domain ( )


 end subroutine atmosphere_up



 subroutine atmosphere_end (Time)
 type (time_type), intent(in) :: Time

  ! initialize domains for writing global physics data
    call set_domain ( domain )

    call get_time (Time, seconds,  days)
    call fv_physics_end(Time)
    call nullify_domain ( )
    call fv_end(Atm)
    deallocate (Atm)

 end subroutine atmosphere_end



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



 subroutine atmosphere_boundary (blon, blat, global)
!---------------------------------------------------------------
!    returns the longitude and latitude grid box edges
!    for either the local PEs grid (default) or the global grid
!---------------------------------------------------------------
    real,    intent(out) :: blon(:,:), blat(:,:)   ! Unit: radian
    logical, intent(in), optional :: global
! Local data:
    logical :: local
    integer i,j

    local = .TRUE.
    if( PRESENT(global) ) local = .NOT.global

    if (local) then
        do j=jsc,jec+1
           do i=isc,iec+1
              blon(i-isc+1,j-jsc+1) = Atm(1)%grid(i,j,1)
              blat(i-isc+1,j-jsc+1) = Atm(1)%grid(i,j,2)
           enddo
        end do
    else
        do j=1,npy
           do i=1,npx
              blon(i,j) = Atm(1)%grid_g(i,j,1)
              blat(i,j) = Atm(1)%grid_g(i,j,2)
           enddo
        end do
    end if

 end subroutine atmosphere_boundary



 subroutine atmosphere_domain ( fv_domain )
 type(domain2d), intent(out) :: fv_domain

!  returns the domain2d variable associated with the coupling grid
!  note: coupling is done using the mass/temperature grid with no halos

   fv_domain = domain

 end subroutine atmosphere_domain



 subroutine get_atmosphere_axes ( axes )
   integer, intent(out) :: axes (:)

!----- returns the axis indices for the atmospheric (mass) grid -----
     if ( size(axes(:)) < 0 .or. size(axes(:)) > 4 ) call error_mesg (    &
                                 'get_atmosphere_axes in atmosphere_mod', &
                                 'size of argument is incorrect', FATAL   )

     axes (1:size(axes(:))) = atmos_axes (1:size(axes(:)))

 
 end subroutine get_atmosphere_axes




 subroutine get_bottom_mass ( t_bot, tr_bot, p_bot, z_bot, p_surf )
!--------------------------------------------------------------
! returns temp, sphum, pres, height at the lowest model level
! and surface pressure
!--------------------------------------------------------------
   real, intent(out), dimension(isc:iec,jsc:jec):: t_bot, p_bot, z_bot, p_surf
   real, intent(out), dimension(isc:iec,jsc:jec,nq):: tr_bot
   integer :: i, j, m
   real rrg

   rrg  = rdgas / grav

     do j=jsc,jec
        do i=isc,iec
           p_surf(i,j) = Atm(1)%ps(i,j)
            t_bot(i,j) = Atm(1)%pt(i,j,npz)
            p_bot(i,j) = Atm(1)%delp(i,j,npz)/(Atm(1)%peln(i,npz+1,j)-Atm(1)%peln(i,npz,j))
            z_bot(i,j) = rrg*t_bot(i,j)*(1.+zvir*Atm(1)%q(i,j,npz,1)) *  &
                        (1. - Atm(1)%pe(i,npz,j)/p_bot(i,j))
        enddo
     enddo

! Copy tracers
     do m=1,nq
        do j=jsc,jec
           do i=isc,iec
              tr_bot(i,j,m) = Atm(1)%q(i,j,npz,m)
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
         u_bot(i,j) = Atm(1)%u_srf(i,j)
         v_bot(i,j) = Atm(1)%v_srf(i,j)
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
              wm(i,j) = wm(i,j) + Atm(1)%delp(i,j,k) * ( Atm(1)%q(i,j,k,1) +    &
                                                         Atm(1)%q(i,j,k,2) +    &
                                                         Atm(1)%q(i,j,k,3) )
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

end module atmosphere_mod
