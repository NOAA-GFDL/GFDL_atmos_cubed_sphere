module fv_nwp_nudge_mod

 use constants_mod,     only: pi, grav, rdgas, cp_air, kappa, radius
 use fms_mod,           only: write_version_number, open_namelist_file, &
                              check_nml_error, file_exist, close_file, read_data
 use fms_io_mod,        only: field_size
 use mpp_mod,           only: mpp_error, FATAL, stdlog, get_unit, input_nml_file
 use mpp_domains_mod,   only: mpp_update_domains
 use time_manager_mod,  only: time_type,  get_time, get_date

 use external_sst_mod,   only: i_sst, j_sst, sst_ncep, sst_anom, forecast_mode
 use fv_control_mod,    only: npx, npy
 use fv_grid_utils_mod, only: vlon, vlat, sina_u, sina_v, da_min, great_circle_dist, ks, intp_great_circle
 use fv_grid_utils_mod, only: latlon2xyz, vect_cross, normalize_vect
 use fv_grid_tools_mod, only: agrid, dx, dy, rdxc, rdyc, rarea, area
 use fv_diagnostics_mod,only: prt_maxmin, fv_time
 use tp_core_mod,       only: copy_corners
 use fv_mapz_mod,       only: mappm
 use fv_mp_mod,         only: is,js,ie,je, isd,jsd,ied,jed, gid, masterproc, domain,    &
                              mp_reduce_sum, mp_reduce_min, mp_reduce_max
 use fv_timing_mod,     only: timing_on, timing_off

!----------------------------------------------------------------------------------------
! Note: using FMS read_data for reading multiple analysis files requires too much memory
! and too many files open (led to fatal errors); 
!----------------------------------------------------------------------------------------
 use sim_nc_mod,        only: open_ncfile, close_ncfile, get_ncdim1, get_var1_double, get_var2_double,   &
                              get_var3_double

 implicit none
 private

 character(len=128) :: version = ''
 character(len=128) :: tagname = ''

 public fv_nwp_nudge, fv_nwp_nudge_init, fv_nwp_nudge_end, breed_slp_inline

 integer im     ! Data x-dimension
 integer jm     ! Data y-dimension
 integer km     ! Data z-dimension
 real, allocatable:: ak0(:), bk0(:)
 real, allocatable:: lat(:), lon(:)

 logical :: module_is_initialized = .false.
 logical :: master
 logical :: no_obs
 real :: deg2rad, rad2deg
 real :: time_nudge = 0.
 integer :: time_interval = 6*3600   ! dataset time interval (seconds)
 integer, parameter :: nfile_max = 125
 integer :: nfile

 integer :: k_breed = 0
 integer :: k_trop = 0
 real    :: p_trop = 200.E2

 real,    allocatable:: s2c(:,:,:)
 integer, allocatable:: id1(:,:), id2(:,:), jdc(:,:)
 real, allocatable :: ps_dat(:,:,:)
 real(KIND=4), allocatable, dimension(:,:,:,:):: u_dat, v_dat, t_dat, q_dat
 real, allocatable:: gz0(:,:)

! Namelist variables:
 character(len=128):: file_names(nfile_max)
 character(len=128):: track_file_name
 integer :: nfile_total = 0       ! =5 for 1-day (if datasets are 6-hr apart)
 real    :: p_wvp = 100.E2        ! cutoff level for specific humidity nudging 
 integer :: kord_data = 8

 logical :: add_bg_wind = .true.
 logical :: pre_test = .false.
 logical :: conserve_mom = .true.
 logical :: tc_mask = .false.
 logical :: strong_mask = .true. 
 logical :: ibtrack = .true. 
 logical :: nudge_debug = .false.
 logical :: nudge_t     = .false.
 logical :: nudge_q     = .false.
 logical :: nudge_winds = .true.
 logical :: nudge_virt  = .false.
 logical :: nudge_hght  = .false.
 logical :: nudge_tpw   = .false.   ! nudge total precipitable water
 logical :: time_varying = .true.
 logical :: time_track   = .false.
 logical :: print_end_breed = .true.
 logical :: print_end_nudge = .true.


! Nudging time-scales (seconds): note, however, the effective time-scale is 2X smaller (stronger) due
! to the use of the time-varying weighting factor
 real :: tau_q      = 86400.       ! 1-day
 real :: tau_tpw    = 86400.       ! 1-day
 real :: tau_winds  = 21600.       !  6-hr
 real :: tau_t      = 86400.
 real :: tau_virt   = 86400. 
 real :: tau_hght   = 86400.

 real :: q_min      = 1.E-8

 integer :: nf_uv = 0 
 integer :: nf_t  = 1 

! starting layer (top layer is sponge layer and is skipped)
 integer :: kstart = 2 

! skip "kbot" layers
 integer :: kbot_winds = 0 
 integer :: kbot_t     = 0 
 integer :: kbot_q     = 1 
 logical :: analysis_time

!-- Tropical cyclones  --------------------------------------------------------------------

! track dataset: 'INPUT/tropical_cyclones.txt'

  logical :: breed_vortex = .false.
  real :: grid_size     = 28.E3
  real :: tau_vt_slp    = 1800.
  real :: tau_vt_wind   = 1200.
  real :: tau_vt_rad    = 4.0  ! TEST.exe: 2.0

  real ::  slp_env = 101010.    ! storm environment pressure (pa)
  real :: pre0_env = 100000.    ! critical storm environment pressure (pa) for size computation
!------------------
  real:: r_lo = 2.0
  real:: r_hi = 5.0    ! try 4.0?
!------------------
  real::  r_fac = 1.2
  real :: r_min = 200.E3
  real :: r_inc =  25.E3
  real, parameter:: del_r = 50.E3
  real:: elapsed_time = 0.0
  real:: nudged_time = 1.E12 ! seconds 
                             ! usage example: set to 43200. to do inline vortex breeding
                             ! for only the first 12 hours
                             ! In addition, specify only 3 analysis files (12 hours)
  integer:: year_track_data
  integer, parameter:: max_storm = 140     ! max # of storms to process
  integer, parameter::  nobs_max = 125     ! Max # of observations per storm

  integer :: nstorms = 0
  integer :: nobs_tc(max_storm)
  real(KIND=4)::     x_obs(nobs_max,max_storm)           ! longitude in degrees
  real(KIND=4)::     y_obs(nobs_max,max_storm)           ! latitude in degrees
  real(KIND=4)::  wind_obs(nobs_max,max_storm)           ! observed 10-m wind speed (m/s)
  real(KIND=4)::  mslp_obs(nobs_max,max_storm)           ! observed SLP in mb
  real(KIND=4)::  mslp_out(nobs_max,max_storm)           ! outer ring SLP in mb
  real(KIND=4)::   rad_out(nobs_max,max_storm)           ! outer ring radius in meters
  real(KIND=4)::   time_tc(nobs_max,max_storm)           ! start time of the track
!------------------------------------------------------------------------------------------

 namelist /fv_nwp_nudge_nml/ nudge_virt, nudge_hght, nudge_t, nudge_q, nudge_winds, nudge_tpw, &
                          tau_winds, tau_t, tau_q, tau_virt, tau_hght,  kstart, kbot_winds, &
                          k_breed, k_trop, p_trop, kord_data, tc_mask, nudge_debug, nf_t,    &
                          nf_uv, breed_vortex, tau_vt_wind, tau_vt_slp, tau_tpw, strong_mask,   &
                          kbot_t, kbot_q, p_wvp, time_varying, time_track, time_interval,  &
                          pre0_env, tau_vt_rad, r_lo, r_hi, add_bg_wind, pre_test, conserve_mom,  &
                          nudged_time, r_fac, r_min, r_inc, ibtrack, track_file_name, file_names

 contains
 

  subroutine fv_nwp_nudge ( Time, dt, npz, ps_dt, u_dt, v_dt, t_dt, q_dt, zvir, &
                            ak, bk, ts, ps, delp, ua, va, pt, nwat, q, phis )

  type(time_type), intent(in):: Time
  integer,         intent(in):: npz           ! vertical dimension
  integer,         intent(in):: nwat
  real,            intent(in):: dt
  real,            intent(in):: zvir
  real, intent(in   ), dimension(npz+1):: ak, bk
  real, intent(in   ), dimension(isd:ied,jsd:jed    ):: phis
  real, intent(inout), dimension(isd:ied,jsd:jed,npz):: pt, ua, va, delp
  real, intent(inout):: q(isd:ied,jsd:jed,npz,nwat)
  real, intent(inout), dimension(isd:ied,jsd:jed):: ps
! Accumulated tendencies
  real, intent(inout), dimension(isd:ied,jsd:jed,npz):: u_dt, v_dt
  real, intent(out):: t_dt(is:ie,js:je,npz)
  real, intent(out):: q_dt(is:ie,js:je,npz)
  real, intent(out), dimension(is:ie,js:je):: ps_dt, ts
! local:
  real:: m_err(is:ie,js:je)         ! height error at specified model interface level
  real:: slp_n(is:ie,js:je)         ! "Observed" SLP
  real:: slp_m(is:ie,js:je)         ! "model" SLP
  real:: tpw_dat(is:ie,js:je)
  real:: tpw_mod(is:ie,js:je)
  real::   mask(is:ie,js:je)
  real:: gz_int(is:ie,js:je), gz(is:ie,npz+1), peln(is:ie,npz+1), pk(is:ie,npz+1)
  real:: pe1(is:ie)
  real:: pkz, ptmp
  real, allocatable :: ps_obs(:,:)
  real, allocatable, dimension(:,:,:):: u_obs, v_obs, t_obs, q_obs
  real, allocatable, dimension(:,:,:):: du_obs, dv_obs
  integer :: seconds, days
  integer :: i,j,k, iq, kht
  real :: factor, rms, bias, co
  real :: q_rat
  real :: dbk, rdt, press(npz), profile(npz), prof_t(npz), prof_q(npz), du, dv


  if ( .not. module_is_initialized ) then 
        call mpp_error(FATAL,'==> Error from fv_nwp_nudge: module not initialized')
  endif

  if ( no_obs ) then
       forecast_mode = .true.
       return
  endif

  call get_time (time, seconds, days)

  do j=js,je
     do i=is,ie
        mask(i,j) = 1.
     enddo
  enddo
  if ( tc_mask )  call get_tc_mask(time, mask)

! The following profile is suitable only for nwp purposes; if the analysis has a good representation
! of the strat-meso-sphere the profile for upper layers should be changed.

  profile(:) = 1.
  do k=1,npz
     press(k) = 0.5*(ak(k) + ak(k+1)) + 0.5*(bk(k)+bk(k+1))*1.E5
     if ( press(k) < 30.E2 ) then
          profile(k) =  max(0.01, press(k)/30.E2) 
     endif
  enddo
  profile(1) = 0.

! Thermodynamics:
  prof_t(:) = 1.
  do k=1,npz
     if ( press(k) < 30.E2 ) then
          prof_t(k) =  max(0.01, press(k)/30.E2) 
     endif
  enddo
  prof_t(1) = 0.
 
! Water vapor:
  prof_q(:) = 1.
  do k=1,npz
     if ( press(k) < 300.E2 ) then
          prof_q(k) =  max(0., press(k)/300.E2) 
     endif
  enddo
  prof_q(1) = 0.

! Height
  if ( k_trop == 0 ) then
       k_trop = 2
       do k=2,npz-1
          ptmp = ak(k+1) + bk(k+1)*1.E5
          if ( ptmp > p_trop ) then
               k_trop = k
               exit              
          endif
       enddo
  endif

  if ( time_varying ) then
       factor = 1. + cos(real(mod(seconds,time_interval))/real(time_interval)*2.*pi)
       factor = max(1.e-5, factor)
  else
       factor = 1.
  endif

  allocate (ps_obs(is:ie,js:je) )
  allocate ( t_obs(is:ie,js:je,npz) )
  allocate ( q_obs(is:ie,js:je,npz) )

  if ( nudge_winds ) then
       allocate (u_obs(is:ie,js:je,npz) )
       allocate (v_obs(is:ie,js:je,npz) )
  endif


  call get_obs(Time, dt, zvir, ak, bk, ps, ts, ps_obs, delp, u_obs, v_obs, t_obs, q_obs,   &
               tpw_dat, phis, gz_int, npz)

  if ( no_obs ) then
       deallocate (ps_obs)
       deallocate (t_obs)
       deallocate (q_obs)
       if ( nudge_winds ) then
            deallocate (u_obs)
            deallocate (v_obs)
       endif
       forecast_mode = .true.
       return
   endif

  if( analysis_time ) then
!     call prt_maxmin('PS_o', ps_obs, is, ie, js, je, 0, 1, 0.01, master)

! Compute RMSE, bias, and correlation of SLP 
      call compute_slp(is, ie, js, je,    pt(is:ie,js:je,npz:npz), ps(is:ie,js:je), phis(is:ie,js:je), slp_m)
      call compute_slp(is, ie, js, je, t_obs(is:ie,js:je,npz:npz), ps_obs, gz0, slp_n)

      call prt_maxmin('SLP_m', slp_m, is, ie, js, je, 0, 1, 0.01, master)
      call prt_maxmin('SLP_o', slp_n, is, ie, js, je, 0, 1, 0.01, master)

      do j=js,je
         do i=is,ie
            if ( phis(i,j)/grav > 500. ) then
! Exclude high terrains region for RMS and bias computation
                 m_err(i,j) = 0.
            else
                 m_err(i,j) = mask(i,j)*(slp_m(i,j) - slp_n(i,j))
            endif
         enddo
      enddo
     
      call rmse_bias(m_err, rms, bias)
      call corr(slp_m, slp_n, co)

      if(master) write(*,*) 'SLP (mb): RMS=', 0.01*rms, ' Bias=', 0.01*bias
      if(master) write(*,*) 'SLP correlation=',co
  endif

  ps_dt(:,:) = 0.

  if ( nudge_winds ) then

       allocate (du_obs(is:ie,js:je,npz) )
       allocate (dv_obs(is:ie,js:je,npz) )

! Compute tendencies:
     rdt = 1. / (tau_winds/factor + dt)
     do k=kstart, npz - kbot_winds
        do j=js,je
           do i=is,ie
              du_obs(i,j,k) = profile(k)*(u_obs(i,j,k)-ua(i,j,k))*rdt
              dv_obs(i,j,k) = profile(k)*(v_obs(i,j,k)-va(i,j,k))*rdt
           enddo
        enddo
     enddo

     if ( nf_uv>0 ) call del2_uv(du_obs, dv_obs, 0.20, npz, nf_uv)

     do k=kstart, npz - kbot_winds
        do j=js,je
           do i=is,ie
! Apply TC mask
              du_obs(i,j,k) = du_obs(i,j,k) * mask(i,j)
              dv_obs(i,j,k) = dv_obs(i,j,k) * mask(i,j)
!
              u_dt(i,j,k) = u_dt(i,j,k) + du_obs(i,j,k)
              v_dt(i,j,k) = v_dt(i,j,k) + dv_obs(i,j,k)
                ua(i,j,k) =   ua(i,j,k) + du_obs(i,j,k)*dt
                va(i,j,k) =   va(i,j,k) + dv_obs(i,j,k)*dt
           enddo
        enddo
     enddo

     deallocate (du_obs)
     deallocate (dv_obs)

  endif

  if ( nudge_t .or. nudge_virt ) then
     if(nudge_debug) call prt_maxmin('T_obs', t_obs, is, ie, js, je, 0, npz, 1., master)
  endif

!---------------------- temp -----------
  if ( nudge_virt .and. nudge_hght ) then
       tau_virt = max(tau_hght, tau_virt)
       kht = k_trop
  else
       kht = npz-kbot_t
  endif
!---------------------- temp -----------

  t_dt(:,:,:) = 0.

  if ( nudge_hght ) then
     if(nudge_debug) call prt_maxmin('H_int', gz_int, is, ie, js, je, 0, 1, 1./grav, master)

        rdt = dt / (tau_hght/factor + dt)

        do j=js,je
 
           do i=is,ie
              pe1(i) = ak(1)
              peln(i,1) = log(pe1(i))
                pk(i,1) = pe1(i)**kappa
           enddo
           do k=2, npz+1
              do i=is,ie
                    pe1(i) = pe1(i) + delp(i,j,k-1)
                 peln(i,k) = log(pe1(i))
                   pk(i,k) = pe1(i)**kappa
              enddo
           enddo

           do i=is,ie
              gz(i,npz+1) = phis(i,j)
           enddo
           do i=is,ie
              do k=npz,k_trop+1,-1
                 gz(i,k) = gz(i,k+1) + rdgas*pt(i,j,k)*(1.+zvir*q(i,j,k,1))*(peln(i,k+1)-peln(i,k))
              enddo
           enddo

           do i=is,ie
              m_err(i,j) = gz(i,k_trop+1)
           enddo

           do i=is,ie
              do k=k_trop+1,npz
! Add constant "virtual potential temperature" increment to correct height at p_interface
                       pkz = (pk(i,k+1)-pk(i,k))/(kappa*(peln(i,k+1)-peln(i,k)))
                 pt(i,j,k) = pt(i,j,k) + mask(i,j)*rdt*pkz*(gz_int(i,j)-gz(i,k_trop+1)) /     &
                            (cp_air*(1.+zvir*q(i,j,k,1))*(pk(i,npz+1)-pk(i,k_trop+1)))
              enddo
           enddo
        enddo   ! j-loop

! Compute RMSE of height
        if( analysis_time ) then
            call corr(m_err, gz_int, co)
            do j=js,je
               do i=is,ie
                  m_err(i,j) = (m_err(i,j) - gz_int(i,j)) / grav
               enddo
            enddo
            call rmse_bias(m_err, rms, bias)
            if(master) write(*,*) 'HGHT: RMSE (m)=', rms, ' Bias (m)=', bias
            if(master) write(*,*) 'HGHT: correlation=', co
        endif
  endif

  if ( nudge_t ) then
       rdt = 1./(tau_t/factor + dt)
     do k=kstart, kht
        do j=js,je
           do i=is,ie
              t_dt(i,j,k) = prof_t(k)*(t_obs(i,j,k)-pt(i,j,k))*rdt
           enddo
        enddo
     enddo
  elseif ( nudge_virt ) then
        rdt = 1./(tau_virt/factor + dt)
     do k=kstart, kht
        do j=js,je
           do i=is,ie
              t_dt(i,j,k) = prof_t(k)*(t_obs(i,j,k)/(1.+zvir*q(i,j,k,1))-pt(i,j,k))*rdt
           enddo
        enddo
     enddo
  endif

  deallocate ( t_obs )

  if ( nudge_t .or. nudge_virt ) then
! Filter t_dt here:
       if ( nf_t>0 ) call del2_scalar(t_dt, 0.20, npz, nf_t)

       do k=kstart, kht
          do j=js,je
             do i=is,ie
                pt(i,j,k) = pt(i,j,k) + t_dt(i,j,k)*dt*mask(i,j)
            enddo
         enddo
       enddo
  endif

  q_dt(:,:,:) = 0.
  if ( nudge_q ) then
       rdt = 1./(tau_q/factor + dt)
     do k=kstart, npz - kbot_q
        if ( press(k) > p_wvp ) then
            do iq=2,nwat
               do j=js,je
                  do i=is,ie
                     q(i,j,k,iq) = q(i,j,k,iq)*delp(i,j,k)
                  enddo
               enddo
            enddo
! Specific humidity:
            do j=js,je
               do i=is,ie
                  delp(i,j,k) = delp(i,j,k)*(1.-q(i,j,k,1))
                  q_dt(i,j,k) = prof_q(k)*(max(q_min,q_obs(i,j,k))-q(i,j,k,1))*rdt*mask(i,j)
                   q(i,j,k,1) = q(i,j,k,1) + q_dt(i,j,k)*dt
                  delp(i,j,k) = delp(i,j,k)/(1.-q(i,j,k,1))
               enddo
            enddo
            do iq=2,nwat
               do j=js,je
                  do i=is,ie
                     q(i,j,k,iq) = q(i,j,k,iq)/delp(i,j,k)
                  enddo
               enddo
            enddo
        endif
     enddo
  elseif ( nudge_tpw ) then
! Compute tpw_model
    tpw_mod(:,:) = 0.
    do k=1,npz
       do j=js,je
          do i=is,ie
             tpw_mod(i,j) = tpw_mod(i,j) + q(i,j,k,1)*delp(i,j,k)
          enddo
       enddo
    enddo

   do j=js,je
      do i=is,ie
         tpw_dat = tpw_dat(i,j) / max(tpw_mod(i,j), q_min)
      enddo
   enddo
   if(nudge_debug) call prt_maxmin('TPW_rat', tpw_dat, is, ie, js, je, 0, 1, 1., master)

   do k=1,npz

      do iq=2,nwat
         do j=js,je
            do i=is,ie
               q(i,j,k,iq) = q(i,j,k,iq)*delp(i,j,k)
            enddo
         enddo
      enddo

      do j=js,je
         do i=is,ie
            delp(i,j,k) = delp(i,j,k)*(1.-q(i,j,k,1))
                  q_rat = max(0.25,  min(tpw_dat(i,j), 4.))
             q(i,j,k,1) = q(i,j,k,1)*(tau_tpw/factor + mask(i,j)*dt*q_rat)/(tau_tpw/factor + mask(i,j)*dt)
            delp(i,j,k) = delp(i,j,k)/(1.-q(i,j,k,1))
         enddo
      enddo

      do iq=2,nwat
         do j=js,je
            do i=is,ie
               q(i,j,k,iq) = q(i,j,k,iq)/delp(i,j,k)
            enddo
         enddo
      enddo

   enddo   ! k-loop
  endif

  deallocate ( q_obs )
  deallocate ( ps_obs )



  if ( breed_vortex )   &
  call breed_srf_winds(Time, dt, npz, u_obs, v_obs, ak, bk, ps, phis, delp, ua, va, u_dt, v_dt, pt, q, nwat, zvir)

  if ( nudge_winds ) then
     deallocate ( u_obs )
     deallocate ( v_obs )
  endif

 end  subroutine fv_nwp_nudge


 subroutine compute_slp(isc, iec, jsc, jec, tm, ps, gz, slp)
 integer, intent(in):: isc, iec, jsc, jec
 real, intent(in), dimension(isc:iec,jsc:jec):: tm, ps, gz
 real, intent(out):: slp(isc:iec,jsc:jec)
 integer:: i,j

    do j=jsc,jec
       do i=isc,iec
          slp(i,j) = ps(i,j) * exp( gz(i,j)/(rdgas*(tm(i,j) + 3.25E-3*gz(i,j)/grav)) )
       enddo
    enddo

 end subroutine compute_slp


 subroutine get_obs(Time, dt, zvir, ak, bk, ps, ts, ps_obs, delp, u_obs, v_obs, t_obs, q_obs,  &
                    tpw_dat, phis, gz_int, npz)
  type(time_type), intent(in):: Time
  integer,         intent(in):: npz           ! vertical dimension
  real,            intent(in):: zvir
  real,            intent(in):: dt
  real, intent(in), dimension(npz+1):: ak, bk
  real, intent(in), dimension(isd:ied,jsd:jed):: phis
  real, intent(in), dimension(isd:ied,jsd:jed,npz):: delp
  real, intent(inout), dimension(isd:ied,jsd:jed):: ps
  real, intent(out), dimension(is:ie,js:je):: ts, ps_obs
  real, intent(out), dimension(is:ie,js:je,npz):: u_obs, v_obs, t_obs, q_obs
  real, intent(out)::  gz_int(is:ie,js:je)
  real, intent(out):: tpw_dat(is:ie,js:je)
! local:
  real(KIND=4), allocatable:: ut(:,:,:), vt(:,:,:)
  real, dimension(is:ie,js:je):: h1, h2
  integer :: seconds, days
  integer :: i,j,k
  real :: alpha, beta

  call get_time (time, seconds, days)

  seconds = seconds - nint(dt)

! Data must be "time_interval" (hr) apart; keep two time levels in memory

  no_obs = .false.
  analysis_time = .false.

  if ( mod(seconds, time_interval) == 0 ) then

    if ( nfile == nfile_total ) then
         no_obs = .true.
         forecast_mode = .true.
         if(print_end_nudge)  then
            print_end_nudge = .false.
            if (master) write(*,*) '*** L-S nudging Ended at', days, seconds
         endif
         return              ! free-running mode
    endif

      ps_dat(:,:,1) = ps_dat(:,:,2)
      if ( nudge_winds ) then
         u_dat(:,:,:,1) = u_dat(:,:,:,2)
         v_dat(:,:,:,1) = v_dat(:,:,:,2)
      endif
      t_dat(:,:,:,1) = t_dat(:,:,:,2)
      q_dat(:,:,:,1) = q_dat(:,:,:,2)

!---------------
! Read next data
!---------------
      nfile = nfile + 1
      call get_ncep_analysis ( ps_dat(:,:,2), u_dat(:,:,:,2), v_dat(:,:,:,2),    &
                              t_dat(:,:,:,2), q_dat(:,:,:,2), zvir,  &
                              ts, nfile, file_names(nfile) )
      analysis_time = .true.
      time_nudge = dt
  else
      time_nudge = time_nudge + dt
  endif

!--------------------
! Time interpolation:
!--------------------

  beta = time_nudge / real(time_interval)

  if ( beta < 0. .or. beta >  (1.+1.E-7) ) then
       call mpp_error(FATAL,'==> Error from get_obs:data out of range')
  endif

  alpha = 1. - beta

! Warning: ps_data are not adjusted for the differences in terrain yet
  ps_obs(:,:)  = alpha*ps_dat(:,:,1) + beta*ps_dat(:,:,2)

  allocate ( ut(is:ie,js:je,npz) )
  allocate ( vt(is:ie,js:je,npz) )

  if ( nudge_winds ) then

       call remap_uv(npz, ak,  bk, ps(is:ie,js:je), delp,  ut,     vt,   &
                     km, ps_dat(is:ie,js:je,1),  u_dat(:,:,:,1), v_dat(:,:,:,1) )

       u_obs(:,:,:) = alpha*ut(:,:,:)
       v_obs(:,:,:) = alpha*vt(:,:,:)

       call remap_uv(npz, ak, bk, ps(is:ie,js:je), delp,   ut,      vt,   &
                     km, ps_dat(is:ie,js:je,2),  u_dat(:,:,:,2), v_dat(:,:,:,2) )

       u_obs(:,:,:) = u_obs(:,:,:) + beta*ut(:,:,:)
       v_obs(:,:,:) = v_obs(:,:,:) + beta*vt(:,:,:)
  endif

! if ( nudge_t .or. nudge_virt .or. nudge_q .or. nudge_hght .or. nudge_tpw ) then

       call remap_tq(npz, ak, bk, ps(is:ie,js:je), delp,  ut,  vt,  &
                     km,  ps_dat(is:ie,js:je,1),  t_dat(:,:,:,1), q_dat(:,:,:,1), zvir)

       t_obs(:,:,:) = alpha*ut(:,:,:)
       q_obs(:,:,:) = alpha*vt(:,:,:)

       call remap_tq(npz, ak, bk, ps(is:ie,js:je), delp,  ut,  vt,  &
                     km,  ps_dat(is:ie,js:je,2),  t_dat(:,:,:,2), q_dat(:,:,:,2), zvir)

       t_obs(:,:,:) = t_obs(:,:,:) + beta*ut(:,:,:)
       q_obs(:,:,:) = q_obs(:,:,:) + beta*vt(:,:,:)

           do j=js,je
              do i=is,ie
                 tpw_dat(i,j) = 0.
              enddo
           enddo
       if ( nudge_tpw ) then
           do k=1,km
           do j=js,je
              do i=is,ie
                 tpw_dat(i,j) = tpw_dat(i,j) + q_obs(i,j,k) *     &
                              ( ak0(k+1)-ak0(k) + (bk0(k+1)-bk0(k))*ps_obs(i,j) )
              enddo
           enddo
           enddo
       endif
! endif

  if ( nudge_hght ) then
       call get_int_hght(h1, npz, ak, bk, ps(is:ie,js:je), delp, ps_dat(is:ie,js:je,1), t_dat(:,:,:,1))
!      if(nudge_debug) call prt_maxmin('H_1', h1, is, ie, js, je, 0, 1, 1./grav, master)

       call get_int_hght(h2, npz, ak, bk, ps(is:ie,js:je), delp, ps_dat(is:ie,js:je,2), t_dat(:,:,:,2))
!      if(nudge_debug) call prt_maxmin('H_2', h2, is, ie, js, je, 0, 1, 1./grav, master)

       gz_int(:,:) = alpha*h1(:,:) + beta*h2(:,:) 
  endif

  deallocate ( ut ) 
  deallocate ( vt ) 

 end subroutine get_obs


 subroutine fv_nwp_nudge_init(npz, zvir, ak, bk, ts, phis)
  integer,  intent(in):: npz           ! vertical dimension 
  real,     intent(in):: zvir
  real, intent(in), dimension(isd:ied,jsd:jed):: phis
  real, intent(in), dimension(npz+1):: ak, bk
  real, intent(out), dimension(is:ie,js:je):: ts
  logical found
  integer tsize(4)
  integer :: i, j, f_unit, unit, io, ierr, nt, k
  integer :: ncid

   master = gid==masterproc

   deg2rad = pi/180.
   rad2deg = 180./pi

   grid_size = 1.E7/real(npx-1)         ! mean grid size

   do nt=1,nfile_max
      file_names(nt) = "No_File_specified"
   enddo

   track_file_name = "No_File_specified"

#ifdef INTERNAL_FILE_NML
    read( input_nml_file, nml = fv_nwp_nudge_nml, iostat = io )
    ierr = check_nml_error(io,'fv_nwp_nudge_nml')
#else
    if( file_exist( 'input.nml' ) ) then
       unit = open_namelist_file ()
       io = 1
       do while ( io .ne. 0 )
          read( unit, nml = fv_nwp_nudge_nml, iostat = io, end = 10 )
          ierr = check_nml_error(io,'fv_nwp_nudge_nml')
       end do
10     call close_file ( unit )
    end if
#endif
    call write_version_number (version, tagname)
    if ( master ) then
         f_unit=stdlog()
         write( f_unit, nml = fv_nwp_nudge_nml )
         write(*,*) 'NWP nudging initialized.'
    endif

    if ( nudge_virt ) then
         nudge_t = .false.
!        nudge_q = .false.
    endif

    if ( nudge_t ) nudge_virt = .false.

    do nt=1,nfile_max
      if ( file_names(nt) == "No_File_specified" ) then
           nfile_total = nt - 1
           if(master) write(*,*) 'Total of NCEP files specified=', nfile_total
           exit
      endif
    enddo


! Initialize remapping coefficients:

!   call field_size(file_names(1), 'T', tsize, field_found=found)
!   if ( found ) then
!        im = tsize(1); jm = tsize(2); km = tsize(3)
!        if(master)  write(*,*) 'NCEP analysis dimensions:', tsize
!   else
!        call mpp_error(FATAL,'==> Error from get_ncep_analysis: T field not found')
!   endif
    call open_ncfile( file_names(1), ncid )        ! open the file
    call get_ncdim1( ncid, 'lon', im )
    call get_ncdim1( ncid, 'lat', jm )
    call get_ncdim1( ncid, 'lev', km )
    if(master)  write(*,*) 'NCEP analysis dimensions:', im, jm, km

    allocate ( s2c(is:ie,js:je,4) )
    allocate ( id1(is:ie,js:je) )
    allocate ( id2(is:ie,js:je) )
    allocate ( jdc(is:ie,js:je) )

    allocate (  lon(im) )
    allocate (  lat(jm) )

!   call read_data (file_names(1), 'LAT', lat, no_domain=.true.)
!   call read_data (file_names(1), 'LON', lon, no_domain=.true.)
    call get_var1_double (ncid, 'lon', im, lon )
    call get_var1_double (ncid, 'lat', jm, lat )

! Convert to radian
    do i=1,im
       lon(i) = lon(i) * deg2rad ! lon(1) = 0.
    enddo
    do j=1,jm
       lat(j) = lat(j) * deg2rad
    enddo
 
    allocate ( ak0(km+1) )
    allocate ( bk0(km+1) )

!   call read_data (file_names(1), 'hyai', ak0, no_domain=.true.)
!   call read_data (file_names(1), 'hybi', bk0, no_domain=.true.)
    call get_var1_double (ncid, 'hyai', km+1, ak0 )
    call get_var1_double (ncid, 'hybi', km+1, bk0 )
    call close_ncfile( ncid )

! Note: definition of NCEP hybrid is p(k) = a(k)*1.E5 + b(k)*ps
    ak0(:) = ak0(:) * 1.E5

! Limiter to prevent NAN at top during remapping 
    ak0(1) = max(1.e-8, ak0(1))

   if ( master ) then
      do k=1,npz
         write(*,*) k, 0.5*(ak(k)+ak(k+1))+0.5*(bk(k)+bk(k+1))*1.E5,  'del-B=', bk(k+1)-bk(k)
      enddo
   endif

   if ( k_breed==0 ) k_breed = ks
!  k_breed = ks

   call slp_obs_init

!-----------------------------------------------------------
! Initialize lat-lon to Cubed bi-linear interpolation coeff:
!-----------------------------------------------------------
    call remap_coef

    allocate ( gz0(is:ie,js:je) )
    allocate (ps_dat(is:ie,js:je,2) )
    allocate ( u_dat(is:ie,js:je,km,2) )
    allocate ( v_dat(is:ie,js:je,km,2) )
    allocate ( t_dat(is:ie,js:je,km,2) )
    allocate ( q_dat(is:ie,js:je,km,2) )


! Get first dataset
    nt = 2
    nfile = 1
    call get_ncep_analysis ( ps_dat(:,:,nt), u_dat(:,:,:,nt), v_dat(:,:,:,nt),     &
                            t_dat(:,:,:,nt), q_dat(:,:,:,nt), zvir,   &
                            ts, nfile, file_names(nfile) )


    module_is_initialized = .true.
    
 end subroutine fv_nwp_nudge_init


 subroutine get_ncep_analysis ( ps, u, v, t, q, zvir, ts, nfile, fname )
  real,     intent(in):: zvir
  character(len=128), intent(in):: fname
  integer,  intent(inout):: nfile
!
  real, intent(out), dimension(is:ie,js:je):: ts
  real, intent(out), dimension(is:ie,js:je):: ps
  real(KIND=4), intent(out), dimension(is:ie,js:je,km):: u, v, t, q
! local:
  real, allocatable:: oro(:,:), wk2(:,:), wk3(:,:,:)
  real tmean
  integer:: i, j, k, npt
  integer:: i1, i2, j1, ncid
  logical found
  logical:: read_ts = .true.
  logical:: land_ts = .false.

  if( .not. file_exist(fname) ) then
     call mpp_error(FATAL,'==> Error from get_ncep_analysis: file not found')
  else
     if(master) write(*,*) 'Reading NCEP anlysis file:', fname 
  endif

!----------------------------------
! remap surface pressure and height:
!----------------------------------
     allocate ( wk2(im,jm) )

!    call read_data (fname, 'PS', wk2, no_domain=.true.)
     call open_ncfile( fname, ncid )        ! open the file
     call get_var2_double( ncid, 'PS', im, jm, wk2 )

     if(gid==0) call pmaxmin( 'PS_ncep', wk2, im,  jm, 0.01)

     do j=js,je
        do i=is,ie
           i1 = id1(i,j)
           i2 = id2(i,j)
           j1 = jdc(i,j)
           ps(i,j) = s2c(i,j,1)*wk2(i1,j1  ) + s2c(i,j,2)*wk2(i2,j1  ) +  &
                     s2c(i,j,3)*wk2(i2,j1+1) + s2c(i,j,4)*wk2(i1,j1+1)
        enddo
     enddo

!    call read_data (fname, 'PHIS', wk2, no_domain=.true.)
     call get_var2_double( ncid, 'PHIS', im, jm, wk2 )

!    if(gid==0) call pmaxmin( 'ZS_ncep', wk2, im,  jm, 1./grav)
     do j=js,je
        do i=is,ie
           i1 = id1(i,j)
           i2 = id2(i,j)
           j1 = jdc(i,j)
           gz0(i,j) = s2c(i,j,1)*wk2(i1,j1  ) + s2c(i,j,2)*wk2(i2,j1  ) +  &
                      s2c(i,j,3)*wk2(i2,j1+1) + s2c(i,j,4)*wk2(i1,j1+1)
        enddo
     enddo
     call prt_maxmin('ZS_ncep', gz0, is, ie, js, je, 0, 1, 1./grav, master)

     if ( read_ts ) then       ! read skin temperature; could be used for SST

!     call read_data (fname, 'TS', wk2, no_domain=.true.)
      call get_var2_double( ncid, 'TS', im, jm, wk2 )

      if ( .not. land_ts ) then
           allocate ( oro(im,jm) )

! Read NCEP ORO (1; land; 0: ocean; 2: sea_ice)
!          call read_data (fname, 'ORO', oro, no_domain=.true.)
           call get_var2_double( ncid, 'ORO', im, jm, oro )

           do j=1,jm
              tmean = 0.
              npt = 0
              do i=1,im
                 if( abs(oro(i,j)-1.) > 0.5 ) then
                     tmean = tmean + wk2(i,j)
                     npt = npt + 1
                 endif
              enddo
!-------------------------------------------------------
! Replace TS over interior land with zonal mean SST/Ice 
!-------------------------------------------------------
              if ( npt /= 0 ) then
                   tmean= tmean / real(npt)
                   do i=1,im
                      if( abs(oro(i,j)-1.) <= 0.5 ) then
                          if ( i==1 ) then
                               i1 = im;     i2 = 2
                          elseif ( i==im ) then
                               i1 = im-1;   i2 = 1
                          else
                               i1 = i-1;    i2 = i+1
                          endif
                          if ( abs(oro(i2,j)-1.)>0.5 ) then     ! east side has priority
                               wk2(i,j) = wk2(i2,j)
                          elseif ( abs(oro(i1,j)-1.)>0.5 ) then ! west side
                               wk2(i,j) = wk2(i1,j)
                          else
                               wk2(i,j) = tmean
                          endif
                      endif
                   enddo
              endif
           enddo
           deallocate ( oro )
      endif   ! land_ts

      if(gid==0) call pmaxmin('SST_ncep', wk2, im,  jm, 1.)
      do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            ts(i,j) = s2c(i,j,1)*wk2(i1,j1  ) + s2c(i,j,2)*wk2(i2,j1  ) +  &
                      s2c(i,j,3)*wk2(i2,j1+1) + s2c(i,j,4)*wk2(i1,j1+1)
         enddo
      enddo
      call prt_maxmin('SST_model', ts, is, ie, js, je, 0, 1, 1., master)

! Perform interp to FMS SST format/grid
      call ncep2fms( wk2 )
      if(gid==0) call pmaxmin( 'SST_ncep', sst_ncep, i_sst, j_sst, 1.)
      if(gid==0) call pmaxmin( 'SST_anom', sst_anom, i_sst, j_sst, 1.)

      endif     ! read_ts

      deallocate ( wk2 ) 

! Read in temperature:
      allocate (  wk3(im,jm,km) )

! Winds:
   if ( nudge_winds ) then

!     call read_data (fname, 'U',  wk3, no_domain=.true.)
      call get_var3_double( ncid, 'U', im, jm, km , wk3 )
      if( master ) call pmaxmin( 'U_ncep',   wk3, im*jm, km, 1.)

      do k=1,km
      do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            u(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                       s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
         enddo
      enddo
      enddo

!     call read_data (fname, 'V',  wk3, no_domain=.true.)
      call get_var3_double( ncid, 'V', im, jm, km , wk3 )

      if( master ) call pmaxmin( 'V_ncep',  wk3, im*jm, km, 1.)
      do k=1,km
      do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            v(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                       s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
         enddo
      enddo
      enddo

   endif

!  if ( nudge_t .or. nudge_virt .or. nudge_q .or. nudge_tpw .or. nudge_hght ) then

! Read in tracers: only sphum at this point
!     call read_data (fname, 'Q', wk3, no_domain=.true.)
      call get_var3_double( ncid, 'Q', im, jm, km , wk3 )

      if(gid==1) call pmaxmin( 'Q_ncep',   wk3, im*jm, km, 1.)
      do k=1,km
      do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            q(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                       s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
         enddo
      enddo
      enddo

!     call read_data (fname, 'T',  wk3, no_domain=.true.)
      call get_var3_double( ncid, 'T', im, jm, km , wk3 )
      call close_ncfile ( ncid )

      if(gid==0) call pmaxmin( 'T_ncep',   wk3, im*jm, km, 1.)

      do k=1,km
      do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            t(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                       s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
! Convert t to virtual temperature:
            t(i,j,k) = t(i,j,k)*(1.+zvir*q(i,j,k))
         enddo
      enddo
      enddo

!  endif

   deallocate ( wk3 ) 

! nfile = nfile + 1

 end subroutine get_ncep_analysis



 subroutine remap_coef

! local:
  real :: rdlon(im)
  real :: rdlat(jm)
  real:: a1, b1
  integer i,j, i1, i2, jc, i0, j0

  do i=1,im-1
     rdlon(i) = 1. / (lon(i+1) - lon(i))
  enddo
     rdlon(im) = 1. / (lon(1) + 2.*pi - lon(im))

  do j=1,jm-1
     rdlat(j) = 1. / (lat(j+1) - lat(j))
  enddo

! * Interpolate to cubed sphere cell center
  do 5000 j=js,je

     do i=is,ie

       if ( agrid(i,j,1)>lon(im) ) then
            i1 = im;     i2 = 1
            a1 = (agrid(i,j,1)-lon(im)) * rdlon(im)
       elseif ( agrid(i,j,1)<lon(1) ) then
            i1 = im;     i2 = 1
            a1 = (agrid(i,j,1)+2.*pi-lon(im)) * rdlon(im)
       else
            do i0=1,im-1
            if ( agrid(i,j,1)>=lon(i0) .and. agrid(i,j,1)<=lon(i0+1) ) then
               i1 = i0;  i2 = i0+1
               a1 = (agrid(i,j,1)-lon(i1)) * rdlon(i0)
               go to 111
            endif
            enddo
       endif
111    continue

       if ( agrid(i,j,2)<lat(1) ) then
            jc = 1
            b1 = 0.
       elseif ( agrid(i,j,2)>lat(jm) ) then
            jc = jm-1
            b1 = 1.
       else
          do j0=1,jm-1
          if ( agrid(i,j,2)>=lat(j0) .and. agrid(i,j,2)<=lat(j0+1) ) then
               jc = j0
               b1 = (agrid(i,j,2)-lat(jc)) * rdlat(jc)
               go to 222
          endif
          enddo
       endif
222    continue

       if ( a1<0.0 .or. a1>1.0 .or.  b1<0.0 .or. b1>1.0 ) then
            write(*,*) 'gid=', gid, i,j,a1, b1
       endif

       s2c(i,j,1) = (1.-a1) * (1.-b1)
       s2c(i,j,2) =     a1  * (1.-b1)
       s2c(i,j,3) =     a1  *     b1
       s2c(i,j,4) = (1.-a1) *     b1
       id1(i,j) = i1
       id2(i,j) = i2
       jdc(i,j) = jc
     enddo   !i-loop
5000 continue   ! j-loop

 end subroutine remap_coef


 subroutine ncep2fms( sst )
  real, intent(in):: sst(im,jm)
! local:
  real :: rdlon(im)
  real :: rdlat(jm)
  real:: a1, b1
  real:: delx, dely
  real:: xc, yc    ! "data" location
  real:: c1, c2, c3, c4
  integer i,j, i1, i2, jc, i0, j0, it, jt

  do i=1,im-1
     rdlon(i) = 1. / (lon(i+1) - lon(i))
  enddo
     rdlon(im) = 1. / (lon(1) + 2.*pi - lon(im))

  do j=1,jm-1
     rdlat(j) = 1. / (lat(j+1) - lat(j))
  enddo

! * Interpolate to "FMS" 1x1 SST data grid
! lon: 0.5, 1.5, ..., 359.5
! lat: -89.5, -88.5, ... , 88.5, 89.5

  delx = 360./real(i_sst) 
  dely = 180./real(j_sst) 

  jt = 1
  do 5000 j=1,j_sst

     yc = (-90. + dely * (0.5+real(j-1)))  * deg2rad
     if ( yc<lat(1) ) then
            jc = 1
            b1 = 0.
     elseif ( yc>lat(jm) ) then
            jc = jm-1
            b1 = 1.
     else
          do j0=jt,jm-1
          if ( yc>=lat(j0) .and. yc<=lat(j0+1) ) then
               jc = j0
               jt = j0
               b1 = (yc-lat(jc)) * rdlat(jc)
               go to 222
          endif
          enddo
     endif
222  continue
     it = 1

     do i=1,i_sst
        xc = delx * (0.5+real(i-1)) * deg2rad
       if ( xc>lon(im) ) then
            i1 = im;     i2 = 1
            a1 = (xc-lon(im)) * rdlon(im)
       elseif ( xc<lon(1) ) then
            i1 = im;     i2 = 1
            a1 = (xc+2.*pi-lon(im)) * rdlon(im)
       else
            do i0=it,im-1
            if ( xc>=lon(i0) .and. xc<=lon(i0+1) ) then
               i1 = i0;  i2 = i0+1
               it = i0
               a1 = (xc-lon(i1)) * rdlon(i0)
               go to 111
            endif
            enddo
       endif
111    continue

!      if ( a1<0.0 .or. a1>1.0 .or.  b1<0.0 .or. b1>1.0 ) then
!           write(*,*) 'gid=', gid, i,j,a1, b1
!      endif
       c1 = (1.-a1) * (1.-b1)
       c2 =     a1  * (1.-b1)
       c3 =     a1  *     b1
       c4 = (1.-a1) *     b1
! Interpolated surface pressure
       sst_ncep(i,j) = c1*sst(i1,jc  ) + c2*sst(i2,jc  ) +    &
                       c3*sst(i2,jc+1) + c4*sst(i1,jc+1)
     enddo   !i-loop
5000 continue   ! j-loop

 end subroutine ncep2fms


 subroutine get_int_hght(h_int, npz, ak, bk, ps, delp, ps0, tv)
  integer, intent(in):: npz
  real,    intent(in):: ak(npz+1), bk(npz+1)
  real,    intent(in), dimension(is:ie,js:je):: ps, ps0
  real, intent(in), dimension(isd:ied,jsd:jed,npz):: delp
  real(KIND=4),  intent(in), dimension(is:ie,js:je,km):: tv
  real,   intent(out), dimension(is:ie,js:je):: h_int  ! g*height
! local:
  real, dimension(is:ie,km+1):: pn0, gz
  real:: logp(is:ie)
  integer i,j,k

  h_int(:,:) = 1.E25

  do 5000 j=js,je

     do k=1,km+1
        do i=is,ie
           pn0(i,k) = log( ak0(k) + bk0(k)*ps0(i,j) )
        enddo
     enddo 
!------
! Model
!------
     do i=is,ie
        logp(i) = ak(1)
     enddo
     do k=1,k_trop
       do i=is,ie
          logp(i) = logp(i) + delp(i,j,k)
       enddo
     enddo
     do i=is,ie
        logp(i) = log( logp(i) )
        gz(i,km+1) = gz0(i,j)   ! Data Surface geopotential
     enddo

! Linear in log-p interpolation
     do i=is,ie
        do k=km,1,-1
           gz(i,k) = gz(i,k+1) + rdgas*tv(i,j,k)*(pn0(i,k+1)-pn0(i,k))
           if ( logp(i)>=pn0(i,k) .and. logp(i)<=pn0(i,k+1) ) then
               h_int(i,j) = gz(i,k+1) + (gz(i,k)-gz(i,k+1))*(pn0(i,k+1)-logp(i))/(pn0(i,k+1)-pn0(i,k))
               goto 400
          endif
       enddo
400    continue
    enddo

5000 continue


 end subroutine get_int_hght



 subroutine remap_tq( npz, ak,  bk,  ps, delp,  t,  q,  &
                      kmd, ps0, ta, qa, zvir)
  integer, intent(in):: npz, kmd
  real,    intent(in):: zvir
  real,    intent(in):: ak(npz+1), bk(npz+1)
  real,    intent(in), dimension(is:ie,js:je):: ps0
  real,    intent(inout), dimension(is:ie,js:je):: ps
  real, intent(in), dimension(isd:ied,jsd:jed,npz):: delp
  real(KIND=4),    intent(in), dimension(is:ie,js:je,kmd):: ta
  real(KIND=4),    intent(in), dimension(is:ie,js:je,kmd):: qa
  real(KIND=4),    intent(out), dimension(is:ie,js:je,npz):: t
  real(KIND=4),    intent(out), dimension(is:ie,js:je,npz):: q
! local:
  real, dimension(is:ie,kmd):: tp, qp
  real, dimension(is:ie,kmd+1):: pe0, pn0
  real, dimension(is:ie,npz):: qn1
  real, dimension(is:ie,npz+1):: pe1, pn1
  integer i,j,k


  do 5000 j=js,je

     do k=1,kmd+1
        do i=is,ie
           pe0(i,k) = ak0(k) + bk0(k)*ps0(i,j)
           pn0(i,k) = log(pe0(i,k))
       enddo
     enddo 
!------
! Model
!------
     do i=is,ie
        pe1(i,1) = ak(1)
     enddo
     do k=1,npz
       do i=is,ie
          pe1(i,k+1) = pe1(i,k) + delp(i,j,k)
       enddo
     enddo
     do i=is,ie
        ps(i,j) = pe1(i,npz+1)
     enddo
     do k=1,npz+1
        do i=is,ie
           pn1(i,k) = log(pe1(i,k))
        enddo
     enddo

   if ( nudge_t .or. nudge_q ) then
        do k=1,kmd
           do i=is,ie
              qp(i,k) = qa(i,j,k)
           enddo
        enddo
        call mappm(kmd, pe0, qp, npz, pe1, qn1, is,ie, 0, kord_data)
        do k=1,npz
           do i=is,ie
              q(i,j,k) = qn1(i,k)
           enddo
        enddo
   endif

   do k=1,kmd
      do i=is,ie
         tp(i,k) = ta(i,j,k)
      enddo
   enddo
   call mappm(kmd, pn0, tp, npz, pn1, qn1, is,ie, 1, kord_data)

   if ( nudge_t ) then
        do k=1,npz
           do i=is,ie
              t(i,j,k) = qn1(i,k)/(1.+zvir*q(i,j,k))
           enddo
        enddo
   else
        do k=1,npz
           do i=is,ie
              t(i,j,k) = qn1(i,k)
           enddo
        enddo
   endif

5000 continue

 end subroutine remap_tq


 subroutine remap_uv(npz, ak, bk, ps, delp, u, v, kmd, ps0, u0, v0)
  integer, intent(in):: npz
  real,    intent(in):: ak(npz+1), bk(npz+1)
  real,    intent(inout):: ps(is:ie,js:je)
  real, intent(in), dimension(isd:ied,jsd:jed,npz):: delp
  real(KIND=4),    intent(inout), dimension(is:ie,js:je,npz):: u, v
!
  integer, intent(in):: kmd
  real,    intent(in):: ps0(is:ie,js:je)
  real(KIND=4),    intent(in), dimension(is:ie,js:je,kmd):: u0, v0
!
! local:
  real, dimension(is:ie,kmd+1):: pe0
  real, dimension(is:ie,npz+1):: pe1
  real, dimension(is:ie,kmd):: qt
  real, dimension(is:ie,npz):: qn1
  integer i,j,k

  do 5000 j=js,je
!------
! Data
!------
     do k=1,kmd+1
       do i=is,ie
          pe0(i,k) = ak0(k) + bk0(k)*ps0(i,j)
       enddo
     enddo
!------
! Model
!------
     do i=is,ie
        pe1(i,1) = ak(1)
     enddo
     do k=1,npz
       do i=is,ie
          pe1(i,k+1) = pe1(i,k) + delp(i,j,k)
       enddo
     enddo
     do i=is,ie
        ps(i,j) = pe1(i,npz+1)
     enddo
!------
! map u
!------
      do k=1,kmd
         do i=is,ie
            qt(i,k) = u0(i,j,k)
         enddo
      enddo
      call mappm(kmd, pe0, qt, npz, pe1, qn1, is,ie, -1, kord_data)
      do k=1,npz
         do i=is,ie
            u(i,j,k) = qn1(i,k)
         enddo
      enddo
!------
! map v
!------
      do k=1,kmd
         do i=is,ie
            qt(i,k) = v0(i,j,k)
         enddo
      enddo
      call mappm(kmd, pe0, qt, npz, pe1, qn1, is,ie, -1, kord_data)
      do k=1,npz
         do i=is,ie
            v(i,j,k) = qn1(i,k)
         enddo
      enddo
5000 continue

 end subroutine remap_uv



 subroutine fv_nwp_nudge_end

    deallocate ( ps_dat )
    deallocate (  t_dat )
    deallocate (  q_dat )

    if ( nudge_winds ) then
         deallocate ( u_dat )
         deallocate ( v_dat )
    endif

    deallocate ( s2c )
    deallocate ( id1 )
    deallocate ( id2 )
    deallocate ( jdc )

    deallocate ( ak0 )
    deallocate ( bk0 )
    deallocate ( lat ) 
    deallocate ( lon ) 

    deallocate ( gz0 ) 

 end subroutine fv_nwp_nudge_end


 subroutine get_tc_mask(time, mask)
      real :: slp_mask = 100900.    ! crtical SLP to apply mask
! Input
      type(time_type), intent(in):: time
      real, intent(inout):: mask(is:ie,js:je)
! local
      real:: pos(2)
      real:: slp_o         ! sea-level pressure (Pa)
      real:: w10_o         ! 10-m wind
      real:: r_vor, p_vor
      real:: dist
      integer n, i, j

    do 5000 n=1,nstorms      ! looop through all storms
!----------------------------------------
! Obtain slp observation
!----------------------------------------
      call get_slp_obs(time, nobs_tc(n), x_obs(1,n), y_obs(1,n), wind_obs(1,n),  mslp_obs(1,n), mslp_out(1,n), rad_out(1,n),   &
                       time_tc(1,n), pos(1), pos(2), w10_o, slp_o, r_vor, p_vor)

      if ( slp_o<880.E2 .or. slp_o>min(slp_env,slp_mask) .or. abs(pos(2))*rad2deg>40. ) goto 5000  ! next storm

      if ( r_vor < 30.E3 ) then
           r_vor = r_min + (slp_env-slp_o)/20.E2*r_inc   ! radius of influence
      endif

      do j=js, je
         do i=is, ie
            dist = great_circle_dist(pos, agrid(i,j,1:2), radius)
            if( dist < 5.*r_vor  ) then 
                if ( strong_mask ) then
                     mask(i,j) = mask(i,j) * ( 1. - exp(-(0.50*dist/r_vor)**2)*min(1.,(slp_env-slp_o)/5.E2) )
                else
! Better analysis data (NCEP 2007 and later) may use a weak mask
                     mask(i,j) = mask(i,j) * ( 1. - exp(-(0.70*dist/r_vor)**2)*min(1.,(slp_env-slp_o)/10.E2) )
                endif
            endif
         enddo             ! i-loop
      enddo                ! end j-loop

5000 continue

 end subroutine get_tc_mask


 subroutine breed_slp_inline(nstep, dt, npz, ak, bk, phis, pe, pk, peln, delp, u, v, pt, q, nwat, zvir)
!------------------------------------------------------------------------------------------
! Purpose:  Vortex-breeding by nudging sea-level-pressure towards single point observations
! Note: conserve water mass, geopotential, and momentum at the expense of dry air mass
!------------------------------------------------------------------------------------------
! Input
      integer, intent(in):: nstep, npz, nwat
      real, intent(in):: dt       ! (small) time step in seconds
      real, intent(in):: zvir
      real, intent(in), dimension(npz+1):: ak, bk
      real, intent(in):: phis(isd:ied,jsd:jed)
! Input/Output
      real, intent(inout):: u(isd:ied,jsd:jed+1,npz)
      real, intent(inout):: v(isd:ied+1,jsd:jed,npz)
      real, intent(inout), dimension(isd:ied,jsd:jed,npz):: delp, pt
      real, intent(inout)::q(isd:ied,jsd:jed,npz,*)

      real, intent(inout):: pk(is:ie,js:je, npz+1)          ! pe**kappa
      real, intent(inout):: pe(is-1:ie+1, npz+1,js-1:je+1)  ! edge pressure (pascal)
      real, intent(out):: peln(is:ie,npz+1,js:je)           ! ln(pe)
! local
      type(time_type):: time
      real:: ps(is:ie,js:je)
      real:: dist(is:ie,js:je)
      real::   tm(is:ie,js:je)
      real::  slp(is:ie,js:je)
      real:: pos(2)
      real:: slp_o         ! sea-level pressure (Pa)
      real:: w10_o, p_env, pre_env
      real:: r_vor
      real:: relx0, relx, f1, pbreed, pbtop, pkz, delp0, dp0
      real:: ratio, p_count, p_sum, a_sum, mass_sink, delps
      real:: p_lo, p_hi, tau_vt
      real:: split_time, fac, pdep, r2, r3
      integer year, month, day, hour, minute, second
      integer n, i, j, k, iq, k0

    if ( nstorms==0 ) then
         if(master) write(*,*) 'NO TC data to process'
         return
    endif

   if ( k_breed==0 ) k_breed = ks
!  k_breed = ks

   k0 = k_breed

! Advance (local) time
    call get_date(fv_time, year, month, day, hour, minute, second)
    if ( year /= year_track_data ) then
        if (master) write(*,*) 'Warning: The year in storm track data is not the same as model year' 
        return
     endif
    time = fv_time   ! fv_time is the time at past time step (set in fv_diag)
    split_time = calday(year, month, day, hour, minute, second) + dt*real(nstep)/86400.

    elapsed_time = elapsed_time + dt
    if ( elapsed_time > nudged_time + 0.1 ) then
         if(print_end_breed)  then
            print_end_breed = .false.
            if (master) write(*,*) '*** Vortext Breeding Ended at', day, hour, minute, second
         endif
         return        !  time to return to forecast mode
    endif

    do j=js,je
! ---- Compute ps
       do i=is,ie
          ps(i,j) = ak(1)
       enddo
       do k=1,npz
          do i=is,ie
             ps(i,j) = ps(i,j) + delp(i,j,k)
          enddo
       enddo
! Compute lowest layer air temperature:
       do i=is,ie
              pkz = (pk(i,j,npz+1)-pk(i,j,npz))/(kappa*log(ps(i,j)/(ps(i,j)-delp(i,j,npz))))
          tm(i,j) = pkz*pt(i,j,npz)/(cp_air*(1.+zvir*q(i,j,npz,1)))
       enddo
    enddo

    do k=k_breed+1,npz

       if ( conserve_mom ) then
       do j=js,je+1
          do i=is,ie
             u(i,j,k) = u(i,j,k) * (delp(i,j-1,k)+delp(i,j,k))
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             v(i,j,k) = v(i,j,k) * (delp(i-1,j,k)+delp(i,j,k))
          enddo
       enddo
       endif

       do j=js,je
          do i=is,ie
             pt(i,j,k) = pt(i,j,k)*(pk(i,j,k+1)-pk(i,j,k))
          enddo
       enddo
! Convert tracer moist mixing ratio to mass
       do iq=1,nwat
          do j=js,je
             do i=is,ie
                q(i,j,k,iq) = q(i,j,k,iq) * delp(i,j,k)
             enddo
          enddo
       enddo

    enddo

    do 5000 n=1,nstorms      ! looop through all storms

!----------------------------------------
! Obtain slp observation
!----------------------------------------
      call get_slp_obs(time, nobs_tc(n), x_obs(1,n), y_obs(1,n), wind_obs(1,n),  mslp_obs(1,n), mslp_out(1,n), rad_out(1,n),   &
                       time_tc(1,n), pos(1), pos(2), w10_o, slp_o, r_vor, p_env, stime=split_time, fact=fac)

      if ( slp_o<87500. .or. slp_o>slp_env .or. abs(pos(2))*rad2deg>45. ) then
           goto 5000         ! next storm
      endif

      if(nudge_debug .and. master)    &
         write(*,*) 'Vortex breeding for TC:', n, ' slp=',slp_o/100.,pos(1)*rad2deg,pos(2)*rad2deg

#ifndef CONST_BREED_DEPTH
! Determine pbtop (top pressure of vortex breeding)

      if ( slp_o > 1000.E2 ) then
!          pbtop = 850.E2
           pbtop = 900.E2
      else
!          pbtop = max(125.E2, 850.E2-30.*(1000.E2-slp_o))
!          pbtop = max(125.E2, 900.E2-30.*(1000.E2-slp_o))
           pbtop = max(120.E2, 900.E2-40.*(1000.E2-slp_o))
      endif

      do k=1,npz
         pbreed = ak(k) + bk(k)*1000.E2
         if ( pbreed>pbtop ) then
              k0 = k
              exit
         endif
      enddo
      k0 = max(k0, k_breed)
#endif

      do j=js, je
         do i=is, ie
            dist(i,j) = great_circle_dist( pos, agrid(i,j,1:2), radius)
         enddo
      enddo

      call compute_slp(is, ie, js, je, tm, ps, phis(is:ie,js:je), slp)

    if ( r_vor < 30.E3 .or. p_env<900.E2 ) then

! Compute r_vor & p_env
         r_vor = r_min + (slp_env-slp_o)/25.E2*r_inc

123   continue
      p_count = 0.
        p_sum = 0.
        a_sum = 0.
      do j=js, je
         do i=is, ie
            if( dist(i,j)<(r_vor+del_r) .and. dist(i,j)>r_vor .and. phis(i,j)<500.*grav ) then 
                p_count = p_count + 1.
                  p_sum = p_sum + slp(i,j)*area(i,j) 
                  a_sum = a_sum + area(i,j) 
            endif
         enddo
      enddo

      call mp_reduce_sum(p_count)

      if ( p_count<32. ) then
           if(nudge_debug .and. master) write(*,*) p_count, 'Skipping obs: too few p_count'
           goto 5000
      endif

      call mp_reduce_sum(p_sum)
      call mp_reduce_sum(a_sum)
      p_env = p_sum / a_sum

      if(nudge_debug .and. master) write(*,*) 'Environmental SLP=', p_env/100., ' computed radius=', r_vor/1.E3

      if ( p_env>1020.E2 .or. p_env<900.E2 ) then
         if( nudge_debug ) then
            if(master)  write(*,*) 'Environmental SLP out of bound; skipping obs. p_count=', p_count, p_sum
            call prt_maxmin('SLP_breeding', slp, is, ie, js, je, 0, 1, 0.01, master)
         endif
         goto 5000
      endif

    endif

      pre_env = pre0_env

      if ( p_env < max(pre_env, slp_o + 250.0) ) then
         if(nudge_debug .and. master) then
            write(*,*) 'Computed environmental SLP too low'
            write(*,*) ' ', p_env/100., slp_o/100.,pos(1)*rad2deg, pos(2)*rad2deg
         endif

         if ( r_vor < 850.E3 ) then
              r_vor = r_vor + del_r
              if(nudge_debug .and. master) write(*,*) 'Vortex radius (km) increased to:', r_vor/1.E3
              goto 123
         else
              p_env = max( slp_o + 250.0, 1000.E2)
         endif
      endif

!     tau_vt = tau_vt_slp + 6.*(980.E2-slp_o)/100.
      tau_vt = tau_vt_slp * (1. + (960.E2-slp_o)/100.E2 )

      tau_vt = max(dt, tau_vt)

      if ( time_track ) then
           relx0  = min(1., fac*dt/tau_vt)
      else
           relx0  = min(1., dt/tau_vt)
      endif

      mass_sink = 0.
      do j=js, je
         do i=is, ie
            if( dist(i,j) < r_vor .and. phis(i,j)<500.*grav ) then
                f1 = dist(i,j)/r_vor
                relx = relx0*exp( -tau_vt_rad*f1**2 )
! Compute p_obs: assuming local radial distributions of slp are Gaussian
                p_hi = p_env - (p_env-slp_o) * exp( -r_hi*f1**2 )    ! upper bound
                p_lo = p_env - (p_env-slp_o) * exp( -r_lo*f1**2 )    ! lower bound

                if ( ps(i,j) > p_hi ) then 
! Under-development:
                     delps = relx*(ps(i,j) - p_hi)   ! Note: ps is used here to prevent
                                                     !       over deepening over terrain
                elseif ( slp(i,j) < p_lo ) then
! Over-development:
                      relx = max(0.5, relx0)
                     delps = relx*(slp(i,j) - p_lo)  ! Note: slp is used here
                else
                     goto 400        ! do nothing; proceed to next storm
                endif 

                if ( delps > 0. ) then
                      pbreed = ak(1)
                      do k=1,k0
                         pbreed = pbreed + delp(i,j,k)
                      enddo
                      f1 = 1. - delps/(ps(i,j)-pbreed)
                      do k=k0+1,npz
                         delp(i,j,k) = delp(i,j,k)*f1
                      enddo
                      mass_sink = mass_sink + delps*area(i,j)
                else
                      dp0 = abs(delps)
                      do k=npz,k0+1,-1
                         if ( abs(delps) < 1. ) then
                              delp(i,j,k) = delp(i,j,k) - delps
                              mass_sink = mass_sink + delps*area(i,j)
                              go to 400
                         else
!                             pdep = max(1.0, min(abs(0.5*delps), 0.2*dp0,  0.02*delp(i,j,k)))
                              pdep = max(1.0, min(abs(0.4*delps), 0.2*dp0,  0.02*delp(i,j,k)))
                              pdep = - min(pdep, abs(delps))
                              delp(i,j,k) = delp(i,j,k) - pdep
                              delps = delps - pdep
                              mass_sink = mass_sink + pdep*area(i,j)
                         endif
                      enddo
                endif

            endif
400     continue 
        enddo        ! end i-loop
      enddo        ! end j-loop

      call mp_reduce_sum(mass_sink)
      if ( abs(mass_sink)<1.E-40 ) goto 5000

      r2 = r_vor + del_r
      r3 = min(2500.E3, 5.*r_vor + del_r)

      p_sum = 0.
      do j=js, je
         do i=is, ie
            if( dist(i,j)<r3 .and. dist(i,j)>r2 ) then
                p_sum = p_sum + area(i,j) 
            endif
         enddo
      enddo

      call mp_reduce_sum(p_sum)
      mass_sink = mass_sink / p_sum ! mean delta pressure to be added back to the environment to conserve mass
      if(master .and. nudge_debug) write(*,*) 'TC#',n, 'Mass tele-ported (pa)=', mass_sink

      do j=js, je
         do i=is, ie
            if( dist(i,j)<r3 .and. dist(i,j)>r2 ) then
                pbreed = ak(1)
                do k=1,k_breed
                   pbreed = pbreed + delp(i,j,k)
                enddo
                f1 = 1. + mass_sink/(ps(i,j)-pbreed)
                do k=k_breed+1,npz
                   delp(i,j,k) = delp(i,j,k)*f1
                enddo
            endif
         enddo
      enddo

! ---- re-compute ps
      do j=js,je
         do i=is,ie
            ps(i,j) = ak(1)
         enddo
         do k=1,npz
            do i=is,ie
               ps(i,j) = ps(i,j) + delp(i,j,k)
            enddo
         enddo
      enddo

5000 continue

!--------------------------
! Update delp halo regions:
!--------------------------
    call mpp_update_domains(delp, domain, complete=.true.)

    do j=js-1,je+1
       do i=is-1,ie+1
          pe(i,1,j) = ak(1)
       enddo
       do k=2,npz+1
          do i=is-1,ie+1
             pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
          enddo
       enddo
    enddo

    do k=k_breed+1,npz+1
      do j=js,je
         do i=is,ie
            peln(i,k,j) = log(pe(i,k,j))
              pk(i,j,k) = pe(i,k,j)**kappa
         enddo
      enddo
    enddo


    do k=k_breed+1,npz

       if ( conserve_mom ) then
       do j=js,je+1
          do i=is,ie
             u(i,j,k) = u(i,j,k) / (delp(i,j-1,k)+delp(i,j,k))
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             v(i,j,k) = v(i,j,k) / (delp(i-1,j,k)+delp(i,j,k))
          enddo
       enddo
       endif

       do j=js,je
          do i=is,ie
             pt(i,j,k) = pt(i,j,k) / (pk(i,j,k+1)-pk(i,j,k))
          enddo
       enddo
    enddo


! Convert tracer mass back to moist mixing ratio
    do iq=1,nwat
       do k=k_breed+1,npz
          do j=js,je
             do i=is,ie
                q(i,j,k,iq) = q(i,j,k,iq) / delp(i,j,k)
             enddo
          enddo
       enddo
    enddo

    call mpp_update_domains(pt, domain, complete=.true.)

  end subroutine breed_slp_inline


 subroutine breed_srf_winds(time, dt, npz, u_obs, v_obs, ak, bk, ps, phis, delp, ua, va, u_dt, v_dt, pt, q, nwat, zvir)
!------------------------------------------------------------------------------------------
! Purpose:  Vortex-breeding by nudging sea-level-pressure towards single point observations
! Note: conserve water mass, geopotential, and momentum at the expense of dry air mass
!------------------------------------------------------------------------------------------
! Input
      type(time_type), intent(in):: time
      integer, intent(in):: npz, nwat
      real, intent(in):: dt       ! time step in seconds
      real, intent(in):: zvir
      real, intent(in), dimension(npz+1):: ak, bk
      real, intent(in):: phis(isd:ied,jsd:jed)
      real, intent(in)::   ps(isd:ied,jsd:jed)
      real, intent(in), dimension(is:ie,js:je,npz):: u_obs, v_obs
! Input/Output
      real, intent(inout), dimension(isd:ied,jsd:jed,npz):: delp, pt, ua, va, u_dt, v_dt
      real, intent(inout)::q(isd:ied,jsd:jed,npz,nwat)
! local
      real:: dist(is:ie,js:je), wind(is:ie,js:je)
      real::  slp(is:ie,js:je)
      real:: pos(2)
      real:: slp_o         ! sea-level pressure (Pa)
      real:: w10_o, p_env
      real:: r_vor, pc, p_count
      real:: r_max, speed, ut, vt, speed_local        ! tangent wind speed
      real:: u_bg, v_bg, mass, t_mass
      real:: relx0, relx, f1, rdt
      real:: z0
      real:: zz = 35.
!     real:: wind_fac = 1.2     ! adjustment factor to account for departure from 10 meter wind
                                ! Computed using Moon et al 2007 
      real:: wind_fac
      integer n, i, j, k, iq

    if ( nstorms==0 ) then
         if(master) write(*,*) 'NO TC data to process'
         return
    endif

       rdt = 1./dt
    relx0  = min(1., dt/tau_vt_wind)

    do j=js, je
       do i=is, ie
           slp(i,j) = ps(i,j)*exp(phis(i,j)/(rdgas*(pt(i,j,npz)+3.25E-3*phis(i,j)/grav)))
          wind(i,j) = sqrt( ua(i,j,npz)**2 + va(i,j,npz)**2 )
       enddo
    enddo

    do 3000 n=1,nstorms      ! looop through all storms

!----------------------------------------
! Obtain slp observation
!----------------------------------------
      call get_slp_obs(time, nobs_tc(n), x_obs(1,n), y_obs(1,n), wind_obs(1,n),  mslp_obs(1,n), mslp_out(1,n), rad_out(1,n),   &
                       time_tc(1,n), pos(1), pos(2), w10_o, slp_o, r_vor, p_env)

      if ( slp_o<87000. .or. slp_o>slp_env .or. abs(pos(2))*rad2deg>35. ) goto 3000         ! next storm

      do j=js, je
         do i=is, ie
            dist(i,j) = great_circle_dist( pos, agrid(i,j,1:2), radius )
         enddo
      enddo

      r_vor = r_min + (slp_env-slp_o)/25.E2*r_inc

!----------------------------------------------------
     if ( w10_o < 0. ) then   ! 10-m wind obs is not available
! Uses Atkinson_Holliday wind-pressure correlation
          w10_o = 3.446778 * (1010.-slp_o/100.)**0.644
     endif

! * Find model's SLP center nearest to the observation
! * Find maximum wind speed at the lowest model level

     speed_local = 0.
           r_max = -999.
              pc = 1013.E2
     do j=js, je
        do i=is, ie
           if( dist(i,j) < r_vor ) then

               pc = min(pc, slp(i,j))

               if ( speed_local < wind(i,j) ) then
                    speed_local = wind(i,j)
                    r_max = dist(i,j)
               endif

           endif
        enddo
     enddo

     speed = speed_local
     call mp_reduce_max(speed)     ! global max wind (near storm)
     call mp_reduce_min(pc)

    if ( speed_local < speed ) then
         r_max = -999.
    endif
    call mp_reduce_max(r_max)
    if( r_max<0. ) call mpp_error(FATAL,'==> Error in r_max')

! ---------------------------------------------------
! Determine surface wind speed and radius for nudging 
! ---------------------------------------------------

! Compute surface roughness z0 from w10, based on Eq (4) & (5) from Moon et al. 2007
     if ( w10_o > 12.5 ) then
          z0 = (0.085*w10_o - 0.58) * 1.E-3
     else
          z0 = 0.0185/grav*(0.001*w10_o**2 + 0.028*w10_o)**2
     endif

! lowest layer height: zz

     wind_fac = log(zz/z0) / log(10./z0)
     if( nudge_debug .and. master ) write(*,*) 'Wind adjustment factor=', wind_fac
     if( wind_fac<1. ) call mpp_error(FATAL,'==> Error in wind_fac')

     if ( pc < slp_o ) then
!--
!         The storm in the model is over developed
!         if ( (pc+3.0E2)>slp_o .or. speed <= w10_o ) go to 3000    ! next storm
!--
! using radius (r_max) as dtermined above;
! What if model's pressure center is very far from the observed?
! using obs wind
          speed = wind_fac*w10_o
     else
!         The storm in the model is under developed; using max wind
          speed = max(wind_fac*w10_o, speed)
          if ( pc>1009.E2 )  r_max = 0.5 * r_vor
     endif

! Some bounds on the radius of maximum wind:
     r_max = max(2.5*grid_size, r_max)      ! at least 2.5X the grid size
     r_max = min(0.75*r_vor, r_max)

     t_mass = 0.
     u_bg = 0.
     v_bg = 0.

     if ( add_bg_wind ) then
       p_count = 0.
       do j=js, je
          do i=is, ie
           if( dist(i,j) <= min(r_vor,2.*r_fac*r_max) .and. phis(i,j)<1.0*grav ) then
               mass = area(i,j)*delp(i,j,npz)
!-- using model winds ----------------------------------
!              u_bg = u_bg + ua(i,j,npz)*mass
!              v_bg = v_bg + va(i,j,npz)*mass
!-------------------------------------------------------
! Using analysis winds
               u_bg = u_bg + u_obs(i,j,npz)*mass
               v_bg = v_bg + v_obs(i,j,npz)*mass
               t_mass = t_mass + mass
               p_count = p_count + 1.
           endif
          enddo
       enddo
       call mp_reduce_sum(p_count)
       if ( p_count<16. ) go to 3000

       call mp_reduce_sum(t_mass)
       call mp_reduce_sum(u_bg)
       call mp_reduce_sum(v_bg)
       u_bg = u_bg / t_mass
       v_bg = v_bg / t_mass
!      if ( master ) write(*,*) pos(2)*rad2deg, 'vortex bg wind=', u_bg, v_bg
     endif

     relx = relx0
     k = npz                 ! lowest layer only
! Nudge wind in the "inner core":
      do j=js, je
         do i=is, ie
            if( dist(i,j) <= min(r_vor, r_fac*r_max) .and. phis(i,j)<1.0*grav ) then
                f1 = dist(i,j)/r_max
!               relx = relx0*exp( -f1**2 )
                relx = relx0*exp( -tau_vt_rad*f1**2 )
                if( dist(i,j)<=r_max ) then
                    speed_local = speed * f1
                else
                    speed_local = speed / f1**0.75
                endif
                call tangent_wind(vlon(i,j,1:3), vlat(i,j,1:3), speed_local, pos, agrid(i,j,1:2), ut, vt)
                ut = ut + u_bg
                vt = vt + v_bg
                u_dt(i,j,k) = u_dt(i,j,k) + relx*(ut-ua(i,j,k)) * rdt
                v_dt(i,j,k) = v_dt(i,j,k) + relx*(vt-va(i,j,k)) * rdt
! Update:
                ua(i,j,k) = ua(i,j,k) + relx*(ut-ua(i,j,k))
                va(i,j,k) = va(i,j,k) + relx*(vt-va(i,j,k))
            endif
400     continue 
        enddo        ! end i-loop
      enddo        ! end j-loop

3000 continue

  end subroutine breed_srf_winds

  subroutine tangent_wind ( elon, elat, speed, po, pp, ut, vt )
  real, intent(in):: speed
  real, intent(in):: po(2), pp(2)
  real, intent(in):: elon(3), elat(3)
  real, intent(out):: ut, vt
! local
  real:: e1(3), eo(3), ep(3), op(3)

  call latlon2xyz(po, eo)
  call latlon2xyz(pp, ep)

  op(:) = ep(:) - eo(:)
  call normalize_vect( op )

  call vect_cross(e1, ep, eo)

  ut = speed * (e1(1)*elon(1) + e1(2)*elon(2) + e1(3)*elon(3))
  vt = speed * (e1(1)*elat(1) + e1(2)*elat(2) + e1(3)*elat(3))

! SH:
  if ( po(2) < 0. ) then
       ut = -ut
       vt = -vt
  endif

  end subroutine tangent_wind


  subroutine get_slp_obs(time, nobs, lon_obs, lat_obs, w10, mslp, slp_out, r_out, time_obs,    &
                         x_o, y_o, w10_o, slp_o, r_vor, p_vor, stime, fact)
! Input
    type(time_type), intent(in):: time
    integer, intent(in)::  nobs   ! number of observations in this particular storm
    real(KIND=4), intent(in)::  lon_obs(nobs)
    real(KIND=4), intent(in)::  lat_obs(nobs)
    real(KIND=4), intent(in)::      w10(nobs)        ! observed 10-m widn speed
    real(KIND=4), intent(in)::     mslp(nobs)        ! observed SLP in pa
    real(KIND=4), intent(in)::  slp_out(nobs)        ! slp at r_out
    real(KIND=4), intent(in)::    r_out(nobs)        ! 
    real(KIND=4), intent(in):: time_obs(nobs)
    real, optional, intent(in):: stime
    real, optional, intent(out):: fact
! Output
    real, intent(out):: x_o , y_o      ! position of the storm center 
    real, intent(out):: w10_o          ! 10-m wind speed
    real, intent(out):: slp_o          ! Observed sea-level-pressure (pa)
    real, intent(out):: r_vor, p_vor
! Internal:
    real:: t_thresh
      real:: p1(2), p2(2)
      real time_model
      real fac
      integer year, month, day, hour, minute, second, n

      t_thresh = 600./86400.  ! unit: days

       w10_o = -100000.
       slp_o = -100000.
         x_o = -100.*pi
         y_o = -100.*pi
       p_vor = -1.E10
       r_vor = -1.E10

   if ( present(stime) ) then
      time_model = stime
   else
      call get_date(time, year, month, day, hour, minute, second)

      if ( year /= year_track_data ) then
           if (master) write(*,*) 'Warning: The year in storm track data is not the same as model year' 
           return
      endif

      time_model = calday(year, month, day, hour, minute, second)
!     if(nudge_debug .and. master) write(*,*) 'Model:', time_model, year, month, day, hour, minute, second
   endif

!-------------------------------------------------------------------------------------------
!     if ( time_model <= time_obs(1)  .or.  time_model >= time_obs(nobs) ) then
!          return
!-------------------------------------------------------------------------------------------

      if ( time_model <= (time_obs(1)-t_thresh)  .or.  time_model >= time_obs(nobs) ) return

      if ( time_model <=  time_obs(1) ) then
!--
! This is an attempt to perform vortex breeding several minutes before the first available observation
!--
                 w10_o =     w10(1)
                 slp_o =    mslp(1)
                   x_o = lon_obs(1)
                   y_o = lat_obs(1)
                 if ( present(fact) )  fact = 1.25
      else
           do n=1,nobs-1
             if( time_model >= time_obs(n) .and. time_model <= time_obs(n+1) ) then
                   fac = (time_model-time_obs(n)) / (time_obs(n+1)-time_obs(n))
                 w10_o =     w10(n) + (    w10(n+1)-    w10(n)) * fac
                 slp_o =    mslp(n) + (   mslp(n+1)-   mslp(n)) * fac
! Trajectory interpolation:
! Linear in (lon,lat) space
!                  x_o = lon_obs(n) + (lon_obs(n+1)-lon_obs(n)) * fac
!                  y_o = lat_obs(n) + (lat_obs(n+1)-lat_obs(n)) * fac
                 p1(1) = lon_obs(n);     p1(2) = lat_obs(n)
                 p2(1) = lon_obs(n+1);   p2(2) = lat_obs(n+1)
                 call intp_great_circle(fac, p1, p2, x_o, y_o)
!----------------------------------------------------------------------
                  if ( present(fact) )   fact = 1. + 0.25*cos(fac*2.*pi)
! Additional data from the extended best track
!                if ( slp_out(n)>0. .and. slp_out(n+1)>0. .and. r_out(n)>0. .and. r_out(n+1)>0. ) then
!                     p_vor = slp_out(n) + ( slp_out(n+1) - slp_out(n)) * fac
!                     r_vor =   r_out(n) + (   r_out(n+1) -   r_out(n)) * fac
!                endif
                 return
             endif
           enddo
      endif

  end subroutine get_slp_obs


  subroutine slp_obs_init
  integer:: unit, n, nobs
  character(len=3):: GMT
  character(len=9):: ts_name
  character(len=19):: comment
  integer:: mmddhh, yr, year, month, day, hour, MPH, islp
  integer:: it, i1, i2, p_ring, d_ring
  real:: lon_deg, lat_deg, cald, slp, mps

  nobs_tc(:) = 0
  time_tc(:,:) = 0.
  wind_obs(:,:) = -100000.
  mslp_obs(:,:) = -100000.
  x_obs(:,:) = - 100.*pi
  y_obs(:,:) = - 100.*pi

  mslp_out(:,:) = -1.E10
   rad_out(:,:) = -1.E10

  if( track_file_name == "No_File_specified" ) then
      if(master) write(*,*) 'No TC track file specified'
      return
  else
      unit = get_unit()
      open( unit, file=track_file_name)
  endif

  read(unit, *) year
  if(master) write(*,*) 'Reading TC track data for YEAR=', year

  year_track_data = year

  nstorms = 0
     nobs = 0
    month = 99

  if ( ibtrack ) then

!---------------------------------------------------------------
! The data format is from Ming Zhoa's processed ibTrack datasets
!---------------------------------------------------------------

    read(unit, *) ts_name, nobs, yr, month, day, hour

    if ( yr /= year ) then
         if(master) write(*, *) 'Year inconsistency found !!!'
         call mpp_error(FATAL,'==> Error in reading best track data')
    endif

    do while ( ts_name=='start' ) 

               nstorms  = nstorms + 1
       nobs_tc(nstorms) = nobs       ! observation count for this storm
       if(nudge_debug .and. master) write(*, *) 'Read Data for TC#', nstorms, nobs

       do it=1, nobs
          read(unit, *) lon_deg, lat_deg, mps, slp, yr, month, day, hour
!         if ( yr /= year ) then
!             if(master) write(*, *) 'Extended to year + 1', yr
!         endif
          cald = calday(yr, month, day, hour, 0, 0)
          time_tc(it,nstorms) = cald
          if(nudge_debug .and. master) write(*, 100) cald, month, day, hour, lon_deg, lat_deg, mps, slp

          wind_obs(it,nstorms) = mps       ! m/s
          mslp_obs(it,nstorms) = 100.*slp
             y_obs(it,nstorms) = lat_deg * deg2rad
             x_obs(it,nstorms) = lon_deg * deg2rad
       enddo

       read(unit, *) ts_name, nobs, yr, month, day, hour
    enddo
100  format(1x, f9.2, 1x, i3, 1x, i2, 1x, i2, 1x, f6.1, 1x, f6.1, 1x, f4.1, 1x, f6.1)

  else

  do while ( month /= 0 )

     read(unit, *) month, day, hour, GMT, lat_deg, lon_deg, MPH, islp, comment

     select case (month)

     case (99)                ! New storm record to start
          nstorms = nstorms + 1
          nobs = 0
          if(master) write(*, *) 'Reading data for TC#', nstorms, comment
     case ( 0)                ! end of record
          if(master) write(*, *) 'End of record reached'
     case default
           nobs = nobs + 1
           cald = calday(year, month, day, hour, 0, 0)
           time_tc(nobs,nstorms) = cald
           nobs_tc(nstorms) = nobs       ! observation count for this storm

          if(master) write(*, 200) nobs, cald,  month, day, hour, GMT, lat_deg, lon_deg, MPH, islp, comment
          mslp_obs(nobs,nstorms) = 100. * real(islp)
             y_obs(nobs,nstorms) = lat_deg * deg2rad
          if ( GMT == 'GMT' ) then
!                                  Transfrom x from (-180 , 180) to (0, 360) then to radian
             if ( lon_deg < 0 ) then 
                  x_obs(nobs,nstorms) = (360.+lon_deg) * deg2rad
             else
                  x_obs(nobs,nstorms) = (360.-lon_deg) * deg2rad
             endif
          elseif ( GMT == 'PAC' ) then   ! Pacific storms
             x_obs(nobs,nstorms) = lon_deg * deg2rad
          endif
     end select

  enddo

  endif

  close(unit)

  if(master) then 
     write(*,*) 'TC vortex breeding: total storms=', nstorms
     if ( nstorms/=0 ) then
          do n=1,nstorms
             write(*, *) 'TC#',n, ' contains ',  nobs_tc(n),' observations'
          enddo
     endif
  endif

200  format(i3, 1x,f9.4, 1x, i2, 1x, i2, 1x, i2, 1x, a3, f5.1, 1x, f5.1, 1x, i3, 1x, i4, 1x, a19)

  end subroutine slp_obs_init


  real function calday(year, month, day, hour, minute, sec)
! For time interpolation; Julian day (0 to 365 for non-leap year)
! input:
    integer, intent(in):: year, month, day, hour
    integer, intent(in):: minute, sec
! Local:
      integer n, m, ds, nday
      real tsec
      integer days(12)
      data days /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

      ds = day - 1

      if( month /= 1 ) then
          do m=1, month-1
            if( m==2  .and. leap_year(year) ) then 
                ds = ds + 29
            else
                ds = ds + days(m)
            endif
          enddo
      endif

      if ( leap_year(year_track_data) ) then
           nday = 366
      else
           nday = 365
      endif

      calday = real((year-year_track_data)*nday + ds)  + real(hour*3600 + minute*60 + sec)/86400.

  end function calday


  logical function leap_year(ny)
  integer, intent(in):: ny
!
! Determine if year ny is a leap year
! Author: S.-J. Lin
   integer ny00
!
! No leap years prior to 0000
!
      parameter ( ny00 = 0000 )   ! The threshold for starting leap-year 

      if( ny >= ny00 ) then
         if( mod(ny,100) == 0. .and. mod(ny,400) == 0. ) then
             leap_year = .true.
         elseif( mod(ny,4) == 0. .and. mod(ny,100) /= 0.  ) then
             leap_year = .true.
         else
             leap_year = .false.
         endif
      else
          leap_year = .false.
      endif

  end function leap_year


 subroutine pmaxmin( qname, a, imax, jmax, fac )

      character(len=*)  qname
      integer imax, jmax
      integer i, j
      real a(imax,jmax)

      real qmin(jmax), qmax(jmax)
      real pmax, pmin
      real fac                     ! multiplication factor

      do j=1,jmax
         pmax = a(1,j)
         pmin = a(1,j)
         do i=2,imax
            pmax = max(pmax, a(i,j))
            pmin = min(pmin, a(i,j))
         enddo
         qmax(j) = pmax
         qmin(j) = pmin
      enddo
!
! Now find max/min of amax/amin
!
            pmax = qmax(1)
            pmin = qmin(1)
         do j=2,jmax
            pmax = max(pmax, qmax(j))
            pmin = min(pmin, qmin(j))
         enddo

      write(*,*) qname, ' max = ', pmax*fac, ' min = ', pmin*fac

 end subroutine pmaxmin


 subroutine del2_uv(du, dv, cd, kmd, ntimes)
! This routine is for filtering the wind tendency
   integer, intent(in):: kmd
   integer, intent(in):: ntimes
   real,    intent(in):: cd            ! cd = K * da_min;   0 < K < 0.25
   real, intent(inout):: du(is:ie,js:je,kmd)
   real, intent(inout):: dv(is:ie,js:je,kmd)
! local:
   real, dimension(is:ie,js:je,kmd):: v1, v2, v3
   integer i,j,k

! transform to 3D Cartesian:
   do k=1,kmd
      do j=js,je
         do i=is,ie
            v1(i,j,k) = du(i,j,k)*vlon(i,j,1) + dv(i,j,k)*vlat(i,j,1)
            v2(i,j,k) = du(i,j,k)*vlon(i,j,2) + dv(i,j,k)*vlat(i,j,2)
            v3(i,j,k) = du(i,j,k)*vlon(i,j,3) + dv(i,j,k)*vlat(i,j,3)
         enddo
      enddo
   enddo

! Filter individual components as scalar:
   call del2_scalar( v1(is,js,1), cd, kmd, ntimes )
   call del2_scalar( v2(is,js,1), cd, kmd, ntimes )
   call del2_scalar( v3(is,js,1), cd, kmd, ntimes )

! Convert back to lat-lon components:
   do k=1,kmd
      do j=js,je
         do i=is,ie
            du(i,j,k) = v1(i,j,k)*vlon(i,j,1) + v2(i,j,k)*vlon(i,j,2) + v3(i,j,k)*vlon(i,j,3)
            dv(i,j,k) = v1(i,j,k)*vlat(i,j,1) + v2(i,j,k)*vlat(i,j,2) + v3(i,j,k)*vlat(i,j,3)
         enddo
      enddo
   enddo

 end subroutine del2_uv

 subroutine del2_scalar(qdt, cd, kmd, ntimes)
! This routine is for filtering the physics tendency
   integer, intent(in):: kmd
   integer, intent(in):: ntimes
   real,    intent(in):: cd            ! cd = K * da_min;   0 < K < 0.25
   real, intent(inout):: qdt(is:ie,js:je,kmd)
! local:
   real::  q(isd:ied,jsd:jed,kmd)
   real:: fx(isd:ied+1,jsd:jed), fy(isd:ied,jsd:jed+1)
   integer i,j,k, n, nt
   real :: damp

   damp = cd * da_min

   do k=1,kmd
      do j=js,je
         do i=is,ie
            q(i,j,k) = qdt(i,j,k)
         enddo
      enddo
   enddo
                     call timing_on('COMM_TOTAL')
   call mpp_update_domains(q, domain, complete=.true.)
                     call timing_off('COMM_TOTAL')

   do n=1,ntimes

   nt = ntimes - n

   do k=1,kmd

      if(nt>0) call copy_corners(q(isd,jsd,k), npx, npy, 1)
      do j=js-nt,je+nt
         do i=is-nt,ie+1+nt
            fx(i,j) = dy(i,j)*sina_u(i,j)*(q(i-1,j,k)-q(i,j,k))*rdxc(i,j)
         enddo
      enddo

      if(nt>0) call copy_corners(q(isd,jsd,k), npx, npy, 2)
      do j=js-nt,je+1+nt
         do i=is-nt,ie+nt
            fy(i,j) = dx(i,j)*sina_v(i,j)*(q(i,j-1,k)-q(i,j,k))*rdyc(i,j)
         enddo
      enddo

      if ( nt==0 ) then
          do j=js,je
             do i=is,ie
                qdt(i,j,k) = q(i,j,k) + damp*rarea(i,j)*(fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))
             enddo
          enddo
      else
          do j=js-nt,je+nt
             do i=is-nt,ie+nt
                q(i,j,k) = q(i,j,k) + damp*rarea(i,j)*(fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))
             enddo
          enddo
      endif
   enddo

   enddo

 end subroutine del2_scalar

 subroutine rmse_bias(a, rms, bias)
   real, intent(in):: a(is:ie,js:je)
   real, intent(out):: rms, bias
   integer:: i,j
   real:: total_area

   total_area = 4.*pi*radius**2

    rms = 0.
   bias = 0.
   do j=js,je
      do i=is,ie
         bias = bias + area(i,j) * a(i,j)
          rms = rms  + area(i,j) * a(i,j)**2
      enddo
   enddo
   call mp_reduce_sum(bias)
   call mp_reduce_sum(rms)

   bias = bias / total_area
    rms = sqrt( rms / total_area )

 end subroutine rmse_bias


 subroutine corr(a, b, co)
 real, intent(in):: a(is:ie,js:je), b(is:ie,js:je)
 real, intent(out):: co
 real:: m_a, m_b, std_a, std_b
 integer:: i,j
 real:: total_area

   total_area = 4.*pi*radius**2

! Compute standard deviation:
   call std(a, m_a, std_a)
   call std(b, m_b, std_b)

! Compute correlation: 
   co = 0.
   do j=js,je
      do i=is,ie
         co = co + area(i,j) * (a(i,j)-m_a)*(b(i,j)-m_b)
      enddo
   enddo
   call mp_reduce_sum(co)
   co = co / (total_area*std_a*std_b )

 end subroutine corr

 subroutine std(a, mean, stdv)
 real,  intent(in):: a(is:ie,js:je)
 real, intent(out):: mean, stdv
 integer:: i,j
 real:: total_area

   total_area = 4.*pi*radius**2

   mean = 0.
   do j=js,je
      do i=is,ie
         mean = mean + area(i,j) * a(i,j)
      enddo
   enddo
   call mp_reduce_sum(mean)
   mean = mean / total_area 

   stdv = 0.
   do j=js,je
      do i=is,ie
         stdv = stdv + area(i,j) * (a(i,j)-mean)**2
      enddo
   enddo
   call mp_reduce_sum(stdv)
   stdv = sqrt( stdv / total_area )

 end subroutine std


end module fv_nwp_nudge_mod
