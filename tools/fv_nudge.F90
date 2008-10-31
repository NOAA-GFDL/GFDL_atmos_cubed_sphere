module fv_nwp_nudge_mod

 use constants_mod,     only: pi, grav, rdgas, rvgas, cp_air, kappa
 use fms_io_mod,        only: field_size
 use fms_mod,           only: write_version_number, open_namelist_file, &
                              check_nml_error, file_exist, close_file,  &
                              error_mesg, FATAL, read_data, field_exist 
 use mpp_mod,           only: mpp_error, FATAL, stdlog
 use time_manager_mod,  only: time_type,  get_time

 use fv_grid_utils_mod, only: i_sst, j_sst, sst_ncep
 use fv_grid_tools_mod, only: agrid
 use fv_diagnostics_mod,only: prt_maxmin
 use fv_mapz_mod,       only: mappm
 use fv_mp_mod,         only: is,js,ie,je, isd,jsd,ied,jed, gid, masterproc

 implicit none
 private

 character(len=128) :: version = ''
 character(len=128) :: tagname = ''

 public fv_nwp_nudge, fv_nwp_nudge_init, fv_nwp_nudge_end

 integer im     ! Data x-dimension
 integer jm     ! Data y-dimension
 integer km     ! Data z-dimension
 real, allocatable:: ak0(:), bk0(:)
 real, allocatable:: lat(:), lon(:)

 logical :: module_is_initialized = .false.
 logical :: master
 logical :: no_obs
 real :: deg2rad
 real :: time_nudge = 0.
 integer :: time_interval = 6*3600   ! dataset time interval (seconds)
 integer, parameter :: nfile_max = 1601
 integer :: nfile = 1

 real,    allocatable:: s2c(:,:,:)
 integer, allocatable:: id1(:,:), id2(:,:), jdc(:,:)
 real, allocatable :: ps_dat(:,:,:)
 real, allocatable, dimension(:,:,:,:):: u_dat, v_dat, t_dat, q_dat

! Namelist variables:
 character(len=128):: file_names(nfile_max)
 integer :: nfile_total = 0       ! =5 for 1-day (datasets are 6-hr apart)
 real    :: p_wvp = 100.E2        ! cutoff level for specific humidity nudging 
 logical :: nudge_ps    = .false. ! actually delp with FV core
 logical :: nudge_t     = .false.
 logical :: nudge_q     = .false.
 logical :: nudge_winds = .true.
 logical :: nudge_virt  = .true.

! Nudging time-scales (seconds):
 real :: tau_winds  = 21600.       !  6-hr
 real :: tau_virt   = 21600.       !  6-hr
 real :: tau_t      = 21600.       !  6-hr
 real :: tau_q      = 43200.       ! 12-hr
 real :: tau_ps     = 86400.       ! 24-hr

! starting layer (top layer is sponge layer and is skipped)
 integer :: kstart = 2 

! skip "kbot" layers
 integer :: kbot_winds = 1 
 integer :: kbot_t     = 1 
 integer :: kbot_q     = 1 

 namelist /fv_nwp_nudge_nml/ nudge_ps, nudge_t, nudge_winds, nudge_virt,tau_ps, &
                          tau_winds, tau_t, tau_q, tau_virt, kstart, kbot_winds, &
                          kbot_t, kbot_q, p_wvp, time_interval, file_names

 contains
 

  subroutine fv_nwp_nudge ( Time, dt, npz, ps_dt, u_dt, v_dt, t_dt, q_dt, zvir, &
                            ak, bk, ts, ps, delp, ua, va, pt, q, phis )

  type(time_type), intent(in):: Time
  integer,         intent(in):: npz           ! vertical dimension
  real,            intent(in):: dt
  real,            intent(in):: zvir
  real, intent(in   ), dimension(npz+1):: ak, bk
  real, intent(in   ), dimension(isd:ied,jsd:jed    ):: phis
  real, intent(inout), dimension(isd:ied,jsd:jed,npz):: pt, q, ua, va
  real, intent(inout), dimension(isd:ied,jsd:jed,npz):: delp    ! p3 not used
  real, intent(inout), dimension(isd:ied,jsd:jed    ):: ps
! Accumulated tendencies
  real, intent(inout), dimension(isd:ied,jsd:jed,npz):: u_dt, v_dt
  real, intent(out):: t_dt(is:ie,js:je,npz)
  real, intent(out):: q_dt(is:ie,js:je,npz)
  real, intent(out), dimension(is:ie,js:je):: ps_dt, ts
! local:
  real, allocatable :: ps_obs(:,:)
  real, allocatable, dimension(:,:,:):: u_obs, v_obs, t_obs, q_obs
  integer :: seconds, days
  integer :: i,j,k
  real :: dbk, rdt, press(npz), profile(npz), du, dv


  if ( .not. module_is_initialized ) then 
        call mpp_error(FATAL,'==> Error from fv_nwp_nudge: module not initialized')
  endif

  do k=kstart,npz
     press(k) = 0.5*(ak(k) + ak(k+1)) + 0.5*(bk(k)+bk(k))*1.E5
     if ( press(k) > 900.E2 ) then
          profile(k) = (1.E5 - press(k)) / 100.E2 
     else
          profile(k) = 1.
     endif
  enddo
  profile(1) = 0.
  profile(2) = 0.5
  profile(3) = 0.75

  allocate (ps_obs(is:ie,js:je) )
  allocate ( t_obs(is:ie,js:je,npz) )
  allocate ( q_obs(is:ie,js:je,npz) )

  if ( nudge_winds ) then
       allocate (u_obs(is:ie,js:je,npz) )
       allocate (v_obs(is:ie,js:je,npz) )
  endif

  call get_obs(Time, dt, zvir, ak, bk, ts, ps_obs, u_obs, v_obs, t_obs, q_obs, phis, npz)

  if ( no_obs ) return

  if ( nudge_ps ) then
     rdt = 1. / (tau_ps + dt)
     do j=js,je
        do i=is,ie
           ps_dt(i,j) = (ps_obs(i,j) - ps(i,j)) * rdt
        enddo
     enddo
     do k=kstart,npz
        dbk = dt * (bk(k+1) - bk(k))
        do j=js,je
           do i=is,ie
              delp(i,j,k) = delp(i,j,k) + dbk*ps_dt(i,j)
           enddo
        enddo
     enddo
  endif
  deallocate ( ps_obs )

  if ( nudge_winds ) then
     rdt = 1. / (tau_winds + dt)
     do k=kstart, npz - kbot_winds
        do j=js,je
           do i=is,ie
                       du = profile(k)*(u_obs(i,j,k)-ua(i,j,k))*rdt
                       dv = profile(k)*(v_obs(i,j,k)-va(i,j,k))*rdt
              u_dt(i,j,k) = u_dt(i,j,k) + du
              v_dt(i,j,k) = v_dt(i,j,k) + dv
! Update (ua, va)
              ua(i,j,k) = ua(i,j,k) + du*dt
              va(i,j,k) = va(i,j,k) + dv*dt
           enddo
        enddo
     enddo
     deallocate ( u_obs )
     deallocate ( v_obs )
  endif

  if ( nudge_t ) then
     rdt = dt / (tau_t + dt)
     do k=kstart, npz - kbot_t
        do j=js,je
           do i=is,ie
              pt(i,j,k) = pt(i,j,k) + profile(k)*(t_obs(i,j,k)-pt(i,j,k))*rdt
           enddo
        enddo
     enddo
  elseif ( nudge_virt ) then
     rdt = dt / (tau_virt + dt)
     do k=kstart, npz - kbot_t
        do j=js,je
           do i=is,ie
              pt(i,j,k) = pt(i,j,k) + profile(k)*(t_obs(i,j,k)/(1.+zvir*q(i,j,k))-pt(i,j,k))*rdt
           enddo
        enddo
     enddo
  endif
  deallocate ( t_obs )

  if ( nudge_q ) then
! Note: Mass is NOT conserved if q is nudged
     rdt = dt / (tau_q + dt)
     do k=kstart, npz - kbot_q
        if ( press(k) > p_wvp ) then
        do j=js,je
           do i=is,ie
              q(i,j,k) = q(i,j,k) + profile(k)*(q_obs(i,j,k)-q(i,j,k))*rdt
           enddo
        enddo
        endif
     enddo
  endif
  deallocate ( q_obs )

 end  subroutine fv_nwp_nudge



 subroutine get_obs(Time, dt, zvir, ak, bk, ts, ps_obs, u_obs, v_obs, t_obs, q_obs, phis, npz)
  type(time_type), intent(in):: Time
  integer,         intent(in):: npz           ! vertical dimension
  real,            intent(in):: zvir
  real,            intent(in):: dt
  real, intent(in), dimension(npz+1):: ak, bk
  real, intent(in), dimension(isd:ied,jsd:jed):: phis
  real, intent(out), dimension(is:ie,js:je):: ts, ps_obs
  real, intent(out), dimension(is:ie,js:je,npz):: u_obs, v_obs, t_obs, q_obs
! local:
  integer :: seconds, days
  integer :: i,j,k
  real :: ratio

  call get_time (time, seconds, days)

  seconds = seconds - nint(dt)

! Data must be 6-hr apart; keep two time levels in memory

  no_obs = .false.

  if ( mod(seconds, time_interval) == 0 ) then

    if ( nfile > nfile_total ) then
         no_obs = .true.
         return              ! free-running mode
    else

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
      call get_ncep_analysis ( ps_dat(:,:,2), u_dat(:,:,:,2), v_dat(:,:,:,2),    &
                              t_dat(:,:,:,2), q_dat(:,:,:,2), npz, zvir, ak, bk, &
                              ts, phis, nfile, file_names(nfile) )
      time_nudge = dt
    endif
  else
      time_nudge = time_nudge + dt
  endif

!--------------------
! Time interpolation:
!--------------------

  ratio = time_nudge / real(time_interval)

  if ( ratio < 0. .or. ratio >  (1.+1.E-7) ) then
       call mpp_error(FATAL,'==> Error from get_obs:data out of range')
  endif

  if ( nudge_ps ) &
       ps_obs(:,:)  = (1.-ratio)*ps_dat(:,:,1) + ratio*ps_dat(:,:,2)

  if ( nudge_winds ) then
       u_obs(:,:,:) = (1.-ratio)*u_dat(:,:,:,1) + ratio*u_dat(:,:,:,2)
       v_obs(:,:,:) = (1.-ratio)*v_dat(:,:,:,1) + ratio*v_dat(:,:,:,2)
  endif

  if ( nudge_t .or. nudge_virt )  &
       t_obs(:,:,:) = (1.-ratio)*t_dat(:,:,:,1) + ratio*t_dat(:,:,:,2)

  if ( nudge_q )  &
       q_obs(:,:,:) = (1.-ratio)*q_dat(:,:,:,1) + ratio*q_dat(:,:,:,2)

 end subroutine get_obs


 subroutine fv_nwp_nudge_init(npz, zvir, ak, bk, ts, phis)
  integer,  intent(in):: npz           ! vertical dimension 
  real,     intent(in):: zvir
  real, intent(in), dimension(isd:ied,jsd:jed):: phis
  real, intent(in), dimension(npz+1):: ak, bk
  real, intent(out), dimension(is:ie,js:je):: ts
  logical found
  integer tsize(4)
  integer :: i, j, unit, io, ierr, nt

   master = gid==masterproc

   deg2rad = pi/180.

   do nt=1,nfile_max
      file_names(nt) = "No_File_specified"
   enddo

    if( file_exist( 'input.nml' ) ) then
       unit = open_namelist_file ()
       io = 1
       do while ( io .ne. 0 )
          read( unit, nml = fv_nwp_nudge_nml, iostat = io, end = 10 )
          ierr = check_nml_error(io,'fv_nwp_nudge_nml')
       end do
10     call close_file ( unit )
    end if
    call write_version_number (version, tagname)
    if ( master ) then
         write( stdlog(), nml = fv_nwp_nudge_nml )
         write(*,*) 'FV: NWP nudging initialized.'
    endif

    if ( nudge_virt ) then
         nudge_t = .false.
         nudge_q = .false.
    endif

    if ( nudge_t ) nudge_virt = .false.

    do nt=1,nfile_max
      if ( file_names(nt) == "No_File_specified" ) then
           nfile_total = nt - 1
           if(master) write(*,*) 'Total of NCEP files specified=', nfile_total
           exit
      endif
    enddo


    allocate (ps_dat(is:ie,js:je,2) )
    ps_dat = 0.

    if ( nudge_winds ) then
         allocate ( u_dat(is:ie,js:je,npz,2) )
         allocate ( v_dat(is:ie,js:je,npz,2) )
         u_dat = 0.
         v_dat = 0.
    endif

    allocate ( t_dat(is:ie,js:je,npz,2) )
    allocate ( q_dat(is:ie,js:je,npz,2) )
    t_dat = 0.
    q_dat = 0.

    allocate ( s2c(is:ie,js:je,4) )
    allocate ( id1(is:ie,js:je) )
    allocate ( id2(is:ie,js:je) )
    allocate ( jdc(is:ie,js:je) )

! Initialize remapping coefficients:

    call field_size(file_names(1), 'T', tsize, field_found=found)

    if ( found ) then
         im = tsize(1); jm = tsize(2); km = tsize(3)
         if(gid==0)  write(*,*) 'NCEP analysis dimensions:', tsize
    else
         call mpp_error(FATAL,'==> Error from get_ncep_analysis:        &
                                   T field not found')
    endif

    allocate (  lon(im) )
    allocate (  lat(jm) )

    call read_data (file_names(1), 'LAT', lat, no_domain=.true.)
    call read_data (file_names(1), 'LON', lon, no_domain=.true.)

! Convert to radian
    do i=1,im
       lon(i) = lon(i) * deg2rad ! lon(1) = 0.
    enddo
    do j=1,jm
       lat(j) = lat(j) * deg2rad
    enddo
 
    allocate ( ak0(km+1) )
    allocate ( bk0(km+1) )

    call read_data (file_names(1), 'hyai', ak0, no_domain=.true.)
    call read_data (file_names(1), 'hybi', bk0, no_domain=.true.)

! Note: definition of NCEP hybrid is p(k) = a(k)*1.E5 + b(k)*ps
    ak0(:) = ak0(:) * 1.E5

! Limiter to prevent NAN at top during remapping 
    ak0(1) = max(1.e-8, ak0(1))

!-----------------------------------------------------------
! Initialize lat-lon to Cubed bi-linear interpolation coeff:
!-----------------------------------------------------------
    call remap_coef

! Get first dataset
    nt = 2
    call get_ncep_analysis ( ps_dat(:,:,nt), u_dat(:,:,:,nt), v_dat(:,:,:,nt),     &
                             t_dat(:,:,:,nt), q_dat(:,:,:,nt), npz, zvir, ak, bk,  &
                             ts, phis, nfile, file_names(nfile) )


    module_is_initialized = .true.
    
 end subroutine fv_nwp_nudge_init


 subroutine get_ncep_analysis ( ps, u, v, t, q, npz, zvir, ak, bk, ts, phis, nfile, fname )
  integer,  intent(in):: npz           ! vertical dimension
  real,     intent(in):: zvir
  real, intent(in), dimension(npz+1):: ak, bk
  real, intent(in), dimension(isd:ied,jsd:jed):: phis
  character(len=128), intent(in):: fname
  integer,  intent(inout):: nfile
!
  real, intent(out), dimension(is:ie,js:je):: ts, ps
  real, intent(out), dimension(is:ie,js:je,npz):: u, v, t, q
! local:
  real, allocatable:: oro(:,:), wk2(:,:), wk3(:,:,:)
  real, allocatable:: tp(:,:,:), qp(:,:,:)
  real, allocatable:: ua(:,:,:), va(:,:,:)
  real psc(is:ie,js:je)
  real gzc(is:ie,js:je)
  real tmean
  integer:: i, j, k, npt
  integer:: i1, i2, j1
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
     call read_data (fname, 'PS', wk2, no_domain=.true.)
     if(gid==0) call pmaxmin( 'PS_ncep', wk2, im,  jm, 0.01)

     do j=js,je
        do i=is,ie
           i1 = id1(i,j)
           i2 = id2(i,j)
           j1 = jdc(i,j)
           psc(i,j) = s2c(i,j,1)*wk2(i1,j1  ) + s2c(i,j,2)*wk2(i2,j1  ) +  &
                      s2c(i,j,3)*wk2(i2,j1+1) + s2c(i,j,4)*wk2(i1,j1+1)
        enddo
     enddo

     call read_data (fname, 'PHIS', wk2, no_domain=.true.)
!    if(gid==0) call pmaxmin( 'ZS_ncep', wk2, im,  jm, 1./grav)
     do j=js,je
        do i=is,ie
           i1 = id1(i,j)
           i2 = id2(i,j)
           j1 = jdc(i,j)
           gzc(i,j) = s2c(i,j,1)*wk2(i1,j1  ) + s2c(i,j,2)*wk2(i2,j1  ) +  &
                      s2c(i,j,3)*wk2(i2,j1+1) + s2c(i,j,4)*wk2(i1,j1+1)
        enddo
     enddo

     if ( read_ts ) then       ! read skin temperature; could be used for SST

      call read_data (fname, 'TS', wk2, no_domain=.true.)

      if ( .not. land_ts ) then
           allocate ( oro(im,jm) )
! Read NCEP ORO (1; land; 0: ocean; 2: sea_ice)
           call read_data (fname, 'ORO', oro, no_domain=.true.)

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

!     if(gid==0) call pmaxmin('SST_ncep', wk2, im,  jm, 1.)
      do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            ts(i,j) = s2c(i,j,1)*wk2(i1,j1  ) + s2c(i,j,2)*wk2(i2,j1  ) +  &
                      s2c(i,j,3)*wk2(i2,j1+1) + s2c(i,j,4)*wk2(i1,j1+1)
         enddo
      enddo
!     call prt_maxmin('SST_model', ts, is, ie, js, je, 0, 1, 1., master)

! Perform interp to FMS SST format/grid
      call ncep2fms( wk2 )
      if(gid==0) call pmaxmin( 'SST_ncep_fms',  sst_ncep, i_sst, j_sst, 1.)

      endif     ! read_ts

      deallocate ( wk2 ) 

! Read in temperature:
      allocate (  wk3(im,jm,km) )
      call read_data (fname, 'T',  wk3, no_domain=.true.)
      if(gid==0) call pmaxmin( 'T_ncep',   wk3, im*jm, km, 1.)
      allocate (  tp(is:ie,js:je,km) )
      do k=1,km
      do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            tp(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                        s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
         enddo
      enddo
      enddo

! Read in tracers: only sphum at this point
      call read_data (fname, 'Q', wk3, no_domain=.true.)
      if(gid==1) call pmaxmin( 'Q_ncep',   wk3, im*jm, km, 1.)
      allocate ( qp(is:ie,js:je,km) )
      do k=1,km
      do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            qp(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                        s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
         enddo
      enddo
      enddo

      call remap_tq(npz, 1, 1, ak, bk, ak0, bk0,   &
                    psc, gzc, tp, qp, ps, phis, t, q, zvir)
      deallocate ( tp ) 
      deallocate ( qp ) 

! Winds:
   if ( nudge_winds ) then

      call read_data (fname, 'U',  wk3, no_domain=.true.)
      if( master ) call pmaxmin( 'U_ncep',   wk3, im*jm, km, 1.)
      allocate ( ua(is:ie,js:je,km) )

      do k=1,km
      do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            ua(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                        s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
         enddo
      enddo
      enddo

      call read_data (fname, 'V',  wk3, no_domain=.true.)
      if( master ) call pmaxmin( 'V_ncep',  wk3, im*jm, km, 1.)
      allocate ( va(is:ie,js:je,km) )
      do k=1,km
      do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            va(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                        s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
         enddo
      enddo
      enddo

      call remap_uv(npz, ak, bk, ak0, bk0, ps, psc, ua, va, u, v)

      deallocate ( ua ) 
      deallocate ( va ) 
   endif

   deallocate ( wk3 ) 

  nfile = nfile + 1

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

 subroutine remap_tq(npz, nq, ncnst, ak, bk, ak0, bk0,     &
                     psc, gzc, ta, qa, ps, phis, t, q, zvir)
  integer, intent(in):: npz, nq, ncnst
  real,    intent(in):: zvir
  real,    intent(in):: ak(npz+1), bk(npz+1)
  real,    intent(in):: ak0(km+1), bk0(km+1)
  real,    intent(in), dimension(is:ie,js:je):: psc, gzc 
  real,    intent(in), dimension(isd:ied,jsd:jed):: phis
  real,    intent(in), dimension(is:ie,js:je,km):: ta
  real,    intent(in), dimension(is:ie,js:je,km,ncnst):: qa
!
  real,    intent(out), dimension(is:ie,js:je):: ps
  real,    intent(out), dimension(is:ie,js:je,npz):: t
  real,    intent(out), dimension(is:ie,js:je,npz,ncnst):: q
! local:
  real, dimension(is:ie,km):: tp
  real, dimension(is:ie,km+1):: pe0, pn0
  real, dimension(is:ie,npz):: qn1
  real, dimension(is:ie,npz+1):: pe1, pn1
  real pt0(km), gz(km+1), pk0(km+1)
  real qp(is:ie,km,ncnst)
  real pst
  integer i,j,k, iq
  integer  sphum

  sphum   = 1

  pk0(1) = ak0(1)**kappa 

  do 5000 j=js,je

     do i=is,ie
        pe0(i,1) = ak0(1)
        pn0(i,1) = log(ak0(1))
     enddo

     do i=is,ie

       do iq=1,ncnst
          do k=1,km
             qp(i,k,iq) = qa(i,j,k,iq)
          enddo
       enddo

       do k=1,km
          tp(i,k) = ta(i,j,k)*(1.+zvir*qp(i,k,sphum))
       enddo
! Tracers:

       do k=2,km+1
          pe0(i,k) = ak0(k) + bk0(k)*psc(i,j)
          pn0(i,k) = log(pe0(i,k))
          pk0(k) = pe0(i,k)**kappa
       enddo

! * Adjust interpolated ps to model terrain
       gz(km+1) = gzc(i,j) 
       do k=km,1,-1
           gz(k) = gz(k+1) + rdgas*tp(i,k)*(pn0(i,k+1)-pn0(i,k)) 
       enddo
! Only lowest layer potential temp is needed
          pt0(km) = tp(i,km)/(pk0(km+1)-pk0(km))*(kappa*(pn0(i,km+1)-pn0(i,km)))
       if( phis(i,j)>gzc(i,j) ) then
           do k=km,1,-1
              if( phis(i,j) <  gz(k)  .and.    &
                  phis(i,j) >= gz(k+1) ) then
                  pst = pk0(k) + (pk0(k+1)-pk0(k))*(gz(k)-phis(i,j))/(gz(k)-gz(k+1))
                  go to 123
              endif
           enddo
       else
! Extrapolation into the ground
           pst = pk0(km+1) + (gzc(i,j)-phis(i,j))/(cp_air*pt0(km))
       endif

123    ps(i,j) = pst**(1./kappa)
     enddo   !i-loop
 
   if ( nudge_t .or. nudge_virt .or. nudge_q ) then

     do i=is,ie
        pe1(i,1) = ak(1)
        pn1(i,1) = log(pe1(i,1))
     enddo
     do k=2,npz+1
       do i=is,ie
          pe1(i,k) = ak(k) + bk(k)*ps(i,j)
          pn1(i,k) = log(pe1(i,k))
       enddo
     enddo

   if ( nudge_t ) then
!---------------
! map tracers
!----------------
      do iq=1,ncnst
         call mappm(km, pe0, qp(is,1,iq), npz, pe1,  qn1, is,ie, 0, 8)
         do k=1,npz
            do i=is,ie
               q(i,j,k,iq) = qn1(i,k)
            enddo
         enddo
      enddo
      call mappm(km, pn0, tp, npz, pn1, qn1, is,ie, 1, 8)
      do k=1,npz
         do i=is,ie
            t(i,j,k) = qn1(i,k)/(1.+zvir*q(i,j,k,sphum))
         enddo
      enddo
   elseif ( nudge_virt ) then
!-------------------------------------------------------------
! map virtual temperature using geopotential conserving scheme.
!-------------------------------------------------------------
      call mappm(km, pn0, tp, npz, pn1, qn1, is,ie, 1, 8)
      do k=1,npz
         do i=is,ie
            t(i,j,k) = qn1(i,k)
         enddo
      enddo
   endif

   endif
5000 continue

! call prt_maxmin('PS_model', ps, is, ie, js, je, ng, 1, 0.01, gid==0)

  if (gid==0) write(*,*) 'done remap_tq'

 end subroutine remap_tq


 subroutine remap_uv(npz, ak, bk, ak0, bk0, ps, psc, ua, va, u, v)
  integer, intent(in):: npz
  real,    intent(in):: ak(npz+1), bk(npz+1)
  real,    intent(in):: ak0(km+1), bk0(km+1)
  real,    intent(in):: ps(is:ie,js:je)
  real,    intent(in):: psc(is:ie,js:je)
  real,    intent(in), dimension(is:ie,js:je,km):: ua, va
!
  real,    intent(out), dimension(is:ie,js:je,npz):: u, v
! local:
  real, dimension(is:ie, km+1):: pe0
  real, dimension(is:ie,npz+1):: pe1
  real, dimension(is:ie,npz):: qn1
  integer i,j,k

  do 5000 j=js,je

     do i=is,ie
        pe0(i,1) = ak0(1)
        pe1(i,1) = ak(1)
     enddo

     do k=2,km+1
        do i=is,ie
           pe0(i,k) = ak0(k) + bk0(k)*psc(i,j)
        enddo
     enddo

     do k=2,npz+1
       do i=is,ie
          pe1(i,k) = ak(k) + bk(k)*ps(i,j)
       enddo
     enddo

!------
! map u
!------
      call mappm(km, pe0, ua(is:ie,j,1:km), npz, pe1, qn1, is,ie, -1, 8)
      do k=1,npz
         do i=is,ie
            u(i,j,k) = qn1(i,k)
         enddo
      enddo
!------
! map v
!------
      call mappm(km, pe0, va(is:ie,j,1:km), npz, pe1, qn1, is,ie, -1, 8)
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

 end subroutine fv_nwp_nudge_end

 subroutine pmaxmin( qname, a, imax, jmax, fac )

      character*(*)  qname
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


end module fv_nwp_nudge_mod
