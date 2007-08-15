module hswf_mod

      use constants_mod,   only: grav, rdgas, RADIAN, kappa
      use mpp_domains_mod, only: mpp_update_domains
      use mp_mod,          only: domain
      use fv_sg_mod,       only: fv_sg_conv, qsmith

      use fv_diagnostics_mod, only: id_prec
      use diag_manager_mod,   only: send_data
      use time_manager_mod,   only: time_type
#ifdef MARS_GCM
      use fms_mod, only: file_exist, read_data, field_size
      use mpp_mod, only: mpp_error, FATAL
      use horiz_interp_mod,  only: horiz_interp
#endif

implicit none
!-----------------------------------------------------------------------
      logical :: rf_initialized = .false.

#ifdef MARS_GCM
      logical :: tmars_initialized = .false.
      real,  allocatable, dimension(:,:,:):: tmars
#endif


private

public :: Held_Suarez_Strat, Held_Suarez_Tend, Sim_phys, age_of_air


contains

!-----------------------------------------------------------------------
 subroutine Held_Suarez_Tend(npx, npy, npz, is, ie, js, je, ng, nq,   &
                              u, v, pt, q, pe, delp, peln, pkz, pdt,  &
                              ua, va, u_dt, v_dt, t_dt, q_dt, agrid,  &
                              delz, hydrostatic, ak, bk, ks,   &
                              strat, rayf, master, Time, time_total)

      integer, INTENT(IN   ) :: npx, npy, npz
      integer, INTENT(IN   ) :: is, ie, js, je, ng, nq
      logical, intent(IN)    :: hydrostatic
      real   , INTENT(IN   ) ::  delz(is:ie,js:je,npz)

      real   , INTENT(INOUT) ::    u(is-ng:ie+  ng,js-ng:je+1+ng,npz)
      real   , INTENT(INOUT) ::    v(is-ng:ie+1+ng,js-ng:je+  ng,npz)
      real   , INTENT(INOUT) ::   pt(is-ng:ie+  ng,js-ng:je+  ng,npz)
      real   , INTENT(INOUT) :: delp(is-ng:ie+  ng,js-ng:je+  ng,npz)
      real   , INTENT(INOUT) ::    q(is-ng:ie+  ng,js-ng:je+  ng,npz, nq)
      real   , INTENT(INOUT) ::   pe(is-1:ie+1 ,1:npz+1,js-1:je+1)
      real   , INTENT(INOUT) :: peln(is  :ie     ,1:npz+1,js   :je     )
      real   , INTENT(INOUT) ::  pkz(is  :ie     ,js   :je     ,1:npz)

      real   , INTENT(INOUT) ::   ua(is-ng:ie+ng,js-ng:je+ng,npz)
      real   , INTENT(INOUT) ::   va(is-ng:ie+ng,js-ng:je+ng,npz)

! Tendencies:
      real, INTENT(INOUT):: u_dt(is-ng:ie+ng,js-ng:je+ng,npz)
      real, INTENT(INOUT):: v_dt(is-ng:ie+ng,js-ng:je+ng,npz)
      real, INTENT(INOUT):: t_dt(is:ie,js:je,npz)
      real, INTENT(INOUT):: q_dt(is:ie,js:je,npz,nq)


      real   , INTENT(IN   ) :: agrid(is-ng:ie+ng,js-ng:je+ng, 2)
      real   , INTENT(IN   ) :: ak(npz+1), bk(npz+1)
      integer, INTENT(IN   ) :: ks

      real   , INTENT(IN   ) :: pdt
      logical, INTENT(IN   ) :: strat, rayf, master

      type(time_type), intent(in) :: Time
      real, INTENT(IN), optional:: time_total

! Local
      real pref(npz)
      integer  i,j,k
      real  ty, tz, akap 
      real  p0, t0, sday, rkv, rka, rks, rkt, sigb, rsgb
      real  tmp
      real  ap0k, algpk
      real  tey, tez, fac, pw, sigl
      real  h0, dz
      real  dt_tropic
      real  rmr, rms
      real  relx, tau
      real  t_st, t_ms
      real  rdt, f1
      real  pc, c1

      real, allocatable, SAVE ::  rf(:)
      real :: frac(is-ng:ie+ng,js-ng:je+ng)
      real ::  pl(is:ie, js:je, 1:npz)
      real :: teq(is:ie, js:je, 1:npz)
      real rdg

      rdg = -rdgas / grav

      ty = 60.0
      tz = 10.0             ! Original value from H-S was 10.
      akap = 2./7.

      p0 = 100000.
      t0 = 200.
      h0 = 7.
      sday = 24.*3600.
      rdt = 1. / pdt

      rkv = pdt/sday
      rka = pdt/ (40.*sday)      ! was 40 days
      rks = pdt/ (4.*sday)       ! was 4 days

! For strat-mesosphere
      t_ms = 10.
      t_st = 40.
      tau = (t_st - t_ms) / log(100.)
      rms = pdt/(t_ms*sday)
      rmr =  1./(1.+rms)

      sigb = 0.7
      rsgb = 1./(1.-sigb)
      ap0k = 1./p0**akap
      algpk = log(ap0k)

      do k=1,npz
         pref(k) = ak(k) + bk(k)*1.E5
      enddo

! Setup...
      if ( rayf .and. (.not. rf_initialized) ) then
          allocate( rf(npz) )
          c1 = 1. / (36.*3600)
          pc = 1.
          if(master) write(6,*) 'HSWF Forcing ...' 
          do k=1,ks
             tmp = (ak(k+1)-ak(k))/log(ak(k+1)/ak(k))
             rf(k) = c1*(1.+tanh(log10(pc/tmp)))
             if(master) write(6,*) k, 0.01*tmp, 1./(rf(k)*sday)
             rf(k) = 1./(1.+pdt*rf(k))
          enddo
          if(master) write(6,*) ' '
          rf_initialized = .true.
      endif

        do k=1,npz
           do j=js,je
              do i=is,ie
                 pl(i,j,k) = delp(i,j,k) / ( peln(i,k+1,j)-peln(i,k,j))
              enddo
           enddo 
        enddo

        if ( .not. hydrostatic ) then
           do k=1,npz
              do j=js,je
                 do i=is,ie
                    pkz(i,j,k) = pl(i,j,k)**akap
                 enddo
              enddo 
           enddo
         endif
        
! Temperature forcing...
        do k=npz,1,-1
           do j=js,je
              do i=is,ie
                 tey = ap0k*( 315.0 - ty*SIN(agrid(i,j,2))*SIN(agrid(i,j,2)) )
                 tez = tz*( ap0k/akap )*COS(agrid(i,j,2))*COS(agrid(i,j,2)) 

                 if (strat .and. pl(i,j,k) < 10000.    &
                           .and. pl(i,j,k) > 100.  )  then
                   dz = h0 * log(pl(i,j,k+1)/pl(i,j,k))
!
! Lapse rate above tropic stratopause is 2.25 deg/km
! Relaxation time is t_st days at 100 mb (as H-S) and gradually
! decreases to t_ms Days at and above the stratopause
!
                   relx =  t_ms + tau*log(0.01*pl(i,j,k))
                   relx = pdt/(relx*sday)
                   dt_tropic = 2.25*COS(agrid(i,j,2)) * dz
                   teq(i,j,k)  =  teq(i,j,k+1) + dt_tropic
                   t_dt(i,j,k) = relx*(teq(i,j,k)-pt(i,j,k))/(1.+relx) * rdt
                 elseif (strat .and. pl(i,j,k) <= 100.)  then
!
! Mesosphere: defined as the region above 1 mb
!
                   dz = h0 * log(pl(i,j,k+1)/pl(i,j,k))
                   dt_tropic = -2.25*COS(agrid(i,j,2)) * dz
                   teq(i,j,k) = teq(i,j,k+1) + dt_tropic
                   t_dt(i,j,k) = ((pt(i,j,k)+rms*teq(i,j,k))*rmr - pt(i,j,k))*rdt
                 else

! Trop:  strictly Held-Suarez

                   sigl = pl(i,j,k)/pe(i,npz+1,j)
                   f1 = max(0., (sigl-sigb) * rsgb )
                   teq(i,j,k) = tey - tez*(log(pkz(i,j,k))+algpk)
                   teq(i,j,k) = max(t0, teq(i,j,k)*pkz(i,j,k))
                   rkt = rka + (rks-rka)*f1*(COS(agrid(i,j,2))**4.0)
                   t_dt(i,j,k) = rkt*(teq(i,j,k)-pt(i,j,k))/(1.+rkt) * rdt
                 endif
              enddo     !i-loop
           enddo     !j-loop
        enddo     !k-loop

        if ( .not. hydrostatic ) then
!------------------------------------------------------------------------
! The nonhydrostatic pressure changes if there is heating (under constant
! volume and mass is locally conserved).
!------------------------------------------------------------------------
           do k=1,npz
              do j=js,je
                 do i=is,ie
                      pkz(i,j,k) = (rdg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k))**akap
                 enddo
              enddo 
           enddo
        endif


! Velocity dissipation damping
      do 2000 k=1,npz

        if (rayf .and. k <= ks) then
! Apply Rayleigh friction
          do j=js,je
             do i=is,ie
                u_dt(i,j,k) = ua(i,j,k)*(rf(k) - 1.) * rdt
                v_dt(i,j,k) = va(i,j,k)*(rf(k) - 1.) * rdt
             enddo
          enddo
        else
! Surface Rayleigh friction according to Held-Suarez
          do j=js,je
             do i=is,ie
                sigl = 0.5*(pe(i,k,j)+pe(i,k+1,j)) / pe(i,npz+1,j)
                frac(i,j) = max(0., (sigl-sigb)*rsgb ) * rkv
                if (frac(i,j) > 0.) then
                    u_dt(i,j,k) = -ua(i,j,k)*frac(i,j)/(1.+frac(i,j)) * rdt
                    v_dt(i,j,k) = -va(i,j,k)*frac(i,j)/(1.+frac(i,j)) * rdt
                endif
             enddo
          enddo
        endif

2000  continue

      if( nq/=0 )     &
      call age_of_air(is, ie, js, je, npz, ng, time_total, pe, q(is-ng,js-ng,1,nq))


 end subroutine Held_Suarez_Tend

 subroutine Held_Suarez_Strat(npx, npy, npz, is, ie, js, je, ng, nq,  &
                              u, v, pt, q, pe, delp, peln, pkz, pdt,  &
                              ua, va, agrid, ak, bk, ks, strat,  &
                              rayf, master, Time, time_total)

      integer, INTENT(IN   ) :: npx, npy, npz
      integer, INTENT(IN   ) :: is, ie, js, je, ng, nq

      real   , INTENT(INOUT) ::    u(is-ng:ie+  ng,js-ng:je+1+ng,npz)
      real   , INTENT(INOUT) ::    v(is-ng:ie+1+ng,js-ng:je+  ng,npz)
      real   , INTENT(INOUT) ::   pt(is-ng:ie+  ng,js-ng:je+  ng,npz)
      real   , INTENT(INOUT) :: delp(is-ng:ie+  ng,js-ng:je+  ng,npz)
      real   , INTENT(INOUT) ::    q(is-ng:ie+  ng,js-ng:je+  ng,npz, nq)
      real   , INTENT(INOUT) ::   pe(is-1:ie+1 ,1:npz+1,js-1:je+1)
      real   , INTENT(INOUT) :: peln(is   :ie     ,1:npz+1,js   :je     )
      real   , INTENT(INOUT) ::  pkz(is   :ie     ,js   :je     ,1:npz)

      real   , INTENT(INOUT) ::   ua(is-ng:ie+ng,js-ng:je+ng,npz)
      real   , INTENT(INOUT) ::   va(is-ng:ie+ng,js-ng:je+ng,npz)

      real   , INTENT(IN   ) :: agrid(is-ng:ie+ng,js-ng:je+ng, 2)
      real   , INTENT(IN   ) :: ak(npz+1), bk(npz+1)
      integer, INTENT(IN   ) :: ks

      real   , INTENT(IN   ) :: pdt
      logical, INTENT(IN   ) :: strat, rayf, master

      type(time_type), intent(in) :: Time
      real, INTENT(IN), optional:: time_total

! Local
      real pref(npz)
      integer  i,j,k
      real  ty, tz, akap 
      real  p0, t0, sday, rkv, rka, rks, rkt, sigb, rsgb
      real  tmp
      real  ap0k, algpk
      real  tey, tez, fac, pw, sigl
      real  h0, dz
      real  dt_tropic
      real  rmr, rms
      real  relx, tau
      real  t_st, t_ms
      real  f1
      real  pc, c1

      real, allocatable, SAVE ::  rf(:)
      real :: frac(is-ng:ie+ng,js-ng:je+ng)
      real ::  pl(is:ie, js:je, 1:npz)
      real :: teq(is:ie, js:je, 1:npz)

      ty = 60.0
      tz = 10.0             ! Original value from H-S was 10.
      akap = 2./7.

      p0 = 100000.
      t0 = 200.
      h0 = 7.
      sday = 24.*3600.

      rkv = 0.5*pdt/sday
      rka = pdt/ (40.*sday)      ! was 40 days
      rks = pdt/ (4.*sday)       ! was 4 days

! For strat-mesosphere
      t_ms = 10.
      t_st = 40.
      tau = (t_st - t_ms) / log(100.)
      rms = pdt/(t_ms*sday)
      rmr =  1./(1.+rms)

      sigb = 0.7
      rsgb = 1./(1.-sigb)
      ap0k = 1./p0**akap
      algpk = log(ap0k)

      do k=1,npz
         pref(k) = ak(k) + bk(k)*1.E5
      enddo

! Setup...
        if ( rayf .and. (.not. rf_initialized) ) then
          allocate( rf(npz) )
          c1 = 1. / (36.*3600)
          pc = 1.
          if(master) write(6,*) 'HSWF Forcing ...' 
          do k=1,ks
             tmp = (ak(k+1)-ak(k))/log(ak(k+1)/ak(k))
             rf(k) = c1*(1.+tanh(log10(pc/tmp)))
             if(master) write(6,*) k, 0.01*tmp, 1./(rf(k)*sday)
             rf(k) = 1./(1.+pdt*rf(k))
          enddo
          if(master) write(6,*) ' '
          rf_initialized = .true.
        endif

        do k=1,npz
           do j=js,je
              do i=is,ie
                 pl(i,j,k) = delp(i,j,k) / ( peln(i,k+1,j)-peln(i,k,j))
              enddo
           enddo 
        enddo

! Temperature forcing...
        do k=npz,1,-1
           do j=js,je
              do i=is,ie
                 tey = ap0k*( 315.0 - ty*SIN(agrid(i,j,2))*SIN(agrid(i,j,2)) )
                 tez = tz*( ap0k/akap )*COS(agrid(i,j,2))*COS(agrid(i,j,2)) 

                 if (strat .and. pl(i,j,k) < 10000.    &
                           .and. pl(i,j,k) > 100.  )  then
                   dz = h0 * log(pl(i,j,k+1)/pl(i,j,k))
!
! Lapse rate above tropic stratopause is 2.25 deg/km
! Relaxation time is t_st days at 100 mb (as H-S) and gradually
! decreases to t_ms Days at and above the stratopause
!
                   relx =  t_ms + tau*log(0.01*pl(i,j,k))
                   relx = pdt/(relx*sday)
                   dt_tropic = 2.25*COS(agrid(i,j,2)) * dz
                   teq(i,j,k)  =  teq(i,j,k+1) + dt_tropic
                   pt(i,j,k) = (pt(i,j,k)+relx*teq(i,j,k))/(1.+relx)
                 elseif (strat .and. pl(i,j,k) <= 100.)  then
!
! Mesosphere: defined as the region above 1 mb
!
                   dz = h0 * log(pl(i,j,k+1)/pl(i,j,k))
                   dt_tropic = -2.25*COS(agrid(i,j,2)) * dz
                   teq(i,j,k) = teq(i,j,k+1) + dt_tropic
                   pt(i,j,k) = (pt(i,j,k)+rms*teq(i,j,k))*rmr
                 else

! Trop:  strictly Held-Suarez

                   sigl = pl(i,j,k)/pe(i,npz+1,j)
                   f1 = max(0., (sigl-sigb) * rsgb )
                   teq(i,j,k) = tey - tez*(log(pkz(i,j,k))+algpk)
                   teq(i,j,k) = max(t0, teq(i,j,k)*pkz(i,j,k))
                   rkt = rka + (rks-rka)*f1*(COS(agrid(i,j,2))**4.0)
                   pt(i,j,k) = (pt(i,j,k)+rkt*teq(i,j,k))/(1.+rkt)
                 endif
              enddo     !i-loop
           enddo     !j-loop
        enddo     !k-loop

! Velocity dissipation damping
      do 2000 k=1,npz

        if (rayf .and. k <= ks) then

! Apply Rayleigh friction
          do j=js,je+1
            do i=is,ie
              u(i,j,k) = u(i,j,k)*rf(k)
            enddo
          enddo
          do j=js,je
            do i=is,ie+1
              v(i,j,k) = v(i,j,k)*rf(k)
            enddo
          enddo

        else

! Surface Rayleigh friction according to Held-Suarez
          do j=js,je
            do i=is,ie
              sigl = 0.5*(pe(i,k,j)+pe(i,k+1,j)) / pe(i,npz+1,j)
              frac(i,j) = max(0., (sigl-sigb)*rsgb )
            enddo
          enddo
#if defined(SPMD)
          call mpp_update_domains( frac, domain )
#endif
! Backward adjustment
          do j=js,je+1
            do i=is,ie+1
              fac = frac(i,j)+frac(i,j-1)
              if (fac > 0.) then
                u(i,j,k) = u(i,j,k)/(1.+rkv*fac)
              endif
            enddo
          enddo
          do j=js,je
            do i=is,ie+1
              fac = frac(i,j)+frac(i-1,j)
              if (fac > 0.) then
                v(i,j,k) = v(i,j,k)/(1.+rkv*fac)
              endif
            enddo
          enddo

        endif

2000  continue

      if( nq/=0 )     &
      call age_of_air(is, ie, js, je, npz, ng, time_total, pe, q(is-ng,js-ng,1,nq))


 end subroutine Held_Suarez_Strat


#ifdef MARS_GCM
! Forcing for MARS_GCM
 subroutine Sim_phys(npx, npy, npz, is, ie, js, je, ng, nq,  &
                     u_dt, v_dt, t_dt, q_dt, u, v, ua, va,   &
                     pt, q, pe, delp, peln, oro, hydrostatic, &
                     pdt, agrid, ak, bk, rayf, p_ref, master, Time, time_total)

      integer, INTENT(IN) :: npx, npy, npz
      integer, INTENT(IN) :: is, ie, js, je, ng, nq
      real   , INTENT(IN) :: pdt
      real   , INTENT(IN) :: agrid(is-ng:ie+ng,js-ng:je+ng, 2)
      real   , INTENT(IN) :: ak(npz+1), bk(npz+1)
      logical, INTENT(IN) :: rayf, master
      real, INTENT(IN):: oro(is:ie,js:je)       ! land fraction
      logical, INTENT(IN):: hydrostatic
      real, INTENT(IN):: p_ref

      type(time_type), intent(in) :: Time
      real, INTENT(IN), optional:: time_total

      real   , INTENT(INOUT) :: u(is-ng:ie+  ng,js-ng:je+1+ng,npz)
      real   , INTENT(INOUT) :: v(is-ng:ie+1+ng,js-ng:je+  ng,npz)

      real, INTENT(INOUT)::   pt(is-ng:ie+ng,js-ng:je+ng,npz)
      real, INTENT(INOUT):: delp(is-ng:ie+ng,js-ng:je+ng,npz)
      real, INTENT(INOUT)::    q(is-ng:ie+ng,js-ng:je+ng,npz, nq)
      real, INTENT(INOUT)::   pe(is-1:ie+1 ,1:npz+1,js-1:je+1)
      real, INTENT(INOUT):: peln(is  :ie   ,1:npz+1,js  :je  )

! Tendencies:
      real, INTENT(INOUT):: u_dt(is-ng:ie+ng,js-ng:je+ng,npz)
      real, INTENT(INOUT):: v_dt(is-ng:ie+ng,js-ng:je+ng,npz)
      real, INTENT(INOUT):: t_dt(is:ie,js:je,npz)
      real, INTENT(INOUT):: q_dt(is:ie,js:je,npz,nq)

      real, INTENT(INOUT):: ua(is-ng:ie+ng,js-ng:je+ng,npz)
      real, INTENT(INOUT):: va(is-ng:ie+ng,js-ng:je+ng,npz)

! Local
      real pedge(npz+1)
      real pref(npz)

      real  sday, rkv, rkt, sigb, rsgb, sigl, cs
      real frac
      real rdt                    ! Reciprocal of dt
      character (len=128) :: filename
      integer  i,j,k

      sday = 24.*3600.

      do k=1,npz+1
         pedge(k) = ak(k) + bk(k)*p_ref
      enddo

      do k=1,npz
         pref(k) = (pedge(k+1)-pedge(k)) / log(pedge(k+1)/pedge(k))
      enddo

      if ( .not. tmars_initialized ) then
          allocate( tmars(is:ie,js:je,npz) )

          filename= 'INPUT/teq.nc' 
          if( file_exist( trim( filename ) ) ) then 
              call read_teq( filename, npz, agrid(is:ie,js:je,1), agrid(is:ie,js:je,2),  &
                             pref, tmars(is:ie,js:je,:)  )
              if(master) write(6,*) 'TEQ for Mars initialized.'
          else
#ifdef FAIL_SAFE
              call mpp_error(FATAL,'Mars_GCM: TEQ data not found')
#else
              do k=1,npz
                 do j=js,je
                    do i=is,ie
                       tmars(i,j,k) = 100.+ 100.*max(0.25, (1.-sin(agrid(i,j,2)))*pref(k)/pedge(npz+1))
!                      tmars(i,j,k) = 120. + 100.*    &
!                      max(0.2, (1.-sin(agrid(i,j,2)))*log(pref(k))/log(pedge(npz+1)))
                    enddo 
                 enddo 
              enddo
              if(master) write(6,*) 'Data for Mars not found; using analytic profile'
#endif
          endif 
          cs = sqrt( rdgas * 273./(1.-kappa) ) 
          if(master) write(6,*) 'Sound speed (T=273)=', cs
          tmars_initialized = .true.
      endif

! ***  Newtonian cooling/relaxation
      rdt = 1. / pdt
!     rkt = pdt / (4.*sday)
      rkt = pdt / (8.*sday)

      do k=1,npz
         do j=js,je
            do i=is,ie
               t_dt(i,j,k) = rkt*(tmars(i,j,k)-pt(i,j,k))/(1.+rkt) * rdt
            enddo 
         enddo 
      enddo

! *** Surface Rayleigh friction according to Held-Suarez
      rkv = pdt / (1.*sday)          ! 1 day
!     rkv = pdt / (2.*sday)          ! 2 day
      sigb = 0.7
      rsgb = 1./(1.-sigb)

      do k=1,npz
          do j=js,je
             do i=is,ie
                sigl = 0.5*(pe(i,k,j)+pe(i,k+1,j)) / pe(i,npz+1,j)
                frac = rkv * (sigl-sigb)*rsgb
                if (frac > 0.) then
                    u_dt(i,j,k) = -ua(i,j,k)*frac/(1.+frac) * rdt
                    v_dt(i,j,k) = -va(i,j,k)*frac/(1.+frac) * rdt
                endif
             enddo
          enddo
      enddo



 end subroutine Sim_phys


 subroutine read_teq ( filename, nlevels, lon, lat, pstd, tout )

!-----------------------------------------------------------------------
!
! routine for initializing the radiative-convective temperature cross-section
!
!-----------------------------------------------------------------------



   character(len=128), intent(in)            :: filename
   integer, intent(in):: nlevels
   real,    intent(in),  dimension(:,:)      :: lon, lat 
   real,    intent(in)       ::  pstd(nlevels)

   real,    intent(out),  dimension(:,:,:)   ::  tout

!-----------------------------------------------------------------------
   integer  unit, io, ierr, i, j, k, klev, kk

   integer  im, jm, km, fld_dims(4)
   real    ::  frac

   real, dimension(:,:,:),  allocatable  ::   teq_inpt
   real, dimension(:,:),    allocatable  ::   txy

   real, dimension(:),  allocatable  ::   lat_inpt, pres_inpt,  presh_inpt
   real, dimension(:),  allocatable  ::   lonb_inpt, latb_inpt 



!       Get input field teq
   call field_size( trim(filename), 'lat', fld_dims )

   allocate( lat_inpt (fld_dims(1)  ) )
   allocate( latb_inpt(fld_dims(1)+1) )

   call read_data( trim(filename), 'lat',  lat_inpt,  no_domain=.true. )
   call read_data( trim(filename), 'latb', latb_inpt, no_domain=.true. )

   call field_size( trim(filename), 'lonb', fld_dims )

   allocate( lonb_inpt (fld_dims(1)  ) )

   call read_data( trim(filename), 'lonb', lonb_inpt, no_domain=.true. )


   call field_size( trim(filename), 'pfull', fld_dims )
   allocate( pres_inpt( fld_dims(1) ) )
   call read_data( trim(filename), 'pfull', pres_inpt, no_domain=.true. )
   print *, 'pfull dims:  ', fld_dims 
   pres_inpt= 100.0 * pres_inpt



   call field_size( trim(filename), 'teq', fld_dims )

   im= fld_dims(1);  jm= fld_dims(2);  km= fld_dims(3) 
       print *, 'Input Teq dims:  ', fld_dims 

   allocate( teq_inpt( im,jm,km ) )
   allocate( txy     ( im,jm )    )

   call read_data( trim(filename), 'teq', teq_inpt, no_domain=.true. )

   latb_inpt(:)= latb_inpt(:)/RADIAN
   lonb_inpt(:)= lonb_inpt(:)/RADIAN


!           If km != nlevels  then require vertical interpolation 
  if( nlevels > km .or. nlevels < km )  then

     DO k= 1, nlevels
        if( pres_inpt(1) > pstd(k) ) then
           klev= 2;  frac= 1.0
        else
          DO kk= 2, km
            if( pres_inpt(kk) > pstd(k) ) then
               frac=  ( pres_inpt(kk) - pstd(k) )/(pres_inpt(kk)-pres_inpt(kk-1))
               klev= kk
               exit
            endif
          ENDDO

          if( kk > km )  then
               klev= km;  frac= 0.0
          endif 
        endif

!             Complete pressure interpolation 
        DO i= 1, im
          DO j= 1, jm
            txy(i,j)= (1.0-frac)*teq_inpt(i,j,klev) + frac *teq_inpt(i,j,klev-1)
          ENDDO
        ENDDO

!             Carry out horizontal interpolation 
         call horiz_interp( txy(:,:), lonb_inpt, latb_inpt, lon, lat,    &
                                      tout(:,:,k), interp_method= 'bilinear' )
     ENDDO      ! -----------  end loop over k 

  else

       DO k= 1, nlevels
         call horiz_interp( teq_inpt(:,:,k), lonb_inpt, latb_inpt, lon, lat, &
                             tout(:,:,k), interp_method= 'bilinear' )
       ENDDO
  endif

   deallocate ( teq_inpt )
   deallocate ( txy )
   deallocate ( lat_inpt )
   deallocate ( latb_inpt )
   deallocate ( lonb_inpt )
   deallocate ( pres_inpt )


 end subroutine read_teq


#else
 subroutine Sim_phys(npx, npy, npz, is, ie, js, je, ng, nq,  &
                     u_dt, v_dt, t_dt, q_dt, u, v, ua, va,   &
                     pt, q, pe, delp, peln, oro, hydrostatic, &
                     pdt, agrid, ak, bk, rayf, p_ref, master, Time, time_total)

      integer, INTENT(IN) :: npx, npy, npz
      integer, INTENT(IN) :: is, ie, js, je, ng, nq
      real   , INTENT(IN) :: pdt
      real   , INTENT(IN) :: agrid(is-ng:ie+ng,js-ng:je+ng, 2)
      real   , INTENT(IN) :: ak(npz+1), bk(npz+1)
      logical, INTENT(IN) :: rayf, master
      real, INTENT(IN):: p_ref
      real, INTENT(IN):: oro(is:ie,js:je)       ! land fraction
      logical, INTENT(IN):: hydrostatic

      type(time_type), intent(in) :: Time
      real, INTENT(IN), optional:: time_total

      real   , INTENT(INOUT) :: u(is-ng:ie+  ng,js-ng:je+1+ng,npz)
      real   , INTENT(INOUT) :: v(is-ng:ie+1+ng,js-ng:je+  ng,npz)

      real, INTENT(INOUT)::   pt(is-ng:ie+ng,js-ng:je+ng,npz)
      real, INTENT(INOUT):: delp(is-ng:ie+ng,js-ng:je+ng,npz)
      real, INTENT(INOUT)::    q(is-ng:ie+ng,js-ng:je+ng,npz, nq)
      real, INTENT(INOUT)::   pe(is-1:ie+1 ,1:npz+1,js-1:je+1)
      real, INTENT(INOUT):: peln(is  :ie   ,1:npz+1,js  :je  )

! Tendencies:
      real, INTENT(INOUT):: u_dt(is-ng:ie+ng,js-ng:je+ng,npz)
      real, INTENT(INOUT):: v_dt(is-ng:ie+ng,js-ng:je+ng,npz)
      real, INTENT(INOUT):: t_dt(is:ie,js:je,npz)
      real, INTENT(INOUT):: q_dt(is:ie,js:je,npz,nq)

      real, INTENT(INOUT):: ua(is-ng:ie+ng,js-ng:je+ng,npz)
      real, INTENT(INOUT):: va(is-ng:ie+ng,js-ng:je+ng,npz)

! Local
      real:: u2(is:ie,npz), v2(is:ie,npz)
      real, allocatable, SAVE ::  rf(:)
      real p2(is:ie,npz), prec(is:ie,js:je)
      real sst(is:ie,js:je)
      real qs(is:ie)
      real pedge(npz+1)
      real pref(npz)
      real  sday, rkv, rka, rks, rkt, sigb, rsgb
      real  rk1
      real  sigl
      real  tmp, rkq
      real  pc, c1
      real cooling
      real frac
      real pi, rdt                    ! Reciprocal of dt
      logical used
      logical sh_con
      integer  i,j,k, iq

      sh_con = .true.

      pi = 4.*atan(1.)
      sday = 24.*3600.

      do k=1,npz+1
         pedge(k) = ak(k) + bk(k)*p_ref
      enddo
      do k=1,npz
         pref(k) = (pedge(k+1)-pedge(k)) / log(pedge(k+1)/pedge(k))
      enddo

        rkt = pdt / (4.*sday)       ! momentum
        rk1 = pdt / (4.*3600.)      ! temp
        rkq = pdt / (4.*3600.)      ! sphum
        rkv = pdt / sday

        rdt = 1. / pdt
        sigb = 0.7
        rsgb = 1./(1.-sigb)
        cooling = pdt / sday

        do k=1,npz
           if ( pref(k) >75.E2 ) then
           do j=js,je
              do i=is,ie
                 tmp = pt(i,j,k) - cooling
                 if ( tmp > 200. ) then
                      t_dt(i,j,k) = tmp
                 else
                      t_dt(i,j,k) = pt(i,j,k)
                 endif
              enddo 
           enddo 
           else
             tmp = 200.
             do j=js,je
                do i=is,ie
                   t_dt(i,j,k) = (pt(i,j,k)+rkt*tmp)/(1.+rkt)
                enddo 
             enddo 
           endif
        enddo


      do j=js,je
      do i=is,ie
#ifdef UNIFORM_SST
         sst(i,j) = 302.
#else
#ifdef FV_LAND
         sst(i,j) = 240. + 63.*cos(agrid(i,j,2))**3
#else
! Neale & Hoskins 2001
         if ( abs(agrid(i,j,2)) <  pi/3. ) then
              sst(i,j) = 273.16 + 27.*(1.-sin(1.5*agrid(i,j,2))**2)
         else
              sst(i,j) = 273.16
         endif
#endif
#endif
      enddo
      enddo

!--- Simple Moist Phys -------------------------------------------------------
     if( nq/=0 )  then

      do j=js,je
         do k=1,npz
            do i=is,ie
               p2(i,k) = delp(i,j,k) / (peln(i,k+1,j)-peln(i,k,j))
               q_dt(i,j,k,1) = q(i,j,k,1)
            enddo
         enddo

         do i=is,ie
#ifdef FV_LAND
            t_dt(i,j,npz) = (t_dt(i,j,npz)+rk1*(1.-oro(i,j))*sst(i,j))/(1.+rk1*(1.-oro(i,j)))
#else
            t_dt(i,j,npz) = (t_dt(i,j,npz)+rk1*sst(i,j))/(1.+rk1)
#endif
         enddo

         call qsmith(ie-is+1, 1, 1, t_dt(is:ie,j,npz), p2(is,npz), qs)

!----------------
! Surface fluxes:
!----------------
         do i=is,ie
            rkt = rkq * max( 0.5, sqrt(ua(i,j,npz)**2+va(i,j,npz)**2) )
!                     * min(5.E2, delp(i,j,npz) ) / delp(i,j,npz)
#ifdef FV_LAND
            q_dt(i,j,npz,1) = (q(i,j,npz,1) + rkt*(1.-oro(i,j))*qs(i))/(1.+rkt*(1.-oro(i,j)))
#else
            q_dt(i,j,npz,1) = (q(i,j,npz,1) + rkt*qs(i))/(1.+rkt)
#endif
#ifdef SURF_DRAG
            u2(i,npz) = u2(i,npz) / (1.+2.*rkt)
            v2(i,npz) = v2(i,npz) / (1.+2.*rkt)
#endif
         enddo

         do k=1,npz
            do i=is,ie
               u2(i,k) = ua(i,j,k)
               v2(i,k) = va(i,j,k)
            enddo
         enddo

         if ( nq > 1 ) then
           do iq=2,nq
             do k=1,npz
                do i=is,ie
                   q_dt(i,j,k,iq) = q(i,j,k,iq)
                enddo
             enddo
           enddo
         endif

         if ( sh_con )     &
         call conv_adj(ie-is+1,   npz, nq,  pdt,   9.8,   1004.64,  p2,   &
                       delp(is:ie,j,1:npz), pe(is:ie,1:npz+1,j), peln(is:ie,1:npz+1,j), &
                       t_dt(is:ie,j,1:npz), q_dt(is:ie,j,1:npz,1), u2, v2, pref, 90.)

         if ( hydrostatic ) then
           call cond_ls(ie-is+1,  npz,  pdt,  9.8,  1004.64,  p2,  prec(is,j),   &
                        delp(is:ie,j,1:npz), pe(is:ie,1:npz+1,j), peln(is:ie,1:npz+1,j), &
                        t_dt(is:ie,j,1:npz), q_dt(is:ie,j,1:npz,1), pref)
         else
! Constant volume heating of the atmosphere:
           call cond_ls_nh(ie-is+1,  npz,  pdt,  9.8,  1004.64,  p2,  prec(is,j),   &
                        delp(is:ie,j,1:npz), pe(is:ie,1:npz+1,j), peln(is:ie,1:npz+1,j), &
                        t_dt(is:ie,j,1:npz), q_dt(is:ie,j,1:npz,1), pref)
         endif

         do k=1,npz
            do i=is,ie
               t_dt(i,j,k  ) = (t_dt(i,j,k  ) - pt(i,j,k) ) * rdt
               q_dt(i,j,k,1) = (q_dt(i,j,k,1) - q(i,j,k,1)) * rdt
! Winds:
                u_dt(i,j,k) = (u2(i,k) - ua(i,j,k)) * rdt
                v_dt(i,j,k) = (v2(i,k) - va(i,j,k)) * rdt
#ifndef SURF_DRAG
!------------------------
! Surface drag as in H-S
!------------------------
                sigl = 0.5*(pe(i,k,j)+pe(i,k+1,j)) / pe(i,npz+1,j)
                frac = rkv * (sigl-sigb)*rsgb
                if (frac > 0.) then
                    u_dt(i,j,k) = u_dt(i,j,k) - u2(i,k)*frac/(1.+frac) * rdt
                    v_dt(i,j,k) = v_dt(i,j,k) - v2(i,k)*frac/(1.+frac) * rdt
                endif
#endif
            enddo
         enddo
      enddo   ! j-loop

      if ( id_prec > 0 ) used = send_data(id_prec, prec, Time)

! Last tracer is the age tracer
         if( nq>1 ) then
             do k=1,npz
                do j=js,je
                   do i=is,ie
                      q_dt(i,j,k,nq) = q(i,j,k,nq)
                   enddo
                enddo
             enddo
             call age_of_air(is, ie, js, je, npz, 0, time_total, pe, q_dt(is,js,1,nq))
             do k=1,npz
                do j=js,je
                   do i=is,ie
                      q_dt(i,j,k,nq) = (q_dt(i,j,k,nq) - q(i,j,k,nq) ) * rdt
                   enddo
                enddo
             enddo
         endif
      endif

  end subroutine Sim_phys

#endif MARS_GCM

  subroutine conv_adj(im, km, nq, dt, gravit, cpair, pm, delp,   &
                      p_edge, pe_ln, t, q, ua, va, pref, tau) 
!----------------------------------------
! Do non-precipitating shallow convection
!----------------------------------------
! Local Physical constants:
      real, parameter:: rdg  = 287.04
      real, parameter:: rvg  = 461.50
      real, parameter:: latv = 2.5104e06
!-------------------------------------------
      integer, intent(in):: im, km, nq
      real, intent(in):: dt, gravit, cpair, tau
      real, intent(in)::   pm(im,km)       ! Pressure at layer midpoints
      real, intent(in):: delp(im,km)       ! Delta p at each model level
      real, intent(in):: p_edge(im,km+1)
      real, intent(in):: pe_ln(im,km+1)
      real, intent(in):: pref(km)          ! Reference pressure (pa)
! Output arguments
      real, intent(inout)::  t(im,km)      ! Temperature
      real, intent(inout)::  q(im,km,nq)   ! Specific humidity & tracers
      real, intent(inout):: ua(im,km)
      real, intent(inout):: va(im,km)
!---------------------------Local variables-----------------------------
      real, dimension(im,km):: t0, u0, v0, gz, qsat, hd
      real  q0(im,km,nq)
      real gzh(im)
      real zvir, dh, dq, tv, tvm, qmix, h0, mc, fra, rdm, rl, rk, rz
      integer mcond, kcond
      integer i, k, n, m, iq

      zvir = rvg/rdg - 1.     ! = 0.607789855
        rz = rvg - rdg        ! rz = zvir * rdgas
        rl = 1. / latv
        rk = cpair/rdg + 1.

! Narrow the convective domain
       mcond = 2
       do k=2,km
          if(pref(k) > 30.E2) then
             mcond = k
             go to 120
          endif
       enddo
120   continue

    fra = dt/tau
    if ( fra < 1. ) then
       m = 3
       do k=mcond,km
          do i=1,im
             t0(i,k) = t(i,k)
             u0(i,k) = ua(i,k)
             v0(i,k) = va(i,k)
          enddo
       enddo
       do iq=1,nq
          do k=mcond,km
             do i=1,im
                q0(i,k,iq) = q(i,k,iq)
             enddo
          enddo
       enddo
    else
       m = 3 + fra
    endif

!------------
! Compute gz
!------------
    do i=1,im
       gzh(i) = 0.
    enddo
    do k=km,mcond,-1
       do i=1,im
              tvm = t(i,k)*(1.+zvir*q(i,k,1))
              tv  = rdg*tvm
          gz(i,k) = gzh(i) + tv*(1.-p_edge(i,k)/pm(i,k))
          hd(i,k) = cpair*tvm + gz(i,k) + 0.5*(ua(i,k)**2+va(i,k)**2)
           gzh(i) = gzh(i) + tv*(pe_ln(i,k+1)-pe_ln(i,k))
       enddo
    enddo

   do n=1,m

!     if ( n==1 ) then
!          kcond = km-2       ! "PBL"
!     else
           kcond = min(km-2, mcond + m - n)
!     endif

      call qsmith(im, km, kcond, t, pm, qsat)
      do k=km,kcond+1,-1
         do i=1,im
            dh = hd(i,k) - hd(i,k-1)
            if ( dh > 0. ) then
! Dry convection:
                rdm = 1. / (delp(i,k-1) + delp(i,k))
                h0 = (delp(i,k-1)*hd(i,k-1)+delp(i,k)*hd(i,k)) * rdm
                hd(i,k-1) = h0
                hd(i,k  ) = h0
                do iq=1,nq
                   qmix = (q(i,k-1,iq)*delp(i,k-1)+q(i,k,iq)*delp(i,k)) * rdm
                   q(i,k-1,iq) = qmix
                   q(i,k  ,iq) = qmix
                enddo
! Momentum mixing:
                qmix = (ua(i,k-1)*delp(i,k-1)+ua(i,k)*delp(i,k)) * rdm
                ua(i,k-1) = qmix
                ua(i,k  ) = qmix
                qmix = (va(i,k-1)*delp(i,k-1)+va(i,k)*delp(i,k)) * rdm
                va(i,k-1) = qmix
                va(i,k  ) = qmix
            else
                dq = q(i,k,1) - qsat(i,k-1) + dh*rl
                if ( dq > 0. ) then
! Moist Mixing/convection
                    qmix = (q(i,k-1,1)*delp(i,k-1)+q(i,k,1)*delp(i,k))/(delp(i,k-1)+delp(i,k))
                    dq = min( q(i,k,1)-qmix, dq )
                  if( dq > 1.E-7 ) then
! Compute equivalent mass flux: mc
                    mc = dq/(q(i,k,1)-q(i,k-1,1)) * delp(i,k)
                    do iq=1,nq
                       h0 = mc*(q(i,k,iq)-q(i,k-1,iq))
                       q(i,k-1,iq) = q(i,k-1,iq) + h0/delp(i,k-1)
                       q(i,k  ,iq) = q(i,k  ,iq) - h0/delp(i,k  )
                    enddo
! h:
                    h0 = mc*dh    ! dh < 0
                    hd(i,k-1) = hd(i,k-1) + h0/delp(i,k-1)
                    hd(i,k  ) = hd(i,k  ) - h0/delp(i,k  )
! u:
                    h0 = mc*(ua(i,k)-ua(i,k-1))
                    ua(i,k-1) = ua(i,k-1) + h0/delp(i,k-1)
                    ua(i,k  ) = ua(i,k  ) - h0/delp(i,k  )
! v:
                    h0 = mc*(va(i,k)-va(i,k-1))
                    va(i,k-1) = va(i,k-1) + h0/delp(i,k-1)
                    va(i,k  ) = va(i,k  ) - h0/delp(i,k  )
                  endif
               endif
            endif
         enddo
      enddo

!-------------
! Retrive Temp
!-------------
      do i=1,im
         gzh(i) = 0.
      enddo
      do k=km,kcond,-1
         do i=1,im
            t(i,k) = (hd(i,k)-gzh(i)-0.5*(ua(i,k)**2+va(i,k)**2))/(rk-p_edge(i,k)/pm(i,k))
            gzh(i) = gzh(i) + t(i,k)*(pe_ln(i,k+1)-pe_ln(i,k))
            t(i,k) = t(i,k) / ( rdg + rz*q(i,k,1) )
         enddo
      enddo

   enddo

   if ( fra < 1. ) then
      do k=mcond,km
         do i=1,im
             t(i,k) = t0(i,k) + ( t(i,k) - t0(i,k))*fra
            ua(i,k) = u0(i,k) + (ua(i,k) - u0(i,k))*fra
            va(i,k) = v0(i,k) + (va(i,k) - v0(i,k))*fra
         enddo
      enddo
      do iq=1,nq
         do k=mcond,km
            do i=1,im
               q(i,k,iq) = q0(i,k,iq) + (q(i,k,iq) - q0(i,k,iq))*fra
            enddo
         enddo
      enddo
   endif

  end subroutine conv_adj


  subroutine cond_ls_nh(plon, plev, dt, gravit, cpair, pmid, prec, pdel,   &
                     p_edge, pe_ln, t, q, pref) 
!
!     cpair = 1004.64
!     epsilo = rdgas/rvgas   ! 0.621971831
!     epsilo = 0.622
!     gravit = 9.80616
!     latvap = 2.5104e06
!     latice = 3.34e5
! Developer: S.-J. Lin

! Input arguments
      integer, intent(in):: plon, plev
      real, intent(in):: dt                   ! Physics time step
      real, intent(in):: gravit, cpair
      real, intent(in):: pmid(plon,plev)       ! Pressure at layer midpoints
      real, intent(in):: pdel(plon,plev)       ! Delta p at each model level
      real, intent(in):: p_edge(plon,plev+1)
      real, intent(in):: pe_ln(plon,plev+1)
      real, intent(in):: pref(plev)
 
! Output arguments
      real, intent(inout):: t(plon,plev)        ! Temperature
      real, intent(inout):: q(plon,plev)        ! Specific humidity
      real, intent(out)  :: prec(plon)          ! Large-scale precipitation rate
!
! Local:
!---------------------------Local variables-----------------------------
      real, parameter:: Tice = 273.16
      real evap(plon,plev)       ! Water evaporation rate
      real dqsdt(plon,plev)      ! Change of qsat with respect to temp.
      real denom                 ! Intermediate quantity
      real omeps                  ! 1 - 0.622
      real qsat(plon,plev)       ! Saturation specific humidity
      real rain(plon)            ! Rain (units of kg/m^2 s)
      real rga                    ! Reciprocal gravitatnl acceleration
      real zqcd(plon)            ! Intermed quantity (several actually)
      real rdt                    ! Reciprocal of dt
      real cndwtr(plon,plev)     ! Water condensation rate (kg/m**2/s)
      real evap_ef
      real relhum                 ! Relative humidity
      real dpovrg                 ! deltap/grav
      real ice(plon)             ! Solid precipation (ice , snow, etc)
      real pice(plon)            ! Temp storage for solid precipatation
      real dice
      real h0, h1
      real fac_frez
      real fac_melt
      real latvap, cld, sl
      real rvgas, zvir, cvair
      integer lcond               ! starting layer for cond. computation
      integer kcond               ! starting layer for evap
      integer i, k, n

! Physical constants:
!     rdgas = 287.04
      rvgas = 461.50

      cvair = cpair - rdgas
       zvir = rvgas/rdgas - 1.     ! = 0.607789855

      latvap = 2.5104e06
!     cld  = latvap/cpair      ! hydrostatic: constant-press heating
      cld  = latvap/cvair      ! constant volume heating of the air

      sl = 3.34E5

! Derived quantities:
      rga  = 1./gravit
      rdt  = 1./dt
      omeps = 1. - rdgas/rvgas

! Tunning parameters:
      evap_ef  = 1.0e-5 ! evaporation efficiency
      fac_frez = 0.5
      fac_melt = 0.5

! Compute starting L-S condensation level
       lcond = 2
       do k=2, plev
          if(pref(k) > 10.E2) then  
             lcond = k
             go to 111
          endif
       enddo
111   continue

!-------------------------------------------------------------------------
   call qsmith(plon, plev, lcond, t, pmid, qsat, dqsdt)

! First diagnose condensation rate due to stable processes
! Update column T and Q (evaporation process is `time-split')
 
      do k=lcond,plev
        do i=1,plon
           cndwtr(i,k) = 0.0
        end do
      end do
 
      do i=1,plon
         ice(i) = 0.
      enddo

      do i=1,plon
         pice(i) = 0.
      enddo

      do k=lcond,plev
! Calculate condensation-rate and new t- and q-values
         do i=1,plon
            zqcd(i) = max(0., qsat(i,k)*(q(i,k)/qsat(i,k)-1.)/(1.+cld*dqsdt(i,k)))
            if (q(i,k) < 0.0) zqcd(i) = 0.
            q(i,k)  = q(i,k) - zqcd(i)
            t(i,k)  = t(i,k) + zqcd(i)*cld
            h1  = pdel(i,k)*rga*rdt
            cndwtr(i,k) = cndwtr(i,k) + zqcd(i)*h1
            h0  = sl / (h1*cvair)
          if(t(i,k) < tice .and. cndwtr(i,k) > 0.) then
! ********************************
! ***   Freezing *****
! ********************************
! The release of latent heat due to freezing should not warm
! the atmosphere above freezing point.
             dice = min(cndwtr(i,k), fac_frez*(tice-t(i,k))/h0)
             pice(i) = pice(i) + dice
             t(i,k) = t(i,k) + dice*h0
             cndwtr(i,k) = cndwtr(i,k) - dice
        elseif(t(i,k) > tice .and. pice(i) > 0.) then

! ********************************
! ***   Melting *****
! ********************************
! The heat needed for melting should not cool the atmosphere
! below frezzing
             dice = min(pice(i), fac_melt*(t(i,k)-tice)/h0)
             cndwtr(i,k) = cndwtr(i,k) + dice
             t(i,k) = t(i,k) - dice*h0
             pice(i) = pice(i) - dice
          endif
          end do
      end do
 
      do i=1,plon
         ice(i) = ice(i) + pice(i)
      enddo

       kcond = lcond
       do k=lcond, plev
          if(pref(k) > 35000.) then   ! rain re-evap starts from 350 mb
             kcond = k
             go to 222
          endif
       enddo
222   continue

!---------------------------------------------------------------

! Initialize rain vector (will be updated as rain falls through column)
      do i=1,plon
         rain(i) = max(cndwtr(i,kcond), 0.0)
      end do

      call qsmith(plon, plev, lcond, t, pmid, qsat)
 
! Evaporate condensate (Sundqvist, 1988: Physically
! Based Modelling ..., pp 433-461, Schlesinger, Ed., Kluwer Academic)
! variable evap has units of 1/s; variable rain has units of kg/m**2/s
! rain is used to accumuluate unevaporated rain water on the way down

      do k=kcond+1,plev
         do i=1,plon
            dpovrg  = pdel(i,k)*rga
            relhum  = 1. - q(i,k)/qsat(i,k)
            evap(i,k) = max(evap_ef*relhum*sqrt(rain(i)), 0.0)
! Nudge factor to prevent over-evaporation: 0.5
            evap(i,k) = min(evap(i,k), 0.5*(qsat(i,k)-q(i,k))*rdt)
            evap(i,k) = min(rain(i)/dpovrg, evap(i,k))
            q(i,k) = q(i,k) + evap(i,k)*dt
            t(i,k) = t(i,k) - evap(i,k)*dt*cld
            rain(i) = max(rain(i) - evap(i,k)*dpovrg + cndwtr(i,k),0.0)
         end do
      end do

      do i=1,plon
         prec(i) = (rain(i)+ice(i)) * 86400.
      end do
 
  end subroutine cond_ls_nh

  subroutine cond_ls(plon, plev, dt, gravit, cpair, pmid, prec, pdel,   &
                     p_edge, pe_ln, t, q, pref) 
!
!     cpair = 1004.64
!     epsilo = rdgas/rvgas   ! 0.621971831
!     epsilo = 0.622
!     gravit = 9.80616
!     latvap = 2.5104e06
!     latice = 3.34e5
! Developer: S.-J. Lin

! Input arguments
      integer, intent(in):: plon, plev
      real, intent(in):: dt                   ! Physics time step
      real, intent(in):: gravit, cpair
      real, intent(in):: pmid(plon,plev)       ! Pressure at layer midpoints
      real, intent(in):: pdel(plon,plev)       ! Delta p at each model level
      real, intent(in):: p_edge(plon,plev+1)
      real, intent(in):: pe_ln(plon,plev+1)
      real, intent(in):: pref(plev)
 
! Output arguments
      real, intent(inout):: t(plon,plev)        ! Temperature
      real, intent(inout):: q(plon,plev)        ! Specific humidity
      real, intent(out)  :: prec(plon)          ! Large-scale precipitation rate
!
! Local:
!---------------------------Local variables-----------------------------
      real, parameter:: Tice = 273.16
      real evap(plon,plev)       ! Water evaporation rate
      real dqsdt(plon,plev)      ! Change of qsat with respect to temp.
      real denom                 ! Intermediate quantity
      real omeps                  ! 1 - 0.622
      real qsat(plon,plev)       ! Saturation specific humidity
      real rain(plon)            ! Rain (units of kg/m^2 s)
      real rga                    ! Reciprocal gravitatnl acceleration
      real zqcd(plon)            ! Intermed quantity (several actually)
      real rdt                    ! Reciprocal of dt
      real cndwtr(plon,plev)     ! Water condensation rate (kg/m**2/s)
      real evap_ef
      real relhum                 ! Relative humidity
      real dpovrg                 ! deltap/grav
      real ice(plon)             ! Solid precipation (ice , snow, etc)
      real pice(plon)            ! Temp storage for solid precipatation
      real dice
      real h0, h1
      real fac_frez
      real fac_melt
      real latvap, cldcp, sl
      real rvgas, zvir
      integer lcond               ! starting layer for cond. computation
      integer kcond               ! starting layer for evap
      integer i, k, n

! Physical constants:
!     rdgas = 287.04
      rvgas = 461.50
      zvir = rvgas/rdgas - 1.     ! = 0.607789855

      latvap = 2.5104e06
      cldcp  = latvap/cpair
      sl = 3.34E5

! Derived quantities:
      rga  = 1./gravit
      rdt  = 1./dt
      omeps = 1. - rdgas/rvgas

! Tunning parameters:
      evap_ef  = 2.5e-5 ! evaporation efficiency
      fac_frez = 0.75
      fac_melt = 0.75

! Compute starting L-S condensation level
       lcond = 2
       do k=2, plev
          if(pref(k) > 10.E2) then  
             lcond = k
             go to 111
          endif
       enddo
111   continue

!-------------------------------------------------------------------------
   call qsmith(plon, plev, lcond, t, pmid, qsat, dqsdt)

! First diagnose condensation rate due to stable processes
! Update column T and Q (evaporation process is `time-split')
 
      do k=lcond,plev
        do i=1,plon
           cndwtr(i,k) = 0.0
        end do
      end do
 
      do i=1,plon
         ice(i) = 0.
      enddo

      do i=1,plon
         pice(i) = 0.
      enddo

      do k=lcond,plev
! Calculate condensation-rate and new t- and q-values
          do i=1,plon
!           zqcd(i) = max(0., qsat(i,k)*(q(i,k)/qsat(i,k)-1.)/(1.+cldcp*dqsdt(i,k)))
            zqcd(i) = max(0., qsat(i,k)*(q(i,k)/qsat(i,k)-1.)   &
                    /    (1.+latvap/(cpair-rdgas*p_edge(i,1)/pmid(i,k))*dqsdt(i,k)))
            if (q(i,k) < 0.0) zqcd(i) = 0.
            q(i,k)  = q(i,k) - zqcd(i)
!           t(i,k)  = t(i,k) + zqcd(i)*cldcp
            t(i,k)  = t(i,k) + zqcd(i)*latvap/(cpair-rdgas*p_edge(i,1)/pmid(i,k))
            h1  = pdel(i,k)*rga*rdt
            cndwtr(i,k) = cndwtr(i,k) + zqcd(i)*h1
            h0  = sl/( h1*(cpair-rdgas*p_edge(i,1)/pmid(i,k)) )
          if(t(i,k) < tice .and. cndwtr(i,k) > 0.) then
! ********************************
! ***   Freezing *****
! ********************************
! The release of latent heat due to freezing should not warm
! the atmosphere above freezing point.
             dice = min(cndwtr(i,k), fac_frez*(tice-t(i,k))/h0)
             pice(i) = pice(i) + dice
             t(i,k) = t(i,k) + dice*h0
             cndwtr(i,k) = cndwtr(i,k) - dice
        elseif(t(i,k) > tice .and. pice(i) > 0.) then

! ********************************
! ***   Melting *****
! ********************************
! The heat needed for melting should not cool the atmosphere
! below frezzing
             dice = min(pice(i), fac_melt*(t(i,k)-tice)/h0)
             cndwtr(i,k) = cndwtr(i,k) + dice
             t(i,k) = t(i,k) - dice*h0
             pice(i) = pice(i) - dice
          endif
          end do
      end do
 
      do i=1,plon
         ice(i) = ice(i) + pice(i)
      enddo

       kcond = lcond
       do k=lcond, plev
          if(pref(k) > 35000.) then   ! rain re-evap starts from 350 mb
             kcond = k
             go to 222
          endif
       enddo
222   continue

!---------------------------------------------------------------

! Initialize rain vector (will be updated as rain falls through column)
      do i=1,plon
         rain(i) = max(cndwtr(i,kcond),0.0)
      end do

      call qsmith(plon, plev, lcond, t, pmid, qsat)
 
! Evaporate condensate (Sundqvist, 1988: Physically
! Based Modelling ..., pp 433-461, Schlesinger, Ed., Kluwer Academic)
! variable evap has units of 1/s; variable rain has units of kg/m**2/s
! rain is used to accumuluate unevaporated rain water on the way down
 

      do k=kcond+1,plev
         do i=1,plon
            dpovrg  = pdel(i,k)*rga
            relhum  = 1. - q(i,k)/qsat(i,k)
            evap(i,k) = max(evap_ef*relhum*sqrt(rain(i)), 0.0)
! Nudge factor to prevent over-evaporation: 0.9
            evap(i,k) = min(evap(i,k), 0.9*(qsat(i,k)-q(i,k))*rdt)
            evap(i,k) = min(rain(i)/dpovrg, evap(i,k))
            q(i,k) = q(i,k) + evap(i,k)*dt
!           t(i,k) = t(i,k) - evap(i,k)*dt*cldcp
            t(i,k) = t(i,k) - evap(i,k)*dt*latvap/(cpair-rdgas*p_edge(i,1)/pmid(i,k))
            rain(i) = max(rain(i) - evap(i,k)*dpovrg + cndwtr(i,k),0.0)
         end do
      end do

      do i=1,plon
         prec(i) = (rain(i)+ice(i)) * 86400.
      end do
 
      end subroutine cond_ls

      subroutine age_of_air(is, ie, js, je, km, ng, time, pe, q)

      integer is, ie, js, je
      integer km
      integer ng

! q is the age tracer
! Need to be converted to mixing ratio (mass of tracer / dry_air-mass)
! Ignore this inconsistency for now.

      real, intent(inout):: pe(is-1:ie+1, km+1, js-1:je+1)
      real, intent(in):: time        ! accumulated time since init
      real, intent(inout):: q(is-ng:ie+ng,js-ng:je+ng,km)

! Local
      integer i, j, k
      real p_source      ! source level (pa)
      real ascale
      real tiny
      parameter ( tiny = 1.e-6 )
      parameter ( p_source = 75000. )
      parameter ( ascale = 5.e-6 / 60. )

!$omp parallel do private(i, j, k)
      do k=1,km
        do j=js,je
            do i=is,ie
               if( time < tiny ) then
                   q(i,j,k) = 0.
               elseif( pe(i,k,j) >= p_source ) then
                   q(i,j,k) = ascale * time
               endif
            enddo
        enddo          ! j-loop
      enddo             ! k-loop

      end subroutine age_of_air

 
end module hswf_mod
