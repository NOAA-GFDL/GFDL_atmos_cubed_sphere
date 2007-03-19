module hswf_mod

      use constants_mod, only: grav, rdgas
      use mpp_domains_mod, only: mpp_update_domains
      use mp_mod,          only: domain

      use fv_diagnostics_mod, only: id_prec
      use diag_manager_mod,   only: send_data
      use time_manager_mod,   only: time_type

implicit none
!-----------------------------------------------------------------------
      real, parameter:: Tice=273.16
      real  table(2621),des(2621)
      logical :: rf_initialized = .false.
      logical :: qs_initialized = .false.

private

public :: Held_Suarez_Strat, Held_Suarez_Tend, Sim_phys, age_of_air, cond_ls


contains

!-----------------------------------------------------------------------
 subroutine Held_Suarez_Tend(npx, npy, npz, is, ie, js, je, ng, nq,   &
                              u, v, pt, q, pe, delp, peln, pkz, pdt,  &
                              ua, va, u_dt, v_dt, t_dt, q_dt, grid,   &
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


      real   , INTENT(IN   ) :: grid(is-ng:ie+ng,js-ng:je+ng, 1:2)
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
                 tey = ap0k*( 315.0 - ty*SIN(grid(i,j,2))*SIN(grid(i,j,2)) )
                 tez = tz*( ap0k/akap )*COS(grid(i,j,2))*COS(grid(i,j,2)) 

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
                   dt_tropic = 2.25*COS(grid(i,j,2)) * dz
                   teq(i,j,k)  =  teq(i,j,k+1) + dt_tropic
                   t_dt(i,j,k) = relx*(teq(i,j,k)-pt(i,j,k))/(1.+relx) * rdt
                 elseif (strat .and. pl(i,j,k) <= 100.)  then
!
! Mesosphere: defined as the region above 1 mb
!
                   dz = h0 * log(pl(i,j,k+1)/pl(i,j,k))
                   dt_tropic = -2.25*COS(grid(i,j,2)) * dz
                   teq(i,j,k) = teq(i,j,k+1) + dt_tropic
                   t_dt(i,j,k) = ((pt(i,j,k)+rms*teq(i,j,k))*rmr - pt(i,j,k))*rdt
                 else

! Trop:  strictly Held-Suarez

                   sigl = pl(i,j,k)/pe(i,npz+1,j)
                   f1 = max(0., (sigl-sigb) * rsgb )
                   teq(i,j,k) = tey - tez*(log(pkz(i,j,k))+algpk)
                   teq(i,j,k) = max(t0, teq(i,j,k)*pkz(i,j,k))
                   rkt = rka + (rks-rka)*f1*(COS(grid(i,j,2))**4.0)
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
                              ua, va, grid, ak, bk, ks, strat,  &
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

      real   , INTENT(IN   ) :: grid(is-ng:ie+ng,js-ng:je+ng, 1:2)
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
                 tey = ap0k*( 315.0 - ty*SIN(grid(i,j,2))*SIN(grid(i,j,2)) )
                 tez = tz*( ap0k/akap )*COS(grid(i,j,2))*COS(grid(i,j,2)) 

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
                   dt_tropic = 2.25*COS(grid(i,j,2)) * dz
                   teq(i,j,k)  =  teq(i,j,k+1) + dt_tropic
                   pt(i,j,k) = (pt(i,j,k)+relx*teq(i,j,k))/(1.+relx)
                 elseif (strat .and. pl(i,j,k) <= 100.)  then
!
! Mesosphere: defined as the region above 1 mb
!
                   dz = h0 * log(pl(i,j,k+1)/pl(i,j,k))
                   dt_tropic = -2.25*COS(grid(i,j,2)) * dz
                   teq(i,j,k) = teq(i,j,k+1) + dt_tropic
                   pt(i,j,k) = (pt(i,j,k)+rms*teq(i,j,k))*rmr
                 else

! Trop:  strictly Held-Suarez

                   sigl = pl(i,j,k)/pe(i,npz+1,j)
                   f1 = max(0., (sigl-sigb) * rsgb )
                   teq(i,j,k) = tey - tez*(log(pkz(i,j,k))+algpk)
                   teq(i,j,k) = max(t0, teq(i,j,k)*pkz(i,j,k))
                   rkt = rka + (rks-rka)*f1*(COS(grid(i,j,2))**4.0)
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

 subroutine Sim_phys(npx, npy, npz, is, ie, js, je, ng, nq,  &
                     u_dt, v_dt, t_dt, q_dt, u, v, ua, va,   &
                     pt, q, pe, delp, peln, oro, &
                     pdt, grid, ak, bk, rayf, master, Time, time_total)

      integer, INTENT(IN) :: npx, npy, npz
      integer, INTENT(IN) :: is, ie, js, je, ng, nq
      real   , INTENT(IN) :: pdt
      real   , INTENT(IN) :: grid(is-ng:ie+ng,js-ng:je+ng, 1:2)
      real   , INTENT(IN) :: ak(npz+1), bk(npz+1)
      logical, INTENT(IN) :: rayf, master
      real, INTENT(IN):: oro(is:ie,js:je)       ! land fraction

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
      real p2(is:ie,npz), prec(is:ie,js:je)
      real sst(is:ie,js:je)
      real qs(is:ie)
      real pedge(npz+1)
      real pref(npz)

      real  sday, rkv, rka, rks, rkt, sigb, rsgb
      real  rk1
      real  tmp, rkq
      real  sigl
      real  pc, c1
      real cooling
      real pi, rdt                    ! Reciprocal of tdt

      real, allocatable, SAVE ::  rf(:)
      real frac
      integer  i,j,k

      logical used

      pi = 4.*atan(1.)
      sday = 24.*3600.

      do k=1,npz+1
         pedge(k) = ak(k) + bk(k)*1.E5
      enddo
      do k=1,npz
         pref(k) = (pedge(k+1)-pedge(k)) / log(pedge(k+1)/pedge(k))
      enddo

      if ( rayf .and. (.not. rf_initialized) ) then

          allocate( rf(npz) )

          c1 = 1. / (10.*sday)
          pc = 1.E2
          do k=1,npz
             rf(k) = c1*(1.+tanh(log10(pc/pref(k))))
             if(master) write(6,*) k, 0.01*pref(k), 1./(rf(k)*sday)
             rf(k) = 1./(1.+pdt*rf(k))
          enddo

          if(master) write(6,*) 'Rayleigh friction initialized.'
          rf_initialized = .true.

      endif

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
!                  t_dt(i,j,k) = pt(i,j,k)
                enddo 
             enddo 
           endif
        enddo


      do j=js,je
      do i=is,ie
#ifdef UNIFORM_SST
         sst(i,j) = 300.
#else
#ifdef FV_LAND
         sst(i,j) = 240. + 63.*cos(grid(i,j,2))**3
#else
! Neale & Hoskins 2001
         if ( abs(grid(i,j,2)) <  pi/3. ) then
              sst(i,j) = 273.16 + 27.*(1.-sin(1.5*grid(i,j,2))**2)
         else
              sst(i,j) = 273.16
         endif
#endif
#endif
      enddo
      enddo

      do 2000 k=1,npz
        if (rayf .and. pref(k) < 50.E2) then
! Stratosphere: Apply Rayleigh friction
#ifdef RAYF_DGRID
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
#else
          do j=js,je
             do i=is,ie
                u_dt(i,j,k) = ua(i,j,k)*(rf(k) - 1.) * rdt
                v_dt(i,j,k) = va(i,j,k)*(rf(k) - 1.) * rdt
             enddo
          enddo
#endif
        else
! Surface Rayleigh friction according to Held-Suarez
          do j=js,je
             do i=is,ie
                sigl = 0.5*(pe(i,k,j)+pe(i,k+1,j)) / pe(i,npz+1,j)
                frac = rkv * (sigl-sigb)*rsgb
                if (frac > 0.) then
                    u_dt(i,j,k) = -ua(i,j,k)*frac/(1.+frac) * rdt
                    v_dt(i,j,k) = -va(i,j,k)*frac/(1.+frac) * rdt
!                   rkt = rkv * ((sigl-sigb)*rsgb)**2
!                   t_dt(i,j,k) = (t_dt(i,j,k)+rkt*sst(i,j))/(1.+rkt)
                endif
             enddo
          enddo
        endif
2000  continue

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

         do i=is,ie
            rkt = rkq * max(0.5, sqrt(ua(i,j,npz)**2+va(i,j,npz)**2))
!                     * min(5.E2, delp(i,j,npz) ) / delp(i,j,npz)
#ifdef FV_LAND
            q_dt(i,j,npz,1) = (q(i,j,npz,1) + rkt*(1.-oro(i,j))*qs(i))/(1.+rkt*(1.-oro(i,j)))
#else
            q_dt(i,j,npz,1) = (q(i,j,npz,1) + rkt*qs(i))/(1.+rkt)
#endif
         enddo

         call cond_ls(ie-is+1,  npz,  pdt,  9.8,  1004.64,  p2,  prec(is,j),   &
                      delp(is:ie,j,1:npz), pe(is:ie,1:npz+1,j), peln(is:ie,1:npz+1,j), &
                      t_dt(is:ie,j,1:npz), q_dt(is:ie,j,1:npz,1), pref, .true.) 

            do k=1,npz
               do i=is,ie
                  t_dt(i,j,k  ) = (t_dt(i,j,k  ) - pt(i,j,k) ) * rdt
                  q_dt(i,j,k,1) = (q_dt(i,j,k,1) - q(i,j,k,1)) * rdt
               enddo
            enddo

         enddo

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
!-----------------------------------------------------------------------------


  end subroutine Sim_phys

  subroutine cond_ls(plon, plev, tdt, gravit, cpair, pmid,   &
                     prec, pdel, p_edge, pe_ln, t, q, pref, shallow_con ) 
!
!     cpair = 1004.64
!     epsilo = 0.622
!     gravit = 9.80616
!     latvap = 2.5104e06
!     latice = 3.34e5
! Developer: S.-J. Lin

! Input arguments
      integer, intent(in):: plon, plev
      real, intent(in):: tdt                   ! Physics time step
      real, intent(in):: gravit, cpair
      real, intent(in):: pmid(plon,plev)       ! Pressure at layer midpoints
      real, intent(in):: pdel(plon,plev)       ! Delta p at each model level
      real, intent(in):: p_edge(plon,plev+1)
      real, intent(in):: pe_ln(plon,plev+1)
      real, intent(in):: pref(plev)
      logical, intent(in):: shallow_con    ! do shallow convection
 
! Output arguments
      real, intent(inout):: t(plon,plev)        ! Temperature
      real, intent(inout):: q(plon,plev)        ! Specific humidity
      real, intent(out)  :: prec(plon)          ! Large-scale precipitation rate
!
! Local:
!---------------------------Local variables-----------------------------
      real evap(plon,plev)       ! Water evaporation rate
      real dqsdt(plon,plev)      ! Change of qsat with respect to temp.
      real denom                 ! Intermediate quantity
      real omeps                  ! 1 - 0.622
      real qsat(plon,plev)       ! Saturation specific humidity
      real rain(plon)            ! Rain (units of kg/m^2 s)
      real rga                    ! Reciprocal gravitatnl acceleration
      real zqcd(plon)            ! Intermed quantity (several actually)
      real rdt                    ! Reciprocal of tdt
      real cndwtr(plon,plev)     ! Water condensation rate (kg/m**2/s)
      real evap_ef
      real relhum                 ! Relative humidity
      real dpovrg                 ! deltap/grav

      real ice(plon)             ! Solid precipation (ice , snow, etc)
      real pice(plon)            ! Temp storage for solid precipatation
      real slcp                   ! Ice-phase latent heat / CP
      real dice
      real h0
      real h1
      real fac_frez
      real fac_melt
      real epsilo
      real latvap, clrv, cldcp
      real rvgas
      real zvir, hq, hs, dh, dq, tvm, qavg
      real gz(plon,plev), gzh(plon)
      integer lcond               ! starting layer for cond. computation
      integer kcond               ! starting layer for evap
      integer mcond               ! starting layer for convection
      integer i, k

! Physical constants:
!     rdgas = 287.04
      rvgas = 461.50
      zvir = rvgas/rdgas - 1.     ! = 0.607789855

      latvap = 2.5104e06
      cldcp  = latvap/cpair
      clrv   = latvap/rvgas
      epsilo = rdgas/rvgas   ! 0.621971831
      slcp = 3.34E5/cpair

! Derived quantities:
      rga  = 1./gravit
      rdt  = 1./tdt
      omeps = 1. - epsilo

! Tunning parameters:
      fac_frez = 0.85
      fac_melt = 0.85

! Compute starting L-S condensation level
       lcond = 2
       do k=2, plev
          if(pref(k) > 10.E2) then  
             lcond = k
             go to 111
          endif
       enddo
111   continue

      call qsmith(plon, plev, lcond, t, pmid, qsat, dqsdt)

!----------------------------------------
! Do non-precipitating shallow convection
!----------------------------------------
    if ( shallow_con ) then

! Narrow the convective domain
       mcond = lcond
       do k=lcond, plev
          if(pref(k) > 100.E2) then
             mcond = k
             go to 120
          endif
       enddo
120   continue

!------------
! Compute gz
!------------
      do i=1,plon
         gzh(i) = 0.
      enddo
      do k=plev,mcond,-1
         do i=1,plon
                 tvm = rdgas*t(i,k)*(1.+zvir*q(i,k))
             gz(i,k) = gzh(i) + tvm*(1.-p_edge(i,k)/pmid(i,k))
              gzh(i) = gzh(i) + tvm*(pe_ln(i,k+1)-pe_ln(i,k))
         enddo
      enddo

      do k=plev,mcond+1,-1
         do i=1,plon
            if ( pmid(i,k)>500.E2 ) then
               hq = cpair*t(i,k  ) + gz(i,k  ) + latvap*q(i,k)
               hs = cpair*t(i,k-1) + gz(i,k-1) + latvap*qsat(i,k-1)
               dh = hq - hs
               if ( dh>0. ) then
! Transport only the water vapor:
                 dq = min( max(0.,q(i,k)), dh/latvap)
                 q(i,k  ) = q(i,k  ) - dq
                 q(i,k-1) = q(i,k-1) + dq*pdel(i,k)/pdel(i,k-1)
               endif
            endif
         enddo
      enddo
    endif

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
            zqcd(i) = max(0., qsat(i,k)*(q(i,k)/qsat(i,k)-1.)/(1.+cldcp*dqsdt(i,k)))
            if (q(i,k) < 0.0) zqcd(i) = 0.
            q(i,k)  = q(i,k) - zqcd(i)
            t(i,k)  = t(i,k) + zqcd(i)*cldcp
            h1  = pdel(i,k)*rga*rdt
            cndwtr(i,k) = cndwtr(i,k) + zqcd(i)*h1
          if(t(i,k) < tice .and. cndwtr(i,k) > 0.) then
! ********************************
! ***   Freezing *****
! ********************************
             h0   = slcp/h1
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
             h0   = slcp/h1
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
 
      evap_ef = 1.5e-5 ! evaporation efficiency

      do k=kcond+1,plev
         do i=1,plon
            dpovrg  = pdel(i,k)*rga
            relhum  = 1. - q(i,k)/qsat(i,k)
            evap(i,k) = max(evap_ef*relhum*sqrt(rain(i)), 0.0)
! Nudge factor to prevent over-evaporation: 0.95
            evap(i,k) = min(evap(i,k), 0.95*(qsat(i,k)-q(i,k))*rdt)
            evap(i,k) = min(rain(i)/dpovrg, evap(i,k))
            q(i,k) = q(i,k) + evap(i,k)*tdt
            t(i,k) = t(i,k) - evap(i,k)*tdt*cldcp
            rain(i) = max(rain(i) - evap(i,k)*dpovrg + cndwtr(i,k),0.0)
         end do
      end do

      do i=1,plon
         prec(i) = (rain(i)+ice(i)) * 86400.
      end do
 
      end subroutine cond_ls


      subroutine qsmith(im, km, k1, t, p, q, dqdt)
! input T in deg K; p (Pa)
      integer, intent(in):: im, km, k1
      real, intent(in),dimension(im,km):: t, p
      real, intent(out),dimension(im,km):: q
      real table(2621),des(2621)
      real, intent(out), optional:: dqdt(im,km)
! Local:
      real es(im,km)
      real ap1
      real Tmin, Tmax, oms
      real, parameter:: esl = 0.621971831
      real:: dt=0.1
      integer i, k, it

      Tmin = Tice-160.
      Tmax = Tice+102.

      oms = 1. - esl
 
      if( .not. qs_initialized ) then
! Generate es table (dT = 0.1 deg. C)
           call esgfdl(table)
           do i=1,2620
              des(i) = table(i+1) - table(i)
           enddo
           des(2621) = des(2620)
           qs_initialized = .true.
      endif
 
      do k=k1,km
         do i=1,im
            ap1 = 10.*DIM(t(i,k), Tmin) + 1
            ap1 = min(2621., ap1)
            it = ap1
            es(i,k) = table(it) + (ap1-it)*des(it)
            q(i,k) = esl*es(i,k) / max(1.E-5, p(i,k)-oms*es(i,k))
            q(i,k) = min(1., q(i,k)) 
         enddo
      enddo

      if ( present(dqdt) ) then
      do k=k1,km
           do i=1,im
              ap1 = 10.*DIM(t(i,k),Tmin) + 1
              ap1 = min(2621., ap1) - 0.5
              it  = ap1
              dqdt(i,k) = 10.*esl*(des(it) + (ap1-it)*(des(it+1)-des(it)))/p(i,k)
           enddo
      enddo
      endif
 
      end subroutine qsmith
 

      subroutine esgfdl(table)
      real table (*)
      real esupc(200)
      real:: dt=0.1
      real esbasw, tbasw, esbasi, tbasi, Tmin, tem, aa, b, c, d, e, esh20 
      real wice, wh2o
      integer i

! Constants
      esbasw = 1013246.0
       tbasw =     373.16
      esbasi =    6107.1
       tbasi =     273.16
! ****************************************************
!  Compute es over ice between -160c and 0 c.
      Tmin = tbasi - 160.
!  see smithsonian meteorological tables page 350.
      do i=1,1600
         tem = Tmin+dt*real(i-1)
         aa  = -9.09718 *(tbasi/tem-1.0)
         b   = -3.56654 *alog10(tbasi/tem)
         c   =  0.876793*(1.0-tem/tbasi)
         e   = alog10(esbasi)
         table(i)=10**(aa+b+c+e)
      enddo
! *****************************************************
!  Compute es over water between -20c and 102c.
!  see smithsonian meteorological tables page 350.
      do  i=1,1221
          tem = 253.16+dt*real(i-1)
          aa  = -7.90298*(tbasw/tem-1)
          b   =  5.02808*alog10(tbasw/tem)
          c   = -1.3816e-07*(10**((1-tem/tbasw)*11.344)-1)
          d   =  8.1328e-03*(10**((tbasw/tem-1)*(-3.49149))-1)
          e   = alog10(esbasw)
          esh20  = 10**(aa+b+c+d+e)
          if (i <= 200) then
              esupc(i) = esh20
          else
              table(i+1400) = esh20
          endif
      enddo
!********************************************************************
!  Derive blended es over ice and supercooled water between -20c and 0c
      do i=1,200
         tem  = 253.16+dt*real(i-1)
         wice = 0.05*(273.16-tem)
         wh2o = 0.05*(tem-253.16)
         table(i+1400) = wice*table(i+1400)+wh2o*esupc(i)
      enddo

      do i=1,2621
         table(i) = table(i)*0.1
      enddo

      end subroutine esgfdl


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
