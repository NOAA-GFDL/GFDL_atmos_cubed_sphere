module fv_sg_mod

!-----------------------------------------------------------------------
! FV sub-grid mixing
!-----------------------------------------------------------------------
  use constants_mod,      only: rdgas, rvgas, cp_air, cp_vapor, hlv, hlf, kappa, grav
  use tracer_manager_mod, only: get_tracer_index
  use field_manager_mod,  only: MODEL_ATMOS

implicit none
private

public  fv_dry_conv, qsmith, neg_adj3

  real, parameter:: esl = 0.621971831
  real, parameter:: tice = 273.16
  real, parameter:: zvir =  rvgas/rdgas - 1.     ! = 0.607789855
  real, allocatable:: table(:),des(:)

!---- version number -----
  character(len=128) :: version = '$Id: fv_sg.F90,v 20.0 2013/12/13 23:04:34 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal_201409 $'

contains

 subroutine fv_dry_conv( isd, ied, jsd, jed, is, ie, js, je, km, nq, dt,    &
                         tau, delp, pe, peln, pkz, ta, qa, ua, va,  &
                         hydrostatic, w, delz, u_dt, v_dt, t_dt, q_dt )
! Dry convective adjustment-mixing
!-------------------------------------------
      integer, intent(in):: is, ie, js, je, km, nq
      integer, intent(in):: isd, ied, jsd, jed
      integer, intent(in):: tau         ! Relaxation time scale
      real, intent(in):: dt             ! model time step
      real, intent(in)::   pe(is-1:ie+1,km+1,js-1:je+1) 
      real, intent(in):: peln(is  :ie,  km+1,js  :je)
      real, intent(in):: delp(isd:ied,jsd:jed,km)      ! Delta p at each model level
      real, intent(in)::  pkz(is:ie,js:je,km)      ! Delta p at each model level
      real, intent(in):: delz(isd:ied,jsd:jed,km)      ! Delta p at each model level
      logical, intent(in)::  hydrostatic
! 
      real, intent(inout):: ua(isd:ied,jsd:jed,km)
      real, intent(inout):: va(isd:ied,jsd:jed,km)
      real, intent(inout)::  w(isd:ied,jsd:jed,km)      ! Delta p at each model level
      real, intent(inout):: ta(isd:ied,jsd:jed,km)      ! Temperature
      real, intent(inout):: qa(isd:ied,jsd:jed,km,nq)   ! Specific humidity & tracers
      real, intent(inout):: u_dt(isd:ied,jsd:jed,km) 
      real, intent(inout):: v_dt(isd:ied,jsd:jed,km) 
      real, intent(inout):: t_dt(is:ie,js:je,km) 
      real, intent(inout):: q_dt(is:ie,js:je,km,nq) 
!---------------------------Local variables-----------------------------
      real, dimension(is:ie,km):: u0, v0, w0, t0, hd, te, gz, tvm, pm
      real q0(is:ie,km,nq), qcon(is:ie,km) 
      real gzh(is:ie)
      real ri, pt1, pt2, ratio, tv, cv
      real qmix, h0, mc, fra, rk, rz, rcv, rdt
      real qs1, qs2, lf, dh, dhs, wz
      integer i, j, k, kk, n, m, iq, km1
      real, parameter:: ustar2 = 1.E-4
      integer :: sphum, liq_wat, rainwat, snowwat, graupel, ice_wat, cld_amt

        rz = rvgas - rdgas          ! rz = zvir * rdgas
        rk = cp_air/rdgas + 1.
        cv = cp_air - rdgas
       rcv = 1./cv

      rdt = 1./ dt

        sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
      liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
      ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
      cld_amt = get_tracer_index (MODEL_ATMOS, 'cld_amt')

      if ( nq.ge.6 ) then
           rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat')
           snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat')
           graupel = get_tracer_index (MODEL_ATMOS, 'graupel')
      endif

!------------------------------------------------------------------------
! The nonhydrostatic pressure changes if there is heating (under constant
! volume and mass is locally conserved).
!------------------------------------------------------------------------
   m = 4
   fra = dt/real(tau)

!$omp parallel do default(shared) private(kk,qcon, q0, t0, u0, v0, w0, h0, pm, gzh, tvm, tv, gz, hd, te, ratio, pt1, pt2,ri, mc)
  do 1000 j=js,je  

    do iq=1, nq
       do k=1,km
          do i=is,ie
             q0(i,k,iq) = qa(i,j,k,iq)
          enddo
       enddo
    enddo

    do k=1,km
       do i=is,ie
          t0(i,k) = ta(i,j,k)
         tvm(i,k) = t0(i,k)*(1.+zvir*q0(i,k,sphum))
          u0(i,k) = ua(i,j,k)
          v0(i,k) = va(i,j,k)
          pm(i,k) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
       enddo
    enddo

    do i=is,ie
       gzh(i) = 0.
    enddo

    if( hydrostatic ) then
       do k=km, 1,-1
          do i=is,ie
                tv  = rdgas*tvm(i,k)
            gz(i,k) = gzh(i) + tv*(1.-pe(i,k,j)/pm(i,k))
            hd(i,k) = cp_air*tvm(i,k)+gz(i,k)+0.5*(u0(i,k)**2+v0(i,k)**2)
             gzh(i) = gzh(i) + tv*(peln(i,k+1,j)-peln(i,k,j))
          enddo
       enddo
    else
       do k=km, 1, -1
          do i=is,ie
             w0(i,k) = w(i,j,k)
             gz(i,k) = gzh(i)  - 0.5*grav*delz(i,j,k)
                 tv  = gz(i,k) + 0.5*(u0(i,k)**2+v0(i,k)**2+w0(i,k)**2)
             hd(i,k) = cp_air*tvm(i,k) + tv
             te(i,k) =     cv*tvm(i,k) + tv
              gzh(i) = gzh(i) - grav*delz(i,j,k)
          enddo
       enddo
    endif

   do n=1,m

      ratio = real(n)/real(m)

      do i=is,ie
         gzh(i) = 0.
      enddo

! Compute total condensate
   if ( nq .le. 3 .or. zvir .lt. 1.e-3 ) then
      do k=1,km
         do i=is,ie
            qcon(i,k) = 0.
         enddo
      enddo
   elseif ( nq .le. 5 ) then
      do k=1,km
         do i=is,ie
            qcon(i,k) = q0(i,k,liq_wat) + q0(i,k,ice_wat)
         enddo
      enddo
   else
      do k=1,km
         do i=is,ie
            qcon(i,k) = q0(i,k,liq_wat) + q0(i,k,ice_wat) + q0(i,k,snowwat) + q0(i,k,rainwat) + q0(i,k,graupel)
         enddo
      enddo
   endif

      do k=km, 2, -1
         km1 = k-1
         do i=is,ie
! Richardson number = g*delz * del_theta/theta / (del_u**2 + del_v**2)
            pt1 = t0(i,km1)/pkz(i,j,km1)*(1.+zvir*q0(i,km1,sphum)-qcon(i,km1))
            pt2 = t0(i,k  )/pkz(i,j,k  )*(1.+zvir*q0(i,k  ,sphum)-qcon(i,k))
            ri = (gz(i,k-1)-gz(i,k))*(pt1-pt2)/( 0.5*(pt1+pt2)*        &
                ((u0(i,k-1)-u0(i,k))**2+(v0(i,k-1)-v0(i,k))**2+ustar2) )
! Dry convective adjustment for K-H instability:
! Compute equivalent mass flux: mc
#ifndef RI_THRESH
            if ( ri < 1. ) then
                 mc = ratio*(1.-max(0.0,ri)) ** 2
#else
            if ( ri < 0.25 ) then
                 mc = (1.-4.*max(0.0,ri)) ** 2
#endif
                 mc = mc*delp(i,j,k-1)*delp(i,j,k)/(delp(i,j,k-1)+delp(i,j,k))
!                if ( .not. hydrostatic ) then
!                   wz = abs(dt*w0(i,k)/delz(i,j,k))
!                   if ( wz > 1. ) then
!                        mc = delp(i,j,k-1)*delp(i,j,k)/(delp(i,j,k-1)+delp(i,j,k))
!                   endif
!                endif
                 do iq=1,nq
                    h0 = mc*(q0(i,k,iq)-q0(i,km1,iq))
                    q0(i,km1,iq) = q0(i,km1,iq) + h0/delp(i,j,km1)
                    q0(i,k  ,iq) = q0(i,k  ,iq) - h0/delp(i,j,k  )
                 enddo
! Recompute qcon
                 if ( nq .le. 3 .or. zvir .lt. 1.e-3 ) then
                    qcon(i,km1) = 0.
                 elseif ( nq .le. 5 ) then
                    qcon(i,km1) = q0(i,km1,liq_wat) + q0(i,km1,ice_wat)
                 else
                    qcon(i,km1) = q0(i,km1,liq_wat) + q0(i,km1,ice_wat) + q0(i,km1,snowwat) + q0(i,km1,rainwat) + q0(i,km1,graupel)
                 endif
! u:
                 h0 = mc*(u0(i,k)-u0(i,k-1))
                 u0(i,k-1) = u0(i,k-1) + h0/delp(i,j,k-1)
                 u0(i,k  ) = u0(i,k  ) - h0/delp(i,j,k  )
! v:
                 h0 = mc*(v0(i,k)-v0(i,k-1))
                 v0(i,k-1) = v0(i,k-1) + h0/delp(i,j,k-1)
                 v0(i,k  ) = v0(i,k  ) - h0/delp(i,j,k  )
              if ( hydrostatic ) then
                 h0 = mc*(hd(i,k)-hd(i,k-1))
                 hd(i,k-1) = hd(i,k-1) + h0/delp(i,j,k-1)
                 hd(i,k  ) = hd(i,k  ) - h0/delp(i,j,k  )
              else
! Total energy
                        h0 = mc*(hd(i,k)-hd(i,k-1))
                 te(i,k-1) = te(i,k-1) + h0/delp(i,j,k-1)
                 te(i,k  ) = te(i,k  ) - h0/delp(i,j,k  )
! w:
                        h0 = mc*(w0(i,k)-w0(i,k-1))
                 w0(i,k-1) = w0(i,k-1) + h0/delp(i,j,k-1)
                 w0(i,k  ) = w0(i,k  ) - h0/delp(i,j,k  )
              endif
            endif    ! Ri condition check
         enddo

!-------------- 
! Retrive Temp:
!--------------
       if ( hydrostatic ) then
         kk = k
         do i=is,ie
            t0(i,kk) = (hd(i,kk)-gzh(i)-0.5*(u0(i,kk)**2+v0(i,kk)**2))  &
                     / ( rk - pe(i,kk,j)/pm(i,kk) )
              gzh(i) = gzh(i) + t0(i,kk)*(peln(i,kk+1,j)-peln(i,kk,j))
            t0(i,kk) = t0(i,kk) / ( rdgas + rz*q0(i,kk,sphum) )
         enddo
         kk = k-1
         do i=is,ie
            t0(i,kk) = (hd(i,kk)-gzh(i)-0.5*(u0(i,kk)**2+v0(i,kk)**2))  &
                     / ((rk-pe(i,kk,j)/pm(i,kk))*(rdgas+rz*q0(i,kk,sphum)))
         enddo
       else
! Non-hydrostatic under constant volume heating/cooling
         do kk=k-1,k
            do i=is,ie
               tv = gz(i,kk) + 0.5*(u0(i,kk)**2+v0(i,kk)**2+w0(i,kk)**2)
               tvm(i,kk) = rcv*(te(i,kk)- tv)
                hd(i,kk) = cp_air*tvm(i,kk) + tv
                t0(i,kk) = tvm(i,kk)/(1.+zvir*q0(i,kk,sphum))
            enddo
         enddo
       endif
      enddo   ! k-loop

   enddo       ! n-loop


!--------------------
   if ( fra < 1. ) then
      do k=1, km
         do i=is,ie
            t0(i,k) = ta(i,j,k) + (t0(i,k) - ta(i,j,k))*fra
            u0(i,k) = ua(i,j,k) + (u0(i,k) - ua(i,j,k))*fra
            v0(i,k) = va(i,j,k) + (v0(i,k) - va(i,j,k))*fra
         enddo
      enddo

      if ( .not. hydrostatic ) then
         do k=1,km
            do i=is,ie
               w0(i,k) = w(i,j,k) + (w0(i,k) - w(i,j,k))*fra
            enddo
         enddo
      endif

      do iq=1,nq
         do k=1,km
            do i=is,ie
               q0(i,k,iq) = qa(i,j,k,iq) + (q0(i,k,iq) - qa(i,j,k,iq))*fra
            enddo
         enddo
      enddo
   endif
!--------------------

   do k=1,km
      do i=is,ie
         u_dt(i,j,k) = rdt*(u0(i,k) - ua(i,j,k))
         v_dt(i,j,k) = rdt*(v0(i,k) - va(i,j,k))
         t_dt(i,j,k) = 0.
           ta(i,j,k) = t0(i,k)   ! *** temperature updated ***
      enddo
      do iq=1,nq
         do i=is,ie
            q_dt(i,j,k,iq) = rdt*(q0(i,k,iq)-qa(i,j,k,iq))
         enddo
      enddo
   enddo

   if ( .not. hydrostatic ) then
      do k=1,km
         do i=is,ie
            w(i,j,k) = w0(i,k)   ! w updated
         enddo
      enddo
   endif

1000 continue


 end subroutine fv_dry_conv



 real function qs1d(t, p, q)
! Based on "moist" mixing ratio, p is the total (dry+vapor) pressure
  real, intent(in):: t, p, q
! Local:
  real es, ap1
  real, parameter:: Tmin=tice - 160.
  integer it

       ap1 = 10.*DIM(t, Tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = table(it) + (ap1-it)*des(it)
      qs1d = esl*es*(1.+zvir*q)/p

  end function qs1d


  subroutine qsmith_init
  integer, parameter:: length=2621 
  integer i

  if( .not. allocated(table) ) then
!                            Generate es table (dT = 0.1 deg. C)

       allocate ( table(length) )
       allocate (  des (length) )

       call qs_table(length, table)

       do i=1,length-1
          des(i) = table(i+1) - table(i)
       enddo
       des(length) = des(length-1)
  endif
 
  end subroutine qsmith_init


  subroutine qsmith(im, km, k1, t, p, q, qs, dqdt)
! input T in deg K; p (Pa)
  integer, intent(in):: im, km, k1
  real, intent(in),dimension(im,km):: t, p, q
  real, intent(out),dimension(im,km):: qs
  real, intent(out), optional:: dqdt(im,km)
! Local:
  real es(im,km)
  real ap1, eps10
  real Tmin
  integer i, k, it

  Tmin = tice-160.
  eps10  = 10.*esl

  if( .not. allocated(table) ) call  qsmith_init
 
      do k=k1,km
         do i=1,im
            ap1 = 10.*DIM(t(i,k), Tmin) + 1.
            ap1 = min(2621., ap1)
            it = ap1
            es(i,k) = table(it) + (ap1-it)*des(it)
            qs(i,k) = esl*es(i,k)*(1.+zvir*q(i,k))/p(i,k)
         enddo
      enddo

      if ( present(dqdt) ) then
      do k=k1,km
           do i=1,im
              ap1 = 10.*DIM(t(i,k), Tmin) + 1.
              ap1 = min(2621., ap1) - 0.5
              it  = ap1
              dqdt(i,k) = eps10*(des(it)+(ap1-it)*(des(it+1)-des(it)))*(1.+zvir*q(i,k))/p(i,k)
           enddo
      enddo
      endif
 
  end subroutine qsmith
 

  subroutine qs_table(n,table)
      integer, intent(in):: n
      real table (n)
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

      do i=1,n
         table(i) = table(i)*0.1
      enddo

 end subroutine qs_table

 subroutine neg_adj3(is, ie, js, je, ng, kbot, hydrostatic,   &
                     peln, delz, pt, dp, qv, ql, qr, qi, qs, qg, qa)

! This is designed for 6-class micro-physics schemes
 use lin_cld_microphys_mod, only: wqsat2_moist
 integer, intent(in):: is, ie, js, je, ng, kbot
 logical, intent(in):: hydrostatic
 real, intent(in):: dp(is-ng:ie+ng,js-ng:je+ng,kbot)  ! total delp-p
 real, intent(in):: delz(is-ng:ie+ng,js-ng:je+ng,kbot)      ! Delta p at each model level
 real, intent(in):: peln(is:ie,kbot+1,js:je)           ! ln(pe)
 real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,kbot)::    &
                                 pt, qv, ql, qr, qi, qs, qg, qa
! Local:
 logical:: sat_adj = .false.
 real, parameter :: c_ice = 2106.
 real, parameter :: c_liq = 4190.
 real, parameter :: t48 = tice - 48.
 real, dimension(is:ie,kbot):: dpk, q2
real, dimension(is:ie,js:je):: t0, pt2, qv2, ql2, qi2, qs2, qr2, qg2, dp2, p2, icpk, lcpk
 real:: cv_air, cv_vap
 real:: dq, qsum, dq1, q_liq, q_sol, cpm, sink, qsw, dwsdt
 integer i, j, k

 cv_air = cp_air - rdgas ! = rdgas * (7/2-1) = 2.5*rdgas=717.68
 cv_vap = 3.*rvgas      ! 1384.5

!$omp parallel do default(shared) private(dq,dq1,qsum,dp2,p2,t0,pt2,qv2,ql2,qi2,qs2,qg2,qr2,lcpk,icpk,qsw,dwsdt,sink,q_liq,q_sol,cpm)
  do k=1, kbot
     do j=js, je
        do i=is, ie
        qv2(i,j) = qv(i,j,k)
        ql2(i,j) = ql(i,j,k)
        qi2(i,j) = qi(i,j,k)
        qs2(i,j) = qs(i,j,k)
        qr2(i,j) = qr(i,j,k)
        qg2(i,j) = qg(i,j,k)
        dp2(i,j) = dp(i,j,k)
        pt2(i,j) = pt(i,j,k)
         t0(i,j) = pt2(i,j)
        enddo
     enddo

     if ( hydrostatic ) then
       do j=js, je
          do i=is, ie
             p2(i,j) = dp2(i,j)/(peln(i,k+1,j)-peln(i,k,j))
             q_liq = ql2(i,j)
             q_sol = qi2(i,j)
             cpm = (1.-(qv2(i,j)+q_liq+q_sol))*cp_air + qv2(i,j)*cp_vapor + q_liq*c_liq + q_sol*c_ice 
             lcpk(i,j) = hlv / cpm
             icpk(i,j) = hlf / cpm
          enddo
       enddo
     else
       do j=js, je
          do i=is, ie
             p2(i,j) = -dp2(i,j)/(grav*delz(i,j,k))*rdgas*pt2(i,j)*(1.+zvir*qv2(i,j))
             q_liq = ql2(i,j)
             q_sol = qi2(i,j)
             cpm = (1.-(qv2(i,j)+q_liq+q_sol))*cv_air + qv2(i,j)*cv_vap + q_liq*c_liq + q_sol*c_ice 
             lcpk(i,j) = hlv / cpm
             icpk(i,j) = hlf / cpm
          enddo
       enddo
     endif

! Fix the negatives:
!-----------
! Ice-phase:
!-----------
    do j=js, je
       do i=is, ie
        qsum = qi2(i,j) + qs2(i,j)
        if ( qsum > 0. ) then
             if ( qi2(i,j) < 0. ) then
                  qi2(i,j) = 0.
                  qs2(i,j) = qsum
             elseif ( qs2(i,j) < 0. ) then
                  qs2(i,j) = 0.
                  qi2(i,j) = qsum
             endif
        else
! borrow froom graupel
             qi2(i,j) = 0.
             qs2(i,j) = 0.
             qg2(i,j) = qg2(i,j) + qsum
        endif

! At this stage qi and qs should be positive definite
! If graupel < 0 then borrow from qs then qi
        if ( qg2(i,j) < 0. ) then
             dq = min( qs2(i,j), -qg2(i,j) )
             qs2(i,j) = qs2(i,j) - dq
             qg2(i,j) = qg2(i,j) + dq
             if ( qg2(i,j) < 0. ) then
! if qg is still negative
                  dq = min( qi2(i,j), -qg2(i,j) )
                  qi2(i,j) = qi2(i,j) - dq
                  qg2(i,j) = qg2(i,j) + dq
             endif
        endif

! If qg is still negative then borrow from rain water: phase change
        if ( qg2(i,j)<0. .and. qr2(i,j)>0. ) then
             dq = min( qr2(i,j), -qg2(i,j) )
             qg2(i,j) = qg2(i,j) + dq
             qr2(i,j) = qr2(i,j) - dq
             pt2(i,j) = pt2(i,j) + dq*icpk(i,j)  ! conserve total energy
        endif
! If qg is still negative then borrow from cloud water: phase change
        if ( qg2(i,j)<0. .and. ql2(i,j)>0. ) then
             dq = min( ql2(i,j), -qg2(i,j) )
             qg2(i,j) = qg2(i,j) + dq
             ql2(i,j) = ql2(i,j) - dq
             pt2(i,j) = pt2(i,j) + dq*icpk(i,j)
        endif
! Last resort; borrow from water vapor
        if ( qg2(i,j)<0. .and. qv2(i,j)>0. ) then
             dq = min( 0.999*qv2(i,j), -qg2(i,j) )
             qg2(i,j) = qg2(i,j) + dq
             qv2(i,j) = qv2(i,j) - dq
             pt2(i,j) = pt2(i,j) + dq*(icpk(i,j)+lcpk(i,j))
        endif

!--------------
! Liquid phase:
!--------------
        qsum = ql2(i,j) + qr2(i,j)
        if ( qsum > 0. ) then
             if ( qr2(i,j) < 0. ) then
                  qr2(i,j) = 0.
                  ql2(i,j) = qsum
             elseif ( ql2(i,j) < 0. ) then
                  ql2(i,j) = 0.
                  qr2(i,j) = qsum
             endif
        else
          ql2(i,j) = 0.
          qr2(i,j) = qsum     ! rain water is still negative
! fill negative rain with qg first
          dq = min( max(0.0, qg2(i,j)), -qr2(i,j) )
          qr2(i,j) = qr2(i,j) + dq
          qg2(i,j) = qg2(i,j) - dq
          pt2(i,j) = pt2(i,j) - dq*icpk(i,j)
          if ( qr(i,j,k) < 0. ) then
! fill negative rain with available qi & qs (cooling)
               dq = min( qi2(i,j)+qs2(i,j), -qr2(i,j) )
               qr2(i,j) = qr2(i,j) + dq
               dq1 = min( dq, qs2(i,j) )
               qs2(i,j) = qs2(i,j) - dq1
               qi2(i,j) = qi2(i,j) + dq1 - dq 
               pt2(i,j) = pt2(i,j) - dq*icpk(i,j)
          endif
! fix negative rain water with available vapor
          if ( qr2(i,j)<0. .and. qv2(i,j)>0. ) then
               dq = min( 0.999*qv2(i,j), -qr2(i,j) )
               qv2(i,j) = qv2(i,j) - dq
               qr2(i,j) = qr2(i,j) + dq
               pt2(i,j) = pt2(i,j) + dq*lcpk(i,j)
          endif
        endif
     enddo
   enddo

!******************************************
! Fast moist physics: Saturation adjustment
!******************************************
 if ( sat_adj ) then

   do j=js, je
     do i=is, ie
! Melting of cloud ice into cloud water ********
        if ( qi2(i,j)>1.e-8 .and. pt2(i,j) > tice ) then
           sink = min( qi2(i,j), (pt2(i,j)-tice)/icpk(i,j) )
           ql2(i,j) = ql2(i,j) + sink
           qi2(i,j) = qi2(i,j) - sink
           pt2(i,j) = pt2(i,j) - sink*icpk(i,j)
        endif

! vapor <---> liquid water --------------------------------
        qsw = wqsat2_moist(pt2(i,j), qv2(i,j), p2(i,j), dwsdt)
        sink = min( ql2(i,j), (qsw-qv2(i,j))/(1.+lcpk(i,j)*dwsdt) )
        qv2(i,j) = qv2(i,j) + sink
        ql2(i,j) = ql2(i,j) - sink
        pt2(i,j) = pt2(i,j) - sink*lcpk(i,j)
!-----------------------------------------------------------

! freezing of cloud water ********
        if( ql2(i,j)>1.e-8 .and. pt2(i,j) < t48 ) then
! Enforce complete freezing below t_00 (-48 C)
            sink = min( ql2(i,j), (t48-pt2(i,j))/icpk(i,j) )
            ql2(i,j) = ql2(i,j) - sink
            qi2(i,j) = qi2(i,j) + sink
            pt2(i,j) = pt2(i,j) + sink*icpk(i,j)
        endif ! significant ql existed
     enddo
   enddo

!----------------------------------------------------------------
! Convert back to moist mixing ratios:
   do j=js, je
     do i=is, ie
        qv(i,j,k) = qv2(i,j)
        ql(i,j,k) = ql2(i,j)
        qi(i,j,k) = qi2(i,j)
        qs(i,j,k) = qs2(i,j)
        qr(i,j,k) = qr2(i,j)
        qg(i,j,k) = qg2(i,j)
        pt(i,j,k) = pt2(i,j)
     enddo
   enddo
 endif

 enddo

!$omp parallel do default(shared) private(dpk, q2)
 do j=js, je
! Graupel:
    do k=1,kbot
       do i=is,ie
          dpk(i,k) = dp(i,j,k)
           q2(i,k) = qg(i,j,k)
       enddo
    enddo
    call fillq(ie-is+1, kbot, q2, dpk)
    do k=1,kbot
       do i=is,ie
          qg(i,j,k) = q2(i,k)
       enddo
    enddo
! Rain water:
    do k=1,kbot
       do i=is,ie
          q2(i,k) = qr(i,j,k)
       enddo
    enddo
    call fillq(ie-is+1, kbot, q2, dpk)
    do k=1,kbot
       do i=is,ie
          qr(i,j,k) = q2(i,k)
       enddo
    enddo
 enddo

!-----------------------------------
! Fix water vapor
!-----------------------------------
! Top layer: borrow from below
    k = 1
!$omp parallel do default(shared)
   do j=js, je
       do i=is, ie
          if( qv(i,j,k) < 0. ) then
              qv(i,j,k+1) = qv(i,j,k+1) + qv(i,j,k)*dp(i,j,k)/dp(i,j,k+1)
              qv(i,j,k  ) = 0.
          endif
     enddo
   enddo

! this OpenMP do-loop cannot be parallelized with recursion on k/k-1
!rab!$omp parallel do default(shared) private(dq)
 do k=2,kbot-1
    do j=js, je
       do i=is, ie
          if( qv(i,j,k) < 0. .and. qv(i,j,k-1) > 0. ) then
              dq = min(-qv(i,j,k)*dp(i,j,k), qv(i,j,k-1)*dp(i,j,k-1))
              qv(i,j,k-1) = qv(i,j,k-1) - dq/dp(i,j,k-1) 
              qv(i,j,k  ) = qv(i,j,k  ) + dq/dp(i,j,k  ) 
          endif
          if( qv(i,j,k) < 0. ) then
              qv(i,j,k+1) = qv(i,j,k+1) + qv(i,j,k)*dp(i,j,k)/dp(i,j,k+1)
              qv(i,j,k  ) = 0.
          endif
       enddo
    enddo
  enddo
 
! Bottom layer; Borrow from above
!$omp parallel do default(shared) private(dq)
  do j=js, je
     do i=is, ie
     if( qv(i,j,kbot) < 0. ) then
         do k=kbot-1,1,-1
            if ( qv(i,j,kbot)>=0. ) goto 123
            if ( qv(i,j,k) > 0. ) then
                 dq = min(-qv(i,j,kbot)*dp(i,j,kbot), qv(i,j,k)*dp(i,j,k))
                 qv(i,j,k   ) = qv(i,j,k   ) - dq/dp(i,j,k) 
                 qv(i,j,kbot) = qv(i,j,kbot) + dq/dp(i,j,kbot) 
            endif
         enddo   ! k-loop
123      continue
     endif
     enddo ! i-loop
 enddo   ! j-loop

!-----------------------------------
! Fix negative cloud fraction
!-----------------------------------
! this OpenMP do-loop cannot be parallelized by the recursion on k/k+1
!rab!$omp parallel do default(shared) 
 do k=1,kbot-1
    do j=js, je
       do i=is, ie
          if( qa(i,j,k) < 0. ) then
              qa(i,j,k+1) = qa(i,j,k+1) + qa(i,j,k)*dp(i,j,k)/dp(i,j,k+1)
              qa(i,j,k  ) = 0.
          endif
     enddo
   enddo
 enddo
 
! Bottom layer; Borrow from above
!$omp parallel do default(shared) private(dq)
  do j=js, je
     do i=is, ie
        if( qa(i,j,kbot) < 0. .and. qa(i,j,kbot-1)>0.) then
            dq = min(-qa(i,j,kbot)*dp(i,j,kbot), qa(i,j,kbot-1)*dp(i,j,kbot-1))
            qa(i,j,kbot-1) = qa(i,j,kbot-1) - dq/dp(i,j,kbot-1) 
            qa(i,j,kbot  ) = qa(i,j,kbot  ) + dq/dp(i,j,kbot  ) 
        endif
! if qa is still < 0
        qa(i,j,kbot) = max(0., qa(i,j,kbot))
   enddo
 enddo


 end subroutine neg_adj3

 subroutine fillq(im, km, q, dp)
! Aggresive 1D filling algorithm for qr and qg
 integer, intent(in):: im, km
 real, intent(inout), dimension(im,km):: q, dp
 integer:: i, k
 real:: sum1, sum2, dq

 do 500 i=1,im
    sum1 = 0.
    do k=1,km
       if ( q(i,k)>0. ) then
            sum1 = sum1 + q(i,k)*dp(i,k)
       endif
    enddo
    if ( sum1<1.E-12  ) goto 500
    sum2 = 0.
    do k=km,1,-1
       if ( q(i,k)<0.0 .and. sum1>0. ) then
            dq = min( sum1, -q(i,k)*dp(i,k) )
            sum1 = sum1 - dq
            sum2 = sum2 + dq
            q(i,k) = q(i,k) + dq/dp(i,k)
       endif
    enddo
    do k=km,1,-1
       if ( q(i,k)>0.0 .and. sum2>0. ) then
            dq = min( sum2, q(i,k)*dp(i,k) )
            sum2 = sum2 - dq
            q(i,k) = q(i,k) - dq/dp(i,k)
       endif
    enddo
500  continue

 end subroutine fillq

end module fv_sg_mod
