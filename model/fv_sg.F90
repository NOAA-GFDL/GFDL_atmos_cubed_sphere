module fv_sg_mod

!-----------------------------------------------------------------------
! FV sub-grid mixing
!-----------------------------------------------------------------------
use constants_mod, only: rdgas, rvgas, cp_air, hlv, hlf, kappa
! use fv_mp_mod,        only: gid

implicit none
private

public  fv_sg_conv, qsmith_init, qsmith, qs1d, neg_adj3

  real, parameter:: Tice = 273.16
  real, allocatable:: table(:),des(:)

contains

 subroutine fv_sg_conv(is, ie, js, je, isd, ied, jsd, jed,    &
                        isc, iec, jsc, jec,  km, nq, dt, tau,  &
                        delp, pe, peln, ta, qa, ua, va,        &
                        u_dt, v_dt, t_dt, q_dt, ak, bk)
! Non-precipitating sub-grid scale convective adjustment-mixing
!-------------------------------------------
      integer, intent(in):: is, ie, js, je, km, nq
      integer, intent(in):: isc, iec, jsc, jec
      integer, intent(in):: isd, ied, jsd, jed
      integer, intent(in):: tau         ! Relaxation time scale
      real, intent(in):: dt             ! model time step
      real, intent(in):: ak(km+1), bk(km+1)
      real, intent(in)::   pe(isc-1:iec+1,km+1,jsc-1:jec+1) 
      real, intent(in):: peln(isc  :iec,  km+1,jsc  :jec)
      real, intent(in):: delp(isd:ied,jsd:jed,km)      ! Delta p at each model level
      real, intent(in)::   ta(isd:ied,jsd:jed,km)      ! Temperature
      real, intent(in)::   qa(isd:ied,jsd:jed,km,nq)   ! Specific humidity & tracers
      real, intent(in)::   ua(isd:ied,jsd:jed,km)
      real, intent(in)::   va(isd:ied,jsd:jed,km)
! Output:
      real, intent(out):: u_dt(isd:ied,jsd:jed,km) 
      real, intent(out):: v_dt(isd:ied,jsd:jed,km) 

      real, intent(out):: t_dt(isc:iec,jsc:jec,km) 
      real, intent(out):: q_dt(isc:iec,jsc:jec,km,nq) 
!---------------------------Local variables-----------------------------
      real, parameter:: qmin1 = 1.E-7
      real, parameter:: qmin2 = 1.E-8
      real, parameter:: pmin  =  5000.
!     real, parameter:: pmin  = 40000.
      real, dimension(is:ie,km):: tvm, u0, v0, t0, gz, qsat, hd, pm, pkz
      real, dimension(is:ie,km+1):: pk
      real q0(is:ie,km,nq) 
      real gzh(is:ie)
      real pbot, ri, pt1, pt2
      real rdt, zvir, dh, dq, tv, qmix, h0, mc, fra, rdm, rl, rk, rz
      integer mcond, kcond
      integer i, j, k, n, m, iq
      real, parameter:: ustar2 = 1.E-8

      zvir = rvgas / rdgas - 1.     ! = 0.607789855
        rz = rvgas - rdgas          ! rz = zvir * rdgas
        rl = 1. / hlv
        rk = cp_air/rdgas + 1.

! Narrow the convective domain
      mcond = 2
      do k=2,km
         pbot = ak(k+1) + bk(k+1)*1.E5
         if( pbot > pmin ) then
             mcond = k
             go to 150
         endif
      enddo
150 continue

    rdt = 1. / dt
    fra = dt/real(tau)


!------------------------------------------------------------------------
! The nonhydrostatic pressure changes if there is heating (under constant
! volume and mass is locally conserved).
!------------------------------------------------------------------------

    m = 3
    if ( fra >= 1. ) m = m + fra

!------------
! Compute gz: center 
!------------
  do 1000 j=js,je       ! this main loop can be OpneMPed in j

    do k=mcond,km+1
       do i=is,ie
          pk(i,k) = pe(i,k,j)**kappa
       enddo
    enddo

    do k=mcond,km
       do i=is,ie
          u0(i,k) = ua(i,j,k)
          v0(i,k) = va(i,j,k)
          t0(i,k) = ta(i,j,k)
       enddo
    enddo
    do iq=1,nq
       do k=mcond,km
          do i=is,ie
             q0(i,k,iq) = qa(i,j,k,iq)
          enddo
       enddo
    enddo

    do i=is,ie
       gzh(i) = 0.
    enddo
    do k=km,mcond,-1
       do i=is,ie
          tvm(i,k) = t0(i,k)*(1.+zvir*q0(i,k,1))
              tv  = rdgas*tvm(i,k)
          pm(i,k) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
          gz(i,k) = gzh(i) + tv*(1.-pe(i,k,j)/pm(i,k))
          hd(i,k) = cp_air*tvm(i,k) + gz(i,k) + 0.5*(u0(i,k)**2+v0(i,k)**2)
           gzh(i) = gzh(i) + tv*(peln(i,k+1,j)-peln(i,k,j))
         pkz(i,k) = (pk(i,k+1)-pk(i,k))/(kappa*(peln(i,k+1,j)-peln(i,k,j)))
       enddo
    enddo

   do n=1,m

      kcond = min(km-2, mcond + m - n)

      call qsmith(ie-is+1, km, kcond, t0, pm, qsat)

      do k=km,kcond+1,-1

!------------------------------------------------------------------------------
         do i=is,ie
! Richardson number = g*delz * theta / ( del_theta * (del_u**2 + del_v**2) )
            pt1 = tvm(i,k-1)/pkz(i,k-1)
            pt2 = tvm(i,k  )/pkz(i,k  )
            ri = (gz(i,k-1)-gz(i,k))*(pt1-pt2)/( 0.5*(pt1+pt2)*        &
                ((u0(i,k-1)-u0(i,k))**2+(v0(i,k-1)-v0(i,k))**2+ustar2) )
            if ( ri < 0.25 ) then
! Dry convective adjustment for K-H instability:
! Compute equivalent mass flux: mc
! Ri=1: no mixing
! Ri=0: complete mixing
                 mc = (1.-4.*max(0.,ri))*delp(i,j,k-1)*delp(i,j,k)/(delp(i,j,k-1)+delp(i,j,k))
                 do iq=1,nq
                    h0 = mc*(q0(i,k,iq)-q0(i,k-1,iq))
                    q0(i,k-1,iq) = q0(i,k-1,iq) + h0/delp(i,j,k-1)
                    q0(i,k  ,iq) = q0(i,k  ,iq) - h0/delp(i,j,k  )
                 enddo
! h:
                 h0 = mc*(hd(i,k)-hd(i,k-1))
                 hd(i,k-1) = hd(i,k-1) + h0/delp(i,j,k-1)
                 hd(i,k  ) = hd(i,k  ) - h0/delp(i,j,k  )
! u:
                 h0 = mc*(u0(i,k)-u0(i,k-1))
                 u0(i,k-1) = u0(i,k-1) + h0/delp(i,j,k-1)
                 u0(i,k  ) = u0(i,k  ) - h0/delp(i,j,k  )
! v:
                 h0 = mc*(v0(i,k)-v0(i,k-1))
                 v0(i,k-1) = v0(i,k-1) + h0/delp(i,j,k-1)
                 v0(i,k  ) = v0(i,k  ) - h0/delp(i,j,k  )
            endif
         enddo
!------------------------------------------------------------------------------

         do i=is,ie
            dh = hd(i,k) - hd(i,k-1)
            if ( dh > 0. ) then
! Dry convective adjustment:
                rdm = 1. / (delp(i,j,k-1) + delp(i,j,k))
                h0 = (delp(i,j,k-1)*hd(i,k-1)+delp(i,j,k)*hd(i,k)) * rdm
                hd(i,k-1) = h0
                hd(i,k  ) = h0
                do iq=1,nq
                   qmix = (q0(i,k-1,iq)*delp(i,j,k-1)+q0(i,k,iq)*delp(i,j,k))*rdm
                   q0(i,k-1,iq) = qmix
                   q0(i,k  ,iq) = qmix
                enddo
! Momentum mixing:
                qmix = (u0(i,k-1)*delp(i,j,k-1)+u0(i,k)*delp(i,j,k))*rdm
                u0(i,k-1) = qmix
                u0(i,k  ) = qmix
                qmix = (v0(i,k-1)*delp(i,j,k-1)+v0(i,k)*delp(i,j,k))*rdm
                v0(i,k-1) = qmix
                v0(i,k  ) = qmix
            else
                dq = q0(i,k,1) - qsat(i,k-1) + dh*rl
                if ( dq > qmin1 ) then
! Moist Mixing/convection
                    qmix = (q0(i,k-1,1)*delp(i,j,k-1)+q0(i,k,1)*delp(i,j,k)) /  &
                           (delp(i,j,k-1)+delp(i,j,k))
                    dq = min( q0(i,k,1) - qmix, dq )
                  if( dq > qmin2 ) then
! Compute equivalent mass flux: mc
                    mc = dq/(q0(i,k,1)-q0(i,k-1,1))*delp(i,j,k)
                    do iq=1,nq
                       h0 = mc*(q0(i,k,iq)-q0(i,k-1,iq))
                       q0(i,k-1,iq) = q0(i,k-1,iq) + h0/delp(i,j,k-1)
                       q0(i,k  ,iq) = q0(i,k  ,iq) - h0/delp(i,j,k  )
                    enddo
! h:
                    h0 = mc*dh    ! dh < 0
                    hd(i,k-1) = hd(i,k-1) + h0/delp(i,j,k-1)
                    hd(i,k  ) = hd(i,k  ) - h0/delp(i,j,k  )
! u:
                    h0 = mc*(u0(i,k)-u0(i,k-1))
                    u0(i,k-1) = u0(i,k-1) + h0/delp(i,j,k-1)
                    u0(i,k  ) = u0(i,k  ) - h0/delp(i,j,k  )
! v:
                    h0 = mc*(v0(i,k)-v0(i,k-1))
                    v0(i,k-1) = v0(i,k-1) + h0/delp(i,j,k-1)
                    v0(i,k  ) = v0(i,k  ) - h0/delp(i,j,k  )
                  endif
               endif
            endif
         enddo
      enddo

!-------------
! Retrive Temp
!-------------
      do i=is,ie
         gzh(i) = 0.
      enddo
      do k=km,kcond,-1
         do i=is,ie
            t0(i,k) = (hd(i,k)-gzh(i)-0.5*(u0(i,k)**2+v0(i,k)**2))   &
                     / ( rk - pe(i,k,j)/pm(i,k) )
            gzh(i) = gzh(i) + t0(i,k)*(peln(i,k+1,j)-peln(i,k,j))
            t0(i,k) = t0(i,k) / ( rdgas + rz*q0(i,k,1) )
         enddo
      enddo

   enddo       ! n-loop

   if ( fra < 1. ) then
      do k=mcond,km
         do i=is,ie
            t0(i,k) = ta(i,j,k) + (t0(i,k) - ta(i,j,k))*fra
            u0(i,k) = ua(i,j,k) + (u0(i,k) - ua(i,j,k))*fra
            v0(i,k) = va(i,j,k) + (v0(i,k) - va(i,j,k))*fra
         enddo
      enddo
      do iq=1,nq
         do k=mcond,km
            do i=is,ie
               q0(i,k,iq) = qa(i,j,k,iq) + (q0(i,k,iq) - qa(i,j,k,iq))*fra
            enddo
         enddo
      enddo
   endif

!--------------------
! Update fields:
!--------------------
   do k=1,mcond-1
      do i=is,ie
         u_dt(i,j,k) = ua(i,j,k)
         v_dt(i,j,k) = va(i,j,k)
         t_dt(i,j,k) = ta(i,j,k)
      enddo
   enddo
   do k=mcond,km
      do i=is,ie
         u_dt(i,j,k) = u0(i,k)
         v_dt(i,j,k) = v0(i,k)
         t_dt(i,j,k) = t0(i,k)
      enddo
   enddo

   do iq=1,nq
      do k=1,mcond-1
         do i=is,ie
            q_dt(i,j,k,iq) = qa(i,j,k,iq)
         enddo
      enddo
      do k=mcond,km
         do i=is,ie
            q_dt(i,j,k,iq) = q0(i,k,iq)
         enddo
      enddo
   enddo

1000 continue

 end subroutine fv_sg_conv

 real function qs1d(t, p, dqdt)
! This routine is for use in mp_lin; which uses dry mixing ratio
  real, intent(in):: t, p
  real, intent(out):: dqdt
! Local:
  real q, es, ap1
  real, parameter:: esl = 0.621971831
  real, parameter:: oms= esl - 1.
  real, parameter:: Tmin=Tice - 160.
  integer it

       ap1 = 10.*DIM(t, Tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = table(it) + (ap1-it)*des(it)
      qs1d = esl*es/max(es, p-es)
! dsq/dt
      it  = ap1 - 0.5
! the following form is consistent with mp_lin
      dqdt = 10.*esl*(des(it) + (ap1-it)*(des(it+1)-des(it)))/max(es, p-es)

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


  subroutine qsmith(im, km, k1, t, p, q, dqdt)
! input T in deg K; p (Pa)
  integer, intent(in):: im, km, k1
  real, intent(in),dimension(im,km):: t, p
  real, intent(out),dimension(im,km):: q
  real, intent(out), optional:: dqdt(im,km)
! Local:
  real es(im,km)
  real ap1
  real Tmin, oms
  real, parameter:: esl = 0.621971831
  real:: dt=0.1
  integer i, k, it

  Tmin = Tice-160.

  oms = 1. - esl
 
  if( .not. allocated(table) ) then
       call  qsmith_init
  endif
 
      do k=k1,km
         do i=1,im
            ap1 = 10.*DIM(t(i,k), Tmin) + 1.
            ap1 = min(2621., ap1)
            it = ap1
            es(i,k) = table(it) + (ap1-it)*des(it)
            q(i,k) = esl*es(i,k) / max(es(i,k), p(i,k)-oms*es(i,k))
            q(i,k) = min(1., q(i,k)) 
         enddo
      enddo

      if ( present(dqdt) ) then
      do k=k1,km
           do i=1,im
              ap1 = 10.*DIM(t(i,k), Tmin) + 1.
              ap1 = min(2621., ap1) - 0.5
              it  = ap1
              dqdt(i,k) = 10.*esl*(des(it) + (ap1-it)*(des(it+1)-des(it)))/p(i,k)
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

 subroutine neg_adj3(is, ie, js, je, ng, kbot,      &
                     pt, dp, qv, ql, qr, qi, qs, qg)

! This is designed for 6-class micro-physics schemes
 integer, intent(in):: is, ie, js, je, ng, kbot
 real, intent(in):: dp(is-ng:ie+ng,js-ng:je+ng,kbot)
 real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,kbot)::    &
                                 pt, qv, ql, qr, qi, qs, qg
! Local:
 real lcp, icp
 real dq
 integer i, j, k

 lcp = hlv / cp_air
 icp = hlf / cp_air

 do k=1, kbot
    do j=js, je
       do i=is, ie
!-----------
! Ice-phase:
!-----------
! if ice<0 borrow from snow
          if( qi(i,j,k) < 0. ) then
              qs(i,j,k) = qs(i,j,k) + qi(i,j,k)
              qi(i,j,k) = 0.
          endif
! if snow<0 borrow from graupel
          if( qs(i,j,k) < 0. ) then
              qg(i,j,k) = qg(i,j,k) + qs(i,j,k)
              qs(i,j,k) = 0.
          endif
! If graupel < 0 then borrow from rain
          if ( qg(i,j,k) < 0. ) then
               qr(i,j,k) = qr(i,j,k) + qg(i,j,k)
               pt(i,j,k) = pt(i,j,k) - qg(i,j,k)*icp   ! heating
               qg(i,j,k) = 0.
          endif

! Liquid phase:
! Fix negative rain by borrowing from cloud water
          if ( qr(i,j,k) < 0. ) then
               ql(i,j,k) = ql(i,j,k) + qr(i,j,k)
               qr(i,j,k) = 0.
          endif
! fix negative cloud water with vapor
          if ( ql(i,j,k) < 0. ) then
               qv(i,j,k) = qv(i,j,k) + ql(i,j,k)
               pt(i,j,k) = pt(i,j,k) - ql(i,j,k)*lcp
               ql(i,j,k) = 0.
          endif
     enddo
   enddo
 enddo

!-----------------------------------
! Fix water vapor; borrow from below
!-----------------------------------
 do k=1,kbot-1
    do j=js, je
       do i=is, ie
          if( qv(i,j,k) < 0. ) then
              qv(i,j,k+1) = qv(i,j,k+1) + qv(i,j,k)*dp(i,j,k)/dp(i,j,k+1)
              qv(i,j,k  ) = 0.
          endif
     enddo
   enddo
 enddo
 
! Bottom layer; Borrow from above
  do j=js, je
     do i=is, ie
        if( qv(i,j,kbot) < 0. .and. qv(i,j,kbot-1)>0.) then
            dq = min(-qv(i,j,kbot)*dp(i,j,kbot), qv(i,j,kbot-1)*dp(i,j,kbot-1))
            qv(i,j,kbot-1) = qv(i,j,kbot-1) - dq/dp(i,j,kbot-1) 
            qv(i,j,kbot  ) = qv(i,j,kbot  ) + dq/dp(i,j,kbot  ) 
        endif
! if qv is still < 0
   enddo
 enddo

 end subroutine neg_adj3

end module fv_sg_mod
