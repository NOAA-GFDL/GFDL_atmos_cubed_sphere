module fv_sg_mod

!-----------------------------------------------------------------------
! FV sub-grid mixing
!-----------------------------------------------------------------------
use constants_mod, only: rdgas, rvgas, cp_air, hlv
! use mp_mod,        only: gid

implicit none
private

public  fv_sg_conv, qsmith

  logical:: qs_initialized = .false.

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
      real, dimension(is:ie,km):: u0, v0, t0, gz, qsat, hd, pm
      real q0(is:ie,km,nq) 
      real gzh(is:ie)
      real pbot
      real rdt, zvir, dh, dq, tv, tvm, qmix, h0, mc, fra, rdm, rl, rk, rz
      integer mcond, kcond
      integer i, j, k, n, m, iq

      zvir = rvgas / rdgas - 1.     ! = 0.607789855
        rz = rvgas - rdgas          ! rz = zvir * rdgas
        rl = 1. / hlv
        rk = cp_air/rdgas + 1.

! Narrow the convective domain
      mcond = 2
      do k=2,km
         pbot = ak(k+1) + bk(k+1)*1.E5
         if( pbot > 50.E2 ) then
             mcond = k
             go to 150
         endif
      enddo
150 continue

    rdt = 1. / dt
    fra = dt/real(tau)

    m = 2
    if ( fra >= 1. ) m = 2 + fra

!------------
! Compute gz
!------------
  do 1000 j=js,je       ! this main loop can be OpneMPed in j

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
              tvm = t0(i,k)*(1.+zvir*q0(i,k,1))
              tv  = rdgas*tvm
          pm(i,k) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
          gz(i,k) = gzh(i) + tv*(1.-pe(i,k,j)/pm(i,k))
          hd(i,k) = cp_air*tvm + gz(i,k) + 0.5*(u0(i,k)**2+v0(i,k)**2)
           gzh(i) = gzh(i) + tv*(peln(i,k+1,j)-peln(i,k,j))
       enddo
    enddo

   do n=1,m

      kcond = min(km-2, mcond + m - n)

      call qsmith(ie-is+1, km, kcond, t0, pm, qsat)

      do k=km,kcond+1,-1
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
                if ( dq > 0. ) then
! Moist Mixing/convection
                    qmix = (q0(i,k-1,1)*delp(i,j,k-1)+q0(i,k,1)*delp(i,j,k)) /  &
                           (delp(i,j,k-1)+delp(i,j,k))
                    dq = min( q0(i,k,1) - qmix, dq )
                  if( dq > 1.E-7 ) then
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
! Compute tendencies:
!--------------------
   do k=1,mcond-1
      do i=is,ie
         u_dt(i,j,k) = 0.
         v_dt(i,j,k) = 0.
         t_dt(i,j,k) = 0.
      enddo
   enddo
   do k=mcond,km
      do i=is,ie
         u_dt(i,j,k) = (u0(i,k) - ua(i,j,k)) * rdt
         v_dt(i,j,k) = (v0(i,k) - va(i,j,k)) * rdt
         t_dt(i,j,k) = (t0(i,k) - ta(i,j,k)) * rdt
      enddo
   enddo
   do iq=1,nq
      do k=1,mcond-1
         do i=is,ie
            q_dt(i,j,k,iq) = 0.
         enddo
      enddo
      do k=mcond,km
         do i=is,ie
            q_dt(i,j,k,iq) = (q0(i,k,iq) - qa(i,j,k,iq))*rdt
         enddo
      enddo
   enddo

1000 continue

  end subroutine fv_sg_conv


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
      real, parameter:: Tice = 273.16
      real:: dt=0.1
      integer i, k, it

      Tmin = Tice-160.
      Tmax = Tice+102.

      oms = 1. - esl
 
      if( .not. qs_initialized ) then
! Generate es table (dT = 0.1 deg. C)
           call qs_table(table)
           do i=1,2620
              des(i) = table(i+1) - table(i)
           enddo
           des(2621) = des(2620)
           qs_initialized = .true.
!          if(gid==0) write(*,*) 'qsmith initialized' 
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
              ap1 = 10.*DIM(t(i,k), Tmin) + 1
              ap1 = min(2621., ap1) - 0.5
              it  = ap1
              dqdt(i,k) = 10.*esl*(des(it) + (ap1-it)*(des(it+1)-des(it)))/p(i,k)
           enddo
      enddo
      endif
 
  end subroutine qsmith
 

  subroutine qs_table(table)
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

  end subroutine qs_table


end module fv_sg_mod
