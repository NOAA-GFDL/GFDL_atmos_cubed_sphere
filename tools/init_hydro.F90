! $Id: init_hydro.F90,v 14.0 2007/03/15 21:58:53 fms Exp $

module init_hydro

      use constants_mod, only: grav, rdgas
      use grid_utils,    only: g_sum
      use grid_tools,    only: area
      use mp_mod,        only: gid, masterproc
!     use fv_diagnostics_mod, only: prt_maxmin

      implicit none
      private

      public :: p_var, hydro_eq

contains

!-------------------------------------------------------------------------------
 subroutine p_var(km, ifirst, ilast, jfirst, jlast, ptop, ptop_min,    &
                  delp, delz, pt, ps,  pe, peln, pk, pkz, cappa, q, ng, nq,    &
                  dry_mass, adjust_dry_mass, mountain, full_phys, hydrostatic, ktop)
               
! Given (ptop, delp) computes (ps, pk, pe, peln, pkz)
! Input:
   integer,  intent(in):: km
   integer,  intent(in):: ifirst, ilast            ! Longitude strip
   integer,  intent(in):: jfirst, jlast            ! Latitude strip
   integer,  intent(in):: nq
   integer,  intent(in):: ng
   integer,  intent(in):: ktop
   logical, intent(in):: adjust_dry_mass, mountain, full_phys, hydrostatic
   real, intent(in):: dry_mass, cappa, ptop, ptop_min
   real, intent(in   )::   pt(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng, km)
   real, intent(in   ):: delz(ifirst   :ilast   ,jfirst   :jlast   , km)
   real, intent(inout):: delp(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng, km)
   real, intent(inout)::    q(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng, km, nq)
! Output:
   real, intent(out) ::   ps(ifirst-ng:ilast+ng, jfirst-ng:jlast+ng)
   real, intent(out) ::   pk(ifirst:ilast, jfirst:jlast, km+1)
   real, intent(out) ::   pe(ifirst-1:ilast+1,km+1,jfirst-1:jlast+1) ! Ghosted Edge pressure
   real, intent(out) :: peln(ifirst:ilast, km+1, jfirst:jlast)    ! Edge pressure
   real, intent(out) ::  pkz(ifirst:ilast, jfirst:jlast, km)

! Local
   real pek
   real lnp
   real ak1, rdg
   integer i, j, k

   pek = ptop ** cappa
   ak1 = (cappa + 1.) / cappa

   do j=jfirst,jlast
      do i=ifirst,ilast
         pe(i,1,j) = ptop
         pk(i,j,1) = pek
      enddo

      do k=2,km+1
         do i=ifirst,ilast
            pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
            peln(i,k,j) = log(pe(i,k,j))
            pk(i,j,k) = pe(i,k,j)**cappa
         enddo
      enddo

      do i=ifirst,ilast
         ps(i,j) = pe(i,km+1,j)
      enddo

      if( ptop < ptop_min ) then
!---- small ptop modification -------------
          do i=ifirst,ilast
             peln(i,1,j) = peln(i,2,j) - ak1
          enddo
      else
             lnp = log( ptop )
          do i=ifirst,ilast
             peln(i,1,j) = lnp
          enddo
      endif

      if ( hydrostatic ) then
         do k=1,km
            do i=ifirst,ilast
               pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k))/(cappa*(peln(i,k+1,j)-peln(i,k,j)))
            enddo
         enddo
      endif
   enddo

!--------
! Caution:
!------------------------------------------------------------------
! The following form is the same as in "update_fv_phys.F90"
! Therefore, restart reproducibility is only enforced in diabatic cases
!------------------------------------------------------------------
! For adiabatic runs, this form is not exactly the same as in mapz_module;
! Therefore, rounding differences will occur with restart!
   if ( .not.hydrostatic ) then

      if ( ktop>1 ) then
! Compute pkz using hydrostatic formular:
         do k=1,ktop-1
            do j=jfirst,jlast
            do i=ifirst,ilast
               pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k))/(cappa*(peln(i,k+1,j)-peln(i,k,j)))
            enddo
            enddo
         enddo
      endif

      rdg = -rdgas / grav
      do k=ktop,km
         do j=jfirst,jlast
            do i=ifirst,ilast
               pkz(i,j,k) = (rdg*delp(i,j,k)*pt(i,j,k)/delz(i,j,k))**cappa
            enddo
         enddo
      enddo
   endif
! Check dry air mass:
   call drymadj(km, ifirst, ilast,  jfirst,  jlast, ng, mountain, cappa, ptop,   &
                ps, delp, pe, pk, peln, pkz, nq, q, dry_mass, adjust_dry_mass, full_phys )

 end subroutine p_var

! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!        
 subroutine drymadj( km,  ifirst, ilast, jfirst,  jlast,  ng, &  
                     moun, cappa,   ptop, ps, delp, pe,  &
                     pk, peln, pkz, nq, q, dry_mass, adjust_dry_mass, full_phys )

! !INPUT PARAMETERS:
      integer km
      integer ifirst, ilast  ! Long strip
      integer jfirst, jlast  ! Latitude strip    
      integer nq             ! Number of tracers         
      integer ng
      logical moun
      real, intent(in):: dry_mass
      logical, intent(in):: adjust_dry_mass
      logical, intent(in):: full_phys
      real, intent(in):: ptop
      real, intent(in):: cappa

! !INPUT/OUTPUT PARAMETERS:     
      real  delp(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng,km)     !
      real  pe(ifirst-1:ilast+1,km+1,jfirst-1:jlast+1)     !
      real   q(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng,km,nq)
      real  ps(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng)        ! surface pressure
      
      real, intent(inout)::   pk(ifirst:ilast, jfirst:jlast, km+1)
      real, intent(inout):: peln(ifirst:ilast, km+1, jfirst:jlast)    ! Edge pressure
      real, intent(inout)::  pkz(ifirst:ilast, jfirst:jlast, km)
      
! Local
      real  psq(ifirst:ilast,jfirst:jlast)
      real  psd(ifirst:ilast,jfirst:jlast)     ! surface pressure  due to dry air mass
!     parameter ( drym = 98288. )       ! setting for USGS
!     parameter ( drym = 98290. )       ! New setting for USGS

      integer   i, j, k
      real   psmo
      real   psdry, qtot
      real   dpd
      integer ic, ip

      ip = min(3,size(q,4))

!$omp  parallel do          &
!$omp  default(shared)      &
!$omp  private(i,j,k)
      do 1000 j=jfirst,jlast
        do i=ifirst,ilast
           psd(i,j) = ptop
           psq(i,j) = 0.
        enddo

        if( full_phys ) then
          do k=1,km
            do i=ifirst,ilast
              psd(i,j) = psd(i,j) + delp(i,j,k)*(1.-q(i,j,k,1))
              psq(i,j) = psq(i,j) + delp(i,j,k)*sum( q(i,j,k,1:ip) )
            enddo
          enddo
        else
          do k=1,km
            do i=ifirst,ilast
              psd(i,j) = psd(i,j) +  delp(i,j,k)
            enddo
          enddo
        endif
1000  continue

! Check global maximum/minimum
       psmo = g_sum( ps(ifirst:ilast,jfirst:jlast), ifirst, ilast, jfirst, jlast, ng, area, 1) 
      psdry = g_sum(psd, ifirst, ilast, jfirst, jlast, ng, area, 1) 
       qtot = g_sum(psq, ifirst, ilast, jfirst, jlast, ng, area, 1) 

      if(gid==masterproc) then
         write(6,*) 'Total surface pressure (mb) = ', 0.01*psmo
         if ( full_phys ) then
         write(6,*) 'mean dry surface pressure = ', 0.01*psdry
         write(6,*) 'TPW-vapor (kg/m**2) =', (psmo-psdry)/GRAV
         write(6,*) 'TPW-total (kg/m**2) =', qtot/GRAV
         endif
      endif

      if( .not. adjust_dry_mass ) return
      if(moun) then
         dpd = dry_mass - psdry
      else
         dpd = 1000.*100. - psdry
      endif

      if(gid==masterproc) write(6,*) 'dry mass to be added (pascals) =', dpd

!$omp  parallel do            &
!$omp  default(shared)        &
!$omp  private(i,j, ic)

      do 2000 j=jfirst,jlast
         do ic=1,nq
            do i=ifirst,ilast
               q(i,j,km,ic) = q(i,j,km,ic)*delp(i,j,km) / (delp(i,j,km)+dpd)
            enddo
         enddo

! Adjust the lowest Lagrangian layer
         do i=ifirst,ilast
            delp(i,j,km) = delp(i,j,km) + dpd
            pe(i,km+1,j) = pe(i,km,j) + delp(i,j,km)
            ps(i,j) = pe(i,km+1,j)
! adjust pk, peln, pkz
            pk(i,j,km+1) = pe(i,km+1,j) ** cappa
            peln(i,km+1,j) = log(pe(i,km+1,j))
            pkz(i,j,km) = (pk(i,j,km+1)-pk(i,j,km)) /      &
                          (cappa*(peln(i,km+1,j)-peln(i,km,j)))
         enddo
2000  continue

      psmo = g_sum(ps(ifirst:ilast,jfirst:jlast), ifirst, ilast, jfirst, jlast, ng, area, 1) 
      if( gid==masterproc ) write(6,*)                                     &
              'Total (moist) surface pressure after adjustment = ', &
               0.01*psmo

 end subroutine drymadj



 subroutine hydro_eq(km, is, ie, js, je, ps, hs, drym, delp, ak, bk,  &
                     pt, delz, ng, mountain, hybrid_z)
! Input: 
  integer, intent(in):: is, ie, js, je, km, ng
  real, intent(in):: ak(km+1), bk(km+1)
  real, intent(in):: hs(is-ng:ie+ng,js-ng:je+ng)
  real, intent(in):: drym
  logical, intent(in):: mountain
  logical, intent(in):: hybrid_z
! Output
  real, intent(out):: ps(is-ng:ie+ng,js-ng:je+ng)
  real, intent(out)::   pt(is-ng:ie+ng,js-ng:je+ng,km)
  real, intent(out):: delp(is-ng:ie+ng,js-ng:je+ng,km)
  real, intent(inout):: delz(is:ie,js:je,km)
! Local
  real   gz(is:ie,km+1)
  real   ph(is:ie,km+1)
  real mslp, z1, t1, p1, t0, a0, psm
  real ztop, c0
#ifdef INIT_4BYTE
  real*4 dps 
#else
  real dps    ! note that different PEs will get differt dps during initialization
              ! this has no effect after cold start
#endif
  real p0, gztop, ptop
  integer  i,j,k

  if ( gid==masterproc ) write(*,*) 'Initializing ATM hydrostatically'

#ifdef MARS_GCM
      p0 = 893.E2         ! need to tune this value
      t0 = 300.
      u = 0;  v = 0.
      pt = t0
! gztop when zs==0
      gztop = rdgas*t0*log(p0/ak(1))

     do j=js,je
        do i=is,ie
           ps(i,j) = ak(1)*exp((gztop-hs(i,j))/(rdgas*t0))
        enddo
     enddo

      do k=1,km
         do j=js,je
            do i=is,ie
               delp(i,j,k) = ak(k+1)-ak(k) + ps(i,j)*(bk(k+1)-bk(k))
            enddo
         enddo
      enddo
#else
  if ( mountain ) then

! Given p1 and z1 (250mb, 10km)
        mslp = 1009.174*100.
        p1 = 25000.
        z1 = 10.E3 * grav
        t1 = 215.
        t0 = 280.            ! sea-level temp.
        a0 = (t1-t0)/z1
        c0 = t0/a0

     if ( hybrid_z ) then
          ptop = 100.   ! *** hardwired model top *** 
     else
          ptop = ak(1)
     endif

     ztop = z1 + (rdgas*t1)*log(p1/ptop)
     if(gid==masterproc) write(6,*) 'ZTOP is computed as', ztop/grav*1.E-3

     do j=js,je
        do i=is,ie
           ps(i,j) = mslp*( c0/(hs(i,j)+c0))**(1./(a0*rdgas))
        enddo
     enddo

     psm = g_sum(ps(is:ie,js:je), is, ie, js, je, ng, area, 1)
     dps = drym - psm
     if(gid==masterproc) write(6,*) 'Computed mean ps=', psm
     if(gid==masterproc) write(6,*) 'Correction delta-ps=', dps


   do j=js,je
      do i=is,ie
         ps(i,j) = ps(i,j) + dps
         gz(i,   1) = ztop
         gz(i,km+1) = hs(i,j)
         ph(i,   1) = ptop                                                     
         ph(i,km+1) = ps(i,j)                                               
      enddo

      if ( hybrid_z ) then
!---------------
! Hybrid Z
!---------------
        do k=km,2,-1
           do i=is,ie
              gz(i,k) = gz(i,k+1) - delz(i,j,k)*grav 
           enddo
        enddo
! Correct delz at the top:
        do i=is,ie
            delz(i,j,1) = (gz(i,2) - ztop) / grav
        enddo
 
        do k=2,km
           do i=is,ie
              if ( gz(i,k) >= z1 ) then
! Isothermal
                 ph(i,k) = ptop*exp( (gz(i,1)-gz(i,k))/(rdgas*t1) )
              else
! Constant lapse rate region (troposphere)
                 ph(i,k) = ps(i,j)*((hs(i,j)+c0)/(gz(i,k)+c0))**(1./(a0*rdgas))
              endif
           enddo
        enddo
      else
!---------------
! Hybrid sigma-p
!---------------
        do k=2,km+1
           do i=is,ie
              ph(i,k) = ak(k) + bk(k)*ps(i,j)
           enddo
        enddo

        do k=2,km
           do i=is,ie
              if ( ph(i,k) <= p1 ) then
! Isothermal
                 gz(i,k) = ztop + (rdgas*t1)*log(ptop/ph(i,k))
              else
! Constant lapse rate region (troposphere)
                 gz(i,k) = (hs(i,j)+c0)/(ph(i,k)/ps(i,j))**(a0*rdgas) - c0
              endif
           enddo
        enddo
      endif  ! end hybrid_z

! Convert geopotential to Temperature
      do k=1,km
         do i=is,ie
              pt(i,j,k) = (gz(i,k)-gz(i,k+1))/(rdgas*(log(ph(i,k+1)/ph(i,k))))
              pt(i,j,k) = max(t1, pt(i,j,k))
            delp(i,j,k) = ph(i,k+1) - ph(i,k)
         enddo
      enddo
   enddo    ! j-loop

   if ( hybrid_z ) then 
!      call prt_maxmin('INIT_hydro: delz', delz, is, ie, js, je,  0, km, 1., gid==masterproc)
!      call prt_maxmin('INIT_hydro: DELP', delp, is, ie, js, je, ng, km, 1., gid==masterproc)
   endif
!  call prt_maxmin('INIT_hydro: PT  ', pt,   is, ie, js, je, ng, km, 1., gid==masterproc)

  else    ! no mountain
!$omp parallel do private (i, j)
        do j=js,je
           do i=is,ie
              ps(i,j) = 1.E5
           enddo
           do k=1,km
               do i=is,ie
                  pt(i,j,k) = 273.0
                  delp(i,j,k) = ak(k+1)-ak(k) + ps(i,j)*(bk(k+1)-bk(k))
               enddo
           enddo
        enddo
  endif
#endif

 end subroutine hydro_eq

! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

end module init_hydro
