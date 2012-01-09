! $Id: init_hydro.F90,v 19.0 2012/01/06 19:59:21 fms Exp $

module init_hydro_mod

      use constants_mod, only: grav, rdgas, rvgas
      use fv_grid_utils_mod,    only: g_sum
      use fv_grid_tools_mod,    only: area
      use fv_mp_mod,        only: gid, masterproc
!     use fv_diagnostics_mod, only: prt_maxmin

      implicit none
      private

      public :: p_var, hydro_eq

!---- version number -----
      character(len=128) :: version = '$Id: init_hydro.F90,v 19.0 2012/01/06 19:59:21 fms Exp $'
      character(len=128) :: tagname = '$Name: siena $'

contains

!-------------------------------------------------------------------------------
 subroutine p_var(km, ifirst, ilast, jfirst, jlast, ptop, ptop_min,    &
                  delp, delz, pt, ps,  pe, peln, pk, pkz, cappa, q, ng, nq,    &
                  dry_mass, adjust_dry_mass, mountain, moist_phys,      &
                  hydrostatic, ktop, nwat, make_nh)
               
! Given (ptop, delp) computes (ps, pk, pe, peln, pkz)
! Input:
   integer,  intent(in):: km
   integer,  intent(in):: ifirst, ilast            ! Longitude strip
   integer,  intent(in):: jfirst, jlast            ! Latitude strip
   integer,  intent(in):: nq, nwat
   integer,  intent(in):: ng
   integer,  intent(in):: ktop
   logical, intent(in):: adjust_dry_mass, mountain, moist_phys, hydrostatic
   real, intent(in):: dry_mass, cappa, ptop, ptop_min
   real, intent(in   )::   pt(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng, km)
   real, intent(inout):: delz(ifirst   :ilast   ,jfirst   :jlast   , km)
   real, intent(inout):: delp(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng, km)
   real, intent(inout)::    q(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng, km, nq)
   logical, optional:: make_nh
! Output:
   real, intent(out) ::   ps(ifirst-ng:ilast+ng, jfirst-ng:jlast+ng)
   real, intent(out) ::   pk(ifirst:ilast, jfirst:jlast, km+1)
   real, intent(out) ::   pe(ifirst-1:ilast+1,km+1,jfirst-1:jlast+1) ! Ghosted Edge pressure
   real, intent(out) :: peln(ifirst:ilast, km+1, jfirst:jlast)    ! Edge pressure
   real, intent(out) ::  pkz(ifirst:ilast, jfirst:jlast, km)

! Local
   real ratio(ifirst:ilast)
   real pek, lnp, ak1, rdg, dpd, zvir
   integer i, j, k


! Check dry air mass & compute the adjustment amount:
   call drymadj(km, ifirst, ilast,  jfirst,  jlast, ng, cappa, ptop, ps, &
                delp, q, nq, nwat, dry_mass, adjust_dry_mass, moist_phys, dpd)

   pek = ptop ** cappa

   do j=jfirst,jlast
      do i=ifirst,ilast
         pe(i,1,j) = ptop
         pk(i,j,1) = pek
      enddo

      if ( adjust_dry_mass ) then
         do i=ifirst,ilast
            ratio(i) = 1. + dpd/(ps(i,j)-ptop)
         enddo 
         do k=1,km
            do i=ifirst,ilast
               delp(i,j,k) = delp(i,j,k) * ratio(i)
            enddo
         enddo
      endif

      do k=2,km+1
         do i=ifirst,ilast
            pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
            peln(i,k,j) = log(pe(i,k,j))
            pk(i,j,k) = exp( cappa*peln(i,k,j) )
!            pk(i,j,k) = pe(i,k,j)**cappa
         enddo
      enddo

      do i=ifirst,ilast
         ps(i,j) = pe(i,km+1,j)
      enddo

      if( ptop < ptop_min ) then
!---- small ptop modification -------------
          ak1 = (cappa + 1.) / cappa
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
      if ( present(make_nh) ) then
          if ( make_nh ) then
             do k=1,km
                do j=jfirst,jlast
                   do i=ifirst,ilast
                      delz(i,j,k) = rdg*pt(i,j,k)*(peln(i,k+1,j)-peln(i,k,j))
                   enddo
                enddo
             enddo
             if(gid==0) write(*,*) 'delz computed from hydrostatic state'
          endif
      endif

     if ( moist_phys ) then
!------------------------------------------------------------------
! The following form is the same as in "fv_update_phys.F90"
! Therefore, restart reproducibility is only enforced in diabatic cases
!------------------------------------------------------------------
! For adiabatic runs, use "-DSOLO_REPRO" to enforce reproducibility
       zvir = rvgas/rdgas - 1.
       do k=ktop,km
          do j=jfirst,jlast
             do i=ifirst,ilast
                pkz(i,j,k) = exp( cappa*log(rdg*delp(i,j,k)*pt(i,j,k)*    &
                                (1.+zvir*q(i,j,k,1))/delz(i,j,k)) )
             enddo
          enddo
       enddo
     else
       do k=ktop,km
          do j=jfirst,jlast
             do i=ifirst,ilast
                pkz(i,j,k) = exp( cappa*log(rdg*delp(i,j,k)*pt(i,j,k)/delz(i,j,k)) )
             enddo
          enddo
       enddo
     endif

   endif

 end subroutine p_var



 subroutine drymadj(km,  ifirst, ilast, jfirst,  jlast,  ng, &  
                    cappa,   ptop, ps, delp, q,  nq,  nwat,  &
                    dry_mass, adjust_dry_mass, moist_phys, dpd)

! !INPUT PARAMETERS:
      integer km
      integer ifirst, ilast  ! Long strip
      integer jfirst, jlast  ! Latitude strip    
      integer nq, ng, nwat
      real, intent(in):: dry_mass
      real, intent(in):: ptop
      real, intent(in):: cappa
      logical, intent(in):: adjust_dry_mass
      logical, intent(in):: moist_phys

! !INPUT/OUTPUT PARAMETERS:     
      real, intent(in)::   q(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng,km,nq)
      real, intent(in)::delp(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng,km)     !
      real, intent(inout):: ps(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng)        ! surface pressure
      real, intent(out):: dpd
! Local
      real  psd(ifirst:ilast,jfirst:jlast)     ! surface pressure  due to dry air mass
      real  psmo, psdry
      integer i, j, k

!$omp parallel do default(shared) 
      do j=jfirst,jlast

         do i=ifirst,ilast
             ps(i,j) = ptop
            psd(i,j) = ptop
         enddo

         do k=1,km
            do i=ifirst,ilast
               ps(i,j) = ps(i,j) + delp(i,j,k)
            enddo
         enddo

       if ( nwat>=1 ) then
          do k=1,km
             do i=ifirst,ilast
                psd(i,j) = psd(i,j) + delp(i,j,k)*(1. - sum(q(i,j,k,1:nwat)))
             enddo
          enddo
        else
          do i=ifirst,ilast
             psd(i,j) = ps(i,j)
          enddo
        endif
      enddo

! Check global maximum/minimum
#ifndef QUICK_SUM
      psdry = g_sum(psd, ifirst, ilast, jfirst, jlast, ng, area, 1, .true.) 
       psmo = g_sum(ps(ifirst:ilast,jfirst:jlast), ifirst, ilast, jfirst, jlast,  &
                     ng, area, 1, .true.) 
#else
      psdry = g_sum(psd, ifirst, ilast, jfirst, jlast, ng, area, 1) 
       psmo = g_sum(ps(ifirst:ilast,jfirst:jlast), ifirst, ilast, jfirst, jlast,  &
                     ng, area, 1) 
#endif

      if(gid==masterproc) then
         write(6,*) 'Total surface pressure (mb) = ', 0.01*psmo
         if ( moist_phys ) then
              write(6,*) 'mean dry surface pressure = ', 0.01*psdry
              write(6,*) 'Total Water (kg/m**2) =', real(psmo-psdry,4)/GRAV
         endif
      endif

      if( adjust_dry_mass ) Then
          dpd = real(dry_mass - psdry,4)
          if(gid==masterproc) write(6,*) 'dry mass to be added (pascals) =', dpd
      endif

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
  real(kind=4) ::  dps 
#else
  real dps    ! note that different PEs will get differt dps during initialization
              ! this has no effect after cold start
#endif
  real p0, gztop, ptop
  integer  i,j,k

  if ( gid==masterproc ) write(*,*) 'Initializing ATM hydrostatically'

#if defined(MARS_GCM)
  if ( gid==masterproc ) write(*,*) 'Initializing Mars'
      p0 = 6.5E2         !
      t0 = 200.0

!         Isothermal temperature
      pt = t0

      gztop = rdgas*t0*log(p0/ak(1))        ! gztop when zs==0

     do j=js,je
        do i=is,ie
           ps(i,j) = ak(1)*exp((gztop-hs(i,j))/(rdgas*t0))
        enddo
     enddo

     psm = g_sum(ps(is:ie,js:je), is, ie, js, je, ng, area, 1, .true.)
     dps = drym - psm

     if(gid==masterproc) write(6,*) 'Initializing:  Computed mean ps=', psm
     if(gid==masterproc) write(6,*) '            Correction delta-ps=', dps

!           Add correction to surface pressure to yield desired
!                globally-integrated atmospheric mass  (drym)
     do j=js,je
        do i=is,ie
           ps(i,j) = ps(i,j) + dps
        enddo
     enddo

      do k=1,km
         do j=js,je
            do i=is,ie
               delp(i,j,k) = ak(k+1)-ak(k) + ps(i,j)*(bk(k+1)-bk(k))
            enddo
         enddo
      enddo

#elif defined(VENUS_GCM)
  if ( gid==masterproc ) write(*,*) 'Initializing Venus'
      p0 = 92.E5         ! need to tune this value
      t0 = 700.
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
  if ( gid==masterproc ) write(*,*) 'Initializing Earth'
! Given p1 and z1 (250mb, 10km)
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

  if ( mountain ) then
     mslp = 100917.4
     do j=js,je
        do i=is,ie
           ps(i,j) = mslp*( c0/(hs(i,j)+c0))**(1./(a0*rdgas))
        enddo
     enddo
     psm = g_sum(ps(is:ie,js:je), is, ie, js, je, ng, area, 1, .true.)
     dps = drym - psm
     if(gid==masterproc) write(6,*) 'Computed mean ps=', psm
     if(gid==masterproc) write(6,*) 'Correction delta-ps=', dps
  else
     mslp = 1000.E2
     do j=js,je
        do i=is,ie
           ps(i,j) = mslp
        enddo
     enddo
     dps = 0.
  endif


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

#endif

 end subroutine hydro_eq


end module init_hydro_mod
