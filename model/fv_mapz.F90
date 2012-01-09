module fv_mapz_mod

  use constants_mod,     only: radius, pi, rvgas, rdgas, grav
  use fv_grid_tools_mod, only: area, dx, dy, rdxa, rdya
  use fv_grid_utils_mod, only: g_sum, ptop, ptop_min, cosa_s, rsin2
  use fv_fill_mod,       only: fillz
  use fv_eta_mod ,       only: compute_dz_L32
  use fv_eta_mod ,       only: compute_dz_L32
  use fv_mp_mod,         only: gid, domain, square_domain
  use mpp_domains_mod,   only: mpp_update_domains
  use mpp_mod,           only: FATAL, mpp_error, get_unit, stdlog, mpp_root_pe, mpp_pe

  implicit none
  real, parameter::  r3 = 1./3., r23 = 2./3., r12 = 1./12.
  real(kind=4) :: E_Flux
  private

  public compute_total_energy, Lagrangian_to_Eulerian,    &
         rst_remap, mappm, E_Flux

!---- version number -----
  character(len=128) :: version = '$Id: fv_mapz.F90,v 19.0 2012/01/06 19:57:44 fms Exp $'
  character(len=128) :: tagname = '$Name: siena $'

contains

 subroutine Lagrangian_to_Eulerian(do_consv, consv, ps, pe, delp, pkz, pk,   &
                      pdt, km, is,ie,js,je, isd,ied,jsd,jed,       &
                      nq, sphum, u, v, w, delz, pt, q, hs, r_vir, cp,  &
                      akap, kord_mt, kord_wz, kord_tr, kord_tm,  peln, te0_2d,        &
                      ng, ua, va, omga, te, fill, reproduce_sum,        &
                      ak, bk, ks, ze0, remap_t, hydrostatic, hybrid_z, do_omega, ktop)
  logical, intent(in):: do_consv
  real,    intent(in):: pdt                   ! phys time step
  integer, intent(in):: km
  integer, intent(in):: nq                    ! number of tracers (including h2o)
  integer, intent(in):: sphum                 ! index for water vapor (specific humidity)
  integer, intent(in):: ng
  integer, intent(in):: is,ie,isd,ied         ! starting & ending X-Dir index
  integer, intent(in):: js,je,jsd,jed         ! starting & ending Y-Dir index
  integer, intent(in):: ks, ktop
  integer, intent(in):: kord_mt               ! Mapping order for the vector winds
  integer, intent(in):: kord_wz               ! Mapping order/option for w
  integer, intent(in):: kord_tr(nq)           ! Mapping order for tracers
  integer, intent(in):: kord_tm               ! Mapping order for thermodynamics

  real, intent(in):: consv                 ! factor for TE conservation
  real, intent(in):: r_vir
  real, intent(in):: cp
  real, intent(in):: akap
  real, intent(in):: hs(isd:ied,jsd:jed)  ! surface geopotential
  real, intent(in):: te0_2d(is:ie,js:je)

  logical, intent(in):: fill                  ! fill negative tracers
  logical, intent(in):: reproduce_sum
  logical, intent(in):: do_omega
  real, intent(in) :: ak(km+1)
  real, intent(in) :: bk(km+1)

! !INPUT/OUTPUT
  real, intent(inout):: pk(is:ie,js:je,km+1) ! pe to the kappa
  real, intent(inout):: q(isd:ied,jsd:jed,km,*)
  real, intent(inout):: delp(isd:ied,jsd:jed,km) ! pressure thickness
  real, intent(inout)::  pe(is-1:ie+1,km+1,js-1:je+1) ! pressure at layer edges
  real, intent(inout):: ps(isd:ied,jsd:jed)      ! surface pressure
  real, intent(inout):: ze0(is:ie,js:je,km+1)    ! Specified height at edges (m)

! u-wind will be ghosted one latitude to the north upon exit
  real, intent(inout)::  u(isd:ied  ,jsd:jed+1,km)   ! u-wind (m/s)
  real, intent(inout)::  v(isd:ied+1,jsd:jed  ,km)   ! v-wind (m/s)
  real, intent(inout)::  w(isd:ied  ,jsd:jed  ,km)   ! vertical velocity (m/s)
  real, intent(inout):: pt(isd:ied  ,jsd:jed  ,km)   ! cp*virtual potential temperature 
                                                     ! as input; output: temperature
  real, intent(inout):: delz(is:ie,js:je,km)   ! delta-height (m)
  logical, intent(in):: remap_t
  logical, intent(in):: hydrostatic
  logical, intent(in):: hybrid_z

  real, intent(inout)::   ua(isd:ied,jsd:jed,km)   ! u-wind (m/s) on physics grid
  real, intent(inout)::   va(isd:ied,jsd:jed,km)   ! v-wind (m/s) on physics grid
  real, intent(inout):: omga(isd:ied,jsd:jed,km)   ! vertical press. velocity (pascal/sec)
  real, intent(inout)::   peln(is:ie,km+1,js:je)     ! log(pe)
  real, intent(out)::    pkz(is:ie,js:je,km)       ! layer-mean pk for converting t to pt
  real, intent(out)::     te(is:ie,js:je,km)

! !DESCRIPTION:
!
! !REVISION HISTORY:
! SJL 03.11.04: Initial version for partial remapping
!
!-----------------------------------------------------------------------
  integer :: i,j,k 
     real q_source(is:ie,js:je,nq)    ! numerical tracer source from surface
                                      ! in case fillz is not sufficient
      real te_2d(is:ie,js:je)
      real zsum0(is:ie,js:je)
      real zsum1(is:ie,js:je)
      real   q2(is:ie,km)
      real  dp2(is:ie,km)
      real  pe1(is:ie,km+1)
      real  pe2(is:ie,km+1)
      real  pk1(is:ie,km+1)
      real  pk2(is:ie,km+1)
      real  pn2(is:ie,km+1)
      real  pe0(is:ie+1,km+1)
      real  pe3(is:ie+1,km+1)
      real phis(is:ie,km+1)
      real   gz(is:ie)
! for nonhydrostatic option with hybrid_z coordinate
      real ze1(is:ie,km+1), ze2(is:ie,km+1), deng(is:ie,km)
      real dz1(km), ztop, z_rat
      real rcp, rg, ak1, tmp, tpe, cv, rgama, rrg
      real bkh, dtmp, dlnp
      integer iq, n, kp, k_next
      logical te_map
      real k1k, kapag

        k1k = akap / (1.-akap)    ! rg/Cv=0.4
      kapag = -akap / grav
         rg = akap * cp
         cv = cp - rg
      rgama = 1. - akap           ! cv/cp
        rcp = 1./ cp
        ak1 = (akap + 1.) / akap
        rrg = - rg/grav

      if ( kord_tm < 0 ) then
           te_map = .false.
          if ( remap_t ) then
! Note: pt at this stage is cp*Theta_v
! Transform virtual pt to virtual Temp
             if ( hydrostatic ) then
!$omp parallel do default(shared)
             do k=1,km
                do j=js,je
                   do i=is,ie
                      pt(i,j,k) = pt(i,j,k) * (pk(i,j,k+1)-pk(i,j,k)) /  &
                                 (rg*(peln(i,k+1,j)-peln(i,k,j)) )
                   enddo
                enddo
             enddo
             else
               if ( ktop>1 ) then
!$omp parallel do default(shared)
                 do k=1,ktop-1
                    do j=js,je
                    do i=is,ie
                       pt(i,j,k) = pt(i,j,k) * (pk(i,j,k+1)-pk(i,j,k)) /  &
                                  (rg*(peln(i,k+1,j)-peln(i,k,j)) )
                    enddo
                    enddo
                 enddo
               endif
!$omp parallel do default(shared)
               do k=ktop,km
                  do j=js,je
                  do i=is,ie
                     pt(i,j,k) = rcp*pt(i,j,k)*exp(k1k*log(kapag*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
                  enddo
                  enddo
               enddo
             endif         ! hydro test
          endif            ! remap_t test
      else
           te_map = .true.
           call pkez(km, is, ie, js, je, pe, pk, akap, peln, pkz)
!           call cubed_to_latlon(u, v, ua, va, dx, dy, rdxa, rdya, km, 1)
! Compute cp*T + KE
!$omp parallel do default(shared) 
           do k=1,km
              do j=js,je
                 do i=is,ie
                    te(i,j,k) = 0.25*rsin2(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                                                 v(i,j,k)**2+v(i+1,j,k)**2 -  &
                               (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s(i,j))  &
                              +  pt(i,j,k)*pkz(i,j,k)
                 enddo
              enddo
           enddo
     endif

     if ( .not.hydrostatic ) then
        if ( hybrid_z ) then
           if ( km==32 ) then
                call compute_dz_L32(km, ztop, dz1)
           else
                call mpp_error(FATAL,'==> Error from fv_mapz: hybrid_z works only with L32')
           endif
        else
!$omp parallel do default(shared)
           do k=1,km
              do j=js,je
                 do i=is,ie
                    delz(i,j,k) = -delz(i,j,k) / delp(i,j,k) ! ="specific volume"/grav
                 enddo
              enddo
           enddo
        endif
     endif

!$omp parallel do default(shared) private(z_rat, kp, k_next, bkh, deng, dp2, pe0, pe1, pe2, pe3, pk1, pk2, pn2, phis, q2, ze1, ze2)
  do 1000 j=js,je+1

        do k=1,km+1
           do i=is,ie
              pe1(i,k) = pe(i,k,j)
           enddo
        enddo

        do i=is,ie
           pe2(i,   1) = ptop
           pe2(i,km+1) = pe(i,km+1,j)
        enddo

  if ( j < (je+1) )  then 
! update ps
        do i=is,ie
            ps(i,j) = pe1(i,km+1)
        enddo

   if ( hybrid_z ) then
!--------------------------
! hybrid z_p coordinate
!--------------------------

        do i=is,ie
           ze1(i,km+1) = ze0(i,j,km+1)
        enddo

        do k=km,1,-1
           do i=is,ie
              ze1(i,k) = ze1(i,k+1) - delz(i,j,k)   ! current height
           enddo
        enddo
!
! Copy ztop; the top layer must be thick enough to prevent numerical problems.
!
        do i=is,ie
           ze2(i,  1) = ze1(i,1)
           ze0(i,j,1) = ze1(i,1)      ! Note: ze0 (top) updated
        enddo

        do k=2,km+1
           do i=is,ie
              ze2(i,k) = ze0(i,j,k)   ! specified height
           enddo
        enddo

! Check/fix monotonicity of top 3 layers thickness
!-----------------------------------------------
        do i=is,ie
           do k=1,3
!          if ( (ze2(i,k)-ze2(i,k+1)) < (ze2(i,k+1)-ze2(i,k+2)) ) then
!                ze2(i,k+1) = 0.5*(ze2(i,k) + ze2(i,k+2))
!          endif
           z_rat = (ze2(i,1)-ze2(i,4)) / (dz1(1)+dz1(2)+dz1(3))
           ze2(i,3) = ze2(i,4) + dz1(3)*z_rat
           ze2(i,2) = ze2(i,3) + dz1(2)*z_rat
           enddo
        enddo
!-----------------------------------------------

        do k=1,km
           do i=is,ie
              deng(i,k) = -delp(i,j,k)/delz(i,j,k)  ! density * grav
           enddo
        enddo

        call remap_z(km, ze1, deng, km, ze2, deng, is, ie, 1, abs(kord_tm))
!-------------
! Update delz
!-------------
        do k=1,km
           do i=is,ie
              delz(i,j,k) = ze2(i,k+1) - ze2(i,k)
           enddo
        enddo

!------------
! update delp
!------------
        do k=1,km-1
           do i=is,ie
               dp2(i,k  ) = -delz(i,j,k) * deng(i,k)
               pe2(i,k+1) =     pe2(i,k) +  dp2(i,k)
           enddo
        enddo

        do i=is,ie
           dp2(i,km) = pe2(i,km+1) - pe2(i,km)  ! to reduce rounding error
        enddo
   else
!
! Hybrid sigma-P coordinate:
!
        do k=2,ks+1
           do i=is,ie
              pe2(i,k) = ak(k)
           enddo
        enddo
        do k=ks+2,km
           do i=is,ie
              pe2(i,k) = ak(k) + bk(k)*pe(i,km+1,j)
           enddo
        enddo
        do k=1,km
           do i=is,ie
              dp2(i,k) = pe2(i,k+1) - pe2(i,k)
           enddo
        enddo
   endif

!------------
! update delp
!------------
      do k=1,km
         do i=is,ie
            delp(i,j,k) = dp2(i,k)
         enddo
      enddo

!----------------
! Map constituents
!----------------
       if( nq /= 0 ) then
!------------------------------------------------------------------
! Do remapping one tracer at a time; seems to be faster on the SGI
! It requires less memory than mapn_ppm
!------------------------------------------------------------------
          do iq=1,nq
             call map1_q2(km, pe1, q(isd,jsd,1,iq),     &
                          km, pe2, q2, dp2,             &
                          is, ie, 0, kord_tr(iq), j, isd, ied, jsd, jed)
!           if (fill) call fillz(ie-is+1, km, 1, q2, dp2, q_source(is,j,iq))
            if (fill) call fillz(ie-is+1, km, 1, q2, dp2)
            do k=1,km
               do i=is,ie
                  q(i,j,k,iq) = q2(i,k)
               enddo
            enddo
          enddo
       endif

!------------------
! Compute p**cappa
!------------------
   do k=1,km+1
      do i=is,ie
         pk1(i,k) = pk(i,j,k)
      enddo
   enddo

   do i=is,ie
      pn2(i,   1) = peln(i,   1,j)
      pn2(i,km+1) = peln(i,km+1,j)
      pk2(i,   1) = pk1(i,   1)
      pk2(i,km+1) = pk1(i,km+1)
   enddo

   do k=2,km
      do i=is,ie
!        pk2(i,k) = pe2(i,k) ** akap
         pn2(i,k) = log(pe2(i,k))
         pk2(i,k) = exp(akap*pn2(i,k))
      enddo
   enddo

   if ( te_map ) then
!---------------------
! Compute Total Energy
!---------------------
        do i=is,ie
           phis(i,km+1) = hs(i,j)
        enddo
        do k=km,1,-1
           do i=is,ie
              phis(i,k) = phis(i,k+1) + pt(i,j,k)*(pk1(i,k+1)-pk1(i,k))
           enddo
        enddo
        do k=1,km+1
           do i=is,ie
              phis(i,k) = phis(i,k) * pe1(i,k)
           enddo
        enddo
        do k=1,km
           do i=is,ie
              te(i,j,k) = te(i,j,k)+(phis(i,k+1)-phis(i,k))/(pe1(i,k+1)-pe1(i,k))
           enddo
        enddo
!----------------
! Map Total Energy
!----------------
        call map1_ppm (km,   pe1,  te,       &
                       km,   pe2,  te,       &
                       is, ie, j, is, ie, js, je, 1, kord_tm)
   else
     if ( remap_t ) then
!----------------------------------
! Map t using ze1 (hybrid_z) or logp 
!----------------------------------
       if ( hybrid_z ) then
           do k=1,km
              do i=is,ie
                 deng(i,k) = pt(i,j,k)
              enddo
           enddo
           call remap_z(km, ze1, deng, km, ze2, deng, is, ie, 2, abs(kord_tm))
           do k=1,km
              do i=is,ie
                 pt(i,j,k) = deng(i,k)
              enddo
           enddo
       else
         call map1_ppm (km,  peln(is,1,j),  pt,    &
                        km,  pn2,           pt,    &
                        is, ie, j, isd, ied, jsd, jed, 2, abs(kord_tm))
       endif
     else
!----------------
! Map pt using pk
!----------------
       call map1_ppm (km,  pk1,  pt,           &
                      km,  pk2,  pt,           &
                      is, ie, j, isd, ied, jsd, jed, 1, abs(kord_tm))
     endif
   endif

   if ( .not. hydrostatic ) then
! Remap vertical wind:
        call map1_ppm (km,   pe1,  w,       &
                       km,   pe2,  w,       &
                       is, ie, j, isd, ied, jsd, jed, -2, kord_wz)
     if ( .not. hybrid_z ) then
! Remap delz for hybrid sigma-p coordinate
        call map1_ppm (km,   pe1, delz,    &
                       km,   pe2, delz,    &
                       is, ie, j, is,  ie,  js,  je,  1, abs(kord_tm))
        do k=1,km
           do i=is,ie
              delz(i,j,k) = -delz(i,j,k)*dp2(i,k)
           enddo
        enddo
     endif
   endif

!----------
! Update pk
!----------
   do k=2,km
      do i=is,ie
         pk(i,j,k) = pk2(i,k)
      enddo
   enddo

!----------------
   if ( do_omega ) then
! Start do_omega
! Copy omega field to pe3
      do i=is,ie
         pe3(i,1) = 0.
      enddo
      do k=2,km+1
         do i=is,ie
            pe3(i,k) = omga(i,j,k-1)
         enddo
      enddo
   endif

   do k=1,km+1
      do i=is,ie
          pe0(i,k)   = peln(i,k,j)
         peln(i,k,j) =  pn2(i,k)
      enddo
   enddo

!------------
! Compute pkz
!------------
   if ( hydrostatic ) then
      do k=1,km
         do i=is,ie
            pkz(i,j,k) = (pk2(i,k+1)-pk2(i,k))/(akap*(peln(i,k+1,j)-peln(i,k,j)))
         enddo
      enddo
   else
      if ( ktop>1 ) then
         do k=1,ktop-1
         do i=is,ie
            pkz(i,j,k) = (pk2(i,k+1)-pk2(i,k))/(akap*(peln(i,k+1,j)-peln(i,k,j)))
         enddo
         enddo
      endif
      if ( remap_t ) then
! Note: pt at this stage is T_v
         do k=ktop,km
         do i=is,ie
            pkz(i,j,k) = exp(akap*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
         enddo
         enddo
      else
! Note: pt at this stage is cp*Theta_v
         do k=ktop,km
         do i=is,ie
            pkz(i,j,k) = exp( k1k*log(kapag*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)) )
         enddo
         enddo
      endif
   endif

! Interpolate omega/pe3 (defined at pe0) to remapped cell center (dp2)
   if ( do_omega ) then
   do k=1,km
      do i=is,ie
         dp2(i,k) = 0.5*(peln(i,k,j) + peln(i,k+1,j))
      enddo
   enddo
   do i=is,ie
       k_next = 1
       do n=1,km
          kp = k_next
          do k=kp,km
             if( dp2(i,n) <= pe0(i,k+1) .and. dp2(i,n) >= pe0(i,k) ) then
                 omga(i,j,n) = pe3(i,k)  +  (pe3(i,k+1) - pe3(i,k)) *    &
                       (dp2(i,n)-pe0(i,k)) / (pe0(i,k+1)-pe0(i,k) )
                 k_next = k
                 exit
             endif
          enddo
       enddo
   enddo
   endif     ! end do_omega

  endif !(j < je+1)

 if ( .not.hybrid_z ) then
      do i=is,ie+1
         pe0(i,1) = pe(i,1,j)
      enddo
!------
! map u
!------
      do k=2,km+1
         do i=is,ie
            pe0(i,k) = 0.5*(pe(i,k,j-1)+pe1(i,k))
         enddo
      enddo


      do k=1,ks+1
         do i=is,ie+1
            pe3(i,k) = ak(k)
         enddo
      enddo

      do k=ks+2,km+1
         bkh = 0.5*bk(k)
         do i=is,ie
            pe3(i,k) = ak(k) + bkh*(pe(i,km+1,j-1)+pe1(i,km+1))
         enddo
      enddo

      call map1_ppm( km, pe0(is:ie,:),   u,       &
                     km, pe3(is:ie,:),   u,       &
                     is, ie, j, isd, ied, jsd, jed+1, -1, kord_mt)

   if (j < je+1) then
!------
! map v
!------
       do k=2,km+1
          do i=is,ie+1
             pe0(i ,k) = 0.5*(pe(i-1,k,j)+pe(i,k,j))
          enddo
       enddo
       do k=ks+2,km+1
          bkh = 0.5*bk(k)
          do i=is,ie+1
             pe3(i,k) = ak(k) + bkh*(pe(i-1,km+1,j)+pe(i,km+1,j))
          enddo
       enddo

       call map1_ppm (km, pe0,  v,              &
                      km, pe3,  v, is, ie+1,    &
                      j, isd, ied+1, jsd, jed, -1, kord_mt)
   endif ! (j < je+1)
 endif    ! end hybrid_z check

     do k=1,km
        do i=is,ie
           ua(i,j,k) = pe2(i,k+1)
        enddo
     enddo

1000  continue

if ( hybrid_z ) then   
!------- Hybrid_z section ---------------
   if ( square_domain ) then
     call mpp_update_domains(ua , domain,  whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
   else
     call mpp_update_domains(ua , domain, complete=.true.)
   endif

! u-wind

!$omp parallel do default(shared) private(pe1, pe2)
   do j=js,je+1
      do i=is,ie
         pe1(i,1) = ptop
         pe2(i,1) = ptop
      enddo
      do k=2,km+1
         do i=is,ie
            pe1(i,k) = 0.5*(pe(i,k,  j-1) + pe(i,k,j  ))
            pe2(i,k) = 0.5*(ua(i,j-1,k-1) + ua(i,j,k-1))
         enddo
      enddo

      call map1_ppm( km, pe1,   u,       &
                     km, pe2,   u,       &
                     is, ie, j, isd, ied, jsd, jed+1, -1, kord_mt)
   enddo

! v-wind
!$omp parallel do default(shared) private(pe0, pe3)
   do j=js,je
      do i=is,ie+1
         pe0(i,1) = ptop
         pe3(i,1) = ptop
      enddo

      do k=2,km+1
         do i=is,ie+1
            pe0(i,k) = 0.5*(pe(i-1,k,j  ) + pe(i,k,j  ))
            pe3(i,k) = 0.5*(ua(i-1,j,k-1) + ua(i,j,k-1))
         enddo
      enddo

      call map1_ppm (km, pe0,  v,              &
                     km, pe3,  v, is, ie+1,    &
                     j, isd, ied+1, jsd, jed, -1, kord_mt)
   enddo
endif 
!------------- Hybrid_z section ----------------------

!$omp parallel do default(shared)
  do k=2,km
     do j=js,je
        do i=is,ie
           pe(i,k,j) = ua(i,j,k-1)
        enddo
     enddo
  enddo

!  call cubed_to_latlon(u,  v, ua, va, dx, dy, rdxa, rdya, km, 1)

  if( do_consv .and. consv > 0. ) then

    if ( te_map ) then
!$omp parallel do default(shared) 
      do j=js,je
          do i=is,ie
             te_2d(i,j) = te(i,j,1)*delp(i,j,1)
          enddo
          do k=2,km
             do i=is,ie
                te_2d(i,j) = te_2d(i,j) + te(i,j,k)*delp(i,j,k)
             enddo
          enddo
      enddo
    else
!$omp parallel do default(shared) private(gz, phis)
      do j=js,je
        if ( remap_t ) then
         if ( hydrostatic ) then
         do i=is,ie
            gz(i) = hs(i,j)
            do k=1,km
               gz(i) = gz(i) + rg*pt(i,j,k)*(peln(i,k+1,j)-peln(i,k,j))
            enddo
         enddo
         do i=is,ie
            te_2d(i,j) = pe(i,km+1,j)*hs(i,j) - pe(i,1,j)*gz(i)
         enddo

         do k=1,km
            do i=is,ie
               te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(cp*pt(i,j,k) +   &
                            0.25*rsin2(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                                             v(i,j,k)**2+v(i+1,j,k)**2 -  &
                           (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s(i,j)))
            enddo
         enddo
         else
! non-hydrostatic & remap_t
           do i=is,ie
              phis(i,km+1) = hs(i,j)
              do k=km,1,-1
                 phis(i,k) = phis(i,k+1) - grav*delz(i,j,k)
              enddo
           enddo
           do i=is,ie
              te_2d(i,j) = 0.
           enddo
           do k=1,km
              do i=is,ie
! KE using 3D winds:
                 te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*( cv*pt(i,j,k) +  &
                              0.5*(phis(i,k)+phis(i,k+1) + w(i,j,k)**2 + 0.5*rsin2(i,j)*( &
                              u(i,j,k)**2+u(i,j+1,k)**2 + v(i,j,k)**2+v(i+1,j,k)**2 -  &
                             (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s(i,j))))
              enddo
           enddo
         endif
        else
         if ( hydrostatic ) then
            do i=is,ie
               gz(i) = hs(i,j)
               do k=1,km
                  gz(i) = gz(i) + pt(i,j,k)*(pk(i,j,k+1)-pk(i,j,k))
               enddo
            enddo

            do i=is,ie
               te_2d(i,j) = pe(i,km+1,j)*hs(i,j) - pe(i,1,j)*gz(i)
            enddo
            do k=1,km
               do i=is,ie
                  te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(pt(i,j,k)*pkz(i,j,k) +   &
                               0.25*rsin2(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                                                v(i,j,k)**2+v(i+1,j,k)**2 -  &
                            (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s(i,j)))
               enddo
            enddo
         else
!-----------------
! Non-hydrostatic:
!-----------------
           do i=is,ie
              phis(i,km+1) = hs(i,j)
              do k=km,1,-1
                 phis(i,k) = phis(i,k+1) - grav*delz(i,j,k)
              enddo
           enddo
           do i=is,ie
              te_2d(i,j) = 0.
           enddo
           do k=1,km
              do i=is,ie
! KE using 3D winds:
                 te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*( rgama*pt(i,j,k)*pkz(i,j,k) +  &
                              0.5*(phis(i,k)+phis(i,k+1) + w(i,j,k)**2 + 0.5*rsin2(i,j)*( &
                              u(i,j,k)**2+u(i,j+1,k)**2 + v(i,j,k)**2+v(i+1,j,k)**2 -  &
                             (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s(i,j))))
              enddo
           enddo
         endif

       endif
      enddo
    endif

!$omp parallel do default(shared)
      do j=js,je
         do i=is,ie
            zsum1(i,j) = pkz(i,j,1)*delp(i,j,1)
         enddo
         do k=2,km
            do i=is,ie
               zsum1(i,j) = zsum1(i,j) + pkz(i,j,k)*delp(i,j,k)
            enddo
         enddo

         do i=is,ie
            zsum0(i,j) = ptop*(pk(i,j,1)-pk(i,j,km+1)) + zsum1(i,j)
            te_2d(i,j) = te0_2d(i,j) - te_2d(i,j)
         enddo
      enddo

         tpe = consv*g_sum(te_2d, is, ie, js, je, ng, area, 0)
      E_Flux = tpe / (grav*pdt*4.*pi*radius**2)    ! unit: W/m**2
                                                   ! Note pdt is "phys" time step

      if ( hydrostatic ) then
           dtmp = tpe / (cp*g_sum(zsum0,  is, ie, js, je, ng, area, 0))
      else
           dtmp = tpe / (cv*g_sum(zsum1,  is, ie, js, je, ng, area, 0))
      endif
!-------------------------------------------------------------------------------
! One may use this quick fix to ensure reproducibility at the expense of a lower
! floating precision; this is fine for the TE correction
!-------------------------------------------------------------------------------
      if ( reproduce_sum ) dtmp = real(dtmp, 4) ! convert to 4-byte real
  else
      dtmp   = 0.
      E_Flux = 0.
  endif        ! end consv check

  if ( te_map ) then
!$omp parallel do default(shared) private(gz, tpe, tmp, dlnp)
      do j=js,je
         do i=is,ie
            gz(i) = hs(i,j)
         enddo
         do k=km,1,-1
            do i=is,ie
               tpe = te(i,j,k) - gz(i) - 0.25*rsin2(i,j)*(    &
                     u(i,j,k)**2+u(i,j+1,k)**2 + v(i,j,k)**2+v(i+1,j,k)**2 -  &
                    (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s(i,j) )
               dlnp = rg*(peln(i,k+1,j) - peln(i,k,j))
#ifdef CONVERT_T
               tmp = tpe / ((cp - pe(i,k,j)*dlnp/delp(i,j,k))*(1.+r_vir*q(i,j,k,sphum)) )
               pt(i,j,k) =  tmp + dtmp*pkz(i,j,k) / (1.+r_vir*q(i,j,k,sphum))
               gz(i) = gz(i) + dlnp*tmp*(1.+r_vir*q(i,j,k,sphum))
#else
               tmp = tpe / (cp - pe(i,k,j)*dlnp/delp(i,j,k))
               pt(i,j,k) = cp*(tmp/pkz(i,j,k) + dtmp)
               gz(i) = gz(i) + dlnp*tmp
#endif
            enddo
         enddo           ! end k-loop
      enddo
  else
    if ( remap_t ) then
!$omp parallel do default(shared) 
      do k=1,km
         do j=js,je
            do i=is,ie
#ifdef CONVERT_T
               pt(i,j,k) = (pt(i,j,k) + dtmp*pkz(i,j,k))/(1.+r_vir*q(i,j,k,sphum))
#else
               pt(i,j,k) = cp*(pt(i,j,k)/pkz(i,j,k) + dtmp)
#endif
            enddo
         enddo   
      enddo
    else
!$omp parallel do default(shared) 
      do k=1,km
         do j=js,je
            do i=is,ie
#ifdef CONVERT_T
               pt(i,j,k) = (rcp*pt(i,j,k) + dtmp)*pkz(i,j,k)/(1.+r_vir*q(i,j,k,sphum))
#else
               pt(i,j,k) = pt(i,j,k) + cp*dtmp
#endif
            enddo
         enddo   
      enddo
    endif
  endif

! call cubed_to_latlon(u, v, ua, va, dx, dy, rdxa, rdya, km, 1)

 end subroutine Lagrangian_to_Eulerian


 subroutine compute_total_energy(is, ie, js, je, isd, ied, jsd, jed, km,       &
                                 u, v, w, delz, pt, delp, q, qc, pe, peln, hs, &
                                 cp, rg, hlv, te_2d, ua, va, teq, &
                                 moist_phys, sphum, hydrostatic, id_te)
!------------------------------------------------------
! Compute vertically integrated total energy per column
!------------------------------------------------------
! !INPUT PARAMETERS:
   integer,  intent(in):: km, is, ie, js, je, isd, ied, jsd, jed, id_te, sphum
   real, intent(inout), dimension(isd:ied,jsd:jed,km):: ua, va
   real, intent(in), dimension(isd:ied,jsd:jed,km):: pt, delp
   real, intent(in), dimension(isd:ied,jsd:jed,km,sphum):: q
   real, intent(in), dimension(is:ie,js:je,km):: qc
   real, intent(inout)::  u(isd:ied,  jsd:jed+1,km)
   real, intent(inout)::  v(isd:ied+1,jsd:jed,  km)
   real, intent(in)::  w(isd:ied,jsd:jed,km)   ! vertical velocity (m/s)
   real, intent(in):: delz(is:ie,js:je,km)
   real, intent(in):: hs(isd:ied,jsd:jed)  ! surface geopotential
   real, intent(in)::   pe(is-1:ie+1,km+1,js-1:je+1) ! pressure at layer edges
   real, intent(in):: peln(is:ie,km+1,js:je)  ! log(pe)
   real, intent(in):: cp, rg, hlv
   logical, intent(in):: moist_phys, hydrostatic
! Output:
   real, intent(out):: te_2d(is:ie,js:je)   ! vertically integrated TE
   real, intent(out)::   teq(is:ie,js:je)   ! Moist TE
! Local
   real, dimension(is:ie,km):: tv
   real  phiz(is:ie,km+1)
   real cv
   integer i, j, k

   cv = cp - rg

!----------------------
! Output lat-lon winds:
!----------------------
!  call cubed_to_latlon(u, v, ua, va, dx, dy, rdxa, rdya, km)

!$omp parallel do default(shared) private(phiz, tv)
  do j=js,je

     if ( hydrostatic ) then

        do i=is,ie
           phiz(i,km+1) = hs(i,j)
        enddo
        do k=km,1,-1
           do i=is,ie
                tv(i,k) = pt(i,j,k)*(1.+qc(i,j,k))
              phiz(i,k) = phiz(i,k+1) + rg*tv(i,k)*(peln(i,k+1,j)-peln(i,k,j))
           enddo
        enddo

        do i=is,ie
           te_2d(i,j) = pe(i,km+1,j)*phiz(i,km+1) - pe(i,1,j)*phiz(i,1)
        enddo

        do k=1,km
           do i=is,ie
              te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(cp*tv(i,k) +            &
                           0.25*rsin2(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +      &
                                            v(i,j,k)**2+v(i+1,j,k)**2 -      &
                       (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s(i,j)))
           enddo
        enddo

     else
!-----------------
! Non-hydrostatic:
!-----------------
     do i=is,ie
        phiz(i,km+1) = hs(i,j)
        do k=km,1,-1
           phiz(i,k) = phiz(i,k+1) - grav*delz(i,j,k)
        enddo
     enddo
     do i=is,ie
        te_2d(i,j) = 0.
     enddo
     do k=1,km
        do i=is,ie
           te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*( cv*pt(i,j,k)*(1.+qc(i,j,k)) +  &
                        0.5*(phiz(i,k)+phiz(i,k+1)+w(i,j,k)**2+0.5*rsin2(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                        v(i,j,k)**2+v(i+1,j,k)**2-(u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s(i,j))))
        enddo
     enddo
     endif
  enddo

!-------------------------------------
! Doganostics computation for moist TE
!-------------------------------------
  if( id_te>0 ) then
!$omp parallel do default(shared) 
      do j=js,je
         do i=is,ie
            teq(i,j) = te_2d(i,j)
         enddo
      enddo
      if ( moist_phys ) then
!$omp parallel do default(shared) 
           do k=1,km
              do j=js,je
                 do i=is,ie
                    teq(i,j) = teq(i,j) + hlv*q(i,j,k,sphum)*delp(i,j,k)
                 enddo
              enddo
           enddo
      endif
!     do j=js,je
!        do i=is,ie
!           teq(i,j) = teq(i,j) / (pe(i,km,j) - pe(i,1,j))
!        enddo
!     enddo
   endif

  end subroutine compute_total_energy


  subroutine pkez(km, ifirst, ilast, jfirst, jlast, &
                  pe, pk, akap, peln, pkz)

! !INPUT PARAMETERS:
   integer, intent(in):: km
   integer, intent(in):: ifirst, ilast        ! Latitude strip
   integer, intent(in):: jfirst, jlast        ! Latitude strip
   real, intent(in):: akap
   real, intent(in):: pe(ifirst-1:ilast+1,km+1,jfirst-1:jlast+1)
   real, intent(in):: pk(ifirst:ilast,jfirst:jlast,km+1)
! !OUTPUT
   real, intent(out):: pkz(ifirst:ilast,jfirst:jlast,km)
   real, intent(inout):: peln(ifirst:ilast, km+1, jfirst:jlast)   ! log (pe)
! Local
   real pk2(ifirst:ilast, km+1)
   real pek
   real lnp
   real ak1
   integer i, j, k

   ak1 = (akap + 1.) / akap

!$omp parallel do default(shared) private(lnp, pek, pk2)
   do j=jfirst, jlast
        pek = pk(ifirst,j,1)
        do i=ifirst, ilast
           pk2(i,1) = pek
        enddo

        do k=2,km+1
           do i=ifirst, ilast
!             peln(i,k,j) =  log(pe(i,k,j))
              pk2(i,k) =  pk(i,j,k)
           enddo
        enddo

!---- GFDL modification
       if( ptop < ptop_min ) then
           do i=ifirst, ilast
               peln(i,1,j) = peln(i,2,j) - ak1
           enddo
       else
           lnp = log( ptop )
           do i=ifirst, ilast
              peln(i,1,j) = lnp
           enddo
       endif
!---- GFDL modification

       do k=1,km
          do i=ifirst, ilast
             pkz(i,j,k) = (pk2(i,k+1) - pk2(i,k) )  /  &
                          (akap*(peln(i,k+1,j) - peln(i,k,j)) )
          enddo
       enddo
    enddo

 end subroutine pkez



 subroutine remap_z(km, pe1, q1, kn, pe2, q2, i1, i2, iv, kord)

! !INPUT PARAMETERS:
      integer, intent(in) :: i1                ! Starting longitude
      integer, intent(in) :: i2                ! Finishing longitude
      integer, intent(in) :: kord              ! Method order
      integer, intent(in) :: km                ! Original vertical dimension
      integer, intent(in) :: kn                ! Target vertical dimension
      integer, intent(in) :: iv

      real, intent(in) ::  pe1(i1:i2,km+1)     ! height at layer edges 
                                               ! (from model top to bottom surface)
      real, intent(in) ::  pe2(i1:i2,kn+1)     ! hieght at layer edges 
                                               ! (from model top to bottom surface)
      real, intent(in) ::  q1(i1:i2,km)        ! Field input

! !INPUT/OUTPUT PARAMETERS:
      real, intent(inout)::  q2(i1:i2,kn)      ! Field output

! !LOCAL VARIABLES:
      real  dp1(  i1:i2,km)
      real   q4(4,i1:i2,km)
      real   pl, pr, qsum, delp, esl
      integer i, k, l, m, k0

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)      ! negative
            q4(1,i,k) = q1(i,k)
         enddo
      enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        call  cs_profile( q4, dp1, km, i1, i2, iv, kord )
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

! Mapping
      do 1000 i=i1,i2
         k0 = 1
      do 555 k=1,kn
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(i,k) <= pe1(i,l) .and. pe2(i,k) >= pe1(i,l+1)) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if(pe2(i,k+1) >= pe1(i,l+1)) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
          else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if(pe2(i,k+1) < pe1(i,m+1) ) then
! Whole layer..
                    qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                    delp = pe2(i,k+1)-pe1(i,m)
                    esl = delp / dp1(i,m)
                    qsum = qsum + delp*(q4(2,i,m)+0.5*esl*               &
                         (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                    k0 = m
                 goto 123
                 endif
              enddo
              goto 123
           endif
      endif
100   continue
123   q2(i,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
1000  continue

 end subroutine remap_z


 subroutine map1_ppm( km,   pe1,    q1,                 &
                      kn,   pe2,    q2,   i1, i2,       &
                      j,    ibeg, iend, jbeg, jend, iv,  kord)
 integer, intent(in) :: i1                ! Starting longitude
 integer, intent(in) :: i2                ! Finishing longitude
 integer, intent(in) :: iv                ! Mode: 0 == constituents  1 == ???
                                          !       2 == remap temp with cs scheme
 integer, intent(in) :: kord              ! Method order
 integer, intent(in) :: j                 ! Current latitude
 integer, intent(in) :: ibeg, iend, jbeg, jend
 integer, intent(in) :: km                ! Original vertical dimension
 integer, intent(in) :: kn                ! Target vertical dimension
 real, intent(in) ::  pe1(i1:i2,km+1)  ! pressure at layer edges 
                                       ! (from model top to bottom surface)
                                       ! in the original vertical coordinate
 real, intent(in) ::  pe2(i1:i2,kn+1)  ! pressure at layer edges 
                                       ! (from model top to bottom surface)
                                       ! in the new vertical coordinate
 real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
 real, intent(inout)::  q2(ibeg:iend,jbeg:jend,kn) ! Field output

! !DESCRIPTION:
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
! !LOCAL VARIABLES:
   real    dp1(i1:i2,km)
   real   q4(4,i1:i2,km)
   real    pl, pr, qsum, dp, esl
   integer i, k, l, m, k0

   do k=1,km
      do i=i1,i2
         dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         q4(1,i,k) = q1(i,j,k)
      enddo
   enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        call  cs_profile( q4, dp1, km, i1, i2, iv, kord )
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

  do i=i1,i2
     k0 = 1
     do 555 k=1,kn
      do l=k0,km
! locate the top edge: pe2(i,k)
      if( pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1) ) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if( pe2(i,k+1) <= pe1(i,l+1) ) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,j,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
         else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if( pe2(i,k+1) > pe1(i,m+1) ) then
! Whole layer
                     qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                     dp = pe2(i,k+1)-pe1(i,m)
                     esl = dp / dp1(i,m)
                     qsum = qsum + dp*(q4(2,i,m)+0.5*esl*               &
                           (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                     k0 = m
                     goto 123
                 endif
              enddo
              goto 123
         endif
      endif
      enddo
123   q2(i,j,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
  enddo

 end subroutine map1_ppm


 subroutine map1_q2(km,   pe1,   q1,            &
                    kn,   pe2,   q2,   dp2,     &
                    i1,   i2,    iv,   kord, j, &
                    ibeg, iend, jbeg, jend )


! !INPUT PARAMETERS:
      integer, intent(in) :: j
      integer, intent(in) :: i1, i2
      integer, intent(in) :: ibeg, iend, jbeg, jend
      integer, intent(in) :: iv                ! Mode: 0 ==  constituents 1 == ???
      integer, intent(in) :: kord
      integer, intent(in) :: km                ! Original vertical dimension
      integer, intent(in) :: kn                ! Target vertical dimension

      real, intent(in) ::  pe1(i1:i2,km+1)     ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the original vertical coordinate
      real, intent(in) ::  pe2(i1:i2,kn+1)     ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the new vertical coordinate
      real, intent(in) ::  q1(ibeg:iend,jbeg:jend,km) ! Field input
      real, intent(in) ::  dp2(i1:i2,kn)
! !INPUT/OUTPUT PARAMETERS:
      real, intent(inout):: q2(i1:i2,kn) ! Field output
! !LOCAL VARIABLES:
      real   dp1(i1:i2,km)
      real   q4(4,i1:i2,km)
      real   pl, pr, qsum, dp, esl

      integer i, k, l, m, k0

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
            q4(1,i,k) = q1(i,j,k)
         enddo
      enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        call  cs_profile( q4, dp1, km, i1, i2, iv, kord )
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

! Mapping
      do 1000 i=i1,i2
         k0 = 1
      do 555 k=1,kn
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1)) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if(pe2(i,k+1) <= pe1(i,l+1)) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
          else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if(pe2(i,k+1) > pe1(i,m+1) ) then
                                                   ! Whole layer..
                    qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                     dp = pe2(i,k+1)-pe1(i,m)
                    esl = dp / dp1(i,m)
                   qsum = qsum + dp*(q4(2,i,m)+0.5*esl*               &
                       (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                   k0 = m
                   goto 123
                 endif
              enddo
              goto 123
          endif
      endif
100   continue
123   q2(i,k) = qsum / dp2(i,k)
555   continue
1000  continue

 end subroutine map1_q2



 subroutine remap_2d(km,   pe1,   q1,        &
                     kn,   pe2,   q2,        &
                     i1,   i2,    iv,   kord )
   integer, intent(in):: i1, i2
   integer, intent(in):: iv               ! Mode: 0 ==  constituents 1 ==others
   integer, intent(in):: kord
   integer, intent(in):: km               ! Original vertical dimension
   integer, intent(in):: kn               ! Target vertical dimension
   real, intent(in):: pe1(i1:i2,km+1)     ! pressure at layer edges 
                                          ! (from model top to bottom surface)
                                          ! in the original vertical coordinate
   real, intent(in):: pe2(i1:i2,kn+1)     ! pressure at layer edges 
                                          ! (from model top to bottom surface)
                                          ! in the new vertical coordinate
   real, intent(in) :: q1(i1:i2,km) ! Field input
   real, intent(out):: q2(i1:i2,kn) ! Field output
! !LOCAL VARIABLES:
   real   dp1(i1:i2,km)
   real   q4(4,i1:i2,km)
   real   pl, pr, qsum, dp, esl
   integer i, k, l, m, k0

   do k=1,km
      do i=i1,i2
          dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         q4(1,i,k) = q1(i,k)
      enddo
   enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        call  cs_profile( q4, dp1, km, i1, i2, iv, kord )
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

   do i=i1,i2
      k0 = 1
      do 555 k=1,kn
         if( pe2(i,k+1) <= pe1(i,1) ) then
! Entire grid above old ptop
             q2(i,k) = q4(2,i,1)
         elseif( pe2(i,k) < pe1(i,1) .and. pe2(i,k+1)>pe1(i,1) ) then
! Partially above old ptop:
             q2(i,k) = q1(i,1)
         else
           do l=k0,km
! locate the top edge: pe2(i,k)
           if( pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1) ) then
               pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
               if(pe2(i,k+1) <= pe1(i,l+1)) then
! entire new grid is within the original grid
                  pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
                  q2(i,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                          *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
                  k0 = l
                  goto 555
               else
! Fractional area...
                 qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                         q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                        (r3*(1.+pl*(1.+pl))))
                 do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                    if(pe2(i,k+1) > pe1(i,m+1) ) then
                                                   ! Whole layer..
                       qsum = qsum + dp1(i,m)*q4(1,i,m)
                    else
                       dp = pe2(i,k+1)-pe1(i,m)
                      esl = dp / dp1(i,m)
                      qsum = qsum + dp*(q4(2,i,m)+0.5*esl*               &
                            (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                      k0 = m
                      goto 123
                    endif
                 enddo
                 goto 123
               endif
           endif
           enddo
123        q2(i,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
         endif
555   continue
   enddo

 end subroutine remap_2d


 subroutine cs_profile(a4, delp, km, i1, i2, iv, kord)
! Optimized vertical profile reconstruction:
! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
 integer, intent(in):: i1, i2
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
 integer, intent(in):: kord
 real, intent(in)   :: delp(i1:i2,km)     ! layer pressure thickness
 real, intent(inout):: a4(4,i1:i2,km)     ! Interpolated values
!-----------------------------------------------------------------------
 logical:: extm(i1:i2,km) 
 real  gam(i1:i2,km)
 real    q(i1:i2,km+1)
 real   d4(i1:i2)
 real   bet, a_bot, grat 
 real   pmp_1, lac_1, pmp_2, lac_2
 integer i, k, im

  do i=i1,i2
         grat = delp(i,2) / delp(i,1)   ! grid ratio
          bet = grat*(grat+0.5)
       q(i,1) = ( (grat+grat)*(grat+1.)*a4(1,i,1) + a4(1,i,2) ) / bet
     gam(i,1) = ( 1. + grat*(grat+1.5) ) / bet
  enddo

  do k=2,km
     do i=i1,i2
           d4(i) = delp(i,k-1) / delp(i,k)
             bet =  2. + d4(i) + d4(i) - gam(i,k-1)
          q(i,k) = ( 3.*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - q(i,k-1) )/bet
        gam(i,k) = d4(i) / bet
     enddo
  enddo
 
  do i=i1,i2
         a_bot = 1. + d4(i)*(d4(i)+1.5)
     q(i,km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1,i,km)+a4(1,i,km-1)-a_bot*q(i,km))  &
               / ( d4(i)*(d4(i)+0.5) - a_bot*gam(i,km) )
  enddo

  do k=km,1,-1
     do i=i1,i2
        q(i,k) = q(i,k) - gam(i,k)*q(i,k+1)
     enddo
  enddo

!------------------
! Apply constraints
!------------------
  im = i2 - i1 + 1

! Apply *large-scale* constraints 
  do i=i1,i2
     q(i,2) = min( q(i,2), max(a4(1,i,1), a4(1,i,2)) )
     q(i,2) = max( q(i,2), min(a4(1,i,1), a4(1,i,2)) )
  enddo

  do k=2,km
     do i=i1,i2
        gam(i,k) = a4(1,i,k) - a4(1,i,k-1)
     enddo
  enddo

! Interior:
  if ( abs(kord)<11 ) then
  do k=3,km-1
     do i=i1,i2
        if ( gam(i,k-1)*gam(i,k+1)>0. ) then
! Apply large-scale constraint to ALL fields if not local max/min
             q(i,k) = min( q(i,k), max(a4(1,i,k-1),a4(1,i,k)) )
             q(i,k) = max( q(i,k), min(a4(1,i,k-1),a4(1,i,k)) )
        else
          if ( gam(i,k-1) > 0. ) then
! There exists a local max
               q(i,k) = max(q(i,k), min(a4(1,i,k-1),a4(1,i,k)))
          else
! There exists a local min
                 q(i,k) = min(q(i,k), max(a4(1,i,k-1),a4(1,i,k)))
               if ( iv==0 ) q(i,k) = max(0., q(i,k))
          endif
        endif
     enddo
  enddo
  else
! abs(kord) >=11
  do k=3,km-1
     do i=i1,i2
        if ( gam(i,k-1)*gam(i,k+1) > 0. ) then
! Apply large-scale constraint to ALL fields if not local max/min
             q(i,k) = min( q(i,k), max(a4(1,i,k-1),a4(1,i,k)) )
             q(i,k) = max( q(i,k), min(a4(1,i,k-1),a4(1,i,k)) )
        else
          if ( gam(i,k-1) > 0. ) then
! There exists a local max
               q(i,k) = max(q(i,k), min(a4(1,i,k-1),a4(1,i,k)))
               q(i,k) = min(q(i,k), a4(1,i,k-1)+0.5*gam(i,k-1),     &
                                    a4(1,i,k  )-0.5*gam(i,k+1) )
          else
! There exists a local min
            if ( iv==0 ) then
                 q(i,k) = min(q(i,k), max(a4(1,i,k-1),a4(1,i,k)))
                 q(i,k) = max(q(i,k), a4(1,i,k-1)+0.5*gam(i,k-1),     &
                                  0., a4(1,i,k  )-0.5*gam(i,k+1) )
            else
                 q(i,k) = min(q(i,k), max(a4(1,i,k-1),a4(1,i,k)))
                 q(i,k) = max(q(i,k), a4(1,i,k-1)+0.5*gam(i,k-1),     &
                                      a4(1,i,k  )-0.5*gam(i,k+1) )
            endif
          endif
        endif
     enddo
  enddo
  endif

! Bottom:
  do i=i1,i2
     q(i,km) = min( q(i,km), max(a4(1,i,km-1), a4(1,i,km)) )
     q(i,km) = max( q(i,km), min(a4(1,i,km-1), a4(1,i,km)) )
  enddo

  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
     enddo
  enddo

  do k=2,km-1
     do i=i1,i2
        if ( gam(i,k)*gam(i,k+1) > 0.0 ) then
             extm(i,k) = .false. 
        else
             extm(i,k) = .true.
        endif
     enddo
  enddo

!---------------------------
! Apply subgrid constraints:
!---------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
! Top 2 and bottom 2 layers always use monotonic mapping

  if ( iv==0 ) then
     do i=i1,i2
        a4(2,i,1) = max(0., a4(2,i,1))
     enddo
  elseif ( iv==-1 ) then 
      do i=i1,i2
         if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
      enddo
  elseif ( abs(iv)==2 ) then
     do i=i1,i2
        a4(2,i,1) = a4(1,i,1)
        a4(3,i,1) = a4(1,i,1)
        a4(4,i,1) = 0.
     enddo
  endif

  if ( abs(iv)/=2 ) then
     do i=i1,i2
        a4(4,i,1) = 3.*(2.*a4(1,i,1) - (a4(2,i,1)+a4(3,i,1)))
     enddo
     call cs_limiters(im, extm(i1,1), a4(1,i1,1), 1)
  endif

! k=2
   do i=i1,i2
      a4(4,i,2) = 3.*(2.*a4(1,i,2) - (a4(2,i,2)+a4(3,i,2)))
   enddo
   call cs_limiters(im, extm(i1,2), a4(1,i1,2), 2)

!-------------------------------------
! Huynh's 2nd constraint for interior:
!-------------------------------------
  do k=3,km-2
     if ( abs(kord)<9 ) then
       do i=i1,i2
! Left  edges
          pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
          lac_1 = pmp_1 + 1.5*gam(i,k+2)
          a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                         max(a4(1,i,k), pmp_1, lac_1) )
! Right edges
          pmp_2 = a4(1,i,k) + 2.*gam(i,k)
          lac_2 = pmp_2 - 1.5*gam(i,k-1)
          a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                         max(a4(1,i,k), pmp_2, lac_2) )

          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo

     elseif ( abs(kord)==9 ) then
       do i=i1,i2
          if ( extm(i,k) .and. (extm(i,k-1).or.extm(i,k+1)) ) then  ! c90_mp122
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==10 ) then
       do i=i1,i2
          if( extm(i,k) ) then
              if( extm(i,k-1) .or. extm(i,k+1) ) then
! grid-scale 2-delta-z wave detected
                   a4(2,i,k) = a4(1,i,k)
                   a4(3,i,k) = a4(1,i,k)
                   a4(4,i,k) = 0.
              else
! True local extremum
                a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
              endif
          else        ! not a local extremum
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     else      ! kord = 11, ...
       do i=i1,i2
         if ( extm(i,k) .and. (extm(i,k-1) .or. extm(i,k+1)) ) then
! Noisy region:
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
         else
              a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         endif
       enddo
     endif

! Additional constraint to ensure positivity
     if ( iv==0 ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 0)

  enddo      ! k-loop

!----------------------------------
! Bottom layer subgrid constraints:
!----------------------------------
  if ( iv==0 ) then
     do i=i1,i2
        a4(3,i,km) = max(0., a4(3,i,km))
     enddo
  elseif ( iv<0 ) then 
      do i=i1,i2
         if ( a4(3,i,km)*a4(1,i,km) <= 0. )  a4(3,i,km) = 0.
      enddo
  endif

  do k=km-1,km
     do i=i1,i2
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
     if(k==(km-1)) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 2)
     if(k== km   ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 1)
  enddo

 end subroutine cs_profile



 subroutine cs_limiters(im, extm, a4, iv)
 integer, intent(in) :: im
 integer, intent(in) :: iv
 logical, intent(in) :: extm(im)
 real , intent(inout) :: a4(4,im)   ! PPM array
! !LOCAL VARIABLES:
 real  da1, da2, a6da
 integer i

 if ( iv==0 ) then
! Positive definite constraint
    do i=1,im
    if( a4(1,i)<=0.) then
        a4(2,i) = a4(1,i)
        a4(3,i) = a4(1,i)
        a4(4,i) = 0.
    else
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
         if( (a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12) < 0. ) then
! local minimum is negative
             if( a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i) ) then
                 a4(3,i) = a4(1,i)
                 a4(2,i) = a4(1,i)
                 a4(4,i) = 0.
             elseif( a4(3,i) > a4(2,i) ) then
                 a4(4,i) = 3.*(a4(2,i)-a4(1,i))
                 a4(3,i) = a4(2,i) - a4(4,i)
             else
                 a4(4,i) = 3.*(a4(3,i)-a4(1,i))
                 a4(2,i) = a4(3,i) - a4(4,i)
             endif
         endif
      endif
    endif
    enddo
 elseif ( iv==1 ) then
    do i=1,im
      if( (a4(1,i)-a4(2,i))*(a4(1,i)-a4(3,i))>=0. ) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
    enddo
 else
! Standard PPM constraint
    do i=1,im
      if( extm(i) ) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
    enddo
 endif
 end subroutine cs_limiters



 subroutine ppm_profile(a4, delp, km, i1, i2, iv, kord)

! !INPUT PARAMETERS:
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
                               ! iv = 2: temp (if remap_t) and w (iv=-2)
 integer, intent(in):: i1      ! Starting longitude
 integer, intent(in):: i2      ! Finishing longitude
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: kord    ! Order (or more accurately method no.):
                               ! 
 real , intent(in):: delp(i1:i2,km)     ! layer pressure thickness

! !INPUT/OUTPUT PARAMETERS:
 real , intent(inout):: a4(4,i1:i2,km)  ! Interpolated values

! DESCRIPTION:
!
!   Perform the piecewise parabolic reconstruction
! 
! !REVISION HISTORY: 
! S.-J. Lin   revised at GFDL 2007
!-----------------------------------------------------------------------
! local arrays:
      real    dc(i1:i2,km)
      real    h2(i1:i2,km)
      real  delq(i1:i2,km)
      real   df2(i1:i2,km)
      real    d4(i1:i2,km)

! local scalars:
      integer i, k, km1, lmt, it
      real  fac
      real  a1, a2, c1, c2, c3, d1, d2
      real  qm, dq, lac, qmp, pmp

      km1 = km - 1
       it = i2 - i1 + 1

      do k=2,km
         do i=i1,i2
            delq(i,k-1) =   a4(1,i,k) - a4(1,i,k-1)
              d4(i,k  ) = delp(i,k-1) + delp(i,k)
         enddo
      enddo

      do k=2,km1
         do i=i1,i2
                 c1  = (delp(i,k-1)+0.5*delp(i,k))/d4(i,k+1)
                 c2  = (delp(i,k+1)+0.5*delp(i,k))/d4(i,k)
            df2(i,k) = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /      &
                                    (d4(i,k)+delp(i,k+1))
            dc(i,k) = sign( min(abs(df2(i,k)),              &
                            max(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))-a4(1,i,k),  &
                  a4(1,i,k)-min(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))), df2(i,k) )
         enddo
      enddo

!-----------------------------------------------------------
! 4th order interpolation of the provisional cell edge value
!-----------------------------------------------------------

      do k=3,km1
         do i=i1,i2
            c1 = delq(i,k-1)*delp(i,k-1) / d4(i,k)
            a1 = d4(i,k-1) / (d4(i,k) + delp(i,k-1))
            a2 = d4(i,k+1) / (d4(i,k) + delp(i,k))
            a4(2,i,k) = a4(1,i,k-1) + c1 + 2./(d4(i,k-1)+d4(i,k+1)) *    &
                      ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -          &
                        delp(i,k-1)*a1*dc(i,k  ) )
         enddo
      enddo

      if(km>8 .and. kord>4) call steepz(i1, i2, km, a4, df2, dc, delq, delp, d4)

! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top
      do i=i1,i2
         d1 = delp(i,1)
         d2 = delp(i,2)
         qm = (d2*a4(1,i,1)+d1*a4(1,i,2)) / (d1+d2)
         dq = 2.*(a4(1,i,2)-a4(1,i,1)) / (d1+d2)
         c1 = 4.*(a4(2,i,3)-qm-d2*dq) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) )
         c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1*d1)
         a4(2,i,2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
! Top edge:
!-------------------------------------------------------
         a4(2,i,1) = d1*(2.*c1*d1**2-c3) + a4(2,i,2)
!-------------------------------------------------------
!        a4(2,i,1) = (12./7.)*a4(1,i,1)-(13./14.)*a4(1,i,2)+(3./14.)*a4(1,i,3)
!-------------------------------------------------------
! No over- and undershoot condition
         a4(2,i,2) = max( a4(2,i,2), min(a4(1,i,1), a4(1,i,2)) )
         a4(2,i,2) = min( a4(2,i,2), max(a4(1,i,1), a4(1,i,2)) )
         dc(i,1) =  0.5*(a4(2,i,2) - a4(1,i,1))
      enddo

! Enforce monotonicity  within the top layer

      if( iv==0 ) then
         do i=i1,i2
            a4(2,i,1) = max(0., a4(2,i,1))
            a4(2,i,2) = max(0., a4(2,i,2))
         enddo 
      elseif( iv==-1 ) then
         do i=i1,i2
            if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
         enddo
      elseif( abs(iv)==2 ) then
         do i=i1,i2
            a4(2,i,1) = a4(1,i,1)
            a4(3,i,1) = a4(1,i,1)
         enddo
      endif

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
      do i=i1,i2
         d1 = delp(i,km)
         d2 = delp(i,km1)
         qm = (d2*a4(1,i,km)+d1*a4(1,i,km1)) / (d1+d2)
         dq = 2.*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2)
         c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
         c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1*d1)
         a4(2,i,km) = qm - c1*d1*d2*(d2+3.*d1)
! Bottom edge:
!-----------------------------------------------------
         a4(3,i,km) = d1*(8.*c1*d1**2-c3) + a4(2,i,km)
!        dc(i,km) = 0.5*(a4(3,i,km) - a4(1,i,km))
!-----------------------------------------------------
!        a4(3,i,km) = (12./7.)*a4(1,i,km)-(13./14.)*a4(1,i,km-1)+(3./14.)*a4(1,i,km-2)
! No over- and under-shoot condition
         a4(2,i,km) = max( a4(2,i,km), min(a4(1,i,km), a4(1,i,km1)) )
         a4(2,i,km) = min( a4(2,i,km), max(a4(1,i,km), a4(1,i,km1)) )
         dc(i,km) = 0.5*(a4(1,i,km) - a4(2,i,km))
      enddo


! Enforce constraint on the "slope" at the surface

#ifdef BOT_MONO
      do i=i1,i2
         a4(4,i,km) = 0
         if( a4(3,i,km) * a4(1,i,km) <= 0. ) a4(3,i,km) = 0.
         d1 = a4(1,i,km) - a4(2,i,km)
         d2 = a4(3,i,km) - a4(1,i,km)
         if ( d1*d2 < 0. ) then
              a4(2,i,km) = a4(1,i,km)
              a4(3,i,km) = a4(1,i,km)
         else
              dq = sign(min(abs(d1),abs(d2),0.5*abs(delq(i,km-1))), d1)
              a4(2,i,km) = a4(1,i,km) - dq
              a4(3,i,km) = a4(1,i,km) + dq
         endif
      enddo
#else
      if( iv==0 ) then
          do i=i1,i2
             a4(2,i,km) = max(0.,a4(2,i,km))
             a4(3,i,km) = max(0.,a4(3,i,km))
          enddo
      elseif( iv<0 ) then
          do i=i1,i2
             if( a4(1,i,km)*a4(3,i,km) <= 0. )  a4(3,i,km) = 0.
          enddo
      endif
#endif

   do k=1,km1
      do i=i1,i2
         a4(3,i,k) = a4(2,i,k+1)
      enddo
   enddo

!-----------------------------------------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
!-----------------------------------------------------------
! Top 2 and bottom 2 layers always use monotonic mapping
      do k=1,2
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

      if(kord >= 7) then
!-----------------------
! Huynh's 2nd constraint
!-----------------------
      do k=2,km1
         do i=i1,i2
! Method#1
!           h2(i,k) = delq(i,k) - delq(i,k-1)
! Method#2 - better
            h2(i,k) = 2.*(dc(i,k+1)/delp(i,k+1) - dc(i,k-1)/delp(i,k-1))  &
                     / ( delp(i,k)+0.5*(delp(i,k-1)+delp(i,k+1)) )        &
                     * delp(i,k)**2 
! Method#3
!!!            h2(i,k) = dc(i,k+1) - dc(i,k-1)
         enddo
      enddo

      fac = 1.5           ! original quasi-monotone

      do k=3,km-2
        do i=i1,i2
! Right edges
!        qmp   = a4(1,i,k) + 2.0*delq(i,k-1)
!        lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
!
         pmp   = 2.*dc(i,k)
         qmp   = a4(1,i,k) + pmp
         lac   = a4(1,i,k) + fac*h2(i,k-1) + dc(i,k)
         a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), qmp, lac)),    &
                                        max(a4(1,i,k), qmp, lac) )
! Left  edges
!        qmp   = a4(1,i,k) - 2.0*delq(i,k)
!        lac   = a4(1,i,k) + fac*h2(i,k+1) - 0.5*delq(i,k)
!
         qmp   = a4(1,i,k) - pmp
         lac   = a4(1,i,k) + fac*h2(i,k+1) - dc(i,k)
         a4(2,i,k) = min(max(a4(2,i,k),  min(a4(1,i,k), qmp, lac)),   &
                     max(a4(1,i,k), qmp, lac))
!-------------
! Recompute A6
!-------------
         a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
! Additional constraint to ensure positivity when kord=7
         if (iv == 0 .and. kord >= 6 )                      &
             call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 2)
      enddo

      else

         lmt = kord - 3
         lmt = max(0, lmt)
         if (iv == 0) lmt = min(2, lmt)

         do k=3,km-2
            if( kord /= 4) then
              do i=i1,i2
                 a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
              enddo
            endif
            if(kord/=6) call ppm_limiters(dc(i1,k), a4(1,i1,k), it, lmt)
         enddo
      endif

      do k=km1,km
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

 end subroutine ppm_profile


 subroutine ppm_limiters(dm, a4, itot, lmt)

! !INPUT PARAMETERS:
      real , intent(in):: dm(*)     ! the linear slope
      integer, intent(in) :: itot      ! Total Longitudes
      integer, intent(in) :: lmt       ! 0: Standard PPM constraint
                                       ! 1: Improved full monotonicity constraint (Lin)
                                       ! 2: Positive definite constraint
                                       ! 3: do nothing (return immediately)
! !INPUT/OUTPUT PARAMETERS:
      real , intent(inout) :: a4(4,*)   ! PPM array
                                           ! AA <-- a4(1,i)
                                           ! AL <-- a4(2,i)
                                           ! AR <-- a4(3,i)
                                           ! A6 <-- a4(4,i)
! !LOCAL VARIABLES:
      real  qmp
      real  da1, da2, a6da
      real  fmin
      integer i

! Developer: S.-J. Lin

      if ( lmt == 3 ) return

      if(lmt == 0) then
! Standard PPM constraint
      do i=1,itot
      if(dm(i) == 0.) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
      enddo

      elseif (lmt == 1) then

! Improved full monotonicity constraint (Lin 2004)
! Note: no need to provide first guess of A6 <-- a4(4,i)
      do i=1, itot
           qmp = 2.*dm(i)
         a4(2,i) = a4(1,i)-sign(min(abs(qmp),abs(a4(2,i)-a4(1,i))), qmp)
         a4(3,i) = a4(1,i)+sign(min(abs(qmp),abs(a4(3,i)-a4(1,i))), qmp)
         a4(4,i) = 3.*( 2.*a4(1,i) - (a4(2,i)+a4(3,i)) )
      enddo

      elseif (lmt == 2) then

! Positive definite constraint
      do i=1,itot
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
      fmin = a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12
         if( fmin < 0. ) then
         if(a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i)) then
            a4(3,i) = a4(1,i)
            a4(2,i) = a4(1,i)
            a4(4,i) = 0.
         elseif(a4(3,i) > a4(2,i)) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         else
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
         endif
      endif
      enddo

      endif

 end subroutine ppm_limiters



 subroutine steepz(i1, i2, km, a4, df2, dm, dq, dp, d4)
 integer, intent(in) :: km, i1, i2
   real , intent(in) ::  dp(i1:i2,km)       ! grid size
   real , intent(in) ::  dq(i1:i2,km)       ! backward diff of q
   real , intent(in) ::  d4(i1:i2,km)       ! backward sum:  dp(k)+ dp(k-1) 
   real , intent(in) :: df2(i1:i2,km)       ! first guess mismatch
   real , intent(in) ::  dm(i1:i2,km)       ! monotonic mismatch
! !INPUT/OUTPUT PARAMETERS:
      real , intent(inout) ::  a4(4,i1:i2,km)  ! first guess/steepened
! !LOCAL VARIABLES:
      integer i, k
      real  alfa(i1:i2,km)
      real     f(i1:i2,km)
      real   rat(i1:i2,km)
      real   dg2

! Compute ratio of dq/dp
      do k=2,km
         do i=i1,i2
            rat(i,k) = dq(i,k-1) / d4(i,k)
         enddo
      enddo

! Compute F
      do k=2,km-1
         do i=i1,i2
            f(i,k) =   (rat(i,k+1) - rat(i,k))                          &
                     / ( dp(i,k-1)+dp(i,k)+dp(i,k+1) )
         enddo
      enddo

      do k=3,km-2
         do i=i1,i2
         if(f(i,k+1)*f(i,k-1)<0. .and. df2(i,k)/=0.) then
            dg2 = (f(i,k+1)-f(i,k-1))*((dp(i,k+1)-dp(i,k-1))**2          &
                   + d4(i,k)*d4(i,k+1) )
            alfa(i,k) = max(0., min(0.5, -0.1875*dg2/df2(i,k)))
         else
            alfa(i,k) = 0.
         endif
         enddo
      enddo

      do k=4,km-2
         do i=i1,i2
            a4(2,i,k) = (1.-alfa(i,k-1)-alfa(i,k)) * a4(2,i,k) +         &
                        alfa(i,k-1)*(a4(1,i,k)-dm(i,k))    +             &
                        alfa(i,k)*(a4(1,i,k-1)+dm(i,k-1))
         enddo
      enddo

 end subroutine steepz



 subroutine rst_remap(km, kn, is,ie,js,je, isd,ied,jsd,jed, nq,  &
                      delp_r, u_r, v_r, w_r, delz_r, pt_r, q_r,      &
                      delp,   u,   v,   w,   delz,   pt,   q,        &
                      ak_r, bk_r, ak, bk, hydrostatic)
!------------------------------------
! Assuming hybrid sigma-P coordinate:
!------------------------------------
! !INPUT PARAMETERS:
  integer, intent(in):: km                    ! Restart z-dimension
  integer, intent(in):: kn                    ! Run time dimension
  integer, intent(in):: nq                    ! number of tracers (including h2o)
  integer, intent(in):: is,ie,isd,ied         ! starting & ending X-Dir index
  integer, intent(in):: js,je,jsd,jed         ! starting & ending Y-Dir index
  logical, intent(in):: hydrostatic
  real, intent(in) :: ak_r(km+1)
  real, intent(in) :: bk_r(km+1)
  real, intent(in) :: ak(kn+1)
  real, intent(in) :: bk(kn+1)
  real, intent(in):: delp_r(is:ie,js:je,km) ! pressure thickness
  real, intent(in)::   u_r(is:ie,  js:je+1,km)   ! u-wind (m/s)
  real, intent(in)::   v_r(is:ie+1,js:je  ,km)   ! v-wind (m/s)
  real, intent(inout)::  pt_r(is:ie,js:je,km)
  real, intent(in)::   w_r(is:ie,js:je,km)
  real, intent(in)::   q_r(is:ie,js:je,km,*)
  real, intent(inout)::delz_r(is:ie,js:je,km)
! Output:
  real, intent(out):: delp(isd:ied,jsd:jed,kn) ! pressure thickness
  real, intent(out)::  u(isd:ied  ,jsd:jed+1,kn)   ! u-wind (m/s)
  real, intent(out)::  v(isd:ied+1,jsd:jed  ,kn)   ! v-wind (m/s)
  real, intent(out)::  w(isd:ied  ,jsd:jed  ,kn)   ! vertical velocity (m/s)
  real, intent(out):: pt(isd:ied  ,jsd:jed  ,kn)   ! temperature
  real, intent(out):: q(isd:ied,jsd:jed,kn,*)
  real, intent(out):: delz(is:ie,js:je,kn)   ! delta-height (m)
!-----------------------------------------------------------------------
  real r_vir
  real ps(isd:ied,jsd:jed)  ! surface pressure
  real  pe1(is:ie,km+1)
  real  pe2(is:ie,kn+1)
  real  pv1(is:ie+1,km+1)
  real  pv2(is:ie+1,kn+1)

  integer i,j,k , iq
  integer, parameter:: kord=4

  r_vir = rvgas/rdgas - 1.

!$omp parallel do default(shared)
  do j=js,je
     do i=is,ie
        ps(i,j) = ak_r(1)
     enddo
  enddo

! this OpenMP do-loop setup cannot work in it's current form....
!rab!$omp parallel do default(shared)
  do k=1,km
     do j=js,je
        do i=is,ie
           ps(i,j) = ps(i,j) + delp_r(i,j,k)
        enddo
     enddo
  enddo

! only one cell is needed
  if ( square_domain ) then
      call mpp_update_domains(ps, domain,  whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
  else
      call mpp_update_domains(ps, domain, complete=.true.)
  endif

! Compute virtual Temp
!$omp parallel do default(shared)
  do k=1,km
     do j=js,je
        do i=is,ie
           pt_r(i,j,k) = pt_r(i,j,k) * (1.+r_vir*q_r(i,j,k,1))
        enddo
     enddo
  enddo

!$omp parallel do default(shared) private(pe1,  pe2, pv1, pv2)
  do 1000 j=js,je+1
!------
! map u
!------
     do k=1,km+1
        do i=is,ie
           pe1(i,k) = ak_r(k) + 0.5*bk_r(k)*(ps(i,j-1)+ps(i,j))
        enddo
     enddo

     do k=1,kn+1
        do i=is,ie
           pe2(i,k) = ak(k) + 0.5*bk(k)*(ps(i,j-1)+ps(i,j))
        enddo
     enddo

     call remap_2d(km, pe1, u_r(is:ie,j:j,1:km),       &
                   kn, pe2,   u(is:ie,j:j,1:kn),       &
                   is, ie, -1, kord)

  if ( j /= (je+1) )  then 

!---------------
! Hybrid sigma-p
!---------------
     do k=1,km+1
        do i=is,ie
           pe1(i,k) = ak_r(k) + bk_r(k)*ps(i,j)
        enddo
     enddo

     do k=1,kn+1
        do i=is,ie
           pe2(i,k) =   ak(k) + bk(k)*ps(i,j)
        enddo
     enddo

!-------------
! Compute delp
!-------------
      do k=1,kn
         do i=is,ie
            delp(i,j,k) = pe2(i,k+1) - pe2(i,k)
         enddo
      enddo

!----------------
! Map constituents
!----------------
      if( nq /= 0 ) then
          do iq=1,nq
             call remap_2d(km, pe1, q_r(is:ie,j:j,1:km,iq:iq),  &
                           kn, pe2,   q(is:ie,j:j,1:kn,iq:iq),  &
                           is, ie, 0, kord)
          enddo
      endif

      if ( .not. hydrostatic ) then
! Remap vertical wind:
         call remap_2d(km, pe1, w_r(is:ie,j:j,1:km),       &
                       kn, pe2,   w(is:ie,j:j,1:kn),       &
                       is, ie, -1, kord)
! Remap delz for hybrid sigma-p coordinate
         do k=1,km
            do i=is,ie
               delz_r(i,j,k) = -delz_r(i,j,k)/delp_r(i,j,k) ! ="specific volume"/grav
            enddo
         enddo
         call remap_2d(km, pe1, delz_r(is:ie,j:j,1:km),       &
                       kn, pe2,   delz(is:ie,j:j,1:kn),       &
                       is, ie, 1, kord)
         do k=1,kn
            do i=is,ie
               delz(i,j,k) = -delz(i,j,k)*delp(i,j,k)
            enddo
         enddo
      endif

! Geopotential conserving remap of virtual temperature:
       do k=1,km+1
          do i=is,ie
             pe1(i,k) = log(pe1(i,k))
          enddo
       enddo
       do k=1,kn+1
          do i=is,ie
             pe2(i,k) = log(pe2(i,k))
          enddo
       enddo

       call remap_2d(km, pe1, pt_r(is:ie,j:j,1:km),       &
                     kn, pe2,   pt(is:ie,j:j,1:kn),       &
                     is, ie, 1, kord)
!------
! map v
!------
       do k=1,km+1
          do i=is,ie+1
             pv1(i,k) = ak_r(k) + 0.5*bk_r(k)*(ps(i-1,j)+ps(i,j))
          enddo
       enddo
       do k=1,kn+1
          do i=is,ie+1
             pv2(i,k) = ak(k) + 0.5*bk(k)*(ps(i-1,j)+ps(i,j))
          enddo
       enddo

       call remap_2d(km, pv1, v_r(is:ie+1,j:j,1:km),       &
                     kn, pv2,   v(is:ie+1,j:j,1:kn),       &
                     is, ie+1, -1, kord)

  endif !(j < je+1)
1000  continue

!$omp parallel do default(shared) 
  do k=1,kn
     do j=js,je
        do i=is,ie
           pt(i,j,k) = pt(i,j,k) / (1.+r_vir*q(i,j,k,1))
        enddo
     enddo   
  enddo

 end subroutine rst_remap



 subroutine mappm(km, pe1, q1, kn, pe2, q2, i1, i2, iv, kord)

! IV = 0: constituents
! IV = 1: potential temp
! IV =-1: winds
 
! Mass flux preserving mapping: q1(im,km) -> q2(im,kn)
 
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate

 integer, intent(in):: i1, i2, km, kn, kord, iv
 real, intent(in ):: pe1(i1:i2,km+1), pe2(i1:i2,kn+1)
 real, intent(in )::  q1(i1:i2,km)
 real, intent(out)::  q2(i1:i2,kn)
! local
      real dp1(i1:i2,km)
      real a4(4,i1:i2,km)
      integer i, k, l
      integer k0, k1
      real pl, pr, tt, delp, qsum, dpsum, esl

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
            a4(1,i,k) = q1(i,k)
         enddo
      enddo

      if ( kord >7 ) then
           call  cs_profile( a4, dp1, km, i1, i2, iv, kord )
      else
           call ppm_profile( a4, dp1, km, i1, i2, iv, kord )
      endif

!------------------------------------
! Lowest layer: constant distribution
!------------------------------------
      do i=i1,i2
         a4(2,i,km) = q1(i,km)
         a4(3,i,km) = q1(i,km)
         a4(4,i,km) = 0.
      enddo

      do 5555 i=i1,i2
         k0 = 1
      do 555 k=1,kn

         if(pe2(i,k+1) .le. pe1(i,1)) then
! Entire grid above old ptop
            q2(i,k) = a4(2,i,1)
         elseif(pe2(i,k) .ge. pe1(i,km+1)) then
! Entire grid below old ps
            q2(i,k) = a4(3,i,km)
         elseif(pe2(i,k  ) .lt. pe1(i,1) .and.   &
                pe2(i,k+1) .gt. pe1(i,1))  then
! Part of the grid above ptop
            q2(i,k) = a4(1,i,1)
         else

         do 45 L=k0,km
! locate the top edge at pe2(i,k)
         if( pe2(i,k) .ge. pe1(i,L) .and.        &
             pe2(i,k) .le. pe1(i,L+1)    ) then
             k0 = L
             PL = (pe2(i,k)-pe1(i,L)) / dp1(i,L)
             if(pe2(i,k+1) .le. pe1(i,L+1)) then

! entire new grid is within the original grid
               PR = (pe2(i,k+1)-pe1(i,L)) / dp1(i,L)
               TT = r3*(PR*(PR+PL)+PL**2)
               q2(i,k) = a4(2,i,L) + 0.5*(a4(4,i,L)+a4(3,i,L)  &
                       - a4(2,i,L))*(PR+PL) - a4(4,i,L)*TT
              goto 555
             else
! Fractional area...
              delp = pe1(i,L+1) - pe2(i,k)
              TT   = r3*(1.+PL*(1.+PL))
              qsum = delp*(a4(2,i,L)+0.5*(a4(4,i,L)+            &
                     a4(3,i,L)-a4(2,i,L))*(1.+PL)-a4(4,i,L)*TT)
              dpsum = delp
              k1 = L + 1
             goto 111
             endif
         endif
45       continue

111      continue
         do 55 L=k1,km
         if( pe2(i,k+1) .gt. pe1(i,L+1) ) then

! Whole layer..

            qsum  =  qsum + dp1(i,L)*q1(i,L)
            dpsum = dpsum + dp1(i,L)
         else
           delp = pe2(i,k+1)-pe1(i,L)
           esl  = delp / dp1(i,L)
           qsum = qsum + delp * (a4(2,i,L)+0.5*esl*            &
                 (a4(3,i,L)-a4(2,i,L)+a4(4,i,L)*(1.-r23*esl)) )
          dpsum = dpsum + delp
           k0 = L
           goto 123
         endif
55       continue
        delp = pe2(i,k+1) - pe1(i,km+1)
        if(delp > 0.) then
! Extended below old ps
           qsum = qsum + delp * a4(3,i,km)
          dpsum = dpsum + delp
        endif
123     q2(i,k) = qsum / dpsum
      endif
555   continue
5555  continue

 end subroutine mappm

end module fv_mapz_mod
