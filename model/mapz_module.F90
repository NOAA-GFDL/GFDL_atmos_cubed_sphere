module mapz_module

  use grid_tools,   only: area, dx, dy, dxa, dya
  use grid_utils,   only: cubed_to_latlon, g_sum, ptop, ptop_min
  use fill_module,  only: fillz
  use mp_mod,       only: gid, domain
  use mpp_domains_mod, only: mpp_update_domains

  implicit none
  real, parameter:: r3 = 1./3., r23 = 2./3.

  private

  public compute_total_energy, Lagrangian_to_Eulerian

CONTAINS

 subroutine Lagrangian_to_Eulerian(consv, ps, pe, delp, pkz, pk,   &
                      mdt, km, is,ie,js,je, isd,ied,jsd,jed,       &
                      nq, u, v, w, delz, pt, q, hs, grav, r_vir, cp,  &
                      akap, kord_mt, kord_tr, kord_tm,  peln, te0_2d,        &
                      ng, ua, va, omga, te, pem, fill, reproduce_sum,        &
                      ak, bk, ks, ze0, hydrostatic, hybrid_z, ktop)
! !INPUT PARAMETERS:
  real,    intent(in):: mdt                   ! mapping time step (same as phys)
  integer, intent(in):: km
  integer, intent(in):: nq                    ! number of tracers (including h2o)
  integer, intent(in):: ng
  integer, intent(in):: is,ie,isd,ied         ! starting & ending X-Dir index
  integer, intent(in):: js,je,jsd,jed         ! starting & ending Y-Dir index
  integer, intent(in):: ks, ktop
  integer, intent(in):: kord_mt               ! Mapping oder for the vector winds
  integer, intent(in):: kord_tr               ! Mapping oder for tracers
  integer, intent(in):: kord_tm               ! Mapping oder for thermodynamics

  real, intent(in):: consv                 ! factor for TE conservation
  real, intent(in):: r_vir
  real, intent(in):: cp, grav
  real, intent(in):: akap
  real, intent(in):: hs(isd:ied,jsd:jed)  ! surface geopotential
  real, intent(in):: te0_2d(is:ie,js:je)

  logical, intent(in):: fill                  ! fill negative tracers
  logical, intent(in):: reproduce_sum
  real, intent(in) :: ak(km+1)
  real, intent(in) :: bk(km+1)

! !INPUT/OUTPUT
  real, intent(inout):: pk(is:ie,js:je,km+1) ! pe to the kappa
  real, intent(inout):: q(isd:ied,jsd:jed,km,*)
  real, intent(inout):: delp(isd:ied,jsd:jed,km) ! pressure thickness
  real, intent(inout)::  pe(is-1:ie+1,km+1,js-1:je+1) ! pressure at layer edges
  real, intent(inout):: pem(is-1:ie+1,km+1,js-1:je+1)
  real, intent(inout):: ps(isd:ied,jsd:jed)      ! surface pressure
  real, intent(inout):: ze0(is:ie,js:je,km+1)    ! Specified height at edges (m)

! u-wind will be ghosted one latitude to the north upon exit
  real, intent(inout)::  u(isd:ied  ,jsd:jed+1,km)   ! u-wind (m/s)
  real, intent(inout)::  v(isd:ied+1,jsd:jed  ,km)   ! v-wind (m/s)
  real, intent(inout)::  w(isd:ied  ,jsd:jed  ,km)   ! vertical velocity (m/s)
  real, intent(inout):: pt(isd:ied  ,jsd:jed  ,km)   ! cp*virtual potential temperature 
                                                     ! as input; output: temperature
  real, intent(inout)::  delz(is:ie,js:je,km)   ! delta-height (m)
  logical, intent(in)::           hydrostatic
  logical, intent(in):: hybrid_z

  real, intent(inout)::   ua(isd:ied,jsd:jed,km)   ! u-wind (m/s) on physics grid
  real, intent(inout)::   va(isd:ied,jsd:jed,km)   ! v-wind (m/s) on physics grid
  real, intent(inout):: omga(isd:ied,jsd:jed,km)   ! vertical press. velocity (pascal/sec)
  real, intent(out)::   peln(is:ie,km+1,js:je)     ! log(pe)
  real, intent(out)::    pkz(is:ie,js:je,km)       ! layer-mean pk for converting t to pt
  real, intent(out)::     te(is:ie,js:je,km)

! !DESCRIPTION:
!
! !REVISION HISTORY:
! SJL 03.11.04: Initial version for partial remapping
!
!-----------------------------------------------------------------------
  integer :: i,j,k 
      real te_2d(is:ie,js:je)
      real  zsum(is:ie,js:je)
      real   q2(is:ie,km)
      real  dp2(is:ie,km)
      real  pe1(is:ie,km+1)
      real  pe2(is:ie,km+1)
      real  pk1(is:ie,km+1)
      real  pk2(is:ie,km+1)
      real  pe0(is:ie+1,km+1)
      real  pe3(is:ie+1,km+1)
      real phis(is:ie,km+1)
      real    gz(is:ie)
! for nonhydrostatic option with hybrid_z coordinate
      real ze1(is:ie,km+1), ze2(is:ie,km+1), deng(is:ie,km)

      real rcp, rg, ak1, tmp, tpe
      real bkh
      real dtmp, pk0
      real dlnp
      real*4 dtmp4
      integer iq, n, kp, k_next
      logical:: diag = .false.
      logical te_map
      real k1k, kapag
  
      k1k   =  akap / (1.-akap)
      kapag = -akap / grav
       rg = akap * cp
      rcp = 1./ cp
      ak1 = (akap + 1.) / akap

      if ( kord_tm < 0 ) then
           te_map = .false.
           do j=js,je
              do k=2,km+1
                 do i=is,ie
                    peln(i,k,j) = log(pe(i,k,j))
                 enddo
              enddo
           enddo
      else
           te_map = .true.
           call pkez(km, is, ie, js, je, pe, pk, akap, peln, pkz)
           call cubed_to_latlon(u, v, ua, va, dx, dy, dxa, dya, km)
! Compute cp*T + KE
!$omp parallel do private(i, j, k)
           do k=1,km
              do j=js,je
                 do i=is,ie
                    te(i,j,k) = 0.5*(ua(i,j,k)**2 + va(i,j,k)**2)  &
                                  +  pt(i,j,k)*pkz(i,j,k)
                 enddo
              enddo
           enddo
     endif

     if ( (.not.hydrostatic) .and. (.not.hybrid_z) ) then
           do k=1,km
              do j=js,je
                 do i=is,ie
                    delz(i,j,k) = -delz(i,j,k) / delp(i,j,k) ! ="specific volume"/grav
                 enddo
              enddo
           enddo
     endif

!$omp parallel do private(i, j, k, pe0, pe1, pe2, pk1, pk2, ......)
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

        do k=2,km+1
           do i=is,ie
              ze2(i,k) = ze0(i,j,k)   ! specified height
           enddo
        enddo
!
! Copy ztop; the top layer must be thick enough to prevent numerical problems.
!
        do i=is,ie
           ze2(i,  1) = ze1(i,1)
           ze0(i,j,1) = ze1(i,1)      ! Note: ze0 updated
        enddo

        do k=1,km
           do i=is,ie
              deng(i,k) = -delp(i,j,k)/delz(i,j,k)  ! density * grav
           enddo
        enddo

        call remap_z(km, ze1, deng, km, ze2, deng, is, ie, abs(kord_tm))
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
               dp2(i,k  ) = -deng(i,k) * delz(i,j,k)
               pe2(i,k+1) =   pe2(i,k) +  dp2(i,k)
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
       if(nq /= 0) then
!------------------------------------------------------------------
! Do remapping one tracer at a time; seems to be faster on the SGI
! It requires less memory than mapn_ppm
!------------------------------------------------------------------
          do iq=1,nq
            call map1_q2(km, pe1, q(isd,jsd,1,iq),     &
                         km, pe2, q2, dp2,             &
                         is, ie, 0, kord_tr, j, isd, ied, jsd, jed)
            if (fill) call fillz(ie-is+1, km, 1, q2, dp2)
            do k=1,km
               do i=is,ie
                  q(i,j,k,iq) = q2(i,k)
               enddo
            enddo
          enddo
       endif

       do k=1,km+1
          do i=is,ie
             pk1(i,k) = pk(i,j,k)
          enddo
       enddo
       do i=is,ie
          pk2(i,   1) = pk1(i,   1)
          pk2(i,km+1) = pk1(i,km+1)
       enddo
       do k=2,km
          do i=is,ie
             pk2(i,k) = pe2(i,k) ** akap
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
!----------------
! Map pt using pk
!----------------
        call map1_ppm (km,  pk1,  pt,           &
                       km,  pk2,  pt,           &
                       is, ie, j, isd, ied, jsd, jed, 1, abs(kord_tm))
   endif

   if ( .not. hydrostatic ) then
! Remap vertical wind:
        call map1_ppm (km,   pe1,  w,       &
                       km,   pe2,  w,       &
                       is, ie, j, isd, ied, jsd, jed, -1, kord_mt)
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

! Copy omega field to pe3
   do i=is,ie
      pe3(i,1) = 0.
   enddo
   do k=2,km+1
      do i=is,ie
         pe3(i,k) = omga(i,j,k-1)
      enddo
   enddo
   do k=1,km+1
      do i=is,ie
         pe0(i,k) = peln(i,k,j)
      enddo
   enddo

!--------------
! Compute peln
!--------------

   if ( hybrid_z ) then
      do k=2,km+1
         do i=is,ie
            peln(i,k,j) = log(pe2(i,k))  ! peln is used anywhere?
         enddo
      enddo
   else
      do k=2,ks+1
         tmp = log(ak(k))
         do i=is,ie
            peln(i,k,j) = tmp
         enddo
      enddo
      do k=ks+2,km+1
         do i=is,ie
            peln(i,k,j) = log(pe2(i,k))
         enddo
      enddo

      if( ptop < ptop_min ) then
          do i=is,ie
             peln(i,1,j) = peln(i,2,j) - ak1
          enddo
      else
          tmp = log( ptop )
          do i=is,ie
             peln(i,1,j) = tmp
          enddo
      endif
   endif

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
      do k=ktop,km
         do i=is,ie
! Note: pt at this stage is cp*Theta_v
            pkz(i,j,k) = ( kapag*delp(i,j,k)*pt(i,j,k) /            &
                          (delz(i,j,k)*(1.+r_vir*q(i,j,k,1))) )**k1k
         enddo
      enddo
   endif

! Interpolate omega/pe3 (defined at pe0) to remapped cell center (dp2)
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

if ( hybrid_z ) then   !------- Hybrid_z section ---------------
     call mpp_update_domains(ua , domain,  whalo=1, ehalo=1,     &
                             shalo=1, nhalo=1, complete=.true.)
! u-wind
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
endif         !------------- Hybrid_z section ----------------------

!$omp parallel do private(i,j,k)
     do k=2,km
        do j=js,je
           do i=is,ie
              pe(i,k,j) = ua(i,j,k-1)
           enddo
        enddo
     enddo

  call cubed_to_latlon(u,  v, ua, va, dx, dy, dxa, dya, km)

  if( consv > 0. ) then

    if ( te_map ) then
!$omp parallel do private(i, j, k)
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
!$omp parallel do private(i,j,k,gz)
      do j=js,je
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
                                      0.5*(ua(i,j,k)**2+va(i,j,k)**2))
            enddo
         enddo
      enddo
    endif

!$omp parallel do private(i, j, k)
      do j=js,je
         do i=is,ie
            te_2d(i,j) = te0_2d(i,j) - te_2d(i,j)
             zsum(i,j) = ptop*(pk(i,j,1)-pk(i,j,km+1))
         enddo

         do k=1,km
            do i=is,ie
               zsum(i,j) = zsum(i,j) + pkz(i,j,k)*delp(i,j,k)
            enddo
         enddo
      enddo

      tpe = g_sum(te_2d, is,ie, js,je, ng, area)
      dtmp = consv * tpe /(cp*g_sum(zsum, is,ie, js,je, ng, area))
!-------------------------------------------------------------------------------
! One may use this quick fix to ensure reproducibility at the expense of a lower
! floating precision; this is fine for the TE correction
!-------------------------------------------------------------------------------
      if ( reproduce_sum ) then
           dtmp4 = dtmp
           dtmp  = dtmp4
      endif

      if( diag ) then
          pk0 = 1.E5 ** akap 
          write(6,*) 'heating (deg/day) =',pk0*dtmp*86400./mdt
      endif
  else
      dtmp = 0.
  endif        ! end consv check

  if ( te_map ) then
!$omp parallel do private(i, j, k, gz, tpe, tmp, dlnp)
      do j=js,je
         do i=is,ie
            gz(i) = hs(i,j)
         enddo
         do k=km,1,-1
            do i=is,ie
               tpe = te(i,j,k) - 0.5*(ua(i,j,k)**2 + va(i,j,k)**2) - gz(i)
               dlnp = rg*(peln(i,k+1,j) - peln(i,k,j))
               tmp = tpe/((cp - pe(i,k,j)*dlnp/delp(i,j,k))*(1.+r_vir*q(i,j,k,1)) )
               pt(i,j,k) =  tmp + dtmp*pkz(i,j,k)/(1.+r_vir*q(i,j,k,1))
               gz(i) = gz(i) + dlnp*tmp*(1.+r_vir*q(i,j,k,1))
            enddo
         enddo           ! end k-loop
      enddo
  else
!$omp parallel do private(i, j, k)
      do k=1,km
         do j=js,je
            do i=is,ie
               pt(i,j,k) = (rcp*pt(i,j,k) + dtmp)*pkz(i,j,k)/(1.+r_vir*q(i,j,k,1))
            enddo
         enddo   
      enddo
  endif

 end subroutine Lagrangian_to_Eulerian


 subroutine compute_total_energy(is, ie, js, je, isd, ied, jsd, jed, km,  &
                                 u, v, pt, delp, q, pe, peln, hs,         &
                                 r_vir,  cp, rg, hlv, te_2d, ua, va, teq, &
                                 moist_phys, id_te)
!------------------------------------------------------
! Compute vertically integrated total energy per column
!------------------------------------------------------
! !INPUT PARAMETERS:
   integer,  intent(in):: km, is, ie, js, je, isd, ied, jsd, jed, id_te
   real, intent(in), dimension(isd:ied,jsd:jed,km):: pt, delp, q
   real, intent(in)::  u(isd:ied,  jsd:jed+1,km)
   real, intent(in)::  v(isd:ied+1,jsd:jed,  km)
   real, intent(in):: hs(isd:ied,jsd:jed)  ! surface geopotential
   real, intent(in)::   pe(is-1:ie+1,km+1,js-1:je+1) ! pressure at layer edges
   real, intent(in):: peln(is:ie,km+1,js:je)  ! log(pe)
   real, intent(in):: cp, rg, r_vir, hlv
   logical, intent(in):: moist_phys
! Output:
   real, intent(out), dimension(isd:ied,jsd:jed,km):: ua, va
   real, intent(out):: te_2d(is:ie,js:je)   ! vertically integrated TE
   real, intent(out)::   teq(is:ie,js:je)   ! Moist TE
! Local
   real  gztop(is:ie)
   integer i, j, k

!----------------------
! Output lat-lon winds:
!----------------------
  call cubed_to_latlon(u, v, ua, va, dx, dy, dxa, dya, km)

!$omp parallel do private(i,j,k,gztop)
  do j=js,je
     do i=is,ie
        gztop(i) = hs(i,j)
        do k=1,km
           gztop(i) = gztop(i) + (peln(i,k+1,j)-peln(i,k,j)) *   &
                      rg*pt(i,j,k)*(1.+r_vir*q(i,j,k))
        enddo
     enddo
     do i=is,ie
        te_2d(i,j) = pe(i,km+1,j)*hs(i,j) - pe(i,1,j)*gztop(i)
     enddo

     do k=1,km
        do i=is,ie
           te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(0.5*(ua(i,j,k)**2+va(i,j,k)**2)  &
                                   + cp*pt(i,j,k)*(1.+r_vir*q(i,j,k)))
        enddo
     enddo
  enddo

!-------------------------------------
! Doganostics computation for moist TE
!-------------------------------------
  if( id_te>0 ) then
      do j=js,je
         do i=is,ie
            teq(i,j) = te_2d(i,j)
         enddo
      enddo
      if ( moist_phys ) then
           do k=1,km
              do j=js,je
                 do i=is,ie
                    teq(i,j) = teq(i,j) + hlv*q(i,j,k)*delp(i,j,k)
                 enddo
              enddo
           enddo
      endif
      do j=js,je
         do i=is,ie
            teq(i,j) = teq(i,j) / (pe(i,km,j) - pe(i,1,j))
         enddo
      enddo
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
   real, intent(out):: peln(ifirst:ilast, km+1, jfirst:jlast)   ! log (pe)
! Local
   real pk2(ifirst:ilast, km+1)
   real pek
   real lnp
   real ak1
   integer i, j, k

   ak1 = (akap + 1.) / akap

!$omp  parallel do default(shared) private(ixj, i1, i2, i, j, k, pek, lnp, pk2)
   do j=jfirst, jlast
        pek = pk(ifirst,j,1)
        do i=ifirst, ilast
           pk2(i,1) = pek
        enddo

        do k=2,km+1
           do i=ifirst, ilast
              peln(i,k,j) =  log(pe(i,k,j))
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



 subroutine remap_z(km, pe1, q1, kn, pe2, q2, i1, i2, kord)

! !INPUT PARAMETERS:
      integer, intent(in) :: i1                ! Starting longitude
      integer, intent(in) :: i2                ! Finishing longitude
      integer, intent(in) :: kord              ! Method order
      integer, intent(in) :: km                ! Original vertical dimension
      integer, intent(in) :: kn                ! Target vertical dimension

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

      integer i, k, l, ll, k0

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)      ! negative
            q4(1,i,k) = q1(i,k)
         enddo
      enddo

! Compute vertical subgrid distribution
      call ppm2m( q4, dp1, km, i1, i2, 0, kord )

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
              do ll=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if(pe2(i,k+1) < pe1(i,ll+1) ) then
! Whole layer..
                    qsum = qsum + dp1(i,ll)*q4(1,i,ll)
                 else
                    delp = pe2(i,k+1)-pe1(i,ll)
                    esl = delp / dp1(i,ll)
                    qsum = qsum + delp*(q4(2,i,ll)+0.5*esl*               &
                         (q4(3,i,ll)-q4(2,i,ll)+q4(4,i,ll)*(1.-r23*esl)))
                    k0 = ll
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
                    
! !INPUT PARAMETERS:
      integer, intent(in) :: i1                ! Starting longitude
      integer, intent(in) :: i2                ! Finishing longitude
      integer, intent(in) :: iv                ! Mode: 0 ==  constituents  1 == ???
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
!
!     Perform piecewise parabolic method on a given latitude    
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
!
! !LOCAL VARIABLES:
      real    dp1(i1:i2,km)
      real   q4(4,i1:i2,km)
      real    pl, pr, qsum, dp, esl

      integer i, k, l, ll, k0

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
            q4(1,i,k) = q1(i,j,k)
         enddo
      enddo

! Compute vertical subgrid distribution
      call ppm2m( q4, dp1, km, i1, i2, iv, kord )

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
            q2(i,j,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
          else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do ll=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if(pe2(i,k+1) > pe1(i,ll+1) ) then
! Whole layer..
                 qsum = qsum + dp1(i,ll)*q4(1,i,ll)
                 else
                   dp = pe2(i,k+1)-pe1(i,ll)
                  esl = dp / dp1(i,ll)
                 qsum = qsum + dp*(q4(2,i,ll)+0.5*esl*               &
                       (q4(3,i,ll)-q4(2,i,ll)+q4(4,i,ll)*(1.-r23*esl)))
                 k0 = ll
                 goto 123
                 endif
              enddo
              goto 123
           endif
      endif
100   continue
123   q2(i,j,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
1000  continue

 end subroutine map1_ppm


subroutine map1_q2( km,   pe1,   q1,            &
                    kn,   pe2,   q2,   dp2,     &
                    i1,   i2,    iv,   kord, j, &
                    ibeg, iend, jbeg, jend )


! !INPUT PARAMETERS:
      integer, intent(in) :: j                 ! Current latitude
      integer, intent(in) :: i1                ! Starting longitude
      integer, intent(in) :: i2                ! Finishing longitude
      integer, intent(in) :: ibeg, iend, jbeg, jend
      integer, intent(in) :: iv                ! Mode: 0 ==  constituents  1 == ???
      integer, intent(in) :: kord              ! Method order
      integer, intent(in) :: km                ! Original vertical dimension
      integer, intent(in) :: kn                ! Target vertical dimension

      real, intent(in) ::  pe1(i1:i2,km+1)  ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the original vertical coordinate
      real, intent(in) ::  pe2(i1:i2,kn+1)  ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the new vertical coordinate
      real, intent(in) ::  q1(ibeg:iend,jbeg:jend,km) ! Field input
      real, intent(in) ::  dp2(i1:i2,kn)
! !INPUT/OUTPUT PARAMETERS:
      real, intent(inout):: q2(i2-i1+1,kn) ! Field output

! !DESCRIPTION:
!
!     Perform piecewise parabolic method on a given latitude    
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
!
! !REVISION HISTORY: 
!    00.04.24   Lin       Last modification
!    01.03.26   Sawyer    Added ProTeX documentation
!    02.04.04   Sawyer    incorporated latest FVGCM version
!    02.06.20   Sawyer    made Q2 inout since the args for Q1/Q2 same
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      real   dp1(i1:i2,km)
      real   q4(4,i1:i2,km)
      real   pl, pr, qsum, dp, esl

      integer i, k, l, ll, k0

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
            q4(1,i,k) = q1(i,j,k)
         enddo
      enddo

! Compute vertical subgrid distribution
      call ppm2m( q4, dp1, km, i1, i2, iv, kord )

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
            q2(i-i1+1,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
          else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do ll=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if(pe2(i,k+1) > pe1(i,ll+1) ) then
! Whole layer..
                 qsum = qsum + dp1(i,ll)*q4(1,i,ll)
                 else
                   dp = pe2(i,k+1)-pe1(i,ll)
                  esl = dp / dp1(i,ll)
                 qsum = qsum + dp*(q4(2,i,ll)+0.5*esl*               &
                       (q4(3,i,ll)-q4(2,i,ll)+q4(4,i,ll)*(1.-r23*esl)))
                 k0 = ll
                 goto 123
                 endif
              enddo
              goto 123
           endif
      endif
100   continue
123   q2(i-i1+1,k) = qsum / dp2(i,k)
555   continue
1000  continue

!EOC
 end subroutine map1_q2

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  ppm2m --- Piecewise parabolic method for fields
!
! !INTERFACE:
 subroutine ppm2m(a4, delp, km, i1, i2, iv, kord)

! !INPUT PARAMETERS:
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
 integer, intent(in):: i1      ! Starting longitude
 integer, intent(in):: i2      ! Finishing longitude
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: kord    ! Order (or more accurately method no.):
                               ! 
 real , intent(in):: delp(i1:i2,km)     ! layer pressure thickness

! !INPUT/OUTPUT PARAMETERS:
 real , intent(inout):: a4(4,i1:i2,km)  ! Interpolated values

! !DESCRIPTION:
!
!   Perform the piecewise parabolic method 
! 
! !REVISION HISTORY: 
!   ??.??.??    Lin        Creation
!   02.04.04    Sawyer     Newest release from FVGCM
! 
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! local arrays:
      real    dc(i1:i2,km)
      real    h2(i1:i2,km)
      real  delq(i1:i2,km)
      real   df2(i1:i2,km)
      real    d4(i1:i2,km)

! local scalars:
      integer i, k, km1, lmt
      integer it
      real  fac
      real  a1, a2, c1, c2, c3, d1, d2
      real  qmax, qmin, cmax, cmin
      real  qm, dq, tmp
      real  qmp, pmp
      real  lac

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
            tmp = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /      &
                                    (d4(i,k)+delp(i,k+1))
            qmax = max(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1)) - a4(1,i,k)
            qmin = a4(1,i,k) - min(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))
             dc(i,k) = sign(min(abs(tmp),qmax,qmin), tmp)
            df2(i,k) = tmp
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

      if(km>8 .and. kord>3) call steepz(i1, i2, km, a4, df2, dc, delq, delp, d4)

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
         a4(2,i,1) = d1*(2.*c1*d1**2-c3) + a4(2,i,2)
!        dc(i,1) =       a4(1,i,1) - a4(2,i,1)
         dc(i,1) =  0.5*(a4(1,i,1) - a4(2,i,1))
! No over- and undershoot condition
         cmax = max(a4(1,i,1), a4(1,i,2))
         cmin = min(a4(1,i,1), a4(1,i,2))
         a4(2,i,2) = max(cmin,a4(2,i,2))
         a4(2,i,2) = min(cmax,a4(2,i,2))
      enddo

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
         a4(3,i,km) = d1*(8.*c1*d1**2-c3) + a4(2,i,km)
!        dc(i,km) =      a4(3,i,km) - a4(1,i,km)
         dc(i,km) = 0.5*(a4(3,i,km) - a4(1,i,km))
! No over- and under-shoot condition
         cmax = max(a4(1,i,km), a4(1,i,km1))
         cmin = min(a4(1,i,km), a4(1,i,km1))
         a4(2,i,km) = max(cmin,a4(2,i,km))
         a4(2,i,km) = min(cmax,a4(2,i,km))
      enddo

      do k=1,km1
         do i=i1,i2
            a4(3,i,k) = a4(2,i,k+1)
         enddo
      enddo

! Enforce monotonicity of the "slope" within the top layer
      do i=i1,i2
         if ( a4(2,i,1) * a4(1,i,1) <= 0. ) then
              a4(2,i,1) = 0.
!               dc(i,1) = a4(1,i,1)
         endif
         if ( dc(i,1) * (a4(2,i,2) - a4(1,i,1)) <= 0. ) then
! Setting DC==0 will force piecewise constant distribution after
! calling kmppm
              dc(i,1) = 0.
         endif
      enddo

! Enforce constraint on the "slope" at the surface

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

!-----------------------------------------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
!-----------------------------------------------------------
! Top 2 and bottom 2 layers always use monotonic mapping
      do k=1,2
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
            call kmppm(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

      if(kord >= 7) then
!-----------------------
! Huynh's 2nd constraint
!-----------------------
      do k=2,km1
         do i=i1,i2
! Method#1
!           h2(i,k) = delq(i,k) - delq(i,k-1)
! Method#2
!           h2(i,k) = 2.*(dc(i,k+1)/delp(i,k+1) - dc(i,k-1)/delp(i,k-1))
!    &               / ( delp(i,k)+0.5*(delp(i,k-1)+delp(i,k+1)) )
!    &               * delp(i,k)**2
! Method#3
            h2(i,k) = dc(i,k+1) - dc(i,k-1)
         enddo
      enddo

      if( kord == 7 ) then
         fac = 1.5           ! original quasi-monotone
      else
         fac = 0.125         ! full monotone
      endif

      do k=3,km-2
        do i=i1,i2
! Right edges
!        qmp   = a4(1,i,k) + 2.0*delq(i,k-1)
!        lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
!
         pmp   = 2.*dc(i,k)
         qmp   = a4(1,i,k) + pmp
         lac   = a4(1,i,k) + fac*h2(i,k-1) + dc(i,k)
         qmin  = min(a4(1,i,k), qmp, lac)
         qmax  = max(a4(1,i,k), qmp, lac)
         a4(3,i,k) = min(max(a4(3,i,k), qmin), qmax)
! Left  edges
!        qmp   = a4(1,i,k) - 2.0*delq(i,k)
!        lac   = a4(1,i,k) + fac*h2(i,k+1) - 0.5*delq(i,k)
!
         qmp   = a4(1,i,k) - pmp
         lac   = a4(1,i,k) + fac*h2(i,k+1) - dc(i,k)
         qmin  = min(a4(1,i,k), qmp, lac)
         qmax  = max(a4(1,i,k), qmp, lac)
         a4(2,i,k) = min(max(a4(2,i,k), qmin), qmax)
! Recompute A6
         a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
! Additional constraint to ensure positivity when kord=7
         if (iv == 0 .and. kord == 7) then
             call kmppm(dc(i1,k), a4(1,i1,k), it, 2)
         endif
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
         call kmppm(dc(i1,k), a4(1,i1,k), it, lmt)
      enddo
      endif

      do i=i1,i2
         a4(4,i,km1) = 3.*(2.*a4(1,i,km1) - (a4(2,i,km1)+a4(3,i,km1)))
      enddo
      call kmppm(dc(i1,km1), a4(1,i1,km1), it, 0)
!EOC
 end subroutine ppm2m
!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  kmppm --- Perform piecewise parabolic method in vertical
!
! !INTERFACE:
 subroutine kmppm(dm, a4, itot, lmt)

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

      real , parameter:: r12 = 1./12.
      real  qmp
      real  da1, da2, a6da
      real  fmin
      integer i

! Developer: S.-J. Lin, NASA-GSFC
! Last modified: Apr 24, 2000

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

 end subroutine kmppm
!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  steepz --- Calculate attributes for PPM
!
! !INTERFACE:
 subroutine steepz(i1, i2, km, a4, df2, dm, dq, dp, d4)

! !INPUT PARAMETERS:
      integer, intent(in) :: km                   ! Total levels
      integer, intent(in) :: i1                   ! Starting longitude
      integer, intent(in) :: i2                   ! Finishing longitude
      real , intent(in) ::  dp(i1:i2,km)       ! grid size
      real , intent(in) ::  dq(i1:i2,km)       ! backward diff of q
      real , intent(in) ::  d4(i1:i2,km)       ! backward sum:  dp(k)+ dp(k-1) 
      real , intent(in) :: df2(i1:i2,km)       ! first guess mismatch
      real , intent(in) ::  dm(i1:i2,km)       ! monotonic mismatch

! !INPUT/OUTPUT PARAMETERS:
      real , intent(inout) ::  a4(4,i1:i2,km)  ! first guess/steepened
! !DESCRIPTION:
!   This is complicated stuff related to the Piecewise Parabolic Method
!   and I need to read the Collela/Woodward paper before documenting
!   thoroughly.
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

end module mapz_module

