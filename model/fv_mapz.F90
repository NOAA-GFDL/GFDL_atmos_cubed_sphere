!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be
!* useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

! SJL: Apr 12, 2012
! This revision may actually produce rounding level differences due to the elimination of KS to compute
! pressure level for remapping.
! Linjiong Zhou: Nov 19, 2019
! Revise the OpenMP code to avoid crash
module fv_mapz_mod

  use constants_mod,     only: pi=>pi_8, rvgas, rdgas, grav, hlv, hlf, cp_air, cp_vapor
  use fv_arrays_mod,     only: radius ! scaled for small earth
  use tracer_manager_mod,only: get_tracer_index, adjust_mass
  use field_manager_mod, only: MODEL_ATMOS
  use fv_grid_utils_mod, only: g_sum, ptop_min, cubed_to_latlon
  use fv_fill_mod,       only: fillz
  use mpp_domains_mod,   only: mpp_update_domains, domain2d
  use mpp_mod,           only: FATAL, NOTE, mpp_error, get_unit, mpp_root_pe, mpp_pe
  use fv_arrays_mod,     only: fv_grid_type, fv_grid_bounds_type, R_GRID, inline_mp_type
  use fv_timing_mod,     only: timing_on, timing_off
  use fv_mp_mod,         only: is_master, mp_reduce_min, mp_reduce_max
  use intermediate_phys_mod, only: intermediate_phys
  use gfdl_mp_mod,       only: c_liq, c_ice

  implicit none
  real, parameter:: consv_min = 0.001   ! below which no correction applies
  real, parameter:: t_min= 184.   ! below which applies stricter constraint
  real, parameter:: r3 = 1./3., r23 = 2./3., r12 = 1./12.
  real, parameter:: cv_vap = 3.*rvgas  ! 1384.5
  real, parameter:: cv_air =  cp_air - rdgas ! = rdgas * (7/2-1) = 2.5*rdgas=717.68
  real, parameter:: cp_vap = cp_vapor   ! 1846.
  real, parameter:: tice = 273.16

  real, parameter :: w_max = 90.
  real, parameter :: w_min = -60.

  real(kind=4) :: E_Flux = 0.
  private

  public compute_total_energy, Lagrangian_to_Eulerian, moist_cv, moist_cp,   &
         rst_remap, mappm, E_Flux, remap_2d, map_scalar, consv_min, map1_q2

contains

 subroutine Lagrangian_to_Eulerian(last_step, consv, ps, pe, delp, pkz, pk,   &
                                   mdt, pdt, npx, npy, km, is,ie,js,je, isd,ied,jsd,jed,       &
                      nq, nwat, sphum, q_con, u, v, w, delz, pt, q, hs, r_vir, cp, te_err, tw_err, &
                      akap, cappa, kord_mt, kord_wz, kord_tr, kord_tm, remap_te, peln, te0_2d,        &
                      ng, ua, va, omga, te, ws, fill, reproduce_sum,      &
                      ptop, ak, bk, pfull, gridstruct, domain, do_sat_adj, &
                      hydrostatic, hybrid_z, adiabatic, do_adiabatic_init, &
                      do_inline_mp, inline_mp, c2l_ord, bd, fv_debug, &
                      w_limiter, do_fast_phys, do_intermediate_phys, consv_checker, adj_mass_vmr)

  logical, intent(in):: last_step
  logical, intent(in):: fv_debug
  logical, intent(in):: w_limiter
  logical, intent(in):: do_fast_phys
  logical, intent(in):: do_intermediate_phys
  logical, intent(in):: consv_checker
  integer, intent(in):: adj_mass_vmr
  real,    intent(in):: mdt                   ! remap time step
  real,    intent(in):: pdt                   ! phys time step
  integer, intent(in):: npx, npy
  integer, intent(in):: km
  integer, intent(in):: nq                    ! number of tracers (including h2o)
  integer, intent(in):: nwat
  integer, intent(in):: sphum                 ! index for water vapor (specific humidity)
  integer, intent(in):: ng
  integer, intent(in):: is,ie,isd,ied         ! starting & ending X-Dir index
  integer, intent(in):: js,je,jsd,jed         ! starting & ending Y-Dir index
  integer, intent(in):: kord_mt               ! Mapping order for the vector winds
  integer, intent(in):: kord_wz               ! Mapping order/option for w
  integer, intent(in):: kord_tr(nq)           ! Mapping order for tracers
  integer, intent(in):: kord_tm               ! Mapping order for thermodynamics
  integer, intent(in):: c2l_ord

  real, intent(in):: consv                 ! factor for TE conservation
  real, intent(in):: r_vir
  real, intent(in):: cp
  real, intent(in):: te_err
  real, intent(in):: tw_err
  real, intent(in):: akap
  real, intent(in):: hs(isd:ied,jsd:jed)  ! surface geopotential
  real, intent(inout):: te0_2d(is:ie,js:je)
  real, intent(in):: ws(is:ie,js:je)

  logical, intent(in):: do_sat_adj
  logical, intent(in):: do_inline_mp
  logical, intent(in):: fill                  ! fill negative tracers
  logical, intent(in):: reproduce_sum
  logical, intent(in):: adiabatic, do_adiabatic_init
  logical, intent(in):: remap_te
  real, intent(in) :: ptop
  real, intent(in) :: ak(km+1)
  real, intent(in) :: bk(km+1)
  real, intent(in):: pfull(km)
  type(fv_grid_type), intent(IN), target :: gridstruct
  type(domain2d), intent(INOUT) :: domain
  type(fv_grid_bounds_type), intent(IN) :: bd

! !INPUT/OUTPUT
  real, intent(inout):: pk(is:ie,js:je,km+1) ! pe to the kappa
  real, intent(inout):: q(isd:ied,jsd:jed,km,*)
  real, intent(inout):: delp(isd:ied,jsd:jed,km) ! pressure thickness
  real, intent(inout)::  pe(is-1:ie+1,km+1,js-1:je+1) ! pressure at layer edges
  real, intent(inout):: ps(isd:ied,jsd:jed)      ! surface pressure

! u-wind will be ghosted one latitude to the north upon exit
  real, intent(inout)::  u(isd:ied  ,jsd:jed+1,km)   ! u-wind (m/s)
  real, intent(inout)::  v(isd:ied+1,jsd:jed  ,km)   ! v-wind (m/s)
  real, intent(inout)::  w(isd:     ,jsd:     ,1:)   ! vertical velocity (m/s)
  real, intent(inout):: pt(isd:ied  ,jsd:jed  ,km)   ! cp*virtual potential temperature
                                                     ! as input; output: temperature
  real, intent(inout), dimension(isd:,jsd:,1:)::q_con, cappa
  real, intent(inout), dimension(is:,js:,1:)::delz
  logical, intent(in):: hydrostatic
  logical, intent(in):: hybrid_z

  real, intent(inout)::   ua(isd:ied,jsd:jed,km)   ! u-wind (m/s) on physics grid
  real, intent(inout)::   va(isd:ied,jsd:jed,km)   ! v-wind (m/s) on physics grid
  real, intent(inout):: omga(isd:ied,jsd:jed,km)   ! vertical press. velocity (pascal/sec)
  real, intent(inout)::   peln(is:ie,km+1,js:je)     ! log(pe)
  real, intent(out)::    pkz(is:ie,js:je,km)       ! layer-mean pk for converting t to pt
  real, intent(out)::     te(isd:ied,jsd:jed,km)

  type(inline_mp_type), intent(inout):: inline_mp

! !DESCRIPTION:
!
! !REVISION HISTORY:
! SJL 03.11.04: Initial version for partial remapping
!
!-----------------------------------------------------------------------
  real, dimension(is:ie,js:je):: te_2d, zsum0, zsum1
  real, dimension(is:ie,km)  :: q2, dp2, t0, w2
  real, dimension(is:ie,km+1):: pe1, pe2, pk1, pk2, pn2, phis
  real, dimension(isd:ied,jsd:jed,km):: pe4
  real, dimension(is:ie+1,km+1):: pe0, pe3
  real, dimension(is:ie+1):: gz, cvm
  real, dimension(isd:ied,jsd:jed,km):: qnl, qni

  real rcp, rg, rrg, bkh, dtmp, k1k, dlnp, tpe
  integer:: i,j,k
  integer:: nt, liq_wat, ice_wat, rainwat, snowwat, cld_amt, graupel, w_diff, iq, n, kmp, kp, k_next
  integer:: ccn_cm3, cin_cm3, aerosol

  k1k = rdgas/cv_air   ! akap / (1.-akap) = rg/Cv=0.4
  rg = rdgas
  rcp = 1./ cp
  rrg = -rdgas/grav

  liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
  ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
  rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat')
  snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat')
  graupel = get_tracer_index (MODEL_ATMOS, 'graupel')
  cld_amt = get_tracer_index (MODEL_ATMOS, 'cld_amt')
  w_diff  = get_tracer_index (MODEL_ATMOS, 'w_diff')
  ccn_cm3 = get_tracer_index (MODEL_ATMOS, 'ccn_cm3')
  cin_cm3 = get_tracer_index (MODEL_ATMOS, 'cin_cm3')
  aerosol = get_tracer_index (MODEL_ATMOS, 'aerosol')

!$OMP parallel do default(none) shared(is,ie,js,je,km,pe,ptop,kord_tm,hydrostatic, &
!$OMP                                  pt,pk,rg,peln,q,nwat,liq_wat,rainwat,ice_wat,snowwat,    &
!$OMP                                  graupel,q_con,sphum,cappa,r_vir,k1k,delp, &
!$OMP                                  delz,akap,pkz,te,u,v,ps, gridstruct, last_step, &
!$OMP                                  ak,bk,nq,isd,ied,jsd,jed,kord_tr,fill, adiabatic, &
!$OMP                                  hs,w,ws,kord_wz,omga,rrg,kord_mt,pe4,w_limiter,cp,remap_te)    &
!$OMP                          private(gz,cvm,kp,k_next,bkh,dp2,dlnp,tpe,   &
!$OMP                                  pe0,pe1,pe2,pe3,pk1,pk2,pn2,phis,q2,w2)

  do j=js,je+1

     !0) Prepare pressure (pe, pk, ps, delp), temperature, density, and energy

     do k=1,km+1
        do i=is,ie
           pe1(i,k) = pe(i,k,j)
        enddo
     enddo

     do i=is,ie
        pe2(i,   1) = ptop
        pe2(i,km+1) = pe(i,km+1,j)
     enddo

     if ( j /= (je+1) ) then

        if (  .not. remap_te ) then
           if (  kord_tm < 0 ) then
              ! Note: pt at this stage is Theta_v
              if ( hydrostatic ) then
                 ! Transform virtual pt to virtual Temp
                 do k=1,km
                    do i=is,ie
                       pt(i,j,k) = pt(i,j,k)*(pk(i,j,k+1)-pk(i,j,k))/(akap*(peln(i,k+1,j)-peln(i,k,j)))
                    enddo
                 enddo
              else
                 ! Transform "density pt" to "density temp"
                 do k=1,km
#ifdef MOIST_CAPPA
                    call moist_cv(is,ie,isd,ied,jsd,jed, km, j, k, nwat, sphum, liq_wat, rainwat,    &
                         ice_wat, snowwat, graupel, q, gz(is:ie), cvm(is:ie))
                    do i=is,ie
                       q_con(i,j,k) = gz(i)
                       cappa(i,j,k) = rdgas / ( rdgas + cvm(i)/(1.+r_vir*q(i,j,k,sphum)) )
                       pt(i,j,k) = pt(i,j,k)*exp(cappa(i,j,k)/(1.-cappa(i,j,k))*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
                    enddo
#else
                    do i=is,ie
                       pt(i,j,k) = pt(i,j,k)*exp(k1k*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
                       ! Using dry pressure for the definition of the virtual potential temperature
                       !                    pt(i,j,k) = pt(i,j,k)*exp(k1k*log(rrg*(1.-q(i,j,k,sphum))*delp(i,j,k)/delz(i,j,k)*    &
                       !                                              pt(i,j,k)/(1.+r_vir*q(i,j,k,sphum))))
                    enddo
#endif
                 enddo
              endif         ! hydro test


           endif !kord_tm
        else
           !----------------------------------
           ! Compute cp*T + KE +phis
           do i=is,ie
              phis(i,km+1) = hs(i,j)
           enddo
           if ( hydrostatic ) then
              call pkez(km, is, ie, js, je, j, pe, pk, akap, peln, pkz, ptop)
              do k=km,1,-1
                 do i=is,ie
                    phis(i,k) = phis(i,k+1) + cp_air*pt(i,j,k)*(pk(i,j,k+1)-pk(i,j,k))  !Pt:theta_v
                 enddo
              enddo
              do k=1,km+1
                 do i=is,ie
                    phis(i,k) = phis(i,k) * pe1(i,k)
                 enddo
              enddo
              do k=1,km
                 do i=is,ie
                    te(i,j,k) = 0.25*gridstruct%rsin2(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                         v(i,j,k)**2+v(i+1,j,k)**2 -  &
                         (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j))  &
                         + cp_air*pt(i,j,k)*pkz(i,j,k) +  (phis(i,k+1)-phis(i,k))/(pe1(i,k+1)-pe1(i,k))
                 enddo
              enddo
           else
              do k=km,1,-1
                 do i=is,ie
                    phis(i,k) = phis(i,k+1) - grav*delz(i,j,k)
                 enddo
#ifdef MOIST_CAPPA
                 call moist_cv(is,ie,isd,ied,jsd,jed, km, j, k, nwat, sphum, liq_wat, rainwat,    &
                      ice_wat, snowwat, graupel, q, gz(is:ie), cvm(is:ie))
                 do i=is,ie
                    q_con(i,j,k) = gz(i)
                    cappa(i,j,k) = rdgas / ( rdgas + cvm(i)/(1.+r_vir*q(i,j,k,sphum)) )
                    pkz(i,j,k) = exp(cappa(i,j,k)/(1.-cappa(i,j,k))*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
                    te(i,j,k) = cvm(i)*pt(i,j,k)*pkz(i,j,k)/((1.+r_vir*q(i,j,k,sphum))*(1.-gz(i))) +     &
                         0.5 * w(i,j,k)**2 + 0.25*gridstruct%rsin2(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                         v(i,j,k)**2+v(i+1,j,k)**2 -  &
                         (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j)) +         &
                         0.5*(phis(i,k+1)+phis(i,k))
                 enddo
#else
                 do i=is,ie
                    pkz(i,j,k) = exp(k1k*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
                    te(i,j,k) = cv_air*pt(i,j,k)*pkz(i,j,k)/(1.+r_vir*q(i,j,k,sphum)) +     &
                         0.5 * w(i,j,k)**2 + 0.25*gridstruct%rsin2(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                         v(i,j,k)**2+v(i+1,j,k)**2 -  &
                         (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j)) +         &
                         0.5*(phis(i,k+1)+phis(i,k))
                 enddo
#endif
              enddo

           endif !end hydrostatic test
        endif ! .not. remap_te

        if ( .not. hydrostatic ) then
           do k=1,km
              do i=is,ie
                 delz(i,j,k) = -delz(i,j,k) / delp(i,j,k) ! ="specific volume"/grav
              enddo
           enddo
        endif

! update ps
        do i=is,ie
           ps(i,j) = pe1(i,km+1)
        enddo
!
! Hybrid sigma-P coordinate:
!
        do k=2,km
           do i=is,ie
              pe2(i,k) = ak(k) + bk(k)*pe(i,km+1,j)
           enddo
        enddo
        do k=1,km
           do i=is,ie
              dp2(i,k) = pe2(i,k+1) - pe2(i,k)
           enddo
        enddo

!------------
! update delp
!------------
      do k=1,km
         do i=is,ie
            delp(i,j,k) = dp2(i,k)
         enddo
      enddo

!------------------
! Compute p**Kappa
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
            pn2(i,k) = log(pe2(i,k))
            pk2(i,k) = exp(akap*pn2(i,k))
         enddo
      enddo

      !1) Remap Tv, thetav, or TE
      if ( remap_te ) then
         if ( kord_tm==0 ) then
!----------------------------------
! map Total Energy using GMAO cubic
!----------------------------------
            call map1_cubic (km,   pe1,  te,       &
                 km,   pe2,  te,       &
                 is, ie, j, isd, ied, jsd, jed, akap, T_VAR=1, conserv=.true.)
         else
            call map_scalar(km,  peln(is,1,j),  te, gz(is:ie),   &
                 km,  pn2,           te,              &
                 is, ie, j, isd, ied, jsd, jed, 1, abs(kord_tm), cp_air*t_min)
         endif

      else
         if ( kord_tm<0 ) then
            ! Map t using logp
            call map_scalar(km,  peln(is,1,j),  pt, gz(is:ie),   &
                 km,  pn2,           pt,              &
                 is, ie, j, isd, ied, jsd, jed, &
                 1, abs(kord_tm), t_min)
         else
            ! Map pt using pe
            call map1_ppm (km,  pe1,  pt,  gz(is:ie),       &
                 km,  pe2,  pt,                  &
                 is, ie, j, isd, ied, jsd, jed, &
                 1, abs(kord_tm))
         endif
      endif

      !2) Map constituents

      if( nq > 5 ) then
           call mapn_tracer(nq, km, pe1, pe2, q, dp2, kord_tr, j,     &
                            is, ie, isd, ied, jsd, jed, 0., fill)
      elseif ( nq > 0 ) then
         ! Remap one tracer at a time
         do iq=1,nq
             call map1_q2(km, pe1, q(isd,jsd,1,iq),     &
                          km, pe2, q2, dp2,             &
                          is, ie, 0, kord_tr(iq), j, &
                          isd, ied, jsd, jed, 0.)
            if (fill) call fillz(ie-is+1, km, 1, q2, dp2)
            do k=1,km
               do i=is,ie
                  q(i,j,k,iq) = q2(i,k)
               enddo
            enddo
         enddo
      endif

      !3) Map W and density; recompute delz; limit w if needed
      if ( .not. hydrostatic ) then
         ! Remap vertical wind:
         if (kord_wz < 0) then
            call map1_ppm (km,   pe1,  w,  ws(is,j),   &
                 km,   pe2,  w,              &
                 is, ie, j, isd, ied, jsd, jed, &
                 -3, abs(kord_wz))
         else
            call map1_ppm (km,   pe1,  w,  ws(is,j),   &
                 km,   pe2,  w,              &
                 is, ie, j, isd, ied, jsd, jed, &
                 -2, abs(kord_wz))
         endif
         ! Remap delz for hybrid sigma-p coordinate
         call map1_ppm (km,   pe1, delz,  gz(is:ie),   & ! works
              km,   pe2, delz,              &
              is, ie, j, is,  ie,  js,  je,  &
              1, abs(kord_tm))
         do k=1,km
            do i=is,ie
               delz(i,j,k) = -delz(i,j,k)*dp2(i,k)
            enddo
         enddo

         !Fix excessive w - momentum conserving --- sjl
         ! gz(:) used here as a temporary array
         if ( w_limiter ) then
            do k=1,km
               do i=is,ie
                  w2(i,k) = w(i,j,k)
               enddo
            enddo
            do k=1, km-1
               do i=is,ie
                  if ( w2(i,k) > w_max ) then
                     gz(i) = (w2(i,k)-w_max) * dp2(i,k)
                     w2(i,k  ) = w_max
                     w2(i,k+1) = w2(i,k+1) + gz(i)/dp2(i,k+1)
                     !print*, ' W_LIMITER down: ', i,j,k, w2(i,k:k+1), w(i,j,k:k+1)
                  elseif ( w2(i,k) < w_min ) then
                     gz(i) = (w2(i,k)-w_min) * dp2(i,k)
                     w2(i,k  ) = w_min
                     w2(i,k+1) = w2(i,k+1) + gz(i)/dp2(i,k+1)
                     !print*, ' W_LIMITER down: ', i,j,k, w2(i,k:k+1), w(i,j,k:k+1)
                  endif
               enddo
            enddo
            do k=km, 2, -1
               do i=is,ie
                  if ( w2(i,k) > w_max ) then
                     gz(i) = (w2(i,k)-w_max) * dp2(i,k)
                     w2(i,k  ) = w_max
                     w2(i,k-1) = w2(i,k-1) + gz(i)/dp2(i,k-1)
                     !print*, ' W_LIMITER up: ', i,j,k, w2(i,k-1:k), w(i,j,k-1:k)
                  elseif ( w2(i,k) < w_min ) then
                     gz(i) = (w2(i,k)-w_min) * dp2(i,k)
                     w2(i,k  ) = w_min
                     w2(i,k-1) = w2(i,k-1) + gz(i)/dp2(i,k-1)
                     !print*, ' W_LIMITER up: ', i,j,k, w2(i,k-1:k), w(i,j,k-1:k)
                  endif
               enddo
            enddo
            do i=is,ie
               if (w2(i,1) > w_max*2. ) then
                  w2(i,1) = w_max*2 ! sink out of the top of the domain
                  !print*, ' W_LIMITER top limited: ', i,j,1, w2(i,1), w(i,j,1)
               elseif (w2(i,1) < w_min*2. ) then
                  w2(i,1) = w_min*2.
                  !print*, ' W_LIMITER top limited: ', i,j,1, w2(i,1), w(i,j,1)
               endif
            enddo
            do k=1,km
               do i=is,ie
                  w(i,j,k) = w2(i,k)
               enddo
            enddo
         endif
      endif

      ! 3.1) Update pressure variables
      do k=1,km+1
         do i=is,ie
            pk(i,j,k) = pk2(i,k)
         enddo
      enddo

      if ( last_step ) then
         ! 4.1) Start do_last_step
         ! save omega field in pe3
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

      ! 3.2) Compute pkz
      if ( .not. remap_te ) then
         if ( hydrostatic ) then
            do k=1,km
               do i=is,ie
                  pkz(i,j,k) = (pk2(i,k+1)-pk2(i,k))/(akap*(peln(i,k+1,j)-peln(i,k,j)))
               enddo
            enddo
         else
            ! Note: pt at this stage is T_v or T_m , unless kord_tm > 0
            do k=1,km
#ifdef MOIST_CAPPA
               call moist_cv(is,ie+1,isd,ied,jsd,jed, km, j, k, nwat, sphum, liq_wat, rainwat,    &
                    ice_wat, snowwat, graupel, q, gz(is:ie+1), cvm(is:ie+1))
               if ( kord_tm < 0 ) then
                  do i=is,ie
                     q_con(i,j,k) = gz(i)
                     cappa(i,j,k) = rdgas / ( rdgas + cvm(i)/(1.+r_vir*q(i,j,k,sphum)) )
                     pkz(i,j,k) = exp(cappa(i,j,k)*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
                  enddo
               else
                  do i=is,ie
                     q_con(i,j,k) = gz(i)
                     cappa(i,j,k) = rdgas / ( rdgas + cvm(i)/(1.+r_vir*q(i,j,k,sphum)) )
                     pkz(i,j,k) = exp(cappa(i,j,k)/(1.-cappa(i,j,k))*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
                  enddo
               endif
#else
               if ( kord_tm < 0 ) then
                  do i=is,ie
                     pkz(i,j,k) = exp(akap*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
                     ! Using dry pressure for the definition of the virtual potential temperature
                     !             pkz(i,j,k) = exp(akap*log(rrg*(1.-q(i,j,k,sphum))*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)/(1.+r_vir*q(i,j,k,sphum))))
                  enddo
               else
                  do i=is,ie
                     pkz(i,j,k) = exp(k1k*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
                     ! Using dry pressure for the definition of the virtual potential temperature
                     !             pkz(i,j,k) = exp(k1k*log(rrg*(1.-q(i,j,k,sphum))*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)/(1.+r_vir*q(i,j,k,sphum))))
                  enddo
               endif
#endif
            enddo
         endif
         if ( kord_tm > 0 ) then
            do k=1,km
               do i=is,ie
                  pt(i,j,k) = pt(i,j,k)*pkz(i,j,k) !Need Tv for energy calculations
               enddo
            enddo
         endif
      endif  ! end not remap_te

      ! 3.3) On last step Interpolate omega/pe3 (defined at pe0) to remapped cell center (dp2)
      if ( last_step ) then
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
      endif     ! end last_step

   endif !(j < je+1)

   do i=is,ie+1
      pe0(i,1) = pe(i,1,j)
   enddo

   !4) Remap winds

   ! 4.1) map u (on STAGGERED grid)
   do k=2,km+1
      do i=is,ie
         pe0(i,k) = 0.5*(pe(i,k,j-1)+pe1(i,k))
      enddo
   enddo

   do k=1,km+1
      bkh = 0.5*bk(k)
      do i=is,ie
         pe3(i,k) = ak(k) + bkh*(pe(i,km+1,j-1)+pe1(i,km+1))
      enddo
   enddo

   call map1_ppm( km, pe0(is:ie,:),   u,   gz(is:ie),   &
                  km, pe3(is:ie,:),   u,               &
                  is, ie, j, isd, ied, jsd, jed+1, &
                  -1, kord_mt)

   if (j < je+1) then

      ! 4.2) map v

      do i=is,ie+1
         pe3(i,1) = ak(1)
      enddo

      do k=2,km+1
         bkh = 0.5*bk(k)
         do i=is,ie+1
            pe0(i,k) =         0.5*(pe(i-1,k,   j)+pe(i,k,   j))
            pe3(i,k) = ak(k) + bkh*(pe(i-1,km+1,j)+pe(i,km+1,j))
         enddo
      enddo

      call map1_ppm (km, pe0,  v, gz(is:ie+1),    &
                     km, pe3,  v, is, ie+1,    &
                     j, isd, ied+1, jsd, jed, -1, kord_mt)

      ! 4a) update Tv and pkz from total energy (if remapping total energy)
      if ( remap_te ) then
         do i=is,ie
            phis(i,km+1) = hs(i,j)
         enddo
         ! calculate Tv from TE
         if ( hydrostatic ) then
            do k=km,1,-1
               do i=is,ie
                  tpe = te(i,j,k) - phis(i,k+1) - 0.25*gridstruct%rsin2(i,j)*(    &
                       u(i,j,k)**2+u(i,j+1,k)**2 + v(i,j,k)**2+v(i+1,j,k)**2 -  &
                       (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j) )
                  dlnp = rg*(peln(i,k+1,j) - peln(i,k,j))
                  pt(i,j,k)= tpe / (cp - pe2(i,k)*dlnp/delp(i,j,k))
                  pkz(i,j,k) = (pk2(i,k+1)-pk2(i,k))/(akap*(peln(i,k+1,j)-peln(i,k,j)))
                  phis(i,k) = phis(i,k+1) + dlnp*pt(i,j,k)
               enddo
            enddo           ! end k-loop
         else
            do k=km,1,-1
#ifdef MOIST_CAPPA
               call moist_cv(is,ie,isd,ied,jsd,jed, km, j, k, nwat, sphum, liq_wat, rainwat,    &
                    ice_wat, snowwat, graupel, q, gz(is:ie), cvm(is:ie))
               do i=is,ie
                  q_con(i,j,k) = gz(i)
                  cappa(i,j,k) = rdgas / ( rdgas + cvm(i)/(1.+r_vir*q(i,j,k,sphum)) )
               enddo
#endif
               do i=is,ie
                  phis(i,k) = phis(i,k+1) - delz(i,j,k)*grav
                  tpe = te(i,j,k) - 0.5*(phis(i,k)+phis(i,k+1)) - 0.5*w(i,j,k)**2 - 0.25*gridstruct%rsin2(i,j)*(    &
                       u(i,j,k)**2+u(i,j+1,k)**2 + v(i,j,k)**2+v(i+1,j,k)**2 -  &
                       (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j) )
#ifdef MOIST_CAPPA
                  pt(i,j,k)= tpe / cvm(i)*(1.+r_vir*q(i,j,k,sphum))*(1.-gz(i))
                  pkz(i,j,k) = exp(cappa(i,j,k)*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
#else
                  pt(i,j,k)= tpe / cv_air *(1.+r_vir*q(i,j,k,sphum))
                  pkz(i,j,k) = exp(akap*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
#endif
               enddo

            enddo           ! end k-loop
         endif
      endif
   endif ! (j < je+1)

   do k=1,km
      do i=is,ie
         pe4(i,j,k) = pe2(i,k+1)
      enddo
   enddo


  enddo !j-loop

  !6) Energy fixer
!$OMP parallel do default(none) shared(is,ie,js,je,km,pe4,pe)
  do k=2,km
     do j=js,je
        do i=is,ie
           pe(i,k,j) = pe4(i,j,k-1)
        enddo
     enddo
  enddo

  dtmp = 0.

  if( last_step .and. (.not.do_adiabatic_init)  ) then

     if ( consv > consv_min ) then

!$OMP parallel do default(none) shared(is,ie,js,je,km,ptop,u,v,pe,isd,ied,jsd,jed,te_2d,delp, &
!$OMP                                  hydrostatic,hs,rg,pt,peln,cp,delz,nwat,rainwat,liq_wat, &
!$OMP                                  ice_wat,snowwat,graupel,q_con,r_vir,sphum,w,pk,pkz,zsum1, &
!$OMP                                  zsum0,te0_2d,gridstruct,q,kord_tm,te,remap_te) &
!$OMP                          private(cvm,gz,phis)
        do j=js,je
           if ( remap_te ) then
              do i=is,ie
                 te_2d(i,j) = te(i,j,1)*delp(i,j,1)
              enddo
              do k=2,km
                 do i=is,ie
                    te_2d(i,j) = te_2d(i,j) + te(i,j,k)*delp(i,j,k)
                 enddo
              enddo
           else
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
                            0.25*gridstruct%rsin2(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                            v(i,j,k)**2+v(i+1,j,k)**2 -  &
                            (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j)))
                    enddo
                 enddo
              else
                 do i=is,ie
                    te_2d(i,j) = 0.
                    phis(i,km+1) = hs(i,j)
                 enddo
                 do k=km,1,-1
                    do i=is,ie
                       phis(i,k) = phis(i,k+1) - grav*delz(i,j,k)
                    enddo
                 enddo

                 do k=1,km
#ifdef USE_COND
                    call moist_cv(is,ie,isd,ied,jsd,jed, km, j, k, nwat, sphum, liq_wat, rainwat,    &
                         ice_wat, snowwat, graupel, q, gz(is:ie), cvm(is:ie))
                    do i=is,ie
                       ! KE using 3D winds:
                       q_con(i,j,k) = gz(i)
                       te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(cvm(i)*pt(i,j,k)/((1.+r_vir*q(i,j,k,sphum))*(1.-gz(i))) + &
                            0.5*(phis(i,k)+phis(i,k+1) + w(i,j,k)**2 + 0.5*gridstruct%rsin2(i,j)*( &
                            u(i,j,k)**2+u(i,j+1,k)**2 + v(i,j,k)**2+v(i+1,j,k)**2 -  &
                            (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j))))
                    enddo
#else
                    do i=is,ie
                       te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(cv_air*pt(i,j,k)/(1.+r_vir*q(i,j,k,sphum)) + &
                            0.5*(phis(i,k)+phis(i,k+1) + w(i,j,k)**2 + 0.5*gridstruct%rsin2(i,j)*( &
                            u(i,j,k)**2+u(i,j+1,k)**2 + v(i,j,k)**2+v(i+1,j,k)**2 -  &
                            (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j))))
                    enddo
#endif
                 enddo   ! k-loop
              endif  ! end non-hydro
           endif  ! end non remapping te

           do i=is,ie
              te_2d(i,j) = te0_2d(i,j) - te_2d(i,j)
              zsum1(i,j) = pkz(i,j,1)*delp(i,j,1)
           enddo
           do k=2,km
              do i=is,ie
                 zsum1(i,j) = zsum1(i,j) + pkz(i,j,k)*delp(i,j,k)
              enddo
           enddo
           if ( hydrostatic ) then
              do i=is,ie
                 zsum0(i,j) = ptop*(pk(i,j,1)-pk(i,j,km+1)) + zsum1(i,j)
              enddo
           endif

        enddo   ! j-loop

        dtmp = consv*g_sum(domain, te_2d, is, ie, js, je, ng, gridstruct%area_64, 0, reproduce=.true.)
        E_Flux = dtmp / (grav*pdt*4.*pi*radius**2)    ! unit: W/m**2
        ! Note pdt is "phys" time step
        if ( hydrostatic ) then !AM4 version multiplies in cp or cv_air to g_sum here
           dtmp = dtmp / g_sum(domain, zsum0, is, ie, js, je, ng, gridstruct%area_64, 0, reproduce=.true.)
        else
           dtmp = dtmp / g_sum(domain, zsum1, is, ie, js, je, ng, gridstruct%area_64, 0, reproduce=.true.)
        endif

     elseif ( consv < -consv_min ) then

!$OMP parallel do default(none) shared(is,ie,js,je,km,pkz,delp,zsum1,zsum0,ptop,pk,hydrostatic)
        do j=js,je
           do i=is,ie
              zsum1(i,j) = pkz(i,j,1)*delp(i,j,1)
           enddo
           do k=2,km
              do i=is,ie
                 zsum1(i,j) = zsum1(i,j) + pkz(i,j,k)*delp(i,j,k)
              enddo
           enddo
           if ( hydrostatic ) then
              do i=is,ie
                 zsum0(i,j) = ptop*(pk(i,j,1)-pk(i,j,km+1)) + zsum1(i,j)
              enddo
           endif
        enddo

        E_Flux = consv
        if ( hydrostatic ) then !AM4 multiplies in cp or cv_air to g_sum here
           dtmp = E_flux*(grav*pdt*4.*pi*radius**2) /    &
                g_sum(domain, zsum0,  is, ie, js, je, ng, gridstruct%area_64, 0, reproduce=.true.)
        else
           dtmp = E_flux*(grav*pdt*4.*pi*radius**2) /    &
                g_sum(domain, zsum1,  is, ie, js, je, ng, gridstruct%area_64, 0, reproduce=.true.)
        endif
     endif        ! end consv check
  endif        ! end last_step check

!-----------------------------------------------------------------------
! Intermediate Physics >>>
! Note: if intemediate physics is disable, cloud fraction will be zero at the first time step.
!-----------------------------------------------------------------------
    if (do_intermediate_phys) then
        call timing_on('INTERMEDIATE_PHYS')
        call intermediate_phys (is, ie, js, je, isd, ied, jsd, jed, km, npx, npy, nq, nwat, &
                 c2l_ord, mdt, consv, akap, ptop, pfull, hs, te0_2d, u, &
                 v, w, pt, delp, delz, q_con, cappa, q, pkz, r_vir, te_err, tw_err, &
                 inline_mp, gridstruct, domain, bd, hydrostatic, do_adiabatic_init, &
                 do_inline_mp, do_sat_adj, last_step, do_fast_phys, consv_checker, adj_mass_vmr)
        call timing_off('INTERMEDIATE_PHYS')
    endif

!-----------------------------------------------------------------------
! <<< Intermediate Physics
!-----------------------------------------------------------------------

  if ( last_step ) then
       ! 9a) Convert T_v/T_m to T if last_step
!!!  if ( is_master() ) write(*,*) 'dtmp=', dtmp, nwat
!$OMP parallel do default(none) shared(is,ie,js,je,km,isd,ied,jsd,jed,hydrostatic,pt,adiabatic,cp, &
!$OMP                                  nwat,rainwat,liq_wat,ice_wat,snowwat,graupel,r_vir,&
!$OMP                                  sphum,pkz,dtmp,q) &
!$OMP                          private(cvm,gz)
     do k=1,km
        do j=js,je
           if (hydrostatic) then !This is re-factored from AM4 so answers may be different
              do i=is,ie
                 pt(i,j,k) = (pt(i,j,k)+dtmp/cp*pkz(i,j,k)) / (1.+r_vir*q(i,j,k,sphum))
              enddo
           else
#ifdef USE_COND
              call moist_cv(is,ie,isd,ied,jsd,jed, km, j, k, nwat, sphum, liq_wat, rainwat,    &
                   ice_wat, snowwat, graupel, q, gz(is:ie), cvm(is:ie))
              do i=is,ie
                 pt(i,j,k) = (pt(i,j,k)+dtmp/cvm(i)*pkz(i,j,k))/((1.+r_vir*q(i,j,k,sphum))*(1.-gz(i)))
              enddo
#else
              if ( .not. adiabatic ) then
                 do i=is,ie
                    pt(i,j,k) = (pt(i,j,k)+dtmp/cv_air*pkz(i,j,k)) / (1.+r_vir*q(i,j,k,sphum))
                 enddo
              endif
#endif
           endif
        enddo   ! j-loop
     enddo  ! k-loop
  else
     ! 9b) not last_step: convert T_v/T_m back to theta_v/theta_m for dyn_core
!$OMP parallel do default(none) shared(is,ie,js,je,km,pkz,pt)
     do k=1,km
        do j=js,je
           do i=is,ie
              pt(i,j,k) = pt(i,j,k)/pkz(i,j,k)
           enddo
        enddo
     enddo
  endif

 end subroutine Lagrangian_to_Eulerian


 subroutine compute_total_energy(is, ie, js, je, isd, ied, jsd, jed, km,       &
                                 u, v, w, delz, pt, delp, q, qc, pe, peln, hs, &
                                 rsin2_l, cosa_s_l, &
                                 r_vir,  cp, rg, hlv, te_2d, ua, va, teq, &
                                 moist_phys, nwat, sphum, liq_wat, rainwat, ice_wat, snowwat, graupel, hydrostatic, id_te)
!------------------------------------------------------
! Compute vertically integrated total energy per column
!------------------------------------------------------
! !INPUT PARAMETERS:
   integer,  intent(in):: km, is, ie, js, je, isd, ied, jsd, jed, id_te
   integer,  intent(in):: sphum, liq_wat, ice_wat, rainwat, snowwat, graupel, nwat
   real, intent(inout), dimension(isd:ied,jsd:jed,km):: ua, va
   real, intent(in), dimension(isd:ied,jsd:jed,km):: pt, delp
   real, intent(in), dimension(isd:ied,jsd:jed,km,*):: q
   real, intent(in), dimension(isd:ied,jsd:jed,km):: qc
   real, intent(inout)::  u(isd:ied,  jsd:jed+1,km)
   real, intent(inout)::  v(isd:ied+1,jsd:jed,  km)
   real, intent(in)::  w(isd:,jsd:,1:)   ! vertical velocity (m/s)
   real, intent(in):: delz(is:,js:,1:)
   real, intent(in):: hs(isd:ied,jsd:jed)  ! surface geopotential
   real, intent(in)::   pe(is-1:ie+1,km+1,js-1:je+1) ! pressure at layer edges
   real, intent(in):: peln(is:ie,km+1,js:je)  ! log(pe)
   real, intent(in):: cp, rg, r_vir, hlv
   real, intent(in) :: rsin2_l(isd:ied, jsd:jed)
   real, intent(in) :: cosa_s_l(isd:ied, jsd:jed)
   logical, intent(in):: moist_phys, hydrostatic
! Output:
   real, intent(out):: te_2d(is:ie,js:je)   ! vertically integrated TE
   real, intent(out)::   teq(is:ie,js:je)   ! Moist TE
! Local
   real, dimension(is:ie,km):: tv
   real  phiz(is:ie,km+1)
   real  cvm(is:ie), qd(is:ie)
   integer i, j, k

!----------------------
! Output lat-lon winds:
!----------------------
!  call cubed_to_latlon(u, v, ua, va, dx, dy, rdxa, rdya, km, flagstruct%c2l_ord)

!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,km,hydrostatic,hs,pt,qc,rg,peln,te_2d, &
!$OMP                                  pe,delp,cp,rsin2_l,u,v,cosa_s_l,delz,moist_phys,w, &
!$OMP                                  q,nwat,liq_wat,rainwat,ice_wat,snowwat,graupel,sphum)   &
!$OMP                          private(phiz, tv, cvm, qd)
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
                           0.25*rsin2_l(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +      &
                                            v(i,j,k)**2+v(i+1,j,k)**2 -      &
                       (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s_l(i,j)))
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
     !TODO moist_phys doesn't seem to make a difference --- lmh 13may21
     if ( moist_phys ) then
     do k=1,km
#ifdef MOIST_CAPPA
        call moist_cv(is,ie,isd,ied,jsd,jed, km, j, k, nwat, sphum, liq_wat, rainwat,    &
                      ice_wat, snowwat, graupel, q, qd, cvm)
#endif
        do i=is,ie
#ifdef MOIST_CAPPA
           te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*( cvm(i)*pt(i,j,k) +  &
#else
           te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*( cv_air*pt(i,j,k) +  &
#endif
                        0.5*(phiz(i,k)+phiz(i,k+1)+w(i,j,k)**2+0.5*rsin2_l(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                        v(i,j,k)**2+v(i+1,j,k)**2-(u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s_l(i,j))))
        enddo
     enddo
     else
       do k=1,km
          do i=is,ie
             te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*( cv_air*pt(i,j,k) +  &
                          0.5*(phiz(i,k)+phiz(i,k+1)+w(i,j,k)**2+0.5*rsin2_l(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                          v(i,j,k)**2+v(i+1,j,k)**2-(u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s_l(i,j))))
          enddo
       enddo
     endif
     endif
  enddo

!-------------------------------------
! Diganostics computation for moist TE
!-------------------------------------
  if( id_te>0 ) then
!$OMP parallel do default(none) shared(is,ie,js,je,teq,te_2d,moist_phys,km,hlv,sphum,q,delp)
      do j=js,je
         do i=is,ie
            teq(i,j) = te_2d(i,j)
         enddo
         if ( moist_phys ) then
           do k=1,km
              do i=is,ie
                 teq(i,j) = teq(i,j) + hlv*q(i,j,k,sphum)*delp(i,j,k)
              enddo
           enddo
         endif
      enddo
   endif

  end subroutine compute_total_energy


  subroutine pkez(km, ifirst, ilast, jfirst, jlast, j, &
                  pe, pk, akap, peln, pkz, ptop)

! !INPUT PARAMETERS:
   integer, intent(in):: km, j
   integer, intent(in):: ifirst, ilast        ! Latitude strip
   integer, intent(in):: jfirst, jlast        ! Latitude strip
   real, intent(in):: akap
   real, intent(in):: pe(ifirst-1:ilast+1,km+1,jfirst-1:jlast+1)
   real, intent(in):: pk(ifirst:ilast,jfirst:jlast,km+1)
   real, intent(IN):: ptop
! !OUTPUT
   real, intent(out):: pkz(ifirst:ilast,jfirst:jlast,km)
   real, intent(inout):: peln(ifirst:ilast, km+1, jfirst:jlast)   ! log (pe)
! Local
   real pk2(ifirst:ilast, km+1)
   real pek
   real lnp
   real ak1
   integer i, k

   ak1 = (akap + 1.) / akap

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

 end subroutine pkez

 subroutine map_scalar( km,   pe1,    q1,   qs,           &
                        kn,   pe2,    q2,   i1, i2,       &
                        j,  ibeg, iend, jbeg, jend,      &
                        iv,  kord, q_min)
! iv=1
 integer, intent(in) :: i1                ! Starting longitude
 integer, intent(in) :: i2                ! Finishing longitude
 integer, intent(in) :: iv                ! Mode: 0 == constituents  1 == temp
                                          !       2 == remap temp with cs scheme
                                          !      -2 or -3 == w with lower bc
 integer, intent(in) :: kord              ! Method order
 integer, intent(in) :: j                 ! Current latitude
 integer, intent(in) :: ibeg, iend, jbeg, jend
 integer, intent(in) :: km                ! Original vertical dimension
 integer, intent(in) :: kn                ! Target vertical dimension
 real, intent(in) ::   qs(i1:i2)       ! bottom BC
 real, intent(in) ::  pe1(i1:i2,km+1)  ! pressure at layer edges
                                       ! (from model top to bottom surface)
                                       ! in the original vertical coordinate
 real, intent(in) ::  pe2(i1:i2,kn+1)  ! pressure at layer edges
                                       ! (from model top to bottom surface)
                                       ! in the new vertical coordinate
 real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
 real, intent(inout)::  q2(ibeg:iend,jbeg:jend,kn) ! Field output
 real, intent(in):: q_min

! !DESCRIPTION:
! IV = 0: constituents: enforce positivity in interface values and reconstruction
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
      call  scalar_profile( qs, q4, dp1, km, i1, i2, iv, kord, q_min )
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

 end subroutine map_scalar


 subroutine map1_ppm( km,   pe1,    q1,   qs,           &
                      kn,   pe2,    q2,   i1, i2,       &
                      j,    ibeg, iend, jbeg, jend,     &
                      iv, kord)
 integer, intent(in) :: i1                ! Starting longitude
 integer, intent(in) :: i2                ! Finishing longitude
 integer, intent(in) :: iv                ! Mode: 0 == constituents  1 == ???
                                          !       2 == remap temp with cs scheme
                                          !      -1 == vertical velocity, with bottom BC
 integer, intent(in) :: kord              ! Method order
 integer, intent(in) :: j                 ! Current latitude
 integer, intent(in) :: ibeg, iend, jbeg, jend
 integer, intent(in) :: km                ! Original vertical dimension
 integer, intent(in) :: kn                ! Target vertical dimension
 real, intent(in) ::   qs(i1:i2)       ! bottom BC (only used if iv == -2 )
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
      call  cs_profile( qs, q4, dp1, km, i1, i2, iv, kord )
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


!Multi-tracer remapping (much faster)
!ONLY supports cubic-spline remapping
 subroutine mapn_tracer(nq, km, pe1, pe2, q1, dp2, kord, j,     &
                        i1, i2, isd, ied, jsd, jed,             &
                        q_min, fill)
! !INPUT PARAMETERS:
      integer, intent(in):: km                ! vertical dimension
      integer, intent(in):: j, nq, i1, i2
      integer, intent(in):: isd, ied, jsd, jed
      integer, intent(in):: kord(nq)
      real, intent(in)::  pe1(i1:i2,km+1)     ! pressure at layer edges
                                              ! (from model top to bottom surface)
                                              ! in the original vertical coordinate
      real, intent(in)::  pe2(i1:i2,km+1)     ! pressure at layer edges
                                              ! (from model top to bottom surface)
                                              ! in the new vertical coordinate
      real, intent(in)::  dp2(i1:i2,km)
      real, intent(in)::  q_min
      logical, intent(in):: fill
      real, intent(inout):: q1(isd:ied,jsd:jed,km,nq) ! Field input
! !LOCAL VARIABLES:
      real:: q4(4,i1:i2,km,nq)
      real:: q2(i1:i2,km,nq) ! Field output
      real:: qsum(nq)
      real:: dp1(i1:i2,km)
      real:: qs(i1:i2)
      real:: pl, pr, dp, esl, fac1, fac2
      integer:: i, k, l, m, k0, iq

      do k=1,km
         do i=i1,i2
            dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         enddo
      enddo

      do iq=1,nq
         do k=1,km
            do i=i1,i2
               q4(1,i,k,iq) = q1(i,j,k,iq)
            enddo
         enddo
         call  scalar_profile( qs, q4(1,i1,1,iq), dp1, km, i1, i2, 0, kord(iq), q_min )
      enddo

! Mapping
      do 1000 i=i1,i2
         k0 = 1
      do 555 k=1,km
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1)) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if(pe2(i,k+1) <= pe1(i,l+1)) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            fac1 = pr + pl
            fac2 = r3*(pr*fac1 + pl*pl)
            fac1 = 0.5*fac1
            do iq=1,nq
               q2(i,k,iq) = q4(2,i,l,iq) + (q4(4,i,l,iq)+q4(3,i,l,iq)-q4(2,i,l,iq))*fac1  &
                                         -  q4(4,i,l,iq)*fac2
            enddo
            k0 = l
            goto 555
          else
! Fractional area...
            dp = pe1(i,l+1) - pe2(i,k)
            fac1 = 1. + pl
            fac2 = r3*(1.+pl*fac1)
            fac1 = 0.5*fac1
            do iq=1,nq
               qsum(iq) = dp*(q4(2,i,l,iq) + (q4(4,i,l,iq)+   &
                              q4(3,i,l,iq) - q4(2,i,l,iq))*fac1 - q4(4,i,l,iq)*fac2)
            enddo
            do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
               if(pe2(i,k+1) > pe1(i,m+1) ) then
                                                   ! Whole layer..
                  do iq=1,nq
                     qsum(iq) = qsum(iq) + dp1(i,m)*q4(1,i,m,iq)
                  enddo
               else
                  dp = pe2(i,k+1)-pe1(i,m)
                  esl = dp / dp1(i,m)
                  fac1 = 0.5*esl
                  fac2 = 1.-r23*esl
                  do iq=1,nq
                     qsum(iq) = qsum(iq) + dp*( q4(2,i,m,iq) + fac1*(         &
                                q4(3,i,m,iq)-q4(2,i,m,iq)+q4(4,i,m,iq)*fac2 ) )
                  enddo
                  k0 = m
                  goto 123
               endif
            enddo
            goto 123
          endif
      endif
100   continue
123   continue
      do iq=1,nq
         q2(i,k,iq) = qsum(iq) / dp2(i,k)
      enddo
555   continue
1000  continue

  if (fill) call fillz(i2-i1+1, km, nq, q2, dp2)

  do iq=1,nq
!    if (fill) call fillz(i2-i1+1, km, 1, q2(i1,1,iq), dp2)
     do k=1,km
        do i=i1,i2
           q1(i,j,k,iq) = q2(i,k,iq)
        enddo
     enddo
  enddo

 end subroutine mapn_tracer


 !This routine remaps a single tracer
 subroutine map1_q2(km,   pe1,   q1,            &
                    kn,   pe2,   q2,   dp2,     &
                    i1,   i2,    iv,   kord, j, &
                    ibeg, iend, jbeg, jend,     &
                    q_min )


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
      real, intent(in) ::  q_min
! !INPUT/OUTPUT PARAMETERS:
      real, intent(inout):: q2(i1:i2,kn) ! Field output
! !LOCAL VARIABLES:
      real   qs(i1:i2)
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
      call  scalar_profile( qs, q4, dp1, km, i1, i2, iv, kord, q_min )
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


 !Currently this routine is only called with kord = 4,
 ! --- lmh 9 june 21
 subroutine remap_2d(km,   pe1,   q1,        &
                     kn,   pe2,   q2,        &
                     i1,   i2,               &
                     iv,   kord)
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
   real   qs(i1:i2)
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
      call  cs_profile( qs, q4, dp1, km, i1, i2, iv, kord )
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

   do i=i1,i2
      k0 = 1
      do 555 k=1,kn
#ifdef OLD_TOP_EDGE
         if( pe2(i,k+1) <= pe1(i,1) ) then
! Entire grid above old ptop
             q2(i,k) = q4(2,i,1)
         elseif( pe2(i,k) < pe1(i,1) .and. pe2(i,k+1)>pe1(i,1) ) then
! Partially above old ptop:
             q2(i,k) = q1(i,1)
#else
         if( pe2(i,k) <= pe1(i,1) ) then
! above old ptop:
             q2(i,k) = q1(i,1)
#endif
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

 !scalar_profile and cs_profile differ ONLY in that scalar_profile
 ! accepts a qmin argument. (Unfortunately I was not able to make
 ! qmin an optional argument in scalar_profile.) --- lmh summer 2020
 subroutine scalar_profile(qs, a4, delp, km, i1, i2, iv, kord, qmin)
! Optimized vertical profile reconstruction:
! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
 integer, intent(in):: i1, i2
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
 integer, intent(in):: kord
 real, intent(in)   ::   qs(i1:i2)
 real, intent(in)   :: delp(i1:i2,km)     ! layer pressure thickness
 real, intent(inout):: a4(4,i1:i2,km)     ! Interpolated values
 real, intent(in):: qmin
!-----------------------------------------------------------------------
 logical, dimension(i1:i2,km):: extm, ext5, ext6
 real  gam(i1:i2,km)
 real    q(i1:i2,km+1)
 real   d4(i1:i2)
 real   bet, a_bot, grat
 real   pmp_1, lac_1, pmp_2, lac_2, x0, x1
 integer i, k, im

 !Compute interface values (\hat{q})
 ! iv=-2 and -3 introduce the lower BC
 ! iv=-2 also uses a simpler calculation
 !       dropping a lot of metric terms
 if ( iv .eq. -2 ) then
      do i=i1,i2
         gam(i,2) = 0.5
           q(i,1) = 1.5*a4(1,i,1)
      enddo
      do k=2,km-1
         do i=i1, i2
                  grat = delp(i,k-1) / delp(i,k)
                   bet =  2. + grat + grat - gam(i,k)
                q(i,k) = (3.*(a4(1,i,k-1)+a4(1,i,k)) - q(i,k-1))/bet
            gam(i,k+1) = grat / bet
         enddo
      enddo
      do i=i1,i2
            grat = delp(i,km-1) / delp(i,km)
         q(i,km) = (3.*(a4(1,i,km-1)+a4(1,i,km)) - grat*qs(i) - q(i,km-1)) /  &
                   (2. + grat + grat - gam(i,km))
         q(i,km+1) = qs(i)
      enddo
      do k=km-1,1,-1
        do i=i1,i2
           q(i,k) = q(i,k) - gam(i,k+1)*q(i,k+1)
        enddo
      enddo
 else
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
 endif

!Perfectly linear scheme
 if ( abs(kord) == 17 ) then
  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
  enddo
  return
 endif

  im = i2 - i1 + 1

  ! Apply *large-scale* constraints to \hat{q}

  !Upper BC for all schemes
  do i=i1,i2
     q(i,2) = min( q(i,2), max(a4(1,i,1), a4(1,i,2)) )
     q(i,2) = max( q(i,2), min(a4(1,i,1), a4(1,i,2)) )
  enddo

  do k=2,km
     do i=i1,i2
        gam(i,k) = a4(1,i,k) - a4(1,i,k-1) !\delta \bar{q}
     enddo
  enddo

! Interior:
  do k=3,km-1
     do i=i1,i2
        if ( abs(kord) >= 14 .or. gam(i,k-1)*gam(i,k+1)>0. ) then
! Apply large-scale constraint to ALL fields if not local max/min
!  first guess interface values cannot exceeed values
!  of adjacent cells
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

! Bottom BC for all schemes:
  do i=i1,i2
     q(i,km) = min( q(i,km), max(a4(1,i,km-1), a4(1,i,km)) )
     q(i,km) = max( q(i,km), min(a4(1,i,km-1), a4(1,i,km)) )
  enddo

  !Set up in-cell reconstruction
  !initially continuous (AL(k) = AR(k-1))
  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
     enddo
  enddo

  !Flags for different extremum/2dz conditions
  ! estimated from first-guess edge values
  do k=1,km
     if ( k==1 .or. k==km ) then
       do i=i1,i2
          extm(i,k) = (a4(2,i,k)-a4(1,i,k)) * (a4(3,i,k)-a4(1,i,k)) > 0.
       enddo
     else
       do i=i1,i2
          extm(i,k) = gam(i,k)*gam(i,k+1) < 0.
       enddo
     endif
     if ( abs(kord) > 9 ) then
       do i=i1,i2
          x0 = 2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k))
          x1 = abs(a4(2,i,k)-a4(3,i,k))
          a4(4,i,k) = 3.*x0
          ext5(i,k) = abs(x0) > x1
          ext6(i,k) = abs(a4(4,i,k)) > x1
       enddo
     endif
  enddo

! Apply subgrid constraints:
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
! Top 2 and bottom 2 layers always use monotonic mapping

  select case (iv)

  case (0)
     do i=i1,i2
        a4(2,i,1) = max(0., a4(2,i,1))
     enddo
  case (-1)
      do i=i1,i2
         if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
      enddo
  case (2)
     do i=i1,i2
        a4(2,i,1) = a4(1,i,1)
        a4(3,i,1) = a4(1,i,1)
        a4(4,i,1) = 0.
     enddo
  end select !iv

  if ( iv/=2 ) then
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
      select case (abs(kord))

      case (0:8)
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

     case (9)
       do i=i1,i2
          if ( extm(i,k) .and. extm(i,k-1) ) then
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else if ( extm(i,k) .and. extm(i,k+1) ) then
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else if ( extm(i,k) .and. a4(1,i,k)<qmin ) then
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
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
              a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
            endif
          endif
       enddo
     case(10) !restored AM4 case 10
       do i=i1,i2
          if( extm(i,k) ) then
              if( a4(1,i,k)<qmin .or. extm(i,k-1) .or. extm(i,k+1) ) then
! grid-scale 2-delta-z wave detected; or q is too small -> ehance vertical mixing
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
    case(11)
       do i=i1,i2
         if ( ext5(i,k) .and. (ext5(i,k-1).or.ext5(i,k+1).or.a4(1,i,k)<qmin) ) then
! Noisy region:
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
         else
              a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         endif
       enddo
    case(12) !post-AM4 case 10
       do i=i1,i2
          if( ext5(i,k) ) then
              if( ext5(i,k-1) .or. ext5(i,k+1) ) then
                   a4(2,i,k) = a4(1,i,k)
                   a4(3,i,k) = a4(1,i,k)
              elseif ( ext6(i,k-1) .or. ext6(i,k+1) ) then
                   pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                   lac_1 = pmp_1 + 1.5*gam(i,k+2)
                   a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                                  max(a4(1,i,k), pmp_1, lac_1) )
                   pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                   lac_2 = pmp_2 - 1.5*gam(i,k-1)
                   a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                                  max(a4(1,i,k), pmp_2, lac_2) )
              endif
          elseif( ext6(i,k) ) then
              if( ext5(i,k-1) .or. ext5(i,k+1) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
                  a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                                 max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
                  a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                                 max(a4(1,i,k), pmp_2, lac_2) )
              endif
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
    case(13) !former 14: no subgrid limiter

       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo

    case (14) !strict monotonicity constraint
       call cs_limiters(im, extm(i1,k), a4(1,i1,k), 2)
    case (15)
       call cs_limiters(im, extm(i1,k), a4(1,i1,k), 1)
    case default
       call mpp_error(FATAL, " kord not implemented")
    end select

! Additional constraint to ensure positivity
     if ( iv==0 .and. abs(kord) <= 13 ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 0)

  enddo      ! k-loop

!----------------------------------
! Bottom layer subgrid constraints:
!----------------------------------
  select case (iv)
  case(0)
     do i=i1,i2
        a4(3,i,km) = max(0., a4(3,i,km))
     enddo
  case(-1)
     do i=i1,i2
        if ( a4(3,i,km)*a4(1,i,km) <= 0. )  a4(3,i,km) = 0.
     enddo
  end select

  do k=km-1,km
     do i=i1,i2
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
     if(k==(km-1)) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 2)
     if(k== km   ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 1)
  enddo

 end subroutine scalar_profile


 subroutine cs_profile(qs, a4, delp, km, i1, i2, iv, kord)
! Optimized vertical profile reconstruction:
! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
 integer, intent(in):: i1, i2
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: iv      ! iv =-2: vertical velocity
                               ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
 integer, intent(in):: kord
 real, intent(in)   ::   qs(i1:i2)
 real, intent(in)   :: delp(i1:i2,km)     ! layer pressure thickness
 real, intent(inout):: a4(4,i1:i2,km)     ! Interpolated values
!-----------------------------------------------------------------------
 logical, dimension(i1:i2,km):: extm, ext5, ext6
 real  gam(i1:i2,km)
 real    q(i1:i2,km+1) ! interface values
 real   d4(i1:i2)
 real   bet, a_bot, grat
 real   pmp_1, lac_1, pmp_2, lac_2, x0, x1
 integer i, k, im

 if ( iv .eq. -2 ) then
      do i=i1,i2
         gam(i,2) = 0.5
           q(i,1) = 1.5*a4(1,i,1)
      enddo
      do k=2,km-1
         do i=i1, i2
                  grat = delp(i,k-1) / delp(i,k)
                   bet =  2. + grat + grat - gam(i,k)
                q(i,k) = (3.*(a4(1,i,k-1)+a4(1,i,k)) - q(i,k-1))/bet
            gam(i,k+1) = grat / bet
         enddo
      enddo
      do i=i1,i2
            grat = delp(i,km-1) / delp(i,km)
         q(i,km) = (3.*(a4(1,i,km-1)+a4(1,i,km)) - grat*qs(i) - q(i,km-1)) /  &
                   (2. + grat + grat - gam(i,km))
         q(i,km+1) = qs(i)
      enddo
      do k=km-1,1,-1
        do i=i1,i2
           q(i,k) = q(i,k) - gam(i,k+1)*q(i,k+1)
        enddo
      enddo

else ! all others
  do i=i1,i2
         grat = delp(i,2) / delp(i,1)   ! grid ratio
          bet = grat*(grat+0.5)
       q(i,1) = ( (grat+grat)*(grat+1.)*a4(1,i,1) + a4(1,i,2) ) / bet
     gam(i,1) = ( 1. + grat*(grat+1.5) ) / bet
  enddo

  if (iv.eq.-3) then !LBC for vertical velocities
    do k=2,km-1
       do i=i1,i2
             d4(i) = delp(i,k-1) / delp(i,k)
               bet =  2. + d4(i) + d4(i) - gam(i,k-1)
            q(i,k) = ( 3.*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - q(i,k-1) )/bet
          gam(i,k) = d4(i) / bet
       enddo
    enddo

    do i=i1,i2
       !    a_bot = 1. + d4(i)*(d4(i)+1.5)
       !q(i,km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1,i,km)+a4(1,i,km-1)-a_bot*q(i,km))  &
       !          / ( d4(i)*(d4(i)+0.5) - a_bot*gam(i,km) )
       d4(i) = delp(i,km-1) / delp(i,km)
       bet =  2. + d4(i) + d4(i) - gam(i,km-1)
       grat = delp(i,km-1) / delp(i,km)
       q(i,km) = ( 3.*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - grat*qs(i) - q(i,k-1) )/bet
       q(i,km+1) = qs(i)
    enddo

  else ! all others
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
  endif

  do k=km,1,-1
     do i=i1,i2
        q(i,k) = q(i,k) - gam(i,k)*q(i,k+1)
     enddo
  enddo
 endif
!----- Perfectly linear scheme --------------------------------
 if ( abs(kord) == 17 ) then
  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
  enddo
  return
 endif
!----- Perfectly linear scheme --------------------------------

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
        gam(i,k) = a4(1,i,k) - a4(1,i,k-1) ! now dq
     enddo
  enddo

! Interior:
  do k=3,km-1
     do i=i1,i2
        if ( abs(kord) >= 14 .or. gam(i,k-1)*gam(i,k+1)>0. ) then
! Apply large-scale constraint to ALL fields if not local max/min
! OR for the strictly monotone schemes
             q(i,k) = min( q(i,k), max(a4(1,i,k-1),a4(1,i,k)) )
             q(i,k) = max( q(i,k), min(a4(1,i,k-1),a4(1,i,k)) )
        else
          if ( gam(i,k-1) > 0. ) then
! There exists a local max
               q(i,k) = max(q(i,k), min(a4(1,i,k-1),a4(1,i,k)))
          else
! There exists a local min
                 q(i,k) = min(q(i,k), max(a4(1,i,k-1),a4(1,i,k)))
               if ( iv==0 ) q(i,k) = max(0., q(i,k)) ! positive-definite
          endif
        endif
     enddo
  enddo

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

  do k=1,km
     if ( k==1 .or. k==km ) then
       do i=i1,i2
          extm(i,k) = (a4(2,i,k)-a4(1,i,k)) * (a4(3,i,k)-a4(1,i,k)) > 0.
       enddo
     else
       do i=i1,i2
          extm(i,k) = gam(i,k)*gam(i,k+1) < 0.
       enddo
     endif
     if ( abs(kord) > 9 ) then
       do i=i1,i2
          x0 = 2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k))
          x1 = abs(a4(2,i,k)-a4(3,i,k))
          a4(4,i,k) = 3.*x0
          ext5(i,k) = abs(x0) > x1
          ext6(i,k) = abs(a4(4,i,k)) > x1
       enddo
     endif
  enddo

!---------------------------
! Apply subgrid constraints:
!---------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
! Top 2 and bottom 2 layers always use monotonic mapping

  select case (iv)
  case (0)
     do i=i1,i2
        a4(2,i,1) = max(0., a4(2,i,1))
     enddo
  case(-1)
      do i=i1,i2
         if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
      enddo
   case(2)
     do i=i1,i2
        a4(2,i,1) = a4(1,i,1)
        a4(3,i,1) = a4(1,i,1)
        a4(4,i,1) = 0.
     enddo
  end select !iv

  if ( iv/=2 ) then
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
     select case (abs(kord))
     case (0:8)
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

    case (9)
       do i=i1,i2
          if ( extm(i,k) .and. extm(i,k-1) ) then  ! c90_mp122
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else if ( extm(i,k) .and. extm(i,k+1) ) then  ! c90_mp122
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
     case(10) !restored AM4 case 10
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
    case (11)
       do i=i1,i2
         if ( ext5(i,k) .and. (ext5(i,k-1) .or. ext5(i,k+1)) ) then
! Noisy region:
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
         else
              a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         endif
       enddo
    case (12) !post-AM4 case 10
       do i=i1,i2
          if( ext5(i,k) ) then
              if( ext5(i,k-1) .or. ext5(i,k+1) ) then
                   a4(2,i,k) = a4(1,i,k)
                   a4(3,i,k) = a4(1,i,k)
              elseif ( ext6(i,k-1) .or. ext6(i,k+1) ) then
                   pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                   lac_1 = pmp_1 + 1.5*gam(i,k+2)
                   a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                                  max(a4(1,i,k), pmp_1, lac_1) )
                   pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                   lac_2 = pmp_2 - 1.5*gam(i,k-1)
                   a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                                  max(a4(1,i,k), pmp_2, lac_2) )
              endif
          elseif( ext6(i,k) ) then
              if( ext5(i,k-1) .or. ext5(i,k+1) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
                  a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                                 max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
                  a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                                 max(a4(1,i,k), pmp_2, lac_2) )
              endif
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
    case (13)  !former 14: no subgrid limiter

       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo

    case (14) !strict monotonicity constraint
       call cs_limiters(im, extm(i1,k), a4(1,i1,k), 2)
    case (15)
       call cs_limiters(im, extm(i1,k), a4(1,i1,k), 1)
    case default
       call mpp_error(FATAL, 'kord not implemented')
    end select

! Additional constraint to ensure positivity
     if ( iv==0 .and. abs(kord) <= 13 ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 0)

  enddo      ! k-loop

!----------------------------------
! Bottom layer subgrid constraints:
!----------------------------------
  select case (iv)
  case (0)
     do i=i1,i2
        a4(3,i,km) = max(0., a4(3,i,km))
     enddo
  case (-1)
      do i=i1,i2
         if ( a4(3,i,km)*a4(1,i,km) <= 0. )  a4(3,i,km) = 0.
      enddo
   end select

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

!     if(km>8 .and. kord>4) call steepz(i1, i2, km, a4, df2, dc, delq, delp, d4)

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


!This routine should be moved to fv_io.F90.
 subroutine rst_remap(km, kn, is,ie,js,je, isd,ied,jsd,jed, nq, ntp, &
                      delp_r, u0_r, v0_r, u_r, v_r, w_r, delz_r, pt_r, q_r, qdiag_r, &
                      delp,   u0,   v0,   u,   v,   w,   delz,   pt,   q,   qdiag,   &
                      ak_r, bk_r, ptop, ak, bk, hydrostatic, make_nh, &
                      domain, square_domain, is_ideal_case)
!------------------------------------
! Assuming hybrid sigma-P coordinate:
!------------------------------------
! !INPUT PARAMETERS:
  integer, intent(in):: km                    ! Restart z-dimension
  integer, intent(in):: kn                    ! Run time dimension
  integer, intent(in):: nq, ntp               ! number of tracers (including h2o)
  integer, intent(in):: is,ie,isd,ied         ! starting & ending X-Dir index
  integer, intent(in):: js,je,jsd,jed         ! starting & ending Y-Dir index
  logical, intent(in):: hydrostatic, make_nh, square_domain, is_ideal_case
  real, intent(IN) :: ptop
  real, intent(in) :: ak_r(km+1)
  real, intent(in) :: bk_r(km+1)
  real, intent(in) :: ak(kn+1)
  real, intent(in) :: bk(kn+1)
  real, intent(in):: delp_r(is:ie,js:je,km) ! pressure thickness
  real, intent(in)::   u0_r(is:ie,  js:je+1,km)   ! initial (t=0) u-wind (m/s)
  real, intent(in)::   v0_r(is:ie+1,js:je  ,km)   ! initial (t=0) v-wind (m/s)
  real, intent(in)::   u_r(is:ie,  js:je+1,km)   ! u-wind (m/s)
  real, intent(in)::   v_r(is:ie+1,js:je  ,km)   ! v-wind (m/s)
  real, intent(inout)::  pt_r(is:ie,js:je,km)
  real, intent(in)::   w_r(is:ie,js:je,km)
  real, intent(in)::   q_r(is:ie,js:je,km,1:ntp)
  real, intent(in)::   qdiag_r(is:ie,js:je,km,ntp+1:nq)
  real, intent(inout)::delz_r(is:ie,js:je,km)
  type(domain2d), intent(INOUT) :: domain
! Output:
  real, intent(out):: delp(isd:ied,jsd:jed,kn) ! pressure thickness
  real, intent(out):: u0(isd:,jsd:,1:)   ! initial (t=0) u-wind (m/s)
  real, intent(out):: v0(isd:,jsd:,1:)   ! initial (t=0) v-wind (m/s)
  real, intent(out)::  u(isd:ied  ,jsd:jed+1,kn)   ! u-wind (m/s)
  real, intent(out)::  v(isd:ied+1,jsd:jed  ,kn)   ! v-wind (m/s)
  real, intent(out)::  w(isd:     ,jsd:     ,1:)   ! vertical velocity (m/s)
  real, intent(out):: pt(isd:ied  ,jsd:jed  ,kn)   ! temperature
  real, intent(out):: q(isd:ied,jsd:jed,kn,1:ntp)
  real, intent(out):: qdiag(isd:ied,jsd:jed,kn,ntp+1:nq)
  real, intent(out):: delz(is:,js:,1:)   ! delta-height (m)
!-----------------------------------------------------------------------
  real r_vir, rgrav
  real ps(isd:ied,jsd:jed)  ! surface pressure
  real  pe1(is:ie,km+1)
  real  pe2(is:ie,kn+1)
  real  pv1(is:ie+1,km+1)
  real  pv2(is:ie+1,kn+1)

  integer i,j,k , iq
  !CS operator replaces original mono PPM 4 --- lmh 19apr23
  integer, parameter:: kord=4 ! 13

#ifdef HYDRO_DELZ_REMAP
  if (is_master() .and. .not. hydrostatic) then
     print*, ''
     print*, ' REMAPPING IC: INITIALIZING DELZ WITH HYDROSTATIC STATE  '
     print*, ''
  endif
#endif

#ifdef HYDRO_DELZ_EXTRAP
  if (is_master() .and. .not. hydrostatic) then
     print*, ''
     print*, ' REMAPPING IC: INITIALIZING DELZ WITH HYDROSTATIC STATE ABOVE INPUT MODEL TOP  '
     print*, ''
  endif
#endif

#ifdef ZERO_W_EXTRAP
  if (is_master() .and. .not. hydrostatic) then
     print*, ''
     print*, ' REMAPPING IC: INITIALIZING W TO ZERO ABOVE INPUT MODEL TOP  '
     print*, ''
  endif
#endif

  r_vir = rvgas/rdgas - 1.
  rgrav = 1./grav

!$OMP parallel do default(none) shared(is,ie,js,je,ps,ak_r)
  do j=js,je
     do i=is,ie
        ps(i,j) = ak_r(1)
     enddo
  enddo

! this OpenMP do-loop setup cannot work in it's current form....
!$OMP parallel do default(none) shared(is,ie,js,je,km,ps,delp_r)
  do j=js,je
     do k=1,km
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
!$OMP parallel do default(none) shared(is,ie,js,je,km,pt_r,r_vir,q_r)
  do k=1,km
     do j=js,je
        do i=is,ie
           pt_r(i,j,k) = pt_r(i,j,k) * (1.+r_vir*q_r(i,j,k,1))
        enddo
     enddo
  enddo

!$OMP parallel do default(none) shared(is,ie,js,je,km,ak_r,bk_r,ps,kn,ak,bk,u0_r,u_r,u0,u,delp, &
!$OMP                                  ntp,nq,hydrostatic,make_nh,w_r,w,delz_r,delp_r,delz, &
!$OMP                                  pt_r,pt,v0_r,v_r,v0,v,q,q_r,qdiag,qdiag_r,is_ideal_case) &
!$OMP                          private(pe1,  pe2, pv1, pv2)
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

     if (is_ideal_case) then
        call remap_2d(km, pe1, u0_r(is:ie,j:j,1:km),      &
                      kn, pe2,   u0(is:ie,j:j,1:kn),      &
                      is, ie, -1, kord)
     endif

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
          do iq=1,ntp
             call remap_2d(km, pe1, q_r(is:ie,j:j,1:km,iq:iq),  &
                           kn, pe2,   q(is:ie,j:j,1:kn,iq:iq),  &
                           is, ie, 0, kord)
          enddo
          do iq=ntp+1,nq
             call remap_2d(km, pe1, qdiag_r(is:ie,j:j,1:km,iq:iq),  &
                           kn, pe2,   qdiag(is:ie,j:j,1:kn,iq:iq),  &
                           is, ie, 0, kord)
          enddo
      endif

      if ( .not. hydrostatic .and. .not. make_nh) then
! Remap vertical wind:
         call remap_2d(km, pe1, w_r(is:ie,j:j,1:km),       &
                       kn, pe2,   w(is:ie,j:j,1:kn),       &
                       is, ie, -1, kord)

#ifdef ZERO_W_EXTRAP
       do k=1,kn
       do i=is,ie
          if (pe2(i,k) < pe1(i,1)) then
             w(i,j,k) = 0.
          endif
       enddo
       enddo
#endif

#ifndef HYDRO_DELZ_REMAP
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
#endif
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

#ifdef HYDRO_DELZ_REMAP
       !initialize delz from the hydrostatic state
       do k=1,kn
       do i=is,ie
          delz(i,j,k) = (rdgas*rgrav)*pt(i,j,k)*(pe2(i,k)-pe2(i,k+1))
       enddo
       enddo
#endif
#ifdef HYDRO_DELZ_EXTRAP
       !initialize delz from the hydrostatic state
       do k=1,kn
       do i=is,ie
          if (pe2(i,k) < pe1(i,1)) then
             delz(i,j,k) = (rdgas*rgrav)*pt(i,j,k)*(pe2(i,k)-pe2(i,k+1))
          endif
       enddo
       enddo
#endif
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

       if (is_ideal_case) then
          call remap_2d(km, pv1, v0_r(is:ie+1,j:j,1:km),      &
                        kn, pv2,   v0(is:ie+1,j:j,1:kn),      &
                        is, ie+1, -1, kord)
       endif

       call remap_2d(km, pv1, v_r(is:ie+1,j:j,1:km),       &
                     kn, pv2,   v(is:ie+1,j:j,1:kn),       &
                     is, ie+1, -1, kord)

  endif !(j < je+1)
1000  continue

!$OMP parallel do default(none) shared(is,ie,js,je,kn,pt,r_vir,q)
  do k=1,kn
     do j=js,je
        do i=is,ie
           pt(i,j,k) = pt(i,j,k) / (1.+r_vir*q(i,j,k,1))
        enddo
     enddo
  enddo

 end subroutine rst_remap


 !This routine is indended to remap between different #
 !     of vertical levels
 subroutine mappm(km, pe1, q1, kn, pe2, q2, i1, i2, iv, kord)

! IV = 0: constituents
! IV = 1: potential temp
! IV =-1: winds
! IV =-2: vertical velocity

! Mass flux preserving mapping: q1(im,km) -> q2(im,kn)

! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate

 integer, intent(in):: i1, i2, km, kn, kord, iv
 real, intent(in ):: pe1(i1:i2,km+1), pe2(i1:i2,kn+1)
 real, intent(in )::  q1(i1:i2,km) ! input field
 real, intent(out)::  q2(i1:i2,kn) ! output field

! local
      real  qs(i1:i2)
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
         call  cs_profile( qs, a4, dp1, km, i1, i2, iv, kord )
      else
           call ppm_profile( a4, dp1, km, i1, i2, iv, kord )
      endif

      do 5555 i=i1,i2
         k0 = 1
      do 555 k=1,kn

         if(pe2(i,k) .le. pe1(i,1)) then
! above old ptop
            q2(i,k) = q1(i,1)
         elseif(pe2(i,k) .ge. pe1(i,km+1)) then
! Entire grid below old ps
            q2(i,k) = q1(i,km)
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
           qsum = qsum + delp * q1(i,km)
          dpsum = dpsum + delp
        endif
123     q2(i,k) = qsum / dpsum
      endif
555   continue
5555  continue

 end subroutine mappm


 subroutine moist_cv(is,ie, isd,ied, jsd,jed, km, j, k, nwat, sphum, liq_wat, rainwat,    &
                     ice_wat, snowwat, graupel, q, qd, cvm, t1)
  integer, intent(in):: is, ie, isd,ied, jsd,jed, km, nwat, j, k
  integer, intent(in):: sphum, liq_wat, rainwat, ice_wat, snowwat, graupel
  real, intent(in), dimension(isd:ied,jsd:jed,km,nwat):: q
  real, intent(out), dimension(is:ie):: cvm, qd  ! qd is q_con
  real, intent(in), optional:: t1(is:ie)
!
  real, parameter:: t_i0 = 15.
  real, dimension(is:ie):: qv, ql, qs
  integer:: i

  select case (nwat)

   case(2)
     if ( present(t1) ) then  ! Special case for GFS physics
        do i=is,ie
           qd(i) = max(0., q(i,j,k,liq_wat))
           if ( t1(i) > tice ) then
                qs(i) = 0.
           elseif ( t1(i) < tice-t_i0 ) then
                qs(i) = qd(i)
           else
                qs(i) = qd(i)*(tice-t1(i))/t_i0
           endif
           ql(i) = qd(i) - qs(i)
           qv(i) = max(0.,q(i,j,k,sphum))
           cvm(i) = (1.-(qv(i)+qd(i)))*cv_air + qv(i)*cv_vap + ql(i)*c_liq + qs(i)*c_ice
        enddo
     else
        do i=is,ie
           qv(i) = max(0.,q(i,j,k,sphum))
           qs(i) = max(0.,q(i,j,k,liq_wat))
           qd(i) = qs(i)
           cvm(i) = (1.-qv(i))*cv_air + qv(i)*cv_vap
        enddo
     endif
  case (3)
     do i=is,ie
        qv(i) = q(i,j,k,sphum)
        ql(i) = q(i,j,k,liq_wat)
        qs(i) = q(i,j,k,ice_wat)
        qd(i) = ql(i) + qs(i)
        cvm(i) = (1.-(qv(i)+qd(i)))*cv_air + qv(i)*cv_vap + ql(i)*c_liq + qs(i)*c_ice
     enddo
  case(4)              ! K_warm_rain with fake ice
     do i=is,ie
        qv(i) = q(i,j,k,sphum)
        qd(i) = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
        cvm(i) = (1.-(qv(i)+qd(i)))*cv_air + qv(i)*cv_vap + qd(i)*c_liq
     enddo
  case(5)
     do i=is,ie
        qv(i) = q(i,j,k,sphum)
        ql(i) = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
        qs(i) = q(i,j,k,ice_wat) + q(i,j,k,snowwat)
        qd(i) = ql(i) + qs(i)
        cvm(i) = (1.-(qv(i)+qd(i)))*cv_air + qv(i)*cv_vap + ql(i)*c_liq + qs(i)*c_ice
     enddo
  case(6)
     do i=is,ie
        qv(i) = q(i,j,k,sphum)
        ql(i) = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
        qs(i) = q(i,j,k,ice_wat) + q(i,j,k,snowwat) + q(i,j,k,graupel)
        qd(i) = ql(i) + qs(i)
        cvm(i) = (1.-(qv(i)+qd(i)))*cv_air + qv(i)*cv_vap + ql(i)*c_liq + qs(i)*c_ice
     enddo
  case default
     !call mpp_error (NOTE, 'fv_mapz::moist_cv - using default cv_air')
     do i=is,ie
         qd(i) = 0.
        cvm(i) = cv_air
     enddo
 end select

 end subroutine moist_cv

 subroutine moist_cp(is,ie, isd,ied, jsd,jed, km, j, k, nwat, sphum, liq_wat, rainwat,    &
                     ice_wat, snowwat, graupel, q, qd, cpm, t1)

  integer, intent(in):: is, ie, isd,ied, jsd,jed, km, nwat, j, k
  integer, intent(in):: sphum, liq_wat, rainwat, ice_wat, snowwat, graupel
  real, intent(in), dimension(isd:ied,jsd:jed,km,nwat):: q
  real, intent(out), dimension(is:ie):: cpm, qd
  real, intent(in), optional:: t1(is:ie)
!
  real, parameter:: t_i0 = 15.
  real, dimension(is:ie):: qv, ql, qs
  integer:: i

  select case (nwat)

  case(2)
     if ( present(t1) ) then  ! Special case for GFS physics
        do i=is,ie
           qd(i) = max(0., q(i,j,k,liq_wat))
           if ( t1(i) > tice ) then
                qs(i) = 0.
           elseif ( t1(i) < tice-t_i0 ) then
                qs(i) = qd(i)
           else
                qs(i) = qd(i)*(tice-t1(i))/t_i0
           endif
           ql(i) = qd(i) - qs(i)
           qv(i) = max(0.,q(i,j,k,sphum))
           cpm(i) = (1.-(qv(i)+qd(i)))*cp_air + qv(i)*cp_vapor + ql(i)*c_liq + qs(i)*c_ice
        enddo
     else
     do i=is,ie
        qv(i) = max(0.,q(i,j,k,sphum))
        qs(i) = max(0.,q(i,j,k,liq_wat))
        qd(i) = qs(i)
        cpm(i) = (1.-qv(i))*cp_air + qv(i)*cp_vapor
     enddo
     endif

  case(3)
     do i=is,ie
        qv(i) = q(i,j,k,sphum)
        ql(i) = q(i,j,k,liq_wat)
        qs(i) = q(i,j,k,ice_wat)
        qd(i) = ql(i) + qs(i)
        cpm(i) = (1.-(qv(i)+qd(i)))*cp_air + qv(i)*cp_vapor + ql(i)*c_liq + qs(i)*c_ice
     enddo
  case(4)    ! K_warm_rain scheme with fake ice
     do i=is,ie
        qv(i) = q(i,j,k,sphum)
        qd(i) = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
        cpm(i) = (1.-(qv(i)+qd(i)))*cp_air + qv(i)*cp_vapor + qd(i)*c_liq
     enddo
  case(5)
     do i=is,ie
        qv(i) = q(i,j,k,sphum)
        ql(i) = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
        qs(i) = q(i,j,k,ice_wat) + q(i,j,k,snowwat)
        qd(i) = ql(i) + qs(i)
        cpm(i) = (1.-(qv(i)+qd(i)))*cp_air + qv(i)*cp_vapor + ql(i)*c_liq + qs(i)*c_ice
     enddo
  case(6)
     do i=is,ie
        qv(i) = q(i,j,k,sphum)
        ql(i) = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
        qs(i) = q(i,j,k,ice_wat) + q(i,j,k,snowwat) + q(i,j,k,graupel)
        qd(i) = ql(i) + qs(i)
        cpm(i) = (1.-(qv(i)+qd(i)))*cp_air + qv(i)*cp_vapor + ql(i)*c_liq + qs(i)*c_ice
     enddo
  case default
     !call mpp_error (NOTE, 'fv_mapz::moist_cp - using default cp_air')
     do i=is,ie
        qd(i) = 0.
        cpm(i) = cp_air
     enddo
  end select

 end subroutine moist_cp
!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  map1_cubic --- Cubic Interpolation for vertical re-mapping
!
! !INTERFACE:
  subroutine map1_cubic( km,   pe1,    q1,                 &
                         kn,   pe2,    q2,   i1, i2,       &
                         j,    ibeg, iend, jbeg, jend, akap, T_VAR, conserv)
      implicit none

! !INPUT PARAMETERS:
      integer, intent(in) :: i1                ! Starting longitude
      integer, intent(in) :: i2                ! Finishing longitude
      real, intent(in) :: akap
      integer, intent(in) :: T_VAR             ! Thermodynamic variable to remap
                                               !     1:TE  2:T  3:PT
      logical, intent(in) :: conserv
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
!     Perform Cubic Interpolation a given latitude
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
!
! !REVISION HISTORY:
!    2005.11.14   Takacs    Initial Code
!    2016.07.20   Putman    Modified to make genaric for any thermodynamic variable
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      real       qx(i1:i2,km)
      real   logpl1(i1:i2,km)
      real   logpl2(i1:i2,kn)
      real   dlogp1(i1:i2,km)
      real    vsum1(i1:i2)
      real    vsum2(i1:i2)
      real   am2,am1,ap0,ap1,P,PLP1,PLP0,PLM1,PLM2,DLP0,DLM1,DLM2

      integer i, k, LM2,LM1,LP0,LP1

! Initialization
! --------------

      select case (T_VAR)
      case(1)
       ! Total Energy Remapping in Log(P)
        do k=1,km
            qx(:,k) = q1(i1:i2,j,k)
        logpl1(:,k) = log( 0.5*(pe1(:,k)+pe1(:,k+1)) )
        enddo
        do k=1,kn
        logpl2(:,k) = log( 0.5*(pe2(:,k)+pe2(:,k+1)) )
        enddo

        do k=1,km-1
        dlogp1(:,k) = logpl1(:,k+1)-logpl1(:,k)
        enddo

      case(2)
       ! Temperature Remapping in Log(P)
        do k=1,km
            qx(:,k) = q1(i1:i2,j,k)
        logpl1(:,k) = log( 0.5*(pe1(:,k)+pe1(:,k+1)) )
        enddo
        do k=1,kn
        logpl2(:,k) = log( 0.5*(pe2(:,k)+pe2(:,k+1)) )
        enddo

        do k=1,km-1
        dlogp1(:,k) = logpl1(:,k+1)-logpl1(:,k)
        enddo

      case(3)
       ! Potential Temperature Remapping in P^KAPPA
        do k=1,km
            qx(:,k) = q1(i1:i2,j,k)
        logpl1(:,k) = exp( akap*log( 0.5*(pe1(:,k)+pe1(:,k+1))) )
        enddo
        do k=1,kn
        logpl2(:,k) = exp( akap*log( 0.5*(pe2(:,k)+pe2(:,k+1))) )
        enddo

        do k=1,km-1
        dlogp1(:,k) = logpl1(:,k+1)-logpl1(:,k)
        enddo

      end select

      if (conserv) then
! Compute vertical integral of Input TE
! -------------------------------------
        vsum1(:) = 0.0
        do i=i1,i2
        do k=1,km
        vsum1(i) = vsum1(i) + qx(i,k)*( pe1(i,k+1)-pe1(i,k) )
        enddo
        vsum1(i) = vsum1(i) / ( pe1(i,km+1)-pe1(i,1) )
        enddo

      endif

! Interpolate TE onto target Pressures
! ------------------------------------
      do i=i1,i2
      do k=1,kn
         LM1 = 1
         LP0 = 1
         do while( LP0.le.km )
            if (logpl1(i,LP0).lt.logpl2(i,k)) then
               LP0 = LP0+1
            else
               exit
            endif
         enddo
         LM1 = max(LP0-1,1)
         LP0 = min(LP0, km)

! Extrapolate Linearly in LogP above first model level
! ----------------------------------------------------
         if( LM1.eq.1 .and. LP0.eq.1 ) then
             q2(i,j,k) = qx(i,1) + ( qx(i,2)-qx(i,1) )*( logpl2(i,k)-logpl1(i,1) ) &
                                                      /( logpl1(i,2)-logpl1(i,1) )

! Extrapolate Linearly in LogP below last model level
! ---------------------------------------------------
         else if( LM1.eq.km .and. LP0.eq.km ) then
             q2(i,j,k) = qx(i,km) + ( qx(i,km)-qx(i,km-1) )*( logpl2(i,k )-logpl1(i,km  ) ) &
                                                           /( logpl1(i,km)-logpl1(i,km-1) )

! Interpolate Linearly in LogP between levels 1 => 2 and km-1 => km
! -----------------------------------------------------------------
         else if( LM1.eq.1 .or. LP0.eq.km ) then
             q2(i,j,k) = qx(i,LP0) + ( qx(i,LM1)-qx(i,LP0) )*( logpl2(i,k  )-logpl1(i,LP0) ) &
                                                            /( logpl1(i,LM1)-logpl1(i,LP0) )
! Interpolate Cubicly in LogP between other model levels
! ------------------------------------------------------
         else
              LP1 = LP0+1
              LM2 = LM1-1
             P    = logpl2(i,k)
             PLP1 = logpl1(i,LP1)
             PLP0 = logpl1(i,LP0)
             PLM1 = logpl1(i,LM1)
             PLM2 = logpl1(i,LM2)
             DLP0 = dlogp1(i,LP0)
             DLM1 = dlogp1(i,LM1)
             DLM2 = dlogp1(i,LM2)

              ap1 = (P-PLP0)*(P-PLM1)*(P-PLM2)/( DLP0*(DLP0+DLM1)*(DLP0+DLM1+DLM2) )
              ap0 = (PLP1-P)*(P-PLM1)*(P-PLM2)/( DLP0*      DLM1 *(     DLM1+DLM2) )
              am1 = (PLP1-P)*(PLP0-P)*(P-PLM2)/( DLM1*      DLM2 *(DLP0+DLM1     ) )
              am2 = (PLP1-P)*(PLP0-P)*(PLM1-P)/( DLM2*(DLM1+DLM2)*(DLP0+DLM1+DLM2) )

             q2(i,j,k) = ap1*qx(i,LP1) + ap0*qx(i,LP0) + am1*qx(i,LM1) + am2*qx(i,LM2)

         endif

      enddo
      enddo
      if (conserv) then

! Compute vertical integral of Output TE
! --------------------------------------
        vsum2(:) = 0.0
        do i=i1,i2
        do k=1,kn
        vsum2(i) = vsum2(i) + q2(i,j,k)*( pe2(i,k+1)-pe2(i,k) )
        enddo
        vsum2(i) = vsum2(i) / ( pe2(i,kn+1)-pe2(i,1) )
        enddo

! Adjust Final TE to conserve
! ---------------------------
        do i=i1,i2
        do k=1,kn
           q2(i,j,k) = q2(i,j,k) + vsum1(i)-vsum2(i)
!          q2(i,j,k) = q2(i,j,k) * vsum1(i)/vsum2(i)
        enddo
        enddo

      endif

      return
!EOC
 end subroutine map1_cubic
end module fv_mapz_mod
