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

! Linjiong Zhou: Nov 19, 2019
! Revise the OpenMP code to avoid crash
module fv_mapz_mod

#ifdef OVERLOAD_R4
  use constantsR4_mod,   only: pi=>pi_8, rdgas, grav, cp_air, cp_vapor
#else
  use constants_mod,     only: pi=>pi_8, rdgas, grav, cp_air, cp_vapor
#endif
  use fv_arrays_mod,     only: radius ! scaled for small earth
  use tracer_manager_mod,only: get_tracer_index
  use field_manager_mod, only: MODEL_ATMOS
  use fv_grid_utils_mod, only: g_sum, ptop_min, cubed_to_latlon
  use fv_fill_mod,       only: fillz
  use mpp_domains_mod,   only: domain2d
  use mpp_mod,           only: FATAL, NOTE, mpp_error
  use fv_arrays_mod,     only: fv_grid_type, fv_grid_bounds_type, R_GRID, inline_mp_type
  use fv_timing_mod,     only: timing_on, timing_off
  use intermediate_phys_mod, only: intermediate_phys
  use fv_operators_mod,  only: map_scalar, map1_ppm, mapn_tracer, map1_q2, map1_cubic
  use fv_thermodynamics_mod, only: moist_cv, fv_thermo_type

  implicit none
  real, parameter:: consv_min = 0.001   ! below which no correction applies
  real, parameter:: t_min= 184.   ! below which applies stricter constraint
  real, parameter:: cv_air =  cp_air - rdgas ! = rdgas * (7/2-1) = 2.5*rdgas=717.68

  real(kind=4) :: E_Flux = 0.
  private

  public Lagrangian_to_Eulerian, E_Flux, consv_min

contains

 subroutine Lagrangian_to_Eulerian(last_step, consv, ps, pe, delp, pkz, pk,   &
                                   mdt, pdt, npx, npy, km, is,ie,js,je, isd,ied,jsd,jed,       &
                      nq, nwat, sphum, q_con, u, v, w, delz, pt, q, hs, r_vir, cp, te_err, tw_err, &
                      akap, cappa, kord_mt, kord_wz, kord_tr, kord_tm, remap_te, peln, te0_2d,        &
                      ng, ua, va, omga, te, ws, fill, reproduce_sum,      &
                      ptop, ak, bk, pfull, gridstruct, thermostruct, domain, do_sat_adj, &
                      hydrostatic, hybrid_z, adiabatic, do_adiabatic_init, &
                      do_inline_mp, inline_mp, bd, fv_debug, &
                      do_fast_phys, do_intermediate_phys, consv_checker, adj_mass_vmr)

  logical, intent(in):: last_step
  logical, intent(in):: fv_debug
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
  type(fv_thermo_type), intent(IN), target :: thermostruct
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
  integer:: nt, liq_wat, ice_wat, rainwat, snowwat, cld_amt, graupel, iq, n, kmp, kp, k_next
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
  ccn_cm3 = get_tracer_index (MODEL_ATMOS, 'ccn_cm3')
  cin_cm3 = get_tracer_index (MODEL_ATMOS, 'cin_cm3')
  aerosol = get_tracer_index (MODEL_ATMOS, 'aerosol')

!$OMP parallel do default(none) shared(is,ie,js,je,km,pe,ptop,kord_tm,hydrostatic, &
!$OMP                                  pt,pk,rg,peln,q,nwat,liq_wat,rainwat,ice_wat,snowwat,    &
!$OMP                                  graupel,q_con,sphum,cappa,r_vir,k1k,delp, &
!$OMP                                  delz,akap,pkz,te,u,v,ps, gridstruct, last_step, &
!$OMP                                  ak,bk,nq,isd,ied,jsd,jed,kord_tr,fill, adiabatic, &
!$OMP                                  hs,w,ws,kord_wz,omga,rrg,kord_mt,pe4,cp,remap_te,thermostruct)&
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
                 ! Transform "density pt" to "density temp". (OK to have flag in outer loop)
                 do k=1,km
                    if (thermostruct%moist_kappa) then
                       call moist_cv(is,ie,isd,ied,jsd,jed, km, j, k, nwat, sphum, liq_wat, rainwat,    &
                            ice_wat, snowwat, graupel, q, gz(is:ie), cvm(is:ie))
                       do i=is,ie
                          q_con(i,j,k) = gz(i)
                          cappa(i,j,k) = rdgas / ( rdgas + cvm(i)/(1.+r_vir*q(i,j,k,sphum)) )
                          pt(i,j,k) = pt(i,j,k)*exp(cappa(i,j,k)/(1.-cappa(i,j,k))*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
                       enddo
                    else
                       do i=is,ie
                          pt(i,j,k) = pt(i,j,k)*exp(k1k*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
                       enddo
                    endif !moist_kappa
                 enddo
              endif         ! hydro test


           endif !kord_tm
        else
           !----------------------------------
           ! remap_te: Compute cp*T + KE +phis
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
                 if (thermostruct%moist_kappa) then
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
                 else
                    do i=is,ie
                       pkz(i,j,k) = exp(k1k*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
                       te(i,j,k) = cv_air*pt(i,j,k)*pkz(i,j,k)/(1.+r_vir*q(i,j,k,sphum)) +     &
                            0.5 * w(i,j,k)**2 + 0.25*gridstruct%rsin2(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                            v(i,j,k)**2+v(i+1,j,k)**2 -  &
                            (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j)) +         &
                            0.5*(phis(i,k+1)+phis(i,k))
                    enddo
                 endif
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
                  pkz(i,j,k) = (pk2(i,k+1)-pk2(i,k))/(akap*(peln(i,k+1,j)-peln(i,k,j))) !MOIST_CAPPA not supported
               enddo
            enddo
         else
            ! Note: pt at this stage is T_v or T_m , unless kord_tm > 0
            do k=1,km
               if (thermostruct%moist_kappa) then
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
               else
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
               endif !moist_kappa
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
         if ( hydrostatic ) then !MOIST_CAPPA not supported
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
               if (thermostruct%moist_kappa) then
                  call moist_cv(is,ie,isd,ied,jsd,jed, km, j, k, nwat, sphum, liq_wat, rainwat,    &
                       ice_wat, snowwat, graupel, q, gz(is:ie), cvm(is:ie))
                  do i=is,ie
                     q_con(i,j,k) = gz(i)
                     cappa(i,j,k) = rdgas / ( rdgas + cvm(i)/(1.+r_vir*q(i,j,k,sphum)) )
                     phis(i,k) = phis(i,k+1) - delz(i,j,k)*grav
                     tpe = te(i,j,k) - 0.5*(phis(i,k)+phis(i,k+1)) - 0.5*w(i,j,k)**2 - 0.25*gridstruct%rsin2(i,j)*(    &
                          u(i,j,k)**2+u(i,j+1,k)**2 + v(i,j,k)**2+v(i+1,j,k)**2 -  &
                          (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j) )
                     pt(i,j,k)= tpe / cvm(i)*(1.+r_vir*q(i,j,k,sphum))*(1.-gz(i))
                     pkz(i,j,k) = exp(cappa(i,j,k)*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
                  enddo
               else
                  do i=is,ie
                     phis(i,k) = phis(i,k+1) - delz(i,j,k)*grav
                     tpe = te(i,j,k) - 0.5*(phis(i,k)+phis(i,k+1)) - 0.5*w(i,j,k)**2 - 0.25*gridstruct%rsin2(i,j)*(    &
                          u(i,j,k)**2+u(i,j+1,k)**2 + v(i,j,k)**2+v(i+1,j,k)**2 -  &
                          (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j) )
                     pt(i,j,k)= tpe / cv_air *(1.+r_vir*q(i,j,k,sphum))
                     pkz(i,j,k) = exp(akap*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
                  enddo
               endif

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
!$OMP                                  zsum0,te0_2d,gridstruct,q,kord_tm,te,remap_te,thermostruct) &
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
                    if (thermostruct%use_cond) then
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
                    else
                       do i=is,ie
                          te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(cv_air*pt(i,j,k)/(1.+r_vir*q(i,j,k,sphum)) + &
                               0.5*(phis(i,k)+phis(i,k+1) + w(i,j,k)**2 + 0.5*gridstruct%rsin2(i,j)*( &
                               u(i,j,k)**2+u(i,j+1,k)**2 + v(i,j,k)**2+v(i+1,j,k)**2 -  &
                               (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j))))
                       enddo
                    endif
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
                 mdt, consv, akap, ptop, pfull, hs, te0_2d, u, &
                 v, w, pt, delp, delz, q_con, cappa, q, pkz, r_vir, te_err, tw_err, &
                 inline_mp, gridstruct, thermostruct, domain, bd, hydrostatic, do_adiabatic_init, &
                 do_inline_mp, do_sat_adj, last_step, do_fast_phys, consv_checker, adj_mass_vmr)
        call timing_off('INTERMEDIATE_PHYS')
    endif

!-----------------------------------------------------------------------
! <<< Intermediate Physics
!-----------------------------------------------------------------------

  if ( last_step ) then
       ! 9a) Convert T_v/T_m to T if last_step
!$OMP parallel do default(none) shared(is,ie,js,je,km,isd,ied,jsd,jed,hydrostatic,pt,adiabatic,cp, &
!$OMP                                  nwat,rainwat,liq_wat,ice_wat,snowwat,graupel,r_vir,&
!$OMP                                  sphum,pkz,dtmp,q,thermostruct) &
!$OMP                          private(cvm,gz)
     do k=1,km
        do j=js,je
           if (hydrostatic) then !This is re-factored from AM4 so answers may be different
              do i=is,ie
                 pt(i,j,k) = (pt(i,j,k)+dtmp/cp*pkz(i,j,k)) / (1.+r_vir*q(i,j,k,sphum)) !use_cond not implemented
              enddo
           else
              if (thermostruct%use_cond) then
                 call moist_cv(is,ie,isd,ied,jsd,jed, km, j, k, nwat, sphum, liq_wat, rainwat,    &
                      ice_wat, snowwat, graupel, q, gz(is:ie), cvm(is:ie)) !gz is q_con
                 do i=is,ie
                    pt(i,j,k) = (pt(i,j,k)+dtmp/cvm(i)*pkz(i,j,k))/((1.+r_vir*q(i,j,k,sphum))*(1.-gz(i)))
                 enddo
              else
                 if ( .not. adiabatic ) then
                    do i=is,ie
                       pt(i,j,k) = (pt(i,j,k)+dtmp/cv_air*pkz(i,j,k)) / (1.+r_vir*q(i,j,k,sphum))
                    enddo
                 endif
              endif
           endif
        enddo   ! j-loop
     enddo  ! k-loop
!!$     if (hydrostatic) then
!!$        thermostruct%pt_is_virtual = .false.
!!$     else
!!$        if (thermostruct%use_cond) then
!!$           thermostruct%pt_is_virtual = .false.
!!$           thermostruct%pt_is_density = .false.
!!$        else
!!$           if (.not. adiabatic) thermostruct%pt_is_virtual = .false.
!!$        endif
!!$     endif
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
     !thermostruct%pt_is_potential = .true.
  endif

 end subroutine Lagrangian_to_Eulerian




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

end module fv_mapz_mod
