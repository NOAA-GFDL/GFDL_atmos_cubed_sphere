module fv_update_phys_mod

  use constants_mod,      only: kappa, rdgas, grav
  use fv_control_mod,     only: npx, npy, npz, ncnst, ch4_chem, k_top, nwat, fv_debug
  use field_manager_mod,  only: MODEL_ATMOS
  use tracer_manager_mod, only: get_tracer_index
  use time_manager_mod,   only: time_type
  use fv_mp_mod,          only: domain, gid
  use fv_eta_mod,         only: get_eta_level
  use mpp_domains_mod,    only: mpp_update_domains
  use mpp_mod,            only: FATAL, mpp_error
  use fv_grid_utils_mod,         only: edge_vect_s,edge_vect_n,edge_vect_w,edge_vect_e,    &
                                es, ew, vlon, vlat
  use fv_grid_tools_mod,         only: grid_type
  use fv_timing_mod,       only: timing_on, timing_off
  use fv_diagnostics_mod,  only: prt_maxmin
#ifdef GFDL_NUDGE
  use atmos_nudge_mod,    only: get_atmos_nudge, do_ps
#endif

  implicit none

  public :: fv_update_phys

  contains

  subroutine fv_update_phys ( dt, is, ie, js, je, isd, ied, jsd, jed, ng, nq,     &
                              u, v, delp, pt, q, ua, va, ps, pe,  peln, pk, pkz,  &
                              ak, bk, u_dt, v_dt, t_dt, q_dt, u_srf, v_srf,       &
                              delz, hydrostatic, full_phys, Time, nudge )
    real, intent(in)   :: dt
    integer, intent(in):: is,  ie,  js,  je, ng
    integer, intent(in):: isd, ied, jsd, jed
    integer, intent(in):: nq            ! tracers modified by physics 
                                        ! ncnst is the total nmber of tracers
    logical, intent(in):: full_phys     ! full physics
    type (time_type), intent(in) :: Time

    logical, intent(in), optional:: nudge

    real, intent(in), dimension(npz+1):: ak, bk
    logical, intent(in):: hydrostatic
    real, intent(inout):: delz(is:ie,js:je,npz)
! Tendencies from Physics:
    real, intent(inout), dimension(isd:ied,jsd:jed,npz):: u_dt, v_dt
    real, intent(inout):: t_dt(is:ie,js:je,npz)
    real, intent(inout):: q_dt(is:ie,js:je,npz,nq)

! Saved Bottom winds for GFDL Physics Interface
    real, intent(out), dimension(is:ie,js:je):: u_srf, v_srf

    real, intent(inout):: u(isd:ied  ,jsd:jed+1,npz)  ! D grid zonal wind (m/s)
    real, intent(inout):: v(isd:ied+1,jsd:jed  ,npz)  ! D grid meridional wind (m/s)
    real, intent(inout), dimension(isd:ied,jsd:jed,npz):: pt, delp
    real, intent(inout):: q(isd:ied,jsd:jed,npz, ncnst)   ! specific humidity and constituents

!-----------------------------------------------------------------------
! Auxilliary pressure arrays:    
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
    real, intent(inout):: ps  (isd:ied  ,jsd:jed)           ! Surface pressure (pascal)
    real, intent(inout):: pe  (is-1:ie+1, npz+1,js-1:je+1)  ! edge pressure (pascal)
    real, intent(inout):: pk  (is:ie,js:je  , npz+1)        ! pe**cappa
    real, intent(inout):: peln(is:ie,npz+1,js:je)           ! ln(pe)
    real, intent(inout):: pkz (is:ie,js:je,npz)             ! finite-volume mean pk

! Winds on lat-lon grid:
    real, intent(inout), dimension(isd:ied,jsd:jed,npz):: ua, va

!***********
! Haloe Data
!***********
    real, parameter::    q1_h2o = 2.2E-6
    real, parameter::    q7_h2o = 3.8E-6
    real, parameter::  q100_h2o = 3.8E-6
    real, parameter:: q1000_h2o = 3.1E-6
    real, parameter:: q2000_h2o = 2.8E-6
    real, parameter::   tau_h2o = 120.*86400.

! Local arrays:
    real  ps_dt(is:ie,js:je)
    real  phalf(npz+1), pfull(npz)

    integer  i, j, k, m
    integer  sphum, liq_wat, ice_wat, cld_amt   ! GFDL AM physics
    integer  rainwat, snowwat, graupel          ! Lin Micro-physics
    real   qstar, dbk, dt5, rdt, rdg


    rdg = -rdgas / grav

    dt5 = 0.5 * dt
    rdt = 1./ dt

   sphum   = 1

   if ( full_phys ) then
        cld_amt = get_tracer_index (MODEL_ATMOS, 'cld_amt')
   else
        cld_amt = 7
   endif

   if ( nwat>=3 ) then
        sphum   = get_tracer_index (MODEL_ATMOS, 'sphum')
        liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
        ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
   endif

   if ( nwat==6 ) then
! Micro-physics:
        rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat')
        snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat')
        graupel = get_tracer_index (MODEL_ATMOS, 'graupel')
        if ( cld_amt<7 ) call mpp_error(FATAL,'Cloud Fraction allocation error') 
   endif

   if ( fv_debug ) then
       if ( gid==0 ) write(*,*) nq, nwat, sphum, liq_wat, ice_wat, rainwat, snowwat, graupel
       call prt_maxmin('ptop_b_update', pe(is:ie,1,js:je), is, ie, js, je, 0, 1, 0.01, gid==0)
       call prt_maxmin('delp_b_update', delp, is, ie, js,  je, ng, npz, 0.01, gid==0)
       do m=1,nq
          call prt_maxmin('q_dt', q_dt(is,js,1,m), is, ie, js, je, 0, npz, 1., gid==0)
       enddo
!         q(:,:,:,graupel) = 0.
!      q_dt(:,:,:,graupel) = 0.
   endif


#ifndef MARS_MGCM
    call get_eta_level(npz, 1.0E5, pfull, phalf, ak, bk)
#endif

!$omp parallel do private (i,j,k,m, qstar)
    do 1000 k=1, npz

! Do idealized Ch4 chemistry
       if ( ch4_chem .and. pfull(k) < 50.E2 ) then

           if ( pfull(k) < 1. ) then
               qstar = q1_h2o
           elseif ( pfull(k) <   7. .and. pfull(k) >=    1. ) then
               qstar = q1_h2o + (q7_h2o-q1_h2o)*log(pfull(k)/1.)/log(7.)
           elseif ( pfull(k) <  100. .and. pfull(k) >=    7. ) then
               qstar = q7_h2o + (q100_h2o-q7_h2o)*log(pfull(k)/7.)/log(100./7.)
           elseif ( pfull(k) < 1000. .and. pfull(k) >=  100. ) then
               qstar = q100_h2o + (q1000_h2o-q100_h2o)*log(pfull(k)/1.E2)/log(10.)
           elseif ( pfull(k) < 2000. .and. pfull(k) >= 1000. ) then
               qstar = q1000_h2o + (q2000_h2o-q1000_h2o)*log(pfull(k)/1.E3)/log(2.)
           endif

           do j=js,je
              do i=is,ie
                 q_dt(i,j,k,sphum) = q_dt(i,j,k,sphum) + (qstar-q(i,j,k,sphum))/tau_h2o
              enddo
           enddo
       endif

       do j=js,je
          do i=is,ie
             ua(i,j,k) = ua(i,j,k) + dt*u_dt(i,j,k)
             va(i,j,k) = va(i,j,k) + dt*v_dt(i,j,k)
          enddo
       enddo
       if ( hydrostatic ) then
          do j=js,je
             do i=is,ie
                pt(i,j,k) = pt(i,j,k) + dt*t_dt(i,j,k)
             enddo
          enddo
       else
! Heating/cooling from physics is assumed to be hydrostatic (constant-pressure)
! i.e., del-T = del-Q / Cp; no need for this if del-T = del-Q / Cv
          do j=js,je
             do i=is,ie
!                 pt(i,j,k) =   pt(i,j,k) + dt*t_dt(i,j,k)
                delz(i,j,k) = delz(i,j,k) / pt(i,j,k)
                  pt(i,j,k) =   pt(i,j,k) + dt*t_dt(i,j,k)
                delz(i,j,k) = delz(i,j,k) * pt(i,j,k)
             enddo
          enddo
       endif

!----------------
! Update tracers:
!----------------
       do m=1,nq
          do j=js,je
             do i=is,ie
                q(i,j,k,m) = q(i,j,k,m) + dt*q_dt(i,j,k,m)
             enddo
          enddo
       enddo

! ncnst == 6
   if ( nwat==6 ) then
        do j=js,je
           do i=is,ie
              ps_dt(i,j)  = 1. + dt * ( q_dt(i,j,k,sphum  ) +    &
                                        q_dt(i,j,k,liq_wat) +    &
                                        q_dt(i,j,k,rainwat) +    &
                                        q_dt(i,j,k,ice_wat) +    &
                                        q_dt(i,j,k,snowwat) +    &
                                        q_dt(i,j,k,graupel) )
              delp(i,j,k) = delp(i,j,k) * ps_dt(i,j)
           enddo
        enddo
   elseif( nwat==3 ) then
! GFDL AM2/3 phys:
        do j=js,je
           do i=is,ie
!--------------------------------------------------------
! Adjust total air mass due to changes in water substance
!--------------------------------------------------------
               ps_dt(i,j) = 1. + dt*(q_dt(i,j,k,sphum  ) +    &
                                     q_dt(i,j,k,liq_wat) +    &
                                     q_dt(i,j,k,ice_wat) )
              delp(i,j,k) = delp(i,j,k) * ps_dt(i,j)
           enddo
        enddo
   elseif ( nwat>0 ) then
        do j=js,je
           do i=is,ie
              ps_dt(i,j)  = 1. + dt*sum(q_dt(i,j,k,1:nwat))
              delp(i,j,k) = delp(i,j,k) * ps_dt(i,j)
           enddo
        enddo
   endif

!-----------------------------------------
! Adjust mass mixing ratio of all tracers 
!-----------------------------------------
   if ( nwat /=0 ) then
      do m=1,ncnst   
      if( m /= cld_amt ) then  ! cloud fraction in GFDL physics
          do j=js,je
             do i=is,ie
                q(i,j,k,m) = q(i,j,k,m) / ps_dt(i,j)
             enddo
          enddo
      endif
      enddo
   endif

1000 continue

! [delp, (ua, va), pt, q] updated. Perform nudging if requested

!------- nudging of atmospheric variables toward specified data --------

#ifdef GFDL_NUDGE
    if (nudge) then
        ps_dt(:,:) = 0.
!--------------------------------------------
! All fields will be updated; tendencies added
!--------------------------------------------
        call get_atmos_nudge ( Time, dt, beglon, endlon, beglat, endlat,    &
             npz, ng, ps(beglon:endlon,:), ua(beglon:endlon,:,:), &
             va(beglon:endlon,:,:), pt(beglon:endlon,:,:), &
             q(beglon:endlon,:,:,:), ps_dt(beglon:endlon,:), u_dt(beglon:endlon,:,:),  & 
             v_dt(beglon:endlon,:,:), t_dt(beglon:endlon,:,:), &
             q_dt(beglon:endlon,:,:,:) )

        if (do_ps) then
!--------------
! Update delp
!--------------
            do k=1,npz
               dbk = dt * (bk(k+1) - bk(k))
               do j=js,je
                  do i=is,ie
                     delp(i,j,k) = delp(i,j,k) + dbk*ps_dt(i,j)
                  enddo
               enddo
            enddo
        endif
    endif
#endif


!----------------------------------------
! Update pe, peln, pkz, and surface winds
!----------------------------------------
  if ( fv_debug ) then
       call prt_maxmin('PS_b_update',     ps, is, ie, js,  je, ng,   1, 0.01, gid==0)
       call prt_maxmin('delp_a_update', delp, is, ie, js,  je, ng, npz, 0.01, gid==0)
  endif

!$omp parallel do private(i, j, k)
   do j=js,je
      do k=2,npz+1                                                                             
         do i=is,ie
            pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
            pk(i,j,k) = pe(i,k,j) ** kappa
            peln(i,k,j) = log(pe(i,k,j))
         enddo
      enddo

      do i=is,ie
            ps(i,j) = pe(i,npz+1,j)
         u_srf(i,j) = ua(i,j,npz)
         v_srf(i,j) = va(i,j,npz)
      enddo

      if ( hydrostatic ) then
         do k=1,npz
            do i=is,ie
               pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k) ) /  &
                            (kappa*(peln(i,k+1,j)-peln(i,k,j)))
            enddo
         enddo
      endif
   enddo      ! j-loop

!-------------------------------------------------------------------------
! Re-compute the full (nonhydrostatic) pressure due to temperature changes
!-------------------------------------------------------------------------
! The assumption here is that "delz" is not changed by the diabatic processes
! and therefore mass & density remain the same within the finite-volume.
! Potential temperature will change becuase (full nonhydrostatic) pressure
! will be changed according to the gas law:
!                                           p = density * R * T
    if ( .not.hydrostatic ) then
      if ( k_top>1 ) then
         do k=1,k_top-1
            do j=js,je
               do i=is,ie
                  pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k)) /     &
                         (kappa*(peln(i,k+1,j)-peln(i,k,j)))
               enddo
            enddo
         enddo
      endif

      do k=k_top,npz
         do j=js,je
            do i=is,ie
               pkz(i,j,k) = ( rdg*delp(i,j,k)*pt(i,j,k)/delz(i,j,k) )**kappa
            enddo
         enddo
      enddo
    endif
                                                    call timing_on(' Update_dwinds')
    call update_dwinds_phys(is, ie, js, je, isd, ied, jsd, jed, dt, u_dt, v_dt, u, v)
                                                    call timing_off(' Update_dwinds')

  if ( fv_debug ) then
       call prt_maxmin('PS_a_update', ps, is, ie, js, je, ng,   1, 0.01, gid==0)
  endif

  end subroutine fv_update_phys


  subroutine update_dwinds_phys(is, ie, js, je, isd, ied, jsd, jed, dt, u_dt, v_dt, u, v)

! Purpose; Transform wind tendencies on A grid to D grid for the final update
 
  integer, intent(in):: is,  ie,  js,  je
  integer, intent(in):: isd, ied, jsd, jed
  real,    intent(in):: dt
  real, intent(inout):: u(isd:ied,  jsd:jed+1,npz)
  real, intent(inout):: v(isd:ied+1,jsd:jed  ,npz)
  real, intent(inout), dimension(isd:ied,jsd:jed,npz):: u_dt, v_dt

! local:
  real v3(is-1:ie+1,js-1:je+1,3)
  real ue(is-1:ie+1,js:je+1,3)    ! 3D winds at edges
  real ve(is:ie+1,js-1:je+1,  3)    ! 3D winds at edges
  real, dimension(is:ie):: ut1, ut2, ut3
  real, dimension(js:je):: vt1, vt2, vt3
  real dt5
  integer i, j, k, im2, jm2


       call timing_on('COMM_TOTAL')
  call mpp_update_domains(u_dt, domain, complete=.false.)
  call mpp_update_domains(v_dt, domain, complete=.true.)
       call timing_off('COMM_TOTAL')

    dt5 = 0.5 * dt
    im2 = (npx-1)/2
    jm2 = (npy-1)/2
!$omp parallel do private (i,j,k)

    do k=1, npz

     if ( grid_type > 3 ) then    ! Local & one tile configurations

       do j=js,je+1
          do i=is,ie
             u(i,j,k) = u(i,j,k) + dt5*(u_dt(i,j-1,k) + u_dt(i,j,k))
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             v(i,j,k) = v(i,j,k) + dt5*(v_dt(i-1,j,k) + v_dt(i,j,k))
          enddo
       enddo

     else
! Compute 3D wind tendency on A grid
       do j=js-1,je+1
          do i=is-1,ie+1
             v3(i,j,1) = u_dt(i,j,k)*vlon(i,j,1) + v_dt(i,j,k)*vlat(i,j,1)
             v3(i,j,2) = u_dt(i,j,k)*vlon(i,j,2) + v_dt(i,j,k)*vlat(i,j,2)
             v3(i,j,3) = u_dt(i,j,k)*vlon(i,j,3) + v_dt(i,j,k)*vlat(i,j,3)
          enddo
       enddo

! A --> D
! Interpolate to cell edges
       do j=js,je+1
          do i=is-1,ie+1
             ue(i,j,1) = v3(i,j-1,1) + v3(i,j,1)
             ue(i,j,2) = v3(i,j-1,2) + v3(i,j,2)
             ue(i,j,3) = v3(i,j-1,3) + v3(i,j,3)
          enddo
       enddo

       do j=js-1,je+1
          do i=is,ie+1
             ve(i,j,1) = v3(i-1,j,1) + v3(i,j,1)
             ve(i,j,2) = v3(i-1,j,2) + v3(i,j,2)
             ve(i,j,3) = v3(i-1,j,3) + v3(i,j,3)
          enddo
       enddo

! --- E_W edges (for v-wind):
     if ( is==1 ) then
       i = 1
       do j=js,je
        if ( j>jm2 ) then
             vt1(j) = edge_vect_w(j)*ve(i,j-1,1)+(1.-edge_vect_w(j))*ve(i,j,1)
             vt2(j) = edge_vect_w(j)*ve(i,j-1,2)+(1.-edge_vect_w(j))*ve(i,j,2)
             vt3(j) = edge_vect_w(j)*ve(i,j-1,3)+(1.-edge_vect_w(j))*ve(i,j,3)
        else
             vt1(j) = edge_vect_w(j)*ve(i,j+1,1)+(1.-edge_vect_w(j))*ve(i,j,1)
             vt2(j) = edge_vect_w(j)*ve(i,j+1,2)+(1.-edge_vect_w(j))*ve(i,j,2)
             vt3(j) = edge_vect_w(j)*ve(i,j+1,3)+(1.-edge_vect_w(j))*ve(i,j,3)
        endif
       enddo
       do j=js,je
          ve(i,j,1) = vt1(j)
          ve(i,j,2) = vt2(j)
          ve(i,j,3) = vt3(j)
       enddo
     endif
     if ( (ie+1)==npx ) then
       i = npx
       do j=js,je
        if ( j>jm2 ) then
             vt1(j) = edge_vect_e(j)*ve(i,j-1,1)+(1.-edge_vect_e(j))*ve(i,j,1)
             vt2(j) = edge_vect_e(j)*ve(i,j-1,2)+(1.-edge_vect_e(j))*ve(i,j,2)
             vt3(j) = edge_vect_e(j)*ve(i,j-1,3)+(1.-edge_vect_e(j))*ve(i,j,3)
        else
             vt1(j) = edge_vect_e(j)*ve(i,j+1,1)+(1.-edge_vect_e(j))*ve(i,j,1)
             vt2(j) = edge_vect_e(j)*ve(i,j+1,2)+(1.-edge_vect_e(j))*ve(i,j,2)
             vt3(j) = edge_vect_e(j)*ve(i,j+1,3)+(1.-edge_vect_e(j))*ve(i,j,3)
        endif
       enddo
       do j=js,je
          ve(i,j,1) = vt1(j)
          ve(i,j,2) = vt2(j)
          ve(i,j,3) = vt3(j)
       enddo
     endif
! N-S edges (for u-wind):
     if ( js==1 ) then
       j = 1
       do i=is,ie
        if ( i>im2 ) then
             ut1(i) = edge_vect_s(i)*ue(i-1,j,1)+(1.-edge_vect_s(i))*ue(i,j,1)
             ut2(i) = edge_vect_s(i)*ue(i-1,j,2)+(1.-edge_vect_s(i))*ue(i,j,2)
             ut3(i) = edge_vect_s(i)*ue(i-1,j,3)+(1.-edge_vect_s(i))*ue(i,j,3)
        else
             ut1(i) = edge_vect_s(i)*ue(i+1,j,1)+(1.-edge_vect_s(i))*ue(i,j,1)
             ut2(i) = edge_vect_s(i)*ue(i+1,j,2)+(1.-edge_vect_s(i))*ue(i,j,2)
             ut3(i) = edge_vect_s(i)*ue(i+1,j,3)+(1.-edge_vect_s(i))*ue(i,j,3)
        endif
       enddo
       do i=is,ie
          ue(i,j,1) = ut1(i)
          ue(i,j,2) = ut2(i)
          ue(i,j,3) = ut3(i)
       enddo
     endif
     if ( (je+1)==npy ) then
       j = npy
       do i=is,ie
        if ( i>im2 ) then
             ut1(i) = edge_vect_n(i)*ue(i-1,j,1)+(1.-edge_vect_n(i))*ue(i,j,1)
             ut2(i) = edge_vect_n(i)*ue(i-1,j,2)+(1.-edge_vect_n(i))*ue(i,j,2)
             ut3(i) = edge_vect_n(i)*ue(i-1,j,3)+(1.-edge_vect_n(i))*ue(i,j,3)
        else
             ut1(i) = edge_vect_n(i)*ue(i+1,j,1)+(1.-edge_vect_n(i))*ue(i,j,1)
             ut2(i) = edge_vect_n(i)*ue(i+1,j,2)+(1.-edge_vect_n(i))*ue(i,j,2)
             ut3(i) = edge_vect_n(i)*ue(i+1,j,3)+(1.-edge_vect_n(i))*ue(i,j,3)
        endif
       enddo
       do i=is,ie
          ue(i,j,1) = ut1(i)
          ue(i,j,2) = ut2(i)
          ue(i,j,3) = ut3(i)
       enddo
     endif

! Update:
       do j=js,je+1
          do i=is,ie
             u(i,j,k) = u(i,j,k) + dt5*( ue(i,j,1)*es(1,i,j,1) +  &
                                         ue(i,j,2)*es(2,i,j,1) +  &
                                         ue(i,j,3)*es(3,i,j,1) )
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             v(i,j,k) = v(i,j,k) + dt5*( ve(i,j,1)*ew(1,i,j,2) +  &
                                         ve(i,j,2)*ew(2,i,j,2) +  &
                                         ve(i,j,3)*ew(3,i,j,2) )
          enddo
       enddo

      endif   ! end grid_type
 
    enddo         ! k-loop

  end subroutine update_dwinds_phys 

end module fv_update_phys_mod
