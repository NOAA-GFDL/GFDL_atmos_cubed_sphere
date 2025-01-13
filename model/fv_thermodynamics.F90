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
module fv_thermodynamics_mod

#ifdef OVERLOAD_R4
  use constantsR4_mod,     only: grav, cp_air, cp_vapor, rvgas, rdgas
#else
  use constants_mod,       only: grav, cp_air, cp_vapor, rvgas, rdgas
#endif
  use gfdl_mp_mod,         only: c_liq, c_ice
  use field_manager_mod,   only: MODEL_ATMOS
  use tracer_manager_mod,  only: get_tracer_name
  use fv_arrays_mod,       only: fv_grid_bounds_type, fv_thermo_type, fv_flags_type
  use mpp_mod,             only: mpp_error, FATAL, input_nml_file
  use fms_mod,             only: check_nml_error

  implicit none
  real, parameter:: cv_vap = 3.*rvgas  ! 1384.5
  real, parameter:: cv_air =  cp_air - rdgas ! = rdgas * (7/2-1) = 2.5*rdgas=717.68
  real, parameter:: tice = 273.16

  public fv_thermo_init, compute_total_energy, moist_cv, moist_cp, compute_q_con

contains


  subroutine fv_thermo_init(thermostruct,flagstruct)

    type(fv_thermo_type), intent(INOUT), target :: thermostruct
    type(fv_flags_type), intent(INOUT) :: flagstruct

    integer :: f_unit, ios, ierr, dum
    logical, pointer :: use_cond, moist_kappa
    namelist /fv_thermo_nml/ use_cond, moist_kappa

    use_cond => thermostruct%use_cond
    moist_kappa => thermostruct%moist_kappa

    !For hydrostatic dynamics, set default to .false. for both
    ! to maintain compatibility with the usual hydrostatic configuration
    if (flagstruct%hydrostatic) then
       use_cond = .false.
       moist_kappa = .false.
    endif


    !read namelist
    read (input_nml_file,fv_thermo_nml,iostat=ios)
    ierr = check_nml_error(ios,'fv_thermo_nml')

    if (moist_kappa .and. .not. use_cond) then
       call mpp_error(FATAL, " moist_kappa = .true. and use_cond = .false. not supported.")
    endif

    if (flagstruct%hydrostatic .and. moist_kappa) then
       call mpp_error(FATAL, " moist_kappa not yet supported for hydrostatic simulation.")
    endif

    if (flagstruct%hydrostatic .and. use_cond) then
       call mpp_error(FATAL, " use_cond not yet supported for hydrostatic simulation.")
    endif

    thermostruct%is_initialized = .true.

  end subroutine fv_thermo_init



 subroutine compute_total_energy(is, ie, js, je, isd, ied, jsd, jed, km,       &
                                 u, v, w, delz, pt, delp, q, qc, q_con, pe, peln, hs, &
                                 rsin2_l, cosa_s_l, &
                                 r_vir,  cp, rg, hlv, te_2d, ua, va, teq, &
                                 moist_phys, nwat, sphum, liq_wat, rainwat, &
                                 ice_wat, snowwat, graupel, hydrostatic, &
                                 moist_kappa, id_te)
!------------------------------------------------------
! Compute vertically integrated total energy per column
!------------------------------------------------------
! !INPUT PARAMETERS:
   integer,  intent(in):: km, is, ie, js, je, isd, ied, jsd, jed, id_te
   integer,  intent(in):: sphum, liq_wat, ice_wat, rainwat, snowwat, graupel, nwat
   real, intent(inout), dimension(isd:ied,jsd:jed,km):: ua, va
   real, intent(in), dimension(isd:ied,jsd:jed,km):: pt, delp
   real, intent(in), dimension(isd:ied,jsd:jed,km,*):: q
   real, intent(in), dimension(isd:ied,jsd:jed,km):: qc, q_con !virtual adjustment zvir*qv
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
   logical, intent(in):: moist_phys, hydrostatic, moist_kappa
! Output:
   real, intent(out):: te_2d(is:ie,js:je)   ! vertically integrated TE
   real, intent(out)::   teq(is:ie,js:je)   ! Moist TE
! Local
   real, dimension(is:ie,km):: tv
   real  phiz(is:ie,km+1)
   real  cvm(is:ie), qd(is:ie)
   integer i, j, k

!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,km,hydrostatic,hs,pt,qc,q_con,rg,peln,te_2d, &
!$OMP                                  pe,delp,cp,rsin2_l,u,v,cosa_s_l,delz,moist_phys,w, &
!$OMP                                  q,nwat,liq_wat,rainwat,ice_wat,snowwat,graupel,sphum,moist_kappa)  &
!$OMP                          private(phiz, tv, cvm, qd)
  do j=js,je

     if ( hydrostatic ) then

        do i=is,ie
           phiz(i,km+1) = hs(i,j)
        enddo
        do k=km,1,-1
           do i=is,ie
#ifdef USE_COND
                tv(i,k) = pt(i,j,k)*(1.+qc(i,j,k))*(1-q_con(i,j,k))
#else
                tv(i,k) = pt(i,j,k)*(1.+qc(i,j,k))
#endif
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
        if (moist_kappa) then
           call moist_cv(is,ie,isd,ied,jsd,jed, km, j, k, nwat, sphum, liq_wat, rainwat,    &
                ice_wat, snowwat, graupel, q, qd, cvm)
           do i=is,ie
              te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*( cvm(i)*pt(i,j,k) +  &
                        0.5*(phiz(i,k)+phiz(i,k+1)+w(i,j,k)**2+0.5*rsin2_l(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                        v(i,j,k)**2+v(i+1,j,k)**2-(u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s_l(i,j))))
           enddo
        else
           do i=is,ie
              te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*( cv_air*pt(i,j,k) +  &
                        0.5*(phiz(i,k)+phiz(i,k+1)+w(i,j,k)**2+0.5*rsin2_l(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                        v(i,j,k)**2+v(i+1,j,k)**2-(u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s_l(i,j))))
           enddo
        endif
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

#ifdef THERMO_PROTOTYPE
 subroutine fv_thermodynamics_init

   !set up heat capacities for each tracer

   do n=1,min(ncnst,nwat)
      tracer_name = ...
      if ( 'sphum' == trim(tracer_name)) then
         dcv(n) = cv_vap - cv_air
         dcp(n) = cp_vap - cp_air
      else if ( ANY( (/'liq_wat','rainwat'/) == trim(tracer_name)) ) then
         dcv(n) = c_liq - cv_air
         dcp(n) = c_liq - cp_air
      else if ( ANY( (/'ice_wat', 'snowwat', 'graupel', 'hailwat'/) == trim(tracer_name)) ) then
         dcv(n) = c_ice - cv_air
         dcp(n) = c_ice - cp_air
      endif
   enddo

 end subroutine fv_thermodynamics_init
#endif


 subroutine moist_cv(is,ie, isd,ied, jsd,jed, km, j, k, nwat, sphum, liq_wat, rainwat,    &
                     ice_wat, snowwat, graupel, q, q_con, cvm, t1)
  integer, intent(in):: is, ie, isd,ied, jsd,jed, km, nwat, j, k
  integer, intent(in):: sphum, liq_wat, rainwat, ice_wat, snowwat, graupel
  real, intent(in), dimension(isd:ied,jsd:jed,km,nwat):: q
  real, intent(out), dimension(is:ie):: cvm, q_con
  real, intent(in), optional:: t1(is:ie)
!
  real, parameter:: t_i0 = 15.
  real, dimension(is:ie):: qv, ql, qs
  integer:: i

  select case (nwat)

   case(2)
     if ( present(t1) ) then  ! Special case for GFS physics
        do i=is,ie
           q_con(i) = max(0., q(i,j,k,liq_wat))
           if ( t1(i) > tice ) then
                qs(i) = 0.
           elseif ( t1(i) < tice-t_i0 ) then
                qs(i) = q_con(i)
           else
                qs(i) = q_con(i)*(tice-t1(i))/t_i0
           endif
           ql(i) = q_con(i) - qs(i)
           qv(i) = max(0.,q(i,j,k,sphum))
           cvm(i) = (1.-(qv(i)+q_con(i)))*cv_air + qv(i)*cv_vap + ql(i)*c_liq + qs(i)*c_ice
        enddo
     else
        do i=is,ie
           qv(i) = max(0.,q(i,j,k,sphum))
           qs(i) = max(0.,q(i,j,k,liq_wat))
           q_con(i) = qs(i)
           cvm(i) = (1.-qv(i))*cv_air + qv(i)*cv_vap
        enddo
     endif
  case (3)
     do i=is,ie
        qv(i) = q(i,j,k,sphum)
        ql(i) = q(i,j,k,liq_wat)
        qs(i) = q(i,j,k,ice_wat)
        q_con(i) = ql(i) + qs(i)
        cvm(i) = (1.-(qv(i)+q_con(i)))*cv_air + qv(i)*cv_vap + ql(i)*c_liq + qs(i)*c_ice
     enddo
  case(4)              ! K_warm_rain with fake ice
     do i=is,ie
        qv(i) = q(i,j,k,sphum)
        q_con(i) = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
        cvm(i) = (1.-(qv(i)+q_con(i)))*cv_air + qv(i)*cv_vap + q_con(i)*c_liq
     enddo
  case(5)
     do i=is,ie
        qv(i) = q(i,j,k,sphum)
        ql(i) = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
        qs(i) = q(i,j,k,ice_wat) + q(i,j,k,snowwat)
        q_con(i) = ql(i) + qs(i)
        cvm(i) = (1.-(qv(i)+q_con(i)))*cv_air + qv(i)*cv_vap + ql(i)*c_liq + qs(i)*c_ice
     enddo
  case(6)
     do i=is,ie
        qv(i) = q(i,j,k,sphum)
        ql(i) = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
        qs(i) = q(i,j,k,ice_wat) + q(i,j,k,snowwat) + q(i,j,k,graupel)
        q_con(i) = ql(i) + qs(i)
        cvm(i) = (1.-(qv(i)+q_con(i)))*cv_air + qv(i)*cv_vap + ql(i)*c_liq + qs(i)*c_ice
     enddo
  case default
     !call mpp_error (NOTE, 'fv_mapz::moist_cv - using default cv_air')
     do i=is,ie
         q_con(i) = 0.
        cvm(i) = cv_air
     enddo
 end select

 end subroutine moist_cv

 subroutine moist_cp(is,ie, isd,ied, jsd,jed, km, j, k, nwat, sphum, liq_wat, rainwat,    &
                     ice_wat, snowwat, graupel, q, q_con, cpm, t1)

  integer, intent(in):: is, ie, isd,ied, jsd,jed, km, nwat, j, k
  integer, intent(in):: sphum, liq_wat, rainwat, ice_wat, snowwat, graupel
  real, intent(in), dimension(isd:ied,jsd:jed,km,nwat):: q
  real, intent(out), dimension(is:ie):: cpm, q_con
  real, intent(in), optional:: t1(is:ie)
!
  real, parameter:: t_i0 = 15.
  real, dimension(is:ie):: qv, ql, qs
  integer:: i

  select case (nwat)

  case(2)
     if ( present(t1) ) then  ! Special case for GFS physics
        do i=is,ie
           q_con(i) = max(0., q(i,j,k,liq_wat))
           if ( t1(i) > tice ) then
                qs(i) = 0.
           elseif ( t1(i) < tice-t_i0 ) then
                qs(i) = q_con(i)
           else
                qs(i) = q_con(i)*(tice-t1(i))/t_i0
           endif
           ql(i) = q_con(i) - qs(i)
           qv(i) = max(0.,q(i,j,k,sphum))
           cpm(i) = (1.-(qv(i)+q_con(i)))*cp_air + qv(i)*cp_vapor + ql(i)*c_liq + qs(i)*c_ice
        enddo
     else
     do i=is,ie
        qv(i) = max(0.,q(i,j,k,sphum))
        qs(i) = max(0.,q(i,j,k,liq_wat))
        q_con(i) = qs(i)
        cpm(i) = (1.-qv(i))*cp_air + qv(i)*cp_vapor
     enddo
     endif

  case(3)
     do i=is,ie
        qv(i) = q(i,j,k,sphum)
        ql(i) = q(i,j,k,liq_wat)
        qs(i) = q(i,j,k,ice_wat)
        q_con(i) = ql(i) + qs(i)
        cpm(i) = (1.-(qv(i)+q_con(i)))*cp_air + qv(i)*cp_vapor + ql(i)*c_liq + qs(i)*c_ice
     enddo
  case(4)    ! K_warm_rain scheme with fake ice
     do i=is,ie
        qv(i) = q(i,j,k,sphum)
        q_con(i) = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
        cpm(i) = (1.-(qv(i)+q_con(i)))*cp_air + qv(i)*cp_vapor + q_con(i)*c_liq
     enddo
  case(5)
     do i=is,ie
        qv(i) = q(i,j,k,sphum)
        ql(i) = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
        qs(i) = q(i,j,k,ice_wat) + q(i,j,k,snowwat)
        q_con(i) = ql(i) + qs(i)
        cpm(i) = (1.-(qv(i)+q_con(i)))*cp_air + qv(i)*cp_vapor + ql(i)*c_liq + qs(i)*c_ice
     enddo
  case(6)
     do i=is,ie
        qv(i) = q(i,j,k,sphum)
        ql(i) = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
        qs(i) = q(i,j,k,ice_wat) + q(i,j,k,snowwat) + q(i,j,k,graupel)
        q_con(i) = ql(i) + qs(i)
        cpm(i) = (1.-(qv(i)+q_con(i)))*cp_air + qv(i)*cp_vapor + ql(i)*c_liq + qs(i)*c_ice
     enddo
  case default
     !call mpp_error (NOTE, 'fv_mapz::moist_cp - using default cp_air')
     do i=is,ie
        q_con(i) = 0.
        cpm(i) = cp_air
     enddo
  end select

 end subroutine moist_cp

 subroutine compute_q_con(bd, npz, nwat, nq, q, q_con)
   type(fv_grid_bounds_type), intent(IN) :: bd
   integer, intent(in):: npz, nwat, nq
   real, intent(in), dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq):: q
   real, intent(out), dimension(bd%is:bd%ie,bd%js:bd%je,npz):: q_con

   integer:: n, dum
   character(len=32) :: tracer_name
   integer :: is,  ie,  js,  je
   integer :: isd, ied, jsd, jed

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed

   !Not optimized for OpenMP yet

   q_con = 0.
   do n=1,nwat
      dum = get_tracer_name(MODEL_ATMOS, n, tracer_name)
      select case (trim(tracer_name))
         case('liq_wat','rainwat')
            q_con = q_con + q(is:ie,js:je,:,n)
         case('ice_wat','snowwat','graupel','hailwat')
            q_con = q_con + q(is:ie,js:je,:,n)
      end select
   enddo

 end subroutine compute_q_con

! subroutine compute_moist_kappa(!!
!
! end subroutine compute_moist_kappa


end module fv_thermodynamics_mod
