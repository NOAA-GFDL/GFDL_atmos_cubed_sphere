module fv_cmp_mod

  use constants_mod, only: pi=>pi_8, rvgas, rdgas, grav, hlv, hlf, cp_air, cp_vapor, R_GRID
  use fv_mp_mod,     only: is_master

  implicit none
  real, parameter:: cv_vap = 3.*rvgas  ! 1384.8
  real, parameter:: cv_air =  cp_air - rdgas ! = rdgas * (7/2-1) = 2.5*rdgas=717.68
  real, parameter:: c_ice = 2106.        ! heat capacity of ice at 0.C
  real, parameter:: c_liq = 4.1855e+3    ! GFS
  real, parameter:: cp_vap = cp_vapor    ! 4*rv_gas=1846.
  real, parameter:: dc_vap = cp_vap - c_liq     ! = -2344.    isobaric heating/cooling
  real, parameter:: dc_ice =  c_liq - c_ice     ! =  2084
  real, parameter:: tice = 273.16
  real, parameter:: t_wfr = tice - 40.
! Values at 0 Deg C
 real, parameter:: hlv0 = 2.5e6
 real, parameter:: hlf0 = 3.3358e5
! Latent heat at absolute zero:
 real, parameter:: Lv0  = hlv0 - dc_vap*tice   ! = 3.141264e6
 real, parameter:: li00 = hlf0 - dc_ice*tice   ! = -2.355446e5
 real(kind=R_GRID), parameter:: e00 = 610.71  ! saturation vapor pressure at T0
 real(kind=R_GRID), parameter:: d2ice  = cp_vap - c_ice
 real(kind=R_GRID), parameter:: Li2 = hlv0+hlf0 - d2ice*tice
! Local:
 real:: ql_gen = 1.0e-3    ! max ql generation during remapping step if fast_sat_adj = .T.
 real:: qi_gen = 1.82E-6
 real:: qi_lim = 2.  ! 2.
 real:: tau_i2s = 1000.
 real:: tau_v2l = 150.
 real:: tau_l2v = 300.
 real:: tau_r  = 600.       ! rain freezing time scale during fast_sat
 real:: tau_s  = 600.       ! snow melt
 real:: tau_mlt = 600.      ! ice melting time-scale
 real, parameter:: tau_l2r = 300.
 real:: sat_adj0 = 0.9  !  0.95
 real:: qi0_max = 1.2e-4    ! Max: ice  --> snow autocon threshold
 real:: ql0_max = 2.0e-3    ! max ql value (auto converted to rain)
 real:: t_sub   = 184.  ! Min temp for sublimation of cloud ice
 real:: cld_min = 0.05
 real:: dw_ocean = 0.12 ! 0.1
 real:: cracw = 3.272
 real:: crevp(5), lat2
 real, allocatable:: table(:), table2(:), tablew(:), des2(:), desw(:)

 logical:: rad_rain = .true.
 logical:: rad_snow = .true.
 logical:: rad_graupel = .true.
 logical:: mp_initialized = .false.

 private
 public fv_sat_adj, qs_init

contains

 subroutine fv_sat_adj(mdt, zvir, is, ie, js, je, ng, hydrostatic, consv_te, &
                       te0, qv, ql, qi, qr, qs, qg, dpln, delz, pt, dp,  &
                       q_con, cappa, area, dtdt, out_dt, last_step, do_qa, qa)
! This is designed for 6-class micro-physics schemes; handles the heat release
! due to in situ phase changes
! input pt is T_vir
 integer, intent(in):: is, ie, js, je, ng
 real, intent(in):: mdt ! remapping time step
 real, intent(in):: zvir
 logical, intent(in):: hydrostatic, consv_te, out_dt
 logical, intent(in):: last_step
 logical, intent(in):: do_qa
 real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng):: dp, delz
 real, intent(in):: dpln(is:ie,js:je)
 real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng):: pt, qv, ql, qi, qr, qs, qg
 real, intent(out):: qa(is-ng:ie+ng,js-ng:je+ng)
 real(kind=R_GRID), intent(in), dimension(is-ng:ie+ng,js-ng:je+ng):: area
 real, intent(inout), dimension(is-ng:,js-ng:):: q_con
 real, intent(inout), dimension(is-ng:,js-ng:):: cappa
 real, intent(inout)::dtdt(is:ie,js:je)
 real, intent(out):: te0(is-ng:ie+ng,js-ng:je+ng)
!---
 real, dimension(is:ie):: wqsat, dq2dt, qpz, cpm, t0, pt1, icp2, lcp2, tcp2, tcp3,    &
                          den, q_liq, q_sol, src, p1, hvar
 real:: sink, qsw, rh, fac_v2l, fac_l2v
 real:: tc, qsi, dqsdt, dq, dq0, pidep, qi_crt, tmp, dtmp
 real:: condensates, tin, qstar, rqi, q_plus, q_minus
 real:: sdt, dt_Bigg, adj_fac, fac_s, fac_r, fac_i2s, fac_mlt, fac_l2r
 real:: factor, qim, tice0
 integer i,j

 sdt = 0.5 * mdt
 dt_Bigg = 3.3333e-10*mdt
 tice0 = tice - 0.01

 fac_i2s = 1. - exp( -mdt/tau_i2s )        !
 fac_v2l = 1. - exp( -sdt/tau_v2l )
 fac_l2v = 1. - exp( -sdt/tau_l2v )        !
 fac_l2v = min(sat_adj0, fac_l2v)
 fac_mlt = 1. - exp( -sdt/tau_mlt )
   fac_r = 1. - exp( -mdt/tau_r )
   fac_s = 1. - exp( -mdt/tau_s )
 fac_l2r = 1. - exp( -mdt/tau_l2r )

!!! hvar(:) = 0.02

 do j=js, je

    do i=is, ie
       q_liq(i) = ql(i,j) + qr(i,j)
       q_sol(i) = qi(i,j) + qs(i,j) + qg(i,j)
         qpz(i) = q_liq(i) + q_sol(i)
#ifdef USE_COND
       q_con(i,j) = qpz(i)
       pt1(i) = pt(i,j) / ((1.+zvir*qv(i,j))*(1-qpz(i)))
#else
       pt1(i) = pt(i,j) / (1.+zvir*qv(i,j))
#endif
        t0(i) = pt1(i)     ! true temperature
        hvar(i) = min(0.2, max(0.01, dw_ocean*sqrt(sqrt(area(i,j)/1.E10))) )
        qpz(i) = qpz(i) + qv(i,j)    ! conserved in this routine
    enddo

    if ( hydrostatic ) then
        do i=is, ie
           den(i) = dp(i,j)/(dpln(i,j)*rdgas*pt(i,j))
           cpm(i) = (1.-qpz(i))*cp_air + qv(i,j)*cp_vap + q_liq(i)*c_liq + q_sol(i)*c_ice
        enddo
    else
        do i=is, ie
           den(i) = -dp(i,j)/(grav*delz(i,j))   ! Moist_air density
           cpm(i) = (1.-qpz(i))*cv_air + qv(i,j)*cv_vap + q_liq(i)*c_liq + q_sol(i)*c_ice
        enddo
    endif

    do i=is, ie
       lcp2(i) = (Lv0 +dc_vap*t0(i)) / cpm(i)
       icp2(i) = (li00+dc_ice*t0(i)) / cpm(i)
       tcp2(i) = lcp2(i) + icp2(i)
! Compute special heat capacity for qv --> ql (dqsdt term)
       tcp3(i) = lcp2(i) + icp2(i)*min(1., dim(tice,t0(i))/40.)
       src(i) = 0.
    enddo

    if ( consv_te ) then 
         do i=is, ie
! Convert to "liquid" form
            te0(i,j) = 0.
!!!!!!!!!!!            te0(i,j) = dp(i,j)*cpm(i)*(icp2(i)*q_sol(i) - lcp2(i)*qv(i,j))
         enddo
    endif

    do i=is, ie
       if( qi(i,j) < 0. ) then
! Fix negative cloud ice with snow
           qs(i,j) = qs(i,j) + qi(i,j)
           qi(i,j) = 0.
       elseif ( qi(i,j)>1.E-8 .and. pt1(i) > tice ) then
! Melting of cloud ice into cloud water ********
           sink = min( qi(i,j), fac_mlt*(pt1(i)-tice)/icp2(i) )
! Maximum amount of melted ice converted to ql
            tmp = min( sink, dim(ql0_max, ql(i,j)) )   ! max ql amount
           ql(i,j) = ql(i,j) + tmp
           qr(i,j) = qr(i,j) + sink - tmp
           qi(i,j) = qi(i,j) - sink
            pt1(i) =  pt1(i) - sink*icp2(i)
       endif
    enddo

    do i=is, ie
       if( qs(i,j) < 0. ) then
! Fix negative snow with graupel
           qg(i,j) = qg(i,j) + qs(i,j)
           qs(i,j) = 0.
       elseif( qg(i,j) < 0. ) then
! Fix negative graupel with available snow
           tmp = min( -qg(i,j), max(0., qs(i,j)) )
           qg(i,j) = qg(i,j) + tmp
           qs(i,j) = qs(i,j) - tmp
       endif
    enddo
! After this point cloud ice & snow are positive definite

    do i=is, ie
       if( ql(i,j) < 0. ) then
! Fix negative cloud water if rain exists
           tmp = min( -ql(i,j), max(0., qr(i,j)) )
           ql(i,j) = ql(i,j) + tmp
           qr(i,j) = qr(i,j) - tmp
       elseif( qr(i,j) < 0. ) then
! Fix negative rain water with available cloud water
           tmp = min( -qr(i,j), max(0., ql(i,j)) )
           ql(i,j) = ql(i,j) - tmp
           qr(i,j) = qr(i,j) + tmp
       endif
    enddo

!-----------------------------------------
! Enforce complete freezing below -48 C
!-----------------------------------------
! (do this before ql & qr to evap)
    do i=is, ie
       dtmp = tice - 48. - pt1(i)
! cloud water --> ice
       if( ql(i,j)>0. .and. dtmp > 0. ) then
           sink = min( ql(i,j),  dtmp/icp2(i) )
           ql(i,j) = ql(i,j) - sink
           qi(i,j) = qi(i,j) + sink
            pt1(i) =  pt1(i) + sink*icp2(i)
       endif
    enddo

 call wqs2_vect(is, ie, pt1, den, wqsat, dq2dt)
 if ( last_step ) then
! Enforce upper (no super_sat) & lower (critical RH) bounds
      do i=is, ie
         dq0 = qv(i,j) - wqsat(i)
         if ( dq0 > 0. ) then ! remove super-saturation
! Prevent super saturation over water:
             src(i) = dq0/(1.+tcp3(i)*dq2dt(i))
         else
! Evaporation of ql
            src(i) = -min( ql(i,j), -dq0/(1.+tcp3(i)*dq2dt(i)) )
         endif
      enddo
      adj_fac = 1.
 else
      adj_fac = sat_adj0
      do i=is, ie
         dq0 = qv(i,j) - wqsat(i)
         if ( dq0 > 0. ) then ! whole grid-box saturated
               tmp = dq0/(1.+tcp3(i)*dq2dt(i))
            src(i) = min(adj_fac*tmp, max(ql_gen-ql(i,j), fac_v2l*tmp))
         else
!                                                     Evaporation of ql
            src(i) = -min( ql(i,j), -fac_l2v*dq0/(1.+tcp3(i)*dq2dt(i)) )
         endif
      enddo
 endif
! Cooling: qi -> (qv, ql, qr); ql -> qv; qr -> qv
! Warming: ql -> (qi, qs); qv -> (ql, qi); qr -> qg

    do i=is, ie
       qv(i,j) = qv(i,j) - src(i)
       ql(i,j) = ql(i,j) + src(i)
        pt1(i) =  pt1(i) + src(i)*lcp2(i)
    enddo

! *********** freezing of cloud water ********
! Enforce complete freezing below -48 C
    do i=is, ie
       dtmp = t_wfr - pt1(i)   ! [-40,-48]
       if( ql(i,j)>0. .and. dtmp > 0. ) then
           sink = min( ql(i,j),  ql(i,j)*dtmp*0.125, dtmp/icp2(i) )
           ql(i,j) = ql(i,j) - sink
           qi(i,j) = qi(i,j) + sink
            pt1(i) =  pt1(i) + sink*icp2(i)
       endif
    enddo

! Bigg mechanism (done only here)
    do i=is, ie
       tc = tice0 - pt1(i)
       if( ql(i,j)>1.E-8 .and. tc > 0. ) then
           sink = dt_Bigg*(exp(0.66*tc)-1.)*den(i)*ql(i,j)**2
           sink = min(ql(i,j), tc/icp2(i), sink)
           ql(i,j) = ql(i,j) - sink
           qi(i,j) = qi(i,j) + sink
            pt1(i) =  pt1(i) + sink*icp2(i)
       endif
    enddo
! *********** freezing of rain water qr-->qg ********
    do i=is, ie
       dtmp = (tice - 1.) - pt1(i)
       if( qr(i,j)>1.E-7 .and. dtmp > 0. ) then
! No limit on freezing below -40 C
            tmp = min( 1., (dtmp*0.025)**2 ) * qr(i,j)
           sink = min( tmp, fac_r*dtmp/icp2(i) )
           qr(i,j) = qr(i,j) - sink
           qg(i,j) = qg(i,j) + sink
            pt1(i) =  pt1(i) + sink*icp2(i)
       endif
    enddo
! *********** Melting of snow qs-->qr ********
    do i=is, ie
       dtmp = pt1(i) - (tice+0.1)
       if( qs(i,j)>1.E-7 .and. dtmp > 0. ) then
! No limter on melting above 10 C
            tmp = min( 1., (dtmp*0.1)**2 ) * qs(i,j)
           sink = min( tmp,  fac_s*dtmp/icp2(i) )
           qs(i,j) = qs(i,j) - sink
           qr(i,j) = qr(i,j) + sink
            pt1(i) =  pt1(i) - sink*icp2(i)
       endif
    enddo

! Update after freezing & before ice-phase adjustment
    do i=is, ie
       q_liq(i) = max(0., ql(i,j) + qr(i,j))
       q_sol(i) = max(0., qi(i,j) + qs(i,j) + qg(i,j))
    enddo

! Enforce upper bounds on ql (can be regarded as autoconversion)
    do i=is, ie
       if ( ql(i,j) > ql0_max ) then
            sink = fac_l2r*(ql(i,j)-ql0_max)
            qr(i,j) = qr(i,j) + sink
            ql(i,j) = ql(i,j) - sink
       endif
    enddo

    if ( hydrostatic ) then
        do i=is, ie
           cpm(i) = (1.-qpz(i))*cp_air + qv(i,j)*cp_vap + q_liq(i)*c_liq + q_sol(i)*c_ice
        enddo
    else
        do i=is, ie
           cpm(i) = (1.-qpz(i))*cv_air + qv(i,j)*cv_vap + q_liq(i)*c_liq + q_sol(i)*c_ice
        enddo
    endif

    do i=is, ie
       lcp2(i) = (Lv0 +dc_vap*pt1(i)) / cpm(i)
       icp2(i) = (li00+dc_ice*pt1(i)) / cpm(i)
       tcp2(i) = lcp2(i) + icp2(i)
       src(i) = 0.
    enddo
 
! Ice-phase
!------------------------------------------
! * pidep: sublimation/deposition of ice:
!------------------------------------------
    do i=is, ie
       p1(i) = dp(i,j)/dpln(i,j)
       if ( pt1(i) < t_sub ) then  ! Too cold to be accurate; freeze qv as a fix
            src(i) = dim(qv(i,j), 1.e-7 )
       elseif ( pt1(i) < tice0 ) then
          qsi = iqs2(pt1(i), den(i), dqsdt)
           dq = qv(i,j) - qsi
         sink = adj_fac*dq/(1.+tcp2(i)*dqsdt)
         if ( qi(i,j) > 1.E-7 ) then
! Eq 9, Hong et al. 2004, MWR; For A and B, see Dudhia 1989: page 3103 Eq (B7)-(B8)
              pidep = sdt*dq*349138.78*exp(0.875*log(qi(i,j)*den(i)))  &
               / (qsi*den(i)*lat2/(0.0243*rvgas*pt1(i)**2) + 4.42478e4)
         else
              pidep = 0.
         endif
         if ( dq > 0. ) then   ! vapor -> ice
             tmp = tice - pt1(i)
             qi_crt = qi_gen*min(qi_lim, 0.1*tmp) / den(i)
             src(i) = min(sink, max(qi_crt-qi(i,j),pidep), tmp/tcp2(i))
         else
             pidep = pidep * min(1., dim(pt1(i),t_sub)*0.2)
             src(i) = max(pidep, sink, -qi(i,j))
         endif
       endif
    enddo
    do i=is, ie
       qv(i,j) = qv(i,j) - src(i)
       qi(i,j) = qi(i,j) + src(i)
        pt1(i) =  pt1(i) + src(i)*tcp2(i)
    enddo

    do i=is, ie
       if( qg(i,j) < 0. ) then
! Fix negative graupel with available ice
           tmp = min( -qg(i,j), max(0., qi(i,j)) )
           qg(i,j) = qg(i,j) + tmp
           qi(i,j) = qi(i,j) - tmp
       endif
    enddo

! Enforce *upper* bounds on qi: consider as autoconversion to qs
    do i=is, ie
       qim = qi0_max / den(i)
       if ( qi(i,j) > qim ) then
            sink = fac_i2s*(qi(i,j)-qim)
            qi(i,j) = qi(i,j) - sink
            qs(i,j) = qs(i,j) + sink
       endif
    enddo
! At this point cloud ice & snow are positive definite
!!! call revap_rac1(hydrostatic, is, ie, sdt, pt1, qv(is,j), ql(is,j), qr(is,j), qi(is,j), qs(is,j), qg(is,j), den, hvar)

! Virtual temp updated !!!
    do i=is, ie
#ifdef USE_COND
       q_liq(i) = ql(i,j) + qr(i,j)
       q_sol(i) = qi(i,j) + qs(i,j) + qg(i,j)
       q_con(i,j) = q_liq(i) + q_sol(i)
       pt(i,j) = pt1(i)*(1.+zvir*qv(i,j))*(1.-q_con(i,j))
       cpm(i) = (1.-qpz(i))*cv_air + qv(i,j)*cv_vap + q_liq(i)*c_liq + q_sol(i)*c_ice
       cappa(i,j) = rdgas/(rdgas + cpm(i)/(1.+zvir*qv(i,j)))
#else
       pt(i,j) = pt1(i)*(1.+zvir*qv(i,j))
#endif
    enddo

    if ( out_dt ) then
         do i=is, ie
            dtdt(i,j) = dtdt(i,j) + pt1(i) - t0(i)
         enddo
    endif

    if ( consv_te ) then 
       if ( hydrostatic ) then
         do i=is, ie
            te0(i,j) = te0(i,j) + cp_air*dp(i,j)*(pt1(i)-t0(i))*(1.+zvir*qv(i,j))
         enddo
       else
         do i=is, ie
#ifdef USE_COND
            te0(i,j) = te0(i,j) + cpm(i)*dp(i,j)*(pt1(i)-t0(i))
#else
            te0(i,j) = te0(i,j) + cv_air*dp(i,j)*(pt1(i)-t0(i))
#endif
         enddo
       endif
!!!      do i=is, ie
!!!         te0(i,j) = te0(i,j) + dp(i,j)*cpm(i)*(lcp2(i)*qv(i,j)-icp2(i)*q_sol(i))
!!!      enddo
    endif

if ( do_qa .and. last_step ) then

    do i=is, ie
       qa(i,j) = 0.
    enddo

    if ( rad_snow ) then
      if ( rad_graupel ) then
       do i=is, ie
          q_sol(i) = qi(i,j) + qs(i,j) + qg(i,j)
       enddo
      else
       do i=is, ie
          q_sol(i) = qi(i,j) + qs(i,j)
       enddo
      endif
    else
       do i=is, ie
          q_sol(i) = qi(i,j)
       enddo
    endif

    do i=is, ie
       if ( rad_rain ) then
            condensates = q_sol(i) + ql(i,j) + qr(i,j)
       else
            condensates = q_sol(i) + ql(i,j)
       endif
! Using the "liquid-frozen water temperature": tin
       tin = pt1(i) - ( lcp2(i)*condensates + icp2(i)*q_sol(i) )  ! minimum  temperature
       if( tin <= t_wfr ) then
           qstar = iqs1(tin, den(i))
       elseif ( tin >= tice ) then
           qstar = wqs1(tin, den(i))
       else
! mixed phase:
           qsi = iqs1(tin, den(i))
           qsw = wqs1(tin, den(i))
           if( condensates > 1.E-6 ) then
               rqi = q_sol(i) / condensates
           else
! Mostly liquid water clouds at initial cloud development stage
               rqi = ( (tice-tin)/(tice-t_wfr) )
           endif
           qstar = rqi*qsi + (1.-rqi)*qsw
       endif

! Partial cloudiness by PDF:
! Assuming subgrid linear distribution in horizontal; this is effectively a smoother for the
! binary cloud scheme;  qa=0.5 if qstar==qpz
 
      rh = qpz(i) / qstar
      if ( rh > 0.75 .and. qpz(i) > 1.E-7 ) then
                dq = hvar(i)*qpz(i)
           q_plus  = qpz(i) + dq
           q_minus = qpz(i) - dq
           if ( qstar < q_minus ) then
                qa(i,j) = 1.
           else
             if ( qstar<q_plus ) then
                qa(i,j) = (q_plus-qstar)/(dq+dq)        ! partial cloud cover:
!               qa(i,j) = sqrt( (q_plus-qstar)/(dq+dq) )
                                                        ! qa = 0 if qstar = q_plus 
                                                        ! qa = 1 if qstar = q_minus
             endif
! Impose minimum cloudiness if substantial condensates exist
             if ( condensates > 1.E-6 ) qa(i,j) = max(cld_min, qa(i,j))
           endif
      endif
   enddo
endif

 enddo

 end subroutine fv_sat_adj


 subroutine revap_rac1(hydrostatic, is, ie, dt, tz, qv, ql, qr, qi, qs, qg, den, hvar)
 logical, intent(in):: hydrostatic
 integer, intent(in):: is, ie
 real, intent(in):: dt       ! time step (s)
 real, intent(in),    dimension(is:ie):: den, hvar, qi, qs, qg
 real, intent(inout), dimension(is:ie):: tz, qv, qr, ql
! local:
 real, parameter:: rh_rain = 0.6
 real, parameter::  sfcrho = 1.2         ! surface air density
 real, dimension(is:ie):: lcp2
 real:: dqv, qsat, dqsdt, evap, qden, q_plus, q_minus, sink
 real:: tin, t2, qpz, dq, dqh, dqm, q_liq, q_sol
 integer i

  if ( hydrostatic ) then
     do i=is, ie
        q_liq = ql(i) + qr(i)
        q_sol = qi(i) + qs(i) + qg(i)
        lcp2(i) = (Lv0+dc_vap*tz(i)) /     &
                  ((1.-(qv(i)+q_liq+q_sol))*cp_air+qv(i)*cp_vap+q_liq*c_liq+q_sol*c_ice)
     enddo
  else
     do i=is, ie
        q_liq = ql(i) + qr(i)
        q_sol = qi(i) + qs(i) + qg(i)
        lcp2(i) = (Lv0+dc_vap*tz(i)) /   &
                  ((1.-(qv(i)+q_liq+q_sol))*cv_air+qv(i)*cv_vap+q_liq*c_liq+q_sol*c_ice)
     enddo
  endif

  do i=is, ie
     if ( qr(i) > 1.E-7 .and. tz(i) > t_wfr ) then
          qpz = qv(i) + ql(i)
          tin = tz(i) - lcp2(i)*ql(i) ! presence of clouds suppresses the rain evap
         qsat = wqs2(tin, den(i), dqsdt)
          dqh = max( ql(i), hvar(i)*max(qpz, 1.E-7) )
          dqh = min( dqh, 0.2*qpz ) ! New limiter
          dqv = qsat - qv(i)
         q_minus = qpz - dqh
         q_plus  = qpz + dqh

! qsat must be > q_minus to activate evaporation
! qsat must be < q_plus  to activate accretion

!-------------------
! * Rain evaporation
!-------------------
         dqm = qsat - q_minus
         if ( dqv > 1.E-8  .and. dqm > 0. ) then
              if ( qsat > q_plus ) then
                   dq = qsat - qpz
              else
! q_minus < qsat < q_plus
! dq == dqh if qsat == q_minus
                   dq = 0.25*dqm*dqm/dqh
              endif
              qden = qr(i)*den(i)
              t2 = tin * tin
              evap =  crevp(1)*t2*dq*(crevp(2)*sqrt(qden)+crevp(3)*exp(0.725*log(qden)))   &
                   / (crevp(4)*t2  +  crevp(5)*qsat*den(i))
              evap = min( qr(i), dt*evap, dqv/(1.+lcp2(i)*dqsdt) )
!             Minimum Evap of rain in dry environmental air
              evap = max( evap, min(qr(i),dim(rh_rain*qsat, qv(i))/(1.+lcp2(i)*dqsdt)) )
             qr(i) = qr(i) - evap
             qv(i) = qv(i) + evap
             tz(i) = tz(i) - evap*lcp2(i)
         endif
         if ( qr(i)>1.E-7 .and. ql(i)>1.E-6 .and. qsat<q_plus ) then
!-------------------
! * Accretion: pracc
!-------------------
               sink = dt*sqrt(sfcrho/den(i))*cracw*exp(0.95*log(qr(i)*den(i)))
               sink = sink/(1.+sink)*ql(i)
              ql(i) = ql(i) - sink
              qr(i) = qr(i) + sink
         endif
     endif
  enddo

 end subroutine revap_rac1

 real function wqs1(ta, den)
! Pure water phase; universal dry/moist formular using air density
! Input "den" can be either dry or moist air density
  real, intent(in):: ta, den
! local:
  real es, ap1
  real, parameter:: tmin=tice - 160.
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = tablew(it) + (ap1-it)*desw(it)
      wqs1 = es / (rvgas*ta*den)

 end function wqs1

 real function iqs1(ta, den)
! water-ice phase; universal dry/moist formular using air density
! Input "den" can be either dry or moist air density
  real, intent(in):: ta, den
! local:
  real es, ap1
  real, parameter:: tmin=tice - 160.
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = table2(it) + (ap1-it)*des2(it)
      iqs1 = es / (rvgas*ta*den)

 end function iqs1


 real function wqs2(ta, den, dqdt)
! Pure water phase; universal dry/moist formular using air density
! Input "den" can be either dry or moist air density
  real, intent(in):: ta, den
  real, intent(out):: dqdt
! local:
  real es, ap1
  real, parameter:: tmin=tice - 160.
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = tablew(it) + (ap1-it)*desw(it)
      wqs2 = es / (rvgas*ta*den)
        it = ap1 - 0.5
! Finite diff, del_T = 0.1:
      dqdt = 10.*(desw(it) + (ap1-it)*(desw(it+1)-desw(it))) / (rvgas*ta*den)

 end function wqs2

 subroutine wqs2_vect(is, ie, ta, den, wqsat, dqdt)
! Pure water phase; universal dry/moist formular using air density
! Input "den" can be either dry or moist air density
  integer, intent(in):: is, ie
  real, intent(in),  dimension(is:ie):: ta, den
  real, intent(out), dimension(is:ie):: wqsat, dqdt
! local:
  real es, ap1
  real, parameter:: tmin=tice - 160.
  integer i, it

   do i=is, ie
       ap1 = 10.*dim(ta(i), tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = tablew(it) + (ap1-it)*desw(it)
      wqsat(i) = es / (rvgas*ta(i)*den(i))
        it = ap1 - 0.5
! Finite diff, del_T = 0.1:
      dqdt(i) = 10.*(desw(it)+(ap1-it)*(desw(it+1)-desw(it)))/(rvgas*ta(i)*den(i))
   enddo

 end subroutine wqs2_vect



 real function iqs2(ta, den, dqdt)
! water-ice phase; universal dry/moist formular using air density
! Input "den" can be either dry or moist air density
  real, intent(in):: ta, den
  real, intent(out):: dqdt
! local:
  real es, ap1
  real, parameter:: tmin=tice - 160.
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = table2(it) + (ap1-it)*des2(it)
      iqs2 = es / (rvgas*ta*den)
        it = ap1 - 0.5
      dqdt = 10.*(des2(it) + (ap1-it)*(des2(it+1)-des2(it))) / (rvgas*ta*den)

 end function iqs2


 subroutine qs_init(kmp)
 integer, intent(in):: kmp
 integer, parameter:: length=2621
 real, parameter:: rhor  = 1.0e3  ! LFO83
 real, parameter:: vdifu = 2.11e-5
 real, parameter:: tcond = 2.36e-2
 real, parameter:: visk = 1.259e-5
 real, parameter:: hltc = 2.5e6
 real, parameter:: gam290 = 1.827363
 real, parameter:: gam380 = 4.694155
 real, parameter:: alin = 842.0
 !Intercept parameters
 real, parameter:: rnzr = 8.0e6
 real, parameter:: c_cracw = 0.9      ! rain accretion efficiency
 real:: scm3, act2
 integer i

 if (  mp_initialized ) return
      if (is_master()) write(*,*) 'Top layer for GFDL_MP=', kmp

      lat2 = (hlv + hlf) ** 2

      scm3 = (visk/vdifu)**(1./3.)
      act2 = pi * rnzr * rhor
      cracw = c_cracw* pi*rnzr*alin*gam380/(4.*act2**0.95)

      crevp(1) = 2.*pi*vdifu*tcond*rvgas*rnzr
      crevp(2) = 0.78/sqrt(act2)
      crevp(3) = 0.31*scm3*gam290*sqrt(alin/visk)/act2**0.725
      crevp(4) = tcond*rvgas
      crevp(5) = hltc**2*vdifu

! generate es table (dt = 0.1 deg. c)
       allocate ( table (length) )
       allocate ( table2(length) )
       allocate ( tablew(length) )
       allocate (   des2(length) )
       allocate (   desw(length) )

       call qs_table (length )
       call qs_table2(length )
       call qs_tablew(length )

       do i=1,length-1
          des2(i) = max(0., table2(i+1) - table2(i))
          desw(i) = max(0., tablew(i+1) - tablew(i))
       enddo
       des2(length) = des2(length-1)
       desw(length) = desw(length-1)

   mp_initialized = .true.

 end subroutine qs_init

 subroutine qs_table(n)
      integer, intent(in):: n
      real(kind=R_GRID):: esupc(200)
      real(kind=R_GRID):: tmin, tem, esh20
      real(kind=R_GRID):: wice, wh2o, t_ice
      real(kind=R_GRID):: delt=0.1
      integer i

! constants
      t_ice = tice

!  compute es over ice between -160c and 0 c.
      tmin = t_ice - 160.
      do i=1,1600
         tem = tmin+delt*real(i-1)
         table(i) = e00*exp((d2ice*log(tem/t_ice)+Li2*(tem-t_ice)/(tem*t_ice))/rvgas)
      enddo

!  compute es over water between -20c and 102c.
      do  i=1,1221
          tem = 253.16+delt*real(i-1)
          esh20 = e00*exp((dc_vap*log(tem/t_ice)+Lv0*(tem-t_ice)/(tem*t_ice))/rvgas)
          if (i <= 200) then
              esupc(i) = esh20
          else
              table(i+1400) = esh20
          endif
      enddo

!  derive blended es over ice and supercooled water between -20c and 0c
      do i=1,200
         tem  = 253.16+delt*real(i-1)
         wice = 0.05*(t_ice-tem)
         wh2o = 0.05*(tem-253.16)
         table(i+1400) = wice*table(i+1400)+wh2o*esupc(i)
      enddo

 end subroutine qs_table

 subroutine qs_tablew(n)
! Over water
   integer, intent(in):: n
   real(kind=R_GRID), parameter:: delt=0.1
   real(kind=R_GRID):: tmin
   real(kind=R_GRID):: tem0, t_ice, fac1
   integer i

! constants
      t_ice = tice
       tmin = t_ice - 160.
     do i=1,n
        tem0 = tmin + delt*real(i-1)
!  compute es over water
        fac1 = Lv0*(tem0-t_ice) / (tem0*t_ice)
        fac1 = (dc_vap*log(tem0/t_ice)+fac1) / rvgas
        fac1 = e00*exp(fac1)
        tablew(i) = fac1
     enddo

 end subroutine qs_tablew


 subroutine qs_table2(n)
! 2-phase table
  integer, intent(in):: n
  real(kind=R_GRID):: delt=0.1
  real(kind=R_GRID):: tmin
  real(kind=R_GRID):: tem0, tem1, t_ice, fac0, fac1, fac2
  integer:: i, i0, i1

! constants
      t_ice = tice
      tmin = t_ice - 160.

! High-precision computation:
     do i=1,n
        tem0 = tmin+delt*real(i-1)
        fac0 = (tem0-t_ice) / (tem0*t_ice)
        if ( i<= 1600 ) then
!  compute es over ice between -160c and 0 c.
             fac1 = fac0*Li2
             fac2 = (d2ice*log(tem0/t_ice)+fac1) / rvgas
        else
!  compute es over water between 0c and 102c.
             fac1 = fac0*Lv0
             fac2 = (dc_vap*log(tem0/t_ice)+fac1) / rvgas
        endif
        fac2 = e00*exp(fac2)
        table2(i) = fac2
     enddo

!----------
! smoother
!----------
      i0 = 1600;  i1 = 1601
      tem0 = 0.25*(table2(i0-1) + 2.*table(i0) + table2(i0+1))
      tem1 = 0.25*(table2(i1-1) + 2.*table(i1) + table2(i1+1))
      table2(i0) = tem0
      table2(i1) = tem1

 end subroutine qs_table2

end module fv_cmp_mod
