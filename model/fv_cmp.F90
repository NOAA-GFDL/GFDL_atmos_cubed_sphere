module fv_cmp_mod

  use constants_mod,         only: pi=>pi_8, rvgas, rdgas, grav, hlv, hlf, cp_air, cp_vapor
  use fv_mp_mod,             only: is_master
  use fv_arrays_mod,         only: R_GRID

  implicit none
  real, parameter:: cv_vap = 3.*rvgas  ! 1384.8
  real, parameter:: cv_air =  cp_air - rdgas ! = rdgas * (7/2-1) = 2.5*rdgas=717.68
! 2050 at 0 deg C; 1972 at -15 C; 1818. at -40 C
! real, parameter:: c_ice = 2106.        ! heat capacity of ice at 0.C (same as IFS)
! real, parameter:: c_liq = 4218.        ! ECMWF-IFS at 0 deg C
  real, parameter:: c_ice = 1972.        ! -15 C
  real, parameter:: c_liq = 4.1855e+3    ! GFS, at 15 deg C
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
! Li (T=113) ~ 0.
!!! real(kind=R_GRID), parameter:: e00 = 610.71  ! saturation vapor pressure at T0
 real(kind=R_GRID), parameter:: e00 = 611.21  ! IFS: saturation vapor pressure at T0
 real(kind=R_GRID), parameter:: d2ice  = cp_vap - c_ice
 real(kind=R_GRID), parameter:: Li2 = hlv0+hlf0 - d2ice*tice
! Local:
 real:: dw_ocean = 0.12 ! This parameter is different from that in major MP
 real:: crevp(5), lat2
 real, allocatable:: table(:), table2(:), tablew(:), des2(:), desw(:)
 real:: d0_vap, lv00

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
 real, dimension(is:ie):: wqsat, dq2dt, qpz, cvm, t0, pt1, icp2, lcp2, tcp2, tcp3,    &
                          den, q_liq, q_sol, src, hvar
 real, dimension(is:ie):: mc_air, lhl, lhi  ! latent heat
 real:: sink, qsw, rh, fac_v2l, fac_l2v
 real:: tc, qsi, dqsdt, dq, dq0, pidep, qi_crt, tmp, dtmp
 real:: condensates, tin, qstar, rqi, q_plus, q_minus
 real:: sdt, dt_Bigg, adj_fac, fac_s, fac_r, fac_i2s, fac_mlt, fac_l2r
 real:: factor, qim, tice0, c_air, c_vap
 integer i,j


end subroutine fv_sat_adj


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
