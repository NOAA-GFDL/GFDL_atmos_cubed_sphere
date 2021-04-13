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
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

#ifndef _FV_DIAG__
#define _FV_DIAG__

 integer ::id_ps, id_slp, id_ua, id_va, id_pt, id_omga, id_vort,  &
      id_tm, id_pv, id_zsurf, id_oro, id_sgh, id_w, &
           id_ke, id_zs, id_ze, id_mq, id_vorts, id_us, id_vs,    &
           id_tq, id_rh, id_c15, id_c25, id_c35, id_c45,          &
                         id_f15, id_f25, id_f35, id_f45, id_ctp,  &
           id_ppt, id_ts, id_tb, id_ctt, id_pmask, id_pmaskv2,    &
           id_delp, id_delz, id_iw, id_lw,                 &
           id_pfhy, id_pfnh, id_ppnh,                             &
           id_qn, id_qn200, id_qn500, id_qn850, id_qp,            &
           id_qdt, id_acly, id_acl, id_acl2,                              &
           id_dbz, id_maxdbz, id_basedbz, id_dbz4km, id_dbztop, id_dbz_m10C, &
           id_ctz, id_w1km, id_wmaxup, id_wmaxdn, id_cape, id_cin

! Selected theta-level fields from 3D variables:
 integer :: id_pv350K, id_pv550K

! Selected p-level fields from 3D variables:
 integer :: id_vort200, id_vort500, id_w500, id_w700
 integer :: id_vort850, id_w850, id_x850, id_srh25, &
            id_uh03, id_uh25, id_theta_e,  &
            id_w200, id_s200, id_sl12, id_sl13, id_w5km, id_rain5km, id_w2500m
 integer :: id_srh1, id_srh3, id_ustm, id_vstm
! NGGPS 31-level diag
 integer, allocatable :: id_u(:), id_v(:), id_t(:), id_h(:), id_q(:), id_omg(:)

 integer:: id_u_plev, id_v_plev, id_t_plev, id_h_plev, id_q_plev, id_omg_plev
 integer:: id_t_plev_ave, id_q_plev_ave, id_qv_dt_gfdlmp_plev_ave, id_t_dt_gfdlmp_plev_ave, id_qv_dt_phys_plev_ave, id_t_dt_phys_plev_ave

 ! IPCC diag
 integer :: id_rh10,  id_rh50,  id_rh100, id_rh200,  id_rh250, id_rh300, &
            id_rh500, id_rh700, id_rh850, id_rh925,  id_rh1000
 integer :: id_dp10,  id_dp50,  id_dp100, id_dp200,  id_dp250, id_dp300, &
            id_dp500, id_dp700, id_dp850, id_dp925,  id_dp1000

 integer :: id_rh1000_cmip, id_rh925_cmip, id_rh850_cmip, id_rh700_cmip, id_rh500_cmip, &
            id_rh300_cmip,  id_rh250_cmip, id_rh100_cmip, id_rh50_cmip,  id_rh10_cmip

 integer :: id_hght3d, id_any_hght
 integer :: id_u100m, id_v100m, id_w100m

     ! For initial conditions:
     integer ic_ps, ic_ua, ic_va, ic_ppt
     integer ic_sphum
     integer, allocatable :: id_tracer(:)
! ESM requested diagnostics  -  dry mass/volume mixing ratios
 integer, allocatable :: id_tracer_dmmr(:)
 integer, allocatable :: id_tracer_dvmr(:)
 real,    allocatable :: w_mr(:)

     real, allocatable :: phalf(:)
     real, allocatable :: zsurf(:,:)
     real, allocatable :: pt1(:)

     integer :: id_prer, id_prei, id_pres, id_preg, id_cond, id_dep, id_reevap, id_sub
     integer :: id_qv_dt_gfdlmp, id_T_dt_gfdlmp, id_ql_dt_gfdlmp, id_qi_dt_gfdlmp
     integer :: id_qr_dt_gfdlmp, id_qg_dt_gfdlmp, id_qs_dt_gfdlmp
     integer :: id_liq_wat_dt_gfdlmp, id_ice_wat_dt_gfdlmp
     integer :: id_u_dt_gfdlmp, id_v_dt_gfdlmp
     integer :: id_t_dt_phys, id_qv_dt_phys, id_ql_dt_phys, id_qi_dt_phys, id_u_dt_phys, id_v_dt_phys
     integer :: id_qr_dt_phys, id_qg_dt_phys, id_qs_dt_phys
     integer :: id_liq_wat_dt_phys, id_ice_wat_dt_phys
     integer :: id_intqv, id_intql, id_intqi, id_intqr, id_intqs, id_intqg

! ESM/CM 3-D diagostics
     integer :: id_uq, id_vq, id_wq, id_iuq, id_ivq, id_iwq,   & ! moisture flux & vertical integral
                id_ut, id_vt, id_wt, id_iut, id_ivt, id_iwt,   & ! heat flux
                id_uu, id_uv, id_vv, id_ww,                    & ! momentum flux
                id_iuu, id_iuv, id_iuw, id_ivv, id_ivw, id_iww   ! vertically integral of momentum flux

     integer :: id_uw, id_vw

     integer :: id_t_dt_nudge, id_ps_dt_nudge, id_delp_dt_nudge, id_u_dt_nudge, id_v_dt_nudge

! EMC additions
     integer :: id_diss, id_zratio, id_hw, id_qvw, id_qlw, id_qiw, id_o3w

#endif _FV_DIAG__
