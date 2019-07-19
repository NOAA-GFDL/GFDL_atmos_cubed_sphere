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
module fv_cmp_mod

  use constants_mod,         only: pi=>pi_8, rvgas, rdgas, grav, hlv, hlf, cp_air, cp_vapor
  use fms_mod,               only: error_mesg, FATAL
  use fv_mp_mod,             only: is_master
  use fv_arrays_mod,         only: R_GRID

  implicit none

  private
  public fv_sat_adj, qs_init

contains

 subroutine fv_sat_adj(mdt, zvir, is, ie, js, je, ng, hydrostatic, consv_te, &
                       te0, qv, ql, qi, qr, qs, qg, dpln, delz, pt, dp,  &
                       q_con, cappa, area, dtdt, out_dt, last_step, do_qa, qa)
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

   call error_mesg('fv_cmp_mod','saturation adjustment is not available.',FATAL)

 end subroutine fv_sat_adj


 subroutine qs_init(kmp)
 integer, intent(in):: kmp

   call error_mesg('fv_cmp_mod','saturation adjustment is not available.',FATAL)

 end subroutine qs_init

end module fv_cmp_mod
