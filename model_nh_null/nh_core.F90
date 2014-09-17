module nh_core_mod

  use fms_mod, only: error_mesg, FATAL
  use fv_arrays_mod,  only: fv_grid_type, fv_grid_bounds_type


  implicit none
  private

  public Riem_Solver3, Riem_Solver_c, update_dz_c, update_dz_d, geopk_halo_nh

!---- version number -----
  character(len=128) :: version = '$Id: nh_core.F90,v 20.0 2013/12/13 23:07:15 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal_201409 $'

contains

  subroutine update_dz_c(is, ie, js, je, km, ng, dt, dp0, zs, area, ut, vt, gz, ws, &
       npx, npy, sw_corner, se_corner, ne_corner, nw_corner, bd, grid_type)
! !INPUT PARAMETERS:
   type(fv_grid_bounds_type), intent(IN) :: bd
   integer, intent(in):: is, ie, js, je, ng, km, npx, npy, grid_type
   logical, intent(IN):: sw_corner, se_corner, ne_corner, nw_corner
   real, intent(in):: dt
   real, intent(in):: dp0(km)
   real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng,km):: ut, vt
   real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng):: area
   real, intent(inout):: gz(is-ng:ie+ng,js-ng:je+ng,km+1)
   real, intent(in   ):: zs(is-ng:ie+ng, js-ng:je+ng)
   real, intent(  out):: ws(is-ng:ie+ng, js-ng:je+ng)

   call error_mesg('update_dz_c','The null version of update_dz_c should not be called.',FATAL)

  end subroutine update_dz_c



  subroutine update_dz_d(ndif, damp, hord, is, ie, js, je, km, ng, npx, npy, area, rarea,   &
                         dp0, zs, zh, crx, cry, xfx, yfx, delz, ws, rdt, gridstruct, bd)
   type(fv_grid_bounds_type), intent(IN) :: bd
   integer, intent(in):: is, ie, js, je, ng, km, npx, npy
   integer, intent(in):: hord
   real, intent(in)   :: rdt
   real, intent(in)   :: dp0(km)
   real, intent(in)   :: area(is-ng:ie+ng,js-ng:je+ng)
   real, intent(in)   :: rarea(is-ng:ie+ng,js-ng:je+ng)
   real,    intent(inout):: damp(km+1)
   integer, intent(inout):: ndif(km+1)
   real, intent(in   ) ::  zs(is-ng:ie+ng,js-ng:je+ng)
   real, intent(inout) ::  zh(is-ng:ie+ng,js-ng:je+ng,km+1)
   real, intent(  out) ::delz(is-ng:ie+ng,js-ng:je+ng,km)
   real, intent(inout), dimension(is:ie+1,js-ng:je+ng,km):: crx, xfx
   real, intent(inout), dimension(is-ng:ie+ng,js:je+1,km):: cry, yfx
   real, intent(inout)   :: ws(is:ie,js:je)
   type(fv_grid_type), intent(IN), target :: gridstruct

   call error_mesg('update_dz_d','The null version of update_dz_d should not be called.',FATAL)

  end subroutine update_dz_d


  subroutine Riem_Solver3(ms, dt,   is,   ie,   js, je, km, ng,    &
                          isd, ied, jsd, jed,        &
                          akap, cp,   ptop, hs, w,  delz, pt,  &
                          delp, zh, gz,  ppe, pk3, pk, pe, peln, &
                          ws, p_fac, a_imp, &
                          scale_z, use_logp, last_call)
   integer, intent(in):: ms, is, ie, js, je, km, ng
   integer, intent(in):: isd, ied, jsd, jed
   real, intent(in):: dt         ! the BIG horizontal Lagrangian time step
   real, intent(in):: akap, cp, ptop, p_fac, a_imp, scale_z
   real, intent(in):: hs(isd:ied,jsd:jed)
   logical, intent(in):: last_call, use_logp
   real, intent(in):: ws(isd:ied,jsd:jed)
   real, intent(inout), dimension(isd:ied,jsd:jed,km+1):: zh
   real, intent(inout):: peln(is:ie,km+1,js:je)          ! ln(pe)
   real, intent(inout), dimension(isd:ied,jsd:jed,km):: w, delp, pt
   real, intent(out), dimension(isd:ied,jsd:jed,km+1):: ppe
   real, intent(out):: delz(is-ng:ie+ng,js-ng:je+ng,km)
   real, intent(out):: pk(is:ie,js:je,km+1)
   real, intent(out):: pe(is-1:ie+1,km+1,js-1:je+1)
   real, intent(out):: pk3(isd:ied,jsd:jed,km+1)
   real, intent(out), dimension(isd:ied,jsd:jed,km+1):: gz

   call error_mesg('Riem_Solver3','The null version of Riem_Solver3 should not be called.',FATAL)

  end subroutine Riem_Solver3



  subroutine Riem_Solver_c(ms,   dt,  is,   ie,   js, je, km,   ng,  &
                           akap, cp,  ptop, hs, w3,  pt,  &
                           delp, gz,  pef,  ws, p_fac, a_imp, scale_z)
   integer, intent(in):: is, ie, js, je, ng, km
   integer, intent(in):: ms
   real, intent(in):: dt,  akap, cp, ptop, p_fac, a_imp, scale_z
   real, intent(in):: ws(is-ng:ie+ng,js-ng:je+ng)
   real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng,km):: pt, delp
   real, intent(in)::   hs(is-ng:ie+ng,js-ng:je+ng)
   real, intent(in):: w3(is-ng:ie+ng,js-ng:je+ng,km)
! OUTPUT PARAMETERS
   real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: gz
   real, intent(  out), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: pef

   call error_mesg('Riem_Solver_c','The null version of Riem_Solver_c should not be called.',FATAL)

  end subroutine Riem_Solver_c



  subroutine geopk_halo_nh(ptop, grav, kappa, cp, delp, delz, pt, phis, pkc, gz, pk3, &
                           npx, npy, npz, nested, pkc_pertn, computepk3, fullhalo, bd)

   integer, intent(in) :: npx, npy, npz
   logical, intent(in) :: pkc_pertn, computepk3, fullhalo, nested
   real, intent(in) :: ptop, kappa, cp, grav
   type(fv_grid_bounds_type), intent(IN) :: bd
   real, intent(in) :: phis(bd%isd:bd%ied,bd%jsd:bd%jed)
   real, intent(in),  dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz):: pt, delp, delz
   real, intent(inout), dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz+1):: gz, pkc, pk3

   call error_mesg('geopk_halo_nh','The null version of geopk_halo_nh should not be called.',FATAL)

  end subroutine geopk_halo_nh

end module nh_core_mod
