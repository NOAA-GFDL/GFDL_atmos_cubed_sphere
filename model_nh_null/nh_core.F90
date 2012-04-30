module nh_core_mod

  use fms_mod, only: error_mesg, FATAL

  implicit none
  private

  public Riem_Solver, update_gz_c, update_dz_d

!---- version number -----
  character(len=128) :: version = '$Id: nh_core.F90,v 1.1.2.1 2012/09/20 19:18:11 rab Exp $'
  character(len=128) :: tagname = '$Name: siena_201204 $'

contains

  subroutine update_gz_c(is, ie, js, je, km, ng, area, ut, vt, gz)

  integer, intent(in):: is, ie, js, je, ng, km
  real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng,km):: ut, vt
  real, intent(in):: area(is-ng:ie+ng,js-ng:je+ng)
  real, intent(inout):: gz(is-ng:ie+ng,js-ng:je+ng,km+1)  ! work array

  call error_mesg('update_gz_c','The null version of update_gz_c should not be called.',FATAL)

  end subroutine update_gz_c


  subroutine update_dz_d(ndif, damp, hord, is, ie, js, je, km, ng, npx, npy, area,    &
                         zs, zh, crx, cry, xfx, yfx, delz, ws, rdt, id_ws)

  integer, intent(in):: is, ie, js, je, ng, km, npx, npy
  integer, intent(in):: ndif(km)
  integer, intent(in):: hord
  integer, intent(in):: id_ws
  real, intent(in)   :: rdt
  real, intent(in)   :: damp(km)
  real, intent(in)   :: area(is-ng:ie+ng,js-ng:je+ng)
  real, intent(inout) ::  zs(is-ng:ie+ng,js-ng:je+ng)
  real, intent(inout) ::  zh(is-ng:ie+ng,js-ng:je+ng,km)
  real, intent(inout) ::delz(is:ie,js:je,km)
  real, intent(inout), dimension(is:ie+1,js-ng:je+ng,km):: crx, xfx
  real, intent(inout), dimension(is-ng:ie+ng,js:je+1,km):: cry, yfx
  real, intent(inout)   :: ws(is:ie,js:je)

  call error_mesg('update_dz_d','The null version of update_dz_d should not be called.',FATAL)

  end subroutine update_dz_d


  subroutine Riem_Solver(dt,   is,   ie,   js, je, km, ng,    &
                         akap, cp,   ptop, hs, w,  delz, pt,  &
                         delp, gz,   pkc, pk, pe, peln,       &
                         da_min, diff_z0, first_call, last_call, compute_pvar)

   integer, intent(in):: is, ie, js, je, km, ng
   real, intent(in):: dt         ! the BIG horizontal Lagrangian time step
   real, intent(in):: akap, cp, ptop
   real, intent(in):: da_min, diff_z0
   real, intent(in):: hs(is-ng:ie+ng,js-ng:je+ng)
   logical, intent(in):: first_call, last_call, compute_pvar
   real, intent(inout):: peln(is:ie,km+1,js:je)          ! ln(pe)
   real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,km):: w, delp, pt
   real, intent(inout):: delz(is:ie,js:je,km)
   real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: pkc
   real, intent(out):: pe(is-1:ie+1,km+1,js-1:je+1)
   real, intent(out):: pk(is:ie,js:je,km+1)
   real, intent(out), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: gz

   call error_mesg('Riem_Solver','The null version of Riem_Solver should not be called.',FATAL)

  end subroutine Riem_Solver

end module nh_core_mod
