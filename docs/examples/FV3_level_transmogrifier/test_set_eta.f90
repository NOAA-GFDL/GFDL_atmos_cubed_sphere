program test_set_eta
use fv_eta_mod

    implicit none

    integer, parameter :: km = 91
    real :: ak(km+1)
    real :: bk(km+1)
    character(24) :: npz_type
    real :: stretch_fac, pint, ptop, s0, s_k_inc, dzmax, stretch_fac_u
    integer :: k_inc, nfilter
    logical :: doenforce_2nd_order_continuity

    doenforce_2nd_order_continuity = .false.
    stretch_fac=1.029
    pint= 10000.
    ptop=40.
    k_inc=25
    s0=0.13
    s_k_inc=1.2
    nfilter = 0
    dzmax = 2000.
    stretch_fac_u = 1.2

    npz_type = 'wrf'
    call set_eta(km, ak, bk, npz_type, doenforce_2nd_order_continuity, dzmax_in = dzmax, stretch_fac_u_in = stretch_fac_u, stretch_fac_in=stretch_fac, pint_in=pint, ptop_in=ptop, k_inc_in=k_inc, s0_in=s0, s_k_inc_in=s_k_inc, nfilter_in=nfilter)
end program
