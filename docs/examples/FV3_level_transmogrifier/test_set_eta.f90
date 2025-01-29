program test_set_eta
use fv_eta_mod

    implicit none

    integer, parameter :: km = 32
    real :: ak(km+1)
    real :: bk(km+1)
    character(24) :: npz_type

    npz_type = '32old'
    call set_eta(km, ak, bk, npz_type)
end program
