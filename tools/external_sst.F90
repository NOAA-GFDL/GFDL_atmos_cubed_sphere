module external_sst_mod

#ifdef NO_GFDL_SHARED
!----------------- Public Data -----------------------------------
integer :: i_sst = -1
integer :: j_sst = -1
logical :: forecast_mode = .false.
logical :: use_ncep_sst  = .false.
real, allocatable, dimension(:,:) ::  sst_ncep, sst_anom
#else
use amip_interp_mod, only: i_sst, j_sst, sst_ncep, sst_anom, &
                           forecast_mode, use_ncep_sst
#endif

public i_sst, j_sst, sst_ncep, sst_anom, forecast_mode, use_ncep_sst

!---- version number -----
character(len=128) :: version = '$Id: external_sst.F90,v 19.0 2012/01/06 19:58:30 fms Exp $'
character(len=128) :: tagname = '$Name: siena $'

end module external_sst_mod
