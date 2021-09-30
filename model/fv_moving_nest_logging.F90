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

!----------------------------------------------------------
! Moving Nest Initial Release    W. Ramstrom - 07/28/2021
!----------------------------------------------------------

!***********************************************************************
!>@brief!   Provides subroutines to debug and log moving nest functionality 
!!>@author W. Ramstrom, AOML/HRD   01/15/2021
!
!   This code is in a separate module so that code review and optimization 
!     can focus on the algorithm code in fv_moving_nest.F90 that implements
!     the core functionality.  These routines will likely be disabled or 
!     removed before operational implementation.
! 
! =======================================================================!
!


! Notes


module fv_moving_nest_logging_mod

#ifdef MOVING_NEST

  use mpp_mod,                  only: mpp_pe
  use fv_arrays_mod,            only: R_GRID
  use fv_moving_nest_utils_mod, only: grid_geometry, fill_grid_from_supergrid
  use fv_arrays_mod,            only: fv_grid_type, fv_nest_type, fv_atmos_type

contains

#endif ! MOVING_NEST

end module fv_moving_nest_logging_mod
