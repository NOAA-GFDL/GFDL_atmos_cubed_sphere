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

!>@brief The module 'fv_eta' contains routine to set up the reference
!! (Eulerian) pressure coordinate

module fv_eta_mod

! <table>
! <tr>
!     <th>Module Name</th>
!     <th>Functions Included</th>
!   </tr>
!   <tr>
!     <td>constants_mod</td>
!     <td>kappa, grav, cp_air, rdgas</td>
!   </tr>
!   <tr>
!     <td>fv_mp_mod</td>
!     <td>is_master</td>
!   </tr>
!   <tr>
!     <td>mpp_mod</td>
!     <td>mpp_error, FATAL</td>
!   </tr>
! </table>

 use constants_mod,  only: kappa, grav, cp_air, rdgas
 use fv_mp_mod,      only: is_master
 use mpp_mod,        only: FATAL, mpp_error
 implicit none
 private
 public set_eta, set_external_eta, get_eta_level, compute_dz_var,  &
        compute_dz_L32, compute_dz_L101, set_hybrid_z, compute_dz, &
        gw_1d, sm1_edge, hybrid_z_dz

 contains

!!!NOTE: USE_VAR_ETA not used in SHiELD
!!! This routine will be kept here
!!! for the time being to not disrupt idealized tests
#ifdef USE_VAR_ETA
 subroutine set_eta(km, ks, ptop, ak, bk, npz_type)
! This is the easy to use version of the set_eta
      integer,  intent(in)::  km           ! vertical dimension
      integer,  intent(out):: ks           ! number of pure p layers
      real:: a60(61),b60(61)
! Thfollowing L63 setting is the same as NCEP GFS's L64 except the top
! 3 layers
      data a60/300.0000,     430.00000,     558.00000,    &
              700.00000,     863.05803,    1051.07995,    &
             1265.75194,    1510.71101,    1790.05098,    &
             2108.36604,    2470.78817,    2883.03811,    &
             3351.46002,    3883.05187,    4485.49315,    &
             5167.14603,    5937.04991,    6804.87379,    &
             7780.84698,    8875.64338,   10100.20534,    &
            11264.35673,   12190.64366,   12905.42546,    &
            13430.87867,   13785.88765,   13986.77987,    &
            14047.96335,   13982.46770,   13802.40331,    &
            13519.33841,   13144.59486,   12689.45608,    &
            12165.28766,   11583.57006,   10955.84778,    &
            10293.60402,    9608.08306,    8910.07678,    &
             8209.70131,    7516.18560,    6837.69250,    &
             6181.19473,    5552.39653,    4955.72632,    &
             4394.37629,    3870.38682,    3384.76586,    &
             2937.63489,    2528.37666,    2155.78385,    &
             1818.20722,    1513.68173,    1240.03585,    &
              994.99144,     776.23591,     581.48797,    &
              408.53400,     255.26520,     119.70243, 0. /

      data b60/0.00000,       0.00000,       0.00000,    &
               0.00000,       0.00000,       0.00000,    &
               0.00000,       0.00000,       0.00000,    &
               0.00000,       0.00000,       0.00000,    &
               0.00000,       0.00000,       0.00000,    &
               0.00000,       0.00000,       0.00000,    &
               0.00000,       0.00000,       0.00000,    &
               0.00201,       0.00792,       0.01755,    &
               0.03079,       0.04751,       0.06761,    &
               0.09097,       0.11746,       0.14690,    &
               0.17911,       0.21382,       0.25076,    &
               0.28960,       0.32994,       0.37140,    &
               0.41353,       0.45589,       0.49806,    &
               0.53961,       0.58015,       0.61935,    &
               0.65692,       0.69261,       0.72625,    &
               0.75773,       0.78698,       0.81398,    &
               0.83876,       0.86138,       0.88192,    &
               0.90050,       0.91722,       0.93223,    &
               0.94565,       0.95762,       0.96827,    &
               0.97771,       0.98608,       0.99347,  1./
      real, intent(out):: ak(km+1)
      real, intent(out):: bk(km+1)
      real, intent(out):: ptop         ! model top (Pa)
      character(24), intent(IN) :: npz_type
      real pint, stretch_fac
      integer  k
      real :: s_rate = -1.0 ! dummy value to not use var_les

      pint = 100.E2

!- Notes ---------------------------------
!  low-top:  ptop = 100.   ! ~45 km
!  mid-top:  ptop = 10.    ! ~60 km
!  hi -top:  ptop = 1.     ! ~80 km
!-----------------------------------------
      select case (km)

! Optimal number = 8 * N -1 (for  8 openMP threads)
!                 16 * M -1 (for 16 openMP threads)

#ifdef HIWPP
#ifdef SUPER_K
        case (20)
        ptop = 56.e2
        pint = ptop
        stretch_fac = 1.03
        case (24)
        ptop = 56.e2
        pint = ptop
        stretch_fac = 1.03
        case (30)
        ptop = 56.e2
        pint = ptop
        stretch_fac = 1.03
        case (40)
        ptop = 56.e2
        pint = ptop
        stretch_fac = 1.03
        case (50)
        ptop = 56.e2
        pint = ptop
        stretch_fac = 1.03
        case (60)
        ptop = 56.e2
        pint = ptop
        stretch_fac = 1.03
        case (80)
        ptop = 56.e2
        pint = ptop
        stretch_fac = 1.03
#else
        case (30)          ! For Baroclinic Instability Test
             ptop = 2.26e2
             pint = 250.E2
             stretch_fac = 1.03
        case (40)
             ptop = 50.e2   ! For super cell test
             pint = 300.E2
             stretch_fac = 1.03
        case (50)          ! Mountain waves?
             ptop = 30.e2
             stretch_fac = 1.025
        case (60)          ! For Baroclinic Instability Test
#ifdef GFSL60
          ks = 20
          ptop = a60(1)
          pint = a60(ks+1)
          do k=1,km+1
            ak(k) = a60(k)
            bk(k) = b60(k)
          enddo
#else
             ptop = 3.e2
!            pint = 250.E2
             pint = 300.E2    ! revised for Moist test
             stretch_fac = 1.03
#endif
#endif
        case (64)
!!!          ptop = 3.e2
             ptop = 2.0e2
             pint = 300.E2
             stretch_fac = 1.03
#else
! *Very-low top: for idealized super-cell simulation:
        case (50)
             ptop = 50.e2
             pint = 250.E2
             stretch_fac = 1.03
        case (60)
             ptop = 40.e2
             pint = 250.E2
             stretch_fac = 1.03
        case (90)          ! super-duper cell
             ptop = 40.e2
             stretch_fac = 1.025
#endif
! Low-top:
        case (31)               ! N = 4, M=2
             ptop = 100.
             stretch_fac = 1.035
        case (32)               ! N = 4, M=2
             ptop = 100.
             stretch_fac = 1.035
        case (39)               ! N = 5
             ptop = 100.
             stretch_fac = 1.035
        case (41)
             ptop = 100.
             stretch_fac = 1.035
        case (47)               ! N = 6, M=3
             ptop = 100.
             stretch_fac = 1.035
        case (51)
             ptop = 100.
             stretch_fac = 1.03
        case (52)  ! very low top
             ptop = 30.e2    ! for special DPM RCE experiments
             stretch_fac = 1.03
! Mid-top:
        case (55)               ! N = 7
             ptop = 10.
             stretch_fac = 1.035
! Hi-top:
        case (63)               ! N = 8, M=4
             ptop = 1.
                                ! c360 or c384
             stretch_fac = 1.035
        case (71)               ! N = 9
             ptop = 1.
             stretch_fac = 1.03
        case (79)               ! N = 10, M=5
             ptop = 1.
             stretch_fac = 1.03
        case (127)               ! N = 10, M=5
             ptop = 1.
             stretch_fac = 1.03
        case (151)
           ptop = 75.e2
           pint = 500.E2
           s_rate = 1.01
        case default
             ptop = 1.
             stretch_fac = 1.03
      end select

#ifdef MOUNTAIN_WAVES
      call mount_waves(km, ak, bk, ptop, ks, pint)
#else
      if (s_rate > 0.) then
           call var_les(km, ak, bk, ptop, ks, pint, s_rate)
      else
         if ( km > 79 ) then
            call var_hi2(km, ak, bk, ptop, ks, pint, stretch_fac)
         elseif (km==5 .or. km==10 ) then
! Equivalent Shallow Water: for NGGPS, variable-resolution testing
            ptop = 500.e2
            ks = 0
            do k=1,km+1
               bk(k) = real(k-1) / real (km)
               ak(k) = ptop*(1.-bk(k))
            enddo
         else
#ifndef GFSL60
            call var_hi(km, ak, bk, ptop, ks, pint, stretch_fac)
#endif
         endif
#endif
      endif

      ptop = ak(1)
      pint = ak(ks+1)

 end subroutine set_eta


#else
 !This is the version of set_eta used in SHiELD and AM4
 subroutine set_eta(km, ks, ptop, ak, bk, npz_type)

!Level definitions are now in this header file
#include <tools/fv_eta.h>

   integer,  intent(in)::  km           ! vertical dimension
   integer,  intent(out):: ks           ! number of pure p layers
   real, intent(out):: ak(km+1)
   real, intent(out):: bk(km+1)
   real, intent(out):: ptop         ! model top (Pa)
   character(24), intent(IN) :: npz_type

   real:: p0=1000.E2
   real:: pc=200.E2

   real pt, lnpe, dlnp
   real press(km+1), pt1(km)
   integer  k
   integer :: var_fn = 0

   real :: pint = 100.E2
   real :: stretch_fac = 1.03
   integer :: auto_routine = 0


   ptop = 1.

   ! Definition: press(i,j,k) = ak(k) + bk(k) * ps(i,j)

   if (trim(npz_type) == 'superC' .or. trim(npz_type) == 'superK') then

      auto_routine = 1
      select case (km)
      case (20)
         ptop = 56.e2
         pint = ptop
         stretch_fac = 1.03
      case (24)
         ptop = 56.e2
         pint = ptop
         stretch_fac = 1.03
      case (30)
         ptop = 56.e2
         pint = ptop
         stretch_fac = 1.03
      case (40)
         ptop = 56.e2
         pint = ptop
         stretch_fac = 1.03
      case (50)
         ptop = 56.e2
         pint = ptop
         stretch_fac = 1.03
      case (60)
         ptop = 56.e2
         pint = ptop
         stretch_fac = 1.03
      case (80)
         ptop = 56.e2
         pint = ptop
         stretch_fac = 1.03
      case (90)          ! super-duper cell
         ptop = 40.e2
         stretch_fac = 1.025
         auto_routine = 2
      end select

   else

      select case (km)

      case (5,10) ! does this work????

         ! Equivalent Shallow Water: for modon test
         ptop = 500.e2
         ks = 0
         do k=1,km+1
            bk(k) = real(k-1) / real (km)
            ak(k) = ptop*(1.-bk(k))
         enddo

      case (24)

         ks = 5
         do k=1,km+1
            ak(k) = a24(k)
            bk(k) = b24(k)
         enddo

      case (26)

         ks = 7
         do k=1,km+1
            ak(k) = a26(k)
            bk(k) = b26(k)
         enddo

      case (30)          ! For Baroclinic Instability Test
         ptop = 2.26e2
         pint = 250.E2
         stretch_fac = 1.03
         auto_routine = 1

      case (31)               ! N = 4, M=2
         if (trim(npz_type) == 'lowtop') then
            ptop = 300.
            pint = 100.E2
            stretch_fac = 1.035
            auto_routine = 5
         else
            ptop = 100.
            stretch_fac = 1.035
            auto_routine = 1
         endif

      case (32)

         if (trim(npz_type) == 'old32') then
            ks = 13              ! high-res trop_32 setup
            do k=1,km+1
               ak(k) = a32old(k)
               bk(k) = b32old(k)
            enddo
         elseif (trim(npz_type) == 'lowtop') then
            ptop = 100.
            stretch_fac = 1.035
            auto_routine = 1
         else
            ks = 7
            do k=1,km+1
               ak(k) = a32(k)
               bk(k) = b32(k)
            enddo
         endif
         !miz
      case (33)
         ks = 7
         do k=1,km+1
            ak(k) = a33(k)
            bk(k) = b33(k)
         enddo
         !miz

      case (39)               ! N = 5
         ptop = 100.
         stretch_fac = 1.035
         auto_routine = 1

      case (40)
         ptop = 50.e2   ! For super cell test
         pint = 300.E2
         stretch_fac = 1.03
         auto_routine = 1

      case (41)
         ptop = 100.
         pint = 100.E2
         stretch_fac = 1.035
         auto_routine = 1

      case (47)

         if (trim(npz_type) == 'lowtop') then
            ptop = 100.
            stretch_fac = 1.035
            auto_routine = 1
         else
            !         ks = 27       ! high-res trop-strat
            ks = 20       ! Oct 23, 2012
            do k=1,km+1
               ak(k) = a47(k)
               bk(k) = b47(k)
            enddo
         endif

      case (48)
         ks = 28
         do k=1,km+1
            ak(k) = a48(k)
            bk(k) = b48(k)
         enddo

      case (50)
         ! ! *Very-low top: for idealized super-cell simulation:
         ! ptop = 50.e2
         ! pint = 250.E2
         ! stretch_fac = 1.03
         ! auto_routine = 1
         ks = 19
         do k=1,km+1
            ak(k) = a50(k)
            bk(k) = b50(k)
         enddo

      case (51)
         if (trim(npz_type) == 'lowtop') then
            ptop = 100.
            stretch_fac = 1.03
            auto_routine = 1
         elseif (trim(npz_type) == 'meso') then
            ptop = 20.E2
            pint = 100.E2
            stretch_fac = 1.05
            auto_routine = 1
         elseif (trim(npz_type) == 'meso2') then
            ptop = 1.E2
            pint = 100.E2
            stretch_fac = 1.05
            auto_routine = 6
         else
            ptop = 100.
            pint = 100.E2
            stretch_fac = 1.035
            auto_routine = 1
         endif

      case (52)

         if (trim(npz_type) == 'rce') then
            ptop = 30.e2    ! for special DPM RCE experiments
            stretch_fac = 1.03
            auto_routine = 1
         else
            ks = 35         ! pint = 223
            do k=1,km+1
               ak(k) = a52(k)
               bk(k) = b52(k)
            enddo
         endif

      case (54)
         ks = 11         ! pint =  109.4
         do k=1,km+1
            ak(k) = a54(k)
            bk(k) = b54(k)
         enddo

         ! Mid-top:
      case (55)               ! N = 7
         ptop = 10.
         pint = 100.E2
         stretch_fac = 1.035
         auto_routine = 1

      case (56)
         ks = 26
         do k=1,km+1
            ak(k) = a56(k)
            bk(k) = b56(k)
         enddo

      case (60)

         if (trim(npz_type) == 'gfs') then
            ks = 20
            do k=1,km+1
               ak(k) = a60gfs(k)
               bk(k) = b60gfs(k)
            enddo
         else if (trim(npz_type) == 'BCwave') then
            ptop = 3.e2
            !            pint = 250.E2
            pint = 300.E2    ! revised for Moist test
            stretch_fac = 1.03
            auto_routine = 1
         else if (trim(npz_type) == 'meso') then

            ptop = 40.e2
            pint = 250.E2
            stretch_fac = 1.03
            auto_routine = 1

         else
            ks = 19
            do k=1,km+1
               ak(k) = a60(k)
               bk(k) = b60(k)
            enddo
         endif

      case (63)
         if (trim(npz_type) == 'meso') then
            ks = 11
            do k=1,km+1
               ak(k) = a63meso(k)
               bk(k) = b63meso(k)
            enddo
         elseif (trim(npz_type) == 'hitop') then
            ptop = 1.   ! high top
            pint = 100.E2
            stretch_fac = 1.035
            auto_routine = 1
         else!if (trim(npz_type) == 'gfs') then
            !Used for SHiELD
            ! GFS L64 equivalent setting
            ks = 23
            do k=1,km+1
               ak(k) = a63(k)
               bk(k) = b63(k)
            enddo
         endif

      case (64)

         if (trim(npz_type) == 'gfs') then
            ks = 23
            do k=1,km+1
               ak(k) = a64gfs(k)
               bk(k) = b64gfs(k)
            enddo

         else

            ks = 46
            do k=1,km+1
               ak(k) = a64(k)
               bk(k) = b64(k)
            enddo

         endif
         !-->cjg
      case (68)
         ks = 27
         do k=1,km+1
            ak(k) = a68(k)
            bk(k) = b68(k)
         enddo

      case (71)               ! N = 9
         ptop = 1.
         stretch_fac = 1.03
         auto_routine = 1
      case (75)   ! HS-SGO test configuration
         pint = 100.E2
         ptop = 10.E2
         stretch_fac = 1.035
         auto_routine = 6
      case (79)               ! N = 10, M=5
         if (trim(npz_type) == 'gcrm') then
           pint = 100.E2
           ptop = 3.E2
           stretch_fac = 1.035
           auto_routine = 6
         else
           ptop = 1.
           stretch_fac = 1.03
           auto_routine = 1
         endif
      case (90)          ! super-duper cell
         ptop = 40.e2
         stretch_fac = 1.025
         auto_routine = 2

         ! NGGPS_GFS
      case (91)
         pint = 100.E2
         ptop = 40.
         stretch_fac = 1.029
         auto_routine = 6

      case (95)
         ! Mid-top settings:
         pint = 100.E2
         ptop = 20.
         stretch_fac = 1.029
         auto_routine = 6

      case (96)
         ks = 27
         do k=1,km+1
            ak(k) = a96(k)
            bk(k) = b96(k)
         enddo
         !<--cjg

      case (100)
         ks = 38
         do k=1,km+1
            ak(k) = a100(k)
            bk(k) = b100(k)
         enddo

      case (104)
         ks = 73
         do k=1,km+1
            ak(k) = a104(k)
            bk(k) = b104(k)
         enddo

         ! IFS-like L125
      case (125)
         ks = 33
         ptop = a125(1)
         pint = a125(ks+1)
         do k=1,km+1
            ak(k) = a125(k)
            bk(k) = b125(k)
         enddo

      case (127)               ! N = 10, M=5
         if (trim(npz_type) == 'hitop') then
            ptop = 1.
            stretch_fac = 1.03
            auto_routine = 2
         else
            ptop = 1.
            pint = 75.E2
            stretch_fac = 1.028
            auto_routine = 6
         endif
      case (151)
         !LES applications
         ptop = 75.e2
         pint = 500.E2
         stretch_fac = 1.01
         auto_routine = 3

      case default

         if(trim(npz_type) == 'hitop') then
            ptop = 1.
            pint = 100.E2
         elseif(trim(npz_type) == 'midtop') then
            ptop = 10.
            pint = 100.E2
         elseif(trim(npz_type) == 'lowtop') then
            ptop = 1.E2
            pint = 100.E2
         endif

         if (trim(npz_type) == 'gfs') then
            auto_routine = 6
         elseif(trim(npz_type) == 'les') then
            auto_routine = 3
         elseif(trim(npz_type) == 'mountain_wave') then
            auto_routine = 4
         elseif (km > 79) then
            auto_routine = 2
         else
            auto_routine = 1
         endif

      end select

   endif ! superC/superK

   select case (auto_routine)

   case (1)
      call var_hi(km, ak, bk, ptop, ks, pint, stretch_fac)
   case (2)
      call var_hi2(km, ak, bk, ptop, ks, pint, stretch_fac)
   case (3)
      call var_les(km, ak, bk, ptop, ks, pint, stretch_fac)
   case (4)
      call mount_waves(km, ak, bk, ptop, ks, pint)
   case (5)
      call var_dz(km, ak, bk, ptop, ks, pint, stretch_fac)
   case (6)
      call var_gfs(km, ak, bk, ptop, ks, pint, stretch_fac)
   end select

   ptop = ak(1)
   pint = ak(ks+1)

   if (is_master()) then
      write(*, '(A4, A13, A13, A11)') 'klev', 'ak', 'bk', 'p_ref'
      do k=1,km+1
         write(*,'(I4, F13.5, F13.5, F11.2)') k, ak(k), bk(k), 1000.E2*bk(k) + ak(k)
      enddo
   endif


 end subroutine set_eta
#endif

!>@brief The subroutine 'set_external_eta' sets 'ptop' (model top) and
!! 'ks' (first level of pure pressure coordinates given the coefficients
!! 'ak' and 'bk'
 subroutine set_external_eta(ak, bk, ptop, ks)
   implicit none
   real,    intent(in)  :: ak(:)
   real,    intent(in)  :: bk(:)
   real,    intent(out) :: ptop         !< model top (Pa)
   integer, intent(out) :: ks           !< number of pure p layers
   !--- local variables
   integer :: k
   real :: eps = 1.d-7

   ptop = ak(1)
   ks = 1
   do k = 1, size(bk(:))
     if (bk(k).lt.eps) ks = k
   enddo
   !--- change ks to layers from levels
   ks = ks - 1
   if (is_master()) write(6,*) ' ptop & ks ', ptop, ks

 end subroutine set_external_eta


 subroutine var_les(km, ak, bk, ptop, ks, pint, s_rate)
 implicit none
  integer, intent(in):: km
  real,    intent(in):: ptop
  real,    intent(in):: s_rate        !< between [1. 1.1]
  real,    intent(out):: ak(km+1), bk(km+1)
  real,    intent(inout):: pint
  integer, intent(out):: ks
! Local
  real, parameter:: p00 = 1.E5
  real, dimension(km+1):: ze, pe1, peln, eta
  real, dimension(km):: dz, s_fac, dlnp, pm, dp, dk
  real ztop, t0, dz0, sum1, tmp1
  real ep, es, alpha, beta, gama
  real, parameter:: akap = 2./7.
!---- Tunable parameters:
  real:: k_inc = 10   !< number of layers from bottom up to near const dz region
  real:: s0 = 0.8     !< lowest layer stretch factor
!-----------------------
  real:: s_inc
  integer  k

     pe1(1) = ptop
     peln(1) = log(pe1(1))
     pe1(km+1) = p00
     peln(km+1) = log(pe1(km+1))

     t0 = 273.
     ztop = rdgas/grav*t0*(peln(km+1) - peln(1))

      s_inc = (1.-s0) / real(k_inc)
      s_fac(km)  = s0

      do k=km-1, km-k_inc, -1
         s_fac(k)  = s_fac(k+1) + s_inc
      enddo

      s_fac(km-k_inc-1) = 0.5*(s_fac(km-k_inc) + s_rate)

      do k=km-k_inc-2, 5, -1
         s_fac(k) = s_rate * s_fac(k+1)
      enddo

      s_fac(4) = 0.5*(1.1+s_rate)*s_fac(5)
      s_fac(3) = 1.1 *s_fac(4)
      s_fac(2) = 1.1 *s_fac(3)
      s_fac(1) = 1.1 *s_fac(2)

      sum1 = 0.
      do k=1,km
         sum1 = sum1 + s_fac(k)
      enddo

      dz0 = ztop / sum1

      do k=1,km
         dz(k) = s_fac(k) * dz0
      enddo

      ze(km+1) = 0.
      do k=km,1,-1
         ze(k) = ze(k+1) + dz(k)
      enddo

! Re-scale dz with the stretched ztop
      do k=1,km
         dz(k) = dz(k) * (ztop/ze(1))
      enddo

      do k=km,1,-1
         ze(k) = ze(k+1) + dz(k)
      enddo
!     ze(1) = ztop

      if ( is_master() ) then
           write(*,*) 'var_les: computed model top (m)=', ztop, ' bottom/top dz=', dz(km), dz(1)
!           do k=1,km
!              write(*,*) k, s_fac(k)
!           enddo
      endif

      call sm1_edge(1, 1, 1, 1, km, 1, 1, ze, 2)

! Given z --> p
      do k=1,km
          dz(k) = ze(k) - ze(k+1)
        dlnp(k) = grav*dz(k) / (rdgas*t0)
        !write(*,*) k, dz(k)
      enddo
      do k=2,km
         peln(k) = peln(k-1) + dlnp(k-1)
          pe1(k) = exp(peln(k))
      enddo


! Pe(k) = ak(k) + bk(k) * PS
! Locate pint and KS
      ks = 0
      do k=2,km
         if ( pint < pe1(k)) then
              ks = k-1
              exit
         endif
      enddo
      if ( is_master() ) then
         write(*,*) 'For (input) PINT=', 0.01*pint, ' KS=', ks, 'pint(computed)=', 0.01*pe1(ks+1)
      endif
      pint = pe1(ks+1)

      do k=1,km+1
         eta(k) = pe1(k) / pe1(km+1)
      enddo

      ep =  eta(ks+1)
      es =  eta(km)
!     es =  1.
      alpha = (ep**2-2.*ep*es) / (es-ep)**2
      beta  = 2.*ep*es**2 / (es-ep)**2
      gama = -(ep*es)**2 / (es-ep)**2

! Pure pressure:
      do k=1,ks+1
         ak(k) = eta(k)*1.e5
         bk(k) = 0.
      enddo

      do k=ks+2, km
         ak(k) = alpha*eta(k) + beta + gama/eta(k)
         ak(k) = ak(k)*1.e5
      enddo
         ak(km+1) = 0.

      do k=ks+2, km
         bk(k) = (pe1(k) - ak(k))/pe1(km+1)
      enddo
         bk(km+1) = 1.

      if ( is_master() ) then
 !         write(*,*) 'KS=', ks, 'PINT (mb)=', pint/100.
 !         do k=1,km
 !            pm(k) = 0.5*(pe1(k)+pe1(k+1))/100.
 !            write(*,*) k, pm(k), dz(k)
 !         enddo
          tmp1 = ak(ks+1)
          do k=ks+1,km
             tmp1 = max(tmp1, (ak(k)-ak(k+1))/max(1.E-5, (bk(k+1)-bk(k))) )
          enddo
          write(*,*) 'Hybrid Sigma-P: minimum allowable surface pressure (hpa)=', tmp1/100.
          write(*,800) (pm(k), k=km,1,-1)
      endif

    do k=1,km
       dp(k) = (pe1(k+1) - pe1(k))/100.
       dk(k) =  pe1(k+1)**akap - pe1(k)**akap
    enddo

800 format(1x,5(1x,f9.4))

 end subroutine var_les



 subroutine var_gfs(km, ak, bk, ptop, ks, pint, s_rate)
  integer, intent(in):: km
  real,    intent(in):: ptop
  real,    intent(in):: s_rate        !< between [1. 1.1]
  real,    intent(out):: ak(km+1), bk(km+1)
  real,    intent(inout):: pint
  integer, intent(out):: ks
! Local
  real, parameter:: p00 = 1.E5
  real, dimension(km+1):: ze, pe1, peln, eta
  real, dimension(km):: dz, s_fac, dlnp
  real ztop, t0, dz0, sum1, tmp1
  real ep, es, alpha, beta, gama
!---- Tunable parameters:
  integer:: k_inc = 25   !< number of layers from bottom up to near const dz region
  real:: s0 = 0.13 !< lowest layer stretch factor
!-----------------------
  real:: s_inc
  integer  k

     pe1(1) = ptop
     peln(1) = log(pe1(1))
     pe1(km+1) = p00
     peln(km+1) = log(pe1(km+1))

     t0 = 270.
     ztop = rdgas/grav*t0*(peln(km+1) - peln(1))

      s_inc = (1.-s0) / real(k_inc)
      s_fac(km)  = s0

      do k=km-1, km-k_inc, -1
         s_fac(k)  = s_fac(k+1) + s_inc
      enddo

      do k=km-k_inc-1, 9, -1
         s_fac(k) = s_rate * s_fac(k+1)
      enddo
      s_fac(8) = 0.5*(1.1+s_rate)*s_fac(9)
      s_fac(7) = 1.10*s_fac(8)
      s_fac(6) = 1.15*s_fac(7)
      s_fac(5) = 1.20*s_fac(6)
      s_fac(4) = 1.26*s_fac(5)
      s_fac(3) = 1.33*s_fac(4)
      s_fac(2) = 1.41*s_fac(3)
      s_fac(1) = 1.60*s_fac(2)

      sum1 = 0.
      do k=1,km
         sum1 = sum1 + s_fac(k)
      enddo

      dz0 = ztop / sum1

      do k=1,km
         dz(k) = s_fac(k) * dz0
      enddo

      ze(km+1) = 0.
      do k=km,1,-1
         ze(k) = ze(k+1) + dz(k)
      enddo

! Re-scale dz with the stretched ztop
      do k=1,km
         dz(k) = dz(k) * (ztop/ze(1))
      enddo

      do k=km,1,-1
         ze(k) = ze(k+1) + dz(k)
      enddo
!     ze(1) = ztop

      if ( is_master() ) then
           write(*,*) 'var_gfs: computed model top (m)=', ztop*0.001, ' bottom/top dz=', dz(km), dz(1)
!          do k=1,km
!             write(*,*) k, s_fac(k)
!          enddo
      endif

!     call sm1_edge(1, 1, 1, 1, km, 1, 1, ze, 3)

! Given z --> p
      do k=1,km
          dz(k) = ze(k) - ze(k+1)
        dlnp(k) = grav*dz(k) / (rdgas*t0)
      enddo
      do k=2,km
         peln(k) = peln(k-1) + dlnp(k-1)
          pe1(k) = exp(peln(k))
      enddo

! Pe(k) = ak(k) + bk(k) * PS
! Locate pint and KS
      ks = 0
      do k=2,km
         if ( pint < pe1(k)) then
              ks = k-1
              exit
         endif
      enddo
      if ( is_master() ) then
         write(*,*) 'For (input) PINT=', 0.01*pint, ' KS=', ks, 'pint(computed)=', 0.01*pe1(ks+1)
         write(*,*) 'ptop =', ptop
      endif
      pint = pe1(ks+1)

#ifdef NO_UKMO_HB
      do k=1,ks+1
         ak(k) = pe1(k)
         bk(k) = 0.
      enddo

      do k=ks+2,km+1
         bk(k) = (pe1(k) - pint) / (pe1(km+1)-pint)  ! bk == sigma
         ak(k) =  pe1(k) - bk(k) * pe1(km+1)
      enddo
      bk(km+1) = 1.
      ak(km+1) = 0.
#else
! Problematic for non-hydrostatic
      do k=1,km+1
         eta(k) = pe1(k) / pe1(km+1)
      enddo

      ep =  eta(ks+1)
      es =  eta(km)
!     es =  1.
      alpha = (ep**2-2.*ep*es) / (es-ep)**2
      beta  = 2.*ep*es**2 / (es-ep)**2
      gama = -(ep*es)**2 / (es-ep)**2

! Pure pressure:
      do k=1,ks+1
         ak(k) = eta(k)*1.e5
         bk(k) = 0.
      enddo

      do k=ks+2, km
         ak(k) = alpha*eta(k) + beta + gama/eta(k)
         ak(k) = ak(k)*1.e5
      enddo
         ak(km+1) = 0.

      do k=ks+2, km
         bk(k) = (pe1(k) - ak(k))/pe1(km+1)
      enddo
         bk(km+1) = 1.
#endif

      if ( is_master() ) then
          write(*,*) 'KS=', ks, 'PINT (mb)=', pint/100.
          do k=1,km
             write(*,*) k, 0.5*(pe1(k)+pe1(k+1))/100., dz(k)
          enddo
          tmp1 = ak(ks+1)
          do k=ks+1,km
             tmp1 = max(tmp1, (ak(k)-ak(k+1))/max(1.E-5, (bk(k+1)-bk(k))) )
          enddo
          write(*,*) 'Hybrid Sigma-P: minimum allowable surface pressure (hpa)=', tmp1/100.
      endif

 end subroutine var_gfs

 subroutine var_hi(km, ak, bk, ptop, ks, pint, s_rate)
  integer, intent(in):: km
  real,    intent(in):: ptop
  real,    intent(in):: s_rate        !< between [1. 1.1]
  real,    intent(out):: ak(km+1), bk(km+1)
  real,    intent(inout):: pint
  integer, intent(out):: ks
! Local
  real, parameter:: p00 = 1.E5
  real, dimension(km+1):: ze, pe1, peln, eta
  real, dimension(km):: dz, s_fac, dlnp
  real ztop, t0, dz0, sum1, tmp1
  real ep, es, alpha, beta, gama
!---- Tunable parameters:
  integer:: k_inc = 15   !<number of layers from bottom up to near const dz region
  real:: s0 = 0.10 !< lowest layer stretch factor
!-----------------------
  real:: s_inc
  integer  k

     pe1(1) = ptop
     peln(1) = log(pe1(1))
     pe1(km+1) = p00
     peln(km+1) = log(pe1(km+1))

     t0 = 270.
     ztop = rdgas/grav*t0*(peln(km+1) - peln(1))

      s_inc = (1.-s0) / real(k_inc)
      s_fac(km)  = s0

      do k=km-1, km-k_inc, -1
         s_fac(k)  = s_fac(k+1) + s_inc
      enddo

      s_fac(km-k_inc-1) = 0.5*(s_fac(km-k_inc) + s_rate)

#ifdef HIWPP
      do k=km-k_inc-2, 4, -1
         s_fac(k) = s_rate * s_fac(k+1)
      enddo
      s_fac(3) = 0.5*(1.15+s_rate)*s_fac(4)
      s_fac(2) = 1.15 *s_fac(3)
      s_fac(1) = 1.3 *s_fac(2)
#else
      do k=km-k_inc-2, 9, -1
         s_fac(k) = s_rate * s_fac(k+1)
      enddo

      s_fac(8) = 0.5*(1.1+s_rate)*s_fac(9)
      s_fac(7) = 1.1 *s_fac(8)
      s_fac(6) = 1.15*s_fac(7)
      s_fac(5) = 1.2 *s_fac(6)
      s_fac(4) = 1.3 *s_fac(5)
      s_fac(3) = 1.4 *s_fac(4)
      s_fac(2) = 1.45 *s_fac(3)
      s_fac(1) = 1.5 *s_fac(2)
#endif

      sum1 = 0.
      do k=1,km
         sum1 = sum1 + s_fac(k)
      enddo

      dz0 = ztop / sum1

      do k=1,km
         dz(k) = s_fac(k) * dz0
      enddo

      ze(km+1) = 0.
      do k=km,1,-1
         ze(k) = ze(k+1) + dz(k)
      enddo

! Re-scale dz with the stretched ztop
      do k=1,km
         dz(k) = dz(k) * (ztop/ze(1))
      enddo

      do k=km,1,-1
         ze(k) = ze(k+1) + dz(k)
      enddo
!     ze(1) = ztop

      if ( is_master() ) then
           write(*,*) 'var_hi: computed model top (m)=', ztop*0.001, ' bottom/top dz=', dz(km), dz(1)
!          do k=1,km
!             write(*,*) k, s_fac(k)
!          enddo
      endif

      call sm1_edge(1, 1, 1, 1, km, 1, 1, ze, 1)

! Given z --> p
      do k=1,km
          dz(k) = ze(k) - ze(k+1)
        dlnp(k) = grav*dz(k) / (rdgas*t0)
      enddo
      do k=2,km
         peln(k) = peln(k-1) + dlnp(k-1)
          pe1(k) = exp(peln(k))
      enddo

! Pe(k) = ak(k) + bk(k) * PS
! Locate pint and KS
      ks = 0
      do k=2,km
         if ( pint < pe1(k)) then
              ks = k-1
              exit
         endif
      enddo
      if ( is_master() ) then
         write(*,*) 'For (input) PINT=', 0.01*pint, ' KS=', ks, 'pint(computed)=', 0.01*pe1(ks+1)
         write(*,*) 'ptop =', ptop
      endif
      pint = pe1(ks+1)

#ifdef NO_UKMO_HB
      do k=1,ks+1
         ak(k) = pe1(k)
         bk(k) = 0.
      enddo

      do k=ks+2,km+1
         bk(k) = (pe1(k) - pint) / (pe1(km+1)-pint)  ! bk == sigma
         ak(k) =  pe1(k) - bk(k) * pe1(km+1)
      enddo
      bk(km+1) = 1.
      ak(km+1) = 0.
#else
! Problematic for non-hydrostatic
      do k=1,km+1
         eta(k) = pe1(k) / pe1(km+1)
      enddo

      ep =  eta(ks+1)
      es =  eta(km)
!     es =  1.
      alpha = (ep**2-2.*ep*es) / (es-ep)**2
      beta  = 2.*ep*es**2 / (es-ep)**2
      gama = -(ep*es)**2 / (es-ep)**2

! Pure pressure:
      do k=1,ks+1
         ak(k) = eta(k)*1.e5
         bk(k) = 0.
      enddo

      do k=ks+2, km
         ak(k) = alpha*eta(k) + beta + gama/eta(k)
         ak(k) = ak(k)*1.e5
      enddo
         ak(km+1) = 0.

      do k=ks+2, km
         bk(k) = (pe1(k) - ak(k))/pe1(km+1)
      enddo
         bk(km+1) = 1.
#endif

      if ( is_master() ) then
          write(*,*) 'KS=', ks, 'PINT (mb)=', pint/100.
          do k=1,km
             write(*,*) k, 0.5*(pe1(k)+pe1(k+1))/100., dz(k)
          enddo
          tmp1 = ak(ks+1)
          do k=ks+1,km
             tmp1 = max(tmp1, (ak(k)-ak(k+1))/max(1.E-5, (bk(k+1)-bk(k))) )
          enddo
          write(*,*) 'Hybrid Sigma-P: minimum allowable surface pressure (hpa)=', tmp1/100.
      endif


 end subroutine var_hi
 subroutine var_hi2(km, ak, bk, ptop, ks, pint, s_rate)
  integer, intent(in):: km
  real,    intent(in):: ptop
  real,    intent(in):: s_rate        !< between [1. 1.1]
  real,    intent(out):: ak(km+1), bk(km+1)
  real,    intent(inout):: pint
  integer, intent(out):: ks
! Local
  real, parameter:: p00 = 1.E5
  real, dimension(km+1):: ze, pe1, peln, eta
  real, dimension(km):: dz, s_fac, dlnp
  real ztop, t0, dz0, sum1, tmp1
  real ep, es, alpha, beta, gama
  integer  k

     pe1(1) = ptop
     peln(1) = log(pe1(1))
     pe1(km+1) = p00
     peln(km+1) = log(pe1(km+1))

     t0 = 270.
     ztop = rdgas/grav*t0*(peln(km+1) - peln(1))

      s_fac(km  ) = 0.15
      s_fac(km-1) = 0.20
      s_fac(km-2) = 0.30
      s_fac(km-3) = 0.40
      s_fac(km-4) = 0.50
      s_fac(km-5) = 0.60
      s_fac(km-6) = 0.70
      s_fac(km-7) = 0.80
      s_fac(km-8) = 0.90
      s_fac(km-9) = 0.95
      s_fac(km-10) = 0.5*(s_fac(km-9) + s_rate)

      do k=km-11, 8, -1
         s_fac(k) = s_rate * s_fac(k+1)
      enddo

      s_fac(7) = 0.5*(1.1+s_rate)*s_fac(9)
      s_fac(6) = 1.05*s_fac(7)
      s_fac(5) = 1.1*s_fac(6)
      s_fac(4) = 1.15*s_fac(5)
      s_fac(3) = 1.2*s_fac(4)
      s_fac(2) = 1.3*s_fac(3)
      s_fac(1) = 1.4*s_fac(2)

      sum1 = 0.
      do k=1,km
         sum1 = sum1 + s_fac(k)
      enddo

      dz0 = ztop / sum1

      do k=1,km
         dz(k) = s_fac(k) * dz0
      enddo

      ze(km+1) = 0.
      do k=km,1,-1
         ze(k) = ze(k+1) + dz(k)
      enddo

! Re-scale dz with the stretched ztop
      do k=1,km
         dz(k) = dz(k) * (ztop/ze(1))
      enddo

      do k=km,1,-1
         ze(k) = ze(k+1) + dz(k)
      enddo
!     ze(1) = ztop

      if ( is_master() ) write(*,*) 'var_hi2: computed model top (m)=', ztop*0.001, ' bottom/top dz=', dz(km), dz(1)
      call sm1_edge(1, 1, 1, 1, km, 1, 1, ze, 1)

! Given z --> p
      do k=1,km
          dz(k) = ze(k) - ze(k+1)
        dlnp(k) = grav*dz(k) / (rdgas*t0)
      enddo
      do k=2,km
         peln(k) = peln(k-1) + dlnp(k-1)
          pe1(k) = exp(peln(k))
      enddo

! Pe(k) = ak(k) + bk(k) * PS
! Locate pint and KS
      ks = 0
      do k=2,km
         if ( pint < pe1(k)) then
              ks = k-1
              exit
         endif
      enddo
      if ( is_master() ) then
         write(*,*) 'For (input) PINT=', 0.01*pint, ' KS=', ks, 'pint(computed)=', 0.01*pe1(ks+1)
      endif
      pint = pe1(ks+1)

#ifdef NO_UKMO_HB
      do k=1,ks+1
         ak(k) = pe1(k)
         bk(k) = 0.
      enddo

      do k=ks+2,km+1
         bk(k) = (pe1(k) - pint) / (pe1(km+1)-pint)  ! bk == sigma
         ak(k) =  pe1(k) - bk(k) * pe1(km+1)
      enddo
      bk(km+1) = 1.
      ak(km+1) = 0.
#else
! Problematic for non-hydrostatic
      do k=1,km+1
         eta(k) = pe1(k) / pe1(km+1)
      enddo

      ep =  eta(ks+1)
      es =  eta(km)
!     es =  1.
      alpha = (ep**2-2.*ep*es) / (es-ep)**2
      beta  = 2.*ep*es**2 / (es-ep)**2
      gama = -(ep*es)**2 / (es-ep)**2

! Pure pressure:
      do k=1,ks+1
         ak(k) = eta(k)*1.e5
         bk(k) = 0.
      enddo

      do k=ks+2, km
         ak(k) = alpha*eta(k) + beta + gama/eta(k)
         ak(k) = ak(k)*1.e5
      enddo
         ak(km+1) = 0.

      do k=ks+2, km
         bk(k) = (pe1(k) - ak(k))/pe1(km+1)
      enddo
         bk(km+1) = 1.
#endif

      if ( is_master() ) then
          write(*,*) 'KS=', ks, 'PINT (mb)=', pint/100.
          do k=1,km
             write(*,*) k, 0.5*(pe1(k)+pe1(k+1))/100., dz(k)
          enddo
          tmp1 = ak(ks+1)
          do k=ks+1,km
             tmp1 = max(tmp1, (ak(k)-ak(k+1))/max(1.E-5, (bk(k+1)-bk(k))) )
          enddo
          write(*,*) 'Hybrid Sigma-P: minimum allowable surface pressure (hpa)=', tmp1/100.
      endif


 end subroutine var_hi2


 subroutine var_dz(km, ak, bk, ptop, ks, pint, s_rate)
  integer, intent(in):: km
  real,    intent(in):: ptop
  real,    intent(in):: s_rate        !< between [1. 1.1]
  real,    intent(out):: ak(km+1), bk(km+1)
  real,    intent(inout):: pint
  integer, intent(out):: ks
! Local
  real, parameter:: p00 = 1.E5
  real, dimension(km+1):: ze, pe1, peln, eta
  real, dimension(km):: dz, s_fac, dlnp
  real ztop, t0, dz0, sum1, tmp1
  real ep, es, alpha, beta, gama
  integer  k

     pe1(1) = ptop
     peln(1) = log(pe1(1))
     pe1(km+1) = p00
     peln(km+1) = log(pe1(km+1))

     t0 = 270.
     ztop = rdgas/grav*t0*(peln(km+1) - peln(1))

      s_fac(km  ) = 0.10
      s_fac(km-1) = 0.20
      s_fac(km-2) = 0.30
      s_fac(km-3) = 0.40
      s_fac(km-4) = 0.50
      s_fac(km-5) = 0.60
      s_fac(km-6) = 0.70
      s_fac(km-7) = 0.80
      s_fac(km-8) = 0.90
      s_fac(km-9) = 0.95
      s_fac(km-10) = 0.5*(s_fac(km-9) + s_rate)

      do k=km-11, 9, -1
         s_fac(k) = min(10.0, s_rate * s_fac(k+1) )
      enddo

      s_fac(8) = 0.5*(1.1+s_rate)*s_fac(9)
      s_fac(7) = 1.1 *s_fac(8)
      s_fac(6) = 1.15*s_fac(7)
      s_fac(5) = 1.2 *s_fac(6)
      s_fac(4) = 1.3 *s_fac(5)
      s_fac(3) = 1.4 *s_fac(4)
      s_fac(2) = 1.5 *s_fac(3)
      s_fac(1) = 1.6 *s_fac(2)

      sum1 = 0.
      do k=1,km
         sum1 = sum1 + s_fac(k)
      enddo

      dz0 = ztop / sum1

      do k=1,km
         dz(k) = s_fac(k) * dz0
      enddo

      ze(km+1) = 0.
      do k=km,1,-1
         ze(k) = ze(k+1) + dz(k)
      enddo

! Re-scale dz with the stretched ztop
      do k=1,km
         dz(k) = dz(k) * (ztop/ze(1))
      enddo

      do k=km,1,-1
         ze(k) = ze(k+1) + dz(k)
      enddo
!     ze(1) = ztop

      if ( is_master() ) write(*,*) 'var_dz: computed model top (m)=', ztop*0.001, ' bottom/top dz=', dz(km), dz(1)
      call sm1_edge(1, 1, 1, 1, km, 1, 1, ze, 1)

! Given z --> p
      do k=1,km
          dz(k) = ze(k) - ze(k+1)
        dlnp(k) = grav*dz(k) / (rdgas*t0)
      enddo
      do k=2,km
         peln(k) = peln(k-1) + dlnp(k-1)
          pe1(k) = exp(peln(k))
      enddo

! Pe(k) = ak(k) + bk(k) * PS
! Locate pint and KS
      ks = 0
      do k=2,km
         if ( pint < pe1(k)) then
              ks = k-1
              exit
         endif
      enddo
      if ( is_master() ) then
         write(*,*) 'For (input) PINT=', 0.01*pint, ' KS=', ks, 'pint(computed)=', 0.01*pe1(ks+1)
         write(*,*) 'ptop =', ptop
      endif
      pint = pe1(ks+1)

#ifdef NO_UKMO_HB
      do k=1,ks+1
         ak(k) = pe1(k)
         bk(k) = 0.
      enddo

      do k=ks+2,km+1
         bk(k) = (pe1(k) - pint) / (pe1(km+1)-pint)  ! bk == sigma
         ak(k) =  pe1(k) - bk(k) * pe1(km+1)
      enddo
      bk(km+1) = 1.
      ak(km+1) = 0.
#else
! Problematic for non-hydrostatic
      do k=1,km+1
         eta(k) = pe1(k) / pe1(km+1)
      enddo

      ep =  eta(ks+1)
      es =  eta(km)
!     es =  1.
      alpha = (ep**2-2.*ep*es) / (es-ep)**2
      beta  = 2.*ep*es**2 / (es-ep)**2
      gama = -(ep*es)**2 / (es-ep)**2

! Pure pressure:
      do k=1,ks+1
         ak(k) = eta(k)*1.e5
         bk(k) = 0.
      enddo

      do k=ks+2, km
         ak(k) = alpha*eta(k) + beta + gama/eta(k)
         ak(k) = ak(k)*1.e5
      enddo
         ak(km+1) = 0.

      do k=ks+2, km
         bk(k) = (pe1(k) - ak(k))/pe1(km+1)
      enddo
         bk(km+1) = 1.
#endif

      if ( is_master() ) then
          write(*,*) 'KS=', ks, 'PINT (mb)=', pint/100.
          do k=1,km
             write(*,*) k, 0.5*(pe1(k)+pe1(k+1))/100., dz(k)
          enddo
          tmp1 = ak(ks+1)
          do k=ks+1,km
             tmp1 = max(tmp1, (ak(k)-ak(k+1))/max(1.E-5, (bk(k+1)-bk(k))) )
          enddo
          write(*,*) 'Hybrid Sigma-P: minimum allowable surface pressure (hpa)=', tmp1/100.
      endif


 end subroutine var_dz


 subroutine var55_dz(km, ak, bk, ptop, ks, pint, s_rate)
  integer, intent(in):: km
  real,    intent(in):: ptop
  real,    intent(in):: s_rate        !< between [1. 1.1]
  real,    intent(out):: ak(km+1), bk(km+1)
  real,    intent(inout):: pint
  integer, intent(out):: ks
! Local
  real, parameter:: p00 = 1.E5
  real, dimension(km+1):: ze, pe1, peln, eta
  real, dimension(km):: dz, s_fac, dlnp
  real ztop, t0, dz0, sum1, tmp1
  real ep, es, alpha, beta, gama
  integer  k

     pe1(1) = ptop
     peln(1) = log(pe1(1))
     pe1(km+1) = p00
     peln(km+1) = log(pe1(km+1))

     t0 = 270.
     ztop = rdgas/grav*t0*(peln(km+1) - peln(1))

      s_fac(km  ) = 0.10
      s_fac(km-1) = 0.15
      s_fac(km-2) = 0.20
      s_fac(km-3) = 0.25
      s_fac(km-4) = 0.30
      s_fac(km-5) = 0.35
      s_fac(km-6) = 0.40
      s_fac(km-7) = 0.45
      s_fac(km-8) = 0.50
      s_fac(km-9) = 0.55
      s_fac(km-10) = 0.60
      s_fac(km-11) = 0.70
      s_fac(km-12) = 0.85
      s_fac(km-13) = 1.00

      do k=km-14, 9, -1
         s_fac(k) = min(10.0, s_rate * s_fac(k+1) )
      enddo

      s_fac(8) = 0.5*(1.1+s_rate)*s_fac(9)
      s_fac(7) = 1.1 *s_fac(8)
      s_fac(6) = 1.15*s_fac(7)
      s_fac(5) = 1.2 *s_fac(6)
      s_fac(4) = 1.3 *s_fac(5)
      s_fac(3) = 1.4 *s_fac(4)
      s_fac(2) = 1.5 *s_fac(3)
      s_fac(1) = 1.6 *s_fac(2)

      sum1 = 0.
      do k=1,km
         sum1 = sum1 + s_fac(k)
      enddo

      dz0 = ztop / sum1

      do k=1,km
         dz(k) = s_fac(k) * dz0
      enddo

      ze(km+1) = 0.
      do k=km,1,-1
         ze(k) = ze(k+1) + dz(k)
      enddo

! Re-scale dz with the stretched ztop
      do k=1,km
         dz(k) = dz(k) * (ztop/ze(1))
      enddo

      do k=km,1,-1
         ze(k) = ze(k+1) + dz(k)
      enddo
!     ze(1) = ztop

      if ( is_master() ) write(*,*) 'var55_dz: computed model top (m)=', ztop*0.001, ' bottom/top dz=', dz(km), dz(1)
      call sm1_edge(1, 1, 1, 1, km, 1, 1, ze, 2)

! Given z --> p
      do k=1,km
          dz(k) = ze(k) - ze(k+1)
        dlnp(k) = grav*dz(k) / (rdgas*t0)
      enddo
      do k=2,km
         peln(k) = peln(k-1) + dlnp(k-1)
          pe1(k) = exp(peln(k))
      enddo

! Pe(k) = ak(k) + bk(k) * PS
! Locate pint and KS
      ks = 0
      do k=2,km
         if ( pint < pe1(k)) then
              ks = k-1
              exit
         endif
      enddo
      if ( is_master() ) then
         write(*,*) 'For (input) PINT=', 0.01*pint, ' KS=', ks, 'pint(computed)=', 0.01*pe1(ks+1)
      endif
      pint = pe1(ks+1)

#ifdef NO_UKMO_HB
      do k=1,ks+1
         ak(k) = pe1(k)
         bk(k) = 0.
      enddo

      do k=ks+2,km+1
         bk(k) = (pe1(k) - pint) / (pe1(km+1)-pint)  ! bk == sigma
         ak(k) =  pe1(k) - bk(k) * pe1(km+1)
      enddo
      bk(km+1) = 1.
      ak(km+1) = 0.
#else
! Problematic for non-hydrostatic
      do k=1,km+1
         eta(k) = pe1(k) / pe1(km+1)
      enddo

      ep =  eta(ks+1)
      es =  eta(km)
!     es =  1.
      alpha = (ep**2-2.*ep*es) / (es-ep)**2
      beta  = 2.*ep*es**2 / (es-ep)**2
      gama = -(ep*es)**2 / (es-ep)**2

! Pure pressure:
      do k=1,ks+1
         ak(k) = eta(k)*1.e5
         bk(k) = 0.
      enddo

      do k=ks+2, km
         ak(k) = alpha*eta(k) + beta + gama/eta(k)
         ak(k) = ak(k)*1.e5
      enddo
         ak(km+1) = 0.

      do k=ks+2, km
         bk(k) = (pe1(k) - ak(k))/pe1(km+1)
      enddo
         bk(km+1) = 1.
#endif

      if ( is_master() ) then
          write(*,*) 'KS=', ks, 'PINT (mb)=', pint/100.
          do k=1,km
             write(*,*) k, 0.5*(pe1(k)+pe1(k+1))/100., dz(k)
          enddo
          tmp1 = ak(ks+1)
          do k=ks+1,km
             tmp1 = max(tmp1, (ak(k)-ak(k+1))/max(1.E-5, (bk(k+1)-bk(k))) )
          enddo
          write(*,*) 'Hybrid Sigma-P: minimum allowable surface pressure (hpa)=', tmp1/100.
      endif


 end subroutine var55_dz

 subroutine hybrid_z_dz(km, dz, ztop, s_rate)
  integer, intent(in):: km
  real,    intent(in):: s_rate        !< between [1. 1.1]
  real,    intent(in):: ztop
  real,    intent(out):: dz(km)
! Local
  real, parameter:: p00 = 1.E5
  real, dimension(km+1):: ze, pe1, peln, eta
  real, dimension(km):: s_fac, dlnp
  real t0, dz0, sum1, tmp1
  real ep, es, alpha, beta, gama
  integer  k

       s_fac(km  ) = 0.12
       s_fac(km-1) = 0.20
       s_fac(km-2) = 0.30
       s_fac(km-3) = 0.40
       s_fac(km-4) = 0.50
       s_fac(km-5) = 0.60
       s_fac(km-6) = 0.70
       s_fac(km-7) = 0.80
       s_fac(km-8) = 0.90
       s_fac(km-9) = 1.

       do k=km-10, 9, -1
          s_fac(k) = min(4.0, s_rate * s_fac(k+1) )
       enddo

       s_fac(8) = 1.05*s_fac(9)
       s_fac(7) = 1.1 *s_fac(8)
       s_fac(6) = 1.15*s_fac(7)
       s_fac(5) = 1.2 *s_fac(6)
       s_fac(4) = 1.3 *s_fac(5)
       s_fac(3) = 1.4 *s_fac(4)
       s_fac(2) = 1.5 *s_fac(3)
       s_fac(1) = 1.6 *s_fac(2)

       sum1 = 0.
       do k=1,km
          sum1 = sum1 + s_fac(k)
       enddo

       dz0 = ztop / sum1

       do k=1,km
          dz(k) = s_fac(k) * dz0
       enddo

       ze(km+1) = 0.
       do k=km,1,-1
          ze(k) = ze(k+1) + dz(k)
       enddo

       ze(1) = ztop

       call sm1_edge(1, 1, 1, 1, km, 1, 1, ze, 2)

       do k=1,km
            dz(k) = ze(k) - ze(k+1)
       enddo

 end subroutine hybrid_z_dz

!>@brief The subroutine 'get_eta_level' returns the interface and
!! layer-mean pressures for reference.
 subroutine get_eta_level(npz, p_s, pf, ph, ak, bk, pscale)
  integer, intent(in) :: npz
  real, intent(in)  :: p_s            !< unit: pascal
  real, intent(in)  :: ak(npz+1)
  real, intent(in)  :: bk(npz+1)
  real, intent(in), optional :: pscale
  real, intent(out) :: pf(npz)
  real, intent(out) :: ph(npz+1)
  integer k

  ph(1) = ak(1)
  do k=2,npz+1
     ph(k) = ak(k) + bk(k)*p_s
  enddo

  if ( present(pscale) ) then
      do k=1,npz+1
         ph(k) = pscale*ph(k)
      enddo
  endif

  if( ak(1) > 1.E-8 ) then
     pf(1) = (ph(2) - ph(1)) / log(ph(2)/ph(1))
  else
     pf(1) = (ph(2) - ph(1)) * kappa/(kappa+1.)
  endif

  do k=2,npz
     pf(k) = (ph(k+1) - ph(k)) / log(ph(k+1)/ph(k))
  enddo

 end subroutine get_eta_level



 subroutine compute_dz(km, ztop, dz)

  integer, intent(in):: km
  real,   intent(in):: ztop        ! try 50.E3
  real,   intent(out):: dz(km)
!------------------------------
  real ze(km+1), dzt(km)
  integer k


! ztop = 30.E3
  dz(1) = ztop / real(km)
  dz(km) = 0.5*dz(1)

  do k=2,km-1
     dz(k) = dz(1)
  enddo

! Top:
  dz(1) = 2.*dz(2)

  ze(km+1) = 0.
  do k=km,1,-1
     ze(k) = ze(k+1) + dz(k)
  enddo

  if ( is_master() ) then
       write(*,*) 'Hybrid_z:  dz, zm'
       do k=1,km
          dzt(k) = 0.5*(ze(k)+ze(k+1)) / 1000.
          write(*,*) k, dz(k), dzt(k)
       enddo
  endif

 end subroutine compute_dz

 subroutine compute_dz_var(km, ztop, dz)

  integer, intent(in):: km
  real,   intent(in):: ztop        ! try 50.E3
  real,   intent(out):: dz(km)
!------------------------------
  real, parameter:: s_rate = 1.0
  real ze(km+1)
  real s_fac(km)
  real sum1, dz0
  integer k

      s_fac(km  ) = 0.125
      s_fac(km-1) = 0.20
      s_fac(km-2) = 0.30
      s_fac(km-3) = 0.40
      s_fac(km-4) = 0.50
      s_fac(km-5) = 0.60
      s_fac(km-6) = 0.70
      s_fac(km-7) = 0.80
      s_fac(km-8) = 0.90
      s_fac(km-9) = 1.

      do k=km-10, 9, -1
         s_fac(k) = s_rate * s_fac(k+1)
      enddo

      s_fac(8) = 1.05*s_fac(9)
      s_fac(7) = 1.1 *s_fac(8)
      s_fac(6) = 1.15*s_fac(7)
      s_fac(5) = 1.2 *s_fac(6)
      s_fac(4) = 1.3 *s_fac(5)
      s_fac(3) = 1.4 *s_fac(4)
      s_fac(2) = 1.5 *s_fac(3)
      s_fac(1) = 1.6 *s_fac(2)

      sum1 = 0.
      do k=1,km
         sum1 = sum1 + s_fac(k)
      enddo

      dz0 = ztop / sum1

      do k=1,km
         dz(k) = s_fac(k) * dz0
      enddo

      ze(1) = ztop
      ze(km+1) = 0.
      do k=km,2,-1
         ze(k) = ze(k+1) + dz(k)
      enddo

! Re-scale dz with the stretched ztop
      do k=1,km
         dz(k) = dz(k) * (ztop/ze(1))
      enddo

      do k=km,2,-1
         ze(k) = ze(k+1) + dz(k)
      enddo

      call sm1_edge(1, 1, 1, 1, km, 1, 1, ze, 2)

      do k=1,km
         dz(k) = ze(k) - ze(k+1)
      enddo

 end subroutine compute_dz_var

 subroutine compute_dz_L32(km, ztop, dz)

  integer, intent(in):: km
  real,   intent(out):: dz(km)
  real,   intent(out):: ztop        ! try 50.E3
!------------------------------
  real dzt(km)
  real ze(km+1)
  real dz0, dz1, dz2
  real z0, z1, z2
  integer k, k0, k1, k2, n

!-------------------
        k2 =  8
        z2 = 30.E3
!-------------------
        k1 = 21
        z1 = 10.0E3
!-------------------
        k0 = 2
        z0 = 0.
        dz0 = 75.   ! meters
!-------------------
! Treat the surface layer as a special layer
        ze(1) = z0
        dz(1) = dz0

        ze(2) = dz(1)
          dz0 = 1.5*dz0
        dz(2) = dz0

        ze(3) = ze(2) + dz(2)

        dz1 = 2.*(z1-ze(3) - k1*dz0) / (k1*(k1-1))

        do k=k0+1,k0+k1
           dz(k) = dz0 + (k-k0)*dz1
           ze(k+1) = ze(k) + dz(k)
        enddo

        dz0 = dz(k1+k0)
        dz2 = 2.*(z2-ze(k0+k1+1)-k2*dz0) / (k2*(k2-1))

        do k=k0+k1+1,k0+k1+k2
           dz(k) = dz0 + (k-k0-k1)*dz2
           ze(k+1) = ze(k) + dz(k)
        enddo

        dz(km) = 2.*dz(km-1)
        ztop = ze(km) + dz(km)
        ze(km+1) = ze(km) + dz(km)

        call zflip (dz, 1, km)

        ze(km+1) = 0.
        do k=km,1,-1
           ze(k) = ze(k+1) + dz(k)
        enddo

!       if ( is_master() ) then
!          write(*,*) 'Hybrid_z:  dz, zm'
!          do k=1,km
!             dzt(k) = 0.5*(ze(k)+ze(k+1)) / 1000.
!             write(*,*) k, dz(k), dzt(k)
!          enddo
!       endif

 end subroutine compute_dz_L32

 subroutine compute_dz_L101(km, ztop, dz)

  integer, intent(in):: km        ! km==101
  real,   intent(out):: dz(km)
  real,   intent(out):: ztop      ! try 30.E3
!------------------------------
  real ze(km+1)
  real dz0, dz1
  real:: stretch_f = 1.16
  integer k, k0, k1

         k1 = 2
         k0 = 25
        dz0 = 40.   ! meters

        ze(km+1) = 0.

        do k=km, k0, -1
           dz(k) = dz0
           ze(k) = ze(k+1) + dz(k)
        enddo

        do k=k0+1, k1, -1
           dz(k) = stretch_f * dz(k+1)
           ze(k) = ze(k+1) + dz(k)
        enddo

        dz(1) = 4.0*dz(2)
        ze(1) = ze(2) + dz(1)
        ztop = ze(1)

        if ( is_master() ) then
           write(*,*) 'Hybrid_z:  dz, ze'
           do k=1,km
              write(*,*) k, 0.001*dz(k), 0.001*ze(k)
           enddo
!  ztop (km) = 20.2859154
           write(*,*) 'ztop (km) =', ztop * 0.001
        endif

 end subroutine compute_dz_L101

 subroutine set_hybrid_z(is, ie, js, je, ng, km, ztop, dz, rgrav, hs, ze, dz3)

 integer,  intent(in):: is, ie, js, je, ng, km
 real, intent(in):: rgrav, ztop
 real, intent(in):: dz(km)       !< Reference vertical resolution for zs=0
 real, intent(in):: hs(is-ng:ie+ng,js-ng:je+ng)
 real, intent(inout)::  ze(is:ie,js:je,km+1)
 real, optional, intent(out):: dz3(is-ng:ie+ng,js-ng:je+ng,km)
! local
 logical:: filter_xy = .false.
 real, allocatable:: delz(:,:,:)
 integer ntimes
 real zint
 real:: z1(is:ie,js:je)
 real:: z(km+1)
 real sig, z_rat
 integer ks(is:ie,js:je)
 integer i, j, k, ks_min, kint

 z(km+1) = 0.
 do k=km,1,-1
    z(k) = z(k+1) + dz(k)
 enddo

  do j=js,je
     do i=is,ie
        ze(i,j,   1) = ztop
        ze(i,j,km+1) = hs(i,j) * rgrav
     enddo
  enddo

 do k=2,km
   do j=js,je
      do i=is,ie
         ze(i,j,k) = z(k)
      enddo
   enddo
 enddo

! Set interface:
#ifndef USE_VAR_ZINT
  zint = 12.0E3
  ntimes = 2
  kint = 2
  do k=2,km
     if ( z(k)<=zint ) then
          kint = k
          exit
     endif
  enddo

  if ( is_master() ) write(*,*) 'Z_coord interface set at k=',kint, ' ZE=', z(kint)

  do j=js,je
     do i=is,ie
        z_rat = (ze(i,j,kint)-ze(i,j,km+1)) / (z(kint)-z(km+1))
        do k=km,kint+1,-1
           ze(i,j,k) = ze(i,j,k+1) + dz(k)*z_rat
        enddo
!--------------------------------------
! Apply vertical smoother locally to dz
!--------------------------------------
        call sm1_edge(is, ie, js, je, km, i, j, ze, ntimes)
    enddo
  enddo
#else
! ZINT is a function of local terrain
  ntimes = 4
  do j=js,je
     do i=is,ie
        z1(i,j) = dim(ze(i,j,km+1), 2500.) + 7500.
     enddo
  enddo

   ks_min = km
   do j=js,je
      do i=is,ie
         do k=km,2,-1
            if ( z(k)>=z1(i,j) ) then
                 ks(i,j) = k
                 ks_min = min(ks_min, k)
                 go to 555
            endif
        enddo
555     continue
      enddo
   enddo

  do j=js,je
     do i=is,ie
        kint = ks(i,j) + 1
        z_rat = (ze(i,j,kint)-ze(i,j,km+1)) / (z(kint)-z(km+1))
        do k=km,kint+1,-1
           ze(i,j,k) = ze(i,j,k+1) + dz(k)*z_rat
        enddo
!--------------------------------------
! Apply vertical smoother locally to dz
!--------------------------------------
        call sm1_edge(is, ie, js, je, km, i, j, ze, ntimes)
    enddo
  enddo
#endif

#ifdef DEV_ETA
    if ( filter_xy ) then
        allocate (delz(isd:ied, jsd:jed, km) )
        ntimes = 2
        do k=1,km
           do j=js,je
              do i=is,ie
                 delz(i,j,k) = ze(i,j,k+1) - ze(i,j,k)
              enddo
           enddo
        enddo
        call del2_cubed(delz, 0.2*da_min, npx, npy, km, ntimes)
        do k=km,2,-1
           do j=js,je
              do i=is,ie
                 ze(i,j,k) = ze(i,j,k+1) - delz(i,j,k)
              enddo
           enddo
        enddo
        deallocate ( delz )
    endif
#endif
    if ( present(dz3) ) then
        do k=1,km
           do j=js,je
              do i=is,ie
                 dz3(i,j,k) = ze(i,j,k+1) - ze(i,j,k)
              enddo
           enddo
        enddo
    endif

  end subroutine set_hybrid_z


  subroutine sm1_edge(is, ie, js, je, km, i, j, ze, ntimes)
  integer, intent(in):: is, ie, js, je, km
  integer, intent(in):: ntimes, i, j
  real, intent(inout):: ze(is:ie,js:je,km+1)
! local:
  real, parameter:: df = 0.25
  real dz(km)
  real flux(km+1)
  integer k, n, k1, k2

      k2 = km-1
      do k=1,km
         dz(k) = ze(i,j,k+1) - ze(i,j,k)
      enddo

   do n=1,ntimes
      k1 = 2 + (ntimes-n)

      flux(k1  ) = 0.
      flux(k2+1) = 0.
      do k=k1+1,k2
         flux(k) = df*(dz(k) - dz(k-1))
      enddo

      do k=k1,k2
         dz(k) = dz(k) - flux(k) + flux(k+1)
      enddo
   enddo

   do k=km,1,-1
      ze(i,j,k) = ze(i,j,k+1) - dz(k)
   enddo

  end subroutine sm1_edge



  subroutine gw_1d(km, p0, ak, bk, ptop, ztop, pt1)
  integer, intent(in):: km
  real,    intent(in):: p0, ztop
  real,    intent(inout):: ptop
  real,    intent(inout):: ak(km+1), bk(km+1)
  real, intent(out):: pt1(km)
! Local
  logical:: isothermal
  real, dimension(km+1):: ze, pe1, pk1
  real, dimension(km):: dz1
  real t0, n2, s0
  integer  k

! Set up vertical coordinare with constant del-z spacing:
       isothermal = .false.
       t0 = 300.

       if ( isothermal ) then
            n2 = grav**2/(cp_air*t0)
       else
            n2 = 0.0001
       endif

       s0 = grav*grav / (cp_air*n2)

       ze(km+1) = 0.
       do k=km,1,-1
          dz1(k) = ztop / real(km)
           ze(k) = ze(k+1) + dz1(k)
       enddo

! Given z --> p
       do k=1,km+1
          pe1(k) = p0*( (1.-s0/t0) + s0/t0*exp(-n2*ze(k)/grav) )**(1./kappa)
       enddo

       ptop = pe1(1)
!      if ( is_master() ) write(*,*) 'GW_1D: computed model top (pa)=', ptop

! Set up "sigma" coordinate
       ak(1) = pe1(1)
       bk(1) = 0.
       do k=2,km
          bk(k) = (pe1(k) - pe1(1)) / (pe1(km+1)-pe1(1))  ! bk == sigma
          ak(k) =  pe1(1)*(1.-bk(k))
       enddo
       ak(km+1) = 0.
       bk(km+1) = 1.

       do k=1,km+1
          pk1(k) = pe1(k) ** kappa
       enddo

! Compute volume mean potential temperature with hydrostatic eqn:
       do k=1,km
          pt1(k) = grav*dz1(k) / ( cp_air*(pk1(k+1)-pk1(k)) )
       enddo

  end subroutine gw_1d

 subroutine mount_waves(km, ak, bk, ptop, ks, pint)
  integer, intent(in):: km
  real,    intent(out):: ak(km+1), bk(km+1)
  real,    intent(out):: ptop, pint
  integer, intent(out):: ks
! Local
  real, parameter:: p00 = 1.E5
  real, dimension(km+1):: ze, pe1, peln, eta
  real, dimension(km):: dz, dlnp
  real ztop, t0, dz0, sum1, tmp1
  real ep, es, alpha, beta, gama, s_fac
  integer  k, k500

  pint = 300.e2
!      s_fac = 1.05
!      dz0 = 500.
  if ( km <= 60 ) then
       s_fac = 1.0
       dz0 = 500.
  else
       s_fac = 1.
       dz0 = 250.
  endif

! Basic parameters for HIWPP mountain waves
   t0 = 300.
! ztop = 20.0e3; 500-m resolution in halft of the vertical domain
! ztop = real(km-1)*500.
!-----------------------
! Compute temp ptop based on isothermal atm
! ptop = p00*exp(-grav*ztop/(rdgas*t0))

! Lowest half has constant resolution
     ze(km+1) = 0.
     do k=km, km-19, -1
        ze(k) = ze(k+1) + dz0
     enddo

! Stretching from 10-km and up:
     do k=km-20, 3,  -1
        dz0 = s_fac * dz0
        ze(k) = ze(k+1) + dz0
     enddo
     ze(2) = ze(3) + sqrt(2.)*dz0
     ze(1) = ze(2) + 2.0*dz0

!    call sm1_edge(1, 1, 1, 1, km, 1, 1, ze, 1)

! Given z --> p
     do k=1,km
         dz(k) = ze(k) - ze(k+1)
       dlnp(k) = grav*dz(k) / (rdgas*t0)
     enddo

      pe1(km+1) = p00
     peln(km+1) = log(p00)
     do k=km,1,-1
        peln(k) = peln(k+1) - dlnp(k)
         pe1(k) = exp(peln(k))
     enddo

! Comnpute new ptop
     ptop = pe1(1)

! Pe(k) = ak(k) + bk(k) * PS
! Locate pint and KS
      ks = 0
      do k=2,km
         if ( pint < pe1(k)) then
              ks = k-1
              exit
         endif
      enddo

      if ( is_master() ) then
         write(*,*) 'For (input) PINT=', 0.01*pint, ' KS=', ks, 'pint(computed)=', 0.01*pe1(ks+1)
         write(*,*) 'Modified ptop =', ptop, ' ztop=', ze(1)/1000.
         do k=1,km
            write(*,*) k, 'ze =', ze(k)/1000.
         enddo
      endif
      pint = pe1(ks+1)

#ifdef NO_UKMO_HB
      do k=1,ks+1
         ak(k) = pe1(k)
         bk(k) = 0.
      enddo

      do k=ks+2,km+1
         bk(k) = (pe1(k) - pint) / (pe1(km+1)-pint)  ! bk == sigma
         ak(k) =  pe1(k) - bk(k) * pe1(km+1)
      enddo
      bk(km+1) = 1.
      ak(km+1) = 0.
#else
! Problematic for non-hydrostatic
      do k=1,km+1
         eta(k) = pe1(k) / pe1(km+1)
      enddo
      ep =  eta(ks+1)
      es =  eta(km)
!     es =  1.
      alpha = (ep**2-2.*ep*es) / (es-ep)**2
      beta  = 2.*ep*es**2 / (es-ep)**2
      gama = -(ep*es)**2 / (es-ep)**2

! Pure pressure:
      do k=1,ks+1
         ak(k) = eta(k)*1.e5
         bk(k) = 0.
      enddo

      do k=ks+2, km
         ak(k) = alpha*eta(k) + beta + gama/eta(k)
         ak(k) = ak(k)*1.e5
      enddo
         ak(km+1) = 0.

      do k=ks+2, km
         bk(k) = (pe1(k) - ak(k))/pe1(km+1)
      enddo
         bk(km+1) = 1.
#endif

      if ( is_master() ) then
          tmp1 = ak(ks+1)
          do k=ks+1,km
             tmp1 = max(tmp1, (ak(k)-ak(k+1))/max(1.E-5, (bk(k+1)-bk(k))) )
          enddo
          write(*,*) 'Hybrid Sigma-P: minimum allowable surface pressure (hpa)=', tmp1/100.
      endif

 end subroutine mount_waves


  subroutine zflip(q, im, km)
  integer, intent(in):: im, km
  real, intent(inout):: q(im,km)
!---
  integer i, k
  real qtmp

    do i = 1, im
       do k = 1, (km+1)/2
          qtmp = q(i,k)
          q(i,k) = q(i,km+1-k)
          q(i,km+1-k) = qtmp
       end do
    end do

  end subroutine zflip

end module fv_eta_mod
