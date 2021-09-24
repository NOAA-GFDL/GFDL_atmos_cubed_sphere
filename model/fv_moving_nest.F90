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


!*************************************************************************
!>@brief!   Provides Moving Nest functionality in FV3 dynamic core.  
!!>@author W. Ramstrom, AOML/HRD  01/15/2021
!
! 
! =======================================================================!
!

! =======================================================================!
!
! Notes
!  
!------------------------------------------------------------------------
! Moving Nest Subroutine Naming Convention
!----------------------------------------------------------------------- 
!
! mn_meta_* subroutines perform moving nest operations for FV3 metadata.
!               These routines will run only once per nest move.
!
! mn_var_*  subroutines perform moving nest operations for an individual FV3 variable.
!               These routines will run many times per nest move.
!
! mn_prog_* subroutines perform moving nest operations for the list of prognostic fields.
!               These routines will run only once per nest move.
!
! mn_phys_* subroutines perform moving nest operations for the list of physics fields.
!               These routines will run only once per nest move.
!
! =======================================================================!

#define REMAP 1 

module fv_moving_nest_mod
#ifdef MOVING_NEST
  
  use block_control_mod,      only : block_control_type
  use fms_mod,                only : mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE, clock_flag_default
  use mpp_mod,                only : mpp_pe, mpp_sync, mpp_sync_self, mpp_send, mpp_error, NOTE, FATAL
  use mpp_domains_mod,        only : mpp_update_domains, mpp_get_data_domain, mpp_get_global_domain
  use mpp_domains_mod,        only : mpp_define_nest_domains, mpp_shift_nest_domains, nest_domain_type, domain2d
  use mpp_domains_mod,        only : mpp_get_C2F_index, mpp_update_nest_fine
  use mpp_domains_mod,        only : mpp_get_F2C_index, mpp_update_nest_coarse
  use mpp_domains_mod,        only : NORTH, SOUTH, EAST, WEST, CORNER, CENTER
  use mpp_domains_mod,        only : NUPDATE, SUPDATE, EUPDATE, WUPDATE, DGRID_NE

#ifdef GFS_TYPES
use GFS_typedefs,           only: IPD_data_type => GFS_data_type, &
                                  IPD_control_type => GFS_control_type, kind_phys
#else
use IPD_typedefs,           only: IPD_data_type, IPD_control_type, kind_phys => IPD_kind_phys
#endif
  use GFS_init,               only: GFS_grid_populate


  use boundary_mod,           only: update_coarse_grid, update_coarse_grid_mpp
  use bounding_box_mod,       only: bbox, bbox_get_C2F_index, fill_bbox, show_bbox
  use constants_mod,          only: cp_air, rdgas, grav, rvgas, kappa, pstd_mks, hlv
  use field_manager_mod,      only: MODEL_ATMOS
  use fms_io_mod,             only: read_data, write_data, get_global_att_value, fms_io_init, fms_io_exit
  use fv_arrays_mod,          only: fv_atmos_type, fv_nest_type, fv_grid_type, R_GRID
  use fv_arrays_mod,          only: fv_moving_nest_prog_type, fv_moving_nest_physics_type
  use fv_arrays_mod,          only: allocate_fv_nest_bc_type, deallocate_fv_nest_bc_type
  use fv_grid_tools_mod,      only: init_grid
  use fv_grid_utils_mod,      only: grid_utils_init, ptop_min, dist2side_latlon
  use fv_mapz_mod,            only: Lagrangian_to_Eulerian, moist_cv, compute_total_energy
  use fv_moving_nest_logging_mod, only: check_array, check_local_array, show_atm, show_atm_grids, show_nest_grid, show_tile_geo, grid_equal
  use fv_nesting_mod,         only: dealloc_nested_buffers
  use fv_nwp_nudge_mod,       only: do_adiabatic_init
  use init_hydro_mod,         only: p_var
  use tracer_manager_mod,     only: get_tracer_index, get_tracer_names
  use fv_moving_nest_utils_mod,  only: alloc_halo_buffer, load_nest_latlons_from_nc, grid_geometry, output_grid_to_nc, find_nest_alignment
  use fv_moving_nest_utils_mod,  only: fill_nest_from_buffer, fill_nest_from_buffer_cell_center, fill_nest_from_buffer_nearest_neighbor
  use fv_moving_nest_utils_mod,  only: fill_nest_halos_from_parent, fill_grid_from_supergrid, fill_weight_grid
  use fv_moving_nest_utils_mod,  only: alloc_read_data

  implicit none

#ifdef NO_QUAD_PRECISION
  ! 64-bit precision (kind=8)
  integer, parameter:: f_p = selected_real_kind(15)
#else
  ! Higher precision (kind=16) for grid geometrical factors:
  integer, parameter:: f_p = selected_real_kind(20)
#endif

#ifdef OVERLOAD_R4
  real, parameter:: real_snan=x'FFBFFFFF'
#else
  real, parameter:: real_snan=x'FFF7FFFFFFFFFFFF'
#endif

  logical :: debug_log = .false.
  logical :: move_physics = .true.
  logical :: move_nsst = .true.

#include <fms_platform.h>

  !! Step 2
  interface mn_var_fill_intern_nest_halos
     module procedure mn_var_fill_intern_nest_halos2D
     module procedure mn_var_fill_intern_nest_halos2D_kindphys
     module procedure mn_var_fill_intern_nest_halos3D
     module procedure mn_var_fill_intern_nest_halos3D_kindphys
     module procedure mn_var_fill_intern_nest_halos4D
     module procedure mn_var_fill_intern_nest_halos4D_kindphys
     module procedure mn_var_fill_intern_nest_halos_wind
  end interface mn_var_fill_intern_nest_halos


  !! Step 6
  interface mn_var_shift_data
     module procedure mn_var_shift_data2D
     module procedure mn_var_shift_data2D_kindphys
     module procedure mn_var_shift_data3D
     module procedure mn_var_shift_data3D_kindphys
     module procedure mn_var_shift_data4D
     module procedure mn_var_shift_data4D_kindphys
  end interface mn_var_shift_data

  !! Step 8
  interface mn_var_dump_to_netcdf
     module procedure mn_var_dump_2d_to_netcdf
     module procedure mn_var_dump_3d_to_netcdf
  end interface mn_var_dump_to_netcdf

contains

  !!===================================================================================== 
  !! Step 1.9 -- Allocate and fill the temporary variable(s)
  !!            This is to manage variables that are not allocated with a halo
  !!            on the Atm structure
  !!=====================================================================================           


  subroutine mn_prog_fill_temp_variables(Atm, n, child_grid_num, is_fine_pe, npz)
    type(fv_atmos_type), allocatable, target, intent(in)     :: Atm(:)
    integer, intent(in)                              :: n, child_grid_num   !  n is the nest level
    logical, intent(in)                              :: is_fine_pe
    integer, intent(in)                              :: npz

    integer :: isd, ied, jsd, jed
    integer :: is, ie, js, je
    integer :: this_pe
    type(fv_moving_nest_prog_type), pointer :: mn_prog

    mn_prog => Atm(n)%mn_prog
    
    this_pe = mpp_pe()

    if (debug_log) print '("[INFO] WDR start mn_prog_fill_temp_variables. npe=",I0," n=",I0)', this_pe, n

    isd = Atm(n)%bd%isd
    ied = Atm(n)%bd%ied
    jsd = Atm(n)%bd%jsd
    jed = Atm(n)%bd%jed
    
    if (debug_log) print '("[INFO] WDR mn_prog_fill_temp_variables. npe=",I0," isd=",I0," ied=",I0," jsd=",I0," jed=",I0)', this_pe, isd, ied, jsd, jed

    is = Atm(n)%bd%is
    ie = Atm(n)%bd%ie
    js = Atm(n)%bd%js
    je = Atm(n)%bd%je
    
    if (debug_log) print '("[INFO] WDR mn_prog_fill_temp_variables. npe=",I0," is=",I0," ie=",I0," js=",I0," je=",I0)', this_pe, is, ie, js, je
    
    ! Reset this to a dummy value, to help flag if the halos don't get updated later.
    mn_prog%delz = +99999.9
    
    mn_prog%delz(is:ie, js:je, 1:npz) =  Atm(n)%delz(is:ie, js:je, 1:npz) 

    if (debug_log) print '("[INFO] WDR Z mn_prog_fill_temp_variables. npe=",I0," npz=",I0," ",I0," ",I0)', this_pe, npz, lbound(Atm(n)%delz,3), ubound(Atm(n)%delz,3)

    if (debug_log) print '("[INFO] WDR end mn_prog_fill_temp_variables. npe=",I0," n=",I0)', this_pe, n

    
  end subroutine mn_prog_fill_temp_variables


  subroutine mn_phys_fill_temp_variables(Atm, Atm_block, IPD_Control, IPD_Data, n, child_grid_num, is_fine_pe, npz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)
    type (block_control_type), intent(in)            :: Atm_block
    type(IPD_control_type), intent(in)               :: IPD_Control
    type(IPD_data_type), intent(inout)               :: IPD_Data(:)
    integer, intent(in)                              :: n, child_grid_num
    logical, intent(in)                              :: is_fine_pe
    integer, intent(in)                              :: npz

    integer :: isd, ied, jsd, jed
    integer :: is, ie, js, je
    integer :: this_pe

    integer :: nb, blen, i, j, k, ix, nv
    type(fv_moving_nest_physics_type), pointer       :: mn_phys

    this_pe = mpp_pe()


    if (debug_log) print '("[INFO] WDR start mn_phys_fill_temp_variables. npe=",I0," n=",I0)', this_pe, n

    isd = Atm(n)%bd%isd
    ied = Atm(n)%bd%ied
    jsd = Atm(n)%bd%jsd
    jed = Atm(n)%bd%jed
    
    if (debug_log) print '("[INFO] WDR mn_phys_fill_temp_variables. npe=",I0," isd=",I0," ied=",I0," jsd=",I0," jed=",I0)', this_pe, isd, ied, jsd, jed

    is = Atm(n)%bd%is
    ie = Atm(n)%bd%ie
    js = Atm(n)%bd%js
    je = Atm(n)%bd%je
    
    if (debug_log) print '("[INFO] WDR mn_phys_fill_temp_variables. npe=",I0," is=",I0," ie=",I0," js=",I0," je=",I0)', this_pe, is, ie, js, je
    
    mn_phys => Atm(n)%mn_phys

    mn_phys%ts(is:ie, js:je) =  Atm(n)%ts(is:ie, js:je) 

    do nb = 1,Atm_block%nblks
      blen = Atm_block%blksz(nb)
      do ix = 1, blen
         ! Get the indices only once, before iterating through vertical levels or number of variables
         !  Was there a different efficiency from having the k loop outside?
         i = Atm_block%index(nb)%ii(ix)  
         j = Atm_block%index(nb)%jj(ix)

         if (move_physics) then
            do k = 1, IPD_Control%lsoil            
               mn_phys%smc(i,j,k) = IPD_Data(nb)%Sfcprop%smc(ix,k)
               mn_phys%stc(i,j,k) = IPD_Data(nb)%Sfcprop%stc(ix,k)
               mn_phys%slc(i,j,k) = IPD_Data(nb)%Sfcprop%slc(ix,k)
            enddo

            mn_phys%u10m(i,j)  = IPD_Data(nb)%IntDiag%u10m(ix)
            mn_phys%v10m(i,j)  = IPD_Data(nb)%IntDiag%v10m(ix)
            mn_phys%tprcp(i,j)  = IPD_Data(nb)%Sfcprop%tprcp(ix)

            do k = 1, IPD_Control%nmtvr
               mn_phys%hprime(i,j,k)  = IPD_Data(nb)%Sfcprop%hprime(ix,k)
            end do
            
            mn_phys%zorl(i,j)  = IPD_Data(nb)%Sfcprop%zorl(ix)
            mn_phys%alvsf(i,j) = IPD_Data(nb)%Sfcprop%alvsf(ix)
            mn_phys%alvwf(i,j) = IPD_Data(nb)%Sfcprop%alvwf(ix)
            mn_phys%alnsf(i,j) = IPD_Data(nb)%Sfcprop%alnsf(ix)
            mn_phys%alnwf(i,j) = IPD_Data(nb)%Sfcprop%alnwf(ix)
            
            mn_phys%facsf(i,j) = IPD_data(nb)%Sfcprop%facsf(ix)   ! fractional coverage for strong zenith angle albedo
            mn_phys%facwf(i,j) = IPD_data(nb)%Sfcprop%facwf(ix)   ! fractional coverage for weak zenith angle albedo

            mn_phys%canopy(i,j) = IPD_Data(nb)%Sfcprop%canopy(ix)
            mn_phys%vegfrac(i,j)= IPD_Data(nb)%Sfcprop%vfrac(ix)
            mn_phys%uustar(i,j) = IPD_Data(nb)%Sfcprop%uustar(ix)
            mn_phys%shdmin(i,j) = IPD_Data(nb)%Sfcprop%shdmin(ix)
            mn_phys%shdmax(i,j) = IPD_Data(nb)%Sfcprop%shdmax(ix)
            mn_phys%zorll(i,j)  = IPD_Data(nb)%Sfcprop%zorll(ix)
            mn_phys%zorlo(i,j)  = IPD_Data(nb)%Sfcprop%zorlo(ix)
            mn_phys%zorlw(i,j)  = IPD_Data(nb)%Sfcprop%zorlw(ix)
            mn_phys%tsfco(i,j)  = IPD_Data(nb)%Sfcprop%tsfco(ix)
            mn_phys%tsfcl(i,j)  = IPD_Data(nb)%Sfcprop%tsfcl(ix)
            mn_phys%tsfc(i,j)   = IPD_Data(nb)%Sfcprop%tsfc(ix)

            do nv = 1, IPD_Control%ntot2d
               mn_phys%phy_f2d(i,j,nv) = IPD_Data(nb)%Tbd%phy_f2d(ix, nv)
            end do
            
            do k = 1, IPD_Control%levs
               do nv = 1, IPD_Control%ntot3d
                  mn_phys%phy_f3d(i,j,k,nv) = IPD_Data(nb)%Tbd%phy_f3d(ix, k, nv)
               end do
            end do
            
            ! Cloud prop data has x,y dimensions
            mn_phys%cv(i,j)  = IPD_Data(nb)%Cldprop%cv(ix)
            mn_phys%cvt(i,j) = IPD_Data(nb)%Cldprop%cvt(ix)
            mn_phys%cvb(i,j) = IPD_Data(nb)%Cldprop%cvb(ix)
         end if
         
         if (move_nsst) then
            mn_phys%tref(i,j)   = IPD_Data(nb)%Sfcprop%tref(ix)
            mn_phys%z_c(i,j)    = IPD_Data(nb)%Sfcprop%z_c(ix)
            mn_phys%c_0(i,j)    = IPD_Data(nb)%Sfcprop%c_0(ix)
            mn_phys%c_d(i,j)    = IPD_Data(nb)%Sfcprop%c_d(ix)
            mn_phys%w_0(i,j)    = IPD_Data(nb)%Sfcprop%w_0(ix)
            mn_phys%w_d(i,j)    = IPD_Data(nb)%Sfcprop%w_d(ix)
            mn_phys%xt(i,j)     = IPD_Data(nb)%Sfcprop%xt(ix) 
            mn_phys%xs(i,j)     = IPD_Data(nb)%Sfcprop%xs(ix)
            mn_phys%xu(i,j)     = IPD_Data(nb)%Sfcprop%xu(ix)
            mn_phys%xv(i,j)     = IPD_Data(nb)%Sfcprop%xv(ix)
            mn_phys%xz(i,j)     = IPD_Data(nb)%Sfcprop%xz(ix)
            mn_phys%zm(i,j)     = IPD_Data(nb)%Sfcprop%zm(ix)
            mn_phys%xtts(i,j)   = IPD_Data(nb)%Sfcprop%xtts(ix)
            mn_phys%xzts(i,j)   = IPD_Data(nb)%Sfcprop%xzts(ix)
            mn_phys%d_conv(i,j) = IPD_Data(nb)%Sfcprop%d_conv(ix)
            !mn_phys%ifd(i,j)    = IPD_Data(nb)%Sfcprop%ifd(ix)
            mn_phys%dt_cool(i,j)= IPD_Data(nb)%Sfcprop%dt_cool(ix)
            mn_phys%qrain(i,j)  = IPD_Data(nb)%Sfcprop%qrain(ix)
         end if
      enddo
    enddo

    if (debug_log) print '("[INFO] WDR end mn_phys_fill_temp_variables. npe=",I0," n=",I0)', this_pe, n
    
  end subroutine mn_phys_fill_temp_variables

  subroutine mn_prog_apply_temp_variables(Atm, n, child_grid_num, is_fine_pe, npz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)
    integer, intent(in)                                      :: n, child_grid_num
    logical, intent(in)                                      :: is_fine_pe
    integer, intent(in)                                      :: npz

    integer :: is, ie, js, je
    integer :: this_pe
    integer :: i,j,k
    integer :: bad_values, good_values
    type(fv_moving_nest_prog_type), pointer :: mn_prog

    mn_prog => Atm(n)%mn_prog

    this_pe = mpp_pe()

    if (debug_log) print '("[INFO] WDR start mn_prog_apply_temp_variables. npe=",I0," n=",I0)', this_pe, n
    
    ! Check if the variables were filled in properly.
    
    if (debug_log) then
       good_values = 0
       bad_values = 0
       
       if (is_fine_pe) then
          do i = Atm(n)%bd%isd, Atm(n)%bd%ied
             do j = Atm(n)%bd%jsd, Atm(n)%bd%jed
                do k = 1, npz
                   if (mn_prog%delz(i,j,k) .gt. 20000.0) then
                      print '("[WARN] WDR BAD NEST mn_prog%delz value. npe=",I0," mn_prog%delz(",I0,",",I0,",",I0,")=",F12.3)', this_pe, i, j, k, mn_prog%delz(i,j,k)
                      bad_values = bad_values + 1
                   else
                      good_values = good_values + 1
                   end if
                end do
             end do
          end do
       else
          do i = Atm(n)%bd%is, Atm(n)%bd%ie
             do j = Atm(n)%bd%js, Atm(n)%bd%je
                do k = 1, npz
                   if (mn_prog%delz(i,j,k) .gt. 20000.0) then
                      print '("[WARN] WDR BAD GLOBAL mn_prog%delz value. npe=",I0," mn_prog%delz(",I0,",",I0,",",I0,")=",F12.3)', this_pe, i, j, k, mn_prog%delz(i,j,k)
                      bad_values = bad_values + 1
                   else
                      good_values = good_values + 1
                   end if
                end do
             end do
          end do
       end if
       

       i = Atm(n)%bd%is
       j = Atm(n)%bd%js
       k = npz
       
       print '("[WARN] WDR Surface mn_prog%delz value. npe=",I0," mn_prog%delz(",I0,",",I0,",",I0,")=",F18.3)', this_pe, i, j, k, mn_prog%delz(i,j,k)

       print '("INFO] WDR mn_prog%delz values. npe=",I0," good_values=",I0," bad_values=",I0)', this_pe, good_values, bad_values
    end if
       

    if (is_fine_pe) then
       is = Atm(n)%bd%is
       ie = Atm(n)%bd%ie
       js = Atm(n)%bd%js
       je = Atm(n)%bd%je

       if (debug_log) print '("[INFO] WDR mn_prog_apply_temp_variables. npe=",I0," is=",I0," ie=",I0," js=",I0," je=",I0)', this_pe, is, ie, js, je


       Atm(n)%delz(is:ie, js:je, 1:npz) =  mn_prog%delz(is:ie, js:je, 1:npz) 
    end if

    if (debug_log) print '("[INFO] WDR end mn_prog_apply_temp_variables. npe=",I0," n=",I0)', this_pe, n

  end subroutine mn_prog_apply_temp_variables

  subroutine mn_phys_apply_temp_variables(Atm, Atm_block, IPD_Control, IPD_Data, n, child_grid_num, is_fine_pe, npz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)
    type (block_control_type), intent(in)            :: Atm_block
    type(IPD_control_type), intent(in)               :: IPD_Control
    type(IPD_data_type), intent(inout)               :: IPD_Data(:)
    integer, intent(in)                              :: n, child_grid_num
    logical, intent(in)                              :: is_fine_pe
    integer, intent(in)                              :: npz

    integer :: is, ie, js, je
    integer :: this_pe
    integer :: nb, blen, i, j ,k, ix, nv
    integer :: bad_values, good_values
    type(fv_moving_nest_physics_type), pointer       :: mn_phys

    this_pe = mpp_pe()
    mn_phys => Atm(n)%mn_phys

    if (debug_log) print '("[INFO] WDR start mn_phys_apply_temp_variables. npe=",I0," n=",I0)', this_pe, n
    
    ! Check if the variables were filled in properly.
    
    if (debug_log) then
       good_values = 0
       bad_values = 0
       
       if (is_fine_pe) then
          do i = Atm(n)%bd%isd, Atm(n)%bd%ied
             do j = Atm(n)%bd%jsd, Atm(n)%bd%jed
                if (mn_phys%ts(i,j) .gt. 20000.0) then
                   print '("[WARN] WDR BAD NEST ts value. npe=",I0," ts(",I0,",",I0,")=",F12.3)', this_pe, i, j, mn_phys%ts(i,j)
                   bad_values = bad_values + 1
                else
                   good_values = good_values + 1
                end if
             end do
          end do
       else
          do i = Atm(n)%bd%is, Atm(n)%bd%ie
             do j = Atm(n)%bd%js, Atm(n)%bd%je
                if (mn_phys%ts(i,j) .gt. 20000.0) then
                   print '("[WARN] WDR BAD GLOBAL ts value. npe=",I0," ts(",I0,",",I0")=",F12.3)', this_pe, i, j, mn_phys%ts(i,j)
                   bad_values = bad_values + 1
                else
                   good_values = good_values + 1
                end if
             end do
          end do
       end if
       

       i = Atm(n)%bd%is
       j = Atm(n)%bd%js
       
       print '("[WARN] WDR Surface ts value. npe=",I0," ts(",I0,",",I0,")=",F18.3)', this_pe, i, j, mn_phys%ts(i,j)

       print '("INFO] WDR ts values. npe=",I0," good_values=",I0," bad_values=",I0)', this_pe, good_values, bad_values
    end if
       
    !  Needed to fill the local grids for parent and nest PEs in order to transmit/interpolate data from parent to nest
    !  But only the nest PE's have changed the values with nest motion, so they are the only ones that need to update the original arrays
    if (is_fine_pe) then
       is = Atm(n)%bd%is
       ie = Atm(n)%bd%ie
       js = Atm(n)%bd%js
       je = Atm(n)%bd%je

       if (debug_log) print '("[INFO] WDR mn_phys_apply_temp_variables. npe=",I0," is=",I0," ie=",I0," js=",I0," je=",I0)', this_pe, is, ie, js, je

       ! SST directly in Atm structure
       Atm(n)%ts(is:ie, js:je) =  mn_phys%ts(is:ie, js:je) 
       
       do nb = 1,Atm_block%nblks
          blen = Atm_block%blksz(nb)
          do ix = 1, blen
             i = Atm_block%index(nb)%ii(ix)
             j = Atm_block%index(nb)%jj(ix)

             if (move_physics) then
                ! Surface properties
                do k = 1, IPD_Control%lsoil
                   IPD_Data(nb)%Sfcprop%smc(ix,k) = mn_phys%smc(i,j,k)
                   IPD_Data(nb)%Sfcprop%stc(ix,k) = mn_phys%stc(i,j,k)
                   IPD_Data(nb)%Sfcprop%slc(ix,k) = mn_phys%slc(i,j,k)
                enddo

                IPD_Data(nb)%IntDiag%u10m(ix) = mn_phys%u10m(i,j)
                IPD_Data(nb)%IntDiag%v10m(ix) = mn_phys%v10m(i,j)
                IPD_Data(nb)%Sfcprop%tprcp(ix) = mn_phys%tprcp(i,j)

                do k = 1, IPD_Control%nmtvr
                   IPD_Data(nb)%Sfcprop%hprime(ix,k) = mn_phys%hprime(i,j,k)
                end do
                
                IPD_Data(nb)%Sfcprop%zorl(ix) = mn_phys%zorl(i,j)
                IPD_Data(nb)%Sfcprop%alvsf(ix) = mn_phys%alvsf(i,j)
                IPD_Data(nb)%Sfcprop%alvwf(ix) = mn_phys%alvwf(i,j)
                IPD_Data(nb)%Sfcprop%alnsf(ix) = mn_phys%alnsf(i,j)
                IPD_Data(nb)%Sfcprop%alnwf(ix) = mn_phys%alnwf(i,j)

                IPD_Data(nb)%Sfcprop%facsf(ix) = mn_phys%facsf(i,j)
                IPD_Data(nb)%Sfcprop%facwf(ix) = mn_phys%facwf(i,j)
                
                IPD_Data(nb)%Sfcprop%canopy(ix) = mn_phys%canopy(i,j)
                IPD_Data(nb)%Sfcprop%vfrac(ix)  = mn_phys%vegfrac(i,j)
                IPD_Data(nb)%Sfcprop%uustar(ix) = mn_phys%uustar(i,j)
                IPD_Data(nb)%Sfcprop%shdmin(ix) = mn_phys%shdmin(i,j)
                IPD_Data(nb)%Sfcprop%shdmax(ix) = mn_phys%shdmax(i,j)
                IPD_Data(nb)%Sfcprop%zorll(ix)  = mn_phys%zorll(i,j)
                IPD_Data(nb)%Sfcprop%zorlo(ix)  = mn_phys%zorlo(i,j)
                IPD_Data(nb)%Sfcprop%zorlw(ix)  = mn_phys%zorlw(i,j)
                IPD_Data(nb)%Sfcprop%tsfco(ix)  = mn_phys%tsfco(i,j)
                IPD_Data(nb)%Sfcprop%tsfcl(ix)  = mn_phys%tsfcl(i,j)
                IPD_Data(nb)%Sfcprop%tsfc(ix)   = mn_phys%tsfc(i,j)
                
                ! Cloud properties
                IPD_Data(nb)%Cldprop%cv(ix) = mn_phys%cv(i,j)
                IPD_Data(nb)%Cldprop%cvt(ix) = mn_phys%cvt(i,j)
                IPD_Data(nb)%Cldprop%cvb(ix) = mn_phys%cvb(i,j)
                
                do nv = 1, IPD_Control%ntot2d
                   IPD_Data(nb)%Tbd%phy_f2d(ix, nv) = mn_phys%phy_f2d(i,j,nv)
                enddo
                
                do k = 1, IPD_Control%levs
                   do nv = 1, IPD_Control%ntot3d
                      IPD_Data(nb)%Tbd%phy_f3d(ix, k, nv) = mn_phys%phy_f3d(i,j,k,nv)
                   enddo
                enddo
             end if

             if (move_nsst) then
                IPD_Data(nb)%Sfcprop%tref(ix)    = mn_phys%tref(i,j)
                IPD_Data(nb)%Sfcprop%z_c(ix)     = mn_phys%z_c(i,j)
                IPD_Data(nb)%Sfcprop%c_0(ix)     = mn_phys%c_0(i,j)
                IPD_Data(nb)%Sfcprop%c_d(ix)     = mn_phys%c_d(i,j)
                IPD_Data(nb)%Sfcprop%w_0(ix)     = mn_phys%w_0(i,j)
                IPD_Data(nb)%Sfcprop%w_d(ix)     = mn_phys%w_d(i,j)
                IPD_Data(nb)%Sfcprop%xt(ix)      = mn_phys%xt(i,j)
                IPD_Data(nb)%Sfcprop%xs(ix)      = mn_phys%xs(i,j)
                IPD_Data(nb)%Sfcprop%xu(ix)      = mn_phys%xu(i,j)
                IPD_Data(nb)%Sfcprop%xv(ix)      = mn_phys%xv(i,j)
                IPD_Data(nb)%Sfcprop%xz(ix)      = mn_phys%xz(i,j)
                IPD_Data(nb)%Sfcprop%zm(ix)      = mn_phys%zm(i,j)
                IPD_Data(nb)%Sfcprop%xtts(ix)    = mn_phys%xtts(i,j)
                IPD_Data(nb)%Sfcprop%xzts(ix)    = mn_phys%xzts(i,j)
                IPD_Data(nb)%Sfcprop%d_conv(ix)  = mn_phys%d_conv(i,j)
                !IPD_Data(nb)%Sfcprop%ifd(ix)    = mn_phys%ifd(i,j)
                IPD_Data(nb)%Sfcprop%dt_cool(ix) = mn_phys%dt_cool(i,j)
                IPD_Data(nb)%Sfcprop%qrain(ix)   = mn_phys%qrain(i,j)
             end if
          enddo
       enddo       
   end if
      
   if (debug_log) print '("[INFO] WDR end mn_phys_apply_temp_variables. npe=",I0," n=",I0)', this_pe, n
   
  end subroutine mn_phys_apply_temp_variables


  !!===================================================================================== 
  !! Step 2 -- Fill the nest edge halos from parent grid before nest motion 
  !!            OR Refill the nest edge halos from parent grid after nest motion   
  !!            Parent and nest PEs need to execute these subroutines
  !!=====================================================================================           

  !  TODO clarify the child_grid_num or child_grid_level to handle multiple levels of nesting
  subroutine mn_prog_fill_nest_halos_from_parent(Atm, n, child_grid_num, is_fine_pe, nest_domain, nz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)
    integer, intent(in)                                      :: n, child_grid_num
    logical, intent(in)                                      :: is_fine_pe
    type(nest_domain_type), intent(inout)                    :: nest_domain
    integer, intent(in)                                      :: nz


    integer  :: position, position_u, position_v
    integer  :: interp_type, interp_type_u, interp_type_v
    integer  :: x_refine, y_refine
    type(fv_moving_nest_prog_type), pointer :: mn_prog

    mn_prog => Atm(n)%mn_prog

    !  TODO: examine how the static nesting code handles this
    !  TODO move terrain and surface parameters, including phis

    !  TODO Rename this from interp_type to stagger_type
    interp_type = 1    ! cell-centered A-grid
    interp_type_u = 4  ! D-grid
    interp_type_v = 4  ! D-grid

    position = CENTER ! CENTER, NORTH, EAST
    position_u = NORTH
    position_v = EAST

    x_refine = Atm(child_grid_num)%neststruct%refinement
    y_refine = x_refine

    !  Fill centered-grid variables
    call fill_nest_halos_from_parent("q_con", Atm(n)%q_con, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
         Atm(child_grid_num)%neststruct%ind_h, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position, nz)

    call fill_nest_halos_from_parent("pt", Atm(n)%pt, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
         Atm(child_grid_num)%neststruct%ind_h, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position, nz)

    call fill_nest_halos_from_parent("w", Atm(n)%w, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
         Atm(child_grid_num)%neststruct%ind_h, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position, nz)

    !call fill_nest_halos_from_parent("omga", Atm(n)%omga, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
    !     Atm(child_grid_num)%neststruct%ind_h, &
    !     x_refine, y_refine, &
    !     is_fine_pe, nest_domain, position, nz)

    call fill_nest_halos_from_parent("delp", Atm(n)%delp, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
         Atm(child_grid_num)%neststruct%ind_h, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position, nz)

    !call fill_nest_halos_from_parent("delz", Atm(n)%delz, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
    !     Atm(child_grid_num)%neststruct%ind_h, &
    !     x_refine, y_refine, &
    !     is_fine_pe, nest_domain, position, nz)

    call fill_nest_halos_from_parent("delz", mn_prog%delz, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
         Atm(child_grid_num)%neststruct%ind_h, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position, nz)


    call fill_nest_halos_from_parent("q", Atm(n)%q, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
         Atm(child_grid_num)%neststruct%ind_h, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position, nz)


    !  move the a-grid winds.  TODO consider recomputing them from D grid instead
    call fill_nest_halos_from_parent("ua", Atm(n)%ua, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
         Atm(child_grid_num)%neststruct%ind_h, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position, nz)

    call fill_nest_halos_from_parent("va", Atm(n)%va, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
         Atm(child_grid_num)%neststruct%ind_h, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position, nz)



    !  Fill staggered D-grid variables
    call fill_nest_halos_from_parent("u", Atm(n)%u, interp_type_u, Atm(child_grid_num)%neststruct%wt_u, &
         Atm(child_grid_num)%neststruct%ind_u, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position_u, nz)

    call fill_nest_halos_from_parent("v", Atm(n)%v, interp_type_v, Atm(child_grid_num)%neststruct%wt_v, &
         Atm(child_grid_num)%neststruct%ind_v, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position_v, nz)


    ! from fv_array.F90
    !allocate (    Atm%u(isd:ied  ,jsd:jed+1,npz) ) ; Atm%u=real_snan
    !allocate (    Atm%v(isd:ied+1,jsd:jed  ,npz) ) ; Atm%v=real_snan
    !allocate (   Atm%pt(isd:ied  ,jsd:jed  ,npz) ) ; Atm%pt=real_snan
    !allocate ( Atm%delp(isd:ied  ,jsd:jed  ,npz) ) ; Atm%delp=real_snan
    !allocate (    Atm%q(isd:ied  ,jsd:jed  ,npz, nq) ) ; Atm%q=real_snan
    !allocate (Atm%qdiag(isd:ied  ,jsd:jed  ,npz, nq+1:ncnst) ) ; Atm%qdiag=real_snan

  end subroutine mn_prog_fill_nest_halos_from_parent


  subroutine mn_phys_fill_nest_halos_from_parent(Atm, IPD_Control, IPD_Data, n, child_grid_num, is_fine_pe, nest_domain, nz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)
    type(IPD_control_type), intent(in)               :: IPD_Control
    type(IPD_data_type), intent(inout)               :: IPD_Data(:)
    integer, intent(in)                              :: n, child_grid_num
    logical, intent(in)                              :: is_fine_pe
    type(nest_domain_type), intent(inout)            :: nest_domain
    integer, intent(in)                              :: nz


    integer  :: position, position_u, position_v
    integer  :: interp_type, interp_type_u, interp_type_v
    integer  :: x_refine, y_refine
    type(fv_moving_nest_physics_type), pointer :: mn_phys

    !  TODO: examine how the static nesting code handles this
    !  TODO move terrain and surface parameters, including phis

    interp_type = 1    ! cell-centered A-grid
    interp_type_u = 4  ! D-grid
    interp_type_v = 4  ! D-grid

    position = CENTER ! CENTER, NORTH, EAST
    position_u = NORTH
    position_v = EAST

    x_refine = Atm(child_grid_num)%neststruct%refinement
    y_refine = x_refine

    mn_phys => Atm(n)%mn_phys

    !  Fill centered-grid variables

    call fill_nest_halos_from_parent("ts", mn_phys%ts, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
         Atm(child_grid_num)%neststruct%ind_h, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position)

    if (move_physics) then
       call fill_nest_halos_from_parent("smc", mn_phys%smc, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position, IPD_Control%lsoil)
       call fill_nest_halos_from_parent("stc", mn_phys%stc, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position, IPD_Control%lsoil)
       call fill_nest_halos_from_parent("slc", mn_phys%slc, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position, IPD_Control%lsoil)
       
       call fill_nest_halos_from_parent("phy_f2d", mn_phys%phy_f2d, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position, IPD_Control%ntot2d)
       
       call fill_nest_halos_from_parent("phy_f3d", mn_phys%phy_f3d, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position, IPD_Control%levs)
       
       !!  Surface variables
       call fill_nest_halos_from_parent("u10m", mn_phys%u10m, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("v10m", mn_phys%v10m, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("tprcp", mn_phys%tprcp, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)

       call fill_nest_halos_from_parent("hprime", mn_phys%hprime, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position, IPD_Control%nmtvr)
       
       call fill_nest_halos_from_parent("zorl", mn_phys%zorl, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       
       call fill_nest_halos_from_parent("alvsf", mn_phys%alvsf, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       
       call fill_nest_halos_from_parent("alvwf", mn_phys%alvwf, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       
       call fill_nest_halos_from_parent("alnsf", mn_phys%alnsf, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       
       call fill_nest_halos_from_parent("alnwf", mn_phys%alnwf, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)


       call fill_nest_halos_from_parent("facsf", mn_phys%facsf, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)

       call fill_nest_halos_from_parent("facwf", mn_phys%facwf, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       
       !!
       
       call fill_nest_halos_from_parent("canopy", mn_phys%canopy, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("vegfrac", mn_phys%vegfrac, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("uustar", mn_phys%uustar, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("shdmin", mn_phys%shdmin, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("shdmax", mn_phys%shdmax, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("zorll", mn_phys%zorll, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("zorlo", mn_phys%zorlo, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("zorlw", mn_phys%zorlw, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)


       call fill_nest_halos_from_parent("tsfco", mn_phys%tsfco, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("tsfcl", mn_phys%tsfcl, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("tsfc", mn_phys%tsfc, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       
       
       call fill_nest_halos_from_parent("cv", mn_phys%cv, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("cvt", mn_phys%cvt, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("cvb", mn_phys%cvb, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
    end if

    if (move_nsst) then

       call fill_nest_halos_from_parent("tref", mn_phys%tref, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("z_c", mn_phys%z_c, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("c_0", mn_phys%c_0, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("c_d", mn_phys%c_d, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("w_0", mn_phys%w_0, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("w_d", mn_phys%w_d, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("xt", mn_phys%xt, interp_type, Atm(child_grid_num)%neststruct%wt_h, & 
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("xs", mn_phys%xs, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("xu", mn_phys%xu, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("xv", mn_phys%xv, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("xz", mn_phys%xz, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("zm", mn_phys%zm, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("xtts", mn_phys%xtts, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("xzts", mn_phys%xzts, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("d_conv", mn_phys%d_conv, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       !call fill_nest_halos_from_parent("ifd", mn_phys%ifd, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
       !     Atm(child_grid_num)%neststruct%ind_h, &
       !     x_refine, y_refine, &
       !     is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("dt_cool", mn_phys%dt_cool, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call fill_nest_halos_from_parent("qrain", mn_phys%qrain, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
            Atm(child_grid_num)%neststruct%ind_h, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)

       end if

  end subroutine mn_phys_fill_nest_halos_from_parent



  !!============================================================================                                            
  !! Step 3 -- Redefine the nest domain to new location                                                                     
  !!   This calls mpp_define_nest_domains.  Following the code in fv_control.F90, only should                               
  !!   be executed on the nest PEs. Operates only on indices.                                                               
  !!  --  Similar to med_nest_configure() from HWRF                                                                         
  !!============================================================================             

  subroutine mn_meta_move_nest(delta_i_c, delta_j_c, pelist, is_fine_pe, extra_halo, &
       nest_domain, domain_fine, domain_coarse, tile_fine, tile_coarse, &
       istart_coarse, iend_coarse, jstart_coarse, jend_coarse,  istart_fine, iend_fine, jstart_fine, jend_fine)

    implicit none

    integer, intent(in)                   :: delta_i_c, delta_j_c
    integer, allocatable, intent(in)      :: pelist(:)
    logical, intent(in)                   :: is_fine_pe
    integer, intent(in)                   :: extra_halo

    type(nest_domain_type), intent(inout) :: nest_domain
    type(domain2d), intent(inout)         :: domain_coarse, domain_fine
    integer, intent(inout)                :: tile_coarse, tile_fine
    integer, intent(inout)                :: istart_coarse, iend_coarse, jstart_coarse, jend_coarse
    integer, intent(in)                   :: istart_fine, iend_fine, jstart_fine, jend_fine

    !  Local variables
    integer   :: num_nest
    integer   :: this_pe

    integer   :: delta_i_coarse(1), delta_j_coarse(1)

    this_pe = mpp_pe()

    if (debug_log) print '("[INFO] WDR start mn_meta_move_nest. npe=",I0)', this_pe


    !  Initial implementation only supports single moving nest.  Update this later.
    !  mpp_shift_nest_domains has a call signature to support multiple moving nests, though has not been tested for correctness.
    delta_i_coarse(1) = delta_i_c
    delta_j_coarse(1) = delta_j_c


    !!===========================================================
    !!
    !! Relocate the nest points
    !!
    !!===========================================================

    istart_coarse = istart_coarse + delta_i_c
    iend_coarse = iend_coarse + delta_i_c

    jstart_coarse = jstart_coarse + delta_j_c
    jend_coarse = jend_coarse + delta_j_c

    ! the fine nest should maintain its indices

    if (debug_log) print '("[INFO] WDR NRD1 about to call mpp_define_nest_domains. npe=",I0," ",I0," ",I0," ",I0," ",I0)', this_pe, istart_coarse, iend_coarse, istart_fine, iend_fine 

    !!===========================================================
    !!
    !! Looks like this is safe to call repeatedly without zapping atmospheric fields; 
    !!   it sets the mapping between coarse and fine tiles
    !!
    !!===========================================================

    !  OLD Dycore form
    !  call mpp_define_nest_domains(nest_domain, domain_fine, domain_coarse, tile_fine, tile_coarse, &
    !       istart_fine, iend_fine, jstart_fine, jend_fine,                  &
    !       istart_coarse, iend_coarse, jstart_coarse, jend_coarse,         &
    !       pelist, extra_halo, name="nest_domain")


    ! New Dycore
    !   type nest_domain_type
    !     character(len=NAME_LENGTH)     :: name
    !     integer                        :: num_level
    !     type(nest_level_type), pointer :: nest(:) => NULL()
    !     integer                        :: num_nest
    !     integer,               pointer :: tile_fine(:), tile_coarse(:)
    !     integer,               pointer :: istart_fine(:), iend_fine(:), jstart_fine(:), jend_fine(:)
    !     integer,               pointer :: istart_coarse(:), iend_coarse(:), jstart_coarse(:), jend_coarse(:)
    !  end type nest_domain_type

    num_nest = nest_domain%num_nest
    ! Which nest are we altering here?


    ! WDR TODO Verify whether rerunning this will cause (small) memory leaks. 

    if (is_fine_pe) then
       call mpp_shift_nest_domains(nest_domain, domain_fine, delta_i_coarse, delta_j_coarse, extra_halo)
    else
       call mpp_shift_nest_domains(nest_domain, domain_coarse, delta_i_coarse, delta_j_coarse, extra_halo)
    end if



    ! New dycore, from fv_control.F90
    !call mpp_define_nest_domains(global_nest_domain, Atm(this_grid)%domain, &
    !     ngrids-1, nest_level=nest_level(2:ngrids) , &
    !     istart_coarse=nest_ioffsets(2:ngrids), jstart_coarse=nest_joffsets(2:ngrids), &
    !     icount_coarse=icount_coarse(2:ngrids), jcount_coarse=jcount_coarse(2:ngrids), &
    !     npes_nest_tile=npes_nest_tile(1:ntiles_nest_all), &
    !     tile_fine=tile_fine(2:ngrids), tile_coarse=tile_coarse(2:ngrids), &
    !     x_refine=nest_refine(2:ngrids), y_refine=nest_refine(2:ngrids), name="global_nest_domain")



    if (debug_log) print '("[INFO] WDR NRD2 after call to mpp_define_nest_domains. npe=",I0)', this_pe

  end subroutine mn_meta_move_nest


  !================================================================================ 
  !! Step 4 --  Updates the internal nest tile halos
  !================================================================================ 

  ! Fill internal nest halos for prognostic variables
  subroutine mn_prog_fill_intern_nest_halos(Atm, domain_fine, is_fine_pe)
    type(fv_atmos_type), target, intent(inout)  :: Atm  
    type(domain2d), intent(inout)               :: domain_fine
    logical, intent(in)                         :: is_fine_pe


    integer :: this_pe
    type(fv_moving_nest_prog_type), pointer :: mn_prog

    mn_prog => Atm%mn_prog


    this_pe = mpp_pe()

    call mn_var_fill_intern_nest_halos(Atm%q_con, domain_fine, is_fine_pe)
    call mn_var_fill_intern_nest_halos(Atm%pt, domain_fine, is_fine_pe)
    call mn_var_fill_intern_nest_halos(Atm%w, domain_fine, is_fine_pe)
    !call mn_var_fill_intern_nest_halos(Atm%omga, domain_fine, is_fine_pe)
    call mn_var_fill_intern_nest_halos(Atm%delp, domain_fine, is_fine_pe)
    !call mn_var_fill_intern_nest_halos(Atm%delz, domain_fine, is_fine_pe)
    call mn_var_fill_intern_nest_halos(mn_prog%delz, domain_fine, is_fine_pe)

    call mn_var_fill_intern_nest_halos(Atm%ua, domain_fine, is_fine_pe)
    call mn_var_fill_intern_nest_halos(Atm%va, domain_fine, is_fine_pe)

    if (debug_log) then
       call check_array(Atm%u, this_pe, "Atm%u", -300.0, 300.0)
       call check_array(Atm%v, this_pe, "Atm%v", -300.0, 300.0)
    end if

    ! The vector form of the subroutine takes care of the staggering of the wind variables internally.
    call mn_var_fill_intern_nest_halos(Atm%u, Atm%v, domain_fine, is_fine_pe)


    call mn_var_fill_intern_nest_halos(Atm%q, domain_fine, is_fine_pe)

  end subroutine mn_prog_fill_intern_nest_halos

  ! Fill internal nest halos for physics variables
  ! 
  subroutine mn_phys_fill_intern_nest_halos(Atm, IPD_Control, IPD_Data, domain_fine, is_fine_pe)
    type(fv_atmos_type), target, intent(inout)       :: Atm  
    type(IPD_control_type), intent(in)               :: IPD_Control
    type(IPD_data_type), intent(inout)               :: IPD_Data(:)
    type(domain2d), intent(inout)                    :: domain_fine
    logical, intent(in)                              :: is_fine_pe

    type(fv_moving_nest_physics_type), pointer :: mn_phys

    mn_phys => Atm%mn_phys

    call mn_var_fill_intern_nest_halos(mn_phys%ts, domain_fine, is_fine_pe)   !! Skin Temp/SST
    if (move_physics) then
       call mn_var_fill_intern_nest_halos(mn_phys%smc, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%stc, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%slc, domain_fine, is_fine_pe)
       
       call mn_var_fill_intern_nest_halos(mn_phys%phy_f2d, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%phy_f3d, domain_fine, is_fine_pe)

       call mn_var_fill_intern_nest_halos(mn_phys%u10m, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%v10m, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%tprcp, domain_fine, is_fine_pe)

       call mn_var_fill_intern_nest_halos(mn_phys%hprime, domain_fine, is_fine_pe)
       
       call mn_var_fill_intern_nest_halos(mn_phys%zorl, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%alvsf, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%alvwf, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%alnsf, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%alnwf, domain_fine, is_fine_pe)

       call mn_var_fill_intern_nest_halos(mn_phys%facsf, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%facwf, domain_fine, is_fine_pe)
       
       call mn_var_fill_intern_nest_halos(mn_phys%canopy, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%vegfrac, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%uustar, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%shdmin, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%shdmax, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%zorll, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%zorlo, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%zorlw, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%tsfco, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%tsfcl, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%tsfc, domain_fine, is_fine_pe)
       
       call mn_var_fill_intern_nest_halos(mn_phys%cv, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%cvt, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%cvb, domain_fine, is_fine_pe)
    end if

    if (move_nsst) then
       call mn_var_fill_intern_nest_halos(mn_phys%tref, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%z_c, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%c_0, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%c_d, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%w_0, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%w_d, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%xt, domain_fine, is_fine_pe) 
       call mn_var_fill_intern_nest_halos(mn_phys%xs, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%xu, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%xv, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%xz, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%zm, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%xtts, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%xzts, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%d_conv, domain_fine, is_fine_pe)
       !call mn_var_fill_intern_nest_halos(mn_phys%ifd, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%dt_cool, domain_fine, is_fine_pe)
       call mn_var_fill_intern_nest_halos(mn_phys%qrain, domain_fine, is_fine_pe)
    end if
    !call check_array(Atm%u, this_pe, "Atm%pt", 100.0, 400.0)

  end subroutine mn_phys_fill_intern_nest_halos



  !================================================================================ 
  !
  !   Step 4 -- Per variable fill internal nest halos
  !
  !================================================================================ 

  subroutine mn_var_fill_intern_nest_halos2D(data_var, domain_fine, is_fine_pe)
    real, allocatable, intent(inout)            :: data_var(:,:)
    type(domain2d), intent(inout)               :: domain_fine
    logical, intent(in)                         :: is_fine_pe

    integer                      :: this_pe
    this_pe = mpp_pe()

    if (is_fine_pe) then
       if (debug_log) print '("[INFO] WDR INH2 before call to mpp_update_domains. npe=",I0)', this_pe
       ! mpp_update_domains fills the halo region of the fine grids for the interior of the nest.  
       ! The fine nest boundary with the coarse grid remains unchanged.
       ! seems that this only performs communication between fine nest PEs
       ! Just transfers halo data between tiles of same resolution -- doesn't perform any interpolation!
       call mpp_update_domains(data_var, domain_fine,  flags=NUPDATE + EUPDATE + SUPDATE + WUPDATE)

       if (debug_log) print '("[INFO] WDR INH2 after call to mpp_update_domains. npe=",I0)', this_pe
    end if

  end subroutine mn_var_fill_intern_nest_halos2D


  subroutine mn_var_fill_intern_nest_halos2D_kindphys(data_var, domain_fine, is_fine_pe)
    real(kind=kind_phys), allocatable, intent(inout)            :: data_var(:,:)
    type(domain2d), intent(inout)               :: domain_fine
    logical, intent(in)                         :: is_fine_pe

    integer                      :: this_pe
    this_pe = mpp_pe()

    if (is_fine_pe) then
       if (debug_log) print '("[INFO] WDR INH2p before call to mpp_update_domains. npe=",I0)', this_pe
       ! mpp_update_domains fills the halo region of the fine grids for the interior of the nest.  
       ! The fine nest boundary with the coarse grid remains unchanged.
       ! seems that this only performs communication between fine nest PEs
       ! Just transfers halo data between tiles of same resolution -- doesn't perform any interpolation!
       call mpp_update_domains(data_var, domain_fine,  flags=NUPDATE + EUPDATE + SUPDATE + WUPDATE)

       if (debug_log) print '("[INFO] WDR INH2p after call to mpp_update_domains. npe=",I0)', this_pe
    end if

  end subroutine mn_var_fill_intern_nest_halos2D_kindphys

  subroutine mn_var_fill_intern_nest_halos3D(data_var, domain_fine, is_fine_pe)
    real, allocatable, intent(inout)            :: data_var(:,:,:)
    type(domain2d), intent(inout)               :: domain_fine
    logical, intent(in)                         :: is_fine_pe

    integer                      :: this_pe
    this_pe = mpp_pe()

    if (is_fine_pe) then
       if (debug_log) print '("[INFO] WDR INH3 before call to mpp_update_domains. npe=",I0)', this_pe
       ! mpp_update_domains fills the halo region of the fine grids for the interior of the nest.  
       ! The fine nest boundary with the coarse grid remains unchanged.
       ! seems that this only performs communication between fine nest PEs
       ! Just transfers halo data between tiles of same resolution -- doesn't perform any interpolation!
       call mpp_update_domains(data_var, domain_fine,  flags=NUPDATE + EUPDATE + SUPDATE + WUPDATE)

       if (debug_log) print '("[INFO] WDR INH3 after call to mpp_update_domains. npe=",I0)', this_pe
    end if

  end subroutine mn_var_fill_intern_nest_halos3D

  subroutine mn_var_fill_intern_nest_halos3D_kindphys(data_var, domain_fine, is_fine_pe)
    real(kind=kind_phys), allocatable, intent(inout)            :: data_var(:,:,:)
    type(domain2d), intent(inout)               :: domain_fine
    logical, intent(in)                         :: is_fine_pe

    integer                      :: this_pe
    this_pe = mpp_pe()

    if (is_fine_pe) then
       if (debug_log) print '("[INFO] WDR INH3p before call to mpp_update_domains. npe=",I0)', this_pe
       ! mpp_update_domains fills the halo region of the fine grids for the interior of the nest.  
       ! The fine nest boundary with the coarse grid remains unchanged.
       ! seems that this only performs communication between fine nest PEs
       ! Just transfers halo data between tiles of same resolution -- doesn't perform any interpolation!
       call mpp_update_domains(data_var, domain_fine,  flags=NUPDATE + EUPDATE + SUPDATE + WUPDATE)

       if (debug_log) print '("[INFO] WDR INH3p after call to mpp_update_domains. npe=",I0)', this_pe
    end if

  end subroutine mn_var_fill_intern_nest_halos3D_kindphys


  subroutine mn_var_fill_intern_nest_halos_wind(u_var, v_var, domain_fine, is_fine_pe)
    real, allocatable, intent(inout)            :: u_var(:,:,:)
    real, allocatable, intent(inout)            :: v_var(:,:,:)
    type(domain2d), intent(inout)               :: domain_fine
    logical, intent(in)                         :: is_fine_pe

    integer                      :: this_pe
    this_pe = mpp_pe()

    if (is_fine_pe) then
       if (debug_log) print '("[INFO] WDR INH3W before call to mpp_update_domains. npe=",I0)', this_pe
       ! mpp_update_domains fills the halo region of the fine grids for the interior of the nest.  
       ! The fine nest boundary with the coarse grid remains unchanged.
       ! seems that this only performs communication between fine nest PEs
       ! Just transfers halo data between tiles of same resolution -- doesn't perform any interpolation!
       call mpp_update_domains(u_var, v_var, domain_fine,  flags=NUPDATE + EUPDATE + SUPDATE + WUPDATE, gridtype=DGRID_NE)

       if (debug_log) print '("[INFO] WDR INH3W after call to mpp_update_domains. npe=",I0)', this_pe
    end if

  end subroutine mn_var_fill_intern_nest_halos_wind


  subroutine mn_var_fill_intern_nest_halos4D(data_var, domain_fine, is_fine_pe)
    real, allocatable, intent(inout)            :: data_var(:,:,:,:)
    type(domain2d), intent(inout)               :: domain_fine
    logical, intent(in)                         :: is_fine_pe

    integer                      :: this_pe
    this_pe = mpp_pe()

    if (is_fine_pe) then
       if (debug_log) print '("[INFO] WDR INH4 before call to mpp_update_domains. npe=",I0)', this_pe
       ! mpp_update_domains fills the halo region of the fine grids for the interior of the nest.  
       ! The fine nest boundary with the coarse grid remains unchanged.
       ! seems that this only performs communication between fine nest PEs
       ! Just transfers halo data between tiles of same resolution -- doesn't perform any interpolation!
       call mpp_update_domains(data_var, domain_fine,  flags=NUPDATE + EUPDATE + SUPDATE + WUPDATE)

       if (debug_log) print '("[INFO] WDR INH4 after call to mpp_update_domains. npe=",I0)', this_pe
    end if

  end subroutine mn_var_fill_intern_nest_halos4D


  subroutine mn_var_fill_intern_nest_halos4D_kindphys(data_var, domain_fine, is_fine_pe)
    real(kind=kind_phys), allocatable, intent(inout)            :: data_var(:,:,:,:)
    type(domain2d), intent(inout)               :: domain_fine
    logical, intent(in)                         :: is_fine_pe

    integer                      :: this_pe
    this_pe = mpp_pe()

    if (is_fine_pe) then
       if (debug_log) print '("[INFO] WDR INH4 before call to mpp_update_domains. npe=",I0)', this_pe
       ! mpp_update_domains fills the halo region of the fine grids for the interior of the nest.  
       ! The fine nest boundary with the coarse grid remains unchanged.
       ! seems that this only performs communication between fine nest PEs
       ! Just transfers halo data between tiles of same resolution -- doesn't perform any interpolation!
       call mpp_update_domains(data_var, domain_fine,  flags=NUPDATE + EUPDATE + SUPDATE + WUPDATE)

       if (debug_log) print '("[INFO] WDR INH4 after call to mpp_update_domains. npe=",I0)', this_pe
    end if

  end subroutine mn_var_fill_intern_nest_halos4D_kindphys



  !!============================================================================                                            
  !! Step 5.1 -- Load the latlon data from NetCDF
  !!             update parent_geo, tile_geo*, p_grid*, n_grid*
  !!============================================================================                                            


  subroutine mn_latlon_load_parent(Atm, n, delta_i_c, delta_j_c, child_grid_num, parent_geo, tile_geo, tile_geo_u, tile_geo_v, fp_super_tile_geo, p_grid, n_grid, p_grid_u, n_grid_u, p_grid_v, n_grid_v)
    type(fv_atmos_type), allocatable, intent(in) :: Atm(:)
    integer, intent(in)                          :: n, delta_i_c, delta_j_c, child_grid_num
    type(grid_geometry), intent(inout)           :: parent_geo, tile_geo, tile_geo_u, tile_geo_v
    type(grid_geometry), intent(in)              :: fp_super_tile_geo
    real(kind=R_GRID), allocatable, intent(out)  :: p_grid(:,:,:), n_grid(:,:,:), p_grid_u(:,:,:), n_grid_u(:,:,:), p_grid_v(:,:,:), n_grid_v(:,:,:)

    logical, save  :: first_nest_move = .true.
    integer, save  :: p_istart_fine, p_iend_fine, p_jstart_fine, p_jend_fine

    integer  :: x, y, fp_i, fp_j
    integer  :: position, position_u, position_v
    integer  :: x_refine, y_refine 
    integer  :: nest_x, nest_y, parent_x, parent_y

    character(len=256) :: res_str

    integer :: this_pe

    this_pe = mpp_pe()

    position = CENTER ! CENTER, NORTH, EAST                                                                                   
    position_u = NORTH
    position_v = EAST

    x_refine = Atm(child_grid_num)%neststruct%refinement
    y_refine = x_refine

    !  Setup parent_geo with the values for the parent tile
    !  Note that lat/lon are stored in the model in RADIANS
    !  Only the netCDF files use degrees

    !parent_geo%lons = Atm(1)%grid_global(:,:,1,6)
    !parent_geo%lats = Atm(1)%grid_global(:,:,2,6)

    write(res_str, '(I0)'), Atm(1)%npx - 1 

    if (first_nest_move) then
       if (debug_log) print '("[INFO] WDR mn_latlon_load_parent READING static coarse file on npe=",I0)', this_pe
       call load_nest_latlons_from_nc(trim(Atm(child_grid_num)%neststruct%surface_dir) //  '/C' // trim(res_str) //  '_grid.tile6.nc', &
            Atm(1)%npx, Atm(1)%npy, 1, &
            parent_geo, &
            p_istart_fine, p_iend_fine, p_jstart_fine, p_jend_fine)
       first_nest_move = .false.
    !else
    !   print '("[INFO] WDR mn_latlon_load_parent SKIPPING static coarse file on npe=",I0)', this_pe
    end if

    parent_geo%nxp = Atm(1)%npx
    parent_geo%nyp = Atm(1)%npy

    parent_geo%nx = Atm(1)%npx - 1
    parent_geo%ny = Atm(1)%npy - 1

    if (debug_log) then
       call show_tile_geo(parent_geo, this_pe, "parent_geo")
       call show_atm_grids(Atm, n)
    end if

    ! Actually, is the nest in grid_global??

    !  Setup tile_geo with the values for the nest
    !  this loses the offset for halos, and starts at 1 instead of at -2
    !tile_geo%lons = Atm(n)%grid_global(:,:,1,1)
    !tile_geo%lats = Atm(n)%grid_global(:,:,2,1)


    !===========================================================
    !  Begin tile_geo per PE.
    !===========================================================

    !------------------------
    ! Grid Definitions 
    !------------------------
    !
    ! tile_geo - lat/lons on A-grid (cell centers) for nest, on data domain (includes halo) for each PE
    ! parent_geo - lat/lons of supergrid for parent
    ! n_grid - lat/lons of cell centers for nest
    ! p_grid - lat/lons of cell centers for parent
    !
    ! gridstruct%agrid - cell centers for each PE
    ! gridstruct%grid - cell corners for each PE


    ! Allocate tile_geo just for this PE, copied from Atm(n)%gridstruct%agrid
    tile_geo%nx = ubound(Atm(n)%gridstruct%agrid, 1) - lbound(Atm(n)%gridstruct%agrid, 1)
    tile_geo%ny = ubound(Atm(n)%gridstruct%agrid, 2) - lbound(Atm(n)%gridstruct%agrid, 2)
    tile_geo%nxp = tile_geo%nx + 1
    tile_geo%nyp = tile_geo%ny + 1

    allocate(tile_geo%lons(lbound(Atm(n)%gridstruct%agrid, 1):ubound(Atm(n)%gridstruct%agrid, 1), lbound(Atm(n)%gridstruct%agrid, 2):ubound(Atm(n)%gridstruct%agrid, 2)))
    allocate(tile_geo%lats(lbound(Atm(n)%gridstruct%agrid, 1):ubound(Atm(n)%gridstruct%agrid, 1), lbound(Atm(n)%gridstruct%agrid, 2):ubound(Atm(n)%gridstruct%agrid, 2)))
    !allocate(tile_geo%area(lbound(Atm(n)%gridstruct%agrid, 1):ubound(Atm(n)%gridstruct%agrid, 1), lbound(Atm(n)%gridstruct%agrid, 2):ubound(Atm(n)%gridstruct%agrid, 2)))

    tile_geo%lats = -999.9
    tile_geo%lons = -999.9
    !tile_geo%area = -999.9

    do x = lbound(Atm(n)%gridstruct%agrid, 1), ubound(Atm(n)%gridstruct%agrid, 1)
       do y = lbound(Atm(n)%gridstruct%agrid, 2), ubound(Atm(n)%gridstruct%agrid, 2)
          tile_geo%lons(x,y) = Atm(n)%gridstruct%agrid(x,y,1)
          tile_geo%lats(x,y) = Atm(n)%gridstruct%agrid(x,y,2)
       end do
    end do

    if (debug_log) call show_tile_geo(tile_geo, this_pe, "tile_geo")
    call find_nest_alignment(tile_geo, fp_super_tile_geo, nest_x, nest_y, parent_x, parent_y)


    ! Allocate tile_geo_u just for this PE, copied from Atm(n)%gridstruct%grid
    ! grid is 1 larger than agrid
    ! u(npx, npy+1)   
    tile_geo_u%nx = ubound(Atm(n)%gridstruct%agrid, 1) - lbound(Atm(n)%gridstruct%agrid, 1)
    tile_geo_u%ny = ubound(Atm(n)%gridstruct%grid, 2) - lbound(Atm(n)%gridstruct%grid, 2)
    tile_geo_u%nxp = tile_geo_u%nx + 1
    tile_geo_u%nyp = tile_geo_u%ny + 1

    allocate(tile_geo_u%lons(lbound(Atm(n)%gridstruct%agrid, 1):ubound(Atm(n)%gridstruct%agrid, 1), lbound(Atm(n)%gridstruct%grid, 2):ubound(Atm(n)%gridstruct%grid, 2)))
    allocate(tile_geo_u%lats(lbound(Atm(n)%gridstruct%agrid, 1):ubound(Atm(n)%gridstruct%agrid, 1), lbound(Atm(n)%gridstruct%grid, 2):ubound(Atm(n)%gridstruct%grid, 2)))

    tile_geo_u%lons = -999.9
    tile_geo_u%lats = -999.9

    do x = lbound(tile_geo_u%lats, 1), ubound(tile_geo_u%lats, 1)
       do y = lbound(tile_geo_u%lats, 2), ubound(tile_geo_u%lats, 2)
          fp_i = (x - nest_x) * 2 + parent_x - 1
          fp_j = (y - nest_y) * 2 + parent_y 

          tile_geo_u%lons(x,y) = fp_super_tile_geo%lons(fp_i, fp_j)
          tile_geo_u%lats(x,y) = fp_super_tile_geo%lats(fp_i, fp_j)
       end do
    end do

    if (debug_log) call show_tile_geo(tile_geo_u, this_pe, "tile_geo_u")


    ! Allocate tile_geo_v just for this PE, copied from Atm(n)%gridstruct%grid
    ! grid is 1 larger than agrid
    ! u(npx, npy+1)   
    tile_geo_v%nx = ubound(Atm(n)%gridstruct%grid, 1) - lbound(Atm(n)%gridstruct%grid, 1)
    tile_geo_v%ny = ubound(Atm(n)%gridstruct%agrid, 2) - lbound(Atm(n)%gridstruct%agrid, 2)
    tile_geo_v%nxp = tile_geo_v%nx + 1
    tile_geo_v%nyp = tile_geo_v%ny + 1

    allocate(tile_geo_v%lons(lbound(Atm(n)%gridstruct%grid, 1):ubound(Atm(n)%gridstruct%grid, 1), lbound(Atm(n)%gridstruct%agrid, 2):ubound(Atm(n)%gridstruct%agrid, 2)))
    allocate(tile_geo_v%lats(lbound(Atm(n)%gridstruct%grid, 1):ubound(Atm(n)%gridstruct%grid, 1), lbound(Atm(n)%gridstruct%agrid, 2):ubound(Atm(n)%gridstruct%agrid, 2)))

    tile_geo_v%lons = -999.9
    tile_geo_v%lats = -999.9

    do x = lbound(tile_geo_v%lats, 1), ubound(tile_geo_v%lats, 1)
       do y = lbound(tile_geo_v%lats, 2), ubound(tile_geo_v%lats, 2)
          fp_i = (x - nest_x) * 2 + parent_x
          fp_j = (y - nest_y) * 2 + parent_y - 1

          tile_geo_v%lons(x,y) = fp_super_tile_geo%lons(fp_i, fp_j)
          tile_geo_v%lats(x,y) = fp_super_tile_geo%lats(fp_i, fp_j)
       end do
    end do

    if (debug_log) call show_tile_geo(tile_geo_v, this_pe, "tile_geo_v")


    !===========================================================
    !  End tile_geo per PE.
    !===========================================================

    allocate(p_grid(1:parent_geo%nxp, 1:parent_geo%nyp,2))
    allocate(n_grid(Atm(child_grid_num)%bd%isd:Atm(child_grid_num)%bd%ied, Atm(child_grid_num)%bd%jsd:Atm(child_grid_num)%bd%jed, 2))
    n_grid = real_snan


    allocate(p_grid_u(1:parent_geo%nxp, 1:parent_geo%nyp+1,2))
    allocate(n_grid_u(Atm(child_grid_num)%bd%isd:Atm(child_grid_num)%bd%ied, Atm(child_grid_num)%bd%jsd:Atm(child_grid_num)%bd%jed+1, 2))
    n_grid_u = real_snan


    allocate(p_grid_v(1:parent_geo%nxp+1, 1:parent_geo%nyp,2))
    allocate(n_grid_v(Atm(child_grid_num)%bd%isd:Atm(child_grid_num)%bd%ied+1, Atm(child_grid_num)%bd%jsd:Atm(child_grid_num)%bd%jed, 2))
    n_grid_v = real_snan

    ! TODO - propagate tile_geo information back to Atm structure
    ! TODO - deallocate tile_geo lat/lons
    ! TODO - ensure the allocation of tile_geo lat/lons is only performed once - outside the loop


    if (debug_log) print '("[INFO] WDR MV_NST2 run step 2 atmosphere.F90 npe=",I0, " tile_geo: nxp=",I0," nyp=",I0," nx=",I0," ny=", I0)', this_pe, tile_geo%nxp, tile_geo%nyp, tile_geo%nx, tile_geo%ny         
    if (debug_log) print *, "[INFO] WDR MV_NST2 run step 2 atmosphere.F90 shape(tile_geo%lats)=", shape(tile_geo%lats)
    if (debug_log) print '("[INFO] WDR MV_NST2 bounds1 (tile_geo%lats)=",I0,"-",I0)', lbound(tile_geo%lats,1), ubound(tile_geo%lats,1)
    if (debug_log) print '("[INFO] WDR MV_NST2 bounds2 (tile_geo%lats)=",I0,"-",I0)', lbound(tile_geo%lats,2), ubound(tile_geo%lats,2)

    call move_nest_geo(tile_geo, tile_geo_u, tile_geo_v, fp_super_tile_geo, delta_i_c, delta_j_c, x_refine, y_refine)


    call assign_n_p_grids(parent_geo, tile_geo, p_grid, n_grid, position)
    call assign_n_p_grids(parent_geo, tile_geo_u, p_grid_u, n_grid_u, position_u)
    call assign_n_p_grids(parent_geo, tile_geo_v, p_grid_v, n_grid_v, position_v)

  end subroutine mn_latlon_load_parent


  subroutine mn_latlon_read_hires_parent(npx, npy, x_refine, fp_super_tile_geo, surface_dir)
    integer, intent(in)                :: npx, npy, x_refine
    type(grid_geometry), intent(inout) :: fp_super_tile_geo
    character(len=120), intent(in)     :: surface_dir

    integer                            :: fp_super_istart_fine, fp_super_jstart_fine,fp_super_iend_fine, fp_super_jend_fine
    integer :: nx_cubic
    character(len=256) :: res_str

    nx_cubic = npx - 1
    write(res_str, '(I0)'), nx_cubic * x_refine

    call load_nest_latlons_from_nc(trim(surface_dir) // '/C' // trim(res_str) // '_grid.tile6.nc', &
         npx, npy, x_refine, &
         fp_super_tile_geo, &
         fp_super_istart_fine, fp_super_iend_fine, fp_super_jstart_fine, fp_super_jend_fine)


  end subroutine mn_latlon_read_hires_parent



  subroutine mn_orog_read_hires_parent(npx, npy, refine, surface_dir, filtered_terrain, orog_grid, orog_std_grid, ls_mask_grid, land_frac_grid)
    integer, intent(in)                :: npx, npy, refine
    character(len=120), intent(in)     :: surface_dir
    logical, intent(in)                :: filtered_terrain
    real, allocatable, intent(out)     :: orog_grid(:,:)
    real, allocatable, intent(out)     :: orog_std_grid(:,:)
    real, allocatable, intent(out)     :: ls_mask_grid(:,:)
    real, allocatable, intent(out)     :: land_frac_grid(:,:)

    integer :: nx_cubic, nx, ny, fp_nx, fp_ny, mid_nx, mid_ny
    integer :: fp_istart_fine, fp_iend_fine, fp_jstart_fine, fp_jend_fine


    character(len=256) :: res_str, parent_str
    character(len=512) :: nc_filename
    character(len=16)  :: orog_var_name

    integer :: parent_tile = 6 ! TODO: Later this will be configurable from namelist
    integer :: this_pe

    this_pe = mpp_pe()

    nx_cubic = npx - 1
    nx = npx - 1
    ny = npy - 1

    fp_istart_fine = 0
    fp_iend_fine = nx * refine
    fp_jstart_fine = 0
    fp_jend_fine = ny * refine

    fp_nx = fp_iend_fine - fp_istart_fine
    fp_ny = fp_jend_fine - fp_jstart_fine

    mid_nx = (fp_iend_fine - fp_istart_fine) / 2 
    mid_ny = (fp_jend_fine - fp_jstart_fine) / 2

    write(res_str, '(I0)'), nx_cubic * refine
    write(parent_str, '(I0)'), parent_tile
    nc_filename = trim(surface_dir) // '/C' // trim(res_str) // '_oro_data.tile' // trim(parent_str) // '.nc'

    if (filtered_terrain) then
       orog_var_name = 'orog_filt'    
    else
       orog_var_name = 'orog_raw'
    end if

    if (debug_log) print '("[INFO] WDR NCREAD LOFC mn_orog_read_hires_parent npe=",I0,I4,I4,I4,I4," ",A12," ",A128)', this_pe, fp_nx, fp_ny, mid_nx,mid_ny, orog_var_name, nc_filename

    call alloc_read_data(nc_filename, orog_var_name, fp_nx, fp_ny, orog_grid)
    call alloc_read_data(nc_filename, 'stddev', fp_nx, fp_ny, orog_std_grid)  ! Not needed

    call alloc_read_data(nc_filename, 'slmsk', fp_nx, fp_ny, ls_mask_grid)
    call alloc_read_data(nc_filename, 'land_frac', fp_nx, fp_ny, land_frac_grid)  ! Not needed

  end subroutine mn_orog_read_hires_parent



  subroutine mn_static_read_hires(npx, npy, refine, surface_dir, file_prefix, var_name, data_grid)
    integer, intent(in)                :: npx, npy, refine
    character(len=*), intent(in)     :: surface_dir, file_prefix
    character(len=*), intent(in)      :: var_name
    real, allocatable, intent(out)     :: data_grid(:,:)

    integer :: nx_cubic, nx, ny, fp_nx, fp_ny
    integer :: fp_istart_fine, fp_iend_fine, fp_jstart_fine, fp_jend_fine

    character(len=256) :: res_str, parent_str
    character(len=512) :: nc_filename

    integer :: parent_tile = 6 ! TODO: Later this will be configurable from namelist
    integer :: this_pe

    this_pe = mpp_pe()

    nx_cubic = npx - 1
    nx = npx - 1
    ny = npy - 1

    fp_istart_fine = 0
    fp_iend_fine = nx * refine
    fp_jstart_fine = 0
    fp_jend_fine = ny * refine

    fp_nx = fp_iend_fine - fp_istart_fine
    fp_ny = fp_jend_fine - fp_jstart_fine

    if (debug_log) print '("[INFO] WDR NCREAD LOFC mn_static_read_hires npe=",I0,I4,I4," ",A128," ",A128)', this_pe, fp_nx, fp_ny, var_name, nc_filename

    write(res_str, '(I0)'), nx_cubic * refine
    write(parent_str, '(I0)'), parent_tile
    
    if (trim(file_prefix) .eq. "oro_data") then
       nc_filename = trim(surface_dir) // '/C' // trim(res_str) // '_' // trim(file_prefix) // '.tile' // trim(parent_str) // '.nc'
    else
       nc_filename = trim(surface_dir) // '/C' // trim(res_str) // '.' // trim(file_prefix) // '.tile' // trim(parent_str) // '.nc'
    end if

    call alloc_read_data(nc_filename, var_name, fp_nx, fp_ny, data_grid)

  end subroutine mn_static_read_hires



  subroutine mn_reset_phys_latlon(Atm, n, tile_geo, fp_super_tile_geo, Atm_block, IPD_control, IPD_data)
    type(fv_atmos_type), intent(in)      :: Atm(:)
    integer, intent(in)                  :: n
    type(grid_geometry), intent(in)      :: tile_geo
    type(grid_geometry), intent(in)      :: fp_super_tile_geo
    type(block_control_type), intent(in) :: Atm_block
    type(IPD_control_type), intent(in) :: IPD_control
    type(IPD_data_type), intent(inout) :: IPD_data(:)
    !type(time_type), intent(in)     :: time_step
    
    integer :: isc, jsc, iec, jec
    integer :: x, y, fp_i, fp_j !, fpa_i, fpa_j
    integer :: nest_x, nest_y, parent_x, parent_y
    integer :: this_pe

    real(kind=kind_phys), allocatable :: lats(:,:), lons(:,:), area(:,:)


    this_pe = mpp_pe()

    isc = Atm(n)%bd%isc
    jsc = Atm(n)%bd%jsc
    iec = Atm(n)%bd%iec
    jec = Atm(n)%bd%jec

    allocate(lats(isc:iec, jsc:jec))
    allocate(lons(isc:iec, jsc:jec))
    allocate(area(isc:iec, jsc:jec))

    !print '("WDR mn_reset_phys_latlon AA npe=",I0)', this_pe

    ! This is going to be slow -- replace with better way to calculate index offsets, or pass them from earlier calculations
    call find_nest_alignment(tile_geo, fp_super_tile_geo, nest_x, nest_y, parent_x, parent_y)
    !print '("WDR mn_reset_phys_latlon AB npe=",I0)', this_pe

    do x = isc, iec
       do y = jsc, jec
          fp_i = (x - nest_x) * 2 + parent_x
          fp_j = (y - nest_y) * 2 + parent_y

          !fpa_i = (x - nest_x) * 2 + parent_x - 1
          !fpa_j = (y - nest_y) * 2 + parent_y - 1

          !print '("WDR mn_reset_phys_latlon BB npe=",I0," ix=",I0," x,y=",I0,",",I0)', this_pe, x,y
          lons(x,y) = fp_super_tile_geo%lons(fp_i, fp_j)
          !print '("WDR mn_reset_phys_latlon CC npe=",I0," ix=",I0," x,y=",I0,",",I0)', this_pe, x,y
          lats(x,y) = fp_super_tile_geo%lats(fp_i, fp_j)
          !print '("WDR mn_reset_phys_latlon DD npe=",I0," ix=",I0," x,y=",I0,",",I0)', this_pe, x,y

          ! Need to add the areas from 4 squares, because the netCDF file has areas calculated for the supergrid cells
          !  We need the area of the whole center of the cell.  
          !  Example dimensions for C288_grid.tile6.nc
          !   longitude -- x(577,577)
          !   latitude  -- y(577,577)
          !   area      -- x(576,576)


          !  Extracting lat/lon/area from Supergrid
          !
          !   1,1----2,1----3,1
          !    |      |      |
          !    | a1,1 | a2,1 |
          !    |      |      |
          !   1,2----2,2----3,2
          !    |      |      |
          !    | a1,2 | a2,2 |
          !    |      |      |
          !   1,3----2,3----3,3
          ! 
          !  The model A-grid cell 1,1 is centered at supergrid location 2,2
          !    The area of the A-grid cell is the sum of the 4 supergrid areas   A = a(1,1) + a(1,2) + a(2,1) + a(2,2)

          area(x,y) = fp_super_tile_geo%area(fp_i - 1, fp_j - 1) + fp_super_tile_geo%area(fp_i - 1, fp_j) + &
               fp_super_tile_geo%area(fp_i, fp_j - 1) + fp_super_tile_geo%area(fp_i, fp_j)   ! TODO make sure these offsets are correct.  
          !print '("WDR mn_reset_phys_latlon EE npe=",I0," ix=",I0," x,y=",I0,",",I0)', this_pe, x,y
       end do
    end do

    !print '("WDR mn_reset_phys_latlon FF npe=",I0," xlon_d=",F15.10," xlat_d=",F15.10," area=",F20.5)', this_pe, IPD_data(1)%Grid%xlon_d(1), IPD_data(1)%Grid%xlat_d(1), IPD_data(1)%Grid%area(1)

    call GFS_grid_populate(IPD_data%Grid, lons, lats, area)
    
    !print '("WDR mn_reset_phys_latlon GG npe=",I0," xlon_d=",F15.10," xlat_d=",F15.10," area=",F20.5)', this_pe, IPD_data(1)%Grid%xlon_d(1), IPD_data(1)%Grid%xlat_d(1), IPD_data(1)%Grid%area(1)

    deallocate(lats)
    deallocate(lons)
    deallocate(area)

  end subroutine mn_reset_phys_latlon



  !!============================================================================                                            
  !! Step 5.2 -- Recalculate nest halo weights
  !!============================================================================                                            

  subroutine mn_meta_recalc( delta_i_c, delta_j_c, x_refine, y_refine, tile_geo, parent_geo, fp_super_tile_geo, &
       is_fine_pe, nest_domain, position, p_grid, n_grid, wt, istart_coarse, jstart_coarse)

    integer, intent(in)                           :: delta_i_c, delta_j_c, x_refine, y_refine
    type(grid_geometry), intent(inout)            :: tile_geo, parent_geo, fp_super_tile_geo
    logical, intent(in)                           :: is_fine_pe
    type(nest_domain_type), intent(in)            :: nest_domain
    real(kind=R_GRID), allocatable, intent(inout) :: p_grid(:,:,:)
    real(kind=R_GRID), allocatable, intent(inout) :: n_grid(:,:,:)
    real, allocatable, intent(inout)              :: wt(:,:,:)
    integer, intent(inout)                        :: position
    integer, intent(in)                           :: istart_coarse, jstart_coarse

    type(bbox) :: wt_fine, wt_coarse
    integer    :: this_pe

    this_pe = mpp_pe()

    ! Update the coarse and fine indices after shifting the nest
    if (is_fine_pe) then

       if (debug_log) print '("[INFO] WDR NRD4 is_fine_pe=TRUE about to call bbox_get_C2F_index. npe=",I0, " position=",I0)', this_pe, position

       !!===========================================================
       !!
       !!  Recalculate halo weights
       !!
       !!===========================================================

       call bbox_get_C2F_index(nest_domain, wt_fine, wt_coarse, EAST,  position)
       call calc_nest_halo_weights(wt_fine, wt_coarse, p_grid, n_grid, wt, istart_coarse, jstart_coarse, x_refine, y_refine)

       call bbox_get_C2F_index(nest_domain, wt_fine, wt_coarse, WEST,  position)
       call calc_nest_halo_weights(wt_fine, wt_coarse, p_grid, n_grid, wt, istart_coarse, jstart_coarse, x_refine, y_refine)

       call bbox_get_C2F_index(nest_domain, wt_fine, wt_coarse, NORTH,  position)
       call calc_nest_halo_weights(wt_fine, wt_coarse, p_grid, n_grid, wt, istart_coarse, jstart_coarse, x_refine, y_refine)

       call bbox_get_C2F_index(nest_domain, wt_fine, wt_coarse, SOUTH,  position)
       call calc_nest_halo_weights(wt_fine, wt_coarse, p_grid, n_grid, wt, istart_coarse, jstart_coarse, x_refine, y_refine)

    end if

  end subroutine mn_meta_recalc




  !!============================================================================                                            
  !! Step 5.3 -- Adjust index by delta_i_c, delta_j_c
  !!============================================================================                                            

  subroutine mn_shift_index(delta_i_c, delta_j_c, ind)
    integer, intent(in)                    :: delta_i_c, delta_j_c
    integer, allocatable, intent(inout)    :: ind(:,:,:)

    ! Shift the index by the delta of this nest move.  
    ! TODO -- validate that we are not moving off the edge of the parent grid.
    integer  :: i, j

    do i = lbound(ind,1), ubound(ind,1)
       do j = lbound(ind,2), ubound(ind,2)
          ind(i,j,1) = ind(i,j,1) + delta_i_c
          ind(i,j,2) = ind(i,j,2) + delta_j_c
       end do
    end do

  end subroutine mn_shift_index



  !================================================================================ 
  !
  !  Prognostic and Physics Variable Nest Motion
  !
  !================================================================================ 

  !!============================================================================                                            
  !! Step 6   Shift the data on each nest PE                                                                                
  !!            -- similar to med_nest_move in HWRF                                                                         
  !!============================================================================      

  subroutine mn_prog_shift_data(Atm, n, child_grid_num, wt_h, wt_u, wt_v, &
       delta_i_c, delta_j_c, x_refine, y_refine, &
       is_fine_pe, nest_domain, nz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)
    integer, intent(in)                                      :: n, child_grid_num
    real, allocatable, intent(in)                            :: wt_h(:,:,:), wt_u(:,:,:), wt_v(:,:,:)
    integer, intent(in)                                      :: delta_i_c, delta_j_c, x_refine, y_refine   
    logical, intent(in)                                      :: is_fine_pe
    type(nest_domain_type), intent(inout)                    :: nest_domain
    integer, intent(in)                                      :: nz

    ! Constants for mpp calls
    integer  :: interp_type   = 1    ! cell-centered A-grid
    integer  :: interp_type_u = 4    ! D-grid
    integer  :: interp_type_v = 4    ! D-grid
    integer  :: position      = CENTER ! CENTER, NORTH, EAST
    integer  :: position_u    = NORTH
    integer  :: position_v    = EAST

    type(fv_moving_nest_prog_type), pointer :: mn_prog

    mn_prog => Atm(n)%mn_prog


    call mn_var_shift_data(Atm(n)%q_con, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
         delta_i_c, delta_j_c, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position, nz)

    call mn_var_shift_data(Atm(n)%pt, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
         delta_i_c, delta_j_c, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position, nz)

    call mn_var_shift_data(Atm(n)%w, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
         delta_i_c, delta_j_c, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position, nz)

    !call mn_var_shift_data(Atm(n)%omga, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
    !     delta_i_c, delta_j_c, &
    !     x_refine, y_refine, &
    !     is_fine_pe, nest_domain, position, nz)

    call mn_var_shift_data(Atm(n)%delp, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
         delta_i_c, delta_j_c, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position, nz)

    !call mn_var_shift_data(Atm(n)%delz, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
    !     delta_i_c, delta_j_c, &
    !     x_refine, y_refine, &
    !     is_fine_pe, nest_domain, position, nz)

    call mn_var_shift_data(mn_prog%delz, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
         delta_i_c, delta_j_c, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position, nz)


    call mn_var_shift_data(Atm(n)%ua, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
         delta_i_c, delta_j_c, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position, nz)

    call mn_var_shift_data(Atm(n)%va, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
         delta_i_c, delta_j_c, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position, nz)

    call mn_var_shift_data(Atm(n)%q, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
         delta_i_c, delta_j_c, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position, nz)

    !if (debug_log) print '("[INFO] WDR MV_NST6 show wt_u run step 6 atmosphere.F90 npe=",I0," n=",I0)', this_pe, n   
    !call check_array(Atm(n)%neststruct%wt_u, this_pe, "Atm(n)%neststruct%wt_u", 0.0, 1.0)
    !call check_array(wt_u, this_pe, "wt_u", 0.0, 1.0)
    !if (debug_log) print '("[INFO] WDR MV_NST6 stagger run step 6 atmosphere.F90 npe=",I0," n=",I0)', this_pe, n

    call mn_var_shift_data(Atm(n)%u, interp_type_u, wt_u, Atm(child_grid_num)%neststruct%ind_u, &
         delta_i_c, delta_j_c, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position_u, nz)

    call mn_var_shift_data(Atm(n)%v, interp_type_v, wt_v, Atm(child_grid_num)%neststruct%ind_v, &
         delta_i_c, delta_j_c, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position_v, nz)

  end subroutine mn_prog_shift_data

  subroutine mn_phys_shift_data(Atm, IPD_Control, IPD_Data, n, child_grid_num, wt_h, wt_u, wt_v, &
       delta_i_c, delta_j_c, x_refine, y_refine, &
       is_fine_pe, nest_domain, nz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)
    type(IPD_control_type), intent(in)               :: IPD_Control
    type(IPD_data_type), intent(inout)               :: IPD_Data(:)
    integer, intent(in)                              :: n, child_grid_num
    real, allocatable, intent(in)                    :: wt_h(:,:,:), wt_u(:,:,:), wt_v(:,:,:)
    integer, intent(in)                              :: delta_i_c, delta_j_c, x_refine, y_refine   
    logical, intent(in)                              :: is_fine_pe
    type(nest_domain_type), intent(inout)            :: nest_domain
    integer, intent(in)                              :: nz

    ! Constants for mpp calls
    integer  :: interp_type   = 1    ! cell-centered A-grid
    integer  :: interp_type_u = 4    ! D-grid
    integer  :: interp_type_v = 4    ! D-grid
    integer  :: position      = CENTER ! CENTER, NORTH, EAST
    integer  :: position_u    = NORTH
    integer  :: position_v    = EAST
    type(fv_moving_nest_physics_type), pointer :: mn_phys

    mn_phys => Atm(n)%mn_phys

    !! Skin temp/SST
    call mn_var_shift_data(mn_phys%ts, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
         delta_i_c, delta_j_c, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position)

    if (move_physics) then
       !! Soil variables
       call mn_var_shift_data(mn_phys%smc, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position, IPD_Control%lsoil)
       call mn_var_shift_data(mn_phys%stc, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position, IPD_Control%lsoil)
       call mn_var_shift_data(mn_phys%slc, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position, IPD_Control%lsoil)
       
       !! Physics arrays
       call mn_var_shift_data(mn_phys%phy_f2d, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position, IPD_control%ntot2d)
       
       call mn_var_shift_data(mn_phys%phy_f3d, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position, IPD_Control%levs)
       
       
       ! Surface variables
       call mn_var_shift_data(mn_phys%u10m, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
         delta_i_c, delta_j_c, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%v10m, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
         delta_i_c, delta_j_c, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%tprcp, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
         delta_i_c, delta_j_c, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position)

       call mn_var_shift_data(mn_phys%hprime, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
         delta_i_c, delta_j_c, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position, IPD_Control%nmtvr)
       

       call mn_var_shift_data(mn_phys%zorl, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
         delta_i_c, delta_j_c, &
         x_refine, y_refine, &
         is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%alvsf, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%alvwf, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%alnsf, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%alnwf, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)

       call mn_var_shift_data(mn_phys%facsf, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%facwf, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       
       
       call mn_var_shift_data(mn_phys%canopy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%vegfrac, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%uustar, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%shdmin, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%shdmax, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%zorll, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%zorlo, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%zorlw, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%tsfco, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%tsfcl, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%tsfc, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       
       call mn_var_shift_data(mn_phys%cv, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%cvt, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       call mn_var_shift_data(mn_phys%cvb, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
            delta_i_c, delta_j_c, &
            x_refine, y_refine, &
            is_fine_pe, nest_domain, position)
       
       end if
       
       if (move_nsst) then
          call mn_var_shift_data(mn_phys%tref, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
               delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
          call mn_var_shift_data(mn_phys%z_c, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
               delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
          call mn_var_shift_data(mn_phys%c_0, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
               delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
          call mn_var_shift_data(mn_phys%c_d, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
               delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
          call mn_var_shift_data(mn_phys%w_0, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
               delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
          call mn_var_shift_data(mn_phys%w_d, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
               delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
          call mn_var_shift_data(mn_phys%xt, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, & 
               delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
          call mn_var_shift_data(mn_phys%xs, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
               delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
          call mn_var_shift_data(mn_phys%xu, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
               delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
          call mn_var_shift_data(mn_phys%xv, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
               delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
          call mn_var_shift_data(mn_phys%xz, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
               delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
          call mn_var_shift_data(mn_phys%zm, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
               delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
          call mn_var_shift_data(mn_phys%xtts, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
               delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
          call mn_var_shift_data(mn_phys%xzts, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
               delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
          call mn_var_shift_data(mn_phys%d_conv, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
               delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
          !call mn_var_shift_data(mn_phys%ifd, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          !     delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
          call mn_var_shift_data(mn_phys%dt_cool, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
               delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
          call mn_var_shift_data(mn_phys%qrain, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
               delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)          
       end if

  end subroutine mn_phys_shift_data




  !!============================================================================                                            
  !! Step 6 - per variable
  !!============================================================================                                            

  !  TODO do we need a version of this for integer data?  land sea

  subroutine mn_var_shift_data2D(data_var, interp_type, wt, ind, delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)

    real, allocatable, intent(inout)            :: data_var(:,:)
    integer, intent(in)                         :: interp_type
    real, allocatable, intent(in)               :: wt(:,:,:)
    integer, allocatable, intent(in)            :: ind(:,:,:)
    integer, intent(in)                         :: delta_i_c, delta_j_c, x_refine, y_refine
    logical, intent(in)                         :: is_fine_pe
    type(nest_domain_type), intent(inout)       :: nest_domain
    integer, intent(in)                         :: position


    real, dimension(:,:), allocatable :: nbuffer, sbuffer, ebuffer, wbuffer
    logical         :: parent_proc, child_proc
    type(bbox)      :: north_fine, north_coarse ! step 4
    type(bbox)      :: south_fine, south_coarse
    type(bbox)      :: east_fine, east_coarse
    type(bbox)      :: west_fine, west_coarse
    integer         :: my_stat
    character(256)  :: my_errmsg
    integer         :: is, ie, js, je
    integer         :: this_pe

    integer         :: nest_level = 1  ! WDR TODO allow to vary


    this_pe = mpp_pe()

    !!===========================================================
    !!
    !! Fill halo buffers
    !!
    !!===========================================================

    if (debug_log) print '("[INFO] WDR NRD5. npe=",I0)', this_pe


    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%tile_fine=",I0," %tile_coarse=",I0)', this_pe, nest_domain%tile_fine, nest_domain%tile_coarse

    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%istart_fine=",I0," %iend_fine=",I0)', this_pe, nest_domain%istart_fine,  nest_domain%iend_fine
    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%jstart_fine=",I0," %jend_fine=",I0)', this_pe, nest_domain%jstart_fine,  nest_domain%jend_fine

    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%istart_coarse=",I0," %iend_coarse=",I0)', this_pe, nest_domain%istart_coarse,  nest_domain%iend_coarse
    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%jstart_coarse=",I0," %jend_coarse=",I0)', this_pe, nest_domain%jstart_coarse,  nest_domain%jend_coarse

    if (debug_log) print '("[INFO] data_var npe=",I0," data_var(",I0,"-",I0,",",I0,"-",I0,")")', this_pe,  lbound(data_var, 1), ubound(data_var, 1), lbound(data_var, 2), ubound(data_var, 2)

    if (debug_log) print '("[INFO] wt npe=",I0," wt(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', this_pe,  lbound(wt, 1), ubound(wt, 1), lbound(wt, 2), ubound(wt, 2), lbound(wt, 3), ubound(wt, 3)

    !====================================================
    if (debug_log) print '("[INFO] WDR ALL1. npe=",I0," position=",I0)', this_pe, position
    call alloc_halo_buffer(nbuffer, north_fine, north_coarse, nest_domain, NORTH,  position)
    call alloc_halo_buffer(sbuffer, south_fine, south_coarse, nest_domain, SOUTH,  position)
    call alloc_halo_buffer(ebuffer, east_fine,  east_coarse,  nest_domain, EAST,   position)
    call alloc_halo_buffer(wbuffer, west_fine,  west_coarse,  nest_domain, WEST,   position)

    if (debug_log) print '("[INFO] WDR allocate_halo_buffers DONE. npe=",I0)', this_pe

    !====================================================

    if (debug_log) print '("[INFO] WDR NRF0.d mn_var_shift_data npe=",I0," data_var(",I0,",",I0,")")', this_pe, size(data_var,1), size(data_var,2)

    if (debug_log) print '("[INFO] WDR NRF1 mn_var_shift_data start. npe=",I0)', this_pe

    ! Passes data from coarse grid to fine grid's halo buffers; requires nest_domain to be intent(inout)
    call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

    if (debug_log) print '("[INFO] WDR NRF2 mn_var_shift_data start. npe=",I0)', this_pe

    if (is_fine_pe) then
       if (debug_log) print '("[INFO] WDR NRF3 mn_var_shift_data start. npe=",I0)', this_pe

       !!===========================================================
       !!
       !! Shift grids internal to each nest PE
       !!
       !!===========================================================

       if ( delta_i_c .ne. 0 ) then
          if (debug_log) print '("[INFO] WDR NREX mn_var_shift_data start. npe=",I0)', this_pe
          data_var = eoshift(data_var, x_refine * delta_i_c, 0.0, 1)
       end if

       if (delta_j_c .ne.  0) then
          if (debug_log) print '("[INFO] WDR NREY mn_var_shift_data start. npe=",I0)', this_pe
          data_var = eoshift(data_var, y_refine * delta_j_c, 0.0, 2)
       end if

       !!===========================================================
       !!
       !! Apply halo data
       !!
       !!===========================================================

       if (debug_log) print '("[INFO] WDR NRFI mn_var_shift_data start. npe=",I0)', this_pe

       call fill_nest_from_buffer(interp_type, data_var, nbuffer, north_fine, north_coarse, NORTH, x_refine, y_refine, wt, ind)
       if (debug_log) print '("[INFO] WDR NRF N mn_var_shift_data start. npe=",I0)', this_pe

       call fill_nest_from_buffer(interp_type, data_var, sbuffer, south_fine, south_coarse, SOUTH, x_refine, y_refine, wt, ind)
       if (debug_log) print '("[INFO] WDR NRF S mn_var_shift_data start. npe=",I0)', this_pe

       call fill_nest_from_buffer(interp_type, data_var, ebuffer, east_fine, east_coarse, EAST, x_refine, y_refine, wt, ind)
       if (debug_log) print '("[INFO] WDR NRF E mn_var_shift_data start. npe=",I0)', this_pe

       call fill_nest_from_buffer(interp_type, data_var, wbuffer, west_fine, west_coarse, WEST, x_refine, y_refine, wt, ind)
       if (debug_log) print '("[INFO] WDR NRF W mn_var_shift_data start. npe=",I0)', this_pe

    end if

    deallocate(nbuffer)
    deallocate(sbuffer)
    deallocate(ebuffer)
    deallocate(wbuffer)

  end subroutine mn_var_shift_data2D


  subroutine mn_var_shift_data2D_kindphys(data_var, interp_type, wt, ind, delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)

    real(kind=kind_phys), allocatable, intent(inout)            :: data_var(:,:)
    integer, intent(in)                         :: interp_type
    real, allocatable, intent(in)               :: wt(:,:,:)
    integer, allocatable, intent(in)            :: ind(:,:,:)
    integer, intent(in)                         :: delta_i_c, delta_j_c, x_refine, y_refine
    logical, intent(in)                         :: is_fine_pe
    type(nest_domain_type), intent(inout)       :: nest_domain
    integer, intent(in)                         :: position


    real(kind=kind_phys), dimension(:,:), allocatable :: nbuffer, sbuffer, ebuffer, wbuffer
    logical         :: parent_proc, child_proc
    type(bbox)      :: north_fine, north_coarse ! step 4
    type(bbox)      :: south_fine, south_coarse
    type(bbox)      :: east_fine, east_coarse
    type(bbox)      :: west_fine, west_coarse
    integer         :: my_stat
    character(256)  :: my_errmsg
    integer         :: is, ie, js, je
    integer         :: this_pe

    integer         :: nest_level = 1  ! WDR TODO allow to vary


    this_pe = mpp_pe()

    !!===========================================================
    !!
    !! Fill halo buffers
    !!
    !!===========================================================

    if (debug_log) print '("[INFO] WDR NRD5. npe=",I0)', this_pe


    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%tile_fine=",I0," %tile_coarse=",I0)', this_pe, nest_domain%tile_fine, nest_domain%tile_coarse

    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%istart_fine=",I0," %iend_fine=",I0)', this_pe, nest_domain%istart_fine,  nest_domain%iend_fine
    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%jstart_fine=",I0," %jend_fine=",I0)', this_pe, nest_domain%jstart_fine,  nest_domain%jend_fine

    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%istart_coarse=",I0," %iend_coarse=",I0)', this_pe, nest_domain%istart_coarse,  nest_domain%iend_coarse
    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%jstart_coarse=",I0," %jend_coarse=",I0)', this_pe, nest_domain%jstart_coarse,  nest_domain%jend_coarse

    if (debug_log) print '("[INFO] data_var npe=",I0," data_var(",I0,"-",I0,",",I0,"-",I0,")")', this_pe,  lbound(data_var, 1), ubound(data_var, 1), lbound(data_var, 2), ubound(data_var, 2)

    if (debug_log) print '("[INFO] wt npe=",I0," wt(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', this_pe,  lbound(wt, 1), ubound(wt, 1), lbound(wt, 2), ubound(wt, 2), lbound(wt, 3), ubound(wt, 3)

    !====================================================
    if (debug_log) print '("[INFO] WDR ALL1. npe=",I0," position=",I0)', this_pe, position
    call alloc_halo_buffer(nbuffer, north_fine, north_coarse, nest_domain, NORTH,  position)
    call alloc_halo_buffer(sbuffer, south_fine, south_coarse, nest_domain, SOUTH,  position)
    call alloc_halo_buffer(ebuffer, east_fine,  east_coarse,  nest_domain, EAST,   position)
    call alloc_halo_buffer(wbuffer, west_fine,  west_coarse,  nest_domain, WEST,   position)

    if (debug_log) print '("[INFO] WDR allocate_halo_buffers DONE. npe=",I0)', this_pe

    !====================================================

    if (debug_log) print '("[INFO] WDR NRF0.d mn_var_shift_data npe=",I0," data_var(",I0,",",I0,")")', this_pe, size(data_var,1), size(data_var,2)

    if (debug_log) print '("[INFO] WDR NRF1 mn_var_shift_data start. npe=",I0)', this_pe

    ! Passes data from coarse grid to fine grid's halo buffers; requires nest_domain to be intent(inout)
    call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

    if (debug_log) print '("[INFO] WDR NRF2 mn_var_shift_data start. npe=",I0)', this_pe

    if (is_fine_pe) then
       if (debug_log) print '("[INFO] WDR NRF3 mn_var_shift_data start. npe=",I0)', this_pe

       !!===========================================================
       !!
       !! Shift grids internal to each nest PE
       !!
       !!===========================================================

       if ( delta_i_c .ne. 0 ) then
          if (debug_log) print '("[INFO] WDR NREX mn_var_shift_data start. npe=",I0)', this_pe
          data_var = eoshift(data_var, x_refine * delta_i_c, 0.0, 1)
       end if

       if (delta_j_c .ne.  0) then
          if (debug_log) print '("[INFO] WDR NREY mn_var_shift_data start. npe=",I0)', this_pe
          data_var = eoshift(data_var, y_refine * delta_j_c, 0.0, 2)
       end if

       !!===========================================================
       !!
       !! Apply halo data
       !!
       !!===========================================================

       if (debug_log) print '("[INFO] WDR NRFI mn_var_shift_data start. npe=",I0)', this_pe

       call fill_nest_from_buffer(interp_type, data_var, nbuffer, north_fine, north_coarse, NORTH, x_refine, y_refine, wt, ind)
       if (debug_log) print '("[INFO] WDR NRF N mn_var_shift_data start. npe=",I0)', this_pe

       call fill_nest_from_buffer(interp_type, data_var, sbuffer, south_fine, south_coarse, SOUTH, x_refine, y_refine, wt, ind)
       if (debug_log) print '("[INFO] WDR NRF S mn_var_shift_data start. npe=",I0)', this_pe

       call fill_nest_from_buffer(interp_type, data_var, ebuffer, east_fine, east_coarse, EAST, x_refine, y_refine, wt, ind)
       if (debug_log) print '("[INFO] WDR NRF E mn_var_shift_data start. npe=",I0)', this_pe

       call fill_nest_from_buffer(interp_type, data_var, wbuffer, west_fine, west_coarse, WEST, x_refine, y_refine, wt, ind)
       if (debug_log) print '("[INFO] WDR NRF W mn_var_shift_data start. npe=",I0)', this_pe

    end if

    deallocate(nbuffer)
    deallocate(sbuffer)
    deallocate(ebuffer)
    deallocate(wbuffer)

  end subroutine mn_var_shift_data2D_kindphys

  subroutine mn_var_shift_data3D(data_var, interp_type, wt, ind, delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)

    real, allocatable, intent(inout)            :: data_var(:,:,:)
    integer, intent(in)                         :: interp_type
    real, allocatable, intent(in)               :: wt(:,:,:)
    integer, allocatable, intent(in)            :: ind(:,:,:)
    integer, intent(in)                         :: delta_i_c, delta_j_c, x_refine, y_refine
    logical, intent(in)                         :: is_fine_pe
    type(nest_domain_type), intent(inout)       :: nest_domain
    integer, intent(in)                         :: position, nz


    real, dimension(:,:,:), allocatable :: nbuffer, sbuffer, ebuffer, wbuffer
    logical         :: parent_proc, child_proc
    type(bbox)      :: north_fine, north_coarse ! step 4
    type(bbox)      :: south_fine, south_coarse
    type(bbox)      :: east_fine, east_coarse
    type(bbox)      :: west_fine, west_coarse
    integer         :: my_stat
    character(256)  :: my_errmsg
    integer         :: is, ie, js, je
    integer         :: this_pe

    integer         :: nest_level = 1  ! WDR TODO allow to vary


    this_pe = mpp_pe()

    !!===========================================================
    !!
    !! Fill halo buffers
    !!
    !!===========================================================

    if (debug_log) print '("[INFO] WDR NRD5. npe=",I0)', this_pe


    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%tile_fine=",I0," %tile_coarse=",I0)', this_pe, nest_domain%tile_fine, nest_domain%tile_coarse

    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%istart_fine=",I0," %iend_fine=",I0)', this_pe, nest_domain%istart_fine,  nest_domain%iend_fine
    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%jstart_fine=",I0," %jend_fine=",I0)', this_pe, nest_domain%jstart_fine,  nest_domain%jend_fine

    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%istart_coarse=",I0," %iend_coarse=",I0)', this_pe, nest_domain%istart_coarse,  nest_domain%iend_coarse
    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%jstart_coarse=",I0," %jend_coarse=",I0)', this_pe, nest_domain%jstart_coarse,  nest_domain%jend_coarse

    if (debug_log) print '("[INFO] data_var npe=",I0," data_var(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', this_pe,  lbound(data_var, 1), ubound(data_var, 1), lbound(data_var, 2), ubound(data_var, 2), lbound(data_var, 3), ubound(data_var, 3)

    if (debug_log) print '("[INFO] wt npe=",I0," wt(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', this_pe,  lbound(wt, 1), ubound(wt, 1), lbound(wt, 2), ubound(wt, 2), lbound(wt, 3), ubound(wt, 3)

    !====================================================
    if (debug_log) print '("[INFO] WDR ALL1. npe=",I0," position=",I0," nz=",I0)', this_pe, position, nz
    call alloc_halo_buffer(nbuffer, north_fine, north_coarse, nest_domain, NORTH,  position, nz)
    call alloc_halo_buffer(sbuffer, south_fine, south_coarse, nest_domain, SOUTH,  position, nz)
    call alloc_halo_buffer(ebuffer, east_fine,  east_coarse,  nest_domain, EAST,   position, nz)
    call alloc_halo_buffer(wbuffer, west_fine,  west_coarse,  nest_domain, WEST,   position, nz)

    if (debug_log) print '("[INFO] WDR allocate_halo_buffers DONE. npe=",I0)', this_pe

    !====================================================

    if (debug_log) print '("[INFO] WDR NRF0.d mn_var_shift_data npe=",I0," data_var(",I0,",",I0,",",I0,")")', this_pe, size(data_var,1), size(data_var,2), size(data_var,3)

    if (debug_log) print '("[INFO] WDR NRF1 mn_var_shift_data start. npe=",I0)', this_pe

    ! Passes data from coarse grid to fine grid's halo buffers; requires nest_domain to be intent(inout)
    call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

    if (debug_log) print '("[INFO] WDR NRF2 mn_var_shift_data start. npe=",I0)', this_pe

    if (is_fine_pe) then
       if (debug_log) print '("[INFO] WDR NRF3 mn_var_shift_data start. npe=",I0)', this_pe

       !!===========================================================
       !!
       !! Shift grids internal to each nest PE
       !!
       !!===========================================================

       if ( delta_i_c .ne. 0 ) then
          if (debug_log) print '("[INFO] WDR NREX mn_var_shift_data start. npe=",I0)', this_pe
          data_var = eoshift(data_var, x_refine * delta_i_c, 0.0, 1)
       end if

       if (delta_j_c .ne.  0) then
          if (debug_log) print '("[INFO] WDR NREY mn_var_shift_data start. npe=",I0)', this_pe
          data_var = eoshift(data_var, y_refine * delta_j_c, 0.0, 2)
       end if

       !!===========================================================
       !!
       !! Apply halo data
       !!
       !!===========================================================

       if (debug_log) print '("[INFO] WDR NRFI mn_var_shift_data start. npe=",I0)', this_pe

       call fill_nest_from_buffer(interp_type, data_var, nbuffer, north_fine, north_coarse, nz, NORTH, x_refine, y_refine, wt, ind)
       if (debug_log) print '("[INFO] WDR NRF N mn_var_shift_data start. npe=",I0)', this_pe

       call fill_nest_from_buffer(interp_type, data_var, sbuffer, south_fine, south_coarse, nz, SOUTH, x_refine, y_refine, wt, ind)
       if (debug_log) print '("[INFO] WDR NRF S mn_var_shift_data start. npe=",I0)', this_pe

       call fill_nest_from_buffer(interp_type, data_var, ebuffer, east_fine, east_coarse, nz, EAST, x_refine, y_refine, wt, ind)
       if (debug_log) print '("[INFO] WDR NRF E mn_var_shift_data start. npe=",I0)', this_pe

       call fill_nest_from_buffer(interp_type, data_var, wbuffer, west_fine, west_coarse, nz, WEST, x_refine, y_refine, wt, ind)
       if (debug_log) print '("[INFO] WDR NRF W mn_var_shift_data start. npe=",I0)', this_pe

    end if

    deallocate(nbuffer)
    deallocate(sbuffer)
    deallocate(ebuffer)
    deallocate(wbuffer)

  end subroutine mn_var_shift_data3D


  subroutine mn_var_shift_data3D_kindphys(data_var, interp_type, wt, ind, delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)

    real(kind=kind_phys), allocatable, intent(inout)            :: data_var(:,:,:)
    integer, intent(in)                         :: interp_type
    real, allocatable, intent(in)               :: wt(:,:,:)
    integer, allocatable, intent(in)            :: ind(:,:,:)
    integer, intent(in)                         :: delta_i_c, delta_j_c, x_refine, y_refine
    logical, intent(in)                         :: is_fine_pe
    type(nest_domain_type), intent(inout)       :: nest_domain
    integer, intent(in)                         :: position, nz


    real(kind=kind_phys), dimension(:,:,:), allocatable :: nbuffer, sbuffer, ebuffer, wbuffer
    logical         :: parent_proc, child_proc
    type(bbox)      :: north_fine, north_coarse ! step 4
    type(bbox)      :: south_fine, south_coarse
    type(bbox)      :: east_fine, east_coarse
    type(bbox)      :: west_fine, west_coarse
    integer         :: my_stat
    character(256)  :: my_errmsg
    integer         :: is, ie, js, je
    integer         :: this_pe

    integer         :: nest_level = 1  ! WDR TODO allow to vary


    this_pe = mpp_pe()

    !!===========================================================
    !!
    !! Fill halo buffers
    !!
    !!===========================================================

    if (debug_log) print '("[INFO] WDR NRD5. npe=",I0)', this_pe


    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%tile_fine=",I0," %tile_coarse=",I0)', this_pe, nest_domain%tile_fine, nest_domain%tile_coarse

    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%istart_fine=",I0," %iend_fine=",I0)', this_pe, nest_domain%istart_fine,  nest_domain%iend_fine
    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%jstart_fine=",I0," %jend_fine=",I0)', this_pe, nest_domain%jstart_fine,  nest_domain%jend_fine

    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%istart_coarse=",I0," %iend_coarse=",I0)', this_pe, nest_domain%istart_coarse,  nest_domain%iend_coarse
    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%jstart_coarse=",I0," %jend_coarse=",I0)', this_pe, nest_domain%jstart_coarse,  nest_domain%jend_coarse

    if (debug_log) print '("[INFO] data_var npe=",I0," data_var(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', this_pe,  lbound(data_var, 1), ubound(data_var, 1), lbound(data_var, 2), ubound(data_var, 2), lbound(data_var, 3), ubound(data_var, 3)

    if (debug_log) print '("[INFO] wt npe=",I0," wt(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', this_pe,  lbound(wt, 1), ubound(wt, 1), lbound(wt, 2), ubound(wt, 2), lbound(wt, 3), ubound(wt, 3)

    !====================================================
    if (debug_log) print '("[INFO] WDR ALL1. npe=",I0," position=",I0," nz=",I0)', this_pe, position, nz
    call alloc_halo_buffer(nbuffer, north_fine, north_coarse, nest_domain, NORTH,  position, nz)
    call alloc_halo_buffer(sbuffer, south_fine, south_coarse, nest_domain, SOUTH,  position, nz)
    call alloc_halo_buffer(ebuffer, east_fine,  east_coarse,  nest_domain, EAST,   position, nz)
    call alloc_halo_buffer(wbuffer, west_fine,  west_coarse,  nest_domain, WEST,   position, nz)

    if (debug_log) print '("[INFO] WDR allocate_halo_buffers DONE. npe=",I0)', this_pe

    !====================================================

    if (debug_log) print '("[INFO] WDR NRF0.d mn_var_shift_data npe=",I0," data_var(",I0,",",I0,",",I0,")")', this_pe, size(data_var,1), size(data_var,2), size(data_var,3)

    if (debug_log) print '("[INFO] WDR NRF1 mn_var_shift_data start. npe=",I0)', this_pe

    ! Passes data from coarse grid to fine grid's halo buffers; requires nest_domain to be intent(inout)
    call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

    if (debug_log) print '("[INFO] WDR NRF2 mn_var_shift_data start. npe=",I0)', this_pe

    if (is_fine_pe) then
       if (debug_log) print '("[INFO] WDR NRF3 mn_var_shift_data start. npe=",I0)', this_pe

       !!===========================================================
       !!
       !! Shift grids internal to each nest PE
       !!
       !!===========================================================

       if ( delta_i_c .ne. 0 ) then
          if (debug_log) print '("[INFO] WDR NREX mn_var_shift_data start. npe=",I0)', this_pe
          data_var = eoshift(data_var, x_refine * delta_i_c, 0.0, 1)
       end if

       if (delta_j_c .ne.  0) then
          if (debug_log) print '("[INFO] WDR NREY mn_var_shift_data start. npe=",I0)', this_pe
          data_var = eoshift(data_var, y_refine * delta_j_c, 0.0, 2)
       end if

       !!===========================================================
       !!
       !! Apply halo data
       !!
       !!===========================================================

       if (debug_log) print '("[INFO] WDR NRFI mn_var_shift_data start. npe=",I0)', this_pe

       call fill_nest_from_buffer(interp_type, data_var, nbuffer, north_fine, north_coarse, nz, NORTH, x_refine, y_refine, wt, ind)
       if (debug_log) print '("[INFO] WDR NRF N mn_var_shift_data start. npe=",I0)', this_pe

       call fill_nest_from_buffer(interp_type, data_var, sbuffer, south_fine, south_coarse, nz, SOUTH, x_refine, y_refine, wt, ind)
       if (debug_log) print '("[INFO] WDR NRF S mn_var_shift_data start. npe=",I0)', this_pe

       call fill_nest_from_buffer(interp_type, data_var, ebuffer, east_fine, east_coarse, nz, EAST, x_refine, y_refine, wt, ind)
       if (debug_log) print '("[INFO] WDR NRF E mn_var_shift_data start. npe=",I0)', this_pe

       call fill_nest_from_buffer(interp_type, data_var, wbuffer, west_fine, west_coarse, nz, WEST, x_refine, y_refine, wt, ind)
       if (debug_log) print '("[INFO] WDR NRF W mn_var_shift_data start. npe=",I0)', this_pe

    end if

    deallocate(nbuffer)
    deallocate(sbuffer)
    deallocate(ebuffer)
    deallocate(wbuffer)

  end subroutine mn_var_shift_data3D_kindphys

  subroutine mn_var_shift_data4D(data_var, interp_type, wt, ind, delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)

    real, allocatable, intent(inout)            :: data_var(:,:,:,:)
    integer, intent(in)                         :: interp_type
    real, allocatable, intent(in)               :: wt(:,:,:)
    integer, allocatable, intent(in)            :: ind(:,:,:)
    integer, intent(in)                         :: delta_i_c, delta_j_c, x_refine, y_refine
    logical, intent(in)                         :: is_fine_pe
    type(nest_domain_type), intent(inout)       :: nest_domain
    integer, intent(in)                         :: position, nz

    real, dimension(:,:,:,:), allocatable         :: nbuffer, sbuffer, ebuffer, wbuffer
    logical         :: parent_proc, child_proc
    type(bbox)      :: north_fine, north_coarse ! step 4
    type(bbox)      :: south_fine, south_coarse
    type(bbox)      :: east_fine, east_coarse
    type(bbox)      :: west_fine, west_coarse
    integer         :: my_stat
    character(256)  :: my_errmsg
    integer         :: n4d
    integer         :: this_pe
    integer         :: is, ie, js, je
    integer         :: nest_level = 1  ! WDR TODO allow to vary


    this_pe = mpp_pe()

    n4d = ubound(data_var, 4)

    !!===========================================================
    !!
    !! Fill halo buffers
    !!
    !!===========================================================

    if (debug_log) print '("[INFO] WDR NRD5. npe=",I0)', this_pe


    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%tile_fine=",I0," %tile_coarse=",I0)', this_pe, nest_domain%tile_fine, nest_domain%tile_coarse

    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%istart_fine=",I0," %iend_fine=",I0)', this_pe, nest_domain%istart_fine,  nest_domain%iend_fine
    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%jstart_fine=",I0," %jend_fine=",I0)', this_pe, nest_domain%jstart_fine,  nest_domain%jend_fine

    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%istart_coarse=",I0," %iend_coarse=",I0)', this_pe, nest_domain%istart_coarse,  nest_domain%iend_coarse
    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%jstart_coarse=",I0," %jend_coarse=",I0)', this_pe, nest_domain%jstart_coarse,  nest_domain%jend_coarse

    if (debug_log) print '("[INFO] data_var npe=",I0," data_var(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', this_pe,  lbound(data_var, 1), ubound(data_var, 1), lbound(data_var, 2), ubound(data_var, 2), lbound(data_var, 3), ubound(data_var, 3)

    if (debug_log) print '("[INFO] wt npe=",I0," wt(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', this_pe,  lbound(wt, 1), ubound(wt, 1), lbound(wt, 2), ubound(wt, 2), lbound(wt, 3), ubound(wt, 3)

    !====================================================
    if (debug_log) print '("[INFO] WDR ALL1. npe=",I0," position=",I0," nz=",I0)', this_pe, position, nz

    call alloc_halo_buffer(nbuffer, north_fine, north_coarse, nest_domain, NORTH,  position, nz, n4d)
    call alloc_halo_buffer(sbuffer, south_fine, south_coarse, nest_domain, SOUTH,  position, nz, n4d)
    call alloc_halo_buffer(ebuffer, east_fine,  east_coarse,  nest_domain, EAST,   position, nz, n4d)
    call alloc_halo_buffer(wbuffer, west_fine,  west_coarse,  nest_domain, WEST,   position, nz, n4d)

    if (debug_log) print '("[INFO] WDR allocate_halo_buffers DONE. npe=",I0)', this_pe


    !====================================================

    !call mpp_sync(full_pelist)

    if (debug_log) print '("[INFO] WDR NRF0.d mn_var_shift_data npe=",I0," data_var(",I0,",",I0,",",I0,",",",I0,")")', this_pe, size(data_var,1), size(data_var,2), size(data_var,3), size(data_var,4)
    if (debug_log) print '("[INFO] WDR NRF1 mn_var_shift_data start. npe=",I0)', this_pe

    ! Passes data from coarse grid to fine grid's halo
    call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

    if (debug_log) print '("[INFO] WDR NRF2 mn_var_shift_data start. npe=",I0)', this_pe

    if (is_fine_pe) then
       if (debug_log) print '("[INFO] WDR NRF3 mn_var_shift_data start. npe=",I0)', this_pe

       !!===========================================================
       !!
       !! Shift grids internal to each nest PE
       !!
       !!===========================================================

       if ( delta_i_c .ne. 0 ) then
          if (debug_log) print '("[INFO] WDR NREX mn_var_shift_data start. npe=",I0)', this_pe
          data_var = eoshift(data_var, x_refine * delta_i_c, 0.0, 1)
       end if

       if (delta_j_c .ne.  0) then
          if (debug_log) print '("[INFO] WDR NREY mn_var_shift_data start. npe=",I0)', this_pe
          data_var = eoshift(data_var, y_refine * delta_j_c, 0.0, 2)
       end if

       !call output_logical_grid("FF", isd_fine, ied_fine, jsd_fine, jed_fine, nz, x)

       !!===========================================================
       !!
       !! Apply halo data
       !!
       !!===========================================================

       if (debug_log) print '("[INFO] WDR NRFI mn_var_shift_data start. npe=",I0)', this_pe

       call fill_nest_from_buffer(interp_type, data_var, nbuffer, north_fine, north_coarse, nz, NORTH, x_refine, y_refine, wt, ind)
       if (debug_log) print '("[INFO] WDR NRF N mn_var_shift_data start. npe=",I0)', this_pe

       call fill_nest_from_buffer(interp_type, data_var, sbuffer, south_fine, south_coarse, nz, SOUTH, x_refine, y_refine, wt, ind)
       if (debug_log) print '("[INFO] WDR NRF S mn_var_shift_data start. npe=",I0)', this_pe

       call fill_nest_from_buffer(interp_type, data_var, ebuffer, east_fine, east_coarse, nz, EAST, x_refine, y_refine, wt, ind)
       if (debug_log) print '("[INFO] WDR NRF E mn_var_shift_data start. npe=",I0)', this_pe

       call fill_nest_from_buffer(interp_type, data_var, wbuffer, west_fine, west_coarse, nz, WEST, x_refine, y_refine, wt, ind)
       if (debug_log) print '("[INFO] WDR NRF W mn_var_shift_data start. npe=",I0)', this_pe

    end if

    deallocate(nbuffer)
    deallocate(sbuffer)
    deallocate(ebuffer)
    deallocate(wbuffer)

  end subroutine mn_var_shift_data4D


  subroutine mn_var_shift_data4D_kindphys(data_var, interp_type, wt, ind, delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, nz)

    real(kind=kind_phys), allocatable, intent(inout)            :: data_var(:,:,:,:)
    integer, intent(in)                         :: interp_type
    real, allocatable, intent(in)               :: wt(:,:,:)
    integer, allocatable, intent(in)            :: ind(:,:,:)
    integer, intent(in)                         :: delta_i_c, delta_j_c, x_refine, y_refine
    logical, intent(in)                         :: is_fine_pe
    type(nest_domain_type), intent(inout)       :: nest_domain
    integer, intent(in)                         :: position, nz

    real(kind=kind_phys), dimension(:,:,:,:), allocatable         :: nbuffer, sbuffer, ebuffer, wbuffer
    logical         :: parent_proc, child_proc
    type(bbox)      :: north_fine, north_coarse ! step 4
    type(bbox)      :: south_fine, south_coarse
    type(bbox)      :: east_fine, east_coarse
    type(bbox)      :: west_fine, west_coarse
    integer         :: my_stat
    character(256)  :: my_errmsg
    integer         :: n4d
    integer         :: this_pe
    integer         :: is, ie, js, je
    integer         :: nest_level = 1  ! WDR TODO allow to vary


    this_pe = mpp_pe()

    n4d = ubound(data_var, 4)

    !!===========================================================
    !!
    !! Fill halo buffers
    !!
    !!===========================================================

    if (debug_log) print '("[INFO] WDR NRD5. npe=",I0)', this_pe


    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%tile_fine=",I0," %tile_coarse=",I0)', this_pe, nest_domain%tile_fine, nest_domain%tile_coarse

    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%istart_fine=",I0," %iend_fine=",I0)', this_pe, nest_domain%istart_fine,  nest_domain%iend_fine
    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%jstart_fine=",I0," %jend_fine=",I0)', this_pe, nest_domain%jstart_fine,  nest_domain%jend_fine

    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%istart_coarse=",I0," %iend_coarse=",I0)', this_pe, nest_domain%istart_coarse,  nest_domain%iend_coarse
    if (debug_log) print '("[INFO] show_nest_domain npe=",I0," nest_domain%jstart_coarse=",I0," %jend_coarse=",I0)', this_pe, nest_domain%jstart_coarse,  nest_domain%jend_coarse

    if (debug_log) print '("[INFO] data_var npe=",I0," data_var(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', this_pe,  lbound(data_var, 1), ubound(data_var, 1), lbound(data_var, 2), ubound(data_var, 2), lbound(data_var, 3), ubound(data_var, 3)

    if (debug_log) print '("[INFO] wt npe=",I0," wt(",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', this_pe,  lbound(wt, 1), ubound(wt, 1), lbound(wt, 2), ubound(wt, 2), lbound(wt, 3), ubound(wt, 3)

    !====================================================
    if (debug_log) print '("[INFO] WDR ALL1. npe=",I0," position=",I0," nz=",I0)', this_pe, position, nz

    call alloc_halo_buffer(nbuffer, north_fine, north_coarse, nest_domain, NORTH,  position, nz, n4d)
    call alloc_halo_buffer(sbuffer, south_fine, south_coarse, nest_domain, SOUTH,  position, nz, n4d)
    call alloc_halo_buffer(ebuffer, east_fine,  east_coarse,  nest_domain, EAST,   position, nz, n4d)
    call alloc_halo_buffer(wbuffer, west_fine,  west_coarse,  nest_domain, WEST,   position, nz, n4d)

    if (debug_log) print '("[INFO] WDR allocate_halo_buffers DONE. npe=",I0)', this_pe


    !====================================================

    !call mpp_sync(full_pelist)

    if (debug_log) print '("[INFO] WDR NRF0.d mn_var_shift_data npe=",I0," data_var(",I0,",",I0,",",I0,",",",I0,")")', this_pe, size(data_var,1), size(data_var,2), size(data_var,3), size(data_var,4)
    if (debug_log) print '("[INFO] WDR NRF1 mn_var_shift_data start. npe=",I0)', this_pe

    ! Passes data from coarse grid to fine grid's halo
    call mpp_update_nest_fine(data_var, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level, position=position)

    if (debug_log) print '("[INFO] WDR NRF2 mn_var_shift_data start. npe=",I0)', this_pe

    if (is_fine_pe) then
       if (debug_log) print '("[INFO] WDR NRF3 mn_var_shift_data start. npe=",I0)', this_pe

       !!===========================================================
       !!
       !! Shift grids internal to each nest PE
       !!
       !!===========================================================

       if ( delta_i_c .ne. 0 ) then
          if (debug_log) print '("[INFO] WDR NREX mn_var_shift_data start. npe=",I0)', this_pe
          data_var = eoshift(data_var, x_refine * delta_i_c, 0.0, 1)
       end if

       if (delta_j_c .ne.  0) then
          if (debug_log) print '("[INFO] WDR NREY mn_var_shift_data start. npe=",I0)', this_pe
          data_var = eoshift(data_var, y_refine * delta_j_c, 0.0, 2)
       end if

       !call output_logical_grid("FF", isd_fine, ied_fine, jsd_fine, jed_fine, nz, x)

       !!===========================================================
       !!
       !! Apply halo data
       !!
       !!===========================================================

       if (debug_log) print '("[INFO] WDR NRFI mn_var_shift_data start. npe=",I0)', this_pe

       call fill_nest_from_buffer(interp_type, data_var, nbuffer, north_fine, north_coarse, nz, NORTH, x_refine, y_refine, wt, ind)
       if (debug_log) print '("[INFO] WDR NRF N mn_var_shift_data start. npe=",I0)', this_pe

       call fill_nest_from_buffer(interp_type, data_var, sbuffer, south_fine, south_coarse, nz, SOUTH, x_refine, y_refine, wt, ind)
       if (debug_log) print '("[INFO] WDR NRF S mn_var_shift_data start. npe=",I0)', this_pe

       call fill_nest_from_buffer(interp_type, data_var, ebuffer, east_fine, east_coarse, nz, EAST, x_refine, y_refine, wt, ind)
       if (debug_log) print '("[INFO] WDR NRF E mn_var_shift_data start. npe=",I0)', this_pe

       call fill_nest_from_buffer(interp_type, data_var, wbuffer, west_fine, west_coarse, nz, WEST, x_refine, y_refine, wt, ind)
       if (debug_log) print '("[INFO] WDR NRF W mn_var_shift_data start. npe=",I0)', this_pe

    end if

    deallocate(nbuffer)
    deallocate(sbuffer)
    deallocate(ebuffer)
    deallocate(wbuffer)

  end subroutine mn_var_shift_data4D_kindphys



  !================================================================================ 
  !
  !  Step 7 -- Gridstruct resetting and reallocation of static buffers
  !      init_grid() also updates the wt arrays
  !================================================================================ 


  subroutine mn_meta_reset_gridstruct(Atm, n, child_grid_num, nest_domain, fp_super_tile_geo, x_refine, y_refine, is_fine_pe, wt_h, wt_u, wt_v, a_step, dt_atmos)
    type(fv_atmos_type), allocatable, intent(inout)  :: Atm(:)
    integer, intent(in)                              :: n, child_grid_num
    type(nest_domain_type),     intent(in)           :: nest_domain
    type(grid_geometry), intent(in)                  :: fp_super_tile_geo
    integer, intent(in)                              :: x_refine, y_refine   
    logical, intent(in)                              :: is_fine_pe
    real, allocatable, intent(in)                    :: wt_h(:,:,:), wt_u(:,:,:), wt_v(:,:,:)   
    integer, intent(in)                              :: a_step
    real, intent(in)                                 :: dt_atmos
    

    integer :: isg, ieg, jsg, jeg
    integer :: ng, pp, nn, parent_tile, refinement, ioffset, joffset
    integer :: this_pe, gid
    integer :: tile_coarse(2)
    integer :: half_x, half_y

    real(kind=R_GRID)   :: pi = 4 * atan(1.0d0)
    real                :: rad2deg, half_lat, half_lon

    logical, save       :: first_time = .true.
    integer, save       :: id_reset1, id_reset2, id_reset3, id_reset4, id_reset5, id_reset6, id_reset7


    if (first_time) then
       id_reset1     = mpp_clock_id ('MN 7 Reset 1',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
       id_reset2     = mpp_clock_id ('MN 7 Reset 2',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
       id_reset3     = mpp_clock_id ('MN 7 Reset 3',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
       id_reset4     = mpp_clock_id ('MN 7 Reset 4',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
       id_reset5     = mpp_clock_id ('MN 7 Reset 5',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
       id_reset6     = mpp_clock_id ('MN 7 Reset 6',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
       id_reset7     = mpp_clock_id ('MN 7 Reset 7',  flags = clock_flag_default, grain=CLOCK_ROUTINE )
    end if

    

    rad2deg = 180.0 / pi

    this_pe = mpp_pe()
    gid = this_pe

    parent_tile = Atm(child_grid_num)%neststruct%parent_tile
    ioffset = Atm(child_grid_num)%neststruct%ioffset
    joffset = Atm(child_grid_num)%neststruct%joffset

    ! Log the bounds of this PE's grid after nest motion.  TODO replace step 4 with timestep
    if (is_fine_pe .and. debug_log) then
       call show_nest_grid(Atm(n), this_pe, 4)
    end if

    !  Reset the gridstruct values for the nest
    if (is_fine_pe) then
       ! Fill in values from high resolution, full panel, supergrid 
       call mpp_clock_begin (id_reset1)

       call fill_grid_from_supergrid(Atm(n)%gridstruct%grid, CORNER, fp_super_tile_geo, ioffset, joffset, &
            x_refine, y_refine)
       call fill_grid_from_supergrid(Atm(n)%gridstruct%agrid, CENTER, fp_super_tile_geo, ioffset, joffset, &
            x_refine, y_refine)
       call fill_grid_from_supergrid(Atm(n)%gridstruct%grid_64, CORNER, fp_super_tile_geo, &
            ioffset, joffset, x_refine, y_refine)
       call fill_grid_from_supergrid(Atm(n)%gridstruct%agrid_64, CENTER, fp_super_tile_geo, &
            ioffset, joffset, x_refine, y_refine)

       ! What's the status of Atm(n)%grid_global?
       if (debug_log) print '("[INFO] WDR Atm(1) GLOBAL npe=",I0," grid_global(",I0,"-",I0",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', this_pe, &
            lbound(Atm(1)%grid_global,1), ubound(Atm(1)%grid_global,1), &
            lbound(Atm(1)%grid_global,2), ubound(Atm(1)%grid_global,2), &
            lbound(Atm(1)%grid_global,3), ubound(Atm(1)%grid_global,3), &
            lbound(Atm(1)%grid_global,4), ubound(Atm(1)%grid_global,4)


       if (debug_log) print '("[INFO] WDR Atm(n) GLOBAL npe=",I0," grid_global(",I0,"-",I0",",I0,"-",I0,",",I0,"-",I0,",",I0,"-",I0,")")', this_pe, &
            lbound(Atm(n)%grid_global,1), ubound(Atm(n)%grid_global,1), &
            lbound(Atm(n)%grid_global,2), ubound(Atm(n)%grid_global,2), &
            lbound(Atm(n)%grid_global,3), ubound(Atm(n)%grid_global,3), &
            lbound(Atm(n)%grid_global,4), ubound(Atm(n)%grid_global,4)

       !! Let this get reset in init_grid()/setup_aligned_nest()
       !call fill_grid_from_supergrid(Atm(n)%grid_global, CORNER, fp_super_tile_geo, &
       !     ioffset, joffset, x_refine, y_refine)

       call mpp_clock_end (id_reset1)
       call mpp_clock_begin (id_reset2)



       ! TODO should these get reset by init_grid instead??
       call fill_weight_grid(Atm(n)%neststruct%wt_h, wt_h)
       call fill_weight_grid(Atm(n)%neststruct%wt_u, wt_u)
       call fill_weight_grid(Atm(n)%neststruct%wt_v, wt_v)
       ! WDR TODO -- Seems like this is not used anywhere, other than being allocated, filled, deallocated
       !call fill_weight_grid(Atm(n)%neststruct%wt_b, wt_b)

       call mpp_clock_end (id_reset2)

    end if

    if (debug_log) print '("[INFO] WDR INIT_GRID AP1 fv_moving_nest.F90 npe=",I0," n=",I0)', this_pe, n

    call mpp_clock_begin (id_reset3)


    ! TODO Write clearer comments on what is happening here.

    ! This code runs several communications steps:
    !  1.  As npe=0, it gets the global_grid domain setup 
    !  2.  sends the global_grid to the other parent PEs
    !  3.  global_grid is received in call to setup_aligned_nest() in fv_grid_tools.F90::init_grid()
    !  Other communication is contained full within setup_aligned_nest().



    ! Sends around data from the parent grids, and recomputes the update indices
    ! This code copied from fv_control.F90
    ! Need to SEND grid_global to any child grids; this is received in setup_aligned_nest in fv_grid_tools
    ! if (Atm(pp)%neststruct%nested) then

    ! TODO phrase this more carefully to choose the parent master PE grid if we are operating in a nested setup. 
    ! Unlike in fv_control.F90, this will be running on Atm(1) when it's on pe=0, so we don't need to navigate to parent_grid.

    first_time = .false.

    ! Seems like we do not need to resend this -- setup_aligned_nest now saves the parent tile information during model initialization,
    !  which happens before we enter the moving nest code.  
    if (this_pe .eq. 0 .and. first_time) then

       ! This is the Atm index for the nest values.  
       pp = child_grid_num

       if (debug_log) print '("[INFO] WDR INIT_GRID AP2 fv_moving_nest.F90 npe=",I0," n=",I0," pp=",I0)', this_pe, n, pp

       refinement = x_refine
       ng = Atm(n)%ng

       call mpp_get_global_domain( Atm(n)%domain, isg, ieg, jsg, jeg)

       !if (debug_log) print '("[INFO] WDR INIT_GRID AP3.1 fv_moving_nest.F90 npe=",I0," gid=",I0," associated(parent_grid)=",L1)', this_pe, gid, associated(Atm(pp)%parent_grid)
       if (debug_log) print '("[INFO] WDR INIT_GRID AP3.1 fv_moving_nest.F90 npe=",I0," gid=",I0," parent_tile=",I0)', this_pe, gid, parent_tile
       if (debug_log) print '("[INFO] WDR INIT_GRID AP3.2 fv_moving_nest.F90 npe=",I0," gid=",I0," size(pelist)=",I0)', this_pe, gid, size(Atm(pp)%pelist)
       if (debug_log) print '("[INFO] WDR INIT_GRID AP3.3 fv_moving_nest.F90 npe=",I0," gid=",I0," pelist1=",I0)', this_pe, gid, Atm(pp)%pelist(1)

       !FIXME: Should replace this by generating the global grid (or at least one face thereof) on the 
       ! nested PEs instead of sending it around.
       !if (gid == Atm(pp)%parent_grid%pelist(1)) then
       if (debug_log) print '("[INFO] WDR INIT_GRID XFER AP4 fv_moving_nest.F90 npe=",I0," send to pe=",I0," size=",I0)', this_pe, Atm(pp)%pelist(1), size(Atm(n)%grid_global(isg-ng:ieg+1+ng,jsg-ng:jeg+1+ng,1:2,parent_tile))

       call mpp_send(Atm(n)%grid_global(isg-ng:ieg+1+ng,jsg-ng:jeg+1+ng,1:2,parent_tile), &
            size(Atm(n)%grid_global(isg-ng:ieg+1+ng,jsg-ng:jeg+1+ng,1:2,parent_tile)), &
            Atm(pp)%pelist(1)) !send to p_ind in setup_aligned_nest
       if (debug_log) print '("[INFO] WDR INIT_GRID AP5 fv_moving_nest.F90 npe=",I0)', this_pe
       call mpp_sync_self()
       if (debug_log) print '("[INFO] WDR INIT_GRID AP6 fv_moving_nest.F90 npe=",I0)', this_pe
       !endif
    endif





    if (debug_log) print '("[INFO] WDR INIT_GRID AP9 fv_moving_nest.F90 npe=",I0)', this_pe

    !if (ngrids > 1) call setup_update_regions   ! Originally from fv_control.F90
    call mn_setup_update_regions(Atm, n, nest_domain)

    call mpp_clock_end (id_reset3)
    call mpp_clock_begin (id_reset4)


    if (Atm(n)%neststruct%nested) then
       if (debug_log) print '("[INFO] WDR INIT_GRID setup_aligned_nestA fv_moving_nest.F90 npe=",I0)', this_pe

       ! New code from fv_control.F90
       ! call init_grid(Atm(this_grid), Atm(this_grid)%flagstruct%grid_name, Atm(this_grid)%flagstruct%grid_file, &
       !    Atm(this_grid)%flagstruct%npx, Atm(this_grid)%flagstruct%npy, Atm(this_grid)%flagstruct%npz, Atm(this_grid)%flagstruct%ndims, Atm(this_grid)%flagstruct%ntiles, Atm(this_grid)%ng, tile_coarse)

       ! Atm(n)%neststruct%parent_tile            = tile_coarse(n)  


       ! Old Code
       !call init_grid(Atm(n), Atm(n)%flagstruct%grid_name, Atm(n)%flagstruct%grid_file, &
       !     Atm(n)%npx, Atm(n)%npy, Atm(n)%npz, Atm(n)%flagstruct%ndims, Atm(n)%flagstruct%ntiles, Atm(n)%ng)


       !tile_coarse(1) = Atm(n)%neststruct%parent_tile
       tile_coarse(1) = parent_tile
       tile_coarse(2) = parent_tile

       call init_grid(Atm(n), Atm(n)%flagstruct%grid_name, Atm(n)%flagstruct%grid_file, &
            Atm(n)%flagstruct%npx, Atm(n)%flagstruct%npy, Atm(n)%flagstruct%npz, &
            Atm(n)%flagstruct%ndims, Atm(n)%flagstruct%ntiles, Atm(n)%ng, tile_coarse)
       if (debug_log) print '("[INFO] WDR INIT_GRID setup_aligned_nestB fv_moving_nest.F90 npe=",I0)', this_pe
    end if

    call mpp_clock_end (id_reset4)
    call mpp_clock_begin (id_reset5)

    !  Reset the gridstruct values for the nest
    if (is_fine_pe) then
       if (debug_log) print '("[INFO] WDR INIT_GRID AA fv_moving_nest.F90 npe=",I0)', this_pe
       if (debug_log) print '("[INFO] WDR INIT_GRID BB fv_moving_nest.F90 npe=",I0)', this_pe

       call grid_utils_init(Atm(n), Atm(n)%npx, Atm(n)%npy, Atm(n)%npz, &
            Atm(n)%flagstruct%non_ortho, Atm(n)%flagstruct%grid_type, Atm(n)%flagstruct%c2l_ord)

       if (debug_log) print '("[INFO] WDR INIT_GRID CC fv_moving_nest.F90 npe=",I0)', this_pe
    end if
    
    call mpp_clock_end (id_reset5)
    call mpp_clock_begin (id_reset6)


    if (debug_log) print '("[INFO] WDR NEST_DOMAIN ZZ fv_moving_nest.F90 npe=",I0)', this_pe

    if (debug_log) print '("[INFO] WDR REINIT1 CT fv_moving_nest.F90. npe=",I0," twowaynest=",L1" Atm(1)%neststruct%parent_tile=",I0)', &
         this_pe, Atm(1)%neststruct%twowaynest, Atm(1)%neststruct%parent_tile

    if (debug_log) print '("[INFO] WDR REINIT2 CT fv_moving_nest.F90. npe=",I0," twowaynest=",L1," Atm(2)%neststruct%parent_tile=",I0," n=",I0)', &
         this_pe, Atm(2)%neststruct%twowaynest, Atm(2)%neststruct%parent_tile, n

    !call mpp_sync(full_pelist)   ! Used to make debugging easier.  Can be removed.

    ! Needs to run for parent and nest Atm(2)
    !    Nest PEs update ind_update_h  -- this now seems obsolete
    !    Parent tile PEs update isu, ieu, jsu, jeu
    !    Global tiles that are not parent have no changes
    if (debug_log) print '("[INFO] WDR REINIT CV fv_moving_nest.F90. npe=",I0, " n=",I0)', this_pe, n


    ! WDR  This is now accomplished with the earlier call to setup_update_regions()
    !call reinit_parent_indices(Atm(2))
    !!call reinit_parent_indices(Atm(n))
    !if (debug_log) print '("[INFO] WDR REINIT CW fv_moving_nest.F90. npe=",I0)', this_pe

    do nn = 1, size(Atm)
       if (debug_log) call show_atm("3", Atm(nn), nn, this_pe)
    end do


    ! Output the center lat/lon of the nest
    !   only the PE that holds the center point will output this information to the logfile
    !   lat = agrid(:,:,2) and lon = agrid(:,:,1), in radians
    if (is_fine_pe) then
       half_x = Atm(child_grid_num)%npx / 2
       half_y = Atm(child_grid_num)%npy / 2

       if (half_x .ge. Atm(child_grid_num)%bd%is .and. half_x .le. Atm(child_grid_num)%bd%ie .and. half_y .ge. Atm(child_grid_num)%bd%js .and. half_y .le. Atm(child_grid_num)%bd%je) then

           half_lat = Atm(child_grid_num)%gridstruct%agrid(half_x, half_y,2) * rad2deg
           half_lon = Atm(child_grid_num)%gridstruct%agrid(half_x, half_y,1) * rad2deg
           if (half_lon .gt. 180.0) half_lon = half_lon - 360.0
           
           print '("[INFO] fv_moving_nest.F90 NEST MOVED to npe=",I0," x=",I0," y=",I0," lat=",F6.2," lon=",F7.2," a_step=",I8," fcst_hr=",F12.3)', this_pe, \
           half_x, half_y, half_lat, half_lon, a_step,  a_step * dt_atmos / 3600.0
        end if

    end if


    ! Reallocate the halo buffers in the neststruct, as some are now the wrong size
    !   Optimization would be to only deallocate the edges that have changed.  


    ! TODO Write comments on the t0 and t1 buffers
    call mpp_clock_end (id_reset6)
    call mpp_clock_begin (id_reset7)

    if (is_fine_pe) then
       !call reallocate_BC_buffers(Atm(child_grid_num))
       call reallocate_BC_buffers(Atm(1))
       if (debug_log) print '("[INFO] WDR INIT_GRID DD fv_moving_nest.F90 npe=",I0)', this_pe

       ! Reallocate buffers that are declared in fv_nesting.F90
       call dealloc_nested_buffers(Atm(1))

       if (debug_log) print '("[INFO] WDR INIT_GRID EE fv_moving_nest.F90 npe=",I0)', this_pe

       ! Set both to true so the call to setup_nested_grid_BCs() (at the beginning of fv_dynamics()) will reset t0 buffers
       ! They will be returned to false by setup_nested_grid_BCs()

       if (debug_log) print '("[INFO] WDR RESET_BCs first_step=.true.  fv_moving_nest.F90 npe=",I0)', this_pe
       Atm(n)%neststruct%first_step = .true.
       !Atm(n)%flagstruct%make_nh= .true. 

       !! Fill in the BC time1 buffers
       !call setup_nested_grid_BCs(npx, npy, npz, zvir, ncnst, &
       !     u, v, w, pt, delp, delz, q, uc, vc, pkz, &
       !     neststruct%nested, flagstruct%inline_q, flagstruct%make_nh, ng, &
       !     gridstruct, flagstruct, neststruct, &
       !     neststruct%nest_timestep, neststruct%tracer_nest_timestep, &
       !     domain, bd, nwat)

       ! Transfer the BC time1 buffers to time0

       !call set_NH_BCs_t0(neststruct)
       !call set_BCs_t0(ncnst, flagstruct%hydrostatic, neststruct)

    end if
    call mpp_clock_end (id_reset7)

  end subroutine mn_meta_reset_gridstruct

  ! WDR Copied and adapted from fv_control.F90; where it is an internal subroutine
  ! Modifications only to pass necessary variables as arguments
  subroutine mn_setup_update_regions(Atm, this_grid, nest_domain)
    type(fv_atmos_type), allocatable, intent(INOUT) :: Atm(:)
    integer, intent(IN)                             :: this_grid
    type(nest_domain_type),     intent(in)          :: nest_domain

    integer :: isu, ieu, jsu, jeu ! update regions
    integer :: isc, jsc, iec, jec
    integer :: upoff
    integer :: ngrids, n, nn

    integer :: this_pe

    this_pe = mpp_pe()

    ! Need to get the following variables from nest_domain
    !   tile_coarse() 
    !   icount_coarse() 
    !         from mpp_define_nest_domains.inc:  iend_coarse(n) = istart_coarse(n) + icount_coarse(n) - 1
    !         rearrange to: iend_coarse(n) - istart_coarse(n) + 1 = icount_coarse(n)
    !   jcount_coarse()
    !   nest_ioffsets()
    !      in fv_control.F90. pass nest_ioffsets as istart_coarse
    !   nest_joffsets()

    isc = Atm(this_grid)%bd%isc
    jsc = Atm(this_grid)%bd%jsc
    iec = Atm(this_grid)%bd%iec
    jec = Atm(this_grid)%bd%jec

    upoff = Atm(this_grid)%neststruct%upoff

    ngrids = size(Atm)

    if (debug_log) print '("[INFO] WDR SUR fv_moving_nest.F90. npe=",I0," ngrids=",I0," nest_domain%tile_coarse(",I0,"-",I0,")")', this_pe, ngrids, lbound(nest_domain%tile_coarse), ubound(nest_domain%tile_coarse)

    if (debug_log) print '("[INFO] WDR tile_coarse fv_moving_nest.F90 npe=",I0," tile_coarse(",I0,"-",I0") ngrids=",I0," tile_coarse(1)=",I0)', this_pe, &
         lbound(nest_domain%tile_coarse,1), ubound(nest_domain%tile_coarse,1), ngrids, nest_domain%tile_coarse(1)

    if (debug_log) print '("[INFO] WDR tile_coarse fv_moving_nest.F90 npe=",I0," istart_coarse(",I0,"-",I0")")', this_pe, &
         lbound(nest_domain%istart_coarse,1), ubound(nest_domain%istart_coarse,1)

    do n=2,ngrids
       nn = n - 1  !  WDR TODO revise this to handle multiple nests.  This adjusts to match fv_control.F90 where these 
       !  arrays are passed in to mpp_define_nest_domains with bounds (2:ngrids)

       if (nest_domain%tile_coarse(nn) == Atm(this_grid)%global_tile) then

          !isu = nest_ioffsets(n)
          isu = nest_domain%istart_coarse(nn)
          !ieu = isu + icount_coarse(n) - 1
          ieu = isu + (nest_domain%iend_coarse(nn) - nest_domain%istart_coarse(nn) + 1) - 1

          !jsu = nest_joffsets(n)
          jsu = nest_domain%jstart_coarse(nn)
          !jeu = jsu + jcount_coarse(n) - 1
          jeu = jsu + (nest_domain%jend_coarse(nn) - nest_domain%jstart_coarse(nn) + 1) - 1

          !update offset adjustment
          isu = isu + upoff
          ieu = ieu - upoff
          jsu = jsu + upoff
          jeu = jeu - upoff

          !restriction to current domain
!!$             !!! DEBUG CODE
!!$             if (Atm(this_grid)%flagstruct%fv_debug) then
!!$                write(*,'(I, A, 4I)') mpp_pe(), 'SETUP_UPDATE_REGIONS  : ', isu, jsu, ieu, jeu
!!$                write(*,'(I, A, 4I)') mpp_pe(), 'SETUP_UPDATE_REGIONS 2: ', isc, jsc, iec, jsc
!!$             endif
!!$             !!! END DEBUG CODE
          if (isu > iec .or. ieu < isc .or. &
               jsu > jec .or. jeu < jsc ) then
             isu = -999 ; jsu = -999 ; ieu = -1000 ; jeu = -1000
          else
             isu = max(isu,isc) ; jsu = max(jsu,jsc)
             ieu = min(ieu,iec) ; jeu = min(jeu,jec)
          endif
!!$             !!! DEBUG CODE
!!$             if (Atm(this_grid)%flagstruct%fv_debug) &
!!$                  write(*,'(I, A, 4I)') mpp_pe(), 'SETUP_UPDATE_REGIONS 3: ', isu, jsu, ieu, jeu
!!$             !!! END DEBUG CODE

          Atm(n)%neststruct%isu = isu
          Atm(n)%neststruct%ieu = ieu
          Atm(n)%neststruct%jsu = jsu
          Atm(n)%neststruct%jeu = jeu
       endif
    enddo

  end subroutine mn_setup_update_regions


  !==================================================================================================
  !
  !  Recalculation Section -- Buffers that have to change size after nest motion
  !
  !==================================================================================================


  ! Deallocate buffers.  Thought they would be reallocated in boundary.F90 nested_grid_BC_recv() when needed, but seem not to.
  ! Meant to be called after the nest has shifted;  some of these buffers will then have the wrong size.
  subroutine reallocate_BC_buffers(Atm)
    type(fv_atmos_type), intent(inout)  :: Atm

    integer :: n, ns
    logical :: dummy = .false. ! same as grids_on_this_pe(n)

    call deallocate_fv_nest_BC_type(Atm%neststruct%delp_BC)
    call deallocate_fv_nest_BC_type(Atm%neststruct%u_BC)
    call deallocate_fv_nest_BC_type(Atm%neststruct%v_BC)
    call deallocate_fv_nest_BC_type(Atm%neststruct%uc_BC)
    call deallocate_fv_nest_BC_type(Atm%neststruct%vc_BC)
    call deallocate_fv_nest_BC_type(Atm%neststruct%divg_BC)

    if (allocated(Atm%neststruct%q_BC)) then
       do n=1,size(Atm%neststruct%q_BC)
          call deallocate_fv_nest_BC_type(Atm%neststruct%q_BC(n))
       enddo
    endif

#ifndef SW_DYNAMICS
    call deallocate_fv_nest_BC_type(Atm%neststruct%pt_BC)
#ifdef USE_COND
    call deallocate_fv_nest_BC_type(Atm%neststruct%q_con_BC)
#ifdef MOIST_CAPPA
    call deallocate_fv_nest_BC_type(Atm%neststruct%cappa_BC)
#endif
#endif
    if (.not.Atm%flagstruct%hydrostatic) then
       call deallocate_fv_nest_BC_type(Atm%neststruct%w_BC)
       call deallocate_fv_nest_BC_type(Atm%neststruct%delz_BC)
    endif
#endif

    ! Reallocate the buffers

    ns = Atm%neststruct%nsponge

    call allocate_fv_nest_BC_type(Atm%neststruct%delp_BC,Atm,ns,0,0,dummy)
    call allocate_fv_nest_BC_type(Atm%neststruct%u_BC,Atm,ns,0,1,dummy)
    call allocate_fv_nest_BC_type(Atm%neststruct%v_BC,Atm,ns,1,0,dummy)
    call allocate_fv_nest_BC_type(Atm%neststruct%uc_BC,Atm,ns,1,0,dummy)
    call allocate_fv_nest_BC_type(Atm%neststruct%vc_BC,Atm,ns,0,1,dummy)
    call allocate_fv_nest_BC_type(Atm%neststruct%divg_BC,Atm,ns,1,1,dummy)

    !  if (ncnst > 0) then
    !     allocate(Atm%neststruct%q_BC(ncnst))
    !     do n=1,ncnst
    !        call allocate_fv_nest_BC_type(Atm%neststruct%q_BC(n),Atm,ns,0,0,dummy)
    !     enddo
    !  endif


    if (allocated(Atm%neststruct%q_BC)) then
       do n=1,size(Atm%neststruct%q_BC)
          call allocate_fv_nest_BC_type(Atm%neststruct%q_BC(n),Atm,ns,0,0,dummy)
       enddo
    endif



#ifndef SW_DYNAMICS
    call allocate_fv_nest_BC_type(Atm%neststruct%pt_BC,Atm,ns,0,0,dummy)
#ifdef USE_COND
    call allocate_fv_nest_BC_type(Atm%neststruct%q_con_BC,Atm,ns,0,0,dummy)
#ifdef MOIST_CAPPA
    call allocate_fv_nest_BC_type(Atm%neststruct%cappa_BC,Atm,ns,0,0,dummy)
#endif
#endif
    if (.not.Atm%flagstruct%hydrostatic) then
       call allocate_fv_nest_BC_type(Atm%neststruct%w_BC,Atm,ns,0,0,dummy)
       call allocate_fv_nest_BC_type(Atm%neststruct%delz_BC,Atm,ns,0,0,dummy)
    endif
#endif

  end subroutine reallocate_BC_buffers

#ifdef OLDCODE

  subroutine reinit_parent_indices(Atm)
    type(fv_atmos_type), intent(inout)   :: Atm

    integer :: i,j
    integer :: isc_p, iec_p, jsc_p, jec_p
    integer :: upoff, ioffset, joffset, refinement, npx, npy
    integer :: jind

    integer :: this_pe
    this_pe = mpp_pe()

    if (debug_log) print '("[INFO] WDR REINIT CT1 fv_moving_nest.F90. npe=",I0)', this_pe

    !if (debug_log) print '("[INFO] WDR REINIT CT fv_moving_nest.F90. npe=",I0," twowaynest=",L1," Atm%parent_grid%tile=",I0," Atm%neststruct%parent_tile=",I0)', &
    !  this_pe, Atm%neststruct%twowaynest, Atm%parent_grid%tile, Atm%neststruct%parent_tile

    !if (debug_log) print '("[INFO] WDR REINIT CT2 fv_moving_nest.F90. npe=",I0," twowaynest=",L1," Atm%neststruct%parent_tile=",I0)', &
    !        this_pe, Atm%neststruct%twowaynest, Atm%neststruct%parent_tile

    !if (debug_log) print '("[INFO] WDR REINIT CT3 fv_moving_nest.F90. npe=",I0," twowaynest=",L1," Atm%parent_grid%tile=",I0," Atm%neststruct%parent_tile=",I0)', &
    !        this_pe, Atm%neststruct%twowaynest, Atm%parent_grid%tile, Atm%neststruct%parent_tile


    ! This code is from fv_control.F90
    if (Atm%neststruct%twowaynest) then

       if (debug_log) print '("[INFO] WDR REINIT AA fv_moving_nest.F90. npe=",I0)', this_pe

       ioffset = Atm%neststruct%ioffset
       joffset = Atm%neststruct%joffset
       refinement = Atm%neststruct%refinement
       npx = Atm%npx
       npy = Atm%npy

       if (debug_log) print '("[INFO] WDR REINIT BB fv_moving_nest.F90. npe=",I0," ioffset=",I0," joffset=",I0)', this_pe, ioffset, joffset
       if (debug_log) print '("[INFO] WDR REINIT CC fv_moving_nest.F90. npe=",I0," isu=",I0," ieu",I0," jsu=",I0," jeu",I0)', &
            this_pe, Atm%neststruct%isu, Atm%neststruct%ieu, Atm%neststruct%jsu, Atm%neststruct%jeu

       !This in reality should be very simple. With the
       ! restriction that only the compute domain data is
       ! sent from the coarse grid, we can compute
       ! exactly which coarse grid cells should use
       ! which nested-grid data. We then don't need to send around p_ind.

       !! WDR 12/17/2020 comment out obsolete ind_update_h
       !Atm%neststruct%ind_update_h = -99999

       !if (debug_log) print '("[INFO] WDR REINIT CT fv_moving_nest.F90. npe=",I0," Atm%parent_grid%tile=",I0," Atm%neststruct%parent_tile=",I0)', &
       !     this_pe, Atm%parent_grid%tile, Atm%neststruct%parent_tile

       if (Atm%parent_grid%tile == Atm%neststruct%parent_tile) then
          if (debug_log) print '("[INFO] WDR REINIT DD fv_moving_nest.F90. npe=",I0)', this_pe

          !! Reinitialize the ind_update_h values 

          isc_p = Atm%parent_grid%bd%isc
          iec_p = Atm%parent_grid%bd%iec
          jsc_p = Atm%parent_grid%bd%jsc
          jec_p = Atm%parent_grid%bd%jec
          upoff = Atm%neststruct%upoff

          if (debug_log) print '("[INFO] WDR REINIT JSU fv_moving_nest.F90. npe=",I0)', this_pe
          !! WDR 12/17/2020 comment out obsolete ind_update_h

          Atm%neststruct%jsu = jsc_p
          Atm%neststruct%jeu = jsc_p-1
          do j=jsc_p,jec_p+1
             if (j < joffset+upoff) then
                !do i=isc_p,iec_p+1
                !   Atm%neststruct%ind_update_h(i,j,2) = -9999
                !enddo
                Atm%neststruct%jsu = Atm%neststruct%jsu + 1
                !elseif (j > joffset + (npy-1)/refinement - upoff) then
                !   !do i=isc_p,iec_p+1
                !   !   Atm%neststruct%ind_update_h(i,j,2) = -9999
                !   !enddo
             else
                jind = (j - joffset)*refinement + 1
                !do i=isc_p,iec_p+1
                !   Atm%neststruct%ind_update_h(i,j,2) = jind
                !enddo
                if ( (j < joffset + (npy-1)/refinement - upoff) .and. j <= jec_p)  Atm%neststruct%jeu = j
             endif
             !write(mpp_pe()+4000,*) j, joffset, upoff, Atm%neststruct%ind_update_h(isc_p,j,2)
          enddo

          Atm%neststruct%isu = isc_p
          Atm%neststruct%ieu = isc_p-1
          do i=isc_p,iec_p+1
             if (i < ioffset+upoff) then
                !Atm%neststruct%ind_update_h(i,:,1) = -9999
                if (debug_log) print '("[INFO] WDR REINIT ISU fv_moving_nest.F90. npe=",I0)', this_pe
                Atm%neststruct%isu = Atm%neststruct%isu + 1
                !elseif (i > ioffset + (npx-1)/refinement - upoff) then
                !   Atm%neststruct%ind_update_h(i,:,1) = -9999
             else
                !Atm%neststruct%ind_update_h(i,:,1) = (i-ioffset)*refinement + 1
                if ( (i < ioffset + (npx-1)/refinement - upoff) .and. i <= iec_p) Atm%neststruct%ieu = i
             end if
             !write(mpp_pe()+5000,*) i, ioffset, upoff, Atm%neststruct%ind_update_h(i,jsc_p,1)
          enddo
       end if
    end if

    if (debug_log) print '("[INFO] WDR REINIT ZZ fv_moving_nest.F90. npe=",I0," isu=",I0," ieu",I0," jsu=",I0," jeu",I0)', &
         this_pe, Atm%neststruct%isu, Atm%neststruct%ieu, Atm%neststruct%jsu, Atm%neststruct%jeu



  end subroutine reinit_parent_indices

#endif




  !!============================================================================                                            
  !!  Step 8 -- Moving Nest Output to NetCDF
  !!============================================================================       


  subroutine mn_prog_dump_to_netcdf(Atm, time_val, file_prefix, is_fine_pe, domain_coarse, domain_fine, nz)
    type(fv_atmos_type), intent(in)            :: Atm
    integer, intent(in)                        :: time_val
    character(len=*), intent(in)               :: file_prefix
    logical, intent(in)                        :: is_fine_pe
    type(domain2d), intent(in)                 :: domain_coarse, domain_fine
    integer, intent(in)                        :: nz

    integer            :: n_moist
    character(len=16)  :: out_var_name
    integer            :: position = CENTER 
    !integer            :: position_u = NORTH
    !integer            :: position_v = EAST

    !call mn_var_dump_to_netcdf(Atm%pt   , is_fine_pe, domain_coarse, domain_fine, position, nz, &
    !     time_val, Atm%global_tile, file_prefix, "tempK")
    !call mn_var_dump_to_netcdf(Atm%delp , is_fine_pe, domain_coarse, domain_fine, position, nz, &
    !     time_val, Atm%global_tile, file_prefix, "DELP")
    call mn_var_dump_to_netcdf(Atm%delz , is_fine_pe, domain_coarse, domain_fine, position, nz, &
         time_val, Atm%global_tile, file_prefix, "DELZ")
    call mn_var_dump_to_netcdf(Atm%q_con, is_fine_pe, domain_coarse, domain_fine, position, nz, &
         time_val, Atm%global_tile, file_prefix, "qcon")

    !call mn_var_dump_to_netcdf(Atm%w    , is_fine_pe, domain_coarse, domain_fine, position, nz, &
    !     time_val, Atm%global_tile, file_prefix, "WWND")
    !call mn_var_dump_to_netcdf(Atm%ua   , is_fine_pe, domain_coarse, domain_fine, position, nz, &
    !     time_val, Atm%global_tile, file_prefix, "UA")
    !call mn_var_dump_to_netcdf(Atm%va   , is_fine_pe, domain_coarse, domain_fine, position, nz, &
    !     time_val, Atm%global_tile, file_prefix, "VA")

    call mn_var_dump_to_netcdf(Atm%ps   , is_fine_pe, domain_coarse, domain_fine, position, 1 , &
         time_val, Atm%global_tile, file_prefix, "PS")

    !! TODO figure out what to do with ze0;  different bounds - only compute domain

    !! TODO Wind worked fine when in its own file.  Can it merge in with the regular file??
    !!call mn_var_dump_to_netcdf(Atm%u, is_fine_pe, domain_coarse, domain_fine, position_u, nz, &
    !!     time_val, Atm%global_tile, "wxvarU", "UWND")
    !!call mn_var_dump_to_netcdf(Atm%v, is_fine_pe, domain_coarse, domain_fine, position_v, nz, &
    !!     time_val, Atm%global_tile, "wxvarU", "VWND")



    ! Latitude and longitude in radians
    call mn_var_dump_to_netcdf( Atm%gridstruct%agrid(:,:,2), is_fine_pe, domain_coarse, domain_fine, position, nz, &
         time_val, Atm%global_tile, file_prefix, "latrad")
    call mn_var_dump_to_netcdf( Atm%gridstruct%agrid(:,:,1), is_fine_pe, domain_coarse, domain_fine, position, nz, &
         time_val, Atm%global_tile, file_prefix, "lonrad")

    !do n_moist = lbound(Atm%q, 4), ubound(Atm%q, 4)
    !   call get_tracer_names(MODEL_ATMOS, n_moist, out_var_name)
    !   call mn_var_dump_to_netcdf( Atm%q(:,:,:,n_moist), is_fine_pe, domain_coarse, domain_fine, position, nz, &
    !        time_val, Atm%global_tile, file_prefix, trim(out_var_name))
    !end do

  end subroutine mn_prog_dump_to_netcdf


  subroutine mn_phys_dump_to_netcdf(Atm, Atm_block, IPD_Control, IPD_Data, time_val, file_prefix, is_fine_pe, domain_coarse, domain_fine, nz)
    type(fv_atmos_type), intent(in)            :: Atm
    type (block_control_type), intent(in)      :: Atm_block
    type(IPD_control_type), intent(in)         :: IPD_Control
    type(IPD_data_type), intent(in)            :: IPD_Data(:)
    integer, intent(in)                        :: time_val
    character(len=*), intent(in)               :: file_prefix
    logical, intent(in)                        :: is_fine_pe
    type(domain2d), intent(in)                 :: domain_coarse, domain_fine
    integer, intent(in)                        :: nz

    integer :: is, ie, js, je
    integer :: nb, blen, i, j, k, ix, nv 
    integer :: this_pe

    integer            :: n_moist
    character(len=16)  :: out_var_name, phys_var_name
    integer            :: position = CENTER 
    !integer            :: position_u = NORTH
    !integer            :: position_v = EAST


    ! TODO PHYSICS Add processing of physics variables, following example in mn_prog_dump_to_netcdf
    !real (kind=kind_phys), allocatable :: smc_pr_local (:,:,:)  !< soil moisture content
    !real (kind=kind_phys), allocatable :: stc_pr_local (:,:,:)  !< soil temperature
    !real (kind=kind_phys), allocatable :: slc_pr_local (:,:,:)  !< soil liquid water content
    

    ! Coerce the high precision variables from physics into regular precision for debugging netCDF output
    ! Does not affect values used in calculations.
    ! TODO do we want to dump these as double precision??
    real, allocatable :: smc_pr_local (:,:,:)  !< soil moisture content
    real, allocatable :: stc_pr_local (:,:,:)  !< soil temperature
    real, allocatable :: slc_pr_local (:,:,:)  !< soil liquid water content

    real, allocatable, dimension(:,:) :: sealand_pr_local, deep_soil_t_pr_local, soil_type_pr_local, veg_type_pr_local, slope_type_pr_local, max_snow_alb_pr_local

    real, allocatable, dimension(:,:) ::  tsfco_pr_local, tsfcl_pr_local, tsfc_pr_local, vegfrac_pr_local, alvsf_pr_local, alvwf_pr_local, facsf_pr_local, facwf_pr_local
    real, allocatable, dimension(:,:) :: tref_pr_local, c_0_pr_local, xt_pr_local,  xu_pr_local,  xv_pr_local, ifd_pr_local

    real, allocatable :: phy_f2d_pr_local (:,:,:)
    real, allocatable :: phy_f3d_pr_local (:,:,:,:)

    this_pe = mpp_pe()

    !  Skin temp/SST
    call mn_var_dump_to_netcdf(Atm%ts   , is_fine_pe, domain_coarse, domain_fine, position, 1, &
         time_val, Atm%global_tile, file_prefix, "SSTK")

    !  Terrain height == phis / grav
    call mn_var_dump_to_netcdf(Atm%phis / grav   , is_fine_pe, domain_coarse, domain_fine, position, 1, &
         time_val, Atm%global_tile, file_prefix, "orog")
    
    ! sgh and oro were only fully allocated if fv_land is True
    !      if false, oro is (1,1), and sgh is not allocated
    if ( Atm%flagstruct%fv_land ) then
       !print '("[INFO] WDR mn_phys_dump_to_netcdf fv_land TRUE npe=",I0," size(oro)=(",I0,",",I0,")")', this_pe, size(Atm%oro, 1), size(Atm%oro, 1)
       !print '("[INFO] WDR mn_phys_dump_to_netcdf fv_land TRUE npe=",I0," size(sgh)=(",I0,",",I0,")")', this_pe, size(Atm%sgh, 1), size(Atm%sgh, 1)
       ! land frac --  called oro in fv_array.F90
       call mn_var_dump_to_netcdf(Atm%oro   , is_fine_pe, domain_coarse, domain_fine, position, 1, &
            time_val, Atm%global_tile, file_prefix, "LFRAC")
       
       ! terrain standard deviation --  called sgh in fv_array.F90
       call mn_var_dump_to_netcdf(Atm%sgh   , is_fine_pe, domain_coarse, domain_fine, position, 1, &
            time_val, Atm%global_tile, file_prefix, "STDDEV")
    !else 
       !print '("[INFO] WDR mn_phys_dump_to_netcdf fv_land FALSE npe=",I0, " size(oro)=(",I0,",",I0,")")', this_pe, size(Atm%oro, 1), size(Atm%oro, 1)
    end if

    ! Latitude and longitude in radians
    !call mn_var_dump_to_netcdf( Atm%gridstruct%agrid(:,:,2), is_fine_pe, domain_coarse, domain_fine, position, nz, &
    !     time_val, Atm%global_tile, file_prefix, "latrad")
    !call mn_var_dump_to_netcdf( Atm%gridstruct%agrid(:,:,1), is_fine_pe, domain_coarse, domain_fine, position, nz, &
    !     time_val, Atm%global_tile, file_prefix, "lonrad")



    is = Atm%bd%is
    ie = Atm%bd%ie
    js = Atm%bd%js
    je = Atm%bd%je
    
    if (debug_log) print '("[INFO] WDR mn_phys_dump_to_netcdf. npe=",I0," is=",I0," ie=",I0," js=",I0," je=",I0)', this_pe, is, ie, js, je


    ! Just allocate compute domain size here for outputs;  the nest moving code also has halos added, but we don't need them here.
    if (move_physics) then
       allocate ( smc_pr_local(is:ie, js:je, IPD_Control%lsoil) )
       allocate ( stc_pr_local(is:ie, js:je, IPD_Control%lsoil) )
       allocate ( slc_pr_local(is:ie, js:je, IPD_Control%lsoil) )
       
       allocate ( sealand_pr_local(is:ie, js:je) )
       
       allocate ( phy_f2d_pr_local(is:ie, js:je, IPD_Control%ntot2d) )
       allocate ( phy_f3d_pr_local(is:ie, js:je, IPD_Control%levs, IPD_Control%ntot3d) )
       
       
       allocate ( tsfco_pr_local(is:ie, js:je) )
       allocate ( tsfcl_pr_local(is:ie, js:je) )
       allocate ( tsfc_pr_local(is:ie, js:je) )
       allocate ( vegfrac_pr_local(is:ie, js:je) )
       allocate ( alvsf_pr_local(is:ie, js:je) )
       allocate ( alvwf_pr_local(is:ie, js:je) )
       
       ! Static file variables
       allocate ( deep_soil_t_pr_local(is:ie, js:je) )
       allocate ( soil_type_pr_local(is:ie, js:je) )
       
       !allocate ( veg_frac_pr_local(is:ie, js:je) )
       allocate ( veg_type_pr_local(is:ie, js:je) )
       
       allocate ( slope_type_pr_local(is:ie, js:je) )
       
       allocate ( max_snow_alb_pr_local(is:ie, js:je) )

       allocate ( facsf_pr_local(is:ie, js:je) )
       allocate ( facwf_pr_local(is:ie, js:je) )

    end if

    if (move_nsst) then
       allocate ( tref_pr_local(is:ie, js:je) )
       allocate ( c_0_pr_local(is:ie, js:je) )
       allocate ( xt_pr_local(is:ie, js:je) )
       allocate ( xu_pr_local(is:ie, js:je) )
       allocate ( xv_pr_local(is:ie, js:je) )
       allocate ( ifd_pr_local(is:ie, js:je) )
    end if

    
    if (move_physics) then
       smc_pr_local = +99999.9
       stc_pr_local = +99999.9
       slc_pr_local = +99999.9
       
       sealand_pr_local = +99999.9
       
       phy_f2d_pr_local = +99999.9
       phy_f3d_pr_local = +99999.9
       
       tsfco_pr_local = +99999.9
       tsfcl_pr_local = +99999.9
       tsfc_pr_local = +99999.9
       vegfrac_pr_local = +99999.9
       alvsf_pr_local = +99999.9
       alvwf_pr_local = +99999.9



    end if
    if (move_nsst) then
       tref_pr_local = +99999.9
       c_0_pr_local = +99999.9
       xt_pr_local = +99999.9
       xu_pr_local = +99999.9
       xv_pr_local = +99999.9
       ifd_pr_local = +99999.9
    end if 

    do nb = 1,Atm_block%nblks
      blen = Atm_block%blksz(nb)
        do ix = 1, blen
          i = Atm_block%index(nb)%ii(ix)
          j = Atm_block%index(nb)%jj(ix)

          if (move_physics) then
             do k = 1, IPD_Control%lsoil
                ! Use real() to lower the precision
                smc_pr_local(i,j,k) = real(IPD_Data(nb)%Sfcprop%smc(ix,k))
                stc_pr_local(i,j,k) = real(IPD_Data(nb)%Sfcprop%stc(ix,k))
                slc_pr_local(i,j,k) = real(IPD_Data(nb)%Sfcprop%slc(ix,k))
             enddo
             
             sealand_pr_local(i,j) = real(IPD_Data(nb)%Sfcprop%slmsk(ix))
             
             deep_soil_t_pr_local(i, j) = IPD_data(nb)%Sfcprop%tg3(ix)
             soil_type_pr_local(i, j) = IPD_data(nb)%Sfcprop%stype(ix)
             
             !veg_frac_pr_local(i, j) = IPD_data(nb)%Sfcprop%vfrac(ix)
             veg_type_pr_local(i, j) = IPD_data(nb)%Sfcprop%vtype(ix)
             
             slope_type_pr_local(i, j) = IPD_data(nb)%Sfcprop%slope(ix)
             
             facsf_pr_local(i, j) = IPD_data(nb)%Sfcprop%facsf(ix)
             facwf_pr_local(i, j) = IPD_data(nb)%Sfcprop%facwf(ix)

             max_snow_alb_pr_local(i, j) = IPD_data(nb)%Sfcprop%snoalb(ix)
             
             
             tsfco_pr_local(i, j) = IPD_data(nb)%Sfcprop%tsfco(ix)
             tsfcl_pr_local(i, j) = IPD_data(nb)%Sfcprop%tsfcl(ix)
             tsfc_pr_local(i, j)  = IPD_data(nb)%Sfcprop%tsfc(ix)
             vegfrac_pr_local(i, j) = IPD_data(nb)%Sfcprop%vfrac(ix)
             alvsf_pr_local(i, j) = IPD_data(nb)%Sfcprop%alvsf(ix)
             alvwf_pr_local(i, j) = IPD_data(nb)%Sfcprop%alvwf(ix)
             
             do nv = 1, IPD_Control%ntot2d
                ! Use real() to lower the precision
                phy_f2d_pr_local(i,j,nv) = real(IPD_Data(nb)%Tbd%phy_f2d(ix, nv))
             end do
             
             do k = 1, IPD_Control%levs
                do nv = 1, IPD_Control%ntot3d
                   ! Use real() to lower the precision
                   phy_f3d_pr_local(i,j,k,nv) = real(IPD_Data(nb)%Tbd%phy_f3d(ix, k, nv))
                end do
             end do
          end if
          
          if (move_nsst) then
             tref_pr_local(i, j) = IPD_data(nb)%Sfcprop%tref(ix)
             c_0_pr_local(i, j) = IPD_data(nb)%Sfcprop%c_0(ix)
             xt_pr_local(i, j) = IPD_data(nb)%Sfcprop%xt(ix)
             xu_pr_local(i, j) = IPD_data(nb)%Sfcprop%xu(ix)
             xv_pr_local(i, j) = IPD_data(nb)%Sfcprop%xv(ix)

             ifd_pr_local(i, j) = IPD_data(nb)%Sfcprop%ifd(ix)
             
          end if
          
 
       enddo
    enddo


    call mn_var_dump_to_netcdf(stc_pr_local   , is_fine_pe, domain_coarse, domain_fine, position, IPD_Control%lsoil, &
         time_val, Atm%global_tile, file_prefix, "SOILT")

    if (move_physics) then
       call mn_var_dump_to_netcdf(smc_pr_local   , is_fine_pe, domain_coarse, domain_fine, position, IPD_Control%lsoil, &
            time_val, Atm%global_tile, file_prefix, "SOILM")
       
       call mn_var_dump_to_netcdf(slc_pr_local   , is_fine_pe, domain_coarse, domain_fine, position, IPD_Control%lsoil, &
            time_val, Atm%global_tile, file_prefix, "SOILL")
       
       call mn_var_dump_to_netcdf(sealand_pr_local   , is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "LMASK")
       
       call mn_var_dump_to_netcdf(deep_soil_t_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "DEEPSOIL")
       call mn_var_dump_to_netcdf(soil_type_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "SOILTP")
       !call mn_var_dump_to_netcdf(veg_frac_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "VEGFRAC")
       call mn_var_dump_to_netcdf(veg_type_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "VEGTYPE")
       call mn_var_dump_to_netcdf(slope_type_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "SLOPE")
       call mn_var_dump_to_netcdf(max_snow_alb_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "SNOWALB")
       
       
       call mn_var_dump_to_netcdf(tsfco_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "TSFCO")
       call mn_var_dump_to_netcdf(tsfcl_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "TSFCL")
       call mn_var_dump_to_netcdf(tsfc_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "TSFC")

       call mn_var_dump_to_netcdf(vegfrac_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "VEGFRAC")
       call mn_var_dump_to_netcdf(alvsf_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "ALVSF")
       call mn_var_dump_to_netcdf(alvwf_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "ALVWF")
       
       call mn_var_dump_to_netcdf(facsf_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "FACSF")
       call mn_var_dump_to_netcdf(facwf_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "FACWF")
       
       
       do nv = 1, IPD_Control%ntot2d
          write (phys_var_name, "(A4,I0.3)")  'PH2D', nv
          call mn_var_dump_to_netcdf(phy_f2d_pr_local(:,:,nv), is_fine_pe, domain_coarse, domain_fine, position, 1, &
               time_val, Atm%global_tile, file_prefix, phys_var_name)       
       end do
       
       do nv = 1, IPD_Control%ntot3d
          write (phys_var_name, "(A4,I0.3)")  'PH3D', nv
          call mn_var_dump_to_netcdf(phy_f3d_pr_local(:,:,:,nv), is_fine_pe, domain_coarse, domain_fine, position, IPD_Control%levs, &
               time_val, Atm%global_tile, file_prefix, phys_var_name)       
       end do
    endif
       
    if (move_nsst) then
       call mn_var_dump_to_netcdf(tref_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "TREF")
       call mn_var_dump_to_netcdf(c_0_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "C_0")
       call mn_var_dump_to_netcdf(xt_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "XT")
       call mn_var_dump_to_netcdf(xu_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "XU")
       call mn_var_dump_to_netcdf(xv_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "XV")

       call mn_var_dump_to_netcdf(ifd_pr_local, is_fine_pe, domain_coarse, domain_fine, position, 1, time_val, Atm%global_tile, file_prefix, "IFD")

    end if

    ! TODO does skin temp need to be deallocated?
    if (move_physics) then
       deallocate(smc_pr_local)
       deallocate(stc_pr_local)
       deallocate(slc_pr_local)
       
       deallocate(sealand_pr_local, deep_soil_t_pr_local, soil_type_pr_local, veg_type_pr_local, max_snow_alb_pr_local)
       
       deallocate(tsfco_pr_local, tsfcl_pr_local, tsfc_pr_local, vegfrac_pr_local, alvsf_pr_local, alvwf_pr_local)
       deallocate(facsf_pr_local, facwf_pr_local)

       deallocate(phy_f2d_pr_local)
       deallocate(phy_f3d_pr_local)
    end if
    if (move_nsst) deallocate(tref_pr_local, c_0_pr_local, xt_pr_local,  xu_pr_local,  xv_pr_local, ifd_pr_local)

    if (debug_log) print '("[INFO] WDR end mn_phys_dump_tp_netcdf npe=",I0)', this_pe


  end subroutine mn_phys_dump_to_netcdf


  !!  Step 8 -- Moving Nest Output Individual Variables

  subroutine mn_var_dump_3d_to_netcdf( data_var, is_fine_pe, domain_coarse, domain_fine, position, nz, time_step, this_tile, file_prefix, var_name)
    implicit none

    !real, allocatable, intent(in)               :: data_var(:,:,:)
    real, intent(in)                            :: data_var(:,:,:)
    logical, intent(in)                         :: is_fine_pe
    type(domain2d), intent(in)                  :: domain_coarse, domain_fine
    integer, intent(in)                         :: position, nz, time_step, this_tile
    character(len=*)                            :: file_prefix, var_name

    integer                      :: isd_coarse, ied_coarse, jsd_coarse, jed_coarse
    integer                      :: isd_fine, ied_fine, jsd_fine, jed_fine
    integer                      :: this_pe
    character(len=64)            :: prefix_fine, prefix_coarse

    this_pe = mpp_pe()


    prefix_fine = trim(file_prefix) // "_fine"
    prefix_coarse = trim(file_prefix) // "_coarse"

    !!===========================================================
    !!
    !! Output the grid data from both nest grids and parent grids to netCDF 
    !!
    !!===========================================================

    if (is_fine_pe) then
       call mpp_get_data_domain(domain_fine, isd_fine, ied_fine, jsd_fine, jed_fine, position=position)

       if (debug_log) print '("[INFO] WDR NRF FG mn_var_dump_to_netcdf start. npe=",I0," size of x=(",I0,",",I0,",",I0")")', this_pe, size(data_var,1), size(data_var,2), size(data_var,3)
       if (debug_log) print '("[INFO] WDR NRF FG mn_var_dump_to_netcdf start. npe=",I0," Data domain i=",I0,"-",I0," j=",I0,"-",I0," (",I0,",",I0,")")', &
            this_pe, isd_fine, ied_fine, jsd_fine, jed_fine, ied_fine - isd_fine + 1,  jed_fine - jsd_fine + 1

       call output_grid_to_nc("GH", isd_fine, ied_fine, jsd_fine, jed_fine, nz, data_var, prefix_fine, var_name, time_step, domain_fine, position)

    else
       if (this_tile == 6) then
          !call mpp_get_compute_domain(domain_coarse, isc_coarse, iec_coarse, jsc_coarse, jec_coarse, position=position)
          call mpp_get_data_domain(domain_coarse, isd_coarse, ied_coarse, jsd_coarse, jed_coarse, position=position)
          !call mpp_get_memory_domain(domain_coarse, ism_coarse, iem_coarse, jsm_coarse, jem_coarse, position=position)

          if (debug_log) print '("[INFO] WDR NRF CG mn_var_dump_to_netcdf start. npe=",I0," size of x=(",I0,",",I0,",",I0")")', this_pe, size(data_var,1), size(data_var,2), size(data_var,3)
          if (debug_log) print '("[INFO] WDR NRF CG mn_var_dump_to_netcdf start. npe=",I0," Data domain i=",I0,"-",I0," j=",I0,"-",I0," (",I0,",",I0,")")', &
               this_pe, isd_coarse, ied_coarse, jsd_coarse, jed_coarse, ied_coarse - isd_coarse + 1,  jed_coarse - jsd_coarse + 1
          !if (debug_log) print '("[INFO] WDR NRF CG mn_var_dump_to_netcdf start. npe=",I0," Compute domain i=",I0,"-",I0," j=",I0,"-",I0," (",I0,",",I0,")")', &
          !     this_pe, isc_coarse, iec_coarse, jsc_coarse, jec_coarse, iec_coarse - isc_coarse + 1,  jec_coarse - jsc_coarse + 1
          !if (debug_log) print '("[INFO] WDR NRF CG mn_var_dump_to_netcdf start. npe=",I0," Memory domain i=",I0,"-",I0," j=",I0,"-",I0," (",I0,",",I0,")")', &
          !     this_pe, ism_coarse, iem_coarse, jsm_coarse, jem_coarse, iem_coarse - ism_coarse + 1,  jem_coarse - jsm_coarse + 1

          call output_grid_to_nc("GH", isd_coarse, ied_coarse, jsd_coarse, jed_coarse, nz, data_var, prefix_coarse, var_name, time_step, domain_coarse, position)

       end if
    end if

  end subroutine mn_var_dump_3d_to_netcdf


  subroutine mn_var_dump_2d_to_netcdf( data_var, is_fine_pe, domain_coarse, domain_fine, position, nz, time_step, this_tile, file_prefix, var_name)
    implicit none

    !real, allocatable, intent(in)               :: data_var(:,:,:)
    real, intent(in)                            :: data_var(:,:)
    logical, intent(in)                         :: is_fine_pe
    type(domain2d), intent(in)                  :: domain_coarse, domain_fine
    integer, intent(in)                         :: position, nz, time_step, this_tile
    character(len=*)                            :: file_prefix, var_name


    integer                      :: isc_coarse, iec_coarse, jsc_coarse, jec_coarse
    integer                      :: isd_coarse, ied_coarse, jsd_coarse, jed_coarse
    integer                      :: isd_fine, ied_fine, jsd_fine, jed_fine
    integer                      :: isc_fine, iec_fine, jsc_fine, jec_fine

    integer                      :: ism_coarse, iem_coarse, jsm_coarse, jem_coarse
    integer                      :: ism_fine, iem_fine, jsm_fine, jem_fine

    integer                      :: this_pe

    character(len=64)            :: prefix_fine, prefix_coarse

    this_pe = mpp_pe()


    prefix_fine = trim(file_prefix) // "_fine"
    prefix_coarse = trim(file_prefix) // "_coarse"

    !!===========================================================
    !!
    !! Output the grid data from both nest grids and parent grids to netCDF 
    !!
    !!===========================================================

    if (is_fine_pe) then
       ! Maybe don't need to call mpp_get_compute_domain here?
       !call mpp_get_compute_domain(domain_fine, isc_fine, iec_fine, jsc_fine, jec_fine, position=position)
       call mpp_get_data_domain(domain_fine, isd_fine, ied_fine, jsd_fine, jed_fine, position=position)
       !call mpp_get_memory_domain(domain_fine, ism_fine, iem_fine, jsm_fine, jem_fine, position=position)


       if (debug_log) print '("[INFO] WDR NRF FG mn_var_dump_to_netcdf start. npe=",I0," size of x=(",I0,",",I0,")")', this_pe, size(data_var,1), size(data_var,2)
       if (debug_log) print '("[INFO] WDR NRF FG mn_var_dump_to_netcdf start. npe=",I0," Data domain i=",I0,"-",I0," j=",I0,"-",I0," (",I0,",",I0,")")', &
            this_pe, isd_fine, ied_fine, jsd_fine, jed_fine, ied_fine - isd_fine + 1,  jed_fine - jsd_fine + 1
       !if (debug_log) print '("[INFO] WDR NRF FG mn_var_dump_to_netcdf start. npe=",I0," Compute domain i=",I0,"-",I0," j=",I0,"-",I0," (",I0,",",I0,")")', &
       !     this_pe, isc_fine, iec_fine, jsc_fine, jec_fine, iec_fine - isc_fine + 1,  jec_fine - jsc_fine + 1
       !if (debug_log) print '("[INFO] WDR NRF FG mn_var_dump_to_netcdf start. npe=",I0," Memory domain i=",I0,"-",I0," j=",I0,"-",I0," (",I0,",",I0,")")', &
       !     this_pe, ism_fine, iem_fine, jsm_fine, jem_fine, iem_fine - ism_fine + 1,  jem_fine - jsm_fine + 1

       call output_grid_to_nc("GH", isd_fine, ied_fine, jsd_fine, jed_fine, nz, data_var, prefix_fine, var_name, time_step, domain_fine, position)

    else


       if (this_tile == 6) then
          !call mpp_get_compute_domain(domain_coarse, isc_coarse, iec_coarse, jsc_coarse, jec_coarse, position=position)
          call mpp_get_data_domain(domain_coarse, isd_coarse, ied_coarse, jsd_coarse, jed_coarse, position=position)
          !call mpp_get_memory_domain(domain_coarse, ism_coarse, iem_coarse, jsm_coarse, jem_coarse, position=position)

          if (debug_log) print '("[INFO] WDR NRF CG mn_var_dump_to_netcdf start. npe=",I0," size of x=(",I0,",",I0,")")', this_pe, size(data_var,1), size(data_var,2)
          if (debug_log) print '("[INFO] WDR NRF CG mn_var_dump_to_netcdf start. npe=",I0," Data domain i=",I0,"-",I0," j=",I0,"-",I0," (",I0,",",I0,")")', &
               this_pe, isd_coarse, ied_coarse, jsd_coarse, jed_coarse, ied_coarse - isd_coarse + 1,  jed_coarse - jsd_coarse + 1
          !if (debug_log) print '("[INFO] WDR NRF CG mn_var_dump_to_netcdf start. npe=",I0," Compute domain i=",I0,"-",I0," j=",I0,"-",I0," (",I0,",",I0,")")', &
          !     this_pe, isc_coarse, iec_coarse, jsc_coarse, jec_coarse, iec_coarse - isc_coarse + 1,  jec_coarse - jsc_coarse + 1
          !if (debug_log) print '("[INFO] WDR NRF CG mn_var_dump_to_netcdf start. npe=",I0," Memory domain i=",I0,"-",I0," j=",I0,"-",I0," (",I0,",",I0,")")', &
          !     this_pe, ism_coarse, iem_coarse, jsm_coarse, jem_coarse, iem_coarse - ism_coarse + 1,  jem_coarse - jsm_coarse + 1

          call output_grid_to_nc("GH", isd_coarse, ied_coarse, jsd_coarse, jed_coarse, nz, data_var, prefix_coarse, var_name, time_step, domain_coarse, position)

       end if
    end if

  end subroutine mn_var_dump_2d_to_netcdf

  !!=========================================================================================                               
  !! Step 9 -- Perform vertical remapping on nest(s) and recalculate auxiliary pressures                                    
  !!           Should help stabilize the fields before dynamics runs                                                        
  !!=========================================================================================                               

  subroutine recalc_aux_pressures(Atm)
    type(fv_atmos_type), intent(inout) :: Atm

    !  Update the auxiliary pressure variables
    !  In nest moving code, we moved delp and delz; this will update ps, pk, pe, peln, and pkz 
    !  Note this routine makes hydrostatic calculations (but has non-hydrostatic branches)
    !  Perhaps not appropriate for a non-hydrostatic run.  
    !  May need to find or write a non-hydrostatic version of this routine

    ! TODO determine if this is the correct way to recalculate the auxiliary pressure variables

    call p_var(Atm%npz, Atm%bd%is, Atm%bd%ie, Atm%bd%js, Atm%bd%je, Atm%ptop, ptop_min,  &
         Atm%delp, Atm%delz, &
         Atm%pt, Atm%ps, &
         Atm%pe, Atm%peln,   &
         Atm%pk,   Atm%pkz, kappa, &
         Atm%q, Atm%ng, Atm%flagstruct%ncnst, Atm%gridstruct%area_64, 0.,  &
         .false.,  .false., & !mountain argument not used
         Atm%flagstruct%moist_phys,  Atm%flagstruct%hydrostatic, &
         Atm%flagstruct%nwat, Atm%domain, .false.)

  end subroutine recalc_aux_pressures

#ifdef REMAP

  subroutine vertical_remap_nest(Atm, dt_atmos, p_split)
#ifdef CCPP
    use mpp_mod,   only: FATAL, mpp_error
    use CCPP_data, only: CCPP_interstitial
#endif
    type(fv_atmos_type), intent(inout) :: Atm
    real, intent(in)    :: dt_atmos
    integer, intent(in) :: p_split


    !! Start Variables for vertical remapping

    logical  :: do_omega
    !logical :: do_adiabatic_init
    logical :: out_dt   
    !integer :: kord_tracer(ncnst)
    integer :: kord_tracer(Atm%ncnst -  Atm%flagstruct%pnats)
    integer :: iq, i, j, k
    integer :: this_pe

    integer :: sphum, liq_wat, rainwat, ice_wat, snowwat, graupel  ! condensate species tracer indices


    integer :: npz, ncnst, pnats, nq, nr


#ifdef CCPP
    integer :: cld_amt
#endif

    real, allocatable :: ws(:,:)
    real, allocatable :: pfull(:) 
    real, allocatable :: cvm(:)

    real, allocatable :: teq(:,:)
    real, allocatable :: te_local(:,:,:)

#ifndef CCPP
    real, allocatable :: cappa(:,:,:) 
    real, allocatable :: dp1(:,:,:) 
    real, allocatable :: dtdt_m(:,:,:) 
    real, allocatable :: te_2d(:,:)
    logical :: last_step
#endif

    real :: akap
    real :: bdt, mdt
    real :: ph1, ph2
    real :: zvir

#ifdef CCPP
    !! Also requires call to "end associate" at end of subroutine
    integer :: ierr

    ccpp_associate:associate(cappa=>CCPP_interstitial%cappa,&
         dp1=>CCPP_interstitial%te0,&
         dtdt_m=>CCPP_interstitial%dtdt,&
         last_step=>CCPP_interstitial%last_step,&
         te_2d=>CCPP_interstitial%te0_2d)

#endif

    this_pe = mpp_pe()

    npz   = Atm%npz
    zvir = rvgas/rdgas - 1.

    ! Moisture species
    ncnst = Atm%ncnst
    pnats = Atm%flagstruct%pnats

    nq = ncnst-pnats
    nr = nq - Atm%flagstruct%dnrts  ! from fv_dynamics.F90


    sphum   = get_tracer_index (MODEL_ATMOS, 'sphum' )
    liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat' )
    ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat' )
    rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat' )
    snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat' )
    graupel = get_tracer_index (MODEL_ATMOS, 'graupel' )
#ifdef CCPP
    cld_amt = get_tracer_index (MODEL_ATMOS, 'cld_amt')
#endif




    !  Update the auxiliary pressure variables
    !  In nest moving code, we moved delp and delz; this will update ps, pk, pe, peln, and pkz 
    !  Note this routine makes hydrostatic calculations (but has non-hydrostatic branches)
    !  Perhaps not appropriate for a non-hydrostatic run.  
    !  May need to find or write a non-hydrostatic version of this routine

    ! TODO determine if this is the correct way to recalculate the auxiliary pressure variables

    call p_var(npz, Atm%bd%is, Atm%bd%ie, Atm%bd%js, Atm%bd%je, Atm%ptop, ptop_min,  &
         Atm%delp, Atm%delz, &
         Atm%pt, Atm%ps, &
         Atm%pe, Atm%peln,   &
         Atm%pk,   Atm%pkz, kappa, &
         Atm%q, Atm%ng, Atm%flagstruct%ncnst, Atm%gridstruct%area_64, 0.,  &
         .false.,  .false., & !mountain argument not used
         Atm%flagstruct%moist_phys,  Atm%flagstruct%hydrostatic, &
         Atm%flagstruct%nwat, Atm%domain, Atm%flagstruct%adiabatic, .false.)




    ! Allocate arrays

    ! TODO verify whether ws needs to be recalculated or shifted
    allocate( ws(Atm%bd%is:Atm%bd%ie, Atm%bd%js:Atm%bd%je) )
    allocate( cvm(Atm%bd%is:Atm%bd%ie))
    allocate( teq(Atm%bd%is:Atm%bd%ie, Atm%bd%js:Atm%bd%je) )
    allocate( te_local(Atm%bd%isd:Atm%bd%ied, Atm%bd%jsd:Atm%bd%jed, npz) )
    allocate( pfull(npz) )


#ifndef CCPP
    allocate( dp1(Atm%bd%isd:Atm%bd%ied, Atm%bd%jsd:Atm%bd%jed, 1:npz) )
    dp1 = 0.0    ! set below in moist_cv loop
    allocate( te_2d(Atm%bd%is:Atm%bd%ie, Atm%bd%js:Atm%bd%je) )
    te_2d = 0.0  ! set 2D total energy in compute_total_energy

    if ( idiag%id_mdt > 0 .and. (.not. do_adiabatic_init) ) then
       allocate ( dtdt_m(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je,npz) )
       dtdt_m = 0.0 ! only used if out_dt is set to .true.
    end if

#ifdef MOIST_CAPPA
    allocate ( cappa(Atm%bd%isd:Atm%bd%ied, Atm%bd%jsd:Atm%bd%jed,npz) )
    call init_ijk_mem(Atm%bd%isd, Atm%bd%ied, Atm%bd%jsd, Atm%bd%jed, npz, cappa, 0.)
#else
    ! Dummy 1x1x1 array
    allocate ( cappa(Atm%bd%isd:Atm%bd%isd, Atm%bd%jsd:Atm%bd%jsd,1) )
    cappa = 0.
#endif

#endif


    ! Set values
    ! ws  - vertical motion w at the surface.  normally returned by dyn_core
    ! te_2d
    ! dtdt_m

    ws = 0.0     ! TODO figure out how to set w surface to a valid value
    cvm = 0.0    ! set in moist_cv loop
    teq = 0.0    ! set in compute_total_energy and then not used
    te_local = 0.0     ! te_local will be set by Lagrangian_to_Eulerian()
    pfull = 0.0  ! set immediately below
    ! cappa is allocated and filled depending on ifdef MOIST_CAPPA
    last_step = .false.
    !if ( n_map==k_split ) last_step = .true.  ! if it's the last step of k_split
    !do_adiabatic_init = .false.   ! TODO Not sure what this should be.

    do k=1,npz
       ph1 = Atm%ak(k  ) + Atm%bk(k  ) * Atm%flagstruct%p_ref
       ph2 = Atm%ak(k+1) + Atm%bk(k+1) * Atm%flagstruct%p_ref
       pfull(k) = (ph2 - ph1) / log(ph2/ph1)
    enddo

    ! Variables needed for vertical remapping Lagrangian_to_Eulerian call

    bdt = dt_atmos/real(abs(p_split))
    mdt = bdt / real(Atm%flagstruct%k_split)

    ! print '("[INFO] WDR REMAP npe=",I0," nr=",I0," nq=",I0)', this_pe, nr, nq    !! TODO Validate which is correct to use here
    
    do iq=1,nq
       kord_tracer(iq) = Atm%flagstruct%kord_tr
       if ( iq==cld_amt )  kord_tracer(iq) = 9      ! monotonic        
    end do

    do_omega = Atm%flagstruct%hydrostatic .and. last_step

    do i=Atm%bd%is, Atm%bd%ie
       do k=1,npz
          do j=Atm%bd%js, Atm%bd%je
             dp1(i,j,k) = zvir * Atm%q(i,j,k,sphum)
          end do
       end do
    end do

    ! cappa initialization code taken from fv_dynamics.F90
#ifdef MOIST_CAPPA
    if (debug_log) print '("[INFO] WDR VERT_REMAP MV_NST L2E computing moist cv before vertical remapping fv_moving_nest.F90 npe=",I0)', this_pe
    do k=1,npz
       do j=Atm%bd%js, Atm%bd%je
          call moist_cv(Atm%bd%is, Atm%bd%ie, Atm%bd%isd, Atm%bd%ied, Atm%bd%jsd, Atm%bd%jed, &
               npz, j, k, Atm%flagstruct%nwat, sphum, liq_wat, rainwat,    &
               ice_wat, snowwat, graupel, Atm%q, Atm%q_con(Atm%bd%is:Atm%bd%ie,j,k), cvm)

          do i=Atm%bd%is, Atm%bd%ie
             cappa(i,j,k) = rdgas/(rdgas + cvm(i)/(1.+dp1(i,j,k)))
          end do
       end do
    end do
#else
    if (debug_log) print '("[INFO] WDR VERT_REMAP MV_NST L2E NOT computing moist cv before vertical remapping fv_moving_nest.F90 npe=",I0)', this_pe  
#endif



#ifdef SW_DYNAMICS
    akap  = 1.
#else
    akap  = kappa
#endif

    ! Returns values for te_2d and teq.  teq is otherwise not used.
    ! The wind variables (u,v,ua,va) in compute_total_energy are intent(inout) - 
    ! but the code that would have modified them is commented out.  They are not changed.
    ! ua,va are not actually referenced in the subroutine

    if ( Atm%flagstruct%consv_te > 0.  .and. (.not. do_adiabatic_init) ) then
       if (debug_log) print '("[INFO] WDR VERT_REMAP MV_NST L2E computing total energy before vertical remapping fv_moving_nest.F90 npe=",I0)', this_pe
       call compute_total_energy(Atm%bd%is, Atm%bd%ie, Atm%bd%js, Atm%bd%je, &
            Atm%bd%isd, Atm%bd%ied, Atm%bd%jsd, Atm%bd%jed, Atm%npz,        &
            Atm%u, Atm%v, Atm%w, Atm%delz, Atm%pt, Atm%delp, Atm%q, &
            dp1, Atm%pe, Atm%peln, Atm%phis, &
            Atm%gridstruct%rsin2, Atm%gridstruct%cosa_s, &
            zvir, cp_air, rdgas, hlv, te_2d, Atm%ua, Atm%va, teq,        &
            Atm%flagstruct%moist_phys, Atm%flagstruct%nwat, sphum, liq_wat, rainwat,   &
            ice_wat, snowwat, graupel, Atm%flagstruct%hydrostatic, Atm%idiag%id_te)
    else
       if (debug_log) print '("[INFO] WDR VERT_REMAP MV_NST L2E NOT computing total energy before vertical remapping fv_moving_nest.F90 npe=",I0," te_2d=",F10.5)', this_pe, te_2d(Atm%bd%is, Atm%bd%js)
    end if

    if (debug_log) print '("[INFO] WDR VERT_REMAP MV_NST L2E before vertical remapping fv_moving_nest.F90 npe=",I0)', this_pe

    ! pkz and te are the only intent(out) arguments
    ! All the rest are intent(in) or intent(inout)

    out_dt = Atm%idiag%id_mdt>0   ! This determines whether dtdt_m is used; is false for current testing configuration

    if (debug_log) call check_local_array(cappa, this_pe, "L2E cappa", -1.0e9,1.0e9)
    if (debug_log) call check_local_array(dp1, this_pe, "L2E dp1", -1.0e9,1.0e9)
    if (debug_log) call check_local_array(dtdt_m, this_pe, "L2E dtdt_m", -1.0e9,1.0e9)
    if (debug_log) call check_local_array(te_2d, this_pe, "L2E te_2d", -1.0e9,1.0e9)

    if (debug_log) print '("[INFO] WDR VERT_REMAP L2E bounds  npe=",I0," ",A32,"(",I4,":",I4,",",I4,":",I4,",",I4,":",I4,")")', this_pe, "cappa", lbound(cappa,1), ubound(cappa,1), lbound(cappa,2), ubound(cappa,2), lbound(cappa,3), ubound(cappa,3)

    if (debug_log) print '("[INFO] WDR VERT_REMAP L2E bounds  npe=",I0," ",A32,"(",I4,":",I4,",",I4,":",I4,",",I4,":",I4,")")', this_pe, "dp1", lbound(dp1,1), ubound(dp1,1), lbound(dp1,2), ubound(dp1,2), lbound(dp1,3), ubound(dp1,3)

    if (debug_log) print '("[INFO] WDR VERT_REMAP L2E bounds  npe=",I0," ",A32,"(",I4,":",I4,",",I4,":",I4,",",I4,":",I4,")")', this_pe, "dtdt_m", lbound(dtdt_m,1), ubound(dtdt_m,1), lbound(dtdt_m,2), ubound(dtdt_m,2), lbound(dtdt_m,3), ubound(dtdt_m,3)

    if (debug_log) print '("[INFO] WDR VERT_REMAP L2E bounds  npe=",I0," ",A32,"(",I4,":",I4,",",I4,":",I4,")")', this_pe, "te_2d", lbound(te_2d,1), ubound(te_2d,1), lbound(te_2d,2), ubound(te_2d,2)


    if (debug_log) call check_array(te_local, this_pe, "L2E te_local", -1.0e9,1.0e9)
    if (debug_log) call check_array(ws, this_pe, "L2E ws", -1.0e9,1.0e9)

    if (debug_log) print '("[INFO] WDR VERT_REMAP MV_NST L2E fv_moving_nest.F90 npe=",I0," out_dt(whether dtdt_m is used)=",L1, " last_step=",L1, " do_omega=",L1)', this_pe, out_dt, last_step, do_omega

    !nq = nq_tot - Atm%flagstruct%dnats



    ! From updated fv_dynamics.F90 on June 24, 2021
!    call Lagrangian_to_Eulerian(last_step, consv_te, &
!         ps, pe, delp,          &
!         pkz, pk, mdt, bdt, npx, npy, npz, &
!         is,ie,js,je, isd,ied,jsd,jed,       &
!         nr, nwat, sphum, &
!         q_con, u,  v, w, &
!         delz, pt, q, phis,    &
!         zvir, cp_air, akap, cappa, flagstruct%kord_mt, flagstruct%kord_wz, &
!         kord_tracer, flagstruct%kord_tm, peln, te_2d,               &
!         ng, ua, va, omga, dp1, ws, &
!         fill, reproduce_sum,             &
!         idiag%id_mdt>0, dtdt_m, &
!         ptop, ak, bk, pfull, gridstruct, domain,   &
!         flagstruct%do_sat_adj, hydrostatic, flagstruct%phys_hydrostatic, &
!         hybrid_z, do_omega,     &
!         flagstruct%adiabatic, do_adiabatic_init, &
!         flagstruct%do_inline_mp, &
!         inline_mp, flagstruct%c2l_ord, bd, flagstruct%fv_debug, &
!         flagstruct%moist_phys)
    



    call Lagrangian_to_Eulerian(last_step, Atm%flagstruct%consv_te,  &
         Atm%ps, Atm%pe, Atm%delp,                              &
         Atm%pkz, Atm%pk, mdt, bdt, Atm%npx, Atm%npy, Atm%npz,  &  
         Atm%bd%is,Atm%bd%ie,Atm%bd%js,Atm%bd%je,               &
         Atm%bd%isd, Atm%bd%ied, Atm%bd%jsd, Atm%bd%jed,        &
         nq, Atm%flagstruct%nwat, sphum,                        &  ! TODO check if nq is the same as nr?
         Atm%q_con, Atm%u,  Atm%v, Atm%w,                       &
         Atm%delz, Atm%pt, Atm%q, Atm%phis,                     &
         zvir, cp_air, akap, cappa, Atm%flagstruct%kord_mt, Atm%flagstruct%kord_wz, &
         kord_tracer, Atm%flagstruct%kord_tm, Atm%peln, te_2d,                      &
         Atm%ng, Atm%ua, Atm%va, Atm%omga, te_local, ws, &  ! TODO check if te_local is the same as dp1?
         Atm%flagstruct%fill, Atm%flagstruct%reproduce_sum,  &
         Atm%idiag%id_mdt>0, dtdt_m, &
         Atm%ptop, Atm%ak, Atm%bk, pfull, Atm%gridstruct, Atm%domain,   &
         Atm%flagstruct%do_sat_adj, Atm%flagstruct%hydrostatic, Atm%flagstruct%phys_hydrostatic, &
         Atm%flagstruct%hybrid_z, do_omega,            &
         Atm%flagstruct%adiabatic, do_adiabatic_init, &
         Atm%flagstruct%do_inline_mp, &
         Atm%inline_mp, Atm%flagstruct%c2l_ord, Atm%bd, Atm%flagstruct%fv_debug, &
         Atm%flagstruct%moist_phys)



    !call vertical_remap(Atm, last_step, mdt, bdt, nq, sphum,       &
    !     zvir, cp_air, akap, cappa, kord_tracer, te_2d,               &
    !     te_local, ws, dtdt_m, pfull, do_omega, do_adiabatic_init)


    !          ! Prognostic variables
    !
    !          Atm%u,  Atm%v, Atm%w,
    !          Atm%pt,
    !          Atm%delz, 
    !          
    !          Atm%delp,              
    !          Atm%q,    
    !          Atm%q_con, 
    !          Atm%phis,
    !
    !          ! Diagnostic variables (must be recalculated or interpolated each timestep)
    !
    !          Atm%ps, 
    !          Atm%pe, 
    !          Atm%pkz, 
    !          Atm%pk, 
    !          Atm%npz, 
    !          Atm%peln, 
    !          Atm%ng, 
    !          Atm%ua, 
    !          Atm%va, 
    !          Atm%omga, 
    !
    !          ! Grid Definitions (require shifts when nest moves)
    !          
    !          Atm%bd%is,Atm%bd%ie,Atm%bd%js,Atm%bd%je,
    !          Atm%bd%isd, Atm%bd%ied, Atm%bd%jsd, Atm%bd%jed,
    !          Atm%gridstruct, 
    !          Atm%domain,
    !
    !          ! Grid Definitions (static for each model run)
    !
    !          Atm%ptop, 
    !          Atm%ak, 
    !          Atm%bk,
    !          
    !          ! Flags
    !          
    !          Atm%flagstruct%consv_te,      
    !          Atm%flagstruct%nwat,
    !          Atm%flagstruct%kord_mt, 
    !          Atm%flagstruct%kord_wz,
    !          Atm%flagstruct%kord_tm,
    !          Atm%flagstruct%fill,
    !          Atm%flagstruct%reproduce_sum,
    !          Atm%idiag%id_mdt
    !          Atm%flagstruct%do_sat_adj,
    !          Atm%flagstruct%hydrostatic,
    !          Atm%flagstruct%hybrid_z,
    !          Atm%flagstruct%adiabatic,




    !! Output some of the data for validation


    !!  This is needed to pair with the associate clause above
#ifdef CCPP
    end associate ccpp_associate
#endif


  end subroutine vertical_remap_nest

#endif REMAP


  !==================================================================================================
  !
  !  Utility Section  -- After Step 9
  !
  !==================================================================================================


  ! copied from dyn_core.F90 to avoid circular dependencies
  subroutine init_ijk_mem(i1, i2, j1, j2, km, array, var)
    integer, intent(in):: i1, i2, j1, j2, km
    real, intent(inout):: array(i1:i2,j1:j2,km)
    real, intent(in):: var
    integer:: i, j, k

    !$OMP parallel do default(none) shared(i1,i2,j1,j2,km,array,var)
    do k=1,km
       do j=j1,j2
          do i=i1,i2
             array(i,j,k) = var
          enddo
       enddo
    enddo

  end subroutine init_ijk_mem


  function almost_equal(a, b)
    logical :: almost_equal
    real, intent(in):: a,b

    real :: tolerance = 0.00001

    if ( abs(a - b) < tolerance ) then
       almost_equal = .true.
    else
       almost_equal = .false.
    end if
  end function almost_equal




  ! Shifts tile_geo values using the data from fp_super_tile_geo
  subroutine move_nest_geo(tile_geo, tile_geo_u, tile_geo_v, fp_super_tile_geo, delta_i_c, delta_j_c, x_refine, y_refine)
    implicit none
    type(grid_geometry), intent(inout)  :: tile_geo
    type(grid_geometry), intent(inout)  :: tile_geo_u
    type(grid_geometry), intent(inout)  :: tile_geo_v
    type(grid_geometry), intent(in)     :: fp_super_tile_geo
    integer, intent(in)                 :: delta_i_c, delta_j_c, x_refine, y_refine


    integer :: nest_x, nest_y, parent_x, parent_y

    type(bbox)  :: tile_bbox, fp_tile_bbox, tile_bbox_u, tile_bbox_v
    integer   :: i, j, fp_i, fp_j

    ! tile_geo is cell-centered, at nest refinement
    ! fp_super_tile_geo is a supergrid, at nest refinement

    call find_nest_alignment(tile_geo, fp_super_tile_geo, nest_x, nest_y, parent_x, parent_y)

    call fill_bbox(tile_bbox, tile_geo%lats)
    call fill_bbox(tile_bbox_u, tile_geo_u%lats)
    call fill_bbox(tile_bbox_v, tile_geo_v%lats)
    call fill_bbox(fp_tile_bbox, fp_super_tile_geo%lats)

    ! Calculate new parent alignment -- supergrid at the refine ratio
    !  delta_{i,j}_c are at the coarse center grid resolution
    parent_x = parent_x + delta_i_c * 2 * x_refine
    parent_y = parent_y + delta_j_c * 2 * y_refine


    ! Brute force repopulation of full tile_geo grids.
    ! Optimization would be to use EOSHIFT and bring in just leading edge
    do i = tile_bbox%is, tile_bbox%ie
       do j = tile_bbox%js, tile_bbox%je
          fp_i = (i - nest_x) * 2 + parent_x
          fp_j = (j - nest_y) * 2 + parent_y

          if (fp_i < fp_tile_bbox%is .or. fp_i > fp_tile_bbox%ie) then
             if (debug_log) print '("[ERROR] WDR move_nest_geo invalid fp_i=",I0," is=",I0," ie=",I0)', fp_i, fp_tile_bbox%is, fp_tile_bbox%ie
             stop  ! replace with a fatal error
          end if
          if (fp_j < fp_tile_bbox%js .or. fp_j > fp_tile_bbox%je) then
             if (debug_log) print '("[ERROR] WDR move_nest_geo invalid fp_j=",I0," js=",I0," je=",I0)', fp_j, fp_tile_bbox%js, fp_tile_bbox%je
             stop  ! replace with a fatal error
          end if

          tile_geo%lats(i,j) = fp_super_tile_geo%lats(fp_i, fp_j)
          tile_geo%lons(i,j) = fp_super_tile_geo%lons(fp_i, fp_j)
       end do
    end do

    do i = tile_bbox_u%is, tile_bbox_u%ie
       do j = tile_bbox_u%js, tile_bbox_u%je
          fp_i = (i - nest_x) * 2 + parent_x
          fp_j = (j - nest_y) * 2 + parent_y - 1

          if (fp_i < fp_tile_bbox%is .or. fp_i > fp_tile_bbox%ie) then
             if (debug_log) print '("[ERROR] WDR move_nest_geo invalid fp_i=",I0," is=",I0," ie=",I0)', fp_i, fp_tile_bbox%is, fp_tile_bbox%ie
             stop  ! replace with a fatal error
          end if
          if (fp_j < fp_tile_bbox%js .or. fp_j > fp_tile_bbox%je) then
             if (debug_log) print '("[ERROR] WDR move_nest_geo invalid fp_j=",I0," js=",I0," je=",I0)', fp_j, fp_tile_bbox%js, fp_tile_bbox%je
             stop  ! replace with a fatal error
          end if

          tile_geo_u%lats(i,j) = fp_super_tile_geo%lats(fp_i, fp_j)
          tile_geo_u%lons(i,j) = fp_super_tile_geo%lons(fp_i, fp_j)
       end do
    end do

    do i = tile_bbox_v%is, tile_bbox_v%ie
       do j = tile_bbox_v%js, tile_bbox_v%je
          fp_i = (i - nest_x) * 2 + parent_x - 1
          fp_j = (j - nest_y) * 2 + parent_y 

          if (fp_i < fp_tile_bbox%is .or. fp_i > fp_tile_bbox%ie) then
             if (debug_log) print '("[ERROR] WDR move_nest_geo invalid fp_i=",I0," is=",I0," ie=",I0)', fp_i, fp_tile_bbox%is, fp_tile_bbox%ie
             stop  ! replace with a fatal error
          end if
          if (fp_j < fp_tile_bbox%js .or. fp_j > fp_tile_bbox%je) then
             if (debug_log) print '("[ERROR] WDR move_nest_geo invalid fp_j=",I0," js=",I0," je=",I0)', fp_j, fp_tile_bbox%js, fp_tile_bbox%je
             stop  ! replace with a fatal error
          end if

          tile_geo_v%lats(i,j) = fp_super_tile_geo%lats(fp_i, fp_j)
          tile_geo_v%lons(i,j) = fp_super_tile_geo%lons(fp_i, fp_j)
       end do
    end do

    ! Validate at the end
    call find_nest_alignment(tile_geo, fp_super_tile_geo, nest_x, nest_y, parent_x, parent_y)


  end subroutine move_nest_geo

  subroutine assign_n_p_grids(parent_geo, tile_geo, p_grid, n_grid, position)
    type(grid_geometry), intent(in)          ::  parent_geo, tile_geo
    real(kind=R_GRID), allocatable, intent(inout)            :: p_grid(:,:,:)
    real(kind=R_GRID), allocatable, intent(inout)            :: n_grid(:,:,:)
    integer, intent(in)                      :: position

    integer :: i,j


    if (position == CENTER) then
       do j = lbound(tile_geo%lats,2), ubound(tile_geo%lats,2)
          do i = lbound(tile_geo%lats,1), ubound(tile_geo%lats,1)
             ! centered grid version
             n_grid(i, j, 1) = tile_geo%lons(i, j)
             n_grid(i, j, 2) = tile_geo%lats(i, j)             
             !if (debug_log) print '("[INFO] WDR populate ngrid npe=",I0, I4,I4, F12.4, F12.4)', this_pe, i, j, n_grid(i,j,1), n_grid(i,j,2)
          end do
       end do

       do j = 1, parent_geo%ny
          do i = 1, parent_geo%nx
             ! centered grid version
             p_grid(i, j, 1) = parent_geo%lons(2*i, 2*j)
             p_grid(i, j, 2) = parent_geo%lats(2*i, 2*j)
          end do
       end do

       ! u(npx, npy+1)   
    elseif (position == NORTH) then  ! u wind on D-stagger
       do j = lbound(tile_geo%lats,2), ubound(tile_geo%lats,2)
          do i = lbound(tile_geo%lats,1), ubound(tile_geo%lats,1)
             ! centered grid version
             n_grid(i, j, 1) = tile_geo%lons(i, j)
             n_grid(i, j, 2) = tile_geo%lats(i, j)             
             !if (debug_log) print '("[INFO] WDR populate ngrid_u npe=",I0, I4,I4, F12.4, F12.4)', this_pe, i, j, n_grid(i,j,1), n_grid(i,j,2)
          end do
       end do


       do j = 1, parent_geo%ny
          do i = 1, parent_geo%nx
             ! centered grid version
             p_grid(i, j, 1) = parent_geo%lons(2*i, 2*j-1)
             p_grid(i, j, 2) = parent_geo%lats(2*i, 2*j-1)
          end do
       end do

       ! v(npx+1, npy)   
    elseif (position == EAST) then  ! v wind on D-stagger
       do j = lbound(tile_geo%lats,2), ubound(tile_geo%lats,2)
          do i = lbound(tile_geo%lats,1), ubound(tile_geo%lats,1)
             ! centered grid version
             n_grid(i, j, 1) = tile_geo%lons(i, j)
             n_grid(i, j, 2) = tile_geo%lats(i, j)             
             !if (debug_log) print '("[INFO] WDR populate ngrid_v npe=",I0, I4,I4, F12.4, F12.4)', this_pe, i, j, n_grid(i,j,1), n_grid(i,j,2)
          end do
       end do

       do j = 1, parent_geo%ny
          do i = 1, parent_geo%nx
             ! centered grid version
             p_grid(i, j, 1) = parent_geo%lons(2*i-1, 2*j)
             p_grid(i, j, 2) = parent_geo%lats(2*i-1, 2*j)
          end do
       end do

    end if



  end subroutine assign_n_p_grids





  !! How should weights for interpolation work?
  !! When grid is aligned, it's much easier - weights are related to fractional differences

  subroutine calc_nest_halo_weights(bbox_fine, bbox_coarse, p_grid, n_grid, wt, istart_coarse, jstart_coarse, x_refine, y_refine)
    implicit none

    type(bbox), intent(in)                       :: bbox_coarse, bbox_fine
    real(kind=R_GRID), allocatable, intent(in)   :: p_grid(:,:,:), n_grid(:,:,:)
    real, allocatable, intent(inout)             :: wt(:,:,:)
    integer, intent(in)                          :: istart_coarse, jstart_coarse, x_refine, y_refine


    integer       :: i,j, ic, jc
    real          :: dist1, dist2, dist3, dist4, sum
    logical       :: verbose = .false.
    !logical       :: verbose = .true.

    integer       :: this_pe

    real(kind=R_GRID)  :: pi = 4 * atan(1.0d0)
    real               :: pi180
    real               :: rad2deg, deg2rad

    pi180 = pi / 180.0
    deg2rad = pi / 180.0
    rad2deg = 1.0 / pi180


    this_pe = mpp_pe()


    if ( bbox_coarse%is == 0 .and. bbox_coarse%ie == -1 ) then
       ! Skip this one
       if (debug_log) print '("[INFO] WDR skip calc weights npe=",I0)', this_pe


    else
       if (debug_log) print '("[INFO] WDR run calc weights npe=",I0)', this_pe

       ! Calculate the bounding parent grid points for the nest grid point
       ! Rely on the nest being aligned
       ! code is from $CUBE/tools/fv_grid_tools.F90
       !

       do j = bbox_fine%js, bbox_fine%je
          ! F90 integer division truncates
          jc = jstart_coarse  + (j + y_refine/2 + 1) / y_refine
          do i = bbox_fine%is, bbox_fine%ie
             ic = istart_coarse  + (i + x_refine/2 + 1) / x_refine

             if (verbose) then
                if (debug_log) print '("[INFO] WDR MAP npe=",I0," istart_coarse, jstart_coarse,   ic,if,jc,jf",I3,I3," ",I3,I3,I3,I3)', this_pe, istart_coarse, jstart_coarse,ic,i,jc,j

                if (debug_log) print '("[INFO] WDR LATLON npe=",I0,"  p_grid(",I3,I3,")",F8.2,F8.2, F8.2)', this_pe, ic, jc, rad2deg*p_grid(ic,jc,1)-360.0 , rad2deg*p_grid(ic,jc,2), rad2deg*p_grid(ic,jc,1)
                if (debug_log) print '("[INFO] WDR LATLON npe=",I0,"  nest n_grid(",I3,I3,") ",F8.2,F8.2, F8.2)', this_pe, i, j, rad2deg*n_grid(i,j,1)-360.0, rad2deg*n_grid(i,j,2), rad2deg*n_grid(i,j,1)


                if (debug_log) print '("[INFO] WDR LOC npe=",I0,"  -------------------")', this_pe
                if (debug_log) print '("[INFO] WDR LOC npe=",I0,"  A p_grid(",I3,I3,")",F12.6,F12.6, F12.6)', this_pe, ic, jc, rad2deg*p_grid(ic,jc,1)-360.0, rad2deg*p_grid(ic,jc,2), rad2deg*p_grid(ic,jc,1)
                if (debug_log) print '("[INFO] WDR LOC npe=",I0,"  B p_grid(",I3,I3,")",F12.6,F12.6, F12.6)', this_pe, ic, jc+1, rad2deg*p_grid(ic,jc+1,1)-360.0, rad2deg*p_grid(ic,jc+1,2), rad2deg*p_grid(ic,jc+1,1)
                if (debug_log) print '("[INFO] WDR LOC npe=",I0,"  C p_grid(",I3,I3,")",F12.6,F12.6, F12.6)', this_pe, ic+1, jc+1, rad2deg*p_grid(ic+1,jc+1,1)-360.0, rad2deg*p_grid(ic+1,jc+1,2), rad2deg*p_grid(ic+1,jc+1,1)
                if (debug_log) print '("[INFO] WDR LOC npe=",I0,"  D p_grid(",I3,I3,")",F12.6,F12.6, F12.6)', this_pe, ic+1, jc, rad2deg*p_grid(ic+1,jc,1)-360.0, rad2deg*p_grid(ic+1,jc,2), rad2deg*p_grid(ic+1,jc,1)
                if (debug_log) print '("[INFO] WDR LOC npe=",I0,"  nest n_grid(",I3,I3,") ",F12.6,F12.6, F12.6)', this_pe, i, j, rad2deg*n_grid(i,j,1)-360.0, rad2deg*n_grid(i,j,2), rad2deg*n_grid(i,j,1)
             end if


             ! dist2side_latlon takes points in longitude-latitude coordinates.
             dist1 = dist2side_latlon(p_grid(ic,jc,:),     p_grid(ic,jc+1,:),   n_grid(i,j,:))

             if (verbose) then
                if (debug_log) print '("[INFO] WDR LATLON npe=",I0," dist1=",F9.4," p_grid(",I3,I3,")=",F9.4,F9.4," p_grid(",I3,I3,")=",F9.4,F9.4," n_grid(",I3,I3,")=",F9.4,F9.4)', this_pe, dist1, ic, jc, rad2deg*p_grid(ic,jc,1)-360.0,  rad2deg*p_grid(ic,jc,2),  ic, jc+1, rad2deg*p_grid(ic,jc+1,1)-360.0,  rad2deg*p_grid(ic,jc+1,2), i, j, rad2deg*n_grid(i,j,1)-360.0,  rad2deg*n_grid(i,j,2)
             end if
             dist2 = dist2side_latlon(p_grid(ic,jc+1,:),   p_grid(ic+1,jc+1,:), n_grid(i,j,:))
             dist3 = dist2side_latlon(p_grid(ic+1,jc+1,:), p_grid(ic+1,jc,:),   n_grid(i,j,:))
             dist4 = dist2side_latlon(p_grid(ic,jc,:),     p_grid(ic+1,jc,:),   n_grid(i,j,:))



             !if (debug_log) print '("[INFO] WDR LATLON npe=",I0,"  dists at (",I3,I3,"): dist: ",F12.4, F12.4, F12.4, F12.4)', this_pe, i, j, dist1*RADIUS, dist2*RADIUS, dist3*RADIUS, dist4*RADIUS
             if (verbose) then
                if (debug_log) print '("[INFO] WDR LATLON npe=",I0,"  dists at (",I3,I3,"): dist: ",F12.4, F12.4, F12.4, F12.4)', this_pe, i, j, dist1, dist2, dist3, dist4
             end if

             wt(i,j,1)=dist2*dist3      ! ic,   jc    weight
             wt(i,j,2)=dist3*dist4      ! ic,   jc+1  weight
             wt(i,j,3)=dist4*dist1      ! ic+1, jc+1  weight
             wt(i,j,4)=dist1*dist2      ! ic+1, jc    weight

             sum=wt(i,j,1)+wt(i,j,2)+wt(i,j,3)+wt(i,j,4)
             wt(i,j,:)=wt(i,j,:)/sum

             if (verbose) then
                if (debug_log) print '("[INFO] WDR LATLON npe=",I0,"  sum (",I3,I3,"): ",F12.2,"  wt: ",F12.6, F12.6, F12.6, F12.6)', this_pe, i, j, sum, wt(i,j,1), wt(i,j,2), wt(i,j,3), wt(i,j,4)
             end if

          end do
       end do
    end if

    if (debug_log) print '("[INFO] WDR DONE calc weights npe=",I0)', this_pe

  end subroutine calc_nest_halo_weights

#endif ! MOVING_NEST

end module fv_moving_nest_mod



