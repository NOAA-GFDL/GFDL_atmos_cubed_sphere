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

!>@brief The module 'fv_nesting' is a collection of routines pertaining to grid nesting 
!! \cite harris2013two.

#undef MULTI_GASES

module fv_nesting_mod

! Modules Included:
! <table>
!   <tr>
!     <th>Module Name</th>
!     <th>Functions Included</th>
!  </tr>
!  <tr>
!     <td>boundary_mod</td>
!     <td>update_coarse_grid,nested_grid_BC_send, nested_grid_BC_recv, nested_grid_BC_save_proc
!         nested_grid_BC, nested_grid_BC_apply_intT</td>
!   </tr>
!   <tr>
!     <td>constants_mod</td>
!     <td>grav, pi=>pi_8, radius, hlv, rdgas, cp_air, rvgas, cp_vapor, kappa</td>
!   </tr>
!   <tr>
!     <td>field_manager_mod</td>
!     <td>MODEL_ATMOS</td>
!   </tr>
!   <tr>
!     <td>fv_arrays_mod</td>
!     <td>fv_grid_type, fv_flags_type, fv_atmos_type, fv_nest_type, fv_diag_type, 
!         fv_nest_BC_type_3D,allocate_fv_nest_BC_type,  fv_atmos_type, fv_grid_bounds_type</td>
!   </tr>
!   <tr>
!     <td>fv_diagnostics_mod</td>
!     <td>sphum_ll_fix, range_check</td>
!   </tr>
!   <tr>
!     <td>fv_grid_utils_mod</td>
!     <td>ptop_min, g_sum, cubed_to_latlon, f_p</td>
!   </tr>
!   <tr>
!     <td>fv_mapz_mod</td>
!     <td>mappm, remap_2d</td>
!   </tr>
!   <tr>
!     <td>fv_mp_mod</td>
!     <td>is, ie, js, je, isd, ied, jsd, jed, isc, iec, jsc, jec,is_master,mp_reduce_sum</td>
!   </tr>
!   <tr>
!     <td>fv_restart_mod</td>
!     <td>d2a_setup, d2c_setup</td>
!   </tr>
!   <tr>
!     <td>fv_sg_mod</td>
!     <td>neg_adj3</td>
!   </tr>
!  <tr>
!     <td>fv_timing_mod</td>
!     <td>timing_on, timing_off</td>
!   </tr>
!   <tr>
!     <td>init_hydro_mod</td>
!     <td>p_var</td>
!   </tr>
!   <tr>
!     <td>mpp_mod/td>
!    <td>mpp_sync_self, mpp_sync, mpp_send, mpp_recv, mpp_error, FATAL</td>
!   </tr>
!   <tr>
!     <td>mpp_domains_mod/td>
!     <td>mpp_update_domains, mpp_global_field,mpp_get_data_domain, mpp_get_compute_domain, 
!         mpp_get_global_domain, DGRID_NE, mpp_update_domains, domain2D, mpp_global_sum, 
!         BITWISE_EFP_SUM, BITWISE_EXACT_SUM</td>
!   </tr>
!   <tr>
!    <td>sw_core_mod</td>
!     <td>divergence_corner, divergence_corner_nest</td>
!   </tr>
!   <tr>
!     <td>tracer_manager_mod</td>
!     <td>get_tracer_index</td>
!   </tr>
! </table>

   use mpp_domains_mod,     only: mpp_update_domains
   use mpp_domains_mod,     only: mpp_global_field
   use field_manager_mod,   only: MODEL_ATMOS
   use tracer_manager_mod,  only: get_tracer_index
   use fv_sg_mod,           only: neg_adj3
   use mpp_domains_mod,     only: mpp_get_data_domain, mpp_get_compute_domain, mpp_get_global_domain
   use mpp_domains_mod,     only: AGRID, CGRID_NE, DGRID_NE, mpp_update_domains, domain2D
   use mpp_mod,             only: mpp_sync_self, mpp_sync, mpp_send, mpp_recv, mpp_error, FATAL, mpp_pe, WARNING, NOTE
   use mpp_domains_mod,     only: mpp_global_sum, BITWISE_EFP_SUM, BITWISE_EXACT_SUM
   use boundary_mod,        only: update_coarse_grid
   use boundary_mod,        only: nested_grid_BC_send, nested_grid_BC_recv, nested_grid_BC_save_proc
   use boundary_mod,        only: nested_grid_BC, nested_grid_BC_apply_intT
   use fv_arrays_mod,       only: fv_grid_type, fv_flags_type, fv_atmos_type, fv_nest_type, fv_diag_type, fv_nest_BC_type_3D
   use fv_arrays_mod,       only: allocate_fv_nest_BC_type, fv_atmos_type, fv_grid_bounds_type, deallocate_fv_nest_BC_type
   use fv_grid_utils_mod,   only: ptop_min, g_sum, cubed_to_latlon, f_p
   use init_hydro_mod,      only: p_var
   use constants_mod,       only: grav, pi=>pi_8, radius, hlv, rdgas, cp_air, rvgas, cp_vapor, kappa
   use fv_mapz_mod,         only: mappm, remap_2d
   use fv_timing_mod,       only: timing_on, timing_off
   use fv_mp_mod,           only: is_master
   use fv_mp_mod,           only: mp_reduce_sum, global_nest_domain
   use fv_diagnostics_mod,  only: sphum_ll_fix, range_check
   use sw_core_mod,         only: divergence_corner, divergence_corner_nest
   use time_manager_mod,    only: time_type

implicit none
   logical :: RF_initialized = .false.
   logical :: bad_range
   real, allocatable ::  rf(:), rw(:)
   integer :: kmax=1
   !Arrays for global grid total energy, used for grid nesting
   real, allocatable :: te_2d_coarse(:,:)
   real, allocatable :: dp1_coarse(:,:,:)

   !For nested grid buffers
   !Individual structures are allocated by nested_grid_BC_recv
   type(fv_nest_BC_type_3d) :: u_buf, v_buf, uc_buf, vc_buf, delp_buf, delz_buf, pt_buf, w_buf, divg_buf, pe_u_buf,pe_v_buf,pe_b_buf
   type(fv_nest_BC_type_3d), allocatable:: q_buf(:)
!#ifdef USE_COND
   real, dimension(:,:,:), allocatable, target :: dum_West, dum_East, dum_North, dum_South
!#endif

private
public :: twoway_nesting, setup_nested_grid_BCs, set_physics_BCs

contains
!>@brief The subroutine 'setup_nested_grid_BCs' fetches data from the coarse grid 	
!! to set up  the nested-grid boundary conditions.
 subroutine setup_nested_grid_BCs(npx, npy, npz, zvir, ncnst,     &
                        u, v, w, pt, delp, delz,q, uc, vc, &
#ifdef USE_COND
                        q_con, &
#ifdef MOIST_CAPPA
                        cappa, &
#endif
#endif
                        nested, inline_q, make_nh, ng, &
                        gridstruct, flagstruct, neststruct, &
                        nest_timestep, tracer_nest_timestep, &
                        domain, parent_grid, bd, nwat, ak, bk)

   
    type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(IN) :: zvir

    integer, intent(IN) :: npx, npy, npz
    integer, intent(IN) :: ncnst, ng, nwat
    logical, intent(IN) :: inline_q, make_nh,nested
    real, intent(IN), dimension(npz) :: ak, bk

    real, intent(inout), dimension(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) :: u !< D grid zonal wind (m/s)
    real, intent(inout), dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) :: v !< D grid meridional wind (m/s)
    real, intent(inout) :: w(   bd%isd:        ,bd%jsd:        ,1:)  !<  W (m/s)
    real, intent(inout) :: pt(  bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  !< temperature (K)
    real, intent(inout) :: delp(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  !< pressure thickness (pascal)
    real, intent(inout) :: delz(bd%is:        ,bd%js:        ,1:)  !< height thickness (m)
    real, intent(inout) :: q(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz, ncnst) !< specific humidity and constituents
    real, intent(inout) :: uc(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) !< (uc,vc) mostly used as the C grid winds
    real, intent(inout) :: vc(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz)
#ifdef USE_COND
    real, intent(inout) :: q_con(  bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)
#ifdef MOIST_CAPPA
    real, intent(inout) :: cappa(  bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)
#endif
#endif
    integer, intent(INOUT) :: nest_timestep, tracer_nest_timestep
    type(fv_atmos_type), pointer, intent(IN) :: parent_grid

    type(fv_grid_type), intent(INOUT) :: gridstruct
    type(fv_flags_type), intent(INOUT) :: flagstruct
    type(fv_nest_type), intent(INOUT), target :: neststruct
    type(domain2d), intent(INOUT) :: domain
    real :: divg(bd%isd:bd%ied+1,bd%jsd:bd%jed+1, npz)
    real :: ua(bd%isd:bd%ied,bd%jsd:bd%jed)
    real :: va(bd%isd:bd%ied,bd%jsd:bd%jed)
    real :: pe_ustag(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz+1)
    real :: pe_vstag(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz+1)
    real :: pe_bstag(bd%isd:bd%ied+1,bd%jsd:bd%jed+1,npz+1)
    real, parameter :: a13 = 1./3.

    integer :: i,j,k,n,p, sphum, npz_coarse, nnest
    logical :: do_pd

    type(fv_nest_BC_type_3d) :: delp_lag_BC, lag_BC, pe_lag_BC, pe_eul_BC
    type(fv_nest_BC_type_3d) :: lag_u_BC, pe_u_lag_BC, pe_u_eul_BC
    type(fv_nest_BC_type_3d) :: lag_v_BC, pe_v_lag_BC, pe_v_eul_BC
    type(fv_nest_BC_type_3d) :: lag_b_BC, pe_b_lag_BC, pe_b_eul_BC
    
    !local pointers
    logical, pointer :: child_grids(:)
    
    integer :: is,  ie,  js,  je
    integer :: isd, ied, jsd, jed
    
    is  = bd%is
    ie  = bd%ie
    js  = bd%js
    je  = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    
    child_grids => neststruct%child_grids
    
    !IF nested, set up nested grid BCs for time-interpolation
    !(actually applying the BCs is done in dyn_core)

    !For multiple grids: Each grid has ONE parent but potentially MULTIPLE nests

    nest_timestep = 0
    if (.not. inline_q) tracer_nest_timestep = 0


    if (neststruct%nested .and. (.not. (neststruct%first_step) .or. make_nh) ) then
       do_pd = .true.
       call set_BCs_t0(ncnst, flagstruct%hydrostatic, neststruct) 
    else
       !On first timestep the t0 BCs are not initialized and may contain garbage
       do_pd = .false.
    end if

    !compute uc/vc for nested-grid BCs
    !!! CLEANUP: if we compute uc/vc here we don't need to do on the first call of c_sw, right?
    if (ANY(neststruct%child_grids)) then
       call timing_on('COMM_TOTAL')
       !!! CLEANUP: could we make this a non-blocking operation?
       !!! Is this needed? it is on the initialization step.
       call mpp_update_domains(delp, domain) !This is needed to make sure delp is updated for pe calculations
       call mpp_update_domains(u, v, &       
            domain, gridtype=DGRID_NE, complete=.true.)
       call timing_off('COMM_TOTAL')
!$OMP parallel do default(none) shared(isd,jsd,ied,jed,is,ie,js,je,npx,npy,npz, &
!$OMP       gridstruct,flagstruct,bd,u,v,uc,vc,nested,divg) &
!$OMP       private(ua,va)
       do k=1,npz
          call d2c_setup(u(isd,jsd,k),  v(isd,jsd,k),   &
               ua, va, &
               uc(isd,jsd,k), vc(isd,jsd,k), flagstruct%nord>0, &
               isd,ied,jsd,jed, is,ie,js,je, npx,npy, &
               gridstruct%grid_type, gridstruct%bounded_domain, &
               gridstruct%se_corner, gridstruct%sw_corner, &
               gridstruct%ne_corner, gridstruct%nw_corner, &
               gridstruct%rsin_u, gridstruct%rsin_v, &
               gridstruct%cosa_s, gridstruct%rsin2 )
          if (nested) then
             call divergence_corner_nest(u(isd,jsd,k), v(isd,jsd,k), ua, va, divg(isd,jsd,k), gridstruct, flagstruct, bd)
          else
             call divergence_corner(u(isd,jsd,k), v(isd,jsd,k), ua, va, divg(isd,jsd,k), gridstruct, flagstruct, bd)
          endif
       end do       
    endif

    nnest = flagstruct%grid_number - 1

!! Nested grid: receive from parent grid (Lagrangian coordinate, npz_coarse)
    if (neststruct%nested) then

       npz_coarse = neststruct%parent_grid%npz

       if (.not. allocated(q_buf)) then
          allocate(q_buf(ncnst))
       endif

       call nested_grid_BC_recv(global_nest_domain, 0, 0,  npz_coarse, bd, &
            delp_buf, nnest)
       do n=1,ncnst
          call nested_grid_BC_recv(global_nest_domain, 0, 0, npz_coarse, bd, &
               q_buf(n), nnest)
       enddo
#ifndef SW_DYNAMICS
       call nested_grid_BC_recv(global_nest_domain, 0, 0, npz_coarse, bd, &
            pt_buf, nnest)

       if (.not. flagstruct%hydrostatic) then
          call nested_grid_BC_recv(global_nest_domain, 0, 0,  npz_coarse, bd, &
               w_buf, nnest)
          call nested_grid_BC_recv(global_nest_domain, 0, 0,  npz_coarse, bd, &
               delz_buf, nnest)
       endif
#endif
       if (neststruct%do_remap_BC(flagstruct%grid_number)) then
          call nested_grid_BC_recv(global_nest_domain, npz_coarse+1, bd, &
               pe_u_buf, pe_v_buf, nnest, gridtype=DGRID_NE)
          call nested_grid_BC_recv(global_nest_domain, 1, 1,  npz_coarse+1, bd, &
               pe_b_buf, nnest)
       endif

       call nested_grid_BC_recv(global_nest_domain, npz_coarse, bd, &
               u_buf, v_buf, nnest, gridtype=DGRID_NE)
       call nested_grid_BC_recv(global_nest_domain, npz_coarse, bd, &
               uc_buf, vc_buf, nnest, gridtype=CGRID_NE)
       call nested_grid_BC_recv(global_nest_domain, 1, 1,  npz_coarse, bd, &
            divg_buf, nnest)
    endif


!! Coarse grid: send to child grids (Lagrangian coordinate, npz_coarse)

    do p=1,size(child_grids)
       if (child_grids(p)) then
          call nested_grid_BC_send(delp, global_nest_domain, 0, 0, p-1)
          do n=1,ncnst
             call nested_grid_BC_send(q(:,:,:,n), global_nest_domain, 0, 0, p-1)
          enddo
#ifndef SW_DYNAMICS
          call nested_grid_BC_send(pt, global_nest_domain, 0, 0, p-1)

          if (.not. flagstruct%hydrostatic) then
             call nested_grid_BC_send(w, global_nest_domain, 0, 0, p-1)
             call nested_grid_BC_send(delz, global_nest_domain, 0, 0, p-1) 
          endif          
#endif

          if (neststruct%do_remap_BC(p)) then 

          !Compute and send staggered pressure
             !u points
!$OMP parallel do default(none) shared(ak,pe_ustag,delp, &
!$OMP                                  is,ie,js,je,npz)
             do j=js,je+1
             do i=is,ie
                pe_ustag(i,j,1) = ak(1)
             enddo
             do k=1,npz
             do i=is,ie
                pe_ustag(i,j,k+1) = pe_ustag(i,j,k) + 0.5*(delp(i,j,k)+delp(i,j-1,k)) 
             enddo
             enddo
             enddo

             !v points
!$OMP parallel do default(none) shared(ak,pe_vstag,delp, &
!$OMP                                  is,ie,js,je,npz)
             do j=js,je
             do i=is,ie+1
                pe_vstag(i,j,1) = ak(1)
             enddo
             do k=1,npz
             do i=is,ie+1
                pe_vstag(i,j,k+1) = pe_vstag(i,j,k) + 0.5*(delp(i,j,k)+delp(i-1,j,k)) 
             enddo
             enddo
             enddo
             call nested_grid_BC_send(pe_ustag, pe_vstag, global_nest_domain, p-1, gridtype=DGRID_NE)

             !b points
!$OMP parallel do default(none) shared(ak,pe_bstag,delp, &
!$OMP                                  is,ie,js,je,npz)
             do j=js,je+1
             do i=is,ie+1
                pe_bstag(i,j,1) = ak(1)
             enddo
             enddo
             !Sets up so 3-point average is automatically done at the corner
             if (is == 1 .and. js == 1) then
                do k=1,npz
                   delp(0,0,k) = a13*(delp(1,1,k) + delp(0,1,k) + delp(1,0,k))
                enddo
             endif
             if (ie == npx-1 .and. js == 1) then
                do k=1,npz
                   delp(npx,0,k) = a13*(delp(npx-1,1,k) + delp(npx,1,k) + delp(npx-1,0,k))
                enddo
             endif
             if (is == 1 .and. je == npy-1) then
                do k=1,npz
                   delp(0,npy,k) = a13*(delp(1,npy-1,k) + delp(0,npy-1,k) + delp(1,npy,k))
                enddo
             endif
             if (ie == npx-1 .and. je == npy-1) then
                do k=1,npz
                   delp(npx,npy,k) = a13*(delp(npx-1,npy-1,k) + delp(npx,npy-1,k) + delp(npx-1,npy,k))
                enddo
             endif
 
!$OMP parallel do default(none) shared(ak,pe_bstag,delp, &
!$OMP                                  is,ie,js,je,npz)
             do j=js,je+1
             do k=1,npz
             do i=is,ie+1
                pe_bstag(i,j,k+1) = pe_bstag(i,j,k) + & 
                     0.25*(delp(i,j,k)+delp(i-1,j,k)+delp(i,j-1,k)+delp(i-1,j-1,k))
             enddo
             enddo
             enddo
             call nested_grid_BC_send(pe_bstag, global_nest_domain, 1, 1, p-1)

          endif

          call nested_grid_BC_send(u, v, global_nest_domain, p-1, gridtype=DGRID_NE)
          call nested_grid_BC_send(uc, vc, global_nest_domain, p-1, gridtype=CGRID_NE)
          call nested_grid_BC_send(divg, global_nest_domain, 1, 1, p-1)
       endif
    enddo
    
    !Nested grid: do computations
    ! Lag: coarse grid, npz_coarse, lagrangian coordinate---receive and use save_proc to copy into lag_BCs
    ! Eul: nested grid, npz, Eulerian (reference) coordinate
    ! Remapping from Lag to Eul
    if (nested) then

       if (neststruct%do_remap_BC(flagstruct%grid_number)) then

          call allocate_fv_nest_BC_type(delp_lag_BC,is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz_coarse,ng,0,0,0,.false.)
          call allocate_fv_nest_BC_type(lag_BC,is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz_coarse,ng,0,0,0,.false.)
          call allocate_fv_nest_BC_type(pe_lag_BC,is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz_coarse+1,ng,0,0,0,.false.)
          call allocate_fv_nest_BC_type(pe_eul_BC,is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz+1,ng,0,0,0,.false.)
          
          call nested_grid_BC_save_proc(global_nest_domain, &
               neststruct%ind_h, neststruct%wt_h, 0, 0,  npx, npy, npz_coarse, bd, &
               delp_lag_BC, delp_buf, pd_in=do_pd)
          !The incoming delp is on the coarse grid's lagrangian coordinate. Re-create the reference coordinate
          call setup_eul_delp_BC(delp_lag_BC, neststruct%delp_BC, pe_lag_BC, pe_eul_BC, ak, bk, npx, npy, npz, npz_coarse, parent_grid%ptop, bd)

       else
          call nested_grid_BC_save_proc(global_nest_domain, &
               neststruct%ind_h, neststruct%wt_h, 0, 0,  npx, npy, npz_coarse, bd, &
               neststruct%delp_BC, delp_buf, pd_in=do_pd)
       endif

!!$       do n=1,ncnst
!!$          call nested_grid_BC_save_proc(global_nest_domain, &
!!$               neststruct%ind_h, neststruct%wt_h, 0, 0, npx,  npy,  npz_coarse, bd, &
!!$               lag_BC, q_buf(n), pd_in=do_pd)
!!$          !This remapping appears to have some trouble with rounding error random noise
!!$          call remap_BC(pe_lag_BC, pe_eul_BC, lag_BC, neststruct%q_BC(n), npx, npy, npz, npz_coarse, bd, 0, 0, 0, flagstruct%kord_tr, 'q')
!!$       enddo
#ifndef SW_DYNAMICS
       if (neststruct%do_remap_BC(flagstruct%grid_number)) then

          call nested_grid_BC_save_proc(global_nest_domain, &
               neststruct%ind_h, neststruct%wt_h, 0, 0, npx,  npy,  npz_coarse, bd, &
               lag_BC, pt_buf) 
          !NOTE: need to remap using peln, not pe
          call remap_BC(pe_lag_BC, pe_eul_BC, lag_BC, neststruct%pt_BC, npx, npy, npz, npz_coarse, bd, 0, 0, 1, abs(flagstruct%kord_tm), 'pt', do_log_pe=.true.)

       else
          call nested_grid_BC_save_proc(global_nest_domain, &
               neststruct%ind_h, neststruct%wt_h, 0, 0, npx,  npy,  npz_coarse, bd, &
               neststruct%pt_BC, pt_buf) 
       endif
          

       !For whatever reason moving the calls for q BC remapping here avoids problems with cross-restart reproducibility. 
       if (neststruct%do_remap_BC(flagstruct%grid_number)) then
          do n=1,ncnst
             call nested_grid_BC_save_proc(global_nest_domain, &
                  neststruct%ind_h, neststruct%wt_h, 0, 0, npx,  npy,  npz_coarse, bd, &
                  lag_BC, q_buf(n), pd_in=do_pd)
             call remap_BC(pe_lag_BC, pe_eul_BC, lag_BC, neststruct%q_BC(n), npx, npy, npz, npz_coarse, bd, 0, 0, 0, flagstruct%kord_tr, 'q2')
          enddo
       else
          do n=1,ncnst
             call nested_grid_BC_save_proc(global_nest_domain, &
                  neststruct%ind_h, neststruct%wt_h, 0, 0, npx,  npy,  npz_coarse, bd, &
                  neststruct%q_BC(n), q_buf(n), pd_in=do_pd)
          enddo
       endif

       sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
       if (flagstruct%hydrostatic) then
          call setup_pt_BC(neststruct%pt_BC, pe_eul_BC, neststruct%q_BC(sphum), npx, npy, npz, zvir, bd)
       else
          if (neststruct%do_remap_BC(flagstruct%grid_number)) then

             call nested_grid_BC_save_proc(global_nest_domain, &
                  neststruct%ind_h, neststruct%wt_h, 0, 0,  npx,  npy,  npz_coarse, bd, &
                  lag_BC, w_buf)
             call remap_BC(pe_lag_BC, pe_eul_BC, lag_BC, neststruct%w_BC, npx, npy, npz, npz_coarse, bd, 0, 0, -1, flagstruct%kord_wz, 'w')
             call nested_grid_BC_save_proc(global_nest_domain, &
                  neststruct%ind_h, neststruct%wt_h, 0, 0,  npx,  npy, npz_coarse, bd, &
                  lag_BC, delz_buf) !Need a negative-definite method? 
             call remap_delz_BC(pe_lag_BC, pe_eul_BC, delp_lag_BC, lag_BC, neststruct%delp_BC, neststruct%delz_BC, npx, npy, npz, npz_coarse, bd, 0, 0, 1, flagstruct%kord_wz)
          
          else
             call nested_grid_BC_save_proc(global_nest_domain, &
                  neststruct%ind_h, neststruct%wt_h, 0, 0,  npx,  npy,  npz_coarse, bd, &
                  neststruct%w_BC, w_buf)
             call nested_grid_BC_save_proc(global_nest_domain, &
                  neststruct%ind_h, neststruct%wt_h, 0, 0,  npx,  npy, npz_coarse, bd, &
                  neststruct%delz_BC, delz_buf) !Need a negative-definite method? 
          endif

          call setup_pt_NH_BC(neststruct%pt_BC, neststruct%delp_BC, neststruct%delz_BC, &
               neststruct%q_BC(sphum), neststruct%q_BC, ncnst, &
#ifdef USE_COND
               neststruct%q_con_BC, &
#ifdef MOIST_CAPPA
               neststruct%cappa_BC, &
#endif
#endif
               npx, npy, npz, zvir, bd)
       endif

#endif

       !!!NOTE: The following require remapping on STAGGERED grids, which requires additional pressure data

       if (neststruct%do_remap_BC(flagstruct%grid_number)) then


          call allocate_fv_nest_BC_type(pe_u_lag_BC,is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz_coarse+1,ng,0,0,1,.false.)
          call allocate_fv_nest_BC_type(pe_u_eul_BC,is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz+1       ,ng,0,0,1,.false.)
          call allocate_fv_nest_BC_type(lag_u_BC,   is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz_coarse  ,ng,0,0,1,.false.)
          call allocate_fv_nest_BC_type(pe_v_lag_BC,is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz_coarse+1,ng,0,1,0,.false.)
          call allocate_fv_nest_BC_type(pe_v_eul_BC,is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz+1       ,ng,0,1,0,.false.)
          call allocate_fv_nest_BC_type(lag_v_BC,   is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz_coarse  ,ng,0,1,0,.false.)
          call allocate_fv_nest_BC_type(pe_b_lag_BC,is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz_coarse+1,ng,0,1,1,.false.)
          call allocate_fv_nest_BC_type(pe_b_eul_BC,is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz+1       ,ng,0,1,1,.false.)
          call allocate_fv_nest_BC_type(lag_b_BC,   is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz_coarse  ,ng,0,1,1,.false.)

          call nested_grid_BC_save_proc(global_nest_domain, &
               neststruct%ind_u, neststruct%wt_u, 0, 1,  npx,  npy,  npz_coarse+1, bd, &
               pe_u_lag_BC, pe_u_buf)
          call setup_eul_pe_BC(pe_u_lag_BC, pe_u_eul_BC, ak, bk, npx, npy, npz, npz_coarse, 0, 1, bd)
          call nested_grid_BC_save_proc(global_nest_domain, &
               neststruct%ind_v, neststruct%wt_v, 1, 0,  npx,  npy,  npz_coarse+1, bd, &
               pe_v_lag_BC, pe_v_buf)
          call setup_eul_pe_BC(pe_v_lag_BC, pe_v_eul_BC, ak, bk, npx, npy, npz, npz_coarse, 1, 0, bd)
          call nested_grid_BC_save_proc(global_nest_domain, &
               neststruct%ind_b, neststruct%wt_b, 1, 1,  npx,  npy,  npz_coarse+1, bd, &
               pe_b_lag_BC, pe_b_buf)
          call setup_eul_pe_BC(pe_b_lag_BC, pe_b_eul_BC, ak, bk, npx, npy, npz, npz_coarse, 1, 1, bd)

          call nested_grid_BC_save_proc(global_nest_domain, &
               neststruct%ind_u, neststruct%wt_u, 0, 1,  npx,  npy,  npz_coarse, bd, &
               lag_u_BC, u_buf)
          call remap_BC(pe_u_lag_BC, pe_u_eul_BC, lag_u_BC, neststruct%u_BC, npx, npy, npz, npz_coarse, bd, 0, 1, -1, flagstruct%kord_mt, 'u')
          call nested_grid_BC_save_proc(global_nest_domain, &
               neststruct%ind_u, neststruct%wt_u, 0, 1,  npx,  npy,  npz_coarse, bd, &
               lag_u_BC, vc_buf)
          call remap_BC(pe_u_lag_BC, pe_u_eul_BC, lag_u_BC, neststruct%vc_BC, npx, npy, npz, npz_coarse, bd, 0, 1, -1, flagstruct%kord_mt, 'vc')
          call nested_grid_BC_save_proc(global_nest_domain, &
               neststruct%ind_v, neststruct%wt_v, 1, 0,  npx,  npy,  npz_coarse, bd, &
               lag_v_BC, v_buf)
          call remap_BC(pe_v_lag_BC, pe_v_eul_BC, lag_v_BC, neststruct%v_BC, npx, npy, npz, npz_coarse, bd, 1, 0, -1, flagstruct%kord_mt, 'v')
          call nested_grid_BC_save_proc(global_nest_domain, &
               neststruct%ind_v, neststruct%wt_v, 1, 0,  npx,  npy,  npz_coarse, bd, &
               lag_v_BC, uc_buf)
          call remap_BC(pe_v_lag_BC, pe_v_eul_BC, lag_v_BC, neststruct%uc_BC, npx, npy, npz, npz_coarse, bd, 1, 0, -1, flagstruct%kord_mt, 'uc')
          call nested_grid_BC_save_proc(global_nest_domain, &
               neststruct%ind_b, neststruct%wt_b, 1, 1,  npx,  npy,  npz_coarse, bd, &
               lag_b_BC, divg_buf)
          call remap_BC(pe_b_lag_BC, pe_b_eul_BC, lag_b_BC, neststruct%divg_BC, npx, npy, npz, npz_coarse, bd, 1, 1, -1, flagstruct%kord_mt, 'divg')

          call deallocate_fv_nest_BC_type(delp_lag_BC)
          call deallocate_fv_nest_BC_type(lag_BC)
          call deallocate_fv_nest_BC_type(pe_lag_BC)
          call deallocate_fv_nest_BC_type(pe_eul_BC)

          call deallocate_fv_nest_BC_type(pe_u_lag_BC)
          call deallocate_fv_nest_BC_type(pe_u_eul_BC)
          call deallocate_fv_nest_BC_type(lag_u_BC)
          call deallocate_fv_nest_BC_type(pe_v_lag_BC)
          call deallocate_fv_nest_BC_type(pe_v_eul_BC)
          call deallocate_fv_nest_BC_type(lag_v_BC)
          call deallocate_fv_nest_BC_type(pe_b_lag_BC)
          call deallocate_fv_nest_BC_type(pe_b_eul_BC)
          call deallocate_fv_nest_BC_type(lag_b_BC)

       else

          call nested_grid_BC_save_proc(global_nest_domain, &
               neststruct%ind_u, neststruct%wt_u, 0, 1,  npx,  npy,  npz_coarse, bd, &
               neststruct%u_BC, u_buf)
          call nested_grid_BC_save_proc(global_nest_domain, &
               neststruct%ind_u, neststruct%wt_u, 0, 1,  npx,  npy,  npz_coarse, bd, &
               neststruct%vc_BC, vc_buf)
          call nested_grid_BC_save_proc(global_nest_domain, &
               neststruct%ind_v, neststruct%wt_v, 1, 0,  npx,  npy,  npz_coarse, bd, &
               neststruct%v_BC, v_buf)
          call nested_grid_BC_save_proc(global_nest_domain, &
               neststruct%ind_v, neststruct%wt_v, 1, 0,  npx,  npy,  npz_coarse, bd, &
               neststruct%uc_BC, uc_buf)
          call nested_grid_BC_save_proc(global_nest_domain, &
               neststruct%ind_b, neststruct%wt_b, 1, 1,  npx,  npy,  npz_coarse, bd, &
               neststruct%divg_BC, divg_buf)

       endif

       !Correct halo values have now been set up for BCs; we can go ahead and apply them too
       call nested_grid_BC_apply_intT(delp, &
            0, 0, npx, npy, npz, bd, 1., 1., &
            neststruct%delp_BC, bctype=neststruct%nestbctype  )
       do n=1,ncnst
          call nested_grid_BC_apply_intT(q(:,:,:,n), &
               0, 0, npx, npy, npz, bd, 1., 1., &
               neststruct%q_BC(n), bctype=neststruct%nestbctype  )          
       enddo
#ifndef SW_DYNAMICS
       call nested_grid_BC_apply_intT(pt, &
            0, 0, npx, npy, npz, bd, 1., 1., &
            neststruct%pt_BC, bctype=neststruct%nestbctype  )
       if (.not. flagstruct%hydrostatic) then
          call nested_grid_BC_apply_intT(w, &
               0, 0, npx, npy, npz, bd, 1., 1., &
               neststruct%w_BC, bctype=neststruct%nestbctype  )
          !Removed halo from delz --- BCs now directly applied in nh_BC --- lmh june 2018
!!$          call nested_grid_BC_apply_intT(delz, &
!!$               0, 0, npx, npy, npz, bd, 1., 1., &
!!$               neststruct%delz_BC, bctype=neststruct%nestbctype  )
       endif
#ifdef USE_COND
       call nested_grid_BC_apply_intT(q_con, &
            0, 0, npx, npy, npz, bd, 1., 1., &
            neststruct%q_con_BC, bctype=neststruct%nestbctype  )            
#ifdef MOIST_CAPPA
       call nested_grid_BC_apply_intT(cappa, &
            0, 0, npx, npy, npz, bd, 1., 1., &
            neststruct%cappa_BC, bctype=neststruct%nestbctype  )            
#endif
#endif
#endif
       call nested_grid_BC_apply_intT(u, &
            0, 1, npx, npy, npz, bd, 1., 1., &
            neststruct%u_BC, bctype=neststruct%nestbctype  )            
       call nested_grid_BC_apply_intT(vc, &
            0, 1, npx, npy, npz, bd, 1., 1., &
            neststruct%vc_BC, bctype=neststruct%nestbctype  )            
       call nested_grid_BC_apply_intT(v, &
            1, 0, npx, npy, npz, bd, 1., 1., &
            neststruct%v_BC, bctype=neststruct%nestbctype  )            
       call nested_grid_BC_apply_intT(uc, &
            1, 0, npx, npy, npz, bd, 1., 1., &
            neststruct%uc_BC, bctype=neststruct%nestbctype  )            
       !!!NOTE: Divg not available here but not needed
       !!! until dyn_core anyway.
!!$       call nested_grid_BC_apply_intT(divg, &
!!$            1, 1, npx, npy, npz, bd, 1., 1., &
!!$            neststruct%divg_BC, bctype=neststruct%nestbctype  )            

       !Update domains needed for Rayleigh damping
       if (.not. flagstruct%hydrostatic) call mpp_update_domains(w, domain) 
       call mpp_update_domains(u, v, domain, gridtype=DGRID_NE, complete=.true.)

    endif

    if (neststruct%first_step) then
       if (neststruct%nested) call set_BCs_t0(ncnst, flagstruct%hydrostatic, neststruct)
       neststruct%first_step = .false.
       if (.not. flagstruct%hydrostatic) flagstruct%make_nh= .false. 
    else if (flagstruct%make_nh) then
       if (neststruct%nested) call set_NH_BCs_t0(neststruct)
       flagstruct%make_nh= .false. 
    endif

    !Unnecessary?
!!$    if ( neststruct%nested .and. .not. neststruct%divg_BC%initialized) then
!!$       neststruct%divg_BC%east_t0  = neststruct%divg_BC%east_t1
!!$       neststruct%divg_BC%west_t0  = neststruct%divg_BC%west_t1
!!$       neststruct%divg_BC%north_t0 = neststruct%divg_BC%north_t1
!!$       neststruct%divg_BC%south_t0 = neststruct%divg_BC%south_t1 
!!$       neststruct%divg_BC%initialized = .true.
!!$    endif


    call mpp_sync_self

 end subroutine setup_nested_grid_BCs

 subroutine set_physics_BCs(ps, u_dt, v_dt, flagstruct, gridstruct, neststruct, npx, npy, npz, ng, ak, bk, bd)

   type(fv_grid_bounds_type), intent(IN) :: bd
   type(fv_flags_type), intent(IN) :: flagstruct
   type(fv_nest_type), intent(INOUT), target :: neststruct
   type(fv_grid_type) :: gridstruct
   integer, intent(IN) :: npx, npy, npz, ng
   real, intent(IN), dimension(npz+1) :: ak, bk
   real, intent(INOUT), dimension(bd%isd:bd%ied,bd%jsd:bd%jed) :: ps
   real, intent(INOUT), dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz) :: u_dt, v_dt
   real, dimension(1,1) :: parent_ps ! dummy variable for nesting
   type(fv_nest_BC_type_3d) :: u_dt_buf, v_dt_buf, pe_src_BC, pe_dst_BC!, var_BC

   integer :: n, npz_coarse, nnest
   integer :: is,  ie,  js,  je
   integer :: isd, ied, jsd, jed
   real    :: dum(1,1,1)

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed

   nnest = flagstruct%grid_number - 1

   if (gridstruct%nested) then
      
      if (neststruct%do_remap_BC(flagstruct%grid_number)) then

         npz_coarse = neststruct%parent_grid%npz

         !Both nested and coarse grids assumed on Eulerian coordinates at this point
         !Only need to fetch ps to form pressure levels
         !Note also u_dt and v_dt are unstaggered
         call nested_grid_BC(ps, parent_ps, global_nest_domain, neststruct%ind_h, neststruct%wt_h, 0, 0, &
              npx, npy, bd, 1, npx-1, 1, npy-1)
         call nested_grid_BC_recv(global_nest_domain, npz_coarse, bd, u_dt_buf, v_dt_buf, nnest, gridtype=AGRID)

         call allocate_fv_nest_BC_type(pe_src_BC, is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz_coarse+1,ng,0,0,0,.false.)
         call allocate_fv_nest_BC_type(pe_dst_BC, is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz+1,ng,0,0,0,.false.)

         call copy_ps_BC(ps, pe_src_BC, npx, npy, npz_coarse, 0, 0, bd)
         call setup_eul_pe_BC(pe_src_BC, pe_dst_BC, ak, bk, npx, npy, npz, npz_coarse, 0, 0, bd, &
              make_src_in=.true., ak_src=neststruct%parent_grid%ak, bk_src=neststruct%parent_grid%bk)

         !Note that iv=-1 is used for remapping winds, which sets the lower reconstructed values to 0 if
         ! there is a 2dx signal. Is this the best for **tendencies** though?? Probably not---so iv=1 here
         call set_BC_direct( pe_src_BC, pe_dst_BC, u_dt_buf, u_dt, neststruct, npx, npy, npz, npz_coarse, ng, bd, 0, 0, 1, flagstruct%kord_mt)
         call set_BC_direct( pe_src_BC, pe_dst_BC, v_dt_buf, v_dt, neststruct, npx, npy, npz, npz_coarse, ng, bd, 0, 0, 1, flagstruct%kord_mt)

         call deallocate_fv_nest_BC_type(pe_src_BC)
         call deallocate_fv_nest_BC_type(pe_dst_BC)

      else
         call nested_grid_BC(u_dt, v_dt, dum, dum, global_nest_domain, neststruct%ind_h, neststruct%ind_h, &
              neststruct%wt_h, neststruct%wt_h, 0, 0, 0, 0, npx, npy, npz, bd, 1, npx-1, 1, npy-1, nnest, gridtype=AGRID)
      endif

   endif
   do n=1,size(neststruct%child_grids)
      if (neststruct%child_grids(n)) then
         if (neststruct%do_remap_BC(n)) &
              call nested_grid_BC(ps, global_nest_domain, 0, 0, n-1)
         call nested_grid_BC_send(u_dt, v_dt, global_nest_domain, n-1, gridtype=AGRID)
      endif
   enddo


 end subroutine set_physics_BCs

 subroutine set_BC_direct( pe_src_BC, pe_dst_BC, buf, var, neststruct, npx, npy, npz, npz_coarse, ng, bd, istag, jstag, iv, kord)

   type(fv_grid_bounds_type), intent(IN) :: bd
   type(fv_nest_type), intent(INOUT) :: neststruct
   integer, intent(IN) :: npx, npy, npz, npz_coarse, ng, istag, jstag, iv, kord
   real, intent(INOUT), dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,npz) :: var
   type(fv_nest_BC_type_3d), intent(INOUT) :: buf, pe_src_BC, pe_dst_BC
   type(fv_nest_BC_type_3d) :: var_BC
   

   call allocate_fv_nest_BC_type(var_BC,bd%is,bd%ie,bd%js,bd%je,bd%isd,bd%ied,bd%jsd,bd%jed,npx,npy,npz_coarse,ng,0,istag,jstag,.false.)

   call nested_grid_BC_save_proc(global_nest_domain, neststruct%ind_h, neststruct%wt_h, istag, jstag, &
        npx, npy, npz_coarse, bd, var_BC, buf)
   call remap_BC_direct(pe_src_BC, pe_dst_BC, var_BC, var, npx, npy, npz, npz_coarse, bd, istag, jstag, iv, kord)
   
   call deallocate_fv_nest_BC_type(var_BC)


 end subroutine set_BC_direct

 subroutine setup_pt_BC(pt_BC, pe_eul_BC, sphum_BC, npx, npy, npz, zvir, bd)

   type(fv_grid_bounds_type), intent(IN)   :: bd
   type(fv_nest_BC_type_3d), intent(IN)    :: pe_eul_BC, sphum_BC
   type(fv_nest_BC_type_3d), intent(INOUT) :: pt_BC
   integer, intent(IN) :: npx, npy, npz
   real, intent(IN) :: zvir

   integer :: istart, iend

   integer :: is,  ie,  js,  je
   integer :: isd, ied, jsd, jed

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed
   
   if (is == 1) then
      call setup_pt_BC_k(pt_BC%west_t1, sphum_BC%west_t1, pe_eul_BC%west_t1, zvir, isd, ied, isd, 0, jsd, jed, npz)
   end if

   if (js == 1) then
      if (is == 1) then
         istart = is
      else
         istart = isd
      end if
      if (ie == npx-1) then
         iend = ie
      else
         iend = ied
      end if

      call setup_pt_BC_k(pt_BC%south_t1, sphum_BC%south_t1, pe_eul_BC%south_t1, zvir, isd, ied, istart, iend, jsd, 0, npz)
   end if


   if (ie == npx-1) then
      call setup_pt_BC_k(pt_BC%east_t1, sphum_BC%east_t1, pe_eul_BC%east_t1, zvir, isd, ied, npx, ied, jsd, jed, npz)
   end if

   if (je == npy-1) then
      if (is == 1) then
         istart = is
      else
         istart = isd
      end if
      if (ie == npx-1) then
         iend = ie
      else
         iend = ied
      end if

      call setup_pt_BC_k(pt_BC%north_t1, sphum_BC%north_t1, pe_eul_BC%north_t1, zvir, isd, ied, istart, iend, npy, jed, npz)
   end if
   
 end subroutine setup_pt_BC


!!!! A NOTE ON NOMENCLATURE
!!!! Originally the BC arrays were bounded by isd and ied in the i-direction.
!!!!   However these were NOT intended to delineate the dimensions of the data domain
!!!!   but instead were of the BC arrays. This is confusing especially in other locations
!!!!   where BCs and data arrays are both present.
 subroutine setup_pt_BC_k(ptBC, sphumBC, peBC, zvir, isd_BC, ied_BC, istart, iend, jstart, jend, npz)

   integer, intent(IN) :: isd_BC, ied_BC, istart, iend, jstart, jend, npz
   real,    intent(IN) :: zvir
   real, intent(INOUT), dimension(isd_BC:ied_BC,jstart:jend,npz) :: ptBC
   real, intent(IN),    dimension(isd_BC:ied_BC,jstart:jend,npz) :: sphumBC
   real, intent(IN),    dimension(isd_BC:ied_BC,jstart:jend,npz+1) :: peBC

   integer :: i,j,k
   real :: pealn, pebln, rpkz

!Assumes dry kappa
!$OMP parallel do default(none) shared(peBC,ptBC,zvir,sphumBC, &
!$OMP                                  istart,iend,jstart,jend,npz) &
!$OMP                           private(pealn,pebln,rpkz)
   do k=1,npz
   do j=jstart,jend
   do i=istart,iend
      pealn = log(peBC(i,j,k))
      pebln = log(peBC(i,j,k+1))

      rpkz =  kappa*(pebln - pealn)/(exp(kappa*pebln)-exp(kappa*pealn) )

      ptBC(i,j,k) = ptBC(i,j,k)*rpkz * &
           (1.+zvir*sphumBC(i,j,k))
   enddo
   enddo
   enddo

 end subroutine setup_pt_BC_k

 subroutine setup_eul_delp_BC(delp_lag_BC, delp_eul_BC, pe_lag_BC, pe_eul_BC, ak_dst, bk_dst, npx, npy, npz, npz_coarse, ptop_src, bd)

   type(fv_grid_bounds_type), intent(IN) :: bd
   type(fv_nest_BC_type_3d), intent(INOUT), target :: delp_lag_BC
   type(fv_nest_BC_type_3d), intent(INOUT), target :: delp_eul_BC, pe_lag_BC, pe_eul_BC
   integer, intent(IN) :: npx, npy, npz, npz_coarse
   real, intent(IN), dimension(npz+1) :: ak_dst, bk_dst
   real, intent(IN) :: ptop_src

   integer :: i,j,k, istart, iend

   integer :: is,  ie,  js,  je
   integer :: isd, ied, jsd, jed

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed
   
   if (is == 1) then
      call setup_eul_delp_BC_k(delp_lag_BC%west_t1, delp_eul_BC%west_t1, pe_lag_BC%west_t1, pe_eul_BC%west_t1, &
           ptop_src, ak_dst, bk_dst, isd, 0, isd, 0, jsd, jed, npz, npz_coarse)
   end if

   if (ie == npx-1) then
      call setup_eul_delp_BC_k(delp_lag_BC%east_t1, delp_eul_BC%east_t1, pe_lag_BC%east_t1, pe_eul_BC%east_t1, &
           ptop_src, ak_dst, bk_dst, npx, ied, npx, ied, jsd, jed, npz, npz_coarse)
   end if

   if (is == 1) then
      istart = is
   else
      istart = isd
   end if
   if (ie == npx-1) then
      iend = ie
   else
      iend = ied
   end if

   if (js == 1) then
      call setup_eul_delp_BC_k(delp_lag_BC%south_t1, delp_eul_BC%south_t1, pe_lag_BC%south_t1, pe_eul_BC%south_t1, &
           ptop_src, ak_dst, bk_dst, isd, ied, istart, iend, jsd, 0, npz, npz_coarse)
   end if

   if (je == npy-1) then
      call setup_eul_delp_BC_k(delp_lag_BC%north_t1, delp_eul_BC%north_t1, pe_lag_BC%north_t1, pe_eul_BC%north_t1, &
           ptop_src, ak_dst, bk_dst, isd, ied, istart, iend, npy, jed, npz, npz_coarse)
   end if
   
 end subroutine setup_eul_delp_BC

 subroutine setup_eul_delp_BC_k(delplagBC, delpeulBC, pelagBC, peeulBC, ptop_src, ak_dst, bk_dst, isd_BC, ied_BC, istart, iend, jstart, jend, npz, npz_coarse)

   integer, intent(IN) :: isd_BC, ied_BC, istart, iend, jstart, jend, npz, npz_coarse
   real, intent(INOUT) :: delplagBC(isd_BC:ied_BC,jstart:jend,npz_coarse), pelagBC(isd_BC:ied_BC,jstart:jend,npz_coarse+1)
   real, intent(INOUT) :: delpeulBC(isd_BC:ied_BC,jstart:jend,npz), peeulBC(isd_BC:ied_BC,jstart:jend,npz+1)
   real, intent(IN) :: ptop_src, ak_dst(npz+1), bk_dst(npz+1)

   integer :: i,j,k

   character(len=120) :: errstring


!$OMP parallel do default(none) shared(istart,iend,jstart,jend,pelagBC,ptop_src)
   do j=jstart,jend
   do i=istart,iend
      pelagBC(i,j,1) = ptop_src 
   enddo
   enddo
!$OMP parallel do default(none) shared(istart,iend,jstart,jend,npz_coarse,pelagBC,delplagBC) 
   do j=jstart,jend
   do k=1,npz_coarse
   do i=istart,iend
      pelagBC(i,j,k+1) = pelagBC(i,j,k) + delplagBC(i,j,k)
   end do
   end do
   end do
!$OMP parallel do default(none) shared(istart,iend,jstart,jend,npz,npz_coarse,peeulBC,pelagBC,ak_dst,bk_dst)
   do k=1,npz+1
   do j=jstart,jend
   do i=istart,iend
      peeulBC(i,j,k) = ak_dst(k) + pelagBC(i,j,npz_coarse+1)*bk_dst(k)
   enddo
   enddo
   enddo
!$OMP parallel do default(none) shared(istart,iend,jstart,jend,npz,peeulBC,delpeulBC)
   do k=1,npz
   do j=jstart,jend
   do i=istart,iend
      delpeulBC(i,j,k) = peeulBC(i,j,k+1) - peeulBC(i,j,k)
   enddo
   enddo
   enddo

!!$!!! DEBUG CODE
!!$   !If more than a few percent difference then log the error
!!$   do k=1,npz
!!$   do j=jstart,jend
!!$   do i=istart,iend
!!$      if (delpeulBC(i,j,k) <= 0.) then
!!$         write(errstring,'(3I5, 3(2x, G))'), i, j, k, pelagBC(i,j,k), peeulBC(i,j,k)
!!$         call mpp_error(WARNING, ' Invalid pressure BC at '//errstring)
!!$      else if (abs( peeulBC(i,j,k) - pelagBC(i,j,k)) > 100.0 ) then
!!$         write(errstring,'(3I5, 3(2x, G))'), i, j, k, pelagBC(i,j,k), peeulBC(i,j,k)
!!$         call mpp_error(WARNING, ' Remap BC: pressure deviation at '//errstring)
!!$      endif
!!$   enddo
!!$   enddo
!!$   enddo
!!$!!! END DEBUG CODE

 end subroutine setup_eul_delp_BC_k

 subroutine copy_ps_BC(ps, pe_BC, npx, npy, npz, istag, jstag, bd)

   integer, intent(IN) :: npx, npy, npz, istag, jstag
   type(fv_grid_bounds_type), intent(IN) :: bd
   real, intent(IN) :: ps(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag)
   type(fv_nest_BC_type_3d), intent(INOUT) :: pe_BC

   integer :: i,j,k, istart, iend

   integer :: is,  ie,  js,  je
   integer :: isd, ied, jsd, jed

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed
   
   if (is == 1) then
!$OMP parallel do default(none) shared(isd,jsd,jed,jstag,npz,pe_BC,ps)
      do j=jsd,jed+jstag
      do i=isd,0
         pe_BC%west_t1(i,j,npz+1) = ps(i,j)
      enddo
      enddo
   end if

   if (ie == npx-1) then
!$OMP parallel do default(none) shared(npx,ied,istag,jsd,jed,jstag,npz,pe_BC,ps)
      do j=jsd,jed+jstag
      do i=npx+istag,ied+istag
         pe_BC%east_t1(i,j,npz+1) = ps(i,j)
      enddo
      enddo
   end if

   if (is == 1) then
      istart = is
   else
      istart = isd
   end if
   if (ie == npx-1) then
      iend = ie
   else
      iend = ied
   end if

   if (js == 1) then
!$OMP parallel do default(none) shared(isd,ied,istag,jsd,npz,pe_BC,ps)
      do j=jsd,0
      do i=isd,ied+istag
         pe_BC%south_t1(i,j,npz+1) = ps(i,j)
      enddo
      enddo
   end if

   if (je == npy-1) then
!$OMP parallel do default(none) shared(isd,ied,istag,npy,jed,jstag,npz,pe_BC,ps)
      do j=npy+jstag,jed+jstag
      do i=isd,ied+istag
         pe_BC%north_t1(i,j,npz+1) = ps(i,j)
      enddo
      enddo
   end if

 end subroutine copy_ps_BC

!In this routine, the pe_*_BC arrays should already have PS filled in on the npz+1 level
 subroutine setup_eul_pe_BC(pe_src_BC, pe_eul_BC, ak_dst, bk_dst, npx, npy, npz, npz_src, istag, jstag, bd, make_src_in, ak_src, bk_src)

   type(fv_grid_bounds_type), intent(IN) :: bd
   type(fv_nest_BC_type_3d), intent(INOUT), target :: pe_src_BC, pe_eul_BC
   integer, intent(IN) :: npx, npy, npz, npz_src, istag, jstag
   real, intent(IN), dimension(npz+1) :: ak_dst, bk_dst
   logical, intent(IN), OPTIONAL :: make_src_in
   real, intent(IN), OPTIONAL :: ak_src(npz_src), bk_src(npz_src)

   logical :: make_src

   integer :: i,j,k, istart, iend

   integer :: is,  ie,  js,  je
   integer :: isd, ied, jsd, jed

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed
   
   make_src = .false.
   if (present(make_src_in)) make_src = make_src_in

   if (is == 1) then
      call setup_eul_pe_BC_k(pe_src_BC%west_t1, pe_eul_BC%west_t1, ak_dst, bk_dst, isd, 0, isd, 0, jsd, jed+jstag, npz, npz_src, &
           make_src, ak_src, bk_src)
   end if

   if (ie == npx-1) then
      call setup_eul_pe_BC_k(pe_src_BC%east_t1, pe_eul_BC%east_t1, ak_dst, bk_dst, npx+istag, ied+istag, npx+istag, ied+istag, jsd, jed+jstag, npz, npz_src, &
           make_src, ak_src, bk_src)
   end if

   if (is == 1) then
      istart = is
   else
      istart = isd
   end if
   if (ie == npx-1) then
      iend = ie
   else
      iend = ied
   end if

   if (js == 1) then
      call setup_eul_pe_BC_k(pe_src_BC%south_t1, pe_eul_BC%south_t1, ak_dst, bk_dst, isd, ied+istag, istart, iend+istag, jsd, 0, npz, npz_src, &
           make_src, ak_src, bk_src)
   end if

   if (je == npy-1) then
      call setup_eul_pe_BC_k(pe_src_BC%north_t1, pe_eul_BC%north_t1, ak_dst, bk_dst, isd, ied+istag, istart, iend+istag, npy+jstag, jed+jstag, npz, npz_src, &
           make_src, ak_src, bk_src)
   end if
   
 end subroutine setup_eul_pe_BC

 subroutine setup_eul_pe_BC_k(pesrcBC, peeulBC, ak_dst, bk_dst, isd_BC, ied_BC, istart, iend, jstart, jend, npz, npz_src, make_src, ak_src, bk_src)

   integer, intent(IN) :: isd_BC, ied_BC, istart, iend, jstart, jend, npz, npz_src
   real, intent(INOUT) :: pesrcBC(isd_BC:ied_BC,jstart:jend,npz_src+1)
   real, intent(INOUT) :: peeulBC(isd_BC:ied_BC,jstart:jend,npz+1)
   real, intent(IN) :: ak_dst(npz+1), bk_dst(npz+1)
   logical, intent(IN) :: make_src
   real, intent(IN) :: ak_src(npz_src+1), bk_src(npz_src+1)

   integer :: i,j,k

   character(len=120) :: errstring

!$OMP parallel do default(none) shared(istart,iend,jstart,jend,npz,npz_src,peeulBC,ak_dst,pesrcBC,bk_dst)
   do k=1,npz+1
   do j=jstart,jend
   do i=istart,iend
      peeulBC(i,j,k) = ak_dst(k) + pesrcBC(i,j,npz_src+1)*bk_dst(k)
   enddo
   enddo
   enddo
   
   if (make_src) then
!$OMP parallel do default(none) shared(istart,iend,jstart,jend,npz_src,pesrcBC,ak_src,bk_src)
   do k=1,npz_src+1
   do j=jstart,jend
   do i=istart,iend
      pesrcBC(i,j,k) = ak_src(k) + pesrcBC(i,j,npz_src+1)*bk_src(k)
   enddo
   enddo
   enddo
   endif


 end subroutine setup_eul_pe_BC_k

 subroutine remap_BC(pe_lag_BC, pe_eul_BC, var_lag_BC, var_eul_BC, npx, npy, npz, npz_coarse, bd, istag, jstag, iv, kord, varname, do_log_pe)

   type(fv_grid_bounds_type), intent(IN) :: bd
   type(fv_nest_BC_type_3d), intent(INOUT), target :: pe_lag_BC, var_lag_BC
   type(fv_nest_BC_type_3d), intent(INOUT), target :: pe_eul_BC, var_eul_BC
   integer, intent(IN) :: npx, npy, npz, npz_coarse, istag, jstag, iv, kord
   character(len=*), intent(IN) :: varname
   logical, intent(IN), OPTIONAL :: do_log_pe

   logical :: log_pe = .false.

   integer :: i,j,k, istart, iend

   integer :: is,  ie,  js,  je
   integer :: isd, ied, jsd, jed

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed
   
   if (present(do_log_pe)) log_pe = do_log_pe
   
   if (is == 1) then
      call remap_BC_k(pe_lag_BC%west_t1, pe_eul_BC%west_t1, var_lag_BC%west_t1, var_eul_BC%west_t1, isd, 0, isd, 0, jsd, jed+jstag, npz, npz_coarse, iv, kord, log_pe)
   end if

   if (ie == npx-1) then
      call remap_BC_k(pe_lag_BC%east_t1, pe_eul_BC%east_t1, var_lag_BC%east_t1, var_eul_BC%east_t1, npx+istag, ied+istag, npx+istag, ied+istag, jsd, jed+jstag, npz, npz_coarse, iv, kord, log_pe)
   end if

   if (is == 1) then
      istart = is
   else
      istart = isd
   end if
   if (ie == npx-1) then
      iend = ie
   else
      iend = ied
   end if

   if (js == 1) then
      call remap_BC_k(pe_lag_BC%south_t1, pe_eul_BC%south_t1, var_lag_BC%south_t1, var_eul_BC%south_t1, isd, ied+istag, istart, iend+istag, jsd, 0, npz, npz_coarse, iv, kord, log_pe)
   end if

   if (je == npy-1) then
      call remap_BC_k(pe_lag_BC%north_t1, pe_eul_BC%north_t1, var_lag_BC%north_t1, var_eul_BC%north_t1, isd, ied+istag, istart, iend+istag, npy+jstag, jed+jstag, npz, npz_coarse, iv, kord, log_pe)
   end if
   
 end subroutine remap_BC

 subroutine remap_BC_direct(pe_lag_BC, pe_eul_BC, var_lag_BC, var, npx, npy, npz, npz_coarse, bd, istag, jstag, iv, kord, do_log_pe)

   type(fv_grid_bounds_type), intent(IN) :: bd
   integer, intent(IN) :: npx, npy, npz, npz_coarse, istag, jstag, iv, kord
   type(fv_nest_BC_type_3d), intent(INOUT), target :: pe_lag_BC, var_lag_BC
   type(fv_nest_BC_type_3d), intent(INOUT), target :: pe_eul_BC
   real, intent(INOUT) ::  var(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,npz)
   logical, intent(IN), OPTIONAL :: do_log_pe

   logical :: log_pe = .false.

   integer :: i,j,k, istart, iend

   integer :: is,  ie,  js,  je
   integer :: isd, ied, jsd, jed

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed
   
   if (present(do_log_pe)) log_pe = do_log_pe
   
   if (is == 1) then
      !I was unable how to do pass-by-memory referencing on parts of the 3D var array,
      ! so instead I am doing an inefficient copy and copy-back. --- lmh 14jun17
      call remap_BC_k(pe_lag_BC%west_t1, pe_eul_BC%west_t1, var_lag_BC%west_t1, var(isd:0,jsd:jed+jstag,:), isd, 0, isd, 0, jsd, jed+jstag, npz, npz_coarse, iv, kord, log_pe)
   end if

   if (ie == npx-1) then
      call remap_BC_k(pe_lag_BC%east_t1, pe_eul_BC%east_t1, var_lag_BC%east_t1, var(npx+istag:ied+istag,jsd:jed+jstag,:), npx+istag, ied+istag, npx+istag, ied+istag, jsd, jed+jstag, npz, npz_coarse, iv, kord, log_pe)
   end if

   if (is == 1) then
      istart = is
   else
      istart = isd
   end if
   if (ie == npx-1) then
      iend = ie
   else
      iend = ied
   end if

   if (js == 1) then
      call remap_BC_k(pe_lag_BC%south_t1, pe_eul_BC%south_t1, var_lag_BC%south_t1, var(isd:ied+istag,jsd:0,:), isd, ied+istag, istart, iend+istag, jsd, 0, npz, npz_coarse, iv, kord, log_pe)
   end if

   if (je == npy-1) then
      call remap_BC_k(pe_lag_BC%north_t1, pe_eul_BC%north_t1, var_lag_BC%north_t1, var(isd:ied+istag,npy+jstag:jed+jstag,:), isd, ied+istag, istart, iend+istag, npy+jstag, jed+jstag, npz, npz_coarse, iv, kord, log_pe)
   end if
   
 end subroutine remap_BC_direct

 subroutine remap_BC_k(pe_lagBC, pe_eulBC, var_lagBC, var_eulBC, isd_BC, ied_BC, istart, iend, jstart, jend, npz, npz_coarse, iv, kord, log_pe)

   integer, intent(IN) :: isd_BC, ied_BC, istart, iend, jstart, jend, npz, npz_coarse, iv, kord
   logical, intent(IN) :: log_pe
   real, intent(INOUT) :: pe_lagBC(isd_BC:ied_BC,jstart:jend,npz_coarse+1), var_lagBC(isd_BC:ied_BC,jstart:jend,npz_coarse)
   real, intent(INOUT) :: pe_eulBC(isd_BC:ied_BC,jstart:jend,npz+1), var_eulBC(isd_BC:ied_BC,jstart:jend,npz)

   integer :: i, j, k
   real peln_lag(istart:iend,npz_coarse+1)
   real peln_eul(istart:iend,npz+1)
   character(120) :: errstring
   
   if (log_pe) then

!$OMP parallel do default(none) shared(istart,iend,jstart,jend,npz,npz_coarse,pe_lagBC,pe_eulBC,var_lagBC,var_eulBC,iv,kord) &
!$OMP                           private(peln_lag,peln_eul)
      do j=jstart,jend

         do k=1,npz_coarse+1
         do i=istart,iend
!!$!!! DEBUG CODE
!!$            if (pe_lagBC(i,j,k) <= 0.) then
!!$               write(errstring,'(3I5, 2x, G)'), i, j, k, pe_lagBC(i,j,k)
!!$               call mpp_error(WARNING, ' Remap BC: invalid pressure at at '//errstring)               
!!$            endif
!!$!!! END DEBUG CODE
            peln_lag(i,k) = log(pe_lagBC(i,j,k))
         enddo
         enddo

         do k=1,npz+1
         do i=istart,iend
!!$!!! DEBUG CODE
!!$            if (pe_lagBC(i,j,k) <= 0.) then
!!$               write(errstring,'(3I5, 2x, G)'), i, j, k, pe_lagBC(i,j,k)
!!$               call mpp_error(WARNING, ' Remap BC: invalid pressure at at '//errstring)               
!!$            endif
!!$!!! END DEBUG CODE
            peln_eul(i,k) = log(pe_eulBC(i,j,k))
         enddo
         enddo

         call mappm(npz_coarse, peln_lag, var_lagBC(istart:iend,j:j,:), &
                    npz, peln_eul, var_eulBC(istart:iend,j:j,:), &
                    istart, iend, iv, kord, pe_eulBC(istart,j,1))

      enddo

   else

!$OMP parallel do default(none) shared(istart,iend,jstart,jend,npz,npz_coarse,pe_lagBC,pe_eulBC,var_lagBC,var_eulBC,iv,kord)
      do j=jstart,jend

         call mappm(npz_coarse, pe_lagBC(istart:iend,j:j,:), var_lagBC(istart:iend,j:j,:), &
                    npz, pe_eulBC(istart:iend,j:j,:), var_eulBC(istart:iend,j:j,:), &
                    istart, iend, iv, kord, pe_eulBC(istart,j,1))
         !!! NEED A FILLQ/FILLZ CALL HERE??

      enddo
   endif

 end subroutine remap_BC_k

 subroutine remap_delz_BC(pe_lag_BC, pe_eul_BC, delp_lag_BC, delz_lag_BC, delp_eul_BC, delz_eul_BC, npx, npy, npz, npz_coarse, bd, istag, jstag, iv, kord)

   type(fv_grid_bounds_type), intent(IN) :: bd
   type(fv_nest_BC_type_3d), intent(INOUT), target :: pe_lag_BC, delp_lag_BC, delz_lag_BC
   type(fv_nest_BC_type_3d), intent(INOUT), target :: pe_eul_BC, delp_eul_BC, delz_eul_BC
   integer, intent(IN) :: npx, npy, npz, npz_coarse, istag, jstag, iv, kord

   integer :: i,j,k, istart, iend

   integer :: is,  ie,  js,  je
   integer :: isd, ied, jsd, jed

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed
   
   if (is == 1) then
      call compute_specific_volume_BC_k(delp_lag_BC%west_t1, delz_lag_BC%west_t1, isd, 0, isd, 0, jsd, jed, npz_coarse)
      call remap_BC_k(pe_lag_BC%west_t1, pe_eul_BC%west_t1, delz_lag_BC%west_t1, delz_eul_BC%west_t1, isd, 0, isd, 0, jsd, jed+jstag, &
           npz, npz_coarse, iv, kord, log_pe=.false.)
      call compute_delz_BC_k(delp_eul_BC%west_t1, delz_eul_BC%west_t1, isd, 0, isd, 0, jsd, jed, npz)
   end if

   if (ie == npx-1) then
      call compute_specific_volume_BC_k(delp_lag_BC%east_t1, delz_lag_BC%east_t1, npx+istag, ied+istag, npx+istag, ied+istag, jsd, jed+jstag, npz_coarse)
      call remap_BC_k(pe_lag_BC%east_t1, pe_eul_BC%east_t1, delz_lag_BC%east_t1, delz_eul_BC%east_t1, npx+istag, ied+istag, npx+istag, ied+istag, jsd, jed+jstag, &
           npz, npz_coarse, iv, kord, log_pe=.false.)
      call compute_delz_BC_k(delp_eul_BC%east_t1, delz_eul_BC%east_t1, npx+istag, ied+istag, npx+istag, ied+istag, jsd, jed+jstag, npz)
   end if

   if (is == 1) then
      istart = is
   else
      istart = isd
   end if
   if (ie == npx-1) then
      iend = ie
   else
      iend = ied
   end if

   if (js == 1) then
      call compute_specific_volume_BC_k(delp_lag_BC%south_t1, delz_lag_BC%south_t1, isd, ied+istag, istart, iend+istag, jsd, 0, npz_coarse)
      call remap_BC_k(pe_lag_BC%south_t1, pe_eul_BC%south_t1, delz_lag_BC%south_t1, delz_eul_BC%south_t1, isd, ied+istag, istart, iend+istag, jsd, 0, npz, npz_coarse, &
           iv, kord, log_pe=.false.)
      call compute_delz_BC_k(delp_eul_BC%south_t1, delz_eul_BC%south_t1, isd, ied+istag, istart, iend+istag, jsd, 0, npz)
   end if

   if (je == npy-1) then
      call compute_specific_volume_BC_k(delp_lag_BC%north_t1, delz_lag_BC%north_t1, isd, ied+istag, istart, iend+istag, npy+jstag, jed+jstag, npz_coarse)
      call remap_BC_k(pe_lag_BC%north_t1, pe_eul_BC%north_t1, delz_lag_BC%north_t1, delz_eul_BC%north_t1, &
           isd, ied+istag, istart, iend+istag, npy+jstag, jed+jstag, npz, npz_coarse, iv, kord, log_pe=.false.)
      call compute_delz_BC_k(delp_eul_BC%north_t1, delz_eul_BC%north_t1, isd, ied+istag, istart, iend+istag, npy+jstag, jed+jstag, npz)
   end if
   
 end subroutine remap_delz_BC

 subroutine compute_specific_volume_BC_k(delpBC, delzBC, isd_BC, ied_BC, istart, iend, jstart, jend, npz)

   integer, intent(IN) :: isd_BC, ied_BC, istart, iend, jstart, jend, npz
   real, intent(IN)    :: delpBC(isd_BC:ied_BC,jstart:jend,npz)
   real, intent(INOUT) :: delzBC(isd_BC:ied_BC,jstart:jend,npz)

   character(len=120) :: errstring
   integer :: i,j,k
   
!$OMP parallel do default(none) shared(istart,iend,jstart,jend,npz,delzBC,delpBC)
   do k=1,npz
   do j=jstart,jend
   do i=istart,iend
      delzBC(i,j,k) = -delzBC(i,j,k)/delpBC(i,j,k)
!!$!!! DEBUG CODE 
!!$      if (delzBC(i,j,k) <= 0. ) then
!!$         write(errstring,'(3I5, 2(2x, G))'), i, j, k, delzBC(i,j,k), delpBC(i,j,k)
!!$         call mpp_error(WARNING, ' Remap BC (sfc volume): invalid delz at '//errstring)               
!!$      endif
!!$!!! END DEBUG CODE
   end do
   end do
   end do

 end subroutine compute_specific_volume_BC_k

 subroutine compute_delz_BC_k(delpBC, delzBC, isd_BC, ied_BC, istart, iend, jstart, jend, npz)

   integer, intent(IN) :: isd_BC, ied_BC, istart, iend, jstart, jend, npz
   real, intent(IN)    :: delpBC(isd_BC:ied_BC,jstart:jend,npz)
   real, intent(INOUT) :: delzBC(isd_BC:ied_BC,jstart:jend,npz)
   
   character(len=120) :: errstring
   integer :: i,j,k
   
!$OMP parallel do default(none) shared(istart,iend,jstart,jend,npz,delzBC,delpBC)
   do k=1,npz
   do j=jstart,jend
   do i=istart,iend
      delzBC(i,j,k) = -delzBC(i,j,k)*delpBC(i,j,k) 
!!$!!! DEBUG CODE
!!$      if (delzBC(i,j,k) >=0. ) then
!!$         write(errstring,'(3I5, 2(2x, G))'), i, j, k, delzBC(i,j,k), delpBC(i,j,k)
!!$         call mpp_error(WARNING, ' Remap BC (compute delz): invalid delz at '//errstring)               
!!$      endif
!!$!!! END DEBUG CODE
   end do
   end do
   end do

 end subroutine compute_delz_BC_k


 subroutine setup_pt_NH_BC(pt_BC, delp_BC, delz_BC, sphum_BC, q_BC, nq, &
#ifdef USE_COND
      q_con_BC, &
#ifdef MOIST_CAPPA
      cappa_BC, &
#endif
#endif
      npx, npy, npz, zvir, bd)

   type(fv_grid_bounds_type), intent(IN) :: bd
   type(fv_nest_BC_type_3d), intent(IN), target    :: delp_BC, delz_BC, sphum_BC
   type(fv_nest_BC_type_3d), intent(INOUT), target :: pt_BC
   integer, intent(IN) :: nq
   type(fv_nest_BC_type_3d), intent(IN), target :: q_BC(nq)
#ifdef USE_COND
   type(fv_nest_BC_type_3d), intent(INOUT), target :: q_con_BC
#ifdef MOIST_CAPPA
   type(fv_nest_BC_type_3d), intent(INOUT), target :: cappa_BC
#endif
#endif
   integer, intent(IN) :: npx, npy, npz
   real, intent(IN) :: zvir

    real, parameter:: c_liq = 4185.5      !< heat capacity of water at 0C
    real, parameter:: c_ice = 1972.       !< heat capacity of ice at 0C: c=c_ice+7.3*(T-Tice) 
    real, parameter:: cv_vap = cp_vapor - rvgas  !< 1384.5

   real, dimension(:,:,:), pointer :: liq_watBC_west, ice_watBC_west, rainwatBC_west, snowwatBC_west, graupelBC_west
   real, dimension(:,:,:), pointer :: liq_watBC_east, ice_watBC_east, rainwatBC_east, snowwatBC_east, graupelBC_east
   real, dimension(:,:,:), pointer :: liq_watBC_north, ice_watBC_north, rainwatBC_north, snowwatBC_north, graupelBC_north
   real, dimension(:,:,:), pointer :: liq_watBC_south, ice_watBC_south, rainwatBC_south, snowwatBC_south, graupelBC_south

   real :: dp1, q_liq, q_sol, q_con = 0., cvm, pkz, rdg, cv_air

   integer :: i,j,k, istart, iend
   integer :: liq_wat, ice_wat, rainwat, snowwat, graupel
   real, parameter:: tice = 273.16 !< For GFS Partitioning
   real, parameter:: t_i0 = 15.

   integer :: is,  ie,  js,  je
   integer :: isd, ied, jsd, jed

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed
   
   rdg = -rdgas / grav
   cv_air =  cp_air - rdgas

   liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
   ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
   rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat')
   snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat')
   graupel = get_tracer_index (MODEL_ATMOS, 'graupel')   

   if (is == 1) then
      if (.not. allocated(dum_West)) then
         allocate(dum_West(isd:0,jsd:jed,npz))
!$OMP parallel do default(none) shared(npz,isd,jsd,jed,dum_West)
         do k=1,npz
         do j=jsd,jed
         do i=isd,0
            dum_West(i,j,k) = 0.
         enddo
         enddo
         enddo
      endif
   endif
   if (js == 1) then
      if (.not. allocated(dum_South)) then
         allocate(dum_South(isd:ied,jsd:0,npz))
!$OMP parallel do default(none) shared(npz,isd,ied,jsd,dum_South)
         do k=1,npz
         do j=jsd,0
         do i=isd,ied
            dum_South(i,j,k) = 0.
         enddo
         enddo
         enddo
      endif
   endif
   if (ie == npx-1) then
      if (.not. allocated(dum_East)) then
         allocate(dum_East(npx:ied,jsd:jed,npz))
!$OMP parallel do default(none) shared(npx,npz,ied,jsd,jed,dum_East)
         do k=1,npz
         do j=jsd,jed
         do i=npx,ied
            dum_East(i,j,k) = 0.
         enddo
         enddo
         enddo
      endif
   endif
   if (je == npy-1) then
      if (.not. allocated(dum_North)) then
         allocate(dum_North(isd:ied,npy:jed,npz))
!$OMP parallel do default(none) shared(npy,npz,isd,ied,jed,dum_North)
         do k=1,npz
         do j=npy,jed
         do i=isd,ied
            dum_North(i,j,k) = 0.
         enddo
         enddo
         enddo
      endif
   endif

   if (liq_wat > 0) then
      liq_watBC_west  => q_BC(liq_wat)%west_t1
      liq_watBC_east  => q_BC(liq_wat)%east_t1
      liq_watBC_north => q_BC(liq_wat)%north_t1
      liq_watBC_south => q_BC(liq_wat)%south_t1
   else
      liq_watBC_west  => dum_west
      liq_watBC_east  => dum_east
      liq_watBC_north => dum_north
      liq_watBC_south => dum_south
   endif
   if (ice_wat > 0) then
      ice_watBC_west  => q_BC(ice_wat)%west_t1
      ice_watBC_east  => q_BC(ice_wat)%east_t1
      ice_watBC_north => q_BC(ice_wat)%north_t1
      ice_watBC_south => q_BC(ice_wat)%south_t1
   else
      ice_watBC_west  => dum_west
      ice_watBC_east  => dum_east
      ice_watBC_north => dum_north
      ice_watBC_south => dum_south
   endif
   if (rainwat > 0) then
      rainwatBC_west  => q_BC(rainwat)%west_t1
      rainwatBC_east  => q_BC(rainwat)%east_t1
      rainwatBC_north => q_BC(rainwat)%north_t1
      rainwatBC_south => q_BC(rainwat)%south_t1
   else
      rainwatBC_west  => dum_west
      rainwatBC_east  => dum_east
      rainwatBC_north => dum_north
      rainwatBC_south => dum_south
   endif
   if (snowwat > 0) then
      snowwatBC_west  => q_BC(snowwat)%west_t1
      snowwatBC_east  => q_BC(snowwat)%east_t1
      snowwatBC_north => q_BC(snowwat)%north_t1
      snowwatBC_south => q_BC(snowwat)%south_t1
   else
      snowwatBC_west  => dum_west
      snowwatBC_east  => dum_east
      snowwatBC_north => dum_north
      snowwatBC_south => dum_south
   endif
   if (graupel > 0) then
      graupelBC_west  => q_BC(graupel)%west_t1
      graupelBC_east  => q_BC(graupel)%east_t1
      graupelBC_north => q_BC(graupel)%north_t1
      graupelBC_south => q_BC(graupel)%south_t1
   else
      graupelBC_west  => dum_west
      graupelBC_east  => dum_east
      graupelBC_north => dum_north
      graupelBC_south => dum_south
   endif

   if (is == 1) then

      call setup_pt_NH_BC_k(pt_BC%west_t1, sphum_BC%west_t1, delp_BC%west_t1, delz_BC%west_t1, &
           liq_watBC_west, rainwatBC_west, ice_watBC_west, snowwatBC_west, graupelBC_west, &
#ifdef USE_COND
           q_con_BC%west_t1, &
#ifdef MOIST_CAPPA
           cappa_BC%west_t1, &
#endif
#endif
           zvir, isd, 0, isd, 0, jsd, jed, npz)
   end if


   if (js == 1) then
      if (is == 1) then
         istart = is
      else
         istart = isd
      end if
      if (ie == npx-1) then
         iend = ie
      else
         iend = ied
      end if

      call setup_pt_NH_BC_k(pt_BC%south_t1, sphum_BC%south_t1, delp_BC%south_t1, delz_BC%south_t1, &
           liq_watBC_south, rainwatBC_south, ice_watBC_south, snowwatBC_south, graupelBC_south, &
#ifdef USE_COND
           q_con_BC%south_t1, &
#ifdef MOIST_CAPPA
           cappa_BC%south_t1, &
#endif
#endif
           zvir, isd, ied, istart, iend, jsd, 0, npz)
   end if


   if (ie == npx-1) then

      call setup_pt_NH_BC_k(pt_BC%east_t1, sphum_BC%east_t1, delp_BC%east_t1, delz_BC%east_t1, &
           liq_watBC_east, rainwatBC_east, ice_watBC_east, snowwatBC_east, graupelBC_east, &
#ifdef USE_COND
           q_con_BC%east_t1, &
#ifdef MOIST_CAPPA
           cappa_BC%east_t1, &
#endif
#endif
           zvir, npx, ied, npx, ied, jsd, jed, npz)
   end if

   if (je == npy-1) then
      if (is == 1) then
         istart = is
      else
         istart = isd
      end if
      if (ie == npx-1) then
         iend = ie
      else
         iend = ied
      end if

      call setup_pt_NH_BC_k(pt_BC%north_t1, sphum_BC%north_t1, delp_BC%north_t1, delz_BC%north_t1, &
           liq_watBC_north, rainwatBC_north, ice_watBC_north, snowwatBC_north, graupelBC_north, &
#ifdef USE_COND
           q_con_BC%north_t1, &
#ifdef MOIST_CAPPA
           cappa_BC%north_t1, &
#endif
#endif
           zvir, isd, ied, istart, iend, npy, jed, npz)
   end if

 end subroutine setup_pt_NH_BC


 subroutine setup_pt_NH_BC_k(ptBC,sphumBC,delpBC,delzBC, &
                             liq_watBC,rainwatBC,ice_watBC,snowwatBC,graupelBC, &
#ifdef USE_COND
                             q_conBC, &
#ifdef MOIST_CAPPA
                             cappaBC, &
#endif
#endif
                             zvir, isd_BC, ied_BC, istart, iend, jstart, jend, npz)

   integer, intent(IN) :: isd_BC, ied_BC, istart, iend, jstart, jend, npz
   real, intent(OUT), dimension(isd_BC:ied_BC,jstart:jend,npz) :: ptBC
   real, intent(IN),  dimension(isd_BC:ied_BC,jstart:jend,npz) :: sphumBC, delpBC, delzBC
   real, intent(IN),  dimension(isd_BC:ied_BC,jstart:jend,npz) :: liq_watBC,rainwatBC,ice_watBC,snowwatBC,graupelBC
#ifdef USE_COND
   real, intent(OUT), dimension(isd_BC:ied_BC,jstart:jend,npz) ::   q_conBC
#ifdef MOIST_CAPPA
   real, intent(OUT), dimension(isd_BC:ied_BC,jstart:jend,npz) ::   cappaBC
#endif
#endif
   real, intent(IN) :: zvir

   integer :: i,j,k
   real :: dp1, q_con, q_sol, q_liq, cvm, pkz, rdg, cv_air

   real, parameter:: c_liq = 4185.5      ! heat capacity of water at 0C
   real, parameter:: c_ice = 1972.       ! heat capacity of ice at 0C: c=c_ice+7.3*(T-Tice) 
   real, parameter:: cv_vap = cp_vapor - rvgas  ! 1384.5
   real, parameter:: tice = 273.16 ! For GFS Partitioning
   real, parameter:: t_i0 = 15.

   rdg = -rdgas / grav
   cv_air =  cp_air - rdgas

!!$!!! DEBUG CODE
!!$   write(*, '(A, 7I5)') 'setup_pt_NH_BC_k', mpp_pe(), isd, ied, istart, iend, lbound(ptBC,1), ubound(ptBC,1)
!!$!!! END DEBUG CODE

!$OMP parallel do default(none) shared(istart,iend,jstart,jend,npz,zvir,ptBC,sphumBC,delpBC,delzBC,liq_watBC,rainwatBC,ice_watBC,snowwatBC,graupelBC, &
#ifdef USE_COND
!$OMP                                  q_conBC, &
#ifdef MOIST_CAPPA
!$OMP                                  cappaBC, &
#endif
#endif
!$OMP                                  rdg, cv_air) &
!$OMP                          private(dp1,q_liq,q_sol,q_con,cvm,pkz) 
   do k=1,npz
   do j=jstart,jend
   do i=istart,iend
         dp1 = zvir*sphumBC(i,j,k)
#ifdef USE_COND
         q_liq = liq_watBC(i,j,k) + rainwatBC(i,j,k)
         q_sol = ice_watBC(i,j,k) + snowwatBC(i,j,k) + graupelBC(i,j,k)
         q_con = q_liq + q_sol
         q_conBC(i,j,k) = q_con
#ifdef MOIST_CAPPA
         cvm = (1.-(sphumBC(i,j,k)+q_con))*cv_air+sphumBC(i,j,k)*cv_vap+q_liq*c_liq+q_sol*c_ice
         cappaBC(i,j,k) = rdgas/(rdgas + cvm/(1.+dp1))
         pkz = exp( cappaBC(i,j,k)*log(rdg*delpBC(i,j,k)*ptBC(i,j,k) * &
              (1.+dp1)*(1.-q_con)/delzBC(i,j,k)))         
#else
         pkz = exp( kappa*log(rdg*delpBC(i,j,k)*ptBC(i,j,k) * &
              (1.+dp1)*(1.-q_con)/delzBC(i,j,k)))
#endif
         ptBC(i,j,k) = ptBC(i,j,k)*(1.+dp1)*(1.-q_con)/pkz
#else
         pkz = exp( kappa*log(rdg*delpBC(i,j,k)*ptBC(i,j,k) * &
              (1.+dp1)/delzBC(i,j,k)))
         ptBC(i,j,k) = ptBC(i,j,k)*(1.+dp1)/pkz
#endif
   end do
   end do
   end do

 end subroutine setup_pt_NH_BC_k

 subroutine set_NH_BCs_t0(neststruct)

   type(fv_nest_type), intent(INOUT) :: neststruct

#ifndef SW_DYNAMICS
   neststruct%delz_BC%east_t0  = neststruct%delz_BC%east_t1
   neststruct%delz_BC%west_t0  = neststruct%delz_BC%west_t1
   neststruct%delz_BC%north_t0 = neststruct%delz_BC%north_t1
   neststruct%delz_BC%south_t0 = neststruct%delz_BC%south_t1

   neststruct%w_BC%east_t0  = neststruct%w_BC%east_t1
   neststruct%w_BC%west_t0  = neststruct%w_BC%west_t1
   neststruct%w_BC%north_t0 = neststruct%w_BC%north_t1
   neststruct%w_BC%south_t0 = neststruct%w_BC%south_t1
#endif

 end subroutine set_NH_BCs_t0

 subroutine set_BCs_t0(ncnst, hydrostatic, neststruct)

   integer, intent(IN) :: ncnst
   logical, intent(IN) :: hydrostatic
   type(fv_nest_type), intent(INOUT) :: neststruct

   integer :: n

   neststruct%delp_BC%east_t0  = neststruct%delp_BC%east_t1
   neststruct%delp_BC%west_t0  = neststruct%delp_BC%west_t1
   neststruct%delp_BC%north_t0 = neststruct%delp_BC%north_t1
   neststruct%delp_BC%south_t0 = neststruct%delp_BC%south_t1
   do n=1,ncnst
      neststruct%q_BC(n)%east_t0  = neststruct%q_BC(n)%east_t1
      neststruct%q_BC(n)%west_t0  = neststruct%q_BC(n)%west_t1
      neststruct%q_BC(n)%north_t0 = neststruct%q_BC(n)%north_t1
      neststruct%q_BC(n)%south_t0 = neststruct%q_BC(n)%south_t1
   enddo
#ifndef SW_DYNAMICS
   neststruct%pt_BC%east_t0    = neststruct%pt_BC%east_t1
   neststruct%pt_BC%west_t0    = neststruct%pt_BC%west_t1
   neststruct%pt_BC%north_t0   = neststruct%pt_BC%north_t1
   neststruct%pt_BC%south_t0   = neststruct%pt_BC%south_t1
   neststruct%pt_BC%east_t0    = neststruct%pt_BC%east_t1
   neststruct%pt_BC%west_t0    = neststruct%pt_BC%west_t1
   neststruct%pt_BC%north_t0   = neststruct%pt_BC%north_t1
   neststruct%pt_BC%south_t0   = neststruct%pt_BC%south_t1

#ifdef USE_COND
   neststruct%q_con_BC%east_t0    = neststruct%q_con_BC%east_t1
   neststruct%q_con_BC%west_t0    = neststruct%q_con_BC%west_t1
   neststruct%q_con_BC%north_t0   = neststruct%q_con_BC%north_t1
   neststruct%q_con_BC%south_t0   = neststruct%q_con_BC%south_t1
#ifdef MOIST_CAPPA
   neststruct%cappa_BC%east_t0    = neststruct%cappa_BC%east_t1
   neststruct%cappa_BC%west_t0    = neststruct%cappa_BC%west_t1
   neststruct%cappa_BC%north_t0   = neststruct%cappa_BC%north_t1
   neststruct%cappa_BC%south_t0   = neststruct%cappa_BC%south_t1
#endif
#endif

   if (.not. hydrostatic) then
      call set_NH_BCs_t0(neststruct)
   endif
#endif
   neststruct%u_BC%east_t0  = neststruct%u_BC%east_t1
   neststruct%u_BC%west_t0  = neststruct%u_BC%west_t1
   neststruct%u_BC%north_t0 = neststruct%u_BC%north_t1
   neststruct%u_BC%south_t0 = neststruct%u_BC%south_t1
   neststruct%v_BC%east_t0  = neststruct%v_BC%east_t1
   neststruct%v_BC%west_t0  = neststruct%v_BC%west_t1
   neststruct%v_BC%north_t0 = neststruct%v_BC%north_t1
   neststruct%v_BC%south_t0 = neststruct%v_BC%south_t1


   neststruct%vc_BC%east_t0  = neststruct%vc_BC%east_t1
   neststruct%vc_BC%west_t0  = neststruct%vc_BC%west_t1
   neststruct%vc_BC%north_t0 = neststruct%vc_BC%north_t1
   neststruct%vc_BC%south_t0 = neststruct%vc_BC%south_t1
   neststruct%uc_BC%east_t0  = neststruct%uc_BC%east_t1
   neststruct%uc_BC%west_t0  = neststruct%uc_BC%west_t1
   neststruct%uc_BC%north_t0 = neststruct%uc_BC%north_t1
   neststruct%uc_BC%south_t0 = neststruct%uc_BC%south_t1

   neststruct%divg_BC%east_t0  = neststruct%divg_BC%east_t1
   neststruct%divg_BC%west_t0  = neststruct%divg_BC%west_t1
   neststruct%divg_BC%north_t0 = neststruct%divg_BC%north_t1
   neststruct%divg_BC%south_t0 = neststruct%divg_BC%south_t1

 end subroutine set_BCs_t0

 subroutine d2c_setup(u, v, &
      ua, va, &
      uc, vc, dord4, &
      isd,ied,jsd,jed, is,ie,js,je, npx,npy, &
      grid_type, bounded_domain, &
      se_corner, sw_corner, ne_corner, nw_corner, &
      rsin_u,rsin_v,cosa_s,rsin2 )

  logical, intent(in):: dord4
  real, intent(in) ::  u(isd:ied,jsd:jed+1)
  real, intent(in) ::  v(isd:ied+1,jsd:jed)
  real, intent(out), dimension(isd:ied  ,jsd:jed  ):: ua
  real, intent(out), dimension(isd:ied  ,jsd:jed  ):: va
  real, intent(out), dimension(isd:ied+1,jsd:jed  ):: uc
  real, intent(out), dimension(isd:ied  ,jsd:jed+1):: vc
  integer, intent(in) :: isd,ied,jsd,jed, is,ie,js,je, npx,npy,grid_type
  logical, intent(in) :: bounded_domain, se_corner, sw_corner, ne_corner, nw_corner
  real, intent(in) :: rsin_u(isd:ied+1,jsd:jed)
  real, intent(in) :: rsin_v(isd:ied,jsd:jed+1)
  real, intent(in) :: cosa_s(isd:ied,jsd:jed)
  real, intent(in) :: rsin2(isd:ied,jsd:jed)

! Local 
  real, dimension(isd:ied,jsd:jed):: utmp, vtmp
  real, parameter:: t11=27./28., t12=-13./28., t13=3./7., t14=6./7., t15=3./28.
  real, parameter:: a1 =  0.5625
  real, parameter:: a2 = -0.0625
  real, parameter:: c1 = -2./14.
  real, parameter:: c2 = 11./14.
  real, parameter:: c3 =  5./14.
  integer npt, i, j, ifirst, ilast, id

  if ( dord4) then
       id = 1
  else
       id = 0
  endif


  if (grid_type < 3 .and. .not. bounded_domain) then
     npt = 4
  else
     npt = -2
  endif

  if ( bounded_domain) then  

     do j=jsd+1,jed-1
        do i=isd,ied
           utmp(i,j) = a2*(u(i,j-1)+u(i,j+2)) + a1*(u(i,j)+u(i,j+1))
        enddo
     enddo
     do i=isd,ied
        j = jsd
        utmp(i,j) = 0.5*(u(i,j)+u(i,j+1))
        j = jed
        utmp(i,j) = 0.5*(u(i,j)+u(i,j+1))
     end do

     do j=jsd,jed
        do i=isd+1,ied-1
           vtmp(i,j) = a2*(v(i-1,j)+v(i+2,j)) + a1*(v(i,j)+v(i+1,j))
        enddo
        i = isd
        vtmp(i,j) = 0.5*(v(i,j)+v(i+1,j)) 
        i = ied
        vtmp(i,j) = 0.5*(v(i,j)+v(i+1,j))
     enddo

     do j=jsd,jed
        do i=isd,ied
           ua(i,j) = (utmp(i,j)-vtmp(i,j)*cosa_s(i,j)) * rsin2(i,j)
           va(i,j) = (vtmp(i,j)-utmp(i,j)*cosa_s(i,j)) * rsin2(i,j)
        enddo
     enddo

  else

     !----------
     ! Interior:
     !----------
     utmp = 0.
     vtmp = 0.


     do j=max(npt,js-1),min(npy-npt,je+1)
        do i=max(npt,isd),min(npx-npt,ied)
           utmp(i,j) = a2*(u(i,j-1)+u(i,j+2)) + a1*(u(i,j)+u(i,j+1))
        enddo
     enddo
     do j=max(npt,jsd),min(npy-npt,jed)
        do i=max(npt,is-1),min(npx-npt,ie+1)
           vtmp(i,j) = a2*(v(i-1,j)+v(i+2,j)) + a1*(v(i,j)+v(i+1,j))
        enddo
     enddo

     !----------
     ! edges:
     !----------
     if (grid_type < 3) then

        if ( js==1 .or. jsd<npt) then
           do j=jsd,npt-1
              do i=isd,ied
                 utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
                 vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
              enddo
           enddo
        endif

        if ( (je+1)==npy .or. jed>=(npy-npt)) then
           do j=npy-npt+1,jed
              do i=isd,ied
                 utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
                 vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
              enddo
           enddo
        endif

        if ( is==1 .or. isd<npt ) then
           do j=max(npt,jsd),min(npy-npt,jed)
              do i=isd,npt-1
                 utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
                 vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
              enddo
           enddo
        endif

        if ( (ie+1)==npx .or. ied>=(npx-npt)) then
           do j=max(npt,jsd),min(npy-npt,jed)
              do i=npx-npt+1,ied
                 utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
                 vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
              enddo
           enddo
        endif

     endif
     do j=js-1-id,je+1+id
        do i=is-1-id,ie+1+id
           ua(i,j) = (utmp(i,j)-vtmp(i,j)*cosa_s(i,j)) * rsin2(i,j)
           va(i,j) = (vtmp(i,j)-utmp(i,j)*cosa_s(i,j)) * rsin2(i,j)
        enddo
     enddo

  end if

! A -> C
!--------------
! Fix the edges
!--------------
! Xdir:
     if( sw_corner ) then
         do i=-2,0
            utmp(i,0) = -vtmp(0,1-i)
         enddo
     endif
     if( se_corner ) then
         do i=0,2
            utmp(npx+i,0) = vtmp(npx,i+1)
         enddo
     endif
     if( ne_corner ) then
         do i=0,2
            utmp(npx+i,npy) = -vtmp(npx,je-i)
         enddo
     endif
     if( nw_corner ) then
         do i=-2,0
            utmp(i,npy) = vtmp(0,je+i)
         enddo
     endif

  if (grid_type < 3 .and. .not. bounded_domain) then
     ifirst = max(3,    is-1)
     ilast  = min(npx-2,ie+2)
  else
     ifirst = is-1
     ilast  = ie+2
  endif
!---------------------------------------------
! 4th order interpolation for interior points:
!---------------------------------------------
     do j=js-1,je+1
        do i=ifirst,ilast
           uc(i,j) = a1*(utmp(i-1,j)+utmp(i,j))+a2*(utmp(i-2,j)+utmp(i+1,j))
        enddo
     enddo

     if (grid_type < 3) then
! Xdir:
     if( is==1 .and. .not. bounded_domain ) then
        do j=js-1,je+1
           uc(0,j) = c1*utmp(-2,j) + c2*utmp(-1,j) + c3*utmp(0,j) 
           uc(1,j) = ( t14*(utmp( 0,j)+utmp(1,j))    &
                     + t12*(utmp(-1,j)+utmp(2,j))    &
                     + t15*(utmp(-2,j)+utmp(3,j)) )*rsin_u(1,j)
           uc(2,j) = c1*utmp(3,j) + c2*utmp(2,j) + c3*utmp(1,j)
        enddo
     endif

     if( (ie+1)==npx .and. .not. bounded_domain ) then
        do j=js-1,je+1
           uc(npx-1,j) = c1*utmp(npx-3,j)+c2*utmp(npx-2,j)+c3*utmp(npx-1,j) 
           uc(npx,j) = (t14*(utmp(npx-1,j)+utmp(npx,j))+      &
                        t12*(utmp(npx-2,j)+utmp(npx+1,j))     &
                      + t15*(utmp(npx-3,j)+utmp(npx+2,j)))*rsin_u(npx,j)
           uc(npx+1,j) = c3*utmp(npx,j)+c2*utmp(npx+1,j)+c1*utmp(npx+2,j) 
        enddo
     endif

     endif

!------
! Ydir:
!------
     if( sw_corner ) then
         do j=-2,0
            vtmp(0,j) = -utmp(1-j,0)
         enddo
     endif
     if( nw_corner ) then
         do j=0,2
            vtmp(0,npy+j) = utmp(j+1,npy)
         enddo
     endif
     if( se_corner ) then
         do j=-2,0
            vtmp(npx,j) = utmp(ie+j,0)
         enddo
     endif
     if( ne_corner ) then
         do j=0,2
            vtmp(npx,npy+j) = -utmp(ie-j,npy)
         enddo
     endif

     if (grid_type < 3) then

     do j=js-1,je+2
      if ( j==1  .and. .not. bounded_domain) then
        do i=is-1,ie+1
           vc(i,1) = (t14*(vtmp(i, 0)+vtmp(i,1))    &
                    + t12*(vtmp(i,-1)+vtmp(i,2))    &
                    + t15*(vtmp(i,-2)+vtmp(i,3)))*rsin_v(i,1)
        enddo
      elseif ( (j==0 .or. j==(npy-1))  .and. .not. bounded_domain) then
        do i=is-1,ie+1
           vc(i,j) = c1*vtmp(i,j-2) + c2*vtmp(i,j-1) + c3*vtmp(i,j)
        enddo
      elseif ( (j==2 .or. j==(npy+1))  .and. .not. bounded_domain) then
        do i=is-1,ie+1
           vc(i,j) = c1*vtmp(i,j+1) + c2*vtmp(i,j) + c3*vtmp(i,j-1)
        enddo
      elseif ( j==npy  .and. .not. bounded_domain) then
        do i=is-1,ie+1
           vc(i,npy) = (t14*(vtmp(i,npy-1)+vtmp(i,npy))    &
                      + t12*(vtmp(i,npy-2)+vtmp(i,npy+1))  &
                      + t15*(vtmp(i,npy-3)+vtmp(i,npy+2)))*rsin_v(i,npy)
        enddo
      else
! 4th order interpolation for interior points:
        do i=is-1,ie+1
           vc(i,j) = a2*(vtmp(i,j-2)+vtmp(i,j+1))+a1*(vtmp(i,j-1)+vtmp(i,j))
        enddo
     endif
     enddo
    else
! 4th order interpolation:
       do j=js-1,je+2
          do i=is-1,ie+1
             vc(i,j) = a2*(vtmp(i,j-2)+vtmp(i,j+1))+a1*(vtmp(i,j-1)+vtmp(i,j))
          enddo
       enddo
    endif

  end subroutine d2c_setup

 subroutine d2a_setup(u, v, ua, va, dord4, &
      isd,ied,jsd,jed, is,ie,js,je, npx,npy, &
      grid_type, bounded_domain, &
      cosa_s,rsin2 )

  logical, intent(in):: dord4
  real, intent(in) ::  u(isd:ied,jsd:jed+1)
  real, intent(in) ::  v(isd:ied+1,jsd:jed)
  real, intent(out), dimension(isd:ied  ,jsd:jed  ):: ua
  real, intent(out), dimension(isd:ied  ,jsd:jed  ):: va
  integer, intent(in) :: isd,ied,jsd,jed, is,ie,js,je, npx,npy,grid_type
  real, intent(in) :: cosa_s(isd:ied,jsd:jed)
  real, intent(in) :: rsin2(isd:ied,jsd:jed)
  logical, intent(in) :: bounded_domain

! Local 
  real, dimension(isd:ied,jsd:jed):: utmp, vtmp
  real, parameter:: t11=27./28., t12=-13./28., t13=3./7., t14=6./7., t15=3./28.
  real, parameter:: a1 =  0.5625
  real, parameter:: a2 = -0.0625
  real, parameter:: c1 = -2./14.
  real, parameter:: c2 = 11./14.
  real, parameter:: c3 =  5./14.
  integer npt, i, j, ifirst, ilast, id

  if ( dord4) then
       id = 1
  else
       id = 0
  endif


  if (grid_type < 3 .and. .not. bounded_domain) then
     npt = 4
  else
     npt = -2
  endif

  if ( bounded_domain) then  

     do j=jsd+1,jed-1
        do i=isd,ied
           utmp(i,j) = a2*(u(i,j-1)+u(i,j+2)) + a1*(u(i,j)+u(i,j+1))
        enddo
     enddo
     do i=isd,ied
        j = jsd
        utmp(i,j) = 0.5*(u(i,j)+u(i,j+1))
        j = jed
        utmp(i,j) = 0.5*(u(i,j)+u(i,j+1))
     end do

     do j=jsd,jed
        do i=isd+1,ied-1
           vtmp(i,j) = a2*(v(i-1,j)+v(i+2,j)) + a1*(v(i,j)+v(i+1,j))
        enddo
        i = isd
        vtmp(i,j) = 0.5*(v(i,j)+v(i+1,j)) 
        i = ied
        vtmp(i,j) = 0.5*(v(i,j)+v(i+1,j))
     enddo

  else

     !----------
     ! Interior:
     !----------

     do j=max(npt,js-1),min(npy-npt,je+1)
        do i=max(npt,isd),min(npx-npt,ied)
           utmp(i,j) = a2*(u(i,j-1)+u(i,j+2)) + a1*(u(i,j)+u(i,j+1))
        enddo
     enddo
     do j=max(npt,jsd),min(npy-npt,jed)
        do i=max(npt,is-1),min(npx-npt,ie+1)
           vtmp(i,j) = a2*(v(i-1,j)+v(i+2,j)) + a1*(v(i,j)+v(i+1,j))
        enddo
     enddo

     !----------
     ! edges:
     !----------
     if (grid_type < 3) then

        if ( js==1 .or. jsd<npt) then
           do j=jsd,npt-1
              do i=isd,ied
                 utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
                 vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
              enddo
           enddo
        endif

        if ( (je+1)==npy .or. jed>=(npy-npt)) then
           do j=npy-npt+1,jed
              do i=isd,ied
                 utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
                 vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
              enddo
           enddo
        endif

        if ( is==1 .or. isd<npt ) then
           do j=max(npt,jsd),min(npy-npt,jed)
              do i=isd,npt-1
                 utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
                 vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
              enddo
           enddo
        endif

        if ( (ie+1)==npx .or. ied>=(npx-npt)) then
           do j=max(npt,jsd),min(npy-npt,jed)
              do i=npx-npt+1,ied
                 utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
                 vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
              enddo
           enddo
        endif

     endif

  end if



  do j=js-1-id,je+1+id
     do i=is-1-id,ie+1+id
        ua(i,j) = (utmp(i,j)-vtmp(i,j)*cosa_s(i,j)) * rsin2(i,j)
        va(i,j) = (vtmp(i,j)-utmp(i,j)*cosa_s(i,j)) * rsin2(i,j)
     enddo
  enddo

end subroutine d2a_setup




!! nestupdate types
!! 1 - Interpolation update on all variables
!! 2 - Conserving update (over areas on cell-
!!     centered variables, over faces on winds) on all variables
!! 3 - Interpolation update on winds only
!! 4 - Interpolation update on all variables except delp (mass conserving)
!! 5 - Remap interpolating update, delp not updated
!! 6 - Remap conserving update, delp not updated
!! 7 - Remap conserving update, delp and q not updated
!! 8 - Remap conserving update, only winds updated

!! Note that nestupdate > 3 will not update delp.

!! "Remap update" remaps updated variables from the nested grid's
!!  vertical coordinate to that of the coarse grid. When delp is not
!!  updated (nestbctype >= 3) the vertical coordinates differ on
!!  the two grids, because the surface pressure will be different
!!  on the two grids.
!! Note: "conserving updates" do not guarantee global conservation
!!  unless flux nested grid BCs are specified, or if a quantity is
!!  not updated at all. This ability has not been implemented.

!>@brief The subroutine'twoway_nesting' performs a two-way update 
!! of nested-grid data onto the parent grid.
subroutine twoway_nesting(Atm, ngrids, grids_on_this_pe, zvir, Time, this_grid)

   type(fv_atmos_type), intent(INOUT) :: Atm(ngrids)
   integer, intent(IN) :: ngrids, this_grid
   logical, intent(IN) :: grids_on_this_pe(ngrids)
   real, intent(IN) :: zvir
   type(time_type), intent(IN) :: Time

   integer :: n, p, sphum

   
   if (ngrids > 1) then

! Re-compute pressures on each grid

      call p_var(Atm(this_grid)%npz, Atm(this_grid)%bd%is, Atm(this_grid)%bd%ie, Atm(this_grid)%bd%js, Atm(this_grid)%bd%je, &
           Atm(this_grid)%ptop, ptop_min, Atm(this_grid)%delp, Atm(this_grid)%delz, Atm(this_grid)%pt, &
           Atm(this_grid)%ps, Atm(this_grid)%pe, Atm(this_grid)%peln, Atm(this_grid)%pk, Atm(this_grid)%pkz, kappa, &
           Atm(this_grid)%q, Atm(this_grid)%ng, Atm(this_grid)%flagstruct%ncnst,  Atm(this_grid)%gridstruct%area_64, 0.,  &
           .false.,  .false., & 
           Atm(this_grid)%flagstruct%moist_phys,  Atm(this_grid)%flagstruct%hydrostatic, &
           Atm(this_grid)%flagstruct%nwat, Atm(this_grid)%domain, Atm(this_grid)%flagstruct%adiabatic, .false.)

      do n=ngrids,2,-1 !loop backwards to allow information to propagate from finest to coarsest grids

         !two-way updating    
         if (Atm(n)%neststruct%twowaynest ) then
            !if  (grids_on_this_pe(n) .or. grids_on_this_pe(Atm(n)%parent_grid%grid_number)) then
            if (n==this_grid .or. Atm(n)%parent_grid%grid_number==this_grid) then
               sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
               call twoway_nest_update(Atm(n)%npx, Atm(n)%npy, Atm(n)%npz, zvir, &
                    Atm(n)%ncnst, sphum, Atm(n)%u, Atm(n)%v, Atm(n)%w, &
                    Atm(n)%pt, Atm(n)%delp, Atm(n)%q, &
                    Atm(n)%pe, Atm(n)%pkz, Atm(n)%delz, Atm(n)%ps, Atm(n)%ptop, Atm(n)%ak, Atm(n)%bk, &
                    Atm(n)%gridstruct, Atm(n)%flagstruct, Atm(n)%neststruct, Atm(n)%domain, &
                    Atm(n)%parent_grid, Atm(N)%bd, n, .false.)
            endif
         endif

      end do

      !NOTE: these routines need to be used with any grid which has been updated to, not just the coarsest grid.
      if (Atm(this_grid)%neststruct%parent_of_twoway .and. grids_on_this_pe(n)) then
            call after_twoway_nest_update( Atm(this_grid)%npx, Atm(this_grid)%npy, Atm(this_grid)%npz, &
                 Atm(this_grid)%ng,     Atm(this_grid)%ncnst,   &
                 Atm(this_grid)%u,      Atm(this_grid)%v,     Atm(this_grid)%w,    Atm(this_grid)%delz, &
                 Atm(this_grid)%pt,     Atm(this_grid)%delp,  Atm(this_grid)%q,   &
                 Atm(this_grid)%ps,     Atm(this_grid)%pe,    Atm(this_grid)%pk,   Atm(this_grid)%peln,  Atm(this_grid)%pkz, &
                 Atm(this_grid)%phis,   Atm(this_grid)%ua,    Atm(this_grid)%va,  &
                 Atm(this_grid)%ptop,   Atm(this_grid)%gridstruct, Atm(this_grid)%flagstruct, &
                 Atm(this_grid)%domain, Atm(this_grid)%bd, Time)
      endif

   endif ! ngrids > 1

  end subroutine twoway_nesting

!!!CLEANUP: this routine assumes that the PARENT GRID has pt = (regular) temperature,
!!!not potential temperature; which may cause problems when updating if this is not the case.

!!! NOTE ALSO: parent_grid%flagstruct is NOT SET UP by default and may be missing much information
!!! Either make sure that parent_grid%flagstruct is filled in fv_control or that proper steps
!!!   are taken to make sure null flags are not used
 subroutine twoway_nest_update(npx, npy, npz, zvir, ncnst, sphum,     &
                        u, v, w, pt, delp, q,   &
                        pe, pkz, delz, ps, ptop, ak, bk, &
                        gridstruct, flagstruct, neststruct, &
                        domain, parent_grid, bd, grid_number, conv_theta_in)

    real, intent(IN) :: zvir, ptop, ak(npz+1), bk(npz+1)

    integer, intent(IN) :: npx, npy, npz
    integer, intent(IN) :: ncnst, sphum, grid_number
    logical, intent(IN), OPTIONAL :: conv_theta_in

    type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(inout), dimension(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) :: u !< D grid zonal wind (m/s)
    real, intent(inout), dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) :: v !< D grid meridional wind (m/s)
    real, intent(inout) :: w(   bd%isd:        ,bd%jsd:        ,1: )  !<  W (m/s)

    real, intent(inout) :: pt(  bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  !< temperature (K)
    real, intent(inout) :: delp(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  !< pressure thickness (pascal)
    real, intent(inout) :: q(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz, ncnst) !< specific humidity and constituents
  
    real, intent(inout) :: pe  (bd%is-1:bd%ie+1,npz+1,bd%js-1:bd%je+1)  !< finite-volume interface p ! NOTE TRANSPOSITION NEEDED
    real, intent(inout) :: pkz (bd%is:bd%ie,bd%js:bd%je,npz)       !< finite-volume mean pk
    real, intent(inout) :: delz(bd%isd:      ,bd%jsd:      ,1: )   !< delta-height (m); non-hydrostatic only
    real, intent(inout) :: ps  (bd%isd:bd%ied  ,bd%jsd:bd%jed)     !< Surface pressure (pascal)

    type(fv_grid_type), intent(INOUT) :: gridstruct
    type(fv_flags_type), intent(INOUT) :: flagstruct
    type(fv_nest_type), intent(INOUT) :: neststruct
    type(domain2d), intent(INOUT) :: domain

    type(fv_atmos_type), pointer, intent(IN) :: parent_grid

    real, allocatable :: t_nest(:,:,:), ps0(:,:)
    integer :: i,j,k,n
    integer :: isd_p, ied_p, jsd_p, jed_p, isc_p, iec_p, jsc_p, jec_p
    integer :: isg, ieg, jsg,jeg, npx_p, npy_p
    integer :: istart, iend
    real :: qmass_b, qmass_a, fix = 1.
    logical :: used, conv_theta=.true.

    real :: qdp(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)
    real, allocatable, dimension(:,:,:) :: qdp_coarse
    real, allocatable, dimension(:,:,:) :: var_src
    real, allocatable, dimension(:,:,:) :: pt_src, w_src, u_src, v_src
    real(kind=f_p), allocatable :: q_diff(:,:,:)
    real :: L_sum_b(npz), L_sum_a(npz), blend_wt(parent_grid%npz)
    real :: pfull, ph1, ph2, rfcut, sgcut
    
    integer :: upoff
    integer :: is,  ie,  js,  je
    integer :: isd, ied, jsd, jed
    integer :: isu, ieu, jsu, jeu
    logical, SAVE :: first_timestep = .true.

    is  = bd%is
    ie  = bd%ie
    js  = bd%js
    je  = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    isu = neststruct%isu
    ieu = neststruct%ieu
    jsu = neststruct%jsu
    jeu = neststruct%jeu

    upoff = neststruct%upoff

    !We update actual temperature, not theta.
    !If pt is actual temperature, set conv_theta to .false.
    if (present(conv_theta_in)) conv_theta = conv_theta_in

    if ((.not. parent_grid%neststruct%parent_proc) .and. (.not. neststruct%child_proc)) return

    call mpp_get_data_domain( parent_grid%domain, &
         isd_p,  ied_p,  jsd_p,  jed_p  )
    call mpp_get_compute_domain( parent_grid%domain, &
         isc_p,  iec_p,  jsc_p,  jec_p  )

    ph2 = parent_grid%ak(1)
    rfcut = max(flagstruct%rf_cutoff, parent_grid%flagstruct%rf_cutoff)
    sgcut = ak(flagstruct%n_sponge+1) + bk(flagstruct%n_sponge+1)*flagstruct%p_ref
    sgcut = max(sgcut, parent_grid%ak(parent_grid%flagstruct%n_sponge+1) + parent_grid%bk(parent_grid%flagstruct%n_sponge+1)*parent_grid%flagstruct%p_ref)
    rfcut = max(rfcut, sgcut)
    do k=1,parent_grid%npz
       ph1 = ph2
       ph2 = parent_grid%ak(k+1) + parent_grid%bk(k+1)*parent_grid%flagstruct%p_ref
       pfull = (ph2 - ph1) / log(ph2/ph1)
       !if above nested-grid ptop or top two nested-grid levels do not remap
       if ( pfull <= ak(3) .or. k <= 2 ) then
          blend_wt(k) = 0.
       !Partial blend of nested-grid's Rayleigh damping region
       !ALSO do not blend n_sponge areas??
       elseif (pfull <= rfcut) then
          blend_wt(k) = 0.
          !blend_wt(k) = neststruct%update_blend*cos(0.5*pi*log(rfcut/pfull)/log(rfcut/ptop))**2
       else
          blend_wt(k) = neststruct%update_blend
       endif
    enddo

    if (parent_grid%neststruct%parent_proc .and. is_master() .and. first_timestep) then
       print*, ' TWO-WAY BLENDING WEIGHTS'
       ph2 = parent_grid%ak(1)
       do k=1,parent_grid%npz
          ph1 = ph2
          ph2 = parent_grid%ak(k+1) + parent_grid%bk(k+1)*parent_grid%flagstruct%p_ref
          pfull = (ph2 - ph1) / log(ph2/ph1)
          print*, k, pfull, blend_wt(k)
       enddo
       first_timestep = .false.
    endif

#ifndef SW_DYNAMICS
   if (neststruct%nestupdate /= 3 .and. neststruct%nestupdate /= 8) then

      if (neststruct%child_proc) then
         call mpp_update_domains(ps, domain, complete=.true.)
         if (.not. flagstruct%hydrostatic) call mpp_update_domains(w, domain)
         ! if (neststruct%child_proc)  call mpp_update_domains(delz, domain)
         call mpp_update_domains(u, v, domain, gridtype=DGRID_NE)
      endif
      allocate(pt_src(isd_p:ied_p,jsd_p:jed_p,npz))
      pt_src = -999.

      if (conv_theta) then

         if (neststruct%child_proc) then
            !pt is potential temperature on the nested grid, but actual
            !temperature on the coarse grid. Compute actual temperature
            !on the nested grid, then gather.
            allocate(t_nest(isd:ied,jsd:jed,1:npz))
!$OMP parallel do default(none) shared(npz,js,je,is,ie,t_nest,pt,pkz,zvir,q,sphum)
            do k=1,npz
               do j=js,je
                  do i=is,ie
                     t_nest(i,j,k) = pt(i,j,k)*pkz(i,j,k)/(1.+zvir*q(i,j,k,sphum))
                  enddo
               enddo
            enddo
            call mpp_update_domains(t_nest, domain, complete=.true.)
         endif

         call update_coarse_grid(pt_src, &
              t_nest, global_nest_domain, &
              gridstruct%dx, gridstruct%dy, gridstruct%area, &
              bd, isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, &
              neststruct%isu, neststruct%ieu, neststruct%jsu, neststruct%jeu, &
              npx, npy, npz, 0, 0, &
              neststruct%refinement, neststruct%nestupdate, upoff, 0, &
              parent_grid%neststruct%parent_proc, neststruct%child_proc, parent_grid, grid_number-1)
         if (neststruct%child_proc)  deallocate(t_nest)
      else
         if (neststruct%child_proc)  call mpp_update_domains(pt, domain, complete=.true.)

         call update_coarse_grid(pt_src, &
              pt, global_nest_domain, &
              gridstruct%dx, gridstruct%dy, gridstruct%area, &
              bd, isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, &
              neststruct%isu, neststruct%ieu, neststruct%jsu, neststruct%jeu, &
              npx, npy, npz, 0, 0, &
              neststruct%refinement, neststruct%nestupdate, upoff, 0, &
              parent_grid%neststruct%parent_proc, neststruct%child_proc, parent_grid, grid_number-1)

      endif !conv_theta

      call mpp_sync!self


      !We don't currently have a good way to communicate all namelist items between
      ! grids (since we cannot assume that we have internal namelists available), so 
      ! we get the clutzy structure here.
      if ( (neststruct%child_proc .and. .not. flagstruct%hydrostatic) .or. &
           (parent_grid%neststruct%parent_proc .and. .not. parent_grid%flagstruct%hydrostatic) ) then

         allocate(w_src(isd_p:ied_p,jsd_p:jed_p,npz))
         w_src = -999.
         call update_coarse_grid(w_src, w, global_nest_domain, &
              gridstruct%dx, gridstruct%dy, gridstruct%area, &
              bd, isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, &
              neststruct%isu, neststruct%ieu, neststruct%jsu, neststruct%jeu, &
              npx, npy, npz, 0, 0, &
              neststruct%refinement, neststruct%nestupdate, upoff, 0, &
              parent_grid%neststruct%parent_proc, neststruct%child_proc, parent_grid, grid_number-1)
         call mpp_sync!self

            !Updating for delz not yet implemented; 
            ! may need to think very carefully how one would do this!!!
            ! consider updating specific volume instead?
!!$            call update_coarse_grid(parent_grid%delz, delz, global_nest_domain, &
!!$                 bd, isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, npz, 0, 0, &
!!$                 neststruct%refinement, neststruct%nestupdate, upoff, 0, parent_grid%neststruct%parent_proc, neststruct%child_proc)

      end if
      
   end if !Neststruct%nestupdate /= 3

#endif

   allocate(u_src(isd_p:ied_p,  jsd_p:jed_p+1,npz))
   allocate(v_src(isd_p:ied_p+1,jsd_p:jed_p,npz))
   u_src = -999.
   v_src = -999.
   call update_coarse_grid(u_src, v_src, u, v, global_nest_domain, &
        gridstruct%dx, gridstruct%dy, gridstruct%area, &
        bd, isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, &
        neststruct%isu, neststruct%ieu, neststruct%jsu, neststruct%jeu, &
        npx, npy, npz, 0, 1, 1, 0, &
        neststruct%refinement, neststruct%nestupdate, upoff, 0, &
        parent_grid%neststruct%parent_proc, neststruct%child_proc, parent_grid, grid_number-1, gridtype=DGRID_NE)

   call mpp_sync()

#ifndef SW_DYNAMICS
   if (neststruct%nestupdate >= 5 .and. npz > 4) then

      !Use PS0 from nested grid, NOT the full delp. Also we assume the same number of levels on both grids.
      !PS0 should be initially set to be ps so that this routine does NOTHING outside of the update region

      !Re-compute nested (AND COARSE) grid ps

      allocate(ps0(isd_p:ied_p,jsd_p:jed_p))
      if (parent_grid%neststruct%parent_proc) then

         parent_grid%ps = parent_grid%ptop
!$OMP parallel do default(none) shared(jsd_p,jed_p,isd_p,ied_p,parent_grid)
         do j=jsd_p,jed_p
            do k=1,parent_grid%npz
               do i=isd_p,ied_p
                  parent_grid%ps(i,j) = parent_grid%ps(i,j) + &
                       parent_grid%delp(i,j,k)
               end do
            end do
         end do

         ps0 = parent_grid%ps
      endif

      if (neststruct%child_proc) then

         ps = ptop
!$OMP parallel do default(none) shared(npz,jsd,jed,isd,ied,ps,delp)
         do j=jsd,jed
            do k=1,npz
               do i=isd,ied
                  ps(i,j) = ps(i,j) + delp(i,j,k)
               end do
            end do
         end do
      endif

      call update_coarse_grid(ps0, ps, global_nest_domain, &
              gridstruct%dx, gridstruct%dy, gridstruct%area, &
              bd, isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, &
              neststruct%isu, neststruct%ieu, neststruct%jsu, neststruct%jeu, &
              npx, npy, 0, 0, &
              neststruct%refinement, neststruct%nestupdate, upoff, 0, parent_grid%neststruct%parent_proc, neststruct%child_proc, parent_grid, grid_number-1)

      !!! The mpp version of update_coarse_grid does not return a consistent value of ps
      !!! across PEs, as it does not go into the haloes of a given coarse-grid PE. This
      !!! update_domains call takes care of the problem.

      if (parent_grid%neststruct%parent_proc) then
         call mpp_update_domains(parent_grid%ps, parent_grid%domain, complete=.false.)
         call mpp_update_domains(ps0, parent_grid%domain, complete=.true.)
      endif

      call mpp_sync!self

      if (parent_grid%global_tile == neststruct%parent_tile) then 

         if (parent_grid%neststruct%parent_proc) then

         !comment out if statement to always remap theta instead of t in the remap-update.
         !(In LtE typically we use remap_t = .true.: remapping t is better (except in
         !idealized simulations with a background uniform theta) since near the top
         !boundary theta is exponential, which is hard to accurately interpolate with a spline
         if (.not. parent_grid%flagstruct%remap_t) then
!$OMP parallel do default(none) shared(jsc_p,jec_p,isc_p,iec_p,parent_grid,zvir,sphum)
            do k=1,parent_grid%npz
               do j=jsc_p,jec_p
                  do i=isc_p,iec_p
                     parent_grid%pt(i,j,k) = &
                          parent_grid%pt(i,j,k)/parent_grid%pkz(i,j,k)*&
                          (1.+zvir*parent_grid%q(i,j,k,sphum))
                  end do
               end do
            end do
         end if
!!$!!!! DEBUG CODE
!!$         do k=1,parent_grid%npz
!!$            write(mpp_pe()+3000,*) 'k = ', k, parent_grid%ak(k), parent_grid%bk(k)
!!$         enddo
!!$         write(mpp_pe()+3000,*) 
!!$         do k=1,npz
!!$            write(mpp_pe()+3000,*) 'k = ', k, ak(k), bk(k)
!!$         enddo
!!$!!!! END DEBUG CODE

         call update_remap_tqw(parent_grid%npz, parent_grid%ak, parent_grid%bk, &
              parent_grid%ps, &
              parent_grid%pt, parent_grid%q, parent_grid%w, &
              parent_grid%flagstruct%hydrostatic, &
              npz, ps0, ak, bk, pt_src, w_src, &
              zvir, parent_grid%ptop, ncnst, &
              parent_grid%flagstruct%kord_tm, parent_grid%flagstruct%kord_tr, &
              parent_grid%flagstruct%kord_wz, &
              isc_p, iec_p, jsc_p, jec_p, isd_p, ied_p, jsd_p, jed_p, .false., &
              neststruct%isu, neststruct%ieu, neststruct%jsu, neststruct%jeu, blend_wt) !neststruct%nestupdate < 7)
         if (.not. parent_grid%flagstruct%remap_t) then
!$OMP parallel do default(none) shared(jsc_p,jec_p,isc_p,iec_p,parent_grid,zvir,sphum)
            do k=1,parent_grid%npz
               do j=jsc_p,jec_p
                  do i=isc_p,iec_p
                     parent_grid%pt(i,j,k) = &
                          parent_grid%pt(i,j,k)*parent_grid%pkz(i,j,k) / &
                          (1.+zvir*parent_grid%q(i,j,k,sphum))
                  end do
               end do
            end do
         end if

         call update_remap_uv(parent_grid%npz, parent_grid%ak, parent_grid%bk, &
              parent_grid%ps, parent_grid%u, parent_grid%v, &
              npz, ak, bk, ps0, u_src, v_src, &
              parent_grid%flagstruct%kord_mt, &
              isc_p, iec_p, jsc_p, jec_p, isd_p, ied_p, jsd_p, jed_p, parent_grid%ptop, &
              neststruct%isu, neststruct%ieu, neststruct%jsu, neststruct%jeu, blend_wt)

         endif !parent_grid%neststruct%parent_proc

      end if

      if (allocated(ps0)) deallocate(ps0)

   end if

#endif



   deallocate(pt_src)
   deallocate(w_src)
   deallocate(u_src)
   deallocate(v_src)


 end subroutine twoway_nest_update

 subroutine level_sum(q, area, domain, bd, npz, L_sum)

    integer, intent(IN) :: npz
    type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(in) :: area(   bd%isd:bd%ied  ,bd%jsd:bd%jed)
    real, intent(in) ::    q(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)
    real, intent(OUT) :: L_sum( npz ) 
    type(domain2d), intent(IN) :: domain
   
    integer :: i, j, k, n
    real :: qA!(bd%is:bd%ie, bd%js:bd%je)

    do k=1,npz
       qA = 0.
       do j=bd%js,bd%je
       do i=bd%is,bd%ie
          !qA(i,j) = q(i,j,k)*area(i,j)
          qA = qA + q(i,j,k)*area(i,j)
       enddo
       enddo
       call mp_reduce_sum(qA)
       L_sum(k) = qA
!       L_sum(k) = mpp_global_sum(domain, qA, flags=BITWISE_EXACT_SUM)
!       L_sum(k) = mpp_global_sum(domain, qA, flags=BITWISE_EFP_SUM) ! doesn't work??
    enddo

 end subroutine level_sum

![ij]start and [ij]end should already take staggering into account
!!! CHECK ARRAY BOUNDS!!
!! Make sure data is in the correct place.
 subroutine remap_up_k(ps_src, ps_dst, ak_src, bk_src, ak_dst, bk_dst, var_src, var_dst, &
      bd, istart, iend, jstart, jend, istag, jstag, npz_src, npz_dst, iv, kord, blend_wt, log_pe)

   !Note here that pe is TRANSPOSED to make loops faster
   type(fv_grid_bounds_type), intent(IN) :: bd
   integer, intent(IN) :: istart, iend, jstart, jend, npz_dst, npz_src, iv, kord, istag, jstag
   logical, intent(IN) :: log_pe
   real, intent(INOUT) :: ps_src(bd%isd:bd%ied,bd%jsd:bd%jed), var_src(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,npz_src)
   real, intent(INOUT) :: ps_dst(bd%isd:bd%ied,bd%jsd:bd%jed), var_dst(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,npz_dst)
   real, intent(IN) :: blend_wt(npz_dst), ak_src(npz_src+1), bk_src(npz_src+1), ak_dst(npz_dst+1), bk_dst(npz_dst+1)

   integer :: i, j, k
   real pe_src(istart:iend,npz_src+1)
   real pe_dst(istart:iend,npz_dst+1)
   real peln_src(istart:iend,npz_src+1)
   real peln_dst(istart:iend,npz_dst+1)
   character(120) :: errstring
   real var_dst_unblend(istart:iend,npz_dst)
   real bw1, bw2
   
   if (iend < istart) return
   if (jend < jstart) return

!!$!!!! DEBUG CODE
!!$      write(debug_unit,*) bd%isd,bd%ied,bd%jsd,bd%jed
!!$      write(debug_unit,*) istart,iend,jstart,jend,istag,jstag
!!$      write(debug_unit,*)
!!$!!! END DEBUG CODE
   
   
   !Compute Eulerian pressures
   !NOTE: assumes that istag + jstag <= 1
   if (istag > 0) then
!$OMP parallel do default(none) shared(istart,iend,jstart,jend,npz_src,npz_dst,pe_src,ak_src,ps_src,bk_src,pe_dst,ak_dst,ps_dst,bk_dst)
      do j=jstart,jend
      do k=1,npz_src+1
      do i=istart,iend
         pe_src(i,k) = ak_src(k) + 0.5*(ps_src(i,j)+ps_src(i-1,j))*bk_src(k)
      enddo
      enddo
      do k=1,npz_dst+1
      do i=istart,iend
         pe_dst(i,k) = ak_dst(k) + 0.5*(ps_dst(i,j)+ps_dst(i-1,j))*bk_dst(k)
      enddo
      enddo
      enddo
   elseif (jstag > 0) then
!$OMP parallel do default(none) shared(istart,iend,jstart,jend,npz_src,npz_dst,pe_src,ak_src,ps_src,bk_src,pe_dst,ak_dst,ps_dst,bk_dst)
      do j=jstart,jend
      do k=1,npz_src+1
      do i=istart,iend
         pe_src(i,k) = ak_src(k) + 0.5*(ps_src(i,j)+ps_src(i,j-1))*bk_src(k)
      enddo
      enddo
      do k=1,npz_dst+1
      do i=istart,iend
         pe_dst(i,k) = ak_dst(k) + 0.5*(ps_dst(i,j)+ps_dst(i,j-1))*bk_dst(k)
      enddo
      enddo
      enddo
   else
!$OMP parallel do default(none) shared(istart,iend,jstart,jend,npz_src,npz_dst,pe_src,ak_src,ps_src,bk_src,pe_dst,ak_dst,ps_dst,bk_dst)
      do j=jstart,jend
      do k=1,npz_src+1
      do i=istart,iend
         pe_src(i,k) = ak_src(k) + ps_src(i,j)*bk_src(k)
      enddo
      enddo
      do k=1,npz_dst+1
      do i=istart,iend
         pe_dst(i,k) = ak_dst(k) + ps_dst(i,j)*bk_dst(k)
      enddo
      enddo
      enddo
   endif
      
   if (log_pe) then

!$OMP parallel do default(none) shared(istart,iend,jstart,jend,npz_src,npz_dst,pe_src,pe_dst,var_src,var_dst,iv,kord,blend_wt) &
!$OMP                          private(peln_src,peln_dst,bw1,bw2,var_dst_unblend)
      do j=jstart,jend

         do k=1,npz_src+1
         do i=istart,iend
            peln_src(i,k) = log(pe_src(i,k))
         enddo
         enddo

         do k=1,npz_dst+1
         do i=istart,iend
            peln_dst(i,k) = log(pe_dst(i,k))
         enddo
         enddo

         !remap_2d seems to have some bugs when doing logp remapping
         call    mappm(npz_src, peln_src, var_src(istart:iend,j:j,:), &
                       npz_dst, peln_dst, var_dst_unblend, &
                       istart, iend, iv, kord, peln_dst(istart,1))

         do k=1,npz_dst
            bw1 = blend_wt(k)
            bw2 = 1. - bw1
         do i=istart,iend
            var_dst(i,j,k) = var_dst(i,j,k)*bw2 + var_dst_unblend(i,k)*bw1
         enddo
         enddo
      enddo

   else

!$OMP parallel do default(none) shared(istart,iend,jstart,jend,npz_src,npz_dst,pe_src,pe_dst,var_src,var_dst,iv,kord,blend_wt) &
!$OMP                          private(bw1,bw2,var_dst_unblend)
      do j=jstart,jend

         call mappm(npz_src, pe_src, var_src(istart:iend,j:j,:), &
                    npz_dst, pe_dst, var_dst_unblend, &
                    istart, iend, iv, kord, pe_dst(istart,1))
         
         do k=1,npz_dst
            bw1 = blend_wt(k)
            bw2 = 1. - bw1
         do i=istart,iend
            var_dst(i,j,k) = var_dst(i,j,k)*bw2 + var_dst_unblend(i,k)*bw1
         enddo
         enddo
      enddo

   endif

 end subroutine remap_up_k

 subroutine after_twoway_nest_update(npx, npy, npz, ng, ncnst,   &
                        u, v, w, delz, pt, delp, q,              &
                        ps, pe, pk, peln, pkz, phis, ua, va,     &
                        ptop, gridstruct, flagstruct,            &
                        domain, bd, Time)

   type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(IN) :: ptop

    integer, intent(IN) :: ng, npx, npy, npz
    integer, intent(IN) :: ncnst

    real, intent(inout), dimension(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) :: u !< D grid zonal wind (m/s)
    real, intent(inout), dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) :: v !< D grid meridional wind (m/s)
    real, intent(inout) :: w(   bd%isd:        ,bd%jsd:        ,1: )  !< W (m/s)
    real, intent(inout) :: pt(  bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  !< temperature (K)
    real, intent(inout) :: delp(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  !< pressure thickness (pascal)
    real, intent(inout) :: q(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz, ncnst) !< specific humidity and constituents
    real, intent(inout) :: delz(bd%is:        ,bd%js:        ,1: )   !< delta-height (m); non-hydrostatic only

!-----------------------------------------------------------------------
! Auxilliary pressure arrays:    
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
    real, intent(inout) :: ps  (bd%isd:bd%ied  ,bd%jsd:bd%jed)           !< Surface pressure (pascal)
    real, intent(inout) :: pe  (bd%is-1:bd%ie+1, npz+1,bd%js-1:bd%je+1)  !< edge pressure (pascal)
    real, intent(inout) :: pk  (bd%is:bd%ie,bd%js:bd%je, npz+1)          !< pe**cappa
    real, intent(inout) :: peln(bd%is:bd%ie,npz+1,bd%js:bd%je)           !< ln(pe)
    real, intent(inout) :: pkz (bd%is:bd%ie,bd%js:bd%je,npz)             !< finite-volume mean pk
    
!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
    real, intent(inout) :: phis(bd%isd:bd%ied,bd%jsd:bd%jed)       !< Surface geopotential (g*Z_surf)

    real, intent(inout), dimension(bd%isd:bd%ied ,bd%jsd:bd%jed ,npz):: ua, va
    type(fv_grid_type), intent(IN) :: gridstruct
    type(fv_flags_type), intent(IN) :: flagstruct
    type(domain2d), intent(INOUT) :: domain
    type(time_type), intent(IN) :: Time

    logical :: bad_range

    integer :: is,  ie,  js,  je
    integer :: isd, ied, jsd, jed
    
    is  = bd%is
    ie  = bd%ie
    js  = bd%js
    je  = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed

    call cubed_to_latlon(u, v, ua, va, &
         gridstruct, npx, npy, npz, &
         1, gridstruct%grid_type, domain, &
         gridstruct%bounded_domain, flagstruct%c2l_ord, bd)

#ifndef SW_DYNAMICS

   !To get coarse grid pkz, etc right after a two-way update so
   !that it is consistent across a restart:
   !(should only be called after doing such an update)

    !! CLEANUP: move to twoway_nest_update??
   call p_var(npz, is, ie, js, je, ptop, ptop_min,  &
        delp, delz, &
        pt, ps, &
        pe, peln,   &
        pk,   pkz, kappa, &
        q, ng, flagstruct%ncnst,  gridstruct%area_64, 0.,  &
        .false.,  .false., & !mountain argument not used
        flagstruct%moist_phys,  flagstruct%hydrostatic, &
        flagstruct%nwat, domain, flagstruct%adiabatic, .false.)

#endif

      if (flagstruct%range_warn) then
         call range_check('TA update', pt, is, ie, js, je, ng, npz, gridstruct%agrid, 130., 350., bad_range, Time)
         call range_check('UA update', ua, is, ie, js, je, ng, npz, gridstruct%agrid, -220., 250., bad_range, Time)
         call range_check('VA update', va, is, ie, js, je, ng, npz, gridstruct%agrid, -220., 220., bad_range, Time)
         if (.not. flagstruct%hydrostatic) then
            call range_check('W update', w, is, ie, js, je, ng, npz, gridstruct%agrid, -50., 100., bad_range, Time)
         endif
      endif



 end subroutine after_twoway_nest_update

!>@brief The subroutine 'update_remap_tqw' remaps (interpolated) nested-grid data 
!! to the coarse-grid's vertical coordinate.
 !This does not yet do anything for the tracers

 subroutine update_remap_tqw( npz, ak_dst,  bk_dst,  ps_dst, t_dst, q_dst, w_dst, &
                      hydrostatic, &
                      kmd, ps_src, ak_src, bk_src, t_src, w_src, &
                      zvir, ptop, nq, kord_tm, kord_tr, kord_wz, &
                      is, ie, js, je, isd, ied, jsd, jed, do_q, &
                      istart, iend, jstart, jend, blend_wt)
  integer, intent(in):: npz, kmd, nq, kord_tm, kord_tr, kord_wz
  real,    intent(in):: zvir, ptop
  real,    intent(in):: ak_src(kmd+1), bk_src(kmd+1)
  real,    intent(in):: ak_dst(npz+1), bk_dst(npz+1), blend_wt(npz)
  real,    intent(in), dimension(isd:ied,jsd:jed):: ps_src
  real,    intent(in), dimension(isd:ied,jsd:jed):: ps_dst
  real,    intent(inout), dimension(isd:ied,jsd:jed,npz):: t_dst, w_dst
  real,    intent(inout), dimension(isd:ied,jsd:jed,npz,nq):: q_dst
  real,    intent(in), dimension(isd:ied,jsd:jed,kmd):: t_src, w_src
  integer,  intent(in) ::  is, ie, js, je, isd, ied, jsd, jed, istart, iend, jstart, jend
  logical,   intent(in) :: hydrostatic, do_q
! local:
  real, dimension(is:ie,kmd):: tp, qp
  real, dimension(is:ie,kmd+1):: pe0, pn0
  real, dimension(is:ie,npz):: qn1
  real, dimension(is:ie,npz+1):: pe1, pn1
  integer i,j,k,iq
  real :: wt1, wt2

  if (do_q) call mpp_error(FATAL, ' update_remap_tqw: q remapping not yet supported') 

  !This line to check if the update region is correctly defined or not is
  ! IMPORTANT. Sometimes one or the other pair of limits will give a
  ! non-empty loop, even though no data was transferred! This is why
  ! I was having so much trouble getting the remap-update to work --- lmh 11jul17
  if (istart > iend .or. jstart > jend) return

!$OMP parallel do default(none) shared(js,je,kmd,is,ie,ak_dst,bk_dst,ps_dst,q_dst,npz,ptop,do_q,&
!$OMP          t_dst,w_dst,t_src,w_src,ak_src,bk_src,ps_src,nq,hydrostatic,kord_tm,kord_tr,kord_wz,istart,iend,jstart,jend,blend_wt) &
!$OMP          private(pe0,pn0,pe1,pn1,qp,tp,qn1,wt1,wt2)
  do 5000 j=jstart,jend

     do k=1,kmd+1
        do i=istart,iend
           pe0(i,k) = ak_src(k) + bk_src(k)*ps_src(i,j)
           pn0(i,k) = log(pe0(i,k))
       enddo
     enddo 
     do k=1,npz+1
        do i=istart,iend
           pe1(i,k) = ak_dst(k) + bk_dst(k)*ps_dst(i,j)
           pn1(i,k) = log(pe1(i,k))
       enddo
     enddo 
     if (do_q) then
        do iq=1,nq
        do k=1,kmd
        do i=istart,iend
           qp(i,k) = q_dst(i,j,k,iq)
        enddo
        enddo
        call mappm(kmd, pe0, qp, npz, pe1,  qn1, is,ie, 0, kord_tr, ptop) !not sure about indices
        do k=1,npz
           do i=istart,iend
              q_dst(i,j,k,iq) = qn1(i,k)
           enddo
        enddo
        enddo
     endif

     do k=1,kmd
        do i=istart,iend
           tp(i,k) = t_src(i,j,k)
        enddo
     enddo
     !Remap T using logp
     call mappm(kmd, pn0(istart:iend,:), tp(istart:iend,:), npz, pn1(istart:iend,:), qn1(istart:iend,:), istart,iend, 1, abs(kord_tm), ptop)
     
     do k=1,npz
        wt1 = blend_wt(k)
        wt2 = 1. - wt1
        do i=istart,iend
           t_dst(i,j,k) = qn1(i,k)*wt1 + t_dst(i,j,k)*wt2
        enddo
     enddo

     if (.not. hydrostatic) then
        do k=1,kmd
           do i=istart,iend
              tp(i,k) = w_src(i,j,k)
           enddo
        enddo
        !Remap w using p
        !Using iv == -1 instead of -2
        call mappm(kmd, pe0(istart:iend,:), tp(istart:iend,:), npz, pe1(istart:iend,:), qn1(istart:iend,:), istart,iend, -1, kord_wz, ptop)

        do k=1,npz
           wt1 = blend_wt(k)
           wt2 = 1. - wt1
           do i=istart,iend
              w_dst(i,j,k) = qn1(i,k)*wt1 + w_dst(i,j,k)*wt2
           enddo
        enddo
     endif

5000 continue

 end subroutine update_remap_tqw

 !remap_uv as-is remaps only a-grid velocities. A new routine has been written to handle staggered grids.
 subroutine update_remap_uv(npz, ak_dst, bk_dst, ps_dst, u_dst, v_dst, &
                            kmd, ak_src, bk_src, ps_src, u_src, v_src, &
                            kord_mt, &
                            is, ie, js, je, isd, ied, jsd, jed, ptop, &
                            istart, iend, jstart, jend, blend_wt)
  integer, intent(in):: npz
  real,    intent(in):: ak_dst(npz+1), bk_dst(npz+1), blend_wt(npz)
  real,    intent(in):: ps_dst(isd:ied,jsd:jed)
  real,    intent(inout), dimension(isd:ied,jsd:jed+1,npz):: u_dst
  real,    intent(inout), dimension(isd:ied+1,jsd:jed,npz):: v_dst
  integer, intent(in):: kmd
  real,    intent(in):: ak_src(kmd+1), bk_src(kmd+1)
  real,    intent(in):: ps_src(isd:ied,jsd:jed)
  real,    intent(inout), dimension(isd:ied,jsd:jed+1,kmd):: u_src
  real,    intent(inout), dimension(isd:ied+1,jsd:jed,kmd):: v_src
!
  integer, intent(in):: kord_mt
  real,    intent(IN) :: ptop
  integer,  intent(in) ::  is, ie, js, je, isd, ied, jsd, jed
  integer,  intent(IN) :: istart, iend, jstart, jend
!
! local:
  real, dimension(is:ie+1,kmd+1):: pe0
  real, dimension(is:ie+1,npz+1):: pe1
  real, dimension(is:ie+1,kmd):: qt
  real, dimension(is:ie+1,npz):: qn1
  integer i,j,k
  real :: wt1, wt2

  !This line to check if the update region is correctly defined or not is
  ! IMPORTANT. Sometimes one or the other pair of limits will give a
  ! non-empty loop, even though no data was transferred!
  if (istart > iend .or. jstart > jend) return

!------
! map u
!------
!$OMP parallel do default(none) shared(js,je,kmd,is,ie,ak_dst,bk_dst,ps_dst,u_dst,v_dst,npz,ak_src,bk_src,ps_src,u_src,v_src,ptop,kord_mt,istart,iend,jstart,jend,blend_wt) &
!$OMP          private(pe0,pe1,qt,qn1,wt1,wt2)
  do j=jstart,jend+1
!------
! Data
!------
     do k=1,kmd+1
       do i=istart,iend
          pe0(i,k) = ak_src(k) + bk_src(k)*0.5*(ps_src(i,j)+ps_src(i,j-1))
       enddo
     enddo
!------
! Model
!------
     do k=1,npz+1
        do i=istart,iend
          pe1(i,k) = ak_dst(k) + bk_dst(k)*0.5*(ps_dst(i,j)+ps_dst(i,j-1))
       enddo
     enddo
!------
!Do map
!------
     qt = 0.
      do k=1,kmd
         do i=istart,iend
            qt(i,k) = u_src(i,j,k)
         enddo
      enddo
      qn1 = 0. 
      call mappm(kmd, pe0(istart:iend,:), qt(istart:iend,:), npz, pe1(istart:iend,:), qn1(istart:iend,:), istart,iend, -1, kord_mt, ptop)
      do k=1,npz
         wt1 = blend_wt(k)
         wt2 = 1. - wt1
         do i=istart,iend
            u_dst(i,j,k) = qn1(i,k)*wt1 + u_dst(i,j,k)*wt2
         enddo
      enddo

   end do

!------
! map v
!------
!$OMP parallel do default(none) shared(js,je,kmd,is,ie,ak_dst,bk_dst,ps_dst,u_dst,v_dst,ak_src,bk_src,ps_src,npz,u_src,v_src,ptop,istart,iend,jstart,jend,blend_wt) &
!$OMP          private(pe0,pe1,qt,qn1,wt1,wt2)
   do j=jstart,jend
!------
! Data
!------
     do k=1,kmd+1
        do i=istart,iend+1
          pe0(i,k) = ak_src(k) + bk_src(k)*0.5*(ps_src(i,j)+ps_src(i-1,j))
       enddo
     enddo
!------
! Model
!------
     do k=1,npz+1
        do i=istart,iend+1
          pe1(i,k) = ak_dst(k) + bk_dst(k)*0.5*(ps_dst(i,j)+ps_dst(i-1,j))
       enddo
     enddo
!------
!Do map
!------
     qt = 0.
      do k=1,kmd
         do i=istart,iend+1
            qt(i,k) = v_src(i,j,k)
         enddo
      enddo
      qn1 = 0.
      call mappm(kmd, pe0(istart:iend+1,:), qt(istart:iend+1,:), npz, pe1(istart:iend+1,:), qn1(istart:iend+1,:), istart,iend+1, -1, 8, ptop)
      do k=1,npz
         wt1 = blend_wt(k)
         wt2 = 1. - wt1
         do i=istart,iend+1
            v_dst(i,j,k) = qn1(i,k)*wt1 + v_dst(i,j,k)*wt2  !Does this kill OMP???
         enddo
      enddo
   end do

 end subroutine update_remap_uv



end module fv_nesting_mod
