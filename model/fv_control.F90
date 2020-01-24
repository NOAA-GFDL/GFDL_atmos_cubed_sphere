
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

!>@brief The module 'FV3_control' is for initialization and termination
!! of the model, and controls namelist parameters in FV3.
!----------------
! FV control panel
!----------------

module fv_control_mod
! Modules Included:
! <table>
! <tr>
!     <th>Module Name</th>
!     <th>Functions Included</th>
!   </tr>
! <table>
!   <tr>
!     <td>constants_mod</td>
!     <td>pi=>pi_8, kappa, radius, grav, rdgas</td>
!   </tr>
!   <tr>
!     <td>field_manager_mod</td>
!     <td>MODEL_ATMOS</td>
!   </tr>
!   <tr>
!     <td>fms_mod</td>
!     <td>write_version_number, open_namelist_file,
!         check_nml_error, close_file, file_exist</td>
!   </tr>
!   <tr>
!     <td>fv_arrays_mod</td>
!     <td>fv_atmos_type, allocate_fv_atmos_type, deallocate_fv_atmos_type,
!          R_GRID</td>
!   </tr>
!   <tr>
!     <td>fv_diagnostics_mod</td>
!     <td>fv_diag_init_gn</td>
!   </tr>
!   <tr>
!     <td>fv_eta_mod</td>
!     <td>set_eta</td>
!   </tr>
!   <tr>
!     <td>fv_grid_tools_mod</td>
!     <td>init_grid</td>
!   </tr>
!   <tr>
!     <td>fv_grid_utils_mod</td>
!     <td>grid_utils_init, grid_utils_end, ptop_min</td>
!   </tr>
!   <tr>
!     <td>fv_mp_mod</td>
!     <td>mp_start, mp_assign_gid, domain_decomp,ng, switch_current_Atm,
!         broadcast_domains, mp_barrier, is_master, setup_master </td>
!   </tr>
!   <tr>
!     <td>fv_io_mod</td>
!     <td>fv_io_exit</td>
!   </tr>
!   <tr>
!     <td>fv_restart_mod</td>
!     <td>fv_restart_init, fv_restart_end</td>
!   </tr>
!   <tr>
!     <td>fv_timing_mod</td>
!     <td>timing_on, timing_off, timing_init, timing_prt</td>
!   </tr>
!   <tr>
!     <td>mpp_mod</td>
!     <td>mpp_send, mpp_sync, mpp_transmit, mpp_set_current_pelist, mpp_declare_pelist, 
!         mpp_root_pe, mpp_recv, mpp_sync_self, mpp_broadcast, read_input_nml,
!         FATAL, mpp_error, mpp_pe, stdlog, mpp_npes, mpp_get_current_pelist, 
!         input_nml_file, get_unit, WARNING, read_ascii_file, INPUT_STR_LENGTH</td>
!   </tr>
!   <tr>
!     <td>mpp_domains_mod</td>
!     <td>mpp_get_data_domain, mpp_get_compute_domain, domain2D, mpp_define_nest_domains, 
!        nest_domain_type, mpp_get_global_domain, mpp_get_C2F_index, mpp_get_F2C_index,
!        mpp_broadcast_domain, CENTER, CORNER, NORTH, EAST, WEST, SOUTH</td>
!   </tr>
!   <tr>
!     <td>mpp_parameter_mod</td>
!     <td>AGRID_PARAM=>AGRID</td>
!   </tr>
!   <tr>
!     <td>test_cases_mod</td>
!     <td>test_case, bubble_do, alpha, nsolitons, soliton_Umax, soliton_size</td>
!   </tr>
!   <tr>
!     <td>tracer_manager_mod</td>
!     <td>tm_get_number_tracers => get_number_tracers,tm_get_tracer_index => get_tracer_index,   
!         tm_get_tracer_indices => get_tracer_indices, tm_set_tracer_profile => set_tracer_profile, 
!         tm_get_tracer_names => get_tracer_names,tm_check_if_prognostic=> check_if_prognostic,
!         tm_register_tracers => register_tracers</td>
!   </tr>
! </table>

   use constants_mod,       only: pi=>pi_8, kappa, radius, grav, rdgas
   use field_manager_mod,   only: MODEL_ATMOS
   use fms_mod,             only: write_version_number, open_namelist_file, &
                                  check_nml_error, close_file, file_exist
   use mpp_mod,             only: FATAL, mpp_error, mpp_pe, stdlog, &
                                  mpp_npes, mpp_get_current_pelist, &
                                  input_nml_file, get_unit, WARNING, &
                                  read_ascii_file, INPUT_STR_LENGTH
   use mpp_domains_mod,     only: mpp_get_data_domain, mpp_get_compute_domain
   use tracer_manager_mod,  only: tm_get_number_tracers => get_number_tracers, &
                                  tm_get_tracer_index   => get_tracer_index,   &
                                  tm_get_tracer_indices => get_tracer_indices, &
                                  tm_set_tracer_profile => set_tracer_profile, &
                                  tm_get_tracer_names   => get_tracer_names,   &
                                  tm_check_if_prognostic=> check_if_prognostic,&
                                  tm_register_tracers   => register_tracers

   use fv_io_mod,           only: fv_io_exit
   use fv_restart_mod,      only: fv_restart_init, fv_restart_end
   use fv_arrays_mod,       only: fv_atmos_type, allocate_fv_atmos_type, deallocate_fv_atmos_type, &
                                  R_GRID
   use fv_grid_utils_mod,   only: grid_utils_init, grid_utils_end, ptop_min
   use fv_eta_mod,          only: set_eta
   use fv_grid_tools_mod,   only: init_grid
   use fv_mp_mod,           only: mp_start, mp_assign_gid, domain_decomp
   use fv_mp_mod,           only: ng, switch_current_Atm
   use fv_mp_mod,           only: broadcast_domains, mp_barrier, is_master, setup_master
!!! CLEANUP: should be replaced by a getter function?
   use test_cases_mod,      only: test_case, bubble_do, alpha, nsolitons, soliton_Umax, soliton_size
   use fv_timing_mod,       only: timing_on, timing_off, timing_init, timing_prt
   use mpp_domains_mod,     only: domain2D
   use mpp_domains_mod,     only: mpp_define_nest_domains, nest_domain_type, mpp_get_global_domain
   use mpp_domains_mod,     only: mpp_get_C2F_index, mpp_get_F2C_index, mpp_broadcast_domain
   use mpp_domains_mod,     only: CENTER, CORNER, NORTH, EAST, WEST, SOUTH
   use mpp_mod,             only: mpp_send, mpp_sync, mpp_transmit, mpp_set_current_pelist, mpp_declare_pelist, mpp_root_pe, mpp_recv, mpp_sync_self, mpp_broadcast, read_input_nml
   use fv_diagnostics_mod,  only: fv_diag_init_gn

#ifdef MULTI_GASES
   use constants_mod,       only: rvgas, cp_air
   use multi_gases_mod,     only: multi_gases_init, &
                                  rilist => ri,     &
                                  cpilist => cpi
#endif

   implicit none
   private

!-----------------------------------------------------------------------
! Grid descriptor file setup
!-----------------------------------------------------------------------
!------------------------------------------
! Model Domain parameters
! See fv_arrays.F90 for descriptions
!------------------------------------------
!CLEANUP module pointers
   character(len=80) , pointer :: grid_name
   character(len=120), pointer :: grid_file
   integer, pointer :: grid_type
   integer , pointer :: hord_mt 
   integer , pointer :: kord_mt 
   integer , pointer :: kord_wz 
   integer , pointer :: hord_vt 
   integer , pointer :: hord_tm 
   integer , pointer :: hord_dp 
   integer , pointer :: kord_tm 
   integer , pointer :: hord_tr 
   integer , pointer :: kord_tr 
   real    , pointer :: scale_z 
   real    , pointer :: w_max 
   real    , pointer :: z_min 
   real    , pointer :: lim_fac

   integer , pointer :: nord
   integer , pointer :: nord_tr
   real    , pointer :: dddmp 
   real    , pointer :: d2_bg 
   real    , pointer :: d4_bg 
   real    , pointer :: vtdm4 
   real    , pointer :: trdm2 
   real    , pointer :: d2_bg_k1 
   real    , pointer :: d2_bg_k2 
   real    , pointer :: d2_divg_max_k1 
   real    , pointer :: d2_divg_max_k2 
   real    , pointer :: damp_k_k1 
   real    , pointer :: damp_k_k2 
   integer , pointer ::    n_zs_filter
   integer , pointer :: nord_zs_filter
   logical , pointer :: full_zs_filter

   logical , pointer :: RF_fast
   logical , pointer :: consv_am
   logical , pointer :: do_sat_adj
   logical , pointer :: do_f3d
   logical , pointer :: no_dycore 
   logical , pointer :: convert_ke 
   logical , pointer :: do_vort_damp 
   logical , pointer :: use_old_omega 
! PG off centering:
   real    , pointer :: beta  
   integer , pointer :: n_sponge 
   real    , pointer :: d_ext 
   integer , pointer :: nwat  
   logical , pointer :: warm_start 
   logical , pointer :: inline_q 
   real , pointer :: shift_fac   
   logical , pointer :: do_schmidt 
   real(kind=R_GRID) , pointer :: stretch_fac 
   real(kind=R_GRID) , pointer :: target_lat  
   real(kind=R_GRID) , pointer :: target_lon  

   logical , pointer :: reset_eta 
   real    , pointer :: p_fac
   real    , pointer :: a_imp
   integer , pointer :: n_split 

   real    , pointer :: fac_n_spl
   real    , pointer :: fhouri
                             ! Default 
   integer , pointer :: m_split 
   integer , pointer :: k_split 
   logical , pointer :: use_logp

   integer , pointer :: q_split 
   integer , pointer :: print_freq 
   logical , pointer :: write_3d_diags

   integer , pointer :: npx           
   integer , pointer :: npy           
   integer , pointer :: npz
   integer , pointer :: npz_rst 
                                      
   integer , pointer :: ncnst 
   integer , pointer :: pnats 
   integer , pointer :: dnats 
   integer , pointer :: ntiles        
   integer , pointer :: nf_omega  
   integer , pointer :: fv_sg_adj 
                                      
   integer , pointer :: na_init 
   logical , pointer :: nudge_dz
   real    , pointer :: p_ref 
   real    , pointer :: dry_mass 
   integer , pointer :: nt_prog 
   integer , pointer :: nt_phys 
   real    , pointer :: tau_h2o 

   real    , pointer :: delt_max
   real    , pointer :: d_con 
   real    , pointer :: ke_bg
   real    , pointer :: consv_te 
   real    , pointer :: tau 
   real    , pointer :: rf_cutoff
   logical , pointer :: filter_phys 
   logical , pointer :: dwind_2d 
   logical , pointer :: breed_vortex_inline 
   logical , pointer :: range_warn 
   logical , pointer :: fill 
   logical , pointer :: fill_dp 
   logical , pointer :: fill_wz 
   logical , pointer :: check_negative
   logical , pointer :: non_ortho 
   logical , pointer :: adiabatic 
   logical , pointer :: moist_phys 
   logical , pointer :: do_Held_Suarez 
   logical , pointer :: do_reed_physics
   logical , pointer :: reed_cond_only
   logical , pointer :: reproduce_sum 
   logical , pointer :: adjust_dry_mass 
   logical , pointer :: fv_debug  
   logical , pointer :: srf_init  
   logical , pointer :: mountain  
   logical , pointer :: remap_t  
   logical , pointer :: z_tracer 

   logical , pointer :: old_divg_damp 
   logical , pointer :: fv_land 
   logical , pointer :: nudge 
   logical , pointer :: nudge_ic
   logical , pointer :: ncep_ic 
   logical , pointer :: nggps_ic 
   logical , pointer :: ecmwf_ic 
   logical , pointer :: gfs_phil
   logical , pointer :: agrid_vel_rst
   logical , pointer :: use_new_ncep 
   logical , pointer :: use_ncep_phy 
   logical , pointer :: fv_diag_ic 
   logical , pointer :: external_ic 
   logical , pointer :: external_eta
   logical , pointer :: read_increment
   character(len=128) , pointer :: res_latlon_dynamics
   character(len=128) , pointer :: res_latlon_tracers 
   logical , pointer :: hydrostatic 
   logical , pointer :: phys_hydrostatic
   logical , pointer :: use_hydro_pressure
   logical , pointer :: do_uni_zfull !miz
   logical , pointer :: adj_mass_vmr ! f1p
   logical , pointer :: hybrid_z    
   logical , pointer :: Make_NH     
   logical , pointer :: make_hybrid_z  
   logical , pointer :: nudge_qv
   real,     pointer :: add_noise

   integer , pointer :: a2b_ord 
   integer , pointer :: c2l_ord 

   integer, pointer :: ndims

  real(kind=R_GRID), pointer :: dx_const
  real(kind=R_GRID), pointer :: dy_const
  real(kind=R_GRID), pointer :: deglon_start, deglon_stop, &  ! boundaries of latlon patch
          deglat_start, deglat_stop
  real(kind=R_GRID), pointer :: deglat

  logical, pointer :: nested, twowaynest
  logical, pointer :: regional
  integer, pointer :: bc_update_interval 
  integer, pointer :: parent_tile, refinement, nestbctype, nestupdate, nsponge, ioffset, joffset
  real, pointer :: s_weight, update_blend

  integer, pointer :: layout(:), io_layout(:)

   integer :: ntilesMe                ! Number of tiles on this process =1 for now

#ifdef OVERLOAD_R4
   real    :: too_big  = 1.E8
#else
   real    :: too_big  = 1.E35
#endif
   public :: fv_init, fv_end

   integer, public :: ngrids = 1
   integer, public, allocatable :: pelist_all(:)
   integer :: commID, max_refinement_of_global = 1.
   integer :: gid

   real :: umax = 350.           !< max wave speed for grid_type>3
   integer :: parent_grid_num = -1

   integer :: halo_update_type = 1 !< 1 for two-interfaces non-block
                                   !< 2 for block
                                   !< 3 for four-interfaces non-block



! version number of this module
! Include variable "version" to be written to log file.
#include<file_version.h>

 contains

!-------------------------------------------------------------------------------
!>@brief The subroutine 'fv_init' initializes FV3.
!>@details It allocates memory, sets up MPI and processor lists,
!! sets up the grid, and controls FV3 namelist parameters.   
 subroutine fv_init(Atm, dt_atmos, grids_on_this_pe, p_split)

   type(fv_atmos_type), allocatable, intent(inout), target :: Atm(:)
   real,                intent(in)    :: dt_atmos
   logical, allocatable, intent(INOUT) :: grids_on_this_pe(:)
   integer, intent(INOUT) :: p_split

   integer :: i, j, k, n, p
   real :: sdt

! tracers
   integer :: num_family          !< output of register_tracers

   integer :: isc_p, iec_p, jsc_p, jec_p, isg, ieg, jsg, jeg, upoff, jind
   integer :: ic, jc

   gid = mpp_pe()
   call init_nesting(Atm, grids_on_this_pe, p_split)

   !This call is needed to set up the pointers for fv_current_grid, even for a single-grid run
   !call switch_current_Atm(Atm(1), .false.)
   call setup_pointers(Atm(1))

! Start up MPI

   !call mp_assign_gid

    ! Initialize timing routines
      call timing_init
      call timing_on('TOTAL')

    ! Setup the run from namelist 
      ntilesMe = size(Atm(:)) !Full number of Atm arrays; one less than number of grids, if multiple grids

      call run_setup(Atm,dt_atmos, grids_on_this_pe, p_split)   ! initializes domain_decomp

      do n=1,ntilesMe

         !In a single-grid run this will still be needed to correctly set the domain
         call switch_current_Atm(Atm(n))
         call setup_pointers(Atm(n))
         
         target_lon = target_lon * pi/180.
         target_lat = target_lat * pi/180.

!--------------------------------------------------
! override number of tracers by reading field_table
!--------------------------------------------------

         !not sure if this works with multiple grids
         call tm_register_tracers (MODEL_ATMOS, ncnst, nt_prog, pnats, num_family)
         if(is_master()) then
            write(*,*) 'ncnst=', ncnst,' num_prog=',nt_prog,' pnats=',pnats,' dnats=',dnats,' num_family=',num_family         
            print*, ''
         endif

         if (grids_on_this_pe(n)) then
            call allocate_fv_atmos_type(Atm(n), Atm(n)%bd%isd, Atm(n)%bd%ied, Atm(n)%bd%jsd, Atm(n)%bd%jed, &
                 Atm(n)%bd%isc, Atm(n)%bd%iec, Atm(n)%bd%jsc, Atm(n)%bd%jec, &
                 npx, npy, npz, ndims, ncnst, ncnst-pnats, ng, .false., grids_on_this_pe(n), ngrids)

            if (grids_on_this_pe(n)) then

               call switch_current_Atm(Atm(n))
               call setup_pointers(Atm(n))

               if ( (Atm(n)%bd%iec-Atm(n)%bd%isc+1).lt.4 .or. (Atm(n)%bd%jec-Atm(n)%bd%jsc+1).lt.4 ) then
                  if (is_master()) write(*,'(6I6)') Atm(n)%bd%isc, Atm(n)%bd%iec, Atm(n)%bd%jsc, Atm(n)%bd%jec, n
                  call mpp_error(FATAL,'Domain Decomposition:  Cubed Sphere compute domain has a &
                       &minium requirement of 4 points in X and Y, respectively')
               end if

            endif

            !!CLEANUP: Convenience pointers
            Atm(n)%gridstruct%nested      => Atm(n)%neststruct%nested
            Atm(n)%gridstruct%grid_type   => Atm(n)%flagstruct%grid_type
            Atm(n)%flagstruct%grid_number => Atm(n)%grid_number
            Atm(n)%gridstruct%regional  => Atm(n)%flagstruct%regional

            call init_grid(Atm(n), grid_name, grid_file, npx, npy, npz, ndims, ntiles, ng)

            ! Initialize the SW (2D) part of the model
            !!!CLEANUP: this call could definitely use some cleaning up
            call grid_utils_init(Atm(n), npx, npy, npz, non_ortho, grid_type, c2l_ord)

            !!!CLEANUP: Are these correctly writing out on all pes?
            if ( is_master() ) then
               sdt =  dt_atmos/real(n_split*k_split*abs(p_split))
               write(*,*) ' '
               write(*,*) 'Divergence damping Coefficients'
               write(*,*) 'For small dt=', sdt
               write(*,*) 'External mode del-2 (m**2/s)=',  d_ext*Atm(n)%gridstruct%da_min_c/sdt
               write(*,*) 'Internal mode del-2 SMAG dimensionless coeff=',  dddmp
               write(*,*) 'Internal mode del-2 background diff=', d2_bg*Atm(n)%gridstruct%da_min_c/sdt

               if (nord==1) then
                   write(*,*) 'Internal mode del-4 background diff=', d4_bg
                   write(*,*) 'Vorticity del-4 (m**4/s)=', (vtdm4*Atm(n)%gridstruct%da_min)**2/sdt*1.E-6
               endif
               if (nord==2) write(*,*) 'Internal mode del-6 background diff=', d4_bg
               if (nord==3) write(*,*) 'Internal mode del-8 background diff=', d4_bg
               write(*,*) 'tracer del-2 diff=', trdm2

               write(*,*) 'Vorticity del-4 (m**4/s)=', (vtdm4*Atm(n)%gridstruct%da_min)**2/sdt*1.E-6
               write(*,*) 'beta=', beta
               write(*,*) ' '
            endif


            Atm(n)%ts   = 300.
            Atm(n)%phis = too_big
            ! The following statements are to prevent the phatom corner regions from
            ! growing instability
            Atm(n)%u  = 0.
            Atm(n)%v  = 0.
            Atm(n)%ua = too_big
            Atm(n)%va = too_big

         else !this grid is NOT defined on this pe

            !Allocate dummy arrays
            call allocate_fv_atmos_type(Atm(n),  Atm(n)%bd%isd, Atm(n)%bd%ied, Atm(n)%bd%jsd, Atm(n)%bd%jed, &
                 Atm(n)%bd%isc, Atm(n)%bd%iec, Atm(n)%bd%jsc, Atm(n)%bd%jec, &
                 npx, npy, npz, ndims, ncnst, ncnst-pnats, ng, .true., .false., ngrids)

            !Need to SEND grid_global to any child grids; this is received in setup_aligned_nest in fv_grid_tools
            if (Atm(n)%neststruct%nested) then

               call mpp_get_global_domain( Atm(n)%parent_grid%domain, &
                    isg, ieg, jsg, jeg)

               !FIXME: Should replace this by generating the global grid (or at least one face thereof) on the 
               ! nested PEs instead of sending it around.
               if (gid == Atm(n)%parent_grid%pelist(1)) then
                     call mpp_send(Atm(n)%parent_grid%grid_global(isg-ng:ieg+1+ng,jsg-ng:jeg+1+ng,1:2,parent_tile), &
                          size(Atm(n)%parent_grid%grid_global(isg-ng:ieg+1+ng,jsg-ng:jeg+1+ng,1:2,parent_tile)), &
                       Atm(n)%pelist(1)) !send to p_ind in setup_aligned_nest
                  call mpp_sync_self()
               endif

               if (Atm(n)%neststruct%twowaynest) then

                  !This in reality should be very simple. With the
                  ! restriction that only the compute domain data is
                  ! sent from the coarse grid, we can compute
                  ! exactly which coarse grid cells should use
                  ! which nested-grid data. We then don't need to send around p_ind.

                  Atm(n)%neststruct%ind_update_h = -99999

                  if (Atm(n)%parent_grid%tile == Atm(n)%neststruct%parent_tile) then

                     isc_p = Atm(n)%parent_grid%bd%isc
                     iec_p = Atm(n)%parent_grid%bd%iec
                     jsc_p = Atm(n)%parent_grid%bd%jsc
                     jec_p = Atm(n)%parent_grid%bd%jec
                     upoff = Atm(n)%neststruct%upoff

                     Atm(n)%neststruct%jsu = jsc_p
                     Atm(n)%neststruct%jeu = jsc_p-1
                     do j=jsc_p,jec_p+1
                        if (j < joffset+upoff) then
                           do i=isc_p,iec_p+1
                              Atm(n)%neststruct%ind_update_h(i,j,2) = -9999
                           enddo
                           Atm(n)%neststruct%jsu = Atm(n)%neststruct%jsu + 1
                        elseif (j > joffset + (npy-1)/refinement - upoff) then
                           do i=isc_p,iec_p+1
                              Atm(n)%neststruct%ind_update_h(i,j,2) = -9999
                           enddo
                        else
                           jind = (j - joffset)*refinement + 1
                           do i=isc_p,iec_p+1
                              Atm(n)%neststruct%ind_update_h(i,j,2) = jind
                           enddo
                           if ( (j < joffset + (npy-1)/refinement - upoff) .and. j <= jec_p)  Atm(n)%neststruct%jeu = j
                        endif
                        !write(mpp_pe()+4000,*) j, joffset, upoff, Atm(n)%neststruct%ind_update_h(isc_p,j,2)
                     enddo

                     Atm(n)%neststruct%isu = isc_p
                     Atm(n)%neststruct%ieu = isc_p-1
                     do i=isc_p,iec_p+1
                        if (i < ioffset+upoff) then
                           Atm(n)%neststruct%ind_update_h(i,:,1) = -9999
                           Atm(n)%neststruct%isu = Atm(n)%neststruct%isu + 1
                        elseif (i > ioffset + (npx-1)/refinement - upoff) then
                           Atm(n)%neststruct%ind_update_h(i,:,1) = -9999
                        else
                           Atm(n)%neststruct%ind_update_h(i,:,1) = (i-ioffset)*refinement + 1
                           if ( (i < ioffset + (npx-1)/refinement - upoff) .and. i <= iec_p) Atm(n)%neststruct%ieu = i
                           end if
                        !write(mpp_pe()+5000,*) i, ioffset, upoff, Atm(n)%neststruct%ind_update_h(i,jsc_p,1)
                     enddo

                  end if


               end if

            endif
         endif
      end do
      
    ! Initialize restart functions
      call fv_restart_init()

!     if ( reset_eta ) then
!         do n=1, ntilesMe
!            call set_eta(npz, Atm(n)%ks, ptop, Atm(n)%ak, Atm(n)%bk)
!         enddo
!         if(is_master()) write(*,*) "Hybrid sigma-p coordinate has been reset"
!     endif

      if (ntilesMe > 1) call switch_current_Atm(Atm(1))
      if (ntilesMe > 1) call setup_pointers(Atm(1))

 end subroutine fv_init
!-------------------------------------------------------------------------------

!>@brief The subroutine 'fv_end' terminates FV3, deallocates memory, 
!! saves restart files, and stops I/O.
 subroutine fv_end(Atm, grids_on_this_pe)

    type(fv_atmos_type), intent(inout) :: Atm(:)
    logical, intent(INOUT) :: grids_on_this_pe(:)

    integer :: n

    call timing_off('TOTAL')
    call timing_prt( gid )

    call fv_restart_end(Atm, grids_on_this_pe)
    call fv_io_exit()

  ! Free temporary memory from sw_core routines

  ! Deallocate
    call grid_utils_end

    do n = 1, ntilesMe
       call deallocate_fv_atmos_type(Atm(n))
    end do


 end subroutine fv_end
!-------------------------------------------------------------------------------

!>@brief The subroutine 'run_setup' initializes the run from a namelist.
 subroutine run_setup(Atm, dt_atmos, grids_on_this_pe, p_split)
   type(fv_atmos_type), intent(inout), target :: Atm(:)
   real, intent(in)                   :: dt_atmos
   logical, intent(INOUT) :: grids_on_this_pe(:)
   integer, intent(INOUT) :: p_split
   !--- local variables ---
   character(len=80) :: tracerName, errString
   character(len=32) :: nested_grid_filename
   integer :: ios, ierr, f_unit, unit
   logical :: exists

   real :: dim0 = 180.           !< base dimension
   real :: dt0  = 1800.          !< base time step
   real :: ns0  = 5.             !< base nsplit for base dimension 
                                 !< For cubed sphere 5 is better
   !real :: umax = 350.          ! max wave speed for grid_type>3 ! Now defined above
   real :: dimx, dl, dp, dxmin, dymin, d_fac

   integer :: n0split
   integer :: n, nn, i

   integer :: pe_counter

!  local version of these variables to allow PGI compiler to compile
   character(len=128) :: res_latlon_dynamics = ''
   character(len=128) :: res_latlon_tracers  = ''
   character(len=80)  :: grid_name = ''
   character(len=120) :: grid_file = ''
!> \defgroup Parameters_List
!! @{
!> ##A.1 Entries in fv\_core\_nml
!!
!>###A.1.1 Required options:
!!
!> \param[in] layout Integer(2): Processor layout on each tile. The number of PEs assigned to a domain must equal layout(1)*layout(2)*ntiles. Must be set. 
!!
!> \param[in] npx  Integer: Number of grid *corners* in the x-direction on one tile of the domain; so one more than the number of grid cells across a tile. On the cubed sphere this is *one* *more* *than* the number of cells across a cube face. Must be set.
!! 
!> \param[in] npy Integer: Number of grid *corners* in the y-direction on one tile of the domain. This value should be identical to npx on a cubed-sphere grid; doubly periodic or nested grids do not have this restriction. Must be set.
!! 
!> \param[in] npz Integer: Number of vertical levels. Each choice of npz comes with a pre-defined set of hybrid sigma pressure levels and model top (see fv\_eta.F90). Must be set.
!! 
!> \param[in]  ntiles Integer: Number of tiles on the domain. For the cubed sphere, this should be 6, one tile for each face of the cubed sphere; normally for most other domains  (including nested grids) this should be set to 1. Must be set. 
!!
!>###A.1.2 Initialization options:

!> \param[in] add\_noise] Real: amplitude of random thermal noise (in K) to add upon startup. Useful for perturbing initial conditions. -1 by default; disabled if 0 or negative. 
!!
!> \param[in] adjust\_dry\_mass Logical: whether to adjust the global dry-air mass to the value set by dry\_mass. This is only done in an initialization step, particularly when using an initial condition from an external dataset, interpolated from another resolution (either horizontal or vertical), or when changing the topography, so that the global mass of the atmosphere matches some estimate of observed value. False by default. It is recommended to only set this to True when initializing the model. 
!!
!> \param[in] breed\_vortex\_inline Logical: whether to bogus tropical cyclones into the model, which are specified by an external file. Options are set in fv\_nwp\_nudge\_nml. False by default. 
!!
!> \param[in] dry\_mass Real: if adjust\_dry\_mass is true, sets the global dry air mass, measured in the globally-averaged surface pressure (Pascals) by adding or removing mass from the lowest layer of the atmosphere as needed. 98290. (Pa) by default. 
!!
!> \param[in] external\_ic Logical: Whether to initialize the models state using the data in an externally specified file, given in res\_latlon\_dynamics. By default this file is assumed to be a legacy lat-lon FV core restart file; set either ncep\_ic or fv\_diag\_ic to true override this behavior..false. by default. Note that external\_ic = true will cause the model to re-initialize the dynamical fields from the input dataset regardless of whether warm\_start is set. 
!!
!> \param[in] full\_zs\_filter Logical: whether to apply the on-line topography filter during initialization. Only active if get\_nggps\_ic = .true. This is so topography filtering can be performed on the initial conditions output by the pre-processing tools, which currently do not support topography filtering for some configurations (such as the nested grid); this also allows the user to easily test changes to the topography filtering on the simulation. Note that for all other initialization methods (if external\_ic = .true.) the on-line topography filter will be applied automatically during the initialization of the topography. .false. by default. 
!!
!> \param[in] mountain Logical: takes topography into account when initializing the model. Set this to true to apply the terrain filter (if n\_zs\_filter = 2 or 4) upon startup; also set to True when cold starting so that the topography can be initialized. Only set this to false if you wish to cold-start without any topography; this value is ignored for the aquaplanet test\_case = 14. True by default. It is highly recommended to not alter this value unless you know what you are doing. 
!!
!> \param[in] na\_init Integer: Number of forward-backward dynamics steps used to initialize adiabatic solver. This is useful for spinning up the nonhydrostatic state from the hydrostatic GFS analyses. 0 by default. Recommended to set this to a non-zero value (1 or 2 is typically sufficient) when initializing from GFS or ECMWF analyses.
!!
!> \param[in] ncep\_ic Logical: If external\_ic =.true., this variable says whether the file in res\_latlon\_dynamics is an NCEP analysis or reanalysis file. This option zeros out all tracer fields except specific humidity..false. by default. 
!!
!> \param[in] nggps\_ic Logical: If external\_ic =.true., reads initial conditions from horizontally-interpolated output from chgres. False by default. Additional options are available through external\_ic\_nml.
!!
!> \param[in] ecmwf\_ic Logical: If external\_ic =.true., reads initial conditions from ECMWF analyses.  .false. by default. 
!!
!> \param[in] external\_eta Logical: If .true., reads the interface coefficients a<sub>k</sub> and b<sub>k</sub> from either the restart file (if restarting) or from the external initial condition file (if nggps\_ic or ecwmf\_ic are .true.). This overrides the hard-coded levels in fv\_eta. .false. by default. 
!!
!> \param[in] nord\_zs\_filter Integer: order of the topography filter applied to n\_zs\_filter. Set to 2 to get a second-order filter, or 4 to get a fourth-order filter; other values do no filtering. 0 by default. This should not be set to a non-zero value on multiple successive simulations; the filter is applied every time the model restarts. This option is useful for testing the terrain filter, and should not be used for regular runs.
!!
!> \param[in] npz\_rst Integer: If using a restart file with a different number of vertical levels, set npz\_rst to be the number of levels in your restart file. The model will then remap the restart file data to the vertical coordinates specified by npz. 0 by default; if 0 or equal to npz no remapping is done. 
!!
!> \param[in] nudge Logical: whether to use the nudging towards the state in some externally-supplied file (such as from reanalysis or another simulation). Further nudging options are set in fv\_nwp\_nudge\_nml. False by default. 
!!
!> \param[in] nudge\_dz Logical: during the adiabatic initialization (na\_init > 0), if set to .true.  delz is nudged back to the value specified in the initial conditions, instead of nudging the temperature back to the initial value. Nudging delz is simpler (faster), doesn't require consideration of the virtual temperature effect, and may be more stable. .false. by default.
!!
!> \param[in] nudge\_ic Logical: same as nudge, but works in adiabatic solo\_core simulations to nudge the field to a single external analysis file. False by default.
!!
!> \param[in] nudge\_qv Logical: during the adiabatic initialization (na\_init > 0), if set to .true., the water vapor profile is nudged to an analytic fit to the HALOE climatology. This is to improve the water vapor concentrations in GFS initial conditions, especially in the stratosphere, where values can be several times higher than observed. This nudging is unnecessary for other ICs, especially the ECMWF initial conditions. .false. by default. 
!!
!> \param[in] n\_zs\_filter Integer: number of times to apply a diffusive filter to the topography upon startup, if mountain is True and the model is not being cold-started. This is applied every time the model is warm-started, so if you want to smooth the topography make sure this is set to 0 after the first simulation. 0 by default. If initializing the model from cold-start the topography is already being filtered by an amount appropriate for the model resolution. 
!!
!> \param[in] res\_latlon\_dynamics character(len=128) If external\_ic =.true. gives the filename of the input IC file. INPUT/fv\_rst.res.nc by default. 
!!
!> \param[in] res\_latlon\_tracers character(len=128) If external\_ic =.true. and both ncep\_ic and fv\_diag\_ic are.false., this variable gives the filename of the initial conditions for the tracers, assumed to be a legacy lat-lon FV core restart file. INPUT/atmos\_tracers.res.nc by default. 
!!
!> \param[in] warm\_start] Logical; whether to start from restart files, instead of cold-starting the model. True by default; if this is set to true and restart files cannot be found the model will stop.
!!
!>###A1.3 I/O and diagnostic options:
!!
!> \param[in] agrid\_vel\_rst Logical: whether to write the unstaggered latitude-longitude winds (u<sub>a</sub> and v<sub>a</sub>) to the restart files. This is useful for data assimilation cycling systems which do not handle staggered winds. .false. by default.
!!
!> \param[in] check\_negative Logical: whether to print the most negative global value of microphysical tracers.
!!
!> \param[in] fv\_debug Logical: whether to turn on additional diagnostics in fv\_dynamics..false. by default. 
!!
!> \param[in] fv\_land Logical: whether to create terrain deviation and land fraction for output to mg\_drag restart files, for use in mg\_drag and in the land model..false. by default;.true. is recommended when, and only when, initializing the model, since the mg\_drag files created provide a much more accurate terrain representation for the mountain gravity wave drag parameterization and for the land surface roughness than either computes internally. This has no effect on the representation of the terrain in the dynamics. 
!!
!> \param[in] io\_layout Integer(2): Layout of output files on each tile. 1,1 by default, which combines all restart and history files on a tile into one file. For 0,0, every process writes out its own restart and history files. If not equal to  1,1, you will have to use mppnccombine to combine these output files prior to post-processing, or if you want to change the number of PEs. Both entries must divide the respective value in layout.
!!
!> \param[in] nf\_omega Integer: number of times to apply second-order smoothing to the diagnosed omega. When 0 the filter is disabled. 1 by default. 
!!
!> \param[in] print\_freq Integer: number of hours between print out of max/min and air/tracer mass diagnostics to standard output. 0 by default, which never prints out any output; set to -1 to see output after every dt\_atmos. Computing these diagnostics requires some computational overhead. 
!!
!> \param[in] range\_warn Logical: checks whether the values of the prognostic variables are within a reasonable range at the end of a dynamics time step, and prints a warning if not. False by default; adds computational overhead, so we only recommend using this when debugging. 
!!
!> \param[in] write\_3d\_diags Logical: whether to write out three-dimensional dynamical diagnostic fields (those defined in fv\_diagnostics.F90). This is useful for runs with multiple grids if you only want very large 3D diagnostics written out for (say) a nested grid, and not for the global grid. False by default. 
!!
!>###A.1.4 Options controlling tracers and interactions with physics
!!
!> \param[in] adiabatic Logical: whether to skip any physics. If true, the physics is not called at all and there is no virtual temperature effect. False by default; this option has no effect if not running solo\_core.
!!
!> \param[in] do\_Held\_Suarez Logical: whether to use Held-Suarez forcing. Requires adiabatic to be false. False by default; this option has no effect if not running solo\_core. 
!!
!> \param[in] do\_uni\_zfull Logical: whether to compute z\_full (the height of each model layer, as opposed to z\_half, the height of each model interface) as the midpoint of the layer, as is done for the nonhydrostatic solver, instead of the height of the location where <SPAN STYLE="text-decoration:overline">p</SPAN>  the mean pressure in the layer. This option is not available for fvGFS or the solo\_core. .false. by default.
!!
!> \param[in] dnats Integer: The number of tracers which are not to be advected by the dynamical core, but still passed into the dynamical core; the last dnats+pnats tracers in field\_table are not advected. 0 by default.
!!
!> \param[in] dwind\_2d Logical: whether to use a simpler \& faster algorithm for interpolating the A-grid (cell-centered) wind tendencies computed from the physics to the D-grid. Typically, the A-grid wind tendencies are first converted in 3D cartesian coordinates and then interpolated before converting back to 2D local coordinates. When this option enabled, a much simpler but less accurate 2D interpolation is used. False by default. 
!!
!> \param[in] fill Logical: Fills in negative tracer values by taking positive tracers from the cells above and below. This option is useful when the physical parameterizations produced negatives. False by default.
!!
!> \param[in] inline\_q Logical: whether to compute tracer transport in-line with the rest of the dynamics instead of sub-cycling, so that tracer transport is done at the same time and on the same time step as is `p` and potential temperature. False by default; if true, q\_split and z\_tracer are ignored. 
!!
!> \param[in] ncnst Integer: Number of tracer species advected by fv\_tracer in the dynamical core. Typically this is set automatically by reading in values from field\_table, but ncnst can be set to a smaller value so only the first ncnst tracers listed in field\_table are not advected. 0 by default, which will use the value from field\_table. 
!!
!> \param[in] nwat Integer: Number of water species to be included in condensate and water vapor loading. The masses of the first nwat tracer species will be added to the dry air mass, so that `p` is the mass of dry air, water vapor, and the included condensate species. The value used depends on the microphysics in the physics package you are using. For GFS physics with only a single condensate species, set to 2. For schemes with prognostic cloud water and cloud ice, such as GFDL AM2/AM3/AM4 Rotsteyn-Klein or Morrison-Gettlean microphysics, set to 3. For warm-rain (Kessler) microphysics set to 4 (with an inactive ice tracer), which only handles three species but uses 4 to avoid interference with the R-K physics. For schemes such as WSM5 or Ferrier that have prognostic rain and snow but not hail, set to 5 (not yet implemented). For six-category schemes that also have prognostic hail or graupel, such as the GFDL, Thompson, or WSM6 microphysics, set to  6. A value of 0 turns off condensate loading. 3 by default.
!!
!> \param[in] phys\_hydrostatic Logical: Option to enable hydrostatic application of heating from the physics in a nonhydrostatic simulation: heating is applied in hydrostatic balance, causing the entire atmospheric column to expand instantaneously. If false, heating from the physics is applied simply as a temperature tendency. True by default; ignored if hydrostatic =.true.
!!
!> \param[in] pnats Integer: The number of tracers not to advect by the dynamical core. Unlike dnats, these tracers are not seen by the dynamical core. The last pnats entries in field\_table are not advected. 0 by default.
!!
!> \param[in] tau\_h2o Real: time-scale (days) for simple methane chemistry to act as a source of water in the stratosphere. Can be useful if your stratosphere dries out too quickly; consider a value between 60 and 120 days if this is the case. 0. by default, which disables the methane chemistry. Values less than zero apply the chemistry above 100 mb; else applied above 30 mb. Requires adiabatic to be false. 
!!
!> \param[in] use\_hydro\_pressure Logical: whether to compute hydrostatic pressure for input to the physics. Currently only enabled for the fvGFS model. Ignored in hydrostatic simulations. False by default.
!!
!> \param[in] z\_tracer Logical: whether to transport sub-cycled tracers layer-by-layer, each with its own computed sub-cycling time step (if q\_split = 0). This may improve efficiency for very large numbers of tracers. False by default; currently not implemented.
!!
!>###A.1.5 Timestep options

!> \param[in] k\_split  Integer: number of vertical remappings per dt\_atmos (physics time step). 1 by default. 
!!
!> \param[in] n\_split  Integer: number of small dynamics (acoustic) time steps between vertical remapping. 0 by default, in which case the model produces a good first guess by examining the resolution, dt\_atmos, and k\_split. 
!!
!> \param[in] umax  Real: for the doubly-periodic grid (grid\_type = 4) an estimate of the maximum wave speed (m/s), used to determine the value of n\_split when n\_split = 0. 350 by default. 
!!
!> \param[in] q\_split]  Integer: number of time steps for sub-cycled tracer advection. 0 by default (recommended), in which case the model determines the number of time steps from the global maximum wind speed at each call to the tracer advection.
!!
!>###A.1.6 Grid options
!!
!> \param[in] deglat  Real: Latitude (in degrees) used to compute the uniform f-plane Coriolis parameter for doubly-periodic simulations (grid\_type = 4). 15. by default. 
!!
!> \param[in] do\_schmidt  Logical: Whether to enable grid stretching and rotation using stretch\_fac, target\_lat, and target\_lon..false. by default. 
!!
!> \param[in] dx\_const  Real: on a doubly-periodic grid (grid\_type = 4) specifies the (uniform) grid-cell-width in the x-direction, in meters. 1000 by default. 
!!
!> \param[in] dy\_const  Real: on a doubly-periodic grid (grid\_type = 4) specifies the (uniform) grid-cell-width in the y-direction, in meters. 1000 by default. 
!!
!> \param[in] grid\_type  Integer: which type of grid to use. If 0, the equidistant gnomonic cubed-sphere will be used. If 4, a doubly-periodic f-plane cartesian grid will be used. If -1, the grid is read from INPUT/grid\_spec.nc. Values 2, 3, 5, 6, and 7 are not supported and will likely not run. 0 by default. 
!!
!> \param[in] hybrid\_z  Logical: whether to use a hybrid-height coordinate, instead of the usual sigma-p coordinate. False by default. (Not currently maintained.)
!!
!> \param[in] make\_hybrid\_z  Logical: Converts the vertical coordinate to a hybrid-height coordinate, instead of the usual sigma-p coordinate. Requires hybrid\_z = True. False by default. 
!!
!> \param[in] p\_ref  Real: surface pressure used to construct a horizontally-uniform reference vertical pressure profile, used in some simple physics packages in the solo\_core and in the Rayleigh damping. *This should not be con- fused with the actual, horizontally-varying pressure levels used for all other dynamical calculations.* 1.e5 by default. *Changing this value is strongly discouraged.*
!!
!> \param[in] shift\_fac  Real: westward zonal rotation (or shift) of cubed-sphere grid from its natural orientation with cube face centers at 0, 90, 180, and 270 degrees longitude. The shift, in degrees, is 180/shift\_fac. This shift does not move the poles. By default this is set to 18, shifting the grid westward 180/18=10 degrees, so that the edges of the cube do not run through the mountains of Japan; all standard CM2.x, AM3, CM3, and HiRAM simulations use this orientation of the grid. Requires do\_schmidt =.false. 
!!
!> \param[in] stretch\_fac  Real: stretching factor for the Schmidt transformation. This is the factor by which tile 6 of the cubed sphere will be shrunk, with the grid size shrinking accordingly. 1 by default, which performs no grid stretching. Requires do\_schmidt  =.true. The model will crash if stretch\_fac is set to zero. Values of up to 40 have been found useful and stable for short-term cloud-scale integrations.
!!
!> \param[in] target\_lat  Real: latitude (in degrees) to which the center of tile 6 will be rotated; if stretching is done with stretch\_fac the center of the high-resolution part of the grid will be at this latitude.  -90 by default, which does no grid rotation (the Schmidt transformation rotates the south pole to the appropriate target). Requires do\_schmidt =.true. 
!!
!> \param[in] target\_lon  Real: longitude to which the center of tile 6 will be rotated. 0 by default. Requires do\_schmidt =.true.
!!
!> \param[in] nested  Logical: whether this is a nested grid. False by default. 
!!
!> \param[in] twowaynest  Logical: whether to use two-way nesting, the process by which the nested-grid solution can feed back onto the coarse-grid solution. False by default. 
!!
!> \param[in] parent\_grid\_num  Integer: Number of the parent to this nested grid. The coarsest grid in a simulation is numbered 1; further nested grids are numbered sequentially. Required to be a positive value if nested = True. Unless you are nesting inside of nested grids or running multiple (coarse) grids that themselves do not interact, this should be set to 1. -1 by default, indicating that this grid does not have a parent grid. 
!!
!> \param[in] parent\_tile  Integer: number of the tile (ie. face) in which this nested grid is found in its parent. Required to be a positive value if nested = true. If the parent grid is not a cubed sphere, or itself is a nested grid, this should be set to 1. If the parent grid has been rotated (using do\_schmidt) with the intent of centering the nested grid at target\_lat and target\_lon, then parent\_tile should be set to 6. 1 by default. 
!!
!> \param[in] refinement  Integer: refinement ratio of the nested grid. This is the number of times that each coarse-grid cell face will be divided into smaller segments on the nested grid. Required to be a positive integer if nested = true. Nested grids are aligned with the coarse grid, so non-integer refinements are not permitted. 3 by default. 
!!
!> \param[in] nestupdate  Integer: type of nested-grid update to use; details are given in model/fv\_nesting.F90. 0 by default.
!!
!>###A.1.7 Solver options
!!
!> \param[in] a2b\_ord  Integer: order of interpolation used by the pressure gradient force to interpolate cell-centered (A-grid) values to the grid corners. 4 by default (recommended), which uses fourth-order interpolation; otherwise second-order interpolation is used. 
!!
!> \param[in] beta  Real: Parameter specifying fraction of time-off-centering for backwards evaluation of the pressure gradient force. 0.0 by default, which produces a fully backwards evaluationthe pressure gradient force is entirely evaluated using the updated (time `n+1` dynamical fields. A value of 0.5 will equally weight the PGF determined at times `n` and `n+1`, but may not be stable; values larger than 0.45 are not recommended. A value of 0.4 is recommended for most hydrostatic simulations, which allows an improved representation of inertia-gravity waves in the tropics. In non-hydrostatic simulations using the semi-implicit solver (a\_imp > 0.5) the values of a\_imp and beta should add to 1, so that the time-centering is consistent between the PGF and the nonhydrostatic solver. Proper range is 0 to 0.45.
!!
!> \param[in] c2l\_ord  Integer: order of interpolation from the solvers native D-grid winds to latitude-longitude A-grid winds, which are used as input to the physics routines and for writing to history files. 4 by default (recommended); fourth-order interpolation is used unless c2l\_ord = 2. 
!!
!> \param[in] consv\_am Logical: whether to enable Angular momentum fixer. False by default.
!!
!> \param[in] consv\_te  Real: fraction of total energy lost during the adiabatic integration between calls of the physics, to be added back globally as heat; essentially the strength of the energy fixer in the physics. Note that this is a global energy fixer and cannot add back energy locally. The default algorithm increments the potential temperature so the pressure gradients are unchanged. 0 by default. Proper range is 0 to 1. 1 will restore the energy completely to its original value before entering the physics; a value of 0.7 *roughly* causes the energy fixer to compensate for the amount of energy changed by the physics in GFDL HiRAM or AM3. 
!!
!> \param[in] convert\_ke  Logical: If true, adds energy dissipated through mechanical damping to heat throughout the *entire* depth of the domain; if false (default) this is only done in the sponge layer at the top of the domain. This option is only enabled if d_con > 1.e-5. 
!!
!> \param[in] d\_con  Real: Fraction of kinetic energy lost to explicit damping to be converted to heat. Acts as a dissipative heating mechanism in the dynamical core. 0. by default. Proper range is 0 to 1. Note that this is a local, physically correct, energy fixer.
!!
!> \param[in] delt\_max  Real: maximum allowed magnitude of the dissipative heating rate, K s$^{-1}$; larger magnitudes are clipped to this amount. This can help avoid instability that can occur due to strong heating when d\_con $> 0$. A value of 0.008 (a rate equivalent to about 800 K/day) is sufficient to stabilize the model at 3-km resolution.  Set to 1 by default, which effectively disables this limitation.
!!
!> \param[in] fill\_dp  Logical: like fill except for  `p`, the hydrostatic pressure thickness. When the filling occurs a diagnostic message is printed out, which is helpful for diagnosing where the problem may be occurring. Typically if the pressure filling is needed a crash is inevitable, and thus this option is often better for debugging than as a safety valve. False by default. 
!!
!> \param[in] fv\_sg\_adj  Integer: timescale (in seconds) at which to remove two-delta-z instability when the local (between two adjacent levels) Richardson number is less than 1. This is achieved by local mixing, which conserves  mass, momentum, and total energy.  Values of 0 or smaller disable this feature. If n\_sponge $< 0$ then the mixing is applied only to the top n\_sponge layers of the domain. Set to -1  (inactive) by default. Proper range is 0 to 3600.
!!
!> \param[in] halo\_update\_type  Integer: which scheme to use for halo updates in multiprocessor simulations. If set to 1 non-blocking updates are used, which can improve simulation efficiency on some systems. Otherwise, blocking halo updates are performed. 1 by default.
!! 
!> \param[in] hord\_mt  Integer: horizontal advection scheme for momentum fluxes. A complete list of kord options is given in the table below. 9 by default, which uses the third-order piecewise-parabolic method with the monotonicity constraint of Huynh, which is less diffusive than other constraints. For hydrostatic simulation, 8 (the L04 monotonicity constraint) is recommended; for nonhydrostatic simulation, the completely unlimited (``linear'' or non-monotone) PPM scheme is recommended. If no monotonicity constraint is applied, enabling the flux damping (do\_vort\_damp = .true.) is highly recommended to control grid-scale noise. It is also recommended that hord\_mt, hord\_vt, hord\_tm, and hord\_dp use the same value, to ensure consistenttransport of all dynamical fields, unless a positivity constraint on mass advection (hord\_dp) is desired.
!!
!> \param[in] hord\_vt  Integer: horizontal advection scheme for absolute vorticity and for vertical velocity in nonhydrostatic simulations. 9 by default.
!!
!> \param[in] hord\_tm  Integer: horizontal advection scheme for potential temperature and layer thickness in nonhydrostatic simulations. 9 by default.
!!
!> \param[in] hord\_dp  Integer: horizontal advection scheme for mass. A positivity constraint may be warranted for hord\_dp but not strictly necessary. 9 by default.
!!
!> \param[in] hord\_tr  Integer: horizontal advection scheme for tracers. 12 by default. This value can differ from the other hord options since tracers are sub-cycled (if inline\_q == False) and require positive-definite advection to control the appearance of non-physical negative masses. 8 (fastest) or 10 (least diffusive) are typically recommended.
!!
!! 
!! hord | Advection Method
!! :----: |           ----
!! 5  | Fastest unlimited fifth-order scheme with built-in 2&Delta;x filter; not recommended for hord_tr. This is also the most accurate and least diffusive FV scheme available within FV&sup3; if monotonicity preservation is not a high priority.
!! 6  | *Developmental* PPM scheme with an intermediate-strength mono- tonicity constraint. More diffusive than 5.
!! 7  | 6, applying a 2&Delta;x filter and a positivity constraint
!! 8  | PPM with the constraint of Lin 2004
!! 9  | PPM with the Hunyh constraint
!! 10  | 9, with a 2&Delta;x filter, and the Huynh constraint applied only if a certain condition is met; otherwise unlimited
!! 
!> \param[in] kord\_mt  Integer: vertical remapping scheme for the winds. 8 by default; 9 is recommended as the safest option, although 10, and 11 can also be useful. See table below for a complete list of kord options.
!!
!> \param[in] kord\_tm  Integer: vertical remapping scheme for temperature. If positive (not recommended), then vertical remapping is performed on total energy instead of temperature (see remap\_t below). -8 by default. 
!!
!> \param[in] kord\_tr  Integer: vertical remapping scheme for tracers. 8 by default. 9 or 11 recommended. It is often recommended to use the same value for kord\_tr as for kord\_tm.
!!
!> \param[in] kord\_wz  Integer: vertical remapping scheme for vertical velocity in nonhydrostatic simulations. 8 by default; 9 recommended. It is also recommended to use the same value for kord\_wz as for kord\_mt. 
!! 
!! kord  | Vertical remapping reconstruction method
!! :---------: | --------------------------------
!!   4  | Monotone PPM 
!!   6  | PPM 
!!   7  | PPM with Hyunh’s second constraint (see L04)
!!   9  | Monotonic cubic spline with 2&Delta;z oscillations removed 
!! 10  | Selectively monotonic cubic spline, where local extrema are retained, with 2&Delta;z oscillations removed 
!! 11  | Non-monotonic (linear) cubic spline with 2&Delta;z oscillations removed; if an invalid value for kord is given, this scheme is used 
!! 13  | Monotonic cubic spline with 2&Delta;z oscillations removed
!! 16  | Non-monotonic cubic spline with a strong 2&Delta;z filter (similar to hord = 6).
!! 
!> \param[in] no\_dycore  Logical: disables execution of the dynamical core, only running the initialization, diagnostic, and I/O routines, and any physics that may be enabled. Essentially turns the model into a column physics model. False by default. 
!!
!> \param[in] remap\_t  Logical: whether the vertical remapping is performed on (virtual) temperature instead of (virtual) potential temperature. Since typically potential temperature increases exponentially from layer to layer near the top boundary, the cubic-spline interpolation in the vertical remapping will have difficulty with the exponential profile. Temperature does not have this problem and will often yield a more accurate result. True by default.
!!
!> \param[in] reproduce\_sum  Logical: uses an exactly-reproducible global sum operation performed when computing the global energy for consv\_te. This is used because the FMS routine mpp\_sum() is not bit-wise reproducible due to its handling of floating point arithmetic, and so can return different answers for (say) different processor layouts. True by default. 
!!
!>###A.1.8 Nonhydrostatic options
!!
!> \param[in] a\_imp  Real: Controls behavior of the non-hydrostatic solver. Values \> 0.5 enable the semi-implicit solver, in which the value of a\_imp controls the time-off-centering: use a\_imp = 1.0 for a fully backward time-stepping. For consistency, the sum of beta and a\_imp should be 1 when the semi-implicit solver is used. The semi-implicit algorithm is substantially more efficient except at very high (km-scale) resolutions with an acoustic time step of a few seconds or less. 0.75 by default. Proper values are 0, or between 0.5 and 1 Only used if hydrostatic =.false.
!!
!> \param[in] hydrostatic  Logical: whether to use the hydrostatic or nonhydrostatic solver. True by default.
!!
!> \param[in] make\_nh  Logical: Whether to re-initialize the nonhydrostatic state, by re-computing dz from hydrostatic balance and setting w to 0. False by default.
!!
!> \param[in] p\_fac  Real: Safety factor for minimum nonhydrostatic pressures, which will be limited so the full pressure is no less than p\_fac times the hydrostatic pressure. This is only of concern in mid-top or high-top models with very low pressures near the model top, and has no effect in most simulations. The pressure limiting activates only when model is in danger of blowup due to unphysical negative total pressures. 0.05 by default. Only used if hydrostatic =.false. and the semi-implicit solver is used. Proper range is 0 to 0.25.
!!
!> \param[in] use\_logp  Logical: Enables a variant of the Lin pressure-gradient force algorithm, which uses the logarithm of pressure instead of the Exner function (as in Lin 1997). This yields more accurate results for regions that are nearly isothermal. Ignored if hydrostatic = true. False by default.
!!
!>###A.1.9 Damping options
!!
!> \param[in] d2\_bg  Real: coefficient for background second-order divergence damping. This option remains active even if nord is nonzero. 0.0 by default. Proper range is 0 to 0.02. 
!!
!> \param[in] d2\_bg\_k1  Real: strength of second-order diffusion in the top sponge layer. 0.16 by default. This value, and d2\_bg\_k2, will be changed appropriately in the model (depending on the height of model top), so the actual damping may be very reduced. See atmos\_cubed\_sphere/model/dyn\_core.F90 for details. Recommended range is 0. to 0.2. Note that since diffusion is converted to heat if d\_con > 0 larger amounts of sponge-layer diffusion may be *less* stable.
!!
!> \param[in] d2\_bg\_k2  Real: strength of second-order diffusion in the second sponge layer from the model top. 0.02 by default. This value should be lower than d2\_bg\_k1.
!!
!> \param[in] d4\_bg  Real: Dimensionless coefficient for background higher-order divergence damping. 0.0 by default. If no second-order divergence damping is used, then values between 0.1 and 0.16 are recommended. Requires nord > 0. Note that the scaling for d4\_bg differs from that of d2\_bg; nord >= 1 and d4\_bg = 0.16 will be less diffusive than nord = 0 and d2\_bg = 0.02.
!!
!> \param[in] dddmp  Real: Dimensionless coefficient for the second-order Smagorinsky-type divergence damping. 0.0 by default. 0.2 (the Smagorinsky constant) is recommended if ICs are noisy.
!!
!> \param[in] d\_ext  Real: coefficient for external (barotropic) mode damping. 0.02 by default. Proper range is 0 to 0.02. A value of 0.01 or 0.02 may help improve the models maximum stable time step in low-resolution (2-degree or poorer) simulations; otherwise a value of 0 is recommended.
!!
!> \param[in] do\_vort\_damp  Logical: whether to apply flux damping (of strength governed by vtdm4) to the fluxes of vorticity, air mass, and nonhydrostatic vertical velocity (there is no dynamically correct way to add explicit diffusion to the tracer fluxes). The form is the same as is used for the divergence damping, including the same order (from nord) damping, unless nord = 0, in which case this damping is fourth-order, or if nord = 3, in which case this damping is sixth-order (instead of eighth-order). We recommend enabling this damping when the linear or non-monotonic horizontal advection schemes are enabled, but is unnecessary and not recommended when using monotonic advection. False by default. 
!!
!> \param[in] n\_sponge  Integer: controls the number of layers at the upper boundary on which the 2&Delta;x filter is applied. *This does not control the sponge layer.* 0 by default. 
!!
!> \param[in] nord  Integer: order of divergence damping: 0 for second-order; 1 for fourth-order (default); 2 for sixth-order; 3 for eighth-order. Sixth-order may yield a better solution for low resolutions (one degree or coarser) by virtue of it being more scale-selective and will not damp moderately-well-resolved disturbances as much as does lower-order damping. 
!!
!> \param[in] nord\_tr  Integer: Order of tracer damping; values mean the same as for nord. 0 by default. 
!!
!> \param[in] rf\_cutoff  Real: pressure below which no Rayleigh damping is applied if tau > 0.
!!
!> \param[in] rf\_fast  Logical: option controlling whether to apply Rayleigh damping (for tau > 0) on the dynamic/acoustic timestep rather than on the physics timestep. This can help stabilize the model by applying the damping more weakly more frequently, so the instantaneous amount of damping (and thereby heat added) is reduced.  .false. by default, which applies the Rayleigh drag every physics timestep.
!!
!> \param[in] tau  Real: time scale (in days) for Rayleigh friction applied to horizontal and vertical winds; lost kinetic energy is converted to heat, except on nested grids. 0.0 by default, which disables damping. Larger values yield less damping. For models with tops at 1mb or lower values between 10 and 30 are useful for preventing overly-strong polar night jets; for higher-top hydrostatic models values between 5 and 15 should be considered; and for non-hydrostatic models values of 10 or less should be considered, with smaller values for higher-resolution.
!!
!> \param[in] vtdm4  Real: coefficient for background other-variable damping. The value of vtdm4 should be less than that of d4\_bg. A good first guess for vtdm4 is about one-third the value of d4\_bg. 0.0 by default. Requires do\_vort\_damp to be .true. Disabled for values less than 1.e-3. Other-variable damping should not be used if a monotonic horizontal advection scheme is used.
!!
!>##A.2 Entries in coupler\_nml
!!
!> \param[in] months, days,hours,minutes,seconds  Integer: length of the model integration in the corresponding units. All are 0 by default, which initializes the model and then immediately writes out the restart files and halts.
!!
!> \param[in] dt\_atmos  Integer: time step for the largest atmosphere model loop, corresponding to the frequency at which the top level routine in the dynamics is called, and the physics timestep. Must be set.
!!
!> \param[in] current\_date  Integer(6): initialization date (in the chosen calendar) for the model, in year, month, day, hour, minute, and second. (0,0,0,0,0,0) by default, a value that is useful for control integrations of coupled models.
!!
!> \param[in] calendar  Character(17): calendar selection; JULIAN is typically recommended, although the other values (THIRTY\_DAY\_MONTHS, NOLEAP, NO\_CALENDAR) have particular uses in idealized models. Must be set.
!!
!> \param[in] force\_date\_from\_namelist  Logical: if .true., will read the initialization date from the namelist variable current\_date, rather than taking the value from a restart file. If the model is being cold-started (such as what is typically but not necessarily done if external\_ic = .true.) then the initialization date must be specified in current\_date, otherwise the model will stop. .false. by default.
!!
!> \param[in] atmos\_nthreads  Integer: number of threads for OpenMP multi-threading. 1 by default.
!!
!> \param[in] use\_hyper\_thread  Logical: indicates that some of the threads in atmos\_nthreads may be hyperthreads. .false. by default.
!!
!> \param[in] ncores\_per\_node  Integer: number of processor codes per physical compute node. Used when setting up hyperthreading to determine number of virtual vs. hardware threads. 0 by default.
!!
!> \param[in] debug\_affinity  Logical: if .true. prints out a message describing cpu affinity characteristics while initializing OpenMP. .false. by default.
!!
!> \param[in] restart\_days, restart\_secs]  Integer: frequency at which to write out "intermediate" restart files, which are useful for checkpointing in the middle of a long run, or to be able to diagnose problems during the model integration. Both are 0 by default, in which case intermediate restarts are not written out.
!!
!>##A.3 Entries in external\_ic\_nml
!!
!> \param[in] filtered\_terrain  Logical: whether to use the terrain filtered by the preprocessing tools rather than the raw terrain. .true. by default. Only active if nggps\_ic = .true. or ecmwf\_ic = .true.
!!
!> \param[in] levp  Integer: number of levels in the input (remapped) initial conditions. 64 by default. Only active if nggps\_ic = .true.
!!
!> \param[in] checker\_tr  Logical: whether to enable the ``checkerboard'' tracer test. .false. by default. Only active if nggps\_ic = .true.
!!
!> \param[in] nt\_checker  Integer: number of tracers (at the end of the list of tracers defined in field\_table) to initialize with an idealized ``checkerboard'' pattern, with values of either 0 or 1. This is intended to test the monotonicity or positivity constraint in the advection scheme. 0 by default. Only active if nggps\_ic = .true.
!!
!>##A.4 Entries in surf\_map\_nml
!!
!> \param[in] surf\_file  Character(len=128): File containing topography data. This file must be in NetCDF format. INPUT/topo1min.nc by default. (Previous versions of the model have used 5 minute USGS data, which is not recommended.) 
!!
!> \param[in] nlon  Integer: Size of the longitude dimension in topography data; not used. 
!!
!> \param[in] nlat  Integer: Size of the latitude dimension in topography data; not used. 
!!
!> \param[in] zero\_ocean Logical: whether to prevent the smoothing from extending topography out into the ocean. True by default. 
!!
!> \param[in] zs\_filter  Logical: whether to apply smoothing to the topography. True by default. 
!!
!>##A.5 Entries in fv\_grid\_nml
!!
!> \param[in] grid\_name Character(len=80): Name of the grid either being read in (if grid\_spec = -1) or being created. This is only used for writing out a binary file in the directory from which the model is run. Gnomonic by default. 
!!
!> \param[in] grid\_file  Character(len=120): If grid\_type = -1 the name of the grid\_spec file to read in. INPUT/grid\_spec.nc by default; other values will not work. 
!!
!>##A.6 Entries in test\_case\_nml
!!
!> \param[in] test\_case  Integer: number of the idealized test case to run. A number of nest cases are defined in tools/test\_cases.F90, of which numbers 19 are intended for the shallow-water model. Requires warm\_start =.false. 11 by default; this creates a resting atmosphere with a very basic thermodynamic profile, with topography. If you wish to initialize an Aquaplanet simulation (no topography) set to 14. 
!!
!> \param[in] alpha  Real: In certain shallow-water test cases specifies the angle (in fractions of a rotation, so 0.25 is a 45-degree rotation) through which the basic state is rotated. 0 by default. 
!!
!>##A.7  Entries in nest\_nml
!!
!> \param[in] ngrids  Integer: number of grids in this simulation. 1 by default.  (The variable ntiles has the same effect, but its use is not recommended as it can be confused with the six tiles of the cubed sphere, which in this case form only one grid.)
!!
!> \param[in] nest\_pes  Integer(100): array carrying number of PEs assigned to each grid, in order. Must be set if ngrids > 1. 
!!
!> \param[in] p\_split  Integer: number of times to sub-cycle dynamics, performing nested-grid BC interpolation and (if twowaynest ==.true.) two-way updating at the end of each set of dynamics calls. If p\_split  > 1 the user should decrease k\_split appropriately so the remapping and dynamics time steps remain the same. 1 by default.
!!
!>##A.8 Entries in nggps\_diag\_nml
!!
!> \param[in] fdiag  Real(1028): Array listing the diagnostic output times (in hours) for the GFS physics. This can either be a list of times after initialization, or an interval if only the first entry is nonzero. All 0 by default, which will result in no outputs. 
!!
!>##A.9 Entries in atmos\_model\_nml (for UFS)
!!
!> \param[in] blocksize  Integer: Number of columns in each ``block'' sent to the physics. OpenMP threading is done over the number of blocks. For best performance this number should divide the number of grid cells per processor ( (npx-1)*(npy-1) /(layout\_x)*(layout\_y) ) and be small enough so the data can best fit into cache?values around 32 appear to be optimal on Gaea. 1 by default
!!
!> \param[in] chksum\_debug  Logical: whether to compute checksums for all variables passed into the GFS physics, before and after each physics timestep. This is very useful for reproducibility checking. .false. by default.
!!
!> \param[in] dycore\_only  Logical: whether only the dynamical core (and not the GFS physics) is executed when running the model, essentially running the model as a solo dynamical core. .false. by default.
!!
!>##A.10 Entries in fms\_nml
!!
!> \param[in] domains\_stack\_size  Integer: size (in bytes) of memory array reserved for domains. For large grids or reduced processor counts this can be large (>10 M); if it is not large enough the model will stop and print a recommended value of the stack size. Default is 0., reverting to the default set in MPP (which is probably not large enough for modern applications).
!!
!>##A.11
!!
!> \param[in] regional Logical: Controls whether this is a regional domain (and thereby needs external BC inputs
!!
!> \param[in] bc_update_interval Integer: Default setting for interval (hours) between external regional BC data files.
!!
!> \param[in] read_increment  Logical: Read in analysis increment and add to restart following are namelist parameters for Stochastic Energy Baskscatter dissipation estimate. This is useful as part of a data-assimilation cycling system or to use native restarts from the six-tile first guess, after which the analysis increment can be applied. 
!!
!> \param[in] do_sat_adj  Logical: The same as fast_sat_adj = .false.  has fast saturation adjustments
!!
!!
!!
!! 
!> @{ 
   namelist /fv_grid_nml/ grid_name, grid_file
   namelist /fv_core_nml/npx, npy, ntiles, npz, npz_rst, layout, io_layout, ncnst, nwat,  &
                         use_logp, p_fac, a_imp, k_split, n_split, m_split, q_split, print_freq, write_3d_diags, do_schmidt,  &
                         hord_mt, hord_vt, hord_tm, hord_dp, hord_tr, shift_fac, stretch_fac, target_lat, target_lon, &
                         kord_mt, kord_wz, kord_tm, kord_tr, fv_debug, fv_land, nudge, do_sat_adj, do_f3d, &
                         external_ic, read_increment, ncep_ic, nggps_ic, ecmwf_ic, use_new_ncep, use_ncep_phy, fv_diag_ic, &
                         external_eta, res_latlon_dynamics, res_latlon_tracers, scale_z, w_max, z_min, lim_fac, &
                         dddmp, d2_bg, d4_bg, vtdm4, trdm2, d_ext, delt_max, beta, non_ortho, n_sponge, &
                         warm_start, adjust_dry_mass, mountain, d_con, ke_bg, nord, nord_tr, convert_ke, use_old_omega, &
                         dry_mass, grid_type, do_Held_Suarez, do_reed_physics, reed_cond_only, &
                         consv_te, fill, filter_phys, fill_dp, fill_wz, consv_am, RF_fast, &
                         range_warn, dwind_2d, inline_q, z_tracer, reproduce_sum, adiabatic, do_vort_damp, no_dycore,   &
                         tau, tau_h2o, rf_cutoff, nf_omega, hydrostatic, fv_sg_adj, breed_vortex_inline,  &
                         na_init, nudge_dz, hybrid_z, Make_NH, n_zs_filter, nord_zs_filter, full_zs_filter, reset_eta,         &
                         pnats, dnats, a2b_ord, remap_t, p_ref, d2_bg_k1, d2_bg_k2,  &
                         c2l_ord, dx_const, dy_const, umax, deglat,      &
                         deglon_start, deglon_stop, deglat_start, deglat_stop, &
                         phys_hydrostatic, use_hydro_pressure, make_hybrid_z, old_divg_damp, add_noise, &
                         nested, twowaynest, parent_grid_num, parent_tile, nudge_qv, &
                         refinement, nestbctype, nestupdate, nsponge, s_weight, &
                         ioffset, joffset, check_negative, nudge_ic, halo_update_type, gfs_phil, agrid_vel_rst,     &
                         do_uni_zfull, adj_mass_vmr, fac_n_spl, fhouri, regional, bc_update_interval

   namelist /test_case_nml/test_case, bubble_do, alpha, nsolitons, soliton_Umax, soliton_size
#ifdef MULTI_GASES
   namelist /multi_gases_nml/ rilist,cpilist
#endif


   pe_counter = mpp_root_pe()

! Make alpha = 0 the default:
   alpha = 0.
   bubble_do = .false.
   test_case = 11   ! (USGS terrain)

#ifdef INTERNAL_FILE_NML
! Read Main namelist
   read (input_nml_file,fv_grid_nml,iostat=ios)
   ierr = check_nml_error(ios,'fv_grid_nml')
#else
   f_unit=open_namelist_file()
   rewind (f_unit)
! Read Main namelist
   read (f_unit,fv_grid_nml,iostat=ios)
   ierr = check_nml_error(ios,'fv_grid_nml')
   call close_file(f_unit)
#endif

   call write_version_number ( 'FV_CONTROL_MOD', version )
   unit = stdlog()
   write(unit, nml=fv_grid_nml)

   do n=1,size(Atm)

      call switch_current_Atm(Atm(n), .false.)
      call setup_pointers(Atm(n))
      Atm(n)%grid_number = n
      if (grids_on_this_pe(n)) then
         call fv_diag_init_gn(Atm(n))
      endif

#ifdef INTERNAL_FILE_NML
   ! Set input_file_nml for correct parent/nest initialization
      if (n > 1) then
         write(nested_grid_filename,'(A4, I2.2)') 'nest', n
         call read_input_nml(nested_grid_filename)
      endif
   ! Read FVCORE namelist 
      read (input_nml_file,fv_core_nml,iostat=ios)
      ierr = check_nml_error(ios,'fv_core_nml')
#ifdef MULTI_GASES
      if( is_master() ) print *,' enter multi_gases: ncnst = ',ncnst
      allocate (rilist(0:ncnst))
      allocate (cpilist(0:ncnst))
      rilist     =    0.0
      cpilist    =    0.0
      rilist(0)  = rdgas
      rilist(1)  = rvgas
      cpilist(0) = cp_air
      cpilist(1) = 4*cp_air
   ! Read multi_gases namelist
      read (input_nml_file,multi_gases_nml,iostat=ios)
      ierr = check_nml_error(ios,'multi_gases_nml')
#endif
   ! Read Test_Case namelist
      read (input_nml_file,test_case_nml,iostat=ios)
      ierr = check_nml_error(ios,'test_case_nml')

   ! Reset input_file_nml to default behavior
      call read_input_nml
#else
      if (size(Atm) == 1) then
         f_unit = open_namelist_file()
      else if (n == 1) then
         f_unit = open_namelist_file('input.nml')
      else 
         write(nested_grid_filename,'(A10, I2.2, A4)') 'input_nest', n, '.nml'
         f_unit = open_namelist_file(nested_grid_filename)
      endif

   ! Read FVCORE namelist 
      read (f_unit,fv_core_nml,iostat=ios)
      ierr = check_nml_error(ios,'fv_core_nml')

#ifdef MULTI_GASES
      if( is_master() ) print *,' enter multi_gases: ncnst = ',ncnst
      allocate (rilist(0:ncnst))
      allocate (cpilist(0:ncnst))
      rilist     =    0.0
      cpilist    =    0.0
      rilist(0)  = rdgas
      rilist(1)  = rvgas
      cpilist(0) = cp_air
      cpilist(1) = 4*cp_air
   ! Read multi_gases namelist
      rewind (f_unit)
      read (f_unit,multi_gases_nml,iostat=ios)
      ierr = check_nml_error(ios,'multi_gases_nml')
#endif
   ! Read Test_Case namelist
      rewind (f_unit)
      read (f_unit,test_case_nml,iostat=ios)
      ierr = check_nml_error(ios,'test_case_nml')
      call close_file(f_unit)
#endif         
      write(unit, nml=fv_core_nml)
      write(unit, nml=test_case_nml)
#ifdef MULTI_GASES
      write(unit, nml=multi_gases_nml)
      call multi_gases_init(ncnst,nwat)
#endif

      if (len_trim(grid_file) /= 0) Atm(n)%flagstruct%grid_file = grid_file
      if (len_trim(grid_name) /= 0) Atm(n)%flagstruct%grid_name = grid_name
      if (len_trim(res_latlon_dynamics) /= 0) Atm(n)%flagstruct%res_latlon_dynamics = res_latlon_dynamics
      if (len_trim(res_latlon_tracers)  /= 0) Atm(n)%flagstruct%res_latlon_tracers = res_latlon_tracers

      !*** single tile for Cartesian grids
      if (grid_type>3) then
         ntiles=1
         non_ortho = .false.
         nf_omega = 0
      endif

      if (.not. (nested .or. regional)) Atm(n)%neststruct%npx_global = npx

      ! Define n_split if not in namelist
      if (ntiles == 6) then
         dimx = 4.0*(npx-1)
         if ( hydrostatic ) then
            if ( npx >= 120 ) ns0 = 6
         else
            if ( npx <= 45 ) then
               ns0 = 6
            elseif ( npx <= 90 ) then
               ns0 = 7
            else
               ns0 = 8
            endif
         endif
      else
         dimx = max ( npx, 2*(npy-1) )
      endif
          
      if (grid_type < 4) then
         n0split = nint ( ns0*abs(dt_atmos)*dimx/(dt0*dim0) + 0.49 )
      elseif (grid_type == 4 .or. grid_type == 7) then
         n0split = nint ( 2.*umax*dt_atmos/sqrt(dx_const**2 + dy_const**2) + 0.49 )
      elseif (grid_type == 5 .or. grid_type == 6) then
         if (grid_type == 6) then
            deglon_start = 0.; deglon_stop  = 360.
         endif
         dl = (deglon_stop-deglon_start)*pi/(180.*(npx-1))
         dp = (deglat_stop-deglat_start)*pi/(180.*(npy-1))

         dxmin=dl*radius*min(cos(deglat_start*pi/180.-ng*dp),   &
                             cos(deglat_stop *pi/180.+ng*dp))
         dymin=dp*radius
         n0split = nint ( 2.*umax*dt_atmos/sqrt(dxmin**2 + dymin**2) + 0.49 )
      endif
      n0split = max ( 1, n0split )

      if ( n_split == 0 ) then
           n_split = nint( real(n0split)/real(k_split*abs(p_split)) * stretch_fac + 0.5 )
           if(is_master()) write(*,*) 'For k_split (remapping)=', k_split
           if(is_master()) write(*,198) 'n_split is set to ', n_split, ' for resolution-dt=',npx,npy,ntiles,dt_atmos
      else
          if(is_master()) write(*,199) 'Using n_split from the namelist: ', n_split
      endif
      if (is_master() .and. n == 1 .and. abs(p_split) > 1) then
         write(*,199) 'Using p_split = ', p_split
      endif

      if (Atm(n)%neststruct%nested) then
         do i=1,n-1
            if (Atm(i)%grid_number == parent_grid_num) then
               Atm(n)%parent_grid => Atm(i)
               exit
            end if
         end do
         if (.not. associated(Atm(n)%parent_grid)) then
            write(errstring,'(2(A,I3))')  "Could not find parent grid #", parent_grid_num, ' for grid #', n
            call mpp_error(FATAL, errstring)
         end if

         !Note that if a gnomonic grid has a parent it is a NESTED gnomonic grid and therefore only has one tile
         if ( Atm(n)%parent_grid%flagstruct%grid_type < 3 .and. &
              .not. associated(Atm(n)%parent_grid%parent_grid)) then
            if (parent_tile > 6 .or. parent_tile < 1) then
               call mpp_error(FATAL, 'parent tile must be between 1 and 6 if the parent is a cubed-sphere grid')
            end if
         else
            if (parent_tile /= 1) then
               call mpp_error(FATAL, 'parent tile must be 1 if the parent is not a cubed-sphere grid')
            end if
         end if

         if ( refinement < 1 ) call mpp_error(FATAL, 'grid refinement must be positive')

         if (nestupdate == 1 .or. nestupdate == 2) then

            if (mod(npx-1,refinement) /= 0 .or. mod(npy-1,refinement) /= 0) then
               call mpp_error(WARNING, 'npx-1 or npy-1 is not evenly divisible by the refinement ratio; averaging update cannot be mass-conservative.')
            end if

         end if

         if ( consv_te > 0.) then
            call mpp_error(FATAL, 'The global energy fixer cannot be used on a nested grid. consv_te must be set to 0.')
         end if

         Atm(n)%neststruct%refinement_of_global = Atm(n)%neststruct%refinement * Atm(n)%parent_grid%neststruct%refinement_of_global
         max_refinement_of_global = max(Atm(n)%neststruct%refinement_of_global,max_refinement_of_global)
         Atm(n)%neststruct%npx_global = Atm(n)%neststruct%refinement * Atm(n)%parent_grid%neststruct%npx_global

      else
         Atm(n)%neststruct%ioffset                = -999
         Atm(n)%neststruct%joffset                = -999   
         Atm(n)%neststruct%parent_tile            = -1      
         Atm(n)%neststruct%refinement             = -1
      end if

      if (Atm(n)%neststruct%nested) then
         if (Atm(n)%flagstruct%grid_type >= 4 .and. Atm(n)%parent_grid%flagstruct%grid_type >= 4) then
            Atm(n)%flagstruct%dx_const = Atm(n)%parent_grid%flagstruct%dx_const / real(Atm(n)%neststruct%refinement)
            Atm(n)%flagstruct%dy_const = Atm(n)%parent_grid%flagstruct%dy_const / real(Atm(n)%neststruct%refinement)
         end if
      end if


!----------------------------------------
! Adjust divergence damping coefficients:
!----------------------------------------
!      d_fac = real(n0split)/real(n_split)
!      dddmp = dddmp * d_fac
!      d2_bg = d2_bg * d_fac
!      d4_bg = d4_bg * d_fac
!      d_ext = d_ext * d_fac
!      vtdm4 = vtdm4 * d_fac
      if (old_divg_damp) then
        if (is_master()) write(*,*) " fv_control: using original values for divergence damping "
        d2_bg_k1 = 6.         ! factor for d2_bg (k=1)  - default(4.)
        d2_bg_k2 = 4.         ! factor for d2_bg (k=2)  - default(2.)
        d2_divg_max_k1 = 0.02 ! d2_divg max value (k=1) - default(0.05)
        d2_divg_max_k2 = 0.01 ! d2_divg max value (k=2) - default(0.02)
        damp_k_k1 = 0.        ! damp_k value (k=1)      - default(0.05)
        damp_k_k2 = 0.        ! damp_k value (k=2)      - default(0.025)
      elseif (n_sponge == 0 ) then
        if ( d2_bg_k1 > 1. ) d2_bg_k1 = 0.20
        if ( d2_bg_k2 > 1. ) d2_bg_k2 = 0.015
      endif

!     if ( beta < 1.e-5 ) beta = 0.   ! beta < 0 is used for non-hydrostatic "one_grad_p"

      if ( .not.hydrostatic ) then
         if ( m_split==0 ) then
            m_split = 1. + abs(dt_atmos)/real(k_split*n_split*abs(p_split))
            if (abs(a_imp) < 0.5) then
               if(is_master()) write(*,199) 'm_split is set to ', m_split
            endif
         endif
         if(is_master()) then
            write(*,*) 'Off center implicit scheme param=', a_imp
            write(*,*) ' p_fac=', p_fac
         endif
      endif

      if(is_master()) then
         if (n_sponge >= 0) write(*,199) 'Using n_sponge : ', n_sponge
         write(*,197) 'Using non_ortho : ', non_ortho
      endif

 197  format(A,l7)
 198  format(A,i2.2,A,i4.4,'x',i4.4,'x',i1.1,'-',f9.3)
 199  format(A,i3.3)

      if (.not. (nested .or. regional)) alpha = alpha*pi


      allocate(Atm(n)%neststruct%child_grids(size(Atm)))
      Atm(N)%neststruct%child_grids = .false.

      !Broadcast data

      !Check layout

   enddo

   !Set pelists
   do n=1,size(Atm)
      if (ANY(Atm(n)%pelist == gid)) then
            call mpp_set_current_pelist(Atm(n)%pelist)
            call mpp_get_current_pelist(Atm(n)%pelist, commID=commID)
            call mp_start(commID,halo_update_type)
      endif

      if (Atm(n)%neststruct%nested) then
         Atm(n)%neststruct%parent_proc = ANY(Atm(n)%parent_grid%pelist == gid)
         Atm(n)%neststruct%child_proc = ANY(Atm(n)%pelist == gid)
      endif
   enddo

   do n=1,size(Atm)
      
      call switch_current_Atm(Atm(n),.false.)
      call setup_pointers(Atm(n))
      !! CLEANUP: WARNING not sure what changes to domain_decomp may cause
      call domain_decomp(npx,npy,ntiles,grid_type,nested,Atm(n),layout,io_layout)
   enddo

   !!! CLEANUP: This sets the pelist to ALL, which is also
   !!! required for the define_nest_domains step in the next loop.
   !!! Later the pelist must be reset to the 'local' pelist.
   call broadcast_domains(Atm)

   do n=1,size(Atm)
      call switch_current_Atm(Atm(n))
      call setup_pointers(Atm(n))

      if (nested) then
         if (mod(npx-1 , refinement) /= 0 .or. mod(npy-1, refinement) /= 0) &
              call mpp_error(FATAL, 'npx or npy not an even refinement of its coarse grid.')
         
         !Pelist needs to be set to ALL (which should have been done
         !in broadcast_domains) to get this to work
         call mpp_define_nest_domains(Atm(n)%neststruct%nest_domain, Atm(n)%domain, Atm(parent_grid_num)%domain, &
              7, parent_tile, &
              1, npx-1, 1, npy-1,                  & !Grid cells, not points
              ioffset, ioffset + (npx-1)/refinement - 1, &
              joffset, joffset + (npy-1)/refinement - 1,         &
              (/ (i,i=0,mpp_npes()-1)  /), extra_halo = 0, name="nest_domain") !What pelist to use?
         call mpp_define_nest_domains(Atm(n)%neststruct%nest_domain, Atm(n)%domain, Atm(parent_grid_num)%domain, &
              7, parent_tile, &
              1, npx-1, 1, npy-1,                  & !Grid cells, not points
              ioffset, ioffset + (npx-1)/refinement - 1, &
              joffset, joffset + (npy-1)/refinement - 1,         &
              (/ (i,i=0,mpp_npes()-1)  /), extra_halo = 0, name="nest_domain") !What pelist to use?
!              (/ (i,i=0,mpp_npes()-1)  /), extra_halo = 2, name="nest_domain_for_BC") !What pelist to use?

         Atm(parent_grid_num)%neststruct%child_grids(n) = .true.

         if (Atm(n)%neststruct%nestbctype > 1) then

            call mpp_error(FATAL, 'nestbctype > 1 not yet implemented')

            !This check is due to a bug which has not yet been identified. Beware.
!            if (Atm(n)%parent_grid%flagstruct%hord_tr == 7) &
!                 call mpp_error(FATAL, "Flux-form nested  BCs (nestbctype > 1) should not use hord_tr == 7 (on parent grid), since there is no guarantee of tracer mass conservation with this option.")

!!$            if (Atm(n)%flagstruct%q_split > 0 .and. Atm(n)%parent_grid%flagstruct%q_split > 0) then
!!$               if (mod(Atm(n)%flagstruct%q_split,Atm(n)%parent_grid%flagstruct%q_split) /= 0) call mpp_error(FATAL, &
!!$                    "Flux-form nested BCs (nestbctype > 1) require q_split on the nested grid to be evenly divisible by that on the coarse grid.")
!!$            endif
!!$            if (mod((Atm(n)%npx-1),Atm(n)%neststruct%refinement) /= 0 .or. mod((Atm(n)%npy-1),Atm(n)%neststruct%refinement) /= 0) call mpp_error(FATAL, &
!!$                 "Flux-form nested BCs (nestbctype > 1) requires npx and npy to be one more than a multiple of the refinement ratio.")
!!$            Atm(n)%parent_grid%neststruct%do_flux_BCs = .true.
!!$            if (Atm(n)%neststruct%nestbctype == 3 .or. Atm(n)%neststruct%nestbctype == 4) Atm(n)%parent_grid%neststruct%do_2way_flux_BCs = .true.
            Atm(n)%neststruct%upoff = 0
         endif

      end if

      do nn=1,size(Atm)
         if (n == 1) allocate(Atm(nn)%neststruct%nest_domain_all(size(Atm)))
         Atm(nn)%neststruct%nest_domain_all(n) = Atm(n)%neststruct%nest_domain
      enddo

   end do

   do n=1,size(Atm)
      if (ANY(Atm(n)%pelist == gid)) then
         call mpp_set_current_pelist(Atm(n)%pelist)
      endif
   enddo

  end subroutine run_setup
  subroutine init_nesting(Atm, grids_on_this_pe, p_split)
    
    type(fv_atmos_type), intent(inout), allocatable :: Atm(:)
   logical, allocatable, intent(INOUT) :: grids_on_this_pe(:)
    integer, intent(INOUT) :: p_split
    character(100) :: pe_list_name
    integer :: nest_pes(100)
    integer :: n, npes, ntiles, pecounter, i
    integer, allocatable :: pelist(:)
    integer :: f_unit, ios, ierr

    !This is an OPTIONAL namelist, that needs to be read before everything else
    namelist /nest_nml/ ngrids, ntiles, nest_pes, p_split

    call mp_assign_gid

    nest_pes = 0
    ntiles = -999

#ifdef INTERNAL_FILE_NML
      read (input_nml_file,nest_nml,iostat=ios)
      ierr = check_nml_error(ios,'nest_nml')
#else
      f_unit=open_namelist_file()
      rewind (f_unit)
      read (f_unit,nest_nml,iostat=ios)
      ierr = check_nml_error(ios,'nest_nml')
      call close_file(f_unit)
#endif

      if (ntiles /= -999) ngrids = ntiles
      if (ngrids > 10) call mpp_error(FATAL, "More than 10 nested grids not supported")

      allocate(Atm(ngrids))
    
      allocate(grids_on_this_pe(ngrids))
      grids_on_this_pe = .false. !initialization

      npes = mpp_npes()

      ! Need to get a global pelist to send data around later?
      allocate( pelist_all(npes) )
      pelist_all = (/ (i,i=0,npes-1) /)
      pelist_all = pelist_all + mpp_root_pe()

      if (ngrids == 1) then

         !Set up the single pelist
         allocate(Atm(1)%pelist(npes))
         Atm(1)%pelist = (/(i, i=0, npes-1)/)
         Atm(1)%pelist = Atm(1)%pelist + mpp_root_pe()
         call mpp_declare_pelist(Atm(1)%pelist)
         call mpp_set_current_pelist(Atm(1)%pelist)
         !Now set in domain_decomp
         !masterproc = Atm(1)%pelist(1)
         call setup_master(Atm(1)%pelist)
         grids_on_this_pe(1) = .true.
         Atm(1)%npes_this_grid = npes

      else

         pecounter = mpp_root_pe()
         do n=1,ngrids
            if (n == 1) then
               pe_list_name = ''
            else
               write(pe_list_name,'(A4, I2.2)') 'nest', n
            endif

            if (nest_pes(n) == 0) then
               if (n < ngrids) call mpp_error(FATAL, 'Only nest_pes(ngrids) in nest_nml can be zero; preceeding values must be nonzero.')
               allocate(Atm(n)%pelist(npes-pecounter))
               Atm(n)%pelist = (/(i, i=pecounter, npes-1)/)
               if (n > 1) then
                  call mpp_declare_pelist(Atm(n)%pelist, trim(pe_list_name))
                  !Make sure nested-grid input file exists
                  if (.not. file_exist('input_'//trim(pe_list_name)//'.nml')) then
                     call mpp_error(FATAL, "Could not find nested grid namelist input_"//trim(pe_list_name)//".nml")
                  endif
               endif
               exit
            else
               allocate(Atm(n)%pelist(nest_pes(n)))
               Atm(n)%pelist = (/ (i, i=pecounter, pecounter+nest_pes(n)-1) /)
               if (Atm(n)%pelist(nest_pes(n)) >= npes) then
                  call mpp_error(FATAL, 'PEs assigned by nest_pes in nest_nml exceeds number of available PEs.')
               endif

               call mpp_declare_pelist(Atm(n)%pelist, trim(pe_list_name))
               !Make sure nested-grid input file exists
               if (n > 1) then
                  if (.not. file_exist('input_'//trim(pe_list_name)//'.nml')) then
                     call mpp_error(FATAL, "Could not find nested grid namelist input_"//trim(pe_list_name)//".nml")
                  endif
               endif
               pecounter = pecounter+nest_pes(n)
            endif
         enddo

         !Set pelists
         do n=1,ngrids
            Atm(n)%npes_this_grid = size(Atm(n)%pelist)
            if (ANY(gid == Atm(n)%pelist)) then
                  call mpp_set_current_pelist(Atm(n)%pelist)
                  !now set in domain_decomp
                  !masterproc = Atm(n)%pelist(1)
                  call setup_master(Atm(n)%pelist)
                  grids_on_this_pe(n) = .true.
                  exit
               endif
         enddo

         if (pecounter /= npes) then
            call mpp_error(FATAL, 'nest_pes in nest_nml does not assign all of the available PEs.')
         endif
      endif

      !Layout is checked later, in fv_control

  end subroutine init_nesting

!>@brief The subroutine 'setup_pointers' associates the MODULE flag pointers
!! with the ARRAY flag variables for the grid active on THIS pe so the flags
!! can be read in from the namelist.
  subroutine setup_pointers(Atm)

    type(fv_atmos_type), intent(INOUT), target :: Atm

    !This routine associates the MODULE flag pointers with the ARRAY flag variables for the grid active on THIS pe so the flags can be read in from the namelist.

     res_latlon_dynamics           => Atm%flagstruct%res_latlon_dynamics
     res_latlon_tracers            => Atm%flagstruct%res_latlon_tracers

     grid_type                     => Atm%flagstruct%grid_type
     grid_name                     => Atm%flagstruct%grid_name
     grid_file                     => Atm%flagstruct%grid_file
     hord_mt                       => Atm%flagstruct%hord_mt
     kord_mt                       => Atm%flagstruct%kord_mt
     kord_wz                       => Atm%flagstruct%kord_wz
     hord_vt                       => Atm%flagstruct%hord_vt
     hord_tm                       => Atm%flagstruct%hord_tm
     hord_dp                       => Atm%flagstruct%hord_dp
     kord_tm                       => Atm%flagstruct%kord_tm
     hord_tr                       => Atm%flagstruct%hord_tr
     kord_tr                       => Atm%flagstruct%kord_tr
     scale_z                       => Atm%flagstruct%scale_z
     w_max                         => Atm%flagstruct%w_max
     z_min                         => Atm%flagstruct%z_min
     lim_fac                       => Atm%flagstruct%lim_fac
     nord                          => Atm%flagstruct%nord
     nord_tr                       => Atm%flagstruct%nord_tr
     dddmp                         => Atm%flagstruct%dddmp
     d2_bg                         => Atm%flagstruct%d2_bg
     d4_bg                         => Atm%flagstruct%d4_bg
     vtdm4                         => Atm%flagstruct%vtdm4
     trdm2                         => Atm%flagstruct%trdm2
     d2_bg_k1                      => Atm%flagstruct%d2_bg_k1
     d2_bg_k2                      => Atm%flagstruct%d2_bg_k2
     d2_divg_max_k1                => Atm%flagstruct%d2_divg_max_k1
     d2_divg_max_k2                => Atm%flagstruct%d2_divg_max_k2
     damp_k_k1                     => Atm%flagstruct%damp_k_k1
     damp_k_k2                     => Atm%flagstruct%damp_k_k2
     n_zs_filter                   => Atm%flagstruct%n_zs_filter
     nord_zs_filter                => Atm%flagstruct%nord_zs_filter
     full_zs_filter                => Atm%flagstruct%full_zs_filter
     RF_fast                       => Atm%flagstruct%RF_fast
     consv_am                      => Atm%flagstruct%consv_am
     do_sat_adj                    => Atm%flagstruct%do_sat_adj
     do_f3d                        => Atm%flagstruct%do_f3d
     no_dycore                     => Atm%flagstruct%no_dycore
     convert_ke                    => Atm%flagstruct%convert_ke
     do_vort_damp                  => Atm%flagstruct%do_vort_damp
     use_old_omega                 => Atm%flagstruct%use_old_omega
     beta                          => Atm%flagstruct%beta
     n_sponge                      => Atm%flagstruct%n_sponge
     d_ext                         => Atm%flagstruct%d_ext
     nwat                          => Atm%flagstruct%nwat
     use_logp                      => Atm%flagstruct%use_logp
     warm_start                    => Atm%flagstruct%warm_start
     inline_q                      => Atm%flagstruct%inline_q
     shift_fac                     => Atm%flagstruct%shift_fac
     do_schmidt                    => Atm%flagstruct%do_schmidt
     stretch_fac                   => Atm%flagstruct%stretch_fac
     target_lat                    => Atm%flagstruct%target_lat
     target_lon                    => Atm%flagstruct%target_lon
     regional                      => Atm%flagstruct%regional
     bc_update_interval            => Atm%flagstruct%bc_update_interval 
     reset_eta                     => Atm%flagstruct%reset_eta
     p_fac                         => Atm%flagstruct%p_fac
     a_imp                         => Atm%flagstruct%a_imp
     n_split                       => Atm%flagstruct%n_split
     fac_n_spl                     => Atm%flagstruct%fac_n_spl
     fhouri                        => Atm%flagstruct%fhouri
     m_split                       => Atm%flagstruct%m_split
     k_split                       => Atm%flagstruct%k_split
     use_logp                      => Atm%flagstruct%use_logp
     q_split                       => Atm%flagstruct%q_split
     print_freq                    => Atm%flagstruct%print_freq
     write_3d_diags                => Atm%flagstruct%write_3d_diags
     npx                           => Atm%flagstruct%npx
     npy                           => Atm%flagstruct%npy
     npz                           => Atm%flagstruct%npz
     npz_rst                       => Atm%flagstruct%npz_rst
     ncnst                         => Atm%flagstruct%ncnst
     pnats                         => Atm%flagstruct%pnats
     dnats                         => Atm%flagstruct%dnats
     ntiles                        => Atm%flagstruct%ntiles
     nf_omega                      => Atm%flagstruct%nf_omega
     fv_sg_adj                     => Atm%flagstruct%fv_sg_adj
     na_init                       => Atm%flagstruct%na_init
     nudge_dz                      => Atm%flagstruct%nudge_dz
     p_ref                         => Atm%flagstruct%p_ref
     dry_mass                      => Atm%flagstruct%dry_mass
     nt_prog                       => Atm%flagstruct%nt_prog
     nt_phys                       => Atm%flagstruct%nt_phys
     tau_h2o                       => Atm%flagstruct%tau_h2o
     delt_max                      => Atm%flagstruct%delt_max
     d_con                         => Atm%flagstruct%d_con
     ke_bg                         => Atm%flagstruct%ke_bg
     consv_te                      => Atm%flagstruct%consv_te
     tau                           => Atm%flagstruct%tau
     rf_cutoff                     => Atm%flagstruct%rf_cutoff
     filter_phys                   => Atm%flagstruct%filter_phys
     dwind_2d                      => Atm%flagstruct%dwind_2d
     breed_vortex_inline           => Atm%flagstruct%breed_vortex_inline
     range_warn                    => Atm%flagstruct%range_warn
     fill                          => Atm%flagstruct%fill
     fill_dp                       => Atm%flagstruct%fill_dp
     fill_wz                       => Atm%flagstruct%fill_wz
     check_negative                => Atm%flagstruct%check_negative
     non_ortho                     => Atm%flagstruct%non_ortho
     adiabatic                     => Atm%flagstruct%adiabatic
     moist_phys                    => Atm%flagstruct%moist_phys
     do_Held_Suarez                => Atm%flagstruct%do_Held_Suarez
     do_reed_physics               => Atm%flagstruct%do_reed_physics
     reed_cond_only                => Atm%flagstruct%reed_cond_only
     reproduce_sum                 => Atm%flagstruct%reproduce_sum
     adjust_dry_mass               => Atm%flagstruct%adjust_dry_mass
     fv_debug                      => Atm%flagstruct%fv_debug
     srf_init                      => Atm%flagstruct%srf_init
     mountain                      => Atm%flagstruct%mountain
     remap_t                       => Atm%flagstruct%remap_t
     z_tracer                      => Atm%flagstruct%z_tracer
     old_divg_damp                 => Atm%flagstruct%old_divg_damp
     fv_land                       => Atm%flagstruct%fv_land
     nudge                         => Atm%flagstruct%nudge
     nudge_ic                      => Atm%flagstruct%nudge_ic
     ncep_ic                       => Atm%flagstruct%ncep_ic
     nggps_ic                      => Atm%flagstruct%nggps_ic
     ecmwf_ic                      => Atm%flagstruct%ecmwf_ic
     gfs_phil                      => Atm%flagstruct%gfs_phil
     agrid_vel_rst                 => Atm%flagstruct%agrid_vel_rst
     use_new_ncep                  => Atm%flagstruct%use_new_ncep
     use_ncep_phy                  => Atm%flagstruct%use_ncep_phy
     fv_diag_ic                    => Atm%flagstruct%fv_diag_ic
     external_ic                   => Atm%flagstruct%external_ic
     external_eta                  => Atm%flagstruct%external_eta
     read_increment                => Atm%flagstruct%read_increment

     hydrostatic                   => Atm%flagstruct%hydrostatic
     phys_hydrostatic              => Atm%flagstruct%phys_hydrostatic
     use_hydro_pressure            => Atm%flagstruct%use_hydro_pressure
     do_uni_zfull                  => Atm%flagstruct%do_uni_zfull !miz
     adj_mass_vmr                  => Atm%flagstruct%adj_mass_vmr !f1p
     hybrid_z                      => Atm%flagstruct%hybrid_z
     Make_NH                       => Atm%flagstruct%Make_NH
     make_hybrid_z                 => Atm%flagstruct%make_hybrid_z
     nudge_qv                      => Atm%flagstruct%nudge_qv
     add_noise                     => Atm%flagstruct%add_noise
     a2b_ord                       => Atm%flagstruct%a2b_ord
     c2l_ord                       => Atm%flagstruct%c2l_ord
     ndims                         => Atm%flagstruct%ndims
    
     dx_const                      => Atm%flagstruct%dx_const
     dy_const                      => Atm%flagstruct%dy_const
     deglon_start                  => Atm%flagstruct%deglon_start
     deglon_stop                   => Atm%flagstruct%deglon_stop
     deglat_start                  => Atm%flagstruct%deglat_start
     deglat_stop                   => Atm%flagstruct%deglat_stop

     deglat                        => Atm%flagstruct%deglat

     nested                        => Atm%neststruct%nested
     twowaynest                    => Atm%neststruct%twowaynest
     parent_tile                   => Atm%neststruct%parent_tile
     refinement                    => Atm%neststruct%refinement
     nestbctype                    => Atm%neststruct%nestbctype
     nestupdate                    => Atm%neststruct%nestupdate
     nsponge                       => Atm%neststruct%nsponge
     s_weight                      => Atm%neststruct%s_weight
     ioffset                       => Atm%neststruct%ioffset
     joffset                       => Atm%neststruct%joffset

     layout                        => Atm%layout
     io_layout                     => Atm%io_layout
  end subroutine setup_pointers

       
end module fv_control_mod
