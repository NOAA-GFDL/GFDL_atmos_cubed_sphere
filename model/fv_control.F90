
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
   use fms_io_mod,          only: set_domain
   use mpp_mod,             only: FATAL, mpp_error, mpp_pe, stdlog, &
                                  mpp_npes, mpp_get_current_pelist, &
                                  input_nml_file, get_unit, WARNING, &
                                  read_ascii_file, INPUT_STR_LENGTH
   use mpp_domains_mod,     only: mpp_get_data_domain, mpp_get_compute_domain, mpp_get_tile_id
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
   use fv_mp_mod,           only: mp_start, domain_decomp, mp_assign_gid, global_nest_domain
   use fv_mp_mod,           only: broadcast_domains, mp_barrier, is_master, setup_master, grids_master_procs, tile_fine
   use fv_mp_mod,           only: MAX_NNEST, MAX_NTILE
   use test_cases_mod,      only: read_namelist_test_case_nml
   use fv_timing_mod,       only: timing_on, timing_off, timing_init, timing_prt
   use mpp_domains_mod,     only: domain2D
   use mpp_domains_mod,     only: mpp_define_nest_domains, nest_domain_type, mpp_get_global_domain
   use mpp_domains_mod,     only: mpp_get_C2F_index, mpp_get_F2C_index
   use mpp_domains_mod,     only: CENTER, CORNER, NORTH, EAST, WEST, SOUTH
   use mpp_mod,             only: mpp_send, mpp_sync, mpp_transmit, mpp_set_current_pelist, &
                                  mpp_declare_pelist, mpp_root_pe, mpp_recv, mpp_sync_self, read_input_nml, &
                                  mpp_max
   use fv_diagnostics_mod,  only: fv_diag_init_gn

#ifdef MULTI_GASES
   use constants_mod,       only: rvgas, cp_air
   use multi_gases_mod,     only: multi_gases_init, &
                                  read_namelist_multi_gases_nml
#endif

   implicit none
   private

#ifdef OVERLOAD_R4
   real    :: too_big  = 1.E8
#else
   real    :: too_big  = 1.E35
#endif
   public :: fv_control_init, fv_end

   integer, public :: ngrids = 1
   integer :: commID, global_commID

   integer :: halo_update_type = 1 ! 1 for two-interfaces non-block
                                   ! 2 for block
                                   ! 3 for four-interfaces non-block

! version number of this module
! Include variable "version" to be written to log file.
#include<file_version.h>

 contains

!-------------------------------------------------------------------------------
         
   subroutine fv_control_init(Atm, dt_atmos, this_grid, grids_on_this_pe, p_split)

     type(fv_atmos_type), allocatable, intent(inout), target :: Atm(:)
     real,                intent(in)    :: dt_atmos
     integer,             intent(OUT)   :: this_grid
     logical, allocatable, intent(OUT) :: grids_on_this_pe(:)

     integer, intent(INOUT) :: p_split
     character(100) :: pe_list_name, errstring
     integer :: n, npes, pecounter, i, num_family, ntiles_nest_all
     integer, allocatable :: global_pelist(:)
     integer, dimension(MAX_NNEST) :: grid_pes = 0
     integer, dimension(MAX_NNEST) :: grid_coarse = -1
     integer, dimension(MAX_NNEST) :: nest_refine = 3 
     integer, dimension(MAX_NNEST) :: nest_ioffsets = -999, nest_joffsets = -999
     integer, dimension(MAX_NNEST) :: all_npx = 0
     integer, dimension(MAX_NNEST) :: all_npy = 0
     integer, dimension(MAX_NNEST) :: all_npz = 0
     integer, dimension(MAX_NNEST) :: all_ntiles = 0
     integer, dimension(MAX_NNEST) :: all_twowaynest = 0 ! > 0 implies two-way
     !integer, dimension(MAX_NNEST) :: tile_fine = 0
     integer, dimension(MAX_NNEST) :: icount_coarse = 1
     integer, dimension(MAX_NNEST) :: jcount_coarse = 1
     integer, dimension(MAX_NNEST) :: nest_level = 0
     integer, dimension(MAX_NNEST) :: tile_coarse = 0
     integer, dimension(MAX_NTILE) :: npes_nest_tile = 0

     real :: sdt
     integer :: unit, ens_root_pe, tile_id(1)

     !!!!!!!!!! POINTERS FOR READING NAMELISTS !!!!!!!!!!

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
     logical , pointer :: do_schmidt, do_cube_transform 
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
     character(len=24), pointer :: npz_type
     integer , pointer :: npz_rst 

     integer , pointer :: ncnst 
     integer , pointer :: pnats 
     integer , pointer :: dnats 
     integer , pointer :: dnrts 
     integer , pointer :: ntiles        
     integer , pointer :: nf_omega  
     integer , pointer :: fv_sg_adj 
     real    , pointer :: sg_cutoff

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
     logical , pointer :: fill_gfs
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
     logical , pointer :: butterfly_effect

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
     integer, pointer :: nrows_blend
     logical, pointer :: regional_bcs_from_gsi
     logical, pointer :: write_restart_with_bcs
     integer, pointer :: parent_tile, refinement, nestbctype, nestupdate, nsponge, ioffset, joffset
     real, pointer :: s_weight, update_blend

     integer, pointer :: layout(:), io_layout(:)

     !!!!!!!!!! END POINTERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!

     this_grid = -1 ! default
     call mp_assign_gid
     ens_root_pe = mpp_root_pe()

     ! 1. read nesting namelists
     call read_namelist_nest_nml
     call read_namelist_fv_nest_nml

     ! 2. Set up Atm and PElists

     ngrids = 1
     do n=2,MAX_NNEST
        if (grid_coarse(n) <= 0) then
           exit
        endif
        ngrids = ngrids + 1
     enddo
     allocate(Atm(ngrids))
     npes = mpp_npes() ! now on global pelist

     allocate(global_pelist(npes))
     call mpp_get_current_pelist(global_pelist, commID=global_commID) ! for commID


     allocate(grids_master_procs(ngrids))
     pecounter = 0
     allocate(grids_on_this_pe(ngrids))
     grids_on_this_pe(:) = .false.

     do n=1,ngrids

        if (ngrids == 1 .or. grid_pes(n) == 0) then
           grid_pes(n) = npes - sum(grid_pes)
           if (grid_pes(n) == 0) then
              if ( n > 1 ) then
                 call mpp_error(FATAL, 'Only one zero entry in grid_pes permitted.') 
              else
                 grid_pes(n) = npes
              endif
           endif
        endif

        allocate(Atm(n)%pelist(grid_pes(n)))
        grids_master_procs(n) = pecounter
        do i=1,grid_pes(n)
           if (pecounter >= npes) then
              if (mpp_pe() == 0) then
                 print*, 'ngrids = ', ngrids, ', grid_pes = ', grid_pes(1:ngrids)
              endif
              call mpp_error(FATAL, 'grid_pes assigns more PEs than are available.')
           endif
           Atm(n)%pelist(i) = pecounter + ens_root_pe !TODO PELIST set up by mpp_define_nest_domains???
           pecounter = pecounter + 1
           Atm(n)%npes_this_grid = grid_pes(n)
        enddo
        Atm(n)%grid_number = n

        !TODO: we are required to use PE name for reading INTERNAL namelist
        ! and the actual file name for EXTERNAL namelists. Need to clean up this code
        if (n == 1) then
           pe_list_name = ''
        else
           write(pe_list_name,'(A4, I2.2)') 'nest', n
        endif
        call mpp_declare_pelist(Atm(n)%pelist, pe_list_name)
        !If nest need to re-initialize internal NML
        if (n > 1) then
           Atm(n)%nml_filename = 'input_'//trim(pe_list_name)//'.nml'
        else
           Atm(n)%nml_filename = 'input.nml'
        endif
        if (.not. file_exist(Atm(n)%nml_filename)) then
           call mpp_error(FATAL, "Could not find nested grid namelist "//Atm(n)%nml_filename)
        endif
     enddo

     do n=1,ngrids        
        !ONE grid per pe
        if (ANY(mpp_pe() == Atm(n)%pelist)) then
           if (this_grid > 0) then
              print*, mpp_pe(), this_grid, n
              call mpp_error(FATAL, " Grid assigned to multiple pes")
           endif
           call mpp_set_current_pelist(Atm(n)%pelist)
           call setup_master(Atm(n)%pelist)
           this_grid = n
           grids_on_this_pe(n) = .true.
        endif
        Atm(n)%neststruct%nested = ( grid_coarse(n) > 0 )

        if (Atm(n)%neststruct%nested) then
           if ( grid_coarse(n) > ngrids .or. grid_coarse(n) == n .or. grid_coarse(n) < 1) then
              write(errstring,'(2(A,I3))')  "Could not find parent grid #", grid_coarse(n), ' for grid #', n
              call mpp_error(FATAL, errstring)
           endif
           Atm(n)%parent_grid => Atm(grid_coarse(n))

           Atm(n)%neststruct%ioffset                = nest_ioffsets(n)
           Atm(n)%neststruct%joffset                = nest_joffsets(n)
           Atm(n)%neststruct%parent_tile            = tile_coarse(n)  
           Atm(n)%neststruct%refinement             = nest_refine(n)

        else

           Atm(n)%neststruct%ioffset                = -999
           Atm(n)%neststruct%joffset                = -999   
           Atm(n)%neststruct%parent_tile            = -1      
           Atm(n)%neststruct%refinement             = -1

        endif

     enddo

     if (pecounter /= npes) then
        if (mpp_pe() == 0) then
           print*, 'npes = ', npes, ', grid_pes = ', grid_pes(1:ngrids)
           call mpp_error(FATAL, 'grid_pes in fv_nest_Nml does not assign all of the available PEs')
        endif
     endif
 
     ! 3pre. 
     call timing_init
     call timing_on('TOTAL')
     
     ! 3. Read namelists, do option processing and I/O

     call set_namelist_pointers(Atm(this_grid))
     call fv_diag_init_gn(Atm(this_grid))
#ifdef INTERNAL_FILE_NML
     if (this_grid .gt. 1) then
        write(Atm(this_grid)%nml_filename,'(A4, I2.2)') 'nest', this_grid
        if (.not. file_exist('input_'//trim(Atm(this_grid)%nml_filename)//'.nml')) then
           call mpp_error(FATAL, "Could not find nested grid namelist "//'input_'//trim(Atm(this_grid)%nml_filename)//'.nml')
        endif
     else
        Atm(this_grid)%nml_filename = ''
     endif
     call read_input_nml(Atm(this_grid)%nml_filename) !re-reads into internal namelist
#endif
     call read_namelist_fv_grid_nml
     call read_namelist_fv_core_nml(Atm(this_grid)) ! do options processing here too?
#ifdef MULTI_GASES
     call read_namelist_multi_gases_nml(Atm(this_grid)%nml_filename, &
           Atm(this_grid)%flagstruct%ncnst,  Atm(this_grid)%flagstruct%nwat)
#endif
     call read_namelist_test_case_nml(Atm(this_grid)%nml_filename)
     call mpp_get_current_pelist(Atm(this_grid)%pelist, commID=commID) ! for commID
     call mp_start(commID,halo_update_type)

     ! 4. Set up domains
     !    This should make use of new fv_nest_nml namelists
     !!!! TODO TEMPORARY location for this code
     if (Atm(this_grid)%neststruct%nested) then

        if ( Atm(this_grid)%flagstruct%consv_te > 0.) then
           call mpp_error(FATAL, 'The global energy fixer cannot be used on a nested grid. consv_te must be set to 0.')
        end if

        if (mod(Atm(this_grid)%flagstruct%npx-1 , Atm(this_grid)%neststruct%refinement) /= 0 .or. &
            mod(Atm(this_grid)%flagstruct%npy-1, Atm(this_grid)%neststruct%refinement) /= 0) then
           call mpp_error(FATAL, 'npx or npy not an even refinement of its coarse grid.')
        endif

     endif

     if (Atm(this_grid)%flagstruct%regional) then
        if ( Atm(this_grid)%flagstruct%consv_te > 0.) then
           call mpp_error(FATAL, 'The global energy fixer cannot be used on a regional grid. consv_te must be set to 0.')
        end if
     endif

     !Now only one call to mpp_define_nest_domains for ALL nests
     ! set up nest_level, tile_fine, tile_coarse
     ! need number of tiles, npx, and npy on each grid
     ! need to define a global PElist

     all_ntiles(this_grid) = ntiles
     call mpp_max(all_ntiles, ngrids, global_pelist)

     all_npx(this_grid) = npx
     call mpp_max(all_npx, ngrids, global_pelist)

     all_npy(this_grid) = npy
     call mpp_max(all_npy, ngrids, global_pelist)

     all_npz(this_grid) = npz
     call mpp_max(all_npz, ngrids, global_pelist)

     if (Atm(this_grid)%neststruct%twowaynest) all_twowaynest(this_grid) = 1
     call mpp_max(all_twowaynest, ngrids, global_pelist)
     ntiles_nest_all = 0
     do n=1,ngrids
        if (n/=this_grid) then
           Atm(n)%flagstruct%npx = all_npx(n)
           Atm(n)%flagstruct%npy = all_npy(n)
           Atm(n)%flagstruct%npz = all_npz(n)
           Atm(n)%flagstruct%ntiles = all_ntiles(n)
           Atm(n)%neststruct%twowaynest = (all_twowaynest(n) > 0) ! disabled
        endif
        npes_nest_tile(ntiles_nest_all+1:ntiles_nest_all+all_ntiles(n)) = &
             Atm(n)%npes_this_grid / all_ntiles(n)
        ntiles_nest_all = ntiles_nest_all + all_ntiles(n)

        if (n > 1) then
           tile_fine(n) = all_ntiles(n) + tile_fine(n-1)
           if (tile_coarse(n) < 1) then !set automatically; only works for single tile parents
              tile_coarse(n) = tile_fine(grid_coarse(n)) 
           endif
           icount_coarse(n) = all_npx(n)/nest_refine(n)
           jcount_coarse(n) = all_npy(n)/nest_refine(n)
           nest_level(n) = nest_level(grid_coarse(n)) + 1
        else
           tile_fine(n) = all_ntiles(n)
           nest_level(n) = 0
        endif
     enddo
         
     if (mpp_pe() == 0 .and. ngrids > 1) then
        print*, ' NESTING TREE'
        do n=1,ngrids
           write(*,'(12i4)') n, nest_level(n), nest_ioffsets(n), nest_joffsets(n), icount_coarse(n), jcount_coarse(n), tile_fine(n), tile_coarse(n), nest_refine(n), all_ntiles(n), all_npx(n), all_npy(n)
           write(*,*)
        enddo
        print*, npes_nest_tile(1:ntiles_nest_all)
        print*, ''
     endif

     ! 5. domain_decomp()
     call domain_decomp(Atm(this_grid)%flagstruct%npx,Atm(this_grid)%flagstruct%npy,Atm(this_grid)%flagstruct%ntiles,&
          Atm(this_grid)%flagstruct%grid_type,Atm(this_grid)%neststruct%nested, &
          Atm(this_grid)%layout,Atm(this_grid)%io_layout,Atm(this_grid)%bd,Atm(this_grid)%tile_of_mosaic, &
          Atm(this_grid)%gridstruct%square_domain,Atm(this_grid)%npes_per_tile,Atm(this_grid)%domain, &
          Atm(this_grid)%domain_for_coupler,Atm(this_grid)%num_contact,Atm(this_grid)%pelist)
     call set_domain(Atm(this_grid)%domain)
     call broadcast_domains(Atm,Atm(this_grid)%pelist,size(Atm(this_grid)%pelist))
     do n=1,ngrids
        tile_id = mpp_get_tile_id(Atm(n)%domain)       
        Atm(n)%global_tile = tile_id(1) ! only meaningful locally
        Atm(n)%npes_per_tile = size(Atm(n)%pelist)/Atm(n)%flagstruct%ntiles ! domain decomp doesn't set this globally
     enddo

     ! 6. Set up domain and Atm structure
     call tm_register_tracers (MODEL_ATMOS, Atm(this_grid)%flagstruct%ncnst, Atm(this_grid)%flagstruct%nt_prog, &
          Atm(this_grid)%flagstruct%pnats, num_family)
     if(is_master()) then
        write(*,*) 'ncnst=', ncnst,' num_prog=',Atm(this_grid)%flagstruct%nt_prog,' pnats=',Atm(this_grid)%flagstruct%pnats,' dnats=',dnats,&
             ' num_family=',num_family         
        print*, ''
     endif
     if (dnrts < 0) dnrts = dnats

     do n=1,ngrids
        !FIXME still setting up dummy structures for other grids for convenience reasons
        !isc, etc. set in domain_decomp
        call allocate_fv_atmos_type(Atm(n), &
             Atm(n)%bd%isd, Atm(n)%bd%ied, &
             Atm(n)%bd%jsd, Atm(n)%bd%jed, &
             Atm(n)%bd%isc, Atm(n)%bd%iec, &
             Atm(n)%bd%jsc, Atm(n)%bd%jec, &
             Atm(n)%flagstruct%npx,    Atm(n)%flagstruct%npy,   Atm(n)%flagstruct%npz, &
             Atm(n)%flagstruct%ndims,  Atm(n)%flagstruct%ncnst, Atm(n)%flagstruct%ncnst-Atm(n)%flagstruct%pnats, &
             n/=this_grid, n==this_grid, ngrids) !TODO don't need both of the last arguments
     enddo
     if ( (Atm(this_grid)%bd%iec-Atm(this_grid)%bd%isc+1).lt.4 .or. (Atm(this_grid)%bd%jec-Atm(this_grid)%bd%jsc+1).lt.4 ) then
        if (is_master()) write(*,'(6I6)') Atm(this_grid)%bd%isc, Atm(this_grid)%bd%iec, Atm(this_grid)%bd%jsc, Atm(this_grid)%bd%jec, this_grid
        call mpp_error(FATAL,'Domain Decomposition:  Cubed Sphere compute domain has a &
             &minium requirement of 4 points in X and Y, respectively')
     end if


     !Tile_coarse is needed to determine which processors are needed to send around their
     ! data for computing the interpolation coefficients
     if (ngrids > 1) then
        !reset to universal pelist
        call mpp_set_current_pelist( global_pelist )
        !Except for npes_nest_tile all arrays should be just the nests and should NOT include the top level
        call mpp_define_nest_domains(global_nest_domain, Atm(this_grid)%domain, &
             ngrids-1, nest_level=nest_level(2:ngrids) , &
             istart_coarse=nest_ioffsets(2:ngrids), jstart_coarse=nest_joffsets(2:ngrids), &
             icount_coarse=icount_coarse(2:ngrids), jcount_coarse=jcount_coarse(2:ngrids), &
             npes_nest_tile=npes_nest_tile(1:ntiles_nest_all), &
             tile_fine=tile_fine(2:ngrids), tile_coarse=tile_coarse(2:ngrids), &
             x_refine=nest_refine(2:ngrids), y_refine=nest_refine(2:ngrids), name="global_nest_domain")
        call mpp_set_current_pelist(Atm(this_grid)%pelist)

     endif

     allocate(Atm(this_grid)%neststruct%child_grids(ngrids))
     do n=1,ngrids
        Atm(this_grid)%neststruct%child_grids(n) = (grid_coarse(n) == this_grid)
        allocate(Atm(n)%neststruct%do_remap_bc(ngrids))
        Atm(n)%neststruct%do_remap_bc(:) = .false.
     enddo
     Atm(this_grid)%neststruct%parent_proc = ANY(Atm(this_grid)%neststruct%child_grids) !ANY(tile_coarse == Atm(this_grid)%global_tile)
     Atm(this_grid)%neststruct%child_proc = ASSOCIATED(Atm(this_grid)%parent_grid) !this means a nested grid

     if (ngrids > 1) call setup_update_regions
     if (Atm(this_grid)%neststruct%nestbctype > 1) then
        call mpp_error(FATAL, 'nestbctype > 1 not yet implemented')
        Atm(this_grid)%neststruct%upoff = 0
     endif

     if (Atm(this_grid)%gridstruct%bounded_domain .and. is_master()) print*, &
          ' Bounded domain: nested = ', Atm(this_grid)%neststruct%nested, ', regional = ', Atm(this_grid)%flagstruct%regional

     ! 7. Init_grid() (including two-way nesting)
     call init_grid(Atm(this_grid), Atm(this_grid)%flagstruct%grid_name, Atm(this_grid)%flagstruct%grid_file, &
          Atm(this_grid)%flagstruct%npx, Atm(this_grid)%flagstruct%npy, Atm(this_grid)%flagstruct%npz, Atm(this_grid)%flagstruct%ndims, Atm(this_grid)%flagstruct%ntiles, Atm(this_grid)%ng, tile_coarse)


     ! 8. grid_utils_init()
     ! Initialize the SW (2D) part of the model
     call grid_utils_init(Atm(this_grid), Atm(this_grid)%flagstruct%npx, Atm(this_grid)%flagstruct%npy, Atm(this_grid)%flagstruct%npz, Atm(this_grid)%flagstruct%non_ortho, Atm(this_grid)%flagstruct%grid_type, Atm(this_grid)%flagstruct%c2l_ord)

     ! Finish up initialization; write damping coefficients dependent upon 

     if ( is_master() ) then
        sdt =  dt_atmos/real(Atm(this_grid)%flagstruct%n_split*Atm(this_grid)%flagstruct%k_split*abs(p_split))
        write(*,*) ' '
        write(*,*) 'Divergence damping Coefficients'
        write(*,*) 'For small dt=', sdt
        write(*,*) 'External mode del-2 (m**2/s)=',  Atm(this_grid)%flagstruct%d_ext*Atm(this_grid)%gridstruct%da_min_c/sdt
        write(*,*) 'Internal mode del-2 SMAG dimensionless coeff=',  Atm(this_grid)%flagstruct%dddmp
        write(*,*) 'Internal mode del-2 background diff=', Atm(this_grid)%flagstruct%d2_bg*Atm(this_grid)%gridstruct%da_min_c/sdt

        if (nord==1) then
           write(*,*) 'Internal mode del-4 background diff=', Atm(this_grid)%flagstruct%d4_bg
           write(*,*) 'Vorticity del-4 (m**4/s)=', (Atm(this_grid)%flagstruct%vtdm4*Atm(this_grid)%gridstruct%da_min)**2/sdt*1.E-6
        endif
        if (Atm(this_grid)%flagstruct%nord==2) write(*,*) 'Internal mode del-6 background diff=', Atm(this_grid)%flagstruct%d4_bg
        if (Atm(this_grid)%flagstruct%nord==3) write(*,*) 'Internal mode del-8 background diff=', Atm(this_grid)%flagstruct%d4_bg
        write(*,*) 'tracer del-2 diff=', Atm(this_grid)%flagstruct%trdm2

        write(*,*) 'Vorticity del-4 (m**4/s)=', (Atm(this_grid)%flagstruct%vtdm4*Atm(this_grid)%gridstruct%da_min)**2/sdt*1.E-6
        write(*,*) 'beta=', Atm(this_grid)%flagstruct%beta
        write(*,*) ' '
     endif


!!$     Atm(this_grid)%ts   = 300.
!!$     Atm(this_grid)%phis = too_big
!!$     ! The following statements are to prevent the phantom corner regions from
!!$     ! growing instability
!!$     Atm(this_grid)%u  = 0.
!!$     Atm(this_grid)%v  = 0.
!!$     Atm(this_grid)%ua = too_big
!!$     Atm(this_grid)%va = too_big
!!$

     !Initialize restart
     call fv_restart_init()
!     if ( reset_eta ) then
!         do n=1, ntilesMe
!            call set_eta(npz, Atm(this_grid)%ks, ptop, Atm(this_grid)%ak, Atm(this_grid)%bk, Atm(this_grid)%flagstruct%npz_type)
!         enddo
!         if(is_master()) write(*,*) "Hybrid sigma-p coordinate has been reset"
!     endif


   contains
!>@brief The subroutine 'setup_namelist_pointers' associates the MODULE flag pointers
!! with the ARRAY flag variables for the grid active on THIS pe so the flags
!! can be read in from the namelist.
     subroutine set_namelist_pointers(Atm)
       type(fv_atmos_type), intent(INOUT), target :: Atm

       !This routine associates the MODULE flag pointers with the ARRAY flag variables for the grid active on THIS pe so the flags can be read in from the namelist.

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
       do_cube_transform             => Atm%flagstruct%do_cube_transform
       stretch_fac                   => Atm%flagstruct%stretch_fac
       target_lat                    => Atm%flagstruct%target_lat
       target_lon                    => Atm%flagstruct%target_lon
       regional                      => Atm%flagstruct%regional
       bc_update_interval            => Atm%flagstruct%bc_update_interval 
       nrows_blend                   => Atm%flagstruct%nrows_blend
       regional_bcs_from_gsi         => Atm%flagstruct%regional_bcs_from_gsi
       write_restart_with_bcs        => Atm%flagstruct%write_restart_with_bcs
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
       npz_type                      => Atm%flagstruct%npz_type
       npz_rst                       => Atm%flagstruct%npz_rst
       ncnst                         => Atm%flagstruct%ncnst
       pnats                         => Atm%flagstruct%pnats
       dnats                         => Atm%flagstruct%dnats
       dnrts                         => Atm%flagstruct%dnrts
       ntiles                        => Atm%flagstruct%ntiles
       nf_omega                      => Atm%flagstruct%nf_omega
       fv_sg_adj                     => Atm%flagstruct%fv_sg_adj
       sg_cutoff                     => Atm%flagstruct%sg_cutoff
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
       fill_gfs                      => Atm%flagstruct%fill_gfs
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
       butterfly_effect              => Atm%flagstruct%butterfly_effect
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
       update_blend                  => Atm%neststruct%update_blend

       layout                        => Atm%layout
       io_layout                     => Atm%io_layout
     end subroutine set_namelist_pointers


     subroutine read_namelist_nest_nml

       integer :: f_unit, ios, ierr, dum
       namelist /nest_nml/ dum ! ngrids, ntiles, nest_pes, p_split !emptied lmh 7may2019

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
       if (ierr > 0) then
          call mpp_error(FATAL, " &nest_nml is depreciated. Please use &fv_nest_nml instead.")
       endif

     end subroutine read_namelist_nest_nml

     subroutine read_namelist_fv_nest_nml

       integer :: f_unit, ios, ierr
       namelist /fv_nest_nml/ grid_pes, grid_coarse, tile_coarse, nest_refine, &
            nest_ioffsets, nest_joffsets, p_split

#ifdef INTERNAL_FILE_NML
       read (input_nml_file,fv_nest_nml,iostat=ios)
       ierr = check_nml_error(ios,'fv_nest_nml')
#else
       f_unit=open_namelist_file()
       rewind (f_unit)
       read (f_unit,fv_nest_nml,iostat=ios)
       ierr = check_nml_error(ios,'fv_nest_nml')
       call close_file(f_unit)
#endif

     end subroutine read_namelist_fv_nest_nml

     subroutine read_namelist_fv_grid_nml

       integer :: f_unit, ios, ierr
       !  local version of these variables to allow PGI compiler to compile
       character(len=80)  :: grid_name = ''
       character(len=120) :: grid_file = ''
       namelist /fv_grid_nml/ grid_name, grid_file

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
       rewind (f_unit)
#endif
       call write_version_number ( 'FV_CONTROL_MOD', version )
       unit = stdlog()
       write(unit, nml=fv_grid_nml)

       !Basic option processing
       if (len_trim(grid_file) /= 0) Atm(this_grid)%flagstruct%grid_file = grid_file
       if (len_trim(grid_name) /= 0) Atm(this_grid)%flagstruct%grid_name = grid_name


     end subroutine read_namelist_fv_grid_nml

     subroutine read_namelist_fv_core_nml(Atm)

       type(fv_atmos_type), intent(inout) :: Atm
       integer :: f_unit, ios, ierr
       real :: dim0 = 180.           ! base dimension
       real :: dt0  = 1800.          ! base time step
       real :: ns0  = 5.             ! base nsplit for base dimension 
       real :: dimx, dl, dp, dxmin, dymin, d_fac
       real :: umax = 350.           ! max wave speed for grid_type>3

       integer :: n0split

       !  local version of these variables to allow PGI compiler to compile
       character(len=128) :: res_latlon_dynamics = ''
       character(len=128) :: res_latlon_tracers  = ''

       namelist /fv_core_nml/npx, npy, ntiles, npz, npz_type, npz_rst, layout, io_layout, ncnst, nwat,  &
            use_logp, p_fac, a_imp, k_split, n_split, m_split, q_split, print_freq, write_3d_diags, &
            do_schmidt, do_cube_transform, &
            hord_mt, hord_vt, hord_tm, hord_dp, hord_tr, shift_fac, stretch_fac, target_lat, target_lon, &
            kord_mt, kord_wz, kord_tm, kord_tr, fv_debug, fv_land, nudge, do_sat_adj, do_f3d, &
            external_ic, read_increment, ncep_ic, nggps_ic, ecmwf_ic, use_new_ncep, use_ncep_phy, fv_diag_ic, &
            external_eta, res_latlon_dynamics, res_latlon_tracers, scale_z, w_max, z_min, lim_fac, &
            dddmp, d2_bg, d4_bg, vtdm4, trdm2, d_ext, delt_max, beta, non_ortho, n_sponge, &
            warm_start, adjust_dry_mass, mountain, d_con, ke_bg, nord, nord_tr, convert_ke, use_old_omega, &
            dry_mass, grid_type, do_Held_Suarez, do_reed_physics, reed_cond_only, &
            consv_te, fill, filter_phys, fill_dp, fill_wz, fill_gfs, consv_am, RF_fast, &
            range_warn, dwind_2d, inline_q, z_tracer, reproduce_sum, adiabatic, do_vort_damp, no_dycore,   &
            tau, tau_h2o, rf_cutoff, nf_omega, hydrostatic, fv_sg_adj, sg_cutoff, breed_vortex_inline,  &
            na_init, nudge_dz, hybrid_z, Make_NH, n_zs_filter, nord_zs_filter, full_zs_filter, reset_eta,         &
            pnats, dnats, dnrts, a2b_ord, remap_t, p_ref, d2_bg_k1, d2_bg_k2,  &
            c2l_ord, dx_const, dy_const, umax, deglat,      &
            deglon_start, deglon_stop, deglat_start, deglat_stop, &
            phys_hydrostatic, use_hydro_pressure, make_hybrid_z, old_divg_damp, add_noise, butterfly_effect, &
            nested, twowaynest, nudge_qv, &
            nestbctype, nestupdate, nsponge, s_weight, &
            check_negative, nudge_ic, halo_update_type, gfs_phil, agrid_vel_rst,     &
            do_uni_zfull, adj_mass_vmr, fac_n_spl, fhouri, update_blend, regional, bc_update_interval,  &
            regional_bcs_from_gsi, write_restart_with_bcs, nrows_blend

#ifdef INTERNAL_FILE_NML
       ! Read FVCORE namelist 
       read (input_nml_file,fv_core_nml,iostat=ios)
       ierr = check_nml_error(ios,'fv_core_nml')
       ! Reset input_file_nml to default behavior (CHECK do we still need this???)
       !call read_input_nml
#else
       f_unit = open_namelist_file(Atm%nml_filename)
       ! Read FVCORE namelist 
       read (f_unit,fv_core_nml,iostat=ios)
       ierr = check_nml_error(ios,'fv_core_nml')
       call close_file(f_unit)
#endif         
       call write_version_number ( 'FV_CONTROL_MOD', version )
       unit = stdlog()
       write(unit, nml=fv_core_nml)

       if (len_trim(res_latlon_dynamics) /= 0) Atm%flagstruct%res_latlon_dynamics = res_latlon_dynamics
       if (len_trim(res_latlon_tracers)  /= 0) Atm%flagstruct%res_latlon_tracers = res_latlon_tracers

       !*** single tile for Cartesian grids
       if (grid_type>3) then
          ntiles=1
          non_ortho = .false.
          nf_omega = 0
       endif

       if (.not. (nested .or. regional)) Atm%neststruct%npx_global = npx

       ! Define n_split if not in namelist
       if (ntiles==6) then
          dimx = 4.0*(npx-1)
          if ( hydrostatic ) then
             if ( npx >= 120 ) ns0 = 6
          else
             if ( npx <= 45 ) then
                ns0 = 6
             elseif ( npx <=90 ) then
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

          dxmin=dl*radius*min(cos(deglat_start*pi/180.-Atm%bd%ng*dp),   &
               cos(deglat_stop *pi/180.+Atm%bd%ng*dp))
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

       if (old_divg_damp) then
          if (is_master()) write(*,*) " fv_control: using AM2/AM3 damping methods "
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

197    format(A,l7)
198    format(A,i2.2,A,i4.4,'x',i4.4,'x',i1.1,'-',f9.3)
199    format(A,i3.3)

       !if (.not. (nested .or. regional)) alpha = alpha*pi  !TODO for test_case_nml

       !allocate(Atm%neststruct%child_grids(size(Atm))) !TODO want to remove
       !Atm(N)%neststruct%child_grids = .false.

       target_lon = target_lon * pi/180.
       target_lat = target_lat * pi/180.

     end subroutine read_namelist_fv_core_nml

     subroutine setup_update_regions

       integer :: isu, ieu, jsu, jeu ! update regions
       integer :: isc, jsc, iec, jec
       integer :: upoff

       isc = Atm(this_grid)%bd%isc
       jsc = Atm(this_grid)%bd%jsc
       iec = Atm(this_grid)%bd%iec
       jec = Atm(this_grid)%bd%jec

       upoff = Atm(this_grid)%neststruct%upoff

       do n=2,ngrids
          if (tile_coarse(n) == Atm(this_grid)%global_tile) then

             isu = nest_ioffsets(n)
             ieu = isu + icount_coarse(n) - 1
             jsu = nest_joffsets(n)
             jeu = jsu + jcount_coarse(n) - 1

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

     end subroutine setup_update_regions

   end subroutine fv_control_init

!-------------------------------------------------------------------------------
         
 !>@brief The subroutine 'fv_end' terminates FV3, deallocates memory, 
!! saves restart files, and stops I/O.
 subroutine fv_end(Atm, this_grid, restart_endfcst)

    type(fv_atmos_type), intent(inout) :: Atm(:)
    integer, intent(IN) :: this_grid
    logical, intent(in) :: restart_endfcst

    integer :: n

    call timing_off('TOTAL')
    call timing_prt( mpp_pe() )

    call fv_restart_end(Atm(this_grid), restart_endfcst)
    call fv_io_exit()

  ! Free temporary memory from sw_core routines
  ! Deallocate
    call grid_utils_end

    do n = 1, ngrids
       call deallocate_fv_atmos_type(Atm(n))
    end do


 end subroutine fv_end
!-------------------------------------------------------------------------------

end module fv_control_mod
