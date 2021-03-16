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
   use coarse_grained_restart_files_mod, only: deallocate_coarse_restart_type

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
     integer :: n, npes, pecounter, i, num_family, ntiles_nest_all, num_tile_top
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
     logical , pointer :: do_inline_mp
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
   logical , pointer :: hrrrv3_ic
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

     character(len=16), pointer :: restart_resolution
     integer, pointer :: layout(:), io_layout(:)
     logical, pointer :: write_coarse_restart_files
     logical, pointer :: write_coarse_diagnostics
     logical, pointer :: write_only_coarse_intermediate_restarts
     logical, pointer :: write_coarse_agrid_vel_rst
     logical, pointer :: write_coarse_dgrid_vel_rst
     !!!!!!!!!! END POINTERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!

     this_grid = -1 ! default
     call mp_assign_gid
     ens_root_pe = mpp_root_pe()

     ! 1. read nesting namelists
     call read_namelist_nest_nml
     call read_namelist_fv_nest_nml

     ! 2. Set up Atm and PElists
     do n=2,MAX_NNEST
        if (tile_coarse(n) > 0) then
           if (tile_coarse(n)<=num_tile_top) then
              grid_coarse(n)=1
           else
              grid_coarse(n)=tile_coarse(n) - num_tile_top + 1
           endif
        endif
     enddo

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
           Atm(n)%neststruct%nlevel=nest_level(n)
           if (n==ngrids) Atm(:)%neststruct%num_nest_level=nest_level(ngrids)
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
     call domain_decomp(this_grid,Atm(this_grid)%flagstruct%npx,Atm(this_grid)%flagstruct%npy,Atm(this_grid)%flagstruct%ntiles,&
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
!!$     Atm(this_grid)%inline_mp%prer = too_big
!!$     Atm(this_grid)%inline_mp%prei = too_big
!!$     Atm(this_grid)%inline_mp%pres = too_big
!!$     Atm(this_grid)%inline_mp%preg = too_big

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
       do_inline_mp                  => Atm%flagstruct%do_inline_mp
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
       hrrrv3_ic                     => Atm%flagstruct%hrrrv3_ic
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

       write_coarse_restart_files    => Atm%coarse_graining%write_coarse_restart_files
       write_coarse_diagnostics      => Atm%coarse_graining%write_coarse_diagnostics
       write_only_coarse_intermediate_restarts => Atm%coarse_graining%write_only_coarse_intermediate_restarts
       write_coarse_agrid_vel_rst    => Atm%coarse_graining%write_coarse_agrid_vel_rst
       write_coarse_dgrid_vel_rst    => Atm%coarse_graining%write_coarse_dgrid_vel_rst
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
       namelist /fv_nest_nml/ grid_pes, num_tile_top, tile_coarse, nest_refine, &
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
       call close_file (f_unit)
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
!> \param[in] butterfly\_effect Logical: whether to flip the least-significant-bit of the lowest level temperature. False by default.
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
!> \param[in] npz\_type Character(24): controls which level setup to use when several vertical level configurations use the same number of levels. These are defined in fv\_eta.F90; check there for details for defaults (ie. empty npz\_type) and alternates. Empty by default.
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
!> \param[in] read_increment  Logical: Read in analysis increment and add to restart following are namelist parameters for Stochastic Energy Baskscatter dissipation estimate. This is useful as part of a data-assimilation cycling system or to use native restarts from the six-tile first guess, after which the analysis increment can be applied.
!!
!> \param[in] res\_latlon\_dynamics character(len=128) If external\_ic =.true. gives the filename of the input IC file. INPUT/fv\_rst.res.nc by default.
!!
!> \param[in] res\_latlon\_tracers character(len=128) If external\_ic =.true. and both ncep\_ic and fv\_diag\_ic are.false., this variable gives the filename of the initial conditions for the tracers, assumed to be a legacy lat-lon FV core restart file. INPUT/atmos\_tracers.res.nc by default.
!!
!> \param[in] warm\_start Logical: whether to start from restart files, instead of cold-starting the model. True by default; if this is set to true and restart files cannot be found the model will stop.
!!
!>###A1.3 I/O and diagnostic options:
!!
!> \param[in] agrid\_vel\_rst Logical: whether to write the unstaggered latitude-longitude winds (u<sub>a</sub> and v<sub>a</sub>) to the restart files. This is useful for data assimilation cycling systems which do not handle staggered winds. .false. by default.
!!
!> \param[in] bc_update_interval Integer: Default setting for interval (hours) between external regional BC data files.
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
!> \param[in] do_sat_adj  Logical: The same as fast_sat_adj = .false.  has fast saturation adjustments
!!
!> \param[in] dnats Integer: The number of tracers which are not to be advected by the dynamical core, but still passed into the dynamical core; the last dnats+pnats tracers in field\_table are not advected. 0 by default.
!!
!> \param[in] dnrts Integer:  the Number of non-remapped consituents. Only makes sense for dnrts <= dnat.
!!
!> \param[in] dwind\_2d Logical: whether to use a simpler \& faster algorithm for interpolating the A-grid (cell-centered) wind tendencies computed from the physics to the D-grid. Typically, the A-grid wind tendencies are first converted in 3D cartesian coordinates and then interpolated before converting back to 2D local coordinates. When this option enabled, a much simpler but less accurate 2D interpolation is used. False by default.
!!
!> \param[in] fill Logical: Fills in negative tracer values by taking positive tracers from the cells above and below. This option is useful when the physical parameterizations produced negatives. False by default.
!!
!> \param[in] gfs_phil Logical: Obsolete - to be removed
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
!> \param[in] do\_cube\_transform logical: applies same transform as when do\_schmidt = True but rotates from the north pole instead of the south pole. If both do\_cube\_transform and do\_schmidt are True the model throws a fatal error. False by default.
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
!> \param[in] nestupdate  Integer: type of nested-grid update to use; details are given in model/fv\_nesting.F90. 0 by default.
!!
!> \param[in] regional Logical: Controls whether this is a regional domain (and thereby needs external BC inputs)
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
!> \param[in] sg\_cutoff Real (in Pa): the pressure above which the 2dz filter in fv_sg is applied, similar to the behavior of rf_cutoff. If this value is set to a non-negative value it overrides the value in n\_sponge. -1 by default, which disables this option and uses n_sponge instead.
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
!> \param[in] d2\_bg\_k2  Real: strength of second-order diffusion in the second sponge layer from the model top. 0.02 by default. This value should be lower than d2\_bg\_k1. If d2\_bg\_k2=0, then d2\_bg\_k1 is applied throughout the depth of the sponge layer (the bottom of the sponge layer is set by rf_cutoff). The amplitude is d2\_bg\_k1 at the top, then decreases downward with the same vertical dependence as the rayleigh damping, going to zero at rf_cutoff.
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
!>###A.1.10 Limited area model (LAM)
!!
!> \param[in] update\_blend Real: Weights to control how much blending is done during two-way nesting update. Default is 1.
!!
!> \param[in] regional\_bcs\_from\_gsi Logical: whether DA-updated BC files are used. Default is false.
!!
!> \param[in] write\_restart\_with\_bcs Logical: whether to write restart files with BC rows.
!!
!> \param[in] nrows\_blend Integer: Number of blending rows in the files.
!!
!>##A.2 Entries in external\_ic\_nml
!!
!> \param[in] filtered\_terrain  Logical: whether to use the terrain filtered by the preprocessing tools rather than the raw terrain. .true. by default. Only active if nggps\_ic = .true. or ecmwf\_ic = .true.
!!
!> \param[in] levp  Integer: number of levels in the input (remapped) initial conditions. 64 by default. Only active if nggps\_ic = .true.
!!
!> \param[in] gfs_dwinds  Logical: obsolete - to be removed
!!
!> \param[in] checker\_tr  Logical: whether to enable the ``checkerboard'' tracer test. .false. by default. Only active if nggps\_ic = .true.
!!
!> \param[in] nt\_checker  Integer: number of tracers (at the end of the list of tracers defined in field\_table) to initialize with an idealized ``checkerboard'' pattern, with values of either 0 or 1. This is intended to test the monotonicity or positivity constraint in the advection scheme. 0 by default. Only active if nggps\_ic = .true.
!!
!>##A.3 Entries in surf\_map\_nml
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
!>##A.4 Entries in fv\_grid\_nml
!!
!> \param[in] grid\_name Character(len=80): Name of the grid either being read in (if grid\_spec = -1) or being created. This is only used for writing out a binary file in the directory from which the model is run. Gnomonic by default.
!!
!> \param[in] grid\_file  Character(len=120): If grid\_type = -1 the name of the grid\_spec file to read in. INPUT/grid\_spec.nc by default; other values will not work.
!!
!>##A.5 Entries in test\_case\_nml
!!
!> \param[in] test\_case  Integer: number of the idealized test case to run. A number of nest cases are defined in tools/test\_cases.F90, of which numbers 19 are intended for the shallow-water model. Requires warm\_start =.false. 11 by default; this creates a resting atmosphere with a very basic thermodynamic profile, with topography. If you wish to initialize an Aquaplanet simulation (no topography) set to 14.
!!
!> \param[in] alpha  Real: In certain shallow-water test cases specifies the angle (in fractions of a rotation, so 0.25 is a 45-degree rotation) through which the basic state is rotated. 0 by default.
!!
!>##A.6  Entries in fv\_nest\_nml
!!
!> \param[in] grid\_pes Integer(:): Number of processor cores (or MPI ranks) assigned to each grid. The sum of the assigned cores in this array must sum to the number of cores allocated to the model. Up to one of the first ngrids entries may be 0, in which case all remaining cores are assigned to it. 0 by default.
!!
!> \param[in] grid_coarse Integer(:): Grid number of parent grid, if any. The first element is ignored; positive values in any successive element indicates that a new nested grid is to be created; the model continues to create grids until it finds a non-positive element. The total number of grids, ngrids, is determined from this array. -1 by default.
!!
!> \param[in] tile\_coarse Integer(:): Absolute index of the tile (sub-component of a grid) a given grid is nested within. The first element is ignored. This is only useful when the parent grid has multiple tiles, such as for the six-tile cubed-sphere grid, or if the parent is a multi-tile nest. If grid_coarse(n) is a single-tile domain the coarse-grid tile can be determined automatically by setting it to a non-positive value. 0 by default.
!!
!> \param[in] nest\_refine Integer(:): Refinement ratio, relative to the parent, for each nest. The first element is ignored. FV3 supports any integer refinement ratio for nesting, although values larger than 6 may cause stability issues. 3 by default.
!!
!> \param[in] nest\_ioffsets integer(:): Index within a tile or grid of the coarse grid cell closest to the local left-hand boundary with a nested grid cell within it. For example, if you have a coarse grid which is 12 grid cells long, and nest\_ioffset for its child is 3, then the nested grid will have its lower-left corner located in the i=3 grid cell.
!!
!> \param[in] nest\_joffsets Integer(:): as for nest\_ioffsets but in the local y-direction.
!!
!>##A.7 Entries in fv\_diag\_column\_nml
!!
!> \param[in] do\_diag\_sonde Logical: whether to enable sounding output specified by the namelist variables diag_sonde* . The output is intended to match the format of text files produced by the University of Wyoming text soundings, except that output is on uninterpolated model levels. False by default.
!!
!> \param[in] sound\_freq integer: Frequency, in hours, of diag_sonde column output. 3 by default.
!!
!> \param[in] runname Character(100): Name of the simulation. This is only to add the runname to the output to enable a user to easily determine the simulation that produced a certain sounding file.
!!
!> \param[in] diag\_sonde\_lon\_in Real (MAX_DIAG_COLUMN): List of longitudes (in degrees) for the desired sounding points. Longitudes may be used as [0, 360] (GFDL style) or [-180,180] (NCAR style).
!!
!> \param[in] diag\_sonde\_lat\_in Real(MAX_DIAG_COLUMN): As for diag_sonde_lon_in except for latitudes.
!!
!> \param[in] diag\_sonde\_names Character(MAX_DIAG_COLUMN): List of names for each sounding point.
!!
!> \param[in] do\_diag\_debug Logical: analogous to the functionality of do\_diag\_sonde, as well as including similar parameters: diag\_debug\_lon\_in, diag\_debug\_lat\_in, and diag\_debug\_names, but outputs different diagnostics at every dt_atmos more appropriate for debugging problems that are known to occur at a specific point in the model.  This functionality is only implemented for the nonhydrostatic solver
!!
!>##A.8 Entries in fms\_nml
!!
!> \param[in] domains\_stack\_size  Integer: size (in bytes) of memory array reserved for domains. For large grids or reduced processor counts this can be large (>10 M); if it is not large enough the model will stop and print a recommended value of the stack size. Default is 0., reverting to the default set in MPP (which is probably not large enough for modern applications).
!!
!!
!!
!> @{
       namelist /fv_core_nml/npx, npy, ntiles, npz, npz_type, npz_rst, layout, io_layout, ncnst, nwat,  &
            use_logp, p_fac, a_imp, k_split, n_split, m_split, q_split, print_freq, write_3d_diags, &
            do_schmidt, do_cube_transform, &
            hord_mt, hord_vt, hord_tm, hord_dp, hord_tr, shift_fac, stretch_fac, target_lat, target_lon, &
            kord_mt, kord_wz, kord_tm, kord_tr, fv_debug, fv_land, nudge, do_sat_adj, do_inline_mp, do_f3d, &
            external_ic, read_increment, ncep_ic, nggps_ic, hrrrv3_ic, ecmwf_ic, use_new_ncep, use_ncep_phy, fv_diag_ic, &
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
            regional_bcs_from_gsi, write_restart_with_bcs, nrows_blend,  &
            write_coarse_restart_files,&
            write_coarse_diagnostics,&
            write_only_coarse_intermediate_restarts, &
            write_coarse_agrid_vel_rst, write_coarse_dgrid_vel_rst


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
          !write(*,'(I, A, 4I)') mpp_pe(), 'SETUP_UPDATE_REGIONS 0: ', mpp_pe(), tile_coarse(n), Atm(this_grid)%global_tile
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
       call deallocate_coarse_restart_type(Atm(n)%coarse_graining%restart)
    end do


 end subroutine fv_end
!-------------------------------------------------------------------------------


end module fv_control_mod
