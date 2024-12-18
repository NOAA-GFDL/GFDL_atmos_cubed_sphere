# RELEASE NOTES for FV3 202411: Summary
FV3-202411-public --- November 2024  
Primary Point of Contact: Lucas Harris, GFDL lucas.harris@noaa.gov  

This version has been tested with:  
SHiELD physics release FV3-202411-public from https://github.com/NOAA-GFDL/SHiELD_physics  
FMS release 2024.03 from https://github.com/NOAA-GFDL/FMS  
FMS Coupler release 2024.03.01 from https://github.com/NOAA-GFDL/FMScoupler  
Atmos Drivers release FV3-202411-public from https://github.com/NOAA-GFDL/atmos_drivers  

This release includes the following:
- Numerics updates (Lucas, Joseph, Linjiong):
  - Removed `USE_COND` and `MOIST_CAPPA` compiler directives and replaced with runtime options `use_cond` and `moist_kappa` in the new namelist `fv_thermo_nml`. Both default to `.T.` in nonhydrostatic simulation and `.F.` in hydrostatic.
  - Added a simple limiter to prevent dissipative heating from creating spurious cooling. Set `prevent_diss_cooling = .F.` to turn off.
  - Fixed bugs in hydrostatic nesting for west and east BCs in `setup_pt_BC_k()` and calculated true pressure in BCs if no BC remapping is done (compute_peBC, compute_peBC_k).
  - Revision of `Lagrangian_to_Eulerian` to fix variable dimension mismatch.
  - Revision of FV3's dissipation heating `diss_est` option to improve numerical consistency with other dissipation options.
  - Fixed edge noise for hord 6 and 7 (suggested by Bill Putman, GMAO).
  - Add mixed precision compilation mode to support 32bit FV3 with other 64bit components (with Uriel Ramirez).
  - New tracers:
    - `w_diff` to allow subgrid mixing of vertical velocity by physics. This requires compiling with the option `-DW_DIFF` to enable.
    - `pbl_age` and `tro_pbl_age` tracers representing the age of air since leaving the PBL and tropical PBL, respectively. 
    - Removed obsolete clock tracers
    - Refer to `docs/HOWTO_tracer-2024.11.md` for more information
- GFDL Microphysics updates (Linjiong)
  - Included fast microphysics tendencies diagnostics 
  - Added two namelist options (`fast_fr_mlt` and `fast_dep_sub`) to control freezing/melting and deposition/sublimation in the fast microphysics.
  - Included a missing term in the energy conservation formula (credit: Tristan Abbott).  May affect  prediction of processes depending strongly on microphysics. Compile the model with `-DENG_CNV_OLD` to revert this change.
  - Added a namelist option, `prog_cin`, to define the source of CIN (cloud ice nuclei) concentration. This is similar to `prog_ccn` but for ice nuclei.
  - Added diagnostics for cloud content and cloud effective radii of all cloud hydrometeors (qc*, re*).
  - Added diagnostics for microphysical process rates (mpp*).
  - Removed unused Keihl et al. (1994) cloud water effective radius diagnosis
- Driver update (Joseph):
  - Implemented a new atmosphere driver to run SHiELD and SHiEMOM with the full FMScoupler.
- Updates to $2\delta z$ filter (fv_sg) (Lucas, Linjiong):
  - Included a missing term in the energy conservation formula (credit: Tristan Abbott).  May affect  prediction of processes depending strongly on microphysics. Compile the model with `-DENG_CNV_OLD` to revert this change.
  - Added option, `fv_sg_adj_weak`, to apply a weaker 2dz filter below sg_cutoff. This may be useful in controlling tropospheric instabilities without interfering with the behavior of the PBL scheme.
  - Renamed routines and eliminated ifdefs for SHiELD vs. AM4 versions.
- Physics interface updates (Linjiong, Kai, Spencer):
  - Fixed negative tracers in the dynamics-physics interface. 
  - Enhanced the fill_gfs function to remove negative tracers. 
  - Enabled data_override for nest domain 
  - Fixed a precipitation diagnostic issue when `ntimes > 1` in the GFDL MP.
  - MPI fix for sedimentation mass transport in GFDL MP.
- Updates to nudging (Lucas):
  - Added an option to turn TC breeding off. 
  - Bugfixes for nudging on a nest/regional domain (in which tendencies in the halo are undefined).
- Coarse-graining updates (Spencer, Kai):
  - Added options `strategy = 'pressure_level_extrapolate' ’blended_area_weighted’` (developed with support from Chris Bretherton, AI2), and simplest `model_level_area_weighted` (like FREgrid first-order conservative scheme). 
  - Renamed `model_level` strategy to `model_level_mass_weighted`.
  - Coarse-grained plev diagnostics for u, v, w, omega, vorticity, height, temperature, tracers, and RH.
  - Coarse-grained plev diagnostics use plevs defined in coarse-grained plev diagnostics for `fv_diag_plevs_nml`
  - OpenMP multi-threaded calculations
- Code refactors (Lucas):
  - Cleaned up `external_ic_nml` and `fv_surf_map_nml`.
  - Cleaned up `fv_mapz.F90` to move vertical remapping operators and thermodynamics/energetics routines into their own modules
- Diagnostics (Lucas, Linjiong, Kai, Spencer):
  - Fixes for nudging and fv_sg diagnostics 
  - Cleaned up fv_diagnostics stdout messages
  - True instantaneous and timestep-mean divergence and dissipative heating.
  - 40 dBz reflectivity height diagnostic.
  - Dissipative heating and dissipation estimate as, even if stochastic physics isn't enabled.
  - Introduced a flag `PRT_LEVEL` (now hard-coded) to control which min/max fields are written to stdout.
  - Fixed a bug for CAPE/CIN/BRN when nonhydrostatic pressure perturbation is also being output.
  - Refactor of plev and standard pressure level diagnostics, added new variables (vort, theta, theta_e, w, RH, dew point) to plevs, and removed unnecessary arguments to cs3_interpolator
- Deprecated/removed options (Lucas):
  - Removed outdated options: scale_z, w_max, w_limiter, z_min, d2_divg_max_k[12], damp_k_k[12], old_divg_damp, do_am4_remap, use_new_ncep, use_ncep_phy, a2b_ord, c2l_ord.
  - Interpolation from cell-means to corner values (a2b) and from local staggered winds to A-grid lat-lon winds, have been hard-coded to be fourth-order, except where it had previously been hard-coded to be second-order. Supporting codes have been cleaned up.
  - Deprecation notice for conserve_ke
  - Added warning messages for poorly-chosen advection scheme options (hord_xx), and a FATAL is thrown for invalid scheme choices.


# RELEASE NOTES for FV3 202305: Summary
FV3-202305-public --- May 2023  
Lucas Harris, GFDL lucas.harris@noaa.gov  

This version has been tested with SHiELD physics release 202305
and with FMS release 2023.01 from https://github.com/NOAA-GFDL/FMS

- Revised Vertical Remapping Operators (Lucas)
  - kord=10 reverted back to AM4 version.
  - Post-AM4 version of kord=10 is now kord=12.
  - do_am4_remap no longer does anything and is deprecated.
  - New strictly-monotone operators kord=14, 15 for improving tracer correlations, and kord=13 without subgrid limiting.
  - kord <= 7 now deprecated; may be removed in a future release.
- New Test Cases: (Joseph, Kun, Lucas)
  - Idealized TC test case with SHiELD physics
  - Zurita-Gotor et al. 2022 Held-Suarez variant
  - New Stable Boundary Layer (Beale at al.) doubly-periodic test case
- New nesting updates: (Joseph)
  - Enable nesting in solo core and add a new idealized test case (58)
  - Enable adding multiple nests in doubly-periodic test cases using absolute coordinates
- Additional idealized capability (Linjiong, Kun, Lucas)
  - Added namelist variable is_ideal_case, which must be used for runs starting (or re-starting) from idealized states.
  - Begin saving the initial wind fields (u0 and v0) to the restart files
- GFDL MP and Integrated Physics (Linjiong):
  - Added options to sub-cycling condensation evaporation (nconds), control timescale or evaporation (do_evap_timescale), and delay condensation and evaporation (delay_cond_evap)
  - Removed unused 3d microphysics diagnostics to save time and memory
  - Optimized the mpp domain updates for fast physics
  - Update gfdl_mp_nml reading code to avoid model crash for absent gfdl_mp_nml
  - Added an option (do_intermediate_phys) to disable intermediate phys
  - Removed grid size in GFDL MP energy and mass calculation
  - Updates to use dry_cp instead of moist_cp in a hydrostatic case
- Added a function to use O3 data from IFS ICs (Jan-Huey)
  - Namelist parameter: “use_gfsO3” with the default value = “false”
  - This function only works when ecmwf_ic = T
  - If the IFS IC does not include O3 data, or the run would like to use GFS O3 with other IFS ICs, set use_gfsO3 = T
- Solver Updates (Lucas)
  - Revised semi-implicit solver to partially linearize vertical sound wave propagation about the hydrostatic state. This removes a specific instability causing deep “columnar” modes in the vertical velocity field due to the equation for the pressure perturbation being updated partially forward-in-time. This removes the spurious modes, reduces vertical velocities, and makes the solver slightly more stable.
  - MPI bug fix for tracer diffusion
  - Fast Rayleigh Damping on w controlled by fast_tau_w_sec.


# RELEASE NOTES for FV3 202210: Summary
FV3-202210-public --- October 2022  
Lucas Harris, GFDL lucas.harris@noaa.gov  

This version has been tested with SHiELD physics release 202210
and with FMS release 2022.03 from https://github.com/NOAA-GFDL/FMS

This release includes the following:
- Release of the GFDL Microphysics Version 3
- Fix pressure-coarse-graining weighting from AI2's fork of FV3GFS
- Add A-grid restart functionality from AI2's fork of FV3GFS
- Fix for telescoping nest and GFS FIX file read
- Total precipitation diag field has changed from prec to pret
- Clean-up of the diagnostic messages to stdout


# RELEASE NOTES for FV3 202204: Summary
FV3-202204-public --- April 2022  
Lucas Harris, GFDL lucas.harris@noaa.gov  

This version has been tested against the current SHiELD physics
and with FMS release 2022.01 from https://github.com/NOAA-GFDL/FMS

This release includes the following:
- Release of stand-alone solo_core functionality with simple physics.
- Updated GFDL Microphysics, used for real-time 2021 C-SHiELD and T-SHiELD.  (L Zhou)
- Merges numerous updates from dev/emc.
- Leverage DA functionality from UFS with additional changes (M Tong).
- Updates to use the latest FMS release, including fms2_io.
- Adds license header to missing files and fixes typo in header.
- Fixes a bug where long_name and units attributes were not being captured in restart files.
- Adds the ability to specify prefix and directory when reading and writing restarts.
- The planetary radius and rotation rate are now re-scalable by a namelist parameter (small_earth_scale) instead of using exclusively the hard-coded FMS constant.
- Removes obsolete driver/SHiELD files.
- Removes unused function fv_diagnostics::max_vorticity_hy1.
- Removes avec timer remnants.
- Removes old style namelist read in favor of read from internal character variable.
- Adds option for a mean wind.
- Addresses GNU warnings.


# RELEASE NOTES for FV3 202107: Summary

FV3-202107-public --- 08 July 2021  
Lucas Harris, GFDL lucas.harris@noaa.gov  

This version has been tested against the current SHiELD physics
and with FMS release 2021.02 from https://github.com/NOAA-GFDL/FMS

This release includes the following:

- Comprehensive documentation in LaTEX format (FV3 team)
- Default changes to some namelist options and updated inline documentation
- Multiple same-level and telescoping nests for the Regional domain (J Mouallem)
- Updated fms2_io functionality (L Chilutti)
- Revised Regional domain code (K-Y Cheng)
- Reproducibility fixes for global+nests and regional+nests (tested for absolute reproducibility across PE counts, restarts)
- Other updates and general cleanup


# RELEASE NOTES for FV3 202101: Summary

FV3-202101-public --- 22 January 2021  
Lucas Harris, GFDL <lucas.harris@noaa.gov>  

This version has been tested against the current SHiELD (formerly fvGFS) physics
and with FMS release candidate 2020.04 from https://github.com/NOAA-GFDL/FMS

This release includes the following:

- Positive-definite advection scheme
- In-line GFDL Microphysics
- Fast-timescale Rayleigh damping
- Updated namelist documentation
- Implemented multiple same-level and telescoping nests for the global system (J Mouallem)
- Updated coarse-graining capabilities (S Clark)
- Re-organized fv_diagnostics, moving the revised fv_diag_column functionality and the declaration of diagnostic IDs to separate files
- and other updates and general cleanup

This version of FV3 is described as component of SHiELD in Harris et al. (2020, JAMES).

## Interface changes in 202101

atmosphere.F90: if using the in-line GFDL microphysics the precipitation rates (available in the structure Atm%inline_mp for rain, ice, snow, and graupel separately) must be passed into the physics and/or land model as appropriate. Here we demonstrate how to do this in SHiELD by copying them into IPD_Data(nb)%Statein%prep (and so on), which are newly defined in the IPD_Data structure within the SHiELD physics.


# RELEASE NOTES for FV3 201912: Summary

FV3-201912-public --- 10 January 2020
Lucas Harris, GFDL <lucas.harris@noaa.gov>

This version has been tested against the current SHiELD (formerly fvGFS) physics
and with FMS release candidate 2020.02 from https://github.com/NOAA-GFDL/FMS

Includes all of the features of the GFDL Release to EMC, as well as:

- Updated 2017 GFDL Microphysics (S-J Lin and L Zhou included in GFSv15)
- Updates for GFSv15 ICs (T Black/J Abeles, EMC)
- Updates to support new nesting capabilities in FMS (Z Liang)
- Re-written grid nesting code for efficiency and parallelization
- Re-organized fv_eta for improved vertical level selection
- 2018 Stand-alone regional capabilities (T Black/J Abeles, EMC)
- Refactored model front-end (fv_control, fv_restart)
- Support for point soundings
- And other updates

## Interface changes

drivers: renamed 'fvGFS' directory to SHiELD

atmosphere.F90: 'mytile' is renamed 'mygrid'

The non-functional gfdl_cloud_microphys.F90 has been removed and replaced with the 2017 public release given to EMC. Also added a proper initialization routine, that includes the use of INTERNAL_FILE_NML and thereby requires the input_nml_file argument. If you do not define the compiler flag INTERNAL_FILE_NML then you can delete this argument.

The namelist nggps_diag_nml has been eliminated. 'fdiag' is no longer handled by the dynamical core, and should be handled by the physics driver.

For a complete technical description see the NOAA Technical Memorandum OAR GFDL: https://repository.library.noaa.gov/view/noaa/23432
