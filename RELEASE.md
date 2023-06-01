# RELEASE NOTES for FV3 202305: Summary
FV3-202305-public --- May 2023
Lucas Harris, GFDL lucas.harris@noaa.gov

This version has been tested with SHiELD physics release 202305
and with FMS release 2023.01 from https://github.com/NOAA-GFDL/FMS

This release includes the following:
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
