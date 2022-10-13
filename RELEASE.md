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
