# RELEASE NOTES for GFDL FV3: Summary
GFDL_2022.04 --- October 2022

This version has been tested against the current SHiELD physics
and with FMS release 2022.04 from https://github.com/NOAA-GFDL/FMS

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


# RELEASE NOTES for GFDL_2021.03.01: Summary
GFDL_2021.03.01 --- August 2021

This version has been tested within the current GFDL Models (AM4+, CM4+, ESM4+, SPEAR, etc.) and requires a
release of the [FMS infrastructure](https://github.com/NOAA-GFDL/FMS) 2021.03 or greater.

This release includes the following:

- Comprehensive documentation in LaTEX format (FV3 team)
- Default changes to some namelist options and updated inline documentation
- Multiple same-level and telescoping nests for the Regional domain (J Mouallem)
- Updated fms2_io functionality (L Chilutti)
- Revised Regional domain code (K-Y Cheng)
- Reproducibility fixes for global+nests and regional+nests (tested for absolute reproducibility across PE counts, restarts)
- Other updates and general cleanup


# RELEASE NOTES for GFDL_2021.01: Summary
GFDL_2021.01 --- April 2021

This version has been tested within the current GFDL Models (AM4+, CM4+, ESM4+, SPEAR, etc.) and requires a
release of the [FMS infrastructure](https://github.com/NOAA-GFDL/FMS) newer than 2020.02.

Includes all of the features from the [202101 Public Release](https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/releases/tag/FV3-202101-public) which include:

- Positive-definite advection scheme
- In-line GFDL Microphysics
- Fast-timescale Rayleigh damping
- Updated namelist documentation
- Implemented multiple same-level and telescoping nests for the global system (J Mouallem)
- Updated coarse-graining capabilities (S Clark)
- Re-organized fv_diagnostics, moving the revised fv_diag_column functionality and the declaration of diagnostic IDs to separate files
- and other updates and general cleanup

This version of FV3 is described as component of SHiELD in Harris et al. (2020, JAMES).


# RELEASE NOTES for 2020.02: Summary
2020.02 --- 22 April 2020

This version has been tested within current GFDL Models (AM4+, CM4+, ESM4+, SPEAR, etc.) and requires the 2020.02 release of the [FMS infrastructure](https://github.com/NOAA-GFDL/FMS).

Includes all of the features from the [201912 Public Release](https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/releases/tag/201912_public_release) which include:

- Updates to support new nesting capabilities in FMS
- Re-written grid nesting code for efficiency and parallelization
- Re-organized fv_eta for improved vertical level selection
- 2018 Stand-alone regional capabilities (from EMC)
- Refactored model front-end (fv_control, fv_restart)
- Support for point soundings
- full non-hydrostatic capability is now included as a runtime option
- And other updates

# Directory structure changes in 2020.02

***drivers/***  (important for those moving from the GFDL internal project)
  - renamed ***fvGFS*** to ***SHiELD***
  - renamed ***coupled*** to ***GFDL***

***model_nh_null/***
  - has been removed

Update your build system as appropriate
