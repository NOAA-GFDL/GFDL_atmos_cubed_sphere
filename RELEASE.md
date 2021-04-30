# RELEASE NOTES for GFDL FV3: Summary

# RELEASE NOTES for GFDL_2021.01: Summary
GFDL_2021.01 --- April 2021

This version has been tested within the current GFDL Models (AM4+, CM4+, ESM4+, SPEAR, etc.) and requires a
release of the [FMS infrastructure](https://github.com/NOAA-GFDL/FMS) newer than 2020.02.

Includes all of the features from the [202101 Public Release](https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/releases/tag/FV3-202101-public) which include:

- Positive-definite advection scheme
- In-line GFDL Microphysics
- Fast-timescale Rayleigh damping
- Updated namelist documentation
- Implemented multiple same-level and telescoping nests for the global system (from J Mouallem)
- Updated coarse-graining capabilities (from S Clark)
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

# Documentation

For a complete technical description see the NOAA Technical Memorandum OAR GFDL: https://repository.library.noaa.gov/view/noaa/23432
