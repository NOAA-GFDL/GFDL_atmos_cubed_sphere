# RELEASE NOTES for GFDL  FV3: Summary

2020.02 --- 22 April 2020

This version has been tested within current GFDL Models (AM, CM, ESM, SPEAR, etc.) and requires the 2020.02 release of the [FMS infrastructure](https://github.com/NOAA-GFDL/FMS).

Includes all of the features of the [201912 Public Release](https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/releases/tag/201912_public_release) which include:

- Updates to support new nesting capabilities in FMS
- Re-written grid nesting code for efficiency and parallelization
- Re-organized fv_eta for improved vertical level selection
- 2018 Stand-alone regional capabilities (from EMC)
- Refactored model front-end (fv_control, fv_restart)
- Support for point soundings
- And other updates

# Directory structure changes

***drivers/*** 
  - renamed ***fvGFS*** to ***SHiELD***
  - renamed ***coupled*** to ***GFDL***

***model_nh_null/*** 
  - has been removed 
  - non-hydrostatic capability is now included as a runtime configurable option

Update your build system as appropriate

# Documentation

For a complete technical description see the NOAA Technical Memorandum OAR GFDL: https://repository.library.noaa.gov/view/noaa/23432
