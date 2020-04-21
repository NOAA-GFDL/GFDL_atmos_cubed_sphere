# RELEASE NOTES for FV3: Summary

FV3-201912-public --- 10 January 2020
Lucas Harris, GFDL <lucas.harris@noaa.gov>

This version has been tested against the current GFDL AM4 phsyics
and with release candidate 2020.02 of the FMS code at https://github.com/NOAA-GFDL/FMS

Includes all of the features of the GFDL Release to EMC, as well as:

- Updated 2017 GFDL Microphysics (from S-J Lin and L Zhou included in GFSv15)
- Updates for GFSv15 ICs (from T Black/J Abeles, EMC)
- Updates to support new nesting capabilities in FMS (from Z Liang)
- Re-written grid nesting code for efficiency and parallelization
- Re-organized fv_eta for improved vertical level selection
- 2018 Stand-alone regional capabilities (from T Black/J Abeles, EMC)
- Refactored model front-end (fv_control, fv_restart)
- Support for point soundings
- And other updates

# Interface changes

drivers: renamed 'fvGFS' directory to SHiELD

atmosphere.F90: 'mytile' is renamed 'mygrid'

The non-functional gfdl_cloud_microphys.F90 has been removed and replaced with the 2017 public release given to EMC. Also added a proper initialization routine, that includes the use of INTERNAL_FILE_NML and thereby requires the input_nml_file argument. If you do not define the compiler flag INTERNAL_FILE_NML then you can delete this argument.

The namelist nggps_diag_nml has been eliminated. 'fdiag' is no longer handled by the dynamical core, and should be handled by the physics driver.

For a complete technical description see the NOAA Technical Memorandum OAR GFDL: https://repository.library.noaa.gov/view/noaa/23432
