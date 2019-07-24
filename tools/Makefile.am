# This is the automake file for GFDL_atmos_cubed_sphere tools
# subdirectory.
# Ed Hartnett 7/24/2019

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libexternal_ic.la libexternal_sst.la		\
libfv_diagnostics.la libfv_eta.la libfv_grid_tools.la libfv_io.la	\
libfv_mp_mod.la libfv_nggps_diag.la libfv_nudge.la libfv_restart.la	\
libfv_surf_map.la libfv_timing.la libinit_hydro.la libsim_nc_mod.la	\
libsorted_index.la libtest_cases.la

#The convenience libries depends on their source.
libtest_cases_la_SOURCES = test_cases.F90
libsorted_index_la_SOURCES = sorted_index.F90
libsim_nc_mod_la_SOURCES = sim_nc_mod.F90
libinit_hydro_la_SOURCES = init_hydro.F90
libfv_timing_la_SOURCES = fv_timing.F90
libfv_surf_map_la_SOURCES = fv_surf_map.F90
libfv_restart_la_SOURCES = fv_restart.F90
libfv_nudge_la_SOURCES = fv_nudge.F90
libfv_nggps_diag_la_SOURCES = fv_nggps_diag.F90
libfv_mp_mod_la_SOURCES = fv_mp_mod.F90
libfv_io_la_SOURCES = fv_io.F90
libfv_grid_tools_la_SOURCES = fv_grid_tools.F90
libfv_eta_la_SOURCES = fv_eta.F90
libfv_diagnostics_la_SOURCES = fv_diagnostics.F90
libexternal_sst_la_SOURCES = external_sst.F90
libexternal_ic_la_SOURCES = external_ic.F90

# Mod file depends on its o file, is built and then installed.
external_ic_mod.mod: external_ic.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
#atmosphere.$(OBJEXT): atmos_co2_mod.mod

MODFILES = external_ic_mod.mod external_sst_mod.mod		\
fv_diagnostics_mod.mod fv_eta_mod.mod fv_grid_tools_mod.mod	\
fv_io_mod.mod fv_mp_mod_mod.mod fv_nggps_diag_mod.mod		\
fv_nudge_mod.mod fv_restart_mod.mod fv_surf_map_mod.mod		\
fv_timing_mod.mod init_hydro_mod.mod sim_nc_mod_mod.mod		\
sorted_index_mod.mod test_cases_mod.mod
BUILT_SOURCES = ${MODFILES}
include_HEADERS = ${MODFILES}

CLEANFILES = *.mod
