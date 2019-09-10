# This is an automake file for the GFDL AM4 model.
# This builds directory ./tools.
#
# Ed Hartnett 2019-09-10

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libfv_surf_map.la libexternal_ic.la libfv_restart.la \
libfv_grid_tools.la libfv_eta.la libinit_hydro.la libsim_nc_mod.la \
libsorted_index.la libfv_mp_mod.la libexternal_sst.la libfv_timing.la \
libfv_diagnostics.la libfv_nudge.la libfv_nggps_diag.la libfv_io.la \
libtest_cases.la

# Each convenience library depends on its source.
libfv_surf_map_la_SOURCES = fv_surf_map.F90
libexternal_ic_la_SOURCES = external_ic.F90
libfv_restart_la_SOURCES = fv_restart.F90
libfv_grid_tools_la_SOURCES = fv_grid_tools.F90
libfv_eta_la_SOURCES = fv_eta.F90
libinit_hydro_la_SOURCES = init_hydro.F90
libsim_nc_mod_la_SOURCES = sim_nc_mod.F90
libsorted_index_la_SOURCES = sorted_index.F90
libfv_mp_mod_la_SOURCES = fv_mp_mod.F90
libexternal_sst_la_SOURCES = external_sst.F90
libfv_timing_la_SOURCES = fv_timing.F90
libfv_diagnostics_la_SOURCES = fv_diagnostics.F90
libfv_nudge_la_SOURCES = fv_nudge.F90
libfv_nggps_diag_la_SOURCES = fv_nggps_diag.F90
libfv_io_la_SOURCES = fv_io.F90
libtest_cases_la_SOURCES = test_cases.F90

# Each mod file depends on the .o file.
fv_surf_map_mod.mod: fv_surf_map.$(OBJEXT)
external_ic_mod.mod: external_ic.$(OBJEXT)
fv_restart_mod.mod: fv_restart.$(OBJEXT)
fv_grid_tools_mod.mod: fv_grid_tools.$(OBJEXT)
fv_eta_mod.mod: fv_eta.$(OBJEXT)
init_hydro_mod.mod: init_hydro.$(OBJEXT)
sim_nc_mod_mod.mod: sim_nc_mod.$(OBJEXT)
sorted_index_mod.mod: sorted_index.$(OBJEXT)
fv_mp_mod_mod.mod: fv_mp_mod.$(OBJEXT)
external_sst_mod.mod: external_sst.$(OBJEXT)
fv_timing_mod.mod: fv_timing.$(OBJEXT)
fv_diagnostics_mod.mod: fv_diagnostics.$(OBJEXT)
fv_nudge_mod.mod: fv_nudge.$(OBJEXT)
fv_nggps_diag_mod.mod: fv_nggps_diag.$(OBJEXT)
fv_io_mod.mod: fv_io.$(OBJEXT)
test_cases_mod.mod: test_cases.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
fv_surf_map.$(OBJEXT): fv_mp_mod_mod.mod fv_timing_mod.mod
external_ic.$(OBJEXT): fv_surf_map_mod.mod fv_eta_mod.mod init_hydro_mod.mod sim_nc_mod_mod.mod fv_mp_mod_mod.mod external_sst_mod.mod fv_timing_mod.mod fv_diagnostics_mod.mod fv_io_mod.mod test_cases_mod.mod
fv_restart.$(OBJEXT): fv_surf_map_mod.mod external_ic_mod.mod fv_eta_mod.mod init_hydro_mod.mod fv_mp_mod_mod.mod fv_timing_mod.mod fv_diagnostics_mod.mod fv_io_mod.mod test_cases_mod.mod
fv_grid_tools.$(OBJEXT): sorted_index_mod.mod fv_mp_mod_mod.mod fv_timing_mod.mod
fv_eta.$(OBJEXT): fv_mp_mod_mod.mod
init_hydro.$(OBJEXT): fv_mp_mod_mod.mod fv_diagnostics_mod.mod
fv_timing.$(OBJEXT): fv_mp_mod_mod.mod
fv_diagnostics.$(OBJEXT): fv_surf_map_mod.mod fv_eta_mod.mod fv_mp_mod_mod.mod
fv_nudge.$(OBJEXT): sim_nc_mod_mod.mod fv_mp_mod_mod.mod external_sst_mod.mod fv_timing_mod.mod fv_diagnostics_mod.mod
fv_nggps_diag.$(OBJEXT): fv_diagnostics_mod.mod
fv_io.$(OBJEXT): fv_eta_mod.mod fv_mp_mod_mod.mod external_sst_mod.mod
test_cases.$(OBJEXT): fv_surf_map_mod.mod fv_grid_tools_mod.mod fv_eta_mod.mod init_hydro_mod.mod fv_mp_mod_mod.mod fv_diagnostics_mod.mod

# Mod files are built and then installed as headers.
MODFILES = fv_surf_map_mod.mod external_ic_mod.mod fv_restart_mod.mod \
fv_grid_tools_mod.mod fv_eta_mod.mod init_hydro_mod.mod \
sim_nc_mod_mod.mod sorted_index_mod.mod fv_mp_mod_mod.mod \
external_sst_mod.mod fv_timing_mod.mod fv_diagnostics_mod.mod \
fv_nudge_mod.mod fv_nggps_diag_mod.mod fv_io_mod.mod \
test_cases_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
