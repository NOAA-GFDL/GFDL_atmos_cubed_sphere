# This is an automake file for the GFDL AM4 model.
# This builds directory ./GFDL_tools.
#
# Ed Hartnett 2019-09-10

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libread_climate_nudge_data.la libfv_cmip_diag.la \
libfv_climate_nudge.la libfv_ada_nudge.la

# Each convenience library depends on its source.
libread_climate_nudge_data_la_SOURCES = read_climate_nudge_data.F90
libfv_cmip_diag_la_SOURCES = fv_cmip_diag.F90
libfv_climate_nudge_la_SOURCES = fv_climate_nudge.F90
libfv_ada_nudge_la_SOURCES = fv_ada_nudge.F90

# Each mod file depends on the .o file.
read_climate_nudge_data_mod.mod: read_climate_nudge_data.$(OBJEXT)
fv_cmip_diag_mod.mod: fv_cmip_diag.$(OBJEXT)
fv_climate_nudge_mod.mod: fv_climate_nudge.$(OBJEXT)
fv_ada_nudge_mod.mod: fv_ada_nudge.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
fv_climate_nudge.$(OBJEXT): read_climate_nudge_data_mod.mod

# Mod files are built and then installed as headers.
MODFILES = read_climate_nudge_data_mod.mod fv_cmip_diag_mod.mod \
fv_climate_nudge_mod.mod fv_ada_nudge_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
