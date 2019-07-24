# This is the automake file for GFDL_atmos_cubed_sphere GFDL_tools
# subdirectory.
# Ed Hartnett 7/24/2019

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libfv_ada_nudge.la libfv_climate_nudge.la	\
libfv_cmip_diag.la

# The convenience library depends on its source.
libfv_climate_nudge_la_SOURCES = fv_climate_nudge.F90
libfv_cmip_diag_la_SOURCES = fv_cmip_diag.F90
libfv_ada_nudge_la_SOURCES = fv_ada_nudge.F90

# Mod file depends on its o file, is built and then installed.
fv_cmip_diag_mod.mod: fv_cmip_diag.$(OBJEXT)
fv_climate_nudge_mod.mod: fv_climate_nudge.$(OBJEXT)
fv_ada_nudge_mod.mod: fv_ada_nudge.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
#atmosphere.$(OBJEXT): atmos_co2_mod.mod

MODFILES = fv_cmip_diag_mod.mod fv_climate_nudge_mod.mod	\
fv_ada_nudge_mod.mod
BUILT_SOURCES = ${MODFILES}
include_HEADERS = ${MODFILES}

CLEANFILES = *.mod
