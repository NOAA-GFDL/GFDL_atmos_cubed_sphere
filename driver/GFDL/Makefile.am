# This is an automake file for the GFDL_atmos_cubed_sphere/driver/GFDL
# directory.

# Ed Hartnett 7/22/19

# Include .h and .mod files.
# AM_CPPFLAGS = -I${top_srcdir}/include
# AM_CPPFLAGS += -I${top_builddir}/fms
# AM_CPPFLAGS += -I${top_builddir}/time_manager
# AM_CPPFLAGS += -I${top_builddir}/constants
# AM_CPPFLAGS += -I${top_builddir}/mpp

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libatmos_co2.la libatmosphere.la

# The convenience library depends on its source.
libatmosphere_la_SOURCES = atmosphere.F90
libatmos_co2_la_SOURCES = atmos_co2.F90

# Mod file depends on its o file, is built and then installed.
atmosphere_mod.mod: atmosphere.$(OBJEXT)
atmos_co2_mod.mod: atmos_co2.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
atmosphere.$(OBJEXT): atmos_co2_mod.mod

BUILT_SOURCES = atmosphere_mod.mod atmos_co2_mod.mod
include_HEADERS = atmosphere_mod.mod atmos_co2_mod.mod

CLEANFILES = *.mod
