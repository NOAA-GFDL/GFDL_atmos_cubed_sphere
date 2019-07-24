# This is the automake file for GFDL_atmos_cubed_sphere driver/SHIELD
# subdirectory.
# Ed Hartnett 7/24/2019

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = liblin_cld_microphys.la libconstants.la	\
libatmosphere.la

# The convenience library depends on its source.
libconstants_la_SOURCES = constants.F90
libatmosphere_la_SOURCES = atmosphere.F90
liblin_cld_microphys_la_SOURCES = lin_cld_microphys.F90

# Mod file depends on its o file, is built and then installed.
atmosphere_mod.mod: atmosphere.$(OBJEXT)
constants_mod.mod: constants.$(OBJEXT)
lin_cloud_microphys_mod.mod: lin_cloud_microphys.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
#atmosphere.$(OBJEXT): atmos_co2_mod.mod

MODFILES = atmosphere_mod.mod constants_mod.mod	\
lin_cloud_microphys_mod.mod
BUILT_SOURCES = ${MODFILES}
include_HEADERS = ${MODFILES}

CLEANFILES = *.mod
