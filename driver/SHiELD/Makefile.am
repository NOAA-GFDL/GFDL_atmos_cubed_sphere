# This is an automake file for the GFDL AM4 model.
# This builds directory ./driver/SHiELD.
#
# Ed Hartnett 2019-09-10

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libatmosphere.la libconstants.la liblin_cloud_microphys.la

# Each convenience library depends on its source.
libatmosphere_la_SOURCES = atmosphere.F90
libconstants_la_SOURCES = constants.F90
liblin_cloud_microphys_la_SOURCES = lin_cloud_microphys.F90

# Each mod file depends on the .o file.
atmosphere_mod.mod: atmosphere.$(OBJEXT)
constants_mod.mod: constants.$(OBJEXT)
lin_cloud_microphys_mod.mod: lin_cloud_microphys.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
atmosphere.$(OBJEXT): constants_mod.mod
lin_cloud_microphys.$(OBJEXT): constants_mod.mod

# Mod files are built and then installed as headers.
MODFILES = atmosphere_mod.mod constants_mod.mod lin_cloud_microphys_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
