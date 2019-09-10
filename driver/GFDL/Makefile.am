# This is an automake file for the GFDL AM4 model.
# This builds directory ./driver/GFDL.
#
# Ed Hartnett 2019-09-10

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libatmosphere.la

# Each convenience library depends on its source.
libatmosphere_la_SOURCES = atmosphere.F90

# Each mod file depends on the .o file.
atmosphere_mod.mod: atmosphere.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = atmosphere_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
