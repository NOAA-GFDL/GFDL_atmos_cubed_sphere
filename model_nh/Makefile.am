# This is an automake file for the GFDL AM4 model.
# This builds directory ./model_nh.
#
# Ed Hartnett 2019-09-10

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libnh_utils.la libnh_core.la

# Each convenience library depends on its source.
libnh_utils_la_SOURCES = nh_utils.F90
libnh_core_la_SOURCES = nh_core.F90

# Each mod file depends on the .o file.
nh_utils_mod.mod: nh_utils.$(OBJEXT)
nh_core_mod.mod: nh_core.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
nh_core.$(OBJEXT): nh_utils_mod.mod

# Mod files are built and then installed as headers.
MODFILES = nh_utils_mod.mod nh_core_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
