# This is an automake file for the GFDL AM4 model.
# This builds directory ./model.
#
# Ed Hartnett 2019-09-10

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libtp_core.la libboundary.la libfv_current_grid.la libfv_arrays.la \
libfv_fill.la libfv_tracer2d.la libsw_core.la libfv_grid_utils.la \
libfv_mapz.la libfv_control.la libfv_nesting.la libdyn_core.la \
libfv_dynamics.la libfv_update_phys.la liba2b_edge.la libfv_sg.la \
libfv_cmp.la

# Each convenience library depends on its source.
libtp_core_la_SOURCES = tp_core.F90
libboundary_la_SOURCES = boundary.F90
libfv_current_grid_la_SOURCES = fv_current_grid.F90
libfv_arrays_la_SOURCES = fv_arrays.F90
libfv_fill_la_SOURCES = fv_fill.F90
libfv_tracer2d_la_SOURCES = fv_tracer2d.F90
libsw_core_la_SOURCES = sw_core.F90
libfv_grid_utils_la_SOURCES = fv_grid_utils.F90
libfv_mapz_la_SOURCES = fv_mapz.F90
libfv_control_la_SOURCES = fv_control.F90
libfv_nesting_la_SOURCES = fv_nesting.F90
libdyn_core_la_SOURCES = dyn_core.F90
libfv_dynamics_la_SOURCES = fv_dynamics.F90
libfv_update_phys_la_SOURCES = fv_update_phys.F90
liba2b_edge_la_SOURCES = a2b_edge.F90
libfv_sg_la_SOURCES = fv_sg.F90
libfv_cmp_la_SOURCES = fv_cmp.F90

# Each mod file depends on the .o file.
tp_core_mod.mod: tp_core.$(OBJEXT)
boundary_mod.mod: boundary.$(OBJEXT)
fv_current_grid_mod.mod: fv_current_grid.$(OBJEXT)
fv_arrays_mod.mod: fv_arrays.$(OBJEXT)
fv_fill_mod.mod: fv_fill.$(OBJEXT)
fv_tracer2d_mod.mod: fv_tracer2d.$(OBJEXT)
sw_core_mod.mod: sw_core.$(OBJEXT)
fv_grid_utils_mod.mod: fv_grid_utils.$(OBJEXT)
fv_mapz_mod.mod: fv_mapz.$(OBJEXT)
fv_control_mod.mod: fv_control.$(OBJEXT)
fv_nesting_mod.mod: fv_nesting.$(OBJEXT)
dyn_core_mod.mod: dyn_core.$(OBJEXT)
fv_dynamics_mod.mod: fv_dynamics.$(OBJEXT)
fv_update_phys_mod.mod: fv_update_phys.$(OBJEXT)
a2b_edge_mod.mod: a2b_edge.$(OBJEXT)
fv_sg_mod.mod: fv_sg.$(OBJEXT)
fv_cmp_mod.mod: fv_cmp.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
tp_core.$(OBJEXT): fv_arrays_mod.mod fv_grid_utils_mod.mod
boundary.$(OBJEXT): fv_arrays_mod.mod
fv_current_grid.$(OBJEXT): fv_arrays_mod.mod
fv_tracer2d.$(OBJEXT): tp_core_mod.mod boundary_mod.mod fv_arrays_mod.mod
sw_core.$(OBJEXT): tp_core_mod.mod fv_arrays_mod.mod a2b_edge_mod.mod
fv_grid_utils.$(OBJEXT): fv_arrays_mod.mod
fv_mapz.$(OBJEXT): fv_arrays_mod.mod fv_fill_mod.mod fv_grid_utils_mod.mod fv_cmp_mod.mod
fv_control.$(OBJEXT): fv_arrays_mod.mod fv_grid_utils_mod.mod
fv_nesting.$(OBJEXT): boundary_mod.mod fv_arrays_mod.mod sw_core_mod.mod fv_grid_utils_mod.mod fv_mapz_mod.mod fv_sg_mod.mod
dyn_core.$(OBJEXT): tp_core_mod.mod boundary_mod.mod fv_arrays_mod.mod sw_core_mod.mod fv_update_phys_mod.mod a2b_edge_mod.mod
fv_dynamics.$(OBJEXT): boundary_mod.mod fv_arrays_mod.mod fv_fill_mod.mod fv_tracer2d_mod.mod fv_grid_utils_mod.mod fv_mapz_mod.mod fv_nesting_mod.mod dyn_core_mod.mod fv_sg_mod.mod
fv_update_phys.$(OBJEXT): boundary_mod.mod fv_arrays_mod.mod fv_grid_utils_mod.mod fv_mapz_mod.mod
a2b_edge.$(OBJEXT): fv_arrays_mod.mod fv_grid_utils_mod.mod
fv_cmp.$(OBJEXT): fv_arrays_mod.mod

# Mod files are built and then installed as headers.
MODFILES = tp_core_mod.mod boundary_mod.mod fv_current_grid_mod.mod \
fv_arrays_mod.mod fv_fill_mod.mod fv_tracer2d_mod.mod sw_core_mod.mod \
fv_grid_utils_mod.mod fv_mapz_mod.mod fv_control_mod.mod \
fv_nesting_mod.mod dyn_core_mod.mod fv_dynamics_mod.mod \
fv_update_phys_mod.mod a2b_edge_mod.mod fv_sg_mod.mod fv_cmp_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
