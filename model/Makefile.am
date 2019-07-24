# This is the automake file for GFDL_atmos_cubed_sphere model
# subdirectory.
# Ed Hartnett 7/24/2019

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = liba2b_edge.la libboundary.la libdyn_core.la	\
libfv_arrays.la libfv_cmp.la libfv_control.la libfv_current_grid.la	\
libfv_dynamics.la libfv_fill.la libfv_grid_utils.la libfv_mapz.la	\
libfv_nesting.la libfv_sg.la libfv_tracer2d.la libfv_update_phys.la	\
libsw_core.la libtp_core.la

# The convenience libries depends on their source.
libtp_core_la_SOURCES = tp_core.F90
libsw_core_la_SOURCES = sw_core.F90
libfv_update_phys_la_SOURCES = fv_update_phys.F90
libfv_tracer2d_la_SOURCES = fv_tracer2d.F90
libfv_sg_la_SOURCES = fv_sg.F90
libfv_nesting_la_SOURCES = fv_nesting.F90
libfv_mapz_la_SOURCES = fv_mapz.F90
libfv_grid_utils_la_SOURCES = fv_grid_utils.F90
libfv_fill_la_SOURCES = fv_fill.F90
libfv_dynamics_la_SOURCES = fv_dynamics.F90
libfv_current_grid_la_SOURCES = fv_current_grid.F90
libfv_control_la_SOURCES = fv_control.F90
libfv_cmp_la_SOURCES = fv_cmp.F90
libfv_arrays_la_SOURCES = fv_arrays.F90
libdyn_core_la_SOURCES = dyn_core.F90
libboundary_la_SOURCES = boundary.F90
liba2b_edge_la_SOURCES = a2b_edge.F90

# Mod file depends on its o file, is built and then installed.
a2b_edge_mod.mod: a2b_edge.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
#atmosphere.$(OBJEXT): atmos_co2_mod.mod

MODFILES = a2b_edge_mod.mod boundary_mod.mod dyn_core_mod.mod	\
fv_arrays_mod.mod fv_cmp_mod.mod fv_control_mod.mod		\
fv_current_grid_mod.mod fv_dynamics_mod.mod fv_fill_mod.mod	\
fv_grid_utils_mod.mod fv_mapz_mod.mod fv_nesting_mod.mod	\
fv_sg_mod.mod fv_tracer2d_mod.mod fv_update_phys_mod.mod	\
sw_core_mod.mod tp_core_mod.mod
BUILT_SOURCES = ${MODFILES}
include_HEADERS = ${MODFILES}

CLEANFILES = *.mod
