#***********************************************************************
#*                   GNU Lesser General Public License
#*
#* This file is part of the FV3 dynamical core.
#*
#* The FV3 dynamical core is free software: you can redistribute it
#* and/or modify it under the terms of the
#* GNU Lesser General Public License as published by the
#* Free Software Foundation, either version 3 of the License, or
#* (at your option) any later version.
#*
#* The FV3 dynamical core is distributed in the hope that it will be
#* useful, but WITHOUT ANYWARRANTY; without even the implied warranty
#* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#* See the GNU General Public License for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with the FV3 dynamical core.
#* If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************
# Copyright (C) 2019 Geophysical Fluid Dynamics Laboratory (GFDL)
#
# Autor(s): Seth Underwood
#***********************************************************************

# This directory stores libtool macros, put there by aclocal.
ACLOCAL_AMFLAGS = -I m4

include $(MK_TEMPLATE)

lib_LTLIBRARIES = libGFDLacs.la

libGFDLacs_la_SOURCES = model/a2b_edge.F90 \
   model/boundary.F90 \
   model/dyn_core.F90 \
   model/fv_arrays.F90 \
   model/fv_cmp.F90 \
   model/fv_control.F90 \
   model/fv_current_grid.F90 \
   model/fv_dynamics.F90 \
   model/fv_fill.F90 \
   model/fv_grid_utils.F90 \
   model/fv_mapz.F90 \
   model/fv_nesting.F90 \
   model/fv_sg.F90 \
   model/fv_tracer2d.F90 \
   model/fv_update_phys.F90 \
   model/sw_core.F90 \
   model/tp_core.F90
libGFDLacs_la_SOURCES += model_nh/nh_core.F90 \
   model_nh/nh_utils.F90
libGFDLacs_la_SOURCES += tools/external_ic.F90 \
   tools/external_sst.F90 \
   tools/fv_diagnostics.F90 \
   tools/fv_eta.F90 \
   tools/fv_grid_tools.F90 \
   tools/fv_io.F90 \
   tools/fv_mp_mod.F90 \
   tools/fv_nggps_diag.F90 \
   tools/fv_nudge.F90 \
   tools/fv_restart.F90 \
   tools/fv_surf_map.F90 \
   tools/fv_timing.F90 \
   tools/init_hydro.F90 \
   tools/sim_nc_mod.F90 \
   tools/sorted_index.F90 \
   tools/test_cases.F90
libGFDLacs_la_SOURCES += GFDL_tools/fv_ada_nudge.F90 \
   GFDL_tools/fv_climate_nudge.F90 \
   GFDL_tools/fv_cmip_diag.F90 \
   GFDL_tools/read_climate_nudge_data.F90
libGFDLacs_la_SOURCES += driver/GFDL/atmosphere.F90
libGFDLacs_la_SOURCES += driver/SHiELD/atmosphere.F90

# Each Fortran mod file (FC_MODEXT) depends on the object file
model/a2b_edge_mod.$(FC_MODEXT) : model/a2b_edge.$(OBJEXT)
model/boundary_mod.$(FC_MODEXT) : model/boundary.$(OBJEXT)
model/dyn_core_mod.$(FC_MODEXT) : model/dyn_core.$(OBJEXT)
model/fv_arrays_mod.$(FC_MODEXT) : model/fv_arrays.$(OBJEXT)
model/fv_cmp_mod.$(FC_MODEXT) : model/fv_cmp.$(OBJEXT)
model/fv_control_mod.$(FC_MODEXT) : model/fv_control.$(OBJEXT)
model/fv_current_grid_mod.$(FC_MODEXT) : model/fv_current_grid.$(OBJEXT)
model/fv_dynamics_mod.$(FC_MODEXT) : model/fv_dynamics.$(OBJEXT)
model/fv_fill_mod.$(FC_MODEXT) : model/fv_fill.$(OBJEXT)
model/fv_grid_utils_mod.$(FC_MODEXT) : model/fv_grid_utils.$(OBJEXT)
model/fv_mapz_mod.$(FC_MODEXT) : model/fv_mapz.$(OBJEXT)
model/fv_nesting_mod.$(FC_MODEXT) : model/fv_nesting.$(OBJEXT)
model/fv_sg_mod.$(FC_MODEXT) : model/fv_sg.$(OBJEXT)
model/fv_tracer2d_mod.$(FC_MODEXT) : model/fv_tracer2d.$(OBJEXT)
model/fv_update_phys_mod.$(FC_MODEXT) : model/fv_update_phys.$(OBJEXT)
model/sw_core_mod.$(FC_MODEXT) : model/sw_core.$(OBJEXT)
model/tp_core_mod.$(FC_MODEXT) : model/tp_core.$(OBJEXT)
model_nh/nh_core_mod.$(FC_MODEXT) : model_nh/nh_core.$(OBJEXT)
model_nh/nh_utils_mod.$(FC_MODEXT) : model_nh/nh_utils.$(OBJEXT)
tools/external_ic_mod.$(FC_MODEXT) : tools/external_ic.$(OBJEXT)
tools/external_sst_mod.$(FC_MODEXT) : tools/external_sst.$(OBJEXT)
tools/fv_diagnostics_mod.$(FC_MODEXT) : tools/fv_diagnostics.$(OBJEXT)
tools/fv_eta_mod.$(FC_MODEXT) : tools/fv_eta.$(OBJEXT)
tools/fv_grid_tools_mod.$(FC_MODEXT) : tools/fv_grid_tools.$(OBJEXT)
tools/fv_io_mod.$(FC_MODEXT) : tools/fv_io.$(OBJEXT)
tools/fv_mp_mod.$(FC_MODEXT) : tools/fv_mp_mod.$(OBJEXT)
tools/fv_nggps_diag_mod.$(FC_MODEXT) : tools/fv_nggps_diag.$(OBJEXT)
tools/fv_nwp_nudge_mod.$(FC_MODEXT) : tools/fv_nudge.$(OBJEXT)
tools/fv_restart_mod.$(FC_MODEXT) : tools/fv_restart.$(OBJEXT)
tools/fv_surf_map_mod.$(FC_MODEXT) : tools/fv_surf_map.$(OBJEXT)
tools/fv_timing_mod.$(FC_MODEXT) : tools/fv_timing.$(OBJEXT)
tools/init_hydro_mod.$(FC_MODEXT) : tools/init_hydro.$(OBJEXT)
tools/sim_nc_mod.$(FC_MODEXT) : tools/sim_nc_mod.$(OBJEXT)
tools/sorted_index_mod.$(FC_MODEXT) : tools/sorted_index.$(OBJEXT)
tools/test_cases_mod.$(FC_MODEXT) : tools/test_cases.$(OBJEXT)
GFDL_tools/fv_ada_nudge_mod.$(FC_MODEXT) : GFDL_tools/fv_ada_nudge.$(OBJEXT)
GFDL_tools/fv_climate_nudge_mod.$(FC_MODEXT) : GFDL_tools/fv_climate_nudge.$(OBJEXT)
GFDL_tools/fv_cmip_diag_mod.$(FC_MODEXT) : GFDL_tools/fv_cmip_diag.$(OBJEXT)
GFDL_tools/read_climate_nudge_data_mod.$(FC_MODEXT) : GFDL_tools/read_climate_nudge_data.$(OBJEXT)
driver/GFDL/atmosphere_mod.$(FC_MODEXT) : driver/GFDL/atmosphere.$(OBJEXT)
driver/SHiELD/atmosphere_mod.$(FC_MODEXT) : driver/SHiELD/atmosphere.$(OBJEXT)

# Fortran modules/objects that depend on other Fortran modules/objects
model/a2b_edge.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT) \
   model/fv_grid_utils_mod.$(FC_MODEXT)
model/boundary.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT) \
   tools/fv_timing_mod.$(FC_MODEXT)
model/dyn_core.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT) \
   tools/fv_timing_mod.$(FC_MODEXT) \
   model/boundary_mod.$(FC_MODEXT) \
   model/tp_core_mod.$(FC_MODEXT) \
   model/a2b_edge_mod.$(FC_MODEXT) \
   tools/fv_diagnostics_mod.$(FC_MODEXT) \
   tools/test_cases_mod.$(FC_MODEXT) \
   model/sw_core_mod.$(FC_MODEXT) \
   model_nh/nh_core_mod.$(FC_MODEXT) \
   GFDL_tools/fv_ada_nudge_mod.$(FC_MODEXT) \
   tools/fv_nwp_nudge_mod.$(FC_MODEXT) \
   model/fv_update_phys_mod.$(FC_MODEXT)
model/fv_cmp.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT)
model/fv_control.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT) \
   tools/fv_timing_mod.$(FC_MODEXT) \
   tools/fv_eta_mod.$(FC_MODEXT) \
   model/fv_grid_utils_mod.$(FC_MODEXT) \
   tools/fv_io_mod.$(FC_MODEXT) \
   tools/fv_grid_tools_mod.$(FC_MODEXT) \
   tools/fv_diagnostics_mod.$(FC_MODEXT) \
   tools/test_cases_mod.$(FC_MODEXT) \
   tools/fv_restart_mod.$(FC_MODEXT)
model/fv_current_grid.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT)
model/fv_dynamics.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT) \
   model/fv_fill_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT) \
   model/fv_sg_mod.$(FC_MODEXT) \
   model/fv_grid_utils_mod.$(FC_MODEXT) \
   tools/fv_timing_mod.$(FC_MODEXT) \
   model/boundary_mod.$(FC_MODEXT) \
   model/fv_tracer2d_mod.$(FC_MODEXT) \
   tools/fv_diagnostics_mod.$(FC_MODEXT) \
   tools/fv_nwp_nudge_mod.$(FC_MODEXT) \
   model/dyn_core_mod.$(FC_MODEXT) \
   model/fv_mapz_mod.$(FC_MODEXT) \
   model/fv_nesting_mod.$(FC_MODEXT)
model/fv_grid_utils.$(OBJEXT) : model/external_sst_mod.$(FC_MODEXT) \
   model/fv_arrays_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT) \
   tools/fv_timing_mod.$(FC_MODEXT) \
   tools/fv_eta_mod.$(FC_MODEXT)
model/fv_mapz.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT) \
   model/fv_fill_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT) \
   model/fv_grid_utils_mod.$(FC_MODEXT) \
   tools/fv_timing_mod.$(FC_MODEXT) \
   model/fv_cmp_mod.$(FC_MODEXT)
model/fv_nesting.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT) \
   model/fv_sg_mod.$(FC_MODEXT) \
   tools/fv_timing_mod.$(FC_MODEXT) \
   model/fv_grid_utils_mod.$(FC_MODEXT) \
   model/boundary_mod.$(FC_MODEXT) \
   model/fv_mapz_mod.$(FC_MODEXT) \
   tools/init_hydro_mod.$(FC_MODEXT) \
   tools/fv_diagnostics_mod.$(FC_MODEXT) \
   model/sw_core_mod.$(FC_MODEXT) \
   tools/fv_restart_mod.$(FC_MODEXT)
model/fv_sg.$(OBJEXT) : tools/fv_mp_mod.$(FC_MODEXT)
model/fv_tracer2d.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT) \
   tools/fv_timing_mod.$(FC_MODEXT) \
   model/boundary_mod.$(FC_MODEXT) \
   model/tp_core_mod.$(FC_MODEXT)
model/fv_update_phys.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT) \
   tools/fv_timing_mod.$(FC_MODEXT) \
   tools/fv_eta_mod.$(FC_MODEXT) \
   model/fv_grid_utils_mod.$(FC_MODEXT) \
   model/boundary_mod.$(FC_MODEXT) \
   model/fv_mapz_mod.$(FC_MODEXT) \
   tools/fv_diagnostics_mod.$(FC_MODEXT) \
   GFDL_tools/fv_climate_nudge_mod.$(FC_MODEXT) \
   GFDL_tools/fv_ada_nudge_mod.$(FC_MODEXT) \
   tools/fv_nwp_nudge_mod.$(FC_MODEXT)
model/sw_core.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT) \
   model/tp_core_mod.$(FC_MODEXT) \
   model/a2b_edge_mod.$(FC_MODEXT) \
   tools/test_cases_mod.$(FC_MODEXT)
model/tp_core.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT) \
   model/fv_grid_utils_mod.$(FC_MODEXT)
model_nh/nh_core.$(OBJEXT) : model/nh_utils_mod.$(FC_MODEXT) \
   model/tp_core_mod.$(FC_MODEXT)
model_nh/nh_utils.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT) \
   model/sw_core_mod.$(FC_MODEXT) \
   model/tp_core_mod.$(FC_MODEXT)
tools/external_ic.$(OBJEXT) : model/external_sst_mod.$(FC_MODEXT) \
   model/fv_arrays_mod.$(FC_MODEXT) \
   model/fv_fill_mod.$(FC_MODEXT) \
   tools/sim_nc_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT) \
   tools/fv_timing_mod.$(FC_MODEXT) \
   tools/fv_eta_mod.$(FC_MODEXT) \
   model/fv_grid_utils_mod.$(FC_MODEXT) \
   model/boundary_mod.$(FC_MODEXT) \
   tools/fv_io_mod.$(FC_MODEXT) \
   tools/init_hydro_mod.$(FC_MODEXT) \
   tools/fv_surf_map_mod.$(FC_MODEXT) \
   tools/fv_diagnostics_mod.$(FC_MODEXT) \
   model/fv_mapz_mod.$(FC_MODEXT) \
   tools/fv_nwp_nudge_mod.$(FC_MODEXT) \
   tools/test_cases_mod.$(FC_MODEXT)
tools/fv_diagnostics.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT) \
   tools/fv_eta_mod.$(FC_MODEXT) \
   model/fv_grid_utils_mod.$(FC_MODEXT) \
   model/a2b_edge_mod.$(FC_MODEXT) \
   model/fv_mapz_mod.$(FC_MODEXT) \
   tools/fv_surf_map_mod.$(FC_MODEXT) \
   model/fv_sg_mod.$(FC_MODEXT)
tools/fv_eta.$(OBJEXT) : tools/fv_mp_mod.$(FC_MODEXT)
tools/fv_grid_tools.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT) \
   model/fv_grid_utils_mod.$(FC_MODEXT) \
   tools/fv_timing_mod.$(FC_MODEXT) \
   tools/sorted_index_mod.$(FC_MODEXT)
tools/fv_io.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT) \
   tools/fv_eta_mod.$(FC_MODEXT) \
   model/fv_mapz_mod.$(FC_MODEXT)
tools/fv_mp_mod.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT)
tools/fv_nggps_diag.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT) \
   tools/fv_diagnostics_mod.$(FC_MODEXT)
tools/fv_nudge.$(OBJEXT) : model/external_sst_mod.$(FC_MODEXT) \
   model/fv_arrays_mod.$(FC_MODEXT) \
   tools/sim_nc_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT) \
   model/fv_grid_utils_mod.$(FC_MODEXT) \
   tools/fv_timing_mod.$(FC_MODEXT) \
   model/tp_core_mod.$(FC_MODEXT) \
   model/fv_mapz_mod.$(FC_MODEXT) \
   tools/fv_diagnostics_mod.$(FC_MODEXT)
tools/fv_restart.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT) \
   tools/fv_timing_mod.$(FC_MODEXT) \
   tools/fv_eta_mod.$(FC_MODEXT) \
   model/fv_grid_utils_mod.$(FC_MODEXT) \
   model/boundary_mod.$(FC_MODEXT) \
   tools/fv_io_mod.$(FC_MODEXT) \
   tools/init_hydro_mod.$(FC_MODEXT) \
   tools/fv_surf_map_mod.$(FC_MODEXT) \
   tools/fv_diagnostics_mod.$(FC_MODEXT) \
   tools/test_cases_mod.$(FC_MODEXT) \
   tools/external_ic_mod.$(FC_MODEXT)
tools/fv_surf_map.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT) \
   model/fv_grid_utils_mod.$(FC_MODEXT) \
   tools/fv_timing_mod.$(FC_MODEXT)
tools/fv_timing.$(OBJEXT) : tools/fv_mp_mod.$(FC_MODEXT)
tools/init_hydro.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT) \
   model/fv_grid_utils_mod.$(FC_MODEXT)
tools/sorted_index.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT)
tools/test_cases.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT) \
   tools/fv_eta_mod.$(FC_MODEXT) \
   model/fv_grid_utils_mod.$(FC_MODEXT) \
   tools/init_hydro_mod.$(FC_MODEXT) \
   tools/fv_surf_map_mod.$(FC_MODEXT) \
   tools/fv_grid_tools_mod.$(FC_MODEXT) \
   model/fv_sg_mod.$(FC_MODEXT) \
   tools/fv_diagnostics_mod.$(FC_MODEXT)
GFDL_tools/fv_ada_nudge.$(OBJEXT) : model/external_sst_mod.$(FC_MODEXT) \
   tools/sim_nc_mod.$(FC_MODEXT) \
   model/fv_arrays_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT) \
   model/fv_grid_utils_mod.$(FC_MODEXT) \
   model/tp_core_mod.$(FC_MODEXT) \
   tools/fv_timing_mod.$(FC_MODEXT) \
   model/fv_mapz_mod.$(FC_MODEXT) \
   tools/fv_diagnostics_mod.$(FC_MODEXT)
GFDL_tools/fv_climate_nudge.$(OBJEXT) : GFDL_tools/read_climate_nudge_data_mod.$(FC_MODEXT) \
   model/fv_mapz_mod.$(FC_MODEXT)
GFDL_tools/fv_cmip_diag.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT) \
   model/fv_mapz_mod.$(FC_MODEXT) \
   tools/fv_diagnostics_mod.$(FC_MODEXT)
driver/GFDL/atmosphere.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT) \
   model/fv_sg_mod.$(FC_MODEXT) \
   tools/fv_timing_mod.$(FC_MODEXT) \
   tools/fv_eta_mod.$(FC_MODEXT) \
   model/boundary_mod.$(FC_MODEXT) \
   model/fv_dynamics_mod.$(FC_MODEXT) \
   tools/fv_io_mod.$(FC_MODEXT) \
   model/fv_nesting_mod.$(FC_MODEXT) \
   model/fv_update_phys_mod.$(FC_MODEXT) \
   GFDL_tools/fv_climate_nudge_mod.$(FC_MODEXT) \
   GFDL_tools/fv_ada_nudge_mod.$(FC_MODEXT) \
   tools/fv_nwp_nudge_mod.$(FC_MODEXT) \
   tools/fv_diagnostics_mod.$(FC_MODEXT) \
   GFDL_tools/fv_cmip_diag_mod.$(FC_MODEXT) \
   tools/fv_restart_mod.$(FC_MODEXT) \
   model/fv_control_mod.$(FC_MODEXT)
driver/SHiELD/atmosphere.$(OBJEXT) : model/fv_arrays_mod.$(FC_MODEXT) \
   model/fv_fill_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT) \
   model/fv_sg_mod.$(FC_MODEXT) \
   tools/fv_timing_mod.$(FC_MODEXT) \
   tools/fv_eta_mod.$(FC_MODEXT) \
   model/boundary_mod.$(FC_MODEXT) \
   model/fv_dynamics_mod.$(FC_MODEXT) \
   model/fv_nesting_mod.$(FC_MODEXT) \
   model/fv_update_phys_mod.$(FC_MODEXT) \
   tools/fv_nggps_diags_mod.$(FC_MODEXT) \
   tools/fv_nwp_nudge_mod.$(FC_MODEXT) \
   tools/fv_diagnostics_mod.$(FC_MODEXT) \
   tools/fv_restart_mod.$(FC_MODEXT) \
   model/fv_control_mod.$(FC_MODEXT)

BUILD_SOURCES = model/a2b_edge_mod.$(FC_MODEXT) \
   model/boundary_mod.$(FC_MODEXT) \
   model/dyn_core_mod.$(FC_MODEXT) \
   model/fv_arrays_mod.$(FC_MODEXT) \
   model/fv_cmp_mod.$(FC_MODEXT) \
   model/fv_control_mod.$(FC_MODEXT) \
   model/fv_current_grid_mod.$(FC_MODEXT) \
   model/fv_dynamics_mod.$(FC_MODEXT) \
   model/fv_fill_mod.$(FC_MODEXT) \
   model/fv_grid_utils_mod.$(FC_MODEXT) \
   model/fv_mapz_mod.$(FC_MODEXT) \
   model/fv_nesting_mod.$(FC_MODEXT) \
   model/fv_sg_mod.$(FC_MODEXT) \
   model/fv_tracer2d_mod.$(FC_MODEXT) \
   model/fv_update_phys_mod.$(FC_MODEXT) \
   model/sw_core_mod.$(FC_MODEXT) \
   model/tp_core_mod.$(FC_MODEXT) \
   model_nh/nh_core_mod.$(FC_MODEXT) \
   model_nh/nh_utils_mod.$(FC_MODEXT) \
   tools/external_ic_mod.$(FC_MODEXT) \
   tools/external_sst_mod.$(FC_MODEXT) \
   tools/fv_diagnostics_mod.$(FC_MODEXT) \
   tools/fv_eta_mod.$(FC_MODEXT) \
   tools/fv_grid_tools_mod.$(FC_MODEXT) \
   tools/fv_io_mod.$(FC_MODEXT) \
   tools/fv_mp_mod.$(FC_MODEXT) \
   tools/fv_nggps_diag_mod.$(FC_MODEXT) \
   tools/fv_nwp_nudge_mod.$(FC_MODEXT) \
   tools/fv_restart_mod.$(FC_MODEXT) \
   tools/fv_surf_map_mod.$(FC_MODEXT) \
   tools/fv_timing_mod.$(FC_MODEXT) \
   tools/init_hydro_mod.$(FC_MODEXT) \
   tools/sim_nc_mod.$(FC_MODEXT) \
   tools/sorted_index_mod.$(FC_MODEXT) \
   tools/test_cases_mod.$(FC_MODEXT) \
   GFDL_tools/fv_ada_nudge_mod.$(FC_MODEXT) \
   GFDL_tools/fv_climate_nudge_mod.$(FC_MODEXT) \
   GFDL_tools/fv_cmip_diag_mod.$(FC_MODEXT) \
   GFDL_tools/read_climate_nudge_data_mod.$(FC_MODEXT) \
   driver/GFDL/atmosphere_mod.$(FC_MODEXT) \
	driver/SHiELD/atmosphere_mod.$(FC_MODEXT)

include_HEADERS = fms_platform.h \
                  fms_version.h
BUILT_SOURCES =

CLEANFILES = *.mod

EXTRA_DIST = README.md \
             LICENSE.md
