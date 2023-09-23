cmake_minimum_required(VERSION 3.18)

project(FV3
  VERSION 1.0.0
  DESCRIPTION  "GFDL Atmosphere Cubed Sphere Dynamical Core"
  HOMEPAGE_URL "https://www.gfdl.noaa.gov/gfdl_atmos_cubed_sphere"
  LANGUAGES Fortran)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(GNUInstallDirs)

if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|Repro|MinSizeRel|RelWithDebInfo)$")
  message(STATUS "No build type specified.")
  set(CMAKE_BUILD_TYPE "Release"
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
                     "Debug" "Release" "Repro" "MinSizeRel" "RelWithDebInfo")
endif()
message(STATUS "Setting build type to '${CMAKE_BUILD_TYPE}'.")

if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|GNU)$")
  message(WARNING "Compiler not officially supported: ${CMAKE_Fortran_COMPILER_ID}")
endif()

option(OPENMP      "Enable OpenMP threading"                               OFF)
option(32BIT       "Enable single precision (r4) arithmetic in FV3 dycore" ON)
option(DEBUG       "Enable compiler definition -DDEBUG"                    OFF)
option(MOVING_NEST "Enable compiler definition -DMOVING_NEST"              OFF)
option(MULTI_GASES "Enable compiler definition -DMULTI_GASES"              OFF)
option(USE_GFSL63  "Enable compiler definition -DUSE_GFSL63"               OFF)
option(GFS_PHYS    "Enable compiler definition -DGFS_PHYS"                 OFF)
option(GFS_TYPES   "Enable compiler definition -DGFS_TYPES"                OFF)
option(use_WRTCOMP "Enable compiler definition -Duse_WRTCOMP"              OFF)
option(INTERNAL_FILE_NML     "Enable compiler definition -DINTERNAL_FILE_NML"      ON)
option(ENABLE_QUAD_PRECISION "Enable compiler definition -DENABLE_QUAD_PRECISION"  ON)

find_package(MPI REQUIRED)
if(OPENMP)
  find_package(OpenMP REQUIRED)
endif()
find_package(NetCDF REQUIRED C Fortran)
if(use_WRTCOMP)
  if(NOT ESMF_FOUND)
  find_package(ESMF REQUIRED)
  endif()
endif()

if(32BIT)
  set(kind "R4")
else()
  set(kind "R8")
endif()

if(NOT FMS_FOUND)
  find_package(FMS REQUIRED COMPONENTS ${kind})
  string(TOLOWER ${kind} kind_lower)
  add_library(fms ALIAS FMS::fms_${kind_lower})
endif()

list(APPEND model_srcs
  model/a2b_edge.F90
  model/multi_gases.F90
  model/molecular_diffusion.F90
  model/boundary.F90
  model/dyn_core.F90
  model/fv_arrays.F90
  model/fv_control.F90
  model/fv_dynamics.F90
  model/fv_fill.F90
  model/fv_grid_utils.F90
  model/fv_mapz.F90
  model/fv_nesting.F90
  model/fv_regional_bc.F90
  model/fv_sg.F90
  model/fv_tracer2d.F90
  model/fv_update_phys.F90
  model/sw_core.F90
  model/tp_core.F90
  model/nh_core.F90
  model/nh_utils.F90)

list(APPEND tools_srcs
  tools/coarse_grained_diagnostics.F90
  tools/coarse_grained_restart_files.F90
  tools/coarse_graining.F90
  tools/external_ic.F90
  tools/external_sst.F90
  tools/fv_diag_column.F90
  tools/fv_diagnostics.F90
  tools/fv_eta.F90
  tools/fv_grid_tools.F90
  tools/fv_io.F90
  tools/fv_mp_mod.F90
  tools/fv_nudge.F90
  tools/fv_treat_da_inc.F90
  tools/fv_restart.F90
  tools/fv_surf_map.F90
  tools/fv_timing.F90
  tools/init_hydro.F90
  tools/module_diag_hailcast.F90
  tools/sim_nc_mod.F90
  tools/sorted_index.F90
  tools/test_cases.F90)

list(APPEND tools_srcs_extra
  tools/fv_iau_mod.F90)

list(APPEND driver_srcs
  driver/fvGFS/DYCORE_typedefs.F90
  driver/fvGFS/fv_nggps_diag.F90
  driver/fvGFS/fv_ufs_restart_io.F90
  driver/fvGFS/atmosphere.F90)

list(APPEND fv3_srcs ${model_srcs}
                     ${tools_srcs})

list(APPEND fv3_defs SPMD
                     MOIST_CAPPA
                     USE_COND)

# Additional (optional) compiler definitions
if(DEBUG)
  list(APPEND fv3_defs DEBUG)
endif()

if(USE_GFSL63)
  list(APPEND fv3_defs USE_GFSL63)
endif()

if(GFS_PHYS)
  list(APPEND fv3_defs GFS_PHYS)
endif()

if(GFS_TYPES)
  list(APPEND fv3_defs GFS_TYPES)
endif()

if(use_WRTCOMP)
  list(APPEND fv3_defs use_WRTCOMP)
  list(APPEND fv3_srcs ${tools_srcs_extra}
                       ${driver_srcs})
endif()

if(MOVING_NEST)
  list(APPEND fv3_defs MOVING_NEST)
endif()

if(MULTI_GASES)
  list(APPEND fv3_defs MULTI_GASES)
endif()

if(32BIT)
  list(APPEND fv3_defs OVERLOAD_R4
                       OVERLOAD_R8)
endif()

if(INTERNAL_FILE_NML)
  list(APPEND fv3_defs INTERNAL_FILE_NML)
endif()

if(ENABLE_QUAD_PRECISION)
  list(APPEND fv3_defs ENABLE_QUAD_PRECISION)
endif()

if(OPENMP)
  list(APPEND fv3_defs OPENMP)
endif()

# Obtain compiler-specific flags
include(fv3_compiler_flags)

# Enable -fPIC compiler flag for all targets
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_library(fv3 ${fv3_srcs})
add_library(FV3::fv3 ALIAS fv3)

set_property(SOURCE model/nh_utils.F90 APPEND_STRING PROPERTY COMPILE_FLAGS "${FAST}")
set_property(SOURCE model/fv_mapz.F90  APPEND_STRING PROPERTY COMPILE_FLAGS "${FAST}")

set_target_properties(fv3 PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/fv3)
target_include_directories(fv3 PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include/fv3>
                                         $<INSTALL_INTERFACE:include/fv3>)

target_include_directories(fv3 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
                                       ${CMAKE_CURRENT_SOURCE_DIR}/tools)

target_compile_definitions(fv3 PRIVATE "${fv3_defs}")

target_link_libraries(fv3 PUBLIC fms)

if(GFS_PHYS)
  target_link_libraries(fv3 PUBLIC fv3ccpp)
endif()

if(use_WRTCOMP)
  target_link_libraries(fv3 PUBLIC esmf)
endif()

if(OPENMP)
  target_link_libraries(fv3 PUBLIC OpenMP::OpenMP_Fortran)
endif()

# Install compiled Fortran module files
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})

install(
  TARGETS fv3
  EXPORT FV3Exports
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

### Package config
include(CMakePackageConfigHelpers)
set(CONFIG_INSTALL_DESTINATION lib/cmake/fv3)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/FV3Config.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/fv3-config.cmake
  INSTALL_DESTINATION ${CONFIG_INSTALL_DESTINATION})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/fv3-config.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/fv3-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/fv3-config-version.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

install(EXPORT FV3Exports
  NAMESPACE FV3::
  FILE fv3-targets.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

