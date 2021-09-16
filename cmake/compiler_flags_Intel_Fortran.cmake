# Precision-based Fortran compiler flags
set(R4_flags "-real-size 32") # Fortran flags for 32BIT precision
set(R8_flags "-real-size 64") # Fortran flags for 64BIT precision
set(R8_flags "${R8_flags} -no-prec-div -no-prec-sqrt")

# Intel Fortran
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback -fpp -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -nowarn -sox -align array64byte -qno-opt-dynamic-align ${${kind}_flags}")

set(CMAKE_Fortran_FLAGS_REPRO "-O2 -debug minimal -fp-model consistent -qoverride-limits")

set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -debug minimal -fp-model source -qoverride-limits -qopt-prefetch=3")

set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fp-stack-check -fstack-protector-all -fpe0 -debug -ftrapuv")

set(CMAKE_Fortran_LINK_FLAGS "")

if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Repro)$")
  set(FAST "-fast-transcendentals")
  if(AVX2)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=core-avx2")
  elseif(SIMDMULTIARCH)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -axSSE4.2,CORE-AVX2")
  endif()
endif()

