# Precision-based Fortran compiler flags
set(r4_flags "-real-size 32") # Fortran flags for 32BIT precision
set(r8_flags "-real-size 64") # Fortran flags for 64BIT precision
set(r8_flags "${r8_flags} -no-prec-div -no-prec-sqrt")

# Intel Fortran
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback -fpp -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -nowarn -sox -align array64byte -qno-opt-dynamic-align")

set(CMAKE_Fortran_FLAGS_REPRO "-O2 -debug minimal -fp-model consistent -qoverride-limits")

set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -debug minimal -fp-model source -qoverride-limits -qopt-prefetch=3")

set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fp-stack-check -fstack-protector-all -fpe0 -debug -ftrapuv")

set(CMAKE_Fortran_LINK_FLAGS "")

set(FAST "-fast-transcendentals")
