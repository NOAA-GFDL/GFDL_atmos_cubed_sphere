# Precision-based Fortran compiler flags
set(R8_flags "-fdefault-real-8 -fdefault-double-8") # Fortran flags for 64BIT precision

# GNU Fortan
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ggdb -fbacktrace -cpp -fcray-pointer -ffree-line-length-none -fno-range-check ${${kind}_flags}")
if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch -fallow-invalid-boz")
endif()

set(CMAKE_Fortran_FLAGS_REPRO "-O2")

set(CMAKE_Fortran_FLAGS_RELEASE "-O2")

set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -ffpe-trap=invalid,zero,overflow -fbounds-check")

set(CMAKE_Fortran_LINK_FLAGS "" )
