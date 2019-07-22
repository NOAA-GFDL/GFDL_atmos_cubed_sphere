# This is the main configure file for the GFDL_atmos_cubed_sphere
# package.
# Ed Hartnett 7/22/2019

AC_PREREQ([2.59])

# Initialize with name, version, and support email address.
AC_INIT([GFDL_atmos_cubed_shere], [2.0-development], [])

# Find out about the host we're building on.
AC_CANONICAL_HOST

# Find out about the target we're building for.
AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE([foreign dist-zip subdir-objects])

# Keep libtool macros in an m4 directory.
AC_CONFIG_MACRO_DIR([m4])

# Set up libtool.
LT_PREREQ([2.4])
LT_INIT()

# Find the C compiler.
AC_PROG_CC
AM_PROG_CC_C_O
AC_C_CONST
AC_PROG_CPP

# Find the Fortran compiler.
AC_PROG_FC
AC_PROG_F77

# Find the install program.
AC_PROG_INSTALL

# Check to see if any macros must be set to enable large (>2GB) files.
AC_SYS_LARGEFILE

# Require MPI.
AC_CHECK_FUNC([MPI_Init], [], [AC_MSG_ERROR([MPI C library required])])

dnl # Check for netCDF C library.
dnl AC_SEARCH_LIBS([nc_create], [netcdf], [],
dnl                             [AC_MSG_ERROR([Can't find or link to the netcdf C library, set CPPFLAGS/LDFLAGS.])])

dnl # Check for netCDF Fortran library.
dnl AC_LANG_PUSH(Fortran)
dnl AC_SEARCH_LIBS([nf_create], [netcdff], [],
dnl                             [AC_MSG_ERROR([Can't find or link to the netcdf Fortran library, set CPPFLAGS/LDFLAGS.])])
dnl AC_LANG_POP(Fortran)

# Require netCDF.
#AC_CHECK_FUNC([nf_open], [], [AC_MSG_ERROR([NetCDF Fortran library required])])

# These defines are required for the build.
AC_DEFINE([use_netCDF], [1])
AC_DEFINE([use_libMPI], [1])

# These files will be created when the configure script is run.
AC_CONFIG_FILES([Makefile
        ])
AC_OUTPUT()
