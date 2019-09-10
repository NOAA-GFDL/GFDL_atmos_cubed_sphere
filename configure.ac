# This is the main configure file for the GFDL AM4
# package.
# Ed Hartnett 7/29/2019

AC_PREREQ([2.59])

# Initialize with name, version, and support email address.
AC_INIT([AM4], [2.0-development], [])

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

# Require FMS. But this does not work because fms_platform.h is
# actually a fortran file, and AC_SEARCH_LIBS does not work with
# fortran mod files. I will swing round to this again, until then
# include the FMS library in FCFLAGS, LDFLAGS and LIBS.
# AC_CHECK_HEADERS([fms_platform.h], [], [AC_MSG_ERROR([Can't find fms_platform.h.])])
# AC_SEARCH_LIBS([__mpp_mod_MOD_mpp_init], [FMS], [],
#                            [AC_MSG_ERROR([Can't find or link to the FMS Fortran library, set CPPFLAGS/LDFLAGS.])])


# Check for netCDF C library.
AC_SEARCH_LIBS([nc_create], [netcdf], [],
                            [AC_MSG_ERROR([Can't find or link to the netcdf C library, set CPPFLAGS/LDFLAGS.])])

# Check for netCDF Fortran library.
AC_LANG_PUSH(Fortran)
AC_SEARCH_LIBS([nf_create], [netcdff], [],
                            [AC_MSG_ERROR([Can't find or link to the netcdf Fortran library, set CPPFLAGS/LDFLAGS.])])
AC_LANG_POP(Fortran)

# These defines are required for the build.
#AC_DEFINE([use_netCDF], [1])
#AC_DEFINE([use_libMPI], [1])
#AC_DEFINE([SPMD], [1])

# These files will be created when the configure script is run.
AC_CONFIG_FILES([Makefile
        GFDL_tools/Makefile
        model_nh/Makefile
        tools/Makefile
        driver/Makefile
        driver/GFDL/Makefile
        driver/SHiELD/Makefile
        model/Makefile
        ])
AC_OUTPUT()
