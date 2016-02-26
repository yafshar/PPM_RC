#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

## Take all the options with the exception of --enable-install-req
AC_PREREQ(2.59)

AC_INIT(PPM_RC, 1.0, [afshar@mpi-cbg.de],[],[https://github.com/yafshar/PPM_RC])
AC_CONFIG_SRCDIR([src/ppm_rc.f])
AM_INIT_AUTOMAKE

# m4 include
m4_ifdef([AX_BOOST],,[m4_include([m4/ax_boost.m4])])
m4_ifdef([PPM_RC_BOOST],,[m4_include([m4/ax_boost.m4])])
m4_ifdef([AX_METIS],,[m4_include([m4/ax_metis.m4])])
m4_ifdef([AC_CHECK_LIBTIFF],,[m4_include([m4/ax_tiff.m4])])
m4_ifdef([AX_F90_INTERNAL_HEADMOD],,[m4_include([m4/ax_f90_internal_headmod.m4])])
m4_ifdef([AX_F90_MODULE],,[m4_include([m4/ax_f90_module.m4])])
m4_ifdef([AX_F90_LIBRARY],,[m4_include([m4/ax_f90_library.m4])])
m4_ifdef([AX_F90_MODULE_EXTENSION],,[m4_include([m4/ax_f90_module_extension.m4])])
m4_ifdef([AX_F90_HEADER],,[m4_include([m4/ax_f90_header.m4])])
m4_ifdef([PPM_LIBRARY_SETUP],,[m4_include([m4/ppm_library_setup.m4])])
m4_ifdef([AX_F90_MODULE_FLAG],,[m4_include([m4/ax_f90_module_flag.m4])])
m4_ifdef([AX_GCC_VERSION],,[m4_include([m4/ax_gcc_version.m4])])
m4_ifdef([AX_GCC_X86_CPUID],,[m4_include([m4/ax_gcc_x86_cpuid.m4])])
m4_ifdef([AX_GCC_ARCHFLAG],,[m4_include([m4/ax_gcc_archflag.m4])])
m4_ifdef([AX_CHECK_COMPILER_FLAGS],,[m4_include([m4/ax_check_compiler_flags.m4])])
m4_ifdef([AX_FC_VERSION],,[m4_include([m4/ax_fc_version.m4])])
m4_ifdef([AX_FC_MODULE_FLAG],,[m4_include([m4/ax_fc_module_flag.m4])])
m4_ifdef([PPM_RC_RUNLOG],,[m4_include([m4/ppm_rc_runlog.m4])])
m4_ifdef([PPM_RC_CPP],,[m4_include([m4/ppm_rc_cpp.m4])])
m4_ifdef([PPM_RC_MPI],,[m4_include([m4/ppm_rc_mpi.m4])])
m4_ifdef([PPM_RC_NOMPI],,[m4_include([m4/ppm_rc_mpi.m4])])

PPM_RC_FCFLAGS=$FCFLAGS
PPM_RC_CFLAGS=$CFLAGS
PPM_RC_CXX_FLAGS=$CXXFLAGS

# Checks for basic compilers 
AC_PROG_FC
AC_PROG_CC
AC_PROG_CXX

AC_PROG_CPP
# PPM_RC needs the preprocessor. The result goes to CPP
PPM_RC_CPP

AC_PROG_SED
AC_PROG_MAKE_SET
AC_PROG_RANLIB

# if MPI is enabled try grabbing those compiler wrappers first
use_mpi=no
AC_ARG_ENABLE([mpi],
      [AS_HELP_STRING([--enable-mpi],
		      [use MPI (default is no), If the MPI implementation of your choice provides
		       compiler wrappers that are in the PATH, I can set them myself,  
		       choose: guess (I will choose the first implementation I can find). 	
		       Else, set CC, CXX and FC to the appropriate compiler wrappers (safest)])],
		       [PPM_RC_MPI
			# compilers should be chosen
			use_mpi=yes],
			[AC_MSG_NOTICE([MPI was not enabled.])])


if [ test x"$use_mpi" = x"no" ]; then
	PPM_RC_NOMPI
fi
if [ test x"$use_mpi" = x"no" ]; then
	echo "***********************************************************************"
	echo "*    MPI was not enabled. Building non-parallel version of PPM_RC.    *"
	echo "***********************************************************************"
fi

AC_SUBST(CFLAGS,["$PPM_RC_CFLAGS"])
AC_SUBST(CXXFLAGS,["$PPM_RC_CXXFLAGS"])
AC_SUBST(FCFLAGS,["$PPM_RC_FCFLAGS"])

# Checks for basic compilers and flags
AC_FC_FREEFORM
AC_FC_LIBRARY_LDFLAGS
# Check architecture
AC_CHECK_SIZEOF([long double])

AX_FC_MODULE_FLAG
AC_SUBST([ax_cv_fc_modflag])

AC_SUBST(CPPFLAGS,[""])

# Check for Fortran2003 support
AC_MSG_CHECKING([Whether your compiler properly supports Fortran 2003])
AC_LANG([Fortran])
use_f03=no 
AC_RUN_IFELSE([
      program procedure_test
        implicit none
        type container
           integer,         pointer         :: var 
           procedure(func), pointer, nopass :: f
        end type container
        abstract interface
           logical function func(var)
             integer, pointer :: var 
           end function func
        end interface
        type(container) :: c
        integer         :: var 
        procedure(func) :: test

        call add(var, test)
        if (c%f(c%var)) stop 0
        stop 1
      contains
        subroutine add(var,f)
          integer,         target    :: var 
          procedure(func), optional  :: f
          c%var => var 
          c%f   => f
        end subroutine add 
      end program procedure_test
    
      logical function test(var)
        integer, pointer :: var 
        test = .true.
      end function test],
               [use_f03=yes],
               [use_f03=no])
if test "x${use_f03}" = xyes; then
    AC_SUBST(CPPFLAGS,["$CPPFLAGS -D__F2003"])
    AC_MSG_RESULT([yes])
else
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([Cannot build PPM_RC without Fortran 2003 support])
fi

if test "x${use_mpi}" = xyes; then
    AC_MSG_CHECKING([if your MPI is usable])
    AC_LINK_IFELSE([
                     program conftest_prog
                     include 'mpif.h'
                     call MPI_Init(ierr)
                     end program conftest_prog],
                     [AC_SUBST(CPPFLAGS,["$CPPFLAGS -D__MPI"])],
                     [AC_MSG_ERROR([Cannot build PPM_RC with MPI support like this])])
    AC_MSG_RESULT([yes])
fi

if test "x${use_mpi}" = xyes; then
    AC_MSG_CHECKING([if your MPI supports MPI3 standard])
    AC_LINK_IFELSE([    
		    program conftest_prog
	  	    include 'mpif.h'
	  	    integer :: rank,s,request
	  	    call MPI_Init(ierr)
	  	    call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
	  	    s=rank
	  	    call MPI_Ibcast(s,1,MPI_INTEGER,0,MPI_COMM_WORLD,request,ierr)
	  	    end program conftest_prog],
		    [AC_SUBST(CPPFLAGS,["$CPPFLAGS -D__MPI3"])],
		    [AC_MSG_ERROR([Cannot build parallel PPM_RC without MPI3 support])])
    AC_MSG_RESULT([yes])
fi

AC_LANG([C])
AC_SUBST(MPI_ENABLED,["false"])
if test "x${use_mpi}" = xyes; then
	AC_SUBST(MPI_ENABLED,["true"])
fi

# Checks for libraries and set preprocessing flags
# Make sure that the libraries are compiled with Fortran
AC_LANG([Fortran])
# Check the compiler vendor.. this list should be extended of course
# but let's hope that this will work in most cases
AC_MSG_CHECKING([your Fortran compiler vendor])
vendor=unknown
for arg in --version -V -v ; do
    rm -f conftest.txt
    PPM_RC_RUNLOG([$FC $arg </dev/null >conftest.txt 2>&1])
    # Ignore the return code, because some compilers set the
    # return code to zero on invalid arguments and some to 
    # non-zero on success (with no files to compile)
    if test -f conftest.txt ; then
        if grep 'Portland Group' conftest.txt >/dev/null 2>&1 ; then
            vendor=pgi
        elif grep 'Sun Workshop' conftest.txt >/dev/null 2>&1 ; then
            vendor=sun
	elif grep 'Sun Fortran 9' conftest.txt >/dev/null 2>&1 ; then 
	    vendor=sun
        elif grep 'Absoft' conftest.txt >/dev/null 2>&1 ; then
            vendor=absoft
        elif grep 'G95' conftest.txt >/dev/null 2>&1 ; then
            vendor=gnu
        elif grep 'GNU Fortran' conftest.txt >/dev/null 2>&1 ; then
            # This is gfortran
            vendor=gnu
        elif grep Intel conftest.txt >/dev/null 2>&1 ; then
            vendor=intel
        fi
    fi
    if test "$vendor" != "unknown" ; then
 	    AC_MSG_RESULT([$vendor]) 
	    break ; 
    fi
done
if test "$vendor" = "unknown" ; then
    # Try to use the compiler name
    if test "$FC" = "ifort" -o "$FC" = "ifc" ; then
        vendor=intel
    elif test "$FC" = "pgf90" ; then
        vendor=pgi
    elif test "$FC" = "xlf90" -o "$FC" = "xlf90_r" ; then
        vendor=ibm
    elif test "$FC" = "xlf95" -o "$FC" = "xlf95_r" ; then
        vendor=ibm
    fi
    AC_MSG_RESULT([$vendor])
fi
rm -f conftest.txt
# End of checking for Fortran compiler vendor

AX_F90_MODULE_EXTENSION
if test x$ax_cv_f90_modext = xunknown ; then
	AC_MSG_ERROR([unable to find f90 modules extension])
fi

AC_CANONICAL_HOST

###### Check for debug compilation
AC_MSG_CHECKING(whether to build with debug information)
debuger=no
AC_ARG_ENABLE(debug,
        AC_HELP_STRING(
            [--enable-debug],
            [enable debug data generation (default is no)]
        ),
        debuger="$enableval"
)

AC_MSG_RESULT($debuger)
if test x"$debuger" = x"yes"; then
        AC_DEFINE([DEBUG_MODE],[],[Debug])
        AC_DEFINE([DEBUG],[],[Debug])

    	if test "x${vendor}" = xintel; then
		AC_SUBST(DEBUG,["-g3 -O0 -debug all -traceback -check all -check bounds -check uninit -check pointers -debug all -traceback -check all -check bounds -check uninit -diag-disable 10120"])
                AC_SUBST(CDEBUG,["-g3 -O0 -debug all -traceback"])
	elif test "x${vendor}" = xgnu; then
		AC_SUBST(DEBUG,["-g3 -O0 -Wall -Wextra -fimplicit-none -fbounds-check -fbacktrace -Wsurprising -Waliasing -Wimplicit-interface -pedantic -ffree-form -ffree-line-length-none -ffree-line-length-none -ftrapv -ffpe-trap=invalid,zero,overflow,underflow -Warray-temporaries -fcheck-array-temporaries  "])
                AC_SUBST(CDEBUG,["-g3 -O0 -Wall -Wextra"])
	else
		AC_SUBST(DEBUG,["-g -O0"])
		AC_SUBST(CDEBUG,["-g -O0"])
	fi

	AC_SUBST(CXXFLAGS,["$CXXFLAGS -O0"])
	AC_SUBST(CFLAGS,["$CFLAGS -O0"])
	AC_SUBST(FCFLAGS,["$FCFLAGS -O0"])

	AC_SUBST(CPPFLAGS,["$CPPFLAGS -D__DEBUG"])
else
    	if test "x${vendor}" = xintel; then
    		AC_SUBST(FCFLAGS,["$FCFLAGS -O3 -diag-disable 5462"])
	elif test "x${vendor}" = xgnu; then
		AC_SUBST(FCFLAGS,["$FCFLAGS -O3 -ffree-line-length-none"])
	else
		AC_SUBST(FCFLAGS,["$FCFLAGS -O3"])
	fi

	AC_SUBST(CXXFLAGS,["$CXXFLAGS -O3"])
	AC_SUBST(CFLAGS,["$CFLAGS -O3"])
fi

# Ciheck for PPM 
PPM_LIBRARY_SETUP(ppm,[],ppm_module_finalize.$ax_cv_f90_modext,libppm*,[
			  use ppm_module_finalize
			  integer :: info
			  call ppm_finalize(info)])  
AC_SUBST(PPMDIR)
MODFLAGS=$ax_cv_f90_modflag
AC_SUBST(MODFLAGS)

AC_LANG([C])
# Check for METIS
have_metis=no
AX_METIS([],[echo "Cannot detect metis, use the --with-metis=path option if it is not installed in the default location"
                    exit 201])
have_metis=$HAVE_METIS
if test "x${have_metis}" = xno; then
  (exit 1); exit 1;
else
  AC_SUBST(CPPFLAGS,["$CPPFLAGS -D__METIS"])
  METISDIR=$METIS_DIR
fi
AC_SUBST(METISDIR)
AC_SUBST(LIBS,["$LIBS -lmetis"])
AC_SUBST(LDFLAGS,["$LDFLAGS $METIS_LIB"])

AC_LANG([C++])
##### CHECK FOR BOOST ##############
PPM_RC_BOOST
BOOSTDIR=$ac_boost_dir_custom
AC_SUBST(BOOSTDIR)
AC_SUBST(CXXFLAGS,["$CXXFLAGS -I$BOOSTDIR"])

AC_LANG([C])
##### CHECK FOR LIBTIFF ##########
dnl tiffdir=
dnl AC_ARG_WITH([tiff],
dnl 	AS_HELP_STRING([--with-tiff=directory],
dnl 			[Directory where libtiff is installed]),
dnl 			[tiffdir="$withval"],
dnl 			[])
dnl 
dnl echo $tiffdir
dnl 
dnl 
dnl oldlibs_tiff=$LIBS
dnl LIBS="-L$tiffdir/lib"
AC_CHECK_LIBTIFF([])
dnl LIBS=$oldlibs_tiff
TIFFDIR=$TIFF_DIR
AC_SUBST(TIFFDIR)
AC_ARG_VAR(CLIBS, [libraries to pass to the linker for C compiler, e.g. -l<library>])
AC_SUBST(CLIBS,["$CLIBS -ltiff -lstdc++"])
AC_SUBST(CFLAGS,["$CFLAGS $TIFF_INCLUDE"])
AC_SUBST(LDFLAGS,["$LDFLAGS $TIFF_LIB"])
AC_SUBST(CPPFLAGS,["$CPPFLAGS -D__TIFF"])
dnl AC_SUBST(LDFLAGS,["$LDFLAGS $METIS_LIB"])

# Check to see if linux was enabled with --enable-linux
AC_ARG_ENABLE([linux],[AS_HELP_STRING([--enable-linux],
	      [compile for linux (default is no)])],
	      [AC_SUBST(CPPFLAGS,["$CPPFLAGS -D__Linux"])]
	      )

AC_SUBST(INCLUDES_PATH)

dnl set src as the default srcdir
if test "$srcdir" = "."; then
   srcdir=src
   AC_SUBST(srcdir)
fi

AC_CONFIG_FILES([Makefile settings])
AC_OUTPUT
echo ""
echo  "**************************************"
if [ test x"$have_bigtiff" = x"yes" ]; then
  echo "*    LIBTIFF supports BIGTIFF        *" 
else
  echo "*    LIBTIFF doesn't support BIGTIFF *" 
fi
echo "*                                    *"
if [ test x"$debuger" = x"yes" ]; then
  echo "*    Debug: Yes                      *"
else
  echo "*    Debug: No                       *"
fi
echo   "**************************************"

