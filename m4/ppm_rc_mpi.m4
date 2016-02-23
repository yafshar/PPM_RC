dnl ===========================================================================
dnl                http://autoconf-archive.cryp.to/acx_mpi.html
dnl ===========================================================================
dnl
dnl SYNOPSIS
dnl
dnl   PPM_RC_MPI([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl DESCRIPTION
dnl
dnl   This macro tries to find out how to compile programs that use MPI
dnl   (Message Passing Interface), a standard API for parallel process
dnl   communication (see http://www-unix.mcs.anl.gov/mpi/)
dnl
dnl   On success, it sets the MPICC, MPICXX, MPIF77, or MPIFC output variable
dnl   to the name of the MPI compiler, depending upon the current language.
dnl   (This may just be $CC/$CXX/$F77/$FC, but is more often something like
dnl   mpicc/mpiCC/mpif77/mpif90.) It also sets MPILIBS to any libraries that
dnl   are needed for linking MPI (e.g. -lmpi or -lfmpi, if a special
dnl   MPICC/MPICXX/MPIF77/MPIFC was not found).
dnl
dnl   If you want to compile everything with MPI, you should set:
dnl
dnl       CC="MPICC" #OR# CXX="MPICXX" #OR# F77="MPIF77" #OR# FC="MPIFC"
dnl       LIBS="$MPILIBS $LIBS"
dnl
dnl   NOTE: The above assumes that you will use $CC (or whatever) for linking
dnl   as well as for compiling. (This is the default for automake and most
dnl   Makefiles.)
dnl
dnl   The user can force a particular library/compiler by setting the
dnl   MPICC/MPICXX/MPIF77/MPIFC and/or MPILIBS environment variables.
dnl
dnl   ACTION-IF-FOUND is a list of shell commands to run if an MPI library is
dnl   found, and ACTION-IF-NOT-FOUND is a list of commands to run if it is not
dnl   found. If ACTION-IF-FOUND is not specified, the default action will
dnl   define HAVE_MPI.
dnl
dnl LAST MODIFICATION
dnl
dnl   2008-04-12
dnl
dnl COPYLEFT
dnl
dnl   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
dnl   Copyright (c) 2008 Julian C. Cummings <cummings@cacr.caltech.edu>
dnl
dnl   This program is free software: you can redistribute it and/or modify it
dnl   under the terms of the GNU General Public License as published by the
dnl   Free Software Foundation, either version 3 of the License, or (at your
dnl   option) any later version.
dnl
dnl   This program is distributed in the hope that it will be useful, but
dnl   WITHOUT ANY WARRANTY; without even the implied warranty of
dnl   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
dnl   Public License for more details.
dnl
dnl   You should have received a copy of the GNU General Public License along
dnl   with this program. If not, see <http://www.gnu.org/licenses/>.
dnl
dnl   As a special exception, the respective Autoconf Macro's copyright owner
dnl   gives unlimited permission to copy, distribute and modify the configure
dnl   scripts that are the output of Autoconf when processing the Macro. You
dnl   need not follow the terms of the GNU General Public License when using
dnl   or distributing such scripts, even though portions of the text of the
dnl   Macro appear in them. The GNU General Public License (GPL) does govern
dnl   all other use of the material that constitutes the Autoconf Macro.
dnl
dnl   This special exception to the GPL applies to versions of the Autoconf
dnl   Macro released by the Autoconf Macro Archive. When you make and
dnl   distribute a modified version of the Autoconf Macro, you may extend this
dnl   special exception to the GPL to apply to your modified version as well.

dnl   Modifications & editions by Yaser Afshar (afshar@mpi-cbg.de)

AC_DEFUN([PPM_RC_MPI], [
AC_PREREQ(2.50) dnl for AC_LANG_CASE

AC_REQUIRE([AC_PROG_CC])
dnl AC_ARG_VAR(MPICC,[MPI C compiler command])
AC_CHECK_PROGS(MPICC, mpicc hcc mpxlc_r mpxlc mpcc cmpicc, [no])
AS_VAR_IF(MPICC,[no],[AC_MSG_ERROR([Could not find mpicc for MPI])],[CC="$MPICC"])
CC="$MPICC"
dnl AC_SUBST(MPICC)

AC_REQUIRE([AC_PROG_CXX])
dnl AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
AC_CHECK_PROGS(MPICXX, mpic++ mpicxx mpiCC hcp mpxlC_r mpxlC mpCC cmpic++, [no])
AS_VAR_IF(MPICXX,[no],[AC_MSG_ERROR([Could not find mpic++ for MPI])],[CXX="$MPICXX"])
CXX="$MPICXX"
dnl AC_SUBST(MPICXX)

if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpi.h])
	AC_TRY_COMPILE([#include <mpi.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
		AC_MSG_RESULT(no)])
fi

AC_REQUIRE([AC_PROG_FC])
dnl AC_ARG_VAR(MPIFC,[MPI Fortran compiler command])
AC_CHECK_PROGS(MPIFC, mpifort mpif90 mpxlf95_r mpxlf90_r mpxlf95 mpxlf90 mpf90 cmpif90c, [no])
AS_VAR_IF(MPIFC,[no],[AC_MSG_ERROR([Could not find mpif90 for MPI])],[FC="$MPIFC"])
FC="$MPIFC"
dnl AC_SUBST(MPIFC)

if test x = x"$MPILIBS"; then
        AC_LANG_CASE([C], [AC_CHECK_FUNC(MPI_Iallgather, [MPILIBS=" "])],
                [C++], [AC_CHECK_FUNC(MPI_Iallgather, [MPILIBS=" "])],
                [Fortran], [AC_MSG_CHECKING([for MPI_Iallgather])
                        AC_TRY_LINK([],[      call MPI_Iallgather], [MPILIBS=" "
                                AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)])])
fi

if test x = x"$MPILIBS"; then
	AC_CHECK_LIB(mpi, MPI_Iallgather, [MPILIBS="-lmpi"])
fi

if test x = x"$MPILIBS"; then
	AC_CHECK_LIB(fmpi, MPI_Iallgather, [MPILIBS="-lfmpi"])
fi

if test x = x"$MPILIBS"; then
	AC_CHECK_LIB(mpichf90, MPI_Iallgather, [MPILIBS="-lmpichf90"])
fi

if test x = x"$MPILIBS"; then
	AC_CHECK_LIB(mpich, MPI_Iallgather, [MPILIBS="-lmpich"])
fi

AC_SUBST(MPILIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x = x"$MPILIBS"; then
        $2
        :
else
        ifelse([$1],,[AC_DEFINE(HAVE_MPI,1,[Define if you have the MPI library.])],[$1])
        :
fi
])dnl PPM_RC_MPI

AC_DEFUN([PPM_RC_NOMPI], [

AS_VAR_IF(CC,[mpicc],[PPM_RC_MPI use_mpi=yes],[])
if [ test x"$use_mpi" = x"no" ]; then
	AS_VAR_IF(CC,[hcc],[PPM_RC_MPI use_mpi=yes],[])
fi
if [ test x"$use_mpi" = x"no" ]; then
	AS_VAR_IF(CC,[mpxlc_r],[PPM_RC_MPI use_mpi=yes],[])
fi
if [ test x"$use_mpi" = x"no" ]; then
	AS_VAR_IF(CC,[mpxlc],[PPM_RC_MPI use_mpi=yes],[])
fi
if [ test x"$use_mpi" = x"no" ]; then
        AS_VAR_IF(CC,[mpcc],[PPM_RC_MPI use_mpi=yes],[])
fi
if [ test x"$use_mpi" = x"no" ]; then
        AS_VAR_IF(CC,[cmpicc],[PPM_RC_MPI use_mpi=yes],[])
fi
if [ test x"$use_mpi" = x"no" ]; then
	AS_VAR_IF(CXX,[mpic++],[PPM_RC_MPI use_mpi=yes],[])
fi
if [ test x"$use_mpi" = x"no" ]; then
	AS_VAR_IF(CXX,[mpicxx],[PPM_RC_MPI use_mpi=yes],[])
fi
if [ test x"$use_mpi" = x"no" ]; then
	AS_VAR_IF(CXX,[mpiCC],[PPM_RC_MPI use_mpi=yes],[])
fi
if [ test x"$use_mpi" = x"no" ]; then
        AS_VAR_IF(CXX,[hcp],[PPM_RC_MPI use_mpi=yes],[])
fi
if [ test x"$use_mpi" = x"no" ]; then
        AS_VAR_IF(CXX,[mpxlC_r],[PPM_RC_MPI use_mpi=yes],[])
fi
if [ test x"$use_mpi" = x"no" ]; then
        AS_VAR_IF(CXX,[mpxlC],[PPM_RC_MPI use_mpi=yes],[])
fi
if [ test x"$use_mpi" = x"no" ]; then
        AS_VAR_IF(CXX,[mpCC],[PPM_RC_MPI use_mpi=yes],[])
fi
if [ test x"$use_mpi" = x"no" ]; then
        AS_VAR_IF(CXX,[cmpic++],[PPM_RC_MPI use_mpi=yes],[])
fi
if [ test x"$use_mpi" = x"no" ]; then
        AS_VAR_IF(FC,[mpifort],[PPM_RC_MPI use_mpi=yes],[])
fi
if [ test x"$use_mpi" = x"no" ]; then
        AS_VAR_IF(FC,[mpif90],[PPM_RC_MPI use_mpi=yes],[])
fi
if [ test x"$use_mpi" = x"no" ]; then
        AS_VAR_IF(FC,[mpxlf95_r],[PPM_RC_MPI use_mpi=yes],[])
fi
if [ test x"$use_mpi" = x"no" ]; then
        AS_VAR_IF(FC,[mpxlf90_r],[PPM_RC_MPI use_mpi=yes],[])
fi
if [ test x"$use_mpi" = x"no" ]; then
        AS_VAR_IF(FC,[mpxlf95],[PPM_RC_MPI use_mpi=yes],[])
fi
if [ test x"$use_mpi" = x"no" ]; then
        AS_VAR_IF(FC,[mpxlf90],[PPM_RC_MPI use_mpi=yes],[])
fi
if [ test x"$use_mpi" = x"no" ]; then
        AS_VAR_IF(FC,[mpf90],[PPM_RC_MPI use_mpi=yes],[])
fi
if [ test x"$use_mpi" = x"no" ]; then
        AS_VAR_IF(FC,[cmpif90c],[PPM_RC_MPI use_mpi=yes],[])
fi
])
