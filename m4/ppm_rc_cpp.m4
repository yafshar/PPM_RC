AC_DEFUN([ACX_GREP_CPP],[
	AC_ARG_VAR(CPP, [Fortran preprocessor])
	echo "$2" > conftest.f
	dnl see if the attempt to preprocess raises an error
	if ! $CPP_base -traditional-cpp -P conftest.f > /dev/null 2>&5 ; then
		$4
	fi
	if (eval "$CPP_base -traditional-cpp -P conftest.f" 2>&5) | grep -q "$1" 2>&5; then :
		$3
	else
		$4
	fi
	rm -f conftest*
])

AC_DEFUN([PPM_RC_CPP],[
	for CPP_base in "$CPP" "`which cpp`"; do
		if test -z "$CPP_base"; then
			continue
		fi
		
		acx_cpp_ok=yes

		ACX_GREP_CPP([anything], AC_LANG_PROGRAM([],[anything]),
			[], [acx_cpp_ok=no; continue])

		ACX_GREP_CPP([hi_s], AC_LANG_PROGRAM([],[
#define  DTYPE(a) a/**/_s
DTYPE(hi)]),
		[], [acx_cpp_ok=no; continue])

		if test x"$acx_cpp_ok" = xyes; then
			AC_MSG_CHECKING([whether $CPP_base is usable for Fortran preprocessing])
			AC_MSG_RESULT([yes])
                        AC_SUBST(CPP,["$CPP_base"])
			break
		fi
	done

	if test x"$acx_cpp_ok" = xno; then
		AC_MSG_ERROR([Could not find preprocessor usable for Fortran! 
		Please set CPP environment variable to the GNU preprocessor. 
		./confgiure CPP=gnucpp 
		The current version of the program for processing Fortran MACROS 
		is only a gnu-compatible. 
		clang preprocessor ond many other ones don't work.])
	fi
])
