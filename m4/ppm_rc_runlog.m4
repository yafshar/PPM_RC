dnl
dnl PPM_RC_RUN_LOG mimics _AC_RUN_LOG which is autoconf internal routine.
dnl
AC_DEFUN([PPM_RC_RUNLOG],[
{ AS_ECHO(["$as_me:$LINENO: $1"]) >&AS_MESSAGE_LOG_FD
  (eval $1) 2>&AS_MESSAGE_LOG_FD
  ac_status=$?
  AS_ECHO(["$as_me:$LINENO: \$? = $ac_status"]) >&AS_MESSAGE_LOG_FD
  test $ac_status = 0; }])
