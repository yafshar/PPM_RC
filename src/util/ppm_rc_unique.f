      !!! [NOTE]
      !!! These routines are not general, they are provided for the especial
      !!! purpose to work on the input array of labels
      !!!
      !!! The unique routines here are for the purpose of unique array of
      !!! Labels which can be in the form of single or double or integer or
      !!! long integer

      SUBROUTINE ppm_rc_uniquers(a,b,info)

        USE ppm_module_util_qsort, ONLY : ppm_util_qsort
        IMPLICIT NONE

        REAL(ppm_kind_single), DIMENSION(:), INTENT(IN   ) :: a
        !input
        INTEGER,               DIMENSION(:), POINTER       :: b
        !unique output
        INTEGER,                             INTENT(  OUT) :: info

        REAL(ppm_kind_double) :: t0

        INTEGER, DIMENSION(SIZE(a))    :: c
        INTEGER, DIMENSION(SIZE(a))    :: MASKc
        INTEGER, DIMENSION(:), POINTER :: d
        INTEGER, DIMENSION(1)          :: n
        INTEGER                        :: i,j,k,iopt

        CHARACTER(LEN=ppm_char) :: caller ='unique'

        CALL substart(caller,t0,info)

        c=INT(a)
        k=SIZE(c,DIM=1)
        IF (k.EQ.0) GOTO 9999

        NULLIFY(d)
        CALL ppm_util_qsort(c,d,info,k)
        or_fail("ppm_util_qsort")

        j=1
        MASKc(1)=c(d(1))
        DO i=2,k
           IF (c(d(i)).NE.c(d(i-1))) THEN
              j=j+1
              MASKc(j)=c(d(i))
           ENDIF
        ENDDO

        iopt=ppm_param_dealloc
        CALL ppm_alloc(d,n,iopt,info)
        or_fail_dealloc("unique d")

        n=j

        iopt=ppm_param_alloc_grow
        CALL ppm_alloc(b,n,iopt,info)
        or_fail_alloc("unique b")

        DO i=1,j
           b(i)=MASKc(i)
        ENDDO
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
      END SUBROUTINE ppm_rc_uniquers

      SUBROUTINE ppm_rc_uniquer(a,b,info)

        USE ppm_module_util_qsort, ONLY : ppm_util_qsort
        IMPLICIT NONE

        REAL(ppm_kind_double), DIMENSION(:), INTENT(IN   ) :: a
        !input
        INTEGER,               DIMENSION(:), POINTER       :: b
        !unique output
        INTEGER,                             INTENT(  OUT) :: info

        REAL(ppm_kind_double) :: t0

        INTEGER, DIMENSION(SIZE(a))    :: c
        INTEGER, DIMENSION(SIZE(a))    :: MASKc
        INTEGER, DIMENSION(:), POINTER :: d
        INTEGER, DIMENSION(1)          :: n
        INTEGER                        :: i,j,k,iopt

        CHARACTER(LEN=ppm_char) :: caller ='unique'

        CALL substart(caller,t0,info)

        c=INT(a)
        k=SIZE(c,DIM=1)
        IF (k.EQ.0) GOTO 9999

        NULLIFY(d)
        CALL ppm_util_qsort(c,d,info,k)
        or_fail("ppm_util_qsort")

        j=1
        MASKc(1)=c(d(1))
        DO i=2,k
           IF (c(d(i)).NE.c(d(i-1))) THEN
              j=j+1
              MASKc(j)=c(d(i))
           ENDIF
        ENDDO

        iopt=ppm_param_dealloc
        CALL ppm_alloc(d,n,iopt,info)
        or_fail_dealloc("unique d")

        n=j

        iopt=ppm_param_alloc_grow
        CALL ppm_alloc(b,n,iopt,info)
        or_fail_alloc("unique b")

        DO i=1,j
           b(i)=MASKc(i)
        ENDDO
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
      END SUBROUTINE ppm_rc_uniquer

      SUBROUTINE ppm_rc_uniquei(a,b,info)

        USE ppm_module_util_qsort, ONLY : ppm_util_qsort
        IMPLICIT NONE

        INTEGER, DIMENSION(:), INTENT(IN   ) :: a
        !input
        INTEGER, DIMENSION(:), POINTER       :: b
        !unique output
        INTEGER,               INTENT(  OUT) :: info

        REAL(ppm_kind_double) :: t0

        INTEGER, DIMENSION(SIZE(a))    :: MASKc
        INTEGER, DIMENSION(:), POINTER :: d
        INTEGER, DIMENSION(1)          :: n
        INTEGER                        :: i,j,k,iopt

        CHARACTER(LEN=ppm_char) :: caller ='unique'

        CALL substart(caller,t0,info)

        k=SIZE(a,DIM=1)
        IF (k.EQ.0) GOTO 9999

        NULLIFY(d)
        CALL ppm_util_qsort(a,d,info,k)
        or_fail("ppm_util_qsort")

        j=1
        MASKc(1)=a(d(1))
        DO i=2,k
           IF (a(d(i)).NE.a(d(i-1))) THEN
              j=j+1
              MASKc(j)=a(d(i))
           ENDIF
        ENDDO

        iopt=ppm_param_dealloc
        CALL ppm_alloc(d,n,iopt,info)
        or_fail_dealloc("unique d")

        n=j

        iopt=ppm_param_alloc_grow
        CALL ppm_alloc(b,n,iopt,info)
        or_fail_alloc("unique b")

        DO i=1,j
           b(i)=MASKc(i)
        ENDDO
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
      END SUBROUTINE ppm_rc_uniquei

      SUBROUTINE ppm_rc_uniqueli(a,b,info)

        USE ppm_module_util_qsort, ONLY : ppm_util_qsort
        IMPLICIT NONE

        INTEGER(ppm_kind_int64), DIMENSION(:), INTENT(IN   ) :: a
        !input
        INTEGER(ppm_kind_int64), DIMENSION(:), POINTER       :: b
        !unique output
        INTEGER,                               INTENT(  OUT) :: info

        REAL(ppm_kind_double) :: t0

        INTEGER(ppm_kind_int64), DIMENSION(SIZE(a))    :: MASKc
        INTEGER,                 DIMENSION(:), POINTER :: d
        INTEGER,                 DIMENSION(1)          :: n
        INTEGER                                        :: i,j,k,iopt

        CHARACTER(LEN=ppm_char) :: caller ='unique'

        CALL substart(caller,t0,info)

        k=SIZE(a,DIM=1)
        IF (k.EQ.0) GOTO 9999

        NULLIFY(d)
        CALL ppm_util_qsort(a,d,info,k)
        or_fail("ppm_util_qsort")

        j=1
        MASKc(1)=a(d(1))
        DO i=2,k
           IF (a(d(i)).NE.a(d(i-1))) THEN
              j=j+1
              MASKc(j)=a(d(i))
           ENDIF
        ENDDO

        iopt=ppm_param_dealloc
        CALL ppm_alloc(d,n,iopt,info)
        or_fail_dealloc("unique d")

        n=j

        iopt=ppm_param_alloc_grow
        CALL ppm_alloc(b,n,iopt,info)
        or_fail_alloc("unique b")

        DO i=1,j
           b(i)=MASKc(i)
        ENDDO
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
      END SUBROUTINE ppm_rc_uniqueli

