      SUBROUTINE uniquers(a,b,info)

        USE ppm_module_util_qsort, ONLY : ppm_util_qsort
        IMPLICIT NONE

        REAL(ppm_kind_single), DIMENSION(:), INTENT(IN   ) :: a
        !input
        INTEGER,               DIMENSION(:), POINTER       :: b
        !unique output
        INTEGER,                             INTENT(  OUT) :: info

        REAL(ppm_kind_double) :: t0

        INTEGER, DIMENSION(SIZE(a))    :: c
        INTEGER, DIMENSION(:), POINTER :: d
        INTEGER, DIMENSION(1)          :: n
        INTEGER                        :: i,j,k,iopt

        CHARACTER(LEN=ppm_char) :: caller ='unique'

        LOGICAL, DIMENSION(SIZE(a)) :: MASK

        CALL substart(caller,t0,info)

        c=INT(a)
        k=SIZE(c,DIM=1)
        IF (k.EQ.0) GOTO 9999

        NULLIFY(d)
        CALL ppm_util_qsort(c,d,info,k)
        or_fail("ppm_util_qsort")

        MASK(d(1))=.TRUE.
        j=1
        DO i=2,k
           IF (c(d(i)).EQ.c(d(i-1))) THEN
              MASK(d(i))=.FALSE.
           ELSE
              MASK(d(i))=.TRUE.
              j=j+1
           ENDIF
        ENDDO

        iopt=ppm_param_dealloc
        CALL ppm_alloc(d,n,iopt,info)
        or_fail_dealloc("unique d")

        n=j
        !COUNT(MASK)

        iopt=ppm_param_alloc_grow
        CALL ppm_alloc(b,n,iopt,info)
        or_fail_alloc("unique b")

        j=1
        DO i=1,k
           IF (MASK(i)) THEN
              b(j)=c(i)
              j=j+1
           ENDIF
        ENDDO
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
      END SUBROUTINE uniquers

      SUBROUTINE uniquer(a,b,info)

        USE ppm_module_util_qsort, ONLY : ppm_util_qsort
        IMPLICIT NONE

        REAL(ppm_kind_double), DIMENSION(:), INTENT(IN   ) :: a
        !input
        INTEGER,               DIMENSION(:), POINTER       :: b
        !unique output
        INTEGER,                             INTENT(  OUT) :: info

        REAL(ppm_kind_double) :: t0

        INTEGER, DIMENSION(SIZE(a))    :: c
        INTEGER, DIMENSION(:), POINTER :: d
        INTEGER, DIMENSION(1)          :: n
        INTEGER                        :: i,j,k,iopt

        CHARACTER(LEN=ppm_char) :: caller ='unique'

        LOGICAL, DIMENSION(SIZE(a)) :: MASK

        CALL substart(caller,t0,info)

        c=INT(a)
        k=SIZE(c,DIM=1)
        IF (k.EQ.0) GOTO 9999

        NULLIFY(d)
        CALL ppm_util_qsort(c,d,info,k)
        or_fail("ppm_util_qsort")

        MASK(d(1))=.TRUE.
        j=1
        DO i=2,k
           IF (c(d(i)).EQ.c(d(i-1))) THEN
              MASK(d(i))=.FALSE.
           ELSE
              MASK(d(i))=.TRUE.
              j=j+1
           ENDIF
        ENDDO

        iopt=ppm_param_dealloc
        CALL ppm_alloc(d,n,iopt,info)
        or_fail_dealloc("unique d")

        n = j
        !COUNT(MASK)

        iopt=ppm_param_alloc_grow
        CALL ppm_alloc(b,n,iopt,info)
        or_fail_alloc("unique b")

        j=1
        DO i=1,k
           IF (MASK(i)) THEN
              b(j)=c(i)
              j=j+1
           ENDIF
        ENDDO
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
      END SUBROUTINE uniquer

      SUBROUTINE uniquei(a,b,info)

        USE ppm_module_util_qsort, ONLY : ppm_util_qsort
        IMPLICIT NONE

        INTEGER, DIMENSION(:), INTENT(IN   ) :: a
        !input
        INTEGER, DIMENSION(:), POINTER       :: b
        !unique output
        INTEGER,               INTENT(  OUT) :: info

        REAL(ppm_kind_double) :: t0

        INTEGER, DIMENSION(:), POINTER :: d
        INTEGER, DIMENSION(1)          :: n
        INTEGER                        :: i,j,k,iopt

        CHARACTER(LEN=ppm_char) :: caller ='unique'

        LOGICAL, DIMENSION(SIZE(a)) :: MASK

        CALL substart(caller,t0,info)

        k=SIZE(a,DIM=1)
        IF (k.EQ.0) GOTO 9999

        NULLIFY(d)
        CALL ppm_util_qsort(a,d,info,k)
        or_fail("ppm_util_qsort")

        MASK(d(1))=.TRUE.
        j=1
        DO i=2,k
           IF (a(d(i)).EQ.a(d(i-1))) THEN
              MASK(d(i))=.FALSE.
           ELSE
              MASK(d(i))=.TRUE.
              j=j+1
           ENDIF
        ENDDO

        iopt=ppm_param_dealloc
        CALL ppm_alloc(d,n,iopt,info)
        or_fail_dealloc("unique d")

        n=j
        !COUNT(MASK)

        iopt=ppm_param_alloc_grow
        CALL ppm_alloc(b,n,iopt,info)
        or_fail_alloc("unique b")

        j=1
        DO i=1,k
           IF (MASK(i)) THEN
              b(j)=a(i)
              j=j+1
           ENDIF
        ENDDO
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
      END SUBROUTINE uniquei

      SUBROUTINE uniqueli(a,b,info)

        USE ppm_module_util_qsort, ONLY : ppm_util_qsort
        IMPLICIT NONE

        INTEGER(ppm_kind_int64), DIMENSION(:), INTENT(IN   ) :: a
        !input
        INTEGER(ppm_kind_int64), DIMENSION(:), POINTER       :: b
        !unique output
        INTEGER,                               INTENT(  OUT) :: info

        REAL(ppm_kind_double) :: t0

        INTEGER(ppm_kind_int64), DIMENSION(SIZE(a))    :: c
        INTEGER,                 DIMENSION(:), POINTER :: d
        INTEGER,                 DIMENSION(1)          :: n
        INTEGER                                        :: i,j,k,iopt

        CHARACTER(LEN=ppm_char) :: caller ='unique'

        LOGICAL, DIMENSION(SIZE(a)) :: MASK

        CALL substart(caller,t0,info)

        c=a
        k=SIZE(c,DIM=1)
        IF (k.EQ.0) GOTO 9999

        NULLIFY(d)
        CALL ppm_util_qsort(c,d,info,k)
        or_fail("ppm_util_qsort")

        MASK(d(1))=.TRUE.
        j=1
        DO i=2,k
           IF (c(d(i)).EQ.c(d(i-1))) THEN
              MASK(d(i))=.FALSE.
           ELSE
              MASK(d(i))=.TRUE.
              j=j+1
           ENDIF
        ENDDO

        iopt=ppm_param_dealloc
        CALL ppm_alloc(d,n,iopt,info)
        or_fail_dealloc("unique d")

        n=j
        !COUNT(MASK)

        iopt=ppm_param_alloc_grow
        CALL ppm_alloc(b,n,iopt,info)
        or_fail_alloc("unique b")

        j=1
        DO i=1,k
           IF (MASK(i)) THEN
              b(j)=c(i)
              j=j+1
           ENDIF
        ENDDO
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
      END SUBROUTINE uniqueli

