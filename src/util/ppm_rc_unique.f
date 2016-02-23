      SUBROUTINE uniquers(a,b,info)

        USE ppm_module_util_qsort, ONLY : ppm_util_qsort
        IMPLICIT NONE

        REAL(ppm_kind_single), DIMENSION(:), INTENT(IN   ) :: a
        !input
        INTEGER,               DIMENSION(:), POINTER       :: b
        !unique output
        INTEGER,                             INTENT(  OUT) :: info

        REAL(ppm_kind_double) :: t0

        INTEGER, DIMENSION(SIZE(a)) :: c
        INTEGER, DIMENSION(1)       :: n
        INTEGER                     :: i,j,iopt

        CHARACTER(LEN=ppm_char) :: caller ='unique'

        LOGICAL, DIMENSION(SIZE(a)) :: MASK

        CALL substart(caller,t0,info)

        c=INT(a)

        n = SIZE(c)
        IF (n(1).EQ.0) GOTO 9999

        CALL ppm_util_qsort(c,info,n(1))
        or_fail("ppm_util_qsort")

!         MASK=.NOT.[(ANY(c(i).EQ.c(1:i-1)),i=1,n(1))]

        MASK(1)=.TRUE.
        j=1
        DO i=2,n(1)
           IF (c(i).EQ.c(i-1)) THEN
              MASK(i)=.FALSE.
           ELSE
              MASK(i)=.TRUE.
              j=j+1
           ENDIF
        ENDDO

        n = j !COUNT(MASK)

        iopt=ppm_param_alloc_grow
        CALL ppm_alloc(b,n,iopt,info)
        or_fail_alloc("unique b")

        b=PACK(c,MASK)
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

        INTEGER, DIMENSION(SIZE(a)) :: c
        INTEGER, DIMENSION(1)       :: n
        INTEGER                     :: i,j,iopt

        CHARACTER(LEN=ppm_char) :: caller ='unique'

        LOGICAL, DIMENSION(SIZE(a)) :: MASK

        CALL substart(caller,t0,info)

        c=INT(a)

        n = SIZE(c)
        IF (n(1).EQ.0) GOTO 9999

        CALL ppm_util_qsort(c,info,n(1))
        or_fail("ppm_util_qsort")

!         MASK=.NOT.[(ANY(c(i).EQ.c(1:i-1)),i=1,n(1))]

        MASK(1)=.TRUE.
        j=1
        DO i=2,n(1)
           IF (c(i).EQ.c(i-1)) THEN
              MASK(i)=.FALSE.
           ELSE
              MASK(i)=.TRUE.
              j=j+1
           ENDIF
        ENDDO

        n = j!COUNT(MASK)

        iopt=ppm_param_alloc_grow
        CALL ppm_alloc(b,n,iopt,info)
        or_fail_alloc("unique b")

        b=PACK(c,MASK)
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

        INTEGER, DIMENSION(SIZE(a)) :: c
        INTEGER, DIMENSION(1)       :: n
        INTEGER                     :: i,j,iopt

        CHARACTER(LEN=ppm_char) :: caller ='unique'

        LOGICAL, DIMENSION(SIZE(a)) :: MASK

        CALL substart(caller,t0,info)

        c=a

        n = SIZE(c)
        IF (n(1).EQ.0) GOTO 9999

        CALL ppm_util_qsort(c,info,n(1))
        or_fail("ppm_util_qsort")

!         MASK=.NOT.[(ANY(a(i).EQ.a(1:i-1)),i=1,n(1))]

        MASK(1)=.TRUE.
        j=1
        DO i=2,n(1)
           IF (c(i).EQ.c(i-1)) THEN
              MASK(i)=.FALSE.
           ELSE
              MASK(i)=.TRUE.
              j=j+1
           ENDIF
        ENDDO

        n = j!COUNT(MASK)

        iopt=ppm_param_alloc_grow
        CALL ppm_alloc(b,n,iopt,info)
        or_fail_alloc("unique b")

        b=PACK(c,MASK)
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

        INTEGER(ppm_kind_int64), DIMENSION(SIZE(a)) :: c
        INTEGER,                 DIMENSION(1)       :: n
        INTEGER                                     :: i,j,iopt

        CHARACTER(LEN=ppm_char) :: caller ='unique'

        LOGICAL, DIMENSION(SIZE(a)) :: MASK

        CALL substart(caller,t0,info)

        c=a

        n = SIZE(c)
        IF (n(1).EQ.0) GOTO 9999

        CALL ppm_util_qsort(c,info,n(1))
        or_fail("ppm_util_qsort")

!         MASK=.NOT.[(ANY(a(i).EQ.a(1:i-1)),i=1,n(1))]
        MASK(1)=.TRUE.
        j=1
        DO i=2,n(1)
           IF (c(i).EQ.c(i-1)) THEN
              MASK(i)=.FALSE.
           ELSE
              MASK(i)=.TRUE.
              j=j+1
           ENDIF
        ENDDO

        n = j!COUNT(MASK)

        iopt=ppm_param_alloc_grow
        CALL ppm_alloc(b,n,iopt,info)
        or_fail_alloc("unique b")

        b=PACK(a,MASK)
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
      END SUBROUTINE uniqueli

