
      SUBROUTINE compute_num_common(a,b,c,n,info)

        IMPLICIT NONE

        INTEGER, DIMENSION(:), POINTER       :: a
        !first input
        INTEGER, DIMENSION(:), POINTER       :: b
        !second input
        INTEGER, DIMENSION(:), POINTER       :: c
        !common output
        INTEGER,               INTENT(  OUT) :: n
        !number of common
        INTEGER,               INTENT(  OUT) :: info

        REAL(ppm_kind_double) :: t0

        INTEGER :: i,iopt

        CHARACTER(LEN=ppm_char) :: caller='compute_num_common'

        LOGICAL, DIMENSION(SIZE(a)) :: MASK

        CALL substart(caller,t0,info)

        MASK=[(ANY(a(i).EQ.b),i=1,SIZE(a))]

        n = COUNT(MASK)

        iopt=ppm_param_alloc_grow
        CALL ppm_alloc(c,(/n/),iopt,info)
        or_fail_alloc("compute_num_common c")

        c=PACK(a,MASK)
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
      END SUBROUTINE compute_num_common

