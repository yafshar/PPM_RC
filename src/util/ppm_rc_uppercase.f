      !    %----------------------------------------------------------
      !    | Converts the string "s1" to upper case
      !    %----------------------------------------------------------

      SUBROUTINE ppm_rc_uppercase(record)
        !    %-----------------------------------------------------------------%
        !    | Author             - y.afshar         January     2010
        !    %-----------------------------------------------------------------%
        IMPLICIT NONE

        CHARACTER(LEN = *), INTENT(INOUT) :: record

        INTEGER :: length, LowerToUpper, i

        CHARACTER(LEN = LEN(record)) :: temp

        LowerToUpper = IACHAR("A") - IACHAR("a")

        temp = ADJUSTL(record)
        length = LEN_TRIM(temp)

        record = ''
        DO i=1,length
           SELECT CASE (temp(i:i))
           CASE ('a' : 'z')
              record(i:i) = ACHAR(IACHAR(temp(i:i)) + LowerToUpper)

           CASE DEFAULT
              record(i:i) = temp(i:i)

           END SELECT
        ENDDO
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END SUBROUTINE ppm_rc_uppercase
