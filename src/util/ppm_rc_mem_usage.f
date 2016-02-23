      SUBROUTINE ppm_rc_mem_usage(valueRSS,info)

!         USE ifport !if on intel compiler
        IMPLICIT NONE

        INTEGER, INTENT(INOUT) :: valueRSS
        !input
        INTEGER, INTENT(  OUT) :: info

        REAL(ppm_kind_double) :: t0

        INTEGER :: gpid

        CHARACTER(LEN=30)       :: dummy
        CHARACTER(LEN=ppm_char) :: rss_use_name,command
        CHARACTER(LEN=ppm_char) :: caller ='ppm_rc_mem_usage'

        CALL substart(caller,t0,info)

        gpid=GETPID()

        WRITE(dummy,'(I0)') gpid
        WRITE(rss_use_name,'(A,I0)')"/tmp/rss_use_",rank

        WRITE(command,'(4A)') "cat /proc/",TRIM(ADJUSTL(dummy)),"/status | grep RSS | grep -Eo '[0-9]{0,9}' > ",TRIM(ADJUSTL(rss_use_name))

        info=SYSTEM(command)

        OPEN(UNIT=777,FILE=TRIM(ADJUSTL(rss_use_name)))
        READ(777,FMT=*,END=7777,ERR=8888) valueRSS
        CLOSE(777,STATUS='DELETE')

        valueRSS=valueRSS/1024-valueRSS0
        !Convert memory use in KB to MB

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        GOTO 9999
      7777 CLOSE(777,STATUS='DELETE')
        valueRSS=-1
        GOTO 9999
      8888 CLOSE(777,STATUS='DELETE')
        valueRSS=-2
      9999 CONTINUE
        CALL substop(caller,t0,info)
      END SUBROUTINE ppm_rc_mem_usage
