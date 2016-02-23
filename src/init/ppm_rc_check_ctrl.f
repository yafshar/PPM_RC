      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_rc_check_ctrl
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Check validity of parameter settings.
      !
      !  Input        :
      !
      !  Input/output :
      !
      !  Output       : info       (I) return status. 0 on success.
      !
      !  Routines     :
      !
      !  Remarks      :
      !
      !  References   :
      !-------------------------------------------------------------------------
      !  MOSAIC Group
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_rc_check_ctrl(info)
        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------

        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER, INTENT(   OUT)  :: info
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        REAL(ppm_kind_double)  :: t0

        INTEGER                :: i

        CHARACTER(LEN=ppm_char) :: caller='ppm_rc_check_ctrl'
        !-------------------------------------------------------------------------
        !  Externals
        !-------------------------------------------------------------------------

        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        !-------------------------------------------------------------------------
        !  Check case name
        !-------------------------------------------------------------------------
        i = LEN_TRIM(casename)
        IF (i.LT.1) THEN
           stdout("CTRL NOTICE: Case name too short. Reset to ppm_rc.")
           casename = 'ppm_rc'
        ENDIF
        IF (i.GT.30) THEN
           stdout("CTRL NOTICE: Case name too long. Truncated to 30 chars.")
           casename = casename(1:30)
        ENDIF

        !-------------------------------------------------------------------------
        !  Check output file name
        !-------------------------------------------------------------------------
        i = LEN_TRIM(outputfile)
        IF (i.LT.1) THEN
           stdout("CTRL NOTICE: Output file name too short. Reset to ppm_rc.")
           outputfile = 'ppm_rc_out'
        ENDIF
        IF (i.GT.30) THEN
           stdout("CTRL NOTICE: Output file name too long.")
        ENDIF

        !-------------------------------------------------------------------------
        !  Check maximum number of iterations
        !-------------------------------------------------------------------------
        IF (maxiter.LT.1) THEN
           stdout("CTRL NOTICE: Maximum number of iterations must be > 0.")
           info = ppm_error_fatal
        ENDIF

        !-------------------------------------------------------------------------
        !  Check input file name
        !-------------------------------------------------------------------------
        i = LEN_TRIM(inputimage)
        IF (i.LT.1) THEN
           stdout("CTRL NOTICE: Image name too short. Reset to ppm_rc.tif")
           inputimage = 'ppm_rc.tif'
        ENDIF
        IF (i.GT.50) THEN
           stdout("CTRL NOTICE: Image name too long.")
        ENDIF

        !-------------------------------------------------------------------------
        !  Check diag file name
        !-------------------------------------------------------------------------
        i = LEN_TRIM(diagfile)
        IF (i.LT.1) THEN
           stdout("CTRL NOTICE: Diag file name too short. Reset to ppm_rc.diag.")
           diagfile = 'ppm_rc_diag'
        ENDIF
        IF (i.GT.30) THEN
           stdout("CTRL NOTICE: Diag file name too long. Truncated to 30.")
           diagfile = diagfile(1:30)
        ENDIF

        !-------------------------------------------------------------------------
        !  Check Abort file name
        !-------------------------------------------------------------------------
        i = LEN_TRIM(abortfile)
        IF (i.LT.1) THEN
           stdout("CTRL NOTICE: Abort file name too short. Reset to ABORT.")
           abortfile = 'ABORT'
        ENDIF
        IF (i.GT.30) THEN
           stdout("CTRL NOTICE: Abort file name too long. Truncated to 30.")
           abortfile = abortfile(1:30)
        ENDIF

        !-------------------------------------------------------------------------
        !  Check output frequency
        !-------------------------------------------------------------------------
        IF (freqoutput.LE.0) THEN
           stdout("CTRL ERROR: Output frequency must be > 0.")
           info = -1
        ENDIF

        !-------------------------------------------------------------------------
        !  Check diagnostics frequency
        !-------------------------------------------------------------------------
        IF (freqdiag.LE.0) THEN
           stdout("CTRL ERROR: Diagnostics frequency must be > 0.")
           info = -1
        ENDIF

!         !-------------------------------------------------------------------------
!         !  Check frame number for validity
!         !-------------------------------------------------------------------------
!         IF (frame.LE.0) THEN
!            stdout("CTRL ERROR: frame number must be > 0.")
!            info = -1
!         ENDIF

        !-------------------------------------------------------------------------
        !  TODO: check that a valid energy has been chosen
        !-------------------------------------------------------------------------

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
      END SUBROUTINE ppm_rc_check_ctrl
