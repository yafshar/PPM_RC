      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_rc_dump_parameters
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Prints some of the parameters to stdout as well as
      !                 the log file.
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
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_rc_dump_parameters(info)

        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------
        USE ppm_rc_module_util, ONLY : ppm_rc_uppercase
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER, INTENT(   OUT) :: info
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        REAL(ppm_kind_double) :: t0

        INTEGER, DIMENSION(8) :: tval

        CHARACTER(LEN = ppm_char) :: caller='ppm_rc_dump_parameters'

        !-------------------------------------------------------------------------
        !  Externals
        !-------------------------------------------------------------------------

        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        !-------------------------------------------------------------------------
        !  Root only
        !-------------------------------------------------------------------------
        IF (rank.NE.0) GOTO 9999

        SELECT CASE (UseMCMC)
        CASE (.FALSE.)
           stdout_f('(A)',"---- PPM_RC REGION COMPETITION SEGMENTATION -------")
        CASE (.TRUE.)
           stdout_f('(A)',"---- PPM_RC DISCRETE REGION SAMPLING -------")
        END SELECT
        CALL ppm_log(caller,cbuf,info)

        !----------------------------------------------------------------------
        !  Get time and date
        !----------------------------------------------------------------------
        CALL DATE_AND_TIME(values=tval)
        stdout_f('(A,I4.4,A,5(I2.2,A))',"   Date and Time of Start:   ", &
        & 'tval(1)',"-",'tval(2)',"-",'tval(3)'," ",'tval(5)',":",       &
        & 'tval(6)',":",'tval(7)'," ")
        CALL ppm_log(caller,cbuf,info)

        !----------------------------------------------------------------------
        !  Control and debug (command line parameters)
        !----------------------------------------------------------------------
        stdout_f('(2A)',"   Control file read:        ",'TRIM(ctrl_file_name)')
        CALL ppm_log(caller,cbuf,info)

        stdout_f('(A,I1.1)',"   Debug level for this run: ",debug)
        CALL ppm_log(caller,cbuf,info)

        !----------------------------------------------------------------------
        !  Run case name
        !----------------------------------------------------------------------
        stdout_f('(A)',"---- CASE DETAILS ---------------------------------")
        CALL ppm_log(caller,cbuf,info)

        stdout_f('(2A)',"   This is problem case      ",'TRIM(casename)')
        CALL ppm_log(caller,cbuf,info)

        !----------------------------------------------------------------------
        !  Restart or input file
        !----------------------------------------------------------------------
        stdout_f('(2A)',"   Image file:               ",'TRIM(inputimage)')
        CALL ppm_log(caller,cbuf,info)

        stdout_f('(A,I0)',"   Case dimensions:          ", ppm_rc_dim)
        CALL ppm_log(caller,cbuf,info)

        stdout_f('(A,3I6)',"   Number of pixels:          ",'Ngrid(1:ppm_rc_dim)')
        CALL ppm_log(caller,cbuf,info)

        stdout_f('(A,I8)',"   Max number of iterations:  ",maxiter)
        CALL ppm_log(caller,cbuf,info)

        !----------------------------------------------------------------------
        !  Domain decomposition
        !----------------------------------------------------------------------
        stdout_f('(A)',"---- DOMAIN DECOMPOSITION -------------------------")
        CALL ppm_log(caller,cbuf,info)

        stdout_f('(2A)',"   Domain decomposition used: bisection")
        CALL ppm_log(caller,cbuf,info)

!         !----------------------------------------------------------------------
!         !  Diagnostics
!         !----------------------------------------------------------------------
!         stdout_f('(A)',"---- DIAGNOSTIC OUTPUT ----------------------------")
!         CALL ppm_log(caller,cbuf,info)
!
!         IF (LEN_TRIM(diagfile).LT.1) THEN
!            stdout_f('(2A)',"   Diagnostics file:         ","disabled")
!            CALL ppm_log(caller,cbuf,info)
!         ELSE
!            stdout_f('(A,I6)',"   Time steps between diag:  ",freqdiag)
!            CALL ppm_log(caller,cbuf,info)
!
!            stdout_f('(2A)',"   Diagnostics file:         ",'TRIM(diagfile)')
!            CALL ppm_log(caller,cbuf,info)
!         ENDIF
        !----------------------------------------------------------------------
        !  Output
        !----------------------------------------------------------------------
        stdout_f('(A)',"---- RESULT DATA OUTPUT ---------------------------")
        CALL ppm_log(caller,cbuf,info)

        stdout_f('(A,I6)',"   Time steps between output:",freqoutput)
        CALL ppm_log(caller,cbuf,info)

        stdout_f('(2A)',"   Output file name:         ",'TRIM(outputfile)')
        CALL ppm_log(caller,TRIM(cbuf),info)

        !----------------------------------------------------------------------
        !  Flags
        !----------------------------------------------------------------------
        cbuf=MERGE('   Probe processor speeds:   yes','   Probe processor speeds:   no ',probeproc)
        stdout_f('(A)','TRIM(cbuf)')
        CALL ppm_log(caller,cbuf,info)

        !----------------------------------------------------------------------
        !  Initialization kind
        !----------------------------------------------------------------------
        stdout_f('(A)',"---- INITIALIZATION MODE --------------------------")
        CALL ppm_log(caller,cbuf,info)

        stdout_f('(A,A)',"   Initialization mode is:    ",'TRIM(init_mode)')
        CALL ppm_log(caller,cbuf,info)

        !-------------------------------------------------------------------------
        !  Initialize init kind
        !-------------------------------------------------------------------------
        stdout_f('(A)',"   Initialization parameters :")
        CALL ppm_log(caller,cbuf,info)

        CALL ppm_rc_uppercase(init_mode)
        SELECT CASE (TRIM(init_mode))
        CASE ("E_RECT","RECT")
           stdout_f('(A,3I5)',"                              ",'INT(init_rd(1:ppm_rc_dim))')
           CALL ppm_log(caller,cbuf,info)
           stdout_f('(A,3I5)',"                              ",'INT(init_sp(1:ppm_rc_dim))')
           CALL ppm_log(caller,cbuf,info)

        CASE ("E_SPHERE","SPHERE")
           stdout_f('(A,3I5)',"                              ",'INT(init_rd(1:ppm_rc_dim))')
           CALL ppm_log(caller,cbuf,info)
           stdout_f('(A,3I5)',"                              ",'INT(init_sp(1:ppm_rc_dim))')
           CALL ppm_log(caller,cbuf,info)

        CASE ("E_OTSU","OTSU")
           stdout_f('(A,4I5)',"                              ",'INT(m_lowerBound)','INT(m_upperBound)',histSize)
           CALL ppm_log(caller,cbuf,info)

        CASE ("E_LOCALMAX","LOCALMAX")
           stdout_f('(A,3I5)',"                              ",'INT(init_rd(1:ppm_rc_dim))')
           CALL ppm_log(caller,cbuf,info)

        CASE DEFAULT
           stdout_f('(A,3I5)',"                              ",'INT(init_rd(1:ppm_rc_dim))')
           CALL ppm_log(caller,cbuf,info)

        END SELECT

        !----------------------------------------------------------------------
        !  Energy related terms
        !----------------------------------------------------------------------
        stdout_f('(A)',"---- ENERGY TERMS ---------------------------------")
        CALL ppm_log(caller,cbuf,info)

        stdout_f('(A,A)',"   External energy term is:   ",'TRIM(energy_ext_name)')
        CALL ppm_log(caller,cbuf,info)

        stdout_f('(A,A)',"   Internal energy term is:   ",'TRIM(energy_int_name)')
        CALL ppm_log(caller,cbuf,info)

        !----------------------------------------------------------------------
        !  Digital topology (control)
        !----------------------------------------------------------------------
        stdout_f('(A)',"---- Digital topology (control) -------------------")
        CALL ppm_log(caller,cbuf,info)

        stdout_f('(A,L4)',"   Allow fusion:  ",AllowFusion)
        CALL ppm_log(caller,cbuf,info)

        stdout_f('(A,L4)',"   Allow fission: ",AllowFission)
        CALL ppm_log(caller,cbuf,info)

        stdout_f('(A,L4)',"   Allow handles: ",AllowHandles)
        CALL ppm_log(caller,cbuf,info)

        stdout_f('(A)',"---------------------------------------------------")
        CALL ppm_log(caller,cbuf,info)

        SELECT CASE (UseMCMC)
        CASE (.TRUE.)
           !----------------------------------------------------------------------
           !  MCMC
           !----------------------------------------------------------------------
           stdout_f('(A)',"---- MCMC algorithm related options ---------------")
           CALL ppm_log(caller,cbuf,info)

           stdout_f('(A,I8)',"   MC step size:                         ",MCMCstepsize)
           CALL ppm_log(caller,cbuf,info)

           stdout_f('(A,g12.5)',"   MC temperature:                       ",MCMCtemperature)
           CALL ppm_log(caller,cbuf,info)

           stdout_f('(A,g12.5)',"   MCMC burn in factor:                  ",MCMCburnInFactor)
           CALL ppm_log(caller,cbuf,info)

!            stdout_f('(A,g12.5)',"   MCMC off-boundary samples percentage: ",MCMCsampleOffBoundaryPercentage)
!            CALL ppm_log(caller,cbuf,info)

           stdout_f('(A,L4)',"   Continue (no relabeling at init):     ",MCMCcontinue)
           CALL ppm_log(caller,cbuf,info)

           stdout_f('(A,L4)',"   Use biased proposals:                 ",MCMCuseBiasedProposal)
           CALL ppm_log(caller,cbuf,info)

           stdout_f('(A,L4)',"   Use pair proposals:                   ",MCMCusePairProposal)
           CALL ppm_log(caller,cbuf,info)

           stdout_f('(A)',"---------------------------------------------------")
           CALL ppm_log(caller,cbuf,info)

        END SELECT

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
      END SUBROUTINE ppm_rc_dump_parameters
