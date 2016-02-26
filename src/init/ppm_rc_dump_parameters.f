      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_rc_dump_parameters
      !-------------------------------------------------------------------------
      ! Copyright (c) 2016 MOSAIC Group (MPI-CBG Dresden)
      !
      !
      ! This file is part of the PPM_RC program.
      !
      ! PPM_RC is free software: you can redistribute it and/or modify
      ! it under the terms of the GNU Lesser General Public License
      ! as published by the Free Software Foundation, either
      ! version 3 of the License, or (at your option) any later
      ! version.
      !
      ! PPM_RC is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! and the GNU Lesser General Public License along with PPM_RC. If not,
      ! see <http://www.gnu.org/licenses/>.
      !-------------------------------------------------------------------------
      !  MOSAIC Group
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

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

        stdout_f('(A,3I8)',"   Number of pixels:          ",'Ngrid(1:ppm_rc_dim)')
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

        stdout_f('(A,I8)',"   Time steps between output:",freqoutput)
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
           stdout_f('(A,3I8)',"                              ",'INT(init_rd(1:ppm_rc_dim))')
           CALL ppm_log(caller,cbuf,info)
           stdout_f('(A,3I8)',"                              ",'INT(init_sp(1:ppm_rc_dim))')
           CALL ppm_log(caller,cbuf,info)

        CASE ("E_SPHERE","SPHERE")
           stdout_f('(A,3I8)',"                              ",'INT(init_rd(1:ppm_rc_dim))')
           CALL ppm_log(caller,cbuf,info)
           stdout_f('(A,3I8)',"                              ",'INT(init_sp(1:ppm_rc_dim))')
           CALL ppm_log(caller,cbuf,info)

        CASE ("E_OTSU","OTSU")
           stdout_f('(A,3I8)',"                              ",'INT(m_lowerBound)','INT(m_upperBound)',histSize)
           CALL ppm_log(caller,cbuf,info)

        CASE ("E_LOCALMAX","LOCALMAX")
           stdout_f('(A,3I8)',"                              ",'INT(init_rd(1:ppm_rc_dim))')
           CALL ppm_log(caller,cbuf,info)

        CASE DEFAULT
           stdout_f('(A,3I8)',"                              ",'INT(init_rd(1:ppm_rc_dim))')
           CALL ppm_log(caller,cbuf,info)

        END SELECT

        !----------------------------------------------------------------------
        !  Energy related terms
        !----------------------------------------------------------------------
        stdout_f('(A)',"---- ENERGY TERMS ---------------------------------")
        CALL ppm_log(caller,cbuf,info)

        stdout_f('(A,A)',"   External energy term is:   ",'TRIM(energy_ext_name)')
        CALL ppm_log(caller,cbuf,info)

        stdout_f('(A)',"")
        CALL ppm_log(caller,cbuf,info)

        stdout_f('(A,f18.8)',"   Data term coefficient:                ",energy_coeff_data)
        CALL ppm_log(caller,cbuf,info)

        stdout_f('(A)',"")
        CALL ppm_log(caller,cbuf,info)

        stdout_f('(A,f18.8)',"   Threshold prior for region merging:   ",energy_region_merge_ths)
        CALL ppm_log(caller,cbuf,info)

        SELECT CASE (TRIM(energy_ext_name))
        CASE ("PS","PSGAUSSIAN","PSPOISSON")
           stdout_f('(A)',"")
           CALL ppm_log(caller,cbuf,info)

           stdout_f('(A,f18.8)',"   Radius of the spherical patches:      ",energy_local_window_radius)
           CALL ppm_log(caller,cbuf,info)

           stdout_f('(A)',"   (Radius of mask for local energy support in pixel)")
           CALL ppm_log(caller,cbuf,info)

           stdout_f('(A)',"")
           CALL ppm_log(caller,cbuf,info)

           stdout_f('(A,f18.8)',"   Outward balloon flow coefficient:     ",energy_coeff_balloon)
           CALL ppm_log(caller,cbuf,info)
        END SELECT

        stdout_f('(A)',"")
        CALL ppm_log(caller,cbuf,info)

        stdout_f('(A,A)',"   Internal energy term is:   ",'TRIM(energy_int_name)')
        CALL ppm_log(caller,cbuf,info)

        stdout_f('(A,f18.8)',"   Length term coefficient:              ",energy_coeff_length)
        CALL ppm_log(caller,cbuf,info)

        SELECT CASE (TRIM(energy_int_name))
        CASE ("CURV")
           stdout_f('(A)',"")
           CALL ppm_log(caller,cbuf,info)

           stdout_f('(A)',"   Curvature mask radius for regularizing")
           CALL ppm_log(caller,cbuf,info)

           stdout_f('(A,f18.8)',"   flows. (See Kybic and Kratky) :       ",energy_curvature_mask_radius)
           CALL ppm_log(caller,cbuf,info)
        END SELECT

        IF (energy_coeff_outward_flow.GT.smallest.OR. &
        &   energy_coeff_outward_flow.LT.-smallest) THEN
           stdout_f('(A)',"")
           CALL ppm_log(caller,cbuf,info)

           stdout_f('(A,f18.8)',"   Coefficient of constant outward flow: ",energy_coeff_outward_flow)
           CALL ppm_log(caller,cbuf,info)
        ENDIF

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

           stdout_f('(A,I8)',"   MCMC step size:                       ",MCMCstepsize)
           CALL ppm_log(caller,cbuf,info)

           stdout_f('(A,f18.8)',"   MCMC temperature:                     ",MCMCtemperature)
           CALL ppm_log(caller,cbuf,info)

           stdout_f('(A,f18.8)',"   MCMC burn in factor:                  ",MCMCburnInFactor)
           CALL ppm_log(caller,cbuf,info)

!            stdout_f('(A,f18.8)',"   MCMC off-boundary samples percentage: ",MCMCsampleOffBoundaryPercentage)
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
