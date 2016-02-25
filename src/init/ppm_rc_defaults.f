      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_rc_defaults
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
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_rc_defaults
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Define default/initial values for some of the
      !                 global variables.
      !
      !  Input        :
      !
      !  Input/output :
      !
      !  Output       : info     (I) return status
      !
      !  Routines     :
      !
      !  Remarks      :
      !
      !  References   :
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_rc_defaults(info)

        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------
        USE ppm_rc_module_mcmc, ONLY : MCMCsampleOffBoundaryPercentage,         &
        &   MCMCFloatingParticlesProposalNormalizer,MCMCZe,MCMCTotalNormalizer, &
        &   MCMCStatsUpdateModulus,MCMCVerbose,MCMCUseForbiddenRegion,          &
        &   MCMCuseSafeSampling,MCMCmarginalFileNamePrefix,MCMCRegionLabelSize, &
        &   MCMCAllParticlesFwdProposals,MCMC_CandidateMove,MCMC_PartnerMove,   &
        &   MCMC_q_A,MCMC_q_B,MCMC_q_A_B,MCMC_q_B_A,MCMC_qb_A,MCMC_qb_B,        &
        &   MCMC_qb_A_B,MCMC_qb_B_A,MCMC_LabelsBeforeJump_A,                    &
        &   MCMC_LabelsBeforeJump_B,MCMC_Particle_Ab_IsFloating,                &
        &   MCMC_Particle_Bb_IsFloating,MCMC_Particle_A_IsFloating,             &
        &   MCMC_CandidateMove_Index,MCMC_PartnerMove_Index
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER, INTENT(  OUT) :: info
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        REAL(ppm_kind_double) :: t0

        CHARACTER(LEN=ppm_char) :: caller = 'ppm_rc_defaults'
        !-------------------------------------------------------------------------
        ! Externals
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        ! Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        !-------------------------------------------------------------------------
        ! Log file units
        !-------------------------------------------------------------------------
        ppm_log_unit = 99

        !-------------------------------------------------------------------------
        ! Default pixel numbers if no input image is given
        !-------------------------------------------------------------------------
        Ngrid = 0
        !Ngrid(1:3)

        !-------------------------------------------------------------------------
        !  Time step
        !-------------------------------------------------------------------------
        istep = 0

        !-------------------------------------------------------------------------
        !  Default energy term
        !-------------------------------------------------------------------------
        SELECT CASE (ppm_rc_dim)
        CASE (2)
           ALLOCATE(init_rd(2),STAT=info)
           or_fail_alloc("init_rd")

           ALLOCATE(init_sp(2),STAT=info)
           or_fail_alloc("init_rd")

           ALLOCATE(ghostsize(2),STAT=info)
           or_fail_alloc("ghostsize")

           ALLOCATE(min_phys(2),max_phys(2),STAT=info)
           or_fail_alloc("min_phys & max_phys")

           ALLOCATE(bcdef(4),STAT=info)
           or_fail_alloc("bcdef")

           IF (AllowFusionZ) AllowFusionZ=.FALSE.

        CASE (3)
           ALLOCATE(init_rd(3), STAT=info)
           or_fail_alloc("blob_radius")

           ALLOCATE(init_sp(3),STAT=info)
           or_fail_alloc("init_rd")

           ALLOCATE(ghostsize(3),STAT=info)
           or_fail_alloc("ghostsize")

           ALLOCATE(min_phys(3),max_phys(3),STAT=info)
           or_fail_alloc("min_phys & max_phys")

           ALLOCATE(bcdef(6),STAT=info)
           or_fail_alloc("bcdef")

           IF (AllowFusionZ) AllowFusion=.TRUE.

        CASE DEFAULT
           fail("NOTICE: Wrong Case dimension.",ppm_error=ppm_error_fatal)

        END SELECT

        !-------------------------------------------------------------------------
        !  Default blob_radius
        !-------------------------------------------------------------------------
        IF (ANY(radius.GT.zero).AND.Sigma.GT.zero) THEN
           fail("You can not have both init_rd & Sigma at the same time!")
        ELSE IF (Sigma.GT.zero) then
           IF (ANY(rect.GT.zero)) THEN
              fail("You can not have two initial modes at the same time!")
           ENDIF
           init_rd=Sigma
        ELSE IF (ANY(radius.GT.zero)) THEN
           IF (ANY(rect.GT.zero)) THEN
              fail("You can not have two initial modes at the same time!")
           ENDIF
           init_rd=radius
        ELSE IF (ANY(rect.GT.zero)) THEN
           IF (ANY(radius.GT.zero).OR.Sigma.GT.zero) THEN
              fail("You can not have two initial modes at the same time!")
           ENDIF
           init_rd=rect
        ELSE
           init_rd=-one
        ENDIF

        IF (ANY(space.GT.zero)) THEN
           init_sp=space
        ELSE
           init_sp=-one
        ENDIF

        !-------------------------------------------------------------------------
        !  Default ghostsize
        !-------------------------------------------------------------------------
        ghostsize = 0

        !-------------------------------------------------------------------------
        !  Default domain size
        !-------------------------------------------------------------------------
        min_phys = zero
        max_phys = one

        !-------------------------------------------------------------------------
        !  Default domain boundary conditions
        !-------------------------------------------------------------------------
        bcdef = ppm_param_bcdef_freespace

        IF (lNormalize) Normalfac=one

        DEALLOCATE(radius,rect,space,STAT=info)
        or_fail_dealloc("radius,rect,space")

        !-------------------------------------------------------------------------
        !  Default debugging time variables
        !-------------------------------------------------------------------------
        IF (debug.GT.0) THEN
           tmove_Simple    = zerod
           tmove_SimpleComm= zerod
           tmove_NotSimple = zerod
           tmove_Part      = zerod
           tmove_ghostfire = zerod
        ENDIF

        !-------------------------------------------------------------------------
        !  Default topo & mesh ID
        !-------------------------------------------------------------------------
        iotopoid=0
        iomeshid=-1
        initopoid=0
        initmeshid=-1
        topoid=0
        meshid=-1

        !-------------------------------------------------------------------------
        !  Default for MCMC
        !-------------------------------------------------------------------------
        SELECT CASE (UseMCMC)
        CASE (.TRUE.)
           IF (MCMCstepsize.GT.1.AND.MCMCusePairProposal) THEN
              fail("Using pair proposal with step size > 1 leads to slightly wrong proposals and hence detailed balance is not guaranteed!", &
              & ppm_error=ppm_error_fatal)
           ENDIF

           MCMCsampleOffBoundaryPercentage         = 0.05_MK
           MCMCFloatingParticlesProposalNormalizer = zerod
           MCMCZe                                  = zero
           MCMCTotalNormalizer                     = zerod

           MCMCStatsUpdateModulus                  = 10000
           MCMCRegionLabelSize                     = 0

           MCMCVerbose                             = .FALSE.
           MCMCUseForbiddenRegion                  = .FALSE.
           MCMCuseSafeSampling                     = .FALSE.

           MCMCmarginalFileNamePrefix              = "marginals_region_"

           ALLOCATE(MCMCAllParticlesFwdProposals(0),STAT=info)
           or_fail_alloc("MCMCAllParticlesFwdProposals")

           ALLOCATE(MCMC_CandidateMove(MCMCstepsize),                  &
           &        MCMC_CandidateMove_Index(ppm_rc_dim,MCMCstepsize), &
           &        MCMC_PartnerMove_Index(ppm_rc_dim,MCMCstepsize),   &
           &        MCMC_PartnerMove(MCMCstepsize),STAT=info)
           or_fail_alloc("MCMCCandidateMove,MCMC_CandidateMove_Index,ppm_rc_dim,MCMCstepsize & MCMCPartnerMove")

           ALLOCATE(MCMC_q_A(MCMCstepsize),   MCMC_q_B(MCMCstepsize),    &
           &        MCMC_q_A_B(MCMCstepsize), MCMC_q_B_A(MCMCstepsize),  &
           &        MCMC_qb_A(MCMCstepsize),  MCMC_qb_B(MCMCstepsize),   &
           &        MCMC_qb_A_B(MCMCstepsize),MCMC_qb_B_A(MCMCstepsize), &
           &        MCMC_LabelsBeforeJump_A(MCMCstepsize),               &
           &        MCMC_LabelsBeforeJump_B(MCMCstepsize),STAT=info)
           or_fail_alloc("MCMC_q_A,q_B,q_A_B,q_B_A,qb_A,qb_B,qb_A_B,qb_B_A,LabelsBeforeJump_A & LabelsBeforeJump_B")

           ALLOCATE(MCMC_Particle_Ab_IsFloating(MCMCstepsize), &
           &        MCMC_Particle_Bb_IsFloating(MCMCstepsize), &
           &        MCMC_Particle_A_IsFloating(MCMCstepsize),STAT=info)
           or_fail_alloc("Particle_Ab_IsFloating,Particle_Bb_IsFloating & Particle_A_IsFloating")

        END SELECT

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
      END SUBROUTINE ppm_rc_defaults
