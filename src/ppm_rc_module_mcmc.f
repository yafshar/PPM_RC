      !-------------------------------------------------------------------------
      !  Module       :                    ppm_rc_module_mcmc
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
      !  Max Planck Institute of Molecular Cell Biology and Genetics
      !  Pfotenhauerstr. 108, 01307 Dresden, Germany
      !
      !  Author           - y.afshar           Nov    2015
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Module       :                    ppm_rc_module_mcmc
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Markov Chain Monte Carlo based sampling of regions in
      !                 an image domain given a fixed number of regions.
      !
      !  Remarks      :
      !                 The probability to minimize is induced by the registered
      !                 energy functionals.
      !                 Any internal or external energy functional might be registered.
      !                  The filter outputs a probability map ...
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      MODULE ppm_rc_module_mcmc
        !----------------------------------------------------------------------
        !  Modules
        !----------------------------------------------------------------------
        USE ppm_module_data, ONLY : ppm_char,ppm_kind_double
        USE ppm_module_inl_hash, ONLY : ppm_htable

        USE ppm_rc_module_global, ONLY : MK,zero,zerod
        USE ppm_rc_module_util, ONLY : ppm_rc_htable,ppm_rc_MCMCParticlehtable, &
        &   ppm_rc_MCMCHistoryParticlehtable,MCMCParticle,MCMCHistoryParticle
        USE ppm_rc_module_energy, ONLY : e_data,e_length
        !NOTE:
        !Intel compiler Error, when I use ppm_rc_module_energy in
        !the corresponding subroutines the compiler complains about internal error
        !which it should not as there is no circular dependency?
        !So as a current solution I put it here
        IMPLICIT NONE

        PRIVATE
        !----------------------------------------------------------------------
        !  MCMC TYPEs
        !----------------------------------------------------------------------
        TYPE(ppm_htable)                                           :: MCMChtable
        !!! hash table for the region labels

        TYPE(ppm_rc_MCMCParticlehtable), DIMENSION(:), ALLOCATABLE, TARGET :: MCMCparents
        TYPE(ppm_rc_MCMCParticlehtable), DIMENSION(:), ALLOCATABLE, TARGET :: MCMCchildren
        TYPE(ppm_rc_MCMCParticlehtable),                            TARGET :: MCMCFloatingParticles
        TYPE(ppm_rc_MCMCParticlehtable)                                    :: MCMCAppliedParticles

        REAL(MK)                                                   :: MCMCsampleOffBoundaryPercentage
        REAL(ppm_kind_double)                                      :: MCMCFloatingParticlesProposalNormalizer
        REAL(ppm_kind_double)                                      :: MCMCProbabilityToProposeAFloatingParticle

        REAL(MK)                                                   :: MCMCZe
        !!! Normalizing constants of the edge image
        REAL(MK),                        DIMENSION(:), ALLOCATABLE :: MCMClengthProposalMask
        !!! For a fast proposal computation
        REAL(MK),                        DIMENSION(:), ALLOCATABLE :: MCMCparentsProposalNormalizer
        REAL(MK),                        DIMENSION(:), ALLOCATABLE :: MCMCchildrenProposalNormalizer
        !!! map containg a label and its float Normalize factor
        REAL(ppm_kind_double)                                      :: MCMCTotalNormalizer

        REAL(MK)                                                   :: MCMCTotEnergyDiff



        TYPE(ppm_rc_MCMCParticlehtable),               POINTER     :: MCMCActiveCandidates
        !!! Temporary pointer to an object

        REAL(MK),                        DIMENSION(:), ALLOCATABLE :: MCMCAllParticlesFwdProposals
        !!! For each particle within the region, calculate the proposal

        INTEGER                                                    :: MCMCActiveCandidatesSize
        !!! Size of the MCMCAllParticlesFwdProposals


        INTEGER                                                    :: MCMCStatsUpdateModulus
        INTEGER,                         DIMENSION(:), ALLOCATABLE :: MCMCVisitedLabels
        !!! Sordered list of visited labels, it is good when you wanna loop through
        !!! Sweep through the label image. If a label didn't exist, register it
        !!! (add an entry in the regionsLabel vector and initialize normalizer
        !!! statistics). And always initialize these for the BG label (=0).
        INTEGER,                         DIMENSION(:), ALLOCATABLE :: MCMCRegionLabel
        !!! vector of labels
        INTEGER                                                    :: MCMCRegionLabelSize




        LOGICAL                                                    :: MCMCVerbose
        LOGICAL                                                    :: MCMCUseForbiddenRegion
        LOGICAL                                                    :: MCMCuseSafeSampling
        !!!
        LOGICAL                                                    :: MCMCHardReject

        CHARACTER(LEN=ppm_char)                                    :: MCMCmarginalFileNamePrefix


        PUBLIC :: MCMChtable

        PUBLIC :: MCMCparents
        PUBLIC :: MCMCchildren
        PUBLIC :: MCMCFloatingParticles
        PUBLIC :: MCMCAppliedParticles

        PUBLIC :: MCMCsampleOffBoundaryPercentage
        PUBLIC :: MCMCFloatingParticlesProposalNormalizer
        PUBLIC :: MCMCProbabilityToProposeAFloatingParticle
        PUBLIC :: MCMCZe

        PUBLIC :: MCMClengthProposalMask
        PUBLIC :: MCMCparentsProposalNormalizer
        PUBLIC :: MCMCchildrenProposalNormalizer

        PUBLIC :: MCMCTotalNormalizer
        PUBLIC :: MCMCTotEnergyDiff

        PUBLIC :: MCMCActiveCandidates
        PUBLIC :: MCMCAllParticlesFwdProposals
        PUBLIC :: MCMCActiveCandidatesSize

        PUBLIC :: MCMCStatsUpdateModulus
        PUBLIC :: MCMCVisitedLabels

        PUBLIC :: MCMCRegionLabel
        PUBLIC :: MCMCRegionLabelSize

        PUBLIC :: MCMCVerbose
        PUBLIC :: MCMCUseForbiddenRegion
        PUBLIC :: MCMCuseSafeSampling
        PUBLIC :: MCMCHardReject

        PUBLIC :: MCMCmarginalFileNamePrefix



        !----------------------------------------------------------------------
        !  Some work memory on the heap
        !----------------------------------------------------------------------
!         REAL(MK), DIMENSION(:), ALLOCATABLE :: tmp1_r

        INTEGER,  DIMENSION(:), ALLOCATABLE :: tmp1_i

        TYPE(MCMCParticle), DIMENSION(:), ALLOCATABLE :: MCMC_CandidateMove
        TYPE(MCMCParticle), DIMENSION(:), ALLOCATABLE :: MCMC_PartnerMove

        REAL(MK), DIMENSION(:), ALLOCATABLE :: MCMC_q_A
        REAL(MK), DIMENSION(:), ALLOCATABLE :: MCMC_q_B
        REAL(MK), DIMENSION(:), ALLOCATABLE :: MCMC_q_A_B
        REAL(MK), DIMENSION(:), ALLOCATABLE :: MCMC_q_B_A
        !!! Forward probability
        REAL(MK), DIMENSION(:), ALLOCATABLE :: MCMC_qb_A
        REAL(MK), DIMENSION(:), ALLOCATABLE :: MCMC_qb_B
        REAL(MK), DIMENSION(:), ALLOCATABLE :: MCMC_qb_A_B
        REAL(MK), DIMENSION(:), ALLOCATABLE :: MCMC_qb_B_A
        !!! Backward probability

        INTEGER, DIMENSION(:),   ALLOCATABLE :: MCMC_LabelsBeforeJump_A
        INTEGER, DIMENSION(:),   ALLOCATABLE :: MCMC_LabelsBeforeJump_B
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: MCMC_CandidateMove_Index
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: MCMC_PartnerMove_Index

        LOGICAL, DIMENSION(:), ALLOCATABLE :: MCMC_Particle_A_IsFloating
        LOGICAL, DIMENSION(:), ALLOCATABLE :: MCMC_Particle_Ab_IsFloating
        LOGICAL, DIMENSION(:), ALLOCATABLE :: MCMC_Particle_Bb_IsFloating

        LOGICAL                            :: MCMC_SingleParticleMoveForPairProposals

        PUBLIC :: MCMC_CandidateMove
        PUBLIC :: MCMC_CandidateMove_Index
        PUBLIC :: MCMC_PartnerMove_Index
        PUBLIC :: MCMC_PartnerMove
        PUBLIC :: MCMC_q_A
        PUBLIC :: MCMC_q_B
        PUBLIC :: MCMC_q_A_B
        PUBLIC :: MCMC_q_B_A
        PUBLIC :: MCMC_qb_A
        PUBLIC :: MCMC_qb_B
        PUBLIC :: MCMC_qb_A_B
        PUBLIC :: MCMC_qb_B_A
        PUBLIC :: MCMC_LabelsBeforeJump_A
        PUBLIC :: MCMC_LabelsBeforeJump_B
        PUBLIC :: MCMC_Particle_Ab_IsFloating
        PUBLIC :: MCMC_Particle_Bb_IsFloating
        PUBLIC :: MCMC_Particle_A_IsFloating
        PUBLIC :: MCMC_SingleParticleMoveForPairProposals

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE MCMCInsertLabelInRegionLabel
          MODULE PROCEDURE MCMCInsertBGLabelInRegionLabel
          MODULE PROCEDURE MCMCInsertLabelInRegionLabel
        END INTERFACE

        INTERFACE MCMCproposal
          MODULE PROCEDURE MCMCproposal_2d
          MODULE PROCEDURE MCMCproposal__2d
          MODULE PROCEDURE MCMCproposal_3d
          MODULE PROCEDURE MCMCproposal__3d
        END INTERFACE

        INTERFACE MCMCIsParticleTopoValid
          MODULE PROCEDURE MCMCIsParticleTopoValid_2d
          MODULE PROCEDURE MCMCIsParticleTopoValid__2d
          MODULE PROCEDURE MCMCIsParticleTopoValid_3d
          MODULE PROCEDURE MCMCIsParticleTopoValid__3d
        END INTERFACE

        INTERFACE MCMCGetIndexFromEdgeDensity
          MODULE PROCEDURE MCMCGetIndexFromEdgeDensity_2d
          MODULE PROCEDURE MCMCGetIndexFromEdgeDensity_3d
        END INTERFACE

        INTERFACE MCMCgetRegularParticlesAtIndex
          MODULE PROCEDURE MCMCgetRegularParticlesAtIndex_2d
          MODULE PROCEDURE MCMCgetRegularParticlesAtIndex__2d
          MODULE PROCEDURE MCMCgetRegularParticlesAtIndex_3d
          MODULE PROCEDURE MCMCgetRegularParticlesAtIndex__3d
        END INTERFACE

        INTERFACE MCMCgetParticlesInBGNeighborhood
          MODULE PROCEDURE MCMCgetParticlesInBGNeighborhood_2d
          MODULE PROCEDURE MCMCgetParticlesInBGNeighborhood_3d
        END INTERFACE

        INTERFACE MCMCgetParticlesInFGNeighborhood
          MODULE PROCEDURE MCMCgetParticlesInFGNeighborhood_2d
          MODULE PROCEDURE MCMCgetParticlesInFGNeighborhood_3d
        END INTERFACE

        INTERFACE MCMCInsertCandidatesToContainers
          MODULE PROCEDURE MCMCInsertCandidatesToContainers_2d
          MODULE PROCEDURE MCMCInsertCandidatesToContainers_3d
        END INTERFACE

        INTERFACE MCMCEraseCandidatesFromContainers
          MODULE PROCEDURE MCMCEraseCandidatesFromContainers_2d
          MODULE PROCEDURE MCMCEraseCandidatesFromContainers_3d
        END INTERFACE

        INTERFACE MCMCEraseFloatingParticle
          MODULE PROCEDURE MCMCEraseFloatingParticle_2d
          !!! first find out whether there is a particle at coords or not
          MODULE PROCEDURE MCMCEraseFloatingParticle__2d
          !!! give the searched tmpParticle at coords as an input
          MODULE PROCEDURE MCMCEraseFloatingParticle_3d
          MODULE PROCEDURE MCMCEraseFloatingParticle__3d
        END INTERFACE

        INTERFACE MCMCInsertFloatingParticle
          MODULE PROCEDURE MCMCInsertFloatingParticle_2d
          MODULE PROCEDURE MCMCInsertFloatingParticle_3d
        END INTERFACE

        INTERFACE MCMCIsRegularParticle
          MODULE PROCEDURE MCMCIsRegularParticle_2d
          MODULE PROCEDURE MCMCIsRegularParticle_3d
        END INTERFACE

        INTERFACE MCMCAddAndRemoveParticlesWhenMove
          MODULE PROCEDURE MCMCAddAndRemoveParticlesWhenMove_2d
          MODULE PROCEDURE MCMCAddAndRemoveParticlesWhenMove_3d
        END INTERFACE

        INTERFACE MCMCupdateProposalsAndFilterTopologyInNeighborhood
          MODULE PROCEDURE MCMCupdateProposalsAndFilterTopologyInNeighborhood_2d
          MODULE PROCEDURE MCMCupdateProposalsAndFilterTopologyInNeighborhood_3d
        END INTERFACE

        INTERFACE MCMCApplyParticle
          MODULE PROCEDURE MCMCApplyParticle_2d
          MODULE PROCEDURE MCMCApplyParticle_3d
        END INTERFACE

        !----------------------------------------------------------------------
        ! Public
        !----------------------------------------------------------------------
        PUBLIC :: MCMCCheckParameters
        PUBLIC :: CreateMCMClengthProposalMask
        PUBLIC :: DestroyMCMClengthProposalMask
        PUBLIC :: MCMCInsertLabelInRegionLabel
        PUBLIC :: MCMCUpdateRegionLabel
        PUBLIC :: MCMCproposal
        PUBLIC :: MCMCIsParticleTopoValid
        PUBLIC :: MCMCGetIndexFromEdgeDensity
        PUBLIC :: MCMCgetRegularParticlesAtIndex
        PUBLIC :: MCMCgetParticlesInBGNeighborhood
        PUBLIC :: MCMCgetParticlesInFGNeighborhood
        PUBLIC :: MCMCInsertCandidatesToContainers
        PUBLIC :: MCMCEraseCandidatesFromContainers
        PUBLIC :: MCMCEraseFloatingParticle
        PUBLIC :: MCMCInsertFloatingParticle
        PUBLIC :: MCMCOffBoundarySample_2d
        PUBLIC :: MCMCOffBoundarySample_3d
        PUBLIC :: MCMCUpdateAllParticlesFwdProposals
        PUBLIC :: MCMCIsRegularParticle
        PUBLIC :: MCMCFindParticleA
        PUBLIC :: MCMCAddAndRemoveParticlesWhenMove
        PUBLIC :: MCMCupdateProposalsAndFilterTopologyInNeighborhood
        PUBLIC :: MCMCApplyParticle
        PUBLIC :: MCMCFindParticleB
!         PUBLIC :: MCMCMove

      CONTAINS

        SUBROUTINE MCMCCheckParameters(info)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_argument

          USE ppm_rc_module_global, ONLY : MCMCstepsize,MCMCusePairProposal
          IMPLICIT NONE

          INTEGER, INTENT(INOUT) :: info
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: t0

          CHARACTER(LEN=ppm_char) :: caller='MCMCCheckParameters'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          IF (MCMCstepsize.GT.1.AND.MCMCusePairProposal) THEN
             fail("Using pair potential with step size > 1 leads to slightly wrong proposals and hence detailed balance is not guaranteed.", &
             & ppm_error=ppm_error_fatal)
          ENDIF
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END SUBROUTINE MCMCCheckParameters

        SUBROUTINE CreateMCMClengthProposalMask(info)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_alloc,ppm_err_dealloc, &
          &   ppm_err_argument

          USE ppm_rc_module_global, ONLy : one,ppm_rc_dim,MCMCuseBiasedProposal
          USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType,BG_ConnectivityType
          IMPLICIT NONE

          INTEGER, INTENT(INOUT) :: info

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: t0
          REAL(MK)              :: vLength

          INTEGER :: i,j,NumberOfNeighbors

          CHARACTER(LEN=ppm_char) :: caller='CreateMCMClengthProposalMask'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          IF (.NOT.MCMCuseBiasedProposal) GOTO 9999

          CALL check

          IF (ALLOCATED(MCMClengthProposalMask)) THEN
             DEALLOCATE(MCMClengthProposalMask,STAT=info)
             or_fail_dealloc("MCMClengthProposalMask")
          ENDIF

          NumberOfNeighbors=BG_ConnectivityType%NumberOfNeighbors

          ALLOCATE(MCMClengthProposalMask(NumberOfNeighbors),STAT=info)
          or_fail_alloc("Failed to allocate MCMClengthProposalMask")

          DO i=1,NumberOfNeighbors
             vLength=REAL(BG_ConnectivityType%NeighborsPoints(1,i),MK)* &
             &       REAL(BG_ConnectivityType%NeighborsPoints(1,i),MK)
             DO j=2,ppm_rc_dim
                vLength=vLength+REAL(BG_ConnectivityType%NeighborsPoints(j,i),MK)* &
                &               REAL(BG_ConnectivityType%NeighborsPoints(j,i),MK)
             ENDDO
             MCMClengthProposalMask(i)=one/vLength
          ENDDO
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        CONTAINS
        SUBROUTINE check
        IMPLICIT NONE
          check_true(<#ASSOCIATED(FG_ConnectivityType)#>, &
          & "Topology Connectivity is not defined yet!",exit_point=8888)
        8888 CONTINUE
        END SUBROUTINE check
        END SUBROUTINE CreateMCMClengthProposalMask

        SUBROUTINE DestroyMCMClengthProposalMask(info)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_dealloc
          IMPLICIT NONE

          INTEGER, INTENT(INOUT) :: info

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: t0

          CHARACTER(LEN=ppm_char) :: caller='DestroyMCMClengthProposalMask'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)
          IF (ALLOCATED(MCMClengthProposalMask)) THEN
             DEALLOCATE(MCMClengthProposalMask,STAT=info)
             or_fail_dealloc("MCMClengthProposalMask")
          ENDIF
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END SUBROUTINE DestroyMCMClengthProposalMask

        FUNCTION MCMCInsertBGLabelInRegionLabel() RESULT(info)
          !!!
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_alloc,ppm_err_argument
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          REAL(ppm_kind_double) :: t0

          CHARACTER(LEN=ppm_char) :: caller='MCMCInsertLabelInRegionLabel'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          IF (MCMCRegionLabelSize.EQ.0) THEN
             MCMCRegionLabelSize=1
             ALLOCATE(MCMCRegionLabel(MCMCRegionLabelSize),STAT=info)
             or_fail_alloc("MCMCRegionLabel",ppm_error=ppm_error_fatal)
             MCMCRegionLabel(1)=0
          ELSE
             fail("MCMCRegionLabel is not allocated yet and can not have MCMCRegionLabelSize bigger than 0!", &
             & ppm_error=ppm_error_fatal)
          ENDIF
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCInsertBGLabelInRegionLabel

        FUNCTION MCMCInsertLabelInRegionLabel(LabelIn) RESULT(info)
          !!!
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_alloc,ppm_err_sub_failed
          USE ppm_module_util_qsort, ONLY : ppm_util_qsort

          USE ppm_rc_module_util, ONLY : label_exist
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, INTENT(IN   ) :: LabelIn
          INTEGER                :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          REAL(ppm_kind_double) :: t0

          INTEGER :: i

          CHARACTER(LEN=ppm_char) :: caller='MCMCInsertLabelInRegionLabel'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          IF (.NOT.label_exist(LabelIn,MCMCRegionLabel,MCMCRegionLabelSize)) THEN
             IF (MCMCRegionLabelSize+1.GT.SIZE(MCMCRegionLabel)) THEN
                ALLOCATE(tmp1_i(MCMCRegionLabelSize*2),STAT=info)
                or_fail_alloc("tmp1_i",ppm_error=ppm_error_fatal)

                FORALL (i=1:MCMCRegionLabelSize) tmp1_i(i)=MCMCRegionLabel(i)

                CALL MOVE_ALLOC(tmp1_i,MCMCRegionLabel)
             ENDIF

             MCMCRegionLabelSize=MCMCRegionLabelSize+1
             MCMCRegionLabel(MCMCRegionLabelSize)=LabelIn

             CALL ppm_util_qsort(MCMCRegionLabel,info,MCMCRegionLabelSize)
             or_fail("ppm_util_qsort",ppm_error=ppm_error_fatal)
          ENDIF
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCInsertLabelInRegionLabel

        FUNCTION MCMCUpdateRegionLabel() RESULT(info)
          !!! Set up a vector MCMCRegionLabel that maps natural numbers, index of
          !!! the vector to the region labels (which are local on this processor).
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_dealloc,ppm_err_sub_failed
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          REAL(ppm_kind_double) :: t0

          INTEGER :: nsize
          INTEGER :: i,RLabel

          CHARACTER(LEN=ppm_char) :: caller='MCMCUpdateRegionLabel'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          DEALLOCATE(MCMCRegionLabel,STAT=info)
          or_fail_dealloc("MCMCRegionLabel",ppm_error=ppm_error_fatal)

          MCMCRegionLabelSize=0

          info=MCMCInsertLabelInRegionLabel()
          or_fail("MCMCInsertLabelInRegionLabel",ppm_error=ppm_error_fatal)

          nsize = SIZE(MCMCparents)-1
          DO i=1,nsize
             IF (MCMCparents(i)%size().GT.0) THEN
                RLabel=e_data%Rlabel(i)

                info=MCMCInsertLabelInRegionLabel(RLabel)
                or_fail("MCMCInsertLabelInRegionLabel",ppm_error=ppm_error_fatal)
             ENDIF
          ENDDO

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCUpdateRegionLabel

        FUNCTION MCMCproposal_2d(coord,labels_)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLy : zero,half,one,MCMCuseBiasedProposal
          USE ppm_rc_module_topologicalnumber, ONLY : BG_ConnectivityType
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: coord
          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_

          REAL(MK)                                           :: MCMCproposal_2d
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER, DIMENSION(:,:), POINTER :: tmplabels
          INTEGER, DIMENSION(2)            :: ll
          INTEGER                          :: cLabel
          INTEGER                          :: LabelNgh
          INTEGER                          :: i

          LOGICAL :: vFloating
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          IF (.NOT.MCMCuseBiasedProposal) THEN
             MCMCproposal_2d=one
             RETURN
          ENDIF
          ! The choice of a data-driven proposal distribution depends a lot on
          ! the contrast in the image (?) and the current segmentation.

          ! L1: 1 / ( 1 + abs(d) )
          ! MCMCproposal_2d=1./(1.+abs(Means[Particle.m_CandidateLabel]-Image(Particle.Index)))
          ! L2: 1 / d^2
          ! MCMCproposal_2d=1./(1.+abs(Means[Particle.m_CandidateLabel]-Image(Particle.Index)))
          ! exp(-d^2)
          ! MCMCproposal_2d=exp((Means[Particle.CandidateLabel] - Image(Particle.Index)) * &
          ! &                   (Means[Particle.CandidateLabel] - Image(Particle.Index)))

          ! Length-prior driven proposal:
          MCMCproposal_2d=zero
          vFloating=.TRUE.

          tmplabels => labels_(coord(1)-1:coord(1)+1,coord(2)-1:coord(2)+1)
          cLabel=ABS(tmplabels(2,2))

          DO i=1,BG_ConnectivityType%NumberOfNeighbors
             ll=BG_ConnectivityType%NeighborsPoints(:,i)+2
             LabelNgh=ABS(tmplabels(ll(1),ll(2)))
             IF (cLabel.NE.LabelNgh) THEN
                MCMCproposal_2d=MCMCproposal_2d+MCMClengthProposalMask(i)
                vFloating=.FALSE.
             ENDIF
          ENDDO

          IF (vFloating) THEN
             ! floating particles need a proposal > 0: We take half of the smallest
             ! element in the mask (hack:we assume the smallest element to be at
             ! position 0).
             MCMCproposal_2d = MCMClengthProposalMask(1)*half
          ENDIF
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCproposal_2d

        FUNCTION MCMCproposal__2d(tmplabels)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLy : zero,half,one,MCMCuseBiasedProposal
          USE ppm_rc_module_topologicalnumber, ONLY : BG_ConnectivityType
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, DIMENSION(:,:), POINTER :: tmplabels
          !!! Temporary label index around coordinate

          REAL(MK)                         :: MCMCproposal__2d
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER, DIMENSION(2) :: ll
          INTEGER               :: cLabel
          INTEGER               :: LabelNgh
          INTEGER               :: i

          LOGICAL :: vFloating
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          IF (.NOT.MCMCuseBiasedProposal) THEN
             MCMCproposal__2d=one
             RETURN
          ENDIF
          ! The choice of a data-driven proposal distribution depends a lot on
          ! the contrast in the image (?) and the current segmentation.

          ! L1: 1 / ( 1 + abs(d) )
          ! MCMCproposal__2d=1./(1.+abs(Means[Particle.m_CandidateLabel]-Image(Particle.Index)))
          ! L2: 1 / d^2
          ! MCMCproposal__2d=1./(1.+abs(Means[Particle.m_CandidateLabel]-Image(Particle.Index)))
          ! exp(-d^2)
          ! MCMCproposal__2d=exp((Means[Particle.CandidateLabel] - Image(Particle.Index)) * &
          ! &                   (Means[Particle.CandidateLabel] - Image(Particle.Index)))

          ! Length-prior driven proposal:
          MCMCproposal__2d=zero
          vFloating=.TRUE.

          cLabel=ABS(tmplabels(2,2))

          DO i=1,BG_ConnectivityType%NumberOfNeighbors
             ll=BG_ConnectivityType%NeighborsPoints(:,i)+2
             LabelNgh=ABS(tmplabels(ll(1),ll(2)))
             IF (cLabel.NE.LabelNgh) THEN
                MCMCproposal__2d=MCMCproposal__2d+MCMClengthProposalMask(i)
                vFloating=.FALSE.
             ENDIF
          ENDDO

          IF (vFloating) THEN
             ! floating particles need a proposal > 0: We take half of the smallest
             ! element in the mask (hack:we assume the smallest element to be at
             ! position 0).
             MCMCproposal__2d = MCMClengthProposalMask(1)*half
          ENDIF
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCproposal__2d

        FUNCTION MCMCproposal_3d(coord,labels_)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLy : zero,half,one,MCMCuseBiasedProposal
          USE ppm_rc_module_topologicalnumber, ONLY : BG_ConnectivityType
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: coord
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_

          REAL(MK)                                             :: MCMCproposal_3d
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER, DIMENSION(:,:,:), POINTER :: tmplabels
          INTEGER, DIMENSION(3)              :: ll
          INTEGER                            :: cLabel
          INTEGER                            :: LabelNgh
          INTEGER                            :: i

          LOGICAL :: vFloating
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          IF (.NOT.MCMCuseBiasedProposal) THEN
             MCMCproposal_3d=one
             RETURN
          ENDIF
          ! The choice of a data-driven proposal distribution depends a lot on
          ! the contrast in the image (?) and the current segmentation.

          ! L1: 1 / ( 1 + abs(d) )
          ! MCMCproposal_3d=1./(1.+abs(Means[Particle.m_CandidateLabel]-Image(Particle.Index)))
          ! L2: 1 / d^2
          ! MCMCproposal_3d=1./(1.+abs(Means[Particle.m_CandidateLabel]-Image(Particle.Index)))
          ! exp(-d^2)
          ! MCMCproposal_3d=exp((Means[Particle.CandidateLabel] - Image(Particle.Index)) * &
          ! &                   (Means[Particle.CandidateLabel] - Image(Particle.Index)))

          ! Length-prior driven proposal:
          MCMCproposal_3d=zero
          vFloating=.TRUE.

          tmplabels => labels_(coord(1)-1:coord(1)+1,coord(2)-1:coord(2)+1,coord(3)-1:coord(3)+1)
          cLabel=ABS(tmplabels(2,2,2))

          DO i=1,BG_ConnectivityType%NumberOfNeighbors
             ll=BG_ConnectivityType%NeighborsPoints(:,i)+2
             LabelNgh=ABS(tmplabels(ll(1),ll(2),ll(3)))
             IF (cLabel.NE.LabelNgh) THEN
                MCMCproposal_3d=MCMCproposal_3d+MCMClengthProposalMask(i)
                vFloating=.FALSE.
             ENDIF
          ENDDO

          IF (vFloating) THEN
             ! floating particles need a proposal > 0: We take half of the smallest
             ! element in the mask (hack:we assume the smallest element to be at
             ! position 0).
             MCMCproposal_3d = MCMClengthProposalMask(1)*half
          ENDIF
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCproposal_3d

        FUNCTION MCMCproposal__3d(tmplabels)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLy : zero,half,one,MCMCuseBiasedProposal
          USE ppm_rc_module_topologicalnumber, ONLY : BG_ConnectivityType
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, DIMENSION(:,:,:), POINTER :: tmplabels
          !!! Temporary label index around coordinate

          REAL(MK)                           :: MCMCproposal__3d
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER, DIMENSION(3) :: ll
          INTEGER               :: cLabel
          INTEGER               :: LabelNgh
          INTEGER               :: i

          LOGICAL :: vFloating
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          IF (.NOT.MCMCuseBiasedProposal) THEN
             MCMCproposal__3d=one
             RETURN
          ENDIF
          ! The choice of a data-driven proposal distribution depends a lot on
          ! the contrast in the image (?) and the current segmentation.

          ! L1: 1 / ( 1 + abs(d) )
          ! MCMCproposal__3d=1./(1.+abs(Means[Particle.m_CandidateLabel]-Image(Particle.Index)))
          ! L2: 1 / d^2
          ! MCMCproposal__3d=1./(1.+abs(Means[Particle.m_CandidateLabel]-Image(Particle.Index)))
          ! exp(-d^2)
          ! MCMCproposal__3d=exp((Means[Particle.CandidateLabel] - Image(Particle.Index)) * &
          ! &                   (Means[Particle.CandidateLabel] - Image(Particle.Index)))

          ! Length-prior driven proposal:
          MCMCproposal__3d=zero
          vFloating=.TRUE.

          cLabel=ABS(tmplabels(2,2,2))

          DO i=1,BG_ConnectivityType%NumberOfNeighbors
             ll=BG_ConnectivityType%NeighborsPoints(:,i)+2
             LabelNgh=ABS(tmplabels(ll(1),ll(2),ll(3)))
             IF (cLabel.NE.LabelNgh) THEN
                MCMCproposal__3d=MCMCproposal__3d+MCMClengthProposalMask(i)
                vFloating=.FALSE.
             ENDIF
          ENDDO

          IF (vFloating) THEN
             ! floating particles need a proposal > 0: We take half of the smallest
             ! element in the mask (hack:we assume the smallest element to be at
             ! position 0).
             MCMCproposal__3d = MCMClengthProposalMask(1)*half
          ENDIF
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCproposal__3d

        ! TOCHECK : The original implementation is wrong!!!!
        ! Returns false if topology is changed when applying this particle. This
        ! is based on the current state of the label image.
        ! TODO: this method should be changed to achieve full topological control.
        !      now the particle is rejected if it changes somehow the topology.
        LOGICAL FUNCTION MCMCIsParticleTopoValid_2d(coord,labels_,CandidateLabel)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLy : AllowFission,AllowHandles
          USE ppm_rc_module_topologicalnumber, ONLY : FGandBGTopoNbPairType, &
          &   TopologicalNumberFunction
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: coord
          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
          INTEGER,                             INTENT(IN   ) :: CandidateLabel
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(FGandBGTopoNbPairType) :: FGBGTNP

          INTEGER :: ContainerLabel
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------

          IF (.NOT.AllowFission.OR..NOT.AllowHandles) THEN
             IF (CandidateLabel.EQ.0) THEN
                ! Get the correct label to access the container of the particle
                ContainerLabel=ABS(labels_(coord(1),coord(2)))

                FGBGTNP = TopologicalNumberFunction%EvaluateFGTNOfLabelAtIndex(coord,labels_,ContainerLabel)

                IF (FGBGTNP%FGNumber.NE.1.OR.FGBGTNP%BGNumber.NE.1) THEN
                   MCMCIsParticleTopoValid_2d=.FALSE.
                   RETURN
                ENDIF
             ELSE
                ! if the both labels are not 0, we must calculate the
                ! topo numbers for the current and the candidate label.
                FGBGTNP = TopologicalNumberFunction%EvaluateFGTNOfLabelAtIndex(coord,labels_,CandidateLabel)

                IF (FGBGTNP%FGNumber.NE.1.OR.FGBGTNP%BGNumber.NE.1) THEN
                   MCMCIsParticleTopoValid_2d=.FALSE.
                   RETURN
                ENDIF

                ContainerLabel=ABS(labels_(coord(1),coord(2)))

                FGBGTNP = TopologicalNumberFunction%EvaluateFGTNOfLabelAtIndex(ContainerLabel)

                IF (FGBGTNP%FGNumber.NE.1.OR.FGBGTNP%BGNumber.NE.1) THEN
                   MCMCIsParticleTopoValid_2d=.FALSE.
                   RETURN
                ENDIF
             ENDIF
          ENDIF !.NOT.AllowFission.OR..NOT.AllowHandles

          MCMCIsParticleTopoValid_2d=.TRUE.
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCIsParticleTopoValid_2d

        LOGICAL FUNCTION MCMCIsParticleTopoValid__2d(tmplabels,CandidateLabel)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLy : AllowFission,AllowHandles
          USE ppm_rc_module_topologicalnumber, ONLY : FGandBGTopoNbPairType, &
          &   TopologicalNumberFunction
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, DIMENSION(:,:), POINTER       :: tmplabels
          INTEGER,                 INTENT(IN   ) :: CandidateLabel
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(FGandBGTopoNbPairType) :: FGBGTNP

          INTEGER :: ContainerLabel
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------

          IF (.NOT.AllowFission.OR..NOT.AllowHandles) THEN
             IF (CandidateLabel.EQ.0) THEN
                ! Get the correct label to access the container of the particle
                ContainerLabel=ABS(tmplabels(2,2))

                FGBGTNP = TopologicalNumberFunction%EvaluateFGTNOfLabelAtIndex(tmplabels,ContainerLabel)

                IF (FGBGTNP%FGNumber.NE.1.OR.FGBGTNP%BGNumber.NE.1) THEN
                   MCMCIsParticleTopoValid__2d=.FALSE.
                   RETURN
                ENDIF
             ELSE
                ! if the both labels are not 0, we must calculate the
                ! topo numbers for the current and the candidate label.
                FGBGTNP = TopologicalNumberFunction%EvaluateFGTNOfLabelAtIndex(tmplabels,CandidateLabel)

                IF (FGBGTNP%FGNumber.NE.1.OR.FGBGTNP%BGNumber.NE.1) THEN
                   MCMCIsParticleTopoValid__2d=.FALSE.
                   RETURN
                ENDIF

                ContainerLabel=ABS(tmplabels(2,2))

                FGBGTNP = TopologicalNumberFunction%EvaluateFGTNOfLabelAtIndex(ContainerLabel)

                IF (FGBGTNP%FGNumber.NE.1.OR.FGBGTNP%BGNumber.NE.1) THEN
                   MCMCIsParticleTopoValid__2d=.FALSE.
                   RETURN
                ENDIF
             ENDIF
          ENDIF !.NOT.AllowFission.OR..NOT.AllowHandles

          MCMCIsParticleTopoValid__2d=.TRUE.
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCIsParticleTopoValid__2d

        LOGICAL FUNCTION MCMCIsParticleTopoValid_3d(coord,labels_,CandidateLabel)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLy : AllowFission,AllowHandles
          USE ppm_rc_module_topologicalnumber, ONLY : FGandBGTopoNbPairType, &
          &   TopologicalNumberFunction
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: coord
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
          INTEGER,                               INTENT(IN   ) :: CandidateLabel
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(FGandBGTopoNbPairType) :: FGBGTNP

          INTEGER :: ContainerLabel

          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------

          IF (.NOT.AllowFission.OR..NOT.AllowHandles) THEN
             IF (CandidateLabel.EQ.0) THEN
                ! Get the correct label to access the container of the particle
                ContainerLabel=ABS(labels_(coord(1),coord(2),coord(3)))

                FGBGTNP = TopologicalNumberFunction%EvaluateFGTNOfLabelAtIndex(coord,labels_,ContainerLabel)

                IF (FGBGTNP%FGNumber.NE.1.OR.FGBGTNP%BGNumber.NE.1) THEN
                   MCMCIsParticleTopoValid_3d=.FALSE.
                   RETURN
                ENDIF
             ELSE
                FGBGTNP = TopologicalNumberFunction%EvaluateFGTNOfLabelAtIndex(coord,labels_,CandidateLabel)

                IF (FGBGTNP%FGNumber.NE.1.OR.FGBGTNP%BGNumber.NE.1) THEN
                   MCMCIsParticleTopoValid_3d=.FALSE.
                   RETURN
                ENDIF

                ! if the both labels are not 0, we must calculate the
                ! topo numbers for the current and the candidate label.
                ContainerLabel=ABS(labels_(coord(1),coord(2),coord(3)))

                FGBGTNP = TopologicalNumberFunction%EvaluateFGTNOfLabelAtIndex(ContainerLabel)

                IF (FGBGTNP%FGNumber.NE.1.OR.FGBGTNP%BGNumber.NE.1) THEN
                   MCMCIsParticleTopoValid_3d=.FALSE.
                   RETURN
                ENDIF
             ENDIF
          ENDIF !.NOT.AllowFission.OR..NOT.AllowHandles

          MCMCIsParticleTopoValid_3d=.TRUE.
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCIsParticleTopoValid_3d

        LOGICAL FUNCTION MCMCIsParticleTopoValid__3d(tmplabels,CandidateLabel)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLy : AllowFission,AllowHandles
          USE ppm_rc_module_topologicalnumber, ONLY : FGandBGTopoNbPairType, &
          &   TopologicalNumberFunction
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, DIMENSION(:,:,:), POINTER       :: tmplabels
          INTEGER,                   INTENT(IN   ) :: CandidateLabel
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(FGandBGTopoNbPairType) :: FGBGTNP

          INTEGER :: ContainerLabel

          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------

          IF (.NOT.AllowFission.OR..NOT.AllowHandles) THEN
             IF (CandidateLabel.EQ.0) THEN
                ! Get the correct label to access the container of the particle
                ContainerLabel=ABS(tmplabels(2,2,2))

                FGBGTNP = TopologicalNumberFunction%EvaluateFGTNOfLabelAtIndex(tmplabels,ContainerLabel)

                IF (FGBGTNP%FGNumber.NE.1.OR.FGBGTNP%BGNumber.NE.1) THEN
                   MCMCIsParticleTopoValid__3d=.FALSE.
                   RETURN
                ENDIF
             ELSE
                FGBGTNP = TopologicalNumberFunction%EvaluateFGTNOfLabelAtIndex(tmplabels,CandidateLabel)

                IF (FGBGTNP%FGNumber.NE.1.OR.FGBGTNP%BGNumber.NE.1) THEN
                   MCMCIsParticleTopoValid__3d=.FALSE.
                   RETURN
                ENDIF

                ! if the both labels are not 0, we must calculate the
                ! topo numbers for the current and the candidate label.
                ContainerLabel=ABS(tmplabels(2,2,2))

                FGBGTNP = TopologicalNumberFunction%EvaluateFGTNOfLabelAtIndex(ContainerLabel)

                IF (FGBGTNP%FGNumber.NE.1.OR.FGBGTNP%BGNumber.NE.1) THEN
                   MCMCIsParticleTopoValid__3d=.FALSE.
                   RETURN
                ENDIF
             ENDIF
          ENDIF !.NOT.AllowFission.OR..NOT.AllowHandles

          MCMCIsParticleTopoValid__3d=.TRUE.
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCIsParticleTopoValid__3d

        FUNCTION MCMCGetIndexFromEdgeDensity_2d(labels_,Nm) RESULT(index_)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          !!!
          USE ppm_module_data, ONLY : ppm_kind_int64

          USE ppm_rc_module_global, ONLY : FORBIDDEN
          USE ppm_rc_module_rnd, ONLY : GetImageDistrIndex
          USE ppm_rc_module_util, ONLY : id_gtl_2d
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: Nm
          INTEGER,             DIMENSION(2)                  :: index_
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER(ppm_kind_int64) :: lindex

          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          DO
             lindex=GetImageDistrIndex()
             CALL id_gtl_2d(Nm,lindex,index_)
             IF (labels_(index_(1),index_(2)).NE.FORBIDDEN) THEN
                EXIT
             ENDIF
          ENDDO
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCGetIndexFromEdgeDensity_2d

        FUNCTION MCMCGetIndexFromEdgeDensity_3d(labels_,Nm) RESULT(index_)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          !!!
          USE ppm_module_data, ONLY : ppm_kind_int64

          USE ppm_rc_module_global, ONLY : FORBIDDEN
          USE ppm_rc_module_rnd, ONLY : GetImageDistrIndex
          USE ppm_rc_module_util, ONLY : id_gtl_3d
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: Nm
          INTEGER,             DIMENSION(3)                    :: index_
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER(ppm_kind_int64) :: lindex

          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          DO
             lindex=GetImageDistrIndex()
             CALL id_gtl_3d(Nm,lindex,index_)
             IF (labels_(index_(1),index_(2),index_(3)).NE.FORBIDDEN) THEN
                EXIT
             ENDIF
          ENDDO
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCGetIndexFromEdgeDensity_3d

        FUNCTION MCMCgetRegularParticlesAtIndex_2d(coord,Nm,labels_,MCMCParticles) RESULT(nsize)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLY : FORBIDDEN
          USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType, &
          &   BG_ConnectivityType
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: coord
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: Nm
          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
          TYPE(MCMCParticle),  DIMENSION(:),   INTENT(INOUT) :: MCMCParticles
          INTEGER                                            :: nsize
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          REAL(MK) :: Proposal

          INTEGER, DIMENSION(:,:), POINTER :: tmplabels
          INTEGER, DIMENSION(2)            :: ll,ld
          INTEGER                          :: cLabel
          INTEGER                          :: LabelNgh
          INTEGER                          :: i

          LOGICAL :: NotParentInserted

          nsize=0

          IF (labels_(coord(1),coord(2)).EQ.FORBIDDEN) RETURN

          tmplabels => labels_(coord(1)-1:coord(1)+1,coord(2)-1:coord(2)+1)
          cLabel=ABS(tmplabels(2,2))

          Proposal=MCMCproposal(tmplabels)

          SELECT CASE (cLabel)
          CASE (0)
             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ld=FG_ConnectivityType%NeighborsPoints(:,i)+coord
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+2
                LabelNgh=ABS(tmplabels(ll(1),ll(2)))
                IF (LabelNgh.NE.0.AND.LabelNgh.NE.FORBIDDEN) THEN
                !!! here there should be a particle since two neighboring pixels
                !!! have different labels.
                   !!! FG labels have a daughter placed at this spot.
                   nsize=nsize+1
                   MCMCParticles(nsize)=MCMCParticle(LabelNgh,Proposal)
                ENDIF
             ENDDO
          CASE DEFAULT
             NotParentInserted=.NOT.ANY(coord.LT.1.OR.coord.GT.Nm)
             !!! If coord is outside the domain, we should not put a parent
             !!! there
             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ld=FG_ConnectivityType%NeighborsPoints(:,i)+coord
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+2
                LabelNgh=ABS(tmplabels(ll(1),ll(2)))
                IF (LabelNgh.NE.cLabel.AND.LabelNgh.NE.FORBIDDEN) THEN
                !!! here there should be a particle since two neighboring pixels
                !!! have different labels.
                   IF (LabelNgh.NE.0) THEN
                      !!! FG labels have a daughter placed at this spot.
                      nsize=nsize+1
                      MCMCParticles(nsize)=MCMCParticle(LabelNgh,Proposal)
                   ENDIF
                   IF (NotParentInserted) THEN
                   !!! this is a non-zero pixel with different neighbors,
                   !!! hence there must be a mother in the list:
                      NotParentInserted=.FALSE.
                      nsize=nsize+1
                      MCMCParticles(nsize)=MCMCParticle(0,Proposal)
                   ENDIF
                ENDIF
             ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors

             !!! Check the BG neighborhood now if we need to insert a parent.
             IF (NotParentInserted) THEN
                DO i=1,BG_ConnectivityType%NumberOfNeighbors
                   ll=BG_ConnectivityType%NeighborsPoints(:,i)+2
                   LabelNgh=ABS(tmplabels(ll(1),ll(2)))
                   IF (cLabel.NE.LabelNgh) THEN
                      !!! This is a FG pixel with a neighbor of a different label.
                      !!! finally, insert a parent particle
                      nsize=nsize+1
                      MCMCParticles(nsize)=MCMCParticle(0,Proposal)
                      EXIT
                   ENDIF
                ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors
             ENDIF !NotParentInserted
          END SELECT !cLabel
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCgetRegularParticlesAtIndex_2d

        FUNCTION MCMCgetRegularParticlesAtIndex__2d(coord1,coord2,Nm,tmplabels,MCMCParticles) RESULT(nsize)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLY : FORBIDDEN
          USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType, &
          &   BG_ConnectivityType
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,                            INTENT(IN   ) :: coord1
          INTEGER,                            INTENT(IN   ) :: coord2
          INTEGER,            DIMENSION(:),   INTENT(IN   ) :: Nm
          INTEGER,            DIMENSION(:,:), POINTER       :: tmplabels
          TYPE(MCMCParticle), DIMENSION(:),   INTENT(INOUT) :: MCMCParticles
          INTEGER                                           :: nsize
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          REAL(MK) :: Proposal

          INTEGER, DIMENSION(2) :: ll,ld
          INTEGER               :: cLabel
          INTEGER               :: LabelNgh
          INTEGER               :: i

          LOGICAL :: NotParentInserted

          nsize=0

          IF (tmplabels(2,2).EQ.FORBIDDEN) RETURN

          cLabel=ABS(tmplabels(2,2))

          proposal=MCMCproposal(tmplabels)

          SELECT CASE (cLabel)
          CASE (0)
             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ld(1)=FG_ConnectivityType%NeighborsPoints(1,i)+coord1
                ld(2)=FG_ConnectivityType%NeighborsPoints(2,i)+coord2
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+2
                LabelNgh=ABS(tmplabels(ll(1),ll(2)))
                IF (LabelNgh.NE.0.AND.LabelNgh.NE.FORBIDDEN) THEN
                !!! here there should be a particle since two neighboring pixels
                !!! have different labels.
                   !!! FG labels have a daughter placed at this spot.
                   nsize=nsize+1
                   MCMCParticles(nsize)=MCMCParticle(LabelNgh,proposal)
                ENDIF
             ENDDO
          CASE DEFAULT
             NotParentInserted=.NOT.(coord1.LT.1.OR.coord2.LT.1.OR.coord1.GT.Nm(1).OR.coord2.GT.Nm(2))
             !!! If coord is outside the domain, we should not put a parent
             !!! there
             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ld(1)=FG_ConnectivityType%NeighborsPoints(1,i)+coord1
                ld(2)=FG_ConnectivityType%NeighborsPoints(2,i)+coord2
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+2
                LabelNgh=ABS(tmplabels(ll(1),ll(2)))
                IF (LabelNgh.NE.cLabel.AND.LabelNgh.NE.FORBIDDEN) THEN
                !!! here there should be a particle since two neighboring pixels
                !!! have different labels.
                   IF (LabelNgh.NE.0) THEN
                      !!! FG labels have a daughter placed at this spot.
                      nsize=nsize+1
                      MCMCParticles(nsize)=MCMCParticle(LabelNgh,proposal)
                   ENDIF
                   IF (NotParentInserted) THEN
                   !!! this is a non-zero pixel with different neighbors,
                   !!! hence there must be a mother in the list:
                      NotParentInserted=.FALSE.
                      nsize=nsize+1
                      MCMCParticles(nsize)=MCMCParticle(0,proposal)
                   ENDIF
                ENDIF
             ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors

             !!! Check the BG neighborhood now if we need to insert a parent.
             IF (NotParentInserted) THEN
                DO i=1,BG_ConnectivityType%NumberOfNeighbors
                   ll=BG_ConnectivityType%NeighborsPoints(:,i)+2
                   LabelNgh=ABS(tmplabels(ll(1),ll(2)))
                   IF (cLabel.NE.LabelNgh) THEN
                      !!! This is a FG pixel with a neighbor of a different label.
                      !!! finally, insert a parent particle
                      nsize=nsize+1
                      MCMCParticles(nsize)=MCMCParticle(0,proposal)
                      EXIT
                   ENDIF
                ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors
             ENDIF !NotParentInserted
          END SELECT !cLabel
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCgetRegularParticlesAtIndex__2d

        FUNCTION MCMCgetRegularParticlesAtIndex_3d(coord,Nm,labels_,MCMCParticles) RESULT(nsize)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLY : FORBIDDEN
          USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType, &
          &   BG_ConnectivityType
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: coord
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: Nm
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
          TYPE(MCMCParticle),  DIMENSION(:),     INTENT(INOUT) :: MCMCParticles
          INTEGER                                              :: nsize
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          REAL(MK) :: proposal

          INTEGER, DIMENSION(:,:,:), POINTER :: tmplabels
          INTEGER, DIMENSION(3)              :: ll,ld
          INTEGER                            :: cLabel
          INTEGER                            :: LabelNgh
          INTEGER                            :: i

          LOGICAL :: NotParentInserted

          nsize=0

          IF (labels_(coord(1),coord(2),coord(3)).EQ.FORBIDDEN) RETURN

          tmplabels => labels_(coord(1)-1:coord(1)+1,coord(2)-1:coord(2)+1,coord(3)-1:coord(3)+1)
          cLabel=ABS(tmplabels(2,2,2))

          proposal=MCMCproposal(tmplabels)

          SELECT CASE (cLabel)
          CASE (0)
             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ld=FG_ConnectivityType%NeighborsPoints(:,i)+coord
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+2
                LabelNgh=ABS(tmplabels(ll(1),ll(2),ll(3)))
                IF (LabelNgh.NE.0.AND.LabelNgh.NE.FORBIDDEN) THEN
                !!! here there should be a particle since two neighboring pixels
                !!! have different labels.
                   !!! FG labels have a daughter placed at this spot.
                   nsize=nsize+1
                   MCMCParticles(nsize)=MCMCParticle(LabelNgh,proposal)
                ENDIF
             ENDDO
          CASE DEFAULT
             NotParentInserted=.NOT.ANY(coord.LT.1.OR.coord.GT.Nm)
             !!! If coord is outside the domain, we should not put a parent
             !!! there
             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ld=FG_ConnectivityType%NeighborsPoints(:,i)+coord
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+2
                LabelNgh=ABS(tmplabels(ll(1),ll(2),ll(3)))
                IF (LabelNgh.NE.cLabel.AND.LabelNgh.NE.FORBIDDEN) THEN
                !!! here there should be a particle since two neighboring pixels
                !!! have different labels.
                   IF (LabelNgh.NE.0) THEN
                      !!! FG labels have a daughter placed at this spot.
                      nsize=nsize+1
                      MCMCParticles(nsize)=MCMCParticle(LabelNgh,proposal)
                   ENDIF
                   IF (NotParentInserted) THEN
                   !!! this is a non-zero pixel with different neighbors,
                   !!! hence there must be a mother in the list:
                      NotParentInserted=.FALSE.
                      nsize=nsize+1
                      MCMCParticles(nsize)=MCMCParticle(0,proposal)
                   ENDIF
                ENDIF
             ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors

             !!! Check the BG neighborhood now if we need to insert a parent.
             IF (NotParentInserted) THEN
                DO i=1,BG_ConnectivityType%NumberOfNeighbors
                   ll=BG_ConnectivityType%NeighborsPoints(:,i)+2
                   LabelNgh=ABS(tmplabels(ll(1),ll(2),ll(3)))
                   IF (cLabel.NE.LabelNgh) THEN
                      !!! This is a FG pixel with a neighbor of a different label.
                      !!! finally, insert a parent particle
                      nsize=nsize+1
                      MCMCParticles(nsize)=MCMCParticle(0,proposal)
                      EXIT
                   ENDIF
                ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors
             ENDIF !NotParentInserted
          END SELECT !cLabel
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCgetRegularParticlesAtIndex_3d

        FUNCTION MCMCgetRegularParticlesAtIndex__3d(coord1,coord2,coord3,Nm,tmplabels,MCMCParticles) RESULT(nsize)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLY : FORBIDDEN
          USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType, &
          &   BG_ConnectivityType
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,                              INTENT(IN   ) :: coord1
          INTEGER,                              INTENT(IN   ) :: coord2
          INTEGER,                              INTENT(IN   ) :: coord3
          INTEGER,            DIMENSION(:),     INTENT(IN   ) :: Nm
          INTEGER,            DIMENSION(:,:,:), POINTER       :: tmplabels
          TYPE(MCMCParticle), DIMENSION(:),     INTENT(INOUT) :: MCMCParticles
          INTEGER                                             :: nsize
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          REAL(MK) :: proposal

          INTEGER, DIMENSION(3) :: ll,ld
          INTEGER               :: cLabel
          INTEGER               :: LabelNgh
          INTEGER               :: i

          LOGICAL :: NotParentInserted

          nsize=0

          IF (tmplabels(2,2,2).EQ.FORBIDDEN) RETURN

          cLabel=ABS(tmplabels(2,2,2))

          proposal=MCMCproposal(tmplabels)

          SELECT CASE (cLabel)
          CASE (0)
             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ld(1)=FG_ConnectivityType%NeighborsPoints(1,i)+coord1
                ld(2)=FG_ConnectivityType%NeighborsPoints(2,i)+coord2
                ld(3)=FG_ConnectivityType%NeighborsPoints(3,i)+coord3
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+2
                LabelNgh=ABS(tmplabels(ll(1),ll(2),ll(3)))
                IF (LabelNgh.NE.0.AND.LabelNgh.NE.FORBIDDEN) THEN
                !!! here there should be a particle since two neighboring pixels
                !!! have different labels.
                   !!! FG labels have a daughter placed at this spot.
                   nsize=nsize+1
                   MCMCParticles(nsize)=MCMCParticle(LabelNgh,proposal)
                ENDIF
             ENDDO
          CASE DEFAULT
             NotParentInserted=.NOT.(coord1.LT.1.OR.coord2.LT.1.OR.coord3.LT.1.OR. &
             &                       coord1.GT.Nm(1).OR.coord2.GT.Nm(2).OR.coord3.GT.Nm(3))
             !!! If coord is outside the domain, we should not put a parent
             !!! there
             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ld(1)=FG_ConnectivityType%NeighborsPoints(1,i)+coord1
                ld(2)=FG_ConnectivityType%NeighborsPoints(2,i)+coord2
                ld(3)=FG_ConnectivityType%NeighborsPoints(3,i)+coord3
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+2
                LabelNgh=ABS(tmplabels(ll(1),ll(2),ll(3)))
                IF (LabelNgh.NE.cLabel.AND.LabelNgh.NE.FORBIDDEN) THEN
                !!! here there should be a particle since two neighboring pixels
                !!! have different labels.
                   IF (LabelNgh.NE.0) THEN
                      !!! FG labels have a daughter placed at this spot.
                      nsize=nsize+1
                      MCMCParticles(nsize)=MCMCParticle(LabelNgh,proposal)
                   ENDIF
                   IF (NotParentInserted) THEN
                   !!! this is a non-zero pixel with different neighbors,
                   !!! hence there must be a mother in the list:
                      NotParentInserted=.FALSE.
                      nsize=nsize+1
                      MCMCParticles(nsize)=MCMCParticle(0,proposal)
                   ENDIF
                ENDIF
             ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors

             !!! Check the BG neighborhood now if we need to insert a parent.
             IF (NotParentInserted) THEN
                DO i=1,BG_ConnectivityType%NumberOfNeighbors
                   ll=BG_ConnectivityType%NeighborsPoints(:,i)+2
                   LabelNgh=ABS(tmplabels(ll(1),ll(2),ll(3)))
                   IF (cLabel.NE.LabelNgh) THEN
                      !!! This is a FG pixel with a neighbor of a different label.
                      !!! finally, insert a parent particle
                      nsize=nsize+1
                      MCMCParticles(nsize)=MCMCParticle(0,proposal)
                      EXIT
                   ENDIF
                ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors
             ENDIF !NotParentInserted
          END SELECT !cLabel
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCgetRegularParticlesAtIndex__3d

        FUNCTION MCMCgetParticlesInBGNeighborhood_2d(coord,Nm,labels_,MCMCParticles,MCMCParticlesCoords) RESULT(nsize)
          !!! The method returns a list of all particles in the neighborhood INCLUSIVE
          !!! the particles at coord, which comes at the end!
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_error_fatal
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error

          USE ppm_rc_module_global, ONLY : AllowFission,AllowHandles
          USE ppm_rc_module_topologicalnumber, ONLY : BG_ConnectivityType
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,             DIMENSION(:),                INTENT(IN   ) :: coord
          INTEGER,             DIMENSION(:),                INTENT(IN   ) :: Nm
          INTEGER, CONTIGUOUS, DIMENSION(:,:),              POINTER       :: labels_
          TYPE(MCMCParticle),  DIMENSION(:),   ALLOCATABLE, INTENT(INOUT) :: MCMCParticles
          INTEGER,             DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: MCMCParticlesCoords
          INTEGER                                                         :: nsize
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle), DIMENSION(:), ALLOCATABLE :: tmpParticles
          TYPE(MCMCParticle), DIMENSION(9)              :: MCMCParticles_

          REAL(MK) :: proposal

          INTEGER, DIMENSION(:,:), ALLOCATABLE :: tmpcoords
          INTEGER, DIMENSION(:,:), POINTER     :: tmplabels_
          INTEGER, DIMENSION(:,:), POINTER     :: tmplabels
          INTEGER, DIMENSION(2)                :: ll,ld
          INTEGER                              :: nsize_
          INTEGER                              :: i,j
          INTEGER                              :: info

          CHARACTER(LEN=*), PARAMETER :: caller="MCMCgetParticlesInBGNeighborhood"
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          nsize=0

          tmplabels => labels_(coord(1)-2:coord(1)+2,coord(2)-2:coord(2)+2)

          DO i=1,BG_ConnectivityType%NumberOfNeighbors
             ld=BG_ConnectivityType%NeighborsPoints(:,i)+coord

             ll=BG_ConnectivityType%NeighborsPoints(:,i)+3

             tmplabels => tmplabels_(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1)

             nsize_=MCMCgetRegularParticlesAtIndex(ld(1),ld(2),Nm,tmplabels,MCMCParticles_)

             IF (nsize_.GT.0) THEN
                IF (.NOT.AllowFission.AND..NOT.AllowHandles) THEN
                   !!! filter the points to meet topological constraints:
                   DO j=nsize_,1,-1
                      !!! Removes all non-simple points from the particle set
                      IF (.NOT.MCMCIsParticleTopoValid(tmplabels,MCMCParticles_(j)%candlabel)) THEN
                         MCMCParticles_(j)=MCMCParticles_(nsize_)
                         nsize_=nsize_-1
                      ENDIF
                   ENDDO
                ENDIF
             ENDIF !nsize_.GT.0

             IF (nsize_.GT.0) THEN
                ALLOCATE(tmpParticles(nsize+nsize_),tmpcoords(2,nsize+nsize_),STAT=info)
                or_fail_alloc("Failed to allocate tmpParticles & tmpcoords!", &
                & ppm_error=ppm_error_fatal,exit_point=no)

                FORALL (j=1:nsize)
                   tmpParticles(j)=MCMCParticles(j)
                   tmpcoords(:,j)=MCMCParticlesCoords(:,j)
                END FORALL
                FORALL (j=1:nsize_)
                   tmpParticles(nsize+j)=MCMCParticles_(j)
                   tmpcoords(:,nsize+j)=ld
                END FORALL

                CALL MOVE_ALLOC(tmpParticles,MCMCParticles)
                CALL MOVE_ALLOC(tmpcoords,MCMCParticlesCoords)

                nsize=nsize+nsize_
             ENDIF
          ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors

          !!! Insert particle A in the set for B such that it is possible
          !!! that we have a single particle move.
          tmplabels => tmplabels_(2:4,2:4)
          proposal=MCMCproposal(tmplabels)

          ALLOCATE(tmpParticles(nsize+1),tmpcoords(2,nsize+1),STAT=info)
          or_fail_alloc("Failed to allocate tmpParticles & tmpcoords!", &
          & ppm_error=ppm_error_fatal,exit_point=no)

          FORALL (j=1:nsize)
             tmpParticles(j)=MCMCParticles(j)
             tmpcoords(:,j)=MCMCParticlesCoords(:,j)
          END FORALL

          tmpParticles(nsize+1)=MCMCParticle(-1,proposal)
          tmpcoords(:,nsize+1)=coord

          nsize=nsize+1

          CALL MOVE_ALLOC(tmpParticles,MCMCParticles)
          CALL MOVE_ALLOC(tmpcoords,MCMCParticlesCoords)

          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCgetParticlesInBGNeighborhood_2d

        FUNCTION MCMCgetParticlesInBGNeighborhood_3d(coord,Nm,labels_,MCMCParticles,MCMCParticlesCoords) RESULT(nsize)
          !!! The method returns a list of all particles in the neighborhood INCLUSIVE
          !!! the particles at coord, which comes at the end!
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_error_fatal
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error

          USE ppm_rc_module_global, ONLY : AllowFission,AllowHandles
          USE ppm_rc_module_topologicalnumber, ONLY : BG_ConnectivityType
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,             DIMENSION(:),                INTENT(IN   ) :: coord
          INTEGER,             DIMENSION(:),                INTENT(IN   ) :: Nm
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:),            POINTER       :: labels_
          TYPE(MCMCParticle),  DIMENSION(:),   ALLOCATABLE, INTENT(INOUT) :: MCMCParticles
          INTEGER,             DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: MCMCParticlesCoords
          INTEGER                                                         :: nsize
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle), DIMENSION(:), ALLOCATABLE :: tmpParticles
          TYPE(MCMCParticle), DIMENSION(27)             :: MCMCParticles_

          REAL(MK) :: proposal

          INTEGER, DIMENSION(:,:),   ALLOCATABLE :: tmpcoords
          INTEGER, DIMENSION(:,:,:), POINTER     :: tmplabels_
          INTEGER, DIMENSION(:,:,:), POINTER     :: tmplabels
          INTEGER, DIMENSION(3)                  :: ll,ld
          INTEGER                                :: nsize_
          INTEGER                                :: i,j
          INTEGER                                :: info

          CHARACTER(LEN=*), PARAMETER :: caller="MCMCgetParticlesInBGNeighborhood"
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          nsize=0

          tmplabels => labels_(coord(1)-2:coord(1)+2,coord(2)-2:coord(2)+2,coord(3)-2:coord(3)+2)

          DO i=1,BG_ConnectivityType%NumberOfNeighbors
             ld=BG_ConnectivityType%NeighborsPoints(:,i)+coord

             ll=BG_ConnectivityType%NeighborsPoints(:,i)+3

             tmplabels => tmplabels_(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1,ll(3)-1:ll(3)+1)

             nsize_=MCMCgetRegularParticlesAtIndex(ld(1),ld(2),ld(3),Nm,tmplabels,MCMCParticles_)

             IF (nsize_.GT.0) THEN
                IF (.NOT.AllowFission.AND..NOT.AllowHandles) THEN
                !!! filter the points to meet topological constraints:
                   DO j=nsize_,1,-1
                      !!! Removes all non-simple points from the particle set
                      IF (.NOT.MCMCIsParticleTopoValid(tmplabels,MCMCParticles_(j)%candlabel)) THEN
                         MCMCParticles_(j)=MCMCParticles_(nsize_)
                         nsize_=nsize_-1
                      ENDIF
                   ENDDO
                ENDIF
             ENDIF !nsize_.GT.0

             IF (nsize_.GT.0) THEN
                ALLOCATE(tmpParticles(nsize+nsize_),tmpcoords(3,nsize+nsize_),STAT=info)
                or_fail_alloc("Failed to allocate tmpParticles & tmpcoords!", &
                & ppm_error=ppm_error_fatal,exit_point=no)

                FORALL (j=1:nsize)
                   tmpParticles(j)=MCMCParticles(j)
                   tmpcoords(:,j)=MCMCParticlesCoords(:,j)
                END FORALL
                FORALL (j=1:nsize_)
                   tmpParticles(nsize+j)=MCMCParticles_(j)
                   tmpcoords(:,nsize+j)=ld
                END FORALL

                CALL MOVE_ALLOC(tmpParticles,MCMCParticles)
                CALL MOVE_ALLOC(tmpcoords,MCMCParticlesCoords)

                nsize=nsize+nsize_
             ENDIF
          ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors

          !!! Insert particle A in the set for B such that it is possible
          !!! that we have a single particle move.
          tmplabels => tmplabels_(2:4,2:4,2:4)
          proposal=MCMCproposal(tmplabels)

          ALLOCATE(tmpParticles(nsize+1),tmpcoords(3,nsize+1),STAT=info)
          or_fail_alloc("Failed to allocate tmpParticles & tmpcoords!", &
          & ppm_error=ppm_error_fatal,exit_point=no)

          FORALL (j=1:nsize)
             tmpParticles(j)=MCMCParticles(j)
             tmpcoords(:,j)=MCMCParticlesCoords(:,j)
          END FORALL

          tmpParticles(nsize+1)=MCMCParticle(-1,proposal)
          tmpcoords(:,nsize+1)=coord

          nsize=nsize+1

          CALL MOVE_ALLOC(tmpParticles,MCMCParticles)
          CALL MOVE_ALLOC(tmpcoords,MCMCParticlesCoords)

          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCgetParticlesInBGNeighborhood_3d

        FUNCTION MCMCgetParticlesInFGNeighborhood_2d(coord,Nm,labels_,MCMCParticles,MCMCParticlesCoords) RESULT(nsize)
          !!! The method returns a list of all particles in the neighborhood INCLUSIVE
          !!! the particles at coord, which comes at the end!
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_error_fatal
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error

          USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,             DIMENSION(:),                INTENT(IN   ) :: coord
          INTEGER,             DIMENSION(:),                INTENT(IN   ) :: Nm
          INTEGER, CONTIGUOUS, DIMENSION(:,:),              POINTER       :: labels_
          TYPE(MCMCParticle),  DIMENSION(:),   ALLOCATABLE, INTENT(INOUT) :: MCMCParticles
          INTEGER,             DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: MCMCParticlesCoords
          INTEGER                                                         :: nsize
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle), DIMENSION(:), ALLOCATABLE :: tmpParticles
          TYPE(MCMCParticle), DIMENSION(9)              :: MCMCParticles_

          REAL(MK) :: proposal

          INTEGER, DIMENSION(:,:), ALLOCATABLE :: tmpcoords
          INTEGER, DIMENSION(:,:), POINTER     :: tmplabels_
          INTEGER, DIMENSION(:,:), POINTER     :: tmplabels
          INTEGER, DIMENSION(2)                :: ll,ld
          INTEGER                              :: nsize_
          INTEGER                              :: i,j
          INTEGER                              :: info

          CHARACTER(LEN=*), PARAMETER :: caller="MCMCgetParticlesInFGNeighborhood"
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          nsize=0

          tmplabels => labels_(coord(1)-2:coord(1)+2,coord(2)-2:coord(2)+2)

          DO i=1,FG_ConnectivityType%NumberOfNeighbors
             ld=FG_ConnectivityType%NeighborsPoints(:,i)+coord

             ll=FG_ConnectivityType%NeighborsPoints(:,i)+3

             tmplabels => tmplabels_(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1)

             nsize_=MCMCgetRegularParticlesAtIndex(ld(1),ld(2),Nm,tmplabels,MCMCParticles_)

             IF (nsize_.GT.0) THEN
                ALLOCATE(tmpParticles(nsize+nsize_),tmpcoords(2,nsize+nsize_),STAT=info)
                or_fail_alloc("Failed to allocate tmpParticles & tmpcoords!", &
                & ppm_error=ppm_error_fatal,exit_point=no)

                FORALL (j=1:nsize)
                   tmpParticles(j)=MCMCParticles(j)
                   tmpcoords(:,j)=MCMCParticlesCoords(:,j)
                END FORALL
                FORALL (j=1:nsize_)
                   tmpParticles(nsize+j)=MCMCParticles_(j)
                   tmpcoords(:,nsize+j)=ld
                END FORALL

                CALL MOVE_ALLOC(tmpParticles,MCMCParticles)
                CALL MOVE_ALLOC(tmpcoords,MCMCParticlesCoords)

                nsize=nsize+nsize_
             ENDIF
          ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors

          !!! Insert particle such that it is possible
          !!! that we have a single particle move.
          tmplabels => tmplabels_(2:4,2:4)
          proposal=MCMCproposal(tmplabels)

          ALLOCATE(tmpParticles(nsize+1),tmpcoords(2,nsize+1),STAT=info)
          or_fail_alloc("Failed to allocate tmpParticles & tmpcoords!", &
          & ppm_error=ppm_error_fatal,exit_point=no)

          FORALL (j=1:nsize)
             tmpParticles(j)=MCMCParticles(j)
             tmpcoords(:,j)=MCMCParticlesCoords(:,j)
          END FORALL

          tmpParticles(nsize+1)=MCMCParticle(-1,proposal)
          tmpcoords(:,nsize+1)=coord

          nsize=nsize+1

          CALL MOVE_ALLOC(tmpParticles,MCMCParticles)
          CALL MOVE_ALLOC(tmpcoords,MCMCParticlesCoords)
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCgetParticlesInFGNeighborhood_2d

        FUNCTION MCMCgetParticlesInFGNeighborhood_3d(coord,Nm,labels_,MCMCParticles,MCMCParticlesCoords) RESULT(nsize)
          !!! The method returns a list of all particles in the neighborhood INCLUSIVE
          !!! the particles at coord, which comes at the end!
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_error_fatal
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error

          USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,             DIMENSION(:),                INTENT(IN   ) :: coord
          INTEGER,             DIMENSION(:),                INTENT(IN   ) :: Nm
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:),            POINTER       :: labels_
          TYPE(MCMCParticle),  DIMENSION(:),   ALLOCATABLE, INTENT(INOUT) :: MCMCParticles
          INTEGER,             DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: MCMCParticlesCoords
          INTEGER                                                         :: nsize
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle), DIMENSION(:), ALLOCATABLE :: tmpParticles
          TYPE(MCMCParticle), DIMENSION(27)             :: MCMCParticles_

          REAL(MK) :: proposal

          INTEGER, DIMENSION(:,:),   ALLOCATABLE :: tmpcoords
          INTEGER, DIMENSION(:,:,:), POINTER     :: tmplabels_
          INTEGER, DIMENSION(:,:,:), POINTER     :: tmplabels
          INTEGER, DIMENSION(3)                  :: ll,ld
          INTEGER                                :: nsize_
          INTEGER                                :: i,j
          INTEGER                                :: info

          CHARACTER(LEN=*), PARAMETER :: caller="MCMCgetParticlesInFGNeighborhood"
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          nsize=0

          tmplabels => labels_(coord(1)-2:coord(1)+2,coord(2)-2:coord(2)+2,coord(3)-2:coord(3)+2)

          DO i=1,FG_ConnectivityType%NumberOfNeighbors
             ld=FG_ConnectivityType%NeighborsPoints(:,i)+coord

             ll=FG_ConnectivityType%NeighborsPoints(:,i)+3

             tmplabels => tmplabels_(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1,ll(3)-1:ll(3)+1)

             nsize_=MCMCgetRegularParticlesAtIndex(ld(1),ld(2),ld(3),Nm,tmplabels,MCMCParticles_)

             IF (nsize_.GT.0) THEN
                ALLOCATE(tmpParticles(nsize+nsize_),tmpcoords(3,nsize+nsize_),STAT=info)
                or_fail_alloc("Failed to allocate tmpParticles & tmpcoords!", &
                & ppm_error=ppm_error_fatal,exit_point=no)

                FORALL (j=1:nsize)
                   tmpParticles(j)=MCMCParticles(j)
                   tmpcoords(:,j)=MCMCParticlesCoords(:,j)
                END FORALL
                FORALL (j=1:nsize_)
                   tmpParticles(nsize+j)=MCMCParticles_(j)
                   tmpcoords(:,nsize+j)=ld
                END FORALL

                CALL MOVE_ALLOC(tmpParticles,MCMCParticles)
                CALL MOVE_ALLOC(tmpcoords,MCMCParticlesCoords)

                nsize=nsize+nsize_
             ENDIF
          ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors

          !!! Insert particle such that it is possible
          !!! that we have a single particle move.
          tmplabels => tmplabels_(2:4,2:4,2:4)
          proposal=MCMCproposal(tmplabels)

          ALLOCATE(tmpParticles(nsize+1),tmpcoords(3,nsize+1),STAT=info)
          or_fail_alloc("Failed to allocate tmpParticles & tmpcoords!", &
          & ppm_error=ppm_error_fatal,exit_point=no)

          FORALL (j=1:nsize)
             tmpParticles(j)=MCMCParticles(j)
             tmpcoords(:,j)=MCMCParticlesCoords(:,j)
          END FORALL

          tmpParticles(nsize+1)=MCMCParticle(-1,proposal)
          tmpcoords(:,nsize+1)=coord

          nsize=nsize+1

          CALL MOVE_ALLOC(tmpParticles,MCMCParticles)
          CALL MOVE_ALLOC(tmpcoords,MCMCParticlesCoords)
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCgetParticlesInFGNeighborhood_3d

        FUNCTION MCMCInsertCandidatesToContainers_2d(ipatch,coord1,coord2,MCMCParticle_,CurrentLabel,DoRecord,ParticleReplaced) RESULT (info)
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error,ppm_err_argument, &
          &   ppm_err_sub_failed

          USE ppm_rc_module_global, ONLY : htable,htable_null
          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list,MCMCParticleInContainerHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: ipatch
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          INTEGER,            INTENT(IN   ) :: CurrentLabel
          LOGICAL,            INTENT(IN   ) :: DoRecord
          LOGICAL,            INTENT(  OUT) :: ParticleReplaced
          INTEGER                           :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle
          TYPE(MCMCHistoryParticle) :: tmpHistoryParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0
          REAL(MK)              :: ProposalDiff

          INTEGER :: label_region

          CHARACTER(LEN=ppm_char) :: caller='MCMCInsertCandidatesToContainers'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          ParticleReplaced=.FALSE.

          IF (MCMCParticle_%candlabel.EQ.0) THEN
             label_region=htable%search(CurrentLabel)
             IF (label_region.EQ.htable_null) THEN
                fail("The region does not exist!!!",ppm_error=ppm_error_fatal)
             ENDIF

             tmpParticle=MCMCparents(label_region)%search(coord1,coord2)

             CALL MCMCparents(label_region)%insert(coord1,coord2,MCMCParticle_,info)
             or_fail("MCMCparents(label_region)%insert")

             IF (tmpParticle%candlabel.EQ.-1) THEN
                MCMCparentsProposalNormalizer(label_region)= &
                & MCMCparentsProposalNormalizer(label_region)+MCMCParticle_%proposal
                MCMCTotalNormalizer=MCMCTotalNormalizer+REAL(MCMCParticle_%proposal,ppm_kind_double)
             ELSE
                !!! this is a replacement:
                ProposalDiff=MCMCParticle_%proposal-tmpParticle%proposal

                MCMCparentsProposalNormalizer(label_region)= &
                & MCMCparentsProposalNormalizer(label_region)+ProposalDiff
                MCMCTotalNormalizer=MCMCTotalNormalizer+REAL(ProposalDiff,ppm_kind_double)

                ParticleReplaced=.TRUE.
             ENDIF
          ELSE
             label_region=htable%search(MCMCParticle_%candlabel)
             IF (label_region.EQ.htable_null) THEN
                fail("The region does not exist!!!",ppm_error=ppm_error_fatal)
             ENDIF

             tmpParticle=MCMCchildren(label_region)%search(coord1,coord2)

             CALL MCMCchildren(label_region)%insert(coord1,coord2,MCMCParticle_,info)
             or_fail("MCMCchildren(label_region)%insert")

             IF (tmpParticle%candlabel.EQ.-1) THEN
                MCMCchildrenProposalNormalizer(label_region)= &
                & MCMCchildrenProposalNormalizer(label_region)+MCMCParticle_%proposal
                MCMCTotalNormalizer=MCMCTotalNormalizer+REAL(MCMCParticle_%proposal,ppm_kind_double)
             ELSE
                !!! this is a replacement:
                ProposalDiff=MCMCParticle_%proposal-tmpParticle%proposal

                MCMCchildrenProposalNormalizer(label_region)= &
                & MCMCchildrenProposalNormalizer(label_region)+ProposalDiff
                MCMCTotalNormalizer=MCMCTotalNormalizer+REAL(ProposalDiff,ppm_kind_double)

                ParticleReplaced=.TRUE.
             ENDIF
          ENDIF

          IF (DoRecord) THEN
             !!! Note that the order here is important: first insert the particle
             !!! that gets replaced. When restoring the state afterwards, the
             !!! particle history is iterated in reverse order.
             IF (ParticleReplaced) THEN
                tmpHistoryParticle=MCMCHistoryParticle(tmpParticle%candlabel, &
                & CurrentLabel,tmpParticle%proposal,.FALSE.)
                ALLOCATE(seed,STAT=info)
                or_fail_alloc("seed")
                CALL seed%add(coord1,coord2,tmpHistoryParticle)
                CALL MCMCParticleInContainerHistory(ipatch)%push(seed,info)
                or_fail("MCMCParticleInContainerHistory(ipatch)%push")
             ENDIF

             tmpHistoryParticle=MCMCHistoryParticle(MCMCParticle_%candlabel, &
             & CurrentLabel,MCMCParticle_%proposal,.TRUE.)
             ALLOCATE(seed,STAT=info)
             or_fail_alloc("seed")
             CALL seed%add(coord1,coord2,tmpHistoryParticle)
             CALL MCMCParticleInContainerHistory(ipatch)%push(seed,info)
             or_fail("MCMCParticleInContainerHistory(ipatch)%push")
          ENDIF

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCInsertCandidatesToContainers_2d

        FUNCTION MCMCInsertCandidatesToContainers_3d(ipatch,coord1,coord2,coord3,MCMCParticle_,CurrentLabel,DoRecord,ParticleReplaced) RESULT(info)
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error,ppm_err_argument, &
          &   ppm_err_sub_failed

          USE ppm_rc_module_global, ONLY : htable,htable_null
          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list,MCMCParticleInContainerHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: ipatch
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          INTEGER,            INTENT(IN   ) :: coord3
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          INTEGER,            INTENT(IN   ) :: CurrentLabel
          LOGICAL,            INTENT(IN   ) :: DoRecord
          LOGICAL,            INTENT(  OUT) :: ParticleReplaced
          INTEGER                           :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle

          TYPE(MCMCHistoryParticle) :: tmpHistoryParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0
          REAL(MK)              :: ProposalDiff

          INTEGER :: label_region

          CHARACTER(LEN=ppm_char) :: caller='MCMCInsertCandidatesToContainers'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          ParticleReplaced=.FALSE.

          IF (MCMCParticle_%candlabel.EQ.0) THEN
             label_region=htable%search(CurrentLabel)
             IF (label_region.EQ.htable_null) THEN
                fail("The region does not exist!!!",ppm_error=ppm_error_fatal)
             ENDIF

             tmpParticle=MCMCparents(label_region)%search(coord1,coord2,coord3)

             CALL MCMCparents(label_region)%insert(coord1,coord2,coord3,MCMCParticle_,info)
             or_fail("MCMCparents(label_region)%insert")

             IF (tmpParticle%candlabel.EQ.-1) THEN
                MCMCparentsProposalNormalizer(label_region)= &
                & MCMCparentsProposalNormalizer(label_region)+MCMCParticle_%proposal
                MCMCTotalNormalizer=MCMCTotalNormalizer+REAL(MCMCParticle_%proposal,ppm_kind_double)
             ELSE
                !!! this is a replacement:
                ProposalDiff=MCMCParticle_%proposal-tmpParticle%proposal

                MCMCparentsProposalNormalizer(label_region)= &
                & MCMCparentsProposalNormalizer(label_region)+ProposalDiff
                MCMCTotalNormalizer=MCMCTotalNormalizer+REAL(ProposalDiff,ppm_kind_double)

                ParticleReplaced=.TRUE.
             ENDIF
          ELSE
             label_region=htable%search(MCMCParticle_%candlabel)
             IF (label_region.EQ.htable_null) THEN
                fail("The region does not exist!!!",ppm_error=ppm_error_fatal)
             ENDIF

             tmpParticle=MCMCchildren(label_region)%search(coord1,coord2,coord3)

             CALL MCMCchildren(label_region)%insert(coord1,coord2,coord3,MCMCParticle_,info)
             or_fail("MCMCchildren(label_region)%insert")

             IF (tmpParticle%candlabel.EQ.-1) THEN
                MCMCchildrenProposalNormalizer(label_region)= &
                & MCMCchildrenProposalNormalizer(label_region)+MCMCParticle_%proposal
                MCMCTotalNormalizer=MCMCTotalNormalizer+REAL(MCMCParticle_%proposal,ppm_kind_double)
             ELSE
                !!! this is a replacement:
                ProposalDiff=MCMCParticle_%proposal-tmpParticle%proposal

                MCMCchildrenProposalNormalizer(label_region)= &
                & MCMCchildrenProposalNormalizer(label_region)+ProposalDiff
                MCMCTotalNormalizer=MCMCTotalNormalizer+REAL(ProposalDiff,ppm_kind_double)

                ParticleReplaced=.TRUE.
             ENDIF
          ENDIF

          IF (DoRecord) THEN
             !!! Note that the order here is important: first insert the particle
             !!! that gets replaced. When restoring the state afterwards, the
             !!! particle history is iterated in reverse order.
             IF (ParticleReplaced) THEN
                tmpHistoryParticle=MCMCHistoryParticle(tmpParticle%candlabel, &
                & CurrentLabel,tmpParticle%proposal,.FALSE.)
                ALLOCATE(seed,STAT=info)
                or_fail_alloc("seed")
                CALL seed%add(coord1,coord2,coord3,tmpHistoryParticle)
                CALL MCMCParticleInContainerHistory(ipatch)%push(seed,info)
                or_fail("MCMCParticleInContainerHistory(ipatch)%push")
             ENDIF

             tmpHistoryParticle=MCMCHistoryParticle(MCMCParticle_%candlabel, &
             & CurrentLabel,MCMCParticle_%proposal,.TRUE.)
             ALLOCATE(seed,STAT=info)
             or_fail_alloc("seed")
             CALL seed%add(coord1,coord2,coord3,tmpHistoryParticle)
             CALL MCMCParticleInContainerHistory(ipatch)%push(seed,info)
             or_fail("MCMCParticleInContainerHistory(ipatch)%push")
          ENDIF

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCInsertCandidatesToContainers_3d

        FUNCTION MCMCEraseCandidatesFromContainers_2d(ipatch,coord1,coord2,MCMCParticle_,CurrentLabel,DoRecord) RESULT (info)
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error,ppm_err_argument, &
          &   ppm_err_sub_failed
          USE ppm_module_inl_hash, ONLY : htable_null

          USE ppm_rc_module_global, ONLY : htable
          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list,MCMCParticleInContainerHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: ipatch
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          INTEGER,            INTENT(IN   ) :: CurrentLabel
          LOGICAL,            INTENT(IN   ) :: DoRecord
          INTEGER                           :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle
          TYPE(MCMCHistoryParticle) :: tmpHistoryParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0

          INTEGER :: label_region

          CHARACTER(LEN=ppm_char) :: caller='MCMCEraseCandidatesFromContainers'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          IF (MCMCParticle_%candlabel.EQ.0) THEN
             label_region=htable%search(CurrentLabel)
             IF (label_region.EQ.htable_null) THEN
                fail("The region does not exist!!!",ppm_error=ppm_error_fatal)
             ENDIF

             tmpParticle=MCMCparents(label_region)%search(coord1,coord2)

             IF (tmpParticle%candlabel.GT.-1) THEN
                !!! Particle exists
                CALL MCMCparents(label_region)%remove(coord1,coord2,info,.TRUE.)
                or_fail("MCMCparents(label_region)%remove")

                MCMCparentsProposalNormalizer(label_region)= &
                & MCMCparentsProposalNormalizer(label_region)-tmpParticle%proposal
                MCMCTotalNormalizer=MCMCTotalNormalizer-REAL(tmpParticle%proposal,ppm_kind_double)
             ENDIF !tmpParticle%candlabel.GT.-1
          ELSE
             label_region=htable%search(MCMCParticle_%candlabel)
             IF (label_region.EQ.htable_null) THEN
                fail("The region does not exist!!!",ppm_error=ppm_error_fatal)
             ENDIF

             tmpParticle=MCMCchildren(label_region)%search(coord1,coord2)

             IF (tmpParticle%candlabel.GT.-1) THEN
                !!! Particle exists
                CALL MCMCchildren(label_region)%remove(coord1,coord2,info,.TRUE.)
                or_fail("MCMCchildren(label_region)%remove")

                MCMCchildrenProposalNormalizer(label_region)= &
                & MCMCchildrenProposalNormalizer(label_region)-tmpParticle%proposal
                MCMCTotalNormalizer=MCMCTotalNormalizer-REAL(tmpParticle%proposal,ppm_kind_double)
             ENDIF
          ENDIF

          IF (DoRecord) THEN
             IF (tmpParticle%candlabel.GT.-1) THEN
                tmpHistoryParticle=MCMCHistoryParticle(tmpParticle%candlabel, &
                & CurrentLabel,tmpParticle%proposal,.FALSE.)
                ALLOCATE(seed,STAT=info)
                or_fail_alloc("seed")
                CALL seed%add(coord1,coord2,tmpHistoryParticle)
                CALL MCMCParticleInContainerHistory(ipatch)%push(seed,info)
                or_fail("MCMCParticleInContainerHistory(ipatch)%push")
             ENDIF
          ENDIF

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCEraseCandidatesFromContainers_2d

        FUNCTION MCMCEraseCandidatesFromContainers_3d(ipatch,coord1,coord2,coord3,MCMCParticle_,CurrentLabel,DoRecord) RESULT (info)
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error,ppm_err_argument, &
          &   ppm_err_sub_failed
          USE ppm_module_inl_hash, ONLY : htable_null

          USE ppm_rc_module_global, ONLY : htable
          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list,MCMCParticleInContainerHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: ipatch
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          INTEGER,            INTENT(IN   ) :: coord3
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          INTEGER,            INTENT(IN   ) :: CurrentLabel
          LOGICAL,            INTENT(IN   ) :: DoRecord
          INTEGER                           :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle
          TYPE(MCMCHistoryParticle) :: tmpHistoryParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0

          INTEGER :: label_region

          CHARACTER(LEN=ppm_char) :: caller='MCMCEraseCandidatesFromContainers'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          IF (MCMCParticle_%candlabel.EQ.0) THEN
             label_region=htable%search(CurrentLabel)
             IF (label_region.EQ.htable_null) THEN
                fail("The region does not exist!!!",ppm_error=ppm_error_fatal)
             ENDIF

             tmpParticle=MCMCparents(label_region)%search(coord1,coord2,coord3)

             IF (tmpParticle%candlabel.GT.-1) THEN
                !!! Particle exists
                CALL MCMCparents(label_region)%remove(coord1,coord2,coord3,info,.TRUE.)
                or_fail("MCMCparents(label_region)%remove")

                MCMCparentsProposalNormalizer(label_region)= &
                & MCMCparentsProposalNormalizer(label_region)-tmpParticle%proposal
                MCMCTotalNormalizer=MCMCTotalNormalizer-REAL(tmpParticle%proposal,ppm_kind_double)
             ENDIF !tmpParticle%candlabel.GT.-1
          ELSE
             label_region=htable%search(MCMCParticle_%candlabel)
             IF (label_region.EQ.htable_null) THEN
                fail("The region does not exist!!!",ppm_error=ppm_error_fatal)
             ENDIF

             tmpParticle=MCMCchildren(label_region)%search(coord1,coord2,coord3)

             IF (tmpParticle%candlabel.GT.-1) THEN
                !!! Particle exists
                CALL MCMCchildren(label_region)%remove(coord1,coord2,coord3,info,.TRUE.)
                or_fail("MCMCchildren(label_region)%remove")

                MCMCchildrenProposalNormalizer(label_region)= &
                & MCMCchildrenProposalNormalizer(label_region)-tmpParticle%proposal
                MCMCTotalNormalizer=MCMCTotalNormalizer-REAL(tmpParticle%proposal,ppm_kind_double)
             ENDIF
          ENDIF

          IF (DoRecord) THEN
             IF (tmpParticle%candlabel.GT.-1) THEN
                tmpHistoryParticle=MCMCHistoryParticle(tmpParticle%candlabel, &
                & CurrentLabel,tmpParticle%proposal,.FALSE.)
                ALLOCATE(seed,STAT=info)
                or_fail_alloc("seed")
                CALL seed%add(coord1,coord2,coord3,tmpHistoryParticle)
                CALL MCMCParticleInContainerHistory(ipatch)%push(seed,info)
                or_fail("MCMCParticleInContainerHistory(ipatch)%push")
             ENDIF
          ENDIF

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCEraseCandidatesFromContainers_3d

        LOGICAL FUNCTION MCMCEraseFloatingParticle_2d(ipatch,coord1,coord2,MCMCParticle_,DoRecord)
          !!! Erase a floating particle only if it exists. If it does not exist, the
          !!! return value will be false.
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_alloc,ppm_err_sub_failed

          USE ppm_rc_module_global, ONLY : small,ten
          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list, &
          &   MCMCFloatingParticleInContainerHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: ipatch
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          LOGICAL,            INTENT(IN   ) :: DoRecord
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle
          TYPE(MCMCHistoryParticle) :: tmpHistoryParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0

          INTEGER :: info

          CHARACTER(LEN=ppm_char) :: caller="MCMCEraseFloatingParticle"
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          tmpParticle=MCMCFloatingParticles%search(coord1,coord2)
          IF (tmpParticle%candlabel.EQ.-1) THEN
             !!! there is no floating particle at this position
             MCMCEraseFloatingParticle_2d=.FALSE.
          ELSE
             IF (DoRecord) THEN
                tmpHistoryParticle=MCMCHistoryParticle(MCMCParticle_%candlabel, &
                & 0,MCMCParticle_%proposal,.FALSE.)
                ALLOCATE(seed,STAT=info)
                or_fail_alloc("seed",ppm_error=ppm_error_fatal)
                CALL seed%add(coord1,coord2,tmpHistoryParticle)
                CALL MCMCFloatingParticleInContainerHistory(ipatch)%push(seed,info)
                or_fail("MCMCFloatingParticleInContainerHistory(ipatch)%push",ppm_error=ppm_error_fatal)
             ENDIF

             IF ((tmpParticle%proposal-MCMCParticle_%proposal).GT.small*ten) THEN
                !!! the particle has still some proposal left
                MCMCFloatingParticlesProposalNormalizer= &
                & MCMCFloatingParticlesProposalNormalizer-REAL(MCMCParticle_%proposal,ppm_kind_double)
                tmpParticle%proposal=tmpParticle%proposal-MCMCParticle_%proposal
                CALL MCMCFloatingParticles%insert(coord1,coord2,tmpParticle,info)
                or_fail("MCMCFloatingParticles%insert",ppm_error=ppm_error_fatal)
             ELSE
                !!! the particle gets deleted. We only remove the amount from the
                !!! normalizer that has been stored in the set.
                MCMCFloatingParticlesProposalNormalizer= &
                & MCMCFloatingParticlesProposalNormalizer-REAL(tmpParticle%proposal,ppm_kind_double)
                CALL MCMCFloatingParticles%remove(coord1,coord2,info,.TRUE.)
                or_fail("MCMCFloatingParticles%remove",ppm_error=ppm_error_fatal)
             ENDIF
             MCMCEraseFloatingParticle_2d=.TRUE.
          ENDIF
        9999 CONTINUE
          CALL substop(caller,t0,info)
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCEraseFloatingParticle_2d

        LOGICAL FUNCTION MCMCEraseFloatingParticle__2d(ipatch,coord1,coord2,searchedParticle,MCMCParticle_,DoRecord)
          !!! Erase a floating particle only if it exists. The particle at coords has been searched as searchedParticle.
          !!! If it does not exist, the return value will be false.
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_alloc,ppm_err_sub_failed

          USE ppm_rc_module_global, ONLY : small,ten
          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list, &
          &   MCMCFloatingParticleInContainerHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: ipatch
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          TYPE(MCMCParticle), INTENT(INOUT) :: searchedParticle
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          LOGICAL,            INTENT(IN   ) :: DoRecord
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCHistoryParticle) :: tmpHistoryParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0

          INTEGER :: info

          CHARACTER(LEN=ppm_char) :: caller="MCMCEraseFloatingParticle"
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          IF (searchedParticle%candlabel.EQ.-1) THEN
             !!! there is no floating particle at this position
             MCMCEraseFloatingParticle__2d=.FALSE.
          ELSE
             IF (DoRecord) THEN
                tmpHistoryParticle=MCMCHistoryParticle(MCMCParticle_%candlabel, &
                & 0,MCMCParticle_%proposal,.FALSE.)
                ALLOCATE(seed,STAT=info)
                or_fail_alloc("seed",ppm_error=ppm_error_fatal)
                CALL seed%add(coord1,coord2,tmpHistoryParticle)
                CALL MCMCFloatingParticleInContainerHistory(ipatch)%push(seed,info)
                or_fail("MCMCFloatingParticleInContainerHistory(ipatch)%push",ppm_error=ppm_error_fatal)
             ENDIF

             IF ((searchedParticle%proposal-MCMCParticle_%proposal).GT.small*ten) THEN
                !!! the particle has still some proposal left
                MCMCFloatingParticlesProposalNormalizer= &
                & MCMCFloatingParticlesProposalNormalizer-REAL(MCMCParticle_%proposal,ppm_kind_double)
                searchedParticle%proposal=searchedParticle%proposal-MCMCParticle_%proposal
                CALL MCMCFloatingParticles%insert(coord1,coord2,searchedParticle,info)
                or_fail("MCMCFloatingParticles%insert",ppm_error=ppm_error_fatal)
             ELSE
                !!! the particle gets deleted. We only remove the amount from the
                !!! normalizer that has been stored in the set.
                MCMCFloatingParticlesProposalNormalizer= &
                & MCMCFloatingParticlesProposalNormalizer-REAL(searchedParticle%proposal,ppm_kind_double)
                CALL MCMCFloatingParticles%remove(coord1,coord2,info,.TRUE.)
                or_fail("MCMCFloatingParticles%remove",ppm_error=ppm_error_fatal)
             ENDIF
             MCMCEraseFloatingParticle__2d=.TRUE.
          ENDIF
        9999 CONTINUE
          CALL substop(caller,t0,info)
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCEraseFloatingParticle__2d

        LOGICAL FUNCTION MCMCEraseFloatingParticle_3d(ipatch,coord1,coord2,coord3,MCMCParticle_,DoRecord)
          !!! Erase a floating particle only if it exists. If it does not exist, the
          !!! return value will be false.
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_alloc,ppm_err_sub_failed

          USE ppm_rc_module_global, ONLY : small,ten
          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list, &
          &   MCMCFloatingParticleInContainerHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: ipatch
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          INTEGER,            INTENT(IN   ) :: coord3
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          LOGICAL,            INTENT(IN   ) :: DoRecord
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle
          TYPE(MCMCHistoryParticle) :: tmpHistoryParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0

          INTEGER :: info

          CHARACTER(LEN=ppm_char) :: caller="MCMCEraseFloatingParticle"
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          tmpParticle=MCMCFloatingParticles%search(coord1,coord2,coord3)
          IF (tmpParticle%candlabel.EQ.-1) THEN
             !!! there is no floating particle at this position
             MCMCEraseFloatingParticle_3d=.FALSE.
          ELSE
             IF (DoRecord) THEN
                tmpHistoryParticle=MCMCHistoryParticle(MCMCParticle_%candlabel, &
                & 0,MCMCParticle_%proposal,.FALSE.)
                ALLOCATE(seed,STAT=info)
                or_fail_alloc("seed",ppm_error=ppm_error_fatal)
                CALL seed%add(coord1,coord2,coord3,tmpHistoryParticle)
                CALL MCMCFloatingParticleInContainerHistory(ipatch)%push(seed,info)
                or_fail("MCMCFloatingParticleInContainerHistory(ipatch)%push",ppm_error=ppm_error_fatal)
             ENDIF

             IF ((tmpParticle%proposal-MCMCParticle_%proposal).GT.small*ten) THEN
                !!! the particle has still some proposal left
                MCMCFloatingParticlesProposalNormalizer= &
                & MCMCFloatingParticlesProposalNormalizer-REAL(MCMCParticle_%proposal,ppm_kind_double)
                tmpParticle%proposal=tmpParticle%proposal-MCMCParticle_%proposal
                CALL MCMCFloatingParticles%insert(coord1,coord2,coord3,tmpParticle,info)
                or_fail("MCMCFloatingParticles%insert",ppm_error=ppm_error_fatal)
             ELSE
                !!! the particle gets deleted. We only remove the amount from the
                !!! normalizer that has been stored in the set.
                MCMCFloatingParticlesProposalNormalizer= &
                & MCMCFloatingParticlesProposalNormalizer-REAL(tmpParticle%proposal,ppm_kind_double)
                CALL MCMCFloatingParticles%remove(coord1,coord2,coord3,info,.TRUE.)
                or_fail("MCMCFloatingParticles%remove",ppm_error=ppm_error_fatal)
             ENDIF
             MCMCEraseFloatingParticle_3d=.TRUE.
          ENDIF
        9999 CONTINUE
          CALL substop(caller,t0,info)
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCEraseFloatingParticle_3d

        LOGICAL FUNCTION MCMCEraseFloatingParticle__3d(ipatch,coord1,coord2,coord3,searchedParticle,MCMCParticle_,DoRecord)
          !!! Erase a floating particle only if it exists. The particle at coords has been searched as searchedParticle.
          !!! If it does not exist, the return value will be false.
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_alloc,ppm_err_sub_failed

          USE ppm_rc_module_global, ONLY : small,ten
          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list, &
          &   MCMCFloatingParticleInContainerHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: ipatch
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          INTEGER,            INTENT(IN   ) :: coord3
          TYPE(MCMCParticle), INTENT(INOUT) :: searchedParticle
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          LOGICAL,            INTENT(IN   ) :: DoRecord
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCHistoryParticle) :: tmpHistoryParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0

          INTEGER :: info

          CHARACTER(LEN=ppm_char) :: caller="MCMCEraseFloatingParticle"
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          IF (searchedParticle%candlabel.EQ.-1) THEN
             !!! there is no floating particle at this position
             MCMCEraseFloatingParticle__3d=.FALSE.
          ELSE
             IF (DoRecord) THEN
                tmpHistoryParticle=MCMCHistoryParticle(MCMCParticle_%candlabel, &
                & 0,MCMCParticle_%proposal,.FALSE.)
                ALLOCATE(seed,STAT=info)
                or_fail_alloc("seed",ppm_error=ppm_error_fatal)
                CALL seed%add(coord1,coord2,coord3,tmpHistoryParticle)
                CALL MCMCFloatingParticleInContainerHistory(ipatch)%push(seed,info)
                or_fail("MCMCFloatingParticleInContainerHistory(ipatch)%push",ppm_error=ppm_error_fatal)
             ENDIF

             IF ((searchedParticle%proposal-MCMCParticle_%proposal).GT.small*ten) THEN
                !!! the particle has still some proposal left
                MCMCFloatingParticlesProposalNormalizer= &
                & MCMCFloatingParticlesProposalNormalizer-REAL(MCMCParticle_%proposal,ppm_kind_double)
                searchedParticle%proposal=searchedParticle%proposal-MCMCParticle_%proposal
                CALL MCMCFloatingParticles%insert(coord1,coord2,coord3,searchedParticle,info)
                or_fail("MCMCFloatingParticles%insert",ppm_error=ppm_error_fatal)
             ELSE
                !!! the particle gets deleted. We only remove the amount from the
                !!! normalizer that has been stored in the set.
                MCMCFloatingParticlesProposalNormalizer= &
                & MCMCFloatingParticlesProposalNormalizer-REAL(searchedParticle%proposal,ppm_kind_double)
                CALL MCMCFloatingParticles%remove(coord1,coord2,coord3,info,.TRUE.)
                or_fail("MCMCFloatingParticles%remove",ppm_error=ppm_error_fatal)
             ENDIF
             MCMCEraseFloatingParticle__3d=.TRUE.
          ENDIF
        9999 CONTINUE
          CALL substop(caller,t0,info)
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCEraseFloatingParticle__3d

        LOGICAL FUNCTION MCMCInsertFloatingParticle_2d(ipatch,coord1,coord2,MCMCParticle_,DoRecord,MCMCParticleExist)
          !!! Insert a floating particle only if it didn't exist. If it existed, the
          !!! return value will be false.
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_alloc,ppm_err_sub_failed

          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list, &
          &   MCMCFloatingParticleInContainerHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: ipatch
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          LOGICAL,            INTENT(IN   ) :: DoRecord
          LOGICAL, OPTIONAL,  INTENT(IN   ) :: MCMCParticleExist
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle
          TYPE(MCMCHistoryParticle) :: tmpHistoryParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0

          INTEGER :: info

          CHARACTER(LEN=ppm_char) :: caller="MCMCInsertFloatingParticle"
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          IF (PRESENT(MCMCParticleExist)) THEN
             IF (MCMCParticleExist) THEN
                !!! the particle already exists
                MCMCInsertFloatingParticle_2d=.FALSE.
                GOTO 9999
             ELSE
                !!! there is no floating particle at this position
                CALL MCMCFloatingParticles%insert(coord1,coord2,MCMCParticle_,info)
                or_fail("MCMCFloatingParticles%insert",ppm_error=ppm_error_fatal)

                IF (DoRecord) THEN
                   tmpHistoryParticle=MCMCHistoryParticle(MCMCParticle_%candlabel, &
                   & 0,MCMCParticle_%proposal,.TRUE.)
                   ALLOCATE(seed,STAT=info)
                   or_fail_alloc("seed",ppm_error=ppm_error_fatal)
                   CALL seed%add(coord1,coord2,tmpHistoryParticle)
                   CALL MCMCFloatingParticleInContainerHistory(ipatch)%push(seed,info)
                   or_fail("MCMCFloatingParticleInContainerHistory(ipatch)%push",ppm_error=ppm_error_fatal)
                ENDIF

                MCMCFloatingParticlesProposalNormalizer= &
                & MCMCFloatingParticlesProposalNormalizer+REAL(MCMCParticle_%proposal,ppm_kind_double)

                MCMCInsertFloatingParticle_2d=.TRUE.
                GOTO 9999
             ENDIF
          ENDIF

          tmpParticle=MCMCFloatingParticles%search(coord1,coord2)
          IF (tmpParticle%candlabel.EQ.-1) THEN
             !!! there is no floating particle at this position
             CALL MCMCFloatingParticles%insert(coord1,coord2,MCMCParticle_,info)
             or_fail("MCMCFloatingParticles%insert",ppm_error=ppm_error_fatal)

             IF (DoRecord) THEN
                tmpHistoryParticle=MCMCHistoryParticle(MCMCParticle_%candlabel, &
                & 0,MCMCParticle_%proposal,.TRUE.)
                ALLOCATE(seed,STAT=info)
                or_fail_alloc("seed",ppm_error=ppm_error_fatal)
                CALL seed%add(coord1,coord2,tmpHistoryParticle)
                CALL MCMCFloatingParticleInContainerHistory(ipatch)%push(seed,info)
                or_fail("MCMCFloatingParticleInContainerHistory(ipatch)%push",ppm_error=ppm_error_fatal)
             ENDIF

             MCMCFloatingParticlesProposalNormalizer= &
             & MCMCFloatingParticlesProposalNormalizer+REAL(MCMCParticle_%proposal,ppm_kind_double)

             MCMCInsertFloatingParticle_2d=.TRUE.
          ELSE
             !!! the particle already exists
             MCMCInsertFloatingParticle_2d=.FALSE.
          ENDIF
        9999 CONTINUE
          CALL substop(caller,t0,info)
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCInsertFloatingParticle_2d

        LOGICAL FUNCTION MCMCInsertFloatingParticle_3d(ipatch,coord1,coord2,coord3,MCMCParticle_,DoRecord,MCMCParticleExist)
          !!! Insert a floating particle only if it didn't exist. If it existed, the
          !!! return value will be false.
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_alloc,ppm_err_sub_failed

          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list, &
          &   MCMCFloatingParticleInContainerHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: ipatch
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          INTEGER,            INTENT(IN   ) :: coord3
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          LOGICAL,            INTENT(IN   ) :: DoRecord
          LOGICAL, OPTIONAL,  INTENT(IN   ) :: MCMCParticleExist
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle
          TYPE(MCMCHistoryParticle) :: tmpHistoryParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0

          INTEGER :: info

          CHARACTER(LEN=ppm_char) :: caller="MCMCInsertFloatingParticle"
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          IF (PRESENT(MCMCParticleExist)) THEN
             IF (MCMCParticleExist) THEN
                !!! the particle already exists
                MCMCInsertFloatingParticle_3d=.FALSE.
                GOTO 9999
             ELSE
                !!! there is no floating particle at this position
                CALL MCMCFloatingParticles%insert(coord1,coord2,coord3,MCMCParticle_,info)
                or_fail("MCMCFloatingParticles%insert",ppm_error=ppm_error_fatal)

                IF (DoRecord) THEN
                   tmpHistoryParticle=MCMCHistoryParticle(MCMCParticle_%candlabel, &
                   & 0,MCMCParticle_%proposal,.TRUE.)
                   ALLOCATE(seed,STAT=info)
                   or_fail_alloc("seed",ppm_error=ppm_error_fatal)
                   CALL seed%add(coord1,coord2,coord3,tmpHistoryParticle)
                   CALL MCMCFloatingParticleInContainerHistory(ipatch)%push(seed,info)
                   or_fail("MCMCFloatingParticleInContainerHistory(ipatch)%push",ppm_error=ppm_error_fatal)
                ENDIF

                MCMCFloatingParticlesProposalNormalizer= &
                & MCMCFloatingParticlesProposalNormalizer+REAL(MCMCParticle_%proposal,ppm_kind_double)

                MCMCInsertFloatingParticle_3d=.TRUE.
                GOTO 9999
             ENDIF
          ENDIF

          tmpParticle=MCMCFloatingParticles%search(coord1,coord2,coord3)
          IF (tmpParticle%candlabel.EQ.-1) THEN
             !!! there is no floating particle at this position
             CALL MCMCFloatingParticles%insert(coord1,coord2,coord3,MCMCParticle_,info)
             or_fail("MCMCFloatingParticles%insert",ppm_error=ppm_error_fatal)

             IF (DoRecord) THEN
                tmpHistoryParticle=MCMCHistoryParticle(MCMCParticle_%candlabel, &
                & 0,MCMCParticle_%proposal,.TRUE.)
                ALLOCATE(seed,STAT=info)
                or_fail_alloc("seed",ppm_error=ppm_error_fatal)
                CALL seed%add(coord1,coord2,coord3,tmpHistoryParticle)
                CALL MCMCFloatingParticleInContainerHistory(ipatch)%push(seed,info)
                or_fail("MCMCFloatingParticleInContainerHistory(ipatch)%push",ppm_error=ppm_error_fatal)
             ENDIF

             MCMCFloatingParticlesProposalNormalizer= &
             & MCMCFloatingParticlesProposalNormalizer+REAL(MCMCParticle_%proposal,ppm_kind_double)

             MCMCInsertFloatingParticle_3d=.TRUE.
          ELSE
             !!! the particle already exists
             MCMCInsertFloatingParticle_3d=.FALSE.
          ENDIF
        9999 CONTINUE
          CALL substop(caller,t0,info)
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCInsertFloatingParticle_3d

        LOGICAL FUNCTION MCMCOffBoundarySample_2d(Growth,LabelIn)
          !!!  For the insertion and deletion of floating particles, we use a
          !!!  constant energy, therefore the MH ratio be equal to the FBR.
          !!!
          !!!  Note that we also could use the energy as if we would do the job as
          !!!  proposal distribution, but we would have to calculate the backward
          !!!  energy as well.
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_sub_failed
          USE ppm_module_mesh_typedef, ONLY : ppm_t_subpatch_

          USE ppm_rc_module_global, ONLY : half,mesh,labels
          USE ppm_rc_module_rnd, ONLY : SaruGetRealVariate
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          LOGICAL, INTENT(IN   ) :: Growth
          INTEGER, INTENT(IN   ) :: LabelIn
          !!! Absolute value label from regionsLabel
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle

          CLASS(ppm_t_subpatch_), POINTER :: sbpitr

          REAL(ppm_kind_double) :: t0

          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER :: labels_2d
          INTEGER,             DIMENSION(:),   POINTER :: Nm
          INTEGER,             DIMENSION(2)            :: index_
          INTEGER                                      :: ipatch
          INTEGER                                      :: info

          LOGICAL :: Child

          CHARACTER(LEN=ppm_char) :: caller='MCMCOffBoundarySample'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)
          !!! For the insertion and deletion of floating particles, we use a
          !!! constant energy, therefore the MH ratio be equal to the FBR.
          !!! Note that we also could use the energy as if we would do the job as
          !!! proposal distribution, but we would have to calculate the backward
          !!! energy as well.

          NULLIFY(labels_2d)

          sbpitr => mesh%subpatch%begin()
          ipatch=1
          DO WHILE (ASSOCIATED(sbpitr))
             NM => sbpitr%nnodes

             CALL sbpitr%get_field(labels,labels_2d,info)
             or_fail("Failed to get field labels_2d data.",ppm_error=ppm_error_fatal)

             index_=MCMCGetIndexFromEdgeDensity(labels_2d,Nm)
             !!! get index from Edge Density

             tmpParticle=MCMCFloatingParticles%search(index_(1),index_(2))
             !!! search for the particle at this index

             IF (tmpParticle%candlabel.EQ.-1) THEN
                !!! Particle does not exist
                IF (Growth) THEN
                   !!! insert particle

                   !!! the target-distribution probability ratio
                   Child=SaruGetRealVariate().GT.half

                   tmpParticle%candlabel=MERGE(LabelIn,0,Child)
                   !!! TO CHECK
                   !!! TO DO
                   !!! This one does not exist in the original code (to compute the proposal)
                   tmpParticle%proposal=MCMCproposal(index_,labels_2d)

                   MCMCOffBoundarySample_2d= &
                   &  MCMCInsertFloatingParticle(ipatch,index_(1),index_(2),tmpParticle,.FALSE.,.FALSE.)
                   GOTO 9999
                ELSE
                   !!! reject
                   MCMCOffBoundarySample_2d=.FALSE.
                   GOTO 9999
                ENDIF
             ELSE
                !!! Particle exists
                IF (Growth) THEN
                   !!! reject
                   MCMCOffBoundarySample_2d=.FALSE.
                   GOTO 9999
                ELSE
                   !!! .NOT.Growth.AND.Particle exists
                   MCMCOffBoundarySample_2d= &
                   & MCMCEraseFloatingParticle(ipatch,index_(1),index_(2),tmpParticle,tmpParticle,.FALSE.)
                   GOTO 9999
                ENDIF
             ENDIF

             sbpitr => mesh%subpatch%next()
             ipatch=ipatch+1
          ENDDO

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCOffBoundarySample_2d

        LOGICAL FUNCTION MCMCOffBoundarySample_3d(Growth,LabelIn)
          !!!  For the insertion and deletion of floating particles, we use a
          !!!  constant energy, therefore the MH ratio be equal to the FBR.
          !!!
          !!!  Note that we also could use the energy as if we would do the job as
          !!!  proposal distribution, but we would have to calculate the backward
          !!!  energy as well.
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_sub_failed
          USE ppm_module_mesh_typedef, ONLY : ppm_t_subpatch_

          USE ppm_rc_module_global, ONLY : half,mesh,labels
          USE ppm_rc_module_rnd, ONLY : SaruGetRealVariate
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          LOGICAL, INTENT(IN   ) :: Growth
          INTEGER, INTENT(IN   ) :: LabelIn
          !!! Absolute value label from regionsLabel
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle

          CLASS(ppm_t_subpatch_), POINTER :: sbpitr

          REAL(ppm_kind_double) :: t0

          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: labels_3d
          INTEGER,             DIMENSION(:),     POINTER :: Nm
          INTEGER,             DIMENSION(3)              :: index_
          INTEGER                                        :: ipatch
          INTEGER                                        :: info

          LOGICAL :: Child

          CHARACTER(LEN=ppm_char) :: caller='MCMCOffBoundarySample'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)
          !!! For the insertion and deletion of floating particles, we use a
          !!! constant energy, therefore the MH ratio be equal to the FBR.
          !!! Note that we also could use the energy as if we would do the job as
          !!! proposal distribution, but we would have to calculate the backward
          !!! energy as well.

          NULLIFY(labels_3d)

          sbpitr => mesh%subpatch%begin()
          ipatch=1
          DO WHILE (ASSOCIATED(sbpitr))
             NM => sbpitr%nnodes

             CALL sbpitr%get_field(labels,labels_3d,info)
             or_fail("Failed to get field labels_3d data.",ppm_error=ppm_error_fatal)

             index_=MCMCGetIndexFromEdgeDensity(labels_3d,Nm)
             !!! get index from Edge Density

             tmpParticle=MCMCFloatingParticles%search(index_(1),index_(2),index_(3))
             !!! search for the particle at this index

             IF (tmpParticle%candlabel.EQ.-1) THEN
                !!! Particle does not exist
                IF (Growth) THEN
                   !!! insert particle

                   !!! the target-distribution probability ratio
                   Child=SaruGetRealVariate().GT.half

                   tmpParticle%candlabel=MERGE(LabelIn,0,Child)
                   !!! TO CHECK
                   !!! TO DO
                   !!! This one does not exist in the original code (to compute the proposal)
                   tmpParticle%proposal=MCMCproposal(index_,labels_3d)

                   MCMCOffBoundarySample_3d= &
                   &  MCMCInsertFloatingParticle(ipatch,index_(1),index_(2),index_(3),tmpParticle,.FALSE.,.FALSE.)
                   GOTO 9999
                ELSE
                   !!! reject
                   MCMCOffBoundarySample_3d=.FALSE.
                   GOTO 9999
                ENDIF
             ELSE
                !!! Particle exists
                IF (Growth) THEN
                   !!! reject
                   MCMCOffBoundarySample_3d=.FALSE.
                   GOTO 9999
                ELSE
                   !!! .NOT.Growth.AND.Particle exists
                   MCMCOffBoundarySample_3d= &
                   & MCMCEraseFloatingParticle(ipatch,index_(1),index_(2),index_(3),tmpParticle,tmpParticle,.FALSE.)
                   GOTO 9999
                ENDIF
             ENDIF

             sbpitr => mesh%subpatch%next()
             ipatch=ipatch+1
          ENDDO

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCOffBoundarySample_3d

!           ! IF (.NOT.Growth) THEN
!           !    !!! we check if we want to delete a floating particle.
!           !    IF (NumberOfFloatingParticles == 0) THEN
!           !       MCMCOffBoundarySample=.FALSE.
!           !       RETURN
!           !    ENDIF
!           !    ParticleIndex = SaruGetIntegerVariate(NumberOfFloatingParticles)
!           !    Particle = MCMCFloatingParticles(ParticleIndex)
!           !    EdgeWeight = m_EdgeImage->GetPixel(vParticle.m_Index)
!           !    !!! the backward probability:
!           !    Pinsertion = EdgeWeight / (MCMCZe + EdgeWeight)
!           !    !!! the forward probability:
!           !    Pdeletion = one / NumberOfFloatingParticles
!           ! ELSE
!           !    !!! sample a location
!           !    Particle.m_Index = MCMCGetIndexFromEdgeDensity()
!           !    !!! the forward probability:
!           !    EdgeWeight = m_EdgeImage->GetPixel(vParticle.m_Index)
!           !    Pinsertion = EdgeWeight / MCMCZe
!           !    !!! the backward probability:
!           !    Pdeletion = one/(vNumberOfFloatingParticles+1)
!           ! ENDIF

        FUNCTION MCMCUpdateAllParticlesFwdProposals(ActiveCandidates) RESULT(info)
          !!!
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_sub_failed

          USE ppm_rc_module_util, ONLY : ppm_rc_MCMCParticlehtable
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          TYPE(ppm_rc_MCMCParticlehtable), INTENT(IN   ) :: ActiveCandidates

          INTEGER                                        :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          REAL(ppm_kind_double) :: t0

          CHARACTER(LEN=ppm_char) :: caller='MCMCUpdateAllParticlesFwdProposals'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          MCMCActiveCandidatesSize=ActiveCandidates%size()

          info=ActiveCandidates%packproposal(MCMCAllParticlesFwdProposals,MCMCActiveCandidatesSize)
          or_fail("ActiveCandidates%packproposal",ppm_error=ppm_error_fatal)

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCUpdateAllParticlesFwdProposals

        LOGICAL FUNCTION MCMCIsRegularParticle_2d (coord1,coord2,MCMCParticle_,LabelIn)
          USE ppm_module_data, ONLY : ppm_error_fatal
          USE ppm_module_error, ONLY : ppm_error,ppm_err_argument
          USE ppm_module_inl_hash, ONLY : htable_null

          USE ppm_rc_module_global, ONLY : htable
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          INTEGER,            INTENT(IN   ) :: LabelIn
          !!! Absolute value label from regionsLabel
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER :: label_region
          INTEGER :: info

          CHARACTER(LEN=*), PARAMETER :: caller='MCMCIsRegularParticle'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          IF (MCMCParticle_%candlabel.EQ.0) THEN
             label_region=htable%search(LabelIn)
             IF (label_region.EQ.htable_null) THEN
                fail("REGION does not exist!",ppm_error=ppm_error_fatal,exit_point=no)
             ENDIF
             MCMCIsRegularParticle_2d=MCMCparents(label_region)%containselement(coord1,coord2,MCMCParticle_)
          ELSE
             label_region=htable%search(MCMCParticle_%candlabel)
             IF (label_region.EQ.htable_null) THEN
                fail("REGION does not exist!",ppm_error=ppm_error_fatal,exit_point=no)
             ENDIF
             MCMCIsRegularParticle_2d=MCMCchildren(label_region)%containselement(coord1,coord2,MCMCParticle_)
          ENDIF
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCIsRegularParticle_2d

        LOGICAL FUNCTION MCMCIsRegularParticle_3d (coord1,coord2,coord3,MCMCParticle_,LabelIn)
          USE ppm_module_data, ONLY : ppm_error_fatal
          USE ppm_module_error, ONLY : ppm_error,ppm_err_argument
          USE ppm_module_inl_hash, ONLY : htable_null

          USE ppm_rc_module_global, ONLY : htable
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          INTEGER,            INTENT(IN   ) :: coord3
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          INTEGER,            INTENT(IN   ) :: LabelIn
          !!! Absolute value label from regionsLabel
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER :: label_region
          INTEGER :: info

          CHARACTER(LEN=*), PARAMETER :: caller='MCMCIsRegularParticle'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          IF (MCMCParticle_%candlabel.EQ.0) THEN
             label_region=htable%search(LabelIn)
             IF (label_region.EQ.htable_null) THEN
                fail("REGION does not exist!",ppm_error=ppm_error_fatal,exit_point=no)
             ENDIF
             MCMCIsRegularParticle_3d=MCMCparents(label_region)%containselement(coord1,coord2,coord3,MCMCParticle_)
          ELSE
             label_region=htable%search(MCMCParticle_%candlabel)
             IF (label_region.EQ.htable_null) THEN
                fail("REGION does not exist!",ppm_error=ppm_error_fatal,exit_point=no)
             ENDIF
             MCMCIsRegularParticle_3d=MCMCchildren(label_region)%containselement(coord1,coord2,coord3,MCMCParticle_)
          ENDIF
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCIsRegularParticle_3d

        FUNCTION MCMCGetProposalNormalizer(LabelIn,CandidateLabel)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_error_fatal
          USE ppm_module_error, ONLY : ppm_error,ppm_err_argument
          USE ppm_module_inl_hash, ONLY : htable_null

          USE ppm_rc_module_global, ONLY : htable
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, INTENT(IN   ) :: LabelIn
          INTEGER, INTENT(IN   ) :: CandidateLabel

          REAL(MK)               :: MCMCGetProposalNormalizer
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER :: label_region
          INTEGER :: info

          CHARACTER(LEN=*), PARAMETER :: caller='MCMCGetProposalNormalizer'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          IF (CandidateLabel.EQ.0) THEN
             label_region=htable%search(LabelIn)
             IF (label_region.EQ.htable_null) THEN
                fail("REGION does not exist!",ppm_error=ppm_error_fatal,exit_point=no)
             ENDIF
             MCMCGetProposalNormalizer=MCMCparentsProposalNormalizer(label_region)
          ELSE
             label_region=htable%search(CandidateLabel)
             IF (label_region.EQ.htable_null) THEN
                fail("REGION does not exist!",ppm_error=ppm_error_fatal,exit_point=no)
             ENDIF
             MCMCGetProposalNormalizer=MCMCchildrenProposalNormalizer(label_region)
          ENDIF
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCGetProposalNormalizer

        FUNCTION MCMCFindParticleA(ParticleAIsFloating,ppmrcdim,rejected) RESULT(info)
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_sub_failed, &
          &   ppm_err_argument
          USE ppm_module_mesh_typedef, ONLY : ppm_t_subpatch_

          USE ppm_rc_module_global, ONLY : one,mesh,labels,MCMCstepsize, &
          &   MCMCuseBiasedProposal
          USE ppm_rc_module_rnd, ONLY : GetPartDistrIndex,SaruGetIntegerVariate,  &
          &   DestroyParticlesDiscrDistr,GenerateParticlesFwdProposalsDiscrDistr
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          LOGICAL, INTENT(IN   ) :: ParticleAIsFloating
          !!!
          INTEGER, INTENT(IN   ) :: ppmrcdim
          !!! Problem dimension
          LOGICAL, INTENT(  OUT) :: rejected

          INTEGER                :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          CLASS(ppm_t_subpatch_), POINTER :: sbpitr

          REAL(ppm_kind_double) :: t0

          INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: labels_2d
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: labels_3d
          INTEGER,             DIMENSION(ppmrcdim)       :: index_
          INTEGER                                        :: ipatch
          INTEGER                                        :: PC
          INTEGER                                        :: ParticleIndex
          INTEGER                                        :: cLabel

          LOGICAL :: DelParticleA

          CHARACTER(LEN=ppm_char) :: caller='MCMCFindParticleA'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)
          !!! For the insertion and deletion of floating particles, we use a
          !!! constant energy, therefore the MH ratio be equal to the FBR.
          !!! Note that we also could use the energy as if we would do the job as
          !!! proposal distribution, but we would have to calculate the backward
          !!! energy as well.

          rejected=.FALSE.

          IF (MCMCuseBiasedProposal.AND..NOT.ParticleAIsFloating) THEN
             info=MCMCUpdateAllParticlesFwdProposals(MCMCActiveCandidates)
             or_fail("MCMCUpdateAllParticlesFwdProposals",ppm_error=ppm_error_fatal)

             !!! Create a discrete distribution over particles
             info=GenerateParticlesFwdProposalsDiscrDistr(MCMCAllParticlesFwdProposals,MCMCActiveCandidatesSize)
             or_fail("GenerateParticlesFwdProposalsDiscrDistr",ppm_error=ppm_error_fatal)
          ELSE
             MCMCActiveCandidatesSize=MCMCActiveCandidates%size()
          ENDIF

          SELECT CASE (ppmrcdim)
          CASE (2)
             NULLIFY(labels_2d)
             sbpitr => mesh%subpatch%begin()
             ipatch=1
             DO WHILE (ASSOCIATED(sbpitr))
                CALL sbpitr%get_field(labels,labels_2d,info)
                or_fail("Failed to get field labels_2d data.",ppm_error=ppm_error_fatal)

                !!! Find particle A:
                DO PC=1,MCMCstepsize
                   IF (MCMCuseBiasedProposal.AND..NOT.ParticleAIsFloating) THEN
                      ParticleIndex=GetPartDistrIndex()
                      MCMC_q_A(PC)=MCMCAllParticlesFwdProposals(ParticleIndex)
                   ELSE
                      ParticleIndex=SaruGetIntegerVariate(MCMCActiveCandidatesSize)
                      MCMC_q_A(PC)=one
                   ENDIF

                   MCMC_CandidateMove(PC)=MCMCActiveCandidates%elementAt(ParticleIndex,index_)
                   IF (MCMC_CandidateMove(PC)%candlabel.EQ.-1) THEN
                      fail("There is no particle in this position!!!",ppm_error=ppm_error_fatal)
                   ENDIF
                   MCMC_CandidateMove_Index(:,PC)=index_
                   cLabel=ABS(labels_2d(index_(1),index_(2)))

                   IF (ParticleAIsFloating) THEN
                      !!! A floating particle was selected.
                      MCMC_Particle_A_IsFloating(PC)=.TRUE.

                      !!! Immediately accept the move (self transition) if the label
                      !!! at the particular position is the same. Since this is a
                      !!! self transition we keep the particle in the proposals.
                      IF (cLabel.EQ.MCMC_CandidateMove(PC)%candlabel) THEN
                         !!! reject
                         rejected=.TRUE.
                         GOTO 8888
                      ENDIF
                      !!! METHOD 5: do it (the regular particle will be deleted) --> 49.06 % / avg. 430 fp

                      !!! METHOD 4: reject (because the backward prob is 0) --> 49.2 %
                      IF (MCMCIsRegularParticle(index_(1),index_(2),MCMC_CandidateMove(PC),cLabel)) THEN
                         !!! reject
                         rejected=.TRUE.
                         GOTO 8888
                      ENDIF

                      !!! if we are here, do the job
                      DelParticleA=MCMCEraseFloatingParticle(ipatch,index_(1),index_(2), &
                      & MCMC_CandidateMove(PC),MCMC_CandidateMove(PC),.TRUE.)

                      !!! METHOD 3 = METHOD 4

                      !!! METHOD 2: apply the floating particle if there is no regular
                      !!! one. If there is a regular one, we apply this and keep the
                      !!! float. (convert the float to regular), i.e. we just don't delete
                      !!! the floating p. --> gets bi-modal with mean 48.75%, avg f.p. 8250

                      !!! METHOD 1: just do it (remove both particles)
                      !!! convert it to a real particle (floating particles always store
                      !!! the region label they belong to). We need to check if we the
                      !!! need to set the candidate label to 0.
                   ENDIF !ParticleAIsFloating

                   !!! Fill in some properties for this particles
                   MCMC_LabelsBeforeJump_A(PC)=cLabel

                   IF (MCMC_Particle_A_IsFloating(PC)) THEN
                      MCMC_q_A(PC)=MCMC_q_A(PC)/REAL(MCMCFloatingParticlesProposalNormalizer,MK)
                   ELSE
                      MCMC_q_A(PC)=MCMC_q_A(PC)/MCMCGetProposalNormalizer(cLabel,MCMC_CandidateMove(PC)%candlabel)
                   ENDIF
                ENDDO !PC=1,MCMCstepsize

                sbpitr => mesh%subpatch%next()
                ipatch=ipatch+1
             ENDDO !ASSOCIATED(sbpitr)
             NULLIFY(labels_2d)
          CASE (3)
             NULLIFY(labels_3d)
             sbpitr => mesh%subpatch%begin()
             ipatch=1
             DO WHILE (ASSOCIATED(sbpitr))
                CALL sbpitr%get_field(labels,labels_3d,info)
                or_fail("Failed to get field labels_3d data.",ppm_error=ppm_error_fatal)

                !!! Find particle A:
                DO PC=1,MCMCstepsize
                   IF (MCMCuseBiasedProposal.AND..NOT.ParticleAIsFloating) THEN
                      ParticleIndex=GetPartDistrIndex()
                      MCMC_q_A(PC)=MCMCAllParticlesFwdProposals(ParticleIndex)
                   ELSE
                      ParticleIndex=SaruGetIntegerVariate(MCMCActiveCandidatesSize)
                      MCMC_q_A(PC)=one
                   ENDIF

                   MCMC_CandidateMove(PC)=MCMCActiveCandidates%elementAt(ParticleIndex,index_)
                   IF (MCMC_CandidateMove(PC)%candlabel.EQ.-1) THEN
                      fail("There is no particle in this position!!!",ppm_error=ppm_error_fatal)
                   ENDIF
                   MCMC_CandidateMove_Index(:,PC)=index_
                   cLabel=ABS(labels_3d(index_(1),index_(2),index_(3)))

                   IF (ParticleAIsFloating) THEN
                      !!! A floating particle was selected.
                      MCMC_Particle_A_IsFloating(PC)=.TRUE.

                      !!! Immediately accept the move (self transition) if the label
                      !!! at the particular position is the same. Since this is a
                      !!! self transition we keep the particle in the proposals.
                      IF (cLabel.EQ.MCMC_CandidateMove(PC)%candlabel) THEN
                         !!! reject
                         rejected=.TRUE.
                         GOTO 8888
                      ENDIF
                      !!! METHOD 5: do it (the regular particle will be deleted) --> 49.06 % / avg. 430 fp

                      !!! METHOD 4: reject (because the backward prob is 0) --> 49.2 %
                      IF (MCMCIsRegularParticle(index_(1),index_(2),index_(3),MCMC_CandidateMove(PC),cLabel)) THEN
                         !!! reject
                         rejected=.TRUE.
                         GOTO 8888
                      ENDIF

                      !!! if we are here, do the job
                      DelParticleA=MCMCEraseFloatingParticle(ipatch,index_(1),index_(2),index_(3), &
                      & MCMC_CandidateMove(PC),MCMC_CandidateMove(PC),.TRUE.)

                      !!! METHOD 3 = METHOD 4

                      !!! METHOD 2: apply the floating particle if there is no regular
                      !!! one. If there is a regular one, we apply this and keep the
                      !!! float. (convert the float to regular), i.e. we just don't delete
                      !!! the floating p. --> gets bi-modal with mean 48.75%, avg f.p. 8250

                      !!! METHOD 1: just do it (remove both particles)
                      !!! convert it to a real particle (floating particles always store
                      !!! the region label they belong to). We need to check if we the
                      !!! need to set the candidate label to 0.
                   ENDIF !ParticleAIsFloating

                   !!! Fill in some properties for this particles
                   MCMC_LabelsBeforeJump_A(PC)=cLabel

                   IF (MCMC_Particle_A_IsFloating(PC)) THEN
                      MCMC_q_A(PC)=MCMC_q_A(PC)/REAL(MCMCFloatingParticlesProposalNormalizer,MK)
                   ELSE
                      MCMC_q_A(PC)=MCMC_q_A(PC)/MCMCGetProposalNormalizer(cLabel,MCMC_CandidateMove(PC)%candlabel)
                   ENDIF
                ENDDO !PC=1,MCMCstepsize

                sbpitr => mesh%subpatch%next()
                ipatch=ipatch+1
             ENDDO
             NULLIFY(labels_3d)
          END SELECT

        8888 CONTINUE
          IF (MCMCuseBiasedProposal.AND..NOT.ParticleAIsFloating) THEN
             !!! make sure that there is no discrete distribution available in memory
             info=DestroyParticlesDiscrDistr()
             or_fail("DestroyParticlesDiscrDistr")
          ENDIF
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCFindParticleA

        FUNCTION MCMCAddAndRemoveParticlesWhenMove_2d(ipatch,coord,Nm,labels_,MCMCParticle_) RESULT(info)
          !!! Updates the DS (label image and the regular particle containers) when
          !!! applying a particle. Note that floating particles will not be updated!
          !!! The method only ensures that L and the regular  particles are correct.
          !!! The method expects the label image NOT to be updated already. Else the
          !!! operations performed will be wrong.
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error,ppm_err_sub_failed

          USE ppm_rc_module_global, ONLY : FORBIDDEN
          USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType, &
          &   BG_ConnectivityType,IsEnclosedByLabel_BGConnectivity
          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list,MCMCLabelImageHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,                             INTENT(IN   ) :: ipatch
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: coord
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: Nm
          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
          TYPE(MCMCParticle),                  INTENT(IN   ) :: MCMCParticle_
          INTEGER                                            :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle
          TYPE(MCMCParticle) :: ReverseParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0
          REAL(MK)              :: proposal

          INTEGER, DIMENSION(:,:), POINTER :: tmplabels_
          INTEGER, DIMENSION(:,:), POINTER :: tmplabels

          INTEGER, DIMENSION(2)            :: ll,ld,dd
          INTEGER                          :: LabelFrom
          INTEGER                          :: LabelTo
          INTEGER                          :: StoreLabel
          INTEGER                          :: LabelNgh
          INTEGER                          :: i,j

          LOGICAL :: ReverseParticleIsFloating
          LOGICAL :: Replaced
          LOGICAL :: HasOtherMother

          CHARACTER(LEN=ppm_char) :: caller='MCMCAddAndRemoveParticlesWhenMove'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          !!! Initilize an neighborhooditerator for fast access to the label image
          !!! in the vicinity of the particle of interest.
          tmplabels_ => labels_(coord(1)-2:coord(1)+2,coord(2)-2:coord(2)+2)

          StoreLabel=tmplabels_(3,3)

          !!! Initialize locals from the particle
          LabelFrom=ABS(StoreLabel)
          LabelTo=MCMCParticle_%candlabel

          !!! We remove the particle and insert the reverse particle:
          !!! Here we cannot decide if the reverse particle is indeed
          !!! a floating particle (this is because we apply B first, then A).
          !!! The solution is that this method only works on the regular particle
          !!! set. Floating particles will be detected and treated outside of
          !!! this method. Here we omit the potential insertion of floating
          !!! particles.
          !!! In order to (maybe) replace the particle with its reverse particle we:
          !!! Simulate the move, calculate the proposal for the backward particle,
          !!! create the backward particle, check if the backward particle is
          !!! floating, and finally restore the label image:

          tmplabels_(3,3)=-LabelTo

          tmplabels => tmplabels_(2:4,2:4)
          proposal=MCMCproposal(tmplabels)

          ReverseParticle=MCMCParticle(LabelFrom,proposal)
          !                                                        (tmplabels,LabelIn,CandidateLabel)
          !                                                         tmplabels,       ,ReverseParticle%candlabel

          tmplabels_(3,3)=StoreLabel

          !!! Insert the reverse particle (if its not floating, see above):
          IF (.NOT.ReverseParticleIsFloating) THEN
             info=MCMCInsertCandidatesToContainers(ipatch,coord(1),coord(2),ReverseParticle,LabelTo,.TRUE.,Replaced)
             or_fail("MCMCInsertCandidatesToContainers")
          ENDIF

          !!! erase the currently applied particle (if its floating this will not hurt
          !!! to call erase here).
          tmpParticle=MCMCParticle(LabelTo,zero)

          info=MCMCEraseCandidatesFromContainers(ipatch,coord(1),coord(2),tmpParticle,LabelFrom,.TRUE.)
          or_fail("MCMCEraseCandidatesFromContainers")

          !!! Since we are storing the parents separately for each region, the
          !!! following patch is needed. We need to shift the mother particle into
          !!! the parents container of each others region:
          IF (LabelFrom.NE.0.AND.LabelTo.NE.0) THEN
             tmpParticle=MCMCParticle(0,proposal)

             info=MCMCEraseCandidatesFromContainers(ipatch,coord(1),coord(2),tmpParticle,LabelFrom,.TRUE.)
             or_fail("MCMCEraseCandidatesFromContainers")

             info=MCMCInsertCandidatesToContainers(ipatch,coord(1),coord(2),tmpParticle,LabelTo,.TRUE.,Replaced)
             or_fail("MCMCInsertCandidatesToContainers")
          ENDIF

          !!! What particle would be added or removed to the contour lists
          !!! of the currently shrinking region:
          !!! (compare to also AddNeighborsAtRemove(..))
          IF (LabelFrom.NE.0) THEN
             !!! a FG region is shrinking
             DO i=1,BG_ConnectivityType%NumberOfNeighbors
                ld=coord+BG_ConnectivityType%NeighborsPoints(:,i)
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=BG_ConnectivityType%NeighborsPoints(:,i)+3
                LabelNgh=tmplabels_(ll(1),ll(2))
                !!! Check if new points enter the contour (internal point becomes a parent):
                !!! Internal points of this region are positive:
                IF (LabelNgh.GT.0.AND.LabelNgh.EQ.LabelFrom) THEN
                   tmplabels => tmplabels_(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1)
                   proposal=MCMCproposal(tmplabels)

                   tmpParticle=MCMCParticle(0,proposal)

                   info=MCMCInsertCandidatesToContainers(ipatch,ld(1),ld(2),tmpParticle,LabelFrom,.TRUE.,Replaced)
                   or_fail("MCMCInsertCandidatesToContainers")

                   IF (.NOT.Replaced) THEN
                      tmplabels_(ll(1),ll(2))=-LabelFrom

                      ALLOCATE(seed,STAT=info)
                      or_fail_alloc("seed")
                      CALL seed%add(ld(1),ld(2),LabelFrom)
                      CALL MCMCLabelImageHistory(ipatch)%push(seed,info)
                      or_fail("MCMCLabelImageHistory(ipatch)%push")
                   ENDIF
                ENDIF !LabelNgh.GT.0.AND.LabelNgh.EQ.LabelFrom
             ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors

             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ld=coord+FG_ConnectivityType%NeighborsPoints(:,i)
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+3
                LabelNgh=ABS(tmplabels_(ll(1),ll(2)))

                !!! check if there are FG-neighbors with no other mother of the
                !!! same label --> orphin.
                !!! we first 'simulate' the move:

                StoreLabel=tmplabels_(3,3)

                tmplabels_(3,3)=-LabelTo

                IF (LabelNgh.NE.LabelFrom) THEN
                   !!! check if this neighbor has other mothers from this label:
                   HasOtherMother=.FALSE.
                   DO j=1,FG_ConnectivityType%NumberOfNeighbors
                      dd=ll+FG_ConnectivityType%NeighborsPoints(:,i)
                      IF (ABS(tmplabels_(dd(1),dd(2))).EQ.LabelFrom) THEN
                         HasOtherMother=.TRUE.
                         EXIT
                      ENDIF
                   ENDDO
                   IF (.NOT.HasOtherMother) THEN
                      !!! The orphin has label equal to what we read from the
                      !!! label image and has a candidate label of the
                      !!! currently shrinking region.
                      tmplabels => tmplabels_(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1)
                      proposal=MCMCproposal(tmplabels)

                      tmpParticle=MCMCParticle(LabelFrom,proposal)

                      info=MCMCEraseCandidatesFromContainers(ipatch,ld(1),ld(2),tmpParticle,LabelNgh,.TRUE.)
                      or_fail("MCMCEraseCandidatesFromContainers")
                   ENDIF !.NOT.HasOtherMother
                ENDIF !LabelNgh.NE.LabelFrom

                tmplabels_(3,3)=StoreLabel
             ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors
          ENDIF ! LabelFrom.NE.0


          !!! Growing: figure out the changes of candidates for the expanding region
          IF (LabelTo.NE.0) THEN
             !!! we are growing

             !!! Neighbors: Figure out what (neighboring)mother points are going
             !!! to be interior points:

             !!! simulate the move
             StoreLabel=tmplabels_(3,3)

             tmplabels_(3,3)=-LabelTo

             DO i=1,BG_ConnectivityType%NumberOfNeighbors
                ld=coord+BG_ConnectivityType%NeighborsPoints(:,i)
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=BG_ConnectivityType%NeighborsPoints(:,i)+3
                LabelNgh=tmplabels_(ll(1),ll(2))

                tmplabels => tmplabels_(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1)

                !!! TODO: the isEnclosedByLabel method could use the above iterator
                IF (LabelNgh.EQ.-LabelTo.AND.IsEnclosedByLabel_BGConnectivity(tmplabels,LabelTo)) THEN
                   !!! Remove the parent that got enclosed; it had a the label
                   !!! of the currently expanding region and a candidate label
                   !!! of 0.
                   tmpParticle=MCMCParticle(0,zero)

                   info=MCMCEraseCandidatesFromContainers(ipatch,ld(1),ld(2),tmpParticle,LabelTo,.TRUE.)
                   or_fail("MCMCEraseCandidatesFromContainers")

                   !!! update the label image (we're not using the update
                   !!! mechanism of the optimizer anymore):
                   tmplabels_(ll(1),ll(2))=LabelTo

                   ALLOCATE(seed,STAT=info)
                   or_fail_alloc("seed")
                   CALL seed%add(ld(1),ld(2),LabelNgh)
                   CALL MCMCLabelImageHistory(ipatch)%push(seed,info)
                   or_fail("MCMCLabelImageHistory(ipatch)%push")
                ENDIF
             ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors

             tmplabels_(3,3)=StoreLabel

             !!! Figure out if a point renders to a candidate. These are
             !!! all the FG-neighbors with a different label that are not yet
             !!! candidates of the currently expanding region.
             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ld=coord+FG_ConnectivityType%NeighborsPoints(:,i)
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+3
                LabelNgh=tmplabels_(ll(1),ll(2))

                IF (ABS(LabelNgh).NE.LabelTo.AND.LabelNgh.NE.FORBIDDEN) THEN
                   !!! check if there is no other mother (hence the particle is
                   !!! not in the container yet). This we could do by checking the
                   !!! neighborhood of the label image or by checking the
                   !!! cointainers.
                   !!! Here: we check the (not yet updated!) label image.
                   HasOtherMother=.FALSE.
                   DO j=1,FG_ConnectivityType%NumberOfNeighbors
                      dd=ll+FG_ConnectivityType%NeighborsPoints(:,i)
                      IF (ABS(tmplabels_(dd(1),dd(2))).EQ.LabelTo) THEN
                         HasOtherMother=.TRUE.
                         EXIT
                      ENDIF
                   ENDDO
                   IF (.NOT.HasOtherMother) THEN
                      !!! This is a new child. It's current label we have to read
                      !!! from the label image, the candidate label is the label of
                      !!! the currently expanding region.
                      tmplabels => tmplabels_(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1)
                      proposal=MCMCproposal(tmplabels)

                      tmpParticle=MCMCParticle(LabelTo,proposal)

                      info=MCMCInsertCandidatesToContainers(ipatch,ld(1),ld(2),tmpParticle,ABS(LabelNgh),.TRUE.,Replaced)
                      or_fail("MCMCInsertCandidatesToContainers")

                      IF (.NOT.Replaced) THEN
                         tmplabels_(ll(1),ll(2))=SIGN(LabelNgh,-1)

                         ALLOCATE(seed,STAT=info)
                         or_fail_alloc("seed")
                         CALL seed%add(ld(1),ld(2),LabelNgh)
                         CALL MCMCLabelImageHistory(ipatch)%push(seed,info)
                         or_fail("MCMCLabelImageHistory(ipatch)%push")
                      ENDIF
                   ENDIF !.NOT.HasOtherMother
                ENDIF !ABS(LabelNgh).NE.LabelTo.AND.LabelNgh.NE.FORBIDDEN
             ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors
          ENDIF ! LabelTo.NE.0
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCAddAndRemoveParticlesWhenMove_2d

        FUNCTION MCMCAddAndRemoveParticlesWhenMove_3d(ipatch,coord,Nm,labels_,MCMCParticle_) RESULT(info)
          !!! Updates the DS (label image and the regular particle containers) when
          !!! applying a particle. Note that floating particles will not be updated!
          !!! The method only ensures that L and the regular  particles are correct.
          !!! The method expects the label image NOT to be updated already. Else the
          !!! operations performed will be wrong.
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error,ppm_err_sub_failed

          USE ppm_rc_module_global, ONLY : FORBIDDEN
          USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType, &
          &   BG_ConnectivityType,IsEnclosedByLabel_BGConnectivity
          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list,MCMCLabelImageHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,                               INTENT(IN   ) :: ipatch
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: coord
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: Nm
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
          TYPE(MCMCParticle),                    INTENT(IN   ) :: MCMCParticle_
          INTEGER                                              :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle
          TYPE(MCMCParticle) :: ReverseParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0
          REAL(MK)              :: proposal

          INTEGER, DIMENSION(:,:,:), POINTER :: tmplabels_
          INTEGER, DIMENSION(:,:,:), POINTER :: tmplabels

          INTEGER, DIMENSION(3)              :: ll,ld,dd
          INTEGER                            :: LabelFrom
          INTEGER                            :: LabelTo
          INTEGER                            :: StoreLabel
          INTEGER                            :: LabelNgh
          INTEGER                            :: i,j

          LOGICAL :: ReverseParticleIsFloating
          LOGICAL :: Replaced
          LOGICAL :: HasOtherMother

          CHARACTER(LEN=ppm_char) :: caller='MCMCAddAndRemoveParticlesWhenMove'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          !!! Initilize an neighborhooditerator for fast access to the label image
          !!! in the vicinity of the particle of interest.
          tmplabels_ => labels_(coord(1)-2:coord(1)+2,coord(2)-2:coord(2)+2,coord(3)-2:coord(3)+2)

          StoreLabel=tmplabels_(3,3,3)

          !!! Initialize locals from the particle
          LabelFrom=ABS(StoreLabel)
          LabelTo=MCMCParticle_%candlabel

          !!! We remove the particle and insert the reverse particle:
          !!! Here we cannot decide if the reverse particle is indeed
          !!! a floating particle (this is because we apply B first, then A).
          !!! The solution is that this method only works on the regular particle
          !!! set. Floating particles will be detected and treated outside of
          !!! this method. Here we omit the potential insertion of floating
          !!! particles.
          !!! In order to (maybe) replace the particle with its reverse particle we:
          !!! Simulate the move, calculate the proposal for the backward particle,
          !!! create the backward particle, check if the backward particle is
          !!! floating, and finally restore the label image:

          tmplabels_(3,3,3)=-LabelTo

          tmplabels => tmplabels_(2:4,2:4,2:4)
          proposal=MCMCproposal(tmplabels)

          ReverseParticle=MCMCParticle(LabelFrom,proposal)
          !                                                        (tmplabels,LabelIn,CandidateLabel)
          !                                                         tmplabels,       ,ReverseParticle%candlabel

          tmplabels_(3,3,3)=StoreLabel

          !!! Insert the reverse particle (if its not floating, see above):
          IF (.NOT.ReverseParticleIsFloating) THEN
             info=MCMCInsertCandidatesToContainers(ipatch,coord(1),coord(2),coord(3), &
             &    ReverseParticle,LabelTo,.TRUE.,Replaced)
             or_fail("MCMCInsertCandidatesToContainers")
          ENDIF

          !!! erase the currently applied particle (if its floating this will not hurt
          !!! to call erase here).
          tmpParticle=MCMCParticle(LabelTo,zero)

          info=MCMCEraseCandidatesFromContainers(ipatch,coord(1),coord(2),coord(3), &
          &    tmpParticle,LabelFrom,.TRUE.)
          or_fail("MCMCEraseCandidatesFromContainers")

          !!! Since we are storing the parents separately for each region, the
          !!! following patch is needed. We need to shift the mother particle into
          !!! the parents container of each others region:
          IF (LabelFrom.NE.0.AND.LabelTo.NE.0) THEN
             tmpParticle=MCMCParticle(0,proposal)

             info=MCMCEraseCandidatesFromContainers(ipatch,coord(1),coord(2),coord(3), &
             &    tmpParticle,LabelFrom,.TRUE.)
             or_fail("MCMCEraseCandidatesFromContainers")

             info=MCMCInsertCandidatesToContainers(ipatch,coord(1),coord(2),coord(3), &
             &    tmpParticle,LabelTo,.TRUE.,Replaced)
             or_fail("MCMCInsertCandidatesToContainers")
          ENDIF

          !!! What particle would be added or removed to the contour lists
          !!! of the currently shrinking region:
          !!! (compare to also AddNeighborsAtRemove(..))
          IF (LabelFrom.NE.0) THEN
             !!! a FG region is shrinking
             DO i=1,BG_ConnectivityType%NumberOfNeighbors
                ld=coord+BG_ConnectivityType%NeighborsPoints(:,i)
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=BG_ConnectivityType%NeighborsPoints(:,i)+3
                LabelNgh=tmplabels_(ll(1),ll(2),ll(3))
                !!! Check if new points enter the contour (internal point becomes a parent):
                !!! Internal points of this region are positive:
                IF (LabelNgh.GT.0.AND.LabelNgh.EQ.LabelFrom) THEN
                   tmplabels => tmplabels_(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1,ll(3)-1:ll(3)+1)
                   proposal=MCMCproposal(tmplabels)

                   tmpParticle=MCMCParticle(0,proposal)

                   info=MCMCInsertCandidatesToContainers(ipatch,ld(1),ld(2),ld(3), &
                   &    tmpParticle,LabelFrom,.TRUE.,Replaced)
                   or_fail("MCMCInsertCandidatesToContainers")

                   IF (.NOT.Replaced) THEN
                      tmplabels_(ll(1),ll(2),ll(3))=-LabelFrom

                      ALLOCATE(seed,STAT=info)
                      or_fail_alloc("seed")
                      CALL seed%add(ld(1),ld(2),ld(3),LabelFrom)
                      CALL MCMCLabelImageHistory(ipatch)%push(seed,info)
                      or_fail("MCMCLabelImageHistory(ipatch)%push")
                   ENDIF
                ENDIF !LabelNgh.GT.0.AND.LabelNgh.EQ.LabelFrom
             ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors

             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ld=coord+FG_ConnectivityType%NeighborsPoints(:,i)
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+3
                LabelNgh=ABS(tmplabels_(ll(1),ll(2),ll(3)))

                !!! check if there are FG-neighbors with no other mother of the
                !!! same label --> orphin.
                !!! we first 'simulate' the move:

                StoreLabel=tmplabels_(3,3,3)

                tmplabels_(3,3,3)=-LabelTo

                IF (LabelNgh.NE.LabelFrom) THEN
                   !!! check if this neighbor has other mothers from this label:
                   HasOtherMother=.FALSE.
                   DO j=1,FG_ConnectivityType%NumberOfNeighbors
                      dd=ll+FG_ConnectivityType%NeighborsPoints(:,i)
                      IF (ABS(tmplabels_(dd(1),dd(2),dd(3))).EQ.LabelFrom) THEN
                         HasOtherMother=.TRUE.
                         EXIT
                      ENDIF
                   ENDDO
                   IF (.NOT.HasOtherMother) THEN
                      !!! The orphin has label equal to what we read from the
                      !!! label image and has a candidate label of the
                      !!! currently shrinking region.
                      tmplabels => tmplabels_(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1,ll(3)-1:ll(3)+1)
                      proposal=MCMCproposal(tmplabels)

                      tmpParticle=MCMCParticle(LabelFrom,proposal)

                      info=MCMCEraseCandidatesFromContainers(ipatch,ld(1),ld(2),ld(3), &
                      &    tmpParticle,LabelNgh,.TRUE.)
                      or_fail("MCMCEraseCandidatesFromContainers")
                   ENDIF !.NOT.HasOtherMother
                ENDIF !LabelNgh.NE.LabelFrom

                tmplabels_(3,3,3)=StoreLabel
             ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors
          ENDIF ! LabelFrom.NE.0


          !!! Growing: figure out the changes of candidates for the expanding region
          IF (LabelTo.NE.0) THEN
             !!! we are growing

             !!! Neighbors: Figure out what (neighboring)mother points are going
             !!! to be interior points:

             !!! simulate the move
             StoreLabel=tmplabels_(3,3,3)

             tmplabels_(3,3,3)=-LabelTo

             DO i=1,BG_ConnectivityType%NumberOfNeighbors
                ld=coord+BG_ConnectivityType%NeighborsPoints(:,i)
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=BG_ConnectivityType%NeighborsPoints(:,i)+3
                LabelNgh=tmplabels_(ll(1),ll(2),ll(3))

                tmplabels => tmplabels_(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1,ll(3)-1:ll(3)+1)

                !!! TODO: the isEnclosedByLabel method could use the above iterator
                IF (LabelNgh.EQ.-LabelTo.AND.IsEnclosedByLabel_BGConnectivity(tmplabels,LabelTo)) THEN
                   !!! Remove the parent that got enclosed; it had a the label
                   !!! of the currently expanding region and a candidate label
                   !!! of 0.
                   tmpParticle=MCMCParticle(0,zero)

                   info=MCMCEraseCandidatesFromContainers(ipatch,ld(1),ld(2),ld(3), &
                   &    tmpParticle,LabelTo,.TRUE.)
                   or_fail("MCMCEraseCandidatesFromContainers")

                   !!! update the label image (we're not using the update
                   !!! mechanism of the optimizer anymore):
                   tmplabels_(ll(1),ll(2),ll(3))=LabelTo

                   ALLOCATE(seed,STAT=info)
                   or_fail_alloc("seed")
                   CALL seed%add(ld(1),ld(2),ld(3),LabelNgh)
                   CALL MCMCLabelImageHistory(ipatch)%push(seed,info)
                   or_fail("MCMCLabelImageHistory(ipatch)%push")
                ENDIF
             ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors

             tmplabels_(3,3,3)=StoreLabel

             !!! Figure out if a point renders to a candidate. These are
             !!! all the FG-neighbors with a different label that are not yet
             !!! candidates of the currently expanding region.
             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ld=coord+FG_ConnectivityType%NeighborsPoints(:,i)
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+3
                LabelNgh=tmplabels_(ll(1),ll(2),ll(3))

                IF (ABS(LabelNgh).NE.LabelTo.AND.LabelNgh.NE.FORBIDDEN) THEN
                   !!! check if there is no other mother (hence the particle is
                   !!! not in the container yet). This we could do by checking the
                   !!! neighborhood of the label image or by checking the
                   !!! cointainers.
                   !!! Here: we check the (not yet updated!) label image.
                   HasOtherMother=.FALSE.
                   DO j=1,FG_ConnectivityType%NumberOfNeighbors
                      dd=ll+FG_ConnectivityType%NeighborsPoints(:,i)
                      IF (ABS(tmplabels_(dd(1),dd(2),dd(3))).EQ.LabelTo) THEN
                         HasOtherMother=.TRUE.
                         EXIT
                      ENDIF
                   ENDDO
                   IF (.NOT.HasOtherMother) THEN
                      !!! This is a new child. It's current label we have to read
                      !!! from the label image, the candidate label is the label of
                      !!! the currently expanding region.
                      tmplabels => tmplabels_(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1,ll(3)-1:ll(3)+1)
                      proposal=MCMCproposal(tmplabels)

                      tmpParticle=MCMCParticle(LabelTo,proposal)

                      info=MCMCInsertCandidatesToContainers(ipatch,ld(1),ld(2),ld(3), &
                      &    tmpParticle,ABS(LabelNgh),.TRUE.,Replaced)
                      or_fail("MCMCInsertCandidatesToContainers")

                      IF (.NOT.Replaced) THEN
                         tmplabels_(ll(1),ll(2),ll(3))=SIGN(LabelNgh,-1)

                         ALLOCATE(seed,STAT=info)
                         or_fail_alloc("seed")
                         CALL seed%add(ld(1),ld(2),ld(3),LabelNgh)
                         CALL MCMCLabelImageHistory(ipatch)%push(seed,info)
                         or_fail("MCMCLabelImageHistory(ipatch)%push")
                      ENDIF
                   ENDIF !.NOT.HasOtherMother
                ENDIF !ABS(LabelNgh).NE.LabelTo.AND.LabelNgh.NE.FORBIDDEN
             ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors
          ENDIF ! LabelTo.NE.0
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCAddAndRemoveParticlesWhenMove_3d

        FUNCTION MCMCupdateProposalsAndFilterTopologyInNeighborhood_2d(ipatch,coord,Nm,labels_) RESULT(info)
          !!! Proposal update and topology control are combined in this method for
          !!! efficiency reason.
          !!! All particles in the 3x3x3 unit cube to be updated. We use the
          !!! label image to find the particles. Their proposal is updated. If the
          !!! particles do not fulfill the topological constraints, they will be
          !!! removed from the containers.
          !!! Note: floating particles will not be updated!
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_sub_failed

          USE ppm_rc_module_topologicalnumber, ONLY : BG_ConnectivityType
          USE ppm_rc_module_util, ONLY : MCMCParticle
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,                             INTENT(IN   ) :: ipatch
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: coord
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: Nm
          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
          INTEGER                                            :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle), DIMENSION(9) :: MCMCParticles

          REAL(ppm_kind_double) :: t0

          INTEGER, DIMENSION(:,:), POINTER :: tmplabels_
          INTEGER, DIMENSION(:,:), POINTER :: tmplabels
          INTEGER, DIMENSION(2)            :: ll,ld
          INTEGER                          :: ContainerLabel
          INTEGER                          :: cLabel
          INTEGER                          :: nsize
          INTEGER                          :: i,j

          LOGICAL :: Replaced

          CHARACTER(LEN=ppm_char) :: caller='MCMCupdateProposalsAndFilterTopologyInNeighborhood'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          tmplabels_ => labels_(coord(1)-2:coord(1)+2,coord(2)-2:coord(2)+2)
          tmplabels  => tmplabels_(2:4,2:4)

          cLabel=ABS(tmplabels(2,2))

          nsize=MCMCgetRegularParticlesAtIndex(coord(1),coord(2),Nm,tmplabels,MCMCParticles)

          DO i=1,nsize
             ContainerLabel=MERGE(cLabel,MCMCParticles(i)%candlabel,MCMCParticles(i)%candlabel.EQ.0)
             IF (MCMCIsParticleTopoValid(tmplabels,MCMCParticles(i)%candlabel)) THEN
                !!! replace the particle with the new propsal
                info=MCMCInsertCandidatesToContainers(ipatch,coord(1),coord(2), &
                &    MCMCParticles(i),ContainerLabel,.TRUE.,Replaced)
                or_fail("MCMCInsertCandidatesToContainers")
             ELSE
                !!! remove the particle
                info=MCMCEraseCandidatesFromContainers(ipatch,coord(1),coord(2), &
                &    MCMCParticles(i),ContainerLabel,.TRUE.)
                or_fail("MCMCEraseCandidatesFromContainers")
             ENDIF
          ENDDO

          DO i=1,BG_ConnectivityType%NumberOfNeighbors
             ld=coord+BG_ConnectivityType%NeighborsPoints(:,i)
             IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
             ll=BG_ConnectivityType%NeighborsPoints(:,i)+3

             tmplabels => tmplabels_(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1)

             cLabel=ABS(tmplabels(2,2))

             nsize=MCMCgetRegularParticlesAtIndex(ld(1),ld(2),Nm,tmplabels,MCMCParticles)

             DO j=1,nsize
                ContainerLabel=MERGE(cLabel,MCMCParticles(j)%candlabel,MCMCParticles(j)%candlabel.EQ.0)
                IF (MCMCIsParticleTopoValid(tmplabels,MCMCParticles(j)%candlabel)) THEN
                   !!! replace the particle with the new propsal
                   info=MCMCInsertCandidatesToContainers(ipatch,ld(1),ld(2), &
                   &    MCMCParticles(j),ContainerLabel,.TRUE.,Replaced)
                   or_fail("MCMCInsertCandidatesToContainers")
                ELSE
                   !!! remove the particle
                   info=MCMCEraseCandidatesFromContainers(ipatch,ld(1),ld(2), &
                   &    MCMCParticles(j),ContainerLabel,.TRUE.)
                   or_fail("MCMCEraseCandidatesFromContainers")
                ENDIF
             ENDDO !j=1,nsize
          ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCupdateProposalsAndFilterTopologyInNeighborhood_2d

        FUNCTION MCMCupdateProposalsAndFilterTopologyInNeighborhood_3d(ipatch,coord,Nm,labels_) RESULT(info)
          !!! Proposal update and topology control are combined in this method for
          !!! efficiency reason.
          !!! All particles in the 3x3x3 unit cube to be updated. We use the
          !!! label image to find the particles. Their proposal is updated. If the
          !!! particles do not fulfill the topological constraints, they will be
          !!! removed from the containers.
          !!! Note: floating particles will not be updated!
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_sub_failed

          USE ppm_rc_module_topologicalnumber, ONLY : BG_ConnectivityType
          USE ppm_rc_module_util, ONLY : MCMCParticle
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,                               INTENT(IN   ) :: ipatch
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: coord
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: Nm
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
          INTEGER                                              :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle), DIMENSION(27) :: MCMCParticles

          REAL(ppm_kind_double) :: t0

          INTEGER, DIMENSION(:,:,:), POINTER :: tmplabels_
          INTEGER, DIMENSION(:,:,:), POINTER :: tmplabels
          INTEGER, DIMENSION(3)              :: ll,ld
          INTEGER                            :: ContainerLabel
          INTEGER                            :: cLabel
          INTEGER                            :: nsize
          INTEGER                            :: i,j

          LOGICAL :: Replaced

          CHARACTER(LEN=ppm_char) :: caller='MCMCupdateProposalsAndFilterTopologyInNeighborhood'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          tmplabels_ => labels_(coord(1)-2:coord(1)+2,coord(2)-2:coord(2)+2,coord(3)-2:coord(3)+2)
          tmplabels  => tmplabels_(2:4,2:4,2:4)

          cLabel=ABS(tmplabels(2,2,2))

          nsize=MCMCgetRegularParticlesAtIndex(coord(1),coord(2),coord(3),Nm,tmplabels,MCMCParticles)

          DO i=1,nsize
             ContainerLabel=MERGE(cLabel,MCMCParticles(i)%candlabel,MCMCParticles(i)%candlabel.EQ.0)
             IF (MCMCIsParticleTopoValid(tmplabels,MCMCParticles(i)%candlabel)) THEN
                !!! replace the particle with the new propsal
                info=MCMCInsertCandidatesToContainers(ipatch,coord(1),coord(2),coord(3), &
                &    MCMCParticles(i),ContainerLabel,.TRUE.,Replaced)
                or_fail("MCMCInsertCandidatesToContainers")
             ELSE
                !!! remove the particle
                info=MCMCEraseCandidatesFromContainers(ipatch,coord(1),coord(2),coord(3), &
                &    MCMCParticles(i),ContainerLabel,.TRUE.)
                or_fail("MCMCEraseCandidatesFromContainers")
             ENDIF
          ENDDO

          DO i=1,BG_ConnectivityType%NumberOfNeighbors
             ld=coord+BG_ConnectivityType%NeighborsPoints(:,i)
             !!! TODO
             !!! TOCHECK: I think this one should be changed in parallel version
             IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
             ll=BG_ConnectivityType%NeighborsPoints(:,i)+3

             tmplabels => tmplabels_(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1,ll(3)-1:ll(3)+1)

             cLabel=ABS(tmplabels(2,2,2))

             nsize=MCMCgetRegularParticlesAtIndex(ld(1),ld(2),ld(3),Nm,tmplabels,MCMCParticles)

             DO j=1,nsize
                ContainerLabel=MERGE(cLabel,MCMCParticles(j)%candlabel,MCMCParticles(j)%candlabel.EQ.0)
                IF (MCMCIsParticleTopoValid(tmplabels,MCMCParticles(j)%candlabel)) THEN
                   !!! replace the particle with the new propsal
                   info=MCMCInsertCandidatesToContainers(ipatch,ld(1),ld(2),ld(3), &
                   &    MCMCParticles(j),ContainerLabel,.TRUE.,Replaced)
                   or_fail("MCMCInsertCandidatesToContainers")
                ELSE
                   !!! remove the particle
                   info=MCMCEraseCandidatesFromContainers(ipatch,ld(1),ld(2),ld(3), &
                   &    MCMCParticles(j),ContainerLabel,.TRUE.)
                   or_fail("MCMCEraseCandidatesFromContainers")
                ENDIF
             ENDDO !j=1,nsize
          ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCupdateProposalsAndFilterTopologyInNeighborhood_3d

        FUNCTION MCMCApplyParticle_2d(image_,ipatch,coord,Nm,labels_,MCMCParticle_,DoSimulate) RESULT(info)
          !!! Updates the DS (label image and the regular particle containers) when
          !!! applying a particle. Note that floating particles will not be updated!
          !!! The method only ensures that L and the regular  particles are correct.
          !!! The method expects the label image NOT to be updated already. Else the
          !!! operations performed will be wrong.
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error,ppm_err_sub_failed, &
          &   ppm_err_argument

          USE ppm_rc_module_global, ONLY : MCMCuseBiasedProposal,AllowFission, &
          &   AllowHandles
          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list,MCMCLabelImageHistory
          USE ppm_rc_module_topologicalnumber, ONLY : IsEnclosedByLabel_BGConnectivity
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          REAL(MK), CONTIGUOUS, DIMENSION(:,:), POINTER       :: image_
          INTEGER,                              INTENT(IN   ) :: ipatch
          INTEGER,              DIMENSION(:),   INTENT(IN   ) :: coord
          INTEGER,              DIMENSION(:),   INTENT(IN   ) :: Nm
          INTEGER,  CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
          TYPE(MCMCParticle),                   INTENT(IN   ) :: MCMCParticle_
          LOGICAL,                              INTENT(IN   ) :: DoSimulate
          INTEGER                                             :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0

          INTEGER :: FromLabel
          INTEGER :: ToLabel
          INTEGER :: cLabel

          CHARACTER(LEN=ppm_char) :: caller='MCMCApplyParticle'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          !!! Maintain the regular particle container and the label image
          info=MCMCAddAndRemoveParticlesWhenMove(ipatch,coord,Nm,labels_,MCMCParticle_)
          or_fail("MCMCAddAndRemoveParticlesWhenMove")

          !!! Update the label image. The new point is either a contour particle or 0,
          !!! therefor the negative label value is set.

          cLabel=labels_(coord(1),coord(2))

          FromLabel=ABS(cLabel)
          ToLabel=MCMCParticle_%candlabel

          check_true(<#FromLabel.NE.ToLabel#>,"particle and its candidate have the same label!")

          ALLOCATE(seed,STAT=info)
          or_fail_alloc("seed")
          CALL seed%add(coord(1),coord(2),cLabel)
          CALL MCMCLabelImageHistory(ipatch)%push(seed,info)
          or_fail("MCMCLabelImageHistory(ipatch)%push")

          IF (IsEnclosedByLabel_BGConnectivity(coord,labels_,ToLabel)) THEN
             !!! we created a floating particle
             labels_(coord(1),coord(2))=ToLabel
          ELSE
             !!! standard sign of the label image of boundary particles
             labels_(coord(1),coord(2))=-ToLabel
          ENDIF

          !!! Update the statistics of the propagating and the loser region.
          IF (.NOT.DoSimulate) THEN
             CALL e_data%UpdateStatisticsWhenJump(image_,coord,FromLabel,ToLabel,info)
             or_fail("e_data%UpdateStatisticsWhenJump")
          ENDIF

          !!! Update the proposals for all particles in the neighborhood (as they
          !!! might have changed).
          IF (MCMCuseBiasedProposal.OR..NOT.AllowFission.OR..NOT.AllowHandles) THEN
             info=MCMCupdateProposalsAndFilterTopologyInNeighborhood(ipatch,coord,Nm,labels_)
             or_fail("MCMCupdateProposalsAndFilterTopologyInNeighborhood")
          ENDIF
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCApplyParticle_2d

        FUNCTION MCMCApplyParticle_3d(image_,ipatch,coord,Nm,labels_,MCMCParticle_,DoSimulate) RESULT(info)
          !!! Updates the DS (label image and the regular particle containers) when
          !!! applying a particle. Note that floating particles will not be updated!
          !!! The method only ensures that L and the regular  particles are correct.
          !!! The method expects the label image NOT to be updated already. Else the
          !!! operations performed will be wrong.
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error,ppm_err_sub_failed, &
          &   ppm_err_argument

          USE ppm_rc_module_global, ONLY : MCMCuseBiasedProposal,AllowFission, &
          &   AllowHandles
          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list,MCMCLabelImageHistory
          USE ppm_rc_module_topologicalnumber, ONLY : IsEnclosedByLabel_BGConnectivity
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: image_
          INTEGER,                                INTENT(IN   ) :: ipatch
          INTEGER,              DIMENSION(:),     INTENT(IN   ) :: coord
          INTEGER,              DIMENSION(:),     INTENT(IN   ) :: Nm
          INTEGER,  CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
          TYPE(MCMCParticle),                     INTENT(IN   ) :: MCMCParticle_
          LOGICAL,                                INTENT(IN   ) :: DoSimulate
          INTEGER                                               :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0

          INTEGER :: FromLabel
          INTEGER :: ToLabel
          INTEGER :: cLabel

          CHARACTER(LEN=ppm_char) :: caller='MCMCApplyParticle'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          !!! Maintain the regular particle container and the label image
          info=MCMCAddAndRemoveParticlesWhenMove(ipatch,coord,Nm,labels_,MCMCParticle_)
          or_fail("MCMCAddAndRemoveParticlesWhenMove")

          !!! Update the label image. The new point is either a contour particle or 0,
          !!! therefor the negative label value is set.

          cLabel=labels_(coord(1),coord(2),coord(3))

          FromLabel=ABS(cLabel)
          ToLabel=MCMCParticle_%candlabel

          check_true(<#FromLabel.NE.ToLabel#>,"particle and its candidate have the same label!")

          ALLOCATE(seed,STAT=info)
          or_fail_alloc("seed")
          CALL seed%add(coord(1),coord(2),coord(3),cLabel)
          CALL MCMCLabelImageHistory(ipatch)%push(seed,info)
          or_fail("MCMCLabelImageHistory(ipatch)%push")

          IF (IsEnclosedByLabel_BGConnectivity(coord,labels_,ToLabel)) THEN
             !!! we created a floating particle
             labels_(coord(1),coord(2),coord(3))=ToLabel
          ELSE
             !!! standard sign of the label image of boundary particles
             labels_(coord(1),coord(2),coord(3))=-ToLabel
          ENDIF

          !!! Update the statistics of the propagating and the loser region.
          IF (.NOT.DoSimulate) THEN
             CALL e_data%UpdateStatisticsWhenJump(image_,coord,FromLabel,ToLabel,info)
             or_fail("e_data%UpdateStatisticsWhenJump")
          ENDIF

          !!! Update the proposals for all particles in the neighborhood (as they
          !!! might have changed).
          IF (MCMCuseBiasedProposal.OR..NOT.AllowFission.OR..NOT.AllowHandles) THEN
             info=MCMCupdateProposalsAndFilterTopologyInNeighborhood(ipatch,coord,Nm,labels_)
             or_fail("MCMCupdateProposalsAndFilterTopologyInNeighborhood")
          ENDIF
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCApplyParticle_3d

        FUNCTION MCMCFindParticleB(ppmrcdim) RESULT(info)
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_sub_failed, &
          &   ppm_err_alloc,ppm_err_dealloc
          USE ppm_module_mesh_typedef, ONLY : ppm_t_subpatch_

          USE ppm_rc_module_global, ONLY : one,mesh,image,labels,MCMCstepsize, &
          &   MCMCuseBiasedProposal
          USE ppm_rc_module_rnd, ONLY : GetPartDistrIndex,SaruGetIntegerVariate,  &
          &   DestroyParticlesDiscrDistr,GenerateParticlesFwdProposalsDiscrDistr
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, INTENT(IN   ) :: ppmrcdim
          !!! Problem dimension

          INTEGER                :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle), DIMENSION(:), ALLOCATABLE :: Parts_Q_BgivenA
          TYPE(MCMCParticle), DIMENSION(:), ALLOCATABLE :: Parts_Qb_AgivenB
          TYPE(MCMCParticle)                            :: Part_A
          TYPE(MCMCParticle)                            :: RPart_A
          TYPE(MCMCParticle)                            :: Part_B

          CLASS(ppm_t_subpatch_), POINTER :: sbpitr

          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER     :: image_2d
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER     :: image_3d
          REAL(ppm_kind_double)                               :: t0
          REAL(MK),             DIMENSION(:),     ALLOCATABLE :: ProposalsVector
          REAL(MK)                                            :: Normalizer_Q_B_A
          REAL(MK)                                            :: Normalizer_Qb_A_B
          REAL(MK)                                            :: qb_A_B_unnorm

          INTEGER,             DIMENSION(:,:),   ALLOCATABLE :: Parts_Q_BgivenA_coords
          INTEGER,             DIMENSION(:,:),   ALLOCATABLE :: Parts_Qb_AgivenB_coords
          INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER     :: labels_2d
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER     :: labels_3d
          INTEGER,             DIMENSION(:),     POINTER     :: Nm
          INTEGER,             DIMENSION(ppmrcdim)           :: index_A
          INTEGER,             DIMENSION(ppmrcdim)           :: index_B
          INTEGER                                            :: nsize
          INTEGER                                            :: ipatch
          INTEGER                                            :: PC,i
          INTEGER                                            :: ParticleIndex
          INTEGER                                            :: CandidateLabel

          CHARACTER(LEN=ppm_char) :: caller='MCMCFindParticleB'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          SELECT CASE (ppmrcdim)
          CASE (2)
             NULLIFY(labels_2d,image_2d)
             sbpitr => mesh%subpatch%begin()
             ipatch=1
             DO WHILE (ASSOCIATED(sbpitr))
                Nm => sbpitr%nnodes

                CALL sbpitr%get_field(image,image_2d,info)
                or_fail("Failed to get field image_2d data.",ppm_error=ppm_error_fatal)

                CALL sbpitr%get_field(labels,labels_2d,info)
                or_fail("Failed to get field labels_2d data.",ppm_error=ppm_error_fatal)

                !!! Find particle A:
                DO PC=1,MCMCstepsize
                   index_A=MCMC_CandidateMove_Index(:,PC)

                   info=MCMCApplyParticle_2d(image_2d,ipatch,index_A,Nm,labels_2d,MCMC_CandidateMove(PC),.TRUE.)
                   or_fail("MCMCApplyParticle")

                   CandidateLabel=MCMC_CandidateMove(PC)%candlabel

                   !!! BTW we have to remember if A' is a floating particle in the
                   !!! state x -> A.

                   !!! Get the particles involved in the second step of the
                   !!! proposal (with updated proposals)
                   !!! nsize is at least 1 to count for particle A
                   nsize=MCMCgetParticlesInBGNeighborhood(index_A,Nm,labels_2d,Parts_Q_BgivenA,Parts_Q_BgivenA_coords)

                   !!! Insert particle A in the set for B such that it is possible
                   !!! that we have a single particle move.
                   Parts_Q_BgivenA(nsize)%candlabel=CandidateLabel

                   Part_A=Parts_Q_BgivenA(nsize)

                   !!! Find B

                   !!! Choose B from Q(B|A) and calculate Q(B|A).
                   IF (MCMCuseBiasedProposal) THEN
                      ALLOCATE(ProposalsVector(nsize),STAT=info)
                      or_fail_alloc("ProposalsVector")

                      FORALL (i=1:nsize) ProposalsVector(i)=Parts_Q_BgivenA(i)%proposal

                      Normalizer_Q_B_A=SUM(ProposalsVector)

                      !!! Create a discrete distribution over particles
                      info=GenerateParticlesFwdProposalsDiscrDistr(ProposalsVector,nsize)
                      or_fail("GenerateParticlesFwdProposalsDiscrDistr", &
                      & ppm_error=ppm_error_fatal)

                      ParticleIndex=GetPartDistrIndex()

                      Part_B=Parts_Q_BgivenA(ParticleIndex)

                      !!! The value m_Proposal of B is currently equal to Q_B_A.
                      MCMC_q_B_A(PC)=Part_B%proposal/Normalizer_Q_B_A
                   ELSE
                      ParticleIndex=SaruGetIntegerVariate(nsize)

                      Part_B=Parts_Q_BgivenA(ParticleIndex)

                      !!! The value m_Proposal of B is currently equal to Q_B_A.
                      MCMC_q_B_A(PC)=one/REAL(nsize,MK)
                   ENDIF !MCMCuseBiasedProposal

                   !!! store B (and its original label).
                   MCMC_PartnerMove(PC)=Part_B

                   index_B=Parts_Q_BgivenA_coords(:,ParticleIndex)

                   !!! TOCHECK
                   MCMC_PartnerMove_Index(:,PC)=index_B

                   DEALLOCATE(Parts_Q_BgivenA,Parts_Q_BgivenA_coords,STAT=info)
                   or_fail_dealloc("Parts_Q_BgivenA,Parts_Q_BgivenA_coords")

                   IF (ALL(index_A.EQ.index_B).AND.Part_A%candlabel.EQ.Part_B%candlabel) THEN
                      MCMC_SingleParticleMoveForPairProposals=.TRUE.
                      MCMC_LabelsBeforeJump_B(PC)=MCMC_LabelsBeforeJump_A(PC)
                   ELSE
                      MCMC_LabelsBeforeJump_B(PC)=ABS(labels_2d(index_B(1),index_B(2)))
                   ENDIF

                   !!! Get the reverse particle of A (without proposal update as it
                   !!! is not necessary):
                   RPart_A=MCMCParticle(MCMC_LabelsBeforeJump_A(PC),MCMC_CandidateMove(PC)%proposal)

                   !!! In case that Part_A == Part_B, we must already undo the simulated
                   !!! move in order to calculate Q'(B'|A') (== Q'(A'|B'))
                   IF (MCMC_SingleParticleMoveForPairProposals) THEN
                      info=MCMCApplyParticle_2d(image_2d,ipatch,index_A,Nm,labels_2d,RPart_A,.TRUE.)
                      or_fail("MCMCApplyParticle")
                   ENDIF

                   !!! Get the reverse particle of B:
                   !!! in the current state of the label image and the
                   !!! containers we can calculate qb_A'_B' as well. We assume now
                   !!! that B' was applied and we calculate the probability for A'.
                   nsize=MCMCgetParticlesInFGNeighborhood(index_B,Nm,labels_2d, &
                   &     Parts_Qb_AgivenB,Parts_Qb_AgivenB_coords)

                   IF (MCMCuseBiasedProposal) THEN
                      Normalizer_Qb_A_B=SUM(Parts_Qb_AgivenB(:)%proposal)
                      qb_A_B_unnorm=MCMCproposal(index_A,labels_2d)
                      MCMC_qb_A_B(PC)=qb_A_B_unnorm/Normalizer_Qb_A_B
                   ELSE
                      MCMC_qb_A_B(PC)=one/REAL(nsize,MK)
                   ENDIF

                   DEALLOCATE(Parts_Qb_AgivenB,Parts_Qb_AgivenB_coords,STAT=info)
                   or_fail_dealloc("Parts_Qb_AgivenB,Parts_Qb_AgivenB_coords")

                   IF (.NOT.MCMC_SingleParticleMoveForPairProposals) THEN
                      !!! undo the simulated move.
                      info=MCMCApplyParticle_2d(image_2d,ipatch,index_A,Nm,labels_2d,RPart_A,.TRUE.)
                      or_fail("MCMCApplyParticle")
                      !!! Now we can calculate Q_B (as we now know B and the original
                      !!! state has been recovered).
                      IF (MCMCIsRegularParticle(index_B(1),index_B(2),Part_B,MCMC_LabelsBeforeJump_B(PC))) THEN
                         IF (MCMCuseBiasedProposal) THEN
                            !!! B Proposal needs to be recalculated here:
                            MCMC_q_B(PC)=MCMCproposal(index_B,labels_2d)/ &
                            &            MCMCGetProposalNormalizer(MCMC_LabelsBeforeJump_B(PC),Part_B%candlabel)
                         ELSE
                            MCMC_q_B(PC)=one/MCMCGetProposalNormalizer(MCMC_LabelsBeforeJump_B(PC),Part_B%candlabel)
                         ENDIF
                      ELSE
                         MCMC_q_B(PC)=zero
                      ENDIF
                   ENDIF

                   IF (MCMCuseBiasedProposal) THEN
                      DEALLOCATE(ProposalsVector,STAT=info)
                      or_fail_dealloc("ProposalsVector")

                      !!! make sure that there is no discrete distribution available in memory
                      info=DestroyParticlesDiscrDistr()
                      or_fail("DestroyParticlesDiscrDistr")
                   ENDIF
                ENDDO !PC=1,MCMCstepsize

                sbpitr => mesh%subpatch%next()
                ipatch=ipatch+1
             ENDDO !ASSOCIATED(sbpitr)
             NULLIFY(labels_2d)
          CASE (3)
             NULLIFY(labels_3d,image_3d)
             sbpitr => mesh%subpatch%begin()
             ipatch=1
             DO WHILE (ASSOCIATED(sbpitr))
                CALL sbpitr%get_field(labels,labels_3d,info)
                or_fail("Failed to get field labels_3d data.",ppm_error=ppm_error_fatal)

                CALL sbpitr%get_field(image,image_3d,info)
                or_fail("Failed to get field labels_3d data.",ppm_error=ppm_error_fatal)

                !!! Find particle A:
                DO PC=1,MCMCstepsize
                   index_A=MCMC_CandidateMove_Index(:,PC)

                   info=MCMCApplyParticle_3d(image_3d,ipatch,index_A, &
                   &    Nm,labels_3d,MCMC_CandidateMove(PC),.TRUE.)
                   or_fail("MCMCApplyParticle")

                   CandidateLabel=MCMC_CandidateMove(PC)%candlabel

                   !!! BTW we have to remember if A' is a floating particle in the
                   !!! state x -> A.

                   !!! Get the particles involved in the second step of the
                   !!! proposal (with updated proposals)
                   !!! nsize is at least 1 to count for particle A
                   nsize=MCMCgetParticlesInBGNeighborhood(index_A,Nm,labels_3d, &
                   &     Parts_Q_BgivenA,Parts_Q_BgivenA_coords)

                   !!! Insert particle A in the set for B such that it is possible
                   !!! that we have a single particle move.
                   Parts_Q_BgivenA(nsize)%candlabel=CandidateLabel

                   Part_A=Parts_Q_BgivenA(nsize)

                   !!! Find B

                   !!! Choose B from Q(B|A) and calculate Q(B|A).
                   IF (MCMCuseBiasedProposal) THEN
                      ALLOCATE(ProposalsVector(nsize),STAT=info)
                      or_fail_alloc("ProposalsVector")

                      FORALL (i=1:nsize) ProposalsVector(i)=Parts_Q_BgivenA(i)%proposal

                      Normalizer_Q_B_A=SUM(ProposalsVector)

                      !!! Create a discrete distribution over particles
                      info=GenerateParticlesFwdProposalsDiscrDistr(ProposalsVector,nsize)
                      or_fail("GenerateParticlesFwdProposalsDiscrDistr", &
                      & ppm_error=ppm_error_fatal)

                      ParticleIndex=GetPartDistrIndex()

                      Part_B=Parts_Q_BgivenA(ParticleIndex)

                      !!! The value m_Proposal of B is currently equal to Q_B_A.
                      MCMC_q_B_A(PC)=Part_B%proposal/Normalizer_Q_B_A
                   ELSE
                      ParticleIndex=SaruGetIntegerVariate(nsize)

                      Part_B=Parts_Q_BgivenA(ParticleIndex)

                      !!! The value m_Proposal of B is currently equal to Q_B_A.
                      MCMC_q_B_A(PC)=one/REAL(nsize,MK)
                   ENDIF !MCMCuseBiasedProposal

                   !!! store B (and its original label).
                   MCMC_PartnerMove(PC)=Part_B

                   index_B=Parts_Q_BgivenA_coords(:,ParticleIndex)

                   !!! TOCHECK
                   MCMC_PartnerMove_Index(:,PC)=index_B

                   DEALLOCATE(Parts_Q_BgivenA,Parts_Q_BgivenA_coords,STAT=info)
                   or_fail_dealloc("Parts_Q_BgivenA,Parts_Q_BgivenA_coords")

                   IF (ALL(index_A.EQ.index_B).AND.Part_A%candlabel.EQ.Part_B%candlabel) THEN
                      MCMC_SingleParticleMoveForPairProposals=.TRUE.
                      MCMC_LabelsBeforeJump_B(PC)=MCMC_LabelsBeforeJump_A(PC)
                   ELSE
                      MCMC_LabelsBeforeJump_B(PC)=ABS(labels_3d(index_B(1),index_B(2),index_B(3)))
                   ENDIF

                   !!! Get the reverse particle of A (without proposal update as it
                   !!! is not necessary):
                   RPart_A=MCMCParticle(MCMC_LabelsBeforeJump_A(PC),MCMC_CandidateMove(PC)%proposal)

                   !!! In case that Part_A == Part_B, we must already undo the simulated
                   !!! move in order to calculate Q'(B'|A') (== Q'(A'|B'))
                   IF (MCMC_SingleParticleMoveForPairProposals) THEN
                      info=MCMCApplyParticle_3d(image_3d,ipatch,index_A,Nm,labels_3d,RPart_A,.TRUE.)
                      or_fail("MCMCApplyParticle")
                   ENDIF

                   !!! Get the reverse particle of B:
                   !!! in the current state of the label image and the
                   !!! containers we can calculate qb_A'_B' as well. We assume now
                   !!! that B' was applied and we calculate the probability for A'.
                   nsize=MCMCgetParticlesInFGNeighborhood(index_B,Nm,labels_3d, &
                   &     Parts_Qb_AgivenB,Parts_Qb_AgivenB_coords)

                   IF (MCMCuseBiasedProposal) THEN
                      Normalizer_Qb_A_B=SUM(Parts_Qb_AgivenB(:)%proposal)
                      qb_A_B_unnorm=MCMCproposal(index_A,labels_3d)
                      MCMC_qb_A_B(PC)=qb_A_B_unnorm/Normalizer_Qb_A_B
                   ELSE
                      MCMC_qb_A_B(PC)=one/REAL(nsize,MK)
                   ENDIF

                   DEALLOCATE(Parts_Qb_AgivenB,Parts_Qb_AgivenB_coords,STAT=info)
                   or_fail_dealloc("Parts_Qb_AgivenB,Parts_Qb_AgivenB_coords")

                   IF (.NOT.MCMC_SingleParticleMoveForPairProposals) THEN
                      !!! undo the simulated move.
                      info=MCMCApplyParticle_3d(image_3d,ipatch,index_A,Nm,labels_3d,RPart_A,.TRUE.)
                      or_fail("MCMCApplyParticle")
                      !!! Now we can calculate Q_B (as we now know B and the original
                      !!! state has been recovered).
                      IF (MCMCIsRegularParticle(index_B(1),index_B(2),index_B(3),Part_B,MCMC_LabelsBeforeJump_B(PC))) THEN
                         IF (MCMCuseBiasedProposal) THEN
                            !!! B Proposal needs to be recalculated here:
                            MCMC_q_B(PC)=MCMCproposal(index_B,labels_3d)/ &
                            &            MCMCGetProposalNormalizer(MCMC_LabelsBeforeJump_B(PC),Part_B%candlabel)
                         ELSE
                            MCMC_q_B(PC)=one/MCMCGetProposalNormalizer(MCMC_LabelsBeforeJump_B(PC),Part_B%candlabel)
                         ENDIF
                      ELSE
                         MCMC_q_B(PC)=zero
                      ENDIF
                   ENDIF

                   IF (MCMCuseBiasedProposal) THEN
                      DEALLOCATE(ProposalsVector,STAT=info)
                      or_fail_dealloc("ProposalsVector")

                      !!! make sure that there is no discrete distribution available in memory
                      info=DestroyParticlesDiscrDistr()
                      or_fail("DestroyParticlesDiscrDistr")
                   ENDIF
                ENDDO !PC=1,MCMCstepsize

                sbpitr => mesh%subpatch%next()
                ipatch=ipatch+1
             ENDDO
             NULLIFY(labels_3d,image_3d)
          END SELECT
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCFindParticleB

        LOGICAL FUNCTION MCMCMove(ppmrcdim,info)
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_sub_failed, &
          &   ppm_err_dealloc
          USE ppm_module_mesh_typedef, ONLY : ppm_t_subpatch_

          USE ppm_rc_module_global, ONLY : one,mesh,image,labels,MCMCstepsize, &
          &   MCMCuseBiasedProposal,MCMCusePairProposal
          USE ppm_rc_module_rnd, ONLY : GetPartDistrIndex,SaruGetIntegerVariate,  &
          &   DestroyParticlesDiscrDistr,GenerateParticlesFwdProposalsDiscrDistr
          USE ppm_rc_module_linkedlist, ONLY : MCMCAppliedParticleOrigLabels, &
          &   ppm_rc_link
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, INTENT(IN   ) :: ppmrcdim
          !!! Problem dimension

          INTEGER, INTENT(  OUT) :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle), DIMENSION(:), ALLOCATABLE :: Parts_Q_AgivenB
          TYPE(MCMCParticle), DIMENSION(:), ALLOCATABLE :: Parts_Qb_BgivenA
          TYPE(MCMCParticle)                            :: Part_A
          TYPE(MCMCParticle)                            :: Part_B
          TYPE(MCMCParticle)                            :: ReverseFloatingP
!
          CLASS(ppm_t_subpatch_), POINTER :: sbpitr
!
! !           TYPE(ppm_rc_link), POINTER :: seedlink
!
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER     :: image_2d
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER     :: image_3d
          REAL(ppm_kind_double)                               :: t0
          REAL(MK)                                            :: Normalizer_Q_A_B
          REAL(MK)                                            :: Normalizer_Qb_B_A
          REAL(MK)                                            :: Normalizer_qb_A
          REAL(MK)                                            :: Normalizer_qb_B
!           REAL(MK)                                            :: ForwardBackwardRatio
          REAL(MK)                                            :: etemp
          REAL(MK)                                            :: proposalA
          REAL(MK)                                            :: proposalBb
!           REAL(ppm_kind_double)                               :: ProposeAFloatInXb
!           REAL(MK)                                            :: HastingsRatio

          INTEGER,             DIMENSION(:,:),   ALLOCATABLE :: Parts_Q_AgivenB_coords
          INTEGER,             DIMENSION(:,:),   ALLOCATABLE :: Parts_Qb_BgivenA_coords
          INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER     :: labels_2d
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER     :: labels_3d
          INTEGER,             DIMENSION(:),     POINTER     :: Nm
!           INTEGER,             DIMENSION(:),     POINTER     :: seedn
          INTEGER,             DIMENSION(ppmrcdim)           :: index_A
          INTEGER,             DIMENSION(ppmrcdim)           :: index_B
          INTEGER,             DIMENSION(ppmrcdim)           :: index_
          INTEGER                                            :: nsize
          INTEGER                                            :: ipatch
          INTEGER                                            :: PC
          INTEGER                                            :: OriginalLabel

          LOGICAL :: e_merge
!           LOGICAL :: Accept

          CHARACTER(LEN=ppm_char) :: caller='MCMCMoves'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          e_merge=.FALSE.

          SELECT CASE (ppmrcdim)
          CASE (2)
             NULLIFY(labels_2d,image_2d)
             sbpitr => mesh%subpatch%begin()
             ipatch=1
             DO WHILE (ASSOCIATED(sbpitr))
                Nm => sbpitr%nnodes

                CALL sbpitr%get_field(image,image_2d,info)
                or_fail("Failed to get field image_2d data.",ppm_error=ppm_error_fatal)

                CALL sbpitr%get_field(labels,labels_2d,info)
                or_fail("Failed to get field labels_2d data.",ppm_error=ppm_error_fatal)

                IF (MCMCusePairProposal.AND..NOT.MCMC_SingleParticleMoveForPairProposals) THEN
                   !!! Iterate the candidates, calculate the energy and perform the moves.
                   DO PC=1,MCMCstepsize
                      !!! apply particle A and B, start with B:
                      !!! it is necessary that we start with particle B as we have
                      !!! to calculate Q(A|B) and Qb(B|A).
                      Part_B =MCMC_PartnerMove(PC)
                      OriginalLabel=MCMC_LabelsBeforeJump_B(PC)
                      index_B=MCMC_PartnerMove_Index(:,PC)
                      !!! We calculate the energy and apply them move if
                      !!! - the move has not been performed beforehand (a particle
                      !!!   was sampled twice)
                      !!! - THE FOLLOWING IS CURRENTLY A NON ISSUE AS B CANNOT BE A':
                      !!!   particle B is not the reverse particle of particle A (in
                      !!!   case of m_usePairProposals, i.e. vN > 1). This is important
                      !!!   because the energy update gets corrupted as particle B is
                      !!!   not a valid particle before A was applied. To resolve this
                      !!!   we just perform a one particle move (we don't apply B).
                      IF (.NOT.MCMCAppliedParticles%containselement(index_B(1),index_B(2),Part_B)) THEN
                         !!! Calculate the energy difference when changing this candidate:
                         etemp=e_data%EvaluateEnergyDifference(image_2d,labels_2d,   &
                         &     index_B,OriginalLabel,Part_B%candlabel,e_merge)       &
                         &    +e_length%EvaluateEnergyDifference(image_2d,labels_2d, &
                         &     index_B,OriginalLabel,Part_B%candlabel,e_merge)

                         MCMCTotEnergyDiff=MCMCTotEnergyDiff+etemp

                         !!! Finally, perform the (particle-)move
                         info=MCMCApplyParticle_2d(image_2d,ipatch,index_B,Nm,labels_2d,Part_B,.FALSE.)
                         or_fail("MCMCApplyParticle!")

                         CALL MCMCAppliedParticles%insert(index_B(1),index_B(2),Part_B,info)
                         or_fail("MCMCAppliedParticles%insert")

                         CALL MCMCAppliedParticleOrigLabels(ipatch)%add(index_B(1),index_B(2),OriginalLabel)
                      ENDIF
                      !!! Calculate Q(A|B) and Qb(B|A) in case we moved B only; this is
                      !!! when vN == 2.
                      !!! Get the neighbors (conditional particles) and sum up
                      !!! their proposal values; this is the normalizer for the
                      !!! discrete probability Q(A|B)
                      nsize=MCMCgetParticlesInFGNeighborhood(index_B,Nm,labels_2d, &
                      &     Parts_Q_AgivenB,Parts_Q_AgivenB_coords)
                      !!! add particle B as this is always a candidate as well

                      index_A=MCMC_CandidateMove_Index(:,PC)
                      !!! Part_A%Proposal is not valid anymore. Part_A
                      !!! got a new proposal when applying particle B.
                      proposalA=MCMCproposal(index_A,labels_2d)

                      IF (MCMCuseBiasedProposal) THEN
                         Normalizer_Q_A_B=SUM(Parts_Q_AgivenB(:)%proposal)
                         MCMC_q_A_B(PC)=proposalA/Normalizer_Q_A_B
                      ELSE
                         MCMC_q_A_B(PC)=one/REAL(nsize,MK)
                      ENDIF !MCMCuseBiasedProposal

                      !!! create A'
                      !!! Calculate Qb(B'|A')
                      nsize=MCMCgetParticlesInFGNeighborhood(index_A,Nm,labels_2d, &
                      &     Parts_Qb_BgivenA,Parts_Qb_BgivenA_coords)

                      IF (MCMCuseBiasedProposal) THEN
                         Normalizer_Qb_B_A=SUM(Parts_Qb_BgivenA(:)%proposal)

                         proposalBb=MCMCproposal(index_B,labels_2d)

                         MCMC_qb_B_A(PC)=proposalBb/Normalizer_Qb_B_A
                      ELSE
                         MCMC_qb_B_A(PC)=one/REAL(nsize,MK)
                      ENDIF !MCMCuseBiasedProposal

                      !!! apply particle A
                      Part_A =MCMC_CandidateMove(PC)
                      OriginalLabel=MCMC_LabelsBeforeJump_A(PC)
                      index_A=MCMC_CandidateMove_Index(:,PC)

                      !!! We calculate the energy and apply them move if
                      !!! - the move has not been performed beforehand (a particle
                      !!!   was sampled twice)
                      !!! - THE FOLLOWING IS CURRENTLY A NON ISSUE AS B CANNOT BE A':
                      !!!   particle B is not the reverse particle of particle A (in
                      !!!   case of m_usePairProposals, i.e. vN > 1). This is important
                      !!!   because the energy update gets corrupted as particle B is
                      !!!   not a valid particle before A was applied. To resolve this
                      !!!   we just perform a one particle move (we don't apply B).
                      IF (.NOT.MCMCAppliedParticles%containselement(index_A(1),index_A(2),Part_A)) THEN
                         !!! Calculate the energy difference when changing this candidate:
                         etemp=e_data%EvaluateEnergyDifference(image_2d,labels_2d,   &
                         &     index_A,OriginalLabel,Part_A%candlabel,e_merge)       &
                         &    +e_length%EvaluateEnergyDifference(image_2d,labels_2d, &
                         &     index_A,OriginalLabel,Part_A%candlabel,e_merge)

                         MCMCTotEnergyDiff=MCMCTotEnergyDiff+etemp
                         !!! Finally, perform the (particle-)move
                         info=MCMCApplyParticle_2d(image_2d,ipatch,index_A,Nm,labels_2d,Part_A,.FALSE.)
                         or_fail("MCMCApplyParticle!")

                         CALL MCMCAppliedParticles%insert(index_A(1),index_A(2),Part_A,info)
                         or_fail("MCMCAppliedParticles%insert")

                         CALL MCMCAppliedParticleOrigLabels(ipatch)%add(index_A(1),index_A(2),OriginalLabel)
                      ENDIF
                   ENDDO !PC=1,MCMCstepsize
                ELSE
                   !!! Iterate the candidates, calculate the energy and perform the moves.
                   DO PC=1,MCMCstepsize
                      Part_A = MCMC_CandidateMove(PC)
                      !!! apply particle A and B, start with B:
                      !!! it is necessary that we start with particle B as we have
                      !!! to calculate Q(A|B) and Qb(B|A).
                      OriginalLabel=MCMC_LabelsBeforeJump_A(PC)
                      index_A=MCMC_CandidateMove_Index(:,PC)

                      !!! We calculate the energy and apply them move if
                      !!! - the move has not been performed beforehand (a particle
                      !!!   was sampled twice)
                      !!! - THE FOLLOWING IS CURRENTLY A NON ISSUE AS B CANNOT BE A':
                      !!!   particle B is not the reverse particle of particle A (in
                      !!!   case of m_usePairProposals, i.e. vN > 1). This is important
                      !!!   because the energy update gets corrupted as particle B is
                      !!!   not a valid particle before A was applied. To resolve this
                      !!!   we just perform a one particle move (we don't apply B).
                      IF (.NOT.MCMCAppliedParticles%containselement(index_A(1),index_A(2),Part_A)) THEN
                         !!! Calculate the energy difference when changing this candidate:
                         etemp=e_data%EvaluateEnergyDifference(image_2d,labels_2d,   &
                         &     index_A,OriginalLabel,Part_A%candlabel,e_merge)       &
                         &    +e_length%EvaluateEnergyDifference(image_2d,labels_2d, &
                         &     index_A,OriginalLabel,Part_A%candlabel,e_merge)

                         MCMCTotEnergyDiff=MCMCTotEnergyDiff+etemp
                         !!! Finally, perform the (particle-)move
                         info=MCMCApplyParticle_2d(image_2d,ipatch,index_A,Nm,labels_2d,Part_A,.FALSE.)
                         or_fail("MCMCApplyParticle!")

                         CALL MCMCAppliedParticles%insert(index_A(1),index_A(2),Part_A,info)
                         or_fail("MCMCAppliedParticles%insert")

                         CALL MCMCAppliedParticleOrigLabels(ipatch)%add(index_A(1),index_A(2),OriginalLabel)
                      ENDIF
                   ENDDO !PC=1,MCMCstepsize
                ENDIF !MCMCusePairProposal.AND..NOT.MCMC_SingleParticleMoveForPairProposals

                !!! Free memory
                IF (ALLOCATED(Parts_Q_AgivenB)) THEN
                   DEALLOCATE(Parts_Q_AgivenB,Parts_Q_AgivenB_coords,STAT=info)
                   or_fail_dealloc("Parts_Q_AgivenB & Parts_Q_AgivenB_coords")
                ENDIF
                IF (ALLOCATED(Parts_Qb_BgivenA)) THEN
                   DEALLOCATE(Parts_Qb_BgivenA,Parts_Qb_BgivenA_coords,STAT=info)
                   or_fail_dealloc("Parts_Qb_BgivenA & Parts_Qb_BgivenA_coords")
                ENDIF

                DO PC=1,MCMCstepsize
                   !!! Correct the containers whenever floating particles were involved:
                   !!! The method moveParticles, for simplicity, only works on the regular
                   !!! particle set.

                   !!! First, figure out if A' or B' is floating:
                   !!! Figure out if the backward particles are floating:
                   IF (MCMCusePairProposal) THEN
                      index_=MCMC_PartnerMove_Index(:,PC)

                   ELSE
                      index_=MCMC_CandidateMove_Index(:,PC)
                      !!! if we're not in pair proposal mode we did not yet check if
                      !!! A's reverse particle is floating (else we did already):
                   ENDIF

                   !!! the first condition is needed when not using pair proposal mode
                   IF (MCMC_Particle_Ab_IsFloating(PC)) THEN
                      ReverseFloatingP=MCMCParticle(MCMC_LabelsBeforeJump_A(PC),MCMCproposal(index_,labels_2d))
                   ENDIF

                   !!! in pair proposal, if A' is floating, B' is as well (they are the
                   !!! same particle):
                   IF (MCMCusePairProposal.AND.MCMC_Particle_Bb_IsFloating(PC)) THEN
                      !!! only possible in pair proposal mode
                      ReverseFloatingP=MCMCParticle(MCMC_LabelsBeforeJump_B(PC),MCMCproposal(index_,labels_2d))
                   ENDIF

                   !!! finally convert the regular particle into a floating particle,
                   !!! i.e. insert it in the floating DS and remove it from the regular:
                   IF (MCMC_Particle_Ab_IsFloating(PC).OR.MCMC_Particle_Bb_IsFloating(PC)) THEN
                      !!! insert the reverse particle in the appropriate container. If
                      !!! there is no space, we reject the move.
                      IF (.NOT.MCMCInsertFloatingParticle(ipatch,index_(1),index_(2),ReverseFloatingP,.TRUE.)) THEN
                         !!! TODO: calling MCMCReject from here invalidates the result.
                         !!! MCMCReject(&vAppliedParticles,&vAppliedParticleOrigLabels);
                         MCMCHardReject=.TRUE.
                      ENDIF
                   ENDIF
                ENDDO !PC=1,MCMCstepsize

                !!! We are now in the state x'.
                !!! Calculate Q'(A) and maybe Q'(B). Note that this has to be done after
                !!! all particles were applied.
                DO PC=1,MCMCstepsize
                   !!! Calculate MCMC_qb_A and MCMC_qb_B
                   IF (.NOT.MCMCuseBiasedProposal) THEN
                      MCMC_qb_A(PC)=one
                      MCMC_qb_B(PC)=one
                   ELSE
                      index_A=MCMC_CandidateMove_Index(:,PC)
                      MCMC_qb_A(PC)=MCMCproposal(index_A,labels_2d)
                      IF (MCMCusePairProposal.AND..NOT.MCMC_SingleParticleMoveForPairProposals) THEN
                         index_B=MCMC_PartnerMove_Index(:,PC)
                         MCMC_qb_B(PC)=MCMCproposal(index_B,labels_2d)
                      ENDIF
                   ENDIF

                   !!! Normalize MCMC_qb_A and MCMC_qb_B
                   IF (MCMC_Particle_Ab_IsFloating(PC)) THEN
                      Normalizer_qb_A=REAL(MCMCFloatingParticlesProposalNormalizer,MK)
                   ELSE
                      Normalizer_qb_A=MCMCGetProposalNormalizer(MCMC_CandidateMove(PC)%candlabel,MCMC_LabelsBeforeJump_A(PC))
                   ENDIF
                   MCMC_qb_A(PC)=MCMC_qb_A(PC)/Normalizer_qb_A

                   IF (MCMCusePairProposal.AND..NOT.MCMC_SingleParticleMoveForPairProposals) THEN
                      IF (MCMC_Particle_Bb_IsFloating(PC)) THEN
                         Normalizer_qb_B=REAL(MCMCFloatingParticlesProposalNormalizer,MK)
                      ELSE
                         Normalizer_qb_B=MCMCGetProposalNormalizer(MCMC_PartnerMove(PC)%candlabel,MCMC_LabelsBeforeJump_B(PC))
                      ENDIF
                      MCMC_qb_B(PC)=MCMC_qb_B(PC)/Normalizer_qb_B
                   ENDIF
                   !!! Finally, we omit half of the calculations if particle A == B
                   IF (MCMC_SingleParticleMoveForPairProposals) THEN
                      MCMC_q_B(PC)   =MCMC_q_A(PC)
                      MCMC_q_A_B(PC) =MCMC_q_B_A(PC);
                      MCMC_qb_B_A(PC)=MCMC_qb_A_B(PC);
                      MCMC_qb_B(PC)  =MCMC_qb_A(PC);
                   ENDIF
                ENDDO !PC=1,MCMCstepsize

!                 !!! Calculate the forward-backward ratio:
!                 ForwardBackwardRatio=one
!                 DO PC=1,MCMCstepsize
!                    IF (MCMCusePairProposal) THEN
!                       IF (MCMC_Particle_Ab_IsFloating(PC).OR.MCMC_Particle_Bb_IsFloating(PC).OR.MCMC_Particle_A_IsFloating(PC)) THEN
!                          ForwardBackwardRatio=ForwardBackwardRatio* &
!                          & (MCMC_qb_B(PC)*MCMC_qb_A_B(PC))/(MCMC_q_A(PC)*MCMC_q_B_A(PC))
!                       ELSE
!                          ForwardBackwardRatio=ForwardBackwardRatio* &
!                          & (MCMC_qb_A(PC)*MCMC_qb_B_A(PC)+MCMC_qb_B(PC)*MCMC_qb_A_B(PC))/(MCMC_q_A(PC)*MCMC_q_B_A(PC)+MCMC_q_B(PC)*MCMC_q_A_B(PC))
!                       ENDIF
!                    ELSE
!                       IF (MCMC_Particle_A_IsFloating(PC)) THEN
!                          !!! we distroy a floating particle, in the next iteration there
!                          !!! will be one floating particle less, hence the probability
!                          !!! in the x' to sample a floating particle is (note that
!                          !!! both normalizers are in state x'):
!                          ProposeAFloatInXb=MCMCTotalNormalizer/(MCMCFloatingParticlesProposalNormalizer+MCMCTotalNormalizer)
!                          ProposeAFloatInXb=halfd*ProposeAFloatInXb/MCMCProbabilityToProposeAFloatingParticle
!
!                          ForwardBackwardRatio=ForwardBackwardRatio*REAL(ProposeAFloatInXb,MK)*MCMC_qb_A(PC)/MCMC_q_A(PC)
!                       ELSE IF (MCMC_Particle_Ab_IsFloating(PC)) THEN
!                          !!! we create a floating particle, in the next iteration there
!                          !!! will be one floating particle more, hence the probability
!                          !!! in the x' to sample a floating particle is (note that
!                          !!! m_MCMCTotalNormalizer is updated to x'):
!                          ProposeAFloatInXb=MCMCFloatingParticlesProposalNormalizer/(MCMCFloatingParticlesProposalNormalizer+MCMCTotalNormalizer)
!                          ProposeAFloatInXb=ProposeAFloatInXb/(halfd*(oned-MCMCProbabilityToProposeAFloatingParticle))
!
!                          ForwardBackwardRatio=ForwardBackwardRatio*REAL(ProposeAFloatInXb,MK)*MCMC_qb_A(PC)/MCMC_q_A(PC)
!                       ELSE
!                          !!! Shrinkage and growth events have the same probability.
!                          !!! We hence only need to compare the individual particle
!                          !!! ratios.
!                          ForwardBackwardRatio=ForwardBackwardRatio*MCMC_qb_A(PC)/MCMC_q_A(PC)
!                       ENDIF
!                    ENDIF
!                 ENDDO !PC=1,MCMCstepsize
!
!                 !!! Compute the Hastingsratio:
!                 HastingsRatio = EXP(-MCMCTotEnergyDiff/MCMCtemperature)*ForwardBackwardRatio
!                 !!! debug: check if the average of the FBR is equal to one.
!                 !!! HastingsRatio = ForwardBackwardRatio
!                 !!! Should I stay or shoud I go; the Metropolis-Hastings algorithm:
!                 IF (HastingsRatio.GE.one) THEN
!                    Accept=.TRUE.
!                 ELSE
!                    IF (HastingsRatio.GT.SaruGetRealVariate()) THEN
!                       Accept=.TRUE.
!                    ELSE
!                       Accept=.FALSE.
!                    ENDIF
!                 ENDIF
!
!                 !!! Register the result (if we accept) or rollback to the previous state.
!                 IF (Accept.AND..NOT.MCMCHardReject) THEN
!                    seedlink => MCMCAppliedParticleOrigLabels(ipatch)%first
!                    DO WHILE (ASSOCIATED(seedlink))
!                       seedn => seedlnk%getValue()
!                       !!! store the results and finish (next iteration).
!                       Part_A=MCMCAppliedParticles%search(seedn(1),seedn(2))
!                       !MCMCStore index,originallabel,step
!                       !!! something changed, particles were inserted and deleted. We need
!                       !!! to update the edge-map constant.
!                       !MCMCUpdateRegularParticleMapInNeighborhood(vAppliedParticles[vM].m_Index);
!                       seedlnk => seedlnk%nextLink()
!                    ENDDO
!                 ELSE
!                    !MCMCReject(&vAppliedParticles,&vAppliedParticleOrigLabels);
!                 ENDIF
!
!                 MCMCMove=Accept

                sbpitr => mesh%subpatch%next()
                ipatch=ipatch+1
             ENDDO !ASSOCIATED(sbpitr)
             NULLIFY(labels_2d,image_2d)
          CASE (3)
             NULLIFY(labels_3d,image_3d)
             sbpitr => mesh%subpatch%begin()
             ipatch=1
             DO WHILE (ASSOCIATED(sbpitr))
                CALL sbpitr%get_field(labels,labels_3d,info)
                or_fail("Failed to get field labels_3d data.",ppm_error=ppm_error_fatal)

                CALL sbpitr%get_field(image,image_3d,info)
                or_fail("Failed to get field labels_3d data.",ppm_error=ppm_error_fatal)

                IF (MCMCusePairProposal.AND..NOT.MCMC_SingleParticleMoveForPairProposals) THEN
                   !!! Iterate the candidates, calculate the energy and perform the moves.
                   DO PC=1,MCMCstepsize
                      !!! apply particle A and B, start with B:
                      Part_B=MCMC_PartnerMove(PC)
                      !!! it is necessary that we start with particle B as we have
                      !!! to calculate Q(A|B) and Qb(B|A).
                      OriginalLabel=MCMC_LabelsBeforeJump_B(PC)
                      index_B=MCMC_PartnerMove_Index(:,PC)
                      !!! We calculate the energy and apply the move if
                      !!! - the move has not been performed beforehand (a particle
                      !!!   was sampled twice)
                      !!! - THE FOLLOWING IS CURRENTLY A NON ISSUE AS B CANNOT BE A':
                      !!!   particle B is not the reverse particle of particle A (in
                      !!!   case of m_usePairProposals, i.e. vN > 1). This is important
                      !!!   because the energy update gets corrupted as particle B is
                      !!!   not a valid particle before A was applied. To resolve this
                      !!!   we just perform a one particle move (we don't apply B).
                      IF (.NOT.MCMCAppliedParticles%containselement(index_B(1),index_B(2),index_B(3),Part_B)) THEN
                         !!! Calculate the energy difference when changing this candidate:
                         etemp=e_data%EvaluateEnergyDifference(image_3d,labels_3d,   &
                         &     index_B,OriginalLabel,Part_B%candlabel,e_merge)       &
                         &    +e_length%EvaluateEnergyDifference(image_3d,labels_3d, &
                         &     index_B,OriginalLabel,Part_B%candlabel,e_merge)

                         MCMCTotEnergyDiff=MCMCTotEnergyDiff+etemp
                         !!! Finally, perform the (particle-)move
                         info=MCMCApplyParticle_3d(image_3d,ipatch,index_B,Nm,labels_3d,Part_B,.FALSE.)
                         or_fail("MCMCApplyParticle")

                         CALL MCMCAppliedParticles%insert(index_B(1),index_B(2),index_B(3),Part_B,info)
                         or_fail("MCMCAppliedParticles%insert")

                         CALL MCMCAppliedParticleOrigLabels(ipatch)%add(index_B(1),index_B(2),index_B(3),OriginalLabel)
                      ENDIF
                      !!! Calculate Q(A|B) and Qb(B|A) in case we moved B only;
                      !!! Get the neighbors (conditional particles) and sum up
                      !!! their proposal values; this is the normalizer for the
                      !!! discrete probability Q(A|B)
                      nsize=MCMCgetParticlesInFGNeighborhood(index_B,Nm,labels_3d, &
                      &     Parts_Q_AgivenB,Parts_Q_AgivenB_coords)
                      Part_B%proposal=Parts_Q_AgivenB(nsize)%proposal

                      index_A=MCMC_CandidateMove_Index(:,PC)
                      !!! Part_A%proposal is not valid anymore. Particle A
                      !!! got a new proposal when applying particle B.
                      proposalA=MCMCproposal(index_A,labels_3d)

                      IF (MCMCuseBiasedProposal) THEN
                         Normalizer_Q_A_B=SUM(Parts_Q_AgivenB(:)%proposal)
                         MCMC_q_A_B(PC)=proposalA/Normalizer_Q_A_B
                      ELSE
                         MCMC_q_A_B(PC)=one/REAL(nsize,MK)
                      ENDIF

                      !!! create A'
                      !!! Calculate Qb(B'|A')
                      nsize=MCMCgetParticlesInFGNeighborhood(index_A,Nm,labels_3d, &
                      &     Parts_Qb_BgivenA,Parts_Qb_BgivenA_coords)

                      IF (MCMCuseBiasedProposal) THEN
                         Normalizer_Qb_B_A=SUM(Parts_Qb_BgivenA(:)%proposal)
                         !!! the proposal of the backward particle (given A) is:
                         proposalBb=MCMCproposal(index_B,labels_3d)

                         MCMC_qb_B_A(PC)=proposalBb/Normalizer_Qb_B_A
                      ELSE
                         MCMC_qb_B_A(PC)=one/REAL(nsize,MK)
                      ENDIF

                      !!! apply particle A
                      Part_A=MCMC_CandidateMove(PC)
                      OriginalLabel=MCMC_LabelsBeforeJump_A(PC)
                      index_A=MCMC_CandidateMove_Index(:,PC)
                      !!! We calculate the energy and apply them move iff
                      !!! - the move has not been performed beforehand (a particle
                      !!!   was sampled twice)
                      !!! - THE FOLLOWING IS CURRENTLY A NON ISSUE AS B CANNOT BE A':
                      !!!   particle B is not the reverse particle of particle A (in
                      !!!   case of m_usePairProposals, i.e. vN > 1). This is important
                      !!!   because the energy update gets corrupted as particle B is
                      !!!   not a valid particle before A was applied. To resolve this
                      !!!   we just perform a one particle move (we don't apply B).
                      IF (.NOT.MCMCAppliedParticles%containselement(index_A(1),index_A(2),index_A(3),Part_A)) THEN
                         !!! Calculate the energy difference when changing this candidate:
                         etemp=e_data%EvaluateEnergyDifference(image_3d,labels_3d,   &
                         &     index_A,OriginalLabel,Part_A%candlabel,e_merge)       &
                         &    +e_length%EvaluateEnergyDifference(image_3d,labels_3d, &
                         &     index_A,OriginalLabel,Part_A%candlabel,e_merge)

                         MCMCTotEnergyDiff=MCMCTotEnergyDiff+etemp
                         !!! Finally, perform the (particle-)move
                         info=MCMCApplyParticle_3d(image_3d,ipatch,index_A,Nm,labels_3d,Part_A,.FALSE.)
                         or_fail("MCMCApplyParticle")

                         CALL MCMCAppliedParticles%insert(index_A(1),index_A(2),index_A(3),Part_A,info)
                         or_fail("MCMCAppliedParticles%insert")

                         CALL MCMCAppliedParticleOrigLabels(ipatch)%add(index_A(1),index_A(2),index_A(3),OriginalLabel)
                      ENDIF
                   ENDDO !PC=1,MCMCstepsize
                ELSE
                   !!! Iterate the candidates, calculate the energy and perform the moves.
                   DO PC=1,MCMCstepsize
                      !!! apply particle A:
                      Part_A=MCMC_CandidateMove(PC)
                      OriginalLabel=MCMC_LabelsBeforeJump_A(PC)
                      index_A=MCMC_CandidateMove_Index(:,PC)
                      !!! We calculate the energy and apply the move if
                      !!! - the move has not been performed beforehand (a particle
                      !!!   was sampled twice)
                      !!! - THE FOLLOWING IS CURRENTLY A NON ISSUE AS B CANNOT BE A':
                      !!!   particle B is not the reverse particle of particle A (in
                      !!!   case of m_usePairProposals, i.e. vN > 1). This is important
                      !!!   because the energy update gets corrupted as particle B is
                      !!!   not a valid particle before A was applied. To resolve this
                      !!!   we just perform a one particle move (we don't apply B).
                      IF (.NOT.MCMCAppliedParticles%containselement(index_A(1),index_A(2),index_A(3),Part_A)) THEN
                         !!! Calculate the energy difference when changing this candidate:
                         etemp=e_data%EvaluateEnergyDifference(image_3d,labels_3d,   &
                         &     index_A,OriginalLabel,Part_A%candlabel,e_merge)       &
                         &    +e_length%EvaluateEnergyDifference(image_3d,labels_3d, &
                         &     index_A,OriginalLabel,Part_A%candlabel,e_merge)

                         MCMCTotEnergyDiff=MCMCTotEnergyDiff+etemp
                         !!! Finally, perform the (particle-)move
                         info=MCMCApplyParticle_3d(image_3d,ipatch,index_A,Nm,labels_3d,Part_A,.FALSE.)
                         or_fail("MCMCApplyParticle")

                         CALL MCMCAppliedParticles%insert(index_A(1),index_A(2),index_A(3),Part_A,info)
                         or_fail("MCMCAppliedParticles%insert")

                         CALL MCMCAppliedParticleOrigLabels(ipatch)%add(index_A(1),index_A(2),index_A(3),OriginalLabel)
                      ENDIF
                   ENDDO !PC=1,MCMCstepsize
                ENDIF !MCMCusePairProposal.AND..NOT.MCMC_SingleParticleMoveForPairProposals

                !!! Free memory
                IF (ALLOCATED(Parts_Q_AgivenB)) THEN
                   DEALLOCATE(Parts_Q_AgivenB,Parts_Q_AgivenB_coords,STAT=info)
                   or_fail_dealloc("Parts_Q_AgivenB & Parts_Q_AgivenB_coords")
                ENDIF
                IF (ALLOCATED(Parts_Qb_BgivenA)) THEN
                   DEALLOCATE(Parts_Qb_BgivenA,Parts_Qb_BgivenA_coords,STAT=info)
                   or_fail_dealloc("Parts_Qb_BgivenA & Parts_Qb_BgivenA_coords")
                ENDIF

                DO PC=1,MCMCstepsize
                   !!! Correct the containers whenever floating particles were involved:
                   !!! The method moveParticles, for simplicity, only works on the regular
                   !!! particle set.

                   !!! First, figure out if A' or B' is floating:
                   !!! Figure out if the backward particles are floating:
                   IF (MCMCusePairProposal) THEN
                      index_=MCMC_PartnerMove_Index(:,PC)

                   ELSE
                      index_=MCMC_CandidateMove_Index(:,PC)
                      !!! if we're not in pair proposal mode we did not yet check if
                      !!! A's reverse particle is floating (else we did already):
                   ENDIF

                   !!! the first condition is needed when not using pair proposal mode
                   IF (MCMC_Particle_Ab_IsFloating(PC)) THEN
                      ReverseFloatingP=MCMCParticle(MCMC_LabelsBeforeJump_A(PC),MCMCproposal(index_,labels_3d))
                   ENDIF

                   !!! in pair proposal, if A' is floating, B' is as well (they are the
                   !!! same particle):
                   IF (MCMCusePairProposal.AND.MCMC_Particle_Bb_IsFloating(PC)) THEN
                      !!! only possible in pair proposal mode
                      ReverseFloatingP=MCMCParticle(MCMC_LabelsBeforeJump_B(PC),MCMCproposal(index_,labels_3d))
                   ENDIF

                   !!! finally convert the regular particle into a floating particle,
                   !!! i.e. insert it in the floating DS and remove it from the regular:
                   IF (MCMC_Particle_Ab_IsFloating(PC).OR.MCMC_Particle_Bb_IsFloating(PC)) THEN
                      !!! insert the reverse particle in the appropriate container. If
                      !!! there is no space, we reject the move.
                      IF (.NOT.MCMCInsertFloatingParticle(ipatch,index_(1),index_(2),index_(3),ReverseFloatingP,.TRUE.)) THEN
                         !!! TODO: calling MCMCReject from here invalidates the result.
                         !!! MCMCReject(&vAppliedParticles,&vAppliedParticleOrigLabels);
                         MCMCHardReject=.TRUE.
                      ENDIF
                   ENDIF
                ENDDO !PC=1,MCMCstepsize

                !!! We are now in the state x'.
                !!! Calculate Q'(A) and maybe Q'(B). Note that this has to be done after
                !!! all particles were applied.
                DO PC=1,MCMCstepsize
                   !!! Calculate MCMC_qb_A and MCMC_qb_B
                   IF (.NOT.MCMCuseBiasedProposal) THEN
                      MCMC_qb_A(PC)=one
                      MCMC_qb_B(PC)=one
                   ELSE
                      index_A=MCMC_CandidateMove_Index(:,PC)
                      MCMC_qb_A(PC)=MCMCproposal(index_A,labels_3d)
                      IF (MCMCusePairProposal.AND..NOT.MCMC_SingleParticleMoveForPairProposals) THEN
                         index_B=MCMC_PartnerMove_Index(:,PC)
                         MCMC_qb_B(PC)=MCMCproposal(index_B,labels_3d)
                      ENDIF
                   ENDIF

                   !!! Normalize MCMC_qb_A and MCMC_qb_B
                   IF (MCMC_Particle_Ab_IsFloating(PC)) THEN
                      Normalizer_qb_A=REAL(MCMCFloatingParticlesProposalNormalizer,MK)
                   ELSE
                      Normalizer_qb_A=MCMCGetProposalNormalizer(MCMC_CandidateMove(PC)%candlabel,MCMC_LabelsBeforeJump_A(PC))
                   ENDIF
                   MCMC_qb_A(PC)=MCMC_qb_A(PC)/Normalizer_qb_A

                   IF (MCMCusePairProposal.AND..NOT.MCMC_SingleParticleMoveForPairProposals) THEN
                      IF (MCMC_Particle_Bb_IsFloating(PC)) THEN
                         Normalizer_qb_B=REAL(MCMCFloatingParticlesProposalNormalizer,MK)
                      ELSE
                         Normalizer_qb_B=MCMCGetProposalNormalizer(MCMC_PartnerMove(PC)%candlabel,MCMC_LabelsBeforeJump_B(PC))
                      ENDIF
                      MCMC_qb_B(PC)=MCMC_qb_B(PC)/Normalizer_qb_B
                   ENDIF
                   !!! Finally, we omit half of the calculations if particle A == B
                   IF (MCMC_SingleParticleMoveForPairProposals) THEN
                      MCMC_q_B(PC)   =MCMC_q_A(PC)
                      MCMC_q_A_B(PC) =MCMC_q_B_A(PC);
                      MCMC_qb_B_A(PC)=MCMC_qb_A_B(PC);
                      MCMC_qb_B(PC)  =MCMC_qb_A(PC);
                   ENDIF
                ENDDO !PC=1,MCMCstepsize

                sbpitr => mesh%subpatch%next()
                ipatch=ipatch+1
             ENDDO
             NULLIFY(labels_3d,image_3d)
          END SELECT

          MCMCMove=.TRUE.
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCMove



      END MODULE ppm_rc_module_mcmc


