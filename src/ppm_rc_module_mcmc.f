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
      !  Please do cite:
      !
      !  Y. Afshar, and I. F. Sbalzarini. A Parallel Distributed-Memory Particle
      !  Method Enables Acquisition-Rate Segmentation of Large Fluorescence
      !  Microscopy Images. PLoS ONE 11(4):e0152528, (2016).
      !
      !  when publishing research data obtained using PPM_RC
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

        USE ppm_rc_module_global, ONLY : MK,zero,zerod,MCMCParticle,MCMCHistoryParticle
        USE ppm_rc_module_energy, ONLY : e_data,e_length
        USE ppm_rc_module_hash, ONLY : ppm_rc_htable,ppm_rc_MCMCParticlehtable, &
        &   ppm_rc_MCMCResultshtable
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
        TYPE(ppm_htable)                                                   :: MCMChtable
        !!! hash table for the region labels

        TYPE(ppm_rc_MCMCParticlehtable), DIMENSION(:), ALLOCATABLE, TARGET :: MCMCRegularParticles
        TYPE(ppm_rc_MCMCParticlehtable),                            TARGET :: MCMCFloatingParticles

        TYPE(ppm_rc_MCMCParticlehtable), DIMENSION(:), ALLOCATABLE, TARGET :: MCMCRegularParticlesInCell
        TYPE(ppm_rc_MCMCParticlehtable), DIMENSION(:), ALLOCATABLE, TARGET :: MCMCFloatingParticlesInCell

        TYPE(ppm_rc_MCMCResultshtable)                                     :: MCMCResults

        REAL(MK)                                                           :: MCMCsampleOffBoundaryPercentage
        REAL(ppm_kind_double)                                              :: MCMCProbabilityToProposeAFloatingParticle

        REAL(ppm_kind_double)                                              :: MCMCZe
        !!! Normalizing constants of the edge image
        REAL(MK),                        DIMENSION(:), ALLOCATABLE         :: MCMClengthProposalMask
        !!! For a fast proposal computation



        REAL(ppm_kind_double)                                              :: MCMCFloatingParticlesProposalNormalizer
        REAL(ppm_kind_double)                                              :: MCMCFloatingParticlesProposalNormalizerlocal

        REAL(MK),                        DIMENSION(:), ALLOCATABLE         :: MCMCRegularParticlesProposalNormalizer
        REAL(MK),                        DIMENSION(:), ALLOCATABLE         :: MCMCRegularParticlesProposalNormalizerlocal
        !!! map containg a label and its float Normalize factor

        REAL(ppm_kind_double)                                              :: MCMCTotalNormalizer
        REAL(ppm_kind_double)                                              :: MCMCTotalNormalizerlocal











        REAL(MK),                        DIMENSION(:), ALLOCATABLE         :: MCMCcellsize
        !!! MCMC Cell size which is using for partitioning the simulation domain
        !!! Into equal size cells
        INTEGER,                         DIMENSION(:), ALLOCATABLE         :: MCMCcellNm
        !!! Number of cells in x,y (& z) direction
        INTEGER,                         DIMENSION(:), ALLOCATABLE         :: MCMCboundarycellIndex
        !!! An array to store all the boundary cell indices
        INTEGER                                                            :: MCMCboundarycellSize
        !!! Size of the MCMCboundarycellIndex array

        !!! Number of cells on each domain which are used for ghost communication

        INTEGER,                         DIMENSION(:), ALLOCATABLE         :: MCMCinteriorcellIndex
        INTEGER,                         DIMENSION(:), ALLOCATABLE         :: MCMCinteriorcellDisp

        TYPE(ppm_rc_MCMCParticlehtable),               POINTER             :: MCMCActiveCandidates
        !!! Temporary pointer to an object



        INTEGER                                                            :: MCMCStatsUpdateModulus
        INTEGER,                         DIMENSION(:), ALLOCATABLE         :: MCMCVisitedLabels
        !!! Ordered list of visited labels, it is good when you wanna loop through
        !!! Sweep through the label image. If a label didn't exist, register it
        !!! (add an entry in the regionsLabel vector and initialize normalizer
        !!! statistics). And always initialize these for the BG label (=0).
        INTEGER,                         DIMENSION(:), ALLOCATABLE         :: MCMCRegionLabel
        !!! vector of labels
        INTEGER                                                            :: MCMCRegionLabelSize



        INTEGER,                         DIMENSION(:), ALLOCATABLE         :: MCMColdghost
        !!! Temporary array for keeping one layer ghost information for later
        !!! comparison with updated ghost and change the particle label accordingly
        !!! This is a provate variable only accessible in MCMC module




        LOGICAL                                                            :: MCMCVerbose
        LOGICAL                                                            :: MCMCUseForbiddenRegion
        LOGICAL                                                            :: MCMCuseSafeSampling
        !!!

        CHARACTER(LEN=ppm_char)                                            :: MCMCmarginalFileNamePrefix

        !----------------------------------------------------------------------
        !  PUBLIC variables
        !----------------------------------------------------------------------
        PUBLIC :: MCMChtable

        PUBLIC :: MCMCRegularParticles
        PUBLIC :: MCMCRegularParticlesInCell
        PUBLIC :: MCMCFloatingParticles
        PUBLIC :: MCMCFloatingParticlesInCell

        PUBLIC :: MCMCsampleOffBoundaryPercentage
        PUBLIC :: MCMCProbabilityToProposeAFloatingParticle
        PUBLIC :: MCMCZe

        PUBLIC :: MCMClengthProposalMask


        PUBLIC :: MCMCRegularParticlesProposalNormalizer
        PUBLIC :: MCMCRegularParticlesProposalNormalizerlocal

        PUBLIC :: MCMCFloatingParticlesProposalNormalizer
        PUBLIC :: MCMCFloatingParticlesProposalNormalizerlocal


        PUBLIC :: MCMCTotalNormalizer
        PUBLIC :: MCMCTotalNormalizerlocal

        PUBLIC :: MCMCActiveCandidates

        PUBLIC :: MCMCStatsUpdateModulus
        PUBLIC :: MCMCVisitedLabels

        PUBLIC :: MCMCRegionLabel
        PUBLIC :: MCMCRegionLabelSize

        PUBLIC :: MCMCVerbose
        PUBLIC :: MCMCUseForbiddenRegion
        PUBLIC :: MCMCuseSafeSampling

        PUBLIC :: MCMCmarginalFileNamePrefix


        PUBLIC :: MCMCcellsize
        PUBLIC :: MCMCcellNm
        PUBLIC :: MCMCboundarycellIndex
        PUBLIC :: MCMCboundarycellSize
        PUBLIC :: MCMCinteriorcellIndex
        PUBLIC :: MCMCinteriorcellDisp


        !----------------------------------------------------------------------
        !  Some work memory on the heap
        !----------------------------------------------------------------------
!         REAL(MK), DIMENSION(:), ALLOCATABLE :: tmp1_r

        INTEGER,  DIMENSION(:), ALLOCATABLE :: tmp1_i

        TYPE(MCMCParticle), DIMENSION(:), ALLOCATABLE :: MCMC_CandidateMove
        TYPE(MCMCParticle), DIMENSION(:), ALLOCATABLE :: MCMC_PartnerMove

        REAL(ppm_kind_double), DIMENSION(:), ALLOCATABLE :: MCMC_q_A
        REAL(ppm_kind_double), DIMENSION(:), ALLOCATABLE :: MCMC_q_B
        REAL(ppm_kind_double), DIMENSION(:), ALLOCATABLE :: MCMC_q_A_B
        REAL(ppm_kind_double), DIMENSION(:), ALLOCATABLE :: MCMC_q_B_A
        !!! Forward probability
        REAL(ppm_kind_double), DIMENSION(:), ALLOCATABLE :: MCMC_qb_A
        REAL(ppm_kind_double), DIMENSION(:), ALLOCATABLE :: MCMC_qb_B
        REAL(ppm_kind_double), DIMENSION(:), ALLOCATABLE :: MCMC_qb_A_B
        REAL(ppm_kind_double), DIMENSION(:), ALLOCATABLE :: MCMC_qb_B_A
        !!! Backward probability

        INTEGER, DIMENSION(:),   ALLOCATABLE :: MCMC_LabelsBeforeJump_A
        INTEGER, DIMENSION(:),   ALLOCATABLE :: MCMC_LabelsBeforeJump_B
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: MCMC_CandidateMove_Index
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: MCMC_PartnerMove_Index

        LOGICAL, DIMENSION(:), ALLOCATABLE :: MCMC_Particle_A_IsFloating
        LOGICAL, DIMENSION(:), ALLOCATABLE :: MCMC_Particle_Ab_IsFloating
        LOGICAL, DIMENSION(:), ALLOCATABLE :: MCMC_Particle_Bb_IsFloating

        LOGICAL                            :: MCMC_SingleParticleMoveForPairProposals

        !----------------------------------------------------------------------
        !  PUBLIC variables
        !----------------------------------------------------------------------
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


!

!
!         INTERFACE MCMCApplyParticle
!           MODULE PROCEDURE MCMCApplyParticle_2d
!           MODULE PROCEDURE MCMCApplyParticle_3d
!         END INTERFACE

        !----------------------------------------------------------------------
        ! Public
        !----------------------------------------------------------------------
        PUBLIC :: MCMCCheckParameters
        PUBLIC :: MCMCFreeMemory
        PUBLIC :: CreateMCMClengthProposalMask
        PUBLIC :: DestroyMCMClengthProposalMask

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE MCMCInsertLabelInRegionLabel
          MODULE PROCEDURE MCMCInsertBGLabelInRegionLabel
          MODULE PROCEDURE MCMCInsertLabelInRegionLabel
        END INTERFACE

        PUBLIC :: MCMCInsertLabelInRegionLabel
        PUBLIC :: MCMCUpdateRegionLabel

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE MCMCproposal
          MODULE PROCEDURE MCMCproposal_2d
          MODULE PROCEDURE MCMCproposal__2d
          MODULE PROCEDURE MCMCproposal_3d
          MODULE PROCEDURE MCMCproposal__3d
        END INTERFACE

        PUBLIC :: MCMCproposal

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE MCMCIsParticleTopoValid
          MODULE PROCEDURE MCMCIsParticleTopoValid_2d
          MODULE PROCEDURE MCMCIsParticleTopoValid__2d
          MODULE PROCEDURE MCMCIsParticleTopoValid_3d
          MODULE PROCEDURE MCMCIsParticleTopoValid__3d
        END INTERFACE

        PUBLIC :: MCMCIsParticleTopoValid

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE MCMCgetRegularParticlesAtIndex
          MODULE PROCEDURE MCMCgetRegularParticlesAtIndex_2d
          MODULE PROCEDURE MCMCgetRegularParticlesAtIndex__2d
          MODULE PROCEDURE MCMCgetRegularParticlesAtIndex_3d
          MODULE PROCEDURE MCMCgetRegularParticlesAtIndex__3d
        END INTERFACE

        PUBLIC :: MCMCgetRegularParticlesAtIndex

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE MCMCgetParticlesInFGNeighborhood
          MODULE PROCEDURE MCMCgetParticlesInFGNeighborhood_2d
          MODULE PROCEDURE MCMCgetParticlesInFGNeighborhood_3d
        END INTERFACE

        PUBLIC :: MCMCgetParticlesInFGNeighborhood

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE MCMCgetParticlesInBGNeighborhood
          MODULE PROCEDURE MCMCgetParticlesInBGNeighborhood_2d
          MODULE PROCEDURE MCMCgetParticlesInBGNeighborhood_3d
        END INTERFACE

        PUBLIC :: MCMCgetParticlesInBGNeighborhood

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE MCMCInsertCandidatesToContainers
          MODULE PROCEDURE MCMCInsertCandidatesToContainers_2d
          MODULE PROCEDURE MCMCInsertCandidatesToContainers_3d
        END INTERFACE

        PUBLIC :: MCMCInsertCandidatesToContainers

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE MCMCEraseCandidatesFromContainers
          MODULE PROCEDURE MCMCEraseCandidatesFromContainers_2d
          MODULE PROCEDURE MCMCEraseCandidatesFromContainers_3d
        END INTERFACE

        PUBLIC :: MCMCEraseCandidatesFromContainers

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE MCMCInsertFloatingParticle
          MODULE PROCEDURE MCMCInsertFloatingParticle_2d
          MODULE PROCEDURE MCMCInsertFloatingParticle_3d
        END INTERFACE

        PUBLIC :: MCMCInsertFloatingParticle

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE MCMCInsertFloatingParticleCumulative
          MODULE PROCEDURE MCMCInsertFloatingParticleCumulative_2d
          MODULE PROCEDURE MCMCInsertFloatingParticleCumulative_3d
        END INTERFACE

        PUBLIC :: MCMCInsertFloatingParticleCumulative

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE MCMCEraseFloatingParticle
          MODULE PROCEDURE MCMCEraseFloatingParticle_2d
          !!! first find out whether there is a particle at coords or not
          MODULE PROCEDURE MCMCEraseFloatingParticle__2d
          !!! give the searched tmpParticle at coords as an input
          MODULE PROCEDURE MCMCEraseFloatingParticle_3d
          MODULE PROCEDURE MCMCEraseFloatingParticle__3d
        END INTERFACE

        PUBLIC :: MCMCEraseFloatingParticle

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE MCMCAddAndRemoveParticlesWhenMove
          MODULE PROCEDURE MCMCAddAndRemoveParticlesWhenMove_2d
          MODULE PROCEDURE MCMCAddAndRemoveParticlesWhenMove_3d
        END INTERFACE

        PUBLIC :: MCMCAddAndRemoveParticlesWhenMove

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE MCMCupdateProposalsAndFilterTopologyInNeighborhood
          MODULE PROCEDURE MCMCupdateProposalsAndFilterTopologyInNeighborhood_2d
          MODULE PROCEDURE MCMCupdateProposalsAndFilterTopologyInNeighborhood_3d
        END INTERFACE

        PUBLIC :: MCMCupdateProposalsAndFilterTopologyInNeighborhood

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE MCMCIsRegularParticle
          MODULE PROCEDURE MCMCIsRegularParticle_2d
          MODULE PROCEDURE MCMCIsRegularParticle_3d
        END INTERFACE

        PUBLIC :: MCMCIsRegularParticle

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE MCMCParticleHasFloatingProperty
          MODULE PROCEDURE MCMCParticleHasFloatingProperty_2d
          MODULE PROCEDURE MCMCParticleHasFloatingProperty__2d
          MODULE PROCEDURE MCMCParticleHasFloatingProperty___2d
          MODULE PROCEDURE MCMCParticleHasFloatingProperty____2d
          MODULE PROCEDURE MCMCParticleHasFloatingProperty_3d
          MODULE PROCEDURE MCMCParticleHasFloatingProperty__3d
          MODULE PROCEDURE MCMCParticleHasFloatingProperty___3d
          MODULE PROCEDURE MCMCParticleHasFloatingProperty____3d
        END INTERFACE

        PUBLIC :: MCMCParticleHasFloatingProperty


        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE MCMCAddAndRemoveParticlesFromGhostMove
          MODULE PROCEDURE MCMCAddAndRemoveParticlesFromGhostMove_2d
          MODULE PROCEDURE MCMCAddAndRemoveParticlesFromGhostMove_3d
        END INTERFACE

        PUBLIC :: MCMCAddAndRemoveParticlesFromGhostMove


        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE MCMCReject
          MODULE PROCEDURE MCMCReject_2d
          MODULE PROCEDURE MCMCReject_3d
        END INTERFACE

        PUBLIC :: MCMCReject




!         PUBLIC :: MCMCOffBoundarySample_2d
!         PUBLIC :: MCMCOffBoundarySample_3d
!         PUBLIC :: MCMCUpdateAllParticlesFwdProposals
!         PUBLIC :: MCMCFindParticleA
!         PUBLIC :: MCMCApplyParticle
!         PUBLIC :: MCMCFindParticleB
!         PUBLIC :: MCMCMove

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
!         INTERFACE MCMCGetIndexFromEdgeDensity
!           MODULE PROCEDURE MCMCGetIndexFromEdgeDensity_2d
!           MODULE PROCEDURE MCMCGetIndexFromEdgeDensity_3d
!         END INTERFACE

!         PUBLIC :: MCMCGetIndexFromEdgeDensity

        PUBLIC :: ppm_rc_MCMCinitmove
        PUBLIC :: ppm_rc_mcmc_move_2d


      CONTAINS

#include "./mcmc/ppm_rc_mcmctypeproc.f"

#define __2D  2
#define __3D  3

#define  DTYPE(a) a/**/_3d
#define __DIME  __3D
#include "./mcmc/ppm_rc_mcmc_move.f"
#undef  __DIME
#undef  DTYPE

#define  DTYPE(a) a/**/_2d
#define __DIME  __2D
#include "./mcmc/ppm_rc_mcmc_move.f"
#undef  __DIME
#undef  DTYPE

#undef  __2D
#undef  __3D


!
! !           ! IF (.NOT.Growth) THEN
! !           !    !!! we check if we want to delete a floating particle.
! !           !    IF (NumberOfFloatingParticles == 0) THEN
! !           !       MCMCOffBoundarySample=.FALSE.
! !           !       RETURN
! !           !    ENDIF
! !           !    ParticleIndex = ppm_rc_Saru_IPRNG(NumberOfFloatingParticles)
! !           !    Particle = MCMCFloatingParticles(ParticleIndex)
! !           !    EdgeWeight = m_EdgeImage->GetPixel(vParticle.m_Index)
! !           !    !!! the backward probability:
! !           !    Pinsertion = EdgeWeight / (MCMCZe + EdgeWeight)
! !           !    !!! the forward probability:
! !           !    Pdeletion = one / NumberOfFloatingParticles
! !           ! ELSE
! !           !    !!! sample a location
! !           !    Particle.m_Index = MCMCGetIndexFromEdgeDensity()
! !           !    !!! the forward probability:
! !           !    EdgeWeight = m_EdgeImage->GetPixel(vParticle.m_Index)
! !           !    Pinsertion = EdgeWeight / MCMCZe
! !           !    !!! the backward probability:
! !           !    Pdeletion = one/(vNumberOfFloatingParticles+1)
! !           ! ENDIF
!








      END MODULE ppm_rc_module_mcmc


