      !-------------------------------------------------------------------------
      !  Module       :                    ppm_rc_module_move
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
      !  Author           - y.afshar           June   2014
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Module       :                    ppm_rc_module_move
      !-------------------------------------------------------------------------
      !
      !  Purpose      :  Module contains read subroutines.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      MODULE ppm_rc_module_move
        !----------------------------------------------------------------------
        !  Modules
        !----------------------------------------------------------------------
        USE ppm_module_inl_hash

        USE ppm_rc_module_global
        USE ppm_rc_module_linkedlist
        USE ppm_rc_module_fire
        USE ppm_rc_module_topologicalnumber
        USE ppm_rc_module_bgl
        USE ppm_rc_module_energy
        IMPLICIT NONE

        PRIVATE
        !----------------------------------------------------------------------
        !  Some work memory on the heap
        !----------------------------------------------------------------------
        TYPE(ppm_htable)                      :: htableg
        TYPE(ppm_rc_list),        POINTER     :: tmp_Candidates => NULL()

        REAL(MK), DIMENSION(:),   ALLOCATABLE :: energyt

        INTEGER, DIMENSION(:),    ALLOCATABLE :: sndbuf
        INTEGER, DIMENSION(:),    ALLOCATABLE :: rcvbuf
        INTEGER, DIMENSION(:),    ALLOCATABLE :: sndrcvbuf

        INTEGER, DIMENSION(:),    POINTER     :: energyrank => NULL()
        INTEGER, DIMENSION(:),    ALLOCATABLE :: bufi
        INTEGER, DIMENSION(:),    ALLOCATABLE :: bufa
        INTEGER, DIMENSION(:),    ALLOCATABLE :: buft
        INTEGER, DIMENSION(:),    ALLOCATABLE :: ltg
        INTEGER, DIMENSION(:),    ALLOCATABLE :: labelt
        INTEGER, DIMENSION(:),    ALLOCATABLE :: candlabelt
        INTEGER, DIMENSION(:),    ALLOCATABLE :: ccandlabelt
        INTEGER, DIMENSION(:),    ALLOCATABLE :: ndaughterst
        INTEGER, DIMENSION(:),    ALLOCATABLE :: nmotherst
        INTEGER, DIMENSION(:,:),  ALLOCATABLE :: daughterst
        INTEGER, DIMENSION(:),    ALLOCATABLE :: acceptedt
        INTEGER, DIMENSION(:),    ALLOCATABLE :: vSortedList
        INTEGER                               :: totNpart

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE ppm_rc_ChangeContourParticleLabelToCandidateLabel_2d
           MODULE PROCEDURE ppm_rc_ChangeContourParticleLabelToCandidateLabel_2d
           MODULE PROCEDURE ppm_rc_ChangeGhostContourParticleLabelToCandidateLabel_2d
        END INTERFACE
        INTERFACE ppm_rc_ChangeContourParticleLabelToCandidateLabel_3d
           MODULE PROCEDURE ppm_rc_ChangeContourParticleLabelToCandidateLabel_3d
           MODULE PROCEDURE ppm_rc_ChangeGhostContourParticleLabelToCandidateLabel_3d
        END INTERFACE
        !----------------------------------------------------------------------
        ! Public
        !----------------------------------------------------------------------
        PUBLIC :: ppm_rc_contour_propagation_2d
        PUBLIC :: ppm_rc_contour_propagation_3d
        PUBLIC :: ppm_rc_FilterCandidates_2d
        PUBLIC :: ppm_rc_FilterCandidates_3d
        PUBLIC :: ppm_rc_FilterCandidatesContainerUsingRanks_2d
        PUBLIC :: ppm_rc_FilterCandidatesContainerUsingRanks_3d
        PUBLIC :: ppm_rc_DetectOscillations
        PUBLIC :: ppm_rc_DetectOscillations2
        PUBLIC :: ppm_rc_ChangeContourParticleLabelToCandidateLabel_2d
        PUBLIC :: ppm_rc_ChangeContourParticleLabelToCandidateLabel_3d
        PUBLIC :: ppm_rc_move_2d
        PUBLIC :: ppm_rc_move_3d

      CONTAINS

#define __2D  2
#define __3D  3

#define  DTYPE(a) a/**/_3d
#define __DIME  __3D
#include "./move/ppm_rc_contour_propagation.f"
#include "./move/ppm_rc_FilterCandidates.f"
#include "./move/ppm_rc_FilterCandidatesContainerUsingRanks.f"
#include "./move/ppm_rc_oscillation.f"
#include "./move/ppm_rc_ChangeContourParticleLabelToCandidateLabel.f"
#include "./move/ppm_rc_ChangeGhostContourParticleLabelToCandidateLabel.f"
#include "./move/ppm_rc_move.f"
#undef  __DIME
#undef  DTYPE

#define  DTYPE(a) a/**/_2d
#define __DIME  __2D
#include "./move/ppm_rc_contour_propagation.f"
#include "./move/ppm_rc_FilterCandidates.f"
#include "./move/ppm_rc_FilterCandidatesContainerUsingRanks.f"
#include "./move/ppm_rc_ChangeContourParticleLabelToCandidateLabel.f"
#include "./move/ppm_rc_ChangeGhostContourParticleLabelToCandidateLabel.f"
#include "./move/ppm_rc_move.f"
#undef  __DIME
#undef  DTYPE

#undef  __2D
#undef  __3D

      END MODULE ppm_rc_module_move
