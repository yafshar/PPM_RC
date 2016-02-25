      !-------------------------------------------------------------------------
      !  Module       :                    ppm_rc_module_linkedlist
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
      !  Module       :                    ppm_rc_module_linkedlist
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains subroutines for VTK outputting
      !
      !  Remarks      :
      !
      !  References   :
      !-------------------------------------------------------------------------
      MODULE ppm_rc_module_linkedlist
        !----------------------------------------------------------------------
        !  Modules
        !----------------------------------------------------------------------
        USE ppm_module_alloc
        USE ppm_module_data
        USE ppm_module_error
        USE ppm_module_write
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_container_typedef
        IMPLICIT NONE

        PRIVATE
        !----------------------------------------------------------------------
        !
        !----------------------------------------------------------------------
#include "./linkedlist/ppm_rc_linkedlistypedef.f"

        !----------------------------------------------------------------------
        !  Some work memory on the heap
        !----------------------------------------------------------------------
        TYPE(ppm_rc_c_list), DIMENSION(:), ALLOCATABLE :: ppm_rc_seeds
        TYPE(ppm_rc_c_list), DIMENSION(:), ALLOCATABLE :: ppm_rc_seeds_to_remove
        TYPE(ppm_rc_c_list), DIMENSION(:), ALLOCATABLE :: InnerContourContainer
        TYPE(ppm_rc_c_list), DIMENSION(:), ALLOCATABLE :: Candidates
        TYPE(ppm_rc_c_list), DIMENSION(:), ALLOCATABLE :: CompetingRegions
        TYPE(ppm_rc_c_list), DIMENSION(:), ALLOCATABLE :: m_Seeds

        TYPE(ppm_rc_c_list), DIMENSION(:), ALLOCATABLE :: MCMCParticleInContainerHistory
        TYPE(ppm_rc_c_list), DIMENSION(:), ALLOCATABLE :: MCMCFloatingParticleInContainerHistory
        TYPE(ppm_rc_c_list), DIMENSION(:), ALLOCATABLE :: MCMCLabelImageHistory

        TYPE(ppm_rc_list),   DIMENSION(:), ALLOCATABLE :: MCMCAppliedParticleOrigLabels

        INTEGER,             DIMENSION(:), ALLOCATABLE :: Candidates_list

        !----------------------------------------------------------------------
        ! Public
        !----------------------------------------------------------------------
        PUBLIC :: ppm_rc_link,   ppm_rc_link_,   ppm_rc_link_2
        PUBLIC :: ppm_rc_list,   ppm_rc_list_,   ppm_rc_list_2
        PUBLIC :: ppm_rc_c_list, ppm_rc_c_list_, ppm_rc_c_list_2

        PUBLIC :: ppm_rc_seeds
        PUBLIC :: ppm_rc_seeds_to_remove
        PUBLIC :: InnerContourContainer
        PUBLIC :: Candidates
        PUBLIC :: CompetingRegions
        PUBLIC :: m_Seeds
        PUBLIC :: Candidates_list

        PUBLIC :: MCMCParticleInContainerHistory
        PUBLIC :: MCMCFloatingParticleInContainerHistory
        PUBLIC :: MCMCLabelImageHistory
        PUBLIC :: MCMCAppliedParticleOrigLabels

      CONTAINS

#include "./linkedlist/ppm_rc_linkedlistypeproc.f"

      END MODULE ppm_rc_module_linkedlist
