      !-------------------------------------------------------------------------
      !  Module       :                 ppm_rc_module_linkedlist
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains subroutines for VTK outputting
      !
      !  Remarks      :
      !
      !  References   :
      !-------------------------------------------------------------------------
      !  MOSAIC Group
      !  Max Planck Institute of Molecular Cell Biology and Genetics
      !  Pfotenhauerstr. 108, 01307 Dresden, Germany
      !
      !  Author           - y.afshar           June   2014
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
