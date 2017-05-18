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
      !  Please do cite:
      !
      !  Y. Afshar, and I. F. Sbalzarini. A Parallel Distributed-Memory Particle
      !  Method Enables Acquisition-Rate Segmentation of Large Fluorescence
      !  Microscopy Images. PLoS ONE 11(4):e0152528, (2016).
      !
      !  when publishing research data obtained using PPM_RC
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

        TYPE(ppm_rc_c_list)                            :: MCMCParticleInContainerHistory
        TYPE(ppm_rc_c_list)                            :: MCMCFloatingParticleInContainerHistory
        TYPE(ppm_rc_c_list)                            :: MCMCLabelImageHistory

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

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------

        INTERFACE IndexHashFunctor_2d
          MODULE PROCEDURE IndexHashFunctor32_2d
          MODULE PROCEDURE IndexHashFunctor32__2d
        END INTERFACE
        INTERFACE IndexHashFunctor_3d
          MODULE PROCEDURE IndexHashFunctor32_3d
          MODULE PROCEDURE IndexHashFunctor32__3d
        END INTERFACE

        !!! I have two versions both for 32 and 64 bits
        !!! as the FORTRAN compiler can not create one interface
        !!! for different output types

        INTERFACE HashIndexFunctor_2d
          MODULE PROCEDURE HashIndexFunctor32_2d
          MODULE PROCEDURE HashIndexFunctor32__2d
          MODULE PROCEDURE HashIndexFunctor64_2d
          MODULE PROCEDURE HashIndexFunctor64__2d
        END INTERFACE
        INTERFACE HashIndexFunctor_3d
          MODULE PROCEDURE HashIndexFunctor32_3d
          MODULE PROCEDURE HashIndexFunctor32__3d
          MODULE PROCEDURE HashIndexFunctor64_3d
          MODULE PROCEDURE HashIndexFunctor64__3d
        END INTERFACE

      CONTAINS
#define __2D  2
#define __3D  3

#define  DTYPE(a) a/**/_3d
#define __DIME  __3D
#include "./linkedlist/ppm_rc_linkedlistypeproc.f"
#undef  __DIME
#undef  DTYPE

#define  DTYPE(a) a/**/_2d
#define __DIME  __2D
#include "./linkedlist/ppm_rc_linkedlistypeproc.f"
#undef  __DIME
#undef  DTYPE
#undef  __2D
#undef  __3D

      END MODULE ppm_rc_module_linkedlist
