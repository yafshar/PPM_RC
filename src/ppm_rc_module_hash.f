      !-------------------------------------------------------------------------
      !  Module       :                    ppm_rc_module_hash
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
      !  Module       :                    ppm_rc_module_hash
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Module contains hash subroutines.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      MODULE ppm_rc_module_hash
        !----------------------------------------------------------------------
        !  Modules
        !----------------------------------------------------------------------
#ifdef __F2003
        USE ISO_C_BINDING
#endif

        USE ppm_rc_module_global
        USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list_2
        IMPLICIT NONE

        PRIVATE

        INTEGER(ppm_kind_int64), PARAMETER :: htable_null_li = -HUGE(1_ppm_kind_int64)
        !!! NULL value for hash table
        INTEGER(ppm_kind_int64), PARAMETER :: seed1 = 738235926_ppm_kind_int64
        !!! Hardcoded seed value taken from MurmurHash
        INTEGER(ppm_kind_int64), PARAMETER :: seed2 = 1243832038_ppm_kind_int64
        !!! Hardcoded seed value taken from MurmurHash

        TYPE ppm_rc_htable
          !---------------------------------------------------------------------
          !  Declaration of arrays
          !---------------------------------------------------------------------
          INTEGER(ppm_kind_int64), DIMENSION(:), ALLOCATABLE :: keys
          !!! Array for keeping hash table keys.
          REAL(MK),                DIMENSION(:), ALLOCATABLE :: borders_pos
          !!! Array for keeping positions of cells on "borders" array.
          !---------------------------------------------------------------------
          !  Declaration of variables
          !--------------------------------------------------------------------
          INTEGER(ppm_kind_int64)                            :: nrow = 0_ppm_kind_int64
          !!! number of rows in hash table
         CONTAINS
          PROCEDURE :: create => create_htable
          PROCEDURE :: destroy => destroy_htable
          PROCEDURE :: h_key
          PROCEDURE :: hash_insert
          PROCEDURE :: hash_insert_32
          PROCEDURE :: hash_insert_2D3D
          PROCEDURE :: hash_insert_2d
          PROCEDURE :: hash_insert_3d
          GENERIC   :: insert =>         &
          &            hash_insert,      &
          &            hash_insert_32,   &
          &            hash_insert_2D3D, &
          &            hash_insert_2d,   &
          &            hash_insert_3d
          PROCEDURE :: hash_search
          PROCEDURE :: hash_search_32
          PROCEDURE :: hash_search_2d
          PROCEDURE :: hash_search_3d
          GENERIC   :: search =>       &
          &            hash_search,    &
          &            hash_search_32, &
          &            hash_search_2d, &
          &            hash_search_3d
          PROCEDURE :: hash_remove
          PROCEDURE :: hash_remove_32
          PROCEDURE :: hash_remove_2D3D
          PROCEDURE :: hash_remove_2d
          PROCEDURE :: hash_remove_3d
          GENERIC   :: remove =>         &
          &            hash_remove,      &
          &            hash_remove_32,   &
          &            hash_remove_2D3D, &
          &            hash_remove_2d,   &
          &            hash_remove_3d
          PROCEDURE :: grow => grow_htable
          PROCEDURE :: shrink => shrink_htable
          PROCEDURE :: size => hash_size
        END TYPE ppm_rc_htable
        !!! This hashtable is for hashing REAL value using the same PPM implemenattion

        TYPE ppm_rc_MCMCParticlehtable
          !---------------------------------------------------------------------
          !  Declaration of arrays
          !---------------------------------------------------------------------
          INTEGER(ppm_kind_int64), DIMENSION(:), ALLOCATABLE :: keys
          !!! Array for keeping hash table keys.
          TYPE(MCMCParticle),      DIMENSION(:), ALLOCATABLE :: borders_pos
          !!! Array for keeping positions of cells on "borders" array.
          !---------------------------------------------------------------------
          !  Declaration of variables
          !--------------------------------------------------------------------
          INTEGER(ppm_kind_int64)                            :: nrow = 0_ppm_kind_int64
          !!! number of rows in hash table
         CONTAINS
          PROCEDURE :: create  => MCMCParticle_create_htable
          PROCEDURE :: destroy => MCMCParticle_destroy_htable
          PROCEDURE :: h_key   => MCMCParticle_h_key
          PROCEDURE :: MCMCParticle_hash_insert
          PROCEDURE :: MCMCParticle_hash_insert_2D3D
          PROCEDURE :: MCMCParticle_hash_insert_2d
          PROCEDURE :: MCMCParticle_hash_insert_3d
          GENERIC   :: insert =>                      &
          &            MCMCParticle_hash_insert,      &
          &            MCMCParticle_hash_insert_2D3D, &
          &            MCMCParticle_hash_insert_2d,   &
          &            MCMCParticle_hash_insert_3d
          PROCEDURE :: MCMCParticle_hash_search
          PROCEDURE :: MCMCParticle_hash_search_2d
          PROCEDURE :: MCMCParticle_hash_search_3d
          GENERIC   :: search =>                    &
          &            MCMCParticle_hash_search,    &
          &            MCMCParticle_hash_search_2d, &
          &            MCMCParticle_hash_search_3d
          PROCEDURE :: MCMCParticle_hash_remove
          PROCEDURE :: MCMCParticle_hash_remove_2D3D
          PROCEDURE :: MCMCParticle_hash_remove_2d
          PROCEDURE :: MCMCParticle_hash_remove_3d
          GENERIC   :: remove =>                      &
          &            MCMCParticle_hash_remove,      &
          &            MCMCParticle_hash_remove_2D3D, &
          &            MCMCParticle_hash_remove_2d,   &
          &            MCMCParticle_hash_remove_3d
          PROCEDURE :: grow => MCMCParticle_grow_htable
          PROCEDURE :: shrink => MCMCParticle_shrink_htable
          PROCEDURE :: size => MCMCParticle_hash_size
          !!! check the elements in the hash table and pack the proposal
          !!! in an allocatable array
          PROCEDURE :: MCMCParticle_hash_packproposal
          PROCEDURE :: MCMCParticle_hash_packproposal_
          GENERIC   :: packproposal =>                 &
          &            MCMCParticle_hash_packproposal, &
          &            MCMCParticle_hash_packproposal_

          !!! This function returns the elementIndex element in the hash table
          !!! If the elementKey is present it gets the coordinates from Hash key
          PROCEDURE :: MCMCParticle_hash_elementAt_2D3D
          PROCEDURE :: MCMCParticle_hash_elementAt_2D
          PROCEDURE :: MCMCParticle_hash_elementAt_3D
          GENERIC   :: elementAt => &
          &            MCMCParticle_hash_elementAt_2D3D, &
          &            MCMCParticle_hash_elementAt_2D,   &
          &            MCMCParticle_hash_elementAt_3D

          PROCEDURE :: MCMCParticle_hash_contains
          PROCEDURE :: MCMCParticle_hash_contains_2d
          PROCEDURE :: MCMCParticle_hash_contains_3d
          GENERIC   :: containselement =>             &
          &            MCMCParticle_hash_contains,    &
          &            MCMCParticle_hash_contains_2d, &
          &            MCMCParticle_hash_contains_3d
        END TYPE ppm_rc_MCMCParticlehtable

        TYPE ppm_rc_HashIndextable
          !---------------------------------------------------------------------
          !  Declaration of arrays
          !---------------------------------------------------------------------
          INTEGER(ppm_kind_int64), DIMENSION(:), ALLOCATABLE :: keys
          !!! Array for keeping hash table keys.
          !---------------------------------------------------------------------
          !  Declaration of variables
          !--------------------------------------------------------------------
          INTEGER(ppm_kind_int64)                            :: nrow = 0_ppm_kind_int64
          !!! number of rows in hash table
         CONTAINS
          PROCEDURE :: create  => HashIndex_create_htable
          PROCEDURE :: destroy => HashIndex_destroy_htable
          PROCEDURE :: h_key   => HashIndex_h_key
          PROCEDURE :: HashIndex_hash_insert
          PROCEDURE :: HashIndex_hash_insert_2D3D
          PROCEDURE :: HashIndex_hash_insert_2d
          PROCEDURE :: HashIndex_hash_insert_3d
          GENERIC   :: insert =>                   &
          &            HashIndex_hash_insert,      &
          &            HashIndex_hash_insert_2D3D, &
          &            HashIndex_hash_insert_2d,   &
          &            HashIndex_hash_insert_3d
          PROCEDURE :: HashIndex_hash_search
          PROCEDURE :: HashIndex_hash_search_2d
          PROCEDURE :: HashIndex_hash_search_3d
          GENERIC   :: search =>                 &
          &            HashIndex_hash_search,    &
          &            HashIndex_hash_search_2d, &
          &            HashIndex_hash_search_3d
          PROCEDURE :: HashIndex_hash_remove
          PROCEDURE :: HashIndex_hash_remove_2D3D
          PROCEDURE :: HashIndex_hash_remove_2d
          PROCEDURE :: HashIndex_hash_remove_3d
          GENERIC   :: remove =>                   &
          &            HashIndex_hash_remove,      &
          &            HashIndex_hash_remove_2D3D, &
          &            HashIndex_hash_remove_2d,   &
          &            HashIndex_hash_remove_3d
          PROCEDURE :: grow => HashIndex_grow_htable
          PROCEDURE :: shrink => HashIndex_shrink_htable
          PROCEDURE :: size => HashIndex_hash_size
        END TYPE ppm_rc_HashIndextable

        !!! Private temporary variables
        TYPE(ppm_rc_MCMCParticlehtable), DIMENSION(:), ALLOCATABLE :: Particlehtabletmp

        TYPE ppm_rc_MCMCResultshtable
          !---------------------------------------------------------------------
          !  Declaration of arrays
          !---------------------------------------------------------------------
          INTEGER(ppm_kind_int64), DIMENSION(:), ALLOCATABLE :: keys
          !!! Array for keeping hash table keys.
          TYPE(ppm_rc_list_2),     DIMENSION(:), ALLOCATABLE :: borders_pos
          !!! Array for keeping positions of cells on "borders" array.
          !---------------------------------------------------------------------
          !  Declaration of variables
          !--------------------------------------------------------------------
          INTEGER(ppm_kind_int64)                            :: nrow = 0_ppm_kind_int64
          !!! number of rows in hash table
         CONTAINS
          PROCEDURE :: create  => MCMCResults_create_htable
          PROCEDURE :: destroy => MCMCResults_destroy_htable
          PROCEDURE :: h_key   => MCMCResults_h_key
          PROCEDURE :: MCMCResults_hash_insert
          PROCEDURE :: MCMCResults_hash_insert_2D3D
          PROCEDURE :: MCMCResults_hash_insert_2d
          PROCEDURE :: MCMCResults_hash_insert_3d
          GENERIC   :: insert =>                     &
          &            MCMCResults_hash_insert,      &
          &            MCMCResults_hash_insert_2D3D, &
          &            MCMCResults_hash_insert_2d,   &
          &            MCMCResults_hash_insert_3d
          PROCEDURE :: MCMCResults_hash_search
          PROCEDURE :: MCMCResults_hash_search_2d
          PROCEDURE :: MCMCResults_hash_search_3d
          GENERIC   :: search =>                   &
          &            MCMCResults_hash_search,    &
          &            MCMCResults_hash_search_2d, &
          &            MCMCResults_hash_search_3d
          PROCEDURE :: MCMCResults_hash_remove
          PROCEDURE :: MCMCResults_hash_remove_2D3D
          PROCEDURE :: MCMCResults_hash_remove_2d
          PROCEDURE :: MCMCResults_hash_remove_3d
          GENERIC   :: remove =>                     &
          &            MCMCResults_hash_remove,      &
          &            MCMCResults_hash_remove_2D3D, &
          &            MCMCResults_hash_remove_2d,   &
          &            MCMCResults_hash_remove_3d
          PROCEDURE :: grow => MCMCResults_grow_htable
          PROCEDURE :: size => MCMCResults_hash_size
          !!! check the elements in the hash table and pack the proposal
          !!! in an allocatable array
          !!! This function returns the elementIndex element in the hash table
          !!! If the elementKey is present it gets the coordinates from Hash key
        END TYPE ppm_rc_MCMCResultshtable

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
#ifdef __F2003
        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE HashFunc
          FUNCTION HashFunc(key,seed,tablenrow) BIND(C,NAME='HashFunc')
            IMPORT             :: C_INT64_T
            INTEGER(C_INT64_T) :: key
            INTEGER(C_INT64_T) :: seed
            INTEGER(C_INT64_T) :: tablenrow
            INTEGER(C_INT64_T) :: HashFunc
          END FUNCTION
        END INTERFACE

        INTERFACE HashFunc_XY
          FUNCTION HashFunc_XY(X,Y,seed,tablenrow,XY) BIND(C,NAME='HashFunc_XY')
            IMPORT             :: C_INT32_T,C_INT64_T
            INTEGER(C_INT32_T) :: X
            !!! X coordinate (less than 32768 on one processor)
            INTEGER(C_INT32_T) :: Y
            !!! Y coordinate (less than 32768 on one processor)
            INTEGER(C_INT64_T) :: seed
            INTEGER(C_INT64_T) :: tablenrow
            INTEGER(C_INT64_T) :: XY
            !!! X & Y in one number
            INTEGER(C_INT64_T) :: HashFunc_XY
          END FUNCTION
        END INTERFACE

        INTERFACE HashFunc_XYLabel
          FUNCTION HashFunc_XYLabel(X,Y,Label,seed,tablenrow,XYLabel) BIND(C,NAME='HashFunc_XYLabel')
            IMPORT             :: C_INT32_T,C_INT64_T
            INTEGER(C_INT32_T) :: X
            INTEGER(C_INT32_T) :: Y
            INTEGER(C_INT32_T) :: Label
            INTEGER(C_INT64_T) :: seed
            INTEGER(C_INT64_T) :: tablenrow
            INTEGER(C_INT64_T) :: XYLabel
            INTEGER(C_INT64_T) :: HashFunc_XYLabel
          END FUNCTION
        END INTERFACE

        INTERFACE HashFunc_XYZ
          FUNCTION HashFunc_XYZ(X,Y,Z,seed,tablenrow,XYZ) BIND(C,NAME='HashFunc_XYZ')
            IMPORT             :: C_INT32_T,C_INT64_T
            INTEGER(C_INT32_T) :: X
            !!! X coordinate (less than 32768 on one processor)
            INTEGER(C_INT32_T) :: Y
            !!! Y coordinate (less than 32768 on one processor)
            INTEGER(C_INT32_T) :: Z
            !!! Z coordinate (less than 32768 on one processor)
            INTEGER(C_INT64_T) :: seed
            INTEGER(C_INT64_T) :: tablenrow
            INTEGER(C_INT64_T) :: XYZ
            INTEGER(C_INT64_T) :: HashFunc_XYZ
          END FUNCTION
        END INTERFACE

        INTERFACE HashFunc_XYZLabel
          FUNCTION HashFunc_XYZLabel(X,Y,Z,Label,seed,tablenrow,XYZLabel) BIND(C,NAME='HashFunc_XYZLabel')
            IMPORT             :: C_INT32_T,C_INT64_T
            INTEGER(C_INT32_T) :: X
            INTEGER(C_INT32_T) :: Y
            INTEGER(C_INT32_T) :: Z
            INTEGER(C_INT32_T) :: Label
            INTEGER(C_INT64_T) :: seed
            INTEGER(C_INT64_T) :: tablenrow
            INTEGER(C_INT64_T) :: XYZLabel
            INTEGER(C_INT64_T) :: HashFunc_XYZLabel
          END FUNCTION
        END INTERFACE
#endif

        INTERFACE IndexHashFunctor_label_2d
          MODULE PROCEDURE IndexHashFunctor64_label0_2d
          MODULE PROCEDURE IndexHashFunctor64_label_2d
        END INTERFACE

        INTERFACE IndexHashFunctor_label_3d
          MODULE PROCEDURE IndexHashFunctor64_label0_3d
          MODULE PROCEDURE IndexHashFunctor64_label_3d
        END INTERFACE

        INTERFACE HashIndexFunctor_label_2d
          MODULE PROCEDURE HashIndexFunctor64_label0_2d
          MODULE PROCEDURE HashIndexFunctor64_label_2d
        END INTERFACE

        INTERFACE HashIndexFunctor_label_3d
          MODULE PROCEDURE HashIndexFunctor64_label0_3d
          MODULE PROCEDURE HashIndexFunctor64_label_3d
        END INTERFACE

        INTERFACE IndexHashFunctor64_2d
          MODULE PROCEDURE IndexHashFunctor64_2d
          MODULE PROCEDURE IndexHashFunctor64__2d
        END INTERFACE

        INTERFACE IndexHashFunctor64_3d
          MODULE PROCEDURE IndexHashFunctor64_3d
          MODULE PROCEDURE IndexHashFunctor64__3d
        END INTERFACE
        !----------------------------------------------------------------------
        ! Public
        !----------------------------------------------------------------------
        PUBLIC :: ppm_rc_htable

        PUBLIC :: ppm_rc_MCMCParticlehtable
        PUBLIC :: ppm_rc_HashIndextable
        PUBLIC :: MCMCParticle_realloc
        PUBLIC :: ppm_rc_MCMCResultshtable

      CONTAINS
#define __2D  2
#define __3D  3

#define  DTYPE(a) a/**/_3d
#define __DIME  __3D
#include "./hash/ppm_rc_hashtypeproc.f"
#undef  __DIME
#undef  DTYPE

#define  DTYPE(a) a/**/_2d
#define __DIME  __2D
#include "./hash/ppm_rc_hashtypeproc.f"
#undef  __DIME
#undef  DTYPE

#undef  __2D
#undef  __3D
      END MODULE ppm_rc_module_hash
