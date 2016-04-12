      !-------------------------------------------------------------------------
      !  Module       :                    ppm_rc_module_util
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
      !  Module       :                    ppm_rc_module_util
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Module contains utility subroutines.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      MODULE ppm_rc_module_util
        !----------------------------------------------------------------------
        !  Modules
        !----------------------------------------------------------------------
        USE ppm_rc_module_global

#ifdef __F2003
        USE ISO_C_BINDING
#endif
        IMPLICIT NONE

        PRIVATE

        INTEGER(ppm_kind_int64), PARAMETER :: htable_null_li = -HUGE(1_ppm_kind_int64)
        !!! NULL value for hash table
        INTEGER,                 PARAMETER :: seed1 = 738235926
        !!! Hardcoded seed value taken from MurmurHash
        INTEGER,                 PARAMETER :: seed2 = 1243832038
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
          INTEGER                                            :: nrow = 0
          !!! number of rows in hash table
         CONTAINS
          PROCEDURE :: create => create_htable
          PROCEDURE :: destroy => destroy_htable
          PROCEDURE :: h_key1
          PROCEDURE :: h_key2
          GENERIC   :: h_key => h_key1,h_key2
          PROCEDURE :: hash_insert
          PROCEDURE :: hash_insert_
          PROCEDURE :: hash_insert__
          PROCEDURE :: hash_insert_2d
          PROCEDURE :: hash_insert_3d
          GENERIC   :: insert =>       &
          &            hash_insert,    &
          &            hash_insert_,   &
          &            hash_insert__,  &
          &            hash_insert_2d, &
          &            hash_insert_3d
          PROCEDURE :: hash_search
          PROCEDURE :: hash_search_
          PROCEDURE :: hash_search_2d
          PROCEDURE :: hash_search_3d
          GENERIC   :: search =>       &
          &            hash_search,    &
          &            hash_search_,   &
          &            hash_search_2d, &
          &            hash_search_3d
          PROCEDURE :: hash_remove
          PROCEDURE :: hash_remove_
          PROCEDURE :: hash_remove__
          PROCEDURE :: hash_remove_2d
          PROCEDURE :: hash_remove_3d
          GENERIC   :: remove =>       &
          &            hash_remove,    &
          &            hash_remove_,   &
          &            hash_remove__,  &
          &            hash_remove_2d, &
          &            hash_remove_3d
          PROCEDURE :: grow => grow_htable
          PROCEDURE :: shrink => shrink_htable
          PROCEDURE :: size => hash_size
        END TYPE ppm_rc_htable
        !!! This hashtable is for hashing REAL value using the same PPM implemenattion

        TYPE MCMCParticle
          INTEGER  :: candlabel
          REAL(MK) :: proposal
        END TYPE MCMCParticle

        TYPE MCMCHistoryParticle
          INTEGER  :: candlabel
          INTEGER  :: orglabel
          REAL(MK) :: proposal
          LOGICAL  :: wasadded
        END TYPE MCMCHistoryParticle

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
          INTEGER                                            :: nrow = 0
          !!! number of rows in hash table
         CONTAINS
          PROCEDURE :: create  => MCMCParticle_create_htable
          PROCEDURE :: destroy => MCMCParticle_destroy_htable
          PROCEDURE :: MCMCParticle_h_key1
          PROCEDURE :: MCMCParticle_h_key2
          GENERIC   :: h_key  => MCMCParticle_h_key1,MCMCParticle_h_key2
          PROCEDURE :: MCMCParticle_hash_insert
          PROCEDURE :: MCMCParticle_hash_insert_
          PROCEDURE :: MCMCParticle_hash_insert__
          PROCEDURE :: MCMCParticle_hash_insert_2d
          PROCEDURE :: MCMCParticle_hash_insert_3d
          GENERIC   :: insert =>                    &
          &            MCMCParticle_hash_insert,    &
          &            MCMCParticle_hash_insert_,   &
          &            MCMCParticle_hash_insert__,  &
          &            MCMCParticle_hash_insert_2d, &
          &            MCMCParticle_hash_insert_3d
          PROCEDURE :: MCMCParticle_hash_search
          PROCEDURE :: MCMCParticle_hash_search_
          PROCEDURE :: MCMCParticle_hash_search_2d
          PROCEDURE :: MCMCParticle_hash_search_3d
          GENERIC   :: search =>                    &
          &            MCMCParticle_hash_search,    &
          &            MCMCParticle_hash_search_,   &
          &            MCMCParticle_hash_search_2d, &
          &            MCMCParticle_hash_search_3d
          PROCEDURE :: MCMCParticle_hash_remove
          PROCEDURE :: MCMCParticle_hash_remove_
          PROCEDURE :: MCMCParticle_hash_remove__
          PROCEDURE :: MCMCParticle_hash_remove_2d
          PROCEDURE :: MCMCParticle_hash_remove_3d
          GENERIC   :: remove =>                    &
          &            MCMCParticle_hash_remove,    &
          &            MCMCParticle_hash_remove_,   &
          &            MCMCParticle_hash_remove__,  &
          &            MCMCParticle_hash_remove_2d, &
          &            MCMCParticle_hash_remove_3d
          PROCEDURE :: grow => MCMCParticle_grow_htable
          PROCEDURE :: shrink => MCMCParticle_shrink_htable
          PROCEDURE :: size => MCMCParticle_hash_size
          !!! check the elements in the hash table and pack the proposal
          !!! in an allocatable array
          PROCEDURE :: packproposal => MCMCParticle_hash_packproposal
          !!! This function returns the elementIndex element in the hash table
          !!! If the elementKey is present it gets the coordinates from Hash key
          PROCEDURE :: elementAt => MCMCParticle_hash_elementAt
          PROCEDURE :: MCMCParticle_hash_contains
          PROCEDURE :: MCMCParticle_hash_contains_
          PROCEDURE :: MCMCParticle_hash_contains_2d
          PROCEDURE :: MCMCParticle_hash_contains_3d
          GENERIC   :: containselement =>             &
          &            MCMCParticle_hash_contains,    &
          &            MCMCParticle_hash_contains_,   &
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
          INTEGER                                            :: nrow = 0
          !!! number of rows in hash table
         CONTAINS
          PROCEDURE :: create  => HashIndex_create_htable
          PROCEDURE :: destroy => HashIndex_destroy_htable
          PROCEDURE :: HashIndex_h_key1
          PROCEDURE :: HashIndex_h_key2
          GENERIC   :: h_key => HashIndex_h_key1,HashIndex_h_key2
          PROCEDURE :: HashIndex_hash_insert
          PROCEDURE :: HashIndex_hash_insert_
          PROCEDURE :: HashIndex_hash_insert__
          PROCEDURE :: HashIndex_hash_insert_2d
          PROCEDURE :: HashIndex_hash_insert_3d
          GENERIC   :: insert =>                 &
          &            HashIndex_hash_insert,    &
          &            HashIndex_hash_insert_,   &
          &            HashIndex_hash_insert__,  &
          &            HashIndex_hash_insert_2d, &
          &            HashIndex_hash_insert_3d
          PROCEDURE :: HashIndex_hash_search
          PROCEDURE :: HashIndex_hash_search_
          PROCEDURE :: HashIndex_hash_search_2d
          PROCEDURE :: HashIndex_hash_search_3d
          GENERIC   :: search =>                 &
          &            HashIndex_hash_search,    &
          &            HashIndex_hash_search_,   &
          &            HashIndex_hash_search_2d, &
          &            HashIndex_hash_search_3d
          PROCEDURE :: HashIndex_hash_remove
          PROCEDURE :: HashIndex_hash_remove_
          PROCEDURE :: HashIndex_hash_remove__
          PROCEDURE :: HashIndex_hash_remove_2d
          PROCEDURE :: HashIndex_hash_remove_3d
          GENERIC   :: remove =>                 &
          &            HashIndex_hash_remove,    &
          &            HashIndex_hash_remove_,   &
          &            HashIndex_hash_remove__,  &
          &            HashIndex_hash_remove_2d, &
          &            HashIndex_hash_remove_3d
          PROCEDURE :: grow => HashIndex_grow_htable
          PROCEDURE :: shrink => HashIndex_shrink_htable
          PROCEDURE :: size => HashIndex_hash_size
        END TYPE ppm_rc_HashIndextable

        !!! Private temporary variables
        TYPE(ppm_rc_MCMCParticlehtable), DIMENSION(:), ALLOCATABLE :: Particlehtabletmp

#ifdef __Linux
        INTEGER :: valueRSS0
#endif

#ifdef __F2003
        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE HashKey
          FUNCTION HashKey(key,seed,tablenrow) BIND(C,NAME='HashKey')
            IMPORT             :: C_INT32_T,C_INT64_T
            INTEGER(C_INT64_T) :: key
            INTEGER(C_INT32_T) :: seed
            INTEGER(C_INT32_T) :: tablenrow
            INTEGER(C_INT32_T) :: HashKey
          END FUNCTION
        END INTERFACE

        INTERFACE HashKey_XY
          FUNCTION HashKey_XY(X,Y,seed,tablenrow,XY) BIND(C,NAME='HashKey_XY')
            IMPORT             :: C_INT32_T,C_INT64_T
            INTEGER(C_INT32_T) :: X
            !!! X coordinate (less than 32768 on one processor)
            INTEGER(C_INT32_T) :: Y
            !!! Y coordinate (less than 32768 on one processor)
            INTEGER(C_INT32_T) :: seed
            INTEGER(C_INT32_T) :: tablenrow
            INTEGER(C_INT64_T) :: XY
            !!! X & Y in one number
            INTEGER(C_INT32_T) :: HashKey_XY
          END FUNCTION
        END INTERFACE

        INTERFACE HashKey_XYLabel
          FUNCTION HashKey_XYLabel(X,Y,Label,seed,tablenrow,XYLabel) BIND(C,NAME='HashKey_XYLabel')
            IMPORT             :: C_INT32_T,C_INT64_T
            INTEGER(C_INT32_T) :: X
            INTEGER(C_INT32_T) :: Y
            INTEGER(C_INT32_T) :: Label
            INTEGER(C_INT32_T) :: seed
            INTEGER(C_INT32_T) :: tablenrow
            INTEGER(C_INT64_T) :: XYLabel
            INTEGER(C_INT32_T) :: HashKey_XYLabel
          END FUNCTION
        END INTERFACE

        INTERFACE HashKey_XYZ
          FUNCTION HashKey_XYZ(X,Y,Z,seed,tablenrow,XYZ) BIND(C,NAME='HashKey_XYZ')
            IMPORT             :: C_INT32_T,C_INT64_T
            INTEGER(C_INT32_T) :: X
            !!! X coordinate (less than 32768 on one processor)
            INTEGER(C_INT32_T) :: Y
            !!! Y coordinate (less than 32768 on one processor)
            INTEGER(C_INT32_T) :: Z
            !!! Z coordinate (less than 32768 on one processor)
            INTEGER(C_INT32_T) :: seed
            INTEGER(C_INT32_T) :: tablenrow
            INTEGER(C_INT64_T) :: XYZ
            INTEGER(C_INT32_T) :: HashKey_XYZ
          END FUNCTION
        END INTERFACE

        INTERFACE HashKey_XYZLabel
          FUNCTION HashKey_XYZLabel(X,Y,Z,Label,seed,tablenrow,XYZLabel) BIND(C,NAME='HashKey_XYZLabel')
            IMPORT             :: C_INT32_T,C_INT64_T
            INTEGER(C_INT32_T) :: X
            INTEGER(C_INT32_T) :: Y
            INTEGER(C_INT32_T) :: Z
            INTEGER(C_INT32_T) :: Label
            INTEGER(C_INT32_T) :: seed
            INTEGER(C_INT32_T) :: tablenrow
            INTEGER(C_INT64_T) :: XYZLabel
            INTEGER(C_INT32_T) :: HashKey_XYZLabel
          END FUNCTION
        END INTERFACE
#endif
        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE ppm_rc_normalize_2d
          MODULE PROCEDURE ppm_rc_normalize1_2d
          MODULE PROCEDURE ppm_rc_normalize2_2d
        END INTERFACE

        INTERFACE ppm_rc_normalize_3d
          MODULE PROCEDURE ppm_rc_normalize1_3d
          MODULE PROCEDURE ppm_rc_normalize2_3d
        END INTERFACE

        INTERFACE ppm_rc_CopyImageAndNormalize_2d
          MODULE PROCEDURE ppm_rc_CopyImageAndNormalize_2d
        END INTERFACE

        INTERFACE ppm_rc_CopyImageAndNormalize_3d
          MODULE PROCEDURE ppm_rc_CopyImageAndNormalize_3d
        END INTERFACE

        INTERFACE id_ltg_2d
          MODULE PROCEDURE ppm_rc_id_LocalToGlobal_i_2d
          MODULE PROCEDURE ppm_rc_id_LocalToGlobal_li_2d
        END INTERFACE

        INTERFACE id_ltg_3d
          MODULE PROCEDURE ppm_rc_id_LocalToGlobal_i_3d
          MODULE PROCEDURE ppm_rc_id_LocalToGlobal_li_3d
        END INTERFACE

        INTERFACE id_gtl_2d
          MODULE PROCEDURE ppm_rc_id_GlobalToLocal_i_2d
          MODULE PROCEDURE ppm_rc_id_GlobalToLocal_li_2d
          MODULE PROCEDURE ppm_rc_id_GlobalToLocal__i_2d
          MODULE PROCEDURE ppm_rc_id_GlobalToLocal__li_2d
        END INTERFACE

        INTERFACE id_gtl_3d
          MODULE PROCEDURE ppm_rc_id_GlobalToLocal_i_3d
          MODULE PROCEDURE ppm_rc_id_GlobalToLocal_li_3d
          MODULE PROCEDURE ppm_rc_id_GlobalToLocal__i_3d
          MODULE PROCEDURE ppm_rc_id_GlobalToLocal__li_3d
        END INTERFACE

        INTERFACE ppm_rc_unique
          MODULE PROCEDURE ppm_rc_uniquers
          MODULE PROCEDURE ppm_rc_uniquer
          MODULE PROCEDURE ppm_rc_uniquei
          MODULE PROCEDURE ppm_rc_uniqueli
        END INTERFACE

        INTERFACE ppm_rc_label_exist
          MODULE PROCEDURE ppm_rc_label_exist_asc
          MODULE PROCEDURE ppm_rc_label_exist_asc2
          MODULE PROCEDURE ppm_rc_label_exist_asc3
          MODULE PROCEDURE ppm_rc_label_exist_asc4
          MODULE PROCEDURE ppm_rc_label_exist_dsc
          MODULE PROCEDURE ppm_rc_label_exist_dsc2
          MODULE PROCEDURE ppm_rc_label_exist_dsc3
          MODULE PROCEDURE ppm_rc_label_exist_dsc4
        END INTERFACE

        INTERFACE IndexHashFunctor_2d
          MODULE PROCEDURE IndexHashFunctor32_2d
          MODULE PROCEDURE IndexHashFunctor32__2d
        END INTERFACE
        INTERFACE IndexHashFunctor_3d
          MODULE PROCEDURE IndexHashFunctor32_3d
          MODULE PROCEDURE IndexHashFunctor32__3d
        END INTERFACE
        INTERFACE IndexHashFunctor64_2d
          MODULE PROCEDURE IndexHashFunctor64_2d
          MODULE PROCEDURE IndexHashFunctor64__2d
        END INTERFACE
        INTERFACE IndexHashFunctor64_3d
          MODULE PROCEDURE IndexHashFunctor64_3d
          MODULE PROCEDURE IndexHashFunctor64__3d
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

        INTERFACE ppm_rc_label_index
          MODULE PROCEDURE ppm_rc_label_index
        END INTERFACE

        INTERFACE ppm_rc_shuffle
          MODULE PROCEDURE ppm_rc_FisherYatesShuffle
          MODULE PROCEDURE ppm_rc_FisherYatesShuffle_
        END INTERFACE

        !----------------------------------------------------------------------
        ! Public
        !----------------------------------------------------------------------
        PUBLIC :: ppm_rc_htable

        PUBLIC :: MCMCParticle
        PUBLIC :: ppm_rc_MCMCParticlehtable
        PUBLIC :: MCMCHistoryParticle
        PUBLIC :: ppm_rc_HashIndextable
        PUBLIC :: MCMCParticle_realloc

        PUBLIC :: ppm_rc_normalize_2d
        PUBLIC :: ppm_rc_normalize_3d
        PUBLIC :: ppm_rc_CopyImageAndNormalize_2d
        PUBLIC :: ppm_rc_CopyImageAndNormalize_3d
        PUBLIC :: id_ltg_2d,id_ltg_3d
        PUBLIC :: id_gtl_2d,id_gtl_3d
        PUBLIC :: IndexHashFunctor_2d
        PUBLIC :: IndexHashFunctor_3d
        PUBLIC :: IndexHashFunctor64_2d
        PUBLIC :: IndexHashFunctor64_3d
        PUBLIC :: HashIndexFunctor_2d
        PUBLIC :: HashIndexFunctor_3d
        PUBLIC :: ppm_rc_unique
        PUBLIC :: ppm_rc_uppercase
        PUBLIC :: ppm_rc_compute_num_common
#ifdef __Linux
        PUBLIC :: ppm_rc_mem_usage
        PUBLIC :: valueRSS0
#endif
        PUBLIC :: ppm_rc_label_exist
        PUBLIC :: ppm_rc_label_index

        PUBLIC :: ppm_rc_shuffle

      CONTAINS
#define __2D  2
#define __3D  3

#define  DTYPE(a) a/**/_3d
#define __DIME  __3D
#include "./util/ppm_rc_CopyImageAndNormalize.f"
#include "./util/ppm_rc_normalize.f"
#include "./util/ppm_rc_id.f"
#undef  __DIME
#undef  DTYPE

#define  DTYPE(a) a/**/_2d
#define __DIME  __2D
#include "./util/ppm_rc_CopyImageAndNormalize.f"
#include "./util/ppm_rc_normalize.f"
#include "./util/ppm_rc_id.f"
#undef  __DIME
#undef  DTYPE

#undef  __2D
#undef  __3D

#include "./util/ppm_rc_uppercase.f"
#include "./util/ppm_rc_num_common.f"
#include "./util/ppm_rc_unique.f"
#ifdef __Linux
#include "./util/ppm_rc_mem_usage.f"
#endif

#include "./util/ppm_rc_inl_hash.f"
#include "./util/ppm_rc_label_exist.f"
#include "./util/ppm_rc_label_index.f"
#include "./util/ppm_rc_FisherYatesShuffle.f"

      END MODULE ppm_rc_module_util
