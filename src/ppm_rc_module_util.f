      !-------------------------------------------------------------------------
      !  Module   :                  ppm_rc_module_util
      !-------------------------------------------------------------------------
      !  MOSAIC Group
      !  Max Planck Institute of Molecular Cell Biology and Genetics
      !  Pfotenhauerstr. 108, 01307 Dresden, Germany
      !
      !  Author           - y.afshar           June   2014
      !-------------------------------------------------------------------------
      MODULE ppm_rc_module_util
        !----------------------------------------------------------------------
        !  Modules
        !----------------------------------------------------------------------
        USE ppm_rc_module_global
        IMPLICIT NONE

        PRIVATE

        INTEGER(ppm_kind_int64), PARAMETER :: htable_null_li = -1_ppm_kind_int64
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
          INTEGER                                            :: nrow = 0
          !!! number of rows in hash table
         CONTAINS
          PROCEDURE :: create => create_htable
          PROCEDURE :: destroy => destroy_htable
          PROCEDURE :: h_func
          PROCEDURE :: h_key
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
          PROCEDURE :: h_func => MCMCParticle_h_func
          PROCEDURE :: h_key  => MCMCParticle_h_key
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
          PROCEDURE :: size => MCMCParticle_hash_size
          !!! check the elements in the hash table and pack the proposal in an allocatable array
          PROCEDURE :: packproposal => MCMCParticle_hash_packproposal
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

        TYPE ppm_rc_MCMCHistoryParticlehtable
          !---------------------------------------------------------------------
          !  Declaration of arrays
          !---------------------------------------------------------------------
          INTEGER(ppm_kind_int64),   DIMENSION(:), ALLOCATABLE :: keys
          !!! Array for keeping hash table keys.
          TYPE(MCMCHistoryParticle), DIMENSION(:), ALLOCATABLE :: borders_pos
          !!! Array for keeping positions of cells on "borders" array.
          !---------------------------------------------------------------------
          !  Declaration of variables
          !--------------------------------------------------------------------
          INTEGER                                              :: nrow = 0
          !!! number of rows in hash table
         CONTAINS
          PROCEDURE :: create  => MCMCHistoryParticle_create_htable
          PROCEDURE :: destroy => MCMCHistoryParticle_destroy_htable
          PROCEDURE :: h_func => MCMCHistoryParticle_h_func
          PROCEDURE :: h_key  => MCMCHistoryParticle_h_key
          PROCEDURE :: MCMCHistoryParticle_hash_insert
          PROCEDURE :: MCMCHistoryParticle_hash_insert_
          PROCEDURE :: MCMCHistoryParticle_hash_insert__
          PROCEDURE :: MCMCHistoryParticle_hash_insert_2d
          PROCEDURE :: MCMCHistoryParticle_hash_insert_3d
          GENERIC   :: insert =>                           &
          &            MCMCHistoryParticle_hash_insert,    &
          &            MCMCHistoryParticle_hash_insert_,   &
          &            MCMCHistoryParticle_hash_insert__,  &
          &            MCMCHistoryParticle_hash_insert_2d, &
          &            MCMCHistoryParticle_hash_insert_3d
          PROCEDURE :: MCMCHistoryParticle_hash_search
          PROCEDURE :: MCMCHistoryParticle_hash_search_
          PROCEDURE :: MCMCHistoryParticle_hash_search_2d
          PROCEDURE :: MCMCHistoryParticle_hash_search_3d
          GENERIC   :: search =>                           &
          &            MCMCHistoryParticle_hash_search,    &
          &            MCMCHistoryParticle_hash_search_,   &
          &            MCMCHistoryParticle_hash_search_2d, &
          &            MCMCHistoryParticle_hash_search_3d
          PROCEDURE :: MCMCHistoryParticle_hash_remove
          PROCEDURE :: MCMCHistoryParticle_hash_remove_
          PROCEDURE :: MCMCHistoryParticle_hash_remove__
          PROCEDURE :: MCMCHistoryParticle_hash_remove_2d
          PROCEDURE :: MCMCHistoryParticle_hash_remove_3d
          GENERIC   :: remove =>                           &
          &            MCMCHistoryParticle_hash_remove,    &
          &            MCMCHistoryParticle_hash_remove_,   &
          &            MCMCHistoryParticle_hash_remove__,  &
          &            MCMCHistoryParticle_hash_remove_2d, &
          &            MCMCHistoryParticle_hash_remove_3d
          PROCEDURE :: grow => MCMCHistoryParticle_grow_htable
          PROCEDURE :: size => MCMCHistoryParticle_hash_size
        END TYPE ppm_rc_MCMCHistoryParticlehtable

        TYPE ppm_rc_HashIndex
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
          PROCEDURE :: h_func => HashIndex_h_func
          PROCEDURE :: h_key  => HashIndex_h_key
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
          PROCEDURE :: size => HashIndex_hash_size
        END TYPE ppm_rc_HashIndex

        !!! Private temporary variables
        TYPE(ppm_rc_MCMCParticlehtable),        DIMENSION(:), ALLOCATABLE :: Particlehtabletmp
        TYPE(ppm_rc_MCMCHistoryParticlehtable), DIMENSION(:), ALLOCATABLE :: HistoryParticlehtabletmp

#ifdef __Linux
        INTEGER :: valueRSS0
#endif

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE normalize_2d
          MODULE PROCEDURE ppm_rc_normalize1_2d
          MODULE PROCEDURE ppm_rc_normalize2_2d
        END INTERFACE

        INTERFACE normalize_3d
          MODULE PROCEDURE ppm_rc_normalize1_3d
          MODULE PROCEDURE ppm_rc_normalize2_3d
        END INTERFACE

        INTERFACE CopyImageAndNormalize_2d
          MODULE PROCEDURE ppm_rc_CopyImageAndNormalize_2d
        END INTERFACE

        INTERFACE CopyImageAndNormalize_3d
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

        INTERFACE unique
          MODULE PROCEDURE uniquers
          MODULE PROCEDURE uniquer
          MODULE PROCEDURE uniquei
          MODULE PROCEDURE uniqueli
        END INTERFACE

        INTERFACE label_exist
          MODULE PROCEDURE ppm_rc_label_exist_asc
          MODULE PROCEDURE ppm_rc_label_exist_asc2
          MODULE PROCEDURE ppm_rc_label_exist_dsc
          MODULE PROCEDURE ppm_rc_label_exist_dsc2
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

        !----------------------------------------------------------------------
        ! Public
        !----------------------------------------------------------------------
        PUBLIC :: ppm_rc_htable

        PUBLIC :: MCMCParticle
        PUBLIC :: ppm_rc_MCMCParticlehtable

        PUBLIC :: MCMCHistoryParticle
        PUBLIC :: ppm_rc_MCMCHistoryParticlehtable

        PUBLIC :: ppm_rc_HashIndex

        PUBLIC :: MCMCParticle_realloc
        PUBLIC :: MCMCHistoryParticle_realloc

        PUBLIC :: normalize_2d
        PUBLIC :: normalize_3d
        PUBLIC :: CopyImageAndNormalize_2d
        PUBLIC :: CopyImageAndNormalize_3d
        PUBLIC :: id_ltg_2d,id_ltg_3d
        PUBLIC :: id_gtl_2d,id_gtl_3d
        PUBLIC :: IndexHashFunctor_2d
        PUBLIC :: IndexHashFunctor_3d
        PUBLIC :: IndexHashFunctor64_2d
        PUBLIC :: IndexHashFunctor64_3d
        PUBLIC :: HashIndexFunctor_2d
        PUBLIC :: HashIndexFunctor_3d
        PUBLIC :: unique
        PUBLIC :: ppm_rc_uppercase
        PUBLIC :: compute_num_common
#ifdef __Linux
        PUBLIC :: ppm_rc_mem_usage
        PUBLIC :: valueRSS0
#endif
        PUBLIC :: label_exist

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

      END MODULE ppm_rc_module_util
