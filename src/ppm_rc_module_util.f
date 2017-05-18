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
        IMPLICIT NONE

        PRIVATE

#ifdef __Linux
        INTEGER :: valueRSS0
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
        PUBLIC :: ppm_rc_normalize_2d
        PUBLIC :: ppm_rc_normalize_3d

        PUBLIC :: ppm_rc_CopyImageAndNormalize_2d
        PUBLIC :: ppm_rc_CopyImageAndNormalize_3d

        PUBLIC :: id_ltg_2d,id_ltg_3d
        PUBLIC :: id_gtl_2d,id_gtl_3d

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

#include "./util/ppm_rc_label_exist.f"
#include "./util/ppm_rc_label_index.f"
#include "./util/ppm_rc_FisherYatesShuffle.f"
#include "./util/ppm_rc_color_vertex.f"

      END MODULE ppm_rc_module_util
