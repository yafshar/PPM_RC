      !-------------------------------------------------------------------------
      !  Module       :                    ppm_rc_module_filter
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
      !  Module       :                    ppm_rc_module_filter
      !-------------------------------------------------------------------------
      !
      !  Purpose      :  Module contains (convolution, filters) subroutines.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      MODULE ppm_rc_module_filter
        !----------------------------------------------------------------------
        !  Modules
        !----------------------------------------------------------------------
        USE ppm_rc_module_global
        IMPLICIT NONE

        PRIVATE
        !----------------------------------------------------------------------
        !  Some work memory on the heap
        !----------------------------------------------------------------------
        CLASS(ppm_t_field_), POINTER :: fld1
        CLASS(ppm_t_field_), POINTER :: fld2

        REAL(MK), DIMENSION(:),         ALLOCATABLE :: tmp1_r
        REAL(MK), DIMENSION(:,:),       ALLOCATABLE :: tmp2_r
        REAL(MK), DIMENSION(:,:,:),     ALLOCATABLE :: tmp3_r
        REAL(MK), DIMENSION(:,:,:,:),   ALLOCATABLE :: tmp4_r

        REAL(ppm_kind_double), DIMENSION(:), ALLOCATABLE :: tmp1_rd

        INTEGER, DIMENSION(:),         ALLOCATABLE :: tmp1_i
        INTEGER, DIMENSION(:,:),       ALLOCATABLE :: tmp2_i
        INTEGER, DIMENSION(:,:,:),     ALLOCATABLE :: tmp3_i
        INTEGER, DIMENSION(:,:,:,:),   ALLOCATABLE :: tmp4_i

        INTEGER,                       DIMENSION(1) :: ldc

#include "./filter/ppm_rc_filtertypedef.f"

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE ppm_rc_convolve
          MODULE PROCEDURE ppm_rc_convolve_i_2d
          MODULE PROCEDURE ppm_rc_convolve_i_3d
          MODULE PROCEDURE ppm_rc_convolve_r_2d
          MODULE PROCEDURE ppm_rc_convolve_r_3d
        END INTERFACE

        INTERFACE ppm_rc_GaussianCurvature_2d
          MODULE PROCEDURE ppm_rc_gc_2d
        END INTERFACE
        INTERFACE ppm_rc_GaussianCurvature_3d
          MODULE PROCEDURE ppm_rc_gc_3d
        END INTERFACE

        INTERFACE ppm_rc_median_2d
          MODULE PROCEDURE ppm_rc_median_2d
        END INTERFACE
        INTERFACE ppm_rc_median_3d
          MODULE PROCEDURE ppm_rc_median_3d
        END INTERFACE

        INTERFACE ppm_rc_GaussianImageFilter_2d
          MODULE PROCEDURE ppm_rc_GaussianImageFilter_2d
        END INTERFACE
        INTERFACE ppm_rc_GaussianImageFilter_3d
          MODULE PROCEDURE ppm_rc_GaussianImageFilter_3d
        END INTERFACE

        INTERFACE ppm_rc_SobelImageFilter_2d
          MODULE PROCEDURE ppm_rc_SobelImageFilter_2d
        END INTERFACE
        INTERFACE ppm_rc_SobelImageFilter_3d
          MODULE PROCEDURE ppm_rc_SobelImageFilter_3d
        END INTERFACE

        INTERFACE ppm_rc_EdgeDetection_2d
          MODULE PROCEDURE ppm_rc_EdgeDetection_2d
        END INTERFACE
        INTERFACE ppm_rc_EdgeDetection_3d
          MODULE PROCEDURE ppm_rc_EdgeDetection_3d
        END INTERFACE
        !----------------------------------------------------------------------
        ! Public
        !----------------------------------------------------------------------
        !!! Gaussian curvature
        PUBLIC :: ppm_rc_GaussianCurvature_2d
        PUBLIC :: ppm_rc_GaussianCurvature_3d
        !!! Convolution or correlation
        PUBLIC :: ppm_rc_convolve
        !!! Median
        PUBLIC :: ppm_rc_median_2d
        PUBLIC :: ppm_rc_median_3d
        !!!
        PUBLIC :: ThresholdImageFilter
        PUBLIC :: OtsuThresholdImageFilter

        PUBLIC :: ppm_rc_GaussianImageFilter_2d
        PUBLIC :: ppm_rc_GaussianImageFilter_3d

        PUBLIC :: ppm_rc_SobelImageFilter_2d
        PUBLIC :: ppm_rc_SobelImageFilter_3d

        PUBLIC :: ppm_rc_EdgeDetection_2d
        PUBLIC :: ppm_rc_EdgeDetection_3d

      CONTAINS

#define __2D  2
#define __3D  3

#define  DTYPE(a) a/**/_3d
#define __DIME  __3D
#include "./filter/ppm_rc_filtertypeproc.f"

#include "./filter/ppm_rc_gc.f"
#include "./filter/ppm_rc_median.f"
#define __TYPE INTEGER
#define  CTYPE(a) a/**/_i
#define __ZERO 0
#include "./filter/ppm_rc_convolve.f"
#undef __ZERO
#undef __TYPE
#undef CTYPE
#define __TYPE REAL(MK)
#define  CTYPE(a) a/**/_r
#define __ZERO zero
#include "./filter/ppm_rc_convolve.f"
#undef __ZERO
#undef __TYPE
#undef CTYPE
#define __ZERO zero
#include "./filter/ppm_rc_GaussianImageFilter.f"
#include "./filter/ppm_rc_SobelImageFilter.f"
#include "./filter/ppm_rc_EdgeDetection.f"
#undef __ZERO
#undef  __DIME
#undef  DTYPE


#define  DTYPE(a) a/**/_2d
#define __DIME  __2D
#include "./filter/ppm_rc_filtertypeproc.f"

#include "./filter/ppm_rc_gc.f"
#include "./filter/ppm_rc_median.f"
#define __TYPE INTEGER
#define __ZERO 0
#define  CTYPE(a) a/**/_i
#include "./filter/ppm_rc_convolve.f"
#undef __ZERO
#undef __TYPE
#undef CTYPE
#define __TYPE REAL(MK)
#define  CTYPE(a) a/**/_r
#define __ZERO zero
#include "./filter/ppm_rc_convolve.f"
#undef __ZERO
#undef __TYPE
#undef CTYPE
#define __ZERO zero
#include "./filter/ppm_rc_GaussianImageFilter.f"
#include "./filter/ppm_rc_SobelImageFilter.f"
#include "./filter/ppm_rc_EdgeDetection.f"
#undef __ZERO
#undef  __DIME
#undef  DTYPE

#undef  __2D
#undef  __3D

      END MODULE ppm_rc_module_filter