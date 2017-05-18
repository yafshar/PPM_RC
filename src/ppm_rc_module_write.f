      !-------------------------------------------------------------------------
      !  Module       :                    ppm_rc_module_write
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
      !  Module       :                    ppm_rc_module_write
      !-------------------------------------------------------------------------
      !
      !  Purpose      :  Module contains write subroutines.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      MODULE ppm_rc_module_write
        !----------------------------------------------------------------------
        !  Modules
        !----------------------------------------------------------------------
        USE ppm_rc_module_global
        IMPLICIT NONE

        PRIVATE
        !---------------------------------------------------------------------
        !  Some work memory on the heap
        !---------------------------------------------------------------------
        INTEGER, DIMENSION(:), ALLOCATABLE :: nlabels
        INTEGER, DIMENSION(:), ALLOCATABLE :: bufi

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        !!! This routine will only write the labels of each regoin
        !!! as whight points in a black background

        !----------------------------------------------------------------------
        ! Public
        !----------------------------------------------------------------------
        PUBLIC :: ppm_rc_write_image_2d
        PUBLIC :: ppm_rc_write_image_3d
        PUBLIC :: ppm_rc_write_image_label_2d
        PUBLIC :: ppm_rc_write_image_label_3d
        PUBLIC :: ppm_rc_read_write_image_2d
        PUBLIC :: ppm_rc_read_write_image_3d
        PUBLIC :: ppm_rc_read_write_unique_label_2d
        PUBLIC :: ppm_rc_read_write_unique_label_3d

      CONTAINS

#define __2D  2
#define __3D  3

#define  DTYPE(a) a/**/_3d
#define __DIME  __3D
#include "./io/ppm_rc_write_image.f"
#include "./io/ppm_rc_read_write_image.f"
#include "./io/ppm_rc_read_write_unique_label.f"
#undef  __DIME
#undef  DTYPE

#define  DTYPE(a) a/**/_2d
#define __DIME  __2D
#include "./io/ppm_rc_write_image.f"
#include "./io/ppm_rc_read_write_image.f"
#include "./io/ppm_rc_read_write_unique_label.f"
#undef  __DIME
#undef  DTYPE
#undef  __2D
#undef  __3D

      END MODULE ppm_rc_module_write