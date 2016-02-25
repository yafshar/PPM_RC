      !-------------------------------------------------------------------------
      !  Module       :                    ppm_rc_module_fire
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
      !  Module       :                    ppm_rc_module_fire
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Spreads fires starting from seed points
      !
      !  Remarks      : The main parallel driving routine is startingfire.
      !                 This routine calls the other subroutines directly or
      !                 indirectly to spread the fire
      !
      !                 Assumptions:
      !                 (1) New label is created based on the global pixel position
      !                     and it would be a unique id.
      !                 (2) The initialization firing assumes that the old
      !                     labels are 1.
      !                 (3) Assumes a mesh ghost layer of 1 grid point
      !                 (4) Need to do a ghost get before calling
      !                     startingfire because it is
      !                     necessary that neighbouring subdomains (along with
      !                     their ghost) form a proper overlapping (through the
      !                     ghost) montage
      !
      !  References   :
      !-------------------------------------------------------------------------
      MODULE ppm_rc_module_fire
        !----------------------------------------------------------------------
        !  Modules
        !----------------------------------------------------------------------
        USE ppm_rc_module_global
        USE ppm_rc_module_linkedlist
        USE ppm_rc_module_stat
        USE ppm_rc_module_topologicalnumber
        USE ppm_rc_module_energy
        USE ppm_rc_module_util
        IMPLICIT NONE

        PRIVATE
        !---------------------------------------------------------------------
        !  Some work memory on the heap
        !---------------------------------------------------------------------
        TYPE(ppm_rc_c_stat), POINTER :: tmp_region_stat => NULL()

        REAL(ppm_kind_double), DIMENSION(:,:), ALLOCATABLE :: tmp_region_stat_compact
        REAL(ppm_kind_double), DIMENSION(:,:), ALLOCATABLE :: tmp_region_stat_aggregate

        INTEGER, DIMENSION(:), ALLOCATABLE :: tmp1_i

        INTEGER, DIMENSION(:), ALLOCATABLE :: seednm
        INTEGER, DIMENSION(:), ALLOCATABLE :: partnm
        INTEGER, DIMENSION(:), ALLOCATABLE :: nlabels
        INTEGER, DIMENSION(:), ALLOCATABLE :: ghost_nlabels

        LOGICAL :: ghostfirenewjob

        !----------------------------------------------------------------------
        ! Public
        !----------------------------------------------------------------------
        PUBLIC :: seednm
        PUBLIC :: partnm
        PUBLIC :: nlabels

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE ppm_rc_floodFill_2d
          MODULE PROCEDURE ppm_rc_floodFillScanline_2d
          MODULE PROCEDURE ppm_rc_floodFillScanline__2d
          MODULE PROCEDURE ppm_rc_floodFillScanlineConditional_2d
        END INTERFACE

        INTERFACE ppm_rc_floodFill_3d
          MODULE PROCEDURE ppm_rc_floodFillScanline_3d
          MODULE PROCEDURE ppm_rc_floodFillScanline__3d
          MODULE PROCEDURE ppm_rc_floodFillScanlineConditional_3d
        END INTERFACE

        !----------------------------------------------------------------------
        ! Public
        !----------------------------------------------------------------------
        PUBLIC :: initforestfire_2d
        PUBLIC :: initforestfire_3d
        PUBLIC :: forestfire_2d
        PUBLIC :: forestfire_3d
        PUBLIC :: ppm_rc_floodFill_2d
        PUBLIC :: ppm_rc_floodFill_3d
        PUBLIC :: ppm_rc_ghost_copy_2d
        PUBLIC :: ppm_rc_ghost_copy_3d
      CONTAINS

#define __2D  2
#define __3D  3

#define  DTYPE(a) a/**/_3d
#define __DIME  __3D
#include "./fire/ppm_rc_forestfire.f"
#include "./fire/ppm_rc_fire.f"
#include "./fire/ppm_rc_ghostfire.f"
#include "./fire/ppm_rc_ghost_copy.f"
#include "./fire/ppm_rc_floodfill.f"
#undef  __DIME
#undef  DTYPE
#define  DTYPE(a) a/**/_2d
#define __DIME  __2D
#include "./fire/ppm_rc_forestfire.f"
#include "./fire/ppm_rc_fire.f"
#include "./fire/ppm_rc_ghostfire.f"
#include "./fire/ppm_rc_ghost_copy.f"
#include "./fire/ppm_rc_floodfill.f"
#undef  __DIME
#undef  DTYPE

#undef  __2D
#undef  __3D

      END MODULE ppm_rc_module_fire
