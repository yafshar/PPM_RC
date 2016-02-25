      !-------------------------------------------------------------------------
      !  Module       :                    ppm_rc_module_stat
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
      !  Module       :                    ppm_rc_module_stat
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains subroutines for stat variable
      !
      !  Remarks      :
      !
      !  References   :
      !-------------------------------------------------------------------------
      MODULE ppm_rc_module_stat
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

        USE ppm_rc_module_global
        IMPLICIT NONE

        PRIVATE
        !----------------------------------------------------------------------
        !
        !----------------------------------------------------------------------
#include "./stat/ppm_rc_statypedef.f"

        !----------------------------------------------------------------------
        ! Public
        !----------------------------------------------------------------------
        PUBLIC :: ppm_rc_stat
        PUBLIC :: ppm_rc_c_stat

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------

      CONTAINS

#include "./stat/ppm_rc_statypeproc.f"

      END MODULE ppm_rc_module_stat
