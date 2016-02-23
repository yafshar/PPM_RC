      !-------------------------------------------------------------------------
      !  Module       :                 ppm_rc_module_stat
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains subroutines for stat variable
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
