      !-------------------------------------------------------------------------
      !  Module       :                    ppm_rc_module_finalize
      !-------------------------------------------------------------------------
      !
      !  Purpose      :
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  MOSAIC Group
      !  Max Planck Institute of Molecular Cell Biology and Genetics
      !  Pfotenhauerstr. 108, 01307 Dresden, Germany
      !
      !  Author           - y.afshar           June   2014
      !-------------------------------------------------------------------------
      MODULE ppm_rc_module_finalize
        !----------------------------------------------------------------------
        !  Modules
        !----------------------------------------------------------------------
        USE ppm_rc_module_global
        IMPLICIT NONE

        PRIVATE
        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        !----------------------------------------------------------------------
        ! Public
        !----------------------------------------------------------------------
        PUBLIC :: ppm_rc_finalize

      CONTAINS

#include "./fin/ppm_rc_finalize.f"

      END MODULE ppm_rc_module_finalize

