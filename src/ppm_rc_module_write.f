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
      !  MOSAIC Group
      !  Max Planck Institute of Molecular Cell Biology and Genetics
      !  Pfotenhauerstr. 108, 01307 Dresden, Germany
      !
      !  Author           - y.afshar           June   2014
      !-------------------------------------------------------------------------

      MODULE ppm_rc_module_write
        !----------------------------------------------------------------------
        !  Modules
        !----------------------------------------------------------------------
        USE ppm_rc_module_global
        IMPLICIT NONE

        PRIVATE
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

      CONTAINS

#define __2D  2
#define __3D  3

#define  DTYPE(a) a/**/_3d
#define __DIME  __3D
#include "./io/ppm_rc_write_image.f"
#include "./io/ppm_rc_read_write_image.f"
#undef  __DIME
#undef  DTYPE

#define  DTYPE(a) a/**/_2d
#define __DIME  __2D
#include "./io/ppm_rc_write_image.f"
#include "./io/ppm_rc_read_write_image.f"
#undef  __DIME
#undef  DTYPE
#undef  __2D
#undef  __3D

      END MODULE ppm_rc_module_write