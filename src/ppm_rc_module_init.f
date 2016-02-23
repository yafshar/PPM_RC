      !-------------------------------------------------------------------------
      !  Module   :                  ppm_rc_module_init
      !-------------------------------------------------------------------------
      !  MOSAIC Group
      !  Max Planck Institute of Molecular Cell Biology and Genetics
      !  Pfotenhauerstr. 108, 01307 Dresden, Germany
      !
      !  Author           - y.afshar           June   2013
      !-------------------------------------------------------------------------
      MODULE ppm_rc_module_init
        !----------------------------------------------------------------------
        !  Modules
        !----------------------------------------------------------------------
        USE ppm_rc_module_global
        USE ppm_rc_module_linkedlist, ONLY : ppm_rc_link,ppm_rc_list,  &
        &   ppm_rc_seeds,InnerContourContainer,Candidates,m_Seeds, &
        &   CompetingRegions
        USE ppm_rc_module_energy
        IMPLICIT NONE

        PRIVATE
        !----------------------------------------------------------------------
        !  Some work memory on the heap
        !----------------------------------------------------------------------
#include "./init/ppm_rc_initypedef.f"

        TYPE(RCInitClass), POINTER :: e_Init => NULL()
        !!! Initialization element

        INTEGER :: vInitKind
        !!! InitializationKindType

        !----------------------------------------------------------------------
        ! Public
        !----------------------------------------------------------------------
        PUBLIC :: vInitKind
        PUBLIC :: e_Init

        !----------------------------------------------------------------------
        ! Public
        !----------------------------------------------------------------------
        PUBLIC :: ppm_rc_defaults
        PUBLIC :: ppm_rc_check_ctrl
        PUBLIC :: ppm_rc_init_arg
        PUBLIC :: ppm_rc_dump_parameters
        PUBLIC :: ppm_rc_init_seg_2d
        PUBLIC :: ppm_rc_init_seg_3d
        PUBLIC :: ppm_rc_init_mcmc_2d
        PUBLIC :: ppm_rc_init_mcmc_3d

      CONTAINS

#include "./init/ppm_rc_initypeproc.f"

#define __2D  2
#define __3D  3
#define  DTYPE(a) a/**/_3d
#define __DIME  __3D
#include "./init/ppm_rc_getoutput.f"
#undef  __DIME
#undef  DTYPE
#define  DTYPE(a) a/**/_2d
#define __DIME  __2D
#include "./init/ppm_rc_getoutput.f"
#undef  __DIME
#undef  DTYPE
#undef  __2D
#undef  __3D

#include "./init/ppm_rc_defaults.f"
#include "./init/ppm_rc_check_ctrl.f"
#include "./init/ppm_rc_init_arg.f"
#include "./init/ppm_rc_dump_parameters.f"

#define __2D  2
#define __3D  3

#define  DTYPE(a) a/**/_3d
#define __DIME  __3D
#include "./init/ppm_rc_init_seg.f"
#include "./init/ppm_rc_init_mcmc.f"
#undef  __DIME
#undef  DTYPE
#define  DTYPE(a) a/**/_2d
#define __DIME  __2D
#include "./init/ppm_rc_init_seg.f"
#include "./init/ppm_rc_init_mcmc.f"
#undef  __DIME
#undef  DTYPE

#undef  __2D
#undef  __3D

      END MODULE ppm_rc_module_init
