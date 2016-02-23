      !-------------------------------------------------------------------------
      !  Module       :                    ppm_rc_module_energy
      !-------------------------------------------------------------------------
      !
      !  Purpose      :  Module contains energy subroutines.
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
      MODULE ppm_rc_module_energy
        !-----------------------------------------------------------------------
        !  Modules
        !-----------------------------------------------------------------------
        USE ppm_rc_module_global
        USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list
        IMPLICIT NONE

        PRIVATE
        !-----------------------------------------------------------------------
        !  Some work memory on the heap
        !-----------------------------------------------------------------------
        REAL(ppm_kind_double), DIMENSION(:),   ALLOCATABLE :: bufr
        REAL(MK),              DIMENSION(:),   ALLOCATABLE :: energyt

        INTEGER,  DIMENSION(:),   ALLOCATABLE :: bufi
        INTEGER,  DIMENSION(:),   ALLOCATABLE :: labelt
        INTEGER,  DIMENSION(:),   ALLOCATABLE :: candlabelt
        INTEGER,  DIMENSION(:),   ALLOCATABLE :: ccandlabelt
        INTEGER,  DIMENSION(:),   ALLOCATABLE :: ndaughterst
        INTEGER,  DIMENSION(:),   ALLOCATABLE :: nmotherst
        INTEGER,  DIMENSION(:,:), ALLOCATABLE :: motherst
        INTEGER,  DIMENSION(:,:), ALLOCATABLE :: daughterst
        INTEGER,  DIMENSION(:),   ALLOCATABLE :: acceptedt
        INTEGER,  DIMENSION(:),   ALLOCATABLE :: missedparticles

#include "./energy/ppm_rc_energytypedef.f"

        CLASS(RCExternalEnergyBaseClass), POINTER :: e_data   => NULL()
        CLASS(RCInternalEnergyBaseClass), POINTER :: e_length => NULL()

        INTEGER :: e_dX,e_dY,e_dZ
        INTEGER :: e_lX,e_lY,e_lZ

#ifdef __MPI
        INTEGER :: requestCount
        INTEGER :: requestSums
#endif

        !----------------------------------------------------------------------
        ! Public
        !----------------------------------------------------------------------
        PUBLIC :: E_Gamma
        PUBLIC :: E_ContourLengthApprox
        PUBLIC :: E_PC
        PUBLIC :: E_PCGaussian
        PUBLIC :: E_PCPoisson
        PUBLIC :: E_PS
        PUBLIC :: E_PSGaussian
        PUBLIC :: E_PSPoisson

        PUBLIC :: e_data
        PUBLIC :: e_length

        PUBLIC :: missedparticles

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        !----------------------------------------------------------------------
        ! Public
        !----------------------------------------------------------------------
        PUBLIC :: ppm_rc_energy_compute_2d
        PUBLIC :: ppm_rc_energy_compute_3d

      CONTAINS

#define __2D  2
#define __3D  3

#define  DTYPE(a) a/**/_3d
#define __DIME  __3D
#include "./energy/ppm_rc_energytypeproc.f"
#include "./energy/ppm_rc_energytype_PC.f"
#include "./energy/ppm_rc_energytype_PS.f"
#undef  __DIME
#undef  DTYPE

#define  DTYPE(a) a/**/_2d
#define __DIME  __2D
#include "./energy/ppm_rc_energytypeproc.f"
#include "./energy/ppm_rc_energytype_PC.f"
#include "./energy/ppm_rc_energytype_PS.f"
#undef  __DIME
#undef  DTYPE

#define  DTYPE(a) a/**/_3d
#define __DIME  __3D
#include "./energy/ppm_rc_energy_compute.f"
#undef  __DIME
#undef  DTYPE

#define  DTYPE(a) a/**/_2d
#define __DIME  __2D
#include "./energy/ppm_rc_energy_compute.f"
#undef  __DIME
#undef  DTYPE

#undef  __2D
#undef  __3D

      END MODULE ppm_rc_module_energy
